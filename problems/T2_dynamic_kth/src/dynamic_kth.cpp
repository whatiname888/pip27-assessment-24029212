#include "dynamic_kth.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <vector>

// T2：动态序列区间第 k 大。
//
// 数据结构：位置维度用块状链表（vector of blocks），每块保存位置序数组与有序
// 数组；值域维度把 [-1e9, 1e9] 均匀划分为 W 个值域桶。
//   * 每块维护 offsets：各值域桶在其有序数组中的起止下标，用于 O(1) 取出
//     "本块中落在某个桶内的有序子段"；
//   * 全局维护桶前缀和 prefix[b][i]：前 i 块中桶 b 的元素个数，用于 O(1) 得到
//     任意块区间内各桶的元素分布（内容修改只标记脏桶，下次查询前惰性重建）。
//
// 查询 [l, r] 第 k 大：
//   1. 两端不足整块的元素直接扫描（同时得到它们的桶分布）；
//   2. 用桶前缀和得到整块部分的桶分布，自大向小累计找到目标桶；
//   3. 收集区间内落在目标桶中的候选元素（整块部分用 offsets 直接取有序子段），
//      对候选做一次 std::nth_element 即可。
// 单次查询代价约 O(B + W + 块数 + 候选数)，插入删除只需 O(B + W) 的内存搬移。
// 候选数异常偏大（数据值域分布极不均匀）时退化为在目标桶内二分答案，保证正确。
namespace {

constexpr std::int64_t kMinValue = -1000000000;
constexpr std::int64_t kMaxValue = 1000000000;

inline int clampInt(int value, int low, int high) {
    return std::max(low, std::min(high, value));
}

}  // namespace

struct DynamicKth::Impl {
    struct Block {
        std::vector<int> values;              // 按位置顺序存放
        std::vector<int> sorted;              // 升序副本
        std::vector<std::uint16_t> offsets;   // 大小 W+1，桶 b 占 sorted[o[b], o[b+1])
    };

    std::vector<Block> blocks_;
    std::size_t length_ = 0;
    int bucket_count_ = 256;
    std::size_t block_target_ = 256;
    // prefix_[b * (block_count + 1) + i] = 前 i 块中桶 b 的元素个数。
    std::vector<std::uint32_t> prefix_;
    std::vector<std::uint8_t> dirty_;
    // 查询期复用的临时缓冲，避免每次查询重新分配。
    std::vector<int> partial_;
    std::vector<std::uint32_t> partial_counts_;
    std::vector<int> candidates_;

    static constexpr std::size_t kMaxCandidates = 262144;

    explicit Impl(const std::vector<int>& initial) {
        // 参数由序列规模标定（实测最优）：块大小取 2.75*sqrt(N)，值域桶数取
        // sqrt(N)，使单次查询的残块扫描、桶扫描、整块遍历与候选数四项均衡。
        // 初始规模决定参数；序列长度最多再增长操作数次，不会明显偏离。
        const std::size_t scale = std::max<std::size_t>(initial.size(), 32768);
        const double root = std::sqrt(static_cast<double>(scale));
        block_target_ = static_cast<std::size_t>(
            clampInt(static_cast<int>(2.75 * root), 256, 8192));
        bucket_count_ = clampInt(static_cast<int>(root), 64, 4096);
        dirty_.assign(static_cast<std::size_t>(bucket_count_), 1);
        partial_counts_.assign(static_cast<std::size_t>(bucket_count_), 0);

        const std::size_t count = initial.size();
        for (std::size_t begin = 0; begin < count; begin += block_target_) {
            const std::size_t end =
                std::min(count, begin + block_target_);
            Block block;
            block.values.assign(initial.begin() + static_cast<std::ptrdiff_t>(begin),
                                initial.begin() + static_cast<std::ptrdiff_t>(end));
            blocks_.push_back(std::move(block));
        }
        for (Block& block : blocks_) {
            rebuildBlock(&block);
        }
        length_ = count;
        resizePrefix();
        ensureClean();
    }

    int bucketOf(int value) const {
        const std::uint64_t offset =
            static_cast<std::uint64_t>(static_cast<std::int64_t>(value) -
                                       kMinValue);
        return static_cast<int>(
            (offset * static_cast<std::uint64_t>(bucket_count_)) >> 31);
    }

    void rebuildBlock(Block* block) {
        block->sorted = block->values;
        std::sort(block->sorted.begin(), block->sorted.end());
        block->offsets.assign(static_cast<std::size_t>(bucket_count_) + 1, 0);
        for (int value : block->sorted) {
            ++block->offsets[static_cast<std::size_t>(bucketOf(value)) + 1];
        }
        for (int b = 0; b < bucket_count_; ++b) {
            block->offsets[static_cast<std::size_t>(b) + 1] +=
                block->offsets[static_cast<std::size_t>(b)];
        }
    }

    void markBucketDirty(int bucket) {
        dirty_[static_cast<std::size_t>(bucket)] = 1;
    }

    void markAllDirty() { std::fill(dirty_.begin(), dirty_.end(), 1); }

    void resizePrefix() {
        const std::size_t rows = blocks_.size() + 1;
        prefix_.assign(static_cast<std::size_t>(bucket_count_) * rows, 0);
        markAllDirty();
    }

    void ensureClean() {
        const std::size_t rows = blocks_.size() + 1;
        for (int b = 0; b < bucket_count_; ++b) {
            if (dirty_[static_cast<std::size_t>(b)] == 0) {
                continue;
            }
            dirty_[static_cast<std::size_t>(b)] = 0;
            std::uint32_t* row =
                &prefix_[static_cast<std::size_t>(b) * rows];
            std::uint32_t accumulated = 0;
            row[0] = 0;
            for (std::size_t i = 0; i < blocks_.size(); ++i) {
                const Block& block = blocks_[i];
                accumulated += static_cast<std::uint32_t>(
                    block.offsets[static_cast<std::size_t>(b) + 1] -
                    block.offsets[static_cast<std::size_t>(b)]);
                row[i + 1] = accumulated;
            }
        }
    }

    // index 允许等于 length_，表示追加位置。
    void locate(std::size_t index, std::size_t* block_index,
                std::size_t* offset) const {
        std::size_t remaining = index;
        for (std::size_t i = 0; i < blocks_.size(); ++i) {
            const std::size_t size = blocks_[i].values.size();
            if (remaining < size) {
                *block_index = i;
                *offset = remaining;
                return;
            }
            remaining -= size;
        }
        *block_index = blocks_.size();
        *offset = 0;
    }

    void blockInsert(Block* block, std::size_t offset, int value) {
        block->values.insert(
            block->values.begin() + static_cast<std::ptrdiff_t>(offset), value);
        const auto position = std::lower_bound(block->sorted.begin(),
                                               block->sorted.end(), value);
        block->sorted.insert(position, value);
        const int bucket = bucketOf(value);
        for (int b = bucket + 1; b <= bucket_count_; ++b) {
            ++block->offsets[static_cast<std::size_t>(b)];
        }
    }

    void blockErase(Block* block, std::size_t offset) {
        const int value = block->values[offset];
        block->values.erase(
            block->values.begin() + static_cast<std::ptrdiff_t>(offset));
        const auto position = std::lower_bound(block->sorted.begin(),
                                               block->sorted.end(), value);
        block->sorted.erase(position);
        const int bucket = bucketOf(value);
        for (int b = bucket + 1; b <= bucket_count_; ++b) {
            --block->offsets[static_cast<std::size_t>(b)];
        }
    }

    void splitBlock(std::size_t index) {
        Block& block = blocks_[index];
        const std::size_t middle = block.values.size() / 2;
        Block tail;
        tail.values.assign(block.values.begin() +
                               static_cast<std::ptrdiff_t>(middle),
                           block.values.end());
        block.values.resize(middle);
        blocks_.insert(blocks_.begin() + static_cast<std::ptrdiff_t>(index) + 1,
                       std::move(tail));
        rebuildBlock(&blocks_[index]);
        rebuildBlock(&blocks_[index + 1]);
        resizePrefix();
    }

    // 块过小时与相邻块合并；无法合并则保留小块，不影响正确性。
    void tryMerge(std::size_t index) {
        const std::size_t size = blocks_[index].values.size();
        if (blocks_.size() < 2) {
            return;
        }
        if (index > 0 &&
            size + blocks_[index - 1].values.size() <= 2 * block_target_) {
            Block& previous = blocks_[index - 1];
            previous.values.insert(previous.values.end(),
                                   blocks_[index].values.begin(),
                                   blocks_[index].values.end());
            blocks_.erase(blocks_.begin() + static_cast<std::ptrdiff_t>(index));
            rebuildBlock(&previous);
            resizePrefix();
            return;
        }
        if (index + 1 < blocks_.size() &&
            size + blocks_[index + 1].values.size() <= 2 * block_target_) {
            Block& block = blocks_[index];
            block.values.insert(block.values.end(),
                                blocks_[index + 1].values.begin(),
                                blocks_[index + 1].values.end());
            blocks_.erase(blocks_.begin() +
                          static_cast<std::ptrdiff_t>(index) + 1);
            rebuildBlock(&block);
            resizePrefix();
        }
    }

    void insert(std::size_t index, int value) {
        std::size_t block_index = 0;
        std::size_t offset = 0;
        locate(index, &block_index, &offset);
        if (block_index == blocks_.size()) {
            if (blocks_.empty()) {
                blocks_.emplace_back();
                rebuildBlock(&blocks_.back());
                resizePrefix();  // 块数变化，桶前缀和需按新行数重建
            }
            block_index = blocks_.size() - 1;
            offset = blocks_[block_index].values.size();
        }
        blockInsert(&blocks_[block_index], offset, value);
        ++length_;
        markBucketDirty(bucketOf(value));
        if (blocks_[block_index].values.size() > 2 * block_target_) {
            splitBlock(block_index);
        }
    }

    void erase(std::size_t index) {
        std::size_t block_index = 0;
        std::size_t offset = 0;
        locate(index, &block_index, &offset);
        Block& block = blocks_[block_index];
        const int value = block.values[offset];
        blockErase(&block, offset);
        --length_;
        markBucketDirty(bucketOf(value));
        if (block.values.size() * 2 < block_target_) {
            tryMerge(block_index);
        }
    }

    int kthLargest(std::size_t left, std::size_t right, std::size_t k) {
        // 一趟扫描同时定位左右端点所在块。
        std::size_t left_block = 0;
        std::size_t left_offset = 0;
        std::size_t right_block = 0;
        std::size_t right_offset = 0;
        {
            std::size_t start = 0;
            for (std::size_t i = 0; i < blocks_.size(); ++i) {
                const std::size_t size = blocks_[i].values.size();
                if (left - start < size) {
                    left_block = i;
                    left_offset = left - start;
                }
                if (right - start < size) {
                    right_block = i;
                    right_offset = right - start;
                    break;
                }
                start += size;
            }
        }

        // 查询区间落在同一块内：直接扫描小段，不经过值域桶。
        if (left_block == right_block) {
            const Block& block = blocks_[left_block];
            candidates_.assign(
                block.values.begin() + static_cast<std::ptrdiff_t>(left_offset),
                block.values.begin() + static_cast<std::ptrdiff_t>(right_offset) + 1);
            const std::size_t index = k - 1;
            std::nth_element(candidates_.begin(),
                             candidates_.begin() +
                                 static_cast<std::ptrdiff_t>(index),
                             candidates_.end(), std::greater<int>());
            return candidates_[index];
        }

        // 两端残块元素直接扫描，同时统计它们的桶分布。
        partial_.clear();
        std::fill(partial_counts_.begin(), partial_counts_.end(), 0);
        {
            const Block& head = blocks_[left_block];
            for (std::size_t i = left_offset; i < head.values.size(); ++i) {
                const int value = head.values[i];
                partial_.push_back(value);
                ++partial_counts_[static_cast<std::size_t>(bucketOf(value))];
            }
            const Block& tail = blocks_[right_block];
            for (std::size_t i = 0; i <= right_offset; ++i) {
                const int value = tail.values[i];
                partial_.push_back(value);
                ++partial_counts_[static_cast<std::size_t>(bucketOf(value))];
            }
        }

        const std::size_t first_full = left_block + 1;
        const std::size_t last_full = right_block - 1;
        const std::size_t rows = blocks_.size() + 1;
        const std::uint32_t* base = prefix_.data();

        // 自大向小累计桶内元素个数，找到包含第 k 大的值域桶。
        std::uint32_t skipped = 0;
        int target_bucket = 0;
        std::uint32_t rank_in_bucket = 0;
        for (int b = bucket_count_ - 1; b >= 0; --b) {
            std::uint32_t count = partial_counts_[static_cast<std::size_t>(b)];
            if (first_full <= last_full) {
                const std::uint32_t* row =
                    base + static_cast<std::size_t>(b) * rows;
                count += row[last_full + 1] - row[first_full];
            }
            if (skipped + count >= k) {
                target_bucket = b;
                rank_in_bucket = static_cast<std::uint32_t>(k) - skipped;
                break;
            }
            skipped += count;
        }

        // 收集候选：整块用 offsets 取出目标桶的有序子段。规模异常偏大时
        // （数据值域分布极不均匀）改为在目标桶内二分答案兜底。
        candidates_.clear();
        for (int value : partial_) {
            if (bucketOf(value) == target_bucket) {
                candidates_.push_back(value);
            }
        }
        for (std::size_t i = first_full; i <= last_full; ++i) {
            const Block& block = blocks_[i];
            const std::size_t begin =
                block.offsets[static_cast<std::size_t>(target_bucket)];
            const std::size_t end =
                block.offsets[static_cast<std::size_t>(target_bucket) + 1];
            candidates_.insert(candidates_.end(),
                               block.sorted.begin() +
                                   static_cast<std::ptrdiff_t>(begin),
                               block.sorted.begin() +
                                   static_cast<std::ptrdiff_t>(end));
            if (candidates_.size() > kMaxCandidates) {
                return kthInBucketByBinarySearch(partial_, first_full,
                                                 last_full, target_bucket,
                                                 rank_in_bucket);
            }
        }

        const std::size_t index = rank_in_bucket - 1;
        std::nth_element(candidates_.begin(),
                         candidates_.begin() + static_cast<std::ptrdiff_t>(index),
                         candidates_.end(), std::greater<int>());
        return candidates_[index];
    }

    // 候选规模过大时的兜底：在目标桶的有序子段上二分答案，只做计数。
    int kthInBucketByBinarySearch(const std::vector<int>& partial,
                                  std::size_t first_full, std::size_t last_full,
                                  int target_bucket,
                                  std::uint32_t rank_in_bucket) const {
        std::vector<int> partial_in_bucket;
        for (int value : partial) {
            if (bucketOf(value) == target_bucket) {
                partial_in_bucket.push_back(value);
            }
        }
        std::int64_t low = kMinValue - 1;
        std::int64_t high = kMaxValue;
        while (high - low > 1) {
            const std::int64_t middle = (low + high) / 2;
            const int probe = static_cast<int>(middle);
            std::uint32_t not_less = 0;
            for (int value : partial_in_bucket) {
                if (value >= probe) {
                    ++not_less;
                }
            }
            for (std::size_t i = first_full; i <= last_full; ++i) {
                const Block& block = blocks_[i];
                const std::size_t begin =
                    block.offsets[static_cast<std::size_t>(target_bucket)];
                const std::size_t end =
                    block.offsets[static_cast<std::size_t>(target_bucket) + 1];
                const auto lower =
                    std::lower_bound(block.sorted.begin() +
                                         static_cast<std::ptrdiff_t>(begin),
                                     block.sorted.begin() +
                                         static_cast<std::ptrdiff_t>(end),
                                     probe);
                not_less += static_cast<std::uint32_t>(
                    (block.sorted.begin() + static_cast<std::ptrdiff_t>(end)) -
                    lower);
            }
            if (not_less >= rank_in_bucket) {
                low = middle;
            } else {
                high = middle;
            }
        }
        return static_cast<int>(low);
    }
};

DynamicKth::DynamicKth(const std::vector<int>& initial)
    : impl_(std::make_unique<Impl>(initial)) {}

DynamicKth::~DynamicKth() = default;

void DynamicKth::insert(std::size_t index, int value) {
    impl_->insert(index, value);
}

void DynamicKth::erase(std::size_t index) { impl_->erase(index); }

void DynamicKth::update(std::size_t index, int value) {
    impl_->erase(index);
    impl_->insert(index, value);
}

int DynamicKth::kthLargest(std::size_t left, std::size_t right,
                           std::size_t k) const {
    Impl& impl = *impl_;
    impl.ensureClean();
    return impl.kthLargest(left, right, k);
}
