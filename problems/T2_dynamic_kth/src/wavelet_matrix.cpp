#include "wavelet_matrix.hpp"

#include <algorithm>
#include <iterator>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace pip27 {
namespace {

inline std::uint32_t lowBit(std::size_t value) {
    return static_cast<std::uint32_t>(value & ~(value - 1));
}

}  // namespace

void WaveletMatrix::rebuildChunk(Chunk* chunk) {
    std::uint32_t running[16] = {0};
    for (int g = 0; g < kRows; ++g) {
        const std::uint32_t begin =
            static_cast<std::uint32_t>(g) * kRowSymbols;
        const std::uint32_t end =
            std::min<std::uint32_t>(chunk->used, begin + kRowSymbols);
        for (std::uint32_t i = begin; i < end; ++i) {
            ++running[chunk->digits[i]];
        }
        for (int d = 0; d < 16; ++d) {
            chunk->rows[d][g] = static_cast<std::uint16_t>(running[d]);
        }
    }
}

void WaveletMatrix::chunkPrefix(const Chunk& chunk, std::uint32_t count,
                                std::uint32_t out[16]) {
    for (int d = 0; d < 16; ++d) {
        out[d] = 0;
    }
    if (count == 0) {
        return;
    }
    const std::uint32_t full_rows = (count - 1) / kRowSymbols;
    if (full_rows > 0) {
        for (int d = 0; d < 16; ++d) {
            out[d] = chunk.rows[d][full_rows - 1];
        }
    }
    for (std::uint32_t i = full_rows * kRowSymbols; i < count; ++i) {
        ++out[chunk.digits[i]];
    }
}

void WaveletMatrix::fenwickRebuild(Level* level) {
    const std::size_t n = level->chunks.size();
    level->fenwick.assign((n + 1) * 17, 0);
    for (std::size_t i = 0; i < n; ++i) {
        std::uint32_t* node = level->fenwick.data() + (i + 1) * 17;
        node[0] = level->chunks[i].used;
        for (int d = 0; d < 16; ++d) {
            node[1 + d] = level->chunks[i].rows[d][kRows - 1];
        }
    }
    for (std::size_t i = 1; i <= n; ++i) {
        const std::size_t parent = i + lowBit(i);
        if (parent <= n) {
            std::uint32_t* target = level->fenwick.data() + parent * 17;
            const std::uint32_t* source = level->fenwick.data() + i * 17;
            for (int k = 0; k < 17; ++k) {
                target[k] += source[k];
            }
        }
    }
}

void WaveletMatrix::fenwickReplaceDigit(Level* level, std::size_t chunk_index,
                                        int old_digit, int new_digit) {
    const std::size_t n = level->chunks.size();
    for (std::size_t i = chunk_index + 1; i <= n; i += lowBit(i)) {
        std::uint32_t* node = level->fenwick.data() + i * 17;
        --node[1 + old_digit];
        ++node[1 + new_digit];
    }
}

void WaveletMatrix::fenwickAdd(Level* level, std::size_t chunk_index,
                               int digit, int delta) {
    const std::size_t n = level->chunks.size();
    for (std::size_t i = chunk_index + 1; i <= n; i += lowBit(i)) {
        std::uint32_t* node = level->fenwick.data() + i * 17;
        node[0] += static_cast<std::uint32_t>(delta);
        node[1 + digit] += static_cast<std::uint32_t>(delta);
    }
}

void WaveletMatrix::locate(const Level& level, std::uint32_t position,
                           std::size_t* chunk_index, std::uint32_t* offset) {
    // Fenwick 上二分定位: 前 idx 个块的符号总数不超过 position。
    const std::size_t n = level.chunks.size();
    std::size_t step = 1;
    while (step <= n) {
        step <<= 1;
    }
    std::size_t index = 0;
    std::uint32_t remaining = position;
    for (step >>= 1; step > 0; step >>= 1) {
        const std::size_t next = index + step;
        if (next <= n && level.fenwick[next * 17] <= remaining) {
            index = next;
            remaining -= level.fenwick[next * 17];
        }
    }
    *chunk_index = index;
    *offset = remaining;
}

std::uint32_t WaveletMatrix::digitAt(const Level& level,
                                     std::uint32_t position) {
    std::size_t chunk_index = 0;
    std::uint32_t offset = 0;
    locate(level, position, &chunk_index, &offset);
    if (chunk_index >= level.chunks.size()) {
        chunk_index = level.chunks.size() - 1;
        offset = level.chunks[chunk_index].used - 1;
    }
    return level.chunks[chunk_index].digits[offset];
}

std::uint32_t WaveletMatrix::locateAndPrefix(int level, std::uint32_t position,
                                             std::size_t* chunk_index,
                                             std::uint32_t* offset,
                                             std::uint32_t out[16]) const {
    const Level& target = levels_[level];
    for (int d = 0; d < 16; ++d) {
        out[d] = 0;
    }
    locate(target, position, chunk_index, offset);
    const std::size_t ci = *chunk_index;
    const std::uint32_t off = *offset;
    if (ci >= target.chunks.size()) {
        // position 等于序列长度: 直方图含全部符号,数字无定义。
        for (int d = 0; d < 16; ++d) {
            out[d] = 0;
        }
        for (std::size_t i = target.chunks.size(); i > 0; i -= lowBit(i)) {
            const std::uint32_t* node = target.fenwick.data() + i * 17;
            for (int d = 0; d < 16; ++d) {
                out[d] += node[1 + d];
            }
        }
        return 0;
    }
    const Chunk& chunk = target.chunks[ci];
    for (std::size_t i = ci; i > 0; i -= lowBit(i)) {
        const std::uint32_t* node = target.fenwick.data() + i * 17;
        for (int d = 0; d < 16; ++d) {
            out[d] += node[1 + d];
        }
    }
    std::uint32_t local[16];
    chunkPrefix(chunk, off, local);
    for (int d = 0; d < 16; ++d) {
        out[d] += local[d];
    }
    return chunk.digits[off];
}

std::uint32_t WaveletMatrix::prefixAndDigit(const Level& level,
                                            std::uint32_t position,
                                            std::uint32_t out[16]) {
    for (int d = 0; d < 16; ++d) {
        out[d] = 0;
    }
    std::size_t chunk_index = 0;
    std::uint32_t offset = 0;
    locate(level, position, &chunk_index, &offset);
    if (chunk_index >= level.chunks.size()) {
        chunk_index = level.chunks.size() - 1;
        offset = level.chunks[chunk_index].used - 1;
    }
    const Chunk& chunk = level.chunks[chunk_index];
    for (std::size_t i = chunk_index; i > 0; i -= lowBit(i)) {
        const std::uint32_t* node = level.fenwick.data() + i * 17;
        for (int d = 0; d < 16; ++d) {
            out[d] += node[1 + d];
        }
    }
    std::uint32_t local[16];
    chunkPrefix(chunk, offset, local);
    for (int d = 0; d < 16; ++d) {
        out[d] += local[d];
    }
    return chunk.digits[offset];
}

void WaveletMatrix::prefixCounts(const Level& level, std::uint32_t position,
                                 std::uint32_t out[16]) {
    for (int d = 0; d < 16; ++d) {
        out[d] = 0;
    }
    if (position == 0 || level.chunks.empty()) {
        return;
    }
    std::size_t chunk_index = 0;
    std::uint32_t offset = 0;
    locate(level, position, &chunk_index, &offset);
    if (chunk_index >= level.chunks.size()) {
        // position 等于总长度: 全部块计入。
        chunk_index = level.chunks.size() - 1;
        offset = level.chunks[chunk_index].used;
    }
    for (std::size_t i = chunk_index; i > 0; i -= lowBit(i)) {
        const std::uint32_t* node = level.fenwick.data() + i * 17;
        for (int d = 0; d < 16; ++d) {
            out[d] += node[1 + d];
        }
    }
    std::uint32_t local[16];
    chunkPrefix(level.chunks[chunk_index], offset, local);
    for (int d = 0; d < 16; ++d) {
        out[d] += local[d];
    }
}

WaveletMatrix::WaveletMatrix(const std::vector<int>& initial) {
    const std::size_t count = initial.size();
    std::vector<std::uint32_t> current(count);
    for (std::size_t i = 0; i < count; ++i) {
        current[i] = keyOf(initial[i]);
    }
    std::vector<std::uint32_t> next(count);
    for (int level = 0; level < kLevels; ++level) {
        Level& target = levels_[level];
        std::size_t bucket_start[17] = {0};
        for (std::size_t begin = 0; begin < count; begin += kChunkSymbols) {
            Chunk chunk;
            chunk.used = static_cast<std::uint16_t>(
                std::min<std::size_t>(kChunkSymbols, count - begin));
            for (std::uint32_t i = 0; i < chunk.used; ++i) {
                chunk.digits[i] = static_cast<std::uint8_t>(
                    digitOf(current[begin + i], level));
            }
            rebuildChunk(&chunk);
            target.chunks.push_back(chunk);
        }
        fenwickRebuild(&target);
        for (std::size_t i = 0; i < count; ++i) {
            const std::uint32_t digit = digitOf(current[i], level);
            ++bucket_start[digit + 1];
            ++target.total[digit];
        }
        for (int d = 0; d < 16; ++d) {
            bucket_start[d + 1] += bucket_start[d];
        }
        std::size_t cursor[16];
        for (int d = 0; d < 16; ++d) {
            cursor[d] = bucket_start[d];
        }
        for (std::size_t i = 0; i < count; ++i) {
            next[cursor[digitOf(current[i], level)]++] = current[i];
        }
        target.size = static_cast<std::uint32_t>(count);
        current.swap(next);
    }
    size_ = count;
}

void WaveletMatrix::replaceDigit(int level, std::size_t chunk_index,
                                 std::uint32_t offset,
                                 std::uint32_t old_digit,
                                 std::uint32_t new_digit) {
    Level& target = levels_[level];
    Chunk& chunk = target.chunks[chunk_index];
    chunk.digits[offset] = static_cast<std::uint8_t>(new_digit);
    // 累计行只需改两种数字的计数。
    const int first_row = static_cast<int>(offset / kRowSymbols);
    for (int g = first_row; g < kRows; ++g) {
        --chunk.rows[old_digit][g];
        ++chunk.rows[new_digit][g];
    }
    fenwickReplaceDigit(&target, chunk_index, static_cast<int>(old_digit),
                        static_cast<int>(new_digit));
    --target.total[old_digit];
    ++target.total[new_digit];
}

void WaveletMatrix::insertDigit(int level, std::size_t chunk_index,
                                std::uint32_t offset, std::uint32_t digit) {
    Level& target = levels_[level];
    if (target.chunks.empty()) {
        Chunk chunk;
        chunk.used = 1;
        chunk.digits[0] = static_cast<std::uint8_t>(digit);
        rebuildChunk(&chunk);
        target.chunks.push_back(chunk);
        fenwickRebuild(&target);
        ++target.total[digit];
        ++target.size;
        return;
    }
    // 块满则对半拆开,保证插入位置所在块总有空位。
    if (target.chunks[chunk_index].used == kChunkSymbols) {
        Chunk tail;
        const std::uint32_t half = kChunkSymbols / 2;
        Chunk& head = target.chunks[chunk_index];
        tail.used = static_cast<std::uint16_t>(head.used - half);
        for (std::uint32_t i = 0; i < tail.used; ++i) {
            tail.digits[i] = head.digits[half + i];
        }
        head.used = half;
        rebuildChunk(&head);
        rebuildChunk(&tail);
        target.chunks.insert(
            target.chunks.begin() + static_cast<std::ptrdiff_t>(chunk_index) + 1,
            std::move(tail));
        fenwickRebuild(&target);
        if (offset >= half) {
            ++chunk_index;
            offset -= half;
        }
    }
    Chunk& chunk = target.chunks[chunk_index];
    const std::uint32_t previous_used = chunk.used;
    for (std::uint32_t i = chunk.used; i > offset; --i) {
        chunk.digits[i] = chunk.digits[i - 1];
    }
    chunk.digits[offset] = static_cast<std::uint8_t>(digit);
    ++chunk.used;
    // 累计行更新: 新数字计入所有包含它的行; 被后移挤出某行上界的那个
    // 符号要从该行扣除,否则行计数会虚增。
    std::uint16_t* digit_rows = chunk.rows[digit];
    const std::uint32_t first_row = offset / kRowSymbols;
    for (int g = static_cast<int>(first_row); g < kRows; ++g) {
        ++digit_rows[g];
    }
    for (int g = 0; g < kRows; ++g) {
        const std::uint32_t boundary =
            (static_cast<std::uint32_t>(g) + 1) * kRowSymbols;
        if (boundary - 1 >= offset && boundary - 1 < previous_used) {
            --chunk.rows[chunk.digits[boundary]][g];
        }
    }
    fenwickAdd(&target, chunk_index, static_cast<int>(digit), 1);
    ++target.total[digit];
    ++target.size;
}

std::uint32_t WaveletMatrix::eraseDigit(int level, std::size_t chunk_index,
                                        std::uint32_t offset) {
    Level& target = levels_[level];
    Chunk& chunk = target.chunks[chunk_index];
    const std::uint32_t digit = chunk.digits[offset];
    const std::uint32_t previous_used = chunk.used;
    for (std::uint32_t i = offset; i + 1 < chunk.used; ++i) {
        chunk.digits[i] = chunk.digits[i + 1];
    }
    --chunk.used;
    // 删除同理: 被前移拉进某行上界内的符号要补回该行。
    std::uint16_t* digit_rows = chunk.rows[digit];
    const std::uint32_t first_row = offset / kRowSymbols;
    for (int g = static_cast<int>(first_row); g < kRows; ++g) {
        --digit_rows[g];
    }
    for (int g = 0; g < kRows; ++g) {
        const std::uint32_t boundary =
            (static_cast<std::uint32_t>(g) + 1) * kRowSymbols;
        if (boundary > offset && boundary < previous_used) {
            ++chunk.rows[chunk.digits[boundary - 1]][g];
        }
    }
    fenwickAdd(&target, chunk_index, static_cast<int>(digit), -1);
    --target.total[digit];
    --target.size;
    // 块过小则与邻块合并,防止小块碎片化。
    if (chunk.used < kMergeSize && target.chunks.size() > 1) {
        if (chunk_index + 1 < target.chunks.size() &&
            chunk.used + target.chunks[chunk_index + 1].used <= kChunkSymbols) {
            Chunk& next = target.chunks[chunk_index + 1];
            for (std::uint32_t i = 0; i < next.used; ++i) {
                chunk.digits[chunk.used + i] = next.digits[i];
            }
            chunk.used += next.used;
            rebuildChunk(&chunk);
            target.chunks.erase(target.chunks.begin() +
                                static_cast<std::ptrdiff_t>(chunk_index) + 1);
            fenwickRebuild(&target);
        } else if (chunk_index > 0 &&
                   chunk.used + target.chunks[chunk_index - 1].used <=
                       kChunkSymbols) {
            Chunk& previous = target.chunks[chunk_index - 1];
            for (std::uint32_t i = 0; i < chunk.used; ++i) {
                previous.digits[previous.used + i] = chunk.digits[i];
            }
            previous.used += chunk.used;
            rebuildChunk(&previous);
            target.chunks.erase(target.chunks.begin() +
                                static_cast<std::ptrdiff_t>(chunk_index));
            fenwickRebuild(&target);
        }
    }
    return digit;
}

void WaveletMatrix::insert(std::size_t position, int value) {
    const std::uint32_t key = keyOf(value);
    std::uint32_t p = static_cast<std::uint32_t>(position);
    for (int level = 0; level < kLevels; ++level) {
        const std::uint32_t digit = digitOf(key, level);
        std::uint32_t low[16];
        std::size_t ci = 0;
        std::uint32_t off = 0;
        if (levels_[level].chunks.empty()) {
            std::fill(std::begin(low), std::end(low), 0);
        } else {
            locateAndPrefix(level, p, &ci, &off, low);
            if (ci >= levels_[level].chunks.size()) {
                ci = levels_[level].chunks.size() - 1;
                off = levels_[level].chunks[ci].used;
            }
        }
        insertDigit(level, ci, off, digit);
        std::uint32_t base = 0;
        for (int d = 0; d < static_cast<int>(digit); ++d) {
            base += levels_[level].total[d];
        }
        p = base + low[digit];
    }
    ++size_;
}

void WaveletMatrix::erase(std::size_t position) {
    std::uint32_t p = static_cast<std::uint32_t>(position);
    for (int level = 0; level < kLevels; ++level) {
        std::uint32_t low[16];
        std::size_t ci = 0;
        std::uint32_t off = 0;
        const std::uint32_t digit = locateAndPrefix(level, p, &ci, &off, low);
        eraseDigit(level, ci, off);
        std::uint32_t base = 0;
        for (int d = 0; d < static_cast<int>(digit); ++d) {
            base += levels_[level].total[d];
        }
        p = base + low[digit];
    }
    --size_;
}

void WaveletMatrix::replace(std::size_t position, int value) {
    // 修改 = 用新符号替换旧符号。自高层向低层处理:
    //   * 两符号数字相同的层完全不动,只更新路由位置;
    //   * 首个数字相异的层两者位置相同,就地换数字即可;
    //   * 再往下的层两者位置已分岔,按"先删旧、再插新"处理,
    //     插入位置要按删除造成的左移调整。
    const std::uint32_t new_key = keyOf(value);
    std::uint32_t old_p = static_cast<std::uint32_t>(position);
    std::uint32_t new_p = static_cast<std::uint32_t>(position);
    int level = 0;
    for (; level < kLevels; ++level) {
        std::uint32_t low[16];
        std::size_t ci = 0;
        std::uint32_t off = 0;
        const std::uint32_t old_digit =
            locateAndPrefix(level, old_p, &ci, &off, low);
        const std::uint32_t new_digit = digitOf(new_key, level);
        std::uint32_t base = 0;
        for (int d = 0; d < static_cast<int>(old_digit); ++d) {
            base += levels_[level].total[d];
        }
        const std::uint32_t routed = base + low[old_digit];
        if (old_digit != new_digit) {
            // 就地替换: 先算好两条路由,再改数字。注意替换会把旧数字的
            // 总数减一,新符号的路由基数要按替换后的总数计算。
            std::uint32_t base_new = 0;
            for (int d = 0; d < static_cast<int>(new_digit); ++d) {
                base_new += levels_[level].total[d];
            }
            if (old_digit < new_digit) {
                --base_new;
            }
            const std::uint32_t routed_new = base_new + low[new_digit];
            replaceDigit(level, ci, off, old_digit, new_digit);
            old_p = routed;
            new_p = routed_new;
            ++level;
            break;
        }
        old_p = routed;
        new_p = routed;
    }
    for (; level < kLevels; ++level) {
        // 先删除旧符号: 路由按删除前状态计算。
        std::uint32_t low_old[16];
        std::size_t ci_old = 0;
        std::uint32_t off_old = 0;
        const std::uint32_t old_digit =
            locateAndPrefix(level, old_p, &ci_old, &off_old, low_old);
        eraseDigit(level, ci_old, off_old);
        std::uint32_t base_old = 0;
        for (int d = 0; d < static_cast<int>(old_digit); ++d) {
            base_old += levels_[level].total[d];
        }
        const std::uint32_t routed_old = base_old + low_old[old_digit];
        // routed_new 就是新符号在本层"最终序列"中的下标,而最终序列恰为
        // "删除旧符号后的序列 + 在此下标插入新符号",故直接在此插入即可。
        const std::uint32_t new_digit = digitOf(new_key, level);
        std::uint32_t low_new[16];
        std::size_t ci_new = 0;
        std::uint32_t off_new = 0;
        if (levels_[level].chunks.empty()) {
            std::fill(std::begin(low_new), std::end(low_new), 0);
        } else {
            locateAndPrefix(level, new_p, &ci_new, &off_new, low_new);
            if (ci_new >= levels_[level].chunks.size()) {
                ci_new = levels_[level].chunks.size() - 1;
                off_new = levels_[level].chunks[ci_new].used;
            }
        }
        insertDigit(level, ci_new, off_new, new_digit);
        std::uint32_t base_new = 0;
        for (int d = 0; d < static_cast<int>(new_digit); ++d) {
            base_new += levels_[level].total[d];
        }
        old_p = routed_old;
        new_p = base_new + low_new[new_digit];
    }
}

int WaveletMatrix::kthLargest(std::size_t left, std::size_t right,
                              std::size_t k) const {
    std::uint32_t l = static_cast<std::uint32_t>(left);
    std::uint32_t r = static_cast<std::uint32_t>(right);
    std::uint32_t rank = static_cast<std::uint32_t>(k);
    std::uint32_t key = 0;
    for (int level = 0; level < kLevels; ++level) {
        std::uint32_t low[16];
        std::uint32_t high[16];
        prefixCounts(levels_[level], l, low);
        prefixCounts(levels_[level], r + 1, high);
        std::uint32_t digit = 0;
        std::uint32_t cumulative = 0;
        for (int d = 15; d >= 0; --d) {
            const std::uint32_t count = high[d] - low[d];
            if (cumulative + count >= rank) {
                digit = static_cast<std::uint32_t>(d);
                rank -= cumulative;
                break;
            }
            cumulative += count;
        }
        key = (key << 4) | digit;
        std::uint32_t base = 0;
        for (int d = 0; d < static_cast<int>(digit); ++d) {
            base += levels_[level].total[d];
        }
        l = base + low[digit];
        r = base + high[digit] - 1;
    }
    return valueOfKey(key);
}

}  // namespace pip27
