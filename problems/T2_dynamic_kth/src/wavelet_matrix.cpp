#include "wavelet_matrix.hpp"

#include <algorithm>
#include <cstring>
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
            chunk->rows[g][d] = static_cast<std::uint16_t>(running[d]);
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
        const std::uint16_t* row = chunk.rows[full_rows - 1];
        for (int d = 0; d < 16; ++d) {
            out[d] = row[d];
        }
    }
    for (std::uint32_t i = full_rows * kRowSymbols; i < count; ++i) {
        ++out[chunk.digits[i]];
    }
}

void WaveletMatrix::fenwickRebuild(Level* level) {
    const std::size_t n = level->chunks.size();
    level->size_fenwick.assign(n + 1, 0);
    level->count_fenwick.assign((n + 1) * 16, 0);
    for (std::size_t i = 0; i < n; ++i) {
        level->size_fenwick[i + 1] = level->chunks[i].used;
        std::uint32_t* node = level->count_fenwick.data() + (i + 1) * 16;
        const std::uint16_t* row = level->chunks[i].rows[kRows - 1];
        for (int d = 0; d < 16; ++d) {
            node[d] = row[d];
        }
    }
    for (std::size_t i = 1; i <= n; ++i) {
        const std::size_t parent = i + lowBit(i);
        if (parent <= n) {
            level->size_fenwick[parent] += level->size_fenwick[i];
            std::uint32_t* target =
                level->count_fenwick.data() + parent * 16;
            const std::uint32_t* source =
                level->count_fenwick.data() + i * 16;
            for (int k = 0; k < 16; ++k) {
                target[k] += source[k];
            }
        }
    }
}

void WaveletMatrix::fenwickReplaceDigit(Level* level, std::size_t chunk_index,
                                        int old_digit, int new_digit) {
    const std::size_t n = level->chunks.size();
    for (std::size_t i = chunk_index + 1; i <= n; i += lowBit(i)) {
        --level->count_fenwick[i * 16 + old_digit];
        ++level->count_fenwick[i * 16 + new_digit];
    }
}

void WaveletMatrix::fenwickAdd(Level* level, std::size_t chunk_index,
                               int digit, int delta) {
    const std::size_t n = level->chunks.size();
    for (std::size_t i = chunk_index + 1; i <= n; i += lowBit(i)) {
        level->size_fenwick[i] += static_cast<std::uint32_t>(delta);
        level->count_fenwick[i * 16 + digit] += static_cast<std::uint32_t>(delta);
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
        if (next <= n && level.size_fenwick[next] <= remaining) {
            index = next;
            remaining -= level.size_fenwick[next];
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

std::uint32_t WaveletMatrix::rankOf(int level, std::uint32_t position,
                                    int digit, std::size_t* chunk_index,
                                    std::uint32_t* offset) const {
    const Level& target = levels_[level];
    locate(target, position, chunk_index, offset);
    const std::size_t ci = *chunk_index;
    const std::uint32_t off = *offset;
    if (ci >= target.chunks.size()) {
        return target.total[digit];
    }
    std::uint32_t count = 0;
    for (std::size_t i = ci; i > 0; i -= lowBit(i)) {
        count += target.count_fenwick[i * 16 + digit];
    }
    const Chunk& chunk = target.chunks[ci];
    const std::uint32_t row = off / kRowSymbols;
    if (row > 0) {
        count += chunk.rows[row - 1][digit];
    }
    for (std::uint32_t i = row * kRowSymbols; i < off; ++i) {
        if (chunk.digits[i] == static_cast<std::uint32_t>(digit)) {
            ++count;
        }
    }
    return count;
}

std::uint32_t WaveletMatrix::locateDigitRank(int level,
                                             std::uint32_t position,
                                             std::uint32_t* same_before,
                                             std::size_t* chunk_index,
                                             std::uint32_t* offset) const {
    const Level& target = levels_[level];
    locate(target, position, chunk_index, offset);
    std::size_t ci = *chunk_index;
    std::uint32_t off = *offset;
    if (ci >= target.chunks.size()) {
        ci = target.chunks.size() - 1;
        off = target.chunks[ci].used - 1;
        *chunk_index = ci;
        *offset = off;
    }
    const Chunk& chunk = target.chunks[ci];
    const std::uint32_t digit = chunk.digits[off];
    std::uint32_t count = 0;
    for (std::size_t i = ci; i > 0; i -= lowBit(i)) {
        count += target.count_fenwick[i * 16 + digit];
    }
    const std::uint32_t row = off / kRowSymbols;
    if (row > 0) {
        count += chunk.rows[row - 1][digit];
    }
    for (std::uint32_t i = row * kRowSymbols; i < off; ++i) {
        if (chunk.digits[i] == digit) {
            ++count;
        }
    }
    *same_before = count;
    return digit;
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
            const std::uint32_t* node = target.count_fenwick.data() + i * 16;
            for (int d = 0; d < 16; ++d) {
                out[d] += node[d];
            }
        }
        return 0;
    }
    const Chunk& chunk = target.chunks[ci];
    for (std::size_t i = ci; i > 0; i -= lowBit(i)) {
        const std::uint32_t* node = target.count_fenwick.data() + i * 16;
        for (int d = 0; d < 16; ++d) {
            out[d] += node[d];
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
        const std::uint32_t* node = level.count_fenwick.data() + i * 16;
        for (int d = 0; d < 16; ++d) {
            out[d] += node[d];
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
        const std::uint32_t* node = level.count_fenwick.data() + i * 16;
        for (int d = 0; d < 16; ++d) {
            out[d] += node[d];
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
    const std::uint32_t first_row = offset / kRowSymbols;
    for (std::uint32_t g = first_row; g < kRows; ++g) {
        std::uint16_t* row = chunk.rows[g];
        --row[old_digit];
        ++row[new_digit];
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
        std::memcpy(tail.digits, head.digits + half, tail.used);
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
    // 块内搬移用 memmove: 逐字节循环在部分编译器上无法向量化,
    // 而这是修改路径最大的内存流量来源。
    std::memmove(chunk.digits + offset + 1, chunk.digits + offset,
                 chunk.used - offset);
    chunk.digits[offset] = static_cast<std::uint8_t>(digit);
    ++chunk.used;
    // 累计行更新: 新数字计入所有包含它的行; 被后移挤出某行上界的那个
    // 符号要从该行扣除。一趟行扫描同时完成两件事,顺序访存。
    const std::uint32_t first_row = offset / kRowSymbols;
    for (std::uint32_t g = first_row; g < kRows; ++g) {
        std::uint16_t* row = chunk.rows[g];
        ++row[digit];
        const std::uint32_t boundary =
            (g + 1) * static_cast<std::uint32_t>(kRowSymbols);
        if (boundary - 1 >= offset && boundary - 1 < previous_used) {
            --row[chunk.digits[boundary]];
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
    std::memmove(chunk.digits + offset, chunk.digits + offset + 1,
                 chunk.used - offset - 1);
    --chunk.used;
    // 删除同理: 被前移拉进某行上界内的符号要补回该行。
    const std::uint32_t first_row = offset / kRowSymbols;
    for (std::uint32_t g = first_row; g < kRows; ++g) {
        std::uint16_t* row = chunk.rows[g];
        --row[digit];
        const std::uint32_t boundary =
            (g + 1) * static_cast<std::uint32_t>(kRowSymbols);
        if (boundary > offset && boundary < previous_used) {
            ++row[chunk.digits[boundary - 1]];
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
            std::memcpy(chunk.digits + chunk.used, next.digits, next.used);
            chunk.used += next.used;
            rebuildChunk(&chunk);
            target.chunks.erase(target.chunks.begin() +
                                static_cast<std::ptrdiff_t>(chunk_index) + 1);
            fenwickRebuild(&target);
        } else if (chunk_index > 0 &&
                   chunk.used + target.chunks[chunk_index - 1].used <=
                       kChunkSymbols) {
            Chunk& previous = target.chunks[chunk_index - 1];
            std::memcpy(previous.digits + previous.used, chunk.digits,
                        chunk.used);
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
        const int digit = static_cast<int>(digitOf(key, level));
        std::size_t ci = 0;
        std::uint32_t off = 0;
        std::uint32_t same_before = 0;
        if (!levels_[level].chunks.empty()) {
            same_before = rankOf(level, p, digit, &ci, &off);
            if (ci >= levels_[level].chunks.size()) {
                ci = levels_[level].chunks.size() - 1;
                off = levels_[level].chunks[ci].used;
            }
        }
        insertDigit(level, ci, off, static_cast<std::uint32_t>(digit));
        std::uint32_t base = 0;
        for (int d = 0; d < digit; ++d) {
            base += levels_[level].total[d];
        }
        p = base + same_before;
    }
    ++size_;
}

void WaveletMatrix::erase(std::size_t position) {
    std::uint32_t p = static_cast<std::uint32_t>(position);
    for (int level = 0; level < kLevels; ++level) {
        std::size_t ci = 0;
        std::uint32_t off = 0;
        std::uint32_t same_before = 0;
        const std::uint32_t digit =
            locateDigitRank(level, p, &same_before, &ci, &off);
        eraseDigit(level, ci, off);
        std::uint32_t base = 0;
        for (int d = 0; d < static_cast<int>(digit); ++d) {
            base += levels_[level].total[d];
        }
        p = base + same_before;
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
        std::size_t ci = 0;
        std::uint32_t off = 0;
        std::uint32_t same_old = 0;
        const std::uint32_t old_digit =
            locateDigitRank(level, old_p, &same_old, &ci, &off);
        const std::uint32_t new_digit = digitOf(new_key, level);
        std::uint32_t base = 0;
        for (int d = 0; d < static_cast<int>(old_digit); ++d) {
            base += levels_[level].total[d];
        }
        const std::uint32_t routed = base + same_old;
        if (old_digit != new_digit) {
            const std::uint32_t same_new =
                rankOf(level, old_p, static_cast<int>(new_digit), &ci, &off);
            std::uint32_t base_new = 0;
            for (int d = 0; d < static_cast<int>(new_digit); ++d) {
                base_new += levels_[level].total[d];
            }
            if (old_digit < new_digit) {
                --base_new;
            }
            const std::uint32_t routed_new = base_new + same_new;
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
        std::size_t ci_old = 0;
        std::uint32_t off_old = 0;
        std::uint32_t same_old = 0;
        const std::uint32_t old_digit =
            locateDigitRank(level, old_p, &same_old, &ci_old, &off_old);
        eraseDigit(level, ci_old, off_old);
        std::uint32_t base_old = 0;
        for (int d = 0; d < static_cast<int>(old_digit); ++d) {
            base_old += levels_[level].total[d];
        }
        const std::uint32_t routed_old = base_old + same_old;
        // routed_new 就是新符号在本层"最终序列"中的下标,而最终序列恰为
        // "删除旧符号后的序列 + 在此下标插入新符号",故直接在此插入即可。
        const int new_digit = static_cast<int>(digitOf(new_key, level));
        std::size_t ci_new = 0;
        std::uint32_t off_new = 0;
        std::uint32_t same_new = 0;
        if (!levels_[level].chunks.empty()) {
            same_new = rankOf(level, new_p, new_digit, &ci_new, &off_new);
            if (ci_new >= levels_[level].chunks.size()) {
                ci_new = levels_[level].chunks.size() - 1;
                off_new = levels_[level].chunks[ci_new].used;
            }
        }
        insertDigit(level, ci_new, off_new,
                    static_cast<std::uint32_t>(new_digit));
        std::uint32_t base_new = 0;
        for (int d = 0; d < new_digit; ++d) {
            base_new += levels_[level].total[d];
        }
        old_p = routed_old;
        new_p = base_new + same_new;
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
