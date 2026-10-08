#pragma once

#include <cstddef>
#include <cstdint>
#include <deque>
#include <vector>

namespace pip27 {

// 动态小波矩阵（wavelet matrix），4 位数字分层、共 8 层。
//
// 值映射为无符号 32 位键（异或符号位，保持大小序），自高位到低位每层取
// 4 位数字；第 L+1 层是第 L 层按该数字稳定划分（小者在前）的结果。区间
// 第 k 大逐层下探，每层只需两次前缀计数即可确定答案的 4 位并收缩区间。
//
// 每层数字序列的存储（这是性能关键）：
//   * 定长块数组：每块 512 个数字，块内每 64 个数字记录累计计数行，
//     块内前缀计数只需一次行查表加至多 63 次扫描；
//   * 块间用扁平 Fenwick 聚合（每结点 16 路计数 + 符号数），前缀计数
//     只走 log(块数)≈12 步且全部命中缓存。
// 单次 rank 约百纳秒量级，插入/删除只在块内搬移一段数字并更新两条
// Fenwick 路径，均摊代价低。存储量 O(N)：N=1e6 时约 12MB。
class WaveletMatrix {
public:
    explicit WaveletMatrix(const std::vector<int>& initial);

    void insert(std::size_t position, int value);
    void erase(std::size_t position);
    void replace(std::size_t position, int value);
    int kthLargest(std::size_t left, std::size_t right, std::size_t k) const;

    std::size_t size() const { return size_; }

private:
    static constexpr int kLevels = 8;
    static constexpr int kChunkSymbols = 2048;   // 每块符号数
    static constexpr int kRowSymbols = 128;      // 累计计数行粒度
    static constexpr int kRows = kChunkSymbols / kRowSymbols;
    static constexpr int kSplitSize = kChunkSymbols * 3 / 4;
    static constexpr int kMergeSize = kChunkSymbols / 4;

    struct Chunk {
        std::uint16_t used = 0;
        std::uint8_t digits[kChunkSymbols] = {0};
        // rows[g][d] = 前 (g+1)*kRowSymbols 个数字中数字 d 的个数(行主序)。
        // 行主序让两处热路径都顺序访存: 查询读一整行是连续 32 字节,
        // 修改的后缀增减与跨行修正沿行号单向扫描——小 L2 机器上散射
        // 访问的缓存缺失是主要开销,顺序访问可被预取器吸收。
        std::uint16_t rows[kRows][16] = {{0}};
    };

    struct Level {
        // deque: 拆块/并块在中间插入删除只搬指针块,不搬数据本体。
        // 小块化后拆块频繁,若用 vector 会反复整层搬移,得不偿失。
        std::deque<Chunk> chunks;
        // 两条扁平 Fenwick: 紧凑的符号数(定位用,常驻一二级缓存)、
        // 16 路数字计数(前缀计数用)。定位走紧凑数组,缓存足迹从
        // 数十 KB 降到数 KB,是查询与更新共同的热路径。
        std::vector<std::uint32_t> size_fenwick;
        std::vector<std::uint32_t> count_fenwick;
        std::uint32_t total[16] = {0};
        std::uint32_t size = 0;

        std::size_t nodeCount() const { return chunks.size(); }
    };

    Level levels_[kLevels];
    std::size_t size_ = 0;

    // ---- 每层基础操作 ----
    static void rebuildChunk(Chunk* chunk);
    static void chunkPrefix(const Chunk& chunk, std::uint32_t count,
                            std::uint32_t out[16]);
    // 定位第 position 个符号所在块与块内偏移。
    static void locate(const Level& level, std::uint32_t position,
                       std::size_t* chunk_index, std::uint32_t* offset);
    // [0, position) 的 16 路直方图。
    static void prefixCounts(const Level& level, std::uint32_t position,
                             std::uint32_t out[16]);
    static std::uint32_t digitAt(const Level& level, std::uint32_t position);
    // 一趟同时给出 [0,position) 直方图与 position 处数字。
    static std::uint32_t prefixAndDigit(const Level& level,
                                        std::uint32_t position,
                                        std::uint32_t out[16]);
    static void fenwickRebuild(Level* level);
    static void fenwickAdd(Level* level, std::size_t chunk_index, int digit,
                           int delta);
    static void fenwickReplaceDigit(Level* level, std::size_t chunk_index,
                                    int old_digit, int new_digit);

    // 以下三个操作均假定 (chunk_index, offset) 已由调用方定位,省去重复查找。
    // 定位并求单个数字的前缀排名(不读取数字)。
    std::uint32_t rankOf(int level, std::uint32_t position, int digit,
                         std::size_t* chunk_index, std::uint32_t* offset) const;

    void insertDigit(int level, std::size_t chunk_index, std::uint32_t offset,
                     std::uint32_t digit);
    std::uint32_t eraseDigit(int level, std::size_t chunk_index,
                             std::uint32_t offset);
    void replaceDigit(int level, std::size_t chunk_index, std::uint32_t offset,
                      std::uint32_t old_digit, std::uint32_t new_digit);
    // 定位并给出 [0,position) 直方图与 position 处数字。
    std::uint32_t locateAndPrefix(int level, std::uint32_t position,
                                  std::size_t* chunk_index,
                                  std::uint32_t* offset,
                                  std::uint32_t out[16]) const;
    // 一趟完成: 定位、取 position 处数字、求单个数字的前缀排名。
    // 插入/删除/修改的路由只需单个排名,免去 16 路直方图的整段累加。
    std::uint32_t locateDigitRank(int level, std::uint32_t position,
                                  std::uint32_t* same_before,
                                  std::size_t* chunk_index,
                                  std::uint32_t* offset) const;

    static std::uint32_t digitOf(std::uint32_t key, int level) {
        return (key >> (28 - 4 * level)) & 0xfu;
    }
    static std::uint32_t keyOf(int value) {
        return static_cast<std::uint32_t>(value) ^ 0x80000000u;
    }
    static int valueOfKey(std::uint32_t key) {
        return static_cast<int>(key ^ 0x80000000u);
    }
};

}  // namespace pip27
