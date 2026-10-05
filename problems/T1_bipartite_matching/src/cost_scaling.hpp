#pragma once

#include "bipartite_matching.hpp"

#include <cstdint>
#include <vector>

namespace pip27 {

// 费用缩放 push-relabel（Goldberg-Tarjan cost scaling）：
// 求不限匹配数量的最大权匹配（允许空匹配）。
//
// 建模为最小费用循环：源点 -> 左部 -> 右部 -> 汇点，容量均为 1，
// 左右之间弧费用取权重的相反数，另加汇->源弧使流量值任意；
// 循环流与匹配一一对应，最小费用循环即最大权匹配。
//
// 算法维护“ε-最优流”：所有残量弧满足归约费用 c^π >= -ε。
// 每轮 ε-细化：先把 c^π < -ε 的残量弧饱和，再用 push-relabel 消化
// 产生的盈余（沿可容许弧推送，无可容许弧时抬高顶点势），恢复 ε-最优；
// 随后 ε 减半。当 ε < 1/(2n) 时收口：
//   * 逐弧校验 c^π >= -ε（线性扫描，极便宜）；
//   * 弧费用全为整数，任意残量环的费用 = Σc^π >= -n·ε > -1/2，
//     但环费用是整数，故 >= 0：残量网络无负环，循环流即全局最优。
//
// 正确性由上述整数环论证 + 线性校验闭环，与实现过程无关；预算耗尽或
// 校验不通过时返回 false，由调用方回退确定性 SSP。
class CostScalingMatcher {
public:
    CostScalingMatcher(std::uint32_t left_count, std::uint32_t right_count,
                       const std::vector<Edge>& edges);

    // 成功时写入长度为 left_count 的匹配数组（-1 表示不匹配）。
    bool solve(std::vector<std::int32_t>* match_left);

private:
    bool refine(double epsilon);
    bool verifyOptimal(double epsilon) const;

    std::uint32_t left_count_;
    std::uint32_t right_count_;
    int source_ = 0;
    int sink_ = 0;
    // 残量网络（CSR）：每条输入弧带一条反向弧，两者尾顶点不同、
    // 在各自尾顶点的区间内连续存放。
    std::vector<std::int32_t> arc_start_;
    std::vector<std::int32_t> arc_to_;
    std::vector<std::int32_t> arc_reverse_;
    std::vector<std::int64_t> arc_cost_;
    std::vector<std::int32_t> arc_capacity_;
    std::vector<std::int64_t> excess_;
    std::vector<double> potential_;
    std::vector<std::int32_t> current_;
    std::int64_t budget_ = 0;
};

}  // namespace pip27
