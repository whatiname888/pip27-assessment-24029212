#include "cost_scaling.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace pip27 {
namespace {

// 推送与重标号的总预算；超出则回退确定性算法。
constexpr std::int64_t kOperationBudget = 400000000;
// 归约费用校验容差；配合 ε < 1/(2n) 的整数环论证仍有充足余量。
constexpr double kTolerance = 1e-5;
// relabel 偏移与可容许判定容差：抵消浮点舍入，否则“恰好 -ε”的弧永远
// 差一个 ulp，relabel 不推进会死循环。
constexpr double kSlack = 1e-7;
constexpr double kAdmissibleTolerance = 1e-6;

}  // namespace

CostScalingMatcher::CostScalingMatcher(std::uint32_t left_count,
                                       std::uint32_t right_count,
                                       const std::vector<Edge>& edges)
    : left_count_(left_count), right_count_(right_count) {
    const int left_begin = 1;
    const int left_end = left_begin + static_cast<int>(left_count);
    source_ = 0;
    sink_ = left_end + static_cast<int>(right_count);
    const int node_count = sink_ + 1;

    struct Spec {
        int from;
        int to;
        std::int64_t cost;
        std::int32_t capacity;
    };
    std::vector<Spec> specs;
    specs.reserve(left_count + right_count + edges.size() + 1);
    for (int left = 0; left < static_cast<int>(left_count); ++left) {
        specs.push_back(Spec{source_, left + left_begin, 0, 1});
    }
    for (const Edge& edge : edges) {
        specs.push_back(
            Spec{left_begin + static_cast<int>(edge.left),
                 left_end + static_cast<int>(edge.right), -edge.weight, 1});
    }
    for (int right = 0; right < static_cast<int>(right_count); ++right) {
        specs.push_back(Spec{left_end + right, sink_, 0, 1});
    }
    // 汇->源弧使流量值任意：最小费用循环 = 不限匹配数的最大权匹配。
    specs.push_back(Spec{sink_, source_, 0,
                         static_cast<std::int32_t>(left_count)});

    arc_start_.assign(static_cast<std::size_t>(node_count) + 1, 0);
    for (const Spec& spec : specs) {
        ++arc_start_[static_cast<std::size_t>(spec.from) + 1];
        ++arc_start_[static_cast<std::size_t>(spec.to) + 1];
    }
    for (int node = 0; node < node_count; ++node) {
        arc_start_[static_cast<std::size_t>(node) + 1] +=
            arc_start_[static_cast<std::size_t>(node)];
    }
    const std::size_t arc_count = arc_start_.back();
    arc_to_.resize(arc_count);
    arc_reverse_.resize(arc_count);
    arc_cost_.resize(arc_count);
    arc_capacity_.resize(arc_count);
    std::vector<std::int32_t> cursor(arc_start_.begin(),
                                     arc_start_.begin() + node_count);
    for (const Spec& spec : specs) {
        const std::int32_t forward = cursor[spec.from]++;
        const std::int32_t backward = cursor[spec.to]++;
        arc_to_[forward] = spec.to;
        arc_reverse_[forward] = backward;
        arc_cost_[forward] = spec.cost;
        arc_capacity_[forward] = spec.capacity;
        arc_to_[backward] = spec.from;
        arc_reverse_[backward] = forward;
        arc_cost_[backward] = -spec.cost;
        arc_capacity_[backward] = 0;
    }

    excess_.assign(static_cast<std::size_t>(node_count), 0);
    potential_.assign(static_cast<std::size_t>(node_count), 0.0);
    current_.assign(static_cast<std::size_t>(node_count), 0);
}

bool CostScalingMatcher::refine(double epsilon) {
    const int node_count = sink_ + 1;
    std::fill(excess_.begin(), excess_.end(), 0);

    // 1) 饱和所有 c^π < -ε 的残量弧，产生待消化的盈余。
    for (int node = 0; node < node_count; ++node) {
        for (std::int32_t a = arc_start_[static_cast<std::size_t>(node)];
             a < arc_start_[static_cast<std::size_t>(node) + 1]; ++a) {
            if (arc_capacity_[static_cast<std::size_t>(a)] == 0) {
                continue;
            }
            const int next = arc_to_[static_cast<std::size_t>(a)];
            const double reduced =
                static_cast<double>(arc_cost_[static_cast<std::size_t>(a)]) -
                potential_[static_cast<std::size_t>(node)] +
                potential_[static_cast<std::size_t>(next)];
            if (reduced < -epsilon) {
                const std::int32_t delta =
                    arc_capacity_[static_cast<std::size_t>(a)];
                arc_capacity_[static_cast<std::size_t>(a)] = 0;
                arc_capacity_[static_cast<std::size_t>(
                    arc_reverse_[static_cast<std::size_t>(a)])] += delta;
                excess_[static_cast<std::size_t>(node)] -= delta;
                excess_[static_cast<std::size_t>(next)] += delta;
            }
        }
    }
    budget_ -= arc_start_[static_cast<std::size_t>(node_count)];
    if (budget_ < 0) {
        return false;
    }

    // 2) push-relabel 消化盈余：沿可容许弧（c^π <= -ε）推送，
    //    无可容许弧时抬高顶点势制造可容许弧。
    std::vector<int> queue;
    queue.reserve(static_cast<std::size_t>(node_count));
    std::vector<bool> in_queue(static_cast<std::size_t>(node_count), false);
    for (int node = 0; node < node_count; ++node) {
        current_[static_cast<std::size_t>(node)] =
            arc_start_[static_cast<std::size_t>(node)];
        if (excess_[static_cast<std::size_t>(node)] > 0) {
            queue.push_back(node);
            in_queue[static_cast<std::size_t>(node)] = true;
        }
    }

    std::size_t head = 0;
    while (head < queue.size()) {
        const int node = queue[head++];
        in_queue[static_cast<std::size_t>(node)] = false;
        while (excess_[static_cast<std::size_t>(node)] > 0) {
            std::int32_t a = current_[static_cast<std::size_t>(node)];
            while (a < arc_start_[static_cast<std::size_t>(node) + 1]) {
                const std::size_t index = static_cast<std::size_t>(a);
                if (arc_capacity_[index] > 0) {
                    const int next = arc_to_[index];
                    const double reduced =
                        static_cast<double>(arc_cost_[index]) -
                        potential_[static_cast<std::size_t>(node)] +
                        potential_[static_cast<std::size_t>(next)];
                    if (reduced <= -epsilon + kAdmissibleTolerance) {
                        break;
                    }
                }
                ++a;
            }
            if (a < arc_start_[static_cast<std::size_t>(node) + 1]) {
                const std::size_t index = static_cast<std::size_t>(a);
                const int next = arc_to_[index];
                const std::int32_t delta = static_cast<std::int32_t>(
                    std::min(excess_[static_cast<std::size_t>(node)],
                             static_cast<std::int64_t>(
                                 arc_capacity_[index])));
                arc_capacity_[index] -= delta;
                arc_capacity_[static_cast<std::size_t>(
                    arc_reverse_[index])] += delta;
                excess_[static_cast<std::size_t>(node)] -= delta;
                excess_[static_cast<std::size_t>(next)] += delta;
                if (excess_[static_cast<std::size_t>(next)] > 0 &&
                    !in_queue[static_cast<std::size_t>(next)]) {
                    in_queue[static_cast<std::size_t>(next)] = true;
                    queue.push_back(next);
                }
                if (arc_capacity_[index] == 0) {
                    current_[static_cast<std::size_t>(node)] = a + 1;
                }
                if (--budget_ < 0) {
                    return false;
                }
            } else {
                // 无可容许弧：抬高势，使最小归约费用弧变成可容许弧。
                double best = std::numeric_limits<double>::infinity();
                for (std::int32_t scan = arc_start_[static_cast<std::size_t>(node)];
                     scan < arc_start_[static_cast<std::size_t>(node) + 1];
                     ++scan) {
                    const std::size_t index = static_cast<std::size_t>(scan);
                    if (arc_capacity_[index] == 0) {
                        continue;
                    }
                    const int next = arc_to_[index];
                    best = std::min(best,
                                    potential_[static_cast<std::size_t>(next)] +
                                        static_cast<double>(arc_cost_[index]));
                }
                if (best == std::numeric_limits<double>::infinity()) {
                    return false;  // 盈余无处可去，交给回退算法
                }
                potential_[static_cast<std::size_t>(node)] =
                    best + epsilon + kSlack;
                current_[static_cast<std::size_t>(node)] =
                    arc_start_[static_cast<std::size_t>(node)];
                if (--budget_ < 0) {
                    return false;
                }
            }
        }
    }
    return true;
}

bool CostScalingMatcher::verifyOptimal(double epsilon) const {
    // 线性扫描残量弧，逐条校验 c^π >= -ε。
    const int node_count = sink_ + 1;
    for (int node = 0; node < node_count; ++node) {
        for (std::int32_t a = arc_start_[static_cast<std::size_t>(node)];
             a < arc_start_[static_cast<std::size_t>(node) + 1]; ++a) {
            const std::size_t index = static_cast<std::size_t>(a);
            if (arc_capacity_[index] == 0) {
                continue;
            }
            const int next = arc_to_[index];
            const double reduced = static_cast<double>(arc_cost_[index]) -
                                   potential_[static_cast<std::size_t>(node)] +
                                   potential_[static_cast<std::size_t>(next)];
            if (reduced < -epsilon - kTolerance) {
                return false;
            }
        }
    }
    return true;
}

bool CostScalingMatcher::solve(std::vector<std::int32_t>* match_left) {
    const int node_count = sink_ + 1;
    if (left_count_ == 0 || right_count_ == 0) {
        match_left->assign(left_count_, -1);
        return true;
    }

    std::int64_t max_cost = 0;
    for (std::size_t a = 0; a < arc_cost_.size(); ++a) {
        max_cost = std::max(max_cost,
                            arc_cost_[a] >= 0 ? arc_cost_[a] : -arc_cost_[a]);
    }
    double epsilon = static_cast<double>(max_cost) + 1.0;
    const double final_epsilon =
        1.0 / (16.0 * static_cast<double>(node_count));
    budget_ = kOperationBudget;
    double last_epsilon = epsilon;
    while (epsilon > final_epsilon) {
        if (!refine(epsilon)) {
            return false;
        }
        last_epsilon = epsilon;
        epsilon *= 0.5;
    }

    // 收口判定：c^π >= -ε 对全部残量弧成立时，任何残量环的费用
    // Σc = Σc^π >= -n·ε > -1/2；环费用为整数故 >= 0，无负环即全局最优。
    if (!verifyOptimal(last_epsilon)) {
        return false;
    }

    // 左部 -> 右部的正向弧剩余容量为 0 表示流量为 1，即该边被选中。
    const int left_end = 1 + static_cast<int>(left_count_);
    match_left->assign(left_count_, -1);
    for (int left = 0; left < static_cast<int>(left_count_); ++left) {
        const int node = 1 + left;
        for (std::int32_t a = arc_start_[static_cast<std::size_t>(node)];
             a < arc_start_[static_cast<std::size_t>(node) + 1]; ++a) {
            const std::size_t index = static_cast<std::size_t>(a);
            const int next = arc_to_[index];
            if (next >= left_end && next < sink_ &&
                arc_capacity_[index] == 0 &&
                arc_capacity_[static_cast<std::size_t>(
                    arc_reverse_[index])] == 1) {
                (*match_left)[static_cast<std::size_t>(left)] =
                    static_cast<std::int32_t>(next - left_end);
                break;
            }
        }
    }
    (void)node_count;
    return true;
}

}  // namespace pip27
