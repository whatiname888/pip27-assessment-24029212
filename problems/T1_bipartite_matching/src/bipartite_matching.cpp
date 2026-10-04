#include "bipartite_matching.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

// T1：带权二分图最大权匹配。
//
// 算法：把最大权匹配转化为最小费用流（源点 -> 左部 -> 右部 -> 汇点，容量均为 1，
// 左右之间边的费用取权重的相反数），用 primal-dual 方法逐相位求增广路：
//   1. 用势函数保证残量网络中所有弧的归约费用非负，再用 Dijkstra 求最短增广路；
//   2. 若最短增广路的实际费用 >= 0，说明继续增广不会更优，立即停止；
//   3. 否则在最短路 DAG 上找出极大顶点不相交的紧弧路径集合，一次全部增广
//      （这些路径两两不交，翻转后新反向弧的归约费用恰为 0，势函数依然可行，
//      等价于若干次顺序增广），再按最短距离抬高势函数，进入下一相位。
// 最小费用关于流量值是凸函数，因此"增广到边际费用非负为止"得到的正是不限
// 匹配数量时的最优匹配（允许空匹配）。最后按题目要求补上非空约束的后处理。
//
// 实现要点：弧用 SoA 分列存储（松弛只读 to/cost/cap，翻转才读 reverse），
// Dijkstra 用 4 叉索引堆降低 sift 深度，顶点数上限约 4002，势函数与距离表
// 常驻缓存。
namespace {

using Cost = std::int64_t;

constexpr Cost kInfCost = std::numeric_limits<Cost>::max() / 4;

// 4 叉索引堆，支持 decrease-key；键值保存在外部数组。
class IndexHeap {
public:
    explicit IndexHeap(const std::vector<Cost>* keys) : keys_(keys) {}

    void reset(std::size_t vertex_count) {
        positions_.assign(vertex_count, -1);
        heap_.clear();
    }

    bool empty() const { return heap_.empty(); }

    void pushOrDecrease(int vertex) {
        const std::size_t index = positions_[static_cast<std::size_t>(vertex)];
        if (static_cast<std::int64_t>(index) < 0) {
            positions_[static_cast<std::size_t>(vertex)] =
                static_cast<std::int32_t>(heap_.size());
            heap_.push_back(vertex);
            siftUp(heap_.size() - 1);
        } else {
            siftUp(static_cast<std::size_t>(index));
        }
    }

    int pop() {
        const int top = heap_.front();
        positions_[static_cast<std::size_t>(top)] = -1;
        const int last = heap_.back();
        heap_.pop_back();
        if (!heap_.empty()) {
            heap_[0] = last;
            positions_[static_cast<std::size_t>(last)] = 0;
            siftDown(0);
        }
        return top;
    }

private:
    Cost keyAt(std::size_t index) const {
        return (*keys_)[static_cast<std::size_t>(heap_[index])];
    }

    void siftUp(std::size_t index) {
        const int vertex = heap_[index];
        const Cost key = (*keys_)[static_cast<std::size_t>(vertex)];
        std::size_t current = index;
        while (current > 0) {
            const std::size_t parent = (current - 1) / 4;
            if (keyAt(parent) <= key) {
                break;
            }
            moveAt(parent, current);
            current = parent;
        }
        heap_[current] = vertex;
        positions_[static_cast<std::size_t>(vertex)] =
            static_cast<std::int32_t>(current);
    }

    void siftDown(std::size_t index) {
        const int vertex = heap_[index];
        const Cost key = (*keys_)[static_cast<std::size_t>(vertex)];
        const std::size_t size = heap_.size();
        std::size_t current = index;
        while (true) {
            const std::size_t child = current * 4 + 1;
            if (child >= size) {
                break;
            }
            std::size_t best = child;
            const std::size_t end = std::min(size, child + 4);
            for (std::size_t candidate = child + 1; candidate < end;
                 ++candidate) {
                if (keyAt(candidate) < keyAt(best)) {
                    best = candidate;
                }
            }
            if (keyAt(best) >= key) {
                break;
            }
            moveAt(best, current);
            current = best;
        }
        heap_[current] = vertex;
        positions_[static_cast<std::size_t>(vertex)] =
            static_cast<std::int32_t>(current);
    }

    void moveAt(std::size_t from, std::size_t to) {
        heap_[to] = heap_[from];
        positions_[static_cast<std::size_t>(heap_[to])] =
            static_cast<std::int32_t>(to);
    }

    const std::vector<Cost>* keys_;
    std::vector<int> heap_;
    std::vector<std::int32_t> positions_;
};

}  // namespace

struct BipartiteMatcher::Impl {
    Impl(std::uint32_t left_count, std::uint32_t right_count,
         const std::vector<Edge>& edges)
        : left_count(left_count), right_count(right_count), edges(edges) {
        const int left_begin = 1;
        left_end = left_begin + static_cast<int>(left_count);
        source = 0;
        sink = left_end + static_cast<int>(right_count);
        const int node_count = sink + 1;

        std::vector<std::pair<int, int>> endpoints;
        endpoints.reserve(edges.size());
        for (const Edge& edge : edges) {
            endpoints.emplace_back(left_begin + static_cast<int>(edge.left),
                                   left_end + static_cast<int>(edge.right));
        }

        // 每条输入边拆成正反两条弧，先统计各顶点出度，再按尾顶点摊平。
        std::vector<int> starts(static_cast<std::size_t>(node_count) + 1, 0);
        auto count_tail = [&starts](int tail) {
            ++starts[static_cast<std::size_t>(tail) + 1];
        };
        for (int left = 0; left < static_cast<int>(left_count); ++left) {
            count_tail(source);
            count_tail(left + left_begin);
        }
        for (const auto& pair : endpoints) {
            count_tail(pair.first);
            count_tail(pair.second);
        }
        for (int right = 0; right < static_cast<int>(right_count); ++right) {
            count_tail(left_end + right);
            count_tail(sink);
        }
        for (int node = 0; node < node_count; ++node) {
            starts[static_cast<std::size_t>(node) + 1] +=
                starts[static_cast<std::size_t>(node)];
        }
        const std::size_t arc_count = starts.back();
        arc_to.assign(arc_count, 0);
        arc_reverse.assign(arc_count, 0);
        arc_cost.assign(arc_count, 0);
        arc_capacity.assign(arc_count, 0);
        std::vector<int> cursor(starts.begin(),
                                starts.begin() + node_count);
        auto add_arc = [this, &cursor](int from, int to, Cost cost) {
            const int forward = cursor[static_cast<std::size_t>(from)]++;
            const int backward = cursor[static_cast<std::size_t>(to)]++;
            arc_to[static_cast<std::size_t>(forward)] = to;
            arc_reverse[static_cast<std::size_t>(forward)] = backward;
            arc_cost[static_cast<std::size_t>(forward)] = cost;
            arc_capacity[static_cast<std::size_t>(forward)] = 1;
            arc_to[static_cast<std::size_t>(backward)] = from;
            arc_reverse[static_cast<std::size_t>(backward)] = forward;
            arc_cost[static_cast<std::size_t>(backward)] = -cost;
            arc_capacity[static_cast<std::size_t>(backward)] = 0;
        };
        for (int left = 0; left < static_cast<int>(left_count); ++left) {
            add_arc(source, left + left_begin, 0);
        }
        for (std::size_t i = 0; i < endpoints.size(); ++i) {
            add_arc(endpoints[i].first, endpoints[i].second, -edges[i].weight);
        }
        for (int right = 0; right < static_cast<int>(right_count); ++right) {
            add_arc(left_end + right, sink, 0);
        }
        this->starts.assign(starts.begin(), starts.end());

        // 初始势函数：源点、左部取 0，右部取所有入边费用的最小值，
        // 汇点取各右部势与 0 的较小者，保证归约费用非负。
        potentials.assign(static_cast<std::size_t>(node_count), 0);
        for (std::size_t i = 0; i < endpoints.size(); ++i) {
            const int right_node = endpoints[i].second;
            potentials[static_cast<std::size_t>(right_node)] =
                std::min(potentials[static_cast<std::size_t>(right_node)],
                         -edges[i].weight);
        }
        Cost sink_potential = 0;
        for (int right = 0; right < static_cast<int>(right_count); ++right) {
            sink_potential =
                std::min(sink_potential,
                         potentials[static_cast<std::size_t>(left_end + right)]);
        }
        potentials[static_cast<std::size_t>(sink)] = sink_potential;

        distances.assign(static_cast<std::size_t>(node_count), 0);
    }

    // 求不限匹配数量时的最优匹配（允许空匹配），结果写入 match_left。
    void solveMatching(std::vector<std::int32_t>* match_left) {
        const int node_count = sink + 1;
        const std::size_t node_size = static_cast<std::size_t>(node_count);
        IndexHeap heap(&distances);
        std::vector<bool> occupied(node_size, false);
        std::vector<int> path;
        path.reserve(node_size);
        Cost* const distance = distances.data();
        const Cost* const potential = potentials.data();
        const int* const arc_start = starts.data();
        const int* const to = arc_to.data();
        const Cost* const cost = arc_cost.data();
        std::int8_t* const capacity = arc_capacity.data();

        while (true) {
            // ---- Dijkstra：归约费用非负，求源点到各点的最短距离 ----
            std::fill(distances.begin(), distances.end(), kInfCost);
            distance[source] = 0;
            heap.reset(node_size);
            heap.pushOrDecrease(source);
            while (!heap.empty()) {
                const int node = heap.pop();
                const Cost base = distance[node];
                for (int i = arc_start[node]; i < arc_start[node + 1]; ++i) {
                    if (capacity[i] == 0) {
                        continue;
                    }
                    const int next = to[i];
                    const Cost candidate =
                        base + cost[i] + potential[node] - potential[next];
                    if (candidate < distance[next]) {
                        distance[next] = candidate;
                        heap.pushOrDecrease(next);
                    }
                }
            }

            if (distance[sink] >= kInfCost) {
                break;
            }
            // 增广路实际费用 = 归约距离 + 势差；非负则当前流量已最优。
            if (distance[sink] + potential[sink] - potential[source] >= 0) {
                break;
            }

            // ---- 在最短路 DAG 上找极大顶点不相交的紧弧路径并增广 ----
            std::fill(occupied.begin(), occupied.end(), false);
            occupied[source] = true;
            while (true) {
                path.clear();
                if (!search(source, &path, occupied)) {
                    break;
                }
                for (int arc_index : path) {
                    --capacity[arc_index];
                    ++capacity[arc_reverse[static_cast<std::size_t>(arc_index)]];
                }
            }

            // ---- 抬高势函数，恢复归约费用非负 ----
            // 可达顶点加真实距离；不可达顶点统一加汇点距离，保持所有弧可行。
            const Cost shift = distance[sink];
            for (int node = 0; node < node_count; ++node) {
                potentials[static_cast<std::size_t>(node)] +=
                    distance[node] < kInfCost ? distance[node] : shift;
            }
        }

        match_left->assign(left_count, -1);
        for (int left = 0; left < static_cast<int>(left_count); ++left) {
            const int node = 1 + left;
            for (int i = arc_start[node]; i < arc_start[node + 1]; ++i) {
                // 左部 -> 右部的正向弧初始容量为 1，变为 0 即被选中。
                if (to[i] >= left_end && to[i] < sink && capacity[i] == 0) {
                    (*match_left)[static_cast<std::size_t>(left)] =
                        static_cast<std::int32_t>(to[i] - left_end);
                    break;
                }
            }
        }
    }

    // 沿紧弧 DFS 找一条源点到汇点的路径（避开已用顶点），弧下标写入 path。
    bool search(int node, std::vector<int>* path,
                std::vector<bool>& occupied) const {
        if (node == sink) {
            return true;
        }
        for (int i = starts[static_cast<std::size_t>(node)];
             i < starts[static_cast<std::size_t>(node) + 1]; ++i) {
            if (arc_capacity[static_cast<std::size_t>(i)] == 0) {
                continue;
            }
            const int next = arc_to[static_cast<std::size_t>(i)];
            if (next == source) {
                continue;
            }
            // 紧弧判定：dist[v] == dist[u] + 归约费用。
            const Cost reduced =
                arc_cost[static_cast<std::size_t>(i)] +
                potentials[static_cast<std::size_t>(node)] -
                potentials[static_cast<std::size_t>(next)];
            if (distances[static_cast<std::size_t>(next)] !=
                distances[static_cast<std::size_t>(node)] + reduced) {
                continue;
            }
            if (next != sink && occupied[static_cast<std::size_t>(next)]) {
                continue;
            }
            if (next != sink) {
                occupied[static_cast<std::size_t>(next)] = true;
            }
            path->push_back(i);
            if (search(next, path, occupied)) {
                return true;
            }
            path->pop_back();
        }
        return false;
    }

    std::uint32_t left_count;
    std::uint32_t right_count;
    std::vector<Edge> edges;
    int source = 0;
    int left_end = 0;
    int sink = 0;
    std::vector<int> arc_to;
    std::vector<int> arc_reverse;
    std::vector<Cost> arc_cost;
    std::vector<std::int8_t> arc_capacity;
    std::vector<int> starts;
    std::vector<Cost> potentials;
    std::vector<Cost> distances;
};

BipartiteMatcher::BipartiteMatcher(std::uint32_t left_count,
                                   std::uint32_t right_count,
                                   const std::vector<Edge>& edges)
    : impl_(std::make_unique<Impl>(left_count, right_count, edges)) {}

BipartiteMatcher::~BipartiteMatcher() = default;

std::vector<std::int32_t> BipartiteMatcher::maximumWeightMatching() const {
    Impl& impl = *impl_;
    if (impl.edges.empty()) {
        // 无边时不存在可行的非空匹配，按题意返回全 -1。
        return std::vector<std::int32_t>(impl.left_count, -1);
    }

    std::vector<std::int32_t> answer;
    impl.solveMatching(&answer);

    // 有边时必须返回非空匹配。允许空匹配时最优解为空，说明所有匹配权重 <= 0：
    // 此时权重最大的单条边就是最优非空匹配（零权边同样满足最优）。
    bool empty = true;
    for (std::int32_t matched : answer) {
        if (matched != -1) {
            empty = false;
            break;
        }
    }
    if (empty) {
        const Edge* best = &impl.edges.front();
        for (const Edge& edge : impl.edges) {
            if (edge.weight > best->weight) {
                best = &edge;
            }
        }
        answer[best->left] = static_cast<std::int32_t>(best->right);
    }
    return answer;
}
