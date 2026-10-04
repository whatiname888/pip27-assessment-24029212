#include "bipartite_matching.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <deque>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <set>
#include <string>
#include <utility>
#include <vector>

// T1 本地测评程序（不参与正式测评）：
//   * 基础档与对拍用例用位压缩 DP 暴力求非空匹配最优值；
//   * 全部用例再用一份独立的 SPFA 最小费用流参考实现重算最优值核对；
//   * 同时检查返回匹配的合法性（边存在、右部不重复、有边时非空、无边全 -1）。
// 用法：./grader [seed]
namespace {

using Cost = std::int64_t;
constexpr Cost kInfCost = std::numeric_limits<Cost>::max() / 4;

struct Graph {
    std::uint32_t left_count = 0;
    std::uint32_t right_count = 0;
    std::vector<Edge> edges;
};

// ---------- 参考实现：SPFA 逐次最短路最小费用流，允许空匹配 ----------
struct RefArc {
    int to;
    int reverse;
    Cost cost;
    int capacity;
};

Cost referenceAllowEmptyWeight(const Graph& graph) {
    const int left_begin = 1;
    const int left_end = left_begin + static_cast<int>(graph.left_count);
    const int sink = left_end + static_cast<int>(graph.right_count);
    const int node_count = sink + 1;

    std::vector<std::vector<RefArc>> network(node_count);
    auto add_edge = [&network](int from, int to, Cost cost) {
        network[from].push_back(
            RefArc{to, static_cast<int>(network[to].size()), cost, 1});
        network[to].push_back(RefArc{
            from, static_cast<int>(network[from].size()) - 1, -cost, 0});
    };
    for (int left = 0; left < static_cast<int>(graph.left_count); ++left) {
        add_edge(0, left + left_begin, 0);
    }
    for (const Edge& edge : graph.edges) {
        add_edge(left_begin + static_cast<int>(edge.left),
                 left_end + static_cast<int>(edge.right), -edge.weight);
    }
    for (int right = 0; right < static_cast<int>(graph.right_count); ++right) {
        add_edge(left_end + right, sink, 0);
    }

    Cost total = 0;
    std::vector<Cost> distance(node_count);
    std::vector<int> previous_node(node_count);
    std::vector<int> previous_arc(node_count);
    std::vector<bool> in_queue(node_count);
    std::deque<int> queue;
    while (true) {
        std::fill(distance.begin(), distance.end(), kInfCost);
        std::fill(previous_node.begin(), previous_node.end(), -1);
        std::fill(in_queue.begin(), in_queue.end(), false);
        queue.clear();
        distance[0] = 0;
        queue.push_back(0);
        in_queue[0] = true;
        while (!queue.empty()) {
            const int node = queue.front();
            queue.pop_front();
            in_queue[node] = false;
            for (std::size_t i = 0; i < network[node].size(); ++i) {
                const RefArc& arc = network[node][i];
                if (arc.capacity == 0) {
                    continue;
                }
                const Cost candidate = distance[node] + arc.cost;
                if (candidate < distance[arc.to]) {
                    distance[arc.to] = candidate;
                    previous_node[arc.to] = node;
                    previous_arc[arc.to] = static_cast<int>(i);
                    if (!in_queue[arc.to]) {
                        in_queue[arc.to] = true;
                        queue.push_back(arc.to);
                    }
                }
            }
        }
        if (distance[sink] >= kInfCost || distance[sink] >= 0) {
            break;
        }
        total += distance[sink];
        int node = sink;
        while (node != 0) {
            RefArc& arc =
                network[previous_node[node]]
                       [static_cast<std::size_t>(previous_arc[node])];
            --arc.capacity;
            ++network[node][static_cast<std::size_t>(arc.reverse)].capacity;
            node = previous_node[node];
        }
    }
    return -total;  // 费用取的是权重的相反数，取反还原最大权重。
}

// ---------- 基础档：位压缩 DP 暴力（只统计非空匹配） ----------
Cost bruteForceNonEmptyWeight(const Graph& graph) {
    const int left = static_cast<int>(graph.left_count);
    const int right = static_cast<int>(graph.right_count);
    const int states = 1 << right;
    const Cost kNone = std::numeric_limits<Cost>::min() / 4;
    std::vector<Cost> best(states, kNone);
    best[0] = 0;
    std::vector<std::vector<Cost>> weights(
        static_cast<std::size_t>(left),
        std::vector<Cost>(static_cast<std::size_t>(right), kNone));
    for (const Edge& edge : graph.edges) {
        weights[edge.left][edge.right] = edge.weight;
    }
    for (int i = 0; i < left; ++i) {
        std::vector<Cost> next = best;
        for (int mask = 0; mask < states; ++mask) {
            if (best[mask] == kNone) {
                continue;
            }
            for (int j = 0; j < right; ++j) {
                const Cost weight =
                    weights[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
                if ((mask & (1 << j)) != 0 || weight == kNone) {
                    continue;
                }
                const int next_mask = mask | (1 << j);
                next[next_mask] = std::max(next[next_mask], best[mask] + weight);
            }
        }
        best.swap(next);
    }
    Cost optimum = kNone;
    for (int mask = 1; mask < states; ++mask) {
        optimum = std::max(optimum, best[mask]);
    }
    return optimum;
}

// 题目期望答案：非空匹配的最大权重。
Cost expectedWeight(const Graph& graph) {
    if (graph.edges.empty()) {
        return 0;
    }
    Cost max_edge = graph.edges.front().weight;
    for (const Edge& edge : graph.edges) {
        max_edge = std::max(max_edge, edge.weight);
    }
    const Cost allow_empty = referenceAllowEmptyWeight(graph);
    // 允许空匹配的最大权为 0 时，非空最优就是权重最大的单条边（<= 0）。
    const Cost reference_non_empty = allow_empty > 0 ? allow_empty : max_edge;
    if (graph.left_count <= 8 && graph.right_count <= 8) {
        const Cost brute = bruteForceNonEmptyWeight(graph);
        if (reference_non_empty != brute) {
            std::cerr << "[自检] 参考实现与暴力解不一致: " << reference_non_empty
                      << " vs " << brute << '\n';
        }
        return brute;
    }
    return reference_non_empty;
}

// 校验返回匹配的合法性，成功时给出总权重。
std::string verify(const Graph& graph, const std::vector<std::int32_t>& answer,
                   Cost* weight) {
    if (answer.size() != graph.left_count) {
        return "返回数组长度错误";
    }
    std::vector<int> used_right(graph.right_count, 0);
    std::vector<std::vector<std::pair<std::int32_t, Cost>>> adjacency(
        graph.left_count);
    for (const Edge& edge : graph.edges) {
        adjacency[edge.left].emplace_back(
            static_cast<std::int32_t>(edge.right), edge.weight);
    }
    Cost total = 0;
    std::size_t matched = 0;
    for (std::size_t i = 0; i < answer.size(); ++i) {
        const std::int32_t matched_right = answer[i];
        if (matched_right == -1) {
            continue;
        }
        if (matched_right < 0 ||
            matched_right >= static_cast<std::int32_t>(graph.right_count)) {
            return "右部顶点编号越界";
        }
        if (used_right[static_cast<std::size_t>(matched_right)] != 0) {
            return "右部顶点重复使用";
        }
        used_right[static_cast<std::size_t>(matched_right)] = 1;
        bool found = false;
        for (const auto& entry : adjacency[i]) {
            if (entry.first == matched_right) {
                total += entry.second;
                found = true;
                break;
            }
        }
        if (!found) {
            return "返回了不存在的边";
        }
        ++matched;
    }
    if (!graph.edges.empty() && matched == 0) {
        return "有边时返回了空匹配";
    }
    if (graph.edges.empty() && matched != 0) {
        return "无边时应返回全 -1";
    }
    *weight = total;
    return "";
}

Graph makeRandom(std::mt19937_64* generator, std::uint32_t left,
                 std::uint32_t right, Cost low, Cost high, double density) {
    Graph graph;
    graph.left_count = left;
    graph.right_count = right;
    std::uniform_int_distribution<Cost> weight(low, high);
    std::vector<std::pair<std::uint32_t, std::uint32_t>> pairs;
    pairs.reserve(static_cast<std::size_t>(left) * right);
    for (std::uint32_t i = 0; i < left; ++i) {
        for (std::uint32_t j = 0; j < right; ++j) {
            pairs.emplace_back(i, j);
        }
    }
    if (density < 1.0) {
        std::shuffle(pairs.begin(), pairs.end(), *generator);
        pairs.resize(static_cast<std::size_t>(
            static_cast<double>(pairs.size()) * density));
    }
    for (const auto& pair : pairs) {
        graph.edges.push_back(Edge{pair.first, pair.second, weight(*generator)});
    }
    return graph;
}

struct TestCase {
    std::string name;
    Graph graph;
};

std::vector<TestCase> buildCases(std::uint32_t seed) {
    std::mt19937_64 generator(seed);
    std::vector<TestCase> cases;

    // ---- 基础档：边界用例 ----
    {
        Graph empty;
        cases.push_back({"基础-空图", empty});

        Graph no_edge;
        no_edge.left_count = 3;
        no_edge.right_count = 3;
        cases.push_back({"基础-无边", no_edge});

        Graph single;
        single.left_count = 1;
        single.right_count = 1;
        single.edges = {Edge{0, 0, -7}};
        cases.push_back({"基础-单条负边", single});

        Graph zero;
        zero.left_count = 2;
        zero.right_count = 2;
        zero.edges = {Edge{0, 0, 0}, Edge{1, 1, -1}};
        cases.push_back({"基础-零权边", zero});

        Graph sample;
        sample.left_count = 3;
        sample.right_count = 3;
        sample.edges = {Edge{0, 0, 5}, Edge{0, 1, 5}, Edge{1, 0, 5},
                        Edge{1, 1, 5}, Edge{2, 2, -7}};
        cases.push_back({"基础-题面样例", sample});

        Graph all_negative;
        all_negative.left_count = 2;
        all_negative.right_count = 2;
        all_negative.edges = {Edge{0, 0, -5}, Edge{1, 1, -3}};
        cases.push_back({"基础-全负边", all_negative});

        Graph single_best;
        single_best.left_count = 2;
        single_best.right_count = 2;
        single_best.edges = {Edge{0, 0, 100}, Edge{0, 1, 1}, Edge{1, 0, 1}};
        cases.push_back({"基础-单边更优", single_best});

        cases.push_back({"基础-随机1",
                         makeRandom(&generator, 8, 8, -20, 20, 0.5)});
        cases.push_back({"基础-随机2",
                         makeRandom(&generator, 6, 7, -100, 100, 1.0)});
        cases.push_back({"基础-随机3",
                         makeRandom(&generator, 8, 5, 0, 5, 0.7)});
    }

    // ---- 对拍：小规模随机图与暴力解逐一核对 ----
    for (std::size_t i = 0; i < 400; ++i) {
        const std::uint32_t left =
            1 + static_cast<std::uint32_t>(generator() % 8);
        const std::uint32_t right =
            1 + static_cast<std::uint32_t>(generator() % 8);
        const double density = 0.15 + static_cast<double>(generator() % 85) / 100.0;
        const Cost span = static_cast<Cost>(generator() % 30) + 1;
        cases.push_back({"对拍-随机" + std::to_string(i),
                         makeRandom(&generator, left, right, -span, span,
                                    density)});
    }

    // ---- 拓展档位 ----
    for (int i = 0; i < 5; ++i) {
        cases.push_back({"拓展1-" + std::to_string(i + 1),
                         makeRandom(&generator, 200, 200, -1000000, 1000000,
                                    1.0)});
    }
    for (int i = 0; i < 5; ++i) {
        cases.push_back({"拓展2-" + std::to_string(i + 1),
                         makeRandom(&generator, 500, 500, -1000000000,
                                    1000000000, 1.0)});
    }
    for (int i = 0; i < 5; ++i) {
        Graph graph;
        graph.left_count = 2000;
        graph.right_count = 2000;
        std::uniform_int_distribution<std::uint32_t> left_index(0, 1999);
        std::uniform_int_distribution<std::uint32_t> right_index(0, 1999);
        std::uniform_int_distribution<Cost> weight(-1000000000, 1000000000);
        std::set<std::uint64_t> seen;
        graph.edges.reserve(50000);
        while (seen.size() < 50000) {
            const std::uint64_t value =
                (static_cast<std::uint64_t>(left_index(generator)) << 32) |
                right_index(generator);
            if (seen.insert(value).second) {
                graph.edges.push_back(
                    Edge{static_cast<std::uint32_t>(value >> 32),
                         static_cast<std::uint32_t>(value & 0xffffffffULL),
                         weight(generator)});
            }
        }
        cases.push_back({"拓展3-" + std::to_string(i + 1), graph});
    }
    return cases;
}

}  // namespace

int main(int argc, char** argv) {
    std::uint32_t seed = 20271001U;
    if (argc > 1) {
        seed = static_cast<std::uint32_t>(std::strtoul(argv[1], nullptr, 10));
    }
    const auto cases = buildCases(seed);
    std::cout << "PIP2027 / T1 二分图最大权匹配 本地测评\n"
              << "用例数: " << cases.size() << "（含 400 组小规模对拍）\n\n";
    std::cout << std::right << std::setw(20) << "Case" << "  " << std::left
              << std::setw(16) << "Verdict" << std::right << std::setw(12)
              << "Solve(ms)" << std::setw(20) << "Weight" << '\n'
              << std::string(72, '-') << '\n'
              << std::fixed << std::setprecision(3);

    std::size_t passed = 0;
    for (const TestCase& test : cases) {
        Cost weight = 0;
        std::string verdict = "PASS";
        double solve_ms = 0.0;
        Cost expected = 0;
        try {
            const auto start = std::chrono::steady_clock::now();
            const BipartiteMatcher matcher(test.graph.left_count,
                                           test.graph.right_count,
                                           test.graph.edges);
            const std::vector<std::int32_t> answer =
                matcher.maximumWeightMatching();
            const auto stop = std::chrono::steady_clock::now();
            solve_ms =
                std::chrono::duration<double, std::milli>(stop - start).count();

            const std::string problem = verify(test.graph, answer, &weight);
            if (!problem.empty()) {
                verdict = problem;
            } else {
                expected = expectedWeight(test.graph);
                if (!test.graph.edges.empty() && weight != expected) {
                    verdict = "WRONG_ANSWER";
                }
            }
        } catch (const std::exception& error) {
            verdict = std::string("EXCEPTION:") + error.what();
        }
        if (verdict == "PASS") {
            ++passed;
        }
        std::cout << std::right << std::setw(20) << test.name << "  "
                  << std::left << std::setw(16) << verdict << std::right
                  << std::setw(12) << solve_ms << std::setw(20) << weight;
        if (verdict == "WRONG_ANSWER") {
            std::cout << "  期望 " << expected;
        }
        std::cout << '\n';
    }
    std::cout << std::string(72, '-') << "\n结果: " << passed << '/'
              << cases.size() << '\n';
    return passed == cases.size() ? 0 : 1;
}
