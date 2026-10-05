#include "dynamic_kth.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <memory>
#include <unordered_set>
#include <utility>
#include <vector>

// T2 本地测评程序（不参与正式测评）：
//   * 小规模、集中值域与拓展 1 用 std::vector + nth_element 暴力逐查询核对；
//   * 拓展 2（N=5e6）同步维护暴力序列做全量操作，按间隔抽验查询结果
//     （暴力单次查询 O(N)，全量核对代价过高），并单独统计实现耗时；
//   * “集中值域”用例把取值压到很窄的区间，专门检验候选退化兜底路径。
// 用法：./grader [seed]
namespace {

class BruteSequence {
public:
    explicit BruteSequence(std::vector<int> initial)
        : data_(std::move(initial)) {}

    void insert(std::size_t index, int value) {
        data_.insert(data_.begin() + static_cast<std::ptrdiff_t>(index), value);
    }

    void erase(std::size_t index) {
        data_.erase(data_.begin() + static_cast<std::ptrdiff_t>(index));
    }

    void update(std::size_t index, int value) { data_[index] = value; }

    int at(std::size_t index) const { return data_[index]; }

    int kthLargest(std::size_t left, std::size_t right, std::size_t k) const {
        std::vector<int> range(
            data_.begin() + static_cast<std::ptrdiff_t>(left),
            data_.begin() + static_cast<std::ptrdiff_t>(right) + 1);
        const std::size_t index = k - 1;
        std::nth_element(range.begin(),
                         range.begin() + static_cast<std::ptrdiff_t>(index),
                         range.end(), std::greater<int>());
        return range[index];
    }

    std::size_t size() const { return data_.size(); }

private:
    std::vector<int> data_;
};

struct Config {
    std::string name;
    std::size_t initial_size;
    std::size_t operations;
    int value_span;               // 取值范围 [-value_span, value_span]
    std::size_t verify_stride;    // 抽验间隔，1 表示每个查询都核对
};

struct CaseResult {
    bool passed = true;
    std::string message;
    double impl_ms = 0.0;
    std::size_t verified_queries = 0;
};

std::vector<int> makeValues(std::mt19937_64* generator, std::size_t count,
                            int span, std::unordered_set<int>* used) {
    std::uniform_int_distribution<int> distribution(-span, span);
    std::vector<int> values;
    values.reserve(count);
    while (values.size() < count) {
        const int value = distribution(*generator);
        if (used->insert(value).second) {
            values.push_back(value);
        }
    }
    return values;
}

CaseResult runCase(const Config& config, std::mt19937_64* generator) {
    CaseResult result;
    std::uniform_int_distribution<std::uint32_t> roll(0, 99);
    std::unordered_set<int> used;
    used.reserve(config.initial_size + config.operations / 2 + 16);

    const std::vector<int> initial =
        makeValues(generator, config.initial_size, config.value_span, &used);
    std::uniform_int_distribution<int> value_distribution(
        -config.value_span, config.value_span);

    double impl_ms = 0.0;
    BruteSequence brute(initial);
    std::unique_ptr<DynamicKth> sequence;
    {
        const auto start = std::chrono::steady_clock::now();
        sequence = std::make_unique<DynamicKth>(initial);
        const auto stop = std::chrono::steady_clock::now();
        impl_ms +=
            std::chrono::duration<double, std::milli>(stop - start).count();
    }

    std::size_t query_seen = 0;
    for (std::size_t op = 0; op < config.operations; ++op) {
        const std::size_t length = brute.size();
        const std::uint32_t action = roll(*generator);
        if (length > 0 && action < 25) {
            // ---- 查询 ----
            std::uniform_int_distribution<std::size_t> pick_left(0, length - 1);
            const std::size_t left = pick_left(*generator);
            std::uniform_int_distribution<std::size_t> pick_right(left,
                                                                 length - 1);
            const std::size_t right = pick_right(*generator);
            std::uniform_int_distribution<std::size_t> pick_rank(
                1, right - left + 1);
            const std::size_t k = pick_rank(*generator);

            const auto start = std::chrono::steady_clock::now();
            const int answer = sequence->kthLargest(left, right, k);
            const auto stop = std::chrono::steady_clock::now();
            impl_ms +=
                std::chrono::duration<double, std::milli>(stop - start).count();

            if (query_seen % config.verify_stride == 0) {
                const int expected = brute.kthLargest(left, right, k);
                ++result.verified_queries;
                if (answer != expected) {
                    result.passed = false;
                    result.message =
                        "查询错误 op=" + std::to_string(op) +
                        " left=" + std::to_string(left) +
                        " right=" + std::to_string(right) +
                        " k=" + std::to_string(k) + " 得到 " +
                        std::to_string(answer) + " 期望 " +
                        std::to_string(expected);
                    result.impl_ms = impl_ms;
                    return result;
                }
            }
            ++query_seen;
        } else if (action < 50) {
            // ---- 插入 ----
            std::uniform_int_distribution<std::size_t> pick_position(0, length);
            const std::size_t index = pick_position(*generator);
            int value = value_distribution(*generator);
            while (!used.insert(value).second) {
                value = value_distribution(*generator);
            }
            const auto start = std::chrono::steady_clock::now();
            sequence->insert(index, value);
            const auto stop = std::chrono::steady_clock::now();
            impl_ms +=
                std::chrono::duration<double, std::milli>(stop - start).count();
            brute.insert(index, value);
        } else if (length > 0 && action < 75) {
            // ---- 删除 ----
            std::uniform_int_distribution<std::size_t> pick_position(0,
                                                                    length - 1);
            const std::size_t index = pick_position(*generator);
            used.erase(brute.at(index));
            const auto start = std::chrono::steady_clock::now();
            sequence->erase(index);
            const auto stop = std::chrono::steady_clock::now();
            impl_ms +=
                std::chrono::duration<double, std::milli>(stop - start).count();
            brute.erase(index);
        } else if (length > 0) {
            // ---- 修改 ----
            std::uniform_int_distribution<std::size_t> pick_position(0,
                                                                    length - 1);
            const std::size_t index = pick_position(*generator);
            used.erase(brute.at(index));
            int value = value_distribution(*generator);
            while (!used.insert(value).second) {
                value = value_distribution(*generator);
            }
            const auto start = std::chrono::steady_clock::now();
            sequence->update(index, value);
            const auto stop = std::chrono::steady_clock::now();
            impl_ms +=
                std::chrono::duration<double, std::milli>(stop - start).count();
            brute.update(index, value);
        }
    }
    result.impl_ms = impl_ms;
    return result;
}

std::vector<Config> buildConfigs() {
    std::vector<Config> configs;
    for (int i = 0; i < 8; ++i) {
        configs.push_back({"对拍-小规模" + std::to_string(i + 1), 30, 200,
                           40, 1});
    }
    for (int i = 0; i < 8; ++i) {
        configs.push_back({"基础-混合" + std::to_string(i + 1), 300, 1000,
                           1000, 1});
    }
    configs.push_back({"基础-空初始", 0, 1000, 1000, 1});
    configs.push_back({"基础-单元素", 1, 1000, 1000, 1});
    for (int i = 0; i < 3; ++i) {
        configs.push_back({"集中值域-" + std::to_string(i + 1), 20000, 3000,
                           30000, 1});
    }
    for (int i = 0; i < 5; ++i) {
        configs.push_back({"拓展1-" + std::to_string(i + 1), 100000, 40000,
                           1000000000, 1});
    }
    for (int i = 0; i < 5; ++i) {
        configs.push_back({"拓展2-" + std::to_string(i + 1), 5000000, 20000,
                           1000000000, 100});
    }
    return configs;
}

}  // namespace

int main(int argc, char** argv) {
    std::uint32_t seed = 20271001U;
    if (argc > 1) {
        seed = static_cast<std::uint32_t>(std::strtoul(argv[1], nullptr, 10));
    }
    std::mt19937_64 generator(seed);
    auto configs = buildConfigs();
    if (argc > 2) {
        configs.resize(std::min(configs.size(),
                                static_cast<std::size_t>(std::strtoul(argv[2], nullptr, 10))));
    }

    std::cout << "PIP2027 / T2 动态序列区间第 k 大 本地测评\n"
              << "用例数: " << configs.size() << "\n\n";
    std::cout << std::right << std::setw(18) << "Case" << "  " << std::left
              << std::setw(14) << "Verdict" << std::right << std::setw(12)
              << "Impl(ms)" << std::setw(14) << "Verified" << '\n'
              << std::string(62, '-') << '\n'
              << std::fixed << std::setprecision(3);

    std::size_t passed = 0;
    for (const Config& config : configs) {
        CaseResult result = runCase(config, &generator);
        if (result.passed) {
            ++passed;
        }
        std::cout << std::right << std::setw(18) << config.name << "  "
                  << std::left << std::setw(14)
                  << (result.passed ? "PASS" : "WRONG_ANSWER") << std::right
                  << std::setw(12) << result.impl_ms << std::setw(14)
                  << result.verified_queries << '\n';
        if (!result.message.empty()) {
            std::cout << "    " << result.message << '\n';
        }
    }
    std::cout << std::string(62, '-') << "\n结果: " << passed << '/'
              << configs.size() << '\n';
    return passed == configs.size() ? 0 : 1;
}
