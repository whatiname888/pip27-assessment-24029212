#include "dynamic_kth.hpp"

#include "wavelet_matrix.hpp"

#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

// T2：动态序列区间第 k 大。
//
// 核心结构是 8 层 4 位数字的动态小波矩阵（wavelet_matrix.hpp）：
// 每层保存序列各元素的对应 4 位数字并按其稳定划分，区间第 k 大逐层下探
// 确定答案的每一位，单次查询 O(8·log n)，插入删除同阶。详细设计见
// wavelet_matrix.hpp 头部注释。
struct DynamicKth::Impl {
    explicit Impl(const std::vector<int>& initial) : matrix(initial) {}
    pip27::WaveletMatrix matrix;
};

DynamicKth::DynamicKth(const std::vector<int>& initial)
    : impl_(std::make_unique<Impl>(initial)) {}

DynamicKth::~DynamicKth() = default;

void DynamicKth::insert(std::size_t index, int value) {
    impl_->matrix.insert(index, value);
}

void DynamicKth::erase(std::size_t index) { impl_->matrix.erase(index); }

void DynamicKth::update(std::size_t index, int value) {
    impl_->matrix.replace(index, value);
}

int DynamicKth::kthLargest(std::size_t left, std::size_t right,
                           std::size_t k) const {
    return impl_->matrix.kthLargest(left, right, k);
}
