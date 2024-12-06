#pragma once
//#include "stdafx.h"
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

template <typename T>
using MatrixType = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template <typename T> using VectorType = Eigen::Matrix<T, Eigen::Dynamic, 1>;

template <typename Scalar, int rank, typename sizeType>
auto Tensor_to_Matrix(const Eigen::Tensor<Scalar, rank> &tensor,
                      const sizeType rows, const sizeType cols) {
  return Eigen::Map<const MatrixType<Scalar>>(tensor.data(), rows, cols);
}

template <typename Scalar>
auto Tensor_to_Matrix(const Eigen::Tensor<Scalar, 2> &tensor) {
  const auto &d = tensor.dimensions();
  return Eigen::Map<const MatrixType<Scalar>>(tensor.data(), d[0], d[1]);
}

template <typename Scalar, typename... Dims>
auto Matrix_to_Tensor(const MatrixType<Scalar> &matrix, Dims... dims) {
  constexpr int rank = sizeof...(Dims);
  return Eigen::TensorMap<Eigen::Tensor<const Scalar, rank>>(matrix.data(),
                                                             {dims...});
}

template <typename Scalar>
auto Matrix_to_Tensor(const MatrixType<Scalar> &matrix) {
  return Eigen::TensorMap<Eigen::Tensor<const Scalar, 2>>(
      matrix.data(), matrix.rows(), matrix.cols());
}

template <typename Scalar>
auto Vector_to_Tensor(Eigen::Matrix<Scalar, -1, 1> &vec) {
  return Eigen::TensorMap<Eigen::Tensor<Scalar, 1>>(vec.data(), vec.size());
}

template <typename Scalar>
auto Tensor_to_Vector(const Eigen::Tensor<Scalar, 1> &tensor) {
  const auto &d = tensor.dimensions();
  return Eigen::Map<const VectorType<Scalar>>(tensor.data(), d[0], 1);
}

template <typename Scalar, int rankA, int rankB>
auto EigenTensorProduct(const Eigen::Tensor<Scalar, rankA> &tensorA,
                        const Eigen::Tensor<Scalar, rankB> &tensorB,
                        const int reduction_index_A,
                        const int reduction_index_B) {
  Eigen::array<Eigen::IndexPair<int>, 1> product_dims = {
      Eigen::IndexPair<int>(reduction_index_A, reduction_index_B)};
  return tensorA.contract(tensorB, product_dims);
}

