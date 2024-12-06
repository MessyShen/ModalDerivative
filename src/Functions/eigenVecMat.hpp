#pragma once
#include <Eigen/Dense>
namespace Lobo
{
    template <typename T>
    inline Eigen::Matrix<T, -1, -1> eigen_vec_2_mat(const Eigen::Matrix<T, -1, 1> &inputmat, int row, int col)
    {
        Eigen::Matrix<T, -1, -1> mat(row, col);
        for (int i = 0; i < row; ++i)
        {
            for (int j = 0; j < col; ++j)
            {
                mat.data()[j * row + i] = inputmat.data()[i * col + j];
            }
        }
        return mat;
    }

    template <typename T>
    inline Eigen::Matrix<T, -1, 1> eigen_mat_2_vec(const Eigen::Matrix<T, -1, -1> &inputmat)
    {
        int row = inputmat.rows();
        int col = inputmat.cols();
        Eigen::Matrix<T, -1, 1> vec(row * col);
        for (int i = 0; i < row; ++i)
        {
            for (int j = 0; j < col; ++j)
            {
                vec.data()[i * col + j] = inputmat.data()[j * row + i];
            }
        }
        return vec;
    }
}