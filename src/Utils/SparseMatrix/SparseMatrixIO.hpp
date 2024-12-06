#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Sparse>

void saveSparseMatrix(const Eigen::SparseMatrix<double>& mat, const std::string& filename) {
    std::ofstream outFile(filename, std::ios::binary);
    if (outFile.is_open()) {
        // Write the matrix dimensions and number of non-zero elements
        int rows = mat.rows();
        int cols = mat.cols();
        int nonZeros = mat.nonZeros();
        outFile.write(reinterpret_cast<const char*>(&rows), sizeof(int));
        outFile.write(reinterpret_cast<const char*>(&cols), sizeof(int));
        outFile.write(reinterpret_cast<const char*>(&nonZeros), sizeof(int));

        // Write the non-zero elements
        for (int k = 0; k < mat.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(mat, k); it; ++it) {
                int row = it.row();
                int col = it.col();
                double value = it.value();
                outFile.write(reinterpret_cast<const char*>(&row), sizeof(int));
                outFile.write(reinterpret_cast<const char*>(&col), sizeof(int));
                outFile.write(reinterpret_cast<const char*>(&value), sizeof(double));
            }
        }
        outFile.close();
    } else {
        std::cerr << "Unable to open file for writing: " << filename << std::endl;
    }
}

Eigen::SparseMatrix<double> loadSparseMatrix(const std::string& filename) {
    std::ifstream inFile(filename, std::ios::binary);
    if (!inFile.is_open()) {
        std::cerr << "Unable to open file for reading: " << filename << std::endl;
        exit(EXIT_FAILURE);
    }

    int rows, cols, nonZeros;
    inFile.read(reinterpret_cast<char*>(&rows), sizeof(int));
    inFile.read(reinterpret_cast<char*>(&cols), sizeof(int));
    inFile.read(reinterpret_cast<char*>(&nonZeros), sizeof(int));

    Eigen::SparseMatrix<double> mat(rows, cols);
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(nonZeros);

    for (int i = 0; i < nonZeros; ++i) {
        int row, col;
        double value;
        inFile.read(reinterpret_cast<char*>(&row), sizeof(int));
        inFile.read(reinterpret_cast<char*>(&col), sizeof(int));
        inFile.read(reinterpret_cast<char*>(&value), sizeof(double));
        tripletList.emplace_back(row, col, value);
    }

    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    inFile.close();
    return mat;
}
