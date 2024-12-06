#pragma once
#define outputprecision 64
#include <fstream>
#include <iostream>

namespace EigenMatrixIO{
	template<class Matrix>
	void write_binary(const char* filename, const Matrix& matrix){
		std::ofstream out(filename, std::ios::out | std::ios::binary | std::ios::trunc);
		out.precision(outputprecision);
		typename Matrix::Index rows = matrix.rows(), cols = matrix.cols();
		out.write((char*)(&rows), sizeof(typename Matrix::Index));
		out.write((char*)(&cols), sizeof(typename Matrix::Index));
		out.write((char*)matrix.data(), rows*cols*sizeof(typename Matrix::Scalar));
		out.close();
	};

	template<class Matrix>
	void write_binary(std::ofstream & out, const Matrix& matrix){
		typename Matrix::Index rows = matrix.rows(), cols = matrix.cols();
		out.write((char*)(&rows), sizeof(typename Matrix::Index));
		out.write((char*)(&cols), sizeof(typename Matrix::Index));
		out.write((char*)matrix.data(), rows*cols*sizeof(typename Matrix::Scalar));
	};

	template<class Matrix>
	bool read_binary(const char* filename, Matrix& matrix){
		std::ifstream in(filename, std::ios::in | std::ios::binary);
		if (!in.good())
		{
			std::cout << "file not open" << std::endl;
			return false;
		}
		typename Matrix::Index rows = 0, cols = 0;
		in.read((char*)(&rows), sizeof(typename Matrix::Index));
		in.read((char*)(&cols), sizeof(typename Matrix::Index));
		matrix.resize(rows, cols);
		in.read((char *)matrix.data(), rows*cols*sizeof(typename Matrix::Scalar));
		in.close();
		return true;
	};

	template <class Matrix>
	void write_vector_eigen(std::ofstream &out,
													std::vector<Matrix> vec) {
		int size_of_vec = vec.size();
		out.write((char *)&size_of_vec, sizeof(int));
		for (auto v : vec) {
			typename Matrix::Index rows = v.rows(), cols = v.cols();
			out.write((char *)(&rows), sizeof(typename Matrix::Index));
			out.write((char *)(&cols), sizeof(typename Matrix::Index));
			out.write((char *)v.data(),
								rows * cols * sizeof(typename Matrix::Scalar));
		}
	}

	template <class Matrix>
  void write_vector_eigen_pure_data(std::ofstream &out, std::vector<Matrix> vec) {
    int size_of_vec = vec.size();
    //out.write((char *)&size_of_vec, sizeof(int));
    for (auto v : vec) {
      typename Matrix::Index rows = v.rows(), cols = v.cols();
      //out.write((char *)(&rows), sizeof(typename Matrix::Index));
      //out.write((char *)(&cols), sizeof(typename Matrix::Index));
      out.write((char *)v.data(),
                rows * cols * sizeof(typename Matrix::Scalar));
    }
  }


	template <class Matrix>
	bool read_vector_eigen(const char *filename, std::vector<Matrix> &vec) {
		std::ifstream in(filename, std::ios::in | std::ios::binary);
		if (!in.good()) {
			std::cout << "file not open" << std::endl;
			return false;
		}
		vec.clear();
		int size_of_vec = 0;
		in.read((char *)&size_of_vec, sizeof(int));

		for (int i = 0; i < size_of_vec; ++i) {
			Matrix matrix;
			typename Matrix::Index rows = 0, cols = 0;
			in.read((char *)(&rows), sizeof(typename Matrix::Index));
			in.read((char *)(&cols), sizeof(typename Matrix::Index));
			matrix.resize(rows, cols);
			in.read((char *)matrix.data(),
							rows * cols * sizeof(typename Matrix::Scalar));
			vec.push_back(matrix);
		}
		in.close();
		return true;
	};

	template<class Matrix>
	bool read_ascii(const char* filename, Matrix& matrix, int nrows = 2070, int ncols = 30) {
      //int nrows = 3;
      //int ncols = 4;
      std::ifstream fin(filename);

      if (fin.is_open()) {
        fin >> nrows >> ncols;
        matrix.resize(nrows, ncols);
        for (int row = 0; row < nrows; row++)
          for (int col = 0; col < ncols; col++) {
            double item = 0.0;
            fin >> item;
            matrix(row, col) = item;
          }
        fin.close();
      }
			return true;
	}

	template<class Matrix>
	void read_binary(std::ifstream &in, Matrix& matrix){
		/*if (!in.good())
		{
			std::cout << "file not open" << std::endl;
			return;
		}*/
		typename Matrix::Index rows = 0, cols = 0;
		in.read((char*)(&rows), sizeof(typename Matrix::Index));
		in.read((char*)(&cols), sizeof(typename Matrix::Index));
		matrix.resize(rows, cols);
		in.read((char *)matrix.data(), rows*cols*sizeof(typename Matrix::Scalar));
	};
};