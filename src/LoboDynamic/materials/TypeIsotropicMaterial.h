#pragma once
#include "Functions/LoboMacros.h"
#include <cmath>
#include <Eigen/Dense>

//#include "AutoDiff/AutoDiffCore.h"

template<class TYPE>
class TypeIsotropicMaterial
{
public:

	typedef Eigen::Matrix<TYPE, 3, 3> Matrix3;
	typedef Eigen::Matrix<TYPE, -1, -1> MatrixX;

	TypeIsotropicMaterial();
	~TypeIsotropicMaterial();

	virtual void updateMaterial() = 0;

	virtual TYPE ComputeEnergy(int elementIndex, TYPE * invariants) = 0;

	virtual void ComputeEnergyGradient(int elementIndex, TYPE * invariants, TYPE * gradient) = 0; // invariants and gradient are 3-vectors

	virtual void ComputeEnergyHessian(int elementIndex, TYPE * invariants, TYPE * hessian) = 0; // invariants is a 3-vector, hessian is a 3x3 symmetric matrix, unrolled into a 6-vector, in the following order: (11, 12, 13, 22, 23, 33).

	virtual void getLambdaLame(int eleid, TYPE &lambda) {};
	virtual void getMuLame(int eleid, TYPE &mu) {};

	double h_CSFD;

};

template <class TYPE> TypeIsotropicMaterial<TYPE>::TypeIsotropicMaterial() {
  //h_CSFD = lobo_h;
}


template<class TYPE>
TypeIsotropicMaterial<TYPE>::~TypeIsotropicMaterial()
{

}


