#pragma once
//#include "AutoDiff/AutoDiffCore.h"
#include "LoboDynamic/LoboVolumtricMesh/LoboTetMesh.h"
#include "TypeMaterialWithCompressionResistance.h"
#include <vector>
#include <iostream>
// adapter between TypeMaterialWithCompressionResistance and material based on
// the tet mesh

template <class TYPE>
class TypeTetElementMaterial
    : public TypeMaterialWithCompressionResistance<TYPE> {
  typedef Eigen::Matrix<TYPE, 3, 3> Matrix3;
  typedef Eigen::Matrix<TYPE, -1, -1> MatrixX;

public:
  //LOBO_MAKE_TYPEDEFS(TYPE, t);

  TypeTetElementMaterial(Lobo::LoboTetMesh *tetmesh,
                         int enableCompressionResistance = 0,
                         TYPE compressionResistance = 0.0);
  ~TypeTetElementMaterial();

  virtual void updateMaterial();

  virtual TYPE ComputeEnergy(int elementIndex, TYPE *invariants) = 0;
  virtual void ComputeEnergyGradient(
      int elementIndex, TYPE *invariants,
      TYPE *gradient) = 0; // invariants and gradient are 3-vectors
  virtual void ComputeEnergyHessian(
      int elementIndex, TYPE *invariants,
      TYPE *hessian) = 0; // invariants is a 3-vector, hessian is a 3x3
                          // symmetric matrix, unrolled into a 6-vector, in the
                          // following order: (11, 12, 13, 22, 23, 33).

  virtual void getLambdaLame(int eleid, TYPE &lambda);
  virtual void getMuLame(int eleid, TYPE &mu);

protected:
  TYPE *lambdaLame;
  TYPE *muLame;

  TYPE compressionResistance;
  TYPE *EdivNuFactor;
  virtual TYPE GetCompressionResistanceFactor(int elementIndex);

  Lobo::LoboTetMesh *tetmesh;

};

template <class TYPE>
void TypeTetElementMaterial<TYPE>::getMuLame(int eleid, TYPE &mu) {
  mu = muLame[eleid];
}

template <class TYPE>
void TypeTetElementMaterial<TYPE>::getLambdaLame(int eleid, TYPE &lambda) {
  lambda = lambdaLame[eleid];
}
