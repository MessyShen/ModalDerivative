#include "TypeStVKMaterial.h"
#include <complex>
// #include <omp.h>
#include "Functions/LoboMacros.h"
//#include "MCSFD/MatrixOp.h"


template <class TYPE>
TypeStVKMaterial<TYPE>::TypeStVKMaterial(Lobo::LoboTetMesh *tetmesh, int enableCompressionResistance /*= 0*/, TYPE compressionResistance /*= 0.0*/) : TypeTetElementMaterial<TYPE>(tetmesh, enableCompressionResistance, compressionResistance)
{
}

template <class TYPE>
TypeStVKMaterial<TYPE>::~TypeStVKMaterial()
{
}

#define ENERGY_FUNCTION(LoboComplex_NAME)                                            \
	E_complex = F_complex.transpose() * F_complex;                                   \
	E_complex.data[0] -= 1;                                                          \
	E_complex.data[4] -= 1;                                                          \
	E_complex.data[8] -= 1;                                                          \
	E_complex = E_complex * (TYPE)0.5;                                               \
	energy = (TYPE)muLame[elementIndex] * lobo::inner_product(E_complex, E_complex); \
	LoboComplex_NAME E_trace = E_complex.trace();                                    \
	E_trace = (TYPE)0.5 * lambdaLame[elementIndex] * E_trace * E_trace;              \
	energy += E_trace;

template <class TYPE>
TYPE TypeStVKMaterial<TYPE>::ComputeEnergy(int elementIndex, TYPE *invariants)
{
	TYPE IC = invariants[0];
	TYPE IIC = invariants[1];
	//double IIIC = invariants[2]; // not needed for StVK

	TYPE energy = 0.125 * lambdaLame[elementIndex] * (IC - 3.0) * (IC - 3.0) + 0.25 * muLame[elementIndex] * (IIC - 2.0 * IC + 3.0);

	this->AddCompressionResistanceEnergy(elementIndex, invariants, &energy);

	return energy;
}

template <class TYPE>
void TypeStVKMaterial<TYPE>::ComputeEnergyGradient(int elementIndex, TYPE *invariants, TYPE *gradient)
{
	TYPE IC = invariants[0];
	gradient[0] = 0.25 * lambdaLame[elementIndex] * (IC - 3.0) - 0.5 * muLame[elementIndex];
	gradient[1] = 0.25 * muLame[elementIndex];
	gradient[2] = 0.0;

	this->AddCompressionResistanceGradient(elementIndex, invariants, gradient);
}

template <class TYPE>
void TypeStVKMaterial<TYPE>::ComputeEnergyHessian(int elementIndex, TYPE *invariants, TYPE *hessian)
{
	hessian[0] = 0.25 * lambdaLame[elementIndex];
	// 12
	hessian[1] = 0.0;
	// 13
	hessian[2] = 0.0;
	// 22
	hessian[3] = 0.0;
	// 23
	hessian[4] = 0.0;
	// 33
	hessian[5] = 0.0;

	this->AddCompressionResistanceHessian(elementIndex, invariants, hessian);
}



LOBO_TEMPLATE_INSTANT(TypeStVKMaterial)