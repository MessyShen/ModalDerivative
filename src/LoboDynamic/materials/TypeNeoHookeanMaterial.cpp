#include "TypeNeoHookeanMaterial.h"
#include <complex>
// #include <omp.h>
#include "Functions/LoboMacros.h"


template <class TYPE>
TypeNeoHookeanMaterial<TYPE>::TypeNeoHookeanMaterial(Lobo::LoboTetMesh* tetmesh, int enableCompressionResistance /*= 0*/, TYPE compressionResistance /*= 0.0*/):TypeTetElementMaterial<TYPE>(tetmesh, enableCompressionResistance, compressionResistance)
{

}

template <class TYPE>
TypeNeoHookeanMaterial<TYPE>::~TypeNeoHookeanMaterial()
{
	
}


#define ENERGY_FUNCTION(LoboComplex_NAME)\
E_complex = F_complex.transpose()*F_complex;\
LoboComplex_NAME I1 = E_complex.trace();\
LoboComplex_NAME J = F_complex.det();\
LoboComplex_NAME A = lobo::pow(J, (TYPE)(-2.0 / 3.0));\
LoboComplex_NAME I1_ = A*I1;\
LoboComplex_NAME J_1_2 = J - (TYPE)1.0;\
LoboComplex_NAME e_j = J_1_2*J_1_2;\
energy = muLame[elementIndex] * (I1_ - (TYPE)3.0) + lambdaLame[elementIndex] * e_j;

//#define ENERGY_FUNCTION(LoboComplex_NAME)\
//E_complex = F_complex.transpose()*F_complex;\
//LoboComplex_NAME det_F = F_complex.det();\
//LoboComplex_NAME log_det_F = lobo::log(det_F);\
//E_complex = F_complex.transpose()*F_complex;\
//LoboComplex_NAME I1 = E_complex.trace();\
//energy = muLame[elementIndex] * (TYPE)0.5*(I1 - (TYPE)3.0);\
//energy -= muLame[elementIndex] * log_det_F;\
//energy += lambdaLame[elementIndex] * (TYPE)0.5 * log_det_F*log_det_F;



template <class TYPE>
TYPE TypeNeoHookeanMaterial<TYPE>::ComputeEnergy(int elementIndex, TYPE * invariants)
{
	TYPE IC = invariants[0];
	TYPE IIIC = invariants[2];
	TYPE J = std::sqrt(IIIC);
	TYPE logJ = std::log(J);


	TYPE energy = 0.5 * muLame[elementIndex] * (IC - 3.0) - muLame[elementIndex] * logJ + 0.5 * lambdaLame[elementIndex] * logJ * logJ;

	this->AddCompressionResistanceEnergy(elementIndex, invariants, &energy);
    /*std::cout << "neohookean energy" << elementIndex  << " " << energy << " "
                  << IC << " " << IIIC
                  << " " << logJ << " " << muLame[elementIndex] << " "
                  << std::endl;*/
	return energy;

}

template <class TYPE>
void TypeNeoHookeanMaterial<TYPE>::ComputeEnergyGradient(int elementIndex, TYPE * invariants, TYPE * gradient)
{
	TYPE IIIC = invariants[2];
	gradient[0] = 0.5 * muLame[elementIndex];
	gradient[1] = 0.0;
	gradient[2] = (-0.5 * muLame[elementIndex] + 0.25 * lambdaLame[elementIndex] * std::log(IIIC)) / IIIC;

	this->AddCompressionResistanceGradient(elementIndex, invariants, gradient);
}

template <class TYPE>
void TypeNeoHookeanMaterial<TYPE>::ComputeEnergyHessian(int elementIndex, TYPE * invariants, TYPE * hessian)
{
	TYPE IIIC = invariants[2];
	// 11
	hessian[0] = 0.0;
	// 12
	hessian[1] = 0.0;
	// 13
	hessian[2] = 0.0;
	// 22
	hessian[3] = 0.0;
	// 23
	hessian[4] = 0.0;
	// 33
	hessian[5] = (0.25 * lambdaLame[elementIndex] + 0.5 * muLame[elementIndex] - 0.25 * lambdaLame[elementIndex] * std::log(IIIC)) / (IIIC * IIIC);

	this->AddCompressionResistanceHessian(elementIndex, invariants, hessian);
}



LOBO_TEMPLATE_INSTANT(TypeNeoHookeanMaterial)