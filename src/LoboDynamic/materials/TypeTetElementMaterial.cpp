#include "TypeTetElementMaterial.h"


template <class TYPE>
TypeTetElementMaterial<TYPE>::TypeTetElementMaterial(
    Lobo::LoboTetMesh *tetmesh, int enableCompressionResistance /*= 0*/,
    TYPE compressionResistance /*= 0.0*/) {
  this->enableCompressionResistance = enableCompressionResistance;
  this->tetmesh = tetmesh;
  int numElements = tetmesh->getNumElements();
  lambdaLame = (TYPE *)malloc(sizeof(TYPE) * numElements);
  muLame = (TYPE *)malloc(sizeof(TYPE) * numElements);

  if (enableCompressionResistance)
    EdivNuFactor = (TYPE *)malloc(sizeof(TYPE) * numElements);
  else
    EdivNuFactor = NULL;

#pragma omp parallel for
  for (int el = 0; el < numElements; el++) {
    Lobo::Material *material = tetmesh->getElementMaterial(el);

    if (material == NULL) {
      printf("Error: NeoHookeanIsotropicMaterial: mesh does not consist of E, "
             "nu materials.\n");
      throw 1;
    }

    lambdaLame[el] = (TYPE)material->getLambda();
    muLame[el] = (TYPE)material->getMu();

    if (enableCompressionResistance) {
      EdivNuFactor[el] = compressionResistance * (TYPE)material->getE() /
                         (1.0 - 2.0 * (TYPE)material->getNu());
      // printf("Setting EdivNuFactor[%d]=%G\n", el, EdivNuFactor[el]);
    }
  }

}

template <class TYPE> TypeTetElementMaterial<TYPE>::~TypeTetElementMaterial() {
  free(EdivNuFactor);
  free(lambdaLame);
  free(muLame);
}

template <class TYPE> void TypeTetElementMaterial<TYPE>::updateMaterial() {
  std::cout << "updateMaterial  " << std::endl;
  int numElements = tetmesh->getNumElements();

  for (int el = 0; el < numElements; el++) {
    Lobo::Material *material = tetmesh->getElementMaterial(el);
    // LoboVolumetricMesh::Material * material =
    // tetmesh->getElementMaterial(el); LoboVolumetricMesh::ENuMaterial *
    // eNuMaterial = downcastENuMaterial(material);
    if (material == NULL) {
      printf("Error: NeoHookeanIsotropicMaterial: mesh does not consist of E, "
             "nu materials.\n");
      throw 1;
    }

    lambdaLame[el] = (TYPE)material->getLambda();
    muLame[el] = (TYPE)material->getMu();

    if (this->enableCompressionResistance) {
      EdivNuFactor[el] = compressionResistance * (TYPE)material->getE() /
                         (1.0 - 2.0 * (TYPE)material->getNu());
      // printf("Setting EdivNuFactor[%d]=%G\n", el, EdivNuFactor[el]);
    }
  }
}

template <class TYPE>
TYPE TypeTetElementMaterial<TYPE>::GetCompressionResistanceFactor(
    int elementIndex) {
  return EdivNuFactor[elementIndex];
}

LOBO_TEMPLATE_INSTANT(TypeTetElementMaterial)