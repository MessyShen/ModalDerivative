#include "DynamicModel.h"

Lobo::DynamicModel::DynamicModel(){
	num_DOFs = 0;
	is_sparse_sovler = true;
	trigger = true;
	row_ = NULL;
	column_ = NULL;
}

Lobo::DynamicModel::~DynamicModel() {}

void Lobo::DynamicModel::setAccelerationIndices(int **row, int **column) {
  this->row_ = row;
  this->column_ = column;
}

void Lobo::DynamicModel::setAccelerationDiagIndices(int *diagonal) {
  this->diagonal_ = diagonal;
}