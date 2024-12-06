#include "LMAModel.hpp"
#include "Utils/SparseMatrix/SparseMatrixRemoveRows.h"
#include <fstream>
#include <iostream>
#include <Spectra/SymGEigsSolver.h>
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <Spectra/MatOp/SparseCholesky.h>
#include "Functions/eigenVecMat.hpp"
#include "Functions/EigenMatrixIO.h"

Lobo::LMAModel::LMAModel()
{
}

Lobo::LMAModel::~LMAModel()
{
}

void Lobo::LMAModel::init(LoboTetMesh *_tetMesh,
                          std::string materialType,
                          double youngsModules,
                          double possionRatio,
                          double density)
{
  tetMesh = _tetMesh;
  num_DOFs = tetMesh->getNumVertex() * 3;
  constrained_num = constraints.size();
  // create hyperelasticmodel
  hyperelasticModel = new HyperelasticModel(tetMesh);
  hyperelasticModel->setMaterial(materialType, youngsModules, possionRatio, density);
  // hyperelasticModel->precompute();

  tetMesh->computeDiagMassMatrix(&mass_matrix);

  hyperelasticModel->computeStiffnessMatrixTopology(
      &stiffness_matrix_topology);
  stiffness_matrix = stiffness_matrix_topology;
  // in this model, the final system matrix has the same topology as stiffness
  // matrix
  SparseMatrixTopologyTYPE<double> sparseMatrixTopology(
      &stiffness_matrix_topology);

  computeAccelerationIndices(&sparseMatrixTopology);

  hyperelasticModel->precompute();
  hyperelasticModel->setAccelerationIndices(row_, column_);
}

void Lobo::LMAModel::resetMaterial(std::string materialType,
                                   double youngsModules,
                                   double possionRatio,
                                   double density)
{
  hyperelasticModel->setMaterial(materialType, youngsModules, possionRatio, density);
  hyperelasticModel->precompute();
}


void Lobo::LMAModel::loadConstraints(const std::string &filename)
{
  std::ifstream file(filename);
  if (!file.is_open())
  {
    std::cerr << "Failed to open file: " << filename << std::endl;
    return;
  }

  file >> constrained_num;

  constraints.resize(constrained_num);
  for (int i = 0; i < constrained_num; ++i)
  {
    file >> constraints[i];
  }

  file.close();
}

void Lobo::LMAModel::SetFulldisplacementFromConstrained(
    Eigen::VectorXd &source, Eigen::VectorXd &target, std::vector<int> &rmap)
{
  target.resize(num_DOFs);
  target.setZero();
  int rmap_size = rmap.size();
  for (int i = 0; i < rmap_size; ++i)
  {
    target[rmap[i]] = source[i];
  }
}

void Lobo::LMAModel::SetConstraineddisplacementFromFull(
    Eigen::VectorXd &source, Eigen::VectorXd &target, std::vector<int> &rmap)
{
  target.resize(subDOFs);
  target.setZero();
  int rmap_size = rmap.size();
  for (int i = 0; i < rmap_size; ++i)
  {
    target[i] = source[rmap[i]];
  }
}

void Lobo::LMAModel::solveLMA()
{
  Eigen::SparseMatrix<double> *hessian = &stiffness_matrix_topology;
  // stiffnessAtDisplacement(hessian);


  int numMode = modeDim;

  Eigen::SparseMatrix<double> subStiffness, subMassMatrix;
  // std::vector<int> constraint;
  // int constraint_num = 0;
  // constrainmodel->getConstraints(constraint, constraint_num);

  std::vector<int> consMap, invMap;

  int sub_size = num_DOFs - constrained_num;
  subDOFs = sub_size;
  subStiffness.resize(sub_size, sub_size);
  subMassMatrix.resize(sub_size, sub_size);

  // q_vel.resize(sub_size);
  // q_vel.setZero();
  // q_1.resize(sub_size);
  // q_1.setZero();

  createReverseMapByConstrains(consMap, invMap, num_DOFs, constrained_num,
                               constraints);

  Eigen::VectorXd zeroDis;
  zeroDis.resize(num_DOFs);
  zeroDis.setZero();
  stiffnessAtDisplacement(&zeroDis, hessian);

  std::cout << "subSparseMatrix..." << std::endl;

  subSparseMatrix(*hessian, subStiffness, consMap);
  subSparseMatrix(mass_matrix, subMassMatrix, consMap);

  // Eigen::MatrixXd Adense = subStiffness;
  double shift = 0.000001;
  Eigen::SparseMatrix<double> K = subStiffness - shift * subMassMatrix;
  Spectra::SparseSymMatProd<double> Aop(subMassMatrix);
  Spectra::SparseCholesky<double> Bop(K);

  std::cout << "Init GEigs..." << std::endl;
  Spectra::SymGEigsSolver<double, Spectra::LARGEST_MAGN,
                          // Spectra::DenseSymMatProd<double>,
                          Spectra::SparseSymMatProd<double>,
                          Spectra::SparseCholesky<double>,
                          Spectra::GEIGS_CHOLESKY>
      geigs(&Aop, &Bop, numMode, numMode * 5);

  geigs.init();

  std::cout << "solving GEigs..." << std::endl;
  int nconv = geigs.compute();

  std::cout << "converge eige" << nconv << std::endl;

  Eigen::VectorXd evalues;
  Eigen::MatrixXd evecs;
  if (geigs.info() == Spectra::SUCCESSFUL)
  {
    evalues = geigs.eigenvalues();
    evecs = geigs.eigenvectors();
  }
  for (int i = 0; i < numMode; ++i)
  {
    evalues(i) = 1.0 / (evalues(i) - shift);
    double evec_mass_norm = evecs.col(i).dot(subMassMatrix * evecs.col(i));
    evecs.col(i) = evecs.col(i) / sqrt(evec_mass_norm);
  }
  LMAphi = evecs;
  std::cout << "evalues\n"
            << evalues << std::endl;

  for (int numID = 0; numID < numMode; ++numID)
  {
    std::cout << "eigenvector " << numID << std::endl;
    Eigen::VectorXd phi1 = evecs.col(numID);
    Eigen::VectorXd fullphi1;
    fullphi1.setZero();
    phi1 = evalues[0] / evalues[numID] * phi1;
    SetFulldisplacementFromConstrained(phi1, fullphi1, invMap);
    LMAPsi.push_back(fullphi1);
    /*for (int i = 0; i < fullphi1.size(); ++i) {
      std::cout << fullphi1[i] << " ";
    }
    std::cout << "\n\n\n\n\n";*/
  }

  for (int psiID = 0; psiID < numMode; ++psiID)
  {

    Eigen::VectorXd PSI = evecs.col(psiID);

    for (int p = 0; p < 9; ++p)
    {
      Eigen::VectorXd basisExt, fullBasisExt;
      basisExt.resizeLike(PSI);

      Eigen::MatrixXd rotMat(3, 3);
      rotMat.setZero();
      rotMat.data()[p] = 1.0;
      for (int i = 0; i < subDOFs / 3; ++i)
      {
        Eigen::VectorXd subvec(3);
        subvec.setZero();
        subvec[0] = PSI[i * 3 + 0];
        subvec[1] = PSI[i * 3 + 1];
        subvec[2] = PSI[i * 3 + 2];
        Eigen::VectorXd rotedVec = rotMat * subvec;
        basisExt[i * 3 + 0] = rotedVec[0];
        basisExt[i * 3 + 1] = rotedVec[1];
        basisExt[i * 3 + 2] = rotedVec[2];
      }

      SetFulldisplacementFromConstrained(basisExt, fullBasisExt, invMap);
      ExtensionBasis.push_back(fullBasisExt);
    }
  }

  std::cout << "geometry collected\n"
            << std::endl;

  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> LDLT;
  Eigen::LDLT<Eigen::MatrixXd> reduceLDLT;
  // LDLT.pardisoParameterArray()[33] = 1;
  // LDLT.pardisoParameterArray()[59] = 0;
  //  auto solver_start = std::chrono::high_resolution_clock::now();
  Eigen::MatrixXd reduce_stiffness = evecs.transpose() * subStiffness * evecs;
  LDLT.compute(subStiffness);
  reduceLDLT.compute(reduce_stiffness);

  // Eigen::VectorXd subGravity;
  // SetConstraineddisplacementFromFull(gravity_force, subGravity, invMap);
  // Eigen::VectorXd sub_ku_g = LDLT.solve(subGravity);
  // Eigen::VectorXd sub_ku_g = evecs.transpose() * subGravity;
  // Eigen::VectorXd reduce_ku_g = reduceLDLT.solve(sub_ku_g);
  // Eigen::VectorXd unProjectKuG = evecs * reduce_ku_g;
  // Eigen::VectorXd full_ku_g;
  // SetFulldisplacementFromConstrained(unProjectKuG, full_ku_g, invMap);
  // LMAPsi.push_back(full_ku_g);

  // check stiffness and mass projection
  // Eigen::MatrixXd project_stiff_mat = LMAphi.transpose() * subStiffness * LMAphi;
  // std::cout << "DEBUG\n" << project_stiff_mat << std::endl;

  // Eigen::MatrixXd project_mass_mat =
  //     LMAphi.transpose() * subMassMatrix * LMAphi;
  // std::cout << "DEBUG\n" << project_mass_mat << std::endl;

  int count_md = 0;
  for (int i = 0; i < numMode; ++i)
  {
    for (int j = i; j < numMode; ++j)
    {
      count_md += 1;

      Eigen::SparseMatrix<double> K_i;
      K_i.resize(sub_size, sub_size);

      Eigen::VectorXd PSI_i = evecs.col(i);
      Eigen::VectorXd PSI_j = evecs.col(j);
      PSI_i = PSI_i / (PSI_i.dot(PSI_i));
      PSI_j = PSI_j / (PSI_j.dot(PSI_j));

      Eigen::VectorXd fullPSI_i, scaled_fullPSI_i;
      SetFulldisplacementFromConstrained(PSI_i, fullPSI_i, invMap);
      scaled_fullPSI_i = scaleCoef * fullPSI_i;
      // fullPSI_i.setZero();
      stiffnessAtDisplacement(&scaled_fullPSI_i, hessian);
      subSparseMatrix(*hessian, K_i, consMap);
      Eigen::VectorXd rhs = K_i * PSI_j;
      Eigen::VectorXd PSI_ij = LDLT.solve(rhs);
      PSI_ij = (PSI_j - PSI_ij) / scaleCoef;

      double norm_psi_ij = PSI_ij.dot(subMassMatrix * PSI_ij);
      PSI_ij = PSI_ij / sqrt(norm_psi_ij);
      norm_psi_ij = PSI_ij.dot(subMassMatrix * PSI_ij);

      PSI_ij = evalues[0] * evalues[0] /
               evalues[i] / evalues[j] * PSI_ij;

      // std::cout << "norm_psi_ij  " << norm_psi_ij << std::endl;

      Eigen::VectorXd fullPSI_ij;
      fullPSI_ij.setZero();
      SetFulldisplacementFromConstrained(PSI_ij, fullPSI_ij, invMap);

      // Eigen::MatrixXd dense_Ki = K_i;
      // std::cout << dense_Ki << std::endl;

      std::cout << i << " " << j << "ID: count_md" << count_md << std::endl;

      MDPhi.push_back(fullPSI_ij);
      /*for (int p = 0; p < fullPSI_ij.size(); ++p) {
        std::cout << fullPSI_ij[p] << " ";
      }
      std::cout << "\n\n\n\n\n";*/
    }
  }

  // solvePCA();
  return;
}

void Lobo::LMAModel::solvePCA()
{
  return;
}

Eigen::MatrixXd Lobo::LMAModel::visualizeModes(int lma_i, int lma_j, double scaleCoef)
{
  Eigen::VectorXd shift;
  shift = LMAPsi[0];
  shift.setZero();


  if (lma_i < 1){
    shift.setZero();
  } 
  else if (lma_j < 1)
    shift = LMAPsi[lma_i - 1];
  else if (lma_j >= lma_i) {
    int numNode = modeDim;
    int ijID = lma_j - lma_i + 1;
    for (int t = 1; t < lma_i; ++t) {
      ijID += numNode + 1 - t;
    }
    shift = MDPhi[ijID - 1];
    std::cout << lma_i << " " << lma_j << " : " << ijID << std::endl; 
  } else if (lma_i >= modeDim) {
    std::cout << "lma_i out of range" << std::endl;
    shift.setZero();
  } else if (lma_j < 1 || lma_j >= modeDim) {
    std::cout << "lma_j out of range" << std::endl;
    shift.setZero();
  } else {
    std::cout << "make sure lma_i <= lma_j" << std::endl;
    shift.setZero();
  }

  //if (lma_i < 1 || lma_i > kineticmodel->modeDim || lma_j < 1 || lma_j >= 10) {
  //  shift.setZero();
  //} else {
  //  int numNode = kineticmodel->modeDim;
  //  shift = kineticmodel->ExtensionBasis[numNode * (lma_i - 1) + lma_j];
  //}
  Eigen::MatrixXd shiftMat = scaleCoef * Lobo::eigen_vec_2_mat(shift, tetMesh->getNumVertex(), 3);
  
  return shiftMat;
}

void Lobo::LMAModel::exportModes(const std::string &filePath)
{
  // make directory for filePath using filesystem
  if (!std::filesystem::exists(filePath))
    std::filesystem::create_directories(filePath);
  // export LMA modes and MD modes
  std::ofstream lma_file(filePath + "/LMAphi.bin", std::ios::binary);
  std::ofstream md_file(filePath + "/MDPhi.bin", std::ios::binary);
  EigenMatrixIO::write_vector_eigen(lma_file, LMAPsi);
  EigenMatrixIO::write_vector_eigen(md_file, MDPhi);
  lma_file.close();
  md_file.close();
  std::cout << "Exported LMA and MD modes" << std::endl;
  std::cout << "LMA modes: " << LMAPsi.size() << std::endl;
  std::cout << "MD modes: " << MDPhi.size() << std::endl;
  std::cout << "at " << filePath << std::endl;
}

void Lobo::LMAModel::exampleLoadModes(const std::string &filePath)
{
  // load LMA and MD modes
  std::string lma_file(filePath + "/LMAphi.bin");
  std::string md_file(filePath + "/MDPhi.bin");
  std::vector<Eigen::VectorXd> example_LMAPsi, example_MDPhi;
  EigenMatrixIO::read_vector_eigen(lma_file.c_str(), example_LMAPsi);
  EigenMatrixIO::read_vector_eigen(md_file.c_str(), example_MDPhi); 
  std::cout << "Loaded LMA and MD modes" << std::endl;
  std::cout << "LMA modes: " << example_LMAPsi.size() << std::endl;
  std::cout << "MD modes: " << example_MDPhi.size() << std::endl;
  if (example_LMAPsi.size() > 0)
    std::cout << "each mode size: " << example_LMAPsi[0].size() << std::endl;
}

void Lobo::LMAModel::stiffnessAtDisplacement(Eigen::VectorXd *free_variables, Eigen::SparseMatrix<double> *hessian)
{
  double energy = 0.0;
  Eigen::VectorXd jacobi;
  jacobi.resize(free_variables->size());

  int flags_all = 0;
  flags_all |= Computeflags_energy | Computeflags_fisrt | Computeflags_second |
               Computeflags_reset;

  hyperelasticModel->computeEnergySparse(free_variables, &energy, &jacobi,
                                         hessian, flags_all);
}

void Lobo::LMAModel::computeAccelerationIndices(
    SparseMatrixTopologyTYPE<double> *sparsetopology)
{
  int numElements = tetMesh->getNumElements();

  if (row_ != NULL)
  {
    for (int i = 0; i < numElements; i++)
      free(row_[i]);
    free(row_);
  }

  if (column_ != NULL)
  {
    for (int i = 0; i < numElements; i++)
      free(column_[i]);
    free(column_);
  }

  row_ = (int **)malloc(sizeof(int *) * numElements);
  column_ = (int **)malloc(sizeof(int *) * numElements);
  diagonal_ = (int *)malloc(sizeof(int) * tetMesh->getNumVertex() * 3); // for diagnal index

  for (int i = 0; i < tetMesh->getNumVertex() * 3; i++)
  {
    diagonal_[i] = sparsetopology->getValueIndex(i, i);
  }

  int numElementVertices = tetMesh->numElementVertices;
  for (int el = 0; el < numElements; el++)
  {
    // std::cout << el << "\r";
    row_[el] = (int *)malloc(sizeof(int) * numElementVertices);
    column_[el] = (int *)malloc(sizeof(int) * numElementVertices *
                                numElementVertices * 3 * 3);

    for (int ver = 0; ver < numElementVertices; ver++)
    {
      row_[el][ver] = tetMesh->tet_indices[el * numElementVertices + ver];
    }

    for (int i = 0; i < numElementVertices; i++)
      for (int j = 0; j < numElementVertices; j++)
      {
        for (int k = 0; k < 3; k++)
        {
          for (int l = 0; l < 3; l++)
          {
            int block_r = i * 3 + k;
            int block_c = j * 3 + l;

            column_[el][3 * numElementVertices * block_c +
                        block_r] =
                sparsetopology->getValueIndex(3 * row_[el][i] + k,
                                              3 * row_[el][j] + l);
          }
        }
      }
  }
}