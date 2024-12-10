#include "Utils/BasicIO/readEleNode.hpp"
#include "Utils/BasicIO/findFileWithExtension.hpp"
#include "ModalDerivative/LMAModel.hpp"
#include "Utils/SparseMatrix/SparseMatrixRemoveRows.h"
#include "Utils/SparseMatrix/SparseMatrixIO.hpp"

#include <iostream>
#include <filesystem>
#include <vector>
#include <string>

int main(int argc, char *argv[])
{
  if (argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " <directory>" << std::endl;
    return 1;
  }

  fs::path directory = argv[1];

  std::string materialType = "StVK";
  double youngsModulus = 5000;
  double possionRatio = 0.4;
  double density = 1.0;
  int numModes = 5;

  for (int i = 2; i < argc; ++i)
  {
    std::string arg = argv[i];
    if (arg == "--materialType" && i + 1 < argc)
    {
      std::string value = argv[++i];
      if (value == "NeoHookean" || value == "StVK")
      {
        materialType = value;
      }
      else
      {
        std::cerr << "Error: materialType must be either 'NeoHookean' or 'StVK'." << std::endl;
        return 1;
      }
    }
    else if (arg == "--youngsModulus" && i + 1 < argc)
    {
      youngsModulus = std::stod(argv[++i]);
    }
    else if (arg == "--possionRatio" && i + 1 < argc)
    {
      possionRatio = std::stod(argv[++i]);
    }
    else if (arg == "--density" && i + 1 < argc)
    {
      density = std::stod(argv[++i]);
    }
    else if (arg == "--numModes" && i + 1 < argc)
    {
      numModes = std::stod(argv[++i]);
    }
  }

  std::vector<std::string> extensions = {".ele", ".node", ".constrainedDoFs"};
  std::string eleFile, nodeFile, consFile, stiffFile, massFile;

  eleFile = findFileWithExtension(directory, extensions[0]);
  nodeFile = findFileWithExtension(directory, extensions[1]);
  consFile = findFileWithExtension(directory, extensions[2]);

  Eigen::MatrixXd V;
  Eigen::MatrixXi T;
  if (!read_node(nodeFile, V) || !read_ele(eleFile, T)) {
        return 1;
  }

  correct_face_order(V, T);

  
  Lobo::LMAModel lma_model;
  if (consFile.empty()) {
    std::cout << "Warning: constrainedDoFs file not found." << std::endl;
  } else {
    lma_model.loadConstraints(consFile);
  }

  lma_model.modeDim = numModes;
  
  Lobo::LoboTetMesh *tetMesh = new Lobo::LoboTetMesh();
  tetMesh->setTetraheralMesh(V, T);
  tetMesh->computeDiagMassMatrix(&lma_model.mass_matrix);
  lma_model.init(tetMesh, materialType, youngsModulus, possionRatio, density);
  int lma_i = 1, lma_j = 1;

  std::string fileDialogResult = directory.string() + "/export";
  lma_model.solveUnconstrainedMD();
  lma_model.exportModes(fileDialogResult);
  std::cout << "Please see LMAModel::exampleLoadModes(const std::string &filePath)" << std::endl;

}