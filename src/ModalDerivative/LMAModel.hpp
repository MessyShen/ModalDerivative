#pragma once
#include <Eigen/Sparse>
#include "LoboDynamic/HyperelasticModel.h"
#include "LoboDynamic/LoboVolumtricMesh/LoboTetMesh.h"
#include "Utils/SparseMatrix/SparseMatrixTopology.h"

namespace Lobo
{
    class LMAModel
    {

    public:
        LMAModel();
        ~LMAModel();

        void init(LoboTetMesh *tetMesh,
                  std::string materialType,
                  double youngsModules,
                  double possionRatio,
                  double density);
        void resetMaterial(std::string materialType, double youngsModules, double possionRatio, double density);
        void precompute();
        void loadConstraints(const std::string &filename);
        void stiffnessAtDisplacement(Eigen::VectorXd *displacement, Eigen::SparseMatrix<double> *hessian);
        void computeAccelerationIndices(SparseMatrixTopologyTYPE<double> *sparsetopology);

        void SetFulldisplacementFromConstrained(Eigen::VectorXd &source, Eigen::VectorXd &target, std::vector<int> &rmap); // SetConstraineddisplacementFromFull
        void SetConstraineddisplacementFromFull(Eigen::VectorXd &source,
                                                Eigen::VectorXd &target,
                                                std::vector<int> &rmap);
        void solveLMA();
        void solveUnconstrainedMD();
        void solvePCA();
        void exportModes(const std::string &filePath);
        void exampleLoadModes(const std::string &filePath);
        Eigen::MatrixXd visualizeModes(int i, int j, double scaleCoef);

        int num_DOFs;

        std::vector<Eigen::VectorXd> LMAPsi;
        std::vector<Eigen::VectorXd> MDPhi;
        std::vector<Eigen::VectorXd> ExtensionBasis;
        int modeDim = 5;

        // constraints, derived from constrained model. USED FOR HARD CONSTRAIN
        std::vector<int> constraints;
        std::vector<int> toFullMap;
        std::vector<int> toSubMap;
        int constrained_num;
        int subDOFs; // num dofs - constrain num
        Eigen::MatrixXd LMAphi;
        double scaleCoef = 0.000001;

        Eigen::SparseMatrix<double> mass_matrix;
        Eigen::SparseMatrix<double> stiffness_matrix;
        Eigen::SparseMatrix<double> stiffness_matrix_topology;
        HyperelasticModel *hyperelasticModel;
        LoboTetMesh *tetMesh;
        double mesh_total_mass = 1.0;

    protected:
        // acceleration indices
        int **row_ = NULL;
        int **column_ = NULL;
        int *diagonal_ = NULL;
    };

}