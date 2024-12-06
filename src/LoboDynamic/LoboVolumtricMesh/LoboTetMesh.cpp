#include "LoboTetMesh.h"
// #include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/project.h>
#include <igl/barycentric_coordinates.h>
#include <igl/in_element.h>
#include <igl/AABB.h>
#include <igl/per_vertex_normals.h>
#include <igl/tet_tet_adjacency.h>

// boost headers

#include "Functions/EigenMatrixIO.h"
#include "Functions/computeTriangle.h"
#include "Functions/findElementInVector.h"
#include "Functions/eigenVecMat.hpp"

// #include "LoboDynamic/LoboDynamicScene.h"
// #include "LoboMesh/LoboMesh.h"
// #include "imgui.h"
// #include "LoboImGui/cpp/imgui_stdlib.h"
// #include "Utils/glm/glmEigenConverter.h"
// #include "Utils/glmMyFunctions.h"


Lobo::LoboTetMesh::LoboTetMesh() {
    initializedGL = false;
    render_normals = false;
    render_cubature_points = false;
    //tetgen_command = "pq1.414";
    tetgen_command ="pq3.0";
    status_flags = 0;

    mesh_total_volume = 0.0;
    numElementVertices = 4;

    clicked_face = 0;
    tet_ids_to_show.clear();
}

Lobo::LoboTetMesh::~LoboTetMesh() {}


void Lobo::LoboTetMesh::setTetraheralMesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &T) {
    tet_vertice = Lobo::eigen_mat_2_vec(V);
    tet_indices = Lobo::eigen_mat_2_vec(T);
    ori_tet_vertice = tet_vertice;
    reinitialTetMesh();
}

void Lobo::LoboTetMesh::reinitialTetMesh() {
    vertices_flags.resize(tet_vertice.size() / 3);
    std::fill(vertices_flags.begin(), vertices_flags.end(), 0);

    tet_vertice_col =
        Lobo::eigen_vec_2_mat(tet_vertice, tet_vertice.size() / 3, 3);
    // tet_faces_col = Lobo::eigen_vec_2_mat(tet_faces, tet_faces.size() / 3, 3);
    tet_indices_col =
        Lobo::eigen_vec_2_mat(tet_indices, tet_indices.size() / 4, 4);

    tet_vertice_force.resize(tet_vertice.size());
    tet_vertice_force.setZero();

    // default material
    setAllMaterial(1.0, 1000.0, 0.4);

    ori_tet_vertice = tet_vertice;

    status_flags &= ~TetMeshStatusFlags_precomputed;

    precomputeElementData();
}


void Lobo::LoboTetMesh::setAllMaterial(double density, double youngsmodulus,
                                       double possionratio) {
    materials.clear();
    materialid.resize(getNumElements());
    std::fill(materialid.begin(), materialid.end(), 0);

    Material m;
    m.density = density;
    m.youngsmodulus = youngsmodulus;
    m.possionratio = possionratio;
    materials.push_back(m);
}

void Lobo::LoboTetMesh::computeDiagMassMatrix(
    Eigen::SparseMatrix<double> *mass) {
    int r = getNumVertex() * 3;
    mass->resize(r, r);
    typedef Eigen::Triplet<double> EIGEN_TRI_T;
    std::vector<EIGEN_TRI_T> coefficients;
    int numElements = getNumElements();
    for (int i = 0; i < numElements; i++) {
        double vol = elements_data[i].volume;
        double density = materials[materialid[i]].density;
        for (int j = 0; j < 4; j++) {
            double mass_value = 0;
            mass_value = density * vol / 4.0;
            int row = tet_indices[i * 4 + j] * 3;
            for (int d = 0; d < 3; d++) {
                coefficients.push_back(
                    EIGEN_TRI_T(row + d, row + d, mass_value));
            }
        }
    }
    mass->setFromTriplets(coefficients.begin(), coefficients.end());
}

void Lobo::LoboTetMesh::computeDiagMassMatrix(
    Eigen::SparseMatrix<double> *mass, Eigen::VectorXd *mass_diag) {
  int r = getNumVertex() * 3;
  mass->resize(r, r);
  mass_diag->resize(r);
  mass_diag->setZero();
  typedef Eigen::Triplet<double> EIGEN_TRI_T;
  std::vector<EIGEN_TRI_T> coefficients;
  int numElements = getNumElements();
  for (int i = 0; i < numElements; i++) {
    double vol = elements_data[i].volume;
    double density = materials[materialid[i]].density;
    for (int j = 0; j < 4; j++) {
      double mass_value = 0;
      mass_value = density * vol / 4.0;
      int row = tet_indices[i * 4 + j] * 3;
      for (int d = 0; d < 3; d++) {
        coefficients.push_back(EIGEN_TRI_T(row + d, row + d, mass_value));
        mass_diag->data()[row+d] += mass_value;
      }
    }
  }
  mass->setFromTriplets(coefficients.begin(), coefficients.end());
}


void Lobo::LoboTetMesh::precomputeElementData() {
    if (status_flags & TetMeshStatusFlags_precomputed) {
        // already precomputed
        return;
    }
    int numElements = getNumElements();
    elements_data.resize(numElements);

    mesh_total_volume = 0.0;
    for (int i = 0; i < numElements; i++) {
        correctElementNodeOrder(i);
        computeElementVolume(i);
        computeElementShapeFunctionDerivate(i);  // precompute dF_du
        mesh_total_volume += elements_data[i].volume;

        // update element center
        elements_data[i].center_p.setZero();
        for (int j = 0; j < 4; j++) {
            int nodeid = tet_indices.data()[i * 4 + j];
            elements_data[i].center_p.data()[0] +=
                tet_vertice.data()[nodeid * 3 + 0];
            elements_data[i].center_p.data()[1] +=
                tet_vertice.data()[nodeid * 3 + 1];
            elements_data[i].center_p.data()[2] +=
                tet_vertice.data()[nodeid * 3 + 2];
        }
        elements_data[i].center_p /= 4;
    }

    // generateBarycentricCoordinate();

    status_flags |= TetMeshStatusFlags_precomputed;
}

void Lobo::LoboTetMesh::getNodeElements(
    std::vector<std::vector<int>> &node_elements) {
    node_elements.resize(getNumVertex());
    int num_elements = getNumElements();
    for (int i = 0; i < num_elements; i++) {
        for (int j = 0; j < numElementVertices; j++) {
            int nodeid = tet_indices.data()[i * 4 + j];
            node_elements[nodeid].push_back(i);
        }
    }
}

void Lobo::LoboTetMesh::getNodeRestPosition(int nodeid, Eigen::Vector3d &p) {
    for (int j = 0; j < 3; j++) {
        p.data()[j] = ori_tet_vertice.data()[nodeid * 3 + j];
    }
}

Eigen::Vector3d Lobo::LoboTetMesh::getNodeRestPosition(int nodeid) {
    Eigen::Vector3d tmp;
    getNodeRestPosition(nodeid, tmp);
    return tmp;
}

Eigen::Vector3d Lobo::LoboTetMesh::getNodeCurPosition(int nodeid)
{
    Eigen::Vector3d tmp;
    tmp.setZero();
    for (int j = 0; j < 3; j++) {
        tmp.data()[j] = tet_vertice.data()[nodeid * 3 + j];
    }
    return tmp;
}

Eigen::Vector3d Lobo::LoboTetMesh::getNodeNormal(int nodeid)
{
    Eigen::Vector3d tmp;
    tmp.setZero();
    if(surface_vertce_map[nodeid] == -1)
    {
        return tmp;
    }
    int sur_node_id = surface_vertce_map[nodeid];

    tmp.data()[0] = tet_sur_vertice_normal.data()[tet_sur_vertice_normal.rows()*0+sur_node_id];
    tmp.data()[1] = tet_sur_vertice_normal.data()[tet_sur_vertice_normal.rows()*1+sur_node_id];
    tmp.data()[2] = tet_sur_vertice_normal.data()[tet_sur_vertice_normal.rows()*2+sur_node_id];

    return tmp;
}

void Lobo::LoboTetMesh::correctElementNodeOrder(int elementid) {
    int ni[4];
    Eigen::Vector3d node_p[4];
    for (int i = 0; i < 4; i++) {
        ni[i] = tet_indices[elementid * 4 + i];
        for (int j = 0; j < 3; j++) {
            getNodeRestPosition(ni[i], node_p[i]);
        }
    }
    Eigen::Vector3d direction_v;
    Lobo::computeTriangleNorm(node_p[0], node_p[1], node_p[2], direction_v);
    Eigen::Vector3d n3n0 = node_p[3] - node_p[0];
    if (n3n0.dot(direction_v) < 0) {
        // element->node_indices[1] = n2;
        // element->node_indices[2] = n1;
        tet_indices[elementid * 4 + 1] = ni[2];
        tet_indices[elementid * 4 + 2] = ni[1];
        // std::cout << "bad order element" << std::endl;
    }
}

void Lobo::LoboTetMesh::computeElementVolume(int elementid) {
    Eigen::Vector3d a =
        this->getNodeRestPosition(tet_indices[elementid * 4 + 0]);
    Eigen::Vector3d b =
        this->getNodeRestPosition(tet_indices[elementid * 4 + 1]);
    Eigen::Vector3d c =
        this->getNodeRestPosition(tet_indices[elementid * 4 + 2]);
    Eigen::Vector3d d =
        this->getNodeRestPosition(tet_indices[elementid * 4 + 3]);

    elements_data[elementid].volume = Lobo::computeTetVolumeABS(a, b, c, d);
}

void Lobo::LoboTetMesh::computeElementShapeFunctionDerivate(int elementid) {
    int ni[4];
    Eigen::Vector3d node_p[4];
    for (int i = 0; i < 4; i++) {
        ni[i] = tet_indices[elementid * 4 + i];
        for (int j = 0; j < 3; j++) {
            getNodeRestPosition(ni[i], node_p[i]);
        }
    }

    TetElementData *te = &elements_data[elementid];

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            te->Dm.data()[j * 3 + i] =
                node_p[j].data()[i] - node_p[3].data()[i];
        }
    }
    te->Dm_inverse = te->Dm.inverse();

    Eigen::Matrix4d referenceShape;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 3; j++) {
            referenceShape.data()[i * 4 + j] = node_p[i].data()[j];
        }
        referenceShape.data()[i * 4 + 3] = 1;
    }

    // need test
    Eigen::Matrix4d inverseShapefunction = referenceShape.inverse();
    te->shape_function_inv = inverseShapefunction;

    te->Phi_derivate.resize(4, 3);
    for (int i = 0; i < 4; i++) {
        te->Phi_derivate.data()[0 * 4 + i] =
            inverseShapefunction.data()[0 * 4 + i];
        te->Phi_derivate.data()[1 * 4 + i] =
            inverseShapefunction.data()[1 * 4 + i];
        te->Phi_derivate.data()[2 * 4 + i] =
            inverseShapefunction.data()[2 * 4 + i];
    }
}

void Lobo::LoboTetMesh::correct_face_order(int eleid, std::vector<int> face_order,Eigen::Vector3i &face_index)
{
    int n[4];
    Eigen::Vector3d nodep[4];

    for(int i=0;i<4;i++)
    {

        n[i] = tet_indices_col.data()[face_order[i] * tet_indices_col.rows() + eleid];
        nodep[i] = this->getNodeRestPosition(n[i]);
    }

    //compute face_normal
    Eigen::Vector3d direction_v;
    Lobo::computeTriangleNorm(nodep[0], nodep[1], nodep[2], direction_v);
    Eigen::Vector3d n3n0 = nodep[3] - nodep[0];
    if (n3n0.dot(direction_v) > 0) {
        // element->node_indices[1] = n2;
        // element->node_indices[2] = n1;
        face_index[1] = n[2];
        face_index[2] = n[1];
    }
}
