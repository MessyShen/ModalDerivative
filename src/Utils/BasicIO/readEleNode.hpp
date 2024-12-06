#include <Eigen/Dense>
#include <fstream>
#include <vector>
#include <iostream>
bool read_node(const std::string &filename, Eigen::MatrixXd &V) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open .node file: " << filename << std::endl;
        return false;
    }

    std::string line;
    std::getline(file, line); 
    std::istringstream iss(line);
    int num_vertices, dimension, num_attributes, num_boundary_markers;
    iss >> num_vertices >> dimension >> num_attributes >> num_boundary_markers;

    V.resize(num_vertices, dimension);
    for (int i = 0; i < num_vertices; ++i) {
        int index;
        file >> index;
        for (int j = 0; j < dimension; ++j) {
            file >> V(i, j);
        }
    }

    file.close();
    return true;
}

bool read_ele(const std::string &filename, Eigen::MatrixXi &F) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open .ele file: " << filename << std::endl;
        return false;
    }

    std::string line;
    std::getline(file, line); 
    std::istringstream iss(line);
    int num_elements, num_nodes_per_element, num_attributes;
    iss >> num_elements >> num_nodes_per_element >> num_attributes;

    F.resize(num_elements, num_nodes_per_element);
    for (int i = 0; i < num_elements; ++i) {
        int index;
        file >> index;
        for (int j = 0; j < num_nodes_per_element; ++j) {
            file >> F(i, j);
        }
    }

    file.close();
    return true;
}

void computeTriangleNorm(const Eigen::Vector3d &v0, const Eigen::Vector3d &v1, const Eigen::Vector3d &v2, Eigen::Vector3d &normal) {
    normal = (v1 - v0).cross(v2 - v0);
    normal.normalize();
}

void correct_face_order(const Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
    for (int i = 0; i < F.rows(); ++i) {
        Eigen::Vector3d nodep[4];
        nodep[0] = V.row(F(i, 0));
        nodep[1] = V.row(F(i, 1));
        nodep[2] = V.row(F(i, 2));
        nodep[3] = V.row(F(i, 3));

        Eigen::Vector3d direction_v;
        computeTriangleNorm(nodep[0], nodep[1], nodep[2], direction_v);
        Eigen::Vector3d n3n0 = nodep[3] - nodep[0];

        if (n3n0.dot(direction_v) > 0) {
            std::swap(F(i, 1), F(i, 2));
        }
    }
}
