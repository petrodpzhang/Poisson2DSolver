#include "Mesh.hpp"
#include "elements/TriangleElement.hpp"
#include "elements/QuadElement.hpp"
#include <algorithm>
#include <stdexcept>

Mesh::Mesh(double lx, double ly, int nx, int ny, const std::string& element_type)
    : lx_(lx), ly_(ly), nx_(nx), ny_(ny), element_type_(element_type) {
    
    generateRectangularMesh();
}

void Mesh::generateRectangularMesh() {
    nodes_.clear();
    elements_.clear();
    
    // generate nodes
    double dx = lx_ / nx_;
    double dy = ly_ / ny_;
    
    int node_id = 0;
    for (int j = 0; j <= ny_; ++j) {
        for (int i = 0; i <= nx_; ++i) {
            Node node;
            node.id = node_id++;
            node.x = i * dx;
            node.y = j * dy;
            node.is_boundary = (i == 0 || i == nx_ || j == 0 || j == ny_);
            nodes_.push_back(node);
        }
    }
    
    // generate elements
    if (element_type_ == "triangle") {
        generateTriangleElements();
    } else if (element_type_ == "quad") {
        generateQuadElements();
    } else {
        throw std::invalid_argument("Unknown element type: " + element_type_);
    }
}

void Mesh::generateTriangleElements() {
    for (int j = 0; j < ny_; ++j) {
        for (int i = 0; i < nx_; ++i) {
            // each rectangle is divided into two triangles
            int n0 = j * (nx_ + 1) + i;
            int n1 = n0 + 1;
            int n2 = n0 + (nx_ + 1);
            int n3 = n2 + 1;
            
            // first triangle (n0, n1, n3)
            elements_.push_back(std::make_unique<TriangleElement>(
                std::vector<int>{n0, n1, n3}));
            
            // second triangle (n0, n3, n2)  
            elements_.push_back(std::make_unique<TriangleElement>(
                std::vector<int>{n0, n3, n2}));
        }
    }
}

void Mesh::generateQuadElements() {
    for (int j = 0; j < ny_; ++j) {
        for (int i = 0; i < nx_; ++i) {
            int n0 = j * (nx_ + 1) + i;
            int n1 = n0 + 1;
            int n2 = n0 + (nx_ + 1) + 1;
            int n3 = n0 + (nx_ + 1);
            
            elements_.push_back(std::make_unique<QuadElement>(
                std::vector<int>{n0, n1, n2, n3}));
        }
    }
}

void Mesh::applyBoundaryConditions(const std::vector<BoundaryCondition>& bcs) {
    for (const auto& bc : bcs) {
        for (auto& node : nodes_) {
            if (isOnBoundary(node, bc.boundary)) {
                node.has_bc = true;
                node.bc_value = bc.value;
            }
        }
    }
}

bool Mesh::isOnBoundary(const Node& node, const std::string& boundary) const {
    if (boundary == "AB") return node.y == 0;
    if (boundary == "BC") return node.x == lx_;
    if (boundary == "CD") return node.y == ly_;
    if (boundary == "AD") return node.x == 0;
    
    // handle combined boundaries
    if (boundary == "ABCD") return node.is_boundary;
    
    throw std::invalid_argument("Unknown boundary: " + boundary);
}

std::vector<std::pair<double, double>> Mesh::getNodeCoordinates() const {
    std::vector<std::pair<double, double>> coords;
    for (const auto& node : nodes_) {
        coords.emplace_back(node.x, node.y);
    }
    return coords;
}

std::vector<std::vector<int>> Mesh::getElementConnectivity() const {
    std::vector<std::vector<int>> connectivity;
    for (const auto& element : elements_) {
        connectivity.push_back(element->getNodeIndices());
    }
    return connectivity;
}