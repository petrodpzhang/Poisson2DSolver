#pragma once
#include <vector>
#include <string>
#include <memory>
#include <Eigen/Dense>
#include "Config.hpp"
#include "elements/Element.hpp"

class Mesh {
private:
    std::vector<Node> nodes_;
    std::vector<std::unique_ptr<Element>> elements_;
    double lx_, ly_;
    int nx_, ny_;
    std::string element_type_;
    
public:
    Mesh(double lx, double ly, int nx, int ny, const std::string& element_type);
    
    void generateRectangularMesh();
    void generateTriangleElements();
    void generateQuadElements();
    void applyBoundaryConditions(const std::vector<BoundaryCondition>& bcs);
    
    const std::vector<Node>& getNodes() const { return nodes_; }
    std::vector<Node>& getNodes() { return nodes_; }
    const std::vector<std::unique_ptr<Element>>& getElements() const { return elements_; }
    int getTotalNodes() const { return nodes_.size(); }
    
    // boundary judgment
    bool isOnBoundary(const Node& node, const std::string& boundary) const;
    bool isOnBoundaryAB(const Node& node) const { return node.y == 0; }
    bool isOnBoundaryBC(const Node& node) const { return node.x == lx_; }
    bool isOnBoundaryCD(const Node& node) const { return node.y == ly_; }
    bool isOnBoundaryAD(const Node& node) const { return node.x == 0; }
    
    std::vector<std::pair<double, double>> getNodeCoordinates() const;
    std::vector<std::vector<int>> getElementConnectivity() const;
};