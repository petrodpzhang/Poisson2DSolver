#pragma once
#include "Element.hpp"
#include <vector>

class TriangleElement : public Element {
private:
    std::vector<int> node_indices_;
    
public:
    TriangleElement(const std::vector<int>& node_indices);
    
    int getNodeCount() const override { return 3; }
    std::vector<int> getNodeIndices() const override { return node_indices_; }
    
    Eigen::MatrixXd computeStiffnessMatrix(
        const std::vector<Node>& nodes,
        double u_avg,
        const std::function<double(double)>& source_deriv) const override;
    
    Eigen::VectorXd computeLoadVector(
        const std::vector<Node>& nodes,
        const Eigen::VectorXd& u_nodal,
        const std::function<double(double)>& source) const override;
    
    double interpolateSolution(
        const std::vector<Node>& nodes,
        const Eigen::VectorXd& u_nodal,
        double xi, double eta) const override;
    
    std::pair<double, double> transformToGlobal(
        const std::vector<Node>& nodes,
        double xi, double eta) const override;
        
    double jacobianDeterminant(
        const std::vector<Node>& nodes,
        double xi, double eta) const override;
        
private:
    double integrateTriangle(const std::function<double(double, double)>& func) const;
};