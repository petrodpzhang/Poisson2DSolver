#pragma once
#include "Element.hpp"
#include <vector>

class QuadElement : public Element {
private:
    std::vector<int> node_indices_;
    
public:
    QuadElement(const std::vector<int>& node_indices);
    
    int getNodeCount() const override { return 4; }
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
    std::tuple<Eigen::Vector4d, Eigen::Vector4d, Eigen::Vector4d> 
    shapeFunctions(double xi, double eta) const;
};