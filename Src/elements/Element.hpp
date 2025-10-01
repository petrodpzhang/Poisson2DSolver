#pragma once
#include <vector>
#include <Eigen/Dense>
#include <memory>

struct Node {
    int id;
    double x, y;
    bool is_boundary = false;
    double bc_value = 0.0;
    bool has_bc = false;
};

class Element {
public:
    virtual ~Element() = default;
    
    virtual int getNodeCount() const = 0;
    virtual std::vector<int> getNodeIndices() const = 0;
    virtual Eigen::MatrixXd computeStiffnessMatrix(
        const std::vector<Node>& nodes,
        double u_avg,
        const std::function<double(double)>& source_deriv) const = 0;
    
    virtual Eigen::VectorXd computeLoadVector(
        const std::vector<Node>& nodes,
        const Eigen::VectorXd& u_nodal,
        const std::function<double(double)>& source) const = 0;
    
    virtual double interpolateSolution(
        const std::vector<Node>& nodes,
        const Eigen::VectorXd& u_nodal,
        double xi, double eta) const = 0;
    
    virtual std::pair<double, double> transformToGlobal(
        const std::vector<Node>& nodes,
        double xi, double eta) const = 0;
        
    virtual double jacobianDeterminant(
        const std::vector<Node>& nodes,
        double xi, double eta) const = 0;
};