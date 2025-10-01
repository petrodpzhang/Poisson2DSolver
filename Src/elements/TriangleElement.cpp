#include "TriangleElement.hpp"
#include <boost/math/quadrature/gauss.hpp>
#include <iostream>

TriangleElement::TriangleElement(const std::vector<int>& node_indices)
    : node_indices_(node_indices) {
    if (node_indices.size() != 3) {
        throw std::invalid_argument("Triangle element requires exactly 3 nodes");
    }
}

Eigen::MatrixXd TriangleElement::computeStiffnessMatrix(
    const std::vector<Node>& nodes,
    double u_avg,
    const std::function<double(double)>& source_deriv) const {
    
    Eigen::MatrixXd ke = Eigen::MatrixXd::Zero(3, 3);
    
    // get node coordinates
    const Node& n0 = nodes[node_indices_[0]];
    const Node& n1 = nodes[node_indices_[1]];
    const Node& n2 = nodes[node_indices_[2]];
    
    // compute Jacobian matrix
    Eigen::Matrix2d J;
    J << n1.x - n0.x, n2.x - n0.x,
         n1.y - n0.y, n2.y - n0.y;
    
    double detJ = std::abs(J.determinant());
    
    if (detJ < 1e-12) {
        throw std::runtime_error("Element Jacobian determinant too small");
    }
    
    Eigen::Matrix2d invJ = J.inverse();
    
    // reference coordinate system shape function gradient
    Eigen::Matrix<double, 2, 3> dN_ref;
    dN_ref << -1, 1, 0,
              -1, 0, 1;
    
    // physical coordinate system shape function gradient
    Eigen::Matrix<double, 2, 3> dN_global = invJ * dN_ref;
    
    // linear stiffness matrix part 
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            ke(i, j) = dN_global.col(i).dot(dN_global.col(j)) * detJ * 0.5;
        }
    }
    
    // nonlinear stiffness matrix part (N_i * N_j * f'(u))
    if (source_deriv) {
        double f_prime = source_deriv(u_avg);
        
        // use analytical integration to compute mass matrix
        Eigen::Matrix3d mass_matrix;
        mass_matrix << 2, 1, 1,
                       1, 2, 1, 
                       1, 1, 2;
        mass_matrix *= detJ / 24.0;
        
        ke -= f_prime * mass_matrix;
    }
    
    return ke;
}

Eigen::VectorXd TriangleElement::computeLoadVector(
    const std::vector<Node>& nodes,
    const Eigen::VectorXd& u_nodal,
    const std::function<double(double)>& source) const {
    
    Eigen::VectorXd fe = Eigen::VectorXd::Zero(3);
    
    const Node& n0 = nodes[node_indices_[0]];
    const Node& n1 = nodes[node_indices_[1]];
    const Node& n2 = nodes[node_indices_[2]];
    
    // compute Jacobian matrix and determinant
    Eigen::Matrix2d J;
    J << n1.x - n0.x, n2.x - n0.x,
         n1.y - n0.y, n2.y - n0.y;
    double detJ = std::abs(J.determinant());
    
    // extract current element node solution
    Eigen::Vector3d u_element;
    u_element << u_nodal[node_indices_[0]], 
                 u_nodal[node_indices_[1]], 
                 u_nodal[node_indices_[2]];
    
    // use 3-point Gauss quadrature to compute load vector
    constexpr int n_points = 3;
    auto rule = boost::math::quadrature::gauss<double, n_points>();
    
    for (int i = 0; i < n_points; ++i) {
        double xi = (rule.abscissa()[i] + 1) / 2;
        for (int j = 0; j < n_points; ++j) {
            double eta = (rule.abscissa()[j] + 1) / 2;
            
            if (xi + eta <= 1.0) {
                double weight = rule.weights()[i] * rule.weights()[j] * 0.25;
                
                // shape function values
                Eigen::Vector3d N;
                N << 1 - xi - eta, xi, eta;
                
                // interpolated solution
                double u_val = N.dot(u_element);
                
                // source term value
                double f_val = source ? source(u_val) : 0.0;
                
                    // Jacobian determinant
                double jac_det = detJ * 2.0; // transformation from [-1,1]^2 to [0,1]^2
                
                fe += weight * jac_det * f_val * N;
            }
        }
    }
    
    return fe;
}

double TriangleElement::interpolateSolution(
    const std::vector<Node>& nodes,
    const Eigen::VectorXd& u_nodal,
    double xi, double eta) const {
    
    Eigen::Vector3d N;
    N << 1 - xi - eta, xi, eta;
    
    double u_val = 0.0;
    for (int i = 0; i < 3; ++i) {
        u_val += N(i) * u_nodal[node_indices_[i]];
    }
    
    return u_val;
}

std::pair<double, double> TriangleElement::transformToGlobal(
    const std::vector<Node>& nodes,
    double xi, double eta) const {
    
    const Node& n0 = nodes[node_indices_[0]];
    const Node& n1 = nodes[node_indices_[1]];
    const Node& n2 = nodes[node_indices_[2]];
    
    double x = n0.x * (1 - xi - eta) + n1.x * xi + n2.x * eta;
    double y = n0.y * (1 - xi - eta) + n1.y * xi + n2.y * eta;
    
    return {x, y};
}

double TriangleElement::jacobianDeterminant(
    const std::vector<Node>& nodes,
    double xi, double eta) const {
    
    const Node& n0 = nodes[node_indices_[0]];
    const Node& n1 = nodes[node_indices_[1]];
    const Node& n2 = nodes[node_indices_[2]];
    
    Eigen::Matrix2d J;
    J << n1.x - n0.x, n2.x - n0.x,
         n1.y - n0.y, n2.y - n0.y;
    
    return std::abs(J.determinant());
}