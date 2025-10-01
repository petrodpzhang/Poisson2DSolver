#include "QuadElement.hpp"
#include <boost/math/quadrature/gauss.hpp>

QuadElement::QuadElement(const std::vector<int>& node_indices)
    : node_indices_(node_indices) {
    if (node_indices.size() != 4) {
        throw std::invalid_argument("Quad element requires exactly 4 nodes");
    }
}

Eigen::MatrixXd QuadElement::computeStiffnessMatrix(
    const std::vector<Node>& nodes,
    double u_avg,
    const std::function<double(double)>& source_deriv) const {
    
    Eigen::MatrixXd ke = Eigen::MatrixXd::Zero(4, 4);
    
    constexpr int n_points = 2; // 2x2 Gauss quadrature
    auto rule = boost::math::quadrature::gauss<double, n_points>();
    
    for (int i = 0; i < n_points; ++i) {
        double xi = rule.abscissa()[i];
        for (int j = 0; j < n_points; ++j) {
            double eta = rule.abscissa()[j];
            double weight = rule.weights()[i] * rule.weights()[j];
            
            // shape function gradient
            auto [N, dN_dxi, dN_deta] = shapeFunctions(xi, eta);
            
            // Jacobian matrix
            Eigen::Matrix2d J = Eigen::Matrix2d::Zero();
            for (int k = 0; k < 4; ++k) {
                const Node& node = nodes[node_indices_[k]];
                J(0,0) += dN_dxi[k] * node.x;
                J(0,1) += dN_dxi[k] * node.y;
                J(1,0) += dN_deta[k] * node.x;
                J(1,1) += dN_deta[k] * node.y;
            }
            
            double detJ = std::abs(J.determinant());
            Eigen::Matrix2d invJ = J.inverse();
            
            // physical coordinate system shape function gradient
            Eigen::Matrix<double, 2, 4> dN_dx;
            for (int k = 0; k < 4; ++k) {
                Eigen::Vector2d dN_ref(dN_dxi[k], dN_deta[k]);
                Eigen::Vector2d dN_global = invJ * dN_ref;
                dN_dx(0, k) = dN_global[0];
                dN_dx(1, k) = dN_global[1];
            }
            
            // assemble stiffness matrix
            for (int m = 0; m < 4; ++m) {
                for (int n = 0; n < 4; ++n) {
                    // linear part
                    ke(m, n) += weight * detJ * dN_dx.col(m).dot(dN_dx.col(n));
                    
                    // nonlinear part
                    if (source_deriv) {
                        double f_prime = source_deriv(u_avg);
                        ke(m, n) -= weight * detJ * f_prime * N[m] * N[n];
                    }
                }
            }
        }
    }
    
    return ke;
}

Eigen::VectorXd QuadElement::computeLoadVector(
    const std::vector<Node>& nodes,
    const Eigen::VectorXd& u_nodal,
    const std::function<double(double)>& source) const {
    
    Eigen::VectorXd fe = Eigen::VectorXd::Zero(4);
    
    constexpr int n_points = 2; // 2x2 Gauss quadrature
    auto rule = boost::math::quadrature::gauss<double, n_points>();
    
    // Extract current element node solution
    Eigen::Vector4d u_element;
    for (int i = 0; i < 4; ++i) {
        u_element(i) = u_nodal[node_indices_[i]];
    }
    
    for (int i = 0; i < n_points; ++i) {
        double xi = rule.abscissa()[i];
        for (int j = 0; j < n_points; ++j) {
            double eta = rule.abscissa()[j];
            double weight = rule.weights()[i] * rule.weights()[j];
            
            auto [N, dN_dxi, dN_deta] = shapeFunctions(xi, eta);
            
            // Interpolate solution at Gauss point
            double u_val = N.dot(u_element);
            
            // Jacobian
            double detJ = jacobianDeterminant(nodes, xi, eta);
            
            // Source term
            double f_val = source ? source(u_val) : 0.0;
            
            fe += weight * detJ * f_val * N;
        }
    }
    
    return fe;
}

double QuadElement::interpolateSolution(
    const std::vector<Node>& nodes,
    const Eigen::VectorXd& u_nodal,
    double xi, double eta) const {
    
    auto [N, dN_dxi, dN_deta] = shapeFunctions(xi, eta);
    
    double u_interp = 0.0;
    for (int i = 0; i < 4; ++i) {
        u_interp += N[i] * u_nodal[node_indices_[i]];
    }
    
    return u_interp;
}

std::pair<double, double> QuadElement::transformToGlobal(
    const std::vector<Node>& nodes,
    double xi, double eta) const {
    
    auto [N, dN_dxi, dN_deta] = shapeFunctions(xi, eta);
    
    double x = 0.0, y = 0.0;
    for (int i = 0; i < 4; ++i) {
        const Node& node = nodes[node_indices_[i]];
        x += N[i] * node.x;
        y += N[i] * node.y;
    }
    
    return {x, y};
}

double QuadElement::jacobianDeterminant(
    const std::vector<Node>& nodes,
    double xi, double eta) const {
    
    auto [N, dN_dxi, dN_deta] = shapeFunctions(xi, eta);
    
    Eigen::Matrix2d J = Eigen::Matrix2d::Zero();
    for (int i = 0; i < 4; ++i) {
        const Node& node = nodes[node_indices_[i]];
        J(0,0) += dN_dxi[i] * node.x;
        J(0,1) += dN_dxi[i] * node.y;
        J(1,0) += dN_deta[i] * node.x;
        J(1,1) += dN_deta[i] * node.y;
    }
    
    return std::abs(J.determinant());
}

std::tuple<Eigen::Vector4d, Eigen::Vector4d, Eigen::Vector4d> 
QuadElement::shapeFunctions(double xi, double eta) const {
    
    Eigen::Vector4d N, dN_dxi, dN_deta;
    
    N[0] = 0.25 * (1 - xi) * (1 - eta);
    N[1] = 0.25 * (1 + xi) * (1 - eta);
    N[2] = 0.25 * (1 + xi) * (1 + eta);
    N[3] = 0.25 * (1 - xi) * (1 + eta);
    
    dN_dxi[0] = -0.25 * (1 - eta);
    dN_dxi[1] =  0.25 * (1 - eta);
    dN_dxi[2] =  0.25 * (1 + eta);
    dN_dxi[3] = -0.25 * (1 + eta);
    
    dN_deta[0] = -0.25 * (1 - xi);
    dN_deta[1] = -0.25 * (1 + xi);
    dN_deta[2] =  0.25 * (1 + xi);
    dN_deta[3] =  0.25 * (1 - xi);
    
    return {N, dN_dxi, dN_deta};
}