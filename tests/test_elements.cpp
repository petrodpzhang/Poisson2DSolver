#include <gtest/gtest.h>
#include "elements/TriangleElement.hpp"
#include "elements/QuadElement.hpp"
#include "elements/Element.hpp"
#include <cmath>

/**
 * @brief Test triangle element shape functions at nodes
 * 
 * Shape functions should satisfy N_i(node_j) = delta_ij
 */
TEST(ElementTest, TriangleShapeFunctionsAtNodes) {
    std::vector<Node> nodes = {
        {0, 0.0, 0.0},
        {1, 1.0, 0.0},
        {2, 0.0, 1.0}
    };
    
    TriangleElement element({0, 1, 2});
    
    // At node 0 (0,0): N0=1, N1=0, N2=0
    Eigen::VectorXd u_nodal(3);
    u_nodal << 1.0, 0.0, 0.0;
    double val0 = element.interpolateSolution(nodes, u_nodal, 0.0, 0.0);
    EXPECT_NEAR(val0, 1.0, 1e-10);
    
    // At node 1 (1,0): N0=0, N1=1, N2=0
    u_nodal << 0.0, 1.0, 0.0;
    double val1 = element.interpolateSolution(nodes, u_nodal, 1.0, 0.0);
    EXPECT_NEAR(val1, 1.0, 1e-10);
    
    // At node 2 (0,1): N0=0, N1=0, N2=1
    u_nodal << 0.0, 0.0, 1.0;
    double val2 = element.interpolateSolution(nodes, u_nodal, 0.0, 1.0);
    EXPECT_NEAR(val2, 1.0, 1e-10);
}

/**
 * @brief Test quadrilateral element shape functions at nodes
 */
TEST(ElementTest, QuadShapeFunctionsAtNodes) {
    std::vector<Node> nodes = {
        {0, 0.0, 0.0},
        {1, 1.0, 0.0},
        {2, 1.0, 1.0},
        {3, 0.0, 1.0}
    };
    
    QuadElement element({0, 1, 2, 3});
    
    // Test at all four nodes in isoparametric coordinates
    std::vector<std::pair<double, double>> iso_coords = {
        {-1.0, -1.0},  // Node 0
        { 1.0, -1.0},  // Node 1
        { 1.0,  1.0},  // Node 2
        {-1.0,  1.0}   // Node 3
    };
    
    for (int i = 0; i < 4; ++i) {
        Eigen::VectorXd u_nodal = Eigen::VectorXd::Zero(4);
        u_nodal(i) = 1.0;
        
        auto [xi, eta] = iso_coords[i];
        double val = element.interpolateSolution(nodes, u_nodal, xi, eta);
        EXPECT_NEAR(val, 1.0, 1e-10) << "Shape function " << i << " should be 1 at node " << i;
    }
}

/**
 * @brief Test triangle element Jacobian determinant
 */
TEST(ElementTest, TriangleJacobianDeterminant) {
    // Right triangle with vertices at (0,0), (1,0), (0,1)
    std::vector<Node> nodes = {
        {0, 0.0, 0.0},
        {1, 1.0, 0.0},
        {2, 0.0, 1.0}
    };
    
    TriangleElement element({0, 1, 2});
    
    // Area of triangle = 0.5
    // Jacobian determinant should be 2 * area = 1.0
    double jac = element.jacobianDeterminant(nodes, 0.33, 0.33);
    EXPECT_NEAR(jac, 1.0, 1e-10);
}

/**
 * @brief Test quad element Jacobian for unit square
 */
TEST(ElementTest, QuadJacobianUnitSquare) {
    // Unit square
    std::vector<Node> nodes = {
        {0, 0.0, 0.0},
        {1, 1.0, 0.0},
        {2, 1.0, 1.0},
        {3, 0.0, 1.0}
    };
    
    QuadElement element({0, 1, 2, 3});
    
    // For unit square in isoparametric space, Jacobian should be constant
    double jac_center = element.jacobianDeterminant(nodes, 0.0, 0.0);
    double jac_corner = element.jacobianDeterminant(nodes, -1.0, -1.0);
    
    // Jacobian should be 0.25 (half of element area in iso space)
    EXPECT_NEAR(jac_center, 0.25, 1e-10);
    EXPECT_NEAR(jac_corner, 0.25, 1e-10);
}

/**
 * @brief Test stiffness matrix symmetry for triangle
 */
TEST(ElementTest, TriangleStiffnessSymmetry) {
    std::vector<Node> nodes = {
        {0, 0.0, 0.0},
        {1, 1.0, 0.0},
        {2, 0.0, 1.0}
    };
    
    TriangleElement element({0, 1, 2});
    
    // Use a simple source derivative
    auto source_deriv = [](double u) { return -1.0; };
    
    auto K = element.computeStiffnessMatrix(nodes, 0.0, source_deriv);
    
    // Check symmetry
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            EXPECT_NEAR(K(i, j), K(j, i), 1e-10) << "Matrix not symmetric at (" << i << "," << j << ")";
        }
    }
}

/**
 * @brief Test stiffness matrix symmetry for quad
 */
TEST(ElementTest, QuadStiffnessSymmetry) {
    std::vector<Node> nodes = {
        {0, 0.0, 0.0},
        {1, 1.0, 0.0},
        {2, 1.0, 1.0},
        {3, 0.0, 1.0}
    };
    
    QuadElement element({0, 1, 2, 3});
    
    // Use a simple source derivative
    auto source_deriv = [](double u) { return -1.0; };
    
    auto K = element.computeStiffnessMatrix(nodes, 0.0, source_deriv);
    
    // Check symmetry
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            EXPECT_NEAR(K(i, j), K(j, i), 1e-10) << "Matrix not symmetric at (" << i << "," << j << ")";
        }
    }
}

/**
 * @brief Test coordinate transformation from isoparametric to global
 */
TEST(ElementTest, QuadCoordinateTransformation) {
    std::vector<Node> nodes = {
        {0, 0.0, 0.0},
        {1, 2.0, 0.0},
        {2, 2.0, 1.0},
        {3, 0.0, 1.0}
    };
    
    QuadElement element({0, 1, 2, 3});
    
    // Center of isoparametric space (0, 0) should map to center of element (1, 0.5)
    auto [x, y] = element.transformToGlobal(nodes, 0.0, 0.0);
    EXPECT_NEAR(x, 1.0, 1e-10);
    EXPECT_NEAR(y, 0.5, 1e-10);
    
    // Corner transformations
    auto [x0, y0] = element.transformToGlobal(nodes, -1.0, -1.0);
    EXPECT_NEAR(x0, 0.0, 1e-10);
    EXPECT_NEAR(y0, 0.0, 1e-10);
    
    auto [x1, y1] = element.transformToGlobal(nodes, 1.0, -1.0);
    EXPECT_NEAR(x1, 2.0, 1e-10);
    EXPECT_NEAR(y1, 0.0, 1e-10);
}

/**
 * @brief Test load vector computation for constant source
 */
TEST(ElementTest, TriangleLoadVectorConstantSource) {
    std::vector<Node> nodes = {
        {0, 0.0, 0.0},
        {1, 1.0, 0.0},
        {2, 0.0, 1.0}
    };
    
    TriangleElement element({0, 1, 2});
    
    // Constant source f = 1.0
    auto source = [](double u) { return 1.0; };
    
    Eigen::VectorXd u_nodal = Eigen::VectorXd::Zero(3);
    auto f = element.computeLoadVector(nodes, u_nodal, source);
    
    // For constant source over triangle, integral should be area * f
    // Area = 0.5, so total integral = 0.5
    // Distributed equally to 3 nodes: ~0.167 each
    double total = f.sum();
    EXPECT_NEAR(total, 0.5, 0.1) << "Total load should equal area times source";
}

/**
 * @brief Test that stiffness matrix has correct properties
 */
TEST(ElementTest, StiffnessMatrixProperties) {
    std::vector<Node> nodes = {
        {0, 0.0, 0.0},
        {1, 1.0, 0.0},
        {2, 1.0, 1.0},
        {3, 0.0, 1.0}
    };
    
    QuadElement element({0, 1, 2, 3});
    
    auto source_deriv = [](double u) { return 0.0; };
    
    auto K = element.computeStiffnessMatrix(nodes, 0.0, source_deriv);
    
    // Check that matrix has non-zero norm
    EXPECT_GT(K.norm(), 0.0);
    
    // For Laplacian, row sums should be close to zero (compatibility condition)
    Eigen::VectorXd row_sums = K.rowwise().sum();
    for (int i = 0; i < 4; ++i) {
        EXPECT_NEAR(row_sums(i), 0.0, 1e-8) << "Row " << i << " sum should be near zero";
    }
}

/**
 * @brief Test element with degenerate geometry (should not crash)
 */
TEST(ElementTest, DegenerateTriangle) {
    // Collinear points (degenerate triangle)
    std::vector<Node> nodes = {
        {0, 0.0, 0.0},
        {1, 1.0, 0.0},
        {2, 2.0, 0.0}  // All on same line
    };
    
    TriangleElement element({0, 1, 2});
    
    // Jacobian should be near zero (or zero)
    double jac = element.jacobianDeterminant(nodes, 0.33, 0.33);
    EXPECT_NEAR(jac, 0.0, 1e-6);
}

/**
 * @brief Test partition of unity for shape functions
 */
TEST(ElementTest, PartitionOfUnity) {
    std::vector<Node> nodes = {
        {0, 0.0, 0.0},
        {1, 1.0, 0.0},
        {2, 1.0, 1.0},
        {3, 0.0, 1.0}
    };
    
    QuadElement element({0, 1, 2, 3});
    
    // Shape functions should sum to 1 at any point
    Eigen::VectorXd ones = Eigen::VectorXd::Ones(4);
    
    // Test at several points
    std::vector<std::pair<double, double>> test_points = {
        {0.0, 0.0}, {0.5, 0.5}, {-0.7, 0.3}, {0.8, -0.8}
    };
    
    for (const auto& [xi, eta] : test_points) {
        double sum = element.interpolateSolution(nodes, ones, xi, eta);
        EXPECT_NEAR(sum, 1.0, 1e-10) << "Partition of unity failed at (" << xi << ", " << eta << ")";
    }
}

