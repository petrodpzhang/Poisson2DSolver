#include <gtest/gtest.h>
#include "core/Solver.hpp"
#include "core/Config.hpp"
#include <cmath>

/**
 * @brief Computes analytical solution for Laplace equation with specific BC
 * 
 * Implements equation 2.2 from the document:
 * u(x,y) = (4*u_hat/pi) * sum_{n=1,3,5,...} sin(n*pi*x/lx) * sinh(n*pi*y/ly) / (n * sinh(n*pi*ly/lx))
 * 
 * @param x X-coordinate
 * @param y Y-coordinate
 * @param lx Domain length in x
 * @param ly Domain length in y
 * @param u_hat Boundary condition value at CD edge
 * @param n_terms Number of terms in series (default 50 for good accuracy)
 * @return Analytical solution value at (x,y)
 */
double analyticalSolution(double x, double y, double lx, double ly, double u_hat, int n_terms = 50) {
    double sum = 0.0;
    const double pi = M_PI;
    
    for (int n = 1; n <= n_terms; n += 2) {  // odd terms only
        double factor = (4.0 * u_hat) / (pi * n);
        double sin_term = std::sin(n * pi * x / lx);
        double sinh_num = std::sinh(n * pi * y / ly);
        double sinh_den = std::sinh(n * pi * ly / lx);
        
        sum += factor * sin_term * sinh_num / sinh_den;
    }
    
    return sum;
}

/**
 * @brief Test case for Laplace equation (f=0) with analytical comparison
 * 
 * Setup: Rectangular domain with f=0, u_CD = 1.0, other boundaries = 0
 * This matches the problem in section 2.11 of the document
 */
TEST(AnalyticalTest, LaplaceSolutionAccuracy) {
    // Setup configuration for Laplace equation
    Config config;
    config.lx = 2.0;
    config.ly = 1.0;
    config.nx = 20;
    config.ny = 10;
    config.mesh_type = "quad";
    
    // Boundary conditions: CD = 1.0, others = 0
    config.boundary_conditions = {
        {"Dirichlet", "CD", 1.0},  // top boundary
        {"Dirichlet", "AB", 0.0},  // bottom boundary
        {"Dirichlet", "AD", 0.0},  // left boundary
        {"Dirichlet", "BC", 0.0}   // right boundary
    };
    
    // Zero source term (Laplace equation)
    config.guess_function = "0";
    config.source_function = "0";
    config.source_derivative_function = "0";
    
    // Tight tolerances
    config.rel_tol = 1e-8;
    config.abs_tol = 1e-10;
    config.max_iterations = 20;
    
    // Solve
    NonlinearPoissonSolver solver(config);
    bool success = solver.solve();
    
    ASSERT_TRUE(success) << "Solver failed to converge";
    
    // Get numerical solution
    auto numerical = solver.getSolution();
    
    // Compare with analytical solution at several points
    double dx = config.lx / config.nx;
    double dy = config.ly / config.ny;
    
    double max_error = 0.0;
    double l2_error_squared = 0.0;
    int num_interior_nodes = 0;
    
    for (int j = 0; j <= config.ny; ++j) {
        for (int i = 0; i <= config.nx; ++i) {
            double x = i * dx;
            double y = j * dy;
            int node_id = j * (config.nx + 1) + i;
            
            double analytical = analyticalSolution(x, y, config.lx, config.ly, 1.0);
            double numerical_val = numerical(node_id);
            double error = std::abs(numerical_val - analytical);
            
            // Skip boundary nodes for error calculation (they are exact)
            if (i > 0 && i < config.nx && j > 0 && j < config.ny) {
                max_error = std::max(max_error, error);
                l2_error_squared += error * error;
                num_interior_nodes++;
            }
        }
    }
    
    double l2_error = std::sqrt(l2_error_squared / num_interior_nodes);
    
    std::cout << "Analytical Solution Test Results:" << std::endl;
    std::cout << "  L2 Error: " << l2_error << std::endl;
    std::cout << "  Max Error: " << max_error << std::endl;
    
    // Error should be small for this mesh resolution
    EXPECT_LT(l2_error, 5e-3) << "L2 error too large";
    EXPECT_LT(max_error, 1e-2) << "Maximum pointwise error too large";
}

/**
 * @brief Test with finer mesh to verify convergence
 */
TEST(AnalyticalTest, MeshRefinementConvergence) {
    std::vector<std::pair<int, int>> grid_sizes = {{10, 5}, {20, 10}, {40, 20}};
    std::vector<double> errors;
    
    for (const auto& [nx, ny] : grid_sizes) {
        Config config;
        config.lx = 2.0;
        config.ly = 1.0;
        config.nx = nx;
        config.ny = ny;
        config.mesh_type = "quad";
        
        config.boundary_conditions = {
            {"Dirichlet", "CD", 1.0},
            {"Dirichlet", "AB", 0.0},
            {"Dirichlet", "AD", 0.0},
            {"Dirichlet", "BC", 0.0}
        };
        
        config.guess_function = "0";
        config.source_function = "0";
        config.source_derivative_function = "0";
        config.rel_tol = 1e-8;
        config.abs_tol = 1e-10;
        config.max_iterations = 20;
        
        NonlinearPoissonSolver solver(config);
        solver.solve();
        auto numerical = solver.getSolution();
        
        // Compute L2 error
        double dx = config.lx / nx;
        double dy = config.ly / ny;
        double l2_error_squared = 0.0;
        int count = 0;
        
        for (int j = 1; j < ny; ++j) {
            for (int i = 1; i < nx; ++i) {
                double x = i * dx;
                double y = j * dy;
                int node_id = j * (nx + 1) + i;
                
                double analytical = analyticalSolution(x, y, config.lx, config.ly, 1.0);
                double error = numerical(node_id) - analytical;
                l2_error_squared += error * error;
                count++;
            }
        }
        
        double l2_error = std::sqrt(l2_error_squared / count);
        errors.push_back(l2_error);
        
        std::cout << "Grid " << nx << "x" << ny << ": L2 error = " << l2_error << std::endl;
    }
    
    // Error should decrease with mesh refinement
    EXPECT_LT(errors[1], errors[0]) << "Error should decrease with refinement";
    EXPECT_LT(errors[2], errors[1]) << "Error should continue decreasing";
}

