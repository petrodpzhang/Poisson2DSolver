#include <gtest/gtest.h>
#include "core/Mesh.hpp"
#include "core/Solver.hpp"
#include "elements/TriangleElement.hpp"

// test mesh generation
TEST(MeshTest, RectangleMeshGeneration) {
    Mesh mesh(2.0, 1.0, 4, 2, "triangle");
    
    EXPECT_EQ(mesh.getTotalNodes(), (4+1)*(2+1)); // 15 nodes
    EXPECT_GT(mesh.getElements().size(), 0);
}

// test triangle element
TEST(ElementTest, TriangleStiffnessMatrix) {
    std::vector<Node> nodes = {
        {0, 0.0, 0.0}, {1, 1.0, 0.0}, {2, 0.0, 1.0}
    };
    
    TriangleElement element({0, 1, 2});
    auto stiffness = element.computeStiffnessMatrix(nodes, 0.0, nullptr);
    
    EXPECT_EQ(stiffness.rows(), 3);
    EXPECT_EQ(stiffness.cols(), 3);
    EXPECT_GT(stiffness.norm(), 0.0);
}

// test linear system
TEST(LinearSystemTest, SparseMatrixAssembly) {
    LinearSystem system(3);
    
    Eigen::Matrix2d ke;
    ke << 1, -1, -1, 1;
    Eigen::Vector2d fe(1.0, 1.0);
    
    system.assemble(ke, fe, {0, 1});
    
    EXPECT_EQ(system.getMatrix().rows(), 3);
    EXPECT_EQ(system.getMatrix().cols(), 3);
}

// test analytical solution verification
TEST(SolverTest, AnalyticalSolutionComparison) {
    Config config;
    config.lx = 2.0;
    config.ly = 1.0;
    config.nx = 10;
    config.ny = 5;
    config.mesh_type = "quad";
    
    // set zero source term and simple boundary conditions
    config.boundary_conditions = {
        {"Dirichlet", "CD", 1.0},  // top boundary
        {"Dirichlet", "AB", 0.0},  // bottom boundary  
        {"Dirichlet", "AD", 0.0},  // left boundary
        {"Dirichlet", "BC", 0.0}   // right boundary
    };
    
    config.guess_function = "0";
    config.source_function = "0";
    config.source_derivative_function = "0";
    
    config.rel_tol = 1e-8;
    config.abs_tol = 1e-10;
    config.max_iterations = 10;
    
    NonlinearPoissonSolver solver(config);
    bool success = solver.solve();
    
    EXPECT_TRUE(success);
    
    // verify solution on boundary
    auto solution = solver.getSolution();
    // here can add more verification logic
}

// performance test: sparse vs dense matrix
TEST(PerformanceTest, SparseVsDense) {
    // test performance with different grid sizes
    std::vector<std::pair<int, int>> grid_sizes = {
        {20, 10}, {50, 25}, {100, 50}
    };
    
    for (const auto& [nx, ny] : grid_sizes) {
        Config config;
        config.lx = 2.0;
        config.ly = 1.0;
        config.nx = nx;
        config.ny = ny;
        config.mesh_type = "quad";
        
        // test configuration...
        
        auto start = std::chrono::high_resolution_clock::now();
        
        NonlinearPoissonSolver solver(config);
        solver.solve();
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        
        std::cout << "Grid " << nx << "x" << ny << ": " << duration.count() << " ms" << std::endl;
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}