#include <gtest/gtest.h>
#include "core/Mesh.hpp"

/**
 * @brief Test basic rectangular mesh generation with triangular elements
 */
TEST(MeshTest, TriangleMeshGeneration) {
    double lx = 2.0, ly = 1.0;
    int nx = 4, ny = 2;
    
    Mesh mesh(lx, ly, nx, ny, "triangle");
    
    // Check number of nodes
    int expected_nodes = (nx + 1) * (ny + 1);  // 5 * 3 = 15 nodes
    EXPECT_EQ(mesh.getTotalNodes(), expected_nodes);
    
    // Check number of elements (each rectangle divided into 2 triangles)
    int expected_elements = nx * ny * 2;  // 4 * 2 * 2 = 16 triangles
    EXPECT_EQ(mesh.getElements().size(), expected_elements);
    
    // Verify all elements have 3 nodes
    for (const auto& element : mesh.getElements()) {
        EXPECT_EQ(element->getNodeCount(), 3);
    }
}

/**
 * @brief Test rectangular mesh generation with quadrilateral elements
 */
TEST(MeshTest, QuadMeshGeneration) {
    double lx = 2.0, ly = 1.0;
    int nx = 4, ny = 2;
    
    Mesh mesh(lx, ly, nx, ny, "quad");
    
    // Check number of nodes
    int expected_nodes = (nx + 1) * (ny + 1);  // 5 * 3 = 15 nodes
    EXPECT_EQ(mesh.getTotalNodes(), expected_nodes);
    
    // Check number of elements
    int expected_elements = nx * ny;  // 4 * 2 = 8 quads
    EXPECT_EQ(mesh.getElements().size(), expected_elements);
    
    // Verify all elements have 4 nodes
    for (const auto& element : mesh.getElements()) {
        EXPECT_EQ(element->getNodeCount(), 4);
    }
}

/**
 * @brief Test node coordinates are correctly generated
 */
TEST(MeshTest, NodeCoordinates) {
    double lx = 2.0, ly = 1.0;
    int nx = 2, ny = 2;
    
    Mesh mesh(lx, ly, nx, ny, "quad");
    
    const auto& nodes = mesh.getNodes();
    double dx = lx / nx;
    double dy = ly / ny;
    
    // Check corner nodes
    EXPECT_DOUBLE_EQ(nodes[0].x, 0.0);
    EXPECT_DOUBLE_EQ(nodes[0].y, 0.0);
    
    EXPECT_DOUBLE_EQ(nodes[nx].x, lx);
    EXPECT_DOUBLE_EQ(nodes[nx].y, 0.0);
    
    EXPECT_DOUBLE_EQ(nodes[(ny) * (nx + 1)].x, 0.0);
    EXPECT_DOUBLE_EQ(nodes[(ny) * (nx + 1)].y, ly);
    
    EXPECT_DOUBLE_EQ(nodes[(ny) * (nx + 1) + nx].x, lx);
    EXPECT_DOUBLE_EQ(nodes[(ny) * (nx + 1) + nx].y, ly);
    
    // Check spacing
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            int node_id = j * (nx + 1) + i;
            EXPECT_NEAR(nodes[node_id].x, i * dx, 1e-10);
            EXPECT_NEAR(nodes[node_id].y, j * dy, 1e-10);
        }
    }
}

/**
 * @brief Test boundary node identification
 */
TEST(MeshTest, BoundaryIdentification) {
    double lx = 2.0, ly = 1.0;
    int nx = 4, ny = 2;
    
    Mesh mesh(lx, ly, nx, ny, "quad");
    const auto& nodes = mesh.getNodes();
    
    // Count boundary nodes
    int boundary_count = 0;
    int interior_count = 0;
    
    for (const auto& node : nodes) {
        if (node.is_boundary) {
            boundary_count++;
            // Boundary nodes should be on the edge
            EXPECT_TRUE(node.x == 0.0 || node.x == lx || 
                       node.y == 0.0 || node.y == ly);
        } else {
            interior_count++;
            // Interior nodes should not be on the edge
            EXPECT_TRUE(node.x > 0.0 && node.x < lx && 
                       node.y > 0.0 && node.y < ly);
        }
    }
    
    // Expected boundary nodes: 2*(nx+1) + 2*(ny-1) = 2*5 + 2*1 = 12
    // Expected interior nodes: (nx-1)*(ny-1) = 3*1 = 3
    EXPECT_EQ(boundary_count, 2 * (nx + 1) + 2 * (ny - 1));
    EXPECT_EQ(interior_count, (nx - 1) * (ny - 1));
}

/**
 * @brief Test boundary condition application
 */
TEST(MeshTest, BoundaryConditionApplication) {
    double lx = 2.0, ly = 1.0;
    int nx = 4, ny = 2;
    
    Mesh mesh(lx, ly, nx, ny, "quad");
    
    // Apply boundary conditions
    std::vector<BoundaryCondition> bcs = {
        {"Dirichlet", "AB", 0.0},  // bottom
        {"Dirichlet", "CD", 1.0}   // top
    };
    
    mesh.applyBoundaryConditions(bcs);
    
    const auto& nodes = mesh.getNodes();
    
    // Check that AB nodes have bc_value = 0.0
    for (const auto& node : nodes) {
        if (node.y == 0.0) {
            EXPECT_TRUE(node.has_bc);
            EXPECT_DOUBLE_EQ(node.bc_value, 0.0);
        }
    }
    
    // Check that CD nodes have bc_value = 1.0
    for (const auto& node : nodes) {
        if (node.y == ly) {
            EXPECT_TRUE(node.has_bc);
            EXPECT_DOUBLE_EQ(node.bc_value, 1.0);
        }
    }
}

/**
 * @brief Test element connectivity for triangle elements
 */
TEST(MeshTest, TriangleElementConnectivity) {
    double lx = 1.0, ly = 1.0;
    int nx = 2, ny = 2;
    
    Mesh mesh(lx, ly, nx, ny, "triangle");
    
    const auto& elements = mesh.getElements();
    
    // Each element should reference valid node indices
    int total_nodes = mesh.getTotalNodes();
    
    for (const auto& element : elements) {
        auto indices = element->getNodeIndices();
        EXPECT_EQ(indices.size(), 3);
        
        for (int idx : indices) {
            EXPECT_GE(idx, 0);
            EXPECT_LT(idx, total_nodes);
        }
        
        // No repeated nodes in an element
        EXPECT_NE(indices[0], indices[1]);
        EXPECT_NE(indices[1], indices[2]);
        EXPECT_NE(indices[0], indices[2]);
    }
}

/**
 * @brief Test element connectivity for quad elements
 */
TEST(MeshTest, QuadElementConnectivity) {
    double lx = 1.0, ly = 1.0;
    int nx = 2, ny = 2;
    
    Mesh mesh(lx, ly, nx, ny, "quad");
    
    const auto& elements = mesh.getElements();
    
    // Each element should reference valid node indices
    int total_nodes = mesh.getTotalNodes();
    
    for (const auto& element : elements) {
        auto indices = element->getNodeIndices();
        EXPECT_EQ(indices.size(), 4);
        
        for (int idx : indices) {
            EXPECT_GE(idx, 0);
            EXPECT_LT(idx, total_nodes);
        }
        
        // No repeated nodes in an element
        for (size_t i = 0; i < indices.size(); ++i) {
            for (size_t j = i + 1; j < indices.size(); ++j) {
                EXPECT_NE(indices[i], indices[j]);
            }
        }
    }
}

/**
 * @brief Test node numbering scheme
 */
TEST(MeshTest, NodeNumberingScheme) {
    double lx = 2.0, ly = 1.0;
    int nx = 4, ny = 2;
    
    Mesh mesh(lx, ly, nx, ny, "quad");
    const auto& nodes = mesh.getNodes();
    
    // Nodes should be numbered in row-major order (bottom to top, left to right)
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            int expected_id = j * (nx + 1) + i;
            EXPECT_EQ(nodes[expected_id].id, expected_id);
        }
    }
}

/**
 * @brief Test invalid element type throws exception
 */
TEST(MeshTest, InvalidElementType) {
    double lx = 2.0, ly = 1.0;
    int nx = 4, ny = 2;
    
    EXPECT_THROW({
        Mesh mesh(lx, ly, nx, ny, "invalid_type");
    }, std::invalid_argument);
}

