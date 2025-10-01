#include "VTKWriter.hpp"
#include "elements/Element.hpp"
#include "elements/TriangleElement.hpp"
#include "elements/QuadElement.hpp"
#include <fstream>
#include <iomanip>

void VTKWriter::writeSolution(const std::string& filename,
                             const std::vector<Node>& nodes,
                             const std::vector<std::unique_ptr<Element>>& elements,
                             const Eigen::VectorXd& solution) {
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    writeVTKHeader(file);
    writePoints(file, nodes);
    writeCells(file, elements);
    writePointData(file, nodes, solution);
    writeCellData(file, elements, nodes, solution);
    
    file.close();
}

void VTKWriter::writeVTKHeader(std::ofstream& file) {
    file << "# vtk DataFile Version 3.0\n";
    file << "2D Nonlinear Poisson Solution\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";
}

void VTKWriter::writePoints(std::ofstream& file, const std::vector<Node>& nodes) {
    file << "POINTS " << nodes.size() << " double\n";
    file << std::scientific << std::setprecision(15);
    
    for (const auto& node : nodes) {
        file << node.x << " " << node.y << " 0.0\n";
    }
    file << "\n";
}

void VTKWriter::writeCells(std::ofstream& file, 
                          const std::vector<std::unique_ptr<Element>>& elements) {
    
    // compute total data size
    size_t total_size = 0;
    for (const auto& element : elements) {
        total_size += 1 + element->getNodeCount(); // cell_type + node_count + node_indices
    }
    
    file << "CELLS " << elements.size() << " " << total_size << "\n";
    
    for (const auto& element : elements) {
        auto node_indices = element->getNodeIndices();
        file << node_indices.size();
        for (int idx : node_indices) {
            file << " " << idx;
        }
        file << "\n";
    }
    file << "\n";
    
    // write cell types
    file << "CELL_TYPES " << elements.size() << "\n";
    for (const auto& element : elements) {
        int vtk_cell_type = getVTKCellType(*element);
        file << vtk_cell_type << "\n";
    }
    file << "\n";
}

void VTKWriter::writePointData(std::ofstream& file, 
                              const std::vector<Node>& nodes,
                              const Eigen::VectorXd& solution) {
    
    file << "POINT_DATA " << nodes.size() << "\n";
    file << "SCALARS solution double 1\n";
    file << "LOOKUP_TABLE default\n";
    
    file << std::scientific << std::setprecision(15);
    for (int i = 0; i < nodes.size(); ++i) {
        file << solution(i) << "\n";
    }
    file << "\n";
    
    // write boundary condition information
    file << "SCALARS boundary_condition double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (const auto& node : nodes) {
        file << (node.has_bc ? node.bc_value : -1.0) << "\n";
    }
    file << "\n";
}

void VTKWriter::writeCellData(std::ofstream& file,
                             const std::vector<std::unique_ptr<Element>>& elements,
                             const std::vector<Node>& nodes,
                             const Eigen::VectorXd& solution) {
    
    file << "CELL_DATA " << elements.size() << "\n";
    file << "SCALARS element_id int 1\n";
    file << "LOOKUP_TABLE default\n";
    
    for (int i = 0; i < elements.size(); ++i) {
        file << i << "\n";
    }
    file << "\n";
}

int VTKWriter::getVTKCellType(const Element& element) {
    if (dynamic_cast<const TriangleElement*>(&element)) {
        return 5; // VTK_TRIANGLE
    } else if (dynamic_cast<const QuadElement*>(&element)) {
        return 9; // VTK_QUAD
    }
    throw std::runtime_error("Unknown element type");
}

void VTKWriter::writeLineData(const std::string& filename,
                             const std::vector<Node>& nodes,
                             const Eigen::VectorXd& solution,
                             double y_slice) {
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    file << "x,u\n";
    file << std::scientific << std::setprecision(15);
    
    // find nodes near y=y_slice
    const double tolerance = 1e-10;
    for (const auto& node : nodes) {
        if (std::abs(node.y - y_slice) < tolerance) {
            file << node.x << "," << solution(node.id) << "\n";
        }
    }
    
    file.close();
}