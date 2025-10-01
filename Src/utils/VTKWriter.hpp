#pragma once
#include <string>
#include <vector>
#include <memory>
#include <Eigen/Dense>
 
struct Node;
class Element;

class VTKWriter {
public:
    static void writeSolution(const std::string& filename,
                             const std::vector<Node>& nodes,
                             const std::vector<std::unique_ptr<Element>>& elements,
                             const Eigen::VectorXd& solution);
    
    static void writeLineData(const std::string& filename,
                             const std::vector<Node>& nodes,
                             const Eigen::VectorXd& solution,
                             double y_slice);
    
private:
    static void writeVTKHeader(std::ofstream& file);
    static void writePoints(std::ofstream& file, const std::vector<Node>& nodes);
    static void writeCells(std::ofstream& file, 
                          const std::vector<std::unique_ptr<Element>>& elements);
    static void writePointData(std::ofstream& file, 
                              const std::vector<Node>& nodes,
                              const Eigen::VectorXd& solution);
    static void writeCellData(std::ofstream& file,
                             const std::vector<std::unique_ptr<Element>>& elements,
                             const std::vector<Node>& nodes,
                             const Eigen::VectorXd& solution);
    
    static int getVTKCellType(const Element& element);
};