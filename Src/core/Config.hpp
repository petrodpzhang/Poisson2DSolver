#pragma once
#include <string>
#include <vector>
#include <map>
#include <Eigen/Dense>

struct BoundaryCondition {
    std::string type;  // "Dirichlet"
    std::string boundary; // "AB", "BC", "CD", "AD"
    double value;
};

struct Config {
    // geometry and mesh
    double lx, ly;
    int nx, ny;
    std::string mesh_type; // "triangle" or "quad"
    
    // boundary conditions
    std::vector<BoundaryCondition> boundary_conditions;
    
    // function expressions
    std::string guess_function;
    std::string source_function;
    std::string source_derivative_function;
    
    // solving parameters
    double rel_tol = 1e-6;
    double abs_tol = 1e-8;
    int max_iterations = 100;
    
    // output
    std::string output_path;
    
    static Config fromJSON(const std::string& filename);
    static Config fromTXT(const std::string& filename);
    
private:
    static void parseConfigValue(Config& config, const std::string& key, const std::string& value);
    static std::string trim(const std::string& str);
};