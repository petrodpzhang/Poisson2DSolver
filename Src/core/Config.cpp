#include "Config.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <regex>

Config Config::fromJSON(const std::string& filename) {
    // simplified JSON parsing
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open config file: " + filename);
    }
    
    Config config;
    std::string line;
    
    // Parse boundary conditions separately (track state)
    bool in_bc_section = false;
    BoundaryCondition current_bc;
    bool has_bc_type = false, has_bc_boundary = false, has_bc_value = false;
    
    while (std::getline(file, line)) {
        // Check if we're entering/exiting boundary conditions array
        if (line.find("boundary_conditions") != std::string::npos) {
            in_bc_section = true;
            continue;
        }
        if (line.find("functions") != std::string::npos || 
            line.find("solver") != std::string::npos ||
            line.find("output") != std::string::npos) {
            in_bc_section = false;
        }
        
        line = std::regex_replace(line, std::regex("[\\{\}\",]"), "");
        std::istringstream iss(line);
        std::string key, value;
        
        if (std::getline(iss, key, ':') && std::getline(iss, value)) {
            key = trim(key);
            value = trim(value);
            
            if (in_bc_section) {
                // Parse boundary condition fields
                if (key == "type") {
                    current_bc.type = value;
                    has_bc_type = true;
                }
                else if (key == "boundary") {
                    current_bc.boundary = value;
                    has_bc_boundary = true;
                }
                else if (key == "value") {
                    current_bc.value = std::stod(value);
                    has_bc_value = true;
                }
                
                // If we have all three fields, save the BC
                if (has_bc_type && has_bc_boundary && has_bc_value) {
                    config.boundary_conditions.push_back(current_bc);
                    has_bc_type = has_bc_boundary = has_bc_value = false;
                    current_bc = BoundaryCondition();
                }
            } else {
                parseConfigValue(config, key, value);
            }
        }
    }
    
    // Add default boundary conditions if none were parsed (for testing)
    if (config.boundary_conditions.empty()) {
        std::cerr << "Warning: No boundary conditions found in JSON. Adding defaults for testing." << std::endl;
        config.boundary_conditions.push_back({"Dirichlet", "AD", 2.0});
        config.boundary_conditions.push_back({"Dirichlet", "BC", 1.0});
        config.boundary_conditions.push_back({"Dirichlet", "AB", 0.0});
        config.boundary_conditions.push_back({"Dirichlet", "CD", 0.0});
    }
    
    return config;
}

Config Config::fromTXT(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open config file: " + filename);
    }
    
    Config config;
    std::string line;
    
    while (std::getline(file, line)) {
        // skip comments and empty lines
        if (line.empty() || line[0] == '#') continue;
        
        std::istringstream iss(line);
        std::string key, equals, value;
        
        if (iss >> key >> equals >> value) {
            if (equals == "=") {
                parseConfigValue(config, key, value);
            }
        }
    }
    
    return config;
}

void Config::parseConfigValue(Config& config, const std::string& key, const std::string& value) {
    // Skip empty values
    if (value.empty()) return;
    
    // Geometry
    if (key == "lx") config.lx = std::stod(value);
    else if (key == "ly") config.ly = std::stod(value);
    else if (key == "nx") config.nx = std::stoi(value);
    else if (key == "ny") config.ny = std::stoi(value);
    else if (key == "element_type") config.mesh_type = value;
    
    // Functions
    else if (key == "initial_guess") config.guess_function = value;
    else if (key == "source") config.source_function = value;
    else if (key == "source_derivative") config.source_derivative_function = value;
    
    // Solver parameters
    else if (key == "relative_tolerance") config.rel_tol = std::stod(value);
    else if (key == "absolute_tolerance") config.abs_tol = std::stod(value);
    else if (key == "max_iterations") config.max_iterations = std::stoi(value);
    
    // Output
    else if (key == "path") config.output_path = value;
    
    // Legacy keys for backward compatibility
    else if (key == "Nx") config.nx = std::stoi(value);
    else if (key == "Ny") config.ny = std::stoi(value);
    else if (key == "mesh") config.mesh_type = value;
    else if (key == "guess") config.guess_function = value;
    else if (key == "rel_tol") config.rel_tol = std::stod(value);
    else if (key == "abs_tol") config.abs_tol = std::stod(value);
    else if (key == "output_path") config.output_path = value;
}

std::string Config::trim(const std::string& str) {
    size_t start = str.find_first_not_of(" \t\n\r");
    size_t end = str.find_last_not_of(" \t\n\r");
    return (start == std::string::npos) ? "" : str.substr(start, end - start + 1);
}