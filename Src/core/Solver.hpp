#pragma once
#include "Mesh.hpp"
#include "LinearSystem.hpp"
#include "Config.hpp"

class NonlinearPoissonSolver {
private:
    Config config_;
    std::unique_ptr<Mesh> mesh_;
    std::unique_ptr<LinearSystem> linear_system_;
    
    // function parser
    std::function<double(double, double)> guess_func_;
    std::function<double(double)> source_func_;
    std::function<double(double)> source_deriv_func_;
    
    Eigen::VectorXd solution_;
    
public:
    NonlinearPoissonSolver(const Config& config);
    
    bool solve();
    void outputResults(const std::string& filename) const;
    
    const Eigen::VectorXd& getSolution() const { return solution_; }
    
private:
    void initializeFunctions();
    void assembleSystem();
    void applyBoundaryConditions();
    bool checkConvergence(const Eigen::VectorXd& delta_u, 
                         const Eigen::VectorXd& residual, 
                         int iteration) const;
    
    Eigen::VectorXd computeInitialGuess() const;
};