#include "Solver.hpp"
#include "utils/Parser.hpp"
#include "utils/VTKWriter.hpp"
#include <iostream>
#include <iomanip>

NonlinearPoissonSolver::NonlinearPoissonSolver(const Config& config)
    : config_(config) {
    
    // create mesh
    mesh_ = std::make_unique<Mesh>(config.lx, config.ly, config.nx, config.ny, config.mesh_type);
    
    // apply boundary conditions
    mesh_->applyBoundaryConditions(config.boundary_conditions);
    
    // initialize linear system
    int total_nodes = mesh_->getTotalNodes();
    linear_system_ = std::make_unique<LinearSystem>(total_nodes);
    
    // initialize function parser
    initializeFunctions();
}

bool NonlinearPoissonSolver::solve() {
    std::cout << "Starting Newton iteration..." << std::endl;
    std::cout << std::setw(6) << "Step" 
              << std::setw(15) << "Abs. Error" 
              << std::setw(15) << "Rel. Error" << std::endl;
    std::cout << std::string(36, '-') << std::endl;
    
    solution_ = computeInitialGuess();
    
    for (int iter = 0; iter < config_.max_iterations; ++iter) {
        // assemble system
        assembleSystem();
        
        // solve linear system
        if (!linear_system_->solve()) {
            std::cerr << "Linear system solve failed at iteration " << iter << std::endl;
            return false;
        }
        
        Eigen::VectorXd delta_u = linear_system_->getSolution();
        Eigen::VectorXd residual = linear_system_->getResidual();
        
        // update solution
        solution_ += delta_u;
        
        // check convergence
        if (checkConvergence(delta_u, residual, iter)) {
            std::cout << "Converged after " << iter + 1 << " iterations" << std::endl;
            return true;
        }
    }
    
    std::cout << "Failed to converge within " << config_.max_iterations << " iterations" << std::endl;
    return false;
}

void NonlinearPoissonSolver::assembleSystem() {
    linear_system_->clear();
    
    const auto& nodes = mesh_->getNodes();
    const auto& elements = mesh_->getElements();
    
    for (const auto& element : elements) {
        // compute element average solution (for linearization)
        double u_avg = 0.0;
        auto node_indices = element->getNodeIndices();
        for (int node_idx : node_indices) {
            u_avg += solution_(node_idx);
        }
        u_avg /= node_indices.size();
        
        // compute element stiffness matrix and load vector
        auto ke = element->computeStiffnessMatrix(nodes, u_avg, source_deriv_func_);
        auto fe = element->computeLoadVector(nodes, solution_, source_func_);
        
        // assemble to global system
        linear_system_->assemble(ke, fe, node_indices);
    }
    
    // apply boundary conditions
    applyBoundaryConditions();
}

void NonlinearPoissonSolver::applyBoundaryConditions() {
    const auto& nodes = mesh_->getNodes();
    
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (nodes[i].has_bc) {
            linear_system_->applyDirichletBC(i, nodes[i].bc_value);
        }
    }
}

bool NonlinearPoissonSolver::checkConvergence(const Eigen::VectorXd& delta_u, 
                                             const Eigen::VectorXd& residual,
                                             int iteration) const {
    
    double delta_norm = delta_u.norm();
    double u_norm = solution_.norm();
    double residual_norm = residual.norm();
    
    double rel_error = (u_norm > 1e-12) ? delta_norm / u_norm : delta_norm;
    
    std::cout << std::setw(6) << iteration
              << std::setw(15) << std::scientific << residual_norm
              << std::setw(15) << rel_error << std::endl;
    
    return (residual_norm < config_.abs_tol) && (rel_error < config_.rel_tol);
}

void NonlinearPoissonSolver::initializeFunctions() {
    // parse user-defined functions
    guess_func_ = Parser::parse2D(config_.guess_function);
    source_func_ = Parser::parse1D(config_.source_function);
    source_deriv_func_ = Parser::parse1D(config_.source_derivative_function);
}

Eigen::VectorXd NonlinearPoissonSolver::computeInitialGuess() const {
    const auto& nodes = mesh_->getNodes();
    Eigen::VectorXd initial_guess(nodes.size());
    
    for (size_t i = 0; i < nodes.size(); ++i) {
        initial_guess(i) = guess_func_(nodes[i].x, nodes[i].y);
    }
    
    return initial_guess;
}

void NonlinearPoissonSolver::outputResults(const std::string& filename) const {
    VTKWriter::writeSolution(filename, mesh_->getNodes(), mesh_->getElements(), solution_);
}