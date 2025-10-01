#include "LinearSystem.hpp"
#include <iostream>

#ifdef USE_SPARSE_MATRIX
LinearSystem::LinearSystem(int size) {
    stiffness_matrix_.resize(size, size);
    load_vector_.resize(size);
    solution_.resize(size);
    
    // preallocate sparse matrix memory (estimated max 9 non-zero elements per row)
    stiffness_matrix_.reserve(Eigen::VectorXi::Constant(size, 9));
}
#else
LinearSystem::LinearSystem(int size) {
    stiffness_matrix_.resize(size, size);
    load_vector_.resize(size);
    solution_.resize(size);
    
    stiffness_matrix_.setZero();
    load_vector_.setZero();
    solution_.setZero();
}
#endif

void LinearSystem::assemble(const Eigen::MatrixXd& element_matrix,
                           const Eigen::VectorXd& element_vector,
                           const std::vector<int>& dof_indices) {
    
    int n_dofs = dof_indices.size();
    
    for (int i = 0; i < n_dofs; ++i) {
        int global_i = dof_indices[i];
        load_vector_(global_i) += element_vector(i);
        
        for (int j = 0; j < n_dofs; ++j) {
            int global_j = dof_indices[j];
#ifdef USE_SPARSE_MATRIX
            stiffness_matrix_.coeffRef(global_i, global_j) += element_matrix(i, j);
#else
            stiffness_matrix_(global_i, global_j) += element_matrix(i, j);
#endif
        }
    }
}

bool LinearSystem::solve() {
#ifdef USE_SPARSE_MATRIX
    // compress sparse matrix
    stiffness_matrix_.makeCompressed();
    
    // use BiCGSTAB solver
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
    solver.compute(stiffness_matrix_);
    
    if (solver.info() != Eigen::Success) {
        std::cerr << "Sparse matrix decomposition failed" << std::endl;
        return false;
    }
    
    solution_ = solver.solve(load_vector_);
    
    // Compute residual
    residual_ = load_vector_ - stiffness_matrix_ * solution_;
    
    return solver.info() == Eigen::Success;
#else
    // use dense matrix solver
    solution_ = stiffness_matrix_.colPivHouseholderQr().solve(load_vector_);
    
    // Compute residual
    residual_ = load_vector_ - stiffness_matrix_ * solution_;
    
    return true;
#endif
}

void LinearSystem::clear() {
#ifdef USE_SPARSE_MATRIX
    stiffness_matrix_.setZero();
    stiffness_matrix_.data().squeeze();  // release unused memory
#else
    stiffness_matrix_.setZero();
#endif
    load_vector_.setZero();
    solution_.setZero();
}

Eigen::VectorXd LinearSystem::computeResidual() const {
#ifdef USE_SPARSE_MATRIX
    return stiffness_matrix_ * solution_ - load_vector_;
#else
    return stiffness_matrix_ * solution_ - load_vector_;
#endif
}

void LinearSystem::applyDirichletBC(int dof_index, double value) {
    // handle stiffness matrix
#ifdef USE_SPARSE_MATRIX
    for (Eigen::SparseMatrix<double>::InnerIterator it(stiffness_matrix_, dof_index); it; ++it) {
        if (it.row() == dof_index) {
            it.valueRef() = 1.0;
        } else {
            it.valueRef() = 0.0;
        }
    }
#else
    stiffness_matrix_.row(dof_index).setZero();
    stiffness_matrix_.col(dof_index).setZero();
    stiffness_matrix_(dof_index, dof_index) = 1.0;
#endif
    
    // handle load vector
    load_vector_(dof_index) = value;
}