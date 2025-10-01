#pragma once
#include <Eigen/Dense>
#include <vector>

#ifdef USE_SPARSE_MATRIX
#include <Eigen/Sparse>
#endif

class LinearSystem {
private:
#ifdef USE_SPARSE_MATRIX
    Eigen::SparseMatrix<double> stiffness_matrix_;
#else
    Eigen::MatrixXd stiffness_matrix_;
#endif
    Eigen::VectorXd load_vector_;
    Eigen::VectorXd solution_;
    Eigen::VectorXd residual_;
    
public:
    LinearSystem(int size);
    
    void clear();
    void assemble(const Eigen::MatrixXd& element_matrix,
                 const Eigen::VectorXd& element_vector,
                 const std::vector<int>& dof_indices);
    
    bool solve();
    void applyDirichletBC(int dof_index, double value);
    
    // Getters
#ifdef USE_SPARSE_MATRIX
    const Eigen::SparseMatrix<double>& getMatrix() const { return stiffness_matrix_; }
#else
    const Eigen::MatrixXd& getMatrix() const { return stiffness_matrix_; }
#endif
    
    const Eigen::VectorXd& getLoadVector() const { return load_vector_; }
    const Eigen::VectorXd& getSolution() const { return solution_; }
    const Eigen::VectorXd& getResidual() const { return residual_; }
    
    Eigen::VectorXd computeResidual() const;
};