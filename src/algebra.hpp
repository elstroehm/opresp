#pragma once

#include <vector>
#include <functional>
#include <iostream>
#include <utility>
#include <random>

#include "../external/eigen-3.4.0/Eigen/Core"
#include "../external/eigen-3.4.0/Eigen/Sparse"

#include "definitions.hpp"

using namespace Eigen;



namespace QSim {


    class TupGen {
    public:
        _ivec lim;
        _ivec state;
        
        TupGen(const _ivec& _lim);

        int operator[](int l);
        int get_nbins();
        _ivec get_state();
        void print_tuple();

        bool is_initial();
        bool is_final();
        void init_initial();
        void init_final();
        void operator++();
        void operator--();
    };

    MatrixXcd liouville_qform_dcomplex_dense(int nlevels);
    MatrixXd liouville_qform_double_dense(int nlevels);
    SparseMatrix<_cd> liouville_qform_dcomplex_sparse(int nlevels);
    SparseMatrix<double> liouville_qform_double_sparse(int nlevels);

    MatrixXcd liouville_lie_algebra_operator_dcomplex_dense(const MatrixXcd& op, double ep);
    MatrixXd liouville_lie_algebra_operator_double_dense(const MatrixXd& op, double ep);
    SparseMatrix<_cd> liouville_lie_algebra_operator_dcomplex_sparse(const MatrixXcd& op, double ep);
    SparseMatrix<double> liouville_lie_algebra_operator_double_sparse(const MatrixXd& op, double ep);

    MatrixXcd liouville_lie_group_operator_dcomplex_dense(const MatrixXcd& op, const MatrixXcd& op_inv, double ep);
    MatrixXd liouville_lie_group_operator_double_dense(const MatrixXd& op, const MatrixXd& op_inv, double ep);
    SparseMatrix<_cd> liouville_lie_group_operator_dcomplex_sparse(const MatrixXcd& op, const MatrixXcd& op_inv, double ep);
    SparseMatrix<double> liouville_lie_group_operator_double_sparse(const MatrixXd& op, const MatrixXd& op_inv, double ep);

    SparseMatrix<_cd> coupled_tls_lindblad_dcomplex(int ntls, const _dvec& omegas, const std::vector<double>& dephasings, const std::vector<double>& relaxations,
            const std::vector<std::pair<int, int>>& connections, const _cdvec& couplings);
    SparseMatrix<_cd> coupled_tls_lindblad_double(int ntls, const _dvec& omegas, const std::vector<double>& dephasings, const std::vector<double>& relaxations,
            const std::vector<std::pair<int, int>>& connections, const _dvec& couplings);

    MatrixXcd coupled_tls_hamiltonian_dcomplex(int ntls, const _dvec& omegas, const std::vector<std::pair<int, int>>& connections, const _cdvec& couplings);
    MatrixXd coupled_tls_hamiltonian_double(int ntls, const _dvec& omegas, const std::vector<std::pair<int, int>>& connections, const _dvec& couplings);

    MatrixXcd coupled_tls_operator_dcomplex(int ntls, const _cdvec& entries);
    MatrixXd coupled_tls_operator_double(int ntls, const _dvec& entries);
    
    SparseMatrix<_cd> eigen_sparse_from_dense_dcomplex(const MatrixXcd& m, double epsilon);
    SparseMatrix<double> eigen_sparse_from_dense_double(const MatrixXd& m, double epsilon);
    
    void eigen_prune_dcomplex(MatrixXcd& m, double epsilon);
    void eigen_prune_double(MatrixXd& m, double epsilon);

    Matrix2cd random_diagonal_state_2();
    Matrix2cd random_pure_state_2();
    MatrixXcd random_diagonal_state_n(int n);
    MatrixXcd random_pure_state_n(int n);


    
    
    // MatrixXcd rho = random_pure_state_n(dim);
    // MatrixXcd rho1 = transform*rho*transform.transpose();
    // MatrixXcd rho2 = (l_transform*rho.reshaped<RowMajor>()).reshaped<RowMajor>(dim, dim);

    // std::cout << "\n\nAusgangsmatrix:\n" << std::endl;
    // std::cout << rho << std::endl;
    // std::cout << "\n\nTransformiert mit O:\n" << std::endl;
    // std::cout << rho1 << std::endl;
    // std::cout << "\n\nTransformiert mit T:\n" << std::endl;
    // std::cout << rho2 << std::endl;
    // std::cout << "\n\nDifferenz:\n" << std::endl;
    // std::cout << rho1 - rho2 << std::endl;
}
