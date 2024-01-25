#include "algebra.hpp"

using namespace QSim;



TupGen::TupGen(const _ivec& _lim) : lim(_lim), state(_ivec(_lim.size(), 0)) {}

int TupGen::operator[](int l) {
    return state[l];
}

int TupGen::get_nbins() {
    return lim.size();
}

_ivec TupGen::get_state() {
    return state;
}

void TupGen::print_tuple() {
    for (int i = 0; i < state.size(); i++) {
        std::cout << state[i] << "  ";
    }
    std::cout << std::endl;
}


bool TupGen::is_initial() {
    for (int i = 0; i < state.size(); i++) {
        if (state[i] != 0) {
            return false;
        }
    }
    return true;
}

bool TupGen::is_final() {
    for (int i = 0; i < state.size(); i++) {
        if (state[i] != lim[i]) {
            return false;
        }
    }
    return true;
}

void TupGen::init_initial() {
    for (int i = 0; i < state.size(); i++) {
        state[i] = 0;
    }
}

void TupGen::init_final() {
    for (int i = 0; i < state.size(); i++) {
        state[i] = lim[i];
    }
}

void TupGen::operator++() {
    int carryover = 1;
    for (int i = 0; i < state.size(); i++) {
        int dummy = state[i] + carryover;
        state[i] = dummy % (lim[i] + 1);
        carryover = dummy / (lim[i] + 1);
        if (carryover == 0) {break;}
    }
}

void TupGen::operator--() {
    int carryover = 1;
    for (int i = 0; i < state.size(); i++) {
        int dummy = state[i] - carryover + (lim[i] + 1);
        state[i] = dummy % (lim[i] + 1);
        carryover = 1 - (dummy / (lim[i] + 1));
        if (carryover == 0) {break;}
    }
}


/*******************************************************************************/


MatrixXcd QSim::liouville_qform_dcomplex_dense(int nlevels) {
    int dim = nlevels*nlevels;
    MatrixXcd result = MatrixXcd::Zero(dim, dim);
    for (int i = 0; i < nlevels; i++) {
        for (int j = 0; j < nlevels; j++) {
            result(i*nlevels+j, j*nlevels+i) = {1.0, 0.0};
        }
    }
    return result;
}


MatrixXd QSim::liouville_qform_double_dense(int nlevels) {
    int dim = nlevels*nlevels;
    MatrixXd result = MatrixXd::Zero(dim, dim);
    for (int i = 0; i < nlevels; i++) {
        for (int j = 0; j < nlevels; j++) {
            result(i*nlevels+j, j*nlevels+i) = 1.0;
        }
    }
    return result;
}


SparseMatrix<_cd> QSim::liouville_qform_dcomplex_sparse(int nlevels) {
    int dim = nlevels*nlevels;
    SparseMatrix<_cd> result = SparseMatrix<_cd>(dim, dim);
    for (int i = 0; i < nlevels; i++) {
        for (int j = 0; j < nlevels; j++) {
            result.insert(i*nlevels+j, j*nlevels+i) = {1.0, 0.0};
        }
    }
    return result;
}


SparseMatrix<double> QSim::liouville_qform_double_sparse(int nlevels) {
    int dim = nlevels*nlevels;
    SparseMatrix<double> result = SparseMatrix<double>(dim, dim);
    for (int i = 0; i < nlevels; i++) {
        for (int j = 0; j < nlevels; j++) {
            result.insert(i*nlevels+j, j*nlevels+i) = 1.0;
        }
    }
    return result;
}


/*******************************************************************************/


MatrixXcd QSim::liouville_lie_algebra_operator_dcomplex_dense(const MatrixXcd& op, double ep) {
    int hdim = op.cols();
    int ldim = hdim*hdim;
    MatrixXcd op_l = MatrixXcd::Zero(ldim, ldim);
    for(int i = 0; i < hdim; i++) {
        for (int j = 0; j < hdim; j++) {
            for(int k = 0; k < hdim; k++) {
                if (std::abs(op(i,k).real()) > ep && std::abs(op(i,k).imag()) > ep) {op_l(i*hdim+j,k*hdim+j) += op(i,k);}
                if (std::abs(op(k,j).real()) > ep && std::abs(op(k,j).imag()) > ep) {op_l(i*hdim+j,i*hdim+k) -= op(k,j);}
            }
        }
    }
    return op_l;
}


MatrixXd QSim::liouville_lie_algebra_operator_double_dense(const MatrixXd& op, double ep) {
    int hdim = op.cols();
    int ldim = hdim*hdim;
    MatrixXd op_l = MatrixXd::Zero(ldim, ldim);
    for(int i = 0; i < hdim; i++) {
        for (int j = 0; j < hdim; j++) {
            for(int k = 0; k < hdim; k++) {
                if (op(i,k) > ep) {op_l(i*hdim+j,k*hdim+j) += op(i,k);}
                if (op(k,j) > ep) {op_l(i*hdim+j,i*hdim+k) -= op(k,j);}
            }
        }
    }
    return op_l;
}


SparseMatrix<_cd> QSim::liouville_lie_algebra_operator_dcomplex_sparse(const MatrixXcd& op, double ep) {
    int hdim = op.cols();
    int ldim = hdim*hdim;
    SparseMatrix<_cd> op_l = SparseMatrix<_cd>(ldim, ldim);
    for(int i = 0; i < hdim; i++) {
        for (int j = 0; j < hdim; j++) {
            for(int k = 0; k < hdim; k++) {
                if (std::abs(op(i,k).real()) > ep && std::abs(op(i,k).imag()) > ep) {op_l.coeffRef(i*hdim+j,k*hdim+j) += op(i,k);}
                if (std::abs(op(k,j).real()) > ep && std::abs(op(k,j).imag()) > ep) {op_l.coeffRef(i*hdim+j,i*hdim+k) -= op(k,j);}
            }
        }
    }
    return op_l;
}


SparseMatrix<double> QSim::liouville_lie_algebra_operator_double_sparse(const MatrixXd& op, double ep) {
    int hdim = op.cols();
    int ldim = hdim*hdim;
    SparseMatrix<double> op_l = SparseMatrix<double>(ldim, ldim);
    for(int i = 0; i < hdim; i++) {
        for (int j = 0; j < hdim; j++) {
            for(int k = 0; k < hdim; k++) {
                if (op(i,k) > ep) {op_l.coeffRef(i*hdim+j,k*hdim+j) += op(i,k);}
                if (op(k,j) > ep) {op_l.coeffRef(i*hdim+j,i*hdim+k) -= op(k,j);}
            }
        }
    }
    return op_l;
}


/*******************************************************************************/


MatrixXcd QSim::liouville_lie_group_operator_dcomplex_dense(const MatrixXcd& op, const MatrixXcd& op_inv, double ep) {
    const int hdim = op.cols();
    const int ldim = hdim*hdim ;
    MatrixXcd op_l = MatrixXcd::Zero(ldim, ldim);
    for (int i = 0; i < hdim; i++) {
        for (int j = 0; j < hdim; j++) {
            for (int k = 0; k < hdim; k++) {
                for (int l = 0; l < hdim; l++) {
                    const _cd val = op(i,k)*op_inv(l,j);
                    if (std::abs(val.real()) > ep && std::abs(val.imag()) > ep) { 
                        op_l(i*hdim+j, k*hdim+l) = val;
                    }
                }
            }
        }
    }
    return op_l;
}


MatrixXd QSim::liouville_lie_group_operator_double_dense(const MatrixXd& op, const MatrixXd& op_inv, double ep) {
    const int hdim = op.cols();
    const int ldim = hdim*hdim;
    MatrixXd op_l = MatrixXd::Zero(ldim, ldim);
    for (int i = 0; i < hdim; i++) {
        for (int j = 0; j < hdim; j++) {
            for (int k = 0; k < hdim; k++) {
                for (int l = 0; l < hdim; l++) {
                    const double val = op(i,k)*op_inv(l,j);
                    if (std::abs(val) > ep) { 
                        op_l(i*hdim+j, k*hdim+l) = val;
                    }
                }
            }
        }
    }
    return op_l;
}


SparseMatrix<_cd> QSim::liouville_lie_group_operator_dcomplex_sparse(const MatrixXcd& op, const MatrixXcd& op_inv, double ep) {
    const int hdim = op.cols();
    const int ldim = hdim*hdim ;
    SparseMatrix<_cd> op_l = SparseMatrix<_cd>(ldim, ldim);
    for (int i = 0; i < hdim; i++) {
        for (int j = 0; j < hdim; j++) {
            for (int k = 0; k < hdim; k++) {
                for (int l = 0; l < hdim; l++) {
                    const _cd val = op(i,k)*op_inv(l,j);
                    if (std::abs(val.real()) > ep && std::abs(val.imag()) > ep) { 
                        op_l.insert(i*hdim+j, k*hdim+l) = val;
                    }
                }
            }
        }
    }
    return op_l;
}


SparseMatrix<double> QSim::liouville_lie_group_operator_double_sparse(const MatrixXd& op, const MatrixXd& op_inv, double ep) {
    const int hdim = op.cols();
    const int ldim = hdim*hdim;
    SparseMatrix<double> op_l = SparseMatrix<double>(ldim, ldim);
    for (int i = 0; i < hdim; i++) {
        for (int j = 0; j < hdim; j++) {
            for (int k = 0; k < hdim; k++) {
                for (int l = 0; l < hdim; l++) {
                    const double val = op(i,k)*op_inv(l,j);
                    if (std::abs(val) > ep) { 
                        op_l.insert(i*hdim+j, k*hdim+l) = val;
                    }
                }
            }
        }
    }
    return op_l;
}


/*******************************************************************************/


SparseMatrix<_cd> QSim::coupled_tls_lindblad_dcomplex(int ntls, const _dvec& omegas, const std::vector<double>& dephasings, const std::vector<double>& relaxations,
        const std::vector<std::pair<int, int>>& connections, const _cdvec& couplings) {
    const int dim = 1 << ntls;
    SparseMatrix<_cd> result = SparseMatrix<_cd>(dim*dim, dim*dim);
    /* Diagonalelemente: Dephasings und Oszillation */
    {
        TupGen tup1 = TupGen(_ivec(ntls, 1));
        TupGen tup2 = TupGen(_ivec(ntls, 1));
        tup1.init_initial();
        for (int i = 0; i < dim; i++) {
            tup2.init_initial();
            for (int j = 0; j < dim; j++) {
                int index1 = 0;
                int index2 = 0;
                _cd component = {0.0, 0.0};
                for (int k = 0; k < ntls; k++) {
                    index1 += tup1[k] * (1 << k);
                    index2 += tup2[k] * (1 << k);
                    if (tup1[k] < tup2[k]) {
                        component += _cd{-dephasings[k], omegas[k]};
                    }
                    else if (tup1[k] > tup2[k]) {
                        component += _cd{-dephasings[k], -omegas[k]};
                    }
                    else if (tup1[k] == 1 && tup2[k] == 1) {
                        component += _cd{-relaxations[k], 0.0};
                    }
                }
                result.insert(index1*dim + index2, index1*dim + index2) = component;
                ++tup2;
            }
            ++tup1;
        }
    }
    /* Populationsfluss |1><1] -> |0><0| */
    {
        TupGen tup1 = TupGen(_ivec(ntls-1, 1));
        TupGen tup2 = TupGen(_ivec(ntls-1, 1));
        tup1.init_initial();
        for (int n = 0; n < ntls; n++) {
            const int shift = (1 << n)*(1 + dim);
            for (int i = 0; i < 1<<(ntls - 1); i++) {
                tup2.init_initial();
                for (int j = 0; j < 1<<(ntls - 1); j++) {
                    int index1 = 0;
                    int index2 = 0;
                    int m = 0;
                    for (int k = 0; k < ntls; k++) {
                        if (k == n) {
                            continue;
                        }
                        index1 += tup1[m] * (1 << k);
                        index2 += tup2[m] * (1 << k);
                        ++m;
                    }
                    result.insert(index1*dim + index2, index1*dim + index2 + shift) = relaxations[n];
                    ++tup2;
                }
                ++tup1;
            }
        }
    }
    /* Kopplungen für den Fall ntls > 2 */
    if (ntls > 2) {
        TupGen tup1 = TupGen(_ivec(ntls-2, 1));
        TupGen tup2 = TupGen(_ivec(ntls, 1));
        for (int i = 0; i < connections.size(); i++) {
            const int i1 = connections[i].first;
            const int i2 = connections[i].second;
            const _cd coupling = couplings[i];
            int shift1 = 1 << i1;
            int shift2 = 1 << i2;
            tup1.init_initial();
            for (int j1 = 0; j1 < (1 << (ntls - 2)); j1++) {
                tup2.init_initial();
                for (int j2 = 0; j2 < dim; j2++) {
                    int index1 = 0;
                    int index2 = 0;
                    int m = 0;
                    for (int k = 0; k < ntls; k++) {
                        if (k != i1 && k != i2) {
                            index1 += tup1[m] * (1 << k);
                            ++m;
                        }
                        index2 += tup2[k] * (1 << k);
                    }
                    // der Faktor -i ist aufgrund des Faktors (-i/hbar) auf der rechten Seite der Neumann-Gleichung für die Dichtematrix benötigt
                    result.insert((index1 + 0*shift1 + 1*shift2)*dim + index2, (index1 + 1*shift1 + 0*shift2)*dim + index2) = -Imd * coupling;
                    result.insert((index1 + 1*shift1 + 0*shift2)*dim + index2, (index1 + 0*shift1 + 1*shift2)*dim + index2) = -Imd * std::conj(coupling);
                    result.insert(index2*dim + index1 + 0*shift1 + 1*shift2, index2*dim + index1 + 1*shift1 + 0*shift2) = -Imd * (- std::conj(coupling));
                    result.insert(index2*dim + index1 + 1*shift1 + 0*shift2, index2*dim + index1 + 0*shift1 + 1*shift2) = -Imd * (- coupling);
                    ++tup2;
                }
                ++tup1;
            }
        }
    }
    /* Kopplungen für den Fall ntls == 2 */
    else if (ntls == 2) {
        TupGen tup = TupGen(_ivec(ntls, 1));
        for (int i = 0; i < connections.size(); i++) {
            const int i1 = connections[i].first;
            const int i2 = connections[i].second;
            const _cd coupling = couplings[i];
            int shift1 = 1 << i1;
            int shift2 = 1 << i2;
            tup.init_initial();
            for (int j = 0; j < dim; j++) {
                int index = 0;
                for (int k = 0; k < ntls; k++) {
                    index += tup[k]*(1 << k);
                }
                result.insert((0*shift1 + 1*shift2)*dim + index, (1*shift1 + 0*shift2)*dim + index) = -Imd * coupling;
                result.insert((1*shift1 + 0*shift2)*dim + index, (0*shift1 + 1*shift2)*dim + index) = -Imd * std::conj(coupling);
                result.insert(index*dim + 0*shift1 + 1*shift2, index*dim + 1*shift1 + 0*shift2) = -Imd * (- std::conj(coupling));
                result.insert(index*dim + 1*shift1 + 0*shift2, index*dim + 0*shift1 + 1*shift2) = -Imd * (- coupling);
                ++tup;
            }
        }
    }
    return result;
}


SparseMatrix<_cd> QSim::coupled_tls_lindblad_double(int ntls, const _dvec& omegas, const std::vector<double>& dephasings, const std::vector<double>& relaxations,
        const std::vector<std::pair<int, int>>& connections, const _dvec& couplings) {
    const int dim = 1 << ntls;
    SparseMatrix<_cd> result = SparseMatrix<_cd>(dim*dim, dim*dim);
    /* Diagonalelemente: Dephasings und Oszillation */
    {
        TupGen tup1 = TupGen(_ivec(ntls, 1));
        TupGen tup2 = TupGen(_ivec(ntls, 1));
        tup1.init_initial();
        for (int i = 0; i < dim; i++) {
            tup2.init_initial();
            for (int j = 0; j < dim; j++) {
                int index1 = 0;
                int index2 = 0;
                _cd component = {0.0, 0.0};
                for (int k = 0; k < ntls; k++) {
                    index1 += tup1[k] * (1 << k);
                    index2 += tup2[k] * (1 << k);
                    if (tup1[k] < tup2[k]) {
                        component += _cd{-dephasings[k], omegas[k]};
                    }
                    else if (tup1[k] > tup2[k]) {
                        component += _cd{-dephasings[k], -omegas[k]};
                    }
                    else if (tup1[k] == 1 && tup2[k] == 1) {
                        component += _cd{-relaxations[k], 0.0};
                    }
                }
                result.insert(index1*dim + index2, index1*dim + index2) = component;
                ++tup2;
            }
            ++tup1;
        }
    }
    /* Populationsfluss |1><1] -> |0><0| */
    {
        TupGen tup1 = TupGen(_ivec(ntls-1, 1));
        TupGen tup2 = TupGen(_ivec(ntls-1, 1));
        tup1.init_initial();
        for (int n = 0; n < ntls; n++) {
            const int shift = (1 << n)*(1 + dim);
            for (int i = 0; i < 1<<(ntls - 1); i++) {
                tup2.init_initial();
                for (int j = 0; j < 1<<(ntls - 1); j++) {
                    int index1 = 0;
                    int index2 = 0;
                    int m = 0;
                    for (int k = 0; k < ntls; k++) {
                        if (k == n) {
                            continue;
                        }
                        index1 += tup1[m] * (1 << k);
                        index2 += tup2[m] * (1 << k);
                        ++m;
                    }
                    result.insert(index1*dim + index2, index1*dim + index2 + shift) = relaxations[n];
                    ++tup2;
                }
                ++tup1;
            }
        }
    }
    /* Kopplungen für den Fall ntls > 2 */
    if (ntls > 2) {
        TupGen tup1 = TupGen(_ivec(ntls-2, 1));
        TupGen tup2 = TupGen(_ivec(ntls, 1));
        for (int i = 0; i < connections.size(); i++) {
            int i1 = connections[i].first;
            int i2 = connections[i].second;
            const double coupling = couplings[i];
            int shift1 = 1 << i1;
            int shift2 = 1 << i2;
            tup1.init_initial();
            for (int j1 = 0; j1 < (1 << (ntls - 2)); j1++) {
                tup2.init_initial();
                for (int j2 = 0; j2 < dim; j2++) {
                    int index1 = 0;
                    int index2 = 0;
                    int m = 0;
                    for (int k = 0; k < ntls; k++) {
                        if (k != i1 && k != i2) {
                            index1 += tup1[m] * (1 << k);
                            ++m;
                        }
                        index2 += tup2[k] * (1 << k);
                    }
                    result.insert((index1 + 0*shift1 + 1*shift2)*dim + index2, (index1 + 1*shift1 + 0*shift2)*dim + index2) = coupling;
                    result.insert((index1 + 1*shift1 + 0*shift2)*dim + index2, (index1 + 0*shift1 + 1*shift2)*dim + index2) = coupling;
                    result.insert(index2*dim + index1 + 0*shift1 + 1*shift2, index2*dim + index1 + 1*shift1 + 0*shift2) = - coupling;
                    result.insert(index2*dim + index1 + 1*shift1 + 0*shift2, index2*dim + index1 + 0*shift1 + 1*shift2) = - coupling;
                    ++tup2;
                }
                ++tup1;
            }
        }
    }
    /* Kopplungen für den Fall ntls == 2 */
    else if (ntls == 2) {
        TupGen tup = TupGen(_ivec(ntls, 1));
        for (int i = 0; i < connections.size(); i++) {
            const int i1 = connections[i].first;
            const int i2 = connections[i].second;
            const double coupling = couplings[i];
            int shift1 = 1 << i1;
            int shift2 = 1 << i2;
            tup.init_initial();
            for (int j = 0; j < dim; j++) {
                int index = 0;
                for (int k = 0; k < ntls; k++) {
                    index += tup[k]*(1 << k);
                }
                result.insert((0*shift1 + 1*shift2)*dim + index, (1*shift1 + 0*shift2)*dim + index) = coupling;
                result.insert((1*shift1 + 0*shift2)*dim + index, (0*shift1 + 1*shift2)*dim + index) = coupling;
                result.insert(index*dim + 0*shift1 + 1*shift2, index*dim + 1*shift1 + 0*shift2) = - coupling;
                result.insert(index*dim + 1*shift1 + 0*shift2, index*dim + 0*shift1 + 1*shift2) = - coupling;
                ++tup;
            }
        }
    }
    return result;
}


/*******************************************************************************/


MatrixXd QSim::coupled_tls_hamiltonian_double(int ntls, const _dvec& omegas, const std::vector<std::pair<int, int>>& connections, const _dvec& couplings) {
    const int dim = 1 << ntls;
    MatrixXd hamiltonian = MatrixXd::Zero(dim, dim);
    // keine Selbstkopplungen im Hamiltonian
    for (int i = 0; i < connections.size(); i++) {
        if (ntls > 2) {
            int ta = connections[i].first;;
            int tb = connections[i].second;
            const int shift_a = 1 << ta;
            const int shift_b = 1 << tb;
            TupGen gen = TupGen(_ivec(ntls-2, 1));
            for (gen.init_initial(); !gen.is_final(); ++gen) {
                int index = 0;
                int m = 0;
                for (int k = 0; k < ntls; k++) {
                    if (k == ta || k == tb) {
                        continue;
                    }
                    index += gen[m] * (1<<k);
                    ++m;
                }
                hamiltonian(index + 0*shift_a + 1*shift_b, index + 1*shift_a + 0*shift_b) = couplings[i];
                hamiltonian(index + 1*shift_a + 0*shift_b, index + 0*shift_a + 1*shift_b) = couplings[i];
            }
            int index = 0;
            int m = 0;
            for (int k = 0; k < ntls; k++) {
                if (k == ta || k == tb) {
                    continue;
                }
                index += gen[m] * (1<<k);
                ++m;
            }
            hamiltonian(index + 0*shift_a + 1*shift_b, index + 1*shift_a + 0*shift_b) = couplings[i];
            hamiltonian(index + 1*shift_a + 0*shift_b, index + 0*shift_a + 1*shift_b) = couplings[i];
        }
        else {
            const int shift_a = 1 << connections[i].first;
            const int shift_b = 1 << connections[i].second;
            hamiltonian(0*shift_a + 1*shift_b, 1*shift_a + 0*shift_b) = couplings[i];
            hamiltonian(1*shift_a + 0*shift_b, 0*shift_a + 1*shift_b) = couplings[i];
        }
    }
    TupGen gen = TupGen(_ivec(ntls, 1));
    for (gen.init_initial(); !gen.is_final(); ++gen) {
        int index = 0;
        double omega = 0.0;
        for (int i = 0; i < ntls; i++) {
            index += gen[i] * (1 << i);
            omega += static_cast<double>(gen[i]) * omegas[i];
        }
        hamiltonian(index, index) = omega;
    }
    {
        int index = 0;
        double omega = 0.0;
        for (int i = 0; i < ntls; i++) {
            index += gen[i] * (1 << i);
            omega += static_cast<double>(gen[i]) * omegas[i];
        }
        hamiltonian(index, index) = omega;
    }
    return hamiltonian;
}


/*******************************************************************************/


MatrixXd QSim::coupled_tls_operator_double(int ntls, const _dvec& entries) {
    const int dim = 1 << ntls;
    MatrixXd tensor_operator = MatrixXd::Zero(dim, dim);
    for (int i = 0; i < ntls; i++) {
        if (ntls > 1) {
            TupGen gen = TupGen(_ivec(ntls-1, 1));
            const int shift = 1 << i;
            for (gen.init_initial(); !gen.is_final(); ++gen) {
                int index = 0;
                int m = 0;
                for (int k = 0; k < ntls; k++) {
                    if (k == i) {
                        continue;
                    }
                    index += gen[m] * (1<<k);
                    ++m;
                }
                tensor_operator(index + 0*shift, index + 1*shift) = entries[i];
                tensor_operator(index + 1*shift, index + 0*shift) = entries[i];
            }
            int index = 0;
            int m = 0;
            for (int k = 0; k < ntls; k++) {
                if (k == i) {
                    continue;
                }
                index += gen[m] * (1<<k);
                ++m;
            }
            tensor_operator(index + 0*shift, index + 1*shift) = entries[i];
            tensor_operator(index + 1*shift, index + 0*shift) = entries[i];
        }
        else {
            tensor_operator(0, 1) = entries[i];
            tensor_operator(1, 0) = entries[i];
        }
    }
    return tensor_operator;
}


/*******************************************************************************/


SparseMatrix<_cd> QSim::eigen_sparse_from_dense_dcomplex(const MatrixXcd& m, double epsilon) {
    const int rows = m.rows();
    const int cols = m.cols();
    SparseMatrix<_cd> result = SparseMatrix<_cd>(rows, cols);
    for (int j = 0; j < cols; j++) {
        for (int i = 0; i < rows; i++) {
            if (std::abs(m(i, j).real()) > epsilon || std::abs(m(i, j).imag()) > epsilon) {
                result.insert(i, j) = m(i, j);
            }
        }
    }
    return result;
}


SparseMatrix<double> QSim::eigen_sparse_from_dense_double(const MatrixXd& m, double epsilon) {
    const int rows = m.rows();
    const int cols = m.cols();
    SparseMatrix<double> result = SparseMatrix<double>(rows, cols);
    for (int j = 0; j < cols; j++) {
        for (int i = 0; i < rows; i++) {
            if (std::abs(m(i, j)) > epsilon) {
                result.insert(i, j) = m(i, j);
            }
        }
    }
    return result;
}


/*******************************************************************************/


void QSim::eigen_prune_dcomplex(MatrixXcd& m, double epsilon) {
    for (int j = 0; j < m.cols(); j++) {
        for (int i = 0; i < m.rows(); i++) {
            if (std::abs(m(i, j).real()) <= epsilon && std::abs(m(i, j).imag()) <= epsilon) {
                m(i, j) = {0.0, 0.0};
            }
        }
    }
}


void QSim::eigen_prune_double(MatrixXd& m, double epsilon) {
    for (int j = 0; j < m.cols(); j++) {
        for (int i = 0; i < m.rows(); i++) {
            if (std::abs(m(i, j)) <= epsilon) {
                m(i, j) = 0.0;
            }
        }
    }
}


/*******************************************************************************/


Matrix2cd QSim::random_diagonal_state_2() {
    std::default_random_engine rand_gen = std::default_random_engine();
    std::uniform_real_distribution<double> rand_dist = std::uniform_real_distribution<double>(0.0, 1.0);
    auto dice = std::bind(rand_dist, rand_gen);
    Matrix2cd random_state = Matrix2cd();
    random_state(0,0) = dice();
    random_state(1,1) = dice();
    random_state(0,1) = 0.0;
    random_state(1,0) = 0.0;
    random_state = (random_state.array() / random_state.trace()).matrix();
    return random_state;
}

Matrix2cd QSim::random_pure_state_2() {
    std::default_random_engine rand_gen = std::default_random_engine();
    std::uniform_real_distribution<double> rand_dist = std::uniform_real_distribution<double>(0.0, 1.0);
    auto dice = std::bind(rand_dist, rand_gen);
    Vector2cd state_vector = Vector2cd();
    state_vector(0) = _cd(dice(), dice());
    state_vector(1) = _cd(dice(), dice());
    _cd norm = std::sqrt( (state_vector.conjugate().transpose() * state_vector)(0,0) );
    state_vector(0) /= norm;
    state_vector(1) /= norm;
    Matrix2cd random_state = state_vector * state_vector.conjugate().transpose();
    return random_state;
}

MatrixXcd QSim::random_diagonal_state_n(int n) {
    std::default_random_engine rand_gen = std::default_random_engine();
    std::uniform_real_distribution<double> rand_dist = std::uniform_real_distribution<double>(0.0, 1.0);
    auto dice = std::bind(rand_dist, rand_gen);
    MatrixXcd random_state = MatrixXcd(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {random_state(i,i) = dice();}
            else {random_state(i,j) = 0.0;}
        }
    }
    random_state = (random_state.array() / random_state.trace()).matrix();
    return random_state;
}

MatrixXcd QSim::random_pure_state_n(int n) {
    std::default_random_engine rand_gen = std::default_random_engine();
    std::uniform_real_distribution<double> rand_dist = std::uniform_real_distribution<double>(0.0, 1.0);
    auto dice = std::bind(rand_dist, rand_gen);
    VectorXcd state_vector = VectorXcd(n);
    for (int i = 0; i < n; i++) {
        state_vector(i) = _cd(dice(), dice());
    }
    _cd norm = std::sqrt( (state_vector.conjugate().transpose() * state_vector)(0,0) );
    for (int i = 0; i < n; i++) {
        state_vector(i) /= norm;
    }
    MatrixXcd random_state = state_vector * state_vector.conjugate().transpose();
    for (int i = 0; i < n; i++) {
        random_state(i,i) = random_state(i,i).real();
    }
    return random_state;
}


