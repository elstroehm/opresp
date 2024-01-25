#include <complex>
#include <vector>
#include <iostream>
#include <fstream>
#include <functional>
#include <memory>
#include <random>
#include <utility>

#include "external/eigen-3.4.0/Eigen/Core"
#include "external/eigen-3.4.0/Eigen/Sparse"
#include "external/eigen-3.4.0/Eigen/Eigenvalues"
#include "external/json/json.hpp"
#include "external/hdf5-1.14.1-2/hdf5/include/H5Cpp.h"

#include "src/definitions.hpp"
#include "src/algebra.hpp"
#include "src/integration.hpp"
#include "src/functions.hpp"
#include "src/propagator.hpp"
#include "src/response_kernels.hpp"
#include "src/response.hpp"

using namespace Eigen;
using json = nlohmann::json;

using namespace QSim;



int main(int argc, char** argv) {



    std::ifstream ifile = std::ifstream(argv[1]);
    json init_json;
    ifile >> init_json;
    ifile.close();



    /*** Systeminitialisierung ***/

    json system_init = init_json["system"];

    const int ntls = system_init["ntls"].get<int>();
    const _dvec tls_omegas = system_init["omegas"].get<_dvec>();
    std::vector<std::pair<int, int>> tls_connections = {};
    _cdvec tls_coupling_strengths = {};
    for (int i = 0; i < system_init["connections"].size(); i++) {
        tls_connections.push_back({system_init["connections"][i]["i"].get<int>(), system_init["connections"][i]["j"].get<int>()});
        tls_coupling_strengths.push_back(system_init["connections"][i]["strength"].get<double>());
    }
    const _dvec dipole_moments = system_init["dipole moments"].get<_dvec>();
    const double eta = system_init["eta"].get<double>();
    const _dvec relaxation_rates = system_init["relaxation rates"].get<_dvec>();

    // const _dvec dephasing_rates = system_init["dephasing rates"].get<_dvec>();

    const _dvec pure_dephasing_rates = system_init["pure dephasing rates"].get<_dvec>();
    _dvec dephasing_rates = {};
    for (int i = 0; i < ntls; i++) {
        dephasing_rates.push_back(0.5*relaxation_rates[i] + pure_dephasing_rates[i]);
    }

    const int dim = 1 << ntls;

    SparseMatrix<_cd> lindblad = coupled_tls_lindblad_dcomplex(ntls, tls_omegas, dephasing_rates, relaxation_rates, tls_connections, tls_coupling_strengths);
    MatrixXcd lindblad_dense = lindblad;
    ComplexEigenSolver<MatrixXcd> solver = ComplexEigenSolver<MatrixXcd>();
    solver.compute(lindblad_dense, true);
    VectorXcd eigenvalues = solver.eigenvalues();
    MatrixXcd new_lindblad = eigenvalues.asDiagonal();

    MatrixXcd transform = solver.eigenvectors();
    MatrixXcd transform_inv = transform.inverse();
    SparseMatrix<_cd> transform_sparse = eigen_sparse_from_dense_dcomplex(transform, 0.0);
    SparseMatrix<_cd> transform_inv_sparse = eigen_sparse_from_dense_dcomplex(transform_inv, 0.0);
    Propagator propagator = Propagator::init_lindblad_diagonal(eigenvalues, eta);

    const MatrixXd dipole_op = coupled_tls_operator_double(ntls, dipole_moments);
    VectorXcd observable = dipole_op.reshaped<RowMajor>().eval();
    VectorXcd observable_new = transform_inv * observable;
    SparseVector<_cd> observable_new_sparse = SparseVector<_cd>(dim*dim);
    for (int i = 0; i < dim*dim; i++) {
        if (observable_new[i] != 0.0) {
            observable_new_sparse.insert(i) = observable_new[i];
        }
    }
    SparseMatrix<double> couplings = liouville_lie_algebra_operator_double_sparse(dipole_op, 0.0);
    SparseMatrix<_cd> couplings_new = (transform_inv_sparse * couplings * transform_sparse).eval();
    MatrixXcd couplings_new_dense = couplings;

    SparseMatrix<_cd> qform = liouville_qform_dcomplex_sparse(dim);
    SparseMatrix<_cd> qform_new = (transform_sparse.transpose() * qform * transform_sparse).eval();
    // SparseMatrix<_cd> qform_new = qform;



    // std::cout << "\n\nLindblad Operator:\n" << std::endl;
    // std::cout << lindblad_dense << std::endl;
    // std::cout << "\n\nBasiswechselabbildung im Liouville-Bild:\n" << std::endl;
    // std::cout << transform << std::endl;
    // std::cout << "\n\nInverse Basiswechselabbildung im Liouville-Bild:\n" << std::endl;
    // std::cout << transform << std::endl;
    // std::cout << "\n\nDiagonalisierter Lindblad Operator:\n" << std::endl;
    // std::cout << new_lindblad << std::endl;
    // std::cout << "\n\nDipoloperator in der alten Basis:\n" << std::endl;
    // std::cout << couplings << std::endl;
    // std::cout << "\n\nDipoloperator in der neuen Basis:\n" << std::endl;
    // std::cout << couplings_new_dense << std::endl;
    // std::cout << "\n\nAlter Dipoloperator im Hilbertraum:\n" << std::endl;
    // std::cout << dipole_op << std::endl;

    // const double teval = 0.2;
    // MatrixXcd rho = random_pure_state_n(4);
    // std::cout << "\n\n\n\n" << std::endl;
    // std::cout << rho << std::endl;
    // std::cout << "\n\n" << std::endl;
    // _cdvec evolve = _cdvec(4, 0.0);
    // MatrixXcd hamiltonian = MatrixXcd::Zero(4, 4);
    // hamiltonian(0,0) = 0.0;
    // hamiltonian(1,1) = tls_omegas[0];
    // hamiltonian(2,2) = tls_omegas[1];
    // hamiltonian(3,3) = tls_omegas[0] + tls_omegas[1];
    // hamiltonian(1,2) = tls_coupling_strengths[0];
    // hamiltonian(2,1) = tls_coupling_strengths[0];
    // SelfAdjointEigenSolver<MatrixXcd> solver2 = SelfAdjointEigenSolver<MatrixXcd>();
    // solver2.compute(hamiltonian);
    // MatrixXcd htransform = solver2.eigenvectors();
    // MatrixXcd hdiag = solver2.eigenvalues().asDiagonal();
    // MatrixXcd Uevolve = MatrixXcd::Zero(4, 4);
    // for (int i = 0; i < 4; i++) {
    //     Uevolve(i, i) = std::exp(-Imd * solver2.eigenvalues()(i) * teval);
    // }
    // MatrixXcd dipo = MatrixXcd::Zero(4, 4);
    // dipo(1,2) = tls_coupling_strengths[0];
    // dipo(2,1) = tls_coupling_strengths[0];
    // MatrixXcd hdipo = htransform * dipo * htransform.conjugate().transpose();
    // MatrixXcd hrho1 = htransform * Uevolve * htransform.conjugate().transpose() * rho * htransform * Uevolve.conjugate().transpose() * htransform.conjugate().transpose();
    // MatrixXcd lrho = htransform * Uevolve * htransform.conjugate().transpose() * rho * htransform * Uevolve.conjugate().transpose() * htransform.conjugate().transpose();
    // SparseMatrix<_cd> lindi_sparse = tls_lindblad_dcomplex(2, tls_omegas, _dvec(2, 0.0), _dvec(2, 0.0), tls_connections, tls_coupling_strengths);
    // MatrixXcd lindi = lindi_sparse;
    // solver.compute(lindi);
    // MatrixXcd ltransform = solver.eigenvectors();
    // MatrixXcd ltransform_inv = ltransform.inverse();
    // MatrixXcd ldiag = solver.eigenvalues().asDiagonal();
    // VectorXcd ldipo = dipo.reshaped<RowMajor>();
    // Propagator Levolve = Propagator::init_lindblad_diagonal(solver.eigenvalues(), eta);
    // MatrixXcd hrho2 = (ltransform * Levolve.eval_t(teval) * ltransform_inv * rho.reshaped<RowMajor>()).reshaped<RowMajor>(4, 4);

    // std::cout << "\n\n" << std::endl;
    // std::cout << hrho1 << std::endl;
    // std::cout << "\n\n" << std::endl;
    // std::cout << hrho2 << std::endl;
    // std::cout << "\n\n" << std::endl;
    // std::cout << hamiltonian << std::endl;
    // std::cout << "\n\n" << std::endl;
    // std::cout << hdiag << std::endl;
    // std::cout << "\n\n" << std::endl;
    // std::cout << lindi << std::endl;
    // std::cout << "\n\n" << std::endl;
    // std::cout << ldiag << std::endl;
    // std::cout << "\n\n" << std::endl;
    // std::cout << (hdipo * hrho1).trace() << "   " << ldipo.dot(hrho2.reshaped<RowMajor>()) << std::endl;
    // std::cout << "\n\n\n\n" << std::endl;




    /* Wenn eine Datei für das Spektrum angegeben wurde, berechne Antwortfunktion, ansonsten nur Konsolenausgabe des diagonalisierten Systems */
    if (argc == 3) {

        /*** Initialisierung der Wechselwirkung ***/

        json i_init = init_json["interaction"];
        std::vector<Env> envelopes = {};
        for (int i = 0; i < i_init["envelopes"].size(); i++) {
            envelopes.push_back( Env(EnvBase::from_json(i_init["envelopes"][i])) );
        }
        Interaction interaction = Interaction {
            i_init["ncontacts"].get<int>(),
            couplings_new,
            envelopes,
            i_init["w0"].get<_dvec>(),
            i_init["phase_matches"].get<int>(),
            i_init["phase_match_signatures"].get<_ivec>(),
        };

        /*** Initialisierung der Berechnungsparameter ***/

        ComputationParams computation_params = init_json["computation params"].get<ComputationParams>();
        CubaCommonParams cuba_common_params = init_json["cuba params"].get<CubaCommonParams>();

        /*** Initialisierung der Parameter der Antwortfunktion ***/
        
        std::vector<Domain> arg_domains = init_json["domains"].get<std::vector<Domain>>();
        ResponseParams response_params = ResponseParams {
            observable_new_sparse,
            qform_new,
            _dvec(interaction.ncontacts, -5.0),
            _dvec(interaction.ncontacts, 5.0),
            propagator,
            interaction,
            arg_domains,
        };

        /*** Initialisierung und Berechnung der Antwortfunktion ***/

        CubaIntegrator integrator = init_kernel(response_params, computation_params, cuba_common_params);
        int ncores = 0;
        int pcores = 0;
        cubacores(&ncores, &pcores);
        _dvec spectral_data = spectrum(integrator, response_params, computation_params);

        /*** Speichere Werte im hdf5 Dateiformat ***/

        const int order = response_params.interaction.ncontacts;
        int rank = order + 1;
        hsize_t* dims = new hsize_t[rank];
        for (int i = 0; i < order; i++) {
            dims[i] = computation_params.nintervals[i] + 1;
        }
        dims[order] = order + 4;

        H5::DataSpace dataspace = H5::DataSpace(rank, dims, NULL);
        H5::H5File outfile = H5::H5File(argv[2], H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        H5::DataSet dataset = outfile.createDataSet("spectrum", H5::PredType::NATIVE_DOUBLE, dataspace);

        dataset.write(&spectral_data[0], H5::PredType::NATIVE_DOUBLE);
        dataset.close();
        dataspace.close();
        outfile.close();

        delete[] dims;

    }

    return 0;

}
