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
    _dvec tls_coupling_strengths = {};
    for (int i = 0; i < system_init["connections"].size(); i++) {
        tls_connections.push_back({system_init["connections"][i]["i"].get<int>(), system_init["connections"][i]["j"].get<int>()});
        tls_coupling_strengths.push_back(system_init["connections"][i]["strength"].get<double>());
    }
    const _dvec dipole_moments = system_init["dipole moments"].get<_dvec>();
    const double eta = system_init["eta"].get<double>();
    const _dvec relaxation_rates = system_init["relaxation rates"].get<_dvec>();

    const _dvec pure_dephasing_rates = system_init["pure dephasing rates"].get<_dvec>();
    _dvec dephasing_rates = {};
    for (int i = 0; i < ntls; i++) {
        dephasing_rates.push_back(0.5*relaxation_rates[i] + pure_dephasing_rates[i]);
    }

    const int dim = 1 << ntls;

    MatrixXd hamiltonian = coupled_tls_hamiltonian_double(ntls, tls_omegas, tls_connections, tls_coupling_strengths);

    SelfAdjointEigenSolver<MatrixXd> solver = SelfAdjointEigenSolver<MatrixXd>();
    solver.compute(hamiltonian);
    MatrixXd hamiltonian_new = solver.eigenvalues().asDiagonal();
    MatrixXd transform = solver.eigenvectors();
    MatrixXd transform_inv = transform.conjugate().transpose();
    VectorXd eigenvalues = solver.eigenvalues();

    _dvec new_omegas = _dvec(dim, 0.0);
    for(int i = 0; i < dim; i++) {
        new_omegas[i] = eigenvalues(i);
    }
    MatrixXd dephasing_mat = MatrixXd::Zero(dim, dim);
    for(int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            dephasing_mat(i, j) = dephasing_rates[0];
        }
    }
    Propagator propagator = Propagator::init_naive(dephasing_mat, new_omegas, eta);


    const MatrixXd dipole_op = coupled_tls_operator_double(ntls, dipole_moments);
    const MatrixXcd dipole_op_new = transform_inv * dipole_op * transform;
    VectorXcd observable_new = dipole_op_new.reshaped<RowMajor>().eval();
    SparseVector<_cd> observable_new_sparse = SparseVector<_cd>(dim*dim);
    for (int i = 0; i < dim*dim; i++) {
        if (observable_new[i] != 0.0) {
            observable_new_sparse.insert(i) = observable_new[i];
        }
    }
    SparseMatrix<_cd> couplings = liouville_from_hilbert_couplings_dcomplex(dipole_op_new, dim);
    
    std::cout << "\n\nHamiltonian in der alten Basis:\n" << std::endl;
    std::cout << hamiltonian << std::endl;
    std::cout << "\n\nHamiltonian in der neuen Basis:\n" << std::endl;
    std::cout << hamiltonian_new << std::endl;
    std::cout << "\n\nBasiswechselabbildung:\n" << std::endl;
    std::cout << transform << std::endl;
    std::cout << "\n\nDipoloperator in der alten Basis:\n" << std::endl;
    std::cout << dipole_op << std::endl;
    std::cout << "\n\nDipoloperator in der neuen Basis:\n" << std::endl;
    std::cout << dipole_op_new << std::endl;


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
            couplings,
            envelopes,
            i_init["center frequencies"].get<_dvec>(),
            i_init["phase matches"].get<int>(),
            i_init["phase match signatures"].get<_ivec>(),
        };

        /*** Initialisierung der Berechnungsparameter ***/

        ComputationParams computation_params = init_json["computation params"].get<ComputationParams>();
        CubaCommonParams cuba_common_params = init_json["cuba params"].get<CubaCommonParams>();

        /*** Initialisierung der Parameter der Antwortfunktion ***/
        
        std::vector<Domain> arg_domains = init_json["domains"].get<std::vector<Domain>>();
        ResponseParams response_params = ResponseParams {
            observable_new_sparse,
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
