#include <complex>
#include <vector>
#include <iostream>
#include <fstream>
#include <functional>
#include <memory>
#include <random>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <json/json.hpp>

#include "src/definitions.hpp"
#include "src/algebra.hpp"
#include "src/functions.hpp"
#include "src/propagator.hpp"
#include "src/response_kernels.hpp"
#include "src/response.hpp"
#include "src/integration.hpp"

using namespace Eigen;
using json = nlohmann::json;

using namespace QSim;



int main(void) {

    // Multipliziere alle Zeiten mit und dividiere alle Frequenzen durch diese Zahl
    const double normalize = 1.0e+15;
    const double f_to_E_conversion = hbar*normalize/e;

    const int order = 3;

    const int nlevels = 4;
    const int dim = nlevels*nlevels;

    const double omega_a = twopi * 365.0e12 / normalize;
    const double omega_b = twopi * 397.0e12 / normalize;
    const double mu_a = -1.1;     // In Debye
    const double mu_b = 1.5;      // In Debye
    const double tls_coupling = hbar / pi * 50.0e12 / normalize;
    const double xi = 0.066 / f_to_E_conversion;
    const double dephasing = 0.02e+15 / normalize;
    const double exciton_decay = 0.1 * dephasing;
    const double eta = 0.01 * dephasing;

    const double w_d = 0.5*(omega_a - omega_b);
    const double w_s = 0.5*(omega_a + omega_b);
    const double rho = std::sqrt(w_d*w_d + xi*xi);

    const double w1 = 0.0;
    const double w2 = w_s - rho;
    const double w3 = w_s + rho;
    const double w4 = omega_a + omega_b;

    _dvec w_q = {w1, w2, w3, w4};

    const double sqrt_p_term = std::sqrt(2.0 * (w_d*w_d + xi*xi + w_d*rho));
    const double sqrt_m_term = std::sqrt(2.0 * (w_d*w_d + xi*xi - w_d*rho));

    const double f_xi_p = xi / sqrt_p_term;
    const double f_xi_m = xi / sqrt_m_term;
    const double f_p = (w_d + rho) / sqrt_p_term;
    const double f_m = (w_d - rho) / sqrt_m_term;

    MatrixXd dipole_moments = MatrixXd::Zero(4, 4);
    dipole_moments(0,1) =   f_xi_p*mu_a + f_xi_m*mu_b;
    dipole_moments(0,2) = - f_p*mu_a    - f_m*mu_b;
    dipole_moments(3,1) =   f_xi_m*mu_a + f_xi_p*mu_b;
    dipole_moments(3,2) = - f_m*mu_a    - f_p*mu_b;
    dipole_moments(1,0) = dipole_moments(0,1);
    dipole_moments(2,0) = dipole_moments(0,2);
    dipole_moments(1,3) = dipole_moments(3,1);
    dipole_moments(2,3) = dipole_moments(3,2);
    SparseMatrix<double> couplings = liouville_from_hilbert_couplings_double(dipole_moments, nlevels);


    _dvec decay_rates = {0.0, exciton_decay, exciton_decay, 2*exciton_decay};
    MatrixXd dephasings = MatrixXd::Zero(nlevels, nlevels);
    for (int i = 0; i < nlevels; i++) {
        for (int j = 0; j < nlevels; j++) {
            if (i != j) {
                dephasings(i,j) = dephasing;
            }
        }
    }
    Propagator propagator = Propagator::init_naive(dephasings, w_q, eta);



    _dvec w_pulse = _dvec(3, w_s);
    double ampl_pulse = 1.0;
    double sigma_pulse = 1.0e-15 * normalize;
    std::vector<Vector3d> k0 = { {1.0, -1.0, 0.0},  {1.0, -1.0, 0.0},  {1.0, 1.0, 0.0} };
    Vector3d ksig = {1.0, 1.0, 0.0};
    VectorXd k_sig = ksig;
    std::vector<VectorXd> k = {};
    for (int i = 0; i < order; i++) {
        k.push_back(k0[i]);
    }

    std::vector<Env> vertex_fns = std::vector<Env>();
    for (int i = 0; i < 3; i++) {
        vertex_fns.emplace_back(new EnvGaussian{ampl_pulse, sigma_pulse});
    }

    Interaction interaction = Interaction {
        order,
        couplings,
        vertex_fns,
        w_pulse,
        k,
        k_sig,
        0.001, // max_k_mismatch
    };




    VectorXd observable_dense = dipole_moments.reshaped<RowMajor>().eval();
    SparseVector<double> observable = SparseVector<double>(dim);
    for (int i = 0; i < dim; i++) {
        if (observable_dense[i] != 0.0) {
            observable.insert(i) = observable_dense[i];
        }
    }



    ResponseParams response_params = ResponseParams {
        observable,
        _dvec(order, -5.0),
        _dvec(order, 5.0),
        propagator,
        interaction,
        {Domain::Freq, Domain::Time, Domain::Freq},
    };




    std::cout << "Energien in der alten Basis:" << std::endl;
    std::cout << "omega_a = " << f_to_E_conversion*omega_a << " eV,   " << omega_a << "*10^15 Hz" << std::endl;
    std::cout << "omega_b = " << f_to_E_conversion*omega_b << " eV,   " << omega_b << "*10^15 Hz" << std::endl;
    std::cout << "Energien in der neuen Basis:" << std::endl;
    for (int i = 0; i < 4; i++) {
        std::cout << f_to_E_conversion*w_q[i] << " eV,   " << w_q[i] << "*10^15 Hz" << std::endl;
    }
    std::cout << "Kopplungsenergie und -frequenz:" << std::endl;
    std::cout << "E_J = " << f_to_E_conversion*xi << " eV,   " << xi << "*10^15 Hz" << std::endl;
    std::cout << "Kopplungen in der alten Basis:" << std::endl;
    std::cout << "mu_a = " << mu_a << std::endl;
    std::cout << "mu_b = " << mu_b << std::endl;
    std::cout << "Kopplungen in der neuen Basis:" << std::endl;
    std::cout << dipole_moments << "\n" << std::endl;





    const double energy_start = 1.300;       // in eV
    const double energy_stop  = 1.800;       // in eV
    const double start  = e*hbarinv * energy_start / normalize; // in 2pi*Hz normalisiert
    const double stop   = e*hbarinv * energy_stop / normalize;  // in 2pi*Hz normalisiert

    int nthreads = 12;
    _dvec lower_borders = _dvec{start, 0.0, start};
    _dvec upper_borders = _dvec{stop, 0.0, stop};
    _ivec nintervals = _ivec{99, 0, 99};
    CubaAlg cuba_alg = CubaAlg::Divonne;
    Representation repr = Representation::DysonFreq;

    ComputationParams computation_params = ComputationParams {
        nthreads, lower_borders, upper_borders, nintervals, cuba_alg, repr
    };


    double epsrel = 1.0e-2;
    double epsabs = 0.0;
    int flags = 0;
    int mineval = 10000;
    int maxeval = 100000;

    CubaCommonParams cuba_common_params = CubaCommonParams {
        epsrel, epsabs, flags, mineval, maxeval
    };


    std::string out_dir = std::string("../sim_instances/template");
    std::string filename_init = out_dir + "/init.json";
    json init_json = json {
        {"response_params", response_params},
        {"computation_params", computation_params},
        {"cuba_common_params", cuba_common_params},
    };

    std::ofstream file;
    std::string init_string = init_json.dump(1, '\t');
    file = std::ofstream(filename_init);
    file.write(&init_string[0], init_string.size());
    file.close();


}
