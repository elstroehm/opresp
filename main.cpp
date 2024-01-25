#include <complex>
#include <vector>
#include <iostream>
#include <functional>
#include <memory>
#include <random>
#include <thread>
#include <fstream>

#include "external/eigen-3.4.0/Eigen/Core"
#include "external/eigen-3.4.0/Eigen/Sparse"
#include "external/json/json.hpp"
#include <hdf5/include/H5Cpp.h>

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

    std::ifstream infile = std::ifstream(argv[1]);
    json init_json;
    infile >> init_json;
    infile.close();

    ResponseParams response_params = init_json["response_params"].get<ResponseParams>();
    ComputationParams computation_params = init_json["computation_params"].get<ComputationParams>();
    CubaCommonParams cuba_common_params = init_json["cuba_common_params"].get<CubaCommonParams>();

    CubaIntegrator integrator = init_kernel(response_params, computation_params, cuba_common_params);
    int ncores = 0;
    int pcores = 0;
    cubacores(&ncores, &pcores);

    // Berechnung der Antwortfunktion
    _dvec spectral_data = spectrum(integrator, response_params, computation_params);

    // Speichere Werte im hdf5 Dateiformat
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

    return 0;

}

