#pragma once

#include <complex>
#include <vector>
#include <iostream>
#include <functional>
#include <memory>
#include <thread>
#include <algorithm>
#include <mutex>

#include "../external/eigen-3.4.0/Eigen/Core"
#include "../external/eigen-3.4.0/Eigen/Sparse"
#include "../external/json/json.hpp"

#include "definitions.hpp"
#include "algebra.hpp"
#include "functions.hpp"
#include "propagator.hpp"
#include "interaction.hpp"
#include "integration.hpp"
#include "response_kernels.hpp"
#include "eigen_serialization.hpp"

using namespace Eigen;
using json = nlohmann::json;



namespace QSim {


    struct ComputationParams {
        int nthreads;
        _dvec lower_borders;
        _dvec upper_borders;
        _ivec nintervals;
        CubaAlg cuba_alg;
        Representation repr;
    };

    void to_json(json& ser, const ComputationParams& val);
    void from_json(const json& ser, ComputationParams& val);


    CubaIntegrator init_kernel(ResponseParams& response_params, ComputationParams computation_params, CubaCommonParams cuba_comm);
    

    void update_progress(int& counter, int max_count, int& prev_progress, std::mutex& mtx);
    void compute_response(CubaIntegrator integrator, ResponseParams response_params, _dvec args, CubaResult* result_ptr, int npoints, std::function<void()> update);
    _dvec spectrum(CubaIntegrator integrator, ResponseParams response_params, ComputationParams computation_params);

}