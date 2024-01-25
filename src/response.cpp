#include "response.hpp"

using namespace QSim;


void QSim::to_json(json& ser, const ComputationParams& val) {
    ser = json {
        {"nthreads", val.nthreads},
        {"lower_borders", val.lower_borders},
        {"upper_borders", val.upper_borders},
        {"nintervals", val.nintervals},
        {"cuba_alg", val.cuba_alg},
        {"repr", val.repr},
    };
}


void QSim::from_json(const json& ser, ComputationParams& val) {
    val = ComputationParams {
        ser["nthreads"].get<int>(),
        ser["lower_borders"].get<_dvec>(),
        ser["upper_borders"].get<_dvec>(),
        ser["nintervals"].get<_ivec>(),
        ser["cuba_alg"].get<CubaAlg>(),
        ser["repr"].get<Representation>(),
    };
}


/********************************************************************************************/


CubaIntegrator QSim::init_kernel(ResponseParams& response_params, ComputationParams computation_params, CubaCommonParams cuba_comm) {
    int order = response_params.interaction.ncontacts;
    int ninsertions = 0;
    for (int i = 0; i < order; i++) {
        if (response_params.arg_domains[i] == Domain::Freq) {
            ++ninsertions;
        }
    }
    int ndim;
    PeakData pdata = PeakData();
    integrand_t integrand;
    peakfinder_t peaksearch_ptr = nullptr;
    _dvec integration_poles = {};
    switch (computation_params.repr) {
        case Representation::DysonFreq: {
            ndim = order - ninsertions;
            std::vector<bool> mask = std::vector<bool>(order, false);
            for (int i = 0; i < order; i++)  {
                if (response_params.arg_domains[i] == Domain::Time) {
                    mask[i] = true;
                }
            }
            pdata = response_params.propagator.peak_data();
            pdata.sort_and_remove_duplicates(true);
            integration_poles = pdata.integration_poles_and_boundaries_freq(response_params.interaction.nphase_matches, mask, response_params.interaction.w_sum, 
                    response_params.x0, response_params.x1, 20);
            integrand = dyson_freq_repr;
            peaksearch_ptr = peaksearch_dyson_freq;
            break;
        }
        case Representation::DysonTime: {
            ndim = order;
            pdata = response_params.interaction.peak_data();
            pdata.sort_and_remove_duplicates(true);
            pdata.integration_boundaries(ndim, response_params.x0, response_params.x1, 5);
            pdata.remove_duplicates_sorted(false);
            integration_poles = pdata.integration_poles(ndim, response_params.x0, response_params.x1);
            integrand = dyson_time_repr;
            peaksearch_ptr = peaksearch_dyson_time;
            break;
        }
        case Representation::SemiImpulsive: {
            ndim = 0;
            integrand = semi_impulsive_repr;
            break;
        }
        // case Representation::SemiImpulsiveTime: {
        //     ndim = ninsertions;
        //     integrand = semi_impulsive_time_repr;
        //     break;
        // }
    }
    response_params.integration_poles = integration_poles;
    response_params.ndim = ndim;
    if (ndim != 0) {
        response_params.npoles = integration_poles.size()/ndim;
    }
    else {
        response_params.npoles = 0;
    }
    CubaIntegrand cuba_integrand = CubaIntegrand {integrand, ndim, 2, 1};
    integrand = nullptr;
    CubaIntegrator integrator;
    switch (computation_params.cuba_alg) {
        case CubaAlg::Vegas: {
            CubaParamsVegas cuba_params = CubaParamsVegas::with_default_args();
            integrator = std::bind(cuba_vegas_integrate_inplace,
                cuba_integrand, std::placeholders::_1, cuba_params, cuba_comm, std::placeholders::_2);
            break;
        }
        case CubaAlg::Suave: {
            CubaParamsSuave cuba_params = CubaParamsSuave::with_default_args();
            integrator = std::bind(cuba_suave_integrate_inplace,
                cuba_integrand, std::placeholders::_1, cuba_params, cuba_comm, std::placeholders::_2);
            break;
        }
        case CubaAlg::Divonne: {
            CubaParamsDivonne cuba_params = CubaParamsDivonne::with_default_args();
            if (ndim != 0) {
                // cuba_params.ngiven = integration_poles.size()/ndim;
                // cuba_params.xgiven = integration_poles;
                cuba_params.nextra = integration_poles.size()/ndim;
                cuba_params.ldxgiven = ndim;
                cuba_params.peakfinder = peaksearch_ptr;
            }
            integrator = std::bind(cuba_divonne_integrate_inplace,
                cuba_integrand, std::placeholders::_1, cuba_params, cuba_comm, std::placeholders::_2);
            break;
        }
        case CubaAlg::Cuhre: {
            CubaParamsCuhre cuba_params = CubaParamsCuhre::with_default_args();
            integrator = std::bind(cuba_cuhre_integrate_inplace,
                cuba_integrand, std::placeholders::_1, cuba_params, cuba_comm, std::placeholders::_2);
            break;
        }
    }
    return integrator;
}


void QSim::update_progress(int& counter, int max_count, int& prev_progress, std::mutex& mtx) {
    mtx.lock();
    ++counter;
    int* progress = new int((counter*100)/max_count);
    if (*progress > prev_progress) {
        std::cout << "Fortschritt: " << *progress << "%" << std::endl;
        prev_progress = *progress;
    }
    mtx.unlock();
}


void QSim::compute_response(CubaIntegrator integrator, ResponseParams response_params, _dvec args, CubaResult* result_ptr, int npoints, std::function<void()> update) {
    const int order = response_params.interaction.ncontacts;
    CubaResult* result = new CubaResult{
        0, 0, 0, 0, _dvec(2, 0.0), _dvec(2, 0.0), _dvec(2, 0.0)
    };
    for (int i = 0; i < npoints; i++) {
        for (int j = 0; j < order; j++) {
            response_params.args[j] = args[i*order + j];
        }
        for (int j = 0; j < response_params.interaction.nphase_matches; j++) {
            response_params.phase_match_state = j;
            integrator(static_cast<void*>(&response_params), result);
            result_ptr[i].ncomp = result->ncomp;
            result_ptr[i].fail = result->fail;
            result_ptr[i].neval = result->neval;
            result_ptr[i].nregions = result->nregions;
            result_ptr[i].integral[0] += result->integral[0];
            result_ptr[i].integral[1] += result->integral[1];
            result_ptr[i].error[0] += result->error[0];
            result_ptr[i].error[1] += result->error[1];
            result_ptr[i].prob[0] += result->prob[0];
            result_ptr[i].prob[1] += result->prob[1];
        }
        update();
    }
    delete result;
}


_dvec QSim::spectrum(CubaIntegrator integrator, ResponseParams response_params, ComputationParams computation_params) {
    int order = response_params.interaction.ncontacts;
    const _ivec nintervals = computation_params.nintervals;
    const _dvec lower_bounds = computation_params.lower_borders;
    const _dvec upper_bounds = computation_params.upper_borders;
    const int nthreads = computation_params.nthreads;
    _dvec stepsizes = _dvec(order, 0.0);
    int npoints = 1;
    for (int i = 0; i < order; i++) {
        npoints *= nintervals[i] + 1;
        if (nintervals[i] != 0) {
            stepsizes[i] = (upper_bounds[i] - lower_bounds[i])/nintervals[i];
        }
    }
    std::vector<CubaResult> results = std::vector<CubaResult>(npoints, CubaResult {
        0, 0, 0, 0, _dvec(2, 0.0), _dvec(2, 0.0), _dvec(2, 0.0)
    });
    TupGen tup_gen = TupGen(nintervals);
    const int ppt2 = npoints / nthreads;
    const int nt1 = npoints - ppt2*nthreads;
    const int nt2 = nthreads - nt1;
    const int ppt1 = ppt2 + 1;
    std::vector<std::thread> workbench = std::vector<std::thread>();
    _dvec args1 = _dvec(ppt1*order, 0.0);
    int index = 0;
    const int data_dimension = order + 4;
    _dvec data = _dvec(npoints*data_dimension);
    std::mutex mtx;
    int counter = 0;
    int prev_progress = 0;
    std::function<void()> update = std::bind(update_progress, std::ref(counter), npoints, std::ref(prev_progress), std::ref(mtx));
    for (int i = 0; i < nt1; i++) {
        int start = index;
        for (int j = 0; j < ppt1; j++) {
            for (int k = 0; k < order; k++) {
                args1[j*order + k] = lower_bounds[k] + tup_gen.state[k]*stepsizes[k];
                data[i*ppt1*data_dimension + j*data_dimension + k] = args1[j*order + k];
            }
            ++tup_gen;
            ++index;
        }
        workbench.push_back(std::thread(compute_response, integrator, response_params, args1, &results[start], ppt1, update));
    }
    const int data_offset = nt1*ppt1*data_dimension;
    _dvec args2 = _dvec(ppt2*order, 0.0);
    for (int i = 0; i < nt2; i++) {
        int start = index;
        for (int j = 0; j < ppt2; j++) {
            for (int k = 0; k < order; k++) {
                args2[j*order + k] = lower_bounds[k] + tup_gen.state[k]*stepsizes[k];
                data[data_offset + i*ppt2*data_dimension + j*data_dimension + k] = args2[j*order + k];
            }
            ++tup_gen;
            ++index;
        }
        workbench.push_back(std::thread(compute_response, integrator, response_params, args2, &results[start], ppt2, update));
    }
    for (int i = 0; i < nthreads; i++) {
        workbench[i].join();
    }
    for (int i = 0; i < npoints; i++) {
        data[i*data_dimension + order + 0] = results[i].integral[0];
        data[i*data_dimension + order + 1] = results[i].integral[1];
        data[i*data_dimension + order + 2] = results[i].error[0];
        data[i*data_dimension + order + 3] = results[i].error[1];
    }
    return data;
}




