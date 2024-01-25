#include "integration.hpp"

using namespace QSim;



CubaParamsVegas CubaParamsVegas::with_default_args() {
    return CubaParamsVegas {
        1, 10000, 10000, 10000, 0
    };
}

CubaParamsSuave CubaParamsSuave::with_default_args() {
    return CubaParamsSuave {
        1, 1000, 1000, 10
    };
}

CubaParamsDivonne CubaParamsDivonne::with_default_args() {
    return CubaParamsDivonne {
        0, -10, -10, 1, 4, 0.0, 1.0, 10.0, 0, 0, _dvec(), 0, nullptr
    };
}

CubaParamsCuhre CubaParamsCuhre::with_default_args() {
    return CubaParamsCuhre {
        10
    };
}



CubaResult QSim::cuba_vegas_integrate(CubaIntegrand integrand, void* userdata,
    const CubaParamsVegas& params, const CubaCommonParams& cuba_args) {

    CubaResult result = CubaResult {
        integrand.ncomp, 0, 0, 0, _dvec(integrand.ncomp), _dvec(integrand.ncomp), _dvec(integrand.ncomp)
    };
    Vegas (integrand.ndim, integrand.ncomp, integrand.integrand, userdata, integrand.nvec,
        cuba_args.epsrel, cuba_args.epsabs, cuba_args.flags, params.seed, cuba_args.mineval, cuba_args.maxeval,
        params.nstart, params.nincrease, params.nbatch, params.gridno,
        nullptr, nullptr,
        &result.neval, &result.fail, 
        &result.integral[0], &result.error[0], &result.prob[0]
    );
    return result;
}

CubaResult QSim::cuba_suave_integrate(CubaIntegrand integrand, void* userdata,
    const CubaParamsSuave& params, const CubaCommonParams& cuba_args) {

    CubaResult result = CubaResult {
        integrand.ncomp, 0, 0, 0, _dvec(integrand.ncomp), _dvec(integrand.ncomp), _dvec(integrand.ncomp)
    };
    Suave (integrand.ndim, integrand.ncomp, integrand.integrand, userdata, integrand.nvec,
        cuba_args.epsrel, cuba_args.epsabs, cuba_args.flags, params.seed, cuba_args.mineval, cuba_args.maxeval,
        params.nnew, params.nmin, params.flatness,
        nullptr, nullptr,
        &result.nregions, &result.neval, &result.fail, 
        &result.integral[0], &result.error[0], &result.prob[0]
    );
    return result;
}

CubaResult QSim::cuba_divonne_integrate(CubaIntegrand integrand, void* userdata,
    const CubaParamsDivonne& params, const CubaCommonParams& cuba_args) {

    CubaResult result = CubaResult {
        integrand.ncomp, 0, 0, 0, _dvec(integrand.ncomp), _dvec(integrand.ncomp), _dvec(integrand.ncomp)
    };
    _dvec xgiven = params.xgiven;
    xgiven.resize(xgiven.size() + params.nextra*params.ldxgiven);
    Divonne (integrand.ndim, integrand.ncomp, integrand.integrand, userdata, integrand.nvec,
        cuba_args.epsrel, cuba_args.epsabs, cuba_args.flags, params.seed, cuba_args.mineval, cuba_args.maxeval,
        params.key1, params.key2, params.key3, params.maxpass, params.border,
        params.maxchisq, params.mindeviation,
        params.ngiven, params.ldxgiven, &xgiven[0],
        params.nextra, params.peakfinder,
        nullptr, nullptr,
        &result.nregions, &result.neval, &result.fail, 
        &result.integral[0], &result.error[0], &result.prob[0]
    );
    return result;
}

CubaResult QSim::cuba_cuhre_integrate(CubaIntegrand integrand, void* userdata,
    const CubaParamsCuhre& params, const CubaCommonParams& cuba_args) {

    CubaResult result = CubaResult {
        integrand.ncomp, 0, 0, 0, _dvec(integrand.ncomp), _dvec(integrand.ncomp), _dvec(integrand.ncomp)
    };
    Cuhre (integrand.ndim, integrand.ncomp, integrand.integrand, userdata, integrand.nvec,
        cuba_args.epsrel, cuba_args.epsabs, cuba_args.flags, cuba_args.mineval, cuba_args.maxeval,
        params.key,
        nullptr, nullptr,
        &result.nregions, &result.neval, &result.fail, 
        &result.integral[0], &result.error[0], &result.prob[0]
    );
    return result;
}



void QSim::cuba_vegas_integrate_inplace(CubaIntegrand integrand, void* userdata,
    CubaParamsVegas params, CubaCommonParams cuba_args, CubaResult* result) {

    result->ncomp = integrand.ncomp;
    result->nregions = 0;
    result->neval = 0;
    result->fail = 0;
    result->integral = _dvec(integrand.ncomp);
    result->prob = _dvec(integrand.ncomp);
    result->error = _dvec(integrand.ncomp);
    if (integrand.ndim == 0) {
        integrand.integrand(&integrand.ndim, nullptr, &integrand.ncomp, &(result->integral[0]), userdata);
        return;
    }
    Vegas (integrand.ndim, integrand.ncomp, integrand.integrand, userdata, integrand.nvec,
        cuba_args.epsrel, cuba_args.epsabs, cuba_args.flags, params.seed, cuba_args.mineval, cuba_args.maxeval,
        params.nstart, params.nincrease, params.nbatch, params.gridno,
        nullptr, nullptr,
        &result->neval, &result->fail, 
        &result->integral[0], &result->error[0], &result->prob[0]
    );
}

void QSim::cuba_suave_integrate_inplace(CubaIntegrand integrand, void* userdata,
    CubaParamsSuave params, CubaCommonParams cuba_args, CubaResult* result) {

    result->ncomp = integrand.ncomp;
    result->nregions = 0;
    result->neval = 0;
    result->fail = 0;
    result->integral = _dvec(integrand.ncomp);
    result->prob = _dvec(integrand.ncomp);
    result->error = _dvec(integrand.ncomp);
    if (integrand.ndim == 0) {
        integrand.integrand(&integrand.ndim, nullptr, &integrand.ncomp, &(result->integral[0]), userdata);
        return;
    }
    Suave (integrand.ndim, integrand.ncomp, integrand.integrand, userdata, integrand.nvec,
        cuba_args.epsrel, cuba_args.epsabs, cuba_args.flags, params.seed, cuba_args.mineval, cuba_args.maxeval,
        params.nnew, params.nmin, params.flatness,
        nullptr, nullptr,
        &result->nregions, &result->neval, &result->fail, 
        &result->integral[0], &result->error[0], &result->prob[0]
    );
}

void QSim::cuba_divonne_integrate_inplace(CubaIntegrand integrand, void* userdata,
    CubaParamsDivonne params, CubaCommonParams cuba_args, CubaResult* result) {

    result->ncomp = integrand.ncomp;
    result->nregions = 0;
    result->neval = 0;
    result->fail = 0;
    result->integral = _dvec(integrand.ncomp);
    result->prob = _dvec(integrand.ncomp);
    result->error = _dvec(integrand.ncomp);
    _dvec xgiven = params.xgiven;
    xgiven.resize(xgiven.size() + params.nextra*params.ldxgiven);
    if (integrand.ndim == 0) {
        integrand.integrand(&integrand.ndim, nullptr, &integrand.ncomp, &(result->integral[0]), userdata);
        return;
    }
    Divonne (integrand.ndim, integrand.ncomp, integrand.integrand, userdata, integrand.nvec,
        cuba_args.epsrel, cuba_args.epsabs, cuba_args.flags, params.seed, cuba_args.mineval, cuba_args.maxeval,
        params.key1, params.key2, params.key3, params.maxpass, params.border,
        params.maxchisq, params.mindeviation,
        params.ngiven, params.ldxgiven, &xgiven[0],
        params.nextra, params.peakfinder,
        nullptr, nullptr,
        &result->nregions, &result->neval, &result->fail, 
        &result->integral[0], &result->error[0], &result->prob[0]
    );
}

void QSim::cuba_cuhre_integrate_inplace(CubaIntegrand integrand, void* userdata,
    CubaParamsCuhre params, CubaCommonParams cuba_args, CubaResult* result) {

    result->ncomp = integrand.ncomp;
    result->nregions = 0;
    result->neval = 0;
    result->fail = 0;
    result->integral = _dvec(integrand.ncomp);
    result->prob = _dvec(integrand.ncomp);
    result->error = _dvec(integrand.ncomp);
    if (integrand.ndim == 0) {
        integrand.integrand(&integrand.ndim, nullptr, &integrand.ncomp, &(result->integral[0]), userdata);
        return;
    }
    Cuhre (integrand.ndim, integrand.ncomp, integrand.integrand, userdata, integrand.nvec,
        cuba_args.epsrel, cuba_args.epsabs, cuba_args.flags, cuba_args.mineval, cuba_args.maxeval,
        params.key,
        nullptr, nullptr,
        &result->nregions, &result->neval, &result->fail, 
        &result->integral[0], &result->error[0], &result->prob[0]
    );
}

