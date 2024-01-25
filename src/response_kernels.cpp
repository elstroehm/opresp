#include "response_kernels.hpp"

using namespace QSim;



ResponseParams::~ResponseParams() = default;

ResponseParams::ResponseParams() = default;

ResponseParams::ResponseParams(const ResponseParams&) = default;

ResponseParams& ResponseParams::operator= (const ResponseParams&) = default;

ResponseParams::ResponseParams(ResponseParams&&) = default;

ResponseParams& ResponseParams::operator= (ResponseParams&&) = default;

ResponseParams::ResponseParams(const SparseVector<_cd> _observable, const SparseMatrix<_cd>& _qform, const _dvec& _x0, const _dvec& _x1, const Propagator& _propagator,
        const Interaction& _interaction, const std::vector<Domain>& _arg_domains)
    : observable(_observable), qform(_qform), x0(_x0), x1(_x1), propagator(_propagator), interaction(_interaction), args(_arg_domains.size()),
        arg_domains(_arg_domains), phase_match_state{0}, integration_poles{}, npoles{0}, ndim{0} {}


void QSim::to_json(json& ser, const ResponseParams& val) {
    ser = json {
        {"qform", val.qform},
        {"observable", val.observable},
        {"x0", val.x0},
        {"x1", val.x1},
        {"propagator", val.propagator},
        {"interaction", val.interaction},
        {"arg_domains", val.arg_domains},
    };
}


void QSim::from_json(const json& ser, ResponseParams& val) {
    val = ResponseParams {
        ser["qform"].get<SparseMatrix<_cd>>(),
        ser["observable"].get<SparseVector<_cd>>(),
        ser["x0"].get<_dvec>(),
        ser["x1"].get<_dvec>(),
        ser["propagator"].get<Propagator>(),
        ser["interaction"].get<Interaction>(),
        ser["arg_domains"].get<std::vector<Domain>>(),
    };
}


/********************************************************************************************/


int QSim::dyson_freq_repr(const int* _ndim, const double i_args[], const int* _ncomp,
    double _val[], void* _params) {

    ResponseParams* params = static_cast<ResponseParams*>(_params);
    const int order = params->interaction.ncontacts;
    SparseVector<_cd> evolve = SparseVector<_cd>(params->propagator.rows);
    evolve.insert(0,0) = {1.0, 0.0};
    _cd phase = 1.0;
    double volume = 1.0;
    double w_i = 0.0;
    double w_f = 0.0;
    int iindex = 0;
    for (int i = 0; i < order; i++) {
        const double w_sum = params->interaction.w_sum[params->phase_match_state*order + i];
        const int sign = params->interaction.phase_match_signatures[params->phase_match_state*order + i];
        if (params->arg_domains[i] == Domain::Freq && sign == 1) {
            w_f = params->args[i];
            phase *= -Imd * params->interaction.envelopes[i]->eval_w(w_f - w_i);
        }
        else if (params->arg_domains[i] == Domain::Freq && sign == -1) {
            w_f = params->args[i];
            phase *= -Imd * std::conj(params->interaction.envelopes[i]->eval_w(w_i - w_f));
        }
        else if (params->arg_domains[i] == Domain::Time && sign == 1) {
            w_f = (params->x1[iindex] - params->x0[iindex])*i_args[iindex] + params->x0[iindex];
            volume *= (params->x1[iindex] - params->x0[iindex])/twopi;
            phase *= -Imd * std::exp(- Imd * params->args[i] * w_f) * params->interaction.envelopes[i]->eval_w(w_f - w_i);
            ++iindex;
        }
        else if (params->arg_domains[i] == Domain::Time && sign == -1) {
            w_f = (params->x1[iindex] - params->x0[iindex])*i_args[iindex] + params->x0[iindex];
            volume *= (params->x1[iindex] - params->x0[iindex])/twopi;
            phase *= -Imd * std::exp(- Imd * params->args[i] * w_f) * std::conj(params->interaction.envelopes[i]->eval_w(w_i - w_f));
            ++iindex;
        }
        evolve =    ( params->propagator.eval_w(w_f + w_sum)
                    * params->interaction.couplings
                    * evolve ).eval();
        w_i = w_f;
    }
    _cd result = Imd * volume * phase * (params->observable.transpose()*params->qform*evolve).eval().coeffRef(0,0);
    if (std::isnan(result.real()) || std::isnan(result.imag()) || std::isinf(result.real()) || std::isinf(result.imag())) {
        _val[0] = 0.0;
        _val[1] = 0.0;
        std::cout << "unzulässiger Wert:   " << result << std::endl;
    }
    else {
        _val[0] = result.real();
        _val[1] = result.imag();
    }
    return 0;
}


int QSim::dyson_time_repr(const int* _ndim, const double i_args[], const int* _ncomp,
    double _val[], void* _params) {

    ResponseParams* params = static_cast<ResponseParams*>(_params);
    const int order = params->interaction.ncontacts;
    SparseVector<_cd> evolve = SparseVector<_cd>(params->propagator.rows);
    evolve.insert(0,0) = {1.0, 0.0};
    _cd phase = 1.0;
    double volume = 1.0;
    double tau_f = 0.0;
    double tau_i = 0.0;
    for (int i = 0; i < order; i++) {
        const double w_sum = params->interaction.w_sum[params->phase_match_state*order + i];
        const int sign = params->interaction.phase_match_signatures[params->phase_match_state*order + i];
        volume *= (params->x1[i] - params->x0[i]);
        tau_i = (params->x1[i] - params->x0[i])*i_args[i] + params->x0[i];
        tau_f = 0.0;
        if (i < order - 1) {
            tau_f = (params->x1[i+1] - params->x0[i+1])*i_args[i+1] + params->x0[i+1];
        }
        if (params->arg_domains[i] == Domain::Freq && sign == 1) {
            phase *= -Imd * std::exp(-Imd*(params->args[i] + w_sum)*(tau_f - tau_i) - Imd*params->interaction.w0[i]*tau_i)
                    * params->interaction.envelopes[i]->eval_t(tau_i);
            evolve = ( params->propagator.eval_w(params->args[i] + w_sum) * params->interaction.couplings * evolve ).eval();
        }
        else if (params->arg_domains[i] == Domain::Freq && sign == -1) {
            phase *= -Imd * std::exp(-Imd*(params->args[i] + w_sum)*(tau_f - tau_i) + Imd*params->interaction.w0[i]*tau_i)
                    * std::conj(params->interaction.envelopes[i]->eval_t(tau_i));
            evolve = ( params->propagator.eval_w(params->args[i] + w_sum) * params->interaction.couplings * evolve ).eval();
        }
        else if (params->arg_domains[i] == Domain::Time && sign == 1) {
            phase *= -Imd * std::exp(Imd*params->args[i]*w_sum - Imd*params->interaction.w0[i]*tau_i)
                    * params->interaction.envelopes[i]->eval_t(tau_i);
            evolve = ( params->propagator.eval_t(params->args[i] + tau_f - tau_i) * params->interaction.couplings * evolve ).eval();
        }
        else if (params->arg_domains[i] == Domain::Time && sign == -1) {
            phase *= -Imd * std::exp(Imd*params->args[i]*w_sum + Imd*params->interaction.w0[i]*tau_i)
                    * std::conj(params->interaction.envelopes[i]->eval_t(tau_i));
            evolve = ( params->propagator.eval_t(params->args[i] + tau_f - tau_i) * params->interaction.couplings * evolve ).eval();
        }
    }
    _cd result = Imd * volume * phase * (params->observable.transpose()*params->qform*evolve).eval().coeffRef(0,0);
    if (std::isnan(result.real()) || std::isnan(result.imag()) || std::isinf(result.real()) || std::isinf(result.imag())) {
        _val[0] = 0.0;
        _val[1] = 0.0;
        std::cout << "unzulässiger Wert:   " << result << std::endl;
    }
    else {
        _val[0] = result.real();
        _val[1] = result.imag();
    }
    return 0;
}


int QSim::semi_impulsive_repr(const int* _ndim, const double i_args[], const int* _ncomp,
    double _val[], void* _params) {
    
    ResponseParams* params = static_cast<ResponseParams*>(_params);
    const int order = params->interaction.ncontacts;
    SparseVector<_cd> evolve = SparseVector<_cd>(params->propagator.rows);
    evolve.insert(0,0) = {1.0, 0.0};
    _cd phase = 1.0;
    for (int i = 0; i < order; i++) {
        const double w_sum = params->interaction.w_sum[params->phase_match_state*order + i];
        phase *= - Imd;
        if (params->arg_domains[i] == Domain::Freq) {
            evolve = (params->propagator.eval_w(params->args[i] + w_sum) * params->interaction.couplings * evolve).eval();
        }
        else if (params->arg_domains[i] == Domain::Time) {
            phase *= std::exp(Imd * params->args[i] * w_sum);
            evolve = (params->propagator.eval_t(params->args[i]) * params->interaction.couplings * evolve).eval();
        }
    }
    _cd result = Imd * phase * (params->observable.transpose()*params->qform*evolve).eval().coeffRef(0,0);
    // _cd result = Imd * phase * params->observable.dot(evolve);
    if (std::isnan(result.real()) || std::isnan(result.imag()) || std::isinf(result.real()) || std::isinf(result.imag())) {
        _val[0] = 0.0;
        _val[1] = 0.0;
        std::cout << "unzulässiger Wert" << std::endl;
    }
    else {
        _val[0] = result.real();
        _val[1] = result.imag();
    }
    return 0;
}


void QSim::peaksearch_dyson_freq(const int* ndim, const double* borders, int* npeaks, double* pos, void* userdata) {
    ResponseParams* params = static_cast<ResponseParams*>(userdata);
    const int pstate = params->phase_match_state;
    const int nintegrals = *ndim;
    const int ppm = params->npoles/params->interaction.nphase_matches;
    const int offs = pstate*ppm*nintegrals;
    int npoles = 0;
    int counter = 0;
    while (npoles < params->npoles && npoles < *npeaks && counter < ppm) {
        for (int j = 0; j < nintegrals; j++) {
            if (params->integration_poles[offs + counter*nintegrals + j] > borders[2*j] 
            && params->integration_poles[offs + counter*nintegrals + j] < borders[2*j + 1]) {
                pos[npoles*nintegrals + j] = params->integration_poles[offs + counter*nintegrals + j];
            }
            else {
                --npoles;
                break;
            }
        }
        ++npoles;
        ++counter;
    }
    *npeaks = npoles;
}


void QSim::peaksearch_dyson_time(const int* ndim, const double* borders, int* npeaks, double* pos, void* userdata) {
    ResponseParams* params = static_cast<ResponseParams*>(userdata);
    const int nintegrals = *ndim;
    int npoles = 0;
    int counter = 0;
    while (npoles < params->npoles && npoles < *npeaks && counter < params->npoles) {
        for (int j = 0; j < nintegrals; j++) {
            if (params->integration_poles[counter*nintegrals + j] > borders[2*j] 
            && params->integration_poles[counter*nintegrals + j] < borders[2*j + 1]) {
                pos[npoles*nintegrals + j] = params->integration_poles[counter*nintegrals + j];
            }
            else {
                --npoles;
                break;
            }
        }
        ++npoles;
        ++counter;
    }
    *npeaks = npoles;
}
