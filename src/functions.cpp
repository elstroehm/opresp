#include "functions.hpp"

using namespace QSim;



bool QSim::cmp_for_pair_sort(std::pair<double*, double*> lhs, std::pair<double*, double*> rhs) {
    return *(lhs.first) < *(rhs.first);
}


void PeakData::sort() {
    std::vector<std::pair<double*, double*>> pdata = std::vector<std::pair<double*, double*>>(pos.size());
    _dvec pos_copy = pos;
    _dvec widths_copy = widths;
    for (int i = 0; i < pos.size(); i++) {
        pdata[i] = std::pair<double*, double*>(&pos_copy[i], &widths_copy[i]);
    }
    std::sort(pdata.begin(), pdata.end(), cmp_for_pair_sort);
    for (int i = 0; i < pos.size(); i++) {
        pos[i] = *(pdata[i].first);
        widths[i] = *(pdata[i].second);
    }
}


void PeakData::remove_duplicates_sorted(bool consider_widths) {
    if (pos.size() == 0) {
        return;
    }
    int index = 0;
    if (consider_widths) {
        for (int i = 1; i < pos.size(); i++) {
            if (pos[index] != pos[i] || widths[index] < widths[i]) {
                ++index;
                pos[index] = pos[i];
                widths[index] = widths[i];
            }
        }
    }
    else {
        for (int i = 1; i < pos.size(); i++) {
            if (pos[index] != pos[i]) {
                ++index;
                pos[index] = pos[i];
                widths[index] = widths[i];
            }
        }
    }
    pos.resize(index+1);
    widths.resize(index+1);
}


void PeakData::sort_and_remove_duplicates(bool consider_widths) {
    this->sort();
    this->remove_duplicates_sorted(consider_widths);
}


_dvec PeakData::integration_poles(int nintegrals, const _dvec& x0, const _dvec& x1) {
    TupGen tup_gen = TupGen(_ivec(nintegrals, pos.size()-1));
    _dvec poles = {};
    for (tup_gen.init_initial(); !tup_gen.is_final(); ++tup_gen) {
        for (int k = 0; k < nintegrals; k++) {
            poles.push_back((pos[tup_gen[k]] - x0[k])/(x1[k] - x0[k]));
        }
    }
    for (int k = 0; k < nintegrals; k++) {
        poles.push_back((pos[tup_gen[k]] - x0[k])/(x1[k] - x0[k]));
    }
    return poles;
}


void PeakData::integration_boundaries(int ndim, _dvec& x0, _dvec& x1, int nwidths) {
    if (pos.size() == 0) {
        return;
    }
    double _x0 = pos[0] - nwidths*widths[0];
    double _x1 = pos[0] + nwidths*widths[0];
    for (int i = 0; i < pos.size(); i++) {
        if (pos[i] - nwidths*widths[i] < _x0) {
            _x0 = pos[i] - nwidths*widths[i];
        }
        else if (pos[i] + nwidths*widths[i] > _x1) {
            _x1 = pos[i] + nwidths*widths[i];
        }
    }
    x0 = _dvec(ndim, _x0);
    x1 = _dvec(ndim, _x1);
}


_dvec PeakData::integration_poles_and_boundaries_freq(int nphase_matches, std::vector<bool> mask, const _dvec& shifts, _dvec& b0, _dvec& b1, int nwidths) {
    int ndim = 0;
    const int order = mask.size();
    _ivec index_map = {};
    for (int i = 0; i < order; i++)  {
        if (mask[i]) {
            ++ndim;
            index_map.push_back(i);
        }
    }
    double maxwidth = widths[0];
    for (int i = 0; i < widths.size(); i++) {
        if (widths[i] > maxwidth) {
            maxwidth = widths[i];
        }
    }
    remove_duplicates_sorted(false);
    const int npoles = pos.size();
    _dvec real_poles = {};
    for (int i = 0; i < nphase_matches; i++) {
        TupGen gen = TupGen(_ivec(ndim, npoles - 1));
        for (gen.init_initial(); !gen.is_final(); ++gen) {
            for (int k = 0; k < ndim; k++) {
                real_poles.push_back(pos[gen[k]] - shifts[i*order + index_map[k]]);
            }
        }
        for (int k = 0; k < ndim; k++) {
            real_poles.push_back(pos[gen[k]] - shifts[i*order + index_map[k]]);
        }
    }
    b0.resize(ndim);
    b0.resize(ndim);
    for (int i = 0; i < ndim; i++) {
        b0[i] = real_poles[i];
        b1[i] = real_poles[i];
    }
    for (int offset = 0; offset < ndim; offset++) {
        for (int i = 0; i < real_poles.size(); i+=ndim) {
            if (real_poles[i*order + offset] < b0[offset]) {
                b0[offset] = real_poles[i*order + offset];
            }
            if (real_poles[i*order + offset] > b1[offset]) {
                b1[offset] = real_poles[i*order + offset];
            }
        }
    }
    for (int i = 0; i < ndim; i++) {
        b0[i] -= nwidths*maxwidth;
        b1[i] -= nwidths*maxwidth;
    }
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < real_poles.size(); j+=ndim) {
            real_poles[j*ndim + i] = (real_poles[j*ndim + i] - b0[i])/(b1[i] - b0[i]);
        }
    }
    return real_poles;
}


/***********************************************************************************************/


FnBase::~FnBase() = default;

FnBase::FnBase(const FnBase&) = default;

FnBase& FnBase::operator= (const FnBase&) = default;

FnBase::FnBase(FnType _fn_type) : fn_type{_fn_type} {}

FnBase* FnBase::from_json(const json& ser) {
    FnType fn_type = ser["fn_type"].get<FnType>();
    switch (fn_type) {
        case FnType::Spectrum: {
            _dvec ampl_real = ser["ampl_real"].get<_dvec>();
            _dvec ampl_imag = ser["ampl_imag"].get<_dvec>();
            _cdvec ampl = _cdvec(ampl_real.size());
            for (int i = 0; i < ampl_real.size(); i++) {
                ampl[i] = {ampl_real[i], ampl_imag[i]};
            }
            return new FnSpectrum {
                ampl,
                ser["w0"].get<_dvec>(),
                ser["gamma"].get<_dvec>(),
                ser["is_causal"].get<std::vector<bool>>(),
                ser["eta"].get<double>(),
            };
        }
        case FnType::Interpolation: {
            std::cout << "FnType::Interpolation does not yet have (de-)serialization, returning nullptr!" << std::endl;
            return nullptr;
        }
    }
}


/***********************************************************************************************/


FnSpectrum::~FnSpectrum() = default;

FnSpectrum::FnSpectrum() : FnBase(FnType::Spectrum), ampl{}, w0{}, gamma{}, is_causal{}, eta{0.0} {}

FnSpectrum::FnSpectrum(const FnSpectrum& other) = default;

FnSpectrum::FnSpectrum(FnSpectrum&& other) : FnBase{FnType::Spectrum}, ampl{}, w0{}, gamma{}, is_causal{}, eta(other.eta) {
    std::swap(ampl, other.ampl);
    std::swap(w0, other.w0);
    std::swap(gamma, other.gamma);
    std::swap(is_causal, other.is_causal);
}

FnSpectrum& FnSpectrum::operator= (const FnSpectrum&) = default;

FnSpectrum& FnSpectrum::operator= (FnSpectrum&& other) {
    std::swap(ampl, other.ampl);
    std::swap(w0, other.w0);
    std::swap(gamma, other.gamma);
    std::swap(is_causal, other.is_causal);
    return *this;
}

FnSpectrum::FnSpectrum(const _cdvec& _ampl, const _dvec& _w0, const _dvec& _gamma, const std::vector<bool>& _is_causal, double _eta)
    : FnBase{FnType::Spectrum}, ampl(_ampl), w0(_w0), gamma(_gamma), is_causal(_is_causal), eta(_eta) {}

json FnSpectrum::to_json() const {
    _dvec ampl_real = _dvec(ampl.size());
    _dvec ampl_imag = _dvec(ampl.size());
    for (int i = 0; i < ampl.size(); i++) {
        ampl_real[i] = ampl[i].real();
        ampl_imag[i] = ampl[i].imag();
    }
    return json {
        {"fn_type", FnType::Spectrum},
        {"eta", eta},
        {"ampl_real", ampl_real},
        {"ampl_imag", ampl_imag},
        {"w0", w0},
        {"gamma", gamma},
        {"is_causal", is_causal},
    };
}

FnBase* FnSpectrum::copy_ptr() const {
    return new FnSpectrum {*this};
}

_cd FnSpectrum::eval_t(double arg) const {
    _cd result = {0.0, 0.0};
    for (int i = 0; i < ampl.size(); i++) {
        if (is_causal[i] && arg < 0.0) {
            continue;
        }
        result += ampl[i]*std::exp(-Imd*w0[i]*arg - gamma[i]*std::abs(arg));
    }
    return result;
}

_cd FnSpectrum::eval_w(double arg) const {
    _cd result = {0.0, 0.0};
    for (int i = 0; i < ampl.size(); i++) {
        const double x = std::max(gamma[i], eta);
        const double w = (arg - w0[i]);
        if (is_causal[i]) {
            result += (Imd*w + x)/(w*w + x*x);
        }
        else {
            result += 2.0*x/(w*w + x*x);
        }
    }
    return result;
}

void FnSpectrum::peak_data(PeakData& pdata) const {
    for (int i = 0; i < w0.size(); i++) {
        pdata.pos.push_back(w0[i]);
        pdata.widths.push_back(gamma[i]);
    }
}

void FnSpectrum::add(const FnBase*& _other) {
    if (_other != this && _other != nullptr) {
        const FnSpectrum* other = dynamic_cast<const FnSpectrum*>(_other);
        for (int i = 0; i < other->ampl.size(); i++) {
            ampl.push_back(other->ampl[i]);
            w0.push_back(other->w0[i]);
            gamma.push_back(other->gamma[i]);
            is_causal.push_back(other->is_causal[i]);
            eta = std::min(eta, other->eta);
        }
    }
}

void FnSpectrum::add(FnBase*&& _other) {
    if (_other != this && _other != nullptr) {
        const FnSpectrum* other = dynamic_cast<const FnSpectrum*>(_other);
        for (int i = 0; i < other->ampl.size(); i++) {
            ampl.push_back(other->ampl[i]);
            w0.push_back(other->w0[i]);
            gamma.push_back(other->gamma[i]);
            is_causal.push_back(other->is_causal[i]);
            eta = std::min(eta, other->eta);
        }
    }
}


/***********************************************************************************************/


FnInterpolation::~FnInterpolation() = default;

FnInterpolation::FnInterpolation() : FnBase(FnType::Interpolation), t_evolve_real{}, t_evolve_imag{}, w_evolve_real{}, w_evolve_imag{} {}

FnInterpolation::FnInterpolation(const FnInterpolation& other) = default;

FnInterpolation::FnInterpolation(FnInterpolation&& other)
    : FnBase{FnType::Interpolation}, t_evolve_real{}, t_evolve_imag{}, w_evolve_real{}, w_evolve_imag{} {
    std::swap(t_evolve_real, other.t_evolve_real);
    std::swap(t_evolve_imag, other.t_evolve_imag);
    std::swap(w_evolve_real, other.w_evolve_real);
    std::swap(w_evolve_imag, other.w_evolve_imag);
}

FnInterpolation& FnInterpolation::operator= (const FnInterpolation&) = default;

FnInterpolation& FnInterpolation::operator= (FnInterpolation&& other) {
    std::swap(t_evolve_real, other.t_evolve_real);
    std::swap(t_evolve_imag, other.t_evolve_imag);
    std::swap(w_evolve_real, other.w_evolve_real);
    std::swap(w_evolve_imag, other.w_evolve_imag);
    return *this;
}

FnInterpolation::FnInterpolation(int _nintervals, double _tstart, double _tstop, double _wstart, double _wstop,
        const spline& _t_evolve_real, const spline& _t_evolve_imag, const spline& _w_evolve_real, const spline& _w_evolve_imag)
    : FnBase{FnType::Interpolation}, nintervals(_nintervals), tstart(_tstart), tstop(_tstop), wstart(_wstart), wstop(_wstop),
        t_evolve_real(_t_evolve_real), t_evolve_imag(_t_evolve_imag), w_evolve_real(_w_evolve_real), w_evolve_imag(_w_evolve_imag) {}

FnBase* FnInterpolation::copy_ptr() const {
    return new FnInterpolation{*this};
}

_cd FnInterpolation::eval_t(double arg) const {
    return {t_evolve_real(arg), t_evolve_imag(arg)};
}

_cd FnInterpolation::eval_w(double arg) const {
    return {w_evolve_real(arg), w_evolve_imag(arg)};
}

void FnInterpolation::peak_data(PeakData& pdata) const {
    std::cout << "   void FnInterpolation::peak_data(PeakData& pdata)   noch nicht implementiert" << std::endl;
}

void FnInterpolation::add(const FnBase*& _other) {
    std::cout << "   void FnInterpolation::add(const FnBase*& _other)   noch nicht implementiert" << std::endl;
    if (_other != this && _other != nullptr) {
        const FnInterpolation* other = dynamic_cast<const FnInterpolation*>(_other);
    }
}

void FnInterpolation::add(FnBase*&& _other) {
    std::cout << "   void FnInterpolation::add(FnBase*&& _other)   noch nicht implementiert" << std::endl;
    if (_other != this && _other != nullptr) {
        const FnInterpolation* other = dynamic_cast<const FnInterpolation*>(_other);
    }
}


/***********************************************************************************************/


// FnBase* QSim::operator+ (const Fn& _lhs, const Fn& _rhs) {
//     FnType lhs_type = _lhs->get_fn_type();
//     FnType rhs_type = _rhs->get_fn_type();
//     switch (lhs_type) {
//     case FnType::Interpolation:
//         switch (rhs_type) {
//         case FnType::Interpolation: {
//             std::cout << "Nicht implementiert" << std::endl;
//             return nullptr;
//         }
//         case FnType::Spectrum: {
//             std::cout << "Nicht implementiert" << std::endl;
//             return nullptr;
//         }
//         }
//     case FnType::Spectrum:
//         switch (rhs_type) {
//         case FnType::Interpolation: {
//             std::cout << "Nicht implementiert" << std::endl;
//             return nullptr;
//         }
//         case FnType::Spectrum: {
//             const FnSpectrum* lhs = dynamic_cast<FnSpectrum*>(_lhs.get());
//             const FnSpectrum* rhs = dynamic_cast<FnSpectrum*>(_rhs.get());
//             FnSpectrum* result = new FnSpectrum {};
//             for (int i = 0; i < lhs->ampl.size(); i++) {
//                 result->add(lhs->copy_ptr());
//             }
//             for (int i = 0; i < rhs->ampl.size(); i++) {
//                 result->add(rhs->copy_ptr());
//             }
//             return result;
//         }
//         }
//     }
// }

// void QSim::operator+= (Fn& lhs, const Fn& rhs) {

// }

// FnBase* QSim::operator* (const Fn& lhs, const Fn& rhs) {

// }

// void QSim::operator*= (Fn& lhs, const Fn& rhs) {

// }


/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/


EnvBase::~EnvBase() = default;

EnvBase::EnvBase(const EnvBase&) = default;

EnvBase& EnvBase::operator= (const EnvBase&) = default;

EnvBase::EnvBase(EnvType _env_type) : env_type{_env_type} {};

EnvType EnvBase::get_fn_type() const {return env_type;}

EnvBase* EnvBase::from_json(const json& ser) {
    EnvType env_type = ser["type"].get<EnvType>();
    switch (env_type) {
        case EnvType::Gaussian: {
            return new EnvGaussian {
                ser["ampl"].get<double>(),
                ser["sigma"].get<double>(),
            };
        }
        case EnvType::Chirped: {
            return new EnvChirped {
                ser["ampl"].get<double>(),
                ser["sigma"].get<double>(),
                ser["chirp"].get<double>(),
            };
        }
    }
}


/***********************************************************************************************/


EnvGaussian::~EnvGaussian() = default;

EnvGaussian::EnvGaussian(): EnvBase(EnvType::Gaussian), ampl(1.0), sigma(1.0) {}

EnvGaussian::EnvGaussian(const EnvGaussian& other) = default;

EnvGaussian::EnvGaussian(EnvGaussian&& other) = default;

EnvGaussian& EnvGaussian::operator= (const EnvGaussian& other) = default;

EnvGaussian& EnvGaussian::operator= (EnvGaussian&& other) = default;

EnvGaussian::EnvGaussian(double _ampl, double _sigma)
    : EnvBase(EnvType::Gaussian), ampl(_ampl), sigma(_sigma) {}

json EnvGaussian::to_json() const {
    return json {
        {"type", get_fn_type()},
        {"ampl", ampl},
        {"sigma", sigma},
    };
}

EnvBase* EnvGaussian::copy_ptr() const {
    return new EnvGaussian{*this};
}

_cd EnvGaussian::eval_t(double arg) const {
    return ampl*std::exp(- (arg/sigma)*(arg/sigma));
}

_cd EnvGaussian::eval_w(double arg) const {
    return ampl*std::sqrt(pi)*sigma*std::exp(- 0.25*sigma*arg*sigma*arg);
}

void EnvGaussian::peak_data(PeakData& pdata) const {
    pdata.pos.push_back(0.0);
    pdata.widths.push_back(sigma);
}


/***********************************************************************************************/


EnvChirped::~EnvChirped() = default;

EnvChirped::EnvChirped(): EnvBase(EnvType::Chirped), ampl(1.0), sigma(1.0), chirp(0.0) {}

EnvChirped::EnvChirped(const EnvChirped& other) = default;

EnvChirped::EnvChirped(EnvChirped&& other) = default;

EnvChirped& EnvChirped::operator= (const EnvChirped& other) = default;

EnvChirped& EnvChirped::operator= (EnvChirped&& other) = default;

EnvChirped::EnvChirped(double _ampl, double _sigma, double _chirp)
    : EnvBase(EnvType::Chirped), ampl(_ampl), sigma(_sigma), chirp(_chirp) {}

json EnvChirped::to_json() const {
    return json {
        {"type", get_fn_type()},
        {"ampl", ampl},
        {"sigma", sigma},
        {"chirp", chirp},
    };
}

EnvBase* EnvChirped::copy_ptr() const {
    return new EnvChirped{*this};
}

_cd EnvChirped::eval_t(double arg) const {
    return ampl*std::exp(- (arg/sigma)*(arg/sigma) - Imd*chirp*arg*arg);
}

_cd EnvChirped::eval_w(double arg) const {
    const double sigma2 = sigma*sigma;
    const double sigma2_inv = 1.0/sigma2;
    _cd b = (1.0 - Imd*sigma2*chirp)/(sigma2_inv + sigma2*chirp*chirp);
    return ampl*std::sqrt(pi*b)*std::exp(- 0.25*b*arg*arg);
}

void EnvChirped::peak_data(PeakData& pdata) const {
    pdata.pos.push_back(0.0);
    pdata.widths.push_back(sigma);
}








