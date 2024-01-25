#pragma once

#include <complex>
#include <vector>
#include <iostream>
#include <memory>
#include <algorithm>
#include <utility>
#include <algorithm>

#include "../external/json/json.hpp"
#include "../external/boost/math/interpolators/cardinal_cubic_b_spline.hpp"

#include "definitions.hpp"
#include "algebra.hpp"

using json = nlohmann::json;
using spline = boost::math::interpolators::cardinal_cubic_b_spline<double>;



namespace QSim {


    bool cmp_for_pair_sort(std::pair<double*, double*> lhs, std::pair<double*, double*> rhs);

    struct PeakData {
        _dvec pos;
        _dvec widths;

        void sort();
        void remove_duplicates_sorted(bool consider_widths);
        void sort_and_remove_duplicates(bool consider_widths);
        _dvec integration_poles(int ndim, const _dvec& x0, const _dvec& x1);
        void integration_boundaries(int ndim, _dvec& x0, _dvec& x1, int nwidths);
        _dvec integration_poles_and_boundaries_freq(int nphase_matches, std::vector<bool> mask, const _dvec& w_sum, _dvec& b0, _dvec& b1, int nwidths);
    };


/***********************************************************************************************/


    enum class FnType {
        Spectrum,
        Interpolation,
    };

    NLOHMANN_JSON_SERIALIZE_ENUM( FnType, {
        {FnType::Spectrum, "Spectrum"},
        {FnType::Interpolation, "Interpolation"},
    })


    class FnBase {
    private:
        FnType fn_type;    // Gibt zwecks Serialisierung den tatsächlichen Funktionstyp an (z.B. FnGaussian)

    public:
        virtual ~FnBase();
        FnBase(const FnBase&);
        FnBase& operator= (const FnBase&);
        FnBase(FnType _fn_type);

        FnType get_fn_type() const;
        static FnBase* from_json(const json& j);

        virtual json to_json() const = 0;
        virtual FnBase* copy_ptr() const = 0;
        virtual _cd eval_t(double t) const = 0;
        virtual _cd eval_w(double w) const = 0;
        virtual void peak_data(PeakData& pdata) const = 0;
        virtual void add(const FnBase*&) = 0;
        virtual void add(FnBase*&&) = 0;
    };


    class FnSpectrum : public FnBase {
    public:
        _cdvec ampl;
        _dvec w0;   // Phase in Time Domain, Shift in Frequency Domain
        _dvec gamma;
        std::vector<bool> is_causal;
        double eta;

        ~FnSpectrum();
        FnSpectrum();
        FnSpectrum(const FnSpectrum&);
        FnSpectrum(FnSpectrum&&);
        FnSpectrum& operator= (const FnSpectrum&);
        FnSpectrum& operator= (FnSpectrum&&);
        FnSpectrum(const _cdvec& _ampl, const _dvec& w0, const _dvec& gamma, const std::vector<bool>& is_causal, double eta);

        virtual json to_json() const override;
        virtual FnBase* copy_ptr() const override;
        virtual _cd eval_t(double arg) const override;
        virtual _cd eval_w(double arg) const override;
        virtual void peak_data(PeakData& pdata) const override;
        virtual void add(const FnBase*&) override;
        virtual void add(FnBase*&&) override;
    };

    
    class FnInterpolation : public FnBase {
    public:
        int nintervals;
        double tstart;
        double tstop;
        double wstart;
        double wstop;
        spline t_evolve_real;
        spline t_evolve_imag;
        spline w_evolve_real;
        spline w_evolve_imag;

        ~FnInterpolation();
        FnInterpolation();
        FnInterpolation(const FnInterpolation&);
        FnInterpolation(FnInterpolation&&);
        FnInterpolation& operator= (const FnInterpolation&);
        FnInterpolation& operator= (FnInterpolation&&);
        FnInterpolation(int _nintervals, double _tstart, double _tstop, double _wstart, double _wstop,
            const spline& _t_evolve_real, const spline& _t_evolve_imag, const spline& _w_evolve_real, const spline& _w_evolve_imag);

        virtual json to_json() const override;
        virtual FnBase* copy_ptr() const override;
        virtual _cd eval_t(double arg) const override;
        virtual _cd eval_w(double arg) const override;
        virtual void peak_data(PeakData& pdata) const override;
        virtual void add(const FnBase*&) override;
        virtual void add(FnBase*&&) override;
    };


    typedef std::unique_ptr<FnBase> Fn;


    // bool operator== (const Fn& _lhs, const Fn& _rhs);
    // FnBase* operator+ (const Fn& lhs, const Fn& rhs);
    // void operator+= (Fn& lhs, const Fn& rhs);
    // FnBase* operator* (const Fn& lhs, const Fn& rhs);
    // void operator*= (Fn& lhs, const Fn& rhs);


/***********************************************************************************************/


    enum class EnvType {
        Gaussian,
        Chirped,
    };

    NLOHMANN_JSON_SERIALIZE_ENUM( EnvType, {
        {EnvType::Gaussian, "Gaussian"},
        {EnvType::Chirped, "Chirped"},
    })


    class EnvBase {
    private:
        EnvType env_type;

    public:
        virtual ~EnvBase();
        EnvBase(const EnvBase&);
        EnvBase& operator= (const EnvBase&);
        EnvBase(EnvType _env_type);

        EnvType get_fn_type() const;
        static EnvBase* from_json(const json& j);

        virtual json to_json() const = 0;
        virtual EnvBase* copy_ptr() const = 0;
        virtual _cd eval_t(double arg) const = 0;
        virtual _cd eval_w(double arg) const = 0;
        virtual void peak_data(PeakData& pdata) const = 0;
    };


    // Represents a Gaussian Pulse without chirp
    class EnvGaussian : public EnvBase {
    public:
        double ampl;
        double sigma;

        ~EnvGaussian();
        EnvGaussian();
        EnvGaussian(const EnvGaussian&);
        EnvGaussian(EnvGaussian&&);
        EnvGaussian& operator= (const EnvGaussian&);
        EnvGaussian& operator= (EnvGaussian&&);
        EnvGaussian(double _ampl, double _sigma);

        virtual json to_json() const override;
        virtual EnvBase* copy_ptr() const override;
        virtual _cd eval_t(double arg) const override;
        virtual _cd eval_w(double arg) const override;
        virtual void peak_data(PeakData& pdata) const override;
    };


    // Represents a chirped Gaussian Pulse
    class EnvChirped : public EnvBase {
    public:
        double ampl;
        double sigma;
        double chirp;

        ~EnvChirped();
        EnvChirped();
        EnvChirped(const EnvChirped&);
        EnvChirped(EnvChirped&&);
        EnvChirped& operator= (const EnvChirped&);
        EnvChirped& operator= (EnvChirped&&);
        EnvChirped(double _ampl, double _sigma, double _chirp);

        virtual json to_json() const override;
        virtual EnvBase* copy_ptr() const override;
        virtual _cd eval_t(double arg) const override;
        virtual _cd eval_w(double arg) const override;
        virtual void peak_data(PeakData& pdata) const override;
    };


    typedef std::unique_ptr<EnvBase> Env;

}
