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
#include "../external/Cuba-4.2.2/cuba.h"

#include "definitions.hpp"
#include "algebra.hpp"
#include "functions.hpp"
#include "propagator.hpp"
#include "interaction.hpp"
#include "integration.hpp"
#include "eigen_serialization.hpp"

using namespace Eigen;
using json = nlohmann::json;



namespace QSim {


    // enum class Representation {
    //     DysonFreq,
    //     DysonTime,
    //     SemiImpulsiveFreq,
    //     SemiImpulsiveTime,
    // };

    // NLOHMANN_JSON_SERIALIZE_ENUM( Representation, {
    //     {Representation::DysonFreq, "DysonFreq"},
    //     {Representation::DysonTime, "DysonTime"},
    //     {Representation::SemiImpulsiveFreq, "SemiImpulsiveFreq"},
    //     {Representation::SemiImpulsiveTime, "SemiImpulsiveTime"},
    // })

    enum class Representation {
        DysonFreq,
        DysonTime,
        SemiImpulsive,
    };

    NLOHMANN_JSON_SERIALIZE_ENUM( Representation, {
        {Representation::DysonFreq, "DysonFreq"},
        {Representation::DysonTime, "DysonTime"},
        {Representation::SemiImpulsive, "SemiImpulsive"},
    })


    enum class Domain {
        Time,
        Freq,
    };

    NLOHMANN_JSON_SERIALIZE_ENUM( Domain, {
        {Domain::Time, "Time"},
        {Domain::Freq, "Freq"},
    })


    // initial state always assumed |0><0|
    class ResponseParams {
    public:
        SparseVector<_cd> observable;
        SparseMatrix<_cd> qform;
        _dvec x0;
        _dvec x1;
        Propagator propagator;
        Interaction interaction;
        _dvec args;
        std::vector<Domain> arg_domains;
        int phase_match_state;
        _dvec integration_poles;
        int npoles;
        int ndim;

        ~ResponseParams();
        ResponseParams();
        ResponseParams(const ResponseParams&);
        ResponseParams& operator= (const ResponseParams&);
        ResponseParams(ResponseParams&&);
        ResponseParams& operator= (ResponseParams&&);
        ResponseParams(const SparseVector<_cd> _observable, const SparseMatrix<_cd>& _qform, const _dvec& _x0, const _dvec& _x1, const Propagator& _propagator,
            const Interaction& _interaction, const std::vector<Domain>& _arg_domains);

        void next_phase_match();
    };

    void to_json(json& ser, const ResponseParams& val);
    void from_json(const json& ser, ResponseParams& val);


    int dyson_freq_repr(const int* _ndim, const double _args[], const int* _ncomp,
        double _val[], void* _params);

    int dyson_time_repr(const int* _ndim, const double i_args[], const int* _ncomp,
        double _val[], void* _params);

    int semi_impulsive_repr(const int* _ndim, const double i_args[], const int* _ncomp,
        double _val[], void* _params);

    // int semi_impulsive_time_repr(const int* _ndim, const double i_args[], const int* _ncomp,
    //     double _val[], void* _params);


    void peaksearch_dyson_freq(const int* ndim, const double* borders, int* npeaks, double* pos, void* userdata);
    void peaksearch_dyson_time(const int* ndim, const double* borders, int* npeaks, double* pos, void* userdata);

}
