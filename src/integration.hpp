#pragma once

#include <complex>
#include <vector>
#include <iostream>
#include <memory>
#include <algorithm>
#include <functional>

#include "../external/eigen-3.4.0/Eigen/Core"
#include "../external/Cuba-4.2.2/cuba.h"
#include "../external/json/json.hpp"

#include "definitions.hpp"
#include "algebra.hpp"

using namespace Eigen;
using json = nlohmann::json;



namespace QSim {


    // Schlüssel für die verschiedenen CUBA Algorithmen
    enum class CubaAlg : int {
        Vegas,
        Suave,
        Divonne,
        Cuhre
    };

    // Manuell festlegbare Parameter, die bei allen CUBA-Routinen gleich sind
    struct CubaCommonParams {
        double epsrel;
        double epsabs;
        int flags;
        int mineval;
        int maxeval;
    };

    // Alle Parameter, welche den Integranden auszeichnen
    struct CubaIntegrand {
        integrand_t integrand;
        int ndim;
        int ncomp;
        int nvec;
    };

    // Vegas-exklusive Parameter
    struct CubaParamsVegas {
        int seed;
        int nstart;
        int nincrease;
        int nbatch;
        int gridno;

        static CubaParamsVegas with_default_args();
    };

    // Suave-exklusive Parameter
    struct CubaParamsSuave {
        int seed;
        int nnew;
        int nmin;
        double flatness;

        static CubaParamsSuave with_default_args();
    };

    // Divonne-exklusive Parameter
    struct CubaParamsDivonne {
        int seed;
        int key1;
        int key2;
        int key3;
        int maxpass;
        double border;
        double maxchisq;
        double mindeviation;
        int ngiven;
        int ldxgiven;
        _dvec xgiven;
        int nextra;
        peakfinder_t peakfinder;

        static CubaParamsDivonne with_default_args();
    };

    // Cuhre-exklusive Parameter
    struct CubaParamsCuhre {
        int key;

        static CubaParamsCuhre with_default_args();
    };

    // Zusammenfassung aller Ausgabegrößen der CUBA-Routinen, für alle gleich (Vegas etwas speziell)
    struct CubaResult {
        int ncomp;
        int nregions;
        int neval;
        int fail;
        _dvec integral;
        _dvec prob;
        _dvec error;
    };



    NLOHMANN_JSON_SERIALIZE_ENUM( CubaAlg, {
        {CubaAlg::Vegas, "Vegas"},
        {CubaAlg::Suave, "Suave"},
        {CubaAlg::Divonne, "Divonne"},
        {CubaAlg::Cuhre, "Cuhre"}
    })

    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(CubaCommonParams, epsrel, epsabs, flags, mineval, maxeval)
    // NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(CubaParamsVegas, seed, nstart, nincrease, nbatch, gridno)
    // NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(CubaParamsSuave, seed, nnew, nmin, flatness)
    // NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(CubaParamsDivonne, seed, key1, key2, key3, maxpass, border, maxchisq, mindeviation, ngiven, ldxgiven, xgiven, nextra)
    // NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(CubaParamsCuhre, key)
    // NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(CubaResult, nregions, ncomp, neval, fail, integral, prob, error)



    // Wrapper für die einzelnen CUBA-Routinen
    
    CubaResult cuba_vegas_integrate(CubaIntegrand integrand, void* userdata,
        const CubaParamsVegas& params, const CubaCommonParams& cuba_args);

    CubaResult cuba_suave_integrate(CubaIntegrand integrand, void* userdata,
        const CubaParamsSuave& params, const CubaCommonParams& cuba_args);

    CubaResult cuba_divonne_integrate(CubaIntegrand integrand, void* userdata,
        const CubaParamsDivonne& params, const CubaCommonParams& cuba_args);

    CubaResult cuba_cuhre_integrate(CubaIntegrand integrand, void* userdata,
        const CubaParamsCuhre& params, const CubaCommonParams& cuba_args);
        

    void cuba_vegas_integrate_inplace(CubaIntegrand integrand, void* userdata,
        CubaParamsVegas params, CubaCommonParams cuba_args, CubaResult* result);

    void cuba_suave_integrate_inplace(CubaIntegrand integrand, void* userdata,
        CubaParamsSuave params, CubaCommonParams cuba_args, CubaResult* result);

    void cuba_divonne_integrate_inplace(CubaIntegrand integrand, void* userdata,
        CubaParamsDivonne params, CubaCommonParams cuba_args, CubaResult* result);

    void cuba_cuhre_integrate_inplace(CubaIntegrand integrand, void* userdata,
        CubaParamsCuhre params, CubaCommonParams cuba_args, CubaResult* result);

    

    // Wrapper für beliebige CUBA-Routine mit bereits definierten Parametern
    typedef std::function<void(void*, CubaCommonParams, CubaResult*)> ExplicitIntegrator;
    typedef std::function<void(void*, CubaResult*)> CubaIntegrator;

   

    
    // #define SOBOL_MINDIM 1
    // #define SOBOL_MAXDIM 40
    // #define KOROBOV_MINDIM 2
    // #define KOROBOV_MAXDIM 33
    // #define MAXDIM 1024

    // #define IsSobol(k) NegQ(k)
    // #define NegQ(a) ((a) >> (sizeof(a)*8 - 1))
    // #define IsRule(k, d) (k == 9 || k == 7 || (k == 11 && d == 3) || (k == 13 && d == 2))

    // inline bool BadDimension(int ndim, int seed, int key)
    // {
    //     if( ndim > MAXDIM ) return true;
    //     if( IsSobol(key) ) return
    //         ndim < SOBOL_MINDIM || (seed == 0 && ndim > SOBOL_MAXDIM);
    //     if( IsRule(key, ndim) ) return ndim < 1;
    //     return ndim < KOROBOV_MINDIM || ndim > KOROBOV_MAXDIM;
    // }

    // inline int BadDimensionBranchTaken(int ndim, int seed, int key)
    // {
    //     if( ndim > MAXDIM ) return 1;
    //     if( IsSobol(key) ) return 2;
    //     if( IsRule(key, ndim) ) return 3;
    //     return 4;
    // }

    // std::cout << "Fehlercode:   " << result.fail << std::endl;
    // if( BadDimension(1, cuba_params.seed, cuba_params.key1) ) {
    //     std::cout << "Key 1, Zweig " << BadDimensionBranchTaken(1, cuba_params.seed, cuba_params.key1) << std::endl;
    // }
    // if (BadDimension(1, cuba_params.seed, cuba_params.key2) ) {
    //     std::cout << "Key 2, Zweig " << BadDimensionBranchTaken(1, cuba_params.seed, cuba_params.key1) << std::endl;
    // }
    // if ((cuba_params.key3 & -2) && BadDimension(1, cuba_params.seed, cuba_params.key3)) {
    //     std::cout << "Key 3, Zweig " << BadDimensionBranchTaken(1, cuba_params.seed, cuba_params.key1) << std::endl;
    // }

}
