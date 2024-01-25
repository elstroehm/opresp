#pragma once

#include <complex>
#include <vector>
#include <iostream>
#include <functional>
#include <memory>

#include "../external/eigen-3.4.0/Eigen/Core"
#include "../external/eigen-3.4.0/Eigen/Sparse"
#include "../external/json/json.hpp"

#include "definitions.hpp"
#include "functions.hpp"
#include "eigen_serialization.hpp"
#include "algebra.hpp"

using namespace Eigen;
using json = nlohmann::json;



namespace QSim {    

    class Interaction {
    public:
        int ncontacts;
        SparseMatrix<_cd> couplings;
        std::vector<Env> envelopes;     // Entspricht den Funktionen a(t) in den Notizen
        std::vector<double> w0;
        int nphase_matches;
        _ivec phase_match_signatures;
        _dvec w_sum;

        ~Interaction();
        Interaction();
        Interaction(const Interaction&);
        Interaction(Interaction&&);
        Interaction& operator= (const Interaction&);
        Interaction& operator= (Interaction&&);
        Interaction(int _ncontacts, const SparseMatrix<_cd>& _couplings, const std::vector<Env>& _envelopes, const _dvec& _w0, 
            int _nphase_matches, const _ivec& _phase_match_signatures);
        
        SparseMatrix<_cd> eval_t(double arg, int signature, int contact);
        SparseMatrix<_cd> eval_w(double arg, int signature, int contact);
        PeakData peak_data();
    };

    void to_json(json& ser, const Interaction& val);
    void from_json(const json& ser, Interaction& val);

}
