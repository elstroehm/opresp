#include "interaction.hpp"

using namespace QSim;



Interaction::~Interaction() = default;

Interaction::Interaction() = default;

Interaction::Interaction(const Interaction& other) : ncontacts(other.ncontacts), couplings(other.couplings), envelopes(other.envelopes.size()),
        w0(other.w0), nphase_matches(other.nphase_matches), phase_match_signatures(other.phase_match_signatures), w_sum{other.w_sum} {
    for (int i = 0; i < other.envelopes.size(); i++) {
        envelopes[i].reset(other.envelopes[i]->copy_ptr());
    }
}

Interaction::Interaction(Interaction&& other) : ncontacts(other.ncontacts), couplings(other.couplings), envelopes{}, w0{},
        nphase_matches(other.nphase_matches), phase_match_signatures{}, w_sum{} {
    std::swap(envelopes, other.envelopes);
    std::swap(w0, other.w0);
    std::swap(phase_match_signatures, other.phase_match_signatures);
    std::swap(w_sum, other.w_sum);
}

Interaction& Interaction::operator= (const Interaction& other) {
    ncontacts = other.ncontacts;
    couplings = other.couplings;
    envelopes = std::vector<Env>(other.envelopes.size());
    nphase_matches = other.nphase_matches;
    phase_match_signatures = other.phase_match_signatures;
    w_sum = other.w_sum;
    for (int i = 0; i < other.envelopes.size(); i++) {
        envelopes[i].reset(other.envelopes[i]->copy_ptr());
    }
    return *this;
}

Interaction& Interaction::operator= (Interaction&& other) {
    ncontacts = other.ncontacts;
    couplings = other.couplings;
    std::swap(envelopes, other.envelopes);
    nphase_matches = other.nphase_matches;
    std::swap(phase_match_signatures, other.phase_match_signatures);
    std::swap(w_sum, other.w_sum);
    return *this;
}

Interaction::Interaction(int _ncontacts, const SparseMatrix<_cd>& _couplings, const std::vector<Env>& _envelopes, const _dvec& _w0, 
            int _nphase_matches, const _ivec& _phase_match_signatures)
    : ncontacts(_ncontacts), couplings(_couplings), envelopes(_envelopes.size()), w0(_w0),
        nphase_matches(_nphase_matches), phase_match_signatures{_phase_match_signatures}, w_sum{} {
    for (int i = 0; i < envelopes.size(); i++) {
        envelopes[i].reset(_envelopes[i]->copy_ptr());
    }
    for (int i = 0; i < nphase_matches; i++) {
        double _w_sum = 0.0;
        for (int j = 0; j < ncontacts; j++) {
            if (phase_match_signatures[i*ncontacts + j] == 1) {
                _w_sum += w0[j];
            }
            else {
                _w_sum -= w0[j];
            }
            w_sum.push_back(_w_sum);
        }
    }
}


// void Interaction::compute_phase_matching() {
//     TupGen tup_gen = TupGen(_ivec(ncontacts, 1));
//     for (tup_gen.init_initial(); !tup_gen.is_final(); ++tup_gen) {
//         VectorXd dk = ksig;
//         for (int i = 0; i < ncontacts; i++) {
//             if (tup_gen[i] == 0) {
//                 dk -= k0[i];
//             }
//             else {
//                 dk += k0[i];
//             }
//         }
//         if (dk.norm() < max_k_mismatch) {
//             ++nphase_matches;
//             double _w_sum = 0.0;
//             for (int i = 0; i < ncontacts; i++) {
//                 if (tup_gen.state[i] == 0) {
//                     phase_match_signatures.push_back(1);
//                     _w_sum += w0[i];
//                 }
//                 else {
//                     phase_match_signatures.push_back(-1);
//                     _w_sum -= w0[i];
//                 }
//                 w_sum.push_back(_w_sum);
//             }
//         }
//         dk = -ksig;
//         for (int i = 0; i < ncontacts; i++) {
//             if (tup_gen[i] == 0) {
//                 dk -= k0[i];
//             }
//             else {
//                 dk += k0[i];
//             }
//         }
//         if (dk.norm() < max_k_mismatch) {
//             ++nphase_matches;
//             double _w_sum = 0.0;
//             for (int i = 0; i < ncontacts; i++) {
//                 if (tup_gen.state[i] == 0) {
//                     phase_match_signatures.push_back(1);
//                     _w_sum += w0[i];
//                 }
//                 else {
//                     phase_match_signatures.push_back(-1);
//                     _w_sum -= w0[i];
//                 }
//                 w_sum.push_back(_w_sum);
//             }
//         }
//     }
// }


SparseMatrix<_cd> Interaction::eval_t(double arg, int signature, int contact) {
    const int sign = phase_match_signatures[signature*ncontacts + contact];
    const _cd value = std::exp(-Imd*w0[contact]*arg)*envelopes[contact]->eval_t(arg);
    if (sign == 1) {
        return value*couplings;
    }
    else {
        return std::conj(value)*couplings;
    }
}


SparseMatrix<_cd> Interaction::eval_w(double arg, int signature, int contact) {
    const int sign = phase_match_signatures[signature*ncontacts + contact];
    if (sign == 1) {
        return envelopes[contact]->eval_w(arg-w0[contact]) * couplings;
    }
    else {
        return std::conj(envelopes[contact]->eval_w(-arg-w0[contact])) * couplings;
    }
}


PeakData Interaction::peak_data() {
    PeakData pdata = PeakData();
    for (int i = 0; i < ncontacts; i++) {
        envelopes[i]->peak_data(pdata);
    }
    return pdata;
}


void QSim::to_json(json& ser, const Interaction& val) {
    json envelopes_ser = json::array();
    for (int i = 0; i < val.envelopes.size(); i++) {
        envelopes_ser.push_back(val.envelopes[i]->to_json());
    }
    ser = json {
        {"ncontacts", val.ncontacts},
        {"couplings", val.couplings},
        {"envelopes", envelopes_ser},
        {"w0", val.w0},
        {"phase_matches", val.nphase_matches},
        {"phase_match_signatures", val.phase_match_signatures},
    };
}


void QSim::from_json(const json& ser, Interaction& val) {
    int ncontacts = ser["ncontacts"].get<int>();
    std::vector<Env> envelopes = {};
    for (int i = 0; i < ncontacts; i++) {
        envelopes.push_back( Env(EnvBase::from_json(ser["envelopes"][i])) );
    }
    val = Interaction {
        ncontacts,
        ser["couplings"].get<SparseMatrix<_cd>>(),
        envelopes,
        ser["w0"].get<_dvec>(),
        ser["phase_matches"].get<int>(),
        ser["phase_match_signatures"].get<_ivec>(),
    };
}