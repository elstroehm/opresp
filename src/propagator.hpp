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

using namespace Eigen;
using json = nlohmann::json;



namespace QSim {    



    class Propagator {
    public:
        int rows;
        int cols;
        std::vector<int> indices;
        std::vector<int> starts;
        std::vector<int> counts;
        std::vector<Fn> fns;

        ~Propagator();
        Propagator();
        Propagator(const Propagator&);
        Propagator(Propagator&&);
        Propagator& operator= (const Propagator&);
        Propagator& operator= (Propagator&&);

        Propagator(int rows, int cols);
        static Propagator init_free(const _dvec& omegas, double eta);
        static Propagator init_naive(const MatrixXd& dephasings, const _dvec& omegas, double eta);
        static Propagator init_simple(_dvec decay_rates, const MatrixXd& dephasings, const _dvec& omegas, double eta);
        static Propagator init_lindblad_diagonal(const VectorXcd& eigenvals, double eta);

        static Propagator from_json(const json& item);
        json to_json() const;

        void insert(int i, int j, const Fn& fn, bool overwrite = true);
        void reserve(const std::vector<int>& entries);

        PeakData peak_data() const;

        SparseMatrix<_cd> eval_t(double t) const;
        SparseMatrix<_cd> eval_w(double w) const;
    };
    
    void to_json(json& ser, const Propagator& mat);
    void from_json(const json& ser, Propagator& mat);

}

