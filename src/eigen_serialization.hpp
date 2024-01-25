#pragma once

#include <complex>
#include <vector>
#include <iostream>

#include "../external/eigen-3.4.0/Eigen/Core"
#include "../external/eigen-3.4.0/Eigen/Sparse"
#include "../external/json/json.hpp"

#include "definitions.hpp"

using namespace Eigen;
using json = nlohmann::json;



namespace QSim {    
    json sparse_mat_to_json_dcomplex(const SparseMatrix<_cd>& mat);
    SparseMatrix<_cd> sparse_mat_from_json_dcomplex(const json& item);

    json sparse_vec_to_json_dcomplex(const SparseVector<_cd>& vec);
    SparseVector<_cd> sparse_vec_from_json_dcomplex(const json& item);

    json dense_mat_to_json_dcomplex(const MatrixXcd& mat);
    MatrixXcd dense_mat_from_json_dcomplex(const json& item);

    json dense_vec_to_json_dcomplex(const VectorXcd& vec);
    VectorXcd dense_vec_from_json_dcomplex(const json& item);

    json sparse_mat_to_json_double(const SparseMatrix<double>& mat);
    SparseMatrix<double> sparse_mat_from_json_double(const json& item);

    json sparse_vec_to_json_double(const SparseVector<double>& vec);
    SparseVector<double> sparse_vec_from_json_double(const json& item);

    json dense_mat_to_json_double(const MatrixXd& mat);
    MatrixXd dense_mat_from_json_double(const json& item);

    json dense_vec_to_json_double(const VectorXd& vec);
    VectorXd dense_vec_from_json_double(const json& item);
}

namespace nlohmann {
    template <>
    struct adl_serializer<Eigen::SparseMatrix<QSim::_cd>> {
        static Eigen::SparseMatrix<QSim::_cd> from_json(const json& j) {
            return QSim::sparse_mat_from_json_dcomplex(j);
        }

        static void to_json(json& j, const Eigen::SparseMatrix<QSim::_cd>& mat) {
            j = QSim::sparse_mat_to_json_dcomplex(mat);
        }
    };

    template <>
    struct adl_serializer<Eigen::SparseVector<QSim::_cd>> {
        static Eigen::SparseVector<QSim::_cd> from_json(const json& j) {
            return QSim::sparse_vec_from_json_dcomplex(j);
        }

        static void to_json(json& j, const Eigen::SparseVector<QSim::_cd>& vec) {
            j = QSim::sparse_vec_to_json_dcomplex(vec);
        }
    };

    template <>
    struct adl_serializer<Eigen::MatrixXcd> {
        static Eigen::MatrixXcd from_json(const json& j) {
            return QSim::dense_mat_from_json_dcomplex(j);
        }

        static void to_json(json& j, const Eigen::MatrixXcd& mat) {
            j = QSim::dense_mat_to_json_dcomplex(mat);
        }
    };

    template <>
    struct adl_serializer<Eigen::VectorXcd> {
        static Eigen::VectorXcd from_json(const json& j) {
            return QSim::dense_vec_from_json_dcomplex(j);
        }

        static void to_json(json& j, const Eigen::VectorXcd& vec) {
            j = QSim::dense_vec_to_json_dcomplex(vec);
        }
    };



    template <>
    struct adl_serializer<Eigen::SparseMatrix<double>> {
        static Eigen::SparseMatrix<double> from_json(const json& j) {
            return QSim::sparse_mat_from_json_double(j);
        }

        static void to_json(json& j, const Eigen::SparseMatrix<double>& mat) {
            j = QSim::sparse_mat_to_json_double(mat);
        }
    };

    template <>
    struct adl_serializer<Eigen::SparseVector<double>> {
        static Eigen::SparseVector<double> from_json(const json& j) {
            return QSim::sparse_vec_from_json_double(j);
        }

        static void to_json(json& j, const Eigen::SparseVector<double>& vec) {
            j = QSim::sparse_vec_to_json_double(vec);
        }
    };

    template <>
    struct adl_serializer<Eigen::MatrixXd> {
        static Eigen::MatrixXd from_json(const json& j) {
            return QSim::dense_mat_from_json_double(j);
        }

        static void to_json(json& j, const Eigen::MatrixXd& mat) {
            j = QSim::dense_mat_to_json_double(mat);
        }
    };

    template <>
    struct adl_serializer<Eigen::VectorXd> {
        static Eigen::VectorXd from_json(const json& j) {
            return QSim::dense_vec_from_json_double(j);
        }

        static void to_json(json& j, const Eigen::VectorXd& vec) {
            j = QSim::dense_vec_to_json_double(vec);
        }
    };
}
