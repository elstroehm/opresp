#include "eigen_serialization.hpp"

using namespace QSim;



json QSim::sparse_mat_to_json_dcomplex(const SparseMatrix<_cd>& mat) {
    json entries = json::array();
    for (int i = 0; i < mat.outerSize(); i++) {
        for (SparseMatrix<_cd>::InnerIterator it(mat, i); it; ++it) {
            entries.push_back( json {
                {"row", it.row()},
                {"col", it.col()},
                {"entry_real", it.value().real()},
                {"entry_imag", it.value().imag()}
            });
        }
    }
    return json {
        {"rows", mat.rows()},
        {"cols", mat.cols()},
        {"entries", entries}
    };
}


SparseMatrix<_cd> QSim::sparse_mat_from_json_dcomplex(const json& item) {
    const int rows = item["rows"].get<int>();
    const int cols = item["cols"].get<int>();
    int nentries = item["entries"].size();
    SparseMatrix<_cd> mat = SparseMatrix<_cd>(rows, cols);
    for (int i = 0; i < nentries; i++) {
        const int row = item["entries"][i]["row"].get<int>();
        const int col = item["entries"][i]["col"].get<int>();
        const _cd entry = {item["entries"][i]["entry_real"].get<double>(), item["entries"][i]["entry_imag"].get<double>()};
        mat.insert(row, col) = entry;
    }
    return mat;
}


json QSim::sparse_vec_to_json_dcomplex(const SparseVector<_cd>& vec) {
    json entries = json::array();
    for (SparseVector<_cd>::InnerIterator it(vec); it; ++it) {
        entries.push_back( json {
            {"row", it.row()},
            {"entry_real", it.value().real()},
            {"entry_imag", it.value().imag()}
        });
    }
    return json {
        {"rows", vec.rows()},
        {"entries", entries}
    };
}


SparseVector<_cd> QSim::sparse_vec_from_json_dcomplex(const json& item) {
    const int rows = item["rows"].get<int>();
    int nentries = item["entries"].size();
    SparseVector<_cd> vec = SparseVector<_cd>(rows);
    for (int i = 0; i < nentries; i++) {
        const int row = item["entries"][i]["row"].get<int>();
        const _cd entry = {item["entries"][i]["entry_real"].get<double>(), item["entries"][i]["entry_imag"].get<double>()};
        vec.insert(row, 0) = entry;
    }
    return vec;
}


json QSim::dense_mat_to_json_dcomplex(const MatrixXcd& mat) {
    json entries_real = json::array();
    json entries_imag = json::array();
    for (int i = 0; i < mat.rows(); i++) {
        entries_real.push_back(json::array());
        entries_imag.push_back(json::array());
        for (int j = 0; j < mat.cols(); j++) {
            entries_real[i].push_back(mat(i,j).real());
            entries_imag[i].push_back(mat(i,j).imag());
        }
    }
    return json {
        {"rows", mat.rows()},
        {"cols", mat.cols()},
        {"entries_real", entries_real},
        {"entries_imag", entries_imag},
    };
}


MatrixXcd QSim::dense_mat_from_json_dcomplex(const json& item) {
    const int rows = item["rows"].get<int>();
    const int cols = item["cols"].get<int>();
    MatrixXcd mat = MatrixXcd(rows, cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < mat.cols(); j++) {
            mat(i, j) = {item["entries_real"][i][j].get<double>(), item["entries_imag"][i][j].get<double>()};
        }
    }
    return mat;
}


json QSim::dense_vec_to_json_dcomplex(const VectorXcd& vec) {
    json entries_real = json::array();
    json entries_imag = json::array();
    for (int i = 0; i < vec.rows(); i++) {
        entries_real.push_back(vec(i).real());
        entries_imag.push_back(vec(i).imag());
    }
    return json {
        {"rows", vec.rows()},
        {"entries_real", entries_real},
        {"entries_imag", entries_imag},
    };
}


VectorXcd QSim::dense_vec_from_json_dcomplex(const json& item) {
    const int rows = item["rows"].get<int>();
    VectorXcd vec = VectorXcd(rows);
    for (int i = 0; i < rows; i++) {
        vec(i) = {item["entries_real"][i].get<double>(), item["entries_imag"][i].get<double>()};
    }
    return vec;
}





json QSim::sparse_mat_to_json_double(const SparseMatrix<double>& mat) {
    json entries = json::array();
    for (int i = 0; i < mat.outerSize(); i++) {
        for (SparseMatrix<double>::InnerIterator it(mat, i); it; ++it) {
            entries.push_back( json {
                {"row", it.row()},
                {"col", it.col()},
                {"entry", it.value()},
            });
        }
    }
    return json {
        {"rows", mat.rows()},
        {"cols", mat.cols()},
        {"entries", entries}
    };
}


SparseMatrix<double> QSim::sparse_mat_from_json_double(const json& item) {
    const int rows = item["rows"].get<int>();
    const int cols = item["cols"].get<int>();
    int nentries = item["entries"].size();
    SparseMatrix<double> mat = SparseMatrix<double>(rows, cols);
    for (int i = 0; i < nentries; i++) {
        const int row = item["entries"][i]["row"].get<int>();
        const int col = item["entries"][i]["col"].get<int>();
        const double entry = item["entries"][i]["entry"].get<double>();
        mat.insert(row, col) = entry;
    }
    return mat;
}


json QSim::sparse_vec_to_json_double(const SparseVector<double>& vec) {
    json entries = json::array();
    for (SparseVector<double>::InnerIterator it(vec); it; ++it) {
        entries.push_back( json {
            {"row", it.row()},
            {"entry", it.value()},
        });
    }
    return json {
        {"rows", vec.rows()},
        {"entries", entries}
    };
}


SparseVector<double> QSim::sparse_vec_from_json_double(const json& item) {
    const int rows = item["rows"].get<int>();
    int nentries = item["entries"].size();
    SparseVector<double> vec = SparseVector<double>(rows);
    for (int i = 0; i < nentries; i++) {
        const int row = item["entries"][i]["row"].get<int>();
        const double entry = item["entries"][i]["entry"].get<double>();
        vec.insert(row, 0) = entry;
    }
    return vec;
}


json QSim::dense_mat_to_json_double(const MatrixXd& mat) {
    json entries = json::array();
    for (int i = 0; i < mat.rows(); i++) {
        entries.push_back(json::array());
        for (int j = 0; j < mat.cols(); j++) {
            entries[i].push_back(mat(i,j));
        }
    }
    return json {
        {"rows", mat.rows()},
        {"cols", mat.cols()},
        {"entries", entries},
    };
}


MatrixXd QSim::dense_mat_from_json_double(const json& item) {
    const int rows = item["rows"].get<int>();
    const int cols = item["cols"].get<int>();
    MatrixXd mat = MatrixXd(rows, cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < mat.cols(); j++) {
            mat(i, j) = item["entries"][i][j].get<double>();
        }
    }
    return mat;
}


json QSim::dense_vec_to_json_double(const VectorXd& vec) {
    json entries = json::array();
    for (int i = 0; i < vec.rows(); i++) {
        entries.push_back(vec(i));
    }
    return json {
        {"rows", vec.rows()},
        {"entries", entries},
    };
}


VectorXd QSim::dense_vec_from_json_double(const json& item) {
    const int rows = item["rows"].get<int>();
    VectorXd vec = VectorXd(rows);
    for (int i = 0; i < rows; i++) {
        vec(i) = item["entries"][i].get<double>();
    }
    return vec;
}




