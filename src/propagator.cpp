#include "propagator.hpp"

using namespace QSim;



Propagator::~Propagator() = default;

Propagator::Propagator() = default;

Propagator::Propagator(const Propagator& other) : rows(other.rows), cols(other.cols), indices(other.indices),
    starts(other.starts), counts(other.counts), fns(other.fns.size()) {
    for (int i = 0; i < cols; i++) {
        for (int j = starts[i]; j < starts[i]+counts[i]; j++) {
            fns[j].reset(other.fns[j]->copy_ptr());
        }
    }
}

Propagator::Propagator(Propagator&& other) : rows(other.rows), cols(other.cols), indices{}, starts{}, counts{}, fns{} {
    std::swap(indices, other.indices);
    std::swap(starts, other.starts);
    std::swap(counts, other.counts);
    std::swap(fns, other.fns);
}

Propagator& Propagator::operator= (const Propagator& other) {
    if (this == &other) {
        return *this;
    }
    else {
        rows = other.rows;
        cols = other.cols;
        indices.resize(other.indices.size());
        starts.resize(other.starts.size());
        counts.resize(other.counts.size());
        fns.resize(other.fns.size());
        for (int i = 0; i < cols; i++) {
            starts[i] = other.starts[i];
            counts[i] = other.counts[i];
        }
        starts[cols] = other.starts[cols];
        for (int i = 0; i < other.indices.size(); i++) {
            indices[i] = other.indices[i];
            fns[i].reset(other.fns[i]->copy_ptr());
        }
        return *this;
    }
}

Propagator& Propagator::operator= (Propagator&& other) {
    if (this == &other) {
        return *this;
    }
    else {
        rows = other.rows;
        cols = other.cols;
        std::swap(indices, other.indices);
        std::swap(starts, other.starts);
        std::swap(counts, other.counts);
        std::swap(fns, other.fns);
        return *this;
    }
}

Propagator::Propagator(int _rows, int _cols) : rows(_rows), cols(_cols), indices{}, starts(_ivec(cols+1, 0)),
    counts(_ivec(cols, 0)), fns{} {
    int pow = 1;
    while (pow*pow < rows) {
        pow *= 2;
    }
    if (pow > 1) {pow /= 2;}
    indices = std::vector<int>(cols*pow, -1);
    fns = std::vector<Fn>(cols*pow);
    for (int i = 0; i < cols*pow; i++) {
        fns.emplace_back(nullptr);
    }
    for (int i = 0; i <= cols; i++) {
        starts[i] = i*pow;
    }
}


Propagator Propagator::init_free(const _dvec& omegas, double eta) {
    int nlevels = omegas.size();
    int dim = nlevels*nlevels;
    Propagator prop = Propagator(dim, dim);
    for (int k = 0; k < nlevels; k++) {
        for (int l = 0; l < nlevels; l++) {
            Fn fn = Fn(new FnSpectrum {
                {_cd{1.0, 0.0}},
                {omegas[k] - omegas[l]},
                {0.0},
                {true},
                eta
            });
            prop.insert(k*nlevels+l, k*nlevels+l, fn);
        }
    }
    return prop;
}


Propagator Propagator::init_naive(const MatrixXd& dephasings, const _dvec& omegas, double eta) {
    int nlevels = omegas.size();
    int dim = nlevels*nlevels;
    Propagator prop = Propagator(dim, dim);
    for (int k = 0; k < nlevels; k++) {
        for (int l = 0; l < nlevels; l++) {
            if (k != l) {
                Fn fn = Fn(new FnSpectrum {
                    {_cd{1.0, 0.0}},
                    {omegas[k] - omegas[l]},
                    {dephasings(k,l)},
                    {true},
                    eta
                });
                prop.insert(k*nlevels+l, k*nlevels+l, fn);
            }
            else {
                Fn fn = Fn(new FnSpectrum {
                    {_cd{1.0, 0.0}},
                    {0.0},
                    {0.0},
                    {true},
                    eta
                });
                prop.insert(k*nlevels+l, k*nlevels+l, fn);
            }
        }
    }
    return prop;
}


Propagator Propagator::init_simple(_dvec decay_rates, const MatrixXd& dephasings, const _dvec& omegas, double eta) {
    int nlevels = omegas.size();
    int dim = nlevels*nlevels;
    Propagator prop = Propagator(dim, dim);
    for (int k = 0; k < nlevels; k++) {
        for (int l = 0; l < nlevels; l++) {
            if (k != l) {
                Fn fn = Fn(new FnSpectrum {
                    {_cd{1.0, 0.0}},
                    {omegas[k] - omegas[l]},
                    {dephasings(k,l)},
                    {true},
                    eta
                });
                prop.insert(k*nlevels+l, k*nlevels+l, fn);
            }
            if (k == 0) {
                Fn fn = Fn(new FnSpectrum {
                    {_cd{1.0, 0.0}, -_cd{1.0, 0.0}},
                    {0.0, 0.0},
                    {0.0, decay_rates[l]},
                    {true, true},
                    eta
                });
                prop.insert(k*(nlevels+1), l*(nlevels+1), fn);
            }
            else if (k == 0 && l == 0) {
                Fn fn = Fn(new FnSpectrum {
                    {_cd{1.0, 0.0}},
                    {0.0},
                    {0.0},
                    {true},
                    eta});
                prop.insert(k*(nlevels+1), l*(nlevels+1), fn);
            }
            else if (k == l) {
                Fn fn = Fn(new FnSpectrum {
                    {_cd{1.0, 0.0}},
                    {0.0},
                    {decay_rates[l]},
                    {true},
                    eta});
                prop.insert(k*(nlevels+1), l*(nlevels+1), fn);
            }
        }
    }
    return prop;
}


Propagator Propagator::init_lindblad_diagonal(const VectorXcd& eigenvals, double eta) {
    const int dim = eigenvals.rows();
    Propagator prop = Propagator(dim, dim);
    for (int k = 0; k < dim; k++) {
        Fn fn = Fn(new FnSpectrum {
            {_cd{1.0, 0.0}},
            {-eigenvals[k].imag()},
            {-eigenvals[k].real()},
            {true},
            eta
        });
        prop.insert(k, k, fn);
    }
    return prop;
}


Propagator Propagator::from_json(const json& item) {
    const int rows = item["rows"].get<int>();
    const int cols = item["cols"].get<int>();
    int nentries = item["entries"].size();
    Propagator mat = Propagator(rows, cols);
    for (int i = 0; i < nentries; i++) {
        const int row = item["entries"][i]["row"].get<int>();
        const int col = item["entries"][i]["col"].get<int>();
        const Fn fn = Fn(FnBase::from_json(item["entries"][i]["function"]));
        mat.insert(row, col, fn);
    }
    return mat;
}


json Propagator::to_json() const {
    json entries = json::array();
    for (int col = 0; col < cols; col++) {
        for (int j = starts[col]; j < starts[col]+counts[col]; j++) {
            const int row = indices[j];
            entries.push_back(json {
                {"row", row},
                {"col", col},
                {"function", fns[j]->to_json()}
            });
        }
    }
    json item = json {
        {"rows", rows},
        {"cols", cols},
        {"entries", entries}
    };
    return item;
}


void Propagator::insert(int row, int col, const Fn& fn, bool overwrite) {
    // Spezialfall, wenn alles noch leer ist
    if (indices.size() == 0) {
        indices.push_back(row);
        fns.emplace_back(fn->copy_ptr());
        for (int i = col; i <= cols; i++) {
            starts[i] = 1;
        }
        counts[col] = 1;
        return;
    }
    // Binäre Suche der Stelle, an der der Eintrag in indices und fns vorgenommen werden muss
    int lower = starts[col];
    int upper = starts[col] + counts[col];
    int step = counts[col]/2;
    while (lower != upper) {
        if (row > indices[lower+step]) {
            lower = lower + step + 1;
        }
        else {
            upper = lower + step;
        }
        step = (upper - lower) / 2;
    }
    const int index = lower;
    // Ist schon ein Eintrag an der entsprechenden Stelle vorhanden, überschreibe oder kehre zurück
    if (indices[index] == row && overwrite) {
        fns[index].reset(fn->copy_ptr());
        return;
    }
    else if (indices[index] == row && !overwrite) {
        return;
    }
    counts[col] += 1;
    int nslots = starts[col+1] - starts[col];
    // Reallokation des Buffers für die entsprechende Spalte, falls Speicher fehlt
    if (nslots < counts[col]) {
        if (nslots == 0) {nslots = 1;}
        else {
            indices.resize(indices.size() + nslots);
            fns.resize(fns.size() + nslots);
        }
        starts[cols] += nslots;
        for (int i = cols-1; i > col; i--) {
            for (int k = starts[i]+counts[i]-1; k >= starts[i]; k--) {
                indices[k+nslots] = indices[k];
                fns[k+nslots].reset(fns[k].release());
            }
            starts[i] += nslots;
        }
    }
    for (int k = starts[col]+counts[col]-1; k > index; k--) {
        indices[k] = indices[k-1];
        fns[k].reset(fns[k-1].release());
    }
    indices[index] = row;
    fns[index].reset(fn->copy_ptr());
    return;
}


void Propagator::reserve(const std::vector<int>& slots) {
    std::vector<int> shifts = std::vector<int>(cols, 0);
    for (int i = 0; i < cols; i++) {
        if (starts[i+1] - starts[i] < slots[i]) {
            shifts[i] = slots[i] - (starts[i+1] - starts[i]);
        }
    }
    for (int i = 0; i < cols-1; i++) {
        shifts[i+1] += shifts[i];
    }
    indices.resize(shifts[cols-1]);
    fns.resize(shifts[cols-1]);
    for (int i = cols-1; i >= 0; i++) {
        for (int j = starts[i]; j < starts[i]+counts[i]; j++) {
            indices[j+shifts[i]] = indices[j];
            fns[j+shifts[i]].reset(fns[j].release());
        }
        starts[i] += shifts[i];
    }
}


PeakData Propagator::peak_data() const {
    PeakData pdata = PeakData {
        {},
        {},
    };
    for (int i = 0; i < cols; i++) {
        for (int j = starts[i]; j < starts[i]+counts[i]; j++) {
            fns[j]->peak_data(pdata);
        }
    }
    return pdata;
}


SparseMatrix<_cd> Propagator::eval_t(double t) const {
    SparseMatrix<_cd> mat = SparseMatrix<_cd>(rows, cols);
    mat.reserve(counts);
    for (int col = 0; col < cols; col++) {
        for (int j = starts[col]; j < starts[col]+counts[col]; j++) {
            const int row = indices[j];
            mat.insert(row, col) = fns[j]->eval_t(t);
        }
    }
    return mat;
}


SparseMatrix<_cd> Propagator::eval_w(double w) const {
    SparseMatrix<_cd> mat = SparseMatrix<_cd>(rows, cols);
    mat.reserve(counts);
    for (int col = 0; col < cols; col++) {
        for (int j = starts[col]; j < starts[col]+counts[col]; j++) {
            const int row = indices[j];
            mat.insert(row, col) = fns[j]->eval_w(w);
        }
    }
    return mat;
}


void QSim::to_json(json& ser, const Propagator& mat) {
    json entries = json::array();
    for (int col = 0; col < mat.cols; col++) {
        for (int j = mat.starts[col]; j < mat.starts[col]+mat.counts[col]; j++) {
            const int row = mat.indices[j];
            entries.push_back(json {
                {"row", row},
                {"col", col},
                {"function", mat.fns[j]->to_json()}
            });
        }
    }
    ser = json {
        {"rows", mat.rows},
        {"cols", mat.cols},
        {"entries", entries}
    };
}


void QSim::from_json(const json& ser, Propagator& mat) {
    const int rows = ser["rows"].get<int>();
    const int cols = ser["cols"].get<int>();
    int nentries = ser["entries"].size();
    mat = Propagator(rows, cols);
    for (int i = 0; i < nentries; i++) {
        const int row = ser["entries"][i]["row"].get<int>();
        const int col = ser["entries"][i]["col"].get<int>();
        const Fn fn = Fn(FnBase::from_json(ser["entries"][i]["function"]));
        mat.insert(row, col, fn);
    }
}

