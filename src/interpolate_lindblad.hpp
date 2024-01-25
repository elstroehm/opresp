// #pragma once

// #include <complex>
// #include <vector>
// #include <iostream>
// #include <memory>
// #include <algorithm>
// #include <thread>

// #include "../external/eigen-3.4.0/Eigen/Core"
// #include "../external/eigen-3.4.0/unsupported/Eigen/MatrixFunctions"
// #include "../external/eigen-3.4.0/unsupported/Eigen/FFT"
// #include "../external/json/json.hpp"
// #include "../external/math/interpolators/cardinal_cubic_b_spline.hpp"

// #include "definitions.hpp"
// #include "functions.hpp"

// using namespace Eigen;
// using json = nlohmann::json;
// using spline = boost::math::interpolators::cardinal_cubic_b_spline<double>;



// namespace QSim {


//     // int nintervals;
//     // double tstart;
//     // double tstop;
//     // double wstart;
//     // double wstop;
//     // spline t_evolve_real;
//     // spline t_evolve_imag;
//     // spline w_evolve_real;
//     // spline w_evolve_imag;


//     std::vector<MatrixXcd> const_ode_solve(MatrixXcd A, double t0, double t1, int nintervals, int nthreads);
//     void const_ode_solve_st(MatrixXcd A, double tstart, double tstep, int nsteps, std::vector<MatrixXcd>& results, int start_index);


//     Fn interpolate_and_fft(const std::vector<MatrixXcd>& X, int row, int col, double tstart, double tstop) {
//         const int npoints = X.size();
//         const int nintervals = npoints - 1;
//         const double tstep = (tstop - tstart)/nintervals;
//         _dvec t_real_data = _dvec(npoints);
//         _dvec t_imag_data = _dvec(npoints);
//         _cdvec t_data = _cdvec(npoints);
//         for (int i = 0; i < npoints; i++) {
//             t_real_data[i] = X[i](row,col).real();
//             t_imag_data[i] = X[i](row,col).imag();
//             t_data[i] = X[i](row,col);
//         }
//         _cdvec w_data = _cdvec(npoints);
//         FFT<_cd> fft;
//         fft.fwd(w_data, t_data);

//         spline t_evolve_real = spline(t_real_data.data(), npoints, tstart, tstep);
//         spline t_evolve_imag = spline(t_imag_data.data(), npoints, tstart, tstep);
//         return Fn { new FnInterpolation {
//             nintervals,
//             tstart,
//             tstop,
//             wstart,
//             wstop,
//             t_evolve_real,
//             t_evolve_imag,
//             w_evolve_real,
//             w_evolve_imag,
//         }};
//     }


// }