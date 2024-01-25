// #include "interpolate_lindblad.hpp"

// using namespace QSim;



// std::vector<MatrixXcd> QSim::const_ode_solve(MatrixXcd A, double t0, double t1, int nintervals, int nthreads) {
//     const int npoints = nintervals + 1;
//     const int ppt2 = npoints / nthreads;
//     const int nt1 = npoints - ppt2*nthreads;
//     const int nt2 = nthreads - nt1;
//     const int ppt1 = ppt2 + 1;
//     const double tstep = (t1 - t0)/nintervals;
//     std::vector<MatrixXcd> results = std::vector<MatrixXcd>(npoints);
//     std::vector<std::thread> workbench = std::vector<std::thread>();
//     int start_index = 0;
//     double tstart = t0;
//     for (int i = 0; i < nt1; i++) {
//         workbench.push_back(std::thread(const_ode_solve_st, A, tstart, tstep, ppt1, results, start_index));
//         start_index += ppt1;
//         tstart += static_cast<double>(ppt1) * tstep;
//     }
//     for (int i = 0; i < nt2; i++) {
//         workbench.push_back(std::thread(const_ode_solve_st, A, tstart, tstep, ppt2, results, start_index));
//         start_index += ppt2;
//         tstart += static_cast<double>(ppt2) * tstep;
//     }
//     for (int i = 0; i < nthreads; i++) {
//         workbench[i].join();
//     }
//     return results;
// }


// void QSim::const_ode_solve_st(MatrixXcd A, double tstart, double tstep, int nsteps, std::vector<MatrixXcd>& results, int start_index) {
//     for (int i = 0; i < nsteps; i++) {
//         results[start_index + i] = ((tstart + i*tstep)*A).eval().exp();
//     }
// }