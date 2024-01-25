#pragma once

#include <complex>
#include <vector>


#define JSON_OBJECT_DERIALIZATION(...) 


namespace QSim {

    typedef unsigned int uint;
    typedef unsigned long ulong;

    typedef std::vector<int> _ivec;
    typedef std::complex<double> _cd;
    typedef std::complex<float> _cf;
    typedef std::vector<std::complex<double>> _cdvec;
    typedef std::vector<std::complex<float>> _cfvec;
    typedef std::vector<double> _dvec;
    typedef std::vector<float> _fvec;

    typedef std::complex<int> _ci;

    constexpr double h = 6.62607015e-13;
    constexpr double hbar = 1.054571817e-34;
    constexpr double hbarinv = 9.482521562467288e+33;
    constexpr double e = 1.602176634e-19;
    constexpr double kB = 1.380649e-23;
    constexpr double rtemp = 293.0;
    constexpr double pi = 3.14159265359;
    constexpr double twopi = 6.28318530718;
    constexpr double twopiinv = 0.1591549430919;
    constexpr std::complex<double> Imd = _cd(0.0, 1.0);
    constexpr std::complex<float> Imf = _cf(0.0, 1.0);
    constexpr std::complex<long> Iml = std::complex<long>(0, 1);
    constexpr std::complex<int> Imi = std::complex<int>(0, 1);

}