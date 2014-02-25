// Visit http://www.johndcook.com/stand_alone_code.html for the source of this code and more like it.

#ifndef _gamma_h_
#define _gamma_h_

#include <cmath>
#include <sstream>
#include <iostream>
#include <stdexcept>

// Note that the functions Gamma and LogGamma are mutually dependent.
double LogGamma(double);
double Gamma(double);

#if defined(_WIN32) || defined(WIN32)
  inline int pow(int x, int y) {
    return (int) pow((double) x, y);
  }
  inline double tgamma(double x) { return Gamma(x); }
  inline double lgamma(double x) { return LogGamma(x); }
#endif //_WIN32

#endif //_gamma_h_
