#ifndef ALIQVECTORS_H
#define ALIQVECTORS_H

/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <array>
#include <complex>
#include <iostream>

template <unsigned int Harmonics, unsigned int Powers>
struct Qv {
  static const unsigned int kMaxH = Harmonics;
  static const unsigned int kMaxP = Powers;
  static constexpr std::size_t kSize = (kMaxH*2+1)*(kMaxP); 
  std::array<std::complex<double>, kSize> qv_;
  std::complex<double> operator()(const int ih, const int ip) const {
    return qv_[index(ih,ip)];
  }
  std::complex<double> &ref(const int ih, const int ip) {
    return qv_[index(ih,ip)];
  }
  unsigned int index(const int ih, const int ip) const {
    auto x = (ih + kMaxH) * (kMaxP) + (ip-1);
    return x;
  }
  void Fill(double phi, double weight) {
    std::array<double, Harmonics+1> cos_arr;
    std::array<double, Harmonics+1> sin_arr;
    for (auto h = 1u; h <= Harmonics; ++h) {
      cos_arr[h] = std::cos(h*phi);
      sin_arr[h] = std::sin(h*phi);
    }
    for (auto p = 1u; p <= Powers; ++p) {
      auto prefactor = std::pow(weight, p);
      ref(0,p) += std::complex<double>(prefactor, 0);
      for (auto h = 1u; h <= Harmonics; ++h) {
        ref(h,p)  += std::complex<double>(prefactor*cos_arr[h],
                                          prefactor*sin_arr[h]);
        ref(-h,p) += std::complex<double>( prefactor*cos_arr[h],
                                          -prefactor*sin_arr[h]);
      }
    }
  }
  void Reset() {
    for (auto &q : qv_) {
      q = std::complex<double>(0., 0.);
    }
  }
  void Print(int nharmonic) {
    int harm = -nharmonic;
    for (auto &q : qv_) {
      std::cout << q << std::endl;
      ++harm;
    }
  }
};

#endif
