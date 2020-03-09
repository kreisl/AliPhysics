#ifndef ALIQVECTOR_H
#define ALIQVECTOR_H

/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <cmath>
#include "TMath.h"

struct AliQvector {

    double x;
    double y;
    double sum;

    enum class WidthRescale {
      YES,
      NO
    };

    inline void Update(double phi, double weight) {
      x += std::cos(phi) * weight;      
      y += std::sin(phi) * weight;
      sum += weight;
    }

    inline void Normalize() {
      if (sum > 0.) {
        x /= sum; 
        y /= sum;
      }
    }

    inline double Psi() {
      return TMath::Pi() + std::atan2(y,x);
    }

    AliQvector Recenter(double x_mean, double y_mean,
                        double x_width, double y_width,
                        bool rescale) const;
};

struct AliZDCQvectors {
  AliQvector za;
  AliQvector zc;
};

#endif
