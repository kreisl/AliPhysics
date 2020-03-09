/**************************************************************************
 * Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliQvector.h"

AliQvector AliQvector::Recenter(const double x_mean, const double y_mean,
                                const double x_width, double y_width, 
                                bool rescale) const {
  if (rescale) {
    if (x_width > 0. && y_width > 0.) {
     return {(x-x_mean)/x_width, (y-y_mean)/y_width, sum};
    }
    return {x-x_mean, y-y_mean, sum};
  } else {
    return {x-x_mean, y-y_mean, sum};
  }
}
