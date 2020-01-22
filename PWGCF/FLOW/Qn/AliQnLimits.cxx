#include "AliQnLimits.h"

Long64_t AliQnLimits::Merge(TCollection *list) {
  if (!list) return 0;
  TIter next(list);
  while (auto limit = (AliQnLimits*) next()) {
    SetNew(limit->Min());
    SetNew(limit->Max());
  }
  return 1;
}
