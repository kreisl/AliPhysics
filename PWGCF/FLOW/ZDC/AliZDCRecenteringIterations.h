#ifndef ALIZDCRECENTERINGITERATIONS_H
#define ALIZDCRECENTERINGITERATIONS_H

#include "AliQvectorCorrectionND.h"

class AliZDCRecenteringIterations : public TObject {
 private:
  AliQvectorCorrectionND fRecenterCoarse4D1; //<
  AliQvectorCorrectionND fRecenterFineCentrality; //<
  AliQvectorCorrectionND fRecenterCoarse4D2; //<
  AliQvectorCorrectionND fRecenterFineVtxX; //<
  AliQvectorCorrectionND fRecenterCoarse4D3; //<
  AliQvectorCorrectionND fRecenterFineVtxY; //<
  AliQvectorCorrectionND fRecenterCoarse4D4; //<
  AliQvectorCorrectionND fRecenterFineVtxZ; //<
 public:
  void Make(TFile *file, int run_number, TList *corrections, TList *qa);
  AliQvector Apply();
  bool IsApplied() {return fIsApplied;}

}

#endif