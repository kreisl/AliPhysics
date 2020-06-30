#ifndef ALIQVECTORCORRECTIONKD_H
#define ALIQVECTORCORRECTIONKD_H

/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <string>
#include "TAxis.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TH2Poly.h"
#include "AliQvector.h"
#include "TFile.h"
#include "AliOADBContainer.h"

class AliQvectorCorrectionKD : public TObject {
 public:
  AliQvectorCorrectionKD() = default;
  void Configure(const std::string &name, const std::vector<std::string> &axes_names,
                 int nbinsxy, int nbinsz, bool equalize);
  AliQvector Correct(const AliQvector &q_vector, const std::vector<double> &vars) const;
  bool IsApplied() const { return fIsApplied; }
  void OpenCorrection(TFile *file, Int_t run);
  void AddCorrectionsToList(TList *hlist, TList *qalist);
 private:
  bool fEqualization = false;
  bool fIsConfigured = false;
  bool fIsApplied = false;
  int fNbinsZ;
  int fNbinsXY;
  int fDim;
  std::string fName;
  std::vector<std::string> fAxesNames;
  std::unique_ptr<TAxis> fAxisZ; //!<!
  std::vector<std::unique_ptr<TH2Poly>> fBinnings; //!<!
  std::unique_ptr<TProfile> fXmeanIn; //!<!
  std::unique_ptr<TProfile> fYmeanIn; //!<!
  TProfile *fXmeanOut = nullptr; //!<!
  TProfile *fYmeanOut = nullptr; //!<!
  TH1D *fCorrectionEntries = nullptr; //!<!

  TProfile* fXmeanQAxyz = nullptr; //!<!
  TProfile* fYmeanQAxyz = nullptr; //!<!
  TH2D* fXmeanQADistributionxyz = nullptr; //!<!
  TH2D* fYmeanQADistributionxyz = nullptr; //!<!
    
  bool CheckBins(TH1 *histo);
  template<typename T>
  T* ReadFromOADB(TFile *file, const std::string &oadb_name, Int_t run_number) {
    T* out = nullptr;
    AliOADBContainer* oadb = dynamic_cast<AliOADBContainer*>(file->Get(oadb_name.data()));
    if (oadb) {
      auto in = dynamic_cast<T*>(oadb->GetObject(run_number));
      if (in) out = in;
    }
    return out;
  }
  
 /// \cond CLASSDEF
 ClassDef(AliQvectorCorrectionKD, 2);
 /// \endcond

};

#endif
