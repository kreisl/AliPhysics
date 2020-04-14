#ifndef ALIQVECTORMAGNITUDE_H
#define ALIQVECTORMAGNITUDE_H

/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <string>

#include "TH2.h"
#include "TSpline.h"
#include "TString.h"

#include "AliQvector.h"

class TList;
class TFile;
class TClonesArray;

class AliQvectorMagnitude : public TObject {
 public:
  void Fill(const AliQvector &qvector, double centrality);
  double GetPercentile(const AliQvector &qvector, double centrality) const;
  void Configure(const std::string &name, int nbins, double max);
  void Initialize();
  void ReadFile(TFile *file);
  void AddToOutputList(TList* list);
 private:
  std::string fName;
  TAxis *fCentralityAxis = nullptr; //!<!
  TH2F *fDistributionOut = nullptr; //!<! 
  TH2F *fDistributionIn = nullptr;  //!<!

  std::vector<std::unique_ptr<TSpline3>> fSplines; //!<!

 /// \cond CLASSDEF
 ClassDef(AliQvectorMagnitude, 2);
 /// \endcond
};

#endif

