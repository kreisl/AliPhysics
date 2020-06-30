#ifndef ALIZDCGAINEQ_H
#define ALIZDCGAINEQ_H

/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <string>
#include <vector>
#include "TProfile.h"
#include "AliQvector.h"

class AliZDCgainEq : public TObject {
  public:
   AliZDCgainEq() = default; 
   AliZDCgainEq(std::string name, UInt_t n_channels, std::vector<UInt_t> groups);
   virtual ~AliZDCgainEq() = default;
   AliQvector ApplyGainEq(const double *channel_energies, const double *channel_phis);
   bool CheckBins(TH1 *histo);
   TProfile* ReadFromOADB(TFile *file, const std::string &oadb_name, Int_t run_number);
   void OpenCorrection(TFile *file, Int_t run_number);
   void AddCorrectionsToList(TList *correction_list, TList *qa_list);
   bool IsApplied() const { return fIsApplied; }
   void Configure(std::string name, UInt_t n_channels, std::vector<UInt_t> groups) {
     fName = name;
     fNChannels = n_channels;
     fGroupEdges = groups;
   }
  private:
   bool fIsApplied = false;
   UInt_t fNChannels;
   std::string fName;
   std::vector<UInt_t> fGroupEdges;
   std::unique_ptr<TProfile> fAverageGainIn; //!<!
   TProfile* fAverageGainOut; //!<!
   TProfile* fGainQAbefore; //!<!
   TProfile* fGainQAafter; //!<!
  /// \cond CLASSDEF
  ClassDef(AliZDCgainEq, 1);
  /// \endcond
};

#endif
