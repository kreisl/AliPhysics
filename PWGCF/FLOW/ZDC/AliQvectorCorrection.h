#ifndef ALIQVECTORCORRECTION_H
#define ALIQVECTORCORRECTION_H

/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <string>
#include "TAxis.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "AliQvector.h"

class AliQvectorCorrection : public TObject {
 public:
  AliQvectorCorrection() = default;
  AliQvectorCorrection(std::string name, bool equalize);
  virtual AliQvector Correct(const AliQvector &q_vector, double var_a, double var_b = 0., double var_c = 0.) const = 0;
  virtual void OpenCorrection(TFile *corrections_file, Int_t run_number) = 0;
  virtual void AddCorrectionsToList(TList *hlist, TList *qalist) = 0;
  bool IsApplied() const { return fIsApplied; }
 protected:
  bool CheckBins(TH1 *histo);
  bool fEqualization = false;
  bool fIsApplied = false;
  bool fIsConfigured = false;
  std::string fName;
};

class AliQvectorCorrection1DInterpolate : public AliQvectorCorrection {
 public:
  AliQvectorCorrection1DInterpolate() = default;
  void Configure(std::string name, std::string var_name,
                 std::vector<Double_t> bin_edges, bool equalize);
  AliQvector Correct(const AliQvector &q_vector, double var_a, double var_b = 0., double var_c = 0.) const final;
  void OpenCorrection(TFile *file, Int_t run_number) final;
  void AddCorrectionsToList(TList *hlist, TList *qalist) final;
 private:
  TF1* ReadFromOADB(TFile *file, const std::string &oadb_name, Int_t run_number);
  std::string fVarName;
  std::vector<Double_t> fBinEdges;
  TProfile *fXmeanOut = nullptr; //!<!
  TProfile *fYmeanOut = nullptr; //!<!
  TProfile *fXmeanQA = nullptr; //!<!
  TProfile *fYmeanQA = nullptr; //!<!
  TH1D *fCorrectionEntries = nullptr; //!<!
  TH2D *fXDistributionQA = nullptr; //!<!
  TH2D *fYDistributionQA = nullptr; //!<!
  std::unique_ptr<TF1> fXmeanIn; //!<!
  std::unique_ptr<TF1> fYmeanIn; //!<!

 /// \cond CLASSDEF
 ClassDef(AliQvectorCorrection1DInterpolate, 2);
 /// \endcond
};

class AliQvectorCorrection1D : public AliQvectorCorrection {
 public:
  AliQvectorCorrection1D() = default;
  void Configure(std::string name, std::string var_name,
                 std::vector<Double_t> bin_edges, bool equalize);
  AliQvector Correct(const AliQvector &q_vector, double var_a, double var_b = 0., double var_c = 0.) const final;
  void OpenCorrection(TFile *file, Int_t run_number) final;
  void AddCorrectionsToList(TList *hlist, TList *qalist) final;
 private:
  TProfile* ReadFromOADB(TFile *file, const std::string &oadb_name, Int_t run_number);
  std::string fVarName;
  std::vector<Double_t> fBinEdges;
  TProfile *fXmeanOut = nullptr; //!<!
  TProfile *fYmeanOut = nullptr; //!<!
  TProfile *fXmeanQA = nullptr; //!<!
  TProfile *fYmeanQA = nullptr; //!<!
  TH1D *fCorrectionEntries = nullptr; //!<!
  TH2D *fXDistributionQA = nullptr; //!<!
  TH2D *fYDistributionQA = nullptr; //!<!
  std::unique_ptr<TProfile> fXmeanIn; //!<!
  std::unique_ptr<TProfile> fYmeanIn; //!<!

 /// \cond CLASSDEF
 ClassDef(AliQvectorCorrection1D, 3);
 /// \endcond
};

class AliQvectorCorrection2D : public AliQvectorCorrection {
 public:
  AliQvectorCorrection2D() = default;
  void Configure(std::string name, 
                 std::string var_a_name, std::vector<Double_t> var_a_bin_edges,
                 std::string var_b_name, std::vector<Double_t> var_b_bin_edges,
                 bool equalize);
  AliQvector Correct(const AliQvector &q_vector, double var_a, double var_b = 0., double var_c = 0.) const final;
  void OpenCorrection(TFile *file, Int_t run_number) final;
  void AddCorrectionsToList(TList *hlist, TList *qalist) final;
 private:
  TProfile2D* ReadFromOADB(TFile *file, const std::string &oadb_name, Int_t run_number);
  std::string fNameVarA;
  std::string fNameVarB;
  std::vector<Double_t> fBinEdgesVarA;
  std::vector<Double_t> fBinEdgesVarB;
  std::unique_ptr<TProfile2D> fXmeanIn; //!<!
  std::unique_ptr<TProfile2D> fYmeanIn; //!<!
  TProfile2D *fXmeanOut = nullptr; //!<!
  TProfile2D *fYmeanOut = nullptr; //!<!
  TProfile2D *fXmeanQA = nullptr; //!<!
  TProfile2D *fYmeanQA = nullptr; //!<!
  TH2D *fCorrectionEntries = nullptr; //!<!
  TH2D *fXDistributionQAvsA = nullptr; //!<!
  TH2D *fYDistributionQAvsA = nullptr; //!<!
  TH2D *fXDistributionQAvsB = nullptr; //!<!
  TH2D *fYDistributionQAvsB = nullptr; //!<!
  TProfile *fXMeanQAvsA = nullptr; //!<!
  TProfile *fYMeanQAvsA = nullptr; //!<!
  TProfile *fXMeanQAvsB = nullptr; //!<!
  TProfile *fYMeanQAvsB = nullptr; //!<!
 
 /// \cond CLASSDEF
 ClassDef(AliQvectorCorrection2D, 3);
 /// \endcond
};

#endif
