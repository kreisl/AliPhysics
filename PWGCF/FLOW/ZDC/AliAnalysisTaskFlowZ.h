#ifndef ALIANALYSISTASKFLOWZ_H
#define ALIANALYSISTASKFLOWZ_H

/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <string>
#include "yaml-cpp/yaml.h"

#include "TTree.h"

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliZDCgainEq.h"
#include "AliZDCdirectedFlow.h"
#include "AliZDCellipticFlow.h"
#include "AliZDCcumulantFlow.h"
#include "AliQvector.h"
#include "AliQvectorMagnitude.h"
#include "AliQvectorCorrection.h"
#include "AliQvectorCorrectionKD.h"
#include "AliQvectorCorrectionND.h"

class AliAnalysisTaskFlowZ : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskFlowZ();
  AliAnalysisTaskFlowZ(const char*);
  virtual ~AliAnalysisTaskFlowZ();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t*);
  virtual void NotifyRun();
  void SetCorrectionFile(std::string name) { fCorrectionFileName = name; }
  void AddDirectedFlowAnalyses(const YAML::Node &node); 
  void AddCumulantFlowAnalyses(const YAML::Node &node); 
  void AddEllipticFlowAnalyses(const YAML::Node &node);
  void ConfigureCorrectionBinning(int nbinsxy, int nbinsz, bool equalize);
  void ConfigureCumulantAnalysis(double eta_gap,
                                    const std::vector<double> &pt_bins,
                                    const std::vector<double> &vtxz_bins,
                                    int n_phi_bins,
                                    int n_eta_bins, double eta_min, double eta_max);
  void ConfigureEllipticAnalysis(const std::vector<double> &pt_bins);
 private:
  void GetSamples() { for (auto &entry : fSamples) entry = fPoisson(fRandomGenerator); }
  void AnalyzeEvent(AliAODEvent* event);
  void ResetTreeValues();

  // QA tree variables
  TTree   *fTree = nullptr;
  Double_t fBcentV0M = -1.;
  Double_t fBcentCL1 = -1.;
  Double_t fBvtxX = -1.;
  Double_t fBvtxY = -1.;
  Double_t fBvtxZ = -1.;
  Long64_t fRunNumber = 0;
  // Q ZNA
  Double_t fBxZNA = 0.;
  Double_t fByZNA = 0.;
  Double_t fBsZNA = 0.;
  // Q ZNC
  Double_t fBxZNC = 0.;
  Double_t fByZNC = 0.;
  Double_t fBsZNC = 0.;
  // Q ZNA iter
  Double_t fBxZNAiter = 0.;
  Double_t fByZNAiter = 0.;
  Double_t fBsZNAiter = 0.;
  // Q ZNC iter
  Double_t fBxZNCiter = 0.;
  Double_t fByZNCiter = 0.;
  Double_t fBsZNCiter = 0.;
  // Q ZNA all
  Double_t fBxZNAall = 0.;
  Double_t fByZNAall = 0.;
  Double_t fBsZNAall = 0.;
  // Q ZNC all
  Double_t fBxZNCall = 0.;
  Double_t fByZNCall = 0.;
  Double_t fBsZNCall = 0.;
  // Q TPC
  Double_t fBxTPC768 = 0.;
  Double_t fByTPC768 = 0.;
  Double_t fBsTPC768 = 0.;
  // Q TPC
  Double_t fBxTPC96 = 0.;
  Double_t fByTPC96 = 0.;
  Double_t fBsTPC96 = 0.;

  AliEventCuts fEventCuts; //< general event cuts
  TList *fCorrelationList; //< output list correlations
  TList *fCorrectionList;  //< output list corrections
  TList *fQAList;          //< output list QA
  std::vector<AliZDCdirectedFlow> fDirectedFlowAnalyses; //< directed flow analysis with zdc
  std::vector<AliZDCellipticFlow> fEllipticFlowAnalyses; //< elliptic flow analysis with zdc
  std::vector<AliZDCcumulantFlow> fCumulantFlowAnalyses; //< elliptic flow analysis with tpc cumulants

  std::string fCorrectionFileName;   //< name of the correction file
  TH1D *fCorrectionStep = nullptr;   //< correction step qa zdc
  TH1D *fCorrectionStepEQ = nullptr; //< correction step qa gain equalized zdc

  AliZDCgainEq fGainEqZNA; //< gain equalization zdc
  AliZDCgainEq fGainEqZNC; //< gain equalization zdc

  Int_t fNsamples = 10;                                  /// Number of samples
  std::vector<Int_t> fSamples = std::vector<Int_t>(fNsamples); /// samples
  std::mt19937 fRandomGenerator{std::random_device{}()}; /// Random number generator
  std::poisson_distribution<> fPoisson{1};               /// distribution of events per sample.

  AliQvectorCorrectionND fCorrectZNAall; //< ND recentering all in one step
  AliQvectorCorrectionND fCorrectZNCall; //< ND recentering all in one step

  AliQvectorCorrection1D fCorrectZNAEQcentStep1; //< recentering corrections gain equalized zdc
  AliQvectorCorrection1D fCorrectZNCEQcentStep1; //< recentering corrections gain equalized zdc
  AliQvectorCorrectionKD fCorrectZNAEQvXYZStep2; //< recentering corrections gain equalized zdc
  AliQvectorCorrectionKD fCorrectZNCEQvXYZStep2; //< recentering corrections gain equalized zdc
  AliQvectorCorrection1D fCorrectZNAEQcentStep3; //< recentering corrections gain equalized zdc
  AliQvectorCorrection1D fCorrectZNCEQcentStep3; //< recentering corrections gain equalized zdc

  AliQvectorCorrection1DInterpolate fCorrectZNAcentInterStep1; //< recentering corrections zdc
  AliQvectorCorrection1DInterpolate fCorrectZNCcentInterStep1; //< recentering corrections zdc
  AliQvectorCorrection1D fCorrectZNAcentStep1; //< recentering corrections zdc
  AliQvectorCorrection1D fCorrectZNCcentStep1; //< recentering corrections zdc
  AliQvectorCorrectionKD fCorrectZNAvXYZStep2; //< recentering corrections zdc
  AliQvectorCorrectionKD fCorrectZNCvXYZStep2; //< recentering corrections zdc
  AliQvectorCorrection1D fCorrectZNAcentStep3; //< recentering corrections zdc
  AliQvectorCorrection1D fCorrectZNCcentStep3; //< recentering corrections zdc

  AliQvectorMagnitude fQZNAmagnitude;
  AliQvectorMagnitude fQZNCmagnitude;

  // QA histograms
  TH1D* fCentralityV0M = nullptr; //!<! centrality QA histogram
  TH1D* fCentralityCL1 = nullptr; //!<! centrality QA histogram
  TH1D* fPsiZA = nullptr; //!<! eventplane angle QA histogram
  TH1D* fPsiZC = nullptr; //!<! eventplane angle QA histogram
  TH1D* fPsiZAEQ = nullptr; //!<! eventplane angle QA histogram
  TH1D* fPsiZCEQ = nullptr; //!<! eventplane angle QA histogram
  TH1D* fPsiTPC1 = nullptr; //!<! eventplane angle QA histogram
  TH1D* fPsiTPC2 = nullptr; //!<! eventplane angle QA histogram
  TH1D* fVertexX = nullptr; //!<! primary vertex QA histogram
  TH1D* fVertexY = nullptr; //!<! primary vertex QA histogram
  TH1D* fVertexZ = nullptr;  //!<! primary vertex QA histogram
  TH2D* fVertexXY = nullptr; //!<! primary vertex QA histogram
  TH2D* fCentralityCL1vsV0M = nullptr; //!<! centrality correlations cl1 vs v0m

  /// \cond CLASSDEF
  ClassDef(AliAnalysisTaskFlowZ, 4);
  /// \endcond
};

#endif
