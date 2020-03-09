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
#include "AliQvectorCorrection.h"
#include "AliQvectorCorrectionKD.h"

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
