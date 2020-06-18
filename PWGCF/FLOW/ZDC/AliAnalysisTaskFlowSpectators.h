#ifndef ALIANALYSISTASKFLOWSPECTATORS_H
#define ALIANALYSISTASKFLOWSPECTATORS_H

/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <string>
#include "yaml-cpp/yaml.h"

#include "TTree.h"

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliZDCgainEq.h"
#include "AliZDCellipticFlow.h"
#include "AliZDCcumulantFlow.h"
#include "AliQvector.h"
#include "AliQvectorMagnitude.h"
#include "AliQvectorCorrectionND.h"
#include "AliQvectorAlignmentND.h"

class AliAnalysisTaskFlowSpectators : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskFlowSpectators();
  AliAnalysisTaskFlowSpectators(const char*);
  virtual ~AliAnalysisTaskFlowSpectators();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t*);
  virtual void NotifyRun();
  void SetCorrectionFile(std::string name) { fCorrectionFileName = name; }
  void AddCumulantFlowAnalyses(const YAML::Node &node); 
  void AddEllipticFlowAnalyses(const YAML::Node &node);
  void ConfigureBinning(int nbinsxy,
                        int nbinsz,
                        bool equalize,
                        double eta_gap,
                        const std::vector<double> &pt_bins,
                        const std::vector<double> &vtxz_bins,
                        int n_phi_bins,
                        int n_eta_bins, 
                        double eta_min, 
                        double eta_max) {
    fDelayedNbinsXY =  nbinsxy;
    fDelayedNbinsZ = nbinsz;
    fDelayedEqualize = equalize;
    fDelayedEtaGap = eta_gap;
    fDelayedPtBins = pt_bins;
    fDelayedVtxZBins = vtxz_bins;
    fDelayedNphiBins = n_phi_bins;
    fDelayedNetaBins = n_eta_bins;
    fDelayedEtaMin = eta_min;
    fDelayedEtaMax = eta_max;
  }
 private:
  void ApplyConfiguration();
  void ConfigureCorrectionBinning(int nbinsxy, int nbinsz, bool equalize);
  void ConfigureCumulantAnalysis(double eta_gap,
                                    const std::vector<double> &pt_bins,
                                    const std::vector<double> &vtxz_bins,
                                    int n_phi_bins,
                                    int n_eta_bins, double eta_min, double eta_max);
  void ConfigureEllipticAnalysis(const std::vector<double> &pt_bins);
  void GetSamples() { for (auto &entry : fSamples) entry = fPoisson(fRandomGenerator); }
  void AnalyzeEvent(AliAODEvent* event);
  void ResetTreeValues();
  bool MultiplicityCut(AliAODEvent* event);

  // delayed configuration
  Int_t fDelayedNbinsXY;
  Int_t fDelayedNbinsZ;
  Bool_t fDelayedEqualize;
  Double_t fDelayedEtaGap;
  std::vector<Double_t> fDelayedPtBins;
  std::vector<Double_t> fDelayedVtxZBins;
  Int_t fDelayedNphiBins;
  Int_t fDelayedNetaBins;
  Double_t fDelayedEtaMin;
  Double_t fDelayedEtaMax;

  AliEventCuts fEventCuts; //< general event cuts
  AliAnalysisUtils *fAnalysisUtils = nullptr; //< analysis utils
  TList *fCorrelationList; //< output list correlations
  TList *fCorrectionList;  //< output list corrections
  TList *fQAList;          //< output list QA
  std::vector<AliZDCellipticFlow> fEllipticFlowAnalyses; //< elliptic flow analysis with zdc
  std::vector<AliZDCcumulantFlow> fCumulantFlowAnalyses; //< elliptic flow analysis with tpc cumulants


  // Pileup cut using multiplicities of TPC only and global tracks
  TH2D *fNtracksESDvsNclsITS = nullptr;
  TH2D *fMultiplicityTPConlyVsGlobal = nullptr; // before cut
  TH2D *fMultiplicityTPConlyVsGlobalCut = nullptr; // after cut
  TH2D *fNSigmaTPConlyVsGlobal = nullptr; // nsigma distribution around mean
  TGraph *fMultMean = nullptr; // linear fit to mean of tpc only tracks
  TGraph *fMult3SigmaPlus = nullptr; // mean + 3sigma
  TGraph *fMult3SigmaMinus = nullptr; // mean - 3 sigma
  TH1D *fCentralityWeightInput = nullptr; // centrality weight after pileup cut

  // Subsampling parameters
  Int_t fNsamples = 10;                                        /// Number of samples
  std::vector<Int_t> fSamples = std::vector<Int_t>(fNsamples); /// samples
  std::mt19937 fRandomGenerator{std::random_device{}()};       /// Random number Gen
  std::poisson_distribution<> fPoisson{1};                     /// distribution of events per sample.

  // Non-uniform acceptance correction for ZDC
  std::string fCorrectionFileName;   //< name of the correction file
  TH1D *fCorrectionStep = nullptr;   //< correction step qa zdc
  TH1D *fCorrectionStepEQ = nullptr; //< correction step qa gain equalized zdc

  // Gain Equalization histograms
  AliZDCgainEq fGainEqualizationZNA; //< gain equalization zdc
  AliZDCgainEq fGainEqualizationZNC; //< gain equalization zdc

  // Recentering correction 4D for VZERO
  AliQvectorCorrectionND fRecenter4DV0A; //< ND recentering all in one step
  AliQvectorCorrectionND fRecenter4DV0C; //< ND recentering all in one step

  // Recentering correction 4D for neutron ZDC
  AliQvectorCorrectionND fRecenter4DZNA; //< ND recentering all in one step
  AliQvectorCorrectionND fRecenter4DZNC; //< ND recentering all in one step

  // Recentering correction 4D for neutron ZDC
  AliQvectorCorrectionND fRecenter4DAfterGainEqZNA; //< ND recentering all in one step
  AliQvectorCorrectionND fRecenter4DAfterGainEqZNC; //< ND recentering all in one step

  // Alignment correction
  AliQvectorAlignmentND fAlignZNA; //< ND alignment all in one step
  AliQvectorAlignmentND fAlignZNC; //< ND alignment all in one step


  // Q-vector magnitude for ESE
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

  // QA tree variables
  TTree   *fTree = nullptr;
  Double_t fBcentV0M = -1.;
  Double_t fBcentCL1 = -1.;
  Double_t fBvtxX = -1.;
  Double_t fBvtxY = -1.;
  Double_t fBvtxZ = -1.;
  Int_t    fBmultFB128 = 0;
  Int_t    fBmultFB768 = 0;
  Int_t    fBnTracksESD;
  std::vector<Int_t> fBnClsITS;
  Long64_t fRunNumber = 0;
  // Q ZNA
  Double_t fBxZNA = 0.;
  Double_t fByZNA = 0.;
  Double_t fBsZNA = 0.;
  // Q ZNC
  Double_t fBxZNC = 0.;
  Double_t fByZNC = 0.;
  Double_t fBsZNC = 0.;
  // Q ZNA all
  Double_t fBxZNAall = 0.;
  Double_t fByZNAall = 0.;
  Double_t fBsZNAall = 0.;
  // Q ZNC all
  Double_t fBxZNCall = 0.;
  Double_t fByZNCall = 0.;
  Double_t fBsZNCall = 0.;
  // Q TPC
  Double_t fBxTPC768_NUA = 0.;
  Double_t fByTPC768_NUA = 0.;
  Double_t fBx2TPC768_NUA = 0.;
  Double_t fBy2TPC768_NUA = 0.;
  Double_t fBsTPC768_NUA = 0.;

  Double_t fBxTPC768 = 0.;
  Double_t fByTPC768 = 0.;
  Double_t fBsTPC768 = 0.;
  Double_t fBx2TPC768 = 0.;
  Double_t fBy2TPC768 = 0.;
  // Q TPC
  Double_t fBxTPC96 = 0.;
  Double_t fByTPC96 = 0.;
  Double_t fBsTPC96 = 0.;

  /// \cond CLASSDEF
  ClassDef(AliAnalysisTaskFlowSpectators, 2);
  /// \endcond
};

#endif
