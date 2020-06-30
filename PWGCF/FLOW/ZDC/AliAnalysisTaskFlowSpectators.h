#ifndef ALIANALYSISTASKFLOWSPECTATORS_H
#define ALIANALYSISTASKFLOWSPECTATORS_H

/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <QnTools/CorrectionManager.hpp>
#include <string>

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliQvector.h"
#include "AliQvectorAlignmentND.h"
#include "AliQvectorCorrectionND.h"
#include "AliQvectorMagnitude.h"
#include "AliZDCResultStorage.h"
#include "AliZDCcumulantFlow.h"
#include "AliZDCellipticFlow.h"
#include "AliZDCgainEq.h"
#include "TTree.h"
#include "yaml-cpp/yaml.h"

class AliAnalysisTaskFlowSpectators : public AliAnalysisTaskSE {
 public:
  static constexpr int kNFilterBits = 10;
  enum Variables {
    kNone,
    kRunNumber,
    kEventNumber,
    kTrigger,
    kTimeStamp,
    kPeriodNumber,
    kOrbitNumber,
    kBunchCrossNumber,
    kCentV0A,
    kCentV0C,
    kCentV0M,
    kCentZNC,
    kCentZNA,
    kCentCL0,
    kCentCL1,
    kNESDTracks,
    kNTracklets,
    kNITSClusters,
    kNTPCTracksHybrid = kNITSClusters + 6,
    kNTPCTracksTPConly,
    kNTPCClusters,
    kV0Mult,
    kVtxX,
    kVtxY,
    kVtxZ,
    kZNSumEnergy,
    kFMDAPhi,
    kFMDAMult = kFMDAPhi + 1200,
    kFMDCPhi = kFMDAMult + 1200,
    kFMDCMult = kFMDCPhi + 1200,
    kV0CChMult = kFMDCMult + 1200,
    kV0CChPhi = kV0CChMult + 32,
    kV0CChRing = kV0CChPhi + 32,
    kV0AChMult = kV0CChRing + 32,
    kV0AChPhi = kV0AChMult + 32,
    kV0AChRing = kV0AChPhi + 32,
    kZDCCChMult = kV0AChRing + 32,
    kZDCCSumMult = kZDCCChMult + 5,
    kZDCCChPhi = kZDCCSumMult + 1,
    kZDCAChMult = kZDCCChPhi + 5,
    kZPAChMult = kZDCAChMult + 5,
    kZPCChMult = kZPAChMult + 5,
    kZPAChOffset = kZPCChMult + 5,
    kZPCChOffset = kZPAChOffset + 5,
    kZPAChPhi = kZPCChOffset + 5,
    kZPCChPhi = kZPAChPhi + 5,
    kZDCASumMult = kZPCChPhi + 5,
    kZDCAChPhi = kZDCASumMult + 1,
    kT0ChMult = kZDCAChPhi + 5,
    kT0ChPhi = kT0ChMult + 24,
    kNEventVariables = kT0ChPhi + 24,
    kPhi = kNEventVariables,
    kPt,
    kPx,
    kPy,
    kPz,
    kEta,
    kDCAxy,
    kDCAz,
    kDCAxySigma,
    kDCAzSigma,
    kCharge,
    kTPCnCls,
    kTPCchi2pCls,
    kFilterBits,
    kNVars = kFilterBits + kNFilterBits
  };

  AliAnalysisTaskFlowSpectators();
  AliAnalysisTaskFlowSpectators(const char *);
  virtual ~AliAnalysisTaskFlowSpectators();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void NotifyRun();
  virtual void Terminate(Option_t *option);
  void SetCorrectionFile(std::string name) { fCorrectionFileName = name; }
  void AddCumulantFlowAnalyses(const YAML::Node &node);
  void AddEllipticFlowAnalyses(const YAML::Node &node);
  void ConfigureBinning(int nbinsxy, int nbinsz, bool equalize, double eta_gap, const std::vector<double> &pt_bins,
                        const std::vector<double> &vtxz_bins, int n_phi_bins, int n_eta_bins, double eta_min,
                        double eta_max) {
    fDelayedNbinsXY = nbinsxy;
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
  Qn::CorrectionManager *GetCorrectionManager() { return fCorrectionManager; }
  void ActivateQnTools() { 
    fActivateQnTools = true; 
    fCorrectionManager = new Qn::CorrectionManager();
    }

 private:
  void ApplyConfiguration();
  void ConfigureCorrectionBinning(int nbinsxy, int nbinsz, bool equalize);
  void ConfigureCumulantAnalysis(double eta_gap, const std::vector<double> &pt_bins,
                                 const std::vector<double> &vtxz_bins, int n_phi_bins, int n_eta_bins, double eta_min,
                                 double eta_max);
  void ConfigureEllipticAnalysis(const std::vector<double> &pt_bins);
  void GetSamples() {
    for (auto &entry : fSamples) entry = fPoisson(fRandomGenerator);
  }
  void AnalyzeEvent(AliAODEvent *event);
  bool MultiplicityCut(AliAODEvent *event);

  Bool_t fActivateQnTools = false;

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

  AliEventCuts fEventCuts;                     //< general event cuts
  AliAnalysisUtils *fAnalysisUtils = nullptr;  //< analysis utils
  AliZDCResultStorage *fStorage = nullptr;     //!<!
  TList *fOutputList = nullptr;
  TList *fCorrelationList;                                //< output list correlations
  TList *fCorrectionList;                                 //< output list corrections
  TList *fQAList;                                         //< output list QA
  std::vector<AliZDCellipticFlow> fEllipticFlowAnalyses;  //< elliptic flow analysis with zdc
  std::vector<AliZDCcumulantFlow> fCumulantFlowAnalyses;  //< elliptic flow analysis with tpc cumulants

  // Pileup cut using multiplicities of TPC only and global tracks
  TH2D *fNtracksESDvsNclsITS = nullptr;
  TH2D *fMultiplicityTPConlyVsGlobal = nullptr;     // before cut
  TH2D *fMultiplicityTPConlyVsGlobalCut = nullptr;  // after cut
  TH2D *fNSigmaTPConlyVsGlobal = nullptr;           // nsigma distribution around mean
  TGraph *fMultMean = nullptr;                      // linear fit to mean of tpc only tracks
  TGraph *fMult3SigmaPlus = nullptr;                // mean + 3sigma
  TGraph *fMult3SigmaMinus = nullptr;               // mean - 3 sigma
  TH1D *fCentralityWeightInput = nullptr;           // centrality weight after pileup cut

  // Subsampling parameters
  Int_t fNsamples = 10;                                         /// Number of samples
  std::vector<Int_t> fSamples = std::vector<Int_t>(fNsamples);  /// samples
  std::mt19937 fRandomGenerator{std::random_device{}()};        /// Random number Gen
  std::poisson_distribution<> fPoisson{1};                      /// distribution of events per sample.

  // Non-uniform acceptance correction for ZDC
  std::string fCorrectionFileName;    //< name of the correction file
  TH1D *fCorrectionStep = nullptr;    //< correction step qa zdc
  TH1D *fCorrectionStepEQ = nullptr;  //< correction step qa gain equalized zdc

  // Gain Equalization histograms
  AliZDCgainEq fGainEqualizationZNA;  //< gain equalization zdc
  AliZDCgainEq fGainEqualizationZNC;  //< gain equalization zdc

  // Recentering correction 4D for VZERO
  AliQvectorCorrectionND fRecenter4DV0A;  //< ND recentering all in one step
  AliQvectorCorrectionND fRecenter4DV0C;  //< ND recentering all in one step

  // Recentering correction 4D for neutron ZDC
  AliQvectorCorrectionND fRecenter4DZNA;  //< ND recentering all in one step
  AliQvectorCorrectionND fRecenter4DZNC;  //< ND recentering all in one step

  // Recentering correction 4D for neutron ZDC
  AliQvectorCorrectionND fRecenter4DAfterGainEqZNA;  //< ND recentering all in one step
  AliQvectorCorrectionND fRecenter4DAfterGainEqZNC;  //< ND recentering all in one step

  // Alignment correction
  AliQvectorAlignmentND fAlignZNA;  //< ND alignment all in one step
  AliQvectorAlignmentND fAlignZNC;  //< ND alignment all in one step

  // Q-vector magnitude for ESE
  AliQvectorMagnitude fQZNAmagnitude;
  AliQvectorMagnitude fQZNCmagnitude;

  // QA histograms
  TH1D *fCentralityV0M = nullptr;       //!<! centrality QA histogram
  TH1D *fCentralityCL1 = nullptr;       //!<! centrality QA histogram
  TH1D *fPsiZA = nullptr;               //!<! eventplane angle QA histogram
  TH1D *fPsiZC = nullptr;               //!<! eventplane angle QA histogram
  TH1D *fPsiZAEQ = nullptr;             //!<! eventplane angle QA histogram
  TH1D *fPsiZCEQ = nullptr;             //!<! eventplane angle QA histogram
  TH1D *fPsiTPC1 = nullptr;             //!<! eventplane angle QA histogram
  TH1D *fPsiTPC2 = nullptr;             //!<! eventplane angle QA histogram
  TH1D *fVertexX = nullptr;             //!<! primary vertex QA histogram
  TH1D *fVertexY = nullptr;             //!<! primary vertex QA histogram
  TH1D *fVertexZ = nullptr;             //!<! primary vertex QA histogram
  TH2D *fVertexXY = nullptr;            //!<! primary vertex QA histogram
  TH2D *fCentralityCL1vsV0M = nullptr;  //!<! centrality correlations cl1 vs v0m

  // QA tree variables
  Qn::CorrectionManager *fCorrectionManager = nullptr;
  TTree *fTree = nullptr;  //!<! values container
  double *fValues;         //!<! values container

  /// \cond CLASSDEF
  ClassDef(AliAnalysisTaskFlowSpectators, 2);
  /// \endcond
};

#endif
