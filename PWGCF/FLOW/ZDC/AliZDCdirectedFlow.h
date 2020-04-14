#ifndef ALIDIRECTEDFLOW_H
#define ALIDIRECTEDFLOW_H

/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <iostream>
#include <array>
#include "TObject.h"
#include "AliQvector.h"
#include "AliZDCflowCuts.h"
#include "AliZDCanalysis.h"

class TList;
class TH1D;
class TProfile;
class TProfile2D;
class AliAODEvent;
class AliAODTrack;

class AliZDCdirectedFlow : public AliZDCanalysis {
  static constexpr Int_t fNcentBinsWide = 3;
  static constexpr Int_t fNptBinsWide = 10;
  static constexpr Int_t fNetaBins = 5;
 public:
  AliZDCdirectedFlow() = default;
  AliZDCdirectedFlow(std::string name) : AliZDCanalysis(name) {}
  virtual ~AliZDCdirectedFlow() = default;
  virtual TList *CreateCorrelations();
  virtual void FindCentralityBin(AliAODEvent *event, Double_t centrality, const std::vector<Int_t> &samples);
  void FillPerTrackCorrelations(AliAODTrack *track);
  void FillPerEventCorrelations();

 private:
  std::array<Bool_t, fNcentBinsWide> fIsInCentralityClassWide;
  Double_t fPtBinsWide[fNptBinsWide+1] = {0.15, 0.35, 0.54, 0.94, 1.33, 1.73, 2.12, 2.51, 2.91, 3.7, 4.88};
  Double_t fEtaBins[fNetaBins+1] = {-0.8,-0.48,-0.16,0.16,0.48,0.8};
  
  TProfile* fXaXcCent = nullptr; //!<! Correlation ZNA ZNC Resolution
  TProfile* fYaYcCent = nullptr; //!<! Correlation ZNA ZNC Resolution
  TProfile* fXaYcCent = nullptr; //!<! Correlation ZNA ZNC Resolution
  TProfile* fYaXcCent = nullptr; //!<! Correlation ZNA ZNC Resolution

  TProfile* fXtXaCentEtaP = nullptr; //!<! Correlation TPC ZNA ZNC
  TProfile* fXtXcCentEtaP = nullptr; //!<! Correlation TPC ZNA ZNC
  TProfile* fYtYaCentEtaP = nullptr; //!<! Correlation TPC ZNA ZNC
  TProfile* fYtYcCentEtaP = nullptr; //!<! Correlation TPC ZNA ZNC
  TProfile* fXtXaCentEtaN = nullptr; //!<! Correlation TPC ZNA ZNC
  TProfile* fXtXcCentEtaN = nullptr; //!<! Correlation TPC ZNA ZNC
  TProfile* fYtYaCentEtaN = nullptr; //!<! Correlation TPC ZNA ZNC
  TProfile* fYtYcCentEtaN = nullptr; //!<! Correlation TPC ZNA ZNC

  TH1D* fPt[fNcentBinsWide]; //!<! pt QA histograms

  std::vector<TProfile2D*> fV1XXaEtaCentBS; //!<! v1 vs eta in centrality classes bs samples
  std::vector<TProfile2D*> fV1XXcEtaCentBS; //!<! v1 vs eta in centrality classes bs samples
  std::vector<TProfile2D*> fV1YYaEtaCentBS; //!<! v1 vs eta in centrality classes bs samples
  std::vector<TProfile2D*> fV1YYcEtaCentBS; //!<! v1 vs eta in centrality classes bs samples
  std::vector<TProfile2D*> fV1XXaPtEtaPCentBS; //!<! v1 vs pt in centrality classes bs samples
  std::vector<TProfile2D*> fV1XXcPtEtaPCentBS; //!<! v1 vs pt in centrality classes bs samples
  std::vector<TProfile2D*> fV1YYaPtEtaPCentBS; //!<! v1 vs pt in centrality classes bs samples
  std::vector<TProfile2D*> fV1YYcPtEtaPCentBS; //!<! v1 vs pt in centrality classes bs samples
  std::vector<TProfile2D*> fV1XXaPtEtaNCentBS; //!<! v1 vs pt in centrality classes bs samples
  std::vector<TProfile2D*> fV1XXcPtEtaNCentBS; //!<! v1 vs pt in centrality classes bs samples
  std::vector<TProfile2D*> fV1YYaPtEtaNCentBS; //!<! v1 vs pt in centrality classes bs samples
  std::vector<TProfile2D*> fV1YYcPtEtaNCentBS; //!<! v1 vs pt in centrality classes bs samples

  TProfile2D* fV1XXaEtaCent = nullptr; //!<! v1 vs eta in wide centrality classes 
  TProfile2D* fV1XXcEtaCent = nullptr; //!<! v1 vs eta in wide centrality classes 
  TProfile2D* fV1YYaEtaCent = nullptr; //!<! v1 vs eta in wide centrality classes 
  TProfile2D* fV1YYcEtaCent = nullptr; //!<! v1 vs eta in wide centrality classes 
  TProfile2D* fV1XYaEtaCent = nullptr; //!<! v1 vs eta in wide centrality classes 
  TProfile2D* fV1XYcEtaCent = nullptr; //!<! v1 vs eta in wide centrality classes 
  TProfile2D* fV1YXaEtaCent = nullptr; //!<! v1 vs eta in wide centrality classes 
  TProfile2D* fV1YXcEtaCent = nullptr; //!<! v1 vs eta in wide centrality classes 

  TProfile2D* fV1XXaPtEtaPCent = nullptr; //!<! v1 vs pt in wide centrality classes 
  TProfile2D* fV1XXcPtEtaPCent = nullptr; //!<! v1 vs pt in wide centrality classes 
  TProfile2D* fV1YYaPtEtaPCent = nullptr; //!<! v1 vs pt in wide centrality classes 
  TProfile2D* fV1YYcPtEtaPCent = nullptr; //!<! v1 vs pt in wide centrality classes 
  TProfile2D* fV1XXaPtEtaNCent = nullptr; //!<! v1 vs pt in wide centrality classes 
  TProfile2D* fV1XXcPtEtaNCent = nullptr; //!<! v1 vs pt in wide centrality classes 
  TProfile2D* fV1YYaPtEtaNCent = nullptr; //!<! v1 vs pt in wide centrality classes 
  TProfile2D* fV1YYcPtEtaNCent = nullptr; //!<! v1 vs pt in wide centrality classes 

  TProfile* fV1XXaEta[fNcentBinsWide]; //!<! v1 vs eta in wide centrality classes 
  TProfile* fV1XXcEta[fNcentBinsWide]; //!<! v1 vs eta in wide centrality classes 
  TProfile* fV1YYaEta[fNcentBinsWide]; //!<! v1 vs eta in wide centrality classes 
  TProfile* fV1YYcEta[fNcentBinsWide]; //!<! v1 vs eta in wide centrality classes 
  TProfile* fV1XYaEta[fNcentBinsWide]; //!<! v1 vs eta in wide centrality classes 
  TProfile* fV1XYcEta[fNcentBinsWide]; //!<! v1 vs eta in wide centrality classes 
  TProfile* fV1YXaEta[fNcentBinsWide]; //!<! v1 vs eta in wide centrality classes 
  TProfile* fV1YXcEta[fNcentBinsWide]; //!<! v1 vs eta in wide centrality classes 

  TProfile* fV1XXaPtEtaP[fNcentBinsWide]; //!<! v1 vs pt in wide centrality classes 
  TProfile* fV1XXcPtEtaP[fNcentBinsWide]; //!<! v1 vs pt in wide centrality classes 
  TProfile* fV1YYaPtEtaP[fNcentBinsWide]; //!<! v1 vs pt in wide centrality classes 
  TProfile* fV1YYcPtEtaP[fNcentBinsWide]; //!<! v1 vs pt in wide centrality classes 
  TProfile* fV1XXaPtEtaN[fNcentBinsWide]; //!<! v1 vs pt in wide centrality classes 
  TProfile* fV1XXcPtEtaN[fNcentBinsWide]; //!<! v1 vs pt in wide centrality classes 
  TProfile* fV1YYaPtEtaN[fNcentBinsWide]; //!<! v1 vs pt in wide centrality classes 
  TProfile* fV1YYcPtEtaN[fNcentBinsWide]; //!<! v1 vs pt in wide centrality classes 

 /// \cond CLASSDEF
 ClassDef(AliZDCdirectedFlow, 1);
 /// \endcond
};

#endif
