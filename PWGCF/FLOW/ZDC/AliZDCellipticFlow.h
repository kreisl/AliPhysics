#ifndef ALIZDCELLIPTICFLOW_H
#define ALIZDCELLIPTICFLOW_H

/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <array>
#include <random>
#include "TObject.h"
#include "AliQvector.h"
#include "AliZDCflowCuts.h"
#include "AliZDCanalysis.h"
#include "AliZDCcumulantFlow.h"

class AliAODEvent;
class AliAODTrack;
class TList;
class TProfile;
class TProfile2D;

class AliZDCellipticFlow : public AliZDCanalysis {
 public:
  static constexpr Int_t fNcentBins = 10;
  AliZDCellipticFlow() = default;
  AliZDCellipticFlow(std::string name) : AliZDCanalysis(name) {}
  virtual ~AliZDCellipticFlow() = default;
  virtual TList *CreateCorrelations();
  virtual void FindCentralityBin(AliAODEvent *event, Double_t centrality, const std::vector<Int_t> &samples);
  void FillPerTrackCorrelations(AliAODTrack *track);
  void FillPerEventCorrelations();
  void SetPtBins(const std::vector<double> &pt_bins) {
    fPtBins = pt_bins;
    auto npt = pt_bins.size()-1;
    fAxisPt = new TAxis(npt, pt_bins.data());
  }
  void CalculateCorrelations(Qv<4,4> qtpc,
                             std::vector<Qv<4,4>> qtpcpt,
                             Qv<4,4> qtpcetap,
                             Qv<4,4> qtpcetan);
 private:
  const Double_t centBins[fNcentBins+1] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.,90.};
  std::vector<double> fPtBins;
  TAxis *fCentralityAxis;
  TAxis* fAxisPt = nullptr;       //!<! p_T axis
  Int_t fCentralityBin = -1;

  TProfile* fVXaZXXCent = nullptr; //!<! Correlation V0 ZNA ZNC Resolution
  TProfile* fVXcZXXCent = nullptr; //!<! Correlation V0 ZNA ZNC Resolution
  TProfile* fVXaVXcCent = nullptr; //!<! Correlation V0 ZNA ZNC Resolution
  TProfile* fVYaZYYCent = nullptr; //!<! Correlation V0 ZNA ZNC Resolution
  TProfile* fVYcZYYCent = nullptr; //!<! Correlation V0 ZNA ZNC Resolution
  TProfile* fVYaVYcCent = nullptr; //!<! Correlation V0 ZNA ZNC Resolution

  TProfile* fTXaZXXCent = nullptr; //!<! Correlation TPC ZNA ZNC 3 subevent Resolution
  TProfile* fTXcZXXCent = nullptr; //!<! Correlation TPC ZNA ZNC 3 subevent Resolution
  TProfile* fTXaTXcCent = nullptr; //!<! Correlation TPC ZNA ZNC 3 subevent Resolution
  TProfile* fTYaZYYCent = nullptr; //!<! Correlation TPC ZNA ZNC 3 subevent Resolution
  TProfile* fTYcZYYCent = nullptr; //!<! Correlation TPC ZNA ZNC 3 subevent Resolution
  TProfile* fTYaTYcCent = nullptr; //!<! Correlation TPC ZNA ZNC 3 subevent Resolution

  TProfile* fXaXcCent = nullptr; //!<! Correlation ZNA ZNC Resolution
  TProfile* fYaYcCent = nullptr; //!<! Correlation ZNA ZNC Resolution
  TProfile* fXaYcCent = nullptr; //!<! Correlation ZNA ZNC Resolution
  TProfile* fYaXcCent = nullptr; //!<! Correlation ZNA ZNC Resolution

  TH2D* fXaXcDistCent = nullptr; //!<! Correlation ZNA ZNC Distribution
  TH2D* fYaYcDistCent = nullptr; //!<! Correlation ZNA ZNC Distribution
  TH2D* fXaYcDistCent = nullptr; //!<! Correlation ZNA ZNC Distribution
  TH2D* fYaXcDistCent = nullptr; //!<! Correlation ZNA ZNC Distribution

  TH2D* fXtXaXcDistCent = nullptr; //!<! Correlation TPC ZNA ZNC Distribution 
  TH2D* fXtYaYcDistCent = nullptr; //!<! Correlation TPC ZNA ZNC Distribution
  TH2D* fYtXaYcDistCent = nullptr; //!<! Correlation TPC ZNA ZNC Distribution
  TH2D* fYtYaXcDistCent = nullptr; //!<! Correlation TPC ZNA ZNC Distribution

  TProfile* fXtXaXcCent = nullptr; //!<! Correlation TPC ZNA ZNC
  TProfile* fXtYaYcCent = nullptr; //!<! Correlation TPC ZNA ZNC
  TProfile* fYtXaYcCent = nullptr; //!<! Correlation TPC ZNA ZNC
  TProfile* fYtYaXcCent = nullptr; //!<! Correlation TPC ZNA ZNC
  TProfile* fYtYaYcCent = nullptr; //!<! Correlation TPC ZNA ZNC
  TProfile* fXtXaYcCent = nullptr; //!<! Correlation TPC ZNA ZNC
  TProfile* fXtYaXcCent = nullptr; //!<! Correlation TPC ZNA ZNC
  TProfile* fYtXaXcCent = nullptr; //!<! Correlation TPC ZNA ZNC

  TH1D* fPt[fNcentBins]; //!<! pt QA histograms

  TH2D*       fPtCent;      //!<! pt in centrality classes 
  TProfile2D* fV2XXXpTcent; //!<! v2 vs pt in centrality classes 
  TProfile2D* fV2XYYpTcent; //!<! v2 vs pt in centrality classes
  TProfile2D* fV2YXYpTcent; //!<! v2 vs pt in centrality classes
  TProfile2D* fV2YYXpTcent; //!<! v2 vs pt in centrality classes
  TProfile2D* fV2YYYpTcent; //!<! v2 vs pt in centrality classes
  TProfile2D* fV2YXXpTcent; //!<! v2 vs pt in centrality classes
  TProfile2D* fV2XXYpTcent; //!<! v2 vs pt in centrality classes
  TProfile2D* fV2XYXpTcent; //!<! v2 vs pt in centrality classes

  std::vector<TProfile2D*> fV2XXXpTcentBS; //!<! v2 vs pt in centrality classes bootstrap samples 
  std::vector<TProfile2D*> fV2XYYpTcentBS; //!<! v2 vs pt in centrality classes bootstrap samples
  std::vector<TProfile2D*> fV2YXYpTcentBS; //!<! v2 vs pt in centrality classes bootstrap samples
  std::vector<TProfile2D*> fV2YYXpTcentBS; //!<! v2 vs pt in centrality classes bootstrap samples

  std::vector<TProfile*> fXaXcCentBS; //!<! Correlation ZNA ZNC Resolution bootstrap samples
  std::vector<TProfile*> fYaYcCentBS; //!<! Correlation ZNA ZNC Resolution bootstrap samples
  std::vector<TProfile*> fXaYcCentBS; //!<! Correlation ZNA ZNC Resolution bootstrap samples
  std::vector<TProfile*> fYaXcCentBS; //!<! Correlation ZNA ZNC Resolution bootstrap samples

  TProfile* fV2XXXpT[fNcentBins]; //!<! v2 vs pt in centrality classes 
  TProfile* fV2XYYpT[fNcentBins]; //!<! v2 vs pt in centrality classes
  TProfile* fV2YXYpT[fNcentBins]; //!<! v2 vs pt in centrality classes
  TProfile* fV2YYXpT[fNcentBins]; //!<! v2 vs pt in centrality classes
  TProfile* fV2YYYpT[fNcentBins]; //!<! v2 vs pt in centrality classes
  TProfile* fV2YXXpT[fNcentBins]; //!<! v2 vs pt in centrality classes
  TProfile* fV2XXYpT[fNcentBins]; //!<! v2 vs pt in centrality classes
  TProfile* fV2XYXpT[fNcentBins]; //!<! v2 vs pt in centrality classes

 /// \cond CLASSDEF
 ClassDef(AliZDCellipticFlow, 1);
 /// \endcond
};

#endif
