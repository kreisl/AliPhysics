#ifndef ALIZDCCUMULANTFLOW_H
#define ALIZDCCUMULANTFLOW_H

/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <array>
#include <memory>
#include <complex>
#include <random>
#include <vector>
#include <iostream>
#include "TObject.h"
#include "TAxis.h"
#include "TH1.h"
#include "TSpline.h"
#include "TF1.h"
#include "TProfile3D.h"
#include "AliQvector.h"
#include "AliQvectors.h"
#include "AliZDCflowCuts.h"
#include "AliZDCanalysis.h"
#include "AliAODTrack.h"

#include "yaml-cpp/yaml.h"

class TList;
class TFile;
class TProfile;
class TProfile2D;
class AliAODEvent;

class AliZDCcumulantFlow : public AliZDCanalysis {
 public:
  AliZDCcumulantFlow() = default;
  AliZDCcumulantFlow(std::string name);
  virtual ~AliZDCcumulantFlow() {
    delete fAxisPt;
    delete fNUAweightsIn;
    for (auto spline : fPtEfficiencySplinesCentralityClasses) {
      delete spline;
    }
  }
  virtual void FindCentralityBin(AliAODEvent *event, std::vector<Double_t> centralities, const std::vector<Int_t> &samples);
  void OpenCorrection(TFile *file, Int_t run_number); 
  void AddCorrectionsToList(TList *hlist, TList *qalist);
  void Configure(double eta_gap, std::vector<double> pt_bins,
                 std::vector<double> vtxz_bins,
                 int n_phi_bins, int n_eta_bins, double eta_min, double eta_max);
  TList *CreateCorrelations();
  void FillPerTrackCorrelations(AliAODTrack *track);
  void CalculateCumulants();
  void FillESE(double qzna, double qznc);
  TSpline3* GetPtSplineIntegrated() { return fPtEfficiencySpline; }

  void SetBootStrapSamples(int n) {
    fNsamples = n;
    fSamples.resize(n);
  }
  void SetESE(bool zdc = true) { fESE = zdc;}

  virtual void ReadYAMLnode(YAML::Node& node) {
    auto analysis_settings = node["analysis_settings"];
    if (analysis_settings.IsDefined()) {
      ParseFlag(fApplyNUA, analysis_settings, "nua_weights");
      ParseFlag(fApplyNUE, analysis_settings, "nue_weights");
      ParseFlag(fUseNUEintegrated, analysis_settings, "nue_weights_integrated");
    }
    fCuts.ReadYAMLnode(node);
  }

  void OpenPtEfficiencies(TFile *file);

  std::vector<Qv<4,4>> GetQTPCpt() {return fP;}
  Qv<4,4> GetQTPC() {return fR;}

  Qv<4,4> GetQTPCEtaNeg() {return fRptCutetaGapN;}
  Qv<4,4> GetQTPCEtaPos() {return fRptCutetaGapP;}

 private:
  unsigned int GetPtEfficiencyCentralityBin(double c) {
    unsigned int i = 9;
    if ( c < 5.) i = 0;
    else if ( 5. < c && c < 10.) i = 1;
    else if (10. < c && c < 20.) i = 2;
    else if (20. < c && c < 30.) i = 3;
    else if (30. < c && c < 40.) i = 4;
    else if (40. < c && c < 50.) i = 5;
    else if (50. < c && c < 60.) i = 6;
    else if (60. < c && c < 70.) i = 7;
    else if (c > 70.) i = 8;
    return i;
  }
  void ResetQvectors();
  TH3D* ReadFromOADB(TFile *file, const std::string&oadb_name, Int_t run_number);
  void ScaleNUA(TH3D* unscaled, TH3D* scaled);

  bool fUseNUEintegrated = true;
  bool fUseEtaGap = false;
  bool fApplyNUA  = true;
  bool fApplyNUE  = true;
  int  fPtEffCentBin = 9;
  double fEtaGap = 0.5;
  double fPtMax = -1.;
  double fPtMin = -1.;
  Double_t fVertexZ = 0.;


  Int_t fNsamples = 10;                                  /// Number of samples
  std::vector<Int_t> fSamples = std::vector<Int_t>(fNsamples);

  Qv<4,4> fRptCut;                //!<!
  Qv<4,4> fRptCutetaGapP;         //!<!
  Qv<4,4> fRptCutetaGapN;         //!<!
  Qv<4,4> fR;                     //!<!
  std::vector<Qv<4,4>> fP;        //!<!
  std::vector<Qv<4,4>> fQ;        //!<!
  Qv<4,4> fRetaGapN;              //!<!
  std::vector<Qv<4,4>> fPetaGapP; //!<!
  std::vector<Qv<4,4>> fQetaGapP; //!<!
  TAxis* fAxisPt = nullptr;       //!<! p_T axis
  std::vector<TF1*> fPtEfficienciesFit; //!<! pt efficiency
  TSpline3* fPtEfficiencySpline = nullptr; //!<! pt efficiency
  std::vector<TSpline3*> fPtEfficiencySplinesCentralityClasses; //!<!
  TH1D* fAreWeightsApplied = nullptr;
  TH3D* fNUAweightsIn = nullptr;
  TH3D* fNUAweightsScaled = nullptr;
  TH3D* fNUAweightsOut = nullptr;
  TH3D* fAfterNUA = nullptr;              //!<! nua weight eta phi vtx z
  TH3D* fBeforeNUA = nullptr;             //!<! nua weight eta phi vtx z
  TH1D* fFilterBit = nullptr;
  TProfile *fC22 = nullptr;              //!<! C_2{2}
  TProfile *fC22EtaGap = nullptr;        //!<! C_2{2}_{#Delta#eta>1.0}
  TProfile *fC24 = nullptr;              //!<! C_2{4}
  TProfile2D *fCdif22 = nullptr;         //!<! C'_2{2}
  TProfile2D *fCdif22EtaGap = nullptr;   //!<! C_2{2}_{#Delta#eta>1.0}
  TProfile2D *fCdif24 = nullptr;         //!<! C'_2{4}

  std::vector<TProfile *> fC22BS;              //!<! C_2{2} bootstrap samples
  std::vector<TProfile *> fC22EtaGapBS;        //!<! C_2{2}_{#Delta#eta>1.0} bootstrap samples
  std::vector<TProfile *> fC24BS;              //!<! C_2{4} bootstrap samples
  std::vector<TProfile2D *> fCdif22BS;         //!<! C'_2{2} bootstrap samples
  std::vector<TProfile2D *> fCdif22EtaGapBS;   //!<! C_2{2}_{#Delta#eta>1.0} bootstrap samples
  std::vector<TProfile2D *> fCdif24BS;         //!<! C'_2{4} bootstrap samples

  bool fESE;
  double fQpercentileZNA;
  double fQpercentileZNC;

  TProfile2D *fC22ESE = nullptr;              //!<! C_2{2} q selection
  TProfile2D *fC22EtaGapESE = nullptr;        //!<! C_2{2}_{#Delta#eta>1.0} q selection
  TProfile2D *fC24ESE = nullptr;              //!<! C_2{4} q selection

  TProfile *fSinTerms = nullptr;         //!<! Sin terms sin(1) | sin(1+2) | sin(1-2-3)
  TProfile *fCosTerms = nullptr;         //!<! Cos terms cos(1) | cos(1+2) | cos(1-2-3)
  TH2D *fMultC24 = nullptr;          //!<! C_2{4} multiplicity QA
  TH2D *fC22Distribution = nullptr;      //!<! C_2{2} distribution
  TH2D *fC22DistributionWhole = nullptr;      //!<! C_2{4} distribution
  TH2D *fC24Distribution = nullptr;      //!<! C_2{4} distribution
  TH2D *fC24DistributionWhole = nullptr;      //!<! C_2{4} distribution

  void FillFilterBitQA(AliAODTrack* track) {
    if (track->TestFilterBit(1))   fFilterBit->Fill(0.);
    if (track->TestFilterBit(2))   fFilterBit->Fill(1.);
    if (track->TestFilterBit(4))   fFilterBit->Fill(2.);
    if (track->TestFilterBit(8))   fFilterBit->Fill(3.);
    if (track->TestFilterBit(16))  fFilterBit->Fill(4.);
    if (track->TestFilterBit(32))  fFilterBit->Fill(5.);
    if (track->TestFilterBit(64))  fFilterBit->Fill(6.);
    if (track->TestFilterBit(128)) fFilterBit->Fill(7.);
    if (track->TestFilterBit(256)) fFilterBit->Fill(8.);
    if (track->TestFilterBit(512)) fFilterBit->Fill(9.);
    if (track->IsHybridGlobalConstrainedGlobal()) fFilterBit->Fill(10.);
  }


  template<typename QQ>
  std::complex<double> TwoParticleRef(const QQ& q, const int n) {
    return q(n,1) * q(-n,1) - q(0,2);
  }

  template<typename QQ>
  std::complex<double> TwoParticleEtaGapRef(const QQ& q1, const QQ& q2, const int n) {
    return q1(n,1) * q2(-n,1);
  }
  
  template<typename QQ>
  std::complex<double> TwoParticleDif(const QQ&q, const QQ&p, const QQ&r, const int n) {
    return p(n,1) * r(-n,1) - q(0,2);
  }
  
  template<typename QQ>
  std::complex<double> FourParticleRef(const QQ& q, const int n) {
    return  2. * q(   0,2) * q(   0,2)
          - 4. * q(   0,2) * q(  -n,1) * q(   n,1)
          - 6. * q(   0,4)               
          -      q(-2*n,2) * q(   n,1) * q(   n,1)
          +      q(-2*n,2) * q( 2*n,2)   
          +      q(  -n,1) * q(  -n,1) * q(   n,1) * q(   n,1)
          -      q(  -n,1) * q(  -n,1) * q( 2*n,2)
          + 4. * q(  -n,1) * q(   n,3) 
          + 4. * q(  -n,3) * q(   n,1);
  }

  template<typename QQ>
  std::complex<double> FourParticleDif(const QQ&q, const QQ&p, const QQ&r, const int n) {
       return  p(   n,1) * r(   n,1) * r(  -n,1) * r(  -n,1)
        -      q( 2*n,2) * r(  -n,1) * r(  -n,1)
        -      p(   n,1) * r(   n,1) * r(-2*n,2)
        - 2. * r(   n,1) * q(   0,2) * r(  -n,1)
        - 2. * p(   n,1) * r(   0,2) * r(  -n,1)
        + 4. * q(   n,3) * r(  -n,1)
        + 2. * r(   0,2) * q(   0,2)
        + 2. * r(   n,1) * q(  -n,3)
        + 2. * p(   n,1) * r(  -n,3)
        +      q( 2*n,2) * r(-2*n,2)
        - 6. * q(   0,4);
  }

 /// \cond CLASSDEF
 ClassDef(AliZDCcumulantFlow, 3);
 /// \endcond
};

#endif
