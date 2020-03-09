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
#include "TProfile3D.h"
#include "AliQvector.h"
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
  template <unsigned int Harmonics, unsigned int Powers>
  struct Qv {
    static const unsigned int kMaxH = Harmonics;
    static const unsigned int kMaxP = Powers;
    static constexpr std::size_t kSize = (kMaxH*2+1)*(kMaxP); 
    std::array<std::complex<double>, kSize> qv_;
    std::complex<double> operator()(const int ih, const int ip) const {
      return qv_[index(ih,ip)];
    }
    std::complex<double> &ref(const int ih, const int ip) {
      return qv_[index(ih,ip)];
    }
    unsigned int index(const int ih, const int ip) const {
      auto x = (ih + kMaxH) * (kMaxP) + (ip-1);
      return x;
    }
    void Fill(double phi, double weight) {
      std::array<double, Harmonics+1> cos_arr;
      std::array<double, Harmonics+1> sin_arr;
      for (auto h = 0u; h <= Harmonics; ++h) {
        cos_arr[h] = std::cos(h*phi);
        sin_arr[h] = std::sin(h*phi);
      }
      for (auto p = 1u; p <= Powers; ++p) {
        for (auto h = 0u; h <= Harmonics; ++h) {
          ref(h,p)  += std::complex<double>(std::pow(weight, p)*cos_arr[h],
                                            std::pow(weight, p)*sin_arr[h]);
          if (h == 0) continue;
          ref(-h,p) += std::complex<double>( std::pow(weight, p)*cos_arr[h],
                                            -std::pow(weight, p)*sin_arr[h]);
        }
      }
    }
    void Reset() {
      for (auto & q : qv_) {
        q = std::complex<double>(0., 0.);
      }
    }
  };
 public:
  AliZDCcumulantFlow() = default;
  AliZDCcumulantFlow(std::string name);
  virtual ~AliZDCcumulantFlow() {
    delete fAxisPt;
    delete fNUAweightsIn;
    delete fPtEfficiencyIn;
  }
  void FindCentralityBin(AliAODEvent *event, Double_t centrality, const std::vector<Int_t> &samples);
  void OpenCorrection(TFile *file, Int_t run_number); 
  void AddCorrectionsToList(TList *hlist, TList *qalist);
  void Configure(double eta_gap, std::vector<double> pt_bins,
                 std::vector<double> vtxz_bins,
                 int n_phi_bins, int n_eta_bins, double eta_min, double eta_max);
  TList *CreateCorrelations();
  void FillPerTrackCorrelations(AliAODTrack *track);
  bool IsApplied() const { return fIsApplied; }
  void CalculateCumulants();
  virtual void ReadYAMLnode(YAML::Node& node) {
    auto analysis_settings = node["analysis_settings"];
    if (analysis_settings.IsDefined()) {
      ParseFlag(fApplyNUA, analysis_settings, "nua_weights");
    }
    fCuts.ReadYAMLnode(node);
  }
  void SetBootStrapSamples(int n) {
    fNsamples = n;
    fSamples.resize(n);
  }

 private:
  TProfile3D* ReadFromOADB(TFile *file, const std::string&oadb_name, Int_t run_number);

  bool fUseEtaGap = false;
  bool fIsApplied = false;
  bool fApplyNUA  = true;
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
  TAxis* fAxisPt = nullptr;              //!<! p_T axis
  TProfile3D* fNUAweightsIn = nullptr;         //!<! nua weight eta phi vtx z
  TH1D* fPtEfficiencyIn = nullptr;       //!<! pt efficiency
  TH1D* fAreNUAapplied = nullptr;
  TH3D* fNUAweightsNtracksTemp = nullptr; //!<! nua weight eta phi vtx z
  TH3D* fAfterNUA = nullptr; //!<! nua weight eta phi vtx z
  TH3D* fBeforeNUA = nullptr; //!<! nua weight eta phi vtx z
  TH1D* fFilterBit = nullptr;
  TProfile3D* fNUAweightsOut = nullptr; //!<! nua weight eta phi vtx z
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

  TProfile *fSinTerms = nullptr;         //!<! Sin terms sin(1) | sin(1+2) | sin(1-2-3)
  TProfile *fCosTerms = nullptr;         //!<! Cos terms cos(1) | cos(1+2) | cos(1-2-3)
  TH2D *fMultC24 = nullptr;          //!<! C_2{4} multiplicity QA
  TH2D *fC22Distribution = nullptr;      //!<! C_2{2} distribution
  TH2D *fC24Distribution = nullptr;      //!<! C_2{4} distribution

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
