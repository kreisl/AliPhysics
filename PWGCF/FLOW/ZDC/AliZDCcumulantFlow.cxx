/**************************************************************************
 * Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include <iostream>
#include <cmath>
#include "TProfile.h"
#include "TList.h"
#include "TFile.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TH3.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliOADBContainer.h"
#include "AliZDCcumulantFlow.h"

AliZDCcumulantFlow::AliZDCcumulantFlow(std::string name) : 
    AliZDCanalysis(name) {
  }

void AliZDCcumulantFlow::FindCentralityBin(AliAODEvent *event, Double_t centrality, const std::vector<Int_t> &samples) {
  fCuts.CheckEventCuts(event); 
  fCentrality = centrality;
  const AliAODVertex* vtx = event->GetPrimaryVertex();
  fVertexZ = vtx->GetZ();
  fCuts.HasPassedCentrality(centrality > 0. && centrality < 90.);
  fSamples = samples;
}

void AliZDCcumulantFlow::FillPerTrackCorrelations(AliAODTrack *track) { 
  if (!(fCuts.PassedEventCuts() && fCuts.CheckTrackCutsNoPtCut(track))) return;
  bool is_reference = false;
  bool is_reference_eta_gap_p = false;
  bool is_reference_eta_gap_n = false;
  bool is_POI = false;
  bool is_POI_eta_gap = false;
  const Double_t phi = track->Phi();
  const Double_t eta = track->Eta();
  const Double_t pt = track->Pt(); 
  double weight = 1.;
  if (pt > fPtMax || pt < fPtMin) return;
  is_reference = true;
  if (fApplyNUA && fNUAweightsIn) {
    auto w = fNUAweightsIn->GetBinContent(fNUAweightsIn->FindBin(phi, eta, fVertexZ));
    if (w > 0.) weight *= 1. / w;
    else return;
  }
  if (fApplyNUA && fPtEfficiencyIn) {
    auto w = fPtEfficiencyIn->GetBinContent(fPtEfficiencyIn->FindBin(pt));
    if (w > 0.) weight *= 1. / w;
    else return;
  }
  if (fUseEtaGap) {
    if      ( eta < -fEtaGap) is_reference_eta_gap_n = true;
    else if ( eta >  fEtaGap) is_reference_eta_gap_p = true;
  }
  auto bin = fAxisPt->FindBin(pt) - 1;
  if (bin >= 0 && bin < fAxisPt->GetNbins()) {
    is_POI = true;
    if (fUseEtaGap && eta > fEtaGap) is_POI_eta_gap = true;
  }
  if (fCuts.CheckTrackCutsPtCutOnly(track)) {
    if (is_reference)                fRptCut.Fill(phi, weight);
    if      (is_reference_eta_gap_p) fRptCutetaGapP.Fill(phi, weight);
    else if (is_reference_eta_gap_n) fRptCutetaGapN.Fill(phi, weight);
  }
  if (is_reference)                             fR.Fill(phi, weight);
  if (is_POI)                                   fP[bin].Fill(phi, weight);
  if (is_reference && is_POI)                   fQ[bin].Fill(phi, weight);
  if (is_reference_eta_gap_n)                   fRetaGapN.Fill(phi, weight);
  if (is_POI_eta_gap)                           fPetaGapP[bin].Fill(phi, weight);
  if (is_reference_eta_gap_n && is_POI_eta_gap) fQetaGapP[bin].Fill(phi, weight);
  fNUAweightsNtracksTemp->Fill(phi, eta, fVertexZ);
  if (fApplyNUA) fBeforeNUA->Fill(phi, eta, fVertexZ, 1.);
  if (fApplyNUA) fAfterNUA->Fill(phi, eta, fVertexZ, weight);
  FillFilterBitQA(track);
}

void AliZDCcumulantFlow::CalculateCumulants() {
  if (!fCuts.PassedEventCuts()) return;
  auto n22 = TwoParticleRef(fRptCut, 2);
  auto n02 = TwoParticleRef(fRptCut, 0);
  if (n02.real() > 0.) {
    fC22->Fill(fCentrality, (n22 / n02).real(), n02.real());
    fC22Distribution->Fill(fCentrality, (n22 / n02).real(), n02.real());
  }
  auto n22eta = TwoParticleEtaGapRef(fRptCutetaGapN, fRptCutetaGapP, 2);
  auto n02eta = TwoParticleEtaGapRef(fRptCutetaGapN, fRptCutetaGapP, 0);
  if (n02eta.real() > 0.) fC22EtaGap->Fill(fCentrality, (n22eta / n02eta).real(), n02eta.real());
  auto n24 = FourParticleRef(fRptCut,2);
  auto n04 = FourParticleRef(fRptCut,0);
  if (n04.real() > 0.) {
    fC24->Fill(fCentrality, (n24 / n04).real(), n04.real());
    fC24Distribution->Fill(fCentrality, (n24 / n04).real(), n04.real());
    fMultC24->Fill(fCentrality, fRptCut(0,1).real());
  }
  for (int i = 0; i < fAxisPt->GetNbins(); ++i) {
    auto pt = fAxisPt->GetBinCenter(i+1);
    auto nd22 = TwoParticleDif(fQ[i], fP[i], fR, 2);
    auto nd02 = TwoParticleDif(fQ[i], fP[i], fR, 0);
    if (nd02.real() > 0.) fCdif22->Fill(fCentrality, pt, (nd22 / nd02).real(), nd02.real());
    auto nd22eta = TwoParticleDif(fQetaGapP[i], fPetaGapP[i], fRetaGapN, 2);
    auto nd02eta = TwoParticleDif(fQetaGapP[i], fPetaGapP[i], fRetaGapN, 0);
    if (nd02eta.real() > 0.) fCdif22EtaGap->Fill(fCentrality, pt, (nd22eta / nd02eta).real(), nd02eta.real());
    auto nd24 = FourParticleDif(fQ[i], fP[i], fR, 2);
    auto nd04 = FourParticleDif(fQ[i], fP[i], fR, 0);
    if (nd04.real() > 0.) fCdif24->Fill(fCentrality, pt, (nd24 / nd04).real(), nd04.real());
    for (int isample = 0; isample < fNsamples; ++isample) {
      for (int i = 0; i < fSamples[isample]; ++i) {
        if (nd02.real() > 0.)    fCdif22BS[isample]->Fill(fCentrality, pt, (nd22 / nd02).real(), nd02.real());
        if (nd02eta.real() > 0.) fCdif22EtaGapBS[isample]->Fill(fCentrality, pt, (nd22eta / nd02eta).real(), nd02eta.real());
        if (nd04.real() > 0.)    fCdif24BS[isample]->Fill(fCentrality, pt, (nd24 / nd04).real(), nd04.real());
      }
    }
  }
  for (int isample = 0; isample < fNsamples; ++isample) {
    for (int i = 0; i < fSamples[isample]; ++i) {
      if (n02.real() > 0.)    fC22BS[isample]->Fill(fCentrality, (n22 / n02).real(), n02.real());
      if (n02eta.real() > 0.) fC22EtaGapBS[isample]->Fill(fCentrality, (n22eta / n02eta).real(), n02eta.real());
      if (n04.real() > 0.)    fC24BS[isample]->Fill(fCentrality, (n24 / n04).real(), n04.real());
    }
  }
  fR.Reset();
  fRptCut.Reset();
  fRetaGapN.Reset();
  fRptCutetaGapP.Reset();
  fRptCutetaGapN.Reset();
  for (int i = 0; i < fAxisPt->GetNbins(); ++i) {
    fP[i].Reset();
    fQ[i].Reset();
    fPetaGapP[i].Reset();
    fQetaGapP[i].Reset();
  }
  auto x_axis = fNUAweightsNtracksTemp->GetXaxis();
  auto y_axis = fNUAweightsNtracksTemp->GetYaxis();
  auto z_axis = fNUAweightsNtracksTemp->GetZaxis();
  for (auto ix = 1; ix <= x_axis->GetNbins(); ++ix) {
    auto x = x_axis->GetBinCenter(ix);
    for (auto iy = 1; iy <= y_axis->GetNbins(); ++iy) {
      auto y = y_axis->GetBinCenter(iy);
      for (auto iz = 1; iz <= z_axis->GetNbins(); ++iz) {
        auto z = z_axis->GetBinCenter(iz);
        auto ntracks = fNUAweightsNtracksTemp->GetBinContent(ix, iy, iz);
        if (ntracks > 0.) fNUAweightsOut->Fill(x, y, z, ntracks);
      }
    }
  }
  fNUAweightsNtracksTemp->Reset();
}

void AliZDCcumulantFlow::Configure(double eta_gap, std::vector<double> pt_bins, std::vector<double> vtxz_bins, int n_phi_bins, int n_eta_bins, double eta_min, double eta_max) {
  fEtaGap = eta_gap / 2;
  if (fEtaGap > 0.) fUseEtaGap = true;
  auto npt = pt_bins.size()-1;
  fPtMax = pt_bins.back();
  fPtMin = pt_bins.front();
  fAxisPt = new TAxis(npt, pt_bins.data());
  fP.resize(npt);
  fQ.resize(npt);
  fPetaGapP.resize(npt);
  fQetaGapP.resize(npt);
  auto nvtxz = vtxz_bins.size()-1;
  std::vector<double> phi_bins;
  auto phi_bin_width = 2*TMath::Pi() / n_phi_bins;
  for (int i = 0; i < n_phi_bins+1; ++i) phi_bins.push_back(i*phi_bin_width);
  std::vector<double> eta_bins;
  auto eta_bin_width = (eta_max - eta_min) / n_eta_bins;
  for (int i = 0; i < n_phi_bins+1; ++i) eta_bins.push_back(eta_min+i*eta_bin_width);
  fNUAweightsNtracksTemp = new TH3D((fName+"_NUA_n_tracks_temp").data(),";#phi;#eta;vertex Z / cm", n_phi_bins, phi_bins.data(), n_eta_bins, eta_bins.data(), nvtxz, vtxz_bins.data());
  fNUAweightsOut = new TProfile3D((fName+"_NUA").data(),";#phi;#eta;vertex Z / cm", n_phi_bins, phi_bins.data(), n_eta_bins, eta_bins.data(), nvtxz, vtxz_bins.data());
  fBeforeNUA = new TH3D((fName+"_BeforeNUAqa").data(),";#phi;#eta;vertex Z / cm", n_phi_bins, phi_bins.data(), n_eta_bins, eta_bins.data(), nvtxz, vtxz_bins.data());
  fAfterNUA = new TH3D((fName+"_AfterNUAqa").data(),";#phi;#eta;vertex Z / cm", n_phi_bins, phi_bins.data(), n_eta_bins, eta_bins.data(), nvtxz, vtxz_bins.data());
  fAreNUAapplied = new TH1D((fName+"_NUA_flag").data(),";are NUA weights appplied?;",2,0.,2.);
  fFilterBit = new TH1D((fName+"FilterBit").data(), ";filter bit;n tracks;" , 11, 0, 11.);
  fMultC24 = new TH2D((fName+"MultC24").data(),";CentralityV0M;multiplicity",100,0.,100.,5000, 0., 5000.);
  fC22Distribution = new TH2D((fName+"C22Distribution").data(),";CentralityV0M;C_{2}{2}",100,0.,100.,1000, -0.02, 0.08);
  fC24Distribution = new TH2D((fName+"C24Distribution").data(),";CentralityV0M;C_{2}{4}",100,0.,100.,1000,-0.0004,0.00030);
}

TList *AliZDCcumulantFlow::CreateCorrelations() {
  auto correlation_list = new TList();
  correlation_list->SetOwner(true);
  correlation_list->SetName(fName.data());
  auto ptbins = fAxisPt->GetXbins()->GetArray();
  auto nptbins = fAxisPt->GetNbins();
  fC22 = new TProfile("C22",";CentralityV0M;C_{2}{2}",100,0.,100.);
  correlation_list->Add(fC22);
  auto etagap = std::to_string(fEtaGap*2.);
  auto eta_gap_axes = std::string(";CentralityV0M;C_{2}{2,#Delta#eta>"+etagap+"}");
  fC22EtaGap = new TProfile("C22EtaGap",eta_gap_axes.data(),100,0.,100.);
  correlation_list->Add(fC22EtaGap);
  fC24 = new TProfile("C24",";CentralityV0M;C_{2}{4}",100,0.,100.);
  correlation_list->Add(fC24);
  fCdif22 = new TProfile2D("Cdif22",";Centrality V0M;#it{p}_{T} / GeV/#it{c};C'_{2}{2}",100,0.,100.,nptbins, ptbins);
  correlation_list->Add(fCdif22);          
  fCdif24 = new TProfile2D("Cdif24",";Centrality V0M;#it{p}_{T} / GeV/#it{c};C'_{2}{4}",100,0.,100.,nptbins, ptbins);
  correlation_list->Add(fCdif24);
  eta_gap_axes += "#it{p}_{T} / GeV/#it{c}";
  fCdif22EtaGap = new TProfile2D("Cdif22EtaGap",eta_gap_axes.data(),100,0.,100.,nptbins, ptbins);
  correlation_list->Add(fCdif22EtaGap);

  auto bootstrap_list = new TList();
  bootstrap_list->SetOwner(true);
  bootstrap_list->SetName("bs_samples");
  for (int i = 0; i < fNsamples; ++i) {
    auto addnum = [i](const std::string &name){ return name+"_bs_"+std::to_string(i); };
    fC22BS.push_back(new TProfile(addnum("C22").data(),";CentralityV0M;C_{2}{2}",100,0.,100.));
    fC22EtaGapBS.push_back(new TProfile(addnum("C22EtaGap").data(),eta_gap_axes.data(),100,0.,100.));
    fC24BS.push_back(new TProfile(addnum("C24").data(),";CentralityV0M;C_{2}{4}",100,0.,100.));
    fCdif22BS.push_back(new TProfile2D(addnum("Cdif22").data(),";Centrality V0M;#it{p}_{T} / GeV/#it{c};C'_{2}{2}",100,0.,100.,nptbins, ptbins));
    fCdif24BS.push_back(new TProfile2D(addnum("Cdif24").data(),";Centrality V0M;#it{p}_{T} / GeV/#it{c};C'_{2}{4}",100,0.,100.,nptbins, ptbins));
    fCdif22EtaGapBS.push_back(new TProfile2D(addnum("Cdif22EtaGap").data(),eta_gap_axes.data(),100,0.,100.,nptbins, ptbins));
    bootstrap_list->Add(fC22BS.back());
    bootstrap_list->Add(fC22EtaGapBS.back());
    bootstrap_list->Add(fC24BS.back());
    bootstrap_list->Add(fCdif22BS.back());
    bootstrap_list->Add(fCdif24BS.back());
    bootstrap_list->Add(fCdif22EtaGapBS.back());
  }
  correlation_list->Add(bootstrap_list);
  return correlation_list;
}

void AliZDCcumulantFlow::AddCorrectionsToList(TList* hlist, TList *qalist) {
  if(fApplyNUA) hlist->Add(fNUAweightsOut);
  auto nua_qa_list = new TList();
  nua_qa_list->SetName((fName+"_QA").data());
  nua_qa_list->SetOwner(true);
  nua_qa_list->Add(fAreNUAapplied);
  nua_qa_list->Add(fBeforeNUA);
  nua_qa_list->Add(fAfterNUA);
  nua_qa_list->Add(fFilterBit);
  nua_qa_list->Add(fMultC24);
  nua_qa_list->Add(fC22Distribution);
  nua_qa_list->Add(fC24Distribution);
  fCuts.AddQAHistograms(nua_qa_list);
  qalist->Add(nua_qa_list);
}

TProfile3D* AliZDCcumulantFlow::ReadFromOADB(TFile *file, const std::string &oadb_name, Int_t run_number) {
  auto check_histo = [](TH1* histo){ 
    for (int i = 0; i < histo->GetNbinsX(); ++i) {
      if (std::isnan(histo->GetBinContent(i)) || std::isnan(histo->GetBinError(i))) return false;
    }
    return true;
  };
  TProfile3D *ret = nullptr;
  fIsApplied = false;
  auto oadb = dynamic_cast<AliOADBContainer*>(file->Get(oadb_name.data()));
  if (oadb) {
    auto profile = dynamic_cast<TProfile3D*>(oadb->GetObject(run_number));
    if (profile && profile->GetEntries() > 0. && check_histo(profile)) {
      ret = profile;
      fIsApplied = true;
    }
  }
  return ret;
}

void AliZDCcumulantFlow::OpenCorrection(TFile *file, Int_t run_number) {
  if (fApplyNUA) {
    if (!file || file->IsZombie()) return;
    fNUAweightsIn = ReadFromOADB(file, fNUAweightsOut->GetName(), run_number);
    if (fNUAweightsIn) fIsApplied = true;
    if (fIsApplied) fAreNUAapplied->Fill(1.);
    else            fAreNUAapplied->Fill(0.);
  }
}
