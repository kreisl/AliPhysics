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

#include <cmath>
#include "TProfile.h"
#include "TProfile2D.h"
#include "TH1D.h"
#include "TList.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliQvector.h"
#include "AliZDCdirectedFlow.h"

void AliZDCdirectedFlow::FindCentralityBin(AliAODEvent *event, const Double_t centrality, const std::vector<Int_t> &samples) {
  fCuts.CheckEventCuts(event); 
  fIsInCentralityClassWide[0] = centrality > 10. && centrality < 20.;
  fIsInCentralityClassWide[1] = centrality > 30. && centrality < 40.;
  fIsInCentralityClassWide[2] = centrality > 10. && centrality < 60.;
  fCentrality = centrality;
  fCuts.HasPassedCentrality(centrality > 0. && centrality < 90.);
  fSamples = samples;
}

void AliZDCdirectedFlow::FillPerEventCorrelations() {
  if (fCuts.PassedEventCuts()) {
    fXaXcCent->Fill(fCentrality, fQZA.x*fQZC.x);
    fYaYcCent->Fill(fCentrality, fQZA.y*fQZC.y);
    fXaYcCent->Fill(fCentrality, fQZA.x*fQZC.y);
    fYaXcCent->Fill(fCentrality, fQZA.y*fQZC.x);
  }
}

void AliZDCdirectedFlow::FillPerTrackCorrelations(AliAODTrack *track) {
  if (!(fCuts.PassedEventCuts() && fCuts.CheckTrackCutsNoPtCut(track))) return;
  const Double_t phi = track->Phi();
  const Double_t eta = track->Eta();
  const Double_t pt = track->Pt(); 

  const auto v1_xxa = std::cos(phi)*fQZA.x;
  const auto v1_yya = std::sin(phi)*fQZA.y;
  const auto v1_yxa = std::sin(phi)*fQZA.x;
  const auto v1_xya = std::cos(phi)*fQZA.y;

  const auto v1_xxc = std::cos(phi)*fQZC.x;
  const auto v1_yyc = std::sin(phi)*fQZC.y;
  const auto v1_yxc = std::sin(phi)*fQZC.x;
  const auto v1_xyc = std::cos(phi)*fQZC.y;
  // wide centrality bins
  const auto n_bins = fIsInCentralityClassWide.size();
  for (int isample = 0; isample < fNsamples; ++isample) {
    for (int i = 0; i < fSamples[isample]; ++i) {
      if (eta > 0.) {
          fV1XXaPtEtaNCentBS[isample]->Fill(fCentrality, pt, v1_xxa);
          fV1XXcPtEtaNCentBS[isample]->Fill(fCentrality, pt, v1_xxc);
          fV1YYaPtEtaNCentBS[isample]->Fill(fCentrality, pt, v1_yya);
          fV1YYcPtEtaNCentBS[isample]->Fill(fCentrality, pt, v1_yyc);
      } else {
          fV1XXaPtEtaPCentBS[isample]->Fill(fCentrality, pt, v1_xxa);
          fV1XXcPtEtaPCentBS[isample]->Fill(fCentrality, pt, v1_xxc);
          fV1YYaPtEtaPCentBS[isample]->Fill(fCentrality, pt, v1_xxa);
          fV1YYcPtEtaPCentBS[isample]->Fill(fCentrality, pt, v1_xxc);
      }
    }
  }
  if (eta > 0.) {
      fV1XXaPtEtaNCent->Fill(fCentrality, pt, v1_xxa);
      fV1XXcPtEtaNCent->Fill(fCentrality, pt, v1_xxc);
      fV1YYaPtEtaNCent->Fill(fCentrality, pt, v1_yya);
      fV1YYcPtEtaNCent->Fill(fCentrality, pt, v1_yyc);
  } else {
      fV1XXaPtEtaPCent->Fill(fCentrality, pt, v1_xxa);
      fV1XXcPtEtaPCent->Fill(fCentrality, pt, v1_xxc);
      fV1YYaPtEtaPCent->Fill(fCentrality, pt, v1_xxa);
      fV1YYcPtEtaPCent->Fill(fCentrality, pt, v1_xxc);
  }
  for (auto i = 0u; i < n_bins; ++i) {
    if (!fIsInCentralityClassWide[i]) continue;
    if (eta > 0.) {
      fV1XXaPtEtaP[i]->Fill(pt,v1_xxa); 
      fV1XXcPtEtaP[i]->Fill(pt,v1_xxc); 
      fV1YYaPtEtaP[i]->Fill(pt,v1_yya); 
      fV1YYcPtEtaP[i]->Fill(pt,v1_yya); 
    } else {
      fV1XXaPtEtaN[i]->Fill(pt,v1_xxa); 
      fV1XXcPtEtaN[i]->Fill(pt,v1_xxc); 
      fV1YYaPtEtaN[i]->Fill(pt,v1_yya); 
      fV1YYcPtEtaN[i]->Fill(pt,v1_yya); 
    }
  }
  if (!fCuts.CheckTrackCutsPtCutOnly(track)) return;
  for (auto i = 0u; i < n_bins; ++i) {
    if (!fIsInCentralityClassWide[i]) continue;
    fV1XXaEta[i]->Fill(eta, v1_xxa);
    fV1XXcEta[i]->Fill(eta, v1_xxc);
    fV1YYaEta[i]->Fill(eta, v1_yya);
    fV1YYcEta[i]->Fill(eta, v1_yyc);
    fV1XYaEta[i]->Fill(eta, v1_xya);
    fV1XYcEta[i]->Fill(eta, v1_xyc);
    fV1YXaEta[i]->Fill(eta, v1_yxa);
    fV1YXcEta[i]->Fill(eta, v1_yxc);
    fPt[i]->Fill(pt);
  }
  fV1XXaEtaCent->Fill(fCentrality, eta, v1_xxa);
  fV1XXcEtaCent->Fill(fCentrality, eta, v1_xxc);
  fV1YYaEtaCent->Fill(fCentrality, eta, v1_yya);
  fV1YYcEtaCent->Fill(fCentrality, eta, v1_yyc);
  fV1XYaEtaCent->Fill(fCentrality, eta, v1_xya);
  fV1XYcEtaCent->Fill(fCentrality, eta, v1_xyc);
  fV1YXaEtaCent->Fill(fCentrality, eta, v1_yxa);
  fV1YXcEtaCent->Fill(fCentrality, eta, v1_yxc);

  for (int isample = 0; isample < fNsamples; ++isample) {
    for (int i = 0; i < fSamples[isample]; ++i) {
      fV1XXaEtaCentBS[isample]->Fill(fCentrality, eta, v1_xxa);
      fV1XXcEtaCentBS[isample]->Fill(fCentrality, eta, v1_xxc);
      fV1YYaEtaCentBS[isample]->Fill(fCentrality, eta, v1_yya);
      fV1YYcEtaCentBS[isample]->Fill(fCentrality, eta, v1_yyc);
    }
  }

  // 1 percent centrality bins
  if (eta > 0.) {
    fXtXaCentEtaP->Fill(fCentrality,v1_xxa);
    fXtXcCentEtaP->Fill(fCentrality,v1_xxc);
    fYtYaCentEtaP->Fill(fCentrality,v1_yya);
    fYtYcCentEtaP->Fill(fCentrality,v1_yya);
  } else {
    fXtXaCentEtaN->Fill(fCentrality,v1_xxa);
    fXtXcCentEtaN->Fill(fCentrality,v1_xxc);
    fYtYaCentEtaN->Fill(fCentrality,v1_yya);
    fYtYcCentEtaN->Fill(fCentrality,v1_yya);
  }
}

TList *AliZDCdirectedFlow::CreateCorrelations() {
  auto correlation_list = new TList();
  correlation_list->SetOwner(true);
  correlation_list->SetName(fName.data());
  auto cent_bins_list = new TList();
  cent_bins_list->SetOwner(true);
  cent_bins_list->SetName("centrality_bins");
  for (unsigned int i = 0; i < fNcentBinsWide; ++i) {
    const auto ibin = std::to_string(i);
    fV1XXaEta[i] = new TProfile((std::string("vnXXaEta_")+ibin).data(),";#eta;#LTx_{1}X_{1}^{ZNA}#GT",fNetaBins,fEtaBins); 
    fV1XXcEta[i] = new TProfile((std::string("vnXXcEta_")+ibin).data(),";#eta;#LTx_{1}X_{1}^{ZNC}#GT",fNetaBins,fEtaBins); 
    fV1YYaEta[i] = new TProfile((std::string("vnYYaEta_")+ibin).data(),";#eta;#LTy_{1}Y_{1}^{ZNA}#GT",fNetaBins,fEtaBins); 
    fV1YYcEta[i] = new TProfile((std::string("vnYYcEta_")+ibin).data(),";#eta;#LTy_{1}Y_{1}^{ZNC}#GT",fNetaBins,fEtaBins); 
    fV1XYaEta[i] = new TProfile((std::string("vnXYaEta_")+ibin).data(),";#eta;#LTx_{1}Y_{1}^{ZNA}#GT",fNetaBins,fEtaBins); 
    fV1XYcEta[i] = new TProfile((std::string("vnXYcEta_")+ibin).data(),";#eta;#LTx_{1}Y_{1}^{ZNC}#GT",fNetaBins,fEtaBins); 
    fV1YXaEta[i] = new TProfile((std::string("vnYXaEta_")+ibin).data(),";#eta;#LTy_{1}X_{1}^{ZNA}#GT",fNetaBins,fEtaBins); 
    fV1YXcEta[i] = new TProfile((std::string("vnYXcEta_")+ibin).data(),";#eta;#LTy_{1}X_{1}^{ZNC}#GT",fNetaBins,fEtaBins); 
    cent_bins_list->Add(fV1XXaEta[i]); 
    cent_bins_list->Add(fV1XXcEta[i]); 
    cent_bins_list->Add(fV1YYaEta[i]); 
    cent_bins_list->Add(fV1YYcEta[i]); 
    cent_bins_list->Add(fV1XYaEta[i]); 
    cent_bins_list->Add(fV1XYcEta[i]); 
    cent_bins_list->Add(fV1YXaEta[i]); 
    cent_bins_list->Add(fV1YXcEta[i]); 
    fV1XXaPtEtaP[i] = new TProfile((std::string("v1XXaPtEtaP")+ibin).data(), ";#it{p}_{T} / GeV/#it{c};#LTx_{1}X_{1}^{ZNA}#GT, #eta > 0", fNptBinsWide, fPtBinsWide); 
    fV1XXcPtEtaP[i] = new TProfile((std::string("v1XXcPtEtaP")+ibin).data(), ";#it{p}_{T} / GeV/#it{c};#LTx_{1}X_{1}^{ZNC}#GT, #eta > 0", fNptBinsWide, fPtBinsWide); 
    fV1YYaPtEtaP[i] = new TProfile((std::string("v1YYaPtEtaP")+ibin).data(), ";#it{p}_{T} / GeV/#it{c};#LTy_{1}Y_{1}^{ZNA}#GT, #eta > 0", fNptBinsWide, fPtBinsWide); 
    fV1YYcPtEtaP[i] = new TProfile((std::string("v1YYcPtEtaP")+ibin).data(), ";#it{p}_{T} / GeV/#it{c};#LTy_{1}Y_{1}^{ZNC}#GT, #eta > 0", fNptBinsWide, fPtBinsWide); 
    cent_bins_list->Add(fV1XXaPtEtaP[i]); 
    cent_bins_list->Add(fV1XXcPtEtaP[i]); 
    cent_bins_list->Add(fV1YYaPtEtaP[i]); 
    cent_bins_list->Add(fV1YYcPtEtaP[i]); 
    fV1XXaPtEtaN[i] = new TProfile((std::string("v1XXaPtEtaN")+ibin).data(), ";#it{p}_{T} / GeV/#it{c};#LTx_{1}X_{1}^{ZNA}#GT, #eta < 0", fNptBinsWide, fPtBinsWide); 
    fV1XXcPtEtaN[i] = new TProfile((std::string("v1XXcPtEtaN")+ibin).data(), ";#it{p}_{T} / GeV/#it{c};#LTx_{1}X_{1}^{ZNC}#GT, #eta < 0", fNptBinsWide, fPtBinsWide); 
    fV1YYaPtEtaN[i] = new TProfile((std::string("v1YYaPtEtaN")+ibin).data(), ";#it{p}_{T} / GeV/#it{c};#LTy_{1}Y_{1}^{ZNA}#GT, #eta < 0", fNptBinsWide, fPtBinsWide); 
    fV1YYcPtEtaN[i] = new TProfile((std::string("v1YYcPtEtaN")+ibin).data(), ";#it{p}_{T} / GeV/#it{c};#LTy_{1}Y_{1}^{ZNC}#GT, #eta < 0", fNptBinsWide, fPtBinsWide); 
    cent_bins_list->Add(fV1XXaPtEtaN[i]); 
    cent_bins_list->Add(fV1XXcPtEtaN[i]); 
    cent_bins_list->Add(fV1YYaPtEtaN[i]); 
    cent_bins_list->Add(fV1YYcPtEtaN[i]); 
    fPt[i] = new TH1D(std::string("Pt"+ibin).data(),";#it{p}_{T} / GeV/#it{c};N",fNptBinsWide, fPtBinsWide);
    cent_bins_list->Add(fPt[i]);
  }
  correlation_list->Add(cent_bins_list);
  fXaXcCent = new TProfile("XXCent",";centrality V0M;#LTX_{1}^{ZNA}X_{1}^{ZNC}#GT",100,0.,100.);
  fYaYcCent = new TProfile("YYCent",";centrality V0M;#LTY_{1}^{ZNA}Y_{1}^{ZNC}#GT",100,0.,100.);
  fXaYcCent = new TProfile("XYCent",";centrality V0M;#LTX_{1}^{ZNA}Y_{1}^{ZNC}#GT",100,0.,100.);
  fYaXcCent = new TProfile("YXCent",";centrality V0M;#LTY_{1}^{ZNA}X_{1}^{ZNC}#GT",100,0.,100.);
  correlation_list->Add(fXaXcCent);
  correlation_list->Add(fYaYcCent);
  correlation_list->Add(fXaYcCent);
  correlation_list->Add(fYaXcCent);
  fXtXaCentEtaP = new TProfile("xxaCentP", ";centrality V0M;#LTx_{1}X_{1}^{ZNA}#GT, #eta > 0", 100, 0., 100.);
  fXtXcCentEtaP = new TProfile("xxcCentP", ";centrality V0M;#LTx_{1}X_{1}^{ZNC}#GT, #eta > 0", 100, 0., 100.);
  fYtYaCentEtaP = new TProfile("yyaCentP", ";centrality V0M;#LTy_{1}Y_{1}^{ZNA}#GT, #eta > 0", 100, 0., 100.);
  fYtYcCentEtaP = new TProfile("yycCentP", ";centrality V0M;#LTy_{1}Y_{1}^{ZNC}#GT, #eta > 0", 100, 0., 100.);
  correlation_list->Add(fXtXaCentEtaP);
  correlation_list->Add(fXtXcCentEtaP);
  correlation_list->Add(fYtYaCentEtaP);
  correlation_list->Add(fYtYcCentEtaP);
  fXtXaCentEtaN = new TProfile("xxaCentN", ";centrality V0M;#LTx_{1}X_{1}^{ZNA}#GT, #eta < 0", 100, 0., 100.);
  fXtXcCentEtaN = new TProfile("xxcCentN", ";centrality V0M;#LTx_{1}X_{1}^{ZNC}#GT, #eta < 0", 100, 0., 100.);
  fYtYaCentEtaN = new TProfile("yyaCentN", ";centrality V0M;#LTy_{1}Y_{1}^{ZNA}#GT, #eta < 0", 100, 0., 100.);
  fYtYcCentEtaN = new TProfile("yycCentN", ";centrality V0M;#LTy_{1}Y_{1}^{ZNC}#GT, #eta < 0", 100, 0., 100.);
  correlation_list->Add(fXtXaCentEtaN);
  correlation_list->Add(fXtXcCentEtaN);
  correlation_list->Add(fYtYaCentEtaN);
  correlation_list->Add(fYtYcCentEtaN);

  fV1XXaEtaCent = new TProfile2D("xxaEtaCent", ";centrality V0M;#eta;#LTx_{1}X_{1}^{ZNA}#GT", 100, 0., 100., fNetaBins, fEtaBins); 
  fV1XXcEtaCent = new TProfile2D("xxcEtaCent", ";centrality V0M;#eta;#LTx_{1}X_{1}^{ZNC}#GT", 100, 0., 100., fNetaBins, fEtaBins); 
  fV1YYaEtaCent = new TProfile2D("yyaEtaCent", ";centrality V0M;#eta;#LTy_{1}Y_{1}^{ZNA}#GT", 100, 0., 100., fNetaBins, fEtaBins); 
  fV1YYcEtaCent = new TProfile2D("yycEtaCent", ";centrality V0M;#eta;#LTy_{1}Y_{1}^{ZNC}#GT", 100, 0., 100., fNetaBins, fEtaBins); 
  fV1XYaEtaCent = new TProfile2D("xyaEtaCent", ";centrality V0M;#eta;#LTx_{1}X_{1}^{ZNA}#GT", 100, 0., 100., fNetaBins, fEtaBins); 
  fV1XYcEtaCent = new TProfile2D("xycEtaCent", ";centrality V0M;#eta;#LTx_{1}X_{1}^{ZNC}#GT", 100, 0., 100., fNetaBins, fEtaBins); 
  fV1YXaEtaCent = new TProfile2D("yxaEtaCent", ";centrality V0M;#eta;#LTy_{1}Y_{1}^{ZNA}#GT", 100, 0., 100., fNetaBins, fEtaBins); 
  fV1YXcEtaCent = new TProfile2D("yxcEtaCent", ";centrality V0M;#eta;#LTy_{1}Y_{1}^{ZNC}#GT", 100, 0., 100., fNetaBins, fEtaBins); 
  correlation_list->Add(fV1XXaEtaCent);
  correlation_list->Add(fV1XXcEtaCent);
  correlation_list->Add(fV1YYaEtaCent);
  correlation_list->Add(fV1YYcEtaCent);
  correlation_list->Add(fV1XYaEtaCent);
  correlation_list->Add(fV1XYcEtaCent);
  correlation_list->Add(fV1YXaEtaCent);
  correlation_list->Add(fV1YXcEtaCent);

  auto bs_list = new TList();
  bs_list->SetOwner(true);
  bs_list->SetName("bs_samples");
  for (int i = 0; i < fNsamples; ++i) {
    auto addnum = [i](const std::string &name){ return name+"_bs_"+std::to_string(i); };
    fV1XXaEtaCentBS.push_back(new TProfile2D(addnum("xxaEtaCent").data(), ";centrality V0M;#eta;#LTx_{1}X_{1}^{ZNA}#GT", 100, 0., 100., fNetaBins, fEtaBins)); 
    fV1XXcEtaCentBS.push_back(new TProfile2D(addnum("xxcEtaCent").data(), ";centrality V0M;#eta;#LTx_{1}X_{1}^{ZNC}#GT", 100, 0., 100., fNetaBins, fEtaBins)); 
    fV1YYaEtaCentBS.push_back(new TProfile2D(addnum("yyaEtaCent").data(), ";centrality V0M;#eta;#LTy_{1}Y_{1}^{ZNA}#GT", 100, 0., 100., fNetaBins, fEtaBins)); 
    fV1YYcEtaCentBS.push_back(new TProfile2D(addnum("yycEtaCent").data(), ";centrality V0M;#eta;#LTy_{1}Y_{1}^{ZNC}#GT", 100, 0., 100., fNetaBins, fEtaBins)); 
    bs_list->Add(fV1XXaEtaCentBS.back());
    bs_list->Add(fV1XXcEtaCentBS.back());
    bs_list->Add(fV1YYaEtaCentBS.back());
    bs_list->Add(fV1YYcEtaCentBS.back());
    fV1XXaPtEtaPCentBS.push_back(new TProfile2D("xxaPtCentP", ";centrality V0M #it{p}_{T} / GeV/#it{c};#LTx_{1}X_{1}^{ZNA}#GT, #eta > 0", 100, 0., 100., fNptBinsWide, fPtBinsWide));
    fV1XXcPtEtaPCentBS.push_back(new TProfile2D("xxcPtCentP", ";centrality V0M #it{p}_{T} / GeV/#it{c};#LTx_{1}X_{1}^{ZNC}#GT, #eta > 0", 100, 0., 100., fNptBinsWide, fPtBinsWide));
    fV1YYaPtEtaPCentBS.push_back(new TProfile2D("yyaPtCentP", ";centrality V0M #it{p}_{T} / GeV/#it{c};#LTy_{1}Y_{1}^{ZNA}#GT, #eta > 0", 100, 0., 100., fNptBinsWide, fPtBinsWide));
    fV1YYcPtEtaPCentBS.push_back(new TProfile2D("yycPtCentP", ";centrality V0M #it{p}_{T} / GeV/#it{c};#LTy_{1}Y_{1}^{ZNC}#GT, #eta > 0", 100, 0., 100., fNptBinsWide, fPtBinsWide));
    fV1XXaPtEtaNCentBS.push_back(new TProfile2D("xxaPtCentN", ";centrality V0M #it{p}_{T} / GeV/#it{c};#LTx_{1}X_{1}^{ZNA}#GT, #eta < 0", 100, 0., 100., fNptBinsWide, fPtBinsWide));
    fV1XXcPtEtaNCentBS.push_back(new TProfile2D("xxcPtCentN", ";centrality V0M #it{p}_{T} / GeV/#it{c};#LTx_{1}X_{1}^{ZNC}#GT, #eta < 0", 100, 0., 100., fNptBinsWide, fPtBinsWide));
    fV1YYaPtEtaNCentBS.push_back(new TProfile2D("yyaPtCentN", ";centrality V0M #it{p}_{T} / GeV/#it{c};#LTy_{1}Y_{1}^{ZNA}#GT, #eta < 0", 100, 0., 100., fNptBinsWide, fPtBinsWide));
    fV1YYcPtEtaNCentBS.push_back(new TProfile2D("yycPtCentN", ";centrality V0M #it{p}_{T} / GeV/#it{c};#LTy_{1}Y_{1}^{ZNC}#GT, #eta < 0", 100, 0., 100., fNptBinsWide, fPtBinsWide));
    bs_list->Add(fV1XXaPtEtaPCentBS.back());
    bs_list->Add(fV1XXcPtEtaPCentBS.back());
    bs_list->Add(fV1YYaPtEtaPCentBS.back());
    bs_list->Add(fV1YYcPtEtaPCentBS.back());
    bs_list->Add(fV1XXaPtEtaNCentBS.back());
    bs_list->Add(fV1XXcPtEtaNCentBS.back());
    bs_list->Add(fV1YYaPtEtaNCentBS.back());
    bs_list->Add(fV1YYcPtEtaNCentBS.back());
  }
  correlation_list->Add(bs_list);

  fV1XXaPtEtaPCent = new TProfile2D("xxaPtCentP", ";centrality V0M #it{p}_{T} / GeV/#it{c};#LTx_{1}X_{1}^{ZNA}#GT, #eta > 0", 100, 0., 100., fNptBinsWide, fPtBinsWide);
  fV1XXcPtEtaPCent = new TProfile2D("xxcPtCentP", ";centrality V0M #it{p}_{T} / GeV/#it{c};#LTx_{1}X_{1}^{ZNC}#GT, #eta > 0", 100, 0., 100., fNptBinsWide, fPtBinsWide);
  fV1YYaPtEtaPCent = new TProfile2D("yyaPtCentP", ";centrality V0M #it{p}_{T} / GeV/#it{c};#LTy_{1}Y_{1}^{ZNA}#GT, #eta > 0", 100, 0., 100., fNptBinsWide, fPtBinsWide);
  fV1YYcPtEtaPCent = new TProfile2D("yycPtCentP", ";centrality V0M #it{p}_{T} / GeV/#it{c};#LTy_{1}Y_{1}^{ZNC}#GT, #eta > 0", 100, 0., 100., fNptBinsWide, fPtBinsWide);
  fV1XXaPtEtaNCent = new TProfile2D("xxaPtCentN", ";centrality V0M #it{p}_{T} / GeV/#it{c};#LTx_{1}X_{1}^{ZNA}#GT, #eta < 0", 100, 0., 100., fNptBinsWide, fPtBinsWide);
  fV1XXcPtEtaNCent = new TProfile2D("xxcPtCentN", ";centrality V0M #it{p}_{T} / GeV/#it{c};#LTx_{1}X_{1}^{ZNC}#GT, #eta < 0", 100, 0., 100., fNptBinsWide, fPtBinsWide);
  fV1YYaPtEtaNCent = new TProfile2D("yyaPtCentN", ";centrality V0M #it{p}_{T} / GeV/#it{c};#LTy_{1}Y_{1}^{ZNA}#GT, #eta < 0", 100, 0., 100., fNptBinsWide, fPtBinsWide);
  fV1YYcPtEtaNCent = new TProfile2D("yycPtCentN", ";centrality V0M #it{p}_{T} / GeV/#it{c};#LTy_{1}Y_{1}^{ZNC}#GT, #eta < 0", 100, 0., 100., fNptBinsWide, fPtBinsWide);
  correlation_list->Add(fV1XXaPtEtaPCent);
  correlation_list->Add(fV1XXcPtEtaPCent);
  correlation_list->Add(fV1YYaPtEtaPCent);
  correlation_list->Add(fV1YYcPtEtaPCent);
  correlation_list->Add(fV1XXaPtEtaNCent);
  correlation_list->Add(fV1XXcPtEtaNCent);
  correlation_list->Add(fV1YYaPtEtaNCent);
  correlation_list->Add(fV1YYcPtEtaNCent);
  return correlation_list;
}
