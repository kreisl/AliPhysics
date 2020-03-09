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
#include <iostream>
#include "TH1.h"
#include "TList.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliQvector.h"
#include "AliZDCellipticFlow.h"

void AliZDCellipticFlow::FindCentralityBin(AliAODEvent *event, Double_t centrality, const std::vector<Int_t> &samples) {
  fCuts.CheckEventCuts(event);
  fCentralityBin = fCentralityAxis->FindBin(centrality) - 1;
  fCentrality = centrality;
  fCuts.HasPassedCentrality(centrality > 0. && centrality < 90.);
  fSamples = samples;
}

void AliZDCellipticFlow::FillPerEventCorrelations() {
  if (fCuts.PassedEventCuts()) {
    fXaXcCent->Fill(fCentrality, fQZA.x*fQZC.x);
    fYaYcCent->Fill(fCentrality, fQZA.y*fQZC.y);
    fXaYcCent->Fill(fCentrality, fQZA.x*fQZC.y);
    fYaXcCent->Fill(fCentrality, fQZA.y*fQZC.x);
    fXaXcDistCent->Fill(fCentrality,fQZA.x*fQZC.x);
    fYaYcDistCent->Fill(fCentrality,fQZA.y*fQZC.y);
    fXaYcDistCent->Fill(fCentrality,fQZA.x*fQZC.y);
    fYaXcDistCent->Fill(fCentrality,fQZA.y*fQZC.x);
    for (int isample = 0; isample < fNsamples; ++isample) {
      for (int i = 0; i < fSamples[isample]; ++i) {
          fXaXcCentBS[isample]->Fill(fCentrality, fQZA.x*fQZC.x);
          fYaYcCentBS[isample]->Fill(fCentrality, fQZA.y*fQZC.y);
          fXaYcCentBS[isample]->Fill(fCentrality, fQZA.x*fQZC.y);
          fYaXcCentBS[isample]->Fill(fCentrality, fQZA.y*fQZC.x);
      }
    }
  }
}

void AliZDCellipticFlow::FillPerTrackCorrelations(AliAODTrack *track) {
  if (!(fCuts.PassedEventCuts() && fCuts.CheckTrackCutsNoPtCut(track))) return;
  const Double_t phi = track->Phi();
  const Double_t pt = track->Pt(); 
  double v2_xxx = std::cos(2.*phi)*fQZA.x*fQZC.x;
  double v2_xyy = std::cos(2.*phi)*fQZA.y*fQZC.y;
  double v2_yyx = std::sin(2.*phi)*fQZA.y*fQZC.x;
  double v2_yxy = std::sin(2.*phi)*fQZA.x*fQZC.y;
  double v2_yyy = std::sin(2.*phi)*fQZA.y*fQZC.y;
  double v2_yxx = std::sin(2.*phi)*fQZA.x*fQZC.x;
  double v2_xyx = std::cos(2.*phi)*fQZA.y*fQZC.x;
  double v2_xxy = std::cos(2.*phi)*fQZA.x*fQZC.y;

  fV2XXXpTcent->Fill(fCentrality, pt, v2_xxx);
  fV2XYYpTcent->Fill(fCentrality, pt, v2_xyy);
  fV2YXYpTcent->Fill(fCentrality, pt, v2_yxy);
  fV2YYXpTcent->Fill(fCentrality, pt, v2_yyx);
  fV2YYYpTcent->Fill(fCentrality, pt, v2_yyy);
  fV2YXXpTcent->Fill(fCentrality, pt, v2_yxx);
  fV2XXYpTcent->Fill(fCentrality, pt, v2_xxy);
  fV2XYXpTcent->Fill(fCentrality, pt, v2_xyx);

  for (int isample = 0; isample < fNsamples; ++isample) {
    for (int i = 0; i < fSamples[isample]; ++i) {
      fV2XXXpTcentBS[isample]->Fill(fCentrality, pt, v2_xxx);
      fV2XYYpTcentBS[isample]->Fill(fCentrality, pt, v2_xyy);
      fV2YXYpTcentBS[isample]->Fill(fCentrality, pt, v2_yxy);
      fV2YYXpTcentBS[isample]->Fill(fCentrality, pt, v2_yyx);
    }
  }

  fV2XXXpT[fCentralityBin]->Fill(pt, v2_xxx);
  fV2XYYpT[fCentralityBin]->Fill(pt, v2_xyy);
  fV2YXYpT[fCentralityBin]->Fill(pt, v2_yxy);
  fV2YYXpT[fCentralityBin]->Fill(pt, v2_yyx);
  fV2YYYpT[fCentralityBin]->Fill(pt, v2_yyy);
  fV2YXXpT[fCentralityBin]->Fill(pt, v2_yxx);
  fV2XXYpT[fCentralityBin]->Fill(pt, v2_xxy);
  fV2XYXpT[fCentralityBin]->Fill(pt, v2_xyx);
  fPt[fCentralityBin]->Fill(pt);
  if (!fCuts.CheckTrackCutsPtCutOnly(track)) return;
  fXtXaXcDistCent->Fill(fCentrality, v2_xxx);
  fXtYaYcDistCent->Fill(fCentrality, v2_xyy);
  fYtXaYcDistCent->Fill(fCentrality, v2_yxy);
  fYtYaXcDistCent->Fill(fCentrality, v2_yyx);

  fXtXaXcCent->Fill(fCentrality, v2_xxx);
  fXtXaYcCent->Fill(fCentrality, v2_xxy);
  fXtYaXcCent->Fill(fCentrality, v2_xyx);
  fYtXaXcCent->Fill(fCentrality, v2_yxx);
  fYtYaYcCent->Fill(fCentrality, v2_yyy);
  fYtYaXcCent->Fill(fCentrality, v2_yyx);
  fXtYaYcCent->Fill(fCentrality, v2_xyy);
  fYtXaYcCent->Fill(fCentrality, v2_yxy);
}

TList *AliZDCellipticFlow::CreateCorrelations() {
  auto nptbins = fPtBins.size() - 1;
  auto ptbins = fPtBins.data();
  fCentralityAxis = new TAxis(fNcentBins, centBins);
  auto correlation_list = new TList();
  correlation_list->SetOwner(true);
  correlation_list->SetName(fName.data());
  auto cent_bins_list = new TList();
  cent_bins_list->SetOwner(true);
  cent_bins_list->SetName("centrality_bins");
  for (unsigned int i = 0; i < fNcentBins; ++i) {
    auto ibin = std::to_string(i);
    fV2XXXpT[i] = new TProfile((std::string("vnXXXpt_")+ibin).data(),";#it{p}_{T} / GeV/#it{c};#LTx_{2}X_{1}^{ZNA}X_{1}^{ZNC}#GT",nptbins,ptbins); 
    fV2XYYpT[i] = new TProfile((std::string("vnXYYpt_")+ibin).data(),";#it{p}_{T} / GeV/#it{c};#LTx_{2}Y_{1}^{ZNA}Y_{1}^{ZNC}#GT",nptbins,ptbins);
    fV2YXYpT[i] = new TProfile((std::string("vnYXYpt_")+ibin).data(),";#it{p}_{T} / GeV/#it{c};#LTy_{2}X_{1}^{ZNA}Y_{1}^{ZNC}#GT",nptbins,ptbins);
    fV2YYXpT[i] = new TProfile((std::string("vnYYXpt_")+ibin).data(),";#it{p}_{T} / GeV/#it{c};#LTy_{2}Y_{1}^{ZNA}X_{1}^{ZNC}#GT",nptbins,ptbins);
    fV2YYYpT[i] = new TProfile((std::string("vnYYYpt_")+ibin).data(),";#it{p}_{T} / GeV/#it{c};#LTy_{2}Y_{1}^{ZNA}Y_{1}^{ZNC}#GT",nptbins,ptbins);
    fV2YXXpT[i] = new TProfile((std::string("vnYXXpt_")+ibin).data(),";#it{p}_{T} / GeV/#it{c};#LTy_{2}X_{1}^{ZNA}X_{1}^{ZNC}#GT",nptbins,ptbins);
    fV2XXYpT[i] = new TProfile((std::string("vnXXYpt_")+ibin).data(),";#it{p}_{T} / GeV/#it{c};#LTx_{2}X_{1}^{ZNA}Y_{1}^{ZNC}#GT",nptbins,ptbins);
    fV2XYXpT[i] = new TProfile((std::string("vnXYXpt_")+ibin).data(),";#it{p}_{T} / GeV/#it{c};#LTx_{2}Y_{1}^{ZNA}X_{1}^{ZNC}#GT",nptbins,ptbins);
    cent_bins_list->Add(fV2XXXpT[i]);
    cent_bins_list->Add(fV2XYYpT[i]);
    cent_bins_list->Add(fV2YXYpT[i]);
    cent_bins_list->Add(fV2YYXpT[i]);
    cent_bins_list->Add(fV2YYYpT[i]);
    cent_bins_list->Add(fV2YXXpT[i]);
    cent_bins_list->Add(fV2XXYpT[i]);
    cent_bins_list->Add(fV2XYXpT[i]);
    fPt[i] = new TH1D(std::string("Pt"+ibin).data(),";#it{p}_{T} / GeV/#it{c};N",nptbins, ptbins);
    cent_bins_list->Add(fPt[i]);
  }
  correlation_list->Add(cent_bins_list);
  fV2XXXpTcent = new TProfile2D("vnXXXptCent",";centrality V0M; #it{p}_{T} / GeV/#it{c};#LTx_{2}X_{1}^{ZNA}X_{1}^{ZNC}#GT",100,0.,100.,nptbins,ptbins); 
  fV2XYYpTcent = new TProfile2D("vnXYYptCent",";centrality V0M; #it{p}_{T} / GeV/#it{c};#LTx_{2}Y_{1}^{ZNA}Y_{1}^{ZNC}#GT",100,0.,100.,nptbins,ptbins);
  fV2YXYpTcent = new TProfile2D("vnYXYptCent",";centrality V0M; #it{p}_{T} / GeV/#it{c};#LTy_{2}X_{1}^{ZNA}Y_{1}^{ZNC}#GT",100,0.,100.,nptbins,ptbins);
  fV2YYXpTcent = new TProfile2D("vnYYXptCent",";centrality V0M; #it{p}_{T} / GeV/#it{c};#LTy_{2}Y_{1}^{ZNA}X_{1}^{ZNC}#GT",100,0.,100.,nptbins,ptbins);
  fV2YYYpTcent = new TProfile2D("vnYYYptCent",";centrality V0M; #it{p}_{T} / GeV/#it{c};#LTy_{2}Y_{1}^{ZNA}Y_{1}^{ZNC}#GT",100,0.,100.,nptbins,ptbins);
  fV2YXXpTcent = new TProfile2D("vnYXXptCent",";centrality V0M; #it{p}_{T} / GeV/#it{c};#LTy_{2}X_{1}^{ZNA}X_{1}^{ZNC}#GT",100,0.,100.,nptbins,ptbins);
  fV2XXYpTcent = new TProfile2D("vnXXYptCent",";centrality V0M; #it{p}_{T} / GeV/#it{c};#LTx_{2}X_{1}^{ZNA}Y_{1}^{ZNC}#GT",100,0.,100.,nptbins,ptbins);
  fV2XYXpTcent = new TProfile2D("vnXYXptCent",";centrality V0M; #it{p}_{T} / GeV/#it{c};#LTx_{2}Y_{1}^{ZNA}X_{1}^{ZNC}#GT",100,0.,100.,nptbins,ptbins);
  correlation_list->Add(fV2XXXpTcent);
  correlation_list->Add(fV2XYYpTcent);
  correlation_list->Add(fV2YXYpTcent);
  correlation_list->Add(fV2YYXpTcent);
  correlation_list->Add(fV2YYYpTcent);
  correlation_list->Add(fV2YXXpTcent);
  correlation_list->Add(fV2XXYpTcent);
  correlation_list->Add(fV2XYXpTcent);
  fXaXcCent = new TProfile("XXCent",";centrality V0M;#LTX_{1}^{ZNA}X_{1}^{ZNC}#GT",100,0.,100.);
  fYaYcCent = new TProfile("YYCent",";centrality V0M;#LTY_{1}^{ZNA}Y_{1}^{ZNC}#GT",100,0.,100.);
  fXaYcCent = new TProfile("XYCent",";centrality V0M;#LTX_{1}^{ZNA}Y_{1}^{ZNC}#GT",100,0.,100.);
  fYaXcCent = new TProfile("YXCent",";centrality V0M;#LTY_{1}^{ZNA}X_{1}^{ZNC}#GT",100,0.,100.);
  correlation_list->Add(fXaXcCent);
  correlation_list->Add(fYaYcCent);
  correlation_list->Add(fXaYcCent);
  correlation_list->Add(fYaXcCent);
  fXaXcDistCent = new TH2D("XXDistributionCent",";centrality V0M;#LTX_{1}^{ZNA}X_{1}^{ZNC}#GT",100,0.,100.,1000,-0.2,0.2);
  fYaYcDistCent = new TH2D("YYDistributionCent",";centrality V0M;#LTY_{1}^{ZNA}Y_{1}^{ZNC}#GT",100,0.,100.,1000,-0.2,0.2);
  fXaYcDistCent = new TH2D("XYDistributionCent",";centrality V0M;#LTX_{1}^{ZNA}Y_{1}^{ZNC}#GT",100,0.,100.,1000,-0.2,0.2);
  fYaXcDistCent = new TH2D("YXDistributionCent",";centrality V0M;#LTY_{1}^{ZNA}X_{1}^{ZNC}#GT",100,0.,100.,1000,-0.2,0.2);
  correlation_list->Add(fXaXcDistCent);
  correlation_list->Add(fYaYcDistCent);
  correlation_list->Add(fXaYcDistCent);
  correlation_list->Add(fYaXcDistCent);
  fXtXaXcDistCent = new TH2D("XXXDistributionCent","; centrality V0M;#LTx_{2}X_{1}^{ZNA}X_{1}^{ZNC}#GT",100,0.,100.,1000,-0.2,0.2);
  fXtYaYcDistCent = new TH2D("XYYDistributionCent","; centrality V0M;#LTx_{2}Y_{1}^{ZNA}Y_{1}^{ZNC}#GT",100,0.,100.,1000,-0.2,0.2);
  fYtXaYcDistCent = new TH2D("YXYDistributionCent","; centrality V0M;#LTy_{2}X_{1}^{ZNA}Y_{1}^{ZNC}#GT",100,0.,100.,1000,-0.2,0.2);
  fYtYaXcDistCent = new TH2D("YYXDistributionCent","; centrality V0M;#LTy_{2}Y_{1}^{ZNA}X_{1}^{ZNC}#GT",100,0.,100.,1000,-0.2,0.2);
  correlation_list->Add(fXtXaXcDistCent);
  correlation_list->Add(fXtYaYcDistCent);
  correlation_list->Add(fYtXaYcDistCent);
  correlation_list->Add(fYtYaXcDistCent);

  fXtXaXcCent = new TProfile("XXXCent","; centrality V0M;#LTx_{2}X_{1}^{ZNA}X_{1}^{ZNC}#GT",100,0.,100.);
  fXtYaYcCent = new TProfile("XYYCent","; centrality V0M;#LTx_{2}Y_{1}^{ZNA}Y_{1}^{ZNC}#GT",100,0.,100.);
  fYtXaYcCent = new TProfile("YXYCent","; centrality V0M;#LTy_{2}X_{1}^{ZNA}Y_{1}^{ZNC}#GT",100,0.,100.);
  fYtYaXcCent = new TProfile("YYXCent","; centrality V0M;#LTy_{2}Y_{1}^{ZNA}X_{1}^{ZNC}#GT",100,0.,100.);
  fYtXaXcCent = new TProfile("XYYCent","; centrality V0M;#LTy_{2}Y_{1}^{ZNA}X_{1}^{ZNC}#GT",100,0.,100.);
  fYtYaYcCent = new TProfile("YYYCent","; centrality V0M;#LTy_{2}Y_{1}^{ZNA}Y_{1}^{ZNC}#GT",100,0.,100.);
  fXtXaYcCent = new TProfile("XXYCent","; centrality V0M;#LTx_{2}X_{1}^{ZNA}Y_{1}^{ZNC}#GT",100,0.,100.);
  fXtYaXcCent = new TProfile("XYXCent","; centrality V0M;#LTx_{2}Y_{1}^{ZNA}X_{1}^{ZNC}#GT",100,0.,100.);
  correlation_list->Add(fXtXaXcCent);
  correlation_list->Add(fXtYaYcCent);
  correlation_list->Add(fYtXaYcCent);
  correlation_list->Add(fYtYaXcCent);
  correlation_list->Add(fYtXaXcCent);
  correlation_list->Add(fYtYaYcCent);
  correlation_list->Add(fXtXaYcCent);
  correlation_list->Add(fXtYaXcCent);

  auto bs_list = new TList();
  bs_list->SetOwner(true);
  bs_list->SetName("bs_samples");
  for (int i = 0; i < fNsamples; ++i) {
    auto addnum = [i](const std::string &name){ return name+"_bs_"+std::to_string(i); };
    fXaXcCentBS.emplace_back(new TProfile(addnum("XXCent").data(),";centrality V0M;#LTX_{1}^{ZNA}X_{1}^{ZNC}#GT",100,0.,100.));
    fYaYcCentBS.emplace_back(new TProfile(addnum("YYCent").data(),";centrality V0M;#LTY_{1}^{ZNA}Y_{1}^{ZNC}#GT",100,0.,100.));
    fXaYcCentBS.emplace_back(new TProfile(addnum("XYCent").data(),";centrality V0M;#LTX_{1}^{ZNA}Y_{1}^{ZNC}#GT",100,0.,100.));
    fYaXcCentBS.emplace_back(new TProfile(addnum("YXCent").data(),";centrality V0M;#LTY_{1}^{ZNA}X_{1}^{ZNC}#GT",100,0.,100.));
    fV2XXXpTcentBS.emplace_back(new TProfile2D(addnum("vnXXXptCent").data(),";centrality V0M; #it{p}_{T} / GeV/#it{c};#LTx_{2}X_{1}^{ZNA}X_{1}^{ZNC}#GT",100,0.,100.,nptbins,ptbins));
    fV2XYYpTcentBS.emplace_back(new TProfile2D(addnum("vnXYYptCent").data(),";centrality V0M; #it{p}_{T} / GeV/#it{c};#LTx_{2}Y_{1}^{ZNA}Y_{1}^{ZNC}#GT",100,0.,100.,nptbins,ptbins));
    fV2YXYpTcentBS.emplace_back(new TProfile2D(addnum("vnYXYptCent").data(),";centrality V0M; #it{p}_{T} / GeV/#it{c};#LTy_{2}X_{1}^{ZNA}Y_{1}^{ZNC}#GT",100,0.,100.,nptbins,ptbins));
    fV2YYXpTcentBS.emplace_back(new TProfile2D(addnum("vnYYXptCent").data(),";centrality V0M; #it{p}_{T} / GeV/#it{c};#LTy_{2}Y_{1}^{ZNA}X_{1}^{ZNC}#GT",100,0.,100.,nptbins,ptbins));
    bs_list->Add(fXaXcCentBS.back());
    bs_list->Add(fYaYcCentBS.back());
    bs_list->Add(fXaYcCentBS.back());
    bs_list->Add(fYaXcCentBS.back());
    bs_list->Add(fV2XXXpTcentBS.back());
    bs_list->Add(fV2XYYpTcentBS.back());
    bs_list->Add(fV2YXYpTcentBS.back());
    bs_list->Add(fV2YYXpTcentBS.back());
  }
  correlation_list->Add(bs_list);
  return correlation_list;
}

