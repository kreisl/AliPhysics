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
#include "AliQvectors.h"
#include "AliZDCellipticFlow.h"

void AliZDCellipticFlow::FindCentralityBin(AliAODEvent *event, Double_t centrality, const std::vector<Int_t> &samples) {
  fCuts.CheckEventCuts(event);
  fCentralityBin = fCentralityAxis->FindBin(centrality) - 1;
  fCentrality = centrality;
  fCuts.HasPassedCentrality(centrality > 0. && centrality < 90.);
  fSamples = samples;
}

void AliZDCellipticFlow::FillPerEventCorrelations() {
}

void AliZDCellipticFlow::CalculateCorrelations(Qv<4,4> qtpc,
                                               std::vector<Qv<4,4>> qtpcpt,
                                               Qv<4,4> qtpcetap,
                                               Qv<4,4> qtpcetan) {
  if (!(fCuts.PassedEventCuts())) return;

  // ZDC and ZDC V0 correlations
  fXaXcCent->Fill(fCentrality, fQZA.x*fQZC.x, fCentralityWeight);
  fYaYcCent->Fill(fCentrality, fQZA.y*fQZC.y, fCentralityWeight);
  fXaYcCent->Fill(fCentrality, fQZA.x*fQZC.y, fCentralityWeight);
  fYaXcCent->Fill(fCentrality, fQZA.y*fQZC.x, fCentralityWeight);
  fXaXcDistCent->Fill(fCentrality,fQZA.x*fQZC.x);
  fYaYcDistCent->Fill(fCentrality,fQZA.y*fQZC.y);
  fXaYcDistCent->Fill(fCentrality,fQZA.x*fQZC.y);
  fYaXcDistCent->Fill(fCentrality,fQZA.y*fQZC.x);
  fVXaZXXCent->Fill(fCentrality, fQVA.x*fQZC.x*fQZA.x, fCentralityWeight);
  fVXcZXXCent->Fill(fCentrality, fQVC.x*fQZC.x*fQZA.x, fCentralityWeight);
  fVXaVXcCent->Fill(fCentrality, fQVA.x*fQVC.x, fCentralityWeight);
  fVYaZYYCent->Fill(fCentrality, fQVA.y*fQZC.y*fQZA.y, fCentralityWeight);
  fVYcZYYCent->Fill(fCentrality, fQVC.y*fQZC.y*fQZA.y, fCentralityWeight);
  fVYaVYcCent->Fill(fCentrality, fQVA.y*fQVC.y, fCentralityWeight);
  for (int isample = 0; isample < fNsamples; ++isample) {
    for (int i = 0; i < fSamples[isample]; ++i) {
        fXaXcCentBS[isample]->Fill(fCentrality, fQZA.x*fQZC.x, fCentralityWeight);
        fYaYcCentBS[isample]->Fill(fCentrality, fQZA.y*fQZC.y, fCentralityWeight);
        fXaYcCentBS[isample]->Fill(fCentrality, fQZA.x*fQZC.y, fCentralityWeight);
        fYaXcCentBS[isample]->Fill(fCentrality, fQZA.y*fQZC.x, fCentralityWeight);
    }
  }

  // ZDC and TPC correlations
  for (int i = 0; i < fAxisPt->GetNbins(); ++i) {
    auto pt = fAxisPt->GetBinCenter(i+1);
    auto qx = qtpcpt[i](2,1).real();
    auto qy = qtpcpt[i](2,1).imag();
    auto m  = qtpcpt[i](0,1).real();
    if (m < 1.) continue;
    qx /= m;
    qy /= m;
    double v2_xxx_pt = qx*fQZA.x*fQZC.x;
    double v2_xyy_pt = qx*fQZA.y*fQZC.y;

    double v2_yxy_pt = qy*fQZA.x*fQZC.y;
    double v2_yyx_pt = qy*fQZA.y*fQZC.x;
    double v2_yyy_pt = qy*fQZA.y*fQZC.y;
    double v2_yxx_pt = qy*fQZA.x*fQZC.x;

    double v2_xyx_pt = qx*fQZA.y*fQZC.x;
    double v2_xxy_pt = qx*fQZA.x*fQZC.y;
    fPtCent->Fill(fCentrality, pt, m);
    fV2XXXpTcent->Fill(fCentrality, pt, v2_xxx_pt, fCentralityWeight);
    fV2XYYpTcent->Fill(fCentrality, pt, v2_xyy_pt, fCentralityWeight);
    fV2YXYpTcent->Fill(fCentrality, pt, v2_yxy_pt, fCentralityWeight);
    fV2YYXpTcent->Fill(fCentrality, pt, v2_yyx_pt, fCentralityWeight);
    fV2YYYpTcent->Fill(fCentrality, pt, v2_yyy_pt, fCentralityWeight);
    fV2YXXpTcent->Fill(fCentrality, pt, v2_yxx_pt, fCentralityWeight);
    fV2XXYpTcent->Fill(fCentrality, pt, v2_xxy_pt, fCentralityWeight);
    fV2XYXpTcent->Fill(fCentrality, pt, v2_xyx_pt, fCentralityWeight);
    for (int isample = 0; isample < fNsamples; ++isample) {
      for (int i = 0; i < fSamples[isample]; ++i) {
        fV2XXXpTcentBS[isample]->Fill(fCentrality, pt, v2_xxx_pt, fCentralityWeight);
        fV2XYYpTcentBS[isample]->Fill(fCentrality, pt, v2_xyy_pt, fCentralityWeight);
        fV2YXYpTcentBS[isample]->Fill(fCentrality, pt, v2_yxy_pt, fCentralityWeight);
        fV2YYXpTcentBS[isample]->Fill(fCentrality, pt, v2_yyx_pt, fCentralityWeight);
      }
    }
    fV2XXXpT[fCentralityBin]->Fill(pt, v2_xxx_pt);
    fV2YYYpT[fCentralityBin]->Fill(pt, v2_yyy_pt);

    fV2XYYpT[fCentralityBin]->Fill(pt, v2_xyy_pt);
    fV2YXYpT[fCentralityBin]->Fill(pt, v2_yxy_pt);
    fV2YYXpT[fCentralityBin]->Fill(pt, v2_yyx_pt);

    fV2YXXpT[fCentralityBin]->Fill(pt, v2_yxx_pt);
    fV2XYXpT[fCentralityBin]->Fill(pt, v2_xyx_pt);
    fV2XXYpT[fCentralityBin]->Fill(pt, v2_xxy_pt);
    fPt[fCentralityBin]->Fill(pt);
  }
  auto qx = qtpc(2,1).real();
  auto qy = qtpc(2,1).imag();
  auto m  = qtpc(0,1).real();
  if (m > 0.) {
    qx /= m;
    qy /= m;
  }
  double v2_xxx = qx*fQZA.x*fQZC.x;
  double v2_xyy = qx*fQZA.y*fQZC.y;
  double v2_yxy = qy*fQZA.x*fQZC.y;
  double v2_yyx = qy*fQZA.y*fQZC.x;
  double v2_yyy = qy*fQZA.y*fQZC.y;
  double v2_yxx = qy*fQZA.x*fQZC.x;
  double v2_xyx = qx*fQZA.y*fQZC.x;
  double v2_xxy = qx*fQZA.x*fQZC.y;

  auto qxetap = qtpcetap(2,1).real();
  auto qyetap = qtpcetap(2,1).imag();
  auto metap  = qtpcetap(0,1).real();
  auto qxetan = qtpcetan(2,1).real();
  auto qyetan = qtpcetan(2,1).imag();
  auto metan  = qtpcetan(0,1).real();
  if (metap > 0. && metan > 0.) {
    qxetap /= metap;
    qyetap /= metap;
    qxetan /= metan;
    qyetan /= metan;
  }

  if (metan > 0. && metap > 0.) {
    fTXaZXXCent->Fill(fCentrality,qxetap*fQZA.x*fQZC.x, fCentralityWeight);
    fTXcZXXCent->Fill(fCentrality,qxetan*fQZA.x*fQZC.x, fCentralityWeight);
    fTXaTXcCent->Fill(fCentrality,qxetan*qxetap       , fCentralityWeight);
    fTYaZYYCent->Fill(fCentrality,qyetap*fQZA.y*fQZC.y, fCentralityWeight);
    fTYcZYYCent->Fill(fCentrality,qyetan*fQZA.y*fQZC.y, fCentralityWeight);
    fTYaTYcCent->Fill(fCentrality,qyetan*qyetap       , fCentralityWeight);
  }

  if (m > 0) {
    fXtYaYcDistCent->Fill(fCentrality, v2_xyy);
    fYtXaYcDistCent->Fill(fCentrality, v2_yxy);
    fYtYaXcDistCent->Fill(fCentrality, v2_yyx);
    fXtXaXcCent->Fill(fCentrality, v2_xxx);
    fYtYaYcCent->Fill(fCentrality, v2_yyy);
    fXtXaYcCent->Fill(fCentrality, v2_xxy);
    fXtYaXcCent->Fill(fCentrality, v2_xyx);
    fYtXaXcCent->Fill(fCentrality, v2_yxx);
    fYtYaXcCent->Fill(fCentrality, v2_yyx);
    fYtXaYcCent->Fill(fCentrality, v2_yxy);
    fXtYaYcCent->Fill(fCentrality, v2_xyy);
  }
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
  fPtCent = new TH2D("PtCent",";centrality V0M; #it{p}_{T} / GeV/#it{c};N",100,0.,100.,nptbins,ptbins);
  correlation_list->Add(fV2XXXpTcent);
  correlation_list->Add(fV2XYYpTcent);
  correlation_list->Add(fV2YXYpTcent);
  correlation_list->Add(fV2YYXpTcent);
  correlation_list->Add(fV2YYYpTcent);
  correlation_list->Add(fV2YXXpTcent);
  correlation_list->Add(fV2XXYpTcent);
  correlation_list->Add(fV2XYXpTcent);
  correlation_list->Add(fPtCent);

  fTXaZXXCent = new TProfile("tXAzXXCent",";centrality V0M;#LTX_{2}^{TPC-P}X_{1}^{ZNA}X_{1}^{ZNC}#GT",100,0.,100.);
  fTXcZXXCent = new TProfile("tXCzXXCent",";centrality V0M;#LTX_{2}^{TPC-N}X_{1}^{ZNA}X_{1}^{ZNC}#GT",100,0.,100.);
  fTXaTXcCent = new TProfile("tXtXCent",";centrality V0M;#LTX_{2}^{TPC-P}X_{2}^{TPC-N}#GT",100,0.,100.);
  fTYaZYYCent = new TProfile("tYAzYYCent",";centrality V0M;#LTX_{2}^{TPC-P}X_{1}^{ZNA}X_{1}^{ZNC}",100,0.,100.);
  fTYcZYYCent = new TProfile("tYCzYYCent",";centrality V0M;#LTX_{2}^{TPC-N}X_{1}^{ZNA}X_{1}^{ZNC}",100,0.,100.);
  fTYaTYcCent = new TProfile("tYtYCent",";centrality V0M;#LTX_{2}^{TPC-P}X_{2}^{TPC-N}#GT",100,0.,100.);

  correlation_list->Add(fTXaZXXCent);
  correlation_list->Add(fTXcZXXCent);
  correlation_list->Add(fTXaTXcCent);
  correlation_list->Add(fTYaZYYCent);
  correlation_list->Add(fTYcZYYCent);
  correlation_list->Add(fTYaTYcCent);

  fVXaZXXCent = new TProfile("vXAzXXCent",";centrality V0M;#LTX_{2}^{V0A}X_{1}^{ZNA}X_{1}^{ZNC}#GT",100,0.,100.);
  fVXcZXXCent = new TProfile("vXCzXXCent",";centrality V0M;#LTX_{2}^{V0C}X_{1}^{ZNA}X_{1}^{ZNC}#GT",100,0.,100.);
  fVXaVXcCent = new TProfile("vXvXCent",";centrality V0M;#LTX_{2}^{V0A}X_{2}^{V0C}#GT",100,0.,100.);
  fVYaZYYCent = new TProfile("vYAzYYCent",";centrality V0M;#LTX_{2}^{V0A}X_{1}^{ZNA}X_{1}^{ZNC}",100,0.,100.);
  fVYcZYYCent = new TProfile("vYCzYYCent",";centrality V0M;#LTX_{2}^{V0C}X_{1}^{ZNA}X_{1}^{ZNC}",100,0.,100.);
  fVYaVYcCent = new TProfile("vYvYCent",";centrality V0M;#LTX_{2}^{V0A}X_{2}^{V0C}#GT",100,0.,100.);

  correlation_list->Add(fVXaZXXCent);
  correlation_list->Add(fVXcZXXCent);
  correlation_list->Add(fVXaVXcCent);
  correlation_list->Add(fVYaZYYCent);
  correlation_list->Add(fVYcZYYCent);
  correlation_list->Add(fVYaVYcCent);

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
  fYtXaXcCent = new TProfile("YXXCent","; centrality V0M;#LTy_{2}X_{1}^{ZNA}X_{1}^{ZNC}#GT",100,0.,100.);
  fYtYaYcCent = new TProfile("YYYCent","; centrality V0M;#LTy_{2}Y_{1}^{ZNA}Y_{1}^{ZNC}#GT",100,0.,100.);
  fXtXaYcCent = new TProfile("XXYCent","; centrality V0M;#LTx_{2}X_{1}^{ZNA}Y_{1}^{ZNC}#GT",100,0.,100.);
  fXtYaXcCent = new TProfile("XYXCent","; centrality V0M;#LTx_{2}Y_{1}^{ZNA}X_{1}^{ZNC}#GT",100,0.,100.);
  correlation_list->Add(fXtXaXcCent);
  correlation_list->Add(fYtYaYcCent);

  correlation_list->Add(fXtYaYcCent);
  correlation_list->Add(fYtXaYcCent);
  correlation_list->Add(fYtYaXcCent);

  correlation_list->Add(fXtXaYcCent);
  correlation_list->Add(fXtYaXcCent);
  correlation_list->Add(fYtXaXcCent);

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

