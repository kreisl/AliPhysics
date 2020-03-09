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
#include "AliAODEvent.h"
#include "AliAODTrack.h"

#include "AliZDCflowCuts.h"

void AliZDCflowCuts::CheckEventCuts(AliAODEvent* event) {
  bool selected = true;
  if(!(std::abs(event->GetPrimaryVertex()->GetZ()) < vtxZcut)) selected = false;
  fEventSelected = selected;
  fEvent = event;
}

void AliZDCflowCuts::AddQAHistograms(TList *list) {
  auto qalist = new TList();
  qalist->SetName("trackQA");
  qalist->SetOwner(true);
  fHistPhi = new TH2D("phi", ";phi;cuts applied;N", 500, 0, 2*TMath::Pi(), 2, 0., 2.);
  fHistEta = new TH2D("eta", ";eta;cuts applied;N", 200, -0.9, 0.9, 2, 0., 2.);
  fHistDCAxy = new TH2D("dcaxy", ";dcaxy;cuts applied;N", 200, 0., 5., 2, 0., 2.);
  fHistDCAz = new TH2D("dcaz", ";dcaz;cuts applied;N", 200, 0., 5., 2, 0., 2.);
  std::vector<Double_t> pt_bins {
    0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.,
    1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 3.,3.25, 3.5,
    3.75, 4., 4.5, 5., 5.5, 6., 7., 8., 9., 10.,
    12., 14., 17., 20., 25., 30., 40., 50.};    
  fHistPt= new TH2D("pt", ";pt; cuts applied;N", pt_bins.size()-1, pt_bins.data(), 2, 0., 2.);
  fHistTPCnCls = new TH1D("TPCnCls",";n cls tpc;N", 160, 0, 160);
  fHistTPCchi2perCls = new TH1D("TPCchi2perCls",";tpc chi2 per cls;N", 100, 0, 4.);
  fHistITSchi2 = new TH1D("ITSchi2",";its chi2;N", 100, 0., 36.);
  fHistTPCchi2CvsGlo = new TH1D("TPCchi2CvsGlo",";tpc chi2 const vs glo;N", 100, -36., 36.);
  fHistTPCSharedClsF = new TH1D("TPCClsSharedFraction",";n cls tpc;N", 100, 0, 1.);
  qalist->Add(fHistPhi);
  qalist->Add(fHistEta);
  qalist->Add(fHistPt);
  qalist->Add(fHistDCAxy);
  qalist->Add(fHistDCAz);
  qalist->Add(fHistTPCnCls);
  qalist->Add(fHistTPCchi2perCls);
  qalist->Add(fHistITSchi2);
  qalist->Add(fHistTPCchi2CvsGlo);
  qalist->Add(fHistTPCSharedClsF);
  list->Add(qalist);
}

bool AliZDCflowCuts::CheckTrackCutsPtCutOnly(AliAODTrack* track) {
  if (!(track->Pt() >= ptMin)) return false;
  if (!(track->Pt() <  ptMax)) return false;
  return true;
}

bool AliZDCflowCuts::CheckTrackCutsNoPtCut(AliAODTrack* track) {
  if (!track->TestFilterBit(filterBit)) return false;
  Double_t track_xyz[3]  = {0.,0.,0.};
  Double_t vertex_xyz[3] = {0.,0.,0.};
  Double_t dca_xyz[3]    = {0.,0.,0.};
  track->GetXYZ(track_xyz);
  const AliAODVertex* vertex = fEvent->GetPrimaryVertex();
  vertex->GetXYZ(vertex_xyz);
  for (auto i = 0u; i < 3; ++i) { dca_xyz[i] = track_xyz[i] - vertex_xyz[i]; }
  const auto dcaxy = std::sqrt(dca_xyz[0]*dca_xyz[0] + dca_xyz[1]*dca_xyz[1]);
  const auto dcaz  = std::fabs(dca_xyz[2]);
  const auto phi = track->Phi();
  const auto eta = track->Eta();
  const auto pt  = track->Pt();
  const auto tpc_ncls        = track->GetTPCncls();
  const auto tpc_ncls_shared = track->GetTPCnclsS();
  const auto tpc_chi2_cls    = track->GetTPCchi2perCluster();
  const auto tpc_chi2_cvg    = track->GetChi2TPCConstrainedVsGlobal();
  const auto its_chi2        = track->GetITSchi2();
  const auto fraction_tpc_shared_cls = static_cast<double>(tpc_ncls_shared) / tpc_ncls;
  // Fill QA Histograms before cuts are applied
  fHistPhi->Fill(phi, 0.);
  fHistEta->Fill(eta, 0.);
  fHistPt->Fill(pt, 0.);
  fHistDCAxy->Fill(dcaxy, 0.);
  fHistDCAz->Fill(dcaz, 0.);
  fHistTPCnCls->Fill(tpc_ncls);
  fHistTPCchi2perCls->Fill(tpc_chi2_cls);
  fHistITSchi2->Fill(its_chi2);
  fHistTPCchi2CvsGlo->Fill(tpc_chi2_cvg);
  fHistTPCSharedClsF->Fill(fraction_tpc_shared_cls);
  // check cuts
  if      (sign == -1 && !(track->GetSign() < 0)) return false;
  else if (sign ==  1 && !(track->GetSign() > 0)) return false;
  if (tpc_ncls < nClsCut)                         return false;
  if (!(tpc_chi2_cls > chi2MinCut))               return false;
  if (!(tpc_chi2_cls < chi2MaxCut))               return false;
  if (!(std::abs(track->Eta()) < etaMax))         return false;
  if (DCAzMax  > 0. && dcaz  > DCAzMax)           return false;
  if (DCAxyMax > 0. && dcaxy > DCAxyMax)          return false;
  // Fill QA Histograms after cuts are applied
  fHistPhi->Fill(phi, 1.);
  fHistEta->Fill(eta, 1.);
  fHistDCAxy->Fill(dcaxy, 1.);
  fHistDCAz->Fill(dcaz, 1.);
  fHistPt->Fill(pt, 1.);
  return true;
}

void AliZDCflowCuts::Print() {
  using std::cout;
  using std::endl;
  cout << "cuts set:" << endl
       << "vtxZcut:   " << vtxZcut    << endl
       << "sign       " << sign       << endl
       << "filterBit  " << filterBit  << endl
       << "nClsCut    " << nClsCut    << endl
       << "chi2MinCut " << chi2MinCut << endl
       << "chi2MaxCut " << chi2MaxCut << endl
       << "ptMin      " << ptMin      << endl
       << "ptMax      " << ptMax      << endl
       << "etaMax     " << etaMax     << endl
       << "DCA Z max  " << DCAzMax    << endl
       << "DCA XY max " << DCAxyMax   << endl;
}
