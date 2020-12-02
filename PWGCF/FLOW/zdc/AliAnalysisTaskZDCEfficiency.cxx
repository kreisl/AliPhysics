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
#include <array>
#include <TFormula.h>
#include <TList.h>

#include "AliAnalysisTaskZDCEfficiency.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliCentrality.h"
#include "AliLog.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliAODMCParticle.h"
#include "AliMultSelection.h"

AliAnalysisTaskZDCEfficiency::AliAnalysisTaskZDCEfficiency() :
  AliAnalysisTaskSE(),
  fMC(true),
  fHistograms(nullptr),
  fHistPtCentrality96(nullptr),
  fHistPtCentrality128(nullptr),
  fHistPtCentrality768(nullptr),
  fHistPtCentralityMC(nullptr),
  fHistPtCentralityMCPrim(nullptr),
  fHistCentrality(nullptr) {
}

AliAnalysisTaskZDCEfficiency::AliAnalysisTaskZDCEfficiency(const char *name) :
  AliAnalysisTaskSE(name),
  fMC(true),
  fHistograms(nullptr),
  fHistPtCentrality96(nullptr),
  fHistPtCentrality128(nullptr),
  fHistPtCentrality768(nullptr),
  fHistPtCentralityMC(nullptr),
  fHistPtCentralityMCPrim(nullptr),
  fHistCentrality(nullptr) {
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskZDCEfficiency::~AliAnalysisTaskZDCEfficiency() {
  delete fAnalysisUtils;
}

void AliAnalysisTaskZDCEfficiency::UserCreateOutputObjects() {
  fAnalysisUtils = new AliAnalysisUtils();
  fAnalysisUtils->SetUseMVPlpSelection(true);
  fAnalysisUtils->SetUseOutOfBunchPileUp(true);
  if (fRun1) {
    fEventCuts.SetupRun1PbPb();
  } else {
    fEventCuts.SetupLHC17n();
  }

  const Int_t npt = 24;
  Double_t pt_bins[npt+1] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.25, 1.5,
      1.75, 2.0, 2.25, 2.5, 3., 3.5, 4., 5., 6., 8., 10., 15., 20., 30.};
  const Int_t ncent = 90;
  Double_t cent_bins[ncent+1];
  for (int i = 0; i < ncent+1; ++i) {
    cent_bins[i] = (double) i;
  }
  fHistograms = new TList();
  fHistograms->SetOwner(true);
  fHistPtCentrality96Primary = new TH2F("PtvsCentrality96Primary", ";#it{p}_{T};centrality", npt, pt_bins, ncent, cent_bins);
  fHistPtCentrality128Primary = new TH2F("PtvsCentrality128Primary", ";#it{p}_{T};centrality", npt, pt_bins, ncent, cent_bins);
  fHistPtCentrality768Primary = new TH2F("PtvsCentrality768Primary", ";#it{p}_{T};centrality", npt, pt_bins, ncent, cent_bins);

  fHistPtCentrality96 = new TH2F("PtvsCentrality96", ";#it{p}_{T};centrality", npt, pt_bins, ncent, cent_bins);
  fHistPtCentrality128 = new TH2F("PtvsCentrality128", ";#it{p}_{T};centrality", npt, pt_bins, ncent, cent_bins);
  fHistPtCentrality768 = new TH2F("PtvsCentrality768", ";#it{p}_{T};centrality", npt, pt_bins, ncent, cent_bins);
  fHistPtCentralityMC = new TH2F("PtvsCentralityMC",";#it{p}_{T};centrality", npt, pt_bins, ncent, cent_bins);
  fHistPtCentralityMCPrim = new TH2F("PtvsCentralityMCPrim",";#it{p}_{T};centrality", npt, pt_bins, ncent, cent_bins);
  fHistCentrality = new TH1F("Centrality",";centrality;N",100,0.,100.);
  fHistograms->Add(fHistPtCentrality96);
  fHistograms->Add(fHistPtCentrality128);
  fHistograms->Add(fHistPtCentrality768);
  fHistograms->Add(fHistPtCentrality96Primary);
  fHistograms->Add(fHistPtCentrality128Primary);
  fHistograms->Add(fHistPtCentrality768Primary);
  fHistograms->Add(fHistPtCentralityMC);
  fHistograms->Add(fHistPtCentralityMCPrim);
  fHistograms->Add(fHistCentrality);
  PostData(1, fHistograms);
}

void AliAnalysisTaskZDCEfficiency::UserExec(Option_t *option) {
  auto event = dynamic_cast<AliAODEvent *>(InputEvent());
  if (!event) {AliDebug(AliLog::kError, "Event not available"); return;}
  auto eventhandler =
      dynamic_cast<AliInputEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!eventhandler->IsEventSelected()) return;
  if (!fEventCuts.AcceptEvent(event)) return;
  if (fAnalysisUtils->IsPileUpEvent(event)) return;
  auto mcevent = MCEvent();
  if (!mcevent) return;

  const auto vertex = event->GetPrimaryVertex();
  double vertex_z = vertex->GetZ();
  if (std::fabs(vertex_z) > 10.) return;

  auto cent_v0m = GetCentrality(event);
  if (cent_v0m > 100. || cent_v0m < 0.) return;

  fHistCentrality->Fill(cent_v0m);
    const unsigned int ntracks = event->GetNumberOfTracks();
  for (unsigned int i = 0; i < ntracks; ++i) {
    auto track = static_cast<AliAODTrack *>(event->GetTrack(i));
    if (!track) continue;
    const auto pt = track->Pt(); 
    const auto eta = track->Eta();
    if (std::fabs(eta) > 0.8) continue;
    auto label = track->GetLabel();
    auto particle = static_cast< AliAODMCParticle* >(mcevent->GetTrack(std::abs(label)));
    bool is_physical_primary = false;
    if (particle) is_physical_primary = particle->IsPhysicalPrimary(); 
    if (track->TestFilterBit(96)) {
      fHistPtCentrality96->Fill(pt, cent_v0m);
      if (is_physical_primary) {
        fHistPtCentrality96Primary->Fill(pt, cent_v0m);
      }
    }
    if (track->TestFilterBit(128)) {
      fHistPtCentrality128->Fill(pt, cent_v0m);
      if (is_physical_primary) {
        fHistPtCentrality128Primary->Fill(pt, cent_v0m);
      }
    }
    if (track->TestFilterBit(768)) {
      fHistPtCentrality768->Fill(pt, cent_v0m);
      if (is_physical_primary) {
        fHistPtCentrality768Primary->Fill(pt, cent_v0m);
      }
    }
  }
  AliDebug(AliLog::kDebug, "Processing MC particles");
  Int_t nTracksMC   = mcevent->GetNumberOfTracks();
  for (Int_t iTr = 0; iTr < nTracksMC; iTr++) {
    auto particle = static_cast< AliAODMCParticle* >(mcevent->GetTrack(iTr));
    if (!particle) continue;
    if (particle->Charge() == 0) continue;
    auto eta = particle->Eta();
    auto pt = particle->Pt();
    if (std::fabs(eta) > 0.8) continue;
    fHistPtCentralityMC->Fill(pt, cent_v0m);
    if (particle->IsPhysicalPrimary()) {
      fHistPtCentralityMCPrim->Fill(pt, cent_v0m);
    }
  }
  PostData(1, fHistograms);
}

void AliAnalysisTaskZDCEfficiency::Terminate(Option_t *option) {}

double AliAnalysisTaskZDCEfficiency::GetCentrality(AliAODEvent* event) {
  bool run_1 = false;
  if (event->GetRunNumber() < 200000) run_1 = true;
  if (run_1) {
    // Run 1
    auto centrality = event->GetCentrality();
    if (centrality) {
      return centrality->GetCentralityPercentile("V0M");
    }
  } else {
    // Run 2
    auto centrality = dynamic_cast<AliMultSelection *>(event->FindListObject("MultSelection"));
    if (centrality) {
      return centrality->GetMultiplicityPercentile("V0M");
    }
  }
  return -1;
}