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

AliAnalysisTaskZDCEfficiency::AliAnalysisTaskZDCEfficiency() :
  AliAnalysisTaskSE(),
  fMC(true),
  fEvent(nullptr),
  fMCEvent(nullptr),
  fHistograms(nullptr),
  fHistPtCentrality96(nullptr),
  fHistPtCentrality128(nullptr),
  fHistPtCentrality768(nullptr),
  fHistPtCentralityMC(nullptr),
  fHistCentrality(nullptr) {
}

AliAnalysisTaskZDCEfficiency::AliAnalysisTaskZDCEfficiency(const char *name) :
  AliAnalysisTaskSE(name),
  fMC(true),
  fEvent(nullptr),
  fMCEvent(nullptr),
  fHistograms(nullptr),
  fHistPtCentrality96(nullptr),
  fHistPtCentrality128(nullptr),
  fHistPtCentrality768(nullptr),
  fHistPtCentralityMC(nullptr),
  fHistCentrality(nullptr) {
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskZDCEfficiency::~AliAnalysisTaskZDCEfficiency() {
  delete fAnalysisUtils;
  delete fHistograms;
}

void AliAnalysisTaskZDCEfficiency::UserCreateOutputObjects() {
  fAnalysisUtils = new AliAnalysisUtils();
  fAnalysisUtils->SetUseMVPlpSelection(true);
  fAnalysisUtils->SetUseOutOfBunchPileUp(true);
  fEventCuts.SetupRun1PbPb();

  const Int_t npt = 24;
  Double_t pt_bins[npt+1] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.25, 1.5,
      1.75, 2.0, 2.25, 2.5, 3., 3.5, 4., 5., 6., 8., 10., 15., 20., 30.};
  const Int_t ncent = 9;
  Double_t cent_bins[ncent+1] = {0., 5., 10., 20., 30., 40., 50., 60., 70., 80.};
  auto fHistograms = new TList();
  fHistograms->SetOwner(true);
  fHistPtCentrality96 = new TH2F("PtvsCentrality96", ";#it{p}_{T};centrality", npt, pt_bins, ncent, cent_bins);
  fHistPtCentrality128 = new TH2F("PtvsCentrality128", ";#it{p}_{T};centrality", npt, pt_bins, ncent, cent_bins);
  fHistPtCentrality768 = new TH2F("PtvsCentrality768", ";#it{p}_{T};centrality", npt, pt_bins, ncent, cent_bins);
  fHistPtCentralityMC = new TH2F("PtvsCentralityMC",";#it{p}_{T};centrality", npt, pt_bins, ncent, cent_bins);
  fHistCentrality = new TH1F("Centrality",";centrality;N",100,0.,100.);
  fHistograms->Add(fHistPtCentrality96);
  fHistograms->Add(fHistPtCentrality128);
  fHistograms->Add(fHistPtCentrality768);
  fHistograms->Add(fHistPtCentralityMC);
  fHistograms->Add(fHistCentrality);
  PostData(1, fHistograms);
}

void AliAnalysisTaskZDCEfficiency::UserExec(Option_t *option) {
  auto event = InputEvent();
  if (!event) { AliDebug(AliLog::kError, "Event not found."); return; }
  
  // MC
  AliHeader* header = 0;
  AliStack* stack = 0;
  if (fMC) {
      fMCEvent = dynamic_cast<AliMCEvent*>(MCEvent());
      if(!fMCEvent) {AliDebug(AliLog::kError, "MCEvent not available"); return; }
      header = fMCEvent->Header();
      if (!header) {AliDebug(AliLog::kError, "Header not available"); return; }
      stack = fMCEvent->Stack();
      if (!stack) { AliDebug(AliLog::kError, "Stack not available");  return; }
	}

  auto event = dynamic_cast<AliAODEvent *>(InputEvent());
  if (!event) {AliDebug(AliLog::kError, "Event not available"); return;}
  auto eventhandler =
      dynamic_cast<AliInputEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!eventhandler->IsEventSelected()) return;
  if (!fEventCuts.AcceptEvent(event)) return;
  if (fAnalysisUtils->IsPileUpEvent(event)) return;

  const auto vertex = fEvent->GetPrimaryVertex();
  double vertex_z = vertex->GetZ();
  if (std::fabs(vertex_z) > 10.) return;

  if (fMC) {
    auto genHeader = header->GenEventHeader();
    if (!genHeader) { AliDebug(AliLog::kError, "Could not retrieve genHeader from Header"); return; }
    TArrayF vtxMC(3);
    genHeader->PrimaryVertex(vtxMC);
    if (std::fabs(vtxMC[2]) > 10.) return;
  }

  auto centrality = fEvent->GetCentrality();
  auto cent_v0m = centrality->GetCentralityPercentile("V0M");
  if (cent_v0m > 100. || cent_v0m < 0.) return;

  fHistCentrality->Fill(cent_v0m);
    const unsigned int ntracks = event->GetNumberOfTracks();
  for (unsigned int i = 0; i < ntracks; ++i) {
    auto track = static_cast<AliAODTrack *>(event->GetTrack(i));
    if (!track) continue;
    const auto pt = track->Pt(); 
    if (track->TestFilterBit(96)) {
      fHistPtCentrality768->Fill(pt, cent_v0m);
    }
    if (track->TestFilterBit(128)) {
      fHistPtCentrality128->Fill(pt, cent_v0m);
    }
    if (track->TestFilterBit(768)) {
      fHistPtCentrality768->Fill(pt, cent_v0m);
    }
  }

  /// mc particle loop
  if (fMC) {
    for (Int_t i = 0; i < stack->GetNtrack(); ++i) {
      TParticle* particle = stack->Particle(i);
      if (!particle) continue;
      if (!particle->GetPDG()) continue;
      Double_t charge = particle->GetPDG()->Charge()/3.;
      if (std::fabs(charge) < 0.001) continue;
      auto is_primary = stack->IsPhysicalPrimary(i);
      Float_t eta = particle->Eta();
      Float_t phi = particle->Phi();
      Float_t pt = particle->Pt();
      if(std::fabs(eta) > 0.8) continue;
      if(pt < 0.2 || pt > 1e15) continue;
      fHistPtCentralityMC->Fill(pt, cent_v0m);
    }
  }
  PostData(1, fHistograms);
}

void AliAnalysisTaskZDCEfficiency::Terminate(Option_t *option) {}