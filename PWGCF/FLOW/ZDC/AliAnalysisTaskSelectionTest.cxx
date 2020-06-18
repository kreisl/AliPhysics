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

#include "TList.h"
#include "TH1F.h"
#include "TString.h"
#include "AliLog.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisTaskSelectionTest.h"

AliAnalysisTaskSelectionTest::AliAnalysisTaskSelectionTest() : 
  AliAnalysisTaskSE() {}

AliAnalysisTaskSelectionTest::AliAnalysisTaskSelectionTest(const char *name) : 
  AliAnalysisTaskSE(name) {
  DefineOutput(1,TList::Class());
}

AliAnalysisTaskSelectionTest::~AliAnalysisTaskSelectionTest() {
  delete fHistogramList;
}

void AliAnalysisTaskSelectionTest::UserCreateOutputObjects() {
  fHistogramList = new TList();
  fHistogramList->SetOwner(true);
  fTPCnCLs = new TH1D("TPCnCls",";TPC number of clusters;number of tracks",
                      160, 0., 160.);
  fHistogramList->Add(fTPCnCls);
  PostData(1,fHistogramList);
}

void AliAnalysisTaskSelectionTest::UserExec(Option_t*) {
  // Event Selection
  auto event = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!event) return;
  auto eventhandler = dynamic_cast<AliInputEventHandler*>(
      AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!eventhandler->IsEventSelected()) return;
  // Track Loop
  const unsigned int ntracks = event->GetNumberOfTracks();
  for (unsigned int i = 0; i < ntracks; ++i) {
    auto track = static_cast<AliAODTrack*>(event->GetTrack(i));
    if (track->TestFilterBit(768)) {
      auto ncls = track->GetTPCncls();
      fTPCnCls->Fill(ncls);
      if (ncls < 70) {
        AliWarning(Form("Number of clusters: %i", ncls))
      }
    }
  }
  PostData(1,fHistogramList);
}
