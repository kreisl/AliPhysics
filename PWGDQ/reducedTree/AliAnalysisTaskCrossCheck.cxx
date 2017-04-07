#include <iostream>

#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TChain.h"

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODForwardMult.h"
#include "AliAnalysisTaskCrossCheck.h"

ClassImp(AliAnalysisTaskCrossCheck)

AliAnalysisTaskCrossCheck::AliAnalysisTaskCrossCheck()
  : AliAnalysisTaskSE("AliAnalysisTaskCrossCheck"),
  fOutputContainer(),
  fHistDcaXY(),
  fHistDcaZ(),
  fHistDcaXYg(),
  fHistDcaZg(),
  fHistSumFmd() {
}

AliAnalysisTaskCrossCheck::AliAnalysisTaskCrossCheck(const char *name)
  : AliAnalysisTaskSE("AliAnalysisTaskCrossCheck"),
  fOutputContainer(0),
  fHistDcaXY(0),
  fHistDcaZ(0),
  fHistDcaXYg(0),
  fHistDcaZg(0),
  fHistSumFmd(0) {
  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}

AliAnalysisTaskCrossCheck::~AliAnalysisTaskCrossCheck() {
}

void AliAnalysisTaskCrossCheck::UserCreateOutputObjects() {
  fHistDcaXY = new TH1D("fHistDcaXY","",400,-10,10);
  fHistDcaZ = new TH1D("fHistDcaZ","",400,-10,10);
  fHistDcaXYg = new TH1D("fHistDcaXYg","global dcaxy",400,-10,10);
  fHistDcaZg   = new TH1D("fHistDcaZg","global dcaz",400,-10,10);
  fHistSumFmd = new TH2D("fHistSumFmd","",200,-4,6,20,0,2*TMath::Pi());
  fOutputContainer = new TList();
  fOutputContainer->SetOwner(kTRUE);
  fOutputContainer->Add(fHistDcaXY);
  fOutputContainer->Add(fHistDcaZ);
  fOutputContainer->Add(fHistDcaXYg);
  fOutputContainer->Add(fHistDcaZg);
  fOutputContainer->Add(fHistSumFmd);
  PostData(1, fOutputContainer);
}

void AliAnalysisTaskCrossCheck::UserExec(Option_t *option) {
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *eventHandler = dynamic_cast<AliInputEventHandler*>(man->GetInputEventHandler());
  AliAODEvent *event = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!event) return;
  AliAODVertex *vertex = event->GetPrimaryVertex();
  if(!vertex) return;
  Int_t nTracks = event->GetNumberOfTracks();
  for (Int_t i = 0; i < nTracks; ++i) {
    AliAODTrack *track = static_cast<AliAODTrack*>(event->GetTrack(i));
    AliExternalTrackParam etp;
    if (!track->TestFilterBit(BIT(0))) continue;
    if (!track->IsOn(AliAODTrack::kTPCrefit)) continue;
    Double_t covar[3] ={-999,-999,-999};
    Double_t dcaGlobal[2] = {-999,-999};
    Double_t dca[2] = {-999,-999};
    dca[0] = track->DCA();
    dca[1] = track->ZAtDCA();
    etp.CopyFromVTrack(track);
    Double_t initialX = etp.GetX();
    Bool_t bProp = kFALSE;
    bProp = etp.PropagateToDCA(vertex,InputEvent()->GetMagneticField(),1e9,dcaGlobal,covar);
    if (track->TestBit(AliAODTrack::kIsDCA)) {
      fHistDcaXY->Fill(dca[0]);
      fHistDcaZ->Fill(dca[1]);
    }
    if (bProp && initialX < 3. && !track->TestBit(AliAODTrack::kIsDCA)) {
      fHistDcaXYg->Fill(dcaGlobal[0]);
      fHistDcaZg->Fill(dcaGlobal[1]);
    }
  }
  TObject *obj = event->FindListObject("Forward");
  if (obj) { 
    AliAODForwardMult *aodForward = static_cast<AliAODForwardMult*>(obj);
    fHistSumFmd->Add(&(aodForward->GetHistogram()));
    std::cout << aodForward->GetHistogram().GetNbinsX() << std::endl;
  } 
  PostData(1, fOutputContainer);
}

void AliAnalysisTaskCrossCheck::Terminate(Option_t*) {
}
