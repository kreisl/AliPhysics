#include <iostream>
#include <cstdlib>

#include "AliAnalysisTaskSE.h"

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAnalysisTaskQn.h"

ClassImp(AliAnalysisTaskQn); // import class inheriting from TObject

using std::cout;
using std::endl;

AliAnalysisTaskQn::AliAnalysisTaskQn() : AliAnalysisTaskSE(),
  fQnManager(nullptr),
  fInputEvent(nullptr),
  fTestHist(nullptr),
  fOutputList(nullptr)
{
  cout<<"Default class constructor called.  Prepare for fun analysis time!"<<endl;
}

AliAnalysisTaskQn::AliAnalysisTaskQn(const char *name) : AliAnalysisTaskSE(name),
  fQnManager(nullptr),
  fInputEvent(nullptr),
  fTestHist(nullptr),
  fOutputList(nullptr)
{
  // Input slot #0 works with a TChain
  DefineInput(0,TChain::Class());
  // Output slot #0 is reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  //DefineOutput(0,TTree::Class());
  DefineOutput(1,TList::Class());
}

AliAnalysisTaskQn::~AliAnalysisTaskQn()
{
  cout<<"Default class destructor called.  Analysis fun time has ended"<<endl;
}

void AliAnalysisTaskQn::UserCreateOutputObjects()
{
  fOutputList = new TList();
  fOutputList->SetName(GetName());
  fOutputList->SetOwner(kTRUE);

  fTestHist = new TH1F("TestHist","TestHist;pt;n",100,0,10);
  fOutputList->Add(fTestHist);

  PostData(1,fOutputList);
}

void AliAnalysisTaskQn::UserExec(Option_t *) 
{
  fInputEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fInputEvent) return;
  int ntracks{fInputEvent->GetNumberOfTracks()};
  for (int i=0;i<ntracks;++i) {
    auto track = static_cast<AliAODTrack*>(fInputEvent->GetTrack(i));
    if (!track) continue;
    fTestHist->Fill(track->Pt());
  }
}      

void AliAnalysisTaskQn::Terminate(Option_t\* option)
{
}

