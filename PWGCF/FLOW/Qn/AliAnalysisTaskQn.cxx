#include <iostream>
#include <cstdlib>

#include "TGrid.h"
#include "TChain.h"

#include "AliAnalysisTaskSE.h"

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliMultSelection.h"
#include "AliMultEstimator.h"
#include "AliAnalysisTaskQn.h"

ClassImp(AliAnalysisTaskQn); // import class inheriting from TObject

using std::cout;
using std::endl;

AliAnalysisTaskQn::AliAnalysisTaskQn() : AliAnalysisTaskSE(),
  fCalibFileType(CalibFile::local),
  fQnManager(nullptr),
  fInputEvent(nullptr),
  fDetectorQAList(nullptr),
  fTestHist(nullptr),
  fQnTree(nullptr),
  fInCalib(nullptr) { }

AliAnalysisTaskQn::AliAnalysisTaskQn(const char *name) : AliAnalysisTaskSE(name),
  fCalibFileType(CalibFile::local),
  fQnManager(new Qn::CorrectionManager()),
  fInputEvent(nullptr),
  fDetectorQAList(nullptr),
  fTestHist(nullptr),
  fQnTree(nullptr),
  fInCalib(nullptr) { 
  // Input slot #0 works with a TChain
  DefineInput(0,TChain::Class());
  // Output slot #0 is reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  // DefineOutput(0,TTree::Class());
  DefineOutput(OutputSlot::QnCalibration,TList::Class());
  DefineOutput(OutputSlot::QnQA,TList::Class());
  DefineOutput(OutputSlot::DetectorQA,TList::Class());
  DefineOutput(OutputSlot::QnTree,TTree::Class());
}

AliAnalysisTaskQn::~AliAnalysisTaskQn() {
  cout<<"Default class destructor called.  Analysis fun time has ended"<<endl;
}

void AliAnalysisTaskQn::UserCreateOutputObjects() {
  fQnTree = new TTree("tree","tree");
  fDetectorQAList= new TList();
  fDetectorQAList->SetName("DetectorQA");
  fDetectorQAList->SetOwner(kTRUE);
  fTestHist = new TH1F("TestHist","TestHist;pt;n",100,0,10);
  fDetectorQAList->Add(fTestHist);
  fQnManager->SetTree(fQnTree);
  fQnManager->Initialize(fInCalib);
  PostData(OutputSlot::QnCalibration,fQnManager->GetCalibrationList());
  PostData(OutputSlot::QnQA,fQnManager->GetCalibrationQAList());
  PostData(OutputSlot::DetectorQA,fDetectorQAList);
  PostData(OutputSlot::QnTree,fQnTree);
}

void AliAnalysisTaskQn::FinishTaskOutput() {
  fQnManager->FillDetectorQAToList(fDetectorQAList);
}

void AliAnalysisTaskQn::UserExec(Option_t *) {
  fQnManager->Reset();
  fInputEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fInputEvent) return;
  fValues = fQnManager->GetVariableContainer();
  auto multselection = (AliMultSelection*)fInputEvent->FindListObject("MultSelection");
  fValues[kCentV0M] = multselection->GetMultiplicityPercentile("V0M");
  //ZDC
  auto zdc = fInputEvent->GetZDCData();
  constexpr std::array<double, 5> zdcX = {{0, 1.75, -1.75, 1.75, -1.75}};
  constexpr std::array<double, 5> zdcY = {{0, -1.75, -1.75, 1.75, 1.75}};
  for (int ich = 0; ich < 5; ++ich) {
    fValues[kZDCCChMult + ich] = zdc->GetZNCTowerEnergy()[ich];
    if (fValues[kZDCCChMult + ich] < 0)  fValues[kZDCCChMult + ich] = 0.f;
    fValues[kZDCCChPhi + ich] = TMath::ATan2(zdcY[ich], zdcX[ich]);
  }
  //case 2015
  for (int ich = 0; ich < 5; ++ich) {
    fValues[kZDCAChPhi + ich] = TMath::ATan2(zdcY[ich], zdcX[ich]);
    if (ich!= 2) fValues[kZDCAChMult + ich] = zdc->GetZNATowerEnergy()[ich];
    if (ich == 2) {
      fValues[kZDCAChMult + ich] = zdc->GetZNATowerEnergy()[0] - zdc->GetZNATowerEnergy()[1] - zdc->GetZNATowerEnergy()[3] - zdc->GetZNATowerEnergy()[4];
    } // special case for ZDC-A channel 2 for 2015 datataking
    if (fValues[kZDCCChMult + ich] < 0)  fValues[kZDCCChMult + ich] = 0.f;
  }
  //case 2010 2011
  for (int ich = 0; ich < 5; ++ich) {
    fValues[kZDCAChPhi + ich] = TMath::ATan2(zdcY[ich], zdcX[ich]);
    fValues[kZDCAChMult + ich] = zdc->GetZNATowerEnergy()[ich];
    if (fValues[kZDCCChMult + ich] < 0)  fValues[kZDCCChMult + ich] = 0.f;
  }
  //V0
  auto vzero = fInputEvent->GetVZEROData();
  if (vzero) {
    // V0-C
    constexpr std::array<double, 8> X = {{0.92388, 0.38268, -0.38268, -0.92388, -0.92388, -0.38268, 0.38268, 0.92388}};
    constexpr std::array<double, 8> Y = {{0.38268, 0.92388, 0.92388, 0.38268, -0.38268, -0.92388, -0.92388, -0.38268}};
    for (int ich = 0; ich < 32; ++ich) {
      fValues[kV0CChPhi+ich] = TMath::ATan2(Y[ich%8], X[ich%8]); 
      fValues[kV0CChMult+ich] = vzero->GetMultiplicityV0C(ich); 
      fValues[kV0CChRing+ich] = ich/8;
    }
    // V0-A
    for (int ich = 0; ich < 32; ++ich) {
      fValues[kV0AChPhi+ich] = TMath::ATan2(Y[ich%8], X[ich%8]); 
      fValues[kV0AChMult+ich] = vzero->GetMultiplicityV0A(ich); 
      std::cout << fValues[kV0AChMult+ich] << std::endl;
      fValues[kV0AChRing+ich] = ich/8;
    }
  }
  fQnManager->ProcessEvent();
  fQnManager->FillChannelDetectors();
  int ntracks = fInputEvent->GetNumberOfTracks();
  for (int i=0;i<ntracks;++i) {
    auto track = static_cast<AliAODTrack*>(fInputEvent->GetTrack(i));
    if (!track) continue;
    if (!(track->TestFilterBit(256) || track->TestFilterBit(512))) continue;
    fTestHist->Fill(track->Pt());
    fValues[kPhi] = track->Phi();
    fValues[kPt] = track->Pt();
    fValues[kEta] = track->Eta();
    fQnManager->FillTrackingDetectors();
  }
  fQnManager->ProcessQnVectors();
}      

void AliAnalysisTaskQn::Terminate(Option_t *option) {
  fQnManager->Finalize();
}

void AliAnalysisTaskQn::SetCalibrationFile(TString name, CalibFile type) {
  TFile *calibfile = nullptr;
  fQnManager->SetProcessName("test");
  switch (type) {
    case CalibFile::local:
      if(name.Contains("alien")) TGrid::Connect("alien://");
      calibfile = TFile::Open(name.Data());
      if (calibfile != nullptr && calibfile->IsOpen()) {
        fInCalib = calibfile;
        std::cout << "setting calibration file to: " << calibfile->GetName() << std::endl;
      }
      break;
    default:
      AliFatal("Calibration file not supported");
      break;
  }
}


