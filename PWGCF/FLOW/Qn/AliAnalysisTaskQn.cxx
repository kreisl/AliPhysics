#include <iostream>
#include <cstdlib>

#include "TGrid.h"
#include "TChain.h"
#include "TProfile.h"

#include "AliAnalysisTaskSE.h"

#include "AliLHCData.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODTracklets.h"
#include "AliMultSelection.h"
#include "AliMultEstimator.h"
#include "AliAnalysisTaskQn.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

ClassImp(AliAnalysisTaskQn);

using std::cout;
using std::endl;

AliAnalysisTaskQn::AliAnalysisTaskQn() : AliAnalysisTaskSE(),
  fCalibFileType(CalibFile::local),
  fQnManager(nullptr),
  fInputEvent(nullptr),
  fInDetectorList(nullptr),
  fQnTree(nullptr),
  fInCalib(nullptr),
  fOCDBPath(),
  fOCDBAvailable(),
  fLHCData(nullptr),
  fEventNumber(0) { }

AliAnalysisTaskQn::AliAnalysisTaskQn(const char *name) : AliAnalysisTaskSE(name),
  fCalibFileType(CalibFile::local),
  fQnManager(new Qn::CorrectionManager()),
  fInputEvent(nullptr),
  fInDetectorList(nullptr),
  fQnTree(nullptr),
  fInCalib(nullptr),
  fOCDBPath(),
  fOCDBAvailable(false),
  fLHCData(nullptr),
  fEventNumber(0) { 

  DefineInput(0,TChain::Class());
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
  fQnManager->SetTree(fQnTree);
  fQnManager->Initialize(fInCalib);
  if (fInCalib) fInDetectorList = (TList*) fInCalib->Get("DetectorQA");
  PostData(OutputSlot::QnCalibration,fQnManager->GetCalibrationList());
  PostData(OutputSlot::QnQA,fQnManager->GetCalibrationQAList());
  PostData(OutputSlot::DetectorQA,fQnManager->GetEventAndDetectorQAList());
  PostData(OutputSlot::QnTree,fQnTree);
}

void AliAnalysisTaskQn::NotifyRun() {
  AliInfo(TString::Format("New run number: %d", this->fCurrentRunNumber).Data());
  fQnManager->SetProcessName(TString::Format("%d", this->fCurrentRunNumber).Data());
  if (!fOCDBPath.empty()) {
    AliCDBManager::Instance()->SetDefaultStorage(fOCDBPath.data());
    AliCDBManager::Instance()->SetRun(this->fCurrentRunNumber);
    fLHCData = (AliLHCData*)((AliCDBEntry*)(AliCDBManager::Instance()->Get("GRP/GRP/LHCData")))->GetObject();
  }
  if (fLHCData) fOCDBAvailable = true;
	fEventNumber = 0;
} 

void AliAnalysisTaskQn::FinishTaskOutput() {
}

void AliAnalysisTaskQn::UserExec(Option_t *) {
  fQnManager->Reset();
  fInputEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fInputEvent) return;
  if(!(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & (AliVEvent::kMB+AliVEvent::kINT7+AliVEvent::kCentral+AliVEvent::kSemiCentral))) return;
  fValues = fQnManager->GetVariableContainer();
  const AliAODVertex* vtxtrk = fInputEvent->GetPrimaryVertex();
  if (vtxtrk->GetZ() > 10 || vtxtrk->GetZ() < -10) return;
  //new vertex selection as Alexandru Dobrin 8.5.19
  const AliAODVertex* vtxspd = fInputEvent->GetPrimaryVertexSPD();
  double covtrk[6], covspd[6];
  vtxtrk->GetCovarianceMatrix(covtrk);
  vtxspd->GetCovarianceMatrix(covspd);
  double dz = vtxtrk->GetZ() - vtxspd->GetZ();
  double errortotal = sqrt(covtrk[5]+covspd[5]);
  double errortracks = sqrt(covtrk[5]);
  double nsigtotal = dz/errortotal;
  double nsigtracks = dz/errortracks;
  if (TMath::Abs(dz) > 0.2 || TMath::Abs(nsigtotal) > 10 || TMath::Abs(nsigtracks) > 20) return; // bad vertexing
  fValues[kVtxX] = vtxtrk->GetX(); 
  fValues[kVtxY] = vtxtrk->GetY(); 
  fValues[kVtxZ] = vtxtrk->GetZ(); 
  bool old_centrality = false;
  if (fInputEvent->GetRunNumber() < 200000) old_centrality = true; // true if runnumber <= run1
  if (!old_centrality) { // run 2
    auto multselection = (AliMultSelection*)fInputEvent->FindListObject("MultSelection");
    if (multselection) {
      fValues[kCentV0A] = multselection->GetMultiplicityPercentile("V0A");
      fValues[kCentV0C] = multselection->GetMultiplicityPercentile("V0C");
      fValues[kCentV0M] = multselection->GetMultiplicityPercentile("V0M");
      fValues[kCentZNC] = multselection->GetMultiplicityPercentile("ZNC");
      fValues[kCentZNA] = multselection->GetMultiplicityPercentile("ZNA");
      fValues[kCentCL0] = multselection->GetMultiplicityPercentile("CL0");
      fValues[kCentCL1] = multselection->GetMultiplicityPercentile("CL1");
    }
  }
  if (old_centrality) {
    auto centrality = fInputEvent->GetCentrality(); // run 1 
    if (centrality) {
      fValues[kCentV0A] = centrality->GetCentralityPercentile("V0A");
      fValues[kCentV0C] = centrality->GetCentralityPercentile("V0C");
      fValues[kCentV0M] = centrality->GetCentralityPercentile("V0M");
      fValues[kCentZNC] = centrality->GetCentralityPercentile("ZNC");
      fValues[kCentZNA] = centrality->GetCentralityPercentile("ZNA");
      fValues[kCentCL0] = centrality->GetCentralityPercentile("CL0");
      fValues[kCentCL1] = centrality->GetCentralityPercentile("CL1");
    }
  }
  if (fValues[kCentV0M] < 0. || fValues[kCentV0M] > 90.) return;
  auto trigger = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  if (trigger & (AliVEvent::kMB+AliVEvent::kINT7)) {
    fValues[kTrigger] = 0;
  }
  if (trigger & AliVEvent::kCentral && !(trigger & (AliVEvent::kMB+AliVEvent::kINT7))) {
    fValues[kTrigger] = 1;
  }
  if (trigger & AliVEvent::kSemiCentral && !(trigger & (AliVEvent::kMB+AliVEvent::kINT7))) {
    fValues[kTrigger] = 2;
  }

	fValues[kRunNumber] = fInputEvent->GetRunNumber();
	fValues[kEventNumber] = fEventNumber;
	++fEventNumber;

  fValues[kNTracklets] = fInputEvent->GetMultiplicity()->GetNumberOfTracklets();
  for (int ilayer = 0; ilayer < 6; ++ilayer) {
    fValues[kNITSClusters + ilayer] = fInputEvent->GetNumberOfITSClusters(ilayer);
  }

  fValues[kNTPCClusters] = fInputEvent->GetNumberOfTPCClusters();

  //ZDC
  auto zdc = fInputEvent->GetZDCData();
  constexpr std::array<double, 5> zdcX = {{0, 1.75, -1.75, 1.75, -1.75}};
  constexpr std::array<double, 5> zdcY = {{0, -1.75, -1.75, 1.75, 1.75}};
  float summult = 0.;
  for (int ich = 0; ich < 5; ++ich) {
    fValues[kZDCCChMult + ich] = zdc->GetZNCTowerEnergy()[ich];
    if (fValues[kZDCCChMult + ich] < 0)  fValues[kZDCCChMult + ich] = 0.f;
    fValues[kZDCCChPhi + ich] = TMath::ATan2(zdcY[ich], zdcX[ich]);
    summult += fValues[kZDCCChMult+ich];
  }
  fValues[kZDCCSumMult] = summult;
  if (fInputEvent->GetRunNumber() >= 240000 && fInputEvent->GetRunNumber() <= 248000) { // in case of 2015 data one channel was offline
    for (int ich = 0; ich < 5; ++ich) {
      fValues[kZDCAChPhi + ich] = TMath::ATan2(zdcY[ich], zdcX[ich]);
      if (ich!= 2) fValues[kZDCAChMult + ich] = zdc->GetZNATowerEnergy()[ich];
      if (ich == 2) {
        fValues[kZDCAChMult + ich] = zdc->GetZNATowerEnergy()[0] - zdc->GetZNATowerEnergy()[1] - zdc->GetZNATowerEnergy()[3] - zdc->GetZNATowerEnergy()[4];
      } // special case for ZDC-A channel 2 for 2015 datataking
      if (fValues[kZDCAChMult + ich] < 0)  fValues[kZDCAChMult + ich] = 0.f;
    }
  } else { // all other data periods
    summult = 0.;
    for (int ich = 0; ich < 5; ++ich) {
      fValues[kZDCAChPhi + ich] = TMath::ATan2(zdcY[ich], zdcX[ich]);
      fValues[kZDCAChMult + ich] = zdc->GetZNATowerEnergy()[ich];
      if (fValues[kZDCAChMult + ich] < 0)  fValues[kZDCAChMult + ich] = 0.f;
      summult += fValues[kZDCAChMult+ich];
    }
  fValues[kZDCASumMult] = summult;
  }
  //V0
  float vzeromult = 0.;
  auto vzero = fInputEvent->GetVZEROData();
  if (vzero) {
    // V0-C
    constexpr std::array<double, 8> X = {{0.92388, 0.38268, -0.38268, -0.92388, -0.92388, -0.38268, 0.38268, 0.92388}};
    constexpr std::array<double, 8> Y = {{0.38268, 0.92388, 0.92388, 0.38268, -0.38268, -0.92388, -0.92388, -0.38268}};
    for (int ich = 0; ich < 32; ++ich) {
      fValues[kV0CChPhi+ich] = TMath::ATan2(Y[ich%8], X[ich%8]); 
      fValues[kV0CChMult+ich] = vzero->GetMultiplicityV0C(ich); 
      fValues[kV0CChRing+ich] = ich/8;
      vzeromult += vzero->GetMultiplicityV0C(ich);
    }
    // V0-A
    for (int ich = 0; ich < 32; ++ich) {
      fValues[kV0AChPhi+ich] = TMath::ATan2(Y[ich%8], X[ich%8]); 
      fValues[kV0AChMult+ich] = vzero->GetMultiplicityV0A(ich); 
      fValues[kV0AChRing+ich] = ich/8;
      vzeromult += vzero->GetMultiplicityV0A(ich);
    }
    fValues[kV0Mult] = vzeromult;
  }
  fQnManager->ProcessEvent();
  fQnManager->FillChannelDetectors();
  int ntracksfilterbit = 0;
  int ntracks = fInputEvent->GetNumberOfTracks();
  for (int i=0;i<ntracks;++i) {
    auto track = static_cast<AliAODTrack*>(fInputEvent->GetTrack(i));
    if (!track) continue;
    if (!track->TestFilterBit(128)) continue;
    ++ntracksfilterbit;
    fValues[kPhi] = track->Phi();
    fValues[kPt] = track->Pt();
    fValues[kEta] = track->Eta();
    fQnManager->FillTrackingDetectors();
  }
  fValues[kNTPCTracks] = ntracksfilterbit;
  fQnManager->ProcessQnVectors();
}      

void AliAnalysisTaskQn::Terminate(Option_t *option) {
  fQnManager->Finalize();
}

void AliAnalysisTaskQn::SetCalibrationFile(TString name, CalibFile type) {
  TFile *calibfile = nullptr;
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

std::array<double, 2> AliAnalysisTaskQn::GetZDCQ(AliAnalysisTaskQn::ZDC zdctype) {
  auto zdc = fInputEvent->GetZDCData();
  constexpr std::array<double, 4> zdcX = {{1.75, -1.75, 1.75, -1.75}};
  constexpr std::array<double, 4> zdcY = {{-1.75, -1.75, 1.75, 1.75}};
  std::array<double,2> qvectorsum = {0.,0.};
  double multsum = 0.;
  double energy = 0.;
  if (zdctype == ZDC::A) {
    if (fInputEvent->GetRunNumber() >= 240000 && fInputEvent->GetRunNumber() <= 248000) { // in case of 2015 data one channel was offline
      for (int ich = 1; ich < 5; ++ich) {
        if (ich!= 2) {
          energy = zdc->GetZNATowerEnergy()[ich];
        }
        if (ich == 2) { 
          energy = zdc->GetZNATowerEnergy()[0] - zdc->GetZNATowerEnergy()[1] - zdc->GetZNATowerEnergy()[3] - zdc->GetZNATowerEnergy()[4]; 
        } // special case for ZDC-A channel 2 for 2015 datataking
          multsum += energy;
          qvectorsum[0] += zdcX[ich-1]*energy;
          qvectorsum[1] += zdcY[ich-1]*energy; 
      }
    } else { // all other data periods
      for (int ich = 1; ich < 5; ++ich) {
        energy = zdc->GetZNATowerEnergy()[ich];
        multsum += energy;
        qvectorsum[0] += zdcX[ich-1]*energy;
        qvectorsum[1] += zdcY[ich-1]*energy; 
      }
    }
  } 
  else if (zdctype== ZDC::C) {
    for (int ich = 1; ich < 5; ++ich) {
      energy = zdc->GetZNCTowerEnergy()[ich];
      multsum += energy;
      qvectorsum[0] += zdcX[ich-1]*energy;
      qvectorsum[1] += zdcY[ich-1]*energy; 
    }
  }
  qvectorsum[0] /= multsum;
  qvectorsum[1] /= multsum;
  return qvectorsum;
}


