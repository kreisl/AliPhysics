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
  fRequireTime(false),
  fQnLimits(nullptr),
  fTimeAxis(nullptr),
  fFilterBit(768),
  fCalibFileType(CalibFile::local),
  fQnManager(nullptr),
  fEvent(nullptr),
  fInDetectorList(nullptr),
  fQnTree(nullptr),
  fInCalib(nullptr),
  fOCDBPath(),
  fOCDBAvailable(),
  fLHCData(nullptr),
  fEventNumber(0) { }

AliAnalysisTaskQn::AliAnalysisTaskQn(const char *name) : AliAnalysisTaskSE(name),
  fRequireTime(false),
  fQnLimits(nullptr),
  fTimeAxis(nullptr),
  fFilterBit(768),
  fCalibFileType(CalibFile::local),
  fQnManager(new Qn::CorrectionManager()),
  fEvent(nullptr),
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
  DefineOutput(OutputSlot::QnTree,TTree::Class());
}

AliAnalysisTaskQn::~AliAnalysisTaskQn() {
}

void AliAnalysisTaskQn::UserCreateOutputObjects() {
  fQnTree = new TTree("tree","tree");
  fQnManager->ConnectOutputTree(fQnTree);
  fQnManager->SetCalibrationInputFile(fInCalib);
  fQnManager->InitializeOnNode();
  if (fInCalib) {
    fInDetectorList = (TList*) fInCalib->Get("CorrectionQAHistograms");
    if (fInDetectorList) fQnLimits = (AliQnLimits*) fInDetectorList->FindObject("TimeStampLimits"); 
  }
  auto qalist = fQnManager->GetCorrectionQAList();
  if (fRequireTime) {
    if (fQnLimits) {
      auto min = fQnLimits->Min();
      auto max = fQnLimits->Max();
      std::cout << std::fixed << "Setting the TimeStamp axis to " << min << " - " << max << std::endl;
      std::cout << std::scientific;
      fTimeAxis = new Qn::Axis<double>("TimeStamp",30,min,max);
      qalist->Add(fQnLimits);
    } else {
      fQnLimits = new AliQnLimits("TimeStampLimits");
      std::cout << "Adding limits to the event qa list" << std::endl;
      qalist->Add(fQnLimits);
    }
  }
  PostData(OutputSlot::QnCalibration,fQnManager->GetCorrectionList());
  PostData(OutputSlot::QnQA,qalist);
  PostData(OutputSlot::QnTree,fQnTree);
}

void AliAnalysisTaskQn::NotifyRun() {
  AliInfo(TString::Format("New run number: %d", this->fCurrentRunNumber).Data());
  fQnManager->SetCurrentRunName(TString::Format("%d", this->fCurrentRunNumber).Data());
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
  fEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fEvent) return;
  auto eventhandler = dynamic_cast<AliInputEventHandler*>(
      AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  auto trigger = eventhandler->IsEventSelected();
  if(!trigger) return;
  fValues = fQnManager->GetVariableContainer();
  // begin vertex selection
  // as of Alexandru Dobrin 8.5.19
    const AliAODVertex* vtxtrk = fEvent->GetPrimaryVertex();
  if (fEvent->GetRunNumber() < 200000) {
    if (vtxtrk->GetZ() > 10 || vtxtrk->GetZ() < -10) return;
    //const AliAODVertex* vtxspd = fEvent->GetPrimaryVertexSPD();
    //double covtrk[6], covspd[6];
    //vtxtrk->GetCovarianceMatrix(covtrk);
    //vtxspd->GetCovarianceMatrix(covspd);
    //double dz = vtxtrk->GetZ() - vtxspd->GetZ();
    //double errortotal = sqrt(covtrk[5]+covspd[5]);
    //double errortracks = sqrt(covtrk[5]);
    //double nstot = dz/errortotal; // n sigmas of deviation with total error
    //double nstrk = dz/errortracks; // n sigmas of deviation with track error
    //if (std::abs(dz) > 0.2 || std::abs(nstot) > 10 || std::abs(nstrk) > 20) return;
    // end vertex selection
  }
  fValues[kVtxX] = vtxtrk->GetX();
  fValues[kVtxY] = vtxtrk->GetY();
  fValues[kVtxZ] = vtxtrk->GetZ();
  bool old_centrality = false;
  if (fEvent->GetRunNumber() < 200000) old_centrality = true; // true if run 1
  if (!old_centrality) { // run 2
    auto mult = (AliMultSelection*)fEvent->FindListObject("MultSelection");
    if (mult) {
      fValues[kCentV0A] = mult->GetMultiplicityPercentile("V0A");
      fValues[kCentV0C] = mult->GetMultiplicityPercentile("V0C");
      fValues[kCentV0M] = mult->GetMultiplicityPercentile("V0M");
      fValues[kCentZNC] = mult->GetMultiplicityPercentile("ZNC");
      fValues[kCentZNA] = mult->GetMultiplicityPercentile("ZNA");
      fValues[kCentCL0] = mult->GetMultiplicityPercentile("CL0");
      fValues[kCentCL1] = mult->GetMultiplicityPercentile("CL1");
    }
  }
  if (old_centrality) {
    auto centrality = fEvent->GetCentrality(); // run 1
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
  // Trigger selection
  // Events marked as kINT7 will not be tagged as kCentral or kSemiCentral
  auto ismb = trigger & (AliVEvent::kMB+AliVEvent::kINT7);
  if (ismb) fValues[kTrigger] = 0;
  if (trigger & AliVEvent::kCentral && !ismb) fValues[kTrigger] = 1;
  if (trigger & AliVEvent::kSemiCentral && !ismb) fValues[kTrigger] = 2;
  if (fRequireTime) {
    if (fTimeAxis) {
      fValues[kTimeStamp] = fTimeAxis->FindBin(fEvent->GetTimeStamp()); 
    } else {
      fQnLimits->SetNew(fEvent->GetTimeStamp());
    }
  } else {
    fValues[kTimeStamp] = -1;
  }
	fValues[kRunNumber] = fEvent->GetRunNumber();
	fValues[kEventNumber] = fEventNumber;
	++fEventNumber;
  fValues[kPeriodNumber] = fEvent->GetPeriodNumber();
  fValues[kOrbitNumber] = fEvent->GetOrbitNumber();
  fValues[kBunchCrossNumber] = fEvent->GetBunchCrossNumber();
  fValues[kNTracklets] = fEvent->GetMultiplicity()->GetNumberOfTracklets();
  for (int ilayer = 0; ilayer < 6; ++ilayer) {
    fValues[kNITSClusters + ilayer] = fEvent->GetNumberOfITSClusters(ilayer);
  }
  fValues[kNTPCClusters] = fEvent->GetNumberOfTPCClusters();

  // |--------------ZNC-----------------|
  // | channel       | coordinatesystem |
  // | configuration |                  |
  // |  _________    | C-Side    A-Side |
  // |  | 1 | 2 |    |  |y        |y    |
  // |  |___|___|    |  |__x      |__-x |
  // |  | 3 | 4 |    | /z        /-z    |
  // |  |_______|    |                  |
  // |     angle from 0 to 2 pi         |
  // |----------------------------------|
  auto zdc = fEvent->GetZDCData();
  double summult = 0.;
  double znsumenergy = 0.;
  const std::array<double, 5> phizna{0, 5.4977871, 3.9269908, 0.78539816, 2.3561945};
  const std::array<double, 5> phiznc{0, 3.9269908, 5.4977871, 2.3561945, 0.78539816}; 
  // fill ZNC
  for (int ich = 0; ich < 5; ++ich) {
    fValues[kZDCCChMult + ich] = zdc->GetZNCTowerEnergy()[ich];
    if (fValues[kZDCCChMult + ich] < 0)  fValues[kZDCCChMult + ich] = 0.f;
    fValues[kZDCCChPhi + ich] = phiznc[ich];
    summult += fValues[kZDCCChMult+ich];
  }
  fValues[kZDCCSumMult] = summult;
  znsumenergy += summult;
  // fill ZNA
  // special case for LHC15o sample where channel 2 was offline
  // it is reconstructed from the other channels using the sum channel
  auto tow = zdc->GetZNATowerEnergy();
  if (fEvent->GetRunNumber() >= 240000 && fEvent->GetRunNumber() <= 248000) {
    for (int ich = 0; ich < 5; ++ich) {
      fValues[kZDCAChPhi + ich] = phizna[ich];
      if (ich!= 2) fValues[kZDCAChMult + ich] = tow[ich];
      if (ich == 2) {
        fValues[kZDCAChMult + ich] = tow[0] - tow[1] - tow[3] - tow[4];
      }
      if (fValues[kZDCAChMult + ich] < 0)  fValues[kZDCAChMult + ich] = 0.f;
    }
    fValues[kZDCASumMult] = -1.;
  } else { // all other data periods
    summult = 0.;
    for (int ich = 0; ich < 5; ++ich) {
      fValues[kZDCAChPhi + ich] = phizna[ich];
      fValues[kZDCAChMult + ich] = zdc->GetZNATowerEnergy()[ich];
      if (fValues[kZDCAChMult + ich] < 0)  fValues[kZDCAChMult + ich] = 0.f;
      summult += fValues[kZDCAChMult+ich];
    }
  fValues[kZDCASumMult] = summult;
  znsumenergy += summult;
  }
  fValues[kZNSumEnergy] = znsumenergy;
  //ZP
  // channels A-side: (BEAM) |4 |3 |2 |1 |
  // channels C-side: (BEAM) |1 |2 |3 |4 |
  constexpr std::array<double, 5> scalingfactor_zpc{1,-1,-1,-1};
  constexpr std::array<double, 5> scalingfactor_zpa{-1,-1,-1,1};
  for (int ich = 0; ich < 5; ++ich) {
    fValues[kZPAChMult + ich] = zdc->GetZPATowerEnergy()[ich];
    fValues[kZPCChMult + ich] = zdc->GetZPCTowerEnergy()[ich];
    fValues[kZPAChOffset + ich] = scalingfactor_zpa[ich];
    fValues[kZPCChOffset + ich] = scalingfactor_zpc[ich];
    fValues[kZPAChPhi + ich] = 0.;
    fValues[kZPCChPhi + ich] = 3.141592;
  }

  // V0
  double vzeromult = 0.;
  auto vzero = fEvent->GetVZEROData();
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
  if (fQnManager->ProcessEvent()) {
    fQnManager->FillChannelDetectors();
    unsigned int ntracksfilterbit = 0;
    unsigned int ntracks = fEvent->GetNumberOfTracks();
    for (unsigned int i = 0; i < ntracks; ++i) {
      auto track = static_cast<AliAODTrack*>(fEvent->GetTrack(i));
      if (!track) continue;
      if (!track->TestFilterBit(fFilterBit)) continue;
      ++ntracksfilterbit;
      fValues[kTPCnCls] = track->GetTPCNcls();
      fValues[kTPCchi2pCls] = track->GetTPCchi2perCluster();
      fValues[kPhi] = track->Phi();
      fValues[kPt] = track->Pt();
      fValues[kEta] = track->Eta();
      fValues[kCharge] = track->GetSign();
      fQnManager->FillTrackingDetectors();
    }
    fValues[kNTPCTracks] = ntracksfilterbit;
    if (fRequireTime) {
      if (fTimeAxis) fQnManager->ProcessCorrections();
    } else {
      fQnManager->ProcessCorrections();
    }
  }
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
        std::cout << "setting calibration file to: " 
          << calibfile->GetName() << std::endl;
      }
      break;
    default:
      std::cout <<"Calibration file not supported" << std::endl;
      break;
  }
}
