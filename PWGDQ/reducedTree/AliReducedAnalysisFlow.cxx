//
// Creation date: 2015/10/01
// Author: Ionut-Cristian Arsene, iarsene@cern.ch, i.c.arsene@fys.uio.no

#include "AliReducedAnalysisFlow.h"

#include <iostream>
using std::cout;
using std::endl;

#include <TClonesArray.h>
#include <TVector2.h>
#include <TList.h>

#include "AliReducedVarManager.h"
#include "AliReducedEventInfo.h"
#include "AliReducedBaseEvent.h"
#include "AliReducedTrackInfo.h"
#include "AliReducedCaloClusterInfo.h"
#include "AliReducedPairInfo.h"
#include "AliHistogramManager.h"

ClassImp(AliReducedAnalysisFlow);


//___________________________________________________________________________
AliReducedAnalysisFlow::AliReducedAnalysisFlow() :
  AliReducedAnalysisTaskSE(),
  fHistosManager(new AliHistogramManager("Histogram Manager", AliReducedVarManager::kNVars)),
  fEventCuts(),
  fTrackCuts(),
  fhVzeroQx(),
  fhVzeroQy()
{
  //
  // default constructor
  //

}


//___________________________________________________________________________
AliReducedAnalysisFlow::AliReducedAnalysisFlow(const Char_t* name, const Char_t* title) :
  AliReducedAnalysisTaskSE(name,title),
  fHistosManager(new AliHistogramManager("Histogram Manager", AliReducedVarManager::kNVars)),
  fEventCuts(),
  fTrackCuts(),
  fhVzeroQx(),
  fhVzeroQy()
{
  //
  // named constructor
  //
   fEventCuts.SetOwner(kTRUE);
   fTrackCuts.SetOwner(kTRUE);
}


//___________________________________________________________________________
AliReducedAnalysisFlow::~AliReducedAnalysisFlow()
{
  //
  // destructor
  //
   fEventCuts.Clear("C"); 
   if(fHistosManager) delete fHistosManager;
}

//___________________________________________________________________________
Bool_t AliReducedAnalysisFlow::IsEventSelected(AliReducedBaseEvent* event) {
  //
  // apply event cuts
  //
  if(fEventCuts.GetEntries()==0) return kTRUE;
  // loop over all the cuts and make a logical and between all cuts in the list
  // TODO: Cut masks or more complicated cut configurations can be implemented here
  for(Int_t i=0; i<fEventCuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*)fEventCuts.At(i);
    if(!cut->IsSelected(event)) return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________________
Bool_t AliReducedAnalysisFlow::IsTrackSelected(AliReducedBaseTrack* track) {
  //
  // apply event cuts
  //
  if(fTrackCuts.GetEntries()==0) return kTRUE;
  // loop over all the cuts and make a logical and between all cuts in the list
  // TODO: Cut masks or more complicated cut configurations can be implemented here
  for(Int_t i=0; i<fTrackCuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*)fTrackCuts.At(i);
    if(!cut->IsSelected(track)) return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________________
void AliReducedAnalysisFlow::Init() {
  //
  // initialize stuff
  //
   AliReducedVarManager::SetDefaultVarNames();
   fHistosManager->SetUseDefaultVariableNames(kTRUE);
   fHistosManager->SetDefaultVarNames(AliReducedVarManager::fgVariableNames,AliReducedVarManager::fgVariableUnits);

   fhVzeroQx = new TProfile("hVzeroQx","hVzeroQx",100,0,100);
   fhVzeroQy = new TProfile("hVzeroQy","hVzeroQy",100,0,100);
   TList *outputlist = new TList();
   outputlist->SetOwner(kTRUE);
   outputlist->SetName("Histograms");
   outputlist->Add(fhVzeroQy);
   outputlist->Add(fhVzeroQx);
   fHistosManager->AddToOutputList(outputlist);
}


//___________________________________________________________________________
void AliReducedAnalysisFlow::Process() {
  //
  // process the current event
  //
  if(!fEvent) return;

  AliReducedVarManager::SetEvent(fEvent);

  AliReducedVarManager::FillEventInfo(fEvent, fValues);
  fHistosManager->FillHistClass("Event_NoCuts", fValues);
  if(fEvent->IsA()==AliReducedEventInfo::Class()) {
    for(UShort_t ibit=0; ibit<64; ++ibit) {
      AliReducedVarManager::FillEventOnlineTrigger(ibit, fValues);
      fHistosManager->FillHistClass("OnlineTriggers_NoCuts", fValues);
    }
  }
  if(!IsEventSelected(fEvent)) return;

  fHistosManager->FillHistClass("Event_AfterCuts", fValues);
  AliReducedEventInfo* eventInfo = NULL;
  if(fEvent->IsA()==AliReducedEventInfo::Class()) eventInfo = (AliReducedEventInfo*)fEvent;

  if (eventInfo) {
    TVector2 q;
    TVector2 qsum;
    Double_t vzeroChannelPhi[8] = {0.3927, 1.1781, 1.9635, 2.7489, -2.7489, -1.9635, -1.1781, -0.3927};
    for (Int_t ich = 0; ich < 64; ++ich) {
      Double_t vzeroMult = eventInfo->MultChannelVZERO(ich);
      q.Set(vzeroMult*TMath::Sin(vzeroChannelPhi[ich%8]),vzeroMult*TMath::Cos(vzeroChannelPhi[ich%8]));
      qsum += q;
    }
    fhVzeroQx->Fill(fEvent->CentralityTPC(),qsum.X());
    fhVzeroQy->Fill(fEvent->CentralityTPC(),qsum.Y());
    fHistosManager->FillHistClass("Flow", fValues);
  }
  
}
  //
  // if(eventInfo) {
  //   for(UShort_t ibit=0; ibit<64; ++ibit) {
  //     AliReducedVarManager::FillEventOnlineTrigger(ibit, fValues);
  //     fHistosManager->FillHistClass("OnlineTriggers_AfterCuts", fValues);
  //     for(UShort_t i=0; i<32; ++i) {
  //       AliReducedVarManager::FillL0TriggerInputs(eventInfo, i, fValues);
  //       fHistosManager->FillHistClass("OnlineTriggers_vs_L0TrigInputs", fValues);
  //     }
  //     for(UShort_t i=0; i<32; ++i) {
  //       AliReducedVarManager::FillL1TriggerInputs(eventInfo, i, fValues);
  //       fHistosManager->FillHistClass("OnlineTriggers_vs_L1TrigInputs", fValues);
  //     }
  //     for(UShort_t i=0; i<16; ++i) {
  //       AliReducedVarManager::FillL2TriggerInputs(eventInfo, i, fValues);
  //       fHistosManager->FillHistClass("OnlineTriggers_vs_L2TrigInputs", fValues);
  //     }
  //   }
  // }
  //
  // for(UShort_t ibit=0; ibit<64; ++ibit) {
  //   AliReducedVarManager::FillEventTagInput(fEvent, ibit, fValues);
  //   fHistosManager->FillHistClass("EvtTags", fValues);
  // }
  //
  // if(eventInfo) {
  //   for(UShort_t ibit=0; ibit<32; ++ibit) {
  //     AliReducedVarManager::FillL0TriggerInputs(eventInfo, ibit, fValues);
  //     fHistosManager->FillHistClass("L0TriggerInput", fValues);
  //   }
  //   for(UShort_t ibit=0; ibit<32; ++ibit) {
  //     AliReducedVarManager::FillL1TriggerInputs(eventInfo, ibit, fValues);
  //     fHistosManager->FillHistClass("L1TriggerInput", fValues);
  //   }
  //   for(UShort_t ibit=0; ibit<16; ++ibit) {
  //     AliReducedVarManager::FillL2TriggerInputs(eventInfo, ibit, fValues);
  //     fHistosManager->FillHistClass("L2TriggerInput", fValues);
  //   }
  //
  //   for(Int_t icl=0; icl<eventInfo->GetNCaloClusters(); ++icl) {
  //     AliReducedVarManager::FillCaloClusterInfo(eventInfo->GetCaloCluster(icl), fValues);
  //     fHistosManager->FillHistClass("CaloClusters", fValues);
  //   }
  // }
//
//   AliReducedBaseTrack* track = 0x0;
//   TClonesArray* trackList = fEvent->GetTracks();
//   TIter nextTrack(trackList);
//   if(trackList) {
//     for(Int_t it=0; it<fEvent->NTracks(); ++it) {
//       track = (AliReducedBaseTrack*)nextTrack();
//
//       if(!IsTrackSelected(track)) continue;
//       AliReducedVarManager::FillTrackInfo(track,fValues);
//       fHistosManager->FillHistClass("TrackQA_AllTracks_ITS_TPC_TRD_TOF_EMCAL", fValues);
//
//       AliReducedTrackInfo* trackInfo = NULL;
//       if(track->IsA()==AliReducedTrackInfo::Class()) trackInfo = (AliReducedTrackInfo*)track;
//
//       if(trackInfo) {
//         for(UInt_t iflag=0; iflag<AliReducedVarManager::kNTrackingFlags; ++iflag) {
//           AliReducedVarManager::FillTrackingFlag(trackInfo, iflag, fValues);
//           fHistosManager->FillHistClass("TrackingFlags", fValues);
//         }
//       }
//       for(UShort_t iflag=0; iflag<64; ++iflag) {
//         AliReducedVarManager::FillTrackQualityFlag(track, iflag, fValues);
//         fHistosManager->FillHistClass("TrackQualityFlags", fValues);
//       }
//       if(trackInfo) {
//         for(Int_t iLayer=0; iLayer<6; ++iLayer) {
//           AliReducedVarManager::FillITSlayerFlag(trackInfo, iLayer, fValues);
//           fHistosManager->FillHistClass("ITSclusterMap", fValues);
//         }
//         for(Int_t iLayer=0; iLayer<8; ++iLayer) {
//           AliReducedVarManager::FillTPCclusterBitFlag(trackInfo, iLayer, fValues);
//           fHistosManager->FillHistClass("TPCclusterMap", fValues);
//         }
//       }
//       if(track->IsGammaLeg()) fHistosManager->FillHistClass("TrackQA_GammaLeg_ITS_TPC_TRD_TOF_EMCAL", fValues);
//       if(track->IsPureGammaLeg()) fHistosManager->FillHistClass("TrackQA_PureGammaLeg_ITS_TPC_TRD_TOF_EMCAL", fValues);
//       if(track->IsK0sLeg()) fHistosManager->FillHistClass("TrackQA_K0sLeg_ITS_TPC_TRD_TOF_EMCAL", fValues);
//       if(track->IsPureK0sLeg()) fHistosManager->FillHistClass("TrackQA_PureK0sLeg_ITS_TPC_TRD_TOF_EMCAL", fValues);
//       if(track->IsLambdaLeg()) {
//         if(track->Charge()>0) fHistosManager->FillHistClass("TrackQA_LambdaPosLeg_ITS_TPC_TRD_TOF_EMCAL", fValues);
//         else fHistosManager->FillHistClass("TrackQA_LambdaNegLeg_ITS_TPC_TRD_TOF_EMCAL", fValues);
//       }
//       if(track->IsPureLambdaLeg()) {
//         if(track->Charge()>0) fHistosManager->FillHistClass("TrackQA_PureLambdaPosLeg_ITS_TPC_TRD_TOF_EMCAL", fValues);
//         else fHistosManager->FillHistClass("TrackQA_PureLambdaNegLeg_ITS_TPC_TRD_TOF_EMCAL", fValues);
//       }
//       if(track->IsALambdaLeg()) {
//         if(track->Charge()>0) fHistosManager->FillHistClass("TrackQA_ALambdaPosLeg_ITS_TPC_TRD_TOF_EMCAL", fValues);
//         else fHistosManager->FillHistClass("TrackQA_ALambdaNegLeg_ITS_TPC_TRD_TOF_EMCAL", fValues);
//       }
//       if(track->IsPureALambdaLeg()) {
//         if(track->Charge()>0) fHistosManager->FillHistClass("TrackQA_PureALambdaPosLeg_ITS_TPC_TRD_TOF_EMCAL", fValues);
//         else fHistosManager->FillHistClass("TrackQA_PureALambdaNegLeg_ITS_TPC_TRD_TOF_EMCAL", fValues);
//       }
//     }  // end loop over tracks
//   }  // end if(trackList)
//   fHistosManager->FillHistClass("Event_AfterCuts", fValues);
// }


//___________________________________________________________________________
void AliReducedAnalysisFlow::Finish() {
  //
  // run stuff after the event loop
  //
}
