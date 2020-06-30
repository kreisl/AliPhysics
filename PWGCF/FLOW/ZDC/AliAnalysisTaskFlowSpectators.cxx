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

#include "AliAnalysisTaskFlowSpectators.h"

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TList.h>

#include <map>

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODTracklets.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliInputEventHandler.h"
#include "AliMultSelection.h"
#include "AliOADBContainer.h"
#include "AliQvectors.h"
#include "AliZDCflowCuts.h"

AliAnalysisTaskFlowSpectators::AliAnalysisTaskFlowSpectators() : AliAnalysisTaskSE() {}

AliAnalysisTaskFlowSpectators::AliAnalysisTaskFlowSpectators(const char *name)
    : AliAnalysisTaskSE(name) {
  DefineOutput(1, AliZDCResultStorage::Class());
  DefineOutput(2, TTree::Class());
}

AliAnalysisTaskFlowSpectators::~AliAnalysisTaskFlowSpectators() {
  delete fCorrelationList;
  delete fCorrectionList;
  delete fQAList;
  delete fTree;
  delete fAnalysisUtils;
  delete fCorrectionManager;
}

void AliAnalysisTaskFlowSpectators::UserCreateOutputObjects() {
  // QnTools
  fTree = new TTree("tree", "tree");
  if (fActivateQnTools) {
    fCorrectionManager->ConnectOutputTree(fTree);
    auto file = TFile::Open(fCorrectionFileName.data());
    fCorrectionManager->SetCalibrationInputFile(file);
    fCorrectionManager->InitializeOnNode();
  }

  fCorrelationList = new TList();
  fCorrelationList->SetName("Correlations");
  fCorrelationList->SetOwner(true);
  fCorrectionList = new TList();
  fCorrectionList->SetName("Corrections");
  fCorrectionList->SetOwner(true);
  fQAList = new TList();
  fQAList->SetName("QA");
  fQAList->SetOwner(true);

  // Apply delayed configuration
  ApplyConfiguration();

  fGainEqualizationZNA.Configure("ZNA", 4, {4});
  fGainEqualizationZNC.Configure("ZNC", 4, {4});

  fGainEqualizationZNA.AddCorrectionsToList(fCorrectionList, fQAList);
  fGainEqualizationZNC.AddCorrectionsToList(fCorrectionList, fQAList);

  auto event_cut_list = new TList();
  event_cut_list->SetOwner(true);
  event_cut_list->SetName("event_cuts");
  fEventCuts.AddQAplotsToList(event_cut_list);
  fQAList->Add(event_cut_list);

  for (auto &analysis : fEllipticFlowAnalyses) analysis.SetBootStrapSamples(fNsamples);
  for (auto &analysis : fCumulantFlowAnalyses) analysis.SetBootStrapSamples(fNsamples);

  auto elliptic_list = new TList();
  elliptic_list->SetOwner(true);
  elliptic_list->SetName("elliptic_flow");
  auto qa_elliptic_list = new TList();
  qa_elliptic_list->SetOwner(true);
  qa_elliptic_list->SetName("elliptic_flow_qa");
  for (auto &analysis : fEllipticFlowAnalyses) {
    elliptic_list->Add(analysis.CreateCorrelations());
    analysis.AddQAHistograms(qa_elliptic_list);
  }
  fQAList->Add(qa_elliptic_list);
  fCorrelationList->Add(elliptic_list);

  auto cumulants_list = new TList();
  cumulants_list->SetOwner(true);
  cumulants_list->SetName("cumulants_flow");
  auto cumulants_qa_list = new TList();
  cumulants_qa_list->SetOwner(true);
  cumulants_qa_list->SetName("cumulants_flow_qa");
  for (auto &analysis : fCumulantFlowAnalyses) {
    cumulants_list->Add(analysis.CreateCorrelations());
    analysis.AddCorrectionsToList(fCorrectionList, cumulants_qa_list);
  }
  fQAList->Add(cumulants_qa_list);
  fCorrelationList->Add(cumulants_list);

  // event QA plots
  fCentralityV0M = new TH1D("CentV0M", ";centrality V0M;N", 100, 0., 100.);
  fCentralityCL1 = new TH1D("CentCL1", ";centrality CL1;N", 100, 0., 100.);
  fCentralityCL1vsV0M = new TH2D("CentCL1vsV0M", ";centrality CL1;centrality V0M;N", 100, 0., 100., 100, 0., 100.);
  fQAList->Add(fCentralityV0M);
  fQAList->Add(fCentralityCL1);
  fQAList->Add(fCentralityCL1vsV0M);
  fVertexX = new TH1D("VertexX", ";vertex x / cm;N", 50, -0.03, 0.02);
  fVertexY = new TH1D("VertexY", ";vertex y / cm;N", 50, 0.14, 0.21);
  fVertexZ = new TH1D("VertexZ", ";vertex z / cm;N", 50, -10., 10.);
  fVertexXY = new TH2D("VertexXVertexY", ";vertex x; vertex y;N", 50, -0.03, 0.02, 50, 0.14, 0.21);
  fQAList->Add(fVertexX);
  fQAList->Add(fVertexY);
  fQAList->Add(fVertexZ);
  fQAList->Add(fVertexXY);

  fMultiplicityTPConlyVsGlobal =
      new TH2D("MultiplicityTPConlyvsGlobal", ";global;tpc only", 1500, 0., 3000., 2000, 0., 4000.);
  fMultiplicityTPConlyVsGlobalCut =
      new TH2D("MultiplicityTPConlyvsGlobalCut", ";global;tpc only", 1500, 0., 3000., 2000, 0., 4000.);
  fNSigmaTPConlyVsGlobal =
      new TH2D("fNSigmaTPConlyVsGlobal", ";global;N sigma ;tpc only", 1500, 0., 3000., 500, -10., 10.);
  fNtracksESDvsNclsITS =
      new TH2D("NtracksESDvsNclsITS", ";esd tracks;ITS clusters", 1500, 0., 14000., 2000, 0., 26000.);
  fQAList->Add(fMultiplicityTPConlyVsGlobal);
  fQAList->Add(fMultiplicityTPConlyVsGlobalCut);
  fQAList->Add(fNSigmaTPConlyVsGlobal);
  fQAList->Add(fNtracksESDvsNclsITS);

  fPsiZA = new TH1D("Psi_ZNA", ";#Psi_{ZNA};N", 50, 0., 2 * TMath::Pi());
  fPsiZC = new TH1D("Psi_ZNC", ";#Psi_{ZNC};N", 50, 0., 2 * TMath::Pi());
  fPsiZAEQ = new TH1D("Psi_ZNA_EQ", ";#Psi_{ZNA};N", 50, 0., 2 * TMath::Pi());
  fPsiZCEQ = new TH1D("Psi_ZNC_EQ", ";#Psi_{ZNC};N", 50, 0., 2 * TMath::Pi());
  fQAList->Add(fPsiZA);
  fQAList->Add(fPsiZC);
  fQAList->Add(fPsiZAEQ);
  fQAList->Add(fPsiZCEQ);
  std::string correctionstepname("CorrectionStep");
  fCorrectionStep = new TH1D(correctionstepname.data(), "current correction step;step;", 4, 0., 4.);
  fCorrectionStepEQ =
      new TH1D((correctionstepname + "_EQ").data(), "current correction step after eq;step;", 4, 0., 4.);
  fCorrectionList->Add(fCorrectionStep);
  fCorrectionList->Add(fCorrectionStepEQ);

  fQZNAmagnitude.Configure("QZNAmagnitude", 1000, 180.);
  fQZNCmagnitude.Configure("QZNCmagnitude", 1000, 180.);
  fQZNAmagnitude.AddToOutputList(fCorrectionList);
  fQZNCmagnitude.AddToOutputList(fCorrectionList);

  fAnalysisUtils = new AliAnalysisUtils();
  fAnalysisUtils->SetUseMVPlpSelection(true);
  fAnalysisUtils->SetUseOutOfBunchPileUp(true);

  fEventCuts.SetupRun1PbPb();

  fOutputList = new TList();
  fOutputList->SetOwner(true);
  fOutputList->Add(fCorrelationList);
  fOutputList->Add(fCorrectionList);
  fOutputList->Add(fQAList);
  if (fActivateQnTools) fOutputList->Add(fCorrectionManager->GetCorrectionList());
  fStorage = new AliZDCResultStorage("ZDC", fOutputList);
  PostData(1, fStorage);
  PostData(2, fTree);
}

// Reads corrections from the OADB file for each new run.
// Logs the status of corrections in the fCorrectionStep histogram.
void AliAnalysisTaskFlowSpectators::NotifyRun() {
  // QnTools
  if (fActivateQnTools) {
    AliInfo(TString::Format("New run number: %d", this->fCurrentRunNumber).Data());
    fCorrectionManager->SetCurrentRunName(TString::Format("%d", this->fCurrentRunNumber).Data());
  }

  // Open Corrections
  auto file = TFile::Open(fCorrectionFileName.data());
  if (file && !file->IsZombie()) {
    fQZNAmagnitude.ReadFile(file);
    fQZNCmagnitude.ReadFile(file);

    // Multiplicity Cut
    auto graph3sigmap = dynamic_cast<TGraph *>(file->Get("TPConlyVsGlobalMeanPlusSigma3"));
    auto graph3sigmam = dynamic_cast<TGraph *>(file->Get("TPConlyVsGlobalMeanMinusSigma3"));
    auto graphmean = dynamic_cast<TGraph *>(file->Get("TPConlyVsGlobalMean"));
    if (graph3sigmap && graph3sigmam && graphmean) {
      std::cout << "multiplicity cut found" << std::endl;
      fMult3SigmaPlus = graph3sigmap;
      fMult3SigmaMinus = graph3sigmam;
      fMultMean = graphmean;
    }

    // Centrality Weight
    fCentralityWeightInput = dynamic_cast<TH1D *>(file->Get("CentralityWeights"));
    if (fCentralityWeightInput) {
      fCentralityWeightInput->SetDirectory(0);
      std::cout << "using Centrality weights" << std::endl;
    }

    // nd
    fAlignZNA.Make(file, fCurrentRunNumber, fCorrectionList, fQAList);
    fAlignZNC.Make(file, fCurrentRunNumber, fCorrectionList, fQAList);

    // Gain Equalization
    fGainEqualizationZNA.OpenCorrection(file, fCurrentRunNumber);
    fGainEqualizationZNC.OpenCorrection(file, fCurrentRunNumber);

    // Recenter V0
    fRecenter4DV0A.Make(file, fCurrentRunNumber, fCorrectionList, fQAList);
    fRecenter4DV0C.Make(file, fCurrentRunNumber, fCorrectionList, fQAList);

    // Recenter ZN
    fRecenter4DZNA.Make(file, fCurrentRunNumber, fCorrectionList, fQAList);
    fRecenter4DZNC.Make(file, fCurrentRunNumber, fCorrectionList, fQAList);

    // Recenter ZN
    fRecenter4DAfterGainEqZNA.Make(file, fCurrentRunNumber, fCorrectionList, fQAList);
    fRecenter4DAfterGainEqZNC.Make(file, fCurrentRunNumber, fCorrectionList, fQAList);

    // cumulants
    for (auto &analysis : fCumulantFlowAnalyses) {
      analysis.OpenCorrection(file, fCurrentRunNumber);
      analysis.OpenPtEfficiencies(file);
    }
    file->Close();
  }
  PostData(1, fStorage);
}

void AliAnalysisTaskFlowSpectators::UserExec(Option_t *) {
  if (fActivateQnTools) fCorrectionManager->Reset();
  auto event = dynamic_cast<AliAODEvent *>(InputEvent());
  if (!event) return;
  auto eventhandler =
      dynamic_cast<AliInputEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!eventhandler->IsEventSelected()) return;
  if (!fEventCuts.AcceptEvent(event)) return;
  if (fAnalysisUtils->IsPileUpEvent(event)) return;
  if (MultiplicityCut(event)) {
    AnalyzeEvent(event);
  }
  fTree->Fill();
  PostData(1, fStorage);
  PostData(2, fTree);
}

bool AliAnalysisTaskFlowSpectators::MultiplicityCut(AliAODEvent *event) {
  bool pass = true;
  double ntracks_tpc_only = 0.;
  double ntracks_global = 0.;
  const unsigned int ntracks = event->GetNumberOfTracks();
  for (unsigned int i = 0; i < ntracks; ++i) {
    auto track = static_cast<AliAODTrack *>(event->GetTrack(i));
    if (!track) continue;
    if (track->TestFilterBit(128)) ++ntracks_tpc_only;
    if (track->TestFilterBit(768)) ++ntracks_global;
  }
  double sum_its_cls = 0.;
  for (int i = 2; i < 6; ++i) {
    sum_its_cls += event->GetNumberOfITSClusters(i);
  }
  fNtracksESDvsNclsITS->Fill(event->GetNumberOfESDTracks(), sum_its_cls);
  fMultiplicityTPConlyVsGlobal->Fill(ntracks_global, ntracks_tpc_only);
  if (fMultMean && fMult3SigmaPlus && fMult3SigmaMinus) {
    auto max = fMult3SigmaPlus->Eval(ntracks_global);
    auto min = fMult3SigmaMinus->Eval(ntracks_global);
    if (ntracks_tpc_only < min || ntracks_tpc_only > max) pass = false;
    auto mean = fMultMean->Eval(ntracks_global);
    auto nsigma = 3. * (ntracks_tpc_only - mean) / (max - mean);
    fNSigmaTPConlyVsGlobal->Fill(ntracks_global, nsigma);
  }
  if (pass) {
    fMultiplicityTPConlyVsGlobalCut->Fill(ntracks_global, ntracks_tpc_only);
  }
  return pass;
}

void AliAnalysisTaskFlowSpectators::AnalyzeEvent(AliAODEvent *event) {
  // QnTools
  if (fActivateQnTools) {
    fValues = fCorrectionManager->GetVariableContainer();
  }
  bool run_1 = false;
  if (event->GetRunNumber() < 200000) run_1 = true;
  float cent_v0m;
  float cent_cl1;
  if (run_1) {
    // Run 1
    auto centrality = event->GetCentrality();
    if (centrality) {
      if (fActivateQnTools) {
        fValues[kCentV0A] = centrality->GetCentralityPercentile("V0A");
        fValues[kCentV0C] = centrality->GetCentralityPercentile("V0C");
        fValues[kCentV0M] = centrality->GetCentralityPercentile("V0M");
        fValues[kCentZNC] = centrality->GetCentralityPercentile("ZNC");
        fValues[kCentZNA] = centrality->GetCentralityPercentile("ZNA");
        fValues[kCentCL0] = centrality->GetCentralityPercentile("CL0");
        fValues[kCentCL1] = centrality->GetCentralityPercentile("CL1");
      }
      cent_v0m = centrality->GetCentralityPercentile("V0M");
      cent_cl1 = centrality->GetCentralityPercentile("CL1");
    }
  } else {
    // Run 2
    auto centrality = dynamic_cast<AliMultSelection *>(event->FindListObject("MultSelection"));
    if (centrality) {
      if (fActivateQnTools) {
        fValues[kCentV0A] = centrality->GetMultiplicityPercentile("V0A");
        fValues[kCentV0C] = centrality->GetMultiplicityPercentile("V0C");
        fValues[kCentV0M] = centrality->GetMultiplicityPercentile("V0M");
        fValues[kCentZNC] = centrality->GetMultiplicityPercentile("ZNC");
        fValues[kCentZNA] = centrality->GetMultiplicityPercentile("ZNA");
        fValues[kCentCL0] = centrality->GetMultiplicityPercentile("CL0");
        fValues[kCentCL1] = centrality->GetMultiplicityPercentile("CL1");
      }
      cent_v0m = centrality->GetMultiplicityPercentile("V0M");
      cent_cl1 = centrality->GetMultiplicityPercentile("CL1");
    }
  }
  if (cent_v0m < 0. || cent_v0m > 90.) return;

  const AliAODVertex *vtx = event->GetPrimaryVertex();
  const auto vtxx = vtx->GetX();
  const auto vtxy = vtx->GetY();
  const auto vtxz = vtx->GetZ();
  if (std::fabs(vtxz) > 10.) return;
  if (fActivateQnTools) {
    fValues[kVtxX] = vtxx;
    fValues[kVtxY] = vtxy;
    fValues[kVtxZ] = vtxz;
    for (int i = 0; i < 6; ++i) fValues[kNITSClusters + i] = event->GetNumberOfITSClusters(i);
    fValues[kRunNumber] = event->GetRunNumber();
    fValues[kPeriodNumber] = event->GetPeriodNumber();
    fValues[kOrbitNumber] = event->GetOrbitNumber();
    fValues[kBunchCrossNumber] = event->GetBunchCrossNumber();
  }
  fVertexX->Fill(vtxx);
  fVertexY->Fill(vtxy);
  fVertexZ->Fill(vtxz);
  fVertexXY->Fill(vtxx, vtxy);

  fCentralityV0M->Fill(cent_v0m);
  fCentralityCL1->Fill(cent_cl1);
  fCentralityCL1vsV0M->Fill(cent_cl1, cent_v0m);
  // Create Bootstrap samples using Poisson
  GetSamples();
  // Set Correction variables for the crosscheck corrections
  std::vector<double> variables = {cent_v0m, vtxx, vtxy, vtxz};

  // Find the centrality weight
  double centrality_weight = 1.;
  if (fCentralityWeightInput) {
    centrality_weight = fCentralityWeightInput->GetBinContent(fCentralityWeightInput->FindBin(cent_v0m));
  }
  for (auto &analysis : fEllipticFlowAnalyses) {
    analysis.FindCentralityBin(event, cent_v0m, fSamples);
    analysis.SetCentralityWeight(centrality_weight);
  }
  for (auto &analysis : fCumulantFlowAnalyses) {
    analysis.FindCentralityBin(event, cent_v0m, fSamples);
    analysis.SetCentralityWeight(centrality_weight);
  }

  // Build Q-vector of ZNA and ZNC
  const auto zdc = event->GetZDCData();
  const std::array<double, 5> phizna{0, 5.4977871, 3.9269908, 0.78539816, 2.3561945};
  const std::array<double, 5> phiznc{0, 3.9269908, 5.4977871, 2.3561945, 0.78539816};
  auto tower_energy_a = zdc->GetZNATowerEnergy();
  auto tower_energy_c = zdc->GetZNCTowerEnergy();
  AliZDCQvectors q_zn_plain;
  AliZDCQvectors q_zn_equalized;
  // Build plain vector
  const auto min_energy = 1e-6;
  for (int ich = 0; ich < 5; ++ich) {
    // QnTools
    if (fActivateQnTools) {
      fValues[kZDCAChMult + ich] = tower_energy_a[ich];
      fValues[kZDCCChMult + ich] = tower_energy_c[ich];
      fValues[kZDCAChPhi + ich] = phizna[ich];
      fValues[kZDCCChPhi + ich] = phiznc[ich];
    }
    if (ich > 0) {
      if (tower_energy_a[ich] > min_energy) q_zn_plain.a.Update(phizna[ich], tower_energy_a[ich]);
      if (tower_energy_c[ich] > min_energy) q_zn_plain.c.Update(phiznc[ich], tower_energy_c[ich]);
    }
  }
  // apply gain equalization
  q_zn_equalized.a = fGainEqualizationZNA.ApplyGainEq(tower_energy_a + 1, phizna.data() + 1);
  q_zn_equalized.c = fGainEqualizationZNC.ApplyGainEq(tower_energy_c + 1, phiznc.data() + 1);
  // Normalize Q-vector for corrections
  const auto min_sum = 1e-3;
  if (q_zn_plain.a.sum < min_sum || q_zn_plain.c.sum < min_sum) return;
  q_zn_equalized.Normalize();
  q_zn_plain.Normalize();
  // Add ZDC vectors for different correction steps
  std::map<std::string, AliZDCQvectors> qvectors_zdc;
  // plain
  qvectors_zdc.emplace("plain", q_zn_plain);
  if (fGainEqualizationZNA.IsApplied() && fGainEqualizationZNC.IsApplied()) {
    // gain equalization
    qvectors_zdc.emplace("equalized", q_zn_equalized);
    AliZDCQvectors q_zn_equalized_recentered;
    // gain equalization + recentering
    q_zn_equalized_recentered.a = fRecenter4DAfterGainEqZNA.Apply(q_zn_equalized.a, variables.data());
    q_zn_equalized_recentered.c = fRecenter4DAfterGainEqZNC.Apply(q_zn_equalized.c, variables.data());
    if (fRecenter4DAfterGainEqZNA.IsApplied() && fRecenter4DAfterGainEqZNC.IsApplied()) {
      qvectors_zdc.emplace("equalized_recentered", q_zn_equalized_recentered);
    }
  }
  // recentering
  AliZDCQvectors q_zn_recentered;
  q_zn_recentered.a = fRecenter4DZNA.Apply(q_zn_plain.a, variables.data());
  q_zn_recentered.c = fRecenter4DZNC.Apply(q_zn_plain.c, variables.data());
  if (fRecenter4DZNA.IsApplied() && fRecenter4DZNC.IsApplied()) {
    qvectors_zdc.emplace("recentered", q_zn_recentered);
  }

  // Build Q-vector of V0A and V0C
  AliZDCQvectors q_v0_plain;
  const double kPiDiv8 = 0.39269908;
  const std::array<double, 8> phiv0 = {1 * kPiDiv8, 3 * kPiDiv8,  5 * kPiDiv8,  7 * kPiDiv8,
                                       9 * kPiDiv8, 11 * kPiDiv8, 13 * kPiDiv8, 15 * kPiDiv8};
  auto vzero = event->GetVZEROData();
  double vzerothreshold = 1e-3;
  for (unsigned int ich = 0; ich < 32; ich++) {
    double wv0a = vzero->GetMultiplicityV0A(ich);
    double wv0c = vzero->GetMultiplicityV0C(ich);
    // QnTools
    if (fActivateQnTools) {
      fValues[kV0AChPhi + ich] = 2. * phiv0[ich % 8];
      fValues[kV0CChPhi + ich] = 2. * phiv0[ich % 8];
      fValues[kV0AChMult + ich] = wv0a;
      fValues[kV0CChMult + ich] = wv0c;
    }
    if (wv0a > vzerothreshold) q_v0_plain.a.Update(2. * phiv0[ich % 8], wv0a);
    if (wv0c > vzerothreshold) q_v0_plain.c.Update(2. * phiv0[ich % 8], wv0c);
  }
  // normalize
  if (q_v0_plain.a.sum < min_sum || q_v0_plain.c.sum < min_sum) return;
  q_v0_plain.Normalize();
  // Add V0 Q vectors for the different correction steps
  std::map<std::string, AliZDCQvectors> qvectors_vzero;
  // plain
  qvectors_vzero.emplace("plain", q_v0_plain);
  // recentering
  AliZDCQvectors q_vzero_recentered;
  q_vzero_recentered.a = fRecenter4DV0A.Apply(q_v0_plain.a, variables.data());
  q_vzero_recentered.c = fRecenter4DV0C.Apply(q_v0_plain.c, variables.data());
  qvectors_vzero.emplace("recentered", q_vzero_recentered);

  // QnTools
  if (fActivateQnTools) {
    fCorrectionManager->ProcessEvent();
    fCorrectionManager->FillChannelDetectors();
  }
  // Track loop
  // builds TPC q-vector in cumulant analysis
  int ntracks_tpconly = 0;
  int ntracks_hybrid = 0;
  AliZDCflowCuts cut;
  cut.CheckEventCuts(event);
  const unsigned int ntracks = event->GetNumberOfTracks();
  for (unsigned int i = 0; i < ntracks; ++i) {
    auto track = static_cast<AliAODTrack *>(event->GetTrack(i));
    if (!track) continue;
    for (auto &analysis : fCumulantFlowAnalyses) analysis.FillPerTrackCorrelations(track);
    // QnTools
    if (fActivateQnTools) {
      if (track->TestFilterBit(128)) ++ntracks_tpconly;
      if (track->TestFilterBit(768)) ++ntracks_hybrid;
      if (!track->TestFilterBit(768)) continue;
      fValues[kTPCnCls] = track->GetTPCNcls();
      fValues[kTPCchi2pCls] = track->GetTPCchi2perCluster();
      fValues[kPhi] = track->Phi();
      fValues[kPt] = track->Pt();
      fValues[kEta] = track->Eta();
      fValues[kCharge] = track->GetSign();
      fCorrectionManager->FillTrackingDetectors();
    }
  }
  if (fActivateQnTools) {
    fValues[kNTPCTracksHybrid] = ntracks_hybrid;
    fValues[kNTPCTracksTPConly] = ntracks_tpconly;
    fCorrectionManager->ProcessCorrections();
  }

  std::map<std::string, Qv<4, 4>> qtpc;
  std::map<std::string, std::vector<Qv<4, 4>>> qtpcpt;
  std::map<std::string, Qv<4, 4>> qtpcp;
  std::map<std::string, Qv<4, 4>> qtpcn;
  for (auto &analysis : fCumulantFlowAnalyses) {
    qtpcpt.emplace(analysis.Name(), analysis.GetQTPCpt());
    qtpc.emplace(analysis.Name(), analysis.GetQTPC());
    qtpcp.emplace(analysis.Name(), analysis.GetQTPCEtaPos());
    qtpcn.emplace(analysis.Name(), analysis.GetQTPCEtaNeg());
  }
  // ZDC alignment
  AliQvector q_tpc = {qtpc["default"](1, 1).real(), qtpc["default"](1, 1).imag(), qtpc["default"](0, 1).real()};
  AliZDCQvectors q_align_c;
  q_align_c.a = fAlignZNA.Apply(q_zn_recentered.a, q_tpc, variables.data());
  q_align_c.c = fAlignZNC.Apply(q_zn_recentered.c, q_tpc, variables.data());
  qvectors_zdc.emplace("aligned", q_align_c);
  // ZDC Event-Shape-Engineering
  fQZNAmagnitude.Fill(q_zn_recentered.a, cent_v0m);
  fQZNCmagnitude.Fill(q_zn_recentered.c, cent_v0m);
  auto qpercentile_zna = fQZNAmagnitude.GetPercentile(q_zn_recentered.a, cent_v0m);
  auto qpercentile_znc = fQZNCmagnitude.GetPercentile(q_zn_recentered.c, cent_v0m);
  for (auto &analysis : fCumulantFlowAnalyses) {
    analysis.FillESE(qpercentile_zna, qpercentile_znc);
  }
  // ZDC event correlations
  for (auto &analysis : fEllipticFlowAnalyses) {
    analysis.SetQvectorsZDC(qvectors_zdc);
    analysis.SetQvectorsV0(qvectors_vzero);
  }
  // Calculate correlations
  for (auto &analysis : fEllipticFlowAnalyses) {
    if (qtpc.find(analysis.Name()) != qtpc.end()) {
      analysis.CalculateCorrelations(qtpc[analysis.Name()], qtpcpt[analysis.Name()], qtpcp[analysis.Name()],
                                     qtpcn[analysis.Name()]);
    } else {
      analysis.CalculateCorrelations(qtpc["default"], qtpcpt["default"], qtpcp["default"], qtpcn["default"]);
    }
  }
  // Calculate cumulants
  for (auto &analysis : fCumulantFlowAnalyses) {
    analysis.CalculateCumulants();
  }
  // Fill Eventplane angle QA histograms.
  fPsiZA->Fill(q_zn_plain.a.Psi());
  fPsiZC->Fill(q_zn_plain.c.Psi());
  fPsiZAEQ->Fill(q_zn_equalized.a.Psi());
  fPsiZCEQ->Fill(q_zn_equalized.c.Psi());
}

void AliAnalysisTaskFlowSpectators::Terminate(Option_t *option) {
  if (fActivateQnTools) fCorrectionManager->Finalize();
}

void AliAnalysisTaskFlowSpectators::AddCumulantFlowAnalyses(const YAML::Node &node) {
  for (YAML::const_iterator it = node.begin(); it != node.end(); ++it) {
    auto config = *it;
    fCumulantFlowAnalyses.emplace_back(config["name"].as<std::string>());
    fCumulantFlowAnalyses.back().ReadYAMLnode(config);
  }
}

void AliAnalysisTaskFlowSpectators::AddEllipticFlowAnalyses(const YAML::Node &node) {
  for (YAML::const_iterator it = node.begin(); it != node.end(); ++it) {
    auto config = *it;
    fEllipticFlowAnalyses.emplace_back(config["name"].as<std::string>());
    fEllipticFlowAnalyses.back().ReadYAMLnode(config);
  }
}

void AliAnalysisTaskFlowSpectators::ApplyConfiguration() {
  ConfigureCorrectionBinning(fDelayedNbinsXY, fDelayedNbinsZ, fDelayedEqualize);
  ConfigureCumulantAnalysis(fDelayedEtaGap, fDelayedPtBins, fDelayedVtxZBins, fDelayedNphiBins, fDelayedNetaBins,
                            fDelayedEtaMin, fDelayedEtaMax);
  ConfigureEllipticAnalysis(fDelayedPtBins);
}

void AliAnalysisTaskFlowSpectators::ConfigureEllipticAnalysis(const std::vector<double> &pt_bins) {
  for (auto &analysis : fEllipticFlowAnalyses) {
    analysis.SetPtBins(pt_bins);
  }
}

void AliAnalysisTaskFlowSpectators::ConfigureCumulantAnalysis(double eta_gap, const std::vector<double> &pt_bins,
                                                              const std::vector<double> &vtxz_bins, int n_phi_bins,
                                                              int n_eta_bins, double eta_min, double eta_max) {
  for (auto &analysis : fCumulantFlowAnalyses) {
    analysis.Configure(eta_gap, pt_bins, vtxz_bins, n_phi_bins, n_eta_bins, eta_min, eta_max);
    analysis.SetESE();
  }
}

void AliAnalysisTaskFlowSpectators::ConfigureCorrectionBinning(int nbinsxy, int nbinsz, bool equalize) {
  // ND correction
  std::vector<TAxis *> axes;
  auto make_axis = [](std::vector<TAxis *> &axes, std::string name, int n, double lo, double hi) {
    auto axis = new TAxis(n, lo, hi);
    axis->SetName(name.c_str());
    axes.push_back(axis);
  };
  make_axis(axes, "centV0M", 100, 0., 100.);
  make_axis(axes, "VtxX", 3, 0., 100.);  // binning configured run by run.
  make_axis(axes, "VtxY", 3, 0., 100.);  // binning configured run by run.
  make_axis(axes, "VtxZ", 3, 0., 100.);  // binning configured run by run.
  std::vector<std::string> rbrvars{"VtxX", "VtxY", "VtxZ"};
  fRecenter4DZNA.Configure("ZNA_4D", axes, rbrvars, equalize, 5);
  fRecenter4DZNC.Configure("ZNC_4D", axes, rbrvars, equalize, 5);

  fRecenter4DV0A.Configure("V0A_4D", axes, rbrvars, equalize, 5);
  fRecenter4DV0C.Configure("V0C_4D", axes, rbrvars, equalize, 5);

  fRecenter4DAfterGainEqZNA.Configure("ZNA_AfterGainEq_4D", axes, rbrvars, equalize, 5);
  fRecenter4DAfterGainEqZNC.Configure("ZNC_AfterGainEq_4D", axes, rbrvars, equalize, 5);

  std::vector<TAxis *> axes_align;
  make_axis(axes_align, "centV0M", 10, 0., 100.);
  fAlignZNA.Configure("ZNA_Align_all", axes_align, {}, 5);
  fAlignZNC.Configure("ZNC_Align_all", axes_align, {}, 5);
}