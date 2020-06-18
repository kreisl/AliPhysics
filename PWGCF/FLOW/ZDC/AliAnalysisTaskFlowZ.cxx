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

#include <map>
#include "TList.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "AliAnalysisTaskSE.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODTracklets.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliOADBContainer.h"
#include "AliZDCflowCuts.h"
#include "AliAnalysisTaskFlowZ.h"
#include "AliAnalysisUtils.h"
#include "AliQvectors.h"

AliAnalysisTaskFlowZ::AliAnalysisTaskFlowZ() : 
  AliAnalysisTaskSE() {}

AliAnalysisTaskFlowZ::AliAnalysisTaskFlowZ(const char *name) : 
  AliAnalysisTaskSE(name) {
  DefineOutput(1,TList::Class());
  DefineOutput(2,TList::Class());
  DefineOutput(3,TList::Class());
  DefineOutput(4,TTree::Class());
}

AliAnalysisTaskFlowZ::~AliAnalysisTaskFlowZ() {
  delete fCorrelationList;
  delete fCorrectionList;
  delete fQAList;
  delete fTree;
  delete fAnalysisUtils;
}

void AliAnalysisTaskFlowZ::UserCreateOutputObjects() {
  fCorrelationList = new TList();
  fCorrelationList->SetOwner(true);
  fCorrectionList = new TList();
  fCorrectionList->SetOwner(true);
  fQAList = new TList();
  fQAList->SetOwner(true);

  // Apply delayed configuration
  ApplyConfiguration();

  fGainEqZNA.Configure("ZNA", 4, {4});
  fGainEqZNC.Configure("ZNC", 4, {4});

  fGainEqZNA.AddCorrectionsToList(fCorrectionList, fQAList);
  fGainEqZNC.AddCorrectionsToList(fCorrectionList, fQAList);

  auto event_cut_list = new TList();
  event_cut_list->SetOwner(true);
  event_cut_list->SetName("event_cuts");
  fEventCuts.AddQAplotsToList(event_cut_list);
  fQAList->Add(event_cut_list);

  for (auto &analysis : fDirectedFlowAnalyses) analysis.SetBootStrapSamples(fNsamples);
  for (auto &analysis : fEllipticFlowAnalyses) analysis.SetBootStrapSamples(fNsamples);
  for (auto &analysis : fCumulantFlowAnalyses) analysis.SetBootStrapSamples(fNsamples);

  auto directed_list = new TList();
  directed_list->SetOwner(true);
  directed_list->SetName("directed_flow");
  auto qa_directed_list = new TList();
  qa_directed_list->SetOwner(true);
  qa_directed_list->SetName("directed_flow");
  for (auto &analysis : fDirectedFlowAnalyses) {
    directed_list->Add(analysis.CreateCorrelations());
    analysis.AddQAHistograms(qa_directed_list);
  }
  fQAList->Add(qa_directed_list);
  fCorrelationList->Add(directed_list);

  auto elliptic_list = new TList();
  elliptic_list->SetOwner(true);
  elliptic_list->SetName("elliptic_flow");
  auto qa_elliptic_list = new TList();
  qa_elliptic_list->SetOwner(true);
  qa_elliptic_list->SetName("elliptic_flow");
  for (auto &analysis : fEllipticFlowAnalyses) {
    elliptic_list->Add(analysis.CreateCorrelations());
    analysis.AddQAHistograms(qa_elliptic_list);
  }
  fQAList->Add(qa_elliptic_list);
  fCorrelationList->Add(elliptic_list);

  auto cumulants_list = new TList();
  cumulants_list->SetOwner(true);
  cumulants_list->SetName("cumulants_flow");
  for (auto &analysis : fCumulantFlowAnalyses) {
    cumulants_list->Add(analysis.CreateCorrelations());
  }
  fCorrelationList->Add(cumulants_list);

  auto cumulants_qa_list = new TList();
  cumulants_qa_list->SetOwner(true);
  cumulants_qa_list->SetName("cumulants_flow");
  for (auto &analysis : fCumulantFlowAnalyses) {
    analysis.AddCorrectionsToList(fCorrectionList, cumulants_qa_list);
  }
  fQAList->Add(cumulants_qa_list);

  // event QA plots
  fCentralityV0M      = new TH1D("CentV0M",";centrality V0M;N",100,0.,100.);
  fCentralityCL1      = new TH1D("CentCL1",";centrality CL1;N",100,0.,100.);
  fCentralityCL1vsV0M = new TH2D("CentCL1vsV0M",";centrality CL1;centrality V0M;N",100,0.,100.,100,0.,100.);
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

  fMultiplicityTPConlyVsGlobal = new TH2D("MultiplicityTPConlyvsGlobal",";global;tpc only",1500, 0., 3000., 2000, 0., 4000.);
  fMultiplicityTPConlyVsGlobalCut = new TH2D("MultiplicityTPConlyvsGlobalCut",";global;tpc only",1500, 0., 3000., 2000, 0., 4000.);
  fNSigmaTPConlyVsGlobal = new TH2D("fNSigmaTPConlyVsGlobal",";global;N sigma ;tpc only",1500, 0., 3000., 500, -10., 10.);
  fNtracksESDvsNclsITS = new TH2D("NtracksESDvsNclsITS",";esd tracks;ITS clusters",1500, 0., 14000., 2000, 0., 26000.);
  fQAList->Add(fMultiplicityTPConlyVsGlobal);
  fQAList->Add(fMultiplicityTPConlyVsGlobalCut);
  fQAList->Add(fNSigmaTPConlyVsGlobal);
  fQAList->Add(fNtracksESDvsNclsITS);

  fPsiZA   = new TH1D("Psi_ZNA", ";#Psi_{ZNA};N", 50, 0., 2*TMath::Pi());
  fPsiZC   = new TH1D("Psi_ZNC", ";#Psi_{ZNC};N", 50, 0., 2*TMath::Pi());
  fPsiZAEQ = new TH1D("Psi_ZNA_EQ", ";#Psi_{ZNA};N", 50, 0., 2*TMath::Pi());
  fPsiZCEQ = new TH1D("Psi_ZNC_EQ", ";#Psi_{ZNC};N", 50, 0., 2*TMath::Pi());
  fQAList->Add(fPsiZA);
  fQAList->Add(fPsiZC);
  fQAList->Add(fPsiZAEQ);
  fQAList->Add(fPsiZCEQ);
  std::string correctionstepname("CorrectionStep");
  fCorrectionStep = new TH1D(correctionstepname.data(),"current correction step;step;", 4, 0., 4.);
  fCorrectionStepEQ = new TH1D((correctionstepname+"_EQ").data(),"current correction step after eq;step;", 4, 0., 4.);
  fCorrectionList->Add(fCorrectionStep);
  fCorrectionList->Add(fCorrectionStepEQ);

  fQZNAmagnitude.Configure("QZNAmagnitude", 1000, 180.);
  fQZNCmagnitude.Configure("QZNCmagnitude", 1000, 180.);
  fQZNAmagnitude.AddToOutputList(fCorrectionList);
  fQZNCmagnitude.AddToOutputList(fCorrectionList);

  // QA tree
  fTree = new TTree("tree","QA Tree of Q vectors");
  fTree->Branch("centV0M", &fBcentV0M);
  fTree->Branch("centCL1", &fBcentCL1);
  fTree->Branch("vtxX", &fBvtxX);
  fTree->Branch("vtxY", &fBvtxY);
  fTree->Branch("vtxZ", &fBvtxZ);
  fTree->Branch("run", &fCurrentRunNumber);
  fTree->Branch("X_ZNA", &fBxZNA);
  fTree->Branch("Y_ZNA", &fByZNA);
  fTree->Branch("M_ZNA", &fBsZNA);
  fTree->Branch("X_ZNC", &fBxZNC);
  fTree->Branch("Y_ZNC", &fByZNC);
  fTree->Branch("M_ZNC", &fBsZNC);
  fTree->Branch("X_ZNA_all", &fBxZNAall);
  fTree->Branch("Y_ZNA_all", &fByZNAall);
  fTree->Branch("M_ZNA_all", &fBsZNAall);
  fTree->Branch("X_ZNC_all", &fBxZNCall);
  fTree->Branch("Y_ZNC_all", &fByZNCall);
  fTree->Branch("M_ZNC_all", &fBsZNCall);

  fTree->Branch("X_TPC768_NUA", &fBxTPC768_NUA);
  fTree->Branch("Y_TPC768_NUA", &fByTPC768_NUA);
  fTree->Branch("X2_TPC768_NUA", &fBx2TPC768_NUA);
  fTree->Branch("Y2_TPC768_NUA", &fBy2TPC768_NUA);
  fTree->Branch("M_TPC768_NUA", &fBsTPC768_NUA);

  fTree->Branch("X_TPC768", &fBxTPC768);
  fTree->Branch("Y_TPC768", &fByTPC768);
  fTree->Branch("M_TPC768", &fBsTPC768);
  fTree->Branch("X2_TPC768", &fBx2TPC768);
  fTree->Branch("Y2_TPC768", &fBy2TPC768);

  fBnClsITS.resize(6);
  fTree->Branch("NclsITS", &fBnClsITS);
  fTree->Branch("NtracksESD", &fBnTracksESD);

  fTree->Branch("X_TPC96", &fBxTPC96);
  fTree->Branch("Y_TPC96", &fByTPC96);
  fTree->Branch("M_TPC96", &fBsTPC96);

  fTree->Branch("MultFB128", &fBmultFB128);
  fTree->Branch("MultFB768", &fBmultFB768);
  
  fAnalysisUtils = new AliAnalysisUtils();
  fAnalysisUtils->SetUseMVPlpSelection(true);
  fAnalysisUtils->SetUseOutOfBunchPileUp(true);
  
  fEventCuts.SetupRun1PbPb();

  PostData(1,fCorrelationList);
  PostData(2,fCorrectionList);
  PostData(3,fQAList);
  PostData(4,fTree);
}

// Reads corrections from the OADB file for each new run.
// Logs the status of corrections in the fCorrectionStep histogram.
void AliAnalysisTaskFlowZ::NotifyRun() {
  double correction_step_za = 0.;
  double correction_step_zc = 0.;
  double correction_step_za_eq = 0.;
  double correction_step_zc_eq = 0.;

  fCorrectZNAcentInterStep1.AddCorrectionsToList(fCorrectionList, fQAList);
  fCorrectZNCcentInterStep1.AddCorrectionsToList(fCorrectionList, fQAList);

  fCorrectZNAEQcentStep1.AddCorrectionsToList(fCorrectionList, fQAList);
  fCorrectZNCEQcentStep1.AddCorrectionsToList(fCorrectionList, fQAList);
  fCorrectZNAEQvXYZStep2.AddCorrectionsToList(fCorrectionList, fQAList);
  fCorrectZNCEQvXYZStep2.AddCorrectionsToList(fCorrectionList, fQAList);
  fCorrectZNAEQcentStep3.AddCorrectionsToList(fCorrectionList, fQAList);
  fCorrectZNCEQcentStep3.AddCorrectionsToList(fCorrectionList, fQAList);

  fCorrectZNAcentStep1.AddCorrectionsToList(fCorrectionList, fQAList);
  fCorrectZNCcentStep1.AddCorrectionsToList(fCorrectionList, fQAList);
  fCorrectZNAvXYZStep2.AddCorrectionsToList(fCorrectionList, fQAList);
  fCorrectZNCvXYZStep2.AddCorrectionsToList(fCorrectionList, fQAList);
  fCorrectZNAcentStep3.AddCorrectionsToList(fCorrectionList, fQAList);
  fCorrectZNCcentStep3.AddCorrectionsToList(fCorrectionList, fQAList);

  // Open Corrections
  auto file = TFile::Open(fCorrectionFileName.data());
  if (file && !file->IsZombie()) {

    fQZNAmagnitude.ReadFile(file);
    fQZNCmagnitude.ReadFile(file);

    auto graph3sigmap = dynamic_cast<TGraph*>(
        file->Get("TPConlyVsGlobalMeanPlusSigma3"));
    auto graph3sigmam = dynamic_cast<TGraph*>(
        file->Get("TPConlyVsGlobalMeanMinusSigma3"));
    auto graphmean = dynamic_cast<TGraph*>(
        file->Get("TPConlyVsGlobalMean"));
    if (graph3sigmap && graph3sigmam && graphmean) {
      std::cout << "multiplicity cut found" << std::endl;
      fMult3SigmaPlus = graph3sigmap;
      fMult3SigmaMinus = graph3sigmam;
      fMultMean = graphmean;
    }

    fCentralityWeightInput = dynamic_cast<TH1D*>(file->Get("CentralityWeights"));
    if (fCentralityWeightInput) {
      fCentralityWeightInput->SetDirectory(0);
      std::cout << "using Centrality weights" << std::endl;
    }

    fGainEqZNA.OpenCorrection(file,fCurrentRunNumber);
    fGainEqZNC.OpenCorrection(file,fCurrentRunNumber);
    // nd
    fAlignZNA.Make(file, fCurrentRunNumber, fCorrectionList, fQAList);
    fAlignZNC.Make(file, fCurrentRunNumber, fCorrectionList, fQAList);

    std::cout << "test" <<std::endl;
    fCorrectV0Aall.Make(file, fCurrentRunNumber, fCorrectionList, fQAList);
    fCorrectV0Call.Make(file, fCurrentRunNumber, fCorrectionList, fQAList);

    fCorrectZNAall.Make(file, fCurrentRunNumber, fCorrectionList, fQAList);
    fCorrectZNCall.Make(file, fCurrentRunNumber, fCorrectionList, fQAList);
    // interpolate
    fCorrectZNAcentInterStep1.OpenCorrection(file, fCurrentRunNumber);
    fCorrectZNAcentInterStep1.IsApplied();
    fCorrectZNCcentInterStep1.OpenCorrection(file, fCurrentRunNumber);
    fCorrectZNCcentInterStep1.IsApplied();
    // plain
    fCorrectZNAcentStep1.OpenCorrection(file, fCurrentRunNumber);
    fCorrectZNAvXYZStep2.OpenCorrection(file, fCurrentRunNumber);
    fCorrectZNAcentStep3.OpenCorrection(file, fCurrentRunNumber);
    if (fCorrectZNAcentStep1.IsApplied()) ++correction_step_za;
    if (fCorrectZNAvXYZStep2.IsApplied()) ++correction_step_za;
    if (fCorrectZNAcentStep3.IsApplied()) ++correction_step_za;
    fCorrectZNCcentStep1.OpenCorrection(file, fCurrentRunNumber);
    fCorrectZNCvXYZStep2.OpenCorrection(file, fCurrentRunNumber);
    fCorrectZNCcentStep3.OpenCorrection(file, fCurrentRunNumber);
    if (fCorrectZNCcentStep1.IsApplied()) ++correction_step_zc;
    if (fCorrectZNCvXYZStep2.IsApplied()) ++correction_step_zc;
    if (fCorrectZNCcentStep3.IsApplied()) ++correction_step_zc;
    // EQ
    fCorrectZNAEQcentStep1.OpenCorrection(file, fCurrentRunNumber);
    fCorrectZNAEQvXYZStep2.OpenCorrection(file, fCurrentRunNumber);
    fCorrectZNAEQcentStep3.OpenCorrection(file, fCurrentRunNumber);
    if (fCorrectZNAEQcentStep1.IsApplied()) ++correction_step_za_eq;
    if (fCorrectZNAEQvXYZStep2.IsApplied()) ++correction_step_za_eq;
    if (fCorrectZNAEQcentStep3.IsApplied()) ++correction_step_za_eq;
    fCorrectZNCEQcentStep1.OpenCorrection(file, fCurrentRunNumber);
    fCorrectZNCEQvXYZStep2.OpenCorrection(file, fCurrentRunNumber);
    fCorrectZNCEQcentStep3.OpenCorrection(file, fCurrentRunNumber);
    if (fCorrectZNCEQcentStep1.IsApplied()) ++correction_step_zc_eq;
    if (fCorrectZNCEQvXYZStep2.IsApplied()) ++correction_step_zc_eq;
    if (fCorrectZNCEQcentStep3.IsApplied()) ++correction_step_zc_eq;
    // cumulants
    for (auto &analysis : fCumulantFlowAnalyses) {
      analysis.OpenCorrection(file, fCurrentRunNumber);
      analysis.OpenPtEfficiencies(file);
    }
  }
  file->Close();
  fCorrectionStep->Fill(correction_step_za);
  fCorrectionStep->Fill(correction_step_zc);
  fCorrectionStepEQ->Fill(correction_step_za_eq);
  fCorrectionStepEQ->Fill(correction_step_zc_eq);
  PostData(2, fCorrectionList);
  PostData(3, fQAList);
}

void AliAnalysisTaskFlowZ::UserExec(Option_t*) {
  ResetTreeValues();
  auto event = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!event) return;
  auto eventhandler = dynamic_cast<AliInputEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!eventhandler->IsEventSelected()) return;
  if (!fEventCuts.AcceptEvent(event))   return;
  if (fAnalysisUtils->IsPileUpEvent(event)) return;
  if (MultiplicityCut(event)) {
    AnalyzeEvent(event);
  }
  fTree->Fill();
  PostData(1, fCorrelationList);
  PostData(2, fCorrectionList);
  PostData(3, fQAList);
  PostData(4, fTree);
}

bool AliAnalysisTaskFlowZ::MultiplicityCut(AliAODEvent* event) {
  bool pass = true;
  double ntracks_tpc_only = 0.;
  double ntracks_global = 0.;
  AliQvector qtpc768 = {0., 0., 0.};
  AliQvector q2tpc768 = {0., 0., 0.};
  AliQvector qtpc96 = {0., 0., 0.};
  // Track QA Qvector in tree
  const unsigned int ntracks = event->GetNumberOfTracks();
  for (unsigned int i = 0; i < ntracks; ++i) {
    auto track = static_cast<AliAODTrack*>(event->GetTrack(i));
    if (!track) continue;
    if (track->TestFilterBit(128)) ++ntracks_tpc_only;
    if (track->TestFilterBit(768)) ++ntracks_global;
    const auto tpc_ncls = track->GetTPCncls();
    const auto tpc_chi2_cls = track->GetTPCchi2perCluster();
    if (track->Pt() > 0.2 &&
        track->Pt() < 30. &&
        std::fabs(track->Eta()) < 0.8 &&
        tpc_ncls > 70 &&
        tpc_chi2_cls > 0.1 &&
        tpc_chi2_cls < 4) {
      if (track->TestFilterBit(768)) {
          qtpc768.Update(2.*track->Phi(), 1.0);
          q2tpc768.Update(4.*track->Phi(), 1.0);
      }
      if (track->TestFilterBit(96)) {
        qtpc96.Update(2.*track->Phi(), 1.0);
      }
    }
  }
  fBxTPC768 = qtpc768.x;
  fByTPC768 = qtpc768.y;
  fBsTPC768 = qtpc768.sum;
  fBx2TPC768 = q2tpc768.x;
  fBy2TPC768 = q2tpc768.y;
  fBxTPC96 = qtpc96.x;
  fByTPC96 = qtpc96.y;
  fBsTPC96 = qtpc96.sum;
  fBmultFB128 = ntracks_tpc_only;
  fBmultFB768 = ntracks_global;
  fBnTracksESD = event->GetNumberOfESDTracks();
  for (int i = 0; i < 6; ++i) {
    fBnClsITS[i] =  event->GetNumberOfITSClusters(i);
  }
  const auto nclsits = fBnClsITS[2]+fBnClsITS[3]+fBnClsITS[4]+fBnClsITS[5];
  fNtracksESDvsNclsITS->Fill(fBnTracksESD, nclsits);
  fMultiplicityTPConlyVsGlobal->Fill(ntracks_global, ntracks_tpc_only);
  if (fMultMean && fMult3SigmaPlus && fMult3SigmaMinus) {
    auto max = fMult3SigmaPlus->Eval(ntracks_global);
    auto min = fMult3SigmaMinus->Eval(ntracks_global);
    if (ntracks_tpc_only < min || ntracks_tpc_only > max) pass = false;
    auto mean = fMultMean->Eval(ntracks_global);
    auto nsigma = 3.*(ntracks_tpc_only - mean)/(max-mean); 
    fNSigmaTPConlyVsGlobal->Fill(ntracks_global, nsigma);
  }
  if (pass) {
    fMultiplicityTPConlyVsGlobalCut->Fill(ntracks_global, ntracks_tpc_only);
  }
  return pass;
}

void AliAnalysisTaskFlowZ::AnalyzeEvent(AliAODEvent* event) {
  //
  float cent_v0m = -1.;
  float cent_cl1 = -1.;
  auto centrality = event->GetCentrality();
  if (centrality) {
    cent_v0m = centrality->GetCentralityPercentile("V0M");
    cent_cl1 = centrality->GetCentralityPercentile("CL1");
  }
  if (cent_v0m < 0. || cent_v0m > 90.) return;

  const AliAODVertex* vtx = event->GetPrimaryVertex();
  const auto vtxx = vtx->GetX();
  const auto vtxy = vtx->GetY();
  const auto vtxz = vtx->GetZ();
  if (std::fabs(vtxz) > 10.) return;

  fVertexX->Fill(vtxx);
  fVertexY->Fill(vtxy);
  fVertexZ->Fill(vtxz);
  fVertexXY->Fill(vtxx, vtxy);

  fBvtxX = vtxx;
  fBvtxY = vtxy;
  fBvtxZ = vtxz;

  fCentralityV0M->Fill(cent_v0m);
  fCentralityCL1->Fill(cent_cl1);
  fCentralityCL1vsV0M->Fill(cent_cl1,cent_v0m);

  fBcentV0M = cent_v0m;
  fBcentCL1 = cent_cl1;

  std::vector<double> vvtx = {vtxx, vtxy, vtxz};
  GetSamples();

  std::vector<double> variables = {cent_v0m, vtxx, vtxy, vtxz};

  double centrality_weight = 1.;
  if (fCentralityWeightInput) {
    centrality_weight = fCentralityWeightInput->GetBinContent(
                          fCentralityWeightInput->FindBin(cent_v0m));
  }

  for (auto &analysis : fDirectedFlowAnalyses) {
    analysis.FindCentralityBin(event, cent_v0m, fSamples);
    analysis.SetCentralityWeight(centrality_weight);
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
  AliQvector qza_eq = {0., 0., 0.};
  AliQvector qzc_eq = {0., 0., 0.};
  AliQvector qza_pl = {0., 0., 0.};
  AliQvector qzc_pl = {0., 0., 0.};

  qza_eq = fGainEqZNA.ApplyGainEq(tower_energy_a + 1, phizna.data() + 1);
  qzc_eq = fGainEqZNC.ApplyGainEq(tower_energy_c + 1, phiznc.data() + 1);
  const auto min_energy = 1e-6;
  for (int ich = 1; ich < 5; ++ich) {
    if (tower_energy_a[ich] > min_energy) qza_pl.Update(phizna[ich], tower_energy_a[ich]);
    if (tower_energy_c[ich] > min_energy) qzc_pl.Update(phiznc[ich], tower_energy_c[ich]);
  }
  const auto min_sum = 1e-3;
  if (qza_pl.sum < min_sum || qzc_pl.sum < min_sum) return;
  qza_eq.Normalize();
  qzc_eq.Normalize();
  qza_pl.Normalize();
  qzc_pl.Normalize();

  std::map<std::string, AliZDCQvectors> qvectors_zdc;
  // Correction in iterations with gain equalization
  AliZDCQvectors q_eq_c;
  if (fGainEqZNA.IsApplied() && fGainEqZNC.IsApplied()) {
    q_eq_c.a = fCorrectZNAEQcentStep1.Correct(qza_eq, cent_v0m);
    q_eq_c.c = fCorrectZNCEQcentStep1.Correct(qzc_eq, cent_v0m);
    if (fCorrectZNAEQcentStep1.IsApplied() && fCorrectZNCEQcentStep1.IsApplied()) {
      qvectors_zdc.emplace("1st iteration gain equalization", q_eq_c);
      q_eq_c.a = fCorrectZNAEQvXYZStep2.Correct(q_eq_c.a, vvtx);
      q_eq_c.c = fCorrectZNCEQvXYZStep2.Correct(q_eq_c.c, vvtx);
    }
    if (fCorrectZNAEQvXYZStep2.IsApplied() && fCorrectZNCEQvXYZStep2.IsApplied()) {
      qvectors_zdc.emplace("2nd iteration gain equalization", q_eq_c);
      q_eq_c.a = fCorrectZNAEQcentStep3.Correct(q_eq_c.a, cent_v0m);
      q_eq_c.c = fCorrectZNCEQcentStep3.Correct(q_eq_c.c, cent_v0m);
    }
    if (fCorrectZNAEQcentStep3.IsApplied() && fCorrectZNCEQcentStep3.IsApplied()) {
      qvectors_zdc.emplace("3rd iteration gain equalization", q_eq_c);
    }
  }
  // Correction with interpolation
  AliZDCQvectors q_it_c;
  q_it_c.a = fCorrectZNAcentInterStep1.Correct(qza_pl, cent_v0m);
  q_it_c.c = fCorrectZNCcentInterStep1.Correct(qzc_pl, cent_v0m);
  qvectors_zdc.emplace("interpolated", q_it_c);

  // Correction in iterations
  AliZDCQvectors q_pl_c;
  q_pl_c.a = fCorrectZNAcentStep1.Correct(qza_pl, cent_v0m);
  q_pl_c.c = fCorrectZNCcentStep1.Correct(qzc_pl, cent_v0m);
  if (fCorrectZNAcentStep1.IsApplied() && fCorrectZNCcentStep1.IsApplied()) {
    qvectors_zdc.emplace("1st iteration", q_pl_c);
    q_pl_c.a = fCorrectZNAvXYZStep2.Correct(q_pl_c.a, vvtx);
    q_pl_c.c = fCorrectZNCvXYZStep2.Correct(q_pl_c.c, vvtx);
  }
  if (fCorrectZNAvXYZStep2.IsApplied() && fCorrectZNCvXYZStep2.IsApplied()) {
    qvectors_zdc.emplace("2nd iteration", q_pl_c);
    q_pl_c.a = fCorrectZNAcentStep3.Correct(q_pl_c.a, cent_v0m);
    q_pl_c.c = fCorrectZNCcentStep3.Correct(q_pl_c.c, cent_v0m);
  }
  if (fCorrectZNAcentStep3.IsApplied() && fCorrectZNCcentStep3.IsApplied()) {
    qvectors_zdc.emplace("3rd iteration", q_pl_c);
  }

  // Simultaneous ND correction  
  AliZDCQvectors q_all_c;
  q_all_c.a = fCorrectZNAall.Apply(qza_pl, variables.data());
  q_all_c.c = fCorrectZNCall.Apply(qzc_pl, variables.data());
  qvectors_zdc.emplace("nd", q_all_c);

  //Fill V0
  AliQvector qva_pl = {0., 0., 0.};
  AliQvector qvc_pl = {0., 0., 0.};
  const double kPiDiv8 = 0.39269908;
  const std::array<double, 8> phiv0 = {1*kPiDiv8, 3*kPiDiv8,   5*kPiDiv8,  7*kPiDiv8,
                                       9*kPiDiv8, 11*kPiDiv8, 13*kPiDiv8, 15*kPiDiv8};
  auto vzero = event->GetVZEROData();
  double vzerothreshold = 1e-3;
  for(unsigned int ich=0; ich<32; ich++){
    double weight = vzero->GetMultiplicityV0A(ich);
    if (weight > vzerothreshold) qva_pl.Update(2.*phiv0[ich%8], weight);
  }
  for(unsigned int ich=0; ich<32; ich++){
    double weight = vzero->GetMultiplicityV0C(ich);
    if (weight > vzerothreshold) qvc_pl.Update(2.*phiv0[ich%8], weight);
  }
  if (qva_pl.sum < min_sum || qvc_pl.sum < min_sum) return;
  qva_pl.Normalize();
  qvc_pl.Normalize();
  std::map<std::string, AliZDCQvectors> qvectors_vzero;
  AliZDCQvectors q_vzero_plain;
  q_vzero_plain.a = qva_pl;
  q_vzero_plain.c = qvc_pl;
  qvectors_vzero.emplace("plain", q_vzero_plain);
  // Simultaneous ND correction  
  AliZDCQvectors q_vzero_all_c;
  q_vzero_all_c.a = fCorrectV0Aall.Apply(qva_pl, variables.data());
  q_vzero_all_c.c = fCorrectV0Call.Apply(qvc_pl, variables.data());
  qvectors_vzero.emplace("nd", q_vzero_all_c);

  // Track loop 
  AliZDCflowCuts cut;
  cut.CheckEventCuts(event);
  const unsigned int ntracks = event->GetNumberOfTracks();
  for (unsigned int i = 0; i < ntracks; ++i) {
    auto track = static_cast<AliAODTrack*>(event->GetTrack(i));
    if (!track) continue;
    for (auto &analysis : fDirectedFlowAnalyses) analysis.FillPerTrackCorrelations(track);
    for (auto &analysis : fCumulantFlowAnalyses) analysis.FillPerTrackCorrelations(track);
  }

  std::map<std::string, Qv<4,4>> qtpc;
  std::map<std::string, std::vector<Qv<4,4>>> qtpcpt;
  std::map<std::string, Qv<4,4>> qtpcp;
  std::map<std::string, Qv<4,4>> qtpcn;
  for (auto &analysis : fCumulantFlowAnalyses) {
    qtpcpt.emplace(analysis.Name(), analysis.GetQTPCpt());
    qtpc.emplace(analysis.Name(), analysis.GetQTPC());
    qtpcp.emplace(analysis.Name(), analysis.GetQTPCEtaPos());
    qtpcn.emplace(analysis.Name(), analysis.GetQTPCEtaNeg());
  }
  //ZDC alignment
  AliQvector q_tpc = {qtpc["default"](1,1).real(), 
                      qtpc["default"](1,1).imag(),
                      qtpc["default"](0,1).real()};
  AliZDCQvectors q_align_c;
  q_align_c.a = fAlignZNA.Apply(q_all_c.a, q_tpc, variables.data());
  q_align_c.c = fAlignZNC.Apply(q_all_c.c, q_tpc, variables.data());
  qvectors_zdc.emplace("aligned", q_align_c);

  //ZNA magnitude
  fQZNAmagnitude.Fill(q_all_c.a, cent_v0m);
  fQZNCmagnitude.Fill(q_all_c.c, cent_v0m);

  auto qpercentile_zna = fQZNAmagnitude.GetPercentile(q_all_c.a, cent_v0m);
  auto qpercentile_znc = fQZNCmagnitude.GetPercentile(q_all_c.c, cent_v0m);

  for (auto &analysis : fCumulantFlowAnalyses) {
    analysis.FillESE(qpercentile_zna, qpercentile_znc);
  }

  // ZDC event correlations
  for (auto &analysis : fDirectedFlowAnalyses) analysis.SetQvectorsZDC(qvectors_zdc);
  for (auto &analysis : fEllipticFlowAnalyses) analysis.SetQvectorsZDC(qvectors_zdc);
  for (auto &analysis : fEllipticFlowAnalyses) analysis.SetQvectorsV0(qvectors_vzero);
  for (auto &analysis : fDirectedFlowAnalyses) analysis.FillPerEventCorrelations();


  //set q vector of ellipticflowanalysis
  for (auto &analysis : fEllipticFlowAnalyses) {
    if (qtpc.find(analysis.Name()) != qtpc.end()) {
      analysis.CalculateCorrelations(qtpc[analysis.Name()],
                                     qtpcpt[analysis.Name()],
                                     qtpcp[analysis.Name()],
                                     qtpcn[analysis.Name()]);
    } else {
      analysis.CalculateCorrelations(qtpc["default"],
                                     qtpcpt["default"],
                                     qtpcp["default"],
                                     qtpcn["default"]);
    }
  }

  // Calculate cumulants
  for (auto &analysis : fCumulantFlowAnalyses) {
    analysis.CalculateCumulants();
  }

  // QA tree variables
  
  fBxTPC768_NUA  = qtpc["default"](2,1).real();
  fByTPC768_NUA  = qtpc["default"](2,1).imag();
  fBx2TPC768_NUA = qtpc["default"](4,1).real();
  fBy2TPC768_NUA = qtpc["default"](4,1).imag();
  fBsTPC768_NUA  = qtpc["default"](0,1).real();
  
  fBxZNA = qza_pl.x;
  fByZNA = qza_pl.y;
  fBsZNA = qza_pl.sum;

  fBxZNC = qzc_pl.x;
  fByZNC = qzc_pl.y;
  fBsZNC = qzc_pl.sum;

  fBxZNAall = q_all_c.a.x;
  fByZNAall = q_all_c.a.y;
  fBsZNAall = q_all_c.a.sum;

  fBxZNCall = q_all_c.c.x;
  fByZNCall = q_all_c.c.y;
  fBsZNCall = q_all_c.c.sum;

  fPsiZA->Fill(q_pl_c.a.Psi());
  fPsiZC->Fill(q_pl_c.c.Psi());
  fPsiZAEQ->Fill(q_eq_c.a.Psi());
  fPsiZCEQ->Fill(q_eq_c.c.Psi());
}

void AliAnalysisTaskFlowZ::AddDirectedFlowAnalyses(const YAML::Node &node) {
  for (YAML::const_iterator it = node.begin(); it != node.end(); ++it) {
    auto config = *it;
    fDirectedFlowAnalyses.emplace_back(config["name"].as<std::string>());
    fDirectedFlowAnalyses.back().ReadYAMLnode(config);
  }
}

void AliAnalysisTaskFlowZ::AddCumulantFlowAnalyses(const YAML::Node &node) {
  for (YAML::const_iterator it = node.begin(); it != node.end(); ++it) {
    auto config = *it;
    fCumulantFlowAnalyses.emplace_back(config["name"].as<std::string>());
    fCumulantFlowAnalyses.back().ReadYAMLnode(config);
  }
}

void AliAnalysisTaskFlowZ::AddEllipticFlowAnalyses(const YAML::Node &node) {
  for (YAML::const_iterator it = node.begin(); it != node.end(); ++it) {
    auto config = *it;
    fEllipticFlowAnalyses.emplace_back(config["name"].as<std::string>());
    fEllipticFlowAnalyses.back().ReadYAMLnode(config);
  }
}

void AliAnalysisTaskFlowZ::ApplyConfiguration() {
  ConfigureCorrectionBinning(fDelayedNbinsXY, fDelayedNbinsZ, fDelayedEqualize);
  ConfigureCumulantAnalysis(fDelayedEtaGap, 
                            fDelayedPtBins,
                            fDelayedVtxZBins,
                            fDelayedNphiBins,
                            fDelayedNetaBins,
                            fDelayedEtaMin,
                            fDelayedEtaMax);
  ConfigureEllipticAnalysis(fDelayedPtBins);
}

void AliAnalysisTaskFlowZ::ConfigureEllipticAnalysis(const std::vector<double> &pt_bins) {
  for (auto &analysis : fEllipticFlowAnalyses) {
    analysis.SetPtBins(pt_bins);
  }
}

void AliAnalysisTaskFlowZ::ConfigureCumulantAnalysis(double eta_gap,
                                                        const std::vector<double> &pt_bins,
                                                        const std::vector<double> &vtxz_bins,
                                                        int n_phi_bins,
                                                        int n_eta_bins, 
                                                        double eta_min,
                                                        double eta_max) {
  for (auto &analysis : fCumulantFlowAnalyses) {
    analysis.Configure(eta_gap, pt_bins, vtxz_bins, n_phi_bins, n_eta_bins, eta_min, eta_max);
    analysis.SetESE();
  }
}

void AliAnalysisTaskFlowZ::ConfigureCorrectionBinning(int nbinsxy, int nbinsz, bool equalize) {
  // ND correction
  std::vector<TAxis*> axes;
  auto make_axis = [](std::vector<TAxis*> &axes, std::string name, int n, double lo, double hi) {
    auto axis = new TAxis(n, lo, hi);
    axis->SetName(name.c_str());
    axes.push_back(axis);
  };
  make_axis(axes,"centV0M", 100, 0., 100.);
  make_axis(axes,"VtxX"   ,   3, 0., 100.); // binning configured run by run.
  make_axis(axes,"VtxY"   ,   3, 0., 100.); // binning configured run by run.
  make_axis(axes,"VtxZ"   ,   3, 0., 100.); // binning configured run by run.
  std::vector<std::string> rbrvars{"VtxX", "VtxY", "VtxZ"};
  fCorrectZNAall.Configure("ZNA_all", axes, rbrvars, equalize, 5);
  fCorrectZNCall.Configure("ZNC_all", axes, rbrvars, equalize, 5);

  fCorrectV0Aall.Configure("V0A_all", axes, rbrvars, equalize, 5);
  fCorrectV0Call.Configure("V0C_all", axes, rbrvars, equalize, 5);

  std::vector<TAxis*> axes_align;
  make_axis(axes_align,"centV0M",10,0.,100.);
  fAlignZNA.Configure("ZNA_Align_all", axes_align, {}, 5);
  fAlignZNC.Configure("ZNC_Align_all", axes_align, {}, 5);

  std::vector<double_t> centrality_bins;
  for (double i = 0; i < 101; ++i) { centrality_bins.push_back(i); }
  // interpolated correction
  fCorrectZNAcentInterStep1.Configure("ZNA_it_step_1", "cent", centrality_bins, equalize);
  fCorrectZNCcentInterStep1.Configure("ZNC_it_step_1", "cent", centrality_bins, equalize);

  // step by step correction and gain equalized q-vectors
  fCorrectZNAcentStep1.Configure("ZNA_step_1", "cent", centrality_bins, equalize);
  fCorrectZNCcentStep1.Configure("ZNC_step_1", "cent", centrality_bins, equalize);
  fCorrectZNAEQcentStep1.Configure("ZNA_eq_step_1", "cent", centrality_bins, equalize);
  fCorrectZNCEQcentStep1.Configure("ZNC_eq_step_1", "cent", centrality_bins, equalize);

  fCorrectZNAvXYZStep2.Configure("ZNA_step_2", {"vx", "vy", "vz"}, nbinsxy, nbinsz, equalize);
  fCorrectZNCvXYZStep2.Configure("ZNC_step_2", {"vx", "vy", "vz"}, nbinsxy, nbinsz, equalize);
  fCorrectZNAEQvXYZStep2.Configure("ZNA_eq_step_2", {"vx", "vy", "vz"}, nbinsxy, nbinsz, equalize);
  fCorrectZNCEQvXYZStep2.Configure("ZNC_eq_step_2", {"vx", "vy", "vz"}, nbinsxy, nbinsz, equalize);

  fCorrectZNAcentStep3.Configure("ZNA_step_3", "cent", centrality_bins, equalize);
  fCorrectZNCcentStep3.Configure("ZNC_step_3", "cent", centrality_bins, equalize);
  fCorrectZNAEQcentStep3.Configure("ZNA_eq_step_3", "cent", centrality_bins, equalize);
  fCorrectZNCEQcentStep3.Configure("ZNC_eq_step_3", "cent", centrality_bins, equalize);

}

void AliAnalysisTaskFlowZ::ResetTreeValues() {
  fBcentV0M = -1.;
  fBcentCL1 = -1.;
  fBvtxX = -999.;
  fBvtxY = -999.;
  fBvtxZ = -999.;
  fBxZNA = 0.;
  fByZNA = 0.;
  fBsZNA = 0.;
  fBxZNC = 0.;
  fByZNC = 0.;
  fBsZNC = 0.;
  fBxZNAall = 0.;
  fByZNAall = 0.;
  fBsZNAall = 0.;
  fBxZNCall = 0.;
  fByZNCall = 0.;
  fBsZNCall = 0.;
  fBxTPC768 = 0.;
  fByTPC768 = 0.;
  fBx2TPC768 = 0.;
  fBy2TPC768 = 0.;
  fBsTPC768 = 0.;
  fBxTPC96 = 0.;
  fByTPC96 = 0.;
  fBsTPC96 = 0.;
  fBmultFB128 = 0;
  fBmultFB768 = 0;
  fBnTracksESD = 0;

  fBxTPC768_NUA  = 0.;
  fByTPC768_NUA  = 0.;
  fBx2TPC768_NUA = 0.;
  fBy2TPC768_NUA = 0.;
  fBsTPC768_NUA  = 0.;

  for (auto &n : fBnClsITS) n = 0;
}
