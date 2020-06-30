#include <iostream>
#include <cstdlib>
#include <cmath>

#include "TProfile.h"
#include "TProfile2D.h"
#include "TH1.h"
#include "TList.h"
#include "TFile.h"

#include "AliAnalysisTaskSE.h"

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODTracklets.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliOADBContainer.h"

#include "AliAnalysisTaskFlowZDC.h"

ClassImp(AliAnalysisTaskFlowZDC);

AliAnalysisTaskFlowZDC::AliAnalysisTaskFlowZDC() : AliAnalysisTaskSE(),
  fFilterBit(768),
  fVtxZcut(10),
  fNclsCut(70),
  fChi2MinCut(0.1),
  fChi2MaxCut(4),
  fPtMin(0.2),
  fPtMax(30.),
  fEtaMax(0.8),
  fPositiveOnly(false),
  fNegativeOnly(false),
  fRun(0) {
  }

AliAnalysisTaskFlowZDC::AliAnalysisTaskFlowZDC(const char *name) : AliAnalysisTaskSE(name),
  fFilterBit(768),
  fVtxZcut(10),
  fNclsCut(70),
  fChi2MinCut(0.1),
  fChi2MaxCut(4),
  fPtMin(0.2),
  fPtMax(30.),
  fEtaMax(0.8),
  fPositiveOnly(false),
  fNegativeOnly(false),
  fRun(0) {
  DefineOutput(1,TList::Class());
  DefineOutput(2,TList::Class());
  DefineOutput(3,TList::Class());
}

AliAnalysisTaskFlowZDC::~AliAnalysisTaskFlowZDC() {
  delete fOutputList;
  delete fCorrectionList;
  delete fQAList;
}

void AliAnalysisTaskFlowZDC::UserCreateOutputObjects() {
  fOutputList = new TList();
  fOutputList->SetOwner(true);
  fCorrectionList = new TList();
  fCorrectionList->SetOwner(true);
  fQAList = new TList();
  fQAList->SetOwner(true);
  const Double_t centBins[fNcentBins+1] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.,90.};
  fCentralityAxis = new TAxis(fNcentBins, centBins);
  const Int_t nPtBins = 24;
  const Double_t ptBins[nPtBins+1] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 3., 3.5, 4., 5., 6., 8., 10., 15., 20., 30.};
  const Int_t nPtBinsWide = 10;
  const Double_t ptBinsWide[nPtBinsWide+1] = {0.15, 0.35, 0.54, 0.94, 1.33, 1.73, 2.12, 2.51, 2.91, 3.7, 4.88};
  const Int_t nVtxXbins = 10;
  const Double_t vtxXbins[nVtxXbins+1] = {-0.03, -0.0128814, -0.0104336, -0.00860641, -0.00706022, -0.00561992, -0.00412592, -0.00256132, -0.00073528, 0.0017518, 0.02};
  const Int_t nVtxYbins = 10;
  const Double_t vtxYbins[nVtxYbins+1] = {0.16, 0.163816, 0.165859, 0.167457, 0.168857, 0.170201, 0.17157, 0.173047, 0.174767, 0.177141, 0.21};
  const Int_t nVtxZbins = 5;
  const Double_t vtxZbins[nVtxZbins+1] = {-10., -4.38032, -0.950887, 2.06604, 5.34186, 10.};
  const Int_t nEtaBins = 5;
  const Double_t etaBins[nEtaBins+1] = {-0.8,-0.48,-0.16,0.16,0.48,0.8};
  for (unsigned int i = 0; i < fNcentBins; ++i) {
    auto ibin = std::to_string(i);
    fV2XXXpT[i] = new TProfile((std::string("fVnXXXpt_")+ibin).data(),";#it{p}_{T} / GeV/#it{c};#LTx_{2}X_{1}^{ZNA}X_{1}^{ZNC}#GT",nPtBins,ptBins); 
    fV2XYYpT[i] = new TProfile((std::string("fVnXYYpt_")+ibin).data(),";#it{p}_{T} / GeV/#it{c};#LTx_{2}Y_{1}^{ZNA}Y_{1}^{ZNC}#GT",nPtBins,ptBins);
    fV2YXYpT[i] = new TProfile((std::string("fVnYXYpt_")+ibin).data(),";#it{p}_{T} / GeV/#it{c};#LTy_{2}X_{1}^{ZNA}Y_{1}^{ZNC}#GT",nPtBins,ptBins);
    fV2YYXpT[i] = new TProfile((std::string("fVnYYXpt_")+ibin).data(),";#it{p}_{T} / GeV/#it{c};#LTy_{2}Y_{1}^{ZNA}X_{1}^{ZNC}#GT",nPtBins,ptBins);
    fV2YYYpT[i] = new TProfile((std::string("fVnYYYpt_")+ibin).data(),";#it{p}_{T} / GeV/#it{c};#LTy_{2}Y_{1}^{ZNA}Y_{1}^{ZNC}#GT",nPtBins,ptBins);
    fV2YXXpT[i] = new TProfile((std::string("fVnYXXpt_")+ibin).data(),";#it{p}_{T} / GeV/#it{c};#LTy_{2}X_{1}^{ZNA}X_{1}^{ZNC}#GT",nPtBins,ptBins);
    fV2XXYpT[i] = new TProfile((std::string("fVnXXYpt_")+ibin).data(),";#it{p}_{T} / GeV/#it{c};#LTx_{2}X_{1}^{ZNA}Y_{1}^{ZNC}#GT",nPtBins,ptBins);
    fV2XYXpT[i] = new TProfile((std::string("fVnXYXpt_")+ibin).data(),";#it{p}_{T} / GeV/#it{c};#LTx_{2}Y_{1}^{ZNA}X_{1}^{ZNC}#GT",nPtBins,ptBins);
    fOutputList->Add(fV2XXXpT[i]);
    fOutputList->Add(fV2XYYpT[i]);
    fOutputList->Add(fV2YXYpT[i]);
    fOutputList->Add(fV2YYXpT[i]);
    fOutputList->Add(fV2YYYpT[i]);
    fOutputList->Add(fV2YXXpT[i]);
    fOutputList->Add(fV2XXYpT[i]);
    fOutputList->Add(fV2XYXpT[i]);
  }
  for (unsigned int i = 0; i < fNcentBinsWide; ++i) {
    auto ibin = std::to_string(i);
    fV1XXaEta[i] = new TProfile((std::string("fVnXXaEta_")+ibin).data(),";#eta;#LTx_{1}X_{1}^{ZNA}#GT",nEtaBins,etaBins); 
    fV1XXcEta[i] = new TProfile((std::string("fVnXXcEta_")+ibin).data(),";#eta;#LTx_{1}X_{1}^{ZNC}#GT",nEtaBins,etaBins); 
    fV1YYaEta[i] = new TProfile((std::string("fVnYYaEta_")+ibin).data(),";#eta;#LTy_{1}Y_{1}^{ZNA}#GT",nEtaBins,etaBins); 
    fV1YYcEta[i] = new TProfile((std::string("fVnYYcEta_")+ibin).data(),";#eta;#LTy_{1}Y_{1}^{ZNC}#GT",nEtaBins,etaBins); 
    fV1XYaEta[i] = new TProfile((std::string("fVnXYaEta_")+ibin).data(),";#eta;#LTx_{1}Y_{1}^{ZNA}#GT",nEtaBins,etaBins); 
    fV1XYcEta[i] = new TProfile((std::string("fVnXYcEta_")+ibin).data(),";#eta;#LTx_{1}Y_{1}^{ZNC}#GT",nEtaBins,etaBins); 
    fV1YXaEta[i] = new TProfile((std::string("fVnYXaEta_")+ibin).data(),";#eta;#LTy_{1}X_{1}^{ZNA}#GT",nEtaBins,etaBins); 
    fV1YXcEta[i] = new TProfile((std::string("fVnYXcEta_")+ibin).data(),";#eta;#LTy_{1}X_{1}^{ZNC}#GT",nEtaBins,etaBins); 
    fOutputList->Add(fV1XXaEta[i]); 
    fOutputList->Add(fV1XXcEta[i]); 
    fOutputList->Add(fV1YYaEta[i]); 
    fOutputList->Add(fV1YYcEta[i]); 
    fOutputList->Add(fV1XYaEta[i]); 
    fOutputList->Add(fV1XYcEta[i]); 
    fOutputList->Add(fV1YXaEta[i]); 
    fOutputList->Add(fV1YXcEta[i]); 
    fV1XXaPtEtaP[i] = new TProfile((std::string("fV1XXaPtEtaP")+ibin).data(), ";#it{p}_{T} / GeV/#it{c};#LTx_{1}X_{1}^{ZNA}#GT, #eta > 0", nPtBinsWide, ptBinsWide); 
    fV1XXcPtEtaP[i] = new TProfile((std::string("fV1XXcPtEtaP")+ibin).data(), ";#it{p}_{T} / GeV/#it{c};#LTx_{1}X_{1}^{ZNC}#GT, #eta > 0", nPtBinsWide, ptBinsWide); 
    fV1YYaPtEtaP[i] = new TProfile((std::string("fV1YYaPtEtaP")+ibin).data(), ";#it{p}_{T} / GeV/#it{c};#LTy_{1}Y_{1}^{ZNA}#GT, #eta > 0", nPtBinsWide, ptBinsWide); 
    fV1YYcPtEtaP[i] = new TProfile((std::string("fV1YYcPtEtaP")+ibin).data(), ";#it{p}_{T} / GeV/#it{c};#LTy_{1}Y_{1}^{ZNC}#GT, #eta > 0", nPtBinsWide, ptBinsWide); 
    fOutputList->Add(fV1XXaPtEtaP[i]); 
    fOutputList->Add(fV1XXcPtEtaP[i]); 
    fOutputList->Add(fV1YYaPtEtaP[i]); 
    fOutputList->Add(fV1YYcPtEtaP[i]); 
    fV1XXaPtEtaN[i] = new TProfile((std::string("fV1XXaPtEtaN")+ibin).data(), ";#it{p}_{T} / GeV/#it{c};#LTx_{1}X_{1}^{ZNA}#GT, #eta < 0", nPtBinsWide, ptBinsWide); 
    fV1XXcPtEtaN[i] = new TProfile((std::string("fV1XXcPtEtaN")+ibin).data(), ";#it{p}_{T} / GeV/#it{c};#LTx_{1}X_{1}^{ZNC}#GT, #eta < 0", nPtBinsWide, ptBinsWide); 
    fV1YYaPtEtaN[i] = new TProfile((std::string("fV1YYaPtEtaN")+ibin).data(), ";#it{p}_{T} / GeV/#it{c};#LTy_{1}Y_{1}^{ZNA}#GT, #eta < 0", nPtBinsWide, ptBinsWide); 
    fV1YYcPtEtaN[i] = new TProfile((std::string("fV1YYcPtEtaN")+ibin).data(), ";#it{p}_{T} / GeV/#it{c};#LTy_{1}Y_{1}^{ZNC}#GT, #eta < 0", nPtBinsWide, ptBinsWide); 
    fOutputList->Add(fV1XXaPtEtaN[i]); 
    fOutputList->Add(fV1XXcPtEtaN[i]); 
    fOutputList->Add(fV1YYaPtEtaN[i]); 
    fOutputList->Add(fV1YYcPtEtaN[i]); 
  }
  fXaXcCent = new TProfile("XXCent",";centrality V0M;#LTX_{1}^{ZNA}X_{1}^{ZNC}#GT",100,0.,100.);
  fYaYcCent = new TProfile("YYCent",";centrality V0M;#LTY_{1}^{ZNA}Y_{1}^{ZNC}#GT",100,0.,100.);
  fXaYcCent = new TProfile("XYCent",";centrality V0M;#LTX_{1}^{ZNA}Y_{1}^{ZNC}#GT",100,0.,100.);
  fYaXcCent = new TProfile("YXCent",";centrality V0M;#LTY_{1}^{ZNA}X_{1}^{ZNC}#GT",100,0.,100.);
  fOutputList->Add(fXaXcCent);
  fOutputList->Add(fYaYcCent);
  fOutputList->Add(fXaYcCent);
  fOutputList->Add(fYaXcCent);

  fXtXaCentEtaP = new TProfile("xxaCentP", ";centrality V0M;#LTx_{1}X_{1}^{ZNA}#GT, #eta > 0", 100, 0., 100.);
  fXtXcCentEtaP = new TProfile("xxcCentP", ";centrality V0M;#LTx_{1}X_{1}^{ZNC}#GT, #eta > 0", 100, 0., 100.);
  fYtYaCentEtaP = new TProfile("yyaCentP", ";centrality V0M;#LTy_{1}Y_{1}^{ZNA}#GT, #eta > 0", 100, 0., 100.);
  fYtYcCentEtaP = new TProfile("yycCentP", ";centrality V0M;#LTy_{1}Y_{1}^{ZNC}#GT, #eta > 0", 100, 0., 100.);
  fOutputList->Add(fXtXaCentEtaP);
  fOutputList->Add(fXtXcCentEtaP);
  fOutputList->Add(fYtYaCentEtaP);
  fOutputList->Add(fYtYcCentEtaP);
  fXtXaCentEtaN = new TProfile("xxaCentN", ";centrality V0M;#LTx_{1}X_{1}^{ZNA}#GT, #eta < 0", 100, 0., 100.);
  fXtXcCentEtaN = new TProfile("xxcCentN", ";centrality V0M;#LTx_{1}X_{1}^{ZNC}#GT, #eta < 0", 100, 0., 100.);
  fYtYaCentEtaN = new TProfile("yyaCentN", ";centrality V0M;#LTy_{1}Y_{1}^{ZNA}#GT, #eta < 0", 100, 0., 100.);
  fYtYcCentEtaN = new TProfile("yycCentN", ";centrality V0M;#LTy_{1}Y_{1}^{ZNC}#GT, #eta < 0", 100, 0., 100.);
  fOutputList->Add(fXtXaCentEtaN);
  fOutputList->Add(fXtXcCentEtaN);
  fOutputList->Add(fYtYaCentEtaN);
  fOutputList->Add(fYtYcCentEtaN);

  fXtXaXcCent = new TProfile("XXXCent","; centrality V0M;#LTx_{2}X_{1}^{ZNA}X_{1}^{ZNC}#GT",100,0.,100.);
  fXtYaYcCent = new TProfile("XYYCent","; centrality V0M;#LTx_{2}Y_{1}^{ZNA}Y_{1}^{ZNC}#GT",100,0.,100.);
  fYtXaYcCent = new TProfile("YXYCent","; centrality V0M;#LTy_{2}X_{1}^{ZNA}Y_{1}^{ZNC}#GT",100,0.,100.);
  fYtYaXcCent = new TProfile("YYXCent","; centrality V0M;#LTy_{2}Y_{1}^{ZNA}X_{1}^{ZNC}#GT",100,0.,100.);
  fYtXaXcCent = new TProfile("XYYCent","; centrality V0M;#LTy_{2}Y_{1}^{ZNA}X_{1}^{ZNC}#GT",100,0.,100.);
  fYtYaYcCent = new TProfile("YYYCent","; centrality V0M;#LTy_{2}Y_{1}^{ZNA}Y_{1}^{ZNC}#GT",100,0.,100.);
  fXtXaYcCent = new TProfile("XXYCent","; centrality V0M;#LTx_{2}X_{1}^{ZNA}Y_{1}^{ZNC}#GT",100,0.,100.);
  fXtYaXcCent = new TProfile("XYXCent","; centrality V0M;#LTx_{2}Y_{1}^{ZNA}X_{1}^{ZNC}#GT",100,0.,100.);
  fOutputList->Add(fXtXaXcCent);
  fOutputList->Add(fXtYaYcCent);
  fOutputList->Add(fYtXaYcCent);
  fOutputList->Add(fYtYaXcCent);
  fOutputList->Add(fYtXaXcCent);
  fOutputList->Add(fYtYaYcCent);
  fOutputList->Add(fXtXaYcCent);
  fOutputList->Add(fXtYaXcCent);
  
  // Corrections
  fXaCent = new TProfile("QxZAcent",";centrality V0M;X_{1}^{ZNA}",100,0.,100.,"s");
  fYaCent = new TProfile("QyZAcent",";centrality V0M;Y_{1}^{ZNA}",100,0.,100.,"s"); 
  fXcCent = new TProfile("QxZCcent",";centrality V0M;X_{1}^{ZNC}",100,0.,100.,"s"); 
  fYcCent = new TProfile("QyZCcent",";centrality V0M;Y_{1}^{ZNC}",100,0.,100.,"s"); 
  fCorrectionList->Add(fXaCent);
  fCorrectionList->Add(fYaCent);
  fCorrectionList->Add(fXcCent);
  fCorrectionList->Add(fYcCent);
  
  fXavXvY = new TProfile2D("QxZAvXvY",";vertex x / cm;vertex y / cm;X_{1}^{ZNA}",nVtxXbins,vtxXbins,nVtxYbins,vtxYbins,"s"); 
  fYavXvY = new TProfile2D("QyZAvXvY",";vertex x / cm;vertex y / cm;Y_{1}^{ZNA}",nVtxXbins,vtxXbins,nVtxYbins,vtxYbins,"s"); 
  fXcvXvY = new TProfile2D("QxZCvXvY",";vertex x / cm;vertex y / cm;X_{1}^{ZNC}",nVtxXbins,vtxXbins,nVtxYbins,vtxYbins,"s"); 
  fYcvXvY = new TProfile2D("QyZCvXvY",";vertex x / cm;vertex y / cm;Y_{1}^{ZNC}",nVtxXbins,vtxXbins,nVtxYbins,vtxYbins,"s"); 
  fCorrectionList->Add(fXavXvY); 
  fCorrectionList->Add(fYavXvY); 
  fCorrectionList->Add(fXcvXvY); 
  fCorrectionList->Add(fYcvXvY); 

  fXaCentVz = new TProfile2D("QxZAcentVz",";centrality V0M;vertex z / cm;X_{1}^{ZNA}",100,0.,100.,nVtxZbins,vtxZbins,"s"); 
  fYaCentVz = new TProfile2D("QyZAcentVz",";centrality V0M;vertex z / cm;Y_{1}^{ZNA}",100,0.,100.,nVtxZbins,vtxZbins,"s"); 
  fXcCentVz = new TProfile2D("QxZCcentVz",";centrality V0M;vertex z / cm;X_{1}^{ZNC}",100,0.,100.,nVtxZbins,vtxZbins,"s"); 
  fYcCentVz = new TProfile2D("QyZCcentVz",";centrality V0M;vertex z / cm;Y_{1}^{ZNC}",100,0.,100.,nVtxZbins,vtxZbins,"s"); 
  fCorrectionList->Add(fXaCentVz); 
  fCorrectionList->Add(fYaCentVz); 
  fCorrectionList->Add(fXcCentVz); 
  fCorrectionList->Add(fYcCentVz); 

  fPt = new TH1D("Pt",";p_{T} / GeV/c; N",nPtBins,ptBins);
  fQAList->Add(fPt);
  fPsiZA   = new TH1D("PsiZA",";#Psi_{ZNA};N",   50,0.,2*TMath::Pi());
  fPsiZC   = new TH1D("PsiZC",";#Psi_{ZNC};N",   50,0.,2*TMath::Pi());
  fPsiTPC1 = new TH1D("PsiTPC1",";#Psi_{TPC1};N",50,0.,2*TMath::Pi());
  fPsiTPC2 = new TH1D("PsiTPC2",";#Psi_{TPC2};N",50,0.,2*TMath::Pi());
  fQAList->Add(fPsiZA);
  fQAList->Add(fPsiZC);
  fQAList->Add(fPsiTPC1);
  fQAList->Add(fPsiTPC2);
  fCentralityV0M      = new TH1D("CentV0M",";centrality V0M;N",100,0.,100.);
  fCentralityCL1      = new TH1D("CentCL1",";centrality CL1;N",100,0.,100.);
  fCentralityCL1vsV0M = new TH2D("CentCL1vsV0M",";centrality CL1;centrality V0M;N",100,0.,100.,100,0.,100.);
  fQAList->Add(fCentralityV0M);
  fQAList->Add(fCentralityCL1);
  fQAList->Add(fCentralityCL1vsV0M);
  fVertexX = new TH1D("VertexX",";vertex x / cm;N", 50, -0.03, 0.02);
  fVertexY = new TH1D("VertexY",";vertex y / cm;N", 50, 0.16, 0.21);
  fVertexZ = new TH1D("VertexZ",";vertex z / cm;N", 50, -10., 10.);
  fQAList->Add(fVertexX);
  fQAList->Add(fVertexY);
  fQAList->Add(fVertexZ);
  fCorrectionStep = new TH1D("CorrectionStep",";correction step;", 4, 0., 4);
  fCorrectionList->Add(fCorrectionStep);

  fEventCuts.AddQAplotsToList(fQAList);

  PostData(1,fOutputList);
  PostData(2,fCorrectionList);
  PostData(3,fQAList);
}

void AliAnalysisTaskFlowZDC::NotifyRun() {
  OpenCorrections();
}

void AliAnalysisTaskFlowZDC::UserExec(Option_t *) {
  auto event = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!event) return;
  auto eventhandler = dynamic_cast<AliInputEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  auto trigger = eventhandler->IsEventSelected();
  if (!trigger) return;
  if (!fEventCuts.AcceptEvent(event)) return;
  const AliAODVertex* vtx = event->GetPrimaryVertex();
  if (std::abs(vtx->GetZ()) < fVtxZcut) {
    RunAnalysis(event);
  }
  PostData(1, fOutputList);
  PostData(2, fCorrectionList);
  PostData(3, fQAList);
}

void AliAnalysisTaskFlowZDC::RunAnalysis(AliAODEvent *event) {
  float cent_v0m = -1.;
  float cent_cl1 = -1.;
  auto centrality = event->GetCentrality();
  if (centrality) {
    cent_v0m = centrality->GetCentralityPercentile("V0M");
    cent_cl1 = centrality->GetCentralityPercentile("CL1");
  }
  if (cent_v0m < 0. || cent_v0m > 90.) return;
  std::vector<unsigned long> cent_bins_wide;
  if (cent_v0m > 10 && cent_v0m < 20) cent_bins_wide.push_back(0);
  if (cent_v0m > 30 && cent_v0m < 40) cent_bins_wide.push_back(1);
  if (cent_v0m > 10 && cent_v0m < 60) cent_bins_wide.push_back(2);

  auto cent_bin = fCentralityAxis->FindBin(cent_v0m) - 1;
  if (cent_bin < 0 && cent_bin >= fNcentBins) return;

  fCentralityV0M->Fill(cent_v0m);
  fCentralityCL1->Fill(cent_cl1);
  fCentralityCL1vsV0M->Fill(cent_cl1,cent_v0m);

  const AliAODVertex* vtx = event->GetPrimaryVertex();
  const auto vtxx = vtx->GetX();
  const auto vtxy = vtx->GetY();
  const auto vtxz = vtx->GetZ();
  fVertexX->Fill(vtxx);
  fVertexY->Fill(vtxy);
  fVertexZ->Fill(vtxz);
  
  // build ZN Q-vector
  auto zdc = event->GetZDCData();
  const std::array<double, 5> phizna{0, 5.4977871, 3.9269908, 0.78539816, 2.3561945};
  const std::array<double, 5> phiznc{0, 3.9269908, 5.4977871, 2.3561945, 0.78539816}; 
  Qv qza_plain = {0., 0., 0.};
  Qv qzc_plain = {0., 0., 0.};
  auto tower_energy_a = zdc->GetZNATowerEnergy();
  auto tower_energy_c = zdc->GetZNCTowerEnergy();
  for (int ich = 1; ich < 5; ++ich) {
    if (tower_energy_a[ich] > 0.) qza_plain.Update(phizna[ich], tower_energy_a[ich]);
    if (tower_energy_c[ich] > 0.) qzc_plain.Update(phiznc[ich], tower_energy_c[ich]);
  }
  if (!(qza_plain.sum > 0. && qzc_plain.sum > 0.)) return;
  // normalize Q vector ZNA and ZNC
  qza_plain.Normalize();
  qzc_plain.Normalize();
  // recentering
  auto qza = qza_plain;
  auto qzc = qzc_plain;
  RecenterZN(qza, qzc, cent_v0m, vtxx, vtxy, vtxz);

  // Fill Correlation ZNA ZNC
  fXaXcCent->Fill(cent_v0m, qza.x*qzc.x);
  fYaYcCent->Fill(cent_v0m, qza.y*qzc.y);
  fXaYcCent->Fill(cent_v0m, qza.x*qzc.y);
  fYaXcCent->Fill(cent_v0m, qza.y*qzc.x);

  // track loop
  Qv qtpc1_plain = {0., 0., 0.};
  Qv qtpc2_plain = {0., 0., 0.};
  unsigned int ntracks = event->GetNumberOfTracks();
  for (unsigned int i = 0; i < ntracks; ++i) {
    auto track = static_cast<AliAODTrack*>(event->GetTrack(i));
    if (!track) continue;
    if (!track->TestFilterBit(fFilterBit)) continue;
    if (track->GetTPCNcls() < fNclsCut) continue;
    if (track->GetTPCchi2perCluster() < fChi2MinCut || track->GetTPCchi2perCluster() > fChi2MaxCut) continue;
    const auto phi  = track->Phi();
    const auto pt   = track->Pt();
    const auto eta  = track->Eta();
    const auto sign = track->GetSign();
    if (std::abs(eta) > fEtaMax) continue;
    if (pt < fPtMin || pt > fPtMax) continue;
    if (fPositiveOnly && sign < 0) continue;
    if (fNegativeOnly && sign > 0) continue;
    fPt->Fill(pt);

    qtpc1_plain.Update(   phi, 1.);
    qtpc2_plain.Update(2.*phi, 1.);
    
    double v1_xxa = std::cos(phi)*qza.x;
    double v1_yya = std::sin(phi)*qza.y;
    double v1_yxa = std::sin(phi)*qza.x;
    double v1_xya = std::cos(phi)*qza.y;

    double v1_xxc = std::cos(phi)*qzc.x;
    double v1_yyc = std::sin(phi)*qzc.y;
    double v1_yxc = std::sin(phi)*qzc.x;
    double v1_xyc = std::cos(phi)*qzc.y;

    double v2_xxx = std::cos(2.*phi)*qza.x*qzc.x;
    double v2_xyy = std::cos(2.*phi)*qza.y*qzc.y;
    double v2_yyx = std::sin(2.*phi)*qza.y*qzc.x;
    double v2_yxy = std::sin(2.*phi)*qza.x*qzc.y;
    double v2_yyy = std::sin(2.*phi)*qza.y*qzc.y;
    double v2_yxx = std::sin(2.*phi)*qza.x*qzc.x;
    double v2_xyx = std::cos(2.*phi)*qza.y*qzc.x;
    double v2_xxy = std::cos(2.*phi)*qza.x*qzc.y;
    
    fV2XXXpT[cent_bin]->Fill(pt, v2_xxx);
    fV2XYYpT[cent_bin]->Fill(pt, v2_xyy);
    fV2YXYpT[cent_bin]->Fill(pt, v2_yxy);
    fV2YYXpT[cent_bin]->Fill(pt, v2_yyx);
    fV2YYYpT[cent_bin]->Fill(pt, v2_yyy);
    fV2YXXpT[cent_bin]->Fill(pt, v2_yxx);
    fV2XXYpT[cent_bin]->Fill(pt, v2_xxy);
    fV2XYXpT[cent_bin]->Fill(pt, v2_xyx);

    fXtXaXcCent->Fill(cent_v0m, std::cos(2.*phi)*qza.x*qzc.x);
    fXtXaYcCent->Fill(cent_v0m, std::cos(2.*phi)*qza.x*qzc.y);
    fXtYaXcCent->Fill(cent_v0m, std::cos(2.*phi)*qza.y*qzc.x);
    fYtXaXcCent->Fill(cent_v0m, std::sin(2.*phi)*qza.x*qzc.x);
    fYtYaYcCent->Fill(cent_v0m, std::sin(2.*phi)*qza.y*qzc.y);
    fYtYaXcCent->Fill(cent_v0m, std::sin(2.*phi)*qza.y*qzc.x);
    fXtYaYcCent->Fill(cent_v0m, std::cos(2.*phi)*qza.y*qzc.y);
    fYtXaYcCent->Fill(cent_v0m, std::sin(2.*phi)*qza.x*qzc.y);

    for (auto ibin : cent_bins_wide) {
      fV1XXaEta[ibin]->Fill(eta, v1_xxa);
      fV1XXcEta[ibin]->Fill(eta, v1_xxc);
      fV1YYaEta[ibin]->Fill(eta, v1_yya);
      fV1YYcEta[ibin]->Fill(eta, v1_yyc);
      fV1XYaEta[ibin]->Fill(eta, v1_xya);
      fV1XYcEta[ibin]->Fill(eta, v1_xyc);
      fV1YXaEta[ibin]->Fill(eta, v1_yxa);
      fV1YXcEta[ibin]->Fill(eta, v1_yxc);
    }
    if (eta > 0.) {
      fXtXaCentEtaP->Fill(cent_v0m,v1_xxa);
      fXtXcCentEtaP->Fill(cent_v0m,v1_xxc);
      fYtYaCentEtaP->Fill(cent_v0m,v1_yya);
      fYtYcCentEtaP->Fill(cent_v0m,v1_yya);
      for (auto ibin : cent_bins_wide) {
        fV1XXaPtEtaP[ibin]->Fill(pt,v1_xxa); 
        fV1XXcPtEtaP[ibin]->Fill(pt,v1_xxc); 
        fV1YYaPtEtaP[ibin]->Fill(pt,v1_yya); 
        fV1YYcPtEtaP[ibin]->Fill(pt,v1_yya); 
      }
    } else {
      fXtXaCentEtaN->Fill(cent_v0m,v1_xxa);
      fXtXcCentEtaN->Fill(cent_v0m,v1_xxc);
      fYtYaCentEtaN->Fill(cent_v0m,v1_yya);
      fYtYcCentEtaN->Fill(cent_v0m,v1_yya);
      for (auto ibin : cent_bins_wide) {
        fV1XXaPtEtaN[ibin]->Fill(pt,v1_xxa); 
        fV1XXcPtEtaN[ibin]->Fill(pt,v1_xxc); 
        fV1YYaPtEtaN[ibin]->Fill(pt,v1_yya); 
        fV1YYcPtEtaN[ibin]->Fill(pt,v1_yya); 
      }
    }
  }

  // normalize TPC Q vector
  qtpc1_plain.Normalize();
  qtpc2_plain.Normalize();

  const auto psi_zna  = TMath::Pi() + std::atan2(qza.y, qza.x);
  const auto psi_znc  = TMath::Pi() + std::atan2(qzc.y, qzc.x);
  const auto psi_tpc1 = TMath::Pi() + std::atan2(qtpc1_plain.y, qtpc1_plain.x);
  const auto psi_tpc2 = TMath::Pi() + std::atan2(qtpc2_plain.y, qtpc2_plain.x);

  fPsiZA->Fill(psi_zna);
  fPsiZC->Fill(psi_znc);
  fPsiTPC1->Fill(psi_tpc1);
  fPsiTPC2->Fill(psi_tpc2);
}

void AliAnalysisTaskFlowZDC::Terminate(Option_t *option) { }

void AliAnalysisTaskFlowZDC::Recenter(Qv &q, TH1 *cx, TH1 *cy, const Double_t coord_x, const Double_t coord_y = 0.) {
  auto ibin = cx->FindBin(coord_x, coord_y);
  if (cx->GetBinError(ibin) > 0. && cy->GetBinError(ibin) > 0.) {
    q.x = (q.x - cx->GetBinContent(ibin)) / cx->GetBinError(ibin);
    q.y = (q.y - cy->GetBinContent(ibin)) / cy->GetBinError(ibin);
  } else {
    q.x = q.x - cx->GetBinContent(ibin);
    q.y = q.y - cy->GetBinContent(ibin);
  }
}


void AliAnalysisTaskFlowZDC::RecenterZN(Qv &qa, Qv &qc, const Double_t cent, const Double_t vtxx, const Double_t vtxy, const Double_t vtxz) {
  fXaCent->Fill(cent, qa.x);
  fYaCent->Fill(cent, qa.y);
  fXcCent->Fill(cent, qc.x);
  fYcCent->Fill(cent, qc.y);
  unsigned int correction_step = 0;
  // step 1: centrality
  if (fMeanXaCent && fMeanYaCent && fMeanXcCent && fMeanYcCent) {
    Recenter(qa, fMeanXaCent, fMeanYaCent, cent);
    Recenter(qc, fMeanXcCent, fMeanYcCent, cent);
    fXavXvY->Fill(vtxx, vtxy, qa.x);
    fYavXvY->Fill(vtxx, vtxy, qa.y);
    fXcvXvY->Fill(vtxx, vtxy, qc.x);
    fYcvXvY->Fill(vtxx, vtxy, qc.y);
    correction_step = 1;
  }
  // step 2: vertex x and vertex y
  if (correction_step==1 && fMeanXavXvY && fMeanYavXvY && fMeanXcvXvY && fMeanYcvXvY) {
    Recenter(qa, fMeanXavXvY, fMeanYavXvY, vtxx, vtxy);
    Recenter(qc, fMeanXcvXvY, fMeanYcvXvY, vtxx, vtxy);
    fXaCentVz->Fill(cent, vtxz, qa.x);
    fYaCentVz->Fill(cent, vtxz, qa.y);
    fXcCentVz->Fill(cent, vtxz, qc.x);
    fYcCentVz->Fill(cent, vtxz, qc.y);
    correction_step = 2;
  }
  // step 3: centrality and vertex z 
  if (correction_step==2 && fMeanXaCentVz && fMeanYaCentVz && fMeanXcCentVz && fMeanYcCentVz) {
    Recenter(qa, fMeanXaCentVz, fMeanYaCentVz, cent, vtxz);
    Recenter(qc, fMeanXcCentVz, fMeanYcCentVz, cent, vtxz);
    correction_step = 3;
  }
  fCorrectionStep->Fill(correction_step);
}

bool AliAnalysisTaskFlowZDC::CheckBins(TH1 *histo) {
  for (int i = 0; i < histo->GetNbinsX(); ++i) {
    if (std::isnan(histo->GetBinContent(i)) || std::isnan(histo->GetBinError(i))) return false;
  }
  return true;
}

void AliAnalysisTaskFlowZDC::ReadFromOADB(TFile *file, const std::string &oadb_name, TProfile** profile) {
  auto oadb = dynamic_cast<AliOADBContainer*>(file->Get(oadb_name.data()));
  if (oadb) {
    auto tprofile = dynamic_cast<TProfile*>(oadb->GetObject(fCurrentRunNumber));
    if (tprofile && tprofile->GetEntries() > 0. && CheckBins(tprofile)) *profile = tprofile;
  }
}

void AliAnalysisTaskFlowZDC::ReadFromOADB(TFile *file, const std::string &oadb_name, TProfile2D** profile) {
  auto oadb = dynamic_cast<AliOADBContainer*>(file->Get(oadb_name.data()));
  if (oadb) {
    auto tprofile = dynamic_cast<TProfile2D*>(oadb->GetObject(fCurrentRunNumber));
    if (tprofile && tprofile->GetEntries() > 0. && CheckBins(tprofile)) *profile = tprofile;
  }
}

void AliAnalysisTaskFlowZDC::OpenCorrections() {
  auto file = TFile::Open(fInputFileName.data());
  if (!file || file->IsZombie()) return;
  // oadb containing first step
  ReadFromOADB(file,"MeanX_ZNA_Cent", &fMeanXaCent);
  ReadFromOADB(file,"MeanY_ZNA_Cent", &fMeanYaCent);
  ReadFromOADB(file,"MeanX_ZNC_Cent", &fMeanXcCent);
  ReadFromOADB(file,"MeanY_ZNC_Cent", &fMeanYcCent);
  // oadb containing first step
  ReadFromOADB(file,"MeanX_ZNA_VxVy", &fMeanXavXvY);
  ReadFromOADB(file,"MeanY_ZNA_VxVy", &fMeanYavXvY);
  ReadFromOADB(file,"MeanX_ZNC_VxVy", &fMeanXcvXvY);
  ReadFromOADB(file,"MeanY_ZNC_VxVy", &fMeanYcvXvY);
  // oadb containing first step
  ReadFromOADB(file,"MeanX_ZNA_CentVz", &fMeanXaCentVz);
  ReadFromOADB(file,"MeanY_ZNA_CentVz", &fMeanYaCentVz);
  ReadFromOADB(file,"MeanX_ZNC_CentVz", &fMeanXcCentVz);
  ReadFromOADB(file,"MeanY_ZNC_CentVz", &fMeanYcCentVz);
}
