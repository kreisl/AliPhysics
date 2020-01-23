#include <iostream>
#include <cstdlib>
#include <cmath>

#include "TProfile.h"
#include "TProfile2D.h"
#include "TH1.h"

#include "AliAnalysisTaskSE.h"

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODTracklets.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliOADBContainer.h"

#include "AliAnalysisTaskFlowZDC.h"

ClassImp(AliAnalysisTaskFlowZDC);

using std::cout;
using std::endl;

AliAnalysisTaskFlowZDC::AliAnalysisTaskFlowZDC() : AliAnalysisTaskSE(),
  fCorrectionInputFile(nullptr),
  fFilterBit(768),
  fVtxZcut(10),
  fNclsCut(70),
  fChi2MinCut(0.1),
  fChi2MaxCut(4),
  fPtMin(0.2),
  fPtMax(30.),
  fEtaMax(0.8)
  fPositiveOnly(false),
  fNegativeOnly(false),
  fRun(0),
  fOutputList(nullptr) {
  }

AliAnalysisTaskFlowZDC::AliAnalysisTaskFlowZDC(const char *name) : AliAnalysisTaskSE(name),
  fCorrectionInputFile(nullptr),
  fFilterBit(768),
  fVtxZcut(10),
  fNclsCut(70),
  fChi2MinCut(0.1),
  fChi2MaxCut(4),
  fPtMin(0.2),
  fPtMax(30.),
  fEtaMax(0.8)
  fPositiveOnly(false),
  fNegativeOnly(false),
  fRun(0),
  fOutputList(nullptr) {
  DefineOutput(1,TList::Class());
}

AliAnalysisTaskFlowZDC::~AliAnalysisTaskFlowZDC() {
  delete fCorrectionInputFile;
  delete fOutputList;
}

void AliAnalysisTaskFlowZDC::UserCreateOutputObjects() {
  fOutputList = new TList();
  fOutputList->SetOwner(true);
  const int fNcentBins = 10;
  float centBins[fNCentBins+1] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.,90.};
  fCentralityAxis = new TAxis(fNcentBins, centBins);
  const int nPtBins = 24;
  float ptBins[nPtBins+1] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 3., 3.5, 4., 5., 6., 8., 10., 15., 20., 30.};
  const int nPtBinsWide = 10;
  float ptBinsWide[nPtBinsWide+1] = {0.15, 0.35, 0.54, 0.94, 1.33, 1.73, 2.12, 2.51, 2.91, 3.7, 4.88};
  const int nVtxXbins = 10;
  float vtxXbins[nVtxXbins+1] = {-0.03, -0.0128814, -0.0104336, -0.00860641, -0.00706022, -0.00561992, -0.00412592, -0.00256132, -0.00073528, 0.0017518, 0.02}
  const int nVtxYbins = 10;
  float vtxYbins[nVtxYbins+1] = {0.16, 0.163816, 0.165859, 0.167457, 0.168857, 0.170201, 0.17157, 0.173047, 0.174767, 0.177141, 0.21} 
  const int nVtxZbins = 5;
  float vtxZbins[nVtxZbins+1] = {-10., -4.38032, -0.950887, 2.06604, 5.34186, 10.} 
  for (unsigned int i = 0; i < fNCentBins; ++i) {
    auto ibin = std::to_string(i);
    fV2XXXpT[i] = new TProfile((std::string("fVnXXXpt_")+ibin).data(),"; p_{t} / GeV/c; v_{2}",nPtBins,ptBins); 
    fOutputList->Add(fV2XXXpT[i]);
    fV2XYYpT[i] = new TProfile((std::string("fVnXYYpt_")+ibin).data(),"; p_{t} / GeV/c; v_{2}",nPtBins,ptBins);
    fOutputList->Add(fV2XYYpT[i]);
    fV2YXYpT[i] = new TProfile((std::string("fVnYXYpt_")+ibin).data(),"; p_{t} / GeV/c; v_{2}",nPtBins,ptBins);
    fOutputList->Add(fV2YXYpT[i]);
    fV2YYXpT[i] = new TProfile((std::string("fVnYYXpt_")+ibin).data(),"; p_{t} / GeV/c; v_{2}",nPtBins,ptBins);
    fOutputList->Add(fV2YYXpT[i]);
    fV2YYYpT[i] = new TProfile((std::string("fVnYYYpt_")+ibin).data(),"; p_{t} / GeV/c; v_{2}",nPtBins,ptBins);
    fOutputList->Add(fV2YYYpT[i]);
    fV2YXXpT[i] = new TProfile((std::string("fVnYXXpt_")+ibin).data(),"; p_{t} / GeV/c; v_{2}",nPtBins,ptBins);
    fOutputList->Add(fV2YXXpT[i]);
    fV2XXYpT[i] = new TProfile((std::string("fVnXXYpt_")+ibin).data(),"; p_{t} / GeV/c; v_{2}",nPtBins,ptBins);
    fOutputList->Add(fV2XXYpT[i]);
    fV2XYXpT[i] = new TProfile((std::string("fVnXYXpt_")+ibin).data(),"; p_{t} / GeV/c; v_{2}",nPtBins,ptBins);
    fOutputList->Add(fV2XYXpT[i]);
  }
  for (unsigned int i = 0; i < fNcentBinsWide; ++i) {
    auto ibin = std::to_string(i);
    fV1XXaEta[i] = new TProfile((std::string("fVnXXaEta_")+ibin).data(),";#eta; v_{1}",nEtaBins,etaBins); 
    fOutputList->Add(fV1XXaEta); 
    fV1XXcEta[i] = new TProfile((std::string("fVnXXcEta_")+ibin).data(),";#eta; v_{1}",nEtaBins,etaBins); 
    fOutputList->Add(fV1XXcEta); 
    fV1YYaEta[i] = new TProfile((std::string("fVnYYaEta_")+ibin).data(),";#eta; v_{1}",nEtaBins,etaBins); 
    fOutputList->Add(fV1YYaEta); 
    fV1YYcEta[i] = new TProfile((std::string("fVnYYcEta_")+ibin).data(),";#eta; v_{1}",nEtaBins,etaBins); 
    fOutputList->Add(fV1YYcEta); 
    fV1XYaEta[i] = new TProfile((std::string("fVnXYaEta_")+ibin).data(),";#eta; v_{1}",nEtaBins,etaBins); 
    fOutputList->Add(fV1XYaEta); 
    fV1XYcEta[i] = new TProfile((std::string("fVnXYcEta_")+ibin).data(),";#eta; v_{1}",nEtaBins,etaBins); 
    fOutputList->Add(fV1XYcEta); 
    fV1YXaEta[i] = new TProfile((std::string("fVnYXaEta_")+ibin).data(),";#eta; v_{1}",nEtaBins,etaBins); 
    fOutputList->Add(fV1YXaEta); 
    fV1YXcEta[i] = new TProfile((std::string("fVnYXcEta_")+ibin).data(),";#eta; v_{1}",nEtaBins,etaBins); 
    fOutputList->Add(fV1YXcEta); 
    fV1XXaPtEtaP[i] = new TProfile((std::string("fV1XXaPtEtaP")+ibin).data(), ";p_{T} / GeV/#it{c}; v_{1}", nPtBinsWide, ptBinsWide); 
    fOutputList->Add(fV1XXaPtEtaP[i]); 
    fV1XXcPtEtaP[i] = new TProfile((std::string("fV1XXcPtEtaP")+ibin).data(), ";p_{T} / GeV/#it{c}; v_{1}", nPtBinsWide, ptBinsWide); 
    fOutputList->Add(fV1XXcPtEtaP[i]); 
    fV1YYaPtEtaP[i] = new TProfile((std::string("fV1YYaPtEtaP")+ibin).data(), ";p_{T} / GeV/#it{c}; v_{1}", nPtBinsWide, ptBinsWide); 
    fOutputList->Add(fV1YYaPtEtaP[i]); 
    fV1YYcPtEtaP[i] = new TProfile((std::string("fV1YYcPtEtaP")+ibin).data(), ";p_{T} / GeV/#it{c}; v_{1}", nPtBinsWide, ptBinsWide); 
    fOutputList->Add(fV1YYcPtEtaP[i]); 
    fV1XXaPtEtaN[i] = new TProfile((std::string("fV1XXaPtEtaN")+ibin).data(), ";p_{T} / GeV/#it{c}; v_{1}", nPtBinsWide, ptBinsWide); 
    fOutputList->Add(fV1XXaPtEtaN[i]); 
    fV1XXcPtEtaN[i] = new TProfile((std::string("fV1XXcPtEtaN")+ibin).data(), ";p_{T} / GeV/#it{c}; v_{1}", nPtBinsWide, ptBinsWide); 
    fOutputList->Add(fV1XXcPtEtaN[i]); 
    fV1YYaPtEtaN[i] = new TProfile((std::string("fV1YYaPtEtaN")+ibin).data(), ";p_{T} / GeV/#it{c}; v_{1}", nPtBinsWide, ptBinsWide); 
    fOutputList->Add(fV1YYaPtEtaN[i]); 
    fV1YYcPtEtaN[i] = new TProfile((std::string("fV1YYcPtEtaN")+ibin).data(), ";p_{T} / GeV/#it{c}; v_{1}", nPtBinsWide, ptBinsWide); 
    fOutputList->Add(fV1YYcPtEtaN[i]); 
  }
  fXznaXzncCent = new TProfile("XXCent",";centrality V0M; xx",100,0.,100.);
  fOutputList->Add(XznaXzncCent);
  fYznaYzncCent = new TProfile("YYCent",";centrality V0M; yy",100,0.,100.);
  fOutputList->Add(YznaYzncCent);
  fXznaYzncCent = new TProfile("XYCent",";centrality V0M; xy",100,0.,100.);
  fOutputList->Add(XznaYzncCent);
  fYznaXzncCent = new TProfile("YXCent",";centrality V0M; yx",100,0.,100.);
  fOutputList->Add(YznaXzncCent);

  fXtpcXznaCentEtaP = new TProfile("xxaCentP", ";centrality V0M; xxa (#eta > 0)", 100, 0., 100.);
  fOutputList->Add(fXtpcXznaCentEtaP);
  fXtpcXzncCentEtaP = new TProfile("xxcCentP", ";centrality V0M; xxc (#eta > 0)", 100, 0., 100.);
  fOutputList->Add(fXtpcXzncCentEtaP);
  fYtpcYznaCentEtaP = new TProfile("yyaCentP", ";centrality V0M; yya (#eta > 0)", 100, 0., 100.);
  fOutputList->Add(fYtpcYznaCentEtaP);
  fYtpcYzncCentEtaP = new TProfile("yycCentP", ";centrality V0M; yyc (#eta > 0)", 100, 0., 100.);
  fOutputList->Add(fYtpcYzncCentEtaP);
  fXtpcXznaCentEtaN = new TProfile("xxaCentN", ";centrality V0M; xxa (#eta < 0)", 100, 0., 100.);
  fOutputList->Add(fXtpcXznaCentEtaN);
  fXtpcXzncCentEtaN = new TProfile("xxcCentN", ";centrality V0M; xxc (#eta < 0)", 100, 0., 100.);
  fOutputList->Add(fXtpcXzncCentEtaN);
  fYtpcYznaCentEtaN = new TProfile("yyaCentN", ";centrality V0M; yya (#eta < 0)", 100, 0., 100.);
  fOutputList->Add(fYtpcYznaCentEtaN);
  fYtpcYzncCentEtaN = new TProfile("yycCentN", ";centrality V0M; yyc (#eta < 0)", 100, 0., 100.);
  fOutputList->Add(fYtpcYzncCentEtaN);

  fXtpcXznaXzncCent = new TProfile("XXXCent","; centrality V0M, xxx",100,0.,100.);
  fOutputList->Add(fXtpcXznaXzncCent);
  fXtpcYznaYzncCent = new TProfile("XYYCent","; centrality V0M, xyy",100,0.,100.);
  fOutputList->Add(fXtpcYznaYzncCent);
  fYtpcXznaYzncCent = new TProfile("YXYCent","; centrality V0M, yxy",100,0.,100.);
  fOutputList->Add(fYtpcXznaYzncCent);
  fYtpcYznaXzncCent = new TProfile("YYXCent","; centrality V0M, yyx",100,0.,100.);
  fOutputList->Add(fYtpcYznaXzncCent);
  fYtpcYznaYzncCent = new TProfile("YYYCent","; centrality V0M, yyy",100,0.,100.);
  fOutputList->Add(fYtpcYznaYzncCent);
  fXtpcXznaYzncCent = new TProfile("XXYCent","; centrality V0M, xxy",100,0.,100.);
  fOutputList->Add(fXtpcXznaYzncCent);
  fXtpcYznaXzncCent = new TProfile("XYXCent","; centrality V0M, xyx",100,0.,100.);
  fOutputList->Add(fXtpcYznaXzncCent);

  fQXZACent = new TProfile("QxZAcent",";centrality V0M; X_{1}{ZNA}",100,0,100);
  fOutputList->Add(fQXZACent);
  fQYZACent = new TProfile("QyZAcent",";centrality V0M; Y_{1}{ZNA}",100,0,100); 
  fOutputList->Add(fQYZACent); 
  fQXZCCent = new TProfile("QxZCcent",";centrality V0M; X_{1}{ZNC}",100,0,100); 
  fOutputList->Add(fQXZCCent); 
  fQYZCCent = new TProfile("QyZCcent",";centrality V0M; Y_{1}{ZNC}",100,0,100); 
  fOutputList->Add(fQYZCCent); 
  
  fQXZAvXvY = new TProfile2D("QxZAvXvY",";vertex X;vertex Y; X_{1}{ZNA}",nVtxXbins,vtxXbins,nVtxYbins,vtxYbins); 
  fOutputList->Add(fQXZAvXvY); 
  fQYZAvXvY = new TProfile2D("QyZAvXvY",";vertex X;vertex Y; Y_{1}{ZNA}",nVtxXbins,vtxXbins,nVtxYbins,vtxYbins); 
  fOutputList->Add(fQYZAvXvY); 
  fQXZCvXvY = new TProfile2D("QxZCvXvY",";vertex X;vertex Y; X_{1}{ZNC}",nVtxXbins,vtxXbins,nVtxYbins,vtxYbins); 
  fOutputList->Add(fQXZCvXvY); 
  fQYZCvXvY = new TProfile2D("QyZCvXvY",";vertex X;vertex Y; Y_{1}{ZNC}",nVtxXbins,vtxXbins,nVtxYbins,vtxYbins); 
  fOutputList->Add(fQYZCvXvY); 

  fQXZACentVz = new TProfile2D("QxZAcentVz",";centrality V0M; vertex Z; X_{1}{ZNA}",100,0,100,nVtxZbins,vtxZbins); 
  fOutputList->Add(fQXZACentVz); 
  fQYZACentVz = new TProfile2D("QyZAcentVz",";centrality V0M; vertex Z; Y_{1}{ZNA}",100,0,100,nVtxZbins,vtxZbins); 
  fOutputList->Add(fQYZACentVz); 
  fQXZCCentVz = new TProfile2D("QxZCcentVz",";centrality V0M; vertex Z; X_{1}{ZNC}",100,0,100,nVtxZbins,vtxZbins); 
  fOutputList->Add(fQXZCCentVz); 
  fQYZCCentVz = new TProfile2D("QyZCcentVz",";centrality V0M; vertex Z; Y_{1}{ZNC}",100,0,100,nVtxZbins,vtxZbins); 
  fOutputList->Add(fQYZCCentVz); 

  fPt = new TH1D("Pt",";p_{T} / GeV/c; N",nPtBins,ptBins);
  fOutputList->Add(fPt);
  fPsiZA = new TH1D("PsiZA","; $Psi_{ZNA}; N",20,0,2*TMath::Pi());
  fOutputList->Add(fPsiZA);
  fPsiZC = new TH1D("PsiZC","; $Psi_{ZNC}; N",20,0,2*TMath::Pi());
  fOutputList->Add(fPsiZC);
  fPsiTPC1 = new TH1D("PsiTPC1","; $Psi_{TPC1}; N",20,0,2*TMath::Pi());
  fOutputList->Add(fPsiTPC1);
  fPsiTPC2 = new TH1D("PsiTPC2","; $Psi_{TPC2}; N",20,0,2*TMath::Pi());
  fOutputList->Add(fPsiTPC2);
  fCentralityV0M = new TH1D("CentralityV0M",";centrality V0M; N",100,0,100);
  fOutputList->Add(fCentralityV0M);
  fCentralityCL1 = new TH1D("CentralityCL1",";centrality CL1; N",100,0,100);
  fOutputList->Add(fCentralityCL1);
  fCentralityCL1vsV0M = new TH2D("CentralityCL1vsV0M",";centrality CL1vsV0M; N",100,0.,100.,100,0.,100.);
  fOutputList->Add(fCentralityCL1vsV0M);
  fVertexX = new TH1D("VertexX",";vertex x / cm; N",50,-0.03,0.02);
  fVertexY = new TH1D("VertexY",";vertex y / cm; N",50,0.16,0.21);
  fVertexZ = new TH1D("VertexZ",";vertex z / cm; N",50,-10,10);
  fCorrectionStep = new TH1D("CorrectionStep",";correction step;",4,0,4);

  PostData(1,fOutputList);
}

void AliAnalysisTaskFlowZDC::UserExec(Option_t *) {
  fEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fEvent) return;
  auto eventhandler = dynamic_cast<AliInputEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  auto trigger = eventhandler->IsEventSelected();
  if(!trigger) return;
  auto run = fEvent->GetRunNumber();
  if (run != fRun) {
    fRun = run;
    OpenCorrections(fInputFile, fRun);
  }
  const AliAODVertex* vtx = fEvent->GetPrimaryVertex();
  if (std::abs(vtxtrk->GetZ()) < fVtxZcut) {
    runAnalysis(fEvent);
  }
  PostData(1, fOutputList);
}

void RunAnalysis(AliAODEvent *event) {
  float cent_v0m = -1.;
  float cent_cl1 = -1.;
  auto centrality = fEvent->GetCentrality(); // run 1
  if (centrality) {
    cent_v0m = centrality->GetCentralityPercentile("V0M");
    cent_cl1 = centrality->GetCentralityPercentile("CL1");
  }
  if (cent_v0m < 0. || cent_v0m > 90.) return;
  auto cent_bin_wide = -1;
  if (cent_v0m > 10 && cent_v0m < 20) cent_bin_wide = 0;
  if (cent_v0m > 30 && cent_v0m < 40) cent_bin_wide = 0;
  if (cent_v0m > 10 && cent_v0m < 60) cent_bin_wide = 0;

  auto cent_bin = fCentralityAxis->FindBin(cent_v0m) - 1;
  if (cent_bin < 0 && cent_bin >= fNcentBins) return;

  fCentralityV0M->Fill(cent_v0m);
  fCentralityCL1->Fill(cent_cl1);

  const AliAODVertex* vtx = fEvent->GetPrimaryVertex();
  fVertexX->Fill(vtxtrk->GetX());
  fVertexY->Fill(vtxtrk->GetY());
  fVertexZ->Fill(vtxtrk->GetZ());
  
  // build ZN Q-vector
  auto zdc = fEvent->GetZDCData();
  const std::array<double, 5> phizna{0, 5.4977871, 3.9269908, 0.78539816, 2.3561945};
  const std::array<double, 5> phiznc{0, 3.9269908, 5.4977871, 2.3561945, 0.78539816}; 
  Qv qza_plain = {0., 0., 0.};
  Qv qzc_plain = {0., 0., 0.};
  auto tower_energy_a = zdc->GetZNATowerEnergy();
  auto tower_energy_c = zdc->GetZNCTowerEnergy();
  for (int ich = 0; ich < 5; ++ich) {
    if(tower_energy_a[ich] > 0.) {
      qza_plain.x += std::cos(phizna[ich]) * tower_energy_a[ich];
      qza_plain.y += std::sin(phizna[ich]) * tower_energy_a[ich];
      qza_plain.sum += tower_energy_a[ich];
    }
    if(tower_energy[ich] > 0.) {
      qzc_plain.x += std::cos(phiznc[ich]) * tower_energy_c[ich];
      qzc_plain.y += std::sin(phiznc[ich]) * tower_energy_c[ich];
      qzc_plain.sum += tower_energy_c[ich];
    }
  }
  if (sumza < 0. || sumzc < 0.) return;
  // normalize Q vector ZNA and ZNC
  qza_plain.x = qza_plain.x / qza_plain.sum;
  qza_plain.y = qza_plain.y / qza_plain.sum;
  qzc_plain.x = qzc_plain.x / qzc_plain.sum;
  qzc_plain.y = qzc_plain.y / qzc_plain.sum;
  // recentering
  auto qza = qza_plain;
  auto qzc = qzc_plain;
  fQXACent->Fill(cent_v0m, qza.x);
  fQYACent->Fill(cent_v0m, qza.y);
  fQXCCent->Fill(cent_v0m, qzc.x);
  fQYCCent->Fill(cent_v0m, qzc.y);
  unsigned int correction_step = 0;
  // step 1: centrality
  if (fMeanQXACent && fMeanQYACent && fMeanQXCCent && fMeanQYCCent) {
    Recenter(qza, fMeanQXACent, fMeanQYACent, cent);
    Recenter(qzc, fMeanQXCCent, fMeanQYCCent, cent);
    fQXAvXvY->Fill(vtxx, vtxy, qza.x);
    fQYAvXvY->Fill(vtxx, vtxy, qza.y);
    fQXCvXvY->Fill(vtxx, vtxy, qzc.x);
    fQYCvXvY->Fill(vtxx, vtxy, qzc.y);
    correction_step = 1;
  }
  // step 2: vertex x and vertex y
  if (correction_step==1 && fMeanQXAvXvY && fMeanQYAvXvY && fMeanQXCvXvY && fMeanQYCvXvY) {
    Recenter(qza, fMeanQXAvXvY, fMeanQYAvXvY, vtxx, vtxy);
    Recenter(qzc, fMeanQXCvXvY, fMeanQYCvXvY, vtxx, vtxy);
    fQXACentVz->Fill(cent_v0m, vtxz, qza.x);
    fQYACentVz->Fill(cent_v0m, vtxz, qza.y);
    fQXCCentVz->Fill(cent_v0m, vtxz, qzc.x);
    fQYCCentVz->Fill(cent_v0m, vtxz, qzc.y);
    correction_step = 2;
  }
  // step 3: centrality and vertex z 
  if (correction_step==2 && fMeanQXACentVz && fMeanQYACentVz && fMeanQXCCentVz && fMeanQYCCentVz) {
    Recenter(qza, fMeanQXACentVz, fMeanQYACentVz, cent, vtxz);
    Recenter(qzc, fMeanQXCCentVz, fMeanQYCCentVz, cent, vtxz);
    correction_step = 3;
  }
  fCorrectionStep->Fill(correction_step);

  // Fill Correlation ZNA ZNC
  fXznaXzncCent->Fill(cent_v0m, qza.x*qzc.x);
  fYznaYzncCent->Fill(cent_v0m, qza.y*qzc.y);
  fXznaYzncCent->Fill(cent_v0m, qza.x*qzc.y);
  fYznaXzncCent->Fill(cent_v0m, qza.y*qzc.x);

  // track loop
  Qv qtpc1_plain = {0., 0., 0.};
  Qv qtpc2_plain = {0., 0., 0.};
  unsigned int ntracks = fEvent->GetNumberOfTracks();
  for (unsigned int i = 0; i < ntracks; ++i) {
    auto track = static_cast<AliAODTrack*>(fEvent->GetTrack(i));
    if (!track) continue;
    if (!track->TestFilterBit(fFilterBit)) continue;
    if (track->GetTPCNcls() < fNclsCut) continue;
    if (track->GetTPCchi2perCluster() < fChi2MinCut || track->GetTPCchi2perCluster() > fChi2MaxCut) continue;
    auto phi = track->Phi();
    auto pt = track->Pt();
    auto eta = track->Eta();
    auto sign = track->GetSign();
    if (std::abs(eta) < fEtaMax) continue;
    if (pt < fPtMin || pt > fPtMax) continue;
    if (fPositiveOnly && sign < 0) continue;
    if (fNegativeOnly && sign > 0) continue;

    qtpc1_plain.x += std::cos(phi);
    qtpc1_plain.y += std::sin(phi);
    qtpc2_plain.x += std::cos(2.*phi);
    qtpc2_plain.y += std::sin(2.*phi);
    qtpc1_plain.sum = qtpc1_plain.sum + 1;
    qtpc2_plain.sum = qtpc2_plain.sum + 1;
    
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
    
    fV1XXXpT[cent_bin]->Fill(pt, v2_xxx);
    fV1XYYpT[cent_bin]->Fill(pt, v2_xyy);
    fV1YXYpT[cent_bin]->Fill(pt, v2_yxy);
    fV1YYXpT[cent_bin]->Fill(pt, v2_yyx);
    fV1YYYpT[cent_bin]->Fill(pt, v2_yyy);
    fV1YXXpT[cent_bin]->Fill(pt, v2_yxx);
    fV1XXYpT[cent_bin]->Fill(pt, v2_xxy);
    fV1XYXpT[cent_bin]->Fill(pt, v2_xyx);

    fXtpcXznaXzncCent->Fill(cent_v0m, std::cos(2.phi)*qza.x*qzc.x);
    fXtpcXznaYzncCent->Fill(cent_v0m, std::cos(2.phi)*qza.x*qzc.y);
    fXtpcYznaXzncCent->Fill(cent_v0m, std::cos(2.phi)*qza.y*qzc.x);
    fYtpcXznaXzncCent->Fill(cent_v0m, std::sin(2.phi)*qza.x*qzc.x);
    fYtpcYznaYzncCent->Fill(cent_v0m, std::sin(2.phi)*qza.y*qzc.y);
    fYtpcYznaXzncCent->Fill(cent_v0m, std::sin(2.phi)*qza.y*qzc.x);
    fXtpcYznaYzncCent->Fill(cent_v0m, std::cos(2.phi)*qza.y*qzc.y);
    fYtpcXznaYzncCent->Fill(cent_v0m, std::sin(2.phi)*qza.x*qzc.y);

    fV1XXaEta[cent_bin_wide]->Fill(eta, v1_xxa);
    fV1XXcEta[cent_bin_wide]->Fill(eta, v1_xxc);
    fV1YYaEta[cent_bin_wide]->Fill(eta, v1_yya);
    fV1YYcEta[cent_bin_wide]->Fill(eta, v1_yyc);
    fV1XYaEta[cent_bin_wide]->Fill(eta, v1_xya);
    fV1XYcEta[cent_bin_wide]->Fill(eta, v1_xyc);
    fV1YXaEta[cent_bin_wide]->Fill(eta, v1_yxa);
    fV1YXcEta[cent_bin_wide]->Fill(eta, v1_yxc);

    if (eta > 0.) {
      fXtpcXznaCentEtaP->Fill(cent_v0m,v1_xxa);
      fXtpcXzncCentEtaP->Fill(cent_v0m,v1_xxc);
      fYtpcYznaCentEtaP->Fill(cent_v0m,v1_yya);
      fYtpcYzncCentEtaP->Fill(cent_v0m,v1_yya);
      fV1XXaPtEtaP[cent_bin_wide]->Fill(pt,v1_xxa); 
      fV1XXcPtEtaP[cent_bin_wide]->Fill(pt,v1_xxc); 
      fV1YYaPtEtaP[cent_bin_wide]->Fill(pt,v1_yya); 
      fV1YYcPtEtaP[cent_bin_wide]->Fill(pt,v1_yya); 
    } else {
      fXtpcXznaCentEtaN->Fill(cent_v0m,v1_xxa);
      fXtpcXzncCentEtaN->Fill(cent_v0m,v1_xxc);
      fYtpcYznaCentEtaN->Fill(cent_v0m,v1_yya);
      fYtpcYzncCentEtaN->Fill(cent_v0m,v1_yya);
      fV1XXaPtEtaN[cent_bin_wide]->Fill(pt,v1_xxa); 
      fV1XXcPtEtaN[cent_bin_wide]->Fill(pt,v1_xxc); 
      fV1YYaPtEtaN[cent_bin_wide]->Fill(pt,v1_yya); 
      fV1YYcPtEtaN[cent_bin_wide]->Fill(pt,v1_yya); 
    }

  }

  // normalize TPC Q vector
  qtpc1_plain.x = qtpc1_plain.x / qtpc1_plain.sum;
  qtpc1_plain.y = qtpc1_plain.y / qtpc1_plain.sum;
  qtpc2_plain.x = qtpc2_plain.x / qtpc2_plain.sum;
  qtpc2_plain.y = qtpc2_plain.y / qtpc2_plain.sum;

  auto psi_zna = TMath::Pi() + std::atan2(qzc.y, qza.x);
  auto psi_znc = TMath::Pi() + std::atan2(qzc.y, qzc.x);
  auto psi_tpc1 = TMath::Pi() + std::atan2(qtpc1_plain.y, qtpc1_plain.x);
  auto psi_tpc2 = TMath::Pi() + std::atan2(qtpc2_plain.y, qtpc2_plain.x);

  fPsiZA->Fill(psi_zna);
  fPsiZC->Fill(psi_znc);
  fPsiTPC1->Fill(psi_tpc1);
  fPsiTPC2->Fill(psi_tpc2);
}

void Recenter(Qv &q, TH1 *cx, TH1 *cy,  float coord_x, float coord_y = 0.) {
  auto ibin = cx->FindBin(coord_x, coord_y);
  q.x = (q.x - cx->GetBinContent(ibin)) / cx->GetBinError(ibin);
  q.y = (q.y - cy->GetBinContent(ibin)) / cy->GetBinError(ibin);
}

void AliAnalysisTaskFlowZDC::Terminate(Option_t *option) {
}

void AliAnalysisTaskFlowZDC::OpenCorrections(std::string file_name, int run) {
  auto file = TFile::Open(name.data());
  // oadb containing first step
  auto oadb_xac = dynamic_cast<AliOADBContainer*>(file->Get("fMeanQxZaCent"));
  if (oadb) fMeanQXZACent = dynamic_cast<TProfile*>(oadb->GetObject(run));
  auto oadb_yac = dynamic_cast<AliOADBContainer*>(file->Get("fMeanQyZaCent"));
  if (oadb) fMeanQYZACent = dynamic_cast<TProfile*>(oadb->GetObject(run));
  auto oadb_xcc = dynamic_cast<AliOADBContainer*>(file->Get("fMeanQxZcCent"));
  if (oadb) fMeanQXZCCent = dynamic_cast<TProfile*>(oadb->GetObject(run));
  auto oadb_ycc = dynamic_cast<AliOADBContainer*>(file->Get("fMeanQyZcCent"));
  if (oadb) fMeanQYZCCent = dynamic_cast<TProfile*>(oadb->GetObject(run));
  // oadb containing first step
  auto oadb_xaxy = dynamic_cast<AliOADBContainer*>(file->Get("fMeanQxZaVxVy"));
  if (oadb) fMeanQXZAvXvY = dynamic_cast<TProfile2D*>(oadb->GetObject(run));
  auto oadb_yaxy = dynamic_cast<AliOADBContainer*>(file->Get("fMeanQyZaVxVy"));
  if (oadb) fMeanQYZAvXvY = dynamic_cast<TProfile2D*>(oadb->GetObject(run));
  auto oadb_xcxy = dynamic_cast<AliOADBContainer*>(file->Get("fMeanQxZcVxVy"));
  if (oadb) fMeanQXZCvXvY = dynamic_cast<TProfile2D*>(oadb->GetObject(run));
  auto oadb_ycxy = dynamic_cast<AliOADBContainer*>(file->Get("fMeanQyZcVxVy"));
  if (oadb) fMeanQYZCvXvY = dynamic_cast<TProfile2D*>(oadb->GetObject(run));
  // oadb containing first step
  auto oadb_xacz = dynamic_cast<AliOADBContainer*>(file->Get("fMeanQxZaCentVz"));
  if (oadb) fMeanQXZACentVz = dynamic_cast<TProfile2D*>(oadb->GetObject(run));
  auto oadb_yacz = dynamic_cast<AliOADBContainer*>(file->Get("fMeanQyZaCentVz"));
  if (oadb) fMeanQYZACentVz = dynamic_cast<TProfile2D*>(oadb->GetObject(run));
  auto oadb_xccz = dynamic_cast<AliOADBContainer*>(file->Get("fMeanQxZcCentVz"));
  if (oadb) fMeanQXZCCentVz = dynamic_cast<TProfile2D*>(oadb->GetObject(run));
  auto oadb_yccz = dynamic_cast<AliOADBContainer*>(file->Get("fMeanQyZcCentVz"));
  if (oadb) fMeanQYZCCentVz = dynamic_cast<TProfile2D*>(oadb->GetObject(run));
}
