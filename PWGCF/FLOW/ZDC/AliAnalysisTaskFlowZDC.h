#ifndef ALIANALYSISTASKFLOWZDC_CXX
#define ALIANALYSISTASKFLOWZDC_CXX

#include <string>
#include "AliAnalysisTaskSE.h"

class TH1;
class TFile;
class AliAODEvent;
class TProfile;
class TProfile2D;
class AliOADBContainer;


class AliAnalysisTaskFlowZDC : public AliAnalysisTaskSE {
 public:
  struct Qv {
    Double_t x;
    Double_t y;
    Double_t sum;
  }
  AliAnalysisTaskFlowZDC();
  AliAnalysisTaskFlowZDC(const char*);
  virtual ~AliAnalysisTaskFlowZDC();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t*);
  void RunAnalysis(AliAODEvent* event);
  void Recenter(Qv &q, TH1 *cx, TH1 *cy, float coord_x, float coord_y);
  virtual void Terminate(Option_t*);
 private:
  const Int_t fNcentBins = 10;
  const Int_t fNcentBinsWide = 3;
  UInt_t fFilterBit;
  Double_t fVtxZcut;
  Double_t fNclsCut
  Double_t fChi2MinCut;
  Double_t fChi2MaxCut;
  Double_t fPtMin;
  Double_t fPtMax;
  Double_t fEtaMax;
  Bool_t fPositiveOnly;
  Bool_t fNegativeOnly;
  Int_t fRun;
  TFile *fCorrectionInputFile  = nullptr;
  TFile *fOutputList = nullptr;

  TAxis *fCentralityAxis = nullptr;
  
  TProfile1D *fMeanQXZACent = nullptr; //! Input to 1st Qx ZNA correction vs centrality 
  TProfile1D *fMeanQYZACent = nullptr; //! Input to 1st Qy ZNA correction vs centrality
  TProfile1D *fMeanQXZCCent = nullptr; //! Input to 1st Qx ZNC correction vs centrality
  TProfile1D *fMeanQYZCCent = nullptr; //! Input to 1st Qy ZNC correction vs centrality
  
  TProfile2D *fMeanQXZAvXvY = nullptr; //! Input to 2nd Qx ZNA correction vs vertex x + vertex y
  TProfile2D *fMeanQYZAvXvY = nullptr; //! Input to 2nd Qy ZNA correction vs vertex x + vertex y
  TProfile2D *fMeanQXZCvXvY = nullptr; //! Input to 2nd Qx ZNC correction vs vertex x + vertex y
  TProfile2D *fMeanQYZCvXvY = nullptr; //! Input to 2nd Qy ZNC correction vs vertex x + vertex y

  TProfile2D *fMeanQXZACentVz = nullptr; //! Input to 3rd Qx ZNA correction vs centrality + vertex z
  TProfile2D *fMeanQYZACentVz = nullptr; //! Input to 3rd Qy ZNA correction vs centrality + vertex z
  TProfile2D *fMeanQXZCCentVz = nullptr; //! Input to 3rd Qx ZNC correction vs centrality + vertex z
  TProfile2D *fMeanQYZCCentVz = nullptr; //! Input to 3rd Qy ZNC correction vs centrality + vertex z
  
  TProfile1D *fQXZACent = nullptr; //! Output of 1st Qx ZNA correction vs centrality 
  TProfile1D *fQYZACent = nullptr; //! Output of 1st Qy ZNA correction vs centrality
  TProfile1D *fQXZCCent = nullptr; //! Output of 1st Qx ZNC correction vs centrality
  TProfile1D *fQYZCCent = nullptr; //! Output of 1st Qy ZNC correction vs centrality
  
  TProfile2D *fQXZAvXvY = nullptr; //! Output of 2nd Qx ZNA correction vs vertex x + vertex y
  TProfile2D *fQYZAvXvY = nullptr; //! Output of 2nd Qy ZNA correction vs vertex x + vertex y
  TProfile2D *fQXZCvXvY = nullptr; //! Output of 2nd Qx ZNC correction vs vertex x + vertex y
  TProfile2D *fQYZCvXvY = nullptr; //! Output of 2nd Qy ZNC correction vs vertex x + vertex y

  TProfile2D *fQXZACentVz = nullptr; //! Output of 3rd Qx ZNA correction vs centrality + vertex z
  TProfile2D *fQYZACentVz = nullptr; //! Output of 3rd Qy ZNA correction vs centrality + vertex z
  TProfile2D *fQXZCCentVz = nullptr; //! Output of 3rd Qx ZNC correction vs centrality + vertex z
  TProfile2D *fQYZCCentVz = nullptr; //! Output of 3rd Qy ZNC correction vs centrality + vertex z
  
  TProfile* fXznaXzncCent = nullptr; //! Correlation ZNA ZNC Resolution
  TProfile* fYznaYzncCent = nullptr; //! Correlation ZNA ZNC Resolution
  TProfile* fXznaYzncCent = nullptr; //! Correlation ZNA ZNC Resolution
  TProfile* fYznaXzncCent = nullptr; //! Correlation ZNA ZNC Resolution

  TProfile* fXtpcXznaXzncCent = nullptr; //! Correlation TPC ZNA ZNC
  TProfile* fXtpcYznaYzncCent = nullptr; //! Correlation TPC ZNA ZNC
  TProfile* fYtpcXznaYzncCent = nullptr; //! Correlation TPC ZNA ZNC
  TProfile* fYtpcYznaXzncCent = nullptr; //! Correlation TPC ZNA ZNC
  TProfile* fYtpcYznaYzncCent = nullptr; //! Correlation TPC ZNA ZNC
  TProfile* fXtpcXznaYzncCent = nullptr; //! Correlation TPC ZNA ZNC
  TProfile* fXtpcYznaXzncCent = nullptr; //! Correlation TPC ZNA ZNC

  TProfile* fXtpcXznaCentEtaP = nullptr; //! Correlation TPC ZNA ZNC
  TProfile* fXtpcXzncCentEtaP = nullptr; //! Correlation TPC ZNA ZNC
  TProfile* fYtpcYznaCentEtaP = nullptr; //! Correlation TPC ZNA ZNC
  TProfile* fYtpcYzncCentEtaP = nullptr; //! Correlation TPC ZNA ZNC
  TProfile* fXtpcXznaCentEtaN = nullptr; //! Correlation TPC ZNA ZNC
  TProfile* fXtpcXzncCentEtaN = nullptr; //! Correlation TPC ZNA ZNC
  TProfile* fYtpcYznaCentEtaN = nullptr; //! Correlation TPC ZNA ZNC
  TProfile* fYtpcYzncCentEtaN = nullptr; //! Correlation TPC ZNA ZNC

  TProfile* fV1XXaEta[fNcentBinsWide]; //! v1 vs eta in wide centrality classes 
  TProfile* fV1XXcEta[fNcentBinsWide]; //! v1 vs eta in wide centrality classes 
  TProfile* fV1YYaEta[fNcentBinsWide]; //! v1 vs eta in wide centrality classes 
  TProfile* fV1YYcEta[fNcentBinsWide]; //! v1 vs eta in wide centrality classes 
  TProfile* fV1XYaEta[fNcentBinsWide]; //! v1 vs eta in wide centrality classes 
  TProfile* fV1XYcEta[fNcentBinsWide]; //! v1 vs eta in wide centrality classes 
  TProfile* fV1YXaEta[fNcentBinsWide]; //! v1 vs eta in wide centrality classes 
  TProfile* fV1YXcEta[fNcentBinsWide]; //! v1 vs eta in wide centrality classes 

  TProfile* fV1XXaPtEtaP[fNcentBinsWide]; //! v1 vs pt in wide centrality classes 
  TProfile* fV1XXcPtEtaP[fNcentBinsWide]; //! v1 vs pt in wide centrality classes 
  TProfile* fV1YYaPtEtaP[fNcentBinsWide]; //! v1 vs pt in wide centrality classes 
  TProfile* fV1YYcPtEtaP[fNcentBinsWide]; //! v1 vs pt in wide centrality classes 
  TProfile* fV1XXaPtEtaN[fNcentBinsWide]; //! v1 vs pt in wide centrality classes 
  TProfile* fV1XXcPtEtaN[fNcentBinsWide]; //! v1 vs pt in wide centrality classes 
  TProfile* fV1YYaPtEtaN[fNcentBinsWide]; //! v1 vs pt in wide centrality classes 
  TProfile* fV1YYcPtEtaN[fNcentBinsWide]; //! v1 vs pt in wide centrality classes 

  TProfile* fV2XXXpT[fNcentBins]; //! v2 vs pt in centrality classes 
  TProfile* fV2XYYpT[fNcentBins]; //! v2 vs pt in centrality classes
  TProfile* fV2YXYpT[fNcentBins]; //! v2 vs pt in centrality classes
  TProfile* fV2YYXpT[fNcentBins]; //! v2 vs pt in centrality classes
  TProfile* fV2YYYpT[fNcentBins]; //! v2 vs pt in centrality classes
  TProfile* fV2YXXpT[fNcentBins]; //! v2 vs pt in centrality classes
  TProfile* fV2XXYpT[fNcentBins]; //! v2 vs pt in centrality classes
  TProfile* fV2XYXpT[fNcentBins]; //! v2 vs pt in centrality classes

  TH1D* fPt = nullptr; //! pt 
  TH1D* fCentralityV0M = nullptr;
  TH1D* fCentralityCL1 = nullptr;
  TH1D* fPsiZA = nullptr;
  TH1D* fPsiZC = nullptr;
  TH1D* fPsiTPC1 = nullptr;
  TH1D* fPsiTPC2 = nullptr;
  TH1D* fVertexX = nullptr;
  TH1D* fVertexY = nullptr;
  TH1D* fVertexZ = nullptr; 
  TH1D* fCorrectionStep = nullptr; 

  TH2D* fCentralityCL1vsV0M = nullptr;

  /// \cond CLASSDEF
  ClassDef(AliAnalysisTaskFlowZDC, 1);
  /// \endcond
};

#endif
