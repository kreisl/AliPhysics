#ifndef ALIANALYSISTASKFLOWZDC_CXX
#define ALIANALYSISTASKFLOWZDC_CXX

#include <string>
#include "AliAnalysisTaskSE.h"

#include "AliEventCuts.h"

class TH1;
class TFile;
class AliAODEvent;
class TProfile;
class TProfile2D;
class AliOADBContainer;
class AliEventCuts;

class AliAnalysisTaskFlowZDC : public AliAnalysisTaskSE {
  static constexpr Int_t fNcentBins = 10;
  static constexpr Int_t fNcentBinsWide = 3;
 public:
  struct Qv {
    Double_t x;
    Double_t y;
    Double_t sum;
    void Update(double phi, double weight) {
      x += std::cos(phi) * weight;      
      y += std::sin(phi) * weight;
      sum += weight;
    }
    void Normalize() {
      if (sum > 0.) {
        x /= sum; 
        y /= sum;
      }
    }
  };
  AliAnalysisTaskFlowZDC();
  AliAnalysisTaskFlowZDC(const char*);
  virtual ~AliAnalysisTaskFlowZDC();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t*);
  virtual void NotifyRun();
  virtual void Terminate(Option_t*);

  void SetCorrectionFile(std::string name) { fInputFileName = name; }

  void SetPtMin(Double_t pt_min) { fPtMin = pt_min; }
  void SetPtMax(Double_t pt_max) { fPtMax = pt_max; }
  void SetEtaMax(Double_t eta_max) { fEtaMax = eta_max; }
  void SetNclsCut(Double_t n_cls) { fNclsCut = n_cls; }
  void SetChi2MinCut(Double_t chi2_min) { fChi2MinCut = chi2_min; }
  void SetChi2MaxCut(Double_t chi2_max) { fChi2MaxCut = chi2_max; }
  void SetVertexZCut(Double_t vtx_z_max) { fVtxZcut = vtx_z_max; }
  void SetFilterBit(UInt_t bit) { fFilterBit = bit; }
  void SetNegativeOnly() { fNegativeOnly = true; }
  void SetPositiveOnly() { fPositiveOnly = true; }

 private:
  void RunAnalysis(AliAODEvent* event);
  void Recenter(Qv &q, TH1 *cx, TH1 *cy, Double_t coord_x, Double_t coord_y);
  void RecenterZN(Qv &qa, Qv &qc, Double_t cent, Double_t vtxx, Double_t vtxy, Double_t vtxz);
  void ReadFromOADB(TFile *file, const std::string &oadb_name, TProfile** profile);
  void ReadFromOADB(TFile *file, const std::string &oadb_name, TProfile2D** profile);
  bool CheckBins(TH1* histo);
  void OpenCorrections();
 private:
  UInt_t fFilterBit;
  Double_t fVtxZcut;
  Double_t fNclsCut;
  Double_t fChi2MinCut;
  Double_t fChi2MaxCut;
  Double_t fPtMin;
  Double_t fPtMax;
  Double_t fEtaMax;
  Bool_t fPositiveOnly;
  Bool_t fNegativeOnly;
  Int_t fRun;
  std::string fInputFileName;
  AliEventCuts fEventCuts;
  TList *fOutputList = nullptr;
  TList *fCorrectionList = nullptr;
  TList *fQAList = nullptr;
  TAxis *fCentralityAxis = nullptr;
  
  TProfile *fMeanXaCent = nullptr; //!<! Input to 1st Qx ZNA correction vs centrality 
  TProfile *fMeanYaCent = nullptr; //!<! Input to 1st Qy ZNA correction vs centrality
  TProfile *fMeanXcCent = nullptr; //!<! Input to 1st Qx ZNC correction vs centrality
  TProfile *fMeanYcCent = nullptr; //!<! Input to 1st Qy ZNC correction vs centrality
  
  TProfile2D *fMeanXavXvY = nullptr; //!<! Input to 2nd Qx ZNA correction vs vertex x + vertex y
  TProfile2D *fMeanYavXvY = nullptr; //!<! Input to 2nd Qy ZNA correction vs vertex x + vertex y
  TProfile2D *fMeanXcvXvY = nullptr; //!<! Input to 2nd Qx ZNC correction vs vertex x + vertex y
  TProfile2D *fMeanYcvXvY = nullptr; //!<! Input to 2nd Qy ZNC correction vs vertex x + vertex y

  TProfile2D *fMeanXaCentVz = nullptr; //!<! Input to 3rd Qx ZNA correction vs centrality + vertex z
  TProfile2D *fMeanYaCentVz = nullptr; //!<! Input to 3rd Qy ZNA correction vs centrality + vertex z
  TProfile2D *fMeanXcCentVz = nullptr; //!<! Input to 3rd Qx ZNC correction vs centrality + vertex z
  TProfile2D *fMeanYcCentVz = nullptr; //!<! Input to 3rd Qy ZNC correction vs centrality + vertex z
  
  TProfile *fXaCent = nullptr; //!<! Output of 1st Qx ZNA correction vs centrality 
  TProfile *fYaCent = nullptr; //!<! Output of 1st Qy ZNA correction vs centrality
  TProfile *fXcCent = nullptr; //!<! Output of 1st Qx ZNC correction vs centrality
  TProfile *fYcCent = nullptr; //!<! Output of 1st Qy ZNC correction vs centrality
  
  TProfile2D *fXavXvY = nullptr; //!<! Output of 2nd Qx ZNA correction vs vertex x + vertex y
  TProfile2D *fYavXvY = nullptr; //!<! Output of 2nd Qy ZNA correction vs vertex x + vertex y
  TProfile2D *fXcvXvY = nullptr; //!<! Output of 2nd Qx ZNC correction vs vertex x + vertex y
  TProfile2D *fYcvXvY = nullptr; //!<! Output of 2nd Qy ZNC correction vs vertex x + vertex y

  TProfile2D *fXaCentVz = nullptr; //!<! Output of 3rd Qx ZNA correction vs centrality + vertex z
  TProfile2D *fYaCentVz = nullptr; //!<! Output of 3rd Qy ZNA correction vs centrality + vertex z
  TProfile2D *fXcCentVz = nullptr; //!<! Output of 3rd Qx ZNC correction vs centrality + vertex z
  TProfile2D *fYcCentVz = nullptr; //!<! Output of 3rd Qy ZNC correction vs centrality + vertex z
  
  TProfile* fXaXcCent = nullptr; //!<! Correlation ZNA ZNC Resolution
  TProfile* fYaYcCent = nullptr; //!<! Correlation ZNA ZNC Resolution
  TProfile* fXaYcCent = nullptr; //!<! Correlation ZNA ZNC Resolution
  TProfile* fYaXcCent = nullptr; //!<! Correlation ZNA ZNC Resolution

  TProfile* fXtXaXcCent = nullptr; //!<! Correlation TPC ZNA ZNC
  TProfile* fXtYaYcCent = nullptr; //!<! Correlation TPC ZNA ZNC
  TProfile* fYtXaYcCent = nullptr; //!<! Correlation TPC ZNA ZNC
  TProfile* fYtYaXcCent = nullptr; //!<! Correlation TPC ZNA ZNC
  TProfile* fYtYaYcCent = nullptr; //!<! Correlation TPC ZNA ZNC
  TProfile* fXtXaYcCent = nullptr; //!<! Correlation TPC ZNA ZNC
  TProfile* fXtYaXcCent = nullptr; //!<! Correlation TPC ZNA ZNC
  TProfile* fYtXaXcCent = nullptr; //!<! Correlation TPC ZNA ZNC

  TProfile* fXtXaCentEtaP = nullptr; //!<! Correlation TPC ZNA ZNC
  TProfile* fXtXcCentEtaP = nullptr; //!<! Correlation TPC ZNA ZNC
  TProfile* fYtYaCentEtaP = nullptr; //!<! Correlation TPC ZNA ZNC
  TProfile* fYtYcCentEtaP = nullptr; //!<! Correlation TPC ZNA ZNC
  TProfile* fXtXaCentEtaN = nullptr; //!<! Correlation TPC ZNA ZNC
  TProfile* fXtXcCentEtaN = nullptr; //!<! Correlation TPC ZNA ZNC
  TProfile* fYtYaCentEtaN = nullptr; //!<! Correlation TPC ZNA ZNC
  TProfile* fYtYcCentEtaN = nullptr; //!<! Correlation TPC ZNA ZNC

  TProfile* fV1XXaEta[fNcentBinsWide]; //!<! v1 vs eta in wide centrality classes 
  TProfile* fV1XXcEta[fNcentBinsWide]; //!<! v1 vs eta in wide centrality classes 
  TProfile* fV1YYaEta[fNcentBinsWide]; //!<! v1 vs eta in wide centrality classes 
  TProfile* fV1YYcEta[fNcentBinsWide]; //!<! v1 vs eta in wide centrality classes 
  TProfile* fV1XYaEta[fNcentBinsWide]; //!<! v1 vs eta in wide centrality classes 
  TProfile* fV1XYcEta[fNcentBinsWide]; //!<! v1 vs eta in wide centrality classes 
  TProfile* fV1YXaEta[fNcentBinsWide]; //!<! v1 vs eta in wide centrality classes 
  TProfile* fV1YXcEta[fNcentBinsWide]; //!<! v1 vs eta in wide centrality classes 

  TProfile* fV1XXaPtEtaP[fNcentBinsWide]; //!<! v1 vs pt in wide centrality classes 
  TProfile* fV1XXcPtEtaP[fNcentBinsWide]; //!<! v1 vs pt in wide centrality classes 
  TProfile* fV1YYaPtEtaP[fNcentBinsWide]; //!<! v1 vs pt in wide centrality classes 
  TProfile* fV1YYcPtEtaP[fNcentBinsWide]; //!<! v1 vs pt in wide centrality classes 
  TProfile* fV1XXaPtEtaN[fNcentBinsWide]; //!<! v1 vs pt in wide centrality classes 
  TProfile* fV1XXcPtEtaN[fNcentBinsWide]; //!<! v1 vs pt in wide centrality classes 
  TProfile* fV1YYaPtEtaN[fNcentBinsWide]; //!<! v1 vs pt in wide centrality classes 
  TProfile* fV1YYcPtEtaN[fNcentBinsWide]; //!<! v1 vs pt in wide centrality classes 

  TProfile* fV2XXXpT[fNcentBins]; //!<! v2 vs pt in centrality classes 
  TProfile* fV2XYYpT[fNcentBins]; //!<! v2 vs pt in centrality classes
  TProfile* fV2YXYpT[fNcentBins]; //!<! v2 vs pt in centrality classes
  TProfile* fV2YYXpT[fNcentBins]; //!<! v2 vs pt in centrality classes
  TProfile* fV2YYYpT[fNcentBins]; //!<! v2 vs pt in centrality classes
  TProfile* fV2YXXpT[fNcentBins]; //!<! v2 vs pt in centrality classes
  TProfile* fV2XXYpT[fNcentBins]; //!<! v2 vs pt in centrality classes
  TProfile* fV2XYXpT[fNcentBins]; //!<! v2 vs pt in centrality classes

  TH1D* fPt = nullptr; //!<! pt QA histogram
  TH1D* fCentralityV0M = nullptr; //!<! centrality QA histogram
  TH1D* fCentralityCL1 = nullptr; //!<! centrality QA histogram
  TH1D* fPsiZA = nullptr; //!<! eventplane angle QA histogram
  TH1D* fPsiZC = nullptr; //!<! eventplane angle QA histogram
  TH1D* fPsiTPC1 = nullptr; //!<! eventplane angle QA histogram
  TH1D* fPsiTPC2 = nullptr; //!<! eventplane angle QA histogram
  TH1D* fVertexX = nullptr; //!<! primary vertex QA histogram
  TH1D* fVertexY = nullptr; //!<! primary vertex QA histogram
  TH1D* fVertexZ = nullptr;  //!<! primary vertex QA histogram
  TH1D* fCorrectionStep = nullptr;  //!<! recentering QA histogram

  TH2D* fCentralityCL1vsV0M = nullptr;

  /// \cond CLASSDEF
  ClassDef(AliAnalysisTaskFlowZDC, 1);
  /// \endcond
};

#endif
