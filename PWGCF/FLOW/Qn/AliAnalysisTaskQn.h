#ifndef AliAnalysisTaskQn_cxx
#define AliAnalysisTaskQn_cxx

#include <string>

#include "AliAnalysisTaskSE.h"
#include "CorrectionManager.h"
#include "AliQnLimits.h"

class AliLHCData;
class AliAODEvent;
class TList;

class AliAnalysisTaskQn : public AliAnalysisTaskSE
{
 public:
  static constexpr int kNFilterBits = 10;
  enum Variables {
   kNone,
	 kRunNumber,
	 kEventNumber,
   kTrigger,
   kTimeStamp,
   kPeriodNumber,
   kOrbitNumber,
   kBunchCrossNumber,
   kCentV0A,
   kCentV0C,
   kCentV0M,
   kCentZNC,
   kCentZNA,
   kCentCL0,
   kCentCL1,
   kNESDTracks,
   kNTracklets,
   kNITSClusters,
   kNTPCTracks = kNITSClusters + 6,
   kNTPCClusters,
   kV0Mult,
   kVtxX,
   kVtxY,
   kVtxZ,
   kZNSumEnergy,
   kFMDAPhi,
   kFMDAMult = kFMDAPhi + 1200,
   kFMDCPhi = kFMDAMult + 1200,
   kFMDCMult = kFMDCPhi + 1200,
   kV0CChMult = kFMDCMult + 1200,
   kV0CChPhi = kV0CChMult + 32,
   kV0CChRing = kV0CChPhi + 32,
   kV0AChMult = kV0CChRing + 32,
   kV0AChPhi = kV0AChMult + 32,
   kV0AChRing = kV0AChPhi + 32,
   kZDCCChMult = kV0AChRing + 32,
   kZDCCSumMult = kZDCCChMult + 5,
   kZDCCChPhi = kZDCCSumMult + 1,
   kZDCAChMult = kZDCCChPhi + 5,
   kZPAChMult = kZDCAChMult + 5,
   kZPCChMult = kZPAChMult + 5,
   kZPAChOffset = kZPCChMult + 5,
   kZPCChOffset = kZPAChOffset + 5,
   kZPAChPhi = kZPCChOffset + 5,
   kZPCChPhi = kZPAChPhi + 5,
   kZDCASumMult = kZPCChPhi + 5,
   kZDCAChPhi = kZDCASumMult + 1,
   kT0ChMult = kZDCAChPhi + 5,
   kT0ChPhi = kT0ChMult + 24,
   kNEventVariables = kT0ChPhi + 24,
   kPhi = kNEventVariables,
   kPt,
   kPx,
   kPy,
   kPz,
   kEta,
   kDCAxy,
   kDCAz,
   kDCAxySigma,
   kDCAzSigma,
   kCharge,
   kTPCnCls,
   kTPCchi2pCls,
   kFilterBits,
   kNVars = kFilterBits + kNFilterBits
  };

  enum CalibFile {
    local,
    alien,
    aodb
  };

  enum OutputSlot {
    QnCalibration = 1,
    QnQA,
    QnTree
  };
  AliAnalysisTaskQn();
  AliAnalysisTaskQn(const char*);
  virtual ~AliAnalysisTaskQn();
  virtual void UserCreateOutputObjects();
  virtual void NotifyRun();
  virtual void UserExec(Option_t*);
  virtual void Terminate(Option_t*);
  virtual void FinishTaskOutput();
  void SetCalibrationFile(TString, CalibFile);
  void SetOCDBPath(std::string path) { fOCDBPath = path;}
  void SetFilterBit(unsigned int filterbit) { fFilterBit = filterbit; }
  Qn::CorrectionManager* GetQnManager() {return fQnManager;}
  void RequireTimeCorrections() { fRequireTime = true; }
 private:
  bool          fRequireTime;
  AliQnLimits*  fQnLimits;           //!<!
  Qn::Axis<double>* fTimeAxis;       //!<!
  unsigned int  fFilterBit;          //!<! AOD track filter bit
  CalibFile     fCalibFileType;      //!<! type of calibration file input
  Qn::CorrectionManager* fQnManager; //!<! Qn correction manager
  AliAODEvent*  fEvent;              //!<! AOD input event
  TList*        fInDetectorList;     //!<! in detector list
  TTree*        fQnTree;             //!<! Qn output tree
  TFile*        fInCalib;            //!<! calibration input file
  double*       fValues;             //!<! values container 
  std::string   fOCDBPath;           //!<! path of OCDB files
  bool          fOCDBAvailable;      //!<! bool if OCDB are available
  AliLHCData*   fLHCData;            //!<! LHC data
	unsigned long fEventNumber;        //!<! LHC number of event in one run

  /// \cond CLASSDEF
  ClassDef(AliAnalysisTaskQn, 2);
  /// \endcond
};

#endif
