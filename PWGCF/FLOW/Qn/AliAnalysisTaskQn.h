#ifndef AliAnalysisTaskQn_cxx
#define AliAnalysisTaskQn_cxx

#include "AliAnalysisTaskSE.h"
#include "CorrectionManager.h"

class AliAODEvent;
class TList;

class AliAnalysisTaskQn : public AliAnalysisTaskSE
{
 public:
  AliAnalysisTaskQn();
  AliAnalysisTaskQn(const char *name);
  virtual ~AliAnalysisTaskQn();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option); 
  virtual void Terminate(Option_t *option); 
  /// \cond CLASSDEF
  ClassDef(AliAnalysisTaskQn, 1);
  /// \endcond
 private:
  Qn::CorrectionManager fQnManager; //!<! Qn correction manager
  AliAODEvent*  fInput;           //!<! input event
  TList*        fOutputList;    //!<! output list
  TH1F*         fTestHist;
};

#endif
