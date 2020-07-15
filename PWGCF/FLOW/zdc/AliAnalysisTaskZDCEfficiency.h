#ifndef ALIANALYSISTASKZDCEFFICIENCY_H
#define ALIANALYSISTASKZDCEFFICIENCY_H

/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskSE.h"

#include "AliEventCuts.h"

class TList;
class TH1F;
class TH2F;
class AliAODEvent;
class AliMCEvent;
class AliAnalysisFilter;
class AliAnalysisUtils;

class AliAnalysisTaskZDCEfficiency : public AliAnalysisTaskSE  {
 public:
  AliAnalysisTaskZDCEfficiency();
  AliAnalysisTaskZDCEfficiency(const char *name);
  virtual ~AliAnalysisTaskZDCEfficiency();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
 private:
  bool fMC;
  AliEventCuts fEventCuts;           //!<!
  AliAnalysisUtils *fAnalysisUtils;  //!<!
  AliAODEvent *fEvent;  //!<!
  AliMCEvent *fMCEvent; //!<!
  TList *fHistograms;   //!<!
  TH2F *fHistPtCentrality96;   //!<!
  TH2F *fHistPtCentrality128;   //!<!
  TH2F *fHistPtCentrality768;   //!<!
  TH2F *fHistPtCentralityMC; //!<!
  TH1F *fHistCentrality;     //!<!

/// \cond CLASSDEF
ClassDef(AliAnalysisTaskZDCEfficiency, 2);
/// \endcond
};

#endif // ALIANALYSISTASKZDCEFFICIENCY_H