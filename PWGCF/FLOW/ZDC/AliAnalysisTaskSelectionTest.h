#ifndef ALIANALYSISTASKSELECTIONTEST_H
#define ALIANALYSISTASKSELECTIONTEST_H

/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskSelectionTest : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskSelectionTest();
  AliAnalysisTaskSelectionTest(const char*);
  virtual ~AliAnalysisTaskSelectionTest();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t*);

  TList *fHistogramList; 
  TH1D *fTPCnCls;

  /// \cond CLASSDEF
  ClassDef(AliAnalysisTaskSelectionTest, 2);
  /// \endcond
};
