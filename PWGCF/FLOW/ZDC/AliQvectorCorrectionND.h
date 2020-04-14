#ifndef ALIQVECTORCORRECTIONND_H
#define ALIQVECTORCORRECTIONND_H

/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <string>
#include <vector>

#include "THn.h"
#include "TFile.h"
#include "TList.h"

#include "AliQvector.h"

class AliQvectorCorrectionND : public TObject {
 public:
  void Configure(std::string name,
                 const std::vector<TAxis*> &axes,
                 const std::vector<std::string> &rbr_axes,
                 bool equalize,
                 int min_entries);
  void Make(TFile *file, int run_number, TList *corrections, TList *qa);
  AliQvector Apply(const AliQvector q_vector, double *variables);

 private:
  bool fEqualization = false; 
  bool fIsApplied    = false;
  int fMinEntries = 0;
  std::string fName;
  std::vector<std::string> fRunByRunAxes;
  TList fAxes;

  THnF *fSumXout; //!<! Correction histogram output
  THnF *fSumYout; //!<! Correction histogram output
  THnF *fSumWout; //!<! Correction histogram output
  std::unique_ptr<THnF> fSumXin; //!<! Correction histogram input
  std::unique_ptr<THnF> fSumYin; //!<! Correction histogram input
  std::unique_ptr<THnF> fSumWin; //!<! Correction histogram input

  THnF *fSumXqa; //!<! Correction histogram output QA
  THnF *fSumYqa; //!<! Correction histogram output QA
  THnF *fSumWqa; //!<! Correction histogram output QA
  THnF *fEntriesQA; //!<! Correction histogram output QA fills events

  void OpenAxes(TFile *clist, int run_number);
  void OpenCorrection(TFile *file, int run_number);
  void AddHistogramsToList(TList *clist);
  void AddQAToList(TList *list);

  void ConfigureTHn(THnF **thn, std::string name, std::string title);

  template<typename T>
  T* ReadFromOADB(TFile *file, const std::string &oadb_name, int run_number) { 
    auto oadb = dynamic_cast<AliOADBContainer*>(file->Get(oadb_name.c_str()));
    if (oadb) return dynamic_cast<T*>(oadb->GetObject(run_number));
    else      return nullptr;
  }

 /// \cond CLASSDEF
 ClassDef(AliQvectorCorrectionND, 2);
 /// \endcond
};

#endif
