#ifndef ALIQVECTORSRECENTERING_H
#define ALIQVECTORSRECENTERING_H

/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <string>
#include <vector>

#include "THn.h"
#include "TFile.h"
#include "TList.h"

#include "AliQvectors.h"

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
  constexpr unsigned int kNharmonics = 4;
  bool fEqualization = false; 
  bool fIsApplied    = false;
  int fMinEntries = 0;
  std::string fName;
  std::vector<std::string> fRunByRunAxes;
  TList fAxes;

  std::array<THnF *fSumXout, kNharmonics>; //!<! Correction histogram output
  std::array<THnF *fSumYout, kNharmonics>; //!<! Correction histogram output
  std::array<THnF *fSumWout, kNharmonics>; //!<! Correction histogram output
  std::array<std::unique_ptr<THnF> fSumXin, kNharmonics>; //!<! Correction histogram input
  std::array<std::unique_ptr<THnF> fSumYin, kNharmonics>; //!<! Correction histogram input
  std::array<std::unique_ptr<THnF> fSumWin, kNharmonics>; //!<! Correction histogram input

  std::array<THnF *fSumXqa, kNharmonics>; //!<! Correction histogram output QA
  std::array<THnF *fSumYqa, kNharmonics>; //!<! Correction histogram output QA
  std::array<THnF *fSumWqa, kNharmonics>; //!<! Correction histogram output QA
  std::array<THnF *fEntriesQA, kNharmonics>; //!<! Correction histogram output QA fills events

  void OpenAxes(TFile *clist, int run_number);
  void OpenCorrection(TFsle *file, int run_number);
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

