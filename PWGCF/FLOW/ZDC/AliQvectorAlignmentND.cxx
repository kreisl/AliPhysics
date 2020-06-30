/**************************************************************************
 * Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include <iostream>
#include "AliOADBContainer.h"
#include "AliQvectorAlignmentND.h"

void AliQvectorAlignmentND::Configure(std::string name,
                                       const std::vector<TAxis*> &axes,
                                       const std::vector<std::string> &rbr_axes,
                                       int min_entries) {
  fName = name;
  fRunByRunAxes = rbr_axes;
  fAxes.SetOwner(true);
  for (auto &axis : axes) fAxes.Add(axis);
}

void AliQvectorAlignmentND::Make(TFile *file, int run_number, TList *corrections, TList *qa) {
  OpenAxes(file, run_number);
  AddHistogramsToList(corrections);
  AddQAToList(qa);
  OpenCorrection(file, run_number);
}

AliQvector AliQvectorAlignmentND::Apply(const AliQvector q, const AliQvector qa, double *variables) {
  AliQvector aligned = q;
  if (!(qa.sum > 0. && q.sum > 0.)) return aligned;
  fSumXXout->Fill(variables, q.x*qa.x/qa.sum);
  fSumYYout->Fill(variables, q.y*qa.y/qa.sum);
  fSumXYout->Fill(variables, q.x*qa.y/qa.sum);
  fSumYXout->Fill(variables, q.y*qa.x/qa.sum);
  fSumWout->Fill(variables, 1.0);
  if (fSumWin && fSumXXin && fSumYYin && fSumXYin && fSumYXin) {
    std::vector<Int_t> bins(fSumWin->GetNdimensions());
    for (auto i = 0; i < fSumWin->GetNdimensions(); ++i) {
      bins[i] = fSumWin->GetAxis(i)->FindBin(variables[i]);
    }
    auto bin  = fSumWin->GetBin(bins.data());
    auto n    = fSumWin->GetBinContent(bin);
    auto xx = fSumXXin->GetBinContent(bin) / n;
    auto yy = fSumYYin->GetBinContent(bin) / n;
    auto xy = fSumXYin->GetBinContent(bin) / n;
    auto yx = fSumYXin->GetBinContent(bin) / n;
    auto exy = fSumXYin->GetBinError(bin) / n;
    auto eyx = fSumYXin->GetBinError(bin) / n;
    auto nharmonic = 1.;
    auto phi = -TMath::ATan2((xy-yx), (xx+yy))/nharmonic;
    auto significance = sqrt((xy-yx)*(xy-yx)/(exy*exy+eyx*eyx));
    if (n > fMinEntries && significance > 2.) {
      aligned.x = std::cos(phi) * q.x;  
      aligned.y = std::sin(phi) * q.y; 
      fSumXXqa->Fill(variables, aligned.x*qa.x);
      fSumYYqa->Fill(variables, aligned.y*qa.y);
      fSumXYqa->Fill(variables, aligned.x*qa.y);
      fSumYXqa->Fill(variables, aligned.y*qa.x);
      fSumWqa->Fill(variables, aligned.sum * qa.sum);
    } else {
      fEntriesQA->Fill(variables);
    }
  }
  return aligned;
}


void AliQvectorAlignmentND::OpenAxes(TFile *file, int run_number) {
  if (!file || file->IsZombie()) return;
  for (auto &name : fRunByRunAxes) {
    auto axis = dynamic_cast<TAxis*>(fAxes.FindObject(name.c_str()));
    auto axisconfig = ReadFromOADB<TAxis>(file, name.c_str(), run_number);
    if (axis && axisconfig) {
      axis->Set(axisconfig->GetNbins(), axisconfig->GetXbins()->GetArray());
    } else {
      std::cout << "false configuration " << name << std::endl;
    }
  }
}

void AliQvectorAlignmentND::OpenCorrection(TFile *file, int run_number) {
  if (!file || file->IsZombie()) return;
  fSumXXin.reset(ReadFromOADB<THnF>(file, fSumXXout->GetName(), run_number));
  fSumYYin.reset(ReadFromOADB<THnF>(file, fSumYYout->GetName(), run_number));
  fSumYXin.reset(ReadFromOADB<THnF>(file, fSumYXout->GetName(), run_number));
  fSumXYin.reset(ReadFromOADB<THnF>(file, fSumXYout->GetName(), run_number));
  fSumWin.reset(ReadFromOADB<THnF>(file, fSumWout->GetName(), run_number));
  if (fSumXXin && fSumYYin && fSumYXin && fSumXYin && fSumWin) {
    std::cout << fName << " ND align is applied" << std::endl;
    fIsApplied = true;
  }
}

void AliQvectorAlignmentND::AddHistogramsToList(TList *list) {
  std::string axestitles;
  for (const auto &&obj : fAxes) axestitles += std::string(";") + obj->GetName();
  ConfigureTHn(&fSumXXout, fName+"_XX", axestitles);
  ConfigureTHn(&fSumYYout, fName+"_YY", axestitles);
  ConfigureTHn(&fSumYXout, fName+"_YX", axestitles);
  ConfigureTHn(&fSumXYout, fName+"_XY", axestitles);
  ConfigureTHn(&fSumWout, fName+"_W", axestitles);
  if (list) {
    list->Add(fSumXXout);
    list->Add(fSumYYout);
    list->Add(fSumXYout);
    list->Add(fSumYXout);
    list->Add(fSumWout);
  }
}

void AliQvectorAlignmentND::AddQAToList(TList *list) {
  std::string axestitles;
  for (const auto &&obj : fAxes) axestitles += std::string(";") + obj->GetName();
  ConfigureTHn(&fSumXXqa, fName+"_XX_QA", axestitles);
  ConfigureTHn(&fSumYYqa, fName+"_YY_QA", axestitles);
  ConfigureTHn(&fSumXYqa, fName+"_XY_QA", axestitles);
  ConfigureTHn(&fSumYXqa, fName+"_XY_QA", axestitles);
  ConfigureTHn(&fSumWqa, fName+"_W_QA", axestitles);
  ConfigureTHn(&fEntriesQA, fName+"_Entries_QA", axestitles);
  if (list) {
    list->Add(fSumXXqa);
    list->Add(fSumYYqa);
    list->Add(fSumXYqa);
    list->Add(fSumYXqa);
    list->Add(fSumWqa);
    list->Add(fEntriesQA);
  }
}

void AliQvectorAlignmentND::ConfigureTHn(THnF **thn, 
                                          std::string name, 
                                          std::string title) {
    std::vector<int> nbins;
    std::vector<double> minvalues;
    std::vector<double> maxvalues;
    for (const auto &&obj : fAxes) {
      auto axis = dynamic_cast<TAxis*>(obj);
      nbins.push_back(axis->GetNbins()); 
      minvalues.push_back(axis->GetXmin()); 
      maxvalues.push_back(axis->GetXmax()); 
    }
    *thn = new THnF(name.c_str(), title.c_str(), fAxes.GetEntries(), 
                    nbins.data(), minvalues.data(), maxvalues.data());
    for (int iaxis = 0; iaxis < fAxes.GetEntries(); ++iaxis) {
      auto axis = dynamic_cast<TAxis*>(fAxes.At(iaxis));
      if (axis->IsVariableBinSize()) {
        (*thn)->GetAxis(iaxis)->Set(axis->GetNbins(), axis->GetXbins()->GetArray());
      } else {
        (*thn)->GetAxis(iaxis)->Set(axis->GetNbins(), axis->GetXmin(), axis->GetXmax());
      }
      (*thn)->GetAxis(iaxis)->SetTitle(fAxes.At(iaxis)->GetName());
    }
    (*thn)->Sumw2();
  }
