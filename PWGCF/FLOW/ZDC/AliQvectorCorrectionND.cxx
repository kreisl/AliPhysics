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
#include "AliQvectorCorrectionND.h"

void AliQvectorCorrectionND::Configure(std::string name,
                                       const std::vector<TAxis*> &axes,
                                       const std::vector<std::string> &rbr_axes,
                                       bool equalize,
                                       int min_entries) {
  fName = name;
  fEqualization = equalize;
  fRunByRunAxes = rbr_axes;
  fAxes.SetOwner(true);
  for (auto &axis : axes) fAxes.Add(axis);
}

void AliQvectorCorrectionND::Make(TFile *file, int run_number, TList *corrections, TList *qa) {
  OpenAxes(file, run_number);
  AddHistogramsToList(corrections);
  AddQAToList(qa);
  OpenCorrection(file, run_number);
  if (fIsApplied) {
    AliDebug(AliLog::kInfo, (fName+": Applying recentering").c_str());
  } else {
    AliDebug(AliLog::kInfo, (fName+": Collecting recentering information").c_str());
  }
}

AliQvector AliQvectorCorrectionND::Apply(const AliQvector q_vector, double *variables) {
  AliQvector rec = q_vector;
  fSumXout->Fill(variables, q_vector.x);
  fSumYout->Fill(variables, q_vector.y);
  fSumWout->Fill(variables, 1.0);
  if (fSumWin && fSumXin && fSumYin) {
    std::vector<Int_t> bins(fSumWin->GetNdimensions());
    for (auto i = 0; i < fSumWin->GetNdimensions(); ++i) {
      bins[i] = fSumWin->GetAxis(i)->FindBin(variables[i]);
    }
    auto bin  = fSumWin->GetBin(bins.data());
    auto n    = fSumWin->GetBinContent(bin);
    auto xsum = fSumXin->GetBinContent(bin);
    auto ysum = fSumYin->GetBinContent(bin);
    if (n > fMinEntries) {
      auto xmean = xsum / n;
      auto ymean = ysum / n;
      auto xwidth = 1.;
      auto ywidth = 1.;
      if (fEqualization) {
        auto xerror = fSumXin->GetBinError2(bin);
        auto yerror = fSumYin->GetBinError2(bin);
        xwidth = sqrt(std::fabs(xerror/n - xmean*xmean));
        ywidth = sqrt(std::fabs(yerror/n - ymean*ymean));
      }
      rec = q_vector.Recenter(xmean, ymean, xwidth, ywidth, fEqualization);
      fSumXqa->Fill(variables, rec.x);
      fSumYqa->Fill(variables, rec.y);
      fSumWqa->Fill(variables, 1.0);
    } else {
      fEntriesQA->Fill(variables);
    }
  }
  return rec;
}


void AliQvectorCorrectionND::OpenAxes(TFile *file, int run_number) {
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

void AliQvectorCorrectionND::OpenCorrection(TFile *file, int run_number) {
  if (!file || file->IsZombie()) return;
  fSumXin.reset(ReadFromOADB<THnF>(file, fSumXout->GetName(), run_number));
  fSumYin.reset(ReadFromOADB<THnF>(file, fSumYout->GetName(), run_number));
  fSumWin.reset(ReadFromOADB<THnF>(file, fSumWout->GetName(), run_number));
  if (fSumXin && fSumYin && fSumWin) {
    fIsApplied = true;
    std::cout << fName << " ND recentering is applied" << std::endl;
  }
}

void AliQvectorCorrectionND::AddHistogramsToList(TList *list) {
  std::string axestitles;
  for (const auto &&obj : fAxes) axestitles += std::string(";") + obj->GetName();
  ConfigureTHn(&fSumXout, fName+"_X", axestitles);
  ConfigureTHn(&fSumYout, fName+"_Y", axestitles);
  ConfigureTHn(&fSumWout, fName+"_W", axestitles);
if (list) {
    list->Add(fSumXout);
    list->Add(fSumYout);
    list->Add(fSumWout);
  }
}

void AliQvectorCorrectionND::AddQAToList(TList *list) {
  std::string axestitles;
  for (const auto &&obj : fAxes) axestitles += std::string(";") + obj->GetName();
  ConfigureTHn(&fSumXqa, fName+"_X_QA", axestitles);
  ConfigureTHn(&fSumYqa, fName+"_Y_QA", axestitles);
  ConfigureTHn(&fSumWqa, fName+"_W_QA", axestitles);
  ConfigureTHn(&fEntriesQA, fName+"_Entries_QA", axestitles);
  if (list) {
    list->Add(fSumXqa);
    list->Add(fSumYqa);
    list->Add(fSumWqa);
    list->Add(fEntriesQA);
  }
}

void AliQvectorCorrectionND::ConfigureTHn(THnF **thn, 
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
