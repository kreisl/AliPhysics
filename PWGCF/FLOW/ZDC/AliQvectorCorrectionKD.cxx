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

#include <cmath>
#include <iostream>
#include <sstream>
#include "TFile.h"
#include "TH2.h"
#include "AliOADBContainer.h"
#include "AliQvectorCorrectionKD.h"

void AliQvectorCorrectionKD::Configure(const std::string &name, 
                                       const std::vector<std::string> &axes_names,
                                       int nbinsxy,
                                       int nbinsz,
                                       bool equalize) {
  fName = name;
  fAxesNames = axes_names;
  fDim = axes_names.size();
  fNbinsXY = nbinsxy;
  fNbinsZ = nbinsz;
  fIsConfigured = true;
}

AliQvector AliQvectorCorrectionKD::Correct(const AliQvector &q_vector,
                                           const std::vector<double> &vtx) const {
  AliQvector ret = q_vector;
  const auto zbin = fAxisZ->FindBin(vtx[2]) - 1; 
  if (zbin < 0 && zbin >= (int) fBinnings.size()) return ret;
  const auto xybin = fBinnings.at(zbin)->FindBin(vtx[0], vtx[1]) - 1;
  const auto bin = (xybin + fNbinsXY * zbin) + 1;
  const auto fbin = (double)bin - 0.5;
  fXmeanOut->Fill(fbin, q_vector.x);
  fYmeanOut->Fill(fbin, q_vector.y);
  fXmeanQADistributionxyz->Fill(fbin, q_vector.x);
  fYmeanQADistributionxyz->Fill(fbin, q_vector.y);
  fCorrectionEntries->Fill(bin);
  if (fXmeanIn && fYmeanIn) {
    const auto x_mean = fXmeanIn->GetBinContent(bin);
    const auto y_mean = fYmeanIn->GetBinContent(bin);
    const auto x_width = fXmeanIn->GetBinError(bin);
    const auto y_width = fYmeanIn->GetBinError(bin);
    ret = q_vector.Recenter(x_mean, y_mean, x_width, y_width, fEqualization);
    fXmeanQAxyz->Fill(fbin, ret.x);
    fYmeanQAxyz->Fill(fbin, ret.y);
  }
  return ret;
}

bool AliQvectorCorrectionKD::CheckBins(TH1 *histo) {
  for (int i = 0; i < histo->GetNbinsX(); ++i) {
    if (std::isnan(histo->GetBinContent(i)) || std::isnan(histo->GetBinError(i))) return false;
  }
  return true;
}

void AliQvectorCorrectionKD::OpenCorrection(TFile *file, Int_t run) {
  if (!file || file->IsZombie()) return;
  std::ostringstream varnames;
  varnames << fAxesNames[0] << "_" << fAxesNames[1] << "_" << fAxesNames[2];
  auto axis_name = varnames.str() + "_z_axis";
  fAxisZ.reset(ReadFromOADB<TAxis>(file, axis_name.c_str(), run));
  fBinnings.resize(fNbinsZ);
  for (auto i = 0; i < fNbinsZ; ++i) {
    auto binning_name = varnames.str() + "_poly_bin_" + std::to_string(i);
    fBinnings.at(i).reset(ReadFromOADB<TH2Poly>(file, binning_name.c_str(), run));
  }
  fXmeanIn.reset(ReadFromOADB<TProfile>(file, fXmeanOut->GetName(), run));
  fYmeanIn.reset(ReadFromOADB<TProfile>(file, fYmeanOut->GetName(), run));
  if (fXmeanIn && CheckBins(fXmeanIn.get()) && fYmeanIn && CheckBins(fYmeanIn.get())) { 
    fIsApplied = true; 
    std::cout << fName << " correction found in file" << file->GetName() << std::endl;
  }
}

void AliQvectorCorrectionKD::AddCorrectionsToList(TList *hlist, TList *qalist) {
  if (!fIsConfigured) return;
  std::ostringstream varnames;
  varnames << fAxesNames[0] << "_" << fAxesNames[1] << "_" << fAxesNames[2];
  auto x_name = fName + "_" +varnames.str() + "_x";
  auto y_name = fName + "_" +varnames.str() + "_y";
  auto axes_base_name = std::string(";correction bins;")+fName;
  auto x_axes = axes_base_name + " X";
  auto y_axes = axes_base_name + " Y";
  auto nbins = fNbinsXY * fNbinsZ;
  fXmeanOut = new TProfile(x_name.c_str(), x_axes.c_str(),
                           nbins, 0., (double) nbins, "s");
  fYmeanOut = new TProfile(y_name.c_str(), y_axes.c_str(), 
                           nbins, 0., (double) nbins, "s");
  hlist->Add(fXmeanOut);
  hlist->Add(fYmeanOut);
  auto correction_qa_list = new TList();
  correction_qa_list->SetName(fName.c_str());
  correction_qa_list->SetOwner(true);
  fCorrectionEntries = new TH1D((fName+"EntriesPerCorrectionBin").c_str(),
                                (axes_base_name+" Entries").c_str(), 
                                nbins, 0., (double) nbins);
  fXmeanQAxyz = new TProfile((x_name+"_after").c_str(), x_axes.c_str(), nbins, 0., (double) nbins, "s");
  fYmeanQAxyz = new TProfile((y_name+"_after").c_str(), y_axes.c_str(), nbins, 0., (double) nbins, "s");
  fXmeanQADistributionxyz = new TH2D((x_name+"_distribution").data(), x_axes.data(), nbins, 0., (double) nbins, 100, -1., 1.);
  fYmeanQADistributionxyz = new TH2D((y_name+"_distribution").data(), y_axes.data(), nbins, 0., (double) nbins, 100, -1., 1.);
  correction_qa_list->Add(fXmeanQAxyz);
  correction_qa_list->Add(fYmeanQAxyz);
  correction_qa_list->Add(fXmeanQADistributionxyz);
  correction_qa_list->Add(fYmeanQADistributionxyz);
  correction_qa_list->Add(fCorrectionEntries);
  qalist->Add(correction_qa_list);
}

