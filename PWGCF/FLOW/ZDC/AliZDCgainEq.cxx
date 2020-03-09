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

#include "TH1.h"
#include "AliOADBContainer.h"
#include "TFile.h"
#include "TList.h"
#include "AliZDCgainEq.h"
#include <iostream>

AliZDCgainEq::AliZDCgainEq(std::string name, UInt_t n_channels, std::vector<UInt_t> groups) :
  TObject(),
  fNChannels(n_channels),
  fName(name),
  fGroupEdges(groups) {
}

AliQvector AliZDCgainEq::ApplyGainEq(const double *channel_energies, const double *channel_phis) {
  AliQvector qvector = {0., 0., 0.};
  for (auto i = 0u; i < fNChannels; ++i) {
    fAverageGainOut->Fill(i, channel_energies[i]);
    fGainQAbefore->Fill(i, channel_energies[i]);
  }
  if (fAverageGainIn) {
    unsigned int i_group = 0;
    for (auto i = 0u; i < fNChannels; ++i) {
      if (i == i_group) i_group = i; 
      auto average = fAverageGainIn->GetBinContent(i+1);
      auto group_average = fAverageGainIn->GetBinContent(i_group+1);
      auto eq_channel_energy = channel_energies[i] / average * group_average;
      if (eq_channel_energy > 0.) {
        qvector.Update(channel_phis[i], eq_channel_energy);
      }
      fGainQAafter->Fill(i, eq_channel_energy);
    }
  }
  return qvector;
}

bool AliZDCgainEq::CheckBins(TH1 *histo) {
  for (int i = 0; i < histo->GetNbinsX(); ++i) {
    if (std::isnan(histo->GetBinContent(i)) || std::isnan(histo->GetBinError(i))) return false;
  }
  return true;
}

TProfile* AliZDCgainEq::ReadFromOADB(TFile *file, const std::string &oadb_name, Int_t run_number) {
  TProfile* ret = nullptr;
  fIsApplied = false;
  auto oadb = dynamic_cast<AliOADBContainer*>(file->Get(oadb_name.data()));
  if (oadb) {
    auto profile = dynamic_cast<TProfile*>(oadb->GetObject(run_number));
    if (profile && profile->GetEntries() > 0. && CheckBins(profile)) {
      ret = profile;
      fIsApplied = true;
    }
  }
  return ret;
}

void AliZDCgainEq::OpenCorrection(TFile *file, Int_t run_number) {
  if (!file || file->IsZombie()) return;
  fAverageGainIn.reset(ReadFromOADB(file, fAverageGainOut->GetName(), run_number));
}

void AliZDCgainEq::AddCorrectionsToList(TList *correction_list, TList *qa_list) {
  auto histo_name = fName+"_GainEQ";
  fAverageGainOut = new TProfile(histo_name.data(),";Channel ID; average Gain",fNChannels, 0, fNChannels);
  correction_list->Add(fAverageGainOut);
  auto correction_qa_list = new TList();
  correction_qa_list->SetName((fName+"_EQ").data());
  correction_qa_list->SetOwner(true);
  auto before_qa_name = fName+"_BeforeGainEQ";
  auto after_qa_name  = fName+"_AfterGainEQ";
  fGainQAbefore = new TProfile(before_qa_name.data(),";Channel ID; average equalized Gain",fNChannels, 0, fNChannels);
  fGainQAafter = new TProfile(after_qa_name.data(),";Channel ID; average equalized Gain",fNChannels, 0, fNChannels);
  correction_qa_list->Add(fGainQAbefore);
  correction_qa_list->Add(fGainQAafter);
  qa_list->Add(correction_qa_list);
}




