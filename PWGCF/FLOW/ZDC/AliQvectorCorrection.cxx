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
#include "TFile.h"
#include "TH2.h"
#include "AliOADBContainer.h"
#include "AliQvectorCorrection.h"

AliQvectorCorrection::AliQvectorCorrection(std::string name, bool equalize) :
  fEqualization(equalize),
  fName(name) {}

bool AliQvectorCorrection::CheckBins(TH1 *histo) {
  for (int i = 0; i < histo->GetNbinsX(); ++i) {
    if (std::isnan(histo->GetBinContent(i)) || std::isnan(histo->GetBinError(i))) return false;
  }
  return true;
}

void AliQvectorCorrection1D::Configure(std::string name, std::string var_name, std::vector<Double_t> bin_edges, bool equalize) {
  fName = name;
  fEqualization = equalize;
  fVarName = var_name;
  fBinEdges = bin_edges;
  fIsConfigured = true;
}

AliQvector AliQvectorCorrection1D::Correct(const AliQvector &q_vector, const double var_a, const double var_b) const {
  AliQvector ret = q_vector;
  fXmeanOut->Fill(var_a, q_vector.x);
  fYmeanOut->Fill(var_a, q_vector.y);
  fCorrectionEntries->Fill(var_a);
  fXDistributionQA->Fill(var_a, q_vector.x);
  fYDistributionQA->Fill(var_a, q_vector.y);
  if (fXmeanIn && fYmeanIn) {
    const auto bin = fXmeanIn->FindBin(var_a);
    const auto x_mean = fXmeanIn->GetBinContent(bin);
    const auto y_mean = fYmeanIn->GetBinContent(bin);
    const auto x_width = fXmeanIn->GetBinError(bin);
    const auto y_width = fYmeanIn->GetBinError(bin);
    ret = q_vector.Recenter(x_mean, y_mean, x_width, y_width, fEqualization);
    fXmeanQA->Fill(var_a, ret.x);
    fYmeanQA->Fill(var_a, ret.y);
  }
  return ret;
}

TProfile* AliQvectorCorrection1D::ReadFromOADB(TFile *file, const std::string &oadb_name, Int_t run_number) {
  TProfile *ret = nullptr;
  fIsApplied = false;
  auto oadb = dynamic_cast<AliOADBContainer*>(file->Get(oadb_name.c_str()));
  if (oadb) {
    auto profile = dynamic_cast<TProfile*>(oadb->GetObject(run_number));
    if (profile && profile->GetEntries() > 0. && CheckBins(profile)) {
      ret = profile;
      fIsApplied = true;
    }
  }
  return ret;
}

void AliQvectorCorrection1D::OpenCorrection(TFile *file, Int_t run_number) {
  if (!file || file->IsZombie()) return;
  fXmeanIn.reset(ReadFromOADB(file, fXmeanOut->GetName(), run_number));
  fYmeanIn.reset(ReadFromOADB(file, fYmeanOut->GetName(), run_number));
  if (fXmeanIn) std::cout << fName << " correction found in file" << file->GetName() << std::endl;
}

void AliQvectorCorrection1D::AddCorrectionsToList(TList *hlist, TList *qalist) {
  if (!fIsConfigured) return;
  auto x_name = fName+"_"+fVarName+"_x";
  auto y_name = fName+"_"+fVarName+"_y";
  auto axes_base_name = std::string(";")+fVarName+";"+fName;
  auto x_axes = axes_base_name+" X";
  auto y_axes = axes_base_name+" Y";
  auto bin_edges = fBinEdges;
  fXmeanOut = new TProfile(x_name.data(), x_axes.data(),
                           bin_edges.size()-1, bin_edges.data(), "s");
  fYmeanOut = new TProfile(y_name.data(), y_axes.data(), 
                           bin_edges.size()-1, bin_edges.data(), "s");
  hlist->Add(fXmeanOut);
  hlist->Add(fYmeanOut);
  auto correction_qa_list = new TList();
  correction_qa_list->SetName(fName.data());
  correction_qa_list->SetOwner(true);
  fCorrectionEntries = new TH1D((fName+"EntriesPerCorrectionBin").data(), (axes_base_name+" Entries").data(), 
                                bin_edges.size()-1, bin_edges.data());
  fXmeanQA = new TProfile((x_name+"_QA").data(), x_axes.data(),
                          bin_edges.size()-1, bin_edges.data(), "s");
  fYmeanQA = new TProfile((y_name+"_QA").data(), y_axes.data(), 
                           bin_edges.size()-1, bin_edges.data(), "s");
  fXDistributionQA = new TH2D((x_name+"_AfterQA").data(), x_axes.data(), bin_edges.size()-1, bin_edges.data(), 100, -1., 1.);
  fYDistributionQA = new TH2D((y_name+"_AfterQA").data(), y_axes.data(), bin_edges.size()-1, bin_edges.data(), 100, -1., 1.);
  correction_qa_list->Add(fCorrectionEntries);
  correction_qa_list->Add(fXDistributionQA);
  correction_qa_list->Add(fYDistributionQA);
  correction_qa_list->Add(fXmeanQA);
  correction_qa_list->Add(fYmeanQA);
  qalist->Add(correction_qa_list);
}

void AliQvectorCorrection1DInterpolate::Configure(std::string name, std::string var_name, std::vector<Double_t> bin_edges, bool equalize) {
  fName = name;
  fEqualization = equalize;
  fVarName = var_name;
  fBinEdges = bin_edges;
  fIsConfigured = true;
}

AliQvector AliQvectorCorrection1DInterpolate::Correct(const AliQvector &q_vector, const double var_a, const double var_b) const {
  AliQvector ret = q_vector;
  fXmeanOut->Fill(var_a, q_vector.x);
  fYmeanOut->Fill(var_a, q_vector.y);
  fCorrectionEntries->Fill(var_a);
  fXDistributionQA->Fill(var_a, q_vector.x);
  fYDistributionQA->Fill(var_a, q_vector.y);
  if (fXmeanIn && fYmeanIn) {
    const auto x_mean = fXmeanIn->Eval(var_a);
    const auto y_mean = fYmeanIn->Eval(var_a);
    const auto x_width = 1;
    const auto y_width = 1;
    ret = q_vector.Recenter(x_mean, y_mean, x_width, y_width, fEqualization);
    fXmeanQA->Fill(var_a, ret.x);
    fYmeanQA->Fill(var_a, ret.y);
  }
  return ret;
}

TF1* AliQvectorCorrection1DInterpolate::ReadFromOADB(TFile *file, const std::string &oadb_name, Int_t run_number) {
  TF1 *ret = nullptr;
  fIsApplied = false;
  auto oadb = dynamic_cast<AliOADBContainer*>(file->Get(oadb_name.c_str()));
  if (oadb) {
    auto profile = dynamic_cast<TProfile*>(oadb->GetObject(run_number));
    if (profile) {
      ret = dynamic_cast<TF1*>(profile->GetFunction("fit"));
      fIsApplied = true;
    }
  }
  return ret;
}

void AliQvectorCorrection1DInterpolate::OpenCorrection(TFile *file, Int_t run_number) {
  if (!file || file->IsZombie()) return;
  fXmeanIn.reset(ReadFromOADB(file, fXmeanOut->GetName(), run_number));
  fYmeanIn.reset(ReadFromOADB(file, fYmeanOut->GetName(), run_number));
  if (fXmeanIn) std::cout << fName << " correction found in file" << file->GetName() << std::endl;
}

void AliQvectorCorrection1DInterpolate::AddCorrectionsToList(TList *hlist, TList *qalist) {
  if (!fIsConfigured) return;
  auto x_name = fName+"_"+fVarName+"_x";
  auto y_name = fName+"_"+fVarName+"_y";
  auto axes_base_name = std::string(";")+fVarName+";"+fName;
  auto x_axes = axes_base_name+" X";
  auto y_axes = axes_base_name+" Y";
  auto bin_edges = fBinEdges;
  fXmeanOut = new TProfile(x_name.data(), x_axes.data(),
                           bin_edges.size()-1, bin_edges.data(), "s");
  fYmeanOut = new TProfile(y_name.data(), y_axes.data(), 
                           bin_edges.size()-1, bin_edges.data(), "s");
  hlist->Add(fXmeanOut);
  hlist->Add(fYmeanOut);
  auto correction_qa_list = new TList();
  correction_qa_list->SetName(fName.data());
  correction_qa_list->SetOwner(true);
  fCorrectionEntries = new TH1D((fName+"EntriesPerCorrectionBin").data(), (axes_base_name+" Entries").data(), 
                                bin_edges.size()-1, bin_edges.data());
  fXmeanQA = new TProfile((x_name+"_QA").data(), x_axes.data(),
                          bin_edges.size()-1, bin_edges.data(), "s");
  fYmeanQA = new TProfile((y_name+"_QA").data(), y_axes.data(), 
                           bin_edges.size()-1, bin_edges.data(), "s");
  fXDistributionQA = new TH2D((x_name+"_AfterQA").data(), x_axes.data(), bin_edges.size()-1, bin_edges.data(), 100, -1., 1.);
  fYDistributionQA = new TH2D((y_name+"_AfterQA").data(), y_axes.data(), bin_edges.size()-1, bin_edges.data(), 100, -1., 1.);
  correction_qa_list->Add(fCorrectionEntries);
  correction_qa_list->Add(fXDistributionQA);
  correction_qa_list->Add(fYDistributionQA);
  correction_qa_list->Add(fXmeanQA);
  correction_qa_list->Add(fYmeanQA);
  qalist->Add(correction_qa_list);
}

// 2D

void AliQvectorCorrection2D::Configure(std::string name, 
                                       std::string var_a_name, 
                                       std::vector<Double_t> var_a_bin_edges,
                                       std::string var_b_name, std::vector<Double_t> var_b_bin_edges,
                                       bool equalize) {
  fName = name;
  fEqualization = equalize;
  fNameVarA = var_a_name;
  fNameVarB = var_b_name;
  fBinEdgesVarA = var_a_bin_edges;
  fBinEdgesVarB = var_b_bin_edges;
}

AliQvector AliQvectorCorrection2D::Correct(const AliQvector &q_vector, const double var_a, const double var_b) const {
  AliQvector ret = q_vector;
  fXmeanOut->Fill(var_a, var_b, q_vector.x);
  fYmeanOut->Fill(var_a, var_b, q_vector.y);
  fCorrectionEntries->Fill(var_a, var_b);
  fXDistributionQAvsA->Fill(var_a, q_vector.x);
  fYDistributionQAvsA->Fill(var_a, q_vector.y);
  fXDistributionQAvsB->Fill(var_b, q_vector.x);
  fYDistributionQAvsB->Fill(var_b, q_vector.y);
  fXMeanQAvsA->Fill(var_a, q_vector.x);
  fYMeanQAvsA->Fill(var_a, q_vector.y);
  fXMeanQAvsB->Fill(var_b, q_vector.x);
  fYMeanQAvsB->Fill(var_b, q_vector.y);
  if (fXmeanIn && fYmeanIn) {
    const auto bin = fXmeanIn->FindBin(var_a, var_b);
    const auto x_mean = fXmeanIn->GetBinContent(bin);
    const auto y_mean = fYmeanIn->GetBinContent(bin);
    const auto x_width = fXmeanIn->GetBinError(bin);
    const auto y_width = fYmeanIn->GetBinError(bin);
    ret = q_vector.Recenter(x_mean, y_mean, x_width, y_width, fEqualization);
    fXmeanQA->Fill(var_a, var_b, ret.x);
    fYmeanQA->Fill(var_a, var_b, ret.y);
  } 
  return ret;
}
  
TProfile2D* AliQvectorCorrection2D::ReadFromOADB(TFile *file, const std::string &oadb_name, Int_t run_number) {
  TProfile2D* ret = nullptr;
  fIsApplied = false;
  AliOADBContainer* oadb = dynamic_cast<AliOADBContainer*>(file->Get(oadb_name.data()));
  if (oadb) {
    auto profile = dynamic_cast<TProfile2D*>(oadb->GetObject(run_number));
    if (profile && profile->GetEntries() > 0. && CheckBins(profile)) {
      ret = profile;
      fIsApplied = true;
    }
  }
  return ret;
}

void AliQvectorCorrection2D::OpenCorrection(TFile *file, Int_t run_number) {
  if (!file || file->IsZombie()) return;
  fXmeanIn.reset(ReadFromOADB(file, fXmeanOut->GetName(), run_number));
  fYmeanIn.reset(ReadFromOADB(file, fYmeanOut->GetName(), run_number));
  if (fXmeanIn) std::cout << fName << " correction found in file" << file->GetName() << std::endl;
}

void AliQvectorCorrection2D::AddCorrectionsToList(TList *hlist, TList *qalist) {
  if (!fIsConfigured) return;
  auto histo_base_name = fName+"_"+fNameVarA+"_"+fNameVarB;
  auto x_name = histo_base_name+"_x";
  auto y_name = histo_base_name+"_y";
  auto axes_base_name = std::string(";")+fNameVarA+";"+fNameVarB+";"+fName;
  auto x_axes = axes_base_name+" X";
  auto y_axes = axes_base_name+" Y";
  auto bin_edges_a = fBinEdgesVarA;
  auto bin_edges_b = fBinEdgesVarB;
  fXmeanOut = new TProfile2D(x_name.data(), x_axes.data(),
                             bin_edges_a.size()-1, bin_edges_a.data(),
                             bin_edges_b.size()-1, bin_edges_b.data(),
                             "s");
  fYmeanOut = new TProfile2D(y_name.data(), y_axes.data(), 
                             bin_edges_a.size()-1, bin_edges_a.data(),
                             bin_edges_b.size()-1, bin_edges_b.data(),
                             "s");
  hlist->Add(fXmeanOut);
  hlist->Add(fYmeanOut);

  auto correction_qa_list = new TList();
  correction_qa_list->SetName(fName.data());
  correction_qa_list->SetOwner(true);
  fCorrectionEntries = new TH2D((fName+"EntriesPerCorrectionBin").data(), (axes_base_name+" Entries").data(), 
                               bin_edges_a.size()-1, bin_edges_a.data(),
                               bin_edges_b.size()-1, bin_edges_b.data());
  fXmeanQA = new TProfile2D((x_name+"_AfterQA").data(), x_axes.data(),
                            bin_edges_a.size()-1, bin_edges_a.data(),
                            bin_edges_b.size()-1, bin_edges_b.data(),
                            "s");
  fYmeanQA = new TProfile2D((y_name+"_AfterQA").data(), y_axes.data(), 
                             bin_edges_a.size()-1, bin_edges_a.data(),
                             bin_edges_b.size()-1, bin_edges_b.data(),
                             "s");
  auto qa_a = fName+"_"+fNameVarA;
  auto qa_b = fName+"_"+fNameVarB;
  auto qa_axes_a = std::string(";")+fNameVarA+";"+fName;
  auto qa_axes_b = std::string(";")+fNameVarB+";"+fName;
  fXMeanQAvsA = new TProfile((qa_a+"_xmean_QA").data(), qa_axes_a.data(), bin_edges_a.size()-1, bin_edges_a.data(),"s");
  fYMeanQAvsA = new TProfile((qa_a+"_ymean_QA").data(), qa_axes_a.data(), bin_edges_a.size()-1, bin_edges_a.data(),"s");
  fXMeanQAvsB = new TProfile((qa_b+"_xmean_QA").data(), qa_axes_b.data(), bin_edges_b.size()-1, bin_edges_b.data(),"s");
  fYMeanQAvsB = new TProfile((qa_b+"_ymean_QA").data(), qa_axes_b.data(), bin_edges_b.size()-1, bin_edges_b.data(),"s");
  fXDistributionQAvsA = new TH2D((qa_a+"_x_QA").data(), (qa_axes_a+"_x").data(), bin_edges_a.size()-1, bin_edges_a.data(), 100, -1., 1);
  fYDistributionQAvsA = new TH2D((qa_a+"_y_QA").data(), (qa_axes_a+"_y").data(), bin_edges_a.size()-1, bin_edges_a.data(), 100, -1., 1);
  fXDistributionQAvsB = new TH2D((qa_b+"_x_QA").data(), (qa_axes_b+"_x").data(), bin_edges_b.size()-1, bin_edges_b.data(), 100, -1., 1);
  fYDistributionQAvsB = new TH2D((qa_b+"_y_QA").data(), (qa_axes_b+"_y").data(), bin_edges_b.size()-1, bin_edges_b.data(), 100, -1., 1);

  correction_qa_list->Add(fCorrectionEntries);
  correction_qa_list->Add(fXmeanQA);
  correction_qa_list->Add(fYmeanQA);
  correction_qa_list->Add(fXMeanQAvsA);
  correction_qa_list->Add(fYMeanQAvsA);
  correction_qa_list->Add(fXMeanQAvsB);
  correction_qa_list->Add(fYMeanQAvsB);
  correction_qa_list->Add(fXDistributionQAvsA);
  correction_qa_list->Add(fYDistributionQAvsA);
  correction_qa_list->Add(fXDistributionQAvsB);
  correction_qa_list->Add(fYDistributionQAvsB);
  qalist->Add(correction_qa_list);
}
