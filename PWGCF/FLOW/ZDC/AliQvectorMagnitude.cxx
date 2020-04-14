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

#include "AliQvectorMagnitude.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TList.h"

void AliQvectorMagnitude::Fill(const AliQvector &qvector, double centrality) {
  fDistributionOut->Fill(qvector.q(), centrality);
}

double AliQvectorMagnitude::GetPercentile(const AliQvector &qvector, double centrality) const {
  auto percentile = -1.;
  if (fDistributionIn) {
    auto ibin = fCentralityAxis->FindBin(centrality) - 1;
    if (ibin > -1 && ibin < fCentralityAxis->GetNbins()) {
      percentile = fSplines.at(ibin)->Eval(qvector.q());
    }
  }
  return percentile;
  
}

void AliQvectorMagnitude::Configure(const std::string &name, int nbins, double max) {
  fName = name+"QVectorMagnitude";
  fCentralityAxis = new TAxis(100, 0., 100.); 
  fDistributionOut = new TH2F((name+"_Distribution").c_str(),";#||{Q};P{Q}",
                           nbins, 0., max, 100, 0., 100.);
}

void AliQvectorMagnitude::Initialize() {
  for (int i = 1; i < fCentralityAxis->GetNbins()+1; ++i) {
    auto lo = fCentralityAxis->GetBinLowEdge(i);
    auto up = fCentralityAxis->GetBinUpEdge(i);
    auto projection = fDistributionIn->ProjectionX((fName+std::to_string(i)).c_str(), lo, up);
    auto percentiles = (TH1F*) projection->Clone((std::string("qpercent")+std::to_string(i)).c_str());
    for (int ib = 1; ib < projection->GetNbinsX()+1; ++ib) {
      double percentile = projection->Integral(1,ib)/projection->Integral();
      percentiles->SetBinContent(ib, percentile);
    }
    fSplines.emplace_back(new TSpline3(percentiles, "sp3"));
    fSplines.back()->SetName((fName+"Spline"+std::to_string(i)).c_str());
    delete percentiles;
    delete projection;
  }
}

void AliQvectorMagnitude::ReadFile(TFile *file) {
  if (file && !file->IsZombie()) {
  fDistributionIn = dynamic_cast<TH2F*>(file->Get(fDistributionOut->GetName()));
  if (fDistributionIn) Initialize();
  }
}

void AliQvectorMagnitude::AddToOutputList(TList *list) {
  list->Add(fDistributionOut);
}
