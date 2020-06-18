#ifndef ALIZDCFLOWCUT_H
#define ALIZDCFLOWCUT_H

/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <limits>
#include "yaml-cpp/yaml.h"

class AliAODTrack;
class AliAODEvent;

class AliZDCflowCuts : public TObject {
  public: 
    virtual ~AliZDCflowCuts() = default;
    void CheckEventCuts(AliAODEvent* event);
    bool CheckTrackCutsNoPtCut(AliAODTrack* track); 
    bool CheckTrackCutsPtCutOnly(AliAODTrack* track); 
    bool PassedEventCuts() const { return fEventSelected; }
    void HasPassedCentrality(bool flag) { fEventSelected &= flag; }
    void AddQAHistograms(TList *list);
    void ReadYAMLnode(YAML::Node& node) {
      auto eventcuts = node["eventcuts"];
      if (eventcuts.IsDefined()) {
        ParseValue(vtxZcut, eventcuts, "vtx_z");
      }
      auto trackcuts = node["trackcuts"];
      if (trackcuts.IsDefined()) {
      ParseValue(sign, trackcuts, "sign");
      ParseValue(filterBit, trackcuts, "filter_bit");
      ParseValue(nClsCut, trackcuts, "ncls_tpc");
      ParseValue(chi2MinCut, trackcuts, "chi2_min");
      ParseValue(chi2MaxCut, trackcuts, "chi2_max");
      ParseValue(ptMin, trackcuts, "pt_min");
      ParseValue(ptMax, trackcuts, "pt_max");
      ParseValue(etaMax, trackcuts, "eta_max");
      }
    }
    void Print();
    // event cuts
    Double_t vtxZcut    = 10.;
    // track cuts
    Int_t          sign =  0;
    UInt_t    filterBit =  768;
    Double_t    nClsCut =  70.;
    Double_t chi2MinCut =  0.1;
    Double_t chi2MaxCut =  4.;
    Double_t      ptMin =  0.2;
    Double_t      ptMax = 30.;
    Double_t     etaMax =  0.8;
    TH2D *fHistPhi;
    TH2D *fHistEta;
    TH2D *fHistDCAxy;
    TH2D *fHistDCAz;
    TH2D *fHistPt;
    TH1D *fHistTPCnCls;
    TH1D *fHistTPCchi2perCls;
    TH1D *fHistITSchi2;
    TH1D *fHistTPCchi2CvsGlo;
    TH1D *fHistTPCSharedClsF;
    bool fQAhistograms = false;
  private:
    template<typename T>
    void ParseValue(T& value, YAML::Node node, std::string name) {
      auto value_node = node[name];
      if (value_node.IsDefined()) {
        value = value_node.as<T>();
      }
    }
    bool fEventSelected = false;
    AliAODEvent *fEvent;

  /// \cond CLASSDEF
  ClassDef(AliZDCflowCuts, 3);
  /// \endcond
};

#endif
