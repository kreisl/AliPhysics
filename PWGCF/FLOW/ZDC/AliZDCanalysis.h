#ifndef ALIZDCANALYSIS_H
#define ALIZDCANALYSIS_H

/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <string>
#include "TObject.h"
#include "AliZDCflowCuts.h"
#include "yaml-cpp/yaml.h"

class AliAODEvent;
class TList;

class AliZDCanalysis : public TObject {
  public:
   AliZDCanalysis() = default;
   AliZDCanalysis(std::string name) : fName(name) { }
   virtual ~AliZDCanalysis() = default;
   virtual TList *CreateCorrelations() = 0;
   virtual void FindCentralityBin(AliAODEvent *event, Double_t centrality, const std::vector<Int_t> &samples) = 0;
   AliZDCflowCuts& Cuts() { return fCuts; }
   void AddQAHistograms(TList* qa) {
     auto qalist = new TList();
     qalist->SetName((fName+"_QA").c_str());
     qalist->SetOwner(true);
     fCuts.AddQAHistograms(qalist);
     qa->Add(qalist);
   }
   inline void SetQvectors(std::map<std::string, AliZDCQvectors> &qv_map) {
     fQZA = qv_map[fQvectorName].za;
     fQZC = qv_map[fQvectorName].zc;
   }
   void ReadYAMLnode(YAML::Node& node) {
     auto analysis_settings = node["analysis_settings"];
     if (analysis_settings.IsDefined()) {
       ParseValue(fQvectorName, analysis_settings, "qvector_name");
     }
     fCuts.ReadYAMLnode(node);
   }
   void SetBootStrapSamples(int n) {
     fNsamples = n;
     fSamples.resize(n);
   }
  protected:
    Int_t fNsamples = 10;                                  /// Number of samples
    std::vector<Int_t> fSamples = std::vector<Int_t>(fNsamples);
    std::string fName = "default";
    std::string fQvectorName = "nd";
    Double_t fCentrality = -1.;
    AliZDCflowCuts fCuts;
    AliQvector fQZA; //!<! 
    AliQvector fQZC; //!<! 

    template<typename T>
    void ParseValue(T& value, YAML::Node node, std::string name) {
      auto value_node = node[name];
      if (value_node.IsDefined()) {
        value = value_node.as<T>();
      }
    }
    
    void ParseFlag(bool& value, YAML::Node node, std::string name) {
      auto value_node = node[name];
      if (value_node.IsDefined()) {
        if (value_node.as<int>() == 1) {
          value = true;
        } else {
          value = false;
        }
      }
    }
};

#endif
