//====================================================================
// 
// Base class for classes that calculate the multiplicity in the
// central region event-by-event
// 
// Inputs: 
//   - AliESDEvent 
//
// Outputs: 
//   - AliAODCentralMult 
// 
// Histograms 
//   
// Corrections used 
#include "AliCentralMultiplicityTask.h"
#include "AliAODForwardMult.h"
#include "AliForwardUtil.h"
#include "AliLog.h"
#include "AliAODHandler.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliMultiplicity.h"
#include <TROOT.h>
#include <TFile.h>
#include <TError.h>
#include <iostream>
#include <iomanip>

//====================================================================
AliCentralMultiplicityTask::AliCentralMultiplicityTask(const char* name) 
  : AliAnalysisTaskSE(name),
    fInspector("centralEventInspector"),
    fData(0),
    fList(0),
    fAODCentral(kFALSE),
    fManager(),
    fUseSecondary(true),
    fUseAcceptance(true),
    fFirstEventSeen(false), 
    fIvz(0)
{
  // 
  // Constructor 
  //   
  DefineOutput(1, TList::Class());
  fBranchNames = 
    "ESD:AliESDRun.,AliESDHeader.,AliMultiplicity.,"
    "SPDVertex.,PrimaryVertex.";
}
//____________________________________________________________________
AliCentralMultiplicityTask::AliCentralMultiplicityTask() 
  : AliAnalysisTaskSE(),
    fInspector(),
    fData(0),
    fList(0),
    fAODCentral(),
    fManager(),
    fUseSecondary(true),
    fUseAcceptance(true),
    fFirstEventSeen(false), 
    fIvz(0)
{
  // 
  // Constructor 
  // 
}
//____________________________________________________________________
AliCentralMultiplicityTask::AliCentralMultiplicityTask(const AliCentralMultiplicityTask& o)
  : AliAnalysisTaskSE(o),
    fInspector(o.fInspector),
    fData(o.fData),
    fList(o.fList),
    fAODCentral(o.fAODCentral),
    fManager(o.fManager),
    fUseSecondary(o.fUseSecondary),
    fUseAcceptance(o.fUseAcceptance),
    fFirstEventSeen(o.fFirstEventSeen), 
    fIvz(0)
{
  //
  // Copy constructor 
  // 
}
//____________________________________________________________________
AliCentralMultiplicityTask&
AliCentralMultiplicityTask::operator=(const AliCentralMultiplicityTask& o)
{
  // 
  // Assignment operator 
  //
  fInspector      = o.fInspector;
  fData           = o.fData;
  fList           = o.fList;
  fAODCentral     = o.fAODCentral;
  fManager        = o.fManager;
  fUseSecondary   = o.fUseSecondary;
  fUseAcceptance  = o.fUseAcceptance;
  fFirstEventSeen = o.fFirstEventSeen;
  fIvz            = 0; 
  return *this;
}
//____________________________________________________________________
void AliCentralMultiplicityTask::UserCreateOutputObjects() 
{
  // 
  // Create output objects 
  // 
  //

  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler*      ah = 
    dynamic_cast<AliAODHandler*>(am->GetOutputEventHandler());
  if (!ah) AliFatal("No AOD output handler set in analysis manager");
  
  
  TObject* obj = &fAODCentral;
  ah->AddBranch("AliAODCentralMult", &obj);

    
  fList = new TList();
  fList->SetOwner();

  fInspector.DefineOutput(fList);

  PostData(1,fList);  
}

//____________________________________________________________________
AliESDEvent*
AliCentralMultiplicityTask::GetESDEvent()
{
  //
  // Get the ESD event. IF this is the first event, initialise
  //
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!esd) {
    AliWarning("No ESD event found for input event");
    return 0;
  }

  if (fFirstEventSeen) return esd;
  if (GetManager().IsInit()) return esd;
  
  fInspector.ReadRunDetails(esd);
  GetManager().Init(fInspector.GetCollisionSystem(),
		    fInspector.GetEnergy(),
		    fInspector.GetField());
  AliCentralCorrSecondaryMap* secMap = GetManager().GetSecMap();
  if (!secMap) 
    AliFatal("No secondary map defined!");
  const TAxis& vaxis = secMap->GetVertexAxis();

  std::cout << "Vertex range is " 
	    << vaxis.GetNbins() << "," 
	    << vaxis.GetXmin() << "," 
	    << vaxis.GetXmax()
	    << std::endl;
  fInspector.Init(vaxis);
  AliInfo("Manager of corrections in AliCentralMultiplicityTask init");
  fFirstEventSeen = kTRUE;
  Print();

  return esd;
}
//____________________________________________________________________
void
AliCentralMultiplicityTask::MarkEventForStore() const
{
  // Make sure the AOD tree is filled 
  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler*      ah = 
    dynamic_cast<AliAODHandler*>(am->GetOutputEventHandler());
  if (!ah)  
    AliFatal("No AOD output handler set in analysis manager");

  ah->SetFillAOD(kTRUE);
}

//____________________________________________________________________
void AliCentralMultiplicityTask::UserExec(Option_t* /*option*/) 
{
  // 
  // Process each event 
  // 
  // Parameters:
  //    option Not used
  //  
  fAODCentral.Clear("");
  fIvz = 0;

  AliESDEvent* esd = GetESDEvent();
  
  Bool_t   lowFlux   = kFALSE;
  UInt_t   triggers  = 0;
  UShort_t ivz       = 0;
  Double_t vz        = 0;
  Double_t cent      = -1;
  UShort_t nClusters = 0;
  UInt_t   found     = fInspector.Process(esd, triggers, lowFlux, 
					  ivz, vz, cent, nClusters);

  // No event or no trigger 
  if (found & AliFMDEventInspector::kNoEvent)    return;
  if (found & AliFMDEventInspector::kNoTriggers) return;
  
  // Make sure AOD is filled
  MarkEventForStore();

  if (found == AliFMDEventInspector::kNoSPD)      return;
  if (found == AliFMDEventInspector::kNoVertex)   return;
  if (triggers & AliAODForwardMult::kPileUp)      return;
  if (found == AliFMDEventInspector::kBadVertex)  return; // Out of range
  
  //Doing analysis
  fIvz = ivz;
  const AliMultiplicity* spdmult = esd->GetMultiplicity();

  TH2D& aodHist = fAODCentral.GetHistogram();

  ProcessESD(aodHist, spdmult);
  CorrectData(aodHist, ivz);

  PostData(1,fList);
}
//____________________________________________________________________
void 
AliCentralMultiplicityTask::ProcessESD(TH2D& aodHist, 
				       const AliMultiplicity* spdmult) const
{
 
  //Filling clusters in layer 1 used for tracklets...
  for(Int_t j = 0; j< spdmult->GetNumberOfTracklets();j++)
    aodHist.Fill(spdmult->GetEta(j),spdmult->GetPhi(j));

  //...and then the unused ones in layer 1 
  for(Int_t j = 0; j< spdmult->GetNumberOfSingleClusters();j++) 
    aodHist.Fill(-TMath::Log(TMath::Tan(spdmult->GetThetaSingle(j)/2.)),
		 spdmult->GetPhiSingle(j));
}

//____________________________________________________________________
void 
AliCentralMultiplicityTask::CorrectData(TH2D& aodHist, UShort_t vtxbin) const
{  
  // Corrections
  TH1D* hAcceptance = fManager.GetAcceptanceCorrection(vtxbin);
  TH2D* hSecMap     = fManager.GetSecMapCorrection(vtxbin);
  
  if (!hSecMap)     AliFatal("No secondary map!");
  if (!hAcceptance) AliFatal("No acceptance!");
    
  if (fUseSecondary && hSecMap) aodHist.Divide(hSecMap);
  
  for(Int_t nx = 1; nx <= aodHist.GetNbinsX(); nx++) {
    Float_t accCor = hAcceptance->GetBinContent(nx);
    Float_t accErr = hAcceptance->GetBinError(nx);
    
    Bool_t etabinSeen = kFALSE;  
    for(Int_t ny = 1; ny <= aodHist.GetNbinsY(); ny++) {
      // Get currrent value 
      Float_t aodValue = aodHist.GetBinContent(nx,ny);
      Float_t aodErr   = aodHist.GetBinError(nx,ny);

      // Set underflow bin
      Float_t secCor   = 0;
      if(hSecMap)       secCor     = hSecMap->GetBinContent(nx,ny);
      if (secCor > 0.5) etabinSeen = kTRUE;
      if (aodValue < 0.000001) { 
	aodHist.SetBinContent(nx,ny, 0); 
	continue; 
      }

      if (!fUseAcceptance) continue; 

      // Acceptance correction 
      if (accCor   < 0.000001) accCor = 1;
      Float_t aodNew   = aodValue / accCor ;
      Float_t error    = aodNew*TMath::Sqrt(TMath::Power(aodErr/aodValue,2) +
					    TMath::Power(accErr/accCor,2) );
      aodHist.SetBinContent(nx,ny, aodNew);
      //test
      aodHist.SetBinError(nx,ny,error);
      aodHist.SetBinError(nx,ny,aodErr);
      
    }
    //Filling underflow bin if we eta bin is in range
    if(etabinSeen) aodHist.SetBinContent(nx,0, 1.);
  }  
}

//____________________________________________________________________
void AliCentralMultiplicityTask::Terminate(Option_t* /*option*/) 
{
  // 
  // End of job
  // 
  // Parameters:
  //    option Not used 
  //
}
//____________________________________________________________________
void
AliCentralMultiplicityTask::Print(Option_t* option) const
{
  // 
  // Print information 
  // 
  // Parameters:
  //    option Not used
  //
  std::cout << ClassName() << ": " << GetName() << "\n" 
	    << std::boolalpha
	    << "  Use secondary correction:  " << fUseSecondary << '\n'
	    << "  Use acceptance correction: " << fUseAcceptance << '\n' 
	    << "  Off-line trigger mask:  0x" 
	    << std::hex     << std::setfill('0') 
	    << std::setw (8) << fOfflineTriggerMask 
	    << std::dec     << std::setfill (' ') 
	    << std::noboolalpha << std::endl;
  gROOT->IncreaseDirLevel();
  fManager.Print(option);
  fInspector.Print(option);
  gROOT->DecreaseDirLevel();
  
}
//====================================================================
AliCentralMultiplicityTask::Manager::Manager() :
  fAcceptancePath("$ALICE_ROOT/PWG2/FORWARD/corrections/CentralAcceptance"),
  fSecMapPath("$ALICE_ROOT/PWG2/FORWARD/corrections/CentralSecMap"),
  fAcceptance(),
  fSecmap(),
  fAcceptanceName("centralacceptance"),
  fSecMapName("centralsecmap"),
  fIsInit(kFALSE)
{
  //
  // Constructor 
  // 
}
//____________________________________________________________________
AliCentralMultiplicityTask::Manager::Manager(const Manager& o) 
  :fAcceptancePath(o.fAcceptancePath),
   fSecMapPath(o.fSecMapPath),
   fAcceptance(o.fAcceptance),
   fSecmap(o.fSecmap),
   fAcceptanceName(o.fAcceptanceName),
   fSecMapName(o.fSecMapName),
   fIsInit(o.fIsInit)
{
  //
  // Copy Constructor 
  // 
}
//____________________________________________________________________
AliCentralMultiplicityTask::Manager&
AliCentralMultiplicityTask::Manager::operator=(const Manager& o)
{
  //
  // Assignment operator  
  // 
  fAcceptancePath = o.fAcceptancePath;
  fSecMapPath     = o.fSecMapPath;
  fAcceptance     = o.fAcceptance;
  fSecmap         = o.fSecmap;
  fAcceptanceName = o.fAcceptanceName;
  fSecMapName     = o.fSecMapName;
  fIsInit         = o.fIsInit;
  return *this;
}

//____________________________________________________________________
const char* 
AliCentralMultiplicityTask::Manager::GetFullFileName(UShort_t what, 
						     UShort_t sys, 
						     UShort_t sNN,  
						     Short_t  field) const
{
  // 
  // Get full path name to object file 
  // 
  // Parameters:
  //    what   What to get 
  //    sys    Collision system
  //    sNN    Center of mass energy 
  //    field  Magnetic field 
  // 
  // Return:
  //    
  //
  return Form("%s/%s",
	      what == 0 ? GetSecMapPath() : GetAcceptancePath(), 
	      GetFileName(what, sys, sNN, field));
}

//____________________________________________________________________
const char* 
AliCentralMultiplicityTask::Manager::GetFileName(UShort_t  what ,
						 UShort_t  sys, 
						 UShort_t  sNN,
						 Short_t   field) const
{
  // 
  // Get the full path name 
  // 
  // Parameters:
  //    what   What to get
  //    sys    Collision system
  //    sNN    Center of mass energy 
  //    field  Magnetic field 
  // 
  // Return:
  //    
  //
  // Must be static - otherwise the data may disappear on return from
  // this member function
  static TString fname = "";
  
  fname = "";
  switch(what) {
  case 0:  fname.Append(fSecMapName.Data());     break;
  case 1:  fname.Append(fAcceptanceName.Data()); break;
  default:
    ::Error("GetFileName", 
	    "Invalid indentifier %d for central object, must be 0 or 1!", what);
    break;
  }
  fname.Append(Form("_%s_%04dGeV_%c%1dkG.root", 
		    AliForwardUtil::CollisionSystemString(sys), 
		    sNN, (field < 0 ? 'm' : 'p'), TMath::Abs(field)));
  
  return fname.Data();
}

//____________________________________________________________________
TH2D* 
AliCentralMultiplicityTask::Manager::GetSecMapCorrection(UShort_t vtxbin) const
{
  // 
  // Get the secondary map
  // 
  // Parameters:
  //    vtxbin 
  // 
  // Return:
  //    
  //
  if (!fSecmap) { 
    ::Warning("GetSecMapCorrection","No secondary map defined");
    return 0;
  }
  return fSecmap->GetCorrection(vtxbin);
}
//____________________________________________________________________
TH1D* 
AliCentralMultiplicityTask::Manager::GetAcceptanceCorrection(UShort_t vtxbin) 
  const 
{
  // 
  // Get the acceptance correction 
  // 
  // Parameters:
  //    vtxbin 
  // 
  // Return:
  //    
  //
  if (!fAcceptance) { 
    ::Warning("GetAcceptanceCorrection","No acceptance map defined");
    return 0;
  }
  return fAcceptance->GetCorrection(vtxbin);
}

//____________________________________________________________________
void 
AliCentralMultiplicityTask::Manager::Init(UShort_t  sys, 
					  UShort_t  sNN,
					  Short_t   field) 
{
  // 
  // Initialize 
  // 
  // Parameters:
  //    sys    Collision system (1: pp, 2: PbPb)
  //    sNN    Center of mass energy per nucleon pair [GeV]
  //    field  Magnetic field [kG]
  //
  if(fIsInit) ::Warning("Init","Already initialised - overriding...");
  
  TFile fsec(GetFullFileName(0,sys,sNN,field));
  fSecmap = 
    dynamic_cast<AliCentralCorrSecondaryMap*>(fsec.Get(fSecMapName.Data()));  
  if(!fSecmap) {
    ::Error("Init", "no central Secondary Map found!") ;
    return;
  }
  TFile facc(GetFullFileName(1,sys,sNN,field));
  fAcceptance = 
    dynamic_cast<AliCentralCorrAcceptance*>(facc.Get(fAcceptanceName.Data()));
  if(!fAcceptance) {
    ::Error("Init", "no central Acceptance found!") ;
    return;
  }
  
  if(fSecmap && fAcceptance) {
    fIsInit = kTRUE;
    ::Info("Init", 
	   "Central Manager initialised for %s, energy %dGeV, field %dkG",
	   sys == 1 ? "pp" : sys == 2 ? "PbPb" : "unknown", sNN,field);
  }  
}

//____________________________________________________________________
void 
AliCentralMultiplicityTask::Manager::Print(Option_t* option) const
{
  std::cout << " AliCentralMultiplicityTask::Manager\n" 
	    << std::boolalpha 
	    << "  Initialized:     " << fIsInit << '\n'
	    << "  Acceptance path: " << fAcceptancePath << '\n'
	    << "  Acceptance name: " << fAcceptanceName << '\n'
	    << "  Acceptance:      " << fAcceptance << '\n'
	    << "  Secondary path:  " << fSecMapPath << '\n'
	    << "  Secondary name:  " << fSecMapName << '\n'
	    << "  Secondary map:   " << fSecmap 
	    << std::noboolalpha << std::endl;
  if (fAcceptance) fAcceptance->Print(option);
  if (fSecmap)     fSecmap->Print(option);
}

//
// EOF
//
