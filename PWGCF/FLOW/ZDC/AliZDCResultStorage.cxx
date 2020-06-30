#include <iostream>

#include <TFile.h>
#include <TObjArray.h>

#include "AliLog.h"
#include "AliZDCResultStorage.h"


/// \cond CLASSIMP
ClassImp(AliZDCResultStorage);
/// \endcond

AliZDCResultStorage::AliZDCResultStorage()
    : TNamed("AliZDCResult", "AliZDCResultStorage"), fObjects(new TObjArray()) {}

AliZDCResultStorage::AliZDCResultStorage(const TString &name)
    : TNamed(name.Data(), "AliZDCResultStorage"), fObjects(new TObjArray()) {}

AliZDCResultStorage::AliZDCResultStorage(const TString &name, TList *list)
    : TNamed(name.Data(), "AliZDCResultStorage"), fObjects(new TObjArray()) {
  fObjects->Add(list);
}

AliZDCResultStorage::~AliZDCResultStorage() { delete fObjects; }

Int_t AliZDCResultStorage::RecursiveDirectoryWrite(const TObject *objects, TDirectory &dir) const {
  Int_t result = 0;
  if (auto *collection = dynamic_cast<const TCollection *>(objects)) {
    TString path = collection->GetName();
    dir.mkdir(path);
    TDirectory *subdir = dir.GetDirectory(path);
    TIter next_object(collection);
    while (TObject *obj = next_object()) {
      result += AliZDCResultStorage::RecursiveDirectoryWrite(obj, *subdir);
    }
  } else {
    result += dir.WriteTObject(objects);
  }
  return result;
}

Int_t AliZDCResultStorage::Write(const char *name, Int_t option, Int_t bufsize) {
  return const_cast<const AliZDCResultStorage *>(this)->Write(name, option, bufsize);
}

Int_t AliZDCResultStorage::Write(const char *name, Int_t option, Int_t bufsize) const {
  TString path(name);
  // override if 'path' is empty
  if (path.IsWhitespace()) {
    path = fName.Strip(TString::kBoth);
  }
  gDirectory->mkdir(path);
  TDirectory *outdir = gDirectory->GetDirectory(path);
  if (!outdir) {
    Error("AliZDCResultStorage::Write", "Could not create path %s", path.Data());
    return 0;
  }
  Int_t result = 0;
  for (TObject *obj : *fObjects) {
    auto *output_object_list = dynamic_cast<TList *>(obj);
    if (!output_object_list) {
      AliWarning(Form("Unexpected type '%s' in output list (name: '%s'). Skipping.",
                      obj->ClassName(), obj->GetName()));
      continue;
    }
    for (TObject *output_object : *output_object_list) {
      result += RecursiveDirectoryWrite(output_object, *outdir);
    }
  }

  return result;
}