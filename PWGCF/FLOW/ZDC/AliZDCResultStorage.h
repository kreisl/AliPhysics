#ifndef ALIZDCRESULTSTORAGE_H
#define ALIZDCRESULTSTORAGE_H

#include <TList.h>
#include <TNamed.h>

/// \class AliZDCResultStorage
/// \brief Objects used as output-container
///
/// Storage is in the form of a TDirectory created upon calling `Write()`
///
class AliZDCResultStorage : public TNamed {
 public:
  /// Construct with defaults
  AliZDCResultStorage();

  /// Construct with name
  AliZDCResultStorage(const TString &name);

  /// Construct with name and List
  AliZDCResultStorage(const TString &name, TList *list);

  /// destroy the objects
  virtual ~AliZDCResultStorage();

  /// Add an output list to storage
  void AddOutputList(TList *list) { fObjects->AddLast(list); }

  /// Called upon by AliAnalysisManager to save data to container
  virtual Int_t Write(const char *name = nullptr, Int_t option = 0, Int_t bufsize = 0);
  virtual Int_t Write(const char *name = nullptr, Int_t option = 0, Int_t bufsize = 0) const;

 private:
  AliZDCResultStorage(AliZDCResultStorage const &) = delete;
  AliZDCResultStorage &operator=(AliZDCResultStorage const &) = delete;
  Int_t RecursiveDirectoryWrite(const TObject *objects, TDirectory &dir) const;
  TObjArray *fObjects;

  /// \cond CLASSIMP
  ClassDef(AliZDCResultStorage, 2);
  /// \endcond
};

#endif