#ifndef ALIQNLIMITS_CXX
#define ALIQNLIMITS_CXX

#include "TNamed.h"
#include "TCollection.h"
#include <limits>

class AliQnLimits : public TNamed {
 public:
  AliQnLimits() = default;
  AliQnLimits(std::string name) : TNamed(name.data(),name.data()) {} 
  virtual ~AliQnLimits() = default;
  double Min() const {return min_;}
  double Max() const {return max_;}
  void SetNew(double lim) {
    if(lim > max_) {
      max_ = lim;
    }
    if(lim < min_) {
      min_ = lim;
    }
  }
  virtual Long64_t Merge(TCollection *list);
  private:
    double min_ = std::numeric_limits<double>::max();
    double max_ = std::numeric_limits<double>::min();

  public:
   //\cond CLASSIMP
   ClassDef(AliQnLimits, 2);
   // \endcond
};

#endif
