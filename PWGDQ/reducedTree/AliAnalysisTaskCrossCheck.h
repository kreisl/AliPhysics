#ifndef ALIANALYSISTASKCROSSCHECK

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskCrossCheck : public AliAnalysisTaskSE {
  public:
    AliAnalysisTaskCrossCheck();
    AliAnalysisTaskCrossCheck(const char *name);
    virtual ~AliAnalysisTaskCrossCheck();

    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t*);

  private:
    TList *fOutputContainer;
    TH1D *fHistDcaXY;
    TH1D *fHistDcaZ;
    TH1D *fHistDcaXYg;
    TH1D *fHistDcaZg;
    TH2D *fHistSumFmd;
 
    AliAnalysisTaskCrossCheck(const AliAnalysisTaskCrossCheck&);
    AliAnalysisTaskCrossCheck &operator = (const AliAnalysisTaskCrossCheck);
    ClassDef(AliAnalysisTaskCrossCheck,1);
};
#endif
