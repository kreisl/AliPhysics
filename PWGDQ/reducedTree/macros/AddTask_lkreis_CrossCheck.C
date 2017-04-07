AliAnalysisTask *AddTask_lkreis_CrossCheck() {
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  if (!man) {
    return 0;
  }
  AliAnalysisTaskCrossCheck *task = new AliAnalysisTaskCrossCheck("lkreisTaskCrossCheck");
  
  AliAnalysisDataContainer *cinput = man->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = man->CreateContainer("QA_histograms",TList::Class(),AliAnalysisManager::kOutputContainer,man->GetCommonFileName());
  man->ConnectInput(task, 0, cinput);
  man->ConnectOutput(task, 1, coutput);
  return task;
}
