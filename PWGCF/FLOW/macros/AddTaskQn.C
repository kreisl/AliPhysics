AliAnalysisTaskQn* AddTaskQn(TString name = "Qn") {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  // resolve the name of the output file
  TString fileName = AliAnalysisManager::GetCommonFileName();
  fileName += ":Qn";      // create a subfolder in the file
  // now we create an instance of your task
  AliAnalysisTaskQn* task = new AliAnalysisTaskQn(name.Data());
  // add your task to the manager
  mgr->AddTask(task);
  // connect the manager to your task
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  // same for the output
  mgr->ConnectOutput(task,1,mgr->CreateContainer("TestHistos", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  // important: return a pointer to your task
  return task;
}
