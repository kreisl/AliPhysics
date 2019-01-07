AliAnalysisTaskQn* AddTaskQn(TString name = "Qn") {
  using QnTask = AliAnalysisTaskQn;
  using VAR = AliAnalysisTaskQn::Variables;
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  // resolve the name of the output file
  TString fileName = AliAnalysisManager::GetCommonFileName();
  fileName += ":Qn";      // create a subfolder in the file
  // now we create an instance of your task
  AliAnalysisTaskQn* task = new AliAnalysisTaskQn(name.Data());
  task->SelectCollisionCandidates(AliVEvent::kAnyINT);
  TString calibname("test.root");
  std::cout << "adding stuff" << std::endl;
  task->SetCalibrationFile(calibname,AliAnalysisTaskQn::CalibFile::local);
  std::cout << "adding stuff" << std::endl;
  auto qnman = task->GetQnManager();
  qnman->AddVariable("CentralityV0M",VAR::kCentV0M,1);
  qnman->AddVariable("V0APhi",VAR::kV0AChPhi,32);
  qnman->AddVariable("V0ARing",VAR::kV0AChRing,32);
  qnman->AddVariable("V0AMult",VAR::kV0AChMult,32);
  qnman->AddVariable("V0CPhi",VAR::kV0CChPhi,32);
  qnman->AddVariable("V0CRing",VAR::kV0CChRing,32);
  qnman->AddVariable("V0CMult",VAR::kV0CChMult,32);
  qnman->AddVariable("ZDCCPhi",VAR::kZDCCChPhi+1,4);
  qnman->AddVariable("ZDCAPhi",VAR::kZDCAChPhi+1,4);
  qnman->AddVariable("ZDCCMult",VAR::kZDCCChMult+1,4);
  qnman->AddVariable("ZDCAMult",VAR::kZDCAChMult+1,4);
  qnman->AddVariable("TPCPhi",VAR::kPhi,1);
  qnman->AddVariable("TPCEta",VAR::kEta,1);
  qnman->AddVariable("TPCPt",VAR::kPt,1);
  qnman->SetEventVariable("CentralityV0M");
  qnman->AddEventCut({"CentralityV0M"}, [](const double &cent) { return 0 < cent && cent < 100; });

  // TPC 
  auto confTPC = [](Qn::DetectorConfiguration *config) { config->SetNormalization(Qn::QVector::Normalization::M); };
  qnman->AddDetector("TPC", Qn::DetectorType::TRACK, "TPCPhi", "Ones", {}, {2, 3, 4});
  qnman->AddCut("TPC", {"TPCEta"}, [](const double &eta) { return -0.8 < eta && eta < 0.8; });
  qnman->AddCut("TPC", {"TPCPt"}, [](const double &pt) { return pt > 0.2 && pt < 10.; });
  qnman->SetCorrectionSteps("TPC", confTPC);

  // ZDC
  auto confZDC = [](Qn::DetectorConfiguration *config) {
  config->SetNormalization(Qn::QVector::Normalization::M);
  auto recenter = new Qn::Recentering();
  config->AddCorrectionOnQnVector(recenter);
  auto channels = new bool[4]{true, true, true, true};
  config->SetChannelsScheme(channels);
  };
  // ZDC A
  qnman->AddDetector("ZDCA", Qn::DetectorType::CHANNEL, "ZDCAPhi", "ZDCAMult", {}, {1});
  qnman->SetCorrectionSteps("ZDCA", confZDC);
  // ZDC C
  qnman->AddDetector("ZDCC", Qn::DetectorType::CHANNEL, "ZDCCPhi", "ZDCCMult", {}, {1});
  qnman->SetCorrectionSteps("ZDCC", confZDC);
  //V0 
  auto confV0 = [](Qn::DetectorConfiguration *config) {
    config->SetNormalization(Qn::QVector::Normalization::M);
    auto recenter = new Qn::Recentering();
    recenter->SetApplyWidthEqualization(true);
    config->AddCorrectionOnQnVector(recenter);
    auto V0Channels = new bool[32];
    auto channelGroups = new int[32];
    for (int ich = 0; ich < 32; ich++) {
      V0Channels[ich] = true;
      channelGroups[ich] = ich/8;
    }
    config->SetChannelsScheme(V0Channels, channelGroups);
    auto equal = new Qn::GainEqualization();
    equal->SetEqualizationMethod(Qn::GainEqualization::Method::AVERAGE);
    equal->SetUseChannelGroupsWeights(true);
    config->AddCorrectionOnInputData(equal);
  };  
  // V0-A
  qnman->AddDetector("V0A", Qn::DetectorType::CHANNEL, "V0APhi", "V0AMult", {}, {2, 3});
  qnman->AddCut("V0A", {"V0AMult"},[](double &mult) { return mult > 0.;});
  qnman->SetCorrectionSteps("V0A", confV0);
  qnman->AddHisto1D("V0A", {{"V0AChannels", 32, 0, 32}}, "V0AMult");
  // V0-C
  qnman->AddDetector("V0C", Qn::DetectorType::CHANNEL, "V0CPhi", "V0CMult", {}, {2, 3});
  qnman->AddCut("V0C", {"V0CMult"},[](double &mult) { return mult > 0.;});
  qnman->SetCorrectionSteps("V0C", confV0);
  qnman->AddHisto1D("V0C", {{"V0CChannels", 32, 0, 32}}, "V0CMult");
  
  TString correctionFile("correction.root");
  // add your task to the manager
  mgr->AddTask(task);
  // connect the manager to your task
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  // same for the output
  mgr->ConnectOutput(task,QnTask::OutputSlot::QnCalibration,mgr->CreateContainer("CalibrationHistograms", TList::Class(), AliAnalysisManager::kOutputContainer, correctionFile.Data()));
  mgr->ConnectOutput(task,QnTask::OutputSlot::QnQA,mgr->CreateContainer("CalibrationQAHistograms", TList::Class(), AliAnalysisManager::kOutputContainer, correctionFile.Data()));
  mgr->ConnectOutput(task,QnTask::OutputSlot::DetectorQA,mgr->CreateContainer("DetectorQA", TList::Class(), AliAnalysisManager::kOutputContainer, correctionFile.Data()));
  mgr->ConnectOutput(task,QnTask::OutputSlot::QnTree,mgr->CreateContainer("tree", TTree::Class(), AliAnalysisManager::kOutputContainer, "tree.root"));
  // important: return a pointer to your task
  return task;
}
