import FWCore.ParameterSet.Config as cms

process = cms.Process("TMVApplyWeights")

process.Options = cms.PSet(

   ## input file name
   InputFileName = cms.vstring(INPUTFILELIST),

   ## list of spectator variables
   InputSpectatorList = cms.vstring("pt","npu"),
   
   ## input tree name
   TreeName = cms.string("chs"),

   ## contains the information about which MVA method has been used
   Label    = cms.string("chs"),

   ## Lepton Type: Muon, Electron ,EleMu and Jets (fully hadronic)
   LeptonType = cms.string("Jets"),

   ## string which is used in the TMVATraining class to define a cut to be applied on the events
   PreselectionCutType = cms.string("basicJetsCutCSA14"),

   ## specify everything we need to take the right weight file
   InputVariableFileParam = cms.VPSet(),

   ## which jet read to fill the output : 0 means leading one, 1 ... etc
   JetToRead = cms.uint32(0),

   ## apply or not training cut when the MVA reading is used .. if yes, events out of the training phase space are put to a default value, otherwise we use the value
   ## given by the TMVAReader
   optionCut = cms.bool(False),

)

# single BDT, hight pT 

process.Options.InputVariableFileParam.insert(False, 
 cms.PSet(
  InputVariableList = cms.vstring("ecf_beta_10[0]"), 
  MethodName = cms.string("BDTG"), 
  JetPTRegion = cms.vdouble(475,600), 
  JetPURegion = cms.vdouble(0,39,39,100),
  inputWeightFiles = cms.vstring("/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/SingleVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU/TMVAWeight_BDTG_NoPruning_ECFbeta10_PTBin_475_600_PU_0_39/chs_ECFbeta10_PTBin_475_600_PU_0_39_BDTG_NoPruning.weights.xml",
                                "/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVAWeight_BDTG_NoPruning_ECFbeta10_PTBin_475_600_PU_39_100/chs_ECFbeta10_PTBin_475_600_PU_39_100_BDTG_NoPruning.weights.xml")))

process.Options.InputVariableFileParam.insert(False,
 cms.PSet(
  InputVariableList = cms.vstring("ecf_beta_15[0]"), 
  MethodName = cms.string("BDTG"), 
  JetPTRegion = cms.vdouble(475,600), 
  JetPURegion = cms.vdouble(0,39,39,100),
  inputWeightFiles = cms.vstring("/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/SingleVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU/TMVAWeight_BDTG_NoPruning_ECFbeta15_PTBin_475_600_PU_0_39/chs_ECFbeta15_PTBin_475_600_PU_0_39_BDTG_NoPruning.weights.xml",
                                "/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVAWeight_BDTG_NoPruning_ECFbeta15_PTBin_475_600_PU_39_100/chs_ECFbeta15_PTBin_475_600_PU_39_100_BDTG_NoPruning.weights.xml")))
 
process.Options.InputVariableFileParam.insert(False,
 cms.PSet(
  InputVariableList = cms.vstring("ecf_beta_20[0]"), 
  MethodName = cms.string("BDTG"), 
  JetPTRegion = cms.vdouble(475,600), 
  JetPURegion = cms.vdouble(0,39,39,100),
  inputWeightFiles = cms.vstring("/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/SingleVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU/TMVAWeight_BDTG_NoPruning_ECFbeta20_PTBin_475_600_PU_0_39/chs_ECFbeta20_PTBin_475_600_PU_0_39_BDTG_NoPruning.weights.xml",
                                "/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVAWeight_BDTG_NoPruning_ECFbeta20_PTBin_475_600_PU_39_100/chs_ECFbeta20_PTBin_475_600_PU_39_100_BDTG_NoPruning.weights.xml")))

process.Options.InputVariableFileParam.insert(False,
 cms.PSet(
  InputVariableList = cms.vstring("QGLikelihood_pr_zcut_010_R_cut_050[0]"), 
  MethodName = cms.string("BDTG"), 
  JetPTRegion = cms.vdouble(475,600), 
  JetPURegion = cms.vdouble(0,39,39,100),
  inputWeightFiles = cms.vstring("/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/SingleVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU/TMVAWeight_BDTG_NoPruning_QGLikelihood_PTBin_475_600_PU_0_39/chs_QGLikelihood_PTBin_475_600_PU_0_39_BDTG_NoPruning.weights.xml",
                                "/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVAWeight_BDTG_NoPruning_QGLikelihood_PTBin_475_600_PU_39_100/chs_QGLikelihood_PTBin_475_600_PU_39_100_BDTG_NoPruning.weights.xml")))

process.Options.InputVariableFileParam.insert(False,
 cms.PSet(
  InputVariableList = cms.vstring("QGLikelihood_pr_sub1_zcut_010_R_cut_050[0]"), 
  MethodName = cms.string("BDTG"), 
  JetPTRegion = cms.vdouble(475,600), 
  JetPURegion = cms.vdouble(0,39,39,100),
  inputWeightFiles = cms.vstring("/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/SingleVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU/TMVAWeight_BDTG_NoPruning_QGLikelihood_sub1_PTBin_475_600_PU_0_39/chs_QGLikelihood_sub1_PTBin_475_600_PU_0_39_BDTG_NoPruning.weights.xml",
                                "/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVAWeight_BDTG_NoPruning_QGLikelihood_sub1_PTBin_475_600_PU_39_100/chs_QGLikelihood_sub1_PTBin_475_600_PU_39_100_BDTG_NoPruning.weights.xml")))

process.Options.InputVariableFileParam.insert(False,
 cms.PSet(
  InputVariableList = cms.vstring("QGLikelihood_pr_sub2_zcut_010_R_cut_050[0]"), 
  MethodName = cms.string("BDTG"), 
  JetPTRegion = cms.vdouble(475,600), 
  JetPURegion = cms.vdouble(0,39,39,100),
  inputWeightFiles = cms.vstring("/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/SingleVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU/TMVAWeight_BDTG_NoPruning_QGLikelihood_sub2_PTBin_475_600_PU_0_39/chs_QGLikelihood_sub2_PTBin_475_600_PU_0_39_BDTG_NoPruning.weights.xml",
                                "/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVAWeight_BDTG_NoPruning_QGLikelihood_sub2_PTBin_475_600_PU_39_100/chs_QGLikelihood_sub2_PTBin_475_600_PU_39_100_BDTG_NoPruning.weights.xml")))

process.Options.InputVariableFileParam.insert(False,
 cms.PSet(
  InputVariableList = cms.vstring("2*QGLikelihood_pr_sub2_zcut_010_R_cut_050[0]+QGLikelihood_pr_sub1_zcut_010_R_cut_050[0]"), 
  MethodName = cms.string("BDTG"), 
  JetPTRegion = cms.vdouble(475,600), 
  JetPURegion = cms.vdouble(0,39,39,100),
  inputWeightFiles = cms.vstring("/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/SingleVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU/TMVAWeight_BDTG_NoPruning_QGLikelihood_comb_PTBin_475_600_PU_0_39/chs_QGLikelihood_comb_PTBin_475_600_PU_0_39_BDTG_NoPruning.weights.xml",
                                "/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVAWeight_BDTG_NoPruning_QGLikelihood_comb_PTBin_475_600_PU_39_100/chs_QGLikelihood_comb_PTBin_475_600_PU_39_100_BDTG_NoPruning.weights.xml")))

process.Options.InputVariableFileParam.insert(False,
 cms.PSet(
  InputVariableList = cms.vstring("Qjets[0]"), 
  MethodName = cms.string("BDTG"), 
  JetPTRegion = cms.vdouble(475,600), 
  JetPURegion = cms.vdouble(0,39,39,100),
  inputWeightFiles = cms.vstring("/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/SingleVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU/TMVAWeight_BDTG_NoPruning_Qjets_PTBin_475_600_PU_0_39/chs_Qjets_PTBin_475_600_PU_0_39_BDTG_NoPruning.weights.xml",
                                "/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVAWeight_BDTG_NoPruning_Qjets_PTBin_475_600_PU_39_100/chs_Qjets_PTBin_475_600_PU_39_100_BDTG_NoPruning.weights.xml")))

process.Options.InputVariableFileParam.insert(False,
 cms.PSet(
  InputVariableList = cms.vstring("mconst[0]"), 
  MethodName = cms.string("BDTG"), 
  JetPTRegion = cms.vdouble(475,600), 
  JetPURegion = cms.vdouble(0,39,39,100),
  inputWeightFiles = cms.vstring("/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/SingleVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU/TMVAWeight_BDTG_NoPruning_mconst_PTBin_475_600_PU_0_39/chs_mconst_PTBin_475_600_PU_0_39_BDTG_NoPruning.weights.xml",
                                "/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVAWeight_BDTG_NoPruning_mconst_PTBin_475_600_PU_39_100/chs_mconst_PTBin_475_600_PU_39_100_BDTG_NoPruning.weights.xml")))

process.Options.InputVariableFileParam.insert(False,
 cms.PSet(
  InputVariableList = cms.vstring("mprunedsafe_zcut_010_R_cut_050[0]"), 
  MethodName = cms.string("BDTG"), 
  JetPTRegion = cms.vdouble(475,600), 
  JetPURegion = cms.vdouble(0,39,39,100),
  inputWeightFiles = cms.vstring("/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/SingleVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU/TMVAWeight_BDTG_NoPruning_mpruned_PTBin_475_600_PU_0_39/chs_mpruned_PTBin_475_600_PU_0_39_BDTG_NoPruning.weights.xml",
                                "/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVAWeight_BDTG_NoPruning_mpruned_PTBin_475_600_PU_39_100/chs_mpruned_PTBin_475_600_PU_39_100_BDTG_NoPruning.weights.xml")))

process.Options.InputVariableFileParam.insert(False,
 cms.PSet(
  InputVariableList = cms.vstring("mraw[0]"), 
  MethodName = cms.string("BDTG"), 
  JetPTRegion = cms.vdouble(475,600), 
  JetPURegion = cms.vdouble(0,39,39,100),
  inputWeightFiles = cms.vstring("/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/SingleVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU/TMVAWeight_BDTG_NoPruning_mraw_PTBin_475_600_PU_0_39/chs_mraw_PTBin_475_600_PU_0_39_BDTG_NoPruning.weights.xml",
                                "/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVAWeight_BDTG_NoPruning_mraw_PTBin_475_600_PU_39_100/chs_mraw_PTBin_475_600_PU_39_100_BDTG_NoPruning.weights.xml")))

process.Options.InputVariableFileParam.insert(False,
 cms.PSet(
  InputVariableList = cms.vstring("msoftdropsafe_beta00[0]"), 
  MethodName = cms.string("BDTG"), 
  JetPTRegion = cms.vdouble(475,600), 
  JetPURegion = cms.vdouble(0,39,39,100),
  inputWeightFiles = cms.vstring("/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/SingleVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU/TMVAWeight_BDTG_NoPruning_msoftbeta00_PTBin_475_600_PU_0_39/chs_msoftbeta00_PTBin_475_600_PU_0_39_BDTG_NoPruning.weights.xml",
                                "/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVAWeight_BDTG_NoPruning_msoftbeta00_PTBin_475_600_PU_39_100/chs_msoftbeta00_PTBin_475_600_PU_39_100_BDTG_NoPruning.weights.xml")))

process.Options.InputVariableFileParam.insert(False,
 cms.PSet(
  InputVariableList = cms.vstring("msoftdropsafe_beta10[0]"), 
  MethodName = cms.string("BDTG"), 
  JetPTRegion = cms.vdouble(475,600), 
  JetPURegion = cms.vdouble(0,39,39,100),
  inputWeightFiles = cms.vstring("/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/SingleVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU/TMVAWeight_BDTG_NoPruning_msoftbeta10_PTBin_475_600_PU_0_39/chs_msoftbeta10_PTBin_475_600_PU_0_39_BDTG_NoPruning.weights.xml",
                                "/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVAWeight_BDTG_NoPruning_msoftbeta10_PTBin_475_600_PU_39_100/chs_msoftbeta10_PTBin_475_600_PU_39_100_BDTG_NoPruning.weights.xml")))

process.Options.InputVariableFileParam.insert(False,
 cms.PSet(
  InputVariableList = cms.vstring("msoftdropsafe_beta20[0]"), 
  MethodName = cms.string("BDTG"), 
  JetPTRegion = cms.vdouble(475,600), 
  JetPURegion = cms.vdouble(0,39,39,100),
  inputWeightFiles = cms.vstring("/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/SingleVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU/TMVAWeight_BDTG_NoPruning_msoftbeta20_PTBin_475_600_PU_0_39/chs_msoftbeta20_PTBin_475_600_PU_0_39_BDTG_NoPruning.weights.xml",
                                "/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVAWeight_BDTG_NoPruning_msoftbeta20_PTBin_475_600_PU_39_100/chs_msoftbeta20_PTBin_475_600_PU_39_100_BDTG_NoPruning.weights.xml")))

process.Options.InputVariableFileParam.insert(False,
 cms.PSet(
  InputVariableList = cms.vstring("mtrimsafe_Rtrim_010_Ptfrac_003[0]"), 
  MethodName = cms.string("BDTG"), 
  JetPTRegion = cms.vdouble(475,600), 
  JetPURegion = cms.vdouble(0,39,39,100),
  inputWeightFiles = cms.vstring("/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/SingleVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU/TMVAWeight_BDTG_NoPruning_mtrim_PTBin_475_600_PU_0_39/chs_mtrim_PTBin_475_600_PU_0_39_BDTG_NoPruning.weights.xml",
                                "/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVAWeight_BDTG_NoPruning_mtrim_PTBin_475_600_PU_39_100/chs_mtrim_PTBin_475_600_PU_39_100_BDTG_NoPruning.weights.xml")))

process.Options.InputVariableFileParam.insert(False,
 cms.PSet(
  InputVariableList = cms.vstring("tau1[0]"), 
  MethodName = cms.string("BDTG"), 
  JetPTRegion = cms.vdouble(475,600), 
  JetPURegion = cms.vdouble(0,39,39,100),
  inputWeightFiles = cms.vstring("/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/SingleVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU/TMVAWeight_BDTG_NoPruning_tau1_PTBin_475_600_PU_0_39/chs_tau1_PTBin_475_600_PU_0_39_BDTG_NoPruning.weights.xml",
                                "/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVAWeight_BDTG_NoPruning_tau1_PTBin_475_600_PU_39_100/chs_tau1_PTBin_475_600_PU_39_100_BDTG_NoPruning.weights.xml")))

process.Options.InputVariableFileParam.insert(False,
 cms.PSet(
  InputVariableList = cms.vstring("tau2[0]"), 
  MethodName = cms.string("BDTG"), 
  JetPTRegion = cms.vdouble(475,600), 
  JetPURegion = cms.vdouble(0,39,39,100),
  inputWeightFiles = cms.vstring("/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/SingleVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU/TMVAWeight_BDTG_NoPruning_tau2_PTBin_475_600_PU_0_39/chs_tau2_PTBin_475_600_PU_0_39_BDTG_NoPruning.weights.xml",
                                "/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVAWeight_BDTG_NoPruning_tau2_PTBin_475_600_PU_39_100/chs_tau2_PTBin_475_600_PU_39_100_BDTG_NoPruning.weights.xml")))

process.Options.InputVariableFileParam.insert(False,
 cms.PSet(
  InputVariableList = cms.vstring("tau2/tau1[0]"), 
  MethodName = cms.string("BDTG"), 
  JetPTRegion = cms.vdouble(475,600), 
  JetPURegion = cms.vdouble(0,39,39,100),
  inputWeightFiles = cms.vstring("/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/SingleVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU/TMVAWeight_BDTG_NoPruning_tau2_tau1_PTBin_475_600_PU_0_39/chs_tau2_tau1_PTBin_475_600_PU_0_39_BDTG_NoPruning.weights.xml",
                                "/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVAWeight_BDTG_NoPruning_tau2_tau1_PTBin_475_600_PU_39_100/chs_tau2_tau1_PTBin_475_600_PU_39_100_BDTG_NoPruning.weights.xml")))
## pair BDT
process.Options.InputVariableFileParam.insert(False,
 cms.PSet(
  InputVariableList = cms.vstring("ecf_beta_10[0]","QGLikelihood_pr_zcut_010_R_cut_050[0]"), 
  MethodName = cms.string("BDTG"), 
  JetPTRegion = cms.vdouble(475,600), 
  JetPURegion = cms.vdouble(0,39,39,100),
  inputWeightFiles = cms.vstring("/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/PairVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU/TMVAWeight_BDTG_NoPruning_ECFbeta10_QGLikelihood_PTBin_475_600_PU_0_39/chs_ECFbeta10_QGLikelihood_PTBin_475_600_PU_0_39_BDTG_NoPruning.weights.xml",
                                "/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/PairVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVAWeight_BDTG_NoPruning_ECFbeta10_QGLikelihood_PTBin_475_600_PU_39_100/chs_ECFbeta10_QGLikelihood_PTBin_475_600_PU_39_100_BDTG_NoPruning.weights.xml")))

process.Options.InputVariableFileParam.insert(False,
 cms.PSet(
  InputVariableList = cms.vstring("ecf_beta_10[0]","2*QGLikelihood_pr_sub2_zcut_010_R_cut_050[0]+QGLikelihood_pr_sub1_zcut_010_R_cut_050[0]"), 
  MethodName = cms.string("BDTG"), 
  JetPTRegion = cms.vdouble(475,600), 
  JetPURegion = cms.vdouble(0,39,39,100),
  inputWeightFiles = cms.vstring("/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/PairVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU/TMVAWeight_BDTG_NoPruning_ECFbeta10_QGLikelihood_comb_PTBin_475_600_PU_0_39/chs_ECFbeta10_QGLikelihood_comb_PTBin_475_600_PU_0_39_BDTG_NoPruning.weights.xml",
                                "/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/PairVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVAWeight_BDTG_NoPruning_ECFbeta10_QGLikelihood_comb_PTBin_475_600_PU_39_100/chs_ECFbeta10_QGLikelihood_comb_PTBin_475_600_PU_39_100_BDTG_NoPruning.weights.xml")))

process.Options.InputVariableFileParam.insert(False,
 cms.PSet(
  InputVariableList = cms.vstring("ecf_beta_10[0]","QGLikelihood_pr_sub1_zcut_010_R_cut_050[0]"), 
  MethodName = cms.string("BDTG"), 
  JetPTRegion = cms.vdouble(475,600), 
  JetPURegion = cms.vdouble(0,39,39,100),
  inputWeightFiles = cms.vstring("/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/PairVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU/TMVAWeight_BDTG_NoPruning_ECFbeta10_QGLikelihood_sub1_PTBin_475_600_PU_0_39/chs_ECFbeta10_QGLikelihood_sub1_PTBin_475_600_PU_0_39_BDTG_NoPruning.weights.xml",
                                "/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/PairVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVAWeight_BDTG_NoPruning_ECFbeta10_QGLikelihood_sub1_PTBin_475_600_PU_39_100/chs_ECFbeta10_QGLikelihood_sub1_PTBin_475_600_PU_39_100_BDTG_NoPruning.weights.xml")))

process.Options.InputVariableFileParam.insert(False,
 cms.PSet(
  InputVariableList = cms.vstring("ecf_beta_10[0]","QGLikelihood_pr_sub2_zcut_010_R_cut_050[0]"), 
  MethodName = cms.string("BDTG"), 
  JetPTRegion = cms.vdouble(475,600), 
  JetPURegion = cms.vdouble(0,39,39,100),
  inputWeightFiles = cms.vstring("/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/PairVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU/TMVAWeight_BDTG_NoPruning_ECFbeta10_QGLikelihood_sub2_PTBin_475_600_PU_0_39/chs_ECFbeta10_QGLikelihood_sub2_PTBin_475_600_PU_0_39_BDTG_NoPruning.weights.xml",
                                "/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/PairVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVAWeight_BDTG_NoPruning_ECFbeta10_QGLikelihood_sub2_PTBin_475_600_PU_39_100/chs_ECFbeta10_QGLikelihood_sub2_PTBin_475_600_PU_39_100_BDTG_NoPruning.weights.xml")))

process.Options.InputVariableFileParam.insert(False,
 cms.PSet(
  InputVariableList = cms.vstring("ecf_beta_10[0]","ecf_beta_15[0]"), 
  MethodName = cms.string("BDTG"), 
  JetPTRegion = cms.vdouble(475,600), 
  JetPURegion = cms.vdouble(0,39,39,100),
  inputWeightFiles = cms.vstring("/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/PairVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU/TMVAWeight_BDTG_NoPruning_ECFbeta15_ECFbeta10_PTBin_475_600_PU_0_39/chs_ECFbeta15_ECFbeta10_PTBin_475_600_PU_0_39_BDTG_NoPruning.weights.xml",
                                "/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/PairVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVAWeight_BDTG_NoPruning_ECFbeta15_ECFbeta10_PTBin_475_600_PU_39_100/chs_ECFbeta15_ECFbeta10_PTBin_475_600_PU_39_100_BDTG_NoPruning.weights.xml")))

process.Options.InputVariableFileParam.insert(False,
 cms.PSet(
  InputVariableList = cms.vstring("ecf_beta_15[0]","QGLikelihood_pr_zcut_010_R_cut_050[0]"), 
  MethodName = cms.string("BDTG"), 
  JetPTRegion = cms.vdouble(475,600), 
  JetPURegion = cms.vdouble(0,39,39,100),
  inputWeightFiles = cms.vstring("/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/PairVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU/TMVAWeight_BDTG_NoPruning_ECFbeta15_QGLikelihood_PTBin_475_600_PU_0_39/chs_ECFbeta15_QGLikelihood_PTBin_475_600_PU_0_39_BDTG_NoPruning.weights.xml",
                                "/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/PairVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVAWeight_BDTG_NoPruning_ECFbeta15_QGLikelihood_PTBin_475_600_PU_39_100/chs_ECFbeta15_QGLikelihood_PTBin_475_600_PU_39_100_BDTG_NoPruning.weights.xml")))

process.Options.InputVariableFileParam.insert(False,
 cms.PSet(
  InputVariableList = cms.vstring("ecf_beta_15[0]","2*QGLikelihood_sub2_pr_zcut_010_R_cut_050[0]+QGLikelihood_sub1_pr_zcut_010_R_cut_050[0]"), 
  MethodName = cms.string("BDTG"), 
  JetPTRegion = cms.vdouble(475,600), 
  JetPURegion = cms.vdouble(0,39,39,100),
  inputWeightFiles = cms.vstring("/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/PairVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU/TMVAWeight_BDTG_NoPruning_ECFbeta15_QGLikelihood_comb_PTBin_475_600_PU_0_39/chs_ECFbeta15_QGLikelihood_comb_PTBin_475_600_PU_0_39_BDTG_NoPruning.weights.xml",
                                "/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/PairVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVAWeight_BDTG_NoPruning_ECFbeta15_QGLikelihood_comb_PTBin_475_600_PU_39_100/chs_ECFbeta15_QGLikelihood_comb_PTBin_475_600_PU_39_100_BDTG_NoPruning.weights.xml")))

process.Options.InputVariableFileParam.insert(False,
 cms.PSet(
  InputVariableList = cms.vstring("ecf_beta_15[0]","QGLikelihood_sub1_pr_zcut_010_R_cut_050[0]"), 
  MethodName = cms.string("BDTG"), 
  JetPTRegion = cms.vdouble(475,600), 
  JetPURegion = cms.vdouble(0,39,39,100),
  inputWeightFiles = cms.vstring("/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/PairVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU/TMVAWeight_BDTG_NoPruning_ECFbeta15_QGLikelihood_sub1_PTBin_475_600_PU_0_39/chs_ECFbeta15_QGLikelihood_sub1_PTBin_475_600_PU_0_39_BDTG_NoPruning.weights.xml",
                                "/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/PairVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVAWeight_BDTG_NoPruning_ECFbeta15_QGLikelihood_sub1_PTBin_475_600_PU_39_100/chs_ECFbeta15_QGLikelihood_sub1_PTBin_475_600_PU_39_100_BDTG_NoPruning.weights.xml")))

process.Options.InputVariableFileParam.insert(False,
 cms.PSet(
  InputVariableList = cms.vstring("ecf_beta_15[0]","QGLikelihood_sub2_pr_zcut_010_R_cut_050[0]"), 
  MethodName = cms.string("BDTG"), 
  JetPTRegion = cms.vdouble(475,600), 
  JetPURegion = cms.vdouble(0,39,39,100),
  inputWeightFiles = cms.vstring("/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/PairVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU/TMVAWeight_BDTG_NoPruning_ECFbeta15_QGLikelihood_sub2_PTBin_475_600_PU_0_39/chs_ECFbeta15_QGLikelihood_sub2_PTBin_475_600_PU_0_39_BDTG_NoPruning.weights.xml",
                                "/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/PairVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVAWeight_BDTG_NoPruning_ECFbeta15_QGLikelihood_sub2_PTBin_475_600_PU_39_100/chs_ECFbeta15_QGLikelihood_sub2_PTBin_475_600_PU_39_100_BDTG_NoPruning.weights.xml")))

process.Options.InputVariableFileParam.insert(False,
 cms.PSet(
  InputVariableList = cms.vstring("QGLikelihood_pr_zcut_010_R_cut_050[0]","QGLikelihood_sub1_pr_zcut_010_R_cut_050[0]"), 
  MethodName = cms.string("BDTG"), 
  JetPTRegion = cms.vdouble(475,600), 
  JetPURegion = cms.vdouble(0,39,39,100),
  inputWeightFiles = cms.vstring("/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/PairVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU/TMVAWeight_BDTG_NoPruning_QGLikelihood_QGLikelihood_sub1_PTBin_475_600_PU_0_39/chs_QGLikelihood_QGLikelihood_sub1_PTBin_475_600_PU_0_39_BDTG_NoPruning.weights.xml",
                                "/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/PairVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVAWeight_BDTG_NoPruning_QGLikelihood_QGLikelihood_sub1_PTBin_475_600_PU_39_100/chs_QGLikelihood_QGLikelihood_sub1_PTBin_475_600_PU_39_100_BDTG_NoPruning.weights.xml")))

process.Options.InputVariableFileParam.insert(False,
 cms.PSet(
  InputVariableList = cms.vstring("QGLikelihood_pr_zcut_010_R_cut_050[0]","QGLikelihood_sub2_pr_zcut_010_R_cut_050[0]"), 
  MethodName = cms.string("BDTG"), 
  JetPTRegion = cms.vdouble(475,600), 
  JetPURegion = cms.vdouble(0,39,39,100),
  inputWeightFiles = cms.vstring("/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/PairVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU/TMVAWeight_BDTG_NoPruning_QGLikelihood_QGLikelihood_sub2_PTBin_475_600_PU_0_39/chs_QGLikelihood_QGLikelihood_sub2_PTBin_475_600_PU_0_39_BDTG_NoPruning.weights.xml",
                                "/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/PairVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVAWeight_BDTG_NoPruning_QGLikelihood_QGLikelihood_sub2_PTBin_475_600_PU_39_100/chs_QGLikelihood_QGLikelihood_sub2_PTBin_475_600_PU_39_100_BDTG_NoPruning.weights.xml")))

process.Options.InputVariableFileParam.insert(False,
 cms.PSet(
  InputVariableList = cms.vstring("QGLikelihood_pr_zcut_010_R_cut_050[0]","2*QGLikelihood_sub2_pr_zcut_010_R_cut_050[0]+QGLikelihood_sub1_pr_zcut_010_R_cut_050[0]"), 
  MethodName = cms.string("BDTG"), 
  JetPTRegion = cms.vdouble(475,600), 
  JetPURegion = cms.vdouble(0,39,39,100),
  inputWeightFiles = cms.vstring("/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/PairVariablesTraining_highPT_lowPU_BDTG/outputTMVATraining_highPT_lowPU/TMVAWeight_BDTG_NoPruning_QGLikelihood_QGLikelihood_comb_PTBin_475_600_PU_0_39/chs_QGLikelihood_QGLikelihood_comb_PTBin_475_600_PU_0_39_BDTG_NoPruning.weights.xml",
                                "/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/PairVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVAWeight_BDTG_NoPruning_QGLikelihood_QGLikelihood_comb_PTBin_475_600_PU_39_100/chs_QGLikelihood_QGLikelihood_comb_PTBin_475_600_PU_39_100_BDTG_NoPruning.weights.xml")))


'''
TMVAWeight_BDTG_NoPruning_QGLikelihood_sub1_QGLikelihood_comb_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_QGLikelihood_sub1_QGLikelihood_sub2_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_QGLikelihood_sub2_QGLikelihood_comb_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_Qjets_ECFbeta10_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_Qjets_ECFbeta15_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_Qjets_QGLikelihood_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_Qjets_QGLikelihood_comb_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_Qjets_QGLikelihood_sub1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_Qjets_QGLikelihood_sub2_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mconst_ECFbeta10_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mconst_ECFbeta15_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mconst_QGLikelihood_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mconst_QGLikelihood_comb_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mconst_QGLikelihood_sub1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mconst_QGLikelihood_sub2_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mconst_Qjets_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mconst_msoftbeta00_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mconst_msoftbeta10_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mconst_mtrim_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mconst_tau1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mconst_tau2_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mconst_tau2tau1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_ECFbeta10_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_ECFbeta15_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_QGLikelihood_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_QGLikelihood_comb_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_QGLikelihood_sub1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_QGLikelihood_sub2_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_Qjets_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_mconst_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_msoftbeta00_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_msoftbeta10_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_mtrim_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_tau1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_tau2_PTBin_475_600_PU_0_39

TMVAWeight_BDTG_NoPruning_msoftbeta00_ECFbeta10_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_msoftbeta00_ECFbeta15_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_msoftbeta00_QGLikelihood_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_msoftbeta00_QGLikelihood_comb_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_msoftbeta00_QGLikelihood_sub1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_msoftbeta00_QGLikelihood_sub2_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_msoftbeta00_Qjets_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_msoftbeta00_tau1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_msoftbeta00_tau2_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_msoftbeta00_tau2tau1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_msoftbeta10_ECFbeta10_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_msoftbeta10_ECFbeta15_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_msoftbeta10_QGLikelihood_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_msoftbeta10_QGLikelihood_comb_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_msoftbeta10_QGLikelihood_sub1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_msoftbeta10_QGLikelihood_sub2_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_msoftbeta10_Qjets_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_msoftbeta10_msoftbeta00_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_msoftbeta10_tau1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_msoftbeta10_tau2_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_msoftbeta10_tau2tau1_PTBin_475_600_PU_0_39

TMVAWeight_BDTG_NoPruning_mtrim_ECFbeta10_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mtrim_ECFbeta15_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mtrim_QGLikelihood_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mtrim_QGLikelihood_comb_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mtrim_QGLikelihood_sub1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mtrim_QGLikelihood_sub2_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mtrim_Qjets_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mtrim_msoftbeta00_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mtrim_msoftbeta10_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mtrim_tau1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mtrim_tau2_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mtrim_tau2tau1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_tau1_ECFbeta10_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_tau1_ECFbeta15_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_tau1_QGLikelihood_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_tau1_QGLikelihood_comb_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_tau1_QGLikelihood_sub1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_tau1_QGLikelihood_sub2_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_tau1_Qjets_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_tau1_tau2_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_tau1_tau2tau1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_tau2_ECFbeta10_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_tau2_ECFbeta15_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_tau2_QGLikelihood_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_tau2_QGLikelihood_comb_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_tau2_QGLikelihood_sub1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_tau2_QGLikelihood_sub2_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_tau2_Qjets_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_tau2_tau2tau1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_tau2tau1_ECFbeta10_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_tau2tau1_ECFbeta15_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_tau2tau1_QGLikelihood_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_tau2tau1_QGLikelihood_comb_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_tau2tau1_QGLikelihood_sub1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_tau2tau1_QGLikelihood_sub2_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_tau2tau1_Qjets_PTBin_475_600_PU_0_39
'''

## Triplet BDT

'''
TMVAWeight_BDTG_NoPruning_mpruned_ECFbeta10_QGLikelihood_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_ECFbeta10_QGLikelihood_comb_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_ECFbeta10_QGLikelihood_sub1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_ECFbeta10_QGLikelihood_sub2_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_ECFbeta10_pullangle_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_ECFbeta15_ECFbeta10_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_ECFbeta15_QGLikelihood_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_ECFbeta15_QGLikelihood_comb_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_ECFbeta15_QGLikelihood_sub1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_ECFbeta15_QGLikelihood_sub2_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_ECFbeta15_pullangle_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_QGLikelihood_QGLikelihood_comb_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_QGLikelihood_QGLikelihood_sub1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_QGLikelihood_QGLikelihood_sub2_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_QGLikelihood_comb_pullangle_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_QGLikelihood_pullangle_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_QGLikelihood_sub1_QGLikelihood_comb_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_QGLikelihood_sub1_QGLikelihood_sub2_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_QGLikelihood_sub1_pullangle_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_QGLikelihood_sub2_QGLikelihood_comb_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_QGLikelihood_sub2_pullangle_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_Qjets_ECFbeta10_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_Qjets_ECFbeta15_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_Qjets_QGLikelihood_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_Qjets_QGLikelihood_comb_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_Qjets_QGLikelihood_sub1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_Qjets_QGLikelihood_sub2_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_Qjets_pullangle_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_mconst_ECFbeta10_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_mconst_ECFbeta15_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_mconst_QGLikelihood_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_mconst_QGLikelihood_comb_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_mconst_QGLikelihood_sub1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_mconst_QGLikelihood_sub2_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_mconst_Qjets_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_mconst_msoftbeta00_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_mconst_msoftbeta10_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_mconst_mtrim_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_mconst_pullangle_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_mconst_tau1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_mconst_tau2_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_mconst_tau2tau1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_msoftbeta00_ECFbeta10_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_msoftbeta00_ECFbeta15_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_msoftbeta00_QGLikelihood_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_msoftbeta00_QGLikelihood_comb_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_msoftbeta00_QGLikelihood_sub1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_msoftbeta00_QGLikelihood_sub2_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_msoftbeta00_Qjets_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_msoftbeta00_pullangle_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_msoftbeta00_tau1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_msoftbeta00_tau2_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_msoftbeta00_tau2tau1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_msoftbeta10_ECFbeta10_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_msoftbeta10_ECFbeta15_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_msoftbeta10_QGLikelihood_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_msoftbeta10_QGLikelihood_comb_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_msoftbeta10_QGLikelihood_sub1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_msoftbeta10_QGLikelihood_sub2_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_msoftbeta10_Qjets_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_msoftbeta10_msoftbeta00_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_msoftbeta10_pullangle_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_msoftbeta10_tau1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_msoftbeta10_tau2_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_msoftbeta10_tau2tau1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_mtrim_ECFbeta10_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_mtrim_ECFbeta15_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_mtrim_QGLikelihood_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_mtrim_QGLikelihood_comb_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_mtrim_QGLikelihood_sub1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_mtrim_QGLikelihood_sub2_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_mtrim_Qjets_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_mtrim_msoftbeta00_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_mtrim_msoftbeta10_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_mtrim_pullangle_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_mtrim_tau1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_mtrim_tau2_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_mtrim_tau2tau1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_tau1_ECFbeta10_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_tau1_ECFbeta15_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_tau1_QGLikelihood_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_tau1_QGLikelihood_comb_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_tau1_QGLikelihood_sub1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_tau1_QGLikelihood_sub2_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_tau1_Qjets_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_tau1_pullangle_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_tau1_tau2_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_tau1_tau2tau1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_tau2_ECFbeta10_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_tau2_ECFbeta15_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_tau2_QGLikelihood_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_tau2_QGLikelihood_comb_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_tau2_QGLikelihood_sub1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_tau2_QGLikelihood_sub2_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_tau2_Qjets_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_tau2_pullangle_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_tau2_tau2tau1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_tau2tau1_ECFbeta10_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_tau2tau1_ECFbeta15_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_tau2tau1_QGLikelihood_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_tau2tau1_QGLikelihood_comb_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_tau2tau1_QGLikelihood_sub1_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_tau2tau1_QGLikelihood_sub2_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_tau2tau1_Qjets_PTBin_475_600_PU_0_39
TMVAWeight_BDTG_NoPruning_mpruned_tau2tau1_pullangle_PTBin_475_600_PU_0_39
'''
