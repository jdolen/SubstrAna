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

# single BDT, hight pT low PU

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
 
