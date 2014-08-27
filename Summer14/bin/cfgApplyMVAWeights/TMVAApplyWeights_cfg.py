import FWCore.ParameterSet.Config as cms

process = cms.Process("")

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
   LeptonType = cms.string("Jets")

   ## string which is used in the TMVATraining class to define a cut to be applied on the events
   PreselectionCutType = cms.string("basicJetsCutCSA14"),

   ## W-jet pt bin for training
   JetPtBinOfTraining = cms.vdouble(475,600),

   ## In time pile-up bin for the training
   PileUpBinOfTraining = cms.vdouble(0,39),

   ## specify everything we need to take the right weight file
   InputWeightFileParam = cms.VPSet(),

   ## which jet read to fill the output : 0 means leading one, 1 ... etc
   JetToRead = cms.uint32(0),

   ## apply or not training cut when the MVA reading is used .. if yes, events out of the training phase space are put to a default value, otherwise we use the value
   ## given by the TMVAReader
   optionCut = cms.bool(False),

)

# single BDT, hight pT low PU


process.Options.InputWeightFileParam.insert(False,cms.PSet(inputWeightFile = cms.string("/afs/cern.ch/user/r/rgerosa/work/JMENtuplesFramework/CMSSW_6_2_8/src/SubstrAna/Summer14/bin/cfgTraining/SingleVariablesTraining_highPT_lowPU/outputTMVATraining_highPT_lowPU/TMVAWeight_CutsSA_ECFbeta10_PTBin_%d_%d_PU_%d_%d_/chs_ECFbeta10_PTBin_ecf_beta_%d_%d_PU_%d_%d_CutsSA.weights.xml"), MethodName = cms.string("BDTG"), InputVariableList = cms.vstring("ecf_beta_10[0]")))
 
