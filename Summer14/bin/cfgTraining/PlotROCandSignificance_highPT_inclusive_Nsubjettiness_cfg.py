import FWCore.ParameterSet.Config as cms

process = cms.Process("PlotROCansSignificance")

process.Options = cms.PSet(

   ## Name of the tree on which perform th training                                                                                                                                
   TreeName           = cms.string("chs"),

   ## Name for ggH signal --> match what is written in the sample list file                                                                                                             
   SignalggHName      = cms.string("RSGWW1000"),

   ## Name for ggH signal --> match what is written in the sample list file                                                                                                             
   SignalqqHName      = cms.string(""),

   ## string which use tree branches to re-weight the events                                                                                                                              
   EventWeight        = cms.string(""),

   ## 0 use both ggH and qqH as signal, 1 use only ggH as signal, 2 use only qqH as signal                                                                                                
   useTypeOfSignal    = cms.uint32(0),

   ## string which is used in the TMVATraining class to define a cut to be applied on the events                                                                                         
   PreselectionCutType = cms.string("basicJetsCutCSA14"),

   ## luminosity in order to  compute signal and bk expectations
   Lumi   = cms.double(19297),

   ## Lepton Type: Muon, Electron ,EleMu and Jets (fully hadronic)                                                                                                                        
   LeptonType         = cms.string("Jets"),

   ## output directory for root and weight file                                                                                                                                           
   outputPlotDirectory  = cms.string("TMVATrainingPlots_CHS_highPT_inclusive/"),

   ##input file and variables   
   InputInformationParam =  cms.VPSet(
    cms.PSet( fileName = cms.string("SingleVariablesTraining_highPT_inclusive_BDTG/outputTMVATraining_highPT/TMVATrainingResult_tau2_tau1_beta_05_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("#tau_{2}/#tau_{1} #beta=0.5"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100), ReducedName = cms.string("tau2tau1beta05")),
    cms.PSet( fileName = cms.string("SingleVariablesTraining_highPT_inclusive_BDTG/outputTMVATraining_highPT/TMVATrainingResult_tau2_tau1_beta_08_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("#tau_{2}/#tau_{1} #beta=0.8"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100), ReducedName = cms.string("tau2tau1beta08")),
    cms.PSet( fileName = cms.string("SingleVariablesTraining_highPT_inclusive_BDTG/outputTMVATraining_highPT/TMVATrainingResult_tau2_tau1_beta_10_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("#tau_{2}/#tau_{1} #beta=1.0"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100), ReducedName = cms.string("tau2tau1beta10")),
    cms.PSet( fileName = cms.string("SingleVariablesTraining_highPT_inclusive_BDTG/outputTMVATraining_highPT/TMVATrainingResult_tau2_tau1_beta_12_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("#tau_{2}/#tau_{1} #beta=1.2"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100), ReducedName = cms.string("tau2tau1beta12")),
    cms.PSet( fileName = cms.string("SingleVariablesTraining_highPT_inclusive_BDTG/outputTMVATraining_highPT/TMVATrainingResult_tau2_tau1_beta_15_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("#tau_{2}/#tau_{1} #beta=1.5"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100), ReducedName = cms.string("tau2tau1beta15")),
    cms.PSet( fileName = cms.string("SingleVariablesTraining_highPT_inclusive_BDTG/outputTMVATraining_highPT/TMVATrainingResult_tau2_tau1_beta_20_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("#tau_{2}/#tau_{1} #beta=2.0"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100), ReducedName = cms.string("tau2tau1beta20")),
    cms.PSet( fileName = cms.string("SingleVariablesTraining_highPT_inclusive_BDTG/outputTMVATraining_highPT/TMVATrainingResult_tau2_tau1_beta_25_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("#tau_{2}/#tau_{1} #beta=2.5"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100), ReducedName = cms.string("tau2tau1beta25")),
    cms.PSet( fileName = cms.string("SingleVariablesTraining_highPT_inclusive_BDTG/outputTMVATraining_highPT/TMVATrainingResult_tau2_tau1_beta_30_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("#tau_{2}/#tau_{1} #beta=3.0"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100), ReducedName = cms.string("tau2tau1beta30")),
#    cms.PSet( fileName = cms.string("AllVariablesTraining_highPT_inclusive_BDTG/outputTMVATraining_highPT/TMVATrainingResult_allvariables_PTBin_475_600_PU_39_100.root"),  inputVariableOrMethodName = cms.vstring("all"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(39,100), ReducedName = cms.string("allvariables")),
  )
)
