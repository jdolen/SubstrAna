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
    cms.PSet( fileName = cms.string("SingleVariablesTraining_highPT_inclusive_BDTG/outputTMVATraining_highPT/TMVATrainingResult_mpruned_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("M_{pruned}"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100), ReducedName = cms.string("mpruned")),
    cms.PSet( fileName = cms.string("SingleVariablesTraining_highPT_inclusive_BDTG/outputTMVATraining_highPT/TMVATrainingResult_msoftz005beta00_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("M_{soft}(z=0.05,#beta=0)"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100), ReducedName = cms.string("msoftz005beta00")),
#    cms.PSet( fileName = cms.string("SingleVariablesTraining_highPT_inclusive_BDTG/outputTMVATraining_highPT/TMVATrainingResult_msoftz005beta10_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("M_{soft}(z=0.05,#beta=1)"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100), ReducedName = cms.string("msoftz005beta10")),
#    cms.PSet( fileName = cms.string("SingleVariablesTraining_highPT_inclusive_BDTG/outputTMVATraining_highPT/TMVATrainingResult_msoftz005beta20_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("M_{soft}(z=0.05,#beta=2)"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100), ReducedName = cms.string("msoftz005beta20")),
#    cms.PSet( fileName = cms.string("SingleVariablesTraining_highPT_inclusive_BDTG/outputTMVATraining_highPT/TMVATrainingResult_msoftz005betam10_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("M_{soft}(z=0.05,#beta=-1)"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100), ReducedName = cms.string("msoftz005betam10")),
#    cms.PSet( fileName = cms.string("SingleVariablesTraining_highPT_inclusive_BDTG/outputTMVATraining_highPT/TMVATrainingResult_msoftz015beta00_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("M_{soft}(z=0.15,#beta=0)"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100), ReducedName = cms.string("msoftz015beta00")),
    cms.PSet( fileName = cms.string("SingleVariablesTraining_highPT_inclusive_BDTG/outputTMVATraining_highPT/TMVATrainingResult_msoftz015beta10_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("M_{soft}(z=0.15,#beta=1)"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100), ReducedName = cms.string("msoftz015beta10")),
#    cms.PSet( fileName = cms.string("SingleVariablesTraining_highPT_inclusive_BDTG/outputTMVATraining_highPT/TMVATrainingResult_msoftz015beta20_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("M_{soft}(z=0.15,#beta=2)"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100), ReducedName = cms.string("msoftz015beta20")),
#    cms.PSet( fileName = cms.string("SingleVariablesTraining_highPT_inclusive_BDTG/outputTMVATraining_highPT/TMVATrainingResult_msoftz015betam10_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("M_{soft}(z=0.15,#beta=-1)"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100), ReducedName = cms.string("msoftz015betam10")),
    cms.PSet( fileName = cms.string("SingleVariablesTraining_highPT_inclusive_BDTG/outputTMVATraining_highPT/TMVATrainingResult_msoftz01beta00_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("M_{soft}(z=0.1,#beta=0)"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100), ReducedName = cms.string("msoftz01beta00")),
    cms.PSet( fileName = cms.string("SingleVariablesTraining_highPT_inclusive_BDTG/outputTMVATraining_highPT/TMVATrainingResult_msoftz01beta10_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("M_{soft}(z=0.1,#beta=1)"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100), ReducedName = cms.string("msoftz01beta10")),
#    cms.PSet( fileName = cms.string("SingleVariablesTraining_highPT_inclusive_BDTG/outputTMVATraining_highPT/TMVATrainingResult_msoftz01beta20_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("M_{soft}(z=0.1,#beta=2)"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100), ReducedName = cms.string("msoftz01beta20")),
#    cms.PSet( fileName = cms.string("SingleVariablesTraining_highPT_inclusive_BDTG/outputTMVATraining_highPT/TMVATrainingResult_msoftz01betam10_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("M_{soft}(z=0.1,#beta=-1)"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(0,100), ReducedName = cms.string("msoftz01betam10")),
#    cms.PSet( fileName = cms.string("AllVariablesTraining_highPT_inclusive_BDTG/outputTMVATraining_highPT/TMVATrainingResult_allvariables_PTBin_475_600_PU_0_100.root"),  inputVariableOrMethodName = cms.vstring("all"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(39,100), ReducedName = cms.string("allvariables")),
  )
)
