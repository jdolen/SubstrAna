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
   outputPlotDirectory  = cms.string("TMVATrainingPlots_CHS_highPT_highPU/"),

   ##input file and variables   
   InputInformationParam =  cms.VPSet(
    cms.PSet( fileName = cms.string("SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVATrainingResult_mraw_PTBin_475_600_PU_39_100.root"),  inputVariableOrMethodName = cms.vstring("M_{raw}"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(39,100), ReducedName = cms.string("mraw")),
    cms.PSet( fileName = cms.string("SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVATrainingResult_mclean_PTBin_475_600_PU_39_100.root"),  inputVariableOrMethodName = cms.vstring("M_{clean}"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(39,100), ReducedName = cms.string("mclean")),
    cms.PSet( fileName = cms.string("SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVATrainingResult_mconst_PTBin_475_600_PU_39_100.root"),  inputVariableOrMethodName = cms.vstring("M_{const}"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(39,100), ReducedName = cms.string("mconst")),
    cms.PSet( fileName = cms.string("SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVATrainingResult_mpruned_PTBin_475_600_PU_39_100.root"),  inputVariableOrMethodName = cms.vstring("M_{pruned}"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(39,100), ReducedName = cms.string("mpruned")),
    cms.PSet( fileName = cms.string("SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVATrainingResult_msoftbeta00_PTBin_475_600_PU_39_100.root"),  inputVariableOrMethodName = cms.vstring("M_{soft}(#beta=0)"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(39,100), ReducedName = cms.string("msoftbeta00")),
    cms.PSet( fileName = cms.string("SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVATrainingResult_msoftbeta10_PTBin_475_600_PU_39_100.root"),  inputVariableOrMethodName = cms.vstring("M_{soft}(#beta=1)"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(39,100), ReducedName = cms.string("msoftbeta10")),
    cms.PSet( fileName = cms.string("SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVATrainingResult_msoftbeta20_PTBin_475_600_PU_39_100.root"),  inputVariableOrMethodName = cms.vstring("M_{soft}(#beta=2)"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(39,100), ReducedName = cms.string("msoftbeta20")),
    cms.PSet( fileName = cms.string("SingleVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVATrainingResult_mtrim_PTBin_475_600_PU_39_100.root"),  inputVariableOrMethodName = cms.vstring("M_{trimmed}"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(39,100), ReducedName = cms.string("mtrim")),
#    cms.PSet( fileName = cms.string("AllVariablesTraining_highPT_highPU_BDTG/outputTMVATraining_highPT_highPU/TMVATrainingResult_allvariables_PTBin_475_600_PU_39_100.root"),  inputVariableOrMethodName = cms.vstring("all"), JetPtBinOfTraining = cms.vdouble(475,600), PileUpBinOfTraining = cms.vdouble(39,100), ReducedName = cms.string("allvariables")),
  )
)
