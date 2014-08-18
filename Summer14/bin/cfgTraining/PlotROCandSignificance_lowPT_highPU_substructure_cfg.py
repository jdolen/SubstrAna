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
   outputPlotDirectory  = cms.string("TMVATrainingPlots_CHS_lowPT_highPU/"),

   ##input file and variables   
   InputInformationParam =  cms.VPSet(
    cms.PSet( fileName = cms.string("SingleVariablesTraining_lowPT_highPU_BDTG/outputTMVATraining_lowPT_highPU/TMVATrainingResult_ECFbeta10_PTBin_300_450_PU_39_100.root"),  inputVariableOrMethodName = cms.vstring("C2(#beta=1)"), JetPtBinOfTraining = cms.vdouble(300,450), PileUpBinOfTraining = cms.vdouble(40,100), ReducedName = cms.string("ECFbeta10")),
    cms.PSet( fileName = cms.string("SingleVariablesTraining_lowPT_highPU_BDTG/outputTMVATraining_lowPT_highPU/TMVATrainingResult_ECFbeta15_PTBin_300_450_PU_39_100.root"),  inputVariableOrMethodName = cms.vstring("C2(#beta=1.5)"), JetPtBinOfTraining = cms.vdouble(300,450), PileUpBinOfTraining = cms.vdouble(40,100), ReducedName = cms.string("ECFbeta15")),
    cms.PSet( fileName = cms.string("SingleVariablesTraining_lowPT_highPU_BDTG/outputTMVATraining_lowPT_highPU/TMVATrainingResult_ECFbeta20_PTBin_300_450_PU_39_100.root"),  inputVariableOrMethodName = cms.vstring("C2(#beta=2)"), JetPtBinOfTraining = cms.vdouble(300,450), PileUpBinOfTraining = cms.vdouble(40,100), ReducedName = cms.string("ECFbeta20")),
    cms.PSet( fileName = cms.string("SingleVariablesTraining_lowPT_highPU_BDTG/outputTMVATraining_lowPT_highPU/TMVATrainingResult_QGLikelihood_PTBin_300_450_PU_39_100.root"),  inputVariableOrMethodName = cms.vstring("QG Likelihood"), JetPtBinOfTraining = cms.vdouble(300,450), PileUpBinOfTraining = cms.vdouble(40,100), ReducedName = cms.string("QGLikelihood")),
    cms.PSet( fileName = cms.string("SingleVariablesTraining_lowPT_highPU_BDTG/outputTMVATraining_lowPT_highPU/TMVATrainingResult_QGLikelihood_sub2_PTBin_300_450_PU_39_100.root"),  inputVariableOrMethodName = cms.vstring("QG Likelihood sub2"), JetPtBinOfTraining = cms.vdouble(300,450), PileUpBinOfTraining = cms.vdouble(40,100), ReducedName = cms.string("QGLikelihood_sub2")),
    cms.PSet( fileName = cms.string("SingleVariablesTraining_lowPT_highPU_BDTG/outputTMVATraining_lowPT_highPU/TMVATrainingResult_Qjets_PTBin_300_450_PU_39_100.root"),  inputVariableOrMethodName = cms.vstring("#Gamma_{Qjet}"), JetPtBinOfTraining = cms.vdouble(300,450), PileUpBinOfTraining = cms.vdouble(40,100), ReducedName = cms.string("qjet")),
    cms.PSet( fileName = cms.string("SingleVariablesTraining_lowPT_highPU_BDTG/outputTMVATraining_lowPT_highPU/TMVATrainingResult_tau1_PTBin_300_450_PU_39_100.root"),  inputVariableOrMethodName = cms.vstring("#tau_{1}"), JetPtBinOfTraining = cms.vdouble(300,450), PileUpBinOfTraining = cms.vdouble(40,100), ReducedName = cms.string("tau1")),
    cms.PSet( fileName = cms.string("SingleVariablesTraining_lowPT_highPU_BDTG/outputTMVATraining_lowPT_highPU/TMVATrainingResult_tau2_PTBin_300_450_PU_39_100.root"),  inputVariableOrMethodName = cms.vstring("#tau_{2}"), JetPtBinOfTraining = cms.vdouble(300,450), PileUpBinOfTraining = cms.vdouble(40,100), ReducedName = cms.string("tau2")),
    cms.PSet( fileName = cms.string("SingleVariablesTraining_lowPT_highPU_BDTG/outputTMVATraining_lowPT_highPU/TMVATrainingResult_tau2_tau1_PTBin_300_450_PU_39_100.root"),  inputVariableOrMethodName = cms.vstring("#tau_{2}/#tau_{1}"), JetPtBinOfTraining = cms.vdouble(300,450), PileUpBinOfTraining = cms.vdouble(40,100), ReducedName = cms.string("tau2tau1")),
  )
)
