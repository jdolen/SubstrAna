import FWCore.ParameterSet.Config as cms

process = cms.Process("MiniNtuplizer")

process.Options = cms.PSet(
    
    isMC            = cms.bool(True),
    maxEvents       = cms.int32(MAXEVENTS),    # maximum events to  run
    minEvents       = cms.int32(MINEVENTS),     # maximum events to  run

    jetR            = cms.double(0.8),  # basic clustering cone size
    jetPtCut        = cms.double(25.0), # pt cut on pf and Gen jets  
    jetAlgo         = cms.string('antikt_algorithm'), # ex: antikt_algorithm, ak, AK, cambridge_algorithm, ca, CA
    doCMSSWJets     = cms.bool(False),  # analyze also basic cmssw reconstructed jets 
    doSoftKillerJets = cms.bool(False),  # produce soft killer jets                                                                                                                        

    puppiConfig     = cms.string("PUPPICONFIG"), # puppi configuration to run

    jetPtTresholdForGroomers      = cms.double(100.),
    genjetPtTresholdForGroomers   = cms.double(50.),
    jetPtTresholdForTopTagging    = cms.double(300.),
    genJetPtTresholdForTopTagging = cms.double(250.),

    L1FastJetJEC    = cms.string("/afs/cern.ch/user/b/bmahakud/public/JEC/POSTLS170_V6_L1FastJet_AK7PF.txt"),
    L2RelativeJEC   = cms.string("/afs/cern.ch/user/b/bmahakud/public/JEC/POSTLS170_V6_L2Relative_AK7PF.txt"),
    L3AbsoluteJEC   = cms.string("/afs/cern.ch/user/b/bmahakud/public/JEC/POSTLS170_V6_L3Absolute_AK7PF.txt"),
    JECUncertainty  = cms.string("/afs/cern.ch/user/b/bmahakud/public/JEC/POSTLS170_V6_Uncertainty_AK7PF.txt"),
    L2L3ResidualJEC = cms.string(""), 

    L1FastJetJEC_CHS    = cms.string("/afs/cern.ch/user/b/bmahakud/public/JEC/POSTLS170_V6_L1FastJet_AK7PFchs.txt"),
    L2RelativeJEC_CHS   = cms.string("/afs/cern.ch/user/b/bmahakud/public/JEC/POSTLS170_V6_L2Relative_AK7PFchs.txt"),
    L3AbsoluteJEC_CHS   = cms.string("/afs/cern.ch/user/b/bmahakud/public/JEC/POSTLS170_V6_L3Absolute_AK7PFchs.txt"),
    JECUncertainty_CHS  = cms.string("/afs/cern.ch/user/b/bmahakud/public/JEC/POSTLS170_V6_Uncertainty_AK7PFchs.txt"),
    L2L3ResidualJEC_CHS = cms.string(""), 
    PhilJEC             = cms.string("/afs/cern.ch/user/p/pharris/public/JEC/PhilJEC.root"),

    QGinputWeightFilePath  = cms.string("SubstrAna/Summer14/data/"),

    computeJetFlavour      = cms.bool(True),

    #mc matching
    DoMatchingToBoson = cms.bool(True), # use this to run on WW, ttbar...
    pdgIdBoson        = cms.int32(24),
    dRMatiching       = cms.double(0.3),
    dRLeptonCleaning  = cms.double(0.8),
  
    #softkiller
    softKiller = cms.PSet( ymax = cms.double(2.5), cell_size = cms.double(0.4) ),

    #softdrop                                                                                                                                                                             
    softDrop = cms.VPSet(
     cms.PSet( beta = cms.double(2.),  symmetry_cut = cms.double(0.1),  R0 = cms.double(1.)),
     cms.PSet( beta = cms.double(0.),  symmetry_cut = cms.double(0.1),  R0 = cms.double(1.)),
     cms.PSet( beta = cms.double(1.),  symmetry_cut = cms.double(0.1),  R0 = cms.double(1.)),
     cms.PSet( beta = cms.double(-1.), symmetry_cut = cms.double(0.1),  R0 = cms.double(1.)),

     cms.PSet( beta = cms.double(2.),  symmetry_cut = cms.double(0.05),  R0 = cms.double(1.)),
     cms.PSet( beta = cms.double(0.),  symmetry_cut = cms.double(0.05),  R0 = cms.double(1.)),
     cms.PSet( beta = cms.double(1.),  symmetry_cut = cms.double(0.05),  R0 = cms.double(1.)),
     cms.PSet( beta = cms.double(-1.), symmetry_cut = cms.double(0.05),  R0 = cms.double(1.)),

     cms.PSet( beta = cms.double(2.),  symmetry_cut = cms.double(0.15),  R0 = cms.double(1.)),
     cms.PSet( beta = cms.double(0.),  symmetry_cut = cms.double(0.15),  R0 = cms.double(1.)),
     cms.PSet( beta = cms.double(1.),  symmetry_cut = cms.double(0.15),  R0 = cms.double(1.)),
     cms.PSet( beta = cms.double(-1.), symmetry_cut = cms.double(0.15),  R0 = cms.double(1.)),
    ),


    # trimming
    trimming = cms.VPSet( 
     #cms.PSet( R_trimming = cms.double(0.2), PtFraction = cms.double(0.05), trimAlgo = cms.string('kt_algorithm')),
     cms.PSet( R_trimming = cms.double(0.1), PtFraction = cms.double(0.03), trimAlgo = cms.string('kt_algorithm')),
     #cms.PSet( R_trimming = cms.double(0.2), PtFraction = cms.double(0.03), trimAlgo = cms.string('kt_algorithm')),
     #cms.PSet( R_trimming = cms.double(0.3), PtFraction = cms.double(0.03), trimAlgo = cms.string('kt_algorithm'))
   ),

    #pruning
    pruning =  cms.VPSet(
     cms.PSet( z_cut = cms.double(0.1),  R_Cut = cms.double(0.5),  R_jet_def_pruning = cms.double(0.8), pruneAlgo = cms.string('cambridge_algorithm')),
     #cms.PSet( z_cut = cms.double(0.05), R_Cut = cms.double(0.5),  R_jet_def_pruning = cms.double(0.8), pruneAlgo = cms.string('cambridge_algorithm')),
     #cms.PSet( z_cut = cms.double(0.05), R_Cut = cms.double(0.75), R_jet_def_pruning = cms.double(0.8), pruneAlgo = cms.string('cambridge_algorithm')),
     #cms.PSet( z_cut = cms.double(0.1),  R_Cut = cms.double(0.75), R_jet_def_pruning = cms.double(0.8), pruneAlgo = cms.string('cambridge_algorithm'))
    ),

    #charge
    jetcharge = cms.vdouble(0.5,0.7,1.0),

    #energy correlator                                                                                                                                                            
    energyCorrelator = cms.VPSet(

     cms.PSet( ecfAlgo = cms.string('antikt_algorithm'), Rparam = cms.double(2.0), nPoint = cms.int32(1), beta = cms.double(1.0), type = cms.int32(0)),

     cms.PSet( ecfAlgo = cms.string('antikt_algorithm'), Rparam = cms.double(2.0), nPoint = cms.int32(2), beta = cms.double(0.5), type = cms.int32(0)),
     cms.PSet( ecfAlgo = cms.string('antikt_algorithm'), Rparam = cms.double(2.0), nPoint = cms.int32(2), beta = cms.double(1.0), type = cms.int32(0)),
     cms.PSet( ecfAlgo = cms.string('antikt_algorithm'), Rparam = cms.double(2.0), nPoint = cms.int32(2), beta = cms.double(1.5), type = cms.int32(0)),
     cms.PSet( ecfAlgo = cms.string('antikt_algorithm'), Rparam = cms.double(2.0), nPoint = cms.int32(2), beta = cms.double(2.0), type = cms.int32(0)),
     cms.PSet( ecfAlgo = cms.string('antikt_algorithm'), Rparam = cms.double(2.0), nPoint = cms.int32(2), beta = cms.double(2.5), type = cms.int32(0)),
     cms.PSet( ecfAlgo = cms.string('antikt_algorithm'), Rparam = cms.double(2.0), nPoint = cms.int32(2), beta = cms.double(3.0), type = cms.int32(0)),

     cms.PSet( ecfAlgo = cms.string('antikt_algorithm'), Rparam = cms.double(2.0), nPoint = cms.int32(3), beta = cms.double(0.5), type = cms.int32(0)),
     cms.PSet( ecfAlgo = cms.string('antikt_algorithm'), Rparam = cms.double(2.0), nPoint = cms.int32(3), beta = cms.double(1.0), type = cms.int32(0)),
     cms.PSet( ecfAlgo = cms.string('antikt_algorithm'), Rparam = cms.double(2.0), nPoint = cms.int32(3), beta = cms.double(1.5), type = cms.int32(0)),
     cms.PSet( ecfAlgo = cms.string('antikt_algorithm'), Rparam = cms.double(2.0), nPoint = cms.int32(3), beta = cms.double(2.0), type = cms.int32(0)),
     cms.PSet( ecfAlgo = cms.string('antikt_algorithm'), Rparam = cms.double(2.0), nPoint = cms.int32(3), beta = cms.double(2.5), type = cms.int32(0)),
     cms.PSet( ecfAlgo = cms.string('antikt_algorithm'), Rparam = cms.double(2.0), nPoint = cms.int32(3), beta = cms.double(3.0), type = cms.int32(0)),

     cms.PSet( ecfAlgo = cms.string('antikt_algorithm'), Rparam = cms.double(2.0), nPoint = cms.int32(4), beta = cms.double(0.5), type = cms.int32(0)),
     cms.PSet( ecfAlgo = cms.string('antikt_algorithm'), Rparam = cms.double(2.0), nPoint = cms.int32(4), beta = cms.double(1.0), type = cms.int32(0)),
     cms.PSet( ecfAlgo = cms.string('antikt_algorithm'), Rparam = cms.double(2.0), nPoint = cms.int32(4), beta = cms.double(1.5), type = cms.int32(0)),
     cms.PSet( ecfAlgo = cms.string('antikt_algorithm'), Rparam = cms.double(2.0), nPoint = cms.int32(4), beta = cms.double(2.0), type = cms.int32(0)),
     cms.PSet( ecfAlgo = cms.string('antikt_algorithm'), Rparam = cms.double(2.0), nPoint = cms.int32(4), beta = cms.double(2.5), type = cms.int32(0)),
     cms.PSet( ecfAlgo = cms.string('antikt_algorithm'), Rparam = cms.double(2.0), nPoint = cms.int32(4), beta = cms.double(3.0), type = cms.int32(0)),
    ),

    Nsubjettiness = cms.VPSet(
     cms.PSet(beta = cms.double(0.5),  R0 = cms.double(0.8), Rcut = cms.double(0.8)),
     cms.PSet(beta = cms.double(0.75), R0 = cms.double(0.8), Rcut = cms.double(0.8)),
     cms.PSet(beta = cms.double(1.0),  R0 = cms.double(0.8), Rcut = cms.double(0.8)),
     cms.PSet(beta = cms.double(1.25), R0 = cms.double(0.8), Rcut = cms.double(0.8)),
     cms.PSet(beta = cms.double(1.5),  R0 = cms.double(0.8), Rcut = cms.double(0.8)),
     cms.PSet(beta = cms.double(2.0),  R0 = cms.double(0.8), Rcut = cms.double(0.8)),
     cms.PSet(beta = cms.double(2.5),  R0 = cms.double(0.8), Rcut = cms.double(0.8)),
     cms.PSet(beta = cms.double(3.0),  R0 = cms.double(0.8), Rcut = cms.double(0.8)),
    ),



)
