// FWCore tools for cfg parsing
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

// FastJet objects
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/contrib/SafeSubtractor.hh"
#include "fastjet/contrib/SoftKiller.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/ConstituentSubtractor.hh"
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/tools/JHTopTagger.hh"
#include "fastjet/Selector.hh"

// SubstrAna objects
#include "../include/MuonLoader.hh"
#include "../include/VTaggingVariables.h"
#include "../include/JetTools.h"
#include "../include/ShapeCorrectionTools.h"
#include "CMSTopTagger.hh"
#include "SubstrAna/Summer14/src/HEPTopTaggerWrapper.hh"

// JEC class
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

// ROOT objects
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TChain.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TGraph.h"

// c++ lybraries
#include <ctime>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cstdlib>

// Bacon Objects
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TVertex.hh"

// some typedef
using namespace std;
using namespace fastjet;
using namespace contrib;

typedef vector<float> vfloat;
typedef vector<bool> vbool;

////////////////////////////////////////
///// Global variables used ////////////
////////////////////////////////////////

// Object Processors
GenLoader       *fGen      = 0; 
PFLoader        *fPFCand   = 0; 
TClonesArray    *fJet      = 0;
TClonesArray    *fGenParticles = 0;

// jet clustering R size and tresholds
double jetR ;
double jetPtTresholdForGroomers, jetPtTresholdForTopTagging, genJetPtTresholdForTopTagging, genjetPtTresholdForGroomers;

// object for VTagging evaluation
VTaggingVariables* vtagger;

// parsing groomers parameter from cfg file
edm::ParameterSet   softKillerParam ;
std::vector<double> chargeParam ;
std::vector<edm::ParameterSet> softDropParam, trimmingParam, pruningParam, ecfParam,NsubjettinessParam ;
fastjet::JetAlgorithm algorithm_Trimming, algorithm_Pruning;

// matching thresholds 
double dRMatching ;
double dRLeptonCleaning;

// random seed
TRandom3 randNumber ;

//QGLikelihood calculator
std::string QGinputWeightFilePath;
QGLikelihoodCalculator* qgLikelihood, *qgLikelihoodCHS ;

//to compute Gen flavour
bool computeJetFlavour;

//////////////////////////////////////////////////////////////
///////// Structure used in order to fill output tree branches
//////////////////////////////////////////////////////////////

class GenJetInfo {

 public :

  GenJetInfo(){};
  ~GenJetInfo(){};

  int npu ;
  int npv ;

  vector<float> pt;
  vector<float> ptcorr;
  vector<float> ptcorrphil;
  vector<float> ptraw;
  vector<float> ptunc;

  vector<float> eta;
  vector<float> phi;
  vector<float> m;      // mass after JEC
  vector<float> mraw;

  vector<float> mclean;
  vector<float> ptclean;
  vector<float> ptconst;
  vector<float> mconst;

  vector<vector<float> > mtrim;
  vector<vector<float> > pttrim;
  vector<vector<float> > pttrimsafe;
  vector<vector<float> > mtrimsafe;

  vector<vector<float> > ptpruned;
  vector<vector<float> > mpruned;
  vector<vector<float> > ptprunedsafe;
  vector<vector<float> > mprunedsafe;

  vector<vector<float> > ptsoftdrop;
  vector<vector<float> > msoftdrop;
  vector<vector<float> > ptsoftdropsafe;
  vector<vector<float> > msoftdropsafe;

  vector<vector<float> > QGLikelihood_pr ;
  vector<vector<float> > QGLikelihood_pr_sub1 ;
  vector<vector<float> > QGLikelihood_pr_sub2 ;

  vector<vector<float> > pullangle;
  vector<vector<float> > pullmagnitude;

  vector<float> sdsymmetry ;
  vector<float> sddeltar ;
  vector<float> sdmu ;
  vector<float> sdenergyloss ;
  vector<float> sdarea ;
  vector<float> sdnconst ;
  vector<float> mfiltsoftdrop ;

  vector<float> hepmass ;
  vector<float> hepwmass ;
  vector<float> hepm01 ;
  vector<float> hepm02 ;
  vector<float> hepm12 ;
  vector<float> hepm12m012 ;
  vector<float> hepatanm02m01;

  // V-tagging base variables
  vector<vector<float> > tau1;
  vector<vector<float> > tau2;
  vector<vector<float> > tau3;

  vector<float> tau1_pr;
  vector<float> tau2_pr;
  vector<float> tau3_pr;

  vector<float> tau1_softdrop;
  vector<float> tau2_softdrop;
  vector<float> tau3_softdrop;

  vector<float> Qjets;

  vector<vector<float> > charge;

  vector<vector<float> > ecf;

  // constituents info
  vector<int>   nparticles;
  vector<int>   nneutrals;
  vector<int>   ncharged;

  vector<float> cmsmass ;
  vector<float> cmsminmass ;
  vector<float> cmshelicity ;
  vector<float> cmsnsubjets ;

  vector<int> jetflavour ;

};

class JetInfo : public GenJetInfo {

 public:

  JetInfo(){};
  ~JetInfo(){};

  // gen level info
  vector<float> ptgen;
  vector<float> etagen;
  vector<float> phigen;
  vector<float> mgen;
  vector<float> mrawgen;
  vector<float> mtrimgen;
  vector<float> mtrimsafegen;
  vector<int>   imatch;
  vector<int>   flavourgen;
  vector<float> msoftdropgen ;
  vector<float> msoftdropsafegen;  
  //matching to the Boson
  vector <bool> is_MatchedToBoson;

};

// ------------------------------------------------------------------------------------------
TTree* load(std::string iName) { 
  TFile *lFile = TFile::Open(iName.c_str());
  TTree *lTree = (TTree*) lFile->FindObjectAny("Events");
  return lTree;
}
 
// ------------------ Matching functions -------------- --> GenJets and RecoJets from GenJetInfo
int matchingIndexFromJetInfo(const PseudoJet& jet, const GenJetInfo& jetInfo) {
  float rmin = 9999.;
  int imatch = -1;
  for(unsigned int i = 0; i < (jetInfo.eta).size(); i++) {
    double dEta = fabs(jet.eta() - (jetInfo.eta)[i]);
    double dPhi = fabs(jet.phi() - (jetInfo.phi)[i]);
    if(dPhi > 2.*TMath::Pi()-dPhi) dPhi =  2.*TMath::Pi()-dPhi;
    float rtemp = sqrt(dEta*dEta+dPhi*dPhi);
    if ( rtemp > dRMatching ) continue;
    if ( rtemp < rmin ){
      rmin =  rtemp;
      imatch = i;
    }
  }
  return (imatch);
}

bool matchingIndex(const PseudoJet & jet, const PseudoJet & genjet, const bool & LeptonCleaning = false) {
  float rtemp = jet.delta_R(genjet);
  if( LeptonCleaning == false ){
    if ( rtemp < dRMatching ) return true;
    else return false;
  }
  else{
    if ( rtemp < dRLeptonCleaning) return true;
    else return false;
  }
}

bool IsMatchedToGenBoson(const vector<float>& eta, const vector<float>& phi, const PseudoJet& Jet) {
  bool IsMatched=false;
  for (unsigned int iGen =0; iGen < eta.size(); ++iGen){
    double dEta = fabs(eta.at(iGen) - (Jet.eta()));
    double dPhi = fabs(phi.at(iGen) - (Jet.phi()));
    if(dPhi > 2.*TMath::Pi()-dPhi) dPhi =  2.*TMath::Pi()-dPhi;
    float rtemp = sqrt(dEta*dEta+dPhi*dPhi);
    if ( rtemp < dRMatching ){
      IsMatched = true;
    }
  }
  return (IsMatched);
}

/// Selectors

Selector SelectorIsPupCharged(){
  return Selector(new SW_IsPupCharged());
}

Selector SelectorIsPupVertex(){
  return Selector(new SW_IsPupVertex());
}

//// Set the output tree structure
void setupGenTree(TTree *iTree, GenJetInfo &iJet, std::string iName) {

  iTree->Branch((iName+"npu"       ).c_str(),&iJet.npu       );
  iTree->Branch((iName+"npv"       ).c_str(),&iJet.npv       );

  iTree->Branch((iName+"pt"        ).c_str(),&iJet.pt        );  
  iTree->Branch((iName+"ptcorr"    ).c_str(),&iJet.ptcorr    );
  iTree->Branch((iName+"ptcorrphil").c_str(),&iJet.ptcorrphil);
  iTree->Branch((iName+"ptraw"     ).c_str(),&iJet.ptraw     );
  iTree->Branch((iName+"ptunc"     ).c_str(),&iJet.ptunc     );

  iTree->Branch((iName+"eta"       ).c_str(),&iJet.eta       );
  iTree->Branch((iName+"phi"       ).c_str(),&iJet.phi       );
  iTree->Branch((iName+"m"         ).c_str(),&iJet.m         );
  iTree->Branch((iName+"mraw"      ).c_str(),&iJet.mraw      );

  iTree->Branch((iName+"ptclean"   ).c_str(),&iJet.ptclean   );
  iTree->Branch((iName+"mclean"    ).c_str(),&iJet.mclean    );

  iTree->Branch((iName+"ptconst"   ).c_str(),&iJet.ptconst   );
  iTree->Branch((iName+"mconst"    ).c_str(),&iJet.mconst    );

  // trimming information   
  std::vector<edm::ParameterSet>::const_iterator itTrim = trimmingParam.begin();
  int iPos = 0 ;
  iJet.pttrim.resize(trimmingParam.size()) ;  
  iJet.mtrim.resize(trimmingParam.size()) ;
  iJet.pttrimsafe.resize(trimmingParam.size());
  iJet.mtrimsafe.resize(trimmingParam.size());
    
  for( ; itTrim != trimmingParam.end() ; ++itTrim){
   TString name ;
   name = Form("_Rtrim_%0.2f_Ptfrac_%0.2f",(*itTrim).getParameter<double>("R_trimming"),(*itTrim).getParameter<double>("PtFraction"));
   name.ReplaceAll(".","");
   iTree->Branch((iName+"pttrim"+std::string(name)     ).c_str(),"vector<float>",&iJet.pttrim[iPos]);
   iTree->Branch((iName+"mtrim"+std::string(name)      ).c_str(),"vector<float>",&iJet.mtrim[iPos]);  
   iTree->Branch((iName+"pttrimsafe"+std::string(name) ).c_str(),"vector<float>",&iJet.pttrimsafe[iPos]);
   iTree->Branch((iName+"mtrimsafe"+std::string(name)  ).c_str(),"vector<float>",&iJet.mtrimsafe[iPos]);
   iPos++ ;
  }
  
  // pruning informations --> QGLikelihood and pull information for each parameter combination  
  std::vector<edm::ParameterSet>::const_iterator itPruned = pruningParam.begin();
  iPos = 0 ;
  iJet.ptpruned.resize(pruningParam.size()) ;  
  iJet.mpruned.resize(pruningParam.size()) ;
  iJet.ptprunedsafe.resize(pruningParam.size());
  iJet.mprunedsafe.resize(pruningParam.size());
  iJet.QGLikelihood_pr.resize(pruningParam.size());
  iJet.QGLikelihood_pr_sub1.resize(pruningParam.size());
  iJet.QGLikelihood_pr_sub2.resize(pruningParam.size());
  iJet.pullangle.resize(pruningParam.size());
  iJet.pullmagnitude.resize(pruningParam.size());

  for( ; itPruned != pruningParam.end() ; ++itPruned){
   TString name ;
   name = Form("_zcut_%0.2f_R_cut_%0.2f",(*itPruned).getParameter<double>("z_cut"),(*itPruned).getParameter<double>("R_Cut"));
   name.ReplaceAll(".","");
   iTree->Branch((iName+"ptpruned"+std::string(name)     ).c_str(),"vector<float>",&iJet.ptpruned[iPos]);
   iTree->Branch((iName+"mpruned"+std::string(name)      ).c_str(),"vector<float>",&iJet.mpruned[iPos]);
   iTree->Branch((iName+"ptprunedsafe"+std::string(name) ).c_str(),"vector<float>",&iJet.ptprunedsafe[iPos]);
   iTree->Branch((iName+"mprunedsafe"+std::string(name)  ).c_str(),"vector<float>",&iJet.mprunedsafe[iPos]);
   iTree->Branch((iName+"QGLikelihood_pr"+std::string(name)  ).c_str(),"vector<float>",&iJet.QGLikelihood_pr[iPos]);
   iTree->Branch((iName+"QGLikelihood_pr_sub1"+std::string(name)  ).c_str(),"vector<float>",&iJet.QGLikelihood_pr_sub1[iPos]);
   iTree->Branch((iName+"QGLikelihood_pr_sub2"+std::string(name)  ).c_str(),"vector<float>",&iJet.QGLikelihood_pr_sub2[iPos]);
   iTree->Branch((iName+"pullangle"+std::string(name)  ).c_str(),"vector<float>",&iJet.pullangle[iPos]);
   iTree->Branch((iName+"pullmagnitude"+std::string(name)  ).c_str(),"vector<float>",&iJet.pullmagnitude[iPos]);
   iPos++ ;
  }

  // soft drop information
  std::vector<edm::ParameterSet>::const_iterator itsoftDrop = softDropParam.begin();
  iPos = 0 ;
  iJet.ptsoftdrop.resize(softDropParam.size()) ;  
  iJet.msoftdrop.resize(softDropParam.size()) ;
  iJet.ptsoftdropsafe.resize(softDropParam.size());
  iJet.msoftdropsafe.resize(softDropParam.size());

  for( ; itsoftDrop != softDropParam.end() ; ++itsoftDrop){
   TString name;
   name = Form("zcut%0.1f_beta%0.1f",(*itsoftDrop).getParameter<double>("symmetry_cut"),(*itsoftDrop).getParameter<double>("beta"));
   name.ReplaceAll(".","");
   iTree->Branch((iName+"ptsoftdrop"+std::string(name)     ).c_str(),"vector<float>",&iJet.ptsoftdrop[iPos]);
   iTree->Branch((iName+"msoftdrop"+std::string(name)      ).c_str(),"vector<float>",&iJet.msoftdrop[iPos]);
   iTree->Branch((iName+"ptsoftdropsafe"+std::string(name) ).c_str(),"vector<float>",&iJet.ptsoftdropsafe[iPos]);
   iTree->Branch((iName+"msoftdropsafe"+std::string(name)  ).c_str(),"vector<float>",&iJet.msoftdropsafe[iPos]);
   iPos++;
  }

  // other branches

  iTree->Branch((iName+"nparticles").c_str(),&iJet.nparticles);
  iTree->Branch((iName+"nneutrals" ).c_str(),&iJet.nneutrals);
  iTree->Branch((iName+"ncharged"  ).c_str(),&iJet.ncharged);


  iTree->Branch((iName+"sdsymmetry" ).c_str(),&iJet.sdsymmetry );
  iTree->Branch((iName+"sddeltar" ).c_str(),&iJet.sddeltar );
  iTree->Branch((iName+"sdmu" ).c_str(),&iJet.sdmu );
  iTree->Branch((iName+"sdenergyloss" ).c_str(),&iJet.sdenergyloss );
  iTree->Branch((iName+"sdarea" ).c_str(),&iJet.sdarea );
  iTree->Branch((iName+"sdnconst" ).c_str(),&iJet.sdnconst );
  iTree->Branch((iName+"mfiltsoftdrop" ).c_str(),&iJet.mfiltsoftdrop );

  // Nsubjettiness
  std::vector<edm::ParameterSet>::const_iterator itNsubjettiness = NsubjettinessParam.begin();
  iPos = 0 ;
  iJet.tau1.resize(NsubjettinessParam.size()) ;  
  iJet.tau2.resize(NsubjettinessParam.size()) ;  
  iJet.tau3.resize(NsubjettinessParam.size()) ;  

  for( ; itNsubjettiness != NsubjettinessParam.end() ; ++itNsubjettiness){
   TString name ;
   name = Form("_beta_%0.1f",(*itNsubjettiness).getParameter<double>("beta"));
   name.ReplaceAll(".","");
   iTree->Branch((iName+"tau1"+std::string(name) ).c_str(),"vector<float>",&iJet.tau1[iPos]);
   iTree->Branch((iName+"tau2"+std::string(name) ).c_str(),"vector<float>",&iJet.tau2[iPos]);
   iTree->Branch((iName+"tau3"+std::string(name) ).c_str(),"vector<float>",&iJet.tau3[iPos]);
   iPos++ ;
  }


  iTree->Branch((iName+"tau1_pr"  ).c_str(),&iJet.tau1_pr);
  iTree->Branch((iName+"tau2_pr"  ).c_str(),&iJet.tau2_pr);
  iTree->Branch((iName+"tau3_pr"  ).c_str(),&iJet.tau3_pr);

  iTree->Branch((iName+"tau1_softdrop"  ).c_str(),&iJet.tau1_softdrop);
  iTree->Branch((iName+"tau2_softdrop"  ).c_str(),&iJet.tau2_softdrop);
  iTree->Branch((iName+"tau3_softdrop"  ).c_str(),&iJet.tau3_softdrop);

  iTree->Branch((iName+"Qjets"  ).c_str(),&iJet.Qjets);

  // jet charge
  std::vector<double>::const_iterator itCharge = chargeParam.begin();
  iPos = 0 ;
  iJet.charge.resize(chargeParam.size()) ;  
  for( ; itCharge !=chargeParam.end() ; ++itCharge){
   TString name;
   name = Form("_k%0.1f",(*itCharge));
   name.ReplaceAll(".","");
   iTree->Branch((iName+"charge"+std::string(name)     ).c_str(),"vector<float>",&iJet.charge[iPos]);
   iPos++;
  }

  // ECF
  std::vector<edm::ParameterSet>::const_iterator itECF = ecfParam.begin();
  iPos = 0 ;
  iJet.ecf.resize(ecfParam.size()) ;  
  for( ; itECF !=ecfParam.end() ; ++itECF){
   TString name;
   name = Form("_nPoint_%d_beta_%0.1f",(*itECF).getParameter<int>("nPoint"),(*itECF).getParameter<double>("beta"));
   name.ReplaceAll(".","");
   iTree->Branch((iName+"ecf"+std::string(name)).c_str(),"vector<float>",&iJet.ecf[iPos]);
   iPos++;
  }

  iTree->Branch((iName+"hepmass" ).c_str(),&iJet.hepmass );
  iTree->Branch((iName+"hepwmass" ).c_str(),&iJet.hepwmass );
  iTree->Branch((iName+"hepm01" ).c_str(),&iJet.hepm01 );
  iTree->Branch((iName+"hepm02" ).c_str(),&iJet.hepm02 );
  iTree->Branch((iName+"hepm12" ).c_str(),&iJet.hepm12 );
  iTree->Branch((iName+"hepm12m012" ).c_str(),&iJet.hepm12m012 );
  iTree->Branch((iName+"hepatanm02m01" ).c_str(),&iJet.hepatanm02m01 );

  iTree->Branch((iName+"cmsmass" ).c_str(),&iJet.cmsmass );
  iTree->Branch((iName+"cmsminmass" ).c_str(),&iJet.cmsminmass );
  iTree->Branch((iName+"cmshelicity" ).c_str(),&iJet.cmshelicity );
  iTree->Branch((iName+"cmsnsubjets" ).c_str(),&iJet.cmsnsubjets );

}

void setupTree(TTree *iTree, JetInfo &iJet, std::string iName) {

  setupGenTree(iTree,iJet,iName); 

  // gen info
  iTree->Branch((iName+"ptgen"       ).c_str(),&iJet.ptgen       );
  iTree->Branch((iName+"etagen"      ).c_str(),&iJet.etagen      );
  iTree->Branch((iName+"phigen"      ).c_str(),&iJet.phigen      );
  iTree->Branch((iName+"mgen"        ).c_str(),&iJet.mgen        );
  iTree->Branch((iName+"mrawgen"     ).c_str(),&iJet.mrawgen     );//needed?
  iTree->Branch((iName+"mtrimgen"    ).c_str(),&iJet.mtrimgen    );//needed?
  iTree->Branch((iName+"mtrimsafegen").c_str(),&iJet.mtrimsafegen);//needed?
  iTree->Branch((iName+"imatch"      ).c_str(),&iJet.imatch      );
  iTree->Branch((iName+"flavourgen"  ).c_str(),&iJet.flavourgen  );  
  iTree->Branch((iName+"msoftdropgen"          ).c_str(),&iJet.msoftdropgen          );
  iTree->Branch((iName+"msoftdropsafegen"      ).c_str(),&iJet.msoftdropsafegen      );  
  //matched to the boson
  iTree->Branch((iName+"is_MatchedToBoson"      ).c_str(),&iJet.is_MatchedToBoson      );
   
}

// clear tree structure content at the beginning of each event
void clear(GenJetInfo &iJet) {

  iJet.npu  = -1;
  iJet.npv  = -1;

  iJet.pt         .clear();
  iJet.ptcorr     .clear();
  iJet.ptcorrphil .clear();
  iJet.ptraw      .clear();
  iJet.ptunc      .clear();

  iJet.jetflavour .clear();

  iJet.eta        .clear();
  iJet.phi        .clear();
  iJet.m          .clear();
  iJet.mraw       .clear();

  iJet.mclean     .clear();
  iJet.ptclean    .clear();
  iJet.ptconst    .clear();
  iJet.mconst     .clear();

  iJet.nparticles .clear();
  iJet.nneutrals  .clear();
  iJet.ncharged   .clear();

  for(unsigned int iTrim = 0; iTrim < iJet.pttrim.size(); iTrim++){
    iJet.pttrim.at(iTrim).clear();
    iJet.mtrim.at(iTrim).clear();
    iJet.pttrimsafe.at(iTrim).clear();
    iJet.mtrimsafe.at(iTrim).clear();
  } 


  for(unsigned int iPruned = 0; iPruned < iJet.ptpruned.size(); iPruned++){
    iJet.ptpruned.at(iPruned).clear();
    iJet.mpruned.at(iPruned).clear();
    iJet.ptprunedsafe.at(iPruned).clear();
    iJet.mprunedsafe.at(iPruned).clear();
    iJet.QGLikelihood_pr.at(iPruned).clear();
    iJet.QGLikelihood_pr_sub1.at(iPruned).clear();
    iJet.QGLikelihood_pr_sub2.at(iPruned).clear();
    iJet.pullangle.at(iPruned).clear();
    iJet.pullmagnitude.at(iPruned).clear();

  } 

  for(unsigned int iSoft = 0; iSoft < iJet.ptsoftdrop.size(); iSoft++){
    iJet.ptsoftdrop.at(iSoft).clear();
    iJet.ptsoftdropsafe.at(iSoft).clear();
    iJet.msoftdrop.at(iSoft).clear();
    iJet.msoftdropsafe.at(iSoft).clear();
  } 


  iJet.sdsymmetry.clear();
  iJet.sddeltar.clear();
  iJet.sdmu.clear();
  iJet.sdenergyloss.clear();
  iJet.sdarea.clear();
  iJet.sdnconst.clear();
  iJet.mfiltsoftdrop.clear();

  for(unsigned int iNsubjettiness = 0; iNsubjettiness < iJet.tau1.size(); iNsubjettiness++){
    iJet.tau1.at(iNsubjettiness).clear();
    iJet.tau2.at(iNsubjettiness).clear();
    iJet.tau3.at(iNsubjettiness).clear();
  }

  iJet.tau1_pr.clear();
  iJet.tau2_pr.clear();
  iJet.tau3_pr.clear();

  iJet.tau1_softdrop.clear();
  iJet.tau2_softdrop.clear();
  iJet.tau3_softdrop.clear();

  iJet.Qjets.clear();

  for( unsigned int iCharge = 0 ; iCharge < iJet.charge.size() ; iCharge++)
    iJet.charge.at(iCharge).clear();

  for( unsigned int iECF = 0 ; iECF < iJet.ecf.size() ; iECF++)
    iJet.ecf.at(iECF).clear();

  iJet.hepmass .clear();
  iJet.hepwmass .clear();
  iJet.hepm01 .clear();
  iJet.hepm02 .clear();
  iJet.hepm12 .clear();
  iJet.hepm12m012 .clear();
  iJet.hepatanm02m01 .clear();
  iJet.cmsmass .clear();
  iJet.cmsminmass .clear();
  iJet.cmshelicity .clear();
  iJet.cmsnsubjets .clear();

}

void clear(JetInfo &iJet) {

  clear((GenJetInfo &)iJet);

  iJet.ptgen       .clear();
  iJet.etagen      .clear();
  iJet.phigen      .clear();
  iJet.mgen        .clear();
  iJet.mrawgen     .clear();
  iJet.mtrimgen    .clear();
  iJet.mtrimsafegen.clear();

  iJet.imatch      .clear();
  iJet.flavourgen  .clear();
  iJet.is_MatchedToBoson.clear();

  iJet.msoftdropgen .clear();
  iJet.msoftdropsafegen.clear();
}

//set the gen jet info in the output tree
void setGenJet( PseudoJet &iJet, // current jet 
                GenJetInfo &iJetI, // output structure
                JetMedianBackgroundEstimator bge_rho, JetMedianBackgroundEstimator bge_rhom, JetMedianBackgroundEstimator bge_rhoC, // info for safe area or constituent subtraction
                vector<JetCleanser> &cleanser_vect, // cleanser objects
                bool is_leadingJet, double rho) {



  // -- area-median subtractor  (safe area subtractor)
  contrib::SafeAreaSubtractor *area_subtractor = 0;
  area_subtractor = new contrib::SafeAreaSubtractor(&bge_rho, &bge_rhom);
  PseudoJet lCorr =  (*area_subtractor)(iJet); // correct the jet for safe are subtraction
  
  // -- constituent subtractor
  contrib::ConstituentSubtractor *const_subtractor = 0;
  const_subtractor = new contrib::ConstituentSubtractor(&bge_rhoC);
  (*const_subtractor).use_common_bge_for_rho_and_rhom(true);
  PseudoJet lConstit = (*const_subtractor)(iJet); // correct the jet for constituent subtraction

  // -- cleansing 
  vector<PseudoJet> neutrals,chargedLV,chargedPU,jetParticles,ghosts;
  SelectorIsPureGhost().sift(iJet.constituents(), ghosts, jetParticles);
  getConstitsForCleansing(jetParticles,neutrals,chargedLV,chargedPU);

  if(is_leadingJet){ // cleansing only for leading jet pt
	for(Int_t i=0; i<Int_t(cleanser_vect.size());i++){
           PseudoJet lClean = cleanser_vect[i](neutrals,chargedLV,chargedPU); // use cleansing
           (iJetI.ptclean   ).push_back(lClean    .pt());
           (iJetI.mclean    ).push_back(lClean    .m());
	}
  }
  

  // -- trimming  
  vector<PseudoJet> lTrim ;
  vector<PseudoJet> lTrimSafe ;
  if(iJet.pt() >= genjetPtTresholdForGroomers){ // if over threshold
   std::vector<edm::ParameterSet>::const_iterator itTrim = trimmingParam.begin();
   for( ; itTrim != trimmingParam.end() ; ++itTrim){
    fastjet::Filter trimmer(fastjet::Filter(fastjet::JetDefinition(get_algo((*itTrim).getParameter<string>("trimAlgo")),(*itTrim).getParameter<double>("R_trimming")), fastjet::SelectorPtFractionMin((*itTrim).getParameter<double>("PtFraction"))));
    lTrim.push_back((trimmer)(iJet));
    trimmer.set_subtractor(area_subtractor);
    lTrimSafe.push_back((trimmer)(iJet));
   }
  }
  
  // -- pruning
  vector<PseudoJet> lPruned ;
  vector<PseudoJet> lPrunedSafe ;
  if(iJet.pt() >= genjetPtTresholdForGroomers){ // if over threshold
   std::vector<edm::ParameterSet>::const_iterator itPruned = pruningParam.begin();
   for( ; itPruned != pruningParam.end() ; ++itPruned){
    JetDefinition jet_def_Pruning(get_algo((*itPruned).getParameter<string>("pruneAlgo")), (*itPruned).getParameter<double>("R_jet_def_pruning"));
    Pruner pruner(jet_def_Pruning,(*itPruned).getParameter<double>("z_cut"), (*itPruned).getParameter<double>("R_Cut"));
    PseudoJet jetTemp = pruner(iJet) ;
    lPruned.push_back(jetTemp);
    lPrunedSafe.push_back((*area_subtractor)(jetTemp));
   }
  }
  
  // -- softdrop
  vector<PseudoJet> lSoftDropped ;
  vector<PseudoJet> lSoftDroppedSafe ;

  double  SoftDropedSymmetry   = -1.0;
  double  SoftDropedDR         = -1.0;
  double  SoftDropedMassDrop   = -1.0;
  double  SoftDropedEnergyLoss = -1.0;
  double  SoftDropedArea       = -1.0;
  double  SoftDropedNconst     = -1.0;
  PseudoJet filtered_softdropped_jet; 
   

  if(iJet.pt() >= genjetPtTresholdForGroomers){ // if over threshold
   std::vector<edm::ParameterSet>::const_iterator itSoft = softDropParam.begin();
   for( ; itSoft != softDropParam.end() ; ++itSoft){
    contrib::SoftDrop softdrop((*itSoft).getParameter<double>("beta"),(*itSoft).getParameter<double>("symmetry_cut"), (*itSoft).getParameter<double>("R0"));
    lSoftDropped.push_back(softdrop(iJet));   
    softdrop.set_subtractor(area_subtractor);
    lSoftDroppedSafe.push_back(softdrop(iJet));
   } 

   if (lSoftDroppedSafe.at(0)!=0 and lSoftDroppedSafe.at(0).m()>0.0){
    SoftDropedSymmetry   = lSoftDroppedSafe.at(0).structure_of<contrib::SoftDrop>().symmetry(); 
    SoftDropedDR         = lSoftDroppedSafe.at(0).structure_of<contrib::SoftDrop>().delta_R();  
    SoftDropedMassDrop   = lSoftDroppedSafe.at(0).structure_of<contrib::SoftDrop>().mu();       
    SoftDropedEnergyLoss = 1-lSoftDroppedSafe.at(0).pt()/iJet.pt();                           
    SoftDropedArea       = lSoftDroppedSafe.at(0).area() ;                                     
    SoftDropedNconst     = lSoftDroppedSafe.at(0).constituents().size() ;                      

    // filter jet dynamically based on deltaR between subjets (arXiv:0802.2470)
    double dyn_Rfilt = min(0.3, SoftDropedDR*0.5);
    int    dyn_nfilt = 3;
    Filter filtersoft(dyn_Rfilt, SelectorNHardest(dyn_nfilt));
    filtered_softdropped_jet = filtersoft(lSoftDroppedSafe.at(0));
   }
  }

  // -- Top Tag
  PseudoJet iJetCA;
  double cmsttJetMass = -1;
  double cmsttMinMass = -1;
  double cmsttHelicity = -1;
  double cmsttNsubjets = -1;
  double hepttJetMass    = -1; 
  double hepttWMass      = -1; 
  double hepttM01        = -1; 
  double hepttM02        = -1; 
  double hepttM12        = -1; 
  double hepttM12M012    = -1; 
  double hepttAtanM02M01 = -1; 

  if (iJet.pt() > genJetPtTresholdForTopTagging){

    // -- recluster jet CA
    JetDefinition jet_def_CA (fastjet::cambridge_algorithm, jetR*10); //large R to cluster all constituents of original jet
    fastjet::ClusterSequence cs_Recluster (jetParticles,jet_def_CA);
    vector<fastjet::PseudoJet> jets_Recluster = sorted_by_pt(cs_Recluster.inclusive_jets());
    iJetCA = jets_Recluster[0];

    // -- HEP Top Tagger 
    double mass_drop_threshold = 0.8;
    double max_subjet_mass     = 30;
    bool use_subjet_mass_cuts  = false;
    HEPTopTagger hep_top_tagger(mass_drop_threshold, max_subjet_mass, use_subjet_mass_cuts);
    
    PseudoJet hep_top_candidate   = hep_top_tagger( iJet );

    if (hep_top_candidate != 0){
      PseudoJet W =     hep_top_candidate.structure_of<HEPTopTagger>().W();
      PseudoJet W1 =    hep_top_candidate.structure_of<HEPTopTagger>().W1();
      PseudoJet W2 =    hep_top_candidate.structure_of<HEPTopTagger>().W2();
      PseudoJet non_W = hep_top_candidate.structure_of<HEPTopTagger>().non_W();

      vector<PseudoJet> all_subjets;
      all_subjets.push_back(W1);
      all_subjets.push_back(W2);
      all_subjets.push_back(non_W);
      all_subjets = sorted_by_pt(all_subjets);

      PseudoJet sum012 = all_subjets[0]+all_subjets[1]+all_subjets[2];
      PseudoJet sum01 = all_subjets[0]+all_subjets[1];
      PseudoJet sum02 = all_subjets[0]+all_subjets[2];
      PseudoJet sum12 = all_subjets[1]+all_subjets[2];

      hepttJetMass       = hep_top_candidate.m();
      hepttWMass         = W.m();
      hepttM01           = sum01.m();
      hepttM02           = sum02.m();
      hepttM12           = sum12.m();
      if ( sum012.m()!=0 ) hepttM12M012     = sum12.m() / sum012.m() ;
      if ( sum01.m()!=0 )  hepttAtanM02M01  = atan( sum02.m() / sum01.m() ) ;
    }

    // -- CMS Top Tagger
    double cms_delta_p = 0.05;
    double cms_delta_r = 0.4;
    double A           = 0.0004;

    CMSTopTagger cms_top_tagger(cms_delta_p, cms_delta_r, A);
    PseudoJet cms_top_candidate = cms_top_tagger( iJetCA );

    if (cms_top_candidate != 0){
        vector<PseudoJet> kept_subjets0 = cms_top_candidate.structure_of<CMSTopTagger>().W().pieces();
        vector<PseudoJet> kept_subjets1 = cms_top_candidate.structure_of<CMSTopTagger>().non_W().pieces();
        vector<PseudoJet> all_subjets = kept_subjets0;
        all_subjets.insert( all_subjets.end(), kept_subjets1.begin(), kept_subjets1.end() );

        cmsttJetMass = cms_top_candidate.m();
        cmsttMinMass = cms_top_candidate.structure_of<CMSTopTagger>().W().m();
        cmsttHelicity = cms_top_candidate.structure_of<CMSTopTagger>().cos_theta_W();
        cmsttNsubjets = all_subjets.size();
    }
  }

         
  // -- Fill jet info
  (iJetI.pt        ).push_back(lCorr     .pt());  // pt after safe area subtraction 
  (iJetI.ptcorr    ).push_back(iJet      .pt());  // jet pt no JEC for gen
  (iJetI.ptcorrphil).push_back(iJet      .pt());  // jet pt no phill JEC for gen
  (iJetI.ptraw     ).push_back(iJet      .pt());  // raw pt
  (iJetI.ptunc     ).push_back(0.);

  (iJetI.eta       ).push_back(iJet      .eta());
  (iJetI.phi       ).push_back(iJet      .phi());
  (iJetI.mraw      ).push_back(iJet      .m());
  (iJetI.m         ).push_back(lCorr     .m());

  (iJetI.ptconst   ).push_back(lConstit  .pt());
  (iJetI.mconst    ).push_back(lConstit  .m());

  if(computeJetFlavour) (iJetI.jetflavour).push_back(computeGenJetFlavour(iJet,fGenParticles,jetR)); // make jet parton association for the flavour
   
  if(iJet.pt() >  genjetPtTresholdForGroomers){ // trimming information
   for( unsigned int iTrim = 0 ; iTrim < lTrim.size() ; iTrim++){
    iJetI.pttrim.at(iTrim).push_back(lTrim.at(iTrim).pt());
    iJetI.mtrim.at(iTrim).push_back(lTrim.at(iTrim).m());
    iJetI.pttrimsafe.at(iTrim).push_back(lTrimSafe.at(iTrim).pt());
    iJetI.mtrimsafe.at(iTrim).push_back(lTrimSafe.at(iTrim).m());
   }

   for( unsigned int iPrun = 0 ; iPrun < lPruned.size() ; iPrun++){// prunung information
    iJetI.mpruned.at(iPrun).push_back(lPruned.at(iPrun).m());
    iJetI.ptpruned.at(iPrun).push_back(lPruned.at(iPrun).pt());
    iJetI.mprunedsafe.at(iPrun).push_back(lPrunedSafe.at(iPrun).m());
    iJetI.ptprunedsafe.at(iPrun).push_back(lPrunedSafe.at(iPrun).pt());
   }

   for( unsigned int iSoft = 0 ; iSoft < lSoftDropped.size() ; iSoft++){ // prunung information
    iJetI.msoftdrop.at(iSoft).push_back(lSoftDropped.at(iSoft).m());
    iJetI.ptsoftdrop.at(iSoft).push_back(lSoftDropped.at(iSoft).pt());
    iJetI.msoftdropsafe.at(iSoft).push_back(lSoftDroppedSafe.at(iSoft).m());
    iJetI.ptsoftdropsafe.at(iSoft).push_back(lSoftDroppedSafe.at(iSoft).pt());
   }
  }

  (iJetI.nparticles).push_back(jetParticles.size()); // filter the ghost taking only real particles
  (iJetI.nneutrals ).push_back(neutrals.size()); 
  (iJetI.ncharged  ).push_back(chargedLV.size()+chargedPU.size());

  (iJetI.sdsymmetry ).push_back( SoftDropedSymmetry );
  (iJetI.sddeltar ).push_back( SoftDropedDR );
  (iJetI.sdmu ).push_back( SoftDropedMassDrop );
  (iJetI.sdenergyloss ).push_back( SoftDropedEnergyLoss );
  (iJetI.sdarea ).push_back( SoftDropedArea );
  (iJetI.sdnconst ).push_back( SoftDropedNconst );
  (iJetI.mfiltsoftdrop ).push_back( filtered_softdropped_jet.m() );

  (iJetI.hepmass ).push_back( hepttJetMass );
  (iJetI.hepwmass ).push_back( hepttWMass );
  (iJetI.hepm01 ).push_back( hepttM01 );
  (iJetI.hepm02 ).push_back( hepttM02 );
  (iJetI.hepm12 ).push_back( hepttM12 );
  (iJetI.hepm12m012 ).push_back( hepttM12M012 );
  (iJetI.hepatanm02m01).push_back( hepttAtanM02M01 );

  (iJetI.cmsmass ).push_back(cmsttJetMass );
  (iJetI.cmsminmass ).push_back(cmsttMinMass );
  (iJetI.cmshelicity ).push_back(cmsttHelicity );
  (iJetI.cmsnsubjets ).push_back(cmsttNsubjets );

  if(iJet.pt() > genjetPtTresholdForGroomers){ // if over treshold calculate other Vtagging observables

   vtagger->setInputJet(iJet); 

   for(unsigned int iNsubjettiness = 0; iNsubjettiness < NsubjettinessParam.size() ; iNsubjettiness++){ // N-subjettines
      (iJetI.tau1 ).at(iNsubjettiness).push_back(vtagger->computeNSubJettines(1,NsubjettinessParam.at(iNsubjettiness).getParameter<double>("beta"),NsubjettinessParam.at(iNsubjettiness).getParameter<double>("R0"),NsubjettinessParam.at(iNsubjettiness).getParameter<double>("Rcut")));
      (iJetI.tau2 ).at(iNsubjettiness).push_back(vtagger->computeNSubJettines(2,NsubjettinessParam.at(iNsubjettiness).getParameter<double>("beta"),NsubjettinessParam.at(iNsubjettiness).getParameter<double>("R0"),NsubjettinessParam.at(iNsubjettiness).getParameter<double>("Rcut")));
      (iJetI.tau3 ).at(iNsubjettiness).push_back(vtagger->computeNSubJettines(3,NsubjettinessParam.at(iNsubjettiness).getParameter<double>("beta"),NsubjettinessParam.at(iNsubjettiness).getParameter<double>("R0"),NsubjettinessParam.at(iNsubjettiness).getParameter<double>("Rcut")));
    }

   (iJetI.Qjets).push_back(vtagger->computeQjets(100,25,randNumber.Uniform(0.,10000),jetR)); // compute Qjets
   
   for( unsigned int iECF = 0; iECF < ecfParam.size() ; iECF++) // compute ECF
      iJetI.ecf.at(iECF).push_back(vtagger->computeECF(get_algo(ecfParam.at(iECF).getParameter<string>("ecfAlgo")),ecfParam.at(iECF).getParameter<double>("Rparam"),ecfParam.at(iECF).getParameter<int>("nPoint"),ecfParam.at(iECF).getParameter<double>("beta"),ecfParam.at(iECF).getParameter<int>("type")));
   

   for( unsigned int iCharge = 0; iCharge < chargeParam.size() ; iCharge++) // compute jet charge
    iJetI.charge.at(iCharge).push_back(vtagger->computeJetChargeReco(chargeParam[iCharge]));

   // using pruned subjets
   vector<PseudoJet> subjets_pruned ;
   for( unsigned int iPrun = 0 ; iPrun < lPruned.size() ; iPrun++){
      if(lPruned.at(iPrun).constituents().size() > 1){
       subjets_pruned = lPruned.at(iPrun).associated_cluster_sequence()->exclusive_subjets(lPruned.at(iPrun),2); // take the pieces
       subjets_pruned = sorted_by_pt(subjets_pruned); // sort by pt
       if(subjets_pruned.at(0).pt()>0 and subjets_pruned.at(1).pt()>0){ // both with pt > 0 --> compute pull
        iJetI.pullangle.at(iPrun).push_back(vtagger->computePullAngle(subjets_pruned,jetR));  
        iJetI.pullmagnitude.at(iPrun).push_back(TMath::Sqrt(iJetI.pullangle.at(iPrun).back()*iJetI.pullangle.at(iPrun).back()+subjets_pruned.at(0).rapidity()*subjets_pruned.at(0).rapidity()));
       }
       else{
	iJetI.pullangle.at(iPrun).push_back(999);
	iJetI.pullmagnitude.at(iPrun).push_back(999);
       }						       
      }
      else{
	iJetI.pullangle.at(iPrun).push_back(999);
	iJetI.pullmagnitude.at(iPrun).push_back(999);
      }
   }

   if(lPruned.at(0).pt() > 0){ // if exsist  the pruned jet with non null pt
     (iJetI.tau1_pr ).push_back(vtagger->computeNSubJettines(1,1.,jetR,jetR));
     (iJetI.tau2_pr ).push_back(vtagger->computeNSubJettines(2,1.,jetR,jetR));
     (iJetI.tau3_pr ).push_back(vtagger->computeNSubJettines(3,1.,jetR,jetR));
   }   
   else{
     (iJetI.tau1_pr ).push_back(999);
     (iJetI.tau2_pr ).push_back(999);
     (iJetI.tau3_pr ).push_back(999);
   }

   // QGLO for pruned jet 
   for( unsigned int iPrun = 0 ; iPrun < lPruned.size() ; iPrun++){
     if(lPruned.at(iPrun).pt() > 0 ) iJetI.QGLikelihood_pr.at(iPrun).push_back(vtagger->computeQGLikelihood(qgLikelihoodCHS,1.,rho));
      else iJetI.QGLikelihood_pr.at(iPrun).push_back(999);      
   }

   subjets_pruned.clear();

   for( unsigned int iPrun = 0 ; iPrun < lPruned.size() ; iPrun++){
     if(lPruned.at(iPrun).constituents().size() > 1){
      subjets_pruned = lPruned.at(iPrun).associated_cluster_sequence()->exclusive_subjets(lPruned.at(iPrun),2);
      subjets_pruned = sorted_by_pt(subjets_pruned);

      if(subjets_pruned.at(0).pt() > 0){
       vtagger->setInputJet(subjets_pruned.at(0));   
       iJetI.QGLikelihood_pr_sub1.at(iPrun).push_back(vtagger->computeQGLikelihood(qgLikelihoodCHS,1.,rho));
      }
      else iJetI.QGLikelihood_pr_sub1.at(iPrun).push_back(999);
     
      if(subjets_pruned.at(1).pt() > 0){
       vtagger->setInputJet(subjets_pruned.at(1));   
       iJetI.QGLikelihood_pr_sub2.at(iPrun).push_back(vtagger->computeQGLikelihood(qgLikelihoodCHS,1.,rho));
      }
      else iJetI.QGLikelihood_pr_sub2.at(iPrun).push_back(999);     
     }
     else{
       iJetI.QGLikelihood_pr_sub1.at(iPrun).push_back(999);
       iJetI.QGLikelihood_pr_sub2.at(iPrun).push_back(999);
    }
   }
   
   // soft drop info
   vtagger->setInputJet(lSoftDropped.at(0)); 
   if(lSoftDropped.at(0).pt() > 0 ){
     (iJetI.tau1_softdrop ).push_back(vtagger->computeNSubJettines(1,1.,jetR,jetR));
     (iJetI.tau2_softdrop ).push_back(vtagger->computeNSubJettines(2,1.,jetR,jetR));
     (iJetI.tau3_softdrop ).push_back(vtagger->computeNSubJettines(3,1.,jetR,jetR));
   }
   else{
     (iJetI.tau1_softdrop ).push_back(999);
     (iJetI.tau2_softdrop ).push_back(999);
     (iJetI.tau3_softdrop ).push_back(999);
   }
  }
  else {

    for( unsigned int iCharge = 0; iCharge < chargeParam.size() ; iCharge++)
      iJetI.charge.at(iCharge).push_back(999.);

    for( unsigned int iECF = 0; iECF < ecfParam.size() ; iECF++)
      iJetI.ecf.at(iECF).push_back(999.);

    for( unsigned int iPrun = 0 ; iPrun < lPruned.size() ; iPrun++){
      iJetI.QGLikelihood_pr.at(iPrun).push_back(999.);
      iJetI.QGLikelihood_pr_sub1.at(iPrun).push_back(999.);
      iJetI.QGLikelihood_pr_sub2.at(iPrun).push_back(999.);
      iJetI.pullangle.at(iPrun).push_back(999.);       
      iJetI.pullmagnitude.at(iPrun).push_back(999.);       
    }

    for(unsigned int iNsubjettiness = 0; iNsubjettiness < NsubjettinessParam.size() ; iNsubjettiness++){
      (iJetI.tau1 ).at(iNsubjettiness).push_back(999.);
      (iJetI.tau2 ).at(iNsubjettiness).push_back(999.);
      (iJetI.tau3 ).at(iNsubjettiness).push_back(999.);
    }

    (iJetI.Qjets).push_back(999.);

    (iJetI.tau1_pr ).push_back(999.);
    (iJetI.tau2_pr ).push_back(999.);
    (iJetI.tau3_pr ).push_back(999.);

    (iJetI.tau1_softdrop ).push_back(999.);
    (iJetI.tau2_softdrop ).push_back(999.);
    (iJetI.tau3_softdrop ).push_back(999.);
  }
  
}

// Set Reco Jet variables 
void setRecoJet(PseudoJet &iJet, // input reco Jet 
                JetInfo &iJetI,  // structure to store reco jet info
                GenJetInfo& iGenJetI, // structure to store reco jet info 
                JetMedianBackgroundEstimator bge_rho, JetMedianBackgroundEstimator bge_rhom, JetMedianBackgroundEstimator bge_rhoC, 
                bool isCHS, // apply CHS or not
                FactorizedJetCorrector *iJetCorr, JetCorrectionUncertainty *iJetUnc, // standard JEC and uncertainty
                std::vector<TGraph*> &iCorr, // Phill corrections
                vector<JetCleanser> &cleanser_vect, // cleansing
                bool is_leadingJet, double rho, 
                vfloat eta_Boson, vfloat phi_Boson, // matched vector boson information
                const bool & isPuppi = false, bool isMC = true) {


  // -- area-median subtractor  ( safe area subtractor )
  contrib::SafeAreaSubtractor *area_subtractor = 0;
  if(!isCHS || fabs(iJet.eta()) > 2.5) area_subtractor = new contrib::SafeAreaSubtractor(&bge_rho, &bge_rhom);
  if( isCHS && fabs(iJet.eta()) < 2.5) area_subtractor = new contrib::SafeAreaSubtractor(&bge_rho, &bge_rhom,SelectorIsPupCharged(),SelectorIsPupVertex()); // within tracker used CHS before apply 4V substraction
  PseudoJet lCorr =  (*area_subtractor)(iJet);
  
  // -- constituent subtractor
  contrib::ConstituentSubtractor *const_subtractor = 0;
  const_subtractor = new contrib::ConstituentSubtractor(&bge_rhoC);
  (*const_subtractor).use_common_bge_for_rho_and_rhom(true);
  PseudoJet lConstit = (*const_subtractor)(iJet);

  // -- cleansing 
  vector<PseudoJet> neutrals,chargedLV,chargedPU,jetParticles,ghosts;
  SelectorIsPureGhost().sift(iJet.constituents(), ghosts, jetParticles);
  getConstitsForCleansing(jetParticles,neutrals,chargedLV,chargedPU);
  if(is_leadingJet){
	for(Int_t i=0; i<Int_t(cleanser_vect.size());i++){
          PseudoJet     lClean = cleanser_vect[i](neutrals,chargedLV,chargedPU); // use cleansing
           (iJetI.ptclean   ).push_back(lClean    .pt());
           (iJetI.mclean    ).push_back(lClean    .m());
	}
  }

  // -- trimming information
  vector<PseudoJet> lTrim ;
  vector<PseudoJet> lTrimSafe ;

  if (iJet.pt() > jetPtTresholdForGroomers){
    std::vector<edm::ParameterSet>::const_iterator itTrim = trimmingParam.begin();
    for( ; itTrim != trimmingParam.end() ; ++itTrim){
      fastjet::Filter trimmer( fastjet::Filter(fastjet::JetDefinition(get_algo((*itTrim).getParameter<string>("trimAlgo")),(*itTrim).getParameter<double>("R_trimming")), fastjet::SelectorPtFractionMin((*itTrim).getParameter<double>("PtFraction"))));
     lTrim.push_back((trimmer)(iJet));
     trimmer.set_subtractor(area_subtractor);
     lTrimSafe.push_back((trimmer)(iJet));
    }
  }

  // -- pruning information
  vector<PseudoJet> lPruned ;
  vector<PseudoJet> lPrunedSafe ;

  if (iJet.pt() > jetPtTresholdForGroomers){
    std::vector<edm::ParameterSet>::const_iterator itPruned = pruningParam.begin();
    for( ; itPruned != pruningParam.end() ; ++itPruned){
     JetDefinition jet_def_Pruning(get_algo((*itPruned).getParameter<string>("pruneAlgo")), (*itPruned).getParameter<double>("R_jet_def_pruning"));
     Pruner pruner(jet_def_Pruning,(*itPruned).getParameter<double>("z_cut"), (*itPruned).getParameter<double>("R_Cut"));
     PseudoJet jetTemp = pruner(iJet) ;
     lPruned.push_back(jetTemp);
     lPrunedSafe.push_back((*area_subtractor)(jetTemp));
    }
  }

  // soft drop information
  vector<PseudoJet> lSoftDropped ;
  vector<PseudoJet> lSoftDroppedSafe ;
  double SoftDropedSymmetry = -1.0;
  double SoftDropedDR = -1.0;
  double SoftDropedMassDrop = -1.0;
  double SoftDropedEnergyLoss = -1.0;
  double SoftDropedArea = -1.0;
  double SoftDropedNconst = -1.0;
  PseudoJet filtered_softdropped_jet;

  if (iJet.pt() > jetPtTresholdForGroomers){
    std::vector<edm::ParameterSet>::const_iterator itSoft = softDropParam.begin();
    for( ; itSoft != softDropParam.end() ; ++itSoft){
     contrib::SoftDrop softdrop((*itSoft).getParameter<double>("beta"),(*itSoft).getParameter<double>("symmetry_cut"), (*itSoft).getParameter<double>("R0"));
     lSoftDropped.push_back(softdrop(iJet));   
     softdrop.set_subtractor(area_subtractor);
     lSoftDroppedSafe.push_back(softdrop(iJet));
    }  

    if (lSoftDroppedSafe.at(0)!=0 and lSoftDroppedSafe.at(0).m()>0.0){
      SoftDropedSymmetry = lSoftDroppedSafe.at(0).structure_of<contrib::SoftDrop>().symmetry();
      SoftDropedDR = lSoftDroppedSafe.at(0).structure_of<contrib::SoftDrop>().delta_R();
      SoftDropedMassDrop = lSoftDroppedSafe.at(0).structure_of<contrib::SoftDrop>().mu();
      SoftDropedEnergyLoss = 1-lSoftDroppedSafe.at(0).pt()/iJet.pt();
      SoftDropedArea = lSoftDroppedSafe.at(0) .area() ;
      SoftDropedNconst = lSoftDroppedSafe.at(0) .constituents().size() ;

      // filter jet dynamically based on deltaR between subjets (arXiv:0802.2470)
      double dyn_Rfilt = min(0.3, SoftDropedDR*0.5);
      int dyn_nfilt = 3;
      Filter filtersoft(dyn_Rfilt, SelectorNHardest(dyn_nfilt));
      filtered_softdropped_jet = filtersoft(lSoftDroppedSafe.at(0));
    }
  }

  // -- apply the JEC --> function from JetTools
  double lJEC     = correction(iJet,iJetCorr,bge_rho.rho());  
  double lUnc     = unc       (iJet,iJetUnc);

  if(isPuppi){
    lJEC = correction(iJet,iJetCorr,1.);
    lUnc = unc       (iJet,iJetUnc);
  }

  int iId = 0; 
  if(isCHS)   iId = 1; 
  if(isPuppi) iId = 2;
  double lJECPhil = correctPhil((iJet.pt())*lJEC,iJet.eta(),iId,iCorr); // take phill JEC 
 
  // -- Top Taggers 
  fastjet::PseudoJet iJetCA;

  double hepttJetMass    = -1; 
  double hepttWMass      = -1; 
  double hepttM01        = -1; 
  double hepttM02        = -1; 
  double hepttM12        = -1; 
  double hepttM12M012    = -1; 
  double hepttAtanM02M01 = -1; 
  double cmsttJetMass     = -1;
  double cmsttMinMass     = -1;
  double cmsttHelicity    = -1;
  double cmsttNsubjets    = -1;


  if (iJet.pt()> jetPtTresholdForTopTagging){

    // -- recluster jet CA
    JetDefinition jet_def_CA (fastjet::cambridge_algorithm, jetR*10); //large R to cluster all constituents of original jet
    fastjet::ClusterSequence cs_Recluster (jetParticles,jet_def_CA);
    vector<fastjet::PseudoJet> jets_Recluster = sorted_by_pt(cs_Recluster.inclusive_jets());
    iJetCA = jets_Recluster[0];

    // -- HEP Top Tagger 
    double mass_drop_threshold = 0.8;
    double max_subjet_mass     = 30;
    bool use_subjet_mass_cuts  = false;
    HEPTopTagger hep_top_tagger(mass_drop_threshold, max_subjet_mass, use_subjet_mass_cuts);
    
    PseudoJet hep_top_candidate   = hep_top_tagger( iJetCA );

    if (hep_top_candidate != 0){
      PseudoJet W =     hep_top_candidate.structure_of<HEPTopTagger>().W();
      PseudoJet W1 =    hep_top_candidate.structure_of<HEPTopTagger>().W1();
      PseudoJet W2 =    hep_top_candidate.structure_of<HEPTopTagger>().W2();
      PseudoJet non_W = hep_top_candidate.structure_of<HEPTopTagger>().non_W();

      vector<PseudoJet> all_subjets;
      all_subjets.push_back(W1);
      all_subjets.push_back(W2);
      all_subjets.push_back(non_W);
      all_subjets = sorted_by_pt(all_subjets);

      PseudoJet sum012 = all_subjets[0]+all_subjets[1]+all_subjets[2];
      PseudoJet sum01 = all_subjets[0]+all_subjets[1];
      PseudoJet sum02 = all_subjets[0]+all_subjets[2];
      PseudoJet sum12 = all_subjets[1]+all_subjets[2];

      hepttJetMass       = hep_top_candidate.m();
      hepttWMass         = W.m();
      hepttM01           = sum01.m();
      hepttM02           = sum02.m();
      hepttM12           = sum12.m();
      if ( sum012.m()!=0 ) hepttM12M012     = sum12.m() / sum012.m() ;
      if ( sum01.m()!=0 )  hepttAtanM02M01  = atan( sum02.m() / sum01.m() ) ;
    }
  
    // -- CMS Top Tagger 
    double cms_delta_p = 0.05;
    double cms_delta_r = 0.4;
    double A           = 0.0004;

    CMSTopTagger cms_top_tagger(cms_delta_p, cms_delta_r, A);
    PseudoJet cms_top_candidate  = cms_top_tagger( iJetCA );

    if (cms_top_candidate != 0){
      vector<PseudoJet> kept_subjets0 = cms_top_candidate.structure_of<CMSTopTagger>().W().pieces();
      vector<PseudoJet> kept_subjets1 = cms_top_candidate.structure_of<CMSTopTagger>().non_W().pieces();
      vector<PseudoJet> all_subjets = kept_subjets0;
      all_subjets.insert( all_subjets.end(), kept_subjets1.begin(), kept_subjets1.end() );

      cmsttJetMass      = cms_top_candidate.m();
      cmsttMinMass      = cms_top_candidate.structure_of<CMSTopTagger>().W().m();
      cmsttHelicity     = cms_top_candidate.structure_of<CMSTopTagger>().cos_theta_W();
      cmsttNsubjets     = all_subjets.size();
    } 
  }

  // -- find the gen jet matched to this reco jet
  int imatch   = -1;
  bool matched = false;
  if (isMC){
    imatch  = matchingIndexFromJetInfo(iJet,iGenJetI); // matching to Genjets
    matched = IsMatchedToGenBoson( eta_Boson, phi_Boson, iJet); // matching to Vector boson at generator level
  }

  float lPtPhil = iJet.pt()*lJEC*float(lJECPhil); // Phill pt corrected
  
  // -- Fil Jet Info
  (iJetI.pt        ).push_back(lCorr     .pt());        // safe 4V jet pT
  (iJetI.ptcorr    ).push_back(iJet      .pt()*lJEC);   // raw pt * JEC (standard)
  (iJetI.ptcorrphil).push_back(lPtPhil);                // raw pt * Phill correction
  (iJetI.ptraw     ).push_back(iJet      .pt());        // raw pt

  (iJetI.eta       ).push_back(iJet      .eta());
  (iJetI.phi       ).push_back(iJet      .phi());
  (iJetI.mraw      ).push_back(iJet      .m());
  (iJetI.m         ).push_back(lCorr     .m());
  (iJetI.ptunc     ).push_back(lUnc);

  (iJetI.ptconst   ).push_back(lConstit  .pt());
  (iJetI.mconst    ).push_back(lConstit  .m());
    
  if(iJet.pt() > jetPtTresholdForGroomers){ // trimming info
   for( unsigned int iTrim = 0 ; iTrim < lTrim.size() ; iTrim++){
     iJetI.pttrim.at(iTrim).push_back(lTrim.at(iTrim).pt());
     iJetI.mtrim.at(iTrim).push_back(lTrim.at(iTrim).m());
     iJetI.pttrimsafe.at(iTrim).push_back(lTrimSafe.at(iTrim).pt());
     iJetI.mtrimsafe.at(iTrim).push_back(lTrimSafe.at(iTrim).m());
   }
  }
  else{
   for( unsigned int iTrim = 0 ; iTrim < trimmingParam.size() ; iTrim++){
     iJetI.pttrim.at(iTrim).push_back(-999.);
     iJetI.mtrim.at(iTrim).push_back(-999.);
     iJetI.pttrimsafe.at(iTrim).push_back(-999.);
     iJetI.mtrimsafe.at(iTrim).push_back(-999.);
   }
  }

  if(iJet.pt() > jetPtTresholdForGroomers){ // pruning info
   for( unsigned int iPruned = 0 ; iPruned < lPruned.size() ; iPruned++){
     iJetI.ptpruned.at(iPruned).push_back(lPruned.at(iPruned).pt());
     iJetI.mpruned.at(iPruned).push_back(lPruned.at(iPruned).m());
     iJetI.ptprunedsafe.at(iPruned).push_back(lPrunedSafe.at(iPruned).pt());
     iJetI.mprunedsafe.at(iPruned).push_back(lPrunedSafe.at(iPruned).m());
   }
  }
  else{
   for( unsigned int iPruned = 0 ; iPruned < pruningParam.size() ; iPruned++){
     iJetI.ptpruned.at(iPruned).push_back(-999.);
     iJetI.mpruned.at(iPruned).push_back(-999.);
     iJetI.ptprunedsafe.at(iPruned).push_back(-999.);
     iJetI.mprunedsafe.at(iPruned).push_back(-999.);
   }
  }

  if(iJet.pt() > jetPtTresholdForGroomers){ // soft drop info
   for( unsigned int iSoft = 0 ; iSoft < lSoftDropped.size() ; iSoft++){
     iJetI.ptsoftdrop.at(iSoft).push_back(lSoftDropped.at(iSoft).pt());
     iJetI.msoftdrop.at(iSoft).push_back(lSoftDropped.at(iSoft).m());
     iJetI.ptsoftdropsafe.at(iSoft).push_back(lSoftDroppedSafe.at(iSoft).pt());
     iJetI.msoftdropsafe.at(iSoft).push_back(lSoftDroppedSafe.at(iSoft).m());
   }
  }
  else{
   for( unsigned int iSoft = 0 ; iSoft < pruningParam.size() ; iSoft++){
     iJetI.ptsoftdrop.at(iSoft).push_back(-999.);
     iJetI.msoftdrop.at(iSoft).push_back(-999.);
     iJetI.ptsoftdropsafe.at(iSoft).push_back(-999.);
     iJetI.msoftdropsafe.at(iSoft).push_back(-999.);
   }
  }

  // other varuables
  (iJetI.sdsymmetry ).push_back( SoftDropedSymmetry );
  (iJetI.sddeltar ).push_back( SoftDropedDR );
  (iJetI.sdmu ).push_back( SoftDropedMassDrop );
  (iJetI.sdenergyloss ).push_back( SoftDropedEnergyLoss );
  (iJetI.sdarea ).push_back( SoftDropedArea );
  (iJetI.sdnconst ).push_back( SoftDropedNconst );
  (iJetI.mfiltsoftdrop ).push_back( filtered_softdropped_jet.m() );

  (iJetI.nparticles).push_back(jetParticles.size());
  (iJetI.nneutrals ).push_back(neutrals.size());
  (iJetI.ncharged  ).push_back(chargedLV.size()+chargedPU.size());

  (iJetI.is_MatchedToBoson ).push_back(matched); //boson matching information

  (iJetI.hepmass ).push_back( hepttJetMass );
  (iJetI.hepwmass ).push_back( hepttWMass );
  (iJetI.hepm01 ).push_back( hepttM01 );
  (iJetI.hepm02 ).push_back( hepttM02 );
  (iJetI.hepm12 ).push_back( hepttM12 );
  (iJetI.hepm12m012 ).push_back( hepttM12M012 );
  (iJetI.hepatanm02m01).push_back( hepttAtanM02M01 );

  (iJetI.cmsmass ).push_back(cmsttJetMass );
  (iJetI.cmsminmass ).push_back(cmsttMinMass );
  (iJetI.cmshelicity ).push_back(cmsttHelicity );
  (iJetI.cmsnsubjets ).push_back(cmsttNsubjets );
 
  // V-tagging variables 
  if (iJet.pt() > jetPtTresholdForGroomers){

    vtagger->setInputJet(iJet); 
    for(unsigned int iNsubjettiness = 0; iNsubjettiness < NsubjettinessParam.size() ; iNsubjettiness++){ // N-subjettiness
      (iJetI.tau1 ).at(iNsubjettiness).push_back(vtagger->computeNSubJettines(1,NsubjettinessParam.at(iNsubjettiness).getParameter<double>("beta"),NsubjettinessParam.at(iNsubjettiness).getParameter<double>("R0"),NsubjettinessParam.at(iNsubjettiness).getParameter<double>("Rcut")));
      (iJetI.tau2 ).at(iNsubjettiness).push_back(vtagger->computeNSubJettines(2,NsubjettinessParam.at(iNsubjettiness).getParameter<double>("beta"),NsubjettinessParam.at(iNsubjettiness).getParameter<double>("R0"),NsubjettinessParam.at(iNsubjettiness).getParameter<double>("Rcut")));
      (iJetI.tau3 ).at(iNsubjettiness).push_back(vtagger->computeNSubJettines(3,NsubjettinessParam.at(iNsubjettiness).getParameter<double>("beta"),NsubjettinessParam.at(iNsubjettiness).getParameter<double>("R0"),NsubjettinessParam.at(iNsubjettiness).getParameter<double>("Rcut")));
    }

    (iJetI.Qjets).push_back(vtagger->computeQjets(100,25,randNumber.Uniform(0,10000),jetR)); // Qjets
    
    for( unsigned int iECF = 0; iECF < ecfParam.size() ; iECF++) // ECF
       iJetI.ecf.at(iECF).push_back(vtagger->computeECF(get_algo(ecfParam.at(iECF).getParameter<string>("ecfAlgo")),ecfParam.at(iECF).getParameter<double>("Rparam"),ecfParam.at(iECF).getParameter<int>("nPoint"),ecfParam.at(iECF).getParameter<double>("beta"),ecfParam.at(iECF).getParameter<int>("type")));

    for( unsigned int iCharge = 0; iCharge < chargeParam.size() ; iCharge++) // Charge
      iJetI.charge.at(iCharge).push_back(vtagger->computeJetChargeReco(chargeParam[iCharge]));
    
    vector<PseudoJet> subjets_pruned ; // pruned subjet
    for( unsigned int iPrun = 0 ; iPrun < lPruned.size() ; iPrun++){
      if(lPruned.at(iPrun).constituents().size() > 1){
       subjets_pruned = lPruned.at(iPrun).associated_cluster_sequence()->exclusive_subjets(lPruned.at(iPrun),2); // take the exclusive subjets
       subjets_pruned = sorted_by_pt(subjets_pruned);
       if(subjets_pruned.at(0).pt()>0){
        iJetI.pullangle.at(iPrun).push_back(vtagger->computePullAngle(subjets_pruned,jetR));  
        iJetI.pullmagnitude.at(iPrun).push_back(TMath::Sqrt(iJetI.pullangle.at(iPrun).back()*iJetI.pullangle.at(iPrun).back()+subjets_pruned.at(0).rapidity()*subjets_pruned.at(0).rapidity()));  
       }
       else{
	iJetI.pullangle.at(iPrun).push_back(999);
	iJetI.pullmagnitude.at(iPrun).push_back(999);
       }						       
      }
      else{
	iJetI.pullangle.at(iPrun).push_back(999);
	iJetI.pullmagnitude.at(iPrun).push_back(999);
      }
    }

    vtagger->setInputJet(lPruned.at(0)); // take the leading parameter combination pruned jet

    if(lPruned.at(0).pt() > 0){
     (iJetI.tau1_pr ).push_back(vtagger->computeNSubJettines(1,1.,jetR,jetR));
     (iJetI.tau2_pr ).push_back(vtagger->computeNSubJettines(2,1.,jetR,jetR));
     (iJetI.tau3_pr ).push_back(vtagger->computeNSubJettines(3,1.,jetR,jetR));
    }   
    else{
     (iJetI.tau1_pr ).push_back(999);
     (iJetI.tau2_pr ).push_back(999);
     (iJetI.tau3_pr ).push_back(999);
    }

    // Calculate QGL for pruned jets 
    for( unsigned int iPrun = 0 ; iPrun < lPruned.size() ; iPrun++){
      if(lPruned.at(iPrun).pt() > 0 ){
	if(isCHS) iJetI.QGLikelihood_pr.at(iPrun).push_back(vtagger->computeQGLikelihood(qgLikelihoodCHS,lJEC,rho));
	else iJetI.QGLikelihood_pr.at(iPrun).push_back(vtagger->computeQGLikelihood(qgLikelihood,lJEC,rho));
      }
      else iJetI.QGLikelihood_pr.at(iPrun).push_back(999);      
    }

    subjets_pruned.clear() ;
    for( unsigned int iPrun = 0 ; iPrun < lPruned.size() ; iPrun++){
      if(lPruned.at(iPrun).constituents().size() > 1){
       subjets_pruned = lPruned.at(iPrun).associated_cluster_sequence()->exclusive_subjets(lPruned.at(iPrun),2);
       subjets_pruned = sorted_by_pt(subjets_pruned);
       if(subjets_pruned.at(0).pt() > 0){
        vtagger->setInputJet(subjets_pruned.at(0));   
        if(isCHS) iJetI.QGLikelihood_pr_sub1.at(iPrun).push_back(vtagger->computeQGLikelihood(qgLikelihoodCHS,lJEC,rho));
        else iJetI.QGLikelihood_pr_sub1.at(iPrun).push_back(vtagger->computeQGLikelihood(qgLikelihood,lJEC,rho));
       }
       else iJetI.QGLikelihood_pr_sub1.at(iPrun).push_back(999);     

       if(subjets_pruned.at(1).pt() > 0){
        vtagger->setInputJet(subjets_pruned.at(1));   
        if(isCHS) iJetI.QGLikelihood_pr_sub2.at(iPrun).push_back(vtagger->computeQGLikelihood(qgLikelihoodCHS,lJEC,rho));
        else iJetI.QGLikelihood_pr_sub2.at(iPrun).push_back(vtagger->computeQGLikelihood(qgLikelihood,lJEC,rho));
       }
       else iJetI.QGLikelihood_pr_sub2.at(iPrun).push_back(999);     
     }
      else{
	iJetI.QGLikelihood_pr_sub1.at(iPrun).push_back(999);
	iJetI.QGLikelihood_pr_sub2.at(iPrun).push_back(999);
      }
    }

    vtagger->setInputJet(lSoftDropped.at(0)); 
    if(lSoftDropped.at(0).pt() > 0 ){
     (iJetI.tau1_softdrop ).push_back(vtagger->computeNSubJettines(1,1.,jetR,jetR));
     (iJetI.tau2_softdrop ).push_back(vtagger->computeNSubJettines(2,1.,jetR,jetR));
     (iJetI.tau3_softdrop ).push_back(vtagger->computeNSubJettines(3,1.,jetR,jetR));
    }
    else{
     (iJetI.tau1_softdrop ).push_back(999);
     (iJetI.tau2_softdrop ).push_back(999);
     (iJetI.tau3_softdrop ).push_back(999);
    }
  }
  else{

    for( unsigned int iCharge = 0; iCharge < chargeParam.size() ; iCharge++)
      iJetI.charge.at(iCharge).push_back(999.);

    for( unsigned int iECF = 0; iECF < ecfParam.size() ; iECF++)
      iJetI.ecf.at(iECF).push_back(999.);

    for( unsigned int iPrun = 0 ; iPrun < lPruned.size() ; iPrun++){
      iJetI.QGLikelihood_pr.at(iPrun).push_back(999.);
      iJetI.QGLikelihood_pr_sub1.at(iPrun).push_back(999.);
      iJetI.QGLikelihood_pr_sub2.at(iPrun).push_back(999.);
      iJetI.pullangle.at(iPrun).push_back(999);
      iJetI.pullmagnitude.at(iPrun).push_back(999);
    }

    for(unsigned int iNsubjettiness = 0; iNsubjettiness < NsubjettinessParam.size() ; iNsubjettiness++){
      (iJetI.tau1 ).at(iNsubjettiness).push_back(999.);
      (iJetI.tau2 ).at(iNsubjettiness).push_back(999.);
      (iJetI.tau3 ).at(iNsubjettiness).push_back(999.);
    }

    (iJetI.Qjets).push_back(999.);
    (iJetI.tau1_pr ).push_back(999.);
    (iJetI.tau2_pr ).push_back(999.);
    (iJetI.tau3_pr ).push_back(999.);

    (iJetI.tau1_softdrop ).push_back(999.);
    (iJetI.tau2_softdrop ).push_back(999.);
    (iJetI.tau3_softdrop ).push_back(999.);

  }

  // store reco-gen matching informations
  if (imatch > -1){

    (iJetI.imatch).push_back(imatch); // store the

    (iJetI.ptgen    ).push_back((iGenJetI.pt)[imatch]);
    (iJetI.etagen   ).push_back((iGenJetI.eta)[imatch]);
    (iJetI.phigen   ).push_back((iGenJetI.phi)[imatch]);
    (iJetI.mgen     ).push_back((iGenJetI.m)[imatch]);
    (iJetI.mrawgen  ).push_back((iGenJetI.mraw)[imatch]);

    if(int(iGenJetI.mtrim.at(0).size()) >= imatch) (iJetI.mtrimgen).push_back((iGenJetI.mtrim.at(0))[imatch]);
    else (iJetI.mtrimgen    ).push_back(-999.);

    if(int(iGenJetI.mtrimsafe.at(0).size())>=imatch) (iJetI.mtrimsafegen).push_back((iGenJetI.mtrimsafe.at(0))[imatch]);
    else (iJetI.mtrimsafegen).push_back(-999.);

    if(int(iGenJetI.msoftdrop.at(0).size()) >= imatch) (iJetI.msoftdropgen).push_back((iGenJetI.msoftdrop.at(0))[imatch]);
    else (iJetI.msoftdropgen).push_back(-999.); 

    if(int(iGenJetI.msoftdropsafe.at(0).size())>= imatch) (iJetI.msoftdropsafegen    ).push_back((iGenJetI.msoftdropsafe.at(0))[imatch]);
    else (iJetI.msoftdropsafegen    ).push_back(-999.);

    if(computeJetFlavour) (iJetI.flavourgen).push_back((iGenJetI.jetflavour)[imatch]);
    else (iJetI.flavourgen).push_back( -999.);
  }
  else { 
    (iJetI.imatch).push_back(imatch);

    (iJetI.ptgen    ).push_back(-999.);
    (iJetI.etagen   ).push_back(-999.);
    (iJetI.phigen   ).push_back(-999.);
    (iJetI.mgen     ).push_back(-999.);
    (iJetI.mrawgen     ).push_back(-999.);
    (iJetI.mtrimgen    ).push_back(-999.);
    (iJetI.mtrimsafegen).push_back(-999.);

    (iJetI.msoftdropgen        ).push_back( -999.);
    (iJetI.msoftdropsafegen    ).push_back( -999.);

    (iJetI.flavourgen          ).push_back( -999.);

  }

}


// ------------------------------------------------------------------------------------------
void fillGenJetsInfo(vector<PseudoJet> &iJets, // set of GenJets in the event
                     vector<PseudoJet> &iParticles, // set of whole gen particles
                     GenJetInfo &iJetInfo, 
                     vector<JetCleanser> &cleanser_vect, // cleansing vector
                     int nPU, int nPV, double rho) { 


  // -- Compute rho, rho_m for SafeAreaSubtraction
  AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(SelectorAbsRapMax(5.0))); // area definition for jet with explicit ghosts
  JetDefinition jet_def_for_rho(kt_algorithm, 0.4); // jet clustering for background estimation
  Selector rho_range =  SelectorAbsRapMax(5.0); // apply just an eta acceptance cut
  ClusterSequenceArea clust_seq_rho(iParticles, jet_def_for_rho, area_def); // cluster the initial particles with the k_t + active ghost for bkg determination
  JetMedianBackgroundEstimator  bge_rho(rho_range, clust_seq_rho);
  JetMedianBackgroundEstimator  bge_rhom(rho_range, clust_seq_rho); // get mediam background estimator
  BackgroundJetPtMDensity m_density;
  bge_rhom.set_jet_density_class(&m_density);
    
  // -- Background estimator for constituents subtractor
  JetMedianBackgroundEstimator bge_rhoC(rho_range,jet_def_for_rho, area_def);
  BackgroundJetScalarPtDensity *scalarPtDensity = new BackgroundJetScalarPtDensity();
  bge_rhoC.set_jet_density_class(scalarPtDensity);
  bge_rhoC.set_particles(iParticles);
    
  // -- Clear jet info for each event                                                                                                                                     
  clear(iJetInfo);  
  iJetInfo.npu = nPU;
  iJetInfo.npv = nPV;

  // -- Loop over jets in the event and set jets variables                                                                                                           
  for (unsigned int j = 0; j < iJets.size(); j++){
     if(j == 0) 
       setGenJet( iJets[j], iJetInfo,  bge_rho, bge_rhom, bge_rhoC, cleanser_vect, 1, rho); // give the original clustered jets, the background estimations and cleansing
     else 
       setGenJet( iJets[j], iJetInfo,  bge_rho, bge_rhom, bge_rhoC, cleanser_vect, 0, rho); // give the original clustered jets, the background estimations and cleansing
  }
}



// ------------------------------------------------------------------------------------------
void fillRecoJetsInfo(vector<PseudoJet> &iJets,  
                      vector<PseudoJet> &iParticles, 
                      vector<PseudoJet> &iAllParticles, 
                      JetInfo &iJetInfo, 
                      GenJetInfo iGenJetInfo, 
                      bool isCHS, 
                      FactorizedJetCorrector *jetCorr, JetCorrectionUncertainty *ijetUnc,
                      std::vector<TGraph*> &iCorr, 
                      vector<JetCleanser> &cleanser_vect, 
                      int nPU, int nPV, double rho, vfloat eta_Boson, vfloat phi_Boson, const bool & isPuppi = false, bool isMC=true){

  
  // -- Compute rho, rho_m for SafeAreaSubtraction -> same procedure is used for GenJets
  AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(SelectorAbsRapMax(5.0))); // jet area definition with explicit ghosts
  JetDefinition jet_def_for_rho(kt_algorithm, 0.4); // jet definition for bkg estimation
  Selector rho_range =  SelectorAbsRapMax(5.0);
  ClusterSequenceArea clust_seq_rho(iAllParticles, jet_def_for_rho, area_def); // cluster sequence
  JetMedianBackgroundEstimator bge_rho(rho_range, clust_seq_rho);
  JetMedianBackgroundEstimator bge_rhom(rho_range, clust_seq_rho);
  BackgroundJetPtMDensity m_density;
  bge_rhom.set_jet_density_class(&m_density);
  
  // -- Background estimator for constituents subtractor
  JetMedianBackgroundEstimator bge_rhoC(rho_range,jet_def_for_rho, area_def);
  BackgroundJetScalarPtDensity *scalarPtDensity = new BackgroundJetScalarPtDensity();
  bge_rhoC.set_jet_density_class(scalarPtDensity);
  bge_rhoC.set_particles(iParticles);

  // -- Compute rho, rho_m for SafeAreaSubtraction -> same procedure is used for GenJets
  AreaDefinition area_def_chs(active_area_explicit_ghosts,GhostedAreaSpec(SelectorAbsRapMax(2.5)));
  JetDefinition jet_def_for_rho_chs(kt_algorithm, 0.4);
  Selector rho_range_chs =  SelectorAbsRapMax(2.5);  
  ClusterSequenceArea clust_seq_rho_chs(iParticles, jet_def_for_rho_chs, area_def_chs); // only for chs particles
  JetMedianBackgroundEstimator  bge_rho_chs  (rho_range_chs, clust_seq_rho_chs);
  JetMedianBackgroundEstimator  bge_rhom_chs (rho_range_chs, clust_seq_rho_chs);
  BackgroundJetPtMDensity m_density_chs;
  bge_rhom_chs.set_jet_density_class(&m_density_chs);
  // -- Background estimator for constituents subtractor
  JetMedianBackgroundEstimator bge_rhoC_chs(rho_range_chs,jet_def_for_rho_chs, area_def_chs);
  BackgroundJetScalarPtDensity *scalarPtDensity_chs = new BackgroundJetScalarPtDensity();
  bge_rhoC_chs.set_jet_density_class(scalarPtDensity_chs);
  bge_rhoC_chs.set_particles(iParticles);

  // -- Clear jet info for each event                                                                                                                                           
  clear(iJetInfo);  
  iJetInfo.npu = nPU;
  iJetInfo.npv = nPV;

  // -- Loop over jets in the event and set jets variables                                                                                                                      
  for (unsigned int j = 0; j < iJets.size(); j++){
    if(fabs(iJets[j].eta()) < 2.5 && isCHS) { // is CHS is applied and jet is inside tracker acceptance used bkg estimated after chs 
      setRecoJet( iJets[j], iJetInfo, iGenJetInfo, bge_rho_chs, bge_rhom_chs, bge_rhoC_chs, isCHS, jetCorr, ijetUnc, iCorr, cleanser_vect, (j==0), rho, eta_Boson, phi_Boson, isPuppi, isMC);
    } else { 
      setRecoJet( iJets[j], iJetInfo, iGenJetInfo, bge_rho, bge_rhom, bge_rhoC, isCHS, jetCorr, ijetUnc ,iCorr, cleanser_vect,(j==0), rho, eta_Boson, phi_Boson, isPuppi, isMC);
    }
  }
}


// -------- Method to setbranch address for TJet in the input file, which have been produced (clustered) inside CMSSW
void setupCMSSWJetReadOut(TTree *iTree, float R ) {  
  cout << "Setting up to read jet collection : " << Form("Jet0%d",int(R*10)) << endl;
  fJet  = new TClonesArray("baconhep::TJet");
  iTree->SetBranchAddress(Form("Jet0%d",int(R*10)), &fJet);
}


// ------------------------------------------------------------------------------------------
void readCMSSWJet(int entry, TTree *iTree, TTree &oTree,  std::vector<fastjet::PseudoJet> genJets, JetInfo &iJetI) {

  // -- Clear jet info for each event
  clear(iJetI);
  // -- Read event and fill jet info
  iTree->GetEntry(entry);

  for (int i = 0; i < fJet->GetEntriesFast(); i++){
    TJet *pJet = (TJet*)((*fJet)[i]);
    
    // -- fill jet info                                                                                                                                                       
    (iJetI.pt        ).push_back(pJet->pt);
    (iJetI.ptcorr    ).push_back(pJet->pt);
    (iJetI.ptcorrphil).push_back(pJet->pt);
    (iJetI.ptraw     ).push_back(pJet->ptRaw);
    (iJetI.eta       ).push_back(pJet->eta);
    (iJetI.phi       ).push_back(pJet->phi);
    (iJetI.m         ).push_back(pJet->mass);
    (iJetI.nparticles).push_back(pJet->nParticles);
    (iJetI.nneutrals ).push_back(pJet->nNeutrals);
    (iJetI.ncharged  ).push_back(pJet->nCharged);
    
    // for now fill this branches with dummy value
    (iJetI.mraw      ).push_back(-999.);
    (iJetI.ptunc     ).push_back(-999.);
    (iJetI.ptclean   ).push_back(-999.);
    (iJetI.mclean   ).push_back(-999.);

    (iJetI.ptconst   ).push_back(-999.);
    (iJetI.mconst    ).push_back(-999.);
    
    for( unsigned int iTrim = 0 ; iTrim != trimmingParam.size() ; iTrim++){
     iJetI.pttrim.at(iTrim).push_back(-999.);
     iJetI.mtrim.at(iTrim).push_back(-999.);
     iJetI.pttrimsafe.at(iTrim).push_back(-999.);
     iJetI.mtrimsafe.at(iTrim).push_back(-999.);
    }

    for( unsigned int iPruned = 0 ; iPruned != pruningParam.size() ; iPruned++){
      (iJetI.ptpruned    ).at(iPruned).push_back(-999.);
      (iJetI.mpruned     ).at(iPruned).push_back(-999.);
      (iJetI.ptprunedsafe).at(iPruned).push_back(-999.);
      (iJetI.mprunedsafe).at(iPruned).push_back(-999.);
    }

    for( unsigned int iSoft = 0 ; iSoft != softDropParam.size() ; iSoft++){
      (iJetI.ptsoftdrop).at(iSoft).push_back(-999.);
      (iJetI.msoftdrop).at(iSoft).push_back(-999.);
      (iJetI.ptsoftdropsafe).at(iSoft).push_back(-999.);
      (iJetI.msoftdropsafe).at(iSoft).push_back(-999.);
    }
    
    //-- gen matching
    int imatch = -1;
    float mindr = dRMatching;
    TLorentzVector *recojet = new TLorentzVector();
    recojet->SetPtEtaPhiM(pJet->pt, pJet->eta, pJet->phi, pJet->mass);
    for (unsigned int ig = 0; ig < genJets.size(); ig++){
      TLorentzVector *genjet = new TLorentzVector();
      genjet->SetPtEtaPhiM(genJets[ig].pt(), genJets[ig].eta(), genJets[ig].phi(), genJets[ig].m());
      double dr = recojet->DeltaR(*genjet);
      if (dr < mindr){
	mindr = dr;
	imatch = ig;
      }
      delete genjet;
    }
    
    delete recojet;

    if (imatch > -1){
      (iJetI.imatch   ).push_back(imatch);
      (iJetI.ptgen    ).push_back(genJets[imatch].pt());
      (iJetI.etagen   ).push_back(genJets[imatch].eta());
      (iJetI.phigen   ).push_back(genJets[imatch].phi());
      (iJetI.mgen     ).push_back(genJets[imatch].m());
      (iJetI.mrawgen     ).push_back(-999.);// dummy val
      (iJetI.mtrimgen    ).push_back(-999.);// dummy val
      (iJetI.mtrimsafegen).push_back(-999.);// dummy val
    }
    else {
      (iJetI.imatch   ).push_back(imatch);
      (iJetI.ptgen    ).push_back(-999.);
      (iJetI.etagen   ).push_back(-999.);
      (iJetI.phigen   ).push_back(-999.);
      (iJetI.mgen     ).push_back(-999.);
      (iJetI.mrawgen     ).push_back(-999.);// dummy val
      (iJetI.mtrimgen    ).push_back(-999.);// dummy val
      (iJetI.mtrimsafegen).push_back(-999.);// dummy val
    }
  }
  // --- fill tree 
  oTree.Fill();  
}



// function to fill a TChain with the list of input files to be processed 
bool FillChain(TChain& chain, const std::string& inputFileList){

  std::ifstream inFile(inputFileList.c_str());
  std::string buffer;

  if(!inFile.is_open()){
      std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
      return false;
  }
  
  while(1){
      inFile >> buffer;
      if(!inFile.good()) break;
      chain.Add(buffer.c_str());
  }

  return true;
}
 

//---------------------------------------------------------------------------------------------------------------
//--- MAIN PROGRAM
//---------------------------------------------------------------------------------------------------------------
int main (int argc, char ** argv) {
    
  // --- args
  if (argc<3){
    cout << "Missing arguments!!!" <<endl;
    cout << "Usage: MiniNtuplizer <config> <input files list> <output file>" <<endl;
  }

  // args 
  std::string inputFilesList = argv[2]; // input file list name
  std::string fOut           = argv[3]; // output file name

  // --- Read configurable parameters from config                                                                                                                        
  std::string configFileName = argv[1];
  boost::shared_ptr<edm::ParameterSet> parameterSet = edm::readConfig(configFileName);  
  edm::ParameterSet Options  = parameterSet -> getParameter<edm::ParameterSet>("Options");
 
  //Global event information  
  bool isMC                  = Options.getParameter<bool>("isMC");           // MC or data
  int maxEvents              = Options.getParameter<int>("maxEvents");       // max num of events to analyze
  int minEvents              = Options.getParameter<int>("minEvents");       // min num of events to analyze -> starting point                                                 
  double jetPtCut            = Options.getParameter<double>("jetPtCut");     //pT cut applied when getting jets from cluster sequence 
  jetR                       = Options.getParameter<double>("jetR");         // jet cone size  
  std::string jetAlgo        = Options.getParameter<std::string>("jetAlgo"); // jet clustering algorithm --> default clustering

  fastjet::JetAlgorithm fatjet_algo = get_algo("jetAlgo"); // take the fastjet definition

  bool doCMSSWJets           = Options.getParameter<bool>("doCMSSWJets");        // analyze also default CMSSW PF jets
  bool doSoftKillerJets      = Options.getParameter<bool>("doSoftKillerJets");   // run soft killer jets
  std::string puppiConfig    = Options.getParameter<std::string>("puppiConfig"); // Puppi congiguration file

  // thresholds for Top and groomed jets 
  jetPtTresholdForGroomers      = Options.getParameter<double>("jetPtTresholdForGroomers");  
  genjetPtTresholdForGroomers   = Options.getParameter<double>("genjetPtTresholdForGroomers");
  jetPtTresholdForTopTagging    = Options.getParameter<double>("jetPtTresholdForTopTagging");
  genJetPtTresholdForTopTagging = Options.getParameter<double>("genJetPtTresholdForTopTagging");

  // JEC set
  std::string L1FastJetJEC    = Options.getParameter<std::string>("L1FastJetJEC");     // L1 JEC 
  std::string L2RelativeJEC   = Options.getParameter<std::string>("L2RelativeJEC");    // L2
  std::string L3AbsoluteJEC   = Options.getParameter<std::string>("L3AbsoluteJEC");    // L3
  std::string L2L3ResidualJEC = Options.getParameter<std::string>("L2L3ResidualJEC");  // L2L3 residual (for data only)
  std::string JECUncertainty  = Options.getParameter<std::string>("JECUncertainty");   // Uncertainty

  std::string L1FastJetJEC_CHS    = Options.getParameter<std::string>("L1FastJetJEC_CHS");    // L1 JEC CHS 
  std::string L2RelativeJEC_CHS   = Options.getParameter<std::string>("L2RelativeJEC_CHS");   // L2 CHS
  std::string L3AbsoluteJEC_CHS   = Options.getParameter<std::string>("L3AbsoluteJEC_CHS");   // L3 CHS
  std::string L2L3ResidualJEC_CHS = Options.getParameter<std::string>("L2L3ResidualJEC_CHS"); // L2L3 residual (for data only) CHS
  std::string JECUncertainty_CHS  = Options.getParameter<std::string>("JECUncertainty_CHS");  // Uncertainty CHS

  std::string JECPhil             = Options.getParameter<std::string>("PhilJEC"); // Phill correction file .root

  // Quark Gluon Likelihood
  QGinputWeightFilePath     = Options.getParameter<std::string>("QGinputWeightFilePath"); // path where find QGL information

  // matching with the truth
  bool DoMatchingToBoson      = Options.getParameter<bool>("DoMatchingToBoson");   // this is relevant for the WW, ttbar etc. samples
  int pdgIdBoson              = Options.getParameter<int>("pdgIdBoson");           // absolute value of pdgId of the boson. Can be used only if the DoMatchingToBoson is set to true.
  dRMatching                  = Options.getParameter<double>("dRMatiching");       // dR matching thresholds with the truth
  dRLeptonCleaning            = Options.getParameter<double>("dRLeptonCleaning");  // cleaning between jets and generated leptons

  //soft killer parameters
  softKillerParam = Options.getParameter<edm::ParameterSet>("softKiller");  
  //softdrop parameters
  softDropParam = Options.getParameter<std::vector<edm::ParameterSet>>("softDrop");
  //trimming
  trimmingParam = Options.getParameter<std::vector<edm::ParameterSet>>("trimming");
  //pruning
  pruningParam  = Options.getParameter<std::vector<edm::ParameterSet>>("pruning");
  //charge param
  chargeParam   = Options.getParameter<std::vector<double>>("jetcharge");
  //ECF param
  ecfParam      = Options.getParameter<std::vector<edm::ParameterSet>>("energyCorrelator");
  //Nsubjettines param
  NsubjettinessParam = Options.getParameter<std::vector<edm::ParameterSet>>("Nsubjettiness");

  //jet flavour for GenJets
  computeJetFlavour =  Options.getParameter<bool>("computeJetFlavour");
  

  // --- Read list of files to be analyzed and fill TChain 
  TChain* lTree = new TChain("Events");
  FillChain(*lTree, inputFilesList); // add all the input files in a TChain
  if (lTree->GetEntries() < maxEvents || maxEvents == -1) maxEvents = lTree->GetEntries(); 

  cout << "This analysis will run on "<< maxEvents << " events" <<endl; 

  fPFCand = new PFLoader (lTree,puppiConfig.c_str()); // Load all the PF candidates collection from baconhep structure
  if (isMC) fGen    = new GenLoader(lTree);           // Load gen particle collection 
  if (doCMSSWJets) setupCMSSWJetReadOut(lTree, jetR); // setup bacon jets if required

  TEventInfo *eventInfo = new TEventInfo(); // read the event info from baconhep input file
  lTree->SetBranchAddress("Info",&eventInfo);

  TClonesArray *PV = new TClonesArray("baconhep::TVertex"); // read vertex info from baconhep input file
  lTree->SetBranchAddress("PV",&PV);

  // --- Setup JEC on the fly  
  std::vector<JetCorrectorParameters> corrParams;
  corrParams.push_back(JetCorrectorParameters(L1FastJetJEC.c_str()));  
  corrParams.push_back(JetCorrectorParameters(L2RelativeJEC.c_str()));  
  corrParams.push_back(JetCorrectorParameters(L3AbsoluteJEC.c_str()));  
  if (L2L3ResidualJEC!="") corrParams.push_back(JetCorrectorParameters(L2L3ResidualJEC.c_str())); // 
  JetCorrectorParameters param(JECUncertainty.c_str());      
  
  FactorizedJetCorrector   *jetCorr = new FactorizedJetCorrector(corrParams); //correction
  JetCorrectionUncertainty *jetUnc  = new JetCorrectionUncertainty(param);    //uncertainty

  vtagger = new VTaggingVariables();
  
  // --- Setup JEC on the fly  for CHS
  std::vector<JetCorrectorParameters> corrParams_CHS;
  corrParams_CHS.push_back(JetCorrectorParameters(L1FastJetJEC_CHS.c_str()));  
  corrParams_CHS.push_back(JetCorrectorParameters(L2RelativeJEC_CHS.c_str()));  
  corrParams_CHS.push_back(JetCorrectorParameters(L3AbsoluteJEC_CHS.c_str()));  
  if (L2L3ResidualJEC_CHS!="") corrParams_CHS.push_back(JetCorrectorParameters(L2L3ResidualJEC_CHS.c_str())); 
  JetCorrectorParameters param_CHS(JECUncertainty_CHS.c_str());      
  
  FactorizedJetCorrector   *jetCorr_CHS = new FactorizedJetCorrector(corrParams_CHS); // correction 
  JetCorrectionUncertainty *jetUnc_CHS  = new JetCorrectionUncertainty(param_CHS); // uncertainty

  std::vector<TGraph*> lCorr;
  loadPhil(JECPhil,lCorr); // load phill corrections

  // Quark Gluon Likelihood
  qgLikelihood    = new QGLikelihoodCalculator(QGinputWeightFilePath,false);  
  qgLikelihoodCHS = new QGLikelihoodCalculator(QGinputWeightFilePath,true);  

  // --- Setup JetAlgos for basic clustering of the event
  JetDefinition jet_def(fatjet_algo,jetR);
  AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(SelectorAbsRapMax(5.0))); // real ghosts in the PseudoJet list 
  
  // --- Setup cleansing
  JetDefinition subjet_def_kt02(kt_algorithm,0.2);
  JetDefinition subjet_def_kt03(kt_algorithm,0.3);

  vector<JetCleanser> cleanser_vect;
  JetCleanser jetcleanser0 = makeJVFCleanser(subjet_def_kt03, "CMS"); cleanser_vect.push_back(jetcleanser0);
  JetCleanser jetcleanser1 = makeJVFCleanser(subjet_def_kt02, "CMS"); cleanser_vect.push_back(jetcleanser1);
  JetCleanser jetcleanser2 = makeLinearCleanser(subjet_def_kt03,0.55, "CMS"); cleanser_vect.push_back(jetcleanser2);
  JetCleanser jetcleanser3 = makeLinearCleanser(subjet_def_kt02,0.55, "CMS"); cleanser_vect.push_back(jetcleanser3);
  JetCleanser jetcleanser4 = makeLinearCleanser(subjet_def_kt03,0.60, "CMS"); cleanser_vect.push_back(jetcleanser4);
  JetCleanser jetcleanser5 = makeLinearCleanser(subjet_def_kt02,0.60, "CMS"); cleanser_vect.push_back(jetcleanser5);
  JetCleanser gsn_cleanser = makeGausCleanser(subjet_def_kt02,0.617,0.62,0.15,0.22, "CMS"); cleanser_vect.push_back(gsn_cleanser);

  // --- Setup soft-killer
  SoftKiller soft_killer (softKillerParam.getParameter<double>("ymax"),softKillerParam.getParameter<double>("cell_size"));
  
  // --- Setup output trees -> one tree for each jet collection type: GenJets, PFJets, PFCHS, Puppi, cmssw and softkiller
  TFile *fout = new TFile(fOut.c_str(),"RECREATE");
  
  TTree *genTree           = new TTree("gen"  , "gen"  );
  TTree *pfTree            = new TTree("pf"   , "pf"   );
  TTree *chsTree           = new TTree("chs"  , "chs"  );
  TTree *puppiTree         = new TTree("puppi", "puppi");
  TTree *softkillerTree    = new TTree("softkiller", "softkiller");
  TTree *cmsswTree         = new TTree("cmsswpf", "cmsswpf");
  
  GenJetInfo JGenInfo;
  JetInfo JPFInfo, JCHSInfo, JPuppiInfo, JSoftKillerInfo, JCMSSWPFInfo; // declare structures to fill the output tree information + make branches
  
  setupGenTree(genTree,   JGenInfo    , "" );
  setupTree(pfTree,    JPFInfo     , "" );
  setupTree(chsTree,   JCHSInfo    , "" );
  setupTree(puppiTree, JPuppiInfo  , "" );

  if(doSoftKillerJets) setupTree(softkillerTree, JSoftKillerInfo  , "" );
  if(doCMSSWJets) setupTree(cmsswTree, JCMSSWPFInfo, "" );
       
  // --- start loop over events
  if (minEvents < 0) minEvents = 0;
  for(int ientry = minEvents; ientry < maxEvents; ientry++) {
    // -- For each event build collections of particles (gen, puppi, etc..) to cluster as a first step
    Long64_t localEntry = lTree->LoadTree(ientry);
    fPFCand->load(localEntry); // load pF information

    // -- nPU and nPV
    lTree->GetEntry(ientry); 
    int nPU    = eventInfo->nPU;
    int nPV    = PV->GetEntries();
    double rho = eventInfo->rhoJet;

    // -- gen info (only if running on MC)
    vector<PseudoJet> genJets;
    vector<PseudoJet> gen_event;
    vector<PseudoJet> genJetsCleaned ; 
    vfloat eta_Boson, phi_Boson; // vector of eta and phi of all the vector bosons at gen level
    PseudoJet leptonVector(0.,0.,0.,0.);

    if (isMC) { 
      fGen->load(localEntry); // load gen information  
      if(fGen->leptonicBosonFilter(leptonVector) < 0) continue; // filter events With W->lnu

      gen_event       = fGen->genFetch();  //gen particles: only status 1 (ME) and user_index set 2
      fGenParticles   = fGen->GetGenParticleArray(); // take the vector of GenParticles, all the status
      ClusterSequenceArea pGen (gen_event,jet_def, area_def);
      genJets     = sorted_by_pt(pGen.inclusive_jets(10.)); // cluster only final state gen particles
     
      if (DoMatchingToBoson){
	fGen -> selectBoson(pdgIdBoson);
	eta_Boson = fGen -> eta_Boson;
	phi_Boson = fGen -> phi_Boson;
      }

      if(leptonVector.pt() > 0){  
	vector<PseudoJet>::iterator itJet = genJets.begin() ;                                                                                                                        
        for( ; itJet != genJets.end() ; ++itJet){                                                                                                                             
	  if( matchingIndex((*itJet),leptonVector,true) == false) genJetsCleaned.push_back((*itJet));                                                                   
	}  	
      }
      else genJetsCleaned = genJets;
      
      fillGenJetsInfo(genJetsCleaned, gen_event, JGenInfo, cleanser_vect, nPU, nPV, rho);          
    }
    
    vector<PseudoJet> pf_event        = fPFCand->pfFetch();       //return all the particles
    vector<PseudoJet> chs_event       = fPFCand->pfchsFetch(-1);  //only chs particles -> user_index set to 1(neutrals) or 2 (chaged from PV)
    vector<PseudoJet> puppi_event     = fPFCand->puppiFetch();    // puppi particles from all pf with puppi weights 

    // -- Cluster jets -> make the clustering
    ClusterSequenceArea pPup    (puppi_event  , jet_def, area_def);   
    ClusterSequenceArea pPF     (pf_event     , jet_def, area_def);
    ClusterSequenceArea pCHS    (chs_event    , jet_def, area_def);

    // -- Order in decreasing pt the final jet collection with an inclusive cut on jets of 25GeV
    vector<PseudoJet> puppiJets   = sorted_by_pt(pPup    .inclusive_jets(jetPtCut));
    vector<PseudoJet> pfJets      = sorted_by_pt(pPF     .inclusive_jets(jetPtCut));    
    vector<PseudoJet> chsJets     = sorted_by_pt(pCHS    .inclusive_jets(jetPtCut));
    
    vector<PseudoJet> puppiJetsCleaned ;
    vector<PseudoJet> pfJetsCleaned ;
    vector<PseudoJet> chsJetsCleaned ;

    
    // clean jets from gen lepton for semi-leptonic events    
    if(isMC && leptonVector.pt() > 0){
      vector<PseudoJet>::iterator itJet = puppiJets.begin() ;
      for( ; itJet != puppiJets.end() ; ++itJet){
        if( matchingIndex((*itJet),leptonVector,true) == false) puppiJetsCleaned.push_back((*itJet));
      }
      itJet = pfJets.begin() ;
      for( ; itJet != pfJets.end() ; ++itJet){
        if( matchingIndex((*itJet),leptonVector,true) == false) pfJetsCleaned.push_back((*itJet)); 
      }
      itJet = chsJets.begin() ;
      for( ; itJet != chsJets.end() ; ++itJet){
        if( matchingIndex((*itJet),leptonVector,true) == false) chsJetsCleaned.push_back((*itJet));
      }
    }
    else{
        puppiJetsCleaned = puppiJets ; 
        pfJetsCleaned = pfJets ; 
	chsJetsCleaned = chsJets ; 
    }   
    
    // save jet info in a tree
    fillRecoJetsInfo(puppiJetsCleaned, puppi_event, puppi_event, JPuppiInfo    , JGenInfo, false, jetCorr, jetUnc,lCorr, cleanser_vect,nPU, nPV, rho, eta_Boson, phi_Boson,true,isMC);
    fillRecoJetsInfo(pfJetsCleaned   , pf_event   , pf_event,    JPFInfo       , JGenInfo, false, jetCorr, jetUnc,lCorr, cleanser_vect,nPU, nPV, rho, eta_Boson, phi_Boson,false,isMC); 
    fillRecoJetsInfo(chsJetsCleaned  , chs_event  , pf_event,    JCHSInfo      , JGenInfo, true , jetCorr_CHS, jetUnc_CHS,lCorr, cleanser_vect, nPU, nPV, rho, eta_Boson, phi_Boson,false,isMC );      
     
    if (isMC) genTree->Fill();        
    puppiTree->Fill();
    pfTree->Fill();
    chsTree->Fill();

    if(doSoftKillerJets){ 
     vector<PseudoJet> soft_event = soft_killer(pf_event);    //retun the list from soft_killer contructor given all pf and the input parameters
     ClusterSequenceArea pSoft   (soft_event   , jet_def, area_def);
     vector<PseudoJet> softJets    = sorted_by_pt(pSoft   .inclusive_jets(jetPtCut));
     vector<PseudoJet> softJetsCleaned ;
     if(isMC && leptonVector.pt() > 0){
      vector<PseudoJet>::iterator itJet = softJets.begin() ;
      for( ; itJet != softJets.end() ; ++itJet){
        if( matchingIndex((*itJet),leptonVector,true) == false) softJetsCleaned.push_back((*itJet));
      }
     }
     else softJetsCleaned = softJets ;    
    fillRecoJetsInfo(softJetsCleaned , soft_event , soft_event,  JSoftKillerInfo , JGenInfo, true , jetCorr, jetUnc,lCorr, cleanser_vect, nPU, nPV, rho, eta_Boson, phi_Boson,false,isMC);
    softkillerTree->Fill();
    }
    
    if (doCMSSWJets) readCMSSWJet(ientry, lTree, *cmsswTree, genJets, JCMSSWPFInfo);        
    
    if (isMC) fGen->reset();         
    fPFCand->reset();
   
    cout << "===> Processed " << ientry << " - Done : " << (float(ientry-minEvents)/float(maxEvents-minEvents))*100 << "%" <<endl ;
        
   }
   
  cout<<"done event loop"<<endl;

  // --- Write trees 
  fout->cd();

  if (isMC) genTree ->Write();  
  pfTree   ->Write();
  chsTree  ->Write();
  puppiTree->Write();

  if (doSoftKillerJets) softkillerTree->Write();
  if (doCMSSWJets)  cmsswTree->Write();

  fout->Close();
  cout<<"done write trees"<<endl;
}  

 
 
