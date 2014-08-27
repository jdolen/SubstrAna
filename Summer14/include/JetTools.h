#ifndef  JETTOOLS_H
#define  JETTOOLS_H

#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"

// system include files                                                                                                                                                                   
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TGraph.h"
#include "TClonesArray.h"

// JEC class                                                                                                                                                                              
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

// Bacon objects
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TVertex.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"

#include <memory>
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <cstdlib>

using namespace std;
using namespace fastjet;


// Function for JetEnergy correction on the fly (Phill corrections)
void loadPhil(std::string iName, std::vector<TGraph*> &iCorr); 
double correctPhil(double iPt,double iEta,int iAlgo,std::vector<TGraph*> &iCorr);

// Apply standard JEC
double correction( PseudoJet &iJet,FactorizedJetCorrector *iJetCorr,double iRho);
double unc( PseudoJet &iJet,JetCorrectionUncertainty *iJetUnc);

// Jet Algorithm definition
fastjet::JetAlgorithm get_algo(string algo);

// Selector definition useful for Puppi
//// Selector for charged particles acting on pseudojets                                                                                                                                  
class SW_IsPupCharged : public SelectorWorker {
 public:
  SW_IsPupCharged(){}
  virtual bool pass(const PseudoJet & jet) const {
    return (fabs(jet.user_index()) >= 1);
  }
};

//// Selector for charged particles from pile up vertexes acting on pseudojets                                                                                                            
class SW_IsPupVertex : public SelectorWorker {
 public:
  SW_IsPupVertex(){}
  virtual bool pass(const PseudoJet & jet) const {
    return (fabs(jet.user_index()) == 2 || fabs(jet.user_index()) == 1);
  }
};


void getConstitsForCleansing(vector<PseudoJet> inputs, vector<PseudoJet> &oNeutrals, vector<PseudoJet> &oChargedLV, vector<PseudoJet> &oChargedPU);

// Matching funtions
Bool_t isMatching( fastjet::PseudoJet j1, fastjet::PseudoJet j2, Double_t deltaR=0.3);

// Jet flavour computation
int computeGenJetFlavour(const PseudoJet & iJet,TClonesArray* fGenParticles, const double & jetR);

#endif
