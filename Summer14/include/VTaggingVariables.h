#ifndef VTaggingVariables_h
#define VTaggingVariables_h

#include <iostream>
#include <math.h>

#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/GhostedAreaSpec.hh"

#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"
#include "fastjet/contrib/Nsubjettiness.hh"

#include "fastjet/contrib/EnergyCorrelator.hh"

#include "SubstrAna/Summer14/src/QjetsPlugin.h"
#include "SubstrAna/Summer14/src/Qjets.h"

#include "GenLoader.hh"
#include "PFLoader.hh"

using namespace fastjet ;

class VTaggingVariables {

 public:

  VTaggingVariables(){};
  ~VTaggingVariables(){};

  VTaggingVariables(const PseudoJet & inputJet);

  ///////////////////////

  double computeJetCharge(const double & jetChargeKappa = 1.0);
  double computeJetChargeReco(const double & jetChargeKappa = 1.0);
  double computeNSubJettines(const int & nJettines, const double & beta, const double & R0, const double & Rcut);
  double computeECF(JetAlgorithm jetAlgoforECF, const double & Rparameter, const int & nPoint, const double & beta);
  double computeQjets(const int & QJetsPreclustering = 35, const int & QJetsN = 25, const int & seed = 1);

  double getQjetVolatility(std::vector<PseudoJet> constits, const int & QJetsN, const int & seed);

  double FindMean(const std::vector< double > & qjetmasses);
  double FindRMS( const std::vector< double > & qjetmasses);
  ///////////////////////
  void   setInputJet(const PseudoJet & inputJet);

 private:

  PseudoJet inputJet_ ;
  double jetChargeKappa_ ;

  TClonesArray  *fGens_;

};


#endif
