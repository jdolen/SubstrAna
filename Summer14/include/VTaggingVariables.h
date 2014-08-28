#ifndef VTaggingVariables_h
#define VTaggingVariables_h

#include <iostream>
#include <math.h>
#include <TVector2.h>

#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

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
#include "QGLikelihoodCalculator.h"

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
  double computeECF(JetAlgorithm jetAlgoforECF, const double & Rparameter, const int & nPoint, const double & beta, const int & type = 0);
  double computeQjets(const int & QJetsPreclustering = 35, const int & QJetsN = 25, const int & seed = 1, const double & jetR = 0.8);

  double getQjetVolatility(std::vector<PseudoJet> constits, const int & QJetsN, const int & seed);

  double computeQGLikelihood(QGLikelihoodCalculator* qgLikelihood, const double & jetCorrection, const double & rho);

  double computePullAngle(std::vector<PseudoJet> subjets, const double & jetR);  
  TVector2 jetPull(const PseudoJet & jet, const int & type = 0);

  double FindMean(std::vector< double > qjetmasses);
  double FindRMS(std::vector< double > qjetmasses);
  ///////////////////////
  void   setInputJet(const PseudoJet & inputJet);

 private:

  PseudoJet inputJet_ ;
  
  std::vector<PseudoJet> particles_;
  std::vector<PseudoJet> ghosts_;

  double jetChargeKappa_ ;

  std::map<std::string,float> QCLikelihoodVariables_ ;

  TClonesArray  *fGens_;

  TVector2 jetPull_;

};


#endif
