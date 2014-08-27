#ifndef ShapeCorrectionTools_h
#define ShapeCorrectionTools_h

// ROOT lybraries
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TMath.h"
#include "TROOT.h"
#include "TLorentzVector.h"

// standard c++
#include <iostream>
#include <cstdio>   
#include <memory>
#include <string>
#include <map>
#include <fstream>

// fasjet inputs
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"

#include "fastjet/contrib/GenericSubtractor.hh"
#include "fastjet/contrib/JetCleanser.hh"

using namespace fastjet;
using namespace std;
using namespace fastjet::contrib;

//Jet Cleansing
JetCleanser makeJVFCleanser   (fastjet::JetDefinition subjet_def, std::string projectmode="CMS", double fcut=-1.0, int nsj=-1 );//projectmode: CMS or ATLAS
JetCleanser makeLinearCleanser(fastjet::JetDefinition subjet_def, double linear_para0,std::string projectmode="CMS", double fcut=-1, int nsj=-1 );//projectmode: CMS or ATLAS
JetCleanser makeGausCleanser  (fastjet::JetDefinition subjet_def, double gaus_para0, double gaus_para1, double gaus_para2, double gaus_para3, std::string projectmode="CMS", double fcut=-1, int nsj=-1 );//projectmode: CMS or ATLAS

#endif
