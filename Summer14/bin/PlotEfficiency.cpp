#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <istream>
#include <sstream>

#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "THStack.h"
#include "TString.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TMath.h"
#include "TF1.h"
#include "TH2F.h"
#include "TList.h"
#include "TRandom3.h"
#include "TTreeFormula.h"
#include "TFormula.h"
#include "TGraphAsymmErrors.h"

#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

std::string ptName  = "pt[0]";
std::string PUName  = "npu";
std::string etaName = "eta[0]";
bool saveRootOutput = false ;

/// Main programme 
int main (int argc, char **argv){
  if(argc<2){ std::cout<<" Not correct number of input parameter --> Need Just one cfg file exit  "<<std::endl; std::exit(EXIT_FAILURE); }

  // Load TTree Lybrary                                                                                                                                                                   
  gSystem->Load("libTree.so");

  // Set Root style from global enviroment path                                                                                                                                           
  std::string ROOTStyle;
  if(getenv ("ROOTStyle")!=NULL){
    ROOTStyle = getenv ("ROOTStyle");
    gROOT->ProcessLine((".x "+ROOTStyle+"/rootLogon.C").c_str());
    gROOT->ProcessLine((".x "+ROOTStyle+"/rootPalette.C").c_str());
    gROOT->ProcessLine((".x "+ROOTStyle+"/rootColors.C").c_str());
    gROOT->ProcessLine((".x "+ROOTStyle+"/setTDRStyle.C").c_str());
  }
  
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadTopMargin(0.09);
  gStyle->SetErrorX(0.5);


  // parse config file parameter --> using CMSSW parser                                                                                                                         
  std::string configFileName = argv[1];
  boost::shared_ptr<edm::ParameterSet> parameterSet = edm::readConfig(configFileName);
  edm::ParameterSet Options  = parameterSet -> getParameter<edm::ParameterSet>("Options");

  // treename -> chs, full pf and puppi
  std::string treeName ;
  if(Options.existsAs<std::string>("TreeName"))
    treeName = Options.getParameter<std::string>("TreeName");
  else{  std::cout<<" Exit from code, no input name tree provided "<<std::endl; std::exit(EXIT_FAILURE) ; }

  // preselectioncut string for Signal
  std::string PreselectionCutSignal ;
  if(Options.existsAs<std::string>("PreselectionCutSignal"))
    PreselectionCutSignal = Options.getParameter<std::string>("PreselectionCutSignal");
  else PreselectionCutSignal = "pt[0] > 200" ;
 
  // preselectioncut string for Background
  std::string PreselectionCutBackground ;
  if(Options.existsAs<std::string>("PreselectionCutBackground"))
    PreselectionCutBackground = Options.getParameter<std::string>("PreselectionCutBackground");
  else PreselectionCutBackground = "pt[0] > 200" ;

  //// BINNING //////////////////////

  // binning in pt for the signal
  std::vector<double> JetPtBinSignal;
  if(Options.existsAs<std::vector<double> >("JetPtBinSignal"))
    JetPtBinSignal = Options.getParameter<std::vector<double> >("JetPtBinSignal");
  else{ std::cout<<" Exit from code, no  jet pT bins signal"<<std::endl; std::exit(EXIT_FAILURE) ; } 

  // binning in pt for the background
  std::vector<double> JetPtBinBackground;
  if(Options.existsAs<std::vector<double> >("JetPtBinBackground"))
    JetPtBinBackground = Options.getParameter<std::vector<double> >("JetPtBinBackground");
  else{ std::cout<<" Exit from code, no  jet pT bins background"<<std::endl; std::exit(EXIT_FAILURE) ; } 

  // binning in PU
  std::vector<double> JetPUBin;
  if(Options.existsAs<std::vector<double> >("JetPUBin"))
    JetPUBin = Options.getParameter<std::vector<double> >("JetPUBin");
  else{ std::cout<<" Exit from code, no  jet PU bins "<<std::endl; std::exit(EXIT_FAILURE) ; } 

  // binning in eta
  std::vector<double> JetEtaBin;
  if(Options.existsAs<std::vector<double> >("JetEtaBin"))
    JetEtaBin = Options.getParameter<std::vector<double> >("JetEtaBin");
  else{ std::cout<<" Exit from code, no  jet Eta bins "<<std::endl; std::exit(EXIT_FAILURE) ; } 

  ////////////////////////////////////////

  // output directory for the plots
  std::string outputFileDirectory ;
  if(Options.existsAs<std::string>("outputFileDirectory"))
    outputFileDirectory = Options.getParameter<std::string>("outputFileDirectory");
  else outputFileDirectory = "cfgTraining/outputEfficiency" ;

  system(("mkdir -p "+outputFileDirectory).c_str());
  system(("rm -r "+outputFileDirectory+"/*").c_str());


  ///// FILES /////////////////////////////

  //background files
  std::vector<std::string> InputBackgroundFiles;
  if(Options.existsAs<std::vector<std::string> >("InputBackgroundFiles"))
    InputBackgroundFiles = Options.getParameter<std::vector<std::string> >("InputBackgroundFiles");
  else{ std::cout<<" Exit from code, no background files "<<std::endl; std::exit(EXIT_FAILURE) ; }

  std::vector<std::string> InputSignalFiles;
  if(Options.existsAs<std::vector<std::string> >("InputSignalFiles"))
    InputSignalFiles = Options.getParameter<std::vector<std::string> >("InputSignalFiles");
  else{ std::cout<<" Exit from code, no signal files "<<std::endl; std::exit(EXIT_FAILURE) ; }

  /////////////////////////////////////////

  /// OBSERVABLES AND CUT ///////////////////

  std::vector<edm::ParameterSet> InputMassVariables;
  if(Options.existsAs<std::vector<edm::ParameterSet> >("InputMassVariables"))
    InputMassVariables = Options.getParameter<std::vector<edm::ParameterSet> >("InputMassVariables");
  else{ std::cout<<" Exit from code, no mass observables "<<std::endl; std::exit(EXIT_FAILURE) ; }

  double MinMassCut = 0 ;
  if(Options.existsAs<double>("MinMassCut"))
    MinMassCut = Options.getParameter<double>("MinMassCut");
  else{ std::cout<<" Exit from code, no min mass cut "<<std::endl; std::exit(EXIT_FAILURE) ; }

  double MaxMassCut = 0 ;
  if(Options.existsAs<double>("MaxMassCut"))
    MaxMassCut = Options.getParameter<double>("MaxMassCut");
  else{ std::cout<<" Exit from code, no max mass cut "<<std::endl; std::exit(EXIT_FAILURE) ; }

  edm::ParameterSet WtagVariable ;
  if(Options.existsAs<edm::ParameterSet>("WtagVariable"))
    WtagVariable = Options.getParameter<edm::ParameterSet>("WtagVariable");
  else{ std::cout<<" Exit from code, no Wtag info "<<std::endl; std::exit(EXIT_FAILURE) ; }

  double WtagCut = 0 ;
  if(Options.existsAs<double>("WtagCut"))
    WtagCut = Options.getParameter<double>("WtagCut");
  else{ std::cout<<" Exit from code, no W-tag cut "<<std::endl; std::exit(EXIT_FAILURE) ; }

  /////////////////////////////////////////

  std::vector<TTree*> inputBackgroundTrees ;  // TTree for input variables inside the TMVA output file
  std::vector<TTree*> inputSignalTrees ;  // TTree for input variables inside the TMVA output file
  TFile* inputFile = NULL ;

  std::vector<std::string>::const_iterator itBack = InputBackgroundFiles.begin();
  for( ; itBack != InputBackgroundFiles.end() ; ++itBack){
    inputFile = TFile::Open((*itBack).c_str());
    if(inputFile == 0 or inputFile == NULL) continue;
    inputBackgroundTrees.push_back((TTree*) inputFile->Get(treeName.c_str())); // make the tree list
  }
   
  std::vector<std::string>::const_iterator itSig = InputSignalFiles.begin();
  for( ; itSig != InputSignalFiles.end() ; ++itSig){
    inputFile = TFile::Open((*itSig).c_str());
    if(inputFile == 0 or inputFile == NULL) continue;
    inputSignalTrees.push_back((TTree*) inputFile->Get(treeName.c_str())); // make the tree list
  }

   
  // histogram definition -> vs PT, vs PU and vs eta
  std::vector<std::vector<TH1F*> > denominatorBackgroundMass_vsPT ; denominatorBackgroundMass_vsPT.resize(inputBackgroundTrees.size());
  std::vector<std::vector<TH1F*> > numeratorBackgroundMass_vsPT ; numeratorBackgroundMass_vsPT.resize(inputBackgroundTrees.size());
  std::vector<std::vector<TH1F*> > numeratorBackgroundWtag_vsPT ; numeratorBackgroundWtag_vsPT.resize(inputBackgroundTrees.size());

  std::vector<std::vector<TH1F*> > denominatorSignalMass_vsPT ; denominatorSignalMass_vsPT.resize(inputSignalTrees.size());
  std::vector<std::vector<TH1F*> > numeratorSignalMass_vsPT ; numeratorSignalMass_vsPT.resize(inputSignalTrees.size());
  std::vector<std::vector<TH1F*> > numeratorSignalWtag_vsPT ; numeratorSignalWtag_vsPT.resize(inputSignalTrees.size());

  std::vector<std::vector<TH1F*> > denominatorBackgroundMass_vsEta ; denominatorBackgroundMass_vsEta.resize(inputBackgroundTrees.size());
  std::vector<std::vector<TH1F*> > numeratorBackgroundMass_vsEta ; numeratorBackgroundMass_vsEta.resize(inputBackgroundTrees.size());
  std::vector<std::vector<TH1F*> > numeratorBackgroundWtag_vsEta ; numeratorBackgroundWtag_vsEta.resize(inputBackgroundTrees.size());

  std::vector<std::vector<TH1F*> > denominatorSignalMass_vsEta ; denominatorSignalMass_vsEta.resize(inputSignalTrees.size());
  std::vector<std::vector<TH1F*> > numeratorSignalMass_vsEta ; numeratorSignalMass_vsEta.resize(inputSignalTrees.size());
  std::vector<std::vector<TH1F*> > numeratorSignalWtag_vsEta ; numeratorSignalWtag_vsEta.resize(inputSignalTrees.size());

  std::vector<std::vector<TH1F*> > denominatorBackgroundMass_vsPU ; denominatorBackgroundMass_vsPU.resize(inputBackgroundTrees.size());
  std::vector<std::vector<TH1F*> > numeratorBackgroundMass_vsPU ; numeratorBackgroundMass_vsPU.resize(inputBackgroundTrees.size());
  std::vector<std::vector<TH1F*> > numeratorBackgroundWtag_vsPU ; numeratorBackgroundWtag_vsPU.resize(inputBackgroundTrees.size());

  std::vector<std::vector<TH1F*> > denominatorSignalMass_vsPU ; denominatorSignalMass_vsPU.resize(inputSignalTrees.size());
  std::vector<std::vector<TH1F*> > numeratorSignalMass_vsPU ; numeratorSignalMass_vsPU.resize(inputSignalTrees.size());
  std::vector<std::vector<TH1F*> > numeratorSignalWtag_vsPU ; numeratorSignalWtag_vsPU.resize(inputSignalTrees.size());

  double* ptBinSignal     = &JetPtBinSignal[0] ;
  double* ptBinBackground = &JetPtBinBackground[0] ;
  double* pUBin  = &JetPUBin[0] ;
  double* etaBin = &JetEtaBin[0] ;

   for(size_t iTree = 0; iTree < inputBackgroundTrees.size() ; iTree++){
    for( size_t iMassVar = 0; iMassVar < InputMassVariables.size(); iMassVar++){

      denominatorBackgroundMass_vsPT.at(iTree).push_back(new TH1F(Form("back_deno_%s_%d_vsPT",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),int(iTree)),"",JetPtBinBackground.size()-1,ptBinBackground));
      numeratorBackgroundMass_vsPT.at(iTree).push_back(new TH1F(Form("back_num_%s_%d_vsPT",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),int(iTree)),"",JetPtBinBackground.size()-1,ptBinBackground));
      numeratorBackgroundWtag_vsPT.at(iTree).push_back(new TH1F(Form("back_num_wtag_%s_%d_vsPT",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),int(iTree)),"",JetPtBinBackground.size()-1,ptBinBackground));
      denominatorBackgroundMass_vsPT.at(iTree).back()->Sumw2();
      numeratorBackgroundMass_vsPT.at(iTree).back()->Sumw2();
      numeratorBackgroundWtag_vsPT.at(iTree).back()->Sumw2();

      denominatorBackgroundMass_vsEta.at(iTree).push_back(new TH1F(Form("back_deno_%s_%d_vsEta",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),int(iTree)),"",JetEtaBin.size()-1,etaBin));
      numeratorBackgroundMass_vsEta.at(iTree).push_back(new TH1F(Form("back_num_%s_%d_vsEta",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),int(iTree)),"",JetEtaBin.size()-1,etaBin));
      numeratorBackgroundWtag_vsEta.at(iTree).push_back(new TH1F(Form("back_num_wtag_%s_%d_vsEta",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),int(iTree)),"",JetEtaBin.size()-1,etaBin));
      denominatorBackgroundMass_vsEta.at(iTree).back()->Sumw2();
      numeratorBackgroundMass_vsEta.at(iTree).back()->Sumw2();
      numeratorBackgroundWtag_vsEta.at(iTree).back()->Sumw2();

      denominatorBackgroundMass_vsPU.at(iTree).push_back(new TH1F(Form("back_deno_%s_%d_vsPU",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),int(iTree)),"",JetPUBin.size()-1,pUBin));
      numeratorBackgroundMass_vsPU.at(iTree).push_back(new TH1F(Form("back_num_%s_%d_vsPU",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),int(iTree)),"",JetPUBin.size()-1,pUBin));
      numeratorBackgroundWtag_vsPU.at(iTree).push_back(new TH1F(Form("back_num_wtag_%s_%d_vsPU",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),int(iTree)),"",JetPUBin.size()-1,pUBin));
      denominatorBackgroundMass_vsPU.at(iTree).back()->Sumw2();
      numeratorBackgroundMass_vsPU.at(iTree).back()->Sumw2();
      numeratorBackgroundWtag_vsPU.at(iTree).back()->Sumw2();
    }
   }   

   for(size_t iTree = 0; iTree < inputSignalTrees.size() ; iTree++){
    for( size_t iMassVar = 0; iMassVar < InputMassVariables.size(); iMassVar++){

      denominatorSignalMass_vsPT.at(iTree).push_back(new TH1F(Form("sig_deno_%s_%d_vsPT",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),int(iTree)),"",JetPtBinSignal.size()-1,ptBinSignal));
      numeratorSignalMass_vsPT.at(iTree).push_back(new TH1F(Form("num_%s_%d_vsPT",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),int(iTree)),"",JetPtBinSignal.size()-1,ptBinSignal));
      numeratorSignalWtag_vsPT.at(iTree).push_back(new TH1F(Form("num_wtag_%s_%d_vsPT",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),int(iTree)),"",JetPtBinSignal.size()-1,ptBinSignal));
      denominatorSignalMass_vsPT.at(iTree).back()->Sumw2();
      numeratorSignalMass_vsPT.at(iTree).back()->Sumw2();
      numeratorSignalWtag_vsPT.at(iTree).back()->Sumw2();

      denominatorSignalMass_vsEta.at(iTree).push_back(new TH1F(Form("sig_deno_%s_%d_vsEta",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),int(iTree)),"",JetEtaBin.size()-1,etaBin));
      numeratorSignalMass_vsEta.at(iTree).push_back(new TH1F(Form("num_%s_%d_vsEta",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),int(iTree)),"",JetEtaBin.size()-1,etaBin));
      numeratorSignalWtag_vsEta.at(iTree).push_back(new TH1F(Form("num_wtag_%s_%d_vsEta",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),int(iTree)),"",JetEtaBin.size()-1,etaBin));
      denominatorSignalMass_vsEta.at(iTree).back()->Sumw2();
      numeratorSignalMass_vsEta.at(iTree).back()->Sumw2();
      numeratorSignalWtag_vsEta.at(iTree).back()->Sumw2();

      denominatorSignalMass_vsPU.at(iTree).push_back(new TH1F(Form("sig_deno_%s_%d_vsPU",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),int(iTree)),"",JetPUBin.size()-1,pUBin));
      numeratorSignalMass_vsPU.at(iTree).push_back(new TH1F(Form("num_%s_%d_vsPU",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),int(iTree)),"",JetPUBin.size()-1,pUBin));
      numeratorSignalWtag_vsPU.at(iTree).push_back(new TH1F(Form("num_wtag_%s_%d_vsPU",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),int(iTree)),"",JetPUBin.size()-1,pUBin));
      denominatorSignalMass_vsPU.at(iTree).back()->Sumw2();
      numeratorSignalMass_vsPU.at(iTree).back()->Sumw2();
      numeratorSignalWtag_vsPU.at(iTree).back()->Sumw2();
    }
   }  
  
   for(size_t iTree = 0; iTree < inputBackgroundTrees.size() ; iTree++){
    for( size_t iMassVar = 0; iMassVar < InputMassVariables.size(); iMassVar++){
      inputBackgroundTrees.at(iTree)->Draw(Form("%s >> %s ",ptName.c_str(),denominatorBackgroundMass_vsPT.at(iTree).at(iMassVar)->GetName()),
                                           Form("%s",PreselectionCutBackground.c_str()),"goff");

      inputBackgroundTrees.at(iTree)->Draw(Form("%s >> %s ",ptName.c_str(),numeratorBackgroundMass_vsPT.at(iTree).at(iMassVar)->GetName()),
                                           Form("%s && %s > %f && %s < %f",PreselectionCutBackground.c_str(),InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),MinMassCut,InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),MaxMassCut),"goff");

      inputBackgroundTrees.at(iTree)->Draw(Form("%s >> %s ",ptName.c_str(),numeratorBackgroundWtag_vsPT.at(iTree).at(iMassVar)->GetName()),
                                           Form("%s && %s > %f && %s < %f && %s < %f ",PreselectionCutBackground.c_str(),InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),MinMassCut,InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),MaxMassCut,WtagVariable.getParameter<std::string>("VariableName").c_str(),WtagCut),"goff");



      inputBackgroundTrees.at(iTree)->Draw(Form("%s >> %s ",etaName.c_str(),denominatorBackgroundMass_vsEta.at(iTree).at(iMassVar)->GetName()),
                                           Form("%s",PreselectionCutBackground.c_str()),"goff");

      inputBackgroundTrees.at(iTree)->Draw(Form("%s >> %s ",etaName.c_str(),numeratorBackgroundMass_vsEta.at(iTree).at(iMassVar)->GetName()),
                                           Form("%s && %s > %f && %s < %f",PreselectionCutBackground.c_str(),InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),MinMassCut,InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),MaxMassCut),"goff");

      inputBackgroundTrees.at(iTree)->Draw(Form("%s >> %s ",etaName.c_str(),numeratorBackgroundWtag_vsEta.at(iTree).at(iMassVar)->GetName()),
                                           Form("%s && %s > %f && %s < %f && %s < %f ",PreselectionCutBackground.c_str(),InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),MinMassCut,InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),MaxMassCut,WtagVariable.getParameter<std::string>("VariableName").c_str(),WtagCut),"goff");


      inputBackgroundTrees.at(iTree)->Draw(Form("%s >> %s ",PUName.c_str(),denominatorBackgroundMass_vsPU.at(iTree).at(iMassVar)->GetName()),
                                           Form("%s",PreselectionCutBackground.c_str()),"goff");

      inputBackgroundTrees.at(iTree)->Draw(Form("%s >> %s ",PUName.c_str(),numeratorBackgroundMass_vsPU.at(iTree).at(iMassVar)->GetName()),
                                           Form("%s && %s > %f && %s < %f",PreselectionCutBackground.c_str(),InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),MinMassCut,InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),MaxMassCut),"goff");

      inputBackgroundTrees.at(iTree)->Draw(Form("%s >> %s ",PUName.c_str(),numeratorBackgroundWtag_vsPU.at(iTree).at(iMassVar)->GetName()),
                                           Form("%s && %s > %f && %s < %f && %s < %f ",PreselectionCutBackground.c_str(),InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),MinMassCut,InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),MaxMassCut,WtagVariable.getParameter<std::string>("VariableName").c_str(),WtagCut),"goff");

    }
   }

   for(size_t iTree = 0; iTree < inputSignalTrees.size() ; iTree++){
    for( size_t iMassVar = 0; iMassVar < InputMassVariables.size(); iMassVar++){

      inputSignalTrees.at(iTree)->Draw(Form("%s >> %s ",ptName.c_str(),denominatorSignalMass_vsPT.at(iTree).at(iMassVar)->GetName()),
                                       Form("%s",PreselectionCutSignal.c_str()),"goff");
  
      inputSignalTrees.at(iTree)->Draw(Form("%s >> %s ",ptName.c_str(),numeratorSignalMass_vsPT.at(iTree).at(iMassVar)->GetName()),
                                       Form("%s && %s > %f && %s < %f",PreselectionCutSignal.c_str(),InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),
                                             MinMassCut,InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),MaxMassCut),"goff");

      inputSignalTrees.at(iTree)->Draw(Form("%s >> %s ",ptName.c_str(),numeratorSignalWtag_vsPT.at(iTree).at(iMassVar)->GetName()),
                                            Form("%s && %s > %f && %s < %f && %s < %f ",PreselectionCutSignal.c_str(),InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),MinMassCut,InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),MaxMassCut,WtagVariable.getParameter<std::string>("VariableName").c_str(),WtagCut),"goff");


      inputSignalTrees.at(iTree)->Draw(Form("%s >> %s ",etaName.c_str(),denominatorSignalMass_vsEta.at(iTree).at(iMassVar)->GetName()),
                                       Form("%s",PreselectionCutSignal.c_str()),"goff");
  
      inputSignalTrees.at(iTree)->Draw(Form("%s >> %s ",etaName.c_str(),numeratorSignalMass_vsEta.at(iTree).at(iMassVar)->GetName()),
                                       Form("%s && %s > %f && %s < %f",PreselectionCutSignal.c_str(),InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),
                                             MinMassCut,InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),MaxMassCut),"goff");

      inputSignalTrees.at(iTree)->Draw(Form("%s >> %s ",etaName.c_str(),numeratorSignalWtag_vsEta.at(iTree).at(iMassVar)->GetName()),
                                            Form("%s && %s > %f && %s < %f && %s < %f ",PreselectionCutSignal.c_str(),InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),MinMassCut,InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),MaxMassCut,WtagVariable.getParameter<std::string>("VariableName").c_str(),WtagCut),"goff");



      inputSignalTrees.at(iTree)->Draw(Form("%s >> %s ",PUName.c_str(),denominatorSignalMass_vsPU.at(iTree).at(iMassVar)->GetName()),
                                       Form("%s",PreselectionCutSignal.c_str()),"goff");
  
      inputSignalTrees.at(iTree)->Draw(Form("%s >> %s ",PUName.c_str(),numeratorSignalMass_vsPU.at(iTree).at(iMassVar)->GetName()),
                                       Form("%s && %s > %f && %s < %f",PreselectionCutSignal.c_str(),InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),
                                             MinMassCut,InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),MaxMassCut),"goff");

      inputSignalTrees.at(iTree)->Draw(Form("%s >> %s ",PUName.c_str(),numeratorSignalWtag_vsPU.at(iTree).at(iMassVar)->GetName()),
                                            Form("%s && %s > %f && %s < %f && %s < %f ",PreselectionCutSignal.c_str(),InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),MinMassCut,InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),MaxMassCut,WtagVariable.getParameter<std::string>("VariableName").c_str(),WtagCut),"goff");
    }
   }

   std::vector<TGraphAsymmErrors*> signalEfficiencyMass_vsPT ; 
   std::vector<TGraphAsymmErrors*> backgroundEfficiencyMass_vsPT ;
   std::vector<TGraphAsymmErrors*> signalEfficiencyWtag_vsPT ; 
   std::vector<TGraphAsymmErrors*> backgroundEfficiencyWtag_vsPT ;

   TH1F* denTempSig_vsPT = new TH1F(Form("numTempMassSig_vsPT"),"",JetPtBinSignal.size()-1,ptBinSignal);   
   TH1F* numTempMassSig_vsPT = new TH1F(Form("numTempWtagSig_vsPT"),"",JetPtBinSignal.size()-1,ptBinSignal);
   TH1F* numTempWtagSig_vsPT = new TH1F(Form("denTempSig_vsPT"),"",JetPtBinSignal.size()-1,ptBinSignal);

   numTempMassSig_vsPT->Sumw2();
   numTempWtagSig_vsPT->Sumw2();
   denTempSig_vsPT->Sumw2();

   TH1F* denTempBack_vsPT = new TH1F(Form("numTempMassBack_vsPT"),"",JetPtBinBackground.size()-1,ptBinBackground);   
   TH1F* numTempMassBack_vsPT = new TH1F(Form("numTempWtagBack_vsPT"),"",JetPtBinBackground.size()-1,ptBinBackground);
   TH1F* numTempWtagBack_vsPT = new TH1F(Form("denTempback_vsPT"),"",JetPtBinBackground.size()-1,ptBinBackground);

   numTempMassBack_vsPT->Sumw2();
   numTempWtagBack_vsPT->Sumw2();
   denTempBack_vsPT->Sumw2();


   std::vector<TGraphAsymmErrors*> signalEfficiencyMass_vsEta ; 
   std::vector<TGraphAsymmErrors*> backgroundEfficiencyMass_vsEta ;
   std::vector<TGraphAsymmErrors*> signalEfficiencyWtag_vsEta ; 
   std::vector<TGraphAsymmErrors*> backgroundEfficiencyWtag_vsEta ;

   TH1F* denTemp_vsEta = new TH1F(Form("numTempMass_vsEta"),"",JetEtaBin.size()-1,etaBin);   
   TH1F* numTempMass_vsEta = new TH1F(Form("numTempWtag_vsEta"),"",JetEtaBin.size()-1,etaBin);
   TH1F* numTempWtag_vsEta = new TH1F(Form("denTemp_vsEta"),"",JetEtaBin.size()-1,etaBin);

   numTempMass_vsEta->Sumw2();
   numTempWtag_vsEta->Sumw2();
   denTemp_vsEta->Sumw2();
    

   std::vector<TGraphAsymmErrors*> signalEfficiencyMass_vsPU ; 
   std::vector<TGraphAsymmErrors*> backgroundEfficiencyMass_vsPU ;
   std::vector<TGraphAsymmErrors*> signalEfficiencyWtag_vsPU ; 
   std::vector<TGraphAsymmErrors*> backgroundEfficiencyWtag_vsPU ;

   TH1F* denTemp_vsPU = new TH1F(Form("numTempMass_vsPU"),"",JetPUBin.size()-1,pUBin);   
   TH1F* numTempMass_vsPU = new TH1F(Form("numTempWtag_vsPU"),"",JetPUBin.size()-1,pUBin);
   TH1F* numTempWtag_vsPU = new TH1F(Form("denTemp_vsPU"),"",JetPUBin.size()-1,pUBin);

   numTempMass_vsPU->Sumw2();
   numTempWtag_vsPU->Sumw2();
   denTemp_vsPU->Sumw2();
    
   for(size_t iMassVar = 0; iMassVar < InputMassVariables.size(); iMassVar++){
     
     numTempMassSig_vsPT->Reset();
     numTempWtagSig_vsPT->Reset();
     denTempSig_vsPT->Reset();

     signalEfficiencyMass_vsPT.push_back(new TGraphAsymmErrors());
     signalEfficiencyMass_vsPT.back()->SetName(Form("signalEfficiencyMass_vsPT_%s",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str()));
     backgroundEfficiencyMass_vsPT.push_back(new TGraphAsymmErrors());
     backgroundEfficiencyMass_vsPT.back()->SetName(Form("backgroundEfficiencyMass_vsPT_%s",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str()));

     signalEfficiencyWtag_vsPT.push_back(new TGraphAsymmErrors());
     signalEfficiencyWtag_vsPT.back()->SetName(Form("signalEfficiencyWtag_vsPT_%s",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str()));
     backgroundEfficiencyWtag_vsPT.push_back(new TGraphAsymmErrors());
     backgroundEfficiencyWtag_vsPT.back()->SetName(Form("backgroundEfficiencyWtag_vsPT_%s",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str()));

     for(size_t iTree = 0; iTree < inputSignalTrees.size() ; iTree++){
       numTempMassSig_vsPT->Add(numeratorSignalMass_vsPT.at(iTree).at(iMassVar));
       numTempWtagSig_vsPT->Add(numeratorSignalWtag_vsPT.at(iTree).at(iMassVar));
       denTempSig_vsPT->Add(denominatorSignalMass_vsPT.at(iTree).at(iMassVar));
     }

     signalEfficiencyMass_vsPT.back()->BayesDivide(numTempMassSig_vsPT,denTempSig_vsPT);
     signalEfficiencyWtag_vsPT.back()->BayesDivide(numTempWtagSig_vsPT,denTempSig_vsPT);
     
     numTempMassBack_vsPT->Reset();
     numTempWtagBack_vsPT->Reset();
     denTempBack_vsPT->Reset();

     for(size_t iTree = 0; iTree < inputBackgroundTrees.size() ; iTree++){
       numTempMassBack_vsPT->Add(numeratorBackgroundMass_vsPT.at(iTree).at(iMassVar));
       numTempWtagBack_vsPT->Add(numeratorBackgroundWtag_vsPT.at(iTree).at(iMassVar));
       denTempBack_vsPT->Add(denominatorBackgroundMass_vsPT.at(iTree).at(iMassVar));
     }

     backgroundEfficiencyMass_vsPT.back()->BayesDivide(numTempMassBack_vsPT,denTempBack_vsPT);
     backgroundEfficiencyWtag_vsPT.back()->BayesDivide(numTempWtagBack_vsPT,denTempBack_vsPT);

     //////////////////////////////////////

     numTempMass_vsEta->Reset();
     numTempWtag_vsEta->Reset();
     denTemp_vsEta->Reset();

     signalEfficiencyMass_vsEta.push_back(new TGraphAsymmErrors());
     signalEfficiencyMass_vsEta.back()->SetName(Form("signalEfficiencyMass_vsEta_%s",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str()));
     backgroundEfficiencyMass_vsEta.push_back(new TGraphAsymmErrors());
     backgroundEfficiencyMass_vsEta.back()->SetName(Form("backgroundEfficiencyMass_vsEta_%s",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str()));

     signalEfficiencyWtag_vsEta.push_back(new TGraphAsymmErrors());
     signalEfficiencyWtag_vsEta.back()->SetName(Form("signalEfficiencyWtag_vsEta_%s",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str()));
     backgroundEfficiencyWtag_vsEta.push_back(new TGraphAsymmErrors());
     backgroundEfficiencyWtag_vsEta.back()->SetName(Form("backgroundEfficiencyWtag_vsEta_%s",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str()));

     for(size_t iTree = 0; iTree < inputSignalTrees.size() ; iTree++){
       numTempMass_vsEta->Add(numeratorSignalMass_vsEta.at(iTree).at(iMassVar));
       numTempWtag_vsEta->Add(numeratorSignalWtag_vsEta.at(iTree).at(iMassVar));
       denTemp_vsEta->Add(denominatorSignalMass_vsEta.at(iTree).at(iMassVar));
     }

     signalEfficiencyMass_vsEta.back()->BayesDivide(numTempMass_vsEta,denTemp_vsEta);
     signalEfficiencyWtag_vsEta.back()->BayesDivide(numTempWtag_vsEta,denTemp_vsEta);
     
     numTempMass_vsEta->Reset();
     numTempWtag_vsEta->Reset();
     denTemp_vsEta->Reset();

     for(size_t iTree = 0; iTree < inputBackgroundTrees.size() ; iTree++){
       numTempMass_vsEta->Add(numeratorBackgroundMass_vsEta.at(iTree).at(iMassVar));
       numTempWtag_vsEta->Add(numeratorBackgroundWtag_vsEta.at(iTree).at(iMassVar));
       denTemp_vsEta->Add(denominatorBackgroundMass_vsEta.at(iTree).at(iMassVar));
     }

     backgroundEfficiencyMass_vsEta.back()->BayesDivide(numTempMass_vsEta,denTemp_vsEta);
     backgroundEfficiencyWtag_vsEta.back()->BayesDivide(numTempWtag_vsEta,denTemp_vsEta);

     //////////////////////////////////////

     numTempMass_vsPU->Reset();
     numTempWtag_vsPU->Reset();
     denTemp_vsPU->Reset();

     signalEfficiencyMass_vsPU.push_back(new TGraphAsymmErrors());
     signalEfficiencyMass_vsPU.back()->SetName(Form("signalEfficiencyMass_vsPU_%s",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str()));
     backgroundEfficiencyMass_vsPU.push_back(new TGraphAsymmErrors());
     backgroundEfficiencyMass_vsPU.back()->SetName(Form("backgroundEfficiencyMass_vsPU_%s",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str()));

     signalEfficiencyWtag_vsPU.push_back(new TGraphAsymmErrors());
     signalEfficiencyWtag_vsPU.back()->SetName(Form("signalEfficiencyWtag_vsPU_%s",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str()));
     backgroundEfficiencyWtag_vsPU.push_back(new TGraphAsymmErrors());
     backgroundEfficiencyWtag_vsPU.back()->SetName(Form("backgroundEfficiencyWtag_vsPU_%s",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str()));

     for(size_t iTree = 0; iTree < inputSignalTrees.size() ; iTree++){
       numTempMass_vsPU->Add(numeratorSignalMass_vsPU.at(iTree).at(iMassVar));
       numTempWtag_vsPU->Add(numeratorSignalWtag_vsPU.at(iTree).at(iMassVar));
       denTemp_vsPU->Add(denominatorSignalMass_vsPU.at(iTree).at(iMassVar));
     }

     signalEfficiencyMass_vsPU.back()->BayesDivide(numTempMass_vsPU,denTemp_vsPU);
     signalEfficiencyWtag_vsPU.back()->BayesDivide(numTempWtag_vsPU,denTemp_vsPU);
     
     numTempMass_vsPU->Reset();
     numTempWtag_vsPU->Reset();
     denTemp_vsPU->Reset();

     for(size_t iTree = 0; iTree < inputBackgroundTrees.size() ; iTree++){
       numTempMass_vsPU->Add(numeratorBackgroundMass_vsPU.at(iTree).at(iMassVar));
       numTempWtag_vsPU->Add(numeratorBackgroundWtag_vsPU.at(iTree).at(iMassVar));
       denTemp_vsPU->Add(denominatorBackgroundMass_vsPU.at(iTree).at(iMassVar));
     }

     backgroundEfficiencyMass_vsPU.back()->BayesDivide(numTempMass_vsPU,denTemp_vsPU);
     backgroundEfficiencyWtag_vsPU.back()->BayesDivide(numTempWtag_vsPU,denTemp_vsPU);


   }

   // DUMP single efficiency plot in a root file
   if(saveRootOutput){
    TFile* output  = new TFile("output.root","RECREATE");

    for(size_t iSig = 0 ; iSig < signalEfficiencyMass_vsPT.size(); iSig++)
     signalEfficiencyMass_vsPT.at(iSig)->Write();

    for(size_t iSig = 0 ; iSig < signalEfficiencyWtag_vsPT.size(); iSig++)
     signalEfficiencyWtag_vsPT.at(iSig)->Write();
   
    for(size_t iSig = 0 ; iSig < backgroundEfficiencyMass_vsPT.size(); iSig++)
     backgroundEfficiencyMass_vsPT.at(iSig)->Write();

    for(size_t iSig = 0 ; iSig < backgroundEfficiencyWtag_vsPT.size(); iSig++)
     backgroundEfficiencyWtag_vsPT.at(iSig)->Write();

    for(size_t iSig = 0 ; iSig < signalEfficiencyMass_vsEta.size(); iSig++)
     signalEfficiencyMass_vsEta.at(iSig)->Write();

    for(size_t iSig = 0 ; iSig < signalEfficiencyWtag_vsEta.size(); iSig++)
     signalEfficiencyWtag_vsEta.at(iSig)->Write();
   
    for(size_t iSig = 0 ; iSig < backgroundEfficiencyMass_vsEta.size(); iSig++)
     backgroundEfficiencyMass_vsEta.at(iSig)->Write();

    for(size_t iSig = 0 ; iSig < backgroundEfficiencyWtag_vsEta.size(); iSig++)
     backgroundEfficiencyWtag_vsEta.at(iSig)->Write();

    for(size_t iSig = 0 ; iSig < signalEfficiencyMass_vsPU.size(); iSig++)
     signalEfficiencyMass_vsPU.at(iSig)->Write();

    for(size_t iSig = 0 ; iSig < signalEfficiencyWtag_vsPU.size(); iSig++)
     signalEfficiencyWtag_vsPU.at(iSig)->Write();
   
    for(size_t iSig = 0 ; iSig < backgroundEfficiencyMass_vsPU.size(); iSig++)
     backgroundEfficiencyMass_vsPU.at(iSig)->Write();

    for(size_t iSig = 0 ; iSig < backgroundEfficiencyWtag_vsPU.size(); iSig++)
     backgroundEfficiencyWtag_vsPU.at(iSig)->Write();
   
    output->Close();
   }

   // make the plot vs PT
   TCanvas *cSignal_vsPT = new TCanvas("cSignal_vsPT","",180,52,550,550);
   
   cSignal_vsPT->SetTicks();
   cSignal_vsPT->SetFillColor(0);
   cSignal_vsPT->SetBorderMode(0);
   cSignal_vsPT->SetBorderSize(2);
   cSignal_vsPT->SetTickx(1);
   cSignal_vsPT->SetTicky(1);
   cSignal_vsPT->SetRightMargin(0.05);
   cSignal_vsPT->SetBottomMargin(0.12);
   cSignal_vsPT->SetFrameBorderMode(0);

   TH2F* frameSignal_vsPT = new TH2F("frameSignal_vsPT","",500,JetPtBinSignal.at(0),JetPtBinSignal.back(),500,0,1);
   frameSignal_vsPT->SetLineWidth(2);
   frameSignal_vsPT->SetMarkerStyle(21);
   frameSignal_vsPT->SetMarkerSize(0.3);
   frameSignal_vsPT->GetXaxis()->SetNdivisions(505);
   frameSignal_vsPT->GetYaxis()->SetNdivisions(505);
   frameSignal_vsPT->GetXaxis()->SetTitle("AK8 jet p_{T} (GeV)");
   frameSignal_vsPT->GetXaxis()->SetLabelOffset(0.012);
   frameSignal_vsPT->GetXaxis()->SetLabelSize(0.038);
   frameSignal_vsPT->GetXaxis()->SetTitleSize(0.05);
   frameSignal_vsPT->GetXaxis()->SetTitleOffset(1.10);
   frameSignal_vsPT->GetYaxis()->SetTitle("efficiency");
   frameSignal_vsPT->GetYaxis()->SetLabelOffset(0.012);
   frameSignal_vsPT->GetYaxis()->SetLabelSize(0.038);
   frameSignal_vsPT->GetYaxis()->SetTitleSize(0.05);
   frameSignal_vsPT->GetYaxis()->SetTitleOffset(1.18);

   TLatex *   tex = new TLatex(0.94,0.92," 13 TeV");
   tex->SetNDC();
   tex->SetTextAlign(31);
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   TLatex * tex2 = new TLatex(0.14,0.92,"CMS");
   tex2->SetNDC();
   tex2->SetTextFont(61);
   tex2->SetTextSize(0.04);
   tex2->SetLineWidth(2);
   TLatex * tex3 = new TLatex(0.236,0.92,"Simulation Preliminary");
   tex3->SetNDC();
   tex3->SetTextFont(52);
   tex3->SetTextSize(0.035);
   tex3->SetLineWidth(2);

   TLegend* SignalLegend = new TLegend(0.175,0.2,0.53,0.4);
   SignalLegend->SetBorderSize(0);
   SignalLegend->SetFillColor(0);
   SignalLegend->SetFillStyle(0);
   SignalLegend->SetTextSize(0.031);
   SignalLegend->SetTextFont(42);

   TLatex *   signalBanner = new TLatex(0.17,0.82,"X #rightarrow W_{T} W_{T} Pythia8");
   signalBanner->SetNDC();
   signalBanner->SetTextAlign(11);
   signalBanner->SetTextFont(42);
   signalBanner->SetTextSize(0.036);
   signalBanner->SetLineWidth(2);

   TLatex *   jetBanner = new TLatex(0.17,0.77,Form("%s jets",treeName.c_str()));
   jetBanner->SetNDC();
   jetBanner->SetTextAlign(11);
   jetBanner->SetTextFont(42);
   jetBanner->SetTextSize(0.03);
   jetBanner->SetLineWidth(2);


   for( size_t iSignal = 0 ; iSignal<signalEfficiencyMass_vsPT.size(); iSignal++){
     cSignal_vsPT->cd();
     frameSignal_vsPT->Draw();
     tex->Draw("same");
     tex2->Draw("same");
     tex3->Draw("same");
     signalBanner->Draw("same");
     jetBanner->Draw("same");

     signalEfficiencyMass_vsPT.at(iSignal)->SetLineColor(kBlack); 
     signalEfficiencyMass_vsPT.at(iSignal)->SetMarkerColor(kBlack); 
     signalEfficiencyMass_vsPT.at(iSignal)->SetLineWidth(2); 
     signalEfficiencyMass_vsPT.at(iSignal)->SetMarkerStyle(20); 
     signalEfficiencyMass_vsPT.at(iSignal)->Draw("psame");

     signalEfficiencyWtag_vsPT.at(iSignal)->SetLineColor(kRed); 
     signalEfficiencyWtag_vsPT.at(iSignal)->SetMarkerColor(kRed); 
     signalEfficiencyWtag_vsPT.at(iSignal)->SetLineWidth(2); 
     signalEfficiencyWtag_vsPT.at(iSignal)->SetMarkerStyle(22); 
     signalEfficiencyWtag_vsPT.at(iSignal)->Draw("psame");

     SignalLegend->Clear();
     SignalLegend->AddEntry(signalEfficiencyMass_vsPT.at(iSignal),Form("%s #in [%d,%d] GeV",InputMassVariables.at(iSignal).getParameter<std::string>("ReducedName").c_str(),int(MinMassCut),int(MaxMassCut)),"pl");
     SignalLegend->AddEntry(signalEfficiencyWtag_vsPT.at(iSignal),Form("%s #in [%d,%d] GeV + %s < %0.1f",InputMassVariables.at(iSignal).getParameter<std::string>("ReducedName").c_str(),int(MinMassCut),int(MaxMassCut),WtagVariable.getParameter<std::string>("ReducedName").c_str(),WtagCut),"pl");
     SignalLegend->Draw("same");

     TString VariableName = Form("%s",InputMassVariables.at(iSignal).getParameter<std::string>("VariableName").c_str());
     if( VariableName.Contains("[0]") ) VariableName.ReplaceAll("[0]","");      
     cSignal_vsPT->SaveAs(std::string(outputFileDirectory+"/"+VariableName+"_signal_vsPT.png").c_str(),"png");
     cSignal_vsPT->SaveAs(std::string(outputFileDirectory+"/"+VariableName+"_signal_vsPT.pdf").c_str(),"pdf");
     cSignal_vsPT->SaveAs(std::string(outputFileDirectory+"/"+VariableName+"_signal_vsPT.root").c_str(),"root");

   }

   ///////////////////
   
   TCanvas *cBackground_vsPT = new TCanvas("cBackground_vsPT","",180,52,550,550);

   cBackground_vsPT->SetTicks();
   cBackground_vsPT->SetFillColor(0);
   cBackground_vsPT->SetBorderMode(0);
   cBackground_vsPT->SetBorderSize(2);
   cBackground_vsPT->SetTickx(1);
   cBackground_vsPT->SetTicky(1);
   cBackground_vsPT->SetRightMargin(0.05);
   cBackground_vsPT->SetBottomMargin(0.12);
   cBackground_vsPT->SetFrameBorderMode(0);

   TH2F* frameBackground_vsPT = new TH2F("frameBackground_vsPT","",500,JetPtBinBackground.at(0),JetPtBinBackground.back(),500,0,0.3);
   frameBackground_vsPT->SetLineWidth(2);
   frameBackground_vsPT->SetMarkerStyle(21);
   frameBackground_vsPT->SetMarkerSize(0.3);
   frameBackground_vsPT->GetXaxis()->SetNdivisions(505);
   frameBackground_vsPT->GetYaxis()->SetNdivisions(10906);
   frameBackground_vsPT->GetXaxis()->SetTitle("AK8 jet p_{T} (GeV)");
   frameBackground_vsPT->GetXaxis()->SetLabelOffset(0.012);
   frameBackground_vsPT->GetXaxis()->SetLabelSize(0.038);
   frameBackground_vsPT->GetXaxis()->SetTitleSize(0.05);
   frameBackground_vsPT->GetXaxis()->SetTitleOffset(1.10);
   frameBackground_vsPT->GetYaxis()->SetTitle("fake rate");
   frameBackground_vsPT->GetYaxis()->SetLabelOffset(0.012);
   frameBackground_vsPT->GetYaxis()->SetLabelSize(0.038);
   frameBackground_vsPT->GetYaxis()->SetTitleSize(0.05);
   frameBackground_vsPT->GetYaxis()->SetTitleOffset(1.38);

   TLegend* BackgroundLegend = new TLegend(0.175,0.54,0.53,0.7);
   BackgroundLegend->SetBorderSize(0);
   BackgroundLegend->SetFillColor(0);
   BackgroundLegend->SetFillStyle(0);
   BackgroundLegend->SetTextSize(0.031);
   BackgroundLegend->SetTextFont(42);

   TLatex *   backgroundBanner = new TLatex(0.17,0.82,"QCD Pythia8");
   backgroundBanner->SetNDC();
   backgroundBanner->SetTextAlign(11);
   backgroundBanner->SetTextFont(42);
   backgroundBanner->SetTextSize(0.036);
   backgroundBanner->SetLineWidth(2);

   for( size_t iBackground = 0 ; iBackground<backgroundEfficiencyMass_vsPT.size(); iBackground++){
     cBackground_vsPT->cd();
     frameBackground_vsPT->Draw("");
     tex->Draw("same");
     tex2->Draw("same");
     tex3->Draw("same");
     backgroundBanner->Draw("same");
     jetBanner->Draw("same");      

     backgroundEfficiencyMass_vsPT.at(iBackground)->SetLineColor(kBlack); 
     backgroundEfficiencyMass_vsPT.at(iBackground)->SetMarkerColor(kBlack); 
     backgroundEfficiencyMass_vsPT.at(iBackground)->SetLineWidth(2); 
     backgroundEfficiencyMass_vsPT.at(iBackground)->SetMarkerStyle(20); 
     backgroundEfficiencyMass_vsPT.at(iBackground)->Draw("psame");

     backgroundEfficiencyWtag_vsPT.at(iBackground)->SetLineColor(kRed); 
     backgroundEfficiencyWtag_vsPT.at(iBackground)->SetMarkerColor(kRed); 
     backgroundEfficiencyWtag_vsPT.at(iBackground)->SetLineWidth(2); 
     backgroundEfficiencyWtag_vsPT.at(iBackground)->SetMarkerStyle(22); 
     backgroundEfficiencyWtag_vsPT.at(iBackground)->Draw("psame");

     BackgroundLegend->Clear();
     BackgroundLegend->AddEntry(signalEfficiencyMass_vsPT.at(iBackground),Form("%s #in [%d,%d] GeV",InputMassVariables.at(iBackground).getParameter<std::string>("ReducedName").c_str(),int(MinMassCut),int(MaxMassCut)),"pl");
     BackgroundLegend->AddEntry(signalEfficiencyWtag_vsPT.at(iBackground),Form("%s #in [%d,%d] GeV + %s < %0.1f",InputMassVariables.at(iBackground).getParameter<std::string>("ReducedName").c_str(),int(MinMassCut),int(MaxMassCut),WtagVariable.getParameter<std::string>("ReducedName").c_str(),WtagCut),"pl");
     BackgroundLegend->Draw("same");
          
     TString VariableName = Form("%s",InputMassVariables.at(iBackground).getParameter<std::string>("VariableName").c_str());
     if( VariableName.Contains("[0]") ) VariableName.ReplaceAll("[0]","");      
     cBackground_vsPT->SaveAs(std::string(outputFileDirectory+"/"+VariableName+"_background_vsPT.png").c_str(),"png");
     cBackground_vsPT->SaveAs(std::string(outputFileDirectory+"/"+VariableName+"_background_vsPT.pdf").c_str(),"pdf");
     cBackground_vsPT->SaveAs(std::string(outputFileDirectory+"/"+VariableName+"_background_vsPT.root").c_str(),"root");

   }

   //////////////////// vs PU
   TCanvas *cSignal_vsPU = new TCanvas("cSignal_vsPU","",180,52,550,550);
   
   cSignal_vsPU->SetTicks();
   cSignal_vsPU->SetFillColor(0);
   cSignal_vsPU->SetBorderMode(0);
   cSignal_vsPU->SetBorderSize(2);
   cSignal_vsPU->SetTickx(1);
   cSignal_vsPU->SetTicky(1);
   cSignal_vsPU->SetRightMargin(0.05);
   cSignal_vsPU->SetBottomMargin(0.12);
   cSignal_vsPU->SetFrameBorderMode(0);

   TH2F* frameSignal_vsPU = new TH2F("frameSignal_vsPU","",500,JetPUBin.at(0),JetPUBin.back(),500,0,1);
   frameSignal_vsPU->SetLineWidth(2);
   frameSignal_vsPU->SetMarkerStyle(21);
   frameSignal_vsPU->SetMarkerSize(0.3);
   frameSignal_vsPU->GetXaxis()->SetNdivisions(405);
   frameSignal_vsPU->GetYaxis()->SetNdivisions(505);
   frameSignal_vsPU->GetXaxis()->SetTitle("N_{PU}");
   frameSignal_vsPU->GetXaxis()->SetLabelOffset(0.012);
   frameSignal_vsPU->GetXaxis()->SetLabelSize(0.038);
   frameSignal_vsPU->GetXaxis()->SetTitleSize(0.05);
   frameSignal_vsPU->GetXaxis()->SetTitleOffset(1.10);
   frameSignal_vsPU->GetYaxis()->SetTitle("efficiency");
   frameSignal_vsPU->GetYaxis()->SetLabelOffset(0.012);
   frameSignal_vsPU->GetYaxis()->SetLabelSize(0.038);
   frameSignal_vsPU->GetYaxis()->SetTitleSize(0.05);
   frameSignal_vsPU->GetYaxis()->SetTitleOffset(1.18);

   for( size_t iSignal = 0 ; iSignal<signalEfficiencyMass_vsPU.size(); iSignal++){
     cSignal_vsPU->cd();
     frameSignal_vsPU->Draw();
     tex->Draw("same");
     tex2->Draw("same");
     tex3->Draw("same");
     signalBanner->Draw("same");
     jetBanner->Draw("same");

     signalEfficiencyMass_vsPU.at(iSignal)->SetLineColor(kBlack); 
     signalEfficiencyMass_vsPU.at(iSignal)->SetMarkerColor(kBlack); 
     signalEfficiencyMass_vsPU.at(iSignal)->SetLineWidth(2); 
     signalEfficiencyMass_vsPU.at(iSignal)->SetMarkerStyle(20); 
     signalEfficiencyMass_vsPU.at(iSignal)->Draw("psame");

     signalEfficiencyWtag_vsPU.at(iSignal)->SetLineColor(kRed); 
     signalEfficiencyWtag_vsPU.at(iSignal)->SetMarkerColor(kRed); 
     signalEfficiencyWtag_vsPU.at(iSignal)->SetLineWidth(2); 
     signalEfficiencyWtag_vsPU.at(iSignal)->SetMarkerStyle(22); 
     signalEfficiencyWtag_vsPU.at(iSignal)->Draw("psame");

     SignalLegend->Clear();
     SignalLegend->AddEntry(signalEfficiencyMass_vsPU.at(iSignal),Form("%s #in [%d,%d] GeV",InputMassVariables.at(iSignal).getParameter<std::string>("ReducedName").c_str(),int(MinMassCut),int(MaxMassCut)),"pl");
     SignalLegend->AddEntry(signalEfficiencyWtag_vsPU.at(iSignal),Form("%s #in [%d,%d] GeV + %s < %0.1f",InputMassVariables.at(iSignal).getParameter<std::string>("ReducedName").c_str(),int(MinMassCut),int(MaxMassCut),WtagVariable.getParameter<std::string>("ReducedName").c_str(),WtagCut),"pl");
     SignalLegend->Draw("same");
          
     TString VariableName = Form("%s",InputMassVariables.at(iSignal).getParameter<std::string>("VariableName").c_str());
     if( VariableName.Contains("[0]") ) VariableName.ReplaceAll("[0]","");      
     cSignal_vsPU->SaveAs(std::string(outputFileDirectory+"/"+VariableName+"_signal_vsPU.png").c_str(),"png");
     cSignal_vsPU->SaveAs(std::string(outputFileDirectory+"/"+VariableName+"_signal_vsPU.pdf").c_str(),"pdf");
     cSignal_vsPU->SaveAs(std::string(outputFileDirectory+"/"+VariableName+"_signal_vsPU.root").c_str(),"root");

   }


   TCanvas *cBackground_vsPU = new TCanvas("cBackground_vsPU","",180,52,550,550);

   cBackground_vsPU->SetTicks();
   cBackground_vsPU->SetFillColor(0);
   cBackground_vsPU->SetBorderMode(0);
   cBackground_vsPU->SetBorderSize(2);
   cBackground_vsPU->SetTickx(1);
   cBackground_vsPU->SetTicky(1);
   cBackground_vsPU->SetRightMargin(0.05);
   cBackground_vsPU->SetBottomMargin(0.12);
   cBackground_vsPU->SetFrameBorderMode(0);

   TH2F* frameBackground_vsPU = new TH2F("frameBackground_vsPU","",500,JetPUBin.at(0),JetPUBin.back(),500,0,0.3);
   frameBackground_vsPU->SetLineWidth(2);
   frameBackground_vsPU->SetMarkerStyle(21);
   frameBackground_vsPU->SetMarkerSize(0.3);
   frameBackground_vsPU->GetXaxis()->SetNdivisions(405);
   frameBackground_vsPU->GetYaxis()->SetNdivisions(505);
   frameBackground_vsPU->GetXaxis()->SetTitle("N_{PU}");
   frameBackground_vsPU->GetXaxis()->SetLabelOffset(0.012);
   frameBackground_vsPU->GetXaxis()->SetLabelSize(0.038);
   frameBackground_vsPU->GetXaxis()->SetTitleSize(0.05);
   frameBackground_vsPU->GetXaxis()->SetTitleOffset(1.10);
   frameBackground_vsPU->GetYaxis()->SetTitle("fake rate");
   frameBackground_vsPU->GetYaxis()->SetLabelOffset(0.012);
   frameBackground_vsPU->GetYaxis()->SetLabelSize(0.038);
   frameBackground_vsPU->GetYaxis()->SetTitleSize(0.05);
   frameBackground_vsPU->GetYaxis()->SetTitleOffset(1.18);

   for( size_t iBackground = 0 ; iBackground<backgroundEfficiencyMass_vsPU.size(); iBackground++){
     cBackground_vsPU->cd();
     frameBackground_vsPU->Draw("");
     tex->Draw("same");
     tex2->Draw("same");
     tex3->Draw("same");
     backgroundBanner->Draw("same");
     jetBanner->Draw("same");      

     backgroundEfficiencyMass_vsPU.at(iBackground)->SetLineColor(kBlack); 
     backgroundEfficiencyMass_vsPU.at(iBackground)->SetMarkerColor(kBlack); 
     backgroundEfficiencyMass_vsPU.at(iBackground)->SetLineWidth(2); 
     backgroundEfficiencyMass_vsPU.at(iBackground)->SetMarkerStyle(20); 
     backgroundEfficiencyMass_vsPU.at(iBackground)->Draw("psame");

     backgroundEfficiencyWtag_vsPU.at(iBackground)->SetLineColor(kRed); 
     backgroundEfficiencyWtag_vsPU.at(iBackground)->SetMarkerColor(kRed); 
     backgroundEfficiencyWtag_vsPU.at(iBackground)->SetLineWidth(2); 
     backgroundEfficiencyWtag_vsPU.at(iBackground)->SetMarkerStyle(22); 
     backgroundEfficiencyWtag_vsPU.at(iBackground)->Draw("psame");

     BackgroundLegend->Clear();
     BackgroundLegend->AddEntry(signalEfficiencyMass_vsPU.at(iBackground),Form("%s #in [%d,%d] GeV",InputMassVariables.at(iBackground).getParameter<std::string>("ReducedName").c_str(),int(MinMassCut),int(MaxMassCut)),"pl");
     BackgroundLegend->AddEntry(signalEfficiencyWtag_vsPU.at(iBackground),Form("%s #in [%d,%d] GeV + %s < %0.1f",InputMassVariables.at(iBackground).getParameter<std::string>("ReducedName").c_str(),int(MinMassCut),int(MaxMassCut),WtagVariable.getParameter<std::string>("ReducedName").c_str(),WtagCut),"pl");
     BackgroundLegend->Draw("same");
          
     TString VariableName = Form("%s",InputMassVariables.at(iBackground).getParameter<std::string>("VariableName").c_str());
     if( VariableName.Contains("[0]") ) VariableName.ReplaceAll("[0]","");      
     cBackground_vsPU->SaveAs(std::string(outputFileDirectory+"/"+VariableName+"_background_vsPU.png").c_str(),"png");
     cBackground_vsPU->SaveAs(std::string(outputFileDirectory+"/"+VariableName+"_background_vsPU.pdf").c_str(),"pdf");
     cBackground_vsPU->SaveAs(std::string(outputFileDirectory+"/"+VariableName+"_background_vsPU.root").c_str(),"root");

   }

   //////////////////// vs PU
   TCanvas *cSignal_vsEta = new TCanvas("cSignal_vsEta","",180,52,550,550);
   
   cSignal_vsEta->SetTicks();
   cSignal_vsEta->SetFillColor(0);
   cSignal_vsEta->SetBorderMode(0);
   cSignal_vsEta->SetBorderSize(2);
   cSignal_vsEta->SetTickx(1);
   cSignal_vsEta->SetTicky(1);
   cSignal_vsEta->SetRightMargin(0.05);
   cSignal_vsEta->SetBottomMargin(0.12);
   cSignal_vsEta->SetFrameBorderMode(0);

   TH2F* frameSignal_vsEta = new TH2F("frameSignal_vsEta","",500,JetEtaBin.at(0),JetEtaBin.back(),500,0,1);
   frameSignal_vsEta->SetLineWidth(2);
   frameSignal_vsEta->SetMarkerStyle(21);
   frameSignal_vsEta->SetMarkerSize(0.3);
   frameSignal_vsEta->GetXaxis()->SetNdivisions(405);
   frameSignal_vsEta->GetYaxis()->SetNdivisions(505);
   frameSignal_vsEta->GetXaxis()->SetTitle("AK8 jet #eta");
   frameSignal_vsEta->GetXaxis()->SetLabelOffset(0.012);
   frameSignal_vsEta->GetXaxis()->SetLabelSize(0.038);
   frameSignal_vsEta->GetXaxis()->SetTitleSize(0.05);
   frameSignal_vsEta->GetXaxis()->SetTitleOffset(1.10);
   frameSignal_vsEta->GetYaxis()->SetTitle("efficiency");
   frameSignal_vsEta->GetYaxis()->SetLabelOffset(0.012);
   frameSignal_vsEta->GetYaxis()->SetLabelSize(0.038);
   frameSignal_vsEta->GetYaxis()->SetTitleSize(0.05);
   frameSignal_vsEta->GetYaxis()->SetTitleOffset(1.18);

   for( size_t iSignal = 0 ; iSignal<signalEfficiencyMass_vsEta.size(); iSignal++){
     cSignal_vsEta->cd();
     frameSignal_vsEta->Draw();
     tex->Draw("same");
     tex2->Draw("same");
     tex3->Draw("same");
     signalBanner->Draw("same");
     jetBanner->Draw("same");

     signalEfficiencyMass_vsEta.at(iSignal)->SetLineColor(kBlack); 
     signalEfficiencyMass_vsEta.at(iSignal)->SetMarkerColor(kBlack); 
     signalEfficiencyMass_vsEta.at(iSignal)->SetLineWidth(2); 
     signalEfficiencyMass_vsEta.at(iSignal)->SetMarkerStyle(20); 
     signalEfficiencyMass_vsEta.at(iSignal)->Draw("psame");

     signalEfficiencyWtag_vsEta.at(iSignal)->SetLineColor(kRed); 
     signalEfficiencyWtag_vsEta.at(iSignal)->SetMarkerColor(kRed); 
     signalEfficiencyWtag_vsEta.at(iSignal)->SetLineWidth(2); 
     signalEfficiencyWtag_vsEta.at(iSignal)->SetMarkerStyle(22); 
     signalEfficiencyWtag_vsEta.at(iSignal)->Draw("psame");

     SignalLegend->Clear();
     SignalLegend->AddEntry(signalEfficiencyMass_vsEta.at(iSignal),Form("%s #in [%d,%d] GeV",InputMassVariables.at(iSignal).getParameter<std::string>("ReducedName").c_str(),int(MinMassCut),int(MaxMassCut)),"pl");
     SignalLegend->AddEntry(signalEfficiencyWtag_vsEta.at(iSignal),Form("%s #in [%d,%d] GeV + %s < %0.1f",InputMassVariables.at(iSignal).getParameter<std::string>("ReducedName").c_str(),int(MinMassCut),int(MaxMassCut),WtagVariable.getParameter<std::string>("ReducedName").c_str(),WtagCut),"pl");
     SignalLegend->Draw("same");
          
     TString VariableName = Form("%s",InputMassVariables.at(iSignal).getParameter<std::string>("VariableName").c_str());
     if( VariableName.Contains("[0]") ) VariableName.ReplaceAll("[0]","");      
     cSignal_vsEta->SaveAs(std::string(outputFileDirectory+"/"+VariableName+"_signal_vsEta.png").c_str(),"png");
     cSignal_vsEta->SaveAs(std::string(outputFileDirectory+"/"+VariableName+"_signal_vsEta.pdf").c_str(),"pdf");
     cSignal_vsEta->SaveAs(std::string(outputFileDirectory+"/"+VariableName+"_signal_vsEta.root").c_str(),"root");

   }


   TCanvas *cBackground_vsEta = new TCanvas("cBackground_vsEta","",180,52,550,550);

   cBackground_vsEta->SetTicks();
   cBackground_vsEta->SetFillColor(0);
   cBackground_vsEta->SetBorderMode(0);
   cBackground_vsEta->SetBorderSize(2);
   cBackground_vsEta->SetTickx(1);
   cBackground_vsEta->SetTicky(1);
   cBackground_vsEta->SetRightMargin(0.05);
   cBackground_vsEta->SetBottomMargin(0.12);
   cBackground_vsEta->SetFrameBorderMode(0);

   TH2F* frameBackground_vsEta = new TH2F("frameBackground_vsEta","",500,JetEtaBin.at(0),JetEtaBin.back(),500,0,0.3);
   frameBackground_vsEta->SetLineWidth(2);
   frameBackground_vsEta->SetMarkerStyle(21);
   frameBackground_vsEta->SetMarkerSize(0.3);
   frameBackground_vsEta->GetXaxis()->SetNdivisions(405);
   frameBackground_vsEta->GetYaxis()->SetNdivisions(505);
   frameBackground_vsEta->GetXaxis()->SetTitle("AK8 jet #eta");
   frameBackground_vsEta->GetXaxis()->SetLabelOffset(0.012);
   frameBackground_vsEta->GetXaxis()->SetLabelSize(0.038);
   frameBackground_vsEta->GetXaxis()->SetTitleSize(0.05);
   frameBackground_vsEta->GetXaxis()->SetTitleOffset(1.10);
   frameBackground_vsEta->GetYaxis()->SetTitle("fake rate");
   frameBackground_vsEta->GetYaxis()->SetLabelOffset(0.012);
   frameBackground_vsEta->GetYaxis()->SetLabelSize(0.038);
   frameBackground_vsEta->GetYaxis()->SetTitleSize(0.05);
   frameBackground_vsEta->GetYaxis()->SetTitleOffset(1.18);

   for( size_t iBackground = 0 ; iBackground<backgroundEfficiencyMass_vsEta.size(); iBackground++){
     cBackground_vsEta->cd();
     frameBackground_vsEta->Draw("");
     tex->Draw("same");
     tex2->Draw("same");
     tex3->Draw("same");
     backgroundBanner->Draw("same");
     jetBanner->Draw("same");      

     backgroundEfficiencyMass_vsEta.at(iBackground)->SetLineColor(kBlack); 
     backgroundEfficiencyMass_vsEta.at(iBackground)->SetMarkerColor(kBlack); 
     backgroundEfficiencyMass_vsEta.at(iBackground)->SetLineWidth(2); 
     backgroundEfficiencyMass_vsEta.at(iBackground)->SetMarkerStyle(20); 
     backgroundEfficiencyMass_vsEta.at(iBackground)->Draw("psame");

     backgroundEfficiencyWtag_vsEta.at(iBackground)->SetLineColor(kRed); 
     backgroundEfficiencyWtag_vsEta.at(iBackground)->SetMarkerColor(kRed); 
     backgroundEfficiencyWtag_vsEta.at(iBackground)->SetLineWidth(2); 
     backgroundEfficiencyWtag_vsEta.at(iBackground)->SetMarkerStyle(22); 
     backgroundEfficiencyWtag_vsEta.at(iBackground)->Draw("psame");

     BackgroundLegend->Clear();
     BackgroundLegend->AddEntry(signalEfficiencyMass_vsEta.at(iBackground),Form("%s #in [%d,%d] GeV",InputMassVariables.at(iBackground).getParameter<std::string>("ReducedName").c_str(),int(MinMassCut),int(MaxMassCut)),"pl");
     BackgroundLegend->AddEntry(signalEfficiencyWtag_vsEta.at(iBackground),Form("%s #in [%d,%d] GeV + %s < %0.1f",InputMassVariables.at(iBackground).getParameter<std::string>("ReducedName").c_str(),int(MinMassCut),int(MaxMassCut),WtagVariable.getParameter<std::string>("ReducedName").c_str(),WtagCut),"pl");
     BackgroundLegend->Draw("same");
          
     TString VariableName = Form("%s",InputMassVariables.at(iBackground).getParameter<std::string>("VariableName").c_str());
     if( VariableName.Contains("[0]") ) VariableName.ReplaceAll("[0]","");      
     cBackground_vsEta->SaveAs(std::string(outputFileDirectory+"/"+VariableName+"_background_vsEta.png").c_str(),"png");
     cBackground_vsEta->SaveAs(std::string(outputFileDirectory+"/"+VariableName+"_background_vsEta.pdf").c_str(),"pdf");
     cBackground_vsEta->SaveAs(std::string(outputFileDirectory+"/"+VariableName+"_background_vsEta.root").c_str(),"root");

   }

   return 0;
}
 
