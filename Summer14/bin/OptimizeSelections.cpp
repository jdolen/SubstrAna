#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <functional>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TMath.h"
#include "TCanvas.h"

#include "TSystem.h"
#include "TROOT.h"

#include "../include/TrainingMVAClass.h"

#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"


/// Main Programme                                                                                                                                                                      
int main (int argc, char** argv){

  if (argc < 2){
    std::cerr << ">>> Usage:   " << argv[1] << "   cfg file" <<  std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Load TTree Lybrary                                                                                                                                                                  
  gSystem->Load("libTree.so");

  TMVA::Tools::Instance();

  // parse config file parameter                                                                                                                                                          
  std::string configFileName = argv[1];
  boost::shared_ptr<edm::ParameterSet> parameterSet = edm::readConfig(configFileName);
  edm::ParameterSet Options  = parameterSet -> getParameter<edm::ParameterSet>("Options");


  std::vector<edm::ParameterSet> InputBackgroundParam ;
  if(Options.existsAs<std::vector<edm::ParameterSet>>("InputBackgroundParam"))
     InputBackgroundParam = Options.getParameter<std::vector<edm::ParameterSet>>("InputBackgroundParam");
  else{ std::cout<<" Exit from code, no input set found for background"<<std::endl; std::exit(EXIT_FAILURE); }

  std::vector<edm::ParameterSet> InputSignalParam ;
  if(Options.existsAs<std::vector<edm::ParameterSet>>("InputSignalParam"))
     InputSignalParam = Options.getParameter<std::vector<edm::ParameterSet>>("InputSignalParam");
  else{ std::cout<<" Exit from code, no input set found for background"<<std::endl; std::exit(EXIT_FAILURE); }


  std::vector<std::string> InputVariableList;
  if(Options.existsAs<std::vector<std::string>>("InputVariableList"))  
     InputVariableList = Options.getParameter<std::vector<std::string>>("InputVariableList");
  else{ std::cout<<" Exit from code, no input variable file list found "<<std::endl; std::exit(EXIT_FAILURE); }

  std::cout << std::endl;
  std::cout << " >>>>> Option::InputVariableList size = " << InputVariableList.size() << std::endl;

  for (unsigned int iVar = 0; iVar < InputVariableList.size(); iVar++){
       std::cout << " " << InputVariableList.at(iVar) << ", ";
  }
  std::cout << std::endl;

  std::vector<std::string> InputSpectatorList;
  if(Options.existsAs<std::vector<std::string>>("InputSpectatorList"))  
     InputSpectatorList = Options.getParameter<std::vector<std::string>>("InputSpectatorList");
  else{ std::cout<<" Exit from code, no input variable file list found "<<std::endl; std::exit(EXIT_FAILURE); }

  std::cout << std::endl;
  std::cout << " >>>>> Option::InputSpectatorList size = " << InputSpectatorList.size() << std::endl;

  for (unsigned int iVar = 0; iVar < InputSpectatorList.size(); iVar++){
       std::cout << " " << InputSpectatorList.at(iVar) << ", ";
  }
  std::cout << std::endl;

  std::string TreeName ; 
  if(Options.existsAs<std::string>("TreeName"))  
   TreeName = Options.getParameter<std::string>("TreeName");
  else{ std::cout<<" Exit from code, no TreeName found "<<std::endl; std::exit(EXIT_FAILURE); }
 
  std::string Label;
  if(Options.existsAs<std::string>("Label"))  
   Label = Options.getParameter<std::string>("Label");
  else{ std::cout<<" Label set to Test "<<std::endl; Label = "Test"; }
  
  std::string LeptonType;
  if(Options.existsAs<std::string>("LeptonType"))
    LeptonType  = Options.getParameter<std::string>("LeptonType");
  else{ std::cout<<" Lepton type set by default to muon "<<std::endl; LeptonType = "Muon"; }


  std::cout<<std::endl;
  std::cout<<" Input TreeName       = "<<TreeName<<std::endl;
  std::cout<<" Input Label          = "<<Label<<std::endl;
  std::cout<<" Input LeptonType     = "<<LeptonType<<std::endl;
  std::cout<<std::endl;

  bool isPrintResultwithTMVA = false ;
  if(Options.existsAs<bool>("isPrintResultwithTMVA"))
    isPrintResultwithTMVA  = Options.getParameter<bool>("isPrintResultwithTMVA");
  else std::cout<<" isPrintResultwithTMVA --> set to false "<<std::endl;

  bool isTrainEachVariable = false ;
  if(Options.existsAs<bool>("TrainEachVariable"))
    isTrainEachVariable  = Options.getParameter<bool>("TrainEachVariable");
  else std::cout<<" TrainEachVariable --> set to false "<<std::endl;

  std::string SignalggHName ;
  if(Options.existsAs<std::string>("SignalggHName"))
    SignalggHName  = Options.getParameter<std::string>("SignalggHName");
  else std::cout<<" SignalggHName --> not found --> empty string "<<std::endl;

  std::cout<<" Option Signal ggH Name = "<<SignalggHName<<std::endl;
  std::cout<<std::endl;

  std::string SignalqqHName ;
  if(Options.existsAs<std::string>("SignalqqHName"))
    SignalqqHName  = Options.getParameter<std::string>("SignalqqHName");
  else std::cout<<" SignalqqHName --> not found --> empty string "<<std::endl;

  std::cout<<" Option Signal qqH Name = "<<SignalqqHName<<std::endl;
  std::cout<<std::endl;

  int useTypeOfSignal = 0;
  if(Options.existsAs<int>("useTypeOfSignal"))
    useTypeOfSignal  = Options.getParameter<int>("useTypeOfSignal");
  else std::cout<<" type of signal to be considered not found --> set to 0 "<<std::endl;

  std::cout<<" Option  useTypeOfSignal = "<<useTypeOfSignal<<std::endl;
  std::cout<<std::endl;

  std::string EventWeight ;
  if(Options.existsAs<std::string>("EventWeight"))
    EventWeight  = Options.getParameter<std::string>("EventWeight");
  else std::cout<<" EventWeight --> not found --> empty string "<<std::endl;

  std::cout<<" Option Event Weight = "<<EventWeight<<std::endl;
  std::cout<<std::endl;


  std::string PreselectionCutType ;
  if(Options.existsAs<std::string>("PreselectionCutType"))
    PreselectionCutType  = Options.getParameter<std::string>("PreselectionCutType");
  else{ std::cout<<" PreselectionCutType --> not found --> empty string "<<std::endl; PreselectionCutType = "none"; }


  std::cout<<" Option Preselection Cut = "<<PreselectionCutType<<std::endl;
  std::cout<<std::endl;

  std::vector<std::string> UseMethodName;
  if(Options.existsAs<std::vector<std::string> >("UseMethodName"))
    UseMethodName  = Options.getParameter<std::vector<std::string>>("UseMethodName");
  else{ std::cout<<" UseMethodName --> not found --> exit from the code "<<std::endl; std::exit(EXIT_FAILURE) ; }

  std::cout << std::endl;
  std::cout << " >>>>> Option::UseMethodName size = " << UseMethodName.size() << std::endl;

  for (unsigned int iCat = 0; iCat < UseMethodName.size(); iCat++){
    std::cout << " " << UseMethodName.at(iCat) << ", ";
  }
  std::cout << std::endl;

  std::vector<double> JetPtBinOfTraining;
  if(Options.existsAs<std::vector<double>>("JetPtBinOfTraining"))
     JetPtBinOfTraining  = Options.getParameter<std::vector<double>>("JetPtBinOfTraining");
  else{ std::cout<<" JetPtBinOfTraining --> not found --> set 1 bin [0,2000] "<<std::endl;
        JetPtBinOfTraining.push_back(0);
        JetPtBinOfTraining.push_back(2000);
  }

  std::cout << std::endl;
  std::cout << " >>>>> Option::JetPtBinOfTraining size = " << JetPtBinOfTraining.size()/2 +1 << std::endl;
  
  for (unsigned int iCat = 0; iCat+1 < JetPtBinOfTraining.size(); iCat++){
    std::cout << " bin min =  " << JetPtBinOfTraining.at(iCat) << " ;  bin max =  "<<JetPtBinOfTraining.at(iCat+1) <<std::endl;
  }
  std::cout << std::endl;

  std::vector<double> PileUpBinOfTraining;
  if(Options.existsAs<std::vector<double>>("PileUpBinOfTraining"))
     PileUpBinOfTraining  = Options.getParameter<std::vector<double>>("PileUpBinOfTraining");
  else{ std::cout<<" PileUpBinOfTraining --> not found --> set 1 bin [0,100] "<<std::endl;
        PileUpBinOfTraining.push_back(0);
        PileUpBinOfTraining.push_back(2000);
  }

  std::cout << std::endl;
  std::cout << " >>>>> Option::PileUpBinOfTraining size = " << PileUpBinOfTraining.size()/2 +1 << std::endl;
  
  for (unsigned int iCat = 0; iCat+1 < PileUpBinOfTraining.size(); iCat++){
    std::cout << " bin min =  " << PileUpBinOfTraining.at(iCat) << " ;  bin max =  "<<PileUpBinOfTraining.at(iCat+1) <<std::endl;
  }
  std::cout << std::endl;

  
  std::string outputFileDirectory ;
  if(Options.existsAs<std::string>("outputFileDirectory"))
   outputFileDirectory = Options.getParameter<std::string>("outputFileDirectory"); 
  else{std::cout<<" Default output directory "<<std::endl; 
       outputFileDirectory = "output/TMVATrainingResult";
  }

  std::string outputFileName ;
  if(Options.existsAs<std::string>("outputFileName"))
   outputFileName = Options.getParameter<std::string>("outputFileName"); 
  else{std::cout<<" Default output file name "<<std::endl; 
       outputFileName = "TMVAoutput";
  }
  
  std::cout<<"                      "<<std::endl;
  std::cout<<" Outout Directory     "<<outputFileDirectory<<std::endl;
  std::cout<<" Input Sample List    "<<outputFileName<<std::endl;
  std::cout<<"                      "<<std::endl;
     
  
  // Import Sample and signal - background collections  
  std::vector <TFile*> signalFileList;
  std::vector <TFile*> backgroundFileList;
  std::vector <TTree*> signalTreeList;
  std::vector <TTree*> backgroundTreeList;

  std::cout<<" Building Tree List for Signal And Background  "<<std::endl;
  std::cout<<std::endl;

  std::vector<int> badBackgroundFiles ;
  std::vector<int> badSignalFiles ;

  std::vector<edm::ParameterSet>::const_iterator itBackground = InputBackgroundParam.begin();
  for( ; itBackground != InputBackgroundParam.end() ; ++itBackground){

    TString NameFile     = Form("%s",(*itBackground).getParameter<std::string>("inputFileName").c_str());
    std::cout<<" Input File Bkg: "<< NameFile.Data()<<std::endl;               
    backgroundFileList.push_back (TFile::Open(NameFile.Data(),"READ") );
    if(not backgroundFileList.back() or backgroundFileList.back() == NULL){ badBackgroundFiles.push_back(backgroundFileList.size()-1); continue; }
    backgroundTreeList.push_back( (TTree*) backgroundFileList.back()->Get(TreeName.c_str())); 	
    if(backgroundTreeList.back() == 0 or backgroundTreeList.back() == NULL) backgroundTreeList.erase(backgroundTreeList.end());            
  }

  backgroundFileList.clear();

  std::vector<edm::ParameterSet>::const_iterator itSignal = InputSignalParam.begin();
  for( ; itSignal != InputSignalParam.end() ; ++itSignal){
    TString NameFile     = Form("%s",(*itSignal).getParameter<std::string>("inputFileName").c_str());
    std::cout<<" Input File Signal: "<< NameFile.Data()<<std::endl;               
    if((*itSignal).getParameter<std::string>("ReducedName") == SignalqqHName and (useTypeOfSignal == 0 or useTypeOfSignal == 1) ){
          signalFileList.push_back (TFile::Open(NameFile.Data(),"READ") );
          if(not signalFileList.back() or signalFileList.back() == NULL){ badSignalFiles.push_back(signalFileList.size()-1); continue; }
	  signalTreeList.push_back( (TTree*) signalFileList.back()->Get(TreeName.c_str())); 
          if(signalTreeList.back() == 0 or signalTreeList.back() == NULL) signalTreeList.erase(signalTreeList.end());            
   }
    else if((*itSignal).getParameter<std::string>("ReducedName") == SignalggHName and (useTypeOfSignal == 0 or useTypeOfSignal == 2)){
          signalFileList.push_back ( TFile::Open(NameFile.Data(),"READ") );
          if(not signalFileList.back() or signalFileList.back() == NULL) { badSignalFiles.push_back(signalFileList.size()-1); continue; }
	  signalTreeList.push_back( (TTree*) signalFileList.back()->Get(TreeName.c_str()));
          if(signalTreeList.back() == 0 or signalTreeList.back() == NULL) signalTreeList.erase(signalTreeList.end());            
    }	
  }

  signalFileList.clear();

  std::cout<<std::endl;
  
  // Book MVA Training Object --> one for each pT bin 
  std::vector<TrainingMVAClass*> WWTrainingVector ;
  // scale factor for W+jet 
  double scaleFactorWjet = 1. ;
  if( argc==3 ) scaleFactorWjet =  std::atof(argv[2]);

  // Loop in order to start the training    
  for(size_t pTBin = 0; pTBin+1 < JetPtBinOfTraining.size() ; pTBin++){

   std::cout<<" pT bin of Training: Min = "<<JetPtBinOfTraining.at(pTBin)<<" Max = "<<JetPtBinOfTraining.at(pTBin+1)<<std::endl;

   for(size_t puBin = 0; puBin+1 < PileUpBinOfTraining.size() ; puBin++){

    std::cout<<" pu bin of Training: Min = "<<PileUpBinOfTraining.at(puBin)<<" Max = "<<PileUpBinOfTraining.at(puBin+1)<<std::endl;
 
    TString Label_ = Form("%s_PTBin_%d_%d_PU_%d_%d",Label.c_str(),int(JetPtBinOfTraining.at(pTBin)),int(JetPtBinOfTraining.at(pTBin+1)),int(PileUpBinOfTraining.at(puBin)),int(PileUpBinOfTraining.at(puBin+1))) ;

    std::string tempLabel; tempLabel = Label_ ;

    if(isTrainEachVariable == true){
     
     for(unsigned int iVariable = 0 ; iVariable < InputVariableList.size() ; iVariable++){
       WWTrainingVector.push_back(new TrainingMVAClass(signalTreeList, backgroundTreeList, TreeName, outputFileDirectory, outputFileName+"_"+InputVariableList.at(iVariable), tempLabel,":TransformationsI,N"));
       // Set Input and Spectator Variables
       std::cout<<std::endl;
       std::cout<<" Start to set input variable  "<<InputVariableList.at(iVariable)<<std::endl;
       std::cout<<std::endl;
 
       WWTrainingVector.back()->AddTrainingVariables(InputVariableList.at(iVariable), InputSpectatorList);

       // Set Global Weight and signal + background Tree for MVA Training
       std::vector<double> signalGlobalWeight (signalTreeList.size(),0.);
       std::vector<double> backgroundGlobalWeight (backgroundTreeList.size(),0.);

       int isSignal = 0;
       int isBackground = 0;

       std::cout<<" Building Global Event Weight  + Add Trees "<<std::endl;
       std::cout<<std::endl; 
       
       std::vector<edm::ParameterSet>::const_iterator itBackground = InputBackgroundParam.begin();
       for( ; itBackground != InputBackgroundParam.end() ; itBackground++){ 
	 if(std::find(badBackgroundFiles.begin(),badBackgroundFiles.end(),isBackground)!=badBackgroundFiles.end()) continue;  
	 if((*itBackground).getParameter<std::string>("ReducedName") == "W+Jets")
	   backgroundGlobalWeight.at(isBackground) = (*itBackground).getParameter<double>("CrossSection")/(double((*itBackground).getParameter<int>("NumberEntriesBefore")))*scaleFactorWjet;
          else backgroundGlobalWeight.at(isBackground) = ((*itBackground).getParameter<double>("CrossSection")/double((*itBackground).getParameter<int>("NumberEntriesBefore")));
          isBackground ++;    
	}

       std::vector<edm::ParameterSet>::const_iterator itSignal = InputSignalParam.begin();
       for( ; itSignal != InputSignalParam.end() ; itSignal++){ 
	 if(std::find(badSignalFiles.begin(),badSignalFiles.end(),isSignal)!=badBackgroundFiles.end()) continue;  
	 if((*itSignal).getParameter<std::string>("ReducedName") == SignalqqHName and (useTypeOfSignal == 0 or useTypeOfSignal == 1)) {       
	   signalGlobalWeight.at(isSignal) = ((*itSignal).getParameter<double>("CrossSection")/double((*itSignal).getParameter<int>("NumberEntriesBefore")));
           isSignal ++; 
         }
         else if((*itSignal).getParameter<std::string>("ReducedName") == SignalggHName and (useTypeOfSignal == 0 or useTypeOfSignal == 2)  ){
	   signalGlobalWeight.at(isSignal) = ((*itSignal).getParameter<double>("CrossSection")/double((*itSignal).getParameter<int>("NumberEntriesBefore")));
           isSignal ++;          
	 }    
       }
       
      WWTrainingVector.back()->BookMVATrees(signalGlobalWeight, backgroundGlobalWeight);  
      // Prepare and Set the MVA Factory
      std::cout<<std::endl;
      std::cout<<" Prepare MVA  "<<std::endl;
      std::cout<<std::endl;
     
      WWTrainingVector.back()->AddPrepareTraining ( LeptonType,PreselectionCutType, EventWeight, EventWeight, &JetPtBinOfTraining, pTBin, &PileUpBinOfTraining, puBin) ;
  
      // Book and Run TMVA Training and testing for the selected methods
      std::cout<<" Loop on the Selected Methods  "<<std::endl;
      std::cout<<std::endl;

      for(size_t iMethod =0; iMethod<UseMethodName.size(); iMethod++){

       // Rectangular Cuts
       if(UseMethodName.at(iMethod) == "CutsMC" )      WWTrainingVector.back()->BookandTrainRectangularCuts("MC",InputVariableList.at(iVariable));
       else if(UseMethodName.at(iMethod) == "CutsGA" ) WWTrainingVector.back()->BookandTrainRectangularCuts("GA",InputVariableList.at(iVariable));
       else if(UseMethodName.at(iMethod) == "CutsSA" ) WWTrainingVector.back()->BookandTrainRectangularCuts("SA",InputVariableList.at(iVariable));
       else { std::cerr<<" Training Method not implemented in the TMVATrainingClass for single variables --> Go to the next one and only rectangluar cuts"<<std::endl; std::cout<<std::endl;}
      }
  
      WWTrainingVector.back()->CloseTrainingAndTesting();

      //Print Output Plots
      std::cout<<" Save Output Image after training and testing ..  "<<std::endl;
      std::cout<<std::endl;
      
      if (isPrintResultwithTMVA) WWTrainingVector.back()->PrintTrainingResults ();
       
     }
       
    }
    
   else{

    WWTrainingVector.push_back(new TrainingMVAClass(signalTreeList, backgroundTreeList, TreeName, outputFileDirectory, outputFileName, tempLabel,":Transformations=I,N:"));

    // Set Input and Spectator Variables
    std::cout<<std::endl;
    std::cout<<" Set Training and Spectator Variables  "<<std::endl;
    std::cout<<std::endl;

    WWTrainingVector.back()->AddTrainingVariables(InputVariableList, InputSpectatorList);
    
    // Set Global Weight and signal + background Tree for MVA Training
    std::vector<double> signalGlobalWeight (signalTreeList.size(),0.);
    std::vector<double> backgroundGlobalWeight (backgroundTreeList.size(),0.);

    int isSignal = 0;
    int isBackground = 0;

    std::cout<<" Building Global Event Weight  + Add Trees "<<std::endl;
    std::cout<<std::endl; 

       std::vector<edm::ParameterSet>::const_iterator itBackground = InputBackgroundParam.begin();
       for( ; itBackground != InputBackgroundParam.end() ; itBackground++){ 
	 if(std::find(badBackgroundFiles.begin(),badBackgroundFiles.end(),isBackground)!=badBackgroundFiles.end()) continue;  
	 if((*itBackground).getParameter<std::string>("ReducedName") == "W+Jets")
	   backgroundGlobalWeight.at(isBackground) = (*itBackground).getParameter<double>("CrossSection")/(double((*itBackground).getParameter<int>("NumberEntriesBefore")))*scaleFactorWjet;
          else backgroundGlobalWeight.at(isBackground) = ((*itBackground).getParameter<double>("CrossSection")/double((*itBackground).getParameter<int>("NumberEntriesBefore")));
          isBackground ++;    
	}

       std::vector<edm::ParameterSet>::const_iterator itSignal = InputSignalParam.begin();
       for( ; itSignal != InputSignalParam.end() ; itSignal++){ 
	 if(std::find(badSignalFiles.begin(),badSignalFiles.end(),isSignal)!=badSignalFiles.end()) continue;  
	 if((*itSignal).getParameter<std::string>("ReducedName") == SignalqqHName and (useTypeOfSignal == 0 or useTypeOfSignal == 1)) {       
	   signalGlobalWeight.at(isSignal) = ((*itSignal).getParameter<double>("CrossSection")/double((*itSignal).getParameter<int>("NumberEntriesBefore")));
           isSignal ++; 
         }
         else if((*itSignal).getParameter<std::string>("ReducedName") == SignalggHName and (useTypeOfSignal == 0 or useTypeOfSignal == 2)  ){
	   signalGlobalWeight.at(isSignal) = ((*itSignal).getParameter<double>("CrossSection")/double((*itSignal).getParameter<int>("NumberEntriesBefore")));
           isSignal ++;          
	 }    
       }
       
    WWTrainingVector.back()->BookMVATrees(signalGlobalWeight, backgroundGlobalWeight);
    
    // Prepare and Set the MVA Factory
    std::cout<<std::endl;
    std::cout<<" Prepare MVA  "<<std::endl;
    std::cout<<std::endl;
     
    WWTrainingVector.back()->AddPrepareTraining ( LeptonType,PreselectionCutType, EventWeight, EventWeight, &JetPtBinOfTraining, pTBin, &PileUpBinOfTraining, puBin) ;
  
    // Book and Run TMVA Training and testing for the selected methods
    std::cout<<" Loop on the Selected Methods  "<<std::endl;
    std::cout<<std::endl;
 
    for(size_t iMethod =0; iMethod<UseMethodName.size(); iMethod++){

     // Rectangular Cuts
     if(UseMethodName.at(iMethod) == "CutsMC" )      WWTrainingVector.back()->BookandTrainRectangularCuts("MC");
     else if(UseMethodName.at(iMethod) == "CutsGA" ) WWTrainingVector.back()->BookandTrainRectangularCuts("GA");
     else if(UseMethodName.at(iMethod) == "CutsSA" ) WWTrainingVector.back()->BookandTrainRectangularCuts("SA");
     else{ 
       //        WWTrainingVector.back()->SetTransformations(":VarTransform=I,N,U,G,P,D,G,D:"); 
       // Likelihood 
       if(UseMethodName.at(iMethod) == "Likelihood")     WWTrainingVector.back()->BookandTrainLikelihood(); 
       else if(UseMethodName.at(iMethod) == "LikelihoodKDE")  WWTrainingVector.back()->BookandTrainLikelihood("LikelihoodKDE"); 
       else if(UseMethodName.at(iMethod) == "PDERS")          WWTrainingVector.back()->BookandTrainLikelihood("PDERS"); 
       else if(UseMethodName.at(iMethod) == "PDEFoam")        WWTrainingVector.back()->BookandTrainLikelihood("PDEFoam"); 
       else if(UseMethodName.at(iMethod) == "PDEFoamBoost")   WWTrainingVector.back()->BookandTrainLikelihood("PDEFoamBoost"); 

       // Fisher Discriminant
       else if(UseMethodName.at(iMethod) == "Fisher")  WWTrainingVector.back()->BookandTrainFisherDiscriminant(); 
    
       // Linear Discriminant
       else if(UseMethodName.at(iMethod) == "LD")      WWTrainingVector.back()->BookandTrainLinearDiscriminant();
    
       // MLP
       else if(UseMethodName.at(iMethod) == "MLP")        WWTrainingVector.back()->BookandTrainMLP();
       else if(UseMethodName.at(iMethod) == "MLPBFG")     WWTrainingVector.back()->BookandTrainMLP(1000,"N+5","sigmoid","BFGS",10,10);
       else if(UseMethodName.at(iMethod) == "CFMlpANN")   WWTrainingVector.back()->BookandTrainCFMlpANN();
       else if(UseMethodName.at(iMethod) == "TMlpANN")    WWTrainingVector.back()->BookandTrainTMlpANN();

       // BDT
       else if(UseMethodName.at(iMethod) == "BDT")     WWTrainingVector.back()->BookandTrainBDT();

       // BDTG
       else if(UseMethodName.at(iMethod) == "BDTG")    WWTrainingVector.back()->BookandTrainBDTG();

       // BDTF
       else if(UseMethodName.at(iMethod) == "BDTF")    WWTrainingVector.back()->BookandTrainBDTF();

       else { std::cerr<<" Training Method not implemented in the TMVATrainingClass >> Go to the next one"<<std::endl; std::cout<<std::endl;}
     }
    }
  
    WWTrainingVector.back()->CloseTrainingAndTesting();

    //Print Output Plots
    std::cout<<" Save Output Image after training and testing ..  "<<std::endl;
    std::cout<<std::endl;

    if (isPrintResultwithTMVA) WWTrainingVector.back()->PrintTrainingResults ();
    }
   }
  }
       
  return 0 ;

}
