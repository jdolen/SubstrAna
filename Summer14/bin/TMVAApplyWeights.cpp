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

#include "../include/TMVAReadingClass.h"

#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

/// Main Programme                                                                                                                                                                 
int main (int argc, char** argv){

  if (argc != 2){
    std::cerr << ">>> Usage:   " << argv[1] << "   cfg file" << std::endl;
    return -1;
  }

  // Load TTree Lybrary
  gSystem->Load("libTree.so");

  TMVA::Tools::Instance();

  // parse config file parameter                                                                                                                                                          
  std::string configFileName = argv[1];
  boost::shared_ptr<edm::ParameterSet> parameterSet = edm::readConfig(configFileName);
  edm::ParameterSet Options  = parameterSet -> getParameter<edm::ParameterSet>("Options");

  std::vector<std::string> InputFileName ;
  if(Options.existsAs<std::vector<std::string>>("InputFileName"))
     InputFileName = Options.getParameter<std::vector<std::string>>("InputFileName");
  else{ std::cout<<" No input directory found for weight TMVA xml file --> exit from the program "<<std::endl; return -1; }

  std::vector<std::string> InputSpectatorList;
  if(Options.existsAs<std::vector<std::string>>("InputSpectatorList"))
    InputSpectatorList = Options.getParameter<std::vector<std::string>>("InputSpectatorList");
  else{ std::cout<<" Exit from code, no input variable file list found "<<std::endl; return -1; }

  std::cout << std::endl;
  std::cout << " >>>>> Option::InputSpectatorList size = " << InputSpectatorList.size() << std::endl;

  for (unsigned int iVar = 0; iVar < InputSpectatorList.size(); iVar++){
    std::cout << " " << InputSpectatorList.at(iVar) << ", ";
  }
  std::cout << std::endl;

  std::string TreeName ;
  if(Options.existsAs<std::string>("TreeName"))
    TreeName = Options.getParameter<std::string>("TreeName");
  else{ std::cout<<" Exit from code, no TreeName found "<<std::endl; return -1; }

  std::string Label;
  if(Options.existsAs<std::string>("Label"))
    Label = Options.getParameter<std::string>("Label");
  else{ std::cout<<" Label set to Test "<<std::endl; Label = "Test"; }

  std::string LeptonType;
  if(Options.existsAs<std::string>("LeptonType"))
    LeptonType  = Options.getParameter<std::string>("LeptonType");
  else{ std::cout<<" Lepton type set by default to muon "<<std::endl; LeptonType = "Jets"; }

  std::cout<<std::endl;
  std::cout<<" Input TreeName       = "<<TreeName<<std::endl;
  std::cout<<" Input Label          = "<<Label<<std::endl;
  std::cout<<" Input LeptonType     = "<<LeptonType<<std::endl;
  std::cout<<std::endl;

  std::vector<edm::ParameterSet> InputVariableFileParam;
  if(Options.existsAs<std::vector<edm::ParameterSet> >("InputVariableFileParam"))
    InputVariableFileParam  = Options.getParameter<std::vector<edm::ParameterSet> >("InputVariableFileParam");
  else{ std::cout<<" No input weight file parameters are provided --> exit from the program "<<std::endl; return -1; }

  std::cout << " >>>>> Option:InputWeightFile size = " << InputVariableFileParam.size() << std::endl;

  for (unsigned int iCat = 0; iCat < InputVariableFileParam.size(); iCat++){
    std::cout<<" input variables "<<std::endl;
    for(unsigned int iVar = 0; iVar < InputVariableFileParam.at(iCat).getParameter<std::vector<std::string> >("InputVariableList").size(); iVar++)
      std::cout<<InputVariableFileParam.at(iCat).getParameter<std::vector<std::string> >("InputVariableList").at(iVar)<<"  "<<std::endl;
    std::cout<<" Method type "<<InputVariableFileParam.at(iCat).getParameter<std::string>("MethodName")<<std::endl;   
    for(unsigned int iWeight = 0; iWeight < InputVariableFileParam.at(iCat).getParameter<std::vector<std::string> >("inputWeightFiles").size(); iWeight++)
      std::cout<<" weight file "<<InputVariableFileParam.at(iCat).getParameter<std::vector<std::string> >("inputWeightFiles").at(iWeight)<<std::endl;    
  }    
   
  std::cout << std::endl;

  std::string PreselectionCutType ;
  if(Options.existsAs<std::string>("PreselectionCutType"))
    PreselectionCutType  = Options.getParameter<std::string>("PreselectionCutType");
  else{ std::cout<<" PreselectionCutType --> not found --> empty string "<<std::endl; PreselectionCutType = "none"; }

  std::cout<<" Option Preselection Cut = "<<PreselectionCutType<<std::endl;
  std::cout<<std::endl;

  // last info to fill the MVA output
  int JetToRead;
  if(Options.existsAs<int>("JetToRead"))
    JetToRead  = Options.getParameter<int>("JetToRead");
  else{ std::cout<<" JetToRead --> not found --> set to 0 "<<std::endl;
        JetToRead = 0 ;
  }

  bool optionCut;
  if(Options.existsAs<bool>("optionCut"))
    optionCut  = Options.getParameter<bool>("optionCut");
  else{ std::cout<<" optionCut --> not found --> set to 0 "<<std::endl;
        optionCut = false ;
  }

  std::cout<<std::endl;
  std::cout<<" jet to read: "<<JetToRead<<std::endl;
  std::cout<<" apply the cut on the readed events "<<optionCut<<std::endl;
  std::cout<<std::endl;

  // Take the input file list to run over --> one file per time is the actual strategy --> more jobs but less cpu time
  std::vector <TFile*> SampleFileList;
  std::vector <TTree*> SampleTreeList;

  std::cout<<std::endl;
  
  std::cout<<" Building Tree List for Signal And Background  "<<std::endl;
  std::cout<<std::endl;

  std::vector<int> badFiles ;

  std::vector<std::string>::const_iterator itFile = InputFileName.begin();
  for( ; itFile != InputFileName.end() ; ++itFile){

    std::cout<<" Input File Bkg: "<<*itFile<<std::endl;
    SampleFileList.push_back (TFile::Open((*itFile).c_str(),"UPDATE") );
    if(not SampleFileList.back() or SampleFileList.back() == NULL){ badFiles.push_back(badFiles.size()-1); continue; }
    SampleTreeList.push_back( (TTree*) SampleFileList.back()->Get(TreeName.c_str()));
    if(SampleTreeList.back() == 0 or SampleTreeList.back() == NULL) SampleTreeList.erase(SampleTreeList.end());
  }

  SampleFileList.clear();

  // Book MVA Reader  Object --> one for each pT bin                                                                                                                        
  std::vector<TMVAReadingClass*> TMVAReaderVector ;


  // loop on the different training performed
  for(size_t iTrain = 0; iTrain < InputVariableFileParam.size() ; iTrain++){

    // loop on the single weights file that can be different since the training can be done in exclusive PT or PU bins
    std::vector<std::string> variableList = InputVariableFileParam.at(iTrain).getParameter<std::vector<std::string> >("InputVariableList");
    std::string methodName = InputVariableFileParam.at(iTrain).getParameter<std::string>("MethodName");
    std::vector<double> JetPTRegion = InputVariableFileParam.at(iTrain).getParameter<std::vector<double> >("JetPTRegion");
    std::vector<double> JetPURegion = InputVariableFileParam.at(iTrain).getParameter<std::vector<double> >("JetPURegion");
    std::vector<std::string> inputWeightFiles = InputVariableFileParam.at(iTrain).getParameter<std::vector<std::string> >("inputWeightFiles");

    TString NameBranch = Form("%s",methodName.c_str());

    for(size_t iVar = 0; iVar < variableList.size(); iVar++){
      TString VarName  = Form("%s",variableList.at(iVar).c_str());
      TString position = Form("[%d]",JetToRead);
      if(VarName.Contains(position)) VarName.ReplaceAll(position.Data(),"");
      NameBranch = Form("%s_%s",NameBranch.Data(),VarName.Data());
    }
    std::cout<<" NameBranch "<<NameBranch<<std::endl;

    TMVAReaderVector.push_back(new TMVAReadingClass(SampleTreeList,TreeName,inputWeightFiles,Label));
    			  
    TMVAReaderVector.back()->AddTrainingVariables(variableList,InputSpectatorList);
      
    TMVAReaderVector.back()->AddPrepareReader(LeptonType,PreselectionCutType,JetPTRegion,JetPURegion);
                
    TMVAReaderVector.back()->BookMVAWeight(methodName,std::string(NameBranch)); 
    
    TMVAReaderVector.back()->FillMVAWeight(JetToRead,optionCut);

  }
            
  //Print Output Plots                                                                                                                                                        
  std::cout<<std::endl;
  std::cout<<" Save Output Root File  ..  "<<std::endl;
  std::cout<<std::endl;

  return 0 ;
    
}
