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

  std::vector<edm::ParameterSet> InputWeightFileParam;
  if(Options.existsAs<std::vector<edm::ParameterSet> >("InputWeightFileParam"))
    InputWeightFileParam  = Options.getParameter<std::vector<edm::ParameterSet> >("InputWeightFileParam");
  else{ std::cout<<" No input weight file parameters are provided --> exit from the program "<<std::endl; return -1; }

  std::cout << " >>>>> Option:InputWeightFile size = " << InputWeightFileParam.size() << std::endl;

  for (unsigned int iCat = 0; iCat < InputWeightFileParam.size(); iCat++){
    std::cout << " " <<InputWeightFileParam.at(iCat).getParameter<std::string>("inputWeightFile")<<"  "<<InputWeightFileParam.at(iCat).getParameter<std::string>("useMethodName")<<std::endl;
  }    
  std::cout << std::endl;


  std::string PreselectionCutType ;
  if(Options.existsAs<std::string>("PreselectionCutType"))
    PreselectionCutType  = Options.getParameter<std::string>("PreselectionCutType");
  else{ std::cout<<" PreselectionCutType --> not found --> empty string "<<std::endl; PreselectionCutType = "none"; }

  std::cout<<" Option Preselection Cut = "<<PreselectionCutType<<std::endl;
  std::cout<<std::endl;

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
  std::string weightFile ;
  
  // loop on the pTbin training 
  for(size_t pTBin = 0; pTBin+1 < JetPtBinOfTraining.size() ; pTBin++){
    std::cout<<" pT bin of Training: Min = "<<JetPtBinOfTraining.at(pTBin)<<" Max = "<<JetPtBinOfTraining.at(pTBin+1)<<std::endl;
    std::cout<<std::endl;
    // loop on the pU training bins 
    for(size_t puBin = 0; puBin+1 < PileUpBinOfTraining.size() ; puBin++){
     std::cout<<" pileup bin of Training: Min = "<<PileUpBinOfTraining.at(puBin)<<" Max = "<<PileUpBinOfTraining.at(puBin+1)<<std::endl;
     std::cout<<std::endl;
     // Loop on the input weight file ... the name should contain some %s for the pTbin and pileUp bin value to be sobstituted
     for( size_t iWeightFile = 0 ; iWeightFile < InputWeightFileParam.size() ; iWeightFile++){ 

      TString fileName = Form("%s",(InputWeightFileParam.at(iWeightFile).getParameter<std::string>("inputWeightFile")).c_str());
      if(fileName.Contains("PTBin_%d_%d")) 
	 fileName.ReplaceAll("PTBin_%d_%d",Form("PTBin_%d_%d",int(JetPtBinOfTraining.at(pTBin)),int(JetPtBinOfTraining.at(pTBin+1))));
      if(fileName.Contains("PU_%d_%d")) 
	 fileName.ReplaceAll("PU_%d_%d",Form("PU_%d_%d",int(PileUpBinOfTraining.at(puBin)),int(PileUpBinOfTraining.at(puBin+1))));
      if(fileName.Contains("%d_%d"))
	 fileName.ReplaceAll("%d_%d",Form("%d_%d",int(PileUpBinOfTraining.at(pTBin)),int(PileUpBinOfTraining.at(pTBin+1))));
  
      std::cout<<" Weight File Name " <<fileName<<std::endl;
      std::cout<<std::endl;
     
      TMVAReaderVector.push_back(new TMVAReadingClass(SampleTreeList,TreeName,std::string(fileName),Label));
            
      TMVAReaderVector.back()->AddTrainingVariables(InputWeightFileParam.at(iWeightFile).getParameter<std::vector<std::string> >("inputVariableList"),InputSpectatorList);
      
      TMVAReaderVector.back()->AddPrepareReader(LeptonType,PreselectionCutType,&JetPtBinOfTraining,pTBin,&PileUpBinOfTraining,puBin);

      TString NameBranch = Form("%s_PTBin_%d_%d_PUBin_%d_%d",InputWeightFileParam.at(iWeightFile).getParameter<std::string>("useMethodName").c_str(),int(JetPtBinOfTraining.at(pTBin)),int(JetPtBinOfTraining.at(pTBin+1)),int(PileUpBinOfTraining.at(puBin)),int(PileUpBinOfTraining.at(puBin+1))) ;

      TMVAReaderVector.back()->BookMVAWeight(InputWeightFileParam.at(iWeightFile).getParameter<std::string>("useMethodName"),std::string(NameBranch)); 

      TMVAReaderVector.back()->FillMVAWeight(JetToRead,optionCut);
    }
   }
 }
     
  //Print Output Plots                                                                                                                                                        
  std::cout<<std::endl;
  std::cout<<" Save Output Root File  ..  "<<std::endl;
  std::cout<<std::endl;

  return 0 ;
    
}
