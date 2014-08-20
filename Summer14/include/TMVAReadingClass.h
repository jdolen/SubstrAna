#ifndef TMVAReadingClass_h
#define TMVAReadingClass_h

#include <iostream>
#include <cstdlib>
#include <string>
#include <map>
#include <vector>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TCut.h"

#include "treeReader.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/MsgLogger.h"
#include "TMVA/Config.h"
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#endif

class TMVAReadingClass {

 public :

  // default constructor
  TMVAReadingClass(){};
  
  // constructor giving input files
  TMVAReadingClass(const std::vector<TFile*> & SampleFileList, const std::string & TreeName,
                      const std::string & InputFilePath , const std::string & Label );

  // constructor giving input trees
  TMVAReadingClass(const std::vector<TTree*> & SampleTreeList, const std::string & TreeName,
                      const std::string & InputFilePath , const std::string & Label);

  // default de-constructor
  ~TMVAReadingClass();

  // Add training variables names 
  void AddTrainingVariables ( const std::vector<std::string> & mapTrainingVariables, const std::vector<std::string> & mapSpectatorVariables);
  
  // prepare the preselection cut and the jet pT of training
  void AddPrepareReader (const std::string & LeptonType, const std::string & preselectionCutType, std::vector<double> * JetPtBinOfTraining = NULL , const int & pTBin = 0, std::vector<double>* PileUpBinOfTraining = NULL, const int & puBin = 0);

  // Book MVA method used in the training + weight file for input TMVA weights and the branch og the name where to store the info
  void BookMVAWeight (const std::string & methodName, const std::string & weightFile, const std::string & nameBranch) ;

  // Fill the MVA weight value
  void FillMVAWeight (const int & iJetToRead = 0, const bool & optionCut = true) ;

  // Set input tree from file
  void SetInputTree (const std::vector<TFile*> & SampleFileList,  const std::string & TreeName = "WJet");

  // Set input tree from tree
  void SetInputTree (const std::vector<TTree*> & SampleTreeList);

  // Set input training variables
  void SetTrainingVariables  (const std::vector<std::string> & mapTrainingVariables);

  // Set input spectator variables
  void SetSpectatorVariables (const std::vector<std::string> & mapSpectatorVariables);

  // Set input file path
  void SetInputFilePath ( const std::string & InputFilePath);

  // Set tree name
  void SetTreeName ( const std::string & TreeName );

  // Set label
  void SetLabel ( const std::string & Label );

  // Set weight file name
  void SetWeightFile ( const std::string & weightFile ) ;

  // Set reader tree 
  void SetReaderTree ();
  void SetReaderTree (const std::vector<TTree*> & SampleTreeList);


 private :

  // pT range of the related training
  double pTJetMin_ ;
  double pTJetMax_ ;
  // PU range of the related training
  double puMin_ ;
  double puMax_ ;

  // sample tree list to fill the weight
  std::vector<TTree*> SampleTreeList_ ;

  // input set of training and spectator variables
  std::vector<std::string> mapTrainingVariables_ ;
  std::vector<std::string> mapSpectatorVariables_ ;

  std::vector<Float_t> *setTrainingVariables_ ;
  std::vector<Float_t> *setSpectatorVariables_ ;

  // weight value
  Float_t weight_ ;

  // tree Name and Label
  std::string TreeName_ ;
  std::string Label_ ;

  // input file path
  std::string inputFilePath_ ;

  // method Name, weight file and name of the output branch
  std::string methodName_ ;
  std::string weightFile_ ;
  std::string nameBranch_ ;
 
  // output pointer to the branch
  TBranch* newBranch_ ;

  // preselection cut string and lepton type
  std::string preselectionCutType_ ;
  std::string LeptonType_ ;

  // Reader Object
  TMVA::Reader*  reader_ ;

  // TreeReader in order to read the inputTree branches
  std::vector<treeReader*> treeReader_ ;

};

#endif
