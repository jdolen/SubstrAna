#include "../include/TMVAReadingClass.h"

// constructor from files
TMVAReadingClass::TMVAReadingClass(const std::vector<TFile*> & SampleFileList, const std::string & TreeName,
                                         const std::string & inputFilePath , const std::string & Label ){

  SetInputTree (SampleFileList, TreeName) ;

  SetReaderTree();

  SetLabel (Label);

  SetInputFilePath (inputFilePath) ;
  
  reader_ = new TMVA::Reader( TreeName_+"_"+Label_);

}

// constructor from trees
TMVAReadingClass::TMVAReadingClass(const std::vector<TTree*> & SampleTreeList, const std::string & TreeName,
		                         const std::string & inputFilePath , const std::string & Label ){


  SetTreeName (TreeName) ;

  SetInputTree (SampleTreeList) ;

  SetReaderTree (SampleTreeList) ;

  SetLabel (Label);

  SetInputFilePath (inputFilePath) ;

  reader_ = new TMVA::Reader( TreeName_+"_"+Label_);

}

// Deconstructor                                                                                                                                                         
TMVAReadingClass::~TMVAReadingClass(){

  for(size_t iTree = 0; iTree < SampleTreeList_.size() ; iTree++) { 
   if(SampleTreeList_.at(iTree)!=0)  delete SampleTreeList_.at(iTree) ; 
  }
  if(reader_!=0) reader_->Delete() ;

  if(setTrainingVariables_ !=0) delete setTrainingVariables_;
  if(setSpectatorVariables_ !=0) delete setSpectatorVariables_;

  for(size_t iTree = 0; iTree < SampleTreeList_.size() ; iTree++) { if(treeReader_.at(iTree)!=0) delete treeReader_.at(iTree) ;}

}

// Set Input variables and spectator 
void TMVAReadingClass::AddTrainingVariables ( const std::vector<std::string> & mapTrainingVariables, const std::vector<std::string> & mapSpectatorVariables){

  SetTrainingVariables(mapTrainingVariables);
  SetSpectatorVariables(mapSpectatorVariables);

  setTrainingVariables_  = new std::vector<Float_t> ( mapTrainingVariables_.size()) ; // create a vector of float values in order to read the input variables for the reader
  setSpectatorVariables_ = new std::vector<Float_t> ( mapSpectatorVariables_.size()) ;

  for( size_t iVar = 0 ; iVar < mapTrainingVariables_.size() ; iVar ++ )
    reader_->AddVariable(mapTrainingVariables_.at(iVar),&setTrainingVariables_->at(iVar)); // add variable to the reader with a right name which should match th input tree
  
  for( size_t iVar = 0 ; iVar < mapSpectatorVariables_.size() ; iVar ++ )
    reader_->AddSpectator(mapSpectatorVariables_.at(iVar),&setSpectatorVariables_->at(iVar));
  
  return;
}

// Prepare the tree -> set the preselection cut
void TMVAReadingClass::AddPrepareReader (const std::string & LeptonType, const std::string & preselectionCutType,
  					    std::vector<double> * JetPtBinOfTraining, const int & pTBin, 
                                            std::vector<double> * PileUpBinOfTraining, const int & puBin){

  if(JetPtBinOfTraining!=NULL || !JetPtBinOfTraining ){
    pTJetMin_ = JetPtBinOfTraining->at(pTBin) ;
    pTJetMax_ = JetPtBinOfTraining->at(pTBin+1) ;
  }
  else{
    pTJetMin_ = 0;
    pTJetMax_ = 2000 ;
  }

  if(PileUpBinOfTraining!=NULL || !PileUpBinOfTraining ){
    puMin_ = PileUpBinOfTraining->at(puBin) ;
    puMax_ = PileUpBinOfTraining->at(puBin+1) ;
  }
  else{
    puMin_ = 0;
    puMax_ = 100 ;
  }

  // store the cut type .. should be a string
  preselectionCutType_ = preselectionCutType ;
  LeptonType_ = LeptonType ;

  return ;
}


// Book the MVA method
void TMVAReadingClass::BookMVAWeight (const std::string & methodName, const std::string & weightFile, const std::string & nameBranch) {

  methodName_ = methodName ;
  weightFile_ = weightFile ;  // destination of weight file
  reader_->BookMVA(methodName_.c_str(),(inputFilePath_+"/"+weightFile_).c_str()) ; // book the MVA reader method
  nameBranch_ = nameBranch ;  // name of the new branch of the MVA output

  return ;

}

// Fill the MVA weight
void TMVAReadingClass::FillMVAWeight (const int & iJetToRead, const bool & optionCut) {

  // loop on the sample list and create the branch
  for(size_t iTree = 0; iTree < SampleTreeList_.size() ; iTree++){

      newBranch_  = SampleTreeList_.at(iTree)->Branch(nameBranch_.c_str(),&weight_,(nameBranch_+"/F").c_str()); // create the new branch as a float

      for( int iEntry = 0; iEntry < SampleTreeList_.at(iTree)->GetEntries(); iEntry++){ // Loop on tree entries

       if(iEntry%10000 == 0 ) std::cout<<" Entry = "<<iEntry<<std::endl;      
       SampleTreeList_.at(iTree)->GetEntry(iEntry);  
       weight_ = 0 ;

       // fill the value of the training and spectator variable for each event
       for( size_t iVar = 0 ; iVar < mapTrainingVariables_.size() ; iVar ++ ) 
         setTrainingVariables_->at(iVar) = treeReader_.at(iTree)->GetFloat(mapTrainingVariables_.at(iVar).c_str())->at(iJetToRead) ; // used treeReader class to fill the value of each input variable  in the vector

      for( size_t iVar = 0 ; iVar < mapSpectatorVariables_.size() ; iVar ++ ) 
	setSpectatorVariables_->at(iVar) = treeReader_.at(iTree)->GetFloat(mapSpectatorVariables_.at(iVar).c_str())->at(iJetToRead) ;

      // preselection type definition --> only the basic cut SB+SR since we want to skip only the events that are never used in the analysis
      bool isGoodEvent = false ;
      
      if(optionCut == true){
	
	if( preselectionCutType_ == "basicJetsCutCSA14" && (LeptonType_ == "Jets" || LeptonType_ == "jets") and  TreeName_ !="gen") {
	  isGoodEvent = treeReader_.at(iTree)->GetFloat("pt")->at(0) > 200 && treeReader_.at(iTree)->GetFloat("imatch")->at(0) >= 0 && (treeReader_.at(iTree)->GetFloat("pt")->at(0) >= pTJetMin_ && treeReader_.at(iTree)->GetFloat("pt")->at(0) <= pTJetMax_) && (treeReader_.at(iTree)->GetFloat("npu")->at(0) >= puMin_ && treeReader_.at(iTree)->GetFloat("npu")->at(0) <= puMax_) && fabs(treeReader_.at(iTree)->GetFloat("eta")->at(0))<2.4; 
	}	
	else if( preselectionCutType_ == "basicJetsCutCSA14" && (LeptonType_ == "Jets" || LeptonType_ == "jets") and  TreeName_ =="gen") {
           isGoodEvent = treeReader_.at(iTree)->GetFloat("pt")->at(0) > 200 && treeReader_.at(iTree)->GetFloat("imatch")->at(0) >= 0 && (treeReader_.at(iTree)->GetFloat("pt")->at(0) >= pTJetMin_ && treeReader_.at(iTree)->GetFloat("pt")->at(0) <= pTJetMax_) && (treeReader_.at(iTree)->GetFloat("npu")->at(0) >= puMin_ && treeReader_.at(iTree)->GetFloat("npu")->at(0) <= puMax_) ;
	}
	else if( preselectionCutType_ == "basicJetsCutCSA14TTbar" && (LeptonType_ == "Jets" || LeptonType_ == "jets") and  TreeName_ !="gen") {
           isGoodEvent = treeReader_.at(iTree)->GetFloat("pt")->at(0) > 200 && treeReader_.at(iTree)->GetFloat("imatch")->at(0) >= 0 && (treeReader_.at(iTree)->GetFloat("pt")->at(0) >= pTJetMin_ && treeReader_.at(iTree)->GetFloat("pt")->at(0) <= pTJetMax_) && (treeReader_.at(iTree)->GetFloat("npu")->at(0) >= puMin_ && treeReader_.at(iTree)->GetFloat("npu")->at(0) <= puMax_) && fabs(treeReader_.at(iTree)->GetFloat("eta")->at(0))<2.4 && treeReader_.at(iTree)->GetFloat("mprunedsafe_zcut_010_R_cut_050")->at(0) <=120;
	}
	else if( preselectionCutType_ == "basicJetsCutCSA14TTbar" && (LeptonType_ == "Jets" || LeptonType_ == "jets") and  TreeName_ =="gen") {
           isGoodEvent = treeReader_.at(iTree)->GetFloat("pt")->at(0) > 200 && treeReader_.at(iTree)->GetFloat("imatch")->at(0) >= 0 && (treeReader_.at(iTree)->GetFloat("pt")->at(0) >= pTJetMin_ && treeReader_.at(iTree)->GetFloat("pt")->at(0) <= pTJetMax_) && (treeReader_.at(iTree)->GetFloat("npu")->at(0) >= puMin_ && treeReader_.at(iTree)->GetFloat("npu")->at(0) <= puMax_) && treeReader_.at(iTree)->GetFloat("mprunedsafe_zcut_010_R_cut_050")->at(0) <=120;
	}
        else {  
           isGoodEvent = treeReader_.at(iTree)->GetFloat("pt")->at(0) > 200 && treeReader_.at(iTree)->GetFloat("imatch")->at(0) >= 0 && (treeReader_.at(iTree)->GetFloat("pt")->at(0) >= pTJetMin_ && treeReader_.at(iTree)->GetFloat("pt")->at(0) <= pTJetMax_) && (treeReader_.at(iTree)->GetFloat("npu")->at(0) >= puMin_ && treeReader_.at(iTree)->GetFloat("npu")->at(0) <= puMax_) ;
	   }	     	
     }
     else isGoodEvent = true ;      
     // if is a good event -> fill the branch with the value of the output      
      if(isGoodEvent) { weight_ = reader_->EvaluateMVA(methodName_.c_str());  // read the weight if the event is good -> re-applying the cut here
 	                newBranch_->Fill() ; // fill the branch
     }
     else { weight_ = -100. ; 
            newBranch_->Fill() ; // fill a default value
     } 
    }       
    SampleTreeList_.at(iTree)->Write("", TObject::kOverwrite); // rewrite the tree
  }
  return ;

}


// Set input tree
void TMVAReadingClass::SetInputTree (const std::vector<TFile*> & SampleFileList,  const std::string & TreeName){

  if (TreeName!="") TreeName_ = TreeName ;
  else TreeName_ = "chs" ;

  for(size_t iFile = 0 ; iFile < SampleFileList.size() ; iFile ++){
    if(SampleFileList.at(iFile)!=0)  SampleTreeList_.push_back((TTree*) SampleFileList.at(iFile)->Get(TreeName_.c_str()));
  }

  return ;

}

void TMVAReadingClass::SetInputTree (const std::vector<TTree*> & SampleTreeList){

  for(unsigned int iTree = 0; iTree< SampleTreeList.size(); iTree++){
    if(SampleTreeList.at(iTree)->GetEntries()>0) SampleTreeList_.push_back(SampleTreeList.at(iTree)) ;
  }

  return ;

}

// Set Training variables name
void TMVAReadingClass::SetTrainingVariables  (const std::vector<std::string> & mapTrainingVariables){

  if(mapTrainingVariables.size()!=0) mapTrainingVariables_ = mapTrainingVariables;

  return ;

}

void TMVAReadingClass::SetSpectatorVariables (const std::vector<std::string> & mapSpectatorVariables){

  if(mapSpectatorVariables.size()!=0) mapSpectatorVariables_ = mapSpectatorVariables;
 
  return ;

}

// Set input path
void TMVAReadingClass::SetInputFilePath ( const std::string & InputFilePath){

  if(!InputFilePath.empty()) inputFilePath_ = InputFilePath ;

  return;
}

void TMVAReadingClass::SetTreeName ( const std::string & TreeName ){

  if(TreeName !="") TreeName_ = TreeName ;
  else TreeName_ = "chs" ;

  return ;
}

void TMVAReadingClass::SetLabel ( const std::string & Label ){
  Label_ = Label ;
  return ;
}

void TMVAReadingClass::SetReaderTree(){

  for(size_t iTree = 0; iTree<SampleTreeList_.size() ; iTree++)
    treeReader_.push_back(new treeReader((TTree*) SampleTreeList_.at(iTree), false));
  return ;
}

void TMVAReadingClass::SetReaderTree(const std::vector<TTree*> & SampleTreeList){

  for(size_t iTree = 0; iTree<SampleTreeList.size() ; iTree++)
    treeReader_.push_back(new treeReader((TTree*) SampleTreeList.at(iTree), false));
  return ;
}
