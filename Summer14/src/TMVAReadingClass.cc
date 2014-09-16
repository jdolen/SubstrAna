#include "../include/TMVAReadingClass.h"

// constructor from files
TMVAReadingClass::TMVAReadingClass(const std::vector<TFile*> & SampleFileList, const std::string & TreeName,
				   const std::vector<std::string> & inputFileWeight , const std::string & Label ){

  SetInputTree (SampleFileList, TreeName) ;
  SetReaderTree();
  SetLabel (Label);
  SetInputFileWeight (inputFileWeight) ;

  for(size_t iWeight = 0 ; iWeight < inputFileWeight.size(); iWeight++){
    TString NameReader = Form("%s_%s_weight_%d",TreeName_.c_str(),Label_.c_str(),int(iWeight));
    reader_.push_back( new TMVA::Reader(NameReader.Data()));
  }
}

// constructor from trees
TMVAReadingClass::TMVAReadingClass(const std::vector<TTree*> & SampleTreeList, const std::string & TreeName,
				   const std::vector<std::string> & inputFileWeight , const std::string & Label ){

  SetTreeName (TreeName) ;
  SetInputTree (SampleTreeList) ;
  SetReaderTree (SampleTreeList) ;
  SetLabel (Label);
  SetInputFileWeight (inputFileWeight) ;

  for(size_t iWeight = 0 ; iWeight < inputFileWeight.size(); iWeight++){
    TString NameReader = Form("%s_%s_weight_%d",TreeName_.c_str(),Label_.c_str(),int(iWeight));
    reader_.push_back( new TMVA::Reader(NameReader.Data()));
  }
}

// Deconstructor                                                                                                                                                         
TMVAReadingClass::~TMVAReadingClass(){

  for(size_t iTree = 0; iTree < SampleTreeList_.size() ; iTree++) { 
   if(SampleTreeList_.at(iTree)!=0)  delete SampleTreeList_.at(iTree) ; 
  }

  for(size_t iReader = 0; iReader < reader_.size(); iReader++){
    if(reader_.at(iReader)!=0) reader_.at(iReader)->Delete() ;
  }

  if(setTrainingVariables_ !=0)  delete setTrainingVariables_;
  if(setSpectatorVariables_ !=0) delete setSpectatorVariables_;

  for(size_t iTree = 0; iTree < SampleTreeList_.size() ; iTree++) { 
     if(treeReader_.at(iTree)!=0) delete treeReader_.at(iTree) ;
  }

}

// Set Input variables and spectator 
void TMVAReadingClass::AddTrainingVariables ( const std::vector<std::string> & mapTrainingVariables, const std::vector<std::string> & mapSpectatorVariables){

  SetTrainingVariables(mapTrainingVariables);
  SetSpectatorVariables(mapSpectatorVariables);

  setTrainingVariables_  = new std::vector<Float_t> ( mapTrainingVariables_.size()) ; // create a vector of float values in order to read the input variables for the reader
  setSpectatorVariables_ = new std::vector<Float_t> ( mapSpectatorVariables_.size()) ;

  for(size_t iWeight = 0; iWeight < reader_.size(); iWeight++){

   for(size_t iVar = 0 ; iVar < mapTrainingVariables_.size() ; iVar ++ )
     reader_.at(iWeight)->AddVariable(mapTrainingVariables_.at(iVar),&setTrainingVariables_->at(iVar)); // add variable to the reader with a right name which should match th input tree
  
   for(size_t iVar = 0 ; iVar < mapSpectatorVariables_.size() ; iVar ++ )
     reader_.at(iWeight)->AddSpectator(mapSpectatorVariables_.at(iVar),&setSpectatorVariables_->at(iVar));
  }

  return;
}

// Prepare the tree -> set the preselection cut
void TMVAReadingClass::AddPrepareReader (const std::string & LeptonType, const std::string & preselectionCutType,
  					 const std::vector<double> & JetPtBinOfTraining, const std::vector<double> & PileUpBinOfTraining){

  if(JetPtBinOfTraining.size()!=0)  JetPtBinOfTraining_ = JetPtBinOfTraining ;
  if(PileUpBinOfTraining.size()!=0) PileUpBinOfTraining_ = PileUpBinOfTraining;
  // store the cut type .. should be a string
  preselectionCutType_ = preselectionCutType ;
  LeptonType_ = LeptonType ;

  return ;
}


// Book the MVA method
void TMVAReadingClass::BookMVAWeight (const std::string & methodName, const std::string & nameBranch) {

  methodName_ = methodName ;
  nameBranch_ = nameBranch ;  // name of the new branch of the MVA output

  for(size_t iWeight = 0; iWeight < inputFileWeight_.size(); iWeight++) 
    reader_.at(iWeight)->BookMVA(methodName_.c_str(),inputFileWeight_.at(iWeight).c_str()) ; // book the MVA reader method

  return ;
}

// Fill the MVA weight
void TMVAReadingClass::FillMVAWeight (const int & iJetToRead, const bool & optionCut) {

       
  // loop on the sample list and create the branch
  for(size_t iTree = 0; iTree < SampleTreeList_.size() ; iTree++){
    
    treeFormulaVar_.clear(); // eliminate the formulas already there
    treeFormulaSpec_.clear(); // eliminate the formulas already there

    // Evaluate the formulas
    for( size_t iVar = 0 ; iVar < mapTrainingVariables_.size() ; iVar ++ ) {
	 TString variableName = Form("%s",mapTrainingVariables_.at(iVar).c_str());
         TString position     = Form("[%d]",iJetToRead);
         if(variableName.Contains(position)) variableName.ReplaceAll(position.Data(),"");
         treeFormulaVar_.push_back(new TTreeFormula (variableName.Data(),variableName.Data(),SampleTreeList_.at(iTree))); 
    }

    for( size_t iVar = 0 ; iVar < mapSpectatorVariables_.size() ; iVar ++ ){
	 TString variableName = Form("%s",mapSpectatorVariables_.at(iVar).c_str());
         TString position     = Form("[%d]",iJetToRead);
         if(variableName.Contains(position)) variableName.ReplaceAll(position.Data(),"");
	 treeFormulaSpec_.push_back(new TTreeFormula (variableName.Data(),variableName.Data(),SampleTreeList_.at(iTree))); 
    }        

    if(SampleTreeList_.at(iTree)->FindBranch(nameBranch_.c_str())) // if already exist the branch take it and fill again (when trainings are done in exclusive pT bins you might want just one branch filled alternatively as a function of the events      
      newBranch_ = SampleTreeList_.at(iTree)->GetBranch(nameBranch_.c_str());
    else // create a new branch from scratch
      newBranch_  = SampleTreeList_.at(iTree)->Branch(nameBranch_.c_str(),&weight_,(nameBranch_+"/F").c_str()); // create the new branch as a float

      for( int iEntry = 0; iEntry < SampleTreeList_.at(iTree)->GetEntries(); iEntry++){ // Loop on tree entries

       if(iEntry%100 == 0 ) std::cout<<" Entry = "<<iEntry<<std::endl;      
       SampleTreeList_.at(iTree)->GetEntry(iEntry);  // load entry information -> treeReader already set all the status to 1
       weight_ = 0 ; // value to be put in the branch
       
       // select only events in the correct pt and pu bin 
       double pTjet = treeReader_.at(iTree)->GetFloat("pt")->at(iJetToRead) ;
       double pUjet = treeReader_.at(iTree)->getInt("npu")[iJetToRead] ;
 
       int pTPosition = -1;
       for(size_t iPT = 0 ; iPT < JetPtBinOfTraining_.size(); iPT = iPT + 2){
	   if(pTjet > JetPtBinOfTraining_.at(iPT) and pTjet <= JetPtBinOfTraining_.at(iPT+1)) {
	     pTPosition = iPT ; 
             continue;
	   }
       }

       int pUPosition = -1;
       for(size_t iPU = 0 ; iPU < PileUpBinOfTraining_.size(); iPU = iPU +2){
	   if(pUjet > PileUpBinOfTraining_.at(iPU) and pUjet <= PileUpBinOfTraining_.at(iPU+1)) {
	     pUPosition = iPU ; 
             continue;
	   }
       }
              
       bool isGoodEvent = false ;      
       if(optionCut == true){	
  	  if( preselectionCutType_ == "basicJetsCutCSA14" && (LeptonType_ == "Jets" || LeptonType_ == "jets") and  TreeName_ !="gen") 
	    isGoodEvent = treeReader_.at(iTree)->GetFloat("pt")->at(iJetToRead) > 200 && treeReader_.at(iTree)->GetInt("imatch")->at(iJetToRead) >= 0 && 
                          fabs(treeReader_.at(iTree)->GetFloat("eta")->at(iJetToRead))<2.4; 	  

	  else if( preselectionCutType_ == "basicJetsCutCSA14" && (LeptonType_ == "Jets" || LeptonType_ == "jets") and  TreeName_ =="gen") 
            isGoodEvent = treeReader_.at(iTree)->GetFloat("pt")->at(iJetToRead) > 200 && treeReader_.at(iTree)->GetInt("imatch")->at(iJetToRead) >= 0  ;	  

	  else if( preselectionCutType_ == "basicJetsCutCSA14NoGen" && (LeptonType_ == "Jets" || LeptonType_ == "jets") and  TreeName_ !="gen") 
            isGoodEvent = treeReader_.at(iTree)->GetFloat("pt")->at(iJetToRead) > 200 && fabs(treeReader_.at(iTree)->GetFloat("eta")->at(iJetToRead))<2.4;	  

  	  else if( preselectionCutType_ == "basicJetsCutCSA14Mass" && (LeptonType_ == "Jets" || LeptonType_ == "jets") and  TreeName_ !="gen") 
	    isGoodEvent = treeReader_.at(iTree)->GetFloat("pt")->at(iJetToRead) > 200 && treeReader_.at(iTree)->GetInt("imatch")->at(iJetToRead) >= 0 && 
                          fabs(treeReader_.at(iTree)->GetFloat("eta")->at(iJetToRead))<2.4 && treeReader_.at(iTree)->GetFloat("mprunedsafe_z_010_R_050")->at(iJetToRead) > 65 && 
                          treeReader_.at(iTree)->GetFloat("mprunedsafe_z_010_R_050")->at(iJetToRead) < 105; 
	  
	  else if( preselectionCutType_ == "basicJetsCutCSA14Mass" && (LeptonType_ == "Jets" || LeptonType_ == "jets") and  TreeName_ =="gen") 
	    isGoodEvent = treeReader_.at(iTree)->GetFloat("pt")->at(iJetToRead) > 200 && treeReader_.at(iTree)->GetInt("imatch")->at(iJetToRead) >= 0 && 
                          treeReader_.at(iTree)->GetFloat("mprunedsafe_z_010_R_050")->at(iJetToRead) > 65 && 
                          treeReader_.at(iTree)->GetFloat("mprunedsafe_z_010_R_050")->at(iJetToRead) < 105  ;
	  
   	  else if( preselectionCutType_ == "basicJetsCutCSA14TTbar" && (LeptonType_ == "Jets" || LeptonType_ == "jets") and  TreeName_ !="gen") 
            isGoodEvent = treeReader_.at(iTree)->GetFloat("pt")->at(iJetToRead) > 200 && treeReader_.at(iTree)->GetInt("imatch")->at(iJetToRead) >= 0 && 
                          fabs(treeReader_.at(iTree)->GetFloat("eta")->at(iJetToRead))<2.4 && treeReader_.at(iTree)->GetFloat("mprunedsafe_z_010_R_050")->at(iJetToRead) <=120;
	  
  	  else if( preselectionCutType_ == "basicJetsCutCSA14TTbar" && (LeptonType_ == "Jets" || LeptonType_ == "jets") and  TreeName_ =="gen")
            isGoodEvent = treeReader_.at(iTree)->GetFloat("pt")->at(iJetToRead) > 200 && treeReader_.at(iTree)->GetInt("imatch")->at(iJetToRead) >= 0 && 
                          treeReader_.at(iTree)->GetFloat("mprunedsafe_z_010_R_050")->at(iJetToRead) <=120;
	  
          else isGoodEvent = treeReader_.at(iTree)->GetFloat("pt")->at(iJetToRead) > 200 && treeReader_.at(iTree)->GetInt("imatch")->at(iJetToRead) >= 0 ;
       }	  	     	       
       else isGoodEvent = true ;      

       //std::cout<<" pt "<<treeReader_.at(iTree)->GetFloat("pt")->at(iJetToRead)<<" puJet "<<pUjet<<" match "<<treeReader_.at(iTree)->GetInt("imatch")->at(iJetToRead)<<" eta "<<fabs(treeReader_.at(iTree)->GetFloat("eta")->at(iJetToRead))<<" isGoodEvent "<<isGoodEvent<<" ptpos "<<pTPosition<<" pile-up "<<pUPosition<<" size "<<JetPtBinOfTraining_.size()<<"  "<<PileUpBinOfTraining_.size()<<std::endl;

       // read input variables
       for( size_t iVar = 0 ; iVar < mapTrainingVariables_.size() ; iVar ++ ) {
         setTrainingVariables_->at(iVar) = treeFormulaVar_.at(iVar)->EvalInstance();
       }
       for( size_t iVar = 0 ; iVar < mapSpectatorVariables_.size() ; iVar ++ ) {
         setSpectatorVariables_->at(iVar) = treeFormulaSpec_.at(iVar)->EvalInstance();
       }
       // take the correct weight file to be read for this event --> look for a string  
       int weightFilePosition = -1 ;
       if( pTPosition !=-1 and pUPosition !=-1){
	 TString weightFileSuffix = Form("PTBin_%d_%d_PU_%d_%d",int(JetPtBinOfTraining_.at(pTPosition)),int(JetPtBinOfTraining_.at(pTPosition+1)),int(PileUpBinOfTraining_.at(pUPosition)),int(PileUpBinOfTraining_.at(pUPosition+1)));
        for(size_t iWeight = 0; iWeight < inputFileWeight_.size(); iWeight++){
         if(TString(inputFileWeight_.at(iWeight).c_str()).Contains(weightFileSuffix.Data())) weightFilePosition = iWeight;
	}
       }
       
       //if(weightFilePosition!=-1) std::cout<<" weight file "<<weightFilePosition<<" name "<<inputFileWeight_.at(weightFilePosition)<<std::endl;
											                     
     // if is a good event -> fill the branch with the value of the output      
      if(isGoodEvent and pTPosition !=-1 and pUPosition !=-1 and weightFilePosition !=-1) { 
  	  weight_ = reader_.at(weightFilePosition)->EvaluateMVA(methodName_.c_str());  // read the weight if the event is good -> re-applying the cut here
	  newBranch_->Fill() ; // fill the branch
     }
     else { weight_ = -100. ; 
            newBranch_->Fill() ; // fill a default value
     }        
    //std::cout<<" weight "<<weight_<<std::endl;  
      
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
void TMVAReadingClass::SetInputFileWeight ( const std::vector<std::string> & InputFileWeight){

  if(!InputFileWeight.empty()) inputFileWeight_ = InputFileWeight ;

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
    treeReader_.push_back(new treeReader((TTree*) SampleTreeList.at(iTree),false));
  return ;
}
