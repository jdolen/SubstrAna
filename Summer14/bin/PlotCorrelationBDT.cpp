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

#include "TKey.h"
#include "TObjArray.h"
#include "TClass.h"
#include "TH2F.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

// To calculate covariance matrix depending on the input
void calculateCovarianceMatrix ( std::vector<TTree*> & inputTrees, TMatrixDSym & correlationMatrixS, TMatrixDSym & correlationMatrixB );
void calculateCovarianceMatrix ( std::vector<TTree*> & inputTrees, TMatrixDSym & correlationMatrix, std::vector<std::string> & variableName, std::string & cut);

// To plot covariance matrix
void plotCovarianceMatrix      ( TCanvas* cCorrelation, TH2D* Matrix, std::string outputDirectory, std::vector<std::string> labelNames);

// diagonalize and input matrix
TMatrixDSym eigenTransform (TMatrixDSym& Matrix);

// dimensional orderign for PCA 
TMatrixDSym dimOrder(TMatrixDSym & Matrix,TMatrixDSym & Eigen,TMatrixD & EigenInv);

/// Main programme 
int main (int argc, char **argv){
  if(argc<3){ std::cout<<" Not correct number of input parameter --> Need Just one cfg file exit and 0 if run on TMVA output, 1 if run on TTree "<<std::endl; std::exit(EXIT_FAILURE); }

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

  if(atoi(argv[3]) == 0){ // analyze TMVA output files .. one file == one trainning + a variables name  

   std::vector<edm::ParameterSet> InputInformationParamLowPU ;
   if(Options.existsAs<std::vector<edm::ParameterSet>>("InputLowPUFiles"))
    InputInformationParamLowPU = Options.getParameter<std::vector<edm::ParameterSet>>("InputLowPUFiles");
   else{ std::cout<<" Exit from code, no input set found for low pile-up files"<<std::endl; std::exit(EXIT_FAILURE) ; }

   std::vector<edm::ParameterSet> InputInformationParamHighPU ;
   if(Options.existsAs<std::vector<edm::ParameterSet>>("InputHighPUFiles"))
   InputInformationParamHighPU = Options.getParameter<std::vector<edm::ParameterSet>>("InputHighPUFiles");
   else{ std::cout<<" Exit from code, no input set found for high pile-up files"<<std::endl; std::exit(EXIT_FAILURE) ; }


   std::string outputDirectory;
   if(Options.existsAs<std::string>("outputDirectory"))
    outputDirectory = Options.getParameter<std::string>("outputDirectory");
   else{ outputDirectory = "output/"; }
 
   system(("mkdir -p "+outputDirectory).c_str());

   std::vector<TTree*> inputLowPileUpTrees ;  // TTree for input variables inside the TMVA output file
   std::vector<std::string> reducedNameLowPileUp ; // reduced name for the variables
   TFile* inputFile = NULL ;

   // take all the TH1 outputs for S and B for  each variable in input-> low Pile Up
   std::vector<edm::ParameterSet>::const_iterator itLowPileUp = InputInformationParamLowPU.begin();
   for( ; itLowPileUp != InputInformationParamLowPU.end() ; ++itLowPileUp){
    inputFile = TFile::Open((*itLowPileUp).getParameter<std::string>("fileName").c_str());
    if(inputFile == 0 or inputFile == NULL) continue;
    inputLowPileUpTrees.push_back((TTree*) inputFile->Get("TestTree")); // make the tree list
    reducedNameLowPileUp.push_back((*itLowPileUp).getParameter<std::string>("variableName"));
   }

   // Build Up the correlation matrix
   TMatrixDSym correlationMatrixS_lowPileUP(reducedNameLowPileUp.size()) ;
   TMatrixDSym correlationMatrixB_lowPileUP(reducedNameLowPileUp.size()) ;
   for( unsigned int irow = 0 ; irow < reducedNameLowPileUp.size() ; irow++){
    for( unsigned int icolum = 0 ; icolum < reducedNameLowPileUp.size() ; icolum++){
     correlationMatrixS_lowPileUP(irow,icolum) = 0;
     correlationMatrixB_lowPileUP(irow,icolum) = 0;
    }
   }

   calculateCovarianceMatrix(inputLowPileUpTrees,correlationMatrixS_lowPileUP,correlationMatrixB_lowPileUP); // calculate mean value for S and B

   TH2D* Matrix_S_lowPileUp = new TH2D(correlationMatrixS_lowPileUP); // make TH2D for plotting reasons
   Matrix_S_lowPileUp->SetName("Matrix_S_lowPileUp");
   TH2D* Matrix_B_lowPileUp = new TH2D(correlationMatrixB_lowPileUP);
   Matrix_B_lowPileUp->SetName("Matrix_B_lowPileUp");

   Matrix_S_lowPileUp->Scale(100); // correlation in %
   Matrix_B_lowPileUp->Scale(100);

   ////////////////////////
   TCanvas* cCorrelationSignal = new TCanvas("CorrelationBDT_S",Form("Correlation Matrix Signal"),180,52,550,550);
   plotCovarianceMatrix(cCorrelationSignal,Matrix_S_lowPileUp,outputDirectory,reducedNameLowPileUp);

   ////////////////////////
   TCanvas* cCorrelationBackground = new TCanvas("CorrelationBDT_B",Form("Correlation Matrix Background"),180,52,550,550);
   plotCovarianceMatrix(cCorrelationBackground,Matrix_B_lowPileUp,outputDirectory,reducedNameLowPileUp);

   //////////////////////////////////////////////////////////////////////
   //// high pile-up correlation
   //////////////////////////////////////////////////////////////////////

   inputLowPileUpTrees.clear() ;
   reducedNameLowPileUp.clear() ;
   // take all the TH1 outputs for S and B for  each variable in input-> low Pile Up
   itLowPileUp = InputInformationParamHighPU.begin();
   for( ; itLowPileUp != InputInformationParamHighPU.end() ; ++itLowPileUp){
    inputFile = TFile::Open((*itLowPileUp).getParameter<std::string>("fileName").c_str());
    if(inputFile == 0 or inputFile == NULL) continue;
    inputLowPileUpTrees.push_back((TTree*) inputFile->Get("TestTree"));
    reducedNameLowPileUp.push_back((*itLowPileUp).getParameter<std::string>("variableName"));
   }
  
   // Build Up the correlation matrix
   TMatrixDSym correlationMatrixS_highPileUP(reducedNameLowPileUp.size()) ;
   TMatrixDSym correlationMatrixB_highPileUP(reducedNameLowPileUp.size()) ;
   for( unsigned int irow = 0 ; irow < reducedNameLowPileUp.size() ; irow++){
    for( unsigned int icolum = 0 ; icolum < reducedNameLowPileUp.size() ; icolum++){
     correlationMatrixS_highPileUP(irow,icolum) = 0;
     correlationMatrixB_highPileUP(irow,icolum) = 0;
    }
   }

   calculateCovarianceMatrix(inputLowPileUpTrees,correlationMatrixS_highPileUP,correlationMatrixB_highPileUP); // calculate mean value for S and B
  
   TH2D* Matrix_S_highPileUp = new TH2D(correlationMatrixS_highPileUP);
   Matrix_S_highPileUp->SetName("Matrix_S_highPileUp");
   TH2D* Matrix_B_highPileUp = new TH2D(correlationMatrixB_highPileUP);
   Matrix_B_highPileUp->SetName("Matrix_B_highPileUp");

   Matrix_S_highPileUp->Scale(100); // correlation in %
   Matrix_B_highPileUp->Scale(100);

   for( int iBinX = 0 ; iBinX < Matrix_S_highPileUp->GetNbinsX() ; iBinX++){
    for( int iBinY = 0 ; iBinY < Matrix_S_highPileUp->GetNbinsX() ; iBinY++){
     Matrix_S_highPileUp->SetBinContent(iBinX+1,iBinY+1,int(Matrix_S_highPileUp->GetBinContent(iBinX+1,iBinY+1)));
     Matrix_B_highPileUp->SetBinContent(iBinX+1,iBinY+1,int(Matrix_B_highPileUp->GetBinContent(iBinX+1,iBinY+1)));
     if(iBinX+1 == iBinY+1 ) {
         Matrix_S_highPileUp->SetBinContent(iBinX+1,iBinY+1,100);
         Matrix_B_highPileUp->SetBinContent(iBinX+1,iBinY+1,100);
     }
    }
   }

   TCanvas* cCorrelationSignal_highPU = new TCanvas("CorrelationBDT_S_highPU",Form("Correlation Matrix Signal highPU"),180,52,550,550);
   plotCovarianceMatrix(cCorrelationSignal_highPU,Matrix_S_highPileUp,outputDirectory,reducedNameLowPileUp);
   TCanvas* cCorrelationBackground_highPU = new TCanvas("CorrelationBDT_B_highPU",Form("Correlation Matrix Background"),180,52,550,550);
   plotCovarianceMatrix(cCorrelationBackground_highPU,Matrix_B_highPileUp,outputDirectory,reducedNameLowPileUp);


   ////////////////////////////
   // Correlation difference among high pile-up bin and low pile-up bin for both S and B
   ////////////////////////////

   TH2D* CorrelationS_difference = new TH2D("CorrelationS_difference","",reducedNameLowPileUp.size(),0,reducedNameLowPileUp.size(),reducedNameLowPileUp.size(),0,reducedNameLowPileUp.size());
   TH2D* CorrelationB_difference = new TH2D("CorrelationB_difference","",reducedNameLowPileUp.size(),0,reducedNameLowPileUp.size(),reducedNameLowPileUp.size(),0,reducedNameLowPileUp.size());

   for(int iBinX = 0; iBinX < Matrix_S_highPileUp->GetNbinsX(); iBinX++){
    for(int iBinY = 0; iBinY < Matrix_S_highPileUp->GetNbinsY(); iBinY++){  
     CorrelationS_difference->SetBinContent(iBinX+1,iBinY+1,Matrix_S_highPileUp->GetBinContent(iBinX+1,iBinY+1)-Matrix_S_lowPileUp->GetBinContent(iBinX+1,iBinY+1));
     CorrelationB_difference->SetBinContent(iBinX+1,iBinY+1,Matrix_B_highPileUp->GetBinContent(iBinX+1,iBinY+1)-Matrix_B_lowPileUp->GetBinContent(iBinX+1,iBinY+1));
    }
   }

   TCanvas* cCorrelationSignal_difference = new TCanvas("cCorrelationSignal_difference",Form("Correlation Matrix Difference Signal"),180,52,550,550);
   plotCovarianceMatrix(cCorrelationSignal_difference,CorrelationS_difference,outputDirectory,reducedNameLowPileUp);
   TCanvas* cCorrelationBackground_difference = new TCanvas("cCorrelationBackground_difference",Form("Correlation Matrix Difference Background"),180,52,550,550);
   plotCovarianceMatrix(cCorrelationBackground_difference,CorrelationS_difference,outputDirectory,reducedNameLowPileUp);

  }

  else{
  
    std::vector<std::string> InputInformationParamLowPUSignal ; // here each variable one root file; S and B. in the input also a cut should be provided and the tree branch name of each variable
    if(Options.existsAs<std::vector<std::string>>("InputLowPUFilesSignal"))
     InputInformationParamLowPUSignal = Options.getParameter<std::vector<std::string>>("InputLowPUFilesSignal");
    else{ std::cout<<" Exit from code, no input set found for low pile-up files Signal"<<std::endl; std::exit(EXIT_FAILURE) ; }

    std::vector<std::string> InputInformationParamHighPUSignal ;
    if(Options.existsAs<std::vector<std::string>>("InputHighPUFilesSignal"))
     InputInformationParamHighPUSignal = Options.getParameter<std::vector<std::string>>("InputHighPUFiles");
    else{ std::cout<<" Exit from code, no input set found for high pile-up files Signal"<<std::endl; std::exit(EXIT_FAILURE) ; }

    std::vector<std::string> InputInformationParamLowPUBackground ;
    if(Options.existsAs<std::vector<std::string>>("InputLowPUFilesBackground"))
     InputInformationParamLowPUBackground = Options.getParameter<std::vector<std::string>>("InputLowPUFilesBackground");
    else{ std::cout<<" Exit from code, no input set found for low pile-up files Background"<<std::endl; std::exit(EXIT_FAILURE) ; }

    std::vector<std::string> InputInformationParamHighPUBackground ;
    if(Options.existsAs<std::vector<std::string>>("InputHighPUFilesBackground"))
     InputInformationParamHighPUBackground = Options.getParameter<std::vector<std::string>>("InputHighPUFilesBackground");
    else{ std::cout<<" Exit from code, no input set found for high pile-up files Background"<<std::endl; std::exit(EXIT_FAILURE) ; }

    std::vector<edm::ParameterSet> InputVariables ;
    if(Options.existsAs<std::vector<edm::ParameterSet>>("InputVariables"))
     InputVariables = Options.getParameter<std::vector<edm::ParameterSet>>("InputVariables");
    else{ std::cout<<" Exit from code, no input variables found"<<std::endl; std::exit(EXIT_FAILURE) ; }


    std::string outputDirectory;
    if(Options.existsAs<std::string>("outputDirectory"))
     outputDirectory = Options.getParameter<std::string>("outputDirectory");
    else{ outputDirectory = "output/"; }

    system(("mkdir -p "+outputDirectory).c_str());

    std::string treeName;
    if(Options.existsAs<std::string>("treeName"))
     treeName = Options.getParameter<std::string>("treeName");
    else{ treeName = ""; }

    std::string cutStringSignal;
    if(Options.existsAs<std::string>("cutStringSignal"))
     cutStringSignal = Options.getParameter<std::string>("cutStringSignal");
    else{ cutStringSignal = ""; }

    std::string cutStringBackground;
    if(Options.existsAs<std::string>("cutStringBackground"))
     cutStringBackground = Options.getParameter<std::string>("cutStringBackground");
    else{ cutStringBackground = ""; }

    std::vector<TTree*> inputLowPileUpTreesSignal ;
    std::vector<TTree*> inputLowPileUpTreesBackground ;
    std::vector<std::string> reducedName ;
    std::vector<std::string> variablesName ;

    TFile* inputFile = NULL ;
    // take all the signal trees low pile up
    std::vector<std::string>::const_iterator itLowPileUp = InputInformationParamLowPUSignal.begin();
    for( ; itLowPileUp != InputInformationParamLowPUSignal.end() ; ++itLowPileUp){
     inputFile = TFile::Open((*itLowPileUp).c_str());
     if(inputFile == 0 or inputFile == NULL) continue;
     inputLowPileUpTreesSignal.push_back((TTree*) inputFile->Get(treeName.c_str())); // tree name should be provided
    }
    // take all the background trees low pile up
    itLowPileUp = InputInformationParamLowPUBackground.begin();
    for( ; itLowPileUp != InputInformationParamLowPUBackground.end() ; ++itLowPileUp){
     inputFile = TFile::Open((*itLowPileUp).c_str());
     if(inputFile == 0 or inputFile == NULL) continue;
     inputLowPileUpTreesBackground.push_back((TTree*) inputFile->Get(treeName.c_str())); // tree name should be provided
    }

    //take the variables
    std::vector<edm::ParameterSet>::const_iterator itVar = InputVariables.begin();
    for( ; itVar != InputVariables.end() ; ++itVar){
      variablesName.push_back((*itVar).getParameter<std::string>("variableName"));
      variablesName.push_back((*itVar).getParameter<std::string>("reducedName"));
    }

    // Build Up the correlation matrix
    TMatrixDSym correlationMatrixS_lowPileUP(reducedName.size()) ;
    TMatrixDSym correlationMatrixB_lowPileUP(reducedName.size()) ;
    for( unsigned int irow = 0 ; irow < reducedName.size() ; irow++){
     for( unsigned int icolum = 0 ; icolum < reducedName.size() ; icolum++){
      correlationMatrixS_lowPileUP(irow,icolum) = 0;
      correlationMatrixB_lowPileUP(irow,icolum) = 0;
     }
    }

    calculateCovarianceMatrix(inputLowPileUpTreesSignal,correlationMatrixS_lowPileUP,variablesName,cutStringSignal); // calculate mean value for S and B
    calculateCovarianceMatrix(inputLowPileUpTreesBackground,correlationMatrixB_lowPileUP,variablesName,cutStringBackground); // calculate mean value for S and B

    TH2D* Matrix_S_lowPileUp = new TH2D(correlationMatrixS_lowPileUP); // make TH2D for plotting reasons
    Matrix_S_lowPileUp->SetName("Matrix_S_lowPileUp");
    TH2D* Matrix_B_lowPileUp = new TH2D(correlationMatrixB_lowPileUP);
    Matrix_B_lowPileUp->SetName("Matrix_B_lowPileUp");

    Matrix_S_lowPileUp->Scale(100); // correlation in %
    Matrix_B_lowPileUp->Scale(100);
  
    for( int iBinX = 0 ; iBinX < Matrix_S_lowPileUp->GetNbinsX() ; iBinX++){
     for( int iBinY = 0 ; iBinY < Matrix_S_lowPileUp->GetNbinsX() ; iBinY++){
      Matrix_S_lowPileUp->SetBinContent(iBinX+1,iBinY+1,int(Matrix_S_lowPileUp->GetBinContent(iBinX+1,iBinY+1)));
      Matrix_B_lowPileUp->SetBinContent(iBinX+1,iBinY+1,int(Matrix_B_lowPileUp->GetBinContent(iBinX+1,iBinY+1)));
      if(iBinX+1 == iBinY+1 ) {
       Matrix_S_lowPileUp->SetBinContent(iBinX+1,iBinY+1,100);
       Matrix_B_lowPileUp->SetBinContent(iBinX+1,iBinY+1,100);
      }
     }
    }

    ////////////////////////
    TCanvas* cCorrelationSignal = new TCanvas("CorrelationBDT_S",Form("Correlation Matrix Signal"),180,52,550,550);
    plotCovarianceMatrix(cCorrelationSignal,Matrix_S_lowPileUp,outputDirectory,reducedName);

    ////////////////////////
    TCanvas* cCorrelationBackground = new TCanvas("CorrelationBDT_B",Form("Correlation Matrix Background"),180,52,550,550);
    plotCovarianceMatrix(cCorrelationBackground,Matrix_B_lowPileUp,outputDirectory,reducedName);

    //////////////////////////////////////////////////////////////////////
    //// high pile-up correlation
    //////////////////////////////////////////////////////////////////////

    std::vector<TTree*> inputHighPileUpTreesSignal ;
    std::vector<TTree*> inputHighPileUpTreesBackground ;

    // take all the signal trees low pile up
    std::vector<std::string>::const_iterator itHighPileUp = InputInformationParamHighPUSignal.begin();
    for( ; itHighPileUp != InputInformationParamHighPUSignal.end() ; ++itHighPileUp){
     inputFile = TFile::Open((*itHighPileUp).c_str());
     if(inputFile == 0 or inputFile == NULL) continue;
     inputHighPileUpTreesSignal.push_back((TTree*) inputFile->Get(treeName.c_str())); // tree name should be provided
    }
    // take all the background trees low pile up
    itHighPileUp = InputInformationParamHighPUBackground.begin();
    for( ; itHighPileUp != InputInformationParamHighPUBackground.end() ; ++itHighPileUp){
     inputFile = TFile::Open((*itHighPileUp).c_str());
     if(inputFile == 0 or inputFile == NULL) continue;
     inputHighPileUpTreesBackground.push_back((TTree*) inputFile->Get(treeName.c_str())); // tree name should be provided
    }

    // Build Up the correlation matrix
    TMatrixDSym correlationMatrixS_highPileUP(reducedName.size()) ;
    TMatrixDSym correlationMatrixB_highPileUP(reducedName.size()) ;
    for( unsigned int irow = 0 ; irow < reducedName.size() ; irow++){
     for( unsigned int icolum = 0 ; icolum < reducedName.size() ; icolum++){
      correlationMatrixS_highPileUP(irow,icolum) = 0;
      correlationMatrixB_highPileUP(irow,icolum) = 0;
     }
    }

    calculateCovarianceMatrix(inputHighPileUpTreesSignal,correlationMatrixS_highPileUP,variablesName,cutStringSignal); // calculate mean value for S and B
    calculateCovarianceMatrix(inputHighPileUpTreesBackground,correlationMatrixB_highPileUP,variablesName,cutStringBackground); // calculate mean value for S and B

    TH2D* Matrix_S_highPileUp = new TH2D(correlationMatrixS_highPileUP);
    Matrix_S_highPileUp->SetName("Matrix_S_highPileUp");
    TH2D* Matrix_B_highPileUp = new TH2D(correlationMatrixB_highPileUP);
    Matrix_B_highPileUp->SetName("Matrix_B_highPileUp");

    Matrix_S_highPileUp->Scale(100); // correlation in %
    Matrix_B_highPileUp->Scale(100);

    for( int iBinX = 0 ; iBinX < Matrix_S_highPileUp->GetNbinsX() ; iBinX++){
     for( int iBinY = 0 ; iBinY < Matrix_S_highPileUp->GetNbinsX() ; iBinY++){
      Matrix_S_highPileUp->SetBinContent(iBinX+1,iBinY+1,int(Matrix_S_highPileUp->GetBinContent(iBinX+1,iBinY+1)));
      Matrix_B_highPileUp->SetBinContent(iBinX+1,iBinY+1,int(Matrix_B_highPileUp->GetBinContent(iBinX+1,iBinY+1)));
      if(iBinX+1 == iBinY+1 ) {
         Matrix_S_highPileUp->SetBinContent(iBinX+1,iBinY+1,100);
         Matrix_B_highPileUp->SetBinContent(iBinX+1,iBinY+1,100);
      }
     }
    }

    TCanvas* cCorrelationSignal_highPU = new TCanvas("CorrelationBDT_S_highPU",Form("Correlation Matrix Signal highPU"),180,52,550,550);
    plotCovarianceMatrix(cCorrelationSignal_highPU,Matrix_S_highPileUp,outputDirectory,reducedName);
    TCanvas* cCorrelationBackground_highPU = new TCanvas("CorrelationBDT_B_highPU",Form("Correlation Matrix Background"),180,52,550,550);
    plotCovarianceMatrix(cCorrelationBackground_highPU,Matrix_B_highPileUp,outputDirectory,reducedName);


    ////////////////////////////
    // Correlation difference among high pile-up bin and low pile-up bin for both S and B
    ////////////////////////////

    TH2D* CorrelationS_difference = new TH2D("CorrelationS_difference","",reducedName.size(),0,reducedName.size(),reducedName.size(),0,reducedName.size());
    TH2D* CorrelationB_difference = new TH2D("CorrelationB_difference","",reducedName.size(),0,reducedName.size(),reducedName.size(),0,reducedName.size());

    for(int iBinX = 0; iBinX < Matrix_S_highPileUp->GetNbinsX(); iBinX++){
     for(int iBinY = 0; iBinY < Matrix_S_highPileUp->GetNbinsY(); iBinY++){  
     CorrelationS_difference->SetBinContent(iBinX+1,iBinY+1,Matrix_S_highPileUp->GetBinContent(iBinX+1,iBinY+1)-Matrix_S_lowPileUp->GetBinContent(iBinX+1,iBinY+1));
     CorrelationB_difference->SetBinContent(iBinX+1,iBinY+1,Matrix_B_highPileUp->GetBinContent(iBinX+1,iBinY+1)-Matrix_B_lowPileUp->GetBinContent(iBinX+1,iBinY+1));
     }
    }

    TCanvas* cCorrelationSignal_difference = new TCanvas("cCorrelationSignal_difference",Form("Correlation Matrix Difference Signal"),180,52,550,550);
    plotCovarianceMatrix(cCorrelationSignal_difference,CorrelationS_difference,outputDirectory,reducedName);
    TCanvas* cCorrelationBackground_difference = new TCanvas("cCorrelationBackground_difference",Form("Correlation Matrix Difference Background"),180,52,550,550);
    plotCovarianceMatrix(cCorrelationBackground_difference,CorrelationS_difference,outputDirectory,reducedName);

  }

  return 0 ;

}

/// Calculate Mean Value
void calculateCovarianceMatrix (std::vector<TTree*> & inputTrees, TMatrixDSym & correlationMatrixS, TMatrixDSym & correlationMatrixB){

  std::vector<float> MeanValue_S; MeanValue_S.resize(inputTrees.size());
  std::vector<float> MeanValue_B; MeanValue_B.resize(inputTrees.size());

  std::vector<int> ientrySignal ; ientrySignal.resize(inputTrees.size());
  std::vector<int> ientryBackground ; ientryBackground.resize(inputTrees.size());

  int classID          = 0;
  float BDTG_NoPruning = 0;
  // compute mean values
  for(unsigned int iTree = 0 ; iTree<inputTrees.size(); iTree++){   
    inputTrees.at(iTree)->SetBranchAddress("classID",&classID);
    inputTrees.at(iTree)->SetBranchAddress("BDTG_NoPruning",&BDTG_NoPruning);
    for(int iEntries = 0; iEntries < inputTrees.at(iTree)->GetEntries(); iEntries++){
      inputTrees.at(iTree)->GetEntry(iEntries);
      if(classID == 0){ MeanValue_S.at(iTree) += BDTG_NoPruning; ientrySignal.at(iTree)++;}
      else{ MeanValue_B.at(iTree) += BDTG_NoPruning; ientryBackground.at(iTree)++;}
   }
  }

  for( unsigned int iVec = 0 ; iVec < MeanValue_S.size(); iVec++){
    MeanValue_S.at(iVec) = MeanValue_S.at(iVec)/ientrySignal.at(iVec);
    MeanValue_B.at(iVec) = MeanValue_B.at(iVec)/ientryBackground.at(iVec);
  }

  // compute distances
  int classID_X          = 0;
  float BDTG_NoPruning_X = 0;
  int classID_Y          = 0;
  float BDTG_NoPruning_Y = 0;

  for(unsigned int iTreeX = 0 ; iTreeX<inputTrees.size(); iTreeX++){   
    for(unsigned int iTreeY = 0 ; iTreeY<inputTrees.size(); iTreeY++){    
     int iEntriesX = 0; int iEntriesY = 0;
     for( ; iEntriesX < inputTrees.at(iTreeX)->GetEntries() and iEntriesY < inputTrees.at(iTreeY)->GetEntries(); iEntriesX++, iEntriesY++){
       inputTrees.at(iTreeX)->SetBranchAddress("classID",&classID_X);
       inputTrees.at(iTreeX)->SetBranchAddress("BDTG_NoPruning",&BDTG_NoPruning_X);
       inputTrees.at(iTreeX)->GetEntry(iEntriesX);
       inputTrees.at(iTreeY)->SetBranchAddress("classID",&classID_Y);
       inputTrees.at(iTreeY)->SetBranchAddress("BDTG_NoPruning",&BDTG_NoPruning_Y);
       inputTrees.at(iTreeY)->GetEntry(iEntriesY);
       if(classID_X == 0 and classID_Y == 0)
         correlationMatrixS(iTreeX,iTreeY) += double((BDTG_NoPruning_X-MeanValue_S.at(iTreeX))*(BDTG_NoPruning_Y-MeanValue_S.at(iTreeY)));
       else if(classID_X == 1 and classID_Y == 1) correlationMatrixB(iTreeX,iTreeY) += double((BDTG_NoPruning_X-MeanValue_B.at(iTreeX))*(BDTG_NoPruning_Y-MeanValue_B.at(iTreeY)));
     }
    }
  }
  
  for( int binX = 0; binX < correlationMatrixS.GetNrows(); binX ++){
    for( int binY = 0; binY < correlationMatrixS.GetNcols(); binY ++){
     correlationMatrixS(binX,binY) /= double(ientrySignal.at(binX));
     correlationMatrixB(binX,binY) /= double(ientryBackground.at(binX));
   }
  }
 
  std::vector<float> Sigma_S ; Sigma_S.resize(inputTrees.size());
  std::vector<float> Sigma_B ; Sigma_B.resize(inputTrees.size());
  for(int binX = 0; binX < correlationMatrixS.GetNrows(); binX ++) {
    Sigma_S.at(binX) = sqrt(correlationMatrixS(binX,binX));
    Sigma_B.at(binX) = sqrt(correlationMatrixB(binX,binX));
  }
  
  for( int binX = 0; binX < correlationMatrixS.GetNrows(); binX ++){
    for( int binY = 0; binY < correlationMatrixS.GetNcols(); binY ++){
     correlationMatrixS(binX,binY) /= double(Sigma_S.at(binX)*Sigma_S.at(binY));
     correlationMatrixB(binX,binY) /= double(Sigma_B.at(binX)*Sigma_B.at(binY));
   }     
  }


}


/// Calculate Mean Value
void calculateCovarianceMatrix ( std::vector<TTree*> & inputTrees, TMatrixDSym & correlationMatrix, std::vector<std::string> & variableName, std::string & cut){

  std::vector<float> MeanValue; MeanValue.resize(variableName.size());
  std::vector<int> ientry ; ientry.resize(variableName.size());

  // signal computation
  std::vector<TTreeFormula*> variablesFormula;
  std::vector<TTreeFormula*> cutFormula;
  
  for(size_t iTree = 0 ; iTree<inputTrees.size(); iTree++){   
    cutFormula.push_back(new TTreeFormula(cut.c_str(),cut.c_str(),inputTrees.at(iTree)));
    for(size_t iVar = 0; iVar < variableName.size(); iVar++){
      variablesFormula.push_back(new TTreeFormula (variableName.at(iVar).c_str(),variableName.at(iVar).c_str(),inputTrees.at(iTree)));
    }
    // loop on signal entries to compute mean values      
    for(int iEntries = 0; iEntries < inputTrees.at(iTree)->GetEntries(); iEntries++){
      inputTrees.at(iTree)->GetEntry(iEntries);
      if(cutFormula.at(iTree)->EvalInstance() == 0) continue;
      for(size_t iVar = 0; iVar < variableName.size(); iVar++){
        ientry.at(iVar)++;
	MeanValue.at(iVar) += variablesFormula.at(inputTrees.size()*iTree+iVar)->EvalInstance();
      }
    }
  }

  // get the mean value
  for( unsigned int iVec = 0 ; iVec < MeanValue.size(); iVec++)
    MeanValue.at(iVec) = MeanValue.at(iVec)/ientry.at(iVec);

  // compute distance for signal
  for(unsigned int iTree = 0 ; iTree<inputTrees.size(); iTree++){     
    // loop on entries 
    for(int iEntry = 0 ; iEntry < inputTrees.at(iTree)->GetEntries(); iEntry++){
      inputTrees.at(iTree)->GetEntry(iEntry);
      if(cutFormula.at(iTree)->EvalInstance() == 0) continue;
      for(unsigned int iVarX = 0 ; iVarX<variableName.size(); iVarX++){    
       for(unsigned int iVarY = 0 ; iVarY<variableName.size(); iVarY++){    
	correlationMatrix(iVarX,iVarY) += double((variablesFormula.at(iVarX+inputTrees.size()*iTree)->EvalInstance()-MeanValue.at(iVarX))*(variablesFormula.at(iVarY+inputTrees.size()*iTree)->EvalInstance()-MeanValue.at(iVarY)));
      }
     }
    }
  }

  for( int binX = 0; binX < correlationMatrix.GetNrows(); binX ++){
    for( int binY = 0; binY < correlationMatrix.GetNcols(); binY ++){
     correlationMatrix(binX,binY) /= double(ientry.at(binX));
   }
  }
 
  std::vector<float> Sigma ; Sigma.resize(variableName.size());
  for(int binX = 0; binX < correlationMatrix.GetNrows(); binX ++) {
    Sigma.at(binX) = sqrt(correlationMatrix(binX,binX));
  }
  
  for( int binX = 0; binX < correlationMatrix.GetNrows(); binX ++){
    for( int binY = 0; binY < correlationMatrix.GetNcols(); binY ++){
     correlationMatrix(binX,binY) /= double(Sigma.at(binX)*Sigma.at(binY));
   }     
  }

}

///////////////
void plotCovarianceMatrix(TCanvas* cCorrelation, TH2D* Matrix, std::string outputDirectory, std::vector<std::string> LabelNames){

  float newMargin1 = 0.13;
  float newMargin2 = 0.15;
  float newMargin3 = 0.20;

  cCorrelation->SetGrid();
  cCorrelation->SetTicks();
  cCorrelation->SetLeftMargin(newMargin3);
  cCorrelation->SetBottomMargin(newMargin2);
  cCorrelation->SetRightMargin(newMargin1);
  cCorrelation->cd();
  gStyle->SetPaintTextFormat("3g");

  Matrix->SetMarkerSize(1.5);
  Matrix->SetMarkerColor(0);
  Matrix->GetXaxis()->SetLabelSize(0.035);
  Matrix->GetYaxis()->SetLabelSize(0.035);
  Matrix->SetLabelOffset(0.011);// label offset on x axis                                                                                                                   

  Matrix->SetMaximum(100);         
  Matrix->SetMinimum(-100);         

  Matrix->Draw("colz");
  for( int binX = 0 ; binX < Matrix->GetNbinsX(); binX++){
    Matrix->GetXaxis()->SetBinLabel(binX+1,LabelNames.at(binX).c_str());
    Matrix->GetYaxis()->SetBinLabel(binX+1,LabelNames.at(binX).c_str());
  }
  Matrix->Draw("textsame");
  Matrix->GetXaxis()->LabelsOption("v");
  Matrix->GetYaxis()->LabelsOption("h");

  TLatex *   tex = new TLatex(0.85,0.92," 13 TeV");
  tex->SetNDC();
  tex->SetTextAlign(31);
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(0.13,0.92,"CMS");
  tex->SetNDC();
  tex->SetTextFont(61);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(0.2324,0.92,"Preliminary Simulation");
  tex->SetNDC();
  tex->SetTextFont(52);
  tex->SetTextSize(0.0304);
  tex->SetLineWidth(2);
  tex->Draw();
  
  cCorrelation->Print((outputDirectory+"/"+std::string(cCorrelation->GetName())+".pdf").c_str(),"pdf");
  cCorrelation->Print((outputDirectory+"/"+std::string(cCorrelation->GetName())+".png").c_str(),"png");
  cCorrelation->Print((outputDirectory+"/"+std::string(cCorrelation->GetName())+".root").c_str(),"root");
}

///////////////////////////////////////////////////

TMatrixDSym eigenTransform (TMatrixDSym& Matrix){

  //copy the original matrix in a new one skipping the last row and columns --> correlation with the full BDT
  TMatrixDSym lSMatrix(Matrix.GetNrows()-1);
  for(int iLine = 0; iLine < Matrix.GetNrows()-1; iLine++) { 
    for(int iRow = 0; iRow < Matrix.GetNcols()-1; iRow++) { 
      lSMatrix(iLine,iRow) = Matrix(iLine,iRow);
    }
  }

  TVectorD lEigVals(Matrix.GetNrows()-1); // make a vector were store the eigenvalues of the diagonalized matrix
  lEigVals = TMatrixDSymEigen(lSMatrix).GetEigenValues();
  TMatrixD lEigen = lSMatrix.EigenVectors(lEigVals); // take the eigenvector matrix

  TMatrixDSym lEigenSym(Matrix.GetNrows()); // new matrix loop on the N-1 x N-1 inner line and columns and copy the eigenVectors
  for(int iLine   = 0; iLine < Matrix.GetNrows()-1; iLine++) { 
    for(int iRow = 0; iRow < Matrix.GetNcols()-1; iRow++) { 
      lEigenSym(iLine,iRow) = lEigen(iLine,iRow);
    }
  }
  for(int iLine   = 0; iLine < Matrix.GetNrows()-1; iLine++) lEigenSym(Matrix.GetNrows()-1,iLine) = lEigVals(iLine); // put in the N-1 line all the eigenvalues
  return lEigenSym;

}

///////////////////////////////////////////////////
TMatrixDSym dimOrder(TMatrixDSym & Matrix,TMatrixDSym & Eigen,TMatrixD & EigenInv) {

  TVectorD lVector(Matrix.GetNrows()-1); 
  for(int iLine = 0; iLine < Matrix.GetNrows()-1; iLine++) lVector(iLine) = Matrix(Matrix.GetNrows()-1,iLine); // take the original values of the N-1 line
  double lMag = 0; 
  for(int iLine = 0; iLine < Matrix.GetNrows()-1; iLine++) lMag += lVector(iLine)*lVector(iLine); 
  for(int iLine = 0; iLine < Matrix.GetNrows()-1; iLine++) lVector(iLine)  = lVector(iLine)/sqrt(lMag); 
  lVector = EigenInv * lVector; // multiply this vector for the inverse of the eigenVector matrix
  std::vector<double> lSum; 
  for(int iLine = 0; iLine < Matrix.GetNrows()-1; iLine++) lSum.push_back(lVector(iLine)*lVector(iLine)); // vector of squared of eigenValues
  lMag = 0; 
  for(int iLine = 0; iLine < Matrix.GetNrows()-1; iLine++) lMag += lVector(iLine)*lVector(iLine); // new magnitude

  std::vector<int>  lRank; // ranking vector
  for(int iLine = 0; iLine < Matrix.GetNrows()-1; iLine++) {
    int pRank = 0;
    for(int iCol = 0; iCol < Matrix.GetNcols()-1; iCol++){
      if(lSum[iCol] > lSum[iLine]) pRank++;
    }
    lRank.push_back(pRank);
  }

  TMatrixDSym lOutMatrix(Matrix.GetNrows());
  for(int iLine = 0; iLine < Matrix.GetNrows(); iLine++) {
    for(int iCol = 0; iCol < Matrix.GetNcols(); iCol++) {
      lOutMatrix(iLine,iCol) = 0;
    }
  }
  for(int iLine = 0; iLine < Matrix.GetNrows()-1; iLine++) {
    for(int iCol = 0; iCol < Matrix.GetNcols()-1; iCol++) {
      lOutMatrix(iLine,lRank[iCol]) = Eigen(iLine,iCol);
    }
  }
  for(int iLine   = 0; iLine < Matrix.GetNrows()-1; iLine++) lOutMatrix(Matrix.GetNrows()-1,lRank[iLine]) = lVector(iLine)/sqrt(lMag);  
  return lOutMatrix;
}


TMatrixDSym eigenCov(TMatrixDSym &iEigen,std::vector<std::string> variableName, std::vector<TTree*> inputTrees, std::string iCut) {

  std::vector<std::string> lOrthoVars; // orthogonal variables

  for(int iLine = 0; iLine < iEigen.GetNrows()-1; iLine++) { 
    std::stringstream pVec;
    for(int iCol = 0; iCol < iEigen.GetNcols(); iCol++) { 
      pVec << iEigen(iCol,iLine) << "*" << variableName.at(iCol);
      if(iCol != iEigen.GetNcols()-2) pVec << " + ";
    } 
    lOrthoVars.push_back(pVec.str());
  }

  lOrthoVars.push_back("bdt_all");
  TMatrixDSym lEigenSym(variableName.size());
  calculateCovarianceMatrix(inputTrees,lEigenSym,lOrthoVars,iCut);

  for(int iLine   = 0; iLine < lEigenSym.GetNrows()-1; iLine++) {
    for(int iCol = 0; iCol < lEigenSym.GetNcols(); iCol++) {
      lEigenSym(iLine,iCol) = iEigen(iLine,iCol);
    }
  }

  for(int iLine = 0; iLine < lEigenSym.GetNrows(); iLine++) lEigenSym(iLine,lEigenSym.GetNrows()-1) = 0; 

  std::vector<int> lRank;
  for(int iLine = 0; iLine < lEigenSym.GetNrows()-1; iLine++) {
    int pRank = 0;
    for(int iCol = 0; iCol < lEigenSym.GetNcols()-1; iCol++) if(fabs(lEigenSym(lEigenSym.GetNrows()-1,iCol)) > fabs(lEigenSym(lEigenSym.GetNrows()-1,iLine))) pRank++;
    lRank.push_back(pRank);
  }

  TMatrixDSym lOutMatrix(lEigenSym.GetNrows());
  for(int iLine     = 0; iLine < lEigenSym.GetNrows(); iLine++) for(int iCol   = 0; iCol < lEigenSym.GetNcols(); iCol++)  lOutMatrix(iLine,iCol) = 0;
  for(int iLine   = 0; iLine < lEigenSym.GetNrows(); iLine++) {
    for(int iCol = 0; iCol < lEigenSym.GetNcols()-1; iCol++) {
      lOutMatrix(iLine,lRank[iCol]) = lEigenSym(iLine,iCol);
    }
  }
  return lOutMatrix;
}

