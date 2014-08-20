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

void calculateCovarianceMatrix ( std::vector<TTree*> & inputTrees, TMatrixDSym & correlationMatrixS, TMatrixDSym & correlationMatrixB );
void plotCovarianceMatrix      ( TCanvas* cCorrelation, TH2D* Matrix, std::string outputDirectory, std::vector<std::string> labelNames);
TMatrixDSym eigenTransform (TMatrixDSym& Matrix);
TMatrixDSym dimOrder(TMatrixDSym & Matrix,TMatrixDSym & Eigen,TMatrixD & EigenInv);

/// Main programme 
int main (int argc, char **argv){
  if(argc<2){ std::cout<<" Not correct number of input parameter --> Need Just one cfg file exit "<<std::endl; return -1; }

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


  // parse config file parameter                                                                                                                                                          
  std::string configFileName = argv[1];
  boost::shared_ptr<edm::ParameterSet> parameterSet = edm::readConfig(configFileName);
  edm::ParameterSet Options  = parameterSet -> getParameter<edm::ParameterSet>("Options");


  std::vector<edm::ParameterSet> InputInformationParamLowPU ;
  if(Options.existsAs<std::vector<edm::ParameterSet>>("InputLowPUFiles"))
    InputInformationParamLowPU = Options.getParameter<std::vector<edm::ParameterSet>>("InputLowPUFiles");
  else{ std::cout<<" Exit from code, no input set found for low pile-up files"<<std::endl; return -1; }

  std::vector<edm::ParameterSet> InputInformationParamHighPU ;
  if(Options.existsAs<std::vector<edm::ParameterSet>>("InputHighPUFiles"))
  InputInformationParamHighPU = Options.getParameter<std::vector<edm::ParameterSet>>("InputHighPUFiles");
  else{ std::cout<<" Exit from code, no input set found for high pile-up files"<<std::endl; return -1; }


  std::string outputDirectory;
  if(Options.existsAs<std::string>("outputDirectory"))
    outputDirectory = Options.getParameter<std::string>("outputDirectory");
  else{ outputDirectory = "output/"; }
 
  system(("mkdir -p "+outputDirectory).c_str());

  std::vector<TTree*> inputLowPileUpTrees ;
  std::vector<std::string> reducedNameLowPileUp ;
  TFile* inputFile = NULL ;
  // take all the TH1 outputs for S and B for  each variable in input-> low Pile Up
  std::vector<edm::ParameterSet>::const_iterator itLowPileUp = InputInformationParamLowPU.begin();
  for( ; itLowPileUp != InputInformationParamLowPU.end() ; ++itLowPileUp){
    inputFile = TFile::Open((*itLowPileUp).getParameter<std::string>("fileName").c_str());
    if(inputFile == 0 or inputFile == NULL) continue;
    inputLowPileUpTrees.push_back((TTree*) inputFile->Get("TestTree"));
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

  TH2D* Matrix_S_lowPileUp = new TH2D(correlationMatrixS_lowPileUP);
  Matrix_S_lowPileUp->SetName("Matrix_S_lowPileUp");
  TH2D* Matrix_B_lowPileUp = new TH2D(correlationMatrixB_lowPileUP);
  Matrix_B_lowPileUp->SetName("Matrix_B_lowPileUp");

  Matrix_S_lowPileUp->Scale(100);
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

  Matrix_S_highPileUp->Scale(100);
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

  /*    
  TH2D* Matrix_S_highPileUp = new TH2D(correlationMatrixS_lowPileUP);
  Matrix_S_highPileUp->SetName("Matrix_S_highPileUp");
  TH2D* Matrix_B_highPileUp = new TH2D(correlationMatrixB_lowPileUP);
  Matrix_B_highPileUp->SetName("Matrix_B_highPileUp");

  Matrix_S_highPileUp->Scale(100);
  Matrix_B_highPileUp->Scale(100);
  
  TRandom3 rand ;
  for( int iBinX = 0 ; iBinX < Matrix_S_highPileUp->GetNbinsX() ; iBinX++){
   for( int iBinY = iBinX ; iBinY < Matrix_S_highPileUp->GetNbinsX() ; iBinY++){
     double svalue = int(rand.Gaus(Matrix_S_highPileUp->GetBinContent(iBinX+1,iBinY+1),5));
     double bvalue = int(rand.Gaus(Matrix_B_highPileUp->GetBinContent(iBinX+1,iBinY+1),5));
     Matrix_S_highPileUp->SetBinContent(iBinX+1,iBinY+1,svalue);
     Matrix_S_highPileUp->SetBinContent(iBinY+1,iBinX+1,svalue);
     Matrix_B_highPileUp->SetBinContent(iBinX+1,iBinY+1,bvalue);
     Matrix_B_highPileUp->SetBinContent(iBinY+1,iBinX+1,bvalue);
     if(iBinX == iBinY){
      Matrix_S_highPileUp->SetBinContent(iBinX+1,iBinY+1,100);
      Matrix_B_highPileUp->SetBinContent(iBinX+1,iBinY+1,100);
     }
   }
  }
  */

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

  //Decomopose Covariance matrix to orthogonal states
  TMatrixDSym lSMatrix(Matrix.GetNrows()-1);
  for(int iLine = 0; iLine < Matrix.GetNrows()-1; iLine++) { 
    for(int iRow = 0; iRow < Matrix.GetNcols()-1; iRow++) { 
      lSMatrix(iLine,iRow) = Matrix(iLine,iRow);
    }
  }

  TVectorD lEigVals(Matrix.GetNrows()-1);
  lEigVals = TMatrixDSymEigen(lSMatrix).GetEigenValues();
  TMatrixD lEigen = lSMatrix.EigenVectors(lEigVals);

  TMatrixDSym lEigenSym(Matrix.GetNrows());
  double pCheck = 0;
  for(int iLine   = 0; iLine < Matrix.GetNrows()-1; iLine++) { 
    for(int iRow = 0; iRow < Matrix.GetNcols()-1; iRow++) { 
      if(iRow == 2) pCheck +=  lEigen(iLine,iRow)*lEigen(iLine,iRow);
      lEigenSym(iLine,iRow) = lEigen(iLine,iRow);
    }
  }
  for(int iLine   = 0; iLine < Matrix.GetNrows()-1; iLine++) lEigenSym(Matrix.GetNrows()-1,iLine) = lEigVals(iLine);
  return lEigenSym;

}

///////////////////////////////////////////////////
TMatrixDSym dimOrder(TMatrixDSym & Matrix,TMatrixDSym & Eigen,TMatrixD & EigenInv) {

  TVectorD lVector(Matrix.GetNrows()-1);
  for(int iLine = 0; iLine < Matrix.GetNrows()-1; iLine++) lVector(iLine) = Matrix(Matrix.GetNrows()-1,iLine); // take the eigenvalues
  double lMag = 0; 
  for(int iLine = 0; iLine < Matrix.GetNrows()-1; iLine++)  lMag += lVector(iLine)*lVector(iLine); // sum of eigenvalues squared
  for(int iLine = 0; iLine < Matrix.GetNrows()-1; iLine++)  lVector(iLine)  = lVector(iLine)/sqrt(lMag); // normalize the eigenvalues to be [0,1] in abs value
  lVector = EigenInv * lVector; // multiply this vector for the inverse of the eigenVector matrix
  std::vector<double> lSum; 
  for(int iLine = 0; iLine < Matrix.GetNrows()-1; iLine++) lSum.push_back(lVector(iLine)*lVector(iLine)); // vector of squared of eigenValues
  lMag = 0; 
  for(int iLine = 0; iLine < Matrix.GetNrows()-1; iLine++)  lMag += lVector(iLine)*lVector(iLine); // new magnitude

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
