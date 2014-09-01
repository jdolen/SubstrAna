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
#include "TGraphErrors.h"

#include "TKey.h"
#include "TObjArray.h"
#include "TClass.h"
#include "TH2F.h"
#include "TMatrixDSym.h"
#include "TProfile2D.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

// fill a TList with all the methods                                                                                                                                                    
int GetListOfMethods( TList & methods, TDirectory *dir = 0);
int GetListOfTitles(  TDirectory*rfdir,TList & titles);
TKey *NextKey(TIter & keyIter, TString className);

struct MatrixEntry{

  std::string binXName ; 
  std::string binYName ;
  double value;
  double error;
};

std::vector<MatrixEntry> rocMatrix(std::vector<TTree*> InputTreeSignal,std::vector<TTree*> InputTreeBackground, std::vector<std::string> branchName, std::vector<std::string> reduceNameX,
                                   std::vector<std::string> reduceNameY, std::string cutSignalLowPileUp, std::string cutBackgroundLowPileUp, double signalEfficiencyTarget);



/// Main programme 
int main (int argc, char **argv){

  if(argc < 3){ std::cout<<" Not correct number of input parameter --> Need Just one cfg and one number to specify the kind of output file --> exit "<<std::endl; return -1; }

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
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetErrorX(0.5);
  gStyle->SetPaintTextFormat("3g");

  // parse config file parameter                                                                                                                                                          
  std::string configFileName = argv[1];
  boost::shared_ptr<edm::ParameterSet> parameterSet = edm::readConfig(configFileName);
  edm::ParameterSet Options  = parameterSet -> getParameter<edm::ParameterSet>("Options");

  if(atoi(argv[2]) == 0){

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

   double signalEfficiencyTarget;
   if(Options.existsAs<double>("signalEfficiencyTarget"))
    signalEfficiencyTarget = Options.getParameter<double>("signalEfficiencyTarget");
   else{ signalEfficiencyTarget = 0.5; }
 
   system(("mkdir -p "+outputDirectory).c_str());

   std::vector<MatrixEntry> performanceValue ;

   TFile* inputFile = NULL ;

   // take all the TH1 outputs for S and B for  each variable in input-> low Pile Up
   TList TrainingMethods;
   TList Titles;

   std::vector<edm::ParameterSet>::const_iterator itLowPileUp = InputInformationParamLowPU.begin();
   for( ; itLowPileUp != InputInformationParamLowPU.end() ; ++itLowPileUp){
     inputFile = TFile::Open((*itLowPileUp).getParameter<std::string>("fileName").c_str());
     if(inputFile == 0 or inputFile == NULL){
      MatrixEntry entry ; 
      entry.binXName = (*itLowPileUp).getParameter<std::string>("variableNameX");
      entry.binYName = (*itLowPileUp).getParameter<std::string>("variableNameY");
      entry.value = 0;   
      entry.error = 0;
      performanceValue.push_back(entry);         
      continue;
     }
     // take the background efficiency related to the target
     inputFile->cd();
     TrainingMethods.Clear();
     int res = GetListOfMethods(TrainingMethods);
     if(res == 0) std::cout<<" No methods found "<<std::endl ;
     TIter next(&TrainingMethods);
     TKey *key = 0, *hkey = 0;
     while ((key = (TKey*)next())) {
      TDirectory *myDir = (TDirectory*)key->ReadObj();
      Titles.Clear();
      int nTitles = GetListOfTitles(myDir,Titles);
      if(nTitles == 0) std::cout<<" No titles found "<<std::endl ;
      TIter nextTitle(&Titles);
      TKey *titkey = 0;
      TDirectory *titDir = 0;
      while ((titkey = NextKey(nextTitle,"TDirectory"))) {
	titDir = (TDirectory*)titkey->ReadObj(); // read each object and take again the method title for each element of the list                                                         
	TString methodTitle;
        methodTitle = titDir->GetName();
	TIter nextKey( titDir->GetListOfKeys() ); // loop and the list of keys                                                                                                
	while ((hkey = NextKey(nextKey,"TH1"))) { // take only the TH1 object type                                                                                                   
	  TH1F* h = (TH1F*) hkey->ReadObj();
	  TString hname = h->GetName();    // only the one which are called rejBvsS            
	  if (hname.Contains("effBvsS") && hname.BeginsWith("MVA_") && not hname.Contains("effBvsSLocal")) {
            MatrixEntry entry ; 
            entry.binXName = (*itLowPileUp).getParameter<std::string>("variableNameX");
            entry.binYName = (*itLowPileUp).getParameter<std::string>("variableNameY");
	    entry.value = (1./h->GetBinContent(h->FindBin(signalEfficiencyTarget)));    
            entry.error = (1/h->GetBinError(h->FindBin(signalEfficiencyTarget)));
            performanceValue.push_back(entry);
	  }
	}
      }
    }    
   }

   //For getting the number of bins is enough to cycle on all the entry and count the diagonal terms
   int numberOfBins = 0;
   for(unsigned int ientry = 0 ; ientry < performanceValue.size(); ientry++){
     if(performanceValue.at(ientry).binXName == performanceValue.at(ientry).binYName) numberOfBins++;
   }

   TH2F* performanceBDT_lowPileUP  = new TH2F("backgroundMatrix_lowPileUP","",numberOfBins,0,numberOfBins,numberOfBins,0,numberOfBins);
   TH2F* performanceBDT_lowPileUP_float  = new TH2F("backgroundMatrix_lowPileUP_float","",numberOfBins,0,numberOfBins,numberOfBins,0,numberOfBins);
   TH2F* performanceBDT_lowPileUP_error  = new TH2F("backgroundMatrix_lowPileUP_error","",numberOfBins,0,numberOfBins,numberOfBins,0,numberOfBins);

   for(int iBinX = 0; iBinX < performanceBDT_lowPileUP->GetNbinsX(); iBinX ++){
    for(int iBinY = 0; iBinY < performanceBDT_lowPileUP->GetNbinsY(); iBinY ++){
     performanceBDT_lowPileUP->SetBinContent(iBinX+1,iBinY+1,0.);
     performanceBDT_lowPileUP_float->SetBinContent(iBinX+1,iBinY+1,0.);
     performanceBDT_lowPileUP_error->SetBinContent(iBinX+1,iBinY+1,0.);
    }
   }

   int vecPos = 0;
   for(int iBinX = 0; iBinX < performanceBDT_lowPileUP->GetNbinsX(); iBinX ++){  
    for(int iBinY = iBinX; iBinY < performanceBDT_lowPileUP->GetNbinsY(); iBinY ++){
     if(iBinY < (performanceBDT_lowPileUP->GetNbinsY()-1) and iBinX < (performanceBDT_lowPileUP->GetNbinsX()-1) and vecPos < int(performanceValue.size()-1)){

       performanceBDT_lowPileUP_float->SetBinContent(iBinX+1,iBinY+1,performanceValue.at(vecPos).value);
       if(performanceValue.at(vecPos).value - int(performanceValue.at(vecPos).value) < 0.5)
        performanceBDT_lowPileUP->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.at(vecPos).value));
       else
        performanceBDT_lowPileUP->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.at(vecPos).value)+1);
       if(performanceValue.at(vecPos).error - int(performanceValue.at(vecPos).error) < 0.5)
        performanceBDT_lowPileUP_error->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.at(vecPos).error));
       else
        performanceBDT_lowPileUP_error->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.at(vecPos).error)+1);

        performanceBDT_lowPileUP->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.at(vecPos).value));
        performanceBDT_lowPileUP_error->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.at(vecPos).error));

       if(std::string(performanceBDT_lowPileUP->GetYaxis()->GetBinLabel(iBinY+1)) == ""){
	 performanceBDT_lowPileUP->GetYaxis()->SetBinLabel(iBinY+1,(performanceValue.at(vecPos).binYName).c_str());
         performanceBDT_lowPileUP_error->GetYaxis()->SetBinLabel(iBinY+1,(performanceValue.at(vecPos).binYName).c_str());
       }
       if(std::string(performanceBDT_lowPileUP->GetXaxis()->GetBinLabel(iBinY+1)) == ""){
         performanceBDT_lowPileUP->GetXaxis()->SetBinLabel(iBinY+1,(performanceValue.at(vecPos).binYName).c_str());
         performanceBDT_lowPileUP_error->GetXaxis()->SetBinLabel(iBinY+1,(performanceValue.at(vecPos).binYName).c_str());
       }
       vecPos++; 
     }    
     else if (iBinY == (performanceBDT_lowPileUP->GetNbinsY()-1) and iBinX < (performanceBDT_lowPileUP->GetNbinsX()-1)) continue ;
     else if (iBinY < (performanceBDT_lowPileUP->GetNbinsY()-1) and iBinX == (performanceBDT_lowPileUP->GetNbinsX()-1)) continue ;
     else if (iBinY == (performanceBDT_lowPileUP->GetNbinsY()-1) and iBinX == (performanceBDT_lowPileUP->GetNbinsX()-1)){
       performanceBDT_lowPileUP_float->SetBinContent(iBinX+1,iBinY+1,performanceValue.at(vecPos).value);
       if(performanceValue.back().value-int(performanceValue.back().value)<0.5)
        performanceBDT_lowPileUP->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.back().value));
       else
        performanceBDT_lowPileUP->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.back().value)+1);

       performanceBDT_lowPileUP->GetXaxis()->SetBinLabel(iBinX+1,(performanceValue.back().binXName).c_str());
       performanceBDT_lowPileUP->GetYaxis()->SetBinLabel(iBinY+1,(performanceValue.back().binXName).c_str());
       if(performanceValue.back().error-int(performanceValue.back().error)<0.5)
        performanceBDT_lowPileUP_error->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.back().error));
       else  
        performanceBDT_lowPileUP_error->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.back().error)+1);
       performanceBDT_lowPileUP_error->GetXaxis()->SetBinLabel(iBinX+1,(performanceValue.back().binXName).c_str());
       performanceBDT_lowPileUP_error->GetYaxis()->SetBinLabel(iBinY+1,(performanceValue.back().binXName).c_str());
     }
    }
   }

   TCanvas* cPerformance_lowPU = new TCanvas("cPerformance_lowPU","",180,52,500,550);

   cPerformance_lowPU->SetGrid();
   cPerformance_lowPU->SetTicks();
   cPerformance_lowPU->cd();

   performanceBDT_lowPileUP->SetMarkerSize(1.5);
   performanceBDT_lowPileUP->SetMarkerColor(0);
   performanceBDT_lowPileUP->GetXaxis()->SetLabelSize(0.035);
   performanceBDT_lowPileUP->GetYaxis()->SetLabelSize(0.035);
   performanceBDT_lowPileUP->SetLabelOffset(0.011);// label offset on x axis                                                                                                  
   performanceBDT_lowPileUP->SetMaximum(performanceBDT_lowPileUP->GetMaximum());
   performanceBDT_lowPileUP->SetMinimum(performanceBDT_lowPileUP->GetMinimum());
   performanceBDT_lowPileUP->Draw("colz");
   performanceBDT_lowPileUP->Draw("textsame");
   performanceBDT_lowPileUP->GetXaxis()->LabelsOption("v");
   performanceBDT_lowPileUP->GetYaxis()->LabelsOption("h");


   TLatex* tex = new TLatex(0.85,0.92," 13 TeV");
   tex->SetNDC();
   tex->SetTextAlign(31);
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
   TLatex* tex2 = new TLatex(0.13,0.92,"CMS");
   tex2->SetNDC();
   tex2->SetTextFont(61);
   tex2->SetTextSize(0.04);
   tex2->SetLineWidth(2);
   tex2->Draw();
   TLatex* tex3 = new TLatex(0.2324,0.92,"Preliminary Simulation");
   tex3->SetNDC();
   tex3->SetTextFont(52);
   tex3->SetTextSize(0.0304);
   tex3->SetLineWidth(2);
   tex3->Draw();

   cPerformance_lowPU->Print((outputDirectory+"/PerformanceBDT_lowPU.pdf").c_str(),"pdf");
   cPerformance_lowPU->Print((outputDirectory+"/PerformanceBDT_lowPU.png").c_str(),"png");
   cPerformance_lowPU->Print((outputDirectory+"/PerformanceBDT_lowPU.root").c_str(),"root");
   cPerformance_lowPU->SetLogz();
   cPerformance_lowPU->Print((outputDirectory+"/PerformanceBDT_lowPU_Log.pdf").c_str(),"pdf");
   cPerformance_lowPU->Print((outputDirectory+"/PerformanceBDT_lowPU_Log.png").c_str(),"png");
   cPerformance_lowPU->Print((outputDirectory+"/PerformanceBDT_lowPU_Log.root").c_str(),"root");

   // Error
   TCanvas* cPerformance_lowPU_error = new TCanvas("cPerformance_lowPU_error","",180,52,500,550);

   cPerformance_lowPU_error->SetGrid();
   cPerformance_lowPU_error->SetTicks();
   cPerformance_lowPU_error->cd();

   performanceBDT_lowPileUP_error->SetMarkerSize(1.5);
   performanceBDT_lowPileUP_error->SetMarkerColor(0);
   performanceBDT_lowPileUP_error->GetXaxis()->SetLabelSize(0.035);
   performanceBDT_lowPileUP_error->GetYaxis()->SetLabelSize(0.035);
   performanceBDT_lowPileUP_error->SetLabelOffset(0.011);// label offset on x axis                                                                                              
   performanceBDT_lowPileUP_error->SetMaximum(performanceBDT_lowPileUP_error->GetMaximum());
   performanceBDT_lowPileUP_error->SetMinimum(performanceBDT_lowPileUP_error->GetMinimum());
   performanceBDT_lowPileUP_error->Draw("colz");
   performanceBDT_lowPileUP_error->Draw("textsame");
   performanceBDT_lowPileUP_error->GetXaxis()->LabelsOption("v");
   performanceBDT_lowPileUP_error->GetYaxis()->LabelsOption("h");

   tex->Draw();
   tex2->Draw();
   tex3->Draw();

   cPerformance_lowPU_error->Print((outputDirectory+"/PerformanceBDT_lowPU_error.pdf").c_str(),"pdf");
   cPerformance_lowPU_error->Print((outputDirectory+"/PerformanceBDT_lowPU_error.png").c_str(),"png");
   cPerformance_lowPU_error->Print((outputDirectory+"/PerformanceBDT_lowPU_error.root").c_str(),"root");
   cPerformance_lowPU_error->SetLogz();
   cPerformance_lowPU_error->Print((outputDirectory+"/PerformanceBDT_lowPU_Log_error.pdf").c_str(),"pdf");
   cPerformance_lowPU_error->Print((outputDirectory+"/PerformanceBDT_lowPU_Log_error.png").c_str(),"png");
   cPerformance_lowPU_error->Print((outputDirectory+"/PerformanceBDT_lowPU_Log_error.root").c_str(),"root");

   // high PU part
   performanceValue.clear();
   std::vector<edm::ParameterSet>::const_iterator itHighPileUp = InputInformationParamHighPU.begin();
   for( ; itHighPileUp != InputInformationParamHighPU.end() ; ++itHighPileUp){
    inputFile = TFile::Open((*itHighPileUp).getParameter<std::string>("fileName").c_str());
    if(inputFile == 0 or inputFile == NULL){
      MatrixEntry entry ; 
      entry.binXName = (*itHighPileUp).getParameter<std::string>("variableNameX");
      entry.binYName = (*itHighPileUp).getParameter<std::string>("variableNameY");
      entry.value = 0;   
      entry.error = 0;
      performanceValue.push_back(entry);         
      continue;
    }
    // take the background efficiency related to the target
    inputFile->cd();
    TrainingMethods.Clear();
    int res = GetListOfMethods(TrainingMethods);
    if(res == 0) std::cout<<" No methods found "<<std::endl ;
    TIter next(&TrainingMethods);
    TKey *key = 0, *hkey = 0;
    while ((key = (TKey*)next())) {
      TDirectory *myDir = (TDirectory*)key->ReadObj();
      Titles.Clear();
      int nTitles = GetListOfTitles(myDir,Titles);
      if(nTitles == 0) std::cout<<" No titles found "<<std::endl ;
      TIter nextTitle(&Titles);
      TKey *titkey = 0;
      TDirectory *titDir = 0;
      while ((titkey = NextKey(nextTitle,"TDirectory"))) {
	titDir = (TDirectory*)titkey->ReadObj(); // read each object and take again the method title for each element of the list                                                         
	TString methodTitle;
        methodTitle = titDir->GetName();
	TIter nextKey( titDir->GetListOfKeys() ); // loop and the list of keys                                                                                                
	while ((hkey = NextKey(nextKey,"TH1"))) { // take only the TH1 object type                                                                                                   
	  TH1F* h = (TH1F*) hkey->ReadObj();
	  TString hname = h->GetName();    // only the one which are called rejBvsS            
	  if (hname.Contains("effBvsS") && hname.BeginsWith("MVA_") && not hname.Contains("effBvsSLocal")) {
            MatrixEntry entry ; 
            entry.binXName = (*itHighPileUp).getParameter<std::string>("variableNameX");
            entry.binYName = (*itHighPileUp).getParameter<std::string>("variableNameY");
	    entry.value = (1./h->GetBinContent(h->FindBin(signalEfficiencyTarget)));   
            entry.error = (1./h->GetBinError(h->FindBin(signalEfficiencyTarget)));
            performanceValue.push_back(entry);         
	  }
	}
      }
    }
   }

   //for getting the number of bins is enough to cycle on all the entry and count the diagonal terms
   numberOfBins = 0;
   for(unsigned int ientry = 0 ; ientry < performanceValue.size(); ientry++){
    if(performanceValue.at(ientry).binXName == performanceValue.at(ientry).binYName) numberOfBins++;
   }

   TH2F* performanceBDT_highPileUP        = new TH2F("backgroundMatrix_highPileUP","",numberOfBins,0,numberOfBins,numberOfBins,0,numberOfBins);
   TH2F* performanceBDT_highPileUP_float  = new TH2F("backgroundMatrix_highPileUP_float","",numberOfBins,0,numberOfBins,numberOfBins,0,numberOfBins);
   TH2F* performanceBDT_highPileUP_error  = new TH2F("backgroundMatrix_highPileUP_error","",numberOfBins,0,numberOfBins,numberOfBins,0,numberOfBins);
   TH2F* performanceDifference            = new TH2F("performanceDifference","",numberOfBins,0,numberOfBins,numberOfBins,0,numberOfBins);
   TH2F* performanceRatio                 = new TH2F("performanceRatio","",numberOfBins,0,numberOfBins,numberOfBins,0,numberOfBins);

   for(int iBinX = 0; iBinX < performanceBDT_highPileUP->GetNbinsX(); iBinX ++){
    for(int iBinY = 0; iBinY < performanceBDT_highPileUP->GetNbinsY(); iBinY ++){
     performanceBDT_highPileUP->SetBinContent(iBinX+1,iBinY+1,0.);
     performanceBDT_highPileUP_error->SetBinContent(iBinX+1,iBinY+1,0.);
    }
   }

   vecPos = 0;
   for(int iBinX = 0; iBinX < performanceBDT_highPileUP->GetNbinsX(); iBinX ++){
    for(int iBinY = iBinX; iBinY < performanceBDT_highPileUP->GetNbinsY(); iBinY ++){
     if(iBinY < (performanceBDT_highPileUP->GetNbinsY()-1) and iBinX < (performanceBDT_highPileUP->GetNbinsX()-1) and vecPos < int(performanceValue.size()-1)){
       performanceBDT_highPileUP_float->SetBinContent(iBinX+1,iBinY+1,performanceValue.at(vecPos).value);
       if(performanceValue.at(vecPos).value-int(performanceValue.at(vecPos).value) < 0.5) 
        performanceBDT_highPileUP->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.at(vecPos).value));
       else 
        performanceBDT_highPileUP->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.at(vecPos).value)+1);
       if( performanceValue.at(vecPos).error-int(performanceValue.at(vecPos).error) <0.5)
        performanceBDT_highPileUP_error->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.at(vecPos).error));
       else
        performanceBDT_highPileUP_error->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.at(vecPos).error)+1);
       if(std::string(performanceBDT_highPileUP->GetYaxis()->GetBinLabel(iBinY+1)) == ""){ 
          performanceBDT_highPileUP->GetYaxis()->SetBinLabel(iBinY+1,(performanceValue.at(vecPos).binYName).c_str());
          performanceBDT_highPileUP_error->GetYaxis()->SetBinLabel(iBinY+1,(performanceValue.at(vecPos).binYName).c_str());
          performanceDifference->GetYaxis()->SetBinLabel(iBinY+1,(performanceValue.at(vecPos).binYName).c_str());
          performanceRatio->GetYaxis()->SetBinLabel(iBinY+1,(performanceValue.at(vecPos).binYName).c_str());
       }
       if(std::string(performanceBDT_highPileUP->GetXaxis()->GetBinLabel(iBinY+1)) == ""){
          performanceBDT_highPileUP->GetXaxis()->SetBinLabel(iBinY+1,(performanceValue.at(vecPos).binYName).c_str());
          performanceBDT_highPileUP_error->GetXaxis()->SetBinLabel(iBinY+1,(performanceValue.at(vecPos).binYName).c_str());
          performanceDifference->GetXaxis()->SetBinLabel(iBinY+1,(performanceValue.at(vecPos).binYName).c_str());
          performanceRatio->GetXaxis()->SetBinLabel(iBinY+1,(performanceValue.at(vecPos).binYName).c_str());
       }
       vecPos++;
     }    
     else if (iBinY == (performanceBDT_highPileUP->GetNbinsY()-1) and iBinX < (performanceBDT_highPileUP->GetNbinsX()-1)) continue ;
     else if (iBinY < (performanceBDT_highPileUP->GetNbinsY()-1) and iBinX == (performanceBDT_highPileUP->GetNbinsX()-1)) continue ;
     else if (iBinY == (performanceBDT_highPileUP->GetNbinsY()-1) and iBinX == (performanceBDT_highPileUP->GetNbinsX()-1)){
       performanceBDT_highPileUP_float->SetBinContent(iBinX+1,iBinY+1,performanceValue.at(vecPos).value);
       if(performanceValue.back().value-int(performanceValue.back().value) < 0.5)
        performanceBDT_highPileUP->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.back().value));
       else 
        performanceBDT_highPileUP->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.back().value)+1);
       if(performanceValue.back().error-int(performanceValue.back().error) < 0.5)
        performanceBDT_highPileUP_error->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.back().error));
       else
        performanceBDT_highPileUP_error->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.back().error)+1);

       performanceBDT_highPileUP->GetXaxis()->SetBinLabel(iBinX+1,(performanceValue.back().binXName).c_str());
       performanceBDT_highPileUP->GetYaxis()->SetBinLabel(iBinY+1,(performanceValue.back().binXName).c_str());
       performanceBDT_highPileUP_error->GetXaxis()->SetBinLabel(iBinX+1,(performanceValue.back().binXName).c_str());
       performanceBDT_highPileUP_error->GetYaxis()->SetBinLabel(iBinY+1,(performanceValue.back().binXName).c_str());
       performanceDifference->GetXaxis()->SetBinLabel(iBinX+1,(performanceValue.back().binXName).c_str());
       performanceDifference->GetYaxis()->SetBinLabel(iBinY+1,(performanceValue.back().binYName).c_str());
       performanceRatio->GetXaxis()->SetBinLabel(iBinX+1,(performanceValue.back().binXName).c_str());
       performanceRatio->GetYaxis()->SetBinLabel(iBinY+1,(performanceValue.back().binYName).c_str());
     }
    }   
   }

   ////
   TCanvas* cPerformance_highPU = new TCanvas("cPerformance_highPU","",180,52,500,550);
   cPerformance_highPU->SetGrid();
   cPerformance_highPU->SetTicks();
   cPerformance_highPU->cd();

   performanceBDT_highPileUP->SetMarkerSize(1.5);
   performanceBDT_highPileUP->SetMarkerColor(0);
   performanceBDT_highPileUP->GetXaxis()->SetLabelSize(0.035);
   performanceBDT_highPileUP->GetYaxis()->SetLabelSize(0.035);
   performanceBDT_highPileUP->SetLabelOffset(0.011);// label offset on x axis                                                                                              
   performanceBDT_highPileUP->SetMaximum(performanceBDT_highPileUP->GetMaximum());
   performanceBDT_highPileUP->SetMinimum(performanceBDT_highPileUP->GetMinimum());
   performanceBDT_highPileUP->Draw("colz");
   performanceBDT_highPileUP->Draw("textsame");
   performanceBDT_highPileUP->GetXaxis()->LabelsOption("v");
   performanceBDT_highPileUP->GetYaxis()->LabelsOption("h");

   tex->Draw();
   tex2->Draw();
   tex3->Draw();

   cPerformance_highPU->Print((outputDirectory+"/PerformanceBDT_highPU.pdf").c_str(),"pdf");
   cPerformance_highPU->Print((outputDirectory+"/PerformanceBDT_highPU.png").c_str(),"png");
   cPerformance_highPU->Print((outputDirectory+"/PerformanceBDT_highPU.root").c_str(),"root");
   cPerformance_highPU->SetLogz();
   cPerformance_highPU->Print((outputDirectory+"/PerformanceBDT_highPU_Log.pdf").c_str(),"pdf");
   cPerformance_highPU->Print((outputDirectory+"/PerformanceBDT_highPU_Log.png").c_str(),"png");
   cPerformance_highPU->Print((outputDirectory+"/PerformanceBDT_highPU_Log.root").c_str(),"root");

   ////
   TCanvas* cPerformance_highPU_error = new TCanvas("cPerformance_highPU_error","",180,52,500,550);

   cPerformance_highPU_error->SetGrid();
   cPerformance_highPU_error->SetTicks();
   cPerformance_highPU_error->cd();

   performanceBDT_highPileUP_error->SetMarkerSize(1.5);
   performanceBDT_highPileUP_error->SetMarkerColor(0);
   performanceBDT_highPileUP_error->GetXaxis()->SetLabelSize(0.035);
   performanceBDT_highPileUP_error->GetYaxis()->SetLabelSize(0.035);
   performanceBDT_highPileUP_error->SetLabelOffset(0.011);// label offset on x axis                                                                                  
   performanceBDT_highPileUP_error->SetMaximum(performanceBDT_highPileUP_error->GetMaximum());
   performanceBDT_highPileUP_error->SetMinimum(performanceBDT_highPileUP_error->GetMinimum());
   performanceBDT_highPileUP_error->Draw("colz");
   performanceBDT_highPileUP_error->Draw("textsame");
   performanceBDT_highPileUP_error->GetXaxis()->LabelsOption("v");
   performanceBDT_highPileUP_error->GetYaxis()->LabelsOption("h");

   tex->Draw();
   tex2->Draw();
   tex3->Draw();

   cPerformance_highPU_error->Print((outputDirectory+"/PerformanceBDT_highPU_error.pdf").c_str(),"pdf");
   cPerformance_highPU_error->Print((outputDirectory+"/PerformanceBDT_highPU_error.png").c_str(),"png");
   cPerformance_highPU_error->Print((outputDirectory+"/PerformanceBDT_highPU_error.root").c_str(),"root");
   cPerformance_highPU_error->SetLogz();
   cPerformance_highPU_error->Print((outputDirectory+"/PerformanceBDT_highPU_Log_error.pdf").c_str(),"pdf");
   cPerformance_highPU_error->Print((outputDirectory+"/PerformanceBDT_highPU_Log_error.png").c_str(),"png");
   cPerformance_highPU_error->Print((outputDirectory+"/PerformanceBDT_highPU_Log_error.root").c_str(),"root");
  
   // Difference Plot
   TCanvas* cPerformance_difference = new TCanvas("cPerformance_difference","",180,52,500,550);

   cPerformance_difference->SetGrid();
   cPerformance_difference->SetTicks();
   cPerformance_difference->cd();
 
   for( int binX = 0 ; binX < performanceBDT_highPileUP->GetNbinsX(); binX++){
    for( int binY = binX ; binY < performanceBDT_highPileUP->GetNbinsY(); binY++){
     if(performanceBDT_highPileUP->GetBinContent(binX+1,binY+1)-performanceBDT_lowPileUP->GetBinContent(binX+1,binY+1) == 0) continue;
     if(performanceBDT_highPileUP->GetBinContent(binX+1,binY+1) == 0 or performanceBDT_lowPileUP->GetBinContent(binX+1,binY+1) == 0) continue;
     else performanceDifference->SetBinContent(binX+1,binY+1,performanceBDT_highPileUP->GetBinContent(binX+1,binY+1)-performanceBDT_lowPileUP->GetBinContent(binX+1,binY+1));
    }
   }

   performanceDifference->SetMarkerSize(1.5);
   performanceDifference->SetMarkerColor(0);
   performanceDifference->GetXaxis()->SetLabelSize(0.035);
   performanceDifference->GetYaxis()->SetLabelSize(0.035);
   performanceDifference->SetLabelOffset(0.011);// label offset on x axis                                                                                                                     
   performanceDifference->SetMaximum(performanceDifference->GetMaximum());
   performanceDifference->SetMinimum(performanceDifference->GetMinimum());

   performanceDifference->Draw("colz");
   performanceDifference->Draw("textsame");
   performanceDifference->GetXaxis()->LabelsOption("v");
   performanceDifference->GetYaxis()->LabelsOption("h");

   tex->Draw();
   tex2->Draw();
   tex3->Draw();

   cPerformance_difference->Print((outputDirectory+"/PerformanceBDT_Difference.pdf").c_str(),"pdf");
   cPerformance_difference->Print((outputDirectory+"/PerformanceBDT_Difference.png").c_str(),"png");
   cPerformance_difference->Print((outputDirectory+"/PerformanceBDT_Difference.root").c_str(),"root");
   cPerformance_difference->SetLogz();
   cPerformance_difference->Print((outputDirectory+"/PerformanceBDT_Difference_Log.pdf").c_str(),"pdf");
   cPerformance_difference->Print((outputDirectory+"/PerformanceBDT_Difference_Log.png").c_str(),"png");
   cPerformance_difference->Print((outputDirectory+"/PerformanceBDT_Difference_Log.root").c_str(),"root");

   // Ratio Plot
   TCanvas* cPerformance_ratio = new TCanvas("cPerformance_ratio","",180,52,500,550);

   cPerformance_ratio->SetGrid();
   cPerformance_ratio->SetTicks();
   cPerformance_ratio->cd();
   int minimum = 100;

   for( int binX = 0 ; binX < performanceBDT_highPileUP->GetNbinsX(); binX++){
    for( int binY = binX ; binY < performanceBDT_highPileUP->GetNbinsY(); binY++){
     if(performanceBDT_highPileUP_float->GetBinContent(binX+1,binY+1)/performanceBDT_lowPileUP_float->GetBinContent(binX+1,binY+1) == 0) continue;
     if(performanceBDT_highPileUP_float->GetBinContent(binX+1,binY+1) == 0 or performanceBDT_lowPileUP_float->GetBinContent(binX+1,binY+1) == 0) continue;     else performanceRatio->SetBinContent(binX+1,binY+1,(performanceBDT_highPileUP_float->GetBinContent(binX+1,binY+1)/performanceBDT_lowPileUP_float->GetBinContent(binX+1,binY+1))*100);    
     if(performanceRatio->GetBinContent(binX+1,binY+1)-int(performanceRatio->GetBinContent(binX+1,binY+1)) < 0.5) performanceRatio->SetBinContent(binX+1,binY+1,int(performanceRatio->GetBinContent(binX+1,binY+1)));
     else performanceRatio->SetBinContent(binX+1,binY+1,int(performanceRatio->GetBinContent(binX+1,binY+1))+1);
     if (performanceRatio->GetBinContent(binX+1,binY+1) > 0 and performanceRatio->GetBinContent(binX+1,binY+1) < minimum ) minimum = performanceRatio->GetBinContent(binX+1,binY+1);
    }
   }

   performanceRatio->SetMarkerSize(1.5);
   performanceRatio->SetMarkerColor(0);
   performanceRatio->GetXaxis()->SetLabelSize(0.035);
   performanceRatio->GetYaxis()->SetLabelSize(0.035);
   performanceRatio->SetLabelOffset(0.011);// label offset on x axis                                                                                                                     

   performanceRatio->SetMinimum(minimum);
   performanceRatio->Draw("colz");
   performanceRatio->Draw("textsame");
   performanceRatio->GetXaxis()->LabelsOption("v");
   performanceRatio->GetYaxis()->LabelsOption("h");

   tex->Draw();
   tex2->Draw();
   tex3->Draw();

  cPerformance_ratio->Print((outputDirectory+"/PerformanceBDT_Ratio.pdf").c_str(),"pdf");
  cPerformance_ratio->Print((outputDirectory+"/PerformanceBDT_Ratio.png").c_str(),"png");
  cPerformance_ratio->Print((outputDirectory+"/PerformanceBDT_Ratio.root").c_str(),"root");
  cPerformance_ratio->SetLogz();
  cPerformance_ratio->Print((outputDirectory+"/PerformanceBDT_Ratio_Log.pdf").c_str(),"pdf");
  cPerformance_ratio->Print((outputDirectory+"/PerformanceBDT_Ratio_Log.png").c_str(),"png");
  cPerformance_ratio->Print((outputDirectory+"/PerformanceBDT_Ratio_Log.root").c_str(),"root");

  }

  else{

   std::vector<std::string> InputTreeSignal ;
   if(Options.existsAs<std::vector<std::string>>("InputTreeSignal"))
    InputTreeSignal = Options.getParameter<std::vector<std::string>>("InputTreeSignal");
   else{ std::cout<<" Exit from code, no input set found for low pile-up files"<<std::endl; return -1; }

   std::vector<std::string> InputTreeBackground ;
   if(Options.existsAs<std::vector<std::string>>("InputTreeBackground"))
    InputTreeBackground = Options.getParameter<std::vector<std::string>>("InputTreeBackground");
   else{ std::cout<<" Exit from code, no input set found for low pile-up files"<<std::endl; return -1; }

   std::vector<edm::ParameterSet> InputVariableName ;
   if(Options.existsAs<std::vector<edm::ParameterSet>>("InputVariableName"))
     InputVariableName = Options.getParameter<std::vector<edm::ParameterSet>>("InputVariableName");
   else{ std::cout<<" Exit from code, no input set of variables found "<<std::endl; return -1; }

   std::string outputDirectory;
   if(Options.existsAs<std::string>("outputDirectory"))
    outputDirectory = Options.getParameter<std::string>("outputDirectory");
   else{ outputDirectory = "output/"; }

   double signalEfficiencyTarget;
   if(Options.existsAs<double>("signalEfficiencyTarget"))
    signalEfficiencyTarget = Options.getParameter<double>("signalEfficiencyTarget");
   else{ signalEfficiencyTarget = 0.5; }

   std::string TreeName;
   if(Options.existsAs<std::string>("TreeName"))
     TreeName = Options.getParameter<std::string>("TreeName");
   else{ TreeName = "chs"; }

   std::string cutSignalLowPileUp;
   if(Options.existsAs<std::string>("cutSignalLowPileUp"))
     cutSignalLowPileUp = Options.getParameter<std::string>("cutSignalLowPileUp");
   else{ cutSignalLowPileUp = "npu < 40"; }

   std::string cutBackgroundLowPileUp;
   if(Options.existsAs<std::string>("cutBackgroundLowPileUp"))
     cutBackgroundLowPileUp = Options.getParameter<std::string>("cutBackgroundLowPileUp");
   else{ cutBackgroundLowPileUp = "npu < 40"; }
 
   system(("mkdir -p "+outputDirectory).c_str());

   std::vector<std::string> branchName ;
   std::vector<std::string> reduceNameX ;
   std::vector<std::string> reduceNameY ;
   
   for(size_t iVar = 0; iVar < InputVariableName.size(); iVar++){
     branchName.push_back(InputVariableName.at(iVar).getParameter<std::string>("branchName"));
     reduceNameX.push_back(InputVariableName.at(iVar).getParameter<std::string>("reduceNameX"));
     reduceNameY.push_back(InputVariableName.at(iVar).getParameter<std::string>("reduceNameY"));
   }

   std::vector<TTree*> signalTree ;
   std::vector<TTree*> backgroundTree ;
 
   TFile* inputFile = NULL ;
   for(size_t iFile = 0; iFile < InputTreeSignal.size(); iFile++){
     inputFile = TFile::Open(InputTreeSignal.at(iFile).c_str(),"READ");
     if(inputFile == 0 or inputFile == NULL) continue;
     signalTree.push_back((TTree*) inputFile->Get(TreeName.c_str())); // make the tree list                                                                                              
   }
   for(size_t iFile = 0; iFile < InputTreeBackground.size(); iFile++){
     inputFile = TFile::Open(InputTreeBackground.at(iFile).c_str(),"READ");
     if(inputFile == 0 or inputFile == NULL) continue;
     backgroundTree.push_back((TTree*) inputFile->Get(TreeName.c_str())); // make the tree list                                                                                              
   }


   // start the performance evaluation
   std::vector<MatrixEntry> performanceValue = rocMatrix(signalTree,backgroundTree,branchName,reduceNameX,reduceNameY,cutSignalLowPileUp,cutBackgroundLowPileUp,signalEfficiencyTarget);

   //For getting the number of bins is enough to cycle on all the entry and count the diagonal terms
   int numberOfBins = 0;
   for(unsigned int ientry = 0 ; ientry < performanceValue.size(); ientry++){
     if(performanceValue.at(ientry).binXName == performanceValue.at(ientry).binYName) numberOfBins++;
   }

   TH2F* performanceBDT_lowPileUP        = new TH2F("backgroundMatrix_lowPileUP",      "",numberOfBins,0,numberOfBins,numberOfBins,0,numberOfBins);
   TH2F* performanceBDT_lowPileUP_float  = new TH2F("backgroundMatrix_lowPileUP_float","",numberOfBins,0,numberOfBins,numberOfBins,0,numberOfBins);
   TH2F* performanceBDT_lowPileUP_error  = new TH2F("backgroundMatrix_lowPileUP_error","",numberOfBins,0,numberOfBins,numberOfBins,0,numberOfBins);

   for(int iBinX = 0; iBinX < performanceBDT_lowPileUP->GetNbinsX(); iBinX ++){
    for(int iBinY = 0; iBinY < performanceBDT_lowPileUP->GetNbinsY(); iBinY ++){
     performanceBDT_lowPileUP->SetBinContent(iBinX+1,iBinY+1,0.);
     performanceBDT_lowPileUP_float->SetBinContent(iBinX+1,iBinY+1,0.);
     performanceBDT_lowPileUP_error->SetBinContent(iBinX+1,iBinY+1,0.);
    }
   }

   int vecPos = 0;
   for(int iBinX = 0; iBinX < performanceBDT_lowPileUP->GetNbinsX(); iBinX ++){  
    for(int iBinY = iBinX; iBinY < performanceBDT_lowPileUP->GetNbinsY(); iBinY ++){
     if(iBinY < (performanceBDT_lowPileUP->GetNbinsY()-1) and iBinX < (performanceBDT_lowPileUP->GetNbinsX()-1) and vecPos < int(performanceValue.size()-1)){

       performanceBDT_lowPileUP_float->SetBinContent(iBinX+1,iBinY+1,performanceValue.at(vecPos).value);
       if(performanceValue.at(vecPos).value - int(performanceValue.at(vecPos).value) < 0.5)
        performanceBDT_lowPileUP->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.at(vecPos).value));
       else
        performanceBDT_lowPileUP->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.at(vecPos).value)+1);
       if(performanceValue.at(vecPos).error - int(performanceValue.at(vecPos).error) < 0.5)
        performanceBDT_lowPileUP_error->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.at(vecPos).error));
       else
        performanceBDT_lowPileUP_error->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.at(vecPos).error)+1);

        performanceBDT_lowPileUP->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.at(vecPos).value));
        performanceBDT_lowPileUP_error->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.at(vecPos).error));

       if(std::string(performanceBDT_lowPileUP->GetYaxis()->GetBinLabel(iBinY+1)) == ""){
	 performanceBDT_lowPileUP->GetYaxis()->SetBinLabel(iBinY+1,(performanceValue.at(vecPos).binYName).c_str());
         performanceBDT_lowPileUP_error->GetYaxis()->SetBinLabel(iBinY+1,(performanceValue.at(vecPos).binYName).c_str());
       }
       if(std::string(performanceBDT_lowPileUP->GetXaxis()->GetBinLabel(iBinY+1)) == ""){
         performanceBDT_lowPileUP->GetXaxis()->SetBinLabel(iBinY+1,(performanceValue.at(vecPos).binYName).c_str());
         performanceBDT_lowPileUP_error->GetXaxis()->SetBinLabel(iBinY+1,(performanceValue.at(vecPos).binYName).c_str());
       }
       vecPos++; 
     }    
     else if (iBinY == (performanceBDT_lowPileUP->GetNbinsY()-1) and iBinX < (performanceBDT_lowPileUP->GetNbinsX()-1)) continue ;
     else if (iBinY < (performanceBDT_lowPileUP->GetNbinsY()-1) and iBinX == (performanceBDT_lowPileUP->GetNbinsX()-1)) continue ;
     else if (iBinY == (performanceBDT_lowPileUP->GetNbinsY()-1) and iBinX == (performanceBDT_lowPileUP->GetNbinsX()-1)){
       performanceBDT_lowPileUP_float->SetBinContent(iBinX+1,iBinY+1,performanceValue.at(vecPos).value);
       if(performanceValue.back().value-int(performanceValue.back().value)<0.5)
        performanceBDT_lowPileUP->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.back().value));
       else
        performanceBDT_lowPileUP->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.back().value)+1);

       performanceBDT_lowPileUP->GetXaxis()->SetBinLabel(iBinX+1,(performanceValue.back().binXName).c_str());
       performanceBDT_lowPileUP->GetYaxis()->SetBinLabel(iBinY+1,(performanceValue.back().binXName).c_str());
       if(performanceValue.back().error-int(performanceValue.back().error)<0.5)
        performanceBDT_lowPileUP_error->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.back().error));
       else  
        performanceBDT_lowPileUP_error->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.back().error)+1);
       performanceBDT_lowPileUP_error->GetXaxis()->SetBinLabel(iBinX+1,(performanceValue.back().binXName).c_str());
       performanceBDT_lowPileUP_error->GetYaxis()->SetBinLabel(iBinY+1,(performanceValue.back().binXName).c_str());
     }
    }
   }

   TCanvas* cPerformance_lowPU = new TCanvas("cPerformance_lowPU","",180,52,500,550);

   cPerformance_lowPU->SetGrid();
   cPerformance_lowPU->SetTicks();
   cPerformance_lowPU->cd();

   performanceBDT_lowPileUP->SetMarkerSize(1.5);
   performanceBDT_lowPileUP->SetMarkerColor(0);
   performanceBDT_lowPileUP->GetXaxis()->SetLabelSize(0.035);
   performanceBDT_lowPileUP->GetYaxis()->SetLabelSize(0.035);
   performanceBDT_lowPileUP->SetLabelOffset(0.011);// label offset on x axis                                                                                                  
   performanceBDT_lowPileUP->SetMaximum(performanceBDT_lowPileUP->GetMaximum());
   performanceBDT_lowPileUP->SetMinimum(performanceBDT_lowPileUP->GetMinimum());
   performanceBDT_lowPileUP->Draw("colz");
   performanceBDT_lowPileUP->Draw("textsame");
   performanceBDT_lowPileUP->GetXaxis()->LabelsOption("v");
   performanceBDT_lowPileUP->GetYaxis()->LabelsOption("h");


   TLatex* tex = new TLatex(0.85,0.92," 13 TeV");
   tex->SetNDC();
   tex->SetTextAlign(31);
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
   TLatex* tex2 = new TLatex(0.13,0.92,"CMS");
   tex2->SetNDC();
   tex2->SetTextFont(61);
   tex2->SetTextSize(0.04);
   tex2->SetLineWidth(2);
   tex2->Draw();
   TLatex* tex3 = new TLatex(0.2324,0.92,"Preliminary Simulation");
   tex3->SetNDC();
   tex3->SetTextFont(52);
   tex3->SetTextSize(0.0304);
   tex3->SetLineWidth(2);
   tex3->Draw();

   cPerformance_lowPU->Print((outputDirectory+"/PerformanceBDT_lowPU.pdf").c_str(),"pdf");
   cPerformance_lowPU->Print((outputDirectory+"/PerformanceBDT_lowPU.png").c_str(),"png");
   cPerformance_lowPU->Print((outputDirectory+"/PerformanceBDT_lowPU.root").c_str(),"root");
   cPerformance_lowPU->SetLogz();
   cPerformance_lowPU->Print((outputDirectory+"/PerformanceBDT_lowPU_Log.pdf").c_str(),"pdf");
   cPerformance_lowPU->Print((outputDirectory+"/PerformanceBDT_lowPU_Log.png").c_str(),"png");
   cPerformance_lowPU->Print((outputDirectory+"/PerformanceBDT_lowPU_Log.root").c_str(),"root");

   // Error
   TCanvas* cPerformance_lowPU_error = new TCanvas("cPerformance_lowPU_error","",180,52,500,550);

   cPerformance_lowPU_error->SetGrid();
   cPerformance_lowPU_error->SetTicks();
   cPerformance_lowPU_error->cd();

   performanceBDT_lowPileUP_error->SetMarkerSize(1.5);
   performanceBDT_lowPileUP_error->SetMarkerColor(0);
   performanceBDT_lowPileUP_error->GetXaxis()->SetLabelSize(0.035);
   performanceBDT_lowPileUP_error->GetYaxis()->SetLabelSize(0.035);
   performanceBDT_lowPileUP_error->SetLabelOffset(0.011);// label offset on x axis                                                                                              
   performanceBDT_lowPileUP_error->SetMaximum(performanceBDT_lowPileUP_error->GetMaximum());
   performanceBDT_lowPileUP_error->SetMinimum(performanceBDT_lowPileUP_error->GetMinimum());
   performanceBDT_lowPileUP_error->Draw("colz");
   performanceBDT_lowPileUP_error->Draw("textsame");
   performanceBDT_lowPileUP_error->GetXaxis()->LabelsOption("v");
   performanceBDT_lowPileUP_error->GetYaxis()->LabelsOption("h");

   tex->Draw();
   tex2->Draw();
   tex3->Draw();

   cPerformance_lowPU_error->Print((outputDirectory+"/PerformanceBDT_lowPU_error.pdf").c_str(),"pdf");
   cPerformance_lowPU_error->Print((outputDirectory+"/PerformanceBDT_lowPU_error.png").c_str(),"png");
   cPerformance_lowPU_error->Print((outputDirectory+"/PerformanceBDT_lowPU_error.root").c_str(),"root");
   cPerformance_lowPU_error->SetLogz();
   cPerformance_lowPU_error->Print((outputDirectory+"/PerformanceBDT_lowPU_Log_error.pdf").c_str(),"pdf");
   cPerformance_lowPU_error->Print((outputDirectory+"/PerformanceBDT_lowPU_Log_error.png").c_str(),"png");
   cPerformance_lowPU_error->Print((outputDirectory+"/PerformanceBDT_lowPU_Log_error.root").c_str(),"root");

   // high PU part
   std::string cutSignalHighPileUp;
   if(Options.existsAs<std::string>("cutSignalHighPileUp"))
     cutSignalHighPileUp = Options.getParameter<std::string>("cutSignalHighPileUp");
   else{ cutSignalHighPileUp = "npu >= 40"; }

   std::string cutBackgroundHighPileUp;
   if(Options.existsAs<std::string>("cutBackgroundHighPileUp"))
     cutBackgroundHighPileUp = Options.getParameter<std::string>("cutBackgroundHighPileUp");
   else{ cutBackgroundHighPileUp = "npu >= 40"; }

   // start the performance evaluation
   performanceValue = rocMatrix(signalTree,backgroundTree,branchName,reduceNameX,reduceNameY,cutSignalHighPileUp,cutBackgroundHighPileUp,signalEfficiencyTarget);

   //for getting the number of bins is enough to cycle on all the entry and count the diagonal terms
   numberOfBins = 0;
   for(unsigned int ientry = 0 ; ientry < performanceValue.size(); ientry++){
    if(performanceValue.at(ientry).binXName == performanceValue.at(ientry).binYName) numberOfBins++;
   }

   TH2F* performanceBDT_highPileUP        = new TH2F("backgroundMatrix_highPileUP","",numberOfBins,0,numberOfBins,numberOfBins,0,numberOfBins);
   TH2F* performanceBDT_highPileUP_float  = new TH2F("backgroundMatrix_highPileUP_float","",numberOfBins,0,numberOfBins,numberOfBins,0,numberOfBins);
   TH2F* performanceBDT_highPileUP_error  = new TH2F("backgroundMatrix_highPileUP_error","",numberOfBins,0,numberOfBins,numberOfBins,0,numberOfBins);
   TH2F* performanceDifference            = new TH2F("performanceDifference","",numberOfBins,0,numberOfBins,numberOfBins,0,numberOfBins);
   TH2F* performanceRatio                 = new TH2F("performanceRatio","",numberOfBins,0,numberOfBins,numberOfBins,0,numberOfBins);

   for(int iBinX = 0; iBinX < performanceBDT_highPileUP->GetNbinsX(); iBinX ++){
    for(int iBinY = 0; iBinY < performanceBDT_highPileUP->GetNbinsY(); iBinY ++){
     performanceBDT_highPileUP->SetBinContent(iBinX+1,iBinY+1,0.);
     performanceBDT_highPileUP_error->SetBinContent(iBinX+1,iBinY+1,0.);
    }
   }

   vecPos = 0;
   for(int iBinX = 0; iBinX < performanceBDT_highPileUP->GetNbinsX(); iBinX ++){
    for(int iBinY = iBinX; iBinY < performanceBDT_highPileUP->GetNbinsY(); iBinY ++){
     if(iBinY < (performanceBDT_highPileUP->GetNbinsY()-1) and iBinX < (performanceBDT_highPileUP->GetNbinsX()-1) and vecPos < int(performanceValue.size()-1)){
       performanceBDT_highPileUP_float->SetBinContent(iBinX+1,iBinY+1,performanceValue.at(vecPos).value);
       if(performanceValue.at(vecPos).value-int(performanceValue.at(vecPos).value) < 0.5) 
        performanceBDT_highPileUP->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.at(vecPos).value));
       else 
        performanceBDT_highPileUP->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.at(vecPos).value)+1);
       if( performanceValue.at(vecPos).error-int(performanceValue.at(vecPos).error) <0.5)
        performanceBDT_highPileUP_error->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.at(vecPos).error));
       else
        performanceBDT_highPileUP_error->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.at(vecPos).error)+1);
       if(std::string(performanceBDT_highPileUP->GetYaxis()->GetBinLabel(iBinY+1)) == ""){ 
          performanceBDT_highPileUP->GetYaxis()->SetBinLabel(iBinY+1,(performanceValue.at(vecPos).binYName).c_str());
          performanceBDT_highPileUP_error->GetYaxis()->SetBinLabel(iBinY+1,(performanceValue.at(vecPos).binYName).c_str());
          performanceDifference->GetYaxis()->SetBinLabel(iBinY+1,(performanceValue.at(vecPos).binYName).c_str());
          performanceRatio->GetYaxis()->SetBinLabel(iBinY+1,(performanceValue.at(vecPos).binYName).c_str());
       }
       if(std::string(performanceBDT_highPileUP->GetXaxis()->GetBinLabel(iBinY+1)) == ""){
          performanceBDT_highPileUP->GetXaxis()->SetBinLabel(iBinY+1,(performanceValue.at(vecPos).binYName).c_str());
          performanceBDT_highPileUP_error->GetXaxis()->SetBinLabel(iBinY+1,(performanceValue.at(vecPos).binYName).c_str());
          performanceDifference->GetXaxis()->SetBinLabel(iBinY+1,(performanceValue.at(vecPos).binYName).c_str());
          performanceRatio->GetXaxis()->SetBinLabel(iBinY+1,(performanceValue.at(vecPos).binYName).c_str());
       }
       vecPos++;
     }    
     else if (iBinY == (performanceBDT_highPileUP->GetNbinsY()-1) and iBinX < (performanceBDT_highPileUP->GetNbinsX()-1)) continue ;
     else if (iBinY < (performanceBDT_highPileUP->GetNbinsY()-1) and iBinX == (performanceBDT_highPileUP->GetNbinsX()-1)) continue ;
     else if (iBinY == (performanceBDT_highPileUP->GetNbinsY()-1) and iBinX == (performanceBDT_highPileUP->GetNbinsX()-1)){
       performanceBDT_highPileUP_float->SetBinContent(iBinX+1,iBinY+1,performanceValue.at(vecPos).value);
       if(performanceValue.back().value-int(performanceValue.back().value) < 0.5)
        performanceBDT_highPileUP->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.back().value));
       else 
        performanceBDT_highPileUP->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.back().value)+1);
       if(performanceValue.back().error-int(performanceValue.back().error) < 0.5)
        performanceBDT_highPileUP_error->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.back().error));
       else
        performanceBDT_highPileUP_error->SetBinContent(iBinX+1,iBinY+1,int(performanceValue.back().error)+1);

       performanceBDT_highPileUP->GetXaxis()->SetBinLabel(iBinX+1,(performanceValue.back().binXName).c_str());
       performanceBDT_highPileUP->GetYaxis()->SetBinLabel(iBinY+1,(performanceValue.back().binXName).c_str());
       performanceBDT_highPileUP_error->GetXaxis()->SetBinLabel(iBinX+1,(performanceValue.back().binXName).c_str());
       performanceBDT_highPileUP_error->GetYaxis()->SetBinLabel(iBinY+1,(performanceValue.back().binXName).c_str());
       performanceDifference->GetXaxis()->SetBinLabel(iBinX+1,(performanceValue.back().binXName).c_str());
       performanceDifference->GetYaxis()->SetBinLabel(iBinY+1,(performanceValue.back().binYName).c_str());
       performanceRatio->GetXaxis()->SetBinLabel(iBinX+1,(performanceValue.back().binXName).c_str());
       performanceRatio->GetYaxis()->SetBinLabel(iBinY+1,(performanceValue.back().binYName).c_str());
     }
    }   
   }

   ////
   TCanvas* cPerformance_highPU = new TCanvas("cPerformance_highPU","",180,52,500,550);
   cPerformance_highPU->SetGrid();
   cPerformance_highPU->SetTicks();
   cPerformance_highPU->cd();

   performanceBDT_highPileUP->SetMarkerSize(1.5);
   performanceBDT_highPileUP->SetMarkerColor(0);
   performanceBDT_highPileUP->GetXaxis()->SetLabelSize(0.035);
   performanceBDT_highPileUP->GetYaxis()->SetLabelSize(0.035);
   performanceBDT_highPileUP->SetLabelOffset(0.011);// label offset on x axis                                                                                              
   performanceBDT_highPileUP->SetMaximum(performanceBDT_highPileUP->GetMaximum());
   performanceBDT_highPileUP->SetMinimum(performanceBDT_highPileUP->GetMinimum());
   performanceBDT_highPileUP->Draw("colz");
   performanceBDT_highPileUP->Draw("textsame");
   performanceBDT_highPileUP->GetXaxis()->LabelsOption("v");
   performanceBDT_highPileUP->GetYaxis()->LabelsOption("h");

   tex->Draw();
   tex2->Draw();
   tex3->Draw();

   cPerformance_highPU->Print((outputDirectory+"/PerformanceBDT_highPU.pdf").c_str(),"pdf");
   cPerformance_highPU->Print((outputDirectory+"/PerformanceBDT_highPU.png").c_str(),"png");
   cPerformance_highPU->Print((outputDirectory+"/PerformanceBDT_highPU.root").c_str(),"root");
   cPerformance_highPU->SetLogz();
   cPerformance_highPU->Print((outputDirectory+"/PerformanceBDT_highPU_Log.pdf").c_str(),"pdf");
   cPerformance_highPU->Print((outputDirectory+"/PerformanceBDT_highPU_Log.png").c_str(),"png");
   cPerformance_highPU->Print((outputDirectory+"/PerformanceBDT_highPU_Log.root").c_str(),"root");

   ////
   TCanvas* cPerformance_highPU_error = new TCanvas("cPerformance_highPU_error","",180,52,500,550);

   cPerformance_highPU_error->SetGrid();
   cPerformance_highPU_error->SetTicks();
   cPerformance_highPU_error->cd();

   performanceBDT_highPileUP_error->SetMarkerSize(1.5);
   performanceBDT_highPileUP_error->SetMarkerColor(0);
   performanceBDT_highPileUP_error->GetXaxis()->SetLabelSize(0.035);
   performanceBDT_highPileUP_error->GetYaxis()->SetLabelSize(0.035);
   performanceBDT_highPileUP_error->SetLabelOffset(0.011);// label offset on x axis                                                                                  
   performanceBDT_highPileUP_error->SetMaximum(performanceBDT_highPileUP_error->GetMaximum());
   performanceBDT_highPileUP_error->SetMinimum(performanceBDT_highPileUP_error->GetMinimum());
   performanceBDT_highPileUP_error->Draw("colz");
   performanceBDT_highPileUP_error->Draw("textsame");
   performanceBDT_highPileUP_error->GetXaxis()->LabelsOption("v");
   performanceBDT_highPileUP_error->GetYaxis()->LabelsOption("h");

   tex->Draw();
   tex2->Draw();
   tex3->Draw();

   cPerformance_highPU_error->Print((outputDirectory+"/PerformanceBDT_highPU_error.pdf").c_str(),"pdf");
   cPerformance_highPU_error->Print((outputDirectory+"/PerformanceBDT_highPU_error.png").c_str(),"png");
   cPerformance_highPU_error->Print((outputDirectory+"/PerformanceBDT_highPU_error.root").c_str(),"root");
   cPerformance_highPU_error->SetLogz();
   cPerformance_highPU_error->Print((outputDirectory+"/PerformanceBDT_highPU_Log_error.pdf").c_str(),"pdf");
   cPerformance_highPU_error->Print((outputDirectory+"/PerformanceBDT_highPU_Log_error.png").c_str(),"png");
   cPerformance_highPU_error->Print((outputDirectory+"/PerformanceBDT_highPU_Log_error.root").c_str(),"root");
  
   // Difference Plot
   TCanvas* cPerformance_difference = new TCanvas("cPerformance_difference","",180,52,500,550);

   cPerformance_difference->SetGrid();
   cPerformance_difference->SetTicks();
   cPerformance_difference->cd();
 
   for( int binX = 0 ; binX < performanceBDT_highPileUP->GetNbinsX(); binX++){
    for( int binY = binX ; binY < performanceBDT_highPileUP->GetNbinsY(); binY++){
     if(performanceBDT_highPileUP->GetBinContent(binX+1,binY+1)-performanceBDT_lowPileUP->GetBinContent(binX+1,binY+1) == 0) continue;
     if(performanceBDT_highPileUP->GetBinContent(binX+1,binY+1) == 0 or performanceBDT_lowPileUP->GetBinContent(binX+1,binY+1) == 0) continue;
     else performanceDifference->SetBinContent(binX+1,binY+1,performanceBDT_highPileUP->GetBinContent(binX+1,binY+1)-performanceBDT_lowPileUP->GetBinContent(binX+1,binY+1));
    }
   }

   performanceDifference->SetMarkerSize(1.5);
   performanceDifference->SetMarkerColor(0);
   performanceDifference->GetXaxis()->SetLabelSize(0.035);
   performanceDifference->GetYaxis()->SetLabelSize(0.035);
   performanceDifference->SetLabelOffset(0.011);// label offset on x axis                                                                                                                     
   performanceDifference->SetMaximum(performanceDifference->GetMaximum());
   performanceDifference->SetMinimum(performanceDifference->GetMinimum());

   performanceDifference->Draw("colz");
   performanceDifference->Draw("textsame");
   performanceDifference->GetXaxis()->LabelsOption("v");
   performanceDifference->GetYaxis()->LabelsOption("h");

   tex->Draw();
   tex2->Draw();
   tex3->Draw();

   cPerformance_difference->Print((outputDirectory+"/PerformanceBDT_Difference.pdf").c_str(),"pdf");
   cPerformance_difference->Print((outputDirectory+"/PerformanceBDT_Difference.png").c_str(),"png");
   cPerformance_difference->Print((outputDirectory+"/PerformanceBDT_Difference.root").c_str(),"root");
   cPerformance_difference->SetLogz();
   cPerformance_difference->Print((outputDirectory+"/PerformanceBDT_Difference_Log.pdf").c_str(),"pdf");
   cPerformance_difference->Print((outputDirectory+"/PerformanceBDT_Difference_Log.png").c_str(),"png");
   cPerformance_difference->Print((outputDirectory+"/PerformanceBDT_Difference_Log.root").c_str(),"root");

   // Ratio Plot
   TCanvas* cPerformance_ratio = new TCanvas("cPerformance_ratio","",180,52,500,550);

   cPerformance_ratio->SetGrid();
   cPerformance_ratio->SetTicks();
   cPerformance_ratio->cd();
   int minimum = 100;

   for( int binX = 0 ; binX < performanceBDT_highPileUP->GetNbinsX(); binX++){
    for( int binY = binX ; binY < performanceBDT_highPileUP->GetNbinsY(); binY++){
     if(performanceBDT_highPileUP_float->GetBinContent(binX+1,binY+1)/performanceBDT_lowPileUP_float->GetBinContent(binX+1,binY+1) == 0) continue;
     if(performanceBDT_highPileUP_float->GetBinContent(binX+1,binY+1) == 0 or performanceBDT_lowPileUP_float->GetBinContent(binX+1,binY+1) == 0) continue;     else performanceRatio->SetBinContent(binX+1,binY+1,(performanceBDT_highPileUP_float->GetBinContent(binX+1,binY+1)/performanceBDT_lowPileUP_float->GetBinContent(binX+1,binY+1))*100);    
     if(performanceRatio->GetBinContent(binX+1,binY+1)-int(performanceRatio->GetBinContent(binX+1,binY+1)) < 0.5) performanceRatio->SetBinContent(binX+1,binY+1,int(performanceRatio->GetBinContent(binX+1,binY+1)));
     else performanceRatio->SetBinContent(binX+1,binY+1,int(performanceRatio->GetBinContent(binX+1,binY+1))+1);
     if (performanceRatio->GetBinContent(binX+1,binY+1) > 0 and performanceRatio->GetBinContent(binX+1,binY+1) < minimum ) minimum = performanceRatio->GetBinContent(binX+1,binY+1);
    }
   }

   performanceRatio->SetMarkerSize(1.5);
   performanceRatio->SetMarkerColor(0);
   performanceRatio->GetXaxis()->SetLabelSize(0.035);
   performanceRatio->GetYaxis()->SetLabelSize(0.035);
   performanceRatio->SetLabelOffset(0.011);// label offset on x axis                                                                                                                     

   performanceRatio->SetMinimum(minimum);
   performanceRatio->Draw("colz");
   performanceRatio->Draw("textsame");
   performanceRatio->GetXaxis()->LabelsOption("v");
   performanceRatio->GetYaxis()->LabelsOption("h");

   tex->Draw();
   tex2->Draw();
   tex3->Draw();

  cPerformance_ratio->Print((outputDirectory+"/PerformanceBDT_Ratio.pdf").c_str(),"pdf");
  cPerformance_ratio->Print((outputDirectory+"/PerformanceBDT_Ratio.png").c_str(),"png");
  cPerformance_ratio->Print((outputDirectory+"/PerformanceBDT_Ratio.root").c_str(),"root");
  cPerformance_ratio->SetLogz();
  cPerformance_ratio->Print((outputDirectory+"/PerformanceBDT_Ratio_Log.pdf").c_str(),"pdf");
  cPerformance_ratio->Print((outputDirectory+"/PerformanceBDT_Ratio_Log.png").c_str(),"png");
  cPerformance_ratio->Print((outputDirectory+"/PerformanceBDT_Ratio_Log.root").c_str(),"root");


  }

  return 0 ;



}

///////////////////////////////
std::vector<MatrixEntry> rocMatrix(std::vector<TTree*> InputTreeSignal,std::vector<TTree*> InputTreeBackground, std::vector<std::string> branchName, std::vector<std::string> reduceNameX,
                                   std::vector<std::string> reduceNameY, std::string cutSignal, std::string cutBackground, double signalEfficiencyTarget){

  //Loop on each branchName and make the signal total histogram
  std::vector<MatrixEntry> performanceValue ;

  for(size_t iVar = 0;  iVar < branchName.size(); iVar++){
    std::vector<TH1F*> signalHisto ;     signalHisto.resize(InputTreeSignal.size());
    std::vector<TH1F*> backgroundHisto ; backgroundHisto.resize(InputTreeBackground.size());
    for(size_t iTree = 0; iTree < InputTreeSignal.size(); iTree++){
      signalHisto.at(iTree) = new TH1F(Form("SHist_%s_%s_%d",reduceNameX.at(iVar).c_str(),reduceNameY.at(iVar).c_str(),int(iTree)),"",500,-1,1);
      InputTreeSignal.at(iTree)->Draw(Form("%s >> %s",branchName.at(iVar).c_str(),signalHisto.at(iTree)->GetName()),cutSignal.c_str(),"goff");
    }
    for(size_t iTree = 0; iTree < InputTreeBackground.size(); iTree++){
      backgroundHisto.at(iTree) = new TH1F(Form("BHist_%s_%s_%d",reduceNameX.at(iVar).c_str(),reduceNameY.at(iVar).c_str(),int(iTree)),"",500,-1,1);
      InputTreeBackground.at(iTree)->Draw(Form("%s >> %s",branchName.at(iVar).c_str(),backgroundHisto.at(iTree)->GetName()),cutBackground.c_str(),"goff");
    }

    TH1F* totalSignal = new TH1F(Form("SHist_%s_%s",reduceNameX.at(iVar).c_str(),reduceNameY.at(iVar).c_str()),"",500,-1,1); totalSignal->Sumw2();
    TH1F* totalBackground = new TH1F(Form("SHist_%s_%s",reduceNameX.at(iVar).c_str(),reduceNameY.at(iVar).c_str()),"",500,-1,1); totalBackground->Sumw2();
    for(size_t iSig = 0; iSig < signalHisto.size(); iSig++) totalSignal->Add(signalHisto.at(iSig));
    for(size_t iBack = 0; iBack < backgroundHisto.size(); iBack++) totalBackground->Add(backgroundHisto.at(iBack));

    // found signal efficiency building first the full roc
    TGraphErrors* rocVar = new TGraphErrors(); rocVar->SetName(Form("Graph_roc_%s_%s",reduceNameX.at(iVar).c_str(),reduceNameY.at(iVar).c_str()));
    for(int iBinX = 0; iBinX < totalSignal->GetNbinsX(); iBinX++){
      double errorNumSig = 0 ;
      double errorDenSig = 0 ;
      double numeratorSig   = totalSignal->IntegralAndError(iBinX,totalSignal->GetNbinsX()+1,errorNumSig);
      double denominatorSig = totalSignal->IntegralAndError(0,totalSignal->GetNbinsX()+1,errorDenSig);

      double effSig   = numeratorSig/denominatorSig;
      double errorSig = 1./(denominatorSig*denominatorSig)*(sqrt(numeratorSig*numeratorSig+errorDenSig*errorDenSig*(numeratorSig*numeratorSig)/(denominatorSig*denominatorSig)));

      double errorNumBack = 0 ;
      double errorDenBack = 0 ;
      double numeratorBack   = totalBackground->IntegralAndError(iBinX,totalBackground->GetNbinsX()+1,errorNumBack);
      double denominatorBack = totalBackground->IntegralAndError(0,totalBackground->GetNbinsX()+1,errorNumBack);

      double effBack   = numeratorBack/denominatorBack;
      double errorBack = 1./(denominatorBack*denominatorBack)*(sqrt(numeratorBack*numeratorBack+errorDenBack*errorDenBack*(numeratorBack*numeratorBack)/(denominatorBack*denominatorBack)));

      rocVar->SetPoint(iBinX+1,effSig,effBack);
      rocVar->SetPointError(iBinX+1,0.,1/(effBack*effBack)*sqrt(errorSig*errorSig+(effSig*effSig)/(effBack*effBack)*errorBack*errorBack));

    }

      MatrixEntry value ;
      value.binXName = reduceNameX.at(iVar);
      value.binYName = reduceNameY.at(iVar);
      value.value = rocVar->Eval(signalEfficiencyTarget);
      value.error = rocVar->GetHistogram()->GetBinError(rocVar->GetHistogram()->FindBin(signalEfficiencyTarget));
      performanceValue.push_back(value);
  }
  
  return performanceValue;					   
}

//////////////////////////////

int GetListOfMethods( TList & methods, TDirectory *dir){

  if (dir==0) dir = gDirectory;
  TIter mnext(dir->GetListOfKeys());
  TKey *mkey;
  methods.Clear();
  methods.SetOwner(kFALSE);
  UInt_t ni=0;
  while ((mkey = (TKey*)mnext())) { // make sure, that we only look at TDirectory with name Method_<xxx>                                                                      
    TString name = mkey->GetClassName();
    TClass *cl = gROOT->GetClass(name);
    if (cl->InheritsFrom("TDirectory")) {
      if (TString(mkey->GetName()).BeginsWith("Method_")) {
	methods.Add(mkey);
	ni++;
      }
    }
  }
  return ni;
}

// get a list of titles (i.e TDirectory) given a method dir                                                                                                           
int GetListOfTitles( TDirectory *rfdir, TList & titles ){
 UInt_t ni=0;
 if (rfdir==0) return 0;
 TList *keys = rfdir->GetListOfKeys();
 if(keys==0) {
    std::cout << "+++ Directory '" << rfdir->GetName() << "' contains no keys" << std::endl;
    return 0;
 }

 TIter rfnext(rfdir->GetListOfKeys());
 TKey *rfkey;
 titles.Clear();
 titles.SetOwner(kFALSE);
 while ((rfkey = (TKey*)rfnext())) { // make sure, that we only look at histograms                                                                                                       
   TClass *cl = gROOT->GetClass(rfkey->GetClassName());
   if (cl->InheritsFrom("TDirectory")) {
   titles.Add(rfkey);
   ni++;
   }
 }
 return ni;
}

// Next key iterator matching the className                                                                                                                    
TKey *NextKey(TIter & keyIter, TString className) {
 TKey *key  = (TKey *) keyIter.Next();
 TKey *rkey = 0;
 Bool_t loop = (key!=0);

 while (loop) {
  TClass *cl = gROOT->GetClass(key->GetClassName());
  if (cl->InheritsFrom(className.Data())) { loop = kFALSE;
       rkey = key;
  }
  else {
     key = (TKey *)keyIter.Next();
     if (key==0) loop = kFALSE;
  }
 }
 return rkey;

}
