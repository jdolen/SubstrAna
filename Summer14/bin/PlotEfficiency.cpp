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

double MinMassCut = 60 ;
double MaxMassCut = 105 ;
double WtagCut    = 0.5 ;
std::string WtagVariable = "tau2_beta_10[0]/tau1_beta_10[0]";

std::string ptName = "pt[0]";
std::string PUName = "pu";
std::string etaName = "eta[0]";


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

  // treename 
  std::string treeName ;
  if(Options.existsAs<std::string>("TreeName"))
    treeName = Options.getParameter<std::string>("TreeName");
  else{  std::cout<<" Exit from code, no input name tree provided "<<std::endl; std::exit(EXIT_FAILURE) ; }

  // preselectioncut string
  std::string PreselectionCut ;
  if(Options.existsAs<std::string>("PreselectionCut"))
    PreselectionCut = Options.getParameter<std::string>("PreselectionCut");
  else PreselectionCut = "pt[0] > 200" ;

  std::vector<double> JetPtBin;
  if(Options.existsAs<std::vector<double> >("JetPtBin"))
    JetPtBin = Options.getParameter<std::vector<double> >("JetPtBin");
  else{ std::cout<<" Exit from code, no  jet pT bins "<<std::endl; std::exit(EXIT_FAILURE) ; } 

  std::vector<double> JetPUBin;
  if(Options.existsAs<std::vector<double> >("JetPUBin"))
    JetPUBin = Options.getParameter<std::vector<double> >("JetPUBin");
  else{ std::cout<<" Exit from code, no  jet PU bins "<<std::endl; std::exit(EXIT_FAILURE) ; } 

  std::vector<double> JetEtaBin;
  if(Options.existsAs<std::vector<double> >("JetEtaBin"))
    JetPUBin = Options.getParameter<std::vector<double> >("JetEtaBin");
  else{ std::cout<<" Exit from code, no  jet Eta bins "<<std::endl; std::exit(EXIT_FAILURE) ; } 

  std::string outputFileDirectory ;
  if(Options.existsAs<std::string>("outputFileDirectory"))
    outputFileDirectory = Options.getParameter<std::string>("outputFileDirectory");
  else outputFileDirectory = "cfgTraining/outputEfficiency" ;

  std::vector<std::string> InputBackgroundFiles;
  if(Options.existsAs<std::vector<std::string> >("InputBackgroundFiles"))
    InputBackgroundFiles = Options.getParameter<std::vector<std::string> >("InputBackgroundFiles");
  else{ std::cout<<" Exit from code, no background files "<<std::endl; std::exit(EXIT_FAILURE) ; }

  std::vector<std::string> InputSignalFiles;
  if(Options.existsAs<std::vector<std::string> >("InputSignalFiles"))
    InputSignalFiles = Options.getParameter<std::vector<std::string> >("InputSignalFiles");
  else{ std::cout<<" Exit from code, no signal files "<<std::endl; std::exit(EXIT_FAILURE) ; }


  std::vector<edm::ParameterSet> InputMassVariables;
  if(Options.existsAs<std::vector<edm::ParameterSet> >("InputMassVariables"))
    InputMassVariables = Options.getParameter<std::vector<edm::ParameterSet> >("InputMassVariables");
  else{ std::cout<<" Exit from code, no mass observables "<<std::endl; std::exit(EXIT_FAILURE) ; }

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

   
   // histogram definition
   std::vector<std::vector<TH1F*> > denominatorBackgroundMass_vsPT ; denominatorBackgroundMass_vsPT.resize(inputBackgroundTrees.size());
   std::vector<std::vector<TH1F*> > numeratorBackgroundMass_vsPT ; numeratorBackgroundMass_vsPT.resize(inputBackgroundTrees.size());
   std::vector<std::vector<TH1F*> > numeratorBackgroundWtag_vsPT ; numeratorBackgroundWtag_vsPT.resize(inputBackgroundTrees.size());

   std::vector<std::vector<TH1F*> > denominatorSignalMass_vsPT ; denominatorSignalMass_vsPT.resize(inputSignalTrees.size());
   std::vector<std::vector<TH1F*> > numeratorSignalMass_vsPT ; numeratorSignalMass_vsPT.resize(inputSignalTrees.size());
   std::vector<std::vector<TH1F*> > numeratorSignalWtag_vsPT ; numeratorSignalWtag_vsPT.resize(inputSignalTrees.size());

   double* ptBin = &JetPtBin[0] ;

   for(size_t iTree = 0; iTree < inputBackgroundTrees.size() ; iTree++){
    for( size_t iMassVar = 0; iMassVar < InputMassVariables.size(); iMassVar++){
      denominatorBackgroundMass_vsPT.at(iTree).push_back(new TH1F(Form("back_deno_%s_%d_vsPT",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),int(iTree)),"",JetPtBin.size()-1,ptBin));
      numeratorBackgroundMass_vsPT.at(iTree).push_back(new TH1F(Form("back_num_%s_%d_vsPT",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),int(iTree)),"",JetPtBin.size()-1,ptBin));
      numeratorBackgroundWtag_vsPT.at(iTree).push_back(new TH1F(Form("back_num_wtag_%s_%d_vsPT",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),int(iTree)),"",JetPtBin.size()-1,ptBin));
      denominatorBackgroundMass_vsPT.at(iTree).back()->Sumw2();
      numeratorBackgroundMass_vsPT.at(iTree).back()->Sumw2();
      numeratorBackgroundWtag_vsPT.at(iTree).back()->Sumw2();
    }
   }   

   for(size_t iTree = 0; iTree < inputSignalTrees.size() ; iTree++){
    for( size_t iMassVar = 0; iMassVar < InputMassVariables.size(); iMassVar++){
      denominatorSignalMass_vsPT.at(iTree).push_back(new TH1F(Form("sig_deno_%s_%d_vsPT",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),int(iTree)),"",JetPtBin.size()-1,ptBin));
      numeratorSignalMass_vsPT.at(iTree).push_back(new TH1F(Form("num_%s_%d_vsPT",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),int(iTree)),"",JetPtBin.size()-1,ptBin));
      numeratorSignalWtag_vsPT.at(iTree).push_back(new TH1F(Form("num_wtag_%s_%d_vsPT",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),int(iTree)),"",JetPtBin.size()-1,ptBin));
      denominatorSignalMass_vsPT.at(iTree).back()->Sumw2();
      numeratorSignalMass_vsPT.at(iTree).back()->Sumw2();
      numeratorSignalWtag_vsPT.at(iTree).back()->Sumw2();
    }
   }  
  
   for(size_t iTree = 0; iTree < inputBackgroundTrees.size() ; iTree++){
    for( size_t iMassVar = 0; iMassVar < InputMassVariables.size(); iMassVar++){
      inputBackgroundTrees.at(iTree)->Draw(Form("%s >> %s ",ptName.c_str(),numeratorBackgroundMass_vsPT.at(iTree).at(iMassVar)->GetName()),Form("%s && %s > %f && %s < %f",PreselectionCut.c_str(),InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),MinMassCut,InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),MaxMassCut),"goff");
      inputBackgroundTrees.at(iTree)->Draw(Form("%s >> %s ",ptName.c_str(),denominatorBackgroundMass_vsPT.at(iTree).at(iMassVar)->GetName()),Form("%s",PreselectionCut.c_str()),"goff");
      inputBackgroundTrees.at(iTree)->Draw(Form("%s >> %s ",ptName.c_str(),numeratorBackgroundWtag_vsPT.at(iTree).at(iMassVar)->GetName()),Form("%s && %s > %f && %s < %f && %s < %f ",PreselectionCut.c_str(),InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),MinMassCut,InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),MaxMassCut,WtagVariable.c_str(),WtagCut),"goff");
    }
   }

   for(size_t iTree = 0; iTree < inputSignalTrees.size() ; iTree++){
    for( size_t iMassVar = 0; iMassVar < InputMassVariables.size(); iMassVar++){
      inputSignalTrees.at(iTree)->Draw(Form("%s >> %s ",ptName.c_str(),numeratorSignalMass_vsPT.at(iTree).at(iMassVar)->GetName()),Form("%s && %s > %f && %s < %f",PreselectionCut.c_str(),InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),MinMassCut,InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),MaxMassCut),"goff");
      inputSignalTrees.at(iTree)->Draw(Form("%s >> %s ",ptName.c_str(),denominatorSignalMass_vsPT.at(iTree).at(iMassVar)->GetName()),Form("%s",PreselectionCut.c_str()),"goff");

      inputSignalTrees.at(iTree)->Draw(Form("%s >> %s ",ptName.c_str(),numeratorSignalWtag_vsPT.at(iTree).at(iMassVar)->GetName()),Form("%s && %s > %f && %s < %f && %s < %f ",PreselectionCut.c_str(),InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),MinMassCut,InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str(),MaxMassCut,WtagVariable.c_str(),WtagCut),"goff");
    }
   }

   std::vector<TGraphAsymmErrors*> signalEfficiencyMass_vsPT ; 
   std::vector<TGraphAsymmErrors*> backgroundEfficiencyMass_vsPT ;
   std::vector<TGraphAsymmErrors*> signalEfficiencyWtag_vsPT ; 
   std::vector<TGraphAsymmErrors*> backgroundEfficiencyWtag_vsPT ;

   TH1F* denTemp ;   
   TH1F* numTempMass ;   
   TH1F* numTempWtag ;   

   numTempMass = new TH1F(Form("numTempMass"),"",JetPtBin.size()-1,ptBin);
   numTempWtag = new TH1F(Form("numTempWtag"),"",JetPtBin.size()-1,ptBin);
   denTemp = new TH1F(Form("denTemp"),"",JetPtBin.size()-1,ptBin);

   numTempMass->Sumw2();
   numTempWtag->Sumw2();
   denTemp->Sumw2();
    

   for(size_t iMassVar = 0; iMassVar < InputMassVariables.size(); iMassVar++){
     
     numTempMass->Reset();
     numTempWtag->Reset();
     denTemp->Reset();

     signalEfficiencyMass_vsPT.push_back(new TGraphAsymmErrors());
     signalEfficiencyMass_vsPT.back()->SetName(Form("signalEfficiencyMass_vsPT_%s",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str()));
     backgroundEfficiencyMass_vsPT.push_back(new TGraphAsymmErrors());
     backgroundEfficiencyMass_vsPT.back()->SetName(Form("backgroundEfficiencyMass_vsPT_%s",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str()));

     signalEfficiencyWtag_vsPT.push_back(new TGraphAsymmErrors());
     signalEfficiencyWtag_vsPT.back()->SetName(Form("signalEfficiencyWtag_vsPT_%s",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str()));
     backgroundEfficiencyWtag_vsPT.push_back(new TGraphAsymmErrors());
     backgroundEfficiencyWtag_vsPT.back()->SetName(Form("backgroundEfficiencyWtag_vsPT_%s",InputMassVariables.at(iMassVar).getParameter<std::string>("VariableName").c_str()));

     for(size_t iTree = 0; iTree < inputSignalTrees.size() ; iTree++){
       numTempMass->Add(numeratorSignalMass_vsPT.at(iTree).at(iMassVar));
       numTempWtag->Add(numeratorSignalWtag_vsPT.at(iTree).at(iMassVar));
       denTemp->Add(denominatorSignalMass_vsPT.at(iTree).at(iMassVar));
     }

     signalEfficiencyMass_vsPT.back()->BayesDivide(numTempMass,denTemp);
     signalEfficiencyWtag_vsPT.back()->BayesDivide(numTempWtag,denTemp);

     numTempMass->Reset();
     numTempWtag->Reset();
     denTemp->Reset();

     for(size_t iTree = 0; iTree < inputBackgroundTrees.size() ; iTree++){
       numTempMass->Add(numeratorBackgroundMass_vsPT.at(iTree).at(iMassVar));
       numTempWtag->Add(numeratorBackgroundWtag_vsPT.at(iTree).at(iMassVar));
       denTemp->Add(denominatorBackgroundMass_vsPT.at(iTree).at(iMassVar));
     }

     backgroundEfficiencyMass_vsPT.back()->BayesDivide(numTempMass,denTemp);
     backgroundEfficiencyWtag_vsPT.back()->BayesDivide(numTempWtag,denTemp);
   }

   TFile* output  = new TFile("output.root","RECREATE");

   for(size_t iSig = 0 ; iSig < signalEfficiencyMass_vsPT.size(); iSig++)
     signalEfficiencyMass_vsPT.at(iSig)->Write();

   for(size_t iSig = 0 ; iSig < signalEfficiencyWtag_vsPT.size(); iSig++)
     signalEfficiencyWtag_vsPT.at(iSig)->Write();

   for(size_t iSig = 0 ; iSig < backgroundEfficiencyMass_vsPT.size(); iSig++)
     backgroundEfficiencyMass_vsPT.at(iSig)->Write();

   for(size_t iSig = 0 ; iSig < backgroundEfficiencyWtag_vsPT.size(); iSig++)
     backgroundEfficiencyWtag_vsPT.at(iSig)->Write();
  
   output->Close();

   /*
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

   TH2F* frame_vsPT = new TH2F("frame_vsPT","",500,JetPtBin.at(0),JetPtBin.back(),500,0,1);
   frame_vsPT->SetLineWidth(2);
   frame_vsPT->SetMarkerStyle(21);
   frame_vsPT->SetMarkerSize(0.3);
   frame_vsPT->GetXaxis()->SetNdivisions(405);
   frame_vsPT->GetYaxis()->SetNdivisions(405);
   frame_vsPT->GetXaxis()->SetTitle("AK8 jet p_{T} (GeV)");
   frame_vsPT->GetXaxis()->SetLabelOffset(0.012);
   frame_vsPT->GetXaxis()->SetLabelSize(0.042);
   frame_vsPT->GetXaxis()->SetTitleSize(0.05);
   frame_vsPT->GetXaxis()->SetTitleOffset(1.05);
   frame_vsPT->GetYaxis()->SetTitle("Efficiency");
   frame_vsPT->GetYaxis()->SetLabelOffset(0.012);
   frame_vsPT->GetYaxis()->SetLabelSize(0.042);
   frame_vsPT->GetYaxis()->SetTitleSize(0.05);
   frame_vsPT->GetYaxis()->SetTitleOffset(1.25);

   TLatex *   tex = new TLatex(0.85,0.92," 13 TeV");
   tex->SetNDC();
   tex->SetTextAlign(31);
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   TLatex * tex2 = new TLatex(0.13,0.92,"CMS");
   tex2->SetNDC();
   tex2->SetTextFont(61);
   tex2->SetTextSize(0.04);
   tex2->SetLineWidth(2);
   TLatex * tex3 = new TLatex(0.2324,0.92,"Simulation Preliminary");
   tex3->SetNDC();
   tex3->SetTextFont(52);
   tex3->SetTextSize(0.0304);
   tex3->SetLineWidth(2);
   tex3->Draw();


   for( size_t iSignal = 0 ; iSignal<signalEfficiencyMass_vsPT.size(); iSignal++){
     cSignal_vsPT->cd();
     signalEfficiencyMass_vsPT.at(iSignal)->SetLineColor(kBlack); 
     signalEfficiencyMass_vsPT.at(iSignal)->SetMarkerColor(kBlack); 
     signalEfficiencyMass_vsPT.at(iSignal)->SetLineWidth(2); 
     signalEfficiencyMass_vsPT.at(iSignal)->SetMarkerStyle(20); 
     frame_vsPT->Draw("");
     signalEfficiencyMass_vsPT.at(iSignal)->Draw("apsame");

     signalEfficiencyWtag_vsPT.at(iSignal)->SetLineColor(kRed); 
     signalEfficiencyWtag_vsPT.at(iSignal)->SetMarkerColor(kRed); 
     signalEfficiencyWtag_vsPT.at(iSignal)->SetLineWidth(2); 
     signalEfficiencyWtag_vsPT.at(iSignal)->SetMarkerStyle(22); 
     signalEfficiencyWtag_vsPT.at(iSignal)->Draw("apsame");
          
     cSignal_vsPT->SaveAs(std::string(outputFileDirectory+"/"+InputMassVariables.at(iSignal).getParameter<std::string>("VariableName")+".png").c_str(),"png");
     cSignal_vsPT->SaveAs(std::string(outputFileDirectory+"/"+InputMassVariables.at(iSignal).getParameter<std::string>("VariableName")+".pdf").c_str(),"pdf");

   }

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
   */
      

   return 0;
}
 
