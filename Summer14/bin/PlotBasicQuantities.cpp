// ROOT objects                                                                                                                                                                       
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TChain.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TString.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TLatex.h"
#include "TLegend.h"
// c++ lybraries                                                                                                                                                                        
#include <ctime>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cstdlib>

using namespace std;

// function to fill a TChain with the list of input files to be processed                                                                                                                 
bool FillChain(TChain* chain, const std::string& inputFileList){
  std::ifstream inFile(inputFileList.c_str());
  std::string buffer;

  if(!inFile.is_open()){
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
    return false;
  }

  while(1){
    inFile >> buffer;
    if(!inFile.good()) break;
    chain->Add(buffer.c_str());
  }
  return true;
}

int main (int argc, char ** argv) {
 
  if (argc<4){
    cout << "Missing arguments!!!" <<endl;
    cout << "Usage: PlotBasicQuantities <input list 1> <input files list 2> <output file>" <<endl;
  }

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

 
  std::string inputFilesList1 = argv[1]; // input file list name                                                                                                                        
  std::string inputFilesList2 = argv[2]; // input file list name                                                                                                                        
  std::string fOut            = argv[3]; // output file name                                                                                                                             

  std::string cut = "pt[0] > 200 && abs(eta[0]) < 2.4 && is_MatchedToBoson[0] == 1";

  // --- Read list of files to be analyzed and fill TChain                                                                                                                             
  TChain* lTree1 = new TChain("puppi");
  TChain* lTree2 = new TChain("puppi");
    
  FillChain(lTree1,inputFilesList1);
  FillChain(lTree2,inputFilesList2);

  std::cout<<" tree1 --> number of entries "<<lTree1->GetEntries()<<std::endl;
  std::cout<<" tree2 --> number of entries "<<lTree2->GetEntries()<<std::endl;

  TFile* output = new TFile(fOut.c_str(),"RECREATE"); 
  output->cd();

  TH1F *npu_1      = new TH1F("npu_1","",20,20,80);  
  TH1F *jet_pt_1   = new TH1F("jet_pt_1","",35,200,600);  
  TH1F *jet_eta_1  = new TH1F("jet_eta_1","",20,-2.4,2.4);  
  TH1F *jet_mraw_1 = new TH1F("jet_mraw_1","",30,20,120);  
  TH1F *jet_mtrim_1         = new TH1F("jet_mtrim_1","",30,20,120);  
  TH1F *jet_mpruned_1       = new TH1F("jet_mpruned_1","",30,20,120);  
  TH1F *jet_msoftdrop_1     = new TH1F("jet_mpsoftdrop_1","",30,20,120);  
  TH1F *jet_QGLikelihood_1  = new TH1F("jet_QGLikelihood_1","",30,-1,1);  
  TH1F *jet_QGLikelihood_sub1_1 = new TH1F("jet_QGLikelihood_sub1_1","",30,-1,1);  
  TH1F *jet_QGLikelihood_sub2_1 = new TH1F("jet_QGLikelihood_sub2_1","",30,-1,1);  
  TH1F *jet_tau1_1 = new TH1F("jet_tau1_1","",25,0,1);  
  TH1F *jet_tau2_1 = new TH1F("jet_tau2_1","",25,0,1);  
  TH1F *jet_tau2tau1_1 = new TH1F("jet_tau2tau1_1","",25,0,1);  
  TH1F *jet_qjets_1 = new TH1F("jet_qjets_1","",25,0,1);  
  TH1F *jet_C2_1 = new TH1F("jet_C2_1","",25,0,1);  
  TH1F *jet_D2_1 = new TH1F("jet_D2_1","",25,0,1);  

  npu_1->Sumw2();
  jet_pt_1->Sumw2();
  jet_eta_1->Sumw2();
  jet_mraw_1->Sumw2();
  jet_mtrim_1->Sumw2();
  jet_mpruned_1->Sumw2();
  jet_msoftdrop_1->Sumw2();
  jet_QGLikelihood_1->Sumw2();
  jet_QGLikelihood_sub1_1->Sumw2();
  jet_QGLikelihood_sub2_1->Sumw2();
  jet_tau1_1->Sumw2();
  jet_tau2_1->Sumw2();
  jet_tau2tau1_1->Sumw2();
  jet_qjets_1->Sumw2();
  jet_C2_1->Sumw2();
  jet_D2_1->Sumw2();

  TH1F *npu_2      = new TH1F("npu_2","",20,20,80);  
  TH1F *jet_pt_2   = new TH1F("jet_pt_2","",35,200,600);  
  TH1F *jet_eta_2  = new TH1F("jet_eta_2","",20,-2.4,2.4);  
  TH1F *jet_mraw_2 = new TH1F("jet_mraw_2","",30,20,120);  
  TH1F *jet_mtrim_2         = new TH1F("jet_mtrim_2","",30,20,120);  
  TH1F *jet_mpruned_2       = new TH1F("jet_mpruned_2","",30,20,120);  
  TH1F *jet_msoftdrop_2     = new TH1F("jet_mpsoftdrop_2","",30,20,120);  
  TH1F *jet_QGLikelihood_2      = new TH1F("jet_QGLikelihood_2","",30,-1,1);  
  TH1F *jet_QGLikelihood_sub1_2 = new TH1F("jet_QGLikelihood_sub1_2","",30,-1,1);  
  TH1F *jet_QGLikelihood_sub2_2 = new TH1F("jet_QGLikelihood_sub2_2","",30,-1,1);  
  TH1F *jet_tau1_2 = new TH1F("jet_tau1_2","",25,0,1);  
  TH1F *jet_tau2_2 = new TH1F("jet_tau2_2","",25,0,1);  
  TH1F *jet_tau2tau1_2 = new TH1F("jet_tau2tau1_2","",25,0,1);  
  TH1F *jet_qjets_2 = new TH1F("jet_qjets_2","",25,0,1);  
  TH1F *jet_C2_2 = new TH1F("jet_C2_2","",25,0,1);  
  TH1F *jet_D2_2 = new TH1F("jet_D2_2","",25,0,1);  

  npu_2->Sumw2();
  jet_pt_2->Sumw2();
  jet_eta_2->Sumw2();
  jet_mraw_2->Sumw2();
  jet_mtrim_2->Sumw2();
  jet_mpruned_2->Sumw2();
  jet_msoftdrop_2->Sumw2();
  jet_QGLikelihood_2->Sumw2();
  jet_QGLikelihood_sub1_2->Sumw2();
  jet_QGLikelihood_sub2_2->Sumw2();
  jet_tau1_2->Sumw2();
  jet_tau2_2->Sumw2();
  jet_tau2tau1_2->Sumw2();
  jet_qjets_2->Sumw2();
  jet_C2_2->Sumw2();
  jet_D2_2->Sumw2();

  //make plot
  lTree1->Draw(Form("npu >> %s",npu_1->GetName()),cut.c_str(),"goff");  
  lTree2->Draw(Form("npu >> %s",npu_2->GetName()),cut.c_str(),"goff");  
  
  lTree1->Draw(Form("pt[0] >> %s",jet_pt_1->GetName()),cut.c_str(),"goff");  
  lTree2->Draw(Form("pt[0] >> %s",jet_pt_2->GetName()),cut.c_str(),"goff");  

  lTree1->Draw(Form("eta[0] >> %s",jet_eta_1->GetName()),cut.c_str(),"goff");  
  lTree2->Draw(Form("eta[0] >> %s",jet_eta_2->GetName()),cut.c_str(),"goff");  

  lTree1->Draw(Form("mraw[0] >> %s",jet_mraw_1->GetName()),cut.c_str(),"goff");  
  lTree2->Draw(Form("mraw[0] >> %s",jet_mraw_2->GetName()),cut.c_str(),"goff");  

  lTree1->Draw(Form("mtrimsafe_R_010_Pt_003[0] >> %s",jet_mtrim_1->GetName()),cut.c_str(),"goff");  
  lTree2->Draw(Form("mtrimsafe_R_010_Pt_003[0] >> %s",jet_mtrim_2->GetName()),cut.c_str(),"goff");  

  lTree1->Draw(Form("mprunedsafe_z_010_R_050[0] >> %s",jet_mpruned_1->GetName()),cut.c_str(),"goff");  
  lTree2->Draw(Form("mprunedsafe_z_010_R_050[0] >> %s",jet_mpruned_2->GetName()),cut.c_str(),"goff");  

  lTree1->Draw(Form("msoftdropsafe_z_010_beta_00[0] >> %s",jet_msoftdrop_1->GetName()),cut.c_str(),"goff");  
  lTree2->Draw(Form("msoftdropsafe_z_010_beta_00[0] >> %s",jet_msoftdrop_2->GetName()),cut.c_str(),"goff");  

  lTree1->Draw(Form("QGLikelihood_pr_z_010_R_050[0] >> %s",jet_QGLikelihood_1->GetName()),cut.c_str(),"goff");  
  lTree2->Draw(Form("QGLikelihood_pr_z_010_R_050[0] >> %s",jet_QGLikelihood_2->GetName()),cut.c_str(),"goff");  

  lTree1->Draw(Form("QGLikelihood_pr_sub1_z_010_R_050[0] >> %s",jet_QGLikelihood_sub1_1->GetName()),cut.c_str(),"goff");  
  lTree2->Draw(Form("QGLikelihood_pr_sub1_z_010_R_050[0] >> %s",jet_QGLikelihood_sub1_2->GetName()),cut.c_str(),"goff");  

  lTree1->Draw(Form("QGLikelihood_pr_sub2_z_010_R_050[0] >> %s",jet_QGLikelihood_sub2_1->GetName()),cut.c_str(),"goff");  
  lTree2->Draw(Form("QGLikelihood_pr_sub2_z_010_R_050[0] >> %s",jet_QGLikelihood_sub2_2->GetName()),cut.c_str(),"goff");  

  lTree1->Draw(Form("tau1_beta_10[0] >> %s",jet_tau1_1->GetName()),cut.c_str(),"goff");  
  lTree2->Draw(Form("tau1_beta_10[0] >> %s",jet_tau1_2->GetName()),cut.c_str(),"goff");  

  lTree1->Draw(Form("tau2_beta_10[0] >> %s",jet_tau2_1->GetName()),cut.c_str(),"goff");  
  lTree2->Draw(Form("tau2_beta_10[0] >> %s",jet_tau2_2->GetName()),cut.c_str(),"goff");  

  lTree1->Draw(Form("tau2_beta_10[0]/tau1_beta_10[0] >> %s",jet_tau2tau1_1->GetName()),cut.c_str(),"goff");  
  lTree2->Draw(Form("tau2_beta_10[0]/tau1_beta_10[0] >> %s",jet_tau2tau1_2->GetName()),cut.c_str(),"goff");  

  lTree1->Draw(Form("Qjets[0] >> %s",jet_qjets_1->GetName()),cut.c_str(),"goff");  
  lTree2->Draw(Form("Qjets[0] >> %s",jet_qjets_2->GetName()),cut.c_str(),"goff");  

  lTree1->Draw(Form("ecf_nP_1_beta_10[0]*ecf_nP_3_beta_10[0]/(ecf_nP_2_beta_10[0]*ecf_nP_2_beta_10[0]) >> %s",jet_C2_1->GetName()),cut.c_str(),"goff");  
  lTree2->Draw(Form("ecf_nP_1_beta_10[0]*ecf_nP_3_beta_10[0]/(ecf_nP_2_beta_10[0]*ecf_nP_2_beta_10[0]) >> %s",jet_C2_2->GetName()),cut.c_str(),"goff");  

  lTree1->Draw(Form("ecf_nP_3_beta_10[0]/(ecf_nP_2_beta_10[0]*ecf_nP_2_beta_10[0]*ecf_nP_2_beta_10[0]) >> %s",jet_D2_1->GetName()),cut.c_str(),"goff");  
  lTree2->Draw(Form("ecf_nP_3_beta_10[0]/(ecf_nP_2_beta_10[0]*ecf_nP_2_beta_10[0]*ecf_nP_2_beta_10[0]) >> %s",jet_D2_2->GetName()),cut.c_str(),"goff");  

  npu_1->Scale(1/npu_1->Integral());
  jet_pt_1->Scale(1/jet_pt_1->Integral());
  jet_eta_1->Scale(1/jet_eta_1->Integral());
  jet_mraw_1->Scale(1/jet_mraw_1->Integral());
  jet_mtrim_1->Scale(1/jet_mtrim_1->Integral());
  jet_mpruned_1->Scale(1/jet_mpruned_1->Integral());
  jet_msoftdrop_1->Scale(1/jet_msoftdrop_1->Integral());
  jet_QGLikelihood_1->Scale(1/jet_QGLikelihood_1->Integral());
  jet_QGLikelihood_sub1_1->Scale(1/jet_QGLikelihood_sub1_1->Integral());
  jet_QGLikelihood_sub2_1->Scale(1/jet_QGLikelihood_sub2_1->Integral());
  jet_tau1_1->Scale(1/jet_tau1_1->Integral());
  jet_tau2_1->Scale(1/jet_tau2_1->Integral());
  jet_tau2tau1_1->Scale(1/jet_tau2tau1_1->Integral());
  jet_qjets_1->Scale(1/jet_qjets_1->Integral());
  jet_C2_1->Scale(1/jet_C2_1->Integral());
  jet_D2_1->Scale(1/jet_D2_1->Integral());

  npu_2->Scale(1/npu_2->Integral());
  jet_pt_2->Scale(1/jet_pt_2->Integral());
  jet_eta_2->Scale(1/jet_eta_2->Integral());
  jet_mraw_2->Scale(1/jet_mraw_2->Integral());
  jet_mtrim_2->Scale(1/jet_mtrim_2->Integral());
  jet_mpruned_2->Scale(1/jet_mpruned_2->Integral());
  jet_msoftdrop_2->Scale(1/jet_msoftdrop_2->Integral());
  jet_QGLikelihood_2->Scale(1/jet_QGLikelihood_2->Integral());
  jet_QGLikelihood_sub1_2->Scale(1/jet_QGLikelihood_sub1_2->Integral());
  jet_QGLikelihood_sub2_2->Scale(1/jet_QGLikelihood_sub2_2->Integral());
  jet_tau1_2->Scale(1/jet_tau1_2->Integral());
  jet_tau2_2->Scale(1/jet_tau2_2->Integral());
  jet_tau2tau1_2->Scale(1/jet_tau2tau1_2->Integral());
  jet_qjets_2->Scale(1/jet_qjets_2->Integral());
  jet_C2_2->Scale(1/jet_C2_2->Integral());
  jet_D2_2->Scale(1/jet_D2_2->Integral());

  output->cd();

  // make the plot vs PT                                                                                                                                                                  
  TCanvas *cPlot = new TCanvas("cPlot","",180,52,550,550);

  cPlot->SetTicks();
  cPlot->SetFillColor(0);
  cPlot->SetBorderMode(0);
  cPlot->SetBorderSize(2);
  cPlot->SetTickx(1);
  cPlot->SetTicky(1);
  cPlot->SetRightMargin(0.05);
  cPlot->SetBottomMargin(0.12);
  cPlot->SetFrameBorderMode(0);

  TLatex *   tex = new TLatex(0.94,0.92," 13 TeV");
  tex->SetNDC();
  tex->SetTextAlign(31);
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  TLatex * tex2 = new TLatex(0.14,0.92,"CMS");
  tex2->SetNDC();
  tex2->SetTextFont(61);
  tex2->SetTextSize(0.04);
  tex2->SetLineWidth(2);
  TLatex * tex3 = new TLatex(0.236,0.92,"Simulation Preliminary");
  tex3->SetNDC();
  tex3->SetTextFont(52);
  tex3->SetTextSize(0.035);
  tex3->SetLineWidth(2);

  TLegend* Legend = new TLegend(0.58,0.62,0.934,0.82);
  Legend->SetBorderSize(0);
  Legend->SetFillColor(0);
  Legend->SetFillStyle(0);
  Legend->SetTextSize(0.031);
  Legend->SetTextFont(42);

  TLatex *   Banner = new TLatex(0.58,0.84,"t#bar{t} Pythia8, puppi jets");
  Banner->SetNDC();
  Banner->SetTextAlign(11);
  Banner->SetTextFont(42);
  Banner->SetTextSize(0.036);
  Banner->SetLineWidth(2);

  
  npu_1->SetLineColor(kBlack);
  npu_1->SetMarkerColor(kBlack);
  npu_1->SetLineWidth(2);
  npu_1->SetMarkerStyle(20);
  npu_1->GetXaxis()->SetTitle("N_{PU}");
  npu_1->GetYaxis()->SetTitle("Normalized Unity");
  npu_1->Draw("pe");

  npu_2->SetLineColor(kRed);
  npu_2->SetMarkerColor(kRed);
  npu_2->SetLineWidth(2);
  npu_2->SetMarkerStyle(22);
  npu_2->Draw("pesame");

  Legend->AddEntry(npu_1,"normal mixing","pl");
  Legend->AddEntry(npu_2,"pre mixing","pl");

  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  Banner->Draw("same");
  Legend->Draw("same");

  cPlot->SetName("npu");
  cPlot->Write();

  /// pt
  cPlot->Clear();

  jet_pt_1->SetLineColor(kBlack);
  jet_pt_1->SetMarkerColor(kBlack);
  jet_pt_1->SetLineWidth(2);
  jet_pt_1->SetMarkerStyle(20);
  jet_pt_1->GetXaxis()->SetTitle("AK8 jet p_{T} (GeV)");
  jet_pt_1->GetYaxis()->SetTitle("Normalized Unity");
  jet_pt_1->Draw("pe");

  jet_pt_2->SetLineColor(kRed);
  jet_pt_2->SetMarkerColor(kRed);
  jet_pt_2->SetLineWidth(2);
  jet_pt_2->SetMarkerStyle(22);
  jet_pt_2->Draw("pesame");


  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  Banner->Draw("same");
  Legend->Draw("same");

  cPlot->SetName("jet_pt");
  cPlot->Write();

  /// eta
  cPlot->Clear();

  jet_eta_1->SetLineColor(kBlack);
  jet_eta_1->SetMarkerColor(kBlack);
  jet_eta_1->SetLineWidth(2);
  jet_eta_1->SetMarkerStyle(20);
  jet_eta_1->GetXaxis()->SetTitle("AK8 jet #eta");
  jet_eta_1->GetYaxis()->SetTitle("Normalized Unity");
  jet_eta_1->Draw("pe");

  jet_eta_2->SetLineColor(kRed);
  jet_eta_2->SetMarkerColor(kRed);
  jet_eta_2->SetLineWidth(2);
  jet_eta_2->SetMarkerStyle(22);
  jet_eta_2->Draw("pesame");


  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  Banner->Draw("same");
  Legend->Draw("same");

  cPlot->SetName("jet_eta");
  cPlot->Write();

  /// mraw
  cPlot->Clear();

  jet_mraw_1->SetLineColor(kBlack);
  jet_mraw_1->SetMarkerColor(kBlack);
  jet_mraw_1->SetLineWidth(2);
  jet_mraw_1->SetMarkerStyle(20);
  jet_mraw_1->GetXaxis()->SetTitle("AK8 jet raw mass (GeV)");
  jet_mraw_1->GetYaxis()->SetTitle("Normalized Unity");
  jet_mraw_1->Draw("pe");

  jet_mraw_2->SetLineColor(kRed);
  jet_mraw_2->SetMarkerColor(kRed);
  jet_mraw_2->SetLineWidth(2);
  jet_mraw_2->SetMarkerStyle(22);
  jet_mraw_2->Draw("pesame");


  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  Banner->Draw("same");
  Legend->Draw("same");

  cPlot->SetName("jet_mraw");
  cPlot->Write();

  /// mtrim
  cPlot->Clear();

  jet_mtrim_1->SetLineColor(kBlack);
  jet_mtrim_1->SetMarkerColor(kBlack);
  jet_mtrim_1->SetLineWidth(2);
  jet_mtrim_1->SetMarkerStyle(20);
  jet_mtrim_1->GetXaxis()->SetTitle("AK8 jet M_{trim} (GeV)");
  jet_mtrim_1->GetYaxis()->SetTitle("Normalized Unity");
  jet_mtrim_1->Draw("pe");

  jet_mtrim_2->SetLineColor(kRed);
  jet_mtrim_2->SetMarkerColor(kRed);
  jet_mtrim_2->SetLineWidth(2);
  jet_mtrim_2->SetMarkerStyle(22);
  jet_mtrim_2->Draw("pesame");


  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  Banner->Draw("same");
  Legend->Draw("same");

  cPlot->SetName("jet_mtrim");
  cPlot->Write();

  /// mpruned
  cPlot->Clear();

  jet_mpruned_1->SetLineColor(kBlack);
  jet_mpruned_1->SetMarkerColor(kBlack);
  jet_mpruned_1->SetLineWidth(2);
  jet_mpruned_1->SetMarkerStyle(20);
  jet_mpruned_1->GetXaxis()->SetTitle("AK8 jet M_{pruned} (GeV)");
  jet_mpruned_1->GetYaxis()->SetTitle("Normalized Unity");
  jet_mpruned_1->Draw("pe");

  jet_mpruned_2->SetLineColor(kRed);
  jet_mpruned_2->SetMarkerColor(kRed);
  jet_mpruned_2->SetLineWidth(2);
  jet_mpruned_2->SetMarkerStyle(22);
  jet_mpruned_2->Draw("pesame");


  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  Banner->Draw("same");
  Legend->Draw("same");

  cPlot->SetName("jet_mpruned");
  cPlot->Write();

  /// msoftdrop
  cPlot->Clear();

  jet_msoftdrop_1->SetLineColor(kBlack);
  jet_msoftdrop_1->SetMarkerColor(kBlack);
  jet_msoftdrop_1->SetLineWidth(2);
  jet_msoftdrop_1->SetMarkerStyle(20);
  jet_msoftdrop_1->GetXaxis()->SetTitle("AK8 jet M_{softdrop} (GeV)");
  jet_msoftdrop_1->GetYaxis()->SetTitle("Normalized Unity");
  jet_msoftdrop_1->Draw("pe");

  jet_msoftdrop_2->SetLineColor(kRed);
  jet_msoftdrop_2->SetMarkerColor(kRed);
  jet_msoftdrop_2->SetLineWidth(2);
  jet_msoftdrop_2->SetMarkerStyle(22);
  jet_msoftdrop_2->Draw("pesame");


  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  Banner->Draw("same");
  Legend->Draw("same");

  cPlot->SetName("jet_msoftdrop");
  cPlot->Write();

  /// QGLikelihood
  cPlot->Clear();

  jet_QGLikelihood_1->SetLineColor(kBlack);
  jet_QGLikelihood_1->SetMarkerColor(kBlack);
  jet_QGLikelihood_1->SetLineWidth(2);
  jet_QGLikelihood_1->SetMarkerStyle(20);
  jet_QGLikelihood_1->GetXaxis()->SetTitle("AK8 jet Q/G Likelihood");
  jet_QGLikelihood_1->GetYaxis()->SetTitle("Normalized Unity");
  jet_QGLikelihood_1->Draw("pe");

  jet_QGLikelihood_2->SetLineColor(kRed);
  jet_QGLikelihood_2->SetMarkerColor(kRed);
  jet_QGLikelihood_2->SetLineWidth(2);
  jet_QGLikelihood_2->SetMarkerStyle(22);
  jet_QGLikelihood_2->Draw("pesame");


  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  Banner->Draw("same");
  Legend->Draw("same");

  cPlot->SetName("jet_QGLikelihood");
  cPlot->Write();

  /// QGLikelihood_sub1
  cPlot->Clear();

  jet_QGLikelihood_sub1_1->SetLineColor(kBlack);
  jet_QGLikelihood_sub1_1->SetMarkerColor(kBlack);
  jet_QGLikelihood_sub1_1->SetLineWidth(2);
  jet_QGLikelihood_sub1_1->SetMarkerStyle(20);
  jet_QGLikelihood_sub1_1->GetXaxis()->SetTitle("AK8 jet Q/G Likelihood sub1");
  jet_QGLikelihood_sub1_1->GetYaxis()->SetTitle("Normalized Unity");
  jet_QGLikelihood_sub1_1->Draw("pe");

  jet_QGLikelihood_sub1_2->SetLineColor(kRed);
  jet_QGLikelihood_sub1_2->SetMarkerColor(kRed);
  jet_QGLikelihood_sub1_2->SetLineWidth(2);
  jet_QGLikelihood_sub1_2->SetMarkerStyle(22);
  jet_QGLikelihood_sub1_2->Draw("pesame");


  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  Banner->Draw("same");
  Legend->Draw("same");

  cPlot->SetName("jet_QGLikelihood_sub1");
  cPlot->Write();

  /// QGLikelihood_sub2
  cPlot->Clear();

  jet_QGLikelihood_sub2_1->SetLineColor(kBlack);
  jet_QGLikelihood_sub2_1->SetMarkerColor(kBlack);
  jet_QGLikelihood_sub2_1->SetLineWidth(2);
  jet_QGLikelihood_sub2_1->SetMarkerStyle(20);
  jet_QGLikelihood_sub2_1->GetXaxis()->SetTitle("AK8 jet Q/G Likelihood sub2");
  jet_QGLikelihood_sub2_1->GetYaxis()->SetTitle("Normalized Unity");
  jet_QGLikelihood_sub2_1->Draw("pe");

  jet_QGLikelihood_sub2_2->SetLineColor(kRed);
  jet_QGLikelihood_sub2_2->SetMarkerColor(kRed);
  jet_QGLikelihood_sub2_2->SetLineWidth(2);
  jet_QGLikelihood_sub2_2->SetMarkerStyle(22);
  jet_QGLikelihood_sub2_2->Draw("pesame");


  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  Banner->Draw("same");
  Legend->Draw("same");

  cPlot->SetName("jet_QGLikelihood_sub2");
  cPlot->Write();

  /// tau1
  cPlot->Clear();

  jet_tau1_1->SetLineColor(kBlack);
  jet_tau1_1->SetMarkerColor(kBlack);
  jet_tau1_1->SetLineWidth(2);
  jet_tau1_1->SetMarkerStyle(20);
  jet_tau1_1->GetXaxis()->SetTitle("AK8 jet #tau_{1}");
  jet_tau1_1->GetYaxis()->SetTitle("Normalized Unity");
  jet_tau1_1->Draw("pe");

  jet_tau1_2->SetLineColor(kRed);
  jet_tau1_2->SetMarkerColor(kRed);
  jet_tau1_2->SetLineWidth(2);
  jet_tau1_2->SetMarkerStyle(22);
  jet_tau1_2->Draw("pesame");


  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  Banner->Draw("same");
  Legend->Draw("same");

  cPlot->SetName("jet_tau1");
  cPlot->Write();

  /// tau2
  cPlot->Clear();

  jet_tau2_1->SetLineColor(kBlack);
  jet_tau2_1->SetMarkerColor(kBlack);
  jet_tau2_1->SetLineWidth(2);
  jet_tau2_1->SetMarkerStyle(20);
  jet_tau2_1->GetXaxis()->SetTitle("AK8 jet #tau_{2}");
  jet_tau2_1->GetYaxis()->SetTitle("Normalized Unity");
  jet_tau2_1->Draw("pe");

  jet_tau2_2->SetLineColor(kRed);
  jet_tau2_2->SetMarkerColor(kRed);
  jet_tau2_2->SetLineWidth(2);
  jet_tau2_2->SetMarkerStyle(22);
  jet_tau2_2->Draw("pesame");


  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  Banner->Draw("same");
  Legend->Draw("same");

  cPlot->SetName("jet_tau2");
  cPlot->Write();

  /// tau2tau1
  cPlot->Clear();

  jet_tau2tau1_1->SetLineColor(kBlack);
  jet_tau2tau1_1->SetMarkerColor(kBlack);
  jet_tau2tau1_1->SetLineWidth(2);
  jet_tau2tau1_1->SetMarkerStyle(20);
  jet_tau2tau1_1->GetXaxis()->SetTitle("AK8 jet #tau_{2}/#tau_{1}");
  jet_tau2tau1_1->GetYaxis()->SetTitle("Normalized Unity");
  jet_tau2tau1_1->Draw("pe");

  jet_tau2tau1_2->SetLineColor(kRed);
  jet_tau2tau1_2->SetMarkerColor(kRed);
  jet_tau2tau1_2->SetLineWidth(2);
  jet_tau2tau1_2->SetMarkerStyle(22);
  jet_tau2tau1_2->Draw("pesame");


  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  Banner->Draw("same");
  Legend->Draw("same");

  cPlot->SetName("jet_tau2tau1");
  cPlot->Write();

  /// qjets
  cPlot->Clear();

  jet_qjets_1->SetLineColor(kBlack);
  jet_qjets_1->SetMarkerColor(kBlack);
  jet_qjets_1->SetLineWidth(2);
  jet_qjets_1->SetMarkerStyle(20);
  jet_qjets_1->GetXaxis()->SetTitle("AK8 jet #Gamma_{qjet}");
  jet_qjets_1->GetYaxis()->SetTitle("Normalized Unity");
  jet_qjets_1->Draw("pe");

  jet_qjets_2->SetLineColor(kRed);
  jet_qjets_2->SetMarkerColor(kRed);
  jet_qjets_2->SetLineWidth(2);
  jet_qjets_2->SetMarkerStyle(22);
  jet_qjets_2->Draw("pesame");


  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  Banner->Draw("same");
  Legend->Draw("same");

  cPlot->SetName("jet_qjets");
  cPlot->Write();

  /// C2
  cPlot->Clear();

  jet_C2_1->SetLineColor(kBlack);
  jet_C2_1->SetMarkerColor(kBlack);
  jet_C2_1->SetLineWidth(2);
  jet_C2_1->SetMarkerStyle(20);
  jet_C2_1->GetXaxis()->SetTitle("AK8 jet C_{2}");
  jet_C2_1->GetYaxis()->SetTitle("Normalized Unity");
  jet_C2_1->Draw("pe");

  jet_C2_2->SetLineColor(kRed);
  jet_C2_2->SetMarkerColor(kRed);
  jet_C2_2->SetLineWidth(2);
  jet_C2_2->SetMarkerStyle(22);
  jet_C2_2->Draw("pesame");


  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  Banner->Draw("same");
  Legend->Draw("same");

  cPlot->SetName("jet_C2");
  cPlot->Write();

  /// D2
  cPlot->Clear();

  jet_D2_1->SetLineColor(kBlack);
  jet_D2_1->SetMarkerColor(kBlack);
  jet_D2_1->SetLineWidth(2);
  jet_D2_1->SetMarkerStyle(20);
  jet_D2_1->GetXaxis()->SetTitle("AK8 jet D_{2}");
  jet_D2_1->GetYaxis()->SetTitle("Normalized Unity");
  jet_D2_1->Draw("pe");

  jet_D2_2->SetLineColor(kRed);
  jet_D2_2->SetMarkerColor(kRed);
  jet_D2_2->SetLineWidth(2);
  jet_D2_2->SetMarkerStyle(22);
  jet_D2_2->Draw("pesame");


  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  Banner->Draw("same");
  Legend->Draw("same");

  cPlot->SetName("jet_D2");
  cPlot->Write();

  return 0 ;

}


