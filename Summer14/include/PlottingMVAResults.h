// global TMVA style settings                                                                                                                                                             
#ifndef TMVA_TMVAGLOB
#define TMVA_TMVAGLOB

#include <iostream>
#include <vector>
#include <istream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <algorithm>

#include "TPad.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TSystem.h"
#include "TImage.h"
#include "TKey.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TObjArray.h"
#include "TClass.h"
#include "TH2F.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TNamed.h"
#include "TFormula.h"

#include "RVersion.h"
#include "TPaletteAxis.h"

// default definition for ROC curve color, style and width                                                                                                                                
static double color []      = {1,2,4,210,6,7,9,12,3,95,28,5,46,40};
static double linestyle []  = {1,1,8,9,7,2,5,9,1,6,10,8,7,1};
static double linewidth []  = {2,2,2,2,2,2,2,2,2,2,2,2,2,2};

static std::vector<double> vec_color(color, color + sizeof(color)/sizeof(double));
static std::vector<double> vec_linewidth(linewidth, linewidth + sizeof(linewidth)/sizeof(double));
static std::vector<double> vec_linestyle(linestyle, linestyle + sizeof(linestyle)/sizeof(double));

// significance Box class 
class significanceBox {
 
 public:

  //default constructor
  significanceBox(){
    Signal_ = NULL ;
    Background_ = NULL ;
    efficiencySignal_ = NULL ;
    efficiencyBackground_ = NULL ;  
    significance_ = NULL;
    maxSig_    = 0; 
    maxSigErr_ = 0; 
    maxbin_    = 0;
  }

  // defailt de-constructor
  ~significanceBox(){
    if(Signal_!=NULL) delete Signal_;
    if(efficiencySignal_!=NULL) delete efficiencySignal_;
    if(Background_!=NULL) delete Background_;
    if(efficiencyBackground_!=NULL) delete efficiencyBackground_;
    if(significance_!=NULL) delete significance_ ;
  } 

  // public members 
 
  // signal and background distribution
  TH1* Signal_;
  TH1* Background_;

  // signal and background efficiency
  TH1* efficiencySignal_;
  TH1* efficiencyBackground_;

  // chosen significance
  TH1F* significance_;  

  // method Name and Title
  TString methodName_;
  TString methodTitle_;

  // maxSignificance and error
  double maxSig_ ;
  double maxSigErr_ ;
  
  // bin related to the maximum
  int    maxbin_ ;

};

// class for global manipulation of TMVA output + ROC, correlation, output and significance plots
class TMVAGlob {

  public:  
  
  // default constructor
  TMVAGlob(); 
  // default deconstructor
  ~TMVAGlob();
 
  // type of plot for TMVA output
  enum TypeOfPlot { kId = 0, kNorm, kDecorrelated, kPCA, kGaussDecorr, kNumOfMethods };

  // type of plot for TMVA output
  enum HistType { MVAType = 0, ProbaType = 1, RarityType = 2, CompareType = 3 };

  // type of significance to plot
  enum SignificanceType { SoverB = 0, SoverSqrtB = 1, SoverSqrtSB = 2, Pvalue = 3};

  //////////////////////////////////////// Input Files

  // open one input file given the name
  void openFileInput( const TString& fin );

  // open a set of input file and store them in the class
  void openFileInput( const std::vector<std::string> & fin );

  // get back a single file from the name
  TFile* GetInputFile (const TString& fin );

  // get back the whole set of file considered
  std::vector<TFile*> GetInputFile();


  ////////////////////////////////////// Canvas methods

  // create canvas and frame for the ROC plot
  void CreateCanvasandFrameROC(const double & ptMin, const double & ptMax, const double & puMin, const double & puMax);
  // destroy exsisting canvas
  void DestroyCanvases(); 
  // Print ROC Curve plot
  void PrintImageROC( TDirectory* dir = 0 , const std::string & outputPlotDirectory = "");    
  // Print a generic plot
  void PrintImage   ( TCanvas* c = 0, const std::string & outputPlotDirectory = "");    


  /////////// Normalize a given histogram
  void NormalizeHist ( TH1F* h );
  void NormalizeHists( TH1F* sig, TH1F* bkg = 0 );

  /////////// Get Method Name and Title
  void GetMethodName ( TString & name, TKey * mkey );
  void GetMethodTitle( TString & name, TKey * ikey );
  void GetMethodName ( TString & name, TDirectory * mdir );
  void GetMethodTitle( TString & name, TDirectory * idir );

  //// GetList of Methods, Title and Jobs
  int GetListOfMethods( TList & methods, TDirectory *dir = 0 );
  int GetListOfTitles ( TDirectory *rfdir, TList & titles );
  int GetListOfTitles ( TString & methodName, TList & titles, TDirectory *dir = 0 );

  ///// Get the name of the input variables and classes
  std::vector<TString> GetInputVariableNames(TDirectory *dir );
  std::vector<TString> GetClassNames(TDirectory *dir );

  // given the iterator on the file and a class name, get the next key
  TKey *NextKey   ( TIter & keyIter, TString className);

  // get list of keys of the root file
  int  GetListOfKeys( TList& keys, TString inherits, TDirectory *dir = 0 );
 
  // get the inputVariable directory and the correlation plot directory
  TDirectory *GetInputVariablesDir  ( TMVAGlob::TypeOfPlot type, TDirectory *dir = 0 );
  TDirectory *GetCorrelationPlotsDir( TMVAGlob::TypeOfPlot type, TDirectory *dir = 0 );

  /// banner inside the ROC plot
  void banner4Plot (const bool & isLabel, const float & ptMin, const float & ptMax, const float & puMin, const float & puMax);
  /// Set the name of the method
  void SetMethodName (const std::vector<std::string> & MethodName);
  /// Set input files
  void SetInputFiles (const std::vector<TFile*> & inputFiles){ inputFiles_ = inputFiles; };

  // Read signal and background histogramm and effciency given an input file
  void ReadHistograms(TFile* file);
  // set the formula for the significance
  void SetFormula(const TString & f) { fFormula_ = f; }

  // Get the formula back for calculation or latex banner
  TString GetFormula();
  TString GetFormulaString() { return fFormula_; }
  TString GetLatexFormula();

  // set if use just signal or bkg efficiency or multiply them for the expected number of events before the training to get the real significance
  void SetSignalType    (const bool & type = false) { signalType_ = type; }
  void SetBackgroundType(const bool & type = false) { backgroundType_ = type; }

  // plot ROC curve
  void plotROCs (TDirectory* dir, const int & = 0, const double & minPTbin = 200, const double & maxPTbin = 1000, const double & minPU = 0, const double & maxPU = 1000, const std::string & outputPlotDirectory = "");
  // plot correlation Matrix
  void plotCorrelationMatrix (TFile* inputFile = 0, const std::string & inputName = "", const std::string & outputPlotDirectory = "");
  // plot output distribution
  void plotMVAs( TFile*inputFile = 0, const std::string & inputName = "", HistType htype = MVAType, const std::string & outputPlotDirectory = "");
  // plot signficance for different formula using or not the expected number of signal and background events
  void plotSignificance (TFile* inputFile = 0, const std::string & inputName = "", SignificanceType stype = SoverSqrtB, const double & numberSignalEvents = 0., 
                         const double & numberBackgroundEvents = 0.,const bool & UseSignalEfficiency = false, const bool & UseBakgroundEfficiency = false, 
                         const std::string & outputPlotDirectory = "");

 private:
  
  // indexes for color, method, mvasplot
  int color_index ;
  int method_index;
  int mvas_index;
  int thisMethod_;

  // list of the input files, inputMethods
  std::vector<TFile*> inputFiles_ ;
  std::vector<std::string> inputMethodName_ ;
  std::vector<std::string> originalMethodName_ ;

  // canvas for ROC, correlation, MVA and significance
  TCanvas* cROC_ ;
  TCanvas* cROCLog_ ;
  TCanvas* cCorrelationSignal_;
  TCanvas* cCorrelationBackground_;
  TCanvas* cMVAs_ ;
  TCanvas* cSignificance_ ;

  TH2F*    frameROC_ ; 
  TH2F*    frameROCLog_ ; 

  TLegend* legROC_ ;
 
  // formula string
  TString fFormula_ ;

  // store if use just efficiency or also the S and B yields
  bool    signalType_ ;
  bool    backgroundType_ ;

  // significance box tool
  std::vector<significanceBox*>*  fInfoList_;

  TH1* histoSignal_;
  TH1* histoBackground_;
  TH1* effBackground_;
  TH1* effSignal_;


};

#endif
