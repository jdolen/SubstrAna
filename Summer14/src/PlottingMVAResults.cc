#include <iomanip>
#include "../include/PlottingMVAResults.h"

// default contructor
TMVAGlob::TMVAGlob(){

  // ROC plots 
  cROC_        = NULL;
  frameROC_    = NULL;
  cROCLog_     = NULL;
  frameROCLog_ = NULL;
  legROC_      = NULL;

  // MVA output, correlation and Significance 
  cMVAs_                  = NULL ;
  cCorrelationSignal_     = NULL ;
  cCorrelationBackground_ = NULL ;
  cSignificance_          = NULL ;

  // method informations 
  color_index  = 0;
  method_index = 0;
  mvas_index   = 0;
  thisMethod_  = 0 ;
  fInfoList_   = NULL ;

  // histograms
  histoSignal_     = NULL ;
  histoBackground_ = NULL ;
  effSignal_       = NULL ;
  effBackground_   = NULL ;

  signalType_     = false ;
  backgroundType_ = false ;

  (*this).DestroyCanvases();

}

TMVAGlob::~TMVAGlob(){

  if(cROC_    !=NULL)     delete cROC_;
  if(frameROC_!=NULL)     delete frameROC_;
  if(cROCLog_    !=NULL)  delete cROCLog_;
  if(frameROCLog_!=NULL)  delete frameROCLog_;
  if(legROC_  !=NULL)     delete legROC_;

  if(cMVAs_ != NULL)                  delete cMVAs_;
  if(cCorrelationSignal_ != NULL)     delete cCorrelationSignal_;
  if(cCorrelationBackground_ != NULL) delete cCorrelationBackground_;
  if(cSignificance_!= NULL)           delete cSignificance_;

  if(histoSignal_!=NULL)     delete histoSignal_ ;
  if(histoBackground_!=NULL) delete histoBackground_;
  if(effSignal_!=NULL)       delete effSignal_ ;
  if(effBackground_!=NULL)   delete effBackground_;

  if(fInfoList_!=NULL){ fInfoList_->clear();  
                        delete fInfoList_; 
  }

  inputFiles_.clear() ;
  inputMethodName_.clear() ;
    
}

// destroy the full list of the exsisting canvas
void TMVAGlob::DestroyCanvases(){

  TList* loc = (TList*)gROOT->GetListOfCanvases();
  TListIter itc(loc);
  TObject *o(0);
  while ((o = itc())) delete o;
  return;

}


// checks if file with name "fin" is already open, and if not opens one
void TMVAGlob::openFileInput( const TString& fin ){

 TFile* file = gDirectory->GetFile();
 if (file ==0 || fin != file->GetName()) {
   if (file != 0) {
     gROOT->cd();
     file->Close();
   }
  std::cout << "--- Opening root file " << fin << " in read mode" << std::endl;
  file = TFile::Open( fin, "READ" );
 }
 else file = gDirectory->GetFile();
      
 inputFiles_.push_back(file);

 return ;

}

// oprn file from a list of file names 
void TMVAGlob::openFileInput( const std::vector<std::string> & fin ){

 for(size_t iFile = 0 ; iFile < fin.size(); iFile++){
  std::cout << "--- Opening root file " << fin.at(iFile) << " in read mode" << std::endl;
  inputFiles_.push_back(new TFile( fin.at(iFile).c_str(), "READ"));
 }

 return ;

}

// Get the input file from the name
TFile* TMVAGlob::GetInputFile(const TString& fin){

  for(size_t iFile = 0 ; iFile < inputFiles_.size(); iFile++){
    if(inputFiles_.at(iFile)->GetName() == fin) 
      return inputFiles_.at(iFile) ;
  }

  return 0 ;
}

// get input file vector and not a single file
std::vector<TFile*> TMVAGlob::GetInputFile(){
  return inputFiles_ ;
}


// Next key iterator matching the className
TKey *TMVAGlob::NextKey(TIter & keyIter, TString className) {

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

// used to create output file for canvas
void TMVAGlob::PrintImageROC(TDirectory* dir, const std::string & outputPlotDirectory){

  TString fname = outputPlotDirectory + "ROCcurve";

  if (TString(dir->GetName()).Contains("multicut")){
    TString fprepend(dir->GetName());
    fprepend.ReplaceAll("multicutMVA_","");
    fname = outputPlotDirectory + fprepend + "_" + "ROCcurve";
  }

  if (cROC_ == NULL){
     std::cout << "*** Error in TMVAGlob::PrintImageROC: canvas is NULL" << std::endl; 
     return ;
  }
  
  // create directory if not existing
  TString f = fname;
  TString dirName = f.Remove( f.Last( '/' ), f.Length() - f.Last( '/' ) );
  gSystem->mkdir( dirName );

  TString pngName = fname + ".png";
  TString pdfName = fname + ".pdf";
  TString rootName = fname + ".root";
  cROC_->cd();
  cROC_->Print(pdfName);
  cROC_->Print(pngName);
  cROC_->Print(rootName);

  TString pngNameLog = fname + "_Log.png";
  TString pdfNameLog = fname + "_Log.pdf";
  TString rootNameLog = fname + "_Log.root";
  cROCLog_->cd();
  cROCLog_->SetLogy();
  cROCLog_->Print(pdfNameLog);
  cROCLog_->Print(pngNameLog);
  cROCLog_->Print(rootNameLog);

  return ;
}

// function to normalize a histo also if it has weights
void TMVAGlob::NormalizeHist(TH1F* h) { 

 if (h==0) return;
 if (h->GetSumw2N() == 0) h->Sumw2();
 if (h->GetSumOfWeights()!=0) {
      Float_t dx = (h->GetXaxis()->GetXmax()-h->GetXaxis()->GetXmin())/h->GetNbinsX();
      h->Scale( 1.0/h->GetSumOfWeights()/dx );
  }
 return ;
}

void TMVAGlob::NormalizeHists( TH1F* sig, TH1F* bkg) {

 if (sig && sig->GetSumw2N() == 0) sig->Sumw2();
 if (bkg && bkg->GetSumw2N() == 0) bkg->Sumw2();
      
 if(sig->GetSumOfWeights()!=0) {
     Float_t dx = (sig->GetXaxis()->GetXmax()-sig->GetXaxis()->GetXmin())/sig->GetNbinsX();
     sig->Scale( 1.0/sig->GetSumOfWeights()/dx );
 }

 if (bkg != 0 && bkg->GetSumOfWeights()!=0) {
     Float_t dx = (bkg->GetXaxis()->GetXmax()-bkg->GetXaxis()->GetXmin())/bkg->GetNbinsX();
     bkg->Scale( 1.0/bkg->GetSumOfWeights()/dx);
 }

 return ;

}

// the following are tools to help handling different methods and titles
void TMVAGlob::GetMethodName( TString & name, TKey * mkey ) {
   if (mkey == 0) return;
   name = mkey->GetName();
   name.ReplaceAll("Method_","");
   return ;
}

void TMVAGlob::GetMethodTitle( TString & name, TKey * ikey ) {
   if (ikey == 0) return;
   name = ikey->GetName();
   return ;
}

void TMVAGlob::GetMethodName( TString & name, TDirectory * mdir ) {
   if (mdir == 0) return;
   name = mdir->GetName();
   name.ReplaceAll("Method_","");
   return ;
}

void TMVAGlob::GetMethodTitle( TString & name, TDirectory * idir ) {
   if (idir == 
0) return;
   name = idir->GetName();
   return ;

}

// get a list of keys with a given inheritance --> the list contains TKey objects
int TMVAGlob::GetListOfKeys( TList& keys, TString inherits, TDirectory *dir){

 if (dir==0) dir = gDirectory;

 TIter mnext(dir->GetListOfKeys());
 TKey *mkey;
 keys.Clear();
 keys.SetOwner(kFALSE);
 UInt_t ni=0;
 
 while ((mkey = (TKey*)mnext())) {         
   // make sure, that we only look at TDirectory with name Method_<xxx>
   TClass *cl = gROOT->GetClass(mkey->GetClassName());
   if (cl->InheritsFrom(inherits)) {
      keys.Add(mkey);
      ni++;
   }
 }
 
 return ni;
}

// Get a vector of string with all the class name object inside a directory   
std::vector<TString> TMVAGlob::GetClassNames(TDirectory *dir ){      
      
 TIter next(dir->GetListOfKeys());
 TKey* key = 0;
 std::vector<TString> names;
      
 while ((key = (TKey*)next())) {
   if (key->GetCycle() != 1) continue;
   TClass *cl = gROOT->GetClass(key->GetClassName());
   if (!cl->InheritsFrom("TH1")) continue;
   TString name(key->GetName());
   name.ReplaceAll("_Deco","");
   name.ReplaceAll("_Gauss","");
   name.ReplaceAll("_PCA","");
   name.ReplaceAll("_Id","");
   name.ReplaceAll("_vs_","");
   char c = '_';
   Int_t pos = name.Last(c);
   name.Remove(0,pos+1);
         
   Bool_t hasname = false;
   std::vector<TString>::const_iterator iter = names.begin();
   while(iter != names.end()){
       if(name.CompareTo(*iter)==0)
          hasname=true;
       iter++;
   }
   if(!hasname) names.push_back(name);
 }
 return names;
}



// fill a TList with all the methods
int TMVAGlob::GetListOfMethods( TList & methods, TDirectory *dir){

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
int TMVAGlob::GetListOfTitles( TDirectory *rfdir, TList & titles ){

 UInt_t ni=0;
 if (rfdir==0) return 0;

 TList *keys = rfdir->GetListOfKeys();

 if (keys==0) {
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


// get the list of all titles for a given method --> if the input dir is 0, gDirectory is used returns a list of keys
int TMVAGlob::GetListOfTitles( TString & methodName, TList & titles, TDirectory *dir){

  if (dir==0) dir = gDirectory;
  TDirectory* rfdir = (TDirectory*)dir->Get( methodName );
  if (rfdir==0) {
    std::cout << "+++ Could not locate directory '" << methodName << std::endl;
    return 0;
  }

  return (*this).GetListOfTitles(rfdir,titles);
}


// return the name of the directory which contains output MVA plots 
TDirectory *TMVAGlob::GetInputVariablesDir( TMVAGlob::TypeOfPlot type, TDirectory *dir){

 // get the InputVariables directory
 const TString directories[TMVAGlob::kNumOfMethods] = { "InputVariables_Id",
                                                        "InputVariables_Deco",
                                                        "InputVariables_PCA",
                                                        "InputVariables_Gauss_Deco" };
 if (dir==0) dir = gDirectory;

 // get top dir containing all hists of the variables
 dir = (TDirectory*)gDirectory->Get(directories[type]);
 if (dir==0) {
   std::cout << "+++ Could not locate input variable directory '" << directories[type] << std::endl;
   return 0;
 }
 
 return dir;

}


// return the name of the directory which correlation plot
TDirectory *TMVAGlob::GetCorrelationPlotsDir( TMVAGlob::TypeOfPlot type, TDirectory *dir){

  if (dir==0) dir = (*this).GetInputVariablesDir(type);
  if (dir==0) return 0;
  TDirectory* corrdir = (TDirectory*)dir->Get( "CorrelationPlots" );
  if (corrdir==0) {
    std::cout << "+++ Could not find CorrelationPlots directory 'CorrelationPlots'" << std::endl;
    return 0;
  }
  return corrdir;
}

// Produce a banner for ROC plots 
void TMVAGlob::banner4Plot (const bool & isLabel, const float & ptMin, const float & ptMax,const float & puMin, const float & puMax ){

  //TPaveText* pt = new TPaveText(.36,0.71,.43,.88,"NDC");
  TPaveText* pt = new TPaveText(.23,0.26,.30,.41,"NDC"); 

  pt->AddText("AK R = 0.8");
  TString BoostLegend ; BoostLegend.Form("%d < p_{T} < %d GeV",int(ptMin),int(ptMax));
  pt->AddText(BoostLegend.Data());
  TString PULegend ; PULegend.Form("%d < N_{PU} < %d ",int(puMin),int(puMax));
  pt->AddText(PULegend.Data());
  pt->AddText("|#eta| < 2.5");

  pt->SetBorderSize(0);
  pt->SetFillColor(0);
  pt->SetTextSize(0.027);
  pt->SetFillStyle(0);
  pt->SetLineColor(1);
  pt->SetLineStyle(1);
  pt->SetLineWidth(0);
  pt->SetMargin(0);
  pt->SetTextFont(42);
  pt->SetShadowColor(0);
  pt->Draw();
}


// create canvas, frame and legend for ROC plot
void TMVAGlob::CreateCanvasandFrameROC(const double & minPTbin, const double & maxPTbin, const double & puMin, const double & puMax){
                               

  cROC_ = new TCanvas((std::string("cROC")+std::string(inputFiles_.at(0)->GetName())).c_str(),"cROC",180,52,550,550);

  cROC_->cd();
  cROC_->SetTicks();
  cROC_->SetFillColor(0);
  cROC_->SetBorderMode(0);
  cROC_->SetBorderSize(2);
  
  cROC_->SetTickx(1);
  cROC_->SetTicky(1);
  cROC_->SetRightMargin(0.05);
  cROC_->SetBottomMargin(0.12);
  cROC_->SetFrameBorderMode(0);


  frameROC_ = new TH2F((std::string("frameROC")+std::string(inputFiles_.at(0)->GetName())).c_str(),"",500,0,1,500,0,1);
  frameROC_->SetLineWidth(2);
  frameROC_->SetMarkerStyle(21);
  frameROC_->SetMarkerSize(0.3);
  frameROC_->GetXaxis()->SetNdivisions(405);
  frameROC_->GetYaxis()->SetNdivisions(405);
  frameROC_->GetXaxis()->SetTitle("#varepsilon_{sig}");
  frameROC_->GetXaxis()->SetLabelOffset(0.012);
  frameROC_->GetXaxis()->SetLabelSize(0.042);
  frameROC_->GetXaxis()->SetTitleSize(0.05);
  frameROC_->GetXaxis()->SetTitleOffset(1.05);
  frameROC_->GetYaxis()->SetTitle("#varepsilon_{bkg}");
  frameROC_->GetYaxis()->SetLabelOffset(0.012);
  frameROC_->GetYaxis()->SetLabelSize(0.042);
  frameROC_->GetYaxis()->SetTitleSize(0.05);
  frameROC_->GetYaxis()->SetTitleOffset(1.25);
  frameROC_->Draw("");

  banner4Plot(false,minPTbin,maxPTbin,puMin,puMax);

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
  tex = new TLatex(0.2324,0.92,"Simulation Preliminary");
  tex->SetNDC();
  tex->SetTextFont(52);
  tex->SetTextSize(0.0304);
  tex->SetLineWidth(2);
  tex->Draw();
  
  cROCLog_ = new TCanvas((std::string("cROCLog")+std::string(inputFiles_.at(0)->GetName())).c_str(),"cROC",180,52,550,550);

  cROCLog_->cd();
  cROCLog_->SetTicks();
  cROCLog_->SetFillColor(0);
  cROCLog_->SetBorderMode(0);
  cROCLog_->SetBorderSize(2);
  
  cROCLog_->SetTickx(1);
  cROCLog_->SetTicky(1);
  cROCLog_->SetRightMargin(0.05);
  cROCLog_->SetBottomMargin(0.12);
  cROCLog_->SetFrameBorderMode(0);

  frameROCLog_ = new TH2F((std::string("frameROCLog")+std::string(inputFiles_.at(0)->GetName())).c_str(),"",500,0.0001,1,500,0.0001,1);
  frameROCLog_->SetLineWidth(2);
  frameROCLog_->SetMarkerStyle(21);
  frameROCLog_->SetMarkerSize(0.3);
  frameROCLog_->GetXaxis()->SetNdivisions(405);
  frameROCLog_->GetYaxis()->SetNdivisions(405);
  frameROCLog_->GetXaxis()->SetTitle("#varepsilon_{sig}");
  frameROCLog_->GetXaxis()->SetLabelOffset(0.012);
  frameROCLog_->GetXaxis()->SetLabelSize(0.042);
  frameROCLog_->GetXaxis()->SetTitleSize(0.05);
  frameROCLog_->GetXaxis()->SetTitleOffset(1.05);
  frameROCLog_->GetYaxis()->SetTitle("#varepsilon_{bkg}");
  frameROCLog_->GetYaxis()->SetLabelOffset(0.012);
  frameROCLog_->GetYaxis()->SetLabelSize(0.042);
  frameROCLog_->GetYaxis()->SetTitleSize(0.05);
  frameROCLog_->GetYaxis()->SetTitleOffset(1.25);
  frameROCLog_->Draw("");

  banner4Plot(false,minPTbin,maxPTbin,puMin,puMax);

  tex = new TLatex(0.85,0.92," 13 TeV");
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
  tex = new TLatex(0.2324,0.92,"Simulation Preliminary");
  tex->SetNDC();
  tex->SetTextFont(52);
  tex->SetTextSize(0.0304);
  tex->SetLineWidth(2);
  tex->Draw();
  
  legROC_ = new TLegend(0.17,0.43,0.56,0.87,NULL,"brNDC");

  legROC_->SetBorderSize(0);
  legROC_->SetTextSize(0.033);
  legROC_->SetTextFont(42);
  legROC_->SetLineColor(1);
  legROC_->SetLineStyle(1);
  legROC_->SetLineWidth(2);
  legROC_->SetFillColor(19);
  legROC_->SetFillStyle(0);
}

// methods for plot the efficiency (ROC) curve for a given inputFile
void TMVAGlob::plotROCs (TDirectory* dir, const int & plotType, const double & minPTbin, const double & maxPTbin, const double & puMin, const double & puMax, const std::string & outputPlotDirectory){

  // Plot the ROC curve with a proper style from root file originated by TMVA                                                                                                         
  if(cROC_==NULL) (*this).CreateCanvasandFrameROC(minPTbin,maxPTbin,puMin,puMax); 
      
  TH1F *h ; 

  for(size_t iFile = 0; iFile <  inputFiles_.size();  iFile++){
 
   inputFiles_.at(iFile)->cd();
   TList TrainingMethods;
   TList hists;
   int res = (*this).GetListOfMethods(TrainingMethods);
   if(res == 0) std::cout<<" TMVAGlob::plotEfficiency --> no methods found "<<std::endl ;
   TIter next(&TrainingMethods);
   TKey *key = 0, *hkey = 0;

   // loop over all methods stored in the TList TrainingMethods                                                                                                              
   while ((key = (TKey*)next())) {
     TDirectory *myDir = (TDirectory*)key->ReadObj();
     TList Titles;
     int nTitles = (*this).GetListOfTitles(myDir,Titles); // get the titles list for eack method                                                                                          
     if(nTitles == 0) std::cout<<" TMVAGlob::plotEfficiency --> no titles found "<<std::endl ;
     TIter nextTitle(&Titles);
     TKey *titkey = 0;
     TDirectory *titDir = 0;
     while ((titkey = (*this).NextKey(nextTitle,"TDirectory"))) {
      titDir = (TDirectory*)titkey->ReadObj(); // read each object and take again the method title for each element of the list
      TString methodTitle;
      (*this).GetMethodTitle(methodTitle,titDir);
      TIter nextKey( titDir->GetListOfKeys() ); // loop and the list of keys
      while ((hkey = (*this).NextKey(nextKey,"TH1"))) { // take only the TH1 object type
        h = (TH1F*)hkey->ReadObj();
        TString hname = h->GetName();    // only the one which are called rejBvsS
        h -> SetName(Form("%s_%s",inputFiles_.at(iFile)->GetName(),h->GetName()));	
        if(plotType == 0){
         if (hname.Contains("effBvsS") && hname.BeginsWith("MVA_") && not hname.Contains("effBvsSLocal")) {
          if(size_t(color_index) <= vec_color.size()){
	    h->SetLineWidth(vec_linewidth[color_index]);
	    h->SetLineColor(vec_color[color_index]);
	    h->SetLineStyle(vec_linestyle[color_index]);
            cROC_->cd();
	    h->Draw("csame");
            cROCLog_->cd();
	    h->Draw("csame");
            hists.Add(h);
            color_index = color_index+1;
          }
          else{
	    h->SetLineWidth(vec_linewidth[color_index-vec_color.size()]);
	    h->SetLineColor(vec_color[color_index-vec_color.size()]);
	    h->SetLineStyle(vec_linestyle[color_index-vec_color.size()]);
            cROC_->cd();
	    h->Draw("csame");
            cROCLog_->cd();
	    h->Draw("csame");
	    hists.Add(h);
            color_index = color_index+1;
          }	
         }
	}
	else{
         if (hname.Contains("rejBvsS") && hname.BeginsWith("MVA_")) {
          if(size_t(color_index) <= vec_color.size()){
	    h->SetLineWidth(vec_linewidth[color_index]);
	    h->SetLineColor(vec_color[color_index]);
	    h->SetLineStyle(vec_linestyle[color_index]);
            cROC_->cd();
	    h->Draw("csame");
            cROCLog_->cd();
	    h->Draw("csame");
            hists.Add(h);
            color_index = color_index+1;
          }
          else{
	    h->SetLineWidth(vec_linewidth[color_index-vec_color.size()]);
	    h->SetLineColor(vec_color[color_index-vec_color.size()]);
	    h->SetLineStyle(vec_linestyle[color_index-vec_color.size()]);
            cROC_->cd();
	    h->Draw("csame");
            cROCLog_->cd();
	    h->Draw("csame");
	    hists.Add(h);
            color_index = color_index+1;
          }	
	 }
	}
      }
     }
   }

  /// Loop on the different histos                                                                                                                                                        
  while (hists.GetSize()) {
    TListIter hIt(&hists); // define an iterator                                                                                                                                          
    TH1F* hist(0);
    Double_t largestInt=-1;
    TH1F* histWithLargestInt(0);
    while ((hist = (TH1F*)hIt())!=0) {
      Double_t integral = hist->Integral(1,hist->GetNbinsX());
      if (integral>largestInt) {
        largestInt = integral;
        histWithLargestInt = hist;
      }
    }
    if (histWithLargestInt == 0) {
      std::cout << "ERROR - unknown hist \"histWithLargestInt\" --> serious problem in ROOT file" << std::endl;
      break;
    }    
      
    // set legend names 
    legROC_->AddEntry(histWithLargestInt,inputMethodName_.at(method_index).c_str(),"l");
    method_index ++ ;
    hists.Remove(histWithLargestInt);
  }
  }  
  

  cROC_->cd();  
  legROC_->Draw("same");
  cROC_->Update();

  cROCLog_->cd();  
  legROC_->Draw("same");
  cROCLog_->Update();
  
  return;

}



// method in order to plot correlation matrix betwenn input variables
void TMVAGlob::plotCorrelationMatrix(TFile* inputFile, const std::string & inputName, const std::string & outputPlotDirectory){
  
  std::string nameCorrelationSignal     = "CorrelationMatrixS";
  std::string nameCorrelationBackground = "CorrelationMatrixB";

  TList TrainingMethods;
  int res = (*this).GetListOfMethods(TrainingMethods);
  if(res == 0) std::cout<<" TMVAGlob::plotCorrelationMatrix --> no methods found "<<std::endl ;

  if(TString(TrainingMethods.At(0)->GetName()).Contains("Cuts")) return ;

  TH2F* hSignal     = (TH2F*) inputFile->Get(nameCorrelationSignal.c_str());
  TH2F* hBackground = (TH2F*) inputFile->Get(nameCorrelationBackground.c_str());

  hSignal->SetName    (std::string(Form("%s_%s",nameCorrelationSignal.c_str(),inputName.c_str())).c_str());
  hBackground->SetName(std::string(Form("%s_%s",nameCorrelationBackground.c_str(),inputName.c_str())).c_str());
  
  if(hSignal == 0 || hBackground == 0){ std::cerr<<" Null Pointer for correlation Matrix --> exit without plot "<<std::endl;  return ; }
  cCorrelationSignal_ = new TCanvas(std::string(Form("c%s_%s",nameCorrelationSignal.c_str(),inputName.c_str())).c_str(),Form("Correlation Matrix Signal"),180,52,550,550);
  float newMargin1 = 0.13;
  float newMargin2 = 0.15;
  float newMargin3 = 0.20;

  cCorrelationSignal_->SetGrid();
  cCorrelationSignal_->SetTicks();
  cCorrelationSignal_->SetLeftMargin(newMargin3);
  cCorrelationSignal_->SetBottomMargin(newMargin2);
  cCorrelationSignal_->SetRightMargin(newMargin1);
  cCorrelationSignal_->SetTopMargin(newMargin1);
  
  gStyle->SetPaintTextFormat("3g");

  hSignal->SetMarkerSize(1.5);
  hSignal->SetMarkerColor(0);
  hSignal->GetXaxis()->SetLabelSize(0.035);
  hSignal->GetYaxis()->SetLabelSize(0.035);
  hSignal->LabelsOption("d");
  hSignal->SetLabelOffset(0.011);// label offset on x axis                                                                                                                                
  
  hBackground->SetMarkerSize(1.5);
  hBackground->SetMarkerColor(0);
  hBackground->GetXaxis()->SetLabelSize(0.035);
  hBackground->GetYaxis()->SetLabelSize(0.035);
  hBackground->LabelsOption("d");
  hBackground->SetLabelOffset( 0.011 );// label offset on x axis                                                                                                                          
 
  // Plot correlation between signal      
  cCorrelationSignal_->cd();
  hSignal->Draw("colz");                                                                                                                                                   

  for( int iBinX = 0 ; iBinX<hSignal->GetNbinsX() ; iBinX++){
     TString LabelBinName = hSignal->GetXaxis()->GetBinLabel(iBinX+1);
     LabelBinName.ReplaceAll("[0]","");
     if(LabelBinName.Contains("mraw")) hSignal->GetXaxis()->SetBinLabel(iBinX+1,"M_{raw}");
     if(LabelBinName.Contains("mconst")) hSignal->GetXaxis()->SetBinLabel(iBinX+1,"M_{const}");
     if(LabelBinName.Contains("mclean")) hSignal->GetXaxis()->SetBinLabel(iBinX+1,"M_{clean}");
     if(LabelBinName.Contains("softdrop") and LabelBinName.Contains("beta00") ) hSignal->GetXaxis()->SetBinLabel(iBinX+1,"M_{SD}^{safe} #beta=0");
     if(LabelBinName.Contains("softdrop") and LabelBinName.Contains("beta10") ) hSignal->GetXaxis()->SetBinLabel(iBinX+1,"M_{SD}^{safe} #beta=1");
     if(LabelBinName.Contains("softdrop") and LabelBinName.Contains("beta20") ) hSignal->GetXaxis()->SetBinLabel(iBinX+1,"M_{SD}^{safe} #beta=2");
     if(LabelBinName.Contains("pruned")) hSignal->GetXaxis()->SetBinLabel(iBinX+1,"M_{pruned}^{safe}");
     if(LabelBinName.Contains("trim")) hSignal->GetXaxis()->SetBinLabel(iBinX+1,"M_{trimmed}^{safe}");
     if(LabelBinName.Contains("ecf_beta_20")) hSignal->GetXaxis()->SetBinLabel(iBinX+1,"C2(#beta=2.0)");
     if(LabelBinName.Contains("ecf_beta_15")) hSignal->GetXaxis()->SetBinLabel(iBinX+1,"C2(#beta=1.5)");
     if(LabelBinName.Contains("ecf_beta_10")) hSignal->GetXaxis()->SetBinLabel(iBinX+1,"C2(#beta=1.0)");
     if(LabelBinName.Contains("tau2/tau1")) hSignal->GetXaxis()->SetBinLabel(iBinX+1,"#tau_{2}/#tau_{1}");
     else if(LabelBinName.Contains("tau2") and not LabelBinName.Contains("tau1")) hSignal->GetXaxis()->SetBinLabel(iBinX+1,"#tau_{2}");
     else if(LabelBinName.Contains("tau1") and not LabelBinName.Contains("tau2")) hSignal->GetXaxis()->SetBinLabel(iBinX+1,"#tau_{1}");
     if(LabelBinName.Contains("QGLikelihood_pr_zcut_010_R_cut_050")) hSignal->GetXaxis()->SetBinLabel(iBinX+1,"QG Like");     
     if(LabelBinName.Contains("QGLikelihood_pr_sub2_zcut_010_R_cut_050") and LabelBinName.Contains("QGLikelihood_pr_sub1_zcut_010_R_cut_050")) hSignal->GetXaxis()->SetBinLabel(iBinX+1,"QG Like combo");     
     else if(LabelBinName.Contains("QGLikelihood_pr_sub1_zcut_010_R_cut_050")) hSignal->GetXaxis()->SetBinLabel(iBinX+1,"QG Like sub 1");     
     else if(LabelBinName.Contains("QGLikelihood_pr_sub2_zcut_010_R_cut_050")) hSignal->GetXaxis()->SetBinLabel(iBinX+1,"QG Like sub 2");     
     if(LabelBinName.Contains("Qjets")) hSignal->GetXaxis()->SetBinLabel(iBinX+1,"#Gamma_{Qjets}");     

     LabelBinName = hSignal->GetYaxis()->GetBinLabel(iBinX+1);
     LabelBinName.ReplaceAll("[0]","");
     if(LabelBinName.Contains("mraw")) hSignal->GetYaxis()->SetBinLabel(iBinX+1,"M_{raw}");
     if(LabelBinName.Contains("mconst")) hSignal->GetYaxis()->SetBinLabel(iBinX+1,"M_{const}");
     if(LabelBinName.Contains("mclean")) hSignal->GetYaxis()->SetBinLabel(iBinX+1,"M_{clean}");
     if(LabelBinName.Contains("softdrop") and LabelBinName.Contains("beta00") ) hSignal->GetYaxis()->SetBinLabel(iBinX+1,"M_{SD}^{safe} #beta=0");
     if(LabelBinName.Contains("softdrop") and LabelBinName.Contains("beta10") ) hSignal->GetYaxis()->SetBinLabel(iBinX+1,"M_{SD}^{safe} #beta=1");
     if(LabelBinName.Contains("softdrop") and LabelBinName.Contains("beta20") ) hSignal->GetYaxis()->SetBinLabel(iBinX+1,"M_{SD}^{safe} #beta=2");
     if(LabelBinName.Contains("pruned")) hSignal->GetYaxis()->SetBinLabel(iBinX+1,"M_{pruned}^{safe}");
     if(LabelBinName.Contains("trim")) hSignal->GetYaxis()->SetBinLabel(iBinX+1,"M_{trimmed}^{safe}");
     if(LabelBinName.Contains("ecf_beta_20")) hSignal->GetYaxis()->SetBinLabel(iBinX+1,"C2(#beta=2.0)");
     if(LabelBinName.Contains("ecf_beta_15")) hSignal->GetYaxis()->SetBinLabel(iBinX+1,"C2(#beta=1.5)");
     if(LabelBinName.Contains("ecf_beta_10")) hSignal->GetYaxis()->SetBinLabel(iBinX+1,"C2(#beta=1.0)");
     if(LabelBinName.Contains("tau2/tau1")) hSignal->GetYaxis()->SetBinLabel(iBinX+1,"#tau_{2}/#tau_{1}");
     else if(LabelBinName.Contains("tau2") and not LabelBinName.Contains("tau1")) hSignal->GetYaxis()->SetBinLabel(iBinX+1,"#tau_{2}");
     else if(LabelBinName.Contains("tau1") and not LabelBinName.Contains("tau2")) hSignal->GetYaxis()->SetBinLabel(iBinX+1,"#tau_{1}");
     if(LabelBinName.Contains("QGLikelihood_pr_zcut_010_R_cut_050")) hSignal->GetYaxis()->SetBinLabel(iBinX+1,"QG Like");     
     if(LabelBinName.Contains("QGLikelihood_pr_sub2_zcut_010_R_cut_050") and LabelBinName.Contains("QGLikelihood_pr_sub1_zcut_010_R_cut_050")) hSignal->GetYaxis()->SetBinLabel(iBinX+1,"QG Like combo");     
     else if(LabelBinName.Contains("QGLikelihood_pr_sub1_zcut_010_R_cut_050")) hSignal->GetYaxis()->SetBinLabel(iBinX+1,"QG Like sub 1");     
     else if(LabelBinName.Contains("QGLikelihood_pr_sub2_zcut_010_R_cut_050")) hSignal->GetYaxis()->SetBinLabel(iBinX+1,"QG Like sub 2");     
     if(LabelBinName.Contains("Qjets")) hSignal->GetYaxis()->SetBinLabel(iBinX+1,"#Gamma_{Qjets}");     
  }


  hSignal->Draw("textsame"); 

  // add comment                                                                                                                                                                       
  TText* text = new TText( 0.53, 0.88, "Linear correlation coefficients in %" );
  text->SetNDC();
  text->SetTextSize( 0.026 );
  text->AppendPad();
  cCorrelationSignal_->Update();

  // add comment                                                                                                                                                                       
  TText* text2 = new TText( 0.3, 0.92, "Correlation Matrix for Signal Events" );
  text2->SetNDC();
  text2->SetTextSize( 0.033 );
  text2->AppendPad();
  cCorrelationSignal_->Update();

  nameCorrelationSignal     = std::string(Form("%s_%s",nameCorrelationSignal.c_str(),inputName.c_str()));
  nameCorrelationBackground = std::string(Form("%s_%s",nameCorrelationBackground.c_str(),inputName.c_str()));

  (*this).PrintImage(cCorrelationSignal_,outputPlotDirectory+"/"+nameCorrelationSignal);
    
  // Background correlation
  cCorrelationBackground_ = new TCanvas(std::string(Form("c%s_%s",nameCorrelationBackground.c_str(),inputName.c_str())).c_str(),Form("Correlation Matrix Signal"),180,52,550,550);
  cCorrelationBackground_->SetGrid();
  cCorrelationBackground_->SetTicks();
  cCorrelationBackground_->SetLeftMargin(newMargin3);
  cCorrelationBackground_->SetBottomMargin(newMargin2);
  cCorrelationBackground_->SetRightMargin(newMargin1);
  cCorrelationBackground_->SetTopMargin(newMargin1);

  cCorrelationBackground_->cd();
  hBackground->Draw("colz");                                                                                                                             

  for( int iBinX = 0 ; iBinX<hBackground->GetNbinsX() ; iBinX++){
     TString LabelBinName = hBackground->GetXaxis()->GetBinLabel(iBinX+1);
     LabelBinName.ReplaceAll("[0]","");
     if(LabelBinName.Contains("mraw")) hBackground->GetXaxis()->SetBinLabel(iBinX+1,"M_{raw}");
     if(LabelBinName.Contains("mconst")) hBackground->GetXaxis()->SetBinLabel(iBinX+1,"M_{const}");
     if(LabelBinName.Contains("mclean")) hBackground->GetXaxis()->SetBinLabel(iBinX+1,"M_{clean}");
     if(LabelBinName.Contains("softdrop") and LabelBinName.Contains("beta00") ) hBackground->GetXaxis()->SetBinLabel(iBinX+1,"M_{SD}^{safe} #beta=0");
     if(LabelBinName.Contains("softdrop") and LabelBinName.Contains("beta10") ) hBackground->GetXaxis()->SetBinLabel(iBinX+1,"M_{SD}^{safe} #beta=1");
     if(LabelBinName.Contains("softdrop") and LabelBinName.Contains("beta20") ) hBackground->GetXaxis()->SetBinLabel(iBinX+1,"M_{SD}^{safe} #beta=2");
     if(LabelBinName.Contains("pruned")) hBackground->GetXaxis()->SetBinLabel(iBinX+1,"M_{pruned}^{safe}");
     if(LabelBinName.Contains("trim")) hBackground->GetXaxis()->SetBinLabel(iBinX+1,"M_{trimmed}^{safe}");
     if(LabelBinName.Contains("ecf_beta_20")) hBackground->GetXaxis()->SetBinLabel(iBinX+1,"C2(#beta=2.0)");
     if(LabelBinName.Contains("ecf_beta_15")) hBackground->GetXaxis()->SetBinLabel(iBinX+1,"C2(#beta=1.5)");
     if(LabelBinName.Contains("ecf_beta_10")) hBackground->GetXaxis()->SetBinLabel(iBinX+1,"C2(#beta=1.0)");
     if(LabelBinName.Contains("tau2/tau1")) hBackground->GetXaxis()->SetBinLabel(iBinX+1,"#tau_{2}/#tau_{1}");
     else if(LabelBinName.Contains("tau2") and not LabelBinName.Contains("tau1")) hBackground->GetXaxis()->SetBinLabel(iBinX+1,"#tau_{2}");
     else if(LabelBinName.Contains("tau1") and not LabelBinName.Contains("tau2")) hBackground->GetXaxis()->SetBinLabel(iBinX+1,"#tau_{1}");
     if(LabelBinName.Contains("QGLikelihood_pr_zcut_010_R_cut_050")) hBackground->GetXaxis()->SetBinLabel(iBinX+1,"QG Like");     
     if(LabelBinName.Contains("QGLikelihood_pr_sub2_zcut_010_R_cut_050") and LabelBinName.Contains("QGLikelihood_pr_sub1_zcut_010_R_cut_050")) hBackground->GetXaxis()->SetBinLabel(iBinX+1,"QG Like combo");     
     else if(LabelBinName.Contains("QGLikelihood_pr_sub1_zcut_010_R_cut_050")) hBackground->GetXaxis()->SetBinLabel(iBinX+1,"QG Like sub 1");     
     else if(LabelBinName.Contains("QGLikelihood_pr_sub2_zcut_010_R_cut_050")) hBackground->GetXaxis()->SetBinLabel(iBinX+1,"QG Like sub 2");     
     if(LabelBinName.Contains("Qjets")) hBackground->GetXaxis()->SetBinLabel(iBinX+1,"#Gamma_{Qjets}");     

     LabelBinName = hBackground->GetYaxis()->GetBinLabel(iBinX+1);
     LabelBinName.ReplaceAll("[0]","");
     if(LabelBinName.Contains("mraw")) hBackground->GetYaxis()->SetBinLabel(iBinX+1,"M_{raw}");
     if(LabelBinName.Contains("mconst")) hBackground->GetYaxis()->SetBinLabel(iBinX+1,"M_{const}");
     if(LabelBinName.Contains("mclean")) hBackground->GetYaxis()->SetBinLabel(iBinX+1,"M_{clean}");
     if(LabelBinName.Contains("softdrop") and LabelBinName.Contains("beta00") ) hBackground->GetYaxis()->SetBinLabel(iBinX+1,"M_{SD}^{safe} #beta=0");
     if(LabelBinName.Contains("softdrop") and LabelBinName.Contains("beta10") ) hBackground->GetYaxis()->SetBinLabel(iBinX+1,"M_{SD}^{safe} #beta=1");
     if(LabelBinName.Contains("softdrop") and LabelBinName.Contains("beta20") ) hBackground->GetYaxis()->SetBinLabel(iBinX+1,"M_{SD}^{safe} #beta=2");
     if(LabelBinName.Contains("pruned")) hBackground->GetYaxis()->SetBinLabel(iBinX+1,"M_{pruned}^{safe}");
     if(LabelBinName.Contains("trim")) hBackground->GetYaxis()->SetBinLabel(iBinX+1,"M_{trimmed}^{safe}");
     if(LabelBinName.Contains("ecf_beta_20")) hBackground->GetYaxis()->SetBinLabel(iBinX+1,"C2(#beta=2.0)");
     if(LabelBinName.Contains("ecf_beta_15")) hBackground->GetYaxis()->SetBinLabel(iBinX+1,"C2(#beta=1.5)");
     if(LabelBinName.Contains("ecf_beta_10")) hBackground->GetYaxis()->SetBinLabel(iBinX+1,"C2(#beta=1.0)");
     if(LabelBinName.Contains("tau2/tau1")) hBackground->GetYaxis()->SetBinLabel(iBinX+1,"#tau_{2}/#tau_{1}");
     else if(LabelBinName.Contains("tau2") and not LabelBinName.Contains("tau1")) hBackground->GetYaxis()->SetBinLabel(iBinX+1,"#tau_{2}");
     else if(LabelBinName.Contains("tau1") and not LabelBinName.Contains("tau2")) hBackground->GetYaxis()->SetBinLabel(iBinX+1,"#tau_{1}");
     if(LabelBinName.Contains("QGLikelihood_pr_zcut_010_R_cut_050")) hBackground->GetYaxis()->SetBinLabel(iBinX+1,"QG Like");     
     if(LabelBinName.Contains("QGLikelihood_pr_sub2_zcut_010_R_cut_050") and LabelBinName.Contains("QGLikelihood_pr_sub1_zcut_010_R_cut_050")) hBackground->GetYaxis()->SetBinLabel(iBinX+1,"QG Like combo");     
     else if(LabelBinName.Contains("QGLikelihood_pr_sub1_zcut_010_R_cut_050")) hBackground->GetYaxis()->SetBinLabel(iBinX+1,"QG Like sub 1");     
     else if(LabelBinName.Contains("QGLikelihood_pr_sub2_zcut_010_R_cut_050")) hBackground->GetYaxis()->SetBinLabel(iBinX+1,"QG Like sub 2");     
     if(LabelBinName.Contains("Qjets")) hBackground->GetYaxis()->SetBinLabel(iBinX+1,"#Gamma_{Qjets}");     
  }


  hBackground->Draw("textsame");          
  // add comment                                                                                                                                                                       
  text->AppendPad();

  // add comment                                                                                                                                                                       
  TText* text3 = new TText( 0.3, 0.92, "Correlation Matrix for Background Events" );
  text3->SetNDC();
  text3->SetTextSize( 0.033 );
  text3->AppendPad();

  cCorrelationBackground_->Update();
  (*this).PrintImage(cCorrelationBackground_,outputPlotDirectory+"/"+nameCorrelationBackground);

  if(hSignal!=0)     delete hSignal;
  if(hBackground!=0) delete hBackground;
  if(text!=0)        delete text;

  return ;

}

// Methods to print image given canvas and a output name
void TMVAGlob::PrintImage(TCanvas* c, const std::string & fname){

  TString pngName = fname+".png";
  TString pdfName = fname+".pdf";
  TString rootName = fname+".root";

  c->Print(pdfName);
  c->Print(pngName);
  c->Print(rootName);

  return ;

}

// plot MVA output, probability and overtraining 
void TMVAGlob::plotMVAs(TFile* inputFile,  const std::string & inputName, HistType htype, const std::string & outputPlotDirectory){

  TIter next(inputFile->GetListOfKeys());
  TKey *key(0);

  while ((key = (TKey*)next())) {
    if (!TString(key->GetName()).BeginsWith("Method_")) continue;
    if (!gROOT->GetClass(key->GetClassName())->InheritsFrom("TDirectory")) continue;

    TString methodName;
    (*this).GetMethodName(methodName,key); // take the method name for each key inside the input file

    TDirectory* mDir = (TDirectory*)key->ReadObj(); // take the list of jey of the selected object

    TIter keyIt(mDir->GetListOfKeys());
    TKey *titkey;

    while ((titkey = (TKey*)keyIt())) { // loop on the second list of keys

      if (!gROOT->GetClass(titkey->GetClassName())->InheritsFrom("TDirectory")) continue; // if is not a TDirectory object skipe

      TDirectory *titDir = (TDirectory *)titkey->ReadObj();
      TString methodTitle;
      (*this).GetMethodTitle(methodTitle,titDir); // get the tutke of the method (another directory found)
      std::cout << "--- Found directory for method: " << methodName << "::" << methodTitle << std::flush;

      TString hname = "MVA_" + methodTitle;
      if      (htype == ProbaType  ) hname += "_Proba";
      else if (htype == RarityType ) hname += "_Rarity";

      // take output distribution for signal and baclground
      histoSignal_     = dynamic_cast<TH1*>(titDir->Get(hname+"_S" ));
      histoBackground_ = dynamic_cast<TH1*>(titDir->Get(hname+"_B"));

      if (histoSignal_==0 || histoBackground_==0) {
	if (htype == MVAType)
	  std::cout << ":\t mva distribution not available (this is normal for Cut classifier)" << std::endl;
	else if(htype == ProbaType)
	  std::cout << ":\t probability distribution not available" << std::endl;
	else if(htype == RarityType)
	  std::cout << ":\t rarity distribution not available" << std::endl;
	else if(htype == CompareType)
	  std::cout << ":\t overtraining check not available" << std::endl;
	else std::cout << std::endl;
	continue;
      }

      if((methodTitle).Contains("Likelihood"))
         methodTitle = Form("Likelhood");
      else if ((methodTitle).Contains("LD"))
         methodTitle = Form("Linear Discriminant");
      else if ((methodTitle).Contains("BDT") && !(methodTitle).Contains("BDTG") && !(methodTitle).Contains("BDTF"))
         methodTitle = Form("Boosted Decision Tree");
      else if ((methodTitle).Contains("MLP"))
         methodTitle = Form("Multi-Layer Perceptron (MLP)");
      else if ((methodTitle).Contains("PDEFoam")) 
         methodTitle = Form("PDEFoam");
      else if ((methodTitle).Contains("SVN")) 
         methodTitle = Form("Supported Vector Machine (SVM)");
      else if ((methodTitle).Contains("Fisher")) 
         methodTitle = Form("Fisher Discriminant");
        
      histoSignal_->SetTitle( Form("Response for classifier: %s", methodTitle.Data()) );
      if(htype == ProbaType)
	histoSignal_->SetTitle( Form("Probability for classifier: %s", methodTitle.Data()) );
      else if (htype == RarityType)
        histoSignal_->SetTitle( Form("Rarity for classifier: %s", methodTitle.Data()) );
      else if (htype == CompareType)
	histoSignal_->SetTitle( Form("Overtraining check for classifier: %s", methodTitle.Data()) );                                                                                


      TString cCanvasTitle = ((htype == MVAType)     ? Form("Response %s",methodTitle.Data()) :
			      (htype == ProbaType)   ? Form("Probability %s",methodTitle.Data()) :
			      (htype == CompareType) ? Form("Comparison %s",methodTitle.Data()) :
			      Form("Rarity %s",methodTitle.Data()));


      cMVAs_ = new TCanvas( Form("cMVAs_%d",(*this).mvas_index), cCanvasTitle,180,100,550,550);
      cMVAs_->cd();
      cMVAs_->SetTicks();
      cMVAs_->SetFillColor(0);
      cMVAs_->SetBorderMode(0);
      cMVAs_->SetBorderSize(2);
  
      cMVAs_->SetTickx(1);
      cMVAs_->SetTicky(1);
      cMVAs_->SetRightMargin(0.05);
      cMVAs_->SetBottomMargin(0.12);
      cMVAs_->SetFrameBorderMode(0);

      // set the histogram style                                                                                                                                                        
      histoSignal_->SetFillColor(4);
      histoSignal_->SetLineColor(4);
      histoSignal_->SetLineWidth(2);
      histoSignal_->SetFillStyle(1001);

      if(htype == CompareType) histoSignal_->SetFillStyle(3001);
         
      histoBackground_->SetFillColor(2);
      histoBackground_->SetLineColor(2);
      histoBackground_->SetLineWidth(2);
      histoBackground_->SetFillStyle(3005);
    
      // normalise both signal and background                                                                                                                                           
      (*this).NormalizeHists((TH1F*)histoSignal_,(TH1F*)histoBackground_);

      // frame limits                                                                                                                
      Float_t nrms = 10;
      std::cout << "--- Mean and RMS (S): " << histoSignal_->GetMean() << ", " << histoSignal_->GetRMS() << std::endl;
      std::cout << "--- Mean and RMS (B): " << histoBackground_->GetMean() << ", " << histoBackground_->GetRMS() << std::endl;

      Float_t xmin = TMath::Max(TMath::Min(histoSignal_->GetMean()-nrms*histoSignal_->GetRMS(),histoBackground_->GetMean()-nrms*histoBackground_->GetRMS()),histoSignal_->GetXaxis()->GetXmin());
      Float_t xmax = TMath::Min(TMath::Max(histoSignal_->GetMean()+nrms*histoSignal_->GetRMS(),histoBackground_->GetMean()+nrms*histoBackground_->GetRMS()),histoSignal_->GetXaxis()->GetXmax());
      Float_t ymin = 0;
      Float_t maxMult = (htype == CompareType) ? 1.3 : 1.2;
      Float_t ymax = TMath::Max(histoSignal_->GetMaximum(),histoBackground_->GetMaximum())*maxMult;

      // build a frame                                                                                                                                                                  
      Int_t nb = 500;
      TString hFrameName(TString("frame") + methodTitle);
      TObject *o = gROOT->FindObject(hFrameName);
      if(o) delete o;

      TH2F* frame = new TH2F( hFrameName, histoSignal_->GetTitle(),nb, xmin, xmax, nb, ymin, ymax );
      frame->GetXaxis()->SetTitle( methodTitle + ((htype == MVAType || htype == CompareType) ? " response" : "") );

      if      (htype == ProbaType  ) frame->GetXaxis()->SetTitle( "Signal probability" );
      else if (htype == RarityType ) frame->GetXaxis()->SetTitle( "Signal rarity" );
      frame->GetYaxis()->SetTitle("(1/N) dN/dx");

      frame->SetLineWidth(2);
      frame->SetMarkerStyle(21);
      frame->SetMarkerSize(0.3);
      frame->GetXaxis()->SetNdivisions(405);
      frame->GetYaxis()->SetNdivisions(405);
      frame->GetXaxis()->SetLabelOffset(0.012);
      frame->GetXaxis()->SetLabelSize(0.03);
      frame->GetXaxis()->SetTitleSize(0.035);
      frame->GetXaxis()->SetTitleOffset(1.05);
      frame->GetYaxis()->SetLabelOffset(0.012);
      frame->GetYaxis()->SetLabelSize(0.03);
      frame->GetYaxis()->SetTitleSize(0.035);
      frame->GetYaxis()->SetTitleOffset(1.10);
      frame->Draw("");
   
      // Draw legend                                                                                                                                                                    
      TLegend *legend= new TLegend( cMVAs_->GetLeftMargin()+0.03, 1-cMVAs_->GetTopMargin()-0.15,
				    cMVAs_->GetLeftMargin()+(htype == CompareType ? 0.40 : 0.3)+0.15, 1-cMVAs_->GetTopMargin()-0.03);
      legend->SetFillStyle(0);
      legend->SetFillColor(0);
      legend->SetTextSize (0.042);
      if(htype == CompareType) legend->SetTextSize (0.032); 
      legend->SetBorderSize(0);
      legend->SetTextFont(42);
      legend->SetLineColor(1);
      legend->SetLineStyle(1);
      legend->SetLineWidth(2);

      legend->AddEntry(histoSignal_,TString("Signal")     + ((htype == CompareType) ? " (testing)" : ""), "F");
      legend->AddEntry(histoBackground_,TString("Background") + ((htype == CompareType) ? " (testing)" : ""), "F");

      legend->SetMargin( (htype == CompareType ? 0.2 : 0.3) );
      legend->Draw("same");

      histoSignal_->Draw("samehist");
      histoBackground_->Draw("samehist");

      TLegend *legend2= new TLegend( 1-cMVAs_->GetRightMargin()-0.4, 1-cMVAs_->GetTopMargin()-0.15,
	        		     1-cMVAs_->GetRightMargin()-0.02,1-cMVAs_->GetTopMargin()-0.03);

      if (htype == CompareType) {
	// if overtraining check, load additional histograms                                                                                                                           
	TH1* sigOv = 0;
	TH1* bgdOv = 0;

	TString ovname = hname += "_Train";
	sigOv = dynamic_cast<TH1*>(titDir->Get( ovname + "_S" ));
	bgdOv = dynamic_cast<TH1*>(titDir->Get( ovname + "_B" ));

	// normalise both signal and background                                                                                                                                        

	if (sigOv == 0 || bgdOv == 0) std::cout << "+++ overtraining check histograms do not exist" << std::endl; return;
	
	
	(*this).NormalizeHists((TH1F*)sigOv,(TH1F*)bgdOv );

	Int_t col = histoSignal_->GetLineColor();
	sigOv->SetMarkerColor(col);
	sigOv->SetMarkerSize(0.7);
	sigOv->SetMarkerStyle(20);
	sigOv->SetLineWidth(1);
	sigOv->SetLineColor(col);
	sigOv->Draw("e1same");

	col = histoBackground_->GetLineColor();
	bgdOv->SetMarkerColor(col);
	bgdOv->SetMarkerSize(0.7);
	bgdOv->SetMarkerStyle(20);
	bgdOv->SetLineWidth(1);
	bgdOv->SetLineColor(col);
	bgdOv->Draw("e1same");

	ymax = TMath::Max( ymax, float(TMath::Max(sigOv->GetMaximum(),bgdOv->GetMaximum())*maxMult));
	frame->GetYaxis()->SetLimits(0,ymax);

        legend2->SetFillStyle(0);
        legend2->SetFillColor(0);
        legend2->SetTextSize (0.03);
        legend2->SetBorderSize(0);
        legend2->SetTextFont(42);
        legend2->SetLineColor(1);
        legend2->SetLineStyle(1);
        legend2->SetLineWidth(2);
 	legend2->AddEntry(sigOv,"Signal (training)","P");
 	legend2->AddEntry(bgdOv,"Background (training)","P");
 	legend2->SetMargin( 0.1 );
 	legend2->Draw("same");

	// perform K-S test                                                                                                                                                            
	std::cout << "--- Perform Kolmogorov-Smirnov tests" << std::endl;
	Double_t kolS = histoSignal_->KolmogorovTest(sigOv);
	Double_t kolB = histoBackground_->KolmogorovTest(bgdOv);
	std::cout << "--- Goodness of signal (background) consistency: " << kolS << " (" << kolB << ")" << std::endl;

	TString probatext = Form("Kolmogorov-Smirnov test:");
        TString probatext2 = Form("signal (background) probability = %5.3g (%5.3g)", kolS, kolB );
	TLatex* tt  = new TLatex(0.25,0.72,probatext);
	TLatex* tt2 = new TLatex(0.25,0.69,probatext2);
	tt->SetNDC(); 
        tt->SetTextSize(0.025); 
        tt->AppendPad();
	tt2->SetNDC(); 
        tt2->SetTextSize(0.025); 
        tt2->AppendPad();

	// chi2 test                                                                                                                                            
	std::cout << "--- Perform Chi2 tests" << std::endl;
	Double_t chi2S = histoSignal_->Chi2Test(sigOv,"WW CHI2/NDF" );
	Double_t chi2B = histoBackground_->Chi2Test(bgdOv,"WW CHI2/NDF" );
	std::cout << "--- Goodness of signal (background) consistency: " << chi2S << " (" << chi2B << ")" << std::endl;

	TString probatext3 = Form("#Chi^{2} test:");
        TString probatext4 = Form("signal (background) #Chi^{2}/ndf = %5.3g (%5.3g)", chi2S, chi2B );
	TLatex* tt3  = new TLatex(0.25,0.65,probatext3);
	TLatex* tt4  = new TLatex(0.25,0.62,probatext4);
	tt3->SetNDC(); 
        tt3->SetTextSize(0.025); 
        tt3->AppendPad();
	tt4->SetNDC(); 
        tt4->SetTextSize(0.025); 
        tt4->AppendPad();

      }

      // redraw axes                                                                                                                                                                    
      frame->Draw("sameaxis");

      // text for overflows                                                                                                                                                             
      Int_t    nbin = histoSignal_->GetNbinsX();
      Double_t dxu  = histoSignal_->GetBinWidth(0);
      Double_t dxo  = histoSignal_->GetBinWidth(nbin+1);
      TString uoflow = Form( "U/O-flow (S,B): (%.1f, %.1f)%% / (%.1f, %.1f)%%",
			     histoSignal_->GetBinContent(0)*dxu*100, histoBackground_->GetBinContent(0)*dxu*100,
			     histoSignal_->GetBinContent(nbin+1)*dxo*100, histoBackground_->GetBinContent(nbin+1)*dxo*100 );

      TText* t = new TText( 0.975, 0.115, uoflow );
      t->SetNDC();
      t->SetTextSize( 0.030 );
      t->SetTextAngle( 90 );
      t->AppendPad();


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
      tex = new TLatex(0.2324,0.92,"Simulation Preliminary");
      tex->SetNDC();
      tex->SetTextFont(52);
      tex->SetTextSize(0.0304);
      tex->SetLineWidth(2);
      tex->Draw();
      cMVAs_->Update();

      methodTitle.ReplaceAll(" ","_");
      methodTitle.ReplaceAll("(","_");
      methodTitle.ReplaceAll(")","_");
      methodTitle.ReplaceAll("[","_");
      methodTitle.ReplaceAll("]","_");
      methodTitle.ReplaceAll("{","_");
      methodTitle.ReplaceAll("}","_");

      if      (htype == MVAType)     (*this).PrintImage(cMVAs_, std::string(Form("%s/mva_output_%s",outputPlotDirectory.c_str(),inputName.c_str())));
      else if (htype == ProbaType)   (*this).PrintImage(cMVAs_, std::string(Form("%s/probability_%s",outputPlotDirectory.c_str(),inputName.c_str())));
      else if (htype == CompareType) (*this).PrintImage(cMVAs_, std::string(Form("%s/overtraining_%s",outputPlotDirectory.c_str(),inputName.c_str())));
      else                           (*this).PrintImage(cMVAs_, std::string(Form("%s/rarity_%s",outputPlotDirectory.c_str(),inputName.c_str())));
    
      mvas_index = mvas_index  +1 ;

      if(frame!=0)   delete frame;
      if(legend!=0)  delete legend;
      if(t!=0)       delete t;
      if(legend2!=0) delete legend2 ;
	       
    }

  }	       

  return ;

}

// methods to read efficiency histograms for each method in order to do significance plot
void TMVAGlob::ReadHistograms (TFile* inputFile){
  
  // if fInfoList is not null cancel it
  if (fInfoList_!=NULL) { fInfoList_->clear();
                          delete fInfoList_;
                          fInfoList_ = NULL ;
  }
  
  // search for the right histograms in full list of keys                                                                                                                                 
  fInfoList_ = new std::vector<significanceBox*> ;
  TIter next(inputFile->GetListOfKeys());
  TKey *key(0);

  while( (key = (TKey*)next()) ) {

    if(!TString(key->GetName()).BeginsWith("Method_")) continue;
    if(!gROOT->GetClass(key->GetClassName())->InheritsFrom("TDirectory") ) continue;
    TDirectory* mDir = (TDirectory*)key->ReadObj();
    
    TIter keyIt(mDir->GetListOfKeys());
    TKey *titkey;
    int maxLenTitle = 0 ;

    while((titkey = (TKey*)keyIt())) {
      significanceBox* significance_ = new significanceBox;
      if(!gROOT->GetClass(titkey->GetClassName())->InheritsFrom("TDirectory") ) continue;
      TDirectory* titDir = (TDirectory *)titkey->ReadObj();      
      (*this).GetMethodName(significance_->methodName_,key); // get method name 
      (*this).GetMethodTitle(significance_->methodTitle_,titDir); // get method title
      if (significance_->methodTitle_.Length() > maxLenTitle) maxLenTitle = significance_->methodTitle_.Length();
      TString hname = "MVA_" + significance_->methodTitle_;

      std::cout << "--- Classifier: " << significance_->methodTitle_ << std::endl;
      
      if(histoBackground_ == NULL ) significance_->Background_ = dynamic_cast<TH1*>(titDir->Get( hname + "_B" ));
      else significance_->Background_ = histoBackground_ ;

      if(histoSignal_ == NULL)      significance_->Signal_ = dynamic_cast<TH1*>(titDir->Get( hname + "_S" ));
      else significance_->Signal_ = histoSignal_ ;

      if(effBackground_ == NULL ) significance_->efficiencyBackground_ = dynamic_cast<TH1*>(titDir->Get( hname + "_effB" ));
      else significance_->efficiencyBackground_ = effBackground_ ;

      if(effSignal_ ==NULL)  significance_->efficiencySignal_  = dynamic_cast<TH1*>(titDir->Get( hname + "_effS" ));
      else significance_->efficiencySignal_ = effSignal_ ;
     
      if (significance_->efficiencyBackground_ == NULL || significance_->efficiencySignal_ == NULL){ delete significance_; continue; }
      fInfoList_->push_back(significance_);
    }
  }
  
 return;

}

// get a different string for the formula in order to be used in a TFormula object
TString TMVAGlob::GetFormula(){

  TString f = fFormula_;

  f.ReplaceAll("epsilonS","x");
  f.ReplaceAll("S","x");
  f.ReplaceAll("epsilonB","y");
  f.ReplaceAll("B","y"); 

  return f;
}

// get the latex string for the formula
TString TMVAGlob::GetLatexFormula(){

  TString f = fFormula_;
  if (signalType_ == true) 
   f.ReplaceAll("S","#epsilon_{S}");
  if (backgroundType_ == true) 
   f.ReplaceAll("B","#epsilon_{B}");   

  f.ReplaceAll("(","{");
  f.ReplaceAll(")","}");
  f.ReplaceAll("sqrt","#sqrt");

  return f;
}


// final method to do significance plot
void TMVAGlob::plotSignificance (TFile* inputFile, const std::string & inputName, SignificanceType stype, 
                                 const double & numberSignalEvents, const double & numberBackgroundEvents,
  		                 const bool & UseSignalEfficiency, const bool & UseBackgroundEfficiency, const std::string & outputPlotDirectory){


  if(inputFile == 0) {std::cerr<<" empty file --> exit "<<std::endl; return ; }

  // read the histo storing the signal and background effciencies and distribution
  (*this).ReadHistograms(inputFile); 
  
  // set the formula for significance inside the available choices
  if(stype == 0)       (*this).SetFormula("S/B"); 
  else if(stype == 1)  (*this).SetFormula("S/sqrt(B)");  
  else if(stype == 2)  (*this).SetFormula("S/sqrt(S+B)");  
  else if(stype == 3)  (*this).SetFormula("2*(sqrt(S+B)-sqrt(B))");  
  else { std::cerr<<" Not known formula --> exit "<<std::endl; return; }

  // set if efficiency have to be used alone in the significance or multiply them for the yields
  (*this).SetSignalType(UseSignalEfficiency);
  (*this).SetBackgroundType(UseBackgroundEfficiency);

  // TFormula for the significance
  TFormula significanceFormula("significanceFormula",(*this).GetFormula());
  
  std::vector<significanceBox*>::iterator itList =  (*this).fInfoList_->begin();
  TString cname = "Classifier";
  int maxLenTitle = 0 ;

  if (cname.Length() >  maxLenTitle)  maxLenTitle = cname.Length();

  TString str = Form( "%*s   (  signal, backgr.)  Optimal-cut  %s      NSig      NBkg   EffSig   EffBkg", maxLenTitle, cname.Data(), GetFormulaString().Data() );
  std::cout << "--- " << std::setfill('=') << std::setw(str.Length()) << "" << std::setfill(' ') << std::endl;
  std::cout << "--- " << str << std::endl;
  std::cout << "--- " << std::setfill('-') << std::setw(str.Length()) << "" << std::setfill(' ') << std::endl;


  // loop on the histo list   
  for ( ;itList!=(*this).fInfoList_->end(); ++itList) {
    if(signalType_)
      (*itList)->significance_ = new TH1F(Form("significance_%s_%s_stype%d",(*itList)->methodTitle_.Data(),inputName.c_str(),stype),"",(*itList)->efficiencySignal_->GetNbinsX(),(*itList)->efficiencySignal_->GetBinLowEdge(1),(*itList)->efficiencySignal_->GetBinLowEdge((*itList)->efficiencySignal_->GetNbinsX()+1));
    else
      (*itList)->significance_ = new TH1F(Form("significance_eff_%s_%s_stype%d",(*itList)->methodTitle_.Data(),inputName.c_str(),stype),"",(*itList)->efficiencySignal_->GetNbinsX(),(*itList)->efficiencySignal_->GetBinLowEdge(1),(*itList)->efficiencySignal_->GetBinLowEdge((*itList)->efficiencySignal_->GetNbinsX()+1));

    for (Int_t iBin=1; iBin<=(*itList)->efficiencySignal_->GetNbinsX(); iBin++) {

      Float_t S = 0;
      if(!signalType_) S = (*itList)->efficiencySignal_->GetBinContent(iBin) * numberSignalEvents;
      else S = (*itList)->efficiencySignal_->GetBinContent(iBin) ;      

      Float_t B = 0;
      if(!backgroundType_) B = (*itList)->efficiencyBackground_->GetBinContent(iBin) * numberBackgroundEvents;
      else B = (*itList)->efficiencyBackground_->GetBinContent(iBin) ;

      Double_t significance = 0. ;
      Double_t sigErr = 0.;

      // evaluate the significance
      if(stype == 0 && B != 0 )        significance = significanceFormula.Eval(S,B);
      else if(stype == 1 && B != 0 )   significance = significanceFormula.Eval(S,B);
      else if(stype == 2 && S+B != 0 ) significance = significanceFormula.Eval(S,B);
      else significance = significanceFormula.Eval(S,B);

      if ((*this).GetFormulaString() == "S/B")              sigErr = sqrt(S/(B*B)+S*S/(B*B*B));	
      else if ((*this).GetFormulaString() == "S/sqrt(B)")   sigErr = significance * sqrt( 1./S + 1./(2.*B));	
      else if ((*this).GetFormulaString() == "S/sqrt(S+B)") sigErr = sqrt(S*(TMath::Power(1-0.5/sqrt(S+B),2))*1/(S+B)+B*S*S/(4*(S+B)));	
      else if ((*this).GetFormulaString() == "2*(sqrt(S+B)-sqrt(B))") sigErr = sqrt(S*TMath::Power(1/sqrt(S+B),2)+B*TMath::Power(1/sqrt(S+B)-1/sqrt(B),2));	

      // set value and error                  
      (*itList)->significance_->SetBinContent(iBin,significance);
      (*itList)->significance_->SetBinError(iBin,sigErr);
    }
    
  }

  
  /// Plot the results  
  std::vector<significanceBox*>::iterator itInfoList = fInfoList_->begin();
  for ( ; itInfoList!=(*this).fInfoList_->end(); ++itInfoList ){
  // create new canvas                                                                                                                                                                 
   if(signalType_)
     cSignificance_ = new TCanvas( Form("cSignificance_%s_%s_type_%d",(*itInfoList)->methodTitle_.Data(),inputName.c_str(),stype),Form("Efficiencies Classifier : %s",(*itInfoList)->methodTitle_.Data()),180,52,550,550);
   else
     cSignificance_ = new TCanvas( Form("cSignificance_eff_%s_%s_type_%d",(*itInfoList)->methodTitle_.Data(),inputName.c_str(),stype),Form("Efficiencies Classifier : %s",(*itInfoList)->methodTitle_.Data()),180,52,550,550);

   cSignificance_->cd();
   cSignificance_->SetTickx(1);
   cSignificance_->SetTicky(0);
   cSignificance_->SetFillColor(0);
   cSignificance_->SetBorderMode(0);
   cSignificance_->SetBorderSize(2);
  
   cSignificance_->SetLeftMargin(0.1);
   cSignificance_->SetRightMargin(0.12);
   cSignificance_->SetBottomMargin(0.12);
   cSignificance_->SetFrameBorderMode(0);

   int bin  = 0 ;   
   for(int iBin = 0; iBin < (*itInfoList)->efficiencyBackground_->GetNbinsX(); iBin++){
     if((*itInfoList)->efficiencyBackground_->GetBinContent(iBin+1)>0.001) continue; // break when the efficiency on the bkg is less than 3% and save the bin
     else{ bin = iBin ; break; }
   }

   // Set X axis limit for significance plots
   (*itInfoList)->efficiencySignal_->SetTitle("Efficiencies and Optimal Cut");
   if ((*itInfoList)->methodTitle_.Contains("Cuts")) {
      (*itInfoList)->efficiencyBackground_->GetXaxis()->SetTitle( "Signal Efficiency (#epsilon_{sig})" );    
      (*itInfoList)->efficiencyBackground_->GetXaxis()->SetRangeUser(0.1,1.); // no less than 10% of signal efficiency for rectangular cut plots 
      (*itInfoList)->significance_->GetXaxis()->SetRangeUser(0.1,1.);
   }
   else if ((*itInfoList)->methodTitle_.Contains("Likelihood")){
      (*itInfoList)->efficiencyBackground_->GetXaxis()->SetTitle( "Likelihood Output" );
      (*itInfoList)->efficiencyBackground_->GetXaxis()->SetRangeUser((*itInfoList)->efficiencyBackground_->GetBinLowEdge(1),(*itInfoList)->efficiencyBackground_->GetBinCenter(bin)-(*itInfoList)->efficiencyBackground_->GetBinWidth(bin));
      (*itInfoList)->significance_->GetXaxis()->SetRangeUser((*itInfoList)->efficiencyBackground_->GetBinLowEdge(1),(*itInfoList)->efficiencyBackground_->GetBinCenter(bin)-(*itInfoList)->efficiencyBackground_->GetBinWidth(bin));
   }
   else if ((*itInfoList)->methodTitle_.Contains("LD")){
      (*itInfoList)->efficiencyBackground_->GetXaxis()->SetTitle( "Linear Discriminant Output" );    
      (*itInfoList)->efficiencyBackground_->GetXaxis()->SetRangeUser((*itInfoList)->efficiencyBackground_->GetBinLowEdge(1),(*itInfoList)->efficiencyBackground_->GetBinCenter(bin)-(*itInfoList)->efficiencyBackground_->GetBinWidth(bin));
      (*itInfoList)->significance_->GetXaxis()->SetRangeUser((*itInfoList)->efficiencyBackground_->GetBinLowEdge(1),(*itInfoList)->efficiencyBackground_->GetBinCenter(bin)-(*itInfoList)->efficiencyBackground_->GetBinWidth(bin));
   }
   else if ((*itInfoList)->methodTitle_.Contains("BDT") && !(*itInfoList)->methodTitle_.Contains("BDTG")){
      (*itInfoList)->efficiencyBackground_->GetXaxis()->SetTitle( "BDT Output" );    
      (*itInfoList)->efficiencyBackground_->GetXaxis()->SetRangeUser((*itInfoList)->efficiencyBackground_->GetBinLowEdge(1),(*itInfoList)->efficiencyBackground_->GetBinCenter(bin)-(*itInfoList)->efficiencyBackground_->GetBinWidth(bin));
      (*itInfoList)->significance_->GetXaxis()->SetRangeUser((*itInfoList)->efficiencyBackground_->GetBinLowEdge(1),(*itInfoList)->efficiencyBackground_->GetBinCenter(bin)-(*itInfoList)->efficiencyBackground_->GetBinWidth(bin));
   }
   else if ((*itInfoList)->methodTitle_.Contains("BDTG")){
      (*itInfoList)->efficiencyBackground_->GetXaxis()->SetTitle( "Gradient BDT Output" );    
      (*itInfoList)->efficiencyBackground_->GetXaxis()->SetRangeUser((*itInfoList)->efficiencyBackground_->GetBinLowEdge(1),(*itInfoList)->efficiencyBackground_->GetBinCenter(bin)-(*itInfoList)->efficiencyBackground_->GetBinWidth(bin));
      (*itInfoList)->significance_->GetXaxis()->SetRangeUser((*itInfoList)->efficiencyBackground_->GetBinLowEdge(1),(*itInfoList)->efficiencyBackground_->GetBinCenter(bin)-(*itInfoList)->efficiencyBackground_->GetBinWidth(bin));
   }
   else if ((*itInfoList)->methodTitle_.Contains("MLP")){
      (*itInfoList)->efficiencyBackground_->GetXaxis()->SetTitle( "MLP Output" );    
      (*itInfoList)->efficiencyBackground_->GetXaxis()->SetRangeUser((*itInfoList)->efficiencyBackground_->GetBinLowEdge(1),(*itInfoList)->efficiencyBackground_->GetBinCenter(bin)-(*itInfoList)->efficiencyBackground_->GetBinWidth(bin));
      (*itInfoList)->significance_->GetXaxis()->SetRangeUser((*itInfoList)->efficiencyBackground_->GetBinLowEdge(1),(*itInfoList)->efficiencyBackground_->GetBinCenter(bin)-(*itInfoList)->efficiencyBackground_->GetBinWidth(bin));
   }
   else if ((*itInfoList)->methodTitle_.Contains("PDEFoam")){
      (*itInfoList)->efficiencyBackground_->GetXaxis()->SetTitle( "PDEFoam Output" );    
      (*itInfoList)->efficiencyBackground_->GetXaxis()->SetRangeUser((*itInfoList)->efficiencyBackground_->GetBinLowEdge(1),(*itInfoList)->efficiencyBackground_->GetBinCenter(bin)-(*itInfoList)->efficiencyBackground_->GetBinWidth(bin));
      (*itInfoList)->significance_->GetXaxis()->SetRangeUser((*itInfoList)->efficiencyBackground_->GetBinLowEdge(1),(*itInfoList)->efficiencyBackground_->GetBinCenter(bin)-(*itInfoList)->efficiencyBackground_->GetBinWidth(bin));
   }
   else if ((*itInfoList)->methodTitle_.Contains("SVM")){
      (*itInfoList)->efficiencyBackground_->GetXaxis()->SetTitle( "SVM Output" );    
      (*itInfoList)->efficiencyBackground_->GetXaxis()->SetRangeUser((*itInfoList)->efficiencyBackground_->GetBinLowEdge(1),(*itInfoList)->efficiencyBackground_->GetBinCenter(bin)-(*itInfoList)->efficiencyBackground_->GetBinWidth(bin));
      (*itInfoList)->significance_->GetXaxis()->SetRangeUser((*itInfoList)->efficiencyBackground_->GetBinLowEdge(1),(*itInfoList)->efficiencyBackground_->GetBinCenter(bin)-(*itInfoList)->efficiencyBackground_->GetBinWidth(bin));
   }
   else if ((*itInfoList)->methodTitle_.Contains("Fisher")){
      (*itInfoList)->efficiencyBackground_->GetXaxis()->SetTitle( "Fisher Output" );    
      (*itInfoList)->efficiencyBackground_->GetXaxis()->SetRangeUser((*itInfoList)->efficiencyBackground_->GetBinLowEdge(1),(*itInfoList)->efficiencyBackground_->GetBinCenter(bin)-(*itInfoList)->efficiencyBackground_->GetBinWidth(bin));
      (*itInfoList)->significance_->GetXaxis()->SetRangeUser((*itInfoList)->efficiencyBackground_->GetBinLowEdge(1),(*itInfoList)->efficiencyBackground_->GetBinCenter(bin)-(*itInfoList)->efficiencyBackground_->GetBinWidth(bin));
   }
   else{
       (*itInfoList)->efficiencyBackground_->GetXaxis()->SetTitle( (*itInfoList)->methodTitle_ + " output" );
       (*itInfoList)->efficiencyBackground_->GetXaxis()->SetRangeUser((*itInfoList)->efficiencyBackground_->GetBinLowEdge(1),(*itInfoList)->efficiencyBackground_->GetBinCenter(bin)-(*itInfoList)->efficiencyBackground_->GetBinWidth(bin));
       (*itInfoList)->significance_->GetXaxis()->SetRangeUser((*itInfoList)->efficiencyBackground_->GetBinLowEdge(1),(*itInfoList)->efficiencyBackground_->GetBinCenter(bin)-(*itInfoList)->efficiencyBackground_->GetBinWidth(bin));
   }

   // some plots style
   (*itInfoList)->efficiencyBackground_->GetYaxis()->SetTitle("Efficiency");
   (*itInfoList)->efficiencyBackground_->GetXaxis()->SetTitleSize(0.035);
   (*itInfoList)->efficiencyBackground_->GetXaxis()->SetLabelSize(0.035);
   (*itInfoList)->efficiencyBackground_->GetXaxis()->SetTitleOffset(1.05);
   (*itInfoList)->efficiencyBackground_->GetYaxis()->SetTitleSize(0.035);
   (*itInfoList)->efficiencyBackground_->GetYaxis()->SetLabelSize(0.035);
   (*itInfoList)->efficiencyBackground_->GetYaxis()->SetTitleOffset(1.20);
   (*itInfoList)->efficiencyBackground_->GetYaxis()->SetRangeUser(0.,(*itInfoList)->efficiencyBackground_->GetMaximum()*1.3);
   (*itInfoList)->efficiencyBackground_->SetLineColor(kRed);
   (*itInfoList)->efficiencyBackground_->SetLineWidth(2);

   (*itInfoList)->efficiencySignal_->SetLineColor(kBlue);
   (*itInfoList)->efficiencySignal_->SetLineWidth(2);


   (*itInfoList)->significance_->SetLineColor(kBlack);
   (*itInfoList)->significance_->SetLineWidth(2);

   (*itInfoList)->efficiencyBackground_->Draw("histl");
   (*itInfoList)->efficiencySignal_->Draw("samehistl");   

    // fill the info for the maximum significance and error and location
    Int_t maxbin = (*itInfoList)->significance_->GetMaximumBin();    
    (*itInfoList)->maxSig_    = (*itInfoList)->significance_->GetBinContent(maxbin) ;
    (*itInfoList)->maxSigErr_ = (*itInfoList)->significance_->GetBinError(maxbin) ;      
    (*itInfoList)->maxbin_    = maxbin ;
    
    TString opt = Form( "%%%is:  (%%8.4g,%%8.4g)    %%9.4g   %%10.6g  %%8.7g  %%8.7g %%8.4g %%8.4g",maxLenTitle);     
    std::cout << "--- "<< Form( opt.Data(), (*itInfoList)->methodTitle_.Data(), numberSignalEvents, numberBackgroundEvents, (*itInfoList)->significance_->GetXaxis()->GetBinCenter(maxbin),(*itInfoList)->maxSig_,(*itInfoList)->efficiencySignal_->GetBinContent(maxbin)*numberSignalEvents, (*itInfoList)->efficiencyBackground_->GetBinContent( maxbin )*numberBackgroundEvents,(*itInfoList)->efficiencySignal_->GetBinContent(maxbin), (*itInfoList)->efficiencyBackground_->GetBinContent(maxbin) ) <<std::endl;


    // scale the significance plot to one in order to be put in the same canvas with efficiencies
    (*itInfoList)->significance_->Scale(1/(*itInfoList)->significance_->GetMaximum());
    (*itInfoList)->significance_->Draw("samehistl");
   
    // Draw legend                                                                                                                                                                       
    TLegend *legend1= new TLegend( cSignificance_->GetLeftMargin()+0.05,1-cSignificance_->GetTopMargin()-0.17,cSignificance_->GetLeftMargin()+0.25,1-cSignificance_->GetTopMargin()-0.02);
    legend1->SetFillStyle(0);
    legend1->SetFillColor(0);
    legend1->SetTextFont(42);
    legend1->SetTextSize(0.032);
    legend1->SetBorderSize(0.);
    legend1->SetMargin(0.3);

    legend1->AddEntry((*itInfoList)->efficiencySignal_,"Signal Efficiency","L");
    legend1->AddEntry((*itInfoList)->efficiencyBackground_,"Background Efficiency","L");
    legend1->Draw("same");
   
    TLegend *legend2= new TLegend(cSignificance_->GetLeftMargin()+0.4,1-cSignificance_->GetTopMargin()-0.09,1-cSignificance_->GetRightMargin()-0.2,1-cSignificance_->GetTopMargin()-0.02);
    legend2->SetFillStyle(0);
    legend2->SetFillColor(0);
    legend2->SetBorderSize(0.);
    legend2->SetTextFont(42);
    legend2->SetTextSize(0.032);
    legend2->SetMargin(0.3);

    legend2->AddEntry((*itInfoList)->significance_,Form("Significance %s",(*this).GetLatexFormula().Data()),"L");
    legend2->Draw("same");

    // line to indicate maximum efficiency                                                                                                                                               
    TLine* effline = 0 ;
    if ((*itInfoList)->methodTitle_.Contains("Cuts")) effline = new TLine(0.1,1,1,1);
    else effline = new TLine((*itInfoList)->efficiencyBackground_->GetBinLowEdge(1),1,(*itInfoList)->efficiencyBackground_->GetBinCenter(bin)-(*itInfoList)->efficiencyBackground_->GetBinWidth(bin),1);

    effline->SetLineWidth(3);
    effline->SetLineColor(6);
    effline->Draw("same");

    cSignificance_->Update();


    // print comments                                                                                                                                                                    
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
    tex = new TLatex(0.2324,0.92,"Simulation Preliminary");
    tex->SetNDC();
    tex->SetTextFont(52);
    tex->SetTextSize(0.0304);
    tex->SetLineWidth(2);
    tex->Draw();

    TString name = Form("For %.4g signal and %.4g background",numberSignalEvents,numberBackgroundEvents);
    TPaveText* line1 = new TPaveText(0.22,0.21,0.5,0.26,"NDC");
    line1->SetBorderSize(0);
    line1->SetFillColor(0);
    line1->SetFillStyle(0);
    line1->AddText(name.Data());
    line1->SetTextSize(0.027);
    line1->SetTextFont(42);
   //   line1->Draw("same");
  
    TPaveText* line2 = new TPaveText(0.22,0.15,0.5,0.2,"NDC");
    name = Form("max significance at %0.4g, cut at %0.4g",(*itInfoList)->maxSig_,(*itInfoList)->significance_->GetXaxis()->GetBinCenter((*itInfoList)->maxbin_));
    line2->AddText(name);
    line2->SetBorderSize(0);
    line2->SetFillColor(0);
    line2->SetFillStyle(0);
    line2->SetTextSize(0.027);
    line2->SetTextFont(42);
    //line2->Draw("same");

    TGaxis* rightAxis = new TGaxis(cSignificance_->GetUxmax(),cSignificance_->GetUymin(),cSignificance_->GetUxmax(),cSignificance_->GetUymax(),0,(*itInfoList)->maxSig_*1.3,510,"+L");
    rightAxis->SetLineColor (kBlack);
    rightAxis->SetLabelColor(kBlack);
    rightAxis->SetTitleColor(kBlack);

    rightAxis->SetLabelSize(0.035);
    rightAxis->SetTitleSize(0.035);
    rightAxis->SetTitleOffset(1.25);
    rightAxis->SetTextFont(42);
    rightAxis->SetLabelFont(42);
    rightAxis->SetTitle("Significance");
    rightAxis->Draw();

    cSignificance_->Update();

    TString baseName ;
    if(UseSignalEfficiency && UseBackgroundEfficiency){
     if(thisMethod_ >= method_index) continue ;
       baseName  = Form("mva_significance_eff_%s_%s",(*itInfoList)->methodTitle_.Data(),inputName.c_str());

     if(stype == 0)      (*this).PrintImage(cSignificance_, std::string(Form("%s/%s_S_over_B_file",outputPlotDirectory.c_str(),baseName.Data())));
     else if(stype == 1) (*this).PrintImage(cSignificance_, std::string(Form("%s/%s_S_over_sqrtB",outputPlotDirectory.c_str(),baseName.Data())));
     else if(stype == 2) (*this).PrintImage(cSignificance_, std::string(Form("%s/%s_S_over_sqrtSB",outputPlotDirectory.c_str(),baseName.Data())));
     else if(stype == 3) (*this).PrintImage(cSignificance_, std::string(Form("%s/%s_pval",outputPlotDirectory.c_str(),baseName.Data())));

     if ((*itInfoList)->methodTitle_.Contains("Cuts")) thisMethod_++;

    }
    else if( !UseSignalEfficiency && !UseBackgroundEfficiency){
     if(thisMethod_ > method_index) continue ;
       baseName  = Form("mva_significance_%s_%s",(*itInfoList)->methodTitle_.Data(),inputName.c_str());

     if(stype == 0) (*this).PrintImage(cSignificance_, std::string(Form("%s/%s_S_over_B",outputPlotDirectory.c_str(),baseName.Data())));
     else if(stype == 1) (*this).PrintImage(cSignificance_, std::string(Form("%s/%s_S_over_sqrtB",outputPlotDirectory.c_str(),baseName.Data())));
     else if(stype == 2) (*this).PrintImage(cSignificance_, std::string(Form("%s/%s_S_over_sqrtSB",outputPlotDirectory.c_str(),baseName.Data())));
     else if(stype == 3) (*this).PrintImage(cSignificance_, std::string(Form("%s/%s_pval",outputPlotDirectory.c_str(),baseName.Data())));

     if ((*itInfoList)->methodTitle_.Contains("Cuts")) thisMethod_++;

    }

    if(thisMethod_ == method_index) thisMethod_ = 0;
   
    //   delete rightAxis ;
    delete cSignificance_;
    legend1->Delete() ;   
    legend2->Delete() ;
    line1->Delete() ;
    line2->Delete() ;
    effline->Delete() ;

  }

  return ;

}


// Set and store the method name in a vector in order to be used in the legend 
void TMVAGlob::SetMethodName(const std::vector<std::string> & SetMethodName){

  inputMethodName_    = SetMethodName;
  originalMethodName_ = SetMethodName;

  for(size_t iMethod = 0 ; iMethod < inputMethodName_.size() ; iMethod++){
    inputMethodName_.at(iMethod) = std::string(TString(inputMethodName_.at(iMethod)).ReplaceAll(":_:","#"));
    inputMethodName_.at(iMethod) = std::string(TString(inputMethodName_.at(iMethod)).ReplaceAll(":__:","_{"));
    inputMethodName_.at(iMethod) = std::string(TString(inputMethodName_.at(iMethod)).ReplaceAll(":___:","}"));
    inputMethodName_.at(iMethod) = std::string(TString(inputMethodName_.at(iMethod)).ReplaceAll("//"," "));
    inputMethodName_.at(iMethod) = std::string(TString(inputMethodName_.at(iMethod)).ReplaceAll("/../","("));
    inputMethodName_.at(iMethod) = std::string(TString(inputMethodName_.at(iMethod)).ReplaceAll("/::/",")"));
    inputMethodName_.at(iMethod) = std::string(TString(inputMethodName_.at(iMethod)).ReplaceAll("/:::/","="));
  }       

  for(size_t iMethod = 0 ; iMethod < originalMethodName_.size() ; iMethod++){
    originalMethodName_.at(iMethod) = std::string(TString(originalMethodName_.at(iMethod)).ReplaceAll(":_:","_"));
    originalMethodName_.at(iMethod) = std::string(TString(originalMethodName_.at(iMethod)).ReplaceAll(":__:","_"));
    originalMethodName_.at(iMethod) = std::string(TString(originalMethodName_.at(iMethod)).ReplaceAll(":___:","_"));
    originalMethodName_.at(iMethod) = std::string(TString(originalMethodName_.at(iMethod)).ReplaceAll("//","_"));
    originalMethodName_.at(iMethod) = std::string(TString(originalMethodName_.at(iMethod)).ReplaceAll("/","_"));
    originalMethodName_.at(iMethod) = std::string(TString(originalMethodName_.at(iMethod)).ReplaceAll("/../","_"));
    originalMethodName_.at(iMethod) = std::string(TString(originalMethodName_.at(iMethod)).ReplaceAll("/::/","_"));
    originalMethodName_.at(iMethod) = std::string(TString(originalMethodName_.at(iMethod)).ReplaceAll("/:::/","_"));
  }       

  return ; 

}
