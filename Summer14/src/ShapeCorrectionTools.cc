#include "../include/ShapeCorrectionTools.h"

//Jet Cleansing
JetCleanser makeJVFCleanser(fastjet::JetDefinition subjet_def, std::string projectmode, double fcut, int nsj ){

   JetCleanser::input_mode tmp_projectmode = JetCleanser::input_nc_separate;
   if(projectmode=="ATLAS") tmp_projectmode = JetCleanser::input_nc_together;
     JetCleanser tmpCleanser(subjet_def, JetCleanser::jvf_cleansing, tmp_projectmode);
     if( fcut>0 && nsj<0){ tmpCleanser.SetTrimming(fcut);}
     else if( fcut<0 && nsj>0 ){ tmpCleanser.SetFiltering(nsj);}
     else if( fcut>0 && nsj>0 ){ tmpCleanser.SetGroomingParameters(fcut,nsj);}

     return tmpCleanser;
}

JetCleanser makeLinearCleanser(fastjet::JetDefinition subjet_def, double linear_para0,std::string projectmode, double fcut, int nsj){

   JetCleanser::input_mode tmp_projectmode=JetCleanser::input_nc_separate;
   if(projectmode=="ATLAS") tmp_projectmode=JetCleanser::input_nc_together;
   JetCleanser tmpCleanser(subjet_def, JetCleanser::linear_cleansing, tmp_projectmode);
   tmpCleanser.SetLinearParameters(linear_para0);

   if( fcut>0 && nsj<0){ tmpCleanser.SetTrimming(fcut);}
   else if( fcut<0 && nsj>0 ){ tmpCleanser.SetFiltering(nsj);}
   else if( fcut>0 && nsj>0 ){ tmpCleanser.SetGroomingParameters(fcut,nsj);}

   return tmpCleanser;
}

JetCleanser makeGausCleanser(fastjet::JetDefinition subjet_def, double gaus_para0, double gaus_para1, double gaus_para2, double gaus_para3, std::string projectmode, double fcut, int nsj ){

  JetCleanser::input_mode tmp_projectmode=JetCleanser::input_nc_separate;
  if(projectmode=="ATLAS") tmp_projectmode=JetCleanser::input_nc_together;
  JetCleanser tmpCleanser(subjet_def, JetCleanser::gaussian_cleansing, tmp_projectmode);
  tmpCleanser.SetGaussianParameters(gaus_para0, gaus_para1, gaus_para2, gaus_para3);

  if( fcut>0 && nsj<0){ tmpCleanser.SetTrimming(fcut);}
  else if( fcut<0 && nsj>0 ){ tmpCleanser.SetFiltering(nsj);}
  else if( fcut>0 && nsj>0 ){ tmpCleanser.SetGroomingParameters(fcut,nsj);}

  return tmpCleanser;
}

