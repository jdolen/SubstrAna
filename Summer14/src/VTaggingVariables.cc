#include "../include/VTaggingVariables.h"

///////////////////////                                                                                                                                                                  
VTaggingVariables::VTaggingVariables(const PseudoJet & inputJet){

  inputJet_ = inputJet ;
  SelectorIsPureGhost().sift(inputJet_.constituents(), ghosts_, particles_);

}

///////////////////////                                                                                                                                                                  

double VTaggingVariables::computeJetCharge(const double & jetChargeKappa){

  jetChargeKappa_ = jetChargeKappa ;

  double  chargeVal = 0.;
  std::vector<fastjet::PseudoJet>::const_iterator itConstituents =  particles_.begin();
  for ( ; itConstituents != particles_.end(); ++itConstituents){
    if(fabs((*itConstituents).user_index()) > 2) continue ;
    chargeVal += (*itConstituents).user_index()*pow((*itConstituents).pt(),jetChargeKappa_);
  }

  chargeVal /= pow(inputJet_.pt(),jetChargeKappa_ );
  return chargeVal;

}

///////////////////////                                                                                                                                                                  
double VTaggingVariables::computeJetChargeReco(const double & jetChargeKappa){

  jetChargeKappa_ = jetChargeKappa ;

  double  chargeVal = 0.;
  std::vector<fastjet::PseudoJet>::const_iterator itConstituents =  particles_.begin();
  for ( ; itConstituents != particles_.end(); ++itConstituents){
    double charge = 0.;        
    
    if(fabs((*itConstituents).user_index()) > 7) continue ;
    if((*itConstituents).user_index() == 0) charge = 0;
    else if(fabs((*itConstituents).user_index() )<=2) charge = (*itConstituents).user_index() ;
    else if((fabs((*itConstituents).user_index() ) > 2 )) charge = (*itConstituents).user_index() -5 ;
    else charge = 0 ;
    
    chargeVal += charge*pow((*itConstituents).pt(),jetChargeKappa_);
  }

  chargeVal /= pow(inputJet_.pt(),jetChargeKappa_ );
  return chargeVal;

}


///////////////////////                                                                                                                                                                  
double VTaggingVariables::computeNSubJettines(const int & nJettines, const double & beta, const double & R0, const double & Rcut){

 if(inputJet_.structure_ptr() == 0 or inputJet_.structure_ptr()== NULL or inputJet_ == 0) return -99;
 contrib::Nsubjettiness nSub(nJettines,contrib::Njettiness::onepass_kt_axes,beta,R0,Rcut);
 double nSubjettinessValues = nSub(inputJet_);

 return nSubjettinessValues ;
}
  

///////////////////////                                                                                                                                                                  
double VTaggingVariables::computeECF(JetAlgorithm jetAlgoforECF, const double & Rparameter, const int & nPoint, const double & beta, const int & type){

   double ECFValues = 0 ;

   if(!jet_def_forECF_) jet_def_forECF_   = new JetDefinition(jetAlgoforECF,Rparameter); // definition for ECF                         
   if(!clust_seq_forECF_) clust_seq_forECF_ = new ClusterSequence(particles_,*(jet_def_forECF_));         
   if(incluisve_jets_forECF_.empty()) incluisve_jets_forECF_ = clust_seq_forECF_->inclusive_jets(0);//take all the set of pseudojets from clustering                                
   
   if(type == 0){ // compute simple energy correlation function for nPoint
     contrib::EnergyCorrelator C2beta(nPoint,beta,contrib::EnergyCorrelator::pt_R);
     ECFValues = C2beta(incluisve_jets_forECF_[0]);  
     return ECFValues ;
   }
   else if(type == 1){
     contrib::EnergyCorrelatorRatio C2beta(nPoint,beta,contrib::EnergyCorrelator::pt_R); // ECF(N+1,beta)/ECF(N,beta)
     ECFValues = C2beta(incluisve_jets_forECF_[0]);  
     return ECFValues ;
   }
   else if(type == 2){
    contrib::EnergyCorrelatorDoubleRatio C2beta(nPoint,beta,contrib::EnergyCorrelator::pt_R); // ECF(N+1,beta)*ECF(N-1,beta)/ECF(N,beta)^2
    ECFValues = C2beta(incluisve_jets_forECF_[0]);  
    return ECFValues ;
   }

   return ECFValues ;
}


///////////////////////                                                                                                                                                                  
double VTaggingVariables::computeQjets(const int & QJetsPreclustering, const int & QJetsN, const int & seed, const double & jetR){

  JetDefinition jet_def_CA (fastjet::cambridge_algorithm, jetR*10); //large R to cluster all constituents of original jet                                                                
  fastjet::ClusterSequence cs_Recluster (particles_, jet_def_CA);
  vector<fastjet::PseudoJet> jets_Recluster = sorted_by_pt(cs_Recluster.inclusive_jets());
 
  std::vector<fastjet::PseudoJet> constits; 
  unsigned int nqjetconstits = jets_Recluster[0].constituents().size(); // take jet constituent size
  if (nqjetconstits < (unsigned int) QJetsPreclustering) constits = jets_Recluster[0].constituents(); // use of subset of particles
  else constits = jets_Recluster[0].associated_cluster_sequence()->exclusive_subjets_up_to(jets_Recluster[0],QJetsPreclustering); // take the exclusive subjets
  
  double qjet_vol = getQjetVolatility(constits, QJetsN, seed*QJetsN) ;
  constits.clear();
  return qjet_vol;
  
}


double VTaggingVariables::getQjetVolatility(std::vector < fastjet::PseudoJet > constits, const int & QJetsN, const int & seed){

    std::vector<double> qjetmasses;
    double zcut(0.1), dcut_fctr(0.5), exp_min(0.), exp_max(0.), rigidity(0.1), truncationFactor(0.01);

    QjetsPlugin qjet_plugin(zcut, dcut_fctr, exp_min, exp_max, rigidity, truncationFactor);
    JetDefinition qjet_def(&qjet_plugin);

    for(unsigned int ii = 0 ; ii < (unsigned int) QJetsN ; ii++){
      qjet_plugin.SetRandSeed(seed+ii); // new feature in Qjets to set the random seed                                                                                                   
      ClusterSequence qjet_seq(constits, qjet_def);
      
      std::vector<PseudoJet> inclusive_jets2 = sorted_by_pt(qjet_seq.inclusive_jets(5.0));
      if (inclusive_jets2.size()>0) { qjetmasses.push_back( inclusive_jets2[0].m() ); }     
    }
    
    // find RMS of a vector                                                                                                                                                               
    double qjetsRMS = FindRMS( qjetmasses );
    // find mean of a vector                                                                                                                                                              
    double qjetsMean = FindMean( qjetmasses );
    double qjetsVolatility = qjetsRMS/qjetsMean;
    return qjetsVolatility;
    
}
  
//Q jets stuff                                                                                                                                                                            
double VTaggingVariables::FindRMS( const std::vector< double > & qjetmasses ){

  double total = 0.;
  double ctr = 0.;

  for (unsigned int i = 0; i < qjetmasses.size(); i++){
      total = total + qjetmasses[i];
      ctr++;
  }

  double mean = total/ctr;
  double totalsquared = 0.;
  for (unsigned int i = 0; i < qjetmasses.size(); i++){
      totalsquared += (qjetmasses[i] - mean)*(qjetmasses[i] - mean) ;
  }
  double RMS = sqrt( totalsquared/ctr );
  return RMS;
}

double VTaggingVariables::FindMean(const std::vector< double > & qjetmasses){

 double total = 0.;
 double ctr = 0.;
 for (unsigned int i = 0; i < qjetmasses.size(); i++){
  total = total + qjetmasses[i];
  ctr++;
 }
  
 return total/ctr;  
}

///////////////////////                                                                                                                                                                  
void VTaggingVariables::setInputJet(const PseudoJet & inputJet){

  inputJet_ = inputJet ;
  SelectorIsPureGhost().sift(inputJet_.constituents(), ghosts_, particles_);
  jet_def_forECF_   = 0 ;
  clust_seq_forECF_ = 0 ;
  incluisve_jets_forECF_.clear() ;


}

//////////////////////
double VTaggingVariables::computeQGLikelihood(QGLikelihoodCalculator* qgLikelihood, const double & jetCorrection, const double & rho){

  QCLikelihoodVariables_.clear();

  float sum_weight = 0., sum_deta = 0., sum_dphi = 0., sum_deta2 = 0., sum_dphi2 = 0., sum_detadphi = 0., sum_pt = 0.;
  QCLikelihoodVariables_["pt"]  = inputJet_.pt()*jetCorrection;
  QCLikelihoodVariables_["eta"] = inputJet_.eta();
  QCLikelihoodVariables_["rhoIso"] = rho;

  //  std::cout<<" jet pt "<<inputJet_.pt()<<" jet eta "<<inputJet_.eta()<<std::endl;

  vector<PseudoJet>::const_iterator itConstituents =  particles_.begin();
  int nChg_QC = 0, nChg_ptCut = 0, nNeutral_ptCut = 0;
  for( ; itConstituents!= particles_.end(); itConstituents++){
    if ((*itConstituents).user_index() == 0 and (*itConstituents).pt() > 1.0 ) nNeutral_ptCut++ ;
    else if (fabs((*itConstituents).user_index()) <= 2 and (*itConstituents).pt() > 1.0){ //leading vertex                                                                        
      nChg_ptCut++;
      nChg_QC++;
    }
    else if (fabs((*itConstituents).user_index()) >= 3 and (*itConstituents).pt() > 1.0) nChg_ptCut++;

    //    std::cout<<" constituent "<<(*itConstituents).user_index()<<" pt "<<(*itConstituents).pt()<<" eta "<<(*itConstituents).eta()<<std::endl;
    sum_pt     += (*itConstituents).pt();
    sum_weight += (*itConstituents).pt()*(*itConstituents).pt();
    sum_deta   += ((*itConstituents).eta()-inputJet_.eta())*(*itConstituents).pt()*(*itConstituents).pt();
    sum_dphi   +=  2*atan(tan((((*itConstituents).phi()-inputJet_.phi()))/2))*(*itConstituents).pt()*(*itConstituents).pt();
    sum_deta2  += pow(((*itConstituents).eta()-inputJet_.eta()),2)*(*itConstituents).pt()*(*itConstituents).pt();
    sum_dphi2  += pow(2*atan(tan((((*itConstituents).phi()-inputJet_.phi()))/2)),2)*(*itConstituents).pt()*(*itConstituents).pt();
    sum_detadphi += ((*itConstituents).eta()-inputJet_.eta())*2*atan(tan((((*itConstituents).phi()-inputJet_.phi()))/2))*(*itConstituents).pt()*(*itConstituents).pt();
  }


  //Calculate axis and ptD                                                                                                                                                                
  float a = 0., b = 0., c = 0.;
  float ave_deta = 0., ave_dphi = 0., ave_deta2 = 0., ave_dphi2 = 0.;
  if(sum_weight > 0){
    QCLikelihoodVariables_["ptD"] = sqrt(sum_weight)/sum_pt;
    //    std::cout<<" ptD : "<<sum_weight<<"  "<<sum_pt<<"  "<<QCLikelihoodVariables_["ptD"]<<std::endl;
    ave_deta = sum_deta/sum_weight;
    ave_dphi = sum_dphi/sum_weight;
    ave_deta2 = sum_deta2/sum_weight;
    ave_dphi2 = sum_dphi2/sum_weight;
    a = ave_deta2 - ave_deta*ave_deta;
    b = ave_dphi2 - ave_dphi*ave_dphi;
    c = -(sum_detadphi/sum_weight - ave_deta*ave_dphi);
  }
  else{ QCLikelihoodVariables_["ptD"] = 0.;
    //    std::cout<<" ptD : "<<sum_weight<<"  "<<sum_pt<<"  "<<QCLikelihoodVariables_["ptD"]<<std::endl;
  }

  float delta = sqrt(fabs((a-b)*(a-b)+4*c*c));
  if(a+b+delta > 0) QCLikelihoodVariables_["axis1"] = sqrt(0.5*(a+b+delta));
  else QCLikelihoodVariables_["axis1"] = 0.;
  if(a+b-delta > 0) QCLikelihoodVariables_["axis2"] = sqrt(0.5*(a+b-delta));
  else QCLikelihoodVariables_["axis2"] = 0.;

  QCLikelihoodVariables_["mult"] = (nChg_QC + nNeutral_ptCut);

  //std::cout<<" axis1 : "<<QCLikelihoodVariables_["axis1"]<<"  "<<QCLikelihoodVariables_["axis2"]<<"  "<<nChg_QC + nNeutral_ptCut<<std::endl;

  return qgLikelihood->QGvalue(QCLikelihoodVariables_);
 
} 

double VTaggingVariables::computePullAngle(const std::vector<PseudoJet> & subjets, const double & jetR){

  PseudoJet subjet0(0,0,0,0) ;
  PseudoJet subjet1(0,0,0,0) ;
  for(unsigned int i0 = 0; i0 < subjets.size(); i0++) {
   double pDR = reco::deltaR(subjets.at(i0).eta(),subjets.at(i0).phi(),inputJet_.eta(),inputJet_.phi());   
   if(pDR > jetR) continue;    
   if(subjet0.pt() < subjets.at(i0).pt()) {
	subjet1 = subjet0;
	subjet0 = subjets.at(i0);            
   } 
   else if(subjet1.pt() < subjets.at(i0).pt()) subjet1 = subjets.at(i0);   
  }  
  if(subjet0.pt() == 0 and subjet1.pt() == 0) return -20;

  // work in dy-dphi space of subjet0
  TVector2 lPull = jetPull(subjet0);
  TVector2 lJet(subjet1.rapidity()-subjet0.rapidity(),reco::deltaPhi(subjet1.phi(),subjet0.phi()));
  double lThetaP = lPull.DeltaPhi(lJet);
  return lThetaP;
}
  
TVector2 VTaggingVariables::jetPull(const PseudoJet & jet, const int & type){

 double dYSum=0, dPhiSum=0;
 const unsigned int nconst = jet.constituents().size();
 
 for(unsigned int iconst = 0; iconst < nconst; iconst++) {

   double pt_i = 0, y_i = 0, phi_i = 0;
   PseudoJet particle = jet.constituents().at(iconst);

   if(type==0) { // PF jet pull
     pt_i  = particle.pt();
     y_i   = particle.rapidity();
     phi_i = particle.phi();    
   } 
   else if(type==1) { // charged jet pull component
     if(particle.user_index()!=0) {
	  pt_i  = particle.pt();
	  y_i   = particle.rapidity();
	  phi_i = particle.phi();
     }
          
   } else if(type==2) { // neutral jet pull component
     if(particle.user_index()==0) {
	  pt_i  = particle.pt();
	  y_i   = particle.rapidity();
	  phi_i = particle.phi();
     }      
   } else assert(0);
      
    
   double dY = y_i - jet.rapidity();
   double dPhi = reco::deltaPhi(phi_i,jet.phi());
   double weight = pt_i*sqrt(dY*dY + dPhi*dPhi);
   dYSum += weight*dY;
   dPhiSum += weight*dPhi;
 }
  
 return TVector2(dYSum/jet.pt(), dPhiSum/jet.pt());

}

