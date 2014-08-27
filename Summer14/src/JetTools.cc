#include "../include/JetTools.h"

//Phill JEC functions to correct on the fly
void loadPhil(std::string iName,std::vector<TGraph*> &iCorr) {
  TFile *lFile = TFile::Open(iName.c_str());
  for(int i0 = 0; i0 < 4; i0++) {
    std::stringstream pSS0,pSS1,pSS2;
    pSS0 << "PFEta" << i0;
    pSS1 << "CHSEta" << i0;
    pSS2 << "PuppiEta" << i0;
    TGraph* lF0 = (TGraph*) lFile->FindObjectAny(pSS0.str().c_str());
    TGraph* lF1 = (TGraph*) lFile->FindObjectAny(pSS1.str().c_str());
    TGraph* lF2 = (TGraph*) lFile->FindObjectAny(pSS2.str().c_str());
    iCorr.push_back(lF0);
    iCorr.push_back(lF1);
    iCorr.push_back(lF2);
  }
  return;
}

double correctPhil(double iPt,double iEta,int iAlgo,std::vector<TGraph*> &iCorr) {
  double lPt = iPt;
  if(lPt > 3000) return 1.;
  int iId = iAlgo;
  if(fabs(iEta) > 1.5) iId += 3;
  if(fabs(iEta) > 2.5) iId += 3;
  if(fabs(iEta) > 3.0) iId += 3;
  if(iPt < 20) lPt = 20;
  double pCorr = iCorr[iId]->Eval(lPt);
  pCorr/=lPt;
  return pCorr;
}

// Get standard JEC information
double correction( PseudoJet &iJet,FactorizedJetCorrector *iJetCorr,double iRho){
  iJetCorr->setJetPt (iJet.pt());
  iJetCorr->setJetEta(iJet.eta());
  iJetCorr->setJetPhi(iJet.phi());
  iJetCorr->setJetE  (iJet.e());
  iJetCorr->setJetA  (iJet.area());
  iJetCorr->setRho(iRho);
  iJetCorr->setJetEMF(-99.0);
  double jetcorr= iJetCorr->getCorrection();
  return jetcorr;
}

double unc( PseudoJet &iJet,JetCorrectionUncertainty *iJetUnc){
  if(fabs(iJet.eta()) > 5. || fabs(iJet.pt()) < 10.) return 1.;
  iJetUnc->setJetPt ( iJet.pt()  );
  iJetUnc->setJetEta( iJet.eta() );
  double jetunc = iJetUnc->getUncertainty(true);
  return jetunc;
}

// Fastjet algorithm definition
fastjet::JetAlgorithm get_algo(string algo){
  fastjet::JetAlgorithm jetalgo;
  if (algo=="kt") jetalgo = fastjet::kt_algorithm ;
  else if (algo=="ca") jetalgo = fastjet::cambridge_algorithm ;
  else if (algo=="ak") jetalgo = fastjet::antikt_algorithm ;
  else if (algo=="KT") jetalgo = fastjet::kt_algorithm ;
  else if (algo=="CA") jetalgo = fastjet::cambridge_algorithm ;
  else if (algo=="AK") jetalgo = fastjet::antikt_algorithm ;
  else if (algo=="kt_algorithm") jetalgo = fastjet::kt_algorithm ;
  else if (algo=="cambridge_algorithm") jetalgo = fastjet::cambridge_algorithm ;
  else if (algo=="antikt_algorithm") jetalgo = fastjet::antikt_algorithm ;
  else if (algo=="0") jetalgo = fastjet::kt_algorithm ;
  else if (algo=="1") jetalgo = fastjet::cambridge_algorithm ;
  else if (algo=="2") jetalgo = fastjet::antikt_algorithm ;
  else jetalgo = fastjet::antikt_algorithm ;
  return jetalgo;
}

///////// divide jet particles after ClusterSequence in neutrals, charged from PV and PU charged                                                                                       
void getConstitsForCleansing(vector<PseudoJet> inputs, vector<PseudoJet> &oNeutrals, vector<PseudoJet> &oChargedLV, vector<PseudoJet> &oChargedPU){
  for (unsigned int i = 0; i < inputs.size(); i++){
    if (inputs[i].user_index() == 0) oNeutrals.push_back(inputs[i]);
    else if (fabs(inputs[i].user_index()) <= 2) oChargedLV.push_back(inputs[i]);
    else if (fabs(inputs[i].user_index()) >= 3) oChargedPU.push_back(inputs[i]);
  }
}

//// Matching functions

Bool_t isMatching( fastjet::PseudoJet j1, fastjet::PseudoJet j2, Double_t deltaR){
	Double_t eta1=j1.eta();
	Double_t eta2=j2.eta();
	Double_t phi1=j1.phi();
	Double_t phi2=j2.phi();
	Double_t tmpR=TMath::Sqrt( (eta1-eta2)*(eta1-eta2) + (phi1-phi2)*(phi1-phi2) );
	if (tmpR<deltaR)return 1;
	else return 0;
}


/// Compute Jet flavour
int computeGenJetFlavour(const PseudoJet & iJet, TClonesArray* fGenParticles, const double & jetR){
  // calculate the flavour of the genJet in order to match it with reco one                                                                                                              
  int tempParticle = -1;
  int tempPartonHighestPt = -1;
  float maxPt = 0.;
  // Loop on the gen particle in order to take the GenPartons                                                                                                                             
  baconhep::TGenParticle *pPartTmp = NULL ;
  baconhep::TGenParticle *pPartTmpD = NULL ;

  for( int iGenParticle = 0; iGenParticle < fGenParticles->GetEntriesFast(); iGenParticle++){
    pPartTmp = (baconhep::TGenParticle*)((*fGenParticles)[iGenParticle]);
    if(!(abs(pPartTmp->pdgId) == 1 || abs(pPartTmp->pdgId) == 2 || abs(pPartTmp->pdgId) == 3 || abs(pPartTmp->pdgId) == 4 || abs(pPartTmp->pdgId) == 5 || abs(pPartTmp->pdgId) == 21)) continue ;
    if(pPartTmp->pt < 0.001) continue; // if gen pt is less than 1MeV, then don't bother matching, the p4 is probably buggy                                                               
    int nDaughters = 0;
    int nPartonDaughters= 0;
    if(pPartTmp->status!=3){
      for( int iGenParticleD = 0; iGenParticleD < fGenParticles->GetEntriesFast(); iGenParticleD++){//9,entries loop,fill the vector particles with PF particles                          
	pPartTmpD = (baconhep::TGenParticle*)((*fGenParticles)[iGenParticleD]);
	if(iGenParticleD!=iGenParticle and pPartTmpD->parent == iGenParticle ){
	  nDaughters++ ;
	  if(abs(pPartTmpD->pdgId) == 1 || abs(pPartTmpD->pdgId) == 2 || abs(pPartTmpD->pdgId) == 3 || abs(pPartTmpD->pdgId) == 4 || abs(pPartTmpD->pdgId) == 5 || abs(pPartTmpD->pdgId)\
	     == 6 || abs(pPartTmpD->pdgId) == 21) nPartonDaughters++;
	}
      }
      if(nDaughters <= 0) continue;
      if(nPartonDaughters > 0) continue ;
    }

    double dPhi = fabs(pPartTmp->phi-iJet.phi());
    if(dPhi > 2.*TMath::Pi()-dPhi) dPhi =  2.*TMath::Pi()-dPhi;
    double deltaR = sqrt(fabs(pPartTmp->eta-iJet.eta())*fabs(pPartTmp->eta-iJet.eta())+dPhi*dPhi);
    if(deltaR > jetR) continue;
    if(tempParticle == -1 && ( abs(pPartTmp->pdgId) == 4 ) ) tempParticle = iGenParticle;
    if(abs(pPartTmp->pdgId) == 5 ) tempParticle = iGenParticle;
    if(pPartTmp->pt > maxPt){
      maxPt = pPartTmp->pt;
      tempPartonHighestPt = iGenParticle;
    }
  }
  if (tempParticle == -1) tempParticle = tempPartonHighestPt;
  if (tempParticle == -1) return 0;
  else { pPartTmp = (baconhep::TGenParticle*)((*fGenParticles)[tempParticle]);
    return int(pPartTmp->pdgId);
  }
}

