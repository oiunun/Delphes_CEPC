#ifdef __CLING__
R__LOAD_LIBRARY(../../libDelphes.so)
#endif
#include "../../classes/DelphesClasses.h"
#include "../../external/ExRootAnalysis/ExRootTreeReader.h"
#include <iostream>
#include "TMath.h"
#include <string>
#include "TFile.h"
#include "TF2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TTree.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TClonesArray.h"

double Inv_Mass(Jet *jet,TLorentzVector &P_pi,TLorentzVector &P_em,TLorentzVector &P_ep);
void InfoInformation(TClonesArray *branchParticle);
void initLatexCode();
string getLatexCode(int pdgid);
map< int, string> LatexCode;
void Bgana()
{
  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add("../../../rootfile/mal/output_pythia6_PID/bb_*.root");
  TFile *f = new TFile("bgana1.root");
  TTree *t1 = (TTree*)f->Get("tree2");
  Int_t id_entry;
  t1->SetBranchAddress("id_entry",&id_entry);

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);

  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchPFC = treeReader->UseBranch("ParticleFlowCandidate");
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TH1F *hist_bg = new TH1F("bg", "bg", 100,5,5.6);
  TH1F *hist_q2 = new TH1F("q2", "q2", 100,0,4);

  double mass_b=5.27934;
  for(Long64_t entry = 0; entry <t1->GetEntries(); ++entry)
  {
    // Load selected branches with data from specified event
    t1->GetEntry(entry);
    treeReader->ReadEntry(id_entry);
    if(branchJet->GetEntriesFast()!= 2) continue;
    Double_t Mass[2]={0.0,0.0};
    Double_t mass=0.0;
    TLorentzVector P_pi,P_em,P_ep;
    TLorentzVector P_pi_1,P_em_1,P_ep_1;
    TLorentzVector P_pi_2,P_em_2,P_ep_2;
    for(Int_t i=0;i<2;i++){
      Jet* jet=(Jet*) branchJet->At(i);
      if(i==0)  Mass[i]=Inv_Mass(jet,P_pi_1,P_em_1,P_ep_1);
      else  Mass[i]=Inv_Mass(jet,P_pi_2,P_em_2,P_ep_2);
    }
    if(abs(Mass[0]-mass_b)<abs(Mass[1]-mass_b)){
      mass=Mass[0];
      P_pi=P_pi_1;
      P_em=P_em_1;
      P_ep=P_ep_1;
    }
    else{
      mass=Mass[1];
      P_pi=P_pi_2;
      P_em=P_em_2;
      P_ep=P_ep_2;
    }
   if(P_pi.Pt()<4 || P_em.Pt()<4 || P_ep.Pt()<4) continue;
   if(mass<5 || mass>5.6) continue;
   double q2=(P_em+P_ep).M();
   if(q2<3.3 && q2>3) continue;
   hist_q2->Fill(q2);
   hist_bg->Fill(mass);
   InfoInformation(branchParticle);
  } 
  TCanvas *c1=new TCanvas("c1","c1",700,700);
  gStyle->SetOptStat("im");
  hist_bg->Draw();
  c1->Print("fig/B_e1e1pi/bgana.png"); 
  c1->Print("fig/B_e1e1pi/bgana.pdf"); 
  hist_q2->Draw();
  c1->Print("fig/B_e1e1pi/bgq2.png"); 
  c1->Clear();
  f->Close();

  delete c1;
  delete treeReader;
  
}

double Inv_Mass(Jet *jet,TLorentzVector &P_pi,TLorentzVector &P_em,TLorentzVector &P_ep){
    int n_Jetp = jet->Constituents.GetEntriesFast();
    if(n_Jetp==0) return 0;
    std::vector < Int_t > id_pi;
    std::vector < Int_t > id_ep;
    std::vector < Int_t > id_em;
    ParticleFlowCandidate *particleflowcandidate=NULL;
    for(int i=0;i<n_Jetp;i++){
       particleflowcandidate = (ParticleFlowCandidate*) jet->Constituents.At(i);
      if((particleflowcandidate->Eem+particleflowcandidate->Ehad)/particleflowcandidate->P>0.4 && (particleflowcandidate->Eem/particleflowcandidate->P)<0.8   && abs(particleflowcandidate->PID)==211 /* &&  particleflowcandidate->P>3 */){
        id_pi.push_back(i);
      }
      if(particleflowcandidate->Eem/particleflowcandidate->P>0.8  /* &&  particleflowcandidate->P>3 */  && particleflowcandidate->Charge == 1 ){
        id_ep.push_back(i);
      }
      if(particleflowcandidate->Eem/particleflowcandidate->P>0.8 /* && particleflowcandidate->P>3 */ && particleflowcandidate->Charge == -1){
        id_em.push_back(i);
      }
      
    }

    float inv_mass = 0;
    if(!id_em.empty() && !id_ep.empty() && !id_pi.empty()) {       
      vector<int>::iterator it_pi;
      vector<int>::iterator it_ep;
      vector<int>::iterator it_em;
      TLorentzVector P_pi_i; 
      TLorentzVector P_ep_i;
      TLorentzVector P_em_i;
      for(it_pi=id_pi.begin();it_pi!=id_pi.end();it_pi++){
        ParticleFlowCandidate* canpi = (ParticleFlowCandidate*) jet->Constituents.At(*it_pi);
	GenParticle *parpi=(GenParticle*)canpi->Particles.At(0);
        P_pi_i=canpi->P4();
        for(it_ep=id_ep.begin();it_ep!=id_ep.end();it_ep++){
         ParticleFlowCandidate* canep = (ParticleFlowCandidate*) jet->Constituents.At(*it_ep);
	 GenParticle *parep=(GenParticle*)canep->Particles.At(0);
          P_ep_i = canep->P4();
          for(it_em=id_em.begin();it_em!=id_em.end();it_em++){
           ParticleFlowCandidate* canem = (ParticleFlowCandidate*) jet->Constituents.At(*it_em);
	   GenParticle *parem=(GenParticle*)canem->Particles.At(0);
            P_em_i = canem->P4();
            Double_t mass = (P_pi_i+P_ep_i+P_em_i).M();
            if(abs(mass-5.27934)<abs(inv_mass-5.27934)) {
              inv_mass = mass;
	      P_pi=P_pi_i;
	      P_em=P_em_i;
	      P_ep=P_ep_i;
   	     if(P_pi.Pt()<4 || P_em.Pt()<4 || P_ep.Pt()<4 || inv_mass>5.6 || inv_mass<5) continue;
             cout<<"emM1 = "<<parem->M1<<"  epM1 = "<<parep->M1<<"  piM1 = "<<parpi->M1<<endl;
            }
          }
        }
      }
    }
    id_em.clear();
    id_ep.clear();
    id_pi.clear();
    return inv_mass;
}
void InfoInformation(TClonesArray *branchParticle){

  initLatexCode();
  GenParticle *particle=NULL;
  for(int i=0;i<branchParticle->GetEntries();i++){
    particle=(GenParticle*)branchParticle->At(i);
    printf("id=%3d, PID =%6d, pdgid=%s, M1=%3d, M2=%3d, D1=%3d, D2=%3d, status=%2d,costheta=%6.3f, phi=%6.3f \n", i,particle->PID,getLatexCode(particle->PID).c_str(),particle->M1,particle->M2,particle->D1,particle->D2,particle->Status,particle->P4().CosTheta(),particle->P4().Phi());

  }
 cout<<endl;
 cout<<"----------------Next Event-------------"<<endl;
 cout<<endl;



}

void initLatexCode(){
	LatexCode.clear();
	LatexCode.insert(map<int,string,less<int> >::value_type(              1, "                 d"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             -1, "            anti-d"));
	LatexCode.insert(map<int,string,less<int> >::value_type(              2, "                 u"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             -2, "            anti-u"));
	LatexCode.insert(map<int,string,less<int> >::value_type(              3, "                 s"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             -3, "            anti-s"));
	LatexCode.insert(map<int,string,less<int> >::value_type(              4, "                 c"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             -4, "            anti-c"));
	LatexCode.insert(map<int,string,less<int> >::value_type(              5, "                 b"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             -5, "            anti-b"));
	LatexCode.insert(map<int,string,less<int> >::value_type(              6, "                 t"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             -6, "            anti-t"));
	LatexCode.insert(map<int,string,less<int> >::value_type(              7, "                b'"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             -7, "           anti-b'"));
	LatexCode.insert(map<int,string,less<int> >::value_type(              8, "                t'"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             -8, "           anti-t'"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             21, "                 g"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             11, "                e-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -11, "                e+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             12, "              nu_e"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -12, "         anti-nu_e"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             13, "               mu-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -13, "               mu+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             14, "             nu_mu"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -14, "        anti-nu_mu"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             15, "              tau-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -15, "              tau+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             16, "            nu_tau"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -16, "       anti-nu_tau"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             17, "                L-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -17, "                L+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             18, "              nu_L"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -18, "         anti-nu_L"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             22, "             gamma"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -22, "          gammaFSR"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10022, "              vpho"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20022, "          Cerenkov"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             23, "                Z0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             24, "                W+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -24, "                W-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             25, "            Higgs0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             28, "           reggeon"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             29, "           pomeron"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             32, "               Z'0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             33, "              Z''0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             34, "               W'+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -34, "               W'-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             35, "           Higgs'0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             36, "                A0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             37, "            Higgs+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -37, "            Higgs-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             40, "                R0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -40, "           anti-R0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             41, "               Xu0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             42, "               Xu+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -42, "               Xu-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             81, "          specflav"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             82, "          rndmflav"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -82, "     anti-rndmflav"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             83, "          phasespa"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             84, "          c-hadron"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -84, "     anti-c-hadron"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             85, "          b-hadron"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -85, "     anti-b-hadron"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             86, "          t-hadron"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -86, "     anti-t-hadron"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             87, "         b'-hadron"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -87, "    anti-b'-hadron"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             88, "         t'-hadron"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -88, "    anti-t'-hadron"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             89, "            Wvirt+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -89, "            Wvirt-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             90, "           diquark"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -90, "      anti-diquark"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             91, "           cluster"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             92, "            string"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             93, "             indep"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             94, "          CMshower"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             95, "          SPHEaxis"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             96, "          THRUaxis"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             97, "           CLUSjet"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             98, "           CELLjet"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             99, "             table"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            111, "               pi0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            211, "               pi+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -211, "               pi-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            210, "          pi_diff+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -210, "          pi_diff-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20111, "           pi(2S)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20211, "           pi(2S)+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -20211, "           pi(2S)-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            221, "               eta"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20221, "           eta(2S)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            331, "              eta'"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            113, "              rho0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            110, "         rho_diff0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            213, "              rho+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -213, "              rho-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          30113, "          rho(2S)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          30213, "          rho(2S)+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -30213, "          rho(2S)-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          40113, "          rho(3S)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          40213, "          rho(3S)+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -40213, "          rho(3S)-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(       -9040213, "        rho(2150)-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9040213, "        rho(2150)+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(       -9040113, "        rho(2150)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            223, "             omega"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            220, "        omega_diff"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          30223, "         omega(2S)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            333, "               phi"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10111, "              a_00"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10211, "              a_0+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10211, "              a_0-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9000221, "          f_0(600)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9010221, "               f_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10221, "              f'_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10113, "              b_10"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10213, "              b_1+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10213, "              b_1-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10223, "               h_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10333, "              h'_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20113, "              a_10"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20213, "              a_1+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -20213, "              a_1-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20223, "               f_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            115, "              a_20"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            215, "              a_2+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -215, "              a_2-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            225, "               f_2"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10221, "         f_0(1370)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          50221, "         f_0(1500)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            335, "              f'_2"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9020221, "         eta(1405)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         100331, "         eta(1475)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10335, "        eta2(1870)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10331, "         f_0(1710)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            229, "         f_4(2050)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20333, "              f'_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        8888888, "         f_0(1790)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9000223, "         f_1(1510)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9050225, "         f_2(1950)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9070225, "         f_2(2150)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9070221, "         f_0(2200)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9020225, "         f_2(1640)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9040225, "         f_2(1910)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9050225, "         f_2(1950)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9080221, "         eta(2225)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9040221, "         eta(1760)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9999999, "           x(1835)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            311, "                K0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -311, "           anti-K0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            310, "              K_S0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            130, "              K_L0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            321, "                K+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -321, "                K-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            313, "               K*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -313, "          anti-K*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            323, "               K*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -323, "               K*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9000311, "             kapa0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(       -9000311, "        anti-kapa0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9000321, "             kapa+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(       -9000321, "             kapa-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10311, "             K_0*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10311, "        anti-K_0*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10321, "             K_0*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10321, "             K_0*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10313, "              K_10"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10313, "         anti-K_10"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10323, "              K_1+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10323, "              K_1-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            315, "             K_2*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -315, "        anti-K_2*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            325, "             K_2*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -325, "             K_2*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20313, "             K'_10"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -20313, "        anti-K'_10"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20323, "             K'_1+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -20323, "             K'_1-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         100313, "              K'*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        -100313, "         anti-K'*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         100323, "              K'*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        -100323, "              K'*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          30313, "             K''*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -30313, "        anti-K''*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          30323, "             K''*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -30323, "             K''*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10315, "              K_20"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10315, "         anti-K_20"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10325, "              K_2+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10325, "              K_2-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            317, "             K_3*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -317, "        anti-K_3*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            327, "             K_3*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -327, "             K_3*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            319, "             K_4*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -319, "        anti-K_4*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            329, "             K_4*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -329, "             K_4*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9060221, "         f_0(2100)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            411, "                D+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -411, "                D-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            421, "                D0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -421, "           anti-D0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            422, "               D0H"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -422, "               D0L"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            413, "               D*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -413, "               D*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            423, "               D*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -423, "          anti-D*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10411, "             D_0*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10411, "             D_0*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10421, "             D_0*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10421, "        anti-D_0*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10413, "              D_1+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10413, "              D_1-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10423, "              D_10"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10423, "         anti-D_10"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            415, "             D_2*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -415, "             D_2*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            425, "             D_2*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -425, "        anti-D_2*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20413, "             D'_1+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -20413, "             D'_1-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20423, "             D'_10"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -20423, "        anti-D'_10"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            431, "              D_s+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -431, "              D_s-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            433, "             D_s*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -433, "             D_s*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10431, "            D_s0*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10431, "            D_s0*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10433, "             D_s1+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10433, "             D_s1-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            435, "            D_s2*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -435, "            D_s2*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20433, "            D'_s1+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -20433, "            D'_s1-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          30411, "            D(2S)+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -30411, "            D(2S)-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          30421, "            D(2S)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -30421, "       anti-D(2S)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          30413, "           D*(2S)+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -30413, "           D*(2S)-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          30423, "           D*(2S)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -30423, "      anti-D*(2S)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            511, "                B0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -511, "           anti-B0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            150, "               B0L"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            510, "               B0H"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            521, "                B+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -521, "                B-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            513, "               B*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -513, "          anti-B*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            523, "               B*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -523, "               B*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10511, "             B_0*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10511, "        anti-B_0*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10521, "             B_0*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10521, "             B_0*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10513, "              B_10"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10513, "         anti-B_10"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10523, "              B_1+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10523, "              B_1-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            515, "             B_2*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -515, "        anti-B_2*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            525, "             B_2*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -525, "             B_2*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20513, "             B'_10"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -20513, "        anti-B'_10"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20523, "             B'_1+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -20523, "             B'_1-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            531, "              B_s0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -531, "         anti-B_s0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            350, "             B_s0L"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            530, "             B_s0H"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            533, "             B_s*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -533, "        anti-B_s*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10531, "            B_s0*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10531, "       anti-B_s0*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10533, "             B_s10"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10533, "        anti-B_s10"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            535, "            B_s2*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -535, "       anti-B_s2*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20533, "            B'_s10"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -20533, "       anti-B'_s10"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            541, "              B_c+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -541, "              B_c-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            543, "             B_c*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -543, "             B_c*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10541, "            B_c0*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10541, "            B_c0*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10543, "             B_c1+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10543, "             B_c1-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            545, "            B_c2*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -545, "            B_c2*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20543, "            B'_c1+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -20543, "            B'_c1-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            441, "             eta_c"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20441, "         eta_c(2S)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            443, "             J/psi"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         100443, "           psi(2S)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          30443, "         psi(3770)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9000443, "         psi(4040)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9010443, "         psi(4160)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9030443, "         psi(4260)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9040443, "         psi(4360)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9020443, "         psi(4415)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10443, "               h_c"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10441, "            chi_c0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20443, "            chi_c1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            445, "            chi_c2"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            551, "             eta_b"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20551, "         eta_b(2S)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          40551, "         eta_b(3S)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            553, "           Upsilon"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          30553, "       Upsilon(2S)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          60553, "       Upsilon(3S)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          70553, "       Upsilon(4S)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          80553, "       Upsilon(5S)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10553, "               h_b"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          40553, "           h_b(2P)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         100553, "           h_b(3P)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10551, "            chi_b0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20553, "            chi_b1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            555, "            chi_b2"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          30551, "        chi_b0(2P)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          50553, "        chi_b1(2P)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10555, "        chi_b2(2P)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          50551, "        chi_b0(3P)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         110553, "        chi_b1(3P)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20555, "        chi_b2(3P)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          40555, "        eta_b2(1D)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          60555, "        eta_b2(2D)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         120553, "     Upsilon_1(1D)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          30555, "     Upsilon_2(1D)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            557, "     Upsilon_3(1D)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         130553, "     Upsilon_1(2D)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          50555, "     Upsilon_2(2D)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10557, "     Upsilon_3(2D)"));
	//
	LatexCode.insert(map<int,string,less<int> >::value_type(          10222, "           sigma_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           1114, "            Delta-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -1114, "       anti-Delta+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           2110, "           n_diffr"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -2110, "      anti-n_diffr"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           2112, "                n0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -2112, "           anti-n0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           2114, "            Delta0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -2114, "       anti-Delta0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           2210, "           p_diff+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -2210, "      anti-p_diff-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           2212, "                p+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -2212, "           anti-p-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           2214, "            Delta+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -2214, "       anti-Delta-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           2224, "           Delta++"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -2224, "      anti-Delta--"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3112, "            Sigma-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3112, "       anti-Sigma+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3114, "           Sigma*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3114, "      anti-Sigma*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3122, "           Lambda0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3122, "      anti-Lambda0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          13122, "     Lambda(1405)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -13122, "anti-Lambda(1405)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3124, "     Lambda(1520)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3124, "anti-Lambda(1520)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          23122, "     Lambda(1600)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -23122, "anti-Lambda(1600)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          33122, "     Lambda(1670)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -33122, "anti-Lambda(1670)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          13124, "     Lambda(1690)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -13124, "anti-Lambda(1690)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          43122, "     Lambda(1800)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -43122, "anti-Lambda(1800)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          53122, "     Lambda(1810)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -53122, "anti-Lambda(1810)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3126, "     Lambda(1820)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3126, "anti-Lambda(1820)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          13126, "     Lambda(1830)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -13126, "anti-Lambda(1830)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          13212, "      Sigma(1660)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -13212, " anti-Sigma(1660)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          13214, "      Sigma(1670)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -13214, " anti-Sigma(1670)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          23212, "      Sigma(1750)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -23212, " anti-Sigma(1750)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3216, "      Sigma(1775)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3216, " anti-Sigma(1775)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3212, "            Sigma0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3212, "       anti-Sigma0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3214, "           Sigma*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3214, "      anti-Sigma*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3222, "            Sigma+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3222, "       anti-Sigma-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3224, "           Sigma*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3224, "      anti-Sigma*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3312, "               Xi-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3312, "          anti-Xi+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3314, "              Xi*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3314, "         anti-Xi*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3322, "               Xi0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3322, "          anti-Xi0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3324, "              Xi*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3324, "         anti-Xi*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3334, "            Omega-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3334, "       anti-Omega+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          14122, "        Lambda_c*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -14122, "   anti-Lambda_c*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          14124, "        Lambda_c*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -14124, "   anti-Lambda_c*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4112, "          Sigma_c0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4112, "     anti-Sigma_c0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4114, "         Sigma_c*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4114, "    anti-Sigma_c*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4212, "          Sigma_c+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4212, "     anti-Sigma_c-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4214, "         Sigma_c*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4214, "    anti-Sigma_c*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4222, "         Sigma_c++"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4222, "    anti-Sigma_c--"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4224, "        Sigma_c*++"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4224, "   anti-Sigma_c*--"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4312, "            Xi'_c0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4312, "       anti-Xi'_c0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4322, "            Xi'_c+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4322, "       anti-Xi'_c-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4324, "            Xi_c*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4324, "       anti-Xi_c*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4122, "         Lambda_c+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4122, "    anti-Lambda_c-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4132, "             Xi_c0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4132, "        anti-Xi_c0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4232, "             Xi_c+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4232, "        anti-Xi_c-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4314, "            Xi_c*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4314, "       anti-Xi_c*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4332, "          Omega_c0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4332, "     anti-Omega_c0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4334, "         Omega_c*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4334, "    anti-Omega_c*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5112, "          Sigma_b-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5112, "     anti-Sigma_b+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5114, "         Sigma_b*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5114, "    anti-Sigma_b*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5122, "         Lambda_b0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5122, "    anti-Lambda_b0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5132, "             Xi_b-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5132, "        anti-Xi_b+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5212, "          Sigma_b0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5212, "     anti-Sigma_b0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5214, "         Sigma_b*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5214, "    anti-Sigma_b*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5222, "          Sigma_b+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5222, "     anti-Sigma_b-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5224, "         Sigma_b*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5224, "    anti-Sigma_b*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5232, "             Xi_b0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5232, "        anti-Xi_b0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5312, "            Xi'_b-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5312, "       anti-Xi'_b+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5314, "            Xi_b*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5314, "       anti-Xi_b*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5322, "            Xi'_b0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5322, "       anti-Xi'_b0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5324, "            Xi_b*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5324, "       anti-Xi_b*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5332, "          Omega_b-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5332, "     anti-Omega_b+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5334, "         Omega_b*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5334, "    anti-Omega_b*+"));
	/*
	LatexCode.insert(map<int,string,less<int> >::value_type(           1101, "              dd_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -1101, "         anti-dd_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           2101, "              ud_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -2101, "         anti-ud_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           2201, "              uu_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -2201, "         anti-uu_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3101, "              sd_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3101, "         anti-sd_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3201, "              su_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3201, "         anti-su_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3301, "              ss_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3301, "         anti-ss_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4101, "              cd_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4101, "         anti-cd_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4201, "              cu_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4201, "         anti-cu_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4301, "              cs_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4301, "         anti-cs_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4401, "              cc_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4401, "         anti-cc_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5101, "              bd_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5101, "         anti-bd_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5201, "              bu_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5201, "         anti-bu_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5301, "              bs_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5301, "         anti-bs_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5401, "              bc_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5401, "         anti-bc_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5501, "              bb_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5501, "         anti-bb_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           1103, "              dd_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -1103, "         anti-dd_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           2103, "              ud_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -2103, "         anti-ud_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           2203, "              uu_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -2203, "         anti-uu_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3103, "              sd_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3103, "         anti-sd_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3203, "              su_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3203, "         anti-su_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3303, "              ss_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3303, "         anti-ss_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4103, "              cd_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4103, "         anti-cd_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4203, "              cu_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4203, "         anti-cu_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4303, "              cs_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4303, "         anti-cs_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4403, "              cc_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4403, "         anti-cc_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5103, "              bd_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5103, "         anti-bd_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5203, "              bu_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5203, "         anti-bu_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5303, "              bs_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5303, "         anti-bs_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5403, "              bc_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5403, "         anti-bc_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5503, "              bb_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5503, "         anti-bb_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           1011, "          deuteron"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -1011, "     anti-deuteron"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           1021, "           tritium"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -1021, "      anti-tritium"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           1012, "               He3"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -1012, "          anti-He3"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           1022, "             alpha"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -1022, "        anti-alpha"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            100, "          geantino"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            101, "   chargedgeantino"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          30343, "               Xsd"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -30343, "          anti-Xsd"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          30353, "               Xsu"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -30353, "          anti-Xsu"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          30373, "               Xdd"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -30373, "          anti-Xdd"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          30383, "               Xdu"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -30383, "          anti-Xdu"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          30363, "               Xss"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -30363, "          anti-Xss"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             51, "         dummy00_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             52, "         dummy10_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             53, "         dummy01_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             54, "         dummy11_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -51, "    anti-dummy00_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -52, "    anti-dummy10_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -53, "    anti-dummy01_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -54, "    anti-dummy11_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             55, "         dummy00_2"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             56, "         dummy10_2"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             57, "         dummy01_2"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             58, "         dummy11_2"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -55, "    anti-dummy00_2"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -56, "    anti-dummy10_2"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -57, "    anti-dummy01_2"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -58, "    anti-dummy11_2"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             61, "           chi_c0p"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             62, "           chi_c1p"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             63, "           chi_c2p"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             64, "             X3940"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             65, "             Y3940"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             66, "             xvpho"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             67, "            eta_c2"));
	*/
}

string getLatexCode(int pdgid){
	string code="           Unknown";
	map< int, string>::iterator pos = LatexCode.find(pdgid); 
	if (pos != LatexCode.end()) 
	{
		code=pos->second;
	}
	return code;
}
