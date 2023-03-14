#include <iostream>
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
#ifdef __CLING__
R__LOAD_LIBRARY(../delphes/libDelphes.so)
#include "../delphes/classes/DelphesClasses.h"
#include "../delphes/external/ExRootAnalysis/ExRootTreeReader.h"
#endif

using namespace std;

void draw(
      const char * inrootfile  = "./root/E240.Pe2e2h_bb.e0.p0.whizard195/e2e2h_bb.e0.p0.00001.root",
      const char * outrootfile = "./root/e2e2h_bb",
      const Int_t lab=1)
{
   TChain chain("Delphes");
   chain.Add(inrootfile);

   ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
   Long64_t numberOfEntries = treeReader->GetEntries();

   TClonesArray *branchPFC  = treeReader->UseBranch("ParticleFlowCandidate");


   printf("%s\n",inrootfile);
   const int np = 128, ne=10000; 
   char outname[100]; 
   sprintf(outname, "%s_%4.4d.root",outrootfile,0);
   printf("%s\n",outname);
   TFile * f1 = TFile::Open (outname,"RECREATE");
   f1->SetCompressionAlgorithm(ROOT::kLZ4);
   f1->SetCompressionLevel(4);
   float En[np], Ch[np], Px[np], Py[np], Pz[np],logEn[np], logPt[np];
   float D0[np], Z0[np], Ma[np], CT[np], TH[np], FI[np], errD0[np], errZ0[np], Ee[np], Eh[np];
   float VP[np];
   int nP = 0, MSK[np]; 
   int lab_bb=0, lab_cc=0, lab_qq=0;


   //cout<<lab<<" "<<LAB<<endl; 
   if( lab>0 && lab<4){
    if(lab == 1) lab_bb=1;
    if(lab == 2) lab_cc=1;
    if(lab == 3) lab_qq=1;
   }
   else{
      printf("wrong label %5d\n", lab);
      exit(1);
   } 

   TTree * t1 = new TTree("tree",outname);
   t1->Branch("nP",  &nP,   "nP/I");
   t1->Branch("Ma",   Ma,   "Ma[128]/F");
   t1->Branch("En",   En,   "En[128]/F");
   t1->Branch("Px",   Px,   "Px[128]/F");
   t1->Branch("Py",   Py,   "Py[128]/F");
   t1->Branch("Pz",   Pz,   "Pz[128]/F");
   t1->Branch("logEn",logEn,"logEn[128]/F");
   t1->Branch("logPt",logPt,"logPt[128]/F");
   t1->Branch("Ch",   Ch,   "Ch[128]/F");
   t1->Branch("D0",   D0,   "D0[128]/F");
   t1->Branch("Z0",   Z0,   "Z0[128]/F");
   t1->Branch("VP",   VP,   "VP[128]/F");
   t1->Branch("Ee",   Ee,   "Ee[128]/F");
   t1->Branch("Eh",   Eh,   "Eh[128]/F");
   t1->Branch("CT",   CT,   "CT[128]/F");
   t1->Branch("TH",   TH,   "TH[128]/F");
   t1->Branch("FI",   FI,   "FI[128]/F");
   t1->Branch("MSK",  MSK,  "MSK[128]/I");

   t1->Branch("lab_bb", &lab_bb,"lab_bb/I");
   t1->Branch("lab_cc", &lab_cc,"lab_cc/I");
   t1->Branch("lab_qq", &lab_qq,"lab_qq/I");

   // Loop over all events
   printf(" Total nubme of events = %8d\n", numberOfEntries);
   for(Long_t entry = 0; entry < numberOfEntries; ++entry)
   {
      for (int i =0; i<np; i++){
         Ma[i] =0.0; 
         En[i] =0.0; 
         Px[i] =0.0; 
         Py[i] =0.0; 
         Pz[i] =0.0; 
         Ch[i] =0.0; 
         D0[i] =0.0; 
         Z0[i] =0.0; 
         VP[i] =0.0; 
         CT[i] =0.0; 
         FI[i] =0.0; 
         TH[i] =0.0; 
         MSK[i]=0  ; 
         logEn[i] = 0.; 
         logPt[i] = 0.; 
      }
      if( entry%10000==0)
      printf(" Event = %8d\n", entry);
      treeReader->ReadEntry(entry);

      Int_t n_PFC = branchPFC->GetEntries();
      ParticleFlowCandidate *PFC = NULL;

      if(n_PFC > 0)
      {
         int NP = 0; 
         for(int k = 0; k < min(n_PFC, np); k++)
         {
            PFC = (ParticleFlowCandidate*) branchPFC->At(k);
            TLorentzVector P4 = PFC->P4(); 
            if ( P4.Rho()>120 || P4.Rho()<0.1 ) continue; 
            Ch   [NP] =( PFC->Charge  ); 
            Ma   [NP] =( P4.M()       ); 
            En   [NP] =( P4.T()/46   ); 
            Px   [NP] =( P4.X()/46   ); 
            Py   [NP] =( P4.Y()/46   ); 
            Pz   [NP] =( P4.Z()/46   ); 
            CT   [NP] =( P4.CosTheta());
            TH   [NP] =( P4.Theta()   );
            FI   [NP] =( P4.Phi()     );
            Ee   [NP] =( PFC->Eem/46    ); 
            Eh   [NP] =( PFC->Ehad/46   ); 
            MSK  [NP] = 1; 
            logEn[NP] = -log(En[NP]); 
            logPt[NP] = -log(Px[NP]*Px[NP] + Py[NP]*Py[NP]); 

            //printf("Mass = %8.3f, Charge = %2.0f, En = %8.3f, D0 = %8.3f,%8.3f\n", Ma[k], Ch[k], En[k], D0[k], PFC->D0 );
            if ( fabs(Ch[NP]) > 0 && fabs(PFC->ErrorDZ)>1e-6 && fabs(PFC->ErrorD0)>1e-6 && fabs(PFC->DZ) < 200 && fabs(PFC->D0) < 200 
                  )
            {
               Z0   [NP] = ( PFC->DZ/10 ); 
               D0   [NP] = ( PFC->D0/10 ); 
               errD0[NP] = ( PFC->ErrorD0/10); 
               errZ0[NP] = ( PFC->ErrorDZ/10); 
               double chid = D0[NP]/errD0[NP];
               double chiz = Z0[NP]/errZ0[NP];
               VP[NP] = (TMath::Prob(chid*chid+chiz*chiz,2));
            }else{
               Z0   [NP] = 0.0; 
               D0   [NP] = 0.0;
               errD0[NP] = 0.0;
               errZ0[NP] = 0.0; 
               VP   [NP] = 0.0;
            }
            NP++; 
         }
         nP = NP; 
         t1->Fill();
      }
      if ( entry % ne==(ne-1) && entry>0 ){
         //f1->Write();
         t1->Write(0, TObject::kWriteDelete, 0);
         f1->Close("R");
         delete f1; 
         //t1->Delete("");

         if ( entry < numberOfEntries ) { 
            sprintf(outname, "%s_%4.4d.root",outrootfile,(entry+1)/ne);
            printf("%s\n",outname);
            f1 = TFile::Open (outname,"RECREATE");
            f1->SetCompressionAlgorithm(ROOT::kLZ4);
            f1->SetCompressionLevel(4);
            t1 = new TTree("tree",outname); 

            t1->Branch("nP",  &nP,   "nP/I");
            t1->Branch("Ma",   Ma,   "Ma[128]/F");
            t1->Branch("En",   En,   "En[128]/F");
            t1->Branch("Px",   Px,   "Px[128]/F");
            t1->Branch("Py",   Py,   "Py[128]/F");
            t1->Branch("Pz",   Pz,   "Pz[128]/F");
            t1->Branch("logEn",logEn,"logEn[128]/F");
            t1->Branch("logPt",logPt,"logPt[128]/F");
            t1->Branch("Ch",   Ch,   "Ch[128]/F");
            t1->Branch("D0",   D0,   "D0[128]/F");
            t1->Branch("Z0",   Z0,   "Z0[128]/F");
            t1->Branch("VP",   VP,   "VP[128]/F");
            t1->Branch("Ee",   Ee,   "Ee[128]/F");
            t1->Branch("Eh",   Eh,   "Eh[128]/F");
            t1->Branch("CT",   CT,   "CT[128]/F");
            t1->Branch("TH",   TH,   "TH[128]/F");
            t1->Branch("FI",   FI,   "FI[128]/F");
            t1->Branch("MSK",  MSK,  "MSK[128]/I");

 
            t1->Branch("lab_bb", &lab_bb,"lab_bb/I");
            t1->Branch("lab_cc", &lab_cc,"lab_cc/I");
            t1->Branch("lab_qq", &lab_qq,"lab_qq/I");

         }
      }
   }
   //f1->Write();
   t1->Write(0, TObject::kWriteDelete, 0);
   f1->Close("R");
   //t1->Delete("");
}

int label(const string file){
   const string prod[ 3]={ "bb","cc", "qq"};

   int lab1=0; 
   for ( int i=0; i<3; i++){
      if ( file.find(prod[i]) != string::npos ) lab1 = i; 

   }
   return lab1; 
}

int main(int argc, char const *argv[])
{
   int is[4] ={0, 0, 4, 10}; 
   for (int a = 0; a < min(4,argc); a++)
   {
      is[a] = atoi(argv[a+1]);
      printf("argv[%d] is %2d\n",a+1,is[a]);
   }

   const string prod[3]={ "bb","cc", "qq"};

   
   char input[100], output[100]; 
   for ( int i=is[0]; i<is[2]; i++){
      for ( int j=is[1]; j<is[3]; j++){
         sprintf( input, "../delphes/rootfile/test_%s.root",prod[i].data()); 
         sprintf(output, "root/%s", prod[i].data()); 
         cout<<output<<endl; 
         int lab = label( input ); 
         draw(input, output, lab); 
      }
   }
   return 0;
}
