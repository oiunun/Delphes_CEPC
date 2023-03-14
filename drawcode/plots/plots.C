// C++ include
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cctype>
#include <string>
#include <sstream>
#include <math.h>
//------ ROOT includes ---------
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TTree.h>
#include <TChain.h>
#include <TNtuple.h>
#include <TBenchmark.h>
#include <TROOT.h>
#include <TLeaf.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TString.h>
#include <cepcplotstyle.h>

using namespace std;
map<string, string>  cut_a_tree;      
void AddHistos( Char_t *bfile, Char_t *cfile, Char_t *lfile)
{
	TKey *key;
	TTree * treeb = NULL, *treec= NULL, *treeq = NULL;
	TFile bbase (bfile);
	TFile cbase (cfile);
	TFile lbase (lfile);
	bool bok = false, cok=false, qok=false;
	if(bbase.IsZombie()) std::cout<< bfile <<" is Zombie"<<std::endl; else bok=true;
	if(cbase.IsZombie()) std::cout<< cfile <<" is Zombie"<<std::endl; else cok=true;
	if(lbase.IsZombie()) std::cout<< lfile <<" is Zombie"<<std::endl; else qok=true;

	TIter next1(cbase.GetListOfKeys());
	//
	while((key=(TKey*)next1())) {
		if(!strcmp(key->GetClassName(),"TTree") ){
			TString tname(key->GetName());
			if(bok)treeb = static_cast<TTree*>(bbase.Get(tname));
			if(cok)treec = static_cast<TTree*>(cbase.Get(tname));
			if(qok)treeq = static_cast<TTree*>(lbase.Get(tname));
			TObjArray* list = treec->GetListOfLeaves();
			TIter next2(list);
			TLeaf* leaf;
			while ( ( leaf = static_cast<TLeaf*>(next2()) ) )
			{
				TString lname(leaf->GetName());
				printf("%s\n", lname.Data());
				TH1D *ha = new TH1D("ha",lname.Data(),1000000,-10000,10000);
				treec->Project("ha",lname);
			   double xmax=-9999,xmin=9999;
			   if(bok)xmax=max(xmax, treeb->GetMaximum(lname));	
			   if(cok)xmax=max(xmax, treec->GetMaximum(lname));	
			   if(qok)xmax=max(xmax, treeq->GetMaximum(lname));	
			   if(bok)xmin=min(xmin, treeb->GetMinimum(lname));	
			   if(cok)xmin=min(xmin, treec->GetMinimum(lname));	
			   if(qok)xmin=min(xmin, treeq->GetMinimum(lname));
			   int nbin=100;
				//
            xmax = int(ha->GetMean()+2.*ha->GetRMS()+1); 
            xmin = int(ha->GetMean()-2.*ha->GetRMS()-1); 
				if(fabs(xmax)<1.1 || fabs(xmin)<1.1) nbin=1000;
				//	
				TH1D *h1 = new TH1D("h1",lname.Data(),nbin,xmin,xmax);
				TH1D *hb = new TH1D("hb",lname.Data(),nbin,xmin,xmax);
				TH1D *hc = new TH1D("hc",lname.Data(),nbin,xmin,xmax);
				TH1D *hq = new TH1D("hq",lname.Data(),nbin,xmin,xmax);
				//	
				if(bok)treeb->Project("hb",lname);
				if(cok)treec->Project("hc",lname);
				if(qok)treeq->Project("hq",lname);
				//
				char filename[256];
				sprintf(filename,"figs/%s", lname.Data());
				//
				PlotDataMC(filename, 
						h1, "",
						hb, "b#bar{b}",
						hc, "c#bar{c}",
						hq, "q#bar{q}",
						 0,         "",
					        false, false, false
						);
				delete h1; 
				delete ha; 
				delete hb; 
				delete hc; 
				delete hq; 
			}
		}
		delete treeb; treeb=NULL;
		delete treec; treec=NULL;
		delete treeq; treeq=NULL;
	}
	bbase.Close();
	cbase.Close();
	lbase.Close();
}

int main(int argc, char *argv[])
{
	// Setbatch
	gROOT->SetBatch();
	AddHistos( "root/h0bq-91.2GeV.root", "root/h0cq-91.2GeV.root", "root/h0lq-91.2GeV.root");   
	return 0;
}
