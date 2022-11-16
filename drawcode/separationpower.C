#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH1F.h"
#include "TMath.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "math.h"
#include <TCanvas.h>
#include <TApplication.h>
#include <TH1F.h>
#include <TH2F.h>
//#include <time.h>
using namespace std;
bool eff(char *a);
bool noeff(char *a);
TMultiGraph *mg = new TMultiGraph();
TLegend* leg = new TLegend(0.65,0.7,0.8,0.8);
bool option = 1;//true is the situation that separation power change with momentum while false is the situation of cos(theta);
void separationpower(){
 //------------seita>45    
    char efile[] = "../rootfile/e2e2h_X_1_9_1.root";//consider cluster counting efficiency
    char nfile[] = "../Noeff_1_100_1.root";//ideal
    char* efffile = efile;
    char* noefffile = nfile;

    eff(efffile);
    // noeff(noefffile);
    //TCanvas *c1 = new TCanvas("c1","kaon PID efficiency with error bars",200,10,700,500);
    mg->SetTitle("k/pi  separation powser");
    if(option){
        mg->GetXaxis()->SetTitle("P(Gev/c)");
    }
    else{
        mg->GetXaxis()->SetTitle("cos#theta");
    }
    mg->GetYaxis()->SetTitle("separation power");
    TCanvas* c = new TCanvas("c","kaon/pi separation power",200,10,700,500);
    mg->Draw("APL");
    // leg->Draw();
    c->Print("fig/separationpower.pdf");
    c->Print("fig/separationpower.png");



}
bool eff(char* a){
    if(option){
        Float_t halfrange = 1;
        Float_t distance = 0.1574/2;//GeV
        int m = 20;
        Float_t sp_a[20];
        Float_t p_a[20];
        for (int n = 0;n<m; n++){
            TFile* f = TFile::Open(a);
            TTree* Delphes = (TTree*)f->Get("Delphes");
            TString cut_common = "Track.L_DC>0 && Track.L>0 && Track.Nclusters>0";
            TString cut_pi_PID = "(Track.Truth_PID==211||Track.Truth_PID==-211)";
            TString cut_k_PID = "(Track.Truth_PID==321||Track.Truth_PID==-321)";
            TString cut_mom_max = "Track.P<="+std::to_string(pow(10,-/*0.3*/0.2+(n+1)*distance));
            TString cut_mom_min = "Track.P>="+std::to_string(pow(10,-/*0.3*/0.2+n*distance));
            TString cut_sita = "abs(Track.CosTheta)>=0&&abs(Track.CosTheta)<=0.854";
            TString cut_pion = cut_common + "&&" + cut_mom_min + "&&" + cut_mom_max + "&&" + cut_sita + "&&" + cut_pi_PID;
            TString cut_kaon = cut_common + "&&" + cut_mom_min + "&&" + cut_mom_max + "&&" + cut_sita + "&&" + cut_k_PID;

            Float_t mean_ncl_pi = 0.;
            Float_t mean_sigma_pi = 0.;

            Float_t mean_ncl_k = 0.;
            Float_t mean_sigma_k = 0.;

            Delphes->Draw("Nclusters>>ncl_pi",cut_pion);
            TH1 *ncl_pi = (TH1*)gDirectory->Get("ncl_pi");
            mean_ncl_pi = ncl_pi->GetMean();
            cout<<"ncl_pi  =  "<<mean_ncl_pi<<endl;
             Delphes->Draw("Nclusters_err>>sigma_pi",cut_pion);
            TH1 *sigma_pi = (TH1*)gDirectory->Get("sigma_pi");
            mean_sigma_pi = sigma_pi->GetMean();
            cout<<"sigma_pi = "<<mean_sigma_pi<<endl;

            Delphes->Draw("Nclusters>>ncl_k",cut_kaon);
            TH1 *ncl_k = (TH1*)gDirectory->Get("ncl_k");
            mean_ncl_k = ncl_k->GetMean();
            cout<<"ncl_k  =  "<<mean_ncl_k<<endl;
            Delphes->Draw("Nclusters_err>>sigma_k",cut_kaon);
            TH1 *sigma_k = (TH1*)gDirectory->Get("sigma_k");
            mean_sigma_k = sigma_k->GetMean();
            cout<<"sigma_k = "<<mean_sigma_k<<endl;

            Float_t sp = 0.;
            Float_t p =n*distance;
            sp = TMath::Abs((mean_ncl_pi-mean_ncl_k)/(mean_sigma_pi+mean_sigma_k)*2);
            cout<< "P: "<<n*distance <<endl;
            cout << "separation power: " << sp << endl;
            sp_a[n]=sp;
            p_a[n]=(pow(10,-/*0.3*/0.0969+n*distance)+pow(10,-/*0.3*/0.2+(n+1)*distance))/2;
            f->Close();
        }
        auto g = new TGraph(20,p_a,sp_a);
        g->SetMarkerStyle(21);
        g->SetMarkerColor(4);
        g->SetLineColor(4);
        leg->AddEntry(g,"eff");
        mg->Add(g);
    }
    else{
        Float_t distance = 0.1;
        int m = 10;
        Float_t sp_a[10];
        Float_t cos_a[10];
    

        for (int n = 0;n<10; n++){
            TFile* f = TFile::Open(a);
            TTree* Delphes = (TTree*)f->Get("Delphes");
            TString cut_common = "Track.L_DC>0 && Track.L>0 && Track.Nclusters>0";
            TString cut_pi_PID = "(Track.Truth_PID==211||Track.Truth_PID==-211)";
            TString cut_k_PID = "(Track.Truth_PID==321||Track.Truth_PID==-321)";
            TString cut_cos_min = "abs(Track.CosTheta)<"+std::to_string(n*distance+0.05);
            TString cut_cos_max = "abs(Track.CosTheta)>"+std::to_string(n*distance);
            TString cut_mom_max = "Track.P<22";
            TString cut_mom_min = "Track.P>18";
            // TString cut_sita = "abs(Track.CosTheta)>=0&&abs(Track.CosTheta)<=0.854";
            // TString cut_pt = "Track.PT > 1";
            TString cut_pion = cut_common + "&&" + cut_mom_min + "&&" + cut_mom_max +  "&&" + cut_cos_max + "&&" + cut_cos_min + "&&" + cut_pi_PID /* + "&&" + cut_pt */;
            TString cut_kaon = /**/cut_common + "&&" + cut_mom_min + "&&" + cut_mom_max +  "&&" + cut_cos_max + "&&" + cut_cos_min + "&&" + cut_k_PID /* + "&&" + cut_pt */;

            Float_t mean_ncl_pi = 0.;
            Float_t mean_sigma_pi = 0.;

            Float_t mean_ncl_k = 0.;
            Float_t mean_sigma_k = 0.;

            Delphes->Draw("Nclusters>>ncl_pi",cut_pion);
            TH1 *ncl_pi = (TH1*)gDirectory->Get("ncl_pi");
            mean_ncl_pi = ncl_pi->GetMean();
            cout<<"ncl_pi  =  "<<mean_ncl_pi<<endl;
            cout<<"n_pi = "<<ncl_pi->GetEntries()<<endl;
            Delphes->Draw("Nclusters_err>>sigma_pi",cut_pion);
            TH1 *sigma_pi = (TH1*)gDirectory->Get("sigma_pi");
            mean_sigma_pi = sigma_pi->GetMean();
            cout<<"sigma_pi = "<<mean_sigma_pi<<endl;

            Delphes->Draw("Nclusters>>ncl_k",cut_kaon);
            TH1 *ncl_k = (TH1*)gDirectory->Get("ncl_k");
            mean_ncl_k = ncl_k->GetMean();
            cout<<"ncl_k  =  "<<mean_ncl_k<<endl;
            cout<<"n_k = "<<ncl_k->GetEntries()<<endl;
            Delphes->Draw("Nclusters_err>>sigma_k",cut_kaon);
            TH1 *sigma_k = (TH1*)gDirectory->Get("sigma_k");
            mean_sigma_k = sigma_k->GetMean();
            cout<<"sigma_k = "<<mean_sigma_k<<endl;

            Float_t sp = 0.;
            Float_t p =((1-n*distance)+(1-(n+1)*distance))/2;
            sp = TMath::Abs((mean_ncl_pi-mean_ncl_k)/(mean_sigma_pi+mean_sigma_k)*2);
            cout<< "cos: "<<n*distance <<endl;
            cout << "separation power: " << sp << endl;
            sp_a[n]=sp;
            cos_a[n]=n*distance;
            f->Close();
        }
        auto g = new TGraph(10,cos_a,sp_a);
        g->SetMarkerStyle(21);
        g->SetMarkerColor(4);
        g->SetLineColor(4);
        leg->AddEntry(g,"eff ");
        mg->Add(g);


    }
    return 0;
}
bool noeff(char* a){
    if(option){
        Float_t halfrange = 1;
        Float_t distance = 0.1574/2;//GeV
        int m = 20;
        Float_t sp_a[20];
        Float_t p_a[20];
    

        for (int n = 0;n<m; n++){
            TFile* f = TFile::Open(a);
            TTree* Delphes = (TTree*)f->Get("Delphes");
            TString cut_common = "Track.L_DC>0 && Track.L>0 && Track.Nclusters>0";
            TString cut_pi_PID = "(Track.Truth_PID==211||Track.Truth_PID==-211)";
            TString cut_k_PID = "(Track.Truth_PID==321||Track.Truth_PID==-321)";
            TString cut_mom_max = "Track.P<="+std::to_string(pow(10,-/*0.3*/0.2+(n+1)*distance));
            TString cut_mom_min = "Track.P>="+std::to_string(pow(10,-/*0.3*/0.2+n*distance));
            TString cut_sita = "abs(Track.CosTheta)>=0&&abs(Track.CosTheta)<=0.854";
            TString cut_pion = cut_common + "&&" + cut_mom_min + "&&" + cut_mom_max + "&&" + cut_sita + "&&" + cut_pi_PID;
            TString cut_kaon = cut_common + "&&" + cut_mom_min + "&&" + cut_mom_max + "&&" + cut_sita + "&&" + cut_k_PID;

            Float_t mean_ncl_pi = 0.;
            Float_t mean_sigma_pi = 0.;

            Float_t mean_ncl_k = 0.;
            Float_t mean_sigma_k = 0.;

            Delphes->Draw("Nclusters>>ncl_pi",cut_pion);
            TH1 *ncl_pi = (TH1*)gDirectory->Get("ncl_pi");
            mean_ncl_pi = ncl_pi->GetMean();
            cout<<"ncl_pi  =  "<<mean_ncl_pi<<endl;
             Delphes->Draw("Nclusters_err>>sigma_pi",cut_pion);
            TH1 *sigma_pi = (TH1*)gDirectory->Get("sigma_pi");
            mean_sigma_pi = sigma_pi->GetMean();
            cout<<"sigma_pi = "<<mean_sigma_pi<<endl;

            Delphes->Draw("Nclusters>>ncl_k",cut_kaon);
            TH1 *ncl_k = (TH1*)gDirectory->Get("ncl_k");
            mean_ncl_k = ncl_k->GetMean();
            cout<<"ncl_k  =  "<<mean_ncl_k<<endl;
            Delphes->Draw("Nclusters_err>>sigma_k",cut_kaon);
            TH1 *sigma_k = (TH1*)gDirectory->Get("sigma_k");
            mean_sigma_k = sigma_k->GetMean();
            cout<<"sigma_k = "<<mean_sigma_k<<endl;

            Float_t sp = 0.;
            Float_t p =n*distance;
            sp = TMath::Abs((mean_ncl_pi-mean_ncl_k)/(mean_sigma_pi+mean_sigma_k)*2);
            cout<< "P: "<<n*distance <<endl;
            cout << "separation power: " << sp << endl;
            sp_a[n]=sp;
            p_a[n]=(pow(10,-/*0.3*/0.0969+n*distance)+pow(10,-/*0.3*/0.2+(n+1)*distance))/2;
            f->Close();

        }
        auto g = new TGraph(20,p_a,sp_a);
        g->SetMarkerStyle(21);
        g->SetMarkerColor(2);
        g->SetLineColor(2);
        leg->AddEntry(g,"ideal");
        mg->Add(g);


    }
    else{
        Float_t distance = 0.1;
        int m = 10;
        Float_t sp_a[10];
        Float_t cos_a[10];
    
  
        for (int n = 0;n<m; n++){
            TFile* f = TFile::Open(a);
            TTree* Delphes = (TTree*)f->Get("Delphes");
            TString cut_common = "Track.L_DC>0 && Track.L>0 && Track.Nclusters>0";
            TString cut_pi_PID = "(Track.Truth_PID==211||Track.Truth_PID==-211)";
            TString cut_k_PID = "(Track.Truth_PID==321||Track.Truth_PID==-321)";
            TString cut_cos_min = "abs(Track.CosTheta)<"+std::to_string(n*distance+0.1);
            TString cut_cos_max = "abs(Track.CosTheta)>"+std::to_string(n*distance);
            TString cut_mom_max = "Track.P<16";
            TString cut_mom_min = "Track.P>12";  
            TString cut_sita = "abs(Track.CosTheta)>=0&&abs(Track.CosTheta)<=0.854";
            // TString cut_pt = "Track.PT > 1";
            TString cut_pion = cut_common + "&&" + cut_mom_min + "&&" + cut_mom_max +  "&&" + cut_cos_max + "&&" + cut_cos_min + "&&" + cut_pi_PID /* + "&&" + cut_pt */;
            TString cut_kaon = cut_common + "&&" + cut_mom_min + "&&" + cut_mom_max +  "&&" + cut_cos_max + "&&" + cut_cos_min + "&&" + cut_k_PID /* + "&&" + cut_pt */;


            Float_t mean_ncl_pi = 0.;
            Float_t mean_sigma_pi = 0.;

            Float_t mean_ncl_k = 0.;
            Float_t mean_sigma_k = 0.;

            Delphes->Draw("Nclusters>>ncl_pi",cut_pion);
            TH1 *ncl_pi = (TH1*)gDirectory->Get("ncl_pi");
            mean_ncl_pi = ncl_pi->GetMean();
            cout<<"ncl_pi  =  "<<mean_ncl_pi<<endl;
             Delphes->Draw("Nclusters_err>>sigma_pi",cut_pion);
            TH1 *sigma_pi = (TH1*)gDirectory->Get("sigma_pi");
            mean_sigma_pi = sigma_pi->GetMean();
            cout<<"sigma_pi = "<<mean_sigma_pi<<endl;

            Delphes->Draw("Nclusters>>ncl_k",cut_kaon);
            TH1 *ncl_k = (TH1*)gDirectory->Get("ncl_k");
            mean_ncl_k = ncl_k->GetMean();
            cout<<"ncl_k  =  "<<mean_ncl_k<<endl;
            Delphes->Draw("Nclusters_err>>sigma_k",cut_kaon);
            TH1 *sigma_k = (TH1*)gDirectory->Get("sigma_k");
            mean_sigma_k = sigma_k->GetMean();
            cout<<"sigma_k = "<<mean_sigma_k<<endl;

            Float_t sp = 0.;
            Float_t p =((1-n*distance)+(1-(n+1)*distance))/2;
            sp = TMath::Abs((mean_ncl_pi-mean_ncl_k)/(mean_sigma_pi+mean_sigma_k)*2);
            cout<< "cos: "<<n*distance <<endl;
            cout << "separation power: " << sp << endl;
            sp_a[n]=sp;
            cos_a[n]=n*distance;
            f->Close();
        }
        auto g = new TGraph(10,cos_a,sp_a);
        g->SetMarkerStyle(21);
        g->SetMarkerColor(2);
        g->SetLineColor(2);
        leg->AddEntry(g,"ideal");
        mg->Add(g);


    }
    return 0;
}
