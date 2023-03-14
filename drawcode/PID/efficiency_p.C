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
// #include <time.h>
using namespace std;
bool eff(char *a);
bool noeff(char *a);
TMultiGraph *mg = new TMultiGraph();
TLegend* leg = new TLegend;
bool particle = 1;//true为kaon；false为pion
void efficiency_p(){
 //------------seita>45   
    char efile[] = "../rootfile/e2e2h_X_1_9_1.root";
    char nfile[] = "../rootfile/Onlynoeff_1_100_1.root"; 
    eff(efile);
    // noeff(nfile);
    if(particle){
        mg->SetTitle("kaon PID efficiency ");
    }
    else{
        mg->SetTitle("pion PID efficiency ");
    }

    mg->GetXaxis()->SetTitle("P(Gev/c)");
    mg->GetYaxis()->SetTitle("efficiency");
    TCanvas* c = new TCanvas("c","kaon PID efficiency",200,10,700,500);
    mg->Draw("APL");
    // leg->Draw();
    c->Print("fig/eff_k_p.pdf");
    c->Print("fig/eff_k_p.png");
}
bool eff(char *a){
    double distance = /*0.178*/0.1574/2;
    int m = 10*2;
    double erry[m];
    double errx[m];
    double accuracy[m];
    double x[m];
    TString cut_truthPID;
    TString cut_totalPID;
    int n_PID=0;
    int Tn_PID=0;
    TFile* f = TFile::Open(a);
    TTree* Delphes = (TTree*)f->Get("Delphes");
    for(int n=0;n<m;n++){
        if(particle){
            cut_truthPID = "((Track.Truth_PID==321||Track.Truth_PID==-321)&&(Track.PID==321||Track.PID==-321))";
            cut_totalPID = "(Track.Truth_PID==321||Track.Truth_PID==-321)";
        }
        else{
            cut_truthPID = "((Track.Truth_PID==211||Track.Truth_PID==-211)&&(Track.PID==211||Track.PID==-211))";
            cut_totalPID = "(Track.Truth_PID==211||Track.Truth_PID==-211)";
        }
        TString cut_mom_min_truth = "Track.Truth_P>="+std::to_string(pow(10,-/*0.3*/0.0969+n*distance))+"&&Track.P>="+std::to_string(pow(10
        ,-/*0.3*/0.0969+n*distance));
        TString cut_mom_max_truth = "Track.Truth_P<="+std::to_string(pow(10,-/*0.3*/0.0969+(n+1)*distance))+"&&Track.P<="+std::to_string(pow(10,-/*0.3*/0.0969+(n+1)*distance));
        TString cut_mom_min_total = "Track.Truth_P>="+std::to_string(pow(10,-/*0.3*/0.0969+n*distance));
        TString cut_mom_max_total = "Track.Truth_P<="+std::to_string(pow(10,-/*0.3*/0.0969+(n+1)*distance));        
        TString cut_sita_truth = "abs(Track.CosTheta)>=0&&abs(Track.CosTheta)<=0.854";
        TString cut_sita_total = "abs(Track.CosTheta)>=0&&abs(Track.CosTheta)<=0.854";
        TString cut_truth = cut_truthPID + "&&" + cut_mom_min_truth + "&&" + cut_mom_max_truth + "&&" + cut_sita_truth;
        TString cut_total = cut_totalPID + "&&" + cut_mom_min_total + "&&" + cut_mom_max_total + "&&" + cut_sita_total;

        Delphes->Draw("Nclusters>>num_truth",cut_truth);
        TH1 *num_truth = (TH1*)gDirectory->Get("num_truth");
        n_PID = num_truth->GetEntries();

        Delphes->Draw("Nclusters>>num_total",cut_total);
        TH1 *num_total = (TH1*)gDirectory->Get("num_total");
        Tn_PID = num_total->GetEntries();

        
        x[n] = (pow(10,-/*0.3*/0.0969+n*distance)+pow(10,-/*0.3*/0.0969+(n+1)*distance))/2;
        accuracy[n] = 1.*n_PID/Tn_PID;
        errx[n] = 0;
        erry[n] = pow(accuracy[n]*(1-accuracy[n])/(1.*Tn_PID),0.5);
        cout << cut_mom_min_total << endl;
        cout << cut_mom_max_total << endl;
        cout<< n_PID << endl;
        cout<< Tn_PID << endl;
        cout<< accuracy[n]  << endl;


    }
    TGraphErrors *gr = new TGraphErrors(m,x,accuracy,errx,erry);
    gr->SetMarkerColor(4);    
    gr->SetLineColor(4);
    gr->SetMarkerStyle(21);
    leg->AddEntry(gr,"eff");
    mg->Add(gr);

    return 0;
}

bool noeff(char *a){
    double distance = /*0.178*/0.1574/2;
    int m = 10*2;
    double erry[m];
    double errx[m];
    double accuracy[m];
    double x[m];
    int n_PID=0;
    int Tn_PID=0;
    TString cut_truthPID;
    TString cut_totalPID;
    TFile* f = TFile::Open(a);
    TTree* Delphes = (TTree*)f->Get("Delphes");
    for(int n=0;n<m;n++){
        if(particle){
            cut_truthPID = "((Track.Truth_PID==321||Track.Truth_PID==-321)&&(Track.PID==321||Track.PID==-321))";
            cut_totalPID = "(Track.Truth_PID==321||Track.Truth_PID==-321)";
        }
        else{
            cut_truthPID = "((Track.Truth_PID==211||Track.Truth_PID==-211)&&(Track.PID==211||Track.PID==-211))";
            cut_totalPID = "(Track.Truth_PID==211||Track.Truth_PID==-211)";
        }
        TString cut_mom_min_truth = "Track.Truth_P>="+std::to_string(pow(10,-/*0.3*/0.0969+n*distance))+"&&Track.P>="+std::to_string(pow(10
        ,-/*0.3*/0.0969+n*distance));
        TString cut_mom_max_truth = "Track.Truth_P<="+std::to_string(pow(10,-/*0.3*/0.0969+(n+1)*distance))+"&&Track.P<="+std::to_string(pow(10,-/*0.3*/0.0969+(n+1)*distance));
        TString cut_mom_min_total = "Track.Truth_P>="+std::to_string(pow(10,-/*0.3*/0.0969+n*distance));
        TString cut_mom_max_total = "Track.Truth_P<="+std::to_string(pow(10,-/*0.3*/0.0969+(n+1)*distance));        
        TString cut_sita_truth = "abs(Track.CosTheta)>=0&&abs(Track.CosTheta)<=0.854";
        TString cut_sita_total = "abs(Track.CosTheta)>=0&&abs(Track.CosTheta)<=0.854";
        TString cut_truth = cut_truthPID + "&&" + cut_mom_min_truth + "&&" + cut_mom_max_truth + "&&" + cut_sita_truth;
        TString cut_total = cut_totalPID + "&&" + cut_mom_min_total + "&&" + cut_mom_max_total + "&&" + cut_sita_total;


        Delphes->Draw("Nclusters>>num_truth",cut_truth);
        TH1 *num_truth = (TH1*)gDirectory->Get("num_truth");
        n_PID = num_truth->GetEntries();

        Delphes->Draw("Nclusters>>num_total",cut_total);
        TH1 *num_total = (TH1*)gDirectory->Get("num_total");
        Tn_PID = num_total->GetEntries();


        x[n] = (pow(10,-/*0.3*/0.0969+n*distance)+pow(10,-/*0.3*/0.0969+(n+1)*distance))/2;
        accuracy[n] = 1.*n_PID/Tn_PID;
        errx[n] = 0;
        erry[n] = pow(accuracy[n]*(1-accuracy[n])/(1.*Tn_PID),0.5);
        cout<< n+1 <<":"<<endl;
        cout << cut_mom_min_total << endl;
        cout << cut_mom_max_total << endl;
        cout<< n_PID << endl;
        cout<< Tn_PID << endl;
        cout<< accuracy[n] << endl;
        cout << erry[n] << endl;


    }
    TGraphErrors *gr = new TGraphErrors(m,x,accuracy,errx,erry);
    gr->SetMarkerColor(2);
    gr->SetLineColor(2);
    gr->SetMarkerStyle(21);
    leg->AddEntry(gr,"ideal");
    mg->Add(gr);
    return 0;
}

