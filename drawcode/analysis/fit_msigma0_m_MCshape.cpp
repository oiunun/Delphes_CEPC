//#include "/besfs5/users/hrqi/fit/Myfitpdf/RooMyRBW/RoodMyRBW_cxx.so"
//#include "/besfs5/users/hrqi/fit/Myfitpdf/RooMyBW/RoodMyBW_cxx.so"
#include "TMath.h"
#include <math.h>
#include <fstream.h>
#include <iomanip.h>
using namespace std;

void fit_msigma0_m_MCshape(){
	gSystem->Load("libRooFit");
	using namespace RooFit;
	gROOT->Reset();
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(0);
	gStyle->SetStripDecimals(kFALSE);
	gStyle->SetPalette(1);
	gStyle->SetOptTitle(0);
	gStyle->SetOptFit(0);

	gStyle->SetPaperSize(20,26);
	gStyle->SetFrameBorderMode(0);
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetPadBorderMode(0);
	gStyle->SetPadColor(0);
	gStyle->SetCanvasColor(0);
	gStyle->SetTitleFillColor(0);
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadBottomMargin(0.15);
	gStyle->SetTitleSize(0.06," ");   // set the pad title size
	gStyle->SetLabelFont(22,"xyz");
	gStyle->SetLabelSize(0.05,"xyz");
	gStyle->SetTextFont(22);  
	gStyle->SetTextSize(0.08);
	gStyle->SetStatFont(22);

	gStyle->SetMarkerStyle(8);
	gStyle->SetHistLineWidth(1.85);
	gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

	RooMsgService::instance().deleteStream(0) ;
	RooMsgService::instance().deleteStream(1) ;

	double x_min=1.12, x_max=1.26;//0.65  0.91
	int Nbins= 40   ;//24
	double cos_min,cos_max;
	ofstream outfile("eps/fit_unbin.dat");  //output nsig  in "data_ang.dat" after fitting

	if{true}{
		TChain *t_data = new TChain ("chi2omg");
		t_data->Add("../data_cut.root");
		Double_t  msigma0_m;
		t_data->SetBranchAddress("msigma0_m",&msigma0_m);

		TFile *newfile = new TFile("./data_combined_tmp.root","recreate");
		TTree *newtree = new TTree("Xi1530","A temp tree");

		newtree->Branch("msigma0_m",&msigma0_m,"msigma0_m/D");

		for(int i=0; i<t_data->GetEntries(); i++){
			t_data->GetEntry(i);
			if(msigma0_m<x_min||msigma0_m>x_max)continue;
			newtree->Fill();
		}

		newtree->Write();
		newfile->Close();
	}//finish cut_opt

	//------fit--------------
	{
		TFile *f_data_t = new TFile("data_combined_tmp.root");
		TTree *tree_t = (TTree*) gDirectory->Get("Xi1530");

    RooRealVar x("msigma0_m","msigma0_m_(GeV/c^2)", x_min, x_max,"GeV/c^{2}");
    RooDataSet data("data","data", tree_t, x);   

    //**************************sig = msigma0_m(mxim1530) MCshape x Gauss**********************/
    TChain *t0 = new TChain("chi2omg");
    t0.Add("/scratchfs/bes/yangrunjia/PublicSources/SSP_neutral/ana/SSKK_cut.root");

    Double_t  msigma0_m;
    t0->SetBranchAddress("msigma0_m",&msigma0_m);

    TFile *newf0 = new TFile("./SD_fit_tmp.root","recreate");
    TTree *newt0 = new TTree("chi2omg","0A temp tree");

    newt0->Branch("msigma0_m",&msigma0_m,"msigma0_m/D");
    for(int i=0; i<t0->GetEntries(); i++){
        t0->GetEntry(i);
        if(msigma0_m<x_min||msigma0_m>x_max)continue;
        newt0->Fill();
    }
    newt0->Write();

    RooDataSet exc1x("msigma0_m","msigma0_m", newt0, x);
    //0819		RooKeysPdf bkg1x ("bkg1x", "bkg1x", x, exc1x, RooKeysPdf::MirrorBoth, 1);
    RooKeysPdf bkg1x ("bkg1x", "bkg1x", x, exc1x, 0);


    //**************************sig = msigma0_m(mxim1530) MCshape x Gauss**********************/
    TChain *t0x = new TChain("chi2omg");
    t0x.Add("/scratchfs/bes/yangrunjia/PublicSources/SSP_neutral/ana/GGLLbPhi_cut.root");

    Double_t  msigma0_m;
    t0x->SetBranchAddress("msigma0_m",&msigma0_m);

    TFile *newf0x = new TFile("./SDxx_fit_tmp.root","recreate");
    TTree *newt0x = new TTree("chi2omg","0A temp tree");

    newt0x->Branch("msigma0_m",&msigma0_m,"msigma0_m/D");
    for(int i=0; i<t0x->GetEntries(); i++){
        t0x->GetEntry(i);
        if(msigma0_m<x_min||msigma0_m>x_max)continue;
        newt0x->Fill();
    }
    newt0x->Write();

    RooDataSet exc0x("msigma0_m","msigma0_m", newt0x, x);
    //0819		RooKeysPdf bkg1x ("bkg1x", "bkg1x", x, exc1x, RooKeysPdf::MirrorBoth, 1);
    RooKeysPdf bkg0x ("bkg0x", "bkg0x", x, exc0x, 0);




    //************************** Double Delta and Omega **********************/
    TChain *tbg = new TChain("chi2omg");
    tbg->Add("/scratchfs/bes/jiaojunkun/psi3686_SSphi/exMC/cut/all/bkg.root");

    Double_t  msigma0_m;
    tbg->SetBranchAddress("msigma0_m",&msigma0_m);

    TFile *newf1 = new TFile("./LLP_chicJ_fit_tmp.root","recreate");
    TTree *newt1 = new TTree("chi2omg","0A temp tree 1");

    newt1->Branch("msigma0_m",&msigma0_m,"msigma0_m/D");
    for(int i=0; i<tbg->GetEntries(); i++){
        tbg->GetEntry(i);
        if(msigma0_m<x_min||msigma0_m>x_max)continue;
        newt1->Fill();
    }
    newt1->Write();

    RooDataSet exc2x("msigma0_m","msigma0_m", newt1, x);
    //0819          RooKeysPdf bkg2x ("bkg2x", "bkg2x", x, exc2x, RooKeysPdf::MirrorBoth, 1);
    RooKeysPdf bkg2x ("bkg2x", "bkg2x", x, exc2x, 0);



    TFile *f1=new TFile("/scratchfs/bes/yangrunjia/PublicSources/SSP_neutral/ana/signal_cut.root");
    TTree *t1 = (TTree*)gDirectory->Get("chi2omg");

    Double_t  msigma0_m;
    t1->SetBranchAddress("msigma0_m",&msigma0_m);

    TFile *newf1 = new TFile("./signal_combined_tmp.root","recreate");
    TTree *newt1 = new TTree("chi2omg","A temp tree");

    newt1->Branch("msigma0_m",&msigma0_m,"msigma0_m/D");
    for(int i=0; i<t1->GetEntries(); i++){
        t1->GetEntry(i);
        if(msigma0_m<x_min||msigma0_m>x_max)continue;
        newt1->Fill();
    }
    newt1->Write();

    // msigma0_m  MCshape
    RooDataSet data1x("msigma0_m","msigma0_m", newt1, x);
    RooKeysPdf sig1x("sig1x", "sig1x", x, data1x, RooKeysPdf::MirrorBoth, 2);
    //                RooKeysPdf sig1x("sig1x", "sig1x", x, data1x, 0);


    x.setBins(10000, "cache_MC");

    //msigma0_m  Gauss
    
    
    RooRealVar GaussM("GaussM","m",2.86939e-03   , -2e-2,  2e-2  );
    RooRealVar GaussS("GaussS","sigma",  5.65007e-04   ,1e-5 ,5e-2);
    RooGaussian Gauss("Gauss","gauss", x, GaussM, GaussS);

    //msigma0_m  MCshape x Gauss
    RooFFTConvPdf sig("sig", "sig", x, sig1x, Gauss);
    //msigma0_m  Gauss
    //	****************************double gauss******************
    /*		RooRealVar  mean1( "mean1" , "mean1" , 2.17595e-03  ,   -1, 1);//函数中得变量，名字，标题，范围，和可选得初值
          RooRealVar  sigma1("sigma1", "sigma1",  0.0008,-1e-2,1e-2);
          RooGaussian gauss1("gauss1", "gauss1", x, mean1, sigma1 );//高斯
    //	*******************************Breit Winger******************************************
    RooRealVar mean("mean","mean",0.783,0.777,0.787);//0.782,0.777,0.787
    RooRealVar width("width","width",0.00849);
    RooBreitWigner BW("BW","BW PDF",x,mean,width);
    RooFFTConvPdf sig("sig", "sig", x,BW, gauss1);
    */

    //*********************************************************/
    RooRealVar p0("p0", "poly 0",  0 ,   -3000.,  3000. );
    RooRealVar p1("p1", "poly 1",   0 ,  -10.,  10. );
    RooRealVar p2("p2", "poly 2",  0.,  -10.,  10. );
    RooRealVar p3("p3", "poly 3",  0.,  -10.,  10. );
    RooChebychev  poly("poly","poly PDF", x, RooArgList(p0,p1));
    //		RooPolynomial  poly("poly","poly PDF", x, RooArgList(p0,p1)  );


    double Nmax=tree_t->GetEntries();
    RooRealVar nsig("N_{sig}","#sig eventas", 0.35*Nmax ,  0, Nmax);//, 0.0, 400000.0);
		//RooRealVar nsig("N_{sig}","#sig events", 0);//, 0.0, 400000.0);
    RooRealVar nbg("N_{poly}","#bg events",   0.43*Nmax, 0, Nmax);//others100000
    RooRealVar nSD("N_{SD}","#gamma eta",  158.4*0.985 );//1242
    RooRealVar nbkg2("N_{chicJ}","#bg events", 70.1*0.866   );//others100000
    RooRealVar nbkg3("N_{sigma0SD}","#bg events",  232);// 0.43*Nmax, 0, Nmax  );//others100000


    RooAddPdf sum("sum","gaus+poly",  RooArgList(sig,poly,bkg0x,bkg1x,bkg2x), RooArgList(nsig, nbg,nSD,nbkg2,nbkg3));
    //		RooAddPdf sum("sum","gaus+poly",  RooArgList(sig,bkg1x), RooArgList(nsig,nSD));

    x.setRange("signal",x_min,x_max);

    //1018	  sum.fitTo(data, "ser",Range("signal"));
    //		RooFitResult *fit_res = sum.fitTo(data, Extended(), Minos(1), Save());
    RooFitResult *fit_res = sum.fitTo(data, Extended(), Save());
    RooPlot* ma0frame = x.frame();

    data.plotOn(ma0frame,Binning(Nbins),LineWidth(3),LineColor(1),MarkerStyle(20-19),MarkerColor(1),MarkerSize(1.3),RooFit::Name("data"));
    data.plotOn(ma0frame,Binning(Nbins),LineWidth(3));
    sum.plotOn(ma0frame,  LineColor(Liver),RooFit::Name("sum"),LineWidth(5));
    sum.plotOn(ma0frame, Components("bkg0x"),LineColor(3),FillColor(3),DrawOption("FL"),FillStyle(1001),RooFit::Name("SD_sigma0"));
    sum.plotOn(ma0frame, Components("bkg1x"),LineColor(parisgreen),FillColor(parisgreen),DrawOption("FL"),RooFit::Name("SD_phi"));
    sum.plotOn(ma0frame, Components("poly"),LineColor(gulan),LineStyle(9),RooFit::Name("poly"),LineWidth(4));
    sum.plotOn(ma0frame, Components("sig"),LineColor(IntOrg),LineStyle(9),RooFit::Name("sig"),LineWidth(4));//signal
    sum.plotOn(ma0frame, Components("bkg2x"),LineColor(suxin),FillColor(suxin),DrawOption("L"),LineStyle(7),LineWidth(4),RooFit::Name("chicJ_LLP"));

    //sum.paramOn(ma0frame,Layout(0.565,0.96,0.97),Format("NEU", AutoPrecision(2)));





    //=============


    x.setRange("signal",0.783-0.022,0.783+0.022);
    RooAbsReal *integral = poly.createIntegral(x,NormSet(x),Range("signal"));
    double normalizedIntegralValue = integral->getVal();

    cout << "normalizedIntegralValue : " << normalizedIntegralValue << endl;

    x.setRange("sideband_1", 0.673,0.717 );
    RooAbsReal *integral = poly.createIntegral(x,NormSet(x),Range("sideband_1"));
    double normalizedIntegralValue_sideband_left = integral->getVal();

    cout << "normalizedIntegralValue_sideband_left : " << normalizedIntegralValue_sideband_left << endl;

    x.setRange("sideband_2", 0.849, 0.893 );
    // x.setRange("sideband_2",0.65 , 0.91 );
    RooAbsReal *integral = poly.createIntegral(x,NormSet(x),Range("sideband_2"));
    double normalizedIntegralValue_sideband_right = integral->getVal();

    cout << "normalizedIntegralValue_sideband_right : " << normalizedIntegralValue_sideband_right << endl;


    cout << "Total SD : " << normalizedIntegralValue_sideband_right+normalizedIntegralValue_sideband_left << endl;



    //=============



    TCanvas *c_msigma0_m = new TCanvas("msigma0_m", "msigma0_m",0,0,800,600);
    //	  c_msigma0_m->SetLogy();	  
    c_msigma0_m->SetFillColor(0);
    c_msigma0_m->Divide(1,1);
    c_msigma0_m->cd(1);
    pady=c_msigma0_m->cd(0);
    pady->SetLeftMargin(0.15);
    pady->SetBottomMargin(0.15);
    pady->SetRightMargin(0.06);
    pady->SetTopMargin(0.03);
    ma0frame->GetXaxis()->SetTitleOffset(0.8);
    ma0frame->GetXaxis()->SetTitleSize(0.08);
    ma0frame->GetXaxis()->CenterTitle();
    ma0frame->GetYaxis()->SetTitleOffset(0.9);
    ma0frame->GetYaxis()->SetTitleSize(0.065);
    //ma0frame->GetYaxis()->SetRangeUser(0,70);
    ma0frame->SetTitle("#bar#Xi(1530)^{+} Fit_data");
    ma0frame->SetXTitle("RM(#Sigma^{0}#phi) [GeV/c^{2}]");	  
    ma0frame->SetYTitle("Events/ 3.5 (MeV/c^{2})");	  
    ma0frame->GetXaxis()->SetTitleFont(22);
    ma0frame->GetYaxis()->CenterTitle();
    ma0frame->GetYaxis()->SetTitleFont(22);
    ma0frame->GetXaxis()->SetLabelFont(22);
    ma0frame->GetYaxis()->SetLabelFont(22);
    ma0frame->Draw();




    Int_t nParsToFit = 0;
    Double_t chi2_ndf=-1; Double_t chi2 = 99999.;
    chi2 = ma0frame->chiSquare("Model_1", "Data", 9);//9
    chi2_ndf = ma0frame->chiSquare(nParsToFit);//拟合优度  
    cout<<"chi2_ndf             "<<chi2_ndf<<endl;



    //		leg = new TLegend(0.7, 0.7, 0.92, 0.92);	  
    leg = new TLegend(0.65, 0.65, 0.92, 0.92);

    leg->AddEntry(ma0frame->findObject("data"), "Data", "lep");
    leg->AddEntry(ma0frame->findObject("sum"), "Total Fitting", "l");
    leg->AddEntry(ma0frame->findObject("sig"), "Signal", "l");
    leg->AddEntry(ma0frame->findObject("poly"), "Smooth backgrounds.", "l");
    leg->AddEntry(ma0frame->findObject("chicJ_LLP"), "#chi_{cJ} #rightarrow #Lambda#bar{#Lambda}#phi", "L");
    leg->AddEntry(ma0frame->findObject("SD_phi"), "#phi sideband", "FL");
    leg->AddEntry(ma0frame->findObject("SD_sigma0"), "#Sigma^{0} Sideband", "fl");

    leg->SetFillColor(0);
    //		leg->SetLineColor(kWhite);
    leg->Draw();

    TLatex lt;
    lt.SetNDC();
    lt.SetTextAngle(0);
    lt.SetTextSize(0.05);
    lt.SetTextColor(2);
    lt.DrawLatex(0.19, 0.90, Form("N_{sig} = %4.1f #pm %4.1f", nsig.getVal(), nsig.getError()));
    //0819		lt.DrawLatex(0.19, 0.82, Form("#Delta M = %4.1f #pm %4.1f", 1000*GaussM.getVal(), 1000*GaussM.getError()));

    string root_name;
    root_name="memo_eps/msigma0_m_MCshape.eps";
    c_msigma0_m->Print(root_name.c_str());

    cout << "total events in this cosThexiRecoil("<<cos_min<<","<<cos_max<<") rang in data: "<<Nmax<<endl;
    outfile << "nsig: "<<nsig.getVal()<<" +- "<<nsig.getError()<<endl;
    outfile << "nbg: "<<nbg.getVal()<<" +- "<<nbg.getError()<<endl;


  }
  cout << "job done!"<< endl;
}
