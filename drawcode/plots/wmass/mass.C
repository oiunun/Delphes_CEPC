{
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   gStyle->SetCanvasBorderMode(0);
   gStyle->SetCanvasBorderSize(0);
   gStyle->SetCanvasColor(10);
   gStyle->SetLabelFont(22,"xyz");
   gStyle->SetLabelSize(0.06,"xyz");
   gStyle->SetLabelOffset(0.04,"xyz");
   gStyle->SetNdivisions(606,"xyz");
   gStyle->SetTitleFont(22,"xyz");
   gStyle->SetTitleColor(1,"xyz");
   gStyle->SetTitleSize(0.06,"xyz");
   gStyle->SetTitleOffset(1.50,"xyz");
   gStyle->SetTitleX(0.40);
   gStyle->SetPadBorderMode(0);
   gStyle->SetPadBorderSize(0);
   gStyle->SetPadColor(10);
   gStyle->SetPadLeftMargin(0.17);
   gStyle->SetPadBottomMargin(0.17);
   gStyle->SetPadRightMargin(0.05);
   gStyle->SetPadTopMargin(0.05);
   gStyle->SetErrorX(0);
   gStyle->SetLegendBorderSize(0);
   gStyle->SetOptDate(0);
   gStyle->SetStatBorderSize(1);
   gStyle->SetStatColor(10);
   gStyle->SetStatFont(22);
   gStyle->SetStatFontSize(0.04);
   gStyle->SetOptFit(1111);
   gStyle->SetLegendFont(22);

   const Int_t n1 = 15, n2=9, n3=1, n4=2; 
   Double_t exx [n1]  = {    0.,    0.,    0.,   0.,   0.,   0.,    0.,    0.,   0.,   0.,     0,      0,     0,       0,      0};
   Double_t xx1 [n1]  = {1983.3,1983.2,1988.9, 1988.9, 1989, 1989,  1990,  1992, 2000, 2002,  2006,   2007,   2009,   2012,   2018};
   Double_t yy1 [n1]  = {    80,    81,   89.,   81.8, 82.7, 80.0, 80.79, 80.84, 81.4, 80.3, 82.87, 80.413, 80.401, 80.367, 80.520};  
   Double_t ey1 [n1]  = {   10.,    5.,   6.7,    6.5,  2.9,  4.1,  0.89,  0.86,  4.7,  2.6,  1.84,  0.048,  0.043,  0.026,  0.115};

   Double_t xx2 [n2]  = { 2000.5, 2002.3, 2006.3, 2005.9, 2005.5, 2008.4, 2012.3,    2014,   2017};
   Double_t yy2 [n2]  = { 80.433, 80.483, 80.440, 80.270, 80.415, 80.336, 80.387,  80.375, 80.370};  
   Double_t ey2 [n2]  = {  0.079,  0.084,  0.051,  0.056,  0.052,  0.067,  0.019,   0.023, 0.0184};

   Double_t xx3 [n2]  = {    2022};
   Double_t yy3 [n2]  = { 80.4335};  
   Double_t ey3 [n2]  = {  0.0094};

   Double_t xx31[n2]  = {    2021};
   Double_t yy31[n2]  = { 80.3540};  
   Double_t ey31[n2]  = {  0.0316};

   Double_t xx4 [n2]  = {    2035,   2043};
   Double_t yy4 [n2]  = { 80.4335,80.4335};  
   Double_t ey4 [n2]  = {  0.0010, 0.0010};
   TGraphErrors gtd("./band.txt",
         "%lg %lg %lg");
   gtd.SetTitle(
         "SM predition;"
         ";"
         "");
   gtd.SetFillColor(kCyan);

   gr1 = new TGraphErrors(n1,xx1,yy1,exx,ey1);
   gr1->SetTitle("Measuring the W mass ");
   gr1->SetMarkerColor(1);
   gr1->SetMarkerStyle(21);

   gr2 = new TGraphErrors(n2,xx2,yy2,exx,ey2);
   gr2->SetTitle("Measuring the W mass ");
   gr2->SetMarkerColor(6);
   gr2->SetMarkerStyle(22);

   gr3 = new TGraphErrors(n3,xx3,yy3,exx,ey3);
   gr3->SetTitle("Measuring the W mass ");
   gr3->SetMarkerColor(3);
   gr3->SetMarkerStyle(20);

   gr31 = new TGraphErrors(n3,xx31,yy31,exx,ey31);
   gr31->SetTitle("Measuring the W mass ");
   gr31->SetMarkerColor(4);
   gr31->SetMarkerStyle(20);

   gr4 = new TGraphErrors(n4,xx4,yy4,exx,ey4);
   gr4->SetTitle("Measuring the W mass ");
   gr4->SetMarkerColor(2);
   gr4->SetMarkerStyle(20);


   TH2D *frm1 = new TH2D("frame1","History of W mass", 100, 1980, 2025, 100, 71.00, 99.00);
   frm1->GetXaxis()->SetTitle("Publication Date");
   frm1->GetYaxis()->SetTitle("W Mass (GeV/c^{2})");
   frm1->GetXaxis()->CenterTitle();
   frm1->GetYaxis()->CenterTitle();
   frm1->GetXaxis()->SetLabelFont(22);
   frm1->GetYaxis()->SetLabelFont(22);
   frm1->GetXaxis()->SetTitleFont(22);
   frm1->GetYaxis()->SetTitleFont(22);

   c1 = new TCanvas("c1","History of W mass measurement",20,10,900,500);
   c1->SetFillColor(0);
   c1->GetFrame()->SetFillColor(21);
   c1->GetFrame()->SetBorderSize(12);
   frm1->Draw("axis");
   gtd.DrawClone("SE3L"); 
   gr1->Draw("P");
   gr2->Draw("P");
   gr3->Draw("P");
   gr31->Draw("P");


   TLegend leg(.70,.73,.94,.93,"");
   leg.SetFillColor(0);
   leg.AddEntry(gr1,"Not used by PDG");
   leg.AddEntry(gr2,"Used by PDG");
   leg.AddEntry(gr3,"CDFII");
   leg.AddEntry(gr31,"LHCb");
   leg.AddEntry(&gtd,"SM prediction");
   leg.DrawClone("Same");


   c1->SaveAs("all.png"); 


   TH2D *frm2 = new TH2D("frame2","History of W mass", 100, 1999, 2045, 100, 80.21, 80.59);
   frm2->GetXaxis()->SetTitle("Publication Date");
   frm2->GetYaxis()->SetTitle("W Mass (GeV/c^{2})");
   frm2->GetXaxis()->CenterTitle();
   frm2->GetYaxis()->CenterTitle();
   frm2->GetXaxis()->SetLabelFont(22);
   frm2->GetYaxis()->SetLabelFont(22);
   frm2->GetXaxis()->SetTitleFont(22);
   frm2->GetYaxis()->SetTitleFont(22);

   c2 = new TCanvas("c2","History of W mass measurement",20,10,900,500);
   c2->SetFillColor(0);
   c2->GetFrame()->SetFillColor(21);
   c2->GetFrame()->SetBorderSize(12);
   frm2->Draw("axis");
   gtd.DrawClone("SE3L"); 
   gr2->Draw("P");
   gr3->Draw("P");
   gr31->Draw("P");
   gr4->Draw("P");


   TLegend leg2(.70,.73,.94,.93,"");
   leg2.SetFillColor(0);
   leg2.AddEntry(gr2,"Used by PDG");
   leg2.AddEntry(gr3,"CDF II");
   leg2.AddEntry(gr31,"LHCb");
   leg2.AddEntry(gr4,"CEPC/FCC-ee");
   leg2.AddEntry(&gtd,"SM prediction");
   leg2.DrawClone("Same");

   c2->SaveAs("zom.png"); 

   return c1;
}
