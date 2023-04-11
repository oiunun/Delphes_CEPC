#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH1F.h"
#include "TMath.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TLine.h"
using namespace std;
void detector(){
  TCanvas *c=new TCanvas("c","c",1300,700);
  c->Range(-0.5,-0.3,4 ,2.2);
  TGaxis *x = new TGaxis(-0.1,-0.1,3.5,-0.1,-0.1,3.5,510,"S");
  x->SetTitle("z(m)");
  x->SetLabelSize(0.03);
  x->SetTickSize(0.01);
  x->Draw();
  TGaxis *y = new TGaxis(-0.1,-0.1,-0.1,2,-0.1,2,510,"S");
  y->SetTitle("R(m)");
  y->SetLabelSize(0.03);
  y->SetTickSize(0.01);
  y->Draw();


  TLegend *le=new TLegend(0.82,0.5 , 0.92, 0.9);
  TLine *VTX1A = new TLine(0,0.016,0.2,0.016);
  TLine *VTX1B = new TLine(0,0.018,0.2,0.018);
  TLine *VTX2A = new TLine(0,0.038,0.2,0.038);
  TLine *VTX2B = new TLine(0,0.040,0.2,0.040);
  TLine *VTX3A = new TLine(0,0.058,0.2,0.058);
  TLine *VTX3B = new TLine(0,0.060,0.2,0.060);
  VTX1A->SetLineWidth(2);
  VTX1B->SetLineWidth(2);
  VTX2A->SetLineWidth(2);
  VTX2B->SetLineWidth(2);
  VTX3A->SetLineWidth(2);
  VTX3B->SetLineWidth(2);

  le->AddEntry(VTX1A,"VERTEX","l");

  TLine *DC1 = new TLine(0,0.605,2.98,0.605);
  TLine *DC2 = new TLine(0,1.79,2.98,1.79);
  TLine *DC3 = new TLine(0,0.605,0,1.79);
  TLine *DC4 = new TLine(2.98,0.605,2.98,1.79);
  DC1->SetLineColor(kMagenta+2);
  DC2->SetLineColor(kMagenta+2);
  DC3->SetLineColor(kMagenta+2);
  DC4->SetLineColor(kMagenta+2);
  DC1->SetLineWidth(2);
  DC2->SetLineWidth(2);
  DC3->SetLineWidth(2);
  DC4->SetLineWidth(2);

  le->AddEntry(DC1,"DC","l");

  TLine *SET  = new TLine(0,1.815,2.98,1.815);
  SET->SetLineColor(kOrange-3);
  SET->SetLineWidth(2);
  le->AddEntry(SET,"SET","l");

  TLine *ETD  = new TLine(3.0,0.6,3.0,1.822);
  ETD->SetLineColor(kCyan);
  ETD->SetLineWidth(2);
  le->AddEntry(ETD,"ETD","l");

  TLine *SIT01 = new TLine(0,0.12,0.241,0.12);
  TLine *SIT02 = new TLine(0,0.27,0.455,0.27);
  TLine *SIT03 = new TLine(0,0.42,0.721,0.42);
  TLine *SIT04 = new TLine(0,0.57,0.988,0.57);
  SIT01->SetLineColor(kBlue);
  SIT02->SetLineColor(kBlue);
  SIT03->SetLineColor(kBlue);
  SIT04->SetLineColor(kBlue);
  SIT01->SetLineWidth(2);
  SIT02->SetLineWidth(2);
  SIT03->SetLineWidth(2);
  SIT04->SetLineWidth(2);
  le->AddEntry(SIT01,"SIT","l");

  TLine *DSK1 = new TLine(0.241,0.0295,0.241,0.12);
  TLine *DSK2 = new TLine(0.455,0.0305,0.455,0.27);
  TLine *DSK3 = new TLine(0.721,0.0325,0.721,0.42);
  TLine *DSK4 = new TLine(0.988,0.0340,0.988,0.57);
  DSK1->SetLineColor(kRed);
  DSK2->SetLineColor(kRed);
  DSK3->SetLineColor(kRed);
  DSK4->SetLineColor(kRed);
  DSK1->SetLineWidth(2);
  DSK2->SetLineWidth(2);
  DSK3->SetLineWidth(2);
  DSK4->SetLineWidth(2);
  le->AddEntry(DSK1,"DSK","l");

  TLine *cos1 = new TLine(0,0,0.408,2);
  TLine *cos2 = new TLine(0,0,0.8729,2);
  TLine *cos3 = new TLine(0,0,1.5,2);
  TLine *cos4 = new TLine(0,0,3.2828,2);
  TLine *cos5 = new TLine(0,0,3.5,0.7107);

  cos2->SetLineColor(kYellow+2);
  cos3->SetLineColor(kBlue);
  cos4->SetLineColor(kRed);
  cos5->SetLineColor(kMagenta+2);
  cos1->SetLineStyle(2);
  cos2->SetLineStyle(2);
  cos3->SetLineStyle(2);
  cos4->SetLineStyle(2);
  cos5->SetLineStyle(2);
  cos1->SetLineWidth(2);
  cos2->SetLineWidth(2);
  cos3->SetLineWidth(2);
  cos4->SetLineWidth(2);
  cos5->SetLineWidth(2);

  le->AddEntry(cos1,"cos#theta = 0.2","l");
  le->AddEntry(cos2,"cos#theta = 0.4","l");
  le->AddEntry(cos3,"cos#theta = 0.6","l");
  le->AddEntry(cos4,"cos#theta = 0.854","l");
  le->AddEntry(cos5,"cos#theta = 0.98","l");


  VTX1A->Draw();
  VTX1A->Draw();
  VTX1B->Draw();
  VTX2A->Draw();
  VTX2B->Draw();
  VTX3A->Draw();
  VTX3B->Draw();
  
  SIT01->Draw();
  SIT02->Draw();
  SIT03->Draw();
  SIT04->Draw();

  DC1->Draw();
  DC2->Draw();  
  DC3->Draw();
  DC4->Draw();

  SET->Draw();

  ETD->Draw();

  DSK1->Draw();
  DSK2->Draw();
  DSK3->Draw();
  DSK4->Draw();

  cos1->Draw();
  cos2->Draw();
  cos3->Draw();
  cos4->Draw();
  cos5->Draw();

  le->Draw();


  // g->Draw();
  c->Print("detector/detector_28.png");
  c->Print("detector/detector.pdf");
}