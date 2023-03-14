/**********************************************************
 *                                                        *
 *         BES III Plotstyle: example 2D data plot        *
 *                                                        *
 *         August 2009, Niklaus Berger                    *
 *         nberger@ihep.ac.cn                             *
 *                                                        *
 *********************************************************/
#include "cepcplotstyle.C"

void Data2DPlot() {
  // Set general style options
  SetStyle();
  // Set options for "final" plots
  SetPrelimStyle();
  // OR: create meeting style plots with stat and fitbox
  //  SetMeetingStyle();

  // Create some dummy histograms - note that you have to call
  // SetStyle() BEFORE you create the histograms
  TH2F * datahist = new TH2F("datahist", "Data", 40, 0, 4, 40, -1, 3);
  datahist->FillRandom("gaus",3000);
  
  /* TH1F * mc1hist = new TH1F("mc1hist", "Some Monte Carlo", 40, 0, 4);
  mc1hist->FillRandom("landau",10000);
  
  TH1F * mc2hist = new TH1F("mc2hist", "Some other Monte Carlo", 40, 0, 4);
  mc2hist->FillRandom("gaus",7000);*/

  // Name the axes of at least the data histogram
  NameAxes(datahist, "p_{#pi^{+}} (GeV/c)", "p_{#pi^{-}} (GeV/c)");
  
  PlotScatter("figs/nicescatter",
	      datahist);



}
