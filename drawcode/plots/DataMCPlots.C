/**********************************************************
 *                                                        *
 *         BES III Plotstyle: example data MC plot        *
 *                                                        *
 *         August 2009, Niklaus Berger                    *
 *         nberger@ihep.ac.cn                             *
 *                                                        *
 *********************************************************/
#include "cepcplotstyle.C"
void DataMCPlots() {
  // Set general style options
  SetStyle();
  // Set options for "final" plots
  SetPrelimStyle();
  // OR: create meeting style plots with stat and fitbox
  //  SetMeetingStyle();

  // Create some dummy histograms - note that you have to call
  // SetStyle() BEFORE you create the histograms
  TH1F * datahist = new TH1F("datahist", "Data:  CUT 2.0 3.0", 40, 0, 4);
  datahist->FillRandom("gaus",10000);
  
  TH1F * mc1hist = new TH1F("mc1hist", "Some Monte Carlo", 40, 0, 4);
  mc1hist->FillRandom("landau",10000);
  
  TH1F * mc2hist = new TH1F("mc2hist", "Some other Monte Carlo", 40, 0, 4);
  mc2hist->FillRandom("gaus",7000);

  // Name the axes of at least the data histogram
  NameAxes(datahist, "p_{#pi} [GeV]", "Events / 0.1 GeV");
  
  PlotDataMC("figs/niceplot",
	     datahist, "CEPC ZH Data",
	     mc1hist, "Some Monte Carlo",
	     mc2hist, "More Monte Carlo",
		  0,"",
		  0,"", 
		  true,
		  false,	
		  true,
		  "Z#rightarrow #mu^{+}#mu^{-}, H#rightarrow b#bar{b}; #int L = 5 ab^{-1}"
		  );
  



}
