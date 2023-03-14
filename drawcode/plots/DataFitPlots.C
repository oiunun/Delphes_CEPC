/**********************************************************
 *                                                        *
 *         CEPC   Plotstyle: example data/fit plot        *
 *                                                        *
 *                                                        *
 *********************************************************/

#include "cepcplotstyle.C"
void DataFitPlots() {
  // Set general style options
  SetStyle();
  // Set options for "final" plots
  SetPrelimStyle();
  // OR: create meeting style plots with stat and fitbox
  //  SetMeetingStyle();

  // Create a dummy histogram - note that you have to call
  // SetStyle() BEFORE you create the histograms
  TH1D * datahist = new TH1D("datahist", "Data", 40, 0, 4);
  datahist->FillRandom("gaus",10000);
  
  //Do some fits - do not forget the "+" option
  datahist->Fit("gaus","+");
  datahist->Fit("landau","+");
  datahist->Fit("pol3","+");


  // Name the axes of the data histogram
  NameAxes(datahist, "p_{#pi} [GeV]", "Events / 0.1 GeV");
  

  // Names for the fits
  char * names[] = {"Good Fit",
		    "Fit 2",
		    "Fit 3",
		    "Bad Fit",
		    "Old Fit",
		    "Another Fit"};


  // Do the plot
  PlotDataFit("figs/nicefit", datahist, "CEPC Simulation", names, true,
		  0.65,                // Left edge of legend   (fraction of canavas width)
		  0.50,                 // Bottom edge of legend (fraction of canavas height)
		  0.90,                // Right edge of legend  (fraction of canavas width)
		  0.70,
		  "Z#rightarrow #mu^{+}#mu^{-}, H#rightarrow b#bar{b}; #int Ldt = 5 ab^{-1}"
		  );

}
