/**********************************************************
 *                                                        *
 *         CEPC Plotstyle: format functions               *
 *                                                        *
 *         October 2014, Gang LI                          *
 *         li.gang@ihep.ac.cn                             *
 *                                                        *
 *********************************************************/

#ifndef CEPCPLOT__H
#define CEPCPLOT__H
#include <TH1.h>
#include <TH2.h>
#include <TH1F.h>
#include <TF1.h>
#include <THStack.h>
#include <TLatex.h>
#include <TList.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TArc.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TArrow.h>
#include <TGraphErrors.h>
#include <RooPlot.h>
#include <RooAbsPdf.h>
#include <RooRealVar.h>
#include <RooAbsData.h>

#include <iostream>
#include <algorithm>
#include <vector>
using namespace std;
using namespace RooFit;

class TH1;
class TH2;
class TH1F;
class TH1D;
class TAxis;
class TGraph;
class RooPlot;

const Double_t makersize = 1.2;
const Double_t h_scale = 1.5;
// Format for data points
void FormatHist(TH1 *datahist, int n = 0);
void FormatData(TH1 *datahist);
void FormatData(TH2 *datahist);
void FormatData(TGraph *datagraph);
// Format Axis
void FormatAxis(TAxis *axis, double offset = 1.1);
// Format for main MC (red line)
void FormatMC1(TH1 *mc1hist);
void FormatMC2(TH1 *mc1hist);
void FormatMC3(TH1 *mc1hist);
void FormatMC4(TH1 *mc1hist);
// Graph Format for second MC or background
void FormatMC1(TGraph *mc1hist);
void FormatMC2(TGraph *mc1hist);
void FormatMC3(TGraph *mc1hist);
void FormatRooFit(RooPlot *rooplot, Char_t *xtitle = (Char_t *)"M(GeV/c^{2})", Char_t *ytitle = (Char_t *)"Entries", double yoffset = 1.2);
void PlotDataRooFit(Char_t *figfilename,
					RooAbsData *data, RooRealVar x,
					vector<RooAbsPdf *> pdfList,
					char *xtitle = (char *)"M(GeV/c^{2})",
					char *ytitle = (char *)"Entries", bool prelim = false,
					char *channel = (char *)"ZH, Z#rightarrow #mu^{+}#mu^{-}; #int L=5 ab^{-1}",
					double maxy = -1);
//
// Name histogram axes
void NameAxes(TH1 *datahist, const char *xname, const char *yname);
void NameAxes(TGraphErrors *datahist, char *xname, char *yname);

void WriteCEPC();
void WritePreliminary();																		// to be used together with WriteCEPC()
void WriteChannel(char *channel = (char *)"ZH, Z#rightarrow #mu^{+}#mu^{-}; #int L=5 ab^{-1}"); // to be used together with WriteCEPC()

// Make a legend;
// position will have to change depending on the data shape
void MakeLegend(TH1 *datahist,		// Histogram with data
				char *dataname,		// Description of data
				TH1 *mc1hist = 0,	// Histogram with first MC
				char *mc1name = 0,	// Description of first MC
				TH1 *mc2hist = 0,	// Histogram with 2nd MC/BG
				char *mc2name = 0,	// Description of second MC/BG
				double xlow = 0.20, // Left edge of legend   (fraction of canavas width)
				double ylow = 0.75, // Bottom edge of legend (fraction of canavas height)
				double xhi = 0.50,	// Right edge of legend  (fraction of canavas width)
				double yhi = 0.94); // Top edge of legend    (fraction of canavas height)

void MakeLegend(vector<TH1D *> datahist, // Histogram with data
				vector<char *> dataname, // Description of data
				double xlow = 0.20,		 // Left edge of legend   (fraction of canavas width)
				double ylow = 0.75,		 // Bottom edge of legend (fraction of canavas height)
				double xhi = 0.50,		 // Right edge of legend  (fraction of canavas width)
				double yhi = 0.94);		 // Top edge of legend    (fraction of canavas height)

// position will have to change depending on the data shape
void MakeLegend4(TH1 *datahist,		 // Histogram with data
				 char *dataname,	 // Description of data
				 TH1 *mc1hist = 0,	 // Histogram with first MC
				 char *mc1name = 0,	 // Description of first MC
				 TH1 *mc2hist = 0,	 // Histogram with 2nd MC/BG
				 char *mc2name = 0,	 // Description of second MC/BG
				 TH1 *mc3hist = 0,	 // Histogram with 2nd MC/BG
				 char *mc3name = 0,	 // Description of second MC/BG
				 double xlow = 0.20, // Left edge of legend   (fraction of canavas width)
				 double ylow = 0.75, // Bottom edge of legend (fraction of canavas height)
				 double xhi = 0.50,	 // Right edge of legend  (fraction of canavas width)
				 double yhi = 0.94); // Top edge of legend    (fraction of canavas height)

void MakeLegend5(
	TH1 *datahist,		// Histogram with data
	char *dataname,		// Description of data
	TH1 *mc1hist = 0,	// Histogram with first MC
	char *mc1name = 0,	// Description of first MC
	TH1 *mc2hist = 0,	// Histogram with 2nd MC/BG
	char *mc2name = 0,	// Description of second MC/BG
	TH1 *mc3hist = 0,	// Histogram with 2nd MC/BG
	char *mc3name = 0,	// Description of second MC/BG
	TH1 *mc4hist = 0,	// Histogram with 2nd MC/BG
	char *mc4name = 0,	// Description of second MC/BG
	double xlow = 0.20, // Left edge of legend   (fraction of canavas width)
	double ylow = 0.75, // Bottom edge of legend (fraction of canavas height)
	double xhi = 0.50,	// Right edge of legend  (fraction of canavas width)
	double yhi = 0.94); // Top edge of legend    (fraction of canavas height)

// Make a legend;
// position will have to change depending on the data shape
void MakeLegend(TGraph *datahist,	 // Graph with data
				char *dataname,		 // Description of data
				TGraph *mc1hist = 0, // Graph with first MC
				char *mc1name = 0,	 // Description of first MC
				TGraph *mc2hist = 0, // Graph with 2nd MC/BG
				char *mc2name = 0,	 // Description of second MC/BG
				TGraph *mc3hist = 0, // Graph with 3rd MC/BG
				char *mc3name = 0,	 // Description of third MC/BG
				double xlow = 0.20,	 // Left edge of legend   (fraction of canavas width)
				double ylow = 0.5,	 // Bottom edge of legend (fraction of canavas height)
				double xhi = 0.50,	 // Right edge of legend  (fraction of canavas width)
				double yhi = 0.7);	 // Top edge of legend    (fraction of canavas height)

// Make a legend (version for fit functions
// position will have to change depending on the data shape
void MakeLegend(TH1 *datahist,		  // Histogram with data
				char *dataname,		  // Description of data
				char **functionnames, // list of function names
				double xlow = 0.20,	  // Left edge of legend   (fraction of canavas width)
				double ylow = 0.5,	  // Bottom edge of legend (fraction of canavas height)
				double xhi = 0.50,	  // Right edge of legend  (fraction of canavas width)
				double yhi = 0.7);	  // Top edge of legend    (fraction of canavas height)

void SetStyle();		// Set the general style options
void SetPrelimStyle();	// Style options for "final" plots      (no stat/fit box)
void SetMeetingStyle(); // Style options for internal meetings  (stat/fit box)
// Plot a data MC plot
void PlotDataMC(
	char *filename,		 // Name for the output files, without extension
	TH1 *datahist,		 // Histogram with data
	char *dataname,		 // Description of data
	TH1 *mc1hist = 0,	 // Histogram with first MC
	char *mc1name = 0,	 // Description of first MC
	TH1 *mc2hist = 0,	 // Histogram with 2nd MC/BG
	char *mc2name = 0,	 // Description of second MC/BG
	TH1 *mc3hist = 0,	 // Histogram with 2nd MC/BG
	char *mc3name = 0,	 // Description of second MC/BG
	TH1 *mc4hist = 0,	 // Histogram with 2nd MC/BG
	char *mc4name = 0,	 // Description of second MC/BG
	bool prelim = false, // Preliminary plot
	bool logy = false,
	bool cut = false,
	char *channel = (char *)"ZH, Z#rightarrow #mu^{+}#mu^{-}; #int L=5 ab^{-1}");
void PlotDataMC(
	char *filename,			 // Name for the output files, without extension
	vector<TH1D *> datahist, // Histogram with data
	vector<char *> dataname, // Description of data
	bool prelim = false,	 // Preliminary plot
	bool logy = false,
	bool cut = false,
	char *channel = (char *)"ZH, Z#rightarrow #mu^{+}#mu^{-}; #int L=5 ab^{-1}");

// Plot data with one or more (fit) functions
// Functions should be part of the data histograms list of functions
// (i.e. perform fits with the "+" option or add other functions via
// datahist->GetListOfFunctions->Add(TF1 * function))
// functionnames should have at least as many elements as the function
// list
void PlotDataFit(
	char *filename,		   // Name for the output files, without extension
	TH1D *datahist,		   // Histogram with data
	char *dataname,		   // Description of data
	char *functionnames[], // Names of associated functions
	bool prelim = false,   // Preliminary plot
	double xlow = 0.70,	   // Left edge of legend   (fraction of canavas width)
	double ylow = 0.5,	   // Bottom edge of legend (fraction of canavas height)
	double xhi = 0.95,	   // Right edge of legend  (fraction of canavas width)
	double yhi = 0.7,
	char *channel = (char *)"ZH, Z#rightarrow #mu^{+}#mu^{-}; #int L=5 ab^{-1}"); // Top edge of legend    (fraction of canavas height)

// Scatter plot
void PlotScatter(char *filename,	// Name for the output files, without extension
				 TH1 *datahist,		// Histogram with data
				 bool prelim = true // Preliminary plot
);

void PlotScatterCircle(char *filename, // Name for the output files, without extension
					   TH1 *datahist,  // Histogram with data
					   double cx = 0,
					   double cy = 0,
					   double cr = 1,
					   bool prelim = true // Preliminary plot
);
#endif
