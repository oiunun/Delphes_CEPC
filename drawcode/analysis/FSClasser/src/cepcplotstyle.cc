/**********************************************************
 *                                                        *
 *         CEPC  Plotstyle: format functions              *
 *                                                        *
 *         October 2014, Gang LI                          *
 *         li.gang@ihep.ac.cn                             *
 *                                                        *
 *********************************************************/

#include "cepcplotstyle.h"

using namespace std;
// Format for data points
void FormatHist(TH1 *datahist, int n)
{

	//int histstyles[] = {1 , 2, 3, 7, 9,10};
	int markstyles[] = {21, 22, 23, 24, 25, 26};
	int histcolors[] = {1, 2, 3, 4, 5, 6};

	datahist->SetLineColor(histcolors[n]);
	datahist->SetFillColor(histcolors[n]);
	datahist->SetMarkerColor(histcolors[n]);
	datahist->SetMarkerStyle(markstyles[n]);
	datahist->SetLineWidth(3);
	datahist->SetMinimum(0);
	datahist->SetMarkerSize(makersize);
	if (datahist->GetXaxis())
		FormatAxis(datahist->GetXaxis());
	if (datahist->GetYaxis())
		FormatAxis(datahist->GetYaxis(), 1.3);
}

void FormatData(TH1 *datahist)
{
	datahist->SetMarkerStyle(20);
	datahist->SetMarkerSize(makersize);
	datahist->SetMarkerColor(1);
	datahist->SetLineWidth(3);
	datahist->SetLineColor(1);
	datahist->SetMinimum(0);
	if (datahist->GetXaxis())
		FormatAxis(datahist->GetXaxis());
	if (datahist->GetYaxis())
		FormatAxis(datahist->GetYaxis(), 1.3);
}

void FormatData(TH2 *datahist)
{
	datahist->SetMarkerStyle(20);
	datahist->SetMarkerSize(makersize);
	datahist->SetLineWidth(2);
	if (datahist->GetXaxis())
		FormatAxis(datahist->GetXaxis());
	if (datahist->GetYaxis())
		FormatAxis(datahist->GetYaxis(), 1.3);
}

// Format for graph data points
void FormatData(TGraph *datahist)
{
	datahist->SetMarkerStyle(20);
	datahist->SetMarkerSize(makersize);
	datahist->SetLineWidth(2);
}

void FormatAxis(TAxis *axis, double offset)
{
	axis->SetLabelFont(22);
	axis->SetLabelSize(0.06);
	axis->SetLabelOffset(0.01);
	axis->SetNdivisions(404);
	axis->SetTitleFont(22);
	axis->SetTitleColor(1);
	axis->SetTitleSize(0.06);
	axis->SetTitleOffset(offset);
	axis->CenterTitle();
}

void NameAxes(TH1 *datahist, const char *xname, const char *yname)
{
	if (xname)
		datahist->GetXaxis()->SetTitle(xname);
	if (yname)
		datahist->GetYaxis()->SetTitle(yname);
}

void NameAxes(TGraphErrors *datahist, char *xname, char *yname)
{
	if (xname)
		datahist->GetXaxis()->SetTitle(xname);
	if (yname)
		datahist->GetYaxis()->SetTitle(yname);
}

// Format for main MC (red line)
void FormatMC1(TH1 *mc1hist)
{
	mc1hist->SetLineColor(2);
	mc1hist->SetFillColor(0);
	mc1hist->SetLineWidth(2);
	mc1hist->SetFillStyle(1001);
}

// Graph Format for main MC (red line)
void FormatMC1(TGraph *mc1hist)
{
	mc1hist->SetLineColor(2);
	mc1hist->SetLineWidth(2);
}

// Format for second MC or background
// (Blue shaded area)
void FormatMC2(TH1 *mc2hist)
{
	mc2hist->SetLineColor(4);
	mc2hist->SetFillColor(0);
	mc2hist->SetLineWidth(2);
	mc2hist->SetFillStyle(1001);
}

// Format for second MC or background
// (purple shaded area)
void FormatMC3(TH1 *mc3hist)
{
	mc3hist->SetLineColor(6);
	mc3hist->SetFillColor(0);
	mc3hist->SetLineWidth(2);
	mc3hist->SetFillStyle(3001);
}
// Format for second MC or background
// (purple shaded area)
void FormatMC4(TH1 *mc4hist)
{
	mc4hist->SetLineColor(3);
	mc4hist->SetFillColor(0);
	mc4hist->SetLineWidth(2);
	mc4hist->SetFillStyle(3001);
}

// Graph Format for second MC or background
// (Blue line)
void FormatMC2(TGraph *mc2hist)
{
	mc2hist->SetLineColor(2);
	mc2hist->SetLineWidth(2);
}

// Graph Format for third MC or background
// (Blue line)
void FormatMC3(TGraph *mc3hist)
{
	mc3hist->SetLineColor(6);
	mc3hist->SetLineWidth(2);
}
//
void FormatRooFit(RooPlot *xframe, Char_t *xtitle, Char_t *ytitle, double yoffset)
{

	xframe->GetXaxis()->SetTitle(xtitle);
	xframe->GetYaxis()->SetTitle(ytitle);
	xframe->GetXaxis()->SetLabelFont(22);
	xframe->GetYaxis()->SetLabelFont(22);
	xframe->GetXaxis()->CenterTitle();
	xframe->GetYaxis()->CenterTitle();
	xframe->GetXaxis()->SetNdivisions(404);
	xframe->GetYaxis()->SetNdivisions(505);
	//xframe->GetYaxis()->SetRangeUser(0,8500);

	xframe->GetYaxis()->SetLabelSize(0.05);
	xframe->GetYaxis()->SetLabelOffset(0.01);
	xframe->GetYaxis()->SetTitleSize(0.06);
	xframe->GetYaxis()->SetTitleFont(22);
	xframe->GetXaxis()->SetTitleOffset(1.00);
	xframe->GetYaxis()->SetTitleOffset(1.10);
	xframe->GetYaxis()->SetTitleOffset(yoffset);
	xframe->GetXaxis()->SetLabelSize(0.05);
	xframe->GetXaxis()->SetLabelOffset(0.015);
	xframe->GetXaxis()->SetTitleSize(0.06);
	xframe->GetXaxis()->SetTitleFont(22);
}
//
// Write "CEPC" in the upper right corner
void WriteCEPC()
{
	TLatex *CEPC = new TLatex(0.94, 0.94, "CEPC");
	CEPC->SetNDC();
	CEPC->SetTextFont(22);
	CEPC->SetTextSize(0.040);
	CEPC->SetTextAlign(33);
	CEPC->Draw();
}

// Write "Preliminary" below CEPC -
// to be used together with WriteCEPC()
void WritePreliminary()
{
	TLatex *prelim = new TLatex(0.94, 0.94, "CEPC Preliminary");
	prelim->SetNDC();
	prelim->SetTextFont(22);
	prelim->SetTextSize(0.040);
	prelim->SetTextAlign(33);
	prelim->Draw();
}

// to be used together with WriteCEPC()
void WriteChannel(char *channel)
{
	TLatex *prelim = new TLatex(0.94, 0.88, channel);
	prelim->SetNDC();
	prelim->SetTextFont(22);
	prelim->SetTextSize(0.030);
	prelim->SetTextAlign(33);
	prelim->Draw();
}
// Make a legend;
// position will have to change depending on the data shape

void MakeLegend(vector<TH1D *> datahist,
				vector<char *> dataname,
				double xlow,
				double ylow,
				double xhi,
				double yhi)
{
	TLegend *leg = new TLegend(xlow, ylow, xhi, yhi);

	for (unsigned int i = 0; i < datahist.size(); i++)
	{
		if (datahist[i]->Integral() > 0 && dataname[i])
		{
			if (i == 0)
				leg->AddEntry(datahist[i], dataname[i], "LEP");
			else
				leg->AddEntry(datahist[i], dataname[i], "LF");
		}
	}

	leg->SetFillColor(0);
	leg->SetTextFont(22);
	leg->Draw();
	//delete leg;
}

void MakeLegend(TH1 *datahist,	// Histogram with data
				char *dataname, // Description of data
				TH1 *mc1hist,	// Histogram with first MC
				char *mc1name,	// Description of first MC
				TH1 *mc2hist,	// Histogram with 2nd MC/BG
				char *mc2name,	// Description of second MC/BG
				double xlow,	// Left edge of legend
								// (fraction of canavas width)
				double ylow,	// Bottom edge of legend
								// (fraction of canavas height)
				double xhi,		// Right edge of legend
								// (fraction of canavas width)
				double yhi)
{	// Top edge of legend
	// (fraction of canavas height)

	TLegend *leg = new TLegend(xlow, ylow, xhi, yhi);
	if (datahist && dataname)
		leg->AddEntry(datahist, dataname, "LEP");
	if (mc1hist && mc1name)
		leg->AddEntry(mc1hist, mc1name, "L");
	if (mc2hist && mc2name)
		leg->AddEntry(mc2hist, mc2name, "LF");

	leg->SetFillColor(0);
	leg->SetTextFont(22);
	leg->Draw();
}

void MakeLegend4(
	TH1 *datahist,	// Histogram with data
	char *dataname, // Description of data
	TH1 *mc1hist,	// Histogram with first MC
	char *mc1name,	// Description of first MC
	TH1 *mc2hist,	// Histogram with 2nd MC/BG
	char *mc2name,	// Description of second MC/BG
	TH1 *mc3hist,	// Histogram with 2nd MC/BG
	char *mc3name,	// Description of second MC/BG
	double xlow,	// Left edge of legend
					// (fraction of canavas width)
	double ylow,	// Bottom edge of legend
					// (fraction of canavas height)
	double xhi,		// Right edge of legend
					// (fraction of canavas width)
	double yhi)
{	// Top edge of legend
	// (fraction of canavas height)

	TLegend *leg = new TLegend(xlow, ylow, xhi, yhi);
	if (datahist && dataname)
		leg->AddEntry(datahist, dataname, "LEP");
	if (mc1hist && mc1name)
		leg->AddEntry(mc1hist, mc1name, "L");
	if (mc2hist && mc2name)
		leg->AddEntry(mc2hist, mc2name, "LF");
	if (mc3hist && mc3name)
		leg->AddEntry(mc3hist, mc3name, "LF");

	leg->SetFillColor(0);
	leg->SetTextFont(22);
	leg->Draw();
}
void MakeLegend5(
	TH1 *datahist,	// Histogram with data
	char *dataname, // Description of data
	TH1 *mc1hist,	// Histogram with first MC
	char *mc1name,	// Description of first MC
	TH1 *mc2hist,	// Histogram with 2nd MC/BG
	char *mc2name,	// Description of second MC/BG
	TH1 *mc3hist,	// Histogram with 2nd MC/BG
	char *mc3name,	// Description of second MC/BG
	TH1 *mc4hist,	// Histogram with 2nd MC/BG
	char *mc4name,	// Description of second MC/BG
	double xlow,	// Left edge of legend
					// (fraction of canavas width)
	double ylow,	// Bottom edge of legend
					// (fraction of canavas height)
	double xhi,		// Right edge of legend
					// (fraction of canavas width)
	double yhi)
{	// Top edge of legend
	// (fraction of canavas height)

	TLegend *leg = new TLegend(xlow, ylow, xhi, yhi);
	if (datahist && dataname)
		leg->AddEntry(datahist, dataname, "LEP");
	if (mc1hist && mc1name != (char *)"")
		leg->AddEntry(mc1hist, mc1name, "L");
	if (mc2hist && mc2name != (char *)"")
		leg->AddEntry(mc2hist, mc2name, "LF");
	if (mc3hist && mc3name != (char *)"")
		leg->AddEntry(mc3hist, mc3name, "LF");
	if (mc4hist && mc4name != (char *)"")
		leg->AddEntry(mc4hist, mc4name, "LF");

	leg->SetFillColor(0);
	leg->SetTextFont(22);
	leg->Draw();
}

// Make a legend;
// position will have to change depending on the data shape
void MakeLegend(TGraph *datahist, // Graph with data
				char *dataname,	  // Description of data
				TGraph *mc1hist,  // Graph with first MC
				char *mc1name,	  // Description of first MC
				TGraph *mc2hist,  // Graph with 2nd MC/BG
				char *mc2name,	  // Description of second MC/BG
				TGraph *mc3hist,  // Graph with 3rd MC/BG
				char *mc3name,	  // Description of third MC/BG
				double xlow,	  // Left edge of legend
								  // (fraction of canavas width)
				double ylow,	  // Bottom edge of legend
								  // (fraction of canavas height)
				double xhi,		  // Right edge of legend
								  // (fraction of canavas width)
				double yhi)
{	// Top edge of legend
	// (fraction of canavas height)

	TLegend *leg = new TLegend(xlow, ylow, xhi, yhi);
	if (datahist && dataname)
		leg->AddEntry(datahist, dataname, "LEP");
	if (mc1hist && mc1name)
		leg->AddEntry(mc1hist, mc1name, "L");
	if (mc2hist && mc2name)
		leg->AddEntry(mc2hist, mc2name, "L");
	if (mc3hist && mc3name)
		leg->AddEntry(mc3hist, mc3name, "L");

	leg->SetFillColor(0);
	leg->SetTextFont(22);
	leg->Draw();
}
// Make a legend (version for fit functions
// position will have to change depending on the data shape
void MakeLegend(TH1 *datahist,		  // Histogram with data
				char *dataname,		  // Description of data
				char **functionnames, // list of function names
				double xlow,		  // Left edge of legend
				//(fraction of canavas width)
				double ylow, // Bottom edge of legend
				//(fraction of canavas height)
				double xhi, // Right edge of legend
				//(fraction of canavas width)
				double yhi)
{ // Top edge of legend
	//(fraction of canavas height)

	TLegend *leg = new TLegend(xlow, ylow, xhi, yhi);
	if (datahist && dataname)
		leg->AddEntry(datahist, dataname, "LEP");

	TList *list = datahist->GetListOfFunctions();
	unsigned int nfun = list->GetEntries();
	//cout<<nfun<<" functions"<<endl;
	for (unsigned int i = 0; i < nfun; i++)
	{
		TF1 *f1 = (TF1 *)(list->At(i));
		leg->AddEntry(f1, functionnames[i], "L");
	}
	leg->SetFillColor(0);
	leg->SetTextFont(22);
	leg->Draw();
}

// Set the general style options
void SetStyle()
{
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	// No Canvas Border
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetCanvasBorderSize(0);
	// White BG
	gStyle->SetCanvasColor(10);
	// Format for axes
	gStyle->SetLabelFont(22, "xyz");
	gStyle->SetLabelSize(0.10, "xyz");
	gStyle->SetLabelOffset(0.04, "xyz");
	gStyle->SetNdivisions(404, "xyz");
	gStyle->SetTitleFont(22, "xyz");
	gStyle->SetTitleColor(1, "xyz");
	gStyle->SetTitleSize(0.20, "xyz");
	gStyle->SetTitleOffset(1.50, "xyz");
	gStyle->SetTitleX(0.50);
	// No pad borders
	gStyle->SetPadBorderMode(0);
	gStyle->SetPadBorderSize(0);
	// White BG
	gStyle->SetPadColor(10);
	// Margins for labels etc.
	gStyle->SetPadLeftMargin(0.17);
	gStyle->SetPadBottomMargin(0.17);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetPadTopMargin(0.05);
	// No error bars in x direction
	gStyle->SetErrorX(0);
	// Format legend
	gStyle->SetLegendBorderSize(0);
}

// Style options for "final" plots
// (no stat/fit box)
void SetPrelimStyle()
{
	gStyle->SetOptDate(0);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	gStyle->SetOptTitle(1);
	gStyle->SetOptTitle(0);
}

// Style options for internal meetings
// (stat/fit box)
void SetMeetingStyle()
{
	gStyle->SetOptDate(0);
	gStyle->SetOptTitle(1);
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(1001111);
	gStyle->SetStatBorderSize(1);
	gStyle->SetStatColor(10);
	gStyle->SetStatFont(22);
	gStyle->SetStatFontSize(0.03);
	gStyle->SetOptFit(1111);
}

void PlotDataMC(
	char *filename,
	vector<TH1D *> datahist,
	vector<char *> dataname,
	bool prelim,
	bool logy,
	bool cut,
	char *channel)
{

	SetStyle();
	if (prelim)
		SetPrelimStyle();
	else
		SetMeetingStyle();
	//
	int nh = datahist.size();
	TString tname, scut, sname;
	if (strstr(datahist[0]->GetTitle(), "CUT"))
	{
		tname = datahist[0]->GetTitle();
		Int_t ipos = tname.First("CUT");
		scut = tname(ipos, tname.Sizeof() - ipos);
		if (ipos > 3)
			sname = tname.Remove(ipos - 2, tname.Sizeof());
		if (ipos < 3)
			sname = tname.Remove(ipos, tname.Sizeof());
		datahist[0]->SetTitle(sname);
	}
	//
	TCanvas *c1 = new TCanvas("CEPCplots", "CEPC Plots", 800, 800);
	FormatData(datahist[0]);
	if (nh > 1)
		FormatMC1(datahist[1]);
	if (nh > 2)
		FormatMC2(datahist[2]);
	if (nh > 3)
		FormatMC3(datahist[3]);
	if (nh > 4)
		FormatMC4(datahist[4]);
	//
	if (nh > 1)
		datahist[0]->SetMaximum(max(datahist[1]->GetMaximum(), datahist[0]->GetMaximum()) * h_scale);
	if (nh > 2)
		datahist[0]->SetMaximum(max(datahist[2]->GetMaximum(), datahist[0]->GetMaximum()) * h_scale);
	if (nh > 3)
		datahist[0]->SetMaximum(max(datahist[3]->GetMaximum(), datahist[0]->GetMaximum()) * h_scale);
	if (nh > 4)
		datahist[0]->SetMaximum(max(datahist[4]->GetMaximum(), datahist[0]->GetMaximum()) * h_scale);
	//
	if (logy)
	{
		gPad->SetLogy();
		datahist[0]->SetMinimum(0.1);
	}
	//
	if (datahist[0]->Integral() > 0)
		datahist[0]->Draw("E");
	else
		datahist[0]->Draw("axis");
	if (nh > 1)
		datahist[1]->Draw("histsame");
	if (nh > 2)
		datahist[2]->Draw("histsame");
	if (nh > 3)
		datahist[3]->Draw("histsame");
	if (nh > 4)
		datahist[4]->Draw("histsame");
	if (datahist[0]->Integral() > 0)
	{
		datahist[0]->Draw("Esame");
		datahist[0]->Draw("axissame");
	}
	if (prelim)
	{
		//WriteCEPC();
		WritePreliminary();
		WriteChannel(channel);
	}
	TArrow *arr1 = 0, *arr2 = 0;
	if (cut && strstr(scut, "CUT"))
	{
		Char_t cuts[10];
		Float_t x0 = 1, x1 = 1;
		//printf("%s\n",scut.Data());
		sscanf(scut.Data(), "%s %f %f", cuts, &x0, &x1);
		{
			Double_t min0 = datahist[0]->GetBinContent(datahist[0]->FindBin(x0));
			Double_t min1 = datahist[0]->GetBinContent(datahist[0]->FindBin(x1));
			Double_t ymin = max(min0, min1);
			Double_t ymax = datahist[0]->GetMaximum();
			Double_t dely = datahist[0]->GetMaximum() - ymin;
			ymin = ymin + dely * 0.05;
			ymax = (ymax - ymin) * 0.5 + ymin;
			//
			arr1 = new TArrow(x0, ymax, x0, ymin, 0.05, ">");
			arr1->SetLineWidth(5);
			arr1->SetLineColor(kMagenta);
			arr1->Draw("");
			//
			if (x0 != x1)
			{
				arr2 = new TArrow(x1, ymax, x1, ymin, 0.05, ">");
				arr2->SetLineWidth(5);
				arr2->SetLineColor(kMagenta);
				arr2->Draw("");
			}
		}
	}
	if (!logy)
	{
		MakeLegend(datahist, dataname);
	}
	char filenameall[256];
	sprintf(filenameall, "%s.eps", filename);
	//c1->SaveAs(filenameall);
	sprintf(filenameall, "%s.pdf", filename);
	c1->SaveAs(filenameall);
	if (arr1)
		delete arr1;
	if (arr2)
		delete arr2;
	delete c1;
}

// Plot a data MC plot
void PlotDataMC(char *filename, // Name for the output files, without extension
				TH1 *datahist,	// Histogram with data
				char *dataname, // Description of data
				TH1 *mc1hist,	// Histogram with 1st MC
				char *mc1name,	// Description of 1st MC
				TH1 *mc2hist,	// Histogram with 2nd MC/BG
				char *mc2name,	// Description of 2nd MC/BG
				TH1 *mc3hist,	// Histogram with 3rd MC/BG
				char *mc3name,	// Description of 3rd MC/BG
				TH1 *mc4hist,	// Histogram with 3rd MC/BG
				char *mc4name,	// Description of 3rd MC/BG
				bool prelim,	// Preliminary plot
				bool logy,		// Preliminary plot
				bool cut,
				char *channel)
{

	SetStyle();
	if (prelim)
		SetPrelimStyle();
	else
		SetMeetingStyle();
	//
	TString tname, scut, sname;
	//cout<<datahist->GetTitle()<<endl;
	if (strstr(datahist->GetTitle(), "CUT"))
	{
		//tname=(TString)datahist->GetTitle();
		tname = datahist->GetTitle();
		Int_t ipos = tname.First("CUT");
		scut = tname(ipos, tname.Sizeof() - ipos);
		//cout<<ipos<< " "<<scut<<endl;
		//sname = tname.Remove(ipos-2,tname.Sizeof());
		if (ipos > 3)
			sname = tname.Remove(ipos - 2, tname.Sizeof());
		if (ipos < 3)
			sname = tname.Remove(ipos, tname.Sizeof());
		//cout<<ipos<< " "<<tname.Remove(ipos,tname.Sizeof())<<" "<<scut<<endl;
		//cout<<ipos<< " "<<sname<<endl;
		datahist->SetTitle(sname);
		if (mc1hist)
			mc1hist->SetTitle(sname);
		if (mc2hist)
			mc2hist->SetTitle(sname);
		if (mc3hist)
			mc3hist->SetTitle(sname);
		if (mc4hist)
			mc4hist->SetTitle(sname);
	}
	//
	TCanvas *c1 = new TCanvas("CEPCplots", "CEPC Plots", 800, 800);
	FormatData(datahist);
	if (mc1hist)
		FormatMC1(mc1hist);
	if (mc2hist)
		FormatMC2(mc2hist);
	if (mc3hist)
		FormatMC3(mc3hist);
	if (mc4hist)
		FormatMC4(mc4hist);
	//
	if (mc1hist)
		datahist->SetMaximum(max(datahist->GetMaximum(), mc1hist->GetMaximum()) * h_scale);
	if (mc2hist)
		datahist->SetMaximum(max(datahist->GetMaximum(), mc2hist->GetMaximum()) * h_scale);
	if (mc3hist)
		datahist->SetMaximum(max(datahist->GetMaximum(), mc3hist->GetMaximum()) * h_scale);
	if (mc4hist)
		datahist->SetMaximum(max(datahist->GetMaximum(), mc4hist->GetMaximum()) * h_scale);
	//
	if (logy)
	{
		gPad->SetLogy();
		datahist->SetMinimum(0.1);
	}
	//
	if (datahist->Integral() > 0)
		datahist->Draw("E");
	else
		datahist->Draw("axis");
	if (mc4hist)
		mc4hist->Draw("histsame");
	if (mc3hist)
		mc3hist->Draw("histsame");
	if (mc2hist)
		mc2hist->Draw("histsame");
	if (mc1hist)
		mc1hist->Draw("histsame");
	if (datahist->Integral() > 0)
	{
		datahist->Draw("Esame");
		datahist->Draw("axissame");
	}
	if (prelim)
	{
		//WriteCEPC();
		WritePreliminary();
		WriteChannel(channel);
	}
	TArrow *arr1 = 0, *arr2 = 0;
	if (cut && strstr(scut, "CUT"))
	{
		Char_t cuts[10];
		Float_t x0 = 1, x1 = 1;
		//printf("%s\n",scut.Data());
		sscanf(scut.Data(), "%s %f %f", cuts, &x0, &x1);
		{
			Double_t min0 = datahist->GetBinContent(datahist->FindBin(x0));
			Double_t min1 = datahist->GetBinContent(datahist->FindBin(x1));
			Double_t ymin = max(min0, min1);
			Double_t ymax = datahist->GetMaximum();
			Double_t dely = datahist->GetMaximum() - ymin;
			ymin = ymin + dely * 0.05;
			ymax = (ymax - ymin) * 0.5 + ymin;
			//
			arr1 = new TArrow(x0, ymax, x0, ymin, 0.05, ">");
			arr1->SetLineWidth(5);
			arr1->SetLineColor(kMagenta);
			arr1->Draw("");
			//
			if (x0 != x1)
			{
				arr2 = new TArrow(x1, ymax, x1, ymin, 0.05, ">");
				arr2->SetLineWidth(5);
				arr2->SetLineColor(kMagenta);
				arr2->Draw("");
			}
		}
	}
	if (!logy)
	{
		if (datahist->Integral() > 0)
		{
			MakeLegend5(
				datahist, dataname,
				mc1hist, mc1name,
				mc2hist, mc2name,
				mc3hist, mc3name,
				mc4hist, mc4name);
		}
		else
		{
		}
		MakeLegend4(
			mc1hist, mc1name,
			mc2hist, mc2name,
			mc3hist, mc3name,
			mc4hist, mc4name);
	}
	char filenameall[256];
	sprintf(filenameall, "%s.eps", filename);
	//c1->SaveAs(filenameall);
	sprintf(filenameall, "%s.pdf", filename);
	c1->SaveAs(filenameall);
	if (arr1)
		delete arr1;
	if (arr2)
		delete arr2;
	delete c1;
}

// Plot data with one or more (fit) functions
// Functions should be part of the data histograms list of functions
// (i.e. perform fits with the "+" option or add other functions via
// datahist->GetListOfFunctions->Add(TF1 * function))
// functionnames should have at least as many elements as the function
// list
void PlotDataFit(char *filename,		// Name for the output files,
										// without extension
				 TH1D *datahist,		// Histogram with data
				 char *dataname,		// Description of data
				 char *functionnames[], // Names of associated functions
				 bool prelim,			// Preliminary plot
				 double xlow,			// Left edge of legend
										// (fraction of canavas width)
				 double ylow,			// Bottom edge of legend
										// (fraction of canavas height)
				 double xhi,			// Right edge of legend
										// (fraction of canavas width)
				 double yhi,			// Top edge of legend
										// (fraction of canavas height)
				 char *channel)
{
	SetStyle();
	if (prelim)
		SetPrelimStyle();
	else
		SetMeetingStyle();

	datahist->SetMaximum(datahist->GetMaximum() * 1.5);
	TCanvas *c1 = new TCanvas("CEPCplots", "CEPC Plots", 800, 800);

	FormatData(datahist);

	int linestyles[] = {
		1,
		1,
		1,
		10,
		9,
		7,
		6,
		5,
	};
	//int linecolors[]   = {2,4,kGreen+2,kOrange+7,kMagenta,2};
	int linecolors[] = {2, 4, 6, 9, 8, 2};

	TList *list = datahist->GetListOfFunctions();
	TH1D *datacopy = new TH1D(*datahist);
	datacopy->Draw("axis");

	unsigned int nfun = list->GetEntries();
	//cout<<"&&&&&&& "<<nfun<<" functions"<<endl;
	if (nfun > 6)
	{
		std::cout << "ERROR: More than six associated functions not forseen" << std::endl;
		return;
	}

	for (unsigned int i = 0; i < nfun; i++)
	{
		TF1 *f1 = (TF1 *)(list->At(i));
		//f1->Print();
		TString title = f1->GetTitle();
		//cout<<i<<" "<<title<<endl;
		if (title.Contains("special TPaveText"))
			continue;
		f1->SetLineColor(linecolors[i]);
		f1->SetLineStyle(linestyles[i]);
		f1->Draw("same");
	}

	MakeLegend(datahist, dataname, functionnames, xlow, ylow, xhi, yhi);

	datacopy->Draw("Esame");
	datacopy->Draw("axissame");

	if (prelim)
	{
		//WriteCEPC();
		WritePreliminary();
		WriteChannel(channel);
	}

	char filenameall[256];
	sprintf(filenameall, "%s.eps", filename);
	//c1->SaveAs(filenameall);
	sprintf(filenameall, "%s.pdf", filename);
	c1->SaveAs(filenameall);

	delete c1;
}

// Scatter plot
void PlotScatter(char *filename, // Name for the output files,
				 // without extension
				 TH1 *datahist, // Histogram with data
				 bool prelim	// preliminary plot
)
{

	SetStyle();
	if (prelim)
		SetPrelimStyle();
	else
		SetMeetingStyle();

	TCanvas *c1 = new TCanvas("CEPCplots", "CEPC Plots", 800, 800);

	FormatData(datahist);

	if (datahist->Integral() > 5000)
		datahist->SetMarkerStyle(1);
	else if (datahist->Integral() > 500)
		datahist->SetMarkerSize(makersize);

	datahist->Draw("BOX");

	if (prelim)
	{
		//WriteCEPC();
		WritePreliminary();
	}

	char filenameall[256];
	sprintf(filenameall, "%s.eps", filename);
	//c1->SaveAs(filenameall);
	sprintf(filenameall, "%s.pdf", filename);
	c1->SaveAs(filenameall);
	delete c1;
}

// Scatter plot
void PlotScatterCircle(char *filename, // Name for the output files,
					   // without extension
					   TH1 *datahist, // Histogram with data
					   double cx,
					   double cy,
					   double cr,
					   bool prelim // preliminary plot
)
{
	SetStyle();
	if (prelim)
		SetPrelimStyle();
	else
		SetMeetingStyle();

	TCanvas *c1 = new TCanvas("CEPCplots", "CEPC Plots", 800, 800);

	FormatData(datahist);

	if (datahist->Integral() > 5000)
		datahist->SetMarkerStyle(1);
	else if (datahist->Integral() > 500)
		datahist->SetMarkerSize(makersize);

	datahist->Draw("BOX");

	if (prelim)
	{
		//WriteCEPC();
		WritePreliminary();
	}

	TArc arc(cx, cy, cr);
	arc.SetLineWidth(4);
	arc.SetLineColor(kRed);
	arc.SetFillColor(0);
	arc.SetFillStyle(0);
	arc.Draw("SAME");

	char filenameall[256];
	sprintf(filenameall, "%s.eps", filename);
	//c1->SaveAs(filenameall);
	sprintf(filenameall, "%s.pdf", filename);
	c1->SaveAs(filenameall);
	delete c1;
}

void PlotDataRooFit(Char_t *figfilename,
					RooAbsData *data, RooRealVar x,
					std::vector<RooAbsPdf *> pdfList,
					char *xtitle, char *ytitle, bool prelim, char *channel, double maxy)
{

	SetStyle();
	if (prelim)
		SetPrelimStyle();
	else
		SetMeetingStyle();
	//
	int linestyles[] = {
		1,
		1,
		1,
		10,
		9,
		7,
		6,
		5,
	};
	int kColor[12] = {
		kBlack, kBlue, kRed, kGreen,
		kViolet,
		kSpring,
		kCyan,
		kOrange,
		kTeal,
		kMagenta,
		kAzure,
		kPink};
	Double_t num = data->sumEntries() / 100;
	Double_t yoffset = 1.0 + 0.1 * log(num) / log(10.);
	//
	RooPlot *xframe = x.frame(100);
	FormatRooFit(xframe, xtitle, ytitle, yoffset);
	if (maxy > 0)
		xframe->SetMaximum(maxy);
	//
	data->plotOn(xframe, MarkerStyle(21));
	for (unsigned int i = 0; i < pdfList.size(); i++)
	{
		if (i == 0)
			pdfList[0]->plotOn(xframe, Components(*pdfList[i]), LineStyle(linestyles[i]), LineWidth(3), LineColor(kColor[i + 1]));
		else
			pdfList[0]->plotOn(xframe, Components(*pdfList[i]), LineStyle(linestyles[i]), LineWidth(3), LineColor(kColor[i + 1]));
	}
	if (!prelim)
	{
		pdfList[0]->paramOn(xframe, Format("NE", AutoPrecision(2)), Layout(0.60, 0.95, 0.95));
		//pdfList[0]->paramOn(xframe,Format("NE",FixedPrecision(2)), Layout(0.60,0.95,0.93));
	}
	//
	TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
	xframe->Draw();
	//
	if (prelim)
	{
		//WriteCEPC();
		WritePreliminary();
		WriteChannel(channel);
	}
	//
	TLegend *leg = new TLegend(0.6, 0.65, 0.9, 0.80);
	Char_t name[100];
	if (data)
		leg->AddEntry(data, "CEPC Simulation", "LEPF");
	for (unsigned int i = 0; i < pdfList.size(); i++)
	{
		TLine *line = new TLine(130, 760, 132, 760);
		line->SetLineWidth(3);
		line->SetLineColor(kColor[i + 1]);
		line->SetLineStyle(linestyles[i]);
		sprintf(name, "%s", pdfList[i]->GetTitle());
		if (i == 0)
			leg->AddEntry(line, name, "LF");
		else
			leg->AddEntry(line, name, "LF");
	}

	leg->SetFillColor(0);
	leg->SetTextFont(22);
	leg->SetTextSize(0.03);
	leg->Draw();

	sprintf(name, "%s.eps", figfilename);
	//c1->Print(name);
	sprintf(name, "%s.pdf", figfilename);
	c1->Print(name);
	sprintf(name, "%s.png", figfilename);
	c1->Print(name);
}
