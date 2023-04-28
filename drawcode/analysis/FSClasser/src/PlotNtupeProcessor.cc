#include "PlotNtupleProcessor.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <utility>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <algorithm>   

// ----- include for verbosity dependend logging ---------
#include <marlin/VerbosityLevels.h>

#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <cepcplotstyle.h>

using namespace std ;
using namespace lcio ;
using namespace marlin ;

PlotNtupleProcessor aPlotNtupleProcessor ;

PlotNtupleProcessor::PlotNtupleProcessor()
	: Processor("PlotNtupleProcessor") 
{

	// Processor description
	_description = "Plot ntuple variables Processor" ;

	registerProcessorParameter( "RootFileNames" , 
			"input root files"  ,
			_rootNames ,
			std::vector<std::string>() ) ;

	registerProcessorParameter( "LegendNames" , 
			"Legend name of samples"  ,
			_legendNames ,
			std::vector<std::string>() ) ;

	registerProcessorParameter( "TreeNames" , 
			"tree names"  ,
			_treeNames ,
			std::vector<std::string>() ) ;

	registerProcessorParameter( "VarsNames" , 
			"variable names"  ,
			_variableNames ,
			std::vector<std::string>() ) ;
	
	registerProcessorParameter("Cut",    "A cut string for the ntuple",    _cut, std::string(""));
	registerProcessorParameter("FigDir", "the dir for the output figs",    _dir, std::string("figs"));
	registerProcessorParameter("OnlyMC", "if true, only MC plotted"   , _onlyMC, 0);
}


void PlotNtupleProcessor::init() { 
	streamlog_out(DEBUG) << "   init called  " << std::endl ;
	printParameters() ;
        if( _treeNames.size()     == 0 ) exit(1);
        if( _rootNames.size()     == 0 ) exit(1);
        if( _variableNames.size() == 0 ) exit(1);
	vector<TFile*> _rootFiles;
	for (unsigned int i=0; i<_rootNames.size(); ++i) {
		TFile* f = new TFile(TString::Format("%s",_rootNames[i].data()));
		_rootFiles.push_back(f);
	}
	//
	for (unsigned int j=0; j<_treeNames.size(); ++j) {
		TTree* ntp = NULL ; 
		vector< TTree* > _trees;
		//
		for (unsigned int i=0; i<_rootNames.size(); ++i) {
			ntp = (TTree*) _rootFiles[i]->Get( _treeNames[j].data() );
			_trees.push_back(ntp);
		}
		//cout<<"Number of trees is "<<_trees.size()<<endl;
		//cout<<_cut<<endl;
		//
		char hname[100], filename[100], title[100];
		vector<char*> _histoname;
		for (unsigned int l=0; l<_legendNames.size(); ++l) {
			_histoname.push_back((char*)_legendNames[l].c_str() );	
		}
		for (unsigned int k=0; k<_variableNames.size(); k+=10) {
			int     nbin  = atoi(_variableNames[k+2].c_str());
			double  xmin  = atof(_variableNames[k+3].c_str());
			double  xmax  = atof(_variableNames[k+4].c_str());
			bool    pril  = atoi(_variableNames[k+5].c_str());
			bool    logy  = atoi(_variableNames[k+6].c_str());
			bool    cut   = atoi(_variableNames[k+7].c_str());
			vector< TH1D* > _histo;
			//
			for (unsigned int l=0; l<_trees.size(); ++l) {
				sprintf(hname,"%s_%d",_variableNames[k+1].data(),l);
				TH1D *hd = new TH1D( hname,  _variableNames[k+1].data(), nbin,xmin,xmax);	
			        NameAxes( hd, _variableNames[k+8].c_str(), _variableNames[k+9].c_str() );
				_histo.push_back(hd);	
			}
			//
			for (unsigned int l=0; l<_trees.size(); ++l) {
				sprintf(hname,"%s_%d",_variableNames[k+1].data(),l);
				//cout<<"hname  "<< hname<<"  variable "<< _variableNames[k].data()<<" cut " <<_cut.data()<<endl;
				if(l==0&&_onlyMC>0)continue; 
				_trees[l]->Project(hname, _variableNames[k].data(), _cut.data());
			}
			//
			sprintf(filename,"%s/%s", _dir.data(), _variableNames[k].data());
			sprintf(title,"%s", _variableNames[k+1].c_str());
			//
			PlotDataMC(     filename, 
					_histo,  _histoname, 
					pril, logy, cut, title
				  );
			//
			for (unsigned int l=0; l<_histo.size(); ++l) {
				delete _histo[l];
			}
		}

	}
	for (unsigned int i=0; i<_rootNames.size(); ++i) {
		delete _rootFiles[i];
	}
}

void PlotNtupleProcessor::processRunHeader( LCRunHeader* run) { 
} 

void PlotNtupleProcessor::processEvent( LCEvent * evt ) { 
}

void PlotNtupleProcessor::check( LCEvent * evt ) { 
}

void PlotNtupleProcessor::end() { 
}
