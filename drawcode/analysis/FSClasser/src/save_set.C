#include <save_set.h>
#include <EVENT/LCCollection.h>
#include <EVENT/ReconstructedParticle.h>
#include <values.h>
#include <string>
#include <iostream>
#include <EVENT/LCFloatVec.h>
#include <EVENT/LCParameters.h>
#include <stdexcept>
#include <TFile.h> 
#include <TTree.h>
#include <Rtypes.h> 
#include <sstream>		
#include <cmath>
#include <math.h>
#include <TVector3.h>
#include <cepcplotstyle.h>

save_set a_save_set_instance;

save_set::save_set()
	: Processor("save_set"),
	_output(0)
{
	_description = "Print MC Truth" ;



	_treeFileName="WWV.root";
	registerProcessorParameter( "TreeOutputFile" , 
			"The name of the file to which the ROOT tree will be written" ,
			_treeFileName ,
			_treeFileName);

	_treeName="WWV";
	registerProcessorParameter( "TreeName" , 
			"The name of the ROOT tree" ,
			_treeName ,
			_treeName);


	_overwrite=0;
	registerProcessorParameter( "OverwriteFile" , 
			"If zero an already existing file will not be overwritten." ,
			_overwrite ,
			_overwrite);

}

void save_set::init() {

	printParameters();

	TFile *tree_file=new TFile(_treeFileName.c_str(),(_overwrite ? "RECREATE" : "UPDATE"));


	if (!tree_file->IsOpen()) {
		delete tree_file;
		tree_file=new TFile(_treeFileName.c_str(),"NEW");
	}

	_outputMCP = new TTree(_treeName.c_str(), "WWV");
	_outputMCP->SetAutoSave(32*1024*1024);  // autosave every 32MB
	_outputMCP->Branch("Run", &_Run, "Run/D");
	_outputMCP->Branch("Evt", &_Evt, "Evt/D");
	_outputMCP->Branch("CosTheta", _CosTheta, "CosTheta[500]/D");
	_outputMCP->Branch("Phi"     ,      _Phi,      "Phi[500]/D");
	_outputMCP->Branch("Mass"    ,     _Mass,     "Mass[500]/D");
	_outputMCP->Branch("posx"    ,     _posx,     "posx[500]/D");
	_outputMCP->Branch("posy"    ,     _posy,     "posy[500]/D");
	_outputMCP->Branch("posz"    ,     _posz,     "posz[500]/D");
	_outputMCP->Branch("dist"    ,     _dist,     "dist[500]/D");
	_outputMCP->Branch("PDGi"    ,     _PDGi,     "PDGi[500]/D");
}

void save_set::processEvent( LCEvent * evtP ) 
{		

	if (evtP) 								
	{
		try{
			//LCCollection* col_SET = evtP->getCollection( "SETCollection" ) ;
			LCCollection* col_SET = evtP->getCollection( "SETSpacePoints" ) ;

			int _nSET=col_SET->getNumberOfElements();
			_Run=evtP->getEventNumber();	
			_Evt=evtP->getEventNumber();	
			for(int j = 0; j < _nSET; j++)
			{
				if (j>500) break; 
				TrackerHit* a_Hit = dynamic_cast<TrackerHit*>(col_SET->getElementAt(j));
				//_PDGi    [j]= (a_Hit->getMCParticle())->getPDG();
				//_Mass    [j]= a_Hit->getMCParticle()->getMass();
				_posx    [j]= a_Hit->getPosition()[0];
				_posy    [j]= a_Hit->getPosition()[1];
				_posz    [j]= a_Hit->getPosition()[2];
				_posr    [j]= pow(_posx[j]*_posx[j]+_posy[j]*_posy[j]+_posz[j]*_posz[j], 0.5 );
				_CosTheta[j]=_posz[j]/_posr[j];
				_Phi     [j]= atan2(_posx[j], _posy[j]);  
			}

			for(int i = 0; i < _nSET; i++)
			{
				if (i>500) break; 
				double d=99999; 
				for(int j = i+1; j < _nSET; j++)
				{
					if (j>500) break; 
					double dd = pow(
							(_posx[j] - _posx[i])*(_posx[j] - _posx[i]) + 
							(_posy[j] - _posy[i])*(_posy[j] - _posy[i]) + 
							(_posz[j] - _posz[i])*(_posz[j] - _posz[i]), 0.5 );
					if ( dd<d ) d =dd;  
				}
				_dist[i] = d; 
			}
			_outputMCP->Fill();
		}catch (lcio::DataNotAvailableException err) { }
	}
}	
//
//******************************************************************
//
TLorentzVector save_set::HPINV(const TLorentzVector &P1, const TLorentzVector &PT){
	const Double_t EPS=1.0e-6;
	//
	Double_t AM0 = PT.M();
	if ( PT.M() < EPS ) {
		printf("AM0=%10.8f\n",AM0);
		printf("Transformation Velocity Is Equal the Speed Of Light\n");
	} 
	Double_t beta  =  PT.Beta(); 
	Double_t the   =  PT.Theta();
	Double_t phi   =  PT.Phi();
	//
	if( beta>EPS ) {	
		//Double_t pp    =  PT.Rho(); 
		Double_t gam   =  PT.Gamma();
		Double_t betx  =  PT.BoostVector().X();
		Double_t bety  =  PT.BoostVector().Y();
		Double_t betz  =  PT.BoostVector().Z();

		Double_t X =  cos(the)*cos(phi)*P1.X()+cos(the)*sin(phi)*P1.Y()-     sin(the)*P1.Z();
		Double_t Y = -         sin(phi)*P1.X()+         cos(phi)*P1.Y()                     ;
		Double_t Z =  gam*betx/beta    *P1.X()+gam*bety/beta    *P1.Y()+gam*betz/beta*P1.Z()-gam*beta*P1.T();
		Double_t T = -gam*betx         *P1.X()-gam*bety         *P1.Y()-gam*betz     *P1.Z()+gam     *P1.T();
		return TLorentzVector(X,Y,Z,T);
	}else{
		return erout4( P1, phi, the );
	}
}
//
//******************************************************************
//
void save_set::GetHelicityAngles( const vector<TLorentzVector>& p4List, Double_t *theta, Double_t *phi )
{
	vector<TLorentzVector> List; List.clear();
	if ( 	p4List.size() < 2 ) {
		cout<<"less than 2 input vectors"<<endl;
		exit(1);	 
	}else{ 
		TLorentzVector mother = p4List[0];
		//printf("mass = %10.6f\n", mother.M());
		if ( p4List.size() >2 ){ 
			for( unsigned int i =1; i< p4List.size(); i++ ){
				List.push_back(HPINV(p4List[i],mother));
			}
			GetHelicityAngles( List, theta, phi );
		}else{
			TLorentzVector p4 = HPINV(p4List[1],mother);
			*theta =  p4.Theta();
			*phi   =  p4.Phi();
		}
	}
	List.clear(); 
}

void save_set::end()
{
	MakePlots();
	if (_outputMCP) {
		TFile *tree_file = _outputMCP->GetCurrentFile(); //just in case we switched to a new file
		tree_file->Write();
		//delete tree_file;
	}
}

int save_set::FindLeptonParent(MCParticle *mcp)
{
	int ret = 0;

	int PDGIDItself = abs(mcp->getPDG()); 
	if ( mcp->getParents().size() == 0 )      // original particle, no parent  
	{
		if ( PDGIDItself==11  ||  PDGIDItself==13 ){
			ret= PDGIDItself;
		}
	} 
	else if ( mcp->getParents().size() == 1 ) // Only one parent  
	{
		MCParticle *mcpa = mcp->getParents()[0];
		int PDGIDParent = abs(mcpa->getPDG()); 
		if ( PDGIDParent >100){
			PDGIDParent %= 1000;
		}
		if ( ( PDGIDItself ==11 || PDGIDItself ==13 ) && (PDGIDParent == PDGIDItself) ){
			ret = FindLeptonParent(mcpa);	
		}else	
			if (     ( PDGIDItself ==11 || PDGIDItself ==13 ) 
					&& ( PDGIDParent ==13 || PDGIDParent ==15 || PDGIDParent ==23 || PDGIDParent ==24 || (PDGIDParent>100&&PDGIDParent<1000) )  
				){
				// parent is W/Z/tau/B/D and itself is a lepton 
				ret = PDGIDParent*100 + PDGIDItself;
			}
	}
	else                // More than one parent  
	{
		printf("There are more than one parent: %2d \n", PDGIDItself);	
	}
	//
	return ret;
}
//
//******************************************************************
//
TLorentzVector save_set::erout4(const TLorentzVector &p1, const double phi, const double theta )
{
	double cp = cos(phi);
	double sp = sin(phi);
	double ct = cos(theta);
	double st = sin(theta);
	double t  = p1.E(), x=p1.X(), y=p1.Y(), z=p1.Z();
	double xp = x*cp*ct+y*sp*ct-z*st, yp= -x*sp+y*cp, zp=x*cp*st+y*sp*st+z*ct;
	return TLorentzVector(xp,yp,zp,t);
}


void save_set::MakePlots(){

	for(int i=0; i<6; i++) {
		if(h_CosTheta[i] && h_CosTheta[i]->Integral()>0){
			TH1D *h2 = 0; 
			char filename[256], title[256];
			sprintf(filename,"figs/%s", (h_CosTheta[i]->GetName()));
			sprintf(title,"%s", (h_CosTheta[i]->GetTitle()));
			//
			h_CosTheta[i]->SetMaximum( h_CosTheta[i]->GetMaximum()*1.2); 
			PlotDataMC(filename, 
					h_CosTheta[i], (char*)"",
					h_CosTheta[i], (char*)"",
					h2       , (char*)"",
					h2       , (char*)"",
					h2       , (char*)"",
					false, false, false, title
					);
		}
		if(h_Phi[i] && h_Phi[i]->Integral()>0){
			TH1D *h2 = 0; 
			char filename[256], title[256];
			sprintf(filename,"figs/%s", (h_Phi[i]->GetName()));
			sprintf(title,"%s", (h_Phi[i]->GetTitle()));
			//
			h_Phi[i]->SetMaximum( h_Phi[i]->GetMaximum()*1.2); 
			PlotDataMC(filename, 
					h_Phi[i] , (char*)"",
					h_Phi[i] , (char*)"",
					h2       , (char*)"",
					h2       , (char*)"",
					h2       , (char*)"",
					false, false, false, title
					);
		}
	}
}
