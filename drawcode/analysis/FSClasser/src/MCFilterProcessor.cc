#include "MCFilterProcessor.h"
#include <iostream>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include "IMPL/ReconstructedParticleImpl.h"
#include "UTIL/LCRelationNavigator.h"
#include "UTIL/PIDHandler.h"
#include "EVENT/LCFloatVec.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include <UTIL/LCTOOLS.h>
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/SISConePlugin.hh"

#include <iostream>
#include <sstream>
#include <fstream>
#include <utility>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <algorithm>   

using namespace lcio ;
using namespace marlin ;
using namespace std ;

MCFilterProcessor aMCFilterProcessor ;

MCFilterProcessor::MCFilterProcessor() : Processor("MCFilterProcessor") {

	// Processor description
	_description = "FastJet Clustering ..." ;


	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			"InputCollection",
			"Collection of reconstructed particles",
			_inputCollection,
			std::string("Unset") );

	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			"InputMCTruthMap",
			"Collection of reconstructed particles",
			_inputMCTruthMap,
			std::string("Unset") );

	registerProcessorParameter("Debug",
			"debug printout",
			_print,
			int(0)); 

	registerProcessorParameter("Algorithm",
			"FastJet algorithm",
			sAlgorithm,
			std::string("antikt_algorithm")); 

	registerProcessorParameter("R",
			"R Parameter",
			_RPar,
			double(0.7)); 
	
	registerProcessorParameter("P",
			"P Parameter",
			_PPar,
			double(0.0)); 

	registerProcessorParameter("EjetMin",
			"Ejet",
			_eJet,
			double(10.0)); 

	registerProcessorParameter("PtMin",
			"PtMin Cut",
			_PtCut,
			double(0.0)); 

	registerProcessorParameter("NJets",
			"max nb of jets",
			_nJetMax,
			int(25)); 

	registerProcessorParameter("nPFOmin",
			"min nb of input PFOs",
			_nPFOmin,
			int(5)); 

	registerProcessorParameter("InclusiveExclusive",
			"Inclusive = 0; Exclusive != 0",
			_InclusiveExclusive,
			int(1)); 

	registerProcessorParameter("EnergyCut",
			"Energy Threshold of particles",
			_EnergyCut,
			double(0)); 
	
	registerProcessorParameter("OutputLcioFile",
			" output file for lcio collections ",
			_lcioOutputFile,
			std::string("output.slcio")); 

	registerProcessorParameter("OutputRejectedLcioFile",
			" output file for lcio collections ",
			_lcioOutputFile_rej,
			std::string("reject.slcio")); 


}

void MCFilterProcessor::init() { 

	printParameters() ;

	if(sAlgorithm=="kt_algorithm"){
		SetAlgorithm(fastjet::kt_algorithm);
	}else if(sAlgorithm=="cambridge_algorithm"){ 
		SetAlgorithm(fastjet::cambridge_algorithm);
	}else if(sAlgorithm=="antikt_algorithm"){ 
		SetAlgorithm(fastjet::antikt_algorithm);
	}else if(sAlgorithm=="ee_kt_algorithm"){ 
		SetAlgorithm(fastjet::ee_kt_algorithm);
	}else if(sAlgorithm=="genkt_algorithm"){ 
		SetAlgorithm(fastjet::genkt_algorithm);
	}else if(sAlgorithm=="ee_genkt_algorithm"){ 
		SetAlgorithm(fastjet::ee_genkt_algorithm);
	}else if(sAlgorithm=="plugin_algorithm"){ 
		SetAlgorithm(fastjet::plugin_algorithm);
	}else{
       printf("please chooose the proper algorithm, now it is %s\n", sAlgorithm.c_str());
		 exit(1);
	}

    _lcWrt = LCFactory::getInstance()->createLCWriter() ;
    _lcWrt->open( _lcioOutputFile ) ;
    _lcWrt_rej = LCFactory::getInstance()->createLCWriter() ;
    _lcWrt_rej->open( _lcioOutputFile_rej ) ;



}

void MCFilterProcessor::processRunHeader( LCRunHeader* run) { 

	_nRun = run->getRunNumber() ;
	_lcWrt->writeRunHeader( run ) ;
	_lcWrt_rej->writeRunHeader( run ) ;

} 

void MCFilterProcessor::processEvent( LCEvent * evt ) { 

	_nRun = evt->getRunNumber();
	_nEvt = evt->getEventNumber();

	if(_print>0)cout <<"Run " << _nRun << " Evt " << _nEvt << endl;

	for(int i1=0;i1<5;i1++){
		for(int i2=0;i2<50;i2++){
			_jetVector[i1][i2]=0.;
		}
	}




	LCCollection* enflowcol  = evt->getCollection(_inputCollection);
	//LCCollection* McPartCol  = evt->getCollection("MCParticle");
	//LCCollection* MCTruthMap = evt->getCollection(_inputMCTruthMap);
	
	int nenflow =  enflowcol->getNumberOfElements(); 
	if(_print>1)cout <<" # of tracks & clusters " << nenflow << endl;


	fastjet::JetAlgorithm algorithm = GetAlgorithm();
	if(_print>1)cout <<" Jet algorithm is  " << algorithm << endl;
	
	vector<fastjet::PseudoJet> input_particles;
	//
	//
	//
	double px, py, pz, E;
	for ( int ienflow=0; ienflow<nenflow ; ienflow++){
		ReconstructedParticle* enflow = dynamic_cast<ReconstructedParticle*>(enflowcol->getElementAt( ienflow ));
		if ( enflow->getEnergy()<_EnergyCut ) continue;;
		if( 1 ){
			px = enflow->getMomentum()[0];
			py = enflow->getMomentum()[1];
			pz = enflow->getMomentum()[2];
			E  = enflow->getEnergy();

			fastjet::PseudoJet thisPtc(px,py,pz,E);
			thisPtc.set_user_index(ienflow);
			if(_print>2) printf("id, px, py, pz, E= %5d, %10.4f, %10.4f, %10.4f, %10.4f\n", enflow->id(),px, py, pz, E);

			input_particles.push_back(thisPtc);
		}
	}

	if(_print>1) printf("%3d PFOs used for jet-clustering\n", (int)input_particles.size());
	double _ymin[8];
	memset(_ymin, 0, sizeof(double)*8 );
	
	vector<ReconstructedParticleImpl*> trash;
	if( (int)input_particles.size() > _nPFOmin ) {
		_rp = _RPar;
		_pp = _PPar;

		fastjet::Strategy strategy = fastjet::Best;

		fastjet::RecombinationScheme recomb_scheme = fastjet::E_scheme;

		float momentum[3], energy, mass;

		fastjet::JetDefinition* jet_def =NULL;

		if(algorithm==fastjet::kt_algorithm || algorithm==fastjet::antikt_algorithm || algorithm==fastjet::cambridge_algorithm){
			jet_def = new fastjet::JetDefinition(algorithm, _rp, recomb_scheme, strategy); 
		}else if(algorithm==fastjet::cambridge_algorithm){
			jet_def = new fastjet::JetDefinition(algorithm, _rp, recomb_scheme, strategy);
		}else if(algorithm==fastjet::ee_kt_algorithm){ 
			jet_def = new fastjet::JetDefinition(algorithm); 
		}else if(algorithm==fastjet::genkt_algorithm || algorithm==fastjet::ee_genkt_algorithm ){ 
			jet_def = new fastjet::JetDefinition(algorithm, _rp, _pp, recomb_scheme, strategy); 
		}else if (algorithm==fastjet::plugin_algorithm ){
			double cone_radius = _rp;
			double overlap_threshold = _pp;
			fastjet::SISConePlugin siscone(cone_radius, overlap_threshold);
			jet_def = new fastjet::JetDefinition(& siscone);
		}

		fastjet::ClusterSequence cs(input_particles, *jet_def);

		vector<fastjet::PseudoJet> JetsVec ;

		if( _InclusiveExclusive ) JetsVec = cs.exclusive_jets(_nJetMax);
		else                      JetsVec = cs.inclusive_jets(_PtCut  );
		
		vector<fastjet::PseudoJet> sortedJets = sorted_by_E(JetsVec);

	
		for(int iy=1; iy<9;iy++){	
			_ymin[iy-1] = cs.exclusive_ymerge (iy); // the ymin corresponding to the recombination that went from iy+1 to iy jets
			if(_print>1)printf("%3d,  ymin = %12.5f\n", iy, _ymin[iy-1] );
		}	
		

		int nmx = sortedJets.size();
		if( _InclusiveExclusive>0 && nmx > _nJetMax) nmx = _nJetMax;

		_nJets = sortedJets.size();

		delete jet_def; 


		if(_print>1)printf("MCFilterProcessor: Nb of Jets %3d\n", (int)sortedJets.size() );
		_nJetsHE=0;
		for(int ij=0; ij<_nJets;ij++){

			momentum[0] = sortedJets[ij].px();
			momentum[1] = sortedJets[ij].py();
			momentum[2] = sortedJets[ij].pz();
			energy      = sortedJets[ij].e();
			mass        = sortedJets[ij].m();

			if(_print>1) printf( "Jet %2d  Energy = %10.4f\n", ij,  energy );

			/*
			if(energy>_eJet && ij < 10 ){
				
				vector<fastjet::PseudoJet> jetConstituents = cs.constituents(sortedJets[ij]);
				//if( jetConstituents.size() < _nJetMax ) continue;

				_jetVector[0][ij] = sortedJets[ij].px();
				_jetVector[1][ij] = sortedJets[ij].py();
				_jetVector[2][ij] = sortedJets[ij].pz();
				_jetVector[3][ij] = sortedJets[ij].e();
				_jetVector[4][ij] = (float)jetConstituents.size();

				_nJetsHE++;

				if( _InclusiveExclusive==0 || ij<_nJetMax){
					ReconstructedParticleImpl* Jets = new ReconstructedParticleImpl;
					for(unsigned ip = 0; ip < jetConstituents.size(); ip++){
						if(jetConstituents[ip].user_index()>=0 && jetConstituents[ip].user_index() < enflowcol->getNumberOfElements()){
							ReconstructedParticle* enflow = dynamic_cast<ReconstructedParticle*>(enflowcol->getElementAt(jetConstituents[ip].user_index()));
							Jets->addParticle(enflow);
							if(_print>2) cout << "add to particle list  " << enflow->id() << endl;
						}
					}
					Jets->setMomentum(momentum);
					Jets->setEnergy(energy);
					Jets->setMass(mass);
					Jets->setType(4);
					trash.push_back(Jets);	
				}
			}
			*/
		}

	}

	for( int iy=1; iy<9; iy++){
		char yname[10];
		sprintf(yname,"y%d%d", iy,iy+1);
		//printf("%s = %10.5f\n", yname, _ymin[iy-1] );
	}

	if(_ymin[3]<0.0001 )
	  _lcWrt->writeEvent( evt ) ;
	else
	  _lcWrt_rej->writeEvent( evt ) ;




	//FreeDelAll(trash);	
}

void  MCFilterProcessor::check( LCEvent * evt ) { 
}


void MCFilterProcessor::end(){ 
	std::cout << "MCFilterProcessor::end()  " 
		<< " processed " << _nEvt << " events in " << _nRun << " runs "
		<< std::endl ;
	  _lcWrt->close() ;
	  _lcWrt_rej->close() ;

}
