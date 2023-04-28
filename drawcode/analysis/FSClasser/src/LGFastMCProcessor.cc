#include "marlin/MarlinConfig.h" // defines MARLIN_CLHEP 

#include "marlin/FastMCParticleType.h"
#include "marlin/ErrorOfSigma.h"


//--- LCIO headers 


#include <iostream>
#include <cmath>

#include "TVector3.h" 
#include "TLorentzVector.h" 
#include "LGFastMCProcessor.h"
#include "LGParticleFactory.h"
#include "LGTrackSmearer.h"
#include "LGClusterSmearer.h"
#include "cepcplotstyle.h"

using namespace lcio ;


namespace marlin{


	LGFastMCProcessor aLGFastMCProcessor ;


	LGFastMCProcessor::LGFastMCProcessor() : Processor("LGFastMCProcessor"),
	_factory(NULL),
	_nRun(-1),
	_nEvt(-1)
	{

		// modify processor description
		_description = "LGFastMCProcessor creates ReconstrcutedParticles from MCParticles " 
			"according to the resolution given in the steering file." ;


		// register steering parameters: name, description, class-variable, default value

		registerInputCollection( LCIO::MCPARTICLE,
				"InputCollectionName" , 
				"Name of the MCParticle input collection"  ,
				_inputCollectionName ,
				std::string("MCParticle") ) ;


		registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
				"RecoParticleCollectionName" , 
				"Name of the ReconstructedParticles output collection"  ,
				_recoParticleCollectionName ,
				std::string("ReconstructedParticles") ) ;

		registerOutputCollection( LCIO::LCRELATION,
				"MCTruthMappingCollectionName" , 
				"Name of the MCTruthMapping output collection"  ,
				_mcTruthCollectionName ,
				std::string("MCTruthMapping") ) ;


		registerProcessorParameter( "MomentumCut" , 
				"No reconstructed particles are produced for smaller momenta (in [GeV])"  ,
				_momentumCut ,
				float( 0.001 ) ) ;

		FloatVec chResDefault ;
		chResDefault.push_back( 5e-5      ) ;
		chResDefault.push_back( 1e-3      ) ;
		chResDefault.push_back( 3.50      ) ;
		chResDefault.push_back( -1.0      ) ;
		chResDefault.push_back(  1.0      ) ;

		registerProcessorParameter( "ChargedResolution" , 
				"Resolution of charged particles in polar angle range:  d(1/P)  th_min  th_max"  ,
				_initChargedRes ,
				chResDefault ,
				chResDefault.size() ) ;

		FloatVec gammaResDefault ;
		gammaResDefault.push_back( 0.01      ) ;
		gammaResDefault.push_back( 0.10      ) ;
		gammaResDefault.push_back( -1.0      ) ;
		gammaResDefault.push_back(  1.0      ) ;

		registerProcessorParameter( "PhotonResolution" , 
				"Resolution dE/E=A+B/sqrt(E/GeV) of photons in polar angle range: A  B th_min  th_max"  ,
				_initPhotonRes ,
				gammaResDefault ,
				gammaResDefault.size() ) ;

		FloatVec hadronResDefault ;
		hadronResDefault.push_back( 0.04      ) ;
		hadronResDefault.push_back( 0.50      ) ;
		hadronResDefault.push_back( -1.0      ) ;
		hadronResDefault.push_back(  1.0      ) ;

		registerProcessorParameter( "NeutralHadronResolution" , 
				"Resolution dE/E=A+B/sqrt(E/GeV) of neutral hadrons in polar angle range: A  B th_min  th_max"  ,
				_initNeutralHadronRes ,
				hadronResDefault ,
				hadronResDefault.size() ) ;

		registerProcessorParameter("MakePlots",     "Make some plots for check",  m_makeplots,       0);
		registerProcessorParameter("Smear",         "modeling the detector res.", m_smear,           1);
		registerProcessorParameter("Luxury",        "save 4 momentum to ntuple ", m_luxury,          0);
		registerProcessorParameter("RejectNeutrino","reject the undetectables  ", m_rejectNeutrino,  1);


	}


	void LGFastMCProcessor::init() { 

		// usually a good idea to
		printParameters() ;

		_nRun = 0 ;
		_nEvt = 0 ;



		_factory = 0 ;
#ifdef MARLIN_CLHEP

		LGParticleFactory* simpleFactory  =  new LGParticleFactory() ; 

		simpleFactory->registerIFourVectorSmearer(  new LGTrackSmearer  ( _initChargedRes       ), CHARGED        ) ;
		simpleFactory->registerIFourVectorSmearer(  new LGClusterSmearer( _initPhotonRes        ), PHOTON         ) ;
		simpleFactory->registerIFourVectorSmearer(  new LGClusterSmearer( _initNeutralHadronRes ), NEUTRAL_HADRON ) ;
		simpleFactory->setMomentumCut( _momentumCut ) ;
		simpleFactory->setSmear( m_smear ) ;
		simpleFactory->setNeutrino( m_rejectNeutrino ) ;

		_factory = simpleFactory ;

		streamlog_out( MESSAGE )  << " LGFastMCProcessor::init() : registering LGParticleFactory " << std::endl ;

#endif // MARLIN_CLHEP
		if( m_luxury>0 ) { 
			char* ntFullName = new char[100];
			sprintf(ntFullName,"nt%s","MCTOPO");
			TTree* outputTree = new TTree(ntFullName,"");
			outputTree->SetAutoSave(32*1024*1024); 
			m_ntp = new NTupleHelper( outputTree  ); 
			delete ntFullName;
		}
		if( m_makeplots>0 ) { 
			SetPrelimStyle();	
			h_Momentum[ 0] = new TH1D("E_gamma"      , "Energy   of #gamma  ", 150,  0.0, 150.0);
			h_Momentum[ 1] = new TH1D("p_electron"   , "Momentum of e^{+}   ", 150,  0.0, 150.0);
			h_Momentum[ 2] = new TH1D("p_positron"   , "Momentum of e^{-}   ", 150,  0.0, 150.0);
			h_Momentum[ 3] = new TH1D("p_muonplus"   , "Momentum of #mu^{+} ", 150,  0.0, 150.0);
			h_Momentum[ 4] = new TH1D("p_muonminus"  , "Momentum of #mu^{-} ", 150,  0.0, 150.0);
			h_Momentum[ 5] = new TH1D("p_pionplus"   , "Momentum of #pi^{+} ", 50,  0.0, 50.0);
			h_Momentum[ 6] = new TH1D("p_pionminus"  , "Momentum of #pi^{-} ", 50,  0.0, 50.0);
			h_Momentum[ 7] = new TH1D("p_kaonplus"   , "Momentum of K^{+}   ", 50,  0.0, 50.0);
			h_Momentum[ 8] = new TH1D("p_kaonminus"  , "Momentum of K^{-}   ", 50,  0.0, 50.0);
			h_Momentum[ 9] = new TH1D("p_proton"     , "Momentum of p       ", 50,  0.0, 50.0);
			h_Momentum[10] = new TH1D("p_antiproton" , "Momentum of #bar{p} ", 50,  0.0, 50.0);
			h_Momentum[11] = new TH1D("p_neutron"    , "Momentum of n       ", 50,  0.0, 50.0);
			h_Momentum[12] = new TH1D("p_antineutron", "Momentum of #bar{n} ", 50,  0.0, 50.0);
			h_Momentum[13] = new TH1D("p_Klong"      , "Momentum of K_{L}   ", 50,  0.0, 50.0);
			h_Momentum[14] = new TH1D("p_pizero"     , "Momentum of #pi^{0} ", 50,  0.0, 50.0);
			h_Momentum[15] = new TH1D("E_charged"    , "Energy   of charged ", 50,  0.0, 50.0);
			h_Momentum[16] = new TH1D("E_neutral"    , "Energy   of neutral ", 50,  0.0, 50.0);
			h_Momentum[17] = new TH1D("E_ISRgamma"   , "Energy   of ISRgamma", 150, 0.0, 150.0);
			h_Momentum[18] = new TH1D("K_Mom_in_B"   , "P of K in B decays  ", 50,  0.0, 50.0);
			h_Momentum[19] = new TH1D("K_Mom_in_D"   , "P of K in D decays  ", 50,  0.0, 50.0);
			h_Momentum[20] = new TH1D("Pi_Mom_in_B"  , "P of Pi in B decays ", 50,  0.0, 50.0);
			h_Momentum[21] = new TH1D("Pi_Mom_in_D"  , "P of Pi in D decays ", 50,  0.0, 50.0);
			//
			h_Mass    [ 0] = new TH1D("M_gammagamma" , "Mass of 2#gamma     ", 350,  0.0, 0.700);
			h_Mass    [ 1] = new TH1D("M_pipi"       , "Mass of 2#pi        ", 350,  0.3, 1.000);
			h_Mass    [ 2] = new TH1D("M_pip"        , "Mass of #pi and p   ", 100,  0.9, 1.200);
			h_Mass    [ 3] = new TH1D("M_KK"         , "Mass of 2Kaon       ", 100,  0.9, 1.100);
			h_Mass    [ 4] = new TH1D("M_Jpsi"       , "Mass of 2lepton     ", 100,  2.9, 3.300);
			h_Mass    [ 5] = new TH1D("M_Upsilon"    , "Mass of 2lepton     ", 100,  9.0, 11.00);
			//
			h_CosTheta[ 0] = new TH1D("all_particles", "Cos of Theta        ", 200, -1.0,  1.00);
			h_CosTheta[ 1] = new TH1D("ISR_Photon"   , "Cos of Theta        ", 200, -1.0,  1.00);
			//
			h_Ntrks   [ 0] = new TH1D("Ntrks"        , "# of tracks         ", 150,  0.0,150.00);
			h_Ntrks   [ 1] = new TH1D("Ngamma"       , "# of photons        ", 150,  0.0,150.00);
			h_Ntrks   [ 2] = new TH1D("Nparticles"   , "# of particles      ", 150,  0.0,150.00);
			h_Ntrks   [ 3] = new TH1D("Nneutralhadron", "# of neutral hadons", 150,  0.0,150.00);

			//
			h_Vertex[ 0]   = new TH1D("D0"           , "Vertex of D0        ", 220, 0.,  2200);
			h_Vertex[ 1]   = new TH1D("D+"           , "Vertex of D+        ", 220, 0.,  2200);
			h_Vertex[ 2]   = new TH1D("Ds"           , "Vertex of Ds        ", 220, 0.,  2200);
			h_Vertex[ 3]   = new TH1D("B0"           , "Vertex of B0        ", 220, 0.,  2200);
			h_Vertex[ 4]   = new TH1D("B+"           , "Vertex of B+        ", 220, 0.,  2200);
			h_Vertex[ 5]   = new TH1D("Bc"           , "Vertex of Bs        ", 220, 0.,  2200);
			h_Vertex[ 6]   = new TH1D("Bs"           , "Vertex of Bc        ", 220, 0.,  2200);
		}
	}


	void LGFastMCProcessor::processRunHeader( LCRunHeader* run) { 
		_nRun++ ;
		//FreeDelAll(_ptrash);
		//FreeDelAll(_tracktrash);
		//FreeDelAll(_clustertrash);
	} 


	void LGFastMCProcessor::processEvent( LCEvent * evt ) { 
		
		//FreeDelAll(_ptrash);
		//FreeDelAll(_tracktrash);
		//FreeDelAll(_clustertrash);

		const LCCollection* mcpCol = evt->getCollection( _inputCollectionName ) ;

		LCCollectionVec * recVec = new LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE ) ;
		LCCollectionVec * parVec = new LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE ) ;
		Cluflag.setBit(LCIO::CHBIT_LONG);
		recVec->setFlag(Cluflag.getFlag());
		parVec->setFlag(Cluflag.getFlag());

		LCRelationNavigator relNav( LCIO::RECONSTRUCTEDPARTICLE , LCIO::MCPARTICLE ) ;


		vector<TLorentzVector> PhotonList, PiMinusList, PiPlusList, KPlusList, KMinusList, 
			PList, AntiPList, ElectronList, PositronList, MuPlusList, MuMinusList;

		vector<TLorentzVector> rawp4list;
		vector<double>         pdgidlist;
		if(m_luxury>0){
			rawp4list.clear();
			pdgidlist.clear();
		}
		double VisEn=0, ntrk=0, ngam=0, npar=0, nNeuHad=0;
		for(int i=0; i<mcpCol->getNumberOfElements() ; i++){

			MCParticle* mcp = dynamic_cast<MCParticle*> ( mcpCol->getElementAt( i ) ) ;

			int status      =  (mcp)->getGeneratorStatus(); 
			int nParents    = ((mcp)->getParents()).size();
			//int pdgid       =  (mcp)->getPDG();
			//int nDaughters  = ((mcp)->getDaughters()).size();
			int idGParent    = 0, nGGParents=0;
			if ( nParents>0) {
				idGParent = ((mcp)->getParents()[0])->getPDG();
				nGGParents= ((mcp)->getParents()[0])->getParents().size();
			}	
			int idGGParent  = 0;
			if ( nGGParents>0) {
				idGGParent = (((mcp)->getParents()[0])->getParents()[0])->getPDG();
			}	
			TVector3 v(mcp->getMomentum());
			double   En(mcp->getEnergy());
			double   theta  = v.Theta() ;  
			double   pmag   = v.Mag();
			double   mass   = mcp->getMass();
			double   charge = mcp->getCharge();
			//
			int    pdgcode = mcp->getPDG();
			int    pdgCODE = abs(mcp->getPDG());
			double gabe = pmag/(mass+1e-8); 
			double vv=0, xx[3]={0,0,0};
			xx[0]= (mcp)->getEndpoint()[0];
			xx[1]= (mcp)->getEndpoint()[1];
			xx[2]= (mcp)->getEndpoint()[2];
			vv=pow(xx[0]*xx[0]+xx[1]*xx[1]+xx[2]*xx[2],0.5)/gabe*1000;
			if( m_makeplots>0 ) {
				h_CosTheta[0]->Fill(cos(theta), 1.0);	
				if(pdgCODE == 421) h_Vertex[0]->Fill(vv,1.0);
				if(pdgCODE == 411) h_Vertex[1]->Fill(vv,1.0);
				if(pdgCODE == 431) h_Vertex[2]->Fill(vv,1.0);
				if(pdgCODE == 511) h_Vertex[3]->Fill(vv,1.0);
				if(pdgCODE == 521) h_Vertex[4]->Fill(vv,1.0);
				if(pdgCODE == 531) h_Vertex[5]->Fill(vv,1.0);
				if(pdgCODE == 541) h_Vertex[6]->Fill(vv,1.0);
				//
			   if(pdgcode == 22 && ( nParents==0 ||  (nParents>0 && idGParent ==22 &&nGGParents==0)) ) h_CosTheta[1]->Fill(cos(theta), 1.0);	
				if(pdgcode ==   211 ) h_Momentum[ 5]->Fill(pmag,1.0); 
				if(pdgcode ==  -211 ) h_Momentum[ 6]->Fill(pmag,1.0); 
				if(pdgcode ==   321 ) h_Momentum[ 7]->Fill(pmag,1.0); 
				if(pdgcode ==  -321 ) h_Momentum[ 8]->Fill(pmag,1.0); 
				if(pdgcode ==  2212 ) h_Momentum[ 9]->Fill(pmag,1.0); 
				if(pdgcode == -2212 ) h_Momentum[10]->Fill(pmag,1.0); 
				if(pdgcode ==  2112 ) h_Momentum[11]->Fill(pmag,1.0); 
				if(pdgcode == -2112 ) h_Momentum[12]->Fill(pmag,1.0); 
				if(pdgcode ==   313 ) h_Momentum[13]->Fill(pmag,1.0); 
				if(pdgcode ==   111 ) h_Momentum[14]->Fill(pmag,1.0); 
				Int_t bc1 = (abs(idGParent )/100)%10; 
				Int_t bc2 = (abs(idGGParent)/100)%10; 
				if(abs(pdgcode) ==  321 ) {
					if ( (nParents>0 && bc1==5) ||(nGGParents>0 && bc2 ==5 ) ) 
						h_Momentum[18]->Fill(pmag,1.0);
				}	
				if(abs(pdgcode) ==  321 ) {
					if ( (nParents>0 && bc1==4) ||(nGGParents>0 && bc2 ==4 ) ) 
						h_Momentum[19]->Fill(pmag,1.0);
				}	
				if(abs(pdgcode) ==  211 ) {
					if ( (nParents>0 && bc1==5) ||(nGGParents>0 && bc2 ==5 ) ) 
						h_Momentum[20]->Fill(pmag,1.0);
				}	
				if(abs(pdgcode) ==  211 ) {
					if ( (nParents>0 && bc1==4) ||(nGGParents>0 && bc2 ==4 ) ) 
						h_Momentum[21]->Fill(pmag,1.0);
				}	
			}
			if( status == 1 ) { // stable particles only 

				if ( fabs(cos(theta) > 1.001) ) continue;
				if ( pmag   <    _momentumCut ) continue;
				if( m_makeplots>0 ) { 
					if(pdgcode ==    22 ) h_Momentum[ 0]->Fill(pmag,1.0); 
					if(pdgcode ==    11 ) h_Momentum[ 1]->Fill(pmag,1.0); 
					if(pdgcode ==   -11 ) h_Momentum[ 2]->Fill(pmag,1.0); 
					if(pdgcode ==    13 ) h_Momentum[ 3]->Fill(pmag,1.0); 
					if(pdgcode ==   -13 ) h_Momentum[ 4]->Fill(pmag,1.0); 
					if(pdgcode != 22 && ( nParents>0 && charge != 0)) h_Momentum[15]->Fill(En, 1.0);	
					if(pdgcode != 22 && ( nParents>0 && charge == 0)) h_Momentum[16]->Fill(En, 1.0);	
					if(pdgcode == 22 && ( nParents==0 ||(nParents>0 && idGParent==22))) h_Momentum[17]->Fill(En, 1.0);	
				}
				ReconstructedParticle*  rec = 0 ;

				if( _factory != 0 ) 
					rec = _factory->createReconstructedParticle( mcp ) ;

				if( rec != 0 ) {

					VisEn += En;
					npar++;
					if ( fabs(charge)>0.01 ) ntrk++;
					if ( pdgcode==22       ) ngam++;
					if ( pdgcode!=22 && fabs(charge)<0.01      ) nNeuHad++;
					recVec->addElement( rec ) ;
					relNav.addRelation( rec , mcp ) ;
					if ( fabs(charge)>0.01 ) relNav.addRelation( rec->getTracks()[0]   , mcp ) ;
					relNav.addRelation( rec->getClusters()[0] , mcp ) ;
					TLorentzVector p4(rec->getMomentum(), rec->getEnergy());
					/*
					EVENT::TrackVec::const_iterator it_trk = (rec->getTracks()).begin();
					for(; it_trk != (rec->getTracks()).end() ; it_trk++)
						_tracktrash.push_back(*it_trk);
					EVENT::ClusterVec::const_iterator it_clu = (rec->getClusters()).begin();
					for(; it_clu != (rec->getClusters()).end() ; it_clu++)
						_clustertrash.push_back(*it_clu);
					*/
					EVENT::ReconstructedParticleVec::const_iterator it_par = (rec->getParticles()).begin();
					for(; it_par != (rec->getParticles()).end() ; it_par++){
						parVec->addElement( *it_par ) ;
						_ptrash.push_back(*it_par);
					}
					//
					if(m_luxury>0){
						rawp4list.push_back(p4);
						pdgidlist.push_back(pdgcode);
					}
					if( m_makeplots>0 ) { 
						if(pdgcode ==   22  ) PhotonList  .push_back(p4) ; 
						if(pdgcode ==   11  ) ElectronList.push_back(p4) ; 
						if(pdgcode ==   13  ) MuMinusList .push_back(p4) ; 
						if(pdgcode ==  -11  ) PositronList.push_back(p4) ; 
						if(pdgcode ==  -13  ) MuPlusList  .push_back(p4) ; 
						if(pdgcode ==   211 ) PiPlusList  .push_back(p4) ; 
						if(pdgcode ==  -211 ) PiMinusList .push_back(p4) ; 
						if(pdgcode ==   321 ) KPlusList   .push_back(p4) ; 
						if(pdgcode ==  -321 ) KMinusList  .push_back(p4) ; 
						if(pdgcode ==  2212 ) PList       .push_back(p4) ; 
						if(pdgcode == -2112 ) AntiPList   .push_back(p4) ; 
					}
				}
			}

		}
		//printf("VisEn = %10.2f\n",VisEn);
		recVec->setDefault   ( true   ) ; // only true works, false has track/cluster but without particle
		recVec->setSubset    ( false  ) ; // can make reasonble slcio file with PFO collection/LCRelation !!! 
		recVec->setTransient ( false  ) ;
		parVec->setDefault   ( true   ) ; // only true works, false has track/cluster but without particle
		parVec->setSubset    ( false  ) ; // can make reasonble slcio file with PFO collection/LCRelation !!! 
		parVec->setTransient ( false  ) ;

		evt->addCollection( recVec, _recoParticleCollectionName ) ;
		evt->addCollection( parVec, "dummyParticles" ) ;
		evt->addCollection( relNav.createLCCollection() , _mcTruthCollectionName ) ;
		//
		if( m_makeplots>0 ) {
         
			//printf("No of Photon is %4d\n", PhotonList.size());
			//printf("No of Kaon+  is %4d\n", KPlusList.size());
			//printf("No of Kaon-  is %4d\n", KMinusList.size());
			h_Ntrks   [ 0]->Fill(ntrk,1.0);
			h_Ntrks   [ 1]->Fill(ngam,1.0);
			h_Ntrks   [ 2]->Fill(npar,1.0);
			h_Ntrks   [ 3]->Fill(nNeuHad,1.0);
			if ( PhotonList.size()>0){
				for(unsigned int i=0; i<PhotonList.size()-1; i++){
					for(unsigned int j=i+1; j<PhotonList.size(); j++){
						h_Mass[0]->Fill((PhotonList[i]+PhotonList[j]).M());
					}	 
				}
			}

			for(unsigned int i=0; i<PiPlusList.size(); i++){
				for(unsigned int j=0; j<PiMinusList.size(); j++){
					h_Mass[1]->Fill((PiPlusList[i]+PiMinusList[j]).M());
				}	 
			} 
		
			for(unsigned int i=0; i<PiPlusList.size(); i++){
				for(unsigned int j=0; j<AntiPList.size(); j++){
					h_Mass[2]->Fill((PiPlusList[i]+AntiPList[j]).M());
				}	 
			}

			for(unsigned int i=0; i<PiMinusList.size(); i++){
				for(unsigned int j=0; j<PList.size(); j++){
					h_Mass[2]->Fill((PiMinusList[i]+PList[j]).M());
				}	 
			} 

			for(unsigned int i=0; i<KPlusList.size(); i++){
				for(unsigned int j=0; j<KMinusList.size(); j++){
					h_Mass[3]->Fill((KPlusList[i]+KMinusList[j]).M());
				}	 
			}

			for(unsigned int i=0; i<ElectronList.size(); i++){
				for(unsigned int j=0; j<PositronList.size(); j++){
					h_Mass[4]->Fill((ElectronList[i]+PositronList[j]).M());
					h_Mass[5]->Fill((ElectronList[i]+PositronList[j]).M());
				}	 
			} 

			for(unsigned int i=0; i<MuPlusList.size(); i++){
				for(unsigned int j=0; j<MuMinusList.size(); j++){
					h_Mass[4]->Fill((MuPlusList[i]+MuMinusList[j]).M());
					h_Mass[5]->Fill((MuPlusList[i]+MuMinusList[j]).M());
				}	 
			} 

		}	
		//
		if(m_luxury>0){
			m_ntp->fill4Momentum("idx1","raw_" , rawp4list, rawp4list.size());
			m_ntp->fillArray    ("idx2","pdgid", pdgidlist, pdgidlist.size());
			m_ntp->write();
		}
		_nEvt ++ ;
	}



	void LGFastMCProcessor::check( LCEvent * evt ) { 
            

	}


	void LGFastMCProcessor::end(){ 

		streamlog_out( MESSAGE4 )  << "LGFastMCProcessor::end()  " << name() 
			<< " processed " << _nEvt << " events in " << _nRun << " runs "
			<< std::endl ;

		if(m_luxury>0){
		}
		if( m_makeplots>0 ) { 
			for( int i=0; i<22; i++){
				h_Momentum[i]->Write();
				//
				TH1D *h1 = new TH1D("h1","",h_Momentum[i]->GetNbinsX(), h_Momentum[i]->GetXaxis()->GetXmin(), h_Momentum[i]->GetXaxis()->GetXmax());
				NameAxes(h1, (char*)"P(GeV/c)", (char*)"Entries/1.0GeV/c");
				TH1D *h2 = 0; 
				char filename[256], title[256];
				//
				sprintf(filename,"figs/%s", (h_Momentum[i]->GetName()));
				sprintf(title,"%s", (h_Momentum[i]->GetTitle()));
				PlotDataMC(filename, 
						h1           , (char*)"",
						h_Momentum[i], (char*)h_Momentum[i]->GetTitle(),
						h2           , (char*)"",
						h2           , (char*)"",
						h2           , (char*)"",
						true, false, false, title
					  );
				//
				sprintf(filename,"figs/%s_log", (h_Momentum[i]->GetName()));
				sprintf(title,"%s", (h_Momentum[i]->GetTitle()));
				PlotDataMC(filename, 
						h1           , (char*)"",
						h_Momentum[i], (char*)h_Momentum[i]->GetTitle(),
						h2           , (char*)"",
						h2           , (char*)"",
						h2           , (char*)"",
						true,  true, false, title
					  );
				delete h1;
				delete h_Momentum[i];
			}
			char filename[256], title[256];
			char * names[] = {(char*)"",(char*)"MC Truth", (char*)"Exp Fit"};
			//
			for( int i=0; i<7; i++){
				NameAxes(h_Vertex[i], (char*)"Vertex(#mu m)", (char*)"Entries / 10 #mu m");
				h_Vertex[i]->Fit("expo", "+");
				sprintf(filename,"figs/Vertex_%s",h_Vertex[i]->GetName());
				sprintf(title,"%s", (h_Vertex[i]->GetTitle()));
				PlotDataFit(filename, h_Vertex[i], (char*)h_Vertex[i]->GetTitle(), names, true, 0.65,0.60,0.90,0.80, title);
			}
			for( int i=0; i<7; i++){
				h_Vertex[i]->Write();
				delete h_Vertex[i];
			}

			for( int i=0; i<6; i++){
				h_Mass[i]->Write();
				//
				TH1D *h1 = new TH1D("h1","",h_Mass[i]->GetNbinsX(), h_Mass[i]->GetXaxis()->GetXmin(), h_Mass[i]->GetXaxis()->GetXmax());
				NameAxes(h1, (char*)"P(GeV/c)", (char*)"Entries/1.0GeV/c");
				TH1D *h2 = 0; 
				char filename[256], title[256];
				//
				sprintf(filename,"figs/%s", (h_Mass[i]->GetName()));
				sprintf(title,"%s", (h_Mass[i]->GetTitle()));
				PlotDataMC(filename, 
						h1           , (char*)"",
						h_Mass[i], (char*)h_Mass[i]->GetTitle(),
						h2           , (char*)"",
						h2           , (char*)"",
						h2           , (char*)"",
						true, false, false, title
					  );
				//
				sprintf(filename,"figs/%s_log", (h_Mass[i]->GetName()));
				sprintf(title,"%s", (h_Mass[i]->GetTitle()));
				PlotDataMC(filename, 
						h1           , (char*)"",
						h_Mass[i], (char*)h_Mass[i]->GetTitle(),
						h2           , (char*)"",
						h2           , (char*)"",
						h2           , (char*)"",
						true,  true, false, title
					  );
				delete h1;
				delete h_Mass[i];
			}
			for( int i=0; i<4; i++){
				h_Ntrks[i]->Write();
				//
				TH1D *h1 = new TH1D("h1","",h_Ntrks[i]->GetNbinsX(), h_Ntrks[i]->GetXaxis()->GetXmin(), h_Ntrks[i]->GetXaxis()->GetXmax());
				NameAxes(h1, (char*)"N_{particles}", (char*)"Entries/1.00");
				TH1D *h2 = 0; 
				char filename[256], title[256];
				//
				sprintf(filename,"figs/%s", (h_Ntrks[i]->GetName()));
				sprintf(title,"%s", (h_Ntrks[i]->GetTitle()));
				PlotDataMC(filename, 
						h1           , (char*)"",
						h_Ntrks[i], (char*)h_Ntrks[i]->GetTitle(),
						h2           , (char*)"",
						h2           , (char*)"",
						h2           , (char*)"",
						true, false, false, title
						);
				//
				sprintf(filename,"figs/%s_log", (h_Ntrks[i]->GetName()));
				sprintf(title,"%s", (h_Ntrks[i]->GetTitle()));
				PlotDataMC(filename, 
						h1           , (char*)"",
						h_Ntrks[i], (char*)h_Ntrks[i]->GetTitle(),
						h2           , (char*)"",
						h2           , (char*)"",
						h2           , (char*)"",
						true,  true, false, title
						);
				delete h1;
				delete h_Ntrks[i];
			}

			for( int i=0; i<2; i++){
				h_CosTheta[i]->Write();
				//
				TH1D *h1 = new TH1D("h1","",h_CosTheta[i]->GetNbinsX(), h_CosTheta[i]->GetXaxis()->GetXmin(), h_CosTheta[i]->GetXaxis()->GetXmax());
				NameAxes(h1, (char*)"cos#theta", (char*)"Entries/0.01");
				TH1D *h2 = 0; 
				char filename[256], title[256];
				//
				sprintf(filename,"figs/%s", (h_CosTheta[i]->GetName()));
				sprintf(title,"%s", (h_CosTheta[i]->GetTitle()));
				PlotDataMC(filename, 
						h1           , (char*)"",
						h_CosTheta[i], (char*)h_CosTheta[i]->GetTitle(),
						h2           , (char*)"",
						h2           , (char*)"",
						h2           , (char*)"",
						true, false, false, title
						);
				//
				sprintf(filename,"figs/%s_log", (h_CosTheta[i]->GetName()));
				sprintf(title,"%s", (h_CosTheta[i]->GetTitle()));
				PlotDataMC(filename, 
						h1           , (char*)"",
						h_CosTheta[i], (char*)h_CosTheta[i]->GetTitle(),
						h2           , (char*)"",
						h2           , (char*)"",
						h2           , (char*)"",
						true,  true, false, title
						);
				delete h1;
				delete h_CosTheta[i];
			}

		}
	}
}
//**********************************************
//**********************************************
int LGFastMCProcessor::FindParton(MCParticle *mcp)
{
	int ret = 0;

	if ( mcp->getParents().size() == 0 )      // original particle, no parent  
	{
		if ( abs(mcp->getPDG())<=6 ){
		  	return mcp->getPDG();
		}
		else return 0;
	} 
	else if ( mcp->getParents().size() == 1 ) // Only one parent  
	{
		MCParticle *mcpa = mcp->getParents()[0];
		int ngp = mcpa->getParents().size();
		//
		if( ngp == 1 )
		{
			if ( mcp->getPDG()==mcpa->getPDG() )
			{  
				// mother is the same particle
				ret = FindParton(mcpa);	
			}
			else if ( abs(mcp->getPDG())<=6 && ( abs(mcpa->getPDG())>=22 && abs(mcpa->getPDG())<=25 ) )
			{
				// parent is gamma/gluon/W/Z and itself is a quark 
				ret = mcp->getPDG();
			}
			else if ( mcp->getPDG()==21  && mcpa->getPDG()==25 )
			{
				// parent is Higgs and itself is a gluon 
				ret = mcp->getPDG();
			}
			else if ( abs(mcp->getPDG())<= 6 && mcpa->getPDG()==21 )
			{
				// parent is gluon and itself is a quark, this kind of object marked by +100   
				ret = 100+mcp->getPDG();
			}	
			else if ( abs(mcpa->getPDG())<= 6 && mcp->getPDG()==21 )
			{
				// parent is quark and itself is a gluon 
				ret = FindParton(mcpa);	
			}	
		} else {
			if ( abs(mcpa->getPDG())>= 81 && abs(mcpa->getPDG())<=100 ) {
				// parent is intermediate particle; need to check the MANY(>1) grandparents
				if ( abs(mcp->getPDG()) <= 6 || ( abs(mcp->getPDG()) >= 11 && abs(mcp->getPDG()) <= 16)  ) 
			{  
				// itself is a quark/lepton, jump over the intermediate paticle, search for the same quark/lepton, then ... 
				for ( unsigned int i=0; i<mcpa->getParents().size(); i++){
					MCParticle *mcpb = mcpa->getParents()[i];
					if ( mcpb->getPDG() == mcp->getPDG() ) {
						//printf( "here: %5d vs. %5d (%5d)\n", mcp->getPDG(),mcpb->getPDG(), mcpb->getParents()[0]->getPDG() );
						ret = FindParton(mcpb);
					}
				}
				// for checking some strangge cases 
				if ( ret == 0  && abs(mcp->getPDG()) <= 6 ) {
					for ( unsigned int i=0; i<mcpa->getParents().size(); i++ ) {
						MCParticle *mcpb = mcpa->getParents()[i];
						printf( "%5d vs. %5d, ", mcp->getPDG(),mcpb->getPDG() );
					}
					printf("\n");	
				}
			} 
			else 
			{  
				// itself is a hadron 
				// first; get all partons from the grandparents
				map<MCParticle *, TVector3 > partonmap; partonmap.clear();
				map<MCParticle *, TVector3>::iterator it;
				for ( unsigned int i=0; i<mcpa->getParents().size(); i++){
					MCParticle *mcpb = mcpa->getParents()[i];
					MCParticle *mcpo = NULL;//FindParton(mcpb);
					TVector3 v(mcpb->getMomentum());
					if ( mcpo != NULL) {
						it = partonmap.find(mcpo);
						if( it != partonmap.end() ) {
							v += it->second; // why add?
							//printf("the mother of %6d is %6d\n", (int)mcpb->getPDG(), (int)mcpo->getPDG());
						}
						partonmap[mcpo] = v;
					}
				}
				// then take the nearest one as the favorite parton
				double nearestangle = 2 * 3.14159;
				int selectedparton = -1, pdgid=0;
				int pSelected = 0;
				int j;
				TVector3 v2(mcp->getMomentum());
				for(j=0,it = partonmap.begin();it != partonmap.end(); it++,j++){
					if ( it->first ){
						double angle = v2.Angle(it->second);
						if(nearestangle > angle)
						{ 
							nearestangle = angle; selectedparton = j; pSelected = (it->first)->getPDG(); pdgid=(it->first)->getPDG();
						}
					}
				}
				ret = pSelected;
				if ( ret == 0) {
					//cout << partonmap.size() << " partons obtained." << endl;
					//cout << "Parton " << selectedparton << " selected. "<< nearestangle  << endl;
				}
				FreeAll(partonmap);
			}
		} else {	
			ret = FindParton(mcpa);	
			if( 0 && ret==0){
				for (unsigned int k=0 ; k < mcpa->getParents().size() ; k++ ) {
					MCParticle 	* part    = mcpa->getParents()[k];
					printf("parent of %10d is %10d (%3d)\n", mcpa->getPDG(), part->getPDG(), (int)part->getParents().size() );
					for (unsigned int l=0 ; l < part->getParents().size() ; l++ ) {
						printf("parent %2d is %10d\n", l, part->getParents()[l]->getPDG() );
					}
				}
			}

		}	
	
		}
	}
	else // More than one parent  
	{ 
		if ( abs(mcp->getPDG())<=6  ){
			int mcpo2 = 0;
			for( unsigned int i=0; i<mcp->getParents().size(); i++ ){
				MCParticle *mcpa = mcp->getParents()[i];
				if( mcp->getParents().size()==2 && abs(mcpa->getPDG()) == 11 ){
					mcpo2 = mcp->getPDG();
					break;
				}else if( mcp->getPDG() == mcp->getPDG()  ){
					int mcpo = FindParton(mcpa);
					mcpo2 = mcpo;
					break;
				}
			}
			ret=mcpo2;
		}
		else
		{	
			int mcpo2 = 0;
			for( unsigned int i=0; i<mcp->getParents().size(); i++ ){
				MCParticle *mcpa = mcp->getParents()[i];
				int mcpo = FindParton(mcpa);
				if(mcpo2 && mcpo != mcpo2)
				{cout << mcp->getPDG()<<" has Multiple parents other than intermidiates with different parton origin!" << endl;}
				mcpo2 = mcpo;
			}
			//
			ret = mcpo2;
		}
	}
	//
	return ret;
}

