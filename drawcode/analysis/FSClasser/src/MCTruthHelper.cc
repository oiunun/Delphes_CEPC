#include "MCTruthHelper.h"

const double MCTruthHelper::Pb[3] = {0.90, 0.09, 0.01};
const double MCTruthHelper::Pc[3] = {0.25, 0.70, 0.05};
const double MCTruthHelper::Pg[3] = {0.03, 0.06, 0.91};

MCTruthHelper::MCTruthHelper( LCCollection* mcParticleCol){

	initLatexCode();
	// initialize flags to -1
        	
	m_allMcPList.clear(); 
	m_EditedList.clear(); 
	m_McTopPList.clear(); 
	m_finalPList.clear();
	m_partonList.clear();
	m_DecayList .clear();
	m_SmearList .clear();
	m_Rawp4List .clear();
	m_PDGTGList .clear();
	m_hDaughters.clear();

	m_nFSRGamma        = -1;
	m_TotalEnergy      = -1;
	m_MissingEnergy    = -1;
	m_nElectron        = -1;
	m_nMuon            = -1;
	m_nTau             = -1;
	m_nNuE             = -1;
	m_nNuMu            = -1;
	m_nNuTau           = -1;
	m_nPhoton          = -1;
	m_nFinalPhoton     = -1;
	TopoEdit(mcParticleCol);
	// make a vector of MCParticle's to use in other methods
	int _nMCP= mcParticleCol->getNumberOfElements();
	if ( 0 ) printf("total MC particles = %8d\n", _nMCP);
	for(int i1 = 0; i1 < _nMCP; i1++)
	{
		MCParticleImpl *mcp = dynamic_cast<MCParticleImpl *>(mcParticleCol->getElementAt(i1));
		//
		int status      =  (mcp)->getGeneratorStatus(); 
		int pdgid       =  (mcp)->getPDG();
		int nParents    = ((mcp)->getParents()).size();
		int nDaughters  = ((mcp)->getDaughters()).size();
		if ( nDaughters == 0 && status ==2     ) continue;
		if ( nParents   == 0 && abs(pdgid)==21 ) continue;
		int ParentID    = 0; 
		if( nParents >0 ) ParentID = ((mcp)->getParents()[0])->getPDG(); 
		//if ( pdgid  == 22 && ( ParentID == 111 || ParentID == 221)  ) continue;
		if ( 0 ) printf("trkidx = %3d, nParents = %1d (%8d), nDaughters= %2d, pdgid = %8d, mass=%9.4f\n", 
				i1, nParents, ParentID, nDaughters, pdgid, (mcp)->getMass());
		// MC particles to process further
		if(ParentID==25&&abs(pdgid)<25) m_hDaughters.push_back(mcp); 
		if(status != 1 && status != 2 ) continue; 
		if(     1              ) m_allMcPList       .push_back(mcp);
		if(  status     == 1   ) m_finalPList       .push_back(mcp);
		if(             ( nParents == 0 && nDaughters == 0 && abs(pdgid) < 20 ) || 
				( nParents == 0 && nDaughters == 1 && abs(pdgid) <= 6 ) || 
				( (abs(ParentID)==23 || abs(ParentID)==24 || abs(ParentID)==25) && abs(pdgid)<23 ) 
		  ) { 
			if( abs(pdgid) <= 6 ) m_partonList       .push_back(mcp);
			if( pdgid ==  21    ) m_partonList       .push_back(mcp);
			if( pdgid ==  21    ) m_GluonPList       .push_back(mcp);
			if( pdgid ==  +1    ) m_UpPList          .push_back(mcp);
			if( pdgid ==  -1    ) m_AntiUpPList      .push_back(mcp);
			if( pdgid ==  +2    ) m_DownPList        .push_back(mcp);
			if( pdgid ==  -2    ) m_AntiDownPList    .push_back(mcp);
			if( pdgid ==  +3    ) m_StrangePList     .push_back(mcp);
			if( pdgid ==  -3    ) m_AntiStrangePList .push_back(mcp);
			if( pdgid ==  +4    ) m_CharmPList       .push_back(mcp);
			if( pdgid ==  -4    ) m_AntiCharmPList   .push_back(mcp);
			if( pdgid ==  +5    ) m_BottomPList      .push_back(mcp);
			if( pdgid ==  -5    ) m_AntiBottomPList  .push_back(mcp);
			if( pdgid ==  +6    ) m_TopPList         .push_back(mcp);
			if( pdgid ==  -6    ) m_AntiTopPList     .push_back(mcp);
			if( pdgid == +11    ) m_ElectronPList    .push_back(mcp);
			if( pdgid == -11    ) m_PositronPList    .push_back(mcp);
			if( pdgid == +13    ) m_MuonPList        .push_back(mcp);
			if( pdgid == -13    ) m_AntiMuonPList    .push_back(mcp);
			if( pdgid == +15    ) m_TauPList         .push_back(mcp);
			if( pdgid == -15    ) m_AntiTauPList     .push_back(mcp);
		}
	}
	m_nFSRGamma        = nParticles(kFSRGamma);
	m_nPhoton          = nParticles(kGamma   );
	m_nFinalPhoton     = nParticles(kGamma   );

	m_nHiggs           = nParticles(kHiggs   );
	m_nZ0              = nParticles(kZ0      );
	m_nWp              = nParticles(kWp      );
	m_nWm              = nParticles(kWm      );
	
	m_nElectron        = nParticles(kEp)    + nParticles(kEm);
	m_nMuon            = nParticles(kMup)   + nParticles(kMum);
	m_nTau             = nParticles(kTaup)  + nParticles(kTaum);
	m_nNuE             = nParticles(kNuE)   + nParticles(kAntiNuE);
	m_nNuMu            = nParticles(kNuMu)  + nParticles(kAntiNuMu);
	m_nNuTau           = nParticles(kNuTau) + nParticles(kAntiNuTau);

	m_nPi              = nParticles(kPip)   + nParticles(kPim  );
	m_nKaon            = nParticles(kKp )   + nParticles(kKm   );
	m_nProton          = nParticles(kPp )   + nParticles(kPm   );
	m_nNeutron         = nParticles(kN  )   + nParticles(kAntiN);
	
	m_nBottom          = m_BottomPList.size() + m_AntiBottomPList.size();
	m_nCharm           = m_CharmPList.size()  + m_AntiCharmPList.size();
	m_nLight     
		= m_StrangePList.size() + m_AntiStrangePList.size()
		+ m_DownPList.size()    + m_AntiDownPList.size()
		+ m_UpPList.size()      + m_AntiUpPList.size()
		+ m_GluonPList.size()  ;

}

double
MCTruthHelper::MCMissingEnergy(){
	if (m_TotalEnergy > 0){
	}
	else {
		MCTotalEnergy();
	}
	return m_MissingEnergy;
}



double
MCTruthHelper::MCTotalEnergy(){
	if (m_TotalEnergy > 0){}
	else { 
		m_MissingEnergy = 0.0;
		m_TotalEnergy   = 0.0;
		for (unsigned int i = 0; i < m_finalPList.size(); i++){
			if (isParticle(m_finalPList[i], kNuE        ) ) {m_MissingEnergy += m_finalPList[i]->getEnergy(); continue;}
			if (isParticle(m_finalPList[i], kNuMu       ) ) {m_MissingEnergy += m_finalPList[i]->getEnergy(); continue;}
			if (isParticle(m_finalPList[i], kNuTau      ) ) {m_MissingEnergy += m_finalPList[i]->getEnergy(); continue;}
			if (isParticle(m_finalPList[i], kAntiNuE    ) ) {m_MissingEnergy += m_finalPList[i]->getEnergy(); continue;}
			if (isParticle(m_finalPList[i], kAntiNuMu   ) ) {m_MissingEnergy += m_finalPList[i]->getEnergy(); continue;}
			if (isParticle(m_finalPList[i], kAntiNuTau  ) ) {m_MissingEnergy += m_finalPList[i]->getEnergy(); continue;}
			m_TotalEnergy += m_finalPList[i]->getEnergy();
		}
	}
	return m_TotalEnergy;
}

TLorentzVector MCTruthHelper::getZP4( LCCollection * mcParticleCol){
	MCParticleImpl *ret =NULL;
	int _nMCP= mcParticleCol->getNumberOfElements();
	int id0=99999999, id1=0;
	for(int i1 = _nMCP-1; i1>=0; i1--)
	{
		MCParticleImpl *mcp = dynamic_cast<MCParticleImpl*>(mcParticleCol->getElementAt(i1));
		int pdgid  =  mcp->getPDG();
		if( i1==0) id1 = mcp->id();
		int idd = mcp->id()-id1;
		if( pdgid==25 && idd < id0){
			ret = mcp;
			id0=idd;
		}
	}

	TLorentzVector p4 = TLorentzVector( 0,0,0,0 ), pnu(0,0,0,0);
	if( ret ) {
		p4 = TLorentzVector(0,0,0,250) - TLorentzVector( ret->getMomentum(), ret->getEnergy() );
		for(int i1 = 0; i1 < _nMCP; i1++)
		{
			MCParticleImpl *mcp = dynamic_cast<MCParticleImpl*>(mcParticleCol->getElementAt(i1));
			int mcPID  = mcp->getPDG();
			int status = mcp->getGeneratorStatus();
			if( status ==1 && ( abs(mcPID)==12 || abs(mcPID)==14 || abs(mcPID)==16) && ! FindHiggsMother(mcp)) {
				p4 -= TLorentzVector( mcp->getMomentum(), mcp->getEnergy() );
				pnu+= TLorentzVector( mcp->getMomentum(), mcp->getEnergy() );
			}	
		}
	}
	_E_nu_Z = pnu.E();
	return p4;
}

TLorentzVector MCTruthHelper::getHiggsP4( LCCollection * mcParticleCol){
	MCParticleImpl *ret =NULL;
	int _nMCP= mcParticleCol->getNumberOfElements();
	int id0=99999999, id1=0;
	for(int i1 = _nMCP-1; i1>=0; i1--)
	{
		MCParticleImpl *mcp = dynamic_cast<MCParticleImpl*>(mcParticleCol->getElementAt(i1));
		int pdgid  =  mcp->getPDG();
		if( i1==0) id1 = mcp->id();
		int idd = mcp->id()-id1;
		if( pdgid==25 && idd < id0){
			ret = mcp;
			id0=idd;
		}
	}

	TLorentzVector p4 = TLorentzVector( 0,0,0,0 ), pnu(0,0,0,0);
	if( ret ) {
		p4 = TLorentzVector( ret->getMomentum(), ret->getEnergy() );
		for(int i1 = 0; i1 < _nMCP; i1++)
		{
			MCParticleImpl *mcp = dynamic_cast<MCParticleImpl*>(mcParticleCol->getElementAt(i1));
			int mcPID  = mcp->getPDG();
			int status = mcp->getGeneratorStatus();
			if( status == 1 && ( abs(mcPID)==12 || abs(mcPID)==14 || abs(mcPID)==16) && FindHiggsMother(mcp)) {
				//cout<<"I am here "<<i1<<" "<<status<<" "<<mcPID<<" "<< FindHiggsMother(mcp)<<endl;
				p4 -= TLorentzVector( mcp->getMomentum(), mcp->getEnergy() );
				pnu+= TLorentzVector( mcp->getMomentum(), mcp->getEnergy() );
			}	
		}
	}
	_E_nu_H = pnu.E();
	return p4;
}
void MCTruthHelper::TopoEdit(LCCollection * mcParticleCol){
	
	int _nMCP= mcParticleCol->getNumberOfElements();
	vector<int> daughterList, parentList; 
	for(int i1 = 0; i1 < _nMCP; i1++)
	{
		MCParticleImpl *mcp = dynamic_cast<MCParticleImpl*>(mcParticleCol->getElementAt(i1));
		int ID = mcp->id();
		unsigned int nParents     = (mcp->getParents()).size();
		if ( nParents<2) continue;
		parentList.clear();
		for( unsigned int id = 0; id< nParents; id++){
			parentList.push_back((mcp->getParents()[id])->getPDG()); 
		}
		unsigned int nDaughters   = (mcp->getDaughters()).size();
		for( unsigned int id = 0; id< nDaughters; id++){
			daughterList.push_back((mcp->getDaughters()[id])->getPDG()); 
		}
	
		int pdgid          =  mcp->getPDG();
		int status         =  mcp->getGeneratorStatus(); 
		TLorentzVector  p4 =  TLorentzVector( mcp->getMomentum(), mcp->getEnergy() );
		if( 0 )
			printf("trkidx = %8d, nParents = %1d,  nDaughters = %2d, pdgid = %8d, status = %2d, mass = %7.3f, E = %7.3f\n", 
					ID, nParents, nDaughters, pdgid, status, fabs(p4.M()), p4.T());
	
	
	}
}

MCTruthHelper::~MCTruthHelper(){
	FreeAll( m_allMcPList      );
	FreeAll( m_EditedList      );
	FreeAll( m_finalPList      );
	FreeAll( m_partonList      );
	FreeAll( m_DecayList       );
	FreeAll( m_SmearList       );
	FreeAll( m_Rawp4List       );
	FreeAll( m_PDGTGList       );
	FreeAll( m_UpPList         );
	FreeAll( m_AntiUpPList     );
	FreeAll( m_DownPList       );
	FreeAll( m_AntiDownPList   );
	FreeAll( m_StrangePList    );
	FreeAll( m_AntiStrangePList);
	FreeAll( m_CharmPList      );
	FreeAll( m_AntiCharmPList  );
	FreeAll( m_BottomPList     );
	FreeAll( m_AntiBottomPList );
	FreeAll( m_TopPList        );
	FreeAll( m_AntiTopPList    );
	FreeAll( m_ElectronPList   );
	FreeAll( m_PositronPList   );
	FreeAll( m_MuonPList       );
	FreeAll( m_AntiMuonPList   );
	FreeAll( m_TauPList        );
	FreeAll( m_AntiTauPList    );
}

bool
MCTruthHelper::isParticle(const MCParticleImpl* mcParticle, int idParticle){
	// if the particle id's don't match, return false
	if (mcParticle->getPDG() != idParticle) return false;
	return true;
}

int
MCTruthHelper::nHiggsDaughters(int idParticle, int idZ){
	int n = 0;
	for (unsigned int i = 0; i < m_hDaughters.size(); i++){
                if( idZ!=0 ){
			if (isParticle(m_hDaughters[i], idParticle)) n += idParticle;
			if (isParticle(m_hDaughters[i], idZ       )) n += idZ;
		}else{
			if (isParticle(m_hDaughters[i], idParticle)) n += idParticle;
			if (isParticle(m_hDaughters[i],-idParticle)) n += idParticle;
		}

	}
	return n;
}

int
MCTruthHelper::nHiggsFinalState(){

	int nHiggs[15];
        //              null  uu    dd    ss	cc    bb     tt     ee   mumu   tautau  gluongluon  gammagamma     gamma Z   ZZ  WW
	int nHSumD[15]={0,    2,    4,    6,     8,   10,    12,    22,     26,     30,         42,         44,          45, 46, 48}; 
	
	nHiggs[ 0] = 0; 
	nHiggs[ 1] = nHiggsDaughters(  1); 	
	nHiggs[ 2] = nHiggsDaughters(  2); 	
	nHiggs[ 3] = nHiggsDaughters(  3); 	
	nHiggs[ 4] = nHiggsDaughters(  4); 	
	nHiggs[ 5] = nHiggsDaughters(  5); 	
	nHiggs[ 6] = nHiggsDaughters(  6); 

	nHiggs[ 7] = nHiggsDaughters( 11); 	
	nHiggs[ 8] = nHiggsDaughters( 13); 	
	nHiggs[ 9] = nHiggsDaughters( 15); 

	nHiggs[10] = nHiggsDaughters( 21); 	
	nHiggs[11] = nHiggsDaughters( 22); 	
	nHiggs[12] = nHiggsDaughters( 22, 23); 	
	nHiggs[13] = nHiggsDaughters( 23); 	
	nHiggs[14] = nHiggsDaughters( 24); 	
	//	
	int n=0;
	for( int i=1; i< 15; i++){
		if( nHiggs[i] == 0        ) continue; 
		if( nHiggs[i] == nHSumD[i]) {
			n=i;
			break;
		}
	}
	//	
	return n;
}




int
MCTruthHelper::nParticles(int idParticle){
	int n = 0;
	for (unsigned int i = 0; i < m_finalPList.size(); i++){
		if (isParticle(m_finalPList[i],idParticle)) n++;
	}
	if (n > 1000){
		cout<< "MCTruthHelper WARNING:  found " << n << " generated particles "
			 << "with ID = " << idParticle << " but returning 1000 " << endl;
		return 1000;
	}
	return n;
}

bool
MCTruthHelper::hasParticle(int idParticle){
	return (nParticles(idParticle) > 0);
}


vector<MCParticle*> 
MCTruthHelper::getDaughters(const MCParticleImpl* P){
	vector<MCParticle*> daughters; daughters.clear(); 
	int isLEAF = (P->getDaughters()).size();
	if ( isLEAF ) return daughters;	
	int  Did  =  abs(P->getPDG());
	bool Dmes = (Did==kD0 || Did==kDp || Did==kDsp);
	//
	vector<MCParticle*> daughtersCheck = P->getDaughters();
	for (unsigned int j = 0; j < daughtersCheck.size(); j++){
		int      leaf = daughtersCheck[j]->getDaughters().size(); 
		int        id = daughtersCheck[j]->getPDG();
		MCParticle* M = daughtersCheck[j]->getParents()[0]; 
		int       idM = M->getPDG();
		if (id == kFSRGamma) continue;
		if( leaf ) {
			daughters.push_back(daughtersCheck[j]); 
			continue;
		}
		// 
		if (id == kGamma && (idM=kPi0 || idM==kEta ))                                  continue;
		if (id == kGamma && Dmes && (idM==kEtaprime && hasDaughters(M, kGamma, kRho0)))continue;
		//
		if ( id==kPi0 && hasDaughters(daughtersCheck[j],kGamma,kGamma)){
			daughters.push_back(daughtersCheck[j]); 
			continue;
		}
		if ( id==kEta && hasDaughters(daughtersCheck[j],kGamma,kGamma)){
			daughters.push_back(daughtersCheck[j]); 
			continue;
		}
		if (id==kEtaprime&& Dmes &&
				(hasDaughters(daughtersCheck[j],kPip,kPim,kEta)||hasDaughters(daughtersCheck[j],kGamma,kRho0))){
			daughters.push_back(daughtersCheck[j]); // only for D tag cases,  potential bug for other usage  
			continue;
		}
		if ( id==kKs && hasDaughters(daughtersCheck[j],kPip,kPim)){
			daughters.push_back(daughtersCheck[j]); 
			continue;
		}
		if ( id==kLambda && hasDaughters(daughtersCheck[j],kPp,kPim)){
			daughters.push_back(daughtersCheck[j]); 
			continue;
		}
		if ( id==kALambda && hasDaughters(daughtersCheck[j],kPm,kPip)){
			daughters.push_back(daughtersCheck[j]); 
			continue;
		}
		//
		vector<MCParticle*> VP=	getDaughters((MCParticleImpl*)daughtersCheck[j]);
		for (unsigned int k = 0; k < VP.size(); k++){
			daughters.push_back(VP[k]); 
		}
	}	
	return daughters;	
}

bool 
MCTruthHelper::hasMothers(const MCParticleImpl* P, const int mid){
	MCParticle* M = ( MCParticle*)P;
	Int_t np = M->getParents().size();
	//
	while( np>0 && M->getPDG() != M->getParents()[0]->getPDG() ){
		int id = abs(M->getParents()[0]->getPDG()); 
		if ( id==mid ) return true;
		np = M->getParents().size();
		if( np>0 ) M=M->getParents()[0];
	        else return false;	
		np = M->getParents().size();
	}
	//
	return false;
}

bool
MCTruthHelper::DhasDaughters(
		const MCParticleImpl* mcParticle, 
		int idDaughter1,
		int idDaughter2,
		int idDaughter3,
		int idDaughter4,
		int idDaughter5,
		int idDaughter6,
		int idDaughter7,
		int idDaughter8,
		int idDaughter9,
		int idDaughter10
		){
	// create a vector of daughter id's to search for
	vector<int> idDaughters;
	if (idDaughter1  != 0)  idDaughters.push_back( idDaughter1);
	if (idDaughter2  != 0)  idDaughters.push_back( idDaughter2);
	if (idDaughter3  != 0)  idDaughters.push_back( idDaughter3);
	if (idDaughter4  != 0)  idDaughters.push_back( idDaughter4);
	if (idDaughter5  != 0)  idDaughters.push_back( idDaughter5);
	if (idDaughter6  != 0)  idDaughters.push_back( idDaughter6);
	if (idDaughter7  != 0)  idDaughters.push_back( idDaughter7);
	if (idDaughter8  != 0)  idDaughters.push_back( idDaughter8);
	if (idDaughter9  != 0)  idDaughters.push_back( idDaughter9);
	if (idDaughter10 != 0)  idDaughters.push_back(idDaughter10);
	if (idDaughters.size() < 2) return false;

	// create a vector of daughter id's from the MCParticleImpl
	vector<int> idDaughtersCheck;
	vector<MCParticle*> daughtersCheck = getDaughters(mcParticle);
	for (unsigned int j = 0; j < daughtersCheck.size(); j++){
		int id = (*(daughtersCheck[j])).getPDG();    
		if (id != kFSRGamma) idDaughtersCheck.push_back(id);
	}

	// if the numbers of daughters don't match return false
	if (idDaughters.size() != idDaughtersCheck.size()) return false;

	// sort the vectors of daughter id's
	for (unsigned int i1 = 0; i1 < idDaughters.size()-1; i1++){
		for (unsigned int i2 = i1+1; i2 < idDaughters.size(); i2++){
			if (idDaughters[i1] > idDaughters[i2]){
				int temp = idDaughters[i1];
				idDaughters[i1] = idDaughters[i2];
				idDaughters[i2] = temp;
			}
			if (idDaughtersCheck[i1] > idDaughtersCheck[i2]){
				int temp = idDaughtersCheck[i1];
				idDaughtersCheck[i1] = idDaughtersCheck[i2];
				idDaughtersCheck[i2] = temp;
			}
		}
	}

	// check if the daughter id's match
	for (unsigned int i1 = 0; i1 < idDaughters.size(); i1++)
	{
		if (idDaughters[i1] != idDaughtersCheck[i1]) return false;
	}
	return true;
}


bool
MCTruthHelper::hasDaughters(
		const MCParticle* mcParticle, 
		int idDaughter1,
		int idDaughter2,
		int idDaughter3,
		int idDaughter4,
		int idDaughter5,
		int idDaughter6,
		int idDaughter7,
		int idDaughter8,
		int idDaughter9,
		int idDaughter10
		){

	// create a vector of daughter id's to search for
	vector<int> idDaughters;
	if (idDaughter1  != 0)  idDaughters.push_back( idDaughter1);
	if (idDaughter2  != 0)  idDaughters.push_back( idDaughter2);
	if (idDaughter3  != 0)  idDaughters.push_back( idDaughter3);
	if (idDaughter4  != 0)  idDaughters.push_back( idDaughter4);
	if (idDaughter5  != 0)  idDaughters.push_back( idDaughter5);
	if (idDaughter6  != 0)  idDaughters.push_back( idDaughter6);
	if (idDaughter7  != 0)  idDaughters.push_back( idDaughter7);
	if (idDaughter8  != 0)  idDaughters.push_back( idDaughter8);
	if (idDaughter9  != 0)  idDaughters.push_back( idDaughter9);
	if (idDaughter10 != 0)  idDaughters.push_back(idDaughter10);
	if (idDaughters.size() < 2) return false;

	// create a vector of daughter id's from the MCParticle
	vector<int> idDaughtersCheck;
	vector<MCParticle*> daughtersCheck = mcParticle->getDaughters();
	for (unsigned int j = 0; j < daughtersCheck.size(); j++){
		int id = ((daughtersCheck[j]))->getPDG();    
		if (id != kFSRGamma) idDaughtersCheck.push_back(id);
	}

	// if the numbers of daughters don't match return false
	if (idDaughters.size() != idDaughtersCheck.size()) return false;

	// sort the vectors of daughter id's
	for (unsigned int i1 = 0; i1 < idDaughters.size()-1; i1++){
		for (unsigned int i2 = i1+1; i2 < idDaughters.size(); i2++){
			if (idDaughters[i1] > idDaughters[i2]){
				int temp = idDaughters[i1];
				idDaughters[i1] = idDaughters[i2];
				idDaughters[i2] = temp;
			}
			if (idDaughtersCheck[i1] > idDaughtersCheck[i2]){
				int temp = idDaughtersCheck[i1];
				idDaughtersCheck[i1] = idDaughtersCheck[i2];
				idDaughtersCheck[i2] = temp;
			}
		}
	}

	// check if the daughter id's match
	for (unsigned int i1 = 0; i1 < idDaughters.size(); i1++)
	{
		if (idDaughters[i1] != idDaughtersCheck[i1]) return false;
	}
	return true;
}


int
MCTruthHelper::nDecays(
		int  idDaughter1,
		int  idDaughter2,
		int  idDaughter3,
		int  idDaughter4,
		int  idDaughter5,
		int  idDaughter6,
		int  idDaughter7,
		int  idDaughter8,
		int  idDaughter9,
		int idDaughter10
		){
	int n = 0;
	for (unsigned int i = 0; i < m_finalPList.size(); i++){
		MCParticleImpl* mcp = m_finalPList[i];
		if (hasDaughters(mcp,    
					idDaughter1,
					idDaughter2,
					idDaughter3,
					idDaughter4,
					idDaughter5,
					idDaughter6,
					idDaughter7,
					idDaughter8,
					idDaughter9,
					idDaughter10
					)) n++;
	}
	if (n > 14){
		cout << "FSClasser WARNING:  found " << n << " generated decays with "
			<< "daughters = " 
			<< idDaughter1  << " "
			<< idDaughter2  << " "
			<< idDaughter3  << " "
			<< idDaughter4  << " "
			<< idDaughter5  << " "
			<< idDaughter6  << " "
			<< idDaughter7  << " "
			<< idDaughter8  << " "
			<< idDaughter9  << " "
			<< idDaughter10 << " "
			<< " but returning 14 " << endl;
		return 14;
	}
	return n;
}


bool
MCTruthHelper::hasDecay(
		int idDaughter1,
		int idDaughter2,
		int idDaughter3,
		int idDaughter4,
		int idDaughter5,
		int idDaughter6,
		int idDaughter7,
		int idDaughter8,
		int idDaughter9,
		int idDaughter10
		){
	return (nDecays(
				idDaughter1,
				idDaughter2,
				idDaughter3,
				idDaughter4,
				idDaughter5,
				idDaughter6,
				idDaughter7,
				idDaughter8,
				idDaughter9,
				idDaughter10
				) > 0);
}



int
MCTruthHelper::nVertices(int idParent, 
		int idDaughter1,
		int idDaughter2,
		int idDaughter3,
		int idDaughter4,
		int idDaughter5,
		int idDaughter6,
		int idDaughter7,
		int idDaughter8,
		int idDaughter9,
		int idDaughter10
		){
	int n = 0;
	for (unsigned int i = 0; i < m_finalPList.size(); i++){
		MCParticleImpl* mcp = m_finalPList[i];
		if (hasMothers(mcp)) continue;
		if (mcp->getPDG() == idParent && 
				hasDaughters(mcp,
					idDaughter1,
					idDaughter2,
					idDaughter3,
					idDaughter4,
					idDaughter5,
					idDaughter6,
					idDaughter7,
					idDaughter8,
					idDaughter9,
					idDaughter10
					)) n++;
	}
	if (n > 14){
		cout << "FSClasser WARNING:  found " << n << " generated vertices with "
			<< "parent = " << idParent << " and "
			<< "daughters = " 
			<< idDaughter1  << " "
			<< idDaughter2  << " "
			<< idDaughter3  << " "
			<< idDaughter4  << " "
			<< idDaughter5  << " "
			<< idDaughter6  << " "
			<< idDaughter7  << " "
			<< idDaughter8  << " "
			<< idDaughter9  << " "
			<< idDaughter10 << " "
			<< " but returning 14" << endl;
		return 14;
	}
	return n;
}

int
MCTruthHelper::dVertices(
		int idParent, 
		int idDaughter1,
		int idDaughter2,
		int idDaughter3,
		int idDaughter4,
		int idDaughter5,
		int idDaughter6,
		int idDaughter7,
		int idDaughter8,
		int idDaughter9,
		int idDaughter10
		){

	int n = 0;
	for (unsigned int i = 0; i < m_finalPList.size(); i++){
		MCParticleImpl* mcp = m_finalPList[i];
		if (mcp->getPDG() == idParent && 
				DhasDaughters(mcp,
					idDaughter1 ,
					idDaughter2 ,
					idDaughter3 ,
					idDaughter4 ,
					idDaughter5 ,
					idDaughter6 ,
					idDaughter7 ,
					idDaughter8 ,
					idDaughter9 ,
					idDaughter10)
				) n++;
	}
	
	if (n > 14){
		cout << "FSClasser WARNING:  found " << n << " generated vertices with "
			<< "parent = " << idParent << " and "
			<< "daughters = " 
			<< idDaughter1  << " "
			<< idDaughter2  << " "
			<< idDaughter3  << " "
			<< idDaughter4  << " "
			<< idDaughter5  << " "
			<< idDaughter6  << " "
			<< idDaughter7  << " "
			<< idDaughter8  << " "
			<< idDaughter9  << " "
			<< idDaughter10 << " "
			<< " but returning 14" << endl;
		return 14;
	}
	return n;
}



bool
MCTruthHelper::hasVertex(
		int idParent, 
		int idDaughter1,
		int idDaughter2,
		int idDaughter3,
		int idDaughter4,
		int idDaughter5,
		int idDaughter6,
		int idDaughter7,
		int idDaughter8,
		int idDaughter9,
		int idDaughter10
		){
	return (nVertices(
				idParent,
				idDaughter1,
				idDaughter2,
				idDaughter3,
				idDaughter4,
				idDaughter5,
				idDaughter6,
				idDaughter7,
				idDaughter8,
				idDaughter9,
				idDaughter10
				) > 0);
}


string
MCTruthHelper::particleType(const MCParticleImpl* mcParticle){
	return particleType(mcParticle->getPDG());
}

string
MCTruthHelper::particleType(int id){

	string name("");

	if      (id == kPsi2S)           name = "psi(2S)";
	else if (id == kPsi3770)         name = "psi(3770)";
	else if (id == kGamma)           name = "gamma";
	else if (id == kFSRGamma)        name = "FSRgamma";
	else if (id == kZ0   )           name = "Z^0";
	else if (id == kCluster)         name = "Cluster";
	else if (id == kString )         name = "String";
	else if (id == kHc)              name = "h_c";
	else if (id == kChic0)           name = "chi_c0";
	else if (id == kChic1)           name = "chi_c1";
	else if (id == kChic2)           name = "chi_c2";
	else if (id == kChic0p)          name = "chi'_c0";
	else if (id == kChic1p)          name = "chi'_c1";
	else if (id == kChic2p)          name = "chi'_c2";
	else if (id == kJpsi)            name = "J/psi";
	else if (id == kEtac)            name = "eta_c";
	else if (id == kPhi)             name = "phi";
	else if (id == kOmega)           name = "omega";
	else if (id == kPi0)             name = "pi0";
	else if (id == kPip)             name = "pi+";
	else if (id == kPim)             name = "pi-";
	else if (id == kRho0)            name = "rho0";
	else if (id == kRhop)            name = "rho+";
	else if (id == kRhom)            name = "rho-";
	else if (id == kA00 )            name = "a_00";
	else if (id == kA0p )            name = "a_0+";
	else if (id == kA0m )            name = "a_0-";
	else if (id == kB10 )            name = "b_10";
	else if (id == kB1p )            name = "b_1+";
	else if (id == kB1m )            name = "b_1-";
	else if (id == kA10 )            name = "a_10";
	else if (id == kA1p )            name = "a_1+";
	else if (id == kA1m )            name = "a_1-";
	else if (id == kF01370)          name = "f0(1370)";
	else if (id == kF01500)          name = "f0(1500)";
	else if (id == kH1  )            name = "h_1";
	else if (id == kH1p )            name = "h'_1";
	else if (id == kA20 )            name = "a_20";
	else if (id == kA2p )            name = "a_2+";
	else if (id == kA2m )            name = "a_2-";
	else if (id == kEtaprime)        name = "etaprime";
	else if (id == kEta)             name = "eta";
	else if (id == kKs)              name = "K_S0";
	else if (id == kKl)              name = "K_L0";
	else if (id == kKp)              name = "K+";
	else if (id == kKm)              name = "K-";
	else if (id == kPp)              name = "p+";
	else if (id == kPm)              name = "p-";
	else if (id == kN)               name = "N";
	else if (id == kAntiN)           name = "anti-N";
	else if (id == kEp)              name = "e+";
	else if (id == kEm)              name = "e-";
	else if (id == kMup)             name = "mu+";
	else if (id == kMum)             name = "mu-";
	else if (id == kTaup)            name = "tau+";
	else if (id == kTaum)            name = "tau-";
	else if (id == kNuE)             name = "nu";
	else if (id == kNuMu)            name = "nu";
	else if (id == kNuTau)           name = "nu";
	else if (id == kAntiNuE)         name = "nu";
	else if (id == kAntiNuMu)        name = "nu";
	else if (id == kAntiNuTau)       name = "nu";
	else if (id == kF0600)           name = "f0(600)";
	else if (id == kK0)              name = "K0";
	else if (id == kAntiK0)          name = "K0";
	else if (id == kKstarp)          name = "K*+";
	else if (id == kKstarm)          name = "K*-";
	else if (id == kKstar0)          name = "K*0";
	else if (id == kAntiKstar0)      name = "K*0";
	else if (id == kK0star0)         name = "K_0*0";
	else if (id == kAntiK0star0)     name = "Anti-K_0*0";
	else if (id == kK0starp)         name = "K_0*+";
	else if (id == kK0starm)         name = "K_0*-";
	else if (id == kK10   )          name = "K_10";
	else if (id == kK1p   )          name = "K_1+";
	else if (id == kK1m   )          name = "K_1-";
	else if (id == kAntiK10)         name = "anti-K_10";
	else if (id == kLambda)          name = "Lambda";
	else if (id == kALambda)         name = "anti-Lambda";
	else if (id == kD0     )         name = "D0";
	else if (id == kD0bar  )         name = "D0bar";
	else if (id == kDp     )         name = "D+";
	else if (id == kDm     )         name = "D-";
	else if (id == kDsp    )         name = "Ds+";
	else if (id == kDsm    )         name = "Ds-";
	else if (id == kDstarP )         name = "D*+";
	else if (id == kDstarM )         name = "D*-";
	else if (id == kDstar )          name = "D*";
	else if (id == kAntiDstar)       name = "anti-D*";
	else if (id == kDsstarP)         name = "D_s*+";
	else if (id == kDsstarM)         name = "D_s*-";
	else if (id == kY4260  )         name = "Y(4260)";
	else if (id == kY4360  )         name = "Y(4360)";
	else if (id == kPsi4040)         name = "psi(4040)";
	else if (id == kPsi4160)         name = "psi(4160)";
	else if (id == kPsi4415)         name = "psi(4415)";
	else if (id == kDeltapp)         name = "Delta++";
	else if (id == kAntiDeltapp)     name = "Delta--";
	else if (id == kSigma0 )         name = "Sigma0";
	else if (id == kAntiSigma0 )     name = "anti-Sigma0";
	else if (id == kSigmastarm )     name = "Sigma*+";
	else if (id == kAntiSigmastarp ) name = "anti-Sigma*-";
	else if (id == kXi0        )     name = "Xi0";
	else if (id == kAntiXi0        ) name = "anti-Xi0";
	else if (id == kXistarp    )     name = "Xi*+";
	else if (id == kAntiXistarm    ) name = "anti-Xi*-";
	else if (abs(id)<9             ) name = "quark";
	else{
		//iostringstream ss;ss << id;
		char ss[20]; sprintf(ss,"%d",id);
		name = ss;
	}
	return name;
}

double
MCTruthHelper::MomentumReso(MCParticleImpl * a_MCP)
{
	double Reso = 0;

	TVector3      currP  = a_MCP->getMomentum();
	int         currPID  = a_MCP->getPDG();
	double   PAmplitude  = currP.Mag();
	double     cosTheta  = fabs(a_MCP->getMomentum()[2]/PAmplitude);
	double     ObjkCoeff = 0; 
	double   Scalefactor = 0;  //Correct EndCap divergence

	if(abs(currPID) == 13 || abs(currPID) == 11 )
	{
		////////mumu_1365/////////////
		//       ObjkCoeff =1.59962e-03*PAmplitude+2.38933e-05*PAmplitude*PAmplitude;
		////////mumu_1465/////////////
		//       ObjkCoeff =1.48628e-03*PAmplitude+2.07520e-05*PAmplitude*PAmplitude;
		////////mumu_1565/////////////
		//       ObjkCoeff =1.32027e-03*PAmplitude+1.87041e-05*PAmplitude*PAmplitude;
		////////mumu_1665/////////////
		//       ObjkCoeff =1.18561e-03*PAmplitude+1.67721e-05*PAmplitude*PAmplitude;
		////////mumu_1808/////////////
		ObjkCoeff =1.05516e-03*PAmplitude+1.41823e-05*PAmplitude*PAmplitude;
		ObjkCoeff /= 2.0; 
      
		if(cosTheta > 0.86)     //Scale as effective R^2
		{
			Scalefactor = (1.0/(cosTheta*cosTheta) - 1)*2.96; //2.96 = (Half_Z/Radius)**2
			ObjkCoeff   = ObjkCoeff*1.0/Scalefactor; 
		}
		// Reso = ObjkCoeff*PAmplitude*PAmplitude;         // To be extended
		Reso = ObjkCoeff;       // To be extended

	}
	else if( abs(currPID) < 6 || abs(currPID) == 21 )
	{
		Reso = 0.06*PAmplitude;           // 6% of jet energy resolution  
	}
	else if( abs(currPID) == 15 )        // of coz need more input/analysis
	{
		Reso = 0.1*PAmplitude;
	}
	else
	{
		Reso = -10000;                    // Protection
	}
	return Reso; 
}

map<MCParticleImpl*, TLorentzVector> 
MCTruthHelper::SmearedList(){
	for (unsigned int i = 0; i < m_allMcPList.size(); i++){
		MCParticleImpl *a_MCP = m_allMcPList[i];
		double PID  = a_MCP->getPDG();
		double _SF  = 1.0;
		TLorentzVector p4(a_MCP->getMomentum()[0],a_MCP->getMomentum()[1],a_MCP->getMomentum()[2],a_MCP->getEnergy());
		if ( a_MCP->getCharge() != 0){
			double CurrPAmpli_p  = p4.Rho();
			double ResoAmpli_p   = MomentumReso(a_MCP);
			_SF  = (gRandom->Gaus(CurrPAmpli_p, ResoAmpli_p))/CurrPAmpli_p; 
		}else{
			double ResoAmpli_p = 0.16*sqrt(p4.E());
			_SF  = (gRandom->Gaus(p4.E(), ResoAmpli_p))/p4.E(); 
		}
		m_Rawp4List.insert( map<MCParticleImpl*, TLorentzVector, less<MCParticleImpl*> >::value_type(a_MCP, p4    ));
		m_SmearList.insert( map<MCParticleImpl*, TLorentzVector, less<MCParticleImpl*> >::value_type(a_MCP, p4*_SF));
		m_PDGTGList.insert( map<MCParticleImpl*, double        , less<MCParticleImpl*> >::value_type(a_MCP, PDGTagging(PID)));
	}
	return m_SmearList;
}	

void
MCTruthHelper::printInformation(int print){

	cout << "--------  TRUTH INFORMATION ------------" << endl;
	//cout << "----------------------------------------" << endl;
	//cout << "----------------------------------------" << endl;
	for (unsigned int i = 0; i < m_allMcPList.size(); i++) {
		//MCParticleImpl* mcp = m_finalPList[i]; 
		MCParticleImpl* mcp    =  m_allMcPList[i]; 
		int ID             =  mcp->id();
		int Did[15] = {-1,-1,-1,-1,-1,  -1,-1,-1,-1,-1, -1,-1,-1,-1,-1};
		int Pid[15] = {-1,-1,-1,-1,-1,  -1,-1,-1,-1,-1, -1,-1,-1,-1,-1};
		int pdgid          =  mcp->getPDG();
		int status         =  mcp->getGeneratorStatus(); 
		const double *vecS =  mcp->getVertex(); 
		const double *vecE =  mcp->getEndpoint();
      double d0 = pow(vecS[0]*vecS[0] + vecS[1]*vecS[1], 0.5);
      double z0 = vecS[2];  
		//double vs          =  pow(vecS[0]*vecS[0] + vecS[1]*vecS[1] + vecS[2]*vecS[2], 0.5);
		double ve          =  pow(vecE[0]*vecE[0] + vecE[1]*vecE[1] + vecE[2]*vecE[2], 0.5);
		TLorentzVector  p4 =  TLorentzVector( mcp->getMomentum(), mcp->getEnergy() );
		int nParents       =  (mcp->getParents()).size();
		int nDaughters     =  (mcp->getDaughters()).size();
		int ParentID[20]   =  {0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; 
		int Daughter[20]   =  {0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; 
		if(  nParents  >0) { Pid[ 0]= (mcp->getParents()[ 0])->id();ParentID[ 0]= (mcp->getParents()[ 0])->getPDG();} 
		if(  nParents  >1) { Pid[ 1]= (mcp->getParents()[ 1])->id();ParentID[ 1]= (mcp->getParents()[ 1])->getPDG();} 
		if(  nParents  >2) { Pid[ 2]= (mcp->getParents()[ 2])->id();ParentID[ 2]= (mcp->getParents()[ 2])->getPDG();} 
		if(  nParents  >3) { Pid[ 3]= (mcp->getParents()[ 3])->id();ParentID[ 3]= (mcp->getParents()[ 3])->getPDG();} 
		if(  nParents  >4) { Pid[ 4]= (mcp->getParents()[ 4])->id();ParentID[ 4]= (mcp->getParents()[ 4])->getPDG();} 
		if(  nParents  >5) { Pid[ 5]= (mcp->getParents()[ 5])->id();ParentID[ 5]= (mcp->getParents()[ 5])->getPDG();} 
		if(  nParents  >6) { Pid[ 6]= (mcp->getParents()[ 6])->id();ParentID[ 6]= (mcp->getParents()[ 6])->getPDG();} 
		if(  nParents  >7) { Pid[ 7]= (mcp->getParents()[ 7])->id();ParentID[ 7]= (mcp->getParents()[ 7])->getPDG();} 
		if(  nParents  >8) { Pid[ 8]= (mcp->getParents()[ 8])->id();ParentID[ 8]= (mcp->getParents()[ 8])->getPDG();} 
		if(  nParents  >9) { Pid[ 9]= (mcp->getParents()[ 9])->id();ParentID[ 9]= (mcp->getParents()[ 9])->getPDG();} 
		if(  nParents >10) { Pid[10]= (mcp->getParents()[10])->id();ParentID[10]= (mcp->getParents()[10])->getPDG();} 
		if(  nParents >11) { Pid[11]= (mcp->getParents()[11])->id();ParentID[11]= (mcp->getParents()[11])->getPDG();} 
		if(  nParents >12) { Pid[12]= (mcp->getParents()[12])->id();ParentID[12]= (mcp->getParents()[12])->getPDG();} 
		if(  nParents >13) { Pid[13]= (mcp->getParents()[13])->id();ParentID[13]= (mcp->getParents()[13])->getPDG();} 
		if(  nParents >14) { Pid[14]= (mcp->getParents()[14])->id();ParentID[14]= (mcp->getParents()[14])->getPDG();} 
		if(  nParents >15) { Pid[15]= (mcp->getParents()[15])->id();ParentID[15]= (mcp->getParents()[15])->getPDG();} 
		if(  nParents >16) { Pid[16]= (mcp->getParents()[16])->id();ParentID[16]= (mcp->getParents()[16])->getPDG();} 
		if(  nParents >17) { Pid[17]= (mcp->getParents()[17])->id();ParentID[17]= (mcp->getParents()[17])->getPDG();} 
		if(  nParents >18) { Pid[18]= (mcp->getParents()[18])->id();ParentID[18]= (mcp->getParents()[18])->getPDG();} 
		if(  nParents >19) { Pid[19]= (mcp->getParents()[19])->id();ParentID[19]= (mcp->getParents()[19])->getPDG();} 

		if( nDaughters> 0) { Did[ 0]= (mcp->getDaughters()[ 0])->id();Daughter[ 0]= (mcp->getDaughters()[ 0])->getPDG();} 
		if( nDaughters> 1) { Did[ 1]= (mcp->getDaughters()[ 1])->id();Daughter[ 1]= (mcp->getDaughters()[ 1])->getPDG();} 
		if( nDaughters> 2) { Did[ 2]= (mcp->getDaughters()[ 2])->id();Daughter[ 2]= (mcp->getDaughters()[ 2])->getPDG();} 
		if( nDaughters> 3) { Did[ 3]= (mcp->getDaughters()[ 3])->id();Daughter[ 3]= (mcp->getDaughters()[ 3])->getPDG();} 
		if( nDaughters> 4) { Did[ 4]= (mcp->getDaughters()[ 4])->id();Daughter[ 4]= (mcp->getDaughters()[ 4])->getPDG();} 
		if( nDaughters> 5) { Did[ 5]= (mcp->getDaughters()[ 5])->id();Daughter[ 5]= (mcp->getDaughters()[ 5])->getPDG();} 
		if( nDaughters> 6) { Did[ 6]= (mcp->getDaughters()[ 6])->id();Daughter[ 6]= (mcp->getDaughters()[ 6])->getPDG();} 
		if( nDaughters> 7) { Did[ 7]= (mcp->getDaughters()[ 7])->id();Daughter[ 7]= (mcp->getDaughters()[ 7])->getPDG();} 
		if( nDaughters> 8) { Did[ 8]= (mcp->getDaughters()[ 8])->id();Daughter[ 8]= (mcp->getDaughters()[ 8])->getPDG();} 
		if( nDaughters> 9) { Did[ 9]= (mcp->getDaughters()[ 9])->id();Daughter[ 9]= (mcp->getDaughters()[ 9])->getPDG();} 
		if( nDaughters>10) { Did[10]= (mcp->getDaughters()[10])->id();Daughter[10]= (mcp->getDaughters()[10])->getPDG();} 
		if( nDaughters>11) { Did[11]= (mcp->getDaughters()[11])->id();Daughter[11]= (mcp->getDaughters()[11])->getPDG();} 
		if( nDaughters>12) { Did[12]= (mcp->getDaughters()[12])->id();Daughter[12]= (mcp->getDaughters()[12])->getPDG();} 
		if( nDaughters>13) { Did[13]= (mcp->getDaughters()[13])->id();Daughter[13]= (mcp->getDaughters()[13])->getPDG();} 
		if( nDaughters>14) { Did[14]= (mcp->getDaughters()[14])->id();Daughter[14]= (mcp->getDaughters()[14])->getPDG();} 
		if( nDaughters>15) { Did[15]= (mcp->getDaughters()[15])->id();Daughter[15]= (mcp->getDaughters()[15])->getPDG();} 
		if( nDaughters>16) { Did[16]= (mcp->getDaughters()[16])->id();Daughter[16]= (mcp->getDaughters()[16])->getPDG();} 
		if( nDaughters>17) { Did[17]= (mcp->getDaughters()[17])->id();Daughter[17]= (mcp->getDaughters()[17])->getPDG();} 
		if( nDaughters>18) { Did[18]= (mcp->getDaughters()[18])->id();Daughter[18]= (mcp->getDaughters()[18])->getPDG();} 
		if( nDaughters>19) { Did[19]= (mcp->getDaughters()[19])->id();Daughter[19]= (mcp->getDaughters()[19])->getPDG();} 



		//if (pdgid == kFSRGamma                    ) continue;
		//if (pdgid == kFSRGamma                    ) continue;
		if (pdgid == kGamma  && ParentID[ 0] == kPi0) continue;
		if (pdgid == kEp     && ParentID[ 0] == kPi0) continue;
		if (pdgid == kEm     && ParentID[ 0] == kPi0) continue;
		if (abs(pdgid)== kEp && ParentID[ 0] == kPi0) continue;
		//if( (nParents == 0 && nDaughters == 1 && abs(pdgid) < 20) || ((ParentID==25 || abs(ParentID)==23 || abs(ParentID)==24)&&pdgid<23) || abs(pdgid)==22 ) 
		//printf("trkidx = %3d, nP = %1d (%8d, %5d, %5d, %5d, %5d), nD= %2d, pdgid = %8d, mass=%9.4f, p4 = (%8.3f %8.3f %8.3f %8.2f)\n", 
		//		i, nParents, ParentID0, ParentID1, ParentID2, ParentID3, ParentID4, nDaughters, pdgid, p4.M(), p4.X(), p4.Y(), p4.Z(), p4.T());
		//if( nParents > 1  )
		if( print  < 0 && status ==1 ){
			printf("trkidx = %8d, nP = %2d(%6d), nD = %2d, pdgid = %10s, status = %2d, mass = %7.3f, E = %7.3f, angles = (%5.2f,%5.2f)\n", 
					ID, nParents, Pid[0], nDaughters, getLatexCode(pdgid).c_str(), status, fabs(p4.M()), p4.T(), p4.CosTheta(), p4.Phi());

		} else if(print == 1 )
			printf("trkidx = %8d, nP = %2d(%6d), nD = %2d, pdgid = %10s, status = %2d, mass = %7.3f, E = %7.3f, angles = (%5.2f,%5.2f)\n", 
					ID, nParents, Pid[0], nDaughters, getLatexCode(pdgid).c_str(), status, fabs(p4.M()), p4.T(), p4.CosTheta(), p4.Phi());
		else if ( print == 2 )
			printf("trkidx = %8d, nPs = %2d (%10s,%3d,%3d,%3d,%3d,%3d,%3d,%3d), nDs = %2d, pdgid = %10s, status = %2d, mass = %7.3f, E = %7.3f\n", 
					ID, nParents, getLatexCode(ParentID[0]).c_str(),ParentID[1],ParentID[2],ParentID[3],ParentID[4],ParentID[5],ParentID[6],ParentID[7],nDaughters, getLatexCode(pdgid).c_str(), status, fabs(p4.M()), p4.T());
		else if ( print == 3 ){
                        if( nParents>0 )
			       if ( (mcp->getParents()[0])->getParents().size()>0 )
				       if ( ((mcp->getParents()[0])->getParents()[0])->getPDG() == 92 ) break;

			if( nParents<=1 && nDaughters >=0)
			printf("trkidx = %8d, nPs = %2d (%41s) --> %8s --> %2d daughters, E = %7.3f, angles = (%5.2f,%5.2f)\n", 
					ID, nParents, getLatexCode(ParentID[0]).c_str(), getLatexCode(pdgid).c_str(),
					nDaughters, p4.T(), p4.CosTheta(), p4.Phi());
			if( nParents >1 && nDaughters >=0)
			printf("trkidx = %8d, nPs = %2d (%18s,%2dg,%18s) --> %8s --> %2d daughters, E = %7.3f, angles = (%5.2f,%5.2f)\n", 
					ID, nParents, getLatexCode(ParentID[0]).c_str(), nParents-2,getLatexCode(ParentID[nParents-1]).c_str(), getLatexCode(pdgid).c_str(),
					nDaughters, p4.T(), p4.CosTheta(), p4.Phi());
		}
		else if ( print == 4 )
			printf("trkidx = %8d, nPs = %2d (%10s,%3d,%3d,%3d,%3d,%3d,%3d,%3d), nDs = %2d, pdgid = %10s, status = %2d, mass = %7.3f, Ep= %7.3f\n", 
					ID, nParents, getLatexCode(ParentID[0]).c_str(),ParentID[1],ParentID[2],ParentID[3],ParentID[4],ParentID[5],ParentID[6],ParentID[7],nDaughters, getLatexCode(pdgid).c_str(), status, fabs(p4.M()), ve);
		else if ( print == 5 ){
			if( status ==1 ){ 
				if( nParents==0 )
					printf("trkidx = %8d, nPs = %2d (%41s) --> %8s --> %2d daughters, E = %7.3f, angles = (%5.2f,%5.2f)\n", 
							ID, nParents, getLatexCode(ParentID[0]).c_str(), getLatexCode(pdgid).c_str(),
							nDaughters, p4.T(), p4.CosTheta(), p4.Phi());
				if( nDaughters ==0)
					printf("trkidx = %8d, nPs = %2d (%18s,%2dg,%18s) --> %8s --> %2d daughters, E = %7.3f, angles = (%5.2f,%5.2f)\n", 
							ID, nParents, getLatexCode(ParentID[0]).c_str(), nParents-2,getLatexCode(ParentID[nParents-1]).c_str(), getLatexCode(pdgid).c_str(),
							nDaughters, p4.T(), p4.CosTheta(), p4.Phi());

			}
		}
		else if ( print >  5 ){
			printf("%8d %10d",ID, pdgid);  
			printf(" {");  

			if ( nParents>0){  
				for( int il =0; il<20; il++){
					if (Pid[il]>0)
						if( 0==il ) printf("%d",Pid[il]);
						else printf(",%d",Pid[il]);
					else break; 
				}
			}else printf("0");

			printf("} {");  

			if ( nDaughters>0){  
				for( int il =0; il<20; il++){
					if (Did[il]>0) 
						if( 0==il ) printf("%d",Did[il]);
						else printf(",%d",Did[il]);
					else break; 
				}
			}else printf(" ");

			printf("} ");  
		   printf("%5.1f %5.3f %5.3f",p4.T(), p4.CosTheta(), p4.Phi()/3.1415926535);
		   printf("%10.2le %10.2le\n",d0, z0);
		}
	}
	//cout << "----------------------------------------" << endl;
}

double 
MCTruthHelper::PDGTagging(double PDGTruth)
{
	double TPDG = -100;

	double Chance = gRandom->Rndm(1);

	if(abs(PDGTruth) == 5)
	{
		if(Chance < Pb[0])
		{
			TPDG = 5;	
		}
		else if(Chance >= Pb[0] && Chance < Pb[1] + Pb[0])
		{
			TPDG = 4; 	//C
		}
		else if(Chance >= Pb[1] + Pb[0])
		{
			TPDG = 1; 	//Soft
		}
	}
	if(abs(PDGTruth) == 4)		//c
	{
		if(Chance < Pc[0])
		{
			TPDG = 5; 
		}
		else if(Chance >= Pc[0] && Chance < Pc[1] + Pc[0])
                {
                        TPDG = 4;       //C
                }
                else if(Chance >= Pc[1] + Pc[0])
                {
                        TPDG = 1;       //Soft
                }
	}
	if(abs(PDGTruth) > 0 && (abs(PDGTruth) < 4 || abs(PDGTruth) == 21))	//No top
	{
		if(Chance < Pg[0])
		{
			TPDG = 5;
		}
		else if(Chance >= Pg[0] && Chance < Pg[1] + Pg[0])
		{
			TPDG = 4;       //C
		}
		else if(Chance >= Pg[1] + Pg[0])
		{
			TPDG = 1;       //Soft
		}
	}

	return TPDG;
}

void
MCTruthHelper::initLatexCode(){
	LatexCode.clear();
	LatexCode.insert(map<int,string,less<int> >::value_type(              1, "                 d"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             -1, "            anti-d"));
	LatexCode.insert(map<int,string,less<int> >::value_type(              2, "                 u"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             -2, "            anti-u"));
	LatexCode.insert(map<int,string,less<int> >::value_type(              3, "                 s"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             -3, "            anti-s"));
	LatexCode.insert(map<int,string,less<int> >::value_type(              4, "                 c"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             -4, "            anti-c"));
	LatexCode.insert(map<int,string,less<int> >::value_type(              5, "                 b"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             -5, "            anti-b"));
	LatexCode.insert(map<int,string,less<int> >::value_type(              6, "                 t"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             -6, "            anti-t"));
	LatexCode.insert(map<int,string,less<int> >::value_type(              7, "                b'"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             -7, "           anti-b'"));
	LatexCode.insert(map<int,string,less<int> >::value_type(              8, "                t'"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             -8, "           anti-t'"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             21, "                 g"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             11, "                e-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -11, "                e+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             12, "              nu_e"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -12, "         anti-nu_e"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             13, "               mu-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -13, "               mu+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             14, "             nu_mu"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -14, "        anti-nu_mu"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             15, "              tau-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -15, "              tau+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             16, "            nu_tau"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -16, "       anti-nu_tau"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             17, "                L-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -17, "                L+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             18, "              nu_L"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -18, "         anti-nu_L"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             22, "             gamma"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -22, "          gammaFSR"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10022, "              vpho"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20022, "          Cerenkov"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             23, "                Z0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             24, "                W+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -24, "                W-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             25, "            Higgs0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             28, "           reggeon"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             29, "           pomeron"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             32, "               Z'0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             33, "              Z''0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             34, "               W'+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -34, "               W'-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             35, "           Higgs'0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             36, "                A0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             37, "            Higgs+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -37, "            Higgs-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             40, "                R0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -40, "           anti-R0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             41, "               Xu0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             42, "               Xu+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -42, "               Xu-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             81, "          specflav"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             82, "          rndmflav"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -82, "     anti-rndmflav"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             83, "          phasespa"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             84, "          c-hadron"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -84, "     anti-c-hadron"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             85, "          b-hadron"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -85, "     anti-b-hadron"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             86, "          t-hadron"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -86, "     anti-t-hadron"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             87, "         b'-hadron"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -87, "    anti-b'-hadron"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             88, "         t'-hadron"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -88, "    anti-t'-hadron"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             89, "            Wvirt+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -89, "            Wvirt-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             90, "           diquark"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -90, "      anti-diquark"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             91, "           cluster"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             92, "            string"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             93, "             indep"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             94, "          CMshower"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             95, "          SPHEaxis"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             96, "          THRUaxis"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             97, "           CLUSjet"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             98, "           CELLjet"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             99, "             table"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            111, "               pi0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            211, "               pi+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -211, "               pi-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            210, "          pi_diff+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -210, "          pi_diff-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20111, "           pi(2S)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20211, "           pi(2S)+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -20211, "           pi(2S)-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            221, "               eta"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20221, "           eta(2S)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            331, "              eta'"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            113, "              rho0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            110, "         rho_diff0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            213, "              rho+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -213, "              rho-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          30113, "          rho(2S)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          30213, "          rho(2S)+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -30213, "          rho(2S)-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          40113, "          rho(3S)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          40213, "          rho(3S)+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -40213, "          rho(3S)-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(       -9040213, "        rho(2150)-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9040213, "        rho(2150)+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(       -9040113, "        rho(2150)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            223, "             omega"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            220, "        omega_diff"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          30223, "         omega(2S)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            333, "               phi"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10111, "              a_00"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10211, "              a_0+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10211, "              a_0-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9000221, "          f_0(600)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9010221, "               f_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10221, "              f'_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10113, "              b_10"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10213, "              b_1+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10213, "              b_1-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10223, "               h_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10333, "              h'_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20113, "              a_10"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20213, "              a_1+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -20213, "              a_1-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20223, "               f_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            115, "              a_20"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            215, "              a_2+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -215, "              a_2-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            225, "               f_2"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10221, "         f_0(1370)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          50221, "         f_0(1500)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            335, "              f'_2"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9020221, "         eta(1405)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         100331, "         eta(1475)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10335, "        eta2(1870)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10331, "         f_0(1710)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            229, "         f_4(2050)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20333, "              f'_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        8888888, "         f_0(1790)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9000223, "         f_1(1510)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9050225, "         f_2(1950)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9070225, "         f_2(2150)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9070221, "         f_0(2200)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9020225, "         f_2(1640)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9040225, "         f_2(1910)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9050225, "         f_2(1950)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9080221, "         eta(2225)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9040221, "         eta(1760)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9999999, "           x(1835)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            311, "                K0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -311, "           anti-K0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            310, "              K_S0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            130, "              K_L0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            321, "                K+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -321, "                K-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            313, "               K*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -313, "          anti-K*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            323, "               K*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -323, "               K*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9000311, "             kapa0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(       -9000311, "        anti-kapa0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9000321, "             kapa+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(       -9000321, "             kapa-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10311, "             K_0*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10311, "        anti-K_0*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10321, "             K_0*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10321, "             K_0*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10313, "              K_10"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10313, "         anti-K_10"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10323, "              K_1+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10323, "              K_1-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            315, "             K_2*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -315, "        anti-K_2*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            325, "             K_2*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -325, "             K_2*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20313, "             K'_10"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -20313, "        anti-K'_10"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20323, "             K'_1+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -20323, "             K'_1-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         100313, "              K'*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        -100313, "         anti-K'*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         100323, "              K'*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        -100323, "              K'*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          30313, "             K''*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -30313, "        anti-K''*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          30323, "             K''*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -30323, "             K''*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10315, "              K_20"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10315, "         anti-K_20"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10325, "              K_2+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10325, "              K_2-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            317, "             K_3*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -317, "        anti-K_3*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            327, "             K_3*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -327, "             K_3*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            319, "             K_4*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -319, "        anti-K_4*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            329, "             K_4*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -329, "             K_4*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9060221, "         f_0(2100)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            411, "                D+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -411, "                D-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            421, "                D0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -421, "           anti-D0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            422, "               D0H"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -422, "               D0L"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            413, "               D*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -413, "               D*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            423, "               D*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -423, "          anti-D*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10411, "             D_0*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10411, "             D_0*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10421, "             D_0*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10421, "        anti-D_0*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10413, "              D_1+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10413, "              D_1-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10423, "              D_10"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10423, "         anti-D_10"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            415, "             D_2*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -415, "             D_2*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            425, "             D_2*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -425, "        anti-D_2*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20413, "             D'_1+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -20413, "             D'_1-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20423, "             D'_10"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -20423, "        anti-D'_10"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            431, "              D_s+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -431, "              D_s-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            433, "             D_s*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -433, "             D_s*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10431, "            D_s0*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10431, "            D_s0*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10433, "             D_s1+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10433, "             D_s1-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            435, "            D_s2*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -435, "            D_s2*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20433, "            D'_s1+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -20433, "            D'_s1-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          30411, "            D(2S)+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -30411, "            D(2S)-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          30421, "            D(2S)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -30421, "       anti-D(2S)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          30413, "           D*(2S)+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -30413, "           D*(2S)-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          30423, "           D*(2S)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -30423, "      anti-D*(2S)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            511, "                B0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -511, "           anti-B0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            150, "               B0L"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            510, "               B0H"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            521, "                B+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -521, "                B-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            513, "               B*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -513, "          anti-B*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            523, "               B*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -523, "               B*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10511, "             B_0*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10511, "        anti-B_0*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10521, "             B_0*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10521, "             B_0*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10513, "              B_10"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10513, "         anti-B_10"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10523, "              B_1+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10523, "              B_1-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            515, "             B_2*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -515, "        anti-B_2*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            525, "             B_2*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -525, "             B_2*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20513, "             B'_10"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -20513, "        anti-B'_10"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20523, "             B'_1+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -20523, "             B'_1-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            531, "              B_s0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -531, "         anti-B_s0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            350, "             B_s0L"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            530, "             B_s0H"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            533, "             B_s*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -533, "        anti-B_s*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10531, "            B_s0*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10531, "       anti-B_s0*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10533, "             B_s10"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10533, "        anti-B_s10"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            535, "            B_s2*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -535, "       anti-B_s2*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20533, "            B'_s10"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -20533, "       anti-B'_s10"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            541, "              B_c+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -541, "              B_c-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            543, "             B_c*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -543, "             B_c*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10541, "            B_c0*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10541, "            B_c0*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10543, "             B_c1+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -10543, "             B_c1-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            545, "            B_c2*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           -545, "            B_c2*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20543, "            B'_c1+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -20543, "            B'_c1-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            441, "             eta_c"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20441, "         eta_c(2S)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            443, "             J/psi"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         100443, "           psi(2S)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          30443, "         psi(3770)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9000443, "         psi(4040)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9010443, "         psi(4160)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9030443, "         psi(4260)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9040443, "         psi(4360)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(        9020443, "         psi(4415)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10443, "               h_c"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10441, "            chi_c0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20443, "            chi_c1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            445, "            chi_c2"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            551, "             eta_b"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20551, "         eta_b(2S)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          40551, "         eta_b(3S)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            553, "           Upsilon"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          30553, "       Upsilon(2S)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          60553, "       Upsilon(3S)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          70553, "       Upsilon(4S)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          80553, "       Upsilon(5S)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10553, "               h_b"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          40553, "           h_b(2P)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         100553, "           h_b(3P)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10551, "            chi_b0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20553, "            chi_b1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            555, "            chi_b2"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          30551, "        chi_b0(2P)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          50553, "        chi_b1(2P)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10555, "        chi_b2(2P)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          50551, "        chi_b0(3P)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         110553, "        chi_b1(3P)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          20555, "        chi_b2(3P)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          40555, "        eta_b2(1D)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          60555, "        eta_b2(2D)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         120553, "     Upsilon_1(1D)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          30555, "     Upsilon_2(1D)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            557, "     Upsilon_3(1D)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         130553, "     Upsilon_1(2D)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          50555, "     Upsilon_2(2D)"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          10557, "     Upsilon_3(2D)"));
	//
	LatexCode.insert(map<int,string,less<int> >::value_type(          10222, "           sigma_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           1114, "            Delta-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -1114, "       anti-Delta+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           2110, "           n_diffr"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -2110, "      anti-n_diffr"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           2112, "                n0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -2112, "           anti-n0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           2114, "            Delta0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -2114, "       anti-Delta0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           2210, "           p_diff+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -2210, "      anti-p_diff-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           2212, "                p+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -2212, "           anti-p-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           2214, "            Delta+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -2214, "       anti-Delta-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           2224, "           Delta++"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -2224, "      anti-Delta--"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3112, "            Sigma-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3112, "       anti-Sigma+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3114, "           Sigma*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3114, "      anti-Sigma*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3122, "           Lambda0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3122, "      anti-Lambda0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          13122, "     Lambda(1405)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -13122, "anti-Lambda(1405)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3124, "     Lambda(1520)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3124, "anti-Lambda(1520)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          23122, "     Lambda(1600)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -23122, "anti-Lambda(1600)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          33122, "     Lambda(1670)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -33122, "anti-Lambda(1670)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          13124, "     Lambda(1690)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -13124, "anti-Lambda(1690)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          43122, "     Lambda(1800)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -43122, "anti-Lambda(1800)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          53122, "     Lambda(1810)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -53122, "anti-Lambda(1810)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3126, "     Lambda(1820)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3126, "anti-Lambda(1820)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          13126, "     Lambda(1830)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -13126, "anti-Lambda(1830)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          13212, "      Sigma(1660)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -13212, " anti-Sigma(1660)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          13214, "      Sigma(1670)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -13214, " anti-Sigma(1670)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          23212, "      Sigma(1750)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -23212, " anti-Sigma(1750)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3216, "      Sigma(1775)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3216, " anti-Sigma(1775)0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3212, "            Sigma0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3212, "       anti-Sigma0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3214, "           Sigma*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3214, "      anti-Sigma*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3222, "            Sigma+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3222, "       anti-Sigma-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3224, "           Sigma*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3224, "      anti-Sigma*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3312, "               Xi-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3312, "          anti-Xi+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3314, "              Xi*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3314, "         anti-Xi*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3322, "               Xi0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3322, "          anti-Xi0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3324, "              Xi*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3324, "         anti-Xi*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3334, "            Omega-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3334, "       anti-Omega+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          14122, "        Lambda_c*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -14122, "   anti-Lambda_c*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          14124, "        Lambda_c*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -14124, "   anti-Lambda_c*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4112, "          Sigma_c0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4112, "     anti-Sigma_c0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4114, "         Sigma_c*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4114, "    anti-Sigma_c*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4212, "          Sigma_c+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4212, "     anti-Sigma_c-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4214, "         Sigma_c*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4214, "    anti-Sigma_c*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4222, "         Sigma_c++"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4222, "    anti-Sigma_c--"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4224, "        Sigma_c*++"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4224, "   anti-Sigma_c*--"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4312, "            Xi'_c0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4312, "       anti-Xi'_c0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4322, "            Xi'_c+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4322, "       anti-Xi'_c-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4324, "            Xi_c*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4324, "       anti-Xi_c*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4122, "         Lambda_c+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4122, "    anti-Lambda_c-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4132, "             Xi_c0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4132, "        anti-Xi_c0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4232, "             Xi_c+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4232, "        anti-Xi_c-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4314, "            Xi_c*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4314, "       anti-Xi_c*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4332, "          Omega_c0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4332, "     anti-Omega_c0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4334, "         Omega_c*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4334, "    anti-Omega_c*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5112, "          Sigma_b-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5112, "     anti-Sigma_b+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5114, "         Sigma_b*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5114, "    anti-Sigma_b*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5122, "         Lambda_b0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5122, "    anti-Lambda_b0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5132, "             Xi_b-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5132, "        anti-Xi_b+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5212, "          Sigma_b0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5212, "     anti-Sigma_b0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5214, "         Sigma_b*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5214, "    anti-Sigma_b*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5222, "          Sigma_b+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5222, "     anti-Sigma_b-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5224, "         Sigma_b*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5224, "    anti-Sigma_b*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5232, "             Xi_b0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5232, "        anti-Xi_b0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5312, "            Xi'_b-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5312, "       anti-Xi'_b+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5314, "            Xi_b*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5314, "       anti-Xi_b*+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5322, "            Xi'_b0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5322, "       anti-Xi'_b0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5324, "            Xi_b*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5324, "       anti-Xi_b*0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5332, "          Omega_b-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5332, "     anti-Omega_b+"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5334, "         Omega_b*-"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5334, "    anti-Omega_b*+"));
	/*
	LatexCode.insert(map<int,string,less<int> >::value_type(           1101, "              dd_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -1101, "         anti-dd_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           2101, "              ud_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -2101, "         anti-ud_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           2201, "              uu_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -2201, "         anti-uu_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3101, "              sd_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3101, "         anti-sd_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3201, "              su_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3201, "         anti-su_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3301, "              ss_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3301, "         anti-ss_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4101, "              cd_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4101, "         anti-cd_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4201, "              cu_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4201, "         anti-cu_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4301, "              cs_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4301, "         anti-cs_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4401, "              cc_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4401, "         anti-cc_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5101, "              bd_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5101, "         anti-bd_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5201, "              bu_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5201, "         anti-bu_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5301, "              bs_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5301, "         anti-bs_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5401, "              bc_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5401, "         anti-bc_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5501, "              bb_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5501, "         anti-bb_0"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           1103, "              dd_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -1103, "         anti-dd_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           2103, "              ud_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -2103, "         anti-ud_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           2203, "              uu_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -2203, "         anti-uu_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3103, "              sd_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3103, "         anti-sd_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3203, "              su_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3203, "         anti-su_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           3303, "              ss_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -3303, "         anti-ss_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4103, "              cd_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4103, "         anti-cd_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4203, "              cu_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4203, "         anti-cu_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4303, "              cs_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4303, "         anti-cs_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           4403, "              cc_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -4403, "         anti-cc_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5103, "              bd_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5103, "         anti-bd_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5203, "              bu_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5203, "         anti-bu_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5303, "              bs_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5303, "         anti-bs_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5403, "              bc_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5403, "         anti-bc_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           5503, "              bb_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -5503, "         anti-bb_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           1011, "          deuteron"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -1011, "     anti-deuteron"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           1021, "           tritium"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -1021, "      anti-tritium"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           1012, "               He3"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -1012, "          anti-He3"));
	LatexCode.insert(map<int,string,less<int> >::value_type(           1022, "             alpha"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          -1022, "        anti-alpha"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            100, "          geantino"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            101, "   chargedgeantino"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          30343, "               Xsd"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -30343, "          anti-Xsd"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          30353, "               Xsu"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -30353, "          anti-Xsu"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          30373, "               Xdd"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -30373, "          anti-Xdd"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          30383, "               Xdu"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -30383, "          anti-Xdu"));
	LatexCode.insert(map<int,string,less<int> >::value_type(          30363, "               Xss"));
	LatexCode.insert(map<int,string,less<int> >::value_type(         -30363, "          anti-Xss"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             51, "         dummy00_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             52, "         dummy10_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             53, "         dummy01_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             54, "         dummy11_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -51, "    anti-dummy00_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -52, "    anti-dummy10_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -53, "    anti-dummy01_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -54, "    anti-dummy11_1"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             55, "         dummy00_2"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             56, "         dummy10_2"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             57, "         dummy01_2"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             58, "         dummy11_2"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -55, "    anti-dummy00_2"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -56, "    anti-dummy10_2"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -57, "    anti-dummy01_2"));
	LatexCode.insert(map<int,string,less<int> >::value_type(            -58, "    anti-dummy11_2"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             61, "           chi_c0p"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             62, "           chi_c1p"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             63, "           chi_c2p"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             64, "             X3940"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             65, "             Y3940"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             66, "             xvpho"));
	LatexCode.insert(map<int,string,less<int> >::value_type(             67, "            eta_c2"));
	*/
}

string
MCTruthHelper::getLatexCode(int pdgid){
	string code="           Unknown";
	map< int, string>::iterator pos = LatexCode.find(pdgid); 
	if (pos != LatexCode.end()) 
	{
		code=pos->second;
	}
	return code;
}
bool MCTruthHelper::FindHiggsMother(MCParticle *mcp)
{
	bool ret = false;

	int PDGIDItself = abs(mcp->getPDG()); 
	if ( PDGIDItself == 25 )      // original particle, no parent  
	{
		ret = true;
	}
	else if ( mcp->getParents().size() == 0 )      // original particle, no parent  
	{
		ret= false; 
	} 
	else if ( mcp->getParents().size() == 1 ) // Only one parent  
	{
		MCParticle *mcpa = mcp->getParents()[0];
		int PDGIDParent  = abs(mcpa->getPDG()); 
		if ( PDGIDParent == 25){
           ret = true; 
		}else{
			ret = FindHiggsMother(mcpa);	
		}
	}else{
		int np = mcp->getParents().size(); 
      for( int i =0; i< np;i++){
			MCParticle *mcpa = mcp->getParents()[i];
			ret = FindHiggsMother(mcpa);
			if( ret )	break;
		}
	}
	//
	return ret;
}


