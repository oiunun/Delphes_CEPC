#include "SimpleMCEvent.h"


SimpleMCEvent::SimpleMCEvent( LCCollection* mcParticleCol ){
	
	int _nMCP= mcParticleCol->getNumberOfElements();
	for(int i1 = 0; i1 < _nMCP; i1++)
	{
		MCParticleImpl *mcp = dynamic_cast<MCParticleImpl *>(mcParticleCol->getElementAt(i1));
		//
		int status      =  (mcp)->getGeneratorStatus(); 
		int pdgid       =  (mcp)->getPDG();
		int nParents    = ((mcp)->getParents()).size();
		int nDaughters  = ((mcp)->getDaughters()).size();
		if ( nParents   == 0 && abs(pdgid)==21) continue;
		if ( nDaughters == 0 && status == 2   ) continue;
		if ( status != 1     && status != 2   ) continue;
		//
	}

}

SimpleMCEvent::~SimpleMCEvent(){

}

