#include "SimpleMCParticle.h"


SimpleMCParticle::SimpleMCParticle( MCParticle* mcp, vector<MCParticle*>& parents, vector<MCParticle*>& daughters  ){
	
	 _status      =  mcp->getGeneratorStatus(); 
	 _id          =  mcp->id();
	 _pdgid       =  mcp->getPDG();
	 _mass        =  mcp->getMass();
	 _charge      =  mcp->getCharge();
	 _nparents    =  parents.size();
	 _ndaughters  =  daughters.size();
	 _p4          =  TLorentzVector( mcp->getMomentum(), mcp->getEnergy() );


		if(  _nparents   >0) {
			for( int i=0; i<_nparents; i++){
			 //_parents.push_back( new SimpleMCParticle( parents[i] ) ); 
			}
		}
		if(  _ndaughters >0){
			for( int i=0; i<_ndaughters; i++){
			 //_daughters.push_back( new SimpleMCParticle( daughters[i] ) ); 
			}
		}
}

SimpleMCParticle::~SimpleMCParticle(){

	FreeDelAll(_parents);
	FreeDelAll(_daughters);

}

