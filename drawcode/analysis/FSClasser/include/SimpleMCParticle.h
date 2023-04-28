#ifndef SimpleMCPARTILCLE_H
#define SimpleMCPARTILCLE_H 1

#include <string>
#include <map>
#include <set>
#include <TLorentzVector.h>
#include <IMPL/MCParticleImpl.h>

#include "Utility.h"

using namespace std; 
using namespace IMPL;
using namespace EVENT;
using namespace Utility;

class SimpleMCParticle  {

	public: 
		SimpleMCParticle( MCParticle* mcp, vector<MCParticle*>& parents, vector<MCParticle*>& daughters );
		~SimpleMCParticle(); 
		int                        getId       () { return _id       ;} 
		int                        getPdgid    () { return _pdgid    ;}
		double                     getCharge   () { return _charge   ;}
		double                     getMass     () { return _mass     ;}
		TLorentzVector             getP4       () { return _p4       ;}
		string                     getName     () { return _name     ;}
		vector<SimpleMCParticle*>  getParents  () { return _parents  ;}
		vector<SimpleMCParticle*>  getDaughters() { return _daughters;}
		
		void                       setId       (int                       id       ) { _id       = id       ;} 
		void                       setPdgid    (int                       pdgid    ) { _pdgid    = pdgid    ;}
		void                       setCharge   (double                    charge   ) { _charge   = charge   ;}
		void                       setMass     (double                    mass     ) { _mass     = mass     ;}
		void                       setP4       (TLorentzVector            p4       ) { _p4       = p4       ;}
		void                       setName     (string                    name     ) { _name     = name     ;}
		void                       setParents  (vector<SimpleMCParticle*> parents  ) { _parents  = parents  ;}
		void                       setDaughters(vector<SimpleMCParticle*> daughters) { _daughters= daughters;}

	protected:

		int                        _id       ;
		int                        _pdgid    ;
		int                        _nparents, _ndaughters, _status;
		double                     _charge   ;
		double                     _mass     ;
		TLorentzVector             _p4       ;
		string                     _name     ;
		
		vector<SimpleMCParticle*>  _parents  ;
		vector<SimpleMCParticle*>  _daughters;

}; 


#endif /* ifndef SimpleMCPARTILCLE__H */
