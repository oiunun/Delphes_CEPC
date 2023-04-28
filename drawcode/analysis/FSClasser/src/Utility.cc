#include "Utility.h"
//
double Utility::Mass(const string name){
	if      ( name == "pi0"    )                       return 0.1349766;
	else if ((name == "pi-"    )||(name == "pi+"    )) return 0.1395702;
	else if ( name == "Ks"     )                       return 0.4976140;
	else if ((name == "K-"     )||(name == "K+"     )) return 0.4936770;
	else if ( name == "gamma"  )                       return 0.0000000;
	else if ( name == "eta"    )                       return 0.5478530;
	else if ((name == "p-"     )||(name == "p+"     )) return 0.9382720;
	else if ((name == "n0"     )||(name == "anti-n0")) return 0.9395653;
	else if ((name == "tau-"   )||(name == "tau+"   )) return 1.7768200;
	else if ((name == "mu-"    )||(name == "mu+"    )) return 0.1056584;
	else if ((name == "e-"     )||(name == "e+"     )) return 0.0005110;
	else if ((name == "ALambda")||(name == "Lambda" )) return 1.1156830;
	else if ((name == "D0"     )||(name == "D0bar"  )) return 1.8648300;
	else if ((name == "D+"     )||(name == "D-"     )) return 1.8696000;
	else if ((name == "Ds+"    )||(name == "Ds-"    )) return 1.9684700;
	else if ( name == "J/psi"  )                       return 3.0969160;
	else if ( name == "h_c"    )                       return 3.5254200;
	else if ( name == "chi_c1" )                       return 3.5106600;
	else if ( name == "psi(2S)")                       return 3.6860900;
	else if ( name == "nu"     )                       return 0.0000000;
	else if ( name == "Z0"     )                       return 91.187600;
	else if ((name == "W+"     )||(name == "W-"     )) return 80.399000;
	return 0.0;
}

int Utility::PdgCode(const string name){

	if      ( name == "e+"    )                       return         -11;
	else if ( name == "e-"    )                       return          11;
	else if ( name == "mu+"   )                       return         -13;
	else if ( name == "mu-"   )                       return          13;
	else if ( name == "tau+"  )                       return         -15;
	else if ( name == "tau-"  )                       return          15;
	else if ( name == "gamma" )                       return          22;
	else if ( name == "jet"   )                       return          5 ;

	return 0;
}

/*
vector<MCParticle*> Utility::getEffectiveParents(MCParticle * mcp){

}

vector<MCParticle*> Utility::getEffectiveDaughters(MCParticle * mcp){

}
*/
