#ifndef LEAKAGE_H
#define LEAKAGE_H
#include <string>
#include <vector>
#include <IMPL/MCParticleImpl.h>
using namespace std;
using namespace IMPL;
using namespace EVENT;

namespace Utility
{
	template <typename T> void FreeDelAll( T & t ) 
	{
		for(unsigned int i=0; i<t.size(); i++) delete t[i];
		T tmp; 	t.swap(tmp);
	}

	template <typename T> void FreeAll( T & t ) 
	{
		T tmp; 	t.swap(tmp);
	}

	double Mass(const string name);

	int    PdgCode(const string name);

	vector<MCParticle*> getEffectiveParent( MCParticle* );
	vector<MCParticle*> getDaughterParent ( MCParticle* );
}
#endif
