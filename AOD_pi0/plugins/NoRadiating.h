#include "AOD_pi0.h"

bool AOD_pi0::NoRadiating(const reco::Candidate * mom, int daug_pdgid)//map_kaons
{
	for(unsigned i = 0; i < mom->numberOfDaughters(); i++)
	{
		const reco::Candidate * prt = mom->daughter(i);
		if (abs(prt->pdgId()) == daug_pdgid) return false;
	}
	return true;
}

bool AOD_pi0::NoRadiating(const reco::Candidate * mom, map <long, string> &map_daug_pdgid)//map_kaons
{
	dout("mom:",mom->pdgId(), "has", mom->numberOfDaughters(), "daughters");
	for(unsigned i = 0; i < mom->numberOfDaughters(); i++)
	{
		const reco::Candidate * prt = mom->daughter(i);
		dout("\tdaughter", i, "is of type", abs(prt->pdgId()));
		if (map_daug_pdgid.find(abs(prt->pdgId())) != map_daug_pdgid.end())
		{
			dout("\t\tfound that daughter", i, "of type", abs(prt->pdgId()), "is in map", map_daug_pdgid.find(abs(prt->pdgId()))->second );
			return false;
		}
	}
	return true;
}