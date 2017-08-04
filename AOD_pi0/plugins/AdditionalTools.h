#include "AOD_pi0.h"

template <typename T>
vector <T*> AOD_pi0::TransformToPointers(vector <T> a, vector <T*> b)
{
	// if b is empty (this will append to the end of b)
	if (a.size() > 0)
	{
		b.reserve(a.size()); // optional, but a good habit
		std::transform(a.begin(), a.end(), std::back_inserter(b), [](T& o){ return &o; });
	}
	return b;
}

template <typename Tc, typename Tr>// T - collection   reco::PFTauRef pftauref(PF_hps_taus, i);
void AOD_pi0::MakeVectorofRef(edm::Handle< Tc > Collection, vector< Tr* > v_of_ref) {}


template <typename Ta, typename Tb>
bool AOD_pi0::Incone(Ta& A, Tb& B, double cone_size)
{
	TLorentzVector a(A.px(), A.py(), A.pz(), A.energy());
	TLorentzVector b(B.px(), B.py(), B.pz(), B.energy());
	double dR = a.DeltaR(b);
	if (dR <= cone_size) return true;
	else return false;
}
