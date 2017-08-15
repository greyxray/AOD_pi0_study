#include "AOD_pi0.h"

void AOD_pi0::ResetBranchesPionTree()
{
	v0_pion.SetXYZM(-1., -1., -1.,-1);
	v0_pion_trackRef.SetXYZM(-1., -1., -1.,-1);
	tau_pion.SetXYZM(-1., -1., -1.,-1);
	tau_pion_besttrack.SetXYZM(-1., -1., -1.,-1);
	// Diff between HPS and HPS best track
	hpsbestTrackDiff_fx = -1.;
	hpsbestTrackDiff_fy = -1.;
	hpsbestTrackDiff_fz = -1.;
	hpsbestTrackDiff_fE = -1.;
	hpsbestTrack_dR = -1.;
	// Diff between V0 and HPS
	hpsv0Diff_fx = -1.;
	hpsv0Diff_fy = -1.;
	hpsv0Diff_fz = -1.;
	hpsv0Diff_fE = -1.;
	hpsv0Diff_dPhi = -1;
	hpsv0Diff_dEta = -1;
	hpsv0_dR = -1.;
	// Diff between V0 and HPS best track
	v0bestTrackDiff_fx = -1.;
	v0bestTrackDiff_fy = -1.;
	v0bestTrackDiff_fz = -1.;
	v0bestTrackDiff_fE = -1.;
	v0bestTrack_dR = -1.;
	// associated mass
	bestTrac_fM = -1.;
	hps_fM = -1.;
	v0_fM = -1.;
}

// ,reco::candidate
template  <typename HPSPion, typename Daughter>
void AOD_pi0::FillPion(HPSPion pion, Daughter daughter)
{
	ResetBranchesPionTree();
	v0_pion.SetPxPyPzE(daughter->px(),daughter->py(),daughter->pz(),daughter->energy());
	v0_pion_trackRef.SetXYZM( ((reco::RecoChargedCandidate *)daughter)->track()->px(), ((reco::RecoChargedCandidate *)daughter)->track()->py(), ((reco::RecoChargedCandidate *)daughter)->track()->pz(), 0.13957018);
	tau_pion.SetPxPyPzE(pion->px(), pion->py(), pion->pz(), pion->energy());
	tau_pion_besttrack.SetXYZM(pion->bestTrack()->px(), pion->bestTrack()->py(), pion->bestTrack()->pz(), 0.13957018);

	dlog("\t\t", "v0_pion_trackRef:", v0_pion_trackRef.X(), v0_pion_trackRef.Y(), v0_pion_trackRef.Z(), v0_pion_trackRef.M(),  v0_pion_trackRef.E());
	dlog("\t\t", "tau_ptau_pion_besttrackion:", tau_pion_besttrack.X(), tau_pion_besttrack.Y(), tau_pion_besttrack.Z(), tau_pion_besttrack.M(), tau_pion_besttrack.E());
	dlog("\t\t", "v0_pion:", v0_pion.X(), v0_pion.Y(), v0_pion.Z(), v0_pion.M(),  v0_pion.E());
	dlog("\t\t", "tau_pion:", tau_pion.X(), tau_pion.Y(), tau_pion.Z(), tau_pion.M(), tau_pion.E());
	dlog();

	hpsbestTrackDiff_fx = tau_pion.X() - tau_pion_besttrack.X();
	hpsbestTrackDiff_fy = tau_pion.Y() - tau_pion_besttrack.Y();
	hpsbestTrackDiff_fz = tau_pion.Z() - tau_pion_besttrack.Z();
	hpsbestTrackDiff_fE = tau_pion.E() - tau_pion_besttrack.E();
	hpsbestTrack_dR = tau_pion.DeltaR(tau_pion_besttrack);

	hpsv0Diff_fx = tau_pion.X() - v0_pion.X();
	hpsv0Diff_fy = tau_pion.Y() - v0_pion.Y();
	hpsv0Diff_fz = tau_pion.Z() - v0_pion.Z();
	hpsv0Diff_fE = tau_pion.E() - v0_pion.E();
	hpsv0Diff_dPhi = tau_pion.DeltaPhi(v0_pion);
	hpsv0Diff_dEta = tau_pion.Eta() - v0_pion.Eta();
	hpsv0_dR = tau_pion.DeltaR(v0_pion);

	v0bestTrackDiff_fx = v0_pion.X() - tau_pion_besttrack.X();
	v0bestTrackDiff_fy = v0_pion.Y() - tau_pion_besttrack.Y();
	v0bestTrackDiff_fz = v0_pion.Z() - tau_pion_besttrack.Z();
	v0bestTrackDiff_fE = v0_pion.E() - tau_pion_besttrack.E();
	v0bestTrack_dR = v0_pion.DeltaR(tau_pion_besttrack);

	bestTrac_fM = tau_pion_besttrack.M();
	hps_fM = tau_pion.M();
	v0_fM = v0_pion.M();

	pions_tree->Fill();	
}