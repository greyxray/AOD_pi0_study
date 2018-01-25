#include "AOD_pi0.h"

template<class T, class P, class R>
float AOD_pi0::getDxy(const T pv, const P p4, const R ref) const
{
	return (
		- (ref.x() - pv.X()) * p4.Y()
		+ (ref.y() - pv.Y()) * p4.X()
	) / sqrtf(p4.Perp2());
}

template<class T, class P, class R>
float AOD_pi0::getDz(const T pv, const P p4, const R ref) const
{
	return ref.z() - pv.Z() - (
			(ref.x() - pv.X()) * p4.X() +
			(ref.y() - pv.Y()) * p4.Y()
		) * p4.Z() / p4.Perp2();
}

void AOD_pi0::K892(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	iEvent.getByToken( TauPiZeroCollectionToken_, Strips);

		///K(892)+- TODO: build using V0 and hps pi0; made: patching to pi+- of hps and V0
	if (false && PF_hps_taus.isValid() )
	{
		vector<reco::CandidateCollection> v0_daughters;
		v0_Ks_count = V0Ks->size();
		// Loop over PF::Taus
		for( unsigned tau_index = 0; tau_index < PF_hps_taus->size(); tau_index++)
		{
			reco::PFTauRef pftauref(PF_hps_taus, tau_index);
			//pi+-
			vector < reco::PFCandidatePtr  > tau_signalPFChargedHadrCands    = pftauref->signalPFChargedHadrCands();//typedef edm::Ptr<PFCandidate> reco::PFCandidatePtr
			vector < reco::PFCandidatePtr  > tau_isolationPFChargedHadrCands = pftauref->isolationPFChargedHadrCands(); // //typedef edm::Ptr<PFCandidate> reco::PFCandidatePtr
			// All pi+-'s
			vector < reco::PFCandidatePtr > tau_picharge = pftauref->signalPFChargedHadrCands(); // vector < edm::Ptr<PFCandidate> >
			tau_picharge.insert(tau_picharge.end(), tau_isolationPFChargedHadrCands.begin(), tau_isolationPFChargedHadrCands.end());

			//dout("tau:", tau_index);
			// Loop over pi+- in Tau
			double E = 0, p_x = 0, p_y = 0, p_z = 0, inv_M = 0;
			for(unsigned pi_i = 0 ; pi_i < tau_picharge.size() ; pi_i++)
			{
				TLorentzVector tau_pion(tau_picharge[pi_i]->px(), tau_picharge[pi_i]->py(), tau_picharge[pi_i]->pz(), tau_picharge[pi_i]->energy());

				// Loop over known Ks
				for(unsigned k0_i = 0 ; k0_i < V0Ks->size() ; k0_i++)// ISSUE: not all K0s's pi+- are matched in one event to the pions of tau
				{
					TLorentzVector v0_ks_pion1((*V0Ks)[k0_i].daughter(0)->px(),(*V0Ks)[k0_i].daughter(0)->py(),(*V0Ks)[k0_i].daughter(0)->pz(),(*V0Ks)[k0_i].daughter(0)->energy());
					TLorentzVector v0_ks_pion2((*V0Ks)[k0_i].daughter(1)->px(),(*V0Ks)[k0_i].daughter(1)->py(),(*V0Ks)[k0_i].daughter(1)->pz(),(*V0Ks)[k0_i].daughter(1)->energy());

					TLorentzVector v0_ks((*V0Ks)[k0_i].px(),(*V0Ks)[k0_i].py(),(*V0Ks)[k0_i].pz(),(*V0Ks)[k0_i].energy());
					double k892_m_inv = (v0_ks + tau_pion).M(); //K892 from the deffinitely K0s from V0 and a pi+- inside the tau
					if (tau_pion.M() < 0.1) break;
					// {
					//   dout("pion mass can not be <135", tau_pion.M(), "is it?", tau_picharge[pi_i]->mass());
					//   exit(1);
					// }
					if (pow(tau_picharge[pi_i]->charge() - (*V0Ks)[k0_i].daughter(0)->charge(), 2) +
										pow(tau_picharge[pi_i]->px() - (*V0Ks)[k0_i].daughter(0)->px(), 2) +
										pow(tau_picharge[pi_i]->py() - (*V0Ks)[k0_i].daughter(0)->py(), 2) +
										pow(tau_picharge[pi_i]->pz() - (*V0Ks)[k0_i].daughter(0)->pz(), 2) +
										pow(tau_picharge[pi_i]->energy() - (*V0Ks)[k0_i].daughter(0)->energy(), 2) < 0.001 ||
									pow(tau_picharge[pi_i]->charge() - (*V0Ks)[k0_i].daughter(1)->charge(), 2) +
										pow(tau_picharge[pi_i]->px() - (*V0Ks)[k0_i].daughter(1)->px(), 2) +
										pow(tau_picharge[pi_i]->py() - (*V0Ks)[k0_i].daughter(1)->py(), 2) +
										pow(tau_picharge[pi_i]->pz() - (*V0Ks)[k0_i].daughter(1)->pz(), 2) +
										pow(tau_picharge[pi_i]->energy() - (*V0Ks)[k0_i].daughter(1)->energy(), 2) < 0.001
										)
					 {
						//HAVE TO PASS THIS PION SINCE IT ALREADY INCLUDED
						//dout("\t\tis the same as pion in Ks:", k0_i, ":");
						//dout("\t\tfirst:", (*V0Ks)[k0_i].daughter(0)->charge(), (*V0Ks)[k0_i].daughter(0)->px(), (*V0Ks)[k0_i].daughter(0)->py(), (*V0Ks)[k0_i].daughter(0)->pz(), (*V0Ks)[k0_i].daughter(0)->energy());
						//dout("\t\tsecond:", (*V0Ks)[k0_i].daughter(1)->charge(), (*V0Ks)[k0_i].daughter(1)->px(), (*V0Ks)[k0_i].daughter(1)->py(), (*V0Ks)[k0_i].daughter(1)->pz(), (*V0Ks)[k0_i].daughter(1)->energy());
					 }
					 else if (//abs(k892_m_inv - 0.892) < 0.02 &&
											(abs(tau_picharge[pi_i]->vz() ) < 0.2))
						{
							dout("FOUND k(892)");
							//TEMP_COMMENTED h_K892->Fill(k892_m_inv);
						}
						else
						{
							//dout(k892_m_inv, "while kaon", v0_ks.M(), "is", (*V0Ks)[k0_i].mass());
						}
				}
				// if (v0ks_piCharged->pt() > max_ks_daughter_pt) max_ks_daughter_pt = v0ks_piCharged->pt();
				// ks_daughter_pt->Fill(v0ks_piCharged->pt());

				// E += v0ks_piCharged->energy();
				// p_x += v0ks_piCharged->px();
				// p_y += v0ks_piCharged->py();
				// p_z += v0ks_piCharged->pz();
			}

			inv_M += sqrt(pow(E, 2) - pow(p_x, 2) - pow(p_y, 2) - pow(p_z, 2));
			// ks_inv_m_pi->Fill(inv_M);
		}
	}
	else if (!PF_hps_taus.isValid()) dout("no valid PF_hps_taus");

	/*
	1) access the standart collection of V0
	2) access the Strips collection
	3) save to root files the two collections separately the Strip and Kaon whose invariant mass is in the range of of K890 (? and oposite direction?)
	*/
	if (V0Ks_standart.isValid() && Strips.isValid() && Strips->size() > 0 && V0Ks_standart->size() > 0)
	{
		int nK892_temp = 0;
		n_K892 = -1;
		for( unsigned v0_index = 0; v0_index < V0Ks_standart->size(); v0_index++)
		{
			//if ( V0Ks_standart[v0_index]->getDxy()<0.045cm and dz<0.2 cm) continue;getDxy
			TLorentzVector v0((*V0Ks_standart)[v0_index].px(), (*V0Ks_standart)[v0_index].py(), (*V0Ks_standart)[v0_index].pz(), (*V0Ks_standart)[v0_index].energy());
			Point ref_v0 = (*V0Ks_standart)[v0_index].vertex();
			K892_K0_dxy = getDxy(pv_position, v0, ref_v0);
			K892_K0_dz = getDz(pv_position, v0, ref_v0);
			for (unsigned int strip_i = 0; strip_i < Strips->size(); strip_i++)
			{
				TLorentzVector strip((*Strips)[strip_i].px(), (*Strips)[strip_i].py(), (*Strips)[strip_i].pz(), (*Strips)[strip_i].energy());
				if (strip.M() < 0) continue;
				lv_K892 = strip + v0;
				if (abs(lv_K892.M() - 0.892) <= 0.4) // mass cut?
				{
					nK892_temp++;

					K892_fx = lv_K892.Px();
					K892_fy = lv_K892.Py();
					K892_fz = lv_K892.Pz();
					K892_fE = lv_K892.E();
					K892_fM = lv_K892.M();

					lv_K892_K0 = v0;
					K892_K0_fx = v0.Px();
					K892_K0_fy = v0.Py();
					K892_K0_fz = v0.Pz();
					K892_K0_fE = v0.E();
					K892_K0_fM = v0.M();

					lv_K892_pi0 = strip;
					K892_pi0_fx = strip.Px();
					K892_pi0_fy = strip.Py();
					K892_pi0_fz = strip.Pz();
					K892_pi0_fE = strip.E();
					K892_pi0_fM = strip.M();

					if (v0_index + 1 == V0Ks_standart->size() && strip_i + 1 == Strips->size()) n_K892 = nK892_temp;
					kaon892_tree->Fill();

				}

			}
		}
	}

}
