#include "AOD_pi0.h"

void AOD_pi0::K892(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
		///K(892)+- TODO: build using V0 and hps pi0; made: patching to pi+- of hps and V0
	if (false && PF_hps_taus.isValid() )
		{
			vector<reco::CandidateCollection> v0_daughters;
			v0_count = V0Ks->size();
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
					//dout("\tpion:", pi_i, tau_picharge[pi_i]->charge(), tau_picharge[pi_i]->px(), tau_picharge[pi_i]->py(), tau_picharge[pi_i]->pz(), tau_picharge[pi_i]->energy());
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

}