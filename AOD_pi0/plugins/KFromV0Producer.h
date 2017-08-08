#include "AOD_pi0.h"

void AOD_pi0::KFromV0Producer(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	// Beam spot
	const reco::BeamSpot* theBeamSpot = TheBeamSpotHandle.product();
	math::XYZPoint BSposition(theBeamSpot->position());

	/// V0's RECO Ks's - all are charged - also matching with pi+- of HPS
	if (V0Ks.isValid())
	{
		dlog("RECO Ks's Particles");
		vector<reco::CandidateCollection> v_daughters;
		v0_count = V0Ks->size();
		dout("Size V0:", v0_count, "vs", V0Ks_standart->size());

		// For events with Kaons
		if (v0_count > 0)
		{
			// Initial
			int number_of_passed_ks = 0;
			int number_of_passed_ks_in_hps_tau = 0;
			int number_of_passed_ks_in_hps_tau_hard_cut = 0;
			int number_of_pion_in_tau_jets_with_good_Ks = 0;

			// For Kaons close to PV (0.2):
				//    h_Ks_v0_vx, h_Ks_v0_vy, h_Ks_v0_vz:               Kaon vertex
				//    h_Ks_v0_PV_dXY, h_Ks_v0_dx, h_Ks_v0_dy, h_Ks_v0_dz:  Kaon distance to PV
				//    h_Ks_v0_pions_dR:                                 dR of Ks pions
				//    h_Ks_v0_BS_dXY:                                   distance from BS of Ks v0
				//    h_Ks_v0_daughter_pt:                              pt of pions of Ks
				//    h_Ks_v0_inv_m_pi:                                 inv mass of pions of Ks
				//    h_Ks_v0_number_per_event:                         Number of pions that passed the dz cut
				//    h_Ks_v0_n_ev_passing_dz_cut:                      Number of Events in which at least 1 Kaon passed the selection
				//    h_Ks_v0_n_Ks_in_jets_per_event:                   number_of_passed_ks_in_hps_tau_hard_cut
				//    h_Ks_v0_count:                                    num of kaons if there are any
				//    h_Ks_v0_found_in_hps_tau:                         number_of_passed_ks_in_hps_tau
				//
				//    h_Ks_v0_found_in_hps_tau_dR:                      dR between matched pions and Ks V0
				//    h_Ks_v0_found_in_hps_tau_m_inv:
				//    h_Ks_v0_found_in_hps_tau_dRcut:                   How many Kaons are kinematicaly in jet of tau
				//    h_Ks_v0_n_pion_in_tau_jets_with_good_Ks:          Number of pions in tau jet                 if assuming that there should be a Kaon matched
				//    h_Ks_v0_hps_pions_combinatoric_dR:                dR of all pairs of tau jets                if assuming that there should be a Kaon matched
				//    h_Ks_v0_pions_and_hps_pions_combined_dR:          dR of all pions of tau jet to Ks V0 pions  if assuming that there should be a Kaon matched
				//    h_Ks_v0_found_in_hps_tau_dR_only_one_pion_left:   if one pion matched dR of Ks0 and tau      if assuming that there should be a Kaon matched
			
			// Loop over V0 K0s with matching of V0 pions to HPS pions
			for(unsigned k0_i = 0; k0_i < V0Ks->size(); k0_i++) 
			{
				// In dz Kaons selected as taus: https://twiki.cern.ch/twiki/bin/view/CMS/HiggsToTauTauWorking2016#Baseline_Selection
				if (abs((*V0Ks)[k0_i].vz() - pv_position.z()) > cDZCut) continue; // from PV to Ks

				// Init
					// K->pipi
						TLorentzVector v_Ks((*V0Ks)[k0_i].px(), (*V0Ks)[k0_i].py(), (*V0Ks)[k0_i].pz(), (*V0Ks)[k0_i].energy());
						TLorentzVector v0_ks_pion1((*V0Ks)[k0_i].daughter(0)->px(),(*V0Ks)[k0_i].daughter(0)->py(),(*V0Ks)[k0_i].daughter(0)->pz(),(*V0Ks)[k0_i].daughter(0)->energy());
						TLorentzVector v0_ks_pion2((*V0Ks)[k0_i].daughter(1)->px(),(*V0Ks)[k0_i].daughter(1)->py(),(*V0Ks)[k0_i].daughter(1)->pz(),(*V0Ks)[k0_i].daughter(1)->energy());
					// Possitioning
						double distMagXY_PV_Ks = sqrt(pow(v_Ks.X() - pv_position.x(), 2) + pow(v_Ks.Y() - pv_position.y(), 2));// distnce from PV to Ks in XY
						double distMagXY_BS_Ks = sqrt(pow(v_Ks.X() - BSposition.x(), 2) + pow(v_Ks.Y() - BSposition.y(), 2));// distnce from BS to Ks in XY
						//h_Ks_v0_BS_dXY->Fill(distMagXY_BS_Ks); will be filled if both pions of K will be matched to the hps pions
						//if (distMagXY_PV_Ks > 4.4) continue; // distance from 0.2 -- something is clearly wrong with this cut.
					// Basic plots for Ks V0 collection
						number_of_passed_ks++; // Number of pions that passed the dz cut
						int num = (*V0Ks)[k0_i].numberOfDaughters();//edm::reco::CompositeCandidate::daughters
						int num_moth = (*V0Ks)[k0_i].numberOfMothers();
						v0_ks_numb++;
						Ks_v0_inv_m_pi = (v0_ks_pion1 + v0_ks_pion2).M();
						pt_1 = v0_ks_pion1.Pt();
						pt_2 = v0_ks_pion2.Pt();
						eta_1 = v0_ks_pion1.Eta();
						eta_2 = v0_ks_pion2.Eta();
						max_ks_daughter_pt = max(v0_ks_pion1.Pt(), v0_ks_pion2.Pt());

				// Basic printout
					dout("distMagXY_PV_Ks:", distMagXY_PV_Ks);
					dlog("\tKs", k0_i, "(", (*V0Ks)[k0_i].charge(), ")", "vertex:",  (*V0Ks)[k0_i].vx(), (*V0Ks)[k0_i].vy(),  (*V0Ks)[k0_i].vz());
					dlog("\t\tp:", (*V0Ks)[k0_i].px(), (*V0Ks)[k0_i].py(), (*V0Ks)[k0_i].pz(), "Energy:", (*V0Ks)[k0_i].energy());
					dlog("\t\tdaughters tot number:", num, ";", " moth number:", num_moth);
					dlog("\t\t\t\tpion 1 (", (*V0Ks)[k0_i].daughter(0)->charge(), "); mass:", v0_ks_pion1.M(), "p:", (*V0Ks)[k0_i].daughter(0)->px(),(*V0Ks)[k0_i].daughter(0)->py(), (*V0Ks)[k0_i].daughter(0)->pz(), "Energy:", (*V0Ks)[k0_i].daughter(0)->energy());
					dlog("\t\t\t\tpion 2 (", (*V0Ks)[k0_i].daughter(1)->charge(), "); mass:", v0_ks_pion2.M(), "p:", (*V0Ks)[k0_i].daughter(1)->px(),(*V0Ks)[k0_i].daughter(1)->py(), (*V0Ks)[k0_i].daughter(1)->pz(), "Energy:", (*V0Ks)[k0_i].daughter(1)->energy());

				// Matching to HPS pions pi+-
				if (match_KsV0_to_HPS && PF_hps_taus.isValid())
				{
					dout("\t\t\tPF taus:", PF_hps_taus->size());

					for( unsigned tau_index = 0; tau_index < PF_hps_taus->size(); tau_index++)
					{
						// Init
							reco::PFTauRef pftauref(PF_hps_taus, tau_index);
							TLorentzVector v_tau(pftauref->px(), pftauref->py(), pftauref->pz(), pftauref->energy());
							dout("\t\t\t\tPF tau i", tau_index, "p:", pftauref->px(), pftauref->py(), pftauref->pz(), "Energy:", pftauref->energy(), "dR:", v_tau.DeltaR(v_Ks));
							// Create pions collection of tau
							vector < reco::PFCandidatePtr  > tau_signalPFChargedHadrCands    = pftauref->signalPFChargedHadrCands();//typedef edm::Ptr<PFCandidate> reco::PFCandidatePtr
							vector < reco::PFCandidatePtr  > tau_isolationPFChargedHadrCands = pftauref->isolationPFChargedHadrCands(); // //typedef edm::Ptr<PFCandidate> reco::PFCandidatePtr
							vector < reco::PFCandidatePtr > tau_picharge = pftauref->signalPFChargedHadrCands(); // vector < edm::Ptr<PFCandidate> >
							tau_picharge.insert(tau_picharge.end(), tau_isolationPFChargedHadrCands.begin(), tau_isolationPFChargedHadrCands.end());
							// Pion to pion matching
							bool both_ks_pions_found_in_tau = true;
							int firstfound = -1, secondfound = -1;
							double deltaR1 = 100, deltaR2 = 100;

						// Check if tau jet kinematicaly matches to Ks of V0: tau^ks dr < 0.5
						// Ks are selected close to PV, so dR to tau jet ensures that the Ks of V0 is in tau jet
						// TODO FIND HERE THE MOST CLOSE JET INSTEAD OF THE FIRST PASSING THE CUT
						if (v_tau.DeltaR(v_Ks) > cKtoTauDR) continue;

						if (false) CombinatoricPairs(tau_picharge);

						// Pion to pion matching - the closest ones in dR
						// TODO: think how to remove already matched from previous kaons matching pions
						// TODO: check how many kaons/pions are not matched with this techniq; what if the first pion matches better to the second
						for(unsigned pi_i = 0; pi_i < tau_picharge.size(); pi_i++)
						{
							if (!tau_picharge[pi_i]->bestTrack()) continue;

							// Init
								TLorentzVector tau_pion(tau_picharge[pi_i]->px(), tau_picharge[pi_i]->py(), tau_picharge[pi_i]->pz(), tau_picharge[pi_i]->energy());
								double firstDaugDR = tau_pion.DeltaR(v0_ks_pion1); //TEMP_COMMENTED h_Ks_v0_pions_and_hps_pions_combined_dR->Fill(tau_pion.DeltaR(v0_ks_pion1));
								double secondDaugDR = tau_pion.DeltaR(v0_ks_pion2); //TEMP_COMMENTED h_Ks_v0_pions_and_hps_pions_combined_dR->Fill(tau_pion.DeltaR(v0_ks_pion2));
								//if (tau_pion.M() < 0.1) break; //ISSUE

							// Best possible matching of tau jets pions to Ks V0 pions
							// if matched to first pion of Ks V0 and it is matched better than the previous candidate
							if (BetterPionMatch(tau_pion, tau_picharge[pi_i]->charge(), v0_ks_pion1, (*V0Ks)[k0_i].daughter(0)->charge(), firstfound, deltaR1))
							{
								firstfound = pi_i;
								dlog("dR:", v_tau.DeltaR(v_Ks), "But first pion found", "\ntau pion1:", tau_picharge[pi_i]->px(), tau_picharge[pi_i]->py(), tau_picharge[pi_i]->pz(), tau_picharge[pi_i]->energy());
							}
							else if (BetterPionMatch(tau_pion, tau_picharge[pi_i]->charge(), v0_ks_pion2, (*V0Ks)[k0_i].daughter(1)->charge(), secondfound, deltaR2))
							{
								secondfound = pi_i;
							}
						}

						// Invariant mass of two matched pions
						if (firstfound != -1 && secondfound != -1)
						{
							// init
								TLorentzVector v_tau_picharge_firstfound  = TLorentzVector(tau_picharge[firstfound]->px(), tau_picharge[firstfound]->py(), tau_picharge[firstfound]->pz(), tau_picharge[firstfound]->energy());
								TLorentzVector v_tau_picharge_secondfound = TLorentzVector(tau_picharge[secondfound]->px(), tau_picharge[secondfound]->py(), tau_picharge[secondfound]->pz(), tau_picharge[secondfound]->energy());
								TLorentzVector v_Ks_from_hps_pions = v_tau_picharge_firstfound + v_tau_picharge_secondfound;
								double invMassHPSpions = v_Ks_from_hps_pions.M();
								number_of_passed_ks_in_hps_tau++;

							// fill histograms
								//TEMP_COMMENTED h_Ks_v0_found_in_hps_tau_dR->Fill(TLorentzVector(tau_picharge[firstfound]->px(), tau_picharge[firstfound]->py(), tau_picharge[firstfound]->pz(), tau_picharge[firstfound]->energy()).DeltaR(v_Ks));//   h_Ks_v0_found_in_hps_tau_dR->Fill(v_tau_picharge_firstfound.DeltaR(v_Ks));
								//TEMP_COMMENTED h_Ks_v0_found_in_hps_tau_dR->Fill(TLorentzVector(tau_picharge[secondfound]->px(), tau_picharge[secondfound]->py(), tau_picharge[secondfound]->pz(), tau_picharge[secondfound]->energy()).DeltaR(v_Ks));//   h_Ks_v0_found_in_hps_tau_dR->Fill(v_tau_picharge_secondfound.DeltaR(v_Ks));
							
							nPionsInJetsWithKs = tau_picharge.size();
							v_HPS_Ks_inv_m_pi.push_back(invMassHPSpions);
							v_HPS_Ks_DR.push_back(v_tau.DeltaR(v_Ks));
							v_HPS_pt_1.push_back(v_tau_picharge_firstfound.Pt());
							v_HPS_pt_2.push_back(v_tau_picharge_secondfound.Pt());
							v_HPS_eta_1.push_back(v_tau_picharge_firstfound.Eta());
							v_HPS_eta_2.push_back(v_tau_picharge_secondfound.Eta());
							v_HPS_Ks_pions_DR.push_back(v_tau_picharge_firstfound.DeltaR(v_tau_picharge_secondfound));
							v_v0_matched_pt_1.push_back(v0_ks_pion1.Pt());
							v_v0_matched_pt_2.push_back(v0_ks_pion2.Pt());
							v_v0_matched_eta_1.push_back(v0_ks_pion1.Eta());
							v_v0_matched_eta_2.push_back(v0_ks_pion2.Eta());
							v_nPionsInJetsWithKs.push_back(tau_picharge.size());
							break; // Since Kaon should be associated to only one jet anyway
						}
						//TEMP_COMMENTED else if (firstfound == -1 || secondfound == -1) h_Ks_v0_found_in_hps_tau_dR_only_one_pion_left->Fill(v_tau.DeltaR(v_Ks));
						
					}
				}

				// Fill branches per-Kaon
				tree->Fill();

				v_v0_Ks_pions_DR.push_back(v0_ks_pion1.DeltaR(v0_ks_pion2));
				v_v0_Ks_inv_m_pi.push_back((v0_ks_pion1 + v0_ks_pion2).M());
				v_v0_Ks_DR.push_back(v0_ks_pion1.DeltaR(v0_ks_pion2));
				v_v0_pt_1.push_back(v0_ks_pion1.Pt());
				v_v0_pt_2.push_back(v0_ks_pion2.Pt());
				v_v0_eta_1.push_back(v0_ks_pion1.Eta());
				v_v0_eta_2.push_back(v0_ks_pion2.Eta());

				// Fill basic hists
					//TEMP_COMMENTED h_Ks_v0_vx->Fill(v_Ks.X());
					//TEMP_COMMENTED h_Ks_v0_vy->Fill(v_Ks.Y());
					//TEMP_COMMENTED h_Ks_v0_vz->Fill(v_Ks.Z());
					//TEMP_COMMENTED h_Ks_v0_dx->Fill(abs(v_Ks.X() - pv_position.x()));
					//TEMP_COMMENTED h_Ks_v0_dy->Fill(abs(v_Ks.Y() - pv_position.y()));
					//TEMP_COMMENTED h_Ks_v0_dz->Fill(abs(v_Ks.Z() - pv_position.z()));
			}

			// Fill histograms
				//TEMP_COMMENTED h_Ks_v0_n_ev_passing_dz_cut->Fill(1 * at_least_one_Ks_passed_dz);
				//TEMP_COMMENTED h_Ks_v0_n_Ks_in_jets_per_event->Fill(number_of_passed_ks_in_hps_tau_hard_cut);
				//TEMP_COMMENTED h_Ks_v0_count->Fill(v0_count); //h_v0_count->Fill(v0_count);
				//TEMP_COMMENTED h_Ks_v0_number_per_event->Fill(number_of_passed_ks);
				// ks from V0 coll, NPE that passed PV and found in HPS tau jets
				//TEMP_COMMENTED if (number_of_passed_ks_in_hps_tau > 0) h_Ks_v0_found_in_hps_tau->Fill(number_of_passed_ks_in_hps_tau);

			dlog("v0_count", v0_count,"number_of_passed_ks", number_of_passed_ks, "number_of_passed_ks_in_hps_tau:", number_of_passed_ks_in_hps_tau, "number_of_passed_ks_in_hps_tau_hard_cut", number_of_passed_ks_in_hps_tau_hard_cut);
		}

		// count the number of the Ks from standart collection in the dz < cDZCut [0.2] cm from the PV 
		if (false)
		{  
			int N_standart_KS = 0;
			if (V0Ks_standart->size() > 0)
				for(unsigned k0_i = 0; k0_i < V0Ks_standart->size(); k0_i++) // over V0 K0s from the PV
				{
					if (abs((*V0Ks_standart)[k0_i].vz() - pv_position.z()) > cDZCut) continue; // 
					N_standart_KS++;
				}
			if (N_standart_KS > 0) dout("\t\tstandart Ks in dZ:", N_standart_KS);
		}

		//  For all the tau jets :
		//    h_Ks_v0_n_pion_in_tau_jets:N of pions in jet
		//    h_Ks_v0_pions_and_hps_pions_combined_dR_with_no_constrain: for each Kaon's pion dR to each tau jet pion
		if (false)
		{
			//TEMP_COMMENTED h_Ks_v0_n_tau_jets_per_event->Fill(PF_hps_taus->size());
			for( unsigned tau_index = 0; tau_index < PF_hps_taus->size(); tau_index++)
			{
				reco::PFTauRef pftauref(PF_hps_taus, tau_index);
				TLorentzVector v_tau(pftauref->px(), pftauref->py(), pftauref->pz(), pftauref->energy());
				double distMagXY = sqrt(pow(v_tau.X() - pv_position.x(),2) + pow(v_tau.Y() - pv_position.y(),2));
				//TEMP_COMMENTED h_tau_v0_dXY->Fill(distMagXY);

				//pi+-
					vector < reco::PFCandidatePtr  > tau_signalPFChargedHadrCands    = pftauref->signalPFChargedHadrCands();//typedef edm::Ptr<PFCandidate> reco::PFCandidatePtr
					vector < reco::PFCandidatePtr  > tau_isolationPFChargedHadrCands = pftauref->isolationPFChargedHadrCands(); // //typedef edm::Ptr<PFCandidate> reco::PFCandidatePtr
				// All pi+-'s
				vector < reco::PFCandidatePtr > tau_picharge = pftauref->signalPFChargedHadrCands(); // vector < edm::Ptr<PFCandidate> >
																				tau_picharge.insert(tau_picharge.end(), tau_isolationPFChargedHadrCands.begin(), tau_isolationPFChargedHadrCands.end());
				//TEMP_COMMENTED h_Ks_v0_n_pion_in_tau_jets->Fill(tau_picharge.size());

				if (v0_count > 0)
				{
					for(unsigned k0_i = 0 ; k0_i < V0Ks->size() ; k0_i++)
					{
						TLorentzVector v0_ks_pion1((*V0Ks)[k0_i].daughter(0)->px(),(*V0Ks)[k0_i].daughter(0)->py(),(*V0Ks)[k0_i].daughter(0)->pz(),(*V0Ks)[k0_i].daughter(0)->energy());
						TLorentzVector v0_ks_pion2((*V0Ks)[k0_i].daughter(1)->px(),(*V0Ks)[k0_i].daughter(1)->py(),(*V0Ks)[k0_i].daughter(1)->pz(),(*V0Ks)[k0_i].daughter(1)->energy());

						for(unsigned pi_i = 0; pi_i < tau_picharge.size(); pi_i++)
						{

							TLorentzVector tau_pion(tau_picharge[pi_i]->px(), tau_picharge[pi_i]->py(), tau_picharge[pi_i]->pz(), tau_picharge[pi_i]->energy());

							//TEMP_COMMENTED h_Ks_v0_pions_and_hps_pions_combined_dR_with_no_constrain->Fill(tau_pion.DeltaR(v0_ks_pion1));
							//TEMP_COMMENTED h_Ks_v0_pions_and_hps_pions_combined_dR_with_no_constrain->Fill(tau_pion.DeltaR(v0_ks_pion2));
						}
					}
				}

				//dR of all tau pions pairs
					for(unsigned pi_i = 0; pi_i < tau_picharge.size() - 1 && tau_picharge.size() > 1; pi_i++)
					{
						TLorentzVector tau_pion(tau_picharge[pi_i]->px(), tau_picharge[pi_i]->py(), tau_picharge[pi_i]->pz(), tau_picharge[pi_i]->energy());
						for(unsigned pi_i2 = pi_i + 1; pi_i2 < tau_picharge.size() ; pi_i2++)
						{
							TLorentzVector tau_pion2(tau_picharge[pi_i2]->px(), tau_picharge[pi_i2]->py(), tau_picharge[pi_i2]->pz(), tau_picharge[pi_i2]->energy());
							//TEMP_COMMENTED h_Ks_v0_hps_pions_combinatoric_dR->Fill((tau_pion).DeltaR(tau_pion2));

							//TEMP_COMMENTED h_tau_comb_pions_m_inv->Fill((tau_pion + tau_pion2).M());
						}
					}
			}
		}
	}
	else dlog("UNVALID Ks's");
}

template  <typename VectoreType >
void AOD_pi0::CombinatoricPairs(vector<VectoreType> pionCollection)
{
	//dR of all tau pions pairs
	for(unsigned pi_i = 0; pi_i < pionCollection.size() - 1 && pionCollection.size() > 1; pi_i++)
	{
		TLorentzVector tau_pion(pionCollection[pi_i]->px(), pionCollection[pi_i]->py(), pionCollection[pi_i]->pz(), pionCollection[pi_i]->energy());
		for(unsigned pi_i2 = pi_i + 1; pi_i2 < pionCollection.size() ; pi_i2++)
		{
			TLorentzVector tau_pion2(pionCollection[pi_i2]->px(), pionCollection[pi_i2]->py(), pionCollection[pi_i2]->pz(), pionCollection[pi_i2]->energy());
			
			v_KsCombinatoricMass.push_back((tau_pion + tau_pion2).M());
			v_KsCombinatoricDR.push_back(tau_pion.DeltaR(tau_pion2));
		}
	}
}

bool AOD_pi0::BetterPionMatch(TLorentzVector tau_pion, int chargeTauPion, TLorentzVector v0_ks_pion, int chargeKPion, int firstfound, double & dR)
{
	if ((chargeTauPion == chargeKPion) &&
		tau_pion.DeltaR(v0_ks_pion) < 0.5 &&// TODO: Can be better tuned
		(firstfound == -1 || dR > tau_pion.DeltaR(v0_ks_pion)))
		{
			dR = tau_pion.DeltaR(v0_ks_pion);
			return true;
		}
	return false;
}
