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
		v0_Ks_count = V0Ks->size();
		dout("Size V0:", v0_Ks_count, "vs", V0Ks_standart->size());
		bool onexit = false;
		// For events with Kaons
		if (v0_Ks_count > 0)
		{
			numOfUnmatchedKaons = 0;
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
				ResetBranchesKaonTree();
				unsigned int v0NegPionIndex(0), v0PosPionIndex(1);
				if ((*V0Ks)[k0_i].daughter(0)->charge()>0)
				{
					v0NegPionIndex = 1; v0PosPionIndex = 0;
				}
				dout( "From v0 pion ", k0_i, "track:", ((reco::RecoChargedCandidate *)((*V0Ks)[k0_i].daughter(v0NegPionIndex)))->track()->px(),
					((reco::RecoChargedCandidate *)((*V0Ks)[k0_i].daughter(v0NegPionIndex)))->track()->py(),
					((reco::RecoChargedCandidate *)((*V0Ks)[k0_i].daughter(v0NegPionIndex)))->track()->pz()
					);
				dout( "\t\t LV:", (*V0Ks)[k0_i].daughter(v0NegPionIndex)->px(),
					(*V0Ks)[k0_i].daughter(v0NegPionIndex)->py(),
					(*V0Ks)[k0_i].daughter(v0NegPionIndex)->pz(),
					(*V0Ks)[k0_i].daughter(v0NegPionIndex)->energy()
					);
				dout();
				dout( "From v0 pion ", k0_i, "track:", ((reco::RecoChargedCandidate *)((*V0Ks)[k0_i].daughter(v0PosPionIndex)))->track()->px(),
					((reco::RecoChargedCandidate *)((*V0Ks)[k0_i].daughter(v0PosPionIndex)))->track()->py(),
					((reco::RecoChargedCandidate *)((*V0Ks)[k0_i].daughter(v0PosPionIndex)))->track()->pz()
					);
				dout( "\t\t LV:", (*V0Ks)[k0_i].daughter(v0PosPionIndex)->px(),
					(*V0Ks)[k0_i].daughter(v0PosPionIndex)->py(),
					(*V0Ks)[k0_i].daughter(v0PosPionIndex)->pz(),
					(*V0Ks)[k0_i].daughter(v0PosPionIndex)->energy()
					);

				// In dz Kaons selected as taus: https://twiki.cern.ch/twiki/bin/view/CMS/HiggsToTauTauWorking2016#Baseline_Selection
				if (abs((*V0Ks)[k0_i].vz() - pv_position.z()) > cDZCut) continue; // from PV to Ks

				// Init
					// K->pipi
						TLorentzVector Ks((*V0Ks)[k0_i].px(), (*V0Ks)[k0_i].py(), (*V0Ks)[k0_i].pz(), (*V0Ks)[k0_i].energy());
						TLorentzVector v0NegPion((*V0Ks)[k0_i].daughter(v0NegPionIndex)->px(),(*V0Ks)[k0_i].daughter(v0NegPionIndex)->py(),(*V0Ks)[k0_i].daughter(v0NegPionIndex)->pz(),(*V0Ks)[k0_i].daughter(v0NegPionIndex)->energy());
						TLorentzVector v0PosPion((*V0Ks)[k0_i].daughter(v0PosPionIndex)->px(),(*V0Ks)[k0_i].daughter(v0PosPionIndex)->py(),(*V0Ks)[k0_i].daughter(v0PosPionIndex)->pz(),(*V0Ks)[k0_i].daughter(v0PosPionIndex)->energy());
					// Possitioning
						double distMagXY_PV_Ks = sqrt(pow(Ks.X() - pv_position.x(), 2) + pow(Ks.Y() - pv_position.y(), 2));// distnce from PV to Ks in XY
						double distMagXY_BS_Ks = sqrt(pow(Ks.X() - BSposition.x(), 2) + pow(Ks.Y() - BSposition.y(), 2));// distnce from BS to Ks in XY
						//h_Ks_v0_BS_dXY->Fill(distMagXY_BS_Ks); will be filled if both pions of K will be matched to the hps pions
						//if (distMagXY_PV_Ks > 4.4) continue; // distance from 0.2 -- something is clearly wrong with this cut.
					// Basic plots for Ks V0 collection
						number_of_passed_ks++; // Number of pions that passed the dz cut
						int num = (*V0Ks)[k0_i].numberOfDaughters();//edm::reco::CompositeCandidate::daughters
						int num_moth = (*V0Ks)[k0_i].numberOfMothers();
						v0_ks_numb++;

				// Basic printout
					dout("\tdistMagXY_PV_Ks:", distMagXY_PV_Ks);
					dlog("\tKs", k0_i, "(", (*V0Ks)[k0_i].charge(), ")", "vertex:",  (*V0Ks)[k0_i].vx(), (*V0Ks)[k0_i].vy(),  (*V0Ks)[k0_i].vz());
					dlog("\t\tp:", (*V0Ks)[k0_i].px(), (*V0Ks)[k0_i].py(), (*V0Ks)[k0_i].pz(), "Energy:", (*V0Ks)[k0_i].energy());
					dlog("\t\tdaughters tot number:", num, ";", " moth number:", num_moth);
					dlog("\t\t\t\tpion 1 (", (*V0Ks)[k0_i].daughter(v0NegPionIndex)->charge(), "); mass:", v0NegPion.M(), "p:", (*V0Ks)[k0_i].daughter(v0NegPionIndex)->px(),(*V0Ks)[k0_i].daughter(v0NegPionIndex)->py(), (*V0Ks)[k0_i].daughter(v0NegPionIndex)->pz(), "Energy:", (*V0Ks)[k0_i].daughter(v0NegPionIndex)->energy());
					dlog("\t\t\t\tpion 2 (", (*V0Ks)[k0_i].daughter(v0PosPionIndex)->charge(), "); mass:", v0PosPion.M(), "p:", (*V0Ks)[k0_i].daughter(v0PosPionIndex)->px(),(*V0Ks)[k0_i].daughter(v0PosPionIndex)->py(), (*V0Ks)[k0_i].daughter(v0PosPionIndex)->pz(), "Energy:", (*V0Ks)[k0_i].daughter(v0PosPionIndex)->energy());

				// Matching to HPS pions pi+-
				if (match_KsV0_to_HPS && PF_hps_taus.isValid())
				{
					dout("\t\t\tPF taus:", PF_hps_taus->size());

					for( unsigned tau_index = 0; tau_index < PF_hps_taus->size(); tau_index++)
					{
						// Init
							reco::PFTauRef pftauref(PF_hps_taus, tau_index);
							TLorentzVector tau(pftauref->px(), pftauref->py(), pftauref->pz(), pftauref->energy());
							dout("\t\t\t\tPF tau i", tau_index, "p:", pftauref->px(), pftauref->py(), pftauref->pz(), "Energy:", pftauref->energy(), "dR:", tau.DeltaR(Ks));
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
						if (tau.DeltaR(Ks) > cKtoTauDR) continue;

						if (false) CombinatoricPairs(tau_picharge);

						// Pion to pion matching - the closest ones in dR
						// TODO: think if to remove already matched from previous kaons matching pions
						// TODO: check how many kaons/pions are not matched with this techniq; what if the first pion matches better to the second - done 50%
						for(unsigned pi_i = 0; pi_i < tau_picharge.size(); pi_i++)
						{
							if (!tau_picharge[pi_i]->bestTrack()) continue;
							if (k0_i == 0) dlog("\ttau pion:", tau_picharge[pi_i]->px(), tau_picharge[pi_i]->py(), tau_picharge[pi_i]->pz(), sqrt(pow(tau_picharge[pi_i]->energy(), 2) + pow(0.13957018, 2)) );
							if (k0_i == 0) dlog("\ttau pion: bestTrack: ", tau_picharge[pi_i]->bestTrack()->px(), tau_picharge[pi_i]->bestTrack()->py(), tau_picharge[pi_i]->bestTrack()->pz());

							// Best possible matching of tau jets pions to Ks V0 pions
							if (matchByReference)
							{
								// Checks by references that the dR of v0 pion associated track and hps track is 0
								if (firstfound != -1 && secondfound != -1) break;

								if (tau_picharge[pi_i]->charge() == (*V0Ks)[k0_i].daughter(v0NegPionIndex)->charge() && PionMatchByRefference(tau_picharge[pi_i], (*V0Ks)[k0_i].daughter(v0NegPionIndex)))
									firstfound = pi_i;
								else if (tau_picharge[pi_i]->charge() == (*V0Ks)[k0_i].daughter(v0PosPionIndex)->charge() && PionMatchByRefference(tau_picharge[pi_i], (*V0Ks)[k0_i].daughter(v0PosPionIndex)))
									secondfound = pi_i;
								else continue;
							}
							else
							{
								// Init
									TLorentzVector tau_pion(tau_picharge[pi_i]->px(), tau_picharge[pi_i]->py(), tau_picharge[pi_i]->pz(), tau_picharge[pi_i]->energy());
									double firstDaugDR = tau_pion.DeltaR(v0NegPion); //TEMP_COMMENTED h_Ks_v0_pions_and_hps_pions_combined_dR->Fill(tau_pion.DeltaR(v0NegPion));
									double secondDaugDR = tau_pion.DeltaR(v0PosPion); //TEMP_COMMENTED h_Ks_v0_pions_and_hps_pions_combined_dR->Fill(tau_pion.DeltaR(v0PosPion));
									//if (tau_pion.M() < 0.1) break; //ISSUE

								// if matched to first pion of Ks V0 and it is matched better than the previous candidate
								if (BetterPionMatch(tau_pion, tau_picharge[pi_i]->charge(), v0NegPion, (*V0Ks)[k0_i].daughter(v0NegPionIndex)->charge(), firstfound, deltaR1))
									firstfound = pi_i;
								else if (BetterPionMatch(tau_pion, tau_picharge[pi_i]->charge(), v0PosPion, (*V0Ks)[k0_i].daughter(v0PosPionIndex)->charge(), secondfound, deltaR2))
									secondfound = pi_i;
								else continue;
							}
						}

						if (firstfound != -1 && secondfound != -1)
						{
							dout("\tFillPion");
							FillPion(tau_picharge[firstfound], (*V0Ks)[k0_i].daughter(v0NegPionIndex));
							FillPion(tau_picharge[secondfound], (*V0Ks)[k0_i].daughter(v0PosPionIndex));
							// init
								TLorentzVector hpsNegPion  = TLorentzVector(tau_picharge[firstfound]->px(), tau_picharge[firstfound]->py(), tau_picharge[firstfound]->pz(), tau_picharge[firstfound]->energy());
								TLorentzVector hpsPosPion = TLorentzVector(tau_picharge[secondfound]->px(), tau_picharge[secondfound]->py(), tau_picharge[secondfound]->pz(), tau_picharge[secondfound]->energy());
								TLorentzVector hpsKs = hpsNegPion + hpsPosPion;

								number_of_passed_ks_in_hps_tau++;
								TLorentzVector v0_di_pion = v0NegPion + v0PosPion;

							dout("\tTEST V0:", v0_di_pion.X(), v0_di_pion.Y(), v0_di_pion.Z(), v0_di_pion.M(), v0_di_pion.E());
							dout("\t\t", sqrt(pow(v0_di_pion.E(), 2) - (pow(v0_di_pion.X(),2) + pow(v0_di_pion.Y(),2) + pow(v0_di_pion.Z(),2) )));
							dout("\tTESTHPS:", hpsKs.X(), hpsKs.Y(), hpsKs.Z(), hpsKs.M(), hpsKs.E());
							dout("\t\t", sqrt(pow(hpsKs.E(), 2) - (pow(hpsKs.X(),2) + pow(hpsKs.Y(),2) + pow(hpsKs.Z(),2) )));
							if (hpsKs.M() < 0.43) onexit = true;

							// tree_kaon
								nPionsInJetsWithKs = tau_picharge.size();
								// V0 comparison to HPS
								v0hps_pions_dPhi_1 = hpsNegPion.DeltaPhi(v0NegPion);
								v0hps_pions_dEta_1 = hpsNegPion.Eta() - v0NegPion.Eta();
								v0hps_pions_dPhi_2 = hpsPosPion.DeltaPhi(v0PosPion);
								v0hps_pions_dEta_2 = hpsPosPion.Eta() - v0PosPion.Eta();
								v0hps_KsTau_DR = tau.DeltaR(Ks);
								v0hps_KsDiPion_DR = hpsKs.DeltaR(Ks);
								// HPS Kaon
								hps_Ks_inv_m_pi = hpsKs.M();
								hps_pions_DR = hpsNegPion.DeltaR(hpsPosPion);
								hps_deta = hpsNegPion.Eta() - hpsPosPion.Eta();
								hps_dphi = hpsNegPion.DeltaPhi(hpsPosPion);
								hps_px = hpsKs.X();
								hps_py = hpsKs.Y();
								hps_pz = hpsKs.Z();
								hps_energy = hpsKs.E();
								// HPS neg pion
								hps_pt_1 = hpsNegPion.Pt();
								hps_eta_1 = hpsNegPion.Eta();
								hps_phi_1 = hpsNegPion.Phi();
								hps_px_1 = hpsNegPion.X();
								hps_py_1 = hpsNegPion.Y();
								hps_pz_1 = hpsNegPion.Z();
								hps_energy_1 = hpsNegPion.E();
								// HPS pos pion
								hps_pt_2 = hpsPosPion.Pt();
								hps_eta_2 = hpsPosPion.Eta();
								hps_phi_2 = hpsPosPion.Phi();
								hps_px_2 = hpsPosPion.X();
								hps_py_2 = hpsPosPion.Y();
								hps_pz_2 = hpsPosPion.Z();
								hps_energy_2 = hpsPosPion.E();
								// V0 kaon
								v0_Ks_inv_m_pi = (v0NegPion + v0PosPion).M();
								v0_Ks_pions_DR = v0NegPion.DeltaR(v0PosPion);
								v0_deta = v0NegPion.Eta() - v0PosPion.Eta();
								v0_dphi = v0NegPion.DeltaPhi(v0PosPion);
								v0_pt = Ks.Pt();
								v0_eta = Ks.Eta();
								v0_phi = Ks.Phi();
								v0_px = Ks.X();
								v0_py = Ks.Y();
								v0_pz = Ks.Z();
								v0_energy = Ks.E();
								// V0 neg pion
								v0_pt_1 = v0NegPion.Pt();
								v0_eta_1 = v0NegPion.Eta();
								v0_phi_1 = v0NegPion.Phi();
								v0_px_1 = v0NegPion.X();
								v0_py_1 = v0NegPion.Y();
								v0_pz_1 = v0NegPion.Z();
								v0_energy_1 = v0NegPion.E();
								// V0 pos pion
								v0_pt_2 = v0PosPion.Pt();
								v0_eta_2 = v0PosPion.Eta();
								v0_phi_2 = v0PosPion.Phi();
								v0_px_2 = v0PosPion.X();
								v0_py_2 = v0PosPion.Y();
								v0_pz_2 = v0PosPion.Z();
								v0_energy_2 = v0PosPion.E();

							// tree_once
								v_hps_Ks_inv_m_pi.push_back(hps_Ks_inv_m_pi);
								v_hps_Ks_DR.push_back(v0hps_KsTau_DR);
								v_hps_pt_1.push_back(hps_pt_1);
								v_hps_pt_2.push_back(hps_pt_2);
								v_hps_eta_1.push_back(hps_eta_1);
								v_hps_eta_2.push_back(hps_eta_2);
								v_hps_Ks_pions_DR.push_back(hps_pions_DR);
								v_v0_matched_pt_1.push_back(v0_pt_1);
								v_v0_matched_pt_2.push_back(v0_pt_2);
								v_v0_matched_eta_1.push_back(v0_eta_1);
								v_v0_matched_eta_2.push_back(v0_eta_2);
								v_v0_Ks_inv_m_pi.push_back(v0_Ks_inv_m_pi);
								v_v0_pions_DR.push_back(v0_Ks_pions_DR);
								v_v0_pt_1.push_back(v0_pt_1);
								v_v0_pt_2.push_back(v0_pt_2);
								v_v0_eta_1.push_back(v0_eta_1);
								v_v0_eta_2.push_back(v0_eta_2);
								v_nPionsInJetsWithKs.push_back(hps_eta_2);

								kaon_tree->Fill();
							break; // Since Kaon should be associated to only one jet anyway
						}
						else numOfUnmatchedKaons++;
					}
				}
			}

			// Fill histograms
				//TEMP_COMMENTED h_Ks_v0_n_ev_passing_dz_cut->Fill(1 * at_least_one_Ks_passed_dz);
				//TEMP_COMMENTED h_Ks_v0_n_Ks_in_jets_per_event->Fill(number_of_passed_ks_in_hps_tau_hard_cut);
				//TEMP_COMMENTED h_Ks_v0_count->Fill(v0_count); //h_v0_count->Fill(v0_count);
				//TEMP_COMMENTED h_Ks_v0_number_per_event->Fill(number_of_passed_ks);
				// ks from V0 coll, NPE that passed PV and found in HPS tau jets
				//TEMP_COMMENTED if (number_of_passed_ks_in_hps_tau > 0) h_Ks_v0_found_in_hps_tau->Fill(number_of_passed_ks_in_hps_tau);

			dlog("v0_Ks_count", v0_Ks_count,"number_of_passed_ks", number_of_passed_ks, "number_of_passed_ks_in_hps_tau:", number_of_passed_ks_in_hps_tau, "number_of_passed_ks_in_hps_tau_hard_cut", number_of_passed_ks_in_hps_tau_hard_cut);
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
				TLorentzVector tau(pftauref->px(), pftauref->py(), pftauref->pz(), pftauref->energy());
				double distMagXY = sqrt(pow(tau.X() - pv_position.x(),2) + pow(tau.Y() - pv_position.y(),2));
				//TEMP_COMMENTED h_tau_v0_dXY->Fill(distMagXY);

				//pi+-
					vector < reco::PFCandidatePtr  > tau_signalPFChargedHadrCands    = pftauref->signalPFChargedHadrCands();//typedef edm::Ptr<PFCandidate> reco::PFCandidatePtr
					vector < reco::PFCandidatePtr  > tau_isolationPFChargedHadrCands = pftauref->isolationPFChargedHadrCands(); // //typedef edm::Ptr<PFCandidate> reco::PFCandidatePtr
				// All pi+-'s
				vector < reco::PFCandidatePtr > tau_picharge = pftauref->signalPFChargedHadrCands(); // vector < edm::Ptr<PFCandidate> >
																				tau_picharge.insert(tau_picharge.end(), tau_isolationPFChargedHadrCands.begin(), tau_isolationPFChargedHadrCands.end());
				//TEMP_COMMENTED h_Ks_v0_n_pion_in_tau_jets->Fill(tau_picharge.size());

				if (v0_Ks_count > 0)
				{
					for(unsigned k0_i = 0 ; k0_i < V0Ks->size() ; k0_i++)
					{
						unsigned int v0NegPionIndex(0), v0PosPionIndex(1);
						if ((*V0Ks)[k0_i].daughter(0)->charge()>0)
						{
							v0NegPionIndex = 1; v0PosPionIndex = 0;
						}
						TLorentzVector v0NegPion((*V0Ks)[k0_i].daughter(v0NegPionIndex)->px(),(*V0Ks)[k0_i].daughter(v0NegPionIndex)->py(),(*V0Ks)[k0_i].daughter(v0NegPionIndex)->pz(),(*V0Ks)[k0_i].daughter(v0NegPionIndex)->energy());
						TLorentzVector v0PosPion((*V0Ks)[k0_i].daughter(v0PosPionIndex)->px(),(*V0Ks)[k0_i].daughter(v0PosPionIndex)->py(),(*V0Ks)[k0_i].daughter(v0PosPionIndex)->pz(),(*V0Ks)[k0_i].daughter(v0PosPionIndex)->energy());

						for(unsigned pi_i = 0; pi_i < tau_picharge.size(); pi_i++)
						{

							TLorentzVector tau_pion(tau_picharge[pi_i]->px(), tau_picharge[pi_i]->py(), tau_picharge[pi_i]->pz(), tau_picharge[pi_i]->energy());

							//TEMP_COMMENTED h_Ks_v0_pions_and_hps_pions_combined_dR_with_no_constrain->Fill(tau_pion.DeltaR(v0NegPion));
							//TEMP_COMMENTED h_Ks_v0_pions_and_hps_pions_combined_dR_with_no_constrain->Fill(tau_pion.DeltaR(v0PosPion));
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
		//if (onexit) exit(1);
		if (onexit) dout("OHNO");
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

template  <typename HPSPion, typename Daughter>
bool AOD_pi0::PionMatchByRefference(HPSPion pion, Daughter daughter)
{
	TLorentzVector tau_pion, v0_pion;
	v0_pion.SetXYZM( ((reco::RecoChargedCandidate *)daughter)->track()->px(), ((reco::RecoChargedCandidate *)daughter)->track()->py(), ((reco::RecoChargedCandidate *)daughter)->track()->pz(), 0.13957018);
	tau_pion.SetXYZM(pion->bestTrack()->px(), pion->bestTrack()->py(), pion->bestTrack()->pz(), 0.13957018);

	if (v0_pion.DeltaR(tau_pion) == 0) return true;
	return false;
}

void AOD_pi0::ResetBranchesKaonTree()
{
	nPionsInJetsWithKs = -999;
	// V0 comparison to HPS
	v0hps_pions_dPhi_1 = -999;
	v0hps_pions_dEta_1 = -999;
	v0hps_pions_dPhi_2 = -999;
	v0hps_pions_dEta_2 = -999;
	v0hps_KsTau_DR = -999;
	v0hps_KsDiPion_DR = -999;
	// HPS Kaon
	hps_Ks_inv_m_pi = -999;
	hps_pions_DR = -999;
	hps_deta = -999;
	hps_dphi = -999;
	hps_px = -999;
	hps_py = -999;
	hps_pz = -999;
	hps_energy = -999;
	// HPS neg pion
	hps_pt_1 = -999;
	hps_eta_1 = -999;
	hps_phi_1 = -999;
	hps_px_1 = -999;
	hps_py_1 = -999;
	hps_pz_1 = -999;
	hps_energy_1 = -999;
	// HPS pos pion
	hps_pt_2 = -999;
	hps_eta_2 = -999;
	hps_phi_2 = -999;
	hps_px_2 = -999;
	hps_py_2 = -999;
	hps_pz_2 = -999;
	hps_energy_2 = -999;
	// V0 kaon
	v0_Ks_inv_m_pi = -999;
	v0_Ks_pions_DR = -999;
	v0_deta = -999;
	v0_dphi = -999;
	v0_pt = -999;
	v0_eta = -999;
	v0_phi = -999;
	v0_px = -999;
	v0_py = -999;
	v0_pz = -999;
	v0_energy = -999;
	// V0 neg pion
	v0_pt_1 = -999;
	v0_eta_1 = -999;
	v0_phi_1 = -999;
	v0_px_1 = -999;
	v0_py_1 = -999;
	v0_pz_1 = -999;
	v0_energy_1 = -999;
	// V0 pos pion
	v0_pt_2 = -999;
	v0_eta_2 = -999;
	v0_phi_2 = -999;
	v0_px_2 = -999;
	v0_py_2 = -999;
	v0_pz_2 = -999;
	v0_energy_2 = -999;
}
