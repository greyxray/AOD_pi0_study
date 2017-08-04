#include "AOD_pi0.h"

void AOD_pi0::KFromV0Producer(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	
	const reco::BeamSpot* theBeamSpot = TheBeamSpotHandle.product();
	math::XYZPoint BSposition(theBeamSpot->position());
	/// V0's RECO Ks's - all are charged - also matching with pi+- of HPS
	if (true && V0Ks.isValid())
	{
		dlog("RECO Ks's Particles");
		vector<reco::CandidateCollection> v_daughters;
		v0_count = V0Ks->size();
		dout("Size V0:", v0_count, "vs", V0Ks_standart->size());

		// For events with Kaons
		if (v0_count > 0)
		{
			int number_of_passed_ks = 0;
			int number_of_passed_ks_in_hps_tau = 0;
			int number_of_passed_ks_in_hps_tau_hard_cut = 0;
			int number_of_pion_in_tau_jets_with_good_Ks = 0;
			bool at_least_one_Ks_passed_dz = false;

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
			
			// loop over V0 K0s with matching of V0 pions to HPS pions
			for(unsigned k0_i = 0; k0_i < V0Ks->size(); k0_i++) 
			{
				// in dz Kaons selected as taus: https://twiki.cern.ch/twiki/bin/view/CMS/HiggsToTauTauWorking2016#Baseline_Selection
					if (abs((*V0Ks)[k0_i].vz() - pv_position.z()) > 0.2) continue; // from PV to Ks
					at_least_one_Ks_passed_dz = true;
					if (at_least_one_Ks_passed_dz) at_least_one_Ks_passed_dz = true;

				// init
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
						double E = 0, p_x = 0, p_y = 0, p_z = 0, inv_M = 0;
					
				// Fill basic hists
					//TEMP_COMMENTED h_Ks_v0_pions_dR->Fill((v0_ks_pion1).DeltaR(v0_ks_pion2));
					//TEMP_COMMENTED h_Ks_v0_PV_dXY->Fill(distMagXY_PV_Ks); 
					//TEMP_COMMENTED h_Ks_v0_vx->Fill(v_Ks.X());
					//TEMP_COMMENTED h_Ks_v0_vy->Fill(v_Ks.Y());
					//TEMP_COMMENTED h_Ks_v0_vz->Fill(v_Ks.Z());
					//TEMP_COMMENTED h_Ks_v0_dx->Fill(abs(v_Ks.X() - pv_position.x()));
					//TEMP_COMMENTED h_Ks_v0_dy->Fill(abs(v_Ks.Y() - pv_position.y()));
					//TEMP_COMMENTED h_Ks_v0_dz->Fill(abs(v_Ks.Z() - pv_position.z()));

				// Basic output
					dout("distMagXY_PV_Ks:", distMagXY_PV_Ks);
					dlog("\tKs", k0_i, "(", (*V0Ks)[k0_i].charge(), ")", "vertex:",  (*V0Ks)[k0_i].vx(), (*V0Ks)[k0_i].vy(),  (*V0Ks)[k0_i].vz());
					dlog("\t\tp:",  (*V0Ks)[k0_i].px(), (*V0Ks)[k0_i].py(),  (*V0Ks)[k0_i].pz(), "Energy:", (*V0Ks)[k0_i].energy());
					dlog("\t\tdaughters tot number:", num, ";", " moth number:", num_moth);
					
					// inv mass of K: Loop over daughters of V0s K - which are basically 2 pions
						for( int j = 0; j < num; j++) // always two
						{
							const reco::Candidate* daug_pion = (*V0Ks)[k0_i].daughter(j); //RecO_Cand_type(daug_pion); is unknown
							TLorentzVector vec_pion(daug_pion->px(),daug_pion->py(),daug_pion->pz(),daug_pion->energy());
								dlog("\t\t\t\tpion", j, "(", daug_pion->charge(), ")", " mass:", vec_pion.M(), "p:", daug_pion->px(),daug_pion->py(), daug_pion->pz(), "Energy:", daug_pion->energy());

							if (daug_pion->pt() > max_ks_daughter_pt) max_ks_daughter_pt = daug_pion->pt();
							//TEMP_COMMENTED h_Ks_v0_daughter_pt->Fill(daug_pion->pt());

							E += daug_pion->energy();
							p_x += daug_pion->px();
							p_y += daug_pion->py();
							p_z += daug_pion->pz();
						}
						inv_M += sqrt(pow(E, 2) - pow(p_x, 2) - pow(p_y, 2) - pow(p_z, 2));
						Ks_v0_inv_m_pi = inv_M;
						pt_1 = v0_ks_pion1.Pt();
						pt_2 = v0_ks_pion2.Pt();
						eta_1 = v0_ks_pion1.Eta();
						eta_2 = v0_ks_pion2.Eta();

						tree->Fill();
						v_Ks_v0_inv_m_pi.push_back(inv_M);
						v_pt_1.push_back(v0_ks_pion1.Pt());
						v_pt_2.push_back(v0_ks_pion2.Pt());
						v_eta_1.push_back(v0_ks_pion1.Eta());
						v_eta_2.push_back(v0_ks_pion2.Eta());
						//TEMP_COMMENTED h_Ks_v0_inv_m_pi->Fill(inv_M);
					
				// Matching to HPS pions pi+-
					if (match_KsV0_to_HPS && PF_hps_taus.isValid() )
					{
						dout("\t\t\tPF taus:", PF_hps_taus->size());

						for( unsigned tau_index = 0; tau_index < PF_hps_taus->size(); tau_index++)
						{
							// init
								reco::PFTauRef pftauref(PF_hps_taus, tau_index);
								TLorentzVector v_tau(pftauref->px(), pftauref->py(), pftauref->pz(), pftauref->energy());
							
							dout("\t\t\t\tPF tau i", tau_index, "p:", pftauref->px(), pftauref->py(), pftauref->pz(), "Energy:", pftauref->energy(), "dR:", v_tau.DeltaR(v_Ks));

							// Check if tau jet kinematicaly matches to Ks of V0: tau^ks dr < 0.5
							// Ks are selected close to PV, so dR to tau jet ensures that the Ks of V0 is in tau jet
								if (v_tau.DeltaR(v_Ks) < 0.5) 
								{
									number_of_passed_ks_in_hps_tau_hard_cut++;
									//TEMP_COMMENTED h_Ks_v0_found_in_hps_tau_dRcut->Fill(v_tau.DeltaR(v_Ks));
								}
								else continue;

							// Create pions collection of tau
								//pi+-
									vector < reco::PFCandidatePtr  > tau_signalPFChargedHadrCands    = pftauref->signalPFChargedHadrCands();//typedef edm::Ptr<PFCandidate> reco::PFCandidatePtr
									vector < reco::PFCandidatePtr  > tau_isolationPFChargedHadrCands = pftauref->isolationPFChargedHadrCands(); // //typedef edm::Ptr<PFCandidate> reco::PFCandidatePtr
								// All pi+-'s
								vector < reco::PFCandidatePtr > tau_picharge = pftauref->signalPFChargedHadrCands(); // vector < edm::Ptr<PFCandidate> >
																								tau_picharge.insert(tau_picharge.end(), tau_isolationPFChargedHadrCands.begin(), tau_isolationPFChargedHadrCands.end());
							
							//TEMP_COMMENTED h_Ks_v0_n_pion_in_tau_jets_with_good_Ks->Fill(tau_picharge.size());
							
							/*
								// //dR of all tau pions pairs
								//   for(unsigned pi_i = 0; pi_i < tau_picharge.size() - 1 && tau_picharge.size() > 1; pi_i++)
								//   {
								//     TLorentzVector tau_pion(tau_picharge[pi_i]->px(), tau_picharge[pi_i]->py(), tau_picharge[pi_i]->pz(), tau_picharge[pi_i]->energy());
								//     for(unsigned pi_i2 = pi_i + 1; pi_i2 < tau_picharge.size() ; pi_i2++)
								//     {
								//       TLorentzVector tau_pion2(tau_picharge[pi_i2]->px(), tau_picharge[pi_i2]->py(), tau_picharge[pi_i2]->pz(), tau_picharge[pi_i2]->energy());
								//       h_Ks_v0_hps_pions_combinatoric_dR->Fill((tau_pion).DeltaR(tau_pion2));

								//       h_tau_comb_pions_m_inv->Fill((tau_pion + tau_pion2).M());
								//     }
								//   }
							*/

							// Pion to pion matching
								// init
									bool both_ks_pions_found_in_tau = true;
									int firstfound = -1, secondfound = -1;
									double deltaR1 = 100, deltaR2 = 100;
								
								for(unsigned pi_i = 0; pi_i < tau_picharge.size(); pi_i++)
								{
									//init
										TLorentzVector tau_pion(tau_picharge[pi_i]->px(), tau_picharge[pi_i]->py(), tau_picharge[pi_i]->pz(), tau_picharge[pi_i]->energy());
										//TEMP_COMMENTED h_Ks_v0_pions_and_hps_pions_combined_dR->Fill(tau_pion.DeltaR(v0_ks_pion1));
										//TEMP_COMMENTED h_Ks_v0_pions_and_hps_pions_combined_dR->Fill(tau_pion.DeltaR(v0_ks_pion2));
										//if (tau_pion.M() < 0.1) break; //ISSUE

									// Best possible matching of tau jets pions to Ks V0 pions
									// if matched to first pion of Ks V0 and it is matched better than the previous candidate
									if (tau_picharge[pi_i]->charge() == (*V0Ks)[k0_i].daughter(0)->charge() &&
											tau_pion.DeltaR(v0_ks_pion1) < 0.5 &&
											(firstfound == -1 || deltaR1 > tau_pion.DeltaR(v0_ks_pion1)))
									{
										deltaR1 = tau_pion.DeltaR(v0_ks_pion1);
										firstfound = pi_i;
										dlog("dR:", v_tau.DeltaR(v_Ks), "But first pion found", "\ntau pion1:", tau_picharge[pi_i]->px(), tau_picharge[pi_i]->py(), tau_picharge[pi_i]->pz(), tau_picharge[pi_i]->energy());

										// check if old first pion fits better to second pion - NO NEED THE CHARGE IS DIFFERENT
										/*
											if (firstfound != -1 &&
													tau_pion.DeltaR(v0_ks_pion2) < 0.5 &&
													tau_picharge[firstfound]->charge() == (*V0Ks)[k0_i].daughter(1)->charge() &&
													(secondfound == -1 ||
														deltaR2 > TLorentzVector(tau_picharge[firstfound]->px(), tau_picharge[firstfound]->py(), tau_picharge[firstfound]->pz(), tau_picharge[firstfound]->energy()).DeltaR(v0_ks_pion2)))
											{
												deltaR2 = TLorentzVector(tau_picharge[firstfound]->px(), tau_picharge[firstfound]->py(), tau_picharge[firstfound]->pz(), tau_picharge[firstfound]->energy()).DeltaR(v0_ks_pion2);
												secondfound = firstfound;
											}
										*/
										
										// find in the rest of the pions the one that is matching the second pion of Ks V0
										/*
											both_ks_pions_found_in_tau = false;
											for(unsigned pi_i2 = pi_i + 1 ; pi_i2 < tau_picharge.size() ; pi_i2++)
											{
												TLorentzVector tau_pion2(tau_picharge[pi_i2]->px(), tau_picharge[pi_i2]->py(), tau_picharge[pi_i2]->pz(), tau_picharge[pi_i2]->energy());
												if (tau_pion2.DeltaR(v0_ks_pion2) < 0.1 &&
														tau_picharge[pi_i2]->charge() == (*V0Ks)[k0_i].daughter(1)->charge())
												{
													dlog("AND THE SECOND FOUND", pi_i2, tau_picharge[pi_i2]->px(),tau_picharge[pi_i2]->py(),tau_picharge[pi_i2]->pz(),tau_picharge[pi_i2]->energy());
													both_ks_pions_found_in_tau = true;
													h_Ks_v0_found_in_hps_tau_dR->Fill(v_tau.DeltaR(v_Ks));
													h_Ks_v0_found_in_hps_tau_m_inv->Fill((tau_pion + tau_pion2).M());
													number_of_passed_ks_in_hps_tau++;
													//exit(1);
													if (v_tau.DeltaR(v_Ks) > 0.5)
													{
														dout("CRITICAL KINEMATICAL ERROR");
														exit(1);
													}
													pi_i = tau_picharge.size() + 1; break;
												}
											}
										*/  
									}
									else if (tau_picharge[pi_i]->charge() == (*V0Ks)[k0_i].daughter(1)->charge() &&
													 tau_pion.DeltaR(v0_ks_pion2) < 0.5 &&
													 (secondfound == -1 || deltaR2 > tau_pion.DeltaR(v0_ks_pion2)))
									{
										deltaR2 = tau_pion.DeltaR(v0_ks_pion2);
										secondfound = pi_i;
										/*
											// dlog("dR:",v_tau.DeltaR(v_Ks),"But first pion found","\ntau pion1:", tau_picharge[pi_i]->px(), tau_picharge[pi_i]->py(), tau_picharge[pi_i]->pz(), tau_picharge[pi_i]->energy());
											// both_ks_pions_found_in_tau = false;
											// for(unsigned pi_i2 = pi_i + 1 ; pi_i2 < tau_picharge.size() ; pi_i2++)
											// {
											//   TLorentzVector tau_pion2(tau_picharge[pi_i2]->px(), tau_picharge[pi_i2]->py(), tau_picharge[pi_i2]->pz(), tau_picharge[pi_i2]->energy());
											//   if( tau_pion2.DeltaR(v0_ks_pion1) < 0.1 && tau_picharge[pi_i2]->charge() == (*V0Ks)[k0_i].daughter(0)->charge())
											//   {
											//     dlog("AND THE SECOND FOUND", pi_i2, tau_picharge[pi_i2]->px(),tau_picharge[pi_i2]->py(),tau_picharge[pi_i2]->pz(),tau_picharge[pi_i2]->energy());
											//     both_ks_pions_found_in_tau = true;
											//     h_Ks_v0_found_in_hps_tau_dR->Fill(v_tau.DeltaR(v_Ks));
											//     h_Ks_v0_found_in_hps_tau_m_inv->Fill((tau_pion + tau_pion2).M());
											//     number_of_passed_ks_in_hps_tau++;
											//     //exit(1);
											//     if (v_tau.DeltaR(v_Ks) > 0.5)
											//     {
											//       dout("CRITICAL KINEMATICAL ERROR");exit(1);
											//     }
											//     pi_i = tau_picharge.size() + 1; break;
											//   }
											// }
										*/
									}
								}

							// Invariant mass of two matched pions
							if (firstfound != -1 && secondfound != -1)
							{
								// init
									TLorentzVector v_tau_picharge_firstfound  = TLorentzVector(tau_picharge[firstfound]->px(), tau_picharge[firstfound]->py(), tau_picharge[firstfound]->pz(), tau_picharge[firstfound]->energy());
									TLorentzVector v_tau_picharge_secondfound = TLorentzVector(tau_picharge[secondfound]->px(), tau_picharge[secondfound]->py(), tau_picharge[secondfound]->pz(), tau_picharge[secondfound]->energy());
									TLorentzVector v_Ks_from_hps_pions = v_tau_picharge_firstfound + v_tau_picharge_secondfound;
									double inv_m_pions = v_Ks_from_hps_pions.M();
									//if (abs(inv_m_pions - 0.498) > 0.07) continue; // cut on the inv m of two matched hps pions 
									// dout("Is it different:", (TLorentzVector(tau_picharge[firstfound]->px(), tau_picharge[firstfound]->py(), tau_picharge[firstfound]->pz(), tau_picharge[firstfound]->energy()) + TLorentzVector(tau_picharge[secondfound]->px(), tau_picharge[secondfound]->py(), tau_picharge[secondfound]->pz(), tau_picharge[secondfound]->energy())).M(), "and", inv_m_pions);
									number_of_passed_ks_in_hps_tau++;


								// Significance to BS or PV
								/*
									dout("Signif of Ks from HPS to BS");
										reco::TransientTrack* posTransTkPtr = 0;//nullptr;
										reco::TransientTrack* negTransTkPtr = 0;//nullptr;
										dout("theMagneticField", theMagneticField);
										//dout("tau_picharge[firstfound]->trackRef()", & tau_picharge[firstfound]->trackRef()); problems on rvalue http://stackoverflow.com/questions/16481490/error-taking-address-of-temporary-fpermissive
										dout("tau_picharge[firstfound]->bestTrack()", tau_picharge[firstfound]->bestTrack());
										reco::TransientTrack tmpTransient1( *tau_picharge[firstfound]->bestTrack(), theMagneticField);
										reco::TransientTrack tmpTransient2( *tau_picharge[secondfound]->bestTrack(), theMagneticField);
										if (tau_picharge[firstfound]->charge() < 0. && tau_picharge[secondfound]->charge() > 0.)
										{
											dout("first <0");
											 negTransTkPtr = &tmpTransient1;
											 posTransTkPtr = &tmpTransient2;
										} else if (tau_picharge[firstfound]->charge() > 0. && tau_picharge[secondfound]->charge() < 0.)
										{
											dout("first >0");
											 negTransTkPtr = &tmpTransient2;//tau_picharge[secondfound]->trackRef();
											 posTransTkPtr = &tmpTransient1;//tau_picharge[firstfound]->trackRef();
										} else
										{
											dout("tau_picharge[firstfound]->charge()", tau_picharge[firstfound]->charge());
										}
									dout("Fill the vector of TransientTracks to send to KVF");
										std::vector<reco::TransientTrack> transTracks;
										transTracks.reserve(2);
										if (negTransTkPtr ==0 || posTransTkPtr ==0) dout("one of the pointers is zero", negTransTkPtr, posTransTkPtr);
										transTracks.push_back(*posTransTkPtr);
										transTracks.push_back(*negTransTkPtr);
									dout("2D decay significance");
									// create the vertex fitter object and vertex the tracks
										 TransientVertex theRecoVertex;
											KalmanVertexFitter theKalmanFitter(true);
											theRecoVertex = theKalmanFitter.vertex(transTracks);
										 reco::Vertex theVtx = theRecoVertex;
										 dout("totalCov");
										SMatrixSym3D totalCov = theBeamSpot->rotatedCovariance3D() + theVtx.covariance();
										//SMatrixSym3D totalCov = BSposition.covariance() + theVtx.covariance();
									dout("distVecXY");
										SVector3 distVecXY(v_Ks_from_hps_pions.X() - BSposition.x(), v_Ks_from_hps_pions.Y() - BSposition.y(), 0.);
										double distMagXY = ROOT::Math::Mag(distVecXY);
										double sigmaDistMagXY = sqrt(ROOT::Math::Similarity(totalCov, distVecXY)) / distMagXY;
										//if (distMagXY/sigmaDistMagXY < vtxDecaySigXYCut_) continue;
									if (sigmaDistMagXY!=0)
									{
										dout("distMagXY/sigmaDistMagXY", distMagXY/sigmaDistMagXY);
										h_Ks_v0_found_in_hps_tau_significance->Fill(distMagXY/sigmaDistMagXY);
									}
									else
									{
										dout("division by zero");
										exit(1);
									}
								*/

								// IP significance
								//const reco::Track* tmpTrack = &(*iTk);

								// fill histograms
									//TEMP_COMMENTED h_Ks_v0_PV_dXY->Fill(distMagXY_PV_Ks);
									//TEMP_COMMENTED h_Ks_v0_BS_dXY->Fill(distMagXY_BS_Ks);
									//TEMP_COMMENTED h_Ks_v0_found_in_hps_tau_m_inv->Fill(inv_m_pions);
									// dR between matched pions of tau and Ks V0
									//TEMP_COMMENTED h_Ks_v0_found_in_hps_tau_dR->Fill(TLorentzVector(tau_picharge[firstfound]->px(), tau_picharge[firstfound]->py(), tau_picharge[firstfound]->pz(), tau_picharge[firstfound]->energy()).DeltaR(v_Ks));//   h_Ks_v0_found_in_hps_tau_dR->Fill(v_tau_picharge_firstfound.DeltaR(v_Ks));
									//TEMP_COMMENTED h_Ks_v0_found_in_hps_tau_dR->Fill(TLorentzVector(tau_picharge[secondfound]->px(), tau_picharge[secondfound]->py(), tau_picharge[secondfound]->pz(), tau_picharge[secondfound]->energy()).DeltaR(v_Ks));//   h_Ks_v0_found_in_hps_tau_dR->Fill(v_tau_picharge_secondfound.DeltaR(v_Ks));
							}
							//TEMP_COMMENTED else if (firstfound == -1 || secondfound == -1) h_Ks_v0_found_in_hps_tau_dR_only_one_pion_left->Fill(v_tau.DeltaR(v_Ks));
						}
					}
			}
			// fill histograms
				//TEMP_COMMENTED h_Ks_v0_n_ev_passing_dz_cut->Fill(1 * at_least_one_Ks_passed_dz);
				//TEMP_COMMENTED h_Ks_v0_n_Ks_in_jets_per_event->Fill(number_of_passed_ks_in_hps_tau_hard_cut);
				//TEMP_COMMENTED h_Ks_v0_count->Fill(v0_count); //h_v0_count->Fill(v0_count);
				//TEMP_COMMENTED h_Ks_v0_number_per_event->Fill(number_of_passed_ks);
				// ks from V0 coll, NPE that passed PV and found in HPS tau jets
				//TEMP_COMMENTED if (number_of_passed_ks_in_hps_tau > 0) h_Ks_v0_found_in_hps_tau->Fill(number_of_passed_ks_in_hps_tau);

			dlog("v0_count", v0_count,"number_of_passed_ks", number_of_passed_ks, "number_of_passed_ks_in_hps_tau:", number_of_passed_ks_in_hps_tau, "number_of_passed_ks_in_hps_tau_hard_cut", number_of_passed_ks_in_hps_tau_hard_cut);
		}

		// count the number of the Ks from standart collection in the dz < 0.2 cm from the PV 
			if (true)
			{  
				int N_standart_KS = 0;

				if (V0Ks_standart->size() > 0)
					for(unsigned k0_i = 0; k0_i < V0Ks_standart->size(); k0_i++) // over V0 K0s from the PV
					{
						if (abs((*V0Ks_standart)[k0_i].vz() - pv_position.z()) > 0.2) continue; // 
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
	// else if (!V0Ks.isValid()) dlog("UNVALID Ks's");
}