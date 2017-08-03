// -*- C++ -*-
//
// Package:    AOD_pi0_study/AOD_pi0
// Class:      AOD_pi0
//
/**\class AOD_pi0 AOD_pi0.cc AOD_pi0_study/AOD_pi0/plugins/AOD_pi0.cc

 Description: build a ks resonan

*/
#include "AOD_pi0.h"

// ------------ method called for each event  ------------
void AOD_pi0::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	v0_count = 0 ;
	pizero_count = 0;
	pt_1 = 0;
	pt_2 = 0;
	eta_1 = 0;
	eta_2 = 0;
	v_Ks_v0_inv_m_pi.clear();
	v_v0_count.clear();
	v_pt_1.clear();
	v_pt_2.clear();
	v_eta_1.clear();
	v_eta_2.clear();
	Ks_v0_inv_m_pi = -1;

	iEvent.getByToken( KshortCollectionToken_, V0Ks); //V0 Ks's token
	iEvent.getByToken( KshortCollectionTag_stand_, V0Ks_standart); //V0 Ks's standart
	iEvent.getByToken( TauPiZeroCollectionToken_, Strips); //actually HPS pi0s
	if (!IsData) iEvent.getByToken( GenParticleCollectionToken_, GenPart);
	iEvent.getByToken(BeamSpotToken_, TheBeamSpotHandle);
	//Tau's token based hps- according to https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookPATDataFormats#PatTau contains the same as PFjets
	//typedef vector< PFTau >  PFTauCollection - the collection of charged pions will be taken from here
	iEvent.getByToken( TauHPSCollectionToken_, PF_hps_taus);
	// Magnetic field
		// edm::ESHandle<MagneticField> theMagneticFieldHandle;
		// iSetup.get<IdealMagneticFieldRecord>().get(theMagneticFieldHandle);
		// const MagneticField* theMagneticField = theMagneticFieldHandle.product();

		//  edm::Handle<reco::TrackCollection> HPSTraHandle;
		// iEvent.getByToken(HPSTrackTagToken_, HPSTraHandle);
		// if (HPSTraHandle.isValid())  dout("HPStrackCollection is fine");
		// else dout("HPStrackCollection is NOT VALID", HPSTraHandle->size());
		// // const reco::TrackCollection* theTrackCollection = HPSTraHandle.product();

	const reco::BeamSpot* theBeamSpot = TheBeamSpotHandle.product();
	math::XYZPoint BSposition(theBeamSpot->position());

	/// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	/// PV -> pv_position
	if (crecprimvertex)
	{
		iEvent.getByToken(PVToken_, Vertex);
		primvertex_count = goodprimvertex_count = 0;

		pv_position = math::XYZPoint(0., 0., 0.);
		if (Vertex.isValid())
		{
			for(unsigned i = 0; i < Vertex->size(); i++)
			{
				primvertex_count++;
				if (i == 0) // pv_position - the first PV
				{
					primvertex_x = (*Vertex)[i].x();
					primvertex_y = (*Vertex)[i].y();
					primvertex_z = (*Vertex)[i].z();
					primvertex_chi2 = (*Vertex)[i].chi2();
					primvertex_ndof = (*Vertex)[i].ndof();
					primvertex_ntracks = (*Vertex)[i].tracksSize();
					primvertex_cov[0] = (*Vertex)[i].covariance(0,0); // xError()
					primvertex_cov[1] = (*Vertex)[i].covariance(0,1);
					primvertex_cov[2] = (*Vertex)[i].covariance(0,2);
					primvertex_cov[3] = (*Vertex)[i].covariance(1,1); // yError()
					primvertex_cov[4] = (*Vertex)[i].covariance(1,2);
					primvertex_cov[5] = (*Vertex)[i].covariance(2,2); // zError()

					Float_t ptq = 0.;
					for(reco::Vertex::trackRef_iterator it = (*Vertex)[i].tracks_begin(); it != (*Vertex)[i].tracks_end(); ++it)
						ptq += (*it)->pt() * (*it)->pt();
					primvertex_ptq = ptq;

					pv_position = (*Vertex)[i].position();
					primvertex = (*Vertex)[i];

					//TEMP_COMMENTED h_primvertex_x->Fill(primvertex_x);
					//TEMP_COMMENTED h_primvertex_y->Fill(primvertex_y);
					//TEMP_COMMENTED h_primvertex_z->Fill(primvertex_z);
					//TEMP_COMMENTED h_primvertex_chi2->Fill(primvertex_chi2);
					//TEMP_COMMENTED h_primvertex_ndof->Fill(primvertex_ndof);
					//TEMP_COMMENTED h_primvertex_ptq->Fill(primvertex_ptq);
					//TEMP_COMMENTED h_primvertex_ntracks->Fill(primvertex_ntracks);
					//TEMP_COMMENTED h_primvertex_cov_x->Fill(primvertex_cov[0]);
					//TEMP_COMMENTED h_primvertex_cov_y->Fill(primvertex_cov[3]);
					//TEMP_COMMENTED h_primvertex_cov_z->Fill(primvertex_cov[5]);
				}

				if ((*Vertex)[i].isValid() &&
						!(*Vertex)[i].isFake() &&
						(*Vertex)[i].ndof() >= 4 &&
						(*Vertex)[i].z() > -24 &&
						(*Vertex)[i].z() < 24 &&
						(*Vertex)[i].position().Rho() < 2.) goodprimvertex_count++;
			}
		}
		//TEMP_COMMENTED h_primvertex_count->Fill(primvertex_count);
		//TEMP_COMMENTED h_goodprimvertex_count->Fill(goodprimvertex_count);
		dout("PrimVertex:", pv_position.X(), pv_position.Y(), pv_position.Z());
	}

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

	if (false) AddGammas(iEvent, iSetup);

	/// HPS Pi0's and taus - with loop among reco::tau
	if (false && Strips.isValid())
		{
			unsigned int matched_pi = 0;
			if (Strips->size() > 0)  dout("Number of HPS piz0's =", Strips->size());
			if (PF_hps_taus->size() > 0) dout("Number of RECO Tau =", PF_hps_taus->size());

			for (unsigned int i = 0; i < Strips->size(); i++)
			{
				if (debug)
				{
					dout("\tPi0_", i, "(", (*Strips)[i].charge(), ")", (*Strips)[i].px(), (*Strips)[i].py(), (*Strips)[i].pz(), (*Strips)[i].p());
					//dout(" \tPI0", i, "px =", (*Strips)[i].px(), "py =", (*Strips)[i].py(), "pz =", (*Strips)[i].pz(), "e =", (*Strips)[i].p(), "pt =", (*Strips)[i].pt(), "; vx =", (*Strips)[i].vx(),  "vy =", (*Strips)[i].vy(),  "vz =", (*Strips)[i].vz());
				}
				pizero_count++;
				if (pizero_count>=1000)
				{
					cerr << "number of pizeros > 1000. They are missing." << endl;
					break;
				}

				// RECO Taus
					bool foundpi = false;
					for (unsigned int i_tau = 0; i_tau < PF_hps_taus->size(); i_tau++) //over all taus
					{
						cout << "\t\tTau_" << i_tau << endl;
						reco::PFTauRef pftauref(PF_hps_taus, i_tau);

						//Signal Pi0
							const vector < reco::RecoTauPiZero > tau_pizeros = pftauref->signalPiZeroCandidates();
							dout("\t\t\tsignal tau_pizeros:", tau_pizeros.size());
							if (tau_pizeros.size() > 0)
								for (unsigned int j_pi = 0; j_pi < tau_pizeros.size(); j_pi++) //over signal pions in tau
								{
									if (pow(tau_pizeros[j_pi].charge() - (*Strips)[i].charge(), 2) +
											pow(tau_pizeros[j_pi].px() - (*Strips)[i].px(), 2) +
											pow(tau_pizeros[j_pi].py() - (*Strips)[i].py(), 2) +
											pow(tau_pizeros[j_pi].pz() - (*Strips)[i].pz(), 2) +
											pow(tau_pizeros[j_pi].energy() - (*Strips)[i].energy(), 2) < 0.001)
									{
										foundpi = true;
										dout("\t\t\t\tPi0_", j_pi, "(", tau_pizeros[j_pi].charge(),")", tau_pizeros[j_pi].px(), tau_pizeros[j_pi].py(), tau_pizeros[j_pi].pz(), tau_pizeros[j_pi].energy());
										break;
									}
								}

						//Isol Pi0
							const vector < reco::RecoTauPiZero > tau_pizeros_isol = pftauref->isolationPiZeroCandidates();
							dout("\t\t\tisolation tau_pizeros_isol :", tau_pizeros_isol.size());
							if (tau_pizeros_isol.size() > 0 && !foundpi)
								for (unsigned int j_pi = 0; j_pi < tau_pizeros_isol.size(); j_pi++) //over isolation pions in tau
								{
									if (pow(tau_pizeros_isol[j_pi].charge() - (*Strips)[i].charge(), 2) +
											pow(tau_pizeros_isol[j_pi].px() - (*Strips)[i].px(), 2) +
											pow(tau_pizeros_isol[j_pi].py() - (*Strips)[i].py(), 2) +
											pow(tau_pizeros_isol[j_pi].pz() - (*Strips)[i].pz(), 2) +
											pow(tau_pizeros_isol[j_pi].energy() - (*Strips)[i].energy(), 2) < 0.001)
									{
										foundpi = true;
										dout("\t\t\t\tPi0_", j_pi, "(", tau_pizeros_isol[j_pi].charge(), ")", tau_pizeros_isol[j_pi].px(), tau_pizeros_isol[j_pi].py(), tau_pizeros_isol[j_pi].pz(), tau_pizeros_isol[j_pi].energy());
										break;
									}
								}

							if (foundpi)
							{
								matched_pi++;
								break;
							}
					}
					if (!foundpi)
					{
						cerr << "PION LOST" << endl;
					}
			}
			if (matched_pi != Strips->size())
			{
				cerr << "===>THE NUMBER OF PIONS DON'T MATCH" << endl;
				num_ev_tau_pi_not_in_hps_pi++;
			}
			else dout("ALL PIONS MATCHED");
		}

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

	/// HPS reco Taus only - general study
	if (false && PF_hps_taus.isValid() )
		{
			if (PF_hps_taus->size() > 0)  dout("Number of Tau = ", PF_hps_taus->size());
			for (unsigned int i = 0; i < PF_hps_taus->size(); i++) // Over Tau's
			{
				reco::PFTauRef pftauref(PF_hps_taus, i); // one Tau instance, typedef edm::Ref<PFTauCollection> reco::PFTauRef
				dlog("tau #", i, "from", PF_hps_taus->size(), "(", pftauref->vx(), pftauref->vy(), pftauref->vz(), ")");
					//reco::PFTau pfTau(PF_hps_taus, i); //same as (*PF_hps_taus)[i]
					//edm::AtomicPtrCache< vector< reco::RecoTauPiZero > >
					//const vector< reco::RecoTauPiZero > a  = (*PF_hps_taus)[i].signalPiZeroCandidates();
					//(*reco::PFTau)pfTau->signalPiZeroCandidates();
					//const vector < RecoTauPiZero > &   isolationPiZeroCandidates () const
					//const vector< RecoTauPiZero > & reco::PFTau::isolationPiZeroCandidates (   ) const
					//const vector < RecoTauPiZero > &   signalPiZeroCandidates () const

				// Lists of objects of the tau
					// Pizeros list
						// Signal pi0's
							vector < reco::RecoTauPiZero > tau_pizeros_sig = pftauref->signalPiZeroCandidates();
						//Isolation pi0's
							const vector < reco::RecoTauPiZero > tau_pizeros_isol = pftauref->isolationPiZeroCandidates(); // pions of the considered Tau
							vector < reco::RecoTauPiZero *> point_tau_pizeros_isol;
						// All pi0's
							vector < reco::RecoTauPiZero > tau_pizeros = pftauref->signalPiZeroCandidates();
																						 tau_pizeros.insert(tau_pizeros.end(), tau_pizeros_isol.begin(), tau_pizeros_isol.end());
							vector < reco::RecoTauPiZero* > point_tau_pizeros;

					// Hadrons list - all are pions charged and neutral
						//pi+-
							vector < reco::PFCandidatePtr  > tau_signalPFChargedHadrCands    = pftauref->signalPFChargedHadrCands();//typedef edm::Ptr<PFCandidate> reco::PFCandidatePtr
							vector < reco::PFCandidatePtr  > tau_isolationPFChargedHadrCands = pftauref->isolationPFChargedHadrCands(); // //typedef edm::Ptr<PFCandidate> reco::PFCandidatePtr
							// All pi+-'s
							vector < reco::PFCandidatePtr > tau_picharge = pftauref->signalPFChargedHadrCands(); // vector < edm::Ptr<PFCandidate> >
																							tau_picharge.insert(tau_picharge.end(), tau_isolationPFChargedHadrCands.begin(), tau_isolationPFChargedHadrCands.end());
							std::vector<reco::PFCandidatePtr * > point_tau_picharge;
						//pi0
							vector < reco::PFCandidatePtr  > tau_signalPFNeutrHadrCands    = pftauref->signalPFNeutrHadrCands(); // //typedef edm::Ptr<PFCandidate> reco::PFCandidatePtr
							vector < reco::PFCandidatePtr  > tau_isolationPFNeutrHadrCands = pftauref->isolationPFNeutrHadrCands(); // //typedef edm::Ptr<PFCandidate> reco::PFCandidatePtr
							// All pi0_had's
							vector < reco::PFCandidatePtr > tau_pizeros_had = pftauref->signalPFNeutrHadrCands(); // //typedef edm::Ptr<PFCandidate> reco::PFCandidatePtr
																							tau_pizeros_had.insert(tau_pizeros_had.end(), tau_isolationPFNeutrHadrCands.begin(), tau_isolationPFNeutrHadrCands.end());
							std::vector<reco::PFCandidatePtr *> point_tau_pizeros_had;

					// NOT USED list
						vector < reco::PFRecoTauChargedHadron > tau_signalTauChargedHadronCandidates = pftauref->signalTauChargedHadronCandidates(); // gives 0 size
						vector < reco::PFRecoTauChargedHadron > tau_isolationTauChargedHadronCandidates = pftauref->isolationTauChargedHadronCandidates();// gives 0 size
						const reco::PFCandidatePtr tau_leadPFChargedHadrCand = pftauref->leadPFChargedHadrCand(); // by output eeror message // //typedef edm::Ptr<PFCandidate> reco::PFCandidatePtr
						float tau_leadPFChargedHadrCandsignedSipt = pftauref->leadPFChargedHadrCandsignedSipt();
						reco::PFRecoTauChargedHadronRef tau_leadTauChargedHadronCandidate = pftauref->leadTauChargedHadronCandidate(); // can not implement

				// Pions of PFRecoTau - all RecoTauPiZero
				dlog("\t...................");
				dlog("\tPions of RecoTauPiZero of PFRecoTau");
					// Signal pi0's
					if (tau_pizeros_sig.size() && true)
					{
						dlog("\t\t signal tau pi0's number:", tau_pizeros_sig.size());
						//TEMP_COMMENTED num_pions->Fill(tau_pizeros_sig.size());
						double inv_M = 0;
						if (tau_pizeros_sig.size() > 0)
						{
							TLorentzVector first(0, 0, 0, 0);
							for (unsigned int j = 0; j < tau_pizeros_sig.size(); j++)
							{
								dout("\t\t\tPi0_", j, "(", tau_pizeros_sig[j].charge(), ")", tau_pizeros_sig[j].px(), tau_pizeros_sig[j].py(), tau_pizeros_sig[j].pz(), tau_pizeros_sig[j].energy());
								TLorentzVector second(tau_pizeros_sig[j].px(), tau_pizeros_sig[j].py(), tau_pizeros_sig[j].pz(), tau_pizeros_sig[j].energy());
								first = (first + second);

								//dlog("\t\t\t signal tau pi0_", j, " vtx:", tau_pizeros_sig[j].vx(), tau_pizeros_sig[j].vy(), tau_pizeros_sig[j].vz(), " || p:", tau_pizeros_sig[j].px(), tau_pizeros_sig[j].py(), tau_pizeros_sig[j].pz());
							}
							inv_M = first.M();
						}
						dout("\tTAU_", i, " WITH TOTAL INVARIANT MASS OF pi0's:", inv_M);
					}
					else dlog("\t\t\tto few pions in signal pions ");

					//Isolation pi0's
					point_tau_pizeros_isol = TransformToPointers(tau_pizeros_isol, point_tau_pizeros_isol);
					CombinatoricOfTwoToNeutralInvM(point_tau_pizeros_isol, "isolat tau pi0", "tau", "pi0", taus_isol_pi0_inv_m_to_ks, taus_isol_pi0_inv_pt);

					//All pi0's loop
					point_tau_pizeros = TransformToPointers(tau_pizeros, point_tau_pizeros);
					CombinatoricOfTwoToNeutralInvM(point_tau_pizeros, "all tau pi0", "tau", "pi0", taus_pi0_inv_m_to_ks, taus_pi0_inv_pt);

				//Charged hadrons of PFRecoTau - all PFCandidatePtr
				dlog("\t...................");
				dlog("\tCharged hadrons of PFRecoTau ");
				point_tau_picharge = TransformToPointers(tau_picharge, point_tau_picharge);
				CombinatoricOfTwoToNeutralInvM(tau_picharge, "all tau pi+-", "tau", "pi+-", taus_pi_charged_inv_m_to_ks, taus_pi_charged_inv_pt);

				//Neutral hadrons of PFRecoTau - all PFCandidatePtr
				dlog("\t...................");
				dlog("\tNeutral hadrons of PFRecoTau "); // this is empty
				point_tau_pizeros_had = TransformToPointers(tau_pizeros_had, point_tau_pizeros_had);
				CombinatoricOfTwoToNeutralInvM(tau_pizeros_had, "all tau pi0_had", "tau", "pi0_had", taus_pi0_had_inv_m_to_ks, taus_pi0_had_inv_pt);
			}
		}
		// else if (!PF_hps_taus.isValid()) dout("no valid PF_hps_taus");

	/// GEN Particles K(892)- gen level study
	if (false && !IsData && GenPart.isValid())
		{
			int gen_count = GenPart->size();
			dlog("Size:", gen_count);
			cout << string(10, '=') << "Size:" << gen_count << endl;
			int NumberOfK892_neutral = 0,  NumberOfK892_charged = 0;
				for(unsigned i = 0 ; i < GenPart->size() ; i++)
				{
					const reco::GenParticle & gen_prt = (*GenPart)[i];
					if (gen_prt.numberOfMothers() != 0) continue;
					//cout << "paritcle # " << i << " (" << (*GenPart)[i].charge() << ") " << "("<< gen_prt.numberOfMothers()<<":"<< gen_prt.numberOfDaughters()<< ")"<< endl;
					//dlog("paritcle #", i, "(", (*GenPart)[i].charge(), ")", "(", gen_prt.numberOfMothers(),":", gen_prt.numberOfDaughters(), ") pdgId:", gen_prt.pdgId() << " status:" <<  gen_prt.status());
					if (abs(gen_prt.pdgId()) == 313 )
						{
							NumberOfK892_neutral++;
							//TEMP_COMMENTED map_K892_gen_histos[313][1]->Fill(gen_prt.vx());
							//TEMP_COMMENTED map_K892_gen_histos[313][2]->Fill(gen_prt.vy());
							//TEMP_COMMENTED map_K892_gen_histos[313][3]->Fill(gen_prt.vz());
						}
					else NumberOfK892_neutral += CountGenPrtID(&gen_prt, 313, 0, 2, &map_K892_gen_histos[313]);

					if (abs(gen_prt.pdgId()) == 323)
						{
							NumberOfK892_charged++;
							//TEMP_COMMENTED map_K892_gen_histos[323][1]->Fill(gen_prt.vx());
							//TEMP_COMMENTED map_K892_gen_histos[323][2]->Fill(gen_prt.vy());
							//TEMP_COMMENTED map_K892_gen_histos[323][3]->Fill(gen_prt.vz());
						}
					else NumberOfK892_charged += CountGenPrtID(&gen_prt, 323, 0, 2, &map_K892_gen_histos[323]);

				}
				dout("NumberOfK892 neutral:",NumberOfK892_neutral);
				//TEMP_COMMENTED map_K892_gen_histos[313][0]->Fill(NumberOfK892_neutral);
				//h_K892_0_gen_number_per_event->Fill(NumberOfK892_neutral);
				dout("NumberOfK892 charged:", NumberOfK892_charged);
				//TEMP_COMMENTED map_K892_gen_histos[323][0]->Fill(NumberOfK892_charged);
		}

	/// GEN Particles - gen level study
	if (false && !IsData && GenPart.isValid())
		{
			dlog("GEN Particles");
			//vector<reco::GenParticleCollection> gen_daughters;
			int gen_count = GenPart->size();
			dlog("Size:", gen_count);
			if (outputGenEvolution)
			{
				ofs << string(10, '=') << "Size:" << gen_count << endl;
				for(unsigned i = 0 ; i < GenPart->size() ; i++)
				{
					const reco::GenParticle & gen_prt = (*GenPart)[i];
					ofs << "paritcle # " << i << " (" << (*GenPart)[i].charge() << ") " << endl;//dlog("paritcle #", i, "(", (*GenPart)[i].charge(), ")");
					ofs << "\t" << gen_prt.pdgId() << " " <<  gen_prt.status() << endl;//dlog("\t", gen_prt.pdgId(), gen_prt.status());
					GenEvolution(&gen_prt, 2);
				}
			}

			v_daughters_k0s_to_pi0.clear();
			BuildTo(GenPart, h_gen_k0_all_to_pi0, v_daughters_k0s_to_pi0, 310, map_kaons, 2, 111, map_pions);// k0s 310
			if (false)//temp
			{
				v_daughters_k0l_to_pi0.clear();
				BuildTo(GenPart, h_gen_k0_all_to_pi0, v_daughters_k0l_to_pi0, 130, map_kaons, 2, 111, map_pions);// k0l 130
				//v_daughters_k0_to_pi0.clear();
				//BuildTo(GenPart, h_gen_k0_all_to_pi0, v_daughters_k0_to_pi0, 311, map_kaons, 2, 111, map_pions);// k0 311 - is empty by definition of converting to 310 or 130
				v_daughters_k0s_to_pic.clear();
				BuildTo(GenPart, h_gen_k0_all_to_pic, v_daughters_k0s_to_pic, 310, map_kaons, 2, 211, map_pions);// k0s 310
				v_daughters_k0l_to_pic.clear();
				BuildTo(GenPart, h_gen_k0_all_to_pic, v_daughters_k0l_to_pic, 130, map_kaons, 2, 211, map_pions);// k0l 130
			}

			v_daughters_k0s_to_pi0_SimToStrips.clear();//vector< vector <const reco::Candidate *>>
			v_daughters_k0s_to_pi0_StipsToSim.clear();//vector< vector <const reco::Candidate *>>
			MakeVectorofRef(Strips, v_strips_ref);//typedef edm::Ref< RecoTauPiZeroCollection >   RecoTauPiZeroRef
			Match(v_daughters_k0s_to_pi0, v_strips_ref, v_daughters_k0s_to_pi0_SimToStrips, v_daughters_k0s_to_pi0_StipsToSim);//Gen to Reco, resulting vactor
			//reco::CandMatchMap map = MCTruthDeltaRMatcher("pions and pizeros")
		}
		// else if (!IsData) dlog("\tno GenPart.isValid()");
		// else dlog("no GenPart in data");

	once_tree->Fill();
}

template <typename Tc, typename Tr>// T - collection   reco::PFTauRef pftauref(PF_hps_taus, i);
void AOD_pi0::MakeVectorofRef(edm::Handle< Tc > Collection, vector< Tr* > v_of_ref) {}

template <typename T>
void AOD_pi0::Match(vector < vector <const reco::Candidate *>> & From,
										vector < T *> & To, //e.g. Strips
										vector < vector< vector <T *>>> & SimToReco,
										vector < vector <const reco::Candidate *>> & RecoToSim
										//multiset < pair< int*, vector <string*> > >
										)
{
	dlog("in the beg To.size()", To.size());
	if (RecoToSim.size() < To.size())
		for(unsigned i = 0; i < To.size(); i++)
			RecoToSim.push_back(vector <const reco::Candidate *> ());

	for(unsigned i = 0; i < From.size(); i++) //Loop on K
	{
		SimToReco.push_back(vector< vector <T *>> ());
		for(unsigned ii = 0; ii < From[i].size(); ii++ ) //Loop on Pions in K
		{
			SimToReco[SimToReco.size() - 1].push_back( vector <T *> ());
			bool found = false;
			dlog("To.size()", To.size());
			for(unsigned j = 0; j < To.size(); j++)
			{
				cout << "\tfrom Kaon num" << i << ", pion number" << ii;
				if (Incone(*To[j], *From[i][ii], 0.2))
				{
					vector< vector <T *>> * K = &(SimToReco[SimToReco.size() - 1]);
					vector <T *> *Pi = &((*K)[K->size() - 1]);
					//SimToReco[SimToReco.back() - 1][last].push_back(& To[j])
					Pi->push_back(To[j]);
					RecoToSim[j].push_back(From[i][ii]);
					cout << "matched to Strip number " << j << " ";
					found = true;
				}
			}

			if (!found) dlog("was not matched");
			else dlog();
		}
	}
}

template <typename Ta, typename Tb>
bool AOD_pi0::Incone(Ta& A, Tb& B, double cone_size)
{
	TLorentzVector a(A.px(), A.py(), A.pz(), A.energy());
	TLorentzVector b(B.px(), B.py(), B.pz(), B.energy());
	double dR = a.DeltaR(b);
	if (dR <= cone_size) return true;
	else return false;
}

void AOD_pi0::BuildTo(edm::Handle<std::vector<reco::GenParticle> >& genPart,
											TH1 *hist,
											vector< vector <const reco::Candidate *>>& daughters,
											int mom_pdgid, map <long, string>& map_mom,
											int num_daugh,
											int daugh_pdgid, map <long, string>& map_daugh)
{
	for(unsigned i = 0; i < genPart->size(); i++)
	{
		const reco::GenParticle & gen_prt = (*genPart)[i];
		BuildTo(&gen_prt, hist, daughters,  mom_pdgid, map_mom,  num_daugh, daugh_pdgid, map_daugh);
	}
	dlog("Found K", mom_pdgid, "to", daugh_pdgid, ":", daughters.size());
	//for(std::vector<vector <const reco::Candidate *>>::iterator it = daughters.begin(); it != daughters.end(); ++it)
	for(unsigned i = 0; i < daughters.size(); i++)
	{
		//dlog("\t\tFound pi", ":", daughters[i].size());
		TLorentzVector temp;//reco::Candidate::LorentzVector
		for(std::vector<const reco::Candidate *>::iterator jt = daughters[i].begin(); jt != daughters[i].end(); ++jt)
			temp += TLorentzVector((*jt)->px(), (*jt)->py(), (*jt)->pz(), (*jt)->energy());
		//dlog("temp:", temp.Px(), temp.Py(), temp.Pz(), temp.E(), temp.M());
		hist->Fill(temp.M());
	}
}

void AOD_pi0::BuildTo(const reco::Candidate *mom,
											TH1 *hist,
											vector< vector <const reco::Candidate *>> & daughters,
											int mom_pdgid, map <long, string>& map_mom,
											int num_daugh,
											int daugh_pdgid, map <long, string>& map_daugh )
{
	for(unsigned i = 0; i < mom->numberOfDaughters(); i++)
	{
		const reco::Candidate * prt = mom->daughter(i);
		if (abs(prt->pdgId()) == mom_pdgid && prt->numberOfDaughters() > 0 && NoRadiating(prt, map_mom))
		{
			dout("K found");
			daughters.push_back(vector <const reco::Candidate *>());
			FindFinalPrt(&daughters.back(), prt, num_daugh, daugh_pdgid);//(*daughters)[daughters.size() - 1]
			if ( (daughters.back()).size() == 0 )
			{
				dout("but no daughers which would fit");
				daughters.erase(daughters.end() - 1);
			}
		}
		else if (prt->numberOfDaughters() > 0) BuildTo(prt, hist, daughters,  mom_pdgid, map_mom, num_daugh, daugh_pdgid, map_daugh);
	}
}

//pi0 111 ;100111 9010111 ...
//pi+ 211 ;100211 9010211 ...
//k0l 130 k0s 310 k0 311...
void AOD_pi0::FindFinalPrt(vector<const reco::Candidate *> *d = 0,
														const reco::Candidate *mom = 0 ,
														int num_of_final_daughters = 2,
														int daug_pdgid = 111)// if charged - of opposite charge
{
	if (d == 0 || mom == 0)
	{
		dlog("AOD_pi0::FindFinalPrt: wrong parameter");
		return;
	}

	int n = mom->numberOfDaughters();
	for(int i = 0; i < n; i++)
	{
		if (num_of_final_daughters == 0) break;
		const reco::Candidate * prt = mom->daughter(i);
		if (abs(prt->pdgId()) == daug_pdgid && NoRadiating(prt, daug_pdgid))// prt->numberOfDaughters() == 0
		{
			dout("\t\tfound one", prt->pdgId());
			d->push_back(prt);
			num_of_final_daughters--;
		}
		else if (prt->numberOfDaughters() > 0) FindFinalPrt(d, prt, num_of_final_daughters, daug_pdgid);
	}
	//dlog("\t\t\tFor this cand not found:", num_of_final_daughters);
}

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

int AOD_pi0::CountGenPrtID(const reco::Candidate * mom, int searched_id = 310, int number = 0, int num_of_tabs = 2, std::vector<TH1D *>* vector_of_histograms=0)
{
	int n = mom->numberOfDaughters();
		for(int i = 0; i < n; i++)
		{
			const reco::Candidate * daughter = mom->daughter(i);
			//if (abs(daughter->pdgId()) == 311/*(abs(daughter->pdgId()) % 1000 == 311 || abs(daughter->pdgId()) % 1000 == 321)*/ && daughter->numberOfDaughters() == 0 ) {dlog(string(num_of_tabs, '\t'), i, ") is Kaon", daughter->pdgId(), daughter->status(), daughter->numberOfDaughters());exit(1);}
			//cout << string(num_of_tabs, '\t') << i << ") "<< number<< " | " << daughter->pdgId() << " status:" << daughter->status() << " nmom:" << daughter->numberOfMothers()  << " ndaugters:" << daughter->numberOfDaughters() << endl;//dlog(string(num_of_tabs, '\t'), i, ") ", daughter->pdgId(), daughter->status(), daughter->numberOfDaughters());
			if (searched_id == abs(daughter->pdgId()))
			{
				number++;
				if (vector_of_histograms->size() > 0)
				{
					(*vector_of_histograms)[1]->Fill(daughter->vx());
					(*vector_of_histograms)[2]->Fill(daughter->vy());
					(*vector_of_histograms)[3]->Fill(daughter->vz());
				}
			}
			else number += CountGenPrtID(daughter, searched_id, 0, num_of_tabs+1, vector_of_histograms);
		}
		return number;
}

void AOD_pi0::GenEvolution(const reco::Candidate * mom, int num_of_tabs = 2)
{
		int n = mom->numberOfDaughters();
		for(int i = 0; i < n; i++)
		{
			const reco::Candidate * daughter = mom->daughter(i);
			//if (abs(daughter->pdgId()) == 311/*(abs(daughter->pdgId()) % 1000 == 311 || abs(daughter->pdgId()) % 1000 == 321)*/ && daughter->numberOfDaughters() == 0 ) {dlog(string(num_of_tabs, '\t'), i, ") is Kaon", daughter->pdgId(), daughter->status(), daughter->numberOfDaughters());exit(1);}
			cout << string(num_of_tabs, '\t') << i << ") " << daughter->pdgId() << " " << daughter->status() << " " << daughter->numberOfDaughters() << endl;//dlog(string(num_of_tabs, '\t'), i, ") ", daughter->pdgId(), daughter->status(), daughter->numberOfDaughters());
			GenEvolution(daughter, num_of_tabs + 1);
		}
}

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

template  <typename VectoreType >
void AOD_pi0::CombinatoricOfTwoToNeutralInvM(vector <VectoreType *> collection,
																		TString typeOfCollection,
																		TString typeOfObjects,
																		TString typeOfConstituences,
																		TH1 * hist_inv_m,
																		TH1 * hist_pt)
{
	dlog();
	dlog("\t\t", typeOfCollection, "'s number:", collection.size());
	if (collection.size() > 1)
	{
		for (unsigned int i = 0; i < collection.size() - 1; i++) //over isolation pions in Tau
		{
			dlog("\t\t\t", typeOfCollection, i, ":", collection[i]->vx(), collection[i]->vy(), collection[i]->vz(),
																					":", collection[i]->px(), collection[i]->py(), collection[i]->pz());

			//Check the vertex position
			 if (false) dout("\t\t\t", typeOfConstituences, i, "(", collection[i]->charge(), "), vertex:", collection[i]->vx(), collection[i]->vy(), collection[i]->vz()); //all from the same vertex
			 TLorentzVector first(collection[i]->px(), collection[i]->py(), collection[i]->pz(), collection[i]->energy());
			//Build the inv mass
				dout("\t\t===>matching to", typeOfConstituences , i, "from", collection.size(), "in this", typeOfObjects);
				for (unsigned int j = i + 1; j < collection.size(); j++)
				{
					if (typeOfConstituences.Contains("pi+-") && collection[i]->charge() * collection[j]->charge() < 0) continue;// combine only to neutral particles K0
					TLorentzVector second(collection[j]->px(), collection[j]->py(), collection[j]->pz(), collection[j]->energy());
					double inv_M = (first + second).M();
					hist_inv_m->Fill(inv_M); //This will give us a combinatiric bg
					dout("\t\t\tm(", typeOfConstituences, i, "+", typeOfConstituences, j, ") =", inv_M);
					if (hist_pt != 0) hist_pt->Fill((first + second).Pt()); //This will give us a combinatiric bg
				}
				//cout,  endl;
				//if (hist_pt != 0) hist_pt->Fill(collection[i]->pt());
		}
		//if (hist_pt != 0) hist_pt->Fill(collection[collection.size() - 1]->pt());
		dlog("\t\t\t", typeOfCollection, collection.size() - 1, ":", collection[collection.size() - 1]->vx(), collection[collection.size() - 1]->vy(), collection[collection.size() - 1]->vz(),
																														":", collection[collection.size() - 1]->px(), collection[collection.size() - 1]->py(), collection[collection.size() - 1]->pz());
	}
	else if (collection.size() == 1) dlog("\t\t\t", typeOfCollection, 0, " vtx:", collection[0]->vx(), collection[0]->vy(), collection[0]->vz(), " || p:", collection[0]->px(), collection[0]->py(), collection[0]->pz());
}

template  <typename VectoreType >
void AOD_pi0::CombinatoricOfTwoToNeutralInvM(vector <VectoreType> collection, // doubt that this function is needed - check it
																		TString typeOfCollection,
																		TString typeOfObjects,
																		TString typeOfConstituences,
																		TH1 * hist_inv_m,
																		TH1 * hist_pt)
{
	dlog();
	dlog("\t\t", typeOfCollection, "'s number:", collection.size());
	if (collection.size() > 1)
	{
		for (unsigned int i = 0; i < collection.size() - 1; i++) //over isolation pions in Tau
		{
			dlog("\t\t\t", typeOfCollection, i, ":", collection[i]->vx(), collection[i]->vy(), collection[i]->vz(),
																					":", collection[i]->px(), collection[i]->py(), collection[i]->pz());

			//Check the vertex position
			 if (false) dout("\t\t\t", typeOfConstituences, i, "(", collection[i]->charge(), "), vertex:", collection[i]->vx(), collection[i]->vy(), collection[i]->vz()); //all from the same vertex
			 TLorentzVector first(collection[i]->px(), collection[i]->py(), collection[i]->pz(), collection[i]->energy());
			//Build the inv mass
				dout("\t\t===>matching to", typeOfConstituences , i, "from", collection.size(), "in this", typeOfObjects);
				for (unsigned int j = i + 1; j < collection.size(); j++)
				{
					if (typeOfConstituences.Contains("pi+-") && collection[i]->charge() * collection[j]->charge() < 0) continue;
					TLorentzVector second(collection[j]->px(), collection[j]->py(), collection[j]->pz(), collection[j]->energy());
					double inv_M = (first + second).M();
					hist_inv_m->Fill(inv_M); //This will give us a combinatiric bg
					dout("\t\t\tm(", typeOfConstituences, i, "+", typeOfConstituences, j, ") =", inv_M);
					if (hist_pt != 0) hist_pt->Fill((first + second).Pt()); //This will give us a combinatiric bg
				}
				// if (hist_pt != 0) hist_pt->Fill(collection[i]->pt());
				//cout,  endl;
		}
		//if (hist_pt != 0) hist_pt->Fill(collection[collection.size() - 1]->pt());
		dlog("\t\t\t", typeOfCollection, collection.size() - 1, " vtx:", collection[collection.size() - 1]->vx(), collection[collection.size() - 1]->vy(), collection[collection.size() - 1]->vz(), " || p:", collection[collection.size() - 1]->px(), collection[collection.size() - 1]->py(), collection[collection.size() - 1]->pz());
	}
	else if (collection.size() == 1) dlog("\t\t\t", typeOfCollection, 0, " vtx:", collection[0]->vx(), collection[0]->vy(), collection[0]->vz(), " || p:", collection[0]->px(), collection[0]->py(), collection[0]->pz());
}

unsigned int AOD_pi0::AddGammas(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	iEvent.getByToken(PFCandidateCollectionToken_, Tracks);

	photon_count = 0;
	std::vector<TLorentzVector> combined_pions;

	if (Tracks.isValid())
	{
		for(unsigned i = 0; i < Tracks->size(); i++)
		{
			if ((*Tracks)[i].pdgId() != 22) continue;
			if ((*Tracks)[i].pt() < cPhotonPtMin) continue;
			if (fabs((*Tracks)[i].eta()) > cPhotonEtaMax) continue;
			//if ((*Tracks)[i].Charge() != 0) continue;

			// photon_px[photon_count]  = (*Tracks)[i].px();
			// photon_py[photon_count]  = (*Tracks)[i].py();
			// photon_pz[photon_count]  = (*Tracks)[i].pz();
			// photon_pt[photon_count]  = (*Tracks)[i].pt();
			// photon_eta[photon_count] = (*Tracks)[i].eta();
			// photon_phi[photon_count] = (*Tracks)[i].phi();
			// photon_e[photon_count]   = (*Tracks)[i].p();
			TLorentzVector photon_1 = TLorentzVector((*Tracks)[i].px(), (*Tracks)[i].py(), (*Tracks)[i].pz(), (*Tracks)[i].energy());
			//dout("M:", photon_1.M());
			//if (abs(photon_1.M() - 0.135)<0.05) exit(1);
			//const reco::SuperClusterRef superCluster = (*Tracks)[i].superClusterRef();//edm::Ref<SuperClusterCollection> reco::SuperClusterRef
			// photon_superClusterX[photon_count] = superCluster->x();
			// photon_superClusterY[photon_count] = superCluster->y();
			// photon_superClusterZ[photon_count] = superCluster->z();
			//nope double mass_1 = superCluster->m();
			photon_count++;
			for(unsigned j = i+1; j < Tracks->size(); j++)
			{
				if (j >= Tracks->size()) break;
				if ((*Tracks)[j].pdgId() != 22) continue;
				if ((*Tracks)[j].pt() < cPhotonPtMin) continue;
				if (fabs((*Tracks)[j].eta()) > cPhotonEtaMax) continue;
				//if ((*Tracks)[j].Charge() != 0) continue;
				TLorentzVector photon_2 = TLorentzVector((*Tracks)[j].px(), (*Tracks)[j].py(), (*Tracks)[j].pz(), (*Tracks)[j].energy());
				combined_pions.push_back((photon_1+photon_2));

				h_ECAL_comb_photons->Fill((photon_2+photon_1).M()); //around zero
			}
			// if (photon_count == M_photonmaxcount)
			// {
			//   cerr << "number of tracks > M_trackmaxcount. They are missing." << endl; errors |= 1<<1;
			//   break;
			// }
		}

		for(std::vector<TLorentzVector>::iterator itv = combined_pions.begin(); itv != combined_pions.end() - 1; ++itv)
		{

			for(std::vector<TLorentzVector>::iterator itv2 = itv + 1; itv2 != combined_pions.end(); ++itv2)
			{
				if (itv2 >= combined_pions.end()) break;
				TLorentzVector kaon_candidate = (*itv) + (*itv2);
				//if (abs(kaon_candidate.M() - 0.492) < 0.07)
				{
					h_ECAL_comb_kaons->Fill(kaon_candidate.M());
					//dout("made ks"); exit(1);
				}
			}
		}
	}

	return photon_count;
}

//define this as a plug-in
DEFINE_FWK_MODULE(AOD_pi0);
