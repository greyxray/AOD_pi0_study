#include "AOD_pi0.h"

void AOD_pi0::GeneralStudy(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
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

}