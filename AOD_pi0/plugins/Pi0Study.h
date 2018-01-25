#include "AOD_pi0.h"

void AOD_pi0::Pi0Study(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	iEvent.getByToken( TauPiZeroCollectionToken_, Strips); //actually HPS pi0s
	/// HPS Pi0's and taus - with loop among reco::tau
	if (Strips.isValid())
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
}