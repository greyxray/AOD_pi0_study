#include "AOD_pi0.h"

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