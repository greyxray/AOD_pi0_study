#include "AOD_pi0.h"

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
