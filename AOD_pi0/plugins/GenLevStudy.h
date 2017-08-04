#include "AOD_pi0.h"

void AOD_pi0::GenLevStudy(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	if (!IsData) iEvent.getByToken(GenParticleCollectionToken_, GenPart);
	else return;
	iEvent.getByToken(TauPiZeroCollectionToken_, Strips); //actually HPS pi0s
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
}