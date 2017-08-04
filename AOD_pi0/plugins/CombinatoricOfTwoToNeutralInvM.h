
#include "AOD_pi0.h"

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
