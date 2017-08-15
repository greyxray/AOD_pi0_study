// -*- C++ -*-
//
// Package:    AOD_pi0_study/AOD_pi0
// Class:      AOD_pi0
//
/**\class AOD_pi0 AOD_pi0.cc AOD_pi0_study/AOD_pi0/plugins/AOD_pi0.cc

 Description: build a ks resonan

*/
#include "AOD_pi0.h"
#include "AddGammas.h"
#include "CombinatoricOfTwoToNeutralInvM.h"
#include "AdditionalTools.h"
#include "NoRadiating.h"
#include "GenLevAdditionalTools.h"
#include "Pi0Study.h"
#include "K892.h"
#include "GeneralStudy.h"
#include "KFromV0Producer.h"
#include "FillPion.h"
#include "FWCore/Framework/interface/MakerMacros.h" //define this as a plug-in

void AOD_pi0::ResetBranchesPerEvent()
{
	// tree_once
	v0_Ks_count = -1;
	primvertex_count = -1;
	goodprimvertex_count = -1;
	numOfUnmatchedKaons = -1;
	v_v0_Ks_inv_m_pi.clear();
	v_v0_pions_DR.clear();
	v_v0_pt_1.clear();
	v_v0_pt_2.clear();
	v_v0_eta_1.clear();
	v_v0_eta_2.clear();
	v_v0_matched_pt_1.clear();
	v_v0_matched_pt_2.clear();
	v_v0_matched_eta_1.clear();
	v_v0_matched_eta_2.clear();
	v_KsCombinatoricMass.clear();
	v_KsCombinatoricDR.clear();
	v_hps_Ks_inv_m_pi.clear();
	v_hps_Ks_DR.clear();
	v_hps_Ks_pions_DR.clear();
	v_hps_pt_1.clear();
	v_hps_pt_2.clear();
	v_hps_eta_1.clear();
	v_hps_eta_2.clear();
	v_nPionsInJetsWithKs.clear();

	pizero_count = -1;

	ResetBranchesKaonTree();
	ResetBranchesPionTree();
}

// ------------ method called for each event  ------------
void AOD_pi0::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	ResetBranchesPerEvent();

	iEvent.getByToken( KshortCollectionToken_, V0Ks); //V0 Ks's token
	iEvent.getByToken( KshortCollectionTag_stand_, V0Ks_standart); //V0 Ks's standart
	iEvent.getByToken( BeamSpotToken_, TheBeamSpotHandle);
	//Tau's token based hps- according to https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookPATDataFormats#PatTau contains the same as PFjets
	//typedef vector< PFTau >  PFTauCollection - the collection of charged pions will be taken from here
	iEvent.getByToken( TauHPSCollectionToken_, PF_hps_taus);
	// Magnetic field
		// edm::ESHandle<MagneticField> theMagneticFieldHandle;
		// iSetup.get<IdealMagneticFieldRecord>().get(theMagneticFieldHandle);
		// const MagneticField* theMagneticField = theMagneticFieldHandle.product();

		//  edm::Handle<reco::TrackCollection> HPSTraHandle;
	iEvent.getByToken( HPSTrackTagToken_, HPSTraHandle);
		// if (HPSTraHandle.isValid())  dout("HPStrackCollection is fine");
		// else dout("HPStrackCollection is NOT VALID", HPSTraHandle->size());
		// // const reco::TrackCollection* theTrackCollection = HPSTraHandle.product();

	/// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	/// PV -> pv_position
	if (crecprimvertex)
	{
		iEvent.getByToken(PVToken_, Vertex);

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

	if (true) KFromV0Producer(iEvent, iSetup);
	if (false) AddGammas(iEvent, iSetup);
	if (false) GenLevStudy(iEvent, iSetup);
	if (false) Pi0Study(iEvent, iSetup);
	if (false) K892(iEvent, iSetup);
	if (false) GeneralStudy(iEvent, iSetup);

	once_tree->Fill();
}

//using PDGMassConstants;
//define this as a plug-in
DEFINE_FWK_MODULE(AOD_pi0);
