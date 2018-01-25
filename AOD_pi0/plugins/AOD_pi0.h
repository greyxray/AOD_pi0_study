#ifndef AOD_pi0_H
#define AOD_pi0_H

#include <memory>
#include <map>
#include <fstream>      // std::ofstream
#include <iostream>
#include <set>
#include <vector>
#include <string>
#include <utility>
#include <TTree.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// Additional include
//error #include "RecoTauTag/RecoTau/interface/PFRecoTauClusterVariables.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

//ROOT
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TString.h"


using namespace std;

// Aleksei
	#include "DataFormats/VertexReco/interface/VertexFwd.h"
	#include "DataFormats/VertexReco/interface/Vertex.h"
	#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
	#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
	//error using namespace reco;

	#include "DataFormats/TauReco/interface/PFTau.h"
	#include "DataFormats/TauReco/interface/PFTauFwd.h"
	#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
	#include "DataFormats/TauReco/interface/RecoTauPiZero.h"
	#include "DataFormats/TauReco/interface/RecoTauPiZeroFwd.h"
	#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

// Alex
	#include "DataFormats/PatCandidates/interface/Tau.h"

//new
	#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
	#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
	#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
	#include "PhysicsTools/HepMCCandAlgos/plugins/MCTruthDeltaRMatcherNew.cc"

// TFileService
	#include "FWCore/ServiceRegistry/interface/Service.h"
	#include "CommonTools/UtilAlgos/interface/TFileService.h"

// From the V0
	#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
	#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
	#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
	#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
	#include "MagneticField/VolumeBasedEngine/interface/VolumeBasedMagneticField.h"

#include "TrackingTools/IPTools/interface/IPTools.h"

//#include "PhysicsTools/HepMCCandAlgos/interface/MCTruthCompositeMatcher.h" // includes missed files
//#include "RecoTauPiZero.h"
/*
	#include <string>
	#include <map>
	#include <vector>
	#include <cstdlib>
	#include <algorithm>

	#include <Math/Vector3D.h>
	#include "Math/LorentzVector.h"
	#include "Math/Point3D.h"

	#include <boost/regex.hpp>
	#include <boost/algorithm/string.hpp>

	#include "FWCore/Utilities/interface/InputTag.h"

	#include "FWCore/Framework/interface/LuminosityBlock.h"
	#include "DataFormats/Luminosity/interface/LumiSummary.h"

	#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
	#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
	#include "DataFormats/Common/interface/TriggerResults.h"
	#include "DataFormats/HLTReco/interface/TriggerEvent.h"
	#include "FWCore/Common/interface/TriggerNames.h"

	#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
	#include "CondFormats/L1TObjects/interface/L1GtPrescaleFactors.h"
	#include "CondFormats/DataRecord/interface/L1GtPrescaleFactorsAlgoTrigRcd.h"
	#include "CondFormats/DataRecord/interface/L1GtPrescaleFactorsTechTrigRcd.h"

	#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
	#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
	//error #include "DataFormats/TrackReco/interface/Track.h"
	#include "DataFormats/TrackReco/interface/TrackFwd.h"
	//error #include "TrackingTools/TransientTrack/interface/TransientTrack.h"
	#include "DataFormats/EgammaReco/interface/SuperCluster.h"
	#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
	#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
	#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
	//error #include "DataFormats/EgammaCandidates/interface/Photon.h"
	#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
	//error #include "DataFormats/EgammaCandidates/interface/Conversion.h"
	#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
	//error #include "DataFormats/MuonReco/interface/Muon.h"
	#include "DataFormats/MuonReco/interface/MuonFwd.h"
	//error #include "DataFormats/MuonReco/interface/MuonSelectors.h"

	//#include "DataFormats/METReco/interface/MET.h"
	//error #include "DataFormats/PatCandidates/interface/MET.h"
	#include "DataFormats/METReco/interface/METFwd.h"
	#include "DataFormats/METReco/interface/CaloMET.h"
	#include "DataFormats/METReco/interface/CaloMETFwd.h"
	#include "DataFormats/METReco/interface/PFMET.h"
	#include "DataFormats/METReco/interface/PFMETFwd.h"
	#include "DataFormats/JetReco/interface/CaloJet.h"
	#include "DataFormats/JetReco/interface/CaloJetCollection.h"
	#include "DataFormats/JetReco/interface/PFJet.h"
	#include "DataFormats/JetReco/interface/PFJetCollection.h"
	#include "DataFormats/JetReco/interface/Jet.h"
	#include "DataFormats/JetReco/interface/JetCollection.h"
	#include "JetMETCorrections/Objects/interface/JetCorrector.h"
	//error #include "DataFormats/VertexReco/interface/Vertex.h"
	#include "DataFormats/VertexReco/interface/VertexFwd.h"


	#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
	#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
	#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
	#include "DataFormats/METReco/interface/GenMET.h"
	#include "DataFormats/METReco/interface/GenMETFwd.h"
	#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

	#include "TrackingTools/Records/interface/TransientTrackRecord.h"
	//error #include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
	//error #include <DataFormats/MuonReco/interface/Muon.h>
	#include <DataFormats/MuonReco/interface/MuonFwd.h>
	#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
	#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"


	#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
	#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
	#include "Geometry/Records/interface/CaloGeometryRecord.h"

	#include "RecoEgamma/EgammaTools/interface/ConversionInfo.h"
	//error #include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"

	#include "TMath.h"
	#include "TLorentzVector.h"
	#include <Math/Functions.h>
	#include <Math/SVector.h>
	#include <Math/SMatrix.h>

	#include "DataFormats/GeometrySurface/interface/SimpleCylinderBounds.h"
	#include "DataFormats/GeometrySurface/interface/SimpleDiskBounds.h"
	#include "DataFormats/GeometrySurface/interface/Cylinder.h"
	#include "DataFormats/GeometrySurface/interface/Plane.h"
	#include "DataFormats/GeometrySurface/interface/BoundCylinder.h"
	#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
	#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
	//error #include "TrackingTools/TransientTrack/interface/TransientTrack.h"
	#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
	#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
	//error #include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
	//error #include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
	//error #include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
	#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

	//error #include "EgammaAnalysis/ElectronTools/interface/EGammaMvaEleEstimatorCSA14.h"
*/

const int MAXNHISTOS = 15;
const int MAXKS = 500;

// pdg mass constants
namespace PDGMassConstants
{
	const double piMass = 0.13957018;
	const double piMassSquared = piMass*piMass;
	const double protonMass = 0.938272046;
	const double protonMassSquared = protonMass*protonMass;
	const double kShortMass = 0.497614;
	const double lambdaMass = 1.115683;
}

typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
typedef ROOT::Math::SVector<double, 3> SVector3;
typedef ROOT::Math::PtEtaPhiM4D<float> RMFLV_Store;
typedef ROOT::Math::LorentzVector<RMFLV_Store> RMFLV;
typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float> > RMPoint;
typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag> Point;

// comment on this: https://twiki.cern.ch/twiki/bin/view/CMSPublic/FWMultithreadedAnalysisEDAnalyzer#Using_TFileService
class AOD_pi0 : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
	public:

		explicit  AOD_pi0(const edm::ParameterSet&);
		~AOD_pi0();

		unsigned int AddGammas(const edm::Event& iEvent, const edm::EventSetup& iSetup);
		void GenLevStudy(const edm::Event& iEvent, const edm::EventSetup& iSetup);
		void Pi0Study(const edm::Event& iEvent, const edm::EventSetup& iSetup);
		void K892(const edm::Event& iEvent, const edm::EventSetup& iSetup);
		void GeneralStudy(const edm::Event& iEvent, const edm::EventSetup& iSetup);
		void KFromV0Producer(const edm::Event& iEvent, const edm::EventSetup& iSetup);

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
		void ResetBranchesPerEvent();
		void ResetBranchesPionTree();
		void ResetBranchesKaonTree();
		static void RecO_Cand_type(const reco::Candidate* cand);
		bool BetterPionMatch(TLorentzVector tau_pion, int chargeTauPion, TLorentzVector v0_ks_pion, int chargeKPion, int firstfound, double & dR);

		template<class T, class P, class R>
		float getDxy(const T pv, const P p4, const R ref) const;
		template<class T, class P, class R>
		float getDz(const T pv, const P p4, const R ref) const;

		template  <typename HPSPion, typename Daughter>
		bool PionMatchByRefference(HPSPion pion, Daughter daughter);

		template  <typename HPSPion, typename Daughter>
		void FillPion(HPSPion pion, Daughter daughter);

		template  <typename VectoreType >
		void CombinatoricPairs(vector<VectoreType> collection);

		template  <typename VectoreType >
		void CombinatoricOfTwoToNeutralInvM(vector<VectoreType*> collection/*reco::RecoTauPiZero */,
																		TString typeOfCollection,
																		TString typeOfObjects,
																		TString typeOfConstituences,
																		TH1 * hist_inv_m,
																		TH1 * hist_pt=0);
		template  <typename VectoreType >
		void CombinatoricOfTwoToNeutralInvM(vector<VectoreType> collection/*reco::RecoTauPiZero */,
																		TString typeOfCollection,
																		TString typeOfObjects,
																		TString typeOfConstituences,
																		TH1 * hist_inv_m,
																		TH1 * hist_pt=0);
		template <typename T>
		vector <T*> TransformToPointers(vector <T> a, vector <T*> b);
		void GenEvolution(const reco::Candidate * , int);
		int CountGenPrtID(const reco::Candidate *, int, int, int, std::vector<TH1D *>*);

		void BuildTo(edm::Handle<std::vector<reco::GenParticle> >& genPart, TH1 *hist, vector< vector <const reco::Candidate *>> &daughters, int mom_pdgid, map <long, string>& , int num_daugh, int daugh_pdgid, map <long, string>& );
		void BuildTo(const reco::Candidate *mom, TH1 *hist, vector< vector <const reco::Candidate *>> &daughters, int mom_pdgid, map <long, string>& , int num_daugh, int daugh_pdgid, map <long, string>& );
		void FindFinalPrt(vector<const reco::Candidate *> *d, const reco::Candidate *mom , int num_of_final_daughters, int daug_pdgid);
		bool NoRadiating(const reco::Candidate * mom, int daug_pdgid);
		bool NoRadiating(const reco::Candidate * mom, map <long, string>& map_daug_pdgid);

		template <typename T>
		void Match(vector< vector <const reco::Candidate *>>& From,
										vector < T *> & To,
										vector < vector< vector <T *>>>& SimToReco,
										vector < vector< const reco::Candidate *>>& RecoToSim
										//multiset < pair< int*, vector <string*> > >
										);
		template <typename Ta, typename Tb>
		bool Incone(Ta& A, Tb& B, double cone_size);
		template <typename Tc, typename Tr>// T - collection   reco::PFTauRef pftauref(PF_taus, i);
		void MakeVectorofRef(edm::Handle< Tc > Collection, vector< Tr* > v_of_ref);

		template<typename T>
		struct is_pointer { static const bool value = false; };
		template<typename T>
		struct is_pointer<T*> { static const bool value = true; };

		static void dout();
		template <typename Head, typename... Tail>
		static void dout(Head, Tail... );
		static void dlog();
		template <typename Head, typename... Tail>
		static void dlog(Head, Tail... );

	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;

	 // ----------member data ---------------------------
		edm::Service<TFileService> fs;
		TTree * kaon_tree; // example: http://www-hep.colorado.edu/~fjensen/temp/trackAnalyzer.cc
		TTree * once_tree;
		TTree * pions_tree;
		TTree * kaon892_tree;

		// Branches variables

			// kaon_tree
				Int_t nPionsInJetsWithKs;

				// V0 comparison to HPS
				Double_t v0hps_pions_dPhi_1;
				Double_t v0hps_pions_dEta_1;
				Double_t v0hps_pions_dPhi_2;
				Double_t v0hps_pions_dEta_2;
				Double_t v0hps_KsTau_DR;
				Double_t v0hps_KsDiPion_DR;
				// HPS Kaon
				Double_t hps_Ks_inv_m_pi;
				Double_t hps_pions_DR;
				Double_t hps_deta;
				Double_t hps_dphi;
				Double_t hps_px;
				Double_t hps_py;
				Double_t hps_pz;
				Double_t hps_energy;
				// HPS neg pion
				Double_t hps_pt_1;
				Double_t hps_eta_1;
				Double_t hps_phi_1;
				Double_t hps_px_1;
				Double_t hps_py_1;
				Double_t hps_pz_1;
				Double_t hps_energy_1;
				// HPS pos pion
				Double_t hps_pt_2;
				Double_t hps_eta_2;
				Double_t hps_phi_2;
				Double_t hps_px_2;
				Double_t hps_py_2;
				Double_t hps_pz_2;
				Double_t hps_energy_2;
				// V0 kaon
				Double_t v0_Ks_inv_m_pi;
				Double_t v0_Ks_pions_DR;
				Double_t v0_deta;
				Double_t v0_dphi;
				Double_t v0_pt;
				Double_t v0_eta;
				Double_t v0_phi;
				Double_t v0_px;
				Double_t v0_py;
				Double_t v0_pz;
				Double_t v0_energy;
				// V0 neg pion
				Double_t v0_pt_1;
				Double_t v0_eta_1;
				Double_t v0_phi_1;
				Double_t v0_px_1;
				Double_t v0_py_1;
				Double_t v0_pz_1;
				Double_t v0_energy_1;
				// V0 pos pion
				Double_t v0_pt_2;
				Double_t v0_eta_2;
				Double_t v0_phi_2;
				Double_t v0_px_2;
				Double_t v0_py_2;
				Double_t v0_pz_2;
				Double_t v0_energy_2;

			// pions_tree
				TLorentzVector v0_pion;
				TLorentzVector v0_pion_trackRef;
				TLorentzVector tau_pion;
				TLorentzVector tau_pion_besttrack;
				// Diff between HPS and HPS best track
				Double_t hpsbestTrackDiff_fx;
				Double_t hpsbestTrackDiff_fy;
				Double_t hpsbestTrackDiff_fz;
				Double_t hpsbestTrackDiff_fE;
				Double_t hpsbestTrack_dR;
				// Diff between V0 and HPS
				Double_t hpsv0Diff_fx;
				Double_t hpsv0Diff_fy;
				Double_t hpsv0Diff_fz;
				Double_t hpsv0Diff_fE;
				Double_t hpsv0Diff_dPhi;
				Double_t hpsv0Diff_dEta;
				Double_t hpsv0_dR;
				// Diff between V0 and HPS best track
				Double_t v0bestTrackDiff_fx;
				Double_t v0bestTrackDiff_fy;
				Double_t v0bestTrackDiff_fz;
				Double_t v0bestTrackDiff_fE;
				Double_t v0bestTrack_dR;
				// associated mass
				Double_t bestTrac_fM;
				Double_t hps_fM;
				Double_t v0_fM;

			// once_tree
				UInt_t primvertex_count;
				UInt_t goodprimvertex_count;
				Int_t v0_Ks_count; // V0s - Ks or Lambdas
				unsigned int numOfUnmatchedKaons;
				std::vector<double> v_v0_Ks_inv_m_pi; //Int_t Ks_v0_inv_m_pi[MAXKS];
				std::vector<double> v_v0_pions_DR;
				std::vector<double> v_v0_Ks_pions_DR;
				std::vector<Double_t> v_v0_pt_1;
				std::vector<Double_t> v_v0_pt_2;
				std::vector<Double_t> v_v0_eta_1;
				std::vector<Double_t> v_v0_eta_2;
				std::vector<Double_t> v_v0_matched_pt_1;
				std::vector<Double_t> v_v0_matched_pt_2;
				std::vector<Double_t> v_v0_matched_eta_1;
				std::vector<Double_t> v_v0_matched_eta_2;
				std::vector<Int_t> v_v0_count;
				std::vector<Int_t> v_nPionsInJetsWithKs;
				std::vector<double> v_KsCombinatoricMass;
				std::vector<double> v_KsCombinatoricDR;
				std::vector<double> v_hps_Ks_inv_m_pi;
				std::vector<double> v_hps_Ks_DR;
				std::vector<double> v_hps_Ks_pions_DR;
				std::vector<Double_t> v_hps_pt_1;
				std::vector<Double_t> v_hps_pt_2;
				std::vector<Double_t> v_hps_eta_1;
				std::vector<Double_t> v_hps_eta_2;

			//kaon892_tree
				TLorentzVector lv_K892;
				double K892_fx;
				double K892_fy;
				double K892_fz;
				double K892_fE;
				double K892_fM;

				TLorentzVector lv_K892_K0;
				double K892_K0_fx;
				double K892_K0_fy;
				double K892_K0_fz;
				double K892_K0_fE;
				double K892_K0_dxy;
				double K892_K0_dz;
				double K892_K0_fM;

				TLorentzVector lv_K892_pi0;
				double K892_pi0_fx;
				double K892_pi0_fy;
				double K892_pi0_fz;
				double K892_pi0_fE;
				double K892_pi0_fM;

				long n_K892;

		// Histograms
			TH1D* taus_isol_pi0_inv_m_to_ks;
			TH1D* taus_isol_pi0_inv_pt;
			TH1D* taus_pi0_inv_m_to_ks;
			TH1D* taus_pi0_inv_pt;
			TH1D* taus_pi_charged_inv_m_to_ks;
			TH1D* taus_pi_charged_inv_pt;
			TH1D* taus_pi0_had_inv_m_to_ks;
			TH1D* taus_pi0_had_inv_pt;
			TH1D* h_ECAL_comb_photons;
			TH1D* h_ECAL_comb_kaons;

		// Tokens for the Collections
			edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> KshortCollectionToken_;
			edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> KshortCollectionTag_stand_;
			edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> LambdaCollectionToken_;
			edm::EDGetTokenT<reco::PFCandidateCollection> PFCandidateCollectionToken_;
			edm::EDGetTokenT<reco::RecoTauPiZeroCollection> TauPiZeroCollectionToken_;
			edm::EDGetTokenT<reco::PFTauCollection> TauHPSCollectionToken_;
			// vector<int>                           "genParticles"              ""                "HLT"
			// vector<reco::GenParticle>             "genParticles"              ""
			edm::EDGetTokenT<reco::GenParticleCollection> GenParticleCollectionToken_;
			edm::EDGetTokenT<reco::VertexCollection> PVToken_;
			edm::EDGetTokenT<reco::BeamSpot> BeamSpotToken_;
			edm::EDGetTokenT<reco::TrackCollection> HPSTrackTagToken_;

		// Handles
			edm::Handle<reco::VertexCompositeCandidateCollection> V0Ks;
			edm::Handle<reco::VertexCompositeCandidateCollection> V0Ks_standart;
			edm::Handle<reco::RecoTauPiZeroCollection> Strips;
			edm::Handle<reco::GenParticleCollection> GenPart;
			edm::Handle<reco::BeamSpot> TheBeamSpotHandle;
			edm::Handle<reco::PFTauCollection> PF_hps_taus;
			edm::Handle<reco::VertexCollection> Vertex;
			edm::Handle<reco::PFCandidateCollection> Tracks;
			edm::Handle<reco::TrackCollection> HPSTraHandle;

		// Variables
			std::ofstream ofs;
				//std::streambuf * buf;
				//std::ostream ofs;//ostream ofstream
				// std::ofstream of;
				// std::ostream & ofs = (!outputGenEvolution) ? cout : of.open("geninfo.txt", std::ofstream::out /*| std::ofstream::app*/);
				// std::ostream* ofs;
				// std::ofstream byname;
				// std::ostream& ofs;

			// pi0's
				UInt_t pizero_count;
				Float_t pizero_px[1000];
				Float_t pizero_py[1000];
				Float_t pizero_pz[1000];
					Float_t pizero_pt[1000];
					Float_t pizero_eta[1000];
					Float_t pizero_phi[1000];
					Float_t pizero_e[1000];
				Float_t pizero_x[1000];
				Float_t pizero_y[1000];
				Float_t pizero_z[1000];

			// primary vertex
				math::XYZPoint pv_position;
				reco::Vertex primvertex;
				Float_t primvertex_x;
				Float_t primvertex_y;
				Float_t primvertex_z;
				Float_t primvertex_chi2;
				Float_t primvertex_ndof;
				Float_t primvertex_ptq;
				Int_t   primvertex_ntracks;
				Float_t primvertex_cov[6];

			// K0
				vector< vector <const reco::Candidate *>> v_daughters_k0s_to_pi0;
				vector< vector <const reco::Candidate *>> v_daughters_k0l_to_pi0;
				vector< vector <const reco::Candidate *>> v_daughters_k0_to_pi0;
				vector< vector <const reco::Candidate *>> v_daughters_k0s_to_pic;
				vector< vector <const reco::Candidate *>> v_daughters_k0l_to_pic;
				vector <reco::RecoTauPiZero *> v_strips_ref;
				vector < vector< vector <reco::RecoTauPiZero *>>> v_daughters_k0s_to_pi0_SimToStrips;
				vector < vector< const reco::Candidate *>> v_daughters_k0s_to_pi0_StipsToSim;

			// maps
				/*static*/ std::map <long, string> map_kaons; //could also be a std::set
				/*static*/ std::map <long, string> map_pions;
				std::map<long, std::vector<TH1D *>> map_K892_gen_histos;

		// Parameters initialised by default
			bool IsData;
			string OutFileName;
			// switches
				bool cgen;
				bool cbeamspot;
				bool crecprimvertex;
				bool crectrack;
				bool crecphoton;
				bool crecpizero;
				bool crecsv;
				bool crecpfjet;
				bool crecv0;
			// tracks
				double cTrackPtMin;
				double cTrackEtaMax;
				double cTrackDxyMax;
				double cTrackDzMax;
				int cTrackNum;
			// photons
				UInt_t photon_count;
				double cPhotonPtMin;
				double cPhotonEtaMax;
				int cPhotonNum;
			// pizeros
				double cPizeroPtMin;
				double cPizeroEtaMax;
				int cPizeroNum;
			// jets
				double cJetPtMin;
				double cJetEtaMax;
				int cJetNum;
			// matched Kaons
				bool match_KsV0_to_HPS;
				double cDZCut;
				double cKtoTauDR;
				bool matchByReference;
			// other
				bool debug;
				bool mute;
				double tkIPSigXYCut;
				double vtxDecaySigXYCut;

			int num_pion_res;
			double max_inv_mass;
			double max_ks_daughter_pt;
			unsigned int num_ev_tau_pi_not_in_hps_pi;
			bool outputGenEvolution;
			double v0_ks_numb;
};


// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

AOD_pi0::AOD_pi0(const edm::ParameterSet& iConfig):
	//the one passed from python configure file
	IsData(iConfig.getUntrackedParameter<bool>("IsData", false)),
	OutFileName(iConfig.getUntrackedParameter<string>("OutFileName", "")),
	 // switches (collections)
	cgen(iConfig.getUntrackedParameter<bool>("GenParticles", false)),
	cbeamspot(iConfig.getUntrackedParameter<bool>("BeamSpot", false)),
	crecprimvertex(iConfig.getUntrackedParameter<bool>("RecPrimVertex", false)),
	crectrack(iConfig.getUntrackedParameter<bool>("RecTrack", false)),
	crecphoton(iConfig.getUntrackedParameter<bool>("RecPhoton", false)),
	crecpizero(iConfig.getUntrackedParameter<bool>("RecPiZero", false)),
	crecsv(iConfig.getUntrackedParameter<bool>("RecSecVertex", false)),
	crecpfjet(iConfig.getUntrackedParameter<bool>("RecJet", false)),
	crecv0(iConfig.getUntrackedParameter<bool>("RecV0", false)),
	// tracks
	cTrackPtMin(iConfig.getUntrackedParameter<double>("RecTrackPtMin", 0.5)),
	cTrackEtaMax(iConfig.getUntrackedParameter<double>("RecTrackEtaMax", 2.4)),
	cTrackDxyMax(iConfig.getUntrackedParameter<double>("RecTrackDxyMax", 2.0)),
	cTrackDzMax(iConfig.getUntrackedParameter<double>("RecTrackDzMax", 2.0)),
	cTrackNum(iConfig.getUntrackedParameter<int>("RecTrackNum", 0)),
	// photons
	cPhotonPtMin(iConfig.getUntrackedParameter<double>("RecPhotonPtMin", 1.)),
	cPhotonEtaMax(iConfig.getUntrackedParameter<double>("RecPhotonEtaMax", 2.5)),
	cPhotonNum(iConfig.getUntrackedParameter<int>("RecPhotonNum", 0)),
	// pzeros
	cPizeroPtMin(iConfig.getUntrackedParameter<double>("RecPiZeroPtMin", 30.)),
	cPizeroEtaMax(iConfig.getUntrackedParameter<double>("RecPiZeroEtaMax", 4.5)),
	cPizeroNum(iConfig.getUntrackedParameter<int>("RecPizeroNum", 0)),
	// jets
	cJetPtMin(iConfig.getUntrackedParameter<double>("RecJetPtMin", 30.)),
	cJetEtaMax(iConfig.getUntrackedParameter<double>("RecJetEtaMax", 4.5)),
	cJetNum(iConfig.getUntrackedParameter<int>("RecJetNum", 0)),
	// matched Kaons
	match_KsV0_to_HPS(iConfig.getUntrackedParameter<bool>("Match_KsV0_to_HPS", true)),
	cDZCut(iConfig.getUntrackedParameter<double>("DZCut", 999.)),
	cKtoTauDR(iConfig.getUntrackedParameter<double>("KtoTauDR", 1.)),
	matchByReference(iConfig.getUntrackedParameter<bool>("MatchByReference", true)),
	//other
	debug(iConfig.getUntrackedParameter<bool>("Debug", true)),
	mute(iConfig.getUntrackedParameter<bool>("Mute", false)),
	tkIPSigXYCut(iConfig.getUntrackedParameter<double>("tkIPSigXYCut", -1)),
	vtxDecaySigXYCut(iConfig.getUntrackedParameter<double>("vtxDecaySigXYCut", -1)),
	num_pion_res(0),
	max_inv_mass(0),
	max_ks_daughter_pt(0),
	num_ev_tau_pi_not_in_hps_pi(0),
	outputGenEvolution(true),
	v0_ks_numb(0)
{
	// TFileService : add the trees
	usesResource("TFileService");
	// create trees
	kaon_tree = fs->make<TTree>("kaon_tree", "kaon_tree");
	once_tree = fs->make<TTree>("once_tree", "once_tree");
	pions_tree = fs->make<TTree>("pions_tree", "pions_tree");
	kaon892_tree = fs->make<TTree>("kaon892_tree", "kaon892_tree");

	map_kaons[311] = "K0";
	map_kaons[310] = "K0s";
	map_kaons[130] = "K0l";
	map_kaons[130] = "K0l";
	map_kaons[313] = "K(892)0";
	map_kaons[323] = "K(892)c";
	map_pions[111] = "pi0";
	map_pions[211] = "pic";

	if (IsData) outputGenEvolution = false;

	if (outputGenEvolution)
		ofs.open("geninfo.txt", std::ofstream::out /*| std::ofstream::app*/);
	/*
		// else ofs  (cout.rdbuf());
		// if(outputGenEvolution)
		// {
		//     ofs.open ("geninfo.txt", std::ofstream::out | std::ofstream::app);
		//     buf = ofs.rdbuf();
		// } else {
		//     buf = std::cout.rdbuf();
		// }
		//std::ostream out(buf);

		// if (outputGenEvolution) {
		//   std::ofstream * os = new std::ofstream("geninfo.txt", std::ofstream::out | std::ofstream::app);
		//   ofs = os;}
		//   else
		//   {
		//     std::ostream& os = std::cout;
		//      ofs = &os;
		//   }

		//  std::ostream& stream = std::cout;
		//  if (outputGenEvolution)
		//  {

		 // ofs(this->byname);
		// }
		//    else
		//    {

		//       ofs = stream;
		//    }
	*/

	// Saved histograms - turned off
		if (OutFileName.empty())
		{
			if (!IsData)
				OutFileName = "simaod_pi0_WHY.root";
			else
				OutFileName = "aod_pi0_WHY.root";
		}

	// Tokens
		//from generalV0Candidates
		KshortCollectionTag_stand_ = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("KshortCollectionTag_stand"));
		KshortCollectionToken_ = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("KshortCollectionTag"));//cms.InputTag("generalV0Candidates","Kshort","RECO"),//vector<reco::VertexCompositeCandidate>    "generalV0Candidates"       "Kshort"          "RECO"
		LambdaCollectionToken_ = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("LambdaCollectionTag"));//cms.InputTag("generalV0Candidates","Lambda","RECO")

		PFCandidateCollectionToken_ = consumes<reco::PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("PFCandidateCollectionTag"));//cms.InputTag("particleFlow")
		TauHPSCollectionToken_ = consumes<reco::PFTauCollection>(edm::InputTag("hpsPFTauProducer","","RECO"));
		TauPiZeroCollectionToken_ = consumes<reco::RecoTauPiZeroCollection>(iConfig.getParameter<edm::InputTag>("TauPiZeroCollectionTag"));//cms.InputTag("hpsPFTauProducer","pizeros","RECO")
		if (!IsData) GenParticleCollectionToken_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("GenParticleCollectionTag"));//GenParticlesToken_ = consumes<reco::GenParticleCollection>(edm::InputTag("genParticles","","")); //typedef std::vector<GenParticle> reco::GenParticleCollection
		PVToken_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("PVCollectionTag")); //offlinePrimaryVerticesToken_ = consumes<reco::VertexCollection>(edm::InputTag("offlinePrimaryVertices","",""));
		BeamSpotToken_ = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("BeamSpotCollectionTag"));
		//SVToken_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("SecVertexCollectionTag"));
		//JetCollectionToken_ = consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("JetCollectionTag"));
		HPSTrackTagToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("HPSTrackTag"));

	// For quick output level control
		if (mute)
		{
			cout.setstate(ios_base::failbit);
			clog.setstate(ios_base::failbit);
		}
		else if (!debug) cout.setstate(ios_base::failbit);
}

//Template functions to simplify output
	void AOD_pi0::dout() { cout << endl; }
	template <typename Head, typename... Tail>
	void AOD_pi0::dout(Head H, Tail... T)
	{
		cout << H << ' ';
		dout(T...);
	}

	void AOD_pi0::dlog() { clog << endl; }
	template <typename Head, typename... Tail>
	void AOD_pi0::dlog(Head H, Tail... T)
	{
		clog << H << ' ';
		dlog(T...);
	}

//Printing on screen functions
void AOD_pi0::RecO_Cand_type(const reco::Candidate* cand)
{
	if (cand->isCaloMuon())
		dout("isCaloMuon():", cand->isCaloMuon());
	else if (cand->isConvertedPhoton())
		dout("isConvertedPhoton():", cand->isConvertedPhoton());
	else if (cand->isElectron())
		dout("isElectron():", cand->isElectron());
	else if (cand->isGlobalMuon())
		dout("isGlobalMuon():", cand->isGlobalMuon());
	else if (cand->isJet())
		dout("isJet():", cand->isJet());
	else if (cand->isMuon())
		dout("isMuon():", cand->isMuon());
	else if (cand->isPhoton())
		dout("isPhoton():", cand->isPhoton());
	else if (cand->isStandAloneMuon())
		dout("isStandAloneMuon():", cand->isStandAloneMuon());
	else if (cand->isTrackerMuon())
		dout("isTrackerMuon():", cand->isTrackerMuon());
	else
		dlog("unknown type of particles");
}

// ------------ method called once each job just before starting event loop  ------------
void AOD_pi0::beginJob()
{
	// initialize branches once per job, for types see https://root.cern.ch/doc/master/classTTree.html
		kaon_tree->Branch("nPionsInJetsWithKs", &nPionsInJetsWithKs, "nPionsInJetsWithKs/i");
		// V0 comparison to HPS
		kaon_tree->Branch("v0hps_pions_dPhi_1", &v0hps_pions_dPhi_1, "v0hps_pions_dPhi_1/D");
		kaon_tree->Branch("v0hps_pions_dEta_1", &v0hps_pions_dEta_1, "v0hps_pions_dEta_1/D");
		kaon_tree->Branch("v0hps_pions_dPhi_2", &v0hps_pions_dPhi_2, "v0hps_pions_dPhi_2/D");
		kaon_tree->Branch("v0hps_pions_dEta_2", &v0hps_pions_dEta_2, "v0hps_pions_dEta_2/D");
		kaon_tree->Branch("v0hps_KsTau_DR", &v0hps_KsTau_DR, "v0hps_KsTau_DR/D");
		kaon_tree->Branch("v0hps_KsDiPion_DR", &v0hps_KsDiPion_DR, "v0hps_KsDiPion_DR/D");
		// HPS Kaon
		kaon_tree->Branch("hps_Ks_inv_m_pi", &hps_Ks_inv_m_pi, "hps_Ks_inv_m_pi/D");
		kaon_tree->Branch("hps_pions_DR", &hps_pions_DR, "hps_pions_DR/D");
		kaon_tree->Branch("hps_deta", &hps_deta, "hps_deta/D");
		kaon_tree->Branch("hps_dphi", &hps_dphi, "hps_dphi/D");
		kaon_tree->Branch("hps_px", &hps_px, "hps_px/D");
		kaon_tree->Branch("hps_py", &hps_py, "hps_py/D");
		kaon_tree->Branch("hps_pz", &hps_pz, "hps_pz/D");
		kaon_tree->Branch("hps_energy", &hps_energy, "hps_energy/D");
		// HPS neg pion
		kaon_tree->Branch("hps_pt_1", &hps_pt_1, "hps_pt_1/D");
		kaon_tree->Branch("hps_eta_1", &hps_eta_1, "hps_eta_1/D");
		kaon_tree->Branch("hps_phi_1", &hps_phi_1, "hps_phi_1/D");
		kaon_tree->Branch("hps_px_1", &hps_px_1, "hps_px_1/D");
		kaon_tree->Branch("hps_py_1", &hps_py_1, "hps_py_1/D");
		kaon_tree->Branch("hps_pz_1", &hps_pz_1, "hps_pz_1/D");
		kaon_tree->Branch("hps_energy_1", &hps_energy_1, "hps_energy_1/D");
		// HPS pos pion
		kaon_tree->Branch("hps_pt_2", &hps_pt_2, "hps_pt_2/D");
		kaon_tree->Branch("hps_eta_2", &hps_eta_2, "hps_eta_2/D");
		kaon_tree->Branch("hps_phi_2", &hps_phi_2, "hps_phi_2/D");
		kaon_tree->Branch("hps_px_2", &hps_px_2, "hps_px_2/D");
		kaon_tree->Branch("hps_py_2", &hps_py_2, "hps_py_2/D");
		kaon_tree->Branch("hps_pz_2", &hps_pz_2, "hps_pz_2/D");
		kaon_tree->Branch("hps_energy_2", &hps_energy_2, "hps_energy_2/D");
		// V0 kaon
		kaon_tree->Branch("v0_Ks_inv_m_pi", &v0_Ks_inv_m_pi, "v0_Ks_inv_m_pi/D");
		kaon_tree->Branch("v0_Ks_pions_DR", &v0_Ks_pions_DR, "v0_Ks_pions_DR/D");
		kaon_tree->Branch("v0_deta", &v0_deta, "v0_deta/D");
		kaon_tree->Branch("v0_dphi", &v0_dphi, "v0_dphi/D");
		kaon_tree->Branch("v0_pt", &v0_pt, "v0_pt/D");
		kaon_tree->Branch("v0_eta", &v0_eta, "v0_eta/D");
		kaon_tree->Branch("v0_phi", &v0_phi, "v0_phi/D");
		kaon_tree->Branch("v0_px", &v0_px, "v0_px/D");
		kaon_tree->Branch("v0_py", &v0_py, "v0_py/D");
		kaon_tree->Branch("v0_pz", &v0_pz, "v0_pz/D");
		kaon_tree->Branch("v0_energy", &v0_energy, "v0_energy/D");
		// V0 neg pion
		kaon_tree->Branch("v0_pt_1", &v0_pt_1, "v0_pt_1/D");
		kaon_tree->Branch("v0_eta_1", &v0_eta_1, "v0_eta_1/D");
		kaon_tree->Branch("v0_phi_1", &v0_phi_1, "v0_phi_1/D");
		kaon_tree->Branch("v0_px_1", &v0_px_1, "v0_px_1/D");
		kaon_tree->Branch("v0_py_1", &v0_py_1, "v0_py_1/D");
		kaon_tree->Branch("v0_pz_1", &v0_pz_1, "v0_pz_1/D");
		kaon_tree->Branch("v0_energy_1", &v0_energy_1, "v0_energy_1/D");
		// V0 pos pion
		kaon_tree->Branch("v0_pt_2", &v0_pt_2, "v0_pt_2/D");
		kaon_tree->Branch("v0_eta_2", &v0_eta_2, "v0_eta_2/D");
		kaon_tree->Branch("v0_phi_2", &v0_phi_2, "v0_phi_2/D");
		kaon_tree->Branch("v0_px_2", &v0_px_2, "v0_px_2/D");
		kaon_tree->Branch("v0_py_2", &v0_py_2, "v0_py_2/D");
		kaon_tree->Branch("v0_pz_2", &v0_pz_2, "v0_pz_2/D");
		kaon_tree->Branch("v0_energy_2", &v0_energy_2, "v0_energy_2/D");

		pions_tree->Branch("v0_pion", &v0_pion, "TLorentzVector");
		pions_tree->Branch("tau_pion", &tau_pion, "TLorentzVector");
		pions_tree->Branch("tau_pion_besttrack", &tau_pion_besttrack, "TLorentzVector");
		pions_tree->Branch("hpsbestTrackDiff_fx", &hpsbestTrackDiff_fx, "hpsbestTrackDiff_fx/D");
		pions_tree->Branch("hpsbestTrackDiff_fy", &hpsbestTrackDiff_fy, "hpsbestTrackDiff_fy/D");
		pions_tree->Branch("hpsbestTrackDiff_fz", &hpsbestTrackDiff_fz, "hpsbestTrackDiff_fz/D");
		pions_tree->Branch("hpsbestTrackDiff_fE", &hpsbestTrackDiff_fE, "hpsbestTrackDiff_fE/D");
		pions_tree->Branch("bestTrac_fM", &bestTrac_fM, "bestTrac_fM/D");
		pions_tree->Branch("hps_fM", &hps_fM, "hps_fM/D");
		pions_tree->Branch("v0_fM", &hps_fM, "hps_fM/D");
		pions_tree->Branch("hpsv0Diff_fx", &hpsv0Diff_fx, "hpsv0Diff_fx/D");
		pions_tree->Branch("hpsv0Diff_fy", &hpsv0Diff_fy, "hpsv0Diff_fy/D");
		pions_tree->Branch("hpsv0Diff_fz", &hpsv0Diff_fz, "hpsv0Diff_fz/D");
		pions_tree->Branch("hpsv0Diff_fE", &hpsv0Diff_fE, "hpsv0Diff_fE/D");
		pions_tree->Branch("v0bestTrackDiff_fx", &v0bestTrackDiff_fx, "v0bestTrackDiff_fx/D");
		pions_tree->Branch("v0bestTrackDiff_fy", &v0bestTrackDiff_fy, "v0bestTrackDiff_fy/D");
		pions_tree->Branch("v0bestTrackDiff_fz", &v0bestTrackDiff_fz, "v0bestTrackDiff_fz/D");
		pions_tree->Branch("v0bestTrackDiff_fE", &v0bestTrackDiff_fE, "v0bestTrackDiff_fE/D");
		pions_tree->Branch("hpsv0Diff_dPhi", &hpsv0Diff_dPhi, "hpsv0Diff_dPhi/D");
		pions_tree->Branch("hpsv0Diff_dEta", &hpsv0Diff_dEta, "hpsv0Diff_dEta/D");
		pions_tree->Branch("v0bestTrack_dR", &v0bestTrack_dR, "v0bestTrack_dR/D");
		pions_tree->Branch("hpsbestTrack_dR", &hpsbestTrack_dR, "hpsbestTrack_dR/D");
		pions_tree->Branch("hpsv0_dR", &hpsv0_dR, "hpsv0_dR/D");

		once_tree->Branch("v0_Ks_count", &v0_Ks_count, "v0_Ks_count/I");
		once_tree->Branch("primvertex_count", &primvertex_count, "primvertex_count/i");
		once_tree->Branch("goodprimvertex_count", &goodprimvertex_count, "goodprimvertex_count/i");
		once_tree->Branch("numOfUnmatchedKaons", &numOfUnmatchedKaons, "numOfUnmatchedKaons/i");
		once_tree->Branch("v_v0_Ks_inv_m_pi", &v_v0_Ks_inv_m_pi);
		once_tree->Branch("v_v0_pions_DR", &v_v0_pions_DR);
		once_tree->Branch("v_v0_pt_1", &v_v0_pt_1);
		once_tree->Branch("v_v0_pt_2", &v_v0_pt_2);
		once_tree->Branch("v_v0_eta_1", &v_v0_eta_1);
		once_tree->Branch("v_v0_eta_2", &v_v0_eta_2);
		once_tree->Branch("v_v0_matched_pt_1", &v_v0_matched_pt_1);
		once_tree->Branch("v_v0_matched_pt_2", &v_v0_matched_pt_2);
		once_tree->Branch("v_v0_matched_eta_1", &v_v0_matched_eta_1);
		once_tree->Branch("v_v0_matched_eta_2", &v_v0_matched_eta_2);
		once_tree->Branch("v_KsCombinatoricMass", &v_KsCombinatoricMass);
		once_tree->Branch("v_KsCombinatoricDR", &v_KsCombinatoricDR);
		once_tree->Branch("v_hps_Ks_inv_m_pi", &v_hps_Ks_inv_m_pi);
		once_tree->Branch("v_hps_Ks_DR", &v_hps_Ks_DR);
		once_tree->Branch("v_hps_Ks_pions_DR", &v_hps_Ks_pions_DR);
		once_tree->Branch("v_hps_pt_1", &v_hps_pt_1);
		once_tree->Branch("v_hps_pt_2", &v_hps_pt_2);
		once_tree->Branch("v_hps_eta_1", &v_hps_eta_1);
		once_tree->Branch("v_hps_eta_2", &v_hps_eta_2);
		once_tree->Branch("v_nPionsInJetsWithKs", &v_nPionsInJetsWithKs);

		//K892
		kaon892_tree->Branch("lv_K892", &lv_K892, "TLorentzVector");
		kaon892_tree->Branch("K892_fx", &K892_fx, "K892_fx/D");
		kaon892_tree->Branch("K892_fy", &K892_fy, "K892_fy/D");
		kaon892_tree->Branch("K892_fz", &K892_fz, "K892_fz/D");
		kaon892_tree->Branch("K892_fE", &K892_fE, "K892_fE/D");
		kaon892_tree->Branch("K892_fM", &K892_fM, "K892_fM/D");

		kaon892_tree->Branch("lv_K892_K0", &lv_K892_K0, "TLorentzVector");
		kaon892_tree->Branch("K892_K0_fx", &K892_K0_fx, "K892_K0_fx/D");
		kaon892_tree->Branch("K892_K0_fy", &K892_K0_fy, "K892_K0_fy/D");
		kaon892_tree->Branch("K892_K0_fz", &K892_K0_fz, "K892_K0_fz/D");
		kaon892_tree->Branch("K892_K0_fE", &K892_K0_fE, "K892_K0_fE/D");
		kaon892_tree->Branch("K892_K0_dxy", &K892_K0_dxy, "K892_K0_dxy/D");
		kaon892_tree->Branch("K892_K0_dz", &K892_K0_dz, "K892_K0_dz/D");
		kaon892_tree->Branch("K892_K0_fM", &K892_K0_fM, "K892_K0_fM/D");

		kaon892_tree->Branch("lv_K892_pi0", &lv_K892_pi0, "TLorentzVector");
		kaon892_tree->Branch("K892_pi0_fx", &K892_pi0_fx, "K892_pi0_fx/D");
		kaon892_tree->Branch("K892_pi0_fy", &K892_pi0_fy, "K892_pi0_fy/D");
		kaon892_tree->Branch("K892_pi0_fz", &K892_pi0_fz, "K892_pi0_fz/D");
		kaon892_tree->Branch("K892_pi0_fE", &K892_pi0_fE, "K892_pi0_fE/D");
		kaon892_tree->Branch("K892_pi0_fM", &K892_pi0_fM, "K892_pi0_fM/D");


		kaon892_tree->Branch("n_K892", &n_K892, "n_K892/I");

}

// ------------ method called once each job just after ending the event loop  ------------
void AOD_pi0::endJob(){}

AOD_pi0::~AOD_pi0()
{
	dlog("Total num of not-matched events:", num_ev_tau_pi_not_in_hps_pi);
	dlog("num_pion_res: ", num_pion_res);
	dlog("max_inv_mass: ", max_inv_mass);
	dlog("max_ks_daughter_pt: ", max_ks_daughter_pt);
	dlog("v0_ks_numb: ", v0_ks_numb);
	if (outputGenEvolution) ofs.close();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void AOD_pi0::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
	//The following says we do not know what parameters are allowed so do no validation
	//Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

#endif
