#include <memory>
#include <map>
#include <fstream>      // std::ofstream
#include <iostream>
#include <set>
#include <vector>
#include <string>
#include <utility>

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
  // #include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
  // #include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
//error using namespace reco;

#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/TauReco/interface/RecoTauPiZero.h"
#include "DataFormats/TauReco/interface/RecoTauPiZeroFwd.h"

// Alex
#include "DataFormats/PatCandidates/interface/Tau.h"

//new
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "PhysicsTools/HepMCCandAlgos/plugins/MCTruthDeltaRMatcherNew.cc"



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
  #include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
  #include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
  #include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

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
  #include "TTree.h"
  #include "TLorentzVector.h"
  #include <Math/Functions.h>
  #include <Math/SVector.h>
  #include <Math/SMatrix.h>

  #include "FWCore/ServiceRegistry/interface/Service.h"
  #include "CommonTools/UtilAlgos/interface/TFileService.h"

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
class AOD_pi0 : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
	 
  public:

    explicit  AOD_pi0(const edm::ParameterSet&);
    ~AOD_pi0();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    static void RecO_Cand_type(const reco::Candidate* cand);
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
    
  
    TFile* outfile;
    TDirectory* hist_directory[4];

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Histograms
      TH1D* pions_inv_m;
      TH1D* num_pions;
      TH1D* taus_isol_pi0_inv_m_to_ks;
      TH1D* taus_isol_pi0_inv_pt;
      TH1D* taus_pi0_inv_m_to_ks;
      TH1D* taus_pi0_inv_pt;
      TH1D* taus_pi_charged_inv_m_to_ks ;
      TH1D* taus_pi_charged_inv_pt;
      TH1D* taus_pi0_had_inv_m_to_ks ;
      TH1D* taus_pi0_had_inv_pt;
      TH1D* h_gen_k0_all_to_pi0;
      TH1D* h_gen_k0_all_to_pic;

      // V0 collection
      TH1D* h_Ks_v0_count;
      TH1D* h_Ks_v0_daughter_pt;
      TH1D* h_Ks_v0_inv_m_pi;
      TH1D* h_Ks_v0_number_per_event;
      TH1D* h_Ks_v0_vx;
      TH1D* h_Ks_v0_vy;
      TH1D* h_Ks_v0_vz;
      TH1D* h_Ks_v0_found_in_hps_tau;
      TH1D* h_Ks_v0_found_in_hps_tau_dR;
      TH1D* h_Ks_v0_found_in_hps_tau_dRcut;
      TH1D* h_Ks_v0_found_in_hps_tau_dR_only_one_pion_left;
      TH1D* h_Ks_v0_found_in_hps_tau_m_inv;

      TH1D* h_Ks_v0_pions_dR;
      TH1D* h_Ks_v0_hps_pions_combinatoric_dR;
      TH1D* h_Ks_v0_pions_and_hps_pions_combined_dR;

      TH1D* h_Ks_v0_n_ev_passing_dz_cut;
      TH1D* h_Ks_v0_n_Ks_in_jets_per_event;
      TH1D* h_Ks_v0_n_tau_jets_per_event;
      TH1D* h_Ks_v0_n_pion_in_tau_jets;
      TH1D* h_Ks_v0_n_pion_in_tau_jets_with_good_Ks;

      

      TH1D* h_Ks_v0_pions_and_hps_pions_combined_dR_with_no_constrain;
      std::vector<TH1D*> map_Ks_v0_histos;

      TH1D* h_K892;

      TH1D* h_K892_0_gen_number_per_event;
      TH1D* h_K892_0_gen_vx;
      TH1D* h_K892_0_gen_vy;
      TH1D* h_K892_0_gen_vz;

      TH1D* h_K892_c_gen_number_per_event;
      TH1D* h_K892_c_gen_vx;
      TH1D* h_K892_c_gen_vy;
      TH1D* h_K892_c_gen_vz;

		// primary vertex
		TH1D* h_primvertex_count;
		TH1D* h_goodprimvertex_count;
		TH1D* h_primvertex_x;
		TH1D* h_primvertex_y;
		TH1D* h_primvertex_z;
		TH1D* h_primvertex_chi2;
		TH1D* h_primvertex_ndof;
		TH1D* h_primvertex_ptq;
		TH1D* h_primvertex_ntracks;
		TH1D* h_primvertex_cov_x;
		TH1D* h_primvertex_cov_y;
		TH1D* h_primvertex_cov_z;
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tokens for the Collections
      edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> KshortCollectionToken_;
      edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> LambdaCollectionToken_;
      edm::EDGetTokenT<reco::PFCandidateCollection> PFCandidateCollectionToken_;
      edm::EDGetTokenT<reco::RecoTauPiZeroCollection> TauPiZeroCollectionToken_;//hpsPFTauProducer
      edm::EDGetTokenT<reco::PFTauCollection> TauHPSCollectionToken_;
      // vector<int>                           "genParticles"              ""                "HLT"   
      // vector<reco::GenJet>                  "ak4GenJets"                ""                "HLT"   
      // vector<reco::GenJet>                  "ak4GenJetsNoNu"            ""                "HLT"   
      // vector<reco::GenJet>                  "ak8GenJets"                ""                "HLT"   
      // vector<reco::GenJet>                  "ak8GenJetsNoNu"            ""                "HLT"   
      // vector<reco::GenMET>                  "genMetCalo"                ""                "HLT"   
      // vector<reco::GenMET>                  "genMetTrue"                ""                "HLT"   
      // vector<reco::GenParticle>             "genParticles"              ""   
      //edm::EDGetTokenT<reco::GenParticleCollection>  GenParticlesToken_;
     // edm::EDGetTokenT<reco::VertexCollection>  offlinePrimaryVerticesToken_ ;
      edm::EDGetTokenT<reco::GenParticleCollection> GenParticleCollectionToken_;
	  edm::EDGetTokenT<reco::VertexCollection> PVToken_;

    //P0's
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

    math::XYZPoint pv_position;
	reco::Vertex primvertex;
     // primary vertex
	  UInt_t  primvertex_count;
	  UInt_t  goodprimvertex_count;
	  Float_t primvertex_x;
	  Float_t primvertex_y;
	  Float_t primvertex_z;
	  Float_t primvertex_chi2;
	  Float_t primvertex_ndof;
	  Float_t primvertex_ptq;
	  Int_t   primvertex_ntracks;
	  Float_t primvertex_cov[6];

    Int_t v0_count; // V0s - Ks or Lambdas

    vector< vector <const reco::Candidate *>> v_daughters_k0s_to_pi0;
    vector< vector <const reco::Candidate *>> v_daughters_k0l_to_pi0;
    vector< vector <const reco::Candidate *>> v_daughters_k0_to_pi0;
    vector< vector <const reco::Candidate *>> v_daughters_k0s_to_pic;
    vector< vector <const reco::Candidate *>> v_daughters_k0l_to_pic;
    vector <reco::RecoTauPiZero *> v_strips_ref;
    vector < vector< vector <reco::RecoTauPiZero *>>> v_daughters_k0s_to_pi0_SimToStrips;
    vector < vector< const reco::Candidate *>> v_daughters_k0s_to_pi0_StipsToSim;


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

		// other
	    bool debug;
	    bool mute;
		bool match_KsV0_to_HPS;

    int num_pion_res;
    double max_inv_mass;
    double max_ks_daughter_pt;
    unsigned int num_ev_tau_pi_not_in_hps_pi;

 //    std::streambuf * buf;
	std::ofstream ofs;
	
    //std::ostream ofs;//ostream ofstream

    // std::ofstream of;
    // std::ostream & ofs = (!outputGenEvolution) ? cout : of.open("geninfo.txt", std::ofstream::out /*| std::ofstream::app*/);

    // std::ostream* ofs;
    bool outputGenEvolution;
    // std::ofstream byname;
    // std::ostream& ofs;


};


// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

AOD_pi0::AOD_pi0(const edm::ParameterSet& iConfig):
  //the one passed from python configure file
  IsData(iConfig.getUntrackedParameter<bool>("IsData", false)),
  OutFileName(iConfig.getUntrackedParameter<string>("OutFileName", "aod_pi0.root")),
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
  //other
  debug(iConfig.getUntrackedParameter<bool>("Debug", true)),
  mute(iConfig.getUntrackedParameter<bool>("Mute", false)),
  match_KsV0_to_HPS(iConfig.getUntrackedParameter<bool>("Match_KsV0_to_HPS", true)),
  num_pion_res(0),
  max_inv_mass(0),
  max_ks_daughter_pt(0),
  num_ev_tau_pi_not_in_hps_pi(0),
  outputGenEvolution(true)
{
  //usesResource("TFileService");
  map_kaons[311] = "K0";
  map_kaons[310] = "K0s";
  map_kaons[130] = "K0l";
  map_kaons[130] = "K0l";
  map_kaons[313] = "K(892)0";
  map_kaons[323] = "K(892)c";

  map_pions[111] = "pi0";
  map_pions[211] = "pic";

  
  //map_K892_c_gen_histos[323] = {h_K892_c_gen_number_per_event, h_K892_c_gen_vx, h_K892_c_gen_vy, h_K892_c_gen_vz};

  if (outputGenEvolution)  
  {
  	ofs.open("geninfo.txt", std::ofstream::out /*| std::ofstream::app*/);
  }

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
 //   	std::ostream& os = std::cout;
 //   	 ofs = &os;
 //   }

// 	std::ostream& stream = std::cout;
// 	if (outputGenEvolution) 
// 	{
	  
	   // ofs(this->byname);
// }
//    else 
//    {
   
//    	 ofs = stream;
//    }
   

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Saved histograms
    if (!IsData) OutFileName = "simaod_pi0_WHY.root";
    else outputGenEvolution = false;
    outfile = new TFile((TString)OutFileName,"RECREATE");

    outfile->cd();
      pions_inv_m  = new TH1D("pions_inv_m","inv mass of pions in taus", 100, 0, 1);
      num_pions = new TH1D("num_of_pios","num of pions", 10, 0, 10);

    hist_directory[1]  = outfile->mkdir("ks_coll", "ks_collection");
    //hist_directory[1]->cd();  //= outfile->mkdir("ks_coll", "ks_collection");
      h_Ks_v0_count = new TH1D("h_Ks_v0_count","v0 count", 10, 0, 10);
      h_Ks_v0_daughter_pt = new TH1D("h_Ks_v0_daughter_pt","ks daughters pt", 1800, 0, 180);
      h_Ks_v0_inv_m_pi = new TH1D("h_Ks_v0_inv_m_pi","ks daughters inv mass", 1000, 0, 10);
      h_Ks_v0_number_per_event = new TH1D("h_Ks_v0_number_per_event","ks from V0 coll, NPE that passed PV", 10, 0, 10);
      h_Ks_v0_vx = new TH1D("h_Ks_v0_vx","ks from V0 coll, x position", 1000, 0, 10);
      h_Ks_v0_vy = new TH1D("h_Ks_v0_vy","ks from V0 coll, y position", 1000, 0, 10);
      h_Ks_v0_vz = new TH1D("h_Ks_v0_vz","ks from V0 coll, z position", 1000, 0, 10);
      h_Ks_v0_found_in_hps_tau = new TH1D("h_Ks_v0_found_in_hps_tau","ks from V0 coll, NPE that passed PV and found in HPS tau jets", 1000, 0, 10);
      h_Ks_v0_found_in_hps_tau_dR = new TH1D("h_Ks_v0_found_in_hps_tau_dR","ks from V0 coll, dR of ks and tau jet, that passed PV and found in HPS tau jets and pions are matched", 1000, 0, 10);
      h_Ks_v0_found_in_hps_tau_dRcut = new TH1D("h_Ks_v0_found_in_hps_tau_dRcut","ks from V0 coll, dR of ks and tau jet, that passed PV and found in HPS tau jets and cut on dR < 0.5", 1000, 0, 10);
      h_Ks_v0_found_in_hps_tau_dR_only_one_pion_left = new TH1D("h_Ks_v0_found_in_hps_tau_dR_only_one_pion_left","ks from V0 coll, dR of ks and tau jet, that passed PV and found in HPS tau jets and only 1 pion is matched", 1000, 0, 10);
      h_Ks_v0_found_in_hps_tau_m_inv = new TH1D("h_Ks_v0_found_in_hps_tau_m_inv","matched tau pions with ks from v0, daughters inv mass", 1000, 0, 10);

      h_Ks_v0_pions_dR = new TH1D("h_Ks_v0_pions_dR","ks from V0 coll, dR for pions of Ks", 1000, 0, 10);
	  h_Ks_v0_hps_pions_combinatoric_dR = new TH1D("h_Ks_v0_hps_pions_combinatoric_dR","combinatoric pions dR for pions of HPS taus jets", 1000, 0, 10);
	  h_Ks_v0_pions_and_hps_pions_combined_dR = new TH1D("h_Ks_v0_pions_and_hps_pions_combined_dR","dr of two pions of KSv V0 and all the pions of the hps tau jet with the kaon", 1000, 0, 10);

	  h_Ks_v0_n_ev_passing_dz_cut = new TH1D("h_Ks_v0_n_ev_passing_dz_cut", "N ev with at least 1 Ks passing dz cut", 10, 0, 10);
      h_Ks_v0_n_Ks_in_jets_per_event = new TH1D("h_Ks_v0_n_Ks_in_jets_per_event","NPE of Ks matched with any tau jet", 10, 0, 10);
      h_Ks_v0_n_tau_jets_per_event = new TH1D("h_Ks_v0_n_tau_jets_per_event","NPE tau jets", 100, 0, 100);
      h_Ks_v0_n_pion_in_tau_jets = new TH1D("h_Ks_v0_n_pion_in_tau_jets","num of pic per jet", 100, 0, 100);
      h_Ks_v0_n_pion_in_tau_jets_with_good_Ks = new TH1D("h_Ks_v0_n_pion_in_tau_jets_with_good_Ks","num of pic per jet with a Ks", 100, 0, 100);

      h_Ks_v0_pions_and_hps_pions_combined_dR_with_no_constrain = new TH1D("h_Ks_v0_pions_and_hps_pions_combined_dR_with_no_constrain", "dr of two pions of KSv V0 and all the pions of all the hps tau jets, no constrains", 1000, 0, 10);
    hist_directory[0]  = outfile->mkdir("Taus_pions_coll", "Taus_pions_collections");
    //hist_directory[0]->cd();  //= outfile->mkdir("Taus_pions_coll", "Taus_pions_collections");
      taus_isol_pi0_inv_m_to_ks = new TH1D("taus_isol_pi0_inv_m_to_ks","all Pairs of tau isolation pions inv mass", 1000, 0, 17);
      taus_isol_pi0_inv_pt = new TH1D("taus_isol_pi0_inv_pt","all Pairs of tau isolation pions int pt", 1000, 0, 10);
      taus_pi0_inv_m_to_ks = new TH1D("taus_pi0_inv_m_to_ks","all Pairs of tau pions inv mass", 1000, 0, 17);
      taus_pi0_inv_pt = new TH1D("taus_pi0_inv_pt","all Pairs of tau pions inv pt", 1000, 0, 10);
    hist_directory[2]  = outfile->mkdir("Taus_charged_had_coll", "Taus_charged_had_coll");
    //hist_directory[2]->cd();  //= outfile->mkdir("Taus_charged_had_coll", "Taus_charged_had_coll");
      taus_pi_charged_inv_m_to_ks = new TH1D("taus_pi_charged_inv_m_to_ks","all Pairs of tau pions from Charged had coll inv mass", 1000, 0, 17);
      taus_pi_charged_inv_pt = new TH1D("taus_pi_charged_inv_pt","all Pairs of tau pions from Charged had coll inv pt", 1000, 0, 10);
    hist_directory[3]  = outfile->mkdir("Taus_neutral_had_coll", "Taus_neutral_had_coll");
    //hist_directory[3]->cd();  //= outfile->mkdir("Taus_neutral_had_coll", "Taus_neutral_had_coll");
      taus_pi0_had_inv_m_to_ks = new TH1D("taus_pi0_had_inv_m_to_ks", "all Pairs of tau pions from Neutal had coll inv mass", 1000, 0, 17);
      taus_pi0_had_inv_pt = new TH1D("taus_pi0_had_inv_pt", "all Pairs of tau pions from Neutal had coll inv pt", 1000, 0, 10);
    h_gen_k0_all_to_pi0 = new TH1D("h_gen_k0_all_to_pi0", "gen. all k0's decaying to pi0", 1000, 0, 10);
    h_gen_k0_all_to_pic = new TH1D("h_gen_k0_all_to_pic", "gen. all k0's decaying to pi+-", 1000, 0, 10);

    h_K892 = new TH1D("h_K892", "m_vis for K(892)", 1000, 0, 10);
    h_K892_0_gen_number_per_event = new TH1D("h_K892_0_gen_number_per_event", "number of generated K(892) neutral", 1000, 0, 1000);
    h_K892_0_gen_vx = new TH1D("h_K892_0_gen_vx", "vx of generated K(892) neutral", 1000, -5, 5);
	h_K892_0_gen_vy = new TH1D("h_K892_0_gen_vy", "vy of generated K(892) neutral", 1000, -5, 5);
	h_K892_0_gen_vz = new TH1D("h_K892_0_gen_vz", "vz of generated K(892) neutral", 1000, -5, 5);

    h_K892_c_gen_number_per_event = new TH1D("h_K892_c_gen_number_per_event", "number of generated K(892) charged", 1000, 0, 1000);
    h_K892_c_gen_vx = new TH1D("h_K892_c_gen_vx", "vx of generated K(892) charged", 1000, -5, 5);
	h_K892_c_gen_vy = new TH1D("h_K892_c_gen_vy", "vy of generated K(892) charged", 1000, -5, 5);
	h_K892_c_gen_vz = new TH1D("h_K892_c_gen_vz", "vz of generated K(892) charged", 1000, -5, 5);

	h_primvertex_count = new TH1D("h_primvertex_count", "h_primvertex_count", 100, 0, 100);
	h_goodprimvertex_count = new TH1D("h_goodprimvertex_count", "h_goodprimvertex_count", 100, 0, 100);
	h_primvertex_x = new TH1D("h_primvertex_x", "h_primvertex_x", 1000, -5, 5);
	h_primvertex_y = new TH1D("h_primvertex_y", "h_primvertex_y", 1000, -5, 5);
	h_primvertex_z = new TH1D("h_primvertex_z", "h_primvertex_z", 1000, -30, 30);
	h_primvertex_chi2 = new TH1D("h_primvertex_chi2", "h_primvertex_chi2", 1000, 0, 100);
	h_primvertex_ndof = new TH1D("h_primvertex_ndof", "h_primvertex_ndof", 1000, 0, 1000);
	h_primvertex_ptq = new TH1D("h_primvertex_ptq", "h_primvertex_ptq", 1000, 0, 1000);
	h_primvertex_ntracks = new TH1D("h_primvertex_ntracks", "h_primvertex_ntracks", 1000, 0, 200);
	h_primvertex_cov_x = new TH1D("h_primvertex_cov_x", "h_primvertex_cov_x", 1000, -5, 5);
	h_primvertex_cov_y = new TH1D("h_primvertex_cov_y", "h_primvertex_cov_y", 1000, -5, 5);
	h_primvertex_cov_z = new TH1D("h_primvertex_cov_z", "h_primvertex_cov_z", 1000, -5, 5);
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAP of histos - only after initialisation!!
		map_K892_gen_histos[313].push_back(h_K892_0_gen_number_per_event);
		map_K892_gen_histos[313].push_back(h_K892_0_gen_vx);
		map_K892_gen_histos[313].push_back(h_K892_0_gen_vy);
		map_K892_gen_histos[313].push_back(h_K892_0_gen_vz);

		map_K892_gen_histos[323].push_back(h_K892_c_gen_number_per_event);
		map_K892_gen_histos[323].push_back(h_K892_c_gen_vx);
		map_K892_gen_histos[323].push_back(h_K892_c_gen_vy);
		map_K892_gen_histos[323].push_back(h_K892_c_gen_vz);


	      TH1D* temp[21] = {h_Ks_v0_number_per_event, h_Ks_v0_vx, h_Ks_v0_vy, h_Ks_v0_vz, h_Ks_v0_count, //0 - 5
	      					h_Ks_v0_daughter_pt, h_Ks_v0_inv_m_pi, h_Ks_v0_found_in_hps_tau, h_Ks_v0_found_in_hps_tau_dR, h_Ks_v0_found_in_hps_tau_dRcut, //6 -10
	      					h_Ks_v0_found_in_hps_tau_dR_only_one_pion_left, h_Ks_v0_found_in_hps_tau_m_inv, h_Ks_v0_pions_dR, h_Ks_v0_hps_pions_combinatoric_dR, h_Ks_v0_pions_and_hps_pions_combined_dR, //11 -15
	      					h_Ks_v0_n_ev_passing_dz_cut, h_Ks_v0_n_Ks_in_jets_per_event, h_Ks_v0_n_tau_jets_per_event, h_Ks_v0_n_pion_in_tau_jets, h_Ks_v0_n_pion_in_tau_jets_with_good_Ks,
	      					h_Ks_v0_pions_and_hps_pions_combined_dR_with_no_constrain};//15-20
	      map_Ks_v0_histos.insert(map_Ks_v0_histos.end(), temp, temp + (sizeof(temp)/sizeof(temp[0])));
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tokens

	//RASP
	    //from generalV0Candidates 
	      KshortCollectionToken_ = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("KshortCollectionTag"));//cms.InputTag("generalV0Candidates","Kshort","RECO"),//vector<reco::VertexCompositeCandidate>    "generalV0Candidates"       "Kshort"          "RECO"
	      LambdaCollectionToken_ = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("LambdaCollectionTag"));//cms.InputTag("generalV0Candidates","Lambda","RECO"),
   		PFCandidateCollectionToken_ = consumes<reco::PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("PFCandidateCollectionTag"));
   		//from hpsPFTauProducer
	     //Pi0 
	      TauPiZeroCollectionToken_ = consumes<reco::RecoTauPiZeroCollection>(iConfig.getParameter<edm::InputTag>("TauPiZeroCollectionTag"));//cms.InputTag("hpsPFTauProducer","pizeros","RECO")
		//SIMAOD
   			if (!IsData) GenParticleCollectionToken_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("GenParticleCollectionTag"));//GenParticlesToken_ = consumes<reco::GenParticleCollection>(edm::InputTag("genParticles","","")); //typedef std::vector<GenParticle> reco::GenParticleCollection
		PVToken_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("PVCollectionTag")); //offlinePrimaryVerticesToken_ = consumes<reco::VertexCollection>(edm::InputTag("offlinePrimaryVertices","",""));
		/*
		  SVToken_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("SecVertexCollectionTag"));
		  JetCollectionToken_ = consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("JetCollectionTag"));
		  BeamSpotToken_ = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("BeamSpotCollectionTag"));
		*/
	//New
    	//Taus
      		TauHPSCollectionToken_ = consumes<reco::PFTauCollection>(edm::InputTag("hpsPFTauProducer","","RECO"));
    //Primary vertex
    	



  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% For quick output level control
    if (mute)
    {
      cout.setstate(ios_base::failbit);
      clog.setstate(ios_base::failbit);
    }
    else if (!debug) cout.setstate(ios_base::failbit);
}


//Template functions to simplify output
  void AOD_pi0::dout()
  {
      cout << endl;
  }
  template <typename Head, typename... Tail>
  void AOD_pi0::dout(Head H, Tail... T)
  {
    cout << H << ' ';
    dout(T...);
  }


  void AOD_pi0::dlog()
  {
      clog << endl;
  }
  template <typename Head, typename... Tail>
  void AOD_pi0::dlog(Head H, Tail... T)
  {
    clog << H << ' ';
    dlog(T...);
  }

//Printing functions
void AOD_pi0::RecO_Cand_type(const reco::Candidate* cand)
{
  if      (cand->isCaloMuon()) dout("isCaloMuon():", cand->isCaloMuon());
  else if (cand->isConvertedPhoton()) dout("isConvertedPhoton():", cand->isConvertedPhoton());
  else if (cand->isElectron()) dout("isElectron():", cand->isElectron());
  else if (cand->isGlobalMuon()) dout("isGlobalMuon():", cand->isGlobalMuon());
  else if (cand->isJet()) dout("isJet():", cand->isJet());
  else if (cand->isMuon()) dout("isMuon():", cand->isMuon());
  else if (cand->isPhoton()) dout("isPhoton():", cand->isPhoton());
  else if (cand->isStandAloneMuon()) dout("isStandAloneMuon():", cand->isStandAloneMuon());
  else if (cand->isTrackerMuon()) dout("isTrackerMuon():", cand->isTrackerMuon());
  else dlog("unknown type of particles");
}

 
// ------------ method called once each job just before starting event loop  ------------
void
AOD_pi0::beginJob()
{
  
}

// ------------ method called once each job just after ending the event loop  ------------
void
AOD_pi0::endJob()
{
	
	for(std::vector<TH1D*>::iterator itv = map_Ks_v0_histos.begin(); itv != map_Ks_v0_histos.end(); ++itv) 
	    (*itv)->Write();
	//h_Ks_v0_count->Write(); 
		h_Ks_v0_count->Print();
	pions_inv_m->Write();
	num_pions->Write();
	taus_isol_pi0_inv_m_to_ks->Write();
	//h_Ks_v0_daughter_pt->Write();
	//h_Ks_v0_inv_m_pi->Write();
	taus_isol_pi0_inv_pt->Write();
	taus_pi0_inv_m_to_ks->Write();
	taus_pi0_inv_pt->Write();
	taus_pi_charged_inv_m_to_ks->Write();
	taus_pi_charged_inv_pt->Write();
	h_gen_k0_all_to_pi0->Write();
	h_gen_k0_all_to_pic->Write();
	h_K892->Write();
	// h_K892_0_gen_number_per_event->Write();
	// h_K892_c_gen_number_per_event->Write();

	for ( std::map<long, std::vector<TH1D *>>::iterator it = map_K892_gen_histos.begin(); it != map_K892_gen_histos.end(); it++ )
	{
	    std::vector<TH1D*> v = it->second;
	    for(std::vector<TH1D*>::iterator itv = v.begin(); itv != v.end(); ++itv) 
	    (*itv)->Write();
	}

	h_primvertex_count->Write();
	h_goodprimvertex_count->Write();
	h_primvertex_x->Write();
	h_primvertex_y->Write();
	h_primvertex_z->Write();
	h_primvertex_chi2->Write();
	h_primvertex_ndof->Write();
	h_primvertex_ptq->Write();
	h_primvertex_ntracks->Write();
	h_primvertex_cov_x->Write();
	h_primvertex_cov_y->Write();
	h_primvertex_cov_z->Write();
	
	outfile->Close();
}

AOD_pi0::~AOD_pi0()
{
  dlog("Total num of not-matched events:", num_ev_tau_pi_not_in_hps_pi);
  dlog("num_pion_res: ", num_pion_res);
  dlog("max_inv_mass: ", max_inv_mass);
  dlog("max_ks_daughter_pt: ", max_ks_daughter_pt);

   if (outputGenEvolution) ofs.close();//ofs.close();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
AOD_pi0::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
