// -*- C++ -*-
//
// Package:    AOD_pi0_study/AOD_pi0
// Class:      AOD_pi0
// 
/**\class AOD_pi0 AOD_pi0.cc AOD_pi0_study/AOD_pi0/plugins/AOD_pi0.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/


// system include files
#include <memory>
#include <map>
#include <fstream>      // std::ofstream

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

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class AOD_pi0 : public edm::one::EDAnalyzer<edm::one::SharedResources>  
{
  public:
    explicit AOD_pi0(const edm::ParameterSet&);
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
    void BuildTo(edm::Handle<std::vector<reco::GenParticle> >& genPart, TH1 *hist, vector< vector <const reco::Candidate *>> &daughters, int mom_pdgid, map <long, string>& , int num_daugh, int daugh_pdgid, map <long, string>& );
    void BuildTo(const reco::Candidate *mom, TH1 *hist, vector< vector <const reco::Candidate *>> &daughters, int mom_pdgid, map <long, string>& , int num_daugh, int daugh_pdgid, map <long, string>& );
    void FindFinalPrt(vector<const reco::Candidate *> *d, const reco::Candidate *mom , int num_of_final_daughters, int daug_pdgid);
    bool NoRadiating(const reco::Candidate * mom, int daug_pdgid);
    bool NoRadiating(const reco::Candidate * mom, map <long, string>& map_daug_pdgid);
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
    std::ofstream ofs;
    
    TFile* outfile;
    TDirectory* hist_directory[4];

    // Histograms
      TH1D* h_v0_count;
      TH1D* pions_inv_m;
      TH1D* num_pions;
      TH1D* taus_isol_pi0_inv_m_to_ks;
      TH1D* taus_isol_pi0_inv_pt;
      TH1D* taus_pi0_inv_m_to_ks;
      TH1D* taus_pi0_inv_pt;
      TH1D* ks_daughter_pt;
      TH1D* ks_inv_m_pi;
      TH1D* taus_pi_charged_inv_m_to_ks ;
      TH1D* taus_pi_charged_inv_pt;
      TH1D* taus_pi0_had_inv_m_to_ks ;
      TH1D* taus_pi0_had_inv_pt;
      TH1D* h_gen_k0_all_to_pi0;

    // Tokens for the Collections 
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
      edm::EDGetTokenT<reco::GenParticleCollection>  GenParticlesToken_;

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


    Int_t v0_count; // V0s - Ks or Lambdas

    vector< vector <const reco::Candidate *>> v_daughters_k0s_to_pi0;
    vector< vector <const reco::Candidate *>> v_daughters_k0l_to_pi0;
    vector< vector <const reco::Candidate *>> v_daughters_k0_to_pi0;
    /*static*/ std::map <long, string> map_kaons;
    /*static*/ std::map <long, string> map_pions;
    //Parameters initialised by default
    bool IsData;
    TString OutFileName;
    bool crecpizero;
    bool debug;
    bool mute;
    int num_pion_res;
    double max_inv_mass;
    double max_ks_daughter_pt;
    unsigned int num_ev_tau_pi_not_in_hps_pi;
    bool outputGenEvolution;
};


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

// constants, enums and typedefs
//

// static data member definitions
//
  
  
  

AOD_pi0::AOD_pi0(const edm::ParameterSet& iConfig):
  //the one passed from python configure file
  IsData(iConfig.getUntrackedParameter<bool>("IsData", false)),
  OutFileName("aod_pi0.root"),
  //other
  crecpizero(iConfig.getUntrackedParameter<bool>("RecPiZero", false)),
  debug(true),
  mute(false),
  num_pion_res(0),
  max_inv_mass(0),
  max_ks_daughter_pt(0),
  num_ev_tau_pi_not_in_hps_pi(0),
  outputGenEvolution(false)
{
  usesResource("TFileService");
  map_kaons[311] = "K0";
  map_kaons[310] = "K0s";
  map_kaons[130] = "K0l";

  map_pions[111] = "pi0";
  map_pions[211] = "pic";
  ofs.open ("geninfo.txt", std::ofstream::out /*| std::ofstream::app*/);
  // Saved histograms 
    if (!IsData) OutFileName = "simaod_pi0.root";
    else outputGenEvolution = false;
    outfile = new TFile(OutFileName,"RECREATE");
    outfile->cd();  
      h_v0_count = new TH1D("v0_count","v0 count", 10, 0, 9);
      pions_inv_m  = new TH1D("pions_inv_m","inv mass of pions in taus", 100, 0, 1);
      num_pions = new TH1D("num_of_pios","num of pions", 10, 0, 9);
    hist_directory[1]  = outfile->mkdir("ks_coll", "ks_collection");
    //hist_directory[1]->cd();  //= outfile->mkdir("ks_coll", "ks_collection");
      ks_daughter_pt = new TH1D("ks_daughter_pt","ks daughters pt", 1000, 0, 10);
      ks_inv_m_pi = new TH1D("ks_inv_m_pi","ks daughters inv mass", 1000, 0, 10);
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
    h_gen_k0_all_to_pi0 = new TH1D("h_gen_k0_all_to_pi0", "gen. all k0's decaying to p0", 1000, 0, 10);
  
  // Tokens
    //Ks's
    KshortCollectionToken_ = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("KshortCollectionTag"));//vector<reco::VertexCompositeCandidate>    "generalV0Candidates"       "Kshort"          "RECO"  
    LambdaCollectionToken_ = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("LambdaCollectionTag"));
    //Pi0
    TauPiZeroCollectionToken_ = consumes<reco::RecoTauPiZeroCollection>(iConfig.getParameter<edm::InputTag>("TauPiZeroCollectionTag"));
    //Taus
    TauHPSCollectionToken_ = consumes<reco::PFTauCollection>(edm::InputTag("hpsPFTauProducer","","RECO"));
    //SIMAOD
    if (!IsData) GenParticlesToken_ = consumes<reco::GenParticleCollection>(edm::InputTag("genParticles","","")); //typedef std::vector<GenParticle> reco::GenParticleCollection

  if (mute) 
  {
    cout.setstate(ios_base::failbit);
    clog.setstate(ios_base::failbit);
  }
  if (debug)
    cout.setstate(ios_base::failbit);
}

void AOD_pi0::RecO_Cand_type(const reco::Candidate* cand)
{
  if (cand->isCaloMuon()) dout("isCaloMuon():", cand->isCaloMuon());
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

// ------------ method called for each event  ------------
void
AOD_pi0::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
    v0_count = 0 ;
    pizero_count = 0;
  //Ks's token
    edm::Handle<reco::VertexCompositeCandidateCollection> Vertices;
    iEvent.getByToken( KshortCollectionToken_, Vertices);

  //HPS Pi0's token
    edm::Handle<reco::RecoTauPiZeroCollection> Strips;
    iEvent.getByToken( TauPiZeroCollectionToken_, Strips);

  //Tau's token based hps
    edm::Handle<reco::PFTauCollection> PF_taus;// typedef vector< PFTau >  PFTauCollection
    iEvent.getByToken( TauHPSCollectionToken_, PF_taus);
  //TauHPSCollectionToken_ = consumes<reco::PFTauCollection>(edm::InputTag("hpsPFTauProducer","","RECO"));

  //SIMAOD's token
    edm::Handle<reco::GenParticleCollection> GenPart;
    if (!IsData) iEvent.getByToken( GenParticlesToken_, GenPart);

  /// RECO Ks's - all are charged
  if (false && Vertices.isValid())
  {
    dlog("RECO Ks's Particles");
    vector<reco::CandidateCollection> v_daughters;
    v0_count = Vertices->size();
    dout("Size:", v0_count);
    for(unsigned i = 0 ; i < Vertices->size() ; i++)
    {
      int num = (*Vertices)[i].numberOfDaughters();//edm::reco::CompositeCandidate::daughters 
      int num_moth = (*Vertices)[i].numberOfMothers();
      dlog("Ks_", i, "(", (*Vertices)[i].charge(), ")");
      dlog("\tdaughters tot number:", num, ";", " moth number:", num_moth); 

      //for (vector<reco::Candidate* >::const_iterator iter = (*Vertices)[i].daughters.begin(); iter != (*Vertices)[i].daughters.end(); ++iter)
      double E = 0, p_x = 0, p_y = 0, p_z = 0, inv_M = 0;
      for( int j = 0; j < num; j++) // Loop over daughters
      {
        const reco::Candidate* daughter = (*Vertices)[i].daughter(j);
        //RecO_Cand_type(daughter); is unknown
        TLorentzVector temp(daughter->px(),daughter->py(),daughter->pz(),daughter->energy());
        dlog("\t\tdaughter_", j, "(", daughter->charge(), ")", " mass:", temp.M());
        if (daughter->pt() > max_ks_daughter_pt) max_ks_daughter_pt = daughter->pt();
        ks_daughter_pt->Fill(daughter->pt());

        E += daughter->energy();
        p_x += daughter->px();
        p_y += daughter->py();
        p_z += daughter->pz();
      }
      inv_M += sqrt(pow(E, 2) - pow(p_x, 2) - pow(p_y, 2) - pow(p_z, 2));
      ks_inv_m_pi->Fill(inv_M);
    }
    // Fill the variables
    h_v0_count->Fill(v0_count);
  }
  else if (!Vertices.isValid()) dlog("UNVALID Ks's");
  
  /// HPS pi0's and taus - with loop among reco::tau
  if (false && Strips.isValid() && false) 
  {
    unsigned int matched_pi = 0;
    if (Strips->size() > 0) dout("Number of HPS Pizeros =", Strips->size());
    if (PF_taus->size() > 0) dout("Number of RECO Tau =", PF_taus->size());

    for (unsigned int i = 0; i < Strips->size(); i++) 
    {
      if (!debug) dout("\tPi0_", i, "(", (*Strips)[i].charge(), ")", (*Strips)[i].px(), (*Strips)[i].py(), (*Strips)[i].pz(), (*Strips)[i].p());
      pizero_px[pizero_count] = (*Strips)[i].px();
      pizero_py[pizero_count] = (*Strips)[i].py();
      pizero_pz[pizero_count] = (*Strips)[i].pz();
      pizero_e[pizero_count]  = (*Strips)[i].p();
      pizero_pt[pizero_count] = (*Strips)[i].pt();
      pizero_eta[pizero_count] = (*Strips)[i].eta();
      pizero_phi[pizero_count] = (*Strips)[i].phi();
      pizero_x[pizero_count] = (*Strips)[i].vx();
      pizero_y[pizero_count] = (*Strips)[i].vy();
      pizero_z[pizero_count] = (*Strips)[i].vz();
      pizero_count++;
      if (!debug) dout(" PI0", i, "px =", pizero_px[pizero_count], "py =", pizero_py[pizero_count], "pz =", pizero_pz[pizero_count], 
                                 "; vx =", pizero_x[pizero_count],  "vy =", pizero_y[pizero_count],  "vz =", pizero_z[pizero_count]);

      if (pizero_count>=1000) 
      {
        cerr << "number of pizeros > 1000. They are missing." << endl; 
        break;
      }

      // RECO Taus
        bool foundpi = false;
        for (unsigned int i_tau = 0; i_tau < PF_taus->size(); i_tau++) //over all taus
        {
          cout << "\t\tTau_" << i_tau << endl;
          reco::PFTauRef pftauref(PF_taus, i_tau);

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

  /// Only hps reco Taus
  if (false && PF_taus.isValid() )
  {
    if (PF_taus->size() > 0)  dout("Number of Tau = ", PF_taus->size());
    for (unsigned int i = 0; i < PF_taus->size(); i++) // Over Tau's
    {
      reco::PFTauRef pftauref(PF_taus, i); // one Tau instance, typedef edm::Ref<PFTauCollection> reco::PFTauRef
      dlog("tau #", i, "from", PF_taus->size(), "(", pftauref->vx(), pftauref->vy(), pftauref->vz(), ")");
        //reco::PFTau pfTau(PF_taus, i); //same as (*PF_taus)[i]
        //edm::AtomicPtrCache< vector< reco::RecoTauPiZero > >
        //const vector< reco::RecoTauPiZero > a  = (*PF_taus)[i].signalPiZeroCandidates();
        //(*reco::PFTau)pfTau->signalPiZeroCandidates();
        //const vector < RecoTauPiZero > &   isolationPiZeroCandidates () const
        //const vector< RecoTauPiZero > & reco::PFTau::isolationPiZeroCandidates (   ) const
        //const vector < RecoTauPiZero > &   signalPiZeroCandidates () const

      // Lists to be used 
        //Pizeros list
          // Signal pi0's
            vector < reco::RecoTauPiZero > tau_pizeros_sig = pftauref->signalPiZeroCandidates();
          //Isolation pi0's
            const vector < reco::RecoTauPiZero > tau_pizeros_isol = pftauref->isolationPiZeroCandidates(); // pions of the considered Tau
            vector < reco::RecoTauPiZero *> point_tau_pizeros_isol;
          // All pi0's
            vector < reco::RecoTauPiZero > tau_pizeros = pftauref->signalPiZeroCandidates();
                                           tau_pizeros.insert(tau_pizeros.end(), tau_pizeros_isol.begin(), tau_pizeros_isol.end());
            vector < reco::RecoTauPiZero* > point_tau_pizeros;
       
        // Hadrons list
          vector < reco::PFCandidatePtr  > tau_signalPFChargedHadrCands    = pftauref->signalPFChargedHadrCands();//typedef edm::Ptr<PFCandidate> reco::PFCandidatePtr
          vector < reco::PFCandidatePtr  > tau_isolationPFChargedHadrCands = pftauref->isolationPFChargedHadrCands(); // //typedef edm::Ptr<PFCandidate> reco::PFCandidatePtr
            // All pi+-'s
            vector < reco::PFCandidatePtr > tau_picharge = pftauref->signalPFChargedHadrCands(); // vector < edm::Ptr<PFCandidate> > 
                                            tau_picharge.insert(tau_picharge.end(), tau_isolationPFChargedHadrCands.begin(), tau_isolationPFChargedHadrCands.end());
            std::vector<reco::PFCandidatePtr * > point_tau_picharge;

          vector < reco::PFCandidatePtr  > tau_signalPFNeutrHadrCands    = pftauref->signalPFNeutrHadrCands(); // //typedef edm::Ptr<PFCandidate> reco::PFCandidatePtr
          vector < reco::PFCandidatePtr  > tau_isolationPFNeutrHadrCands = pftauref->isolationPFNeutrHadrCands(); // //typedef edm::Ptr<PFCandidate> reco::PFCandidatePtr
            // All pi0_had's
            vector < reco::PFCandidatePtr > tau_pizeros_had = pftauref->signalPFNeutrHadrCands(); // //typedef edm::Ptr<PFCandidate> reco::PFCandidatePtr
                                            tau_pizeros_had.insert(tau_pizeros_had.end(), tau_isolationPFNeutrHadrCands.begin(), tau_isolationPFNeutrHadrCands.end());
            std::vector<reco::PFCandidatePtr *> point_tau_pizeros_had;
        
        //NOT USED list
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
          num_pions->Fill(tau_pizeros_sig.size());
          double inv_M = 0;
          if (tau_pizeros_sig.size() > 0)
          {
            TLorentzVector first(0, 0, 0, 0);
            for (unsigned int j = 0; j < tau_pizeros_sig.size(); j++)
            {
              dout("\t\t\tPi0_", j, "(", tau_pizeros_sig[j].charge(), ")", tau_pizeros_sig[j].px(), tau_pizeros_sig[j].py(), tau_pizeros_sig[j].pz(), tau_pizeros_sig[j].energy());
              TLorentzVector second(tau_pizeros_sig[j].px(), tau_pizeros_sig[j].py(), tau_pizeros_sig[j].pz(), tau_pizeros_sig[j].energy());
              first = (first + second);
              
              dlog("\t\t\t signal tau pi0_", j, ":", tau_pizeros_sig[j].vx(), tau_pizeros_sig[j].vy(), tau_pizeros_sig[j].vz(), 
                                                ":", tau_pizeros_sig[j].px(), tau_pizeros_sig[j].py(), tau_pizeros_sig[j].pz());
            }
            inv_M = first.M();
          }
          dout("TAU_", i, " WITH TOTAL INVARIANT MASS OF PIONS:", inv_M);
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
  else if (!PF_taus.isValid()) dout("no valid PF_taus");

  /// GEN Particles
  dlog("GEN Particles");
  if (!IsData && GenPart.isValid())
  {
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
    v_daughters_k0l_to_pi0.clear();
    BuildTo(GenPart, h_gen_k0_all_to_pi0, v_daughters_k0l_to_pi0, 130, map_kaons, 2, 111, map_pions);// k0l 130
    //v_daughters_k0_to_pi0.clear();
    //BuildTo(GenPart, h_gen_k0_all_to_pi0, v_daughters_k0_to_pi0, 311, map_kaons, 2, 111, map_pions);// k0 311 - is empty by definition 

    //reco::CandMatchMap map = MCTruthDeltaRMatcher("pions and pizeros")
  }
  else if (!IsData) dlog("\tno GenPart.isValid()");
  else dlog("no GenPart in data");
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
  dlog("Found K", mom_pdgid, ":", daughters.size());
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

void AOD_pi0::GenEvolution(const reco::Candidate * mom, int num_of_tabs = 2)
{
    int n = mom->numberOfDaughters();
    for(int i = 0; i < n; i++) 
    {
      const reco::Candidate * daughter = mom->daughter(i);
      //if (abs(daughter->pdgId()) == 311/*(abs(daughter->pdgId()) % 1000 == 311 || abs(daughter->pdgId()) % 1000 == 321)*/ && daughter->numberOfDaughters() == 0 ) {dlog(string(num_of_tabs, '\t'), i, ") is Kaon", daughter->pdgId(), daughter->status(), daughter->numberOfDaughters());exit(1);}
      ofs << string(num_of_tabs, '\t') << i << ") " << daughter->pdgId() << " " << daughter->status() << " " << daughter->numberOfDaughters() << endl;//dlog(string(num_of_tabs, '\t'), i, ") ", daughter->pdgId(), daughter->status(), daughter->numberOfDaughters()); 
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
          if (typeOfConstituences.Contains("pi+-") && collection[i]->charge() * collection[j]->charge() < 0) continue;
          TLorentzVector second(collection[j]->px(), collection[j]->py(), collection[j]->pz(), collection[j]->energy());
          double inv_M = (first + second).M();
          hist_inv_m->Fill(inv_M);
          dout("\t\t\tm(", typeOfConstituences, i, "+", typeOfConstituences, j, ") =", inv_M);
          if (hist_pt != 0) hist_pt->Fill((first + second).Pt());
        }
        //cout,  endl;
        //if (hist_pt != 0) hist_pt->Fill(collection[i]->pt());
    }
    //if (hist_pt != 0) hist_pt->Fill(collection[collection.size() - 1]->pt());
    dlog("\t\t\t", typeOfCollection, collection.size() - 1, ":", collection[collection.size() - 1]->vx(), collection[collection.size() - 1]->vy(), collection[collection.size() - 1]->vz(), 
                                                            ":", collection[collection.size() - 1]->px(), collection[collection.size() - 1]->py(), collection[collection.size() - 1]->pz());
  }
  else if (collection.size() == 1) dlog("\t\t\t", typeOfCollection, 0, ":", collection[0]->vx(), collection[0]->vy(), collection[0]->vz(), 
                                           ":", collection[0]->px(), collection[0]->py(), collection[0]->pz());
}
template  <typename VectoreType >
void AOD_pi0::CombinatoricOfTwoToNeutralInvM(vector <VectoreType> collection, 
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
          hist_inv_m->Fill(inv_M);
          dout("\t\t\tm(", typeOfConstituences, i, "+", typeOfConstituences, j, ") =", inv_M);
          if (hist_pt != 0) hist_pt->Fill((first + second).Pt());
        }
        // if (hist_pt != 0) hist_pt->Fill(collection[i]->pt());
        //cout,  endl;
    }
    //if (hist_pt != 0) hist_pt->Fill(collection[collection.size() - 1]->pt());
    dlog("\t\t\t", typeOfCollection, collection.size() - 1, ":", collection[collection.size() - 1]->vx(), collection[collection.size() - 1]->vy(), collection[collection.size() - 1]->vz(), 
                                                            ":", collection[collection.size() - 1]->px(), collection[collection.size() - 1]->py(), collection[collection.size() - 1]->pz());
  }
  else if (collection.size() == 1) dlog("\t\t\t", typeOfCollection, 0, ":", collection[0]->vx(), collection[0]->vy(), collection[0]->vz(), 
                                        ":", collection[0]->px(), collection[0]->py(), collection[0]->pz());
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
  h_v0_count->Write();
  pions_inv_m->Write();
  num_pions->Write();
  taus_isol_pi0_inv_m_to_ks->Write();
  ks_daughter_pt->Write();
  ks_inv_m_pi->Write();
  taus_isol_pi0_inv_pt->Write();
  taus_pi0_inv_m_to_ks->Write();
  taus_pi0_inv_pt->Write();
  taus_pi_charged_inv_m_to_ks->Write();
  taus_pi_charged_inv_pt->Write();
  h_gen_k0_all_to_pi0->Write();
}

AOD_pi0::~AOD_pi0()
{
  dlog("Total num of not-matched events:", num_ev_tau_pi_not_in_hps_pi);
  dlog("num_pion_res: ", num_pion_res);
  dlog("max_inv_mass: ", max_inv_mass);
  dlog("max_ks_daughter_pt: ", max_ks_daughter_pt);
  outfile->Close();
  ofs.close();
  //delete v_daughters_k0s_to_pi0;
  //delete v_daughters_k0l_to_pi0;
  //delete v_daughters_k0_to_pi0;
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

//define this as a plug-in
DEFINE_FWK_MODULE(AOD_pi0);
