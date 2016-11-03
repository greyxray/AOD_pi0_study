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

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"

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
    void CombinatoricOfTwoInvM(vector<VectoreType*> collection/*reco::RecoTauPiZero */, 
                                    TString typeOfCollection, 
                                    TString typeOfObjects,
                                    TString typeOfConstituences, 
                                    TH1 * hist_inv_m, 
                                    TH1 * hist_pt=0);
    template  <typename VectoreType >
    void CombinatoricOfTwoInvM(vector<VectoreType> collection/*reco::RecoTauPiZero */, 
                                    TString typeOfCollection, 
                                    TString typeOfObjects,
                                    TString typeOfConstituences, 
                                    TH1 * hist_inv_m, 
                                    TH1 * hist_pt=0);
    template <typename T>
    vector <T*> TransformToPointers(vector <T> a, vector <T*> b);

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

    // Histograms
      TH1D* h_v0_count;
      TH1D* pions_inv_m;
      TH1D* num_pions;
      TH1D* taus_isol_pi0_inv_m_to_ks;
      TH1D* taus_isol_pi0_pt;
      TH1D* taus_pi0_inv_m_to_ks;
      TH1D* taus_pi0_pt;
      TH1D* ks_daughter_pt;
      TH1D* ks_inv_m_pi;
      TH1D* taus_pi_charged_inv_m_to_ks ;
      TH1D* taus_pi_charged_pt;
      TH1D* taus_pi0_had_inv_m_to_ks ;
      TH1D* taus_pi0_had_pt;

    // Tokens for the Collections 
      edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> KshortCollectionToken_;
      edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> LambdaCollectionToken_;
      edm::EDGetTokenT<reco::PFCandidateCollection> PFCandidateCollectionToken_;

      edm::EDGetTokenT<reco::RecoTauPiZeroCollection> TauPiZeroCollectionToken_;//hpsPFTauProducer

      edm::EDGetTokenT<reco::PFTauCollection> TauHPSCollectionToken_;

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
    bool crecpizero;
    bool debug;
    bool mute;
    int num_pion_res;
    double max_inv_mass;
    double max_ks_daughter_pt;
    unsigned int num_ev_tau_pi_not_in_hps_pi;
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
  crecpizero(iConfig.getUntrackedParameter<bool>("RecPiZero", false)),
  debug(true),
  mute(false),
  num_pion_res(0),
  max_inv_mass(0),
  max_ks_daughter_pt(0),
  num_ev_tau_pi_not_in_hps_pi(0)
{
  usesResource("TFileService");

  // Saved histograms 
    outfile = new TFile("aod_pi0.root","RECREATE");
    outfile->cd();  
      h_v0_count = new TH1D("v0_count","v0 count", 10, 0, 9);
      pions_inv_m  = new TH1D("pions_inv_m","inv mass of pions in taus", 100, 0, 1);
      num_pions = new TH1D("num_of_pios","num of pions", 10, 0, 9);
    hist_directory[1]  = outfile->mkdir("ks_coll", "ks_collection");
    hist_directory[1]->cd();  //= outfile->mkdir("ks_coll", "ks_collection");
      ks_daughter_pt = new TH1D("ks_daughter_pt","ks daughters pt", 1000, 0, 10);
      ks_inv_m_pi = new TH1D("ks_inv_m_pi","ks daughters inv mass", 1000, 0, 10);
    hist_directory[0]  = outfile->mkdir("Taus_pions_coll", "Taus_pions_collections");
    hist_directory[0]->cd();  //= outfile->mkdir("Taus_pions_coll", "Taus_pions_collections");
      taus_isol_pi0_inv_m_to_ks = new TH1D("taus_isol_pi0_inv_m_to_ks","all pairs of tau isolation pions inv mass", 1000, 0, 17);
      taus_isol_pi0_pt = new TH1D("taus_isol_pi0_pt","all pairs of tau isolation pionspt", 1000, 0, 10);
      taus_pi0_inv_m_to_ks = new TH1D("taus_pi0_inv_m_to_ks","all pairs of tau pions inv mass", 1000, 0, 17);
      taus_pi0_pt = new TH1D("taus_pi0_pt","all pairs of tau pions pt", 1000, 0, 10);
    hist_directory[2]  = outfile->mkdir("Taus_charged_had_coll", "Taus_charged_had_coll");
    hist_directory[2]->cd();  //= outfile->mkdir("Taus_charged_had_coll", "Taus_charged_had_coll");
      taus_pi_charged_inv_m_to_ks = new TH1D("taus_pi_charged_inv_m_to_ks","all pairs of tau charged pions inv mass", 1000, 0, 17);
      taus_pi_charged_pt = new TH1D("taus_pi_charged_pt","all pairs of tau charged pions pt", 1000, 0, 10);
    hist_directory[3]  = outfile->mkdir("Taus_neutral_had_coll", "Taus_neutral_had_coll");
    hist_directory[3]->cd();  //= outfile->mkdir("Taus_neutral_had_coll", "Taus_neutral_had_coll");
      taus_pi0_had_inv_m_to_ks = new TH1D("taus_pi0_had_inv_m_to_ks","all pairs of tau pions 0 had inv mass", 1000, 0, 17);
      taus_pi0_had_pt = new TH1D("taus_pi0_had_pt","all pairs of tau pions 0 had pt", 1000, 0, 10);

  // Tokens
    //Ks's
    KshortCollectionToken_ = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("KshortCollectionTag"));//vector<reco::VertexCompositeCandidate>    "generalV0Candidates"       "Kshort"          "RECO"  
    LambdaCollectionToken_ = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("LambdaCollectionTag"));
    //Pi0
    TauPiZeroCollectionToken_ = consumes<reco::RecoTauPiZeroCollection>(iConfig.getParameter<edm::InputTag>("TauPiZeroCollectionTag"));
    //Taus
    TauHPSCollectionToken_ = consumes<reco::PFTauCollection>(edm::InputTag("hpsPFTauProducer","","RECO"));

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

  //Tau's token
    edm::Handle<reco::PFTauCollection> pf_taus;// typedef vector< PFTau >  PFTauCollection
    iEvent.getByToken( TauHPSCollectionToken_, pf_taus);
  //TauHPSCollectionToken_ = consumes<reco::PFTauCollection>(edm::InputTag("hpsPFTauProducer","","RECO"));


  /// RECO Ks's
  if (Vertices.isValid())
  {
    vector<reco::CandidateCollection> v_daughters;
    v0_count = Vertices->size();
    dout("Size:", v0_count);
    for(unsigned i = 0 ; i < Vertices->size() ; i++)
    {
      int num = (*Vertices)[i].numberOfDaughters();//edm::reco::CompositeCandidate::daughters 
      int num_moth = (*Vertices)[i].numberOfMothers();
      dout("Ks_", i, "(", (*Vertices)[i].charge(), ")");
      dlog("\tdaughters number:", num, ";", " moth number:", num_moth); 

      //for (vector<reco::Candidate* >::const_iterator iter = (*Vertices)[i].daughters.begin(); iter != (*Vertices)[i].daughters.end(); ++iter)
      double E = 0, p_x = 0, p_y = 0, p_z = 0, inv_M = 0;
      for( int j = 0; j < num; j++) // Loop over daughters
      {
        const reco::Candidate* daughter = (*Vertices)[i].daughter(j);
        RecO_Cand_type(daughter);

        dlog("daughter_", j, "(", daughter->charge(), ")", " pt:", daughter->pt());
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
  else if (!Vertices.isValid()) dout("UNVALID VERTICES");

  /// HPS pi0's - with loop among reco::tau
  if (Strips.isValid() && false) 
  {
    unsigned int matched_pi = 0;
    if (Strips->size() > 0) dout("Number of HPS Pizeros =", Strips->size());
    if (pf_taus->size() > 0) dout("Number of RECO Tau =", pf_taus->size());

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
        for (unsigned int i_tau = 0; i_tau < pf_taus->size(); i_tau++) //over all taus
        {
          cout << "\t\tTau_" << i_tau << endl;
          reco::PFTauRef pftauref(pf_taus, i_tau);

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

  /// RECO Taus
  if (pf_taus.isValid() )
  {
    if (pf_taus->size() > 0)  dout("Number of Tau = ", pf_taus->size());
    for (unsigned int i = 0; i < pf_taus->size(); i++) // Over Tau's
    {
      reco::PFTauRef pftauref(pf_taus, i); // one Tau instance, typedef edm::Ref<PFTauCollection> reco::PFTauRef
      dlog("tau #", i, "from", pf_taus->size(), "(", pftauref->vx(), pftauref->vy(), pftauref->vz(), ")");
        //reco::PFTau pfTau(pf_taus, i); //same as (*pf_taus)[i]
        //edm::AtomicPtrCache< vector< reco::RecoTauPiZero > >
        //const vector< reco::RecoTauPiZero > a  = (*pf_taus)[i].signalPiZeroCandidates();
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
      dlog("\t-----------------------------------------");
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
        CombinatoricOfTwoInvM(point_tau_pizeros_isol, "isolat tau pi0", "tau", "pi0", taus_isol_pi0_inv_m_to_ks, taus_isol_pi0_pt);
        
        //All pi0's loop
        point_tau_pizeros = TransformToPointers(tau_pizeros, point_tau_pizeros);
        CombinatoricOfTwoInvM(point_tau_pizeros, "all tau pi0", "tau", "pi0", taus_pi0_inv_m_to_ks, taus_pi0_pt);

      //Charged hadrons of PFRecoTau - all PFCandidatePtr
      dlog("\t-----------------------------------------");
      dlog("\tCharged hadrons of PFRecoTau ");
      point_tau_picharge = TransformToPointers(tau_picharge, point_tau_picharge);
      CombinatoricOfTwoInvM(tau_picharge, "all tau pi+-", "tau", "pi+-", taus_pi_charged_inv_m_to_ks, taus_pi_charged_pt);
      // if (tau_picharge.size() > 1)
      // {
      //   dlog("\t\t all tau pi+-'s number:", tau_picharge.size());

      //   for (unsigned int j_pi = 0; j_pi < tau_picharge.size() - 1; j_pi++) //over isolation pions in Tau
      //   {
      //     dlog("\t\t\t all tau pi+_", j_pi, ":", tau_picharge[j_pi]->vx(), tau_picharge[j_pi]->vy(), tau_picharge[j_pi]->vz(),  
      //                                       ":", tau_picharge[j_pi]->px(), tau_picharge[j_pi]->py(), tau_picharge[j_pi]->pz());

      //     //Check the vertex position
      //      if (false ) dout("\t\t\tPi+_", j_pi, "(", tau_picharge[j_pi]->charge(), "), vertex:", tau_picharge[j_pi]->vx(), tau_picharge[j_pi]->vy(), tau_picharge[j_pi]->vz());//all from the same vertex
      //       TLorentzVector first(tau_picharge[j_pi]->px(), tau_picharge[j_pi]->py(), tau_picharge[j_pi]->pz(), tau_picharge[j_pi]->energy());
      //     //Build the inv mass
      //       dout("\t\t===>matching to Pi+_", j_pi, "from", tau_picharge.size(), "in this tau");
      //       for (unsigned int j2_pi = j_pi + 1; j2_pi < tau_picharge.size(); j2_pi++)
      //       {
      //         TLorentzVector second(tau_picharge[j2_pi]->px(), tau_picharge[j2_pi]->py(), tau_picharge[j2_pi]->pz(), tau_picharge[j2_pi]->energy());
      //         double inv_M = (first + second).M();//(*(tau_picharge[j_pi].get()) + *(tau_picharge[j2_pi].get()))->mass();// 
              
      //         taus_pi_charged_inv_m_to_ks->Fill(inv_M);
      //         dout("\t\t\tm( Pi+_", j_pi, "+ Pi+_", j2_pi, ") =", inv_M);
      //       }
      //       dout();

      //       taus_pi_charged_pt->Fill(tau_picharge[j_pi]->pt());
      //   }
      //   taus_pi_charged_pt->Fill(tau_picharge[tau_picharge.size() - 1]->pt());
      //   dlog("\t\t\t all tau pi+_", tau_picharge.size() - 1, ":", tau_picharge[tau_picharge.size() - 1]->vx(), tau_picharge[tau_picharge.size() - 1]->vy(), tau_picharge[tau_picharge.size() - 1]->vz(), 
      //                                                        ":", tau_picharge[tau_picharge.size() - 1]->px(), tau_picharge[tau_picharge.size() - 1]->py(), tau_picharge[tau_picharge.size() - 1]->pz());
      // }

      //Neutral hadrons of PFRecoTau - all PFCandidatePtr
      dlog("\t-----------------------------------------");
      dlog("\tNeutral hadrons of PFRecoTau ");
      point_tau_pizeros_had = TransformToPointers(tau_pizeros_had, point_tau_pizeros_had);
      CombinatoricOfTwoInvM(tau_pizeros_had, "all tau pi0_had", "tau", "pi0_had", taus_pi0_had_inv_m_to_ks, taus_pi0_had_pt);
      // if (tau_pizeros_had.size() > 1)
      // {
      //   dlog("\t\t all tau pi0_had's number:", tau_pizeros_had.size());

      //   for (unsigned int j_pi = 0; j_pi < tau_pizeros_had.size() - 1; j_pi++) //over isolation pions in Tau
      //   {
      //     dlog("\t\t\t all tau pi0_had", j_pi, ":", tau_pizeros_had[j_pi]->vx(), tau_pizeros_had[j_pi]->vy(), tau_pizeros_had[j_pi]->vz(),  
      //                                       ":", tau_pizeros_had[j_pi]->px(), tau_pizeros_had[j_pi]->py(), tau_pizeros_had[j_pi]->pz());

      //     //Check the vertex position
      //      if (false ) dout("\t\t\tPi0_had", j_pi, "(", tau_pizeros_had[j_pi]->charge(), "), vertex:", tau_pizeros_had[j_pi]->vx(), tau_pizeros_had[j_pi]->vy(), tau_pizeros_had[j_pi]->vz());//all from the same vertex
      //       TLorentzVector first(tau_pizeros_had[j_pi]->px(), tau_pizeros_had[j_pi]->py(), tau_pizeros_had[j_pi]->pz(), tau_pizeros_had[j_pi]->energy());
      //     //Build the inv mass
      //       dout("\t\t===>matching to Pi0_had", j_pi, "from", tau_pizeros_had.size(), "in this tau");
      //       for (unsigned int j2_pi = j_pi + 1; j2_pi < tau_pizeros_had.size(); j2_pi++)
      //       {
      //         TLorentzVector second(tau_pizeros_had[j2_pi]->px(), tau_pizeros_had[j2_pi]->py(), tau_pizeros_had[j2_pi]->pz(), tau_pizeros_had[j2_pi]->energy());
      //         double inv_M = (first + second).M();//(*(tau_pizeros_had[j_pi].get()) + *(tau_pizeros_had[j2_pi].get()))->mass();// 
              
      //         taus_pi0_had_inv_m_to_ks->Fill(inv_M);
      //         dout("\t\t\tm( Pi0_had_", j_pi, "+ Pi0_had_", j2_pi, ") =", inv_M);
      //       }
      //       dout();

      //       taus_pi0_had_pt->Fill(tau_pizeros_had[j_pi]->pt());
      //   }
      //   taus_pi0_had_pt->Fill(tau_pizeros_had[tau_pizeros_had.size() - 1]->pt());
      //   dlog("\t\t\t all tau pi0_had_", tau_pizeros_had.size() - 1, ":", tau_pizeros_had[tau_pizeros_had.size() - 1]->vx(), tau_pizeros_had[tau_pizeros_had.size() - 1]->vy(), tau_pizeros_had[tau_pizeros_had.size() - 1]->vz(), 
      //                                                               ":", tau_pizeros_had[tau_pizeros_had.size() - 1]->px(), tau_pizeros_had[tau_pizeros_had.size() - 1]->py(), tau_pizeros_had[tau_pizeros_had.size() - 1]->pz());
      // }
    }
  } 
  else if (!pf_taus.isValid()) dout("no valid pf_taus");
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
void AOD_pi0::CombinatoricOfTwoInvM(vector <VectoreType *> collection, 
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
          TLorentzVector second(collection[j]->px(), collection[j]->py(), collection[j]->pz(), collection[j]->energy());
          double inv_M = (first + second).M();
          hist_inv_m->Fill(inv_M);
          dout("\t\t\tm(", typeOfConstituences, i, "+", typeOfConstituences, j, ") =", inv_M);
        }
        //cout,  endl;
        if (hist_pt != 0) hist_pt->Fill(collection[i]->pt());
    }
    if (hist_pt != 0) hist_pt->Fill(collection[collection.size() - 1]->pt());
    dlog("\t\t\t", typeOfCollection, collection.size() - 1, ":", collection[collection.size() - 1]->vx(), collection[collection.size() - 1]->vy(), collection[collection.size() - 1]->vz(), 
                                                            ":", collection[collection.size() - 1]->px(), collection[collection.size() - 1]->py(), collection[collection.size() - 1]->pz());
  }
  else if (collection.size() == 1) dlog("\t\t\t", typeOfCollection, 0, ":", collection[0]->vx(), collection[0]->vy(), collection[0]->vz(), 
                                           ":", collection[0]->px(), collection[0]->py(), collection[0]->pz());
}
      

template  <typename VectoreType >
void AOD_pi0::CombinatoricOfTwoInvM(vector <VectoreType> collection, 
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
          TLorentzVector second(collection[j]->px(), collection[j]->py(), collection[j]->pz(), collection[j]->energy());
          double inv_M = (first + second).M();
          hist_inv_m->Fill(inv_M);
          dout("\t\t\tm(", typeOfConstituences, i, "+", typeOfConstituences, j, ") =", inv_M);
        }
        //cout,  endl;
        if (hist_pt != 0) hist_pt->Fill(collection[i]->pt());
    }
    if (hist_pt != 0) hist_pt->Fill(collection[collection.size() - 1]->pt());
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
  taus_isol_pi0_pt->Write();
  taus_pi0_inv_m_to_ks->Write();
  taus_pi0_pt->Write();
  taus_pi_charged_inv_m_to_ks->Write();
  taus_pi_charged_pt->Write();
}

AOD_pi0::~AOD_pi0()
{
  dlog("Total num of not-matched events:", num_ev_tau_pi_not_in_hps_pi);
  dlog("num_pion_res: ", num_pion_res);
  dlog("max_inv_mass: ", max_inv_mass);
  dlog("max_ks_daughter_pt: ", max_ks_daughter_pt);
  outfile->Close();
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
/*FYI
  isolationPFChargedHadrCands() const reco::PFTau 
  isolationTauChargedHadronCandidates() const reco::PFTau 
  //isolationPFChargedHadrCandsPtSum() const  reco::PFTau 
    isolationPFNeutrHadrCands() const reco::PFTau 


  signalPFChargedHadrCands() const  reco::PFTau 
  signalTauChargedHadronCandidates() const  reco::PFTau 
    signalPFNeutrHadrCands() const  reco::PFTau 


  leadPFChargedHadrCand() const reco::PFTau 
  leadPFChargedHadrCandsignedSipt() const reco::PFTau 
  leadTauChargedHadronCandidate() const reco::PFTau 
========
   PFCandidatePtr   isolationPFChargedHadrCands() const reco::PFTau 
  PFRecoTauChargedHadron  isolationTauChargedHadronCandidates() const reco::PFTau 
  //isolationPFChargedHadrCandsPtSum() const  reco::PFTau 
  PFCandidatePtr  isolationPFNeutrHadrCands() const reco::PFTau 

  const PFCandidatePtr & leadPFChargedHadrCand() const reco::PFTau 
  
float  leadPFChargedHadrCandsignedSipt() const reco::PFTau 
 PFRecoTauChargedHadronRef =>typedef edm::Ref<PFRecoTauChargedHadronCollection> =>typedef vector<PFRecoTauChargedHadron> with #include <PFRecoTauChargedHadron.h>  leadTauChargedHadronCandidate() const reco::PFTau 

*/
