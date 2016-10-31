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

//
// class declaration
//

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


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    // ----------member data ---------------------------
    TFile* outfile;

    // Histograms
      TH1D* h_v0_count;
      TH1D* pions_inv_m;
      TH1D* num_pions;
      TH1D* isol_pi_inv_m_to_ks;
      TH1D* ks_daughter_pt;
      TH1D* ks_inv_m_pi;
      TH1D* isol_pi_pt;

    // Tokens for the Collections 
      edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> KshortCollectionToken_;
      edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> LambdaCollectionToken_;
      edm::EDGetTokenT<reco::PFCandidateCollection> PFCandidateCollectionToken_;

      edm::EDGetTokenT<reco::RecoTauPiZeroCollection> TauPiZeroCollectionToken_;//hpsPFTauProducer

      edm::EDGetTokenT<reco::PFTauCollection> TauHPSCollectionToken_;


    
    Int_t v0_count; // V0s - Ks or Lambdas
    bool crecpizero;
    bool debug;
    int num_pion_res;
    double max_inv_mass;
    double max_ks_daughter_pt;

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

    //Tau
    Float_t tau_px[1000];
      Float_t tau_py[1000];
      Float_t tau_pz[1000];
      Float_t tau_e[1000];
      Float_t tau_pt[1000];
      Float_t tau_eta[1000];
      Float_t tau_phi[1000];
      Float_t tau_x[1000];
      Float_t tau_y[1000];
      Float_t tau_z[1000];
      UInt_t tau_count;

      unsigned int num_ev_tau_pi_not_in_hps_pi;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
AOD_pi0::AOD_pi0(const edm::ParameterSet& iConfig):
  crecpizero(iConfig.getUntrackedParameter<bool>("RecPiZero", false))
{
  debug = true;
  num_ev_tau_pi_not_in_hps_pi = 0;
  num_pion_res = 0;
  max_inv_mass = 0;
  max_ks_daughter_pt = 0;
  //now do what ever initialization is needed
  usesResource("TFileService");

  outfile = new TFile("aod_pi0.root","RECREATE");
  // Saved brunches 
  h_v0_count = new TH1D("v0_count","v0 count", 10, 0, 9);
  pions_inv_m  = new TH1D("pions_inv_m","inv mass of pions in taus", 100, 0, 1);
  num_pions = new TH1D("num_of_pios","num of pions", 10, 0, 9);
  isol_pi_inv_m_to_ks = new TH1D("isol_pi_inv_m_to_ks","all pairs of pions inv mass", 1000, 0, 17);
  ks_daughter_pt = new TH1D("ks_daughter_pt","ks daughters pt", 1000, 0, 10);
  ks_inv_m_pi = new TH1D("ks_inv_m_pi","ks daughters inv mass", 1000, 0, 10);
  isol_pi_pt = new TH1D("isol_pi_pt","isol_pi_pt", 1000, 0, 10);

  // Tokens
    //Ks's
    //vector<reco::VertexCompositeCandidate>    "generalV0Candidates"       "Kshort"          "RECO"  
    KshortCollectionToken_ = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("KshortCollectionTag"));
    LambdaCollectionToken_ = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("LambdaCollectionTag"));
    
    //Pi0
    TauPiZeroCollectionToken_ = consumes<reco::RecoTauPiZeroCollection>(iConfig.getParameter<edm::InputTag>("TauPiZeroCollectionTag"));
    
    //Taus
    TauHPSCollectionToken_ = consumes<reco::PFTauCollection>(edm::InputTag("hpsPFTauProducer","","RECO"));

}




//
// member functions
//

// ------------ method called for each event  ------------
void
AOD_pi0::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  bool ev_tau_pi_not_in_hps_pi = false;
  unsigned int matched_pi = 0;
  using namespace edm;
    v0_count = 0 ;
    pizero_count = 0;
    tau_count = 0;
  //Ks's token
    edm::Handle<reco::VertexCompositeCandidateCollection> Vertices;
    iEvent.getByToken( KshortCollectionToken_, Vertices);

  //HPS Pi0's token
    edm::Handle<reco::RecoTauPiZeroCollection> Strips;
    iEvent.getByToken( TauPiZeroCollectionToken_, Strips);

  //Tau's token
    edm::Handle<reco::PFTauCollection> pf_taus;// typedef std::vector< PFTau >  PFTauCollection
    iEvent.getByToken( TauHPSCollectionToken_, pf_taus);
  //TauHPSCollectionToken_ = consumes<reco::PFTauCollection>(edm::InputTag("hpsPFTauProducer","","RECO"));


  /// RECO Ks's
  if (Vertices.isValid())
  {
    std::vector<reco::CandidateCollection> v_daughters;
    v0_count = Vertices->size();
    if (!debug) cout << "Size: " << v0_count << endl;
    for(unsigned i = 0 ; i < Vertices->size() ; i++)
    {
      //edm::reco::CompositeCandidate::daughters 
      int num = (*Vertices)[i].numberOfDaughters();
      int num_moth = (*Vertices)[i].numberOfMothers();
      cout << "Ks_"<< i << "(" << (*Vertices)[i].charge() << ")" << endl;
      if (debug) cout << "\tdaughters number: " << num << ";" << " moth number: " << num_moth << endl;
      //reco::Candidate* a = (*Vertices)[i].daughter(0);//reco::CompositeCandidate::daughters

      //for (std::vector<reco::Candidate* >::const_iterator iter = (*Vertices)[i].daughters.begin(); iter != (*Vertices)[i].daughters.end(); ++iter)
      double E = 0;
      double p_x = 0;
      double p_y = 0;
      double p_z = 0;
      double inv_M = 0;
      for( int j = 0; j < num; j++) // Loop over daughters
      {
        const reco::Candidate* daughter = (*Vertices)[i].daughter(j);
        if (daughter->isCaloMuon()) cout << "isCaloMuon(): " << daughter->isCaloMuon() << endl;
          else if (daughter->isConvertedPhoton()) cout << "isConvertedPhoton(): " << daughter->isConvertedPhoton() << endl;
          else if (daughter->isElectron()) cout << "isElectron(): " << daughter->isElectron() << endl;
          else if (daughter->isGlobalMuon()) cout << "isGlobalMuon(): " << daughter->isGlobalMuon() << endl;
          else if (daughter->isJet()) cout << "isJet(): " << daughter->isJet() << endl;
          else if (daughter->isMuon()) cout << "isMuon(): " << daughter->isMuon() << endl;
          else if (daughter->isPhoton()) cout << "isPhoton(): " << daughter->isPhoton() << endl;
          else if (daughter->isStandAloneMuon()) cout << "isStandAloneMuon(): " << daughter->isStandAloneMuon() << endl;
          else if (daughter->isTrackerMuon()) cout << "isTrackerMuon(): " << daughter->isTrackerMuon() << endl;
        //else cout << "daughter of unknown type" << endl;

        cout << "daughter_" << j << "(" << daughter->charge() << ")" << " pt: " << daughter->pt() << endl;
        if (daughter->pt() > max_ks_daughter_pt) max_ks_daughter_pt = daughter->pt();
        ks_daughter_pt->Fill(daughter->pt());

        E += daughter->energy();
        p_x += daughter->px();
        p_y += daughter->py();
        p_z += daughter->pz();
        // if (inv_M > max_inv_mass) max_inv_mass = inv_M;
        // isol_pi_inv_m_to_ks->Fill(inv_M);
      }
      inv_M += sqrt(pow(E, 2) - pow(p_x, 2) - pow(p_y, 2) - pow(p_z, 2));
      ks_inv_m_pi->Fill(inv_M);

      /* or
        for (size_t index = 0; index != vec.size(); ++index)
          // do something with vec[index]

        // as of C++11
        for (const auto& item: vec)
          // do something with item
      */
    
      //int a = (*Vertices)[i].daughters()->size();
      //reco::CompositeCandidate a = (*Vertices)[i].daughters;

    }
    // Fill the variables
    h_v0_count->Fill(v0_count);
  }
  else if (!Vertices.isValid() && !debug) cout << "UNVALID VERTICES" << endl;

  /// HPS pi0's - with loop among reco::tau
  if (Strips.isValid() && false) 
  {
    if (Strips->size() > 0) if (!debug) cout << "Number of HPS Pizeros = " << Strips->size() << std::endl;
    if (pf_taus->size() > 0)  if (!debug) cout << "Number of RECO Tau = " << pf_taus->size() << endl;

    for (unsigned int i = 0; i < Strips->size(); i++) 
    {
      if (!debug) cout << "\tPi0_" << i << " (" << (*Strips)[i].charge() <<") " <<   
            (*Strips)[i].px() << " " <<
            (*Strips)[i].py() << " " <<
            (*Strips)[i].pz() << " " <<
            (*Strips)[i].p()  << " " <<
            endl;
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
      if (!debug) cout << " PI0 " << i 
                        << "  px = " << pizero_px[pizero_count]
                        << "  py = " << pizero_py[pizero_count]
                        << "  pz = " << pizero_pz[pizero_count]
                        << "  vx = " << pizero_x[pizero_count]
                        << "  vy = " << pizero_y[pizero_count]
                        << "  vz = " << pizero_z[pizero_count] << endl;

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
            const std::vector < reco::RecoTauPiZero > tau_pizeros = pftauref->signalPiZeroCandidates();
            if (!debug) cout << "\t\t\tsignal tau_pizeros: " << tau_pizeros.size() << endl;
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

                  if (!debug) cout << "\t\t\t\tPi0_" << j_pi << " (" << tau_pizeros[j_pi].charge() <<") " <<   
                    tau_pizeros[j_pi].px() << " " <<
                    tau_pizeros[j_pi].py() << " " <<
                    tau_pizeros[j_pi].pz() << " " <<
                    tau_pizeros[j_pi].energy()  << " " <<
                    endl;

                  break;
                }
              }

          //Isol Pi0
            const std::vector < reco::RecoTauPiZero > tau_pizeros_isol = pftauref->isolationPiZeroCandidates();
            if (!debug) cout << "\t\t\tisolation tau_pizeros_isol :" << tau_pizeros_isol.size() << endl;
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

                   if (!debug) cout << "\t\t\t\tPi0_" << j_pi << " (" << tau_pizeros_isol[j_pi].charge() <<") " <<   
                    tau_pizeros_isol[j_pi].px() << " " <<
                    tau_pizeros_isol[j_pi].py() << " " <<
                    tau_pizeros_isol[j_pi].pz() << " " <<
                    tau_pizeros_isol[j_pi].energy()  << " " <<
                    endl;
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
      else if (!debug) cout <<"ALL PIONS MATCHED" << endl;
  } 

  /// RECO Taus
  if (pf_taus.isValid() )
  {
    if (pf_taus->size() > 0 && !debug)  cout << "Number of Tau = " << pf_taus->size() << std::endl;
    for (unsigned int i = 0; i < pf_taus->size(); i++) // Over Tau's
    {
      if (!debug) cout << "\ttau #" << i << " from " << pf_taus->size() ;//<< endl;
      //reco::PFTau pfTau(pf_taus, i); //same as (*pf_taus)[i]
      //edm::AtomicPtrCache< std::vector< reco::RecoTauPiZero > >
      //const std::vector< reco::RecoTauPiZero > a  = (*pf_taus)[i].signalPiZeroCandidates();
      //(*reco::PFTau)pfTau->signalPiZeroCandidates();
      //const std::vector < RecoTauPiZero > &   isolationPiZeroCandidates () const
      //const std::vector< RecoTauPiZero > & reco::PFTau::isolationPiZeroCandidates (   ) const
      //const std::vector < RecoTauPiZero > &   signalPiZeroCandidates () const
      reco::PFTauRef pftauref(pf_taus, i); // one Tau instance, typedef edm::Ref<PFTauCollection> reco::PFTauRef


      // Signal pi0's
      if (false)
      {
        const std::vector < reco::RecoTauPiZero > tau_pizeros = pftauref->signalPiZeroCandidates();
        if (debug) cout << "\ttau_pizeros num:" << tau_pizeros.size() << endl;
        num_pions->Fill(tau_pizeros.size());
        double inv_M = 0;
        double p_x(0), p_y(0), p_z(0), E(0);
        if (tau_pizeros.size() > 0)
          for (unsigned int j = 0; j < tau_pizeros.size(); j++)
          {
            if (!debug) cout << "\t\tPi0_" << j << " (" << tau_pizeros[j].charge() <<") " <<   
              tau_pizeros[j].px() << " " <<
              tau_pizeros[j].py() << " " <<
              tau_pizeros[j].pz() << " " <<
              tau_pizeros[j].energy()  << " " <<
              endl;
              p_x += tau_pizeros[j].px();
              p_y += tau_pizeros[j].py();
              p_z += tau_pizeros[j].pz();
              E += tau_pizeros[j].energy();
          }
        inv_M = sqrt(pow(E, 2) - pow(p_x, 2) - pow(p_y, 2) - pow(p_z, 2));
        if (!debug) cout << "TAU_" << i << " WITH INVARIANT MASS OF PIONS:" <<  inv_M << endl;
        if (false && inv_M > 0.140) 
        {
          if (debug) cout << "TAU_" << i << " WITH INVARIANT MASS OF PIONS:" <<  inv_M << endl;
          pions_inv_m->Fill(inv_M);
          num_pion_res++;
        }

        // Fill arrays with tau's properties
        {
          tau_px[tau_count] = (*pf_taus)[i].px();
          tau_py[tau_count] = (*pf_taus)[i].py();
          tau_pz[tau_count] = (*pf_taus)[i].pz();
          tau_e[tau_count]  = (*pf_taus)[i].p();
          tau_pt[tau_count] = (*pf_taus)[i].pt();
          tau_eta[tau_count] = (*pf_taus)[i].eta();
          tau_phi[tau_count] = (*pf_taus)[i].phi();
          tau_x[tau_count] = (*pf_taus)[i].vx();
          tau_y[tau_count] = (*pf_taus)[i].vy();
          tau_z[tau_count] = (*pf_taus)[i].vz();
          tau_count++;

          if (!debug) cout << " Tau " << i 
            << "  px = " << tau_px[tau_count]
            << "  py = " << tau_py[tau_count]
            << "  pz = " << tau_pz[tau_count]
            << "  vx = " << tau_x[tau_count]
            << "  vy = " << tau_y[tau_count]
            << "  vz = " << tau_z[tau_count] << endl;
        }
      }
      
      //Matching Isolation Pions to obtain Ks peak
      const std::vector < reco::RecoTauPiZero > tau_pizeros_isol = pftauref->isolationPiZeroCandidates(); // pions of the considered Tau
      if (tau_pizeros_isol.size() > 1)
      {
        cout << endl;
        for (unsigned int j_pi = 0; j_pi < tau_pizeros_isol.size() - 1; j_pi++) //over isolation pions in Tau
        {
            //Check the vertex position
             if (false && !debug) cout << "\t\t\tPi0_" << j_pi << " (" << tau_pizeros_isol[j_pi].charge() <<"), vertex:  " <<   
              tau_pizeros_isol[j_pi].vx() << " " <<
              tau_pizeros_isol[j_pi].vy() << " " <<
              tau_pizeros_isol[j_pi].vz() << " " <<
              endl; //all from the same vertex

            //Build the inv mass
              if ( !debug) cout << "\t\t===>matching to Pi0_"<< j_pi << " from " << tau_pizeros_isol.size() << " in this tau" << endl;
              for (unsigned int j2_pi = j_pi + 1; j2_pi < tau_pizeros_isol.size(); j2_pi++)
              {
                double E = tau_pizeros_isol[j_pi].energy() +  tau_pizeros_isol[j2_pi].energy();
                double p_x = tau_pizeros_isol[j_pi].px() +  tau_pizeros_isol[j2_pi].px();
                double p_y = tau_pizeros_isol[j_pi].py() +  tau_pizeros_isol[j2_pi].py();
                double p_z = tau_pizeros_isol[j_pi].pz() +  tau_pizeros_isol[j2_pi].pz();
                double inv_M = sqrt(pow(E, 2) - pow(p_x, 2) - pow(p_y, 2) - pow(p_z, 2));
                if (inv_M > max_inv_mass) max_inv_mass = inv_M;
                isol_pi_inv_m_to_ks->Fill(inv_M);
                if ( !debug) cout << "\t\t\tm( Pi0_" << j_pi << " + ";
                if ( !debug) cout <<  "Pi0_" << j2_pi << " ) = " << inv_M << endl;
              }
              cout << endl;

              isol_pi_pt->Fill(tau_pizeros_isol[j_pi].pt());
        }
              isol_pi_pt->Fill(tau_pizeros_isol[tau_pizeros_isol.size() - 1].pt());


      }
      else if (!debug) cout << "\t\t\tto few pions " << endl;

      if (tau_count >= 1000) 
      {
        cerr << "number of taus > 1000. They are missing." << endl; 
        break;
      }
    }
  } 
  else if (!pf_taus.isValid() && !debug ) cout << "no valid pf_taus"  << endl;

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
  isol_pi_inv_m_to_ks->Write();
  ks_daughter_pt->Write();
  ks_inv_m_pi->Write();
  isol_pi_pt->Write();
}

AOD_pi0::~AOD_pi0()
{
  cout << "Total num of not-matched events:" << num_ev_tau_pi_not_in_hps_pi << endl;
  cout << "num_pion_res: " << num_pion_res << endl;
  cout << "max_inv_mass: " << max_inv_mass << endl;
  cout << "max_ks_daughter_pt: " << max_ks_daughter_pt << endl;
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
