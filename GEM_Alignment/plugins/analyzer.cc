#include <assert.h>
#include <memory>
#include <cmath>
#include <iostream>
#include <sstream>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/CSCRecHit/interface/CSCRecHit2D.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include <DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h>
#include <DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigi.h>
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"

#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"


using namespace std;
using namespace edm;


struct MuonData
{
  void init();
  TTree* book(TTree *t, int prop_type);
  //Muon Info//////////////////////////////////////////////////////
  int muon_charge; float muon_pt; float muon_eta; float muon_momentum;
  unsigned long long  evtNum; unsigned long long  lumiBlock; int muonIdx;
  int runNum;
  //Propagation Info//////////////////////////////////////////////////////
  float prop_GP[3]; float prop_LP[3]; float prop_startingPoint_GP[3];
  float prop_yroll; float prop_localphi_rad; float prop_localphi_deg;
  bool has_prop; bool has_fidcut;
  int prop_location[5];
  //Track Info//////////////////////////////////////////////////////
  float track_chi2; float track_ndof; int n_ME11_segment; int which_track;
  int hasME11; int hasME11RecHit; int hasME11A; int hasME11ARecHit;
  int nCSCSeg; int nDTSeg; int nME11RecHits; float ME11_BunchX; int ME11_strip;
  int ME11_location[5]; int inner_or_outer_mom;
  //Rechit Info//////////////////////////////////////////////////////
  float rechit_GP[3]; float rechit_LP[3];
  float rechit_yroll; float rechit_localphi_rad; float rechit_localphi_deg;
  bool has_rechit;
  int rechit_first_strip; int rechit_CLS; int rechit_BunchX;
  float RdPhi; float RdPhi_Corrected; int rechit_detId;
  int nRecHitsTot; int nRecHits5; int nRecHits2;
  int rechit_location[5];
  int nRecHitsRpos1L1; int nRecHitsRpos1L2;
  int nRecHitsRneg1L1; int nRecHitsRneg1L2;
  //Sim info for MC
  float sim_GP[3]; float sim_LP[3];
  float simDy; float sim_yroll; int nSim;
};

void MuonData::init()
{
  //Muon Info//////////////////////////////////////////////////////
  muon_charge = 9999; muon_pt = 9999; muon_eta = 9999; muon_momentum = 9999;
  evtNum = 99999999; lumiBlock = 99999999; muonIdx = 99999999; runNum = 99999999;
  //Propagation Info//////////////////////////////////////////////////////
  for(int i=0; i<3; ++i){
    prop_GP[i] = 99999; prop_LP[i] = 99999; prop_startingPoint_GP[i] = 99999;
  }
  prop_yroll = 99999; prop_localphi_rad = 99999; prop_localphi_deg = 99999;
  has_prop = false; has_fidcut = false;
  for(int i=0; i<5; ++i){
    prop_location[i] = 99999;
  }
  //Track Info//////////////////////////////////////////////////////
  track_chi2 = 999999; track_ndof = 999999; n_ME11_segment = 999999; which_track = 999999;
  hasME11 = 0; hasME11RecHit = 0; hasME11A = 0; hasME11ARecHit = 0;
  nCSCSeg = 999999; nDTSeg = 999999; nME11RecHits = 999999; ME11_BunchX = 999999; ME11_strip = 999999;
  for(int i=0; i<5; ++i){
    ME11_location[i] = 999999;
  }
  inner_or_outer_mom = 99999;
  //Rechit Info//////////////////////////////////////////////////////
  for(int i=0; i<3; ++i){
    rechit_GP[i] = 999999; rechit_LP[i] = 999999;
  }
  rechit_yroll = 999999; rechit_localphi_rad = 999999; rechit_localphi_deg = 999999;
  has_rechit = false;
  rechit_first_strip = 999999; rechit_CLS = 999999; rechit_BunchX = 999999;
  RdPhi = 999999; RdPhi_Corrected = 999999; rechit_detId = 999999;
  nRecHitsTot = 999999; nRecHits5 = 999999; nRecHits2 = 999999;
  for(int i=0; i<5; ++i){
    rechit_location[i] = 999999;
  }
  nRecHitsRpos1L1 = 999999; nRecHitsRpos1L2 = 999999;
  nRecHitsRneg1L1 = 999999; nRecHitsRneg1L2 = 999999;
  //Sim info for MC
  for(int i=0; i<3; ++i){
    sim_GP[i] = 9999999; sim_LP[i] = 9999999;
  }
  simDy = 9999999; sim_yroll = 9999999; nSim = 9999999;
}

TTree* MuonData::book(TTree *t, int prop_type){
  edm::Service< TFileService > fs;
  if(prop_type == 1){
    t = fs->make<TTree>("CSC_Prop", "CSC_Prop");
  }
  else if(prop_type == 2){
    t = fs->make<TTree>("Inner_Prop", "Inner_Prop");
  }
  else if(prop_type == 3){
    t = fs->make<TTree>("ME11Seg_Prop", "ME11Seg_Prop");
  }
  else{
    std::cout << "Bad prop type, failure" << std::endl;
  }
  //Muon Info//////////////////////////////////////////////////////
  t->Branch("muon_charge", &muon_charge); t->Branch("muon_pt", &muon_pt);
  t->Branch("muon_eta", &muon_eta); t->Branch("muon_momentum", &muon_momentum);
  t->Branch("evtNum", &evtNum); t->Branch("lumiBlock", &lumiBlock); t->Branch("muonIdx", &muonIdx);
  t->Branch("runNum", &runNum);
  //Propagation Info//////////////////////////////////////////////////////
  t->Branch("prop_GP", &prop_GP, "prop_GP[3] (x,y,z)/F");
  t->Branch("prop_LP", &prop_LP, "prop_LP[3] (x,y,z)/F");
  t->Branch("prop_startingPoint_GP", &prop_startingPoint_GP, "prop_startingPoint_GP[3] (x,y,z)/F");
  t->Branch("prop_yroll", &prop_yroll);
  t->Branch("prop_localphi_rad", &prop_localphi_rad);
  t->Branch("prop_localphi_deg", &prop_localphi_deg);
  t->Branch("has_prop", &has_prop);
  t->Branch("has_fidcut", &has_fidcut);
  t->Branch("prop_location", &prop_location, "prop_location[5] (reg, sta, cha, lay, rol)/I");
  //Track Info//////////////////////////////////////////////////////
  t->Branch("track_chi2", &track_chi2); t->Branch("track_ndof", &track_ndof);
  t->Branch("n_ME11_segment", &n_ME11_segment); t->Branch("which_track", &which_track);
  t->Branch("hasME11", &hasME11); t->Branch("hasME11RecHit", &hasME11RecHit);
  t->Branch("hasME11A", &hasME11A); t->Branch("hasME11ARecHit", &hasME11ARecHit);
  t->Branch("nCSCSeg", &nCSCSeg); t->Branch("nDTSeg", &nDTSeg);
  t->Branch("nME11RecHits", &nME11RecHits); t->Branch("ME11_BunchX", &ME11_BunchX);
  t->Branch("ME11_strip", &ME11_strip);
  t->Branch("ME11_location", &ME11_location, "ME11_location[5] (end, sta, ring, cha, lay)/I");
  t->Branch("inner_or_outer_mom", &inner_or_outer_mom, "inner_or_outer_mom (0 = inner, 1 = outer)/I");
  //Rechit Info//////////////////////////////////////////////////////
  t->Branch("rechit_GP", &rechit_GP, "rechit_GP[3] (x,y,z)/F");
  t->Branch("rechit_LP", &rechit_LP, "rechit_LP[3] (x,y,z)/F");
  t->Branch("rechit_yroll", &rechit_yroll);
  t->Branch("rechit_localphi_rad", &rechit_localphi_rad);
  t->Branch("rechit_localphi_deg", &rechit_localphi_deg);
  t->Branch("has_rechit", &has_rechit);
  t->Branch("rechit_first_strip", &rechit_first_strip);
  t->Branch("rechit_CLS", &rechit_CLS);
  t->Branch("rechit_BunchX", &rechit_BunchX);
  t->Branch("RdPhi", &RdPhi);
  t->Branch("RdPhi_Corrected", &RdPhi_Corrected);
  t->Branch("rechit_detId", &rechit_detId);
  t->Branch("nRecHitsTot", &nRecHitsTot);
  t->Branch("nRecHits2", &nRecHits2);
  t->Branch("nRecHits5", &nRecHits5);
  t->Branch("rechit_location", &rechit_location, "rechit_location[5] (reg, sta, cha, lay, rol)/I");
  t->Branch("nRecHitsRpos1L1", &nRecHitsRpos1L1);
  t->Branch("nRecHitsRpos1L2", &nRecHitsRpos1L2);
  t->Branch("nRecHitsRneg1L1", &nRecHitsRneg1L1);
  t->Branch("nRecHitsRneg1L2", &nRecHitsRneg1L2);
  //Sim info for MC
  t->Branch("sim_GP", &sim_GP, "sim_GP[3] (x,y,z)/F");
  t->Branch("sim_LP", &sim_LP, "sim_LP[3] (x,y,z)/F");
  t->Branch("simDy", &simDy);
  t->Branch("sim_yroll", &sim_yroll);
  t->Branch("nSim", &nSim);
  return t;
}


class analyzer : public edm::one::EDAnalyzer<> {
public:
  explicit analyzer(const edm::ParameterSet&);
  ~analyzer(){};

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob() ;
  virtual void endJob() ;

  void propagate(const reco::Muon* mu, int prop_type, const edm::Event& iEvent, int i);
  void CSCSegmentCounter(const reco::Muon* mu, MuonData& data_);
  void propagate_to_GEM(const reco::Muon* mu, const GEMEtaPartition* ch, int prop_type, bool &tmp_has_prop, GlobalPoint &pos_GP, MuonData& data_);
  void GEM_rechit_matcher(const GEMEtaPartition* ch, LocalPoint prop_LP, MuonData& data_);
  void GEM_simhit_matcher(const GEMEtaPartition* ch, GlobalPoint prop_GP, MuonData& data_);
  float RdPhi_func(float stripAngle, const edm::OwnVector<GEMRecHit, edm::ClonePolicy<GEMRecHit> >::const_iterator rechit, float prop_localx, float prop_localy, const GEMEtaPartition* ch);
  bool fidcutCheck(float local_y, float localphi_deg, const GEMEtaPartition* ch);

  edm::EDGetTokenT<GEMRecHitCollection> gemRecHits_;
  edm::Handle<GEMRecHitCollection> gemRecHits;
  edm::EDGetTokenT<vector<PSimHit> > gemSimHits_;
  edm::Handle<vector<PSimHit> > gemSimHits;
  edm::EDGetTokenT<edm::View<reco::Muon> > muons_;
  edm::EDGetTokenT<CSCSegmentCollection> cscSegments_;

  edm::Service<TFileService> fs;

  MuonServiceProxy* theService_;
  edm::ESHandle<Propagator> propagator_;
  edm::ESHandle<TransientTrackBuilder> ttrackBuilder_;

  ESHandle<GlobalTrackingGeometry> theTrackingGeometry;

  edm::ESHandle<GEMGeometry> GEMGeometry_;
  edm::ESHandle<CSCGeometry> CSCGeometry_;

  bool CSC_prop; bool tracker_prop; bool Segment_prop;
  vector<int> prop_list;
  bool debug;
  bool isCosmic;

  MuonData data_;
  TTree* CSC_tree; TTree* Tracker_tree; TTree* Segment_tree;
  TH2D* nME11_col_vs_matches = new TH2D("nME11_test", "nME11_test", 5, 0, 5, 5, 0, 5);

  bool isMC;
  const CSCSegment *ME11_segment;

  const edm::ESGetToken<GEMGeometry, MuonGeometryRecord> gemGeomToken_;
  const edm::ESGetToken<CSCGeometry, MuonGeometryRecord> cscGeomToken_;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttkToken_;
  const edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> geomToken_;
};


analyzer::analyzer(const edm::ParameterSet& iConfig)
  : gemGeomToken_(esConsumes()),
    cscGeomToken_(esConsumes()),
    ttkToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
    geomToken_(esConsumes())
{
  cout << "Begin analyzer" << endl;
  edm::ParameterSet serviceParameters = iConfig.getParameter<edm::ParameterSet>("ServiceParameters");
  theService_ = new MuonServiceProxy(serviceParameters, consumesCollector());

  muons_ = consumes<View<reco::Muon> >(iConfig.getParameter<InputTag>("muons"));
  gemRecHits_ = consumes<GEMRecHitCollection>(iConfig.getParameter<edm::InputTag>("gemRecHits"));
  gemSimHits_ = consumes<vector<PSimHit> >(iConfig.getParameter<edm::InputTag>("gemSimHits"));
  cscSegments_ = consumes<CSCSegmentCollection>(edm::InputTag("cscSegments"));


  tracker_prop = iConfig.getParameter<bool>("tracker_prop");
  CSC_prop = iConfig.getParameter<bool>("CSC_prop");
  Segment_prop = iConfig.getParameter<bool>("Segment_prop");
  debug = iConfig.getParameter<bool>("debug");
  isCosmic = iConfig.getParameter<bool>("isCosmic");
  cout << "tracker_prop " << tracker_prop << " CSC_prop " << CSC_prop << " debug " << debug << std::endl;

  if(CSC_prop){CSC_tree = data_.book(CSC_tree, 1); prop_list.push_back(1);}
  if(tracker_prop){Tracker_tree = data_.book(Tracker_tree, 2); prop_list.push_back(2);}
  if(Segment_prop){Segment_tree = data_.book(Segment_tree, 3); prop_list.push_back(3);}
}


void
analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  //iSetup.get<MuonGeometryRecord>().get(GEMGeometry_);
  //iSetup.get<MuonGeometryRecord>().get(CSCGeometry_);
  //iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttrackBuilder_);
  //iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);

  GEMGeometry_ = &iSetup.getData(gemGeomToken_);
  CSCGeometry_ = &iSetup.getData(cscGeomToken_);
  ttrackBuilder_ = &iSetup.getData(ttkToken_);
  theTrackingGeometry = &iSetup.getData(geomToken_);

  theService_->update(iSetup);
  auto propagator = theService_->propagator("SteppingHelixPropagatorAny");

  isMC = false;
  if (! iEvent.eventAuxiliary().isRealData()) isMC = true;
  iEvent.getByToken(gemRecHits_, gemRecHits);
  if (isMC) {
    iEvent.getByToken(gemSimHits_, gemSimHits);
  }
  edm::Handle<View<reco::Muon> > muons;
  if (! iEvent.getByToken(muons_, muons)) return;
  if (muons->size() == 0) return;

  edm::Handle<CSCSegmentCollection> cscSegments;
  if (! iEvent.getByToken(cscSegments_, cscSegments)){std::cout << "Bad segments" << std::endl;}

  if (debug) cout << "New! EvtNumber = " << iEvent.eventAuxiliary().event() << " LumiBlock = " << iEvent.eventAuxiliary().luminosityBlock() << " RunNumber = " << iEvent.run() << endl;

  for (size_t i = 0; i < muons->size(); ++i){
    edm::RefToBase<reco::Muon> muRef = muons->refAt(i);
    const reco::Muon* mu = muRef.get();
    //if (mu->pt() < 2.0) continue;  //can apply a pt cut later
    if (not mu->standAloneMuon()) continue;
    if (debug) cout << "new standalone" << endl;
    for(auto it = std::begin(prop_list); it != std::end(prop_list); ++it){
      if (debug) std::cout << "prop " << *it << "about to start propagate" << std::endl;
      int prop_type = *it;
      propagate(mu, prop_type, iEvent, i);
    }
  }

  if (debug){
    int muon_size = muons->size();
    int CSCSegments_size = cscSegments->size();
    int muon_STA_counter = 0;
    int ME11_counter = 0;
    int Matches_counter_norepeat = 0;
    int ME11_counter_norepeat = 0;
    int tot_muon_ME11_counter = 0;

    //Commented out for github, may be many print statements
    //cout << "Starting loop over all gemRecHits " << gemRecHits->size() << std::endl;
    //for(GEMRecHitCollection::const_iterator gemRecHit = gemRecHits->begin(); gemRecHit != gemRecHits->end(); gemRecHit++){
    //  cout << gemRecHit->gemId() << endl;
    //}

    std::cout << "MUON MATCHES" << std::endl;
    std::vector<int> list_of_matches = {};
    for (size_t i = 0; i < muons->size(); ++i){
      edm::RefToBase<reco::Muon> muRef = muons->refAt(i);
      const reco::Muon* mu = muRef.get();
      if (mu->isStandAloneMuon()){
        muon_STA_counter ++;
        int muon_ME11_counter_tmp = 0;
        auto matches = mu->matches();
        for (auto MCM : matches){
          if(MCM.detector() != 2) continue;
          for(auto MSM : MCM.segmentMatches){
            auto cscSegRef = MSM.cscSegmentRef;
            auto cscDetID = cscSegRef->cscDetId();
            std::cout << cscDetID.endcap() << cscDetID.station() << cscDetID.ring() << cscDetID.chamber() << std::endl;
            if(cscDetID.station() == 1 and (cscDetID.ring() == 1 or cscDetID.ring() == 4)){
              muon_ME11_counter_tmp++;
              int tmp_ID = cscDetID.endcap()*10000+cscDetID.station()*1000+cscDetID.ring()*100+cscDetID.chamber();
              if (std::find(list_of_matches.begin(), list_of_matches.end(), tmp_ID) == list_of_matches.end()){
                list_of_matches.push_back(tmp_ID);
                std::cout << "New Det Triggered " << tmp_ID << std::endl;
                Matches_counter_norepeat++;
              }
            }
          }
        }
        std::cout << "Muon had " << muon_ME11_counter_tmp << " segments" << std::endl;
        tot_muon_ME11_counter += muon_ME11_counter_tmp;
      }
    }

    std::cout << "SEGMENT COLLECTION" << std::endl;
    std::vector<int> list_of_collection = {};
    for(CSCSegmentCollection::const_iterator segmentCSC = cscSegments->begin(); segmentCSC != cscSegments->end(); segmentCSC++){
      auto cscDetID = (segmentCSC)->cscDetId();
      std::cout << cscDetID.endcap() << cscDetID.station() << cscDetID.ring() << cscDetID.chamber() << std::endl;
      if(cscDetID.station() == 1 and (cscDetID.ring() == 1 or cscDetID.ring() == 4)){
        ME11_counter++;
        int tmp_ID = cscDetID.endcap()*10000+cscDetID.station()*1000+cscDetID.ring()*100+cscDetID.chamber();
        if (std::find(list_of_collection.begin(), list_of_collection.end(), tmp_ID) == list_of_collection.end()){
          list_of_collection.push_back(tmp_ID);
          std::cout << "New Det Triggered " << tmp_ID << std::endl;
          ME11_counter_norepeat++;
        }
      }
    }
  std::cout << "Muon size           = " << muon_size << std::endl;
  std::cout << "Segment Size        = " << CSCSegments_size << std::endl;
  std::cout << "muon STA  counter   = " << muon_STA_counter << std::endl;
  std::cout << "Muon Matches ME11   = " << tot_muon_ME11_counter << std::endl;
  std::cout << "Seg Collection ME11 = " << ME11_counter << std::endl;
  std::cout << "Muon Matches no rep = " << Matches_counter_norepeat << std::endl;
  std::cout << "Coll Matches no rep = " << ME11_counter_norepeat << std::endl;
  nME11_col_vs_matches->Fill(Matches_counter_norepeat, ME11_counter_norepeat);
  }
}


float analyzer::RdPhi_func(float stripAngle, const edm::OwnVector<GEMRecHit, edm::ClonePolicy<GEMRecHit> >::const_iterator rechit, float prop_localx, float prop_localy, const GEMEtaPartition* ch){
  GEMDetId gemid((rechit)->geographicalId());
  const auto& etaPart = GEMGeometry_->etaPartition(gemid);
  const auto& etaPart_ch = GEMGeometry_->etaPartition(ch->id());
  float deltay_roll =  etaPart_ch->toGlobal(etaPart_ch->centreOfStrip(etaPart_ch->nstrips()/2)).perp() - etaPart->toGlobal(etaPart->centreOfStrip(etaPart->nstrips()/2)).perp();
  return cos(stripAngle) * (prop_localx - (rechit)->localPosition().x()) - sin(stripAngle) * (prop_localy + deltay_roll);
}
void analyzer::CSCSegmentCounter(const reco::Muon* mu, MuonData& data_){
  const reco::Track* Track = mu->outerTrack().get();
  int tmp_CSC_counter = 0; int tmp_DT_counter = 0; int tmp_ME11_counter = 0; int tmp_ME11RecHit_counter = 0; float tmp_ME11_BunchX = 99999; int tmp_ME11_strip = 99999; bool tmp_hasME11A = 0;
  if(isCosmic){
    tmp_CSC_counter = mu->numberOfSegments(1,2) + mu->numberOfSegments(2,2) + mu->numberOfSegments(3,2) + mu->numberOfSegments(4,2);
    tmp_DT_counter = mu->numberOfSegments(1,1) + mu->numberOfSegments(2,1) + mu->numberOfSegments(3,1) + mu->numberOfSegments(4,1);
    auto matches = mu->matches();
    for (auto MCM : matches){
      if(MCM.detector() != 2) continue;
      for(auto MSM : MCM.segmentMatches){
        auto cscSegRef = MSM.cscSegmentRef;
        auto cscDetID = cscSegRef->cscDetId();
        if(cscDetID.station() == 1 and (cscDetID.ring() == 1 or cscDetID.ring() == 4)){
          if(cscDetID.ring() == 4){tmp_hasME11A = 1;}
          tmp_ME11_counter++;
          if (debug){std::cout << "isCosmic = True! Getting ME11 Segment" << std::endl;}
          ME11_segment = cscSegRef.get();
          tmp_ME11RecHit_counter = (cscSegRef.get())->nRecHits(); // Find the real function for this. Bad if multiple segments.
          tmp_ME11_BunchX = ME11_segment->time();
          auto cscDetID_FAKE = CSCDetId(cscDetID.endcap(), cscDetID.station(), cscDetID.ring(), cscDetID.chamber(), 3);
          const CSCLayer* tmp_ME11_layer = CSCGeometry_->layer(cscDetID_FAKE);
          const CSCLayerGeometry* tmp_ME11_layer_geo = tmp_ME11_layer->geometry();
          tmp_ME11_strip = tmp_ME11_layer_geo->nearestStrip(ME11_segment->localPosition());
          data_.ME11_location[0] = cscDetID.endcap(); data_.ME11_location[1] = cscDetID.station(); data_.ME11_location[2] = cscDetID.ring(); data_.ME11_location[3] = cscDetID.chamber(); data_.ME11_location[4] = cscDetID.layer();
        }
      }
    }
  }
  else{
    for (size_t RecHit_iter = 0; RecHit_iter != Track->recHitsSize(); RecHit_iter++){
      const TrackingRecHit* RecHit = (Track->recHit(RecHit_iter)).get();
      DetId RecHitId = RecHit->geographicalId();
      uint16_t RecHitDetId = RecHitId.det();
      if (RecHitDetId == DetId::Muon){
        uint16_t RecHitSubDet = RecHitId.subdetId();
        if (RecHitSubDet == (uint16_t)MuonSubdetId::CSC){
          if (CSCDetId(RecHitId).station() == 1 and CSCDetId(RecHitId).ring() == 1 and RecHit->dimension() == 4){
            tmp_ME11_counter++;
            if (debug){std::cout << "isCosmic = False! Getting ME11 Segment" << std::endl;}
            RecSegment* Rec_segment = (RecSegment*)RecHit;
            ME11_segment = (CSCSegment*)Rec_segment;
            tmp_ME11_BunchX = ((CSCRecHit2D*)RecHit)->wgroupsBX();
            auto cscDetID_FAKE = CSCDetId(CSCDetId(RecHitId).endcap(), CSCDetId(RecHitId).station(), CSCDetId(RecHitId).ring(), CSCDetId(RecHitId).chamber(), 3);
            const CSCLayer* tmp_ME11_layer = CSCGeometry_->layer(cscDetID_FAKE);
            const CSCLayerGeometry* tmp_ME11_layer_geo = tmp_ME11_layer->geometry();
            tmp_ME11_strip = tmp_ME11_layer_geo->nearestStrip(ME11_segment->localPosition());
            data_.ME11_location[0] = CSCDetId(RecHitId).endcap(); data_.ME11_location[1] = CSCDetId(RecHitId).station(); data_.ME11_location[2] = CSCDetId(RecHitId).ring(); data_.ME11_location[3] = CSCDetId(RecHitId).chamber(); data_.ME11_location[4] = CSCDetId(RecHitId).layer();
          }
          if (CSCDetId(RecHitId).station() == 1 and CSCDetId(RecHitId).ring() == 1){tmp_ME11RecHit_counter++;}
          if (RecHit->dimension() == 4){tmp_CSC_counter++;}
        }
        if (RecHitSubDet == (uint16_t)MuonSubdetId::DT){
          if (RecHit->dimension() > 1){tmp_DT_counter++;}
        }
      }
    }
  }
  data_.nCSCSeg = tmp_CSC_counter; data_.nDTSeg = tmp_DT_counter;
  data_.n_ME11_segment = tmp_ME11_counter;
  data_.nME11RecHits = tmp_ME11RecHit_counter;
  data_.ME11_BunchX = tmp_ME11_BunchX;
  data_.ME11_strip = tmp_ME11_strip;
  data_.hasME11A =  tmp_hasME11A;
  if(data_.n_ME11_segment >= 1 and data_.n_ME11_segment < 1000){data_.hasME11 = 1;}
}
void analyzer::propagate_to_GEM(const reco::Muon* mu, const GEMEtaPartition* ch, int prop_type, bool &tmp_has_prop, GlobalPoint &pos_GP, MuonData& data_){
  const reco::Track* Track;
  reco::TransientTrack ttrack;
  tmp_has_prop = false;
  int tmp_inner_or_outer_mom = 99999;
  const BoundPlane& bps(ch->surface());
  auto propagator = theService_->propagator("SteppingHelixPropagatorAny");
  const auto& etaPart_ch = GEMGeometry_->etaPartition(ch->id());
  TrajectoryStateOnSurface tsos_ch; TrajectoryStateOnSurface tsos_seg;
  GlobalPoint pos_startingPoint_GP;
  if(prop_type == 1 or prop_type == 2){
    if(prop_type == 1){
      Track = mu->outerTrack().get();
      ttrack = ttrackBuilder_->build(Track);
    }
    if(prop_type == 2){
      Track = mu->track().get();
      ttrack = ttrackBuilder_->build(Track);
    }
    float inner_delta = abs(ttrack.innermostMeasurementState().globalPosition().z() - GEMGeometry_->etaPartition(ch->id())->toGlobal(etaPart_ch->centreOfStrip(etaPart_ch->nstrips()/2)).z());
    float outer_delta = abs(ttrack.outermostMeasurementState().globalPosition().z() - GEMGeometry_->etaPartition(ch->id())->toGlobal(etaPart_ch->centreOfStrip(etaPart_ch->nstrips()/2)).z());
    float used_delta = 0;
    if (inner_delta < outer_delta){
      tsos_seg = ttrack.innermostMeasurementState(); tsos_ch = propagator->propagate(tsos_seg, ch->surface()); used_delta = inner_delta;
      if(prop_type == 1){data_.which_track = 1;}
      else{data_.which_track = 0;}
    }
    else{
      tsos_seg = ttrack.outermostMeasurementState(); tsos_ch = propagator->propagate(tsos_seg, ch->surface()); used_delta = outer_delta;
      if(prop_type == 1){data_.which_track = 0;}
      else{data_.which_track = 1;}
    }
    if (tsos_ch.isValid()){
      const LocalPoint pos_local_ch = ch->toLocal(tsos_ch.globalPosition());
      const LocalPoint pos2D_local_ch(pos_local_ch.x(), pos_local_ch.y(), 0);
      if (!(tsos_ch.globalPosition().z() * tsos_seg.globalPosition().z() < 0) and bps.bounds().inside(pos2D_local_ch) and ch->id().station() == 1 and ch->id().ring() == 1){
        tmp_has_prop = true;
        if(debug){std::cout << "Delta to GEM!!! = " << used_delta << " prop " << prop_type << std::endl;}
        pos_GP = tsos_ch.globalPosition();
        pos_startingPoint_GP = tsos_seg.globalPosition();
      }
    }
  }
  if(prop_type == 3){
    Track = mu->track().get();
    DetId segDetId = ME11_segment->geographicalId();
    const GeomDet* segDet = theTrackingGeometry->idToDet(segDetId);
    LocalVector momentum_at_surface = ME11_segment->localDirection(); //No momentum for segments
    if (Track != 0){
      tmp_inner_or_outer_mom = 0;
      momentum_at_surface = momentum_at_surface*(Track->outerP()); //If innerTrack exists, use momentum
    }
    else{
      Track = mu->outerTrack().get();
      if (Track != 0){
        tmp_inner_or_outer_mom = 1;
        momentum_at_surface = momentum_at_surface*(Track->outerP()); //If no innerTrack, use outerTrack
      }
      else{
        if(debug){std::cout << "No tracks!" << std::endl;}
        return;
      }
    }
    LocalTrajectoryParameters param(ME11_segment->localPosition(), momentum_at_surface, mu->charge());
    AlgebraicSymMatrix mat(5,0);
    mat = ME11_segment->parametersError().similarityT( ME11_segment->projectionMatrix() );
    LocalTrajectoryError error(asSMatrix<5>(mat));
    TrajectoryStateOnSurface tsos_seg(param, error, segDet->surface(), &*theService_->magneticField());
    TrajectoryStateOnSurface tsos_ch = propagator->propagate(tsos_seg, ch->surface());
    if (tsos_ch.isValid()){
      const LocalPoint pos_local_ch = ch->toLocal(tsos_ch.globalPosition());
      const LocalPoint pos2D_local_ch(pos_local_ch.x(), pos_local_ch.y(), 0);
      if (!(tsos_ch.globalPosition().z() * tsos_seg.globalPosition().z() < 0) and bps.bounds().inside(pos2D_local_ch) and ch->id().station() == 1 and ch->id().ring() == 1){
        tmp_has_prop = true;
        pos_GP = tsos_ch.globalPosition();
        pos_startingPoint_GP = tsos_seg.globalPosition();
      }
    }
  }
  if(tmp_has_prop){
    const auto& etaPart_ch = GEMGeometry_->etaPartition(ch->id());
    const float prop_y_to_center = etaPart_ch->toGlobal(etaPart_ch->centreOfStrip(etaPart_ch->nstrips()/2)).perp(); //y distance to the current eta part
    const float prop_y_to_chamber = (GEMGeometry_->chamber(ch->id()))->toLocal(etaPart_ch->toGlobal(etaPart_ch->centreOfStrip(etaPart_ch->nstrips()/2))).y();
    LocalPoint tmp_prop_LP = ch->toLocal(pos_GP);
    data_.prop_GP[0] = pos_GP.x(); data_.prop_GP[1] = pos_GP.y(); data_.prop_GP[2] = pos_GP.z();
    data_.prop_LP[0] = tmp_prop_LP.x(); data_.prop_LP[1] = tmp_prop_LP.y() + prop_y_to_chamber; data_.prop_LP[2] = tmp_prop_LP.z();
    data_.prop_startingPoint_GP[0] = pos_startingPoint_GP.x(); data_.prop_startingPoint_GP[1] = pos_startingPoint_GP.y(); data_.prop_startingPoint_GP[2] = pos_startingPoint_GP.z();
    data_.prop_yroll = tmp_prop_LP.y();
    LocalPoint local_to_center(tmp_prop_LP.x(), tmp_prop_LP.y() + prop_y_to_center, 0);
    float local_phi = local_to_center.phi();
    data_.prop_localphi_rad = (3.14159265/2.) - local_phi;
    data_.prop_localphi_deg = ((3.14159265/2.) - local_phi)*(180./3.14159265);
    data_.has_prop = tmp_has_prop;
    data_.has_fidcut = fidcutCheck(tmp_prop_LP.y(), ((3.14159265/2.) - local_phi)*(180./3.14159265), ch);
    data_.prop_location[0] = ch->id().region(); data_.prop_location[1] = ch->id().station(); data_.prop_location[2] = ch->id().chamber(); data_.prop_location[3] = ch->id().layer(); data_.prop_location[4] = ch->id().roll();
    data_.inner_or_outer_mom = tmp_inner_or_outer_mom;
  }
}
void analyzer::GEM_rechit_matcher(const GEMEtaPartition* ch, LocalPoint prop_LP, MuonData& data_){
  float tmp_rechit_GP_x; float tmp_rechit_GP_y; float tmp_rechit_GP_z;
  float tmp_rechit_LP_x; float tmp_rechit_LP_y; float tmp_rechit_LP_z;
  float tmp_rechit_yroll; float tmp_rechit_localphi_rad; float tmp_rechit_localphi_deg;
  bool tmp_has_rechit = false;
  int tmp_rechit_first_strip; int tmp_rechit_CLS; int tmp_rechit_BunchX;
  float tmp_RdPhi = 9999.; float tmp_RdPhi_Corrected; int tmp_rechit_detId;
  int tmp_nRecHitsTot = 0; int tmp_nRecHits5 = 0; int tmp_nRecHits2 = 0;
  int tmp_rechit_region; int tmp_rechit_station; int tmp_rechit_chamber; int tmp_rechit_layer; int tmp_rechit_roll;
  int tmp_nRecHitsRpos1L1 = 0; int tmp_nRecHitsRpos1L2 = 0; int tmp_nRecHitsRneg1L1 = 0; int tmp_nRecHitsRneg1L2 = 0;
  for(auto hit = gemRecHits->begin(); hit != gemRecHits->end(); hit++){
    if((hit)->geographicalId().det() == DetId::Detector::Muon && (hit)->geographicalId().subdetId() == MuonSubdetId::GEM){
      GEMDetId gemid((hit)->geographicalId());
      if(gemid.region() == 1){
        if(gemid.layer() == 1){tmp_nRecHitsRpos1L1++;}
        if(gemid.layer() == 2){tmp_nRecHitsRpos1L2++;}
      }
      if(gemid.region() == -1){
        if(gemid.layer() == 1){tmp_nRecHitsRneg1L1++;}
        if(gemid.layer() == 2){tmp_nRecHitsRneg1L2++;}
      }
      if(gemid.station() == ch->id().station() and gemid.chamber() == ch->id().chamber() and gemid.layer() == ch->id().layer() and abs(gemid.roll() - ch->id().roll()) <= 1 and gemid.region() == ch->id().region()){
        const auto& etaPart = GEMGeometry_->etaPartition(gemid);
        float strip = etaPart->strip(hit->localPosition());
        float stripAngle = etaPart->specificTopology().stripAngle(strip);
        float rechit_y_to_center = etaPart->toGlobal(etaPart->centreOfStrip(etaPart->nstrips()/2)).perp();
        float rechit_y_to_chamber = (GEMGeometry_->chamber(ch->id()))->toLocal(etaPart->toGlobal(etaPart->centreOfStrip(etaPart->nstrips()/2))).y();
        LocalPoint local_to_center((hit)->localPosition().x(), rechit_y_to_center + (hit)->localPosition().y(), 0);
        if (ch->id().station() == 1 and ch->id().ring() == 1 and fabs((hit)->localPosition().x() - prop_LP.x()) < 999.0){
          tmp_nRecHitsTot++;
          if(abs(RdPhi_func(stripAngle, hit, prop_LP.x(), prop_LP.y(), ch)) < 5){tmp_nRecHits5++;}
          if(abs(RdPhi_func(stripAngle, hit, prop_LP.x(), prop_LP.y(), ch)) < 2){tmp_nRecHits2++;}
          if(abs(tmp_RdPhi) > abs(RdPhi_func(stripAngle, hit, prop_LP.x(), prop_LP.y(), ch))){
            tmp_rechit_GP_x = etaPart->toGlobal((hit)->localPosition()).x(); tmp_rechit_GP_y = etaPart->toGlobal((hit)->localPosition()).y(); tmp_rechit_GP_z = etaPart->toGlobal((hit)->localPosition()).z();
            tmp_rechit_LP_x = (hit)->localPosition().x(); tmp_rechit_LP_y = rechit_y_to_chamber + (hit)->localPosition().y(); tmp_rechit_LP_z = (hit)->localPosition().z();
            tmp_rechit_yroll = (hit)->localPosition().y();
            float local_phi = local_to_center.phi();
            tmp_rechit_localphi_rad = (3.14159265/2.) - local_phi;
            tmp_rechit_localphi_deg = ((3.14159265/2.) - local_phi)*(180./3.14159265);
            tmp_has_rechit = true;
            tmp_rechit_first_strip = (hit)->firstClusterStrip();
            tmp_rechit_CLS = (hit)->clusterSize();
            tmp_rechit_BunchX = (hit)->BunchX();
            tmp_RdPhi = RdPhi_func(stripAngle, hit, prop_LP.x(), prop_LP.y(), ch);
            tmp_RdPhi_Corrected = tmp_RdPhi;
            if((gemid.region() == 1 and gemid.chamber()%2 == 1) || (gemid.region() == -1 && gemid.chamber()%2 == 0)){
              tmp_RdPhi_Corrected = -1.0*tmp_RdPhi_Corrected;
            }
            tmp_rechit_detId = gemid.region()*(gemid.station()*100 + gemid.chamber());
            tmp_rechit_region = gemid.region(); tmp_rechit_station = gemid.station(); tmp_rechit_chamber = gemid.chamber(); tmp_rechit_layer = gemid.layer(); tmp_rechit_roll = gemid.roll();
          }
        }
      }
    }
  }
  if(tmp_has_rechit){
    data_.rechit_GP[0] = tmp_rechit_GP_x; data_.rechit_GP[1] = tmp_rechit_GP_y; data_.rechit_GP[2] = tmp_rechit_GP_z;
    data_.rechit_LP[0] = tmp_rechit_LP_x; data_.rechit_LP[1] = tmp_rechit_LP_y; data_.rechit_LP[2] = tmp_rechit_LP_z;
    data_.rechit_yroll = tmp_rechit_yroll;
    data_.rechit_localphi_rad = tmp_rechit_localphi_rad;
    data_.rechit_localphi_deg = tmp_rechit_localphi_deg;
    data_.has_rechit = tmp_has_rechit;
    data_.rechit_first_strip = tmp_rechit_first_strip;
    data_.rechit_CLS = tmp_rechit_CLS;
    data_.rechit_BunchX = tmp_rechit_BunchX;
    data_.RdPhi = tmp_RdPhi;
    data_.RdPhi_Corrected = tmp_RdPhi_Corrected;
    data_.rechit_detId = tmp_rechit_detId;
    data_.nRecHitsTot = tmp_nRecHitsTot; data_.nRecHits5 = tmp_nRecHits5; data_.nRecHits2 = tmp_nRecHits2;
    data_.rechit_location[0] = tmp_rechit_region; data_.rechit_location[1] = tmp_rechit_station; data_.rechit_location[2] = tmp_rechit_chamber; data_.rechit_location[3] = tmp_rechit_layer; data_.rechit_location[4] = tmp_rechit_roll;
    data_.nRecHitsRpos1L1 = tmp_nRecHitsRpos1L1; data_.nRecHitsRpos1L2 = tmp_nRecHitsRpos1L2; data_.nRecHitsRneg1L1 = tmp_nRecHitsRneg1L1; data_.nRecHitsRneg1L2 = tmp_nRecHitsRneg1L2;
  }
}
void analyzer::GEM_simhit_matcher(const GEMEtaPartition* ch, GlobalPoint prop_GP, MuonData& data_){
  float tmpDy = 999.; float tmpDr = 999.; int tmpSimCounter = 0;
  float tmp_sim_GP_x; float tmp_sim_GP_y; float tmp_sim_GP_z;
  float tmp_sim_LP_x; float tmp_sim_LP_y; float tmp_sim_LP_z;
  bool has_tmp = false;
  for (const auto& simHit:*gemSimHits.product()){
    GEMDetId gemid((simHit).detUnitId());
    if (gemid.station() == ch->id().station() and gemid.chamber() == ch->id().chamber() and gemid.layer() == ch->id().layer() and abs(gemid.roll() - ch->id().roll()) <= 1 and gemid.region() == ch->id().region()){
      tmpSimCounter++;
      const auto& etaPart = GEMGeometry_->etaPartition(gemid);
      float dy = prop_GP.y() - etaPart->toGlobal(simHit.localPosition()).y();
      float dx = prop_GP.x() - etaPart->toGlobal(simHit.localPosition()).x();
      if (dy < tmpDy) tmpDy = dy;
      if (pow(pow(dy, 2) + pow(dx, 2), 0.5) < tmpDr){
        tmp_sim_GP_x = etaPart->toGlobal(simHit.localPosition()).x();
        tmp_sim_GP_y = etaPart->toGlobal(simHit.localPosition()).y();
        tmp_sim_GP_z = etaPart->toGlobal(simHit.localPosition()).z();
        tmp_sim_LP_x = simHit.localPosition().x();
        tmp_sim_LP_y = (GEMGeometry_->chamber(ch->id()))->toLocal(etaPart->toGlobal(etaPart->centreOfStrip(etaPart->nstrips()/2))).y() + simHit.localPosition().y();
        tmp_sim_LP_z = simHit.localPosition().z();
        tmpDr = pow(pow(dy, 2) + pow(dx, 2), 0.5);
        has_tmp = true;
      }
    }
  }
  if(has_tmp){
    data_.sim_GP[0] = tmp_sim_GP_x; data_.sim_GP[1] = tmp_sim_GP_y; data_.sim_GP[2] = tmp_sim_GP_z;
    data_.sim_LP[0] = tmp_sim_LP_x; data_.sim_LP[1] = tmp_sim_LP_y; data_.sim_LP[2] = tmp_sim_LP_z;
    data_.simDy = (tmpDy);
    data_.nSim = (tmpSimCounter);
  }
}
void analyzer::propagate(const reco::Muon* mu, int prop_type, const edm::Event& iEvent, int i){
  const reco::Track* Track;
  reco::TransientTrack ttTrack;
  TTree* tree;
  if (debug) std::cout << "Getting tree, track, tttrack" << std::endl;
  if(prop_type == 1){ //If want to swith to global, use mu->globalTrack().get()
    tree = CSC_tree;
    if(!(mu->outerTrack().isNonnull())){return;}
    Track = mu->outerTrack().get();
    ttTrack = ttrackBuilder_->build(Track);
  }
  else if (prop_type == 2){
    tree = Tracker_tree;
    if(!(mu->track().isNonnull())){return;}
    Track = mu->track().get();
    ttTrack = ttrackBuilder_->build(Track);
  }
  else if (prop_type == 3){
    tree = Segment_tree;
    if(isCosmic){
      if(!(mu->outerTrack().isNonnull())){return;}
      Track = mu->outerTrack().get();
      ttTrack = ttrackBuilder_->build(Track);
    }
    else{
      if(!(mu->track().isNonnull())){return;}
      Track = mu->track().get();
      ttTrack = ttrackBuilder_->build(Track);
    }
  }
  else{
    std::cout << "Bad prop type, failure." << std::endl; return;
  }
  if(!ttTrack.isValid()){std::cout << "BAD EVENT! NO TRACK" << std::endl;}
  if (debug) std::cout << "Got track, init data" << std::endl;
  data_.init();
  //Muon Info//////////////////////////////////////////////////////
  data_.muon_charge = mu->charge(); data_.muon_pt = mu->pt(); data_.muon_eta = mu->eta(); data_.muon_momentum = mu->momentum().mag2();
  data_.evtNum = iEvent.eventAuxiliary().event(); data_.lumiBlock = iEvent.eventAuxiliary().luminosityBlock(); data_.muonIdx = data_.evtNum*100 + i;
  data_.runNum = iEvent.run();
  //Track Info//////////////////////////////////////////////////////
  data_.track_chi2 = Track->chi2(); data_.track_ndof = Track->ndof();
  if (debug) std::cout << "Counting segments" << std::endl;
  CSCSegmentCounter(mu, data_);
  if(prop_type == 3 and data_.hasME11 != 1){return;}
  //which_track
  //Propagation Info//////////////////////////////////////////////////////
  if (debug) std::cout << "Starting chamber loop" << std::endl;
  for (const auto& ch : GEMGeometry_->etaPartitions()) {
    if (ch->id().station() != 1) continue; //Only takes GE1/1
    GlobalPoint tmp_prop_GP; bool tmp_has_prop = 0;
    propagate_to_GEM(mu, ch, prop_type, tmp_has_prop, tmp_prop_GP, data_);
    if(tmp_has_prop){
      LocalPoint tmp_prop_LP = ch->toLocal(tmp_prop_GP);
      //Rechit Info//////////////////////////////////////////////////////
      GEM_rechit_matcher(ch, tmp_prop_LP, data_);
      if(isMC){
        GEM_simhit_matcher(ch, tmp_prop_GP, data_);
      }
      tree->Fill();
    }
  }
  //tree->Fill();
}



bool analyzer::fidcutCheck(float local_y, float localphi_deg, const GEMEtaPartition* ch){
  const float fidcut_angle = 1.0;
  const float cut_chamber = 5.0;
  const float cut_angle = 5.0 - fidcut_angle;
  auto& parameters(ch->specs()->parameters());
  float height(parameters[2]);
  if ((abs(localphi_deg) < cut_angle) && ((local_y < (height - cut_chamber) && ch->id().roll() == 1) || (local_y > -1.0*(height - cut_chamber) && ch->id().roll() == 8) || (ch->id().roll() != 1 && ch->id().roll() != 8))){return 1;}
  else{return 0;}
}



void analyzer::beginJob(){}
void analyzer::endJob(){
  if(debug){nME11_col_vs_matches->Write();}
  }

DEFINE_FWK_MODULE(analyzer);
