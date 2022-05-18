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



struct CSC_tbma_Data_ChamberLevel{
  void init();
  TTree* book(TTree *t);
  //Only keep required data for 6DOF RdPhi Fitter//////////////////
  float res_x; float res_slope_x;
  float pos_x; float pos_y;
  float angle_x; float angle_y;
  float pz; float pt; float q;
  int detId;
};

void CSC_tbma_Data_ChamberLevel::init(){
  //Only keep required data for 6DOF RdPhi Fitter//////////////////
  res_x = 9999; res_slope_x = 9999;
  pos_x = 9999; pos_y = 9999;
  angle_x = 9999; angle_y = 9999;
  pz = 9999; pt = 9999; q = 9999;
  detId = 9999;
}

TTree* CSC_tbma_Data_ChamberLevel::book(TTree *t){
  edm::Service< TFileService > fs;
  t = fs->make<TTree>("Inner_Prop_ChamberLevel", "Inner_PropChamberLevel");
  //Only keep required data for 6DOF RdPhi Fitter//////////////////
  t->Branch("res_x", &res_x);
  t->Branch("res_slope_x", &res_slope_x);
  t->Branch("pos_x", &pos_x);
  t->Branch("pos_y", &pos_y);
  t->Branch("angle_x", &angle_x);
  t->Branch("angle_y", &angle_y);
  t->Branch("pz", &pz);
  t->Branch("pt", &pt);
  t->Branch("q", &q);
  t->Branch("detId", &detId);
  return t;
}



struct CSC_tbma_Data{
  void init();
  TTree* book(TTree *t);
  //Muon Info//////////////////////////////////////////////////////
  int muon_charge; float muon_pt; float muon_eta; float muon_momentum;
  unsigned long long  evtNum; unsigned long long  lumiBlock; int muonIdx;
  int runNum;
  //Propagation Info//////////////////////////////////////////////////////
  float prop_GP[3]; float prop_LP[3]; float prop_startingPoint_GP[3];
  float prop_localphi_rad; float prop_localphi_deg;
  bool has_prop; bool has_fidcut;
  int prop_location[5];
  //Track Info//////////////////////////////////////////////////////
  float track_chi2; float track_ndof; int which_track;
  //Rechit Info//////////////////////////////////////////////////////
  float rechit_GP[3]; float rechit_LP[3];
  float rechit_localphi_rad; float rechit_localphi_deg;
  bool has_rechit;
  float RdPhi; int rechit_detId;
  int rechit_location[5];
};

void CSC_tbma_Data::init(){
  //Muon Info//////////////////////////////////////////////////////
  muon_charge = 9999; muon_pt = 9999; muon_eta = 9999; muon_momentum = 9999;
  evtNum = 99999999; lumiBlock = 99999999; muonIdx = 99999999; runNum = 99999999;
  //Propagation Info//////////////////////////////////////////////////////
  for(int i=0; i<3; ++i){
    prop_GP[i] = 99999; prop_LP[i] = 99999; prop_startingPoint_GP[i] = 99999;
  }
  prop_localphi_rad = 99999; prop_localphi_deg = 99999;
  has_prop = false; has_fidcut = false;
  for(int i=0; i<5; ++i){
    prop_location[i] = 99999;
  }
  //Track Info//////////////////////////////////////////////////////
  track_chi2 = 999999; track_ndof = 999999; which_track = 999999;
  //Rechit Info//////////////////////////////////////////////////////
  for(int i=0; i<3; ++i){
    rechit_GP[i] = 999999; rechit_LP[i] = 999999;
  }
  rechit_localphi_rad = 999999; rechit_localphi_deg = 999999;
  has_rechit = false;
  RdPhi = 999999; rechit_detId = 999999;
  for(int i=0; i<5; ++i){
    rechit_location[i] = 999999;
  }
}

TTree* CSC_tbma_Data::book(TTree *t){
  edm::Service< TFileService > fs;
  t = fs->make<TTree>("Inner_Prop", "Inner_Prop");
  //Muon Info//////////////////////////////////////////////////////
  t->Branch("muon_charge", &muon_charge); t->Branch("muon_pt", &muon_pt);
  t->Branch("muon_eta", &muon_eta); t->Branch("muon_momentum", &muon_momentum);
  t->Branch("evtNum", &evtNum); t->Branch("lumiBlock", &lumiBlock); t->Branch("muonIdx", &muonIdx);
  t->Branch("runNum", &runNum);
  //Propagation Info//////////////////////////////////////////////////////
  t->Branch("prop_GP", &prop_GP, "prop_GP[3] (x,y,z)/F");
  t->Branch("prop_LP", &prop_LP, "prop_LP[3] (x,y,z)/F");
  t->Branch("prop_startingPoint_GP", &prop_startingPoint_GP, "prop_startingPoint_GP[3] (x,y,z)/F");
  t->Branch("prop_localphi_rad", &prop_localphi_rad);
  t->Branch("prop_localphi_deg", &prop_localphi_deg);
  t->Branch("has_prop", &has_prop);
  t->Branch("has_fidcut", &has_fidcut);
  t->Branch("prop_location", &prop_location, "prop_location[5] (reg, sta, ring, cha, lay)/I");
  //Track Info//////////////////////////////////////////////////////
  t->Branch("track_chi2", &track_chi2); t->Branch("track_ndof", &track_ndof);
  t->Branch("which_track", &which_track);
  //Rechit Info//////////////////////////////////////////////////////
  t->Branch("rechit_GP", &rechit_GP, "rechit_GP[3] (x,y,z)/F");
  t->Branch("rechit_LP", &rechit_LP, "rechit_LP[3] (x,y,z)/F");
  t->Branch("rechit_localphi_rad", &rechit_localphi_rad);
  t->Branch("rechit_localphi_deg", &rechit_localphi_deg);
  t->Branch("has_rechit", &has_rechit);
  t->Branch("RdPhi", &RdPhi);
  t->Branch("rechit_detId", &rechit_detId);
  t->Branch("rechit_location", &rechit_location, "rechit_location[5] (reg, sta, ring, cha, lay)/I");
  return t;
}


class CSC_tbma : public edm::one::EDAnalyzer<> {
public:
  explicit CSC_tbma(const edm::ParameterSet&);
  ~CSC_tbma(){};

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob() ;
  virtual void endJob() ;

  void propagate(const reco::Muon* mu, const edm::Event& iEvent, int i);

  edm::EDGetTokenT<edm::View<reco::Muon> > muons_;
  edm::Handle<View<reco::Muon> > muons;

  edm::EDGetTokenT<CSCSegmentCollection> cscSegments_;
  edm::Handle<CSCSegmentCollection> cscSegments;

  edm::EDGetTokenT<CSCRecHit2DCollection> csc2DRecHits_;
  edm::Handle<CSCRecHit2DCollection> csc2DRecHits;


  edm::Service<TFileService> fs;

  MuonServiceProxy* theService_;
  edm::ESHandle<Propagator> propagator_;
  edm::ESHandle<TransientTrackBuilder> ttrackBuilder_;

  ESHandle<GlobalTrackingGeometry> theTrackingGeometry;

  edm::ESHandle<CSCGeometry> CSCGeometry_;

  bool debug;
  bool isCosmic;

  CSC_tbma_Data data_;
  TTree* Tracker_tree;

  CSC_tbma_Data_ChamberLevel data_ChamberLevel_;
  TTree* Tracker_tree_ChamberLevel;

  bool isMC;

  const edm::ESGetToken<CSCGeometry, MuonGeometryRecord> cscGeomToken_;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttkToken_;
  const edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> geomToken_;
};

CSC_tbma::CSC_tbma(const edm::ParameterSet& iConfig)
  : cscGeomToken_(esConsumes()),
    ttkToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
    geomToken_(esConsumes())
{
  cout << "Begin analyzer" << endl;
  edm::ParameterSet serviceParameters = iConfig.getParameter<edm::ParameterSet>("ServiceParameters");
  theService_ = new MuonServiceProxy(serviceParameters, consumesCollector());

  muons_ = consumes<View<reco::Muon> >(iConfig.getParameter<InputTag>("muons"));
  cscSegments_ = consumes<CSCSegmentCollection>(edm::InputTag("cscSegments"));
  csc2DRecHits_ = consumes<CSCRecHit2DCollection>(iConfig.getParameter<edm::InputTag>("csc2DRecHits"));

  debug = iConfig.getParameter<bool>("debug");
  isCosmic = iConfig.getParameter<bool>("isCosmic");

  Tracker_tree = data_.book(Tracker_tree);
  Tracker_tree_ChamberLevel = data_ChamberLevel_.book(Tracker_tree_ChamberLevel);

}


void CSC_tbma::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  //iSetup.get<MuonGeometryRecord>().get(CSCGeometry_);
  //iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttrackBuilder_);
  //iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);

  CSCGeometry_ = &iSetup.getData(cscGeomToken_);
  ttrackBuilder_ = &iSetup.getData(ttkToken_);
  theTrackingGeometry = &iSetup.getData(geomToken_);

  theService_->update(iSetup);
  auto propagator = theService_->propagator("SteppingHelixPropagatorAny");

  isMC = false;
  if (! iEvent.eventAuxiliary().isRealData()) isMC = true;
  if (! iEvent.getByToken(muons_, muons)) return;
  if (muons->size() == 0) return;


  iEvent.getByToken(cscSegments_, cscSegments);
  iEvent.getByToken(csc2DRecHits_, csc2DRecHits);


  if (debug) cout << "New! EvtNumber = " << iEvent.eventAuxiliary().event() << " LumiBlock = " << iEvent.eventAuxiliary().luminosityBlock() << " RunNumber = " << iEvent.run() << endl;

  for (size_t i = 0; i < muons->size(); ++i){
    edm::RefToBase<reco::Muon> muRef = muons->refAt(i);
    const reco::Muon* mu = muRef.get();
    //if (mu->pt() < 2.0) continue;  //can apply a pt cut later
    if (not mu->standAloneMuon()) continue;
    if (debug) cout << "new standalone" << endl;
    propagate(mu, iEvent, i);
  }
}


void CSC_tbma::propagate(const reco::Muon* mu, const edm::Event& iEvent, int i){
  const reco::Track* Track;
  reco::TransientTrack ttTrack;
  TTree* tree;
  tree = Tracker_tree;
  TTree* tree_ChamberLevel;
  tree_ChamberLevel = Tracker_tree_ChamberLevel;
  if(!(mu->track().isNonnull())){return;}
  Track = mu->track().get();
  ttTrack = ttrackBuilder_->build(Track);
  if(!ttTrack.isValid()){std::cout << "BAD EVENT! NO TRACK" << std::endl;}
  data_.init();
  data_ChamberLevel_.init();
  //Muon Info//////////////////////////////////////////////////////
  data_.muon_charge = mu->charge(); data_.muon_pt = mu->pt(); data_.muon_eta = mu->eta(); data_.muon_momentum = mu->momentum().mag2();
  data_.evtNum = iEvent.eventAuxiliary().event(); data_.lumiBlock = iEvent.eventAuxiliary().luminosityBlock(); data_.muonIdx = data_.evtNum*100 + i;
  data_.runNum = iEvent.run();
  //Track Info//////////////////////////////////////////////////////
  data_.track_chi2 = Track->chi2(); data_.track_ndof = Track->ndof();
  //Track tsos//////////////////////////////////////////////////////
  auto propagator = theService_->propagator("SteppingHelixPropagatorAny");
  TrajectoryStateOnSurface tsos_from_tracker; TrajectoryStateOnSurface tsos_on_chamber;
  //Find outer-most tracker point
  int incoming_or_outgoing = 0;
  std::cout << "TTTrack has " << ttTrack.recHitsSize() << " hits" << std::endl;
  if (ttTrack.outermostMeasurementState().globalPosition().perp() > ttTrack.innermostMeasurementState().globalPosition().perp()){
    tsos_from_tracker = ttTrack.outermostMeasurementState();
    incoming_or_outgoing = 1;
  }
  else{
    tsos_from_tracker = ttTrack.innermostMeasurementState();
    incoming_or_outgoing = -1;
  }
  data_.which_track = incoming_or_outgoing;
  //Get CSC Segments////////////////////////////////////////////////
  std::vector<const CSCSegment*> AllCSCSegments_Matches;
  std::vector<CSCDetId> Chambers_Matches;
  std::vector<CSCDetId> Chambers_Matches_repeats;
  for (auto muon_match : mu->matches()){
    if (muon_match.detector() != 2) continue;
    for (auto muon_seg_match : muon_match.segmentMatches){
      auto CSCSeg = muon_seg_match.cscSegmentRef;
      AllCSCSegments_Matches.push_back(CSCSeg.get());
      if (std::find(Chambers_Matches.begin(), Chambers_Matches.end(), CSCSeg.get()->cscDetId()) != Chambers_Matches.end()){
        Chambers_Matches_repeats.push_back(CSCSeg.get()->cscDetId());
      }
      else{Chambers_Matches.push_back(CSCSeg.get()->cscDetId());}
    }
  }
  std::vector<const CSCSegment*> AllCSCSegments_Collection;
  std::vector<CSCDetId> Chambers_Collection;
  std::vector<CSCDetId> Chambers_Collection_repeats;
  for (CSCSegmentCollection::const_iterator itCSCSeg = cscSegments->begin(); itCSCSeg != cscSegments->end(); itCSCSeg++){
    const CSCSegment* CSCSeg = itCSCSeg->clone();
    AllCSCSegments_Collection.push_back(CSCSeg);
    if (std::find(Chambers_Collection.begin(), Chambers_Collection.end(), CSCSeg->cscDetId()) != Chambers_Collection.end()){
      Chambers_Collection_repeats.push_back(CSCSeg->cscDetId());
    }
    else{Chambers_Collection.push_back(CSCSeg->cscDetId());}
  }
  std::cout << "All Segs Size from Matches    = " << AllCSCSegments_Matches.size() << std::endl;
  std::cout << "All Segs Size from Collection = " << AllCSCSegments_Collection.size() << std::endl;

  //Choose segments and loop over rechits///////////////////////////
  auto AllCSCSegments = AllCSCSegments_Matches;
  auto AllChamberRepeats = Chambers_Matches_repeats;
  for (auto CSCSeg : AllCSCSegments){
    if (debug)std::cout << "Starting loop over CSC Segment, DetID = " << CSCSeg->cscDetId() << std::endl;
    if (debug)std::cout << "Segment rechits = " << CSCSeg->nRecHits() << std::endl;
    if (std::find(AllChamberRepeats.begin(), AllChamberRepeats.end(), CSCSeg->cscDetId()) != AllChamberRepeats.end()){std::cout << "REPEATED CHAMBER" << std::endl; continue;}
    if (CSCSeg->nRecHits() < 5){continue;} //TBMA skips segments with less than 5 rechits
    float sum_weight = 0; float sum_weight_propz = 0; float sum_weight_residual = 0; float sum_weight_propz_propz = 0; float sum_weight_propz_residual = 0; int nRecHits_in_smoothing = 0; // Testing the TBMA method of CHAMBER RESIDUALS including slope
    float sum_weight_propx = 0; float sum_weight_propz_propx = 0; float sum_weight_propy = 0; float sum_weight_propz_propy = 0;
    for (auto CSCRecHit : CSCSeg->recHits()){
      DetId hitId = CSCRecHit->geographicalId();
      CSCDetId cscDetId = CSCDetId(CSCRecHit->geographicalId());
      if (cscDetId.layer() == 0){continue;} //TBMA skips rechits on layer 0 (This makes sense, rechits should always have a layer)
      //Got Hit, Propgate to same layer/////////////////////////////////
      tsos_on_chamber = propagator->propagate(tsos_from_tracker, CSCGeometry_->idToDet(hitId)->surface());
      if (!(tsos_on_chamber.isValid())) continue;
      if (debug)std::cout << "SUCCESS!!! ID " << cscDetId << std::endl;
      //Time to start filling TTree. We have a propagation (tsos_on_chamber) and we have a recHit (CSCRecHit)
      LocalPoint hit_LP = CSCRecHit->localPosition(); GlobalPoint hit_GP = CSCGeometry_->idToDet(hitId)->toGlobal(CSCRecHit->localPosition());
      LocalPoint prop_LP = tsos_on_chamber.localPosition(); GlobalPoint prop_GP = tsos_on_chamber.globalPosition();
      if (debug)std::cout << "prop lp = " << prop_LP << " prop gp = " << prop_GP << std::endl;
      if (debug)std::cout << "hit  lp = " << hit_LP << " hit  gp = " << hit_GP << std::endl;

      data_.prop_GP[0] = prop_GP.x(); data_.prop_GP[1] = prop_GP.y(); data_.prop_GP[2] = prop_GP.z();
      data_.prop_LP[0] = prop_LP.x(); data_.prop_LP[1] = prop_LP.y(); data_.prop_LP[2] = prop_LP.z();
      data_.prop_startingPoint_GP[0] = tsos_from_tracker.globalPosition().x();
      data_.prop_startingPoint_GP[1] = tsos_from_tracker.globalPosition().y();
      data_.prop_startingPoint_GP[2] = tsos_from_tracker.globalPosition().z();
      data_.prop_location[0] = cscDetId.zendcap();
      data_.prop_location[1] = cscDetId.station();
      data_.prop_location[2] = cscDetId.ring();
      data_.prop_location[3] = cscDetId.chamber();
      data_.prop_location[4] = cscDetId.layer();

      data_.rechit_GP[0] = hit_GP.x(); data_.rechit_GP[1] = hit_GP.y(); data_.rechit_GP[2] = hit_GP.z();
      data_.rechit_LP[0] = hit_LP.x(); data_.rechit_LP[1] = hit_LP.y(); data_.rechit_LP[2] = hit_LP.z();
      data_.rechit_location[0] = cscDetId.zendcap();
      data_.rechit_location[1] = cscDetId.station();
      data_.rechit_location[2] = cscDetId.ring();
      data_.rechit_location[3] = cscDetId.chamber();
      data_.rechit_location[4] = cscDetId.layer();

      int strip = CSCGeometry_->layer(cscDetId)->geometry()->nearestStrip(hit_LP);
      double angle = CSCGeometry_->layer(cscDetId)->geometry()->stripAngle(strip) - M_PI/2.;
      double sinAngle = sin(angle);
      double cosAngle = cos(angle);
      double residual = cosAngle * (prop_LP.x() - hit_LP.x()) + sinAngle * (prop_LP.y() - hit_LP.y());  // yes, that's +sin()

      data_.RdPhi = residual;
      //DetId for CSC :    +-   :    _    :    _    :   __   :    _    :
      //              :  Region : Station :  Ring   : Chamber:  Layer  :
      //              :   sign  : *10,000 : *1,000  :  *100  :   *1    :
      data_.rechit_detId = data_.rechit_location[0]*(data_.rechit_location[1]*10000 + data_.rechit_location[2]*1000 + data_.rechit_location[3]*100 + data_.rechit_location[4]);
      if(debug) std::cout << "Success, set variables, example RdPhi = " << residual << std::endl;
      tree->Fill();

      //Adding the TBMA CHAMBER RESIDUALS values
      float xx = CSCRecHit->localPositionError().xx(); float xy = CSCRecHit->localPositionError().xy(); float yy = CSCRecHit->localPositionError().yy();
      std::cout << "Error mat xx = " << xx << " yy = " << yy << " xy = " << xy << std::endl;
      float weight = 1.0/(xx*cosAngle*cosAngle + 2.0*xy*sinAngle*cosAngle + yy*sinAngle*sinAngle);
      float prop_local_z = CSCGeometry_->idToDet(CSCSeg->cscDetId())->toLocal((tsos_on_chamber.globalPosition())).z();
      sum_weight += weight;
      //TBMA residual
      sum_weight_propz += weight*prop_local_z;
      sum_weight_residual += weight*residual;
      sum_weight_propz_propz += weight*prop_local_z*prop_local_z;
      sum_weight_propz_residual += weight*prop_local_z*residual;
      nRecHits_in_smoothing += 1;
      std::cout << "RecHit " << nRecHits_in_smoothing << "Using w/propz/res " << weight << "/" << prop_local_z << "/" << residual << std::endl;

      //TBMA prop x
      float prop_x = tsos_on_chamber.localPosition().x();
      sum_weight_propx += weight*prop_x;
      sum_weight_propz_propx += weight*prop_local_z*prop_x;

      //TBMA prop y
      float prop_y = tsos_on_chamber.localPosition().y();
      sum_weight_propy += weight*prop_y;
      sum_weight_propz_propy += weight*prop_local_z*prop_y;

    }
    std::cout << "Begin smoothing! On Chamber " << CSCSeg->cscDetId() << std::endl;
    float delta_res = (sum_weight * sum_weight_propz_propz) - (sum_weight_propz * sum_weight_propz);
    float real_residual = ((sum_weight_propz_propz * sum_weight_residual) - (sum_weight_propz * sum_weight_propz_residual)) / delta_res;
    float real_slope = ((sum_weight * sum_weight_propz_residual) - (sum_weight_propz * sum_weight_residual)) / delta_res;
    std::cout << "Smoothed delta/res/slope " << delta_res << "/" << real_residual << "/" << real_slope << std::endl;

    float delta_propx = (sum_weight * sum_weight_propz_propz) - (sum_weight_propz * sum_weight_propz);
    float real_x = ((sum_weight_propz_propz * sum_weight_propx) - (sum_weight_propz * sum_weight_propz_propx)) / delta_propx;
    float real_angle_x = ((sum_weight * sum_weight_propz_propx) - (sum_weight_propz * sum_weight_propx)) / delta_propx;

    float delta_propy = (sum_weight * sum_weight_propz_propz) - (sum_weight_propz * sum_weight_propz);
    float real_y = ((sum_weight_propz_propz * sum_weight_propy) - (sum_weight_propz * sum_weight_propz_propy)) / delta_propy;
    float real_angle_y = ((sum_weight * sum_weight_propz_propy) - (sum_weight_propz * sum_weight_propy)) / delta_propy;


    data_ChamberLevel_.res_x = real_residual;
    data_ChamberLevel_.res_slope_x = real_slope;
    data_ChamberLevel_.pos_x = real_x;
    data_ChamberLevel_.pos_y = real_y;
    data_ChamberLevel_.angle_x = real_angle_x;
    data_ChamberLevel_.angle_y = real_angle_y;
    data_ChamberLevel_.pz = mu->pz();
    data_ChamberLevel_.pt = mu->pt();
    data_ChamberLevel_.q = mu->charge();
    //DetId for CSC :    +-   :    _    :    _    :   __   :
    //Chamber Level :  Region : Station :  Ring   : Chamber:
    //              :   sign  :  *1,000 :   *100  :    *1  :
    data_ChamberLevel_.detId = CSCSeg->cscDetId().zendcap()*(CSCSeg->cscDetId().station()*1000 + CSCSeg->cscDetId().ring()*100 + CSCSeg->cscDetId().chamber());
    std::cout << "On a chamber level!" << std::endl;
    std::cout << "RdPhi    " << data_ChamberLevel_.res_x << std::endl;
    std::cout << "RdPhiDz  " << data_ChamberLevel_.res_slope_x << std::endl;
    std::cout << "pos x    " << data_ChamberLevel_.pos_x << std::endl;
    std::cout << "pos y    " << data_ChamberLevel_.pos_y << std::endl;
    std::cout << "dxdz     " << data_ChamberLevel_.angle_x << std::endl;
    std::cout << "dydz     " << data_ChamberLevel_.angle_y << std::endl;
    std::cout << "pZ       " << data_ChamberLevel_.pz << std::endl;
    std::cout << "pT       " << data_ChamberLevel_.pt << std::endl;
    std::cout << "Charge   " << data_ChamberLevel_.q << std::endl;
    std::cout << "DetID    " << data_ChamberLevel_.detId << std::endl;


    tree_ChamberLevel->Fill();
  }


  /*
  std::cout << "Muon matches length = " << mu->numberOfMatches() << std::endl;
  std::cout << "CSC Segments collection = " << cscSegments->size() << std::endl;


  //std::cout << typeid(muon_segment_matches).name() << std::endl;
  if (isCosmic){
    std::cout << "Cosmic Data" << std::endl;
    for (auto muon_match : mu->matches()){
      if (muon_match.detector() != 2) continue;
      for (auto muon_seg_match : muon_match.segmentMatches){
        auto CSCSeg = muon_seg_match.cscSegmentRef;
        std::cout << "Starting loop over CSC Segment, DetID = " << CSCSeg->cscDetId() << std::endl;
        std::cout << "Segment rechits = " << CSCSeg->nRecHits() << std::endl;
        if (CSCSeg->nRecHits() < 5){continue;} //TBMA skips segments with less than 5 rechits
        for (auto CSCRecHit : CSCSeg->recHits()){
          DetId hitId = CSCRecHit->geographicalId();
          CSCDetId cscDetId = CSCDetId(CSCRecHit->geographicalId());
          if (cscDetId.layer() == 0){continue;} //TBMA skips rechits on layer 0 (This makes sense, rechits should always have a layer)
          //Got Hit, Propgate to same layer/////////////////////////////////
          tsos_on_chamber = propagator->propagate(tsos_from_tracker, CSCGeometry_->idToDet(hitId)->surface());
          if (tsos_on_chamber.isValid()){
            std::cout << "SUCCESS!!! ID " << cscDetId << std::endl;
            //Time to start filling TTree. We have a propagation (tsos_on_chamber) and we have a recHit (CSCRecHit)
            LocalPoint hit_LP = CSCRecHit->localPosition(); GlobalPoint hit_GP = CSCGeometry_->idToDet(hitId)->toGlobal(CSCRecHit->localPosition());
            LocalPoint prop_LP = tsos_on_chamber.localPosition(); GlobalPoint prop_GP = tsos_on_chamber.globalPosition();
            std::cout << "prop lp = " << prop_LP << " prop gp = " << prop_GP << std::endl;
            std::cout << "hit  lp = " << hit_LP << " hit  gp = " << hit_GP << std::endl;
          }
        }
      }
    }
  }
  //else{
  if (isCosmic){ //Testing both cases
    std::cout << "PP Data" << std::endl;
    for (CSCSegmentCollection::const_iterator CSCSeg = cscSegments->begin(); CSCSeg != cscSegments->end(); CSCSeg++){
      std::cout << "Starting loop over CSC Segment, DetID = " << CSCSeg->cscDetId() << std::endl;
      std::cout << "Segment rechits = " << CSCSeg->nRecHits() << std::endl;
      if (CSCSeg->nRecHits() < 5){continue;} //TBMA skips segments with less than 5 rechits
      for (auto CSCRecHit : CSCSeg->recHits()){
        DetId hitId = CSCRecHit->geographicalId();
        CSCDetId cscDetId = CSCDetId(CSCRecHit->geographicalId());
        if (cscDetId.layer() == 0){continue;} //TBMA skips rechits on layer 0 (This makes sense, rechits should always have a layer)
        //Got Hit, Propgate to same layer/////////////////////////////////
        tsos_on_chamber = propagator->propagate(tsos_from_tracker, CSCGeometry_->idToDet(hitId)->surface());
        if (tsos_on_chamber.isValid()){
          std::cout << "SUCCESS!!! ID " << cscDetId << " and position was " << tsos_on_chamber.globalPosition() << std::endl;
          //Time to start filling TTree. We have a propagation (tsos_on_chamber) and we have a recHit (CSCRecHit)
          LocalPoint hit_LP = CSCRecHit->localPosition(); GlobalPoint hit_GP = CSCGeometry_->idToDet(hitId)->toGlobal(CSCRecHit->localPosition());
          std::cout << "hit lp = " << hit_LP << " hit gp = " << hit_GP << std::endl;
        }
      }
    }
  }
  */
  //tree->Fill();
}


void CSC_tbma::beginJob(){}
void CSC_tbma::endJob(){}

DEFINE_FWK_MODULE(CSC_tbma);
