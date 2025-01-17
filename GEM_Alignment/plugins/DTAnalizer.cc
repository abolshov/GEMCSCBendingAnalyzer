#include <assert.h>
#include <memory>
#include <cmath>
#include <iostream>
#include <sstream>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

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

#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment2D.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment2DCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4D.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecHit1DPair.h"
#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/MuonDetId/interface/DTBtiId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"

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

struct DT_tbma_data_sectorLevel {
    void init();
    TTree* book(TTree* t);
    float res_dx;
    float res_dy;
    float res_dxdz;
    float res_dydz;
    float prop_local_x;
    float prop_local_y;

    // next two are needed for the fitter
    float sectorLevel_LP[3];
    float sectorLevel_GP[3];

    float global_z;
    float prop_phi;
    int location[3];
    int charge;
    float pz;
    float pt;
    int detId;
};

void DT_tbma_data_sectorLevel::init() {
    res_dx = 9999.0;
    res_dy = 9999.0;
    res_dxdz = 9999.0;
    res_dydz = 9999.0;
    prop_local_x = 9999.0;
    prop_local_y = 9999.0;
    for (size_t i = 0; i < 3; i++) {
        sectorLevel_LP[i] = 9999.0;
        sectorLevel_GP[i] = 9999.0;
    }
    prop_phi = 9999.0;
    for (size_t i = 0; i < 3; i++) {
        location[i] = 9999;
    }
    global_z = 9999.0;
    charge = 9999;
    pz = 9999.0;
    pt = 9999.0;
    detId = 9999;
}

TTree* DT_tbma_data_sectorLevel::book(TTree* t) {
    edm::Service< TFileService > fs;
    t = fs->make<TTree>("Inner_Prop_ChamberLevel", "Inner_PropChamberLevel");
    t->Branch("res_dx", &res_dx);
    t->Branch("res_dy", &res_dy);
    t->Branch("res_dxdz", &res_dxdz);
    t->Branch("res_dydz", &res_dydz);
    t->Branch("prop_local_x", &prop_local_x);
    t->Branch("prop_local_y", &prop_local_y);
    t->Branch("prop_phi", &prop_phi);
    t->Branch("location", &location, "location[3] (wheel, station, sector)/I");
    t->Branch("global_z", &global_z);
    t->Branch("sectorLevel_LP", &sectorLevel_LP, "sectorLevel_LP[3] (x,y,z)/F");
    t->Branch("sectorLevel_GP", &sectorLevel_GP, "sectorLevel_GP[3] (x,y,z)/F");
    t->Branch("charge", &charge);
    t->Branch("pz", &pz);
    t->Branch("pt", &pt);
    t->Branch("detId", &detId);
    return t;
}

struct DT_tbma_data {
    void init();
    TTree* book(TTree *t);
    //Muon Info//////////////////////////////////////////////////////
    int muon_charge; float muon_pt; float muon_eta; float muon_momentum;
    unsigned long long  evtNum; unsigned long long  lumiBlock; int muonIdx;
    int runNum;
    //Propagation Info//////////////////////////////////////////////////////
    float prop_GP[3]; float prop_LP[3]; float prop_startingPoint_GP[3];
    int prop_location[4]; float prop_global_phi;
    //Track Info//////////////////////////////////////////////////////
    float track_chi2; float track_ndof; int which_track;
    //Rechit Info//////////////////////////////////////////////////////
    float rechit_GP[3]; float rechit_LP[3];
    bool has_rechit;
    float dx; float dy; float dxdz; float dydz; int rechit_detId;
    int rechit_location[4];
};

void DT_tbma_data::init() {
    //Muon Info//////////////////////////////////////////////////////
    muon_charge = 9999; muon_pt = 9999; muon_eta = 9999; muon_momentum = 9999;
    evtNum = 99999999; lumiBlock = 99999999; muonIdx = 99999999; runNum = 99999999;
    //Propagation Info//////////////////////////////////////////////////////
    for(int i=0; i<3; ++i){
        prop_GP[i] = 99999; prop_LP[i] = 99999; prop_startingPoint_GP[i] = 99999;
    }
    for(int i=0; i<4; ++i){
        prop_location[i] = 99999;
    }
    prop_global_phi = 99999;
    //Track Info//////////////////////////////////////////////////////
    track_chi2 = 999999; track_ndof = 999999; which_track = 999999;
    //Rechit Info//////////////////////////////////////////////////////
    for(int i=0; i<3; ++i){
        rechit_GP[i] = 999999; rechit_LP[i] = 999999;
    }
    has_rechit = false;
    dx = dy = dxdz = dydz = 99999;
    for(int i=0; i<4; ++i){
        rechit_location[i] = 999999;
    }
    rechit_detId = 99999;
}

TTree* DT_tbma_data::book(TTree *t){
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
    t->Branch("prop_location", &prop_location, "prop_location[4] (wheel, station, sector, supLay)/I");
    t->Branch("prop_global_phi", &prop_global_phi);
    //Track Info//////////////////////////////////////////////////////
    t->Branch("track_chi2", &track_chi2); t->Branch("track_ndof", &track_ndof);
    t->Branch("which_track", &which_track);
    //Rechit Info//////////////////////////////////////////////////////
    t->Branch("rechit_GP", &rechit_GP, "rechit_GP[3] (x,y,z)/F");
    t->Branch("rechit_LP", &rechit_LP, "rechit_LP[3] (x,y,z)/F");
    t->Branch("dx", &dx);
    t->Branch("dy", &dy);
    t->Branch("dxdz", &dxdz);
    t->Branch("dydz", &dydz);
    t->Branch("has_rechit", &has_rechit);
    t->Branch("rechit_detId", &rechit_detId);
    t->Branch("rechit_location", &rechit_location, "rechit_location[4] (wheel, station, sector, supLay)/I");
  return t;
}

class DT_tbma : public edm::EDAnalyzer {
public:
  explicit DT_tbma(const edm::ParameterSet&);
  ~DT_tbma(){};

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob() ;
  virtual void endJob() ;

  void propagate(const reco::Muon* mu, const edm::Event& iEvent, int i);
  //void propagate_to_CSC(const reco::Muon* mu, const CSCLayer* ch, bool &tmp_has_prop, GlobalPoint &pos_GP, CSC_tbma_Data& data_, reco::TransientTrack track);
  //void CSC_rechit_matcher(const CSCLayer* ch, LocalPoint prop_LP, CSC_tbma_Data& data_);
  //bool fidcutCheck(float local_y, float localphi_deg);

  edm::EDGetTokenT<edm::View<reco::Muon> > muons_;

  edm::EDGetTokenT<DTRecSegment4DCollection> dt4DSegments_;
  edm::Handle<DTRecSegment4DCollection> dt4DSegments;

  edm::EDGetTokenT<DTRecHitCollection> dt1DRecHits_;
  edm::Handle<DTRecHitCollection> dt1DRecHits;


  edm::Service<TFileService> fs;

  MuonServiceProxy* theService_;
  edm::ESHandle<Propagator> propagator_;
  edm::ESHandle<TransientTrackBuilder> ttrackBuilder_;

  ESHandle<GlobalTrackingGeometry> theTrackingGeometry;

  edm::ESHandle<DTGeometry> DTGeometry_;

  bool debug;
  bool isCosmic;

  DT_tbma_data data_;
  TTree* Tracker_tree;

  DT_tbma_data_sectorLevel data_sectorLevel_;
  TTree* Tracker_tree_sectorLevel;

  bool isMC;

  // new stuff required for CMSSW update
  const edm::ESGetToken<DTGeometry, MuonGeometryRecord> dtGeomToken_;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttkToken_;
  const edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> geomToken_;
};

DT_tbma::DT_tbma(const edm::ParameterSet& iConfig)
    // for CMSSW update
    : dtGeomToken_(esConsumes()),
    ttkToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
    geomToken_(esConsumes())
{
    cout << "Begin analyzer" << endl;
    edm::ParameterSet serviceParameters = iConfig.getParameter<edm::ParameterSet>("ServiceParameters");
    theService_ = new MuonServiceProxy(serviceParameters, consumesCollector());

    muons_ = consumes<View<reco::Muon> >(iConfig.getParameter<InputTag>("muons"));
    dt4DSegments_ = consumes<DTRecSegment4DCollection>(edm::InputTag("dt4DSegments"));
    dt1DRecHits_ = consumes<DTRecHitCollection>(iConfig.getParameter<edm::InputTag>("dt1DRecHits"));

    debug = iConfig.getParameter<bool>("debug");
    isCosmic = iConfig.getParameter<bool>("isCosmic");


    Tracker_tree = data_.book(Tracker_tree);
    Tracker_tree_sectorLevel = data_sectorLevel_.book(Tracker_tree_sectorLevel);

}

void DT_tbma::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  // iSetup.get<MuonGeometryRecord>().get(DTGeometry_);
  // iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttrackBuilder_);
  // iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);

  // for CMSSW update
  DTGeometry_ = &iSetup.getData(dtGeomToken_);
  ttrackBuilder_ = &iSetup.getData(ttkToken_);
  theTrackingGeometry = &iSetup.getData(geomToken_);

  theService_->update(iSetup);
  auto propagator = theService_->propagator("SteppingHelixPropagatorAny");

  isMC = false;
  if (! iEvent.eventAuxiliary().isRealData()) isMC = true;
  edm::Handle<View<reco::Muon> > muons;
  if (! iEvent.getByToken(muons_, muons)) return;
  if (muons->size() == 0) return;


  iEvent.getByToken(dt4DSegments_, dt4DSegments);
  iEvent.getByToken(dt1DRecHits_, dt1DRecHits);


  // if (debug) cout << "New! EvtNumber = " << iEvent.eventAuxiliary().event() << " LumiBlock = "
  //                 << iEvent.eventAuxiliary().luminosityBlock() << " RunNumber = " << iEvent.run() << endl;

  for (size_t i = 0; i < muons->size(); ++i){
    edm::RefToBase<reco::Muon> muRef = muons->refAt(i);
    const reco::Muon* mu = muRef.get();
    //if (mu->pt() < 2.0) continue;  //can apply a pt cut later
    if (not mu->standAloneMuon()) continue;
    // if (debug) cout << "new standalone" << endl;
    propagate(mu, iEvent, i);
  }
}

void DT_tbma::propagate(const reco::Muon* mu, const edm::Event& iEvent, int i){
  const reco::Track* Track;
  reco::TransientTrack ttTrack;
  TTree* tree;
  tree = Tracker_tree;

  // sector level tree for slope calculation
  TTree* tree_sectorLevel;
  tree_sectorLevel = Tracker_tree_sectorLevel;

  if(!(mu->track().isNonnull())){return;}
  Track = mu->track().get();
  ttTrack = ttrackBuilder_->build(Track);
  // if(!ttTrack.isValid()){std::cout << "BAD EVENT! NO TRACK" << std::endl;}
  data_.init();
  //Muon Info//////////////////////////////////////////////////////
  data_.muon_charge = mu->charge();
  data_.muon_pt = mu->pt();
  data_.muon_eta = mu->eta();
  data_.muon_momentum = mu->momentum().mag2();
  data_.evtNum = iEvent.eventAuxiliary().event();
  data_.lumiBlock = iEvent.eventAuxiliary().luminosityBlock();
  data_.muonIdx = data_.evtNum*100 + i;
  data_.runNum = iEvent.run();
  //Track Info//////////////////////////////////////////////////////
  data_.track_chi2 = Track->chi2(); data_.track_ndof = Track->ndof();
  //Track tsos//////////////////////////////////////////////////////
  auto propagator = theService_->propagator("SteppingHelixPropagatorAny");
  TrajectoryStateOnSurface tsos_from_tracker;
  TrajectoryStateOnSurface tsos_on_chamber;
  //Find outer-most tracker point
  int incoming_or_outgoing = 0;
  // std::cout << "TTTrack has " << ttTrack.recHitsSize() << " hits" << std::endl;
  if (ttTrack.outermostMeasurementState().globalPosition().perp() > ttTrack.innermostMeasurementState().globalPosition().perp()){
    tsos_from_tracker = ttTrack.outermostMeasurementState();
    incoming_or_outgoing = 1;
  }
  else{
    tsos_from_tracker = ttTrack.innermostMeasurementState();
    incoming_or_outgoing = -1;
  }
  data_.which_track = incoming_or_outgoing;
  // std::cout << "Global position is " << tsos_from_tracker.globalPosition() << '\n';
  //tree->Fill();

  //Get DT Segments////////////////////////////////////////////////
  std::vector<const DTRecSegment4D*> AllDTSegments_Matches;

  // loop over muon objects, extract rechits in DTs matching track of this muon
  for (auto muon_match : mu->matches()){
      if (muon_match.detector() != 1) continue;
      for (auto muon_seg_match: muon_match.segmentMatches) {
          auto DTSeg = muon_seg_match.dtSegmentRef;
          AllDTSegments_Matches.push_back(DTSeg->clone());
      }
  }
  // std::cout << "Size of DT segment is " << AllDTSegments_Matches.size() << std::endl;

// 4D DT segment = 2 2D segments (x, slope), (y, slope)
// choose 4D segment and loop over recHits
// looping over 4D segments
    for (auto dtSeg4D: AllDTSegments_Matches) {
        DTChamberId sectorId = dtSeg4D->chamberId();
        // cout << "chamber ID from 4D segment " << sectorId.wheel() << "/" << sectorId.station() << "/" << sectorId.sector() << std::endl;

        // initialize sums for accumulation
        float sum_weight_13 = 0.0; // for x layers
        float sum_weight_2 = 0.0; // for y layers
        float sum_weight_propz_13 = 0.0;
        float sum_weight_propz_2 = 0.0;
        float sum_weight_residual_dx = 0.0;
        float sum_weight_residual_dy = 0.0;
        float sum_weight_propz_propz_13 = 0.0;
        float sum_weight_propz_propz_2 = 0.0;
        float sum_weight_propz_residual_dx = 0;
        float sum_weight_propz_residual_dy = 0;
        float sum_weight_propx_13 = 0.0;
        float sum_weight_propx_2 = 0.0;
        float sum_weight_propz_propx_13 = 0.0;
        float sum_weight_propz_propx_2 = 0.0;
        float sum_weight_propy_13 = 0.0;
        float sum_weight_propy_2 = 0.0;
        float sum_weight_propz_propy_13 = 0.0;
        float sum_weight_propz_propy_2 = 0.0;

        // get 2D segment from 4D segment
        auto vDTSeg2D = dtSeg4D->recHits();
        // loop over contets of the 4D segment - 2D segments
        for(auto itDTSeg2D : vDTSeg2D){
    		  // get 1D segment from 2D segment
              auto vDTSeg1D = itDTSeg2D->recHits();

              // looping over 1d hits; propagate them and calulate residuals
              for (auto itDTHits1D: vDTSeg1D) {

                DetId hitId  = itDTHits1D->geographicalId(); // get global position of the hit
                const DTSuperLayerId superLayerId(hitId.rawId());
    		    const DTLayerId layerId(hitId.rawId());
                DTChamberId chamberId = DTChamberId(hitId);

                //Extrapolate the propagated position for residual
    			TrajectoryStateOnSurface tsos_on_chamber;
    			tsos_on_chamber = propagator->propagate( tsos_from_tracker, DTGeometry_->idToDet(hitId)->surface() );

                // filling TTree
                // LocalPoint hit_LP = itDTHits1D->localPosition();
                GlobalPoint hit_GP = DTGeometry_->idToDet(hitId)->toGlobal(itDTHits1D->localPosition());
                // LocalPoint prop_LP = tsos_on_chamber.localPosition();
                GlobalPoint prop_GP = tsos_on_chamber.globalPosition();

                LocalPoint hit_LP = DTGeometry_->idToDet(sectorId)->toLocal(DTGeometry_->idToDet(hitId)->toGlobal(itDTHits1D->localPosition()));
                LocalPoint prop_LP = DTGeometry_->idToDet(sectorId)->toLocal(DTGeometry_->idToDet(hitId)->toGlobal(tsos_on_chamber.localPosition()));

                data_.prop_GP[0] = prop_GP.x();
                data_.prop_GP[1] = prop_GP.y();
                data_.prop_GP[2] = prop_GP.z();
                data_.prop_LP[0] = prop_LP.x();
                data_.prop_LP[1] = prop_LP.y();
                data_.prop_LP[2] = prop_LP.z();
                data_.prop_global_phi = prop_GP.phi();
                data_.prop_startingPoint_GP[0] = tsos_from_tracker.globalPosition().x();
                data_.prop_startingPoint_GP[1] = tsos_from_tracker.globalPosition().y();
                data_.prop_startingPoint_GP[2] = tsos_from_tracker.globalPosition().z();
                data_.prop_location[0] = chamberId.wheel();
                data_.prop_location[1] = chamberId.station();
                data_.prop_location[2] = chamberId.sector();
                // data_.prop_location[3] = dtDetId.superlayer();

                // cout << "superlayer ID " << superLayerId.superLayer() << " " << layerId.layer() << std::endl;

                data_.rechit_GP[0] = hit_GP.x();
                data_.rechit_GP[1] = hit_GP.y();
                data_.rechit_GP[2] = hit_GP.z();
                data_.rechit_LP[0] = hit_LP.x();
                data_.rechit_LP[1] = hit_LP.y();
                data_.rechit_LP[2] = hit_LP.z();
                data_.rechit_location[0] = chamberId.wheel();
                data_.rechit_location[1] = chamberId.station();
                data_.rechit_location[2] = chamberId.sector();
                // data_.rechit_location[3] = dtDetId.superlayer

                if ( superLayerId.superlayer() == 2  && vDTSeg1D.size() >= 3 ) {
                    float residual_dy = prop_LP.y() - hit_LP.y();
                    data_.dx = 99999;
                    data_.dy = residual_dy;
                    float xx = itDTHits1D->localPositionError().xx();
                    float weight_2 = 1/xx;
                    sum_weight_2 += weight_2;

                    // float prop_local_z = DTGeometry_->idToDet(hitId)->toLocal((tsos_on_chamber.globalPosition())).z();
                    float prop_local_z = prop_LP.z();
                    sum_weight_propz_2 += weight_2*prop_local_z;
                    sum_weight_propz_propz_2 += weight_2*prop_local_z*prop_local_z;
                    sum_weight_propz_residual_dy += weight_2*prop_local_z*residual_dy;
                    sum_weight_residual_dy += weight_2*residual_dy;

                    // float prop_x = tsos_on_chamber.localPosition().x();
                    float prop_x = prop_LP.x();
                    sum_weight_propx_2 += weight_2*prop_x;
                    sum_weight_propz_propx_2 += weight_2*prop_local_z*prop_x;

                    // float prop_y = tsos_on_chamber.localPosition().y();
                    float prop_y = prop_LP.y();
                    sum_weight_propy_2 += weight_2*prop_y;
                    sum_weight_propz_propy_2 += weight_2*prop_local_z*prop_y;

                    tree->Fill();
                }

                if ((superLayerId.superlayer() == 1 || superLayerId.superlayer() == 3) && vDTSeg1D.size() >= 6) {
                    float residual_dx = prop_LP.x() - hit_LP.x();
                    data_.dx = residual_dx;
                    data_.dy = 99999;
                    float xx = itDTHits1D->localPositionError().xx();
                    float weight_13 = 1/xx;
                    sum_weight_13 += weight_13;

                    // float prop_local_z = DTGeometry_->idToDet(hitId)->toLocal((tsos_on_chamber.globalPosition())).z();
                    float prop_local_z = prop_LP.z();
                    sum_weight_propz_13 += weight_13*prop_local_z;
                    sum_weight_propz_propz_13 += weight_13*prop_local_z*prop_local_z;
                    sum_weight_propz_residual_dx += weight_13*prop_local_z*residual_dx;
                    sum_weight_residual_dx += weight_13*residual_dx;

                    // float prop_x = tsos_on_chamber.localPosition().x();
                    float prop_x = prop_LP.x();
                    sum_weight_propx_13 += weight_13*prop_x;
                    sum_weight_propz_propx_13 += weight_13*prop_local_z*prop_x;

                    // float prop_y = tsos_on_chamber.localPosition().y();
                    float prop_y = prop_LP.y();
                    sum_weight_propy_13 += weight_13*prop_y;
                    sum_weight_propz_propy_13 += weight_13*prop_local_z*prop_y;

                    if (chamberId.wheel() == 2 && chamberId.station() == 4 && chamberId.sector() == 9) {
                        cout << "residual dx is " << residual_dx << std::endl;
                        cout << "prop_local_z = " << prop_local_z << std::endl;
                        cout << "prop local position is (" << prop_LP.x() << ", " << prop_LP.y() << ")" << std::endl;
                        cout << "hit local position is (" << hit_LP.x() << ", " << hit_LP.y() << ")" << std::endl;
                        cout << "weight is " << weight_13 << std::endl;
                    }

                    tree->Fill();
                }
            }
        }
        // calculate smoothened quantities here
        float delta_2 = (sum_weight_2 * sum_weight_propz_propz_2) - (sum_weight_propz_2 * sum_weight_propz_2);
        float delta_13 = (sum_weight_13 * sum_weight_propz_propz_13) - (sum_weight_propz_13 * sum_weight_propz_13);
        float real_residual_dx = ((sum_weight_propz_propz_13 * sum_weight_residual_dx) - (sum_weight_propz_13 * sum_weight_propz_residual_dx)) / delta_13;
        float real_residual_dy = ((sum_weight_propz_propz_2 * sum_weight_residual_dy) - (sum_weight_propz_2 * sum_weight_propz_residual_dy)) / delta_2;
        float real_slope_dxdz = ((sum_weight_13 * sum_weight_propz_residual_dx) - (sum_weight_propz_13 * sum_weight_residual_dx)) / delta_13;
        float real_slope_dydz = ((sum_weight_2 * sum_weight_propz_residual_dy) - (sum_weight_propz_2 * sum_weight_residual_dy)) / delta_2;

        float real_x = ((sum_weight_propz_propz_13 * sum_weight_propx_13) - (sum_weight_propz_13 * sum_weight_propz_propx_13)) / delta_13;
        float real_y = 0.0;
        if (sectorId.station() == 4) {
            real_y = ((sum_weight_propz_propz_13 * sum_weight_propy_13) - (sum_weight_propz_13 * sum_weight_propz_propy_13)) / delta_13;
        }
        else {
            real_y = ((sum_weight_propz_propz_2 * sum_weight_propy_2) - (sum_weight_propz_2 * sum_weight_propz_propy_2)) / delta_2;
        }

        LocalPoint sectorLevelPosLocal(real_x, real_y, 0.0);
        GlobalPoint sectorLevelPosGlobal = DTGeometry_->idToDet(sectorId)->toGlobal(sectorLevelPosLocal);

        if (sectorId.wheel() == 2 && sectorId.station() == 4 && sectorId.sector() == 9) {
                cout << "chamber " << sectorId.wheel() << "/" << sectorId.station() << "/" << sectorId.sector() << std::endl;
                cout << "sum_weight_13 = " << sum_weight_13 << std::endl;
                cout << "sum_weight_propz_propz_13 = " << sum_weight_propz_propz_13 << std::endl;
                cout << "sum_weight_propz_13 = " << sum_weight_propz_13 << std::endl;
                cout << "sum_weight_residual_dx = " << sum_weight_residual_dx << std::endl;
                cout << "sum_weight_propz_residual_dx = " << sum_weight_propz_residual_dx << std::endl;
                cout << "smoothened variables:" << std::endl;
                cout << "real residual dx is " << real_residual_dx << std::endl;
                cout << "(real_x, real_y, real_phi, real_global_z) = " << "(" << real_x << ", " << real_y << ", " << sectorLevelPosGlobal.phi() << ", " << sectorLevelPosGlobal.z() << ")" << std::endl;
                cout << "*****************************************" << std::endl;
            }


        // filling the sector level tree
        data_sectorLevel_.res_dx = real_residual_dx;
        if (sectorId.station() != 4) {
            data_sectorLevel_.res_dy = real_residual_dy;
        }
        else {
            data_sectorLevel_.res_dy = 9999.0;
        }
        data_sectorLevel_.res_dxdz = real_slope_dxdz;
        data_sectorLevel_.res_dydz = real_slope_dydz;
        data_sectorLevel_.prop_local_x = real_x;
        data_sectorLevel_.prop_local_y = real_y;
        data_sectorLevel_.sectorLevel_LP[0] = real_x;
        data_sectorLevel_.sectorLevel_LP[1] = real_y;
        data_sectorLevel_.sectorLevel_LP[2] = 0.0;
        data_sectorLevel_.sectorLevel_GP[0] = sectorLevelPosGlobal.x();
        data_sectorLevel_.sectorLevel_GP[1] = sectorLevelPosGlobal.y();
        data_sectorLevel_.sectorLevel_GP[2] = sectorLevelPosGlobal.z();
        data_sectorLevel_.prop_phi = sectorLevelPosGlobal.phi();
        data_sectorLevel_.global_z = sectorLevelPosGlobal.z();
        data_sectorLevel_.pz = mu->pz();
        data_sectorLevel_.pt = mu->pt();
        data_sectorLevel_.charge = mu->charge();
        data_sectorLevel_.location[0] = sectorId.wheel();
        data_sectorLevel_.location[1] = sectorId.station();
        data_sectorLevel_.location[2] = sectorId.sector();
        tree_sectorLevel->Fill();
    }

}







void DT_tbma::beginJob(){}
void DT_tbma::endJob(){}

DEFINE_FWK_MODULE(DT_tbma);
