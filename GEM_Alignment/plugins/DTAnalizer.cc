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

  bool isMC;
};

DT_tbma::DT_tbma(const edm::ParameterSet& iConfig)
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

}

void DT_tbma::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  iSetup.get<MuonGeometryRecord>().get(DTGeometry_);

  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttrackBuilder_);

  iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);

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

// for (auto dtSeg: AllDTSegments_Matches) {
//     std::cout << "rechit size is " << (dtSeg->recHits()).size() << std::endl;
//     for (auto recHit: dtSeg->recHits()) {
//         std::cout << "rechit local position is " << recHit->localPosition() << std::endl;
//         for (auto recHit2 : recHit->recHits()){
//             std::cout << "Yup i'm right" << std::endl;
//         }
//     }
// }

// 4D DT segment = 2 2D segments (x, slope), (y, slope)
// choose segment and loop over recHits
for (auto dtSeg: AllDTSegments_Matches) {
    //std::vector<const TrackingRecHit*> vDTSeg2D = (*dtSeg).recHits(); // get 2D segment from 4D segment
    auto vDTSeg2D = dtSeg->recHits();

    //for ( std::vector<TrackingRecHit*>::const_iterator itDTSeg2D =  vDTSeg2D.begin(); itDTSeg2D != vDTSeg2D.end(); ++itDTSeg2D ) {
    for(auto itDTSeg2D : vDTSeg2D){
		  // std::vector<const TrackingRecHit*> vDTHits1D =  itDTSeg2D->recHits(); // get 1D segment from 2D segment
          auto vDTSeg1D = itDTSeg2D->recHits();

          // float sum_weight = 0;
          // float sum_weight_propz = 0;
          // float sum_weight_residual = 0;
          // float sum_weight_propz_propz = 0;
          // float sum_weight_propz_residual = 0;
          // int nRecHits_in_smoothing = 0; // Testing the TBMA method of CHAMBER RESIDUALS including slope
          // float sum_weight_propx = 0;
          // float sum_weight_propz_propx = 0;
          // float sum_weight_propy = 0;
          // float sum_weight_propz_propy = 0;

          // looping over 1d hits; propagate them and calulate residuals
          for (auto itDTHits1D: vDTSeg1D) {

          // for ( std::vector<TrackingRecHit*>::const_iterator itDTHits1D =  vDTHits1D.begin(); itDTHits1D != vDTHits1D.end(); ++itDTHits1D ) {
		    // TrackingRecHit* hit = *itDTHits1D;

            DetId hitId  = itDTHits1D->geographicalId(); // get global position of the hit
            const DTSuperLayerId superLayerId(hitId.rawId());
		    const DTLayerId layerId(hitId.rawId());
            DTChamberId chamberId = DTChamberId(hitId);
            // cout << "Chamber ID is " << chamberId.wheel() << " " << chamberId.station() << " " << chamberId.sector() << std::endl;

            //Extrapolate the propagated position for residual
			TrajectoryStateOnSurface tsos_on_chamber;
			tsos_on_chamber = propagator->propagate( tsos_from_tracker, DTGeometry_->idToDet(hitId)->surface() );

            // if ( tsos_on_chamber.isValid() ) {
			//     std::cout << " extrapolation localPosition()"
			// 	<< " x: " << tsos_on_chamber.localPosition().x()
			// 	<< " y: " << tsos_on_chamber.localPosition().y()
			// 	<< " z: " << tsos_on_chamber.localPosition().z() << std::endl;
		    // }

            // filling TTree
            LocalPoint hit_LP = itDTHits1D->localPosition();
            GlobalPoint hit_GP = DTGeometry_->idToDet(hitId)->toGlobal(itDTHits1D->localPosition());
            LocalPoint prop_LP = tsos_on_chamber.localPosition();
            GlobalPoint prop_GP = tsos_on_chamber.globalPosition();

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
            // data_.rechit_location[3] = dtDetId.superlayer();

            if ( superLayerId.superlayer() == 2  && vDTSeg1D.size() >= 3 ) {
                double residual_dy = prop_LP.y() - hit_LP.y();
                data_.dx = 99999;
                data_.dy = residual_dy;
                // cout << "residual dy is" << " " << residual_dy << std::endl;
                if (chamberId.station() == 4) {
                    cout << "WARNING! STATION 4!" << std::endl;
                    cout << "dy WAS CALCULATED!" << std::endl;
                }
                tree->Fill();
            }

            if ((superLayerId.superlayer() == 1 || superLayerId.superlayer() == 3) && vDTSeg1D.size() >= 6) {
                double residual_dx = prop_LP.x() - hit_LP.x();
                data_.dx = residual_dx;
                data_.dy = 99999;
                // cout << "residual dx is" << " " << residual_dx << std::endl;
                if (chamberId.station() == 4 && chamberId.sector() == 9 && chamberId.wheel() == 0) {
                    cout << "**** WARNING! STATION 4! ****" << std::endl;
                    cout << "superlayer ID is " << superLayerId.superlayer() << std::endl;
                    cout << "Chamber ID is " << chamberId.wheel() << " " << chamberId.station() << " " << chamberId.sector() << std::endl;
                    cout << "local positions of hits are " << hit_LP.x() << " " << hit_LP.y() << " " << hit_LP.z() << std::endl;
                    cout << "local prop positions are " << prop_LP.x() << " " << prop_LP.y() << " " << prop_LP.z() << std::endl;
                    cout << "residual is " << residual_dx << std::endl;
                    cout << "muon index is " << data_.muonIdx << std::endl;
                    cout << "**************************" << std::endl;
                }
                tree->Fill();
            }
            // tree->Fill();
        }
    }
// fill chamber-level tree here
}

}







void DT_tbma::beginJob(){}
void DT_tbma::endJob(){}

DEFINE_FWK_MODULE(DT_tbma);
