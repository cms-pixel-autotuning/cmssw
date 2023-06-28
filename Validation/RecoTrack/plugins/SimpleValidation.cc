// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "SimTracker/TrackerHitAssociation/interface/ClusterTPAssociation.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;

class SimpleValidation : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit SimpleValidation(const edm::ParameterSet&);
  ~SimpleValidation() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  
  std::vector<edm::InputTag> trackLabels_;
  edm::EDGetTokenT<ClusterTPAssociation> tpMap_;
//   edm::EDGetTokenT<std::vector<PileupSummaryInfo>>  infoPileUp_;
  std::vector<edm::EDGetTokenT<edm::View<reco::Track>>> trackTokens_;
  edm::EDGetTokenT<reco::TrackToTrackingParticleAssociator> trackAssociatorToken_;
  edm::EDGetTokenT<TrackingParticleCollection> trackingParticleToken_;

//   const double sharingFraction_;
//   const double sharingFractionForTriplets_;


};

SimpleValidation::SimpleValidation(const edm::ParameterSet& iConfig)
    : trackLabels_(iConfig.getParameter<std::vector<edm::InputTag>>("trackLabels")),
      // tpMap_(consumes(iConfig.getParameter<edm::InputTag>("tpMap"))),
    //   infoPileUp_(consumes(iConfig.getParameter< edm::InputTag >("infoPileUp"))),
      trackAssociatorToken_(consumes<reco::TrackToTrackingParticleAssociator>(iConfig.getUntrackedParameter<edm::InputTag>("trackAssociator"))),
      trackingParticleToken_(consumes<TrackingParticleCollection>(iConfig.getParameter< edm::InputTag >("trackingParticles")))
    //   sharingFraction_(iConfig.getUntrackedParameter<double>("sharingFraction")),
    //   sharingFractionForTriplets_(iConfig.getUntrackedParameter<double>("sharingFractionForTriplets"))
{

  for (auto& itag : trackLabels_) {
    trackTokens_.push_back(consumes<edm::View<reco::Track>>(itag));
  }

  //now do what ever initialization is needed
}

SimpleValidation::~SimpleValidation() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void SimpleValidation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
   

//   auto const& tpClust = iEvent.get(tpMap_);
  auto const& associatorByHits = iEvent.get(trackAssociatorToken_);
  
  TrackingParticleRefVector tpCollection;
  edm::Handle<TrackingParticleCollection> TPCollectionH;
  iEvent.getByToken(trackingParticleToken_, TPCollectionH);
//   auto const& tp = iEvent.get(trackingParticleToken_);

  for (size_t i = 0, size = TPCollectionH->size(); i < size; ++i) {
   tpCollection.push_back(TrackingParticleRef(TPCollectionH, i));
  }
  
  for (const auto& trackToken : trackTokens_) 
  {

    edm::Handle<edm::View<reco::Track>> tracksHandle;
    iEvent.getByToken(trackToken, tracksHandle);
    const edm::View<reco::Track>& tracks = *tracksHandle;
    
    edm::RefToBaseVector<reco::Track> trackRefs;
    for (edm::View<reco::Track>::size_type i = 0; i < tracks.size(); ++i) {
        trackRefs.push_back(tracks.refAt(i));
    }

    reco::RecoToSimCollection recSimColl = associatorByHits.associateRecoToSim(trackRefs, tpCollection);
    int rt = 0;
    int at = 0;
    for (const auto& track : trackRefs) {
        rt++;
        // int charge = track.charge();
        // float pt = track.pt();
        // std::cout << pt << std::endl;
        // int nSimHIts = 0;
        // std::vector<int> tpIdx;
        // std::vector<float> sharedFraction;
        // std::vector<float> tpChi2;
        // bool isSimMatched = false;
        auto foundTPs = recSimColl.find(track);
        if (foundTPs != recSimColl.end()) {
          if (!foundTPs->val.empty()) {
            at++;
            // nSimHits = foundTPs->val[0].first->numberOfTrackerHits();
            // isSimMatched = true;
          }
        }
        // //tP Matching
        // auto rangeIn = tpClust->equal_range(hits[0]->firstClusterRef());
        // auto rangeOut = tpClust->equal_range(hits[1]->firstClusterRef());
    }
    LogPrint("TrackValidator") << "Associated tracks" << at << "Total reconstructed" << rt <<'\n';
  }

}

// ------------ method called once each job just before starting event loop  ------------
void SimpleValidation::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void SimpleValidation::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void SimpleValidation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SimpleValidation);
