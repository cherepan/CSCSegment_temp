/** \file CSCSegmentProducer.cc
 *
 */


#include "CSCRecoConditions.h"
#include <RecoLocalMuon/CSCSegment/src/CSCSegmentProducer.h>
#include <RecoLocalMuon/CSCSegment/src/CSCSegmentBuilder.h>


#include <DataFormats/Common/interface/Handle.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/Utilities/interface/InputTag.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h> 

#include <Geometry/Records/interface/MuonGeometryRecord.h>

#include <DataFormats/CSCRecHit/interface/CSCSegmentCollection.h>
#include <DataFormats/CSCRecHit/interface/CSCSegment.h>

CSCSegmentProducer::CSCSegmentProducer(const edm::ParameterSet& pas) : iev(0) {
	
    m_token         = consumes<CSCRecHit2DCollection>( pas.getParameter<edm::InputTag>("inputObjects") );
    m_token_wire    = consumes<CSCWireHitCollection> ( pas.getParameter<edm::InputTag>("inputObjects") );
    m_token_strip   = consumes<CSCStripHitCollection>( pas.getParameter<edm::InputTag>("inputObjects") );
    segmentBuilder_ = new CSCSegmentBuilder(pas); // pass on the PS
    m_cscGeometryToken = esConsumes<CSCGeometry, MuonGeometryRecord>();

    recoConditions_    = new CSCRecoConditions( pas,  consumesCollector() ); // access to conditions data
    segmentBuilder_    ->setConditions( recoConditions_ ); // pass down to who needs access

  	// register what this produces
    produces<CSCSegmentCollection>();
}

CSCSegmentProducer::~CSCSegmentProducer() {

    LogDebug("CSCSegment|CSC") << "deleting CSCSegmentBuilder after " << iev << " events w/csc data.";
    delete segmentBuilder_;
    delete recoConditions_;

}

void CSCSegmentProducer::produce(edm::Event& ev, const edm::EventSetup& setup) {

    LogDebug("CSCSegment|CSC") << "start producing segments for " << ++iev << "th event with csc data";

    std::cout<<" ===> Event:  "<< ev.id().event() <<std::endl;
    // find the geometry (& conditions?) for this event & cache it in the builder
  
    edm::ESHandle<CSCGeometry> h = setup.getHandle(m_cscGeometryToken);
    const CSCGeometry* pgeom = &*h;
    segmentBuilder_->setGeometry(pgeom);



    //    edm::ESHandle<CSCGeometry> h;
    //    setup.get<MuonGeometryRecord>().get(h);
    //    const CSCGeometry* pgeom = &*h;

    segmentBuilder_->setGeometry(pgeom);
    recoConditions_->initializeEvent( setup );
	
    // get the collection of CSCRecHit2D
    edm::Handle<CSCRecHit2DCollection> cscRecHits;
    edm::Handle<CSCWireHitCollection>  cscWireHits;
    edm::Handle<CSCStripHitCollection> cscStripHits;

    
    ev.getByToken( m_token,       cscRecHits);
    ev.getByToken( m_token_wire,  cscWireHits);
    ev.getByToken( m_token_strip, cscStripHits);


    // create empty collection of Segments
    auto oc = std::make_unique<CSCSegmentCollection>();


    segmentBuilder_->build(cscRecHits.product(), cscWireHits.product(), cscStripHits.product(), *oc); //@@ FILL oc


    // put collection in event
    ev.put(std::move(oc));
}

