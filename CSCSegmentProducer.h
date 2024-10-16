#ifndef CSCSegment_CSCSegmentProducer_h
#define CSCSegment_CSCSegmentProducer_h

/** \class CSCSegmentProducer 
 * Produces a collection of CSCSegment's in endcap muon CSCs. 
 *
 */

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/CSCRecHit/interface/CSCRecHit2DCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCWireHitCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCStripHitCollection.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"


class CSCSegmentBuilder; 
class CSCRecoConditions;

class CSCSegmentProducer : public edm::stream::EDProducer<> {
public:
    /// Constructor
    explicit CSCSegmentProducer(const edm::ParameterSet&);
    /// Destructor
    ~CSCSegmentProducer();
    /// Produce the CSCSegment collection
    virtual void produce(edm::Event&, const edm::EventSetup&) override;


private:
    int iev; // events through
    CSCSegmentBuilder* segmentBuilder_;
    CSCRecoConditions* recoConditions_;
    edm::EDGetTokenT<CSCRecHit2DCollection> m_token;
    edm::EDGetTokenT<CSCWireHitCollection> m_token_wire;
    edm::EDGetTokenT<CSCStripHitCollection> m_token_strip;
    edm::ESGetToken<CSCGeometry, MuonGeometryRecord> m_cscGeometryToken;

};

#endif
