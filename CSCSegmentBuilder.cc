
/** \file CSCSegmentBuilder.cc
 *
 *
 */

#include <RecoLocalMuon/CSCSegment/src/CSCSegmentBuilder.h>
#include <Geometry/CSCGeometry/interface/CSCChamberSpecs.h>
#include <Geometry/CSCGeometry/interface/CSCLayer.h>
#include <Geometry/CSCGeometry/interface/CSCGeometry.h>
#include <DataFormats/MuonDetId/interface/CSCDetId.h>
#include <DataFormats/CSCRecHit/interface/CSCRecHit2D.h>
#include <DataFormats/CSCRecHit/interface/CSCRangeMapAccessor.h>
#include <DataFormats/CSCRecHit/interface/CSCWireHit.h>
#include <DataFormats/CSCRecHit/interface/CSCStripHit.h>

#include <RecoLocalMuon/CSCSegment/src/CSCSegmentAlgorithm.h>
#include <RecoLocalMuon/CSCSegment/src/CSCSegmentBuilderPluginFactory.h>

#include <FWCore/Utilities/interface/Exception.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h> 



CSCSegmentBuilder::CSCSegmentBuilder(const edm::ParameterSet& ps) : geom_(0) {
    
    // The algo chosen for the segment building
    int chosenAlgo = ps.getParameter<int>("algo_type") - 1;
    
    // Find appropriate ParameterSets for each algo type
    std::vector<edm::ParameterSet> algoPSets = ps.getParameter<std::vector<edm::ParameterSet> >("algo_psets");

    // Now load the right parameter set
    // Algo name
    std::string algoName = algoPSets[chosenAlgo].getParameter<std::string>("algo_name");
    //    std::cout<<" >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> algo_name : " << algoName << std::endl;
    LogDebug("CSCSegment|CSC")<< "CSCSegmentBuilder algorithm name: " << algoName;

    // SegAlgo parameter set
    std::vector<edm::ParameterSet> segAlgoPSet = algoPSets[chosenAlgo].getParameter<std::vector<edm::ParameterSet> >("algo_psets");

    // Chamber types to handle
    std::vector<std::string> chType = algoPSets[chosenAlgo].getParameter<std::vector<std::string> >("chamber_types");
    LogDebug("CSCSegment|CSC")<< "No. of chamber types to handle: " << chType.size();

    // Algo to chamber type 
    std::vector<int> algoToType = algoPSets[chosenAlgo].getParameter<std::vector<int> >("parameters_per_chamber_type");

    // Trap if we don't have enough parameter sets or haven't assigned an algo to every type   
    if (algoToType.size() !=  chType.size()) 
      {
        throw cms::Exception("ParameterSetError") << 
	  "#dim algosToType=" << algoToType.size() << ", #dim chType=" << chType.size() << std::endl;
      }

    // Ask factory to build this algorithm, giving it appropriate ParameterSet
            
    for (size_t j=0; j<chType.size(); ++j) 
      {
	//	algoMap[chType[j]] = CSCSegmentBuilderPluginFactory::get()->create(algoName, segAlgoPSet[algoToType[j]-1]);  
	algoMap.emplace(chType[j], CSCSegmentBuilderPluginFactory::get()->create(algoName, segAlgoPSet[algoToType[j] - 1])); // CMSSW_12_4_4 
	edm::LogVerbatim("CSCSegment|CSC")<< "using algorithm #" << algoToType[j] << " for chamber type " << chType[j];
      }
}

CSCSegmentBuilder::~CSCSegmentBuilder() = default; 


/*
CSCSegmentBuilder::~CSCSegmentBuilder() {
  //
  // loop on algomap and delete them
  //


  //  for (std::map<std::string, std::unique_ptr<CSCSegmentAlgorithm> >::iterator it = algoMap.begin();it != algoMap.end(); it++)    delete *(*it).second; 
  for (std::map<std::string, std::unique_ptr<CSCSegmentAlgorithm> >::iterator it = algoMap.begin();it != algoMap.end(); it++)    delete it->second; 
  
}
*/






void CSCSegmentBuilder::build(const CSCRecHit2DCollection* recHits,
                              const CSCWireHitCollection*  wireHits,
                              const CSCStripHitCollection* stripHits, CSCSegmentCollection& oc) {
  	
  //  edm::LogVerbatim("CSCSegment|CSC")<< "Total number of rechits in this event: " << recHits->size();
  //  edm::LogVerbatim("CSCSegment|CSC")<< "Total number of wirehits in this event: " << wireHits->size();
  //  edm::LogVerbatim("CSCSegment|CSC")<< "Total number of striphits in this event: " << stripHits->size();

  //  std::cout << "Total number of rechits in this event: " << recHits->size() << std::endl;
  //  std::cout << "Total number of wirehits in this event: " << wireHits->size() << std::endl;
  //  std::cout << "Total number of striphits in this event: " << stripHits->size() << std::endl;

    std::vector<CSCDetId> chambers;
    std::vector<CSCDetId>::const_iterator chIt;
    
//    for(CSCRecHit2DCollection::const_iterator it2 = recHits->begin(); it2 != recHits->end(); it2++) {

    // find chambers with strip hits
    for(CSCStripHitCollection::const_iterator it2 = stripHits->begin(); it2 != stripHits->end(); it2++) 
      {
	
        bool insert = true;
        for(chIt=chambers.begin(); chIt != chambers.end(); ++chIt) 
	  if (((*it2).cscDetId().chamber() == (*chIt).chamber()) &&
	      ((*it2).cscDetId().station() == (*chIt).station()) &&
	      ((*it2).cscDetId().ring()    == (*chIt).ring())    &&
	      ((*it2).cscDetId().endcap()  == (*chIt).endcap()))
	    insert = false;
	
        if (insert)
	  chambers.push_back((*it2).cscDetId().chamberId());
      }
    
    // loop over chamber with strip hits, use it for segment building if it also has wire hits
    for(chIt=chambers.begin(); chIt != chambers.end(); ++chIt) 
      {

        std::vector<const CSCRecHit2D*> cscRecHits;
        std::vector<const CSCWireHit *> cscWireHits;
        std::vector<const CSCStripHit*> cscStripHits;

        const CSCChamber* chamber = geom_->chamber(*chIt);
        
        CSCRangeMapAccessor acc;
        CSCRecHit2DCollection::range range   = recHits->get(acc.cscChamber(*chIt));
        CSCWireHitCollection::range  range_w = wireHits->get(acc.cscChamber(*chIt));
        CSCStripHitCollection::range range_s = stripHits->get(acc.cscChamber(*chIt));
        
        std::vector<int> hitPerLayer(6);
        for(CSCRecHit2DCollection::const_iterator rechit = range.first; rechit != range.second; rechit++) {
            
            hitPerLayer[(*rechit).cscDetId().layer()-1]++;
            cscRecHits.push_back(&(*rechit));
        }    

        std::vector<int> wHitPerLayer(6);
        for(CSCWireHitCollection::const_iterator wirehit = range_w.first; wirehit != range_w.second; wirehit++) {

            wHitPerLayer[(*wirehit).cscDetId().layer()-1]++;
            cscWireHits.push_back(&(*wirehit));
        }

        std::vector<int> sHitPerLayer(6);
        for(CSCStripHitCollection::const_iterator striphit = range_s.first; striphit != range_s.second; striphit++) {

            sHitPerLayer[(*striphit).cscDetId().layer()-1]++;
            cscStripHits.push_back(&(*striphit));
        }

        if (cscWireHits.size() == 0) continue;
        
        LogDebug("CSCSegment|CSC") << "found " << cscRecHits.size() << " rechits in chamber " << *chIt;

//        std::cout << std::endl;
//        std::cout << "found " << cscRecHits.size() << " rechits in chamber " << *chIt << std::endl;
//        std::cout << "found " << cscWireHits.size() << " wirehits in chamber " << *chIt << std::endl;
//        std::cout << "found " << cscStripHits.size() << " striphits in chamber " << *chIt << std::endl;
        // given the chamber select the appropriate algo... and run it
//        if (chamber->specs()->chamberTypeName() == "CSCSegAlgoUF")
//           algoMap[chamber->specs()->chamberTypeName()]->setConditions(recoConditions_);
        
        std::vector<CSCSegment> segv = algoMap[chamber->specs()->chamberTypeName()]->run(chamber, cscRecHits, cscWireHits, cscStripHits, recoConditions_);

	//        LogDebug("CSCSegment|CSC") << "found " << segv.size() << " segments in chamber " << *chIt;
	std::cout<<"  ===================== FOUND   " << segv.size() << " segments in chamber " << *chIt << std::endl<<std::endl<<std::endl;
        // Add the segments to master collection
        oc.put((*chIt), segv.begin(), segv.end());
    }
}

void CSCSegmentBuilder::setGeometry(const CSCGeometry* geom) {
	geom_ = geom;
}

void CSCSegmentBuilder::setConditions( CSCRecoConditions* reco ) {
  recoConditions_ = reco;
}

