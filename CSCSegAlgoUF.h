#ifndef CSCSegment_CSCSegAlgoUF_h
#define CSCSegment_CSCSegAlgoUF_h

/**
 * \class CSCSegAlgoUF
 *

 */


#include "CSCSegmentAlgorithm.h"

#include <DataFormats/CSCRecHit/interface/CSCRecHit2D.h>
#include <DataFormats/CSCRecHit/interface/CSCWireHit.h>
#include <DataFormats/CSCRecHit/interface/CSCStripHit.h>
#include "RecoLocalMuon/CSCSegment/src/CSCRecoConditions.h"
#include "RecoLocalMuon/CSCSegment/src/CSCMake2DRecHit.h"
#include "RecoLocalMuon/CSCSegment/src/CSCXonStrip_MatchGatti.h"
#include "RecoLocalMuon/CSCSegment/src/CSCFindPeakTime.h"


#include "CSCSegFit.h"
#include "CSCWireSegment.h"
#include "CSCStripSegment.h"

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>

#include <vector>

#include "TH2F.h"

class CSCSegFit;
class CSCRecoConditions;
class CSCMake2DRecHit;

class CSCSegAlgoUF : public CSCSegmentAlgorithm {

public:


    typedef std::vector<int> LayerIndex;
    typedef std::vector<const CSCRecHit2D*> ChamberHitContainer;
    typedef std::vector<const CSCWireHit*>  ChamberWireHitContainer;
    typedef std::vector<const CSCStripHit*> ChamberStripHitContainer;
    typedef std::vector<const CSCRecHit2D*>::const_iterator ChamberHitContainerCIt;

    

    explicit CSCSegAlgoUF(const edm::ParameterSet& ps);

    ~CSCSegAlgoUF();


    std::vector<CSCSegment> buildSegments(const ChamberHitContainer      & rechits);
    std::vector<CSCSegment> buildSegments(const ChamberWireHitContainer  & wirehits, 
                                          const ChamberStripHitContainer & striphits);

    void FillWireMatrix( TH2F* whitsMatrix, ChamberWireHitContainer  whits);
    void FillStripMatrix(TH2F* shitsMatrix, ChamberStripHitContainer shits);

  
    void ScanForWireSeg( TH2F* wireHitsInChamber, std::list<CSCWireSegment>& wireSegments, std::vector<TH2F*>& wireSegmentsTH2F, std::vector<int>& wireSegements_rank/* <---- for debugging, to be removed */,  int nLayer);
  void ScanForStripSeg(TH2F* stripHitsInChamber, std::list<CSCStripSegment>& stripSegments, std::vector<TH2F*>& stripSegmentsTH2F, std::vector<int>& stripSegments_rank/* <---- for debugging, to be removed */, int nLayer);

  
    void GetWireHitFromWireSeg  (CSCWireSegment  wireSeg,  ChamberWireHitContainer  whits, int* wireHitIndex);
    void GetStripHitFromStripSeg(CSCStripSegment stripSeg, ChamberStripHitContainer shits, int* stripHitIndex);


    std::vector< CSCSegment > prune_bad_hits(const CSCChamber* aChamber, std::vector< CSCSegment > & segments);


  
    void PrintTH2F(TH2F* hist);

  

    std::vector<CSCSegment> run(const CSCChamber* aChamber, const ChamberHitContainer& rechits, 
                                                            const ChamberWireHitContainer& wirehits,
                                                            const ChamberStripHitContainer& striphits,
                                                            CSCRecoConditions* reco); 

    CSCMake2DRecHit*       make2DHits_;

    void setConditions ( CSCRecoConditions* reco );

private:

  
    bool condpass1, condpass2;  // 
    ChamberHitContainer protoSegment;
    const edm::ParameterSet& pset(void) const { return ps_;}
    const edm::ParameterSet ps_;
  
    double chi2Norm_3D_;            /// Chi^2 normalization for the initial fit
    bool prePrun_;                  /// Allow to prune a (rechit in a) segment in segment buld method
                                    /// once it passed through Chi^2-X and  chi2Correction is big.
    double prePrunLimit_;           //  what the hell is this !!!!!!!??????????? It is not even initialised but somehow used in the code!!!!!!!!!!!!!

  
    int nWireGroups;
    int nStrips;



    const CSCChamber* theChamber;
    int theStation;
    int theRing;
    bool isME11;
  
    const std::string myName; 
		
    double theChi2;
    LocalPoint theOrigin;
    LocalVector theDirection;
  
    std::unique_ptr<CSCSegFit> sfit_;
    CSCRecoConditions* recoConditions_;

};

#endif
