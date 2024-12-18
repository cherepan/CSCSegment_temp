
#ifndef CSCSegment_CSCWireSegment_H
#define CSCSegment_CSCWireSegment_H

#include "TH2F.h"

class CSCWireSegment
{

public:

  CSCWireSegment();
  CSCWireSegment(int wg, int nLayer, TH2F* wireSegHist);

  int keyWG() const { return theKeyWG; }
  int nLayersWithHits() const { return nlayersWithHits; }
  
  double* wireHitsPosition()   {return HitPosition;}
  int*    nLayerHits() {return nHitsInLayer;}

  void   updateWHits(double* wHits2, int* nHits2);
  void   printWireSegment();
  double LowestHitInLayer();
  double HighestHitInLayer();
  bool   SegmentWithMissingLayers();
  TH2F*  getPatternMatrix() { return wirePattern; }
  
  ~CSCWireSegment();

private:

  int theKeyWG;
  int nlayersWithHits;
  double HitPosition[6] ;           // wire hit position in each layer
  int    nHitsInLayer[6];               // number of wire hits in each layer

  //  double GetMean(TH1D* h1);  // to be removed, unused
  TH2F* wirePattern; // unused so far

};
#endif
