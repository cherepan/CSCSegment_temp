
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
  
  double* wireHits()   {return wireHitPosition;}
  int*    nLayerHits() {return nHitsInLayer;}

  void   updateWHits(double* wHits2, int* nHits2);
  double comHitLow();
  double comHitHigh();
  TH2F*  getPatternMatrix() { return wirePattern; }
  
  ~CSCWireSegment();

private:

  int theKeyWG;
  int nlayersWithHits;
  double wireHitPosition[6] ;    // wire hit position in each layer
  int    nHitsInLayer[6];               // number of wire hits in each layer

  //  double GetMean(TH1D* h1);  // to be removed, unused
  TH2F* wirePattern; // unused so far

};
#endif
