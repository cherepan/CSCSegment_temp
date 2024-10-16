#ifndef CSCSegment_CSCStripSegment_H
#define CSCSegment_CSCStripSegment_H

#include "TH2F.h"

class CSCStripSegment 
{

public:

  CSCStripSegment();
  CSCStripSegment(int halfStrip, int nLayer, int patternRank, TH2F* stripSegHist);

  int keyHalfStrip() const {return theKeyHalfStrip;}
  int nLayersWithHits() const {return nlayersWithHits;}
  int patternRank() const {return thePatternRank;}
  double* stripHits() {return sHits;}
  int* nLayerHits() {return nHits;}

  void updateSHits(double* sHits2, int* nHits2);

  double comHitLow(bool isME11);
  double comHitHigh(bool isME11);

  ~CSCStripSegment();

private:
  
  int theKeyHalfStrip;
  int nlayersWithHits;
  int thePatternRank;
  double sHits[6];
  int nHits[6]; // number of shits in each layer

  double GetMean(TH1D* h1);
  // TH2F* stripPattern

};

#endif
