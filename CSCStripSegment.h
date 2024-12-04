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
  double* stripHits()  {return stripHitsPosition;}
  int*    nLayerHits() {return numberOfHitsInLayer;}
  void printStripSegment();
  void updateSHits(double* sHits2, int* nHits2);

  double LowestHitInLayer(bool isME11);
  double HighestHitInLayer(bool isME11);
  
  bool   SegmentWithMissingLayers();

  
  ~CSCStripSegment();

private:
  
  int theKeyHalfStrip;
  int nlayersWithHits;
  int thePatternRank;
  double stripHitsPosition[6];
  int numberOfHitsInLayer[6]; // number of shits in each layer

  
  double GetMean(TH1D* h1);
  // TH2F* stripPattern

};

#endif
