#include "RecoLocalMuon/CSCSegment/src/CSCStripSegment.h"
#include <iostream>

CSCStripSegment::CSCStripSegment() {}

CSCStripSegment::CSCStripSegment(int halfStrip,
                               int nLayer,
                               int patternRank,
                               TH2F* stripSegHist) :
   theKeyHalfStrip( halfStrip   ),
   nlayersWithHits( nLayer      ),
   thePatternRank(  patternRank )

{


   TH1D* sHitHists[6];
   TH1D* nHitHists = stripSegHist->ProjectionY("nLayerHits", 0, -1);
   for (int i = 0; i < 6; i++)
     {

       sHitHists[i] = stripSegHist->ProjectionX("layer"+TString(i+1),i+1,i+1);
       stripHitsPosition[i] = sHitHists[i]->GetMean()+0.5;
       //       sHits[i] = GetMean(sHitHists[i]);
       if (sHitHists[i]->GetMean() == 0) stripHitsPosition[i] = 0;
       numberOfHitsInLayer[i] = nHitHists->GetBinContent(i+1);
       //std::cout << "sHitHists[i]->GetMean(): " << sHitHists[i]->GetMean() << std::endl;
     }

}


CSCStripSegment::~CSCStripSegment() {}



void CSCStripSegment::updateSHits(double* NextSegmentStripHitsPosition, int* NextSegmentNHits)
{

   for (int i = 0; i < 6; i++)
     {
       
       if (numberOfHitsInLayer[i] + NextSegmentNHits[i] > 0) 
	 {
	   stripHitsPosition[i] = (stripHitsPosition[i]*numberOfHitsInLayer[i] + NextSegmentStripHitsPosition[i]*NextSegmentNHits[i] ) / (numberOfHitsInLayer[i] + NextSegmentNHits[i]);
	 }
       numberOfHitsInLayer[i] = numberOfHitsInLayer[i] + NextSegmentNHits[i];
       
     }

}


double CSCStripSegment::comHitLow(bool isME11)
{

  double low = 161;

  for (int i = 0; i < 6; i++) {

      double tmpHit = stripHitsPosition[i];
      if (!isME11 && (i==0 || i==2 || i==4)) tmpHit -= 1;
      if (tmpHit < low && stripHitsPosition[i] > 0) low = tmpHit;

      }

  return low;

}


double CSCStripSegment::comHitHigh(bool isME11)
{

  double high = -1;

  for (int i = 0; i < 6; i++) {

      double tmpHit = stripHitsPosition[i];
      if (!isME11 && (i==0 || i==2 || i==4)) tmpHit -= 1;
      if (tmpHit > high && stripHitsPosition[i] > 0) high = tmpHit;

      }

  return high;

}



/*
double CSCStripSegment::GetMean(TH1D* h1)
{

  double mean = 0;
  double count = 0;
  double sum = 0;
  for (int i = 0; i < h1->GetNbinsX(); i++) {

      if (h1->GetBinContent(i+1) > 0) {
         count += 1; sum += (h1->GetBinLowEdge(i+1) + 1);
         }
      }

  if (count>0) mean = sum/count;
  return mean;
}

*/
