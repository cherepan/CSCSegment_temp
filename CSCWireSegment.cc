#include <iostream>
#include "RecoLocalMuon/CSCSegment/src/CSCWireSegment.h"

CSCWireSegment::CSCWireSegment() {}

CSCWireSegment::CSCWireSegment(int   wg,
                               int   nLayer,
                               TH2F* wireSegHist) :
   theKeyWG( wg ),
   nlayersWithHits( nLayer )
   //   wirePattern(wireSegHist)
  
{

   TH1D* wHitHists[6];
   TH1D* nHitHists = wireSegHist->ProjectionY("nLayerHits", 0, -1);
   
   for (int i = 0; i < 6; i++) 
     {

       wHitHists[i] = wireSegHist->ProjectionX("layer"+TString(i+1), i+1, i+1);
       wireHitPosition[i] = wHitHists[i]->GetMean() + 0.5;
       std::cout<<"CSCWireSegment:  csc wire segment  "<< wHitHists[i]->GetMean() + 0.5  << std::endl;
       
       if (wHitHists[i]->GetMean() == 0) wireHitPosition[i] = 0;

       
       // wHits[i] = GetMean(wHitHists[i]);
       // first number is layer, second number is wirePosition in unit of wire group
       // +0.5 is because of ROOT histogram property, GetMean() returns 4.5 for example if fill only one entry at 4 ??
       // Vladimir: That's not true;

       
       nHitsInLayer[i] = nHitHists->GetBinContent(i+1);
       //       std::cout<<"CSCWireSegment:  nHits   "<<  nHitHists->GetBinContent(i+1)  << std::endl;
       
     }

}


CSCWireSegment::~CSCWireSegment() {}





void CSCWireSegment::updateWHits(double* NextSegmentWireHitsPosition, int* NextSegmentNHits)
// This way of merging segments to be tested.
// If second segment is 1 WG apart merge two in a big segment, computing in each layer center of gravity of two segments
{

   for (int i = 0; i < 6; i++)
     {
       if( (nHitsInLayer[i] + NextSegmentNHits[i]) != 0)  // if layer is not empty
	 {
	   
	   wireHitPosition[i] = ( wireHitPosition[i]*nHitsInLayer[i]  +  NextSegmentWireHitsPosition[i]* NextSegmentNHits [i] ) / ( nHitsInLayer[i] + NextSegmentNHits[i] );
	   
	 }
       
       nHitsInLayer[i] = nHitsInLayer[i] + NextSegmentNHits [i];

     }
}



double CSCWireSegment::comHitLow()
{

  double low = 120;

  for (int i = 0; i < 6; i++)
    {
      double tmpHit = wireHitPosition[i];
      if (tmpHit < low && wireHitPosition[i] > 0) low = tmpHit;
      
    }

  return low;

}


double CSCWireSegment::comHitHigh()
{

  double high = -1;

  for (int i = 0; i < 6; i++)
    {
      
      double tmpHit = wireHitPosition[i];
      
      if (tmpHit > high && wireHitPosition[i] > 0) high = tmpHit;
      
    }

  return high;

}



/*
double CSCWireSegment::GetMean(TH1D* h1)
{

  double mean = 0;
  double count = 0;
  double sum = 0;

//std::cout << h1->GetNbinsX() << std::endl;
  for (int i = 0; i < h1->GetNbinsX(); i++) {

      if (h1->GetBinContent(i+1) > 0) {
         count += 1;
	 sum += (h1->GetBinLowEdge(i+1) + 1);
//std::cout << "count: " << count << ", sum: " << sum << std::endl;
         }
      }

  if (count > 0) mean = sum/count;
//std::cout << mean << std::endl;
  return mean;
}
*/
