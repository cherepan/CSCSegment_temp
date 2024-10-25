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

       wHitHists[i] = wireSegHist->ProjectionX("layer" + TString(i+1), i+1, i+1);
       HitPosition[i] = wHitHists[i]->GetMean() + 0.5;

       if (wHitHists[i]->GetMean() == 0) HitPosition[i] = 0;
       nHitsInLayer[i] = nHitHists->GetBinContent(i+1);
       
     }

}


CSCWireSegment::~CSCWireSegment() {}



void CSCWireSegment::updateWHits(double* NextSegmentWireHitsPosition, int* NextSegmentNHits)
// This way of merging segments to be tested
// If second segment is 1 WG apart merge two in a big segment, computing in each layer center of gravity of two segments
{

   for (int i = 0; i < 6; i++)
     {
       if( (nHitsInLayer[i] + NextSegmentNHits[i]) != 0)  // if layer is not empty
	 {
	   
	   HitPosition[i] = ( HitPosition[i]*nHitsInLayer[i]  +  NextSegmentWireHitsPosition[i]* NextSegmentNHits [i] ) / ( nHitsInLayer[i] + NextSegmentNHits[i] );
	   
	 }
       
       nHitsInLayer[i] = nHitsInLayer[i] + NextSegmentNHits [i];

     }
}



double CSCWireSegment::LowestHitInLayer()
{

  double low = 121; // why 120 ??    Low is low, why this idiotic logic????   Max N WG = 120

  for (int i = 0; i < 6; i++)
    {
      
      double tmpHit = HitPosition[i];
      
      if (tmpHit < low && HitPosition[i] > 0) low = tmpHit;
      
    }

  return low;

}


double CSCWireSegment::HighestHitInLayer()
{

  double high = -1;

  for (int i = 0; i < 6; i++)
    {
      
      double tmpHit = HitPosition[i];
      
      if (tmpHit > high && HitPosition[i] > 0) high = tmpHit;
      
    }

  return high;

}

void CSCWireSegment::printWireSegment()
{
  
  std::cout<<"  ===>>>>>>>> Print Wire Segment with key WG  "<< theKeyWG << std::endl;
  for (int i = 0; i < 6; i++)
    {
      std::cout<<"  Layer:   "<< i  << "  position:  "<< HitPosition[i] << "  n hits  "<< nHitsInLayer[i] <<std::endl;
    }
  
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
