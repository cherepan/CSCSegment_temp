#include <iostream>
#include "RecoLocalMuon/CSCSegment/src/CSCWireSegment.h"

CSCWireSegment::CSCWireSegment() {}

CSCWireSegment::CSCWireSegment(int wg,
                               int nLayer,
                               TH2F* wireSegHist) :
   theKeyWG( wg ),
   nlayersWithHits( nLayer )
   //   wirePattern(wireSegHist)
  
{

   TH1D* wHitHists[6];
   TH1D* nHitHists = wireSegHist->ProjectionY("nLayerHits",0,-1);
   
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
       std::cout<<"CSCWireSegment:  nHits   "<<  nHitHists->GetBinContent(i+1)  << std::endl;
       
     }

}


CSCWireSegment::~CSCWireSegment() {}





void CSCWireSegment::updateWHits(double* wHits2, int* nHits2)  // no clue what it's doing
{

   for (int i = 0; i < 6; i++)
     {
       //       std::cout<<"  nHits[i]+nHits2[i]  " << nHitsInLayer[i] + nHits2[i] <<  "   true or false  "<< (nHitsInLayer[i] + nHits2[i]) << "   wireHitPosition    " << wireHitPosition[i] <<std::endl;

       std::cout<<"============  BEFORE   layer    "<< i << "  nHits    " << nHitsInLayer[i] << "  positions   "<< wireHitPosition[i]  <<std::endl;
       std::cout<<"============  adj      layer    "<< i << "  nHits    " << nHits2[i] << "  positions   "<< wHits2[i]  <<std::endl;
       
       if( (nHitsInLayer[i] + nHits2[i]) != 0)  // if layer is not empty
	 {
	   //	   if((nHitsInLayer[i] + nHits2[i]) == 1) std::cout<<"  nHitsInLayer   "<< nHitsInLayer[i] << " adjacent    "<< nHits2[i] <<std::endl;
	   
	   wireHitPosition[i] = ( wireHitPosition[i]*nHitsInLayer[i]  +  wHits2[i]*nHits2[i] ) / ( nHitsInLayer[i] + nHits2[i] );

	   
	 }
       
       nHitsInLayer[i] = nHitsInLayer[i] + nHits2[i];
       std::cout<<"============  AFTER  layer    "<< i << "  nHits    " << nHitsInLayer[i] << "  positions   "<< wireHitPosition[i]  <<std::endl;
       //       std::cout<<" UPDATED   nHits[i]+nHits2[i]  " << (nHitsInLayer[i]+nHits2[i]) << "   wireHitPosition    " << wireHitPosition[i] <<std::endl;
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
