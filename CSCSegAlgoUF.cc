/**
 * \file CSCSegAlgUF.cc
 *
 * 

 */

#include "CSCSegAlgoUF.h"
#include "CSCSegFit.h"
#include "CSCWireSegmentPattern.h"
#include "CSCStripSegmentPattern.h"
#include "Geometry/CSCGeometry/interface/CSCLayer.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>




#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "CSCSegAlgoST.h"
#include "CSCCondSegFit.h"

CSCSegAlgoUF::CSCSegAlgoUF(const edm::ParameterSet& ps)
  : CSCSegmentAlgorithm(ps), ps_(ps), myName("CSCSegAlgoUF"), sfit_(nullptr)
{

  
  make2DHits_ = new CSCMake2DRecHit( ps );

  
}


CSCSegAlgoUF::~CSCSegAlgoUF()
{

  delete make2DHits_;

}


void CSCSegAlgoUF::setConditions( CSCRecoConditions* reco )
{
  recoConditions_ = reco;
  make2DHits_->setConditions( reco );
}


std::vector<CSCSegment> CSCSegAlgoUF::run(const CSCChamber* Chamber, const ChamberHitContainer& rechits,
					  const ChamberWireHitContainer& wirehits,
					  const ChamberStripHitContainer& striphits,
					  CSCRecoConditions* reco)
{
  setConditions(reco);
  theChamber = Chamber;
  
  //  std::cout << theChamber->id() << std::endl;
  //  get number of wire groups and strips for this chamber
  nWireGroups = theChamber->layer(1)->geometry()->numberOfWireGroups();
  nStrips     = theChamber->layer(1)->geometry()->numberOfStrips();


  theStation = theChamber->id().station();
  theRing    = theChamber->id().ring();


  
  if (theStation==1 && (theRing == 1 || theRing == 4) ) // 
    //if (theStation==1 && theRing == 1)  
    {isME11 = true;}
  else {isME11 = false;}

  
  return buildSegments(wirehits, striphits);
  //return buildSegments(rechits);

}


std::vector<CSCSegment> CSCSegAlgoUF::buildSegments(const ChamberWireHitContainer&  uwirehits,
                                                    const ChamberStripHitContainer& ustriphits) {

  std::vector<CSCSegment> segments;

  //  int e = theChamber->id().endcap();
  //  int s = theChamber->id().station();
  //  int r = theChamber->id().ring();
  //  int c = theChamber->id().chamber();
  //  if (!(e == 1 && s == 2 && r == 1 && c == 16)) return segments;



  
  ChamberWireHitContainer  wirehits   = uwirehits;   //  input wire hits
  ChamberStripHitContainer striphits  = ustriphits;  //  input strip hits

  
  TH2F* wireHitsInChamber = new TH2F("wireHitsInChamber", "", nWireGroups, 0, nWireGroups, 6 ,0, 6);
  TH2F* sHitsPerChamber = new TH2F("sHitsPerChamber", "", 2*nStrips+1, 0, 2*nStrips+1, 6 ,0, 6); // half strip stagger at 1,3,5 layer

  
  FillWireMatrix(wireHitsInChamber,  wirehits);   TH2F* wireHitsInChamber_clone = (TH2F*)wireHitsInChamber->Clone("wireHitsInChamber_clone");
  FillStripMatrix(sHitsPerChamber, striphits);  //TH2F* sHitsPerChamber_clone = (TH2F*)sHitsPerChamber->Clone("sHitsPerChamber_clone");


  
  std::cout<<">>>>>>>>>>>>>>>>>>>>>>>  Before Scan for a Segment  "<< std::endl;  WriteTH2F(wireHitsInChamber_clone);
  
  std::list<CSCWireSegment> wireSegs;
  ScanForWireSeg(wireHitsInChamber, wireSegs, 6, true);
  //  ScanForWireSeg(wireHitsInChamber, wireSegs, 5, false);
  //  ScanForWireSeg(wireHitsInChamber, wireSegs, 4, false);
  //  ScanForWireSeg(wireHitsInChamber, wireSegs, 3, false);

  //  std::cout<<">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  End SCan for a wire  Segment  "<< "   nSegments   "  << wireSegs.size() << std::endl;  WriteTH2F(wHitsPerChamber_clone);

  
  std::list<CSCStripSegment> stripSegs;
  ScanForStripSeg(sHitsPerChamber, stripSegs, 6);
  ScanForStripSeg(sHitsPerChamber, stripSegs, 5);
  ScanForStripSeg(sHitsPerChamber, stripSegs, 4);
  ScanForStripSeg(sHitsPerChamber, stripSegs, 3);


  
  std::cout << theChamber->id() <<  ", nWireSegments: " << wireSegs.size() << 
    ", nStripSegments: " << stripSegs.size() << 
    ", nSegments: " << wireSegs.size() * stripSegs.size() << std::endl;

  
  
  for (auto i_wire = wireSegs.begin(); i_wire != wireSegs.end(); i_wire++)
    {
    for (auto i_strip = stripSegs.begin(); i_strip != stripSegs.end(); i_strip++)
      {

	CSCWireSegment  wireSeg  = *i_wire;
	CSCStripSegment stripSeg = *i_strip;

	
	//	std::cout << "keyWG: " << i_wire->keyWG() << ", nLayer: " << i_wire->nLayersWithHits() << std::endl;
	//	std::cout << "keyHS: " << i_strip->keyHalfStrip() << ", nLayer: " << i_strip->nLayersWithHits() << std::endl;

	
	int wHitsFromWSeg[6] = {};        GetWireHitFromWireSeg(wireSeg,    wirehits,  wHitsFromWSeg); 
	int sHitsFromSSeg[6] = {};        GetStripHitFromStripSeg(stripSeg, striphits, sHitsFromSSeg);

	
	std::vector<CSCRecHit2D> csc2DRecHits;

	
	for (int k = 0; k < 6; k++) // loop over layers
	    {
              if (wHitsFromWSeg[k] == -1 || sHitsFromSSeg[k] == -1) continue;
              const CSCWireHit* cscwirehit   = wirehits[wHitsFromWSeg[k]]; 
              const CSCStripHit* cscstriphit = striphits[sHitsFromSSeg[k]];
	      
              const CSCWireHit&  wirehit  = *cscwirehit;
              const CSCStripHit& striphit = *cscstriphit;
              const CSCDetId&    detId    = CSCDetId(theChamber->id().endcap(),
                                               theChamber->id().station(),
                                               theChamber->id().ring(),
                                               theChamber->id().chamber(), k+1);
	      
              const CSCLayer* cscLayer = theChamber->layer(k+1);


	      
              CSCRecHit2D rechit = make2DHits_->hitFromStripAndWire(detId, cscLayer, wirehit, striphit );   //to be reviewed


	      
	      //	      std::cout << " ============================================   CSC UF Algo    print out the rechit              " <<std::endl;
	      //	      cscwirehit->print();
              if (make2DHits_->isHitInFiducial( cscLayer, rechit ))    csc2DRecHits.push_back(rechit);
	    }
	  
	//	  std::cout << csc2DRecHits.size() << " RHs survived fiducial cut" << std::endl;
	  ChamberHitContainer csc2DRecHits_p;
	  
	  for (std::vector<CSCRecHit2D>::const_iterator it = csc2DRecHits.begin(); it != csc2DRecHits.end(); it++)     csc2DRecHits_p.push_back(&(*it)); 
	  
	  if ( int(csc2DRecHits_p.size() ) < 3 ) continue; // why number of hits is 0 ??? or smaller than3, wire and strip hit not found ???
	  







// borrow from ST  //  What this part of code is doing ???
	  
    CSCCondSegFit* segfit = new CSCCondSegFit( pset(), theChamber, csc2DRecHits_p );
    condpass1 = false;
    condpass2 = false;
    segfit->setScaleXError( 1.0 );
    segfit->fit(condpass1, condpass2);

    if( true ){
      if(segfit->chi2()/segfit->ndof()>chi2Norm_3D_){
        condpass1 = true;
        segfit->fit(condpass1, condpass2);
      }
      
      if(segfit->scaleXError() < 1.00005){
        LogTrace("CSCWeirdSegment") << "[CSCSegAlgoST::buildSegments] Segment ErrXX scaled and refit " << std::endl;
        if(segfit->chi2()/segfit->ndof()>chi2Norm_3D_){
          LogTrace("CSCWeirdSegment") << "[CSCSegAlgoST::buildSegments] Segment ErrXY changed to match cond. number and refit " << std::endl;
          condpass2 = true;
          segfit->fit(condpass1, condpass2);
        }
      }

      
      if(prePrun_ && (sqrt(segfit->scaleXError())>prePrunLimit_) &&
         (segfit->nhits()>3)){   
        LogTrace("CSCWeirdSegment") << "[CSCSegAlgoST::buildSegments] Scale factor chi2uCorrection too big, pre-Prune and refit " << std::endl;
        csc2DRecHits_p.erase(csc2DRecHits_p.begin() + segfit->worstHit(),
                           csc2DRecHits_p.begin() + segfit->worstHit()+1 );
        double tempcorr = segfit->scaleXError(); // save current value
        delete segfit;    
        segfit = new CSCCondSegFit( pset(), theChamber, csc2DRecHits_p );
        segfit->setScaleXError( tempcorr ); // reset to previous value (rather than init to 1)
        segfit->fit(condpass1, condpass2);
      }
    }

    CSCSegment temp(csc2DRecHits_p, segfit->intercept(), segfit->localdir(),
                       segfit->covarianceMatrix(), segfit->chi2() );
    delete segfit;
    segments.push_back(temp); 

    
      }
    }

//  std::cout << "return " << segments.size() << " 2D segments" << std::endl;
//  std::cout << std::endl;

  std::cout << "n2DSeg before prune: " << segments.size() <<  std::endl;
  for(auto iseg : segments)std::cout<<"  nRecHits  "<< (iseg.recHits()).size() << std::endl;
  std::vector<CSCSegment> segments_prune = prune_bad_hits(theChamber, segments); 

  std::cout << "n2DSeg after prune: " << segments_prune.size() << std::endl;
  for(auto iseg : segments)std::cout<<"  nRecHits  "<< (iseg.recHits()).size() << std::endl;
  return segments_prune;

}





void CSCSegAlgoUF::WriteTH2F(TH2F* hist) {

     for (int i = 1; i < hist->GetNbinsY()+1; i++)
       {
         for (int j = 1; j < hist->GetNbinsX()+1; j++)
	   {

             if (hist->GetBinContent(j,i)==0)
  	       {
		 std::cout << "-";
	       }
	     else
	       {
		 //		 std::cout << hist->GetBinContent(j,i);
		 std::cout << "+"; hist->GetBinContent(j,i);
	       }

	 }
         std::cout << std::endl;
       }
     
}





void CSCSegAlgoUF::FillWireMatrix(TH2F* whitsMatrix, ChamberWireHitContainer whits)
{
  
  std::vector<double> rows_v;
  std::vector<double> cols_v;
  std::vector<double> data_v;
  
  for (unsigned int i = 0; i < whits.size(); i++)
    {
	 
      const CSCWireHit* whit = whits[i];
      int wLayer = whit->cscDetId().layer();
	 
      for (unsigned int j = 0; j < whit->wgroups().size(); j++)
	{
	     
	  int wg = (whit->wgroups())[j];
	  cols_v.push_back(wg-1);

	  rows_v.push_back(wLayer-1);

	  data_v.push_back(1);
	  
	}
    }
     
  double* rows_a = &rows_v[0];
  double* cols_a = &cols_v[0];
  double* data_a = &data_v[0];
  
  
  whitsMatrix->FillN(int(rows_v.size()), cols_a, rows_a, data_a);


  
}






void CSCSegAlgoUF::FillStripMatrix(TH2F* shitsMatrix, ChamberStripHitContainer shits) {

     std::vector<double> rows_v;
     std::vector<double> cols_v;
     std::vector<double> data_v;

     
     for (unsigned int i = 0; i < shits.size(); i++)
       {
         const CSCStripHit* shit = shits[i];
         int sLayer = shit->cscDetId().layer();


	 
         if (int(shit->strips().size()) == 1) 
	   {

            int sp = (shit->strips())[0];
            if (isME11 || (!isME11 && (sLayer == 2 || sLayer == 4 || sLayer == 6) ) )
	      {
		cols_v.push_back(2*sp-2); //cols_v.push_back(2*sp-1);
		//		std::cout << "layer: " << sLayer << " strip    "<< sp <<", half strip: " << 2*sp-2 << std::endl;
	      }
            if (!isME11 && (sLayer == 1 || sLayer == 3 || sLayer == 5) ) 
	      {
		
		cols_v.push_back(2*sp-1); //cols_v.push_back(2*sp);
		
	      }

            rows_v.push_back(sLayer-1);
            data_v.push_back(1);

	   }

	 

         if (int(shit->strips().size()) == 3)
	   {

            int sp = (shit->strips())[1];
	    //	    int sp0 = (shit->strips())[0];            // not used
	    //	    int sp2 = (shit->strips())[2];            // not used
	    
//            double leftC = (shit->s_adc())[1]; double rightC = (shit->s_adc())[9];
//	    shit->print();
//	    std::cout << (shit->s_adc())[0] << " " << (shit->s_adc())[4] << " " << (shit->s_adc())[8] << std::endl;
//	    std::cout << (shit->s_adc())[1] << " " << (shit->s_adc())[5] << " " << (shit->s_adc())[9] << std::endl;
//	    std::cout << (shit->s_adc())[2] << " " << (shit->s_adc())[6] << " " << (shit->s_adc())[10] << std::endl;
//	    std::cout << (shit->s_adc())[3] << " " << (shit->s_adc())[7] << " " << (shit->s_adc())[11] << std::endl;
	    
//	    std::cout<<"   s_adc  size     " <<     (shit->s_adc()).size() << "  maximum time bin    " <<shit->tmax() << "    sp 1,2,3     "
//		     << sp0 << "  " << sp <<"   " << sp2 << "    Layer      "<< sLayer  <<std::endl;

            double leftC  = (shit->s_adc())[0] + (shit->s_adc())[1] + (shit->s_adc())[2];
            double rightC = (shit->s_adc())[8] + (shit->s_adc())[9] + (shit->s_adc())[10];

	    
            // above line ugly hardcode, can be improved!
            // reason is for each strip hit, 4 time bins of adc for 3 strips in one strip cluster are saved: CSCHitFromStripOnly.cc

	    
            rows_v.push_back(sLayer-1);
	    data_v.push_back(1);
	    
            if (isME11 || (!isME11 && (sLayer == 2 || sLayer == 4 || sLayer == 6) ) )
	      {
//if (sLayer==3) cols_v.push_back(52);
/*else*/       leftC  >=  rightC ? cols_v.push_back(2*sp-2) : cols_v.push_back(2*sp-1);
		//std::cout << "layer: " << sLayer << ", half strip: " << (leftC >= rightC ? (2*sp-2) : (2*sp-1)) << std::endl;
		//std::cout << "layer: " << sLayer << ", leftC: " << leftC << ", half strip: " <<  (2*sp-2) << std::endl;
		//std::cout << "layer: " << sLayer << ", rightC: " << rightC << ", half strip: " <<  (2*sp-1) << std::endl;
               }
            if (!isME11 && (sLayer == 1 || sLayer == 3 || sLayer == 5) ) {
	      //	      std::cout<<" 2*sp - 1   "<< 2*sp-1 << "  2*sp   "<< 2*sp << std::endl;
               leftC >= rightC ? cols_v.push_back(2*sp-1) : cols_v.push_back(2*sp);
               } 

            }
     }
     
     double* rows_a = &rows_v[0];
     double* cols_a = &cols_v[0];
     double* data_a = &data_v[0];

     
     shitsMatrix->FillN(int(rows_v.size()), cols_a, rows_a, data_a);

}






void CSCSegAlgoUF::ScanForWireSeg(TH2F* wireHitsInChamber, std::list<CSCWireSegment>& wireSegs, int NumberOfLayers, bool debug) {


     TH1D* nonEmptyWG = wireHitsInChamber->ProjectionX();
     
     TH2F* wireHitsInChamber_AtStartingPoint  =  wireHitsInChamber;
     
     std::vector<TH2F*> SegmentsTH2F;

     
     for (int iWireGroup = 0; iWireGroup < nWireGroups; iWireGroup++)
       {
	 int thisKeyWG = nWireGroups - iWireGroup;   // scan from the wider chamber side


	 //	 if (nonEmptyWG->GetBinContent(thisKeyWG) == 0 && // only scan non empty area, what about under/over flow ? need to check !
	 //	     nonEmptyWG->GetBinContent(thisKeyWG+1) == 0 && nonEmptyWG->GetBinContent(thisKeyWG+2) == 0 ) continue;
	 

       for (int iPattern = 0; iPattern < nDefinedWirePatterns; iPattern++)  //  loop over 4 predefined wire patterns
	 {
	   
	   if(wireHitsInChamber->Integral()==0) continue;  // Don't proceed if either there are no hits or if no hits left after creating a segment

	   
	   if(iPattern!=0) continue; // remove me later

	   
	   //	   std::cout<<"   keyWG   "<< thisKeyWG << "   wire rank pattern   "<< iPattern  << "   nWireGroups   " << nWireGroups <<"  WG  loop counter      " << iWireGroup <<std::endl;
	   //	   std::cout<<"  wireHitsInChamber    total INtegral "<< wireHitsInChamber->Integral() << std::endl;

	   WriteTH2F(wireHitsInChamber);
	   //	   TH1D* hitsLayer = wireHitsInChamber->ProjectionY();
	   //	   wireHitsInChamber
	   
	   
	   // ========================== fill up the pattern matrix =================================
	   int nWGsInPattern = nWGsInPatterns[iPattern];

	   double w_cols_scan[nWGsInPattern] = {};
	   
	   for (int k = 0; k < nWGsInPattern; k++) // loop over pattern
	       {
		 w_cols_scan[k] = w_cols[iPattern][k] + thisKeyWG - 1;
	       }
	   double w_rows_scan[nWGsInPattern] = {};

	   
	   if (theStation == 3 || theStation == 4) //  invert the order for these chambers as layers are counted in opposite Z direction
	     for (int k = 0; k < nWGsInPattern; k++){w_rows_scan[k] = 5-w_rows[iPattern][k];}
	   else
	     for (int k = 0; k < nWGsInPattern; k++){ w_rows_scan[k] = w_rows[iPattern][k];}
	   
	   TH2F* wirePattern = new TH2F("wirePattern","", nWireGroups, 0, nWireGroups, 6, 0, 6);
	   wirePattern->FillN(nWGsInPattern, w_cols_scan, w_rows_scan, w_data[iPattern]);
	   

	   std::cout<<">>>>>>>>>>>>>>>>>> starting point    "<< std::endl; WriteTH2F(wireHitsInChamber_AtStartingPoint);
	   std::cout<<"|\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ <<<<<<<<<<<<<<<<<< pattern   "<< std::endl;  WriteTH2F(wirePattern);	   


	   


	   wirePattern->Multiply(wireHitsInChamber); // Multiply pattern and actual hits in the chamber
	   TH2F* actualSegment  = wirePattern;            // intersection between pattern and hits;  after multiplied wirePattern becomes an actual segment

	   
	   TH1D* hitsLayer = wirePattern->ProjectionY(); //<- this histo is Y projection, so it's layers with hits
	   hitsLayer->Divide(hitsLayer);                 //  normalised to itself
	   int thisNLayer = hitsLayer->Integral();       //  how many layers with hits 




	   

	   

	   if(actualSegment->Integral()!=0)
	     {
	       std::cout<<"================================== actual segment "<< std::endl;
	       WriteTH2F(actualSegment);
	     }



	   


	   
	   
	   if (thisNLayer !=  NumberOfLayers)  // break the pattern loop if N layers with hits != to requested: 6,5,4,3
	     {
	       
	       delete wirePattern;
	       continue;
	       
	     } 








	   

	   
	   

	   if (  int(wireSegs.size()) == 0 )  // if the first segment
	       {

		 wireSegs.push_back( CSCWireSegment(thisKeyWG, thisNLayer, wirePattern) ); 
		 SegmentsTH2F.push_back(actualSegment);

		 wireHitsInChamber->Add( wirePattern, -1 );  // remove "segment" from chamber hits
		 
		 delete wirePattern;
		 continue;  // continue loop over patterns
		 
	       }
	 
	   
             CSCWireSegment lastSegment                         =  wireSegs.back();
             int            lastKeyWG                           =  lastSegment.keyWG();
	     
	     CSCWireSegment potentialAnotherSegment             = CSCWireSegment(thisKeyWG, thisNLayer, wirePattern);
	     






	     

             if (abs(lastKeyWG - thisKeyWG) >  1)
	       {

		 wireSegs.push_back(potentialAnotherSegment); //  if the WG difference between segment > 1 create the new segment
		 SegmentsTH2F.push_back(actualSegment);

	       }








	     
             if (abs(lastKeyWG - thisKeyWG) == 1)  // combine two segments if the "keyWG difference is 1"; In text it's misleading;
	       {

		 //		 std::cout<<" Origianla    Before UpdateHits    N Layers Hits   "<< (wireSegs.back().nLayerHits()) << "   nWire Hits    " << (wireSegs.back().wireHitsPosition())  <<std::endl;
		 //		 std::cout<<"   nHits Before UpdateHits    N Layers Hits   "<< (*potentialAnotherSegment.nLayerHits()) << "   nWire Hits    " << (potentialAnotherSegment.wireHitsPosition())  <<std::endl;
		 wireSegs.back().updateWHits( potentialAnotherSegment.wireHitsPosition(), potentialAnotherSegment.nLayerHits());
		 SegmentsTH2F.push_back(actualSegment);
		 //	 std::cout<<"   nHits After  UpdateHits    N Layers Hits   "<< (wireSegs.back().nLayerHits()) << "   nWire Hits   in layer 1  " << (wireSegs.back().wireHitsPosition()[1])  <<std::endl;

	       }
	     //	     std::cout<<"   -------- lastKeyWG    "<< lastKeyWG <<"  thisKeyWG    "  << thisKeyWG <<"    thisNlayer   " << thisNLayer <<"  and wire segment size   " << wireSegs.size()  <<std::endl;
	     //	     std::cout<<"  Let's check the starting point     "<< std::endl; WriteTH2F(wireHitsInChamberCopy_removeMeLater);std::cout << std::endl;







	     
             wireHitsInChamber->Add(wirePattern, -1); // delete wire group being used
	     //	     std::cout << "after:" << std::endl; WriteTH2F(wireHitsInChamber);
	     //	     std::cout<<" actual Segment:   "<< std::endl; WriteTH2F(actualSegment);

	     //	     potentialAnotherSegment.printWireSegment();
             delete wirePattern;
	     
	     std::cout<<"    end of pattern loop   " << std::endl;
             } // pattern loop
       std::cout<<"    end of WG loop    last WG:  " << iWireGroup << std::endl;
       } // WG loop
     //     std::cout<<"====================================  wieSegsSize at the end  layer  "  << NumberOfLayers  << "  size      "<<     wireSegs.size() << "  vec size  "  << SegmentsTH2F.size() <<std::endl;

     
     //     if(SegmentsTH2F.size()!=0) std::cout<<"  bins check "<<  SegmentsTH2F.at(0)->GetNbinsX() << std::endl;
     //     for(auto Wiresegment : wireSegs)
       //       Wiresegment.printWireSegment();
       //     std::cout<<"------------------------------------------------"<< std::endl;
     //     for(auto segmentHist : SegmentsTH2F)
     //       {
	 //	 std::cout<<" ??? "<< std::endl;
	 //	 std::cout<< "  binx   "<< segmentHist->GetNbinsX() << "   biny  " << segmentHist->GetNbinsY() <<std::endl;
     //	 WriteTH2F(segmentHist);
     //       }
}







void CSCSegAlgoUF::ScanForStripSeg(TH2F* sHitsPerChamber, std::list<CSCStripSegment>& stripSegs, int nLayer) {

     TH1D* nonEmptyStrip = sHitsPerChamber->ProjectionX();
     
     for (int iStrip = 0; iStrip < nStrips*2+1; iStrip++) {

       
         int thisKeyHalfStrip = iStrip+1;

         if (nonEmptyStrip->GetBinContent(iStrip+1) == 0 && // only scan non empty area, what about under/over flow ? need to check !
             nonEmptyStrip->GetBinContent(iStrip) == 0 && nonEmptyStrip->GetBinContent(iStrip+2) == 0 ) continue;


         for (int iStripPattern = 0; iStripPattern < nDefinedStripPatterns;   iStripPattern++   )
	   {
	     // scan patterns sequentially, according to ranks
	     // scan from high rank (high pt) patterns to low pt patterns
	     if(sHitsPerChamber->Integral()==0) continue;  // dont continue if there are no hits or no hits left when strip segment is formed
	     
             int nhalfStrips = nHalfStrips[iStripPattern]; 
             int thisRank = patternRanks[iStripPattern];
             double s_cols_scan[nhalfStrips] = {};

	     
	     
             if (isME11)
	       {
		 for (int k = 0; k < nhalfStrips; k++) s_cols_scan[k] = s_cols[iStripPattern][k] + thisKeyHalfStrip -1;
	       }
	     else
	       {
		 for (int k = 0; k < nhalfStrips; k++) s_cols_scan[k] = s_cols[iStripPattern][k] + thisKeyHalfStrip ;
	       }



	     double s_rows_scan[nhalfStrips] = {};
             if (theStation == 3 || theStation == 4)   // invert layers numbering for these stations
                for (int k = 0; k < nhalfStrips; k++) s_rows_scan[k] = 5-s_rows[iStripPattern][k];
             else
                for (int k = 0; k < nhalfStrips; k++) s_rows_scan[k] = s_rows[iStripPattern][k];


	     
             TH2F* stripPattern = new TH2F("stripPattern","", nStrips*2+1, 0, nStrips*2+1, 6, 0, 6);

             stripPattern->FillN(nhalfStrips,s_cols_scan,s_rows_scan,/*s_rows[j]*/s_data[iStripPattern]); 

             stripPattern->Multiply(sHitsPerChamber); // scan is done here

             TH1D* hitsLayer = stripPattern->ProjectionY();
             hitsLayer->Divide(hitsLayer); 
             int thisNLayer = hitsLayer->Integral(); // find layers with s hits

             if (thisNLayer != nLayer)  // 3 or 4 is hardcode

	       {
		 
		 delete stripPattern;
		 continue;	       // here continue loop over patterns
		 
	       }


	     
             if (int(stripSegs.size()) == 0 )
	       {

		 stripSegs.push_back( CSCStripSegment(thisKeyHalfStrip, thisNLayer, thisRank, stripPattern) ); 
		 sHitsPerChamber->Add(stripPattern,-1);
		 
		 delete stripPattern;
		 continue;
	       }

             CSCStripSegment lastSegment = stripSegs.back();
             int lastKeyHalfStrip = lastSegment.keyHalfStrip();
	     //std::cout << "updated sPattern: " << std::endl; WriteTH2F(stripPattern);
	     //std::cout << std::endl;

             CSCStripSegment tmpStripSeg = CSCStripSegment(thisKeyHalfStrip,thisNLayer,thisRank,stripPattern);
	     
             if (abs(thisKeyHalfStrip - lastKeyHalfStrip) > 1) stripSegs.push_back(tmpStripSeg);
	     
             if (abs(thisKeyHalfStrip - lastKeyHalfStrip) == 1)
	       {
		 
		 stripSegs.back().updateSHits(tmpStripSeg.stripHits(), tmpStripSeg.nLayerHits());

		 
	       }

             sHitsPerChamber->Add(stripPattern,-1); 
             delete stripPattern;


             }
         }
}



void CSCSegAlgoUF::GetWireHitFromWireSeg(CSCWireSegment wireSeg, ChamberWireHitContainer WireHitsInChamber, int* wireHitIndex) {


  //  std::cout<<"||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||    GetWireHitFromWireSeg  ||||||||||||||||  "   << std::endl;
  double lowerWireGroup  = wireSeg.LowestHitInLayer();
  double higherWireGroup = wireSeg.HighestHitInLayer();
  
//std::cout << "keyWH: " << keyWH << std::endl;
//  std::cout << " ============================  low: " << lowerWireGroup << ", high: " << higherWireGroup << std::endl;
  //  wireSeg.printWireSegment();

  
  std::vector<int> addBackCounter; // a layer w/o hits

  for (int i = 0; i < 6; i++)     // loop over 6 layers
    { 
      double wireHitPos_from_segment = (wireSeg.wireHitsPosition())[i];
      //      std::cout<<"  ----------------- wHit Pos       "<< wireHitPos_from_segment << "  layer      "<< i <<std::endl;
      int wireHitIndex_from_hits_collection = -999;
      
      double wireHitPosDelta = 113; //112 is max one can get
      //std::cout << "wpos: " << wireHitPos_from_segment << std::endl;
      bool hitMissed = false; // ??????????
      
      //    std::cout<<"||||||||||||||  N wire hits   "<< int(WireHitsInChamber.size()) << std::endl;
      for (int hit_index = 0; hit_index < int(WireHitsInChamber.size()); hit_index++) // loop over all hits 
	{
	  
          const CSCWireHit* tmpHit = WireHitsInChamber[hit_index]; //  this is the Chamber Wire Hit
          int wireHitLayer         = tmpHit->cscDetId().layer();
          double wireHitPos_from_hits_collection = tmpHit->wHitPos();

	  //if (i==0) std::cout << "wHitPos2: " << wireHitPos_from_hits_collection << std::endl;
          if (wireHitLayer != (i+1) ) continue;  // skip if not in this layer
	  //  std::cout<<"   layer     "<< i << "   wire hit position from segment    " <<wireHitPos_from_segment <<"   hit position from collection    " << wireHitPos_from_hits_collection <<std::endl;	  
          if (wireHitPos_from_segment == 0 && /*abs(keyWH-  wireHitPos_from_hits_collection) <= 2 &&*/ (wireHitPos_from_hits_collection >= lowerWireGroup-1) && (wireHitPos_from_hits_collection <= higherWireGroup+1) )
	    {
	      //
	      // std::cout << "here:   whits Pos  " << wireHitPos_from_hits_collection <<"   wireSeg Wire Posi     " << wHitPos << "   dif to keyWH     " <<  abs(keyWH-wHitPos2) <<" lo   " <<lowerWireGroup <<" hi  " << higherWireGroup <<std::endl;
	      //	      std::cout<<" -------------------------------------------------------------   add backcounter " << i << std::endl;
	      //	      std::cout<<"  ever here ??? "  << std::endl;
	      addBackCounter.push_back(i);   // store layer number where there no hit from Wire segment  and no wire hit
	      wireHitIndex_from_hits_collection = hit_index;
	      hitMissed = true;
	      break;
	     
	    }




          if (abs( wireHitPos_from_segment  -  wireHitPos_from_hits_collection ) < wireHitPosDelta)
	    {
	      wireHitPosDelta = abs(wireHitPos_from_segment  -  wireHitPos_from_hits_collection);
	      wireHitIndex_from_hits_collection = hit_index; // find the closest hit 
	    }



	  
	} // loop over hits collection
      //      std::cout<<"  delta   " << wireHitPosDelta << "   backcounetr  size  "<< addBackCounter.size() <<std::endl;


      
      if (  (  wireHitPosDelta < 113 && wireHitPos_from_segment >= 1  )    ||    ( wireHitPos_from_segment ==0 && hitMissed )   )
	// if pattern doesn't cover some hit, and is close, to keyWH, include it // <-- It's not true!
	{
	  
	  wireHitIndex[i] = wireHitIndex_from_hits_collection;
	  if(hitMissed)
	    std::cout<<"  Missing hits =================            " <<wireHitIndex[i] << "  position " << WireHitsInChamber[wireHitIndex_from_hits_collection]->wHitPos()<< "  and delta   " << wireHitPosDelta  <<std::endl;
	  
	}
      else
	{
	  
	  wireHitIndex[i] = -1;
	  
	}
      
    }// loop over layers


  
  if (int(addBackCounter.size())>0) std::cout << "Wire add back " << addBackCounter.size() << " times" << std::endl;



  
  if (wireSeg.nLayersWithHits() == 3 && int(addBackCounter.size()) > 0) std::cout << "add back recover ALCT" << addBackCounter.size() << " times" << std::endl;

}




void CSCSegAlgoUF::GetStripHitFromStripSeg(CSCStripSegment stripSeg, ChamberStripHitContainer shits, int* stripHitIndex)
{

//  double keySH = ceil(stripSeg.keyHalfStrip()/2.0);
//std::cout << "keySH: " << keySH << std::endl;
  double LowestHitInLayer  = stripSeg.LowestHitInLayer(isME11);
  double HighestHitInLayer = stripSeg.HighestHitInLayer(isME11);
//std::cout << "low: " << LowestHitInLayer << ", high: " << HighestHitInLayer << std::endl;

  std::vector<int> addBackCounter;

  for (int i = 0; i < 6; i++) {
      double sHitPos_ = (stripSeg.stripHits())[i]; // comparator position in unit of half strip
      double sHitPos = -1;                         // strip hit position in unit of strip
 

      
      if (!isME11 && (i==0 || i==2 || i==4) )              sHitPos = ceil((sHitPos_-1)/2.0); // convert comparator number to strip number
      if (isME11 || (!isME11 && (i==1 || i==3 || i==5) ) ) sHitPos = ceil(sHitPos_/2.0);
      if (sHitPos_==1 || sHitPos_==2)                      sHitPos = 1;

      
      int sHitIndex = -999;
      double sPosDiff = 81;
      bool hitMissed = false;
      
//std::cout << "sHitPos: " << sHitPos << std::endl;
      for (int j = 0; j < int(shits.size()); j++) {
          const CSCStripHit* tmpHit = shits[j];
          int sHitLayer = tmpHit->cscDetId().layer();
          double sHitPos2 = tmpHit->sHitPos();
//if (i==0) std::cout << "sHitPos2: " << sHitPos2 << std::endl;

          if ((i+1) != sHitLayer) continue;

          if (sHitPos == 0 /*&& abs(keySH-sHitPos2)<2*/ && ((2*sHitPos2) >= LowestHitInLayer-2) && (2*sHitPos2-1 <= HighestHitInLayer+2) ) {
             addBackCounter.push_back(i);
             sHitIndex = j; hitMissed = true; break;
             }

          if (abs(sHitPos-sHitPos2) < sPosDiff) {
             sPosDiff = abs(sHitPos-sHitPos2);
             sHitIndex = j;

             }
          }


      if ((sPosDiff < 81 && sHitPos >= 1) || (sHitPos==0&&hitMissed) ) {  // this line needs some simplification ! this is really ugly ...
         stripHitIndex[i] = sHitIndex;
//std::cout << "layer:" << i+1 << ", sHitIndex: " << sHitIndex << std::endl;

         } else {
                stripHitIndex[i] = -1;
                }

      }

      if (int(addBackCounter.size())>0) std::cout << "Strip add back " << addBackCounter.size() << " times" << std::endl;
      if (int(addBackCounter.size())>0 && stripSeg.nLayersWithHits() == 5 && stripSeg.patternRank() == 1) std::cout << "5 hits rank1 strip add back " << addBackCounter.size() << " times" << std::endl;
      if (int(addBackCounter.size())>0 && stripSeg.nLayersWithHits() == 5 ) std::cout << "5 hits all ranks strip add back " << addBackCounter.size() << " times" << std::endl;
      if (stripSeg.nLayersWithHits()==3 && int(addBackCounter.size())>0) std::cout << "add back recover CLCT" << addBackCounter.size() << " times" << std::endl;
      
}








std::vector<CSCSegment> CSCSegAlgoUF::prune_bad_hits(const CSCChamber* Chamber, std::vector<CSCSegment> & segments) {
  
//     std::cout<<"*************************************************************"<<std::endl;
//     std::cout<<"Called prune_bad_hits in Chamber "<< theChamber->specs()->chamberTypeName()<<std::endl;
//     std::cout<<"*************************************************************"<<std::endl;
 
 
  std::vector<CSCSegment>          segments_temp;
  std::vector<ChamberHitContainer> rechits_clusters; // this is a collection of groups of rechits
  
  const float chi2ndfProbMin = 1.0e-4;
  bool   use_brute_force = true;//BrutePruning; hm

  int hit_nr = 0;
  int hit_nr_worst = -1;
  //int hit_nr_2ndworst = -1;
  
  for(std::vector<CSCSegment>::iterator it=segments.begin(); it != segments.end(); ++it) {
    
    // do nothing for nhit <= minHitPerSegment
    if( (*it).nRecHits() <= 3/*minHitsPerSegment hm*/) continue;
    
    if( !use_brute_force ) {// find worst hit
      
      float chisq    = (*it).chi2();
      int nhits      = (*it).nRecHits();
      LocalPoint localPos = (*it).localPosition();
      LocalVector segDir = (*it).localDirection();
      const CSCChamber* cscchamber = theChamber;
      float globZ;
	  
      GlobalPoint globalPosition = cscchamber->toGlobal(localPos);
      globZ = globalPosition.z();
      
      
      if( ChiSquaredProbability((double)chisq,(double)(2*nhits-4)) < chi2ndfProbMin  ) {

	// find (rough) "residuals" (NOT excluding the hit from the fit - speed!) of hits on segment
	std::vector<CSCRecHit2D> theseRecHits = (*it).specificRecHits();
	std::vector<CSCRecHit2D>::const_iterator iRH_worst;
	//float xdist_local       = -99999.;

	float xdist_local_worst_sig = -99999.;
	float xdist_local_2ndworst_sig = -99999.;
	float xdist_local_sig       = -99999.;

	hit_nr = 0;
	hit_nr_worst = -1;
	//hit_nr_2ndworst = -1;

	for ( std::vector<CSCRecHit2D>::const_iterator iRH = theseRecHits.begin(); iRH != theseRecHits.end(); ++iRH) {
	  //mark "worst" hit:
	  
 	  //float z_at_target ;
	  //float radius      ;
	  float loc_x_at_target;
	  //float loc_y_at_target;
	  //float loc_z_at_target;

	  //z_at_target  = 0.;
	  //radius       = 0.;
	  
	  // set the z target in CMS global coordinates:
	  const CSCLayer* csclayerRH = theChamber->layer((*iRH).cscDetId().layer());
	  LocalPoint localPositionRH = (*iRH).localPosition();
	  GlobalPoint globalPositionRH = csclayerRH->toGlobal(localPositionRH);	
	  
	  LocalError rerrlocal = (*iRH).localPositionError();  
	  float xxerr = rerrlocal.xx();
	  
	  float target_z     = globalPositionRH.z();  // target z position in cm (z pos of the hit)
	  
	  if(target_z > 0.) {
	    loc_x_at_target = localPos.x() + (segDir.x()/fabs(segDir.z())*( target_z - globZ ));
	    //loc_y_at_target = localPos.y() + (segDir.y()/fabs(segDir.z())*( target_z - globZ ));
	    //loc_z_at_target = target_z;
	  }
	  else {
	    loc_x_at_target = localPos.x() + ((-1)*segDir.x()/fabs(segDir.z())*( target_z - globZ ));
	    //loc_y_at_target = localPos.y() + ((-1)*segDir.y()/fabs(segDir.z())*( target_z - globZ ));
	    //loc_z_at_target = target_z;
	  }
	  // have to transform the segments coordinates back to the local frame... how?!!!!!!!!!!!!
	  
	  //xdist_local  = fabs(localPositionRH.x() - loc_x_at_target);
	  xdist_local_sig  = fabs((localPositionRH.x()  -  loc_x_at_target)/(xxerr));
	  
	  if( xdist_local_sig > xdist_local_worst_sig ) {
	    xdist_local_2ndworst_sig = xdist_local_worst_sig;
	    xdist_local_worst_sig    = xdist_local_sig;
	    iRH_worst            = iRH;
	    //hit_nr_2ndworst = hit_nr_worst;
	    hit_nr_worst = hit_nr;
	  }
	  else if(xdist_local_sig > xdist_local_2ndworst_sig) {
	    xdist_local_2ndworst_sig = xdist_local_sig;
	    //hit_nr_2ndworst = hit_nr;
	  }
	  ++hit_nr;
	}

	// reset worst hit number if certain criteria apply.
	// Criteria: 2nd worst hit must be at least a factor of
	// 1.5 better than the worst in terms of sigma:
	if ( xdist_local_worst_sig / xdist_local_2ndworst_sig < 1.5 ) {
	  hit_nr_worst    = -1;
	  //hit_nr_2ndworst = -1;
	}
      }
    }

    // if worst hit was found, refit without worst hit and select if considerably better than original fit.
    // Can also use brute force: refit all n-1 hit segments and choose one over the n hit if considerably "better"
   
    std::vector< CSCRecHit2D > buffer;
    std::vector< std::vector< CSCRecHit2D > > reduced_segments;
    std::vector< CSCRecHit2D > theseRecHits = (*it).specificRecHits();
    float best_red_seg_prob = 0.0;
    // usefor chi2 1 diff   float best_red_seg_prob = 99999.;
    buffer.clear();

    if( ChiSquaredProbability((double)(*it).chi2(),(double)((2*(*it).nRecHits())-4)) < chi2ndfProbMin  ) {
	
      buffer = theseRecHits;

      // Dirty switch: here one can select to refit all possible subsets or just the one without the 
      // tagged worst hit:
      if( use_brute_force ) { // Brute force method: loop over all possible segments:
	for(size_t bi = 0; bi < buffer.size(); ++bi) {
	  reduced_segments.push_back(buffer);
	  reduced_segments[bi].erase(reduced_segments[bi].begin()+(bi),reduced_segments[bi].begin()+(bi+1));
	}
      }
      else { // More elegant but still biased: erase only worst hit
	// Comment: There is not a very strong correlation of the worst hit with the one that one should remove... 
	if( hit_nr_worst >= 0 && hit_nr_worst <= int(buffer.size())  ) {
	  // fill segment in buffer, delete worst hit
	  buffer.erase(buffer.begin()+(hit_nr_worst),buffer.begin()+(hit_nr_worst+1));
	  reduced_segments.push_back(buffer);
	}
	else {
	  // only fill segment in array, do not delete anything
	  reduced_segments.push_back(buffer);
	}
      }
    }
      
    // Loop over the subsegments and fit (only one segment if "use_brute_force" is false):
    for(size_t iSegment=0; iSegment<reduced_segments.size(); ++iSegment) {
      // loop over hits on given segment and push pointers to hits into protosegment
      protoSegment.clear();
      for(size_t m = 0; m<reduced_segments[iSegment].size(); ++m ) {
	protoSegment.push_back(&reduced_segments[iSegment][m]);
      }

      // Create fitter object
      CSCCondSegFit* segfit = new CSCCondSegFit( pset(), theChamber/*chamber() hm*/, protoSegment );
      condpass1 = false;
      condpass2 = false;
      segfit->setScaleXError( 1.0 );
      segfit->fit(condpass1, condpass2);

      // Attempt to handle numerical instability of the fit;
      // The same as in the build method;
      // Preprune is not applied;
      if( true/*adjustCovariance() hm*/){
	if(segfit->chi2()/segfit->ndof()>10/*chi2Norm_3D_*/){
	  condpass1 = true;
	  segfit->fit(condpass1, condpass2);
	}
	if( (segfit->scaleXError() < 1.00005)&&(segfit->chi2()/segfit->ndof()>10/*chi2Norm_3D_*/) ){
	  condpass2 = true;
	  segfit->fit(condpass1, condpass2);
	}
      }
    
      // calculate error matrix 
      //      AlgebraicSymMatrix temp2 = segfit->covarianceMatrix();

      // build an actual segment
      CSCSegment temp(protoSegment, segfit->intercept(), segfit->localdir(), 
                         segfit->covarianceMatrix(), segfit->chi2() );

      // and finished with this fit
      delete segfit;

      // n-hit segment is (*it)
      // (n-1)-hit segment is temp
      // replace n-hit segment with (n-1)-hit segment if segment probability is BPMinImprovement better
      double oldchi2 = (*it).chi2();
      double olddof  =  2 * (*it).nRecHits() - 4;
      double newchi2 = temp.chi2();
      double newdof  = 2 * temp.nRecHits() - 4;
      if( ( ChiSquaredProbability(oldchi2,olddof) < (1./10000.0/*BPMinImprovement hm*/)*
	    ChiSquaredProbability(newchi2,newdof) ) 
	  && ( ChiSquaredProbability(newchi2,newdof) > best_red_seg_prob )
	  && ( ChiSquaredProbability(newchi2,newdof) > 1e-10 )
	) {
	best_red_seg_prob = ChiSquaredProbability(newchi2,newdof);
        // The (n-1)- segment is better than the n-hit segment.
	// If it has at least  minHitsPerSegment  hits replace the n-hit segment
        // with this  better (n-1)-hit segment:
        if( temp.nRecHits() >= 3/*minHitsPerSegment hm*/) {
          (*it) = temp;
        }
      }
    } // end of loop over subsegments (iSegment)
    
  } // end loop over segments (it)
  
  return segments;
  
}






