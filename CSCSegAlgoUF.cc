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
  : CSCSegmentAlgorithm(ps), ps_(ps), myName("CSCSegAlgoUF")
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
  
  std::cout << theChamber->id() << std::endl;

  
  //  get number of wire groups and strips for this chamber
  nWireGroups = theChamber->layer(1)->geometry()->numberOfWireGroups();
  nStrips     = theChamber->layer(1)->geometry()->numberOfStrips();
  

  MaxStripNumber = 81;
  MaxWGNumber    = 113;


  theStation = theChamber->id().station();
  theRing    = theChamber->id().ring();


  
  if (theStation==1 && (theRing == 1 || theRing == 4) ) 
    {
      isME11 = true;
    }
  else
    {
      isME11 = false;
    }

  
  return buildSegments(wirehits, striphits);
}




std::vector<CSCSegment> CSCSegAlgoUF::buildSegments(const ChamberWireHitContainer&  uwirehits,  const ChamberStripHitContainer& ustriphits)
{
  std::cout<<"|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"<<std::endl;
  std::cout<<"|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"<<std::endl;
  std::cout<<"|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"<<std::endl;
  
  std::vector<CSCSegment> segments;

  //  int e = theChamber->id().endcap();
  //  int s = theChamber->id().station();
  //  int r = theChamber->id().ring();
  //  int c = theChamber->id().chamber();
  //  if (!(e == 1 && s == 2 && r == 1 && c == 16)) return segments;



  
  ChamberWireHitContainer  wirehits   = uwirehits;   //  input wire hits
  ChamberStripHitContainer striphits  = ustriphits;  //  input strip hits

  
  TH2F* wireHitsInChamber  = new TH2F("wireHitsInChamber",  "", nWireGroups, 0, nWireGroups, 6 ,0, 6);
  TH2F* stripHitsInChamber = new TH2F("stripHitsInChamber", "", 2*nStrips+1, 0, 2*nStrips+1, 6 ,0, 6); // half strip stagger at 1,3,5 layer

  
  FillWireMatrix(wireHitsInChamber,  wirehits);      TH2F* wireHitsInChamber_clone  = (TH2F*)wireHitsInChamber->Clone("wireHitsInChamber_clone");
  FillStripMatrix(stripHitsInChamber, striphits);    TH2F* stripHitsInChamber_clone = (TH2F*)stripHitsInChamber->Clone("stripHitsInChamber_clone");


  

  
  std::list<CSCWireSegment> wireSegments;
  std::vector<TH2F*>        wireSegmentsTH2F;
  std::vector<int>          wireSegments_rank; //

  
  //  Call Scans for segment 4 times, requesting 6,5,4 and 3 layers; wireSegments, wireSegmentsTH2F, wireSegments_rank
  //  are updated after every call; wireHitsInChamber -> actual wire hits are updated as well, hits that have already
  //  been used are removed sequentially. 
  ScanForWireSegment(wireHitsInChamber, wireSegments, wireSegmentsTH2F, wireSegments_rank, 6);
  ScanForWireSegment(wireHitsInChamber, wireSegments, wireSegmentsTH2F, wireSegments_rank, 5);
  ScanForWireSegment(wireHitsInChamber, wireSegments, wireSegmentsTH2F, wireSegments_rank, 4);
  ScanForWireSegment(wireHitsInChamber, wireSegments, wireSegmentsTH2F, wireSegments_rank, 3);



  
  std::list<CSCStripSegment> stripSegments;
  std::vector<TH2F*>         stripSegmentsTH2F;
  std::vector<int>           stripSegments_rank;

  //  Call scans for segment 4 times, requesting 6,5,4 and 3 layers; stripSegments, stripSegmentsTH2F, stripSegments_rank
  //  are updated after every call; wireHitsInChamber -> actual wire hits are updated as well, hits that have already
  //  been used are removed sequentially.
  
  ScanForStripSegment(stripHitsInChamber, stripSegments, stripSegmentsTH2F, stripSegments_rank,  6);
  ScanForStripSegment(stripHitsInChamber, stripSegments, stripSegmentsTH2F, stripSegments_rank,  5);
  ScanForStripSegment(stripHitsInChamber, stripSegments, stripSegmentsTH2F, stripSegments_rank,  4);
  ScanForStripSegment(stripHitsInChamber, stripSegments, stripSegmentsTH2F, stripSegments_rank,  3);







  /*
  
  // print out segment for debugging purposes
  std::cout<<"  Wire  Hits in Chamber:  " << std::endl;
  PrintTH2F(wireHitsInChamber_clone);

  std::cout<<"  N built segments: "<< wireSegmentsTH2F.size() << std::endl;
  for(auto s : wireSegments)
    {
      s.printWireSegment();  std::cout<<std::endl;
    }
  std::cout<<" Built Wire Segments:   "<< std::endl;
  for(auto sH : wireSegmentsTH2F)
    {
      PrintTH2F(sH); std::cout<<std::endl;
    }
  for(auto sR : wireSegments_rank)
    {
      std::cout<<"strip  segment rank  "<< sR << std::endl;
    }

  */
  
  
  
  // print out segment for debugging purposes
  std::cout<<"  Strip  Hits in Chamber:  " << std::endl;
  PrintTH2F(stripHitsInChamber_clone);
  for(auto s : stripSegments)
    {
      s.printStripSegment();  std::cout<<std::endl;
    }
  std::cout<<" Built Strip Segments:   "<< std::endl;
  for(auto sH : stripSegmentsTH2F)
    {
      PrintTH2F(sH); std::cout<<std::endl;
    }
  for(auto sR : stripSegments_rank)
    {
      std::cout<<" segment rank  "<< sR << std::endl;
    }

  



  
  
  std::cout << theChamber->id() <<  ", nWireSegments: " << wireSegments.size() << 
    ", nStripSegments: " << stripSegments.size() << 
    ", nSegments: " << wireSegments.size() * stripSegments.size() << std::endl;


  

  //  first loop over wireSegments
  
  for (auto i_wire = wireSegments.begin(); i_wire != wireSegments.end(); i_wire++)
    {


      int wireHitsFromWireSegment[6]   = {};        GetWireHitFromWireSegment(*i_wire,    wirehits,  wireHitsFromWireSegment);     // misleading naming; Those are matched hits, not actual hits from segment found by patterns
      CSCWireSegment s = *i_wire;
      s.printWireSegment();
      
      std::cout<<"  print out ========   "<<wireHitsFromWireSegment[0] << std::endl;
      std::cout<<"  print out ========   "<<wireHitsFromWireSegment[1] << std::endl;
      
       ///////////////////////////////// Fill TH2F to dispaly; to be commented out 
       ChamberWireHitContainer FinalSegmentWireHits;
       for(int iLayer =0; iLayer < 6; iLayer++)
	 {
	   if(wireHitsFromWireSegment[iLayer]!=-1)
	     {
	       FinalSegmentWireHits.push_back(wirehits[wireHitsFromWireSegment[iLayer]]);
	     }
	 }
       TH2F* SegmentwireHitsInChamber  = new TH2F("SegmentwireHitsInChamber",  "", nWireGroups, 0, nWireGroups, 6 ,0, 6);
       FillWireMatrix(SegmentwireHitsInChamber,  FinalSegmentWireHits);
       //       std::cout<<"  FINAL WIRE  SEGMENT     "<< std::endl;PrintTH2F(SegmentwireHitsInChamber);
       ///////////////////////////////////////////////////////////////////////////


       
       

	
       //  loop over strip segments
       for (auto i_strip = stripSegments.begin(); i_strip != stripSegments.end(); i_strip++)
	 {

	   int stripHitsFromStripSegment[6] = {};        GetStripHitFromStripSegment(*i_strip, striphits, stripHitsFromStripSegment);   // misleading naming; Those are matched hits, not actual hits from segment found by patterns 

	   

	   ///////////////////////////////// Fill TH2F to dispaly; to be commented out 
	   ChamberStripHitContainer FinalSegmentStripsHits;
	   for(int iLayer =0; iLayer < 6; iLayer++)
	     {
	       if(stripHitsFromStripSegment[iLayer]!=-1)
		 {
		   FinalSegmentStripsHits.push_back(striphits[stripHitsFromStripSegment[iLayer]]);
		 }
	     }
	   TH2F* SegmentstripHitsInChamber = new TH2F("SegmentstripHitsInChamber", "", 2*nStrips+1, 0, 2*nStrips+1, 6 ,0, 6); // half strip stagger at 1,3,5 layer
	   FillStripMatrix(SegmentstripHitsInChamber, FinalSegmentStripsHits);    
	   std::cout<<"  FINAL STRIP  SEGMENT     "<< std::endl;PrintTH2F(SegmentstripHitsInChamber);
	   ///////////////////////////////////////////////////////////////////////////
	


	
	   std::vector<CSCRecHit2D> csc2DRecHits;

	
	   for (unsigned int iLayer = 0; iLayer < 6; iLayer++) // loop over layers
	     {
	       if (wireHitsFromWireSegment[iLayer] == -1 || stripHitsFromStripSegment[iLayer] == -1) continue;
	       
	       const CSCWireHit*  cscwirehit   = wirehits[wireHitsFromWireSegment[iLayer]];  // 
	       const CSCStripHit* cscstriphit  = striphits[stripHitsFromStripSegment[iLayer]]; // 
	       
	       const CSCWireHit&  wirehit  = *cscwirehit;
	       const CSCStripHit& striphit = *cscstriphit;
	       
	       const CSCDetId&    detId    = CSCDetId(theChamber->id().endcap(),
						      theChamber->id().station(),
						      theChamber->id().ring(),
						      theChamber->id().chamber(), iLayer+1);
	       const CSCLayer* cscLayer = theChamber->layer(iLayer+1);
	       
	       CSCRecHit2D rechit = make2DHits_->hitFromStripAndWire(detId, cscLayer, wirehit, striphit );   //to be reviewed

	       
	       //	       std::cout << csc2DRecHits.size() << " RHs before  fiducial cut" << std::endl;
	       if( make2DHits_->isHitInFiducial( cscLayer, rechit )  )
		 csc2DRecHits.push_back(rechit);
	       
	     }
	   



	
	   //	   std::cout << csc2DRecHits.size() << " RHs survived fiducial cut" << std::endl;
	   ChamberHitContainer csc2DRecHits_input_to_build_segment;                                         // this action here is only to reconile formats
	   for (std::vector<CSCRecHit2D>::const_iterator it = csc2DRecHits.begin(); it != csc2DRecHits.end(); it++)
	     csc2DRecHits_input_to_build_segment.push_back(&(*it));


	   
	   if ( int(csc2DRecHits_input_to_build_segment.size() ) < 3 ) continue; // proceed only if >=3 rechits
	  




	  



	  
	   // borrow from ST  //  -------------------------------------- to debug .... 
	   
	   CSCCondSegFit* segment_fit = new CSCCondSegFit( pset(), theChamber, csc2DRecHits_input_to_build_segment );
	   condpass1 = false;
	   condpass2 = false;
	   
	   segment_fit->setScaleXError( 1.0 );
	   segment_fit->fit(condpass1, condpass2);
	   
	   if(segment_fit->chi2()/segment_fit->ndof() > chi2Norm_3D_)
	     {
	       condpass1 = true;
	       segment_fit->fit(condpass1, condpass2);
	     }
	   
	   if(segment_fit->scaleXError() < 1.00005)
	     {
	       LogTrace("CSCWeirdSegment") << "[CSCSegAlgoST::buildSegments] Segment ErrXX scaled and refit " << std::endl;
	       if(segment_fit->chi2()/segment_fit->ndof()>chi2Norm_3D_)
		 {
		   LogTrace("CSCWeirdSegment") << "[CSCSegAlgoST::buildSegments] Segment ErrXY changed to match cond. number and refit " << std::endl;
		   condpass2 = true;
		   segment_fit->fit(condpass1, condpass2);
		 }
	     }
	   
	   
	   if(prePrun_ && (sqrt(segment_fit->scaleXError())>prePrunLimit_) &&   (segment_fit->nhits()>3))
	     {
	       
	       LogTrace("CSCWeirdSegment") << "[CSCSegAlgoST::buildSegments] Scale factor chi2uCorrection too big, pre-Prune and refit " << std::endl;
	       
	       csc2DRecHits_input_to_build_segment.erase(csc2DRecHits_input_to_build_segment.begin() + segment_fit->worstHit(), csc2DRecHits_input_to_build_segment.begin() + segment_fit->worstHit()+1 );
	       
	       double tempcorr = segment_fit->scaleXError(); // save current value
	       delete segment_fit;    
	       segment_fit = new CSCCondSegFit( pset(), theChamber, csc2DRecHits_input_to_build_segment );
	       segment_fit->setScaleXError( tempcorr ); // reset to previous value (rather than init to 1)
	       segment_fit->fit(condpass1, condpass2);
	       
	     }
	   
	  
	   CSCSegment temp(csc2DRecHits_input_to_build_segment,  segment_fit->intercept(),  segment_fit->localdir(),  segment_fit->covarianceMatrix(),  segment_fit->chi2() );
	   delete segment_fit;
	   segments.push_back(temp);
	   
	  
	 } // end strip hit loop
    } // end wire hit loop







  
  std::cout << "n2DSeg before prune: " << segments.size() <<  std::endl;
  for(auto iseg : segments)std::cout<<"  nRecHits  "<< (iseg.recHits()).size() << std::endl;
  

  
  std::vector<CSCSegment> segments_prune = prune_bad_hits(theChamber, segments); 
  
  
  std::cout << "n2DSeg after prune: " << segments_prune.size() << std::endl;
  for(auto iseg : segments)std::cout<<"  nRecHits  "<< (iseg.recHits()).size() << std::endl;
  return segments_prune;
  
}





void CSCSegAlgoUF::PrintTH2F(TH2F* hist)
{

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
	  cols_v.push_back(wg - 1);

	  rows_v.push_back(wLayer - 1);

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






void
CSCSegAlgoUF::ScanForWireSegment(TH2F* wireHitsInChamber, std::list<CSCWireSegment>& wireSegments, std::vector<TH2F*>& wireSegmentsTH2F, std::vector<int> &wireSegments_rank/* <---- for debugging, to be removed */, int NumberOfLayers)
{


  TH1D* nonEmptyWG = wireHitsInChamber->ProjectionX();
     
  for (int iWireGroup = 0; iWireGroup < nWireGroups; iWireGroup++)
    {
      int thisKeyWG = nWireGroups - iWireGroup;   // scan from the wider chamber side

      
      if (nonEmptyWG->GetBinContent(thisKeyWG)   == 0 && 
	  nonEmptyWG->GetBinContent(thisKeyWG+1) == 0 &&
	  nonEmptyWG->GetBinContent(thisKeyWG+2) == 0 ) continue; // only scan non empty area, what about under/over flow ? need to check ! (Vladimir: what are the under/overflow bins? )

	 

      for (int iPattern = 0; iPattern < nDefinedWirePatterns; iPattern++)  //  loop over 4 predefined wire patterns
	{
	   
	  if(wireHitsInChamber->Integral() == 0) continue;  // Do not proceed if either there are no hits or if no hits left after creating a segment

	  
	  // ========================== fill up the pattern matrix =================================
	  int nWGsInPattern = nWGsInPatterns[iPattern];
	  double w_cols_scan[nWGsInPattern] = {};
	  for (int k = 0; k < nWGsInPattern; k++) // loop over pattern
	    {
	      w_cols_scan[k] = w_cols[iPattern][k] + thisKeyWG - 1;
	    }
	  double w_rows_scan[nWGsInPattern] = {};

	  if (theStation == 3 || theStation == 4) //  invert the order for these chambers whose layers are counted in opposite Z direction
	    for (int k = 0; k < nWGsInPattern; k++){w_rows_scan[k] = 5-w_rows[iPattern][k];}
	  else
	    for (int k = 0; k < nWGsInPattern; k++){ w_rows_scan[k] = w_rows[iPattern][k];}
	     
	  TH2F* wirePattern = new TH2F("wirePattern","", nWireGroups, 0, nWireGroups, 6, 0, 6);
	  wirePattern->FillN(nWGsInPattern, w_cols_scan, w_rows_scan, w_data[iPattern]);
	  //=========================================================================================



	   


	  wirePattern->Multiply(wireHitsInChamber);           // Multiply pattern and actual hits in the chamber
	  TH2F* actualSegment  = wirePattern;                 // intersection between pattern and hits;  after multiplied wirePattern becomes an actual segment

	  
	  TH1D* hitsLayer = wirePattern->ProjectionY();       //<- this histo is Y projection, so it's layers with hits
	  hitsLayer->Divide(hitsLayer);                       //  normalised to itself
	  int nLayersOfSegment = hitsLayer->Integral();       //  how many layers with hits 


	   
	  if (nLayersOfSegment !=  NumberOfLayers)  continue; // break the pattern loop if N layers with hits != to requested: 6,5,4,3

	   

	  if (  wireSegments.size() == 0  )                   // if the first segment; (that's a bit weird logic)
	    {
	      wireSegments.push_back( CSCWireSegment(thisKeyWG, nLayersOfSegment, wirePattern) );  // create a segment
	      wireSegmentsTH2F.push_back(actualSegment);                                           // keep TH2F for control and debug
	      wireSegments_rank.push_back(iPattern);
	      wireHitsInChamber->Add( wirePattern, -1 );                                           // remove hits belonging to the segment from chamber hits
	      
	      continue;                                                                            // continue loop over patterns
	    }

	     
	     CSCWireSegment lastSegment                         =  wireSegments.back();            // Last segment in vector
	     int            lastKeyWG                           =  lastSegment.keyWG();            // its keyWG; To be checked, in text it says that it's WG in the 3rd layer, but i dont see how.
	     
	     CSCWireSegment potentialAnotherSegment             =  CSCWireSegment(thisKeyWG, nLayersOfSegment, wirePattern);
	     

	     
	     
	     if (abs(lastKeyWG - thisKeyWG) >  1)                 // if two segments are far away by more than 1 WG -> create a new segment;
	       {
		 
		 wireSegments.push_back(potentialAnotherSegment); //  if the WG difference between segment > 1 create  new segment
		 wireSegmentsTH2F.push_back(actualSegment);
		 wireSegments_rank.push_back(iPattern);
		 
	       }


	   
	     /*	     /// remove for now; Here he tried to merge two segments into a wider segment if by keyWG they dont differ by 1; See comment above how the keyWG is defined;
	     if (abs(lastKeyWG -  thisKeyWG) == 1)
	       {
		 
		 wireSegments.back().updateWHits( potentialAnotherSegment.wireHitsPosition(), potentialAnotherSegment.nLayerHits());
	       
	       }
	     */

	     
	     
	     wireHitsInChamber->Add(wirePattern, -1); // delete wire group being used
	     
	   } // pattern loop
       } // WG loop
     // std::cout<<" ===>  End  ScanForWireSegment    "<< std::endl;
}






void CSCSegAlgoUF::ScanForStripSegment(TH2F* stripHitsInChamber, std::list<CSCStripSegment>& stripSegments, std::vector<TH2F*>& stripSegmentsTH2F, std::vector<int>& stripSegments_rank/* <---- for debugging, to be removed */, int nLayer)
{

     TH1D* nonEmptyStrip = stripHitsInChamber->ProjectionX();
     for (int iStrip = 0; iStrip < nStrips*2+1; iStrip++)
       {

       
	 int thisKeyHalfStrip = iStrip + 1; // why +1 ? 
       
	 if (nonEmptyStrip->GetBinContent(iStrip+1) == 0 && 
	     nonEmptyStrip->GetBinContent(iStrip)   == 0 &&
	     nonEmptyStrip->GetBinContent(iStrip+2) == 0 ) continue; // only scan non empty area, what about under/over flow ? need to check !


	 for (unsigned int iStripPattern = 0; iStripPattern < nDefinedStripPatterns;   iStripPattern++ )
	   {
	     // scan patterns sequentially, according to ranks
	     // scan from high rank (high pt) patterns to low pt patterns //
	   
	     if( stripHitsInChamber->Integral()==0)       continue;  // dont continue if there are no hits or no hits left when strip segment is formed 



	     // ========================== fill up the pattern matrix =================================
	     int nhalfStrips     = nHalfStrips[iStripPattern]; 
	     int thisRank        = patternRanks[iStripPattern];
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
	     if (theStation == 3 || theStation == 4)         // invert layers numbering for these stations
	       for (int k = 0; k < nhalfStrips; k++) s_rows_scan[k] = 5-s_rows[iStripPattern][k];
	     else
	       for (int k = 0; k < nhalfStrips; k++) s_rows_scan[k] = s_rows[iStripPattern][k];
	     
	     TH2F* stripPattern = new TH2F("stripPattern","", nStrips*2+1, 0, nStrips*2+1, 6, 0, 6);
	     stripPattern->FillN(nhalfStrips, s_cols_scan, s_rows_scan, s_data[iStripPattern]);
	     // =========================================================================================
	     



	     
	     stripPattern->Multiply(stripHitsInChamber);    // Multiply pattern and actual hits in the chamber
	     TH2F* actualSegment  = stripPattern;


	     
	     TH1D* hitsLayer = stripPattern->ProjectionY(); //<- this histo is Y projection, so it's layers with hits
	     hitsLayer->Divide(hitsLayer);                  //  normalised to itself
	     int thisNLayer = hitsLayer->Integral();        //  how many layers with hits 


	     

	     if (thisNLayer != nLayer) continue;            // continue the pattern loop if N layers with hits != to requested: 6,5,4,3

	     

	     
	     if ( stripSegments.size() == 0 )
	       {
		 
		 stripSegments.push_back( CSCStripSegment(thisKeyHalfStrip, thisNLayer, thisRank, stripPattern) ); 
		 stripHitsInChamber->Add(stripPattern,-1);
		 stripSegmentsTH2F.push_back(actualSegment);
		 stripSegments_rank.push_back(iStripPattern);
	       
		 continue;
	       }

	     
	     CSCStripSegment lastSegment = stripSegments.back();
	     int lastKeyHalfStrip        = lastSegment.keyHalfStrip();


	     CSCStripSegment potentialAnotherSegment = CSCStripSegment(thisKeyHalfStrip,thisNLayer,thisRank,stripPattern);


	     
	     if (abs(thisKeyHalfStrip - lastKeyHalfStrip) > 1)
	       {
		 stripSegments.push_back(potentialAnotherSegment);
		 stripSegmentsTH2F.push_back(actualSegment);
		 stripSegments_rank.push_back(iStripPattern);
		 
	       }


	     /* ///  remove for now 
		if (abs(thisKeyHalfStrip - lastKeyHalfStrip) == 1)
		  {
		
		  stripSegments.back().updateSHits(tmpStripSeg.stripHits(), tmpStripSeg.nLayerHits());
		
		  }
	     */

	     
	     
             stripHitsInChamber->Add(stripPattern,-1); 
	     

	   } // pattern loop
       } // WG loop
     //std::cout<<" ===>  End  ScanForStripSegment    "<< std::endl;
}










void CSCSegAlgoUF::GetWireHitFromWireSegment(CSCWireSegment wireSegment, ChamberWireHitContainer WireHitsInChamber, int* wireHitIndex) {



  
  
  double lowerWireGroup  = wireSegment.LowestHitInLayer();
  double higherWireGroup = wireSegment.HighestHitInLayer();
  double segment_keyWG   = wireSegment.keyWG();     // to check what it really returns

  
  //std::cout << "keyWH: " << keyWH << std::endl;
  //std::cout << "low: " << wgLow << ", high: " << wgHigh << std::endl;



  std::vector<int> addBackCounter;               // it counts layers w/o hits

  for (int iLayer = 0; iLayer < 6; iLayer++)     // loop over 6 layers
    { 
      double wireHitPosition_from_segment = (wireSegment.wireHitsPosition())[iLayer];
      int wireHitIndex_from_hits_collection = -1; // to be found by matching, but why??
      
      double wireHitPositionDelta = MaxWGNumber;         //112 is max one can get; Max WG number  = 112 in ME21

      bool hitMissed_in_segment = false; // check if the segment does not have hit in this layer

      

      
      for (unsigned int hit_index = 0; hit_index < WireHitsInChamber.size(); hit_index++) // loop over all hits 
	{
	  
          const CSCWireHit* tmpHit = WireHitsInChamber[hit_index];        //  this is the Chamber Wire Hit
          int  wireHitLayer = tmpHit->cscDetId().layer();                 // 
          float wireHitPosition_from_hits_collection = tmpHit->wHitPos(); // why is it float and may take non-integer values ?; That's should be a wire group; 


          if (wireHitLayer != (iLayer+1) ) continue;  // skip if not in this layer

	  
	  // check if there is a missing layer in the segment but there are actual hits  within low-high segment boundaries	  
          if (wireHitPosition_from_segment == 0 && (wireHitPosition_from_hits_collection >= lowerWireGroup  - 1) &&
	                                           (wireHitPosition_from_hits_collection <= higherWireGroup + 1) )

	    {
	      addBackCounter.push_back(iLayer);   // store layer number where there no hit from Wire segment  and no wire hit
	      wireHitIndex_from_hits_collection = hit_index;
	      hitMissed_in_segment = true;
	      
	      break;
	     
	    }


	  

	  // look for the hit from Wire Hit Collection closest to the hit in segment
          if (abs( wireHitPosition_from_segment  -  wireHitPosition_from_hits_collection ) < wireHitPositionDelta)
	    {

	      
	      wireHitPositionDelta = abs( wireHitPosition_from_segment  -  wireHitPosition_from_hits_collection );
	      wireHitIndex_from_hits_collection = hit_index; // find the closest hit

	      
	    }
	} // loop over hits collection



      if( wireHitPositionDelta < MaxWGNumber && wireHitPosition_from_segment > 0 ) // if there is a hit in a segment in the given layer and any (!) hit in chamber return hit from hit collection
	{

	  wireHitIndex[iLayer] = wireHitIndex_from_hits_collection;
	  
	}
      else if(hitMissed_in_segment && wireHitPosition_from_segment == 0)
	{
	  
	  wireHitIndex[iLayer] = wireHitIndex_from_hits_collection;
	  
	}
      else
	{
	  
	  wireHitIndex[iLayer] = -1;
	  
	}


      /*
      if (  (  wireHitPositionDelta < 113 && wireHitPosition_from_segment > 0  )    ||    ( wireHitPosition_from_segment == 0 && hitMissed_in_segment )   ) // this requires more debugging, delta is just < 113 ??? 
	// if pattern doesn't cover some hit, and is close, to keyWH, include it // <-- It's not true!
	{
	  
	  wireHitIndex[iLayer] = wireHitIndex_from_hits_collection;
	  //	  std::cout<<"    Layer         "<<iLayer<<"     " <<wireHitIndex[iLayer] << "  position " << WireHitsInChamber[wireHitIndex_from_hits_collection]->wHitPos()<< "  and delta   " << wireHitPositionDelta  <<std::endl;
	  //	  if(hitMissed_in_segment)
	    //	    std::cout<<"  Missing hits =================    Layer         "<< iLayer <<"     " << "  position " << WireHitsInChamber[wireHitIndex_from_hits_collection]->wHitPos()<< "  and delta   " << wireHitPositionDelta  <<std::endl;
	  
	}
      
      else
	
	{
	  
	  wireHitIndex[iLayer] = -1;
	  
	}
      */

      
      
    }// loop over layers


  
  if (addBackCounter.size() > 0 ) std::cout << "Wire add back " << addBackCounter.size() << " times" << std::endl; 

  if (wireSegment.nLayersWithHits() == 3 && addBackCounter.size() > 0 ) std::cout << "add back recover ALCT:   " << addBackCounter.size() << " times" << std::endl;

}










void CSCSegAlgoUF::GetStripHitFromStripSegment(CSCStripSegment stripSegment, ChamberStripHitContainer stripHitsInChamber, int* stripHitIndex)
{


  double LowestHitInLayer  = stripSegment.LowestHitInLayer(isME11);
  double HighestHitInLayer = stripSegment.HighestHitInLayer(isME11);
  
  //  double keySH = ceil(stripSegment.keyHalfStrip()/2.0);
  //  std::cout << "keySH: " << keySH << std::endl;

  std::vector<int> addBackCounter;

  for (int iLayer = 0; iLayer < 6; iLayer++) {
    double stripHitPosition_from_segment_ = (stripSegment.stripHits())[iLayer]; // comparator position in unit of half strip
    double stripHitPosition_from_segment = -1;                             // strip hit position in unit of strip
 


    if (!isME11                   && (iLayer==0 || iLayer==2 || iLayer==4) )         stripHitPosition_from_segment = ceil((stripHitPosition_from_segment_ - 1)/2.0); // convert comparator number to strip number
    if (isME11        || (!isME11 && (iLayer==1 || iLayer==3 || iLayer==5) ) )       stripHitPosition_from_segment = ceil(stripHitPosition_from_segment_/2.0);
    if (stripHitPosition_from_segment_ == 1 ||  stripHitPosition_from_segment_==2)   stripHitPosition_from_segment = 1;

      

    int    stripHitIndex_from_hits_collection = -1;
    double stripHitPositionDelta = MaxStripNumber;     //MaxStripNumber = 81; max strip number is 80
    bool   hitMissed_in_segment = false;
      

    for (unsigned int hit_index = 0; hit_index < stripHitsInChamber.size(); hit_index++)
	{
	  
          const CSCStripHit* tmpHit = stripHitsInChamber[hit_index];
          int stripHitLayer         = tmpHit->cscDetId().layer();
          double stripHitPosition_from_hits_collection = tmpHit->sHitPos();
	  

          if ( stripHitLayer != iLayer + 1) continue;

	  
          if (stripHitPosition_from_segment == 0 /*&& abs(keySH - stripHitPosition_from_hits_collection)<2*/ && ((2*stripHitPosition_from_hits_collection) >= LowestHitInLayer-2) && (2*stripHitPosition_from_hits_collection-1 <= HighestHitInLayer+2) )
	    {
	      addBackCounter.push_back(iLayer);
	      stripHitIndex_from_hits_collection = hit_index;
	      
	      hitMissed_in_segment = true;
	      
	      break;
	    }

	  
	  //  find the closest hit in hits collection
          if (abs(  stripHitPosition_from_segment -  stripHitPosition_from_hits_collection) < stripHitPositionDelta)
	    {
	      
	      stripHitPositionDelta = abs(stripHitPosition_from_segment-stripHitPosition_from_hits_collection);
	      stripHitIndex_from_hits_collection = hit_index;
	    
	    }
	}



    if( stripHitPositionDelta < MaxStripNumber   &&   stripHitPosition_from_segment >= 1)
      {
	stripHitIndex[iLayer] = stripHitIndex_from_hits_collection;
      }
    
    else if(stripHitPosition_from_segment == 0 && hitMissed_in_segment)
      {
	stripHitIndex[iLayer] = stripHitIndex_from_hits_collection;
      }
    
    else
      {
	stripHitIndex[iLayer] = -1;
      }
      
  }

  if (int(addBackCounter.size())>0) std::cout << "Strip add back " << addBackCounter.size() << " times" << std::endl;
  if (int(addBackCounter.size())>0 && stripSegment.nLayersWithHits() == 5 && stripSegment.patternRank() == 1) std::cout << "5 hits rank1 strip add back " << addBackCounter.size() << " times" << std::endl;
  if (int(addBackCounter.size())>0 && stripSegment.nLayersWithHits() == 5 ) std::cout << "5 hits all ranks strip add back " << addBackCounter.size() << " times" << std::endl;
  if (stripSegment.nLayersWithHits()==3 && int(addBackCounter.size()) >0  ) std::cout << "add back recover CLCT" << addBackCounter.size() << " times" << std::endl;
  
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
      CSCCondSegFit* segment_fit = new CSCCondSegFit( pset(), theChamber/*chamber() hm*/, protoSegment );
      condpass1 = false;
      condpass2 = false;
      segment_fit->setScaleXError( 1.0 );
      segment_fit->fit(condpass1, condpass2);

      // Attempt to handle numerical instability of the fit;
      // The same as in the build method;
      // Preprune is not applied;

      
      if( true/*adjustCovariance() hm*/){ // if true? :) Seriously ? 
	if(segment_fit->chi2()/segment_fit->ndof()>10/*chi2Norm_3D_*/){
	  condpass1 = true;
	  segment_fit->fit(condpass1, condpass2);
	}
	if( (segment_fit->scaleXError() < 1.00005)&&(segment_fit->chi2()/segment_fit->ndof()>10/*chi2Norm_3D_*/) ){
	  condpass2 = true;
	  segment_fit->fit(condpass1, condpass2);
	}
      }
    
      // calculate error matrix 
      //      AlgebraicSymMatrix temp2 = segfit->covarianceMatrix();

      // build an actual segment

      
      CSCSegment temp(protoSegment, segment_fit->intercept(), segment_fit->localdir(), 
                         segment_fit->covarianceMatrix(), segment_fit->chi2() );

      // and finished with this fit
      delete segment_fit;

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






