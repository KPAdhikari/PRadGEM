#include "GEMZeroHitDecoder.h"

// split large cluster or not
#define SPLITCLUSTER
#undef SPLITCLUSTER

#define CheckChargeRatio
#undef CheckChargeRatio

GEMZeroHitDecoder::GEMZeroHitDecoder(uint32_t *rbuf, int Size)
{
  // if no more sigma cut, use this constructor
  // no sigma cut -> so no need to load pedestal data
  fIsHitMaxOrTotalADCs = "signalPeak"; 

  buf = rbuf;
  fSize = Size;
  //cout<<"GEMZeroHitDecoder:: event size:  "<<fSize<<endl;
  cut_buf = NULL;
  cut_fSize = 0;

  fMapping = PRDMapping::GetInstance();
  ped = NULL; // pass pedestal file

  fFECID = 0;
  fADCChannel = 0;
  fAPVID = 0;
  fAPVKey = 0;
  fZeroSupCut = 5;
  nTimeBin = 3;

  fListOfHitsZero.clear();
  fListOfHitsZeroFromPlane.clear();
  
  // cluster information 
  fMinClusterSize = 1;
  fMaxClusterSize = 20;
  fIsClusterMaxOrTotalADCs = "totalADCs";
  //fIsClusterMaxOrTotalADCs = "maximumADCs";
  fIsGoodClusterEvent = kFALSE;
  fListOfClustersZeroFromPlane.clear();

  ProcessEvent();
  //ped->Delete();
  //
}

GEMZeroHitDecoder::GEMZeroHitDecoder(uint32_t *rbuf, int Size, GEMPedestal *pPed)
{
   // If need to do one more sigma cut, use this constructor.
   // need pass pedestal to do sigma cut
  fIsHitMaxOrTotalADCs = "signalPeak"; 

  cut_buf = rbuf;
  cut_fSize = Size;
  buf = NULL;
  fSize = 0;
  //cout<<"GEMZeroHitDecoder:: event size:  "<<fSize<<endl;

  fMapping = PRDMapping::GetInstance();
  ped = pPed; // pass pedestal file

  fFECID = 0;
  fADCChannel = 0;
  fAPVID = 0;
  fAPVKey = 0;
  fZeroSupCut = 5;
  nTimeBin = 3;

  fListOfHitsZero.clear();
  fListOfHitsZeroFromPlane.clear();
  
  // cluster information   
  fMinClusterSize = 1;
  fMaxClusterSize = 20;
  fIsClusterMaxOrTotalADCs = "totalADCs";
  //fIsClusterMaxOrTotalADCs = "maximumADCs";
  fIsGoodClusterEvent = kFALSE;
  fListOfClustersZeroFromPlane.clear();

  ProcessEvent();
  //ped->Delete();
}

GEMZeroHitDecoder::~GEMZeroHitDecoder()
{
  //clear hits
  map<int, GEMHit *>::iterator it = fListOfHitsZero.begin();
  for(;it!=fListOfHitsZero.end();++it)
  {
    delete (it->second);
  }

  //clear clusters
  if( fListOfClustersZeroFromPlane.size() > 0)
  {
    map<TString, list<GEMCluster*> >::iterator itt = fListOfClustersZeroFromPlane.begin();
    for(;itt!=fListOfClustersZeroFromPlane.end();++itt)
    {
      list<GEMCluster*>::iterator itc = (itt->second).begin();
      for(;itc!=(itt->second).end();++itc)
        delete *itc;
      (itt->second).clear();
    }
    
  }
  fListOfClustersZeroFromPlane.clear();
  fListOfHitsZero.clear();
  fListOfHitsZeroFromPlane.clear();

  DeleteClustersInPlaneMap();
  //clear buffer
  if(cut_buf != NULL) delete[] buf;
}

void GEMZeroHitDecoder::ProcessEvent( )
{
  Cut();
  EventHandler( );
}

void GEMZeroHitDecoder::Cut()
{
  if(cut_buf == NULL) return;

  vector<uint32_t> strip;
  strip.clear();

  int pol = 0;
  int adc = 0;
  int TS = 0;
  int chNo = 0;

  int totalEntry = 0;

  for(int i=0;i<cut_fSize;i++)
  {
    //GEM data identifier by Sergey, no use for now.
    if(cut_buf[i] == 0xfecfec99) 
    {
        totalEntry++; 
	continue;
    } 
    // check overflow adc value
    if( (cut_buf[i] & 0x7ffff000) == 0x7ffff000 ) 
    {
        totalEntry++; 
	cout<<"7ffff: "<<cut_buf[i]
	    <<endl;
	continue;
    } 

    fFECID = ( (cut_buf[i]>>26) & 0xf);
    fADCChannel = ( (cut_buf[i]>>22) & 0xf);
    chNo = ( (cut_buf[i]>>15) & 0x7f );
    TS = ( (cut_buf[i]>>12) & 0x7 );
    pol = ( (cut_buf[i]>>11) & 0x1);
    adc = ( (cut_buf[i]) & 0x7ff );
    if(pol == 1) 
        adc = -adc;
 
    int cnt = 1;
    int avg_adc = adc; 

     // Note: If Time Sample 0 missing, temporarily throw away the currrent event
    if(TS == 0)
    {
      totalEntry ++;
      vector<uint32_t> s;
      s.clear();
      s.push_back(cut_buf[i]);
      for(int j=i+1;j<cut_fSize;j++)
      {
        if( ( ( (cut_buf[j]>>26) & 0xf) == fFECID ) && 
	    (fADCChannel == ( (cut_buf[j]>>22) & 0xf) ) && 
	    (chNo == ( (cut_buf[j]>>15) & 0x7f )) &&  
	    ( (( (cut_buf[j]>>12) & 0x7 ) == 1 ) || 
	    (( (cut_buf[j]>>12) & 0x7 ) == 2 ) ))
	{
	    s.push_back(cut_buf[j]);
	    int p = ( (cut_buf[j]>>11) & 0x1);
	    int val = ( (cut_buf[j]) & 0x7ff );
            if(p == 1) 
	        val = -val;
	    avg_adc += val;
	    cnt ++;
	    totalEntry++;
	    if(cnt == 4) 
	        break;
	}
      }
      avg_adc/=cnt;
      fAPVID = (fFECID<<4)|fADCChannel;
      fPedestalNoises = ped ->GetAPVNoises(fAPVID); 
      if(avg_adc > fZeroSupCut * fPedestalNoises[chNo] ) 
      {
        //strip.insert(strip.end(), s.begin(), s.end());
	//vector size < 12, push_back faster than insert; vector size > 12, insert faster than push_back
	int N = s.size();
	for(int i=0;i<N;i++)
	{
	  strip.push_back(s[i]);
	}
      }
    }

  }
  fSize = strip.size();
  buf = new uint32_t[fSize];
  for(int i=0;i<fSize;i++)
  {
    buf[i] = strip[i];
  }
  /* 
  //debug
  cout<<"after 5 sigma cut: "<<strip.size()<<"  before: "<<fSize<<endl;
  cout<<"totalEntry:  "<<totalEntry<<endl;
  //assert(totalEntry == fSize);
   
  //debug
  for(int i=0;i<fSize;i++)
  {
    if(cut_buf[i] == 0xfecfec99) {continue;} //GEM data identifier by Sergey, no use for now.
    if( (cut_buf[i] & 0x7ffff000) == 0x7ffff000 ) { cout<<"7ffff: "<<cut_buf[i]<<endl;continue;} // check overflow adc value

    fFECID = ( (cut_buf[i]>>26) & 0xf);
    fADCChannel = ( (cut_buf[i]>>22) & 0xf);
    chNo = ( (cut_buf[i]>>15) & 0x7f );
    TS = ( (cut_buf[i]>>12) & 0x7 );
    pol = ( (cut_buf[i]>>11) & 0x1);
    adc = ( (cut_buf[i]) & 0x7ff );
    if(pol == 1) adc = -adc;
    cout<<"fec: "<<fFECID<<" adc: "<<fADCChannel<<" strip: "<<chNo<<" timebin: "<<TS<<" adc: "<<adc<<endl;
  }
  cout<<"============"<<endl;
  for(int i=0;i<strip.size();i++)
  {
    fFECID = ( (strip[i]>>26) & 0xf);
    fADCChannel = ( (strip[i]>>22) & 0xf);
    chNo = ( (strip[i]>>15) & 0x7f );
    TS = ( (strip[i]>>12) & 0x7 );
    pol = ( (strip[i]>>11) & 0x1);
    adc = ( (strip[i]) & 0x7ff );
    if(pol == 1) adc = -adc;
    cout<<"fec: "<<fFECID<<" adc: "<<fADCChannel<<" strip: "<<chNo<<" timebin: "<<TS<<" adc: "<<adc<<endl;
    
  }
  */
}

void GEMZeroHitDecoder::EventHandler( )
{
  //cout<<"GEMHitDecoder::EventHandler begin..."<<endl;
  //Cut();

  /* hit structure
   * det: 1 bit
   * plane: 1 bit
   * fec: 4 bit
   * adcchannel: 4 bit
   * strip: 7 bit
   * time sample: 3 bit
   * polarity: 1 bit
   * val: 11 bit
   */

  //Fill Hits
  int pol = 0;
  int adc = 0;
  int TS = 0;
  int chNo = 0;

  for(int i=0;i<fSize;i++)
  {
    //GEM data identifier by Sergey, no use for now.
    if(buf[i] == 0xfecfec99) 
        continue; 

    // Temporary solution
    if( (buf[i] & 0x7ffff000) == 0x7ffff000 ) 
    {
        cout<<"7ffff: "<<buf[i]
	    <<endl;
	continue;
    } 
    fFECID = ( (buf[i]>>26) & 0xf);
    fADCChannel = ( (buf[i]>>22) & 0xf);
    chNo = ( (buf[i]>>15) & 0x7f );
    TS = ( (buf[i]>>12) & 0x7 );
    pol = ( (buf[i]>>11) & 0x1);
    adc = ( (buf[i]) & 0x7ff );

    if(pol == 1) 
        adc = -adc;

    fAPVID = (fFECID<<4)|fADCChannel;
    fAPVKey = fMapping->GetAPVNoFromID(fAPVID);
    int hitID = (fAPVKey << 8) | chNo ;

    //debug
    //cout<<"fec: "<<fFECID<<" adc: "<<fADCChannel<<" strip: "<<chNo<<" timebin: "<<TS<<" adc: "<<adc<<endl;

    if( !fListOfHitsZero[hitID])
    {
      //cout<<hitID<<"  "<<fAPVID<<"  "<<chNo<<"  "<<fZeroSupCut<<"  "<<fIsHitMaxOrTotalADCs<<endl;
      GEMHit * hit = new GEMHit(hitID, fAPVID, chNo, fZeroSupCut, fIsHitMaxOrTotalADCs);
      fListOfHitsZero[hitID] = hit;
      //cout<<"GEMZeroHitDecoder:: EventHandler:  "<<endl;
    }
    fListOfHitsZero[hitID] -> AddTimeBinADCs(TS, adc);
  }


  GetListOfHitsZeroFromPlanes();

  // cluster computing       
#ifndef SPLITCLUSTER
  ComputeClusters();
#endif
#ifdef SPLITCLUSTER
  ComputeClustersFineTune();
#endif
}

map < TString, list <GEMHit * > > GEMZeroHitDecoder::GetListOfHitsZeroFromPlanes() {
  map < Int_t, GEMHit * >::const_iterator hit_itr ;
  for(hit_itr = fListOfHitsZero.begin(); hit_itr != fListOfHitsZero.end(); ++hit_itr) { 
    GEMHit * hit = (* hit_itr).second ;
    TString planename = hit->GetPlane() ;
    fListOfHitsZeroFromPlane[planename].push_back(hit) ;
    int adc = hit->GetHitADCs() ;
  }
  return fListOfHitsZeroFromPlane ;
}

TH1F* GEMZeroHitDecoder::GetZeroHit(TString str)
{
  TH1F * h1;

  int NN = fMapping->GetNbOfPlane();

  int nbDetector = fMapping->GetNbOfDetectors();

  for(int i=0;i<nbDetector;i++)
  {
    TString detectorName = fMapping->GetDetectorFromID(i);
    list<TString> planeList = fMapping->GetPlaneListFromDetector(detectorName);
    int nbPlane = planeList.size();

    list<TString>::iterator it;
    for(it=planeList.begin();it!=planeList.end();++it)
    {
      if(*it == str)
      {
        TString hh = detectorName+"_"+(*it)+"_hit_distribution_zero_suppression";
	h1 = new TH1F(hh, hh, 2000, -fMapping->GetPlaneSize(*it)/2-100, fMapping->GetPlaneSize(*it)/2+100 );
        list< GEMHit* > hitList = fListOfHitsZeroFromPlane[ *it  ];
        list< GEMHit* >::iterator hit_it;
        for(hit_it=hitList.begin(); hit_it!=hitList.end();++hit_it)
        {
          Float_t pos = (*hit_it) -> GetStripPosition();
	  Float_t adc = (*hit_it) -> GetHitADCs(); 
	  h1 -> Fill(pos, adc);
        }
      }
    }
  }
  return h1;
}

/*****************************************************************************
 * **                    Compute Clusters Information                     ** *
 *****************************************************************************/
static Bool_t CompareStripNo( TObject *obj1, TObject *obj2) {
  Bool_t compare ;
  if ( ( (GEMHit*) obj1 )->GetStripNo() < ( ( GEMHit*) obj2 )->GetStripNo() ) compare = kTRUE ;
  else compare = kFALSE ;
  return compare ;
}

static Bool_t CompareHitADCs( TObject *obj1, TObject *obj2) {
  Bool_t compare ;
  if ( ( (GEMHit*) obj1 )->GetHitADCs() > ( ( GEMHit*) obj2 )->GetHitADCs()) compare = kTRUE ;
  else compare = kFALSE ;
  return compare ;
}

static Bool_t CompareClusterADCs( TObject *obj1, TObject *obj2) {
  Bool_t compare ;
  if ( ( (GEMCluster*) obj1 )->GetClusterADCs() > ( ( GEMCluster*) obj2 )->GetClusterADCs()) compare = kTRUE ;
  else compare = kFALSE ;
  return compare ;
}

// rearrange cluster by its size
static Bool_t CompareClusterSize( TObject *obj1, TObject *obj2) {
  Bool_t compare ;
  if ( ( (GEMCluster*) obj1 )->GetNbOfHits() > ( ( GEMCluster*) obj2 )->GetNbOfHits()) compare = kTRUE ;
  else compare = kFALSE ;
  return compare ;
}

void GEMZeroHitDecoder::ComputeClusters()
{
  map < TString, list <GEMHit*> >::const_iterator  hitsFromPlane_itr ;

  for (hitsFromPlane_itr = fListOfHitsZeroFromPlane.begin(); hitsFromPlane_itr != fListOfHitsZeroFromPlane.end(); ++hitsFromPlane_itr) 
    {
      TString plane =  (*hitsFromPlane_itr).first ;
      list <GEMHit*> hitsFromPlane = (*hitsFromPlane_itr).second ; 
      hitsFromPlane.sort(CompareStripNo) ;
      Int_t listSize = hitsFromPlane.size() ;

      if (listSize < fMinClusterSize) 
	{
	  fIsGoodClusterEvent = kFALSE ;
	  continue ;
	}

      Int_t previousStrip = -2 ;
      Int_t clusterNo = -1 ;
      map<Int_t, GEMCluster *> clustersMap ;
      list <GEMHit *>::const_iterator hit_itr ;

      for (hit_itr = hitsFromPlane.begin(); hit_itr != hitsFromPlane.end(); hit_itr++) 
	{
	  GEMHit * hit =  * hit_itr ; 
	  
	  Int_t currentStrip = hit->GetStripNo() ;
	  //cout<<"ComputeClusters:  "<<currentStrip<<endl;
	  
	  // remove first 16 strips (apv index 0 on X side) and last 16 strips (apv index 10 on X side)
	  if( plane.Contains("X") && ( (currentStrip<16) || (currentStrip > 1391) )  ) continue;

	  if(currentStrip != (previousStrip + 1)) 
	    {
	      clusterNo++ ;
	    }
	  if(!clustersMap[clusterNo]) 
	    {
	      clustersMap[clusterNo] = new GEMCluster(fMinClusterSize, fMaxClusterSize, fIsClusterMaxOrTotalADCs) ;
	      clustersMap[clusterNo]->SetNbAPVsFromPlane(hit->GetNbAPVsFromPlane());
	      clustersMap[clusterNo]->SetAPVIndexOnPlane(hit->GetAPVIndexOnPlane());
	      clustersMap[clusterNo]->SetPlaneSize(hit->GetPlaneSize());
	      clustersMap[clusterNo]->SetPlane(hit->GetPlane());
	    }
	  clustersMap[clusterNo]->AddHit(hit) ;
	  previousStrip = currentStrip;
	}

      map<Int_t, GEMCluster *>::const_iterator  cluster_itr ;
      for (cluster_itr = clustersMap.begin(); cluster_itr != clustersMap.end(); cluster_itr++) 
	{
	  GEMCluster * cluster = ( * cluster_itr ).second ;
	  if (!cluster->IsGoodCluster()) 
	    {
	      delete cluster ;
	      continue ;
	    }
	  cluster->ComputeClusterPosition() ;
	  fListOfClustersZeroFromPlane[plane].push_back(cluster) ;
	}
       
      fListOfClustersZeroFromPlane[plane].sort(CompareClusterADCs) ;
      //xb : sort cluster by their size
      //fListOfClustersZeroFromPlane[plane].sort(CompareClusterSize) ;
      hitsFromPlane.clear() ;
      clustersMap.clear() ;
    }

}


void GEMZeroHitDecoder::ComputeClustersFineTune()
{
  map < TString, list <GEMHit*> >::const_iterator  hitsFromPlane_itr ;

  for (hitsFromPlane_itr = fListOfHitsZeroFromPlane.begin(); hitsFromPlane_itr != fListOfHitsZeroFromPlane.end(); ++hitsFromPlane_itr) 
    {
      TString plane =  (*hitsFromPlane_itr).first ;
      list <GEMHit*> hitsFromPlane = (*hitsFromPlane_itr).second ; 
      hitsFromPlane.sort(CompareStripNo) ;
      Int_t listSize = hitsFromPlane.size() ;

      if (listSize < fMinClusterSize) 
	{
	  fIsGoodClusterEvent = kFALSE ;
	  continue ;
	}

      Int_t previousStrip = -2 ;
      Int_t clusterNo = -1 ;
      map<Int_t, GEMCluster *> clustersMap ;
      list <GEMHit *>::const_iterator hit_itr ;

      //split big clusters 
      int leading=0;
      int trailing = 0;
      int newHitFlag = 0;
      GEMHit * previousHit = NULL;
      GEMHit * currentHit = NULL;

      for (hit_itr = hitsFromPlane.begin(); hit_itr != hitsFromPlane.end(); hit_itr++) 
	{
	  GEMHit * hit =  * hit_itr ; 

	  //split big clusters
	  if(currentHit == NULL) previousHit = *hit_itr;
	  else previousHit = currentHit;
	  currentHit = *hit_itr;
	  
	  Int_t currentStrip = hit->GetStripNo() ;

          //split big clusters
	  Float_t previousADC = previousHit->GetHitADCs();
	  Float_t currentADC = currentHit->GetHitADCs();
	  if(currentADC >= previousADC) 
	    { 
	      if(trailing == 1) newHitFlag = 1; 
	      else newHitFlag = 0;
	      leading = 1; trailing=0;
	    }
	  else if(currentADC < previousADC) 
	    { 
	      leading=0; 
	      trailing =1;
	    }

	  if(currentStrip != (previousStrip + 1)) 
	    {
	      clusterNo++ ;
	    }
	  //split big clusters
	  else
	    {
	      if(newHitFlag==1) clusterNo++;
	    }

	  if(!clustersMap[clusterNo]) 
	    {
	      clustersMap[clusterNo] = new GEMCluster(fMinClusterSize, fMaxClusterSize, fIsClusterMaxOrTotalADCs) ;
	      clustersMap[clusterNo]->SetNbAPVsFromPlane(hit->GetNbAPVsFromPlane());
	      clustersMap[clusterNo]->SetAPVIndexOnPlane(hit->GetAPVIndexOnPlane());
	      clustersMap[clusterNo]->SetPlaneSize(hit->GetPlaneSize());
	      clustersMap[clusterNo]->SetPlane(hit->GetPlane());
	    }
	  clustersMap[clusterNo]->AddHit(hit) ;
	  previousStrip = currentStrip;
	}

      map<Int_t, GEMCluster *>::const_iterator  cluster_itr ;
      for (cluster_itr = clustersMap.begin(); cluster_itr != clustersMap.end(); cluster_itr++) 
	{
	  GEMCluster * cluster = ( * cluster_itr ).second ;
	  if (!cluster->IsGoodCluster()) 
	    {
	      delete cluster ;
	      continue ;
	    }
	  cluster->ComputeClusterPosition() ;
	  fListOfClustersZeroFromPlane[plane].push_back(cluster) ;
	}

      fListOfClustersZeroFromPlane[plane].sort(CompareClusterADCs) ;
      // xb: sort cluster by their size
      //fListOfClustersZeroFromPlane[plane].sort(CompareClusterSize) ;
      hitsFromPlane.clear() ;
      clustersMap.clear() ;
    }

}


TH1F* GEMZeroHitDecoder::GetCluster(TString str)
{
  TH1F * hc1;

  int NN = fMapping->GetNbOfPlane();

  int nbDetector = fMapping->GetNbOfDetectors();

  for(int i=0;i<nbDetector;i++)
  {
    TString detectorName = fMapping->GetDetectorFromID(i);
    list<TString> planeList = fMapping->GetPlaneListFromDetector(detectorName);
    int nbPlane = planeList.size();

    list<TString>::iterator it;
    for(it=planeList.begin();it!=planeList.end();++it)
    {
      if(*it == str)
      {
        TString hh = detectorName+"_"+(*it)+"_Cluster_Distribution_zero_suppression";
	hc1 = new TH1F(hh, hh, 2000, -fMapping->GetPlaneSize(*it)/2-100, fMapping->GetPlaneSize(*it)/2+100 );
        list< GEMCluster* > clusterList = fListOfClustersZeroFromPlane[ *it  ];
        list< GEMCluster* >::iterator cluster_it;
        for(cluster_it=clusterList.begin(); cluster_it!=clusterList.end();++cluster_it)
        {
          Float_t pos = (*cluster_it) -> GetClusterPosition();
	  Float_t adc = (*cluster_it) -> GetClusterADCs(); 
	  hc1 -> Fill(pos, adc);
        }
      }
    }
  }
  return hc1;
}

void GEMZeroHitDecoder::DeleteClustersInPlaneMap() {
  map < TString, list <GEMCluster *> >::const_iterator itr ;
  for (itr = fListOfClustersZeroFromPlane.begin(); itr != fListOfClustersZeroFromPlane.end(); itr++) {
    list <GEMCluster *> listOfCluster = (*itr).second ;
    listOfCluster.clear() ;
  }
  fListOfClustersZeroFromPlane.clear() ;
}

void GEMZeroHitDecoder::GetClusterHyCalCutMode(vector<GEMClusterStruct> &gem1, vector<GEMClusterStruct> &gem2)
{
  list<GEMCluster*> cluster_x1 = fListOfClustersZeroFromPlane["pRadGEM1X"];
  list<GEMCluster*> cluster_y1 = fListOfClustersZeroFromPlane["pRadGEM1Y"];
  list<GEMCluster*> cluster_x2 = fListOfClustersZeroFromPlane["pRadGEM2X"];
  list<GEMCluster*> cluster_y2 = fListOfClustersZeroFromPlane["pRadGEM2Y"];

  int s1 = cluster_x1.size();
  int s2 = cluster_y1.size();  
  int nbCluster1 = (s1<s2)?s1:s2;
  s1 = cluster_x2.size();
  s2 = cluster_y2.size();
  int nbCluster2 = (s1<s2)?s1:s2;

  /*
   * Do the coordinate convert here
   * Convert GEM coordinate to HyCal coordinate
   *
   * **************************************
   *    Move GEM coordinate to Beam Hole center
   *    overlapping area: hole diameter (44mm)
   *    origin shift: 550.4/2 - (44-pitch)/2 = 253.2
   *
   *    gem1: x = x-253.2; y = y
   *    gem2: x = -x+253.2; y=-y
   *    
   *    right-hand coordinate, HyCal Y axis
   *    must pointing downward
   *
   *    beam downstream is z axis direction
   * **************************************
   *
   */
  double O_Transfer = 253.2;
  double OverlapLength = 44;
  double z_gem1 = 5300; //mm
  double z_gem2 = 5260; //mm

  // cutting edge
  float edge1 = 0;
  float edge2 = 0;

  double xoffset = -0.3618;
  double yoffset = 0.1792;
  
  //the above offsets are from the projection on GEM2
  //real GEM1 offsets should be projected back
  xoffset = xoffset*z_gem1/z_gem2;
  yoffset = yoffset*z_gem1/z_gem2;

  if(nbCluster1>0)
  {
    list<GEMCluster*>::iterator itx = cluster_x1.begin();
    list<GEMCluster*>::iterator ity = cluster_y1.begin();
    for(int i = 0;i<nbCluster1;i++)
    {
      /*
       * xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
       * overlapping area: according to frame design, overlapping area is in fact 44mm
       * xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
       */
      if(((*itx)->GetClusterPosition() -xoffset -O_Transfer) <= edge1) 
      // remove overlapping area on GEM1, use the corresponding area on GEM2
      // do not use edge as the cut line. choose some arbitrary line and find the corresponding line on GEM2
      {
        float c_x = (*itx)->GetClusterADCs();
	float c_y = (*ity)->GetClusterADCs();
        float x = (*itx++)->GetClusterPosition() -O_Transfer - xoffset;
	float y = (*ity++)->GetClusterPosition() -yoffset; 
        gem1.push_back( GEMClusterStruct(x, y, c_x, c_y) ) ;
      }
    }
  }

   if(nbCluster2>0)
  {
    list<GEMCluster*>::iterator itx2 = cluster_x2.begin();
    list<GEMCluster*>::iterator ity2 = cluster_y2.begin();
    for(int i = 0;i<nbCluster2;i++)
    {
      if( ( O_Transfer - ((*itx2)->GetClusterPosition())  ) > edge2)
      {
        float c_x = (*itx2)->GetClusterADCs();
	float c_y = (*ity2)->GetClusterADCs();
        float x =  O_Transfer - ((*itx2++)->GetClusterPosition());
	float y =  -(*ity2++)->GetClusterPosition() ; 
	gem2.push_back(GEMClusterStruct(x, y, c_x, c_y));
      }
    }
  }

}

void GEMZeroHitDecoder::GetClusterHyCalPlusMode(vector<GEMClusterStruct> &gem1, vector<GEMClusterStruct> &gem2)
{
  list<GEMCluster*> cluster_x1 = fListOfClustersZeroFromPlane["pRadGEM1X"];
  list<GEMCluster*> cluster_y1 = fListOfClustersZeroFromPlane["pRadGEM1Y"];
  list<GEMCluster*> cluster_x2 = fListOfClustersZeroFromPlane["pRadGEM2X"];
  list<GEMCluster*> cluster_y2 = fListOfClustersZeroFromPlane["pRadGEM2Y"];

  /*
   * Do the coordinate convert here
   * Convert GEM coordinate to HyCal coordinate
   *
   * **************************************
   *    Move GEM coordinate to HyCal Hole center
   *    overlapping area: hole diameter (44mm)
   *    origin shift: 550.4/2 - (44-pitch)/2 = 253.2
   *
   *    gem1: x = x-253.2; y = y
   *    gem2: x = -x+253.2; y=-y
   *    
   *    right-hand coordinate, HyCal Y axis
   *    must pointing downward
   *
   *    beam downstream is z axis direction
   * **************************************
   *
   */
  double O_Transfer = 253.2;
  double OverlapLength = 44;
  double z_gem1 = 5300; //mm
  double z_gem2 = 5260; //mm

  // offset from data
  double xoffset = -0.3618;
  double yoffset = 0.1792;

  //the above offsets are from the projection on GEM2
  //real GEM1 offsets should be projected back
  xoffset = xoffset*z_gem1/z_gem2;
  yoffset = yoffset*z_gem1/z_gem2;

  for( auto &i : cluster_x1)
  {
      for(auto &j : cluster_y1)
      {
          float c_x = i->GetClusterADCs();
	  float c_y = j->GetClusterADCs();
	  float x = i->GetClusterPosition() - O_Transfer-xoffset;
	  float y = j->GetClusterPosition() - yoffset;
	  gem1.push_back(GEMClusterStruct(x,y,c_x,c_y));
      }
  }
  for( auto &i : cluster_x2)
  {
      for(auto &j : cluster_y2)
      {
          float c_x = i->GetClusterADCs();
	  float c_y = j->GetClusterADCs();
	  float x = O_Transfer - i->GetClusterPosition();
	  float y = - j->GetClusterPosition() ;
	  gem2.push_back(GEMClusterStruct(x,y,c_x,c_y));
      }
  }

}


void GEMZeroHitDecoder::GetClusterBeamLine(vector<GEMClusterStruct> &gem1, vector<GEMClusterStruct> &gem2)
{
  list<GEMCluster*> cluster_x1 = fListOfClustersZeroFromPlane["pRadGEM1X"];
  list<GEMCluster*> cluster_y1 = fListOfClustersZeroFromPlane["pRadGEM1Y"];
  list<GEMCluster*> cluster_x2 = fListOfClustersZeroFromPlane["pRadGEM2X"];
  list<GEMCluster*> cluster_y2 = fListOfClustersZeroFromPlane["pRadGEM2Y"];

  int s1 = cluster_x1.size();
  int s2 = cluster_y1.size();  
  int nbCluster1 = (s1<s2)?s1:s2;
  s1 = cluster_x2.size();
  s2 = cluster_y2.size();
  int nbCluster2 = (s1<s2)?s1:s2;

  /*
   * Do the coordinate convert here
   * Convert GEM coordinate to HyCal coordinate
   *
   * **************************************
   *    Move GEM coordinate to Beam Hole center
   *    overlapping area: hole diameter (44mm)
   *    origin shift: 550.4/2 - (44-pitch)/2 = 253.2
   *
   *    gem1: x = x-253.2; y = y
   *    gem2: x = -x+253.2; y=-y
   *    
   *    right-hand coordinate, HyCal Y axis
   *    must pointing downward
   *
   *    beam downstream is z axis direction
   * **************************************
   *
   */
  double O_Transfer = 253.2;
  double OverlapLength = 44;
 
  // cutting edge
  float edge1 = 0;
  float edge2 = 0;
  double z_gem1 = 5300; //mm
  double z_gem2 = 5260; //mm

  double xoffset_gem = -0.3618;
  double yoffset_gem = 0.1792;

  //the above offsets are from the projection on GEM2
  //real GEM1 offsets should be projected back
  xoffset_gem = xoffset_gem*z_gem1/z_gem2;
  yoffset_gem = yoffset_gem*z_gem1/z_gem2;


  double xoffset_beam = 1.631;
  double yoffset_beam = -0.366;

  if(nbCluster1>0)
  {
    list<GEMCluster*>::iterator itx = cluster_x1.begin();
    list<GEMCluster*>::iterator ity = cluster_y1.begin();
    for(int i = 0;i<nbCluster1;i++)
    {
      /*
       * xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
       * overlapping area: according to frame design, overlapping area is in fact 44mm
       * xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
       */
      if(((*itx)->GetClusterPosition() -xoffset_gem -O_Transfer) <= edge1) 
      // remove overlapping area on GEM1, use the corresponding area on GEM2
      // do not use edge as the cut line. choose some arbitrary line and find the corresponding line on GEM2
      {
        float c_x = (*itx)->GetClusterADCs();
	float c_y = (*ity)->GetClusterADCs();
        float x = (*itx++)->GetClusterPosition() -O_Transfer - xoffset_gem - xoffset_beam ;
	float y =  (*ity++)->GetClusterPosition() -yoffset_gem - yoffset_beam;
	gem1.push_back(GEMClusterStruct(x, y, c_x, c_y));
      }
    }
  }

   if(nbCluster2>0)
  {
    list<GEMCluster*>::iterator itx2 = cluster_x2.begin();
    list<GEMCluster*>::iterator ity2 = cluster_y2.begin();
    for(int i = 0;i<nbCluster2;i++)
    {
      if( ( O_Transfer - ((*itx2)->GetClusterPosition())  ) < edge2) continue;
      {
        float c_x = (*itx2)->GetClusterADCs();
	float c_y = (*ity2)->GetClusterADCs();
        float x =  O_Transfer - ((*itx2++)->GetClusterPosition()) - xoffset_beam;
	float y = -(*ity2++)->GetClusterPosition()  - yoffset_beam;
	gem2.push_back(GEMClusterStruct(x, y, c_x, c_y) );
      }
    }
  }

}

void GEMZeroHitDecoder::GetClusterGEM(vector<GEMClusterStruct> &gem1, vector<GEMClusterStruct> &gem2)
{
  list<GEMCluster*> cluster_x1 = fListOfClustersZeroFromPlane["pRadGEM1X"];
  list<GEMCluster*> cluster_y1 = fListOfClustersZeroFromPlane["pRadGEM1Y"];
  list<GEMCluster*> cluster_x2 = fListOfClustersZeroFromPlane["pRadGEM2X"];
  list<GEMCluster*> cluster_y2 = fListOfClustersZeroFromPlane["pRadGEM2Y"];

  int s1 = cluster_x1.size();
  int s2 = cluster_y1.size();  
  int nbCluster1 = (s1<s2)?s1:s2;
  s1 = cluster_x2.size();
  s2 = cluster_y2.size();
  int nbCluster2 = (s1<s2)?s1:s2;

  if(nbCluster1>0)
  {
    list<GEMCluster*>::iterator itx = cluster_x1.begin();
    list<GEMCluster*>::iterator ity = cluster_y1.begin();
    for(int i = 0;i<nbCluster1;i++)
    {
      float c_x = (*itx)->GetClusterADCs();
      float c_y = (*ity)->GetClusterADCs();
      float x = (*itx++)->GetClusterPosition();
      float y = (*ity++)->GetClusterPosition();
      gem1.push_back(GEMClusterStruct(x, y, c_x, c_y));
    }
  }

   if(nbCluster2>0)
  {
    list<GEMCluster*>::iterator itx2 = cluster_x2.begin();
    list<GEMCluster*>::iterator ity2 = cluster_y2.begin();
    for(int i = 0;i<nbCluster2;i++)
    {
      float c_x = (*itx2)->GetClusterADCs();
      float c_y = (*ity2)->GetClusterADCs();
      float x =  (*itx2++)->GetClusterPosition(); 
      float y =  (*ity2++)->GetClusterPosition(); 
      gem2.push_back( GEMClusterStruct(x, y, c_x, c_y) );
    }
  }

}

void GEMZeroHitDecoder::FillHistos(TH1F* hNbClusterPerPlaneX[], TH1F* hNbClusterPerPlaneY[], TH1F* hClusterDistX[], TH1F* hClusterDistY[])
{
  list<GEMCluster*> cluster_x1 = fListOfClustersZeroFromPlane["pRadGEM1X"];
  list<GEMCluster*> cluster_y1 = fListOfClustersZeroFromPlane["pRadGEM1Y"];
  list<GEMCluster*> cluster_x2 = fListOfClustersZeroFromPlane["pRadGEM2X"];
  list<GEMCluster*> cluster_y2 = fListOfClustersZeroFromPlane["pRadGEM2Y"];

  int s1 = cluster_x1.size();
  hNbClusterPerPlaneX[0]->Fill(s1);
  int s2 = cluster_y1.size();
  hNbClusterPerPlaneY[0]->Fill(s2);
  for(int i=0;i<s1;i++)
  {
    list<GEMCluster*>::iterator itx = cluster_x1.begin();
    hClusterDistX[0] -> Fill( ( *itx++)->GetClusterPosition() );
  }
  for(int i=0;i<s2;i++)
  {
    list<GEMCluster*>::iterator ity = cluster_y1.begin();
    hClusterDistY[0]->Fill( (*ity++)->GetClusterPosition() );
  }

  s1 = cluster_x2.size();
  hNbClusterPerPlaneX[1] -> Fill(s1);
  s2 = cluster_y2.size();
  hNbClusterPerPlaneY[1] -> Fill(s2);
  for(int i=0;i<s1;i++)
  {
    list<GEMCluster*>::iterator itx = cluster_x2.begin();
    hClusterDistX[1] -> Fill( ( *itx++)->GetClusterPosition() );
  }
  for(int i=0;i<s2;i++)
  {
    list<GEMCluster*>::iterator ity = cluster_y2.begin();
    hClusterDistY[1]->Fill( (*ity++)->GetClusterPosition() );
  }
}

