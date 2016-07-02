/************************************************
 *
 * Xinzhan Bai 
 * xb4zp@virginia.edu
 *
 ************************************************/

#include <arpa/inet.h>
#include <assert.h>
#include <utility>
#include "GEMRawDecoder.h"

#include <fstream>
#include <iostream>
#include <utility>

#include <TH1F.h>

using namespace std;

GEMRawDecoder::GEMRawDecoder(unsigned int * buffer, int n)
{ 
  mapping = PRDMapping::GetInstance();
  mAPVRawSingleEvent.clear();
  mAPVRawHisto.clear();
  vActiveAdcChannels.clear();

  fBuf = n;
  buf = new unsigned int[fBuf];
  for(int i=0;i<fBuf;i++)
  {
    buf[i] = buffer[i];
  }
  Decode();
}

GEMRawDecoder::GEMRawDecoder(vector<int> buffer)
{
  mapping = PRDMapping::GetInstance();
  mAPVRawSingleEvent.clear();
  mAPVRawHisto.clear();
  vActiveAdcChannels.clear();

  fBuf = buffer.size();
  buf = new unsigned int[fBuf];
  for(int i=0;i<fBuf;i++)
  {
    buf[i] = buffer[i];
  }
  if(fBuf <= 0)
  {
    cout<<"empty vector passed in..."<<endl;
    return;
  }
  Decode();
};

GEMRawDecoder::~GEMRawDecoder()
{
  //free buf
  delete[] buf;

  //clear maps
  map<int, map<int, vector<int> > >::iterator it;
  for(it=mAPVRawSingleEvent.begin(); it!=mAPVRawSingleEvent.end(); ++it)
  {
    map<int, vector<int> > temp = it->second;
    map<int, vector<int> >::iterator it1;
    for(it1=temp.begin(); it1 != temp.end(); ++it1)
    {
      it1->second.clear();
    } 
    it->second.clear();
  }
  mAPVRawSingleEvent.clear();

  //clear APV raw histos
  map<int, map<int, TH1F*> >::iterator it_apv;
  for(it_apv=mAPVRawHisto.begin(); it_apv!=mAPVRawHisto.end(); ++it_apv)
  {
    map<int, TH1F*> temp = it_apv->second;
    map<int, TH1F*>::iterator itt_apv;

    for(itt_apv=temp.begin(); itt_apv!=temp.end(); ++itt_apv)
    {
      TH1F *h = itt_apv->second;
      h->Delete();
    }
    temp.clear();
  }
  mAPVRawHisto.clear();
  vActiveAdcChannels.clear();
}

void GEMRawDecoder::Decode()
{
  unsigned int word32bit;
  unsigned int word16bit1;
  unsigned int word16bit2;

  int fecid;
  int adc_ch;

  vector<int> apv_event;
  vector<int> apv_margin;
  vector<int> fec_margin;

  //for debug
  //fstream fs;
  //fs.open("raw.txt",fstream::out|fstream::app);
  for(int i=0;i<fBuf;i++)
  {
    //endianness switch ...
    //buf[i] = ntohl(buf[i]);
    //fs<<hex<<buf[i]<<endl;
  }
  //cout<<"fBuf  :"<<fBuf<<endl;
  //fs<<"#################"<<endl;
 
  //for debug
  assert( ((buf[0] >> 8) & 0xffffff) != 0x414443  );

  fec_margin.clear();
  apv_margin.clear();
  
  //find fec margin
  //... to be implemented...
  //fecid=0;
  fec_margin.push_back(-1);
  for(int i=0;i<fBuf;i++)
  {
    if(buf[i]==0xfafafafa)
    {
      fec_margin.push_back(i);
      //cout<<dec<<i<<endl;
      //cout<<hex<<buf[i]<<"  "<<buf[i+3]<<endl;
    }
  }
  //cout<<dec<<fBuf<<endl;
  //cout<<"very end "<<hex<<buf[fBuf]<<endl;
  int nFecMargin = fec_margin.size();
  //cout<<"GEMRawDecoder: Nb of FECs:"<<nFecMargin-1<<endl;
  if(nFecMargin < 2)
  {
    cout<<"GEM Raw Decoder: No SRS Data Receved..."<<endl;
  }

  assert( buf[fec_margin[nFecMargin-1]]==0xfafafafa);

  for(int ii = 0;ii<=nFecMargin-2;ii++)
  {
    //find apv margin
    apv_margin.clear();
    for(int i=fec_margin[ii]+1;i<=fec_margin[ii+1];i++)
    {
      if( ( (buf[i] >> 8) & 0xffffff ) == 0x414443 ) 
      {
        apv_margin.push_back(i-1);
      }
      if(buf[i]==0xfafafafa) apv_margin.push_back(i);
    }

    int nMargin = apv_margin.size();
    //cout<<"Number of APVs:   "<<nMargin-1<<endl;
    assert( buf[apv_margin[nMargin-1]] == 0xfafafafa  );

    for(int i=0;i<=nMargin-2;i++)
    {
      //fecid = buf[ margin[i] ] & 0xff;
      apv_event.clear();
      adc_ch = buf[ apv_margin[i]+1 ] & 0xff;
      
      //If use Bryan's setup, following three lines should be true
      fecid = ( buf[ apv_margin[i]+2 ] >> 16 )&0xff;
      //int fecid_bank = buf[ fec_margin[ii+1] + 1 ];
      //assert(fecid == fecid_bank);

      for(int k=apv_margin[i]+3; k<apv_margin[i+1]; k++)
      {
        word32bit = buf[k];
      
        Word32ToWord16(&word32bit, &word16bit1, &word16bit2);
        //for debug
        //if(word16bit1>3500 || word16bit2>3500) continue;
        apv_event.push_back( (int)word16bit1 );
        apv_event.push_back( (int)word16bit2 );
      }
      //fec_event.insert( pair<int, vector<int> >(adc_ch, apv_event)  );
      //mAPVRawSingleEvent.insert( pair<int, map<int, vector<int> > > (fecid, fec_event) );
      //mAPVRawSingleEvent[fecid].insert( pair<int, vector<int> >(adc_ch, apv_event)  );
      //cout<<"fec id: "<<fecid<<"  adc_ch:"<<adc_ch<<endl;

      /////////////////////////////////////////////////////////////////////////////////////
      //check first 100 words, to make sure no data lost in case using the wrong mapping //
      //by checking if this channel has APV header                                       //
      //int apvid = fecid << 4 | adc_ch;                                                 //
      //int header = mapping->GetAPVHeaderLevelFromID(apvid);                            //
      //cout<<"header: "<<header<<endl;                                                  //
      /////////////////////////////////////////////////////////////////////////////////////
      int header = 1500;
      if(! IsAdcChannelActive(fecid, adc_ch) )
      {
        //cout<<"header: "<<header<<endl;
        for(int i=0;i<100;i++)
	{
	  if(apv_event[i] < header)
	  {
	    i++;
	    if(apv_event[i] < header)
	    {
	      i++;
	      if(apv_event[i] < header )
	      {
	        i++;
		cout<<"## Deocder WARNNING: ##  Found meaningful data in INACTIVE channels..."<<endl;
		cout<<"                     ##  FEC: "<<fecid<<"  Channel:  "<<adc_ch<<endl;
		cout<<"                     ##  check you mapping file, or check the raw data..."<<endl;
		break;
	      }
	    }
	  }
	}
        //cout<<"INACTIVE Channel: FEC: "<<fecid<<"  Channel:  "<<adc_ch<<endl;
	continue;
      }

      mAPVRawSingleEvent[fecid][adc_ch] = apv_event;
    }
  }  
}

map<int, map<int, vector<int> > >  GEMRawDecoder::GetDecoded()
{
  return mAPVRawSingleEvent;
}

map<int, map<int, TH1F* > > GEMRawDecoder::GetAPVRawHisto()
{
  int fec_id=0;
  int adc_ch=0;

  map<int, TH1F*> ch_apv;
  map<int, map<int, vector<int> > >::iterator it;
  for(it = mAPVRawSingleEvent.begin(); it!=mAPVRawSingleEvent.end(); ++it)
  {
    fec_id = it->first;
    map<int, vector<int> > temp = it->second; 

    map<int, vector<int> >::iterator itt;
    for(itt=temp.begin(); itt!=temp.end(); ++itt)
    {
      adc_ch = itt->first;
      vector<int> adc_temp = itt->second;

      int N = adc_temp.size();
      TH1F* h = new TH1F(Form("fec_%d_ch_%d",fec_id, adc_ch), Form("fec_%d_ch_%d_raw_data",fec_id, adc_ch), N+128, 0, N+127);
      for(int i=0;i<N;i++)
      {
        h->Fill(i+1, (Float_t) adc_temp[i]);
      }
      mAPVRawHisto[fec_id][adc_ch] = h;
    }
  }
  return mAPVRawHisto;
}

void GEMRawDecoder::Word32ToWord16(unsigned int *word, unsigned int *word1, unsigned int *word2)
{
  unsigned int data1=0;
  unsigned int data2=0;
  unsigned int data3=0;
  unsigned int data4=0;

  data1 = ( (*word)>>24 ) & 0xff;
  data2 = ( (*word)>>16 ) & 0xff;
  data3 = ( (*word)>>8  ) & 0xff;
  data4 = (*word) & 0xff;

  (*word1) = (data2 << 8) | data1;
  (*word2) = (data4 << 8) | data3;
}

int GEMRawDecoder::IsAdcChannelActive(int fecid, int ch)
{
  int isActive = 0;
  vActiveAdcChannels = mapping->GetActiveADCchannels(fecid);
  if( find( vActiveAdcChannels.begin(), vActiveAdcChannels.end(), ch) != vActiveAdcChannels.end() )
  {
    isActive =1;
  }
  return isActive;
}
