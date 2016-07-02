/*******************************************************************************************
 * Xinzhan Bai
 *
 * Usage:
 *
 * Input:  1), A buffer contains purely SRS data, and the buffer size in [int] unit
 *
 *         2), Or, SRS data filled in a vector<int>
 *
 * Output: 1), map<fec_id, map<ch_id, vector<int> > >, vector: adc values (all time samples)
 *                map<int, map<int, vector<int> > > GetDecoded();
 *
 *         2), map<fec_id, map<ch_id, TH1F* > > , TH1F*: adc values filled in histograms
 *             No need to worry memory leakage, this class owns every object it produced.
 *                map<int, map<int, TH1F* > > GetAPVRawHisto();
 * *****************************************************************************************/

#ifndef __RAW_DECODER_H__
#define __RAW_DECODER_H__

#include <map>
#include <vector>
#include <algorithm>

#include <TH1F.h>

#include "PRDMapping.h"

using namespace std;

class GEMRawDecoder
{
public:
  GEMRawDecoder( vector<int> );
  GEMRawDecoder( unsigned int *, int);
  ~GEMRawDecoder();

  void Decode();
  map<int, map<int, vector<int> > > GetDecoded();
  map<int, map<int, TH1F* > > GetAPVRawHisto();
  void Word32ToWord16(unsigned int *, unsigned int*, unsigned int*);

  int IsAdcChannelActive(int, int );

private:
  unsigned int * buf;
  int fBuf;
  map<int, map<int, vector<int> > > mAPVRawSingleEvent;
  map<int, map<int, TH1F* > > mAPVRawHisto;

  PRDMapping* mapping;

  vector<int> vActiveAdcChannels;
};

#endif
