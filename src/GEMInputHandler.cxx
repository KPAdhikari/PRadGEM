#include "GEMInputHandler.h"

//to test gem online show
#include "GEMOnlineHitDecoder.h"
//to test zero suppression
#include "GEMZeroHitDecoder.h"

using namespace std;
using namespace evio;

//#define CHECKRAWDATA
#define HITCONSTRUCTION

//_______________________________________________
GEMInputHandler::GEMInputHandler(string str)
{
  configure.LoadConfigure();
  filename = str;
  //file.open(filename.c_str(), ios::in|ios::binary);
  fRawDecoder = 0;

  mAPVRawHistos.clear();
  mAPVRawTSs.clear();
  vSRSSingleEventData.clear();
  vROCZeroSupData.clear();
  FECs.clear();

  PRDMapping *mapping = PRDMapping::GetInstance();
  FECs = mapping->GetBankIDSet();

  ped = new GEMPedestal(configure.GetLoadPedPath());
  ped->LoadPedestal();
  //for debug
  cout<<"FECs Found:  ";
  set<int>::iterator it;
  for(it=FECs.begin(); it!=FECs.end(); ++it)
  {
    cout<<(*it)<<"  ";
  }
  cout<<endl;
}

//______________________________________________
GEMInputHandler::~GEMInputHandler()
{
  mAPVRawHistos.clear();
  mAPVRawTSs.clear();
  ped->Delete();
}

//________________________________________________
int GEMInputHandler::ProcessAllEvents(int evtID)
{
  int entry = 0;
  try{

    evioFileChannel chan(filename.c_str(), "r");

    chan.open();

    while(chan.read())
    {
      //vector<uint32_t> srsSingleEvent;
      vSRSSingleEventData.clear();
      vROCZeroSupData.clear();

      evioDOMTree event(chan);
      evioDOMNodeListP fecEventList = event.getNodeList( isLeaf() );
      cout<<"total number of all banks: "<<fecEventList->size()<<endl;
      //int n_bank = fecEventList->size();
      //if(n_bank)

      evioDOMNodeList::iterator iter;
      for(iter=fecEventList->begin(); iter!=fecEventList->end(); ++iter)
      {
       // cout<<"bank #:"<< (*iter)->tag<<endl;

	//( ( (*iter)->tag == 5) || ( (*iter)->tag == 6  ) || ( (*iter)->tag ==7  ) || ((*iter)->tag ==8) ) 
	//if( FECs.find( ((*iter)->tag)-9 ) != FECs.end()  )
        
	if( ( (*iter)->tag == 57631) && ( ((int)((*iter)->num)==5) || ((int)((*iter)->num)==6) || ((int)((*iter)->num)==7) || ((int)((*iter)->num)==8) || ((int)((*iter)->num)==9) || ((int)((*iter)->num)==10) || ((int)((*iter)->num)==11) || ((int)((*iter)->num)==12) ) ) // new bank
	{
          //for debug
          //cout<<"GEMInputHandler() Processing FEC: "<< ((*iter)->tag)-5<<endl;
          vector<uint32_t> *vec = (*iter)->getVector<uint32_t>();
	  if(vec!=NULL  && vec->size()>0)
	  {
	    //for debug
	    //for(unsigned int i=0;i<vec->size();i++)
	    //{
	      //cout<<hex<<(*vec)[i]<<endl;
	      //cout<<(vec)->size()<<endl;
	    //}
	    //
	    //vec->pop_back(); //remove the last 0xfafafafa word
	    //int fec_id =( (*iter)->tag )-5; // one fec occupys one bank; one to one mapping
	    //cout<<"fec_id   bank:  "<<fec_id<<endl;
	    //vec->push_back(htonl(fec_id));
	    vSRSSingleEventData.reserve(vSRSSingleEventData.size() + vec->size() );
	    vSRSSingleEventData.insert(vSRSSingleEventData.end(), vec->begin(), vec->end() );
	  }
	  else
	  {
	    cout<<"found NULL contents in fec.."<<endl;
	  }
	}
        
	
        // ROC Level Zero Suppression Data
        // ------------------------------------------------------------------------------------------------
        //( ( (*iter)->tag == 5) || ( (*iter)->tag == 6  ) || ( (*iter)->tag ==7  ) || ((*iter)->tag ==8) ) 
	//if( FECs.find( ((*iter)->tag)-9 ) != FECs.end()  )
	if( ((*iter)->tag == 57631) && ((int)((*iter)->num)==99) ) // new bank
	{
          //for debug
          //cout<<"GEMInputHandler() Processing FEC: "<< ((*iter)->tag)-5<<endl;
          vector<uint32_t> *vec = (*iter)->getVector<uint32_t>();
	  if(vec!=NULL  && vec->size()>0)
	  {
	    //for debug
	    //for(unsigned int i=0;i<vec->size();i++)
	    //{
	      //cout<<hex<<(*vec)[i]<<endl;
	      //cout<<(vec)->size()<<endl;
	    //}
	    //
	    //vec->pop_back(); //remove the last 0xfafafafa word
	    //int fec_id =( (*iter)->tag )-5; // one fec occupys one bank; one to one mapping
	    //cout<<"fec_id   bank:  "<<fec_id<<endl;
	    //vec->push_back(htonl(fec_id));
	    vROCZeroSupData.reserve(vROCZeroSupData.size() + vec->size() );
	    vROCZeroSupData.insert(vROCZeroSupData.end(), vec->begin(), vec->end() );
	  }
	  else
	  {
	    cout<<"ROC Zero Sup: found NULL contents in fec.."<<endl;
	  }
	}
        
      }
      
      ++entry;
      cout<<"entry: "<<entry<<endl;
      cout<<"SRS Event Size [uint_32]:"<<vSRSSingleEventData.size()<<endl;

      if (vSRSSingleEventData.size() == 0 ) continue; // if no srs event found, go to next event
      
      //vSRSSingleEventData.push_back(0xfafafafa); // add back srs event trailor
      GEMRawDecoder raw_decoder(vSRSSingleEventData);

      //==============================================================================//
      //test online hit show
#ifdef HITCONSTRUCTION      
      int size = vSRSSingleEventData.size();
      uint32_t buf[size];
      for(int i=0;i<size;i++)
      {
        buf[i]=vSRSSingleEventData[i];
      }
      int roc_size = vROCZeroSupData.size();
      uint32_t roc_buf[roc_size];
      for(int i=0;i<roc_size;i++)
      {
        roc_buf[i] = vROCZeroSupData[i];
      }
      GEMOnlineHitDecoder online_hit(buf, size, ped);
      // non-zero suppression
      TH1F* ho1x = online_hit.GetHit("pRadGEM1X");
      TH1F* ho1y = online_hit.GetHit("pRadGEM1Y");
      TH1F* ho2x = online_hit.GetHit("pRadGEM2X");
      TH1F* ho2y = online_hit.GetHit("pRadGEM2Y");
      // offline zero suppression
      TH1F* hz1x = online_hit.GetCleanHit("pRadGEM1X");
      TH1F* hz1y = online_hit.GetCleanHit("pRadGEM1Y");
      TH1F* hz2x = online_hit.GetCleanHit("pRadGEM2X");
      TH1F* hz2y = online_hit.GetCleanHit("pRadGEM2Y");
      TH1F* hc = online_hit.GetCluster("pRadGEM1X");
      // roc level online zero suppression
      GEMZeroHitDecoder zero_hit(roc_buf, roc_size);
      TH1F* hr1x = zero_hit.GetZeroHit("pRadGEM1X");
      TH1F* hr1y = zero_hit.GetZeroHit("pRadGEM1Y");
      TH1F* hr2x = zero_hit.GetZeroHit("pRadGEM2X");
      TH1F* hr2y = zero_hit.GetZeroHit("pRadGEM2Y");

      // plot histograms
      int total = hz1x->GetEntries() + hz1y->GetEntries() + hz2x->GetEntries() + hz2y->GetEntries();
      cout<<"GEMInputHandler offline hits found: "<<total<<endl;
      TCanvas *ca = new TCanvas("ca", "ca", 2000, 1000);
      ca->Divide(2,2);
      ca->cd(1);
      ho1x->Draw();
      ca->cd(2);
      ho1y->Draw();
      ca->cd(3);
      ho2x->Draw();
      ca->cd(4);
      ho2y->Draw();
      ca->Update();

/*
      TCanvas *ca = new TCanvas("ca", "ca", 2000, 1000);
      ca->Divide(2,2);
      ca->cd(1);
      hr1x->Draw();
      ca->cd(2);
      hr1y->Draw();
      ca->cd(3);
      hr2x->Draw();
      ca->cd(4);
      hr2y->Draw();
      ca->Update();
*/
      /*
      // 1x
      TCanvas *ccc = new TCanvas("ccc", "1x", 1000, 1000);
      TCanvas *ccc1 = new TCanvas("ccc1", "1y", 1000, 1000);
      TCanvas *ccc2 = new TCanvas("ccc2", "2x", 1000, 1000);
      TCanvas *ccc3 = new TCanvas("ccc3", "2y", 1000, 1000);
      ccc->Divide(1,3);
      ccc1->Divide(1,3);
      ccc2->Divide(1,3);
      ccc3->Divide(1,3);

      ccc->cd(1);
      //h->Print("all");
      ho1x->Draw();
      ccc->cd(2);
      hz1x->Draw();
      ccc->cd(3);
      hr1x->Draw();
      ccc->Update();

      ccc1->cd(1);
      //h->Print("all");
      ho1y->Draw();
      ccc1->cd(2);
      hz1y->Draw();
      ccc1->cd(3);
      hr1y->Draw();
      ccc1->Update();

      ccc2->cd(1);
      //h->Print("all");
      ho2x->Draw();
      ccc2->cd(2);
      hz2x->Draw();
      ccc2->cd(3);
      hr2x->Draw();
      ccc2->Update();

      ccc3->cd(1);
      //h->Print("all");
      ho2y->Draw();
      ccc3->cd(2);
      hz2y->Draw();
      ccc3->cd(3);
      hr2y->Draw();
      ccc3->Update();
      */

      // next event
      getchar();
      ho1x->Delete();
      ho1y->Delete();
      ho2x->Delete();
      ho2y->Delete();

      hz1x->Delete();
      hz1y->Delete();
      hz2x->Delete();
      hz2y->Delete();
      hc->Delete();

      hr1x->Delete();
      hr1y->Delete();
      hr2x->Delete();
      hr2y->Delete();

      //delete ccc;

      vector<GEMClusterStruct> gem1, gem2;
      online_hit.GetClusterGEM(gem1, gem2);
      if(gem1.size() > 0)
      {
        for(int i=0;i<gem1.size();i++)
        {
          cout<<"x1: "<<gem1[i].x<<"  y1: "<<gem1[i].y<<endl;
        }
      }
      if(gem2.size()>0)
      {
        for(int i=0;i<gem2.size();i++)
	{
	  cout<<"x2: "<<gem2[i].x<<"  gem2: "<<gem2[i].y<<endl;
	}
      }
#endif
      //==============================================================================//
      //test zero suppression hit show
      /*
      uint32_t zbuf[12];
      zbuf[0] = 0x10008008;
      zbuf[1] = 0x10009008;
      zbuf[2] = 0x1000a008;
      zbuf[3] = 0x48008008;
      zbuf[4] = 0x48009008;
      zbuf[5] = 0x4800a008;
      zbuf[6] = 0x80008008;
      zbuf[7] = 0x80009008;
      zbuf[8] = 0x8000a008;
      zbuf[9] = 0xc4008008;
      zbuf[10] = 0xc4009008;
      zbuf[11] = 0xc400a008;

      GEMZeroHitDecoder zero_decoder(zbuf, 12);
      TH1F* h4 = zero_decoder.GetZeroHit("pRadGEM1Y");
      TH1F* h4h = zero_decoder.GetCluster("pRadGEM1Y");
      TCanvas *c4 = new TCanvas("c4", "c4", 800,800);
      c4->Divide(1,2);
      c4->cd(1);
      h4->Draw();
      c4->cd(2);
      h4h->Draw();
      c4->Update();
      getchar();
      h4->Delete();
      h4h->Delete();
      delete c4;
      */

      /*
      //test GEMPedestal
      GEMPedestal * ped;
      ped = new GEMPedestal("./pedestal.root");
      ped->LoadPedestal();
      ped->Delete();
      */
      //==============================================================================//

      //Draw raw histos
      //
      //debug
      //GEMRawPedestal pedestal(raw_decoder.GetDecoded());
      //pedestal.PrintEventPedestal();
      //end debug

#ifdef CHECKRAWDATA
      vector<TH1F*> RawAPVs;
      TCanvas *c = new TCanvas("c", "APV Raw Signal", 10, 10, 2400, 1000);
      c->Divide(9,8);

     
      mAPVRawHistos = raw_decoder.GetAPVRawHisto();
      map<int, map<int, TH1F*> >::iterator it;
      RawAPVs.clear();
      TFile *f = new TFile("histos.root", "recreate");
      for(it=mAPVRawHistos.begin(); it!=mAPVRawHistos.end();++it)
      {
        map<int, TH1F*> temp = it->second;
        map<int, TH1F*>::iterator itt;
        for(itt=temp.begin();itt!=temp.end();++itt)
        {
          TH1F* h = itt->second;
          RawAPVs.push_back(h);
        }
      }
      int nn = RawAPVs.size();
      for(int i=0; i<nn;i++)
      {
       c->cd(i+1);
       RawAPVs[i]->Draw();
       RawAPVs[i]->Write();
       //if(RawAPVs[i]->Get
       //cout<<RawAPVs[i]->GetName()<<endl;
      }
      c->Update();
      c->Print("signal.pdf");
      c->Print("signal.png");
      getchar();
      f->Close();
      c -> Delete(); 
      delete c;
#endif

    }

    chan.close();

  } catch (evioException e) {
    cerr <<endl <<e.toString() <<endl <<endl;
    exit(EXIT_FAILURE);
  }

  //exit(EXIT_SUCCESS);
  return entry;
}
