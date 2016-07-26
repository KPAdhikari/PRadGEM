#include "GEMPedestal.h"

using namespace std;
using namespace evio;

GEMPedestal::GEMPedestal(string str)
{
  config.LoadConfigure();
  filename = str;
  pedestal_file = str;
  mapping = PRDMapping::GetInstance();
  NCH = 128;
  nNbofAPVs = mapping->GetNbOfAPVs();

  //mAPVRawHistos.clear();
  mAPVRawTSs.clear();
  vSRSSingleEventData.clear();
  FECs.clear();

  PRDMapping *mapping = PRDMapping::GetInstance();
  FECs = mapping->GetBankIDSet();
}

GEMPedestal::~GEMPedestal()
{
  //mAPVRawHistos.clear();
  vSRSSingleEventData.clear();
  FECs.clear();
  mAPVRawTSs.clear();
}

void GEMPedestal::Delete()
{
  int N = vStripOffsetHistos.size();
  for(int i=0;i<N;i++)
  {
    //vStripOffsetHistos[i]->Delete();
    delete vStripOffsetHistos[i];
    //vStripNoiseHistos[i]->Delete();
    delete vStripNoiseHistos[i];
  }
  int M = vApvPedestalOffset.size();
  for(int i=0;i<M;i++)
  {
    //vApvPedestalOffset[i]->Delete();
    delete vApvPedestalOffset[i];
    //vApvPedestalNoise[i]->Delete();
    delete vApvPedestalNoise[i];
  }
}

void GEMPedestal::BookHistos()
{
  //this function mainly book all needed histograms
  
  //book histograms for each strip
  for(int chNo = 0; chNo < NCH; chNo++)
  {
    for(int apvKey = 0; apvKey < nNbofAPVs; apvKey++)
    {
      stringstream out;
      stringstream apv;
      apv << apvKey;
      out << chNo;
      TString chNoStr = out.str();
      TString apvStr = apv.str();
      TString noise = "hNoise_" + apvStr+ "_" + chNoStr;
      TString offset = "hOffset_" +apvStr+ "_"+ chNoStr;
      vStripNoiseHistos.push_back( new TH1F(noise, noise, 8097, -4048, 4048) ); 
      vStripOffsetHistos.push_back( new TH1F(offset, offset, 8097, -4048, 4048)  );
    }
  }

  // book histograms for each APV
  for(int apvKey = 0; apvKey < nNbofAPVs; apvKey++)
  {
    stringstream out;
    out << apvKey;
    TString outStr = out.str();
    TString offset = "offset_APV_#:" + outStr;
    TString noise = "noise_APV_#:" + outStr;
    int fAPVID = mapping->GetAPVIDFromAPVNo(apvKey);
    int fFECID = (fAPVID>>4);
    int fADCCh = (fAPVID)&0xf;
    TString offset_title = offset;
    TString noise_title = noise;
    out.str("");
    out << fFECID;
    outStr = out.str();
    offset_title = offset_title+"_FEC_"+outStr;
    noise_title = noise_title+"_FEC_"+outStr;
    out.str("");
    out << fADCCh;
    outStr = out.str();
    offset_title = offset_title+"_ADCCh_"+outStr;
    noise_title = noise_title+"_ADCCh_"+outStr;
    vApvPedestalOffset.push_back( new TH1F(offset, offset_title, 128, -0.5, 127.5)  );
    vApvPedestalNoise.push_back(  new TH1F(noise, noise_title, 128, -0.5, 127.5)  );
  }

  // book histograms for overall distribution
  hAllStripNoise = new TH1F("hAllStripNoise", "Overall Noise Distribution", 1000, 0, 1000);
  hAllXStripNoise = new TH1F("hAllXStripNoise", "Overall X Direction Noise Distribution", 1000, 0, 1000);
  hAllYStripNoise = new TH1F("hAllYStripNoise", "Overall Y Direction Noise Distribution", 1000, 0, 1000);
  
}

int GEMPedestal::ProcessAllEvents(int evtID)
{
  int entry = 0;

  try{

    evioFileChannel chan(filename.c_str(), "r");

    chan.open();

    while(chan.read())
    {
      vSRSSingleEventData.clear();

      evioDOMTree event(chan);
      evioDOMNodeListP fecEventList = event.getNodeList( isLeaf() );
      int n_bank = fecEventList->size();

      evioDOMNodeList::iterator iter;
      for(iter=fecEventList->begin(); iter!=fecEventList->end(); ++iter)
      {

	if (  ((*iter)->tag == 57631)  && ( ((int)((*iter)->num)==5) || ((int)((*iter)->num)==6)|| ((int)((*iter)->num)==7) || ((int)((*iter)->num)==8) || ((int)((*iter)->num)==9)  ||((int)((*iter)->num)==10) ||((int)((*iter)->num)==11) ||((int)((*iter)->num)==12)  )  )
	{
          vector<uint32_t> *vec = (*iter)->getVector<uint32_t>();
	  if(vec!=NULL  && vec->size()>0)
	  {
	    vSRSSingleEventData.reserve(vSRSSingleEventData.size() + vec->size() );
	    vSRSSingleEventData.insert(vSRSSingleEventData.end(), vec->begin(), vec->end() );
	  }
	  else
	  {
	    cout<<"found NULL contents in fec.."<<endl;
	  }
	}
      }
      
      if (vSRSSingleEventData.size() == 0 ) continue; // if no srs event found, go to next event
       ++entry; //++nPedestalEntry;
      if( evtID!=-1 && entry > evtID) return entry;
      cout<<"entry: "<<entry<<endl;
      cout<<"GEMPedestal::SRS Event Size [uint_32]:"<<vSRSSingleEventData.size()<<endl;

     
      GEMRawDecoder raw_decoder(vSRSSingleEventData);

      mAPVRawTSs = raw_decoder.GetDecoded();
      GEMRawPedestal raw_pedestal(mAPVRawTSs);

      map<int, map<int, vector<int> > >::iterator it;
      for(it=mAPVRawTSs.begin(); it!=mAPVRawTSs.end();++it)
      {
        int fecid = it->first;
        map<int, vector<int> > fecdata = it->second;
        map<int, vector<int> >::iterator itt;
        for(itt=fecdata.begin();itt!=fecdata.end();++itt)
        {
	  int adc_ch = itt->first;

	  for(int chNo=0; chNo < NCH; chNo++)
	  {
	    int apvID = fecid << 4 | adc_ch;
	    int apvKey = mapping->GetAPVNoFromID(apvID);
	    vStripNoiseHistos[apvKey*NCH + chNo] ->  Fill( raw_pedestal.GetStripNoise(fecid, adc_ch, chNo)   );
	    vStripOffsetHistos[apvKey*NCH + chNo] -> Fill( raw_pedestal.GetStripOffset(fecid, adc_ch, chNo)  );
	  }
        }
      }

    }

    chan.close();

  } catch (evioException e) {
    cerr <<endl <<e.toString() <<endl <<endl;
    exit(EXIT_FAILURE);
  }

  return entry;
}

void GEMPedestal::ComputePedestal()
{
  for(int apvKey = 0; apvKey < nNbofAPVs; apvKey++)
  {
    for(int chNo = 0; chNo < NCH; chNo++)
    {
      Float_t offset = vStripOffsetHistos[apvKey*NCH + chNo] -> GetMean();
      Float_t noise  = vStripNoiseHistos[apvKey*NCH + chNo] -> GetRMS();
      vApvPedestalOffset[apvKey]->SetBinContent(chNo, offset);
      vApvPedestalNoise[apvKey]->SetBinContent(chNo, noise);
      hAllStripNoise->Fill(noise);
      //X and Y to be implemented
      Int_t fAPVID = mapping->GetAPVIDFromAPVNo(apvKey);
      TString plane = mapping->GetPlaneFromAPVID(fAPVID);
      if(plane.Contains("X"))
      {
        hAllXStripNoise->Fill(noise);
      }
      else if( plane.Contains("Y"))
      {
        hAllYStripNoise->Fill(noise);
      }
      else
      {
        cout<<"GEMPedestal::ComputePedestal: Error: Unrecongnized plane name..."<<endl;
      }
    }
  }
}

void GEMPedestal::SavePedestalFile()
{
  BookHistos();
  TFile *file = new TFile(config.GetSavePedPath().c_str(), "recreate");
  int NN = config.GetNumEvtForPed();
  
  cout<<"Number of Events for Pedestal: "<<NN<<endl;
  ProcessAllEvents(NN);
  ComputePedestal();

  for(int apvKey=0;apvKey<nNbofAPVs;apvKey++)
  {
    vApvPedestalOffset[apvKey]->Write();
    vApvPedestalNoise[apvKey]->Write();
  }
  hAllStripNoise->Write();
  hAllXStripNoise->Write();
  hAllYStripNoise->Write();
}

void GEMPedestal::LoadPedestal()
{
  TFile *file = new TFile( pedestal_file.c_str(), "READ" );
  if(file->IsZombie() )
  {
    cout<<"#### Cannot Load pedestal file... ####"<<endl;
    return;
  }
  //stringstream out;
  Int_t nAPVs = mapping->GetNbOfAPVs();
  for(int i=0;i<nAPVs;i++)
  {
    stringstream out;
    out << i;
    TString outStr = out.str();
    TString offset = "offset_APV_#:" + outStr;
    TString noise =  "noise_APV_#:" + outStr;
    vApvPedestalOffset.push_back( (TH1F*) file->Get( offset )  );
    vApvPedestalOffset[i]->SetDirectory(0);
    vApvPedestalNoise.push_back( (TH1F*) file->Get( noise )  ); 
    vApvPedestalNoise[i]->SetDirectory(0);
  }

  // Cannot close file while pedestal histograms are still being used...
  file -> Close();
}

vector<Float_t> GEMPedestal::GetAPVNoises(Int_t apvid)
{
  vector<Float_t> noises;
  Int_t apvNo = mapping->GetAPVNoFromID(apvid);
  for(int i=0;i<NCH;i++)
  {
    noises.push_back( vApvPedestalNoise[apvNo]->GetBinContent(i)  );
  }
  return noises;
}

vector<Float_t> GEMPedestal::GetAPVOffsets(Int_t apvid)
{
  vector<Float_t> offsets;
  Int_t apvNo = mapping->GetAPVNoFromID(apvid);
  for(int i=0;i<NCH;i++)
  {
    offsets.push_back( vApvPedestalOffset[apvNo]->GetBinContent(i)  );
  }
  return offsets;
}
