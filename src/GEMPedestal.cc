#include "GEMPedestal.h"
#include "GEMRawDecoder.h"
#include "GEMMapping.h"
#include "GEMRawPedestal.h"
#include "GEMConfigure.h"
#include <stdio.h>
#include <TCanvas.h>
#include <iostream>

using namespace std;

GEMPedestal::GEMPedestal()
{
    config = new GEMConfigure();
    config->LoadConfigure();
    eventPedestal = new GEMRawPedestal();
    mapping = GEMMapping::GetInstance();
    NCH = 128;
    nNbofAPVs = mapping->GetNbOfAPVs();

    FECs.clear();

    FECs = mapping->GetBankIDSet();
    SetFile();
}

GEMPedestal::~GEMPedestal()
{
    FECs.clear();
}

void GEMPedestal::SetFile()
{
    filename = config->GetLoadPedPath();
}

void GEMPedestal::Delete()
{
    int N = vStripOffsetHistos.size();
    for(int i=0;i<N;i++)
    {
	delete vStripOffsetHistos[i];
	delete vStripNoiseHistos[i];
    }
    int M = vApvPedestalOffset.size();
    for(int i=0;i<M;i++)
    {
	delete vApvPedestalOffset[i];
	delete vApvPedestalNoise[i];
    }
}

void GEMPedestal::BookHistos()
{
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

int GEMPedestal::AccumulateEvent(int evtID, unordered_map<int, vector<int> > & mAPVRawTSs)
{
    int entry = evtID;
    eventPedestal -> ComputeEventPedestal(mAPVRawTSs);

    for(auto &i: mAPVRawTSs)
    {
        int apvid = i.first;
	int apvKey = mapping->GetAPVNoFromID(apvid);
	for(int chNo=0;chNo<NCH;chNo++)
	{
	    vStripNoiseHistos[apvKey*NCH + chNo] ->  Fill( eventPedestal->GetStripNoise(apvid, chNo));
	    vStripOffsetHistos[apvKey*NCH + chNo] -> Fill( eventPedestal->GetStripOffset(apvid, chNo));
	}
    }
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
		cout<<"GEMPedestal::ComputePedestal: Error: Unrecongnized plane name..."
		    <<endl;
	    }
	}
    }
}

void GEMPedestal::SavePedestalFile()
{
    TFile *file = new TFile(config->GetSavePedPath().c_str(), "recreate");
    int NN = config->GetNumEvtForPed();
 
    cout<<"Number of Events for Pedestal: "<<NN<<endl;
    ComputePedestal();

    for(int apvKey=0;apvKey<nNbofAPVs;apvKey++)
    {
	vApvPedestalOffset[apvKey]->Write();
	vApvPedestalNoise[apvKey]->Write();
    }
    hAllStripNoise->Write();
    hAllXStripNoise->Write();
    hAllYStripNoise->Write();
    file->Write();
}

void GEMPedestal::LoadPedestal()
{
    cout<<" loading pedestal from: "
        <<filename<<endl;
    TFile *_file = new TFile( filename.c_str(), "READ" );
    if(_file->IsZombie() )
    {
	cout<<"#### Cannot Load pedestal file... ####"<<endl;
	return;
    }
    Int_t nAPVs = mapping->GetNbOfAPVs();
    for(int i=0;i<nAPVs;i++)
    {
	stringstream out;
	out << i;
	TString outStr = out.str();
	TString offset = "offset_APV_#:" + outStr;
	TString noise =  "noise_APV_#:" + outStr;
	vApvPedestalOffset.push_back( (TH1F*) _file->Get( offset )  );
	vApvPedestalOffset[i]->SetDirectory(0);
	vApvPedestalNoise.push_back( (TH1F*) _file->Get( noise )  ); 
	vApvPedestalNoise[i]->SetDirectory(0);
    }

    // Cannot close file while pedestal histograms are still being used...
    _file -> Close();
}

vector<Float_t>  GEMPedestal::GetAPVNoises(Int_t apvid)
{
    vector<Float_t> noises;
    Int_t apvNo = mapping->GetAPVNoFromID(apvid);
    for(int i=0;i<NCH;i++)
    {
	noises.push_back( vApvPedestalNoise[apvNo]->GetBinContent(i)  );
    }
    return noises;
}

vector<Float_t>  GEMPedestal::GetAPVOffsets(Int_t apvid)
{
    vector<Float_t> offsets;
    Int_t apvNo = mapping->GetAPVNoFromID(apvid);
    for(int i=0;i<NCH;i++)
    {
	offsets.push_back( vApvPedestalOffset[apvNo]->GetBinContent(i)  );
    }
    return offsets;
}
