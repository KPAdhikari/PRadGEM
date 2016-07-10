#include <iostream>
#include <fstream>
#include <map>
#include <string>

#include <stdio.h> // for getchar()

#include <TH1F.h>
#include <TApplication.h>
#include <TCanvas.h>

#include "GEMRawDecoder.h"
#include "GEMInputHandler.h"
#include "GEMPhysHandler.h"
#include "PRDMapping.h"
#include "GEMPedestal.h"
#include "GEMRawPedestal.h"
#include "GEMHitDecoder.h"
#include "GEMConfigure.h"

int main(int argc, char* argv[])
{
  //TApplication theApp("app", &argc, argv);
  TApplication theApp("app", NULL, NULL);
  //TCanvas *c = new TCanvas("c", "APV Raw Signal", 10, 10, 1200, 1000);
  //c->Divide(3,3);
 
  //Analysis Configure
  GEMConfigure configure;
  configure.LoadConfigure();
 
  // Load data file
  string file_name;
  if((configure.GetRunType() != "PHYSICS")  )
  {
    if(argc < 2)
    {
      cout<<"Non Physis analysis: Usage:"<<endl;
      cout<<"                          ./main <file_name>"<<endl;
      return 0;
    }
    file_name=argv[1];
    cout<<"Analyzing file:  "<<file_name<<endl;
  }

  //TApplication theApp("app", &argc, argv);
  //  Load mapping fiel
  //string mapping_file("./mapping.cfg");
  string mapping_file(configure.GetMapping());
  cout<<"Loading Mapping File from:  "<<mapping_file.c_str()<<endl;
  PRDMapping* mapping = PRDMapping::GetInstance();
  mapping->LoadMapping(mapping_file.c_str());
  
  //---------------------------------------------------------
  // RAW Data Run
  if( configure.GetRunType().compare("RAWDATA") == 0)
  {
    cout<<"RAWDATA RUN ... "<<endl;
    //Process all Events
    GEMInputHandler inputHandler(file_name);
    int entries = inputHandler.ProcessAllEvents();
    cout<<"Total Events: "<<entries<<endl;
  }
  //---------------------------------------------------------
  // PEDESTAL RUN
  else if( configure.GetRunType().compare("PEDESTAL") == 0)
  {
    cout<<"PEDESTAL RUN ..."<<endl;
    GEMPedestal pedestal(file_name);
    pedestal.SavePedestalFile();
    //pedestal.LoadPedestal();
    vector<Float_t> noise = pedestal.GetAPVNoises(0);
    vector<Float_t> offset = pedestal.GetAPVOffsets(0);

    cout<< noise.size()<<endl;
    cout<<offset.size()<<endl;
    for(int i=0;i<128;i++)
    {
      cout<<"offset: "<<offset[i]<<"  noise:  "<<noise[i]<<endl;
    }
  }
  //---------------------------------------------------------
  // PHYSICS RUN
  else if( configure.GetRunType().compare("PHYSICS" ) == 0)
  {
    cout<<"PHYSICS RUN ..."<<endl;
    GEMPhysHandler phys_handler;
    //phys_handler.ProcessAllEvents(-1);
    phys_handler.ProcessAllFiles();
  }
  //---------------------------------------------------------
  // HIT RUN
  else if( configure.GetRunType().compare("HIT" ) ==0 )
  {
    cout<<"HITRUN"<<endl;
    GEMHitDecoder hitDecoder(file_name);
    //hitDecoder.ShowHits();
    hitDecoder.ProcessAllEvents();
    //GEMPedestal pedestal("./pedestal.root");
    //pedestal.LoadPedestal();
  }
  //---------------------------------------------------------
  // DEFAULT RUN
  else
  {
    cout<<"please specify a RUNTYPE ..."<<endl;
    cout<<mapping->GetFECIPFromFECID(1)<<endl;
    cout<<mapping->GetFECIPFromFECID(3)<<endl;
    cout<<mapping->GetFECIDFromFECIP(TString("10.0.1.2"))<<endl;
    cout<<mapping->GetFECIDFromFECIP(TString("10.0.8.2"))<<endl;
  }
  //theApp.Run(true);
  return 0;
}

