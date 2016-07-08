#include <iostream>
#include <fstream>
#include <map>
#include <string>

#include <stdio.h>

#include <TH1F.h>
#include <TApplication.h>
#include <TCanvas.h>

#include "GEMPhysHandler.h"
#include "GEMConfigure.h"

int main(int argc, char* argv[])
{
  //TApplication theApp("app", &argc, argv);
  TApplication theApp("app", NULL, NULL);
 
  //Analysis Configure
  GEMConfigure configure;
 
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

  // PHYSICS RUN
  if( configure.GetRunType().compare("PHYSICS" ) == 0)
  {
    cout<<"PHYSICS RUN ..."<<endl;
    GEMPhysHandler phys_handler;
    //phys_handler.ProcessAllEvents(-1);
    phys_handler.ProcessAllFiles();
  }

  // DEFAULT RUN
  else
  {
    cout<<"please specify a RUNTYPE ..."<<endl;
  }
  //theApp.Run(true);
  return 0;
}

