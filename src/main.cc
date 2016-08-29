#include <iostream>
#include <fstream>
#include <map>
#include <string>

#include <TH1F.h>
#include <TApplication.h>
#include <TCanvas.h>

#include "GEMRawDecoder.h"
#include "GEMAnalyzer.h"
#include "GEMMapping.h"
#include "GEMConfigure.h"

int main(int argc, char* argv[])
{
  TApplication theApp("app", NULL, NULL);
 
  //Analysis Configure
  GEMConfigure configure;
  configure.LoadConfigure();
  string mapping_file(configure.GetMapping());
  cout<<"Loading Mapping File from:  "<<mapping_file.c_str()<<endl;
  GEMMapping* mapping = GEMMapping::GetInstance();
  mapping->LoadMapping(mapping_file.c_str());

  GEMAnalyzer *gem_analyzer = new GEMAnalyzer();
  gem_analyzer->ProcessFiles();
  delete gem_analyzer;
  //theApp.Run(true);
  return 0;
}

