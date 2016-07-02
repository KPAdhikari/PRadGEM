//======================================//
//                                      //
// Xinzhan Bai                          //
// 03/20/2016                           //
// xb4zp@virginia.edu                   //
//======================================//

#ifndef __GEMCONFIGURE_H__
#define __GEMCONFIGURE_H__

#include <string>
#include <cstring>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cassert>

using namespace std;

class GEMConfigure
{
public:
  GEMConfigure(const char * file );
  GEMConfigure();
  ~GEMConfigure();

  std::string GetRunType();
  string GetMapping();
  string GetSavePedPath();
  string GetLoadPedPath();
  int GetNumEvtForPed();

  void LoadConfigure();

  // make following variables public
  // input files to be analyzed...
  string fileList[1000];
  string fileHeader;
  int nFile;
  int evioStart;
  int evioEnd;
  // physics cut
  int UseScinTDCCut;
  int UseHyCalTimingCut;
  int UseHyCalPosMatch;

  string TDC_Channel;
  double TDC_Start;
  double TDC_End;
  double Hycal_Energy;
  double Hycal_Timing_Cut_Start;
  double Hycal_Timing_Cut_End;
  
  // file to save physics analysis results
  string phys_results_path;

private:
  string configure_file;

  string runType;
  int maxCluster;
  int minCluster;
  int nPedestal;
  
  string mapping;
  string save_pedestal;
  string load_pedestal;

};

#endif
