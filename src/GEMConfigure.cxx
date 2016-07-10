#include "GEMConfigure.h"

using namespace std;

GEMConfigure::GEMConfigure()
{
  nFile = 0;
  TDC_Channel ="";
  TDC_Quan = 0;
  TDC_Start = 0;
  TDC_End = 0;
  Hycal_Timing_Cut_Start = 0;
  Hycal_Timing_Cut_End = 0;
  Hycal_Energy = 0;
  configure_file = "./gem.cfg";
  //LoadConfigure();
}

GEMConfigure::GEMConfigure(const char * file)
{
  configure_file = file;
  //LoadConfigure();
}

GEMConfigure::~GEMConfigure()
{
}

void GEMConfigure::LoadConfigure()
{
  ifstream filestream( configure_file.c_str(), ifstream::in);
  string line;
  int kk_file = 0;
  while( getline(filestream, line) )
  {
    line.erase( remove_if(line.begin(), line.end(), ::isspace), line.end() ); // strip spaces
    if(line.find("#") == 0) continue; //skip comments

    char* tokens = strtok( (char*)line.data(), " :,");
    while( tokens != NULL)
    {
      string s = tokens;
      if(s == "RUNTYPE")
      {
        tokens = strtok(NULL, " :,");
	runType = tokens;
      }
      else if(s == "MAPPING")
      {
        tokens = strtok(NULL, " :,");
	mapping = tokens;
      }
      else if(s=="SAVEPED")
      {
        tokens = strtok(NULL, " :,");
	save_pedestal = tokens;
      }
      else if(s=="NEVENTPED")
      {
        tokens = strtok(NULL, " :,");
	nPedestal = atoi(tokens);
      }
      else if(s=="LOADPED")
      {
        tokens = strtok(NULL, " :,");
        load_pedestal = tokens;
      }
      else if(s=="NFILE")
      {
        tokens = strtok(NULL, " :,");
	nFile = atoi(tokens);
      }
     else if(s=="INPUTFILE")
      {
	tokens = strtok(NULL, " :,");
	fileHeader = tokens;
        tokens = strtok(NULL, " :,");
	evioStart = atoi(tokens);
	tokens = strtok(NULL, " :,");
	evioEnd = atoi(tokens);

	for(int i=evioStart; i<= evioEnd;i++)
	{
	  fileList[kk_file++] = fileHeader+"."+to_string(i);
	}
      }
      else if(s=="HyCalTDCGroupQuantity")
      {
        tokens = strtok(NULL, " :,");
	TDC_Quan = atoi(tokens);
        cout<<"TDC_Quan: "<<TDC_Quan<<endl;
      }
      else if(s=="HyCalTDCGroup")
      {
        int ii = TDC_Quan;
	while(--ii>=0)
	{
	 tokens = strtok(NULL, " :,");
	 TDC[ii] = tokens;
         cout<<TDC[ii]<<endl;
	}
      }
      else if(s == "TDCChannel")
      {
        tokens = strtok(NULL, " :,");
	TDC_Channel = tokens;
      }
      else if(s == "TDCValueStart")
      {
        tokens = strtok(NULL, " :,");
        TDC_Start = atof(tokens);
      }
      else if(s == "TDCValueEnd")
      {
        tokens = strtok(NULL, " :,");
	TDC_End = atof(tokens);
      }
      else if( s=="HycalEnergy")
      {
        tokens = strtok(NULL, " :,");
	Hycal_Energy = atof(tokens);
      }
      else if(s == "PHYSRESULTS")
      {
        tokens = strtok(NULL, " :,");
	phys_results_path = tokens;
      }
      else if(s == "HyCalTimingCutStart")
      {
        tokens = strtok(NULL, " :,");
        Hycal_Timing_Cut_Start = atof(tokens);
      }
      else if(s == "HyCalTimingCutEnd")
      {
        tokens = strtok(NULL, " :,");
	Hycal_Timing_Cut_End = atof(tokens);
      }
      else if(s == "UseScinTDC")
      {
        tokens = strtok(NULL, " :,");
	UseScinTDCCut = atoi(tokens);
      }
      else if(s == "UseHyCalPosMatch")
      {
        tokens = strtok(NULL, " :,");
	UseHyCalPosMatch = atoi(tokens);
      }
      else if(s == "UseHycalTimingCut")
      {
        tokens = strtok(NULL, " :,");
	UseHyCalTimingCut = atoi(tokens);
      }
      else
      {
        tokens = strtok(NULL, " :,");
      }
    }
  }
  
  assert(kk_file == nFile);
}

string GEMConfigure::GetRunType()
{
  return runType;
}

string GEMConfigure::GetMapping()
{
  return mapping;
}

string GEMConfigure::GetSavePedPath()
{
  return save_pedestal;
}

string GEMConfigure::GetLoadPedPath()
{
  return load_pedestal;
}


int GEMConfigure::GetNumEvtForPed()
{
  return nPedestal;
}
