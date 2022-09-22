#include <fstream>
#include <iostream>
#include <climits>
#include <string>
#include <vector>
#include "TGraph.h"
#include "TFile.h"
#include "Run.h"
#include "Config.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TRootCanvas.h"

using namespace std;


int main(int argc, char* argv[]){
  string cfg_fname = "";
  string rundir = "";
  bool display = false;
  bool noise = false;
  bool TPC = false;
  bool generate = false;
  bool SiPM = false;
  int nevts = -1;
  for(int i_arg = 1; i_arg < argc; i_arg++){
    string arg = argv[i_arg];
    if(arg == "-c" && i_arg+1 < argc)        cfg_fname = argv[i_arg+1];
    if(arg == "-d" && i_arg+1 < argc)        rundir    = argv[i_arg+1];
    if(arg == "--display")                   display  = true;
    if(arg == "--noise")                     noise = true;
    if(arg == "--SiPM")                      SiPM = true;
    if(arg == "--TPC")                       TPC = true;
    if(arg == "--generate")                  generate = true;
    if(arg == "--nevents" && i_arg+1 < argc) nevts     = atoi(argv[i_arg+1]);
  } // for i_arg                                                               

  
  if(rundir == "") {
    cout << "ERROR: Enter a run directory" << endl;
    return 1;
  }                                                                                                                                                                       
  if(cfg_fname == "") {
    cout << "ERROR: Enter a config file name" << endl;
    return 1;
  } 


  Config::Get()->SetParameter("nevents",nevts);
  if(TPC) Config::Get()->SetParameter("TPC","y");
  else Config::Get()->SetParameter("TPC","n");
  Config::Get()->AddFile(cfg_fname);
  

  Run* run = new Run(rundir);
  if (noise) run->ProcessNoise();
  else if (TPC && !display) run->ProcessTPC();
  //else if (display_deconv) {TApplication *myApp = new TApplication("myApp", &argc, argv); run->EventViewerDeconv(false); myApp->Run();}
  else if(generate) run->GenerateWaveforms();  
  else if (display)  {
    TApplication *myApp = new TApplication("myApp", &argc, argv);
    cout << "EVENT VIEWER (enter q to exit)" << endl;
    run->EventViewer();
    myApp->Run();
  }
  else if (SiPM && !display) run->ProcessSiPM();     
  else {
    cout << "Please specify a software mode. " << endl;
  }

  return 0;
}
