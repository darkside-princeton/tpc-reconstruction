/*---------------------------------------------------------------------*
| Output.h                                                             |
|                                                                      |
| Description: The Output class handles output TTrees and TGraphs for  |
|              each mode of running the software.                      |
|                                                                      |
*----------------------------------------------------------------------*/


#ifndef output_h
#define output_h

#include <algorithm>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <TFile.h>
#include <TObject.h>
#include <TTree.h>
#include "Event.h"
#include "Config.h"
#include "TH1D.h"


using namespace std;


class Output{
  static Output *_oh;

  Output(){
    _ofile = 0;
    total_S1_pulses = 0;
    for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){
      total_pulses.push_back(0);
      total_SPE_pulses.push_back(0);
    }
  }

  ~Output(){ ; }

 public:

  static Output *Get() {
    if(!_oh)
      _oh = new Output;
    return _oh;
  }

  //General Variables 
  TFile *_ofile;
  int eventID;


  // Noise Characterization
  void     SetupNoise();
  void     WriteNoise();
  
  TTree**  _noiseTree;
  double   baseline_mean;
  double   baseline_stdev;
  TH1D**   hAvgFFT; //for each channel
  TGraph** gHitsThreshold; //for each channel
  TGraph** gHitsToT; //for each channel        


  //SiPM characterization
  void SetupSiPM();
  void WriteSiPM();

  vector<int> total_pulses; // For average pulse
  vector<int> total_SPE_pulses; // For average SPE pulse
  
  TTree*                           _pulseTree;
  vector<int>                      npulses;
  vector<vector<int> >             nafterpulses;
  vector<vector<double> >          start_time;
  vector<vector<double> >          end_time;
  vector<vector<double> >          integral;
  vector<vector<double> >          max_val;
  vector<vector<double> >          min_val;
  vector<vector<double> >          height;
  vector<vector<double> >          filtered_height;
  vector<vector<double> >          peak_time;
  vector<vector<double> >          half_time;
  vector<vector<double> >          filtered_integral;
  vector<vector<vector<double> > > ap_peak_times;
  vector<vector<vector<double> > > ap_heights;
  vector<vector<double> >          deconv_integral;
  TGraph**                         gAvgPulse;
  TGraph**                         gAvgSPEPulse;
  TGraph**                         gAvgWF;


  //Generating waveforms
  void SetupGenerateWF();
  void WriteGenerateWF();

  TTree**                  _trueTree;
  vector<int>              nPE;
  vector<double>           rTime;
  vector<vector<double> >  rEventID;
  vector<vector<double> >  rPulseID;
  vector<vector<double> >  rStart_time;
  vector<vector<double> >  rIntegral;
  vector<vector<double> >  rHeight;
  vector<vector<double> >  rAfterPulses;


  // TPC runs
  void SetupTPC();
  void WriteTPC();


  int total_S1_pulses; // for average S1

  TTree*                           _clusterTree;
  int                              nclusters;
  int                              overlap;
  int                              S1_pulses;
  double                           S1_charge;
  int                              S2_pulses;
  double                           S2_charge;
  double                           S1_S2_dt;
  vector<double>                   cl_start_time;
  vector<double>                   cl_end_time;
  vector<double>                   cl_integral;
  vector<double>                   cl_max_val;
  vector<double>                   cl_min_val;
  vector<double>                   cl_height;
  vector<double>                   cl_peak_time;
  vector<double>                   cl_half_time;
  vector<double>                   cl_symmetry;
  vector<vector<double> >          cl_fprompt;
  vector<vector<vector<double> > > pulse_peak_time;
  vector<vector<vector<double> > > pulse_height;
  vector<vector<double> >          cl_deconv_integral;
  vector<vector<double> >          cl_deconv_fprompt;
  TGraph*                          gAvgS1;


};



#endif
