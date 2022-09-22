/*---------------------------------------------------------------------* 
| Run.h                                                                |
|                                                                      |
| Description: The Run class organizes different modes of running the  |
|              software and calculates averages across channels.       |
|                                                                      |
*----------------------------------------------------------------------*/

#ifndef Run_h
#define Run_h

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include "Event.h"
#include "Config.h"
#include "TGraph.h"
#include "TFile.h"
#include "Waveform.h"
#include "Channel.h"
#include "TRandom.h"

using namespace std;

class Run{
  public:
  
  Run(string rundir) {
    int NChannel = Config::Get()->GetParameterI("num_chans");
    ifstream** fin = new ifstream*[NChannel];
    for(int i_ch = 0; i_ch < NChannel; i_ch++){
	fin[i_ch] = new ifstream();
	vector<float> ch_wf;
	vector<float> ch_pulse;
	string data = rundir + "/wave%d.dat";
	fin[i_ch]->open(Form(data.c_str(),Config::Get()->GetParameterIvec("ch_ids")[i_ch]), std::fstream::binary);
	if(!fin[i_ch]->is_open()){
	  cout<<"file not found"<<endl;
	}
	for(int i_s = 0; i_s < Config::Get()->GetParameterI("num_samps"); i_s++){
	  ch_wf.push_back(0);
	  if(i_s < Config::Get()->GetParameterI("pulse_avg_samps") + Config::Get()->GetParameterI("pulse_avg_buf")) ch_pulse.push_back(0);
	  if(i_s < Config::Get()->GetParameterI("S1_avg_samps") + Config::Get()->GetParameterI("S1_avg_buf")) _avgS1.push_back(0);
	}
	_avgWF.push_back(ch_wf);
	_avgPulseAll.push_back(ch_pulse);
	_avgPulse.push_back(ch_pulse);
	_ifiles.push_back(fin[i_ch]);
    }    
  }
  ~Run() { ; }
  

  // Software modes
  int ProcessSiPM();
  void ProcessNoise();
  void ProcessTPC();
  void EventViewer();
  void GenerateWaveforms();

  // Averages
  void CalcAvgWF();
  void CalcAvgPulse();
  void CalcAvgS1();

 private:
  vector<ifstream*>      _ifiles; // one file for each channel                                                                                                                    
  vector<Event*>         _events;
  vector<vector<float> > _avgWF;
  vector<vector<float> > _avgPulseAll;
  vector<vector<float> > _avgPulse;
  vector<float>          _avgS1;
  
  Waveform* ReadBinary(int channel, bool header);
  Waveform* ReadBinary2(int channel, bool header);
  void AddRandPulses(TRandom* r1, double* wf, int npulses,  double tau, int offset, int sig_num);
  void AddTemplatePulses(TGraph* avg, TRandom* r1, double* wf, int npulses,  double tau, int offset, int sig_num);
  int Finish();

};


#endif
