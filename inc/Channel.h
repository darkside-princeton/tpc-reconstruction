/*---------------------------------------------------------------------*
| Channel.h                                                            |
|                                                                      |
| Description: A class that holds waveforms and low-level processing   |
|              functions for each channel, including various methods   |
|	      of pulse-finding. Also holds vectors of pulse objects    |
|	      found through each pulse-finding method, and display     |
|             functions for the event viewrer.                         |
|                                                                      |
*----------------------------------------------------------------------*/


#ifndef Channel_h
#define Channel_h

#include <iostream>
#include <map>
#include "Config.h"
#include "Pulse.h"
#include "Waveform.h"

using namespace std;

class Channel{
 public:
  Channel() {
    Channel(-1);
  }
  Channel(int ch_id){
    _ch_id = ch_id;
    _gain = 0;
    _normalized = false;
    _npulses = 0;
    _total_integral = 0;
    _raw_wf = NULL;
    _integral_wf = NULL;
    _bl_sub_wf = NULL;
    _filtered_wf_fft = NULL;
    _filtered_wf_low = NULL;
    _filtered_wf_high = NULL;
    _filtered_wf_high2 = NULL;
    _filtered_wf_high3 = NULL;
    _diff_wf = NULL; 
    _deconv_wf = NULL;
    if(Config::Get()->GetParameterIvec("ch_pos")[ch_id] == 1) _top = true;
    else _top = false;
  }

  ~Channel() { ; }

  // Get Functions
  int            GetChannelID()     {return _ch_id;}
  double         GetGain()          {return _gain;}
  int            GetNPulses()       {return _pulses.size();}
  Waveform*      GetRawWF()         {return _raw_wf; }
  Waveform*      GetIntegralWF()    {return _integral_wf;}
  Waveform*      GetBlSubWF()       {return _bl_sub_wf;}
  Waveform*      GetFilteredWF()    {return _filtered_wf_low;}
  Waveform*      GetDeconvWF()      {return _deconv_wf;}
  Waveform*      GetFFTFilteredWF() {return _filtered_wf_fft;}
  Waveform*      GetDiffWF()        {return _diff_wf;}
  vector<Pulse*> GetPulses()        {return _pulses;}
  Pulse*         GetPulse(int p)    {return _pulses[p];}
  vector<Pulse*> GetDeconvPulses()  {return _deconv_pulses;}
  bool           GetTop()           {return _top;}
  bool           GetNormalized()    {return _normalized;}

  // Set Functions
  void SetChannelID(int ch_id) {_ch_id = ch_id; }
  void SetRawWF(Waveform *wf);
  void SetBlSubWF() {_bl_sub_wf = _raw_wf->SubtractBaseline();}
  void SetDeconvWF(Waveform* response);

  // Processing Functions
  void ProcessChannel();
  void ProcessNoise(); // function for noise-characterization
  void Normalize(); // divides all waveforms by _gain except _raw_wf and recalculates bl_sub_wf baseline
  void PulseDeconv(Waveform * response);

  // Display Functions
  void ChannelPrint();
  void ChannelDisplay(double min_x = 0, double max_x = Config::Get()->GetParameterI("num_samps") * 1e-3 / Config::Get()->GetParameterD("sampling_rate"), double min_y = 0, double max_y = 0);


  void Reset() {
    _gain = 0;
    _normalized = false;
    _npulses = 0;
    _total_integral = 0;
    _raw_wf = NULL;
    _integral_wf = NULL;
    _bl_sub_wf = NULL;
    _filtered_wf_fft = NULL;
    _filtered_wf_low = NULL;
    _filtered_wf_high = NULL;
    _filtered_wf_high2 = NULL;
    _filtered_wf_high3 = NULL;
    _diff_wf = NULL;
    _deconv_wf = NULL;
  }

  void ChannelFree();

 private:
  int            _ch_id;
  double         _gain;
  bool           _top;
  bool           _normalized;
  int            _npulses;
  double         _total_integral;
  Waveform*      _raw_wf;
  Waveform*      _integral_wf; //cumulative integral
  Waveform*      _bl_sub_wf; // baseline_subtracted waveform
  Waveform*      _filtered_wf_low; // LP filtered waveform
  Waveform*      _filtered_wf_high; // 1-pole HP filtered waveform 
  Waveform*      _filtered_wf_high2; // 2-pole HP filtered waveform
  Waveform*      _filtered_wf_high3; // 3-pole HP filtered waveform
  Waveform*      _filtered_wf_fft; //filtered waveform using fft
  Waveform*      _diff_wf; // derivative of LP filtered waveform
  Waveform*      _deconv_wf; // deconvolved waveform
  Waveform*      _deconv_diff; //derivative of deconvolved waveform
  vector<Pulse*> _pulses; // pulses identified through threshold on _diff_wf
  vector<Pulse*> _deconv_pulses; // pulses identified through threshold on _deconv_diff

  // Set Functions
  void SetWaveforms();
  void SetGain();

  // Processing Functions
  void ClearPulses() {_pulses.clear();}
  void FindPulses();
  void FindAfterPulses(Pulse* p, int pulse_start, int pulse_end);

  // Display Functions
  void PlotPulseTimes();
  void PlotLPFilter(double min_x = 0, double max_x = Config::Get()->GetParameterI("num_samps") * 1e-3 / Config::Get()->GetParameterD("sampling_rate"), double min_y = 0, double max_y = 0);
  void PlotHPFilter(double min_x = 0, double max_x = Config::Get()->GetParameterI("num_samps") * 1e-3 / Config::Get()->GetParameterD("sampling_rate"), double min_y = 0, double max_y = 0);
  void PlotBaseline(double min_x = 0, double max_x = Config::Get()->GetParameterI("num_samps") * 1e-3 / Config::Get()->GetParameterD("sampling_rate"), double min_y = 0, double max_y = 0);
  void PlotDerivative(double min_x = 0, double max_x = Config::Get()->GetParameterI("num_samps") * 1e-3 / Config::Get()->GetParameterD("sampling_rate"), double min_y = 0, double max_y = 0);
  void PlotIntegral(double min_x = 0, double max_x = Config::Get()->GetParameterI("num_samps") * 1e-3 / Config::Get()->GetParameterD("sampling_rate"), double min_y = 0, double max_y = 0);
  void PlotDeconv(double min_x = 0, double max_x = Config::Get()->GetParameterI("num_samps") * 1e-3 / Config::Get()->GetParameterD("sampling_rate"), double min_y = 0, double max_y = 0);
  void PlotFFTFilter(double min_x = 0, double max_x = Config::Get()->GetParameterI("num_samps") * 1e-3 / Config::Get()->GetParameterD("sampling_rate"), double min_y = 0, double max_y = 0);

};

#endif
