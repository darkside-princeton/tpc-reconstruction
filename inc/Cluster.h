/*---------------------------------------------------------------------*
| Cluster.h                                                            |
|                                                                      |
| Description: A Cluster object holds time and charge information      |
|              about correlated arrivals of PEs across channels. Each  |
|              Cluster is also classified as an S1/S2 signal.          |
|                                                                      |
*----------------------------------------------------------------------*/

#ifndef Cluster_h
#define Cluster_h

#include <iostream>
#include <map>
#include "Config.h"
#include "Pulse.h"

using namespace std;

class Cluster{
 public:
  Cluster(){
    _start_samp = -1;
    _start_time = -1;
    _end_samp = -1;
    _end_time = -1;
    _integral = -1;
    _height = -1;
    _max_val = -1;
    _min_val = -1;
    _peak_time = -1;
    _peak_samp = -1;
    _symmetry = -1;
    for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){
      vector<Pulse*> pulse_vec;
      vector<double> deconv_vec;
      _pulses.push_back(pulse_vec);
      _deconv_integral.push_back(0);
    }
    for(unsigned int i_fp = 0; i_fp < Config::Get()->GetParameterDvec("prompt_windows").size(); i_fp++){
      _deconv_fprompt.push_back(-1);
      _fprompt.push_back(-1);
    }
    
  }
  ~Cluster() {;}

  //Get Functions
  int            GetStartSamp()              {return _start_samp; }
  double         GetStartTime()              {return _start_time; }
  int            GetEndSamp()                {return _end_samp;}
  double         GetEndTime()                {return _end_time;}
  double         GetIntegral()               {return _integral;}
  double         GetMaxVal()                 {return _max_val;}
  double         GetMinVal()                 {return _min_val;}
  double         GetHeight()                 {return _height;}
  double         GetPeakTime()               {return _peak_time;}
  int            GetPeakSamp()               {return _peak_samp;}
  vector<double> GetFprompt()                {return _fprompt;}
  double         GetFprompt(int i)           {return _fprompt[i];}
  double         GetHalfTime()               {return _half_time;}
  vector<Pulse*> GetPulses(int i_ch)         {return _pulses[i_ch];}
  Pulse*         GetPulse(int i_ch, int i_p) {return _pulses[i_ch][i_p];}
  double         GetDeconvIntegral(int i_ch) {return _deconv_integral[i_ch];}
  vector<double> GetDeconvFprompts(int i_ch) {return _deconv_fprompt;}
  double         GetDeconvFprompt(int i_fp)  {return _deconv_fprompt[i_fp];}
  bool           GetS1()                     {return _S1;}
  double         GetSymmetry()               {return _symmetry;}

  //Set Functions
  void SetStartSamp(int start_samp)                        {_start_samp = start_samp;}
  void SetStartTime(double start_time)                     {_start_time = start_time;}
  void SetEndSamp(int end_samp)                            {_end_samp = end_samp;}
  void SetEndTime(double end_time)                         {_end_time = end_time;}
  void SetIntegral(double integral)                        {_integral = integral;}
  void SetMaxVal(double max_val)                           {_max_val = max_val;}
  void SetMinVal(double min_val)                           {_min_val = min_val;}
  void SetHeight(double height)                            {_height = height;}
  void SetPeakTime(double peak_time)                       {_peak_time = peak_time;}
  void SetPeakSamp(int peak_samp)                          {_peak_samp = peak_samp;}
  void SetFprompt(vector<double> fprompt)                  {_fprompt = fprompt;}
  void SetHalfTime(double half_time)                       {_half_time = half_time;}
  void AddPulse(int i_ch, Pulse* pulse)                    {_pulses[i_ch].push_back(pulse);}
  void SetDeconvIntegral(int i_ch, double deconv_integral) {_deconv_integral[i_ch] = deconv_integral;}
  void SetDeconvFprompts(vector<double> fprompts)          {_deconv_fprompt = fprompts;} 
  void SetS1(bool S1)                                      {_S1 = S1;}
  void SetSymmetry(double symmetry)                        {_symmetry = symmetry;}


  void ClearPulses() {for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){ _pulses[i_ch].clear();}}


  void Print(){
    cout << "start_time = " << _start_time << " , end_time = " << _end_samp * 1e-3 / Config::Get()->GetParameterD("sampling_rate") <<
      " , max_val = " << _max_val << " , min_val = " << _min_val << " , integral = " << _integral << " , height = " << _height << 
      " , half_time = " << _half_time << " , symmetry = " << _symmetry << endl;
    for(unsigned int i_fp = 0; i_fp < Config::Get()->GetParameterDvec("prompt_windows").size(); i_fp++){
      cout << Form("fprompt %d = ", i_fp+1) << _fprompt[i_fp];
      if(i_fp !=  Config::Get()->GetParameterDvec("prompt_windows").size()) cout << " , ";
    }
    cout << endl;
    if(Config::Get()->GetParameterS("plot_deconv") == "y" && Config::Get()->GetParameterS("deconvolve") == "y") {
      cout << "Deconvolved Integrals: " << endl;
      for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++) {
	cout << "Channel " << i_ch << ": integral = " << _deconv_integral[i_ch] << endl;
      }
      for(unsigned int i_fp = 0; i_fp < Config::Get()->GetParameterDvec("prompt_windows").size(); i_fp++){
	cout << Form("fprompt %d = ", i_fp+1) << _deconv_fprompt[i_fp];
	if(i_fp !=  Config::Get()->GetParameterDvec("prompt_windows").size()) cout << " , ";
      }
      cout << endl;
    }
  }

  void ClusterFree() {
    for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++) {
      for(unsigned int i_p = 0; i_p < _pulses[i_ch].size(); i_p++) {
	_pulses[i_ch][i_p]->PulseFree();
	delete _pulses[i_ch][i_p];
      }
      _pulses[i_ch].clear(); vector<Pulse*>().swap(_pulses[i_ch]);
    }
    _deconv_integral.clear(); vector<double>().swap(_deconv_integral);
    _deconv_fprompt.clear(); vector<double>().swap(_deconv_fprompt);
    _pulses.clear(); vector<vector<Pulse*> >().swap(_pulses);
   
  }

 private:
  int                     _start_samp;
  int                     _end_samp;
  double                  _start_time;
  double                  _end_time;
  double                  _integral;
  double                  _max_val;
  double                  _min_val;
  double                  _height;
  int                     _nafter_pulses;
  double                  _peak_time;
  int                     _peak_samp;
  vector<double>          _fprompt;
  double                  _half_time;
  double                  _symmetry;
  vector<vector<Pulse*> > _pulses;
  vector<double>          _deconv_integral;
  vector<double>          _deconv_fprompt;
  bool _S1; //true if the pulse is an S1 pulse, false otherwise

};

#endif

