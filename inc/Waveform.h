/*---------------------------------------------------------------------*
| Waveform.h                                                           |
|                                                                      |
| Description: A class that holds an array of times, samples,          |
|              amplitudes, and baseline and FFT information for each   |
|              waveform. Contains corresponding functions for waveform |
|              processing.                                             |
|                                                                      |
*----------------------------------------------------------------------*/


#ifndef Waveform_h
#define Waveform_h

#include <iostream>
#include <map>
#include "Config.h"
#include "TH1D.h"

using namespace std;


class Waveform{
 public:
  Waveform(){
    _max = 0;
    _max_samp = 0;
    _min = 0;
    _min_samp = 0;
    for(int i_s = 0; i_s < Config::Get()->GetParameterI("num_samps"); i_s++){
        _samples.push_back(i_s);
        _times.push_back(i_s * 1e-3 / Config::Get()->GetParameterD("sampling_rate"));
	_amplitudes.push_back(0);
	_baseline.push_back(0);
	_baseline_stdev.push_back(0);
	_baseline_range.push_back(0);
    } 
  }
  ~Waveform() {;}

  // Get Functions
  float GetMin()                   {return _min;}
  int GetMinSamp()                  {return _min_samp;}
  float GetMax()                   {return _max;}
  int GetMaxSamp()                  {return _max_samp;}
  vector<int> GetSamples()          {return _samples;}
  vector<float> GetTimes()         {return _times;}
  vector<float> GetAmp()           {return _amplitudes;}
  float GetAmp(int s)              {return _amplitudes[s];}
  vector<float> GetBaseline()      {return _baseline;}
  float GetBaseline(int s)         {return _baseline[s];}
  vector<float> GetBaselineStdev()   {return _baseline_stdev;}
  float GetBaselineStdev(int s)    {return _baseline_stdev[s];}
  vector<float> GetBaselineRange() {return _baseline_range;}
  float GetBaselineRange(int s)    {return _baseline_range[s];}
  vector<float>  GetFFT()          {return _FFTMag;}
  float GetFFT(int s)              {return _FFTMag[s];}
  
  //Set Functions
  void CalcExtrema();
  void CalcBaseline(double threshold = Config::Get()->GetParameterD("pulse_threshold"));
  void SetAmp(int s,float x) {_amplitudes[s] = x;}
  void CalcFFT();


  // Waveform Analysis
  double Integrate(int s0 = 0, int s1 = Config::Get()->GetParameterI("num_samps"));
  double Integrate(double t0 = 0, double t1 = Config::Get()->GetParameterI("num_samps") * 1e-3 / Config::Get()->GetParameterD("sampling_rate"));
  Waveform* CalcRollingIntegral();  
  Waveform* SubtractBaseline();  
  Waveform* LowPassFilter();   
  Waveform* HighPassFilter(); 
  Waveform* Deconvolve(Waveform* response); 
  Waveform* WienerDeconvolve(Waveform* response);
  Waveform* FFTFilter(); 

  void WaveformFree();

 private:
  float          _max;
  int            _max_samp;
  float          _min;
  int            _min_samp;
  vector<int>    _samples;
  vector<float> _times;
  vector<float> _amplitudes;
  vector<float> _baseline;
  vector<float> _baseline_stdev;
  vector<float> _baseline_range;
  vector<float> _FFTMag;

};

#endif
