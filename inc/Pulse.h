/*---------------------------------------------------------------------* 
| Pulse.h                                                              |
|                                                                      |
| Description: A Pulse object holds time and charge information        |
|              about the arrivals of individual PE pulses. Each        |
|              Pulse also contains a list of afterpulses identified    |
|              within its recharge time.                               |
|                                                                      | 
*----------------------------------------------------------------------*/


#ifndef Pulse_h
#define Pulse_h

#include <iostream>
#include <map>
#include "Config.h"

using namespace std;

class Pulse{
 public:
  Pulse(){
    _start_samp = -1;
    _start_time = -1;
    _end_samp = -1;
    _end_time = -1;
    _integral = -1;
    _height = -1;
    _filtered_height = -1;
    _max_val = -1;
    _min_val = -1;
    _peak_time = -1;
    _peak_samp = -1;
    _deconv_integral = -1;
  }
  ~Pulse() {;}


  //Get Functions
  int    GetStartSamp()      {return _start_samp;}
  double GetStartTime()      {return _start_time;}
  int    GetEndSamp()        {return _end_samp;}
  double GetEndTime()        {return _end_time;}
  double GetIntegral()       {return _integral;}
  double GetMaxVal()         {return _max_val;}
  double GetMinVal()         {return _min_val;}
  double GetHeight()         {return _height;}
  double GetFilteredHeight() {return _filtered_height;}
  double GetPeakTime()       {return _peak_time;}
  int    GetPeakSamp()       {return _peak_samp;}
  double GetHalfTime()       {return _half_time;}
  double GetDeconvIntegral() {return _deconv_integral;}
  int    GetNAfterPulses()   {return _after_pulses.size();}
  
  vector<Pulse*> GetAfterPulses() {return _after_pulses;}
  Pulse* GetAfterPulse(int p) {return _after_pulses[p];}


  //Set Functions
  void SetStartSamp(int start_samp)              {_start_samp = start_samp;}
  void SetStartTime(double start_time)           {_start_time = start_time;}
  void SetEndSamp(int end_samp)                  {_end_samp = end_samp;}
  void SetEndTime(double end_time)               {_end_time = end_time;}
  void SetIntegral(double integral)              {_integral = integral;}
  void SetMaxVal(double max_val)                 {_max_val = max_val;}
  void SetMinVal(double min_val)                 {_min_val = min_val;}
  void SetHeight(double height)                  {_height = height;}
  void SetFilteredHeight(double height)          {_filtered_height = height;}
  void SetPeakTime(double peak_time)             {_peak_time = peak_time;}
  void SetPeakSamp(int peak_samp)                {_peak_samp = peak_samp;}
  void SetHalfTime(double half_time)             {_half_time = half_time;}
  void SetDeconvIntegral(double deconv_integral) {_deconv_integral = deconv_integral;}
  void AddAfterPulse(Pulse* p)                   {_after_pulses.push_back(p);}
  void ClearAfterPulses()                        {_after_pulses.clear();}


  void Print(){
    cout << "start_time = " << _start_time << " , end_time = " << _end_samp * 1e-3 / Config::Get()->GetParameterD("sampling_rate") << 
" , max_val = " << _max_val << " , min_val = " << _min_val << " , integral = " << _integral << " , height = " << _height <<  " , filtered height = " <<  _filtered_height << " , half_time = " << _half_time << endl;
    for(unsigned int i = 0; i < _after_pulses.size(); i++) {
      cout << "     AfterPulse " << i << ": ";
      cout << "start_time = " << _after_pulses[i]->GetStartTime() << " , end_time = " << _after_pulses[i]->GetEndTime() << " , peak_time = " << _after_pulses[i]->GetPeakTime() 
	   << " , height = " << _after_pulses[i]->GetHeight() << " , integral = " << _after_pulses[i]->GetIntegral() <<  endl;
    }
  }

  void PulseFree() {
    for(unsigned int i = 0; i < _after_pulses.size(); i++) delete _after_pulses[i];
    _after_pulses.clear(); vector<Pulse*>().swap(_after_pulses);
  }


 private:
  int _start_samp;
  int _end_samp;
  double _start_time;
  double _end_time;
  double _integral;
  double _max_val;
  double _min_val;
  double _height;
  double _filtered_height;
  int _nafter_pulses;
  double _peak_time;
  int _peak_samp;
  double _half_time; // time at half-height
  double _deconv_integral;
  vector<Pulse*> _after_pulses; //vector of pulses that occur within rise time of pulse, not including primary pulse
};

#endif
