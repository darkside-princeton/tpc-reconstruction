#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <climits>
#include <vector>
#include "Event.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TMath.h"
#include "Pulse.h"
#include "Waveform.h"
#include "Channel.h"
#include "Output.h"

using namespace std;

/*--------------------------------------------------------------------*/

void Channel::SetRawWF(Waveform *wf) {
  _raw_wf = wf;
  int num_samps = Config::Get()->GetParameterI("num_samps");
  if(Config::Get()->GetParameterS("polarity") == "neg") {
    for(int i_s = 0 ; i_s < num_samps; i_s++){
      _raw_wf->SetAmp(i_s, _raw_wf->GetAmp(i_s) * -1);
    }
  }
}

/*--------------------------------------------------------------------*/

void Channel::SetGain() {
  ifstream gain_file("Config/gain.txt");
  if(!gain_file.is_open()){
    cout << "ERROR: Could not open file gain.txt" << endl;
    return;
  }
  string buff;
  while(getline(gain_file,buff)){
    if(buff == "") continue;
    stringstream ss(buff);
    string runID, chID, val;
    ss >> runID >> chID >> val;
    if(runID != Config::Get()->GetParameterS("gain_runID")) continue;
    if(stoi(chID) != _ch_id) continue;
    _gain = stod(val);
  }
}


/*--------------------------------------------------------------------*/

void Channel::ProcessChannel(){
  SetWaveforms();
  ClearPulses();
  FindPulses();
}


/*--------------------------------------------------------------------*/
void Channel::ProcessNoise() {
  _raw_wf->CalcBaseline();
  SetBlSubWF();
  _bl_sub_wf->CalcFFT();

  // Determining the total hits as a function of threshold and time over threshold
  int peak_start, peak_end;
  int num_samps = Config::Get()->GetParameterI("num_samps");
  double sampling_rate = Config::Get()->GetParameterD("sampling_rate");
  int nevents = Config::Get()->GetParameterI("nevents");

  for(int i_t = 0; i_t <  Output::Get()->gHitsThreshold[_ch_id]->GetN(); i_t++) {
    for(int i_s = 0 ; i_s < num_samps; i_s++){
      if(_bl_sub_wf->GetAmp(i_s) > Output::Get()->gHitsThreshold[_ch_id]->GetX()[i_t] * _raw_wf->GetBaselineStdev(i_s)) {
	Output::Get()->gHitsThreshold[_ch_id]->GetY()[i_t] +=  
	  1/(nevents * num_samps * 1e-3/ sampling_rate);
	peak_start = i_s;
	peak_end = i_s;
	while(_bl_sub_wf->GetAmp(peak_end) > Output::Get()->gHitsThreshold[_ch_id]->GetX()[i_t] * _raw_wf->GetBaselineStdev(peak_end) && peak_end < num_samps) 
	  peak_end++;
	i_s = peak_end;
	if(Output::Get()->gHitsThreshold[_ch_id]->GetX()[i_t] == 1) {
	  for(int i_ti = 0; i_ti < Output::Get()->gHitsToT[_ch_id]->GetN();i_ti++){
	    if(peak_end - peak_start > Output::Get()->gHitsToT[_ch_id]->GetX()[i_ti])
	      Output::Get()->gHitsToT[_ch_id]->GetY()[i_ti] += 
		1/(nevents * num_samps * 1e-3/ sampling_rate);
	  }
	}
      }
    }
  }
}

/*--------------------------------------------------------------------*/

void Channel::SetWaveforms() {
  _raw_wf->CalcBaseline();
  _bl_sub_wf = _raw_wf->SubtractBaseline();
  float max_val_raw = INT_MIN; int max_samp_raw = 0;
  float min_val_raw = INT_MAX; int min_samp_raw = 0;
  float max_val_sub = INT_MIN; int max_samp_sub = 0;
  float min_val_sub = INT_MAX; int min_samp_sub = 0;
  _total_integral = 0;
  _integral_wf = new Waveform();
  _integral_wf->SetAmp(0,0); int int_samp = 0;
  _filtered_wf_low = new Waveform; int filt_samp = 0;
  double cutoff_low = Config::Get()->GetParameterD("cutoff_freq_low");
  double RC_low = 1/ (2 * 3.1415 * cutoff);
  double alpha_low = (1e-3 /  Config::Get()->GetParameterD("sampling_rate")) / (RC_low + 1e-3 /  Config::Get()->GetParameterD("sampling_rate"));
  int num_samps = Config::Get()->GetParameterI("num_samps");
  int derivative_offset = Config::Get()->GetParameterI("derivative_offset");

  _filtered_wf_high = new Waveform;
  double cutoff_high = Config::Get()->GetParameterD("cutoff_freq_high");
  double RC_high = 1/ (2 * 3.1415 * cutoff);
  double alpha_high = RC_high / (RC_high + (1e-3 /  Config::Get()->GetParameterD("sampling_rate")));
  
  _filtered_wf_low->SetAmp(0, _bl_sub_wf->GetAmp(0));
  _filtered_wf_high->SetAmp(0, _bl_sub_wf->GetAmp(0));

  _diff_wf = new Waveform();
  for(int i_s = 0; i_s < num_samps; i_s++) {
    
    // Extrema for raw waveform
    if(_raw_wf->GetAmp(i_s) > max_val_raw){
      max_val_raw = _raw_wf->GetAmp(i_s);
      max_samp_raw = i_s;
    }
    if(_raw_wf->GetAmp(i_s) < min_val_raw){
      min_val_raw = _raw_wf->GetAmp(i_s);
      min_samp_raw = i_s;
    }

    // Extrema for subtracted waveform
    if(_bl_sub_wf->GetAmp(i_s) > max_val_sub){
      max_val_sub = _bl_sub_wf->GetAmp(i_s);
      max_samp_sub = i_s;
    }
    if(_bl_sub_wf->GetAmp(i_s) < min_val_sub){
      min_val_sub = _bl_sub_wf->GetAmp(i_s);
      min_samp_sub = i_s;
    }

    // Calculating total integral
    _total_integral += _bl_sub_wf->GetAmp(i_s);
  
    // Calculating rolling integral
    int_samp = min(i_s + 1, num_samps-1);
    _integral_wf->SetAmp(int_samp , _bl_sub_wf->GetAmp(int_samp ) + _integral_wf->GetAmp(int_samp -1));

    // Calculating low-pass filtered waveform
    filt_samp = min(i_s, num_samps-2);
    _filtered_wf_low->SetAmp(filt_samp + 1, (1-alpha_low) * _filtered_wf_low->GetAmp(filt_samp) + alpha_low * _bl_sub_wf->GetAmp(filt_samp ));  

    // Calculating high-pass filtered waveform
    _filtered_wf_high->SetAmp(filt_samp + 1, alpha_high * _filtered_wf_high->GetAmp(i_s) + alpha_high * (_bl_sub_wf->GetAmp(filt_samp  + 1) - _bl_sub_wf->GetAmp(filt_samp )));

    // Calculating derivative waveform
    if(i_s > derivative_offset) {
      _diff_wf->SetAmp(i_s,_filtered_wf_low->GetAmp(i_s) - _filtered_wf_low->GetAmp(i_s-derivative_offset));
    }

  }

  _raw_wf->SetExtrema(max_val_raw, max_samp_raw, min_val_raw, min_samp_raw);
  _bl_sub_wf->SetExtrema(max_val_sub, max_samp_sub, min_val_sub, min_samp_sub);
}

/*--------------------------------------------------------------------*/


void Channel::Normalize() {
  if(_normalized == true) return;
  if (_gain == 0) SetGain(); 
  if(_bl_sub_wf == NULL) {
    _raw_wf->CalcBaseline();
    _bl_sub_wf = _raw_wf->SubtractBaseline();
  }

  int num_samps = Config::Get()->GetParameterI("num_samps");
  for(int i_s = 0; i_s < num_samps; i_s++) {
    _bl_sub_wf->SetAmp(i_s, _bl_sub_wf->GetAmp(i_s) / _gain);
    if(_integral_wf != NULL) _integral_wf->SetAmp(i_s, _integral_wf->GetAmp(i_s) / _gain);
    if(_filtered_wf_low != NULL) _filtered_wf_low->SetAmp(i_s, _filtered_wf_low->GetAmp(i_s)/ _gain);
    if(_filtered_wf_high!= NULL) _filtered_wf_high->SetAmp(i_s, _filtered_wf_high->GetAmp(i_s)/ _gain);
    if(_filtered_wf_high2!= NULL) _filtered_wf_high2->SetAmp(i_s, _filtered_wf_high2->GetAmp(i_s)/ _gain);
    if(_filtered_wf_high3!= NULL) _filtered_wf_high3->SetAmp(i_s, _filtered_wf_high3->GetAmp(i_s)/ _gain);
    if(_filtered_wf_fft != NULL) _filtered_wf_fft->SetAmp(i_s, _filtered_wf_fft->GetAmp(i_s)/ _gain);
    if(_diff_wf != NULL) _diff_wf->SetAmp(i_s, _diff_wf->GetAmp(i_s) / _gain);
  }
  _bl_sub_wf->CalcBaseline();
  _bl_sub_wf->CalcExtrema();
  _total_integral = _bl_sub_wf->Integrate(0, num_samps);

  _normalized = true;
}

/*--------------------------------------------------------------------*/

void Channel::FindPulses() {
  int num_samps = Config::Get()->GetParameterI("num_samps");
  double sampling_rate = Config::Get()->GetParameterD("sampling_rate");
  int sample_size = Config::Get()->GetParameterI("pulse_end_size");
  double pulse_threshold = Config::Get()->GetParameterD("pulse_threshold");
  int pulse_width = Config::Get()->GetParameterI("pulse_width");

  for(int i_s = 0; i_s < num_samps; i_s++) {
    //Finding Main Pulse
    int pulse_start = i_s; int pulse_end = i_s; int peak_time = i_s;
    int half_samp;
    double half_time;
    double baseline_end;
    int temp_s;
    double min_val = _bl_sub_wf->GetAmp(i_s); double max_val = _bl_sub_wf->GetAmp(i_s);
    bool oob = false;
    if (_bl_sub_wf->GetAmp(i_s) >= pulse_threshold * _raw_wf->GetBaselineStdev(i_s)) {
      while(_bl_sub_wf->GetAmp(pulse_start) > 0) {pulse_start--;} // shifting backwards to find the start of the pulse
      while (_bl_sub_wf->GetAmp(peak_time+1) > _bl_sub_wf->GetAmp(peak_time)) peak_time++; //Finding peak time
      half_samp = pulse_start;
      while(_bl_sub_wf->GetAmp(half_samp) < 0.5 * _bl_sub_wf->GetAmp(peak_time)) half_samp++; //finding time of 1/2 peak height
      // Interpolating the true half-time
      double slope = (_bl_sub_wf->GetAmp(half_samp) - _bl_sub_wf->GetAmp(half_samp - 1)) / (1e-3 / sampling_rate);
      double intercept = _bl_sub_wf->GetAmp(half_samp) - (slope * half_samp * 1e-3 / sampling_rate);
      half_time = (0.5 * _bl_sub_wf->GetAmp(peak_time) - intercept) / slope;
      // Finding the end of the pulse
      while(true) {
	      pulse_end++;
	      if (pulse_end > num_samps) {oob = true; break;}
	      if(_bl_sub_wf->GetAmp(pulse_end) > max_val) {max_val = _bl_sub_wf->GetAmp(pulse_end);}
	      if(_bl_sub_wf->GetAmp(pulse_end) < min_val) {min_val = _bl_sub_wf->GetAmp(pulse_end);}
	      // End of pulse is fixed when waveform is below zero and the baseline over a given number of samples is within
	      // one sigma of the rolling baseline
	      if (_bl_sub_wf->GetAmp(pulse_end) <= 0){
	        baseline_end = 0;
	        temp_s = pulse_end;
	        for(int i = 0 ; i < sample_size; i++, temp_s++) {
	          if (temp_s > num_samps) {oob = true; break;}
	          baseline_end += _bl_sub_wf->GetAmp(temp_s)/sample_size;
	        }
	        if (abs(baseline_end) < _raw_wf->GetBaselineStdev(pulse_end)) break;
	      }
      }
      // If pulse length and integral exceed minimum, sets all pulse parameters
      float integral = _bl_sub_wf->Integrate(pulse_start, pulse_end);
      if (integral > max_val && pulse_end-pulse_start > pulse_width && !oob) {
	      Pulse* newPulse = new Pulse();
	      newPulse->SetStartSamp(pulse_start); newPulse->SetStartTime(pulse_start * 1e-3 / sampling_rate);
	      newPulse->SetPeakSamp(peak_time); newPulse->SetPeakTime(peak_time * 1e-3 / sampling_rate);
        newPulse->SetHalfTime(half_time);
	      newPulse->SetEndSamp(pulse_end); newPulse->SetMinVal(min_val); newPulse->SetMaxVal(max_val); newPulse->SetHeight(_bl_sub_wf->GetAmp(peak_time));
        newPulse->SetEndTime(pulse_end * 1e-3 / sampling_rate);
        newPulse->SetIntegral(integral);
	      _pulses.push_back(newPulse);
	      i_s = pulse_end;

        // Finding all pulses between the start-time and end-time of the pulse
	      FindAfterPulses(newPulse, pulse_start, pulse_end); 
      }
    }
  }
  _npulses = _pulses.size();
}


/*--------------------------------------------------------------------*/
    
void Channel::FindAfterPulses(Pulse* p, int pulse_start, int pulse_end) {
  int counter = 0;
  double pulse_threshold = Config::Get()->GetParameterD("pulse_threshold");
  double sampling_rate = Config::Get()->GetParameterD("sampling_rate");
  for(int i_s = pulse_start; i_s < pulse_end; i_s++) {
    int afterpulse_time = i_s; int afterpulse_end = i_s;
    // If derivative exceeds the threshold:
    if (_diff_wf->GetAmp(i_s) > (pulse_threshold * _raw_wf->GetBaselineStdev(i_s))) { 
      //Finding peak time at nearest maximum
      if (_bl_sub_wf->GetAmp(i_s + 1) > _bl_sub_wf->GetAmp(i_s - 1)) {while (_bl_sub_wf->GetAmp(afterpulse_time+1) > _bl_sub_wf->GetAmp(afterpulse_time-1)) afterpulse_time++;}
      else if (_bl_sub_wf->GetAmp(i_s + 1) < _bl_sub_wf->GetAmp(i_s-1)) {while (_bl_sub_wf->GetAmp(afterpulse_time+1) < _bl_sub_wf->GetAmp(afterpulse_time-1)) afterpulse_time--;}
      while (_diff_wf->GetAmp(afterpulse_end) > 0) afterpulse_end++; //advancing until decreasing edge of the pulse
      i_s = afterpulse_end;
      counter++;
      // Setting after-pulse parameters and adding to afterpulse array (excluding the first pulse)
      if (counter > 1) {
	      Pulse* newPulse = new Pulse();
	      newPulse->SetPeakSamp(afterpulse_time); newPulse->SetPeakTime(afterpulse_time * 1e-3 / sampling_rate);
	      newPulse->SetHeight(_bl_sub_wf->GetAmp(afterpulse_time));
	      p->AddAfterPulse(newPulse);
      }
    }
  }
}


/*--------------------------------------------------------------------*/

void Channel::PulseDeconv(Waveform* response) {
  if(_deconv_wf == NULL) SetDeconvWF(response);
  
  for(unsigned int i_p = 0; i_p < _pulses.size(); i_p++) {
    int deconv_start = _pulses[i_p]->GetStartSamp();
    while(_deconv_wf->GetAmp(deconv_start) > 0) deconv_start--; // Finding the time at which the deconvolved waveform crosses 0
    _pulses[i_p]->SetDeconvIntegral(_deconv_wf->Integrate(deconv_start, _pulses[i_p]->GetEndSamp()));
  }
}

/*--------------------------------------------------------------------*/

void Channel::SetDeconvWF(Waveform* response) {

  //Rescaling the response function
  if (_gain == 0) SetGain();
  Waveform * response_scale = new Waveform();
  double scale = (double)_gain / response->Integrate(0,  Config::Get()->GetParameterI("num_samps"));

  for(int i_s = 0 ; i_s < Config::Get()->GetParameterI("num_samps"); i_s++) {
    response_scale->SetAmp(i_s, response->GetAmp(i_s) * scale);
    if(_normalized == true) response_scale->SetAmp(i_s, response_scale->GetAmp(i_s) / _gain); 
  }

  _deconv_wf = _bl_sub_wf->Deconvolve(response_scale);
  _deconv_wf->CalcExtrema();

  response_scale->WaveformFree();
  delete response_scale;
}

/*--------------------------------------------------------------------*/

void Channel::ChannelPrint() {
  cout << "Total Integral: " << _total_integral << endl;

  if(Config::Get()->GetParameterS("print_pulses") == "y") {
    for (unsigned int i = 0; i < _pulses.size(); i++){
      cout << "Pulse " << i+1 << ": ";
      _pulses[i]->Print();
    }
  }
}


/*--------------------------------------------------------------------*/

void Channel::ChannelDisplay(double min_x, double max_x, double min_y, double max_y) {
  // Plotting Waveform
  TGraph *g = new TGraph(Config::Get()->GetParameterI("num_samps"),&_raw_wf->GetTimes()[0], &_bl_sub_wf->GetAmp()[0]);
  g->SetTitle(Form("Channel %d; Time [#mus]",_ch_id));
  if(_normalized == true) g->GetYaxis()->SetTitle(Form("PE / %g ns", 1 /  Config::Get()->GetParameterD("sampling_rate")));
  else  g->GetYaxis()->SetTitle("ADC Counts");
  g->GetXaxis()->SetRangeUser(min_x,max_x);
  if (min_y != 0 || max_y != 0) g->GetYaxis()->SetRangeUser(min_y, max_y);
  g->Draw();
  gPad->Update();

  if(Config::Get()->GetParameterS("plot_pulses") == "y") PlotPulseTimes();
  if(Config::Get()->GetParameterS("plot_baseline") == "y") PlotBaseline();
  if(Config::Get()->GetParameterS("plot_LP_filter") == "y") PlotLPFilter();
  if(Config::Get()->GetParameterS("plot_FFT_filter") == "y") PlotFFTFilter();
  if(Config::Get()->GetParameterS("plot_HP_filter") == "y") PlotHPFilter();
  if(Config::Get()->GetParameterS("plot_derivative") == "y") PlotDerivative();
  if(Config::Get()->GetParameterS("plot_integral") == "y") PlotIntegral();
  if(Config::Get()->GetParameterS("plot_deconv") == "y") PlotDeconv();

}

/*--------------------------------------------------------------------*/

void Channel::PlotPulseTimes() {
  double rate =  Config::Get()->GetParameterD("sampling_rate");
  for (unsigned int i = 0; i < _pulses.size(); i++) {
    TLine *pulse_start = new TLine(_pulses[i]->GetStartSamp() * 1e-3 / rate , _bl_sub_wf->GetMax(), _pulses[i]->GetStartSamp() * 1e-3/ rate, _bl_sub_wf->GetMin());
    pulse_start->SetLineColor(4);
    pulse_start->SetLineWidth(2);
    pulse_start->Draw("same");

    TLine *pulse_end = new TLine(_pulses[i]->GetEndSamp() * 1e-3 / rate , _bl_sub_wf->GetMax(), _pulses[i]->GetEndSamp() * 1e-3/ rate, _bl_sub_wf->GetMin());
    pulse_end->SetLineColor(4);
    pulse_end->SetLineWidth(2);
    pulse_end->Draw("same");

    TLine *pulse_time = new TLine(_pulses[i]->GetPeakTime() , _bl_sub_wf->GetMax(), _pulses[i]->GetPeakTime(), _bl_sub_wf->GetMin());
    pulse_time->SetLineColor(9);                                                                   
    pulse_time->Draw("same");

    TLine *half_time = new TLine(_pulses[i]->GetHalfTime() , _bl_sub_wf->GetMax(), _pulses[i]->GetHalfTime(), _bl_sub_wf->GetMin());
    half_time->SetLineColor(9);
    half_time->Draw("same");

    for (int j = 0; j < _pulses[i]->GetNAfterPulses(); j++) {
      TLine *after_pulse_start = new TLine(_pulses[i]->GetAfterPulse(j)->GetPeakTime(),
                                           _bl_sub_wf->GetMax(), _pulses[i]->GetAfterPulse(j)->GetPeakTime(), _bl_sub_wf->GetMin());
      after_pulse_start->SetLineColor(9);
      after_pulse_start->Draw("same");
    }
  }
}

/*--------------------------------------------------------------------*/

void Channel::PlotBaseline(double min_x, double max_x, double min_y, double max_y) { 
  TGraph *dev = new TGraph(Config::Get()->GetParameterI("num_samps"),&_raw_wf->GetTimes()[0], &_bl_sub_wf->GetBaselineStdev()[0]);
  dev->GetXaxis()->SetRangeUser(min_x,max_x);
  if (min_y != 0 || max_y != 0) dev->GetYaxis()->SetRangeUser(min_y, max_y);
  dev->SetLineColor(3);
  dev->Draw();
  gPad->Update();
}

/*--------------------------------------------------------------------*/

void Channel::PlotLPFilter(double min_x, double max_x, double min_y, double max_y) {
  TGraph *filter = new TGraph(Config::Get()->GetParameterI("num_samps"),&_raw_wf->GetTimes()[0], &_filtered_wf_low->GetAmp()[0]);
  filter->GetXaxis()->SetRangeUser(min_x,max_x);
  if (min_y != 0 || max_y != 0) filter->GetYaxis()->SetRangeUser(min_y, max_y);
  filter->SetLineColor(7);
  filter->Draw();
  gPad->Update();
}

/*--------------------------------------------------------------------*/

void Channel::PlotFFTFilter(double min_x, double max_x, double min_y, double max_y) {
  if(_filtered_wf_fft == NULL) _filtered_wf_fft = _bl_sub_wf->FFTFilter();

  TGraph *filter_fft = new TGraph(Config::Get()->GetParameterI("num_samps"),&_raw_wf->GetTimes()[0], &_filtered_wf_fft->GetAmp()[0]);
  filter_fft->GetXaxis()->SetRangeUser(min_x,max_x);
  if (min_y != 0 || max_y != 0) filter_fft->GetYaxis()->SetRangeUser(min_y, max_y);
  filter_fft->SetLineColor(51);
  filter_fft->Draw();
  gPad->Update();
}


/*--------------------------------------------------------------------*/

void Channel::PlotHPFilter(double min_x, double max_x, double min_y, double max_y) {

  //Plotting high-pass filtered wf
  TGraph *filter_high = new TGraph(Config::Get()->GetParameterI("num_samps"),&_raw_wf->GetTimes()[0], &_filtered_wf_high->GetAmp()[0]);
  filter_high->GetXaxis()->SetRangeUser(min_x,max_x);
  if (min_y != 0 || max_y != 0) filter_high->GetYaxis()->SetRangeUser(min_y, max_y);
  filter_high->SetLineColor(8);
  filter_high->Draw(); 
  gPad->Update();

  //Plotting high-pass filtered wf 2 
  TGraph *filter_high2 = new TGraph(Config::Get()->GetParameterI("num_samps"),&_raw_wf->GetTimes()[0], &_filtered_wf_high2->GetAmp()[0]);                                                         
  filter_high2->GetXaxis()->SetRangeUser(min_x,max_x);
  if (min_y != 0 || max_y != 0) filter_high2->GetYaxis()->SetRangeUser(min_y, max_y);
  filter_high2->SetLineColor(9);
  filter_high2->Draw();
  gPad->Update();

  //Plotting high-pass filtered wf 3                                                                                                                                                              
  TGraph *filter_high3 = new TGraph(Config::Get()->GetParameterI("num_samps"),&_raw_wf->GetTimes()[0], &_filtered_wf_high3->GetAmp()[0]);
  filter_high3->GetXaxis()->SetRangeUser(min_x,max_x);
  if (min_y != 0 || max_y != 0) filter_high3->GetYaxis()->SetRangeUser(min_y, max_y);
  filter_high3->SetLineColor(93);
  filter_high3->Draw();                                                                                                                                                                      
  gPad->Update();
}

/*--------------------------------------------------------------------*/

void Channel::PlotDerivative(double min_x, double max_x, double min_y, double max_y) {
  TGraph *diff = new TGraph(Config::Get()->GetParameterI("num_samps"),&_raw_wf->GetTimes()[0], &_diff_wf->GetAmp()[0]);
  diff->GetXaxis()->SetRangeUser(min_x,max_x);
  if (min_y != 0 || max_y != 0) diff->GetYaxis()->SetRangeUser(min_y, max_y);
  diff->SetLineColor(6);
  diff->Draw();
  gPad->Update();
}

/*--------------------------------------------------------------------*/

void Channel::PlotIntegral(double min_x, double max_x, double min_y, double max_y) {
  TGraph *gint = new TGraph(Config::Get()->GetParameterI("num_samps"),&_raw_wf->GetTimes()[0], &_integral_wf->GetAmp()[0]);
  Float_t rightmax = 1.1*TMath::MaxElement(Config::Get()->GetParameterI("num_samps"),gint->GetY());
  Float_t rightmin = 1.1*TMath::MinElement(Config::Get()->GetParameterI("num_samps"),gint->GetY());
  Float_t scale;
  TGaxis *axis;
  if (abs(rightmax) < abs(rightmin)) {
    scale = abs(gPad->GetUymin()/rightmin);
    axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
                      gPad->GetUxmax(), gPad->GetUymax(),rightmin,gPad->GetUymax()/scale,510,"+L");
  }
  else {
    scale = abs(gPad->GetUymax()/rightmax);
    axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
                      gPad->GetUxmax(), gPad->GetUymax(),gPad->GetUymin()/scale,rightmax,510,"+L");
  }

  gint->SetLineColor(2);
  for (int i=0;i<gint->GetN();i++) gint->GetY()[i] *= scale;
  gint->GetXaxis()->SetRangeUser(min_x,max_x);
  if (min_y != 0 || max_y != 0) gint->GetYaxis()->SetRangeUser(min_y, max_y);
  gint->SetLineWidth(2);
  gint->Draw("same");
  axis->SetLineColor(2);
  axis->SetLabelSize(0.03);
  if(_normalized == true) axis->SetTitle("PE");
  axis->SetLabelColor(2);
  axis->Draw();       
}

/*--------------------------------------------------------------------*/

void Channel::PlotDeconv(double min_x, double max_x, double min_y, double max_y) {

  if(_deconv_wf == NULL) {
    Waveform *response = new Waveform();
    TFile *f = new TFile(Config::Get()->GetParameterS("avgPulse_file").c_str());
    TGraph* avgPulse_input = (TGraph*)f->Get(Config::Get()->GetParameterS("avgPulse_name").c_str());
    for(int i_s = 0 ; i_s < Config::Get()->GetParameterI("pulse_avg_samps"); i_s++) {
      response->SetAmp(i_s, avgPulse_input->GetY()[Config::Get()->GetParameterI("pulse_avg_buf") + i_s]);
    }
    delete f;
    delete avgPulse_input;
    SetDeconvWF(response);
    response->WaveformFree();
    delete response;
  }

  TGraph *deconv = new TGraph(Config::Get()->GetParameterI("num_samps"),&_raw_wf->GetTimes()[0], &_deconv_wf->GetAmp()[0]);
  deconv->GetXaxis()->SetRangeUser(min_x,max_x);
  if (min_y != 0 || max_y != 0) deconv->GetYaxis()->SetRangeUser(min_y, max_y);
  deconv->SetLineColor(95);
  deconv->Draw("same");
  gPad->Update();  
}


/*--------------------------------------------------------------------*/

void Channel::ChannelFree() {
  if(_raw_wf!= NULL) _raw_wf->WaveformFree();
  if(_bl_sub_wf!= NULL) _bl_sub_wf->WaveformFree();
  if(_integral_wf!= NULL) _integral_wf->WaveformFree();
  if(_diff_wf!= NULL) _diff_wf->WaveformFree();
  if(_filtered_wf_low!= NULL) _filtered_wf_low->WaveformFree();
  if(_filtered_wf_high != NULL) _filtered_wf_high->WaveformFree();
  if(_filtered_wf_high2 != NULL) _filtered_wf_high2->WaveformFree();
  if(_filtered_wf_high3 != NULL) _filtered_wf_high3->WaveformFree();
  if(_deconv_wf != NULL) _deconv_wf->WaveformFree();
  if(_filtered_wf_fft!= NULL) _filtered_wf_fft->WaveformFree();

  delete _raw_wf;
  delete _bl_sub_wf;
  delete _integral_wf;
  delete _filtered_wf_low;
  delete _filtered_wf_high;
  delete _filtered_wf_high2;
  delete _filtered_wf_high3;
  delete _diff_wf;
  delete _deconv_wf;
  delete _filtered_wf_fft;

  for (unsigned int i = 0; i < _pulses.size(); i++) {
    _pulses[i]->PulseFree();
    delete _pulses[i];
  }
  _pulses.clear(); vector<Pulse*>().swap(_pulses);

}
