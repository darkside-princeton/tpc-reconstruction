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
#include "Channel.h"
#include "Output.h"


using namespace std;


/*--------------------------------------------------------------------*/

void Event::SetRawWF(int channel, Waveform* wf) {
  _channels[channel]->SetRawWF(wf);
}

/*--------------------------------------------------------------------*/

void Event::ProcessEvent() {
  for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){
    _channels[i_ch]->ProcessChannel();
    vector<Pulse*> pulses = _channels[i_ch]->GetPulses();
    Output::Get()->total_pulses[i_ch] += _channels[i_ch]->GetNPulses(); // Counting total pulses
    for(int i_p = 0; i_p < _channels[i_ch]->GetNPulses(); i_p++) {
      if(pulses[i_p]->GetNAfterPulses() == 0 &&
	    pulses[i_p]->GetIntegral() > Config::Get()->GetParameterD("min_SPE_charge") &&
	    pulses[i_p]->GetIntegral() < Config::Get()->GetParameterD("max_SPE_charge") &&
	    pulses[i_p]->GetHeight() > Config::Get()->GetParameterD("min_SPE_height") &&
	    pulses[i_p]->GetHeight() < Config::Get()->GetParameterD("max_SPE_height") &&
	    pulses[i_p]->GetEndTime() - pulses[i_p]->GetStartTime() > Config::Get()->GetParameterD("min_SPE_length") &&
	    pulses[i_p]->GetEndTime() - pulses[i_p]->GetStartTime() < Config::Get()->GetParameterD("max_SPE_length")) 
	    Output::Get()->total_SPE_pulses[i_ch]++; // Counting total SPE pulses
    }
  }
  
  if(Config::Get()->GetParameterS("TPC") == "y") {
    CalcSumWaveform();
    ClearClusters();
    FindClusters();
    FindClusterPulses();
    if(Config::Get()->GetParameterS("deconvolve") == "y") {
      Waveform *response = new Waveform();
      TFile *f = new TFile(Config::Get()->GetParameterS("avgPulse_file").c_str());
      TGraph* avgPulse_input = (TGraph*)f->Get(Config::Get()->GetParameterS("avgPulse_name").c_str());
      for(int i_s = 0 ; i_s < Config::Get()->GetParameterI("pulse_avg_samps"); i_s++) {
	      response->SetAmp(i_s, avgPulse_input->GetY()[Config::Get()->GetParameterI("pulse_avg_buf") + i_s]);
      }
      delete f;
      delete avgPulse_input;
      ClusterDeconv(response);
      response->WaveformFree();
      delete response;
      for(unsigned int i_c = 0; i_c < _clusters.size(); i_c++) {
	      if(_clusters[i_c]->GetS1() == true) Output::Get()->total_S1_pulses++; // Counting total S1 pulses
      }
    }
  }
}

/*--------------------------------------------------------------------*/

void Event::PulseOutput() {
  Output::Get()->eventID = _eventID;
  Output::Get()->npulses.clear();
  for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){
    vector<Pulse*> pulses = _channels[i_ch]->GetPulses();
    Output::Get()->start_time[i_ch].clear();
    Output::Get()->end_time[i_ch].clear();
    Output::Get()->integral[i_ch].clear();
    Output::Get()->max_val[i_ch].clear();
    Output::Get()->min_val[i_ch].clear();
    Output::Get()->height[i_ch].clear();
    Output::Get()->peak_time[i_ch].clear();
    Output::Get()->half_time[i_ch].clear();
    Output::Get()->nafterpulses[i_ch].clear();
    Output::Get()->deconv_integral[i_ch].clear();
    Output::Get()->filtered_integral[i_ch].clear();
    Output::Get()->ap_peak_times[i_ch].clear();
    Output::Get()->ap_heights[i_ch].clear(); 
    Output::Get()->npulses[i_ch] = pulses.size();
    for(int i = 0; i < _channels[i_ch]->GetNPulses(); i++) {
      Output::Get()->start_time[i_ch].push_back(pulses[i]->GetStartTime());
      Output::Get()->end_time[i_ch].push_back(pulses[i]->GetEndTime());
      Output::Get()->integral[i_ch].push_back(pulses[i]->GetIntegral());
      Output::Get()->max_val[i_ch].push_back(pulses[i]->GetMaxVal());
      Output::Get()->min_val[i_ch].push_back(pulses[i]->GetMinVal());
      Output::Get()->height[i_ch].push_back(pulses[i]->GetHeight());
      Output::Get()->peak_time[i_ch].push_back(pulses[i]->GetPeakTime());
      Output::Get()->half_time[i_ch].push_back(pulses[i]->GetHalfTime());
      Output::Get()->nafterpulses[i_ch].push_back(pulses[i]->GetNAfterPulses());
      Output::Get()->deconv_integral[i_ch].push_back(pulses[i]->GetDeconvIntegral());
      vector<double> ap_peak_time;
      vector<double> ap_height;
      for(int j = 0; j < pulses[i]->GetNAfterPulses(); j++) {
	      ap_peak_time.push_back(pulses[i]->GetAfterPulse(j)->GetPeakTime());
	      ap_height.push_back(pulses[i]->GetAfterPulse(j)->GetHeight());
      }
      Output::Get()->ap_peak_times[i_ch].push_back(ap_peak_time);
      Output::Get()->ap_heights[i_ch].push_back(ap_height);
    }
  }

  Output::Get()->_pulseTree->Fill();
}

/*--------------------------------------------------------------------*/

void Event::ClusterOutput() {
  Output::Get()->eventID = _eventID;
  Output::Get()->overlap = _cluster_overlap;
  Output::Get()->S1_pulses = 0;
  Output::Get()->S2_pulses = 0;
  Output::Get()->S1_charge = 0;
  Output::Get()->S2_charge = 0;
  Output::Get()->S1_S2_dt = -1;
  Output::Get()->nclusters = _clusters.size();
  Output::Get()->cl_start_time.clear();
  Output::Get()->cl_end_time.clear();
  Output::Get()->cl_integral.clear();
  Output::Get()->cl_max_val.clear();
  Output::Get()->cl_min_val.clear();
  Output::Get()->cl_height.clear();
  Output::Get()->cl_peak_time.clear();
  Output::Get()->cl_half_time.clear();
  Output::Get()->cl_symmetry.clear();
  for(unsigned int i_fp = 0; i_fp < Config::Get()->GetParameterDvec("prompt_windows").size(); i_fp++){
    Output::Get()->cl_fprompt[i_fp].clear();
    Output::Get()->cl_deconv_fprompt[i_fp].clear();
  }
  for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){
    Output::Get()->pulse_peak_time[i_ch].clear();
    Output::Get()->pulse_height[i_ch].clear();
    Output::Get()->cl_deconv_integral[i_ch].clear();
  }
  double first_s1_time = -1;
  double first_s2_time = -1;
  for(unsigned int i_c = 0; i_c < _clusters.size(); i_c++) {
    Output::Get()->cl_start_time.push_back(_clusters[i_c]->GetStartTime());
    Output::Get()->cl_end_time.push_back(_clusters[i_c]->GetEndTime());
    Output::Get()->cl_integral.push_back(_clusters[i_c]->GetIntegral());
    Output::Get()->cl_min_val.push_back(_clusters[i_c]->GetMinVal());
    Output::Get()->cl_max_val.push_back(_clusters[i_c]->GetMaxVal());
    Output::Get()->cl_height.push_back(_clusters[i_c]->GetHeight());
    Output::Get()->cl_peak_time.push_back(_clusters[i_c]->GetPeakTime());
    Output::Get()->cl_half_time.push_back(_clusters[i_c]->GetHalfTime());
    Output::Get()->cl_symmetry.push_back(_clusters[i_c]->GetSymmetry());
    if(_clusters[i_c]->GetS1() == true) {
      Output::Get()->S1_pulses++;
      Output::Get()->S1_charge += _clusters[i_c]->GetIntegral();
      if(first_s1_time < 0) first_s1_time = _clusters[i_c]->GetHalfTime();
    }
    else {
      Output::Get()->S2_pulses++;
      Output::Get()->S2_charge += _clusters[i_c]->GetIntegral();
      if(first_s2_time < 0) first_s2_time = _clusters[i_c]->GetHalfTime();
    } 
    for(unsigned int i_fp = 0; i_fp < Config::Get()->GetParameterDvec("prompt_windows").size(); i_fp++){
      Output::Get()->cl_fprompt[i_fp].push_back(_clusters[i_c]->GetFprompt(i_fp));
      Output::Get()->cl_deconv_fprompt[i_fp].push_back(_clusters[i_c]->GetDeconvFprompt(i_fp));
    }
    for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){
      vector<double> cl_pulse_peak_time;
      vector<double> cl_pulse_height;
      for(unsigned int i_p = 0; i_p < _clusters[i_c]->GetPulses(i_ch).size(); i_p++) {
	      cl_pulse_peak_time.push_back(_clusters[i_c]->GetPulse(i_ch, i_p)->GetPeakTime());
	      cl_pulse_height.push_back(_clusters[i_c]->GetPulse(i_ch, i_p)->GetHeight());
      }
      Output::Get()->pulse_peak_time[i_ch].push_back(cl_pulse_peak_time);
      Output::Get()->pulse_height[i_ch].push_back(cl_pulse_height);
      Output::Get()->cl_deconv_integral[i_ch].push_back(_clusters[i_c]->GetDeconvIntegral(i_ch));
    }
  }
  if(first_s2_time > first_s1_time) Output::Get()->S1_S2_dt = first_s2_time - first_s1_time;
  Output::Get()->_clusterTree->Fill();
}

/*--------------------------------------------------------------------*/

void Event::ProcessEventNoise() {
  for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){
    _channels[i_ch]->ProcessNoise();
    Output::Get()->baseline_mean = _channels[i_ch]->GetRawWF()->GetBaseline(Config::Get()->GetParameterI("num_samps") - Config::Get()->GetParameterI("baseline_samps")/2);
    Output::Get()->baseline_stdev = _channels[i_ch]->GetRawWF()->GetBaselineStdev(Config::Get()->GetParameterI("num_samps") - Config::Get()->GetParameterI("baseline_samps")/2);
    Output::Get()->eventID = _eventID;
    Output::Get()->_noiseTree[i_ch]->Fill();
  }
}

/*--------------------------------------------------------------------*/

void Event::CalcSumWaveform() {
  _sum_wf = new Waveform();
  for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){
    _channels[i_ch]->Normalize();
    for(int i_s = 0; i_s < Config::Get()->GetParameterI("num_samps"); i_s++){
      _sum_wf->SetAmp(i_s, _sum_wf->GetAmp(i_s) + _channels[i_ch]->GetBlSubWF()->GetAmp(i_s));
    }
  }
  _sum_wf->CalcBaseline();
  _sum_wf->CalcExtrema();
}

/*--------------------------------------------------------------------*/

void Event::FindClusters() {
  // Performs pulse-finding on the summed-waveform
  if(_sum_wf == NULL) CalcSumWaveform();

  for(int i_s = 0; i_s < Config::Get()->GetParameterI("num_samps"); i_s++) {
    int sample_size = Config::Get()->GetParameterI("cluster_end_size"); 
    int cluster_start = i_s; int cluster_end = i_s; int peak_time = i_s;
    int half_samp;
    double half_time;
    double baseline_end;
    int temp_s;
    double min_val = _sum_wf->GetAmp(i_s); double max_val = _sum_wf->GetAmp(i_s);
    bool oob = false;
    bool second_pulse = false;
    int pulse_match = 0;
    int baseline_samps = Config::Get()->GetParameterI("baseline_samps");
    if (_sum_wf->GetAmp(i_s) >= Config::Get()->GetParameterD("cluster_threshold") * _sum_wf->GetBaselineStdev(max(i_s-baseline_samps, baseline_samps/2))) {
      while(_sum_wf->GetAmp(cluster_start) > 0) {cluster_start--;} // determining cluster start time
      while (_sum_wf->GetAmp(peak_time+1) > _sum_wf->GetAmp(peak_time)) peak_time++; // determining cluster peak time
      // Interpolating time at half-height
      half_samp = cluster_start;
      while(_sum_wf->GetAmp(half_samp) < 0.5 * _sum_wf->GetAmp(peak_time)) half_samp++;
      double slope = (_sum_wf->GetAmp(half_samp) - _sum_wf->GetAmp(half_samp - 1)) / (1e-3 / Config::Get()->GetParameterD("sampling_rate"));
      double intercept = _sum_wf->GetAmp(half_samp) - (slope * half_samp * 1e-3 / Config::Get()->GetParameterD("sampling_rate"));
      half_time = (0.5 * _sum_wf->GetAmp(peak_time) - intercept) / slope;
      while(true) {
        cluster_end++;
        if (cluster_end > Config::Get()->GetParameterI("num_samps")) {oob = true; break;}
        if(_sum_wf->GetAmp(cluster_end) > max_val) {max_val = _sum_wf->GetAmp(cluster_end);}
        if(_sum_wf->GetAmp(cluster_end) < min_val) {min_val = _sum_wf->GetAmp(cluster_end);}
	      // Checking for pulse overlap (NEEDS WORK!!)
	      int second_pulse_length = 0;
	      if (_sum_wf->GetAmp(cluster_end) - min_val >= Config::Get()->GetParameterD("double_cluster_threshold") * _sum_wf->GetBaselineStdev(max(i_s-baseline_samps, baseline_samps/2))) {
	        while(cluster_end + second_pulse_length < Config::Get()->GetParameterI("num_samps") - 1 &&_sum_wf->GetAmp(cluster_end + second_pulse_length) - min_val >= 
		      Config::Get()->GetParameterD("double_cluster_threshold") * _sum_wf->GetBaselineStdev(max(i_s-baseline_samps, baseline_samps/2))) second_pulse_length++;
	      }
	      if(second_pulse_length > Config::Get()->GetParameterI("double_cluster_width")) second_pulse = true; 
	      if (_sum_wf->GetAmp(cluster_end) <= 0){
	        baseline_end = 0;
	        temp_s = cluster_end;
	        for(int i = 0 ; i < sample_size; i++, temp_s++) {
	          if (temp_s > Config::Get()->GetParameterI("num_samps")) {oob = true; break;}
	          baseline_end += _sum_wf->GetAmp(temp_s)/sample_size;
	        }
	        if (abs(baseline_end) < 1 * _sum_wf->GetBaselineStdev(cluster_end)) {i_s = cluster_end; cluster_end = temp_s; break; } // cluster_end = temp_s is buffer
	      }
      }
      // Checking if each channel has a corresponding pulse in the cluster region
      for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){
	      if(_channels[i_ch]->GetBlSubWF()->Integrate(cluster_start, cluster_end) >= 1) pulse_match++;
      }	
      // Setting all cluster parameters
      float integral = _sum_wf->Integrate(cluster_start, cluster_end);
      if (second_pulse) _cluster_overlap = true;
      if (pulse_match >=  Config::Get()->GetParameterI("num_chans") / 2 && integral > max_val 
	      && cluster_end-cluster_start > Config::Get()->GetParameterI("cluster_width") && !oob) {
        Cluster* newCluster = new Cluster();
        newCluster->SetStartSamp(cluster_start); newCluster->SetStartTime(cluster_start * 1e-3 / Config::Get()->GetParameterD("sampling_rate"));
        newCluster->SetPeakSamp(peak_time); newCluster->SetPeakTime(peak_time * 1e-3 / Config::Get()->GetParameterD("sampling_rate"));
        newCluster->SetHalfTime(half_time);
        newCluster->SetEndSamp(cluster_end); newCluster->SetMinVal(min_val); newCluster->SetMaxVal(max_val); newCluster->SetHeight(_sum_wf->GetAmp(peak_time));
        newCluster->SetEndTime(cluster_end * 1e-3 / Config::Get()->GetParameterD("sampling_rate"));
        newCluster->SetIntegral(integral);
	      // Calculating fprompt values
	      vector<double> fprompts;
	      for(unsigned int i_fp = 0; i_fp < Config::Get()->GetParameterDvec("prompt_windows").size(); i_fp++){
	        double fp = _sum_wf->Integrate(newCluster->GetStartTime(), newCluster->GetStartTime() + Config::Get()->GetParameterDvec("prompt_windows")[i_fp]);
	        if(integral != 0) fprompts.push_back(fp / integral);
	        else fprompts.push_back(0);
	      }
	      newCluster->SetFprompt(fprompts);
	      if(Config::Get()->GetParameterS("deconvolve") == "n") {
	        if(fprompts[0] > Config::Get()->GetParameterD("S1_S2_threshold")) newCluster->SetS1(true);
	        else newCluster->SetS1(false);
	      }
	      // Calculating asymmetry
	      double top_frac = 0;
	      for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){
	        if(_channels[i_ch]->GetTop() == true) top_frac += _channels[i_ch]->GetBlSubWF()->Integrate(cluster_start, cluster_end);
	      }
	      if(integral != 0) newCluster->SetSymmetry(top_frac/integral);
	      else newCluster->SetSymmetry(0);
	      _clusters.push_back(newCluster);
        i_s = cluster_end;
      }
    }
  }
}

/*--------------------------------------------------------------------*/

void Event::FindClusterPulses() {
  // Finding pulses between cluster start and end for each channel
  for(unsigned int i_c = 0; i_c < _clusters.size(); i_c++) {
    _clusters[i_c]->ClearPulses();
    for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){
      for(int i_s = _clusters[i_c]->GetStartSamp(); i_s < _clusters[i_c]->GetEndSamp(); i_s++) {
	      int pulse_time = i_s; int pulse_end = i_s;
	      if (_channels[i_ch]->GetDiffWF()->GetAmp(i_s) > (Config::Get()->GetParameterD("pulse_threshold") * _channels[i_ch]->GetBlSubWF()->GetBaselineStdev(i_s))) {
	        if (_channels[i_ch]->GetBlSubWF()->GetAmp(i_s + 1) > _channels[i_ch]->GetBlSubWF()->GetAmp(i_s - 1)) {
	          while (_channels[i_ch]->GetBlSubWF()->GetAmp(pulse_time+1) > _channels[i_ch]->GetBlSubWF()->GetAmp(pulse_time-1)) pulse_time++;
	        }
	        else if (_channels[i_ch]->GetBlSubWF()->GetAmp(i_s + 1) < _channels[i_ch]->GetBlSubWF()->GetAmp(i_s-1)) {
	          while (_channels[i_ch]->GetBlSubWF()->GetAmp(pulse_time+1) < _channels[i_ch]->GetBlSubWF()->GetAmp(pulse_time-1)) pulse_time--;
          }
	        while (_channels[i_ch]->GetDiffWF()->GetAmp(pulse_end) > 0) pulse_end++; //advancing until decreasing edge of the pulse
	        i_s = pulse_end;
	        Pulse* newPulse = new Pulse();
	        newPulse->SetPeakSamp(pulse_time); newPulse->SetPeakTime(pulse_time * 1e-3 / Config::Get()->GetParameterD("sampling_rate"));
	        newPulse->SetHeight(_channels[i_ch]->GetBlSubWF()->GetAmp(pulse_time));
	        _clusters[i_c]->AddPulse(i_ch, newPulse);
	      }
      }
    }
  }
}

/*--------------------------------------------------------------------*/

void Event::ClusterDeconv(Waveform* response) {
  
  //Setting up fprompt calculation
  vector<double> total_integral;
  vector<vector<double> > total_fprompts;

  for(unsigned int i_c = 0; i_c < _clusters.size(); i_c++) {
    total_integral.push_back(0);
    vector<double> cl_fprompts;
    total_fprompts.push_back(cl_fprompts);
    for(unsigned int i_fp = 0; i_fp < Config::Get()->GetParameterDvec("prompt_windows").size(); i_fp++){
      total_fprompts[i_c].push_back(0);
    }
  }

  for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){
    _channels[i_ch]->SetDeconvWF(response);
    for(unsigned int i_c = 0; i_c < _clusters.size(); i_c++) {
      int deconv_start = _clusters[i_c]->GetStartSamp();
      while(_channels[i_ch]->GetDeconvWF()->GetAmp(deconv_start) > 0) deconv_start--; // finding start-time of deconvolved waveform
     
      _clusters[i_c]->SetDeconvIntegral(i_ch, _channels[i_ch]->GetDeconvWF()->Integrate(deconv_start, _clusters[i_c]->GetEndSamp()));
      total_integral[i_c] += _channels[i_ch]->GetDeconvWF()->Integrate(deconv_start, _clusters[i_c]->GetEndSamp());

      // Summming integrals in prompt-windows of each channel
      for(unsigned int i_fp = 0; i_fp < Config::Get()->GetParameterDvec("prompt_windows").size(); i_fp++){
	double fp = _channels[i_ch]->GetDeconvWF()->Integrate(deconv_start * 1e-3 / Config::Get()->GetParameterD("sampling_rate") ,
							      deconv_start * 1e-3 / Config::Get()->GetParameterD("sampling_rate") + Config::Get()->GetParameterDvec("prompt_windows")[i_fp]);
	total_fprompts[i_c][i_fp] += fp;
      }
    }
  }
 
  // Calculating total fprompt
  for(unsigned int i_c = 0; i_c < _clusters.size(); i_c++) {
    for(unsigned int i_fp = 0; i_fp < Config::Get()->GetParameterDvec("prompt_windows").size(); i_fp++){
      total_fprompts[i_c][i_fp] /= total_integral[i_c];
    }
    _clusters[i_c]->SetDeconvFprompts(total_fprompts[i_c]);
    
    // Classifying S1/S2
    if(total_fprompts[i_c][0] > Config::Get()->GetParameterD("S1_S2_threshold")) _clusters[i_c]->SetS1(true);
    else _clusters[i_c]->SetS1(false);
  }
}


/*--------------------------------------------------------------------*/

void Event::EventPrint(int channel){
  // cout << "EVENT ID: " << _eventID << endl;
  //cout << "====================================================================" << endl;
  //cout << "====================================================================" << endl;

  if(Config::Get()->GetParameterS("plot_deconv") == "y" && _channels[channel]->GetDeconvWF() == NULL) {
    Waveform *response = new Waveform();
    TFile *f = new TFile(Config::Get()->GetParameterS("avgPulse_file").c_str());
    TGraph* avgPulse_input = (TGraph*)f->Get(Config::Get()->GetParameterS("avgPulse_name").c_str());
    for(int i_s = 0 ; i_s < Config::Get()->GetParameterI("pulse_avg_samps"); i_s++) {
      response->SetAmp(i_s, avgPulse_input->GetY()[Config::Get()->GetParameterI("pulse_avg_buf") + i_s]);
    }
    delete f;
    delete avgPulse_input;
    ClusterDeconv(response);
    response->WaveformFree();
    delete response;
  }

  if(channel == 0) {
    for (unsigned int i = 0; i < _clusters.size(); i++){
      cout << "Cluster " << i+1 << ": ";
      _clusters[i]->Print();
    }
  }

  cout << "====================================================================" << endl;
  cout << "Channel " << channel << endl;
  _channels[channel]->ChannelPrint();

}

/*--------------------------------------------------------------------*/

void Event::DisplayEvent(int channel, double min_x, double max_x, double min_y, double max_y) {
  _channels[channel]->ChannelDisplay(min_x,max_x, min_y, max_y);

  if(Config::Get()->GetParameterS("plot_cluster_times") == "y") PlotClusterTimes(channel, min_x, max_x, min_y, max_y);
  if(Config::Get()->GetParameterS("plot_cluster_pulses") == "y") PlotClusterPulses(channel, min_x, max_x, min_y, max_y);
  

  if(Config::Get()->GetParameterS("plot_sum_wf") == "y" && Config::Get()->GetParameterS("TPC") == "y") {
    if(channel == 0) {
      TCanvas* cSum = new TCanvas("Sum Waveform", "Sum Waveform", 1400, 800);
      cSum->cd();
      PlotSumWaveform(min_x, max_x, min_y, max_y);
    }
  }
}


/*--------------------------------------------------------------------*/

void Event::PlotClusterTimes(int channel, double min_x, double max_x, double min_y, double max_y) {
  TLine *cluster_start;
  TLine *cluster_end;
  for (unsigned int i = 0; i < _clusters.size(); i++) {
    cluster_start = new TLine(_clusters[i]->GetStartTime() , _channels[channel]->GetBlSubWF()->GetMax(), _clusters[i]->GetStartTime(),  _channels[channel]->GetBlSubWF()->GetMin());
    cluster_start->SetLineColor(4);
    cluster_start->SetLineWidth(2);
    cluster_start->Draw("same");

    cluster_end = new TLine(_clusters[i]->GetEndTime() ,  _channels[channel]->GetBlSubWF()->GetMax(), _clusters[i]->GetEndTime(),  _channels[channel]->GetBlSubWF()->GetMin());
    cluster_end->SetLineColor(4);
    cluster_end->SetLineWidth(2);
    cluster_end->Draw("same");
  }
  gPad->Update();
}

/*--------------------------------------------------------------------*/

void Event::PlotClusterPulses(int channel, double min_x, double max_x, double min_y, double max_y) {
  for (unsigned int i = 0; i < _clusters.size(); i++) {    
    for (unsigned int i_p = 0; i_p < _clusters[i]->GetPulses(channel).size(); i_p++) {
      TLine *pulse_start = new TLine(_clusters[i]->GetPulse(channel,i_p)->GetPeakTime(),
                                     _channels[channel]->GetBlSubWF()->GetMax(), _clusters[i]->GetPulse(channel,i_p)->GetPeakTime(),  _channels[channel]->GetBlSubWF()->GetMin());
      pulse_start->SetLineColor(9);
      pulse_start->Draw("same");
    }
  }
}

/*--------------------------------------------------------------------*/

void Event::PlotSumWaveform(double min_x, double max_x, double min_y, double max_y) {
  TGraph *gSum = new TGraph(Config::Get()->GetParameterI("num_samps"),&_sum_wf->GetTimes()[0], &_sum_wf->GetAmp()[0]);
  gSum->SetTitle(Form("Sum Waveform; Time [#mus]; PE / %g ns", 1 / Config::Get()->GetParameterD("sampling_rate")));
  gSum->GetXaxis()->SetRangeUser(min_x,max_x);
  if (min_y != 0 || max_y != 0) gSum->GetYaxis()->SetRangeUser(min_y, max_y);
  gSum->Draw();

  if(Config::Get()->GetParameterS("plot_baseline") == "y") {
    TGraph *dev = new TGraph(Config::Get()->GetParameterI("num_samps"),&_sum_wf->GetTimes()[0], &_sum_wf->GetBaselineStdev()[0]);
    dev->GetXaxis()->SetRangeUser(min_x,max_x);
    if (min_y != 0 || max_y != 0) dev->GetYaxis()->SetRangeUser(min_y, max_y);
    dev->SetLineColor(3);
    dev->Draw("same");
    gPad->Update();
  }

  if(Config::Get()->GetParameterS("plot_cluster_times") == "y") {
    TLine *cluster_start;
    TLine *cluster_end;
    for (unsigned int i = 0; i < _clusters.size(); i++) {                                                                                                
      cluster_start = new TLine(_clusters[i]->GetStartTime() , _sum_wf->GetMax(), _clusters[i]->GetStartTime(),  _sum_wf->GetMin());
      cluster_start->SetLineColor(4);
      cluster_start->SetLineWidth(2);
      cluster_start->Draw("same");

      cluster_end = new TLine(_clusters[i]->GetEndTime() ,  _sum_wf->GetMax(), _clusters[i]->GetEndTime(),  _sum_wf->GetMin());
      cluster_end->SetLineColor(4);
      cluster_end->SetLineWidth(2);
      cluster_end->Draw("same");
    }
    gPad->Update();
  }

  if(Config::Get()->GetParameterS("plot_integral") == "y") {
    Waveform* _sum_int_wf = _sum_wf->CalcRollingIntegral();
    TGraph *gint = new TGraph(Config::Get()->GetParameterI("num_samps"),&_sum_int_wf->GetTimes()[0], &_sum_int_wf->GetAmp()[0]);    
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
    axis->SetLabelColor(2);
    axis->SetLabelSize(0.03);
    axis->SetTitle("PE");
    axis->Draw();
  }
}


/*--------------------------------------------------------------------*/

void Event::EventFree() {
  if(_sum_wf!=NULL) _sum_wf->WaveformFree();
  delete _sum_wf;
  for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){
    _channels[i_ch]->ChannelFree();
    delete _channels[i_ch];
  }
  for(unsigned int i_c = 0; i_c < _clusters.size(); i_c++) {
    _clusters[i_c]->ClusterFree();
    delete _clusters[i_c];
  }

}
