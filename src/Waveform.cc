#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <climits>
#include <vector>
#include <deque>
#include "TMath.h"
#include "Waveform.h"
#include "TVirtualFFT.h"
#include "TSpectrum.h"
#include "TComplex.h"
#include "TFile.h"
#include "TGraph.h"
#include "TF1.h"

using namespace std;


/*--------------------------------------------------------------------*/

void Waveform::CalcBaseline(double threshold) {
  int sample_size = Config::Get()->GetParameterI("baseline_samps");
  int num_samps = Config::Get()->GetParameterI("num_samps");
  deque<int> used_samps; //array to hold samples used to calculate the baseline

  if (Config::Get()->GetParameterS("constant_baseline") == "n") {

    //Initializing first sample
    double avg = 0, nrolled = 0, stdev = 0;
    float minval = INT_MAX, maxval = INT_MIN;

    // Calculates an average over the first sample_size/2 samples
    for(int i_roll = 0; i_roll <= sample_size/2; i_roll++) {
      avg += _amplitudes[i_roll];
      stdev += pow(_amplitudes[i_roll],2);
      if (_amplitudes[i_roll] > maxval) maxval = _amplitudes[i_roll];
      if (_amplitudes[i_roll] < minval) minval = _amplitudes[i_roll]; 
      nrolled++;
      used_samps.push_back(i_roll);
    }

    _baseline[0] = (float)(avg / nrolled);
    _baseline_stdev[0] = (float)(sqrt(stdev / nrolled - pow(_baseline[0],2)));
    _baseline_range[0] = (float)(maxval - minval);


    //Now for the rest of the samples
    int peak_start = num_samps;
    int peak_end = 0;
    // First sample_size/2 samples
    for(int i_s = 1; i_s < num_samps; i_s++) {
      if(i_s < sample_size/2) {
        avg += _amplitudes[used_samps.back()+1];
        stdev += pow(_amplitudes[used_samps.back()+1],2);
        maxval = max(maxval, _amplitudes[used_samps.back()+1]);
        minval = min(minval, _amplitudes[used_samps.back()+1]);
        nrolled++;
        used_samps.push_back(used_samps.back()+1);
      }
      // Last sample_size/2 samples
      else if(i_s > num_samps - sample_size/2 -1) {
        avg -= _amplitudes[used_samps[0]];
        stdev -= pow(_amplitudes[used_samps[0]],2);
        nrolled--;
        if(_amplitudes[used_samps[0]] == maxval || _amplitudes[used_samps[0]] == minval) {
	        maxval = _amplitudes[used_samps[1]]; minval =_amplitudes[used_samps[1]];
	        for(unsigned int i = 1; i < used_samps.size(); i++) {
	          if(_amplitudes[used_samps[i]] > maxval) maxval = _amplitudes[used_samps[i]];
	          if(_amplitudes[used_samps[i]] < minval) minval = _amplitudes[used_samps[i]];
	        }
        }
        used_samps.pop_front();
      }
      // All other samples
      else {
        int newsamp = used_samps[sample_size-1]+1;
        int buff = 5;
        // Ignoring samples that exceed a threshold * baseline_stdev, where the standard deviation is calculated from the 
        // previous non-overlapping sample_size range
        if(abs(_amplitudes[newsamp + buff] - _baseline[i_s- sample_size/2]) > threshold * _baseline_stdev[i_s-sample_size/2]){
	        peak_start = newsamp;
          // If positively exceeding threshold, continue until first sample is below baseline
	        if(_amplitudes[newsamp + buff] - _baseline[i_s- sample_size/2] > threshold * _baseline_stdev[i_s-sample_size/2]) { 
	          while(true) {
	            newsamp++;
	            if (newsamp + buff > num_samps) break;
	            // End of pulse is fixed when waveform is below the baseline and the baseline over a given number of samples is within
	            // one sigma of the rolling baseline
	            if (_amplitudes[newsamp + buff] < _baseline[i_s-sample_size/2] ){
	              double baseline_end = 0; int temp_s;
	              temp_s = newsamp + buff;
	              for(int i = 0 ; i < sample_size; i++, temp_s++) {
	  	            if (temp_s > num_samps) break;
	  	            baseline_end += _amplitudes[temp_s]/sample_size;
	              }
	              if (abs(baseline_end -  _baseline[i_s- sample_size/2]) < _baseline_stdev[i_s-sample_size/2]) break;
	            }
	          }
	        }
          // If negatively exceeding threshold, continue until first sample is above baseline
	        else if (_amplitudes[newsamp + buff] - _baseline[i_s-sample_size/2] < threshold * _baseline_stdev[i_s-sample_size/2]) {
	          while(true) {
              newsamp++;
              if (newsamp + buff > num_samps) break;
              // End of pulse is fixed when waveform is below the baseline and the baseline over a given number of samples is within
              // one sigma of the rolling baseline
              if (_amplitudes[newsamp + buff] > _baseline[i_s-sample_size/2] ){
                double baseline_end = 0; int temp_s;
                temp_s = newsamp + buff;
                for(int i = 0 ; i < sample_size; i++, temp_s++) {
                  if (temp_s > num_samps) break;
                  baseline_end += _amplitudes[temp_s]/sample_size;
                }
                if (abs(baseline_end -  _baseline[i_s- sample_size/2]) < _baseline_stdev[i_s-sample_size/2]) break;
              }
            }
	        }
	        newsamp += buff;
	        peak_end = newsamp;
        }
        if((i_s < peak_start || i_s > peak_end) && newsamp < num_samps) {
	        avg = avg - _amplitudes[used_samps[0]] + _amplitudes[newsamp]; // subtracting first sample and adding new sample
	        stdev = stdev - pow(_amplitudes[used_samps[0]],2) + pow(_amplitudes[newsamp],2);
	        used_samps.push_back(newsamp); // updating used samples list
	        if(_amplitudes[newsamp] > maxval) maxval = _amplitudes[newsamp];
	        if(_amplitudes[newsamp] < minval) minval = _amplitudes[newsamp];
	        else if(_amplitudes[used_samps[0]] == maxval || _amplitudes[used_samps[0]] == minval) {
	          maxval = _amplitudes[used_samps[1]]; minval = _amplitudes[used_samps[1]];
	          for(int i = 1; i <= sample_size; i++) {
	            if(_amplitudes[used_samps[i]] > maxval) maxval = _amplitudes[used_samps[i]];
	            if(_amplitudes[used_samps[i]] < minval) minval = _amplitudes[used_samps[i]];
	          }
	        }
	        used_samps.pop_front();
        }
      }
      _baseline[i_s] = (float)(avg / nrolled);
      _baseline_stdev[i_s] = (float)(sqrt(stdev / nrolled - pow(_baseline[i_s],2)));
      _baseline_range[i_s] = (float)(maxval - minval);
    }
  }
}


/*--------------------------------------------------------------------*/


Waveform* Waveform::SubtractBaseline(){
  Waveform* subtracted_wf = new Waveform();
  int num_samps = Config::Get()->GetParameterI("num_samps");
  //if(_baseline[0] == 0) CalcBaseline();
  for(int i_s = 0; i_s < num_samps; i_s++){
    subtracted_wf->SetAmp(i_s, _amplitudes[i_s] -_baseline[i_s]);
  }
  return subtracted_wf;
}


/*--------------------------------------------------------------------*/


double Waveform::Integrate(int s0, int s1){
  double sum = 0;
  s0 = max(s0,0);
  s1 = min(s1,Config::Get()->GetParameterI("num_samps"));
  for(int i_s = s0; i_s < s1; i_s++){
    sum += _amplitudes[i_s];
  }
  return sum;
}

/*--------------------------------------------------------------------*/

double Waveform::Integrate(double t0, double t1) {
  double sum = 0;
  int s0 = t0*1e3*Config::Get()->GetParameterD("sampling_rate");
  int s1 = t1*1e3*Config::Get()->GetParameterD("sampling_rate");
  sum = Integrate(int(s0),int(s1));
  return sum;
}


/*--------------------------------------------------------------------*/


void Waveform::CalcExtrema(){
  _max = INT_MIN;
  _min = INT_MAX;
  int num_samps = Config::Get()->GetParameterI("num_samps");
  for(int i_s = 0; i_s < num_samps; i_s++){
    if(_amplitudes[i_s] > _max){
      _max = _amplitudes[i_s];
      _max_samp = i_s;
    }
    if(_amplitudes[i_s] < _min){
      _min = _amplitudes[i_s];
      _min_samp = i_s;
    }
  }
}


/*--------------------------------------------------------------------*/


Waveform* Waveform::CalcRollingIntegral(){
  Waveform* roll_int = new Waveform();
  int num_samps = Config::Get()->GetParameterI("num_samps");
  roll_int->SetAmp(0,0);
  for(int i_s = 1; i_s < num_samps; i_s++){
    roll_int->SetAmp(i_s, _amplitudes[i_s] + roll_int->GetAmp(i_s-1));
  }
  return roll_int;
}


/*--------------------------------------------------------------------*/


Waveform* Waveform::LowPassFilter(){
  Waveform* filtered = new Waveform();
  double cutoff = Config::Get()->GetParameterD("cutoff_freq_low");
  double RC = 1/ (2 * 3.1415 * cutoff);
  double sampling_rate = Config::Get()->GetParameterD("sampling_rate");
  int num_samps = Config::Get()->GetParameterI("num_samps");
  double alpha = (1e-3 /  sampling_rate) / (RC + 1e-3 /  sampling_rate);
  // double beta = exp(-cutoff * 1e-3 / Config::Get()->GetParameterD("sampling_rate"));
  filtered->SetAmp(0, _amplitudes[0]);
  for(int i_s = 0; i_s < num_samps-1; i_s++){
    filtered->SetAmp(i_s + 1, (1-alpha) * filtered->GetAmp(i_s) + alpha * _amplitudes[i_s]);  
  }
  return filtered;
}


/*--------------------------------------------------------------------*/

Waveform* Waveform::HighPassFilter(){
  Waveform* filtered = new Waveform();
  double cutoff = Config::Get()->GetParameterD("cutoff_freq_high");
  double RC = 1/ (2 * 3.1415 * cutoff);
  double alpha = RC / (RC + (1e-3 /  Config::Get()->GetParameterD("sampling_rate")));
  int num_samps = Config::Get()->GetParameterI("num_samps");
  filtered->SetAmp(0, _amplitudes[0]);
  for(int i_s = 0; i_s < num_samps-1; i_s++){
    filtered->SetAmp(i_s + 1, alpha * filtered->GetAmp(i_s) + alpha * (_amplitudes[i_s + 1] - _amplitudes[i_s]));
  }
  return filtered;
}


/*--------------------------------------------------------------------*/


void Waveform::CalcFFT() {

  TH1F* hRawWF = new TH1F("hRawWF", "hRawWF", Config::Get()->GetParameterI("num_samps") + 1, 0, 
			  Config::Get()->GetParameterI("num_samps")*1e-3 / Config::Get()->GetParameterD("sampling_rate"));

  for(int i_s = 0 ; i_s < Config::Get()->GetParameterI("num_samps"); i_s++){
    hRawWF->SetBinContent(i_s+1, _amplitudes[i_s]);
  }

  TH1 *hFFT = 0;
  TVirtualFFT::SetTransform(0);
  hFFT = hRawWF->FFT(hFFT, "MAG");
  for (int bin = 1; bin <= Config::Get()->GetParameterI("num_samps")/2; bin++){
    _FFTMag.push_back(2*hFFT->GetBinContent(bin)/ Config::Get()->GetParameterI("num_samps"));
  }
  delete hFFT;
  delete hRawWF;
}


/*--------------------------------------------------------------------*/

Waveform* Waveform::Deconvolve(Waveform* response) {
  Waveform* deconv_wf = new Waveform();
  
  TH1D* hResponse = new TH1D("hResponse", "hResponse", Config::Get()->GetParameterI("num_samps"), 0, 
			     Config::Get()->GetParameterI("num_samps") * 1e-3 / Config::Get()->GetParameterD("sampling_rate"));
  TH1D* hSource = new TH1D("hSource", "hSource", Config::Get()->GetParameterI("num_samps"), 0, 
			   Config::Get()->GetParameterI("num_samps") * 1e-3 / Config::Get()->GetParameterD("sampling_rate"));

  //Setting up histograms
  for (int i = 0; i < Config::Get()->GetParameterI("num_samps"); i++) { 
    hResponse->SetBinContent(i + 1, response->GetAmp(i));
    hSource->SetBinContent(i + 1, (double)_amplitudes[i]);
  }

  //Calculating FFT for source and response functions
  int N = Config::Get()->GetParameterI("num_samps");
  TH1 *hresponse_m = 0;
  TVirtualFFT::SetTransform(0);
  hresponse_m = hResponse->FFT(hresponse_m, "RE");
  TVirtualFFT *fft_response = TVirtualFFT::GetCurrentTransform();
  delete hresponse_m;
  TH1 *hsource_m = 0;
  TVirtualFFT::SetTransform(0);
  hsource_m = hSource->FFT(hsource_m, "RE");
  TVirtualFFT *fft_source = TVirtualFFT::GetCurrentTransform();
  delete hsource_m;

  //Setting up IFFT
  TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &N, "C2R K");
  for(int i_s = 0; i_s < Config::Get()->GetParameterI("num_samps"); i_s++) {
    Double_t re, im;
    fft_response->GetPointComplex(i_s, re, im);
    TComplex cResponse(re, im); cResponse += 1e-30;
    fft_source->GetPointComplex(i_s, re, im);
    TComplex cSource(re, im); 
    TComplex div = cSource / cResponse;
    fft_back->SetPointComplex(i_s, div);
    // Filtering
    if(i_s >  Config::Get()->GetParameterD("fft_filter_freq") * Config::Get()->GetParameterI("num_samps") / (Config::Get()->GetParameterD("sampling_rate") * 1e3) && 
       i_s < Config::Get()->GetParameterI("num_samps") - 
       Config::Get()->GetParameterD("fft_filter_freq") * Config::Get()->GetParameterI("num_samps") / (Config::Get()->GetParameterD("sampling_rate") * 1e3)) {
      TComplex filter(0,0);
      fft_back->SetPointComplex(i_s, filter);
    }
  }

  fft_back->Transform();
  
  
  TH1 *hResult = 0;
  hResult = TH1::TransformHisto(fft_back, hResult,"Re");
  hResult->Scale((double)1/ Config::Get()->GetParameterI("num_samps"));

  for(int i_s = 0 ; i_s < Config::Get()->GetParameterI("num_samps"); i_s++) {                                                                                                                 
    deconv_wf->SetAmp(i_s, hResult->GetBinContent(i_s + 1));
  }

  delete hResponse;
  delete hSource;
  delete hResult;
  delete fft_source;
  delete fft_response;
  delete fft_back;
  return deconv_wf;
}

/*--------------------------------------------------------------------*/

Waveform* Waveform::FFTFilter() {
  Waveform* filtered = new Waveform();

  TH1D* hSource = new TH1D("hSource", "hSource", Config::Get()->GetParameterI("num_samps"), 0,
                           Config::Get()->GetParameterI("num_samps") * 1e-3 / Config::Get()->GetParameterD("sampling_rate"));
  for (int i = 0; i < Config::Get()->GetParameterI("num_samps"); i++) {
    hSource->SetBinContent(i + 1, (double)_amplitudes[i]);
  }

  int N = Config::Get()->GetParameterI("num_samps");

  TH1 *hsource_m = 0;
  TVirtualFFT::SetTransform(0);
  hsource_m = hSource->FFT(hsource_m, "RE");
  TVirtualFFT *fft_source = TVirtualFFT::GetCurrentTransform();
  delete hsource_m;

  TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &N, "C2R K");
  for(int i_s = 0; i_s < Config::Get()->GetParameterI("num_samps"); i_s++) {
    Double_t re, im;
    fft_source->GetPointComplex(i_s, re, im);
    TComplex cSource(re, im);
    fft_back->SetPointComplex(i_s, cSource);
    if(i_s >  Config::Get()->GetParameterD("fft_filter_freq") * Config::Get()->GetParameterI("num_samps") / (Config::Get()->GetParameterD("sampling_rate") * 1e3) &&
       i_s < Config::Get()->GetParameterI("num_samps") -  
       Config::Get()->GetParameterD("fft_filter_freq") * Config::Get()->GetParameterI("num_samps") / (Config::Get()->GetParameterD("sampling_rate") * 1e3)) {
      TComplex filter(0,0);
      fft_back->SetPointComplex(i_s, filter);
    }
  }

  fft_back->Transform();

  TH1 *hResult = 0;
  hResult = TH1::TransformHisto(fft_back, hResult,"Re");
  hResult->Scale((double)1/ Config::Get()->GetParameterI("num_samps"));

  for(int i_s = 0 ; i_s < Config::Get()->GetParameterI("num_samps"); i_s++) {
    filtered->SetAmp(i_s, hResult->GetBinContent(i_s + 1));
  }

  delete hSource;
  delete hResult;
  delete fft_source;
  delete fft_back;
  
  return filtered;
  

}

/*--------------------------------------------------------------------*/

Waveform* Waveform::WienerDeconvolve(Waveform* response) {
  Waveform* deconv_wf = new Waveform();


  TFile *fNoise = new TFile(Config::Get()->GetParameterS("noise_samp_file").c_str());
  TGraph* gNoise = (TGraph*)fNoise->Get(Config::Get()->GetParameterS("noise_samp_name").c_str());
  delete fNoise;

  TF1 *fSignal = new TF1("fSignal","pow([0],2) / (2* TMath::Pi()) + pow([0],2) * (exp(2 * [1] / [0]) - 2 * exp(2 * [1] / [0]) * (cos([1]*x) - [0]*x*sin([1] * x))) / (2 * TMath::Pi() * pow(1+pow([0],2) * pow(x,2),2))",-3, 3);
  fSignal->SetParameters(1.4, 1);


  TH1D* hResponse = new TH1D("hResponse", "hResponse", Config::Get()->GetParameterI("num_samps"), 0,
                             Config::Get()->GetParameterI("num_samps") * 1e-3 / Config::Get()->GetParameterD("sampling_rate"));
  TH1D* hSource = new TH1D("hSource", "hSource", Config::Get()->GetParameterI("num_samps"), 0,
                           Config::Get()->GetParameterI("num_samps") * 1e-3 / Config::Get()->GetParameterD("sampling_rate"));
  TH1D* hNoise = new TH1D("hNoise", "hNoise", Config::Get()->GetParameterI("num_samps"), 0,
                           Config::Get()->GetParameterI("num_samps") * 1e-3 / Config::Get()->GetParameterD("sampling_rate"));


  //Setting up histograms 
  for (int i = 0; i < Config::Get()->GetParameterI("num_samps"); i++) {
    hResponse->SetBinContent(i + 1, response->GetAmp(i));
    hSource->SetBinContent(i + 1, (double)_amplitudes[i]);
    hNoise->SetBinContent(i + 1, gNoise->GetY()[i]);
  }

  //Calculating relevant Fourier transforms
  int N = Config::Get()->GetParameterI("num_samps");
  TH1 *hresponse_m = new TH1D("hresponse_m", "hresponse_m", Config::Get()->GetParameterI("num_samps"), 0, 
			      Config::Get()->GetParameterI("num_samps") * 1e-3 / Config::Get()->GetParameterD("sampling_rate"));
  TVirtualFFT::SetTransform(0);
  hresponse_m = hResponse->FFT(hresponse_m, "MAG");
  TVirtualFFT *fft_response = TVirtualFFT::GetCurrentTransform();
  TH1 *hsource_m = 0;
  TVirtualFFT::SetTransform(0);
  hsource_m = hSource->FFT(hsource_m, "MAG");
  TVirtualFFT *fft_source = TVirtualFFT::GetCurrentTransform();
  delete hsource_m;
  TH1 *hnoise_m = new TH1D("hnoise_m", "hnoise_m", Config::Get()->GetParameterI("num_samps"), 0,
			   Config::Get()->GetParameterI("num_samps") * 1e-3 / Config::Get()->GetParameterD("sampling_rate"));
  TVirtualFFT::SetTransform(0);
  hnoise_m = hNoise->FFT(hnoise_m, "MAG");
  TVirtualFFT *fft_noise = TVirtualFFT::GetCurrentTransform();

  TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &N, "C2R K");
  for(int i_s = 0; i_s < Config::Get()->GetParameterI("num_samps"); i_s++) {
    Double_t re, im;
    fft_response->GetPointComplex(i_s, re, im);
    TComplex cResponse(re, im);
    fft_source->GetPointComplex(i_s, re, im);
    TComplex cSource(re, im);
    TComplex one(1,0);

    //Calculating the SNR
    // cout << i_s << endl;
    double SNR;
    if(i_s <= Config::Get()->GetParameterI("num_samps")/2) SNR = fSignal->Eval(2 * TMath::Pi() * Config::Get()->GetParameterD("sampling_rate") * 1e9 * i_s 
									    / Config::Get()->GetParameterI("num_samps")) / pow(hnoise_m->GetBinContent(i_s + 1),2);
    else SNR = fSignal->Eval(2 * TMath::Pi() *  Config::Get()->GetParameterD("sampling_rate") * 1e9 * (Config::Get()->GetParameterI("num_samps") - i_s) 
			     / Config::Get()->GetParameterI("num_samps")) / pow(hnoise_m->GetBinContent(i_s + 1),2);
    
    TComplex g = (one / cResponse) * (1/ (1 + 1/(pow(hresponse_m->GetBinContent(i_s + 1),2) * SNR)));
    TComplex x = g * cSource;
    fft_back->SetPointComplex(i_s, x);
  }

  fft_back->Transform();

  TH1 *hResult = 0;
  hResult = TH1::TransformHisto(fft_back, hResult,"Re");
  hResult->Scale((double)1/ Config::Get()->GetParameterI("num_samps"));
  
  for(int i_s = 0 ; i_s < Config::Get()->GetParameterI("num_samps"); i_s++) {
    deconv_wf->SetAmp(i_s, hResult->GetBinContent(i_s + 1));
  }

  delete hResponse;
  delete hresponse_m;  
  delete hSource;
  delete hNoise;
  delete hnoise_m; 
  delete gNoise;
  delete fSignal;
  delete hResult;
  delete fft_source;
  delete fft_response;
  delete fft_noise;
  delete fft_back;
  return deconv_wf;

}

/*--------------------------------------------------------------------*/

void Waveform::WaveformFree() {
  _times.clear(); vector<float>().swap(_times);
  _samples.clear(); vector<int>().swap(_samples);
  _amplitudes.clear(); vector<float>().swap(_amplitudes);
  _baseline.clear(); vector<float>().swap(_baseline);
  _baseline_stdev.clear(); vector<float>().swap(_baseline_stdev);
  _baseline_range.clear();vector<float>().swap(_baseline_range); 
}
