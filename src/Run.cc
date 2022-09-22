#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include "Run.h"
#include "Event.h"
#include "TGraph.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TAxis.h"
#include "TTree.h"
#include "TH1F.h"
#include "TVirtualFFT.h"
#include "TSpectrum.h"
#include "Output.h"
#include "TLine.h"
#include "TRandom.h"
#include "TMath.h"
#include "Waveform.h"

using namespace std;

/*--------------------------------------------------------------------*/

int Run::ProcessSiPM() {
  _ifiles[0]->seekg(0, _ifiles[0]->end);
  long int rawFileLength = _ifiles[0]->tellg();
  _ifiles[0]->seekg(0, _ifiles[0]->beg);
  int tot_events = rawFileLength / (Config::Get()->GetParameterI("num_samps")*2);
  int Nevents = Config::Get()->GetParameterI("nevents");
  if ((Nevents == -1) || (Nevents > tot_events)) {Nevents = tot_events;}
  Config::Get()->SetParameter("nevents",Nevents);
  
  // Setting up output
  Output::Get()->SetupSiPM();  

  for (int i_ev =  0; i_ev < Nevents; i_ev++){
    if(i_ev % 500 == 499) cout << "event" << i_ev << endl;
    Event* anEvent = new Event(i_ev);
    for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){
      Waveform* wf = ReadBinary(i_ch, false);
      anEvent->SetRawWF(i_ch, wf);
    }
    anEvent->ProcessEvent();
    _events.push_back(anEvent);

    // When the event number is a multiple of avg_frequency, recalculate the average pulse and deconvolve 
    if ((i_ev + 1) % Config::Get()->GetParameterI("avg_frequency") == 0) {
      CalcAvgPulse();
      CalcAvgWF();
      
      if(Config::Get()->GetParameterS("deconvolve") == "y") {
	Waveform* response = new Waveform();
	if(Config::Get()->GetParameterS("calc_avg_pulse") == "y") {
	  int deconv_channel = Config::Get()->GetParameterI("avg_channel");
	  for(int i_s = 0 ; i_s < Config::Get()->GetParameterI("pulse_avg_samps"); i_s++) {
	      response->SetAmp(i_s, _avgPulse[deconv_channel][Config::Get()->GetParameterI("pulse_avg_buf") + i_s] / Output::Get()->total_SPE_pulses[deconv_channel]);
	  }
	}
	else {
	  // Reading response funcion from file
	  TFile *f = new TFile(Config::Get()->GetParameterS("avgPulse_file").c_str());
	  TGraph* avgPulse_input = (TGraph*)f->Get(Config::Get()->GetParameterS("avgPulse_name").c_str());
	  for(int i_s = 0 ; i_s < Config::Get()->GetParameterI("pulse_avg_samps"); i_s++) {
	      response->SetAmp(i_s, avgPulse_input->GetY()[Config::Get()->GetParameterI("pulse_avg_buf") + i_s]);
	  }
	  delete f;
	  delete avgPulse_input;
	}
	cout << "Deconvolving.." << endl;
	for (unsigned int i_d =  0; i_d < _events.size(); i_d++){
	  for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){
	    if(_events[i_d]->GetChannel(i_ch)->GetNPulses() > 0) {
	      Output::Get()->eventID = i_ev - Config::Get()->GetParameterI("avg_frequency") + i_d + 1;
	      _events[i_d]->GetChannel(i_ch)->PulseDeconv(response);
	    }
	  }
	}
	response->WaveformFree();
	delete response;
      }
      //Writing to tree
      for (unsigned int i_d =  0; i_d < _events.size(); i_d++){
	_events[i_d]->PulseOutput();
      }

      // Freeing events
      if(i_ev < Nevents - 1) {
	for (unsigned int i_d =  0; i_d < _events.size(); i_d++){ 
	  _events[i_d]->EventFree();
	  delete _events[i_d];
	}
	_events.clear();
      }
    }
  }

  // Calculating average pulses
   for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){
     for(int i = 0; i < Config::Get()->GetParameterI("pulse_avg_samps") + Config::Get()->GetParameterI("pulse_avg_buf"); i++) {
       _avgPulseAll[i_ch][i] /= Output::Get()->total_pulses[i_ch]; 
       _avgPulse[i_ch][i] /=  Output::Get()->total_SPE_pulses[i_ch]; 
     }
   }

  //Create Graphs of Average Waveforms and Pulses
  for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){
    Output::Get()->gAvgWF[i_ch] = new TGraph(Config::Get()->GetParameterI("num_samps"),&_events[0]->GetChannel(i_ch)->GetRawWF()->GetTimes()[0], &_avgWF[i_ch][0]); 
    Output::Get()->gAvgWF[i_ch]->SetName(Form("gAvgWF_ch%d", Config::Get()->GetParameterIvec("ch_ids")[i_ch]));
    Output::Get()->gAvgPulse[i_ch] = new TGraph(Config::Get()->GetParameterI("pulse_avg_samps") + Config::Get()->GetParameterI("pulse_avg_buf"),
					  &_events[0]->GetChannel(i_ch)->GetRawWF()->GetTimes()[0], &_avgPulseAll[i_ch][0]); 
    Output::Get()->gAvgPulse[i_ch]->SetName(Form("gAvgPulse_ch%d", Config::Get()->GetParameterIvec("ch_ids")[i_ch]));
    Output::Get()->gAvgSPEPulse[i_ch] = new TGraph(Config::Get()->GetParameterI("pulse_avg_samps") + Config::Get()->GetParameterI("pulse_avg_buf"),
					      &_events[0]->GetChannel(i_ch)->GetRawWF()->GetTimes()[0],&_avgPulse[i_ch][0]);  
    Output::Get()->gAvgSPEPulse[i_ch]->SetName(Form("gAvgSPEPulse_ch%d", Config::Get()->GetParameterIvec("ch_ids")[i_ch]));
  }
  Output::Get()->WriteSiPM();
  Output::Get()->_ofile->Close();
  return 0;
}

/*--------------------------------------------------------------------*/

void Run::ProcessTPC() {
  _ifiles[0]->seekg(0, _ifiles[0]->end);
  long int rawFileLength = _ifiles[0]->tellg();
  _ifiles[0]->seekg(0, _ifiles[0]->beg);
  int tot_events = rawFileLength / (Config::Get()->GetParameterI("num_samps")*2);
  int Nevents = Config::Get()->GetParameterI("nevents");
  if ((Nevents == -1) || (Nevents > tot_events)) {Nevents = tot_events;}
  Config::Get()->SetParameter("nevents",Nevents);

  Output::Get()->SetupTPC();          

  for (int i_ev =  0; i_ev < Nevents; i_ev++){
    if(i_ev % 500 == 0) cout << "event" << i_ev << endl;
    Event* anEvent = new Event(i_ev);
    for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){
      Waveform* wf = ReadBinary(i_ch, false);
      anEvent->SetRawWF(i_ch, wf);
    }
    anEvent->ProcessEvent();
    anEvent->ClusterOutput();
    _events.push_back(anEvent);
    CalcAvgS1();
    anEvent->EventFree();
    delete anEvent;
    _events.clear();
  }

  for(int i = 0; i < Config::Get()->GetParameterI("S1_avg_samps") + Config::Get()->GetParameterI("S1_avg_buf"); i++) {
    _avgS1[i] /= Output::Get()->total_S1_pulses;
  }

  Waveform* temp_time = new Waveform();
  Output::Get()->gAvgS1 = new TGraph(Config::Get()->GetParameterI("S1_avg_samps") + Config::Get()->GetParameterI("S1_avg_buf"),&temp_time->GetTimes()[0], &_avgS1[0]);
  Output::Get()->gAvgS1->SetName("gAvgS1");

 
  Output::Get()->WriteTPC();
  Output::Get()->_ofile->Close();
}

/*--------------------------------------------------------------------*/

void Run::ProcessNoise() {
  _ifiles[0]->seekg(0, _ifiles[0]->end);
  long int rawFileLength = _ifiles[0]->tellg();
  _ifiles[0]->seekg(0, _ifiles[0]->beg);
  int tot_events = rawFileLength / (Config::Get()->GetParameterI("num_samps")*2);
  int Nevents = Config::Get()->GetParameterI("nevents");
  if ((Nevents == -1) || (Nevents > tot_events)) {Nevents = tot_events;}
  Config::Get()->SetParameter("nevents",Nevents);

  Output::Get()->SetupNoise();

  for (int i_ev =  0; i_ev < Nevents; i_ev++){
    if(i_ev % 500 == 0) cout << "event" << i_ev << endl;
    Event* anEvent = new Event(i_ev);
    for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){
      Waveform* wf = ReadBinary(i_ch, false);
      anEvent->SetRawWF(i_ch, wf);
    }
    anEvent->ProcessEventNoise();
    //Build average FFT
    for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){
      for (int bin = 1; bin <= Config::Get()->GetParameterI("num_samps")/2; bin++) {
	double old_bin_content =  Output::Get()->hAvgFFT[i_ch]->GetBinContent(bin);  
	Output::Get()->hAvgFFT[i_ch]->SetBinContent(bin, old_bin_content + anEvent->GetChannel(i_ch)->GetBlSubWF()->GetFFT(bin) / Nevents);
      }
    }
  }
   Output::Get()->WriteNoise();
}

/*--------------------------------------------------------------------*/


Waveform* Run::ReadBinary(int channel, bool header){
  int Ns = Config::Get()->GetParameterI("num_samps");
  char* eventBinary = new char[Ns*2];
  short int* waveform = new short int[Ns];
  Waveform* wf = new Waveform();
  
  _ifiles[channel]->read((char*) waveform, Ns*2);
  
  for(int is=0; is<Ns; is++){  
    wf->SetAmp(is, waveform[is]);
  }

  delete[] eventBinary;
  delete[] waveform;

  return wf;
}

/*--------------------------------------------------------------------*/


Waveform* Run::ReadBinary2(int channel, bool header){
  int Ns = Config::Get()->GetParameterI("num_samps");
  char* eventBinary = new char[Ns*2];

  Waveform* wf = new Waveform();

  _ifiles[channel]->read(eventBinary, Ns*2+header*24);

  for(int is=0; is<Ns; is++){
    wf->SetAmp(is,(((unsigned short)((eventBinary[is*2+1 + header*24]&0x7F) + (eventBinary[is*2+1 + header*24]<0 ? 128 : 0)))<<8)
               + ((eventBinary[is*2 + header*24]&0x7F) + (eventBinary[is*2 + header*24]<0 ? 128 : 0) )) ;
  }

  return wf;
}


/*--------------------------------------------------------------------*/


void Run::CalcAvgWF(){
  Waveform*  wf;
  for (unsigned int i_ev =  0; i_ev < _events.size(); i_ev++){
    for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){
      wf =_events[i_ev]->GetChannel(i_ch)->GetBlSubWF();
      for(int i_s = 0 ; i_s < Config::Get()->GetParameterI("num_samps"); i_s++){
	_avgWF[i_ch][i_s] += wf->GetAmp(i_s) / Config::Get()->GetParameterI("nevents");
      }
    }
  }
}


/*--------------------------------------------------------------------*/ 

void Run::CalcAvgPulse() {
  Waveform* wf;
  vector<Pulse*> pulses;
  int total_pulses = 0;
  for (unsigned int i_ev =  0; i_ev < _events.size(); i_ev++){
    for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){
      wf =_events[i_ev]->GetChannel(i_ch)->GetBlSubWF();
      pulses = _events[i_ev]->GetChannel(i_ch)->GetPulses();
      for(unsigned int i_p = 0; i_p < pulses.size(); i_p++) {
	total_pulses++;
	int start_samp = pulses[i_p]->GetStartSamp() - Config::Get()->GetParameterI("pulse_avg_buf");
	for(int i = 0; i < Config::Get()->GetParameterI("pulse_avg_samps") + Config::Get()->GetParameterI("pulse_avg_buf"); i++, start_samp++) {
	  _avgPulseAll[i_ch][i] += wf->GetAmp(start_samp); 
	  if(pulses[i_p]->GetNAfterPulses() == 0 && 
	     pulses[i_p]->GetIntegral() > Config::Get()->GetParameterD("min_SPE_charge") && 
	     pulses[i_p]->GetIntegral() < Config::Get()->GetParameterD("max_SPE_charge") &&
	     pulses[i_p]->GetHeight() > Config::Get()->GetParameterD("min_SPE_height") &&
	     pulses[i_p]->GetHeight() < Config::Get()->GetParameterD("max_SPE_height") &&
	     pulses[i_p]->GetEndTime() - pulses[i_p]->GetStartTime() > Config::Get()->GetParameterD("min_SPE_length") &&
	     pulses[i_p]->GetEndTime() - pulses[i_p]->GetStartTime() < Config::Get()->GetParameterD("max_SPE_length")) { 
	    _avgPulse[i_ch][i] += wf->GetAmp(start_samp);
	  }
	}
      }
    }
  }
}


/*--------------------------------------------------------------------*/

void Run::CalcAvgS1() {
  Waveform* wf;
  vector<Cluster*> clusters;
  for (unsigned int i_ev =  0; i_ev < _events.size(); i_ev++){
    wf = _events[i_ev]->GetSumWF();
    clusters = _events[i_ev]->GetClusters();
    for(unsigned int i_c = 0; i_c < clusters.size(); i_c++) {
      int start_samp = clusters[i_c]->GetStartSamp() -  Config::Get()->GetParameterI("S1_avg_buf");
      for(int i = 0; i < Config::Get()->GetParameterI("S1_avg_samps") + Config::Get()->GetParameterI("S1_avg_buf"); i++, start_samp++) {
	if(clusters[i_c]->GetS1()) _avgS1[i] += wf->GetAmp(start_samp);
      }
    }
  }
}


/*--------------------------------------------------------------------*/

void Run::EventViewer() {
  
  int prev_evID = 0;
  string input;
  int input_evID = 0;
  int skip_ID = 0;
  bool header = false; 
  double min_range_x = 0; double max_range_x = 0;
  double min_range_y = 0; double max_range_y = 0;

  TCanvas **c= new TCanvas*[Config::Get()->GetParameterI("num_chans")];
  for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){
    c[i_ch] = new TCanvas(Form("c%d",i_ch),"Event 0",1400,800);
  }


  Event* anEvent = new Event(0);
  cout << "Displaying the first event..." << endl;
  
  //Set up Event
  for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){
    Waveform* wf = ReadBinary(i_ch, false);
    anEvent->SetRawWF(i_ch, wf);
  }
  anEvent->ProcessEvent();

  // Plot Event
  cout << "EVENT ID: " << anEvent->GetEventID() << endl;
  cout << "====================================================================" << endl;
  cout << "====================================================================" << endl;
  for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){
    c[i_ch]->cd();
    anEvent->EventPrint(i_ch);
    anEvent->DisplayEvent(i_ch);
    c[i_ch]->Modified();
    c[i_ch]->Update();
    gSystem->ProcessEvents();
  }

  while(!_ifiles[0]->eof()) {
    
      cout << "====================================================================" << endl;
      cout << "Type \"zoom <min_xvalue> <max_xvalue> <min_yvalue> <max_yvalue>\" to zoom in on the waveform  OR " << endl;
      cout << "Type \"print\" to save waveform to file OR " << endl;
      cout << "Press enter to see next event or input an event ID to skip forward: ";
      
      getline(cin,input);

      // Validate input
      while(true) {
	bool digit = true;
	if (input == "" || input == "print" || input == "q") break;
	if (input.substr(0,4) == "zoom") {
	  istringstream ss(input);
	  ss.ignore(5);
	  string min_x, max_x, min_y, max_y;
	  getline(ss,min_x,' ');
	  for (unsigned int i = 0; i < min_x.length(); i++) {if (isdigit(min_x[i]) == false && min_x[i] != '.') digit = false;}
	  getline(ss,max_x,' ');
	  for (unsigned int i = 0; i < max_x.length(); i++) {if (isdigit(max_x[i]) == false && max_x[i] != '.') digit = false;}
	  getline(ss,min_y,' ');
	  for (unsigned int i = 0; i < min_y.length(); i++) {if (isdigit(min_y[i]) == false && min_y[i] != '.' && min_y[0] != '-') digit = false;}
	  getline(ss,max_y,' ');
	  for (unsigned int i = 0; i < max_y.length(); i++) {if (isdigit(max_y[i]) == false && max_y[i] != '.' && min_y[0] != '-') digit = false;}
	  if (digit && (atof(max_x.c_str()) > atof(min_x.c_str())) && (atof(max_y.c_str()) >= atof(min_y.c_str()))){
	    min_range_x = atof(min_x.c_str());
	    max_range_x = atof(max_x.c_str());
	    min_range_y = atof(min_y.c_str());
            max_range_y = atof(max_y.c_str());
	    break;
	  } 
	}
	for (unsigned int i = 0; i < input.length(); i++) {
	  if (isdigit(input[i]) == false) digit = false;
        }
	if (digit && (atoi(input.c_str()) > prev_evID)) {
	    input_evID = atoi(input.c_str());
	    break;
	}
	cout << "Invalid input. Please try again: ";
	getline(cin, input);
      }

      if (input == "q") { 
	for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){
	  c[i_ch]->Close();
	} 
	c[0]->Connect("TCanvas", "Closed()", "TApplication", gApplication, "Terminate()");
	return;
      }
      if (input == "print") {
        for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){
          TFile * fout  = new TFile(Form("output/eventView%d.root", i_ch),"UPDATE");
	  c[i_ch]->Write();
          fout->Close();
        }
        continue;
      } // if print write to file
      else if (input.substr(0,4) == "zoom") {
	for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){
	  c[i_ch]->Clear(); 
	  c[i_ch]->cd();
	  anEvent->DisplayEvent(i_ch,min_range_x,max_range_x, min_range_y, max_range_y);
	  c[i_ch]->Modified();
	  c[i_ch]->Update();
	  gSystem->ProcessEvents();
	}
	continue;
      } // if zoom, modify canvas
      else if (input == "") {
	prev_evID++;
	skip_ID = 0;
	input_evID = prev_evID;
      } // if enter advance to next event
      else {
	skip_ID = input_evID - prev_evID - 1;
	prev_evID = input_evID;	
      } // if valid integer skip ahead
      
      anEvent->SetEventID(input_evID);
      // for each channel, skip ahead in the files corresponding to skipID, then display event
      for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){
	_ifiles[i_ch]->ignore(skip_ID*Config::Get()->GetParameterI("num_samps")*2+header*24);
        Waveform* wf = ReadBinary(i_ch, false);
	anEvent->GetChannel(i_ch)->Reset();
	anEvent->SetRawWF(i_ch, wf);
      }
      anEvent->ProcessEvent();
      
      cout << "EVENT ID: " << anEvent->GetEventID() << endl;
      cout << "====================================================================" << endl;
      cout << "====================================================================" << endl;
      for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){
	c[i_ch]->Clear();
        c[i_ch]->cd();
	c[i_ch]->SetTitle(Form("Event %d", anEvent->GetEventID()));
	anEvent->EventPrint(i_ch);
	anEvent->DisplayEvent(i_ch);
	c[i_ch]->Modified();
	c[i_ch]->Update();
	gSystem->ProcessEvents();
      }    
  }
  
}

/*--------------------------------------------------------------------*/


void Run::AddTemplatePulses(TGraph* avg, TRandom* r1, double* wf, int npulses, double tau, int offset, int sig_num) {
  double start_time = 0;
  int start_samp = 0;
  int num_samps = Config::Get()->GetParameterI("waveform_length") * 1e3 * Config::Get()->GetParameterD("sampling_rate");

  vector<double> pulses, events, afterpulses;
  vector<double> start_times, integrals, heights;


  for(int i = 0; i < npulses; i++) {
    Waveform* avg_template = new Waveform();
    for(int i_s = 0 ; i_s < Config::Get()->GetParameterI("num_samps"); i_s++) {
      if(i_s < Config::Get()->GetParameterI("pulse_avg_samps")) {
	avg_template->SetAmp(i_s, avg->GetY()[Config::Get()->GetParameterI("pulse_avg_buf") + i_s]);
      }
    }
    // Randomly determining the multiplicative factor
    double charge_diff = r1->Gaus(0, 25.54);
    double charge_ratio = (367 + charge_diff) / 367;
    
    for(int i_s = 0 ; i_s < Config::Get()->GetParameterI("num_samps"); i_s++) {
      avg_template->SetAmp(i_s, avg_template->GetAmp(i_s) * charge_ratio);
    }
    avg_template->CalcExtrema();

    start_time = r1->Exp(tau);
    start_samp = offset + start_time * 1e3 * Config::Get()->GetParameterD("sampling_rate");
    for(int j = start_samp; j < start_samp + Config::Get()->GetParameterI("pulse_avg_samps") && j < num_samps; j++) {
      wf[j] += avg_template->GetAmp(j - start_samp);
    }

    Output::Get()->rPulseID[sig_num].push_back(-1);
    Output::Get()->rEventID[sig_num].push_back(-1);
    Output::Get()->rStart_time[sig_num].push_back(start_samp * 1e-3 / Config::Get()->GetParameterD("sampling_rate"));
    Output::Get()->rIntegral[sig_num].push_back(avg_template->Integrate(0, Config::Get()->GetParameterI("pulse_avg_samps")));
    Output::Get()->rHeight[sig_num].push_back(avg_template->GetMax());
    Output::Get()->rAfterPulses[sig_num].push_back(0);
    
    avg_template->WaveformFree();
    delete avg_template;
  }
}


/*--------------------------------------------------------------------*/


void Run::AddRandPulses(TRandom* r1, double* wf, int npulses, double tau, int offset, int sig_num) {
  int event_number = 0, pulse_number = 0, start_samp = 0;
  double start_time = 0;
  Pulse* rand_pulse;
  int num_samps = Config::Get()->GetParameterI("waveform_length") * 1e3 * Config::Get()->GetParameterD("sampling_rate");
  int source_chID = Config::Get()->GetParameterI("source_channel");

  vector<double> pulses, events, afterpulses;
  vector<double> start_times, integrals, heights;

  for(int i = 0; i < npulses; i++) {
    event_number = 0; pulse_number = 0;
    while (pulse_number == 0) {
      event_number = (int)r1->Uniform(0,Config::Get()->GetParameterI("nevents"));
      pulse_number = (int)r1->Uniform(1, _events[event_number]->GetChannel(source_chID)->GetNPulses());
      if(pulse_number == 0) continue;
      if( _events[event_number]->GetChannel(source_chID)->GetPulse(pulse_number-1)->GetNAfterPulses() == 0 &&
    	  _events[event_number]->GetChannel(source_chID)->GetPulse(pulse_number-1)->GetIntegral() > Config::Get()->GetParameterD("min_SPE_charge") &&
    	  _events[event_number]->GetChannel(source_chID)->GetPulse(pulse_number-1)->GetIntegral() < Config::Get()->GetParameterD("max_SPE_charge") &&
    	  _events[event_number]->GetChannel(source_chID)->GetPulse(pulse_number-1)->GetHeight() > Config::Get()->GetParameterD("min_SPE_height") &&
    	  _events[event_number]->GetChannel(source_chID)->GetPulse(pulse_number-1)->GetHeight() < Config::Get()->GetParameterD("max_SPE_height") &&
    	  _events[event_number]->GetChannel(source_chID)->GetPulse(pulse_number-1)->GetEndTime() -
    	  _events[event_number]->GetChannel(source_chID)->GetPulse(pulse_number-1)->GetStartTime() > Config::Get()->GetParameterD("min_SPE_length") &&
    	  _events[event_number]->GetChannel(source_chID)->GetPulse(pulse_number-1)->GetEndTime() -
    	  _events[event_number]->GetChannel(source_chID)->GetPulse(pulse_number-1)->GetStartTime() < Config::Get()->GetParameterD("max_SPE_length")) break;
      else pulse_number = 0;
    }
    
    rand_pulse = _events[event_number]->GetChannel(source_chID)->GetPulse(pulse_number-1);
    start_time = r1->Exp(tau);
    start_samp = offset + start_time * 1e3 * Config::Get()->GetParameterD("sampling_rate");
    for(int j = start_samp; j < start_samp + (rand_pulse->GetEndSamp() - rand_pulse->GetStartSamp()) && j < num_samps; j++) {
      wf[j] += _events[event_number]->GetChannel(source_chID)->GetBlSubWF()->GetAmp(rand_pulse->GetStartSamp() + j - start_samp);
    }

    // Filling output
    Output::Get()->rPulseID[sig_num].push_back(pulse_number-1);
    Output::Get()->rEventID[sig_num].push_back(event_number);
    Output::Get()->rStart_time[sig_num].push_back(start_samp * 1e-3 / Config::Get()->GetParameterD("sampling_rate"));
    Output::Get()->rIntegral[sig_num].push_back(_events[event_number]->GetChannel(source_chID)->GetPulse(pulse_number-1)->GetIntegral());
    Output::Get()->rHeight[sig_num].push_back(_events[event_number]->GetChannel(source_chID)->GetPulse(pulse_number-1)->GetHeight());
    Output::Get()->rAfterPulses[sig_num].push_back(_events[event_number]->GetChannel(source_chID)->GetPulse(pulse_number-1)->GetNAfterPulses());

  }
}


/*--------------------------------------------------------------------*/

void Run::GenerateWaveforms() {
  _ifiles[0]->seekg(0, _ifiles[0]->end);
  long int rawFileLength = _ifiles[0]->tellg();
  _ifiles[0]->seekg(0, _ifiles[0]->beg);
  int tot_events = rawFileLength / (Config::Get()->GetParameterI("num_samps")*2);
  int Nevents = Config::Get()->GetParameterI("nevents");
  if ((Nevents == -1) || (Nevents > tot_events)) {Nevents = tot_events;}
  Config::Get()->SetParameter("nevents",Nevents);

  Output::Get()->SetupSiPM();
  Output::Get()->SetupGenerateWF();

  int num_samps = Config::Get()->GetParameterI("waveform_length") * 1e3 * Config::Get()->GetParameterD("sampling_rate");
  int buffer = 250; //time before first pulse in samples
  int source_chID = Config::Get()->GetParameterI("source_channel");


  //Processing all events
  cout << "Processing events..." << endl;
  for (int i_ev =  0; i_ev < Nevents; i_ev++){
    if(i_ev % 500 == 499) cout << "event" << i_ev << endl;
    Event* anEvent = new Event(i_ev);
    for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){
      Waveform* wf = ReadBinary(i_ch, false);
      anEvent->SetRawWF(i_ch, wf);
    }
    anEvent->ProcessEvent();
    _events.push_back(anEvent);
  }


  //Find S2 amplification information
  TFile *fS2 = new TFile(Config::Get()->GetParameterS("s2_amp_file").c_str());
  TGraph* gS2 = (TGraph*)fS2->Get(Config::Get()->GetParameterS("s2_amp_name").c_str());
  delete fS2;
  

  TRandom *r1 = new TRandom(0);

  //Calculating PEs for each scatter event
  vector<vector<int> > nPE_s1;
  vector<vector<int> > nPE_s2;
  vector<vector<double> > s1_s2_dt;
  vector<double> fprompt;


  for(int i_s2 = 0; i_s2 < Config::Get()->GetParameterI("num_S2"); i_s2++) {
    vector<int> s2_vecI;
    vector<double> s2_vecD;
    nPE_s1.push_back(s2_vecI);
    nPE_s2.push_back(s2_vecI);
    s1_s2_dt.push_back(s2_vecD);

    double ep = 0.0135 * Config::Get()->GetParameterDvec("energy")[i_s2];
    double sn = log(1 + 1.1383 * 0.953 * ep) / (2 * (0.953 * ep + 0.01321 * pow(0.953 * ep, 0.21226) + 0.19593 * pow(0.953 * ep, 0.5)));
    double se = 0.145 * sqrt(ep); 
    double n_electrons = 6.8 * 1e3 * ep * se / (se + sn);
  
    if(Config::Get()->GetParameterS("recoil") == "ER") fprompt.push_back(0.205563 * exp(-Config::Get()->GetParameterDvec("energy")[i_s2] / 7.16924) + 0.279701);
    else {
      double NR_energy[] = {10.3, 14.8, 16.9, 20.5, 25.4, 28.7, 36.1, 49.7, 57.3};
      double NR_fprompt[] = {0.527, 0.567, 0.582, 0.606, 0.641, 0.643, 0.664, 0.704, 0.715};
      TGraph* gFpromptNR = new TGraph(9, NR_energy, NR_fprompt);
      fprompt.push_back(gFpromptNR->Eval(Config::Get()->GetParameterDvec("energy")[i_s2], 0, "S"));
    }


    // For each event, randomly determining the S1 and S2 response
    for (int n = 0; n < Config::Get()->GetParameterI("generate_events"); n++) {
      int nS1, nS2;
      double dt;
      double radius;
      if(Config::Get()->GetParameterS("recoil") == "ER") {
	nS2 = (200/9.2 + 0.11 * pow(Config::Get()->GetParameterDvec("energy")[i_s2],1.71)) * log(1 + (9.2/200) * 54.4 * Config::Get()->GetParameterDvec("energy")[i_s2]);
	nS1 = (int)(Config::Get()->GetParameterDvec("energy")[i_s2] * 1e3 / Config::Get()->GetParameterD("work")) - nS2;
	nS1 = r1->Binomial(nS1, Config::Get()->GetParameterD("s1_detection_eff"));   
	radius = r1->Uniform(1,Config::Get()->GetParameterD("max_radius"));
	nS2 = r1->Poisson(gS2->Eval(radius, 0, "") * nS2);
      }
      else if(Config::Get()->GetParameterS("recoil") == "NR") {
	nS2 = 6.8 * 1e3 * ep * se / (se + sn);
	nS2 = 200/(8.1 * n_electrons) * log(1 + 8.1 * n_electrons/200) * n_electrons; 
	nS1 = (int)(Config::Get()->GetParameterDvec("energy")[i_s2] * 1e3 / Config::Get()->GetParameterD("work")) - nS2;
	nS1 = r1->Binomial(nS1, Config::Get()->GetParameterD("s1_detection_eff"));
	radius = r1->Uniform(1,Config::Get()->GetParameterD("max_radius"));
	nS2 = r1->Poisson(gS2->Eval(radius, 0, "") * nS2);
      }

      nPE_s1[i_s2].push_back(nS1);
      nPE_s2[i_s2].push_back(nS2);
      dt = r1->Uniform(0, Config::Get()->GetParameterD("max_dz")) * 10 / Config::Get()->GetParameterD("drift_velocity");
      s1_s2_dt[i_s2].push_back(dt);
    }
  }

  // Setting channel probability distribution
  vector<double> channel_prob;
  if(Config::Get()->GetParameterS("uniform_channels") == "y") {
    for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("generate_channels"); i_ch++) {
      channel_prob.push_back((double)1/Config::Get()->GetParameterI("generate_channels"));
    }
  }
  else {
    double sum = 0;
    for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("generate_channels"); i_ch++) {
      double temp = r1->Uniform(1,10);
      channel_prob.push_back(temp);
      sum += temp;
    }
    for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("generate_channels"); i_ch++) {
      channel_prob[i_ch] /= sum;
    }
  }

  TFile *f = new TFile(Config::Get()->GetParameterS("avgPulse_file").c_str());
  TGraph* avgPulse_input = (TGraph*)f->Get(Config::Get()->GetParameterS("avgPulse_name").c_str());
  delete f;

  string output_dir = Config::Get()->GetParameterS("generate_output_dir") + "/wave%d.dat";

  //Repeating for a give number of channels
  for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("generate_channels"); i_ch++) {
    cout << "Generating Channel " << i_ch << "..." << endl;

    ofstream writef(Form(output_dir.c_str(), i_ch), ios::out | ios::binary);
    
    // Repeating for each event
    for (int n = 0; n < Config::Get()->GetParameterI("generate_events"); n++) { 
      double* wf = new double[num_samps];
      
      //Initializing the array that will hold the waveform
      for(int k = 0; k < num_samps; k++) {
	wf[k] = 0;
      }

      // Dividing PEs per channel
      vector<int> nS1ScatterPulses;
      int nS1TotalPulses = 0;
      for(int i_s2 = 0; i_s2 < Config::Get()->GetParameterI("num_S2"); i_s2++) {
	nS1ScatterPulses.push_back(r1->Binomial(nPE_s1[i_s2][n], channel_prob[i_ch]));
	nS1TotalPulses += nS1ScatterPulses[i_s2];
      }

      // Finding the number of PEs in the fast component of the S1 pulse
      int nS1FastPulses = 0;
      for(int i_s2 = 0; i_s2 < Config::Get()->GetParameterI("num_S2"); i_s2++) {
	if(Config::Get()->GetParameterS("recoil") == "ER") nS1FastPulses += r1->Binomial(nS1ScatterPulses[i_s2], fprompt[i_s2]);
	else nS1FastPulses += r1->Binomial(nS1ScatterPulses[i_s2], fprompt[i_s2]);
      }

      vector<double> noise;
	
      // Setting up output
      Output::Get()->rPulseID[0].clear();
      Output::Get()->rEventID[0].clear();
      Output::Get()->rStart_time[0].clear();
      Output::Get()->rIntegral[0].clear();
      Output::Get()->rHeight[0].clear();
      Output::Get()->rAfterPulses[0].clear();
      Output::Get()->nPE.clear();
      Output::Get()->rTime.clear();
      Output::Get()->nPE.push_back(nS1TotalPulses);
      Output::Get()->rTime.push_back(buffer * 1e-3 / Config::Get()->GetParameterD("sampling_rate"));
      

      // Adding pulses to construct an S1 signal
      if(Config::Get()->GetParameterS("use_templates") == "y") {
	AddTemplatePulses(avgPulse_input, r1, wf, nS1FastPulses, 0.007, buffer,0);
	AddTemplatePulses(avgPulse_input, r1, wf, nS1TotalPulses - nS1FastPulses, 1.4, buffer,0);
      }
      else {
	AddRandPulses(r1, wf, nS1FastPulses, 0.007, buffer,0); 
	AddRandPulses(r1, wf, nS1TotalPulses - nS1FastPulses, 1.4, buffer,0);
      }
  
      
      // Repeating for each S2 pulse
      for(int i_s2 = 0; i_s2 < Config::Get()->GetParameterI("num_S2"); i_s2++) { 

	// Dividing pulses per channel
	int nS2TotalPulses = r1->Binomial(nPE_s2[i_s2][n], channel_prob[i_ch]);

	Output::Get()->rPulseID[i_s2+1].clear();
	Output::Get()->rEventID[i_s2+1].clear();
	Output::Get()->rStart_time[i_s2+1].clear();
	Output::Get()->rIntegral[i_s2+1].clear();
	Output::Get()->rHeight[i_s2+1].clear();
	Output::Get()->rAfterPulses[i_s2+1].clear();
	Output::Get()->nPE.push_back(nS2TotalPulses);
	Output::Get()->rTime.push_back(buffer * 1e-3 / Config::Get()->GetParameterD("sampling_rate") + s1_s2_dt[i_s2][n]);

	// Adding pulses to construct an S2 signal
	if(Config::Get()->GetParameterS("use_templates") == "y") 
	  AddTemplatePulses(avgPulse_input, r1, wf, nS2TotalPulses, 3.8, buffer + s1_s2_dt[i_s2][n] * 1e3 * Config::Get()->GetParameterD("sampling_rate"), i_s2+1);
	else AddRandPulses(r1, wf, nS2TotalPulses, 3.8, buffer + s1_s2_dt[i_s2][n] * 1e3 * Config::Get()->GetParameterD("sampling_rate"), i_s2+1);

      }

      Output::Get()->_trueTree[i_ch]->Fill();
     
      //Generating noise array
      while ((int)noise.size() < num_samps) {
	int event_number = 0; int pulse_number = 0;
	while (pulse_number == 0) {
	  event_number = (int)r1->Uniform(0,Nevents);
	  pulse_number = (int)r1->Uniform(1, _events[event_number]->GetChannel(source_chID)->GetNPulses());
	}
	int nNoise = 0; int pulse_count = 0;                                                   
	int pulse_start = _events[event_number]->GetChannel(source_chID)->GetPulse(pulse_count)->GetStartSamp();
	int pulse_end = _events[event_number]->GetChannel(source_chID)->GetPulse(pulse_count)->GetEndSamp();
	for(int k = 0; k < Config::Get()->GetParameterI("num_samps") && nNoise < 500 ; k++) {                                                                                 
	  if(k >= pulse_start) {
	    k = pulse_end;
	    if(pulse_count + 1 == _events[event_number]->GetChannel(source_chID)->GetNPulses()) {
	      pulse_start = Config::Get()->GetParameterI("num_samps");
	      pulse_end = Config::Get()->GetParameterI("num_samps");                                                                                                                           
	    }                                                                                                                              
	    else {
	      pulse_count++;                         
	      pulse_start = _events[event_number]->GetChannel(source_chID)->GetPulse(pulse_count)->GetStartSamp();
	      pulse_end = _events[event_number]->GetChannel(source_chID)->GetPulse(pulse_count)->GetEndSamp();
	    }                                                                                             
	  }                                                                                                                                                                                      
	  noise.push_back(_events[event_number]->GetChannel(source_chID)->GetBlSubWF()->GetAmp(k));
	  nNoise++;
	}
      }                                                                                                                                                                                 
              
      // Writing to file
      vector<short int> amplitudes;
      for(int k = 0; k < num_samps; k++) { 
	amplitudes.push_back((short int)(wf[k] + noise[k] + 5));
      }
      writef.write((char*) &amplitudes[0], sizeof(short int) * num_samps);
    }
  }

  Output::Get()->WriteGenerateWF();
  Output::Get()->_ofile->Close();
}
