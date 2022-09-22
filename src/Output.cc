#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <TFile.h>
#include <TTree.h>
#include "Config.h"
#include "Run.h"
#include "Event.h"
#include "Output.h"
#include "TH1D.h"


using namespace std;
Output *Output::_oh = NULL;

/*--------------------------------------------------------------------*/

void Output::SetupNoise() {
  _ofile = new TFile(Config::Get()->GetParameterS("output_file").c_str(),"recreate");
  _noiseTree = new TTree*[Config::Get()->GetParameterI("num_chans")];
  gHitsThreshold = new TGraph*[Config::Get()->GetParameterI("num_chans")];
  gHitsToT = new TGraph*[Config::Get()->GetParameterI("num_chans")];
  hAvgFFT = new TH1D*[Config::Get()->GetParameterI("num_chans")];

  for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){
    _noiseTree[i_ch] = new TTree(Form("Noise_Tree_%d",i_ch), Form("tree for channel %d",i_ch));
    _noiseTree[i_ch]->Branch("baseline", &baseline_mean);
    _noiseTree[i_ch]->Branch("baseline_dev", &baseline_stdev);
    _noiseTree[i_ch]->Branch("event_ID", &eventID);
    hAvgFFT[i_ch] = new TH1D("avg_FFT", "Average FFT", Config::Get()->GetParameterI("num_samps")/2, 0, Config::Get()->GetParameterD("sampling_rate")*1e3/2);
    for (int bin = 1; bin <= Config::Get()->GetParameterI("num_samps")/2; bin++){hAvgFFT[i_ch]->SetBinContent(bin,0); }
    gHitsThreshold[i_ch] = new TGraph(); gHitsThreshold[i_ch]->SetName("gHitsThreshold");
    gHitsToT[i_ch] = new TGraph(); gHitsToT[i_ch]->SetName("gHitsToT");
    for (int i_t = 0; i_t <= Config::Get()->GetParameterI("max_threshold")/ Config::Get()->GetParameterD("step_size_threshold"); i_t++){
      gHitsThreshold[i_ch]->SetPoint(i_t,i_t * Config::Get()->GetParameterD("step_size_threshold"),0); 
    }
    for (int i_ti = 0; i_ti <= Config::Get()->GetParameterI("max_time")/ Config::Get()->GetParameterD("step_size_time"); i_ti++){
      gHitsToT[i_ch]->SetPoint(i_ti, i_ti *  Config::Get()->GetParameterD("step_size_time"),0);
    }
  }
}

/*--------------------------------------------------------------------*/

void Output::WriteNoise() {
  _ofile->cd();
  for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){
    _noiseTree[i_ch]->Write();
    gHitsThreshold[i_ch]->Write();
    gHitsToT[i_ch]->Write();
    hAvgFFT[i_ch]->Write();
  }
}

/*--------------------------------------------------------------------*/

void Output::SetupSiPM() {
  _ofile = new TFile(Config::Get()->GetParameterS("output_file").c_str(),"recreate");
  _pulseTree = new TTree("Pulse_Tree", "Pulse_Tree");
  gAvgPulse = new TGraph*[Config::Get()->GetParameterI("num_chans")];
  gAvgSPEPulse = new TGraph*[Config::Get()->GetParameterI("num_chans")];
  gAvgWF = new TGraph*[Config::Get()->GetParameterI("num_chans")];


  for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){
    npulses.push_back(-1);
    vector<vector<double> > ap_ch_vec;
    ap_peak_times.push_back(ap_ch_vec);
    ap_heights.push_back(ap_ch_vec);
    vector<double> ch_vec;
    vector<int> ch_vecI;
    start_time.push_back(ch_vec);
    end_time.push_back(ch_vec);
    integral.push_back(ch_vec);
    max_val.push_back(ch_vec);
    min_val.push_back(ch_vec);
    height.push_back(ch_vec);
    peak_time.push_back(ch_vec);
    half_time.push_back(ch_vec);
    filtered_integral.push_back(ch_vec);
    deconv_integral.push_back(ch_vec);
    nafterpulses.push_back(ch_vecI);
  }

  for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){
    _pulseTree->Branch(Form("start_time_ch%d", Config::Get()->GetParameterIvec("ch_ids")[i_ch]),&start_time[i_ch]);
    _pulseTree->Branch(Form("end_time_ch%d", Config::Get()->GetParameterIvec("ch_ids")[i_ch]),&end_time[i_ch]);
    _pulseTree->Branch(Form("integral_ch%d", Config::Get()->GetParameterIvec("ch_ids")[i_ch]),&integral[i_ch]);
    _pulseTree->Branch(Form("max_val_ch%d", Config::Get()->GetParameterIvec("ch_ids")[i_ch]),&max_val[i_ch]);
    _pulseTree->Branch(Form("min_val_ch%d", Config::Get()->GetParameterIvec("ch_ids")[i_ch]),&min_val[i_ch]);
    _pulseTree->Branch(Form("height_ch%d", Config::Get()->GetParameterIvec("ch_ids")[i_ch]),&height[i_ch]);
    _pulseTree->Branch(Form("peak_time_ch%d", Config::Get()->GetParameterIvec("ch_ids")[i_ch]),&peak_time[i_ch]);
    _pulseTree->Branch(Form("half_time_ch%d", Config::Get()->GetParameterIvec("ch_ids")[i_ch]),&half_time[i_ch]);
    _pulseTree->Branch(Form("nafterpulses_ch%d", Config::Get()->GetParameterIvec("ch_ids")[i_ch]),&nafterpulses[i_ch]);
    _pulseTree->Branch(Form("npulses_ch%d", Config::Get()->GetParameterIvec("ch_ids")[i_ch]),&npulses[i_ch]);
    _pulseTree->Branch("eventID", &eventID);
    _pulseTree->Branch(Form("filtered_integral_ch%d", Config::Get()->GetParameterIvec("ch_ids")[i_ch]), &filtered_integral[i_ch]);
    _pulseTree->Branch(Form("ap_peak_times_ch%d", Config::Get()->GetParameterIvec("ch_ids")[i_ch]), &ap_peak_times[i_ch]);
    _pulseTree->Branch(Form("ap_heights_ch%d", Config::Get()->GetParameterIvec("ch_ids")[i_ch]), &ap_heights[i_ch]);
    _pulseTree->Branch(Form("deconv_integral_ch%d", Config::Get()->GetParameterIvec("ch_ids")[i_ch]), &deconv_integral[i_ch]);
  }
}

/*--------------------------------------------------------------------*/

void Output::WriteSiPM() {
  _ofile->cd();
  _pulseTree->Write(_pulseTree->GetName(), TTree::kOverwrite);
  for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){
    gAvgPulse[i_ch]->Write();
    gAvgSPEPulse[i_ch]->Write();
    gAvgWF[i_ch]->Write();
  }
}

/*--------------------------------------------------------------------*/

void Output::SetupTPC() {
  _ofile = new TFile(Config::Get()->GetParameterS("output_file").c_str(),"recreate");
  _clusterTree = new TTree("Cluster_Tree", "Cluster_Tree");
  for(unsigned int i_fp = 0; i_fp < Config::Get()->GetParameterDvec("prompt_windows").size(); i_fp++){
    vector<double> fp_vec;
    cl_fprompt.push_back(fp_vec);
    cl_deconv_fprompt.push_back(fp_vec);
  }
  for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){
    vector<vector<double> > ch_vec;
    pulse_peak_time.push_back(ch_vec);
    pulse_height.push_back(ch_vec);
    vector<double> int_vec;
    cl_deconv_integral.push_back(int_vec);
  } 
  
  _clusterTree->Branch("start_time",&cl_start_time);
  _clusterTree->Branch("end_time",&cl_end_time);
  _clusterTree->Branch("integral",&cl_integral);
  _clusterTree->Branch("max_val",&cl_max_val);
  _clusterTree->Branch("min_val",&cl_min_val);
  _clusterTree->Branch("height",&cl_height);
  _clusterTree->Branch("peak_time",&cl_peak_time);
  _clusterTree->Branch("half_time",&cl_half_time);
  _clusterTree->Branch("eventID", &eventID);
  _clusterTree->Branch("overlap", &overlap);
  _clusterTree->Branch("S1_pulses", &S1_pulses);
  _clusterTree->Branch("S1_charge", &S1_charge);
  _clusterTree->Branch("S2_pulses", &S2_pulses);
  _clusterTree->Branch("S2_charge", &S2_charge);
  _clusterTree->Branch("S1_S2_dt", &S1_S2_dt);
  _clusterTree->Branch("nclusters",&nclusters);
  _clusterTree->Branch("symmetry",&cl_symmetry);
  for(unsigned int i_fp = 0; i_fp < Config::Get()->GetParameterDvec("prompt_windows").size(); i_fp++){
    _clusterTree->Branch(Form("frompt%d", i_fp+1), &cl_fprompt[i_fp]);
    _clusterTree->Branch(Form("deconv_frompt%d", i_fp+1), &cl_deconv_fprompt[i_fp]);
  } 
  for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("num_chans"); i_ch++){  
    _clusterTree->Branch(Form("pulse_peak_time_ch%d", Config::Get()->GetParameterIvec("ch_ids")[i_ch]), &pulse_peak_time[i_ch]);
    _clusterTree->Branch(Form("pulse_height_ch%d", Config::Get()->GetParameterIvec("ch_ids")[i_ch]), &pulse_height[i_ch]);
    _clusterTree->Branch(Form("deconv_integral_ch%d", Config::Get()->GetParameterIvec("ch_ids")[i_ch]), &cl_deconv_integral[i_ch]);
  }
}

/*--------------------------------------------------------------------*/

void Output::WriteTPC() {
  _ofile->cd();
  _clusterTree->Write(_clusterTree->GetName(),TTree::kOverwrite);
  gAvgS1->Write(); 
}


/*--------------------------------------------------------------------*/

void Output::SetupGenerateWF() {
  _ofile = new TFile(Config::Get()->GetParameterS("output_file").c_str(),"recreate");

  for(int i_s2 = 0; i_s2 < Config::Get()->GetParameterI("num_S2") + 1; i_s2++) {
    vector<double> s2_vec;
    rPulseID.push_back(s2_vec);
    rEventID.push_back(s2_vec);
    rStart_time.push_back(s2_vec);
    rIntegral.push_back(s2_vec);
    rHeight.push_back(s2_vec);
    rAfterPulses.push_back(s2_vec);  
  }

  _trueTree = new TTree*[Config::Get()->GetParameterI("generate_channels")];
  for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("generate_channels"); i_ch++){
    _trueTree[i_ch] = new TTree(Form("True_Tree_ch%d", i_ch),Form("True_Tree_ch%d", i_ch));
    _trueTree[i_ch]->Branch("rPulseID", &rPulseID);
    _trueTree[i_ch]->Branch("rEventID", &rEventID);
    _trueTree[i_ch]->Branch("rStart_time", &rStart_time);
    _trueTree[i_ch]->Branch("rIntegral", &rIntegral);
    _trueTree[i_ch]->Branch("rHeight", &rHeight);
    _trueTree[i_ch]->Branch("rAfterPulses", &rAfterPulses);
    _trueTree[i_ch]->Branch("nPE", &nPE);
    _trueTree[i_ch]->Branch("rTime", &rTime);
  }
}

/*--------------------------------------------------------------------*/

void Output::WriteGenerateWF() {
  _ofile->cd();
  for(int i_ch = 0; i_ch < Config::Get()->GetParameterI("generate_channels"); i_ch++){
    _trueTree[i_ch]->Write();
  }
}

