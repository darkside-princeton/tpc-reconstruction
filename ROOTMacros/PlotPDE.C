#include "Tfile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TF1.h"
#include "TPaveStats.h"
#include <vector>
#include <cmath>
#include <climits>

using namespace std;


void PlotPDE() {

  double voltages[] = {44, 44.5, 45, 45.5, 46, 46.5, 47};
  double min_height[] = {20,30,35,40, 50, 60 , 80};
  double max_height[] = {40,60,70, 80, 90,100, 120};

  vector<double> vDark_count_rate;

  vector<double> vNo_pulses;

  string voltage_fnames[] = {"44", "44dot5", "45", "45dot5", "46", "46dot5", "47"};

  for(int v = 0; v < sizeof(voltages) / sizeof(voltages[0]); v++) {

    TFile *f = new TFile(Form("FinalData/Laser070821/Laser%sv_070821.root", voltage_fnames[v].c_str()));
    TTree *tPulse = (TTree*)f->Get("Pulse_Tree");

    int eventID;
    vector<int> npulses;
    vector<vector<int> *> nafterpulses;
    vector<vector<double> *> start_time;
    vector<vector<double> *> end_time;
    vector<vector<double> *> integral;
    vector<vector<double> *> deconv_integral;
    vector<vector<double> *>  max_val;
    vector<vector<double> *>  min_val;
    vector<vector<double> *>  height;
    vector<vector<double> *>  peak_time;
    vector<vector<double> *>  half_time;
    vector<vector<double> *>  filtered_integral;
    vector<vector<vector<double> > *> ap_peak_times;
    vector<vector<vector<double> > *> ap_heights;

    for(int i_ch = 0; i_ch < 2; i_ch++) {
      npulses.push_back(-1);
      vector<vector<double> > *ap_ch_vec = 0;
      ap_peak_times.push_back(ap_ch_vec);
      ap_heights.push_back(ap_ch_vec);
      vector<double> *ch_vec = 0;
      vector<int> *ch_vecI = 0;
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


    for(int i_ch = 0; i_ch < 2; i_ch++){
      tPulse->SetBranchAddress(Form("start_time_ch%d", i_ch),&start_time[i_ch]);
      tPulse->SetBranchAddress(Form("end_time_ch%d", i_ch),&end_time[i_ch]);
      tPulse->SetBranchAddress(Form("integral_ch%d",i_ch),&integral[i_ch]);
      tPulse->SetBranchAddress(Form("max_val_ch%d", i_ch),&max_val[i_ch]);
      tPulse->SetBranchAddress(Form("min_val_ch%d", i_ch),&min_val[i_ch]);
      tPulse->SetBranchAddress(Form("height_ch%d", i_ch),&height[i_ch]);
      tPulse->SetBranchAddress(Form("peak_time_ch%d", i_ch),&peak_time[i_ch]);
      tPulse->SetBranchAddress(Form("half_time_ch%d", i_ch),&half_time[i_ch]);
      tPulse->SetBranchAddress(Form("nafterpulses_ch%d", i_ch),&nafterpulses[i_ch]);
      tPulse->SetBranchAddress(Form("npulses_ch%d", i_ch),&npulses[i_ch]);
      tPulse->SetBranchAddress("eventID", &eventID);
      tPulse->SetBranchAddress(Form("filtered_integral_ch%d", i_ch), &filtered_integral[i_ch]);
      tPulse->SetBranchAddress(Form("ap_peak_times_ch%d", i_ch), &ap_peak_times[i_ch]);
      tPulse->SetBranchAddress(Form("ap_heights_ch%d", i_ch), &ap_heights[i_ch]);
      tPulse->SetBranchAddress(Form("deconv_integral_ch%d", i_ch), &deconv_integral[i_ch]);
   }

    Int_t nentries = (Int_t)tPulse->GetEntries();
    int dark_counts = 0;

    int total_pulses = 0;
    for (Int_t i = 0; i < nentries; i++) {
      tPulse->GetEntry(i);
      if (npulses[1] > 0) {
        total_pulses++;
        if(peak_time[1]->at(0) < 1.5) dark_counts++;
      }
    }

    vDark_count_rate.push_back((double)dark_counts * 1000000 / (1.5 * total_pulses));
 }

 for(int v = 0; v < sizeof(voltages) / sizeof(voltages[0]); v++) {
   TFile *f = new TFile(Form("FinalData/Laser070821/Laser%sv_pde_070821.root", voltage_fnames[v].c_str()));
   TTree *tPulse = (TTree*)f->Get("Pulse_Tree");

   int eventID;
   vector<int> npulses;
   vector<vector<int> *> nafterpulses;
   vector<vector<double> *> start_time;
   vector<vector<double> *> end_time;
   vector<vector<double> *> integral;
   vector<vector<double> *> deconv_integral;
   vector<vector<double> *>  max_val;
   vector<vector<double> *>  min_val;
   vector<vector<double> *>  height;
   vector<vector<double> *>  peak_time;
   vector<vector<double> *>  half_time;
   vector<vector<double> *>  filtered_integral;
   vector<vector<vector<double> > *> ap_peak_times;
   vector<vector<vector<double> > *> ap_heights;

   for(int i_ch = 0; i_ch < 2; i_ch++) {
     npulses.push_back(-1);
     vector<vector<double> > *ap_ch_vec = 0;
     ap_peak_times.push_back(ap_ch_vec);
     ap_heights.push_back(ap_ch_vec);
     vector<double> *ch_vec = 0;
     vector<int> *ch_vecI = 0;
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


   for(int i_ch = 0; i_ch < 2; i_ch++){
     tPulse->SetBranchAddress(Form("start_time_ch%d", i_ch),&start_time[i_ch]);
     tPulse->SetBranchAddress(Form("end_time_ch%d", i_ch),&end_time[i_ch]);
     tPulse->SetBranchAddress(Form("integral_ch%d",i_ch),&integral[i_ch]);
     tPulse->SetBranchAddress(Form("max_val_ch%d", i_ch),&max_val[i_ch]);
     tPulse->SetBranchAddress(Form("min_val_ch%d", i_ch),&min_val[i_ch]);
     tPulse->SetBranchAddress(Form("height_ch%d", i_ch),&height[i_ch]);
     tPulse->SetBranchAddress(Form("peak_time_ch%d", i_ch),&peak_time[i_ch]);
     tPulse->SetBranchAddress(Form("half_time_ch%d", i_ch),&half_time[i_ch]);
     tPulse->SetBranchAddress(Form("nafterpulses_ch%d", i_ch),&nafterpulses[i_ch]);
     tPulse->SetBranchAddress(Form("npulses_ch%d", i_ch),&npulses[i_ch]);
     tPulse->SetBranchAddress("eventID", &eventID);
     tPulse->SetBranchAddress(Form("filtered_integral_ch%d", i_ch), &filtered_integral[i_ch]);
     tPulse->SetBranchAddress(Form("ap_peak_times_ch%d", i_ch), &ap_peak_times[i_ch]);
     tPulse->SetBranchAddress(Form("ap_heights_ch%d", i_ch), &ap_heights[i_ch]);
     tPulse->SetBranchAddress(Form("deconv_integral_ch%d", i_ch), &deconv_integral[i_ch]);
   }



   Int_t nentries = (Int_t)tPulse->GetEntries();
   int missed = 0;
   int total_pulses = 0;
   for (Int_t i = 0; i < nentries; i++) {
     tPulse->GetEntry(i);
     if (npulses[1] > 0) {
       total_pulses++;
       if(half_time[1]->at(0) > 2) missed++;
     }
   }
   vNo_pulses.push_back(-log((double)(100000 - total_pulses + missed) / 100000));
 }

 TGraph* gDark_counts = new TGraph(sizeof(voltages) / sizeof(voltages[0]), voltages, &vDark_count_rate[0]);
 TGraph* gNoPulse = new TGraph(sizeof(voltages) / sizeof(voltages[0]), voltages, &vNo_pulses[0]);

 gDark_counts->SetMarkerStyle(20);
 gDark_counts->SetTitle("Dark Count Rate; Bias[V]; Counts / s");
 gNoPulse->SetTitle("Efficiency; Bias[V]; #mu");
 gNoPulse->SetMarkerStyle(20);

 TCanvas* c = new TCanvas("PDE", "PDE", 800, 400);
 c->Divide(2,1);
 c->cd(1);
 //gAP_count_full->GetYaxis()->SetRangeUser(100, 800);
 gDark_counts->Draw();
 c->cd(2);
 gNoPulse->Draw("ALP");






}
