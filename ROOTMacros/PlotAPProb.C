// This file plot afterpulsing probability estimates from Laser runs


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



void PlotAPProb() {
  double voltages[] = {44, 44.5, 45, 45.5, 46, 46.5, 47};
  double min_height[] = {20,30,35,40, 50, 60 , 80};
  double max_height[] = {40,60,70, 80, 90,100, 120};

  vector<double> vAP_count;
  vector<double> vNAP_count;

  vector<double> vAP_count_full;
  vector<double> vNAP_count_full;

  // vector<double> vDark_count_rate;
  //
  // vector<double> vNo_pulses;


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

   double noAP_count = 0;
   int AP_count = 0;
   double total_SPEs = 0;
   int AP_count_full = 0;
   double noAP_count_full = 0;

   for (Int_t i = 0; i < nentries; i++) {
     tPulse->GetEntry(i);
     for(int p = 0; p < npulses[1]; p++) { // for each of the pulses in the event
       if(height[1]->at(0) > min_height[v] &&  height[1]->at(0) < max_height[v]) {
         total_SPEs++; // counting total SPEs
         if(nafterpulses[1]->at(0) == 0) {noAP_count++; noAP_count_full++;}
         AP_count += nafterpulses[1]->at(p); AP_count_full += nafterpulses[1]->at(p);
         if (p == 1 && peak_time[1]->at(p) - peak_time[1]->at(0) < 1 && nafterpulses[1]->at(0) == 0) noAP_count--;
         if (p == 1 && nafterpulses[1]->at(0) == 0) noAP_count_full--;
         if (p != 0 && peak_time[1]->at(p) -  peak_time[1]->at(0) < 1) AP_count++;
         if (p != 0) AP_count_full++;

       }
     }
   }

   vAP_count.push_back((double)AP_count / (double)total_SPEs);
   vAP_count_full.push_back((double)AP_count_full / (double)total_SPEs);
   vNAP_count.push_back(1-noAP_count/ total_SPEs);
   vNAP_count_full.push_back(1-noAP_count_full/total_SPEs);
 }

   TGraph* gAP_count = new TGraph(sizeof(voltages) / sizeof(voltages[0]), voltages, &vAP_count[0]);
   TGraph* gNAP_count = new TGraph(sizeof(voltages) / sizeof(voltages[0]), voltages, &vNAP_count[0]);
   TGraph* gAP_count_full = new TGraph(sizeof(voltages) / sizeof(voltages[0]), voltages, &vAP_count_full[0]);
   TGraph* gNAP_count_full = new TGraph(sizeof(voltages) / sizeof(voltages[0]), voltages, &vNAP_count_full[0]);

   gAP_count->SetMarkerStyle(20); gAP_count->SetMarkerColor(2); gAP_count->SetLineColor(2);
   gNAP_count->SetMarkerStyle(20); gNAP_count->SetMarkerColor(2); gNAP_count->SetLineColor(2);

   gAP_count->SetTitle("Number of Afterpulses"); gAP_count->SetName("< 1 #mus");
   gNAP_count->SetTitle("Fraction of Pulses with Afterpulses"); gAP_count->SetName("< 1 #mus");

   gAP_count_full->SetMarkerStyle(20); gAP_count_full->SetMarkerColor(4); gAP_count_full->SetLineColor(4);
   gNAP_count_full->SetMarkerStyle(20); gNAP_count_full->SetMarkerColor(4); gNAP_count_full->SetLineColor(4);

   gAP_count_full->SetTitle("Number of Afterpulses; Bias [V]; Number of APs / Total SPEs"); gAP_count->SetName("Full Waveform");
   gNAP_count_full->SetTitle("Fraction of Pulses with Afterpulses; Bias [V]; Fraction"); gAP_count->SetName("Full Waveform");

   auto legend = new TLegend(0.1,0.7,0.4,0.9);
   legend->AddEntry(gAP_count,"< 1 #mus");
   legend->AddEntry(gAP_count_full,"Full Waveform");


   TCanvas* c = new TCanvas("Afterpulse Counts", "Afterpulse Counts", 800, 400);
   c->Divide(2,1);
   c->cd(1);
   gAP_count_full->Draw("ALP"); gAP_count->Draw("LP same"); legend->Draw();
   c->cd(2);
   gNAP_count_full->GetYaxis()->SetRangeUser(0.1, 0.6);
   gNAP_count_full->Draw("ALP"); gNAP_count->Draw("LP same"); legend->Draw();

}
