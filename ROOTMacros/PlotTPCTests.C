// Plots a comparison of reconstructed time and charge information from
// the cluster finding algorithm and true information


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

void PlotTPCTests() {
  TFile *f = new TFile("FinalData/ER50keV/clusterER50keV.root");
  TTree *tCluster_electron = (TTree*)f->Get("Cluster_Tree");

  TFile *fTrue = new TFile("FinalData/ER50keV/trueTreeER50keV.root");
  TTree** trueTree = new TTree*[8];

  double fprompts[] = {0.09, 0.135, 0.180, 0.225, 0.270, 0.315};


  int nclusters;
  int eventID;
  int overlap;
  int S1_pulses;
  double S1_charge;
  int S2_pulses;
  double S2_charge;
  double S1_S2_dt;
  vector<double> *cl_start_time = 0;
  vector<double> *cl_end_time = 0;
  vector<double> *cl_integral = 0;
  vector<double> *cl_filtered_integral = 0;
  vector<double> *cl_max_val = 0;
  vector<double> *cl_min_val = 0;
  vector<double> *cl_height = 0;
  vector<double> *cl_peak_time = 0;
  vector<double> *cl_half_time = 0;
  vector<double> *cl_symmetry = 0;
  vector<vector<double> *> cl_fprompt;
  vector<vector<vector<double> > *> pulse_peak_time;
  vector<vector<vector<double> > *> pulse_height;
  vector<vector<double> *> cl_deconv_integral;
  vector<vector<double> *> cl_deconv_fprompt;
  vector<vector<double> *> cl_filtered_fprompt;
  for(unsigned int i_fp = 0; i_fp < 6; i_fp++){
    vector<double> *fp_vec = 0;
    cl_fprompt.push_back(fp_vec);
    cl_deconv_fprompt.push_back(fp_vec);
    cl_filtered_fprompt.push_back(fp_vec);
  }
  for(int i_ch = 0; i_ch < 8; i_ch++){
   vector<vector<double> > *ch_vec = 0;
   pulse_peak_time.push_back(ch_vec);
   pulse_height.push_back(ch_vec);
   vector<double> *int_vec = 0;
   cl_deconv_integral.push_back(int_vec);
  }


  tCluster_electron->SetBranchAddress("start_time",&cl_start_time);
  tCluster_electron->SetBranchAddress("end_time",&cl_end_time);
  tCluster_electron->SetBranchAddress("integral",&cl_integral);
  tCluster_electron->SetBranchAddress("max_val",&cl_max_val);
  tCluster_electron->SetBranchAddress("min_val",&cl_min_val);
  tCluster_electron->SetBranchAddress("height",&cl_height);
  tCluster_electron->SetBranchAddress("peak_time",&cl_peak_time);
  tCluster_electron->SetBranchAddress("half_time",&cl_half_time);
  tCluster_electron->SetBranchAddress("eventID", &eventID);
  tCluster_electron->SetBranchAddress("S1_pulses", &S1_pulses);
  tCluster_electron->SetBranchAddress("S1_charge", &S1_charge);
  tCluster_electron->SetBranchAddress("S2_pulses", &S2_pulses);
  tCluster_electron->SetBranchAddress("S2_charge", &S2_charge);
  tCluster_electron->SetBranchAddress("S1_S2_dt", &S1_S2_dt);
  tCluster_electron->SetBranchAddress("nclusters",&nclusters);
  tCluster_electron->SetBranchAddress("symmetry",&cl_symmetry);
  tCluster_electron->SetBranchAddress("overlap",&overlap);

  for(unsigned int i_fp = 0; i_fp < 6; i_fp++){
    tCluster_electron->SetBranchAddress(Form("frompt%d", i_fp+1), &cl_fprompt[i_fp]);
    tCluster_electron->SetBranchAddress(Form("deconv_frompt%d", i_fp+1), &cl_deconv_fprompt[i_fp]);
  }
  for(int i_ch = 0; i_ch < 8; i_ch++){
    tCluster_electron->SetBranchAddress(Form("pulse_peak_time_ch%d", i_ch), &pulse_peak_time[i_ch]);
    tCluster_electron->SetBranchAddress(Form("pulse_height_ch%d", i_ch), &pulse_height[i_ch]);
    tCluster_electron->SetBranchAddress(Form("deconv_integral_ch%d", i_ch), &cl_deconv_integral[i_ch]);
  }


  vector<int> *nPE = 0;
  vector<double> *rTime = 0;
  vector<vector<double> >* rEventID = 0;
  vector<vector<double> >* rPulseID = 0;
  vector<vector<double> >* rStart_time = 0;
  vector<vector<double> >* rIntegral = 0;
  vector<vector<double> >* rHeight = 0;
  vector<vector<double> >* rAfterPulses = 0;

  for(int i_ch = 0; i_ch < 8; i_ch++){
    trueTree[i_ch] = (TTree*)fTrue->Get(Form("True_Tree_ch%d", i_ch));
    trueTree[i_ch]->SetBranchAddress("rPulseID", &rPulseID);
    trueTree[i_ch]->SetBranchAddress("rEventID", &rEventID);
    trueTree[i_ch]->SetBranchAddress("rStart_time", &rStart_time);
    trueTree[i_ch]->SetBranchAddress("rIntegral", &rIntegral);
    trueTree[i_ch]->SetBranchAddress("rHeight", &rHeight);
    trueTree[i_ch]->SetBranchAddress("rAfterPulses", &rAfterPulses);
    trueTree[i_ch]->SetBranchAddress("nPE", &nPE);
    trueTree[i_ch]->SetBranchAddress("rTime", &rTime);

  }

  TH1D *hS1Charge   = new TH1D("hS1Charge","S1 Charge; Charge[PE]; Counts",100,0,0);
  TH1D *hS2Charge   = new TH1D("hS2Charge","S2 Charge; Charge[PE]; Counts",100,0,0);

  TH1D *hS1Charge_deconv   = new TH1D("hS1Charge_deconv","S1 Charge Deconv.; Charge[PE]; Counts",100,0,0);
  TH1D *hS2Charge_deconv   = new TH1D("hS2Charge_deconv","S2 Charge Deconv.; Charge[PE]; Counts",100,0,0);

  TH1D *hS1TrueCharge = new TH1D("hS1TrueCharge","S1 True Charge",100,0,0);
  TH1D *hS2TrueCharge = new TH1D("hS2TrueCharge","S2 True Charge",100,0,0);

  TH2D *hNs1PE = new TH2D("hNs1PE","Reconstructed PEs S1; Reconstructed PEs; True PEs",75,0,0,75,0,0);
  TH2D *hNs2PE = new TH2D("hNs2PE","Reconstructed PEs S2; Reconstructed PEs; True PEs",75,0,0,75,0,0);

  TH2D *hNs1PE_deconv = new TH2D("hNs1PE_deconv","Reconstructed PEs S1 Deconv.; Reconstructed PEs; True PEs",75,0,0,75,0,0);
  TH2D *hNs2PE_deconv = new TH2D("hNs2PE_deconv","Reconstructed PEs S2 Deconv.; Reconstructed PEs; True PEs",75,0,0,75,0,0);

  TH2D *hZPos = new TH2D("hZPos","Reconstructed z; Reconstructed z [cm]; True z [cm]",75,0,0,75,0,0);

  TH1D **hFPrompt_s1 = new TH1D*[6];
  TH1D **hFPrompt_s2 = new TH1D*[6];

  TH1D **hFPrompt_s1_deconv = new TH1D*[6];
  TH1D **hFPrompt_s2_deconv = new TH1D*[6];

  TH1D **hFPrompt_s1_true = new TH1D*[6];
  TH1D **hFPrompt_s2_true = new TH1D*[6];

  TH1D *hFPrompt_diff   = new TH1D("hFPrompt_diff","Fprompt mean difference",6,0,0);
  TH1D *hFPrompt_diff_deconv   = new TH1D("hFPrompt_diff_deconv","Fprompt mean difference",6,0,0);
  TH1D *hFPrompt_diff_true   = new TH1D("hFPrompt_diff_true","Fprompt mean difference",6,0,0);

  for(int i = 0; i < 6; i++) {
    hFPrompt_s1[i] = new TH1D(Form("hFromptS1%d", i), Form("%d #tau; Fprompt; Counts", i + 2), 200, 0, 1);
    hFPrompt_s2[i] = new TH1D(Form("hFromptS2%d", i), Form("%d #tau; Fprompt; Counts", i + 2), 200, 0, 1);

    hFPrompt_s1_deconv[i] = new TH1D(Form("hFromptS1_deconv%d", i), Form("%d #tau; Fprompt; Counts", i + 2), 200, 0, 1);
    hFPrompt_s2_deconv[i] = new TH1D(Form("hFromptS2_deconv%d", i), Form("%d #tau; Fprompt; Counts", i + 2), 500, 0, 1);

    hFPrompt_s1_true[i] = new TH1D(Form("hFromptS1_true%d", i), Form("%d #tau; Fprompt; Counts", i + 2), 200, 0, 1);
    hFPrompt_s2_true[i] = new TH1D(Form("hFromptS2_true%d", i), Form("%d #tau; Fprompt; Counts", i + 2), 500, 0, 1);
  }


  // Filling histograms
  Int_t nentries = (Int_t)tCluster_electron->GetEntries();

  for (Int_t i = 0; i < nentries; i++) {
    tCluster_electron->GetEntry(i);
    if (nclusters == 2 && !overlap) {
      hS1Charge->Fill(cl_integral->at(0));
      hS2Charge->Fill(cl_integral->at(1));
      double totalTruePEs_S1 = 0;
      double totalTruePEs_S2 = 0;
      double totalS1charge = 0;
      double totalS2charge = 0;
      double min_s1_time = 100;
      double min_s2_time = 100;
      double S1_deconv_charge = 0;
      double S2_deconv_charge = 0;
      for(int i_ch = 0; i_ch < 8; i_ch++){
        S1_deconv_charge += cl_deconv_integral[i_ch]->at(0);
        S2_deconv_charge += cl_deconv_integral[i_ch]->at(1);
        trueTree[i_ch]->GetEntry(i);
        totalTruePEs_S1 += nPE->at(0);
        totalTruePEs_S2 += nPE->at(1);
        for (unsigned int j = 0; j < rStart_time->at(0).size(); j++) {
          totalS1charge += rIntegral->at(0).at(j);
        }
        for (unsigned int j = 0; j < rEventID->at(1).size(); j++) {
          totalS2charge += rIntegral->at(1).at(j);
        }
      }
      // Calculating true fprompts
      for(int i_fp = 0; i_fp < 6; i_fp++){
        double total_fprompt_s1 = 0;
        double total_fprompt_s2 = 0;
        for(int i_ch = 0; i_ch < 8; i_ch++){
          trueTree[i_ch]->GetEntry(i);
          for (unsigned int j = 0; j < rEventID->at(0).size(); j++) {
            if(rStart_time->at(0).at(j) - rTime->at(0) < fprompts[i_fp]) total_fprompt_s1++;
          }
          for (unsigned int j = 0; j < rEventID->at(1).size(); j++) {
            if(rStart_time->at(1).at(j) - rTime->at(1) < fprompts[i_fp]) total_fprompt_s2 ++;
          }
        }
        hFPrompt_s1_true[i_fp]->Fill(total_fprompt_s1 / totalTruePEs_S1);
        hFPrompt_s2_true[i_fp]->Fill(total_fprompt_s2 / totalTruePEs_S2);
      }

      hS1Charge_deconv->Fill(S1_deconv_charge);
      hS2Charge_deconv->Fill(S2_deconv_charge);
      hS1TrueCharge->Fill(totalS1charge / 367);
      hS2TrueCharge->Fill(totalS2charge / 367);
      hNs1PE->Fill(cl_integral->at(0), totalTruePEs_S1);
      hNs2PE->Fill(cl_integral->at(1), totalTruePEs_S2);
      hNs1PE_deconv->Fill(S1_deconv_charge, totalTruePEs_S1);
      hNs2PE_deconv->Fill(S2_deconv_charge, totalTruePEs_S2);
      hZPos->Fill(S1_S2_dt/10, (rTime->at(1) -rTime->at(0))/10);

      //Filling Fprompts
      double* ER_s1 = new double[6];
      for(int i = 0; i < 6; i++) {
        hFPrompt_s1[i]->Fill(cl_fprompt[i]->at(0)); ER_s1[i] = cl_fprompt[i]->at(0);
        hFPrompt_s2[i]->Fill(cl_fprompt[i]->at(1));
        hFPrompt_s1_deconv[i]->Fill(cl_deconv_fprompt[i]->at(0));
        hFPrompt_s2_deconv[i]->Fill(cl_deconv_fprompt[i]->at(1));
      }
     }
  }

  TCanvas *cTruthComp = new TCanvas("Truth Comparison", "Truth Comparison", 800, 800);
  TLine* lNs1PE = new TLine(500, 500, 700, 700); lNs1PE->SetLineColor(2);
  TLine* lNs2PE = new TLine(12800, 12800, 13800, 13800); lNs2PE->SetLineColor(2);
  TLine* lZPos = new TLine(1, 1, 7, 7); lZPos->SetLineColor(2);

  cTruthComp->Divide(2,2);
  cTruthComp->cd(1); hS1Charge->Draw();
  hS1TrueCharge->SetLineColor(2); hS1TrueCharge->Draw("same");
  cTruthComp->cd(2); hS2Charge->Draw(); hS2TrueCharge->SetLineColor(2); hS2TrueCharge->Draw("same");
  cTruthComp->cd(3); hNs1PE->Draw("Colz"); lNs1PE->Draw("same");
  cTruthComp->cd(4); hNs2PE->Draw("Colz"); lNs2PE->Draw("same");

  TCanvas *cZPos = new TCanvas("Position Reconstruction", "Position Reconstruction", 400, 400);
  hZPos->Draw("Colz"); lZPos->Draw("same");

  TCanvas *cFPrompt_deconv = new TCanvas("fprompt Comparison Deconv", "fprompt Comparison Deconv", 1400, 800);
  cFPrompt_deconv->Divide(3,2);
  for(int i = 0; i < 6; i++) {
    cFPrompt_deconv->cd(i+1);
    hFPrompt_s1_deconv[i]->SetStats(0); hFPrompt_s2_deconv[i]->SetStats(0);
    hFPrompt_s1_deconv[i]->Draw(); hFPrompt_s2_deconv[i]->Draw("sames");
    hFPrompt_s2_true[i]->SetStats(0); hFPrompt_s1_true[i]->SetStats(0);
    hFPrompt_s2_true[i]->SetLineColor(2);
    hFPrompt_s1_true[i]->SetLineColor(2); hFPrompt_s1_true[i]->Draw("same"); hFPrompt_s2_true[i]->Draw("sames");
  }

  TCanvas *cFPrompt = new TCanvas("fprompt Comparison", "fprompt Comparison", 1400, 800);
  cFPrompt->Divide(3,2);
  for(int i = 0; i < 6; i++) {
    cFPrompt->cd(i+1);
    hFPrompt_s1[i]->SetStats(0); hFPrompt_s2[i]->SetStats(0);
    hFPrompt_s1[i]->Draw(); hFPrompt_s2[i]->Draw("same");
    hFPrompt_s2_true[i]->SetStats(0); hFPrompt_s1_true[i]->SetStats(0);
    hFPrompt_s2_true[i]->SetLineColor(2);
    hFPrompt_s1_true[i]->SetLineColor(2); hFPrompt_s1_true[i]->Draw("same"); hFPrompt_s2_true[i]->Draw("sames");
  }

  TCanvas *cFPromptDiff = new TCanvas("fprompt Diff", "fprompt Diff", 400, 400);
  const char*fpromptdiff_labels[18] = {"2 #tau", "3 #tau", "4 #tau", "5 #tau", "6 #tau", "7 #tau"};

    for(int i = 0; i < 6; i++) {
      hFPrompt_diff->Fill(fpromptdiff_labels[i], (hFPrompt_s1[i]->GetMean(1) - hFPrompt_s2[i]->GetMean(1))/ sqrt(pow(hFPrompt_s1[i]->GetRMS(1),2) + pow(hFPrompt_s2[i]->GetRMS(1),2)));
      hFPrompt_diff_deconv->Fill(fpromptdiff_labels[i], (hFPrompt_s1_deconv[i]->GetMean(1) - hFPrompt_s2_deconv[i]->GetMean(1))/ sqrt(pow(hFPrompt_s1_deconv[i]->GetRMS(1),2) + pow(hFPrompt_s2_deconv[i]->GetRMS(1),2)));
      hFPrompt_diff_true->Fill(fpromptdiff_labels[i], (hFPrompt_s1_true[i]->GetMean(1) - hFPrompt_s2_true[i]->GetMean(1))/ sqrt(pow(hFPrompt_s1_true[i]->GetRMS(1),2) + pow(hFPrompt_s2_true[i]->GetRMS(1),2)));
    }

  hFPrompt_diff->Draw("Hist");
  hFPrompt_diff_deconv->SetLineColor(1);
  hFPrompt_diff_deconv->Draw("Hist same");
  hFPrompt_diff_true->SetLineColor(2);
  hFPrompt_diff_true->SetLineStyle(9); hFPrompt_diff_true->Draw("Hist same");

  TCanvas *cTruthCompDeconv = new TCanvas("Truth Comparison Deconv", "Truth Comparison Deconv", 800, 800);
  cTruthCompDeconv->Divide(2,2);
  cTruthCompDeconv->cd(1); hS1Charge_deconv->Draw(); hS1TrueCharge->SetLineColor(2); hS1TrueCharge->Draw("same");
  cTruthCompDeconv->cd(2); hS2Charge_deconv->Draw(); hS2TrueCharge->SetLineColor(2); hS2TrueCharge->Draw("same");
  cTruthCompDeconv->cd(3); hNs1PE_deconv->Draw("Colz"); lNs1PE->Draw("same");
  cTruthCompDeconv->cd(4); hNs2PE_deconv->Draw("Colz"); lNs2PE->Draw("same");


}
