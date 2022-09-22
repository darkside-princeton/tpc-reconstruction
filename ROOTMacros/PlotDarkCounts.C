// Generates distributions for dark count data taken with channel 0

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


static Double_t findMaxPosition(TGraph *G) {
    Double_t x, y;
    G->GetPoint(0, x, y);
    Double_t max_x = x, max_y = y;
    for(int i = 1; i < G->GetN(); i++)
    {
        G->GetPoint(i, x, y);
        if(y > max_y) {
           max_x = x;
           max_y = y;
        }
    }
    return max_x;
}


void PlotDarkCounts() {
  TFile *f = new TFile("FinalData/DarkCounts/WarmDarkCounts53v_310721.root");
  TTree *tPulse = (TTree*)f->Get("Pulse_Tree");

  TGraph* avg_pulse = (TGraph*)f->Get("gAvgSPEPulse_ch0");

  vector<int> npulses;
  int eventID;
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

  for(int i_ch = 0; i_ch < 1; i_ch++){
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

  for(int i_ch = 0; i_ch < 1; i_ch++){
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

 //Pulses with no afterpulses
 TH1D *hStartTime   = new TH1D("hStartTime","Start Time Distribution; Time [#mus]; Counts",300,0,10);
 TH1D *hLength   = new TH1D("hLength","Pulse Length Distribution; Time [#mus]; Counts",100,0,1);
 TH1D *hIntegral   = new TH1D("hIntegral","Charge Distribution; Charge [ADC Counts]; Counts",300,0,0);
 TH1D *hFilteredIntegral   = new TH1D("hFilteredIntegral","Charge Distribution",300,0,0);
 TH1D *hHeight   = new TH1D("hHeight","Height Distribution; Height [ADC Counts]; Counts",300,0,0);
 TH2D *hHeightIntegral = new TH2D("hHeightIntegral","Height vs. Charge; Height [ADC Counts]; Charge [ADC Counts]",300,0,0,300,0,0);
 TH2D *hLengthIntegral = new TH2D("hLengthIntegral","Length vs. Charge; Pulse Length [#mus]; Charge [ADC Counts]", 200,0,0,200,0,0);


 //TH1D *hPeakTime = new TH1D("hPeakTime","Start Time Distribution",100,0,0);

 //Pulses with afterpulses
 TH1D *hStartTime_ap   = new TH1D("hStartTime_ap","Start Time Distribution",300,0,10);
 TH1D *hLength_ap   = new TH1D("hLength_ap","Pulse Length Distribution",100,0,3);
 TH1D *hIntegral_ap   = new TH1D("hIntegral_ap","Charge Distribution",300,0,0);
 TH1D *hFilteredIntegral_ap   = new TH1D("hFilteredIntegral_ap","Charge Distribution",300,0,0);
 TH1D *hHeight_ap   = new TH1D("hHeight_ap","Height Distribution",300,0,0);
 TH2D *hHeightIntegral_ap = new TH2D("hHeightIntegral_ap","Height vs Integral",300,0,0,300,0,0);
 TH2D *hLengthIntegral_ap = new TH2D("hLengthIntegral_ap","Length vs Integral",200,0,0,200,0,0);

 //After-pulse Plots
 TH1D *hAPTimeDiff = new TH1D("hAPTimeDiff","Afterpulse Time Distribution; Time Since Previous Pulse [#mus]; Counts",100,0,1);
 TH1D *hAPHeight = new TH1D("hAPHeight","Afterpulse Height Distribution; Height [ADC Counts]; Counts",100,0,0);
 TH2D *hAPTimeHeight = new TH2D("hAPTimeHeight","Afterpulse Time vs. Height Distribution; Time Since Previous Pulse [#mus]; Height [ADC Counts]",200,0,0,200,0,0);

 Int_t nentries = (Int_t)tPulse->GetEntries();

 for (Int_t i = 0; i < nentries; i++) {
   tPulse->GetEntry(i);

   for(int p = 0; p < npulses[0]; p++) { // for each of the pulses in the event
     if(nafterpulses[0]->at(p) == 0) {
       hStartTime->Fill(start_time[0]->at(p));
       hLength->Fill(end_time[0]->at(p) - start_time[0]->at(p));
       hIntegral->Fill(integral[0]->at(p));
       hHeight->Fill(height[0]->at(p));
       hHeightIntegral->Fill(height[0]->at(p), integral[0]->at(p));
       hLengthIntegral->Fill(end_time[0]->at(p) - start_time[0]->at(p), integral[0]->at(p));
     }

     if(nafterpulses[0]->at(p) > 0) {
       hStartTime_ap->Fill(start_time[0]->at(p));
       hLength_ap->Fill(end_time[0]->at(p) - start_time[0]->at(p));
       hIntegral_ap->Fill(integral[0]->at(p));
       hHeight_ap->Fill(height[0]->at(p));
       hHeightIntegral_ap->Fill(height[0]->at(p), integral[0]->at(p));
       hLengthIntegral_ap->Fill(end_time[0]->at(p) - start_time[0]->at(p), integral[0]->at(p));
     }
     for(int j = 0; j < ap_peak_times[0]->at(p).size(); j++) {
       if(j == 0) {
         hAPTimeDiff->Fill(ap_peak_times[0]->at(p).at(j) - peak_time[0]->at(p));
         hAPTimeHeight->Fill(ap_peak_times[0]->at(p).at(j) - peak_time[0]->at(p), ap_heights[0]->at(p).at(j));
       }
       else if (j > 0) {
         hAPTimeDiff->Fill(ap_peak_times[0]->at(p).at(j) - peak_time[0]->at(p));
         hAPTimeHeight->Fill(ap_peak_times[0]->at(p).at(j) - peak_time[0]->at(p), ap_heights[0]->at(p).at(j));
       }
       hAPHeight->Fill(ap_heights[0]->at(p).at(j));
     }
   }
 }


   TCanvas* c = new TCanvas("w/o After-Pulses","w/o After-Pulses",1400,800);
   c->Divide(3,2);

   c->cd(1); hStartTime->Draw();
   c->cd(2); hLength->Draw();
   c->cd(3); hHeight->Draw();
   c->cd(4); hIntegral->Fit("gaus", "", "", 0, 500); gStyle->SetOptFit(1); hIntegral->Draw();
   c->cd(5); hLengthIntegral->Draw("Colz");
   c->cd(6); hHeightIntegral->Draw("Colz");

   TCanvas* c1 = new TCanvas("w/ After-Pulses","w/ After-Pulses",1400,800);
   c1->Divide(3,2);

   c1->cd(1); hStartTime_ap->Draw();
   c1->cd(2); hLength_ap->Draw();
   c1->cd(3); hHeight_ap->Draw();
   c1->cd(4); hIntegral_ap->Draw();
   c1->cd(5); hLengthIntegral_ap->Draw("Colz");
   c1->cd(6); hHeightIntegral_ap->Draw("Colz");

   TCanvas* c2 = new TCanvas("After-Pulses","After-Pulses",1400,400);
   c2->Divide(3,1);

   c2->cd(2); hAPTimeDiff->Draw();
   c2->cd(1); hAPHeight->Draw();
   c2->cd(3); hAPTimeHeight->Draw("Colz");

   TH1F* havgPulse = new TH1F("havgPulse", "havgPulse", avg_pulse->GetN() + 1, 0, avg_pulse->GetXaxis()->GetXmax());
   for(int i_s = 0 ; i_s < avg_pulse->GetN(); i_s++){
     Double_t x, y;
     avg_pulse->GetPoint(i_s, x, y);
     havgPulse->SetBinContent(i_s+1, y);
  }

  TCanvas* c3 = new TCanvas("Avg Pulse Shape", "Avg Pulse Shape", 400,400);
  avg_pulse->SetTitle("Average Pulse Shape");
  TF1 *f1 = new TF1("f1","expo",findMaxPosition(avg_pulse) + 0.01, avg_pulse->GetPointX(avg_pulse->GetN()-1));
  avg_pulse->Fit("f1", "R");
  double time_constant = -1/f1->GetParameter(1) * 1000;
  avg_pulse->Draw();
  c3->Update();

  //Adding the time constant to the stat box:
  TPaveStats *ps = (TPaveStats*)c3->GetPrimitive("stats");
  TList *listOfLines = ps->GetListOfLines();
  TLatex *myt = new TLatex(0,0,TString::Format("#tau = %g ns", time_constant));
  myt ->SetTextFont(42);
  myt ->SetTextSize(0.035);
  myt ->SetTextColor(kRed);
  listOfLines->Add(myt);
  avg_pulse->GetListOfFunctions()->Remove(ps); ps->Draw();
  avg_pulse->PaintStats(0);
  c3->Modified();
  c3->Update();










}
