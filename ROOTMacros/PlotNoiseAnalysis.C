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

static void rescaleaxis(TGraph *g,double scale) {
  int N=g->GetN();
  double *y=g->GetX();
  int i=0;
  while(i<N) {y[i]=y[i]*scale; i=i+1;}
  g->GetHistogram()->Delete();
  g->SetHistogram(0);
}

void PlotNoiseAnalysis() {
  TFile *f = new TFile("NoiseData/noise_warm.root");
  TTree *tNoise = (TTree*)f->Get("Noise_Tree_0");
  TGraph* gHitThreshold = (TGraph*)f->Get("gHitsThreshold;1");
  TGraph* gHitTime = (TGraph*)f->Get("gHitsToT;1");
  TH1D* hAvgFFT = (TH1D*)f->Get("avg_FFT;1");

  double baseline_mean;
  double baseline_stdev;
  int event_ID;

  tNoise->SetBranchAddress("baseline", &baseline_mean);
  tNoise->SetBranchAddress("baseline_dev", &baseline_stdev);
  tNoise->SetBranchAddress("event_ID", &event_ID);

  TGraph* gEventBaseline = new TGraph();
  TGraph* gEventBaselineDev = new TGraph();
  TH1D *hBaseline   = new TH1D("hBaseline","Baseline Mean",50,0,0);
  TH1D *hBaselineDev   = new TH1D("hBaselineDev","Baseline Std. Dev.",50,0,0);

  Int_t nentries = (Int_t)tNoise->GetEntries();

  for (Int_t i = 0; i < nentries; i++) {
    tNoise->GetEntry(i);
    gEventBaseline->SetPoint(i, event_ID, baseline_mean);
    gEventBaselineDev->SetPoint(i, event_ID, baseline_stdev);
    hBaseline->Fill(baseline_mean);
    hBaselineDev->Fill(baseline_stdev);
  }

  //Plot Formatting
  gEventBaseline->SetTitle("Baseline Mean; Event Number; Baseline");
  gEventBaselineDev->SetTitle("Baseline Std. Dev.; Event Number; Baseline Std. Dev.");
  gHitThreshold->SetTitle("Hit Rate vs. Threshold; Threshold; Hit Rate [Hits / #mus]");
  rescaleaxis(gHitTime, 1e-3 / 0.25);
  gHitTime->GetXaxis()->SetRangeUser(0,0.1);
  gHitTime->SetTitle("Hit Rate vs. Time over Threshold; Time [#mus]; Hit Rate [Hits / #mus]");
  hAvgFFT->SetTitle("Avg. FFT; Frequency [MHz]; Amplitude");
  hAvgFFT->SetStats(FALSE);

  TCanvas* c = new TCanvas("Noise Analysis","Noise Analysis",1400,800);
  c->Divide(3,2);

  c->cd(1); gEventBaseline->Draw("AL");
  c->cd(2); gEventBaselineDev->Draw("AL");
  c->cd(4); hBaseline->Draw();
  c->cd(5); hBaselineDev->Draw();
  c->cd(3); gHitThreshold->Draw("AL");
  c->cd(6); gHitTime->Draw("AL");

  TCanvas* c1 = new TCanvas("FFT","FFT",400,400);
  hAvgFFT->Draw();


}
