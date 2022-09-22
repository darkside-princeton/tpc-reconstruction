// Generates a plot of bias voltage vs. SPE charge, as determined from
// laser runs. Used to find SiPM breakdown voltage.


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


void PlotBreakdown() {

  // Charges found from fits to SPE distribution in Laser runs
  double voltages[] = {44, 44.5, 45, 45.5, 46, 46.5, 47};
  double charges[] = {367, 496.3, 637.8, 787.8, 943.9, 1096, 1263};

  TCanvas* cOVCharge = new TCanvas("Overvoltage vs. Charge", "Overvoltage vs. Charge", 400, 400);
  TGraph* gOverV = new TGraph(7, voltages, charges);
  cOVCharge->cd();
  gOverV->SetTitle("Breakdown Voltage; Bias [V]; SPE Charge [ADC Counts]");
  gOverV->SetMarkerStyle(20);
  gOverV->Fit("pol1","", "", 40, 50);
  gOverV->GetXaxis()->SetLimits(0, 48);
  gOverV->GetYaxis()->SetLimits(0, 1400);
  gOverV->Draw();


}
