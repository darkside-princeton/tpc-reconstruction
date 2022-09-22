// get radius. this should be anywhere between 0 and 25.

void get_radius(double * x, double * y, double * radius, TRandom3 * tr3, double upper_rad) {
  *x = tr3->Uniform(0, upper_rad);
  *y = tr3->Uniform(0, upper_rad);
  *radius = sqrt(*x * *x + *y * *y);
  if (*radius > upper_rad) {
    get_radius(x, y, radius, tr3, upper_rad);
  }
}

// find Y for a given E. E_Vcm in V/cm, output in mm^-1 electron^-1.
double get_Y(double E_Vcm) {
  double E = E_Vcm / 1000.;
  double Y = 858.7 - 572.6 * E + 117.6 * pow(E, 2) - 6.632 * pow(E, 3);
  return Y / 10.;
}

// upload files
void upload_files(string pfn, TGraph2D * gr, int side) {
  for (int i = 1; i <= 25; i++) {
    ifstream file;
    string filename = pfn + Form("ex_el_field_s%d_r%d.txt", side, i);
    file.open(filename);
    // distance in mm, field in V/m (converted to V/cm)
    double dist, field;
    while (file >> dist >> field) {
      gr->AddPoint(i, dist, field/100.);
    }
  }
}

// given a radius, run a trial and output the number of photons produced.
double run_trial(TGraph2D * gr, double radius, const double step_size) {
  double n = 0;

  // for the whole gas pocket, get the field, increment photons produced based on Y, step up, and repeat
  for (double z = 6; z < 17.9; z += step_size) {
    double E = (gr->Interpolate(radius, z));
    double Y = get_Y(E);
    n += Y * step_size;
  }

  return n;
}

void draw_field(string pfn, int side, TGraph2D * gr, TPad * pad) {
  string t = "Electric Field Strength, Side " + to_string(side) + ";Radius (mm);Distance from grid (mm);Field strength (V/cm)";
  const char * title = t.c_str();
  gr->SetTitle(title);
  upload_files(pfn, gr, side);
  pad->cd();
  gr->Draw("surf1");
  pad->Update();
  gr->GetZaxis()->SetRangeUser(4500, 5500);
  pad->Modified();
}


void s2_amplification_sim(const string folder_path, const int n_trials, const double step_size, const double upper_rad) {
  // set up tree, radius
  TRandom3 * tr3 = new TRandom3(0);
  double xval = 0;
  double * x = &xval;
  double yval = 0;
  double * y = &yval;
  double radval = 0;
  double * radius = &radval;

  get_radius(x, y, radius, tr3, upper_rad);

  double n1, n2, n3, n4;
  int n5, n6, n7, n8;

  // GRAPH FIELDS
  TCanvas * fields  = new TCanvas("fields", "Fields");
  TPad * p1 = new TPad("p1", "Electric Field, Side 1", 0, 0.5, 0.5, 1);
  p1->Draw();
  TPad * p2 = new TPad("p2", "Electric Field, Side 2", 0.5, 0.5, 1, 1);
  p2->Draw();
  TPad * p3 = new TPad("p3", "Electric Field, Side 3", 0, 0, 0.5, 0.5);
  p3->Draw();
  TPad * p4 = new TPad("p4", "Electric Field, Side 4", 0.5, 0, 1, 0.5);
  p4->Draw();
  
  TGraph2D * gr1 = new TGraph2D();
  draw_field(folder_path, 1, gr1, p1);

  TGraph2D * gr2 = new TGraph2D();
  draw_field(folder_path, 2, gr2, p2);
 
  TGraph2D * gr3 = new TGraph2D();
  draw_field(folder_path, 3, gr3, p3);

  TGraph2D * gr4 = new TGraph2D();
  draw_field(folder_path, 4, gr4, p4);


  // make trees to store the data  
  TTree * tree1 = new TTree("tree1", "S2 Signal, Side 1");
  tree1->Branch("x", x);
  tree1->Branch("y", y);
  tree1->Branch("ideal", &n1);
  tree1->Branch("S2", &n5);
  TTree * tree2 = new TTree("tree2", "S2 Signal, Side 2");
  tree2->Branch("x", x);
  tree2->Branch("y", y);
  tree2->Branch("ideal", &n2);
  tree2->Branch("S2", &n6);
  TTree * tree3 = new TTree("tree3", "S2 Signal, Side 3");
  tree3->Branch("x", x);
  tree3->Branch("y", y);
  tree3->Branch("ideal", &n3);
  tree3->Branch("S2", &n7);
  TTree * tree4 = new TTree("tree4", "S2 Signal, Side 4");
  tree4->Branch("x", x);
  tree4->Branch("y", y);
  tree4->Branch("ideal", &n4);
  tree4->Branch("S2", &n8);

  for (int i = 0; i < n_trials; i++) {
    // run trial and enter into the histogram
    get_radius(x, y, radius, tr3, upper_rad);

    n1 = run_trial(gr1, *radius, step_size);
    n2 = run_trial(gr2, *radius, step_size);
    n3 = run_trial(gr3, *radius, step_size);
    n4 = run_trial(gr4, *radius, step_size);

    // adjust using poisson statistics, binomial distribution
    n5 = tr3->Poisson(n1);
    n5 = tr3->Binomial(n5, .219);
    n6 = tr3->Poisson(n2);
    n6 = tr3->Binomial(n6, .219);
    n7 = tr3->Poisson(n3);
    n7 = tr3->Binomial(n7, .219);
    n8 = tr3->Poisson(n4);
    n8 = tr3->Binomial(n8, .219);

    tree1->Fill();
    tree2->Fill();
    tree3->Fill();
    tree4->Fill();
  }

  // draw the histogram
  TCanvas * c2 = new TCanvas("c2", "S2 Amplification Expectation Value");
  TPad * pad1 = new TPad("pad1", "S2 Amplification, Side 1", 0, 0.5, 0.5, 1);
  pad1->Draw();
  TPad * pad2 = new TPad("pad2", "S2 Amplification, Side 2", 0.5, 0.5, 1, 1);
  pad2->Draw();
  TPad * pad3 = new TPad("pad3", "S2 Amplification, Side 3", 0, 0, 0.5, 0.5);
  pad3->Draw();
  TPad * pad4 = new TPad("pad4", "S2 Amplifciation, Side 4", 0.5, 0, 1, 0.5);
  pad4->Draw();

  pad1->cd();
  tree1->Draw("ideal>>h1(300, 110, 140)");
  TH1F *h1 = (TH1F*)gDirectory->Get("h1");
  h1->SetTitle("S2 Amplification, Side 1;Number of Photons;Count");
  pad2->cd();
  tree2->Draw("ideal>>h2(300, 110, 140)");
  TH1F *h2 = (TH1F*)gDirectory->Get("h2");
  h2->SetTitle("S2 Amplification, Side 2;Number of Photons;Count");
  pad3->cd();
  tree3->Draw("ideal>>h3(300, 110, 140)");
  TH1F *h3 = (TH1F*)gDirectory->Get("h3");
  h3->SetTitle("S2 Amplification, Side 3;Number of Photons;Count");
  pad4->cd();
  tree4->Draw("ideal>>h4(300, 110, 140)");
  TH1F *h4 = (TH1F*)gDirectory->Get("h4");
  h4->SetTitle("S2 Amplification, Side 4;Number of Photons;Count");

  TCanvas * c3 = new TCanvas("c3", "Photoelectrons Detected");
  TPad * pad9 = new TPad("pad9", "Photoelectrons Detected, Side 1", 0, 0.5, 0.5, 1);
  pad9->Draw();
  TPad * pad10 = new TPad("pad10", "Photoelectrons Detected, Side 2", 0.5, 0.5, 1, 1);
  pad10->Draw();
  TPad * pad11 = new TPad("pad11", "Photoelectrons Detected, Side 3", 0, 0, 0.5, 0.5);
  pad11->Draw();
  TPad * pad12 = new TPad("pad12", "Photoelectrons Detected, Side 4", 0.5, 0, 1, 0.5);
  pad12->Draw();

  pad9->cd();
  tree1->Draw("S2:sqrt(x*x+y*y)>>h5", "S2 < 100", "COLZ");
  TH1F *h5 = (TH1F*)gDirectory->Get("h5");
  h5->SetTitle("Photoelectrons Detected, Side 1;Radius;Number of Photoelectrons;Count");
  pad10->cd();
  tree2->Draw("S2:sqrt(x*x+y*y)>>h6", "S2 < 100", "COLZ");
  TH1F *h6 = (TH1F*)gDirectory->Get("h6");
  h6->SetTitle("Photoelectrons Detected, Side 2;Radius;Number of Photoelectrons;Count");
  pad11->cd();
  tree3->Draw("S2:sqrt(x*x+y*y)>>h7", "S2 < 100", "COLZ");
  TH1F *h7 = (TH1F*)gDirectory->Get("h7");
  h7->SetTitle("Photoelectrons Detected, Side 3;Radius;Number of Photoelectrons;Count");
  pad12->cd();
  tree4->Draw("S2:sqrt(x*x+y*y)>>h8", "S2 < 100", "COLZ");
  TH1F *h8 = (TH1F*)gDirectory->Get("h8");
  h8->SetTitle("Photoelectrons Detected, Side 4;Radius;Number of Photoelectrons;Count");

  TGraph * mean1 = new TGraph();
  mean1->SetLineColor(1);
  TGraph * std1 = new TGraph();
  std1->SetLineColor(1);
  TGraph * stdmean1 = new TGraph();
  stdmean1->SetLineColor(1);
  TGraph * mean2 = new TGraph();
  mean2->SetLineColor(2);
  TGraph * std2 = new TGraph();
  std2->SetLineColor(2);
  TGraph * stdmean2 = new TGraph();
  stdmean2->SetLineColor(2);
  TGraph * mean3 = new TGraph();
  mean3->SetLineColor(3);
  TGraph * std3 = new TGraph();
  std3->SetLineColor(3);
  TGraph * stdmean3 = new TGraph();
  stdmean3->SetLineColor(3);
  TGraph * mean4 = new TGraph();
  mean4->SetLineColor(4);
  TGraph * std4 = new TGraph();
  std4->SetLineColor(4);
  TGraph * stdmean4 = new TGraph();
  stdmean4->SetLineColor(4);

  for (int i = 1; i < 25; i++) {
    TCanvas * temp = new TCanvas(Form("c1_%d", i), Form("Photoelectrons Detected, Radius = %d", i));
    string str = "S2 < 100 && sqrt(x*x + y*y) >= " + to_string(i) + " && sqrt(x*x + y*y) < (" + to_string(i) + "+1)";
    const char * c = str.c_str();
    TPad * pad13 = new TPad(Form("pad1_%d", i), "", 0, 0.5, 0.5, 1);
    pad13->Draw();
    TPad * pad14 = new TPad(Form("pad2_%d", i), "", 0.5, 0.5, 1, 1);
    pad14->Draw();
    TPad * pad15 = new TPad(Form("pad3_%d", i), "", 0, 0, 0.5, 0.5);
    pad15->Draw();
    TPad * pad16 = new TPad(Form("pad4_%d", i), "", 0.5, 0, 1, 0.5);
    pad16->Draw();

    pad13->cd();
    tree1->Draw(Form("S2>>hist1_%d(100, 0, 100)", i), c);
    pad14->cd();
    tree2->Draw(Form("S2>>hist2_%d(100, 0, 100)", i), c);
    pad15->cd();
    tree3->Draw(Form("S2>>hist3_%d(100, 0, 100)", i), c);
    pad16->cd();
    tree4->Draw(Form("S2>>hist4_%d(100, 0, 100)", i), c);

    TH1F *hist1 = (TH1F*)gDirectory->Get(Form("hist1_%d", i));
    hist1->SetTitle(Form("Photoelectrons Detected at r = %d mm, Side 1;Number of Photoelectrons;Count", i));
    mean1->AddPoint(i, hist1->GetMean());
    std1->AddPoint(i, hist1->GetStdDev());
    stdmean1->AddPoint(i, (hist1->GetStdDev()/hist1->GetMean()));
    TH1F *hist2 = (TH1F*)gDirectory->Get(Form("hist2_%d", i));
    hist2->SetTitle(Form("Photoelectrons Detected at r = %d mm, Side 2;Number of Photoelectrons;Count", i));
    mean2->AddPoint(i, hist2->GetMean());
    std2->AddPoint(i, hist2->GetStdDev());
    stdmean2->AddPoint(i, (hist2->GetStdDev()/hist2->GetMean()));
    TH1F *hist3 = (TH1F*)gDirectory->Get(Form("hist3_%d", i));
    hist3->SetTitle(Form("Photoelectrons Detected at r = %d mm, Side 3;Number of Photoelectrons;Count", i));
    mean3->AddPoint(i, hist3->GetMean());
    std3->AddPoint(i, hist3->GetStdDev());
    stdmean3->AddPoint(i, (hist3->GetStdDev()/hist3->GetMean()));
    TH1F *hist4 = (TH1F*)gDirectory->Get(Form("hist4_%d", i));
    hist4->SetTitle(Form("Photoelectrons Detected at r = %d mm, Side 4;Number of Photoelectrons;Count", i));
    mean4->AddPoint(i, hist4->GetMean());
    std4->AddPoint(i, hist4->GetStdDev());
    stdmean4->AddPoint(i, (hist4->GetStdDev()/hist4->GetMean()));
  }
  TCanvas * c4 = new TCanvas("c4", "Mean");
  TMultiGraph *meansss = new TMultiGraph();
  meansss->Add(mean1);
  meansss->Add(mean2);
  meansss->Add(mean3);
  meansss->Add(mean4);
  meansss->SetTitle("Mean Number of Photoelectrons Detected by Radius;Radius (mm);Mean number of photoelectrons");
  meansss->Draw("A");

  TCanvas * c5 = new TCanvas("c5", "Standard Deviation");
  TMultiGraph *stdstd = new TMultiGraph();
  stdstd->Add(std1);
  stdstd->Add(std2);
  stdstd->Add(std3);
  stdstd->Add(std4);
  stdstd->SetTitle("Standard Deviation in Number of Photoelectrons Detected by Radius;Radius (mm); Std deviation of number of photoelectrons detected");
  stdstd->Draw("A");

  TCanvas * c6 = new TCanvas("c6", "StdDev/Mean");
  TMultiGraph *stdmeans = new TMultiGraph();
  stdmeans->Add(stdmean1);
  stdmeans->Add(stdmean2);
  stdmeans->Add(stdmean3);
  stdmeans->Add(stdmean4);
  stdmeans->SetTitle("Standard Deviation / Mean by Radius;Radius (mm);Std deviation/mean");
  stdmeans->Draw("A");
}

void s2_amplification_sim(const string partial_filename, const int n_trials, const double step_size) {
  s2_amplification_sim(partial_filename, n_trials, step_size, 25);
}

void s2_amplification_sim(const string partial_filename, const int n_trials) {
  s2_amplification_sim(partial_filename, n_trials, 0.1, 25);
}
