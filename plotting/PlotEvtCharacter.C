#include"BlairUtils.C"
#include <vector>
#include <string>
#include <cmath>
#include <yaml-cpp/yaml.h>

const char *fdata = "../histMakers/results/evt_characterdata_aa.root";
const char *fphoton = "../histMakers/results/evt_characterphoton12_aa.root";
const char *fjet = "../histMakers/results/evt_characterjet20_aa.root";
const char *plot_outdir = "figs/EvtCharacter";

const char *legend_label_data = "Data";
const char *legend_label_photon = "Pythia #gamma Overlay";
const char *legend_label_jet = "Pythia jet Overlay";
const char *config_yaml_path = "../histMakers/nom.yaml";


void PlotEvtCharacter()
{
  SetsPhenixStyle();
  gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");

  std::vector<float> pT_bins;
  try {
    YAML::Node configYaml = YAML::LoadFile(config_yaml_path);
    if (configYaml["analysis"] && configYaml["analysis"]["pT_bins"]) {
      pT_bins = configYaml["analysis"]["pT_bins"].as<std::vector<float>>();
    }
    std::cout << "Loaded config: " << config_yaml_path << std::endl;
  } catch (const std::exception &e) {
    std::cerr << "Warning: failed to read " << config_yaml_path << " (" << e.what() << ")" << std::endl;
    std::cerr << "Proceeding without pT-bin labels from config." << std::endl;
  }

  TFile *fd = TFile::Open(fdata, "READ");
  TFile *fp = TFile::Open(fphoton, "READ");
  TFile *fj = TFile::Open(fjet, "READ");

  if (!fd || fd->IsZombie()) { std::cerr << "Cannot open data file: " << fdata << std::endl; return; }
  if (!fp || fp->IsZombie()) { std::cerr << "Cannot open photon10_aa file: " << fphoton << std::endl; return; }
  if (!fj || fj->IsZombie()) { std::cerr << "Cannot open jet10_aa file: " << fjet << std::endl; return; }

  TH1 *hd_cent   = (TH1 *)fd->Get("h_centrality");
  TH1 *hd_vtx    = (TH1 *)fd->Get("h_vertexz");
  TH1 *hd_cent10 = (TH1 *)fd->Get("h_centrality_with_cluster10");
  TH1 *hd_pt     = (TH1 *)fd->Get("h_cluster_pt_above10");

  TH1 *hp_cent   = (TH1 *)fp->Get("h_centrality");
  TH1 *hp_vtx    = (TH1 *)fp->Get("h_vertexz");
  TH1 *hp_cent10 = (TH1 *)fp->Get("h_centrality_with_cluster10");
  TH1 *hp_pt     = (TH1 *)fp->Get("h_cluster_pt_above10");

  TH1 *hj_cent   = (TH1 *)fj->Get("h_centrality");
  TH1 *hj_vtx    = (TH1 *)fj->Get("h_vertexz");
  TH1 *hj_cent10 = (TH1 *)fj->Get("h_centrality_with_cluster10");
  TH1 *hj_pt     = (TH1 *)fj->Get("h_cluster_pt_above10");

  TH2D *h2d_eta_phi_d = (TH2D *)fd->Get("h_cluster_eta_phi");
  TH2D *h2d_eta_phi_p = (TH2D *)fp->Get("h_cluster_eta_phi");
  TH2D *h2d_eta_phi_j = (TH2D *)fj->Get("h_cluster_eta_phi");

  TH2D *h2d_time_et_d = (TH2D *)fd->Get("h_cluster_time_et");
  TH2D *h2d_time_et_p = (TH2D *)fp->Get("h_cluster_time_et");
  TH2D *h2d_time_et_j = (TH2D *)fj->Get("h_cluster_time_et");

  TH2D *h2d_ieta_iphi_d = (TH2D *)fd->Get("h_cluster_ieta_iphi");
  TH2D *h2d_ieta_iphi_p = (TH2D *)fp->Get("h_cluster_ieta_iphi");
  TH2D *h2d_ieta_iphi_j = (TH2D *)fj->Get("h_cluster_ieta_iphi");

  TH1 *hd_scaler_per_run = (TH1 *)fd->Get("h_scaler_per_run");
  TH1 *hd_clusters_per_run = (TH1 *)fd->Get("h_clusters_per_run");
  TH1 *hd_avg_totalEMCal_energy_per_run = (TH1 *)fd->Get("h_avg_totalEMCal_energy_per_run");

  TH1 *hd_jet_pt = (TH1 *)fd->Get("h_jet_pt");
  TH1 *hp_jet_pt = (TH1 *)fp->Get("h_jet_pt");
  TH1 *hj_jet_pt = (TH1 *)fj->Get("h_jet_pt");
  TH2D *h2d_jet_eta_phi_d = (TH2D *)fd->Get("h_jet_eta_phi");
  TH2D *h2d_jet_eta_phi_p = (TH2D *)fp->Get("h_jet_eta_phi");
  TH2D *h2d_jet_eta_phi_j = (TH2D *)fj->Get("h_jet_eta_phi");

  if (!hd_cent || !hp_cent || !hj_cent) { std::cerr << "Missing h_centrality in one or more files." << std::endl; return; }
  if (!hd_vtx || !hp_vtx || !hj_vtx) { std::cerr << "Missing h_vertexz in one or more files." << std::endl; return; }
  if (!hd_cent10 || !hp_cent10 || !hj_cent10) { std::cerr << "Missing h_centrality_with_cluster10 in one or more files." << std::endl; return; }
  if (!hd_pt || !hp_pt || !hj_pt) { std::cerr << "Missing h_cluster_pt_above10 in one or more files." << std::endl; return; }
  if (!h2d_eta_phi_d || !h2d_eta_phi_p || !h2d_eta_phi_j) { std::cerr << "Missing h_cluster_eta_phi in one or more files." << std::endl; return; }
  if (!h2d_time_et_d || !h2d_time_et_p || !h2d_time_et_j) { std::cerr << "Missing h_cluster_time_et in one or more files." << std::endl; return; }
  if (!h2d_ieta_iphi_d || !h2d_ieta_iphi_p || !h2d_ieta_iphi_j) { std::cerr << "Missing h_cluster_ieta_iphi in one or more files." << std::endl; return; }

  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(1);
  gSystem->mkdir(plot_outdir, true);

  auto draw_compare = [&](TH1 *hd, TH1 *hp, TH1 *hj, const char *title, const char *outname, bool logy = false) {
    TH1 *h1 = (TH1 *)hd->Clone();
    TH1 *h2 = (TH1 *)hp->Clone();
    TH1 *h3 = (TH1 *)hj->Clone();

    TCanvas *c = new TCanvas("c", title, 700, 600);
    if (logy) c->SetLogy();

    float max = h1->GetMaximum();
    if(!logy) h1->GetYaxis()->SetRangeUser(0, max*1.15);

    h1->SetMarkerStyle(20);
    h1->SetMarkerSize(1.0);
    h1->SetMarkerColor(kBlack);
    h2->SetLineColor(kBlue);
    h2->SetLineWidth(2);
    h3->SetLineColor(kRed);
    h3->SetLineWidth(2);

    double idata = h1->Integral();
    if (idata > 0 && h2->Integral() > 0) h2->Scale(idata / h2->Integral());
    if (idata > 0 && h3->Integral() > 0) h3->Scale(idata / h3->Integral());

    double ymax = std::max(h1->GetMaximum(), std::max(h2->GetMaximum(), h3->GetMaximum()));
    if (hd->GetNbinsX() > 0) ymax = std::max(ymax, h1->GetBinContent(h1->GetMaximumBin()) + h1->GetBinError(h1->GetMaximumBin()));
    h1->SetMaximum(ymax * 1.15);
    h1->Draw("E");
    h2->Draw("HIST SAME");
    h3->Draw("HIST SAME");
    h1->Draw("E SAME");

    myText      (0.65, 0.92, 1, "#bf{#it{sPHENIX}} Internal", 0.04);
    myMarkerText(0.65, 0.87, kBlack, 20, legend_label_data, 1, 0.04);
    myMarkerText(0.65, 0.82, kBlue, 33, legend_label_photon, 1, 0.04);
    myMarkerText(0.65, 0.77, kRed, 34, legend_label_jet, 1, 0.04);

    c->SaveAs(Form("%s/%s.pdf", plot_outdir, outname));
    delete c;
    delete h1;
    delete h2;
    delete h3;
  };

  hd_cent->Scale(1.0,"width");
  hp_cent->Scale(1.0,"width");
  hj_cent->Scale(1.0,"width");

  
  draw_compare(hd_cent, hp_cent, hj_cent, ";Centrality [%];Events / Bin Width", "centrality");
  draw_compare(hd_vtx, hp_vtx, hj_vtx, ";Vertex Z [cm];Events", "vertexz");
  draw_compare(hd_cent10, hp_cent10, hj_cent10, ";Centrality [%];Events", "centrality");
  draw_compare(hd_pt, hp_pt, hj_pt, ";Cluster E_{T} [GeV];Clusters", "cluster_pt", true);



  // TH2: cluster yield eta vs phi (COLZ), one canvas per data type
  auto draw_th2_colz = [&](TH2 *h2, const char *label, const char *outname, const char *sublabel = nullptr) {
    TCanvas *c = new TCanvas("c", "", 700, 600);
    h2->Draw("COLZ");
    gPad->SetRightMargin(0.15);
    //gPad->SetLogz();

    myText      (0.65, 0.92, 1, "#bf{#it{sPHENIX}} Internal", 0.04);
    myText      (0.65, 0.87, 1,label , 0.04);
    if (sublabel) myText(0.65, 0.82, 1, sublabel, 0.04);
    c->SaveAs(Form("%s/%s.pdf", plot_outdir, outname));
    delete c;
  };

  std::string pt10_label = "#it{p}_{T}^{clus} > 10 GeV";

  draw_th2_colz(h2d_eta_phi_d, legend_label_data, "cluster_eta_phi_data",pt10_label.c_str());
  draw_th2_colz(h2d_eta_phi_p, legend_label_photon, "cluster_eta_phi_photon",pt10_label.c_str());
  draw_th2_colz(h2d_eta_phi_j, legend_label_jet, "cluster_eta_phi_jet",pt10_label.c_str());

  // TH2: cluster yield ieta vs iphi (COLZ), one canvas per data type
  draw_th2_colz(h2d_ieta_iphi_d, legend_label_data, "cluster_ieta_iphi_data",pt10_label.c_str());
  draw_th2_colz(h2d_ieta_iphi_p, legend_label_photon, "cluster_ieta_iphi_photon",pt10_label.c_str());
  draw_th2_colz(h2d_ieta_iphi_j, legend_label_jet, "cluster_ieta_iphi_jet",pt10_label.c_str());

  // TH2: cluster time vs ET (COLZ), one canvas per data type
  draw_th2_colz(h2d_time_et_d, legend_label_data, "cluster_time_et_data");
  draw_th2_colz(h2d_time_et_p, legend_label_photon, "cluster_time_et_photon");
  draw_th2_colz(h2d_time_et_j, legend_label_jet, "cluster_time_et_jet");

  // 1D cluster time: project TH2 onto time axis, all three on same plot
  TH1D *hd_time = h2d_time_et_d->ProjectionY("hd_time", 1, -1);
  TH1D *hp_time = h2d_time_et_p->ProjectionY("hp_time", 1, -1);
  TH1D *hj_time = h2d_time_et_j->ProjectionY("hj_time", 1, -1);
  draw_compare(hd_time, hp_time, hj_time, ";Cluster time [ns];Clusters", "cluster_time");
  delete hd_time;
  delete hp_time;
  delete hj_time;



  ////////////////////////////////////////
  // jet histograms
  draw_compare(hd_jet_pt, hp_jet_pt, hj_jet_pt, ";Jet p_{T} [GeV];Jets", "jet_pt", true);
  draw_th2_colz(h2d_jet_eta_phi_d, legend_label_data, "jet_eta_phi_data");
  draw_th2_colz(h2d_jet_eta_phi_p, legend_label_photon, "jet_eta_phi_photon");
  draw_th2_colz(h2d_jet_eta_phi_j, legend_label_jet, "jet_eta_phi_jet");

  ////////////////////////////////////////
  // iso ET vs centrality
  /////////////////////////////////////////
  for (int ipt = 0; ; ++ipt) {
    TString hname = Form("h_iso_vs_cent_pt%d", ipt);
    TH2D *h2_d = (TH2D *)fd->Get(hname);
    TH2D *h2_p = (TH2D *)fp->Get(hname);
    TH2D *h2_j = (TH2D *)fj->Get(hname);

    if (!h2_d || !h2_p || !h2_j) break;  // stop when no further pT bin exists

    h2_d->RebinX(2);
    h2_p->RebinX(2);
    h2_j->RebinX(2);

    std::string pt_label = Form("p_{T} bin %d", ipt);
    if (ipt + 1 < (int)pT_bins.size()) {
      pt_label = Form("%.0f < #it{p}_{T}^{clus} < %.0f GeV", pT_bins[ipt], pT_bins[ipt + 1]);
    }
    
    // Draw TH2 iso vs centrality for each sample
    draw_th2_colz(h2_d, legend_label_data, Form("iso_vs_cent_pt%d_data", ipt), pt_label.c_str());
    draw_th2_colz(h2_p, legend_label_photon, Form("iso_vs_cent_pt%d_photon", ipt), pt_label.c_str());
    draw_th2_colz(h2_j, legend_label_jet, Form("iso_vs_cent_pt%d_jet", ipt), pt_label.c_str());

    // TProfile <iso ET> vs centrality: data + MC on same plot
    TProfile *p_d = h2_d->ProfileX(Form("%s_prof_d", hname.Data()), 1, -1, "s");
    TProfile *p_p = h2_p->ProfileX(Form("%s_prof_p", hname.Data()), 1, -1, "s");
    TProfile *p_j = h2_j->ProfileX(Form("%s_prof_j", hname.Data()), 1, -1, "s");


    TGraphErrors *g_p = new TGraphErrors(p_p);
    TGraphErrors *g_j = new TGraphErrors(p_j);
    for (int ip = 0; ip < g_p->GetN(); ++ip) {
       g_p->GetX()[ip]-=0.5;
       g_j->GetX()[ip]+=0.5;
    }

    TCanvas *c = new TCanvas("c", "", 700, 600);
    p_d->SetLineColor(kBlack);
    p_d->SetMarkerColor(kBlack);
    p_d->SetMarkerStyle(20);
    p_d->SetMarkerSize(1.0);
    p_d->GetXaxis()->SetTitle("Centrality [%]");
    p_d->GetYaxis()->SetTitle("<Iso E_{T}> [GeV]");
    p_d->GetYaxis()->SetRangeUser(-20,20);
    p_d->Draw("ex0");

    g_p->SetLineColor(kBlue);
    g_p->SetLineWidth(2);
    g_p->SetMarkerColor(kBlue);
    g_p->SetMarkerStyle(23);
    g_p->Draw("p");

    g_j->SetLineColor(kRed);
    g_j->SetLineWidth(2);
    g_j->SetMarkerColor(kRed);
    g_j->SetMarkerStyle(22);
    g_j->Draw("p");

    p_d->Fit("pol2");

    myText      (0.65, 0.92, 1, "#bf{#it{sPHENIX}} Internal", 0.04);
    myText      (0.65, 0.87, 1, pt_label.c_str(), 0.035);
    myMarkerText(0.65, 0.82, kBlack, 20, legend_label_data, 1, 0.04);
    myMarkerText(0.65, 0.77, kBlue, 23, legend_label_photon, 1, 0.04);
    myMarkerText(0.65, 0.72, kRed, 22, legend_label_jet, 1, 0.04);

    c->SaveAs(Form("%s/iso_prof_vs_cent_pt%d.pdf", plot_outdir, ipt));
    delete c;
    delete p_d;
    if (p_p) delete p_p;
    if (p_j) delete p_j;
  }

  ////////////////////////////////////////
  // Cluster yield vs centrality (5% bins), per analysis pT bin
  ////////////////////////////////////////
  for (int ipt = 0;; ++ipt) {
    TString hname = Form("h_cluster_yield_vs_cent_pt%d", ipt);
    TH1 *hd = (TH1 *)fd->Get(hname);
    TH1 *hp = (TH1 *)fp->Get(hname);
    TH1 *hj = (TH1 *)fj->Get(hname);
    if (!hd || !hp || !hj)
      break;

    draw_compare(hd, hp, hj, ";Centrality [%];Clusters", Form("cluster_yield_vs_cent_pt%d", ipt));
  }

  ////////////////////////////////////////
  // run-by-run histograms
  // Data only: per-run scaler and cluster count histograms
  auto draw_data_only_1d = [&](TH1 *h, const char *outname, bool logy = false) {
    if (!h) return;
    TH1 *h1 = (TH1 *)h->Clone();
    TCanvas *c = new TCanvas("c", "", 1500, 600);
    if (logy) c->SetLogy();
    h1->SetMarkerStyle(20);
    h1->SetMarkerSize(1.0);
    h1->SetMarkerColor(kBlack);

    // Set y-range to [min, max] with 10% headroom
    double ymin = 0;
    double ymax = 0;
    bool first = true;
    for (int ib = 1; ib <= h1->GetNbinsX(); ++ib) {
      double v = h1->GetBinContent(ib);
      if (v <= 0) continue;
      if (first) {
        ymin = ymax = v;
        first = false;
      } else {
        if (v < ymin) ymin = v;
        if (v > ymax) ymax = v;
      }
    }
    if (!first) {
      double pad = 0.1 * (ymax - ymin);
      if (pad <= 0) pad = 0.1 * ymax;
      if (pad <= 0) pad = 1.0;
      if (!logy) {
        h1->SetMinimum(std::max(0.0, ymin - pad));
      } else {
        h1->SetMinimum(std::max(1e-6, ymin * 0.9));
      }
      h1->SetMaximum(ymax + pad);
    }

    h1->Draw("E");
    h1->GetXaxis()->SetNdivisions(505);
    
    myText      (0.65, 0.92, 1, "#bf{#it{sPHENIX}} Internal", 0.04);
    myMarkerText(0.65, 0.87, kBlack, 20, legend_label_data, 1, 0.04);
    c->SaveAs(Form("%s/%s.pdf", plot_outdir, outname));
    delete c;
    delete h1;
  };
  
  if (!hd_scaler_per_run || !hd_clusters_per_run || !hd_avg_totalEMCal_energy_per_run) {
    std::cerr << "Missing h_scaler_per_run, h_clusters_per_run, or h_avg_totalEMCal_energy_per_run in data file." << std::endl;
  } else {
  long total_scaler = 0;
  for (int ib = 1; ib <= hd_scaler_per_run->GetNbinsX(); ++ib) {
    double n = hd_scaler_per_run->GetBinContent(ib);
    hd_scaler_per_run->SetBinError(ib, std::sqrt(n));
    total_scaler += n;
  }
  std::cout << "Total scaler: " <<  total_scaler << "   lumi " <<  (float) total_scaler / (6.0 * 1e9) << " [nb^-1]" << std::endl;
  for (int ib = 1; ib <= hd_clusters_per_run->GetNbinsX(); ++ib) {
    double x = hd_clusters_per_run->GetBinContent(ib);
    if (x == 0) continue;
    float run = hd_clusters_per_run->GetBinLowEdge(ib);
    if (x > 200) 
    cout << "run " << run << " x " << x << endl;
    hd_clusters_per_run->SetBinError(ib, std::sqrt(x));
  }

  // Non-overlapping groups of kRunsPerGroup adjacent run bins (includes runs with 0 entries).
  const int kRunsPerGroup = 5;
  auto make_run_group_hist = [&](const TH1 *hPerRun, const char *name, const char *xtitle,
                                 const char *ytitle) -> TH1D * {
    const int nb = hPerRun->GetNbinsX();
    if (nb <= 0)
      return nullptr;
    const int ngrp = (nb + kRunsPerGroup - 1) / kRunsPerGroup;
    const double x0 = hPerRun->GetXaxis()->GetXmin();
    TH1D *hg =
        new TH1D(name, Form(";%s;%s", xtitle, ytitle), ngrp, x0, x0 + static_cast<double>(ngrp * kRunsPerGroup));
    hg->SetDirectory(nullptr);
    hg->Sumw2();
    for (int ig = 0; ig < ngrp; ++ig) {
      double sum = 0.;
      double err2 = 0.;
      for (int k = 0; k < kRunsPerGroup; ++k) {
        const int ib = ig * kRunsPerGroup + k + 1;
        if (ib > nb)
          break;
        sum += hPerRun->GetBinContent(ib);
        const double e = hPerRun->GetBinError(ib);
        err2 += e * e;
      }
      hg->SetBinContent(ig + 1, sum);
      hg->SetBinError(ig + 1, std::sqrt(err2));
    }
    return hg;
  };

  draw_data_only_1d(hd_scaler_per_run, "scaler_per_run_data", false);
  draw_data_only_1d(hd_clusters_per_run, "clusters_per_run_data", false);
  draw_data_only_1d(hd_avg_totalEMCal_energy_per_run, "avg_totalEMCal_energy_per_run_data", false);


  TH1F* hd_clusters_per_scaler = (TH1F*)hd_clusters_per_run->Clone("hd_clusters_per_scaler");
  hd_clusters_per_scaler->Divide(hd_scaler_per_run);

  draw_data_only_1d(hd_clusters_per_scaler, "clusters_per_scaler", false);

  TH1D *hd_scaler_5runs = make_run_group_hist(
      hd_scaler_per_run, "hd_scaler_per_5runs",
      "First run # in 5-run group", "Total MB scaler (sum of 5 runs)");
  TH1D *hd_clusters_5runs = make_run_group_hist(
      hd_clusters_per_run, "hd_clusters_per_5runs",
      "First run # in 5-run group", "N clusters (sum of 5 runs)");
  if (hd_scaler_5runs && hd_clusters_5runs) {
    draw_data_only_1d(hd_scaler_5runs, "scaler_per_5runs_data", false);
    draw_data_only_1d(hd_clusters_5runs, "clusters_per_5runs_data", false);
    TH1D *hd_clusters_per_scaler_5runs = (TH1D *)hd_clusters_5runs->Clone("hd_clusters_per_scaler_5runs");
    hd_clusters_per_scaler_5runs->Divide(hd_scaler_5runs);
    draw_data_only_1d(hd_clusters_per_scaler_5runs, "clusters_per_scaler_5runs", false);
    delete hd_clusters_per_scaler_5runs;
    delete hd_clusters_5runs;
    delete hd_scaler_5runs;
  }
  }

  
  fd->Close();
  fp->Close();
  fj->Close();

  std::cout << "Plots saved to " << plot_outdir << std::endl;
}
