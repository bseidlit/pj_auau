#include"BlairUtils.C"
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <yaml-cpp/yaml.h>
#include <TF1.h>

const char *fdata = "../histMakers/results/evt_characterdata_aa.root";
const char *fphoton = "../histMakers/results/evt_characterphoton20_aa.root";
const char *fjet = "../histMakers/results/evt_characterjet20_aa.root";
const char *plot_outdir = "figs/";
const char *reweight_outfile = "output/centrality_reweighting.root";

const char *legend_label_data = "Data";
const char *legend_label_photon = "Pythia #gamma Overlay";
const char *legend_label_jet = "Pythia jet Overlay";
const char *config_yaml_path = "../histMakers/nom.yaml";


void reweightingDer()
{
  SetsPhenixStyle();
  gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
  gStyle->SetOptStat(0);

  std::vector<float> pT_bins;
  std::vector<float> centrality_bins;
  try {
    YAML::Node configYaml = YAML::LoadFile(config_yaml_path);
    if (configYaml["analysis"] && configYaml["analysis"]["pT_bins"]) {
      pT_bins = configYaml["analysis"]["pT_bins"].as<std::vector<float>>();
    }
    if (configYaml["analysis"] && configYaml["analysis"]["centrality_bins"]) {
      centrality_bins = configYaml["analysis"]["centrality_bins"].as<std::vector<float>>();
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
  gSystem->mkdir(plot_outdir, true);
  gSystem->mkdir("output", true);
  TFile *fout = TFile::Open(reweight_outfile, "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "Cannot open output ROOT file: " << reweight_outfile << std::endl;
    return;
  }
    
  ////////////////////////////////////////
  // Cluster yield vs centrality (5% bins), per analysis pT bin
  ////////////////////////////////////////

  TF1 *nom_cent_rw = nullptr;
  TH1 *nom_cent_rw_hist = nullptr;

  for (int ipt = 0;; ++ipt) {
    TString hname = Form("h_cluster_yield_vs_cent_pt%d", ipt);
    TH1 *hd = (TH1 *)fd->Get(hname);
    TH1 *hp = (TH1 *)fp->Get(hname);
    TH1 *hj = (TH1 *)fj->Get(hname);
    if (!hd || !hp || !hj)
      break;

    TH1D *h_data = (TH1D *)hd->Clone(Form("h_data_cent_pt%d", ipt));
    TH1D *h_photon = (TH1D *)hp->Clone(Form("h_photon_cent_pt%d", ipt));
    TH1D *h_jet = (TH1D *)hj->Clone(Form("h_jet_cent_pt%d", ipt));
    h_data->SetDirectory(nullptr);
    h_photon->SetDirectory(nullptr);
    h_jet->SetDirectory(nullptr);

    // Self-normalize each centrality spectrum before taking the ratio.
    const double i_data = h_data->Integral();
    const double i_photon = h_photon->Integral();
    const double i_jet = h_jet->Integral();
    if (i_data <= 0 || i_photon <= 0 || i_jet <= 0) {
      std::cerr << "Skipping pT bin " << ipt << " due to non-positive integral(s)." << std::endl;
      delete h_data;
      delete h_photon;
      delete h_jet;
      continue;
    }
    h_data->Scale(1.0 / i_data);
    h_photon->Scale(1.0 / i_photon);
    h_jet->Scale(1.0 / i_jet);

    TH1D *h_ratio_photon = (TH1D *)h_data->Clone(Form("h_ratio_data_over_photon_cent_pt%d", ipt));
    TH1D *h_ratio_jet = (TH1D *)h_data->Clone(Form("h_ratio_data_over_jet_cent_pt%d", ipt));
    h_ratio_photon->SetDirectory(nullptr);
    h_ratio_jet->SetDirectory(nullptr);
    h_ratio_photon->Divide(h_photon);
    h_ratio_jet->Divide(h_jet);

    TF1 *f_pol2_photon = new TF1(Form("f_cent_reweight_photon_pt%d", ipt), "pol2", 0.0, 100.0);
    TF1 *f_pol2_jet = new TF1(Form("f_cent_reweight_jet_pt%d", ipt), "pol2", 0.0, 100.0);
    h_ratio_photon->Fit(f_pol2_photon, "Q0R");
    h_ratio_jet->Fit(f_pol2_jet, "Q0R+");

    if (ipt == 0) {
      nom_cent_rw = dynamic_cast<TF1 *>(f_pol2_photon->Clone("nom_cent_rw"));
      nom_cent_rw_hist = dynamic_cast<TH1 *>(h_ratio_photon->Clone("nom_cent_rw_hist"));
      if (nom_cent_rw) {
        nom_cent_rw->SetTitle(
            "Nominal centrality reweight (photon MC, p_{T} bin 0);Centrality [%];w");
      }
    }

    h_ratio_photon->SetMarkerStyle(20);
    h_ratio_photon->SetMarkerSize(1.0);
    h_ratio_photon->SetMarkerColor(kBlue + 1);
    h_ratio_photon->SetLineColor(kBlue + 1);
    h_ratio_jet->SetMarkerStyle(21);
    h_ratio_jet->SetMarkerSize(1.0);
    h_ratio_jet->SetMarkerColor(kRed + 1);
    h_ratio_jet->SetLineColor(kRed + 1);
    f_pol2_photon->SetLineColor(kBlue + 1);
    f_pol2_photon->SetLineWidth(2);
    f_pol2_jet->SetLineColor(kRed + 1);
    f_pol2_jet->SetLineWidth(2);
    f_pol2_jet->SetLineStyle(2);

    double ymax = std::max(h_ratio_photon->GetMaximum(), h_ratio_jet->GetMaximum());
    double ymin = std::min(h_ratio_photon->GetMinimum(), h_ratio_jet->GetMinimum());
    if (!std::isfinite(ymin)) ymin = 0.0;
    if (!std::isfinite(ymax) || ymax <= 0.0) ymax = 2.0;
    h_ratio_photon->SetMinimum(std::max(0.0, ymin * 0.9));
    h_ratio_photon->SetMaximum(ymax * 1.2);
    h_ratio_photon->GetXaxis()->SetTitle("Centrality [%]");
    h_ratio_photon->GetYaxis()->SetTitle("Data / Sim (self-normalized)");

    TCanvas *c = new TCanvas(Form("c_reweight_cent_pt%d", ipt), "", 750, 620);
    h_ratio_photon->Draw("E");
    h_ratio_jet->Draw("E SAME");
    f_pol2_photon->Draw("SAME");
    f_pol2_jet->Draw("SAME");

    TLine *line1 = new TLine(0.0, 1.0, 100.0, 1.0);
    line1->SetLineColor(kGray + 2);
    line1->SetLineStyle(7);
    line1->Draw("SAME");

    std::string pt_label = Form("p_{T} bin %d", ipt);
    if (ipt + 1 < (int)pT_bins.size()) {
      pt_label = Form("%.0f < #it{E}_{T}^{clus} < %.0f GeV", pT_bins[ipt], pT_bins[ipt + 1]);
    }
    myText(0.60, 0.92, 1, "#bf{#it{sPHENIX}} Internal", 0.04);
    myText(0.60, 0.87, 1, pt_label.c_str(), 0.035);
    myMarkerText(0.60, 0.82, kBlue + 1, 20, "Data / #gamma MC", 1, 0.035);
    myMarkerText(0.60, 0.77, kRed + 1, 21, "Data / Jet MC", 1, 0.035);
    myText(0.60, 0.72, kBlue + 1, "pol2 fit (Data / #gamma MC)", 0.032);
    myText(0.60, 0.68, kRed + 1, "pol2 fit (Data / Jet MC)", 0.032);

    c->SaveAs(Form("%s/cent_reweight_ratio_fit_pt%d.pdf", plot_outdir, ipt));
    delete line1;
    delete c;

    fout->cd();
    h_ratio_photon->Write();
    h_ratio_jet->Write();
    f_pol2_photon->Write();
    f_pol2_jet->Write();
    nom_cent_rw->Write();
    nom_cent_rw_hist->Write();
    
    delete h_data;
    delete h_photon;
    delete h_jet;
    delete h_ratio_photon;
    delete h_ratio_jet;
    delete f_pol2_photon;
    delete f_pol2_jet;

  }

  fout->Close();
  fd->Close();
  fp->Close();
  fj->Close();
  std::cout << "Wrote reweight functions to " << reweight_outfile << std::endl;

}
