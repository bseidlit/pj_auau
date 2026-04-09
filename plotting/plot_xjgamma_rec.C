#include <yaml-cpp/yaml.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <iostream>
#include <string>
#include <vector>

#include "plotcommon.h"

namespace {
TH1D *project_xjgamma_ptbin(TH2D *h2, int yBin, const char *name)
{
  if (!h2)
    return nullptr;
  TH1D *h1 = h2->ProjectionX(name, yBin, yBin);
  if (h1)
    h1->SetDirectory(nullptr);
  return h1;
}

void style_mc_hist(TH1D *h, int color, float alpha = 0.35f)
{
  if (!h)
    return;
  h->SetFillColorAlpha(color, alpha);
  h->SetFillStyle(1001);
  h->SetLineColor(color - 1);
  h->SetLineWidth(2);
}

void style_data_hist(TH1D *h)
{
  if (!h)
    return;
  h->SetMarkerStyle(20);
  h->SetMarkerColor(kBlack);
  h->SetLineColor(kBlack);
}
} // namespace




void plot_xjgamma_rec(std::string tune = "nom")
{
  init_plot();
  std::string savePath = "./figs/xjgamma";
  gSystem->mkdir(savePath.c_str(), true);

  gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
  std::string yamlPath = Form("../histMakers/%s.yaml", tune.c_str());
  YAML::Node config = YAML::LoadFile(yamlPath.c_str());
  std::vector<double> pT_bins = config["analysis"]["pT_bins"].as<std::vector<double>>();
  std::vector<double> centrality_bins = config["analysis"]["centrality_bins"].as<std::vector<double>>();

  const int jet_cone_size = config["analysis"]["jet_cone_size"].as<int>(3);
  const bool use_jet_r04 = (jet_cone_size == 4);
  const float jet_eta = config["analysis"]["jet_eta"].as<float>();
  const float b2bjet_pT_min = config["analysis"]["b2bjet_pT_min"].as<float>();

  std::cout << "[plot_xjgamma_rec] jet_cone_size=" << jet_cone_size
            << " -> Anti-kT r0" << (use_jet_r04 ? "4" : "3")
            << ", jet_eta=" << jet_eta << ", b2bjet_pT_min=" << b2bjet_pT_min << " GeV" << std::endl;

  // Matches RecoEffCalculator_TTreeReader.C back-to-back jet cuts (reco jets)
  const std::string str_jet_cone =
      Form("#it{R}=0.%i, |#it{#eta}^{j}|<%.1f, #it{p}_{T}^{calib}>%.0f GeV",
           jet_cone_size, jet_eta, b2bjet_pT_min);

  const int nPtBins = static_cast<int>(pT_bins.size()) - 1;
  const int nCentBins = static_cast<int>(centrality_bins.size()) - 1;

  TFile *fdata = TFile::Open("../histMakers/results/data_histo_showershape.root", "READ");
  TFile *fsig = TFile::Open("../histMakers/results/MC_efficiency_photon20_aa_showershape.root", "READ");
  TFile *fbkg = TFile::Open("../histMakers/results/MC_efficiency_jet20_aa_showershape.root", "READ");



  for (int icent = 0; icent < nCentBins; ++icent)
  {
    const std::string str_cent =
        Form("%.0f-%.0f%% centrality", centrality_bins[icent], centrality_bins[icent + 1]);

    TH2D *h2_data = (TH2D *)fdata->Get(Form("h_tight_iso_xjgamma_cent%d", icent));
    TH2D *h2_sig  = (TH2D *)fsig->Get(Form("h_tight_iso_xjgamma_signal_cent%d", icent));
    TH2D *h2_bkg  = (TH2D *)fbkg->Get(Form("h_tight_iso_xjgamma_background_cent%d", icent));

    TH2D *h2_sig_reco = (TH2D *)fsig->Get(Form("h_tight_iso_xjgamma_signal_cent%d", icent));
    TH2D *h2_sig_truth = (TH2D *)fsig->Get(Form("h_tight_iso_truthjet_xjgamma_signal_cent%d", icent));
    TH2D *h2_sig_truthmatchreco =
        (TH2D *)fsig->Get(Form("h_tight_iso_truthmatchreco_xjgamma_signal_cent%d", icent));

    for (int ipt = 0; ipt < nPtBins; ++ipt)
    {
      const std::string str_photon_et =
          Form("%.0f < #it{E}_{T}^{#gamma} < %.0f GeV", pT_bins[ipt], pT_bins[ipt + 1]);

      TH1D *h_data = project_xjgamma_ptbin(h2_data, ipt + 1, Form("h_data_xj_tightiso_cent%d_pt%d", icent, ipt));
      TH1D *h_sig = project_xjgamma_ptbin(h2_sig, ipt + 1, Form("h_sig_xj_tightiso_cent%d_pt%d", icent, ipt));
      TH1D *h_bkg = project_xjgamma_ptbin(h2_bkg, ipt + 1, Form("h_bkg_xj_tightiso_cent%d_pt%d", icent, ipt));

      h_data->Scale(1., "width");
      h_sig->Scale(1., "width");
      h_bkg->Scale(1., "width");

      const double int_data = h_data->Integral("width");
      const double int_sig = h_sig->Integral("width");
      const double int_bkg = h_bkg->Integral("width");
      h_sig->Scale(int_data / int_sig);
      h_bkg->Scale(int_data / int_bkg);

      h_data->SetXTitle("#it{x}_{J#gamma} = #it{p}_{T}^{jet,reco} / #it{E}_{T}^{#gamma}");
      h_data->SetYTitle("Counts / Bin width");

      style_data_hist(h_data);
      style_mc_hist(h_sig, kBlue, 0.35f);
      style_mc_hist(h_bkg, kRed, 0.35f);

      TCanvas *c1 = new TCanvas(Form("c_xj_data_mc_cent%d_pt%d", icent, ipt), Form("c_xj_data_mc_cent%d_pt%d", icent, ipt), 700, 560);
      h_data->GetXaxis()->SetTitleOffset(1.15);
      h_data->Draw("e0");
      h_sig->Draw("hist same");
      h_bkg->Draw("hist same");
      h_data->Draw("same e0");
      h_data->Draw("same axis");

      const float xpos = 0.6f, ypos = 0.875f, dy = 0.054f, fontsize = 0.046f, fontsize1 = 0.048f;
      const float mx = 0.75f;
      const double tleg = 0.04;
      myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize1);
      myText(xpos, ypos - 1 * dy, 1, strleg2.c_str(), fontsize);
      myText(xpos, ypos - 2 * dy, 1, str_cent.c_str(), fontsize);
      myText(xpos, ypos - 3 * dy, 1, strlegphotonjet.c_str(), fontsize);
      myText(xpos, ypos - 4 * dy, 1, str_photon_et.c_str(), fontsize);
      myText(xpos, ypos - 5 * dy, 1, strdijet.c_str(), fontsize);
      myText(xpos, ypos - 6 * dy, 1, str_jet_cone.c_str(), fontsize);
      myMarkerText(mx, ypos - 7 * dy - 0.02f, kBlack, 20, "Data", 1, tleg);
      myMarkerText(mx, ypos - 7 * dy - 0.07f, kBlue + 1, 22, "Signal MC (#gamma), norm. to data", 1, tleg);
      myMarkerText(mx, ypos - 7 * dy - 0.12f, kRed + 1, 21, "Jet MC (bkg), norm. to data", 1, tleg);

      c1->SaveAs(Form("%s/h1D_xjgamma_tightiso_data_mc_cent%d_pt%d.pdf", savePath.c_str(), icent, ipt));




      // --- MC only: reco jet vs truth jet (signal), same selection ---
      TH1D *h_mc_reco = project_xjgamma_ptbin(h2_sig_reco, ipt + 1, Form("h_mc_reco_xj_cent%d_pt%d", icent, ipt));
      TH1D *h_mc_truth = project_xjgamma_ptbin(h2_sig_truth, ipt + 1, Form("h_mc_truth_xj_cent%d_pt%d", icent, ipt));
      TH1D *h_mc_truthmatchreco =
          project_xjgamma_ptbin(h2_sig_truthmatchreco, ipt + 1, Form("h_mc_truthmatchreco_xj_cent%d_pt%d", icent, ipt));

      h_mc_reco->Scale(1., "width");
      h_mc_truth->Scale(1., "width");
      h_mc_truthmatchreco->Scale(1., "width");
      //const double ir = h_mc_reco->Integral("width");
      //const double it = h_mc_truth->Integral("width");
      //h_mc_reco->Scale(1. / ir);
      //h_mc_truth->Scale(1. / it);

      h_mc_reco->SetLineColor(kBlue + 1);
      h_mc_reco->SetLineWidth(2);
      h_mc_reco->SetFillStyle(0);
      h_mc_truth->SetLineColor(kRed + 1);
      h_mc_truth->SetLineStyle(2);
      h_mc_truth->SetLineWidth(2);
      h_mc_truth->SetFillStyle(0);

      h_mc_truthmatchreco->SetLineColor(kGreen + 2);
      h_mc_truthmatchreco->SetLineStyle(1);
      h_mc_truthmatchreco->SetLineWidth(2);
      h_mc_truthmatchreco->SetFillStyle(0);

      h_mc_reco->SetXTitle("#it{x}_{J#gamma}");
      h_mc_reco->SetYTitle("Pairs / bin width");
      h_mc_truth->SetXTitle("#it{x}_{J#gamma}");
      h_mc_truth->SetYTitle("Pairs / bin width");

      float max_y = std::max(h_mc_reco->GetMaximum(), h_mc_truth->GetMaximum());
      if (h_mc_truthmatchreco)
        max_y = std::max(max_y, static_cast<float>(h_mc_truthmatchreco->GetMaximum()));
      h_mc_reco->GetYaxis()->SetRangeUser(0, max_y * 1.1);
      h_mc_truth->GetYaxis()->SetRangeUser(0, max_y * 1.1);
      h_mc_truthmatchreco->GetYaxis()->SetRangeUser(0, max_y * 1.1);

      TCanvas *c2 = new TCanvas(Form("c_xj_mc_truth_reco_cent%d_pt%d", icent, ipt), Form("c_xj_mc_truth_reco_cent%d_pt%d", icent, ipt), 700, 560);
      h_mc_truth->Draw("hist");
      h_mc_reco->Draw("hist same");
      h_mc_truthmatchreco->Draw("hist same");

      myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize1);
      myText(xpos, ypos - 1 * dy, 1, strleg2.c_str(), fontsize);
      myText(xpos, ypos - 2 * dy, 1, str_cent.c_str(), fontsize);
      myText(xpos, ypos - 3 * dy, 1, strSigMC.c_str(), fontsize);
      myText(xpos, ypos - 4 * dy, 1, str_photon_et.c_str(), fontsize);
      myText(xpos, ypos - 5 * dy, 1, str_jet_cone.c_str(), fontsize);
      myMarkerText(mx, ypos - 6 * dy - 0.02f, kBlue + 1, 22, "Reco jet", 1, tleg);
      myMarkerText(mx, ypos - 6 * dy - 0.07f, kRed + 1, 21, "Truth jet", 1, tleg);
      myMarkerText(mx, ypos - 6 * dy - 0.12f, kGreen + 2, 22, "Truth-matched reco jet (#Delta#it{R}<0.2)", 1, tleg);

      c2->SaveAs(Form("%s/h1D_xjgamma_mc_reco_vs_truth_cent%d_pt%d.pdf", savePath.c_str(), icent, ipt));
    
    }
  }

  fbkg->Close();
  fsig->Close();
  fdata->Close();
}
