#include <yaml-cpp/yaml.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLegend.h>
#include <TStyle.h>
#include <iostream>
#include <string>
#include <vector>

#include "plotcommon.h"

namespace {
TH1D *make_rms_vs_truth_pt(TH2D *h2, const char *name)
{
  if (!h2)
    return nullptr;
  const int nx = h2->GetNbinsX();
  TH1D *hRms = new TH1D(name, "", nx, h2->GetXaxis()->GetXmin(), h2->GetXaxis()->GetXmax());
  hRms->SetDirectory(nullptr);
  for (int ix = 1; ix <= nx; ++ix)
  {
    TH1D *py = h2->ProjectionY(Form("%s_py_%d", name, ix), ix, ix);
    py->SetDirectory(nullptr);
    hRms->SetBinContent(ix, py->GetRMS());
    hRms->SetBinError(ix, py->GetRMSError());
    delete py;
  }
  return hRms;
}
} // namespace

void plot_jet(std::string tune = "nom")
{
  init_plot();
  gStyle->SetOptStat(0);

  const std::string savePath = "./figs/jet";
  gSystem->mkdir(savePath.c_str(), true);

  gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
  const std::string yamlPath = Form("../histMakers/%s.yaml", tune.c_str());
  YAML::Node config = YAML::LoadFile(yamlPath.c_str());
  const std::vector<double> centrality_bins =
      config["analysis"]["centrality_bins"].as<std::vector<double>>();
  const int nCentBins = static_cast<int>(centrality_bins.size()) - 1;

  const int jet_cone_size = config["analysis"]["jet_cone_size"].as<int>(3);
  const bool use_jet_r04 = (jet_cone_size == 4);
  const float jet_eta = config["analysis"]["jet_eta"].as<float>();
  const float b2bjet_pT_min = config["analysis"]["b2bjet_pT_min"].as<float>();

  std::cout << "[plot_jet] jet_cone_size=" << jet_cone_size
            << " -> Anti-kT r0" << (use_jet_r04 ? "4" : "3")
            << ", jet_eta=" << jet_eta << ", b2bjet_pT_min=" << b2bjet_pT_min << " GeV" << std::endl;

  const std::string str_jet_cone =
      Form("#it{R}=0.%i, |#it{#eta}^{jet}| < %.1f",
           jet_cone_size, jet_eta);

  TFile *fin = TFile::Open("../histMakers/results/MC_efficiency_photon20_aa_showershape.root", "READ");
  if (!fin || fin->IsZombie())
  {
    std::cerr << "Could not open jet-response input file." << std::endl;
    return;
  }

  const int kMaxCent = 12;
  const int colors[kMaxCent] = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1,
                                kOrange + 7, kCyan + 1, kPink + 1, kYellow + 2, kAzure + 1,
                                kViolet + 1, kSpring + 5};

  for (int icent = 0; icent < nCentBins; ++icent)
  {
    TH2D *h2 = (TH2D *)fin->Get(Form("h_jet_pT_response_cent%d", icent));
    if (!h2)
    {
      std::cerr << "Missing histogram h_jet_pT_response_cent" << icent << std::endl;
      continue;
    }

    TH2D *h2_draw = (TH2D *)h2->Clone(Form("h2_draw_cent%d", icent));
    h2_draw->SetDirectory(nullptr);
    h2_draw->SetTitle(Form("Jet response %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]));
    h2_draw->GetXaxis()->SetTitle("Truth jet #it{p}_{T} [GeV]");
    h2_draw->GetYaxis()->SetTitle("#it{p}_{T}^{reco} / #it{p}_{T}^{truth}");

    TProfile *pfx = h2_draw->ProfileX(Form("pfx_jet_response_cent%d", icent));
    pfx->SetDirectory(nullptr);
    pfx->SetLineColor(kBlack);
    pfx->SetMarkerColor(kBlack);
    pfx->SetMarkerStyle(20);
    pfx->SetMarkerSize(0.9);
    pfx->SetLineWidth(2);

    TCanvas *c = new TCanvas(Form("c_jet_response_cent%d", icent), "", 760, 620);
    c->SetRightMargin(0.14);
    h2_draw->Draw("COLZ");
    pfx->Draw("EP SAME");

    TLine line_unity(5.0, 1.0, 50.0, 1.0);
    line_unity.SetLineStyle(7);
    line_unity.SetLineWidth(2);
    line_unity.SetLineColor(kBlack);
    line_unity.Draw("SAME");

    myText(0.18, 0.93, 1, strleg1.c_str(), 0.045);
    myText(0.18, 0.88, 1, strleg2.c_str(), 0.040);
    myText(0.18, 0.83, 1, Form("%.0f-%.0f%% centrality", centrality_bins[icent], centrality_bins[icent + 1]), 0.040);
    myText(0.18, 0.78, 1, str_jet_cone.c_str(), 0.038);

    c->SaveAs(Form("%s/h2_jet_response_with_profile_cent%d.pdf", savePath.c_str(), icent));
  }

  // Jet energy scale: mean response vs truth jet pT — one curve per centrality
  TCanvas *c_jes = new TCanvas("c_jes_vs_truthpt", "", 760, 620);
  TLegend *leg_jes = new TLegend(0.50, 0.15, 0.90, 0.45);
  leg_jes->SetBorderSize(0);
  leg_jes->SetFillStyle(0);

  bool first_jes = true;
  for (int icent = 0; icent < nCentBins; ++icent)
  {
    TH2D *h2 = (TH2D *)fin->Get(Form("h_jet_pT_response_cent%d", icent));
    if (!h2)
      continue;

    TProfile *pfx = h2->ProfileX(Form("pfx_jes_cent%d", icent), 1, -1, "");
    pfx->SetDirectory(nullptr);
    const int col = colors[icent % kMaxCent];
    pfx->SetLineColor(col);
    pfx->SetMarkerColor(col);
    pfx->SetMarkerStyle(20);
    pfx->SetMarkerSize(0.85);
    pfx->SetLineWidth(2);
    pfx->SetTitle(";Truth jet #it{p}_{T} [GeV];#LT#it{p}_{T}^{reco} / #it{p}_{T}^{truth}#GT");

    if (first_jes)
    {
      pfx->GetXaxis()->SetRangeUser(5.0, 50.0);
      pfx->GetYaxis()->SetRangeUser(0.0, 1.5);
      pfx->Draw("E1");
      first_jes = false;
    }
    else
    {
      pfx->Draw("E1 SAME");
    }

    leg_jes->AddEntry(pfx, Form("%.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), "lep");
  }

  TLine jes_unity(5.0, 1.0, 50.0, 1.0);
  jes_unity.SetLineStyle(7);
  jes_unity.SetLineWidth(2);
  jes_unity.SetLineColor(kGray + 2);
  jes_unity.Draw("SAME");

  myText(0.20, 0.93, 1, strleg1.c_str(), 0.045);
  myText(0.20, 0.88, 1, strleg2.c_str(), 0.040);
  myText(0.20, 0.83, 1, str_jet_cone.c_str(), 0.038);
  leg_jes->Draw();
  c_jes->SaveAs(Form("%s/h1_jet_energy_scale_vs_truthpt_by_cent.pdf", savePath.c_str()));



  // Resolution proxy: RMS of response vs truth jet pT — one curve per centrality
  TCanvas *c_jer = new TCanvas("c_jer_vs_truthpt", "", 760, 620);
  TLegend *leg_jer = new TLegend(0.50, 0.50, 0.90, 0.88);
  leg_jer->SetBorderSize(0);
  leg_jer->SetFillStyle(0);

  bool first_jer = true;
  for (int icent = 0; icent < nCentBins; ++icent)
  {
    TH2D *h2 = (TH2D *)fin->Get(Form("h_jet_pT_response_cent%d", icent));
    if (!h2)
      continue;

    TH1D *hRms = make_rms_vs_truth_pt(h2, Form("hjer_vs_truthpt_cent%d", icent));
    if (!hRms)
      continue;
    hRms->SetTitle(";Truth jet #it{p}_{T} [GeV];RMS(#it{p}_{T}^{reco} / #it{p}_{T}^{truth})");

    const int col = colors[icent % kMaxCent];
    hRms->SetLineColor(col);
    hRms->SetMarkerColor(col);
    hRms->SetMarkerStyle(21);
    hRms->SetMarkerSize(0.85);
    hRms->SetLineWidth(2);

    if (first_jer)
    {
      hRms->GetXaxis()->SetRangeUser(5.0, 50.0);
      hRms->GetYaxis()->SetRangeUser(0.0, 0.6);
      hRms->Draw("E1");
      first_jer = false;
    }
    else
    {
      hRms->Draw("E1 SAME");
    }

    leg_jer->AddEntry(hRms, Form("%.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), "lep");
  }

  myText(0.20, 0.93, 1, strleg1.c_str(), 0.045);
  myText(0.20, 0.88, 1, strleg2.c_str(), 0.040);
  myText(0.20, 0.83, 1, str_jet_cone.c_str(), 0.038);
  leg_jer->Draw();
  c_jer->SaveAs(Form("%s/h1_jet_resolution_vs_truthpt_by_cent.pdf", savePath.c_str()));

  fin->Close();
}
