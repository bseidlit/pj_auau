#include <yaml-cpp/yaml.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPad.h>
#include <TLine.h>

#include <string>
#include <vector>
#include <iostream>

#include "plotcommon.h"

void plot_truth(std::string tune = "nom")
{
  init_plot();
  const std::string savePath = "./figs/truth";
  gSystem->mkdir(savePath.c_str(), true);

  gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
  const std::string yamlPath = Form("../histMakers/%s.yaml", tune.c_str());
  YAML::Node config = YAML::LoadFile(yamlPath.c_str());

  std::vector<double> centrality_bins = config["analysis"]["centrality_bins"].as<std::vector<double>>();
  const std::string var_type = config["output"]["var_type"].as<std::string>();
  const int nCentBins = static_cast<int>(centrality_bins.size()) - 1;

  const std::string base_eff = config["output"]["eff_outfile"].as<std::string>();
  const std::string f10_name = base_eff + "_photon10_aa_" + var_type + ".root";
  const std::string f20_name = base_eff + "_photon20_aa_" + var_type + ".root";

  TFile *f10 = TFile::Open(f10_name.c_str(), "READ");
  TFile *f20 = TFile::Open(f20_name.c_str(), "READ");
  if (!f10 || f10->IsZombie())
  {
    std::cerr << "[plot_truth] Could not open photon10 file: " << f10_name << std::endl;
    return;
  }
  if (!f20 || f20->IsZombie())
  {
    std::cerr << "[plot_truth] Could not open photon20 file: " << f20_name << std::endl;
    return;
  }

  for (int icent = 0; icent < nCentBins; ++icent)
  {
    TH1D *h10_in = dynamic_cast<TH1D *>(f10->Get(Form("h_truth_pT_%d", icent)));
    TH1D *h20_in = dynamic_cast<TH1D *>(f20->Get(Form("h_truth_pT_%d", icent)));
    if (!h10_in || !h20_in)
    {
      std::cerr << "[plot_truth] Missing h_truth_pT_" << icent << " in photon10/photon20 files." << std::endl;
      continue;
    }

    TH1D *h10 = dynamic_cast<TH1D *>(h10_in->Clone(Form("h_truth_pT_photon10_cent%d", icent)));
    TH1D *h20 = dynamic_cast<TH1D *>(h20_in->Clone(Form("h_truth_pT_photon20_cent%d", icent)));
    if (!h10 || !h20)
      continue;
    h10->SetDirectory(nullptr);
    h20->SetDirectory(nullptr);

    h10->Scale(1., "width");
    h20->Scale(1., "width");

    h10->SetLineColor(kRed + 1);
    h10->SetLineWidth(2);
    h10->SetFillStyle(0);
    h20->SetLineColor(kBlue + 1);
    h20->SetLineWidth(2);
    h20->SetFillStyle(0);

    const std::string st_centbin =
        Form("%.0f-%.0f%% centrality", centrality_bins[icent], centrality_bins[icent + 1]);

    TH1D *h_ratio = dynamic_cast<TH1D *>(h10->Clone(Form("h_truth_pT_ratio_10_over_20_cent%d", icent)));
    if (!h_ratio)
      continue;
    h_ratio->SetDirectory(nullptr);
    h_ratio->Divide(h20);
    h_ratio->SetLineColor(kBlack);
    h_ratio->SetLineWidth(2);
    h_ratio->SetMarkerStyle(20);
    h_ratio->SetMarkerSize(0.6);
    h_ratio->SetFillStyle(0);
    h_ratio->SetTitle("");
    h_ratio->SetYTitle("#gamma10 / #gamma20");
    h_ratio->SetXTitle("#it{p}_{T}^{#gamma,truth} [GeV]");
    // Larger text in the small lower pad (sizes are fractions of pad height).
    const float ratioLabelSize = 0.065f;
    const float ratioTitleSize = 0.072f;
    h_ratio->GetXaxis()->SetLabelSize(ratioLabelSize);
    h_ratio->GetXaxis()->SetTitleSize(ratioTitleSize);
    h_ratio->GetXaxis()->SetTitleOffset(1.25);
    h_ratio->GetXaxis()->SetTickLength(0.03);
    h_ratio->GetYaxis()->SetLabelSize(ratioLabelSize);
    h_ratio->GetYaxis()->SetTitleSize(ratioTitleSize);
    h_ratio->GetYaxis()->SetTitleOffset(1.45);
    h_ratio->GetYaxis()->SetNdivisions(505);
    h_ratio->GetYaxis()->SetTickLength(0.02);
    h_ratio->GetYaxis()->SetRangeUser(0.5,1.5);

    TCanvas *c = new TCanvas(Form("c_truth_pt_cent%d", icent), Form("c_truth_pt_cent%d", icent), 700, 700);
    TPad *padUp = new TPad("padUp", "padUp", 0.0, 0.35, 1.0, 1.0);
    TPad *padDn = new TPad("padDn", "padDn", 0.0, 0.0, 1.0, 0.35);
    padUp->SetBottomMargin(0.02);
    padUp->SetLeftMargin(0.12);
    padUp->SetRightMargin(0.05);
    padDn->SetTopMargin(0.03);
    padDn->SetBottomMargin(0.32); // room for larger x-axis title/labels
    padDn->SetLeftMargin(0.12);
    padDn->SetRightMargin(0.05);
    padUp->Draw();
    padDn->Draw();

    padUp->cd();
    h10->Draw("hist");
    h10->SetXTitle("");
    h10->SetYTitle("Counts / Bin width");
    h10->GetXaxis()->SetLabelSize(0);
    h10->GetXaxis()->SetTitleSize(0);
    h10->GetYaxis()->SetTitleOffset(1.35);
    h20->Draw("hist same");
    gPad->SetLogy();

    const float xpos = 0.60f, ypos = 0.875f, dy = 0.054f, fontsize = 0.046f, fontsize1 = 0.048f;
    myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize1);
    myText(xpos, ypos - 1 * dy, 1, strleg2.c_str(), fontsize);
    myText(xpos, ypos - 2 * dy, 1, st_centbin.c_str(), fontsize);
    myText(xpos, ypos - 3 * dy, 1, strleg3.c_str(), fontsize);

    TLegend *leg = new TLegend(0.58, 0.62, 0.90, 0.76);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h20, "#gamma20 sample", "l");
    leg->AddEntry(h10, "#gamma10 sample", "l");
    leg->Draw("same");

    padDn->cd();
    h_ratio->Draw("e0");
    const double xmin = h_ratio->GetXaxis()->GetXmin();
    const double xmax = h_ratio->GetXaxis()->GetXmax();
    TLine lineOne(xmin, 1.0, xmax, 1.0);
    lineOne.SetLineColor(kGray + 2);
    lineOne.SetLineStyle(2);
    lineOne.SetLineWidth(1);
    lineOne.Draw("same");

    c->SaveAs(Form("%s/h1D_truth_pT_photon10_vs_photon20_cent%d_%s.pdf",
                   savePath.c_str(), icent, tune.c_str()));

    delete leg;
    delete c;
    delete h10;
    delete h20;
    delete h_ratio;
  }

  f10->Close();
  f20->Close();
  delete f10;
  delete f20;

  // Leading truth jet pT (RecoEffCalculator): filled for jet MC only — compare jet10 vs jet20
  const std::string fj10_name = base_eff + "_jet10_aa_" + var_type + ".root";
  const std::string fj20_name = base_eff + "_jet20_aa_" + var_type + ".root";
  TFile *fj10 = TFile::Open(fj10_name.c_str(), "READ");
  TFile *fj20 = TFile::Open(fj20_name.c_str(), "READ");
  if (!fj10 || fj10->IsZombie() || !fj20 || fj20->IsZombie())
  {
    std::cerr << "[plot_truth] Skip leading truth jet plot (open jet10/jet20 files): "
              << fj10_name << " / " << fj20_name << std::endl;
    if (fj10)
    {
      fj10->Close();
      delete fj10;
    }
    if (fj20)
    {
      fj20->Close();
      delete fj20;
    }
    return;
  }

  TH1 *hj10_in = dynamic_cast<TH1 *>(fj10->Get("h_max_truth_jet_pT"));
  TH1 *hj20_in = dynamic_cast<TH1 *>(fj20->Get("h_max_truth_jet_pT"));
  if (!hj10_in || !hj20_in)
  {
    std::cerr << "[plot_truth] Missing h_max_truth_jet_pT in jet10/jet20 files." << std::endl;
    fj10->Close();
    fj20->Close();
    delete fj10;
    delete fj20;
    return;
  }

  TH1 *hj10 = dynamic_cast<TH1 *>(hj10_in->Clone("h_max_truth_jet_pT_jet10"));
  TH1 *hj20 = dynamic_cast<TH1 *>(hj20_in->Clone("h_max_truth_jet_pT_jet20"));
  if (!hj10 || !hj20)
  {
    fj10->Close();
    fj20->Close();
    delete fj10;
    delete fj20;
    return;
  }
  hj10->SetDirectory(nullptr);
  hj20->SetDirectory(nullptr);

  //hj10->Scale(1., "width");
  //hj20->Scale(1., "width");

  hj10->SetLineColor(kRed + 1);
  hj10->SetLineWidth(2);
  hj10->SetFillStyle(0);
  hj20->SetLineColor(kBlue + 1);
  hj20->SetLineWidth(2);
  hj20->SetFillStyle(0);

  TH1 *h_ratio_j = dynamic_cast<TH1 *>(hj10->Clone("h_max_truth_jet_pT_ratio_10_over_20"));
  if (!h_ratio_j)
  {
    delete hj10;
    delete hj20;
    fj10->Close();
    fj20->Close();
    delete fj10;
    delete fj20;
    return;
  }
  h_ratio_j->SetDirectory(nullptr);
  h_ratio_j->Divide(hj20);
  h_ratio_j->SetLineColor(kBlack);
  h_ratio_j->SetLineWidth(2);
  h_ratio_j->SetMarkerStyle(20);
  h_ratio_j->SetMarkerSize(0.6);
  h_ratio_j->SetFillStyle(0);
  h_ratio_j->SetTitle("");
  h_ratio_j->SetYTitle("jet10 / jet20");
  h_ratio_j->SetXTitle("Leading truth jet #it{p}_{T} [GeV]");
  const float ratioLabelSizeJ = 0.065f;
  const float ratioTitleSizeJ = 0.072f;
  h_ratio_j->GetXaxis()->SetLabelSize(ratioLabelSizeJ);
  h_ratio_j->GetXaxis()->SetTitleSize(ratioTitleSizeJ);
  h_ratio_j->GetXaxis()->SetTitleOffset(1.25);
  h_ratio_j->GetXaxis()->SetTickLength(0.03);
  h_ratio_j->GetYaxis()->SetLabelSize(ratioLabelSizeJ);
  h_ratio_j->GetYaxis()->SetTitleSize(ratioTitleSizeJ);
  h_ratio_j->GetYaxis()->SetTitleOffset(1.45);
  h_ratio_j->GetYaxis()->SetNdivisions(505);
  h_ratio_j->GetYaxis()->SetTickLength(0.02);
  h_ratio_j->GetYaxis()->SetRangeUser(0.5, 1.5);

  TCanvas *cj = new TCanvas("c_max_truth_jet_pT", "c_max_truth_jet_pT", 700, 700);
  TPad *padUpJ = new TPad("padUpJ", "padUpJ", 0.0, 0.35, 1.0, 1.0);
  TPad *padDnJ = new TPad("padDnJ", "padDnJ", 0.0, 0.0, 1.0, 0.35);
  padUpJ->SetBottomMargin(0.02);
  padUpJ->SetLeftMargin(0.12);
  padUpJ->SetRightMargin(0.05);
  padDnJ->SetTopMargin(0.03);
  padDnJ->SetBottomMargin(0.32);
  padDnJ->SetLeftMargin(0.12);
  padDnJ->SetRightMargin(0.05);
  padUpJ->Draw();
  padDnJ->Draw();

  padUpJ->cd();
  hj10->Draw("hist");
  hj10->SetXTitle("p_{T}^{truth, leading jet}} [GeV]");
  hj10->GetYaxis()->SetRangeUser(10, 1e15);
  hj10->SetYTitle("Counts / Bin width");
  hj10->GetXaxis()->SetLabelSize(0);
  hj10->GetXaxis()->SetTitleSize(0);
  hj10->GetYaxis()->SetTitleOffset(1.35);
  hj20->Draw("hist same");
  gPad->SetLogy();

  const float xposJ = 0.60f, yposJ = 0.875f, dyJ = 0.054f, fontsizeJ = 0.046f, fontsize1J = 0.048f;
  myText(xposJ, yposJ - 0 * dyJ, 1, strleg1.c_str(), fontsize1J);
  myText(xposJ, yposJ - 1 * dyJ, 1, strleg2.c_str(), fontsizeJ);
  myText(xposJ, yposJ - 2 * dyJ, 1, "Leading truth jet (jet MC)", fontsizeJ);

  TLegend *legJ = new TLegend(0.58, 0.62, 0.90, 0.76);
  legJ->SetBorderSize(0);
  legJ->SetFillStyle(0);
  legJ->AddEntry(hj20, "jet20 sample", "l");
  legJ->AddEntry(hj10, "jet10 sample", "l");
  legJ->Draw("same");

  padDnJ->cd();
  h_ratio_j->Draw("e0");
  const double xminJ = h_ratio_j->GetXaxis()->GetXmin();
  const double xmaxJ = h_ratio_j->GetXaxis()->GetXmax();
  TLine lineOneJ(xminJ, 1.0, xmaxJ, 1.0);
  lineOneJ.SetLineColor(kGray + 2);
  lineOneJ.SetLineStyle(2);
  lineOneJ.SetLineWidth(1);
  lineOneJ.Draw("same");

  cj->SaveAs(Form("%s/h1D_max_truth_jet_pT_jet10_vs_jet20_%s.pdf", savePath.c_str(), tune.c_str()));

  delete legJ;
  delete cj;
  delete hj10;
  delete hj20;
  delete h_ratio_j;

  fj10->Close();
  fj20->Close();
  delete fj10;
  delete fj20;
}
