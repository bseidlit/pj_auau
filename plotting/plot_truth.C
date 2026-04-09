#include <yaml-cpp/yaml.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>

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

    TCanvas *c = new TCanvas(Form("c_truth_pt_cent%d", icent), Form("c_truth_pt_cent%d", icent), 700, 560);
    h20->Draw("hist");
    h20->SetXTitle("#it{p}_{T}^{#gamma,truth} [GeV]");
    h20->SetYTitle("Counts / Bin width");
    h20->GetXaxis()->SetTitleOffset(1.15);
    h10->Draw("hist same");

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

    c->SaveAs(Form("%s/h1D_truth_pT_photon10_vs_photon20_cent%d_%s.pdf",
                   savePath.c_str(), icent, tune.c_str()));

    delete leg;
    delete c;
    delete h10;
    delete h20;
  }

  f10->Close();
  f20->Close();
  delete f10;
  delete f20;
}
