#include <yaml-cpp/yaml.h>
#include <TSystem.h>
#include <TH2D.h>
#include <TFile.h>
#include <TLegend.h>
#include "plotcommon.h"

// jet20 MC: reco #it{E}_{T}^{iso} from h_iso_truth_reco vs h_background_iso_truth_reco (ProjectionY)
void plot_iso_jet20_sigbg(const std::string &tune = "nom", const char *rootpath = nullptr)
{
  init_plot();
  const std::string savePath = "./figs/iso";

  gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
  YAML::Node config = YAML::LoadFile("../histMakers/nom.yaml");
  std::vector<double> pT_bins = config["analysis"]["pT_bins"].as<std::vector<double>>();
  std::vector<double> centrality_bins = config["analysis"]["centrality_bins"].as<std::vector<double>>();
  const std::string var_type = config["output"]["var_type"].as<std::string>();
  std::string mcfile = "../histMakers/results/MC_efficiency_jet20_aa_showershape.root";////config["output"]["eff_outfile"].as<std::string>() + "_jet20_aa_" + var_type + ".root";
  if (rootpath && rootpath[0] != '\0')
    mcfile = rootpath;

  const int nPtBins = static_cast<int>(pT_bins.size()) - 1;
  const int nCentBins = static_cast<int>(centrality_bins.size()) - 1;

  TFile *fjet = TFile::Open(mcfile.c_str(), "READ");
  if (!fjet || fjet->IsZombie())
  {
    std::cerr << "[plot_iso_jet20_sigbg] Could not open jet20 MC file: " << mcfile << std::endl;
    return;
  }

  const int rbf = 10;
  const float xpos = 0.60, ypos = 0.92, dy = 0.05, fontsize = 0.042, fontsize1 = 0.045;

  for (int icent = 0; icent < nCentBins; ++icent)
  {
    for (int ipt = 0; ipt < nPtBins; ++ipt)
    {
      TH2D *h2sig = dynamic_cast<TH2D *>(fjet->Get(Form("h_iso_truth_reco_cent%d_pt%d", icent, ipt)));
      TH2D *h2bg = dynamic_cast<TH2D *>(fjet->Get(Form("h_background_iso_truth_reco_cent%d_pt%d", icent, ipt)));
      if (!h2sig || !h2bg)
      {
        std::cerr << "[plot_iso_jet20_sigbg] Missing TH2 for cent " << icent << " pT " << ipt << std::endl;
        continue;
      }

      TH1D *h_sig = h2sig->ProjectionY(Form("h_jet20_sig_cent%d_pt%d", icent, ipt));
      TH1D *h_bg = h2bg->ProjectionY(Form("h_jet20_bg_cent%d_pt%d", icent, ipt));
      if (!h_sig || !h_bg)
        continue;
      h_sig->SetDirectory(nullptr);
      h_bg->SetDirectory(nullptr);

      h_sig->Scale(1., "width");
      h_bg->Scale(1., "width");
      h_sig->Rebin(rbf);
      h_bg->Rebin(rbf);

      const double int_sig = h_sig->Integral("width");
      const double int_bg = h_bg->Integral("width");
      if (int_sig <= 0.0 || int_bg <= 0.0)
      {
        std::cerr << "[plot_iso_jet20_sigbg] Empty histogram cent " << icent << " pT " << ipt << std::endl;
        delete h_sig;
        delete h_bg;
        continue;
      }
      h_sig->Scale(1.0 / int_sig);
      h_bg->Scale(1.0 / int_bg);

      h_sig->SetLineColor(kBlue + 1);
      h_sig->SetLineWidth(2);
      h_sig->SetFillStyle(0);
      h_bg->SetLineColor(kRed + 1);
      h_bg->SetLineWidth(2);
      h_bg->SetFillStyle(0);

      TCanvas *c = new TCanvas(Form("c_jet20_sigbg_cent%d_pt%d", icent, ipt), "", 600, 560);
      h_sig->Draw("hist");
      h_sig->GetXaxis()->SetRangeUser(-10, 30);
      h_sig->SetYTitle("Self-normalized counts");
      h_sig->SetXTitle("#it{E}_{T}^{iso} [GeV]");
      h_sig->GetXaxis()->SetTitleOffset(1.2);
      h_bg->Draw("hist same");

      const std::string st_etbin =
          Form("%.0f < #it{E}_{T}^{#gamma} < %.0f GeV", pT_bins[ipt], pT_bins[ipt + 1]);
      const std::string st_centbin =
          Form("%.0f-%.0f%%", centrality_bins[icent], centrality_bins[icent + 1]);

      myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize1);
      myText(xpos, ypos - 1 * dy, 1, strleg2.c_str(), fontsize);
      myText(xpos, ypos - 2 * dy, 1, st_etbin.c_str(), fontsize);
      myText(xpos, ypos - 3 * dy, 1, st_centbin.c_str(), fontsize);
      myText(xpos, ypos - 4 * dy, 1, "jet20 MC", fontsize);

      TLegend *leg = new TLegend(0.50, 0.55, 0.90, 0.72);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->AddEntry(h_sig, "Signal (truth #gamma)", "l");
      leg->AddEntry(h_bg, "Background (jet)", "l");
      leg->Draw("same");

      c->SaveAs(Form("%s/h1D_iso_jet20_sigbg_%s_cent%d_pt%d.pdf", savePath.c_str(), tune.c_str(), icent, ipt));

      delete h_sig;
      delete h_bg;
      delete c;
    }
  }
  fjet->Close();
  delete fjet;
}

void plot_iso(string tune = "nom", int lowbin = 0, int highbin = 1)
{
  (void)lowbin;
  (void)highbin;
  init_plot();
  string savePath = "./figs/iso";

  gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
  YAML::Node config = YAML::LoadFile("../histMakers/nom.yaml");
  std::vector<double> pT_bins = config["analysis"]["pT_bins"].as<std::vector<double>>();
  std::vector<double> centrality_bins = config["analysis"]["centrality_bins"].as<std::vector<double>>();

  const int nPtBins = static_cast<int>(pT_bins.size()) - 1;
  const int nCentBins = static_cast<int>(centrality_bins.size()) - 1;


  TFile *fdata = TFile::Open("../histMakers/results/data_histo_showershape.root", "READ");
  TFile *fmc = TFile::Open("../histMakers/results/MC_efficiency_photon20_aa_showershape.root", "READ");
  if (!fdata || fdata->IsZombie())
  {
    std::cerr << "Could not open data file." << std::endl;
    return;
  }
  if (!fmc || fmc->IsZombie())
  {
    std::cerr << "Could not open signal MC file." << std::endl;
    return;
  }

  for (int icent = 0; icent < nCentBins; ++icent)
  {
    for (int ipt = 0; ipt < nPtBins; ++ipt)
    {
      TH1D *h_tight_isoET = (TH1D *)fdata->Get(Form("h_tight_isoET_cent%d_pt%d", icent, ipt));
      TH1D *h_tight_isoET_mcSig = (TH1D *)fmc->Get(Form("h_tight_isoET_cent%d_pt%d", icent, ipt));
      TH1D *h_nontight_isoET = (TH1D *)fdata->Get(Form("h_nontight_isoET_cent%d_pt%d", icent, ipt));
      if (!h_tight_isoET || !h_tight_isoET_mcSig || !h_nontight_isoET)
      {
        std::cerr << "Missing iso histogram for cent " << icent << ", pT bin " << ipt << std::endl;
        continue;
      }

      h_tight_isoET = (TH1D *)h_tight_isoET->Clone(Form("h_tight_isoET_cent%d_pt%d", icent, ipt));
      h_tight_isoET_mcSig = (TH1D *)h_tight_isoET_mcSig->Clone(Form("h_tight_isoET_mcSig_cent%d_pt%d", icent, ipt));
      h_nontight_isoET = (TH1D *)h_nontight_isoET->Clone(Form("h_nontight_isoET_cent%d_pt%d", icent, ipt));

    h_tight_isoET->Scale(1., "width");
    h_tight_isoET_mcSig->Scale(1., "width");
    h_nontight_isoET->Scale(1., "width");

    const int rbf = 20;
    h_tight_isoET->Rebin(rbf);
    h_tight_isoET_mcSig->Rebin(rbf);
    h_nontight_isoET->Rebin(rbf);


    float meanIsoET = h_tight_isoET->GetMean();
    float rmsIsoET = h_tight_isoET->GetRMS();
    const float ptcutnorm = meanIsoET + 1.75* rmsIsoET;
    
    const int ptB = h_tight_isoET->FindBin(ptcutnorm);
    const int ptBs = h_tight_isoET->GetNbinsX();
    const float normTight = h_tight_isoET->Integral(ptB, ptBs);
    const float normNontight = h_nontight_isoET->Integral(ptB, ptBs);
    if (normNontight > 0)
      h_nontight_isoET->Scale(normTight / normNontight);
    const float dataSig = h_tight_isoET->Integral() - h_nontight_isoET->Integral();
    const float mcInt = h_tight_isoET_mcSig->Integral();
    if (mcInt > 0)
      h_tight_isoET_mcSig->Scale(dataSig / mcInt);

      TCanvas *c1 = new TCanvas(Form("c1_cent%d_pt%d", icent, ipt), Form("c1_cent%d_pt%d", icent, ipt), 600, 560);
      h_tight_isoET->Draw("ex0");
      float ymax = h_tight_isoET->GetMaximum();
      h_tight_isoET->GetXaxis()->SetRangeUser(-20, 50);
      h_tight_isoET->GetYaxis()->SetRangeUser(0, ymax * 1.3);

      h_tight_isoET_mcSig->SetFillColorAlpha(kBlue, 0.3);
      h_tight_isoET_mcSig->SetFillStyle(1001);
      h_tight_isoET_mcSig->SetLineColor(kBlue - 1);
      h_tight_isoET_mcSig->SetLineWidth(2);

      h_nontight_isoET->SetFillColorAlpha(kRed, 0.3);
      h_nontight_isoET->SetFillStyle(1001);
      h_nontight_isoET->SetLineColor(kRed - 1);
      h_nontight_isoET->SetLineWidth(2);

      THStack *hs = new THStack(Form("hs_cent%d_pt%d", icent, ipt), "");
      hs->Add(h_nontight_isoET);
      hs->Add(h_tight_isoET_mcSig);
      hs->Draw("hist same");
      h_tight_isoET->Draw("same ex0");
      h_tight_isoET->Draw("same axis");

      h_tight_isoET->SetYTitle("Counts / Bin Width");
      h_tight_isoET->SetXTitle("#it{E}_{T}^{iso} [GeV]");
      h_tight_isoET->GetXaxis()->SetTitleOffset(1.2);

      string st_etbin = Form("%.0f < #it{E}_{T}^{#gamma} < %.0f GeV", pT_bins[ipt], pT_bins[ipt + 1]);
      string st_centbin = Form("%.0f-%.0f%% ", centrality_bins[icent], centrality_bins[icent + 1]);
      const float xpos = 0.6, ypos = 0.875, dy = 0.054, dy1 = 0.06, fontsize = 0.046, fontsize1 = 0.048;
      myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize1);
      myText(xpos, ypos - 1 * dy, 1, strleg2.c_str(), fontsize);
      myText(xpos, ypos - 2 * dy, 1, st_etbin.c_str(), fontsize);
      myText(xpos, ypos - 3 * dy, 1, st_centbin.c_str(), fontsize);

      TLegend *l1 = new TLegend(0.51, ypos - (3 + 4) * dy1, 0.9, ypos - 4 * dy + 0.03);
      l1->AddEntry(h_tight_isoET, "Data (Tight)", "pe");
      l1->AddEntry(h_nontight_isoET, "Data (Non-tight)", "f");
      l1->AddEntry(h_tight_isoET_mcSig, "Signal MC (Tight)", "f");
      l1->SetBorderSize(0);
      l1->SetFillStyle(0);
      l1->Draw("same");

      c1->SaveAs(Form("%s/h1D_iso_%s_cent%d_pt%d.pdf", savePath.c_str(), tune.c_str(), icent, ipt));





      TH1D *h_tight_isoET_shape = (TH1D *)h_tight_isoET->Clone(Form("h_tight_isoET_shape_cent%d_pt%d", icent, ipt));
      TH1D *h_nontight_isoET_shape = (TH1D *)h_nontight_isoET->Clone(Form("h_nontight_isoET_shape_cent%d_pt%d", icent, ipt));
      TH1D *h_tight_isoET_mcSig_shape = (TH1D *)h_tight_isoET_mcSig->Clone(Form("h_tight_isoET_mcSig_shape_cent%d_pt%d", icent, ipt));

      const double int_tight = h_tight_isoET_shape->Integral("width");
      if (int_tight > 0)
        h_tight_isoET_shape->Scale(1.0 / int_tight);
      const double int_nt = h_nontight_isoET_shape->Integral("width");
      if (int_nt > 0)
        h_nontight_isoET_shape->Scale(1.0 / int_nt);
      const double int_mc = h_tight_isoET_mcSig_shape->Integral("width");
      if (int_mc > 0)
        h_tight_isoET_mcSig_shape->Scale(1.0 / int_mc);

      h_tight_isoET_shape->SetMarkerStyle(20);
      h_tight_isoET_shape->SetMarkerColor(kBlack);
      h_tight_isoET_shape->SetLineColor(kBlack);

      h_nontight_isoET_shape->SetLineColor(kRed + 1);
      h_nontight_isoET_shape->SetLineWidth(2);
      h_nontight_isoET_shape->SetFillStyle(0);

      h_tight_isoET_mcSig_shape->SetLineColor(kBlue + 1);
      h_tight_isoET_mcSig_shape->SetLineWidth(2);
      h_tight_isoET_mcSig_shape->SetFillStyle(0);

      TCanvas *c2 = new TCanvas(Form("c2_cent%d_pt%d", icent, ipt), Form("c2_cent%d_pt%d", icent, ipt), 600, 560);
      h_tight_isoET_shape->Draw("ex0");
      h_tight_isoET_shape->GetXaxis()->SetRangeUser(-50, 50);
      h_tight_isoET_shape->SetYTitle("Self-normalized counts");
      h_tight_isoET_shape->SetXTitle("#it{E}_{T}^{iso} [GeV]");
      h_tight_isoET_shape->GetXaxis()->SetTitleOffset(1.2);
      h_nontight_isoET_shape->Draw("hist same");
      h_tight_isoET_mcSig_shape->Draw("hist same");
      h_tight_isoET_shape->Draw("same ex0");

      myText      (0.60, 0.92, 1, strleg1.c_str(), fontsize1);
      myText      (0.60, 0.87, 1, strleg2.c_str(), fontsize);
      myText      (0.60, 0.82, 1, st_etbin.c_str(), fontsize);
      myText      (0.80, 0.82, 1, st_centbin.c_str(), fontsize);
      myMarkerText(0.60, 0.77, kBlack, 20, "Data (Tight)", 1, 0.04);
      myMarkerText(0.60, 0.72, kRed + 1, 21, "Data (Non-tight)", 1, 0.04);
      myMarkerText(0.60, 0.67, kBlue + 1, 22, "Signal MC (Tight)", 1, 0.04);

      c2->SaveAs(Form("%s/h1D_iso_shape_%s_cent%d_pt%d.pdf", savePath.c_str(), tune.c_str(), icent, ipt));
    }
  }

  plot_iso_jet20_sigbg(tune, nullptr);
}

