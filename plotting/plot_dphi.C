#include "BlairUtils.C"
#include <TSystem.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

// RecoEffCalculator_TTreeReader.C output (var_type from nom.yaml, e.g. showershape)
const char *fdata = "../histMakers/results/data_histo_showershape.root";
const char *fphoton = "../histMakers/results/MC_efficiency_photon12_aa_showershape.root";
const char *fjet = "../histMakers/results/MC_efficiency_jet20_aa_showershape.root";
const char *plot_outdir = "figs/dphi";
const char *config_yaml_path = "../histMakers/nom.yaml";

const char *legend_label_data = "Data";
const char *legend_label_photon = "Pythia #gamma Overlay";
const char *legend_label_jet = "Pythia jet Overlay";

// Histograms: h_dphi_clusterJets_tight_cent{icent}_pt{ipt} (tight clusters, RecoEff)
// Per-(tight) cluster: scale each Δφ spectrum by 1 / yield in the same (cent, pT) bin of
// h_tight_cluster_{icent} from the same file (RecoEff: tight cluster E_T spectrum).
//
// Event-level: h_totalEMCal_energy_tight_weight, h_totalIHCal_energy_tight_weight,
// h_totalOHCal_energy_tight_weight — total calorimeter energy for events with at least one tight
// cluster (weighted); shape-normalized overlays (data vs photon vs jet MC).
// h_centrality_tight_weight — same selection, weighted centrality [%] (shape-normalized overlay).
// h_vertexz — weighted primary vertex z [cm] for events in the RecoEff event loop (shape-normalized).

void plot_dphi()
{
  SetsPhenixStyle();
  gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");

  std::vector<float> pT_bins;
  std::vector<float> centrality_bins;
  int jet_cone_size = 3;
  float b2bjet_pT_min = 8.0f;
  try {
    YAML::Node configYaml = YAML::LoadFile(config_yaml_path);
    if (configYaml["analysis"] && configYaml["analysis"]["pT_bins"])
      pT_bins = configYaml["analysis"]["pT_bins"].as<std::vector<float>>();
    if (configYaml["analysis"] && configYaml["analysis"]["centrality_bins"])
      centrality_bins = configYaml["analysis"]["centrality_bins"].as<std::vector<float>>();
    if (configYaml["analysis"] && configYaml["analysis"]["jet_cone_size"])
      jet_cone_size = configYaml["analysis"]["jet_cone_size"].as<int>();
    if (configYaml["analysis"] && configYaml["analysis"]["b2bjet_pT_min"])
      b2bjet_pT_min = configYaml["analysis"]["b2bjet_pT_min"].as<float>();
    std::cout << "[plot_dphi] Loaded config: " << config_yaml_path
              << "  jet_cone_size=" << jet_cone_size << "  b2bjet_pT_min=" << b2bjet_pT_min << " GeV"
              << std::endl;
  } catch (const std::exception &e) {
    std::cerr << "[plot_dphi] Warning: failed to read " << config_yaml_path << " (" << e.what() << ")"
              << std::endl;
  }

  TFile *fd = TFile::Open(fdata, "READ");
  TFile *fp = TFile::Open(fphoton, "READ");
  TFile *fj = TFile::Open(fjet, "READ");

  if (!fd || fd->IsZombie()) {
    std::cerr << "[plot_dphi] Cannot open data file: " << fdata << std::endl;
    return;
  }
  if (!fp || fp->IsZombie()) {
    std::cerr << "[plot_dphi] Cannot open photon MC file: " << fphoton << std::endl;
    return;
  }
  if (!fj || fj->IsZombie()) {
    std::cerr << "[plot_dphi] Cannot open jet MC file: " << fjet << std::endl;
    return;
  }

  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(1);
  gSystem->mkdir(plot_outdir, true);

  if (pT_bins.size() < 2 || centrality_bins.size() < 2) {
    std::cerr << "[plot_dphi] centrality_bins or pT_bins missing/too short in config; abort." << std::endl;
    fd->Close();
    fp->Close();
    fj->Close();
    return;
  }

  const int nPtBins = static_cast<int>(pT_bins.size()) - 1;
  const int nCentBins = static_cast<int>(centrality_bins.size()) - 1;

  for (int icent = 0; icent < nCentBins; ++icent) {
    TString tight_name = Form("h_tight_cluster_%d", icent);
    TH1 *htd = dynamic_cast<TH1 *>(fd->Get(tight_name));
    TH1 *htp = dynamic_cast<TH1 *>(fp->Get(tight_name));
    TH1 *htj = dynamic_cast<TH1 *>(fj->Get(tight_name));
    if (!htd || !htp || !htj) {
      std::cerr << "[plot_dphi] Missing " << tight_name.Data()
                << " in one or more files." << std::endl;
      continue;
    }
    if (htd->GetNbinsX() < nPtBins || htp->GetNbinsX() < nPtBins || htj->GetNbinsX() < nPtBins) {
      std::cerr << "[plot_dphi] " << tight_name.Data() << " has fewer bins than pT_bins; skipping cent "
                << icent << std::endl;
      continue;
    }

    for (int ipt = 0; ipt < nPtBins; ++ipt) {
      TString hname = Form("h_dphi_clusterJets_tight_cent%d_pt%d", icent, ipt);
      TH1 *hd = dynamic_cast<TH1 *>(fd->Get(hname));
      TH1 *hp = dynamic_cast<TH1 *>(fp->Get(hname));
      TH1 *hj = dynamic_cast<TH1 *>(fj->Get(hname));

      if (!hd || !hp || !hj) {
        std::cerr << "[plot_dphi] Missing histogram " << hname.Data()
                  << " in one or more files (run RecoEffCalculator_TTreeReader on data + MC)." << std::endl;
        continue;
      }

      const int pt_bin = ipt + 1;
      const double n_tight_d = htd->GetBinContent(pt_bin);
      const double n_tight_p = htp->GetBinContent(pt_bin);
      const double n_tight_j = htj->GetBinContent(pt_bin);
      if (n_tight_d <= 0.0 || n_tight_p <= 0.0 || n_tight_j <= 0.0) {
        std::cerr << "[plot_dphi] Non-positive N_tight in " << tight_name.Data() << " bin " << pt_bin
                  << " for cent " << icent << "; skipping pt " << ipt << std::endl;
        continue;
      }

      TH1 *h1 = (TH1 *)hd->Clone();
      TH1 *h2 = (TH1 *)hp->Clone();
      TH1 *h3 = (TH1 *)hj->Clone();
      h1->SetDirectory(nullptr);
      h2->SetDirectory(nullptr);
      h3->SetDirectory(nullptr);

      const int dphi_rb = 4;
      h1->Rebin(dphi_rb);
      h2->Rebin(dphi_rb);
      h3->Rebin(dphi_rb);

      h1->Scale(1.0 / n_tight_d);
      h2->Scale(1.0 / n_tight_p);
      h3->Scale(1.0 / n_tight_j);

      TCanvas *c = new TCanvas(Form("c_dphi_%d_%d", icent, ipt), "", 700, 600);

      const float ymax = static_cast<float>(
          std::max(h1->GetMaximum(), std::max(h2->GetMaximum(), h3->GetMaximum())));
      h1->GetYaxis()->SetRangeUser(0, ymax * 1.15f);
      h1->SetXTitle("|#Delta#phi_{clus-jet}|");
      h1->GetYaxis()->SetTitle("(1/#it{N}_{tight}) d#it{N}_{pairs}/d|#Delta#phi|");

      h1->SetMarkerStyle(20);
      h1->SetMarkerSize(1.0);
      h1->SetMarkerColor(kBlack);
      h2->SetLineColor(kBlue);
      h2->SetLineWidth(2);
      h3->SetLineColor(kRed);
      h3->SetLineWidth(2);

      h1->Draw("E");
      h2->Draw("HIST SAME");
      h3->Draw("HIST SAME");
      h1->Draw("E SAME");

      myText(      0.55, 0.92, 1, "#bf{#it{sPHENIX}} Internal", 0.04);
      myMarkerText(0.60, 0.87, kBlack, 20, legend_label_data, 1, 0.04);
      myMarkerText(0.60, 0.82, kBlue, 33, legend_label_photon, 1, 0.04);
      myMarkerText(0.60, 0.77, kRed, 34, legend_label_jet, 1, 0.04);

      // Match RecoEffCalculator: Anti-kT Tower jets, R = 0.{cone}; reco jet pT uses calib pT > b2bjet_pT_min
      const std::string cone_label =
          Form("Anti-#it{k}_{T} tower jets, #it{R} = 0.%i", jet_cone_size);
      const std::string jetpt_label =
          Form("#it{p}_{T}^{calib} > %.1f GeV (reco jet)", static_cast<double>(b2bjet_pT_min));
      const std::string pt_label =
          Form("%.0f < #it{E}_{T}^{clus} < %.0f GeV", pT_bins[ipt], pT_bins[ipt + 1]);
      const std::string cent_label =
          Form("%.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]);
      myText(0.20, 0.88, 1, pt_label.c_str(), 0.035);
      myText(0.20, 0.83, 1, cent_label.c_str(), 0.035);
      myText(0.20, 0.78, 1, cone_label.c_str(), 0.032);
      myText(0.20, 0.73, 1, jetpt_label.c_str(), 0.032);
      myText(0.20, 0.68, 1, "Tight clusters (RecoEff)", 0.030);

      c->SaveAs(Form("%s/dphi_clusterJet_cent%d_pt%d.pdf", plot_outdir, icent, ipt));
      delete c;
      delete h1;
      delete h2;
      delete h3;
    }
  }

  // --- Total EMCal energy (events with ≥1 tight cluster), weighted ---
  {
    TH1 *hd_em = dynamic_cast<TH1 *>(fd->Get("h_totalEMCal_energy_tight_weight"));
    TH1 *hp_em = dynamic_cast<TH1 *>(fp->Get("h_totalEMCal_energy_tight_weight"));
    TH1 *hj_em = dynamic_cast<TH1 *>(fj->Get("h_totalEMCal_energy_tight_weight"));
    if (!hd_em || !hp_em || !hj_em) {
      std::cerr << "[plot_dphi] Missing h_totalEMCal_energy_tight_weight in one or more files; "
                   "skip total EMCal plot (re-run RecoEffCalculator_TTreeReader if needed)."
                << std::endl;
    } else {
      TH1 *e1 = (TH1 *)hd_em->Clone("emcal_d");
      TH1 *e2 = (TH1 *)hp_em->Clone("emcal_p");
      TH1 *e3 = (TH1 *)hj_em->Clone("emcal_j");
      e1->SetDirectory(nullptr);
      e2->SetDirectory(nullptr);
      e3->SetDirectory(nullptr);

      const double id = e1->Integral();
      const double ip = e2->Integral();
      const double ij = e3->Integral();
      if (id <= 0.0 || ip <= 0.0 || ij <= 0.0) {
        std::cerr << "[plot_dphi] Non-positive integral on h_totalEMCal_energy_tight_weight; skip EMCal plot."
                  << std::endl;
      } else {
        e1->Scale(1.0 / id);
        e2->Scale(1.0 / ip);
        e3->Scale(1.0 / ij);

        TCanvas *c_em = new TCanvas("c_totalEMCal_energy_tight", "", 700, 600);
        const float ymax_em = static_cast<float>(
            std::max(e1->GetMaximum(), std::max(e2->GetMaximum(), e3->GetMaximum())));
        e1->GetYaxis()->SetRangeUser(0, ymax_em * 1.15f);
        e1->SetXTitle("E_{tot}^{EMCal} [GeV]");
        e1->GetYaxis()->SetTitle("normalized yield");
        e1->GetYaxis()->SetTitleOffset(1.2);
        e1->GetXaxis()->SetNdivisions(505);

        e1->SetMarkerStyle(20);
        e1->SetMarkerSize(1.0);
        e1->SetMarkerColor(kBlack);
        e2->SetLineColor(kBlue);
        e2->SetLineWidth(2);
        e3->SetLineColor(kRed);
        e3->SetLineWidth(2);

        e1->Draw("E");
        e2->Draw("HIST SAME");
        e3->Draw("HIST SAME");
        e1->Draw("E SAME");

        myText(0.55, 0.92, 1, "#bf{#it{sPHENIX}} Internal", 0.04);
        myMarkerText(0.60, 0.87, kBlack, 20, legend_label_data, 1, 0.04);
        myMarkerText(0.60, 0.82, kBlue, 33, legend_label_photon, 1, 0.04);
        myMarkerText(0.60, 0.77, kRed, 34, legend_label_jet, 1, 0.04);


        c_em->SaveAs(Form("%s/totalEMCal_energy_tight_weight.pdf", plot_outdir));
        delete c_em;
      }
      delete e1;
      delete e2;
      delete e3;
    }
  }

  // --- Total IHCal energy (events with ≥1 tight cluster), weighted ---
  {
    TH1 *hd_ih = dynamic_cast<TH1 *>(fd->Get("h_totalIHCal_energy_tight_weight"));
    TH1 *hp_ih = dynamic_cast<TH1 *>(fp->Get("h_totalIHCal_energy_tight_weight"));
    TH1 *hj_ih = dynamic_cast<TH1 *>(fj->Get("h_totalIHCal_energy_tight_weight"));
    if (!hd_ih || !hp_ih || !hj_ih) {
      std::cerr << "[plot_dphi] Missing h_totalIHCal_energy_tight_weight in one or more files; "
                   "skip total IHCal plot (re-run RecoEffCalculator_TTreeReader if needed)."
                << std::endl;
    } else {
      TH1 *i1 = (TH1 *)hd_ih->Clone("ihcal_d");
      TH1 *i2 = (TH1 *)hp_ih->Clone("ihcal_p");
      TH1 *i3 = (TH1 *)hj_ih->Clone("ihcal_j");
      i1->SetDirectory(nullptr);
      i2->SetDirectory(nullptr);
      i3->SetDirectory(nullptr);

      const double iid = i1->Integral();
      const double iip = i2->Integral();
      const double iij = i3->Integral();
      if (iid <= 0.0 || iip <= 0.0 || iij <= 0.0) {
        std::cerr << "[plot_dphi] Non-positive integral on h_totalIHCal_energy_tight_weight; skip IHCal plot."
                  << std::endl;
      } else {
        i1->Scale(1.0 / iid);
        i2->Scale(1.0 / iip);
        i3->Scale(1.0 / iij);

        TCanvas *c_ih = new TCanvas("c_totalIHCal_energy_tight", "", 700, 600);
        const float ymax_ih = static_cast<float>(
            std::max(i1->GetMaximum(), std::max(i2->GetMaximum(), i3->GetMaximum())));
        i1->GetYaxis()->SetRangeUser(0, ymax_ih * 1.15f);
        i1->SetXTitle("E_{tot}^{IHCal} [GeV]");
        i1->GetYaxis()->SetTitle("normalized yield");
        i1->GetYaxis()->SetTitleOffset(1.2);
        i1->GetXaxis()->SetNdivisions(505);

        i1->SetMarkerStyle(20);
        i1->SetMarkerSize(1.0);
        i1->SetMarkerColor(kBlack);
        i2->SetLineColor(kBlue);
        i2->SetLineWidth(2);
        i3->SetLineColor(kRed);
        i3->SetLineWidth(2);

        i1->Draw("E");
        i2->Draw("HIST SAME");
        i3->Draw("HIST SAME");
        i1->Draw("E SAME");

        myText(0.55, 0.92, 1, "#bf{#it{sPHENIX}} Internal", 0.04);
        myMarkerText(0.60, 0.87, kBlack, 20, legend_label_data, 1, 0.04);
        myMarkerText(0.60, 0.82, kBlue, 33, legend_label_photon, 1, 0.04);
        myMarkerText(0.60, 0.77, kRed, 34, legend_label_jet, 1, 0.04);

        c_ih->SaveAs(Form("%s/totalIHCal_energy_tight_weight.pdf", plot_outdir));
        delete c_ih;
      }
      delete i1;
      delete i2;
      delete i3;
    }
  }

  // --- Total OHCal energy (events with ≥1 tight cluster), weighted ---
  {
    TH1 *hd_oh = dynamic_cast<TH1 *>(fd->Get("h_totalOHCal_energy_tight_weight"));
    TH1 *hp_oh = dynamic_cast<TH1 *>(fp->Get("h_totalOHCal_energy_tight_weight"));
    TH1 *hj_oh = dynamic_cast<TH1 *>(fj->Get("h_totalOHCal_energy_tight_weight"));
    if (!hd_oh || !hp_oh || !hj_oh) {
      std::cerr << "[plot_dphi] Missing h_totalOHCal_energy_tight_weight in one or more files; "
                   "skip total OHCal plot (re-run RecoEffCalculator_TTreeReader if needed)."
                << std::endl;
    } else {
      TH1 *o1 = (TH1 *)hd_oh->Clone("ohcal_d");
      TH1 *o2 = (TH1 *)hp_oh->Clone("ohcal_p");
      TH1 *o3 = (TH1 *)hj_oh->Clone("ohcal_j");
      o1->SetDirectory(nullptr);
      o2->SetDirectory(nullptr);
      o3->SetDirectory(nullptr);

      const double iod = o1->Integral();
      const double iop = o2->Integral();
      const double ioj = o3->Integral();
      if (iod <= 0.0 || iop <= 0.0 || ioj <= 0.0) {
        std::cerr << "[plot_dphi] Non-positive integral on h_totalOHCal_energy_tight_weight; skip OHCal plot."
                  << std::endl;
      } else {
        o1->Scale(1.0 / iod);
        o2->Scale(1.0 / iop);
        o3->Scale(1.0 / ioj);

        TCanvas *c_oh = new TCanvas("c_totalOHCal_energy_tight", "", 700, 600);
        const float ymax_oh = static_cast<float>(
            std::max(o1->GetMaximum(), std::max(o2->GetMaximum(), o3->GetMaximum())));
        o1->GetYaxis()->SetRangeUser(0, ymax_oh * 1.15f);
        o1->GetXaxis()->SetRangeUser(0, 700.0);
        o1->SetXTitle("E_{tot}^{OHCal} [GeV]");
        o1->GetYaxis()->SetTitle("normalized yield");
        o1->GetYaxis()->SetTitleOffset(1.2);
        o1->GetXaxis()->SetNdivisions(505);

        o1->SetMarkerStyle(20);
        o1->SetMarkerSize(1.0);
        o1->SetMarkerColor(kBlack);
        o2->SetLineColor(kBlue);
        o2->SetLineWidth(2);
        o3->SetLineColor(kRed);
        o3->SetLineWidth(2);

        o1->Draw("E");
        o2->Draw("HIST SAME");
        o3->Draw("HIST SAME");
        o1->Draw("E SAME");

        myText(0.55, 0.92, 1, "#bf{#it{sPHENIX}} Internal", 0.04);
        myMarkerText(0.60, 0.87, kBlack, 20, legend_label_data, 1, 0.04);
        myMarkerText(0.60, 0.82, kBlue, 33, legend_label_photon, 1, 0.04);
        myMarkerText(0.60, 0.77, kRed, 34, legend_label_jet, 1, 0.04);

        c_oh->SaveAs(Form("%s/totalOHCal_energy_tight_weight.pdf", plot_outdir));
        delete c_oh;
      }
      delete o1;
      delete o2;
      delete o3;
    }
  }

  // --- Centrality (events with ≥1 tight cluster), weighted ---
  {
    TH1 *hd_c = dynamic_cast<TH1 *>(fd->Get("h_centrality_tight_weight"));
    TH1 *hp_c = dynamic_cast<TH1 *>(fp->Get("h_centrality_tight_weight"));
    TH1 *hj_c = dynamic_cast<TH1 *>(fj->Get("h_centrality_tight_weight"));
    if (!hd_c || !hp_c || !hj_c) {
      std::cerr << "[plot_dphi] Missing h_centrality_tight_weight in one or more files; "
                   "skip centrality plot (re-run RecoEffCalculator_TTreeReader if needed)."
                << std::endl;
    } else {
      TH1 *c1 = (TH1 *)hd_c->Clone("cent_d");
      TH1 *c2 = (TH1 *)hp_c->Clone("cent_p");
      TH1 *c3 = (TH1 *)hj_c->Clone("cent_j");
      c1->SetDirectory(nullptr);
      c2->SetDirectory(nullptr);
      c3->SetDirectory(nullptr);

      const double icd = c1->Integral();
      const double icp = c2->Integral();
      const double icj = c3->Integral();
      if (icd <= 0.0 || icp <= 0.0 || icj <= 0.0) {
        std::cerr << "[plot_dphi] Non-positive integral on h_centrality_tight_weight; skip centrality plot."
                  << std::endl;
      } else {
        c1->Scale(1.0 / icd);
        c2->Scale(1.0 / icp);
        c3->Scale(1.0 / icj);

        TCanvas *c_cent = new TCanvas("c_centrality_tight", "", 700, 600);
        const float ymax_c = static_cast<float>(
            std::max(c1->GetMaximum(), std::max(c2->GetMaximum(), c3->GetMaximum())));
        c1->GetYaxis()->SetRangeUser(0, ymax_c * 1.15f);
        c1->SetXTitle("Centrality [%]");
        c1->GetYaxis()->SetTitle("normalized yield");

        c1->SetMarkerStyle(20);
        c1->SetMarkerSize(1.0);
        c1->SetMarkerColor(kBlack);
        c2->SetLineColor(kBlue);
        c2->SetLineWidth(2);
        c3->SetLineColor(kRed);
        c3->SetLineWidth(2);

        c1->Draw("E");
        c2->Draw("HIST SAME");
        c3->Draw("HIST SAME");
        c1->Draw("E SAME");

        myText(0.55, 0.92, 1, "#bf{#it{sPHENIX}} Internal", 0.04);
        myMarkerText(0.60, 0.87, kBlack, 20, legend_label_data, 1, 0.04);
        myMarkerText(0.60, 0.82, kBlue, 33, legend_label_photon, 1, 0.04);
        myMarkerText(0.60, 0.77, kRed, 34, legend_label_jet, 1, 0.04);
        myText(0.20, 0.88, 1, "Events with #geq1 tight cluster (RecoEff)", 0.035);
        myText(0.20, 0.83, 1, "Shape-normalized", 0.032);

        c_cent->SaveAs(Form("%s/centrality_tight_weight.pdf", plot_outdir));
        delete c_cent;
      }
      delete c1;
      delete c2;
      delete c3;
    }
  }

  // --- Vertex z (weighted, all events in RecoEff loop) ---
  {
    TH1 *hd_vz = dynamic_cast<TH1 *>(fd->Get("h_vertexz"));
    TH1 *hp_vz = dynamic_cast<TH1 *>(fp->Get("h_vertexz"));
    TH1 *hj_vz = dynamic_cast<TH1 *>(fj->Get("h_vertexz"));
    if (!hd_vz || !hp_vz || !hj_vz) {
      std::cerr << "[plot_dphi] Missing h_vertexz in one or more files; skip vertex z plot." << std::endl;
    } else {
      TH1 *v1 = (TH1 *)hd_vz->Clone("vz_d");
      TH1 *v2 = (TH1 *)hp_vz->Clone("vz_p");
      TH1 *v3 = (TH1 *)hj_vz->Clone("vz_j");
      v1->SetDirectory(nullptr);
      v2->SetDirectory(nullptr);
      v3->SetDirectory(nullptr);

      const double ivd = v1->Integral();
      const double ivp = v2->Integral();
      const double ivj = v3->Integral();
      if (ivd <= 0.0 || ivp <= 0.0 || ivj <= 0.0) {
        std::cerr << "[plot_dphi] Non-positive integral on h_vertexz; skip vertex z plot." << std::endl;
      } else {
        v1->Scale(1.0 / ivd);
        v2->Scale(1.0 / ivp);
        v3->Scale(1.0 / ivj);

        TCanvas *c_vz = new TCanvas("c_vertexz", "", 700, 600);
        const float ymax_vz = static_cast<float>(
            std::max(v1->GetMaximum(), std::max(v2->GetMaximum(), v3->GetMaximum())));
        v1->GetYaxis()->SetRangeUser(0, ymax_vz * 1.15f);
        v1->SetXTitle("Vertex #it{z} [cm]");
        v1->GetYaxis()->SetTitle("normalized yield");

        v1->SetMarkerStyle(20);
        v1->SetMarkerSize(1.0);
        v1->SetMarkerColor(kBlack);
        v2->SetLineColor(kBlue);
        v2->SetLineWidth(2);
        v3->SetLineColor(kRed);
        v3->SetLineWidth(2);

        v1->Draw("E");
        v2->Draw("HIST SAME");
        v3->Draw("HIST SAME");
        v1->Draw("E SAME");

        myText(0.55, 0.92, 1, "#bf{#it{sPHENIX}} Internal", 0.04);
        myMarkerText(0.60, 0.87, kBlack, 20, legend_label_data, 1, 0.04);
        myMarkerText(0.60, 0.82, kBlue, 33, legend_label_photon, 1, 0.04);
        myMarkerText(0.60, 0.77, kRed, 34, legend_label_jet, 1, 0.04);
        myText(0.20, 0.88, 1, "Weighted vertex #it{z} (RecoEff)", 0.035);
        myText(0.20, 0.83, 1, "Shape-normalized", 0.032);

        c_vz->SaveAs(Form("%s/vertexz.pdf", plot_outdir));
        delete c_vz;
      }
      delete v1;
      delete v2;
      delete v3;
    }
  }

  fd->Close();
  fp->Close();
  fj->Close();
  std::cout << "[plot_dphi] Plots saved to " << plot_outdir << std::endl;
}
