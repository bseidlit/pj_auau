#include <yaml-cpp/yaml.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TSystem.h>
#include <cmath>
#include <memory>
#include <vector>
#include <map>
#include <string>

// Cluster ET threshold (GeV) for "cluster above 10 GeV" event selection and pT spectrum
const float CLUSTER_ET_MIN_GEV = 10.0;
// EMCAL time sample period for cluster time (ns)
const float TIME_SAMPLE_NS = 17.6;
// Index in currentscaler_raw for minimum-bias trigger (12th element = index 11)
const int MB_SCALER_INDEX = 12;

void EvtCharacter(const std::string &configname = "nom.yaml", const std::string filetype = "jet10")
{
  gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
  YAML::Node configYaml = YAML::LoadFile(configname);

  bool issim = true;
  if (filetype == "data" || filetype == "data_pp" || filetype == "data_aa")
    issim = false;

  // Parse filetype into base and collision system (pp vs Au+Au)
  std::string filetype_base = filetype;
  bool is_pp = true;
  if (filetype.size() > 3 && filetype.compare(filetype.size() - 3, 3, "_aa") == 0) {
    filetype_base = filetype.substr(0, filetype.size() - 3);
    is_pp = false;
  } else if (filetype.size() > 3 && filetype.compare(filetype.size() - 3, 3, "_pp") == 0) {
    filetype_base = filetype.substr(0, filetype.size() - 3);
    is_pp = true;
  }
  std::cout << "filetype: " << filetype << " -> base: " << filetype_base << ", system: " << (is_pp ? "pp" : "Au+Au (aa)") << std::endl;

  std::string infilename_root_dir = configYaml["input"]["photon_jet_file_root_dir"].as<std::string>();
  std::string infilename_branch_dir = configYaml["input"]["photon_jet_file_branch_dir"].as<std::string>();
  std::string infilename = infilename_root_dir + filetype + infilename_branch_dir;

  if (!issim) {
    if (filetype == "data_aa" && configYaml["input"]["data_file_aa"])
      infilename = configYaml["input"]["data_file_aa"].as<std::string>();
    else if ((filetype == "data_pp" || filetype == "data") && configYaml["input"]["data_file_pp"])
      infilename = configYaml["input"]["data_file_pp"].as<std::string>();
    else
      infilename = configYaml["input"]["data_file"].as<std::string>();
  }
  std::cout << "infilename: " << infilename << std::endl;

  float vertexcut = configYaml["analysis"]["vertex_cut"].as<float>();
  std::vector<float> eta_bins = configYaml["analysis"]["eta_bins"].as<std::vector<float>>();
  float jet_eta_max = configYaml["analysis"]["jet_eta"].as<float>(0.6f);
  std::vector<float> centrality_bins = configYaml["analysis"]["centrality_bins"].as<std::vector<float>>();
  int n_cent_bins = centrality_bins.size() - 1;
  std::vector<double> cent_edges(centrality_bins.begin(), centrality_bins.end());

  std::vector<float> pT_bins = configYaml["analysis"]["pT_bins"].as<std::vector<float>>();
  int n_pT_bins = pT_bins.size() - 1;

  const int npb_cut_on = configYaml["analysis"]["npb_cut_on"].as<int>(0);
  const float npb_score_cut = configYaml["analysis"]["npb_score_cut"].as<float>(0.5f);

  // Jet collection: Anti-kT tower jets; jet_cone_size 3 -> R=0.3, 4 -> R=0.4
  const int jet_cone_size = configYaml["analysis"]["jet_cone_size"].as<int>(3);
  const bool use_jet_r04 = (jet_cone_size == 4);
  std::cout << "[EvtCharacter] jet_cone_size=" << jet_cone_size << " -> Anti-kT Tower r0"
            << (use_jet_r04 ? "4" : "3") << "_Sub1" << std::endl;
  std::cout << "[EvtCharacter] npb_cut_on=" << npb_cut_on << ", npb_score_cut=" << npb_score_cut << std::endl;

  std::vector<int> trigger_used;
  {
    YAML::Node trigNode = configYaml["analysis"]["trigger_used"];
    if (trigNode && trigNode.IsSequence())
      trigger_used = trigNode.as<std::vector<int>>();
    else
      trigger_used.push_back(configYaml["analysis"]["trigger_used"].as<int>());
  }

  // Vertex reweighting for simulation
  TH1 *h_vertex_reweight = nullptr;
  int vertex_reweight_on = 1;
  std::string vertex_reweight_file = "results/vertex_reweight.root";
  if (issim) {
    vertex_reweight_on = configYaml["analysis"]["vertex_reweight_on"].as<int>(1);
    vertex_reweight_file = configYaml["analysis"]["vertex_reweight_file"].as<std::string>("results/vertex_reweight.root");
    if (vertex_reweight_on) {
      TFile *fvtx = TFile::Open(vertex_reweight_file.c_str(), "READ");
      if (!fvtx || fvtx->IsZombie()) {
        std::cerr << "[VertexReweight] ERROR: cannot open " << vertex_reweight_file << std::endl;
        return;
      }
      TH1 *htmp = dynamic_cast<TH1 *>(fvtx->Get("h_vertexz_ratio_data_over_mccombined"));
      if (!htmp) {
        std::cerr << "[VertexReweight] ERROR: histogram not found in " << vertex_reweight_file << std::endl;
        fvtx->Close();
        delete fvtx;
        return;
      }
      h_vertex_reweight = dynamic_cast<TH1 *>(htmp->Clone("h_vertexz_ratio_data_over_mccombined_clone"));
      h_vertex_reweight->SetDirectory(nullptr);
      h_vertex_reweight->Sumw2();
      std::cout << "[VertexReweight] Using " << vertex_reweight_file << std::endl;
    }
  }

  std::string outfilename = configYaml["output"]["evt_character_outfile"].as<std::string>() + filetype + ".root";

  std::string treename = configYaml["input"]["tree"].as<std::string>();
  std::string clusternodename = configYaml["input"]["cluster_node_name"].as<std::string>();

  TChain chain(treename.c_str());
  chain.Add(infilename.c_str());
  if (chain.GetListOfFiles()->GetEntries() == 0) {
    std::cerr << "Error: No files found for path: " << infilename << std::endl;
    return;
  }



  TTreeReader reader(&chain);

  TTreeReaderValue<int> runnumber(reader, "runnumber");
  TTreeReaderValue<float> vertexz(reader, "vertexz");
  TTreeReaderValue<float> Psi2(reader, "Psi2");
  TTreeReaderValue<int> ncluster(reader, Form("ncluster_%s", clusternodename.c_str()));
  TTreeReaderArray<Bool_t> scaledtrigger(reader, "scaledtrigger");

  // currentscaler_raw[64]/L in tree: 12th element (index 11) = MB trigger
  TTreeReaderArray<Long64_t> *currentscaler_raw = nullptr;
  if (!issim && chain.FindBranch("currentscaler_raw")) {
    currentscaler_raw = new TTreeReaderArray<Long64_t>(reader, "currentscaler_raw");
    std::cout << "Using currentscaler_raw (element " << MB_SCALER_INDEX << " = MB) for collisions per run." << std::endl;
  }

  std::unique_ptr<TTreeReaderValue<float>> cent_reader;
  if (!is_pp)
    cent_reader = std::make_unique<TTreeReaderValue<float>>(reader, "cent");

  TTreeReaderArray<float> cluster_Et(reader, Form("cluster_Et_%s", clusternodename.c_str()));
  TTreeReaderArray<float> cluster_Eta(reader, Form("cluster_Eta_%s", clusternodename.c_str()));
  TTreeReaderArray<float> cluster_Phi(reader, Form("cluster_Phi_%s", clusternodename.c_str()));
  TTreeReaderArray<float> cluster_ietacent(reader, Form("cluster_ietacent_%s", clusternodename.c_str()));
  TTreeReaderArray<float> cluster_iphicent(reader, Form("cluster_iphicent_%s", clusternodename.c_str()));
  TTreeReaderArray<float> cluster_e_array(reader, Form("cluster_e_array_%s", clusternodename.c_str()));
  TTreeReaderArray<float> cluster_time_array(reader, Form("cluster_time_array_%s", clusternodename.c_str()));
  TTreeReaderArray<int> cluster_ownership_array(reader, Form("cluster_ownership_array_%s", clusternodename.c_str()));
  TTreeReaderArray<float> cluster_npb_score(reader, Form("cluster_npb_score_%s", clusternodename.c_str()));

  //TTreeReaderArray<float> cluster_iso_03_70_emcal(reader, Form("cluster_iso_03_70_emcal_%s", clusternodename.c_str()));
  //TTreeReaderArray<float> cluster_iso_03_70_hcalin(reader, Form("cluster_iso_03_70_hcalin_%s", clusternodename.c_str()));
  //TTreeReaderArray<float> cluster_iso_03_70_hcalout(reader, Form("cluster_iso_03_70_hcalout_%s", clusternodename.c_str()));

  TTreeReaderArray<float> cluster_iso_03_70_emcal(reader, Form("cluster_iso_04_sub1_emcal_%s", clusternodename.c_str()));
  TTreeReaderArray<float> cluster_iso_03_70_hcalin(reader, Form("cluster_iso_04_sub1_hcalin_%s", clusternodename.c_str()));
  TTreeReaderArray<float> cluster_iso_03_70_hcalout(reader, Form("cluster_iso_04_sub1_hcalout_%s", clusternodename.c_str()));


  TTreeReaderValue<int> njet_R3(reader, "njet_AntiKt_Tower_r03_Sub1");
  TTreeReaderArray<float> jet_Pt_R3(reader, "jet_Pt_AntiKt_Tower_r03_Sub1");
  TTreeReaderArray<float> jet_Eta_R3(reader, "jet_Eta_AntiKt_Tower_r03_Sub1");
  TTreeReaderArray<float> jet_Phi_R3(reader, "jet_Phi_AntiKt_Tower_r03_Sub1");
  
  TTreeReaderValue<int> njet_R4(reader, "njet_AntiKt_Tower_r04_Sub1");
  TTreeReaderArray<float> jet_Pt_R4(reader, "jet_Pt_AntiKt_Tower_r04_Sub1");
  TTreeReaderArray<float> jet_Eta_R4(reader, "jet_Eta_AntiKt_Tower_r04_Sub1");
  TTreeReaderArray<float> jet_Phi_R4(reader, "jet_Phi_AntiKt_Tower_r04_Sub1");

  TTreeReaderValue<int> &njet = use_jet_r04 ? njet_R4 : njet_R3;
  TTreeReaderArray<float> &jet_Pt = use_jet_r04 ? jet_Pt_R4 : jet_Pt_R3;
  TTreeReaderArray<float> &jet_Eta = use_jet_r04 ? jet_Eta_R4 : jet_Eta_R3;
  TTreeReaderArray<float> &jet_Phi = use_jet_r04 ? jet_Phi_R4 : jet_Phi_R3;

  TFile *fout = new TFile(outfilename.c_str(), "RECREATE");
  fout->cd();

  // Event-level histograms
  TH1D *h_centrality = new TH1D("h_centrality", "Centrality distribution;Centrality [%];Events", n_cent_bins, cent_edges.data());
  TH1D *h_vertexz = new TH1D("h_vertexz", "Vertex Z;Vertex Z [cm];Events", 200, -100, 100);
  TH1D *h_Psi2 = new TH1D("h_Psi2", "Event plane #Psi_{2};#Psi_{2} [rad];Events", 100, -M_PI, M_PI);
  TH1D *h_centrality_with_cluster10 = new TH1D("h_centrality_with_cluster10", "Centrality (events with at least one cluster E_{T} > 10 GeV);Centrality [%];Events", n_cent_bins, cent_edges.data());
  TH1D *h_cluster_pt_above10 = new TH1D("h_cluster_pt_above10", "Cluster E_{T} (E_{T} > 10 GeV);Cluster E_{T} [GeV];Clusters", 80, 10, 50);
  TH2D *h_cluster_eta_phi = new TH2D("h_cluster_eta_phi", ";#eta;#phi [rad]", 50, -1.0, 1.0, 64, -M_PI, M_PI);
  TH2D *h_cluster_ieta_iphi = new TH2D("h_cluster_ieta_iphi", ";i#eta;i#phi", 96, -0.5, 95.5, 256, -0.5, 255.5);
  TH2D *h_cluster_time_et = new TH2D("h_cluster_time_et", "Cluster time vs E_{T} (E_{T} > 10 GeV);Cluster E_{T} [GeV];Cluster time [ns]", 80, 10, 50, 100, -30, 30);
  TH1D *h_jet_pt = new TH1D("h_jet_pt", "Jet p_{T};Jet p_{T} [GeV];Jets", 50, 0, 100);
  TH2D *h_jet_eta_phi = new TH2D("h_jet_eta_phi", "Jet yield (p_{T} > 10 GeV);#eta;#phi [rad]", 10, -1.0, 1.0, 10, -M_PI, M_PI);

  TH2D *h_Psi2_vs_jetPt = new TH2D(
      "h_Psi2_vs_jetPt",
      "#Psi_{2} vs jet p_{T} (jets with p_{T} > 10 GeV, |#eta| < jet cut);Jet p_{T} [GeV];#Psi_{2} [rad]",
      50, 0, 100, 100, -M_PI, M_PI);

  // Iso ET vs centrality (1% bins), one TH2 per pT bin
  std::vector<TH2D *> h_iso_vs_cent_pt;
  for (int ipt = 0; ipt < n_pT_bins; ipt++) {
    h_iso_vs_cent_pt.push_back(new TH2D(Form("h_iso_vs_cent_pt%d", ipt),
      Form("Iso E_{T} vs centrality (%.0f < E_{T} < %.0f GeV);Centrality [%%];Iso E_{T} [GeV]", pT_bins[ipt], pT_bins[ipt + 1]),
      100, 0, 100, 100, -50, 50));
    h_iso_vs_cent_pt.back()->Sumw2();
  }

  // Cluster yield vs centrality (5% bins), one TH1 per analysis pT bin
  std::vector<TH1D *> h_cluster_yield_vs_cent_pt;
  for (int ipt = 0; ipt < n_pT_bins; ipt++) {
    h_cluster_yield_vs_cent_pt.push_back(new TH1D(
        Form("h_cluster_yield_vs_cent_pt%d", ipt),
        Form("Cluster yield vs centrality (%.0f < E_{T} < %.0f GeV);Centrality [%%];Clusters",
             pT_bins[ipt], pT_bins[ipt + 1]),
        20, 0., 100.));
    h_cluster_yield_vs_cent_pt.back()->Sumw2();
  }

  const float dR_min = 0.2f;
  std::vector<std::vector<TH1D *>> h_dphi_clusterJets;
  h_dphi_clusterJets.resize(n_cent_bins);
  for (int icent = 0; icent < n_cent_bins; ++icent)
  {
    h_dphi_clusterJets[icent].resize(n_pT_bins, nullptr);
    for (int ipt = 0; ipt < n_pT_bins; ++ipt)
    {
      h_dphi_clusterJets[icent][ipt] = new TH1D(
          Form("h_dphi_clusterJets_cent%d_pt%d", icent, ipt),
          Form("|#Delta#phi|(cluster, jet) for %.0f-%.0f%% cent and %.0f < E_{T}^{clus} < %.0f GeV;|#Delta#phi|;Pairs",
               centrality_bins[icent], centrality_bins[icent + 1], pT_bins[ipt], pT_bins[ipt + 1]),
          64, 0, M_PI);
      h_dphi_clusterJets[icent][ipt]->Sumw2();
    }
  }

  // Per-run: fill these from accumulated run_max_scaler and run_nclusters_above10
  TH1D *h_scaler_per_run = new TH1D("h_scaler_per_run", ";Runs;Max MB scaler", 14000, 66000, 80000);
  h_scaler_per_run->Sumw2();
  TH1D *h_clusters_per_run = new TH1D("h_clusters_per_run", ";Runs;N clusters", 14000, 66000, 80000);
  h_clusters_per_run->Sumw2();

  // Per-run accumulation (runs mixed in tree)
  std::map<int, long long> run_max_scaler;
  std::map<int, int> run_nclusters_above10;

  int nentries = chain.GetEntries();
  int ientry = 0;

  while (reader.Next()) {
    if (ientry % 50000 == 0)
      std::cout << "Processing entry " << ientry << " / " << nentries << std::endl;

    int run = *runnumber;
    if (!issim && currentscaler_raw && currentscaler_raw->GetSize() > (unsigned int)MB_SCALER_INDEX) {
      Long64_t val = (*currentscaler_raw)[MB_SCALER_INDEX];
      if (run_max_scaler.count(run))
        run_max_scaler[run] = std::max(run_max_scaler[run], (long long)val);
      else
        run_max_scaler[run] = (long long)val;
    }

    if (!issim) {
      bool any_trigger_fired = false;
      const auto ntrig = scaledtrigger.GetSize();
      for (int itrig : trigger_used) {
        if (itrig >= 0 && (unsigned int)itrig < ntrig && scaledtrigger[itrig] != 0) {
          any_trigger_fired = true;
          break;
        }
      }
      if (!any_trigger_fired)
        continue;
    }

    float weight = 1.0;
    if (issim && vertex_reweight_on && h_vertex_reweight) {
      int bin = h_vertex_reweight->FindBin(*vertexz);
      if (bin < 1) bin = 1;
      if (bin > h_vertex_reweight->GetNbinsX()) bin = h_vertex_reweight->GetNbinsX();
      weight = h_vertex_reweight->GetBinContent(bin);
      if (!std::isfinite(weight) || weight <= 0.0)
        weight = 1.0;
    }

    if (std::abs(*vertexz) > vertexcut)
      continue;

    int centbin = -1;
    float cent_percent = 50.0f;
    if (is_pp) {
      centbin = 0;
    } else {
      cent_percent = (*(*cent_reader)) * 100.0f;
      for (int icent = 0; icent < n_cent_bins; icent++) {
        if (cent_percent >= centrality_bins[icent] && cent_percent < centrality_bins[icent + 1]) {
          centbin = icent;
          break;
        }
      }
    }

    if (centbin < 0)
      continue;

    h_vertexz->Fill(*vertexz, weight);
    h_centrality->Fill(cent_percent, weight);
    if (std::isfinite(*Psi2) && *Psi2 > -9000.f)
      h_Psi2->Fill(*Psi2, weight);

    int n_clusters_this_evt = 0;
    bool has_cluster_above_10 = false;
    for (int icluster = 0; icluster < *ncluster; icluster++) {
      if (cluster_Et[icluster] < CLUSTER_ET_MIN_GEV)
        continue;
      if (cluster_npb_score[icluster] < npb_score_cut)
        continue;
      
      float eta = cluster_Eta[icluster];
      if (eta < eta_bins[0] || eta > eta_bins[eta_bins.size() - 1])
        continue;

      has_cluster_above_10 = true;
      n_clusters_this_evt++;
      h_cluster_pt_above10->Fill(cluster_Et[icluster], weight);
      h_cluster_eta_phi->Fill(eta, cluster_Phi[icluster], weight);
      h_cluster_ieta_iphi->Fill(cluster_ietacent[icluster], cluster_iphicent[icluster], weight);

      float recoisoET = cluster_iso_03_70_emcal[icluster] + cluster_iso_03_70_hcalin[icluster] + cluster_iso_03_70_hcalout[icluster];
      int pTbin = -1;
      float clusterET = cluster_Et[icluster];
      for (int ipt = 0; ipt < n_pT_bins; ipt++) {
        if (clusterET >= pT_bins[ipt] && clusterET < pT_bins[ipt + 1]) { pTbin = ipt; break; }
      }
      if (pTbin >= 0)
      {
        h_iso_vs_cent_pt[pTbin]->Fill(cent_percent, recoisoET, weight);
        h_cluster_yield_vs_cent_pt[pTbin]->Fill(cent_percent, weight);

        // Fill |DeltaPhi|(cluster, jet) for this cluster pT bin.
        for (int ijet = 0; ijet < *njet; ijet++)
        {
          if (jet_Pt[ijet] < 10.0f)
            continue;
          float jeta = jet_Eta[ijet];
          if (std::abs(jeta) > jet_eta_max)
            continue;

          const float deta = cluster_Eta[icluster] - jeta;
          float dphi = cluster_Phi[icluster] - jet_Phi[ijet];
          while (dphi > M_PI)
            dphi -= 2 * M_PI;
          while (dphi < -M_PI)
            dphi += 2 * M_PI;

          const float dR = std::sqrt(deta * deta + dphi * dphi);
          if (dR <= dR_min)
            continue;

          h_dphi_clusterJets[centbin][pTbin]->Fill(std::abs(dphi), weight);
        }
      }

      float clusteravgtime = 0;
      float cluster_total_e = 0;
      for (int i = 0; i < 49; i++) {
        int idx = icluster * 49 + i;
        if (cluster_ownership_array[idx] == 1) {
          clusteravgtime += cluster_time_array[idx] * cluster_e_array[idx];
          cluster_total_e += cluster_e_array[idx];
        }
      }
      clusteravgtime = cluster_total_e > 0 ? clusteravgtime / cluster_total_e : 0;
      clusteravgtime *= TIME_SAMPLE_NS;
      h_cluster_time_et->Fill(cluster_Et[icluster], clusteravgtime, weight);
    }

    if (has_cluster_above_10)
      h_centrality_with_cluster10->Fill(cent_percent, weight);

    if (has_cluster_above_10)
    {
      const bool psi2_ok = std::isfinite(*Psi2) && *Psi2 > -9000.f;
      for (int ijet = 0; ijet < *njet; ijet++) {
        h_jet_pt->Fill(jet_Pt[ijet], weight);
        if (jet_Pt[ijet] < 8.0f) continue;
        float jeta = jet_Eta[ijet];
        if (std::abs(jeta) > jet_eta_max) continue;
        h_jet_eta_phi->Fill(jeta, jet_Phi[ijet], weight);
        if (psi2_ok){
          float dphi = fabs(*Psi2 - jet_Phi[ijet]);
          dphi = fmod(dphi, M_PI);
          if (dphi < 0) dphi += M_PI;
          if (dphi > M_PI/2) dphi = M_PI - dphi;
          h_Psi2_vs_jetPt->Fill(jet_Pt[ijet],dphi, weight);
        }
      }
    }

    run_nclusters_above10[run] += n_clusters_this_evt;

    ientry++;
  }

  for (auto &p : run_max_scaler){
    h_scaler_per_run->Fill(p.first, p.second);
    std::cout << "h_scaler_per_run: " << p.first << " " << p.second << std::endl;
  }
    
  for (auto &p : run_nclusters_above10){
    h_clusters_per_run->Fill(p.first, p.second);
    std::cout << "h_clusters_per_run: " << p.first << " " << p.second << std::endl;
  }

  // Non-overlapping sums over 5 adjacent run bins (same grouping as PlotEvtCharacter.C).
  // Errors: quadrature sum of sqrt(count) per run bin (matches plotting macro).
  const int kRunsPerGroup = 5;
  auto build_and_write_run_groups = [&](const TH1D *h, const char *name, const char *xtitle, const char *ytitle) {
    const int nb = h->GetNbinsX();
    if (nb <= 0)
      return;
    const int ngrp = (nb + kRunsPerGroup - 1) / kRunsPerGroup;
    const double x0 = h->GetXaxis()->GetXmin();
    auto *hg = new TH1D(name, Form(";%s;%s", xtitle, ytitle), ngrp, x0,
                        x0 + static_cast<double>(ngrp * kRunsPerGroup));
    hg->Sumw2();
    for (int ig = 0; ig < ngrp; ++ig) {
      double sum = 0.;
      double err2 = 0.;
      for (int k = 0; k < kRunsPerGroup; ++k) {
        const int ib = ig * kRunsPerGroup + k + 1;
        if (ib > nb)
          break;
        const double c = h->GetBinContent(ib);
        sum += c;
        const double sigma = std::sqrt(std::fmax(0.0, c));
        err2 += sigma * sigma;
      }
      hg->SetBinContent(ig + 1, sum);
      hg->SetBinError(ig + 1, std::sqrt(err2));
    }
    fout->cd();
    hg->Write();
    delete hg;
  };

  build_and_write_run_groups(h_scaler_per_run, "h_scaler_per_5runs", "First run # in 5-run group",
                             "Total MB scaler (sum of 5 runs)");
  build_and_write_run_groups(h_clusters_per_run, "h_clusters_per_5runs", "First run # in 5-run group",
                             "N clusters (sum of 5 runs)");

  fout->cd();
  fout->Write();
  fout->Close();
  std::cout << "EvtCharacter done. Output: " << outfilename << std::endl;
}
