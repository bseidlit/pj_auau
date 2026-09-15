#pragma once
// Frozen regression reference. Normal analysis uses rdf/PhotonJetRDF.C.
//
//   root -l -b -q 'tests/reference/PhotonJetHistMaker.C("configs/pp/config_pj_pp_nom.yaml","data","chunks/pj_pp_nom/data_000.list","data_000")'
//
//   product   : data | sim_signal | sim_inclusive
//   chunk_list: text file, one part per line, "index<space>path" where index is
//               the 0-based line number of the part in the product's files.txt
//               (needed for the sample map). Bare chunk paths work without a sample map.
//               Empty -> configured master files.txt, whose bare paths get canonical indices.
//   chunk_tag : empty -> write the final PPG12 file names directly,
//               otherwise write <results>/chunks/<name>.<tag>.root for hadd.
//
// Output names (identical to RecoEffCalculator_TTreeReader.C):
//   data          -> {output.data_outfile}_{var_type}.root
//   sim_signal    -> {output.eff_outfile}_{var_type}.root + {output.response_outfile}_{var_type}.root
//   sim_inclusive -> {output.eff_outfile}_jet_{var_type}.root
//
// Only eta index 0 is produced (eta_bins must have exactly one bin), because
// CalculatePhotonYield.C reads `_0` only. Au+Au centrality classes are
// separate configs / var_types.
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TString.h>

#include <RooUnfoldResponse.h>

#include <yaml-cpp/yaml.h>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "PhotonJetSelection.h"
#include "PhotonJetWeights.h"
#include "PhotonJetIO.h"

namespace
{
TH1D *Book1(const char *name, const char *title, const std::vector<double> &edges)
{
    TH1D *h = new TH1D(name, title, (int)edges.size() - 1, edges.data());
    h->Sumw2();
    return h;
}
} // namespace

static void MakeHistograms(const std::string &configname = "",
                        const std::string &product = "data",
                        const std::string &chunk_list = "",
                        const std::string &chunk_tag = "")
{
    const PJ::Product kind = PJ::ParseProduct(product);
    const PJ::ConfigFile config = PJ::ReadConfig(configname);
    const YAML::Node &cfg = config.yaml;
    const PJ::Cuts c = PJ::LoadCuts(cfg);
    const PJ::InputMap imap = PJ::InputMapFor(c.system);
    const bool issim = kind != PJ::Product::Data;
    const bool is_signal = kind == PJ::Product::Signal;
    const bool need_sample = issim && c.weight_mode != "stored";
    if (issim && c.system == "auau" && c.weight_mode != "sample_map")
        throw std::runtime_error("Au+Au simulation requires sample_map weights");

    // ---- inputs ------------------------------------------------------------
    std::string list_path = chunk_list;
    if (list_path.empty())
    {
        list_path = PJ::RequiredString(cfg["photonjet"], PJ::InputListKey(kind));
    }
    const std::vector<PJ::ChunkEntry> parts = PJ::ReadParts(list_path, chunk_list.empty());
    std::cout << "[PJ] product=" << product << " system=" << c.system << " parts=" << parts.size()
              << " threshold_source=" << c.threshold_source << " weight_mode=" << c.weight_mode << std::endl;

    // per-product sample map override (Au+Au: signal and inclusive products differ)
    std::string sample_map_file = c.sample_map_file;
    if (issim)
    {
        const char *pkey = is_signal ? "sample_map_signal" : "sample_map_inclusive";
        if (cfg["photonjet"] && cfg["photonjet"][pkey]) sample_map_file = cfg["photonjet"][pkey].as<std::string>();
    }
    const std::vector<PJ::SampleRange> sample_map = need_sample ? PJ::LoadSampleMap(sample_map_file) : std::vector<PJ::SampleRange>{};
    if (need_sample && sample_map.empty())
        throw std::runtime_error("weight_mode " + c.weight_mode + " needs photonjet.sample_map");
    for (const auto &part : parts)
        if (need_sample && !PPG12::GetSampleConfig(PJ::SampleForPart(sample_map, part.index)).valid)
            throw std::runtime_error("no sample for part index " + std::to_string(part.index) +
                "; use the master files.txt or an indexed chunk list");
    PJ::VertexWeight vtxw = need_sample ? PJ::LoadVertexWeight(c.vertex_weight_file) : PJ::VertexWeight{};
    std::vector<std::string> dependencies;
    if (!issim && !c.run_list_file.empty()) dependencies.push_back(c.run_list_file);
    if (need_sample) { dependencies.push_back(sample_map_file); dependencies.push_back(c.vertex_weight_file); }
    YAML::Node provenance = PJ::RunProvenance(config, kind, c, list_path, dependencies);

    // ---- outputs -----------------------------------------------------------
    const PJ::OutputNames names = PJ::OutputPaths(cfg, kind, chunk_tag);
    PJ::RootOutput main_file(names.main);
    std::unique_ptr<PJ::RootOutput> response_file;
    if (is_signal) response_file = std::make_unique<PJ::RootOutput>(names.response);
    TFile *fout = main_file.get();
    TFile *fresp = response_file ? response_file->get() : nullptr;
    fout->cd();

    // ---- booking (names/binning = RecoEff, eta index 0) ---------------------
    const std::vector<double> &pb = c.pT_bins;
    const std::vector<double> &tb = c.pT_bins_truth;
    const double pTmin = pb.front(), pTmax = pb.back();
    const double pTmin_truth = tb.front(), pTmax_truth = tb.back();
    const char *et = Form("%.1f < eta < %.1f", c.eta_bins[0], c.eta_bins[1]);

    TH1D *h_tight_iso = Book1("h_tight_iso_cluster_0", Form("Tight Iso Cluster %s", et), pb);
    TH1D *h_tight_noniso = Book1("h_tight_noniso_cluster_0", Form("Tight Non-Iso Cluster %s", et), pb);
    TH1D *h_nontight_iso = Book1("h_nontight_iso_cluster_0", Form("Non-Tight Iso Cluster %s", et), pb);
    TH1D *h_nontight_noniso = Book1("h_nontight_noniso_cluster_0", Form("Non-Tight Non-Iso Cluster %s", et), pb);
    TH1D *h_common = Book1("h_common_cluster_0", Form("Common Cluster %s", et), pb);
    TH1D *h_all = Book1("h_all_cluster_0", Form("All Cluster %s", et), pb);
    TH1D *h_tight = Book1("h_tight_cluster_0", Form("Tight Cluster %s", et), pb);
    TH1D *h_sig[4] = {Book1("h_tight_iso_cluster_signal_0", Form("Tight Iso Cluster %s", et), pb),
                      Book1("h_tight_noniso_cluster_signal_0", Form("Tight Non-Iso Cluster %s", et), pb),
                      Book1("h_nontight_iso_cluster_signal_0", Form("Non-Tight Iso Cluster %s", et), pb),
                      Book1("h_nontight_noniso_cluster_signal_0", Form("Non-Tight Non-Iso Cluster %s", et), pb)};
    TH1D *h_all_sig = Book1("h_all_cluster_signal_0", Form("All Cluster %s", et), pb);
    TH1D *h_tight_sig = Book1("h_tight_cluster_signal_0", Form("Tight Cluster %s", et), pb);
    TH1D *h_nm[4] = {Book1("h_tight_iso_cluster_notmatch_0", Form("Tight Iso Cluster %s", et), pb),
                     Book1("h_tight_noniso_cluster_notmatch_0", Form("Tight Non-Iso Cluster %s", et), pb),
                     Book1("h_nontight_iso_cluster_notmatch_0", Form("Non-Tight Iso Cluster %s", et), pb),
                     Book1("h_nontight_noniso_cluster_notmatch_0", Form("Non-Tight Non-Iso Cluster %s", et), pb)};
    TH1D *h_abcd[4] = {h_tight_iso, h_tight_noniso, h_nontight_iso, h_nontight_noniso};

    TH1D *h_truth = Book1("h_truth_pT_0", Form("Truth pT %s", et), tb);
    TH1D *h_truth_novtx = Book1("h_truth_pT_novtx_0", Form("Truth pT (no vtx wt) %s", et), tb);
    TH1D *h_truth_vc = Book1("h_truth_pT_vertexcut_0", Form("Truth pT %s", et), tb);
    TH1D *h_truth_vc_mbd = Book1("h_truth_pT_vertexcut_mbd_cut_0", Form("Truth pT %s", et), tb);
    TH1D *h_truth_vc_n = Book1("h_truth_pT_vertexcut_mbd_north_cut_0", Form("Truth pT %s", et), tb);
    TH1D *h_truth_vc_s = Book1("h_truth_pT_vertexcut_mbd_south_cut_0", Form("Truth pT %s", et), tb);
    TH1D *h_truth_vc_on = Book1("h_truth_pT_vertexcut_mbd_only_north_0", Form("Truth pT %s", et), tb);
    TH1D *h_truth_vc_os = Book1("h_truth_pT_vertexcut_mbd_only_south_0", Form("Truth pT %s", et), tb);
    TH1D *h_truth_vc_ne = Book1("h_truth_pT_vertexcut_mbd_neither_0", Form("Truth pT %s", et), tb);

    const TEfficiency::EStatOption effopt = TEfficiency::kBUniform;
    auto BookEff = [&](const char *name, const char *title) {
        TEfficiency *e = new TEfficiency(name, Form("%s %s", title, et), (int)tb.size() - 1, tb.data());
        e->SetStatisticOption(effopt);
        e->SetUseWeightedEvents();
        e->SetDirectory(fout);
        return e;
    };
    TEfficiency *eff_reco = BookEff("eff_reco_eta_0", "Reco Efficiency");
    TEfficiency *eff_iso = BookEff("eff_iso_eta_0", "Iso Efficiency");
    TEfficiency *eff_id = BookEff("eff_id_eta_0", "ID Efficiency");
    TEfficiency *eff_all = BookEff("eff_all_eta_0", "All Efficiency");
    TEfficiency *eff_converts = BookEff("eff_converts_eta_0", "Conversion Prob");

    // response (signal only, written to the response file)
    TH1D *h_pT_truth_resp = nullptr, *h_pT_reco_resp = nullptr;
    TH2D *h_resp_full = nullptr;
    RooUnfoldResponse *response = nullptr;
    if (is_signal)
    {
        fresp->cd();
        h_pT_truth_resp = Book1("h_pT_truth_response_0", Form("Truth pT %s", et), tb);
        h_pT_reco_resp = Book1("h_pT_reco_response_0", Form("Reco pT %s", et), pb);
        h_resp_full = new TH2D("h_response_full_0", Form("Response Matrix %s", et), (int)pb.size() - 1, pb.data(), (int)tb.size() - 1, tb.data());
        h_resp_full->Sumw2();
        response = new RooUnfoldResponse((const TH1 *)h_pT_reco_resp, (const TH1 *)h_pT_truth_resp, h_resp_full, "response_matrix_full_0", "", false);
        fout->cd();
    }
    std::unique_ptr<TF1> f_reweight;
    if (is_signal && c.unfold_reweight)
    {
        f_reweight = std::make_unique<TF1>("f_reweight", c.trw_formula.c_str(), c.trw_xmin, c.trw_xmax);
        if (!f_reweight->IsValid() || f_reweight->GetNpar() != static_cast<int>(c.trw_params.size()))
            throw std::runtime_error("invalid response prior formula or parameter count");
        for (size_t i = 0; i < c.trw_params.size(); ++i) f_reweight->SetParameter(i, c.trw_params[i]);
    }

    // PJ bookkeeping / QA
    TH1D *h_cutflow = new TH1D("h_pj_cutflow", "PJ cutflow", 12, 0, 12);
    const char *cf_labels[12] = {"events read", "events selected", "photons read", "photons ET>=min", "photons common",
                                 "photons tight", "photons A", "photons B", "photons C", "photons D",
                                 "truth fiducial", "truth matched"};
    for (int i = 0; i < 12; ++i) h_cutflow->GetXaxis()->SetBinLabel(i + 1, cf_labels[i]);
    TH1D *h_flag = new TH1D("h_pj_flagcheck", "stored-flag cross-check (data, formula mode)", 4, 0, 4);
    const char *fl_labels[4] = {"photons checked", "tight mismatch", "nontight mismatch", "iso mismatch"};
    for (int i = 0; i < 4; ++i) h_flag->GetXaxis()->SetBinLabel(i + 1, fl_labels[i]);
    TH1D *h_vz = new TH1D("h_pj_vertex_z", "selected events vertex z", 200, -100, 100);
    TH1D *h_w = new TH1D("h_pj_event_weight_log10", "log10(event weight) selected events", 200, -4, 8);
    TH1D *h_cent = new TH1D("h_pj_centrality", "centrality of selected events", 101, -1.5, 99.5);
    const double run_lo = (c.system == "auau") ? 70000 : 47000, run_hi = (c.system == "auau") ? 80000 : 55000;
    TH1D *h_run = new TH1D("h_pj_run_events", "selected events per run", (int)(run_hi - run_lo), run_lo, run_hi);
    TH2D *h_iso_et = new TH2D("h_pj_iso_vs_et", "iso (used) vs ET, common photons", 80, 0, 40, 200, -20, 40);
    TH2D *h_score_et = new TH2D("h_pj_score_vs_et", "BDT score vs ET, common photons", 80, 0, 40, 100, 0, 1);
    TH1D *h_anom = new TH1D("h_pj_auau_threshold_recipe", "Au+Au iso threshold recipe (0 table, 1 small)", 2, 0, 2);

    // ---- event loop ----------------------------------------------------------
    TRandom3 rng(c.random_seed);
    std::vector<PJ::Event> events;
    std::vector<PJ::Photon> photons;
    std::vector<PJ::TruthPhoton> truths;
    std::unordered_map<PJ::EventKey, double, PJ::EventKeyHash> leadjet;
    std::unordered_map<PJ::LinkKey, Int_t, PJ::LinkKeyHash> links;
    const bool use_links = issim && (c.match_source == "links");
    std::cout << "[PJ] match_source=" << c.match_source << std::endl;
    long n_parts_done = 0;
    long n_below_floor = 0;

    for (const auto &pe : parts)
    {
        auto fin = PJ::OpenRoot(pe.path);
        PJ::LoadEvents(fin.get(), events, !issim && c.system == "auau");
        PJ::LoadPhotons(fin.get(), photons);
        std::string sample;
        PPG12::SampleConfig sc;
        if (issim)
        {
            PJ::LoadTruthPhotons(fin.get(), truths);
            if (use_links) PJ::LoadPhotonTruthLinks(fin.get(), links);
            if (need_sample)
            {
                sample = PJ::SampleForPart(sample_map, pe.index);
                sc = PPG12::GetSampleConfig(sample);
                if (!sc.valid)
                    throw std::runtime_error(Form("part index %d (%s) has no valid sample in the sample map", pe.index, pe.path.c_str()));
                if (sc.isbackground && c.weight_mode == "sample_map") PJ::LoadLeadingTruthJetPt(fin.get(), leadjet);
            }
        }
        YAML::Node input;
        input["index"] = pe.index;
        input["path"] = pe.path;
        input["ROOT_UUID"] = fin->GetUUID().AsString();
        input["events"] = events.size();
        input["photons"] = photons.size();
        input["truth_photons"] = truths.size();
        provenance["inputs"].push_back(input);
        fin->Close();

        const auto ph_by_ev = PJ::GroupByEvent(photons);
        const auto tr_by_ev = PJ::GroupByEvent(truths);
        PJ::ValidateEventKeys(events, ph_by_ev, tr_by_ev);
        const auto photon_truth_row = issim
            ? PJ::MatchPhotons(photons, truths, ph_by_ev, tr_by_ev, links, use_links)
            : std::vector<int>(photons.size(), -1);

        for (const auto &ev : events)
        {
            h_cutflow->Fill(0.5);
            double weight = 1.0;
            if (!issim)
            {
                if (!PJ::PassEventData(ev, c)) continue;
            }
            else if (!PJ::SimulationWeight(ev, c, sc, vtxw, truths, tr_by_ev, leadjet, weight)) continue;
            h_cutflow->Fill(1.5);
            h_vz->Fill(ev.vertex_z);
            h_w->Fill(std::log10(std::max(weight, 1e-12)));
            h_cent->Fill(ev.centrality);
            h_run->Fill(ev.run);

            // ---- photons ---------------------------------------------------
            struct Cand { int idx; PJ::Decision d; double w; };
            std::vector<Cand> cands;
            auto pit = ph_by_ev.find(ev.key());
            if (pit != ph_by_ev.end())
            {
                for (int ip : pit->second)
                {
                    const PJ::Photon &p = photons[ip];
                    h_cutflow->Fill(2.5);
                    const double etmod = PJ::CalibratedET(c, p.et, issim);
                    if (!std::isfinite(etmod)) throw std::runtime_error("nonfinite calibrated photon ET");
                    if (!issim && c.min_photon_et > 0 && p.et < c.min_photon_et) ++n_below_floor; // data skim floor check only
                    if (etmod < c.reco_min_ET) continue;
                    if (issim && need_sample && sc.isbackground && etmod > sc.cluster_ET_upper) continue;
                    if (!(p.eta > c.eta_bins[0] && p.eta < c.eta_bins[1])) continue;
                    h_cutflow->Fill(3.5);
                    PJ::Decision d = PJ::Classify(p, c, imap, issim, etmod);
                    const double w = issim ? weight : weight * PJ::TriggerEffWeight(c, etmod);
                    if (c.system == "auau")
                    {
                        // small-threshold recipe seen in embedded MC (t = 0.081 + 0.019 ET, gap 0.8)
                        const bool small = std::fabs(p.iso4_thr - (0.081 + 0.019 * p.et)) < 1e-6;
                        h_anom->Fill(small ? 1.5 : 0.5);
                    }
                    if (!issim && c.threshold_source == "formula" && c.flag_check != "off")
                    {
                        h_flag->Fill(0.5);
                        if (d.flag_mismatch_tight) h_flag->Fill(1.5);
                        if (d.flag_mismatch_nontight) h_flag->Fill(2.5);
                        if (d.flag_mismatch_iso) h_flag->Fill(3.5);
                    }
                    h_all->Fill(d.et, w);
                    if (d.common)
                    {
                        h_cutflow->Fill(4.5);
                        h_common->Fill(d.et, w);
                        h_iso_et->Fill(d.et, d.iso, w);
                        h_score_et->Fill(d.et, p.score, w);
                    }
                    if (d.tight) { h_cutflow->Fill(5.5); h_tight->Fill(d.et, w); }
                    if (d.region >= 0)
                    {
                        h_cutflow->Fill(6.5 + d.region);
                        h_abcd[d.region]->Fill(d.et, w);
                    }
                    cands.push_back({ip, d, w});
                }
            }

            if (!issim) continue;

            // ---- truth loop (signal templates, efficiencies, response) -----
            std::vector<char> matched_to_fiducial(cands.size(), 0);
            auto tit = tr_by_ev.find(ev.key());
            if (tit != tr_by_ev.end())
            {
                for (int ti : tit->second)
                {
                    const PJ::TruthPhoton &t = truths[ti];
                    if (!(t.eta > c.eta_bins[0] && t.eta < c.eta_bins[1])) continue;
                    if (!PJ::IsFiducialTruth(t, c)) continue;
                    h_cutflow->Fill(10.5);
                    h_truth->Fill(t.pt, weight);
                    h_truth_novtx->Fill(t.pt, weight);
                    h_truth_vc->Fill(t.pt, weight);
                    h_truth_vc_mbd->Fill(t.pt, weight);
                    h_truth_vc_n->Fill(t.pt, weight);
                    h_truth_vc_s->Fill(t.pt, weight);

                    bool photon_reco = false, photon_iso = false, photon_tight_iso = false;
                    for (size_t ic = 0; ic < cands.size(); ++ic)
                    {
                        if (photon_truth_row[cands[ic].idx] != ti) continue;
                        matched_to_fiducial[ic] = 1;
                        photon_reco = true;
                        const PJ::Decision &d = cands[ic].d;
                        if (d.iso_pass) photon_iso = true;
                        if (d.tight && d.iso_pass) photon_tight_iso = true;
                        const bool in_range = (t.pt > pTmin_truth && t.pt < pTmax_truth && d.et > pTmin && d.et < pTmax);
                        if (in_range)
                        {
                            h_all_sig->Fill(d.et, weight);
                            if (d.tight) h_tight_sig->Fill(d.et, weight);
                            if (d.region >= 0) h_sig[d.region]->Fill(d.et, weight);
                        }
                        if (is_signal && d.tight && d.iso_pass && t.pt > pTmin_truth && t.pt < pTmax_truth)
                        {
                            const double response_et = PJ::ResponseET(c, d.et, t.pt, rng);
                            if (!(response_et > pTmin && response_et < pTmax)) continue;
                            double rw = 1.0;
                            if (c.unfold_reweight)
                                rw = f_reweight->Eval(std::min(c.trw_clamp_max, std::max(c.trw_clamp_min, t.pt)));
                            if (!std::isfinite(rw) || rw <= 0) throw std::runtime_error("invalid response prior weight");
                            h_pT_truth_resp->Fill(t.pt, weight * rw);
                            h_pT_reco_resp->Fill(response_et, weight * rw);
                            response->Fill(response_et, t.pt, weight * rw);
                            h_resp_full->Fill(response_et, t.pt, weight * rw);
                        }
                    }
                    if (photon_reco) h_cutflow->Fill(11.5);
                    eff_reco->FillWeighted(photon_reco, weight, t.pt);
                    eff_all->FillWeighted(photon_reco && photon_iso && photon_tight_iso, weight, t.pt);
                    if (photon_reco)
                    {
                        eff_iso->FillWeighted(photon_iso, weight, t.pt);
                        if (photon_iso) eff_id->FillWeighted(photon_tight_iso, weight, t.pt);
                    }
                }
            }
            for (size_t ic = 0; ic < cands.size(); ++ic)
            {
                if (matched_to_fiducial[ic]) continue;
                const PJ::Decision &d = cands[ic].d;
                if (d.region >= 0) h_nm[d.region]->Fill(d.et, cands[ic].w);
            }
        }
        ++n_parts_done;
        if (n_parts_done % 20 == 0) std::cout << "[PJ] parts done: " << n_parts_done << "/" << parts.size() << std::endl;
    }

    // ---- summary -----------------------------------------------------------
    std::cout << "[PJ] cutflow:";
    for (int i = 1; i <= 12; ++i) std::cout << " " << cf_labels[i - 1] << "=" << (long long)h_cutflow->GetBinContent(i);
    std::cout << std::endl;
    if (!issim && c.threshold_source == "formula" && c.flag_check != "off")
    {
        std::cout << "[PJ] flag check: photons=" << (long long)h_flag->GetBinContent(1)
                  << " tight mismatches=" << (long long)h_flag->GetBinContent(2)
                  << " nontight mismatches=" << (long long)h_flag->GetBinContent(3)
                  << " iso mismatches=" << (long long)h_flag->GetBinContent(4) << std::endl;
    }
    if (c.system == "auau")
        std::cout << "[PJ] Au+Au iso-threshold recipe: table=" << (long long)h_anom->GetBinContent(1)
                  << " small(0.081+0.019ET)=" << (long long)h_anom->GetBinContent(2) << std::endl;
    if (n_below_floor > 0)
        std::cout << "[PJ] WARNING: " << n_below_floor << " data photons below photonjet.min_photon_et=" << c.min_photon_et
                  << " (expected 0 for the PhotonJetTrees skim)" << std::endl;
    std::cout << "[PJ] A=" << h_tight_iso->Integral() << " B=" << h_tight_noniso->Integral()
              << " C=" << h_nontight_iso->Integral() << " D=" << h_nontight_noniso->Integral() << std::endl;

    if (!issim && c.threshold_source == "formula" && c.flag_check == "strict" &&
        h_flag->GetBinContent(2) + h_flag->GetBinContent(3) + h_flag->GetBinContent(4) > 0)
        throw std::runtime_error("stored flags disagree with configured thresholds; see h_pj_flagcheck counts above");
    if (n_below_floor > 0)
        throw std::runtime_error("data violate the configured PhotonJetTrees minimum photon ET");

    // Close both complete files before publishing either final filename.
    PJ::WriteRunMetadata(fout, config, provenance, parts, n_parts_done);
    if (fresp)
    {
        fresp->cd();
        if (response->Write() <= 0) throw std::runtime_error("cannot write RooUnfold response");
        PJ::WriteRunMetadata(fresp, config, provenance, parts, n_parts_done);
        response_file->Close();
    }
    main_file.Close();
    if (response_file) response_file->Commit();
    main_file.Commit();
    std::cout << "[PJ] wrote " << names.main << (fresp ? (" and " + names.response) : "") << std::endl;
}

void PhotonJetHistMaker(const std::string &configname = "",
                        const std::string &product = "data",
                        const std::string &chunk_list = "",
                        const std::string &chunk_tag = "")
{
    try { MakeHistograms(configname, product, chunk_list, chunk_tag); }
    catch (const std::exception &error)
    {
        std::cerr << "[PJ] ERROR: " << error.what() << std::endl;
        gSystem->Exit(1);
    }
}
