// Typed histogram ownership and booking. Physics fills stay in the maker.
#pragma once
#include "PhotonJetConfig.h"
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <RooUnfoldResponse.h>
#include <array>
#include <memory>

namespace PJ {
using H1Ptr = std::unique_ptr<TH1D>;
using H2Ptr = std::unique_ptr<TH2D>;
using EffPtr = std::unique_ptr<TEfficiency>;
struct RecoHistograms {
    H1Ptr all, common, tight, signal_all, signal_tight;
    std::array<H1Ptr, 4> abcd, signal, unmatched;
};
struct PhotonQA { H2Ptr isolation, score; };
struct TruthHistograms {
    H1Ptr spectrum, novtx, vertexcut, mbd, north, south, only_north, only_south, neither;
    EffPtr reco, iso, id, all, converts;
};
struct ResponseHistograms {
    H1Ptr reco, truth;
    H2Ptr matrix;
    std::unique_ptr<RooUnfoldResponse> unfold;
};
struct HistogramSet {
    RecoHistograms reco;
    PhotonQA qa;
    TruthHistograms truth;
    ResponseHistograms response;
};
using HistogramGrid = std::vector<std::vector<HistogramSet>>;
struct JobQA { H1Ptr cutflow, flags, vertex, weight, centrality, run, recipe; };

// One inventory home for publication and Python merge checks. An observable
// addition updates this per-cell count beside its member, booking and writing.
inline int ExpectedRootObjectCount(int cells, bool response_file, bool imported_mbd = false)
{
    constexpr int metadata = 4; // config, provenance, input manifest, input-parts QA
    constexpr int job_qa = 7, main_per_cell = 33, response_per_cell = 4;
    if (response_file) return metadata + response_per_cell * cells;
    return metadata + job_qa + main_per_cell * cells + (imported_mbd ? 2 : 0);
}

inline H1Ptr BookSpectrum(const std::string &name, const std::string &title, const std::vector<double> &edges)
{
    auto h = std::make_unique<TH1D>(name.c_str(), title.c_str(), static_cast<int>(edges.size()) - 1, edges.data());
    h->SetDirectory(nullptr); h->Sumw2();
    return h;
}
inline EffPtr BookEfficiency(const std::string &name, const std::string &title, const std::vector<double> &edges)
{
    auto h = std::make_unique<TEfficiency>(name.c_str(), title.c_str(), static_cast<int>(edges.size()) - 1, edges.data());
    h->SetDirectory(nullptr); h->SetStatisticOption(TEfficiency::kBUniform); h->SetUseWeightedEvents();
    return h;
}
inline JobQA BookJobQA(const Cuts &c)
{
    TDirectory::TContext detached(nullptr);
    JobQA q;
    q.cutflow = std::make_unique<TH1D>("h_pj_cutflow", "PJ cutflow", 12, 0, 12);
    const char *labels[] = {"events read", "events selected", "photons read", "photons ET>=min", "photons common", "photons tight",
                            "photons A", "photons B", "photons C", "photons D", "truth fiducial", "truth matched"};
    for (int i = 0; i < 12; ++i) q.cutflow->GetXaxis()->SetBinLabel(i + 1, labels[i]);
    q.flags = std::make_unique<TH1D>("h_pj_flagcheck", "stored-flag cross-check (data, formula mode)", 4, 0, 4);
    const char *flags[] = {"photons checked", "tight mismatch", "nontight mismatch", "iso mismatch"};
    for (int i = 0; i < 4; ++i) q.flags->GetXaxis()->SetBinLabel(i + 1, flags[i]);
    q.vertex = std::make_unique<TH1D>("h_pj_vertex_z", "selected events vertex z", 200, -100, 100);
    q.weight = std::make_unique<TH1D>("h_pj_event_weight_log10", "log10(event weight) selected events", 200, -4, 8);
    q.centrality = std::make_unique<TH1D>("h_pj_centrality", "centrality of selected events", 101, -1.5, 99.5);
    const double lo = c.system == "auau" ? 70000 : 47000, hi = c.system == "auau" ? 80000 : 55000;
    q.run = std::make_unique<TH1D>("h_pj_run_events", "selected events per run", static_cast<int>(hi - lo), lo, hi);
    q.recipe = std::make_unique<TH1D>("h_pj_auau_threshold_recipe", "Au+Au iso threshold recipe (0 table, 1 small)", 2, 0, 2);
    return q;
}
inline HistogramSet BookHistogramSet(const Cuts &c, const BinInfo &bin, bool signal)
{
    TDirectory::TContext detached(nullptr);
    HistogramSet h;
    const std::string title = " " + bin.Title();
    auto reco = [&](const std::string &base, const std::string &label) {
        return BookSpectrum(bin.Name(base), label + title, c.pT_bins);
    };
    auto truth = [&](const std::string &base, const std::string &label = "Truth pT") {
        return BookSpectrum(bin.Name(base), label + title, c.pT_bins_truth);
    };
    auto efficiency = [&](const std::string &base, const std::string &label) {
        return BookEfficiency(bin.Name(base, "_eta_0"), label + title, c.pT_bins_truth);
    };
    h.reco.all = reco("h_all_cluster", "All Cluster");
    h.reco.common = reco("h_common_cluster", "Common Cluster");
    h.reco.tight = reco("h_tight_cluster", "Tight Cluster");
    h.reco.signal_all = reco("h_all_cluster_signal", "All Cluster");
    h.reco.signal_tight = reco("h_tight_cluster_signal", "Tight Cluster");
    const char *regions[] = {"tight_iso", "tight_noniso", "nontight_iso", "nontight_noniso"};
    const char *labels[] = {"Tight Iso Cluster", "Tight Non-Iso Cluster", "Non-Tight Iso Cluster", "Non-Tight Non-Iso Cluster"};
    for (int r = 0; r < 4; ++r) {
        const std::string base = std::string("h_") + regions[r] + "_cluster";
        h.reco.abcd[r] = reco(base, labels[r]);
        h.reco.signal[r] = reco(base + "_signal", labels[r]);
        h.reco.unmatched[r] = reco(base + "_notmatch", labels[r]);
    }
    const std::string qa_title = bin.single ? "" : title;
    h.qa.isolation = std::make_unique<TH2D>(bin.Name("h_pj_iso_vs_et", "").c_str(),
        ("iso (used) vs ET, common photons" + qa_title).c_str(), 80, 0, 40, 200, -20, 40);
    h.qa.score = std::make_unique<TH2D>(bin.Name("h_pj_score_vs_et", "").c_str(),
        ("BDT score vs ET, common photons" + qa_title).c_str(), 80, 0, 40, 100, 0, 1);
    h.truth.spectrum = truth("h_truth_pT");
    h.truth.novtx = truth("h_truth_pT_novtx", "Truth pT (no vtx wt)");
    h.truth.vertexcut = truth("h_truth_pT_vertexcut");
    h.truth.mbd = truth("h_truth_pT_vertexcut_mbd_cut");
    h.truth.north = truth("h_truth_pT_vertexcut_mbd_north_cut");
    h.truth.south = truth("h_truth_pT_vertexcut_mbd_south_cut");
    h.truth.only_north = truth("h_truth_pT_vertexcut_mbd_only_north");
    h.truth.only_south = truth("h_truth_pT_vertexcut_mbd_only_south");
    h.truth.neither = truth("h_truth_pT_vertexcut_mbd_neither");
    h.truth.reco = efficiency("eff_reco", "Reco Efficiency");
    h.truth.iso = efficiency("eff_iso", "Iso Efficiency");
    h.truth.id = efficiency("eff_id", "ID Efficiency");
    h.truth.all = efficiency("eff_all", "All Efficiency");
    h.truth.converts = efficiency("eff_converts", "Conversion Prob");
    if (signal) {
        h.response.reco = reco("h_pT_reco_response", "Reco pT");
        h.response.truth = truth("h_pT_truth_response");
        h.response.matrix = std::make_unique<TH2D>(bin.Name("h_response_full").c_str(), ("Response Matrix" + title).c_str(),
            static_cast<int>(c.pT_bins.size()) - 1, c.pT_bins.data(),
            static_cast<int>(c.pT_bins_truth.size()) - 1, c.pT_bins_truth.data());
        h.response.matrix->Sumw2();
        h.response.unfold = std::make_unique<RooUnfoldResponse>(h.response.reco.get(), h.response.truth.get(),
            h.response.matrix.get(), bin.Name("response_matrix_full").c_str(), bin.single ? "" : bin.Title().c_str(), false);
    }
    return h;
}
inline HistogramGrid BookHistograms(const Cuts &c, const BinLayout &bins, bool signal)
{
    HistogramGrid histograms(bins.nCentrality());
    for (int centrality = 0; centrality < bins.nCentrality(); ++centrality) {
        auto &row = histograms[centrality];
        row.reserve(bins.nEta());
        for (int eta = 0; eta < bins.nEta(); ++eta)
            row.push_back(BookHistogramSet(c, bins.cell(centrality, eta), signal));
    }
    return histograms;
}
inline void WriteHistogramObject(TObject &object)
{
    if (object.Write(object.GetName(), TObject::kOverwrite) <= 0)
        throw std::runtime_error("cannot write " + std::string(object.GetName()));
}
inline void WriteHistograms(TFile &main, TFile *response, const JobQA &q, const HistogramGrid &histograms)
{
    TDirectory::TContext directory(&main);
    for (auto *h : {q.cutflow.get(), q.flags.get(), q.vertex.get(), q.weight.get(), q.centrality.get(), q.run.get(), q.recipe.get()})
        WriteHistogramObject(*h);
    for (const auto &row : histograms) for (const auto &h : row) {
        main.cd();
        for (auto *spectrum : {h.reco.all.get(), h.reco.common.get(), h.reco.tight.get(), h.reco.signal_all.get(), h.reco.signal_tight.get(),
                              h.truth.spectrum.get(), h.truth.novtx.get(), h.truth.vertexcut.get(), h.truth.mbd.get(), h.truth.north.get(),
                              h.truth.south.get(), h.truth.only_north.get(), h.truth.only_south.get(), h.truth.neither.get()})
            WriteHistogramObject(*spectrum);
        for (int region = 0; region < 4; ++region)
            for (auto *spectrum : {h.reco.abcd[region].get(), h.reco.signal[region].get(), h.reco.unmatched[region].get()})
                WriteHistogramObject(*spectrum);
        WriteHistogramObject(*h.qa.isolation); WriteHistogramObject(*h.qa.score);
        for (auto *efficiency : {h.truth.reco.get(), h.truth.iso.get(), h.truth.id.get(), h.truth.all.get(), h.truth.converts.get()})
            WriteHistogramObject(*efficiency);
        if (response) {
            response->cd();
            WriteHistogramObject(*h.response.reco); WriteHistogramObject(*h.response.truth);
            WriteHistogramObject(*h.response.matrix); WriteHistogramObject(*h.response.unfold);
        }
    }
}

// ROOT's empty-response constructor contributes projection entries. Rebuild once
// after merging so those offsets occur once per complete sample, independently per cell.
inline std::unique_ptr<RooUnfoldResponse> FinalizeResponse(TH1D &measured, TH1D &truth, TH2D &matrix,
    const std::string &name = "response_matrix_full_0", const std::string &title = "")
{
    for (TH1 *h : {static_cast<TH1 *>(&measured), static_cast<TH1 *>(&truth), static_cast<TH1 *>(&matrix)})
        if (!h->GetSumw2N()) h->Sumw2();
    TDirectory::TContext detached(nullptr);
    TH1D empty_m(measured), empty_t(truth); TH2D empty_r(matrix);
    empty_m.Reset(); empty_t.Reset(); empty_r.Reset();
    auto result = std::make_unique<RooUnfoldResponse>(&empty_m, &empty_t, &empty_r, name.c_str(), title.c_str(), false);
    const double measured_setup = result->Hmeasured()->GetEntries(), truth_setup = result->Htruth()->GetEntries();
    measured.Copy(*result->Hmeasured()); truth.Copy(*result->Htruth()); matrix.Copy(*result->Hresponse());
    result->Hmeasured()->SetEntries(measured.GetEntries() + measured_setup);
    result->Htruth()->SetEntries(truth.GetEntries() + truth_setup);
    result->ClearCache();
    return result;
}
inline void FinalizeMergedResponses(TFile &file, const Cuts &cuts)
{
    TDirectory::TContext directory(&file);
    const BinLayout bins(cuts);
    for (int centrality = 0; centrality < bins.nCentrality(); ++centrality) for (int eta = 0; eta < bins.nEta(); ++eta) {
        const auto bin = bins.cell(centrality, eta);
        auto *measured = dynamic_cast<TH1D *>(file.Get(bin.Name("h_pT_reco_response").c_str()));
        auto *truth = dynamic_cast<TH1D *>(file.Get(bin.Name("h_pT_truth_response").c_str()));
        auto *matrix = dynamic_cast<TH2D *>(file.Get(bin.Name("h_response_full").c_str()));
        if (!measured || !truth || !matrix) throw std::runtime_error("missing response histograms for cell " + bin.Name(""));
        const auto name = bin.Name("response_matrix_full");
        auto result = FinalizeResponse(*measured, *truth, *matrix, name, bin.single ? "" : bin.Title());
        file.Delete((name + ";*").c_str());
        WriteHistogramObject(*result);
    }
}
} // namespace PJ
