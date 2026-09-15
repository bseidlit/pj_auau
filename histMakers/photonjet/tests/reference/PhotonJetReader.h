// PhotonJetTrees_v1 rows and checked schema loading. Physics choices live in
// PhotonJetConfig.h / PhotonJetSelection.h / PhotonJetWeights.h.
#pragma once
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace PJ
{
// ---------------------------------------------------------------------------
// Event identity within one part: (source_file_index, event_id_lo).
// event_id_hi is the constant 2 in the whole release, it is carried but not
// used as a key.
// ---------------------------------------------------------------------------
struct EventKey
{
    Int_t src = -1;
    ULong64_t lo = 0;
    bool operator==(const EventKey &o) const { return src == o.src && lo == o.lo; }
};
struct EventKeyHash
{
    size_t operator()(const EventKey &k) const
    {
        return std::hash<ULong64_t>()(k.lo * 1000003ULL + (ULong64_t)(k.src + 7));
    }
};

// ---------------------------------------------------------------------------
// Row structs (one per tree). Leaf types of the release: D, I, L, l only.
// ---------------------------------------------------------------------------
struct Event
{
    Int_t src = -1;
    Long64_t source_entry = -1;
    ULong64_t id_hi = 0, id_lo = 0;
    Int_t run = -1;
    ULong64_t trigger_bits = 0, live_trigger_bits = 0, scaled_trigger_bits = 0;
    Int_t scaled_bit30 = -1;
    Double_t vertex_z = 0, centrality = -1, event_weight = 1, total_calo_energy = 0;
    Int_t terminal_status = -1;
    // Au+Au data only (bound when present, defaults otherwise)
    Int_t scaled_bit22 = -1, minimum_bias_pass = -1, nominal_event_selection_pass = -1;
    EventKey key() const { return EventKey{src, id_lo}; }
};

struct Photon
{
    Int_t src = -1;
    ULong64_t id_lo = 0;
    Int_t run = -1, scaled_bit30 = -1, terminal_status = -1;
    Double_t vertex_z = 0, centrality = -1, event_weight = 1;
    Int_t scaled_bit22 = -1, minimum_bias_pass = -1;
    Int_t ordinal = -1;
    Double_t et = 0, eta = 0, phi = 0;
    Double_t score = 0, thr_tight = 0, thr_nt_low = 0, thr_nt_high = 0;
    Int_t is_tight = 0, is_nontight = 0, is_not_tight = 0, input_count = 0;
    Double_t iso3 = 0, iso3_thr = 0, iso3_nonthr = 0;
    Double_t iso4 = 0, iso4_thr = 0, iso4_nonthr = 0;
    Int_t iso3_pass = 0, iso4_pass = 0;
    Int_t truth_matched = -1, truth_barcode = -1;
    Double_t native_weta = 0, native_wphi = 0, native_weta33 = 0, native_wphi33 = 0;
    Double_t native_e11e33 = 0, native_e32e35 = 0, native_et1 = 0;
    Double_t in[14] = {0};
    EventKey key() const { return EventKey{src, id_lo}; }
};

struct TruthPhoton
{
    Int_t src = -1;
    ULong64_t id_lo = 0;
    Double_t pt = 0, eta = 0, phi = 0, iso = 0;
    Int_t prompt_class = -1, source_role = -1, barcode = -1;
    EventKey key() const { return EventKey{src, id_lo}; }
};

// ---------------------------------------------------------------------------
// Checked branch binding. Any missing branch or leaf-type mismatch is fatal.
// ---------------------------------------------------------------------------
template <class T>
inline void Bind(TTree *t, const char *name, T *addr, std::vector<std::string> &problems)
{
    if (!t->GetBranch(name))
    {
        problems.push_back(std::string("missing branch: ") + name);
        return;
    }
    t->SetBranchStatus(name, 1);
    const Int_t rc = t->SetBranchAddress(name, addr);
    if (rc < 0)
    {
        problems.push_back(std::string("type mismatch (rc=") + std::to_string(rc) + ") on branch: " + name);
    }
}

template <class T>
inline bool BindOptional(TTree *t, const char *name, T *addr, std::vector<std::string> &problems)
{
    if (!t->GetBranch(name)) return false;
    Bind(t, name, addr, problems);
    return true;
}

inline void AbortOnProblems(const std::vector<std::string> &problems, const std::string &tree, const std::string &file)
{
    if (problems.empty()) return;
    std::cerr << "[PhotonJetReader] FATAL: schema problems in tree '" << tree << "' of " << file << "\n";
    for (const auto &p : problems) std::cerr << "    " << p << "\n";
    throw std::runtime_error("PhotonJetReader schema mismatch");
}

inline TTree *GetTree(TFile *f, const char *name)
{
    TTree *t = dynamic_cast<TTree *>(f->Get(name));
    if (!t)
    {
        std::cerr << "[PhotonJetReader] FATAL: tree '" << name << "' not found in " << f->GetName() << std::endl;
        throw std::runtime_error("PhotonJetReader missing tree");
    }
    t->SetBranchStatus("*", 0);
    return t;
}

// ---------------------------------------------------------------------------
// Loaders: read a whole per-part tree into memory (parts are small: <= 14k
// events, <= 80k photons, <= 10k truth photons).
// ---------------------------------------------------------------------------
inline void LoadEvents(TFile *f, std::vector<Event> &out, bool require_mb = false)
{
    TTree *t = GetTree(f, "events");
    Event e;
    std::vector<std::string> pr;
    Bind(t, "source_file_index", &e.src, pr);
    Bind(t, "source_entry", &e.source_entry, pr);
    Bind(t, "event_id_hi", &e.id_hi, pr);
    Bind(t, "event_id_lo", &e.id_lo, pr);
    Bind(t, "run", &e.run, pr);
    Bind(t, "trigger_bits", &e.trigger_bits, pr);
    Bind(t, "live_trigger_bits", &e.live_trigger_bits, pr);
    Bind(t, "scaled_trigger_bits", &e.scaled_trigger_bits, pr);
    Bind(t, "scaled_bit30", &e.scaled_bit30, pr);
    Bind(t, "vertex_z", &e.vertex_z, pr);
    Bind(t, "centrality", &e.centrality, pr);
    Bind(t, "event_weight", &e.event_weight, pr);
    Bind(t, "total_calo_energy", &e.total_calo_energy, pr);
    Bind(t, "terminal_status", &e.terminal_status, pr);
    if (require_mb)
    {
        Bind(t, "scaled_bit22", &e.scaled_bit22, pr);
        Bind(t, "minimum_bias_pass", &e.minimum_bias_pass, pr);
    }
    else
    {
        BindOptional(t, "scaled_bit22", &e.scaled_bit22, pr);
        BindOptional(t, "minimum_bias_pass", &e.minimum_bias_pass, pr);
    }
    BindOptional(t, "nominal_event_selection_pass", &e.nominal_event_selection_pass, pr);
    AbortOnProblems(pr, "events", f->GetName());
    const Long64_t n = t->GetEntries();
    out.clear();
    out.reserve(n);
    for (Long64_t i = 0; i < n; ++i)
    {
        if (t->GetEntry(i) <= 0)
            throw std::runtime_error(std::string("failed reading ") + t->GetName() + " entry " + std::to_string(i) + " in " + f->GetName());
        out.push_back(e);
    }
}

inline void LoadPhotons(TFile *f, std::vector<Photon> &out)
{
    TTree *t = GetTree(f, "photons");
    Photon p;
    std::vector<std::string> pr;
    Bind(t, "source_file_index", &p.src, pr);
    Bind(t, "event_id_lo", &p.id_lo, pr);
    Bind(t, "run", &p.run, pr);
    Bind(t, "scaled_bit30", &p.scaled_bit30, pr);
    Bind(t, "terminal_status", &p.terminal_status, pr);
    Bind(t, "vertex_z", &p.vertex_z, pr);
    Bind(t, "centrality", &p.centrality, pr);
    Bind(t, "event_weight", &p.event_weight, pr);
    BindOptional(t, "scaled_bit22", &p.scaled_bit22, pr);
    BindOptional(t, "minimum_bias_pass", &p.minimum_bias_pass, pr);
    Bind(t, "photon_encounter_ordinal", &p.ordinal, pr);
    Bind(t, "photon_et", &p.et, pr);
    Bind(t, "photon_eta", &p.eta, pr);
    Bind(t, "photon_phi", &p.phi, pr);
    Bind(t, "bdt_score", &p.score, pr);
    Bind(t, "bdt_tight_threshold", &p.thr_tight, pr);
    Bind(t, "bdt_nontight_low_threshold", &p.thr_nt_low, pr);
    Bind(t, "bdt_nontight_high_threshold", &p.thr_nt_high, pr);
    Bind(t, "bdt_is_tight", &p.is_tight, pr);
    Bind(t, "bdt_is_nontight", &p.is_nontight, pr);
    Bind(t, "bdt_is_not_tight", &p.is_not_tight, pr);
    Bind(t, "bdt_input_count", &p.input_count, pr);
    Bind(t, "iso_r03", &p.iso3, pr);
    Bind(t, "iso_r03_threshold", &p.iso3_thr, pr);
    Bind(t, "iso_r03_nonisolated_threshold", &p.iso3_nonthr, pr);
    Bind(t, "iso_r03_pass", &p.iso3_pass, pr);
    Bind(t, "iso_r04", &p.iso4, pr);
    Bind(t, "iso_r04_threshold", &p.iso4_thr, pr);
    Bind(t, "iso_r04_nonisolated_threshold", &p.iso4_nonthr, pr);
    Bind(t, "iso_r04_pass", &p.iso4_pass, pr);
    Bind(t, "truth_matched", &p.truth_matched, pr);
    Bind(t, "truth_barcode", &p.truth_barcode, pr);
    Bind(t, "native_weta_cogx", &p.native_weta, pr);
    Bind(t, "native_wphi_cogx", &p.native_wphi, pr);
    Bind(t, "native_weta33_cogx", &p.native_weta33, pr);
    Bind(t, "native_wphi33_cogx", &p.native_wphi33, pr);
    Bind(t, "native_e11_over_e33", &p.native_e11e33, pr);
    Bind(t, "native_e32_over_e35", &p.native_e32e35, pr);
    Bind(t, "native_et1", &p.native_et1, pr);
    for (int k = 0; k < 14; ++k)
    {
        Bind(t, Form("bdt_input_%02d", k), &p.in[k], pr);
    }
    AbortOnProblems(pr, "photons", f->GetName());
    const Long64_t n = t->GetEntries();
    out.clear();
    out.reserve(n);
    for (Long64_t i = 0; i < n; ++i)
    {
        if (t->GetEntry(i) <= 0)
            throw std::runtime_error(std::string("failed reading ") + t->GetName() + " entry " + std::to_string(i) + " in " + f->GetName());
        out.push_back(p);
    }
}

inline void LoadTruthPhotons(TFile *f, std::vector<TruthPhoton> &out)
{
    TTree *t = GetTree(f, "truthPhotons");
    TruthPhoton tp;
    std::vector<std::string> pr;
    Bind(t, "source_file_index", &tp.src, pr);
    Bind(t, "event_id_lo", &tp.id_lo, pr);
    Bind(t, "truth_photon_pt", &tp.pt, pr);
    Bind(t, "truth_photon_eta", &tp.eta, pr);
    Bind(t, "truth_photon_phi", &tp.phi, pr);
    Bind(t, "prompt_class", &tp.prompt_class, pr);
    Bind(t, "source_role", &tp.source_role, pr);
    Bind(t, "generator_barcode", &tp.barcode, pr);
    Bind(t, "truth_isolation", &tp.iso, pr);
    AbortOnProblems(pr, "truthPhotons", f->GetName());
    const Long64_t n = t->GetEntries();
    out.clear();
    out.reserve(n);
    for (Long64_t i = 0; i < n; ++i)
    {
        if (t->GetEntry(i) <= 0)
            throw std::runtime_error(std::string("failed reading ") + t->GetName() + " entry " + std::to_string(i) + " in " + f->GetName());
        out.push_back(tp);
    }
}

// Photon -> truth-photon links from recoTruthLinks: class (1,1,0) rows only.
// Key = (source_file_index, event_id_lo, reco_index), value = truth_index, both
// event-local positions in file order of the `photons` / `truthPhotons` trees
// (verified: candidate_id_lo of photons[reco_index] == reco_id_lo of the link).
// Needed for Au+Au, where photons.truth_barcode is -1 for ~80 % of linked photons.
struct LinkKey
{
    Int_t src = -1;
    ULong64_t lo = 0;
    Int_t reco = -1;
    bool operator==(const LinkKey &o) const { return src == o.src && lo == o.lo && reco == o.reco; }
};
struct LinkKeyHash
{
    size_t operator()(const LinkKey &k) const
    {
        return std::hash<ULong64_t>()(k.lo * 1000003ULL + (ULong64_t)(k.src + 7)) ^ (std::hash<Int_t>()(k.reco) << 1);
    }
};
inline void LoadPhotonTruthLinks(TFile *f, std::unordered_map<LinkKey, Int_t, LinkKeyHash> &out)
{
    TTree *t = GetTree(f, "recoTruthLinks");
    Int_t src = -1, reco_type = -1, truth_type = -1, link_class = -1, reco_index = -1, truth_index = -1;
    ULong64_t lo = 0;
    std::vector<std::string> pr;
    Bind(t, "source_file_index", &src, pr);
    Bind(t, "event_id_lo", &lo, pr);
    Bind(t, "reco_type", &reco_type, pr);
    Bind(t, "truth_type", &truth_type, pr);
    Bind(t, "link_class", &link_class, pr);
    Bind(t, "reco_index", &reco_index, pr);
    Bind(t, "truth_index", &truth_index, pr);
    AbortOnProblems(pr, "recoTruthLinks", f->GetName());
    const Long64_t n = t->GetEntries();
    out.clear();
    for (Long64_t i = 0; i < n; ++i)
    {
        if (t->GetEntry(i) <= 0)
            throw std::runtime_error(std::string("failed reading ") + t->GetName() + " entry " + std::to_string(i) + " in " + f->GetName());
        if (reco_type != 1 || truth_type != 1 || link_class != 0 || reco_index < 0 || truth_index < 0) continue;
        if (!out.emplace(LinkKey{src, lo, reco_index}, truth_index).second)
            throw std::runtime_error("duplicate photon truth link in " + std::string(f->GetName()));
    }
}

// Leading truth-jet pT per event (needed only for the sample_map window of
// jet samples in products without stitching weights, i.e. Au+Au). The
// truthJets tree is large (4e5-6e5 rows per part), only three branches are read.
inline void LoadLeadingTruthJetPt(TFile *f, std::unordered_map<EventKey, double, EventKeyHash> &out)
{
    TTree *t = GetTree(f, "truthJets");
    Int_t src = -1;
    ULong64_t lo = 0;
    Double_t pt = 0;
    std::vector<std::string> pr;
    Bind(t, "source_file_index", &src, pr);
    Bind(t, "event_id_lo", &lo, pr);
    Bind(t, "truth_jet_pt", &pt, pr);
    AbortOnProblems(pr, "truthJets", f->GetName());
    const Long64_t n = t->GetEntries();
    out.clear();
    for (Long64_t i = 0; i < n; ++i)
    {
        if (t->GetEntry(i) <= 0)
            throw std::runtime_error(std::string("failed reading ") + t->GetName() + " entry " + std::to_string(i) + " in " + f->GetName());
        auto &v = out[EventKey{src, lo}];
        if (pt > v) v = pt;
    }
}


using EventRows = std::unordered_map<EventKey, std::vector<int>, EventKeyHash>;
template <class Row>
inline EventRows GroupByEvent(const std::vector<Row> &rows)
{
    EventRows groups;
    for (size_t i = 0; i < rows.size(); ++i) groups[rows[i].key()].push_back(static_cast<int>(i));
    return groups;
}
inline void ValidateEventKeys(const std::vector<Event> &events, const EventRows &photons, const EventRows &truths)
{
    std::unordered_set<EventKey, EventKeyHash> keys;
    for (const auto &event : events)
        if (!keys.insert(event.key()).second) throw std::runtime_error("duplicate event key in a part");
    for (const auto &group : photons)
        if (!keys.count(group.first)) throw std::runtime_error("photon references an absent event");
    for (const auto &group : truths)
        if (!keys.count(group.first)) throw std::runtime_error("truth photon references an absent event");
}
inline std::vector<int> MatchPhotons(const std::vector<Photon> &photons, const std::vector<TruthPhoton> &truths,
                                    const EventRows &photon_rows, const EventRows &truth_rows,
                                    const std::unordered_map<LinkKey, Int_t, LinkKeyHash> &links, bool use_links)
{
    std::vector<int> matches(photons.size(), -1);
    for (const auto &group : photon_rows)
    {
        auto found = truth_rows.find(group.first);
        const std::vector<int> empty;
        const auto &truth = found == truth_rows.end() ? empty : found->second;
        for (size_t local = 0; local < group.second.size(); ++local)
        {
            const int photon = group.second[local];
            if (use_links)
            {
                auto link = links.find(LinkKey{group.first.src, group.first.lo, static_cast<int>(local)});
                if (link == links.end()) continue;
                if (link->second < 0 || link->second >= static_cast<int>(truth.size()))
                    throw std::runtime_error("photon link has an invalid truth index");
                matches[photon] = truth[link->second];
            }
            else if (photons[photon].truth_barcode >= 0)
            {
                for (int index : truth)
                    if (truths[index].barcode == photons[photon].truth_barcode)
                    {
                        if (matches[photon] >= 0) throw std::runtime_error("ambiguous truth barcode within an event");
                        matches[photon] = index;
                    }
            }
        }
    }
    if (use_links)
        for (const auto &link : links)
        {
            auto found = photon_rows.find(EventKey{link.first.src, link.first.lo});
            if (found == photon_rows.end() || link.first.reco >= static_cast<int>(found->second.size()))
                throw std::runtime_error("photon link references an absent reconstructed candidate");
        }
    return matches;
}

} // namespace PJ
