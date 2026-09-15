// Checked PhotonJetTrees_v1 input, held one original part at a time.
// Records contain input data only; cuts and weights belong to PhotonJetPhysics.h.
#pragma once
#include <TFile.h>
#include <TLeaf.h>
#include <TTree.h>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace PJ {
using EventKey = std::tuple<Int_t, ULong64_t, ULong64_t>;
struct Identity {
    Int_t source_file_index = -1;
    ULong64_t event_id_hi = 0, event_id_lo = 0;
    EventKey key() const { return {source_file_index, event_id_hi, event_id_lo}; }
};

// These fixed lists keep the raw scalar declarations and checked bindings in
// agreement. They are local to this header and do not define persisted columns.
#define PJ_EVENT_FIELDS(X) \
    X(Long64_t, source_entry) X(Int_t, run) \
    X(Long64_t, event_sequence) X(Long64_t, physical_event_sequence) \
    X(ULong64_t, trigger_bits) X(ULong64_t, live_trigger_bits) X(ULong64_t, scaled_trigger_bits) \
    X(Int_t, scaled_bit30) X(Double_t, vertex_z) X(Double_t, centrality) \
    X(Double_t, event_weight) X(Double_t, total_calo_energy) X(Int_t, terminal_status)
#define PJ_MB_FIELDS(X) \
    X(Int_t, scaled_bit22) X(Int_t, minimum_bias_pass) X(Int_t, nominal_event_selection_pass)
#define PJ_PHOTON_FIELDS(X) \
    X(ULong64_t, candidate_id_hi) X(ULong64_t, candidate_id_lo) \
    X(Double_t, photon_et) X(Double_t, photon_eta) \
    X(Double_t, bdt_score) X(Double_t, bdt_tight_threshold) \
    X(Double_t, bdt_nontight_low_threshold) X(Double_t, bdt_nontight_high_threshold) \
    X(Int_t, bdt_is_tight) X(Int_t, bdt_is_nontight) \
    X(Double_t, iso_r03) X(Double_t, iso_r03_threshold) X(Double_t, iso_r03_nonisolated_threshold) X(Int_t, iso_r03_pass) \
    X(Double_t, iso_r04) X(Double_t, iso_r04_threshold) X(Double_t, iso_r04_nonisolated_threshold) X(Int_t, iso_r04_pass) \
    X(Int_t, truth_barcode)
#define PJ_TRUTH_FIELDS(X) \
    X(ULong64_t, truth_photon_id_hi) X(ULong64_t, truth_photon_id_lo) \
    X(Double_t, truth_photon_pt) X(Double_t, truth_photon_eta) \
    X(Int_t, prompt_class) X(Int_t, generator_barcode) X(Double_t, truth_isolation)
#define PJ_DECLARE(type, name) type name{};
struct Event : Identity {
    PJ_EVENT_FIELDS(PJ_DECLARE)
    Int_t scaled_bit22 = -1, minimum_bias_pass = -1, nominal_event_selection_pass = -1;
};
struct Photon : Identity {
    PJ_PHOTON_FIELDS(PJ_DECLARE)
    // Resolved producer shower-shape inputs, with no selection applied.
    Double_t weta = 0, wphi = 0, e11e33 = 0, et1 = 0, et2 = 0, et3 = 0, et4 = 0, e32e35 = 0;
    // The RNG uses this original position within the event, never the producer
    // encounter ordinal or a position after selection.
    Int_t original_index = -1;
};
struct TruthPhoton : Identity {
    PJ_TRUTH_FIELDS(PJ_DECLARE)
    Int_t original_index = -1;
};
#undef PJ_DECLARE

struct PartData {
    std::vector<Event> events;
    std::vector<Photon> photons;
    std::vector<TruthPhoton> truths;
    std::map<EventKey, size_t> event_index;
    // Entries refer to the original part's photon/truth row order. links[e][p]
    // is a truth index within event e, or -1, preserving zero-reco events.
    std::vector<std::vector<Long64_t>> reco_rows, truth_rows;
    std::vector<std::vector<Int_t>> links;
    std::vector<Double_t> leading_jet;
    std::vector<Long64_t> jet_count;
    bool links_available = false, jets_available = false;
    Long64_t link_rows = 0, photon_links = 0, jet_rows = 0;
};

namespace ReaderDetail {
inline TTree *Tree(TFile &file, const char *name, bool required = true)
{
    auto *object = file.Get(name);
    auto *tree = dynamic_cast<TTree *>(object);
    if ((!tree && required) || (object && !tree))
        throw std::runtime_error(std::string("missing or wrong-type tree ") + name + " in " + file.GetName());
    if (tree) { tree->SetBranchStatus("*", 0); tree->SetImplicitMT(false); }
    return tree;
}
// Also release stack-bound addresses when an invalid row throws during loading.
struct ResetBranches {
    TTree *tree;
    ~ResetBranches() { if (tree) tree->ResetBranchAddresses(); }
};
template<class T>
inline void Bind(TTree &tree, const char *name, T &value, const char *type, bool required = true)
{
    auto *leaf = tree.GetLeaf(name);
    if (!leaf && !required && !tree.GetBranch(name)) return;
    if (!leaf || std::string(leaf->GetTypeName()) != type || leaf->GetLeafCount() || leaf->GetLenStatic() != 1 ||
        leaf->GetBranch()->GetListOfLeaves()->GetEntries() != 1)
        throw std::runtime_error(std::string("missing or wrong-type scalar ") + tree.GetName() + "." + name);
    tree.SetBranchStatus(name, 1);
    if (tree.SetBranchAddress(name, &value) != TTree::kMatch)
        throw std::runtime_error(std::string("cannot bind ") + tree.GetName() + "." + name);
}
inline void Read(TTree &tree, Long64_t entry)
{
    if (tree.GetEntry(entry) <= 0)
        throw std::runtime_error(std::string("read failed: ") + tree.GetName() + " entry " + std::to_string(entry));
}
inline void BindIdentity(TTree &tree, Identity &row)
{
    Bind(tree, "source_file_index", row.source_file_index, "Int_t");
    Bind(tree, "event_id_hi", row.event_id_hi, "ULong64_t");
    Bind(tree, "event_id_lo", row.event_id_lo, "ULong64_t");
}
inline void BindEvent(TTree &tree, Event &row, bool require_mb)
{
    BindIdentity(tree, row);
#define PJ_BIND(type, name) Bind(tree, #name, row.name, #type);
    PJ_EVENT_FIELDS(PJ_BIND)
#undef PJ_BIND
#define PJ_BIND(type, name) Bind(tree, #name, row.name, #type, require_mb);
    PJ_MB_FIELDS(PJ_BIND)
#undef PJ_BIND
}
inline bool Equal(Double_t a, Double_t b) { return a == b || (std::isnan(a) && std::isnan(b)); }
template<class T> inline bool Equal(T a, T b) { return a == b; }
inline void CheckEventCopy(const Event &event, const Event &copy)
{
#define PJ_COMPARE(type, name) if (!Equal(event.name, copy.name)) \
    throw std::runtime_error("photon/event disagreement: " #name);
    PJ_EVENT_FIELDS(PJ_COMPARE)
    PJ_MB_FIELDS(PJ_COMPARE)
#undef PJ_COMPARE
}
inline size_t EventIndex(const PartData &data, const Identity &row, const char *what)
{
    const auto found = data.event_index.find(row.key());
    if (found == data.event_index.end()) throw std::runtime_error(std::string(what) + " references an absent event");
    return found->second;
}
inline Int_t ObjectIndex(size_t size)
{
    if (size > static_cast<size_t>(std::numeric_limits<Int_t>::max()))
        throw std::runtime_error("too many objects within one event");
    return static_cast<Int_t>(size);
}
} // namespace ReaderDetail

inline PartData LoadPart(TFile &file, bool require_mb, bool simulation, bool auau = false)
{
    using namespace ReaderDetail;
    if (file.IsZombie() || file.TestBit(TFile::kRecovered))
        throw std::runtime_error(std::string("corrupt or recovered ROOT file: ") + file.GetName());
    PartData data;
    auto *events = Tree(file, "events");
    ResetBranches reset_events{events};
    Event event;
    BindEvent(*events, event, require_mb);
    data.events.reserve(events->GetEntries());
    for (Long64_t i = 0; i < events->GetEntries(); ++i) {
        Read(*events, i);
        if (event.source_file_index < 0 || event.source_entry < 0 || !data.event_index.emplace(event.key(), data.events.size()).second)
            throw std::runtime_error("invalid or duplicate event identity");
        data.events.push_back(event);
    }
    events->ResetBranchAddresses();
    const auto n = data.events.size();
    data.reco_rows.resize(n); data.truth_rows.resize(n); data.links.resize(n);
    data.leading_jet.assign(n, 0); data.jet_count.assign(n, 0);

    auto *photons = Tree(file, "photons");
    ResetBranches reset_photons{photons};
    Event event_copy;
    Photon photon;
    BindEvent(*photons, event_copy, require_mb);
#define PJ_BIND(type, name) Bind(*photons, #name, photon.name, #type);
    PJ_PHOTON_FIELDS(PJ_BIND)
#undef PJ_BIND
    // Validate unused raw diagnostics without retaining one copy per photon.
    Double_t photon_phi = 0;
    Int_t photon_encounter_ordinal = 0, bdt_is_not_tight = 0, bdt_input_count = 0, truth_matched = 0;
    Bind(*photons, "photon_phi", photon_phi, "Double_t");
    Bind(*photons, "photon_encounter_ordinal", photon_encounter_ordinal, "Int_t");
    Bind(*photons, "bdt_is_not_tight", bdt_is_not_tight, "Int_t");
    Bind(*photons, "bdt_input_count", bdt_input_count, "Int_t");
    Bind(*photons, "truth_matched", truth_matched, "Int_t");
    // Only the resolved eight shape values survive in each stored Photon.
    Double_t native[7]{}, inputs[14]{};
    const char *native_names[] = {"native_weta_cogx", "native_wphi_cogx", "native_weta33_cogx", "native_wphi33_cogx",
                                 "native_e11_over_e33", "native_e32_over_e35", "native_et1"};
    for (int j = 0; j < 7; ++j) Bind(*photons, native_names[j], native[j], "Double_t");
    for (int j = 0; j < 14; ++j) {
        const std::string name = "bdt_input_" + std::string(j < 10 ? "0" : "") + std::to_string(j);
        Bind(*photons, name.c_str(), inputs[j], "Double_t");
    }
    data.photons.reserve(photons->GetEntries());
    std::set<std::tuple<EventKey, ULong64_t, ULong64_t>> candidate_ids;
    for (Long64_t i = 0; i < photons->GetEntries(); ++i) {
        Read(*photons, i);
        const auto e = EventIndex(data, event_copy, "photon");
        CheckEventCopy(data.events[e], event_copy);
        static_cast<Identity &>(photon) = event_copy;
        if (!candidate_ids.emplace(photon.key(), photon.candidate_id_hi, photon.candidate_id_lo).second)
            throw std::runtime_error("duplicate photon identity within an event");
        photon.original_index = ObjectIndex(data.reco_rows[e].size());
        photon.weta = inputs[1]; photon.wphi = inputs[2];
        const int offset = auau ? 7 : 5;
        photon.e11e33 = inputs[offset]; photon.et1 = inputs[offset + 1]; photon.et2 = inputs[offset + 2];
        photon.et3 = inputs[offset + 3]; photon.et4 = inputs[offset + 4]; photon.e32e35 = inputs[offset + 5];
        data.reco_rows[e].push_back(i);
        data.links[e].push_back(-1);
        data.photons.push_back(photon);
    }
    photons->ResetBranchAddresses();

    auto *truths = Tree(file, "truthPhotons", simulation);
    ResetBranches reset_truths{truths};
    if (truths) {
        TruthPhoton truth;
        BindIdentity(*truths, truth);
#define PJ_BIND(type, name) Bind(*truths, #name, truth.name, #type);
        PJ_TRUTH_FIELDS(PJ_BIND)
#undef PJ_BIND
        Double_t truth_photon_phi = 0;
        Int_t source_role = 0;
        Bind(*truths, "truth_photon_phi", truth_photon_phi, "Double_t");
        Bind(*truths, "source_role", source_role, "Int_t");
        data.truths.reserve(truths->GetEntries());
        std::set<std::tuple<EventKey, ULong64_t, ULong64_t>> truth_ids;
        for (Long64_t i = 0; i < truths->GetEntries(); ++i) {
            Read(*truths, i);
            const auto e = EventIndex(data, truth, "truth photon");
            if (!truth_ids.emplace(truth.key(), truth.truth_photon_id_hi, truth.truth_photon_id_lo).second)
                throw std::runtime_error("duplicate truth identity within an event");
            truth.original_index = ObjectIndex(data.truth_rows[e].size());
            data.truth_rows[e].push_back(i);
            data.truths.push_back(truth);
        }
        truths->ResetBranchAddresses();
    }

    // Stream side trees without retaining their rows. Check every event key,
    // even for link classes unused by this histogram analysis.
    auto *links = Tree(file, "recoTruthLinks", false);
    ResetBranches reset_links{links};
    data.links_available = links != nullptr;
    if (links) {
        Identity key;
        BindIdentity(*links, key);
        Int_t reco_type = -1, truth_type = -1, link_class = -1, reco_index = -1, truth_index = -1;
        ULong64_t reco_id_hi = 0, reco_id_lo = 0, truth_id_hi = 0, truth_id_lo = 0;
#define PJ_BIND(type, name) Bind(*links, #name, name, #type);
        PJ_BIND(Int_t, reco_type) PJ_BIND(Int_t, truth_type) PJ_BIND(Int_t, link_class)
        PJ_BIND(Int_t, reco_index) PJ_BIND(Int_t, truth_index)
        PJ_BIND(ULong64_t, reco_id_hi) PJ_BIND(ULong64_t, reco_id_lo)
        PJ_BIND(ULong64_t, truth_id_hi) PJ_BIND(ULong64_t, truth_id_lo)
#undef PJ_BIND
        data.link_rows = links->GetEntries();
        for (Long64_t i = 0; i < data.link_rows; ++i) {
            Read(*links, i);
            const auto e = EventIndex(data, key, "link");
            if (reco_type != 1 || truth_type != 1 || link_class != 0) continue;
            if (reco_index < 0 || truth_index < 0 || reco_index >= static_cast<Long64_t>(data.reco_rows[e].size()) ||
                truth_index >= static_cast<Long64_t>(data.truth_rows[e].size())) throw std::runtime_error("invalid photon link index");
            const auto &reco = data.photons[data.reco_rows[e][reco_index]];
            const auto &truth = data.truths[data.truth_rows[e][truth_index]];
            if (reco.candidate_id_hi != reco_id_hi || reco.candidate_id_lo != reco_id_lo ||
                truth.truth_photon_id_hi != truth_id_hi || truth.truth_photon_id_lo != truth_id_lo)
                throw std::runtime_error("photon link index/identity disagreement");
            if (data.links[e][reco_index] != -1) throw std::runtime_error("duplicate photon link");
            data.links[e][reco_index] = truth_index;
            ++data.photon_links;
        }
        links->ResetBranchAddresses();
    }
    auto *jets = Tree(file, "truthJets", false);
    ResetBranches reset_jets{jets};
    data.jets_available = jets != nullptr;
    if (jets) {
        Identity key;
        BindIdentity(*jets, key);
        Double_t pt = 0;
        Bind(*jets, "truth_jet_pt", pt, "Double_t");
        data.jet_rows = jets->GetEntries();
        for (Long64_t i = 0; i < data.jet_rows; ++i) {
            Read(*jets, i);
            const auto e = EventIndex(data, key, "truth jet");
            if (!std::isfinite(pt)) throw std::runtime_error("nonfinite truth jet pT");
            if (!data.jet_count[e] || pt > data.leading_jet[e]) data.leading_jet[e] = pt;
            ++data.jet_count[e];
        }
        jets->ResetBranchAddresses();
    }
    if (!simulation && (!data.truths.empty() || data.jet_rows || data.link_rows))
        throw std::runtime_error("data product contains simulation objects");
    return data;
}
#undef PJ_EVENT_FIELDS
#undef PJ_MB_FIELDS
#undef PJ_PHOTON_FIELDS
#undef PJ_TRUTH_FIELDS
} // namespace PJ
