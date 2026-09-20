// Reads the flat trees of one part, each once and front to back, and rebuilds the events:
//   events -> truthPhotons -> photons -> jets -> truthJets -> recoTruthLinks.
// Every row carries the key of its event, and all pairing happens inside one event: the link
// tree gives, for a photon, the position of its truth photon within the same event.
// The flow is the same for pp and Au+Au. No selection is applied here.
#pragma once
#include "PhotonJetObjects.h"
#include "PhotonJetConfig.h"
#include "CrossSectionWeights.h"
#include <TFile.h>
#include <TObjString.h>
#include <TLeaf.h>
#include <TTree.h>
#include <cmath>
#include <cstring>
#include <map>

namespace PJ {
namespace ReaderDetail {
inline TTree &Tree(TFile &file, const char *name)
{
    TTree *tree = dynamic_cast<TTree *>(file.Get(name));
    if (!tree) throw std::runtime_error(std::string("missing tree ") + name + " in " + file.GetName());
    tree->SetBranchStatus("*", 0);   // only the bound branches are read
    return *tree;
}
// A missing branch or one of another type is refused by ROOT, which makes a schema change fail loudly.
template<class T> inline void Bind(TTree &tree, const char *name, T &value)
{
    tree.SetBranchStatus(name, 1);
    if (tree.SetBranchAddress(name, &value) != TTree::kMatch)
        throw std::runtime_error(std::string("missing or wrong-type branch ") + tree.GetName() + "." + name);
}
inline void Read(TTree &tree, Long64_t entry)
{
    if (tree.GetEntry(entry) <= 0) throw std::runtime_error(std::string("read failed: ") + tree.GetName());
}
// A C-string branch is read into `buffer`. The producer declared these leaves 5 characters long and
// ROOT truncates a string to the declared length, so the length is raised to the buffer size first.
inline void BindText(TTree &tree, const char *name, std::vector<char> &buffer)
{
    TLeaf *leaf = tree.GetLeaf(name);
    if (!leaf || std::string(leaf->GetTypeName()) != "Char_t")
        throw std::runtime_error(std::string("missing or wrong-type text ") + tree.GetName() + "." + name);
    buffer.assign(64, '\0');
    leaf->SetLen(buffer.size());
    tree.SetBranchStatus(name, 1);
    if (tree.SetBranchAddress(name, buffer.data()) < 0) throw std::runtime_error(std::string("cannot bind ") + name);
}
inline void BindIdentity(TTree &tree, Identity &row)
{
    Bind(tree, "source_file_index", row.source_file_index);
    Bind(tree, "event_id_hi", row.event_id_hi); Bind(tree, "event_id_lo", row.event_id_lo);
}
// Bind the branch to the member of `row` that has the same name.
#define PJ_BIND(name) Bind(tree, #name, row.name)
inline void BindEvent(TTree &tree, Event &row, bool simulation)
{
    BindIdentity(tree, row);
    PJ_BIND(run); PJ_BIND(scaled_trigger_bits); PJ_BIND(scaled_bit30);
    PJ_BIND(vertex_z); PJ_BIND(centrality);
    PJ_BIND(total_calo_energy); PJ_BIND(emcal_total_energy);
    if (simulation) { PJ_BIND(truth_vertex_z); PJ_BIND(truth_mb_vertex_z); }
}
// One recoTruthLinks row. Class 0 with both types 1 links a photon to a truth photon.
struct Link : Identity {
    Int_t reco_type = 0, truth_type = 0, link_class = 0, reco_index = -1, truth_index = -1;
};
// The generator sample of one source file, from its name in the part's source_files metadata:
// run28_photonjet10 -> photon10, run28_jet12_double -> jet12_double, as CrossSectionWeights.h knows them.
inline PPG12::SampleConfig SampleOf(std::string name)
{
    if (name.rfind("run28_", 0) == 0) name = name.substr(6);
    if (name.rfind("photonjet", 0) == 0) name = "photon" + name.substr(9);
    return PPG12::GetSampleConfig(name);
}
// Tag every event with the weight and truth-pT window of its sample.
inline void AssignSamples(TFile &file, std::vector<Event> &events)
{
    TObjString *text = dynamic_cast<TObjString *>(file.Get("source_files"));
    if (!text) throw std::runtime_error(std::string("missing source_files metadata in ") + file.GetName());
    struct Source { PPG12::SampleConfig sample; bool double_interaction = false; };
    std::map<Int_t, Source> sources;
    for (const YAML::Node &node : YAML::Load(text->GetString().Data())) {
        const YAML::Node occurrence = node["source_occurrences"][0];
        Source source;
        source.sample = SampleOf(occurrence["sample"].as<std::string>());
        source.double_interaction = occurrence["si_di_role"].as<std::string>() == "DI";
        if (!source.sample.valid)
            throw std::runtime_error("unknown generator sample " + occurrence["sample"].as<std::string>() + " in " + file.GetName());
        sources[node["source_file_index"].as<Int_t>()] = source;
    }
    for (Event &event : events) {
        const Source &source = sources.at(event.source_file_index);
        event.sample_weight = source.sample.weight;
        event.jet_sample = source.sample.isbackground;
        event.window_low = source.sample.isbackground ? source.sample.jet_pt_lower : source.sample.photon_pt_lower;
        event.window_high = source.sample.isbackground ? source.sample.jet_pt_upper : source.sample.photon_pt_upper;
        event.cluster_et_upper = source.sample.isbackground ? source.sample.cluster_ET_upper : 1e9;
        event.double_interaction = source.double_interaction;
    }
}
// Position of the event that `row` belongs to. `event_index` maps an event key to that position.
inline size_t EventIndex(const std::map<EventKey, size_t> &event_index, const Identity &row)
{
    const std::map<EventKey, size_t>::const_iterator found = event_index.find(row.key());
    if (found == event_index.end()) throw std::runtime_error("object references an absent event");
    return found->second;
}
} // namespace ReaderDetail

// The trees are read here and nowhere else, and the file is closed afterwards.
inline std::vector<Event> LoadPart(TFile &file, const Config &config, bool simulation)
{
    using namespace ReaderDetail;
    std::vector<Event> events;
    // Event key -> position in `events`. The keys are unique within one file only, so this is
    // rebuilt for every part. Rows of one event need not be contiguous in the other trees.
    std::map<EventKey, size_t> event_index;

    /////////////////////////////////////
    // events: one row per event
    /////////////////////////////////////
    {
        TTree &tree = Tree(file, "events");
        Event row;
        BindEvent(tree, row, simulation);
        events.reserve(tree.GetEntries());
        for (Long64_t i = 0; i < tree.GetEntries(); ++i) {
            Read(tree, i);
            if (!event_index.emplace(row.key(), events.size()).second) throw std::runtime_error("duplicate event identity");
            events.push_back(row);
        }
    }
    if (simulation && config.sample_weights) AssignSamples(file, events);

    /////////////////////////////////////
    // truthPhotons, simulation only: appended to their event
    /////////////////////////////////////
    if (simulation) {
        TTree &tree = Tree(file, "truthPhotons");
        TruthPhoton row; BindIdentity(tree, row);
        PJ_BIND(truth_photon_pt); PJ_BIND(truth_photon_eta); PJ_BIND(truth_photon_phi);
        PJ_BIND(prompt_class);
        // Truth isolation uses the same cone as the reconstructed photon.
        Bind(tree, config.cone_size == 3 ? "truth_isolation_r03" : "truth_isolation_r04", row.truth_isolation);
        for (Long64_t i = 0; i < tree.GetEntries(); ++i) {
            Read(tree, i);
            events[EventIndex(event_index, row)].truths.push_back(row);
        }
    }

    /////////////////////////////////////
    // photons: appended to their event
    /////////////////////////////////////
    {
        TTree &tree = Tree(file, "photons");
        Photon row; BindIdentity(tree, row);
        PJ_BIND(photon_et); PJ_BIND(photon_eta); PJ_BIND(photon_phi);
        PJ_BIND(bdt_score); PJ_BIND(bdt_tight_threshold);
        PJ_BIND(bdt_nontight_low_threshold); PJ_BIND(bdt_nontight_high_threshold);
        PJ_BIND(iso_r03); PJ_BIND(iso_r03_threshold); PJ_BIND(iso_r03_nonisolated_threshold);
        PJ_BIND(iso_r04); PJ_BIND(iso_r04_threshold); PJ_BIND(iso_r04_nonisolated_threshold);
        Bind(tree, "shower_weta_cogx", row.weta); Bind(tree, "shower_wphi_cogx", row.wphi);
        Bind(tree, "shower_e11_over_e33", row.e11e33); Bind(tree, "shower_e32_over_e35", row.e32e35);
        Bind(tree, "shower_et1", row.et1); Bind(tree, "shower_et2", row.et2);
        Bind(tree, "shower_et3", row.et3); Bind(tree, "shower_et4", row.et4);
        PJ_BIND(npb_score); PJ_BIND(shower_center_eta_index); PJ_BIND(shower_center_phi_index);
        for (Long64_t i = 0; i < tree.GetEntries(); ++i) {
            Read(tree, i);
            events[EventIndex(event_index, row)].photons.push_back(row);
        }
    }

    /////////////////////////////////////
    // jets and truthJets. The sample windows of jet samples need the truth jets.
    /////////////////////////////////////
    {
        TTree &tree = Tree(file, "jets");
        Identity key; BindIdentity(tree, key);
        Jet row;
        std::vector<char> input;
        Bind(tree, "jet_pt", row.pt); Bind(tree, "jet_eta", row.eta); Bind(tree, "jet_phi", row.phi);
        Bind(tree, "jet_radius", row.radius);
        BindText(tree, "jet_input_identity", input);
        Long64_t usable = 0;
        for (Long64_t i = 0; i < tree.GetEntries(); ++i) {
            Read(tree, i);
            row.configured_input = std::strcmp(input.data(), config.jet_input.c_str()) == 0;
            if (row.configured_input && std::fabs(row.radius - config.jet_radius) < 1e-6) ++usable;
            events[EventIndex(event_index, key)].jets.push_back(row);
        }
        // A wrong radius or input name would otherwise give jet histograms that are silently empty.
        if (tree.GetEntries() > 0 && usable == 0)
            throw std::runtime_error("no jets with radius " + std::to_string(config.jet_radius) + " and input '" + config.jet_input + "' in " + file.GetName());
    }
    if (simulation) {
        TTree &tree = Tree(file, "truthJets");
        Identity key; BindIdentity(tree, key);
        Jet row;
        Bind(tree, "truth_jet_pt", row.pt); Bind(tree, "truth_jet_eta", row.eta); Bind(tree, "truth_jet_phi", row.phi);
        Bind(tree, "truth_jet_radius", row.radius);
        for (Long64_t i = 0; i < tree.GetEntries(); ++i) {
            Read(tree, i);
            events[EventIndex(event_index, key)].truth_jets.push_back(row);
        }
    }

    /////////////////////////////////////
    // recoTruthLinks, simulation only: the producer's photon-to-truth-photon match, Photon::truth_index
    /////////////////////////////////////
    if (simulation) {
        TTree &tree = Tree(file, "recoTruthLinks");
        Link row; BindIdentity(tree, row);
        PJ_BIND(reco_type); PJ_BIND(truth_type); PJ_BIND(link_class);
        PJ_BIND(reco_index); PJ_BIND(truth_index);
        for (Long64_t i = 0; i < tree.GetEntries(); ++i) {
            Read(tree, i);
            Event &event = events[EventIndex(event_index, row)];
            if (row.link_class != 0 || row.reco_type != 1 || row.truth_type != 1) continue;
            // reco_index and truth_index count the photons and truth photons of the event in row order.
            if (row.reco_index < 0 || row.truth_index < 0 ||
                row.reco_index >= static_cast<Long64_t>(event.photons.size()) ||
                row.truth_index >= static_cast<Long64_t>(event.truths.size()))
                throw std::runtime_error("invalid photon link index");
            Photon &photon = event.photons[row.reco_index];
            if (photon.truth_index != -1) throw std::runtime_error("a photon with two truth links");
            photon.truth_index = row.truth_index;
        }
    }
    return events;
}
#undef PJ_BIND
} // namespace PJ
