// Sample stitching and optional vertex weighting.
#pragma once
#include "../../support/CrossSectionWeights.h"
#include "PhotonJetSelection.h"
#include <TFile.h>
#include <TH1.h>
#include <yaml-cpp/yaml.h>
#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace PJ
{
struct SampleRange { int first, last; std::string sample; };
inline std::vector<SampleRange> LoadSampleMap(const std::string &path)
{
    std::vector<SampleRange> v;
    if (path.empty()) return v;
    YAML::Node n = YAML::LoadFile(path);
    const YAML::Node r = n["ranges"] ? n["ranges"] : n;
    for (const auto &x : r)
    {
        v.push_back(SampleRange{x["first"].as<int>(), x["last"].as<int>(), x["sample"].as<std::string>()});
    }
    std::sort(v.begin(), v.end(), [](const SampleRange &a, const SampleRange &b) { return a.first < b.first; });
    for (size_t i = 0; i < v.size(); ++i)
        if (v[i].first < 0 || v[i].last < v[i].first ||
            (i && v[i].first <= v[i - 1].last) || !PPG12::GetSampleConfig(v[i].sample).valid)
            throw std::runtime_error("invalid or overlapping sample-map range in " + path);
    return v;
}
inline std::string SampleForPart(const std::vector<SampleRange> &m, int part)
{
    for (const auto &r : m) if (part >= r.first && part <= r.last) return r.sample;
    return "";
}

// Optional reco-vertex reweight histogram (data/MC ratio vs z), used by
// weight_mode recompute/sample_map. Returns 1 when no file is configured.
struct VertexWeight
{
    std::unique_ptr<TH1> h;
    double operator()(double z) const
    {
        if (!h) return 1.0;
        int b = h->FindBin(z);
        b = std::max(1, std::min(h->GetNbinsX(), b));
        const double weight = h->GetBinContent(b);
        if (!std::isfinite(weight) || weight <= 0) throw std::runtime_error("invalid vertex weight");
        return weight;
    }
};
inline VertexWeight LoadVertexWeight(const std::string &path, const std::string &hname = "h_vertex_reweight")
{
    VertexWeight vw;
    if (path.empty()) return vw;
    std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
    if (!f || f->IsZombie()) throw std::runtime_error("PhotonJetReader: cannot open vertex_weight_file " + path);
    TH1 *h = dynamic_cast<TH1 *>(f->Get(hname.c_str()));
    if (!h || h->GetDimension() != 1) throw std::runtime_error("missing or invalid vertex histogram " + hname + " in " + path);
    vw.h.reset(dynamic_cast<TH1 *>(h->Clone("pj_vertex_weight")));
    vw.h->SetDirectory(nullptr);
    f->Close();
    return vw;
}


// Return false for events outside the configured analysis/sample window.
inline bool SimulationWeight(const Event &event, const Cuts &cuts, const PPG12::SampleConfig &sample,
                              const VertexWeight &vertex, const std::vector<TruthPhoton> &truth,
                              const EventRows &truth_rows,
                              const std::unordered_map<EventKey, double, EventKeyHash> &leading_jet,
                              double &weight)
{
    if (!std::isfinite(event.vertex_z) || !std::isfinite(event.event_weight))
        throw std::runtime_error("nonfinite simulation vertex or event weight");
    if (std::fabs(event.vertex_z) > cuts.vertex_cut) return false;
    if (cuts.system == "auau" && !InCentralityClass(event.centrality, cuts)) return false;
    if (cuts.weight_mode == "stored" || cuts.weight_mode == "recompute")
    {
        if (event.event_weight == 1.0) return false; // Producer's out-of-window/no-vertex sentinel.
        weight = cuts.weight_mode == "stored" ? event.event_weight : sample.weight * vertex(event.vertex_z);
    }
    else
    {
        double leading = -1;
        if (sample.isbackground)
        {
            auto found = leading_jet.find(event.key());
            if (found != leading_jet.end()) leading = found->second;
            if (!(leading >= sample.jet_pt_lower && leading < sample.jet_pt_upper)) return false;
        }
        else
        {
            auto found = truth_rows.find(event.key());
            if (found != truth_rows.end())
                for (int index : found->second)
                    if (truth[index].prompt_class == 1 || truth[index].prompt_class == 2)
                        leading = std::max(leading, truth[index].pt);
            if (!(leading >= sample.photon_pt_lower && leading < sample.photon_pt_upper)) return false;
        }
        if (event.event_weight == 1.0 && event.vertex_z == 0.0) return false;
        weight = sample.weight * vertex(event.vertex_z);
    }
    if (!std::isfinite(weight) || weight <= 0) throw std::runtime_error("simulation weight must be finite and positive");
    return true;
}

} // namespace PJ
