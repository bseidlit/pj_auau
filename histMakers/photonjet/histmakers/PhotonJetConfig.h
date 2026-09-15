// Photonjet configuration: existing YAML keys, checked once before processing.
#pragma once
#include <yaml-cpp/yaml.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <initializer_list>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace PJ
{
enum class Product { Data, Signal, Inclusive };
inline Product ParseProduct(const std::string &name)
{
    if (name == "data") return Product::Data;
    if (name == "sim_signal") return Product::Signal;
    if (name == "sim_inclusive") return Product::Inclusive;
    throw std::runtime_error("unknown product '" + name + "' (expected data, sim_signal, sim_inclusive)");
}
inline const char *ProductName(Product p)
{
    return p == Product::Data ? "data" : p == Product::Signal ? "sim_signal" : "sim_inclusive";
}
inline const char *InputListKey(Product p)
{
    return p == Product::Data ? "data_files_list" : p == Product::Signal ? "sim_signal_files_list" : "sim_inclusive_files_list";
}

struct Cuts
{
    // photonjet block
    std::string system = "pp";           // pp | auau
    std::string threshold_source = "formula"; // formula | tree
    std::string weight_mode = "stored";  // stored | recompute | sample_map
    std::string match_source = "barcode"; // barcode (photons.truth_barcode) | links (recoTruthLinks class 1,1,0)
    std::string sample_map_file;
    std::string vertex_weight_file;      // optional TH1 "h_vertex_reweight" for recompute mode
    int apply_shower_shape_windows = 1;
    double min_photon_et = 0;
    std::vector<double> centrality_bins; // AuAu centrality edges; each class is [lo, hi)
    int require_prompt_class_2 = 1;      // treat fragmentation (class 2) as signal like PPG12

    // event level
    double vertex_cut = 30.0;
    int run_min = -1, run_max = -1;
    std::string run_list_file;
    std::set<int> runs;
    std::string flag_check = "strict"; // strict | report | off (data formula mode)
    unsigned int random_seed = 42;
    std::vector<int> trigger_used;       // pp bit list (30)
    int apply_trigger_eff_correction = 0;
    double trigger_eff_p0 = 1.0, trigger_eff_mu = -2.04, trigger_eff_beta = 3.33;

    // kinematics
    double reco_min_ET = 5.0;
    std::vector<double> eta_bins, pT_bins, pT_bins_truth;

    // isolation (PPG12 nominal: topo R=0.4, use_topo_iso 2)
    int use_topo_iso = 2;                // 1 -> iso_r03, 2 -> iso_r04
    double reco_iso_min = -20, reco_iso_max_b = 0.49, reco_iso_max_s = 0.037;
    double reco_noniso_min_shift = 0.8, reco_noniso_max = 20;
    double mc_iso_scale = 1.0, mc_iso_shift = 0.0;
    double truth_iso_max = 4.0;

    // Selection uses scale/NL only; extra resolution is applied to the response.
    double cluster_escale = 1.0, cluster_escale_nl_slope = 0.0, cluster_escale_nl_ref_ET = 10.0;
    double eres_data_p0 = 0, eres_data_p1 = 0, eres_data_p2 = 0, eres_mc_p0 = 0, eres_mc_p1 = 0, eres_mc_p2 = 0;

    // common
    double common_e11_over_e33_min = 0, common_e11_over_e33_max = 0.98, common_wr_cogx_bound = 0, common_weta_cogx_bound = 2.0;

    // tight
    double t_weta_min = 0, t_weta_max_b = 1.0, t_weta_max_s = 0;
    double t_wphi_min = 0, t_wphi_max_b = 1.0, t_wphi_max_s = 0;
    double t_et1_min_b = 0.5, t_et1_min_s = 0, t_et1_max = 1;
    double t_et2_min = 0, t_et2_max = 1, t_et3_min = 0, t_et3_max = 1, t_et4_min = 0, t_et4_max = 1;
    double t_e11e33_min = 0, t_e11e33_max = 1, t_e32e35_min = 0.8, t_e32e35_max = 1;
    double t_bdt_max = 1.0, t_bdt_min_intercept = 0.0, t_bdt_min_slope = 0.0;

    // non-tight
    double nt_weta_min = 0, nt_weta_max = 1.0, nt_wphi_min = 0, nt_wphi_max = 1.0;
    double nt_et1_min = 0.6, nt_et1_max = 1, nt_et4_min = 0, nt_et4_max = 1;
    double nt_e11e33_min = 0, nt_e11e33_max = 1, nt_e32e35_min = 0.8, nt_e32e35_max = 1;
    double nt_bdt_min_intercept = 0.0, nt_bdt_min_slope = 0.0;
    double nt_bdt_max_intercept = 1.0, nt_bdt_max_slope = 0.0;
    int n_nt_fail = 1;
    int weta_on = 1, wphi_on = 1, et1_on = 1, et2_on = 1, et3_on = 1, e11_to_e33_on = 1, e32_to_e35_on = 1, et4_on = 1, bdt_on = 1;
    int weta_fail = 0, wphi_fail = 0, et1_fail = 0, e11_to_e33_fail = 0, e32_to_e35_fail = 0, bdt_fail = 0;

    // response prior reweight
    int unfold_reweight = 0;
    std::string trw_formula = "1";
    std::vector<double> trw_params;
    double trw_xmin = 0, trw_xmax = 100, trw_clamp_min = 0, trw_clamp_max = 1e9;
};

// Defaults apply only to absent keys. A malformed supplied value is an error.
template <class T>
inline T Y(const YAML::Node &n, const char *key, T def)
{
    if (!n || !n[key]) return def;
    try
    {
        T value = n[key].as<T>();
        if constexpr (std::is_floating_point<T>::value)
            if (!std::isfinite(value)) throw std::runtime_error("must be finite");
        return value;
    }
    catch (const std::exception &e)
    {
        throw std::runtime_error(std::string("invalid config key '") + key + "': " + e.what());
    }
}
template <class T>
inline std::vector<T> YV(const YAML::Node &n, const char *key)
{
    std::vector<T> v;
    if (!n || !n[key]) return v;
    if (!n[key].IsSequence()) throw std::runtime_error(std::string(key) + " must be a sequence");
    for (const auto &x : n[key])
    {
        T value = x.as<T>();
        if constexpr (std::is_floating_point<T>::value)
            if (!std::isfinite(value)) throw std::runtime_error(std::string(key) + " must be finite");
        v.push_back(value);
    }
    return v;
}

inline void RequireChoice(const std::string &value, const char *key,
                          std::initializer_list<const char *> allowed)
{
    std::string choices;
    for (const char *choice : allowed)
    {
        if (value == choice) return;
        if (!choices.empty()) choices += ", ";
        choices += choice;
    }
    throw std::runtime_error(std::string(key) + " must be one of: " + choices + " (got '" + value + "')");
}

inline void ValidateEdges(const std::vector<double> &edges, const char *name, size_t count = 0)
{
    if (edges.size() < 2 || (count && edges.size() != count))
        throw std::runtime_error(std::string(name) + " has the wrong number of edges");
    for (size_t i = 0; i < edges.size(); ++i)
        if (!std::isfinite(edges[i]) || (i && edges[i] <= edges[i - 1]))
            throw std::runtime_error(std::string(name) + " must contain finite, strictly increasing edges");
}


// One cell's identity. Naming is centralized; a one-cell run keeps legacy keys.
struct BinInfo
{
    int centrality = 0, eta = 0;
    bool auau = false, single = true;
    double centrality_low = 0, centrality_high = 0, eta_low = 0, eta_high = 0;
    std::string Name(const std::string &base, const std::string &legacy_suffix = "_0") const
    {
        return base + (single ? legacy_suffix : "_cent" + std::to_string(centrality) + "_eta" + std::to_string(eta));
    }
    std::string Title() const
    {
        char text[160];
        if (single) std::snprintf(text, sizeof(text), "%.1f < eta < %.1f", eta_low, eta_high);
        else if (auau) std::snprintf(text, sizeof(text), "%g <= centrality < %g%%, %g %s eta < %g",
                                    centrality_low, centrality_high, eta_low, eta == 0 ? "<" : "<=", eta_high);
        else std::snprintf(text, sizeof(text), "pp inclusive, %g %s eta < %g", eta_low, eta == 0 ? "<" : "<=", eta_high);
        return text;
    }
};

class BinLayout
{
public:
    explicit BinLayout(const Cuts &cuts) : auau_(cuts.system == "auau"), centrality_(cuts.centrality_bins), eta_(cuts.eta_bins)
    {
        ValidateEdges(eta_, "eta_bins");
        if (auau_) {
            ValidateEdges(centrality_, "photonjet.centrality_bins");
            if (centrality_.front() < 0 || centrality_.back() > 100)
                throw std::runtime_error("centrality must lie in [0,100]");
        } else if (!centrality_.empty()) throw std::runtime_error("centrality_bins apply only to Au+Au");
    }
    int nCentrality() const { return auau_ ? static_cast<int>(centrality_.size()) - 1 : 1; }
    int nEta() const { return static_cast<int>(eta_.size()) - 1; }
    bool singleCell() const { return nCentrality() == 1 && nEta() == 1; }
    int CentralityBin(double value) const
    {
        if (!auau_) return 0; // pp has an inclusive slot, without centrality selection.
        if (!std::isfinite(value) || value < centrality_.front() || value >= centrality_.back()) return -1;
        return static_cast<int>(std::upper_bound(centrality_.begin(), centrality_.end(), value) - centrality_.begin()) - 1;
    }
    int EtaBin(double value) const
    {
        // Preserve strict exterior cuts; shared interior edges belong to the higher bin.
        if (!std::isfinite(value) || value <= eta_.front() || value >= eta_.back()) return -1;
        return static_cast<int>(std::upper_bound(eta_.begin(), eta_.end(), value) - eta_.begin()) - 1;
    }
    BinInfo cell(int centrality, int eta) const
    {
        if (centrality < 0 || centrality >= nCentrality() || eta < 0 || eta >= nEta())
            throw std::out_of_range("histogram cell outside configured layout");
        return {centrality, eta, auau_, singleCell(), auau_ ? centrality_[centrality] : 0,
                auau_ ? centrality_[centrality + 1] : 0, eta_[eta], eta_[eta + 1]};
    }
    YAML::Node Metadata() const
    {
        YAML::Node node;
        node["format_version"] = 1;
        node["centrality_edges"] = centrality_; node["eta_edges"] = eta_;
        node["n_centrality"] = nCentrality(); node["n_eta"] = nEta();
        node["centrality_policy"] = auau_ ? "low-inclusive-high-exclusive" : "inclusive-pp";
        node["eta_policy"] = "outer-exclusive-interior-to-higher";
        return node;
    }
private:
    bool auau_;
    std::vector<double> centrality_, eta_;
};

inline std::set<int> ReadRunList(const std::string &path)
{
    std::set<int> runs;
    if (path.empty()) return runs;
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open run_list_file " + path);
    std::string line;
    while (std::getline(input, line))
    {
        line = line.substr(0, line.find('#'));
        std::istringstream row(line);
        row >> std::ws;
        if (row.eof()) continue;
        int run;
        std::string extra;
        if (!(row >> run) || run <= 0 || (row >> extra))
            throw std::runtime_error("invalid run-list row in " + path + ": " + line);
        runs.insert(run);
    }
    if (input.bad() || runs.empty()) throw std::runtime_error("empty or unreadable run list " + path);
    return runs;
}

inline void ValidateCuts(const Cuts &c)
{
    RequireChoice(c.system, "photonjet.system", {"pp", "auau"});
    RequireChoice(c.threshold_source, "photonjet.threshold_source", {"formula", "tree"});
    RequireChoice(c.weight_mode, "photonjet.weight_mode", {"stored", "recompute", "sample_map"});
    RequireChoice(c.match_source, "photonjet.match_source", {"barcode", "links"});
    RequireChoice(c.flag_check, "photonjet.flag_check", {"strict", "report", "off"});
    ValidateEdges(c.eta_bins, "eta_bins");
    ValidateEdges(c.pT_bins, "pT_bins");
    ValidateEdges(c.pT_bins_truth, "pT_bins_truth");
    if (c.pT_bins.front() < 0 || c.pT_bins_truth.front() < 0 ||
        c.pT_bins_truth.front() > c.pT_bins.front() || c.pT_bins_truth.back() < c.pT_bins.back())
        throw std::runtime_error("truth pT bins must cover the nonnegative reco pT range");
    if (c.system == "auau")
    {
        ValidateEdges(c.centrality_bins, "photonjet.centrality_bins");
        if (c.centrality_bins.front() < 0 || c.centrality_bins.back() > 100)
            throw std::runtime_error("centrality must lie in [0,100]");
    }
    else if (!c.centrality_bins.empty())
        throw std::runtime_error("centrality_bins apply only to Au+Au");
    if (c.use_topo_iso != 1 && c.use_topo_iso != 2)
        throw std::runtime_error("PhotonJetTrees supports use_topo_iso 1 (R=0.3) or 2 (R=0.4) only");
    if (!(c.vertex_cut > 0) || c.reco_min_ET < 0 || c.min_photon_et < 0 || c.cluster_escale <= 0)
        throw std::runtime_error("invalid vertex cut, minimum ET or energy scale");
    if (c.run_min > 0 && c.run_max > 0 && c.run_min > c.run_max)
        throw std::runtime_error("run_min exceeds run_max");
    for (int bit : c.trigger_used)
        if (bit < 0 || bit >= 64) throw std::runtime_error("trigger_used must contain bits in [0,63]");
    if (c.random_seed == 0) throw std::runtime_error("random_seed must be nonzero for reproducible ROOT draws");
    if (c.apply_trigger_eff_correction && (!(c.trigger_eff_p0 > 0 && c.trigger_eff_p0 <= 1) || c.trigger_eff_beta <= 0))
        throw std::runtime_error("trigger efficiency requires 0 < p0 <= 1 and beta > 0");
    if (c.n_nt_fail < 0 || c.reco_noniso_min_shift < 0 || c.reco_iso_min >= c.reco_noniso_max)
        throw std::runtime_error("invalid isolation window or n_nt_fail");
    if (c.trw_xmin >= c.trw_xmax || c.trw_clamp_min > c.trw_clamp_max)
        throw std::runtime_error("invalid response prior range");
    for (double value : {c.eres_data_p0, c.eres_data_p1, c.eres_data_p2, c.eres_mc_p0, c.eres_mc_p1, c.eres_mc_p2})
        if (value < 0) throw std::runtime_error("resolution coefficients must be nonnegative");
    for (int flag : {c.apply_shower_shape_windows, c.require_prompt_class_2, c.apply_trigger_eff_correction,
                    c.unfold_reweight, c.weta_on, c.wphi_on, c.et1_on, c.et2_on, c.et3_on, c.et4_on,
                    c.e11_to_e33_on, c.e32_to_e35_on, c.bdt_on, c.weta_fail, c.wphi_fail, c.et1_fail,
                    c.e11_to_e33_fail, c.e32_to_e35_fail, c.bdt_fail})
        if (flag != 0 && flag != 1) throw std::runtime_error("analysis switches must be 0 or 1");
}

inline Cuts LoadCuts(const YAML::Node &cfg)
{
    Cuts c;
    if (!cfg.IsMap() || !cfg["analysis"].IsMap() || !cfg["photonjet"].IsMap())
        throw std::runtime_error("configuration needs analysis and photonjet mappings");
    const YAML::Node a = cfg["analysis"];
    const YAML::Node pj = cfg["photonjet"];
    const YAML::Node co = a["common"], ti = a["tight"], nt = a["non_tight"];

    c.system = Y<std::string>(pj, "system", "pp");
    c.threshold_source = Y<std::string>(pj, "threshold_source", c.system == "auau" ? "tree" : "formula");
    c.weight_mode = Y<std::string>(pj, "weight_mode", c.system == "auau" ? "sample_map" : "stored");
    c.match_source = Y<std::string>(pj, "match_source", c.system == "auau" ? "links" : "barcode");
    c.sample_map_file = Y<std::string>(pj, "sample_map", "");
    c.vertex_weight_file = Y<std::string>(pj, "vertex_weight_file", "");
    c.apply_shower_shape_windows = Y<int>(pj, "apply_shower_shape_windows", c.system == "auau" ? 0 : 1);
    c.min_photon_et = Y<double>(pj, "min_photon_et", 0.0);
    c.centrality_bins = YV<double>(pj, "centrality_bins");
    c.require_prompt_class_2 = Y<int>(pj, "count_fragmentation_as_signal", 1);

    c.vertex_cut = Y<double>(a, "vertex_cut", 30.0);
    c.run_list_file = Y<std::string>(a, "run_list_file", "");
    c.runs = ReadRunList(c.run_list_file);
    c.flag_check = Y<std::string>(pj, "flag_check", "strict");
    const int seed = Y<int>(pj, "random_seed", 42);
    if (seed <= 0) throw std::runtime_error("random_seed must be positive");
    c.random_seed = static_cast<unsigned int>(seed);
    c.run_min = Y<int>(a, "run_min", -1);
    c.run_max = Y<int>(a, "run_max", -1);
    if (a["trigger_used"])
    {
        if (a["trigger_used"].IsSequence()) c.trigger_used = YV<int>(a, "trigger_used");
        else c.trigger_used.push_back(a["trigger_used"].as<int>());
    }
    c.apply_trigger_eff_correction = Y<int>(a, "apply_trigger_eff_correction", 0);
    c.trigger_eff_p0 = Y<double>(a, "trigger_eff_p0", 1.0);
    c.trigger_eff_mu = Y<double>(a, "trigger_eff_mu", -2.04);
    c.trigger_eff_beta = Y<double>(a, "trigger_eff_beta", 3.33);

    c.reco_min_ET = Y<double>(a, "reco_min_ET", 5.0);
    c.eta_bins = YV<double>(a, "eta_bins");
    c.pT_bins = YV<double>(a, "pT_bins");
    c.pT_bins_truth = YV<double>(a, "pT_bins_truth");

    c.use_topo_iso = Y<int>(a, "use_topo_iso", 2);
    c.reco_iso_min = Y<double>(a, "reco_iso_min", -20.0);
    c.reco_iso_max_b = Y<double>(a, "reco_iso_max_b", 0.49);
    c.reco_iso_max_s = Y<double>(a, "reco_iso_max_s", 0.037);
    c.reco_noniso_min_shift = Y<double>(a, "reco_noniso_min_shift", 0.8);
    c.reco_noniso_max = Y<double>(a, "reco_noniso_max", 20.0);
    c.mc_iso_scale = Y<double>(a, "mc_iso_scale", 1.0);
    c.mc_iso_shift = Y<double>(a, "mc_iso_shift", 0.0);
    c.truth_iso_max = Y<double>(a, "truth_iso_max", 4.0);

    c.cluster_escale = Y<double>(a, "cluster_escale", 1.0);
    c.cluster_escale_nl_slope = Y<double>(a, "cluster_escale_nl_slope", 0.0);
    c.cluster_escale_nl_ref_ET = Y<double>(a, "cluster_escale_nl_ref_ET", 10.0);
    c.eres_data_p0 = Y<double>(a, "cluster_eres_data_p0", 0.0);
    c.eres_data_p1 = Y<double>(a, "cluster_eres_data_p1", 0.0);
    c.eres_data_p2 = Y<double>(a, "cluster_eres_data_p2", 0.0);
    c.eres_mc_p0 = Y<double>(a, "cluster_eres_mc_p0", 0.0);
    c.eres_mc_p1 = Y<double>(a, "cluster_eres_mc_p1", 0.0);
    c.eres_mc_p2 = Y<double>(a, "cluster_eres_mc_p2", 0.0);

    c.common_e11_over_e33_min = Y<double>(co, "e11_over_e33_min", 0.0);
    c.common_e11_over_e33_max = Y<double>(co, "e11_over_e33_max", 0.98);
    c.common_wr_cogx_bound = Y<double>(co, "wr_cogx_bound", 0.0);
    c.common_weta_cogx_bound = Y<double>(co, "cluster_weta_cogx_bound", 2.0);

    c.t_weta_min = Y<double>(ti, "weta_cogx_min", 0.0);
    c.t_weta_max_b = Y<double>(ti, "weta_cogx_max_b", Y<double>(ti, "weta_cogx_max", 1.0));
    c.t_weta_max_s = Y<double>(ti, "weta_cogx_max_s", 0.0);
    c.t_wphi_min = Y<double>(ti, "wphi_cogx_min", 0.0);
    c.t_wphi_max_b = Y<double>(ti, "wphi_cogx_max_b", Y<double>(ti, "wphi_cogx_max", 1.0));
    c.t_wphi_max_s = Y<double>(ti, "wphi_cogx_max_s", 0.0);
    c.t_et1_min_b = Y<double>(ti, "et1_min_b", Y<double>(ti, "et1_min", 0.5));
    c.t_et1_min_s = Y<double>(ti, "et1_min_s", 0.0);
    c.t_et1_max = Y<double>(ti, "et1_max", 1.0);
    c.t_et2_min = Y<double>(ti, "et2_min", 0.0); c.t_et2_max = Y<double>(ti, "et2_max", 1.0);
    c.t_et3_min = Y<double>(ti, "et3_min", 0.0); c.t_et3_max = Y<double>(ti, "et3_max", 1.0);
    c.t_et4_min = Y<double>(ti, "et4_min", 0.0); c.t_et4_max = Y<double>(ti, "et4_max", 1.0);
    c.t_e11e33_min = Y<double>(ti, "e11_over_e33_min", 0.0); c.t_e11e33_max = Y<double>(ti, "e11_over_e33_max", 1.0);
    c.t_e32e35_min = Y<double>(ti, "e32_over_e35_min", 0.8); c.t_e32e35_max = Y<double>(ti, "e32_over_e35_max", 1.0);
    c.t_bdt_max = Y<double>(ti, "bdt_max", 1.0);
    c.t_bdt_min_intercept = Y<double>(ti, "bdt_min_intercept", Y<double>(ti, "bdt_min", 0.0));
    c.t_bdt_min_slope = Y<double>(ti, "bdt_min_slope", 0.0);

    c.nt_weta_min = Y<double>(nt, "weta_cogx_min", 0.0); c.nt_weta_max = Y<double>(nt, "weta_cogx_max", 1.0);
    c.nt_wphi_min = Y<double>(nt, "wphi_cogx_min", 0.0); c.nt_wphi_max = Y<double>(nt, "wphi_cogx_max", 1.0);
    c.nt_et1_min = Y<double>(nt, "et1_min", 0.6); c.nt_et1_max = Y<double>(nt, "et1_max", 1.0);
    c.nt_et4_min = Y<double>(nt, "et4_min", 0.0); c.nt_et4_max = Y<double>(nt, "et4_max", 1.0);
    c.nt_e11e33_min = Y<double>(nt, "e11_over_e33_min", 0.0); c.nt_e11e33_max = Y<double>(nt, "e11_over_e33_max", 1.0);
    c.nt_e32e35_min = Y<double>(nt, "e32_over_e35_min", 0.8); c.nt_e32e35_max = Y<double>(nt, "e32_over_e35_max", 1.0);
    c.nt_bdt_min_intercept = Y<double>(nt, "bdt_min_intercept", Y<double>(nt, "bdt_min", 0.0));
    c.nt_bdt_min_slope = Y<double>(nt, "bdt_min_slope", 0.0);
    c.nt_bdt_max_intercept = Y<double>(nt, "bdt_max_intercept", Y<double>(nt, "bdt_max", 1.0));
    c.nt_bdt_max_slope = Y<double>(nt, "bdt_max_slope", 0.0);
    c.n_nt_fail = Y<int>(a, "n_nt_fail", 1);
    c.weta_on = Y<int>(a, "weta_on", 1); c.wphi_on = Y<int>(a, "wphi_on", 1); c.et1_on = Y<int>(a, "et1_on", 1);
    c.et2_on = Y<int>(a, "et2_on", 1); c.et3_on = Y<int>(a, "et3_on", 1); c.e11_to_e33_on = Y<int>(a, "e11_to_e33_on", 1);
    c.e32_to_e35_on = Y<int>(a, "e32_to_e35_on", 1); c.et4_on = Y<int>(a, "et4_on", 1); c.bdt_on = Y<int>(a, "bdt_on", 1);
    c.weta_fail = Y<int>(a, "weta_fail", 0); c.wphi_fail = Y<int>(a, "wphi_fail", 0); c.et1_fail = Y<int>(a, "et1_fail", 0);
    c.e11_to_e33_fail = Y<int>(a, "e11_to_e33_fail", 0); c.e32_to_e35_fail = Y<int>(a, "e32_to_e35_fail", 0); c.bdt_fail = Y<int>(a, "bdt_fail", 0);

    const YAML::Node un = a["unfold"];
    c.unfold_reweight = Y<int>(un, "reweight", 0);
    if (un && un["truth_reweight"])
    {
        const YAML::Node tr = un["truth_reweight"];
        c.trw_formula = Y<std::string>(tr, "formula", "1");
        c.trw_params = YV<double>(tr, "params");
        c.trw_xmin = Y<double>(tr, "xmin", 0.0);
        c.trw_xmax = Y<double>(tr, "xmax", 100.0);
        c.trw_clamp_min = Y<double>(tr, "clamp_pT_min", 0.0);
        c.trw_clamp_max = Y<double>(tr, "clamp_pT_max", 1e9);
    }
    if (Y<int>(a, "tower_mask_on", 0) || Y<int>(co, "npb_cut_on", 0) ||
        Y<int>(a, "truth_vertex_reweight_on", 0) || Y<double>(a, "cluster_eres", 0.0) != 0.0)
        throw std::runtime_error("tower mask, NPB, truth-vertex weights and legacy multiplicative smearing are not supported by PhotonJetTrees");
    const std::set<std::string> keys = {"release", "system", "threshold_source", "weight_mode", "match_source",
        "data_files_list", "sim_signal_files_list", "sim_inclusive_files_list", "sample_map", "sample_map_signal",
        "sample_map_inclusive", "vertex_weight_file", "centrality_bins", "apply_shower_shape_windows",
        "count_fragmentation_as_signal", "min_photon_et", "external_mbd_eff_file", "flag_check", "random_seed"};
    for (const auto &entry : pj)
        if (!keys.count(entry.first.as<std::string>()))
            throw std::runtime_error("unknown photonjet key: " + entry.first.as<std::string>());
    ValidateCuts(c);
    return c;
}


} // namespace PJ
