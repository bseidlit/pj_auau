// Photonjet configuration: existing YAML keys, checked once before processing.
#pragma once
#include <yaml-cpp/yaml.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
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

struct Config
{
    // photonjet block
    std::string system = "pp";           // pp | auau
    std::string threshold_source = "formula"; // formula | tree
    int apply_shower_shape_windows = 1;
    int require_prompt_class_2 = 1;      // count_fragmentation_as_signal: class 2 is signal like class 1, as in PPG12
    std::string jet_input;               // jets.jet_input_identity to keep: towerinfo_calib (pp), towerinfo_sub1 (Au+Au)
    int sample_weights = 1;              // simulation: PPG12 sample weights and truth-pT windows from CrossSectionWeights.h
    struct Period { std::string name; double lumi_fraction = 1, single_fraction = 1; std::string reweight_file; };
    std::vector<Period> periods;         // simulation: the run periods of the PPG12 blend, see Job::Weight

    // event level
    double vertex_cut = 30.0;
    int run_min = -1, run_max = -1;
    std::string run_list_file;
    std::set<int> runs;
    // Jets, with the pj_auau reference keys: jet_cone_size 3 is R = 0.3.
    double jet_radius = 0.3, jet_eta_max = 0.8;
    double b2bjet_pT_min = 8.0, b2bjet_dphi = 2.749;   // back-to-back jets of x_Jgamma
    std::vector<double> xjgamma_bins, xjgamma_bins_truth;
    int npb_cut_on = 1;                  // EvtCharacter clusters; default 1 (pp), 0 (Au+Au, where the tree has no score)
    double npb_score_cut = 0.5;
    unsigned int random_seed = 42;
    std::vector<int> trigger_used;       // data: any one of these bits of scaled_trigger_bits
    int apply_trigger_eff_correction = 0;
    double trigger_eff_p0 = 1.0, trigger_eff_mu = -2.04, trigger_eff_beta = 3.33;

    // kinematics
    double reco_min_ET = 5.0;
    std::vector<double> eta_bins, pT_bins, pT_bins_truth;
    std::vector<double> centrality_bins; // Au+Au centrality classes [lo, hi); read but unused in pp

    // isolation
    int cone_size = 4;                   // 3 -> iso_r03, 4 -> iso_r04, for reconstructed and truth isolation
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

// Keys read so far. The analysis block is shared with the PPG12 and pj_auau configs, and the keys
// of theirs that this code does not read are listed once at start-up instead of failing.
inline std::set<std::string> &ReadKeys()
{
    static std::set<std::string> keys;
    return keys;
}
// Read one key. An absent key returns the default, which LoadConfig takes from the Config struct
// so that every default is written once. A malformed supplied value is an error.
template <class T>
inline T Get(const YAML::Node &n, const char *key, T def)
{
    ReadKeys().insert(key);
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
inline std::vector<T> GetList(const YAML::Node &n, const char *key)
{
    ReadKeys().insert(key);
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

inline void ValidateEdges(const std::vector<double> &edges, const char *name)
{
    if (edges.size() < 2) throw std::runtime_error(std::string(name) + " needs at least two edges");
    for (size_t i = 0; i < edges.size(); ++i)
        if (!std::isfinite(edges[i]) || (i && edges[i] <= edges[i - 1]))
            throw std::runtime_error(std::string(name) + " must contain finite, strictly increasing edges");
}


// One centrality/eta cell; naming is independent of the number of cells.
struct BinInfo
{
    int centrality = 0, eta = 0;
    bool auau = false, single = true, split_eta = false;
    double centrality_low = 0, centrality_high = 0, eta_low = 0, eta_high = 0;
    std::string Name(const std::string &base) const
    {
        return base + "_" + std::to_string(centrality) + (split_eta ? "_eta" + std::to_string(eta) : "");
    }
    // The reference names its x_Jgamma histograms base_cent<c>.
    std::string CentName(const std::string &base) const
    {
        return base + "_cent" + std::to_string(centrality) + (split_eta ? "_eta" + std::to_string(eta) : "");
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
    // The edges were checked by ValidateConfig when the configuration was loaded.
    explicit BinLayout(const Config &config) : auau_(config.system == "auau"), centrality_(config.centrality_bins), eta_(config.eta_bins) {}
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
        BinInfo info;
        info.centrality = centrality;
        info.eta = eta;
        info.auau = auau_;
        info.single = singleCell();
        info.split_eta = nEta() > 1;
        if (auau_) {
            info.centrality_low = centrality_[centrality];
            info.centrality_high = centrality_[centrality + 1];
        }
        info.eta_low = eta_[eta];
        info.eta_high = eta_[eta + 1];
        return info;
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

// The values that select a code path or size a histogram. A physics value is the user's business.
inline void ValidateConfig(const Config &c)
{
    if (c.system != "pp" && c.system != "auau") throw std::runtime_error("photonjet.system must be pp or auau");
    if (c.threshold_source != "formula" && c.threshold_source != "tree") throw std::runtime_error("photonjet.threshold_source must be formula or tree");
    ValidateEdges(c.eta_bins, "eta_bins");
    ValidateEdges(c.pT_bins, "pT_bins");
    ValidateEdges(c.pT_bins_truth, "pT_bins_truth");
    if (c.pT_bins_truth.front() > c.pT_bins.front() || c.pT_bins_truth.back() < c.pT_bins.back())
        throw std::runtime_error("truth pT bins must cover the reco pT range");
    if (c.system == "auau")
    {
        ValidateEdges(c.centrality_bins, "centrality_bins");
        if (c.centrality_bins.front() < 0 || c.centrality_bins.back() > 100)
            throw std::runtime_error("centrality must lie in [0,100]");
    }
    if (c.cone_size != 3 && c.cone_size != 4)
        throw std::runtime_error("PhotonJetTrees has isolation cones 3 (R=0.3) and 4 (R=0.4) only");
    ValidateEdges(c.xjgamma_bins, "xjgamma_bins");
    ValidateEdges(c.xjgamma_bins_truth, "xjgamma_bins_truth");
    for (const Config::Period &period : c.periods)
        if (period.lumi_fraction <= 0 || period.single_fraction <= 0 || period.single_fraction > 1 || period.reweight_file.empty())
            throw std::runtime_error("period " + period.name + " needs lumi_fraction > 0, 0 < single_fraction <= 1 and a reweight file");
    if (c.run_min > 0 && c.run_max > 0 && c.run_min > c.run_max)
        throw std::runtime_error("run_min exceeds run_max");
    for (int bit : c.trigger_used)
        if (bit < 0 || bit >= 64) throw std::runtime_error("trigger_used must contain bits in [0,63]");
    if (c.jet_input.empty()) throw std::runtime_error("photonjet.jet_input must name the jet collection");
    if (c.apply_trigger_eff_correction && (!(c.trigger_eff_p0 > 0 && c.trigger_eff_p0 <= 1) || c.trigger_eff_beta <= 0))
        throw std::runtime_error("trigger efficiency requires 0 < p0 <= 1 and beta > 0");
    if (c.trw_xmin >= c.trw_xmax || c.trw_clamp_min > c.trw_clamp_max)
        throw std::runtime_error("invalid response prior range");
    for (int flag : {c.apply_shower_shape_windows, c.require_prompt_class_2, c.apply_trigger_eff_correction,
                    c.unfold_reweight, c.npb_cut_on, c.sample_weights, c.weta_on, c.wphi_on, c.et1_on, c.et2_on, c.et3_on, c.et4_on,
                    c.e11_to_e33_on, c.e32_to_e35_on, c.bdt_on, c.weta_fail, c.wphi_fail, c.et1_fail,
                    c.e11_to_e33_fail, c.e32_to_e35_fail, c.bdt_fail})
        if (flag != 0 && flag != 1) throw std::runtime_error("analysis switches must be 0 or 1");
}

inline Config LoadConfig(const YAML::Node &cfg)
{
    Config c;
    if (!cfg.IsMap() || !cfg["analysis"].IsMap() || !cfg["photonjet"].IsMap())
        throw std::runtime_error("configuration needs analysis and photonjet mappings");
    const YAML::Node a = cfg["analysis"];
    const YAML::Node pj = cfg["photonjet"];
    const YAML::Node co = a["common"], ti = a["tight"], nt = a["non_tight"];

    c.system = Get(pj, "system", c.system);
    c.threshold_source = Get(pj, "threshold_source", c.threshold_source);
    c.apply_shower_shape_windows = Get(pj, "apply_shower_shape_windows", c.apply_shower_shape_windows);
    c.require_prompt_class_2 = Get(pj, "count_fragmentation_as_signal", c.require_prompt_class_2);
    c.jet_input = Get(pj, "jet_input", c.jet_input);
    c.sample_weights = Get(pj, "sample_weights", c.sample_weights);
    for (const YAML::Node &node : pj["periods"]) {
        Config::Period period;
        period.name = Get<std::string>(node, "name", "");
        period.lumi_fraction = Get(node, "lumi_fraction", period.lumi_fraction);
        period.single_fraction = Get(node, "single_fraction", period.single_fraction);
        period.reweight_file = Get<std::string>(node, "truth_vertex_reweight_file", "");
        c.periods.push_back(period);
    }

    c.vertex_cut = Get(a, "vertex_cut", c.vertex_cut);
    c.run_list_file = Get(a, "run_list_file", c.run_list_file);
    c.runs = ReadRunList(c.run_list_file);
    // Jets and x_Jgamma, with the pj_auau reference keys.
    c.jet_radius = 0.1 * Get(a, "jet_cone_size", 3);
    c.jet_eta_max = Get(a, "jet_eta", c.jet_eta_max);
    c.b2bjet_pT_min = Get(a, "b2bjet_pT_min", c.b2bjet_pT_min);
    c.b2bjet_dphi = Get(a, "b2bjet_dphi", c.b2bjet_dphi);
    c.xjgamma_bins = GetList<double>(a, "xjgamma_bins");
    c.xjgamma_bins_truth = GetList<double>(a, "xjgamma_bins_truth");
    if (c.xjgamma_bins_truth.empty()) c.xjgamma_bins_truth = c.xjgamma_bins;
    c.centrality_bins = GetList<double>(a, "centrality_bins");
    c.npb_cut_on = Get(a, "npb_cut_on", c.npb_cut_on);
    c.npb_score_cut = Get(a, "npb_score_cut", c.npb_score_cut);
    const int seed = Get(pj, "random_seed", static_cast<int>(c.random_seed));
    if (seed <= 0) throw std::runtime_error("random_seed must be positive: a zero seed would not be reproducible");
    c.random_seed = static_cast<unsigned int>(seed);
    c.run_min = Get(a, "run_min", c.run_min);
    c.run_max = Get(a, "run_max", c.run_max);
    // One bit or a list of bits.
    if (a["trigger_used"].IsSequence()) c.trigger_used = GetList<int>(a, "trigger_used");
    else if (a["trigger_used"]) c.trigger_used.push_back(Get(a, "trigger_used", 0));
    c.apply_trigger_eff_correction = Get(a, "apply_trigger_eff_correction", c.apply_trigger_eff_correction);
    c.trigger_eff_p0 = Get(a, "trigger_eff_p0", c.trigger_eff_p0);
    c.trigger_eff_mu = Get(a, "trigger_eff_mu", c.trigger_eff_mu);
    c.trigger_eff_beta = Get(a, "trigger_eff_beta", c.trigger_eff_beta);

    c.reco_min_ET = Get(a, "reco_min_ET", c.reco_min_ET);
    c.eta_bins = GetList<double>(a, "eta_bins");
    c.pT_bins = GetList<double>(a, "pT_bins");
    c.pT_bins_truth = GetList<double>(a, "pT_bins_truth");

    c.cone_size = Get(a, "cone_size", c.cone_size);
    c.reco_iso_min = Get(a, "reco_iso_min", c.reco_iso_min);
    c.reco_iso_max_b = Get(a, "reco_iso_max_b", c.reco_iso_max_b);
    c.reco_iso_max_s = Get(a, "reco_iso_max_s", c.reco_iso_max_s);
    c.reco_noniso_min_shift = Get(a, "reco_noniso_min_shift", c.reco_noniso_min_shift);
    c.reco_noniso_max = Get(a, "reco_noniso_max", c.reco_noniso_max);
    c.mc_iso_scale = Get(a, "mc_iso_scale", c.mc_iso_scale);
    c.mc_iso_shift = Get(a, "mc_iso_shift", c.mc_iso_shift);
    c.truth_iso_max = Get(a, "truth_iso_max", c.truth_iso_max);

    c.cluster_escale = Get(a, "cluster_escale", c.cluster_escale);
    c.cluster_escale_nl_slope = Get(a, "cluster_escale_nl_slope", c.cluster_escale_nl_slope);
    c.cluster_escale_nl_ref_ET = Get(a, "cluster_escale_nl_ref_ET", c.cluster_escale_nl_ref_ET);
    c.eres_data_p0 = Get(a, "cluster_eres_data_p0", c.eres_data_p0);
    c.eres_data_p1 = Get(a, "cluster_eres_data_p1", c.eres_data_p1);
    c.eres_data_p2 = Get(a, "cluster_eres_data_p2", c.eres_data_p2);
    c.eres_mc_p0 = Get(a, "cluster_eres_mc_p0", c.eres_mc_p0);
    c.eres_mc_p1 = Get(a, "cluster_eres_mc_p1", c.eres_mc_p1);
    c.eres_mc_p2 = Get(a, "cluster_eres_mc_p2", c.eres_mc_p2);

    c.common_e11_over_e33_min = Get(co, "e11_over_e33_min", c.common_e11_over_e33_min);
    c.common_e11_over_e33_max = Get(co, "e11_over_e33_max", c.common_e11_over_e33_max);
    c.common_wr_cogx_bound = Get(co, "wr_cogx_bound", c.common_wr_cogx_bound);
    c.common_weta_cogx_bound = Get(co, "cluster_weta_cogx_bound", c.common_weta_cogx_bound);

    c.t_weta_min = Get(ti, "weta_cogx_min", c.t_weta_min);
    c.t_weta_max_b = Get(ti, "weta_cogx_max_b", c.t_weta_max_b);
    c.t_weta_max_s = Get(ti, "weta_cogx_max_s", c.t_weta_max_s);
    c.t_wphi_min = Get(ti, "wphi_cogx_min", c.t_wphi_min);
    c.t_wphi_max_b = Get(ti, "wphi_cogx_max_b", c.t_wphi_max_b);
    c.t_wphi_max_s = Get(ti, "wphi_cogx_max_s", c.t_wphi_max_s);
    c.t_et1_min_b = Get(ti, "et1_min_b", c.t_et1_min_b);
    c.t_et1_min_s = Get(ti, "et1_min_s", c.t_et1_min_s);
    c.t_et1_max = Get(ti, "et1_max", c.t_et1_max);
    c.t_et2_min = Get(ti, "et2_min", c.t_et2_min);
    c.t_et2_max = Get(ti, "et2_max", c.t_et2_max);
    c.t_et3_min = Get(ti, "et3_min", c.t_et3_min);
    c.t_et3_max = Get(ti, "et3_max", c.t_et3_max);
    c.t_et4_min = Get(ti, "et4_min", c.t_et4_min);
    c.t_et4_max = Get(ti, "et4_max", c.t_et4_max);
    c.t_e11e33_min = Get(ti, "e11_over_e33_min", c.t_e11e33_min);
    c.t_e11e33_max = Get(ti, "e11_over_e33_max", c.t_e11e33_max);
    c.t_e32e35_min = Get(ti, "e32_over_e35_min", c.t_e32e35_min);
    c.t_e32e35_max = Get(ti, "e32_over_e35_max", c.t_e32e35_max);
    c.t_bdt_max = Get(ti, "bdt_max", c.t_bdt_max);
    c.t_bdt_min_intercept = Get(ti, "bdt_min_intercept", c.t_bdt_min_intercept);
    c.t_bdt_min_slope = Get(ti, "bdt_min_slope", c.t_bdt_min_slope);

    c.nt_weta_min = Get(nt, "weta_cogx_min", c.nt_weta_min);
    c.nt_weta_max = Get(nt, "weta_cogx_max", c.nt_weta_max);
    c.nt_wphi_min = Get(nt, "wphi_cogx_min", c.nt_wphi_min);
    c.nt_wphi_max = Get(nt, "wphi_cogx_max", c.nt_wphi_max);
    c.nt_et1_min = Get(nt, "et1_min", c.nt_et1_min);
    c.nt_et1_max = Get(nt, "et1_max", c.nt_et1_max);
    c.nt_et4_min = Get(nt, "et4_min", c.nt_et4_min);
    c.nt_et4_max = Get(nt, "et4_max", c.nt_et4_max);
    c.nt_e11e33_min = Get(nt, "e11_over_e33_min", c.nt_e11e33_min);
    c.nt_e11e33_max = Get(nt, "e11_over_e33_max", c.nt_e11e33_max);
    c.nt_e32e35_min = Get(nt, "e32_over_e35_min", c.nt_e32e35_min);
    c.nt_e32e35_max = Get(nt, "e32_over_e35_max", c.nt_e32e35_max);
    c.nt_bdt_min_intercept = Get(nt, "bdt_min_intercept", c.nt_bdt_min_intercept);
    c.nt_bdt_min_slope = Get(nt, "bdt_min_slope", c.nt_bdt_min_slope);
    c.nt_bdt_max_intercept = Get(nt, "bdt_max_intercept", c.nt_bdt_max_intercept);
    c.nt_bdt_max_slope = Get(nt, "bdt_max_slope", c.nt_bdt_max_slope);
    c.n_nt_fail = Get(a, "n_nt_fail", c.n_nt_fail);
    c.weta_on = Get(a, "weta_on", c.weta_on);
    c.wphi_on = Get(a, "wphi_on", c.wphi_on);
    c.et1_on = Get(a, "et1_on", c.et1_on);
    c.et2_on = Get(a, "et2_on", c.et2_on);
    c.et3_on = Get(a, "et3_on", c.et3_on);
    c.e11_to_e33_on = Get(a, "e11_to_e33_on", c.e11_to_e33_on);
    c.e32_to_e35_on = Get(a, "e32_to_e35_on", c.e32_to_e35_on);
    c.et4_on = Get(a, "et4_on", c.et4_on);
    c.bdt_on = Get(a, "bdt_on", c.bdt_on);
    c.weta_fail = Get(a, "weta_fail", c.weta_fail);
    c.wphi_fail = Get(a, "wphi_fail", c.wphi_fail);
    c.et1_fail = Get(a, "et1_fail", c.et1_fail);
    c.e11_to_e33_fail = Get(a, "e11_to_e33_fail", c.e11_to_e33_fail);
    c.e32_to_e35_fail = Get(a, "e32_to_e35_fail", c.e32_to_e35_fail);
    c.bdt_fail = Get(a, "bdt_fail", c.bdt_fail);

    const YAML::Node un = a["unfold"];
    c.unfold_reweight = Get(un, "reweight", c.unfold_reweight);
    if (un && un["truth_reweight"])
    {
        const YAML::Node tr = un["truth_reweight"];
        c.trw_formula = Get(tr, "formula", c.trw_formula);
        c.trw_params = GetList<double>(tr, "params");
        c.trw_xmin = Get(tr, "xmin", c.trw_xmin);
        c.trw_xmax = Get(tr, "xmax", c.trw_xmax);
        c.trw_clamp_min = Get(tr, "clamp_pT_min", c.trw_clamp_min);
        c.trw_clamp_max = Get(tr, "clamp_pT_max", c.trw_clamp_max);
    }
    // The photonjet block is ours, so a misspelt key is an error there.
    const std::set<std::string> keys = {"system", "threshold_source", "data_files_list", "sim_signal_files_list",
        "sim_inclusive_files_list", "apply_shower_shape_windows", "count_fragmentation_as_signal",
        "random_seed", "jet_input", "sample_weights", "periods"};
    for (const auto &entry : pj)
        if (!keys.count(entry.first.as<std::string>()))
            throw std::runtime_error("unknown photonjet key: " + entry.first.as<std::string>());
    // The analysis block may come from a PPG12 or pj_auau config. Say once which of its keys are not read.
    std::string ignored;
    for (const std::pair<std::string, YAML::Node> block : {std::pair<std::string, YAML::Node>{"analysis", a}, {"analysis.tight", ti},
                                                            {"analysis.non_tight", nt}, {"analysis.common", co}, {"analysis.unfold", un}})
        for (const auto &entry : block.second)
            if (!ReadKeys().count(entry.first.as<std::string>()) && !entry.second.IsMap())
                ignored += " " + block.first + "." + entry.first.as<std::string>();
    if (!ignored.empty()) std::cout << "photonjet: config keys not read:" << ignored << "\n";
    ValidateConfig(c);
    return c;
}


} // namespace PJ
