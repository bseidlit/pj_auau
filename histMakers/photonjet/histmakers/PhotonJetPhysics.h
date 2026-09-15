// Named physics decisions for the direct event loop; no histogram mutation.
#pragma once
#include "PhotonJetConfig.h"
#include "PhotonJetReader.h"
#include "../support/CrossSectionWeights.h"
#include <TRandomGen.h>
#include <algorithm>
#include <cstdint>
#include <map>

namespace PJ {
constexpr const char *BackendVersion = "photonjet-cpp-v1";
constexpr const char *RNGVersion = "fnv1a64-splitmix64-mixmax17-gaus-v1";
struct Shapes
{
    double weta, wphi, e11e33, et1, et2, et3, et4, e32e35;
};
struct Decision
{
    double et = 0;        // calibrated selection ET; never resolution-smeared
    double iso = 0;       // isolation value used (after the MC fudge)
    bool common = false, tight = false, nontight = false, iso_pass = false, noniso_pass = false;
    int region = -1;      // 0 A tight-iso, 1 B tight-noniso, 2 C nontight-iso, 3 D nontight-noniso
    // stored-flag cross-check (formula mode only): true when our decision differs
    bool flag_mismatch_tight = false, flag_mismatch_nontight = false, flag_mismatch_iso = false;
};

inline double SigmaExtraFrac(const PJ::Cuts &c, double pt)
{
    if (pt <= 0) return 0.0;
    const double sd2 = (c.eres_data_p0 * c.eres_data_p0) / pt + (c.eres_data_p1 * c.eres_data_p1) / (pt * pt) + c.eres_data_p2 * c.eres_data_p2;
    const double sm2 = (c.eres_mc_p0 * c.eres_mc_p0) / pt + (c.eres_mc_p1 * c.eres_mc_p1) / (pt * pt) + c.eres_mc_p2 * c.eres_mc_p2;
    return (sd2 > sm2) ? std::sqrt(sd2 - sm2) : 0.0;
}

// Calibration affects selection. Resolution affects the response only.
inline double CalibratedET(const PJ::Cuts &c, double et, bool simulation)
{
    if (!simulation) return et;
    et *= c.cluster_escale;
    return et * (1.0 + c.cluster_escale_nl_slope * (et - c.cluster_escale_nl_ref_ET));
}
inline double TriggerEffWeight(const PJ::Cuts &c, double et)
{
    if (!c.apply_trigger_eff_correction) return 1.0;
    const double eps = c.trigger_eff_p0 * std::exp(-std::exp(-(et - c.trigger_eff_mu) / c.trigger_eff_beta));
    if (!std::isfinite(eps) || eps <= 1e-6)
        throw std::runtime_error("trigger efficiency is zero or invalid at ET=" + std::to_string(et));
    return 1.0 / eps;
}

inline Decision Classify(const PJ::Cuts &c, bool issim, double et_modified, const Shapes &s,
                         double score, double tight_threshold, double nt_low_threshold, double nt_high_threshold,
                         double isolation, double iso_threshold, double noniso_threshold,
                         int stored_tight, int stored_nontight, int stored_iso)
{
    Decision d;
    const double et = et_modified;
    d.et = et;
    const bool from_tree = (c.threshold_source == "tree");

    // ---- isolation ---------------------------------------------------------
    double iso = isolation;
    if (issim) iso = iso * c.mc_iso_scale + c.mc_iso_shift;
    d.iso = iso;
    double iso_max, noniso_min;
    if (from_tree)
    {
        iso_max = iso_threshold;
        noniso_min = noniso_threshold;
    }
    else
    {
        iso_max = c.reco_iso_max_b + c.reco_iso_max_s * et;
        noniso_min = iso_max + c.reco_noniso_min_shift;
    }
    d.iso_pass = (iso > c.reco_iso_min && iso < iso_max);
    d.noniso_pass = (iso > noniso_min && iso < c.reco_noniso_max);

    double thr_tight, nt_low, nt_high;
    if (from_tree)
    {
        thr_tight = tight_threshold;
        nt_low = nt_low_threshold;
        nt_high = nt_high_threshold;
    }
    else
    {
        thr_tight = c.t_bdt_min_intercept + c.t_bdt_min_slope * et;
        nt_low = c.nt_bdt_min_intercept + c.nt_bdt_min_slope * et;
        nt_high = c.nt_bdt_max_intercept + c.nt_bdt_max_slope * et;
    }
    const bool w_bdt = score > thr_tight && score < c.t_bdt_max;
    const bool in_band = from_tree ? (score >= nt_low && score < nt_high)
                                   : (score > nt_low && score < nt_high);
    // ---- stored-flag cross-check (pure BDT / iso decisions, no windows) ------
    if (!from_tree && !issim && c.flag_check != "off")
    {
        const bool tree_tight = stored_tight != 0;
        const bool tree_nt = stored_nontight != 0;
        const bool tree_iso = stored_iso != 0;
        d.flag_mismatch_tight = (w_bdt != tree_tight);
        d.flag_mismatch_nontight = (in_band != tree_nt);
        d.flag_mismatch_iso = ((iso < iso_max) != tree_iso);
    }

    // ---- common --------------------------------------------------------------
    // prob and NPB are not in the release (documented drop); wr_cogx = wphi/weta.
    const double wr_cogx = (s.weta != 0) ? s.wphi / s.weta : 0.0;
    bool common = true;
    if (c.apply_shower_shape_windows)
    {
        common = s.e11e33 > c.common_e11_over_e33_min && s.e11e33 < c.common_e11_over_e33_max &&
                 wr_cogx > c.common_wr_cogx_bound && s.weta < c.common_weta_cogx_bound;
    }
    d.common = common;
    if (!common) return d;

    // ---- tight ---------------------------------------------------------------
    const double t_weta_max = c.t_weta_max_b + c.t_weta_max_s * et;
    const double t_wphi_max = c.t_wphi_max_b + c.t_wphi_max_s * et;
    const double t_et1_min = c.t_et1_min_b + c.t_et1_min_s * et;
    bool w_weta = true, w_wphi = true, w_et1 = true, w_et2 = true, w_et3 = true, w_e11 = true, w_e32 = true, w_et4 = true;
    if (c.apply_shower_shape_windows)
    {
        w_weta = s.weta > c.t_weta_min && s.weta < t_weta_max;
        w_wphi = s.wphi > c.t_wphi_min && s.wphi < t_wphi_max;
        w_et1 = s.et1 > t_et1_min && s.et1 < c.t_et1_max;
        w_et2 = s.et2 > c.t_et2_min && s.et2 < c.t_et2_max;
        w_et3 = s.et3 > c.t_et3_min && s.et3 < c.t_et3_max;
        w_e11 = s.e11e33 > c.t_e11e33_min && s.e11e33 < c.t_e11e33_max;
        w_e32 = s.e32e35 > c.t_e32e35_min && s.e32e35 < c.t_e32e35_max;
        w_et4 = s.et4 > c.t_et4_min && s.et4 < c.t_et4_max;
    }
    d.tight = w_weta && w_wphi && w_et1 && w_et2 && w_et3 && w_e11 && w_e32 && w_et4 && w_bdt;

    // ---- non-tight -----------------------------------------------------------
    bool nt_windows = true;
    if (c.apply_shower_shape_windows)
    {
        nt_windows = s.weta > c.nt_weta_min && s.weta < c.nt_weta_max &&
                     s.wphi > c.nt_wphi_min && s.wphi < c.nt_wphi_max &&
                     s.e11e33 > c.nt_e11e33_min && s.e11e33 < c.nt_e11e33_max &&
                     s.e32e35 > c.nt_e32e35_min && s.e32e35 < c.nt_e32e35_max &&
                     s.et1 > c.nt_et1_min && s.et1 < c.nt_et1_max &&
                     s.et4 > c.nt_et4_min && s.et4 < c.nt_et4_max;
    }
    // Tree semantics: is_nontight == (low <= score < high). PPG12 formula path
    // uses strict inequalities on both sides (identical on continuous scores).
    if (nt_windows && in_band)
    {
        int nfail = 0;
        if (!w_weta) nfail += c.weta_on;
        if (!w_wphi) nfail += c.wphi_on;
        if (!w_et1) nfail += c.et1_on;
        if (!w_et2) nfail += c.et2_on;
        if (!w_et3) nfail += c.et3_on;
        if (!w_e11) nfail += c.e11_to_e33_on;
        if (!w_e32) nfail += c.e32_to_e35_on;
        if (!w_et4) nfail += c.et4_on;
        if (!w_bdt) nfail += c.bdt_on;
        if (nfail > c.n_nt_fail)
        {
            bool all_flags_fail = true;
            if (c.weta_fail && w_weta) all_flags_fail = false;
            if (c.wphi_fail && w_wphi) all_flags_fail = false;
            if (c.et1_fail && w_et1) all_flags_fail = false;
            if (c.e11_to_e33_fail && w_e11) all_flags_fail = false;
            if (c.e32_to_e35_fail && w_e32) all_flags_fail = false;
            if (c.bdt_fail && w_bdt) all_flags_fail = false;
            if (all_flags_fail) d.nontight = true;
        }
    }
    if (d.tight) d.nontight = false; // disjoint by construction

    if (d.tight && d.iso_pass) d.region = 0;
    else if (d.tight && d.noniso_pass) d.region = 1;
    else if (d.nontight && d.iso_pass) d.region = 2;
    else if (d.nontight && d.noniso_pass) d.region = 3;

    return d;
}

// Explicit byte encoding and avalanche; independent of std::hash, host endianness,
// input paths and processing order. ULong_t must retain all 64 seed bits.
inline std::uint64_t CandidateSeed(std::uint64_t seed, const std::string &part,
                                  int source, std::uint64_t hi, std::uint64_t lo, int candidate)
{
    std::uint64_t hash = UINT64_C(14695981039346656037);
    auto byte = [&](unsigned char value) { hash = (hash ^ value) * UINT64_C(1099511628211); };
    auto word = [&](std::uint64_t value) { for (int i = 0; i < 8; ++i) { byte(value & 255); value >>= 8; } };
    for (unsigned char ch : std::string("PJRDF-RNG-v1")) byte(ch);
    word(seed); word(part.size());
    for (unsigned char ch : part) byte(ch);
    word(static_cast<std::uint64_t>(source)); word(hi); word(lo); word(static_cast<std::uint64_t>(candidate));
    hash += UINT64_C(0x9e3779b97f4a7c15);
    hash = (hash ^ (hash >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    hash = (hash ^ (hash >> 27)) * UINT64_C(0x94d049bb133111eb);
    hash ^= hash >> 31;
    return hash ? hash : UINT64_C(0x9e3779b97f4a7c15);
}
inline double Gaussian(std::uint64_t seed)
{
    static_assert(sizeof(ULong_t) == 8, "candidate RNG requires a 64-bit ROOT seed type");
    if (!seed) throw std::runtime_error("zero candidate RNG seed");
    TRandomMixMax17 rng(static_cast<ULong_t>(seed));
    return rng.Gaus(0., 1.);
}
inline double ResponseET(const PJ::Cuts &c, double reco, double truth, double gaussian)
{
    return reco + SigmaExtraFrac(c, truth) * truth * gaussian;
}

struct VertexTable {
    std::vector<double> edges, weights;
    double operator()(double z) const {
        if (weights.empty()) return 1.;
        auto bin = std::upper_bound(edges.begin(), edges.end(), z) - edges.begin() - 1;
        bin = std::max<decltype(bin)>(0, std::min<decltype(bin)>(weights.size() - 1, bin));
        const double weight = weights[bin];
        if (!std::isfinite(weight) || weight <= 0) throw std::runtime_error("invalid vertex weight");
        return weight;
    }
};
struct EventDecision { bool pass = false; double weight = 1.; int centrality_bin = -1; };
inline EventDecision SelectEvent(const Cuts &c, const BinLayout &bins, bool simulation,
                                const PPG12::SampleConfig &sample, const VertexTable &vertex,
                                const PartData &data, size_t event_index)
{
    const Event &event = data.events[event_index];
    const double z = event.vertex_z, raw_weight = event.event_weight;
    if (simulation && (!std::isfinite(z) || !std::isfinite(raw_weight)))
        throw std::runtime_error("nonfinite simulation vertex or event weight");
    if (!std::isfinite(z) || std::fabs(z) > c.vertex_cut) return {};
    const int centrality = bins.CentralityBin(event.centrality);
    if (centrality < 0) return {};
    if (!simulation) {
        if (event.run <= 0 || (!c.runs.empty() && !c.runs.count(event.run)) ||
            (c.run_min > 0 && event.run < c.run_min) || (c.run_max > 0 && event.run > c.run_max)) return {};
        bool pass_trigger = c.trigger_used.empty();
        for (int bit : c.trigger_used)
            if ((bit == 30 && event.scaled_bit30 == 1) || ((event.scaled_trigger_bits >> bit) & 1ULL)) pass_trigger = true;
        if (c.system == "auau") pass_trigger = event.scaled_bit22 == 1 && event.minimum_bias_pass == 1;
        return {pass_trigger, 1., centrality};
    }
    double weight = raw_weight;
    if (c.weight_mode == "stored" || c.weight_mode == "recompute") {
        if (raw_weight == 1.) return {}; // Producer sentinel; sample-map mode differs below.
        if (c.weight_mode == "recompute") weight = sample.weight * vertex(z);
    } else {
        double leading = -1.;
        if (sample.isbackground) {
            if (data.jet_count[event_index] > 0) leading = std::max(0., data.leading_jet[event_index]);
            if (!(leading >= sample.jet_pt_lower && leading < sample.jet_pt_upper)) return {};
        } else {
            for (const Long64_t row : data.truth_rows[event_index]) {
                const TruthPhoton &truth = data.truths[row];
                if (truth.prompt_class == 1 || truth.prompt_class == 2) leading = std::max(leading, truth.truth_photon_pt);
            }
            if (!(leading >= sample.photon_pt_lower && leading < sample.photon_pt_upper)) return {};
        }
        if (raw_weight == 1. && z == 0.) return {};
        weight = sample.weight * vertex(z);
    }
    if (!std::isfinite(weight) || weight <= 0) throw std::runtime_error("simulation weight must be finite and positive");
    return {true, weight, centrality};
}

// This validation runs before event selection, including for rejected events.
inline std::vector<int> MatchPhotons(const PartData &data, size_t event, bool simulation, bool use_links)
{
    const std::vector<Long64_t> &reco = data.reco_rows[event];
    std::vector<int> match(reco.size(), -1);
    if (!simulation) return match;
    if (use_links) return data.links[event];
    std::map<int, int> barcodes;
    const std::vector<Long64_t> &truth = data.truth_rows[event];
    for (size_t t = 0; t < truth.size(); ++t) {
        const int barcode = data.truths[truth[t]].generator_barcode;
        if (!barcodes.emplace(barcode, static_cast<int>(t)).second) barcodes[barcode] = -2;
    }
    for (size_t p = 0; p < reco.size(); ++p) {
        const int barcode = data.photons[reco[p]].truth_barcode;
        if (barcode < 0) continue;
        const auto found = barcodes.find(barcode);
        if (found == barcodes.end()) continue;
        if (found->second == -2) throw std::runtime_error("ambiguous truth barcode within an event");
        match[p] = found->second;
    }
    return match;
}

struct Candidate
{
    bool accepted = false, matched_to_fiducial = false;
    int original_index = -1, truth_index = -1, eta_bin = -1;
    Decision decision;
    double weight = 1.;
};
inline Candidate SelectPhoton(const Photon &photon, const Cuts &c, const BinLayout &bins,
                              bool simulation, const PPG12::SampleConfig &sample,
                              double event_weight, int truth_index)
{
    const double et = CalibratedET(c, photon.photon_et, simulation);
    if (!std::isfinite(et)) throw std::runtime_error("nonfinite calibrated photon ET");
    if (!simulation && c.min_photon_et > 0 && photon.photon_et < c.min_photon_et)
        throw std::runtime_error("data violate the configured PhotonJetTrees minimum photon ET");
    const int eta = bins.EtaBin(photon.photon_eta);
    if (et < c.reco_min_ET || eta < 0) return {};
    if (simulation && c.weight_mode != "stored" && sample.isbackground && et > sample.cluster_ET_upper) return {};
    const bool cone3 = c.use_topo_iso == 1;
    const Decision decision = Classify(c, simulation, et,
        {photon.weta, photon.wphi, photon.e11e33, photon.et1, photon.et2, photon.et3, photon.et4, photon.e32e35},
        photon.bdt_score, photon.bdt_tight_threshold, photon.bdt_nontight_low_threshold, photon.bdt_nontight_high_threshold,
        cone3 ? photon.iso_r03 : photon.iso_r04, cone3 ? photon.iso_r03_threshold : photon.iso_r04_threshold,
        cone3 ? photon.iso_r03_nonisolated_threshold : photon.iso_r04_nonisolated_threshold,
        photon.bdt_is_tight, photon.bdt_is_nontight, cone3 ? photon.iso_r03_pass : photon.iso_r04_pass);
    const double weight = simulation ? event_weight : event_weight * TriggerEffWeight(c, et);
    return {true, false, photon.original_index, truth_index, eta, decision, weight};
}
inline bool IsFiducialTruth(const TruthPhoton &truth, const Cuts &c)
{
    return truth.truth_isolation < c.truth_iso_max &&
           (truth.prompt_class == 1 || (c.require_prompt_class_2 && truth.prompt_class == 2));
}
struct ResponsePoint
{
    bool accepted = false;
    std::uint64_t seed = 0;
    double gaussian = 0, reco = 0, truth = 0, weight = 0;
};
template<class Prior>
inline ResponsePoint MakeResponse(const Cuts &c, const Candidate &candidate, const TruthPhoton &truth,
                                  const Event &event, const std::string &part, double event_weight, const Prior &prior)
{
    const Decision &d = candidate.decision;
    const double pt = truth.truth_photon_pt;
    if (!d.tight || !d.iso_pass || !(pt > c.pT_bins_truth.front() && pt < c.pT_bins_truth.back())) return {};
    const std::uint64_t seed = CandidateSeed(c.random_seed, part, event.source_file_index, event.event_id_hi, event.event_id_lo,
                                    candidate.original_index);
    const double gaussian = Gaussian(seed), reco = ResponseET(c, d.et, pt, gaussian);
    if (!(reco > c.pT_bins.front() && reco < c.pT_bins.back())) return {};
    const double reweight = prior(pt);
    if (!std::isfinite(reweight) || reweight <= 0) throw std::runtime_error("invalid response prior weight");
    return {true, seed, gaussian, reco, pt, event_weight * reweight};
}
} // namespace PJ
