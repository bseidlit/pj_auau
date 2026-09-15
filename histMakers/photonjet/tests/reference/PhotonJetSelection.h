// Photonjet event and photon decisions. No file opening or histogram booking.
#pragma once
#include "PhotonJetConfig.h"
#include "PhotonJetReader.h"
#include <TRandom3.h>
#include <cmath>

namespace PJ
{
struct InputMap
{
    int et = 0, weta = 1, wphi = 2, vz = 3, eta = 4, e11e33 = 5, et1 = 6, et2 = 7, et3 = 8, et4 = 9, e32e35 = 10;
    int weta33 = -1, wphi33 = -1, cent = -1;
};
inline InputMap InputMapFor(const std::string &system)
{
    InputMap m;
    if (system == "auau")
    {
        m.et = 0; m.weta = 1; m.wphi = 2; m.weta33 = 3; m.wphi33 = 4; m.vz = 5; m.eta = 6;
        m.e11e33 = 7; m.et1 = 8; m.et2 = 9; m.et3 = 10; m.et4 = 11; m.e32e35 = 12; m.cent = 13;
    }
    return m;
}

struct Shapes
{
    double weta, wphi, e11e33, et1, et2, et3, et4, e32e35;
};
inline Shapes ShapesOf(const Photon &p, const InputMap &m)
{
    return Shapes{p.in[m.weta], p.in[m.wphi], p.in[m.e11e33], p.in[m.et1], p.in[m.et2], p.in[m.et3], p.in[m.et4], p.in[m.e32e35]};
}

struct Decision
{
    double et = 0;        // calibrated selection ET; never resolution-smeared
    double iso = 0;       // isolation value used (after the MC fudge)
    bool common = false, tight = false, nontight = false, iso_pass = false, noniso_pass = false;
    int region = -1;      // 0 A tight-iso, 1 B tight-noniso, 2 C nontight-iso, 3 D nontight-noniso
    // stored-flag cross-check (formula mode only): true when our decision differs
    bool flag_mismatch_tight = false, flag_mismatch_nontight = false, flag_mismatch_iso = false;
};

inline double SigmaExtraFrac(const Cuts &c, double pt)
{
    if (pt <= 0) return 0.0;
    const double sd2 = (c.eres_data_p0 * c.eres_data_p0) / pt + (c.eres_data_p1 * c.eres_data_p1) / (pt * pt) + c.eres_data_p2 * c.eres_data_p2;
    const double sm2 = (c.eres_mc_p0 * c.eres_mc_p0) / pt + (c.eres_mc_p1 * c.eres_mc_p1) / (pt * pt) + c.eres_mc_p2 * c.eres_mc_p2;
    return (sd2 > sm2) ? std::sqrt(sd2 - sm2) : 0.0;
}

// Calibration affects selection. Resolution affects the response only.
inline double CalibratedET(const Cuts &c, double et, bool simulation)
{
    if (!simulation) return et;
    et *= c.cluster_escale;
    return et * (1.0 + c.cluster_escale_nl_slope * (et - c.cluster_escale_nl_ref_ET));
}
inline double ResponseET(const Cuts &c, double reco_et, double truth_pt, TRandom3 &rng)
{
    const double sigma = SigmaExtraFrac(c, truth_pt) * truth_pt;
    return reco_et + (sigma > 0 ? rng.Gaus(0, sigma) : 0.0);
}

inline double TriggerEffWeight(const Cuts &c, double et)
{
    if (!c.apply_trigger_eff_correction) return 1.0;
    const double eps = c.trigger_eff_p0 * std::exp(-std::exp(-(et - c.trigger_eff_mu) / c.trigger_eff_beta));
    if (!std::isfinite(eps) || eps <= 1e-6)
        throw std::runtime_error("trigger efficiency is zero or invalid at ET=" + std::to_string(et));
    return 1.0 / eps;
}

inline Decision Classify(const Photon &p, const Cuts &c, const InputMap &m, bool issim, double et_modified)
{
    Decision d;
    const double et = et_modified;
    d.et = et;
    const Shapes s = ShapesOf(p, m);
    const bool from_tree = (c.threshold_source == "tree");

    // ---- isolation ---------------------------------------------------------
    double iso = (c.use_topo_iso == 1) ? p.iso3 : p.iso4;
    if (issim) iso = iso * c.mc_iso_scale + c.mc_iso_shift;
    d.iso = iso;
    double iso_max, noniso_min;
    if (from_tree)
    {
        iso_max = (c.use_topo_iso == 1) ? p.iso3_thr : p.iso4_thr;
        noniso_min = (c.use_topo_iso == 1) ? p.iso3_nonthr : p.iso4_nonthr;
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
        thr_tight = p.thr_tight;
        nt_low = p.thr_nt_low;
        nt_high = p.thr_nt_high;
    }
    else
    {
        thr_tight = c.t_bdt_min_intercept + c.t_bdt_min_slope * et;
        nt_low = c.nt_bdt_min_intercept + c.nt_bdt_min_slope * et;
        nt_high = c.nt_bdt_max_intercept + c.nt_bdt_max_slope * et;
    }
    const bool w_bdt = p.score > thr_tight && p.score < c.t_bdt_max;
    const bool in_band = from_tree ? (p.score >= nt_low && p.score < nt_high)
                                   : (p.score > nt_low && p.score < nt_high);
    // ---- stored-flag cross-check (pure BDT / iso decisions, no windows) ------
    if (!from_tree && !issim && c.flag_check != "off")
    {
        const bool tree_tight = p.is_tight != 0;
        const bool tree_nt = p.is_nontight != 0;
        const bool tree_iso = ((c.use_topo_iso == 1) ? p.iso3_pass : p.iso4_pass) != 0;
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

// ---------------------------------------------------------------------------
// Event-level selection.
// ---------------------------------------------------------------------------
inline bool PassTrigger(const Event &e, const Cuts &c)
{
    if (c.system == "auau")
    {
        // Used for data only; simulation does not require an MB trigger.
        return e.scaled_bit22 == 1 && e.minimum_bias_pass == 1;
    }
    if (c.trigger_used.empty()) return true;
    for (int bit : c.trigger_used)
    {
        if (bit == 30 && e.scaled_bit30 == 1) return true;
        if (bit >= 0 && bit < 64 && ((e.scaled_trigger_bits >> bit) & 1ULL)) return true;
    }
    return false;
}

inline bool InCentralityClass(double cent, const Cuts &c)
{
    if (c.centrality_bins.size() < 2) return true;
    return cent >= c.centrality_bins.front() && cent < c.centrality_bins.back();
}

inline bool PassEventData(const Event &e, const Cuts &c)
{
    if (e.run <= 0 || !std::isfinite(e.vertex_z)) return false;
    if (!c.runs.empty() && !c.runs.count(e.run)) return false;
    if (c.run_min > 0 && e.run < c.run_min) return false;
    if (c.run_max > 0 && e.run > c.run_max) return false;
    if (!PassTrigger(e, c)) return false;
    if (std::fabs(e.vertex_z) > c.vertex_cut) return false;
    if (c.system == "auau" && !InCentralityClass(e.centrality, c)) return false;
    return true;
}

// Truth fiducial definition (PPG12: pid 22, photonclass in {1,2}, iso_R03 < 4).
inline bool IsFiducialTruth(const TruthPhoton &t, const Cuts &c)
{
    if (t.prompt_class != 1 && !(c.require_prompt_class_2 && t.prompt_class == 2)) return false;
    if (!(t.iso < c.truth_iso_max)) return false;
    return true;
}


} // namespace PJ
