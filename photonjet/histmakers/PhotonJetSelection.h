// Named physics decisions for the event loops. No histogram is touched here.
#pragma once
#include "PhotonJetConfig.h"
#include "PhotonJetObjects.h"
#include <TRandom3.h>

namespace PJ {
// ABCD regions. The values index the per-region histogram arrays.
enum Region { kNoRegion = -1, kTightIso = 0, kTightNonIso = 1, kNonTightIso = 2, kNonTightNonIso = 3, kNRegions = 4 };

// One reconstructed photon after SelectPhoton. Only `et` and `eta_bin` are meaningful when not accepted.
struct Candidate
{
    bool accepted = false;      // inside the ET and eta acceptance
    int truth_index = -1;       // copy of Photon::truth_index: position in the event's `truths`, -1 if none
    int eta_bin = -1;
    double et = 0;              // calibrated selection ET, never resolution-smeared
    double eta = 0, phi = 0;
    double iso = 0;             // isolation after the simulation scale and shift
    double weight = 1.;         // event weight, times the trigger-efficiency weight in data
    double smearing_draw = 0;   // copy of Photon::smearing_draw
    bool common = false, tight = false, nontight = false, iso_pass = false, noniso_pass = false;
    int region = kNoRegion;
    bool matched_to_fiducial = false, signal = false;   // set by the truth loop of the spectrum macro
    bool leading = false;                                // the highest-ET candidate of the event, the photon of x_Jgamma
    std::vector<double> b2b_jets, b2b_matched_jets;      // pT of its back-to-back jets, and of those matched to a truth jet
};

inline double SigmaExtraFrac(const Config &c, double pt)
{
    if (pt <= 0) return 0.0;
    const double sd2 = (c.eres_data_p0 * c.eres_data_p0) / pt + (c.eres_data_p1 * c.eres_data_p1) / (pt * pt) + c.eres_data_p2 * c.eres_data_p2;
    const double sm2 = (c.eres_mc_p0 * c.eres_mc_p0) / pt + (c.eres_mc_p1 * c.eres_mc_p1) / (pt * pt) + c.eres_mc_p2 * c.eres_mc_p2;
    return (sd2 > sm2) ? std::sqrt(sd2 - sm2) : 0.0;
}

// Calibration affects selection. Resolution affects the response only.
inline double CalibratedET(const Config &c, double et, bool simulation)
{
    if (!simulation) return et;
    et *= c.cluster_escale;
    return et * (1.0 + c.cluster_escale_nl_slope * (et - c.cluster_escale_nl_ref_ET));
}
inline double TriggerEffWeight(const Config &c, double et)
{
    if (!c.apply_trigger_eff_correction) return 1.0;
    const double eps = c.trigger_eff_p0 * std::exp(-std::exp(-(et - c.trigger_eff_mu) / c.trigger_eff_beta));
    if (!std::isfinite(eps) || eps <= 1e-6)
        throw std::runtime_error("trigger efficiency is zero or invalid at ET=" + std::to_string(et));
    return 1.0 / eps;
}

// Azimuthal difference a - b in (-pi, pi].
inline double DeltaPhi(double a, double b)
{
    const double pi = 3.141592653589793;
    double dphi = a - b;
    while (dphi > pi) dphi -= 2 * pi;
    while (dphi < -pi) dphi += 2 * pi;
    return dphi;
}

// Strictly inside the first and last edge, as the binned spectra are.
inline bool InRange(double value, const std::vector<double> &edges)
{
    return value > edges.front() && value < edges.back();
}
// Bin of `value` in `edges`, counted from 0, or -1 outside. Lower edge inclusive, as the reference counts.
inline int FindBin(double value, const std::vector<double> &edges)
{
    if (value < edges.front() || value >= edges.back()) return -1;
    return std::upper_bound(edges.begin(), edges.end(), value) - edges.begin() - 1;
}
inline double DeltaR(double eta1, double phi1, double eta2, double phi2)
{
    return std::hypot(eta1 - eta2, DeltaPhi(phi1, phi2));
}

// Fill the identification and isolation decisions of `k` from the photon `p`,
// at the calibrated ET already stored in `k.et`.
inline void Classify(const Config &c, bool simulation, const Photon &p, Candidate &k)
{
    const double et = k.et;
    const bool from_tree = (c.threshold_source == "tree");
    const bool cone3 = (c.cone_size == 3);

    // ---- isolation ---------------------------------------------------------
    double iso = cone3 ? p.iso_r03 : p.iso_r04;
    if (simulation) iso = iso * c.mc_iso_scale + c.mc_iso_shift;
    k.iso = iso;
    double iso_max, noniso_min;
    if (from_tree)
    {
        iso_max = cone3 ? p.iso_r03_threshold : p.iso_r04_threshold;
        noniso_min = cone3 ? p.iso_r03_nonisolated_threshold : p.iso_r04_nonisolated_threshold;
    }
    else
    {
        iso_max = c.reco_iso_max_b + c.reco_iso_max_s * et;
        noniso_min = iso_max + c.reco_noniso_min_shift;
    }
    k.iso_pass = (iso > c.reco_iso_min && iso < iso_max);
    k.noniso_pass = (iso > noniso_min && iso < c.reco_noniso_max);

    // ---- BDT ---------------------------------------------------------------
    const double score = p.bdt_score;
    double thr_tight, nt_low, nt_high;
    if (from_tree)
    {
        thr_tight = p.bdt_tight_threshold;
        nt_low = p.bdt_nontight_low_threshold;
        nt_high = p.bdt_nontight_high_threshold;
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
    // ---- common --------------------------------------------------------------
    // prob and NPB are not in the release (documented drop); wr_cogx = wphi/weta.
    const double wr_cogx = (p.weta != 0) ? p.wphi / p.weta : 0.0;
    bool common = true;
    if (c.apply_shower_shape_windows)
    {
        common = p.e11e33 > c.common_e11_over_e33_min && p.e11e33 < c.common_e11_over_e33_max &&
                 wr_cogx > c.common_wr_cogx_bound && p.weta < c.common_weta_cogx_bound;
    }
    k.common = common;
    if (!common) return;

    // ---- tight ---------------------------------------------------------------
    const double t_weta_max = c.t_weta_max_b + c.t_weta_max_s * et;
    const double t_wphi_max = c.t_wphi_max_b + c.t_wphi_max_s * et;
    const double t_et1_min = c.t_et1_min_b + c.t_et1_min_s * et;
    bool w_weta = true, w_wphi = true, w_et1 = true, w_et2 = true, w_et3 = true, w_e11 = true, w_e32 = true, w_et4 = true;
    if (c.apply_shower_shape_windows)
    {
        w_weta = p.weta > c.t_weta_min && p.weta < t_weta_max;
        w_wphi = p.wphi > c.t_wphi_min && p.wphi < t_wphi_max;
        w_et1 = p.et1 > t_et1_min && p.et1 < c.t_et1_max;
        w_et2 = p.et2 > c.t_et2_min && p.et2 < c.t_et2_max;
        w_et3 = p.et3 > c.t_et3_min && p.et3 < c.t_et3_max;
        w_e11 = p.e11e33 > c.t_e11e33_min && p.e11e33 < c.t_e11e33_max;
        w_e32 = p.e32e35 > c.t_e32e35_min && p.e32e35 < c.t_e32e35_max;
        w_et4 = p.et4 > c.t_et4_min && p.et4 < c.t_et4_max;
    }
    k.tight = w_weta && w_wphi && w_et1 && w_et2 && w_et3 && w_e11 && w_e32 && w_et4 && w_bdt;

    // ---- non-tight -----------------------------------------------------------
    bool nt_windows = true;
    if (c.apply_shower_shape_windows)
    {
        nt_windows = p.weta > c.nt_weta_min && p.weta < c.nt_weta_max &&
                     p.wphi > c.nt_wphi_min && p.wphi < c.nt_wphi_max &&
                     p.e11e33 > c.nt_e11e33_min && p.e11e33 < c.nt_e11e33_max &&
                     p.e32e35 > c.nt_e32e35_min && p.e32e35 < c.nt_e32e35_max &&
                     p.et1 > c.nt_et1_min && p.et1 < c.nt_et1_max &&
                     p.et4 > c.nt_et4_min && p.et4 < c.nt_et4_max;
    }
    // Tree semantics: is_nontight == (low <= score < high). PPG12 formula path
    // uses strict inequalities on both sides (identical on continuous scores).
    // Non-tight: more than n_nt_fail of the counted tight cuts fail, and every cut flagged *_fail fails.
    if (nt_windows && in_band)
    {
        const int nfail = !w_weta * c.weta_on + !w_wphi * c.wphi_on + !w_et1 * c.et1_on + !w_et2 * c.et2_on + !w_et3 * c.et3_on +
                          !w_e11 * c.e11_to_e33_on + !w_e32 * c.e32_to_e35_on + !w_et4 * c.et4_on + !w_bdt * c.bdt_on;
        const bool flagged_fail = !(c.weta_fail && w_weta) && !(c.wphi_fail && w_wphi) && !(c.et1_fail && w_et1) &&
                                  !(c.e11_to_e33_fail && w_e11) && !(c.e32_to_e35_fail && w_e32) && !(c.bdt_fail && w_bdt);
        k.nontight = nfail > c.n_nt_fail && flagged_fail;
    }
    if (k.tight) k.nontight = false; // disjoint by construction

    if (k.tight && k.iso_pass) k.region = kTightIso;
    else if (k.tight && k.noniso_pass) k.region = kTightNonIso;
    else if (k.nontight && k.iso_pass) k.region = kNonTightIso;
    else if (k.nontight && k.noniso_pass) k.region = kNonTightNonIso;
}

// Random initialization of one part, done once right after it is read: every reconstructed photon
// gets one standard normal (Photon::smearing_draw), event by event in file order. The response
// smearing uses it later. Because the numbers are fixed before any cut, they depend on the seed
// and the part number only: not on the cuts, and not on how the list is split into jobs.
inline void InitRandom(std::vector<Event> &events, unsigned int random_seed, int part_number)
{
    TRandom3 rng(random_seed * 1000003u + static_cast<UInt_t>(part_number) + 1u);   // never 0: TRandom3(0) seeds from the clock
    for (Event &event : events)
        for (Photon &photon : event.photons) photon.smearing_draw = rng.Gaus(0., 1.);
}

// Data only: the run selection and the trigger. Any one bit of analysis.trigger_used is enough.
inline bool PassesRunAndTrigger(const Config &c, const Event &event)
{
    if (event.run <= 0 || (!c.runs.empty() && !c.runs.count(event.run)) ||
        (c.run_min > 0 && event.run < c.run_min) || (c.run_max > 0 && event.run > c.run_max)) return false;
    if (c.trigger_used.empty()) return true;
    for (int bit : c.trigger_used)
        if ((bit == 30 && event.scaled_bit30 == 1) || ((event.scaled_trigger_bits >> bit) & 1ULL)) return true;
    return false;
}

// Simulation only: is the event inside the truth-pT window of its generator sample? The samples are
// stitched by these windows (CrossSectionWeights.h), so an event outside belongs to another sample and
// is dropped. Photon samples cut on the leading prompt truth photon, jet samples on the leading R = 0.4
// truth jet, as PPG12 does. An event with no such object fails, as in PPG12.
inline bool InSampleWindow(const Event &event)
{
    double leading = -1;
    if (event.jet_sample) {
        for (const Jet &jet : event.truth_jets)
            if (std::fabs(jet.radius - 0.4) < 1e-6) leading = std::max(leading, jet.pt);
    } else {
        for (const TruthPhoton &truth : event.truths)
            if (truth.prompt_class == 1 || truth.prompt_class == 2) leading = std::max(leading, truth.truth_photon_pt);
    }
    return leading >= event.window_low && leading < event.window_high;
}

// The event selection is PassesRunAndTrigger (data) or InSampleWindow (simulation), a centrality bin
// (BinLayout::CentralityBin) and this vertex cut. The event weight comes from Job::Weight.
inline bool HasVertex(const Config &c, const Event &event)
{
    return std::isfinite(event.vertex_z) && std::fabs(event.vertex_z) <= c.vertex_cut;
}

// Acceptance, calibration, weight and classification of one reconstructed photon.
inline Candidate SelectPhoton(const Photon &photon, const Config &c, const BinLayout &bins,
                              bool simulation, double event_weight)
{
    Candidate k;
    k.et = CalibratedET(c, photon.photon_et, simulation);
    k.eta = photon.photon_eta;
    k.phi = photon.photon_phi;
    k.eta_bin = bins.EtaBin(photon.photon_eta);
    if (k.et < c.reco_min_ET || k.eta_bin < 0) return k;
    k.accepted = true;
    k.truth_index = photon.truth_index;
    k.smearing_draw = photon.smearing_draw;
    k.weight = simulation ? event_weight : event_weight * TriggerEffWeight(c, k.et);
    Classify(c, simulation, photon, k);
    return k;
}
inline bool IsFiducialTruth(const TruthPhoton &truth, const Config &c)
{
    return truth.truth_isolation < c.truth_iso_max &&
           (truth.prompt_class == 1 || (c.require_prompt_class_2 && truth.prompt_class == 2));
}
// One entry of the matched response: a tight isolated candidate and its fiducial truth photon.
// `prior_weight` is the truth-spectrum reweight of the response.
struct ResponsePoint { bool accepted = false; double reco = 0, truth = 0, weight = 0; };
inline ResponsePoint MakeResponse(const Config &c, const Candidate &k, const TruthPhoton &truth,
                                  double prior_weight, double event_weight)
{
    ResponsePoint point;
    const double pt = truth.truth_photon_pt;
    if (!k.tight || !k.iso_pass || !InRange(pt, c.pT_bins_truth)) return point;
    const double reco = k.et + SigmaExtraFrac(c, pt) * pt * k.smearing_draw;   // data-simulation resolution difference
    if (!InRange(reco, c.pT_bins)) return point;
    if (!std::isfinite(prior_weight) || prior_weight <= 0) throw std::runtime_error("invalid response prior weight");
    point.accepted = true;
    point.reco = reco;
    point.truth = pt;
    point.weight = event_weight * prior_weight;
    return point;
}
// The jet collection of the analysis: the radius of jet_cone_size and the input of photonjet.jet_input.
inline bool IsConfiguredJet(const Jet &jet, const Config &config)
{
    return jet.configured_input && std::fabs(jet.radius - config.jet_radius) < 1e-6;
}
// pT, smallest first, of the jets of the configured collection that recoil against the direction `phi`:
// |eta| < jet_eta, |dphi| > b2bjet_dphi, pT > b2bjet_pT_min. `matched` keeps only jets with a truth match.
// Works on truth jets too. The jet pT enters as stored: the reference multiplies it by 1/0.65 here, next
// to a comment saying the calibration was removed. That factor is deliberately not applied.
inline std::vector<double> BackToBackJets(const Config &c, const std::vector<Jet> &jets, double phi, bool matched = false)
{
    std::vector<double> pts;
    for (const Jet &jet : jets)
        if (IsConfiguredJet(jet, c) && std::fabs(jet.eta) < c.jet_eta_max && std::fabs(DeltaPhi(phi, jet.phi)) > c.b2bjet_dphi &&
            jet.pt > c.b2bjet_pT_min && (!matched || jet.truth_index >= 0)) pts.push_back(jet.pt);
    std::sort(pts.begin(), pts.end());
    return pts;
}
// Simulation: the reference's geometric jet matching. A truth jet within dR < 0.2 of a truth photon is
// the photon itself and is marked. A reconstructed jet of the configured collection takes the nearest
// other truth jet of the same radius within dR < 0.2 as Jet::truth_index, a position in event.truth_jets.
inline void MatchJets(const Config &c, Event &event)
{
    for (Jet &truth : event.truth_jets)
        for (const TruthPhoton &photon : event.truths)
            if (DeltaR(truth.eta, truth.phi, photon.truth_photon_eta, photon.truth_photon_phi) < 0.2) truth.near_photon = true;
    for (Jet &jet : event.jets) {
        if (!IsConfiguredJet(jet, c)) continue;
        double best = 0.2;
        for (size_t t = 0; t < event.truth_jets.size(); ++t) {
            const Jet &truth = event.truth_jets[t];
            if (truth.near_photon || std::fabs(truth.radius - jet.radius) > 1e-6) continue;
            const double dr = DeltaR(jet.eta, jet.phi, truth.eta, truth.phi);
            if (dr < best) { best = dr; jet.truth_index = t; }
        }
    }
}
} // namespace PJ
