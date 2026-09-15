#ifndef CROSSSECTIONWEIGHTS_H
#define CROSSSECTIONWEIGHTS_H

// Single source of truth for MC sample cross-section weights (pb) and
// per-sample kinematic windows.
// Include this header instead of copy-pasting constants per file.

#include <string>

namespace PPG12 {

// Photon-jet samples (Pythia)
constexpr float photon5cross  = 146359.3f;   // canonical window 0-14 GeV (truth max photon pT)
constexpr float photon10cross = 6944.675f;    // canonical window 14-22 GeV
constexpr float photon20cross = 130.4461f;    // canonical window 22+ GeV (photon20 sample is fully on by pT=22)

// Photon-jet samples (Herwig, signal-only cross-check)
// photonjet5 (run 28 type 45) is omitted — generator-level shape mismatch with
// photonjet10 in their overlap region and an unresolved xsec issue. The
// effective truth-photon-pT threshold of photon10_herwig (~7 GeV after FSR
// softening) covers the full PPG12 analysis pT range from 8 GeV upward.
// Truth-pT split point: 22 GeV (mirrors Pythia photon10/photon20 split).
// Cross-check: see reports/herwig_xsec_stitch_check.py + herwig_shape_check.py.
constexpr float herwig_photon10cross = 3.62808e+02f;  // run 28 type 46 (Photonjet pTmin=10 GeV)
constexpr float herwig_photon20cross = 5.34010e+01f;  // run 28 type 47 (Photonjet pTmin=20 GeV)

// QCD jet samples
constexpr float jet5cross  = 1.3878e+08f;
constexpr float jet8cross  = 1.15e+07f;    // restored 2026-09-08: 4.929e6 (used 2026-05-14..09-08) left a x2.5 step
                                             // in the leading-truth-jet spectrum at the jet8/jet12 boundary
constexpr float jet10cross = 3.997e+06f;
constexpr float jet12cross = 1.4903e+06f;
constexpr float jet15cross = 4.073e+05f;
constexpr float jet20cross = 6.2623e+04f;
constexpr float jet30cross = 2.5298e+03f;
constexpr float jet40cross = 1.3553e+02f;
constexpr float jet50cross = 7.3113f;

// Number of generated events per sample (default)
constexpr float default_nsimevents = 1e7f;

// ---------------------------------------------------------------------------
// SampleConfig: per-filetype kinematic windows and cross-section weight.
// Replaces the ~120-line if/else chains duplicated across analysis macros.
// ---------------------------------------------------------------------------
struct SampleConfig {
    float photon_pt_lower = 0;     // truth photon pT window lower bound
    float photon_pt_upper = 0;     // truth photon pT window upper bound
    float jet_pt_lower    = 0;     // truth jet pT window lower bound
    float jet_pt_upper    = 100;   // truth jet pT window upper bound
    float cluster_ET_upper = 100;  // reco cluster ET upper bound
    float weight          = 1.0f;  // cross-section weight (relative to reference sample)
    bool  isbackground    = false; // true for QCD jet samples
    bool  valid           = false; // false if filetype was not recognized
};

// Return kinematic windows and weight for a given filetype string.
// Photon weights are normalized to photon20cross; jet weights to jet50cross.
inline SampleConfig GetSampleConfig(const std::string &filetype)
{
    SampleConfig c;
    c.valid = true;

    // --- Photon samples ---
    if (filetype == "photon5") {
        c.photon_pt_lower = 0;   c.photon_pt_upper = 14;
        c.weight = photon5cross / photon20cross;
    }
    else if (filetype == "photon10") {
        c.photon_pt_lower = 14;  c.photon_pt_upper = 22;
        c.weight = photon10cross / photon20cross;
    }
    else if (filetype == "photon10_double" || filetype == "photon10_nom") {
        c.photon_pt_lower = 14;  c.photon_pt_upper = 22;
        c.weight = photon10cross / photon20cross;
    }
    else if (filetype == "photon20") {
        c.photon_pt_lower = 22;  c.photon_pt_upper = 200;
        c.weight = 1.0f;
    }
    else if (filetype == "photon5_double" || filetype == "photon5_nom") {
        c.photon_pt_lower = 0;   c.photon_pt_upper = 14;
        c.weight = photon5cross / photon20cross;
    }
    else if (filetype == "photon20_double" || filetype == "photon20_nom") {
        c.photon_pt_lower = 22;  c.photon_pt_upper = 200;
        c.weight = 1.0f;
    }
    // --- Herwig photon-jet samples (signal-only cross-check, no Pythia photon5 analog) ---
    // photon10_herwig fills truth pT [0, 22], photon20_herwig fills [22, 200].
    // Weights normalized to herwig_photon20cross (HERWIG-internal reference) so
    // a HERWIG-only run is self-consistent.
    else if (filetype == "photon10_herwig") {
        c.photon_pt_lower = 0;   c.photon_pt_upper = 22;
        c.weight = herwig_photon10cross / herwig_photon20cross;
    }
    else if (filetype == "photon20_herwig") {
        c.photon_pt_lower = 22;  c.photon_pt_upper = 200;
        c.weight = 1.0f;
    }
    // --- Jet samples ---
    else if (filetype == "jet5") {
        c.jet_pt_lower = 7;   c.jet_pt_upper = 9;   c.cluster_ET_upper = 10;
        c.weight = jet5cross / jet50cross;   c.isbackground = true;
    }
    else if (filetype == "jet8") {
        c.jet_pt_lower = 9;   c.jet_pt_upper = 14;  c.cluster_ET_upper = 15;
        c.weight = jet8cross / jet50cross;   c.isbackground = true;
    }
    else if (filetype == "jet8_double" || filetype == "jet8_nom") {
        c.jet_pt_lower = 9;   c.jet_pt_upper = 14;  c.cluster_ET_upper = 15;
        c.weight = jet8cross / jet50cross;   c.isbackground = true;
    }
    else if (filetype == "jet10") {
        c.jet_pt_lower = 10;  c.jet_pt_upper = 15;  c.cluster_ET_upper = 18;
        c.weight = jet10cross / jet50cross;  c.isbackground = true;
    }
    else if (filetype == "jet12") {
        c.jet_pt_lower = 14;  c.jet_pt_upper = 21;  c.cluster_ET_upper = 23;
        c.weight = jet12cross / jet50cross;  c.isbackground = true;
    }
    else if (filetype == "jet12_double" || filetype == "jet12_nom") {
        c.jet_pt_lower = 14;  c.jet_pt_upper = 21;  c.cluster_ET_upper = 23;
        c.weight = jet12cross / jet50cross;  c.isbackground = true;
    }
    else if (filetype == "jet15") {
        c.jet_pt_lower = 15;  c.jet_pt_upper = 21;  c.cluster_ET_upper = 23;
        c.weight = jet15cross / jet50cross;  c.isbackground = true;
    }
    else if (filetype == "jet20") {
        c.jet_pt_lower = 21;  c.jet_pt_upper = 32;  c.cluster_ET_upper = 35;
        c.weight = jet20cross / jet50cross;  c.isbackground = true;
    }
    else if (filetype == "jet20_double" || filetype == "jet20_nom") {
        c.jet_pt_lower = 21;  c.jet_pt_upper = 32;  c.cluster_ET_upper = 35;
        c.weight = jet20cross / jet50cross;  c.isbackground = true;
    }
    else if (filetype == "jet30") {
        c.jet_pt_lower = 32;  c.jet_pt_upper = 42;  c.cluster_ET_upper = 45;
        c.weight = jet30cross / jet50cross;  c.isbackground = true;
    }
    else if (filetype == "jet30_double" || filetype == "jet30_nom") {
        c.jet_pt_lower = 32;  c.jet_pt_upper = 42;  c.cluster_ET_upper = 45;
        c.weight = jet30cross / jet50cross;  c.isbackground = true;
    }
    else if (filetype == "jet40") {
        c.jet_pt_lower = 42;  c.jet_pt_upper = 100; c.cluster_ET_upper = 100;
        c.weight = jet40cross / jet50cross;  c.isbackground = true;
    }
    else if (filetype == "jet40_double" || filetype == "jet40_nom") {
        c.jet_pt_lower = 42;  c.jet_pt_upper = 100; c.cluster_ET_upper = 100;
        c.weight = jet40cross / jet50cross;  c.isbackground = true;
    }
    else if (filetype == "jet50") {
        c.jet_pt_lower = 52;  c.jet_pt_upper = 100;
        c.weight = jet50cross / jet50cross;  c.isbackground = true;
    }
    else {
        c.valid = false;
    }
    return c;
}

}  // namespace PPG12

#endif  // CROSSSECTIONWEIGHTS_H
