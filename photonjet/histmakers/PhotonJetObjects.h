// What the analysis sees: one Event that owns its photons, truth photons and jets.
// The input stores these in separate flat trees. PhotonJetReader.h rebuilds the events.
#pragma once
#include <RtypesCore.h>
#include <limits>
#include <tuple>
#include <vector>

namespace PJ {
// Every row of every tree carries the key of the event it belongs to.
using EventKey = std::tuple<Int_t, ULong64_t, ULong64_t>;
struct Identity {
    Int_t source_file_index = -1;
    ULong64_t event_id_hi = 0, event_id_lo = 0;
    EventKey key() const { return {source_file_index, event_id_hi, event_id_lo}; }
};

// Members named after a branch are read from that branch.
struct Photon : Identity {
    Double_t photon_et = 0, photon_eta = 0, photon_phi = 0;
    Double_t bdt_score = 0, bdt_tight_threshold = 0;
    Double_t bdt_nontight_low_threshold = 0, bdt_nontight_high_threshold = 0;
    Double_t iso_r03 = 0, iso_r03_threshold = 0, iso_r03_nonisolated_threshold = 0;
    Double_t iso_r04 = 0, iso_r04_threshold = 0, iso_r04_nonisolated_threshold = 0;
    Double_t weta = 0, wphi = 0, e11e33 = 0, et1 = 0, et2 = 0, et3 = 0, et4 = 0, e32e35 = 0;   // shower_* branches
    Double_t npb_score = 0;                                          // filled in pp, NaN in Au+Au
    Int_t shower_center_eta_index = -1, shower_center_phi_index = -1;   // tower indices of the cluster centre
    // Not branches. Set after the rows are read.
    Int_t truth_index = -1;       // position of the matched truth photon in the event's `truths`, -1 if none
    Double_t smearing_draw = 0;   // standard normal for the response smearing, set by InitRandom in signal simulation
};
struct TruthPhoton : Identity {
    Double_t truth_photon_pt = 0, truth_photon_eta = 0, truth_photon_phi = 0;
    Double_t truth_isolation = 0;   // truth_isolation_r03 or _r04, the cone of the reconstructed isolation
    Int_t prompt_class = 0;
};
// A reconstructed jet (jets tree) or a truth jet (truthJets tree). The trees hold every radius and, in
// Au+Au, two inputs per radius. All of them are kept, in row order, because the link indices count all
// of them. The analysis picks its collection with IsConfiguredJet.
struct Jet {
    Double_t pt = 0, eta = 0, phi = 0, radius = 0;
    // Not branches. See MatchJets.
    Int_t truth_index = -1;         // reconstructed jet: position of the matched truth jet in the event's `truth_jets`, -1 if none
    bool near_photon = false;       // truth jet: within dR < 0.2 of a truth photon, so it is the photon itself
    bool configured_input = true;   // reconstructed jet: its jet_input_identity equals photonjet.jet_input
};

struct Event : Identity {
    Int_t run = 0;
    ULong64_t scaled_trigger_bits = 0;
    Int_t scaled_bit30 = -1;
    Double_t vertex_z = 0, centrality = 0, total_calo_energy = 0;
    Double_t emcal_total_energy = std::numeric_limits<double>::quiet_NaN();
    // Simulation only: the truth vertices of the hard scattering and, in double-interaction samples, of
    // the second collision. The producer's weight columns are not read, see Job::Weight.
    Double_t truth_vertex_z = std::numeric_limits<double>::quiet_NaN();
    Double_t truth_mb_vertex_z = std::numeric_limits<double>::quiet_NaN();
    // Simulation only, not branches: the generator sample of the event's source file, from the part's
    // source_files metadata and CrossSectionWeights.h. Unknown samples get weight 1 and an open window.
    Double_t sample_weight = 1;                       // cross section relative to the reference sample
    Double_t window_low = -1e9, window_high = 1e9;    // truth-pT window that keeps the samples from overlapping
    Double_t cluster_et_upper = 1e9;                  // jet samples: reconstructed photons above this belong to the next sample
    bool jet_sample = false;                          // the window is on the leading truth jet, else on the leading prompt truth photon
    bool double_interaction = false;                  // a double-interaction sample

    // The objects of this event, in the order of their rows in the input trees.
    std::vector<Photon> photons;
    std::vector<TruthPhoton> truths;   // simulation only
    std::vector<Jet> jets;
    std::vector<Jet> truth_jets;       // simulation only
};
} // namespace PJ
