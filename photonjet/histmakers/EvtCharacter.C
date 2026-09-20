// Event characterization. Follows pj_auau/histMakers/EvtCharacter.C: the same histograms, binning and
// selections. A comment starting with "v21:" marks a difference that the input trees force.
#include "PhotonJetIO.h"
#include "PhotonJetHistograms.h"
#include <TString.h>

// Thresholds of the reference macro.
namespace Evt {
const double kClusterEtMin = 10.0;      // GeV: every cluster histogram
const double kJetPtMinDphi = 10.0;      // GeV: jets in |dphi|(cluster, jet)
const double kJetPtMinEtaPhi = 8.0;     // GeV: jets in h_jet_eta_phi
const double kClusterJetDrMin = 0.2;    // a jet closer than this to the cluster is the cluster itself
const int kRunsPerGroup = 5;            // h_clusters_per_5runs
}

/////////////////////////////////////
// Histograms
/////////////////////////////////////
struct EventHistograms
{
    // Events
    TH1D *vertex_before = nullptr;               // h_vertexz_no_vertexcut_unweighted
    TH1D *vertex = nullptr;                      // h_vertexz
    TH1D *centrality = nullptr;                  // h_centrality
    TH1D *centrality_fine = nullptr;             // h_cent_fine
    TH1D *centrality_with_cluster = nullptr;     // h_centrality_with_cluster10
    TH1D *psi2 = nullptr;                        // h_Psi2, booked but not filled
    TH1D *emcal = nullptr;                       // h_totalEMCal_energy, unweighted
    TH1D *emcal_weighted = nullptr;              // h_totalEMCal_energy_weight
    TH2D *emcal_vs_centrality = nullptr;         // h_totalEMCal_energy_vs_cent, unweighted
    TH1D *mc_weight = nullptr;                   // h_mc_full_weight
    // Clusters with ET >= 10 GeV
    TH1D *cluster_et = nullptr;                  // h_cluster_pt_above10
    TH2D *cluster_eta_phi = nullptr;             // h_cluster_eta_phi
    TH2D *cluster_ieta_iphi = nullptr;           // h_cluster_ieta_iphi
    std::vector<TH2D *> iso_vs_centrality;            // h_iso_vs_cent_pt<i>, one per pT bin
    std::vector<TH1D *> cluster_yield_vs_centrality;  // h_cluster_yield_vs_cent_pt<i>
    std::vector<std::vector<TH1D *>> dphi_cluster_jets;   // h_dphi_clusterJets_cent<c>_pt<i>
    // Jets, in events with such a cluster
    TH1D *jet_pt = nullptr;                      // h_jet_pt
    TH2D *jet_eta_phi = nullptr;                 // h_jet_eta_phi
    TH2D *psi2_vs_jet_pt = nullptr;              // h_Psi2_vs_jetPt, booked but not filled
    // Per run. The run axis comes from analysis.run_min and analysis.run_max, one bin per run.
    TH1D *clusters_per_run = nullptr;            // h_clusters_per_run
    TH1D *clusters_per_5runs = nullptr;          // h_clusters_per_5runs
    TProfile *emcal_per_run = nullptr;           // h_avg_totalEMCal_energy_per_run
    TH1D *scaler_per_run = nullptr;              // h_scaler_per_run, booked but not filled
    TH1D *scaler_per_5runs = nullptr;            // h_scaler_per_5runs, booked but not filled
    // Not in the reference
    TH1D *calo_energy = nullptr;                 // calo_energy: stored total calorimeter energy
    TH1D *photon_multiplicity = nullptr;         // photon_multiplicity: clusters above per event
    TH2D *isolation = nullptr;                   // isolation: isolation ET vs cluster ET
    TH2D *bdt_score = nullptr;                   // bdt_score: BDT score vs cluster ET
};

EventHistograms BookEventHistograms(PJ::HistogramList &booked, const PJ::Config &config, const PJ::BinLayout &bins)
{
    using PJ::UniformBins;
    const double pi = 3.141592653589793;
    const std::vector<double> &pt = config.pT_bins;
    // v21: pp has no centrality bins. The reference books them from its config for both systems.
    const std::vector<double> centrality_edges = config.system == "auau" ? config.centrality_bins : UniformBins(100, 0, 100);
    if (config.run_min <= 0 || config.run_max < config.run_min)
        throw std::runtime_error("EvtCharacter needs analysis.run_min and analysis.run_max for the per-run histograms");
    const int n_runs = config.run_max - config.run_min + 1;
    const int n_groups = (n_runs + Evt::kRunsPerGroup - 1) / Evt::kRunsPerGroup;
    const std::vector<double> runs = UniformBins(n_runs, config.run_min, config.run_max + 1);
    const std::vector<double> run_groups = UniformBins(n_groups, config.run_min, config.run_min + n_groups * Evt::kRunsPerGroup);

    EventHistograms h;
    h.vertex_before = booked.Book1D("h_vertexz_no_vertexcut_unweighted", "Vertex Z (before |z| vertex cut, unweighted);Vertex Z [cm];Events", UniformBins(200, -100, 100));
    h.vertex = booked.Book1D("h_vertexz", "Vertex Z;Vertex Z [cm];Events", UniformBins(200, -100, 100));
    h.centrality = booked.Book1D("h_centrality", "Centrality distribution;Centrality [%];Events", centrality_edges);
    h.centrality_fine = booked.Book1D("h_cent_fine", " ", UniformBins(50, 0.001, 100.001));
    h.centrality_with_cluster = booked.Book1D("h_centrality_with_cluster10", "Centrality (events with at least one cluster E_{T} > 10 GeV);Centrality [%];Events", centrality_edges);
    h.psi2 = booked.Book1D("h_Psi2", "Event plane #Psi_{2};#Psi_{2} [rad];Events", UniformBins(100, -pi, pi));
    h.emcal = booked.Book1D("h_totalEMCal_energy", "Total EMCal energy;E_{tot}^{EMCal} [GeV];events", UniformBins(200, 0, 2000));
    h.emcal_weighted = booked.Book1D("h_totalEMCal_energy_weight", "Total EMCal energy;E_{tot}^{EMCal} [GeV];events", UniformBins(200, 0, 2000));
    h.emcal_vs_centrality = booked.Book2D("h_totalEMCal_energy_vs_cent", "Total EMCal energy vs centrality (unweighted);Centrality [%];E_{tot}^{EMCal} [GeV]", UniformBins(100, 0, 100), UniformBins(200, 0, 2000));
    h.mc_weight = booked.Book1D("h_mc_full_weight", "Per-event MC weight;w_{MC};events", UniformBins(1000, 0, 1e4));

    h.cluster_et = booked.Book1D("h_cluster_pt_above10", "Cluster E_{T} (E_{T} > 10 GeV);Cluster E_{T} [GeV];Clusters", UniformBins(80, 10, 50));
    h.cluster_eta_phi = booked.Book2D("h_cluster_eta_phi", ";#eta;#phi [rad]", UniformBins(50, -1, 1), UniformBins(64, -pi, pi));
    h.cluster_ieta_iphi = booked.Book2D("h_cluster_ieta_iphi", ";i#eta;i#phi", UniformBins(96, -0.5, 95.5), UniformBins(256, -0.5, 255.5));
    h.dphi_cluster_jets.resize(bins.nCentrality());
    for (size_t i = 0; i + 1 < pt.size(); ++i) {
        const std::string range = Form("%.0f < E_{T} < %.0f GeV", pt[i], pt[i + 1]);
        h.iso_vs_centrality.push_back(booked.Book2D("h_iso_vs_cent_pt" + std::to_string(i), "Iso E_{T} vs centrality (" + range + ");Centrality [%];Iso E_{T} [GeV]",
                                                    UniformBins(100, 0, 100), UniformBins(100, -50, 50)));
        h.cluster_yield_vs_centrality.push_back(booked.Book1D("h_cluster_yield_vs_cent_pt" + std::to_string(i), "Cluster yield vs centrality (" + range + ");Centrality [%];Clusters",
                                                              UniformBins(50, 0.001, 100.001)));
        for (int c = 0; c < bins.nCentrality(); ++c)
            h.dphi_cluster_jets[c].push_back(booked.Book1D("h_dphi_clusterJets_cent" + std::to_string(c) + "_pt" + std::to_string(i),
                                                           "|#Delta#phi|(cluster, jet), " + range + ";|#Delta#phi|;Pairs", UniformBins(64, 0, pi)));
    }

    h.jet_pt = booked.Book1D("h_jet_pt", "Jet p_{T};Jet p_{T} [GeV];Jets", UniformBins(50, 0, 100));
    h.jet_eta_phi = booked.Book2D("h_jet_eta_phi", "Jet yield (p_{T} > 8 GeV);#eta;#phi [rad]", UniformBins(10, -1, 1), UniformBins(10, -pi, pi));
    h.psi2_vs_jet_pt = booked.Book2D("h_Psi2_vs_jetPt", "#Psi_{2} vs jet p_{T};Jet p_{T} [GeV];#Psi_{2} [rad]", UniformBins(50, 0, 100), UniformBins(100, -pi, pi));

    h.clusters_per_run = booked.Book1D("h_clusters_per_run", ";Runs;N clusters", runs);
    h.clusters_per_5runs = booked.Book1D("h_clusters_per_5runs", ";First run # in 5-run group;N clusters (sum of 5 runs)", run_groups);
    h.emcal_per_run = booked.BookProfile("h_avg_totalEMCal_energy_per_run", ";Run;#LT E_{tot}^{EMCal} #GT [GeV]", runs);
    h.scaler_per_run = booked.Book1D("h_scaler_per_run", ";Runs;Max MB scaler", runs);
    h.scaler_per_5runs = booked.Book1D("h_scaler_per_5runs", ";First run # in 5-run group;Total MB scaler (sum of 5 runs)", run_groups);

    h.calo_energy = booked.Book1D("calo_energy", ";Stored total calorimeter energy;Weighted selected events", UniformBins(200, 0, 2000));
    h.photon_multiplicity = booked.Book1D("photon_multiplicity", ";Clusters with E_{T} > 10 GeV;Weighted selected events", UniformBins(51, -.5, 50.5));
    h.isolation = booked.Book2D("isolation", ";Cluster E_{T} [GeV];Isolation E_{T} [GeV]", pt, UniformBins(100, -20, 30));
    h.bdt_score = booked.Book2D("bdt_score", ";Cluster E_{T} [GeV];BDT score", pt, UniformBins(100, 0, 1));
    return h;
}

// Sums over 5 adjacent runs with Poisson errors, as the reference builds them. Additive when chunks are merged.
void SumRunGroups(const TH1D &per_run, TH1D &per_group)
{
    for (int group = 0; group < per_group.GetNbinsX(); ++group) {
        double sum = 0;
        for (int k = 0; k < Evt::kRunsPerGroup; ++k) {
            const int bin = group * Evt::kRunsPerGroup + k + 1;
            if (bin > per_run.GetNbinsX()) break;
            sum += per_run.GetBinContent(bin);
        }
        per_group.SetBinContent(group + 1, sum);
        per_group.SetBinError(group + 1, std::sqrt(std::max(0., sum)));
    }
}

void EvtCharacter(const std::string &config_path, const std::string &product,
                  const std::string &input = "", const std::string &tag = "")
{
    try {
        /////////////////////////////////////
        // Configuration, inputs and histograms
        /////////////////////////////////////
        PJ::Job job(config_path, product, input, tag, "events");
        const PJ::Config &config = job.config;
        const PJ::BinLayout bins(config);
        PJ::HistogramList booked;
        const EventHistograms h = BookEventHistograms(booked, config, bins);
        const PJ::EventCounter counter(booked, bins);
        const std::vector<double> &pt = config.pT_bins;

        for (const PJ::InputPart &part : job.parts) {
            std::unique_ptr<TFile> file = PJ::OpenRoot(part.path);
            const std::vector<PJ::Event> events = PJ::LoadPart(*file, config, job.simulation);

            for (const PJ::Event &event : events) {
                const int centrality_bin = bins.CentralityBin(event.centrality);
                counter.Read(centrality_bin);

                /////////////////////////////////////
                // Event selection
                /////////////////////////////////////
                if (!job.simulation && !PJ::PassesRunAndTrigger(config, event)) continue;
                if (job.simulation && (!job.OwnsSample(event) || !PJ::InSampleWindow(event))) continue;   // the sample stitching of PPG12
                if (std::isfinite(event.vertex_z)) h.vertex_before->Fill(event.vertex_z);
                if (centrality_bin < 0 || !PJ::HasVertex(config, event)) continue;
                // v21: PPG12's weight, sample cross section times period luminosity, SI fraction and truth-vertex
                // reweight (Job::Weight). The reference uses its own vertex reweight only.
                const double weight = job.Weight(event);
                counter.Selected(centrality_bin, weight);
                // v21: pp has no centrality. The reference uses the placeholder 50 there.
                const double centrality = config.system == "auau" ? event.centrality : 50.;

                /////////////////////////////////////
                // Event quantities
                /////////////////////////////////////
                if (std::isfinite(event.emcal_total_energy)) {
                    h.emcal->Fill(event.emcal_total_energy);
                    h.emcal_weighted->Fill(event.emcal_total_energy, weight);
                    h.emcal_vs_centrality->Fill(centrality, event.emcal_total_energy);
                    h.emcal_per_run->Fill(event.run, event.emcal_total_energy);
                }
                h.mc_weight->Fill(weight);
                h.vertex->Fill(event.vertex_z, weight);
                h.centrality->Fill(centrality, weight);
                h.centrality_fine->Fill(centrality, weight);
                // v21: events.jet_background_psi2 is NaN in pp and identically 0 in Au+Au data, so h_Psi2 stays empty.
                // if (std::isfinite(event.jet_background_psi2)) h.psi2->Fill(event.jet_background_psi2, weight);

                /////////////////////////////////////
                // Clusters with ET >= 10 GeV
                /////////////////////////////////////
                int n_clusters = 0;
                for (const PJ::Photon &photon : event.photons) {
                    const double et = photon.photon_et, eta = photon.photon_eta;
                    if (et < Evt::kClusterEtMin) continue;
                    // v21: the NPB score exists in pp only, so the cut is a switch. The reference always applies it.
                    if (config.npb_cut_on && !(photon.npb_score >= config.npb_score_cut)) continue;
                    if (eta < config.eta_bins.front() || eta > config.eta_bins.back()) continue;
                    ++n_clusters;
                    h.cluster_et->Fill(et, weight);
                    h.cluster_eta_phi->Fill(eta, photon.photon_phi, weight);
                    h.cluster_ieta_iphi->Fill(photon.shower_center_eta_index, photon.shower_center_phi_index, weight);

                    const double iso = config.cone_size == 3 ? photon.iso_r03 : photon.iso_r04;
                    h.isolation->Fill(et, iso, weight);
                    h.bdt_score->Fill(et, photon.bdt_score, weight);
                    const int pt_bin = PJ::FindBin(et, pt);
                    if (pt_bin < 0) continue;
                    h.iso_vs_centrality[pt_bin]->Fill(centrality, iso, weight);
                    h.cluster_yield_vs_centrality[pt_bin]->Fill(centrality, weight);
                    for (const PJ::Jet &jet : event.jets) {
                        if (!PJ::IsConfiguredJet(jet, config) || jet.pt < Evt::kJetPtMinDphi || std::fabs(jet.eta) > config.jet_eta_max) continue;
                        const double dphi = PJ::DeltaPhi(photon.photon_phi, jet.phi);
                        if (std::hypot(eta - jet.eta, dphi) <= Evt::kClusterJetDrMin) continue;
                        h.dphi_cluster_jets[centrality_bin][pt_bin]->Fill(std::fabs(dphi), weight);
                    }
                }
                if (n_clusters > 0) h.clusters_per_run->Fill(event.run, n_clusters);
                h.photon_multiplicity->Fill(n_clusters, weight);
                if (std::isfinite(event.total_calo_energy)) h.calo_energy->Fill(event.total_calo_energy, weight);

                /////////////////////////////////////
                // Jets, in events with such a cluster
                /////////////////////////////////////
                if (n_clusters == 0) continue;
                h.centrality_with_cluster->Fill(centrality, weight);
                for (const PJ::Jet &jet : event.jets) {
                    if (!PJ::IsConfiguredJet(jet, config)) continue;
                    h.jet_pt->Fill(jet.pt, weight);
                    if (jet.pt < Evt::kJetPtMinEtaPhi || std::fabs(jet.eta) > config.jet_eta_max) continue;
                    h.jet_eta_phi->Fill(jet.eta, jet.phi, weight);
                    // h_Psi2_vs_jetPt: not filled, see h_Psi2 above.
                }
            }
        }

        /////////////////////////////////////
        // Per-run sums
        /////////////////////////////////////
        SumRunGroups(*h.clusters_per_run, *h.clusters_per_5runs);
        // v21: h_scaler_per_run stays empty. The reference takes the maximum of a per-event scaler in each run,
        // which the trees do not store, and a maximum cannot be added up when chunks are merged. The run totals
        // are in the triggerRunInfo tree of the input parts and belong to the luminosity step.

        PJ::WriteOutput(job, booked);
    } catch (const std::exception &error) { PJ::Fail(error); }
}
