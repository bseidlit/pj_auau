// Photon spectra, ABCD regions, truth-efficiency counts, the matched response and the photon-jet
// observable x_Jgamma = jet pT / photon ET of pj_auau/histMakers/RecoEffCalculator_TTreeReader.C.
// x_Jgamma is an observable of the leading photon of the event and its back-to-back jets.
#include "PhotonJetIO.h"
#include "PhotonJetHistograms.h"
#include <TString.h>

/////////////////////////////////////
// Histograms of one centrality x eta cell. Names follow the reference:
// <base>_<c> for the spectra, <base>_cent<c> for everything with jets.
/////////////////////////////////////
struct SpectrumHistograms
{
    // Reconstructed photons vs calibrated ET. Data and simulation.
    TH1D *all = nullptr, *common = nullptr, *tight = nullptr;
    TH1D *abcd[PJ::kNRegions] = {};                                     // A, B, C, D
    TH2D *xj_all = nullptr, *xj_common = nullptr, *xj_tight = nullptr;  // x_Jgamma vs photon ET, one entry per back-to-back jet
    TH2D *xj_abcd[PJ::kNRegions] = {};
    std::vector<TH1D *> iso_tight, iso_nontight;                        // isolation ET, one per photon ET bin
    std::vector<TH1D *> dphi_tight;                                     // |dphi|(tight photon, jet), one per photon ET bin
    // Simulation: every ABCD region split by truth, so that abcd = signal + background.
    TH1D *signal[PJ::kNRegions] = {};                                   // matched to a fiducial truth photon, both inside the binned ranges
    TH1D *background[PJ::kNRegions] = {};                               // all other candidates
    TH1D *unmatched[PJ::kNRegions] = {};                                // no fiducial truth photon at all
    TH1D *all_signal = nullptr, *tight_signal = nullptr;
    TH2D *all_signal_max_b2b = nullptr;                                 // photon ET vs the largest back-to-back jet pT
    TH1D *fake = nullptr;                                               // tight isolated, matched to a prompt truth photon that is not isolated
    TH2D *xj_signal[PJ::kNRegions] = {}, *xj_all_signal = nullptr, *xj_tight_signal = nullptr;
    TH2D *xj_matched_signal[PJ::kNRegions] = {}, *xj_matched_all_signal = nullptr, *xj_matched_tight_signal = nullptr;  // truth-matched jets only
    TH2D *xj_truth_signal[PJ::kNRegions] = {}, *xj_truth_all_signal = nullptr, *xj_truth_tight_signal = nullptr;        // truth jets around the truth photon
    TH2D *xj_background[PJ::kNRegions] = {};
    TH2D *jet_response = nullptr;                                       // reco / truth jet pT vs truth jet pT
    // Simulation: fiducial truth photons vs truth pT.
    TH1D *truth_before = nullptr;                                       // every event: denominator of the vertex efficiency
    TH1D *truth_vertex = nullptr;                                       // events with a reconstructed vertex inside the cut
    TH1D *truth = nullptr;                                              // in selected events: denominator of the photon efficiencies
    TH1D *reconstructed = nullptr, *isolated = nullptr, *identified = nullptr;
    // Signal simulation: the matched response in ET and, for the leading photon, in (x_Jgamma, ET).
    TH1D *response_reco = nullptr, *response_truth = nullptr;
    TH2D *response = nullptr;
    TH2D *xj_response = nullptr;                                        // flattened: truth global bin vs reco global bin, pT slow and x_Jgamma fast
    TH2D *xj_response_reco = nullptr, *xj_response_truth = nullptr;     // every reco entry (fakes included), every truth entry (misses included)
};
using SpectrumGrid = std::vector<std::vector<SpectrumHistograms>>;   // [centrality][eta]

SpectrumGrid BookSpectrumHistograms(PJ::HistogramList &booked, const PJ::Config &config, const PJ::BinLayout &bins,
                                    bool simulation, bool signal)
{
    using PJ::UniformBins;
    const char *regions[PJ::kNRegions] = {"tight_iso", "tight_noniso", "nontight_iso", "nontight_noniso"};
    const std::vector<double> &pt = config.pT_bins, &pt_truth = config.pT_bins_truth;
    const std::vector<double> &xj = config.xjgamma_bins, &xj_truth = config.xjgamma_bins_truth;
    const int n_reco = (pt.size() - 1) * (xj.size() - 1), n_truth = (pt_truth.size() - 1) * (xj_truth.size() - 1);
    const double pi = 3.141592653589793;
    SpectrumGrid grid(bins.nCentrality());
    for (int ic = 0; ic < bins.nCentrality(); ++ic) {
        for (int ie = 0; ie < bins.nEta(); ++ie) {
            const PJ::BinInfo cell = bins.cell(ic, ie);
            const std::string reco_title = cell.Title() + ";Photon E_{T} [GeV];Weighted photons";
            const std::string truth_title = cell.Title() + ";Truth photon p_{T} [GeV];Weighted photons";
            const std::string xj_title = cell.Title() + ";x_{J#gamma};Photon E_{T} [GeV]";
            const std::string xj_truth_title = cell.Title() + ";Truth x_{J#gamma};Photon E_{T} [GeV]";
            SpectrumHistograms h;
            h.all = booked.Book1D(cell.Name("h_all_cluster"), reco_title, pt);
            h.common = booked.Book1D(cell.Name("h_common_cluster"), reco_title, pt);
            h.tight = booked.Book1D(cell.Name("h_tight_cluster"), reco_title, pt);
            h.xj_all = booked.Book2D(cell.CentName("h_all_xjgamma"), xj_title, xj, pt);
            h.xj_common = booked.Book2D(cell.CentName("h_common_xjgamma"), xj_title, xj, pt);
            h.xj_tight = booked.Book2D(cell.CentName("h_tight_xjgamma"), xj_title, xj, pt);
            for (int r = 0; r < PJ::kNRegions; ++r) {
                const std::string region = std::string("h_") + regions[r];
                h.abcd[r] = booked.Book1D(cell.Name(region + "_cluster"), reco_title, pt);
                h.xj_abcd[r] = booked.Book2D(cell.CentName(region + "_xjgamma"), xj_title, xj, pt);
                if (!simulation) continue;
                h.signal[r] = booked.Book1D(cell.Name(region + "_cluster_signal"), reco_title, pt);
                h.background[r] = booked.Book1D(cell.Name(region + "_cluster_background"), reco_title, pt);
                h.unmatched[r] = booked.Book1D(cell.Name(region + "_cluster_notmatch"), reco_title, pt);
                h.xj_signal[r] = booked.Book2D(cell.CentName(region + "_xjgamma_signal"), xj_title, xj, pt);
                h.xj_matched_signal[r] = booked.Book2D(cell.CentName(region + "_truthmatchreco_xjgamma_signal"), xj_title, xj, pt);
                h.xj_truth_signal[r] = booked.Book2D(cell.CentName(region + "_truthjet_xjgamma_signal"), xj_truth_title, xj_truth, pt);
                h.xj_background[r] = booked.Book2D(cell.CentName(region + "_xjgamma_background"), xj_title, xj, pt);
            }
            for (size_t i = 0; i + 1 < pt.size(); ++i) {
                const std::string bin = "_pt" + std::to_string(i), range = Form(", %g < E_{T} < %g GeV", pt[i], pt[i + 1]);
                h.iso_tight.push_back(booked.Book1D(cell.CentName("h_tight_isoET") + bin, cell.Title() + range + ";Isolation E_{T} [GeV];Weighted photons", UniformBins(1000, -50, 50)));
                h.iso_nontight.push_back(booked.Book1D(cell.CentName("h_nontight_isoET") + bin, cell.Title() + range + ";Isolation E_{T} [GeV];Weighted photons", UniformBins(1000, -50, 50)));
                h.dphi_tight.push_back(booked.Book1D(cell.CentName("h_dphi_clusterJets_tight") + bin, cell.Title() + range + ";|#Delta#phi|(tight photon, jet);Weighted pairs", UniformBins(64, 0, pi)));
            }
            if (simulation) {
                h.all_signal = booked.Book1D(cell.Name("h_all_cluster_signal"), reco_title, pt);
                h.tight_signal = booked.Book1D(cell.Name("h_tight_cluster_signal"), reco_title, pt);
                h.all_signal_max_b2b = booked.Book2D(cell.Name("h_all_cluster_Et_max_b2bjet"), cell.Title() + ";Photon E_{T} [GeV];Largest back-to-back jet p_{T} [GeV]", pt, UniformBins(100, 0, 100));
                h.fake = booked.Book1D(cell.Name("h_pT_reco_fake"), reco_title, pt);
                h.xj_all_signal = booked.Book2D(cell.CentName("h_all_xjgamma_signal"), xj_title, xj, pt);
                h.xj_tight_signal = booked.Book2D(cell.CentName("h_tight_xjgamma_signal"), xj_title, xj, pt);
                h.xj_matched_all_signal = booked.Book2D(cell.CentName("h_all_truthmatchreco_xjgamma_signal"), xj_title, xj, pt);
                h.xj_matched_tight_signal = booked.Book2D(cell.CentName("h_tight_truthmatchreco_xjgamma_signal"), xj_title, xj, pt);
                h.xj_truth_all_signal = booked.Book2D(cell.CentName("h_all_truthjet_xjgamma_signal"), xj_truth_title, xj_truth, pt);
                h.xj_truth_tight_signal = booked.Book2D(cell.CentName("h_tight_truthjet_xjgamma_signal"), xj_truth_title, xj_truth, pt);
                h.jet_response = booked.Book2D(cell.CentName("h_jet_pT_response"), cell.Title() + ";Truth jet p_{T} [GeV];Reco jet p_{T} / truth jet p_{T}", UniformBins(45, 5, 50), UniformBins(100, 0, 2));
                h.truth_before = booked.Book1D(cell.Name("h_truth_pT"), truth_title, pt_truth);
                h.truth_vertex = booked.Book1D(cell.Name("h_truth_pT_vertexcut"), truth_title, pt_truth);
                h.truth = booked.Book1D(cell.Name("h_truth_pT_selected"), truth_title, pt_truth);
                h.reconstructed = booked.Book1D(cell.Name("h_truth_pT_reconstructed"), truth_title, pt_truth);
                h.isolated = booked.Book1D(cell.Name("h_truth_pT_isolated"), truth_title, pt_truth);
                h.identified = booked.Book1D(cell.Name("h_truth_pT_identified"), truth_title, pt_truth);
            }
            if (signal) {
                h.response_reco = booked.Book1D(cell.Name("h_pT_reco_response"), reco_title, pt);
                h.response_truth = booked.Book1D(cell.Name("h_pT_truth_response"), truth_title, pt_truth);
                h.response = booked.Book2D(cell.Name("h_response_full"), cell.Title() + ";Reco E_{T} [GeV];Truth p_{T} [GeV]", pt, pt_truth);
                h.xj_response = booked.Book2D(cell.CentName("h_response_xjgamma_global"), cell.Title() + ";Truth global bin (p_{T}, x_{J#gamma});Reco global bin (E_{T}, x_{J#gamma})", UniformBins(n_truth, 0.5, n_truth + 0.5), UniformBins(n_reco, 0.5, n_reco + 0.5));
                h.xj_response_reco = booked.Book2D(cell.CentName("h_xjgamma_reco_response"), cell.Title() + ";Reco x_{J#gamma};Reco E_{T} [GeV]", xj, pt);
                h.xj_response_truth = booked.Book2D(cell.CentName("h_xjgamma_truth_response"), cell.Title() + ";Truth x_{J#gamma};Truth p_{T} [GeV]", xj_truth, pt_truth);
            }
            grid[ic].push_back(h);
        }
    }
    return grid;
}

// One entry per back-to-back jet at x_Jgamma = jet pT / photon ET, or one underflow entry when there is none,
// so that the ET projection still counts every photon.
void FillXJ(TH2D *h, const std::vector<double> &jets, double et, double weight)
{
    if (jets.empty()) h->Fill(h->GetXaxis()->GetBinLowEdge(1) - 1e-4, et, weight);
    for (double pt : jets) h->Fill(pt / et, et, weight);
}

void RecoEffCalculator_TTreeReader(const std::string &config_path, const std::string &product,
                                  const std::string &input = "", const std::string &tag = "")
{
    try {
        /////////////////////////////////////
        // Configuration, inputs and histograms
        /////////////////////////////////////
        PJ::Job job(config_path, product, input, tag, "spectrum");
        const PJ::Config &config = job.config;
        const PJ::BinLayout bins(config);
        PJ::HistogramList booked;
        SpectrumGrid histograms = BookSpectrumHistograms(booked, config, bins, job.simulation, job.signal);
        TH1D *h_vertexz = booked.Book1D("h_vertexz", "Vertex Z;Vertex Z [cm];Weighted events", PJ::UniformBins(200, -100, 100));
        const PJ::EventCounter counter(booked, bins);
        const int n_xj = config.xjgamma_bins.size() - 1, n_xj_truth = config.xjgamma_bins_truth.size() - 1;

        for (const PJ::InputPart &part : job.parts) {
            /////////////////////////////////////
            // Initialize one part: read its events, then give every photon its random number
            /////////////////////////////////////
            std::unique_ptr<TFile> file = PJ::OpenRoot(part.path);
            std::vector<PJ::Event> events = PJ::LoadPart(*file, config, job.simulation);
            // The random number is used by the response smearing only, so only signal simulation needs it.
            if (job.signal) PJ::InitRandom(events, config.random_seed, part.number);

            for (PJ::Event &event : events) {
                const int centrality = bins.CentralityBin(event.centrality);
                counter.Read(centrality);
                if (!job.simulation && !PJ::PassesRunAndTrigger(config, event)) continue;
                if (job.simulation && (!job.OwnsSample(event) || !PJ::InSampleWindow(event))) continue;   // the sample stitching of PPG12
                const double weight = job.Weight(event);

                /////////////////////////////////////
                // Simulation: truth photons before the reconstructed-event selection, for the vertex efficiency
                /////////////////////////////////////
                if (job.simulation && centrality >= 0) {
                    for (const PJ::TruthPhoton &truth : event.truths) {
                        const int eta = bins.EtaBin(truth.truth_photon_eta);
                        if (eta < 0 || !PJ::IsFiducialTruth(truth, config)) continue;
                        SpectrumHistograms &h = histograms[centrality][eta];
                        h.truth_before->Fill(truth.truth_photon_pt, weight);
                        if (PJ::HasVertex(config, event)) h.truth_vertex->Fill(truth.truth_photon_pt, weight);
                    }
                }

                /////////////////////////////////////
                // Event selection
                /////////////////////////////////////
                if (centrality < 0 || !PJ::HasVertex(config, event)) continue;
                counter.Selected(centrality, weight);
                h_vertexz->Fill(event.vertex_z, weight);
                std::vector<SpectrumHistograms> &cell = histograms[centrality];

                /////////////////////////////////////
                // Simulation: jet matching and the jet pT response. Every truth jet in the acceptance takes
                // the nearest reconstructed jet in the acceptance within dR < 0.2.
                /////////////////////////////////////
                if (job.simulation) {
                    PJ::MatchJets(config, event);
                    for (const PJ::Jet &truth : event.truth_jets) {
                        if (!PJ::IsConfiguredJet(truth, config) || truth.near_photon || truth.pt <= 0 || std::fabs(truth.eta) >= config.jet_eta_max) continue;
                        const PJ::Jet *nearest = nullptr;
                        double best = 0.2;
                        for (const PJ::Jet &jet : event.jets) {
                            if (!PJ::IsConfiguredJet(jet, config) || std::fabs(jet.eta) >= config.jet_eta_max) continue;
                            const double dr = PJ::DeltaR(jet.eta, jet.phi, truth.eta, truth.phi);
                            if (dr < best) { best = dr; nearest = &jet; }
                        }
                        // Jets have no eta cell: they go with the first eta bin.
                        if (nearest) cell[0].jet_response->Fill(truth.pt, nearest->pt / truth.pt, weight);
                    }
                }

                /////////////////////////////////////
                // Reconstructed photons: acceptance, classification and the back-to-back jets of each
                /////////////////////////////////////
                std::vector<PJ::Candidate> candidates;
                for (const PJ::Photon &photon : event.photons) {
                    PJ::Candidate k = PJ::SelectPhoton(photon, config, bins, job.simulation, weight);
                    if (!k.accepted || k.et > event.cluster_et_upper) continue;   // jet samples: photons above belong to the next sample
                    k.b2b_jets = PJ::BackToBackJets(config, event.jets, k.phi);
                    if (job.simulation) k.b2b_matched_jets = PJ::BackToBackJets(config, event.jets, k.phi, true);
                    candidates.push_back(k);
                }
                if (!candidates.empty())
                    std::max_element(candidates.begin(), candidates.end(), [](const PJ::Candidate &a, const PJ::Candidate &b) { return a.et < b.et; })->leading = true;

                /////////////////////////////////////
                // Spectra, x_Jgamma of the leading photon, isolation and photon-jet |dphi| per ET bin
                /////////////////////////////////////
                for (const PJ::Candidate &k : candidates) {
                    SpectrumHistograms &h = cell[k.eta_bin];
                    auto fill = [&](TH1D *spectrum, TH2D *xj) {
                        spectrum->Fill(k.et, k.weight);
                        if (k.leading) FillXJ(xj, k.b2b_jets, k.et, k.weight);
                    };
                    fill(h.all, h.xj_all);
                    if (k.common) fill(h.common, h.xj_common);
                    if (k.tight) fill(h.tight, h.xj_tight);
                    if (k.region != PJ::kNoRegion) fill(h.abcd[k.region], h.xj_abcd[k.region]);
                    const int pt_bin = PJ::FindBin(k.et, config.pT_bins);
                    if (pt_bin < 0) continue;
                    if (k.tight) h.iso_tight[pt_bin]->Fill(k.iso, k.weight);
                    if (k.nontight) h.iso_nontight[pt_bin]->Fill(k.iso, k.weight);
                    if (!k.tight) continue;
                    // Every jet of the collection in the acceptance and above the x_Jgamma threshold, except the photon itself.
                    for (const PJ::Jet &jet : event.jets)
                        if (PJ::IsConfiguredJet(jet, config) && std::fabs(jet.eta) < config.jet_eta_max && jet.pt > config.b2bjet_pT_min &&
                            PJ::DeltaR(k.eta, k.phi, jet.eta, jet.phi) > 0.2)
                            h.dphi_tight[pt_bin]->Fill(std::fabs(PJ::DeltaPhi(k.phi, jet.phi)), k.weight);
                }
                if (!job.simulation) continue;

                /////////////////////////////////////
                // Simulation: efficiencies, signal in the ABCD regions and the responses
                /////////////////////////////////////
                // These efficiencies are conditional on the selected reconstructed event.
                for (size_t t = 0; t < event.truths.size(); ++t) {
                    const PJ::TruthPhoton &truth = event.truths[t];
                    const int eta = bins.EtaBin(truth.truth_photon_eta);
                    if (eta < 0 || !PJ::IsFiducialTruth(truth, config)) continue;
                    SpectrumHistograms &h = cell[eta];
                    const double pt = truth.truth_photon_pt;
                    h.truth->Fill(pt, weight);
                    // The truth jets back-to-back with this truth photon, for the truth x_Jgamma.
                    const std::vector<double> truth_b2b = PJ::BackToBackJets(config, event.truth_jets, truth.truth_photon_phi);

                    bool reconstructed = false, isolated = false, identified = false;
                    // The candidates matched to this truth photon: truth_index is a position in event.truths.
                    for (PJ::Candidate &k : candidates) {
                        if (k.truth_index != static_cast<int>(t) || k.eta_bin != eta) continue;
                        k.matched_to_fiducial = true;
                        reconstructed = true;
                        isolated |= k.iso_pass;
                        identified |= k.tight && k.iso_pass;
                        // Signal: both inside the binned ranges. Reco x_Jgamma with all and with truth-matched jets, and truth x_Jgamma.
                        if (PJ::InRange(pt, config.pT_bins_truth) && PJ::InRange(k.et, config.pT_bins)) {
                            k.signal = true;
                            auto fill = [&](TH1D *spectrum, TH2D *xj, TH2D *xj_matched, TH2D *xj_truth) {
                                spectrum->Fill(k.et, weight);
                                if (!k.leading) return;
                                FillXJ(xj, k.b2b_jets, k.et, weight);
                                FillXJ(xj_matched, k.b2b_matched_jets, k.et, weight);
                                FillXJ(xj_truth, truth_b2b, pt, weight);
                            };
                            fill(h.all_signal, h.xj_all_signal, h.xj_matched_all_signal, h.xj_truth_all_signal);
                            if (!k.b2b_jets.empty()) h.all_signal_max_b2b->Fill(k.et, k.b2b_jets.back(), weight);
                            if (k.tight) fill(h.tight_signal, h.xj_tight_signal, h.xj_matched_tight_signal, h.xj_truth_tight_signal);
                            if (k.region != PJ::kNoRegion) fill(h.signal[k.region], h.xj_signal[k.region], h.xj_matched_signal[k.region], h.xj_truth_signal[k.region]);
                        }
                        if (!job.signal) continue;

                        /////////////////////////////////////
                        // Signal simulation: the response in ET, and in (x_Jgamma, ET) for the leading photon. The reco
                        // side is the smeared ET, in range after smearing as in PPG12. Reco and truth x_Jgamma values are
                        // paired largest with largest, the rest are fakes and misses.
                        /////////////////////////////////////
                        const PJ::ResponsePoint point = PJ::MakeResponse(config, k, truth, job.Prior(pt), weight);
                        if (!point.accepted) continue;
                        h.response_reco->Fill(point.reco, point.weight);
                        h.response_truth->Fill(point.truth, point.weight);
                        h.response->Fill(point.reco, point.truth, point.weight);
                        if (!k.leading) continue;
                        std::vector<double> reco_xj, truth_xj;   // largest first
                        for (auto jet = k.b2b_jets.rbegin(); jet != k.b2b_jets.rend(); ++jet) reco_xj.push_back(*jet / point.reco);
                        for (auto jet = truth_b2b.rbegin(); jet != truth_b2b.rend(); ++jet) truth_xj.push_back(*jet / pt);
                        const int reco_pt_bin = PJ::FindBin(point.reco, config.pT_bins), truth_pt_bin = PJ::FindBin(pt, config.pT_bins_truth);
                        const size_t pairs = std::min(reco_xj.size(), truth_xj.size());
                        for (size_t i = 0; i < pairs; ++i) {
                            const int reco_bin = PJ::FindBin(reco_xj[i], config.xjgamma_bins), truth_bin = PJ::FindBin(truth_xj[i], config.xjgamma_bins_truth);
                            if (reco_bin < 0 || truth_bin < 0) continue;
                            h.xj_response->Fill(truth_pt_bin * n_xj_truth + truth_bin + 1, reco_pt_bin * n_xj + reco_bin + 1, point.weight);
                            h.xj_response_reco->Fill(reco_xj[i], point.reco, point.weight);
                            h.xj_response_truth->Fill(truth_xj[i], pt, point.weight);
                        }
                        for (size_t i = pairs; i < reco_xj.size(); ++i) h.xj_response_reco->Fill(reco_xj[i], point.reco, point.weight);   // fakes
                        for (size_t i = pairs; i < truth_xj.size(); ++i) h.xj_response_truth->Fill(truth_xj[i], pt, point.weight);      // misses
                    }
                    if (reconstructed) h.reconstructed->Fill(pt, weight);
                    if (isolated) h.isolated->Fill(pt, weight);
                    if (identified) h.identified->Fill(pt, weight);
                }

                /////////////////////////////////////
                // Simulation: candidates that are not signal, and fakes
                /////////////////////////////////////
                for (const PJ::Candidate &k : candidates) {
                    SpectrumHistograms &h = cell[k.eta_bin];
                    // A prompt truth photon that is not fiducial (outside the truth isolation), reconstructed as tight and isolated.
                    if (k.tight && k.iso_pass && k.truth_index >= 0) {
                        const PJ::TruthPhoton &truth = event.truths[k.truth_index];
                        if ((truth.prompt_class == 1 || truth.prompt_class == 2) && !PJ::IsFiducialTruth(truth, config)) h.fake->Fill(k.et, k.weight);
                    }
                    if (k.region == PJ::kNoRegion) continue;
                    if (!k.matched_to_fiducial) h.unmatched[k.region]->Fill(k.et, k.weight);
                    if (k.signal) continue;
                    h.background[k.region]->Fill(k.et, k.weight);
                    if (k.leading) FillXJ(h.xj_background[k.region], k.b2b_jets, k.et, k.weight);
                }
            }
        }
        PJ::WriteOutput(job, booked);
    } catch (const std::exception &error) { PJ::Fail(error); }
}
