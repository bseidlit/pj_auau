#pragma once
// Direct PhotonJetTrees_v1 maker: setup, booking, event/reco/truth loops, writing.
// root -l -b -q 'histmakers/PhotonJetHistMaker.C("config.yaml","sim_signal","parts.list","tag")'
#include "PhotonJetIO.h"
#include "PhotonJetHistograms.h"
#include <chrono>
#include <iomanip>
#include <iostream>

namespace PJ {
inline void Run(const std::string &config_path, const std::string &product,
                const std::string &chunk_list = "", const std::string &tag = "",
                const std::string &audit_path = "")
{
    // ---- configuration and checked original inputs --------------------------
    const auto start = std::chrono::steady_clock::now();
    const auto config = ReadConfig(config_path);
    const auto kind = ParseProduct(product);
    const auto cuts = LoadCuts(config.yaml);
    const BinLayout bins(cuts);
    auto inputs = ReadInputs(chunk_list, config, kind, cuts);
    if (!inputs.job_tag.empty() && tag != inputs.job_tag) throw std::runtime_error("chunk tag differs from sealed job plan");
    const auto model = Configure(config, kind, inputs);
    auto provenance = Provenance(config, kind, inputs, model);
    const auto names = OutputPaths(config.yaml, kind, tag);
    CheckOutputDestinations(config, inputs, names);
    if (model.signal && OutputIdentityPath(names.main) == OutputIdentityPath(names.response))
        throw std::runtime_error("main and response destinations must differ");
    // Optional accepted-candidate evidence used by the regression driver.
    std::ofstream audit;
    if (!audit_path.empty()) {
        if (std::filesystem::exists(audit_path) || OutputIdentityPath(audit_path) == OutputIdentityPath(names.main) ||
            (model.signal && OutputIdentityPath(audit_path) == OutputIdentityPath(names.response)))
            throw std::runtime_error("response audit requires a distinct new file");
        audit.open(audit_path);
        if (!audit) throw std::runtime_error("cannot open response audit " + audit_path);
    }
    const auto configured = std::chrono::steady_clock::now();

    // ---- histogram booking -------------------------------------------------
    auto qa = BookJobQA(cuts);
    auto histograms = BookHistograms(cuts, bins, model.signal);
    std::vector<ULong64_t> selected_by_centrality(bins.nCentrality(), 0);
    ULong64_t events_processed = 0;
    size_t parts_processed = 0;
    double loading_seconds = 0, processing_seconds = 0;
    const bool use_links = model.simulation && cuts.match_source == "links";
    const auto prior = [&](double pt) {
        return model.prior ? model.prior->Eval(std::min(cuts.trw_clamp_max, std::max(cuts.trw_clamp_min, pt))) : 1.;
    };

    // ---- one original part at a time ---------------------------------------
    for (const auto &part : inputs.originals) {
        const auto loading_start = std::chrono::steady_clock::now();
        auto file = OpenRoot(part.path);
        const auto data = LoadPart(*file, !model.simulation && cuts.system == "auau", model.simulation, cuts.system == "auau");
        const auto &record = inputs.records[parts_processed];
        if (data.events.size() != record["events"].as<size_t>() || data.photons.size() != record["photons"].as<size_t>() ||
            data.truths.size() != record["truth_photons"].as<size_t>() || data.link_rows != record["source_link_rows"].as<Long64_t>() ||
            data.jet_rows != record["truth_jet_rows"].as<Long64_t>() ||
            data.links_available != record["links_available"].as<bool>() || data.jets_available != record["truth_jets_available"].as<bool>() ||
            file->GetUUID().AsString() != record["source_uuid"].as<std::string>())
            throw std::runtime_error("original input inventory changed while reading " + part.path);
        file->Close();
        const auto loaded = std::chrono::steady_clock::now();
        loading_seconds += std::chrono::duration<double>(loaded - loading_start).count();
        const auto &sample = model.samples.at(part.index);

        // ---- event selection and job-level accounting ----------------------
        for (size_t e = 0; e < data.events.size(); ++e) {
            const auto &event = data.events[e];
            const auto matches = MatchPhotons(data, e, model.simulation, use_links);
            ++events_processed;
            qa.cutflow->Fill(.5);
            const auto selected = SelectEvent(cuts, bins, model.simulation, sample, model.vertex, data, e);
            if (!selected.pass) continue;
            const int icent = selected.centrality_bin;
            const double event_weight = selected.weight;
            ++selected_by_centrality[icent];
            qa.cutflow->Fill(1.5);
            qa.vertex->Fill(event.vertex_z);
            qa.weight->Fill(std::log10(std::max(event_weight, 1e-12)));
            qa.centrality->Fill(event.centrality);
            qa.run->Fill(event.run);

            // ---- reco photons: the decision and its fills stay together -----
            std::vector<Candidate> candidates;
            candidates.reserve(data.reco_rows[e].size());
            for (const auto row : data.reco_rows[e]) {
                const auto &photon = data.photons[row];
                qa.cutflow->Fill(2.5);
                auto candidate = SelectPhoton(photon, cuts, bins, model.simulation, sample, event_weight, matches[photon.original_index]);
                if (!candidate.accepted) continue;
                auto &h = histograms[icent][candidate.eta_bin];
                const auto &d = candidate.decision;
                const double photon_weight = candidate.weight;
                qa.cutflow->Fill(3.5);
                if (cuts.system == "auau")
                    qa.recipe->Fill(std::fabs(photon.iso_r04_threshold - (.081 + .019 * photon.photon_et)) < 1e-6 ? 1.5 : .5);
                if (!model.simulation && cuts.threshold_source == "formula" && cuts.flag_check != "off") {
                    qa.flags->Fill(.5);
                    if (d.flag_mismatch_tight) qa.flags->Fill(1.5);
                    if (d.flag_mismatch_nontight) qa.flags->Fill(2.5);
                    if (d.flag_mismatch_iso) qa.flags->Fill(3.5);
                }
                h.reco.all->Fill(d.et, photon_weight);
                if (d.common) {
                    qa.cutflow->Fill(4.5);
                    h.reco.common->Fill(d.et, photon_weight);
                    h.qa.isolation->Fill(d.et, d.iso, photon_weight);
                    h.qa.score->Fill(d.et, photon.bdt_score, photon_weight);
                }
                if (d.tight) { qa.cutflow->Fill(5.5); h.reco.tight->Fill(d.et, photon_weight); }
                if (d.region >= 0) { qa.cutflow->Fill(6.5 + d.region); h.reco.abcd[d.region]->Fill(d.et, photon_weight); }
                candidates.push_back(candidate);
            }
            if (!model.simulation) continue;

            // ---- truth, efficiencies and response: each truth counted once --
            for (const auto row : data.truth_rows[e]) {
                const auto &truth = data.truths[row];
                const int ieta = bins.EtaBin(truth.truth_photon_eta);
                if (ieta < 0 || !IsFiducialTruth(truth, cuts)) continue;
                auto &h = histograms[icent][ieta];
                const double pt = truth.truth_photon_pt;
                qa.cutflow->Fill(10.5);
                h.truth.spectrum->Fill(pt, event_weight);
                h.truth.novtx->Fill(pt, event_weight);
                h.truth.vertexcut->Fill(pt, event_weight);
                h.truth.mbd->Fill(pt, event_weight);
                h.truth.north->Fill(pt, event_weight);
                h.truth.south->Fill(pt, event_weight);
                bool reconstructed = false, isolated = false, tight_isolated = false;
                for (auto &candidate : candidates) {
                    if (candidate.truth_index != truth.original_index || candidate.eta_bin != ieta) continue;
                    candidate.matched_to_fiducial = true;
                    reconstructed = true;
                    const auto &d = candidate.decision;
                    isolated |= d.iso_pass;
                    tight_isolated |= d.tight && d.iso_pass;
                    if (pt > cuts.pT_bins_truth.front() && pt < cuts.pT_bins_truth.back() &&
                        d.et > cuts.pT_bins.front() && d.et < cuts.pT_bins.back()) {
                        h.reco.signal_all->Fill(d.et, event_weight);
                        if (d.tight) h.reco.signal_tight->Fill(d.et, event_weight);
                        if (d.region >= 0) h.reco.signal[d.region]->Fill(d.et, event_weight);
                    }
                    if (model.signal) {
                        const auto point = MakeResponse(cuts, candidate, truth, event, inputs.ids.at(part.index), event_weight, prior);
                        if (!point.accepted) continue;
                        const double response_weight = point.weight;
                        h.response.reco->Fill(point.reco, response_weight);
                        h.response.truth->Fill(point.truth, response_weight);
                        h.response.matrix->Fill(point.reco, point.truth, response_weight);
                        h.response.unfold->Fill(point.reco, point.truth, response_weight);
                        if (audit.is_open())
                            audit << inputs.ids.at(part.index) << '/' << event.source_file_index << '/' << event.event_id_hi << '/'
                                  << event.event_id_lo << '/' << candidate.original_index << ' ' << point.seed << ' ' << std::hexfloat
                                  << point.gaussian << ' ' << point.reco << ' ' << point.truth << ' ' << point.weight << '\n';
                    }
                }
                if (reconstructed) qa.cutflow->Fill(11.5);
                h.truth.reco->FillWeighted(reconstructed, event_weight, pt);
                h.truth.all->FillWeighted(tight_isolated, event_weight, pt);
                if (reconstructed) {
                    h.truth.iso->FillWeighted(isolated, event_weight, pt);
                    if (isolated) h.truth.id->FillWeighted(tight_isolated, event_weight, pt);
                }
            }
            for (const auto &candidate : candidates) {
                const auto &d = candidate.decision;
                if (!candidate.matched_to_fiducial && d.region >= 0)
                    histograms[icent][candidate.eta_bin].reco.unmatched[d.region]->Fill(d.et, candidate.weight);
            }
        }
        processing_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - loaded).count();
        ++parts_processed;
    }

    // ---- validate, write and publish one complete output generation ---------
    if (events_processed != inputs.events) throw std::runtime_error("incomplete event processing");
    if (!model.simulation && cuts.threshold_source == "formula" && cuts.flag_check == "strict" &&
        qa.flags->GetBinContent(2) + qa.flags->GetBinContent(3) + qa.flags->GetBinContent(4) > 0)
        throw std::runtime_error("stored flags disagree with configured thresholds; see h_pj_flagcheck");
    if (audit.is_open()) { audit.flush(); if (!audit) throw std::runtime_error("response audit write failed"); audit.close(); }
    provenance["events_processed"] = events_processed; provenance["parts_processed"] = parts_processed;
    provenance["workers"] = 1; provenance["event_loops"] = 1;
    provenance["selected_events_by_centrality"] = selected_by_centrality;
    provenance["preflight_seconds"] = std::chrono::duration<double>(configured - start).count();
    provenance["loading_seconds"] = loading_seconds; provenance["processing_seconds"] = processing_seconds;
    provenance["analysis_seconds"] = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    RecheckInputs(inputs, provenance);
    RootOutput main(names.main);
    std::unique_ptr<RootOutput> response;
    if (model.signal) response = std::make_unique<RootOutput>(names.response);
    WriteHistograms(*main.get(), response ? response->get() : nullptr, qa, histograms);
    WriteRunMetadata(main.get(), config, provenance, inputs.originals, parts_processed);
    if (response) WriteRunMetadata(response->get(), config, provenance, inputs.originals, parts_processed);
    Publish(main, response.get(), names, provenance);
    std::cout << "[PJ] parts=" << parts_processed << " events=" << events_processed << " cells=" << bins.nCentrality() << 'x'
              << bins.nEta() << " wrote " << names.main << std::endl;
}
} // namespace PJ

void PhotonJetHistMaker(const std::string &config = "", const std::string &product = "data",
                       const std::string &chunk_list = "", const std::string &chunk_tag = "")
{
    try { PJ::Run(config, product, chunk_list, chunk_tag); }
    catch (const std::exception &error) { std::cerr << "[PJ] ERROR: " << error.what() << std::endl; gSystem->Exit(1); }
}
