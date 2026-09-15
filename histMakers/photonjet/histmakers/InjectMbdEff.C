#pragma once
// Import the external MBD/vertex correction after merging signal MC.
// Validate the complete set before updating a temporary copy of the target.
// Empty external_mbd_eff_file explicitly keeps the maker's unit placeholders.
#include "PhotonJetIO.h"
#include <TNamed.h>
#include <iostream>
#include <memory>
#include <vector>

namespace PJ
{
inline void InjectCorrection(const std::string &configname, const std::string &target_override = "")
{
    const PJ::ConfigFile config = PJ::ReadConfig(configname);
    const auto &cfg = config.yaml;
    RequireSingleCellLayout(BinLayout(LoadCuts(cfg)).Metadata(), "MBD injection");
    const std::string source = PJ::Y<std::string>(cfg["photonjet"], "external_mbd_eff_file", "");
    if (source.empty())
    {
        std::cout << "[InjectMbdEff] no external source configured; target unchanged (maker outputs use unit placeholders)" << std::endl;
        return;
    }
    const std::string target = target_override.empty() ? PJ::MainOutputPath(cfg, PJ::Product::Signal) : target_override;
    const auto source_digest = PJ::FileDigest(source);
    auto fs = PJ::OpenRoot(source);
    auto original = PJ::OpenRoot(target); // READ must succeed; never create a missing merged target.
    for (auto *file : {fs.get(), original.get()})
        if (auto *saved = dynamic_cast<TObjString *>(file->Get("photonjet_provenance"))) {
            const auto info = YAML::Load(saved->GetString().Data());
            if (info["layout"]) RequireSingleCellLayout(info["layout"], "MBD injection");
        }
    if (std::filesystem::equivalent(source, target)) throw std::runtime_error("MBD source and target must differ");
    const auto edges = PJ::YV<double>(cfg["analysis"], "pT_bins_truth");
    PJ::ValidateEdges(edges, "pT_bins_truth");
    std::vector<std::unique_ptr<TH1>> replacements;
    for (const auto &name : PJ::MbdHistogramNames())
    {
        auto *source_hist = PJ::RequireHistogram(fs.get(), name, edges);
        PJ::RequireHistogram(original.get(), name, edges);
        replacements.emplace_back(static_cast<TH1 *>(source_hist->Clone(name.c_str())));
        replacements.back()->SetDirectory(nullptr);
    }
    const TH1 &denominator = *replacements.front();
    for (size_t i = 1; i < replacements.size(); ++i)
        for (int bin = 0; bin < denominator.GetNcells(); ++bin)
            if (replacements[i]->GetBinContent(bin) > denominator.GetBinContent(bin) +
                1e-12 * std::max(1.0, denominator.GetBinContent(bin)))
                throw std::runtime_error("MBD numerator exceeds denominator: " + std::string(replacements[i]->GetName()));

    // Slimtree sources normally archive their config. Check the selections that
    // define this correction, independently of BDT/reco selections and luminosity.
    bool source_selection_checked = false;
    if (auto *saved = dynamic_cast<TObjString *>(fs->Get("config")))
    {
        const auto analysis = YAML::Load(saved->GetString().Data())["analysis"];
        bool complete = true;
        for (const char *key : {"vertex_cut", "truth_iso_max", "vertex_cut_truth"})
        {
            if (!analysis[key] || !cfg["analysis"][key]) { complete = false; continue; }
            if (analysis[key].as<double>() != cfg["analysis"][key].as<double>())
                throw std::runtime_error(std::string("MBD source has a different ") + key);
        }
        if (PJ::YV<double>(analysis, "eta_bins") != PJ::YV<double>(cfg["analysis"], "eta_bins"))
            throw std::runtime_error("MBD source has different eta bins");
        source_selection_checked = complete;
    }
    YAML::Node provenance;
    provenance["source"] = source;
    if (PJ::FileDigest(source) != source_digest) throw std::runtime_error("MBD source changed during import");
    provenance["source_md5"] = source_digest;
    provenance["target_before_md5"] = PJ::FileDigest(target);
    provenance["source_selection_checked"] = source_selection_checked;
    provenance["histograms"] = PJ::MbdHistogramNames();
    original.reset();

    PJ::RootOutput output(target, true);
    for (auto &hist : replacements)
    {
        output.get()->cd();
        output.get()->Delete((std::string(hist->GetName()) + ";*").c_str());
        hist->SetDirectory(output.get());
        if (hist->Write(hist->GetName(), TObject::kOverwrite) <= 0)
            throw std::runtime_error("cannot write MBD correction");
        hist.release(); // TFile owns it until Close().
    }
    output.get()->cd();
    TNamed source_name("mbd_eff_source", source.c_str());
    if (source_name.Write("mbd_eff_source", TObject::kOverwrite) <= 0)
        throw std::runtime_error("cannot write MBD source provenance");
    PJ::WriteText(output.get(), "mbd_eff_provenance", YAML::Dump(provenance));
    output.Close();
    output.Commit();
    std::cout << "[InjectMbdEff] replaced all " << replacements.size() << " correction histograms in " << target << std::endl;
}
} // namespace PJ

void InjectMbdEff(const std::string &configname = "", const std::string &target_override = "")
{
    try { PJ::InjectCorrection(configname, target_override); }
    catch (const std::exception &error)
    {
        std::cerr << "[InjectMbdEff] ERROR: " << error.what() << std::endl;
        gSystem->Exit(1);
    }
}
