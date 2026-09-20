// Small job setup and ROOT output handling.
#pragma once
#include "PhotonJetHistograms.h"
#include "PhotonJetReader.h"
#include "PhotonJetSelection.h"
#include "TruthVertexReweightLoader.h"
#include <TF1.h>
#include <iostream>
#include <TFile.h>
#include <TObjString.h>
#include <TSystem.h>
#include <filesystem>
#include <fstream>
#include <memory>
#include <set>
#include <sstream>
#include <unistd.h>

namespace PJ {
inline std::string ReadText(const std::string &path)
{
    std::ifstream file(path);
    if (!file) throw std::runtime_error("cannot read " + path);
    std::ostringstream text; text << file.rdbuf();
    if (file.bad()) throw std::runtime_error("read error in " + path);
    return text.str();
}
struct ConfigFile { std::string text; YAML::Node yaml; };
inline ConfigFile ReadConfigFile(const std::string &path)
{
    static const int status = gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so", "", true);
    if (status < 0) throw std::runtime_error("cannot load yaml-cpp; source the analysis environment");
    const std::string text = ReadText(path);
    return {text, YAML::Load(text)};
}
inline std::string RequiredString(const YAML::Node &node, const char *key)
{
    const std::string value = Get<std::string>(node, key, "");
    if (value.empty()) throw std::runtime_error(std::string("missing config key ") + key);
    return value;
}
inline std::unique_ptr<TFile> OpenRoot(const std::string &path, const char *mode = "READ")
{
    std::unique_ptr<TFile> file(TFile::Open(path.c_str(), mode));
    if (!file || file->IsZombie() || file->TestBit(TFile::kRecovered))
        throw std::runtime_error("cannot open intact ROOT file " + path);
    return file;
}

/////////////////////////////////////
// Input parts
/////////////////////////////////////
// part_001234.root -> 1234. The number is the part's identity and seeds its response smearing,
// so it must not depend on where the file sits in a list.
inline int PartNumber(const std::string &path)
{
    const std::string name = std::filesystem::path(path).stem().string();
    if (name.rfind("part_", 0) != 0) throw std::runtime_error("input file name must look like part_NNNNNN.root: " + path);
    return std::stoi(name.substr(5));
}
struct InputPart { int number; std::string path; };
inline InputPart MakePart(const std::string &path)
{
    InputPart part;
    part.path = std::filesystem::canonical(path).string();
    part.number = PartNumber(part.path);
    return part;
}
// One path per line, as the producer's files.txt.
inline std::vector<InputPart> ReadParts(const std::string &list)
{
    std::istringstream input(ReadText(list));
    std::vector<InputPart> parts;
    for (std::string path; std::getline(input, path);)
        if (!path.empty()) parts.push_back(MakePart(path));
    if (parts.empty()) throw std::runtime_error("empty input list " + list);
    return parts;
}
// No input argument: the configured list of the product. Otherwise one ROOT file, or a list of them.
inline std::vector<InputPart> InputParts(const YAML::Node &photonjet, Product product, const std::string &input)
{
    std::vector<InputPart> parts;
    if (input.empty()) parts = ReadParts(RequiredString(photonjet, InputListKey(product)));
    else if (std::filesystem::path(input).extension() == ".root") parts.push_back(MakePart(input));
    else parts = ReadParts(input);
    // Two exports may share part numbers (the DI parts are listed with the photon-jet and the inclusive parts).
    std::set<std::string> paths;
    for (const InputPart &part : parts)
        if (!paths.insert(part.path).second) throw std::runtime_error("part listed twice: " + part.path);
    return parts;
}

/////////////////////////////////////
// Output file name and response prior
/////////////////////////////////////
// <directory>/<core>_<product>_<var_type>[.<tag>].root. An existing file is never replaced.
inline std::string OutputPath(const YAML::Node &output, const std::string &core, const std::string &product, const std::string &tag)
{
    const std::string variation = RequiredString(output, "var_type");
    const std::string name = core + "_" + product + "_" + variation + (tag.empty() ? "" : "." + tag) + ".root";
    const std::string path = (std::filesystem::path(RequiredString(output, "directory")) / name).string();
    if (std::filesystem::exists(path)) throw std::runtime_error("output already exists: " + path);
    return path;
}
inline std::unique_ptr<TF1> MakePrior(const Config &config)
{
    std::unique_ptr<TF1> prior = std::make_unique<TF1>("photonjet_prior", config.trw_formula.c_str(), config.trw_xmin, config.trw_xmax);
    prior->AddToGlobalList(false);
    if (!prior->IsValid() || prior->GetNpar() != int(config.trw_params.size()))
        throw std::runtime_error("invalid response prior formula or parameter count");
    for (size_t i = 0; i < config.trw_params.size(); ++i) prior->SetParameter(i, config.trw_params[i]);
    return prior;
}

/////////////////////////////////////
// One macro call
/////////////////////////////////////
struct Job {
    ConfigFile config_file;
    Config config;
    Product product;
    bool simulation, signal;
    std::string core, output;
    std::vector<InputPart> parts;
    std::unique_ptr<TF1> prior;   // truth-spectrum reweight of the response, signal simulation only
    struct Period { double lumi_fraction, single_fraction; std::unique_ptr<TH1> reweight; };
    std::vector<Period> periods;  // the run periods of the PPG12 blend, simulation only

    Job(const std::string &config_path, const std::string &product_name,
        const std::string &input, const std::string &tag, const std::string &core_name)
        : config_file(ReadConfigFile(config_path)), config(LoadConfig(config_file.yaml)), product(ParseProduct(product_name)),
          simulation(product != Product::Data), signal(product == Product::Signal), core(core_name)
    {
        parts = InputParts(config_file.yaml["photonjet"], product, input);
        output = OutputPath(config_file.yaml["output"], core, product_name, tag);
        if (signal && config.unfold_reweight) prior = MakePrior(config);
        for (const Config::Period &period : config.periods) {
            TH1 *reweight = simulation ? LoadTruthVertexReweight(period.reweight_file) : nullptr;
            if (simulation && !reweight) throw std::runtime_error("cannot load the truth-vertex reweight of period " + period.name);
            periods.push_back({period.lumi_fraction, period.single_fraction, std::unique_ptr<TH1>(reweight)});
        }
    }
    // PPG12's simulation weight: the sample's cross-section ratio times, summed over the run periods,
    // the period's luminosity fraction, its single- or double-interaction fraction and the truth-vertex
    // reweight of that period. Data events have weight 1. Without periods the sample weight alone is used.
    double Weight(const Event &event) const
    {
        if (!simulation) return 1.;
        if (periods.empty()) return event.sample_weight;
        double blend = 0;
        for (const Period &period : periods) {
            const double fraction = event.double_interaction ? 1 - period.single_fraction : period.single_fraction;
            blend += period.lumi_fraction * fraction * TruthVertexWeight(period.reweight.get(), event.truth_vertex_z, event.truth_mb_vertex_z);
        }
        return event.sample_weight * blend;
    }
    double Prior(double pt) const { return prior ? prior->Eval(std::clamp(pt, config.trw_clamp_min, config.trw_clamp_max)) : 1.; }
    // Simulation with sample weights: the DI export mixes photon and jet samples in one part and is listed for both
    // products, so each product keeps the events of its own kind of sample.
    bool OwnsSample(const Event &event) const
    {
        return !config.sample_weights || event.jet_sample == (product == Product::Inclusive);
    }
};

inline void WriteText(TFile &file, const char *name, const std::string &text)
{
    file.cd();
    TObjString object(text.c_str());
    if (object.Write(name, TObject::kOverwrite) <= 0) throw std::runtime_error(std::string("cannot write ") + name);
}
// Histograms plus the configuration and the list of inputs that made them.
inline void WriteOutput(const Job &job, const HistogramList &booked)
{
    const std::filesystem::path parent = std::filesystem::path(job.output).parent_path();
    if (!parent.empty()) std::filesystem::create_directories(parent);
    const std::string temporary = job.output + ".tmp." + std::to_string(getpid());
    try {
        std::unique_ptr<TFile> file = OpenRoot(temporary, "CREATE");
        booked.Write(*file);
        WriteText(*file, "config", job.config_file.text);
        WriteText(*file, "core", job.core);
        WriteText(*file, "product", ProductName(job.product));
        std::ostringstream manifest;
        for (const InputPart &part : job.parts) manifest << part.number << ' ' << part.path << '\n';
        WriteText(*file, "input_manifest", manifest.str());
        file->Close();
        if (file->TestBit(TFile::kWriteError)) throw std::runtime_error("ROOT output write failed");
        // Publish without replacing an existing result, including another job's result.
        if (link(temporary.c_str(), job.output.c_str()) != 0) throw std::runtime_error("cannot publish " + job.output);
        std::filesystem::remove(temporary);
    } catch (...) { std::filesystem::remove(temporary); throw; }
}
inline void Fail(const std::exception &error)
{
    std::cerr << "photonjet: " << error.what() << '\n';
    gSystem->Exit(1);
}
} // namespace PJ
