// Shared file handling for histogram production and correction import.
#pragma once
#include "PhotonJetConfig.h"
#include "PhotonJetPhysics.h"
#include "PhotonJetHistograms.h"
#include <TF1.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TMD5.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TUUID.h>
#include <cctype>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <unistd.h>
#include <sys/stat.h>

namespace PJ
{
inline std::string ReadText(const std::string &path)
{
    std::ifstream file(path);
    if (!file) throw std::runtime_error("cannot read " + path);
    std::ostringstream text;
    text << file.rdbuf();
    if (file.bad()) throw std::runtime_error("read error in " + path);
    return text.str();
}

struct ConfigFile
{
    std::string path, text;
    YAML::Node yaml;
};
inline ConfigFile ReadConfig(const std::string &path)
{
    // Permanent, once-per-process loading also supports several runs in PyROOT.
    // Re-loading a version-suffixed .so as a macro can unload Cling declarations.
    static const int yaml_status = gSystem->Load("libyaml-cpp", "", true);
    if (yaml_status < 0)
        throw std::runtime_error("cannot load yaml-cpp; source the sPHENIX setup first");
    ConfigFile config;
    config.path = path.empty() ? std::string(gSystem->DirName(__FILE__)) + "/../configs/pp/config_pj_pp_nom.yaml" : path;
    config.text = ReadText(config.path);
    config.yaml = YAML::Load(config.text);
    if (!config.yaml.IsMap()) throw std::runtime_error("configuration must be a YAML mapping");
    return config;
}

inline std::string RequiredString(const YAML::Node &node, const char *key)
{
    const std::string value = Y<std::string>(node, key, "");
    if (value.empty()) throw std::runtime_error(std::string("missing or empty config key ") + key);
    return value;
}
inline void ValidateLabel(const std::string &label, const char *name)
{
    if (label.empty() || label == "." || label == "..") throw std::runtime_error(std::string("invalid ") + name);
    for (unsigned char ch : label)
        if (!std::isalnum(ch) && ch != '_' && ch != '-' && ch != '.')
            throw std::runtime_error(std::string(name) + " may contain only letters, digits, _, - and .");
}
inline std::string TaggedPath(const std::string &name, const std::string &tag)
{
    if (tag.empty()) return name;
    ValidateLabel(tag, "chunk tag");
    std::filesystem::path path(name);
    return (path.parent_path() / "chunks" / (path.stem().string() + "." + tag + ".root")).string();
}
inline std::string MainOutputPath(const YAML::Node &cfg, Product product, const std::string &tag = "")
{
    const YAML::Node out = cfg["output"];
    const std::string variation = RequiredString(out, "var_type");
    ValidateLabel(variation, "var_type");
    const char *key = product == Product::Data ? "data_outfile" : "eff_outfile";
    const std::string jet = product == Product::Inclusive ? "_jet" : "";
    return TaggedPath(RequiredString(out, key) + jet + "_" + variation + ".root", tag);
}
struct OutputNames { std::string main, response; };
inline OutputNames OutputPaths(const YAML::Node &cfg, Product product, const std::string &tag = "")
{
    OutputNames names{MainOutputPath(cfg, product, tag), ""};
    if (product == Product::Signal)
    {
        const YAML::Node out = cfg["output"];
        names.response = TaggedPath(RequiredString(out, "response_outfile") + "_" + RequiredString(out, "var_type") + ".root", tag);
        if (std::filesystem::absolute(names.main).lexically_normal() == std::filesystem::absolute(names.response).lexically_normal())
            throw std::runtime_error("main and response output paths must differ");
    }
    return names;
}

struct ChunkEntry { int index; std::string path; };
// Master files.txt: derive canonical zero-based indices for bare paths.
// Explicit chunks: preserve supplied indices; bare paths remain usable when
// sample-map lookup is not required. Never renumber an arbitrary subset.
inline std::vector<ChunkEntry> ReadParts(const std::string &path, bool master_list)
{
    std::istringstream input(ReadText(path));
    std::vector<ChunkEntry> parts;
    std::set<std::string> paths;
    std::set<int> indices;
    std::string line;
    int line_number = 0;
    while (std::getline(input, line))
    {
        ++line_number;
        line = line.substr(0, line.find('#'));
        std::istringstream row(line);
        std::string first, second, extra;
        if (!(row >> first)) continue;
        int index = master_list ? static_cast<int>(parts.size()) : -1;
        std::string file = first;
        if (row >> second)
        {
            size_t used = 0;
            try { index = std::stoi(first, &used); }
            catch (...) { throw std::runtime_error("invalid part index at " + path + ":" + std::to_string(line_number)); }
            if (used != first.size() || index < 0 || (row >> extra))
                throw std::runtime_error("expected 'index path' at " + path + ":" + std::to_string(line_number));
            file = second;
        }
        if (!paths.insert(file).second || (index >= 0 && !indices.insert(index).second))
            throw std::runtime_error("duplicate input path or index in " + path);
        parts.push_back({index, file});
    }
    if (parts.empty()) throw std::runtime_error("empty input list " + path);
    return parts;
}
inline std::string InputManifest(const std::vector<ChunkEntry> &parts)
{
    std::ostringstream text;
    for (const auto &part : parts) text << part.index << ' ' << part.path << '\n';
    return text.str();
}

inline std::unique_ptr<TFile> OpenRoot(const std::string &path, const char *mode = "READ")
{
    std::unique_ptr<TFile> file(TFile::Open(path.c_str(), mode));
    if (!file || file->IsZombie() || file->TestBit(TFile::kRecovered))
        throw std::runtime_error("cannot open intact ROOT file " + path + " (" + mode + ")");
    return file;
}

// Write next to the final file and rename only after a successful close.
// Exceptions remove our temporary file; an existing final file is preserved.
class RootOutput
{
public:
    explicit RootOutput(const std::string &path, bool copy_existing = false) : final_(path)
    {
        const auto parent = std::filesystem::path(path).parent_path();
        if (!parent.empty()) std::filesystem::create_directories(parent);
        std::string pattern = path + ".tmp.XXXXXX";
        std::vector<char> buffer(pattern.begin(), pattern.end());
        buffer.push_back('\0');
        const int fd = mkstemp(buffer.data());
        if (fd < 0) throw std::runtime_error("cannot create temporary output next to " + path);
        close(fd);
        temporary_ = buffer.data();
        try
        {
            if (copy_existing) std::filesystem::copy_file(path, temporary_, std::filesystem::copy_options::overwrite_existing);
            file_ = OpenRoot(temporary_, copy_existing ? "UPDATE" : "RECREATE");
            if (!file_->IsWritable()) throw std::runtime_error("ROOT output is not writable: " + temporary_);
        }
        catch (...)
        {
            std::remove(temporary_.c_str());
            throw;
        }
    }
    ~RootOutput()
    {
        if (file_) file_->Close();
        if (!temporary_.empty()) std::remove(temporary_.c_str());
    }
    TFile *get() const { return file_.get(); }
    void Close()
    {
        // Explicit writers may have already flushed every object: zero additional
        // bytes is valid for a nonempty file. Write/close error bits remain fatal.
        if (file_->Write("", TObject::kOverwrite) < 0 || file_->TestBit(TFile::kWriteError) ||
            file_->GetListOfKeys()->GetEntries() == 0)
            throw std::runtime_error("ROOT write failed: " + temporary_);
        file_->Close();
        if (file_->TestBit(TFile::kWriteError)) throw std::runtime_error("ROOT close failed: " + temporary_);
    }
    void Commit()
    {
        if (file_->IsOpen()) throw std::runtime_error("close ROOT output before publishing it");
        if (std::rename(temporary_.c_str(), final_.c_str()) != 0)
            throw std::runtime_error("cannot publish " + final_);
        temporary_.clear();
    }
private:
    std::string final_, temporary_;
    std::unique_ptr<TFile> file_;
};

inline void WriteText(TFile *file, const char *name, const std::string &text)
{
    file->cd();
    TObjString object(text.c_str());
    if (object.Write(name, TObject::kOverwrite) <= 0) throw std::runtime_error(std::string("cannot write ") + name);
}
inline std::string FileDigest(const std::string &path)
{
    // ROOT's FileChecksum restores timestamps with second precision. Stream the
    // bytes ourselves so hashing never changes the source's nanosecond mtime.
    std::ifstream file(path, std::ios::binary);
    if (!file) throw std::runtime_error("cannot fingerprint " + path);
    TMD5 digest;
    char buffer[64 * 1024];
    while (file) {
        file.read(buffer, sizeof(buffer));
        digest.Update(reinterpret_cast<const UChar_t *>(buffer), static_cast<UInt_t>(file.gcount()));
    }
    if (file.bad() || !file.eof()) throw std::runtime_error("read error while fingerprinting " + path);
    digest.Final();
    return digest.AsString();
}
inline YAML::Node CodeDigests()
{
    YAML::Node code;
    const auto base = std::filesystem::path(__FILE__).parent_path().parent_path();
    for (const char *name : {"histmakers/PhotonJetHistMaker.C", "histmakers/PhotonJetReader.h",
                            "histmakers/PhotonJetPhysics.h", "histmakers/PhotonJetHistograms.h",
                            "histmakers/PhotonJetConfig.h", "histmakers/PhotonJetIO.h", "support/CrossSectionWeights.h"})
        code[name] = FileDigest((base / name).string());
    return code;
}
inline void WriteRunMetadata(TFile *file, const ConfigFile &config, const YAML::Node &info,
                             const std::vector<ChunkEntry> &parts, size_t processed)
{
    if (processed != parts.size()) throw std::runtime_error("refusing to publish an incomplete input set");
    WriteText(file, "config", config.text); // The text actually parsed at startup, never reread at completion.
    WriteText(file, "photonjet_provenance", YAML::Dump(info));
    WriteText(file, "input_manifest", InputManifest(parts));
    file->cd();
    TH1D accounting("h_pj_input_parts", "Input completeness", 3, 0, 3);
    accounting.GetXaxis()->SetBinLabel(1, "expected");
    accounting.GetXaxis()->SetBinLabel(2, "processed");
    accounting.GetXaxis()->SetBinLabel(3, "failed");
    accounting.SetBinContent(1, parts.size());
    accounting.SetBinContent(2, processed);
    if (accounting.Write("h_pj_input_parts", TObject::kOverwrite) <= 0)
        throw std::runtime_error("cannot write input accounting");
}

inline const std::vector<std::string> &MbdHistogramNames()
{
    static const std::vector<std::string> names = {
        "h_truth_pT_vertexcut_0", "h_truth_pT_vertexcut_mbd_cut_0",
        "h_truth_pT_vertexcut_mbd_north_cut_0", "h_truth_pT_vertexcut_mbd_south_cut_0",
        "h_truth_pT_vertexcut_mbd_only_north_0", "h_truth_pT_vertexcut_mbd_only_south_0",
        "h_truth_pT_vertexcut_mbd_neither_0"};
    return names;
}
inline TH1 *RequireHistogram(TFile *file, const std::string &name, const std::vector<double> &edges)
{
    auto *hist = dynamic_cast<TH1 *>(file->Get(name.c_str()));
    if (!hist || hist->GetDimension() != 1 || hist->GetNbinsX() + 1 != static_cast<int>(edges.size()))
        throw std::runtime_error("missing or incompatible " + name + " in " + file->GetName());
    for (int i = 0; i <= hist->GetNbinsX(); ++i)
        if (hist->GetXaxis()->GetBinLowEdge(i + 1) != edges[i])
            throw std::runtime_error("axis mismatch for " + name + " in " + file->GetName());
    for (int i = 0; i < hist->GetNcells(); ++i)
        if (!std::isfinite(hist->GetBinContent(i)) || !std::isfinite(hist->GetBinError(i)) || hist->GetBinContent(i) < 0)
            throw std::runtime_error("invalid correction bins in " + name);
    return hist;
}
// Original sources and their identities are checked before creating output files.
inline constexpr const char *OriginalSchemaVersion = "PhotonJetTrees_v1";
inline constexpr int PlanFormatVersion = 4, ProvenanceFormatVersion = 4, CompletionFormatVersion = 3;
inline constexpr const char *InputVerification = "worker-md5-stat-v1";
struct Inputs {
    std::vector<ChunkEntry> originals;
    std::map<int, std::string> ids;
    YAML::Node records, dependencies, job_plan;
    std::string release, product, job_tag;
    ULong64_t events = 0;
};
inline std::string ReleaseBasename(const std::string &value) {
    // Match pathlib.Path.name without resolving directory aliases or '..'.
    auto path = std::filesystem::path(value);
    while (!path.empty() && (path.filename().empty() || path.filename() == ".")) {
        const auto parent = path.parent_path();
        if (parent == path) break; // Root has no release label; ValidateLabel rejects it.
        path = parent;
    }
    return path.filename().string();
}
inline std::string ReleaseProduct(Product kind, const std::string &system) {
    if (kind == Product::Data) return "data/" + system;
    return "simulation/" + (system == "auau" ? std::string("auau_embedded_") : std::string("pp_")) +
           (kind == Product::Signal ? "photonjet" : "inclusive");
}
inline bool LayoutMatches(const YAML::Node &a, const YAML::Node &b) {
    if (!a.IsMap() || !b.IsMap()) return false;
    for (const char *key : {"format_version", "n_centrality", "n_eta"})
        if (!a[key] || !b[key] || a[key].as<int>() != b[key].as<int>()) return false;
    for (const char *key : {"centrality_policy", "eta_policy"})
        if (!a[key] || !b[key] || a[key].as<std::string>() != b[key].as<std::string>()) return false;
    for (const char *key : {"centrality_edges", "eta_edges"})
        if (!a[key] || !b[key] || a[key].as<std::vector<double>>() != b[key].as<std::vector<double>>()) return false;
    return true;
}
inline void RequireSingleCellLayout(const YAML::Node &layout, const std::string &consumer) {
    if (!layout || layout["n_centrality"].as<int>() != 1 || layout["n_eta"].as<int>() != 1)
        throw std::runtime_error(consumer + " does not support multiple centrality/eta cells; explicit cell selection and normalization are required");
}
// Match Python's stat().st_mtime_ns; filesystem::file_time_type has a different epoch.
struct SourceStat { uintmax_t bytes; int64_t mtime_ns; };
inline SourceStat StatOriginal(const std::string &path) {
    struct stat value{};
    if (::stat(path.c_str(), &value) != 0 || !S_ISREG(value.st_mode))
        throw std::runtime_error("original source is missing or not a regular file: " + path);
    return {static_cast<uintmax_t>(value.st_size),
            int64_t(value.st_mtim.tv_sec) * 1000000000LL + value.st_mtim.tv_nsec};
}
inline void CheckOriginalStat(const YAML::Node &record) {
    const auto path = record["source_path"].as<std::string>();
    const auto current = StatOriginal(path);
    if (current.bytes != record["source_bytes"].as<uintmax_t>() ||
        current.mtime_ns != record["source_mtime_ns"].as<int64_t>())
        throw std::runtime_error("original source size or timestamp changed: " + path);
}
// Plans contain assignment/stat fields only. Workers add hashes, UUIDs and counts.
inline bool MatchesPlannedOriginal(const YAML::Node &record, const YAML::Node &planned) {
    for (const char *key : {"original_part_index", "source_bytes", "source_mtime_ns"})
        if (!record[key] || !planned[key] || record[key].as<int64_t>() != planned[key].as<int64_t>()) return false;
    for (const char *key : {"original_part_id", "source_path"})
        if (!record[key] || !planned[key] || record[key].as<std::string>() != planned[key].as<std::string>()) return false;
    return true;
}
inline YAML::Node ReadCheckedPlan(const std::string &path, const ConfigFile &config, Product kind) {
    const auto plan = YAML::Load(ReadText(path));
    if (plan["format_version"].as<int>() != PlanFormatVersion ||
        plan["input_verification"].as<std::string>() != InputVerification ||
        plan["backend"].as<std::string>() != BackendVersion ||
        plan["schema_version"].as<std::string>() != OriginalSchemaVersion ||
        plan["config_md5"].as<std::string>() != FileDigest(config.path) ||
        !LayoutMatches(plan["layout"], BinLayout(LoadCuts(config.yaml)).Metadata()))
        throw std::runtime_error("incompatible or stale job plan; regenerate chunks");
    const auto code = CodeDigests();
    if (plan["code"].size() != code.size()) throw std::runtime_error("job plan source inventory differs");
    for (const auto &entry : code)
        if (plan["code"][entry.first.as<std::string>()].as<std::string>() != entry.second.as<std::string>())
            throw std::runtime_error("maker code changed since chunk planning");
    for (const auto &entry : plan["list_md5"])
        if (FileDigest(entry.first.as<std::string>()) != entry.second.as<std::string>())
            throw std::runtime_error("job list changed since planning");
    const auto product = plan["products"][ProductName(kind)];
    const auto master = RequiredString(config.yaml["photonjet"], InputListKey(kind));
    if (product["master_md5"].as<std::string>() != FileDigest(master) ||
        !std::filesystem::equivalent(product["master"].as<std::string>(), master))
        throw std::runtime_error("master list changed since planning");
    for (const auto &dependency : product["dependencies"])
        if (FileDigest(dependency.first.as<std::string>()) != dependency.second.as<std::string>())
            throw std::runtime_error("dependency changed since planning: " + dependency.first.as<std::string>());
    return plan;
}
inline YAML::Node InspectOriginal(const ChunkEntry &part, const std::string &id, bool simulation) {
    YAML::Node record;
    record["original_part_index"] = part.index; record["original_part_id"] = id;
    record["source_path"] = part.path;
    const auto stat = StatOriginal(part.path);
    record["source_bytes"] = stat.bytes; record["source_mtime_ns"] = stat.mtime_ns;
    record["source_md5"] = FileDigest(part.path);
    auto file = OpenRoot(part.path);
    record["source_uuid"] = file->GetUUID().AsString();
    auto count = [&](const char *name, bool required) -> Long64_t {
        auto *object = file->Get(name);
        auto *tree = dynamic_cast<TTree *>(object);
        if ((required && !tree) || (object && !tree))
            throw std::runtime_error(std::string("missing or wrong-type original tree ") + name + ": " + part.path);
        return tree ? tree->GetEntries() : 0;
    };
    record["events"] = count("events", true); record["photons"] = count("photons", true);
    record["truth_photons"] = count("truthPhotons", simulation);
    record["source_link_rows"] = count("recoTruthLinks", false); record["truth_jet_rows"] = count("truthJets", false);
    record["links_available"] = file->Get("recoTruthLinks") != nullptr;
    record["truth_jets_available"] = file->Get("truthJets") != nullptr;
    CheckOriginalStat(record);
    return record;
}
// This full content check runs once at worker completion, never in status/merge/yield.
inline void CheckOriginalContents(const YAML::Node &record) {
    const auto path = record["source_path"].as<std::string>();
    CheckOriginalStat(record);
    if (FileDigest(path) != record["source_md5"].as<std::string>())
        throw std::runtime_error("original source fingerprint changed: " + path);
    CheckOriginalStat(record);
}
inline Inputs ReadInputs(const std::string &requested, const ConfigFile &config, Product kind, const Cuts &cuts) {
    Inputs out;
    const auto master = RequiredString(config.yaml["photonjet"], InputListKey(kind));
    const auto master_parts = ReadParts(master, true);
    for (size_t i=0; i<master_parts.size(); ++i)
        if (master_parts[i].index != int(i) || !std::filesystem::path(master_parts[i].path).is_absolute())
            throw std::runtime_error("master list needs canonical zero-based indices and absolute original paths");
    const auto path = requested.empty() ? master : requested;
    if (std::filesystem::path(path).extension() == ".json")
        throw std::runtime_error("prepared manifests are retired; provide an indexed original ROOT list");
    const bool use_master = std::filesystem::equivalent(path, master);
    auto selected = ReadParts(path, use_master);
    out.product = ReleaseProduct(kind, cuts.system);
    const char *release_override = gSystem->Getenv("PJ_RELEASE_ID");
    out.release = release_override && *release_override ? release_override :
        ReleaseBasename(RequiredString(config.yaml["photonjet"], "release"));
    out.dependencies[path] = FileDigest(path); out.dependencies[master] = FileDigest(master);
    out.dependencies[config.path] = FileDigest(config.path);
    if (!cuts.run_list_file.empty()) out.dependencies[cuts.run_list_file] = FileDigest(cuts.run_list_file);
    if (ReadText(config.path) != config.text) throw std::runtime_error("config changed after it was read");
    YAML::Node assignment(YAML::NodeType::Undefined);
    const auto plan_path = std::filesystem::path(path).parent_path() / "jobs.json";
    if (!use_master && std::filesystem::exists(plan_path)) {
        const auto plan = ReadCheckedPlan(plan_path.string(), config, kind);
        const auto product_plan = plan["products"][ProductName(kind)];
        for (const auto &dependency : product_plan["dependencies"])
            out.dependencies[dependency.first.as<std::string>()] = dependency.second;
        for (const auto &job : product_plan["jobs"])
            if (std::filesystem::weakly_canonical(job["input"].as<std::string>()) == std::filesystem::weakly_canonical(path)) {
                if (assignment) throw std::runtime_error("duplicate input assignment in job plan");
                assignment = YAML::Clone(job); out.job_tag = job["tag"].as<std::string>();
            }
        if (!assignment) throw std::runtime_error("input list is not assigned in the sealed job plan");
        out.release = plan["release_id"].as<std::string>();
        out.dependencies[plan_path.string()] = FileDigest(plan_path.string());
        out.job_plan["path"] = std::filesystem::weakly_canonical(plan_path).string();
        out.job_plan["md5"] = out.dependencies[plan_path.string()];
        YAML::Node binding;
        binding["input"] = assignment["input"]; binding["tag"] = assignment["tag"];
        out.job_plan["assignments"].push_back(binding);
    }
    ValidateLabel(out.release, "release identity");
    std::set<std::string> uuids, canonical_paths;
    for (auto part : selected) {
        if (part.index < 0 || part.index >= int(master_parts.size()) || master_parts[part.index].index != part.index ||
            std::filesystem::weakly_canonical(part.path) != std::filesystem::weakly_canonical(master_parts[part.index].path))
            throw std::runtime_error("input must use canonical original indices and paths from its master list");
        part.path = master_parts[part.index].path;
        const auto id = out.release + "/" + out.product + "/part/" + std::to_string(part.index);
        if (!out.ids.emplace(part.index, id).second || !canonical_paths.insert(std::filesystem::weakly_canonical(part.path).string()).second)
            throw std::runtime_error("duplicate original input coverage");
        if (assignment) {
            size_t matches = 0;
            for (const auto &sealed : assignment["parts"]) if (sealed["original_part_index"].as<int>() == part.index) {
                CheckOriginalStat(sealed); ++matches;
            }
            if (matches != 1) throw std::runtime_error("input is not uniquely assigned in job plan");
        }
        const auto record = InspectOriginal(part, id, kind != Product::Data);
        if (!uuids.insert(record["source_uuid"].as<std::string>()).second) throw std::runtime_error("duplicate original ROOT UUID");
        if (kind != Product::Data && cuts.match_source == "links" && !record["links_available"].as<bool>())
            throw std::runtime_error("requested photon links are unavailable");
        if (assignment) {
            size_t matches = 0;
            for (const auto &sealed : assignment["parts"]) if (MatchesPlannedOriginal(record, sealed)) ++matches;
            if (matches != 1) throw std::runtime_error("original input differs from its sealed job plan");
        }
        out.originals.push_back(part); out.records.push_back(record); out.events += record["events"].as<ULong64_t>();
    }
    if (assignment && assignment["parts"].size() != out.originals.size()) throw std::runtime_error("incomplete assigned input coverage");
    return out;
}
inline void RecordDependency(Inputs &inputs, const std::string &path) {
    const auto digest = FileDigest(path);
    if (inputs.dependencies[path] && inputs.dependencies[path].as<std::string>() != digest)
        throw std::runtime_error("dependency changed during setup: " + path);
    inputs.dependencies[path] = digest;
}
struct Model {
    PJ::Cuts cuts;
    bool simulation, signal;
    std::map<int, PPG12::SampleConfig> samples;
    VertexTable vertex;
    std::unique_ptr<TF1> prior;
};
inline Model Configure(const PJ::ConfigFile &config, PJ::Product kind, Inputs &inputs) {
    for (const auto &dependency : inputs.dependencies)
        if (FileDigest(dependency.first.as<std::string>()) != dependency.second.as<std::string>())
            throw std::runtime_error("dependency changed before configuration");
    Model model{}; model.cuts = PJ::LoadCuts(config.yaml);
    model.simulation = kind != PJ::Product::Data; model.signal = kind == PJ::Product::Signal;
    const auto &c = model.cuts;
    if (model.simulation && c.system == "auau" && c.weight_mode != "sample_map")
        throw std::runtime_error("Au+Au simulation requires sample_map weights");
    const bool need_sample = model.simulation && c.weight_mode != "stored";
    if (!c.run_list_file.empty()) RecordDependency(inputs, c.run_list_file);
    std::string sample_path = c.sample_map_file;
    const auto override = config.yaml["photonjet"][model.signal ? "sample_map_signal" : "sample_map_inclusive"];
    if (model.simulation && override) sample_path = override.as<std::string>();
    YAML::Node ranges;
    if (need_sample) {
        if (sample_path.empty()) throw std::runtime_error("sample_map weights need a sample map");
        RecordDependency(inputs, sample_path);
        auto doc = YAML::LoadFile(sample_path); ranges = doc["ranges"] ? doc["ranges"] : doc;
        std::set<int> covered;
        for (const auto &r : ranges) {
            const int first = r["first"].as<int>(), last = r["last"].as<int>();
            if (first < 0 || last < first || !PPG12::GetSampleConfig(r["sample"].as<std::string>()).valid)
                throw std::runtime_error("invalid sample map range");
            for (int i=first; i<=last; ++i) if (!covered.insert(i).second) throw std::runtime_error("overlapping sample map ranges");
        }
        if (!c.vertex_weight_file.empty()) {
            RecordDependency(inputs, c.vertex_weight_file);
            auto file = OpenRoot(c.vertex_weight_file);
            auto *h = dynamic_cast<TH1 *>(file->Get("h_vertex_reweight"));
            if (!h || h->GetDimension() != 1) throw std::runtime_error("missing or invalid vertex histogram");
            for (int i=1; i<=h->GetNbinsX(); ++i) { model.vertex.edges.push_back(h->GetXaxis()->GetBinLowEdge(i)); model.vertex.weights.push_back(h->GetBinContent(i)); }
            model.vertex.edges.push_back(h->GetXaxis()->GetBinUpEdge(h->GetNbinsX()));
        }
    }
    for (const auto &part : inputs.ids) {
        PPG12::SampleConfig sample;
        if (need_sample) for (const auto &r : ranges)
            if (part.first >= r["first"].as<int>() && part.first <= r["last"].as<int>()) sample = PPG12::GetSampleConfig(r["sample"].as<std::string>());
        if (need_sample && !sample.valid) throw std::runtime_error("no sample for canonical part " + std::to_string(part.first));
        if (need_sample && sample.isbackground && c.weight_mode == "sample_map")
            for (const auto &record : inputs.records) if (record["original_part_index"].as<int>() == part.first && !record["truth_jets_available"].as<bool>())
                throw std::runtime_error("sample stitching requires truth jets");
        model.samples.emplace(part.first, sample);
    }
    if (model.signal && c.unfold_reweight) {
        model.prior = std::make_unique<TF1>("pj_prior", c.trw_formula.c_str(), c.trw_xmin, c.trw_xmax);
        model.prior->AddToGlobalList(false);
        if (!model.prior->IsValid() || model.prior->GetNpar() != int(c.trw_params.size()))
            throw std::runtime_error("invalid response prior formula or parameter count");
        for (size_t p=0; p<c.trw_params.size(); ++p) model.prior->SetParameter(p, c.trw_params[p]);
    }
    return model;
}
inline YAML::Node Provenance(const PJ::ConfigFile &config, PJ::Product kind, const Inputs &inputs, const Model &model) {
    const auto &c = model.cuts;
    YAML::Node info;
    info["format_version"] = ProvenanceFormatVersion; info["input_verification"] = InputVerification; info["backend"] = BackendVersion; info["schema_version"] = OriginalSchemaVersion;
    info["rng_version"] = RNGVersion;
    info["stage"] = "histograms";
    info["run_id"] = TUUID().AsString(); info["product"] = PJ::ProductName(kind);
    info["config_path"] = config.path; info["ROOT_version"] = gROOT->GetVersion(); info["system"] = c.system;
    info["threshold_source"] = c.threshold_source; info["weight_mode"] = c.weight_mode; info["match_source"] = c.match_source;
    info["flag_check"] = c.flag_check; info["random_seed"] = c.random_seed;
    info["smearing"] = "response_only_truth_pt"; info["combined_efficiency"] = "any_single_candidate_tight_and_isolated";
    info["fingerprint_algorithm"] = "MD5"; info["release_id"] = inputs.release; info["release_product"] = inputs.product;
    info["layout"] = BinLayout(c).Metadata();
    if (inputs.job_plan.IsMap()) info["job_plan"] = inputs.job_plan;
    info["inputs"] = inputs.records; info["dependencies"] = inputs.dependencies; info["code"] = CodeDigests();
    return info;
}
// Directory aliases name the same destination. Keep the final filename literal:
// publication replaces that directory entry even if an older file is a symlink.
inline std::filesystem::path OutputIdentityPath(const std::string &name) {
    const auto path = std::filesystem::absolute(name).lexically_normal();
    return std::filesystem::weakly_canonical(path.parent_path()) / path.filename();
}
inline bool SameDestination(const std::string &a, const std::string &b) {
    if (OutputIdentityPath(a) == OutputIdentityPath(b)) return true;
    return std::filesystem::exists(a) && std::filesystem::exists(b) && std::filesystem::equivalent(a, b);
}
inline void CheckOutputDestinations(const ConfigFile &config, const Inputs &inputs, const OutputNames &names) {
    std::vector<std::string> destinations{names.main};
    if (!names.response.empty()) destinations.push_back(names.response);
    // The marker is published too; it must not replace a configured input.
    destinations.push_back(names.main + ".complete.yaml");
    for (size_t i=0; i<destinations.size(); ++i) for (size_t j=0; j<i; ++j)
        if (SameDestination(destinations[i], destinations[j])) throw std::runtime_error("output destinations must differ");
    std::vector<std::string> protected_paths{config.path};
    for (const auto &part : inputs.originals) protected_paths.push_back(part.path);
    for (const auto &dependency : inputs.dependencies) protected_paths.push_back(dependency.first.as<std::string>());
    const auto mbd = Y<std::string>(config.yaml["photonjet"], "external_mbd_eff_file", "");
    if (!mbd.empty()) protected_paths.push_back(mbd);
    for (const auto &destination : destinations) for (const auto &source : protected_paths)
        if (SameDestination(destination, source))
            throw std::runtime_error("output destination collides with an input or dependency: " + destination);

    // Check configured product outputs when those products are configured. A
    // data-only direct config need not supply unrelated simulation output keys.
    std::vector<std::string> product_outputs;
    const auto output = config.yaml["output"];
    if (!Y<std::string>(output, "data_outfile", "").empty())
        product_outputs.push_back(MainOutputPath(config.yaml, Product::Data));
    if (!Y<std::string>(output, "eff_outfile", "").empty()) {
        product_outputs.push_back(MainOutputPath(config.yaml, Product::Signal));
        product_outputs.push_back(MainOutputPath(config.yaml, Product::Inclusive));
        if (!Y<std::string>(output, "response_outfile", "").empty())
            product_outputs.push_back(OutputPaths(config.yaml, Product::Signal).response);
    }
    for (size_t i=0; i<product_outputs.size(); ++i) for (size_t j=0; j<i; ++j)
        if (SameDestination(product_outputs[i], product_outputs[j]))
            throw std::runtime_error("configured output destinations collide across products");
}
inline std::string CompletionPath(const std::string &main) { return main + ".complete.yaml"; }
inline void Publish(PJ::RootOutput &main, PJ::RootOutput *response, const PJ::OutputNames &names,
                    const YAML::Node &provenance, const std::string &stage = "histograms") {
    main.Close(); if (response) response->Close();
    YAML::Node marker;
    marker["format_version"] = CompletionFormatVersion; marker["input_verification"] = InputVerification; marker["run_id"] = provenance["run_id"]; marker["stage"] = stage;
    marker["product"] = provenance["product"]; marker["backend"] = BackendVersion;
    marker["layout"] = provenance["layout"];
    for (const auto &part : provenance["inputs"]) marker["original_part_ids"].push_back(part["original_part_id"]);
    auto record = [&](PJ::RootOutput &output, const std::string &path, const std::string &kind) {
        const std::string temp = output.get()->GetName();
        auto file = OpenRoot(temp);
        for (const char *key : {"config", "photonjet_provenance", "input_manifest", "h_pj_input_parts"})
            if (!file->Get(key)) throw std::runtime_error(std::string("missing output metadata ") + key);
        const int cells = provenance["layout"]["n_centrality"].as<int>() * provenance["layout"]["n_eta"].as<int>();
        if (file->GetListOfKeys()->GetEntries() != ExpectedRootObjectCount(cells, kind == "response"))
            throw std::runtime_error("incomplete output inventory");
        YAML::Node identity;
        identity["kind"] = kind; identity["path"] = OutputIdentityPath(path).string();
        identity["md5"] = PJ::FileDigest(temp); identity["bytes"] = std::filesystem::file_size(temp); identity["ROOT_UUID"] = file->GetUUID().AsString();
        marker["files"].push_back(identity);
    };
    record(main, names.main, "main"); if (response) record(*response, names.response, "response");
    const auto final = CompletionPath(names.main);
    const auto temporary = final + "." + provenance["run_id"].as<std::string>() + ".tmp";
    {
        std::ofstream stream(temporary);
        if (!(stream << YAML::Dump(marker) << "\n") || !stream.flush())
            throw std::runtime_error("cannot write completion marker " + temporary);
    }
    try {
        // Invalidate first. A failed second rename can never authorize a mixed pair.
        std::filesystem::remove(final);
        if (response) response->Commit(); main.Commit();
        std::filesystem::rename(temporary, final);
    } catch (...) { std::filesystem::remove(temporary); throw; }
}
// Verify worker-discovered fingerprints cover exactly the stat-only plan assignments.
inline void CheckPlanBinding(const YAML::Node &info, const ConfigFile &config, Product kind,
                             const std::string &main_path) {
    const auto stage = info["stage"].as<std::string>();
    if (stage != "histograms" && stage != "merged") throw std::runtime_error("unknown histogram stage");
    const auto binding = info["job_plan"];
    if (!binding) {
        if (stage == "merged") throw std::runtime_error("merged output has no job plan binding");
        return; // A direct unplanned macro invocation is also supported.
    }
    const auto path = binding["path"].as<std::string>();
    if (FileDigest(path) != binding["md5"].as<std::string>()) throw std::runtime_error("published job plan changed");
    const auto plan = ReadCheckedPlan(path, config, kind);
    if (plan["release_id"].as<std::string>() != info["release_id"].as<std::string>())
        throw std::runtime_error("published release differs from job plan");
    const auto jobs = plan["products"][ProductName(kind)]["jobs"];
    const auto assignments = binding["assignments"];
    if (!assignments.IsSequence() || assignments.size() == 0 ||
        (stage == "histograms" && assignments.size() != 1) ||
        (stage == "merged" && assignments.size() != jobs.size()))
        throw std::runtime_error("incomplete published job assignments");
    std::map<int, YAML::Node> records;
    for (const auto &record : info["inputs"])
        if (!records.emplace(record["original_part_index"].as<int>(), record).second)
            throw std::runtime_error("duplicate published original index");
    std::set<std::string> lists, tags;
    std::set<int> indices;
    for (const auto &assignment : assignments) {
        const auto input = std::filesystem::weakly_canonical(assignment["input"].as<std::string>()).string();
        const auto tag = assignment["tag"].as<std::string>();
        if (!lists.insert(input).second || !tags.insert(tag).second) throw std::runtime_error("duplicate published job assignment");
        size_t matches = 0;
        for (const auto &job : jobs) {
            if (std::filesystem::weakly_canonical(job["input"].as<std::string>()).string() != input ||
                job["tag"].as<std::string>() != tag) continue;
            ++matches;
            if (stage == "histograms" && OutputIdentityPath(main_path) != OutputIdentityPath(OutputPaths(config.yaml, kind, tag).main))
                throw std::runtime_error("chunk destination differs from its job assignment");
            for (const auto &part : job["parts"]) {
                if (!indices.insert(part["original_part_index"].as<int>()).second) throw std::runtime_error("duplicate planned original coverage");
                const auto found = records.find(part["original_part_index"].as<int>());
                if (found == records.end() || !MatchesPlannedOriginal(found->second, part))
                    throw std::runtime_error("published original differs from planned assignment");
            }
        }
        if (matches != 1) throw std::runtime_error("unplanned published job assignment");
    }
    if (indices.size() != info["inputs"].size()) throw std::runtime_error("extra published original coverage");
    if (stage == "merged" && OutputIdentityPath(main_path) != OutputIdentityPath(OutputPaths(config.yaml, kind).main))
        throw std::runtime_error("merged destination differs from configuration");
}
// Validate an exact published generation before local merge or yield consumes it.
inline YAML::Node CheckCompletion(const std::string &main_path, const PJ::ConfigFile &config, PJ::Product product) {
    const auto marker = YAML::Load(PJ::ReadText(CompletionPath(main_path)));
    if (marker["format_version"].as<int>() != CompletionFormatVersion ||
        marker["input_verification"].as<std::string>() != InputVerification || marker["backend"].as<std::string>() != BackendVersion ||
        marker["product"].as<std::string>() != PJ::ProductName(product)) throw std::runtime_error("incompatible completion marker");
    const auto layout = BinLayout(LoadCuts(config.yaml)).Metadata();
    if (!LayoutMatches(marker["layout"], layout)) throw std::runtime_error("completion histogram layout differs");
    const std::set<std::string> kinds = product == PJ::Product::Signal ? std::set<std::string>{"main","response"} : std::set<std::string>{"main"};
    std::set<std::string> seen_kinds, seen_paths;
    YAML::Node main_info;
    std::string common_provenance;
    const auto code = CodeDigests();
    for (const auto &entry : marker["files"]) {
        const auto kind = entry["kind"].as<std::string>(), path = entry["path"].as<std::string>();
        if (!kinds.count(kind) || !seen_kinds.insert(kind).second || !seen_paths.insert(OutputIdentityPath(path).string()).second ||
            (kind == "main" && OutputIdentityPath(main_path) != OutputIdentityPath(path)))
            throw std::runtime_error("invalid completion file membership");
        if (std::filesystem::file_size(path) != entry["bytes"].as<uintmax_t>() || PJ::FileDigest(path) != entry["md5"].as<std::string>())
            throw std::runtime_error("published file fingerprint mismatch: " + path);
        auto file = OpenRoot(path);
        if (file->GetUUID().AsString() != entry["ROOT_UUID"].as<std::string>()) throw std::runtime_error("published ROOT UUID mismatch");
        auto *saved = dynamic_cast<TObjString *>(file->Get("config"));
        auto *record = dynamic_cast<TObjString *>(file->Get("photonjet_provenance"));
        auto *parts = dynamic_cast<TH1 *>(file->Get("h_pj_input_parts"));
        auto *manifest = dynamic_cast<TObjString *>(file->Get("input_manifest"));
        if (!saved || saved->GetString().Data() != config.text || !record || !parts || !manifest)
            throw std::runtime_error("config differs or published metadata is missing");
        const std::string text = record->GetString().Data();
        if (!common_provenance.empty() && common_provenance != text) throw std::runtime_error("mixed main/response generations");
        common_provenance = text;
        const auto info = YAML::Load(text);
        if (info["run_id"].as<std::string>() != marker["run_id"].as<std::string>() ||
            info["product"].as<std::string>() != PJ::ProductName(product) || info["format_version"].as<int>() != ProvenanceFormatVersion ||
            info["input_verification"].as<std::string>() != InputVerification ||
            info["release_product"].as<std::string>() != ReleaseProduct(product, LoadCuts(config.yaml).system) ||
            info["ROOT_version"].as<std::string>() != gROOT->GetVersion() ||
            info["backend"].as<std::string>() != BackendVersion || info["schema_version"].as<std::string>() != OriginalSchemaVersion ||
            info["rng_version"].as<std::string>() != RNGVersion)
            throw std::runtime_error("incompatible published generation/backend/schema/RNG");
        if (!LayoutMatches(info["layout"], layout)) throw std::runtime_error("published histogram layout differs");
        if (info["stage"].as<std::string>() != marker["stage"].as<std::string>()) throw std::runtime_error("completion stage differs");
        if (info["code"].size() != code.size()) throw std::runtime_error("incompatible direct-maker source inventory");
        for (const auto &source : code)
            if (info["code"][source.first.as<std::string>()].as<std::string>() != source.second.as<std::string>())
                throw std::runtime_error("code differs from published generation: " + source.first.as<std::string>());
        for (const auto &dependency : info["dependencies"])
            if (PJ::FileDigest(dependency.first.as<std::string>()) != dependency.second.as<std::string>())
                throw std::runtime_error("dependency changed since production: " + dependency.first.as<std::string>());
        std::set<std::string> ids, uuids;
        std::set<int> indices;
        std::vector<PJ::ChunkEntry> originals;
        for (const auto &input : info["inputs"]) {
            const int index = input["original_part_index"].as<int>();
            const auto id = input["original_part_id"].as<std::string>();
            if (index < 0 || !indices.insert(index).second || !ids.insert(id).second || !uuids.insert(input["source_uuid"].as<std::string>()).second ||
                id != info["release_id"].as<std::string>() + "/" + info["release_product"].as<std::string>() + "/part/" + std::to_string(index))
                throw std::runtime_error("invalid published original-part identities");
            const auto digest = input["source_md5"].as<std::string>();
            if (digest.size() != 32 || digest.find_first_not_of("0123456789abcdef") != std::string::npos ||
                input["source_uuid"].as<std::string>().empty()) throw std::runtime_error("invalid worker source fingerprint");
            for (const char *key : {"events", "photons", "truth_photons", "source_link_rows", "truth_jet_rows"})
                if (input[key].as<Long64_t>() < 0) throw std::runtime_error("invalid worker source count");
            input["links_available"].as<bool>(); input["truth_jets_available"].as<bool>();
            originals.push_back({index, input["source_path"].as<std::string>()});
            if (kind == "main") CheckOriginalStat(input);
        }
        std::set<std::string> sealed;
        for (const auto &id : marker["original_part_ids"])
            if (!sealed.insert(id.as<std::string>()).second) throw std::runtime_error("duplicate completion coverage");
        if (ids.empty() || ids != sealed || parts->GetNbinsX() != 3 || parts->GetBinContent(1) != double(ids.size()) ||
            parts->GetBinContent(2) != double(ids.size()) || parts->GetBinContent(3) != 0 ||
            manifest->GetString().Data() != PJ::InputManifest(originals)) throw std::runtime_error("incomplete published input accounting");
        ULong64_t expected_events = 0;
        for (const auto &input : info["inputs"]) expected_events += input["events"].as<ULong64_t>();
        if (info["workers"].as<unsigned>() != 1 || info["events_processed"].as<ULong64_t>() != expected_events ||
            info["parts_processed"].as<size_t>() != ids.size()) throw std::runtime_error("incomplete serial processing accounting");
        const auto selected = info["selected_events_by_centrality"].as<std::vector<ULong64_t>>();
        ULong64_t selected_total = 0;
        for (auto count : selected) selected_total += count;
        if (selected.size() != size_t(layout["n_centrality"].as<int>()) || selected_total > expected_events ||
            !info["event_loops"] || info["event_loops"].as<ULong64_t>() == 0 ||
            (info["stage"].as<std::string>() == "histograms" && info["event_loops"].as<ULong64_t>() != 1))
            throw std::runtime_error("invalid selected-event or loop accounting");
        if (kind == "main") { CheckPlanBinding(info, config, product, main_path); main_info = info; }
    }
    if (seen_kinds != kinds || !main_info.IsMap()) throw std::runtime_error("incomplete published file pair");
    return main_info;
}
inline std::string CheckCompletionText(const std::string &main, const std::string &config, const std::string &product) {
    return YAML::Dump(CheckCompletion(main, PJ::ReadConfig(config), PJ::ParseProduct(product)));
}

inline void RecheckInputs(const Inputs &inputs, const YAML::Node &provenance) {
    for (const auto &record : inputs.records) CheckOriginalContents(record);
    for (const auto &dependency : inputs.dependencies)
        if (FileDigest(dependency.first.as<std::string>()) != dependency.second.as<std::string>())
            throw std::runtime_error("dependency changed during analysis: " + dependency.first.as<std::string>());
    if (YAML::Dump(CodeDigests()) != YAML::Dump(provenance["code"])) throw std::runtime_error("code changed during analysis");
}
} // namespace PJ
