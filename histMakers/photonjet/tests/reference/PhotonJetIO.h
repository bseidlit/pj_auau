// Shared file handling for histogram production and correction import.
#pragma once
#include "PhotonJetConfig.h"
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
    // Permanent, once-per-process loading also supports several graphs in PyROOT.
    // Re-loading a version-suffixed .so as a macro can unload Cling declarations.
    static const int yaml_status = gSystem->Load("libyaml-cpp", "", true);
    if (yaml_status < 0)
        throw std::runtime_error("cannot load yaml-cpp; source the sPHENIX setup first");
    ConfigFile config;
    config.path = path.empty() ? std::string(gSystem->DirName(__FILE__)) + "/../../configs/pp/config_pj_pp_nom.yaml" : path;
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
    if (!file || file->IsZombie()) throw std::runtime_error("cannot open ROOT file " + path + " (" + mode + ")");
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
    std::unique_ptr<TMD5> digest(TMD5::FileChecksum(path.c_str()));
    if (!digest) throw std::runtime_error("cannot fingerprint " + path);
    return digest->AsString();
}
inline YAML::Node CodeDigests()
{
    YAML::Node code;
    const std::string directory = gSystem->DirName(__FILE__);
    for (const char *name : {"PhotonJetHistMaker.C", "PhotonJetReader.h", "PhotonJetConfig.h",
                             "PhotonJetSelection.h", "PhotonJetWeights.h", "PhotonJetIO.h", "../../support/CrossSectionWeights.h"})
        code[name] = FileDigest(directory + "/" + name);
    return code;
}
inline YAML::Node RunProvenance(const ConfigFile &config, Product product, const Cuts &cuts,
                                const std::string &list, const std::vector<std::string> &dependencies)
{
    YAML::Node info;
    info["format_version"] = 1;
    info["run_id"] = TUUID().AsString();
    info["product"] = ProductName(product);
    info["config_path"] = config.path;
    info["ROOT_version"] = gROOT->GetVersion();
    info["system"] = cuts.system;
    info["threshold_source"] = cuts.threshold_source;
    info["weight_mode"] = cuts.weight_mode;
    info["match_source"] = cuts.match_source;
    info["flag_check"] = cuts.flag_check;
    info["random_seed"] = cuts.random_seed;
    info["smearing"] = "response_only_truth_pt";
    info["combined_efficiency"] = "any_single_candidate_tight_and_isolated";
    info["fingerprint_algorithm"] = "MD5"; // ROOT's built-in content fingerprint; not a security signature.
    info["input_list"] = list;
    info["dependencies"][list] = FileDigest(list);
    for (const auto &path : dependencies) if (!path.empty()) info["dependencies"][path] = FileDigest(path);
    info["code"] = CodeDigests();
    return info;
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
} // namespace PJ
