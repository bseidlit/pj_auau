// Plot consumers currently support the legacy single-cell histogram layout.
#pragma once
#include <TFile.h>
#include <TObjString.h>
#include <TSystem.h>
#include <yaml-cpp/yaml.h>
#include <memory>
#include <stdexcept>
#include <string>

namespace PJPlot {
inline TFile *Open(const char *path, const char *mode = "READ")
{
    std::unique_ptr<TFile> file(TFile::Open(path, mode));
    if (!file || file->IsZombie()) return file.release(); // Preserve caller error handling.
    auto *saved = dynamic_cast<TObjString *>(file->Get("photonjet_provenance"));
    if (saved) {
        static const int yaml_status = gSystem->Load("libyaml-cpp", "", true);
        if (yaml_status < 0) throw std::runtime_error("cannot load yaml-cpp to validate plotting input");
        const auto info = YAML::Load(saved->GetString().Data());
        const auto layout = info["layout"];
        if (layout && (layout["n_centrality"].as<int>() != 1 || layout["n_eta"].as<int>() != 1))
            throw std::runtime_error(std::string("plotting does not support multiple centrality/eta cells; explicit cell selection is required: ") + path);
    }
    return file.release();
}
} // namespace PJPlot
