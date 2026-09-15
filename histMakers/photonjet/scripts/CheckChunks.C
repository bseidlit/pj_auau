#pragma once
// Enforce complete generations, direct-maker code/schema/RNG, and exact original coverage.
#include "../histmakers/PhotonJetIO.h"
#include <iostream>

namespace PJ {
inline void ValidateChunks(const std::string &manifest, const std::string &config_path, const std::string &product,
                 const std::string &job_plan = "") {
        const auto kind = PJ::ParseProduct(product);
        const auto config = PJ::ReadConfig(config_path);
        const auto master = PJ::RequiredString(config.yaml["photonjet"], PJ::InputListKey(kind));
        const auto master_parts = PJ::ReadParts(master, true);
        std::set<int> expected, seen;
        std::set<std::string> uuids, paths;
        YAML::Node jobs(YAML::NodeType::Undefined);
        if (!job_plan.empty()) {
            const auto plan = PJ::ReadCheckedPlan(job_plan, config, kind);
            for (const auto &index : plan["products"][product]["indices"])
                if (!expected.insert(index.as<int>()).second) throw std::runtime_error("duplicate planned part");
            jobs = plan["products"][product]["jobs"];
        } else for (const auto &part : master_parts) expected.insert(part.index);
        std::istringstream files(PJ::ReadText(manifest)); std::string path, release;
        size_t count=0;
        while (std::getline(files,path)) {
            if (path.empty()) continue;
            if (!paths.insert(OutputIdentityPath(path).string()).second) throw std::runtime_error("duplicate chunk path");
            const auto info = PJ::CheckCompletion(path,config,kind);
            if (info["event_loops"].as<unsigned>() != 1 || info["workers"].as<unsigned>() != 1 ||
                info["stage"].as<std::string>() != "histograms") throw std::runtime_error("chunk must have one serial original-input pass");
            if (!release.empty() && info["release_id"].as<std::string>() != release) throw std::runtime_error("incompatible release identities");
            release=info["release_id"].as<std::string>();
            std::set<int> chunk_indices;
            for (const auto &input : info["inputs"]) {
                const int index=input["original_part_index"].as<int>();
                if (index < 0 || index >= int(master_parts.size()) || master_parts[index].path != input["source_path"].as<std::string>() ||
                    !seen.insert(index).second || !uuids.insert(input["source_uuid"].as<std::string>()).second) throw std::runtime_error("duplicate or noncanonical chunk coverage");
                chunk_indices.insert(index);
            }
            if (jobs) {
                const auto binding = info["job_plan"];
                if (!binding || binding["assignments"].size() != 1 ||
                    std::filesystem::weakly_canonical(binding["path"].as<std::string>()) != std::filesystem::weakly_canonical(job_plan) ||
                    binding["md5"].as<std::string>() != FileDigest(job_plan))
                    throw std::runtime_error("chunk is not bound to this job plan");
                bool found=false;
                for (const auto &job : jobs) if (OutputIdentityPath(PJ::OutputPaths(config.yaml,kind,job["tag"].as<std::string>()).main) == OutputIdentityPath(path)) {
                    if (found) throw std::runtime_error("duplicate job output assignment");
                    found=true; std::set<int> planned;
                    const auto assignment = info["job_plan"]["assignments"][0];
                    if (assignment["tag"].as<std::string>() != job["tag"].as<std::string>() ||
                        std::filesystem::weakly_canonical(assignment["input"].as<std::string>()) != std::filesystem::weakly_canonical(job["input"].as<std::string>()))
                        throw std::runtime_error("chunk job assignment differs from plan");
                    for (const auto &part : job["parts"]) {
                        if (!planned.insert(part["original_part_index"].as<int>()).second) throw std::runtime_error("duplicate assigned original part");
                        bool same=false;
                        for (const auto &input : info["inputs"]) if (input["original_part_index"].as<int>() == part["original_part_index"].as<int>())
                            same = MatchesPlannedOriginal(input, part);
                        if (!same) throw std::runtime_error("chunk original assignment differs from job plan");
                    }
                    if (planned != chunk_indices) throw std::runtime_error("chunk does not cover its assigned original parts");
                }
                if (!found) throw std::runtime_error("unplanned chunk output");
            }
            ++count;
        }
        if (!count || seen != expected || (jobs && jobs.size() != count)) throw std::runtime_error("missing or extra original-part coverage");
        std::cout << "[CheckChunks] verified " << count << " complete " << product << " chunks; " << seen.size() << " disjoint original parts" << std::endl;
    }
} // namespace PJ

void CheckChunks(const std::string &manifest, const std::string &config_path, const std::string &product,
                 const std::string &job_plan = "") {
    try { PJ::ValidateChunks(manifest, config_path, product, job_plan); }
    catch (const std::exception &error) {
        std::cerr << "[CheckChunks] ERROR: " << error.what() << std::endl; gSystem->Exit(1);
    }
}
