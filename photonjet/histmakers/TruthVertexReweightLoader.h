#pragma once

#include <cmath>
#include <iostream>
#include <string>
#include <TFile.h>
#include <TH1.h>
#include <TNamed.h>

// Sentinel written by CaloAna24 when a truth vertex is unavailable; weight
// of 1 means "do not modify the event". Values strictly above the threshold
// are real coordinates and get the full reweight treatment.
constexpr float TRUTH_VTX_SENTINEL_THRESHOLD = -9000.0f;

// Per-vertex truth-vertex weight: 1 for sentinels, 0 for out-of-range
// (drops the event), interpolated value otherwise. Falls back to 0 if the
// interpolation returns a non-finite or non-positive value.
inline float TruthVertexWeight(TH1* h, float z)
{
    if (!std::isfinite(z) || z <= TRUTH_VTX_SENTINEL_THRESHOLD) return 1.0f;
    if (z < h->GetXaxis()->GetXmin() || z > h->GetXaxis()->GetXmax()) return 0.0f;
    float w = h->Interpolate(z);
    return (std::isfinite(w) && w > 0.0f) ? w : 0.0f;
}

// Double-MC overload: w(z_hard) * w(z_mb).
inline float TruthVertexWeight(TH1* h, float z_hard, float z_mb)
{
    return TruthVertexWeight(h, z_hard) * TruthVertexWeight(h, z_mb);
}

// Opens `file`, clones h_w_iterative (detached from any directory), logs the
// period/run_min/run_max/z_cut TNamed metadata, closes the file, and returns
// the clone. Returns nullptr on any failure (empty path, missing file,
// missing histogram); callers must check and bail out.
inline TH1* LoadTruthVertexReweight(const std::string& file)
{
    if (file.empty())
    {
        std::cerr << "[TruthVertexReweight] FATAL: truth_vertex_reweight_file empty" << std::endl;
        return nullptr;
    }
    TFile* f = TFile::Open(file.c_str(), "READ");
    if (!f || f->IsZombie())
    {
        std::cerr << "[TruthVertexReweight] FATAL: cannot open " << file << std::endl;
        return nullptr;
    }
    TH1* hw = dynamic_cast<TH1*>(f->Get("h_w_iterative"));
    if (!hw)
    {
        std::cerr << "[TruthVertexReweight] FATAL: h_w_iterative not found in " << file << std::endl;
        f->Close(); delete f;
        return nullptr;
    }
    TH1* clone = dynamic_cast<TH1*>(hw->Clone("h_truth_vtx_reweight_clone"));
    clone->SetDirectory(nullptr);
    for (const char* k : {"period", "run_min", "run_max", "z_cut"})
        if (auto* tn = dynamic_cast<TNamed*>(f->Get(k)))
            std::cout << "[TruthVertexReweight] " << k << "=" << tn->GetTitle() << std::endl;
    std::cout << "[TruthVertexReweight] loaded h_w_iterative from " << file
              << " (" << clone->GetNbinsX() << " bins, "
              << clone->GetBinWidth(1) << " cm/bin)" << std::endl;
    f->Close(); delete f;
    return clone;
}
