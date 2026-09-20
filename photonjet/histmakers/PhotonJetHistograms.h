// Booking and writing of histograms. Each macro declares and books its own histograms at the top.
#pragma once
#include "PhotonJetConfig.h"
#include <TDirectory.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <memory>

namespace PJ {
inline std::vector<double> UniformBins(int n, double lo, double hi)
{
    std::vector<double> edges(n + 1);
    for (int i = 0; i <= n; ++i) edges[i] = lo + (hi - lo) * i / n;
    return edges;
}

// Owns every histogram it books, so a macro writes them all with one call.
// A booked histogram has Sumw2 and belongs to no ROOT directory.
class HistogramList
{
public:
    TH1D *Book1D(const std::string &name, const std::string &title, const std::vector<double> &edges)
    {
        return Keep(new TH1D(name.c_str(), title.c_str(), edges.size() - 1, edges.data()));
    }
    TH2D *Book2D(const std::string &name, const std::string &title,
                 const std::vector<double> &x, const std::vector<double> &y)
    {
        return Keep(new TH2D(name.c_str(), title.c_str(), x.size() - 1, x.data(), y.size() - 1, y.data()));
    }
    // A profile keeps ROOT's default errors, as in the reference.
    TProfile *BookProfile(const std::string &name, const std::string &title, const std::vector<double> &edges)
    {
        TProfile *profile = new TProfile(name.c_str(), title.c_str(), edges.size() - 1, edges.data());
        profile->SetDirectory(nullptr);
        histograms_.emplace_back(profile);
        return profile;
    }
    void Write(TDirectory &output) const
    {
        TDirectory::TContext context(&output);
        for (const std::unique_ptr<TH1> &histogram : histograms_)
            if (histogram->Write() <= 0) throw std::runtime_error(std::string("cannot write histogram ") + histogram->GetName());
    }
private:
    template<class Histogram> Histogram *Keep(Histogram *histogram)
    {
        histogram->SetDirectory(nullptr);
        histogram->Sumw2();
        histograms_.emplace_back(histogram);
        return histogram;
    }
    std::vector<std::unique_ptr<TH1>> histograms_;
};

// Event bookkeeping shared by both macros: "cutflow", and the same counts per centrality bin.
constexpr int kRead = 0, kSelected = 1, kSelectedWeighted = 2, kNCutflowBins = 3;   // cutflow bins
struct EventCounter
{
    TH1D *cutflow = nullptr;
    TH2D *by_centrality = nullptr;   // x: centrality bin, events outside the configured classes go to underflow

    EventCounter(HistogramList &booked, const BinLayout &bins)
    {
        const char *labels[kNCutflowBins] = {"read", "selected", "selected weighted"};
        cutflow = booked.Book1D("cutflow", ";Selection;Events", UniformBins(kNCutflowBins, 0, kNCutflowBins));
        by_centrality = booked.Book2D("h_event_counts", ";Centrality bin;Selection",
                                      UniformBins(bins.nCentrality(), 0, bins.nCentrality()),
                                      UniformBins(kNCutflowBins, 0, kNCutflowBins));
        for (int i = 0; i < kNCutflowBins; ++i) {
            cutflow->GetXaxis()->SetBinLabel(i + 1, labels[i]);
            by_centrality->GetYaxis()->SetBinLabel(i + 1, labels[i]);
        }
    }
    void Read(int centrality_bin) const
    {
        cutflow->Fill(kRead + .5);
        by_centrality->Fill(centrality_bin + .5, kRead + .5);
    }
    void Selected(int centrality_bin, double weight) const
    {
        cutflow->Fill(kSelected + .5);
        cutflow->Fill(kSelectedWeighted + .5, weight);
        by_centrality->Fill(centrality_bin + .5, kSelected + .5);
        by_centrality->Fill(centrality_bin + .5, kSelectedWeighted + .5, weight);
    }
};
} // namespace PJ
