// Run once on a hadd-merged file. Checks that the merge kept every event, records the full
// input manifest, and builds the efficiencies and the response from the complete counts.
#include "PhotonJetIO.h"
#include <RooUnfoldResponse.h>
#include <TEfficiency.h>
#include <TH2.h>
#include <TKey.h>

namespace {
std::string Text(TFile &file, const char *name)
{
    TObjString *text = dynamic_cast<TObjString *>(file.Get(name));
    if (!text) throw std::runtime_error(std::string("missing ") + name + " in " + file.GetName());
    return text->GetString().Data();
}
double EventsRead(TFile &file)
{
    TH1 *cutflow = dynamic_cast<TH1 *>(file.Get("cutflow"));
    if (!cutflow) throw std::runtime_error(std::string("missing cutflow in ") + file.GetName());
    return cutflow->GetBinContent(1);
}
TH1 *Counts(TFile &file, const std::string &name)
{
    TH1 *histogram = dynamic_cast<TH1 *>(file.Get(name.c_str()));
    if (!histogram) throw std::runtime_error("missing " + name);
    return histogram;
}
TH2 *Counts2D(TFile &file, const std::string &name)
{
    TH2 *histogram = dynamic_cast<TH2 *>(Counts(file, name));
    if (!histogram) throw std::runtime_error(name + " is not a TH2");
    return histogram;
}
// The reference's flattened x_Jgamma matrix holds the truth global bin on x and the reco global bin on y.
// RooUnfold wants reco on x and truth on y. Both count x_Jgamma fast and pT slow, from 1.
TH2D Transposed(const TH2 &matrix)
{
    TH2D result("transposed", "", matrix.GetNbinsY(), 0.5, matrix.GetNbinsY() + 0.5, matrix.GetNbinsX(), 0.5, matrix.GetNbinsX() + 0.5);
    result.Sumw2();
    for (int x = 1; x <= matrix.GetNbinsX(); ++x)
        for (int y = 1; y <= matrix.GetNbinsY(); ++y) {
            result.SetBinContent(y, x, matrix.GetBinContent(x, y));
            result.SetBinError(y, x, matrix.GetBinError(x, y));
        }
    return result;
}
void WriteResponse(TFile &file, const std::string &name, const TH1 *reco, const TH1 *truth, const TH2 *matrix)
{
    RooUnfoldResponse response(reco, truth, matrix, name.c_str(), "", false);
    if (file.WriteObject(&response, name.c_str()) <= 0) throw std::runtime_error("cannot write " + name);
}
void WriteEfficiency(TFile &file, const std::string &name, TH1 *passed, TH1 *total)
{
    if (!TEfficiency::CheckConsistency(*passed, *total, "w")) throw std::runtime_error("inconsistent efficiency counts: " + name);
    TEfficiency efficiency(*passed, *total);
    efficiency.SetName(name.c_str());
    efficiency.SetUseWeightedEvents();
    efficiency.SetStatisticOption(TEfficiency::kBUniform);
    if (file.WriteObject(&efficiency, name.c_str()) <= 0) throw std::runtime_error("cannot write " + name);
}
} // namespace

void FinalizeMerged(const std::string &merged_path, const std::string &chunk_list)
{
    try {
        /////////////////////////////////////
        // What the chunks hold
        /////////////////////////////////////
        std::istringstream chunks(PJ::ReadText(chunk_list));
        std::string path, manifest, core, product;
        double events_read = 0;
        while (std::getline(chunks, path)) {
            if (path.empty()) continue;
            std::unique_ptr<TFile> chunk = PJ::OpenRoot(path);
            if (core.empty()) { core = Text(*chunk, "core"); product = Text(*chunk, "product"); }
            else if (core != Text(*chunk, "core") || product != Text(*chunk, "product"))
                throw std::runtime_error("chunks mix cores or products: " + path);
            events_read += EventsRead(*chunk);
            manifest += Text(*chunk, "input_manifest");
        }
        if (core.empty()) throw std::runtime_error("empty chunk list " + chunk_list);

        /////////////////////////////////////
        // The merge must hold every event; hadd keeps only the first chunk's manifest
        /////////////////////////////////////
        std::unique_ptr<TFile> merged = PJ::OpenRoot(merged_path, "UPDATE");
        if (EventsRead(*merged) != events_read)
            throw std::runtime_error("merged file holds " + std::to_string(EventsRead(*merged)) + " events, chunks hold " + std::to_string(events_read));
        PJ::WriteText(*merged, "input_manifest", manifest);

        /////////////////////////////////////
        // Efficiencies and responses, once per centrality/eta cell
        /////////////////////////////////////
        if (core == "spectrum" && product != "data") {
            const std::string prefix = "h_truth_pT_selected_";
            std::vector<std::string> cells;
            for (TObject *key : *merged->GetListOfKeys())
                if (std::string(key->GetName()).rfind(prefix, 0) == 0) cells.push_back(std::string(key->GetName()).substr(prefix.size()));
            if (cells.empty()) throw std::runtime_error("no truth counts in " + merged_path);
            for (const std::string &cell : cells) {
                TH1 *truth = Counts(*merged, prefix + cell);
                TH1 *reconstructed = Counts(*merged, "h_truth_pT_reconstructed_" + cell);
                TH1 *isolated = Counts(*merged, "h_truth_pT_isolated_" + cell);
                TH1 *identified = Counts(*merged, "h_truth_pT_identified_" + cell);
                WriteEfficiency(*merged, "eff_reco_cent" + cell, reconstructed, truth);
                WriteEfficiency(*merged, "eff_iso_cent" + cell, isolated, reconstructed);
                WriteEfficiency(*merged, "eff_id_cent" + cell, identified, isolated);
                WriteEfficiency(*merged, "eff_all_cent" + cell, identified, truth);
                WriteEfficiency(*merged, "eff_vertex_cent" + cell, Counts(*merged, "h_truth_pT_vertexcut_" + cell), Counts(*merged, "h_truth_pT_" + cell));
                if (product != "sim_signal") continue;
                WriteResponse(*merged, "response_matrix_full_" + cell, Counts(*merged, "h_pT_reco_response_" + cell),
                              Counts(*merged, "h_pT_truth_response_" + cell), Counts2D(*merged, "h_response_full_" + cell));
                const TH2D matrix = Transposed(*Counts2D(*merged, "h_response_xjgamma_global_cent" + cell));
                WriteResponse(*merged, "response_matrix_xjgamma_2d_cent" + cell, Counts2D(*merged, "h_xjgamma_reco_response_cent" + cell),
                              Counts2D(*merged, "h_xjgamma_truth_response_cent" + cell), &matrix);
            }
        }
        merged->Close();
        if (merged->TestBit(TFile::kWriteError)) throw std::runtime_error("write failed: " + merged_path);
    } catch (const std::exception &error) { PJ::Fail(error); }
}
