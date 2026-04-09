#include <yaml-cpp/yaml.h>
#include <TSystem.h>
#include "plotcommon.h"

const int col[] = {kBlack, kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1, kAzure + 7, kOrange + 7, kViolet + 1};
const int mkStyle[] = {20, 21, 22, 33, 34, 24, 25, 27};
const float mkSize[] = {1.2, 1.2, 1.2, 1.1, 1.1, 1.0, 1.0, 1.0};

void plot_efficiency()
{
    init_plot();
    string savePath = "./figs/efficiency";
    gSystem->mkdir(savePath.c_str(), true);

    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    YAML::Node config = YAML::LoadFile("../histMakers/nom.yaml");
    std::vector<double> pT_bins_truth = config["analysis"]["pT_bins_truth"].as<std::vector<double>>();
    std::vector<double> centrality_bins = config["analysis"]["centrality_bins"].as<std::vector<double>>();

    const int nCentBins = static_cast<int>(centrality_bins.size()) - 1;

    TFile *fmc = TFile::Open("../histMakers/results/MC_efficiency_photon20_aa_showershape.root", "READ");
    if (fmc->IsZombie())
    {
        std::cerr << "Could not open MC efficiency file." << std::endl;
        return;
    }

    auto drawEffOverlayByCent = [&](const char *basename, const char *ytitle, const char *outfile) {
        TCanvas *c = new TCanvas(Form("c_%s", basename), Form("c_%s", basename), 650, 600);
        frame_et_truth->SetYTitle(ytitle);
        frame_et_truth->SetXTitle("#it{E}_{T}^{#gamma, truth} [GeV]");
        frame_et_truth->GetXaxis()->SetRangeUser(pT_bins_truth.front(), pT_bins_truth.back());
        frame_et_truth->GetYaxis()->SetRangeUser(0.0, 1.4);
        frame_et_truth->Draw("axis");

        TLegend *leg = new TLegend(0.50, 0.20, 0.90, 0.20 + std::min(0.06 * nCentBins, 0.52));
        legStyle(leg, 0.20, 0.04);

        int nDrawn = 0;
        for (int icent = 0; icent < nCentBins; ++icent)
        {
            TEfficiency *eff = (TEfficiency *)fmc->Get(Form("%s%d", basename, icent));
            if (!eff)
            {
                std::cerr << "Missing TEfficiency: " << Form("%s%d", basename, icent) << std::endl;
                continue;
            }

            const int istyle = icent % 8;
            eff->SetMarkerColor(col[istyle]);
            eff->SetMarkerStyle(mkStyle[istyle]);
            eff->SetMarkerSize(mkSize[istyle]);
            eff->SetLineColor(col[istyle]);
            eff->SetLineWidth(2);
            eff->Draw("same");


            leg->AddEntry(eff, Form("%.0f-%.0f%%", centrality_bins[icent], centrality_bins[icent + 1]), "pl");
            ++nDrawn;
        }

        myText(0.20, 0.90, 1, strleg1.c_str(), 0.048);
        myText(0.20, 0.85, 1, strleg2.c_str(), 0.046);
        myText(0.20, 0.80, 1, strMC.c_str(), 0.046);
        myText(0.20, 0.75, 1, strleg3.c_str(), 0.046);

        LINE1->Draw("l");
        
        leg->Draw("same");

        c->SaveAs(Form("%s/%s", savePath.c_str(), outfile));
    };

    drawEffOverlayByCent("eff_reco_cent", "Reconstruction Efficiency", "eff_reco.pdf");
    drawEffOverlayByCent("eff_iso_cent", "Isolation Efficiency", "eff_iso.pdf");
    drawEffOverlayByCent("eff_id_cent", "Identification Efficiency", "eff_id.pdf");
    drawEffOverlayByCent("eff_all_cent", "Total Efficiency", "eff_total.pdf");
    drawEffOverlayByCent("eff_converts_cent", "Conversion Probability", "eff_converts.pdf");


}
