#include <TFile.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <string>
#include <yaml-cpp/yaml.h>
#include <TSystem.h>

#include "plotcommon.h"

void plot_showershapes()
{
    init_plot();
    std::string savePath = "figs/showershapes";
    gSystem->mkdir(savePath.c_str(), true);

    // Load pT and centrality bins from config file
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    YAML::Node config = YAML::LoadFile("../histMakers/nom.yaml");
    std::vector<double> pT_bin_edges = config["analysis"]["pT_bins"].as<std::vector<double>>();
    std::vector<double> centrality_bins = config["analysis"]["centrality_bins"].as<std::vector<double>>();

    int nCentBins = centrality_bins.size() - 1;
    int nPtBins = pT_bin_edges.size() - 1;
    std::cout << "Using " << nPtBins << " pT bins and " << nCentBins << " centrality bins from config file" << std::endl;
    int nCuts = 4; // e.g. 0 or 1 (two cut options)

    //------------------------------------------------------------------------------
    // 2) Open your three ROOT files: data, signal, background
    //    Paths updated to match new histMakers layout (see PlotEvtCharacter.C)
    //------------------------------------------------------------------------------
    TFile *f_data = TFile::Open("../histMakers/results/data_histoshower_shape_data_aa.root", "READ");
    TFile *f_sig  = TFile::Open("../histMakers/results/MC_efficiencyshower_shape_photon10_aa_inclusive.root", "READ");
    TFile *f_bkg  = TFile::Open("../histMakers/results/MC_efficiencyshower_shape_jet10_aa_inclusive.root", "READ");


    std::string leg_sig = "Sig Pythia overlay";
    std::string leg_bkg = "Jet 10 overlay";
    std::string leg_data = "Data AuAu";

    if (!f_data || f_data->IsZombie())
    {
        std::cerr << "Error: Could not open data file!" << std::endl;
        return;
    }
    if (!f_sig || f_sig->IsZombie())
    {
        std::cerr << "Error: Could not open signal file!" << std::endl;
        return;
    }
    if (!f_bkg || f_bkg->IsZombie())
    {
        std::cerr << "Error: Could not open background file!" << std::endl;
        return;
    }

    //------------------------------------------------------------------------------
    // 3) List of base histogram names (the keys in your map).
    //    They should match exactly what you used to create them (e.g. "h2d_prob").
    //------------------------------------------------------------------------------
    std::vector<std::string> histNames = {
        //"h2d_prob",
        //"h2d_CNN_prob",
        "h2d_e17_to_e77",
        //"h2d_e37_to_e77",
        "h2d_e32_to_e35",
        "h2d_e33_to_e35",
        "h2d_e11_to_e33",
        "h2d_e11_to_E",
        "h2d_e33_to_E",
        "h2d_hcalet33_to_ettot",
        "h2d_ihcalet33_to_ettot",
        "h2d_ohcalet33_to_ettot",
        "h2d_hcalet22_to_ettot",
        "h2d_ihcalet22_to_ettot",
        "h2d_ohcalet22_to_ettot",
        //"h2d_detamax",
        //"h2d_dphimax",
        //"h2d_e1",
        //"h2d_e2",
        //"h2d_e3",
        //"h2d_e4",
        "h2d_et1",
        "h2d_et2",
        "h2d_et3",
        "h2d_et4",
        "h2d_weta",
        "h2d_wphi",
        "h2d_w32",
        "h2d_w52",
        "h2d_w72",
        "h2d_wr",
        "h2d_wrr",
        //"h2d_weta_cog",
        //"h2d_wphi_cog",
        "h2d_weta_cogx",
        "h2d_wphi_cogx",
        "h2d_bdt",
        "h2d_npb_score"
    };

    //------------------------------------------------------------------------------
    // Helper lambda to scale a TH1D to unit area
    //------------------------------------------------------------------------------
    auto scaleToUnit = [](TH1D *h)
    {
        if (!h)
            return;
        double integral = h->Integral();
        if (integral > 1e-12)
            h->Scale(1.0 / integral);
    };

    //------------------------------------------------------------------------------
    // 4) Loop over each histogram base name, then over cuts, centrality-bins, and pT-bins
    //------------------------------------------------------------------------------
    for (const auto &hbase : histNames)
    {
        //set

        for (int icut = 0; icut < nCuts; ++icut)
        {
            for (int icent = 0; icent < nCentBins; ++icent)
            {
                for (int ipt = 0; ipt < nPtBins; ++ipt)
                {
                    // x axis name is the base name after first 4 characters
                    TString xaxisname = hbase.substr(4, hbase.size() - 4);

                    float xaxismax = 1.0;
                    float xaxismin = 0;
                    int nrebin = 5;

                    //if the first character of xaxisname is 'w', then the xaxismax is 2.0
                    if (xaxisname[0] == 'w')
                    {
                        xaxismax = 2.0;
                    }
                    //if first is h , then xaxismax is 0.4
                    if (xaxisname[0] == 'h')
                    {
                        xaxismax = 0.3;
                        nrebin = 1;
                    }
                    if (xaxisname.CompareTo("et4") ==0)
                    {
                        xaxismax = 0.25;
                        nrebin = 1;
                    }
                    if (xaxisname.CompareTo("et2") ==0 || xaxisname.CompareTo("et3") ==0)
                    {
                        xaxismax = 1.0;
                        xaxismin = -0.5;
                    }
                    if (xaxisname.CompareTo("e32_to_e35") ==0)
                    {
                        xaxismax = 1;
                        xaxismin = 0.6;
                        nrebin = 1;
                    }
                    if (xaxisname.CompareTo("et1") ==0)
                    {
                        xaxismax = 1;
                        xaxismin = 0.5;
                        nrebin = 2;
                    }

                    // BDT handling
                    if (xaxisname.CompareTo("bdt") == 0)
                    {
                        xaxismax = 1.0;
                        xaxismin = 0.0;
                        nrebin = 2;
                    }

                    TString histNameFull = Form(
                        "%s_cent%d_pt%d_cut%d",
                        hbase.c_str(),
                        icent,
                        ipt,
                        icut);

                    TString histNamesave = Form(
                        "%s_cent%d_pt%d_cut%d",
                        xaxisname.Data(),
                        icent,
                        ipt,
                        icut);

                    TString histNameFull_pp = Form(
                        "%s_cent%d_pt%d_cut%d",
                        hbase.c_str(),
                        0,
                        ipt,
                        icut);

                    TString histNameFull_pp_c = Form(
                        "%s_cent%d_pt%d_cut%d_c",
                        hbase.c_str(),
                        icent,
                        ipt,
                        icut);

                    // Retrieve TH2F from each file
                    TH2F *h2_data = dynamic_cast<TH2F *>(f_data->Get(histNameFull));
                    TH2F *h2_sig = dynamic_cast<TH2F *>(f_sig->Get(histNameFull));
                    TH2F *h2_bkg = dynamic_cast<TH2F *>(f_bkg->Get(histNameFull_pp)->Clone(histNameFull_pp_c));
                    if(!h2_bkg || !h2_sig || !h2_data){
                      cout << "ERROR: cannot find " << histNameFull << endl; 
                      continue;
                    }
                    h2_bkg->RebinX(nrebin);
                    h2_sig->RebinX(nrebin);
                    h2_data->RebinX(nrebin);

                    h2_bkg ->GetXaxis()->SetRangeUser(xaxismin, xaxismax);
                    h2_sig ->GetXaxis()->SetRangeUser(xaxismin, xaxismax);
                    h2_data->GetXaxis()->SetRangeUser(xaxismin, xaxismax);

                    // Check if they exist (skip if not found)
                    if (!h2_data || !h2_sig || !h2_bkg)
                    {
                        // Possibly comment this out if you expect many missing
                        std::cerr << "Warning: Could not retrieve "
                                  << histNameFull << " from one or more files!\n";
                        continue;
                    }

                    TH1D *proj_data = h2_data->ProjectionX(
                        Form("%s_px_data", histNameFull.Data()));
                    TH1D *proj_sig = h2_sig->ProjectionX(
                        Form("%s_px_sig", histNameFull.Data()));
                    TH1D *proj_bkg = h2_bkg->ProjectionX(
                        Form("%s_px_bkg", histNameFull_pp.Data()));

                    scaleToUnit(proj_data);
                    scaleToUnit(proj_sig);
                    scaleToUnit(proj_bkg);

                    float maxy = std::max({proj_data->GetMaximum(), proj_sig->GetMaximum(), proj_bkg->GetMaximum()});

                    TCanvas *c_proj = new TCanvas(
                        Form("c_proj_%s", histNameFull.Data()),
                        Form("ProjectionX - %s", histNameFull.Data()),
                        600, 600);
                    c_proj->cd();

                    // Set colors
                    proj_data->SetLineColor(kBlack);
                    proj_data->SetMarkerColor(kBlack);
                    proj_sig->SetLineColor(kRed);
                    proj_sig->SetMarkerColor(kRed);
                    proj_bkg->SetLineColor(kBlue);
                    proj_bkg->SetMarkerColor(kBlue);

                    // Draw them
                    proj_sig->SetYTitle("normalized counts");
                    proj_sig->GetYaxis()->SetTitleOffset(1.5);
                    proj_sig->SetXTitle(xaxisname.Data());
                    proj_sig->GetYaxis()->SetRangeUser(0, maxy * 1.3);
                    proj_sig->GetXaxis()->SetNdivisions(505);
                    proj_sig->Draw("HIST");

                    proj_bkg->Draw("HIST SAME");

                    proj_data->Draw("ex0 SAME");

                    myText(0.20, 0.90, 1, strleg1.c_str(), 0.04);
                    myText(0.20, 0.85, 1, strleg2.c_str(), 0.04);
                    myText(0.20, 0.80, 1, strleg3.c_str(), 0.04);
                    float pTlow = pT_bin_edges[ipt];
                    float pThigh = pT_bin_edges[ipt + 1];
                    float centlow = centrality_bins[icent];
                    float centhigh = centrality_bins[icent + 1];
                    std::string bgcut;
                    if (icut == 0) bgcut = "w/o nbkg cut";
                    if (icut == 1) bgcut = "w/  nbkg cut";
                    if (icut == 2) bgcut = "w/ tight cut";
                    if (icut == 3) bgcut = "w/ nontight cut";
                    myText(0.2, 0.75, 1, Form("%.0f<p_{T}<%.0f GeV, %.0f-%.0f%% cent, %s", pTlow, pThigh, centlow, centhigh, bgcut.c_str()), 0.04);

                    myMarkerLineText(0.6, 0.90, 1.5, kBlack, 20, kBlack, 1,
                                     leg_data.c_str(), 0.05, true);
                    myMarkerLineText(0.6, 0.85, 0, kRed, 0, kRed, 1,
                                     leg_sig.c_str(), 0.05, true);
                    myMarkerLineText(0.6, 0.80, 0, kBlue, 0, kBlue, 1,
                                     leg_bkg.c_str(), 0.05, true);

                    c_proj->SaveAs(Form("%s/dis_%s.pdf", savePath.c_str(), histNamesave.Data()));


                    
                    TProfile *pfx_data = h2_data->ProfileX(
                        Form("%s_pfx_data", histNameFull.Data()),
                        1,  // first y-bin
                        -1, // last y-bin (use all bins)
                        "s" // 's' option to store RMS in the bin error
                    );

                    TProfile *pfx_sig = h2_sig->ProfileX(
                        Form("%s_pfx_sig", histNameFull.Data()),
                        1, -1, "s");

                    TProfile *pfx_bkg = h2_bkg->ProfileX(
                        Form("%s_pfx_bkg", histNameFull.Data()),
                        1, -1, "");
                    TCanvas *c_prof = new TCanvas(
                        Form("c_prof_%s", histNameFull.Data()),
                        Form("ProfileX - %s", histNameFull.Data()),
                        600, 600);
                    c_prof->cd();

                    // Style
                    if (pfx_data)
                    {
                        pfx_data->SetLineColor(kBlack);
                        pfx_data->SetMarkerColor(kBlack);
                    }
                    if (pfx_sig)
                    {
                        pfx_sig->SetLineColor(kRed);
                        pfx_sig->SetMarkerColor(kRed);
                    }
                    if (pfx_bkg)
                    {
                        pfx_bkg->SetLineColor(kBlue);
                        pfx_bkg->SetMarkerColor(kBlue);
                    }

                    // Draw
                    if (pfx_data)
                        pfx_data->Draw("E");
                    if (pfx_sig)
                        pfx_sig->Draw("E SAME");
                    if (pfx_bkg){
                        pfx_bkg->Draw("hist SAME");
                        pfx_bkg->Draw("ex0 SAME");
                        pfx_bkg->SetMarkerSize(0);
                    }

                    myText(0.50, 0.90, 1, strleg1.c_str(), 0.04);
                    myText(0.50, 0.85, 1, strleg2.c_str(), 0.04);
                    myText(0.50, 0.80, 1, strleg3.c_str(), 0.04);

                    myMarkerLineText(0.55, 0.75, 0, kBlack, 0, kBlack, 1,
                                     leg_data.c_str(), 0.05, true);
                    myMarkerLineText(0.55, 0.70, 0, kRed, 0, kRed, 1,
                                     leg_sig.c_str(), 0.05, true);
                    myMarkerLineText(0.55, 0.65, 0, kBlue, 0, kBlue, 1,
                                     leg_bkg.c_str(), 0.05, true);

                    TCanvas *c_bkg = new TCanvas(
                        Form("c_bkg_%s", histNameFull.Data()),
                        Form("Bkg only - %s", histNameFull.Data()),
                        600, 600);
                    c_bkg->cd();

                    TH1D *proj_bkg_clone = (TH1D *)pfx_bkg->Clone(
                        Form("%s_px_bkgOnly", histNameFull.Data()));
                    proj_bkg_clone->SetLineColor(kBlue);

                    proj_bkg_clone->SetYTitle("<#it{E}_{T}^{iso}> [GeV]");
                    proj_bkg_clone->SetXTitle(xaxisname.Data());
                    float max = proj_bkg_clone->GetMaximum();
                    proj_bkg_clone->GetYaxis()->SetRangeUser(0,max*1.3);

                    proj_bkg_clone->Draw("HIST");
                    proj_bkg_clone->Draw("same ex0");
                    proj_bkg_clone->SetMarkerSize(0);

                    float corr = h2_bkg->GetCorrelationFactor();
                   
                    myText(0.20, 0.90, 1, strleg1.c_str(), 0.04);
                    myText(0.20, 0.85, 1, strleg2.c_str(), 0.04);
                    myText(0.20, 0.80, 1, strleg3.c_str(), 0.04);
                    myText(0.55, 0.90, 1, Form("%.0f<p_{T}<%.0f GeV, %.0f-%.0f%% cent, %s", pTlow, pThigh, centlow, centhigh, bgcut.c_str()), 0.04);
                    myText(0.55, 0.85, 1, "Background MC", 0.04);
                    myText(0.55, 0.80, 1, Form("Correlation: %.3f", corr), 0.04);

                    // Optionally save
                    //c_bkg->SaveAs(Form("%s/pfx_%s.pdf", savePath.c_str(), histNamesave.Data()));

                } // end ipt
            } // end icent
        } // end icut
    } // end loop over histNames

    //------------------------------------------------------------------------------
    // 5) Close files (if you want).
    //------------------------------------------------------------------------------
    f_data->Close();
    f_sig->Close();
    f_bkg->Close();
}
