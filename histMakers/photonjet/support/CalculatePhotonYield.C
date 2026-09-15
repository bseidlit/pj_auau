#include <RooUnfoldResponse.h>
#include <RooUnfoldBayes.h>
#include <yaml-cpp/yaml.h>
#include <cmath>
#include "CrossSectionWeights.h"
using namespace PPG12;

double myfunc(double *x, double *params)
{
    double NsigA = x[0]; // Independent variable (N_A^sig)

    // Parameters passed to the function
    double NA = params[0];
    double NB = params[1];
    double NC = params[2];
    double ND = params[3];
    double cB = params[4];
    double cC = params[5];
    double cD = params[6];
    double R = params[7];

    // The equation:
    double numerator = NB - cB * NsigA;
    double denominator = ND - cD * NsigA;

    // Avoid division by zero
    if (denominator == 0)
    {
        denominator = 1e-12;
    }

    double fraction = (NC - cC * NsigA) / denominator;

    // The function to find the root of
    double result = NsigA - (NA - R * numerator * fraction);
    return result;
}
void scale_histogram(TH1 *h, float lumi)
{
    for (int ibin = 1; ibin <= h->GetNbinsX(); ibin++)
    {
        float binwidth = h->GetBinWidth(ibin);
        float bincenter = h->GetBinCenter(ibin);

        float scale = 1.0 / binwidth / lumi;

        h->SetBinContent(ibin, h->GetBinContent(ibin) * scale);
        h->SetBinError(ibin, h->GetBinError(ibin) * scale);
    }
}

void CalculatePhotonYield(const std::string &configname = "config_bdt_purity_pade.yaml", bool isMC = false)
{
    // isMC=true runs the MC closure: the inclusive jet MC (which already contains
    // the prompt-photon processes at the PYTHIA rate) takes the place of the data,
    // while the photon MC provides the signal leakage and efficiencies as in the data run.
    float mbdcorr = 25.2/42 / 0.57;
    float solid_angle = 2 * M_PI * 0.7 * 2;
    // lumi times cross section is events
    float nsimevents = 1E7;
    // float nsimevents = 2417664.0;
    //const float photon20cross = 1.571e+05 * 0.000673448;
    float simluminosity = nsimevents / photon20cross;

    // [PPG12] Herwig cross-check: photon20_herwig has 2.787M generated events,
    // and the reference cross-section is herwig_photon20cross. Detected by
    // var_type substring; loaded after configYaml below.

    // float jetevents = 2194879.0;
    float jetevents = 0.3555 * 1E7;
    const float jetcross = jet50cross;

    float jetluminosity = jetevents / jetcross;
    bool fit_purity = true;
    bool fit_purity_dis = true;

    gSystem->Load("libyaml-cpp");
    YAML::Node configYaml = YAML::LoadFile(configname);

    // [PPG12] Herwig cross-check: signal MC is photon{10,20}_herwig, with the
    // reference cross-section herwig_photon20cross and N_events=2.787M (run-28
    // type 47). Override the Pythia-tuned nsimevents/simluminosity here so the
    // downstream MC luminosity normalization is correct. Pythia path unchanged.
    {
        std::string vt_check = configYaml["output"]["var_type"].as<std::string>();
        if (vt_check.find("herwig") != std::string::npos) {
            nsimevents = 2.787e+06f;
            simluminosity = nsimevents / herwig_photon20cross;
            std::cout << "[CalculatePhotonYield] HERWIG cross-check (var_type=" << vt_check
                      << "): nsimevents=" << nsimevents
                      << " simluminosity=" << simluminosity << " pb-1" << std::endl;
        }
    }

    // luminosity: read from config (pb^-1), default 49.562
    float luminosity = configYaml["analysis"]["lumi"].as<float>(49.562);

    // for mc we can check the actual purity
    if (isMC)
    {
        // the inclusive jet MC is the "data", so use its luminosity
        luminosity = jetluminosity;
    }

    int fittingerror = configYaml["analysis"]["fittingerror"].as<int>(0);
    bool mc_purity_correction = configYaml["analysis"]["mc_purity_correction"].as<int>(0) != 0;
    std::string mc_purity_correction_file = configYaml["analysis"]["mc_purity_correction_file"].as<std::string>("");
    // Polynomial order of the smooth fit to the MC purity nonclosure ratio.
    // Default 1 (pol1, linear); set to 2 to recover the legacy pol2 smoothing.
    int mc_purity_corr_fitorder = configYaml["analysis"]["mc_purity_corr_fitorder"].as<int>(1);

    std::string var_type = configYaml["output"]["var_type"].as<std::string>();

    std::string mcstring = isMC ? "_mc" : "";

    std::string outfilename = configYaml["output"]["final_outfile"].as<std::string>() + "_" + var_type + mcstring + ".root";
    std::string mc_outfilename = configYaml["output"]["final_outfile"].as<std::string>() + "_" + var_type + "_mc.root";

    // calculate and unfold the isolated photon spectrum given result from RecoEffCalculator.C

    std::string datainput = configYaml["output"]["data_outfile"].as<std::string>() + "_" + var_type + ".root";
    if (isMC)
    {
        datainput = configYaml["output"]["eff_outfile"].as<std::string>() + "_jet_" + var_type + ".root";
    }
    std::string siminput = configYaml["output"]["eff_outfile"].as<std::string>() + "_" + var_type + ".root";

    std::string unfoldinput = configYaml["output"]["response_outfile"].as<std::string>() + "_" + var_type + ".root";

    TFile *f_mc_purity_corr = nullptr;
    TGraphErrors *g_mc_purity_fit_ratio = nullptr;
    TF1 *f_mc_purity_corr_fit = nullptr;  // smooth fit to bin-by-bin truth/ABCD ratio
    if (!isMC && mc_purity_correction)
    {
        std::string mc_ratio_source = mc_purity_correction_file.empty() ? mc_outfilename : mc_purity_correction_file;
        f_mc_purity_corr = new TFile(mc_ratio_source.c_str(), "READ");
        if (!f_mc_purity_corr || f_mc_purity_corr->IsZombie())
        {
            std::cout << "WARNING: cannot open MC purity correction file " << mc_ratio_source
                      << ". Skip MC purity correction." << std::endl;
            mc_purity_correction = false;
        }
        else
        {
            g_mc_purity_fit_ratio = (TGraphErrors *)f_mc_purity_corr->Get("g_mc_purity_fit_ratio");
            f_mc_purity_corr_fit = (TF1 *)f_mc_purity_corr->Get("f_mc_purity_corr_fit");
            if (!g_mc_purity_fit_ratio || !f_mc_purity_corr_fit)
            {
                std::cout << "WARNING: g_mc_purity_fit_ratio or f_mc_purity_corr_fit not found in " << mc_ratio_source
                          << ". Skip MC purity correction." << std::endl;
                mc_purity_correction = false;
            }
            else
            {
                std::cout << "Apply MC purity correction (fit-of-bin-by-bin truth/ABCD ratio) using " << mc_ratio_source << std::endl;
            }
        }
    }

    int fitoption = configYaml["analysis"]["fit_option"].as<int>(0);

    float mbd_eff_scale = configYaml["analysis"]["mbd_eff_scale"].as<float>(0.0);

    std::string histogram_postfix = "_0";

    std::string tight_iso_cluster_name = "h_tight_iso_cluster"; // A

    std::string tight_noniso_cluster_name = "h_tight_noniso_cluster"; // B

    std::string nontight_iso_cluster_name = "h_nontight_iso_cluster"; // C

    std::string nontight_noniso_cluster_name = "h_nontight_noniso_cluster"; // D

    std::string common_cluster_name = "h_common_cluster";

    std::string tight_iso_cluster_signal_name = tight_iso_cluster_name + "_signal" + histogram_postfix; // A

    std::string tight_noniso_cluster_signal_name = tight_noniso_cluster_name + "_signal" + histogram_postfix; // B

    std::string nontight_iso_cluster_signal_name = nontight_iso_cluster_name + "_signal" + histogram_postfix; // C

    std::string nontight_noniso_cluster_signal_name = nontight_noniso_cluster_name + "_signal" + histogram_postfix; // D

    std::string tight_iso_cluster_notmatch_name     = tight_iso_cluster_name     + "_notmatch" + histogram_postfix; // A
    std::string tight_noniso_cluster_notmatch_name  = tight_noniso_cluster_name  + "_notmatch" + histogram_postfix; // B
    std::string nontight_iso_cluster_notmatch_name  = nontight_iso_cluster_name  + "_notmatch" + histogram_postfix; // C
    std::string nontight_noniso_cluster_notmatch_name = nontight_noniso_cluster_name + "_notmatch" + histogram_postfix; // D

    std::string reco_efficiency_name = "eff_reco_eta" + histogram_postfix;

    std::string reco_efficiency_iso_name = "eff_iso_eta" + histogram_postfix;

    std::string reco_efficiency_id_name = "eff_id_eta" + histogram_postfix;

    std::string truth_pythia_name = "h_truth_pT" + histogram_postfix;
    std::string truth_pythia_novtx_name = "h_truth_pT_novtx" + histogram_postfix;

    std::string truth_with_vertex_name = "h_truth_pT_vertexcut" + histogram_postfix;

    std::string truth_with_vertex_mbd_name = "h_truth_pT_vertexcut_mbd_cut" + histogram_postfix;
    std::string truth_with_vertex_mbd_north_name = "h_truth_pT_vertexcut_mbd_north_cut" + histogram_postfix;
    std::string truth_with_vertex_mbd_south_name = "h_truth_pT_vertexcut_mbd_south_cut" + histogram_postfix;

    // calculate the signal leakage from simulation

    TFile *fsimin = new TFile(siminput.c_str(), "READ");

    TH1F *h_tight_iso_cluster_signal = (TH1F *)fsimin->Get(tight_iso_cluster_signal_name.c_str());

    TH1F *h_tight_noniso_cluster_signal = (TH1F *)fsimin->Get(tight_noniso_cluster_signal_name.c_str());

    TH1F *h_nontight_iso_cluster_signal = (TH1F *)fsimin->Get(nontight_iso_cluster_signal_name.c_str());

    TH1F *h_nontight_noniso_cluster_signal = (TH1F *)fsimin->Get(nontight_noniso_cluster_signal_name.c_str());

    TH1F *h_tight_iso_cluster_notmatch     = (TH1F *)fsimin->Get(tight_iso_cluster_notmatch_name.c_str());
    TH1F *h_tight_noniso_cluster_notmatch  = (TH1F *)fsimin->Get(tight_noniso_cluster_notmatch_name.c_str());
    TH1F *h_nontight_iso_cluster_notmatch  = (TH1F *)fsimin->Get(nontight_iso_cluster_notmatch_name.c_str());
    TH1F *h_nontight_noniso_cluster_notmatch = (TH1F *)fsimin->Get(nontight_noniso_cluster_notmatch_name.c_str());

    TH1F *h_truth_pT_vertexcut = (TH1F *)fsimin->Get(truth_with_vertex_name.c_str());

    TH1F *h_truth_pT_vertexcut_mbd_cut = (TH1F *)fsimin->Get(truth_with_vertex_mbd_name.c_str());
    TH1F *h_truth_pT_vertexcut_mbd_north_cut = (TH1F *)fsimin->Get(truth_with_vertex_mbd_north_name.c_str());
    TH1F *h_truth_pT_vertexcut_mbd_south_cut = (TH1F *)fsimin->Get(truth_with_vertex_mbd_south_name.c_str());

    TH1F *h_leak_B = (TH1F *)h_tight_noniso_cluster_signal->Clone("h_leak_B");
    TH1F *h_leak_C = (TH1F *)h_nontight_iso_cluster_signal->Clone("h_leak_C");
    TH1F *h_leak_D = (TH1F *)h_nontight_noniso_cluster_signal->Clone("h_leak_D");

    TEfficiency *eff_reco = (TEfficiency *)fsimin->Get(reco_efficiency_name.c_str());
    TEfficiency *eff_iso = (TEfficiency *)fsimin->Get(reco_efficiency_iso_name.c_str());
    TEfficiency *eff_id = (TEfficiency *)fsimin->Get(reco_efficiency_id_name.c_str());

    TH1F *h_pythia_truth = (TH1F *)fsimin->Get(truth_pythia_name.c_str());
    TH1F *h_pythia_truth_novtx = (TH1F *)fsimin->Get(truth_pythia_novtx_name.c_str());

    // TH1F* h_truth_with_cut = (TH1F *)fsimin->Get(truth_with_cut_name.c_str());

    // h_tight_iso_cluster_signal->Sumw2();
    // h_leak_B->Sumw2();
    // h_leak_C->Sumw2();
    // h_leak_D->Sumw2();

    h_leak_B->Divide(h_tight_iso_cluster_signal);
    h_leak_C->Divide(h_tight_iso_cluster_signal);
    h_leak_D->Divide(h_tight_iso_cluster_signal);
    // calculate error manually
    /*
    for (int ibin = 1; ibin <= h_leak_B->GetNbinsX(); ibin++)
    {
        double A = h_tight_iso_cluster_signal->GetBinContent(ibin);
        double B = h_tight_noniso_cluster_signal->GetBinContent(ibin);
        double C = h_nontight_iso_cluster_signal->GetBinContent(ibin);
        double D = h_nontight_noniso_cluster_signal->GetBinContent(ibin);


        double errA = h_tight_iso_cluster_signal->GetBinError(ibin);
        double errB = h_tight_noniso_cluster_signal->GetBinError(ibin);
        double errC = h_nontight_iso_cluster_signal->GetBinError(ibin);
        double errD = h_nontight_noniso_cluster_signal->GetBinError(ibin);

        double errBA = h_leak_B->GetBinContent(ibin) * sqrt(pow(errB / B, 2) + pow(errA / A, 2));
        double errCA = h_leak_C->GetBinContent(ibin) * sqrt(pow(errC / C, 2) + pow(errA / A, 2));
        double errDA = h_leak_D->GetBinContent(ibin) * sqrt(pow(errD / D, 2) + pow(errA / A, 2));
        std::cout << "errBA: " << errBA << " errCA: " << errCA << " errDA: " << errDA << std::endl;
        //nan gard
        if (isnan(errBA))
        {
            errBA = 0;
        }
        if (isnan(errCA))
        {
            errCA = 0;
        }
        if (isnan(errDA))
        {
            errDA = 0;
        }
        h_leak_B->SetBinError(ibin, errBA);
        h_leak_C->SetBinError(ibin, errCA);
        h_leak_D->SetBinError(ibin, errDA);
    }
    */

    // get the raw spectrum from data

    std::string tight_iso_cluster_signal_name_data = tight_iso_cluster_name + histogram_postfix; // A
    // std::cout<<"tight_iso_cluster_signal_name_data: "<<tight_iso_cluster_signal_name_data<<std::endl;

    std::string tight_noniso_cluster_signal_name_data = tight_noniso_cluster_name + histogram_postfix; // B

    std::string nontight_iso_cluster_signal_name_data = nontight_iso_cluster_name + histogram_postfix; // C

    std::string nontight_noniso_cluster_signal_name_data = nontight_noniso_cluster_name + histogram_postfix; // D

    std::string common_cluster_name_data = common_cluster_name + histogram_postfix;

    TFile *fdatain = new TFile(datainput.c_str(), "READ");

    // Defensive guard (F2): catch corrupted/recovered files from transient
    // gpfs I/O during heavy condor batches. ROOT auto-recovers truncated
    // files and silently leaves TH1 streamer state inconsistent — downstream
    // Clone()/Add() then segfaults. Abort cleanly with an explicit message.
    if (!fdatain || fdatain->IsZombie())
    {
        std::cerr << "[FATAL] CalculatePhotonYield: cannot open data input file "
                  << datainput << " (zombie). Aborting." << std::endl;
        return;
    }
    if (fdatain->TestBit(TFile::kRecovered))
    {
        std::cerr << "[FATAL] CalculatePhotonYield: data input file " << datainput
                  << " was auto-recovered by ROOT (likely corrupted/truncated). "
                  << "Streamer state of TH1 objects is unreliable. Aborting." << std::endl;
        fdatain->Close();
        return;
    }

    TH1F *h_tight_iso_cluster_signal_data = (TH1F *)fdatain->Get(tight_iso_cluster_signal_name_data.c_str());

    TH1F *h_tight_noniso_cluster_signal_data = (TH1F *)fdatain->Get(tight_noniso_cluster_signal_name_data.c_str());

    TH1F *h_nontight_iso_cluster_signal_data = (TH1F *)fdatain->Get(nontight_iso_cluster_signal_name_data.c_str());

    TH1F *h_nontight_noniso_cluster_signal_data = (TH1F *)fdatain->Get(nontight_noniso_cluster_signal_name_data.c_str());

    if (!h_tight_iso_cluster_signal_data || !h_tight_noniso_cluster_signal_data ||
        !h_nontight_iso_cluster_signal_data || !h_nontight_noniso_cluster_signal_data)
    {
        std::cerr << "[FATAL] CalculatePhotonYield: missing ABCD histogram(s) in "
                  << datainput << ". Aborting." << std::endl;
        fdatain->Close();
        return;
    }

    TH1F *h_tight_iso_cluster_background_data = (TH1F *) h_tight_iso_cluster_signal_data->Clone("h_tight_iso_cluster_background_data");

    TH1F *h_tight_noniso_cluster_background_data = (TH1F *) h_tight_noniso_cluster_signal_data->Clone("h_tight_noniso_cluster_background_data");

    TH1F *h_nontight_iso_cluster_background_data = (TH1F *) h_nontight_iso_cluster_signal_data->Clone("h_nontight_iso_cluster_background_data");

    TH1F *h_nontight_noniso_cluster_background_data = (TH1F *) h_nontight_noniso_cluster_signal_data->Clone("h_nontight_noniso_cluster_background_data");

    TH1F *h_tight_iso_cluster_signal_inclusive = (TH1F *)fdatain->Get(tight_iso_cluster_signal_name.c_str());

    TH1F *h_tight_noniso_cluster_signal_inclusive = (TH1F *)fdatain->Get(tight_noniso_cluster_signal_name.c_str());

    TH1F *h_nontight_iso_cluster_signal_inclusive = (TH1F *)fdatain->Get(nontight_iso_cluster_signal_name.c_str());

    TH1F *h_nontight_noniso_cluster_signal_inclusive = (TH1F *)fdatain->Get(nontight_noniso_cluster_signal_name.c_str());

    h_tight_iso_cluster_background_data->Add(h_tight_iso_cluster_signal_inclusive, -1);
    h_tight_noniso_cluster_background_data->Add(h_tight_noniso_cluster_signal_inclusive, -1);
    h_nontight_iso_cluster_background_data->Add(h_nontight_iso_cluster_signal_inclusive, -1);
    h_nontight_noniso_cluster_background_data->Add(h_nontight_noniso_cluster_signal_inclusive, -1);

    // R = A/B*D/C  (data-driven: background = data - signal leak)
    TH1F *h_R = (TH1F *)h_tight_iso_cluster_background_data->Clone("h_R");
    h_R->Divide(h_tight_noniso_cluster_background_data);
    h_R->Multiply(h_nontight_noniso_cluster_background_data);
    h_R->Divide(h_nontight_iso_cluster_background_data);

    // R = A/B*D/C  (MC truth: background = non-truth-matched clusters directly from simulation)
    TH1F *h_R_notmatch = (TH1F *)h_tight_iso_cluster_notmatch->Clone("h_R_notmatch");
    h_R_notmatch->Divide(h_tight_noniso_cluster_notmatch);
    h_R_notmatch->Multiply(h_nontight_noniso_cluster_notmatch);
    h_R_notmatch->Divide(h_nontight_iso_cluster_notmatch);
    

    TH1F *h_common_cluster_data = (TH1F *)fdatain->Get(common_cluster_name_data.c_str());

    TH1F *h_truth_iso_cluster_data = (TH1F *)fdatain->Get(tight_iso_cluster_signal_name.c_str());

    // a histogram with same binning as h_tight_iso_cluster_signal_data
    // TH1F *h_purity = (TH1F *)h_tight_iso_cluster_signal_data->Clone("h_purity");
    // h_purity->Reset();
    TH1F *h_data_sub = (TH1F *)h_tight_iso_cluster_signal_data->Clone("h_data_sub");
    h_data_sub->Reset();
    TGraphErrors *gpurity = new TGraphErrors(h_tight_iso_cluster_signal_data);
    gpurity->SetName("gpurity");

    // purity histogram with leakage correction
    // TH1F *h_purity_leak = (TH1F *)h_tight_iso_cluster_signal_data->Clone("h_purity_leak");
    // h_purity_leak->Reset();
    TH1F *h_data_sub_leak = (TH1F *)h_tight_iso_cluster_signal_data->Clone("h_data_sub_leak");
    h_data_sub_leak->Reset();
    TGraphErrors *gpurity_leak = new TGraphErrors(h_tight_iso_cluster_signal_data);
    gpurity_leak->SetName("gpurity_leak");

    TGraphAsymmErrors *g_purity_truth = new TGraphAsymmErrors(h_truth_iso_cluster_data, h_tight_iso_cluster_signal_data);
    g_purity_truth->SetName("g_purity_truth");
    TGraphErrors *g_mc_purity_fit_ratio_out = nullptr;
    TF1 *f_purity_truth_fit = nullptr;

    TGraphAsymmErrors *g_mbd_eff = new TGraphAsymmErrors(h_truth_pT_vertexcut_mbd_cut, h_truth_pT_vertexcut);
    g_mbd_eff->SetName("g_mbd_eff");
    TGraphAsymmErrors *g_mbd_eff_north = new TGraphAsymmErrors(h_truth_pT_vertexcut_mbd_north_cut, h_truth_pT_vertexcut);
    g_mbd_eff_north->SetName("g_mbd_eff_north");
    TGraphAsymmErrors *g_mbd_eff_south = new TGraphAsymmErrors(h_truth_pT_vertexcut_mbd_south_cut, h_truth_pT_vertexcut);
    g_mbd_eff_south->SetName("g_mbd_eff_south");

    // calculate the purity
    //std::vector<TH1F *> h_NA_sig_list;
    std::vector<TH1F *> h_purity_list;
    std::vector<TH1F *> h_purity_leak_list;
    std::vector<TF1 *> f_purity_list;
    std::vector<TF1 *> f_purity_leak_list;
    // Fixed seed (42) so the toy-throw step in the ABCD purity uncertainty
    // (Poisson re-throws for raw yields + Gaussian re-throws for the
    // signal-leakage corrections cB/cC/cD) is bit-reproducible across
    // pipeline runs. Seed 0 (CPU clock) had been used previously, which
    // made the per-bin ±2-sigma fluctuation envelope drift between
    // re-runs and contaminated the systematic comparison vs nominal.
    // Matches the global analysis seed used elsewhere in the pipeline.
    TRandom3 *rand = new TRandom3(42);
    h_tight_iso_cluster_signal_data->Sumw2();
    h_tight_noniso_cluster_signal_data->Sumw2();
    h_nontight_iso_cluster_signal_data->Sumw2();
    h_nontight_noniso_cluster_signal_data->Sumw2();
    for (int ibin = 1; ibin <= h_tight_iso_cluster_signal_data->GetNbinsX(); ibin++)
    {
        std::cout << "ibin: " << ibin << "bin center: " << h_tight_iso_cluster_signal_data->GetBinCenter(ibin) << std::endl;
        int nsamples = 20000;

        double A = h_tight_iso_cluster_signal_data->GetBinContent(ibin);
        double B = h_tight_noniso_cluster_signal_data->GetBinContent(ibin);
        double C = h_nontight_iso_cluster_signal_data->GetBinContent(ibin);
        double D = h_nontight_noniso_cluster_signal_data->GetBinContent(ibin);

        double A_w2 = h_tight_iso_cluster_signal_data->GetSumw2()->At(ibin);
        double B_w2 = h_tight_noniso_cluster_signal_data->GetSumw2()->At(ibin);
        double C_w2 = h_nontight_iso_cluster_signal_data->GetSumw2()->At(ibin);
        double D_w2 = h_nontight_noniso_cluster_signal_data->GetSumw2()->At(ibin);

        //std::cout<< "A: " << A << " B: " << B << " C: " << C << " D: " << D << " A_w2: " << A_w2 << " B_w2: " << B_w2 << " C_w2: " << C_w2 << " D_w2: " << D_w2 << std::endl;

        double A_eff = (A != 0) ? (A * A) / A_w2 : 0;
        double B_eff = (B != 0) ? (B * B) / B_w2 : 0;
        double C_eff = (C != 0) ? (C * C) / C_w2 : 0;
        double D_eff = (D != 0) ? (D * D) / D_w2 : 0;

        double errA = h_tight_iso_cluster_signal_data->GetBinError(ibin);
        double errB = h_tight_noniso_cluster_signal_data->GetBinError(ibin);
        double errC = h_nontight_iso_cluster_signal_data->GetBinError(ibin);
        double errD = h_nontight_noniso_cluster_signal_data->GetBinError(ibin);

        double cB = h_leak_B->GetBinContent(ibin);
        // double cC = h_leak_C->GetBinContent(ibin);
        double cC = h_leak_C->GetBinContent(ibin);
        double cD = h_leak_D->GetBinContent(ibin);

        double errcb = h_leak_B->GetBinError(ibin);
        double errcc = h_leak_C->GetBinError(ibin);
        double errcd = h_leak_D->GetBinError(ibin);

        TH1D *h_result = new TH1D("h_result", "h_result", 1000, -A, A * 10);
        TH1D *h_result_purity = new TH1D("h_result_purity", "h_result_purity", 1000, -1, 2);
        
        TH1D *h_result_leak = new TH1D("h_result_leak", "h_result_leak", 1000, -A, A * 10);
        TH1D *h_result_purity_leak = new TH1D("h_result_purity_leak", "h_result_purity_leak", 1000, -1, 2);

        for (int isample = 0; isample < nsamples; isample++)
        {
            
            //use poisson to sample the data
            double NA_eff = rand->PoissonD(A_eff);
            double NB_eff = rand->PoissonD(B_eff);
            double NC_eff = rand->PoissonD(C_eff);
            double ND_eff = rand->PoissonD(D_eff);

            double NA = A * NA_eff / A_eff;
            double NB = B * NB_eff / B_eff;
            double NC = C * NC_eff / C_eff;
            double ND = D * ND_eff / D_eff;

            //std::cout << "NA: " << NA << " NB: " << NB << " NC: " << NC << " ND: " << ND << std::endl;

            /*
            double NA = rand->Gaus(A, 0);
            double NB = rand->Gaus(B, 0);
            double NC = rand->Gaus(C, 0);
            double ND = rand->Gaus(D, 0);
            */
            double CB = 0;
            double CC = 0;
            double CD = 0;

            TF1 *f = new TF1("myfunc", myfunc, 0.0, A, 8);
            f->SetParameter(0, NA);
            f->SetParameter(1, NB);
            f->SetParameter(2, NC);
            f->SetParameter(3, ND);
            f->SetParameter(4, CB);
            f->SetParameter(5, CC);
            f->SetParameter(6, CD);
            f->SetParameter(7, 1);

            double root = f->GetX(0, -0.5*NA, 2*NA);

            double root_manual = NA - NB * NC / ND;

            // fill if the root is not nan nor inf
            if (root == root && abs(root) != std::numeric_limits<double>::infinity())
            {
                h_result->Fill(root);
                h_result_purity->Fill(root / NA);
            }

            // std::cout<< "NA: " << NA << " NB: " << NB << " NC: " << NC << " ND: " << ND << std::endl;
            // std::cout << "root: " << root << " NA: " << NA << std::endl;
            // std::cout << "root_manual: " << root_manual << " NA: " << NA << std::endl;

            CB = rand->Gaus(cB, errcb);
            CC = rand->Gaus(cC, errcc);
            CD = rand->Gaus(cD, errcd);

            f->SetParameter(4, CB);
            f->SetParameter(5, CC);
            f->SetParameter(6, CD);

            double root_leak = f->GetX(0, -0.5*NA, 2*NA);

            if (root_leak == root_leak && abs(root_leak) != std::numeric_limits<double>::infinity())
            {
                h_result_leak->Fill(root_leak);
                h_result_purity_leak->Fill(root_leak / NA);
            }

            // std::cout << "root_leak: " << root_leak << " NA: " << NA << std::endl;
        }
        //get the nomilal centered value
        double NA = A;
        double NB = B;
        double NC = C;
        double ND = D;

        double CB = cB;
        double CC = cC;
        double CD = cD;
        //no leakage correction
        TF1 *f = new TF1("myfunc", myfunc, 0.0, A, 8);
        f->SetParameter(0, NA);
        f->SetParameter(1, NB);
        f->SetParameter(2, NC);
        f->SetParameter(3, ND);
        f->SetParameter(4, 0);
        f->SetParameter(5, 0);
        f->SetParameter(6, 0);
        f->SetParameter(7, 1);

        double root = f->GetX(0, 0, 2*NA);

        float purity_result = 0;
        if (root == root && abs(root) != std::numeric_limits<double>::infinity())
        {
            purity_result = root / NA;
        }

        //leakage correction
        f->SetParameter(4, CB);
        f->SetParameter(5, CC);
        f->SetParameter(6, CD);

        double root_leak = f->GetX(0, 0, 2*NA);

        float purity_result_leak = 0;
        if (root_leak == root_leak && abs(root_leak) != std::numeric_limits<double>::infinity())
        {
            purity_result_leak = root_leak / NA;
        }

        std::cout << "NA: " << A << " NB: " << B << " NC: " << C << " ND: " << D << std::endl;
        std::cout << "result mean: " << h_result->GetMean() << " result std: " << h_result->GetRMS() << std::endl;
        std::cout << "result leak mean: " << h_result_leak->GetMean() << " result leak std: " << h_result_leak->GetRMS() << std::endl;

        h_data_sub->SetBinContent(ibin, h_result->GetMean());
        h_data_sub->SetBinError(ibin, h_result->GetRMS());

        h_data_sub_leak->SetBinContent(ibin, h_result_leak->GetMean());
        h_data_sub_leak->SetBinError(ibin, h_result_leak->GetRMS());
        
        TH1F *h_purity = (TH1F *)h_result_purity->Clone(Form("h_purity_%d", ibin));
        h_purity->SetName(Form("h_purity_%d", ibin));


        TH1F *h_purity_leak = (TH1F *)h_result_purity_leak->Clone(Form("h_purity_leak_%d", ibin));
        h_purity_leak->SetName(Form("h_purity_leak_%d", ibin));
   
        // calculate the purity
        float purity = h_result_purity->GetMean();
        float purity_err = h_result_purity->GetRMS();

        float purity_leak = h_result_purity_leak->GetMean();
        float purity_err_leak = h_result_purity_leak->GetRMS();



        //gaussian fit the purity
        if (fit_purity)
        {
            float fitlower = h_result_purity->GetMean() - 1 * h_result_purity->GetRMS();
            float fitupper = h_result_purity->GetMean() + 1.5 * h_result_purity->GetRMS();
            TF1 *f_purity = new TF1(Form("f_purity_%d", ibin), "gaus", fitlower, fitupper);
            h_result_purity->Fit(f_purity, "REMQN", "", fitlower, fitupper);
            f_purity_list.push_back(f_purity);
            purity = f_purity->GetParameter(1);
            purity_err = f_purity->GetParameter(2);

            float fitlower_leak = h_result_purity_leak->GetMean() - 1 * h_result_purity_leak->GetRMS();
            float fitupper_leak = h_result_purity_leak->GetMean() + 1.5 * h_result_purity_leak->GetRMS();
            TF1 *f_purity_leak = new TF1(Form("f_purity_leak_%d", ibin), "gaus", fitlower_leak, fitupper_leak);
            h_result_purity_leak->Fit(f_purity_leak, "REMQN", "", fitlower_leak, fitupper_leak);
            f_purity_leak_list.push_back(f_purity_leak);
            purity_leak = f_purity_leak->GetParameter(1);
            purity_err_leak = f_purity_leak->GetParameter(2);
        }

        h_purity_leak_list.push_back(h_purity_leak);
        h_purity_list.push_back(h_purity);
        
        //gpurity->SetPoint(ibin - 1, h_tight_iso_cluster_signal_data->GetBinCenter(ibin), purity_result);
        gpurity->SetPoint(ibin - 1, h_tight_iso_cluster_signal_data->GetBinCenter(ibin), purity);
        gpurity->SetPointError(ibin - 1, h_tight_iso_cluster_signal_data->GetBinWidth(ibin) / 2, purity_err);

        if (h_result->GetMean() == 0)
        {
            gpurity->SetPoint(ibin - 1, h_tight_iso_cluster_signal_data->GetBinCenter(ibin), 0);
            gpurity->SetPointError(ibin - 1, h_tight_iso_cluster_signal_data->GetBinWidth(ibin) / 2, 0);
        }

        //gpurity_leak->SetPoint(ibin - 1, h_tight_iso_cluster_signal_data->GetBinCenter(ibin), purity_result_leak);
        gpurity_leak->SetPoint(ibin - 1, h_tight_iso_cluster_signal_data->GetBinCenter(ibin), purity_leak);
        gpurity_leak->SetPointError(ibin - 1, h_tight_iso_cluster_signal_data->GetBinWidth(ibin) / 2, purity_err_leak);

        if (h_result_leak->GetMean() == 0)
        {
            gpurity_leak->SetPoint(ibin - 1, h_tight_iso_cluster_signal_data->GetBinCenter(ibin), 0);
            gpurity_leak->SetPointError(ibin - 1, h_tight_iso_cluster_signal_data->GetBinWidth(ibin) / 2, 0);
        }
    }



    int nFinePoints = 1000; // Number of points for granular intervals
    double xMin = 10;       // Fit range start (tightened 2026-04-24; was 8)
    double xMax = 36;       // Fit range end (extended 2026-04-27; was 32)

    // Padé init seeded from the nom-converged values so variant fits don't
    // fall into the runaway minimum where (p0+p1*x) and (1+p2^2*x) both blow
    // up but their ratio mimics the data. Per-variant gpurity values are
    // ~1% from nom, so starting near the nom minimum keeps the optimizer
    // there. Bare TGraphErrors->Fit() previously diverged for npb03 (chi2
    // 224 vs 5.7 nominal, params at 1e+27, fit eval 1.25 at low pT
    // when the per-bin gpurity value was 0.52). See session 2026-05-03.
    constexpr double pade_purity_p0      =  0.129795; // f_purity_fit nom
    constexpr double pade_purity_p1      =  0.032936;
    constexpr double pade_purity_p2      =  0.151076;
    constexpr double pade_purity_leak_p0 =  0.085956; // f_purity_leak_fit nom
    constexpr double pade_purity_leak_p1 =  0.052318;
    constexpr double pade_purity_leak_p2 = -0.164694;

    TF1 *f_purity_fit = new TF1("f_purity_fit", "[0]*TMath::Erf((x - [1])/[2])", xMin, xMax);
    f_purity_fit->SetParameters(1.0, 15.0, 5.0);
    if(fitoption ==1)
    {
        f_purity_fit = new TF1("f_purity_fit", "([0] + [1]*x) / (1 + [2]*[2]*x)", xMin, xMax);
        f_purity_fit->SetParameters(pade_purity_p0, pade_purity_p1, pade_purity_p2);
    }
    

    gpurity->Fit(f_purity_fit, "REMN", "", xMin, xMax);

    TGraphErrors *grFineConf = new TGraphErrors(nFinePoints);
    grFineConf->SetName("grFineConf");

    // Fill the graph with x-values spaced uniformly
    for (int i = 0; i < nFinePoints; ++i) {
        double x = xMin + i * (xMax - xMin) / (nFinePoints - 1);
        grFineConf->SetPoint(i, x, 0); 
    }


    TGraphErrors *confInt = new TGraphErrors(gpurity->GetN());
    confInt->SetName("confInt");
    for (int i = 0; i < gpurity->GetN(); ++i) {
        confInt->SetPoint(i, gpurity->GetX()[i], 0); // y=0 is a placeholder
    }

    //Compute confidence intervals (1 sigma = 68.3%)
    TVirtualFitter *fitter = TVirtualFitter::GetFitter();
    fitter->GetConfidenceIntervals(confInt, 0.683);
    fitter->GetConfidenceIntervals(grFineConf, 0.683);



    for (int i = 0; i < gpurity->GetN(); ++i) {
        double x = gpurity->GetX()[i];
        double y = gpurity->GetY()[i];
        double ey = gpurity->GetErrorY(i);
        double err_low = confInt->GetErrorYlow(i);
        double err_high = confInt->GetErrorYhigh(i);
        //std::cout << "x: " << x << " y: " << y << " ey: " << ey << " err_low: " << err_low << " err_high: " << err_high << std::endl;
    }

    TF1 *f_purity_leak_fit = new TF1("f_purity_leak_fit", "[0]*TMath::Erf((x - [1])/[2])", xMin, xMax);
    f_purity_leak_fit->SetParameters(1.0, 15.0, 5.0);

    if(fitoption ==1)
    {
        f_purity_leak_fit = new TF1("f_purity_leak_fit", "([0] + [1]*x) / (1 + [2]*[2]*x)", xMin, xMax);
        f_purity_leak_fit->SetParameters(pade_purity_leak_p0, pade_purity_leak_p1, pade_purity_leak_p2);
    }

    //TF1 *f_purity_leak_fit = new TF1("f_purity_leak_fit", "pol3", 8, 35);
    //TF1 *f_purity_leak_fit = new TF1("f_purity_leak_fit", "([0] + [1]*x) / (1 + [2]*[2]*x)", xMin, xMax);
    //    f_purity_leak_fit->SetParameters(0.5, 0.5, 0.5);




    gpurity_leak->Fit(f_purity_leak_fit, "REMN","", xMin, xMax);
    gpurity_leak->Fit(f_purity_leak_fit, "REMN","", xMin, xMax);
double chi2 = f_purity_leak_fit->GetChisquare();
int ndf = f_purity_leak_fit->GetNDF();
double pvalue = f_purity_leak_fit->GetProb();

std::cout << "chi2 = " << chi2 << "\n";
std::cout << "ndf = " << ndf << "\n";
std::cout << "chi2/ndf = " << chi2/ndf << "\n";
std::cout << "p-value = " << pvalue << "\n";

    TGraphErrors *grFineConf_leak = new TGraphErrors(nFinePoints);
    grFineConf_leak->SetName("grFineConf_leak");

    // Fill the graph with x-values spaced uniformly
    for (int i = 0; i < nFinePoints; ++i) {
        double x = xMin + i * (xMax - xMin) / (nFinePoints - 1);
        grFineConf_leak->SetPoint(i, x, 0); 
    }

    TGraphErrors *confInt_leak = new TGraphErrors(gpurity_leak->GetN());
    confInt_leak->SetName("confInt_leak");
    for (int i = 0; i < gpurity_leak->GetN(); ++i) {
        confInt_leak->SetPoint(i, gpurity_leak->GetX()[i], 0); // y=0 is a placeholder
    }

    //Compute confidence intervals (1 sigma = 68.3%)
    TVirtualFitter *fitter_leak = TVirtualFitter::GetFitter();
    fitter_leak->GetConfidenceIntervals(confInt_leak, 0.683);
    fitter_leak->GetConfidenceIntervals(grFineConf_leak, 0.683);

    for (int i = 0; i < gpurity_leak->GetN(); ++i) {
        double x = gpurity_leak->GetX()[i];
        double y = gpurity_leak->GetY()[i];
        double ey = gpurity_leak->GetErrorY(i);
        double err_low = confInt_leak->GetErrorYlow(i);
        double err_high = confInt_leak->GetErrorYhigh(i);
        //std::cout << "x: " << x << " y: " << y << " ey: " << ey << " err_low: " << err_low << " err_high: " << err_high << std::endl;
    }

    if (isMC)
    {
        // Truth-purity fit (kept for plotting; not used in the closure correction).
        f_purity_truth_fit = new TF1("f_purity_truth_fit", "[0]*TMath::Erf((x - [1])/[2])", xMin, xMax);
        f_purity_truth_fit->SetParameters(1.0, 15.0, 5.0);
        if (fitoption == 1)
        {
            f_purity_truth_fit = new TF1("f_purity_truth_fit", "([0] + [1]*x) / (1 + [2]*[2]*x)", xMin, xMax);
            f_purity_truth_fit->SetParameters(0.5, 0.5, 0.5);
        }
        g_purity_truth->Fit(f_purity_truth_fit, "REMN", "", xMin, xMax);
        g_purity_truth->Fit(f_purity_truth_fit, "REMN", "", xMin, xMax);

        // Bin-by-bin closure correction = (raw truth purity) / (raw leakage-corrected ABCD purity),
        // computed at the same per-ET-bin sampling as the underlying TGraphs. Then fit with a
        // smooth function (order = analysis.mc_purity_corr_fitorder, default pol1) over
        // [xMin, xMax] so applying the correction on data is robust against per-bin noise.
        g_mc_purity_fit_ratio_out = new TGraphErrors(g_purity_truth->GetN());
        g_mc_purity_fit_ratio_out->SetName("g_mc_purity_fit_ratio");
        for (int i = 0; i < g_purity_truth->GetN(); ++i)
        {
            double x = 0, y_truth = 0;
            g_purity_truth->GetPoint(i, x, y_truth);
            double x_leak = 0, y_leak = 0;
            gpurity_leak->GetPoint(i, x_leak, y_leak);
            double ratio = 1.0;
            double err = 0.0;
            if (y_leak != 0 && !std::isnan(y_leak) && !std::isinf(y_leak)
                            && !std::isnan(y_truth) && !std::isinf(y_truth))
            {
                ratio = y_truth / y_leak;
                double truth_err = 0.5 * (g_purity_truth->GetErrorYlow(i) + g_purity_truth->GetErrorYhigh(i));
                double leak_err  = gpurity_leak->GetErrorY(i);
                double rel2 = 0.0;
                if (y_truth > 0) rel2 += (truth_err / y_truth) * (truth_err / y_truth);
                if (y_leak  > 0) rel2 += (leak_err  / y_leak ) * (leak_err  / y_leak );
                err = std::abs(ratio) * std::sqrt(rel2);
            }
            g_mc_purity_fit_ratio_out->SetPoint(i, x, ratio);
            g_mc_purity_fit_ratio_out->SetPointError(i, 0, err);
        }
        // Smooth fit over the analysis range; saved separately (named f_mc_purity_corr_fit)
        // and used as the actual correction on the data side. Polynomial order is
        // config-driven via analysis.mc_purity_corr_fitorder (default 1 = pol1).
        TString purity_corr_formula = TString::Format("pol%d", mc_purity_corr_fitorder);
        TF1 *f_mc_purity_corr_fit_out = new TF1("f_mc_purity_corr_fit", purity_corr_formula, xMin, xMax);
        f_mc_purity_corr_fit_out->SetParameter(0, 1.0);
        for (int ip = 1; ip <= mc_purity_corr_fitorder; ++ip)
            f_mc_purity_corr_fit_out->SetParameter(ip, 0.0);
        g_mc_purity_fit_ratio_out->Fit(f_mc_purity_corr_fit_out, "REMN", "", xMin, xMax);
        g_mc_purity_fit_ratio_out->Fit(f_mc_purity_corr_fit_out, "REMN", "", xMin, xMax);
        // Stash the fit in a global so the writer block can persist it.
        gROOT->GetListOfFunctions()->Add(f_mc_purity_corr_fit_out);
    }

    if(fit_purity_dis)
    {
        //apply the confint to h_tight_iso_cluster_signal_data
        h_data_sub = (TH1F*) h_tight_iso_cluster_signal_data->Clone("h_data_sub");

        for (int ibin = 1; ibin <= h_tight_iso_cluster_signal_data->GetNbinsX(); ibin++)
        {
            double NA_count = h_tight_iso_cluster_signal_data->GetBinContent(ibin);
            double NA_err = h_tight_iso_cluster_signal_data->GetBinError(ibin);
            //double NA_purity = gpurity->GetY()[ibin - 1];
            double NA_purity = f_purity_fit->Eval(h_tight_iso_cluster_signal_data->GetBinCenter(ibin));
            //get the upper and lower error for confint
            double pusity_fit_low = f_purity_fit->Eval(h_tight_iso_cluster_signal_data->GetBinCenter(ibin)) - confInt->GetErrorYlow(ibin - 1);
            double pusity_fit_high = f_purity_fit->Eval(h_tight_iso_cluster_signal_data->GetBinCenter(ibin)) + confInt->GetErrorYhigh(ibin - 1);
            if(fittingerror == -1)
            {
                NA_purity = pusity_fit_low;
            }
            if(fittingerror == 1)
            {
                NA_purity = pusity_fit_high;
            }
            // Purity is a signal fraction and cannot exceed unity by
            // construction. Some systematic variations push the fitted ABCD
            // purity (or its +1sigma CI edge) slightly above 1, which would
            // imply more signal than data in region A. Clamp to the physical
            // ceiling so the variant cross-section (and hence the purity
            // systematic) is not inflated by an unphysical value. No-op for
            // nominal, whose purity is below 1 in every bin.
            if (NA_purity > 1.0)
            {
                std::cout << "[purity-clamp][noleak] ibin=" << ibin
                          << " ET=" << h_tight_iso_cluster_signal_data->GetBinCenter(ibin)
                          << " raw purity=" << NA_purity << " -> 1.0" << std::endl;
                NA_purity = 1.0;
            }
            double NA_sig_count = NA_count * NA_purity;
            double NA_sig_err = NA_err * NA_purity;

            h_data_sub->SetBinContent(ibin, NA_sig_count);
            h_data_sub->SetBinError(ibin, NA_sig_err);
        }
        
        h_data_sub_leak = (TH1F*) h_tight_iso_cluster_signal_data->Clone("h_data_sub_leak");
        
        for (int ibin = 1; ibin <= h_tight_iso_cluster_signal_data->GetNbinsX(); ibin++)
        {
            double NA_count = h_tight_iso_cluster_signal_data->GetBinContent(ibin);
            double NA_err = h_tight_iso_cluster_signal_data->GetBinError(ibin);
            //double NA_purity = gpurity_leak->GetY()[ibin - 1];
            double x_center = h_tight_iso_cluster_signal_data->GetBinCenter(ibin);
            double mc_corr_ratio = 1.0;
            if (!isMC && mc_purity_correction && f_mc_purity_corr_fit)
            {
                // Evaluate the smooth fit (config-driven order, default pol1) of the bin-by-bin
                // truth/ABCD ratio, NOT the raw TGraph. This avoids interpolation noise.
                mc_corr_ratio = f_mc_purity_corr_fit->Eval(x_center);
                if (mc_corr_ratio <= 0 || std::isnan(mc_corr_ratio) || std::isinf(mc_corr_ratio))
                {
                    mc_corr_ratio = 1.0;
                }
            }
            double NA_purity = f_purity_leak_fit->Eval(x_center) * mc_corr_ratio;
            //get the upper and lower error for confint
            double pusity_fit_low = (f_purity_leak_fit->Eval(x_center) - confInt_leak->GetErrorYlow(ibin - 1)) * mc_corr_ratio;
            double pusity_fit_high = (f_purity_leak_fit->Eval(x_center) + confInt_leak->GetErrorYhigh(ibin - 1)) * mc_corr_ratio;
            if(fittingerror == -1)
            {
                NA_purity = pusity_fit_low;
            }
            if(fittingerror == 1)
            {
                NA_purity = pusity_fit_high;
            }

            // Purity is a signal fraction and cannot exceed unity by
            // construction. Some systematic variations push the leakage-
            // corrected ABCD purity (or its +1sigma CI edge, optionally scaled
            // by the MC closure ratio) slightly above 1, which would imply more
            // signal than data in region A. Clamp to the physical ceiling so the
            // variant cross-section (and hence the purity systematic) is not
            // inflated by an unphysical value. No-op for nominal, whose purity
            // is below 1 in every bin. This is the leakage-corrected path that
            // feeds the nominal cross-section.
            if (NA_purity > 1.0)
            {
                std::cout << "[purity-clamp][leak] ibin=" << ibin
                          << " ET=" << x_center
                          << " raw purity=" << NA_purity << " -> 1.0" << std::endl;
                NA_purity = 1.0;
            }

            double NA_sig_count = NA_count * NA_purity;
            double NA_sig_err = NA_err * NA_purity;

            h_data_sub_leak->SetBinContent(ibin, NA_sig_count);
            h_data_sub_leak->SetBinError(ibin, NA_sig_err);
        }

    }
    


    TH1F *h_data_sub_copy = (TH1F *)h_data_sub->Clone("h_data_sub_copy");
    TH1F *h_data_sub_leak_copy = (TH1F *)h_data_sub_leak->Clone("h_data_sub_leak_copy");

    scale_histogram(h_data_sub_copy, luminosity);
    scale_histogram(h_data_sub_leak_copy, luminosity);

    TH1F *h_tight_iso_cluster_signal_copy = (TH1F *)h_tight_iso_cluster_signal->Clone("h_tight_iso_cluster_signal_copy");

    scale_histogram(h_tight_iso_cluster_signal_copy, (simluminosity * 0.3555));

    std::cout << "n bins: " << h_data_sub_copy->GetNbinsX() << "nbins sim: " << h_tight_iso_cluster_signal_copy->GetNbinsX() << std::endl;

    // h_data_sub_copy->Divide(h_tight_iso_cluster_signal_copy);
    // h_data_sub_leak_copy->Divide(h_tight_iso_cluster_signal_copy);

    /*
    // purity correction with and without leakage correction
    //TGraph for the purity
    TGraphAsymmErrors *gpurity = new TGraphAsymmErrors(h_data_sub, h_tight_iso_cluster_signal_data);
    gpurity->SetName("gpurity");
    gpurity->SetTitle("Purity");
    gpurity->GetXaxis()->SetTitle("E_{T} [GeV]");

    TGraphAsymmErrors *gpurity_leak = new TGraphAsymmErrors(h_data_sub_leak, h_tight_iso_cluster_signal_data);
    gpurity_leak->SetName("gpurity_leak");
    gpurity_leak->SetTitle("Purity with leakage correction");
    gpurity_leak->GetXaxis()->SetTitle("E_{T} [GeV]");
    */

    // place holder for reweighting before unfolding
    TFile *funfold = new TFile(unfoldinput.c_str(), "READ");

    std::string responsematrixname = "response_matrix_full" + histogram_postfix;
    RooUnfoldResponse *response = (RooUnfoldResponse *)funfold->Get(responsematrixname.c_str());

    RooUnfoldResponse *response_reweighted = nullptr;

    RooUnfoldResponse *response_leak_reweighted = nullptr;

    std::string truth_with_cut_name = "h_pT_truth_response" + histogram_postfix;

    TH1D *h_truth_with_cut = (TH1D *)funfold->Get(truth_with_cut_name.c_str());

    int niterations_total = 10;
    int resultit = configYaml["analysis"]["unfold"]["resultit"].as<int>();
    int resultleak = configYaml["analysis"]["unfold"]["resultleak"].as<int>(); // 0 for no leakage correction, 1 for leakage correction
    int reweight = configYaml["analysis"]["unfold"]["reweight"].as<int>();     // 0 for no reweighting, 1 for reweighting
    //we move the reweight handling to updtream turn it off here
    reweight = 0;
    // T2 (unfold_prior_flat): when set, replace the response prior with a
    // uniform truth distribution before Bayes unfolding. With finite iter
    // (resultit=2 nominal), prior choice biases the unfolded result; the
    // difference vs the MC-shape prior is the prior systematic.
    int flat_prior = configYaml["analysis"]["unfold"]["flat_prior"].as<int>(0);

    std::vector<TH1D *> h_unfold_sub_list;
    std::vector<TH1D *> h_unfold_sub_list_copy;

    std::vector<TH1D *> h_unfold_sub_leak_list;
    std::vector<TH1D *> h_unfold_sub_leak_list_copy;

    TH1D *h_reweight_factor;
    TH1D *h_reweight_factor_leak;

    // do reweighting
    if (reweight == 1)
    {
        TH1F *h_reco = (TH1F *)response->Hmeasured();
        TH1F *h_truth = (TH1F *)response->Htruth();
        TH2F *h_response = (TH2F *)response->Hresponse();

        h_reweight_factor = (TH1D *)h_reco->Clone("h_reweight_factor");
        h_reweight_factor->Reset();
        h_reweight_factor_leak = (TH1D *)h_reco->Clone("h_reweight_factor_leak");
        h_reweight_factor_leak->Reset();

        TH2F *h_response_reweighted = (TH2F *)h_response->Clone("h_response_reweighted");
        h_response_reweighted->Reset();

        TH2F *h_response_leak_reweighted = (TH2F *)h_response->Clone("h_response_leak_reweighted");
        h_response_leak_reweighted->Reset();
        // looping over Bins and reweight
        // first check if h_data_sub and h_data_sub_leak have the same binning with h_reco
        if (h_data_sub->GetNbinsX() != h_reco->GetNbinsX())
        {
            std::cout << "h_data_sub and h_reco have different binning" << std::endl;
            return;
        }
        if (h_data_sub_leak->GetNbinsX() != h_reco->GetNbinsX())
        {
            std::cout << "h_data_sub_leak and h_reco have different binning" << std::endl;
            return;
        }

        for (int ibin = 1; ibin <= h_reco->GetNbinsX(); ibin++)
        {
            float reco_val = h_reco->GetBinContent(ibin);
            float truth_val = h_truth->GetBinContent(ibin);
            float data_sub_val = h_data_sub->GetBinContent(ibin);
            float data_sub_leak_val = h_data_sub_leak->GetBinContent(ibin);

            float normcount = h_reco->Integral() / h_data_sub->Integral();

            float reweight_val = data_sub_val / reco_val * normcount;
            float reweight_leak_val = data_sub_leak_val / reco_val * normcount;
            if (reco_val == 0)
            {
                reweight_val = 0;
                reweight_leak_val = 0;
            }
            h_reweight_factor->SetBinContent(ibin, reweight_val);
            h_reweight_factor->SetBinError(ibin, 0);
            h_reweight_factor_leak->SetBinContent(ibin, reweight_leak_val);
            h_reweight_factor_leak->SetBinError(ibin, 0);

            std::cout << "ibin: " << ibin << " reweight_val: " << reweight_val << " reweight_leak_val: " << reweight_leak_val << std::endl;

            for (int jbin = 1; jbin <= h_response->GetNbinsY(); jbin++)
            {
                float response_val = h_response->GetBinContent(ibin, jbin);
                float response_leak_val = h_response->GetBinContent(ibin, jbin);

                h_response_reweighted->SetBinContent(ibin, jbin, response_val * reweight_val);
                h_response_reweighted->SetBinError(ibin, jbin, h_response->GetBinError(ibin, jbin) * reweight_val);
                h_response_leak_reweighted->SetBinContent(ibin, jbin, response_leak_val * reweight_leak_val);
                h_response_leak_reweighted->SetBinError(ibin, jbin, h_response->GetBinError(ibin, jbin) * reweight_leak_val);
            }
        }
        TH1F *h_reweighted_reco = (TH1F *)h_response_reweighted->ProjectionX("h_reweighted_reco");
        TH1F *h_reweighted_truth = (TH1F *)h_response_reweighted->ProjectionY("h_reweighted_truth");
        response_reweighted = new RooUnfoldResponse(h_reweighted_reco, h_reweighted_truth, h_response_reweighted, "response_reweighted", "response_reweighted", false);

        TH1F *h_reweighted_leak_reco = (TH1F *)h_response_leak_reweighted->ProjectionX("h_reweighted_leak_reco");
        TH1F *h_reweighted_leak_truth = (TH1F *)h_response_leak_reweighted->ProjectionY("h_reweighted_leak_truth");

        response_leak_reweighted = new RooUnfoldResponse(h_reweighted_leak_reco, h_reweighted_leak_truth, h_response_leak_reweighted, "response_leak_reweighted", "response_leak_reweighted", false);
    }
    else
    {
        response_reweighted = response;
        response_leak_reweighted = response;
    }

    // T2: flatten response truth prior in place. RooUnfoldBayes consumes
    // Htruth as the iter-0 prior in absolute terms, so we set each bin to
    // (Integral / N_bins) to preserve total normalization while making the
    // SHAPE uniform — that is the actual prior systematic (shape, not norm).
    if (flat_prior)
    {
        std::cout << "[unfold] flat_prior=1 — flattening response truth prior shape" << std::endl;
        auto flatten_truth = [](RooUnfoldResponse *resp) {
            if (!resp) return;
            TH1 *htruth = (TH1 *)resp->Htruth();
            int nbins = htruth->GetNbinsX();
            double total = htruth->Integral();
            double per_bin = (nbins > 0) ? (total / nbins) : 1.0;
            std::cout << "  [flat_prior] nbins=" << nbins
                      << " total=" << total
                      << " per_bin=" << per_bin << std::endl;
            for (int ib = 1; ib <= nbins; ++ib) {
                htruth->SetBinContent(ib, per_bin);
                htruth->SetBinError(ib, 0.0);
            }
        };
        flatten_truth(response_reweighted);
        if (response_leak_reweighted != response_reweighted)
            flatten_truth(response_leak_reweighted);
    }

    for (int i = 0; i < niterations_total; i++)
    {
        RooUnfoldBayes *unfold = new RooUnfoldBayes(response_reweighted, h_data_sub, i + 1);
        //  RooUnfoldBayes *unfold = new RooUnfoldBayes(response, h_data_sub, i + 1);
        //  RooUnfoldBayes *unfold = new RooUnfoldBayes(response, h_tight_iso_cluster_signal_data, i + 1);
        // RooUnfoldBayes *unfold = new RooUnfoldBayes(response_reweighted, h_truth_iso_cluster_data, i + 1);

        TH1D *h_unfold_sub = (TH1D *)unfold->Hunfold()->Clone(Form("h_unfold_sub_%d", i + 1));
        h_unfold_sub->SetName(Form("h_unfold_sub_%d", i + 1));

        TH1D *h_reweight_response = (TH1D *)unfold->Hmeasured();
        TH1D *h_reweight_truth = (TH1D *)response_reweighted->Htruth();
        for (int ibin = 1; ibin <= h_reweight_response->GetNbinsX(); ibin++)
        {
            std::cout << "iteration: " << i + 1 << " h_truth_iso_cluster_data->GetBinContent(1): " << h_truth_iso_cluster_data->GetBinContent(ibin) << std::endl;

            std::cout << "iteration: " << i + 1 << " h_reweight_response->GetBinContent(1): " << h_reweight_response->GetBinContent(ibin) << std::endl;

            std::cout << "iteration: " << i + 1 << " h_reweight_truth->GetBinContent(1): " << h_reweight_truth->GetBinContent(ibin) << std::endl;
            // h_unfold_sub

            std::cout << "iteration: " << i + 1 << " h_unfold_sub->GetBinContent(1): " << h_unfold_sub->GetBinContent(ibin) << std::endl;
        }

        // h_unfold_sub = h_reweight_truth;

        // scale by the lumi and cross section
        // loop over the bins
        /*
        for (int ibin = 1; ibin <= h_unfold_sub->GetNbinsX(); ibin++)
        {
            float binwidth = h_unfold_sub->GetBinWidth(ibin);
            float bincenter = h_unfold_sub->GetBinCenter(ibin);

            float scale = 1.0 / binwidth / luminosity;

            h_unfold_sub->SetBinContent(ibin, h_unfold_sub->GetBinContent(ibin) * scale);
            h_unfold_sub->SetBinError(ibin, h_unfold_sub->GetBinError(ibin) * scale);
        }
        */
        scale_histogram(h_unfold_sub, luminosity);
        h_unfold_sub_list.push_back(h_unfold_sub);
        h_unfold_sub_list_copy.push_back((TH1D *)h_unfold_sub->Clone(Form("h_unfold_sub_%d_copy", i + 1)));

        RooUnfoldBayes *unfold_leak = new RooUnfoldBayes(response_leak_reweighted, h_data_sub_leak, i + 1);
        //   RooUnfoldBayes *unfold_leak = new RooUnfoldBayes(response, h_data_sub_leak, i + 1);
        //   RooUnfoldBayes *unfold_leak = new RooUnfoldBayes(response, h_tight_iso_cluster_signal_data, i + 1);
        // RooUnfoldBayes *unfold_leak = new RooUnfoldBayes(response_leak_reweighted, h_truth_iso_cluster_data, i + 1);
        TH1D *h_unfold_sub_leak = (TH1D *)unfold_leak->Hunfold()->Clone(Form("h_unfold_sub_leak_%d", i + 1));
        h_unfold_sub_leak->SetName(Form("h_unfold_sub_leak_%d", i + 1));
        // h_unfold_sub_leak = h_reweight_truth;
        //  scale by the lumi and cross section
        //  loop over the bins
        /*
        for (int ibin = 1; ibin <= h_unfold_sub_leak->GetNbinsX(); ibin++)
        {
            float binwidth = h_unfold_sub_leak->GetBinWidth(ibin);
            float bincenter = h_unfold_sub_leak->GetBinCenter(ibin);

            float scale = 1.0 / binwidth / luminosity;

            h_unfold_sub_leak->SetBinContent(ibin, h_unfold_sub_leak->GetBinContent(ibin) * scale);
            h_unfold_sub_leak->SetBinError(ibin, h_unfold_sub_leak->GetBinError(ibin) * scale);
        }
        */
        scale_histogram(h_unfold_sub_leak, luminosity);
        h_unfold_sub_leak_list.push_back(h_unfold_sub_leak);
        h_unfold_sub_leak_list_copy.push_back((TH1D *)h_unfold_sub_leak->Clone(Form("h_unfold_sub_leak_%d_copy", i + 1)));
    }

    // efficiency correction

    for (int i = 0; i < niterations_total; i++)
    {
        TH1D *h_unfold_sub = h_unfold_sub_list[i];
        TH1D *h_unfold_sub_leak = h_unfold_sub_leak_list[i];

        for (int ibin = 1; ibin <= h_unfold_sub->GetNbinsX(); ibin++)
        {
            float eff_reco_val = eff_reco->GetEfficiency(ibin);
            float eff_iso_val = eff_iso->GetEfficiency(ibin);
            float eff_id_val = eff_id->GetEfficiency(ibin);

            float vertex_eff_val = h_truth_pT_vertexcut_mbd_cut->GetBinContent(ibin) / h_truth_pT_vertexcut->GetBinContent(ibin) + mbd_eff_scale;

            float total_eff = eff_reco_val * eff_iso_val * eff_id_val * vertex_eff_val;

            // std::cout << "ibin: " << ibin << " eff_reco_val: " << eff_reco_val << " eff_iso_val: " << eff_iso_val << " eff_id_val: " << eff_id_val << std::endl;

            h_unfold_sub->SetBinContent(ibin, h_unfold_sub->GetBinContent(ibin) / total_eff);
            h_unfold_sub->SetBinError(ibin, h_unfold_sub->GetBinError(ibin) / total_eff);

            h_unfold_sub_leak->SetBinContent(ibin, h_unfold_sub_leak->GetBinContent(ibin) / total_eff);
            h_unfold_sub_leak->SetBinError(ibin, h_unfold_sub_leak->GetBinError(ibin) / total_eff);
        }
    }

    // loop over bins of h_pythia_truth and scale by the lumi and cross section
    /*
    for (int ibin = 1; ibin <= h_pythia_truth->GetNbinsX(); ibin++)
    {
        float binwidth = h_pythia_truth->GetBinWidth(ibin);
        float bincenter = h_pythia_truth->GetBinCenter(ibin);

        float scale = 1.0 / binwidth / simluminosity;

        h_pythia_truth->SetBinContent(ibin, h_pythia_truth->GetBinContent(ibin) * scale);
        h_pythia_truth->SetBinError(ibin, h_pythia_truth->GetBinError(ibin) * scale);
    }
    */
    scale_histogram(h_pythia_truth, simluminosity);
    scale_histogram(h_pythia_truth_novtx, simluminosity);
    // h_truth_with_cut
    /*
    for (int ibin = 1; ibin <= h_truth_with_cut->GetNbinsX(); ibin++)
    {
        float binwidth = h_truth_with_cut->GetBinWidth(ibin);
        float bincenter = h_truth_with_cut->GetBinCenter(ibin);

        float scale = 1.0 / binwidth / simluminosity;

        h_truth_with_cut->SetBinContent(ibin, h_truth_with_cut->GetBinContent(ibin) * scale);
        h_truth_with_cut->SetBinError(ibin, h_truth_with_cut->GetBinError(ibin) * scale);
    }
    */
    scale_histogram(h_truth_with_cut, simluminosity);

    // h_data_sub and h_data_sub_leak
    /*
    for (int ibin = 1; ibin <= h_data_sub->GetNbinsX(); ibin++)
    {
        float binwidth = h_data_sub->GetBinWidth(ibin);
        float bincenter = h_data_sub->GetBinCenter(ibin);

        float scale = 1.0 / binwidth / luminosity;

        h_data_sub->SetBinContent(ibin, h_data_sub->GetBinContent(ibin) * scale);
        h_data_sub->SetBinError(ibin, h_data_sub->GetBinError(ibin) * scale);

        h_data_sub_leak->SetBinContent(ibin, h_data_sub_leak->GetBinContent(ibin) * scale);
        h_data_sub_leak->SetBinError(ibin, h_data_sub_leak->GetBinError(ibin) * scale);
    }
    */
    scale_histogram(h_data_sub, luminosity);
    scale_histogram(h_data_sub_leak, luminosity);

    scale_histogram(h_common_cluster_data, luminosity);

    scale_histogram(h_tight_iso_cluster_signal_data, luminosity);

    TFile *fout = new TFile(outfilename.c_str(), "RECREATE");

    h_pythia_truth->Write();
    h_pythia_truth_novtx->Write();
    h_truth_with_cut->Write();

    for (int i = 0; i < niterations_total; i++)
    {
        std::cout << "i: " << i << std::endl;
        h_unfold_sub_list[i]->Write();
        h_unfold_sub_list_copy[i]->Write();
        h_unfold_sub_leak_list[i]->Write();
        h_unfold_sub_leak_list_copy[i]->Write();

        if (i + 1 == resultit)
        {
            if (resultleak == 0)
            {
                TH1F *h_unfold_sub_result = (TH1F *)h_unfold_sub_list[i]->Clone("h_unfold_sub_result");
                h_unfold_sub_result->Write("h_unfold_sub_result");
                TH1F *h_unfold_sub_result_woeff = (TH1F *)h_unfold_sub_list_copy[i]->Clone("h_unfold_sub_result_woeff");
                h_unfold_sub_result_woeff->Write("h_unfold_sub_result_woeff");
            }
            else
            {
                TH1F *h_unfold_sub_result = (TH1F *)h_unfold_sub_leak_list[i]->Clone("h_unfold_sub_result");
                h_unfold_sub_result->Write("h_unfold_sub_result");
                TH1F *h_unfold_sub_result_woeff = (TH1F *)h_unfold_sub_leak_list_copy[i]->Clone("h_unfold_sub_result_woeff");
                h_unfold_sub_result_woeff->Write("h_unfold_sub_result_woeff");
            }
        }
    }
    gpurity->Write();
    gpurity_leak->Write();
    g_purity_truth->Write();
    if (f_purity_truth_fit)
    {
        f_purity_truth_fit->Write();
    }
    if (g_mc_purity_fit_ratio_out)
    {
        g_mc_purity_fit_ratio_out->Write();
        TF1 *f_mc_purity_corr_fit_persist = (TF1 *)gROOT->GetListOfFunctions()->FindObject("f_mc_purity_corr_fit");
        if (f_mc_purity_corr_fit_persist)
        {
            f_mc_purity_corr_fit_persist->Write();
            gROOT->GetListOfFunctions()->Remove(f_mc_purity_corr_fit_persist);
        }
    }
    g_mbd_eff->Write();
    g_mbd_eff_north->Write();
    g_mbd_eff_south->Write();
    std::cout << "saving histograms" << std::endl;

    h_data_sub_copy->Write();
    h_data_sub_leak_copy->Write();
    h_tight_iso_cluster_signal_copy->Write();
    h_R->Write();
    h_R_notmatch->Write();
    h_leak_B->Write();
    h_leak_C->Write();
    h_leak_D->Write();
    if (reweight == 1)
    {
        h_reweight_factor->Write();
        h_reweight_factor_leak->Write();
    }
    h_data_sub->Write();
    h_data_sub_leak->Write();
    h_tight_iso_cluster_signal_data->Write();
    h_common_cluster_data->Write();
    for(auto h_purity: h_purity_list){
        h_purity->Write();
    }
    for(auto h_purity_leak: h_purity_leak_list){
        h_purity_leak->Write();
    }
    for(auto f_purity: f_purity_list){
        f_purity->Write();
    }
    for(auto f_purity_leak: f_purity_leak_list){
        f_purity_leak->Write();
    }
    f_purity_fit->Write();
    f_purity_leak_fit->Write();
    confInt->Write();
    confInt_leak->Write();
    grFineConf->Write();
    grFineConf_leak->Write();
    if (reweight == 1)
    {
        response_reweighted->Hresponse()->Write();
        response_leak_reweighted->Hresponse()->Write();
    }
    else
    {
        std::cout << "response->Hresponse()->Write()" << std::endl;
        response->Hresponse()->Write();
    }

    // ---------------------------------------------------------------------
    // Tower-index (ietacent x iphicent) acceptance maps:
    //   - Read the 4 TH2s from the jet merge. The jet MC is the inclusive sample
    //     (it already contains the prompt-photon processes at the PYTHIA rate), so
    //     it is the inclusive MC tower map on its own. The photon MC is not added.
    //   - Read the 4 TH2s from data (fdatain) as-is.
    //   - Write both with _mc_inclusive / _data suffix into Photon_final_{var_type}.root
    //     for direct downstream plotting (no re-hadd needed).
    // ---------------------------------------------------------------------
    fout->cd();
    const std::vector<std::string> tower_levels = {
        "preselect", "common", "tight", "tight_iso"
    };
    std::string jetinput = configYaml["output"]["eff_outfile"].as<std::string>()
                         + "_jet_" + var_type + ".root";
    TFile *fjetin = TFile::Open(jetinput.c_str(), "READ");
    for (const auto &lvl : tower_levels)
    {
        std::string hname = "h_etaphi_tower_" + lvl;
        TH2F *h_jet = (fjetin && !fjetin->IsZombie()) ? (TH2F *) fjetin->Get(hname.c_str()) : nullptr;
        if (!h_jet) {
            std::cout << "WARNING: missing " << hname << " in " << jetinput << std::endl;
            continue;
        }
        TH2F *h_inc = (TH2F *) h_jet->Clone((hname + "_mc_inclusive").c_str());
        h_inc->SetDirectory(fout);
        h_inc->Write((hname + "_mc_inclusive").c_str(), TObject::kOverwrite);

        // Pass data through too (only if fdatain holds real data, not the
        // isMC-closure-mode reuse of the MC file).
        if (!isMC)
        {
            TH2F *h_dat = (TH2F *) fdatain->Get(hname.c_str());
            if (h_dat) {
                TH2F *h_dat_out = (TH2F *) h_dat->Clone((hname + "_data").c_str());
                h_dat_out->SetDirectory(fout);
                h_dat_out->Write((hname + "_data").c_str(), TObject::kOverwrite);
            }
        }
    }
    if (fjetin) fjetin->Close();

    fout->Write();
    fout->Close();
}
