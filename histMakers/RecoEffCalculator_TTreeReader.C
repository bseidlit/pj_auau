#include <iostream>
#include <string>
#include <fstream>
#include <iterator>
#include <vector>
#include <sstream>
#include <cmath>
#include <limits>
#include <algorithm>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TF1.h>
#include <TSystem.h>
#include <TEfficiency.h>
#include <TObjString.h>
#include <TRandom3.h>
#include <yaml-cpp/yaml.h>
// unfolding
#include <RooUnfoldResponse.h>
#include <RooUnfoldBayes.h>

// R__LOAD_LIBRARY(/sphenix/user/egm2153/calib_study/JetValidation/analysis/roounfold/libRooUnfold.so)
void SaveYamlToRoot(TFile *f, const char *yaml_filename)
{
    // Read YAML file into a string
    std::ifstream yaml_file(yaml_filename);
    std::string yaml_content((std::istreambuf_iterator<char>(yaml_file)),
                             std::istreambuf_iterator<char>());

    // Create a ROOT file and save the YAML string
    TObjString yaml_obj(yaml_content.c_str());
    f->cd();
    yaml_obj.Write("config");
}




void RecoEffCalculator_TTreeReader(const std::string &configname = "config_bdt_nom.yaml",
                                   const std::string &filetype_in = "jet10_aa")
{
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so"); 
    YAML::Node configYaml = YAML::LoadFile(configname);


    /////////////////////////////////////
    // Parse filetype  and input / output filenames
    /////////////////////////////////////
    bool issim = true;
    std::string filetype = filetype_in;
    std::string filetype_base = filetype;
    bool is_pp = true;
    if (filetype.size() > 3 && filetype.compare(filetype.size() - 3, 3, "_aa") == 0) {
        filetype_base = filetype.substr(0, filetype.size() - 3);
        is_pp = false;
    } else if (filetype.size() > 3 && filetype.compare(filetype.size() - 3, 3, "_pp") == 0) {
        filetype_base = filetype.substr(0, filetype.size() - 3);
        is_pp = true;
    }

    if (filetype == "data" || filetype == "data_pp" || filetype == "data_aa")
    {
        issim = false;
    }

    // bool isbackground = (bool)configYaml["input"]["isbackground"].as<int>();
    static const bool isbackground = filetype.find("jet")!=string::npos ? true : false;

    std::string infilename_root_dir   = configYaml["input"]["photon_jet_file_root_dir"].as<std::string>();
    std::string infilename_branch_dir = configYaml["input"]["photon_jet_file_branch_dir"].as<std::string>();
    std::string infilename            = infilename_root_dir + filetype + infilename_branch_dir;

    // For data, allow per-system files: data_file_pp, data_file_aa; else fall back to data_file
    if (!issim)
    {
        if (filetype == "data_aa" && configYaml["input"]["data_file_aa"])
            infilename = configYaml["input"]["data_file_aa"].as<std::string>();
        else if ((filetype == "data_pp" || filetype == "data") && configYaml["input"]["data_file_pp"])
            infilename = configYaml["input"]["data_file_pp"].as<std::string>();
        else
            infilename = configYaml["input"]["data_file"].as<std::string>();
    }

    std::cout << "infilename: " << infilename << std::endl;


    /////////////////////////////////////
    // Define cross section values
    /////////////////////////////////////
    float max_photon_lower = 0;
    float max_photon_upper = 100;
    // unit in pb
    //9.26915e+10 146359.3
    //const float photon5cross = 2.017e+08 * 0.000442571;
    const float photon5cross = 146359.3;
    //1.2613e+08 6944.675
    //const float photon10cross = 3.688e+07 * 0.000181474;
    const float photon10cross = 6944.675;
    //5.2244e+06 130.4461
    //const float photon20cross = 1.571e+05 * 0.000673448;
    const float photon20cross = 130.4461;

    // Hanpu uses unit in b
    const float jet10cross = 3.997e+06;
    const float jet15cross = 4.073e+05;
    const float jet5cross = 1.3878e+08;
    const float jet8cross = 1.3013e+07*0.7;
    const float jet12cross = 1.4903e+06;
    const float jet20cross = 6.2623e+04;
    const float jet30cross = 2.5298e+03;
    const float jet40cross = 1.3553e+02;
    const float jet50cross = 7.3113;

    float max_jet_lower = 0;
    float max_jet_upper = 100;

    float energy_scale_lower = 0;
    float energy_scale_upper = 100;

    float cluster_ET_upper = 100;

    float weight = 1.0;
    float vertex_weight = 1.0;
    float cent_weight = 1.0;
    float cross_weight = 1.0;

    if (filetype_base == "photon10")
    {
        max_photon_lower = 10;
        max_photon_upper = 20;
        weight = photon10cross / photon20cross /0.045;
    }
    else if (filetype_base == "photon12")
    {
        max_photon_lower = 12;
        max_photon_upper = 30;
        weight = 1.0;
    }
    else if (filetype_base == "photon20")
    {
        max_photon_lower = 20;
        max_photon_upper = 100;
        weight = 1.0;
    }
    else if (filetype_base == "jet10")
    {
        max_jet_lower = 10;
        max_jet_upper = 20;
        cluster_ET_upper = 23;
        weight = jet10cross / jet20cross;
        //isbackground = true;
    }
    else if (filetype_base == "jet20")
    {
        max_jet_lower = 20;
        max_jet_upper = 100;
        cluster_ET_upper = 100;
        weight = 1.0;
        //isbackground = true;
    }

    cross_weight = weight;



    /////////////////////////////////////
    // Vertex reweighting
    /////////////////////////////////////
    TH1* h_vertex_reweight = nullptr;
    int vertex_reweight_on = 1;
    std::string vertex_reweight_file = "results/vertex_reweight_bdt_none.root";
    std::string vtx_scan_data_file = "";
    if (issim)
    {
        // optional config knobs (safe defaults)
        vertex_reweight_on = configYaml["analysis"]["vertex_reweight_on"].as<int>(1);
        vertex_reweight_file =
            configYaml["analysis"]["vertex_reweight_file"].as<std::string>("results/vertex_reweight.root");
        vtx_scan_data_file = configYaml["analysis"]["vertex_scan_data_file"].as<std::string>("");
        //vertex_reweight_on = false;
        if (vertex_reweight_on && vtx_scan_data_file.empty())
        {
            TFile* fvtx = TFile::Open(vertex_reweight_file.c_str(), "READ");
            if (!fvtx || fvtx->IsZombie())
            {
                std::cerr << "[VertexReweight] ERROR: cannot open vertex reweight file: " << vertex_reweight_file << std::endl;
                return;
            }

            TH1* htmp = dynamic_cast<TH1*>(fvtx->Get("data_over_MC_ratios/h_zvtx_ratio_data_over_photonJet"));
            if (!htmp)
            {
                std::cerr << "[VertexReweight] ERROR: cannot find histogram " << std::endl;
                return;
            }

            std::string vtx_histname = htmp->GetName();  // save before Close() invalidates htmp
	        h_vertex_reweight = dynamic_cast<TH1*>(htmp->Clone("h_zvtx_ratio_data_over_photonJet_clone"));

            h_vertex_reweight->SetDirectory(nullptr);
            h_vertex_reweight->Sumw2();
            std::cout << "[VertexReweight] Using histogram weights from " << vertex_reweight_file << " : " << vtx_histname << std::endl;
        }
    }


    /////////////////////////////////////
    // Centrality reweighting
    /////////////////////////////////////
    TH1 *h_cent_reweight = nullptr;
    int cent_reweight_on = 0;
    std::string cent_reweight_file = "../reweightingDer/output/centrality_reweighting.root";
    if (issim)
    {
        // preferred keys
        cent_reweight_on = configYaml["analysis"]["centrality_reweight_on"].as<int>(0);
        cent_reweight_file =
            configYaml["analysis"]["centrality_reweight_file"].as<std::string>(cent_reweight_file);
        // backward-compatible key
        if (configYaml["analysis"]["cent_reweight_on"])
            cent_reweight_on = configYaml["analysis"]["cent_reweight_on"].as<int>(cent_reweight_on);
        if (cent_reweight_on)
        {
            TFile *fcent = TFile::Open(cent_reweight_file.c_str(), "READ");
            if (!fcent || fcent->IsZombie())
            {
                std::cerr << "[CentralityReweight] ERROR: cannot open centrality reweight file: " << cent_reweight_file << std::endl;
                return;
            }
            TH1 *htmp = dynamic_cast<TH1 *>(fcent->Get("nom_cent_rw_hist"));
            if (!htmp)
            {
                std::cerr << "[CentralityReweight] ERROR: cannot find histogram 'nom_cent_rw_hist' in " << cent_reweight_file << std::endl;
                fcent->Close();
                delete fcent;
                return;
            }
            h_cent_reweight = dynamic_cast<TH1 *>(htmp->Clone("nom_cent_rw_hist_clone"));
            if (!h_cent_reweight)
            {
                std::cerr << "[CentralityReweight] ERROR: failed to clone histogram 'nom_cent_rw_hist' from "
                          << cent_reweight_file << std::endl;
                fcent->Close();
                delete fcent;
                return;
            }
            h_cent_reweight->SetDirectory(nullptr);
            h_cent_reweight->Sumw2();
            fcent->Close();
            delete fcent;
            std::cout << "[CentralityReweight] Using histogram weights from "
                      << cent_reweight_file << " : nom_cent_rw_hist" << std::endl;
        }
    }



    // std::string infilename = "/sphenix/tg/tg01/commissioning/CaloCalibWG/sli/ppg12/ana450/condorout/combine.root";
    // Build input chain and connect reader
    std::string treename = configYaml["input"]["tree"].as<std::string>();
    TChain chain(treename.c_str());
    chain.Add(infilename.c_str());
    std::string var_type = configYaml["output"]["var_type"].as<std::string>();

    std::string outfilename = configYaml["output"]["eff_outfile"].as<std::string>() + "_" + filetype + "_" + var_type + ".root";

    std::string responsefilename = configYaml["output"]["response_outfile"].as<std::string>() + "_" + filetype + "_" + var_type + ".root";

    if (!issim)
    {
        outfilename = configYaml["output"]["data_outfile"].as<std::string>() + "_" + var_type + ".root";
        // unfolding is only for sim
        responsefilename = "bla.root";
    }

    std::cout << "outfilename: " << outfilename << std::endl;
    std::cout << "responsefilename: " << responsefilename << std::endl;

    // TChain is used instead of a single TTree

    std::string clusternodename = configYaml["input"]["cluster_node_name"].as<std::string>();


    /////////////////////////////////////
    // BDT model
    /////////////////////////////////////
    std::string bdt_model_name = configYaml["input"]["bdt_model_name"].as<std::string>("base");
    //std::string bdt_model_name = "base";

    std::vector<float> bdt_et_bin_edges;
    std::vector<std::string> bdt_et_bin_models;
    bool use_et_binned_bdt = false;

    if (configYaml["input"]["bdt_et_bin_edges"] && configYaml["input"]["bdt_et_bin_models"]) {
        for (auto v : configYaml["input"]["bdt_et_bin_edges"])
            bdt_et_bin_edges.push_back(v.as<float>());
        for (auto v : configYaml["input"]["bdt_et_bin_models"])
            bdt_et_bin_models.push_back(v.as<std::string>());
        use_et_binned_bdt = (bdt_et_bin_models.size() == bdt_et_bin_edges.size() - 1);
        if (!use_et_binned_bdt)
            std::cout << "WARNING: bdt_et_bin_edges/bdt_et_bin_models size mismatch; falling back to single model" << std::endl;
    }


    /////////////////////////////////////
    // Iso threshold
    /////////////////////////////////////
    int iso_threshold = configYaml["analysis"]["iso_threshold"].as<int>(0);
    int iso_hcalonly = configYaml["analysis"]["iso_hcalonly"].as<int>(0);
    int use_topo_iso = configYaml["analysis"]["use_topo_iso"].as<int>(0);
    // we have 0, 0.05, 0.1, 0.2 options
    float iso_emcalinnerr = configYaml["analysis"]["iso_emcalinnerr"].as<float>(0.0);
    std::cout<<"iso_emcalinnerr: "<<iso_emcalinnerr<<std::endl;



 
    float truthisocut = configYaml["analysis"]["truth_iso_max"].as<float>();

    float recoiso_min = configYaml["analysis"]["reco_iso_min"].as<float>();
    float recoiso_max_b = configYaml["analysis"]["reco_iso_max_b"].as<float>();
    float recoiso_max_s = configYaml["analysis"]["reco_iso_max_s"].as<float>();
    float reco_iso_max_cent_c0 = 0.0f, reco_iso_max_cent_c1 = 0.0f, reco_iso_max_cent_c2 = 0.0f;
    if (configYaml["analysis"]["reco_iso_max_cent_coeff"] && configYaml["analysis"]["reco_iso_max_cent_coeff"].size() >= 3)
    {
        reco_iso_max_cent_c0 = configYaml["analysis"]["reco_iso_max_cent_coeff"][0].as<float>();
        reco_iso_max_cent_c1 = configYaml["analysis"]["reco_iso_max_cent_coeff"][1].as<float>();
        reco_iso_max_cent_c2 = configYaml["analysis"]["reco_iso_max_cent_coeff"][2].as<float>();
    }

    float recononiso_min_shift = configYaml["analysis"]["reco_noniso_min_shift"].as<float>();
    float recononiso_max = configYaml["analysis"]["reco_noniso_max"].as<float>();

    float vertexcut = configYaml["analysis"]["vertex_cut"].as<float>();
    std::vector<float> eta_bins = configYaml["analysis"]["eta_bins"].as<std::vector<float>>();
    float eta_acceptance_min = -0.7f;
    float eta_acceptance_max = 0.7f;
    if (!eta_bins.empty())
    {
        eta_acceptance_min = eta_bins.front();
        eta_acceptance_max = eta_bins.back();
    }
    int n_eta_bins = eta_bins.size() - 1;
    double eta_bin_edges[n_eta_bins + 1];
    std::copy(eta_bins.begin(), eta_bins.end(), eta_bin_edges);


    /////////////////////////////////////
    // Bining
    /////////////////////////////////////
    std::vector<float> centrality_bins;
    if (configYaml["analysis"]["centrality_bins"])
      centrality_bins = configYaml["analysis"]["centrality_bins"].as<std::vector<float>>();
    if (centrality_bins.empty())
      centrality_bins = {0.0f, 100.0f};  // default single bin for pp/no-centrality
    int n_cent_bins = centrality_bins.size() - 1;

    std::vector<float> pT_bins = configYaml["analysis"]["pT_bins"].as<std::vector<float>>();
    int n_pT_bins = pT_bins.size() - 1;
    double pT_bin_edges[n_pT_bins + 1];
    double pTmin = pT_bins[0];
    double pTmax = pT_bins[n_pT_bins];

    float b2bjet_dphi = 3 * M_PI / 4;
    if (configYaml["analysis"]["b2bjet_dphi"])
        b2bjet_dphi = configYaml["analysis"]["b2bjet_dphi"].as<float>();

    float  jet_eta = configYaml["analysis"]["jet_eta"].as<float>();

    float b2bjet_pT_min = configYaml["analysis"]["b2bjet_pT_min"].as<float>();
    // Anti-kT tower jets: 3 -> R=0.3, 4 -> R=0.4 (same as EvtCharacter.C)
    const int jet_cone_size = configYaml["analysis"]["jet_cone_size"].as<int>(3);
    const bool use_jet_r04 = (jet_cone_size == 4);
    std::cout << "[RecoEffCalculator] jet_cone_size=" << jet_cone_size
              << " -> Anti-kT Tower r0" << (use_jet_r04 ? "4" : "3") << "_Sub1"
              << ", truth jets AntiKt_Truth_r0" << (use_jet_r04 ? "4" : "3") << std::endl;

    std::vector<float> pT_bins_truth = configYaml["analysis"]["pT_bins_truth"].as<std::vector<float>>();
    int n_pT_bins_truth = pT_bins_truth.size() - 1;
    double pT_bin_edges_truth[n_pT_bins_truth + 1];
    double pTmin_truth = pT_bins_truth[0];
    double pTmax_truth = pT_bins_truth[n_pT_bins_truth];

    std::cout << "n_pT_bins_truth: " << n_pT_bins_truth << std::endl;
    for (int i = 0; i < n_pT_bins_truth + 1; i++)
    {
        pT_bin_edges_truth[i] = pT_bins_truth[i];
        std::cout << "pT_bin_edges_truth: " << pT_bin_edges_truth[i] << std::endl;
    }

    std::copy(pT_bins.begin(), pT_bins.end(), pT_bin_edges);

    // xjgamma binning (variable bins from YAML; default provided)
    std::vector<float> xjgamma_bins;
    if (configYaml["analysis"]["xjgamma_bins"])
        xjgamma_bins = configYaml["analysis"]["xjgamma_bins"].as<std::vector<float>>();
    if (xjgamma_bins.empty())
        xjgamma_bins = {0.0f, 0.2f, 0.4f, 0.6f, 0.8f, 1.0f, 1.2f, 1.4f, 1.6f, 1.8f};
    int n_xj_bins = xjgamma_bins.size() - 1;
    double xj_bin_edges[n_xj_bins + 1];
    for (int i = 0; i < n_xj_bins + 1; ++i)
    {
        xj_bin_edges[i] = xjgamma_bins[i];
    }

    // truth xjgamma binning (variable bins from YAML; default to reco xjgamma bins)
    std::vector<float> xjgamma_bins_truth = configYaml["analysis"]["xjgamma_bins_truth"].as<std::vector<float>>();
    if (xjgamma_bins_truth.empty())
        xjgamma_bins_truth = xjgamma_bins;
    int n_xj_bins_truth = xjgamma_bins_truth.size() - 1;
    double xj_bin_edges_truth[n_xj_bins_truth + 1];
    for (int i = 0; i < n_xj_bins_truth + 1; ++i)
    {
        xj_bin_edges_truth[i] = xjgamma_bins_truth[i];
    }
    const int n_global_xjgamma_truth_bins = n_pT_bins_truth * n_xj_bins_truth;
    const int n_global_xjgamma_reco_bins = n_pT_bins * n_xj_bins;

    int conesize = configYaml["analysis"]["cone_size"].as<int>();

    float reco_min_ET = configYaml["analysis"]["reco_min_ET"].as<float>();

    float eff_dR = configYaml["analysis"]["eff_dR"].as<float>();

    // trigger_used can be either a scalar int or a YAML sequence of ints.
    // Example: trigger_used: [26, 29, 30, 31, 36, 37, 38]
    std::vector<int> trigger_used;
    {
        YAML::Node trigNode = configYaml["analysis"]["trigger_used"];
        if (trigNode && trigNode.IsSequence())
        {
            trigger_used = trigNode.as<std::vector<int>>();
        }
        else
        {
            // backward-compatible: allow single int
            trigger_used.push_back(configYaml["analysis"]["trigger_used"].as<int>());
        }
    }

    float mc_iso_shift = configYaml["analysis"]["mc_iso_shift"].as<float>(0.0);
    float mc_iso_scale = configYaml["analysis"]["mc_iso_scale"].as<float>(1.0);
    /////////////////////////////////////
    // Photon ID cuts
    /////////////////////////////////////
    int n_nt_fail = configYaml["analysis"]["n_nt_fail"].as<int>(1);

    int weta_fail = configYaml["analysis"]["weta_fail"].as<int>(0);
    int wphi_fail = configYaml["analysis"]["wphi_fail"].as<int>(0);
    int e11_to_e33_fail = configYaml["analysis"]["e11_to_e33_fail"].as<int>(0);
    int e32_to_e35_fail = configYaml["analysis"]["e32_to_e35_fail"].as<int>(0);
    int et1_fail = configYaml["analysis"]["et1_fail"].as<int>(0);
    int bdt_fail = configYaml["analysis"]["bdt_fail"].as<int>(0);

    int weta_on = configYaml["analysis"]["weta_on"].as<int>(1);
    int wphi_on = configYaml["analysis"]["wphi_on"].as<int>(1);
    int e11_to_e33_on = configYaml["analysis"]["e11_to_e33_on"].as<int>(1);
    int e32_to_e35_on = configYaml["analysis"]["e32_to_e35_on"].as<int>(1);
    int et1_on = configYaml["analysis"]["et1_on"].as<int>(1);
    int et2_on = configYaml["analysis"]["et2_on"].as<int>(1);
    int et3_on = configYaml["analysis"]["et3_on"].as<int>(1);
    int et4_on = configYaml["analysis"]["et4_on"].as<int>(1);
    int bdt_on = configYaml["analysis"]["bdt_on"].as<int>(1);
    int nosat = configYaml["analysis"]["nosat"].as<int>(0);
    int common_b2bjet_cut = configYaml["analysis"]["common_b2bjet_cut"].as<int>(0);
    float common_b2bjet_pt_min = configYaml["analysis"]["common_b2bjet_pt_min"].as<float>(7.0);

    std::cout << "tight cuts" << std::endl;
    float tight_reta77_min = configYaml["analysis"]["tight"]["reta77_min"].as<float>();
    float tight_reta77_max = configYaml["analysis"]["tight"]["reta77_max"].as<float>();

    float tight_rhad33_max = configYaml["analysis"]["tight"]["rhad33_max"].as<float>();
    float tight_rhad33_min = configYaml["analysis"]["tight"]["rhad33_min"].as<float>();

    float tight_w72_max = configYaml["analysis"]["tight"]["w72_max"].as<float>();
    float tight_w72_min = configYaml["analysis"]["tight"]["w72_min"].as<float>();

    float tight_re11_E_max = configYaml["analysis"]["tight"]["re11_E_max"].as<float>();
    float tight_re11_E_min = configYaml["analysis"]["tight"]["re11_E_min"].as<float>();

    float tight_CNN_min = configYaml["analysis"]["tight"]["CNN_min"].as<float>();
    float tight_CNN_max = configYaml["analysis"]["tight"]["CNN_max"].as<float>();

    float tight_weta_cogx_max = configYaml["analysis"]["tight"]["weta_cogx_max"].as<float>();
    float tight_weta_cogx_min = configYaml["analysis"]["tight"]["weta_cogx_min"].as<float>();
    float tight_weta_cogx_max_b = configYaml["analysis"]["tight"]["weta_cogx_max_b"].as<float>();
    float tight_weta_cogx_max_s = configYaml["analysis"]["tight"]["weta_cogx_max_s"].as<float>();

    float tight_wphi_cogx_max = configYaml["analysis"]["tight"]["wphi_cogx_max"].as<float>();
    float tight_wphi_cogx_min = configYaml["analysis"]["tight"]["wphi_cogx_min"].as<float>();
    float tight_wphi_cogx_max_b = configYaml["analysis"]["tight"]["wphi_cogx_max_b"].as<float>();
    float tight_wphi_cogx_max_s = configYaml["analysis"]["tight"]["wphi_cogx_max_s"].as<float>();

    float tight_e11_over_e33_max = configYaml["analysis"]["tight"]["e11_over_e33_max"].as<float>();
    float tight_e11_over_e33_min = configYaml["analysis"]["tight"]["e11_over_e33_min"].as<float>();

    float tight_et1_max = configYaml["analysis"]["tight"]["et1_max"].as<float>();
    float tight_et1_min = configYaml["analysis"]["tight"]["et1_min"].as<float>();
    float tight_et1_min_b = configYaml["analysis"]["tight"]["et1_min_b"].as<float>();
    float tight_et1_min_s = configYaml["analysis"]["tight"]["et1_min_s"].as<float>();

    float tight_et2_max = configYaml["analysis"]["tight"]["et2_max"].as<float>(1.0);
    float tight_et2_min = configYaml["analysis"]["tight"]["et2_min"].as<float>(0.0);

    float tight_et3_max = configYaml["analysis"]["tight"]["et3_max"].as<float>(1.0);
    float tight_et3_min = configYaml["analysis"]["tight"]["et3_min"].as<float>(0.0);

    float tight_e32_over_e35_max = configYaml["analysis"]["tight"]["e32_over_e35_max"].as<float>();
    float tight_e32_over_e35_min = configYaml["analysis"]["tight"]["e32_over_e35_min"].as<float>();

    float tight_prob_max = configYaml["analysis"]["tight"]["prob_max"].as<float>();
    float tight_prob_min = configYaml["analysis"]["tight"]["prob_min"].as<float>();

    float tight_et4_max = configYaml["analysis"]["tight"]["et4_max"].as<float>();
    float tight_et4_min = configYaml["analysis"]["tight"]["et4_min"].as<float>();

    float tight_w32_max = configYaml["analysis"]["tight"]["w32_max"].as<float>();
    float tight_w32_min = configYaml["analysis"]["tight"]["w32_min"].as<float>();

    float tight_bdt_max = configYaml["analysis"]["tight"]["bdt_max"].as<float>(1);
    float tight_bdt_min = configYaml["analysis"]["tight"]["bdt_min"].as<float>(0);
    float tight_bdt_min_slope = configYaml["analysis"]["tight"]["bdt_min_slope"].as<float>(0);
    float tight_bdt_min_intercept = configYaml["analysis"]["tight"]["bdt_min_intercept"].as<float>(tight_bdt_min);

    // non tight cuts
    std::cout << "non tight cuts" << std::endl;
    float non_tight_reta77_min = configYaml["analysis"]["non_tight"]["reta77_min"].as<float>();
    float non_tight_reta77_max = configYaml["analysis"]["non_tight"]["reta77_max"].as<float>();

    float non_tight_rhad33_max = configYaml["analysis"]["non_tight"]["rhad33_max"].as<float>();
    float non_tight_rhad33_min = configYaml["analysis"]["non_tight"]["rhad33_min"].as<float>();

    float non_tight_w72_max = configYaml["analysis"]["non_tight"]["w72_max"].as<float>();
    float non_tight_w72_min = configYaml["analysis"]["non_tight"]["w72_min"].as<float>();

    float non_tight_re11_E_max = configYaml["analysis"]["non_tight"]["re11_E_max"].as<float>();
    float non_tight_re11_E_min = configYaml["analysis"]["non_tight"]["re11_E_min"].as<float>();

    float non_tight_CNN_min = configYaml["analysis"]["non_tight"]["CNN_min"].as<float>();
    float non_tight_CNN_max = configYaml["analysis"]["non_tight"]["CNN_max"].as<float>();

    float non_tight_weta_cogx_max = configYaml["analysis"]["non_tight"]["weta_cogx_max"].as<float>();
    float non_tight_weta_cogx_min = configYaml["analysis"]["non_tight"]["weta_cogx_min"].as<float>();
    float non_tight_weta_cogx_max_b = configYaml["analysis"]["non_tight"]["weta_cogx_max_b"].as<float>();
    float non_tight_weta_cogx_max_s = configYaml["analysis"]["non_tight"]["weta_cogx_max_s"].as<float>();

    float non_tight_wphi_cogx_max = configYaml["analysis"]["non_tight"]["wphi_cogx_max"].as<float>();
    float non_tight_wphi_cogx_min = configYaml["analysis"]["non_tight"]["wphi_cogx_min"].as<float>();
    float non_tight_wphi_cogx_max_b = configYaml["analysis"]["non_tight"]["wphi_cogx_max_b"].as<float>();
    float non_tight_wphi_cogx_max_s = configYaml["analysis"]["non_tight"]["wphi_cogx_max_s"].as<float>();

    float non_tight_prob_max = configYaml["analysis"]["non_tight"]["prob_max"].as<float>();
    float non_tight_prob_min = configYaml["analysis"]["non_tight"]["prob_min"].as<float>();

    float non_tight_et1_max = configYaml["analysis"]["non_tight"]["et1_max"].as<float>();
    float non_tight_et1_min = configYaml["analysis"]["non_tight"]["et1_min"].as<float>();

    float non_tight_e11_over_e33_max = configYaml["analysis"]["non_tight"]["e11_over_e33_max"].as<float>();
    float non_tight_e11_over_e33_min = configYaml["analysis"]["non_tight"]["e11_over_e33_min"].as<float>();

    float non_tight_e32_over_e35_max = configYaml["analysis"]["non_tight"]["e32_over_e35_max"].as<float>();
    float non_tight_e32_over_e35_min = configYaml["analysis"]["non_tight"]["e32_over_e35_min"].as<float>();

    float non_tight_et4_max = configYaml["analysis"]["non_tight"]["et4_max"].as<float>();
    float non_tight_et4_min = configYaml["analysis"]["non_tight"]["et4_min"].as<float>();

    float non_tight_w32_max = configYaml["analysis"]["non_tight"]["w32_max"].as<float>();
    float non_tight_w32_min = configYaml["analysis"]["non_tight"]["w32_min"].as<float>();

    float non_tight_bdt_max = configYaml["analysis"]["non_tight"]["bdt_max"].as<float>(1);
    float non_tight_bdt_min = configYaml["analysis"]["non_tight"]["bdt_min"].as<float>(0);
    float non_tight_bdt_max_slope = configYaml["analysis"]["non_tight"]["bdt_max_slope"].as<float>(0);
    float non_tight_bdt_max_intercept = configYaml["analysis"]["non_tight"]["bdt_max_intercept"].as<float>(non_tight_bdt_max);

    // common cuts for both tight and non tight

    float common_prob_max = configYaml["analysis"]["common"]["prob_max"].as<float>();
    float common_prob_min = configYaml["analysis"]["common"]["prob_min"].as<float>();

    float common_e11_over_e33_max = configYaml["analysis"]["common"]["e11_over_e33_max"].as<float>();
    float common_e11_over_e33_min = configYaml["analysis"]["common"]["e11_over_e33_min"].as<float>();

    float common_wr_cogx_bound = configYaml["analysis"]["common"]["wr_cogx_bound"].as<float>();
    float common_cluster_weta_cogx_bound = configYaml["analysis"]["common"]["cluster_weta_cogx_bound"].as<float>();

    int common_npb_cut_on = configYaml["analysis"]["common"]["npb_cut_on"].as<int>(0);
    float common_npb_score_cut = configYaml["analysis"]["common"]["npb_score_cut"].as<float>(0.5);

    int reweight = configYaml["analysis"]["unfold"]["reweight"].as<int>(); // 0 for no reweighting, 1 for reweighting

    float clusterescale = configYaml["analysis"]["cluster_escale"].as<float>(1.0);
    float clustereres = configYaml["analysis"]["cluster_eres"].as<float>(0.0);

    // polynomial 3 for the reweighting for response matrix
    TF1 *f_reweight = new TF1("f_reweight", "([0] + [1]*x + [3]*x*x) / (1 + [2]*x + [4]*x*x)", 0, 100);
    //f_reweight->SetParameters(0.714962, -0.0856443, -0.125383, 0.00345831, 0.00462972);
    f_reweight->SetParameters(1, 0, 0, 0, 0);

    /////////////////////////////////////
    // TTreeReader 
    /////////////////////////////////////   
    TTreeReader reader(&chain);
    
    // Basic event variables
    TTreeReaderValue<int> pythiaid(reader, "pythiaid");
    TTreeReaderValue<int> nparticles(reader, "nparticles");
    TTreeReaderValue<int> ncluster(reader, Form("ncluster_%s", clusternodename.c_str()));
    TTreeReaderValue<int> runnumber(reader, "runnumber");
    TTreeReaderArray<Bool_t> scaledtrigger(reader, "scaledtrigger");
    TTreeReaderArray<Bool_t> livetrigger(reader, "livetrigger");
    TTreeReaderValue<float> energy_scale(reader, "energy_scale");
    TTreeReaderValue<float> vertexz(reader, "vertexz");
    TTreeReaderValue<float> vertexz_truth(reader, "vertexz_truth");
    // Centrality branch only for Au+Au; pp has no cent branch and uses a single bin (centbin=0)
    std::unique_ptr<TTreeReaderValue<float>> cent_reader;
    if (!is_pp && chain.FindBranch("cent"))
      cent_reader = std::make_unique<TTreeReaderValue<float>>(reader, "cent");
    else if (!is_pp && !chain.FindBranch("cent"))
      std::cout << "RecoEffCalculator_TTreeReader: WARNING data_aa but no 'cent' branch in tree; all events will be skipped (centbin invalid)" << std::endl;
    TTreeReaderValue<float> totalEMCal_energy(reader, "totalEMCal_energy");
    TTreeReaderValue<float> totalIHCal_energy(reader, "totalIHCal_energy");
    TTreeReaderValue<float> totalOHCal_energy(reader, "totalOHCal_energy");
    TTreeReaderArray<float> trigger_prescale(reader, "trigger_prescale");

    // Particle arrays
    TTreeReaderArray<float> particle_E(reader, "particle_E");
    TTreeReaderArray<float> particle_Pt(reader, "particle_Pt");
    TTreeReaderArray<float> particle_Eta(reader, "particle_Eta");
    TTreeReaderArray<float> particle_Phi(reader, "particle_Phi");
    TTreeReaderArray<float> particle_truth_iso_02(reader, "particle_truth_iso_02");
    TTreeReaderArray<float> particle_truth_iso_03(reader, "particle_truth_iso_03");
    TTreeReaderArray<float> particle_truth_iso_04(reader, "particle_truth_iso_04");
    TTreeReaderArray<int> particle_pid(reader, "particle_pid");
    TTreeReaderArray<int> particle_trkid(reader, "particle_trkid");
    TTreeReaderArray<int> particle_photonclass(reader, "particle_photonclass");
    TTreeReaderArray<int> particle_converted(reader, "particle_converted");
    // Cluster arrays
    TTreeReaderArray<float> cluster_E(reader, Form("cluster_E_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_Et(reader, Form("cluster_Et_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_Eta(reader, Form("cluster_Eta_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_Phi(reader, Form("cluster_Phi_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_prob(reader, Form("cluster_prob_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_CNN_prob(reader, Form("cluster_CNN_prob_%s", clusternodename.c_str()));
    TTreeReaderArray<int> cluster_truthtrkID(reader, Form("cluster_truthtrkID_%s", clusternodename.c_str()));
    TTreeReaderArray<int> cluster_pid(reader, Form("cluster_pid_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_02(reader, Form("cluster_iso_02_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03(reader, Form("cluster_iso_03_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_04(reader, Form("cluster_iso_04_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e1(reader, Form("cluster_e1_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e2(reader, Form("cluster_e2_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e3(reader, Form("cluster_e3_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e4(reader, Form("cluster_e4_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_et1(reader, Form("cluster_et1_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_et2(reader, Form("cluster_et2_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_et3(reader, Form("cluster_et3_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_et4(reader, Form("cluster_et4_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_weta(reader, Form("cluster_weta_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_wphi(reader, Form("cluster_wphi_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_ietacent(reader, Form("cluster_ietacent_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iphicent(reader, Form("cluster_iphicent_%s", clusternodename.c_str()));
    TTreeReaderArray<int> cluster_detamax(reader, Form("cluster_detamax_%s", clusternodename.c_str()));
    TTreeReaderArray<int> cluster_dphimax(reader, Form("cluster_dphimax_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_weta_cogx(reader, Form("cluster_weta_cogx_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_wphi_cogx(reader, Form("cluster_wphi_cogx_%s", clusternodename.c_str()));
    TTreeReaderArray<int> cluster_nsaturated(reader, Form("cluster_nsaturated_%s", clusternodename.c_str()));

    // Cluster energy arrays
    TTreeReaderArray<float> cluster_e11(reader, Form("cluster_e11_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e22(reader, Form("cluster_e22_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e13(reader, Form("cluster_e13_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e15(reader, Form("cluster_e15_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e17(reader, Form("cluster_e17_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e31(reader, Form("cluster_e31_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e51(reader, Form("cluster_e51_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e71(reader, Form("cluster_e71_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e33(reader, Form("cluster_e33_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e35(reader, Form("cluster_e35_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e37(reader, Form("cluster_e37_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e53(reader, Form("cluster_e53_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e73(reader, Form("cluster_e73_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e55(reader, Form("cluster_e55_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e57(reader, Form("cluster_e57_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e75(reader, Form("cluster_e75_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e77(reader, Form("cluster_e77_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_w32(reader, Form("cluster_w32_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e32(reader, Form("cluster_e32_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_w72(reader, Form("cluster_w72_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e72(reader, Form("cluster_e72_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_w52(reader, Form("cluster_w52_%s", clusternodename.c_str()));

    // Cluster arrays for 2D data - using regular TTreeReaderArray for 2D arrays
    TTreeReaderArray<float> cluster_e_array(reader, Form("cluster_e_array_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_adc_array(reader, Form("cluster_adc_array_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_time_array(reader, Form("cluster_time_array_%s", clusternodename.c_str()));
    TTreeReaderArray<int> cluster_status_array(reader, Form("cluster_status_array_%s", clusternodename.c_str()));
    TTreeReaderArray<int> cluster_ownership_array(reader, Form("cluster_ownership_array_%s", clusternodename.c_str()));
 
    // Cluster isolation arrays
    TTreeReaderArray<float> cluster_iso_03_emcal(reader, Form("cluster_iso_03_emcal_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_hcalin(reader, Form("cluster_iso_03_hcalin_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_hcalout(reader, Form("cluster_iso_03_hcalout_%s", clusternodename.c_str()));

    TTreeReaderArray<float> cluster_iso_02_70_emcal(reader, Form("cluster_iso_02_70_emcal_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_005_70_emcal(reader, Form("cluster_iso_005_70_emcal_%s", clusternodename.c_str()));
  
    TTreeReaderArray<float> cluster_iso_03_70_emcal(reader, Form("cluster_iso_04_sub1_emcal_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_70_hcalin(reader, Form("cluster_iso_04_sub1_hcalin_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_70_hcalout(reader, Form("cluster_iso_04_sub1_hcalout_%s", clusternodename.c_str()));
  
    //TTreeReaderArray<float> cluster_iso_03_70_emcal(reader, Form("cluster_iso_03_70_emcal_%s", clusternodename.c_str()));
    //TTreeReaderArray<float> cluster_iso_03_70_hcalin(reader, Form("cluster_iso_03_70_hcalin_%s", clusternodename.c_str()));
    //TTreeReaderArray<float> cluster_iso_03_70_hcalout(reader, Form("cluster_iso_03_70_hcalout_%s", clusternodename.c_str()));


    TTreeReaderArray<float> cluster_iso_01_70_emcal(reader, Form("cluster_iso_01_70_emcal_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_120_emcal(reader, Form("cluster_iso_03_120_emcal_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_120_hcalin(reader, Form("cluster_iso_03_120_hcalin_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_120_hcalout(reader, Form("cluster_iso_03_120_hcalout_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_60_emcal(reader, Form("cluster_iso_03_60_emcal_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_60_hcalin(reader, Form("cluster_iso_03_60_hcalin_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_60_hcalout(reader, Form("cluster_iso_03_60_hcalout_%s", clusternodename.c_str()));

    // Cluster HCal arrays
    TTreeReaderArray<float> cluster_ihcal_et(reader, Form("cluster_ihcal_et_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_ohcal_et(reader, Form("cluster_ohcal_et_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_ihcal_et22(reader, Form("cluster_ihcal_et22_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_ohcal_et22(reader, Form("cluster_ohcal_et22_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_ihcal_et33(reader, Form("cluster_ihcal_et33_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_ohcal_et33(reader, Form("cluster_ohcal_et33_%s", clusternodename.c_str()));
    TTreeReaderArray<int> cluster_ihcal_ieta(reader, Form("cluster_ihcal_ieta_%s", clusternodename.c_str()));
    TTreeReaderArray<int> cluster_ihcal_iphi(reader, Form("cluster_ihcal_iphi_%s", clusternodename.c_str()));
    TTreeReaderArray<int> cluster_ohcal_ieta(reader, Form("cluster_ohcal_ieta_%s", clusternodename.c_str()));
    TTreeReaderArray<int> cluster_ohcal_iphi(reader, Form("cluster_ohcal_iphi_%s", clusternodename.c_str()));

    // BDT score — build one TTreeReaderArray per unique model name
    std::vector<std::string> all_bdt_models = {bdt_model_name};
    if (use_et_binned_bdt)
        all_bdt_models.insert(all_bdt_models.end(), bdt_et_bin_models.begin(), bdt_et_bin_models.end());
    std::sort(all_bdt_models.begin(), all_bdt_models.end());
    all_bdt_models.erase(std::unique(all_bdt_models.begin(), all_bdt_models.end()), all_bdt_models.end());

    std::map<std::string, TTreeReaderArray<float>*> bdt_arrays;
    for (auto& mname : all_bdt_models) {
        bdt_arrays[mname] = new TTreeReaderArray<float>(reader,
            Form("cluster_bdt_%s_%s", clusternodename.c_str(), mname.c_str()));
    }

    // NPB score
    TTreeReaderArray<float> cluster_npb_score(reader, Form("cluster_npb_score_%s", clusternodename.c_str()));

    // Truth jet arrays (Anti-kT truth; radius matches analysis.jet_cone_size, same as reco tower jets)
    TTreeReaderValue<int> njet_truth_R3(reader, "njet_truth_AntiKt_Truth_r03");
    TTreeReaderArray<float> jet_truth_E_R3(reader, "jet_truth_E_AntiKt_Truth_r03");
    TTreeReaderArray<float> jet_truth_Et_R3(reader, "jet_truth_Et_AntiKt_Truth_r03");
    TTreeReaderArray<float> jet_truth_Pt_R3(reader, "jet_truth_Pt_AntiKt_Truth_r03");
    TTreeReaderArray<float> jet_truth_Eta_R3(reader, "jet_truth_Eta_AntiKt_Truth_r03");
    TTreeReaderArray<float> jet_truth_Phi_R3(reader, "jet_truth_Phi_AntiKt_Truth_r03");

    TTreeReaderValue<int> njet_truth_R4(reader, "njet_truth_AntiKt_Truth_r04");
    TTreeReaderArray<float> jet_truth_E_R4(reader, "jet_truth_E_AntiKt_Truth_r04");
    TTreeReaderArray<float> jet_truth_Et_R4(reader, "jet_truth_Et_AntiKt_Truth_r04");
    TTreeReaderArray<float> jet_truth_Pt_R4(reader, "jet_truth_Pt_AntiKt_Truth_r04");
    TTreeReaderArray<float> jet_truth_Eta_R4(reader, "jet_truth_Eta_AntiKt_Truth_r04");
    TTreeReaderArray<float> jet_truth_Phi_R4(reader, "jet_truth_Phi_AntiKt_Truth_r04");

    TTreeReaderValue<int> &njet_truth = use_jet_r04 ? njet_truth_R4 : njet_truth_R3;
    TTreeReaderArray<float> &jet_truth_E = use_jet_r04 ? jet_truth_E_R4 : jet_truth_E_R3;
    TTreeReaderArray<float> &jet_truth_Et = use_jet_r04 ? jet_truth_Et_R4 : jet_truth_Et_R3;
    TTreeReaderArray<float> &jet_truth_Pt = use_jet_r04 ? jet_truth_Pt_R4 : jet_truth_Pt_R3;
    TTreeReaderArray<float> &jet_truth_Eta = use_jet_r04 ? jet_truth_Eta_R4 : jet_truth_Eta_R3;
    TTreeReaderArray<float> &jet_truth_Phi = use_jet_r04 ? jet_truth_Phi_R4 : jet_truth_Phi_R3;

    // Reco jet arrays (Anti-kT tower; radius from analysis.jet_cone_size)
    TTreeReaderValue<int> njet_R3(reader, "njet_AntiKt_Tower_r03_Sub1");
    TTreeReaderArray<float> jet_Pt_R3(reader, "jet_Pt_AntiKt_Tower_r03_Sub1");
    TTreeReaderArray<float> jet_Eta_R3(reader, "jet_Eta_AntiKt_Tower_r03_Sub1");
    TTreeReaderArray<float> jet_Phi_R3(reader, "jet_Phi_AntiKt_Tower_r03_Sub1");

    TTreeReaderValue<int> njet_R4(reader, "njet_AntiKt_Tower_r04_Sub1");
    TTreeReaderArray<float> jet_Pt_R4(reader, "jet_Pt_AntiKt_Tower_r04_Sub1");
    TTreeReaderArray<float> jet_Eta_R4(reader, "jet_Eta_AntiKt_Tower_r04_Sub1");
    TTreeReaderArray<float> jet_Phi_R4(reader, "jet_Phi_AntiKt_Tower_r04_Sub1");

    TTreeReaderValue<int> &njet = use_jet_r04 ? njet_R4 : njet_R3;
    TTreeReaderArray<float> &jet_Pt = use_jet_r04 ? jet_Pt_R4 : jet_Pt_R3;
    TTreeReaderArray<float> &jet_Eta = use_jet_r04 ? jet_Eta_R4 : jet_Eta_R3;
    TTreeReaderArray<float> &jet_Phi = use_jet_r04 ? jet_Phi_R4 : jet_Phi_R3;



    /////////////////////////////////////
    // Output file
    /////////////////////////////////////
    TFile *fout = new TFile(outfilename.c_str(), "RECREATE");
    TH1::SetDefaultSumw2(kTRUE);
    TH1F *h_max_photon_pT = new TH1F("h_max_photon_pT", "Max Photon pT", 1000, 0, 100);
    TH1F *h_photon_pT = new TH1F("h_photon_pT", "Photon pT", 1000, 0, 100);
    TH1F *h_max_direct_pT = new TH1F("h_max_direct_pT", "Max Direct Photon pT", 1000, 0, 100);
    TH1F *h_direct_pT = new TH1F("h_direct_pT", "Direct Photon pT", 1000, 0, 100);
    TH1F *h_max_frag_pT = new TH1F("h_max_frag_pT", "Max Fragmentation Photon pT", 1000, 0, 100);
    TH1F *h_frag_pT = new TH1F("h_frag_pT", "Fragmentation Photon pT", 1000, 0, 100);
    TH1F *h_max_decay_pT = new TH1F("h_max_decay_pT", "Max Decay Photon pT", 1000, 0, 100);
    TH1F *h_decay_photon_pT = new TH1F("h_decay_photon_pT", "Decay Photon pT", 1000, 0, 100);
    TH1F *h_vertexz = new TH1F("h_vertexz", "Vertex z", 200, -100, 100);
    TH1D *h_totalEMCal_energy_tight_weight = new TH1D(
        "h_totalEMCal_energy_tight_weight",
        "Total EMCal energy (events with at least one tight cluster);E_{tot}^{EMCal} [GeV];events",
        200, 0., 2e3);
    TH1D *h_totalIHCal_energy_tight_weight = new TH1D(
        "h_totalIHCal_energy_tight_weight",
        "Total IHCal energy (events with at least one tight cluster);E_{tot}^{IHCal} [GeV];events",
        200, 0., 2e3);
    TH1D *h_totalOHCal_energy_tight_weight = new TH1D(
        "h_totalOHCal_energy_tight_weight",
        "Total OHCal energy (events with at least one tight cluster);E_{tot}^{OHCal} [GeV];events",
        200, 0., 2e3);
    std::vector<double> cent_hist_edges(centrality_bins.begin(), centrality_bins.end());
    TH1D *h_centrality_tight_weight = new TH1D(
        "h_centrality_tight_weight",
        ";Centrality [%];events",
        100,0,100);
    TH1F *h_cluster_common_Et = new TH1F("h_cluster_common_E", "Cluster Common E", 1000, 0, 100);
    TH1F *h_cluster_common_leading_Et = new TH1F("h_cluster_common_leading_E", "Cluster Common Leading E", 1000, 0, 100);

    TH1F *h_max_truth_jet_pT = new TH1F("h_max_truth_jet_pT", "Max Truth Jet pT", 1000, 0, 100);

    TH1F *h_max_photon_pT_vertexcut = new TH1F("h_max_photon_pT_vertexcut", "Max Photon pT Vertex Cut", 1000, 0, 100);


    // TEfficiency for conversion and reco
    TEfficiency::EStatOption effopt = TEfficiency::kBUniform;
    TEfficiency *eff_reco = new TEfficiency("eff_reco", "Reco Efficiency", 40, 10, 50, 50, -1, 1);
    eff_reco->SetStatisticOption(effopt);
    TEfficiency *eff_id = new TEfficiency("eff_id", "ID Efficiency", 40, 10, 50, 50, -1, 1);
    eff_id->SetStatisticOption(effopt);
    TEfficiency *eff_converts = new TEfficiency("eff_converts", "Conversion Prob", 40, 10, 50, 50, -1, 1);
    eff_converts->SetStatisticOption(effopt);

    // Per-centrality TEfficiencies 
    std::vector<TEfficiency *> eff_reco_cent;      
    std::vector<TEfficiency *> eff_iso_cent;       
    std::vector<TEfficiency *> eff_id_cent;        
    std::vector<TEfficiency *> eff_converts_cent;  
    std::vector<TEfficiency *> eff_all_cent;       

    // truth pythia
    std::vector<TH1D *> h_truth_pT;

    std::vector<TH1D *> h_truth_pT_vertexcut;


    std::vector<TH1D *> h_tight_iso_cluster_signal;
    std::vector<TH1D *> h_tight_noniso_cluster_signal;
    std::vector<TH1D *> h_nontight_iso_cluster_signal;
    std::vector<TH1D *> h_nontight_noniso_cluster_signal;
    std::vector<TH1D *> h_tight_iso_cluster_notmatch;
    std::vector<TH1D *> h_tight_noniso_cluster_notmatch;
    std::vector<TH1D *> h_nontight_iso_cluster_notmatch;
    std::vector<TH1D *> h_nontight_noniso_cluster_notmatch;
    std::vector<TH1D *> h_all_cluster_signal;
    std::vector<TH2D *> h_all_cluster_Et_max_b2bjet; // max cluster Et vs max backtobjets Et
    std::vector<TH1D *> h_tight_cluster_signal;

    std::vector<TH1D *> h_tight_noniso_cluster_background;
    std::vector<TH1D *> h_nontight_iso_cluster_background;
    std::vector<TH1D *> h_nontight_noniso_cluster_background;

    std::vector<TH2D *> h_singal_reco_isoET;
    std::vector<TH2D *> h_singal_truth_isoET;
    std::vector<TH2D *> h_background_truth_isoET;

    // here are for the plots we gonna make for both data and simulation
    std::vector<TH1D *> h_tight_iso_cluster;
    std::vector<TH1D *> h_tight_noniso_cluster;
    std::vector<TH1D *> h_nontight_iso_cluster;
    std::vector<TH1D *> h_nontight_noniso_cluster;
    std::vector<TH1D *> h_common_cluster;
    std::vector<TH1D *> h_all_cluster;
    std::vector<TH1D *> h_tight_cluster;


    // unfold response matrix
    std::vector<RooUnfoldResponse *> responses_full;
    std::vector<RooUnfoldResponse *> responses_half;
    // vector for the response matrix th2
    std::vector<TH2D *> h_response_full_list;
    std::vector<TH2D *> h_response_half_list;
    // Flattened 2D response matrix for unfolding in (pT^gamma, xjgamma):
    // x-axis = truth global bin, y-axis = reco global bin
    std::vector<TH2D *> h_response_xjgamma_global_list;
    // Native 2D RooUnfold response in (xjgamma, pT^gamma)
    std::vector<TH2D *> h_xjgamma_truth_response;
    std::vector<TH2D *> h_xjgamma_reco_response;
    std::vector<RooUnfoldResponse *> responses_xjgamma_2d;
    // id histogram for unfolding
    std::vector<TH1D *> h_pT_truth_response;
    std::vector<TH1D *> h_pT_reco_response;

    std::vector<TH1D *> h_pT_reco_fake;

    std::vector<TH1D *> h_pT_truth_half_response;
    std::vector<TH1D *> h_pT_reco_half_response;

    std::vector<TH1D *> h_pT_truth_secondhalf_response;
    std::vector<TH1D *> h_pT_reco_secondhalf_response;

    // direct and fragmentation photon pT vs truth iso ET
    std::vector<TH2D *> h_direct_pT_truth_isoET;
    std::vector<TH2D *> h_frag_pT_truth_isoET;

    // n cluster per photon
    std::vector<TH2D *> h_ncluster_truth;

    // energy resolution
    std::vector<TH2D *> h_pT_truth_reco;
    // reco-jet / truth-jet response vs truth-jet pT (matched with DeltaR < 0.2)
    std::vector<TH2D *> h_jet_pT_response;

    // xjgamma, why not?
    // this is for sim only
    std::vector<TH2D *> h_tight_iso_xjgamma_signal;
    std::vector<TH2D *> h_tight_noniso_xjgamma_signal;
    std::vector<TH2D *> h_nontight_iso_xjgamma_signal;
    std::vector<TH2D *> h_nontight_noniso_xjgamma_signal;
    std::vector<TH2D *> h_all_xjgamma_signal;
    std::vector<TH2D *> h_tight_xjgamma_signal;
    // using reco jets that are matched to truth jets (DeltaR < 0.2), sim only
    std::vector<TH2D *> h_tight_iso_truthmatchreco_xjgamma_signal;
    std::vector<TH2D *> h_tight_noniso_truthmatchreco_xjgamma_signal;
    std::vector<TH2D *> h_nontight_iso_truthmatchreco_xjgamma_signal;
    std::vector<TH2D *> h_nontight_noniso_truthmatchreco_xjgamma_signal;
    std::vector<TH2D *> h_all_truthmatchreco_xjgamma_signal;
    std::vector<TH2D *> h_tight_truthmatchreco_xjgamma_signal;

    // using truth jet
    std::vector<TH2D *> h_tight_iso_truthjet_xjgamma_signal;
    std::vector<TH2D *> h_tight_noniso_truthjet_xjgamma_signal;
    std::vector<TH2D *> h_nontight_iso_truthjet_xjgamma_signal;
    std::vector<TH2D *> h_nontight_noniso_truthjet_xjgamma_signal;
    std::vector<TH2D *> h_all_truthjet_xjgamma_signal;
    std::vector<TH2D *> h_tight_truthjet_xjgamma_signal;

    std::vector<TH2D *> h_tight_iso_xjgamma_background;
    std::vector<TH2D *> h_tight_noniso_xjgamma_background;
    std::vector<TH2D *> h_nontight_iso_xjgamma_background;
    std::vector<TH2D *> h_nontight_noniso_xjgamma_background;

    // for both data and sim
    std::vector<TH2D *> h_tight_iso_xjgamma;
    std::vector<TH2D *> h_tight_noniso_xjgamma;
    std::vector<TH2D *> h_nontight_iso_xjgamma;
    std::vector<TH2D *> h_nontight_noniso_xjgamma;
    std::vector<TH2D *> h_common_xjgamma;
    std::vector<TH2D *> h_all_xjgamma;
    std::vector<TH2D *> h_tight_xjgamma;

    // isolation profile for debugging reasons
    std::vector<std::vector<TH1D *>> h_tight_cluster_pT;
    h_tight_cluster_pT.resize(n_cent_bins);
    std::vector<std::vector<TH1D *>> h_nontight_cluster_pT;
    h_nontight_cluster_pT.resize(n_cent_bins);
    // truth iso vs reco iso for different pT bins and eta bins
    std::vector<std::vector<TH2D *>> h_iso_truth_reco;
    h_iso_truth_reco.resize(n_cent_bins);
    std::vector<std::vector<TH2D *>> h_background_iso_truth_reco;
    h_background_iso_truth_reco.resize(n_cent_bins);
    // response vs. isoET
    std::vector<std::vector<TH2D *>> h_response_isoET;
    h_response_isoET.resize(n_cent_bins);
    // |DeltaPhi|(cluster, jet), tight clusters only — same kinematics as EvtCharacter.C dphi block
    std::vector<std::vector<TH1D *>> h_dphi_clusterJets_tight;
    h_dphi_clusterJets_tight.resize(n_cent_bins);

    // Standard efficiencies and spectra binned by centrality (eta only as acceptance cut)
    for (int icent = 0; icent < n_cent_bins; icent++)
    {   
        // efficiencies
        eff_reco_cent.push_back(new TEfficiency(Form("eff_reco_cent%d", icent), Form("Reco Efficiency %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_pT_bins_truth, pT_bin_edges_truth));
        eff_reco_cent[icent]->SetStatisticOption(effopt);
        eff_reco_cent[icent]->SetWeight(weight);
        eff_id_cent.push_back(new TEfficiency(Form("eff_id_cent%d", icent), Form("ID Efficiency %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_pT_bins_truth, pT_bin_edges_truth));
        eff_id_cent[icent]->SetStatisticOption(effopt);
        eff_id_cent[icent]->SetWeight(weight);
        eff_converts_cent.push_back(new TEfficiency(Form("eff_converts_cent%d", icent), Form("Conversion Prob %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_pT_bins_truth, pT_bin_edges_truth));
        eff_converts_cent[icent]->SetStatisticOption(effopt);
        eff_converts_cent[icent]->SetWeight(weight);
        eff_iso_cent.push_back(new TEfficiency(Form("eff_iso_cent%d", icent), Form("Iso Efficiency %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_pT_bins_truth, pT_bin_edges_truth));
        eff_iso_cent[icent]->SetStatisticOption(effopt);
        eff_iso_cent[icent]->SetWeight(weight);
        eff_all_cent.push_back(new TEfficiency(Form("eff_all_cent%d", icent), Form("All Efficiency %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_pT_bins_truth, pT_bin_edges_truth));
        eff_all_cent[icent]->SetStatisticOption(effopt);
        eff_all_cent[icent]->SetWeight(weight);
        h_truth_pT.push_back(new TH1D(Form("h_truth_pT_%d", icent), Form("Truth pT %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_pT_bins_truth, pT_bin_edges_truth));
        // vertex eff
        h_truth_pT_vertexcut.push_back(new TH1D(Form("h_truth_pT_vertexcut_%d", icent), Form("Truth pT %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_pT_bins_truth, pT_bin_edges_truth));
        h_tight_iso_cluster_signal.push_back(new TH1D(Form("h_tight_iso_cluster_signal_%d", icent), Form("Tight Iso Cluster %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_pT_bins, pT_bin_edges));
        h_tight_noniso_cluster_signal.push_back(new TH1D(Form("h_tight_noniso_cluster_signal_%d", icent), Form("Tight Non-Iso Cluster %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_pT_bins, pT_bin_edges));
        h_nontight_iso_cluster_signal.push_back(new TH1D(Form("h_nontight_iso_cluster_signal_%d", icent), Form("Non-Tight Iso Cluster %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_pT_bins, pT_bin_edges));
        h_nontight_noniso_cluster_signal.push_back(new TH1D(Form("h_nontight_noniso_cluster_signal_%d", icent), Form("Non-Tight Non-Iso Cluster %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_pT_bins, pT_bin_edges));
        h_tight_iso_cluster_notmatch.push_back(new TH1D(Form("h_tight_iso_cluster_notmatch_%d", icent), Form("Tight Iso Not Matched %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_pT_bins, pT_bin_edges));
        h_tight_noniso_cluster_notmatch.push_back(new TH1D(Form("h_tight_noniso_cluster_notmatch_%d", icent), Form("Tight Non-Iso Not Matched %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_pT_bins, pT_bin_edges));
        h_nontight_iso_cluster_notmatch.push_back(new TH1D(Form("h_nontight_iso_cluster_notmatch_%d", icent), Form("Non-Tight Iso Not Matched %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_pT_bins, pT_bin_edges));
        h_nontight_noniso_cluster_notmatch.push_back(new TH1D(Form("h_nontight_noniso_cluster_notmatch_%d", icent), Form("Non-Tight Non-Iso Not Matched %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_pT_bins, pT_bin_edges));
        h_all_cluster_signal.push_back(new TH1D(Form("h_all_cluster_signal_%d", icent), Form("All Cluster %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_pT_bins, pT_bin_edges));
        h_all_cluster_Et_max_b2bjet.push_back(new TH2D(Form("h_all_cluster_Et_max_b2bjet_%d", icent), Form("All Cluster Et Max Backtobjets %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_pT_bins, pT_bin_edges, 100, 0, 100));
        h_tight_cluster_signal.push_back(new TH1D(Form("h_tight_cluster_signal_%d", icent), Form("Tight Cluster %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_pT_bins, pT_bin_edges));
        h_tight_noniso_cluster_background.push_back(new TH1D(Form("h_tight_noniso_cluster_background_%d", icent), Form("Tight Non-Iso Cluster %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_pT_bins, pT_bin_edges));
        h_nontight_iso_cluster_background.push_back(new TH1D(Form("h_nontight_iso_cluster_background_%d", icent), Form("Non-Tight Iso Cluster %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_pT_bins, pT_bin_edges));
        h_nontight_noniso_cluster_background.push_back(new TH1D(Form("h_nontight_noniso_cluster_background_%d", icent), Form("Non-Tight Non-Iso Cluster %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_pT_bins, pT_bin_edges));
        h_singal_reco_isoET.push_back(new TH2D(Form("h_singal_reco_isoET_%d", icent), Form("Signal Reco Iso ET %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), 400, 0, 50, 4400, -5, 50));
        h_singal_truth_isoET.push_back(new TH2D(Form("h_singal_truth_isoET_%d", icent), Form("Signal Truth Iso ET %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), 400, 0, 50, 4400, -5, 50));
        h_background_truth_isoET.push_back(new TH2D(Form("h_background_truth_isoET_%d", icent), Form("Background Truth Iso ET %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), 400, 0, 50, 4400, -5, 50));
        h_tight_iso_cluster.push_back(new TH1D(Form("h_tight_iso_cluster_%d", icent), Form("Tight Iso Cluster %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_pT_bins, pT_bin_edges));
        h_tight_noniso_cluster.push_back(new TH1D(Form("h_tight_noniso_cluster_%d", icent), Form("Tight Non-Iso Cluster %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_pT_bins, pT_bin_edges));
        h_nontight_iso_cluster.push_back(new TH1D(Form("h_nontight_iso_cluster_%d", icent), Form("Non-Tight Iso Cluster %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_pT_bins, pT_bin_edges));
        h_nontight_noniso_cluster.push_back(new TH1D(Form("h_nontight_noniso_cluster_%d", icent), Form("Non-Tight Non-Iso Cluster %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_pT_bins, pT_bin_edges));
        h_common_cluster.push_back(new TH1D(Form("h_common_cluster_%d", icent), Form("Common Cluster %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_pT_bins, pT_bin_edges));
        h_all_cluster.push_back(new TH1D(Form("h_all_cluster_%d", icent), Form("All Cluster %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_pT_bins, pT_bin_edges));
        h_tight_cluster.push_back(new TH1D(Form("h_tight_cluster_%d", icent), Form("Tight Cluster %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_pT_bins, pT_bin_edges));
        h_pT_truth_response.push_back(new TH1D(Form("h_pT_truth_response_%d", icent), Form("Truth pT %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_pT_bins_truth, pT_bin_edges_truth));
        h_pT_reco_response.push_back(new TH1D(Form("h_pT_reco_response_%d", icent), Form("Reco pT %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_pT_bins, pT_bin_edges));
        h_pT_reco_fake.push_back(new TH1D(Form("h_pT_reco_fake_%d", icent), Form("Reco Fake pT %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_pT_bins, pT_bin_edges));

        std::cout << h_pT_reco_response[icent]->GetNbinsX() << std::endl;
        TH2D *h_response_full = new TH2D(Form("h_response_full_%d", icent), Form("Response Matrix %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_pT_bins, pT_bin_edges, n_pT_bins_truth, pT_bin_edges_truth);
        responses_full.push_back(new RooUnfoldResponse((const TH1 *)h_pT_reco_response[icent], (const TH1 *)h_pT_truth_response[icent], h_response_full, Form("response_matrix_full_%d", icent), "", false));
        h_response_full_list.push_back(h_response_full);
        h_pT_truth_half_response.push_back(new TH1D(Form("h_pT_truth_half_response_%d", icent), Form("Truth pT %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_pT_bins_truth, pT_bin_edges_truth));
        h_pT_reco_half_response.push_back(new TH1D(Form("h_pT_reco_half_response_%d", icent), Form("Reco pT %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_pT_bins, pT_bin_edges));
        h_pT_reco_secondhalf_response.push_back(new TH1D(Form("h_pT_reco_secondhalf_response_%d", icent), Form("Reco pT %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_pT_bins, pT_bin_edges));
        h_pT_truth_secondhalf_response.push_back(new TH1D(Form("h_pT_truth_secondhalf_response_%d", icent), Form("Truth pT %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_pT_bins_truth, pT_bin_edges_truth));
        TH2D *h_response_half = new TH2D(Form("h_response_half_%d", icent), Form("Response Matrix %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_pT_bins, pT_bin_edges, n_pT_bins_truth, pT_bin_edges_truth);
        responses_half.push_back(new RooUnfoldResponse(h_pT_reco_half_response[icent], h_pT_truth_half_response[icent], h_response_half, Form("response_matrix_half_%d", icent), ""));
        h_response_half_list.push_back(h_response_half);
        h_response_xjgamma_global_list.push_back(new TH2D(
            Form("h_response_xjgamma_global_cent%d", icent),
            Form("Flattened response in (p_{T}^{#gamma},x_{J#gamma}) %.0f-%.0f%% cent;truth global bin (p_{T}^{#gamma},x_{J#gamma});reco global bin (p_{T}^{#gamma},x_{J#gamma})",
                 centrality_bins[icent], centrality_bins[icent + 1]),
            n_global_xjgamma_truth_bins, 0.5, n_global_xjgamma_truth_bins + 0.5,
            n_global_xjgamma_reco_bins, 0.5, n_global_xjgamma_reco_bins + 0.5));
        h_xjgamma_truth_response.push_back(new TH2D(
            Form("h_xjgamma_truth_response_cent%d", icent),
            Form("Truth (x_{J#gamma},p_{T}^{#gamma}) %.0f-%.0f%% cent;x_{J#gamma}^{truth};p_{T,truth}^{#gamma} [GeV]",
                 centrality_bins[icent], centrality_bins[icent + 1]),
            n_xj_bins_truth, xj_bin_edges_truth, n_pT_bins_truth, pT_bin_edges_truth));
        h_xjgamma_reco_response.push_back(new TH2D(
            Form("h_xjgamma_reco_response_cent%d", icent),
            Form("Reco (x_{J#gamma},p_{T}^{#gamma}) %.0f-%.0f%% cent;x_{J#gamma}^{reco};p_{T,reco}^{#gamma} [GeV]",
                 centrality_bins[icent], centrality_bins[icent + 1]),
            n_xj_bins, xj_bin_edges, n_pT_bins, pT_bin_edges));
        responses_xjgamma_2d.push_back(new RooUnfoldResponse(
            h_xjgamma_reco_response[icent], h_xjgamma_truth_response[icent],
            Form("response_matrix_xjgamma_2d_cent%d", icent), ""));
        h_ncluster_truth.push_back(new TH2D(Form("h_ncluster_truth_%d", icent), Form("N Cluster From Truth %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), 400, 0, 40, 10, 0, 10));
        h_pT_truth_reco.push_back(new TH2D(Form("h_pT_truth_reco_%d", icent), Form("Truth Reco %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), 400, 5, 40, 150, 0, 1.5));
        h_jet_pT_response.push_back(new TH2D(
            Form("h_jet_pT_response_cent%d", icent),
            Form("Jet p_{T} response;Truth jet p_{T} [GeV];Reco jet p_{T} / Truth jet p_{T} (%.0f-%.0f%% cent)", centrality_bins[icent], centrality_bins[icent + 1]),
            45, 5, 50,   // 1 GeV bins from 5 to 50 GeV
            100, 0, 2)); // response axis
        h_direct_pT_truth_isoET.push_back(new TH2D(Form("h_direct_pT_truth_isoET_%d", icent), Form("Direct Photon pT vs Truth Iso ET %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), 400, 0, 40, 400, 0, 40));
        h_frag_pT_truth_isoET.push_back(new TH2D(Form("h_frag_pT_truth_isoET_%d", icent), Form("Fragmentation Photon pT vs Truth Iso ET %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), 400, 0, 40, 400, 0, 40));

        // xjgamma (variable binning from xj_bin_edges on x-axis)
        h_tight_iso_xjgamma_signal.push_back(new TH2D(Form("h_tight_iso_xjgamma_signal_cent%d", icent), Form("Tight Iso XJGamma Signal %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_xj_bins, xj_bin_edges, n_pT_bins, pT_bin_edges));
        h_tight_noniso_xjgamma_signal.push_back(new TH2D(Form("h_tight_noniso_xjgamma_signal_cent%d", icent), Form("Tight Non-Iso XJGamma Signal %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_xj_bins, xj_bin_edges, n_pT_bins, pT_bin_edges));
        h_nontight_iso_xjgamma_signal.push_back(new TH2D(Form("h_nontight_iso_xjgamma_signal_cent%d", icent), Form("Non-Tight Iso XJGamma Signal %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_xj_bins, xj_bin_edges, n_pT_bins, pT_bin_edges));
        h_nontight_noniso_xjgamma_signal.push_back(new TH2D(Form("h_nontight_noniso_xjgamma_signal_cent%d", icent), Form("Non-Tight Non-Iso XJGamma Signal %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_xj_bins, xj_bin_edges, n_pT_bins, pT_bin_edges));
        h_all_xjgamma_signal.push_back(new TH2D(Form("h_all_xjgamma_signal_cent%d", icent), Form("All XJGamma Signal %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_xj_bins, xj_bin_edges, n_pT_bins, pT_bin_edges));
        h_tight_xjgamma_signal.push_back(new TH2D(Form("h_tight_xjgamma_signal_cent%d", icent), Form("Tight XJGamma Signal %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_xj_bins, xj_bin_edges, n_pT_bins, pT_bin_edges));
        h_tight_iso_truthmatchreco_xjgamma_signal.push_back(new TH2D(Form("h_tight_iso_truthmatchreco_xjgamma_signal_cent%d", icent), Form("Tight Iso TruthMatchedRecoJet XJGamma Signal %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_xj_bins, xj_bin_edges, n_pT_bins, pT_bin_edges));
        h_tight_noniso_truthmatchreco_xjgamma_signal.push_back(new TH2D(Form("h_tight_noniso_truthmatchreco_xjgamma_signal_cent%d", icent), Form("Tight Non-Iso TruthMatchedRecoJet XJGamma Signal %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_xj_bins, xj_bin_edges, n_pT_bins, pT_bin_edges));
        h_nontight_iso_truthmatchreco_xjgamma_signal.push_back(new TH2D(Form("h_nontight_iso_truthmatchreco_xjgamma_signal_cent%d", icent), Form("Non-Tight Iso TruthMatchedRecoJet XJGamma Signal %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_xj_bins, xj_bin_edges, n_pT_bins, pT_bin_edges));
        h_nontight_noniso_truthmatchreco_xjgamma_signal.push_back(new TH2D(Form("h_nontight_noniso_truthmatchreco_xjgamma_signal_cent%d", icent), Form("Non-Tight Non-Iso TruthMatchedRecoJet XJGamma Signal %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_xj_bins, xj_bin_edges, n_pT_bins, pT_bin_edges));
        h_all_truthmatchreco_xjgamma_signal.push_back(new TH2D(Form("h_all_truthmatchreco_xjgamma_signal_cent%d", icent), Form("All TruthMatchedRecoJet XJGamma Signal %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_xj_bins, xj_bin_edges, n_pT_bins, pT_bin_edges));
        h_tight_truthmatchreco_xjgamma_signal.push_back(new TH2D(Form("h_tight_truthmatchreco_xjgamma_signal_cent%d", icent), Form("Tight TruthMatchedRecoJet XJGamma Signal %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_xj_bins, xj_bin_edges, n_pT_bins, pT_bin_edges));
        h_tight_iso_truthjet_xjgamma_signal.push_back(new TH2D(Form("h_tight_iso_truthjet_xjgamma_signal_cent%d", icent), Form("Tight Iso TruthJet XJGamma Signal %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_xj_bins_truth, xj_bin_edges_truth, n_pT_bins, pT_bin_edges));
        h_tight_noniso_truthjet_xjgamma_signal.push_back(new TH2D(Form("h_tight_noniso_truthjet_xjgamma_signal_cent%d", icent), Form("Tight Non-Iso TruthJet XJGamma Signal %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_xj_bins_truth, xj_bin_edges_truth, n_pT_bins, pT_bin_edges));
        h_nontight_iso_truthjet_xjgamma_signal.push_back(new TH2D(Form("h_nontight_iso_truthjet_xjgamma_signal_cent%d", icent), Form("Non-Tight Iso TruthJet XJGamma Signal %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_xj_bins_truth, xj_bin_edges_truth, n_pT_bins, pT_bin_edges));
        h_nontight_noniso_truthjet_xjgamma_signal.push_back(new TH2D(Form("h_nontight_noniso_truthjet_xjgamma_signal_cent%d", icent), Form("Non-Tight Non-Iso TruthJet XJGamma Signal %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_xj_bins_truth, xj_bin_edges_truth, n_pT_bins, pT_bin_edges));
        h_all_truthjet_xjgamma_signal.push_back(new TH2D(Form("h_all_truthjet_xjgamma_signal_cent%d", icent), Form("All TruthJet XJGamma Signal %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_xj_bins_truth, xj_bin_edges_truth, n_pT_bins, pT_bin_edges));
        h_tight_truthjet_xjgamma_signal.push_back(new TH2D(Form("h_tight_truthjet_xjgamma_signal_cent%d", icent), Form("Tight TruthJet XJGamma Signal %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_xj_bins_truth, xj_bin_edges_truth, n_pT_bins, pT_bin_edges));
        h_tight_iso_xjgamma_background.push_back(new TH2D(Form("h_tight_iso_xjgamma_background_cent%d", icent), Form("Tight Iso XJGamma Background %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_xj_bins, xj_bin_edges, n_pT_bins, pT_bin_edges));
        h_tight_noniso_xjgamma_background.push_back(new TH2D(Form("h_tight_noniso_xjgamma_background_cent%d", icent), Form("Tight Non-Iso XJGamma Background %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_xj_bins, xj_bin_edges, n_pT_bins, pT_bin_edges));
        h_nontight_iso_xjgamma_background.push_back(new TH2D(Form("h_nontight_iso_xjgamma_background_cent%d", icent), Form("Non-Tight Iso XJGamma Background %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_xj_bins, xj_bin_edges, n_pT_bins, pT_bin_edges));
        h_nontight_noniso_xjgamma_background.push_back(new TH2D(Form("h_nontight_noniso_xjgamma_background_cent%d", icent), Form("Non-Tight Non-Iso XJGamma Background %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_xj_bins, xj_bin_edges, n_pT_bins, pT_bin_edges));
        h_tight_iso_xjgamma.push_back(new TH2D(Form("h_tight_iso_xjgamma_cent%d", icent), Form("Tight Iso XJGamma %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_xj_bins, xj_bin_edges, n_pT_bins, pT_bin_edges));
        h_tight_noniso_xjgamma.push_back(new TH2D(Form("h_tight_noniso_xjgamma_cent%d", icent), Form("Tight Non-Iso XJGamma %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_xj_bins, xj_bin_edges, n_pT_bins, pT_bin_edges));
        h_nontight_iso_xjgamma.push_back(new TH2D(Form("h_nontight_iso_xjgamma_cent%d", icent), Form("Non-Tight Iso XJGamma %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_xj_bins, xj_bin_edges, n_pT_bins, pT_bin_edges));
        h_nontight_noniso_xjgamma.push_back(new TH2D(Form("h_nontight_noniso_xjgamma_cent%d", icent), Form("Non-Tight Non-Iso XJGamma %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_xj_bins, xj_bin_edges, n_pT_bins, pT_bin_edges));
        h_common_xjgamma.push_back(new TH2D(Form("h_common_xjgamma_cent%d", icent), Form("Common XJGamma %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_xj_bins, xj_bin_edges, n_pT_bins, pT_bin_edges));
        h_all_xjgamma.push_back(new TH2D(Form("h_all_xjgamma_cent%d", icent), Form("All XJGamma %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_xj_bins, xj_bin_edges, n_pT_bins, pT_bin_edges));
        h_tight_xjgamma.push_back(new TH2D(Form("h_tight_xjgamma_cent%d", icent), Form("Tight XJGamma %.0f-%.0f%% cent", centrality_bins[icent], centrality_bins[icent + 1]), n_xj_bins, xj_bin_edges, n_pT_bins, pT_bin_edges));
        for (int ipt = 0; ipt < n_pT_bins; ipt++)
        {
            h_tight_cluster_pT[icent].push_back(new TH1D(Form("h_tight_isoET_cent%d_pt%d", icent, ipt), Form("Tight Iso ET %.0f-%.0f%% cent %.1f < pT < %.1f", centrality_bins[icent], centrality_bins[icent + 1], pT_bin_edges[ipt], pT_bin_edges[ipt + 1]), 1000, -50, 50));
            h_nontight_cluster_pT[icent].push_back(new TH1D(Form("h_nontight_isoET_cent%d_pt%d", icent, ipt), Form("Non-Tight Iso ET %.0f-%.0f%% cent %.1f < pT < %.1f", centrality_bins[icent], centrality_bins[icent + 1], pT_bin_edges[ipt], pT_bin_edges[ipt + 1]), 1000, -50, 50));
            h_iso_truth_reco[icent].push_back(new TH2D(Form("h_iso_truth_reco_cent%d_pt%d", icent, ipt), Form("Iso Truth Reco %.0f-%.0f%% cent %.1f < pT < %.1f", centrality_bins[icent], centrality_bins[icent + 1], pT_bin_edges[ipt], pT_bin_edges[ipt + 1]), 300, 0, 30, 1000, -50, 50));
            h_background_iso_truth_reco[icent].push_back(new TH2D(Form("h_background_iso_truth_reco_cent%d_pt%d", icent, ipt), Form("Background Truth Iso ET %.0f-%.0f%% cent %.1f < pT < %.1f", centrality_bins[icent], centrality_bins[icent + 1], pT_bin_edges[ipt], pT_bin_edges[ipt + 1]), 300, 0, 30, 1000, -50, 50));
            h_response_isoET[icent].push_back(new TH2D(Form("h_response_isoET_cent%d_pt%d", icent, ipt), Form("Response Iso ET %.0f-%.0f%% cent %.1f < pT < %.1f", centrality_bins[icent], centrality_bins[icent + 1], pT_bin_edges[ipt], pT_bin_edges[ipt + 1]), 150, 0, 1.5, 1000, -50, 50));
            h_dphi_clusterJets_tight[icent].push_back(new TH1D(
                Form("h_dphi_clusterJets_tight_cent%d_pt%d", icent, ipt),
                Form("Tight |#Delta#phi|(cluster, jet) (xj#gamma reco jets: |#eta|<cut, |#Delta#phi|>b2b, p_{T}^{cal}>min) %.0f-%.0f%% cent, %.0f < E_{T}^{clus} < %.0f GeV;|#Delta#phi|;Pairs",
                     centrality_bins[icent], centrality_bins[icent + 1], pT_bin_edges[ipt], pT_bin_edges[ipt + 1]),
                64, 0, M_PI));
            h_dphi_clusterJets_tight[icent][ipt]->Sumw2();
        }
    }

    TRandom3 *rand = new TRandom3(0);
    std::set<int> skiprunnumbers = {0};
    int nentries = chain.GetEntries();
    int ientry = 0;








    /////////////////////////////////////
    // Event loop
    /////////////////////////////////////
    while (reader.Next())
    {
        // Event-by-event MC vertex reweighting.
        // We directly update `weight` so all existing Fill(..., weight) calls use the per-event weight.
        weight = cross_weight;
        vertex_weight = 1.0;
        cent_weight = 1.0;
        if (issim)
        {
            if (vertex_reweight_on)
            {
                int bin = h_vertex_reweight->FindBin(*vertexz);
                if (bin < 1) bin = 1;
                if (bin > h_vertex_reweight->GetNbinsX()) bin = h_vertex_reweight->GetNbinsX();
                vertex_weight = h_vertex_reweight->GetBinContent(bin);
            }

            if (!std::isfinite(vertex_weight) || vertex_weight <= 0.0)
            {
                std::cout << "Warning: vertex weight is nan or inf" << std::endl;
                std::cout << "vertexz: " << *vertexz << std::endl;
                vertex_weight = 1.0;
            }
        }

        if (ientry < 0)
        {
            ientry++;
            continue;
        }

        if (ientry % 10000 == 0)  std::cout << "Processing entry " << ientry << " / " << nentries << std::endl;


        if (!issim)
        {
            // Accept event if ANY trigger in trigger_used fired.
            bool any_trigger_fired = false;
            const auto ntrig = scaledtrigger.GetSize();
            for (int itrig : trigger_used)
            {
                if (itrig < 0 || (unsigned int)itrig >= ntrig)
                {
                    continue;
                }
                if (scaledtrigger[itrig] != 0)
                {
                    any_trigger_fired = true;
                    break;
                }
            }

            if (!any_trigger_fired)
            {
                ientry++;
                continue;
            }
        }

        // Centrality bin for this event: pp uses single bin (0); Au+Au uses branch "cent"
        int centbin = -1;
        float cent_percent = 50.0f;
        if (is_pp)
        {
            centbin = 0;
        }
        else if (cent_reader)
        {
            cent_percent = (*(*cent_reader)) * 100.0f;
            for (int icent = 0; icent < n_cent_bins; icent++)
            {
                if (cent_percent >= centrality_bins[icent] && cent_percent < centrality_bins[icent + 1])
                {
                    centbin = icent;
                    break;
                }
            }
        }
        // Skip events outside analysis centrality bins before applying centrality reweighting.
        if (centbin < 0 || centbin >= n_cent_bins)
        {
            ientry++;
            continue;
        }
        if (issim && !is_pp && cent_reweight_on)
        {
            if (!h_cent_reweight)
            {
                std::cerr << "[CentralityReweight] ERROR: centrality reweighting is enabled but histogram is not loaded."
                          << std::endl;
                return;
            }
            int bin = h_cent_reweight->FindBin(cent_percent);
            if (bin < 1) bin = 1;
            if (bin > h_cent_reweight->GetNbinsX()) bin = h_cent_reweight->GetNbinsX();
            cent_weight = h_cent_reweight->GetBinContent(bin);
            if (!std::isfinite(cent_weight) || cent_weight <= 0.0)
            {
                std::cout << "Warning: centrality weight is nan/inf/non-positive, reset to 1.0" << std::endl;
                std::cout << "cent_percent: " << cent_percent << std::endl;
                cent_weight = 1.0;
            }
        }
        if (issim)
        {
            weight *= vertex_weight * cent_weight;
        }
        std::map<int, int> particle_trkidmap;
        // map for photon reco eff
        std::map<int, bool> photon_converts;
        // map for photon reco eff
        std::map<int, bool> photon_reco;
        // map for iso eff
        std::map<int, bool> photon_iso;
        std::map<int, float> photon_iso_ET;
        std::map<int, float> all_iso_ET;
        // map for id eff
        std::map<int, bool> photon_id;
        // map for n cluster per truth photon
        std::map<int, int> photon_ncluster;
        // event-level mask: truth jets matched to truth photons (DeltaR < 0.2)
        std::vector<bool> truthjet_is_photon_matched;
        // per reco jet: true if matched to some truth jet with DeltaR < 0.2 (filled after cent/vertex cuts)
        std::vector<bool> reco_jet_truth_matched;
        if (issim)
        {
            truthjet_is_photon_matched.assign(*njet_truth, false);
            float maxphotonpT = 0;
            int maxphotonclass = 0;
            for (int iparticle = 0; iparticle < *nparticles; iparticle++)
            {
                particle_trkidmap[particle_trkid[iparticle]] = iparticle;

                // Mark truth jets that geometrically match truth photons; these are excluded
                // from the jet-response histograms.
                if (particle_pid[iparticle] == 22)
                {
                    for (int ijet = 0; ijet < *njet_truth; ++ijet)
                    {
                        float dphi = particle_Phi[iparticle] - jet_truth_Phi[ijet];
                        while (dphi > M_PI) dphi -= 2 * M_PI;
                        while (dphi < -M_PI) dphi += 2 * M_PI;
                        const float deta = particle_Eta[iparticle] - jet_truth_Eta[ijet];
                        const float dr = std::sqrt(deta * deta + dphi * dphi);
                        if (dr < 0.2f)
                            truthjet_is_photon_matched[ijet] = true;
                    }
                }
                // if (!isbackground)
                float truthisoET = 0;
                if (conesize == 4)
                {
                    truthisoET = particle_truth_iso_04[iparticle];
                }
                else if (conesize == 3)
                {
                    truthisoET = particle_truth_iso_03[iparticle];
                }
                else if (conesize == 2)
                {
                    truthisoET = particle_truth_iso_02[iparticle];
                }
                else
                {
                    std::cout << "Error: conesize not supported" << std::endl;
                    continue;
                }

                if (particle_pid[iparticle] == 22)
                {
                    if (particle_Pt[iparticle] > maxphotonpT)
                    {
                        maxphotonpT = particle_Pt[iparticle];
                        maxphotonclass = particle_photonclass[iparticle];
                    }


                    if (particle_photonclass[iparticle] < 3) // direct or fragmentation
                    {
                        if (truthisoET < truthisocut)
                        {
                            // we check for conversion here
                            photon_converts[iparticle] = (particle_converted[iparticle] > 0);
                            // initialize reco to false
                            photon_reco[iparticle] = false;
                            // initialize id to false
                            photon_id[iparticle] = false;
                            // initialize iso to false
                            photon_iso[iparticle] = false;

                            photon_ncluster[iparticle] = 0;

                            // fill the truth histogram
                            float photonpT = particle_Pt[iparticle];
                            float photon_eta = particle_Eta[iparticle];

                            if (photon_eta <= eta_acceptance_min || photon_eta >= eta_acceptance_max)
                            {
                                continue;
                            }
                        }
                        photon_iso_ET[iparticle] = truthisoET;
                    }
                }
                all_iso_ET[iparticle] = truthisoET;
            }
            if (!isbackground)
            {
                if ((maxphotonpT > max_photon_upper) || (maxphotonpT < max_photon_lower))
                {
                    ientry++;
                    continue;
                }
            }
            // if (abs(particle_Eta[iparticle]) < 0.7)
            {
                h_max_photon_pT->Fill(maxphotonpT, weight);

                if (maxphotonclass == 1)
                {
                    h_max_direct_pT->Fill(maxphotonpT, weight);
                }
                else if (maxphotonclass == 2)
                {
                    h_max_frag_pT->Fill(maxphotonpT, weight);
                }
                else if (maxphotonclass == 3)
                {
                    h_max_decay_pT->Fill(maxphotonpT, weight);
                }
            }

            // loop over truth_jets for the background sample
            if (isbackground)
            {
                float maxjetpT = 0;
                for (int ijet = 0; ijet < *njet_truth; ijet++)
                {
                    if (jet_truth_Pt[ijet] > maxjetpT)
                    {
                        maxjetpT = jet_truth_Pt[ijet];
                    }
                }
                if (maxjetpT == 0) continue;

                if ((maxjetpT > max_jet_upper) )
                {
                    ientry++;
                    continue;
                }
                if (maxjetpT < max_jet_lower){
                    ientry++;
                    continue;
                }
                //std::cout<<"accepted maxjetpT: "<<maxjetpT<<std::endl;
                h_max_truth_jet_pT->Fill(maxjetpT, weight);
            }

            if (!isbackground)
            {
                for (int iparticle = 0; iparticle < *nparticles; iparticle++)
                {
                    particle_trkidmap[particle_trkid[iparticle]] = iparticle;

                    if (particle_pid[iparticle] == 22)
                    {
                        if (particle_Pt[iparticle] > maxphotonpT)
                        {
                            maxphotonpT = particle_Pt[iparticle];
                            maxphotonclass = particle_photonclass[iparticle];
                        }

                        float truthisoET = 0;
                        if (conesize == 4)
                        {
                            truthisoET = particle_truth_iso_04[iparticle];
                        }
                        else if (conesize == 3)
                        {
                            truthisoET = particle_truth_iso_03[iparticle];
                        }
                        else if (conesize == 2)
                        {
                            truthisoET = particle_truth_iso_02[iparticle];
                        }
                        else
                        {
                            std::cout << "Error: conesize not supported" << std::endl;
                            continue;
                        }

                        if (particle_Eta[iparticle] > eta_acceptance_min && particle_Eta[iparticle] < eta_acceptance_max && centbin >= 0 && centbin < n_cent_bins)
                        {
                            if (particle_photonclass[iparticle] == 1)
                            {
                                h_direct_pT_truth_isoET[centbin]->Fill(particle_Pt[iparticle], truthisoET, weight);
                                h_direct_pT->Fill(particle_Pt[iparticle], weight);
                            }
                            else if (particle_photonclass[iparticle] == 2)
                            {
                                h_frag_pT_truth_isoET[centbin]->Fill(particle_Pt[iparticle], truthisoET, weight);
                                h_frag_pT->Fill(particle_Pt[iparticle], weight);
                            }
                            else if (particle_photonclass[iparticle] == 3)
                            {
                                h_decay_photon_pT->Fill(particle_Pt[iparticle], weight);
                            }

                            h_photon_pT->Fill(particle_Pt[iparticle], weight);
                        }

                        if (particle_photonclass[iparticle] < 3) // direct or fragmentation
                        {
                            if (truthisoET < truthisocut)
                            {
                                float photonpT = particle_Pt[iparticle];
                                float photon_eta = particle_Eta[iparticle];

                                if (photon_eta <= eta_acceptance_min || photon_eta >= eta_acceptance_max)
                                {
                                    continue;
                                }

                                if (centbin >= 0 && centbin < n_cent_bins)
                                {
                                    h_truth_pT[centbin]->Fill(photonpT, weight);
                                }

                                // check if truth vertex is within the vertex cut
                                if (abs(*vertexz_truth) < vertexcut)
                                {
                                    if (centbin >= 0 && centbin < n_cent_bins)
                                        h_truth_pT_vertexcut[centbin]->Fill(photonpT, weight);
                                }
                            }
                        }
                    }
                }
            }
        }

        if (abs(*vertexz) > vertexcut)
        {
            ientry++;
            continue;
        }

        if (centbin < 0 || centbin >= n_cent_bins)
        {
            ientry++;
            continue;
        }

        constexpr float match_dr_reco_truth = 0.2f;
        if (issim)
        {
            reco_jet_truth_matched.assign(*njet, false);
            for (int ijet_reco = 0; ijet_reco < *njet; ++ijet_reco)
            {
                for (int ijet_truth = 0; ijet_truth < *njet_truth; ++ijet_truth)
                {
                    if (truthjet_is_photon_matched[ijet_truth])
                        continue;
                    float dphi = jet_Phi[ijet_reco] - jet_truth_Phi[ijet_truth];
                    while (dphi > M_PI)
                        dphi -= 2 * M_PI;
                    while (dphi < -M_PI)
                        dphi += 2 * M_PI;
                    const float deta = jet_Eta[ijet_reco] - jet_truth_Eta[ijet_truth];
                    const float dr = std::sqrt(deta * deta + dphi * dphi);
                    if (dr < match_dr_reco_truth)
                    {
                        reco_jet_truth_matched[ijet_reco] = true;
                        break;
                    }
                }
            }
        }
        else
        {
            reco_jet_truth_matched.assign(*njet, false);
        }

        std::vector<float> jetphi;

        for (int ijet = 0; ijet < *njet; ijet++)
        {
            if (jet_Pt[ijet] < 10)
                continue;
            jetphi.push_back(jet_Phi[ijet]);
        }

        // Jet pT response: match each truth jet to nearest reco jet within DeltaR < 0.2
        // and fill (truth jet pT, reco/truth pT) per centrality bin.
        if (issim)
        {
            for (int ijet_truth = 0; ijet_truth < *njet_truth; ++ijet_truth)
            {
                if (truthjet_is_photon_matched[ijet_truth])
                    continue;
                const float truth_pt = jet_truth_Pt[ijet_truth];
                if (truth_pt <= 0.f)
                    continue;
                if (std::abs(jet_truth_Eta[ijet_truth]) >= jet_eta)
                    continue;

                float best_dr = match_dr_reco_truth;
                int best_reco = -1;
                for (int ijet_reco = 0; ijet_reco < *njet; ++ijet_reco)
                {
                    if (std::abs(jet_Eta[ijet_reco]) >= jet_eta)
                        continue;

                    float dphi = jet_Phi[ijet_reco] - jet_truth_Phi[ijet_truth];
                    while (dphi > M_PI)
                        dphi -= 2 * M_PI;
                    while (dphi < -M_PI)
                        dphi += 2 * M_PI;
                    const float deta = jet_Eta[ijet_reco] - jet_truth_Eta[ijet_truth];
                    const float dr = std::sqrt(deta * deta + dphi * dphi);
                    if (dr < best_dr)
                    {
                        best_dr = dr;
                        best_reco = ijet_reco;
                    }
                }

                if (best_reco < 0)
                    continue;

                const float response = jet_Pt[best_reco] / truth_pt;
                h_jet_pT_response[centbin]->Fill(truth_pt, response, weight);
            }
        }

        h_vertexz->Fill(*vertexz, weight);
        float leading_common_cluster_ET = 0;
        int leading_cluster_ET_index = -1;
        float leading_cluster_ET = 0;
        for (int icluster = 0; icluster < *ncluster; icluster++)
        {
            if (issim)
            {
                cluster_Et[icluster] = cluster_Et[icluster] * clusterescale;
                if (clustereres > 0)
                {
                    cluster_Et[icluster] = cluster_Et[icluster] * rand->Gaus(1, clustereres);
                }
            }
            // need ET > 10 GeV
            if (cluster_Et[icluster] < reco_min_ET)
                continue;
            if (nosat)
            {
                if (cluster_nsaturated[icluster] > 0)
                    continue;
            }
            if (cluster_Et[icluster] > leading_common_cluster_ET)
            {
                leading_cluster_ET = cluster_Et[icluster];
                leading_cluster_ET_index = icluster;
            }
        }

        for (int icluster = 0; icluster < *ncluster; icluster++)
        {
            bool is_leading_cluster = (icluster == leading_cluster_ET_index);
            // need ET > 10 GeV
            if (cluster_Et[icluster] < reco_min_ET)
                continue;
            
            if (isbackground)
            {
                if (cluster_Et[icluster] > cluster_ET_upper)
                {
                    continue;
                }
            }
            if (nosat)
            {
                if (cluster_nsaturated[icluster] > 0)
                    continue;
            }
            float rhad22 = (cluster_ihcal_et22[icluster] + cluster_ohcal_et22[icluster]) / (cluster_Et[icluster] + (cluster_ihcal_et22[icluster] + cluster_ohcal_et22[icluster]));
            float rhad33 = (cluster_ihcal_et33[icluster] + cluster_ohcal_et33[icluster]) / (cluster_Et[icluster] + (cluster_ihcal_et22[icluster] + cluster_ohcal_et22[icluster]));

            float reta77 = cluster_e37[icluster] / cluster_e77[icluster];
            float rphi77 = cluster_e73[icluster] / cluster_e77[icluster];

            float reta55 = cluster_e35[icluster] / cluster_e55[icluster];
            float rphi55 = cluster_e53[icluster] / cluster_e55[icluster];

            float reta = cluster_e33[icluster] / cluster_e73[icluster];

            float rphi = cluster_e33[icluster] / cluster_e37[icluster];

            float re11_E = cluster_e11[icluster] / cluster_E[icluster];

            float e11_over_e33 = cluster_e11[icluster] / cluster_e33[icluster];

            float e32_over_e35 = cluster_e32[icluster] / cluster_e35[icluster];

            float wr_cogx = cluster_wphi_cogx[icluster] / cluster_weta_cogx[icluster];

            // reco cut (cone-based isolation only; topo cluster iso not used)
            float recoisoET = -999;
            if (conesize == 4)
            {
                recoisoET = cluster_iso_04[icluster];
            }
            else if (conesize == 3)
            {
                recoisoET = cluster_iso_03[icluster];
            }
            else if (conesize == 2)
            {
                recoisoET = cluster_iso_02[icluster];
            }
            else
            {
                std::cout << "Error: conesize not supported" << std::endl;
                continue;
            }

            if (!use_topo_iso && iso_threshold)
            {
                if (iso_hcalonly)
                {
                    recoisoET = cluster_iso_03_70_hcalin[icluster] + cluster_iso_03_70_hcalout[icluster];
                }
                else
                {
                    recoisoET = cluster_iso_03_70_emcal[icluster] + cluster_iso_03_70_hcalin[icluster] + cluster_iso_03_70_hcalout[icluster];
                    //recoisoET = cluster_iso_03_120_emcal[icluster] + cluster_iso_03_120_hcalin[icluster] + cluster_iso_03_120_hcalout[icluster];
                    //recoisoET = cluster_iso_03_60_emcal[icluster] + cluster_iso_03_60_hcalin[icluster] + cluster_iso_03_60_hcalout[icluster];
                    if (iso_emcalinnerr > 0.04 && iso_emcalinnerr < 0.06)
                    {
                        recoisoET -= cluster_iso_005_70_emcal[icluster];
                    }
                    else if (iso_emcalinnerr > 0.09 && iso_emcalinnerr < 0.11)
                    {
                        recoisoET -= cluster_iso_01_70_emcal[icluster];
                    }
                    else if (iso_emcalinnerr > 0.19 && iso_emcalinnerr < 0.21)
                    {
                        recoisoET -= cluster_iso_02_70_emcal[icluster];
                    }
                }
            }
            // fudge the MC isoET
            if (issim)
            {
                recoisoET = recoisoET * mc_iso_scale;
                recoisoET += mc_iso_shift;
            }

            float cluster_eta = cluster_Eta[icluster];
            if (cluster_eta <= eta_acceptance_min || cluster_eta >= eta_acceptance_max)
            {
                continue;
            }

            float pTbin = -1;
            float clusterET = cluster_Et[icluster];
            for (int ipt = 0; ipt < n_pT_bins; ipt++)
            {
                if (clusterET > pT_bins[ipt] && clusterET < pT_bins[ipt + 1])
                {
                    pTbin = ipt;
                    break;
                }
            }

            bool otherside_jet = false;
            vector<float> b2bjet_pT = {};
            vector<float> b2bjet_truthmatched_pT = {};
            vector<float> b2btruthjet_pT = {};

            for (int ijet = 0; ijet < *njet; ijet++)
            {
                float dphi = cluster_Phi[icluster] - jet_Phi[ijet];

                while (dphi > M_PI)
                    dphi = dphi - 2 * M_PI;
                while (dphi < -M_PI)
                    dphi = dphi + 2 * M_PI;

                if (abs(dphi) > (M_PI / 2))
                {
                    otherside_jet = true;
                    // break;
                }
                float calibrated_jet_pT = 1/0.65*jet_Pt[ijet];  // JES calibration removed – jets already calibrated
                if (abs(jet_Eta[ijet]) < jet_eta)
                {
                    if (abs(dphi) > b2bjet_dphi)
                    {
                        if (calibrated_jet_pT > b2bjet_pT_min)
                        {
                            b2bjet_pT.push_back(calibrated_jet_pT);
                            if (issim && ijet < static_cast<int>(reco_jet_truth_matched.size()) && reco_jet_truth_matched[ijet])
                            {
                                b2bjet_truthmatched_pT.push_back(calibrated_jet_pT);
                            }
                        }
                    }
                }
            }
            std::sort(b2bjet_pT.begin(), b2bjet_pT.end());
            std::sort(b2bjet_truthmatched_pT.begin(), b2bjet_truthmatched_pT.end());
            const float max_b2bjet_pT = b2bjet_pT.empty() ? -1.f : b2bjet_pT.back();

            bool passes_common_b2bjet = (!common_b2bjet_cut) || (max_b2bjet_pT >= common_b2bjet_pt_min);
            // One xj entry per back-to-back jet above threshold; only for leading cluster in the event
            const bool fill_xjgamma = is_leading_cluster && leading_cluster_ET > 0.f;
            
            bool common_pass = false;
            bool tight = false;
            bool nontight = false;
            bool iso = false;
            bool noniso = false;

            float cent_shift = reco_iso_max_cent_c0 + reco_iso_max_cent_c1 * cent_percent + reco_iso_max_cent_c2 * cent_percent * cent_percent;
            recoisoET -= cent_shift;
            float recoiso_max = recoiso_max_b + recoiso_max_s * clusterET;// + cent_shift;   
            float recononiso_min = recoiso_max + recononiso_min_shift;

            if (recoisoET > recoiso_min && recoisoET < recoiso_max)
            {
                iso = true;
            }
            if (recoisoET > recononiso_min && recoisoET < recononiso_max)
            {
                noniso = true;
            }

            // common cuts
            bool passes_common_shape =
                cluster_prob[icluster] > common_prob_min &&
                cluster_prob[icluster] < common_prob_max &&
                e11_over_e33 > common_e11_over_e33_min &&
                e11_over_e33 < common_e11_over_e33_max &&
                //(!(wr_cogx < common_wr_cogx_bound && cluster_weta_cogx[icluster] > common_cluster_weta_cogx_bound))
                (cluster_weta_cogx[icluster] < common_cluster_weta_cogx_bound) &&
                (!common_npb_cut_on || cluster_npb_score[icluster] > common_npb_score_cut);

            if (passes_common_shape && passes_common_b2bjet)
            {
                common_pass = true;
            }

            if (common_pass)
            {
                {
                    tight_weta_cogx_max = tight_weta_cogx_max_b + tight_weta_cogx_max_s * clusterET;
                    tight_wphi_cogx_max = tight_wphi_cogx_max_b + tight_wphi_cogx_max_s * clusterET;
                    tight_et1_min = tight_et1_min_b + tight_et1_min_s * clusterET;
                }
                // need to update to a function to find tight and non tight
                bool is_cluster_weta_cogx_tight =
                    (cluster_weta_cogx[icluster] > tight_weta_cogx_min) &&
                    (cluster_weta_cogx[icluster] < tight_weta_cogx_max);

                bool is_cluster_wphi_cogx_tight =
                    (cluster_wphi_cogx[icluster] > tight_wphi_cogx_min) &&
                    (cluster_wphi_cogx[icluster] < tight_wphi_cogx_max);

                bool is_cluster_et1_tight =
                    (cluster_et1[icluster] > tight_et1_min) &&
                    (cluster_et1[icluster] < tight_et1_max);

                bool is_cluster_et2_tight =
                    (cluster_et2[icluster] > tight_et2_min) &&
                    (cluster_et2[icluster] < tight_et2_max);

                bool is_cluster_et3_tight =
                    (cluster_et3[icluster] > tight_et3_min) &&
                    (cluster_et3[icluster] < tight_et3_max);

                bool is_e11_over_e33_tight =
                    (e11_over_e33 > tight_e11_over_e33_min) &&
                    (e11_over_e33 < tight_e11_over_e33_max);

                bool is_e32_over_e35_tight =
                    (e32_over_e35 > tight_e32_over_e35_min) &&
                    (e32_over_e35 < tight_e32_over_e35_max);

                bool is_cluster_et4_tight =
                    (cluster_et4[icluster] > tight_et4_min) &&
                    (cluster_et4[icluster] < tight_et4_max);

                bool is_cluster_prob_tight =
                    (cluster_prob[icluster] > tight_prob_min) &&
                    (cluster_prob[icluster] < tight_prob_max);
                
                // Select BDT model based on cluster ET (ET-binned or fallback)
                std::string selected_bdt_model = bdt_model_name;
                if (use_et_binned_bdt) {
                    float et = cluster_Et[icluster];
                    for (int ib = 0; ib < (int)bdt_et_bin_models.size(); ++ib) {
                        if (et >= bdt_et_bin_edges[ib] && et < bdt_et_bin_edges[ib + 1]) {
                            selected_bdt_model = bdt_et_bin_models[ib];
                            break;
                        }
                    }
                }
                float bdt_score = (*bdt_arrays[selected_bdt_model])[icluster];

                //std::cout<<"bdt_score: "<<bdt_score<<std::endl;
                float tight_bdt_min_et = tight_bdt_min_slope * cluster_Et[icluster] + tight_bdt_min_intercept;
                bool is_bdt_tight =
                    (bdt_score > tight_bdt_min_et) &&
                    (bdt_score < tight_bdt_max);

                // Combined condition
                if (is_cluster_weta_cogx_tight &&
                    is_cluster_wphi_cogx_tight &&
                    is_cluster_et1_tight &&
                    is_cluster_et2_tight &&
                    is_cluster_et3_tight &&
                    is_e11_over_e33_tight &&
                    is_e32_over_e35_tight &&
                    is_cluster_et4_tight &&
                    is_cluster_prob_tight &&
                    is_bdt_tight)
                {
                    tight = true;

                    // |DeltaPhi|(cluster, jet) for tight clusters — reco jets as xjγ (b2bjet_pT); exclude ΔR(cluster,jet) ≤ 0.2
                    int pTbin_dphi = -1;
                    for (int ipt = 0; ipt < n_pT_bins; ipt++)
                    {
                        if (clusterET >= pT_bins[ipt] && clusterET < pT_bins[ipt + 1])
                        {
                            pTbin_dphi = ipt;
                            break;
                        }
                    }

                    h_totalEMCal_energy_tight_weight->Fill(*totalEMCal_energy, weight);
                    h_totalIHCal_energy_tight_weight->Fill(*totalIHCal_energy, weight);
                    h_totalOHCal_energy_tight_weight->Fill(*totalOHCal_energy, weight);
                    h_centrality_tight_weight->Fill(cent_percent, weight);

                    const float dphi_cluster_jet_dR_min = 0.2f;
                    if (pTbin_dphi >= 0)
                    {
                        for (int ijet = 0; ijet < *njet; ijet++)
                        {
                            float jdphi = cluster_Phi[icluster] - jet_Phi[ijet];
                            while (jdphi > M_PI)
                                jdphi = jdphi - 2 * M_PI;
                            while (jdphi < -M_PI)
                                jdphi = jdphi + 2 * M_PI;

                            const float deta = cluster_Eta[icluster] - jet_Eta[ijet];
                            const float dR = std::sqrt(deta * deta + jdphi * jdphi);
                            if (dR <= dphi_cluster_jet_dR_min)
                                continue;

                            const float calibrated_jet_pT = 1.f / 0.65f * jet_Pt[ijet];
                            if (!(std::abs(jet_Eta[ijet]) < jet_eta))
                                continue;
                            if (!(calibrated_jet_pT > b2bjet_pT_min))
                                continue;

                            h_dphi_clusterJets_tight[centbin][pTbin_dphi]->Fill(std::abs(jdphi), weight);
                        }
                    }
                }
                if (
                    cluster_weta_cogx[icluster] > non_tight_weta_cogx_min &&
                    cluster_weta_cogx[icluster] < non_tight_weta_cogx_max &&
                    cluster_wphi_cogx[icluster] > non_tight_wphi_cogx_min &&
                    cluster_wphi_cogx[icluster] < non_tight_wphi_cogx_max &&
                    cluster_prob[icluster] > non_tight_prob_min &&
                    cluster_prob[icluster] < non_tight_prob_max &&
                    e11_over_e33 > non_tight_e11_over_e33_min &&
                    e11_over_e33 < non_tight_e11_over_e33_max &&
                    e32_over_e35 > non_tight_e32_over_e35_min &&
                    e32_over_e35 < non_tight_e32_over_e35_max &&
                    cluster_et1[icluster] > non_tight_et1_min &&
                    cluster_et1[icluster] < non_tight_et1_max &&
                    cluster_et4[icluster] > non_tight_et4_min &&
                    cluster_et4[icluster] < non_tight_et4_max &&
                    bdt_score > non_tight_bdt_min &&
                    bdt_score < non_tight_bdt_max_slope * cluster_Et[icluster] + non_tight_bdt_max_intercept)
                {
                    // fail at least one of the tight cuts with small correlation
                    int nfail = 0;
                    if (!is_cluster_weta_cogx_tight)  nfail += weta_on;
                    if (!is_cluster_wphi_cogx_tight)  nfail += wphi_on;
                    if (!is_cluster_et1_tight)        nfail += et1_on;
                    if (!is_cluster_et2_tight)        nfail += et2_on;
                    if (!is_cluster_et3_tight)        nfail += et3_on;
                    if (!is_e11_over_e33_tight)       nfail += e11_to_e33_on;
                    if (!is_e32_over_e35_tight)       nfail += e32_to_e35_on;
                    if (!is_cluster_et4_tight)        nfail += et4_on;
                    if (!is_cluster_prob_tight)       nfail++;
                    if (!is_bdt_tight)                nfail += bdt_on;

                    if ((nfail > n_nt_fail))
                    {
                        bool all_flags_fail = true;
                        if (weta_fail)
                        {
                            if (is_cluster_weta_cogx_tight)
                                all_flags_fail = false;
                        }
                        if (wphi_fail)
                        {
                            if (is_cluster_wphi_cogx_tight)
                                all_flags_fail = false;
                        }
                        if (et1_fail)
                        {
                            if (is_cluster_et1_tight)
                                all_flags_fail = false;
                        }
                        if (e11_to_e33_fail)
                        {
                            if (is_e11_over_e33_tight)
                                all_flags_fail = false;
                        }
                        if (e32_to_e35_fail)
                        {
                            if (is_e32_over_e35_tight)
                                all_flags_fail = false;
                        }
                        if (bdt_fail)
                        {
                            if (is_bdt_tight)
                                all_flags_fail = false;
                        }

                        if (all_flags_fail)
                        {
                            nontight = true;
                        }
                    }
                }
            }

            auto xj_underflow_x = [](TH2D *h) {
                return h->GetXaxis()->GetBinLowEdge(1) - 1e-4;
            };
            auto fill_reco_xjgamma = [&](TH2D *h) {
                if (!fill_xjgamma || !h) return;
                if (b2bjet_pT.empty())
                    h->Fill(xj_underflow_x(h), cluster_Et[icluster], weight);
                else
                    for (float jpt : b2bjet_pT)
                        h->Fill(jpt / leading_cluster_ET, cluster_Et[icluster], weight);
            };
            auto fill_truth_xjgamma = [&](TH2D *h, float truth_photon_pT) {
                if (!fill_xjgamma || !issim || !h) return;
                if (truth_photon_pT <= 0.f) return;
                if (b2btruthjet_pT.empty())
                    h->Fill(xj_underflow_x(h), truth_photon_pT, weight);
                else
                    for (float tpt : b2btruthjet_pT)
                        h->Fill(tpt / truth_photon_pT, truth_photon_pT, weight);
            };
            auto fill_truthmatched_reco_xjgamma = [&](TH2D *h) {
                if (!fill_xjgamma || !issim || !h) return;
                if (b2bjet_truthmatched_pT.empty())
                    h->Fill(xj_underflow_x(h), cluster_Et[icluster], weight);
                else
                    for (float jpt : b2bjet_truthmatched_pT)
                        h->Fill(jpt / leading_cluster_ET, cluster_Et[icluster], weight);
            };
            auto find_bin = [](float value, const std::vector<float> &edges) {
                for (int ib = 0; ib + 1 < static_cast<int>(edges.size()); ++ib)
                {
                    if (value >= edges[ib] && value < edges[ib + 1])
                        return ib;
                }
                return -1;
            };
            auto in_range = [](float value, const std::vector<float> &edges) {
                if (edges.size() < 2)
                    return false;
                return (value >= edges.front() && value < edges.back());
            };

            if (tight && iso)
            {
                h_tight_iso_cluster[centbin]->Fill(cluster_Et[icluster], weight);
                fill_reco_xjgamma(h_tight_iso_xjgamma[centbin]);
                // h_pT_reco_response[centbin]->Fill(cluster_Et[icluster]);
            }
            if (tight)
            {
                h_tight_cluster[centbin]->Fill(cluster_Et[icluster], weight);
                fill_reco_xjgamma(h_tight_xjgamma[centbin]);
            }
            if (tight && noniso)
            {
                h_tight_noniso_cluster[centbin]->Fill(cluster_Et[icluster], weight);
                fill_reco_xjgamma(h_tight_noniso_xjgamma[centbin]);
            }
            if (nontight && iso)
            {
                h_nontight_iso_cluster[centbin]->Fill(cluster_Et[icluster], weight);
                fill_reco_xjgamma(h_nontight_iso_xjgamma[centbin]);
            }
            if (nontight && noniso)
            {
                h_nontight_noniso_cluster[centbin]->Fill(cluster_Et[icluster], weight);
                fill_reco_xjgamma(h_nontight_noniso_xjgamma[centbin]);
            }
            if (common_pass)
            {
                h_common_cluster[centbin]->Fill(cluster_Et[icluster], weight);
                fill_reco_xjgamma(h_common_xjgamma[centbin]);
                h_cluster_common_Et->Fill(cluster_Et[icluster], weight);
                if (cluster_Et[icluster] > leading_common_cluster_ET)
                {
                    leading_common_cluster_ET = cluster_Et[icluster];
                }
            }
            h_all_cluster[centbin]->Fill(cluster_Et[icluster], weight);
            fill_reco_xjgamma(h_all_xjgamma[centbin]);
            if (pTbin != -1)
            {
                if (tight)
                {
                    h_tight_cluster_pT[centbin][pTbin]->Fill(recoisoET, weight);
                }
                if (nontight)
                {
                    h_nontight_cluster_pT[centbin][pTbin]->Fill(recoisoET, weight);
                }
            }


            if (issim)
            {

                if (particle_trkidmap.find(cluster_truthtrkID[icluster]) == particle_trkidmap.end())
                {
                    continue;
                }
                int iparticle = particle_trkidmap[cluster_truthtrkID[icluster]];

                b2btruthjet_pT.clear();
                if (photon_reco.find(iparticle) != photon_reco.end())
                {
                    for (int ijet = 0; ijet < *njet_truth; ++ijet)
                    {
                        float truth_dphi = particle_Phi[iparticle] - jet_truth_Phi[ijet];
                        while (truth_dphi > M_PI)
                            truth_dphi -= 2 * M_PI;
                        while (truth_dphi < -M_PI)
                            truth_dphi += 2 * M_PI;

                        if (!(std::abs(jet_truth_Eta[ijet]) < jet_eta))
                            continue;
                        if (!(std::abs(truth_dphi) > b2bjet_dphi))
                            continue;
                        if (!(jet_truth_Pt[ijet] > b2bjet_pT_min))
                            continue;

                        b2btruthjet_pT.push_back(jet_truth_Pt[ijet]);
                    }
                    std::sort(b2btruthjet_pT.begin(), b2btruthjet_pT.end());
                }

                // delta R cut
                float deta = cluster_Eta[icluster] - particle_Eta[iparticle];
                float dphi = cluster_Phi[icluster] - particle_Phi[iparticle];
                if (dphi > M_PI)
                {
                    dphi -= 2 * M_PI;
                }
                else if (dphi < -M_PI)
                {
                    dphi += 2 * M_PI;
                }
                float dR = sqrt(deta * deta + dphi * dphi);
                // truth iso vs. reco iso

                if (photon_iso_ET.find(iparticle) != photon_iso_ET.end())
                {
                    if (pTbin != -1)
                    {
                        h_iso_truth_reco[centbin][pTbin]->Fill(photon_iso_ET[iparticle], recoisoET, weight);
                    }
                    if (photon_reco.find(iparticle) == photon_reco.end())
                    {
                        // then it is non truth signal, if it pass the reco, iso, and tight cuts, then it is a fake
                        if (iso && tight && (dR < eff_dR))
                        {
                            h_pT_reco_fake[centbin]->Fill(cluster_Et[icluster], weight);
                        }
                    }
                }

                // iparticle has to be in the photon map
                if (photon_reco.find(iparticle) == photon_reco.end())
                {
                    // this is non-photon background
                    if (pTbin != -1)
                    {
                        h_background_iso_truth_reco[centbin][pTbin]->Fill(all_iso_ET[iparticle], recoisoET, weight);
                    }
                    h_background_truth_isoET[centbin]->Fill(particle_Pt[iparticle], recoisoET, weight);
                    // not matched to a direct or fragmentation photon — fill ABCD background yield
                    if (tight && iso)      h_tight_iso_cluster_notmatch[centbin]->Fill(cluster_Et[icluster], weight);
                    if (tight && noniso)   h_tight_noniso_cluster_notmatch[centbin]->Fill(cluster_Et[icluster], weight);
                    if (nontight && iso)   h_nontight_iso_cluster_notmatch[centbin]->Fill(cluster_Et[icluster], weight);
                    if (nontight && noniso) h_nontight_noniso_cluster_notmatch[centbin]->Fill(cluster_Et[icluster], weight);
                    continue;
                }

                photon_ncluster[iparticle]++;

                if (dR < eff_dR)
                {

                    // if(photon_converts[iparticle]) continue;

                    photon_reco[iparticle] = true;

                    h_pT_truth_reco[centbin]->Fill(particle_Pt[iparticle], cluster_Et[icluster] / particle_Pt[iparticle], weight);
                    if (pTbin != -1)
                    {
                        h_response_isoET[centbin][pTbin]->Fill(cluster_Et[icluster] / particle_Pt[iparticle], recoisoET, weight);
                    }
                    if (iso)
                    {
                        photon_iso[iparticle] = true;
                    }

                    if (tight)
                    {
                        photon_id[iparticle] = true;
                    }

                    if (particle_Pt[iparticle] > pTmin_truth && particle_Pt[iparticle] < pTmax_truth && cluster_Et[icluster] > pTmin && cluster_Et[icluster] < pTmax)
                    {
                        h_all_cluster_signal[centbin]->Fill(cluster_Et[icluster], weight);
                        fill_truthmatched_reco_xjgamma(h_all_truthmatchreco_xjgamma_signal[centbin]);
                        // fill the max backtobjets
                        if (max_b2bjet_pT >= 0.f)
                            h_all_cluster_Et_max_b2bjet[centbin]->Fill(cluster_Et[icluster], max_b2bjet_pT, weight);
                    }

                    if (tight)
                    {
                        if (particle_Pt[iparticle] > pTmin_truth && particle_Pt[iparticle] < pTmax_truth && cluster_Et[icluster] > pTmin && cluster_Et[icluster] < pTmax)
                        {
                            h_tight_cluster_signal[centbin]->Fill(cluster_Et[icluster], weight);
                            fill_reco_xjgamma(h_tight_xjgamma_signal[centbin]);
                            fill_truthmatched_reco_xjgamma(h_tight_truthmatchreco_xjgamma_signal[centbin]);
                            fill_truth_xjgamma(h_tight_truthjet_xjgamma_signal[centbin], particle_Pt[iparticle]);
                            if (iso)
                            {
                                h_tight_iso_cluster_signal[centbin]->Fill(cluster_Et[icluster], weight);
                                fill_reco_xjgamma(h_tight_iso_xjgamma_signal[centbin]);
                                fill_truthmatched_reco_xjgamma(h_tight_iso_truthmatchreco_xjgamma_signal[centbin]);
                                fill_truth_xjgamma(h_tight_iso_truthjet_xjgamma_signal[centbin], particle_Pt[iparticle]);
                                // fill the response matrix

                                float response_reweight = 1.0;
                                if (reweight)
                                {
                                    response_reweight = cluster_Et[icluster] > 30 ? f_reweight->Eval(30) : f_reweight->Eval(cluster_Et[icluster]);
                                }
                                h_pT_truth_response[centbin]->Fill(particle_Pt[iparticle], weight*response_reweight);
                                h_pT_reco_response[centbin]->Fill(cluster_Et[icluster], weight*response_reweight);
                                responses_full[centbin]->Fill(cluster_Et[icluster], particle_Pt[iparticle], weight * response_reweight);
                                h_response_full_list[centbin]->Fill(cluster_Et[icluster], particle_Pt[iparticle], weight * response_reweight);

                                // Build matched pairs in xjgamma by sorting both reco and truth xj lists
                                // in descending order and pairing highest-with-highest.
                                std::vector<float> reco_xj_values;
                                std::vector<float> truth_xj_values;
                                if (leading_cluster_ET > 0.f)
                                {
                                    for (float jpt : b2bjet_pT)
                                        reco_xj_values.push_back(jpt / leading_cluster_ET);
                                }
                                if (particle_Pt[iparticle] > 0.f)
                                {
                                    for (float tpt : b2btruthjet_pT)
                                        truth_xj_values.push_back(tpt / particle_Pt[iparticle]);
                                }

                                std::sort(reco_xj_values.begin(), reco_xj_values.end(), std::greater<float>());
                                std::sort(truth_xj_values.begin(), truth_xj_values.end(), std::greater<float>());

                                const int reco_pt_bin = find_bin(cluster_Et[icluster], pT_bins);
                                const int truth_pt_bin = find_bin(particle_Pt[iparticle], pT_bins_truth);
                                if (reco_pt_bin >= 0 && truth_pt_bin >= 0)
                                {
                                    const int npairs = std::min(static_cast<int>(reco_xj_values.size()), static_cast<int>(truth_xj_values.size()));
                                    for (int im = 0; im < npairs; ++im)
                                    {
                                        const float reco_xj = reco_xj_values[im];
                                        const float truth_xj = truth_xj_values[im];
                                        const int reco_xj_bin = find_bin(reco_xj_values[im], xjgamma_bins);
                                        const int truth_xj_bin = find_bin(truth_xj_values[im], xjgamma_bins_truth);
                                        if (reco_xj_bin < 0 || truth_xj_bin < 0)
                                            continue;

                                        const int reco_global_bin = reco_pt_bin * n_xj_bins + reco_xj_bin + 1;
                                        const int truth_global_bin = truth_pt_bin * n_xj_bins_truth + truth_xj_bin + 1;
                                        h_response_xjgamma_global_list[centbin]->Fill(truth_global_bin, reco_global_bin, weight * response_reweight);
                                        responses_xjgamma_2d[centbin]->Fill(
                                            reco_xj, cluster_Et[icluster],
                                            truth_xj, particle_Pt[iparticle],
                                            weight * response_reweight);
                                    }

                                    // Fill explicit fake/miss entries for unpaired xj values.
                                    for (int im = npairs; im < static_cast<int>(reco_xj_values.size()); ++im)
                                    {
                                        const float reco_xj = reco_xj_values[im];
                                        if (!in_range(reco_xj, xjgamma_bins))
                                            continue;
                                        responses_xjgamma_2d[centbin]->Fake(reco_xj, cluster_Et[icluster], weight * response_reweight);
                                    }
                                    for (int im = npairs; im < static_cast<int>(truth_xj_values.size()); ++im)
                                    {
                                        const float truth_xj = truth_xj_values[im];
                                        if (!in_range(truth_xj, xjgamma_bins_truth))
                                            continue;
                                        responses_xjgamma_2d[centbin]->Miss(truth_xj, particle_Pt[iparticle], weight * response_reweight);
                                    }
                                }
                                if (ientry < (nentries / 2))
                                {
                                    h_pT_truth_half_response[centbin]->Fill(particle_Pt[iparticle], weight);
                                    h_pT_reco_half_response[centbin]->Fill(cluster_Et[icluster], weight);
                                    responses_half[centbin]->Fill(cluster_Et[icluster], particle_Pt[iparticle], weight * response_reweight);
                                    h_response_half_list[centbin]->Fill(cluster_Et[icluster], particle_Pt[iparticle], weight * response_reweight);
                                }
                                else
                                {
                                    h_pT_truth_secondhalf_response[centbin]->Fill(particle_Pt[iparticle], weight);
                                    h_pT_reco_secondhalf_response[centbin]->Fill(cluster_Et[icluster], weight);
                                }
                            }
                        }
                    }
                    if (tight && noniso)
                    {
                        h_tight_noniso_cluster_signal[centbin]->Fill(cluster_Et[icluster], weight);
                        fill_reco_xjgamma(h_tight_noniso_xjgamma_signal[centbin]);
                        fill_truthmatched_reco_xjgamma(h_tight_noniso_truthmatchreco_xjgamma_signal[centbin]);
                        fill_truth_xjgamma(h_tight_noniso_truthjet_xjgamma_signal[centbin], particle_Pt[iparticle]);
                    }
                    if (nontight && iso)
                    {
                        h_nontight_iso_cluster_signal[centbin]->Fill(cluster_Et[icluster], weight);
                        fill_reco_xjgamma(h_nontight_iso_xjgamma_signal[centbin]);
                        fill_truthmatched_reco_xjgamma(h_nontight_iso_truthmatchreco_xjgamma_signal[centbin]);
                        fill_truth_xjgamma(h_nontight_iso_truthjet_xjgamma_signal[centbin], particle_Pt[iparticle]);
                    }
                    if (nontight && noniso)
                    {
                        h_nontight_noniso_cluster_signal[centbin]->Fill(cluster_Et[icluster], weight);
                        fill_reco_xjgamma(h_nontight_noniso_xjgamma_signal[centbin]);
                        fill_truthmatched_reco_xjgamma(h_nontight_noniso_truthmatchreco_xjgamma_signal[centbin]);
                        fill_truth_xjgamma(h_nontight_noniso_truthjet_xjgamma_signal[centbin], particle_Pt[iparticle]);
                    }

                    h_singal_reco_isoET[centbin]->Fill(cluster_Et[icluster], recoisoET, weight);
                    h_singal_truth_isoET[centbin]->Fill(particle_Pt[iparticle], recoisoET, weight);
                }
            }
        } // end of cluster loop
        if (leading_common_cluster_ET > 0)
        {
            h_cluster_common_leading_Et->Fill(leading_common_cluster_ET, weight);
        }

        // go over the map and fill the TEfficiency
        for (auto it = photon_reco.begin(); it != photon_reco.end(); ++it)
        {
            float photon_pT = particle_Pt[it->first];
            float photon_eta = particle_Eta[it->first];
            if (photon_eta <= eta_acceptance_min || photon_eta >= eta_acceptance_max)
            {
                continue;
            }
            if (centbin >= 0 && centbin < n_cent_bins)
                eff_converts_cent[centbin]->Fill(photon_converts[it->first], photon_pT);

            h_ncluster_truth[centbin]->Fill(photon_pT, photon_ncluster[it->first], weight);

            if (centbin >= 0 && centbin < n_cent_bins)
                eff_reco_cent[centbin]->Fill(photon_reco[it->first], photon_pT);


            bool totalpass = photon_reco[it->first] && photon_iso[it->first] && photon_id[it->first];
            if (centbin >= 0 && centbin < n_cent_bins)
                eff_all_cent[centbin]->Fill(totalpass, photon_pT);

            if (photon_reco[it->first])
            {
                // if(!photon_converts[it->first])
                {
                    if (centbin >= 0 && centbin < n_cent_bins)
                        eff_iso_cent[centbin]->Fill(photon_iso[it->first], photon_pT);
                }

                if (photon_iso[it->first])
                {
                    if (centbin >= 0 && centbin < n_cent_bins)
                        eff_id_cent[centbin]->Fill(photon_id[it->first], photon_pT);
                }
            }
        }
        ientry++;
    }

    TFile *fresponse = new TFile(responsefilename.c_str(), "RECREATE");
    for (int icent = 0; icent < n_cent_bins; icent++)
    {

        responses_full[icent]->Write();

        responses_half[icent]->Write();
        responses_xjgamma_2d[icent]->Write();

        h_pT_truth_response[icent]->Write();
        h_pT_reco_response[icent]->Write();
        h_xjgamma_truth_response[icent]->Write();
        h_xjgamma_reco_response[icent]->Write();

        h_pT_truth_half_response[icent]->Write();
        h_pT_reco_half_response[icent]->Write();

        h_pT_truth_secondhalf_response[icent]->Write();
        h_pT_reco_secondhalf_response[icent]->Write();
        h_response_xjgamma_global_list[icent]->Write();
    }

    fout->Write();
    if (!issim)
      SaveYamlToRoot(fout, configname.c_str());
    fout->Close();

    fresponse->Write();
    fresponse->Close();
}
