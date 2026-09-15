// plot_pj_summary.C -- summary figures for the PhotonJetTrees analysis outputs.
//
//   root -l -b -q 'plot_pj_summary.C("auau")'   // three centrality classes
//   root -l -b -q 'plot_pj_summary.C("pp")'     // pj_pp_nom vs PPG12 bdt_pjlike and bdt_nom
//
// Au+Au: figures/auau_abcd_spectra.pdf (data A/B/C/D per class), auau_purity.pdf,
//        auau_efficiency.pdf (eps_reco/iso/id per class), auau_leakage.pdf (cB/cC/cD),
//        auau_qa.pdf (centrality, vertex, weights, threshold recipe of selected data/MC).
// pp:    pp_xsec_compare.pdf (final result + ratio), pp_purity_compare.pdf, pp_eff_compare.pdf.
// Absolute Au+Au normalization uses the placeholder lumi = 1, so only shapes and
// purities/efficiencies are meaningful there.
#include "PhotonJetPlotInput.h"
#include <TH1.h>
#include <TEfficiency.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <cstdio>
#include <string>
#include <vector>

namespace
{
const char *R = "results/";
const char *FIG = "figures/";
const int cols[4] = {kBlack, kRed + 1, kBlue + 1, kGreen + 2};
const int mks[4] = {20, 21, 22, 23};

TObject *Get(const std::string &file, const char *name)
{
    TFile *f = PJPlot::Open(file.c_str(), "READ");
    if (!f || f->IsZombie()) { printf("cannot open %s\n", file.c_str()); return nullptr; }
    TObject *o = f->Get(name);
    if (!o) { printf("missing %s in %s\n", name, file.c_str()); f->Close(); return nullptr; }
    TObject *c = o->Clone(Form("%s_%s", name, gSystem->BaseName(file.c_str())));
    if (TH1 *h = dynamic_cast<TH1 *>(c)) h->SetDirectory(nullptr);
    f->Close();
    return c;
}
TH1 *GetH(const std::string &file, const char *name) { return dynamic_cast<TH1 *>(Get(file, name)); }

TLatex *Label(double x, double y, const char *txt, double size = 0.045)
{
    TLatex *t = new TLatex(x, y, txt); t->SetNDC(); t->SetTextSize(size); t->Draw(); return t;
}

TGraphAsymmErrors *EffGraph(const std::string &file, const char *name, int col, int mk)
{
    TEfficiency *e = dynamic_cast<TEfficiency *>(Get(file, name));
    if (!e) return nullptr;
    TGraphAsymmErrors *g = e->CreateGraph();
    g->SetLineColor(col); g->SetMarkerColor(col); g->SetMarkerStyle(mk); g->SetMarkerSize(1.1);
    return g;
}
} // namespace

void plot_pj_auau()
{
    const std::vector<std::string> cls = {"pj_auau_c00_20", "pj_auau_c20_50", "pj_auau_c50_80"};
    const std::vector<std::string> lab = {"0-20%", "20-50%", "50-80%"};
    const char *sys = "Au+Au #sqrt{s_{NN}} = 200 GeV, PhotonJetTrees_v1, |z_{vtx}| < 30 cm";

    // ---- ABCD data spectra per class ----------------------------------------
    TCanvas *c1 = new TCanvas("c_auau_abcd", "", 1500, 500);
    c1->Divide(3, 1, 0.002, 0.002);
    const char *reg[4] = {"h_tight_iso_cluster_0", "h_tight_noniso_cluster_0", "h_nontight_iso_cluster_0", "h_nontight_noniso_cluster_0"};
    const char *regl[4] = {"A: tight, isolated", "B: tight, non-isolated", "C: non-tight, isolated", "D: non-tight, non-isolated"};
    for (size_t i = 0; i < cls.size(); ++i)
    {
        c1->cd(i + 1); gPad->SetLogy(); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
        TLegend *leg = new TLegend(0.45, 0.62, 0.92, 0.88); leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.035);
        for (int r = 0; r < 4; ++r)
        {
            TH1 *h = GetH(Form("%sdata_histo_%s.root", R, cls[i].c_str()), reg[r]);
            if (!h) continue;
            h->SetLineColor(cols[r]); h->SetMarkerColor(cols[r]); h->SetMarkerStyle(mks[r]); h->SetLineWidth(2);
            h->SetTitle(Form("Au+Au data %s;E_{T}^{#gamma} [GeV];photon candidates", lab[i].c_str()));
            h->GetYaxis()->SetTitleOffset(1.5); h->SetMinimum(0.5); h->SetMaximum(3e5);
            h->Draw(r == 0 ? "p e" : "p e same");
            leg->AddEntry(h, Form("%s (%.0f)", regl[r], h->Integral()), "pl");
        }
        leg->Draw();
        if (i == 0) Label(0.16, 0.92, sys, 0.035);
    }
    c1->SaveAs(Form("%sauau_abcd_spectra.pdf", FIG));

    // ---- purity per class ------------------------------------------------------
    TCanvas *c2 = new TCanvas("c_auau_purity", "", 700, 600);
    gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
    TH1 *fr = new TH1F("fr_pur", ";E_{T}^{#gamma} [GeV];purity (leakage corrected)", 10, 16, 36); fr->SetMinimum(0); fr->SetMaximum(1.3); fr->Draw();
    TLegend *lp = new TLegend(0.16, 0.7, 0.55, 0.88); lp->SetBorderSize(0); lp->SetFillStyle(0);
    for (size_t i = 0; i < cls.size(); ++i)
    {
        TGraphErrors *g = dynamic_cast<TGraphErrors *>(Get(Form("%sPhoton_final_%s.root", R, cls[i].c_str()), "gpurity_leak"));
        if (!g) continue;
        g->SetLineColor(cols[i + 1]); g->SetMarkerColor(cols[i + 1]); g->SetMarkerStyle(mks[i]); g->SetMarkerSize(1.2);
        g->Draw("p e same"); lp->AddEntry(g, lab[i].c_str(), "pl");
    }
    TLine *l1 = new TLine(16, 1, 36, 1); l1->SetLineStyle(2); l1->Draw();
    lp->Draw(); Label(0.15, 0.92, sys, 0.035);
    c2->SaveAs(Form("%sauau_purity.pdf", FIG));

    // ---- efficiencies per class ---------------------------------------------------
    TCanvas *c3 = new TCanvas("c_auau_eff", "", 1500, 500);
    c3->Divide(3, 1, 0.002, 0.002);
    const char *en[3] = {"eff_reco_eta_0", "eff_iso_eta_0", "eff_id_eta_0"};
    const char *et[3] = {"#varepsilon_{reco}", "#varepsilon_{iso} (given reco)", "#varepsilon_{ID} (given reco, iso)"};
    for (int k = 0; k < 3; ++k)
    {
        c3->cd(k + 1); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
        TH1 *f = new TH1F(Form("fr_eff%d", k), Form(";p_{T}^{#gamma, truth} [GeV];%s", et[k]), 10, 8, 45); f->SetMinimum(0); f->SetMaximum(1.15); f->Draw();
        TLegend *le = new TLegend(0.18, 0.18, 0.5, 0.38); le->SetBorderSize(0); le->SetFillStyle(0);
        for (size_t i = 0; i < cls.size(); ++i)
        {
            TGraphAsymmErrors *g = EffGraph(Form("%sMC_efficiency_%s.root", R, cls[i].c_str()), en[k], cols[i + 1], mks[i]);
            if (!g) continue;
            g->Draw("p e same"); le->AddEntry(g, lab[i].c_str(), "pl");
        }
        le->Draw();
        if (k == 0) Label(0.16, 0.92, "embedded photon MC, thresholds from tree", 0.035);
    }
    c3->SaveAs(Form("%sauau_efficiency.pdf", FIG));

    // ---- signal leakage per class --------------------------------------------------
    TCanvas *c4 = new TCanvas("c_auau_leak", "", 1500, 500);
    c4->Divide(3, 1, 0.002, 0.002);
    const char *ln[3] = {"h_tight_noniso_cluster_signal_0", "h_nontight_iso_cluster_signal_0", "h_nontight_noniso_cluster_signal_0"};
    const char *lt[3] = {"c_{B} = B_{sig}/A_{sig}", "c_{C} = C_{sig}/A_{sig}", "c_{D} = D_{sig}/A_{sig}"};
    for (int k = 0; k < 3; ++k)
    {
        c4->cd(k + 1); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
        TH1 *f = new TH1F(Form("fr_leak%d", k), Form(";E_{T}^{#gamma} [GeV];%s", lt[k]), 10, 16, 36); f->SetMinimum(0); f->SetMaximum(k == 1 ? 1.0 : 0.4); f->Draw();
        TLegend *le = new TLegend(0.18, 0.65, 0.5, 0.88); le->SetBorderSize(0); le->SetFillStyle(0);
        for (size_t i = 0; i < cls.size(); ++i)
        {
            TH1 *a = GetH(Form("%sMC_efficiency_%s.root", R, cls[i].c_str()), "h_tight_iso_cluster_signal_0");
            TH1 *x = GetH(Form("%sMC_efficiency_%s.root", R, cls[i].c_str()), ln[k]);
            if (!a || !x) continue;
            x->Divide(a); x->SetLineColor(cols[i + 1]); x->SetMarkerColor(cols[i + 1]); x->SetMarkerStyle(mks[i]);
            x->Draw("p e same"); le->AddEntry(x, lab[i].c_str(), "pl");
        }
        le->Draw();
    }
    c4->SaveAs(Form("%sauau_leakage.pdf", FIG));

    // ---- QA: centrality / vertex / weights / threshold recipe ---------------------------
    TCanvas *c5 = new TCanvas("c_auau_qa", "", 1500, 900);
    c5->Divide(3, 2, 0.002, 0.002);
    const char *qa[6] = {"h_pj_centrality", "h_pj_vertex_z", "h_pj_cutflow", "h_pj_centrality", "h_pj_event_weight_log10", "h_pj_auau_threshold_recipe"};
    const char *qt[6] = {"data: centrality of selected events", "data: vertex z of selected events", "data: cutflow", "signal MC: centrality of selected events", "signal MC: log10(event weight)", "signal MC: iso-threshold recipe (0 table, 1 small)"};
    for (int k = 0; k < 6; ++k)
    {
        c5->cd(k + 1); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
        TLegend *le = new TLegend(0.55, 0.65, 0.9, 0.88); le->SetBorderSize(0); le->SetFillStyle(0);
        bool first = true;
        for (size_t i = 0; i < cls.size(); ++i)
        {
            const std::string file = (k < 3) ? Form("%sdata_histo_%s.root", R, cls[i].c_str()) : Form("%sMC_efficiency_%s.root", R, cls[i].c_str());
            TH1 *h = GetH(file, qa[k]);
            if (!h) continue;
            h->SetLineColor(cols[i + 1]); h->SetLineWidth(2); h->SetTitle(qt[k]);
            if (k == 2) { gPad->SetLogy(); h->GetXaxis()->SetLabelSize(0.04); }
            h->Draw(first ? "hist" : "hist same"); first = false; le->AddEntry(h, lab[i].c_str(), "l");
        }
        le->Draw();
    }
    c5->SaveAs(Form("%sauau_qa.pdf", FIG));
}

void plot_pj_pp()
{
    // all three at the same 30 cm vertex cut and 16+ GeV binning: PPG19 (SI by necessity),
    // PPG12 with SI-only MC, PPG12 with the full nominal treatment (SI+DI blend + reweight)
    const std::vector<std::string> vars = {"pj_pp_nom", "bdt_pjcert", "bdt_pjlike"};
    const std::vector<std::string> lab = {"PPG19 tree", "PPG12 tree, SI only", "PPG12 tree, full treatment"};
    // ---- final result + ratio -------------------------------------------------------
    TCanvas *c1 = new TCanvas("c_pp_xsec", "", 700, 800);
    TPad *p1 = new TPad("p1", "", 0, 0.35, 1, 1); TPad *p2 = new TPad("p2", "", 0, 0, 1, 0.35);
    p1->SetBottomMargin(0.02); p2->SetTopMargin(0.03); p2->SetBottomMargin(0.3); p1->Draw(); p2->Draw();
    p1->cd(); gPad->SetLogy(); gPad->SetLeftMargin(0.14);
    TLegend *lg = new TLegend(0.50, 0.66, 0.93, 0.89); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.032);
    std::vector<TH1 *> hs;
    // bdt_pjcert corrects with the unreweighted-SI eps_MBD+vtx (vertexcut_mbd_cut/vertexcut,
    // 0.23 at 17 GeV) while PPG19 and bdt_pjlike share the truth-vertex-reweighted one (0.40).
    // The raw-sim vertex profile is not data-like inside 30 cm, so for an apples-to-apples
    // overlay the SI-only series is rescaled per truth bin onto the common (pjlike) factor.
    TH1 *e_own_n = GetH(Form("%sMC_efficiency_bdt_pjcert.root", R), "h_truth_pT_vertexcut_mbd_cut_0");
    TH1 *e_own_d = GetH(Form("%sMC_efficiency_bdt_pjcert.root", R), "h_truth_pT_vertexcut_0");
    TH1 *e_com_n = GetH(Form("%sMC_efficiency_bdt_pjlike.root", R), "h_truth_pT_vertexcut_mbd_cut_0");
    TH1 *e_com_d = GetH(Form("%sMC_efficiency_bdt_pjlike.root", R), "h_truth_pT_vertexcut_0");
    for (size_t i = 0; i < vars.size(); ++i)
    {
        TH1 *h = GetH(Form("%sPhoton_final_%s.root", R, vars[i].c_str()), "h_unfold_sub_result");
        hs.push_back(h);
        if (!h) continue;
        if (vars[i] == "bdt_pjcert" && e_own_n && e_own_d && e_com_n && e_com_d)
            for (int b = 1; b <= h->GetNbinsX(); ++b)
            {
                const double x = h->GetBinCenter(b);
                const int bo = e_own_d->GetXaxis()->FindBin(x), bc = e_com_d->GetXaxis()->FindBin(x);
                if (e_own_d->GetBinContent(bo) <= 0 || e_com_d->GetBinContent(bc) <= 0 || e_com_n->GetBinContent(bc) <= 0) continue;
                const double eps_own = e_own_n->GetBinContent(bo) / e_own_d->GetBinContent(bo);
                const double eps_com = e_com_n->GetBinContent(bc) / e_com_d->GetBinContent(bc);
                h->SetBinContent(b, h->GetBinContent(b) * eps_own / eps_com);
                h->SetBinError(b, h->GetBinError(b) * eps_own / eps_com);
            }
        h->SetLineColor(cols[i]); h->SetMarkerColor(cols[i]); h->SetMarkerStyle(mks[i]); h->SetMarkerSize(1.1);
        // common comparable window: maker/pjlike acceptance starts at 16 (tree skim), their
        // unfolded truth bins below 16 are migration artifacts, not measurements
        h->SetTitle(";E_{T}^{#gamma} [GeV];d^{2}#sigma/dE_{T}d#eta [pb/GeV]"); h->GetXaxis()->SetRangeUser(16, 32);
        h->Draw(i == 0 ? "p e" : "p e same"); lg->AddEntry(h, lab[i].c_str(), "pl");
    }
    lg->Draw(); Label(0.16, 0.92, "#it{p}+#it{p} #sqrt{s} = 200 GeV, |#eta^{#gamma}| < 0.7, isolated photons", 0.035);
    Label(0.5, 0.62, "SI-only rescaled to common #varepsilon_{MBD+vtx}", 0.028);
    // y range must cover all three series (bdt_nom starts at 10 GeV and reaches ~600 pb/GeV),
    // not just the first-drawn maker histogram; top pad keeps no x labels (ratio pad has them)
    double gmax = 0, gmin = 1e30;
    for (auto *h : hs)
    {
        if (!h) continue;
        for (int b = h->GetXaxis()->FindBin(16.1); b <= h->GetXaxis()->FindBin(31.9); ++b)
        {
            const double v = h->GetBinContent(b);
            if (v > 0 && v > gmax) gmax = v;
            if (v > 0 && v < gmin) gmin = v;
        }
    }
    if (hs[0] && gmax > 0) { hs[0]->SetMaximum(gmax * 4); hs[0]->SetMinimum(gmin * 0.3); hs[0]->GetXaxis()->SetLabelSize(0); }
    p2->cd(); gPad->SetLeftMargin(0.14); gPad->SetGridy();
    for (size_t i = 1; i < vars.size(); ++i)
    {
        if (!hs[0] || !hs[i]) continue;
        TH1 *r = (TH1 *)hs[0]->Clone(Form("r_%zu", i)); r->Divide(hs[i]);
        r->SetLineColor(cols[i]); r->SetMarkerColor(cols[i]); r->SetMarkerStyle(mks[i]);
        r->SetTitle(";E_{T}^{#gamma} [GeV];PPG19 / PPG12"); r->SetMinimum(0.55); r->SetMaximum(1.4);
        r->GetXaxis()->SetRangeUser(16, 32); r->GetXaxis()->SetLabelSize(0.08); r->GetXaxis()->SetTitleSize(0.09); r->GetYaxis()->SetLabelSize(0.07); r->GetYaxis()->SetTitleSize(0.08); r->GetYaxis()->SetTitleOffset(0.7);
        r->Draw(i == 1 ? "p e" : "p e same");
    }
    TLine *l1 = new TLine(16, 1, 32, 1); l1->SetLineStyle(2); l1->Draw();
    c1->SaveAs(Form("%spp_xsec_compare.pdf", FIG));

    // ---- purity ---------------------------------------------------------------------
    TCanvas *c2 = new TCanvas("c_pp_pur", "", 700, 600); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
    TH1 *fr = new TH1F("fr_ppur", ";E_{T}^{#gamma} [GeV];purity (leakage corrected)", 10, 10, 32); fr->SetMinimum(0.3); fr->SetMaximum(1.3); fr->Draw();
    TLegend *lp = new TLegend(0.15, 0.7, 0.7, 0.88); lp->SetBorderSize(0); lp->SetFillStyle(0); lp->SetTextSize(0.035);
    for (size_t i = 0; i < vars.size(); ++i)
    {
        TGraphErrors *g = dynamic_cast<TGraphErrors *>(Get(Form("%sPhoton_final_%s.root", R, vars[i].c_str()), "gpurity_leak"));
        if (!g) continue;
        g->SetLineColor(cols[i]); g->SetMarkerColor(cols[i]); g->SetMarkerStyle(mks[i]); g->Draw("p e same"); lp->AddEntry(g, lab[i].c_str(), "pl");
    }
    lp->Draw(); c2->SaveAs(Form("%spp_purity_compare.pdf", FIG));

    // ---- efficiencies -----------------------------------------------------------------
    TCanvas *c3 = new TCanvas("c_pp_eff", "", 1500, 500); c3->Divide(3, 1, 0.002, 0.002);
    const char *en[3] = {"eff_reco_eta_0", "eff_iso_eta_0", "eff_id_eta_0"};
    const char *et[3] = {"#varepsilon_{reco}", "#varepsilon_{iso} (given reco)", "#varepsilon_{ID} (given reco, iso)"};
    for (int k = 0; k < 3; ++k)
    {
        c3->cd(k + 1); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
        TH1 *f = new TH1F(Form("fr_ppeff%d", k), Form(";p_{T}^{#gamma, truth} [GeV];%s", et[k]), 10, 8, 32); f->SetMinimum(0.3); f->SetMaximum(1.15); f->Draw();
        TLegend *le = new TLegend(0.16, 0.16, 0.7, 0.36); le->SetBorderSize(0); le->SetFillStyle(0); le->SetTextSize(0.033);
        for (size_t i = 0; i < vars.size(); ++i)
        {
            TGraphAsymmErrors *g = EffGraph(Form("%sMC_efficiency_%s.root", R, vars[i].c_str()), en[k], cols[i], mks[i]);
            if (!g) continue;
            g->Draw("p e same"); le->AddEntry(g, lab[i].c_str(), "pl");
        }
        le->Draw();
    }
    c3->SaveAs(Form("%spp_eff_compare.pdf", FIG));
}

void plot_pj_summary(const char *system = "auau")
{
    gStyle->SetOptStat(0); gStyle->SetOptTitle(1); gStyle->SetTitleFontSize(0.045);
    gSystem->mkdir(FIG, true);
    if (std::string(system) == "auau") plot_pj_auau();
    else plot_pj_pp();
}
