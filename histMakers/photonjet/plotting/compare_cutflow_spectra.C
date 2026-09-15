// compare_cutflow_spectra.C -- photon-cluster ET spectra at every selection level for two
// var_types (e.g. the PhotonJetTrees maker vs PPG12 RecoEff on the slimtree), with ratios.
//
//   root -l -b -q 'compare_cutflow_spectra.C("pj_pp_cert","bdt_pjcert")'
//
// Levels (data, signal MC, inclusive MC): all clusters, common, tight, A=tight+iso,
// B=tight+noniso, C=nontight+iso, D=nontight+noniso, plus the truth-matched `_signal`
// versions for the signal MC. Prints per-bin ratios a/b and the step efficiencies
// common/all, tight/common, A/tight for each code. Figures -> photonjet/figures/.
#include "PhotonJetPlotInput.h"
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>
#include <cstdio>
#include <string>
#include <vector>

namespace
{
const char *R = "results/";
struct Level { const char *label; const char *hname; };

TH1 *Get(TFile *f, const char *name)
{
    TH1 *h = f ? dynamic_cast<TH1 *>(f->Get(name)) : nullptr;
    if (h) h->SetDirectory(nullptr);
    return h;
}

void PrintTable(const char *title, TFile *fa, TFile *fb, const std::vector<Level> &levels, const char *a, const char *b)
{
    printf("\n== %s: %s / %s ==\n", title, a, b);
    TH1 *ref = Get(fa, levels[0].hname);
    if (!ref) { printf("  (missing %s)\n", levels[0].hname); return; }
    printf("%-9s", "bin");
    for (const auto &l : levels) printf(" | %-26s", l.label);
    printf("\n");
    for (int bin = 1; bin <= ref->GetNbinsX(); ++bin)
    {
        printf("%3.0f-%-5.0f", ref->GetBinLowEdge(bin), ref->GetBinLowEdge(bin + 1));
        for (const auto &l : levels)
        {
            TH1 *ha = Get(fa, l.hname), *hb = Get(fb, l.hname);
            if (!ha || !hb) { printf(" | %-26s", "missing"); continue; }
            const double va = ha->GetBinContent(bin), vb = hb->GetBinContent(bin);
            printf(" | %9.4g %9.4g r=%.3f", va, vb, vb > 0 ? va / vb : 0.0);
        }
        printf("\n");
    }
    printf("%-9s", "16-36 sum");
    for (const auto &l : levels)
    {
        TH1 *ha = Get(fa, l.hname), *hb = Get(fb, l.hname);
        if (!ha || !hb) { printf(" | %-26s", "missing"); continue; }
        const double va = ha->Integral(), vb = hb->Integral();
        printf(" | %9.4g %9.4g r=%.3f", va, vb, vb > 0 ? va / vb : 0.0);
    }
    printf("\n");
}

void PrintSteps(const char *title, TFile *f, const char *tag)
{
    TH1 *all = Get(f, "h_all_cluster_0"), *com = Get(f, "h_common_cluster_0"), *ti = Get(f, "h_tight_cluster_0"), *A = Get(f, "h_tight_iso_cluster_0");
    if (!all || !com || !ti || !A) return;
    printf("  %-10s %-11s common/all  tight/common  A/tight   (per bin)\n", title, tag);
    for (int bin = 1; bin <= all->GetNbinsX(); ++bin)
        printf("  %-10s %3.0f-%-5.0f   %.4f      %.4f        %.4f\n", "", all->GetBinLowEdge(bin), all->GetBinLowEdge(bin + 1),
               all->GetBinContent(bin) > 0 ? com->GetBinContent(bin) / all->GetBinContent(bin) : 0,
               com->GetBinContent(bin) > 0 ? ti->GetBinContent(bin) / com->GetBinContent(bin) : 0,
               ti->GetBinContent(bin) > 0 ? A->GetBinContent(bin) / ti->GetBinContent(bin) : 0);
}

void Draw(const char *title, TFile *fa, TFile *fb, const std::vector<Level> &levels, const char *a, const char *b, const char *outpdf)
{
    const int n = (int)levels.size();
    TCanvas *c = new TCanvas(Form("c_%s", title), title, 420 * n, 760);
    c->Divide(n, 2, 0.002, 0.002);
    std::vector<TObject *> keep;
    for (int i = 0; i < n; ++i)
    {
        TH1 *ha = Get(fa, levels[i].hname), *hb = Get(fb, levels[i].hname);
        if (!ha || !hb) continue;
        ha = (TH1 *)ha->Clone(Form("a_%s_%d", title, i)); hb = (TH1 *)hb->Clone(Form("b_%s_%d", title, i));
        const int lastclean = ha->GetXaxis()->FindBin(31.9); // base_E boundary at 35 touches the 32-36 bin
        ha->GetXaxis()->SetRange(1, lastclean); hb->GetXaxis()->SetRange(1, lastclean);
        keep.push_back(ha); keep.push_back(hb);
        c->cd(i + 1); gPad->SetLogy(); gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.13);
        hb->SetLineColor(kBlack); hb->SetMarkerColor(kBlack); hb->SetMarkerStyle(24); hb->SetLineWidth(2);
        ha->SetLineColor(kRed + 1); ha->SetMarkerColor(kRed + 1); ha->SetMarkerStyle(20);
        hb->SetTitle(Form("%s;E_{T} [GeV];counts (weighted)", levels[i].label));
        hb->GetYaxis()->SetTitleOffset(1.6);
        const double mx = std::max(ha->GetMaximum(), hb->GetMaximum());
        hb->SetMaximum(mx * 5); hb->SetMinimum(std::max(1e-3, std::min(ha->GetMinimum(1e-9), hb->GetMinimum(1e-9)) * 0.3));
        hb->Draw("hist e"); ha->Draw("p e same");
        if (i == 0)
        {
            TLegend *leg = new TLegend(0.35, 0.72, 0.92, 0.88); leg->SetBorderSize(0); leg->SetFillStyle(0);
            leg->AddEntry(ha, "PPG19 tree", "pl"); leg->AddEntry(hb, "PPG12 tree", "l"); leg->Draw(); keep.push_back(leg);
        }
        c->cd(n + i + 1); gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.13); gPad->SetGridy();
        TH1 *r = (TH1 *)ha->Clone(Form("r_%s_%d", title, i)); r->Divide(hb); keep.push_back(r);
        r->GetXaxis()->SetRange(1, lastclean);
        for (int bb = 0; bb <= r->GetNbinsX() + 1; ++bb) r->SetBinError(bb, 0); // same events on both sides, correlated
        r->SetTitle(";E_{T} [GeV];PPG19 / PPG12"); r->SetMinimum(0.8); r->SetMaximum(1.2);
        r->SetLineColor(kRed + 1); r->SetMarkerColor(kRed + 1); r->GetYaxis()->SetTitleOffset(1.6);
        r->Draw("p");
        TLine *l1 = new TLine(r->GetXaxis()->GetXmin(), 1, r->GetXaxis()->GetXmax(), 1); l1->SetLineStyle(2); l1->Draw(); keep.push_back(l1);
    }
    c->SaveAs(outpdf);
    printf("  wrote %s\n", outpdf);
}
} // namespace

void compare_cutflow_spectra(const char *a = "pj_pp_cert", const char *b = "bdt_pjcert")
{
    gStyle->SetOptStat(0); gStyle->SetTitleFontSize(0.06);
    gSystem->mkdir("/sphenix/user/shuhangli/ppg12/efficiencytool/photonjet/figures", true);
    const std::vector<Level> reco = {{"all clusters", "h_all_cluster_0"}, {"common", "h_common_cluster_0"}, {"tight", "h_tight_cluster_0"},
                                     {"A tight+iso", "h_tight_iso_cluster_0"}, {"B tight+noniso", "h_tight_noniso_cluster_0"},
                                     {"C nontight+iso", "h_nontight_iso_cluster_0"}, {"D nontight+noniso", "h_nontight_noniso_cluster_0"}};
    // figures show all clusters, tight, and the signal region A; common (differs by
    // construction, prob cut) and the B/C/D control regions stay in the printed tables only
    const std::vector<Level> reco_draw = {reco[0], reco[2], reco[3]};
    const std::vector<Level> sig = {{"signal all", "h_all_cluster_signal_0"}, {"signal tight", "h_tight_cluster_signal_0"},
                                    {"signal A", "h_tight_iso_cluster_signal_0"}, {"signal B", "h_tight_noniso_cluster_signal_0"},
                                    {"signal C", "h_nontight_iso_cluster_signal_0"}, {"signal D", "h_nontight_noniso_cluster_signal_0"}};
    struct Set { const char *title; std::string fa, fb; std::vector<Level> levels; const char *pdf; };
    std::vector<Set> sets = {
        {"data", Form("%sdata_histo_%s.root", R, a), Form("%sdata_histo_%s.root", R, b), reco, "data"},
        {"signal_MC_reco", Form("%sMC_efficiency_%s.root", R, a), Form("%sMC_efficiency_%s.root", R, b), reco, "signal_reco"},
        {"signal_MC_matched", Form("%sMC_efficiency_%s.root", R, a), Form("%sMC_efficiency_%s.root", R, b), sig, "signal_matched"},
        {"inclusive_MC", Form("%sMC_efficiency_jet_%s.root", R, a), Form("%sMC_efficiency_jet_%s.root", R, b), reco, "inclusive"}};
    for (auto &s : sets)
    {
        TFile *fa = PJPlot::Open(s.fa.c_str()), *fb = PJPlot::Open(s.fb.c_str());
        if (!fa || fa->IsZombie() || !fb || fb->IsZombie()) { printf("cannot open %s or %s\n", s.fa.c_str(), s.fb.c_str()); continue; }
        PrintTable(s.title, fa, fb, s.levels, a, b);
        const bool isreco = s.levels[0].hname == std::string("h_all_cluster_0");
        if (isreco) { PrintSteps(s.title, fa, a); PrintSteps(s.title, fb, b); }
        Draw(s.title, fa, fb, isreco ? reco_draw : s.levels, a, b, Form("figures/cutflow_%s_%s_vs_%s.pdf", s.pdf, a, b));
        fa->Close(); fb->Close();
    }
}
