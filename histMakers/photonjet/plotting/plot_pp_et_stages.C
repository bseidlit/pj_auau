// plot_pp_et_stages.C -- pp data photon ET spectrum, ONE figure per selection stage,
// new-tree maker (pj_pp_nom, filled red) vs PPG12 RecoEff on the slimtree
// (bdt_pjlike, open black), with a maker/PPG12 ratio panel.
//   root -l -b -q plot_pp_et_stages.C
// Output: figures/pp_et_stage_{all,common,tight,A_tight_iso,B_tight_noniso,
//         C_nontight_iso,D_nontight_noniso}.pdf
#include "PhotonJetPlotInput.h"
#include <TH1.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <cstdio>

void plot_pp_et_stages(const char *a = "pj_pp_nom", const char *b = "bdt_pjlike")
{
  gStyle->SetOptStat(0);
  const char *R = "results/";
  const char *FIG = "figures/";
  TFile *fa = PJPlot::Open(Form("%sdata_histo_%s.root", R, a));
  TFile *fb = PJPlot::Open(Form("%sdata_histo_%s.root", R, b));
  const int N = 7;
  const char *hn[N] = {"h_all_cluster_0", "h_common_cluster_0", "h_tight_cluster_0", "h_tight_iso_cluster_0",
                       "h_tight_noniso_cluster_0", "h_nontight_iso_cluster_0", "h_nontight_noniso_cluster_0"};
  const char *tag[N] = {"all", "common", "tight", "A_tight_iso", "B_tight_noniso", "C_nontight_iso", "D_nontight_noniso"};
  const char *ttl[N] = {"all clusters (E_{T} #geq 15 GeV)", "common cuts", "tight",
                        "A: tight + isolated", "B: tight + non-isolated", "C: non-tight + isolated", "D: non-tight + non-isolated"};
  for (int i = 0; i < N; ++i)
  {
    TH1 *ha = (TH1 *)fa->Get(hn[i]);
    TH1 *hb = (TH1 *)fb->Get(hn[i]);
    if (!ha || !hb) { printf("missing %s\n", hn[i]); continue; }
    ha = (TH1 *)ha->Clone(); hb = (TH1 *)hb->Clone();
    // comparison shown only below 32 GeV: PPG12 switches to the base_E model at
    // ET >= 35 GeV (not stored in the tree), which touches the 32-36 bin
    const int lastclean = ha->GetXaxis()->FindBin(31.9);
    ha->GetXaxis()->SetRange(1, lastclean); hb->GetXaxis()->SetRange(1, lastclean);
    TCanvas *c = new TCanvas(Form("c_%s", tag[i]), "", 700, 750);
    TPad *p1 = new TPad("p1", "", 0, 0.32, 1, 1);
    TPad *p2 = new TPad("p2", "", 0, 0, 1, 0.32);
    p1->SetBottomMargin(0.025); p1->SetLeftMargin(0.14);
    p2->SetTopMargin(0.04); p2->SetBottomMargin(0.32); p2->SetLeftMargin(0.14); p2->SetGridy();
    p1->Draw(); p2->Draw();

    p1->cd(); p1->SetLogy();
    hb->SetLineColor(kBlack); hb->SetMarkerColor(kBlack); hb->SetMarkerStyle(24); hb->SetMarkerSize(1.2); hb->SetLineWidth(2);
    ha->SetLineColor(kRed + 1); ha->SetMarkerColor(kRed + 1); ha->SetMarkerStyle(20); ha->SetMarkerSize(1.1);
    hb->SetTitle(Form(";;photon candidates / bin (weighted by 1/#varepsilon_{L1})"));
    hb->GetXaxis()->SetLabelSize(0); hb->GetYaxis()->SetTitleOffset(1.35);
    double mx = TMath::Max(ha->GetMaximum(), hb->GetMaximum());
    double mn = TMath::Min(ha->GetMinimum(1e-9), hb->GetMinimum(1e-9));
    hb->SetMaximum(mx * 8); hb->SetMinimum(TMath::Max(0.5, mn * 0.4));
    hb->Draw("p e"); ha->Draw("p e same");
    TLegend *leg = new TLegend(0.45, 0.62, 0.93, 0.80); leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.038);
    leg->AddEntry(ha, Form("PPG19 tree (%.0f)", ha->Integral(1, lastclean)), "pl");
    leg->AddEntry(hb, Form("PPG12 tree (%.0f)", hb->Integral(1, lastclean)), "pl");
    leg->Draw();
    TLatex t; t.SetNDC();
    t.SetTextSize(0.045); t.SetTextFont(62); t.DrawLatex(0.17, 0.86, ttl[i]);
    t.SetTextSize(0.036); t.SetTextFont(42);
    t.DrawLatex(0.17, 0.81, "#it{p}+#it{p} #sqrt{s} = 200 GeV data, |#eta^{#gamma}| < 0.7, |z_{vtx}| < 30 cm");
    if (i == 1) { t.SetTextColor(kGray + 2); t.DrawLatex(0.17, 0.25, "#splitline{expected to differ: PPG12 common includes prob #in (0,1),}{cluster prob is not stored in PhotonJetTrees}"); t.SetTextColor(kBlack); }

    p2->cd();
    TH1 *r = (TH1 *)ha->Clone(Form("r_%s", tag[i]));
    r->Divide(hb);
    r->GetXaxis()->SetRange(1, lastclean);
    // same underlying events in both inputs -> fully correlated, no meaningful ratio error
    for (int bb = 0; bb <= r->GetNbinsX() + 1; ++bb) r->SetBinError(bb, 0);
    r->SetTitle(";E_{T}^{#gamma} [GeV];PPG19 / PPG12");
    r->SetMinimum(i == 1 ? 0.0 : 0.75); r->SetMaximum(i == 1 ? 4.5 : 1.25);
    r->GetXaxis()->SetLabelSize(0.10); r->GetXaxis()->SetTitleSize(0.12); r->GetXaxis()->SetTitleOffset(1.1);
    r->GetYaxis()->SetLabelSize(0.09); r->GetYaxis()->SetTitleSize(0.11); r->GetYaxis()->SetTitleOffset(0.55); r->GetYaxis()->SetNdivisions(505);
    r->Draw("p");
    TLine *l1 = new TLine(r->GetXaxis()->GetXmin(), 1, r->GetXaxis()->GetXmax(), 1); l1->SetLineStyle(2); l1->Draw();
    c->SaveAs(Form("%spp_et_stage_%s.pdf", FIG, tag[i]));
  }
}
