// plot_pp_final_corrected.C -- fully corrected pp spectrum (ABCD purity + Bayesian
// unfolding + eps_reco*eps_iso*eps_id*eps_MBDvtx + luminosity), PPG19 (tree maker) vs PPG12, with
// ratio and a per-bin decomposition printed (purity-corrected ratio x 1/efficiency
// ratio ~ final ratio, the rest is unfolding). Ratio without error bars (same data).
//   root -l -b -q 'plot_pp_final_corrected.C("pj_pp_cert","bdt_pjcert")'
#include "PhotonJetPlotInput.h"
#include <TH1.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TEfficiency.h>
#include <TStyle.h>
#include <cstdio>

void plot_pp_final_corrected(const char *a = "pj_pp_cert", const char *b = "bdt_pjcert")
{
  gStyle->SetOptStat(0);
  const char *R = "results/";
  const char *FIG = "figures/";
  TFile *fa = PJPlot::Open(Form("%sPhoton_final_%s.root", R, a));
  TFile *fb = PJPlot::Open(Form("%sPhoton_final_%s.root", R, b));
  TFile *ma = PJPlot::Open(Form("%sMC_efficiency_%s.root", R, a));
  TFile *mb = PJPlot::Open(Form("%sMC_efficiency_%s.root", R, b));
  TH1 *ha = (TH1 *)((TH1 *)fa->Get("h_unfold_sub_result"))->Clone("fa_h");
  TH1 *hb = (TH1 *)((TH1 *)fb->Get("h_unfold_sub_result"))->Clone("fb_h");
  TH1 *pa = (TH1 *)fa->Get("h_data_sub_leak"); TH1 *pb = (TH1 *)fb->Get("h_data_sub_leak");
  // h_unfold_sub_result is in the TRUTH binning (8-45); show only the reco-supported 16-32 range
  const int firstb = ha->GetXaxis()->FindBin(16.1), lastclean = ha->GetXaxis()->FindBin(31.9);
  ha->GetXaxis()->SetRange(firstb, lastclean); hb->GetXaxis()->SetRange(firstb, lastclean);

  auto effprod = [](TFile *m, TFile *f, double et) {
    TEfficiency *er = (TEfficiency *)m->Get("eff_reco_eta_0"), *ei = (TEfficiency *)m->Get("eff_iso_eta_0"), *ed = (TEfficiency *)m->Get("eff_id_eta_0");
    TH1 *v = (TH1 *)m->Get("h_truth_pT_vertexcut_0"), *mm = (TH1 *)m->Get("h_truth_pT_vertexcut_mbd_cut_0");
    const int bt = er->GetTotalHistogram()->GetXaxis()->FindBin(et);
    return er->GetEfficiency(bt) * ei->GetEfficiency(bt) * ed->GetEfficiency(bt) * (mm->GetBinContent(bt) / v->GetBinContent(bt));
  };
  printf("== fully corrected spectrum: %s / %s ==\n", a, b);
  printf("%-9s | final ratio | purity-corr ratio | eff-product ratio | residual (unfolding+lumi)\n", "bin");
  for (int i = firstb; i <= lastclean; ++i)
  {
    const double et = ha->GetBinCenter(i);
    const int ip = pa->GetXaxis()->FindBin(et); // reco binning differs from the truth binning
    const double rf = hb->GetBinContent(i) > 0 ? ha->GetBinContent(i) / hb->GetBinContent(i) : 0;
    const double rp = pb->GetBinContent(ip) > 0 ? pa->GetBinContent(ip) / pb->GetBinContent(ip) : 0;
    const double re = effprod(ma, fa, et) / effprod(mb, fb, et);
    printf("%3.0f-%-5.0f | %11.3f | %17.3f | %17.3f | %9.3f\n", ha->GetBinLowEdge(i), ha->GetBinLowEdge(i + 1), rf, rp, re, (rp > 0 && re > 0) ? rf / (rp / re) : 0);
  }

  TCanvas *c = new TCanvas("c_fc", "", 700, 750);
  TPad *p1 = new TPad("p1", "", 0, 0.32, 1, 1); TPad *p2 = new TPad("p2", "", 0, 0, 1, 0.32);
  p1->SetBottomMargin(0.025); p1->SetLeftMargin(0.14);
  p2->SetTopMargin(0.04); p2->SetBottomMargin(0.32); p2->SetLeftMargin(0.14); p2->SetGridy();
  p1->Draw(); p2->Draw();
  p1->cd(); p1->SetLogy();
  hb->SetLineColor(kBlack); hb->SetMarkerColor(kBlack); hb->SetMarkerStyle(24); hb->SetMarkerSize(1.2); hb->SetLineWidth(2);
  ha->SetLineColor(kRed + 1); ha->SetMarkerColor(kRed + 1); ha->SetMarkerStyle(20); ha->SetMarkerSize(1.1);
  hb->SetTitle(";;d^{2}#sigma/dE_{T}d#eta [pb/GeV]");
  hb->GetXaxis()->SetLabelSize(0); hb->GetYaxis()->SetTitleOffset(1.35);
  hb->SetMaximum(TMath::Max(ha->GetMaximum(), hb->GetMaximum()) * 8); hb->SetMinimum(0.3);
  hb->Draw("p e"); ha->Draw("p e same");
  TLegend *leg = new TLegend(0.42, 0.62, 0.93, 0.80); leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.038);
  leg->AddEntry(ha, "PPG19 tree", "pl");
  leg->AddEntry(hb, "PPG12 tree", "pl");
  leg->Draw();
  TLatex t; t.SetNDC();
  t.SetTextSize(0.045); t.SetTextFont(62); t.DrawLatex(0.17, 0.86, "purity + unfolding + efficiency corrected");
  t.SetTextSize(0.034); t.SetTextFont(42);
  t.DrawLatex(0.17, 0.81, "#it{p}+#it{p} #sqrt{s} = 200 GeV, |#eta^{#gamma}| < 0.7, |z_{vtx}| < 30 cm, isolated photons");
  p2->cd();
  TH1 *r = (TH1 *)ha->Clone("r_fc"); r->Divide(hb);
  for (int bb = 0; bb <= r->GetNbinsX() + 1; ++bb) r->SetBinError(bb, 0);
  r->GetXaxis()->SetRange(firstb, lastclean);
  r->SetTitle(";E_{T}^{#gamma} [GeV];PPG19 / PPG12");
  r->SetMinimum(0.8); r->SetMaximum(1.2);
  r->GetXaxis()->SetLabelSize(0.10); r->GetXaxis()->SetTitleSize(0.12); r->GetXaxis()->SetTitleOffset(1.1);
  r->GetYaxis()->SetLabelSize(0.09); r->GetYaxis()->SetTitleSize(0.11); r->GetYaxis()->SetTitleOffset(0.55); r->GetYaxis()->SetNdivisions(505);
  r->Draw("p");
  TLine *l1 = new TLine(16, 1, 32, 1); l1->SetLineStyle(2); l1->Draw();
  c->SaveAs(Form("%spp_final_corrected_%s_vs_%s.pdf", FIG, a, b));
}
