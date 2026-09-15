// plot_pp_purity_corrected.C -- purity-corrected (ABCD background-subtracted, signal-leakage
// corrected) pp data ET spectrum, maker vs PPG12, with ratio. Uses the UNSCALED
// h_data_sub_leak from Photon_final_*.root (raw corrected counts, no lumi/binwidth scaling).
// The data entering both sides is identical; the purity differs only through the MC
// leakage templates, so the ratio directly shows the sim-sample effect on a data-level
// quantity. Ratio drawn without error bars (fully correlated data component).
//   root -l -b -q 'plot_pp_purity_corrected.C("pj_pp_cert","bdt_pjcert")'
#include "PhotonJetPlotInput.h"
#include <TH1.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <cstdio>

void plot_pp_purity_corrected(const char *a = "pj_pp_cert", const char *b = "bdt_pjcert")
{
  gStyle->SetOptStat(0);
  const char *R = "results/";
  const char *FIG = "figures/";
  TFile *fa = PJPlot::Open(Form("%sPhoton_final_%s.root", R, a));
  TFile *fb = PJPlot::Open(Form("%sPhoton_final_%s.root", R, b));
  TH1 *ha = (TH1 *)fa->Get("h_data_sub_leak"); TH1 *hb = (TH1 *)fb->Get("h_data_sub_leak");
  TGraphErrors *pa = (TGraphErrors *)fa->Get("gpurity_leak"); TGraphErrors *pb = (TGraphErrors *)fb->Get("gpurity_leak");
  ha = (TH1 *)ha->Clone("pa_h"); hb = (TH1 *)hb->Clone("pb_h");
  const int lastclean = ha->GetXaxis()->FindBin(31.9);
  ha->GetXaxis()->SetRange(1, lastclean); hb->GetXaxis()->SetRange(1, lastclean);

  printf("== purity-corrected data spectrum (A x purity, leakage-corrected counts): %s vs %s ==\n", a, b);
  printf("%-9s | %10s %10s %7s | purity %s / %s\n", "bin", a, b, "ratio", a, b);
  for (int i = 1; i <= lastclean; ++i)
  {
    double xa, ya, xb, yb; pa->GetPoint(i - 1, xa, ya); pb->GetPoint(i - 1, xb, yb);
    printf("%3.0f-%-5.0f | %10.1f %10.1f %7.3f | %.3f / %.3f\n", ha->GetBinLowEdge(i), ha->GetBinLowEdge(i + 1),
           ha->GetBinContent(i), hb->GetBinContent(i), hb->GetBinContent(i) > 0 ? ha->GetBinContent(i) / hb->GetBinContent(i) : 0, ya, yb);
  }

  TCanvas *c = new TCanvas("c_pc", "", 700, 750);
  TPad *p1 = new TPad("p1", "", 0, 0.32, 1, 1); TPad *p2 = new TPad("p2", "", 0, 0, 1, 0.32);
  p1->SetBottomMargin(0.025); p1->SetLeftMargin(0.14);
  p2->SetTopMargin(0.04); p2->SetBottomMargin(0.32); p2->SetLeftMargin(0.14); p2->SetGridy();
  p1->Draw(); p2->Draw();
  p1->cd(); p1->SetLogy();
  hb->SetLineColor(kBlack); hb->SetMarkerColor(kBlack); hb->SetMarkerStyle(24); hb->SetMarkerSize(1.2); hb->SetLineWidth(2);
  ha->SetLineColor(kRed + 1); ha->SetMarkerColor(kRed + 1); ha->SetMarkerStyle(20); ha->SetMarkerSize(1.1);
  hb->SetTitle(";;(A #times purity) / (#DeltaE_{T} #upoint L) [pb/GeV]");
  hb->GetXaxis()->SetLabelSize(0); hb->GetYaxis()->SetTitleOffset(1.35);
  hb->SetMaximum(TMath::Max(ha->GetMaximum(), hb->GetMaximum()) * 8); hb->SetMinimum(0.05);
  hb->Draw("p e"); ha->Draw("p e same");
  TLegend *leg = new TLegend(0.42, 0.62, 0.93, 0.80); leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.038);
  leg->AddEntry(ha, Form("%s (maker)", a), "pl");
  leg->AddEntry(hb, Form("%s (PPG12)", b), "pl");
  leg->Draw();
  TLatex t; t.SetNDC();
  t.SetTextSize(0.045); t.SetTextFont(62); t.DrawLatex(0.17, 0.86, "A #times purity (leakage corrected)");
  t.SetTextSize(0.034); t.SetTextFont(42);
  t.DrawLatex(0.17, 0.81, "#it{p}+#it{p} #sqrt{s} = 200 GeV data, |#eta^{#gamma}| < 0.7, |z_{vtx}| < 30 cm");
  t.DrawLatex(0.17, 0.765, "identical data, purity differs only via the MC leakage templates");
  p2->cd();
  TH1 *r = (TH1 *)ha->Clone("r_pc"); r->Divide(hb);
  for (int bb = 0; bb <= r->GetNbinsX() + 1; ++bb) r->SetBinError(bb, 0);
  r->GetXaxis()->SetRange(1, lastclean);
  r->SetTitle(";E_{T}^{#gamma} [GeV];maker / PPG12");
  r->SetMinimum(0.85); r->SetMaximum(1.15);
  r->GetXaxis()->SetLabelSize(0.10); r->GetXaxis()->SetTitleSize(0.12); r->GetXaxis()->SetTitleOffset(1.1);
  r->GetYaxis()->SetLabelSize(0.09); r->GetYaxis()->SetTitleSize(0.11); r->GetYaxis()->SetTitleOffset(0.55); r->GetYaxis()->SetNdivisions(505);
  r->Draw("p");
  TLine *l1 = new TLine(16, 1, 32, 1); l1->SetLineStyle(2); l1->Draw();
  c->SaveAs(Form("%spp_data_puritycorr_%s_vs_%s.pdf", FIG, a, b));
}
