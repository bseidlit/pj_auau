#include "PhotonJetPlotInput.h"
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
void plot_pp_et_levels()
{
  gStyle->SetOptStat(0);
  const char* R="results/";
  TFile* fa=PJPlot::Open(Form("%sdata_histo_pj_pp_nom.root",R));
  TFile* fb=PJPlot::Open(Form("%sdata_histo_bdt_pjlike.root",R));
  const char* hn[6]={"h_all_cluster_0","h_tight_cluster_0","h_tight_iso_cluster_0","h_tight_noniso_cluster_0","h_nontight_iso_cluster_0","h_nontight_noniso_cluster_0"};
  const char* lb[6]={"all clusters (E_{T} #geq 15)","tight","A: tight + isolated","B: tight + non-isolated","C: non-tight + isolated","D: non-tight + non-isolated"};
  const int col[6]={kGray+2,kBlack,kRed+1,kOrange+7,kBlue+1,kGreen+2};
  TCanvas* c=new TCanvas("c","",750,650); gPad->SetLogy(); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12);
  TLegend* leg=new TLegend(0.45,0.55,0.93,0.88); leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.03);
  for(int i=0;i<6;i++){
    TH1* ha=(TH1*)fa->Get(hn[i]); TH1* hb=(TH1*)fb->Get(hn[i]);
    ha->GetXaxis()->SetRange(1, ha->GetXaxis()->FindBin(31.9)); hb->GetXaxis()->SetRange(1, hb->GetXaxis()->FindBin(31.9)); // base_E boundary at 35 touches the 32-36 bin
    ha->SetLineColor(col[i]); ha->SetMarkerColor(col[i]); ha->SetMarkerStyle(20); ha->SetLineWidth(2);
    hb->SetLineColor(col[i]); hb->SetMarkerColor(col[i]); hb->SetMarkerStyle(24); hb->SetLineWidth(1);
    ha->SetTitle(";E_{T}^{#gamma} [GeV];photon candidates / bin (weighted by 1/#varepsilon_{L1})");
    ha->SetMinimum(1); ha->SetMaximum(3e5); ha->GetYaxis()->SetTitleOffset(1.4);
    ha->Draw(i==0?"p e":"p e same"); hb->Draw("p e same");
    leg->AddEntry(ha,lb[i],"pl");
  }
  leg->AddEntry((TObject*)0,"filled: new tree, open: PPG12","");
  leg->Draw();
  TLatex t; t.SetNDC(); t.SetTextSize(0.035); t.DrawLatex(0.15,0.92,"#it{p}+#it{p} #sqrt{s} = 200 GeV data, |#eta^{#gamma}| < 0.7, |z_{vtx}| < 30 cm");
  c->SaveAs("figures/pp_data_et_spectrum_levels.pdf");
}
