// compare_pj_vs_ppg12.C -- side-by-side of a PhotonJetTrees maker result with a PPG12 slimtree
// result (raw ABCD counts, signal leakage fractions, TEfficiencies, MBD eff, purity, final result).
//   root -l -b -q 'compare_pj_vs_ppg12.C("pj_pp_nom","bdt_pjcert")'
#include "PhotonJetPlotInput.h"
#include <TH1.h>
#include <TEfficiency.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <cstdio>
void compare_pj_vs_ppg12(const char* pj="pj_pp_nom", const char* nom="bdt_nom")
{
  const char* R="results/";
  TFile *fd1=PJPlot::Open(Form("%sdata_histo_%s.root",R,pj)), *fd0=PJPlot::Open(Form("%sdata_histo_%s.root",R,nom));
  TFile *fm1=PJPlot::Open(Form("%sMC_efficiency_%s.root",R,pj)), *fm0=PJPlot::Open(Form("%sMC_efficiency_%s.root",R,nom));
  TFile *ff1=PJPlot::Open(Form("%sPhoton_final_%s.root",R,pj)), *ff0=PJPlot::Open(Form("%sPhoton_final_%s.root",R,nom));
  const char* reg[4]={"tight_iso","tight_noniso","nontight_iso","nontight_noniso"};
  printf("\n== raw data ABCD counts (weighted by 1/eps_L1), %s vs %s ==\n", pj, nom);
  printf("%-10s","bin"); for(int r=0;r<4;r++) printf(" | %14s %14s",(TString("pj_")+TString(reg[r]).Remove(5)).Data(),"nom"); printf("\n");
  TH1* hA1=(TH1*)fd1->Get("h_tight_iso_cluster_0"); TH1* hA0=(TH1*)fd0->Get("h_tight_iso_cluster_0");
  for(int b=1;b<=hA1->GetNbinsX();b++){
    double lo=hA1->GetBinLowEdge(b), hi=hA1->GetBinLowEdge(b+1);
    int b0=hA0->FindBin(0.5*(lo+hi));
    printf("%4.0f-%-5.0f",lo,hi);
    for(int r=0;r<4;r++){ TH1* h1=(TH1*)fd1->Get(Form("h_%s_cluster_0",reg[r])); TH1* h0=(TH1*)fd0->Get(Form("h_%s_cluster_0",reg[r]));
      printf(" | %14.1f %14.1f", h1->GetBinContent(b), h0 ? h0->GetBinContent(b0):-1); }
    printf("\n");
  }
  printf("\n== MC signal leakage fractions cB=B/A cC=C/A cD=D/A (signal MC) ==\n");
  for(int b=1;b<=hA1->GetNbinsX();b++){
    double lo=hA1->GetBinLowEdge(b), hi=hA1->GetBinLowEdge(b+1); int b0=hA0->FindBin(0.5*(lo+hi));
    TH1 *a1=(TH1*)fm1->Get("h_tight_iso_cluster_signal_0"), *a0=(TH1*)fm0->Get("h_tight_iso_cluster_signal_0");
    printf("%4.0f-%-5.0f",lo,hi);
    for(int r=1;r<4;r++){ TH1* h1=(TH1*)fm1->Get(Form("h_%s_cluster_signal_0",reg[r])); TH1* h0=(TH1*)fm0->Get(Form("h_%s_cluster_signal_0",reg[r]));
      printf(" | c%c pj=%.4f nom=%.4f", "ABCD"[r], h1->GetBinContent(b)/a1->GetBinContent(b), h0->GetBinContent(b0)/a0->GetBinContent(b0)); }
    printf("\n");
  }
  printf("\n== efficiencies vs truth pT (TEfficiency, %s vs %s) ==\n", pj, nom);
  const char* en[3]={"eff_reco_eta_0","eff_iso_eta_0","eff_id_eta_0"};
  TEfficiency* e1=(TEfficiency*)fm1->Get(en[0]); const TH1* tot=e1->GetTotalHistogram();
  printf("%-10s | %-23s | %-23s | %-23s | %-15s\n","truth bin","reco pj/nom","iso pj/nom","id pj/nom","MBDeff pj/nom");
  TH1 *v1=(TH1*)fm1->Get("h_truth_pT_vertexcut_0"), *m1=(TH1*)fm1->Get("h_truth_pT_vertexcut_mbd_cut_0");
  TH1 *v0=(TH1*)fm0->Get("h_truth_pT_vertexcut_0"), *m0=(TH1*)fm0->Get("h_truth_pT_vertexcut_mbd_cut_0");
  for(int b=1;b<=tot->GetNbinsX();b++){
    printf("%4.0f-%-5.0f", tot->GetBinLowEdge(b), tot->GetBinLowEdge(b+1));
    for(int k=0;k<3;k++){ TEfficiency* a=(TEfficiency*)fm1->Get(en[k]); TEfficiency* c=(TEfficiency*)fm0->Get(en[k]);
      printf(" | %.4f / %.4f", a->GetEfficiency(b), c->GetEfficiency(b)); }
    printf(" | %.3f / %.3f\n", m1->GetBinContent(b)/v1->GetBinContent(b), m0->GetBinContent(b)/v0->GetBinContent(b));
  }
  printf("\n== purity (gpurity_leak) and final result, %s vs %s ==\n", pj, nom);
  TGraphErrors *p1=(TGraphErrors*)ff1->Get("gpurity_leak"), *p0=(TGraphErrors*)ff0->Get("gpurity_leak");
  TH1 *r1=(TH1*)ff1->Get("h_unfold_sub_result"), *r0=(TH1*)ff0->Get("h_unfold_sub_result");
  for(int i=0;i<p1->GetN();i++){ double x,y,x0,y0; p1->GetPoint(i,x,y); int j=-1; for(int k=0;k<p0->GetN();k++){p0->GetPoint(k,x0,y0); if(fabs(x0-x)<0.5){j=k;break;}}
    double y0v=-1; if(j>=0) p0->GetPoint(j,x0,y0v);
    int br=r1?r1->FindBin(x):-1, br0=r0?r0->FindBin(x):-1;
    printf("pT=%5.1f  purity pj=%.3f +- %.3f  nom=%.3f   xsec pj=%.4g nom=%.4g ratio=%.3f \n", x, y, p1->GetErrorY(i), y0v,
      r1?r1->GetBinContent(br):0, r0?r0->GetBinContent(br0):0, (r0&&r1&&r0->GetBinContent(br0)>0)? r1->GetBinContent(br)/r0->GetBinContent(br0):0);
  }
}
