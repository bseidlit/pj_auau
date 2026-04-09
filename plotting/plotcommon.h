#include "BlairUtils.C"
#include "sPhenixStyle.C"


const int NptBins = 2;
const float ptRanges[NptBins + 1] = {15, 20, 30};

const float pTmin = 10;
const float pTmax = 26;

TH1F *frame_et_rec;
TH1F *frame_isoET;
TH1F *frame_et_truth;
TH1F *frame_iteration;
TH2F *frame_response;
TGraph *lineone;
TGraph *linezero;

string strleg1 = "#bf{#it{sPHENIX}} Internal";
string strleg2 = "AuAu #kern[-0.1]{#sqrt{#it{s}}_{NN} = 200 GeV}";
string strleg2_pp = "#it{p}+#it{p} #kern[-0.05]{#sqrt{#it{s}} = 200 GeV, 16.6 pb^{-1}}";
string strleg3 = "|#it{#eta^{#gamma}}| < 0.7";
string strlegphotonjet = "|#it{#eta^{#gamma}}| < 0.7, |#it{#eta^{jet}}| < 0.6";
string strleg4 = "#it{E}_{T}^{iso, #kern[-0.2]{#it{R}=0.3}}< 4 GeV";
string strleg5 = "         = 16.6 pb^{-1}";
string strdijet = "#Delta#it{#phi} > 3#pi/4";
string strSigMC = "PYTHIA Signal";
string strMC = "PYTHIA8 #gamma + data";
string strIncMC = "PYTHIA Inclusive Jet";

TGraph* LINE1 = new TGraph();
LINE1->SetPoint(0,-100000,1);
LINE1->SetPoint(1,100000,1);
LINE1->SetLineStyle(7);

void init_plot()
{
  SetsPhenixStyle();

  gStyle->SetHatchesLineWidth(4);

  const float et_low = 7;
  const float et_high = 50;
  frame_et_rec = new TH1F("frame_et", "", 43, et_low, et_high);
  frame_et_rec->SetXTitle("#it{E}_{T}^{#gamma,rec} [GeV]");
  frame_et_rec->GetXaxis()->SetRangeUser(8, 40);
  frame_et_rec->SetYTitle("Counts");
  frame_et_rec->GetYaxis()->SetRangeUser(5, 5e5);

  frame_et_truth = new TH1F("frame_et_truth", "", 430, et_low, et_high);
  frame_et_truth->SetXTitle("#it{E}_{T}^{#gamma, truth} [GeV]");
  frame_et_truth->GetXaxis()->SetRangeUser(8, 40);
  frame_et_truth->SetYTitle("Efficiency");
  frame_et_truth->GetYaxis()->SetRangeUser(0.2, 1.1);

  const float isoet_low = -5;
  const float isoet_high = 15;
  frame_isoET = new TH1F("frame_et", "", 43, isoet_low, isoet_high);
  frame_isoET->SetXTitle("#it{E}_{T}^{iso} [GeV]");
  frame_isoET->GetXaxis()->SetRangeUser(-3,15);
  frame_isoET->SetYTitle("scaled counts");
  frame_isoET->GetYaxis()->SetTitleOffset(1.5);
  // frame_isoET->GetYaxis()->SetRangeUser(5,5e5);

  frame_iteration = new TH1F("frame_iteration", "", 10, 0, 10);
  frame_iteration->SetXTitle("Iteration");
  frame_iteration->SetYTitle("Relative changes");
  frame_iteration->GetYaxis()->SetRangeUser(0, 0.5);
  frame_iteration->GetXaxis()->SetRangeUser(0, 10);

  frame_response = new TH2F("frame_response", "", 43, et_low, et_high, 43, et_low, et_high);
  frame_response->SetXTitle("#it{E}_{T}^{#gamma,rec} [GeV]");
  frame_response->SetYTitle("#it{E}_{T}^{#gamma,truth} [GeV]");
  frame_response->GetXaxis()->SetRangeUser(7, 40);
  frame_response->GetYaxis()->SetRangeUser(7, 40);

  lineone = new TGraph();
  lineone->SetPoint(0, -1e6, 1);
  lineone->SetPoint(1, 1e6, 1);
  lineone->SetLineColor(kBlack);
  lineone->SetLineStyle(7);

  linezero = new TGraph();
  linezero->SetPoint(0, -1e6, 0);
  linezero->SetPoint(1, 1e6, 0);
  linezero->SetLineColor(kBlack);
  linezero->SetLineStyle(7);
}

std::pair<TH1F *, TH1F *> calcDelta(TH1F *h1, TH1F *h2, std::string name)
{
  TH1F *hDelta = (TH1F *)h1->Clone(name.c_str());
  hDelta->Add(h2, -1);
  TH1F *hDeltaRel = (TH1F *)hDelta->Clone((name + "_rel").c_str());
  hDeltaRel->Divide(h2);
  return {hDelta, hDeltaRel};
}
