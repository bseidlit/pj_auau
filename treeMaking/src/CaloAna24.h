// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef CALOANA24_H
#define CALOANA24_H

#include <fun4all/SubsysReco.h>

#include <phool/onnxlib.h>

#include <TFile.h>
#include <TH3.h>
#include <TH2.h>
#include <TH1.h>
#include <TProfile2D.h>
#include <TTree.h>
#include <algorithm>
#include <string>
#include <TLorentzVector.h>

class PHCompositeNode;
class CaloEvalStack;
class CaloRawClusterEval;
class CaloTruthEval;
class RawTowerGeom;
class RawTowerGeomContainer;
class RawClusterContainer;
class TowerInfoContainer;
class PHG4TruthInfoContainer;
class CaloEvalStack;

namespace HepMC
{
  class GenEvent;
}

class CaloAna24 : public SubsysReco
{
public:
  CaloAna24(const std::string &name = "CaloAna24", const std::string &filename = "caloana.root");

  ~CaloAna24() override;

  int Init(PHCompositeNode *topNode) override;

  int InitRun(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  /// Reset
  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;

  void set_isMC(bool isMC_) { isMC = isMC_; }

  void set_isSingleParticle(bool isSingleParticle_)
  {
    isSingleParticle = isSingleParticle_;
    if (isSingleParticle)
      isMC = true;
  }

  void set_using_trigger_bits(std::vector<int> trigger_bits)
  {
    using_trigger_bits = trigger_bits;
  }

  void set_clusterpTmin(float pTmin) { clusterpTmin = pTmin; }

  void set_particlepTmin(float pTmin) { particlepTmin = pTmin; }

  void set_jet_node_name(const std::string &jet_node_name);
  void set_jet_node_names(const std::vector<std::string> &jet_node_names);
  void set_truth_jet_node_name(const std::string &truth_jet_node_name);
  void set_truth_jet_node_names(const std::vector<std::string> &truth_jet_node_names);

private:
  int ievent = 0;
  Ort::Session *onnxmodule{nullptr};
  std::string m_modelPath{"/sphenix/u/shuhang98/core_patch/coresoftware/offline/packages/CaloReco/functional_model_single.onnx"};
  std::string m_outputFileName{"caloana.root"};
  TFile *fout;

  TTree *slimtree;

  int runnumber{0};

  bool isMC{true};
  bool isSingleParticle{false};
  int m_scaledtrigger[64] = {0};
  bool initilized = false;
  long long initscaler[64][3] = {0};
  long long currentscaler[64][3] = {0};
  long long currentscaler_raw[64] = {0};
  long long currentscaler_live[64] = {0};
  long long currentscaler_scaled[64] = {0};
  bool scaledtrigger[64] = {false};
  bool livetrigger[64] = {false};
  int nscaledtrigger[64] = {0};
  int nlivetrigger[64] = {0};
  float trigger_prescale[64] = {-1};
  int m_eventnumber{0};
  float m_totalEMCal_energy{0};
  float m_totalIHCal_energy{0};
  float m_totalOHCal_energy{0};
  int nevent{0};

  float m_cent{-9999};
  float vertexz{-9999};
  int mbdnorthhit{0};
  int mbdsouthhit{0};
  float mbd_time{-9999};
  float mbd_north_time{-9999};
  float mbd_south_time{-9999};
  float mbdnorthq[64] = {0};
  float mbdsouthq[64] = {0};
  float mbdnortht[64] = {0};
  float mbdsoutht[64] = {0};
  float mbdnorthqsum{0};
  float mbdsouthqsum{0};
  float mbdnorthtmean{0};
  float mbdsouthtmean{0};
  float _Psi2{0};
  float vertexz_truth{-9999};
  int m_pythiaid{-9999};
  float m_energyscale{-1.0};
  float particlepTmin{1};
  static const int nparticlesmax = 10000;
  int nparticles{0};
  float particle_E[nparticlesmax] = {0};
  float particle_Pt[nparticlesmax] = {0};
  float particle_Eta[nparticlesmax] = {0};
  float particle_Phi[nparticlesmax] = {0};
  int particle_pid[nparticlesmax] = {0};
  int particle_trkid[nparticlesmax] = {0};
  int particle_photonclass[nparticlesmax] = {0};
  int particle_photon_mother_pid[nparticlesmax] = {0};
  float particle_truth_iso_02[nparticlesmax] = {0};
  float particle_truth_iso_03[nparticlesmax] = {0};
  float particle_truth_iso_04[nparticlesmax] = {0};

  int particle_converted[nparticlesmax] = {0};

  static const int ndaughtermax = 100;
  int ndaughter{0};
  int daughter_pid[ndaughtermax] = {0};
  int daughter_parent_trackid[ndaughtermax] = {0};
  float daughter_E[ndaughtermax] = {0};
  float daughter_Pt[ndaughtermax] = {0};
  float daughter_Eta[ndaughtermax] = {0};
  float daughter_Phi[ndaughtermax] = {0};
  float daughter_vtx_x[ndaughtermax] = {0};
  float daughter_vtx_y[ndaughtermax] = {0};
  float daughter_vtx_z[ndaughtermax] = {0};

  std::vector<std::string> clusternamelist = {"CLUSTERINFO_CEMC_NO_SPLIT", "CLUSTERINFO_CEMC"};
  static const int nclustercontainer = 2;
  // cluster wise stuff
  float clusterpTmin{12};
  static const int nclustermax = 10000;
  int ncluster[nclustercontainer] = {0};
  float cluster_E[nclustercontainer][nclustermax] = {0};
  float cluster_ecore[nclustercontainer][nclustermax] = {0};
  float cluster_Et[nclustercontainer][nclustermax] = {0};
  float cluster_Eta[nclustercontainer][nclustermax] = {0};
  float cluster_Phi[nclustercontainer][nclustermax] = {0};
  float cluster_prob[nclustercontainer][nclustermax] = {0};
  float cluster_merged_prob[nclustercontainer][nclustermax] = {0};
  float cluster_CNN_prob[nclustercontainer][nclustermax] = {0};
  int cluster_truthtrkID[nclustercontainer][nclustermax] = {0};
  int cluster_pid[nclustercontainer][nclustermax] = {0};
  float cluster_iso_02[nclustercontainer][nclustermax] = {0};
  float cluster_iso_03[nclustercontainer][nclustermax] = {0};
  float cluster_iso_04[nclustercontainer][nclustermax] = {0};
  float cluster_iso_03_emcal[nclustercontainer][nclustermax] = {0};
  float cluster_iso_03_hcalin[nclustercontainer][nclustermax] = {0};
  float cluster_iso_03_hcalout[nclustercontainer][nclustermax] = {0};
  float cluster_iso_03_60_emcal[nclustercontainer][nclustermax] = {0};
  float cluster_iso_03_60_hcalin[nclustercontainer][nclustermax] = {0};
  float cluster_iso_03_60_hcalout[nclustercontainer][nclustermax] = {0};
  float cluster_iso_03_70_emcal[nclustercontainer][nclustermax] = {0};
  float cluster_iso_03_70_hcalin[nclustercontainer][nclustermax] = {0};
  float cluster_iso_03_70_hcalout[nclustercontainer][nclustermax] = {0};
  float cluster_iso_005_70_emcal[nclustercontainer][nclustermax] = {0};
  float cluster_iso_01_70_emcal[nclustercontainer][nclustermax] = {0};
  float cluster_iso_02_70_emcal[nclustercontainer][nclustermax] = {0};
  float cluster_iso_03_120_emcal[nclustercontainer][nclustermax] = {0};
  float cluster_iso_03_120_hcalin[nclustercontainer][nclustermax] = {0};
  float cluster_iso_03_120_hcalout[nclustercontainer][nclustermax] = {0};
  float cluster_iso_04_emcal[nclustercontainer][nclustermax] = {0};
  float cluster_iso_04_hcalin[nclustercontainer][nclustermax] = {0};
  float cluster_iso_04_hcalout[nclustercontainer][nclustermax] = {0};
  float cluster_iso_03_sub1_emcal[nclustercontainer][nclustermax] = {0};
  float cluster_iso_03_sub1_hcalin[nclustercontainer][nclustermax] = {0};
  float cluster_iso_03_sub1_hcalout[nclustercontainer][nclustermax] = {0};
  float cluster_iso_04_sub1_emcal[nclustercontainer][nclustermax] = {0};
  float cluster_iso_04_sub1_hcalin[nclustercontainer][nclustermax] = {0};
  float cluster_iso_04_sub1_hcalout[nclustercontainer][nclustermax] = {0};

  float cluster_iso_topo_03[nclustercontainer][nclustermax] = {0};
  float cluster_iso_topo_04[nclustercontainer][nclustermax] = {0};
  float cluster_iso_topo_soft_03[nclustercontainer][nclustermax] = {0};
  float cluster_iso_topo_soft_04[nclustercontainer][nclustermax] = {0};

  static const int arrayntower = 49;
  // shower shapes
  float cluster_e_array[nclustercontainer][nclustermax][arrayntower] = {0};
  float cluster_adc_array[nclustercontainer][nclustermax][arrayntower] = {0};
  float cluster_time_array[nclustercontainer][nclustermax][arrayntower] = {0};
  int cluster_e_array_idx[nclustercontainer][nclustermax][arrayntower] = {0};
  int cluster_status_array[nclustercontainer][nclustermax][arrayntower] = {0};
  int cluster_ownership_array[nclustercontainer][nclustermax][arrayntower] = {0};
  float cluster_e1[nclustercontainer][nclustermax] = {0};
  float cluster_e2[nclustercontainer][nclustermax] = {0};
  float cluster_e3[nclustercontainer][nclustermax] = {0};
  float cluster_e4[nclustercontainer][nclustermax] = {0};
  float cluster_et1[nclustercontainer][nclustermax] = {0};
  float cluster_et2[nclustercontainer][nclustermax] = {0};
  float cluster_et3[nclustercontainer][nclustermax] = {0};
  float cluster_et4[nclustercontainer][nclustermax] = {0};
  float cluster_ietacent[nclustercontainer][nclustermax] = {0};
  float cluster_iphicent[nclustercontainer][nclustermax] = {0};
  float cluster_weta[nclustercontainer][nclustermax] = {0};
  float cluster_wphi[nclustercontainer][nclustermax] = {0};
  float cluster_weta_cog[nclustercontainer][nclustermax] = {0};
  float cluster_wphi_cog[nclustercontainer][nclustermax] = {0};
  float cluster_weta_cogx[nclustercontainer][nclustermax] = {0};
  float cluster_wphi_cogx[nclustercontainer][nclustermax] = {0};
  int cluster_detamax[nclustercontainer][nclustermax] = {0};
  int cluster_dphimax[nclustercontainer][nclustermax] = {0};
  int cluster_nsaturated[nclustercontainer][nclustermax] = {0};

  float cluster_e11[nclustercontainer][nclustermax] = {0};
  float cluster_e22[nclustercontainer][nclustermax] = {0};
  float cluster_e13[nclustercontainer][nclustermax] = {0};
  float cluster_e15[nclustercontainer][nclustermax] = {0};
  float cluster_e17[nclustercontainer][nclustermax] = {0};
  float cluster_e31[nclustercontainer][nclustermax] = {0};
  float cluster_e51[nclustercontainer][nclustermax] = {0};
  float cluster_e71[nclustercontainer][nclustermax] = {0};
  float cluster_e33[nclustercontainer][nclustermax] = {0};
  float cluster_e35[nclustercontainer][nclustermax] = {0};
  float cluster_e37[nclustercontainer][nclustermax] = {0};
  float cluster_e53[nclustercontainer][nclustermax] = {0};
  float cluster_e73[nclustercontainer][nclustermax] = {0};
  float cluster_e55[nclustercontainer][nclustermax] = {0};
  float cluster_e57[nclustercontainer][nclustermax] = {0};
  float cluster_e75[nclustercontainer][nclustermax] = {0};
  float cluster_e77[nclustercontainer][nclustermax] = {0};
  float cluster_w32[nclustercontainer][nclustermax] = {0};
  float cluster_e32[nclustercontainer][nclustermax] = {0};
  float cluster_w52[nclustercontainer][nclustermax] = {0};
  float cluster_e52[nclustercontainer][nclustermax] = {0};
  float cluster_w72[nclustercontainer][nclustermax] = {0};
  float cluster_e72[nclustercontainer][nclustermax] = {0};

  float cluster_ihcal_et[nclustercontainer][nclustermax] = {0};
  float cluster_ohcal_et[nclustercontainer][nclustermax] = {0};
  float cluster_ihcal_et22[nclustercontainer][nclustermax] = {0};
  float cluster_ohcal_et22[nclustercontainer][nclustermax] = {0};
  float cluster_ihcal_et33[nclustercontainer][nclustermax] = {0};
  float cluster_ohcal_et33[nclustercontainer][nclustermax] = {0};
  int cluster_ihcal_ieta[nclustercontainer][nclustermax] = {0};
  int cluster_ihcal_iphi[nclustercontainer][nclustermax] = {0};
  int cluster_ohcal_ieta[nclustercontainer][nclustermax] = {0};
  int cluster_ohcal_iphi[nclustercontainer][nclustermax] = {0};
  // Number of radii

  static const int njetcontainermax = 8;
  std::vector<std::string> jetnamelist = {"AntiKt_unsubtracted_r04"};
  int njetcontainer{1};
  static const int njetmax = 1000;
  int njet[njetcontainermax] = {0};
  float jet_E[njetcontainermax][njetmax] = {0};
  float jet_Et[njetcontainermax][njetmax] = {0};
  float jet_Pt[njetcontainermax][njetmax] = {0};
  float jet_Eta[njetcontainermax][njetmax] = {0};
  float jet_Phi[njetcontainermax][njetmax] = {0};
  float jet_time[njetcontainermax][njetmax] = {0};
  float jet_emcal_calo_E[njetcontainermax][njetmax] = {0};
  float jet_ihcal_calo_E[njetcontainermax][njetmax] = {0};
  float jet_ohcal_calo_E[njetcontainermax][njetmax] = {0};

  static const int njet_truthcontainermax = 8;
  std::vector<std::string> truthjetnamelist = {"AntiKt_Truth_r04"};
  int njet_truthcontainer{1};
  static const int njet_truthmax = 1000;
  int njet_truth[njet_truthcontainermax] = {0};
  float jet_truth_E[njet_truthcontainermax][njet_truthmax] = {0};
  float jet_truth_Et[njet_truthcontainermax][njet_truthmax] = {0};
  float jet_truth_Pt[njet_truthcontainermax][njet_truthmax] = {0};
  float jet_truth_Eta[njet_truthcontainermax][njet_truthmax] = {0};
  float jet_truth_Phi[njet_truthcontainermax][njet_truthmax] = {0};

  static const int nRadii = 3;

  TH3F *h_tracking_radiograph;

  TH1I *h_sim_cross_counting;

  const int truthisocut = 4;

  int process_cluster(std::vector<TLorentzVector> goodcluster);

  std::pair<int, int> photon_type(int barcode);

  void shift_tower_index(int &ieta, int &iphi, int maxeta, int maxphi)
  {
    if (ieta < 0)
      ieta = -1;
    if (ieta >= maxeta)
      ieta = -1;
    if (iphi < 0)
      iphi += maxphi;
    if (iphi >= maxphi)
      iphi -= maxphi;
  }

  double getTowerEta(RawTowerGeom *tower_geom, double vx, double vy, double vz);

  float calculateET(float eta, float phi, float dR, int layer, float min_E, bool use_subtracted = false); // layer: 0 EMCal, 1 IHCal, 2 OHCal

  float calculateET_topo(float eta, float phi, float dR, RawClusterContainer *clusterContainer);

  std::vector<int> find_closest_hcal_tower(float eta, float phi, RawTowerGeomContainer *rawtowergeom, TowerInfoContainer *towercontainer, float vertex_z, bool isihcal);

  inline /*const*/ float deltaR(float eta1, float eta2, float phi1, float phi2)
  {
    float deta = eta1 - eta2;
    float dphi = phi1 - phi2;
    if (dphi > M_PI)
      dphi -= 2 * M_PI; // corrects to keep range -pi to pi
    if (dphi < -1 * M_PI)
      dphi += 2 * M_PI; // corrects to keep range -pi to pi
    return sqrt(deta * deta + dphi * dphi);
  }

  float DeltaR(TLorentzVector photon1, TLorentzVector photon2)
  {
    float deta = photon1.PseudoRapidity() - photon2.PseudoRapidity();
    float dphi = abs(photon1.Phi() - photon2.Phi());

    if (dphi > M_PI)
      dphi = 2 * M_PI - dphi;
    float dr = sqrt(deta * deta + dphi * dphi);
    return dr;
  }
  std::vector<int> using_trigger_bits{12, 13, 22, 24, 25, 26, 27, 30, 31, 36, 37, 38};
  std::unique_ptr<CaloEvalStack> m_caloevalstack;
  float m_vertex_cut{200.0};
  float m_shower_shape_min_tower_E{0.07};
  CaloRawClusterEval *clustereval{nullptr};
  CaloTruthEval *trutheval{nullptr};
  HepMC::GenEvent *singal_event{nullptr};
  PHG4TruthInfoContainer *truthinfo{nullptr};
  CaloEvalStack *caloevalstack{nullptr};

  RawTowerGeomContainer *geomEM{nullptr};
  RawTowerGeomContainer *geomIH{nullptr};
  RawTowerGeomContainer *geomOH{nullptr};

  TowerInfoContainer *emcTowerContainer{nullptr};
  TowerInfoContainer *emcRawTowerContainer{nullptr};
  TowerInfoContainer *ihcalTowerContainer{nullptr};
  TowerInfoContainer *ihcalRawTowerContainer{nullptr};
  TowerInfoContainer *ohcalTowerContainer{nullptr};
  TowerInfoContainer *ohcalRawTowerContainer{nullptr};
  TowerInfoContainer *m_emc_sub1_tower_container{nullptr};
  TowerInfoContainer *m_ihcal_sub1_tower_container{nullptr};
  TowerInfoContainer *m_ohcal_sub1_tower_container{nullptr};

  RawClusterContainer *topoClusterContainer{nullptr};
  RawClusterContainer *topoClusterContainer_soft{nullptr};

  PHCompositeNode *topNodeptr{nullptr};
};

#endif // CALOANA24_H
