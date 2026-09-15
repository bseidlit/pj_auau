"""Tiny PhotonJetTrees-format files for behavioral regression tests."""
from array import array
import ROOT

def write_tree(f, name, schema, rows):
    f.cd()
    t=ROOT.TTree(name,name)
    buffers={key:array({"I":"i","L":"q","l":"Q","D":"d"}[kind],[0]) for key,kind in schema.items()}
    for key,kind in schema.items(): t.Branch(key,buffers[key],key+"/"+kind)
    for row in rows:
        for key in schema: buffers[key][0]=row.get(key,0)
        t.Fill()
    t.Write()

event_schema=dict(source_file_index="I",source_entry="L",event_id_hi="l",event_id_lo="l",run="I",trigger_bits="l",live_trigger_bits="l",scaled_trigger_bits="l",scaled_bit30="I",vertex_z="D",centrality="D",event_weight="D",total_calo_energy="D",terminal_status="I")
photon_ints="source_file_index run scaled_bit30 terminal_status photon_encounter_ordinal bdt_is_tight bdt_is_nontight bdt_is_not_tight bdt_input_count iso_r03_pass iso_r04_pass truth_matched truth_barcode".split()
photon_doubles="vertex_z centrality event_weight photon_et photon_eta photon_phi bdt_score bdt_tight_threshold bdt_nontight_low_threshold bdt_nontight_high_threshold iso_r03 iso_r03_threshold iso_r03_nonisolated_threshold iso_r04 iso_r04_threshold iso_r04_nonisolated_threshold native_weta_cogx native_wphi_cogx native_weta33_cogx native_wphi33_cogx native_e11_over_e33 native_e32_over_e35 native_et1".split()
photon_schema={**{k:"I" for k in photon_ints},**{k:"D" for k in photon_doubles},"event_id_lo":"l",**{"bdt_input_%02d"%i:"D" for i in range(14)}}
truth_schema=dict(source_file_index="I",event_id_lo="l",truth_photon_pt="D",truth_photon_eta="D",truth_photon_phi="D",prompt_class="I",source_role="I",generator_barcode="I",truth_isolation="D")
ev=dict(source_file_index=0,source_entry=0,event_id_hi=2,event_id_lo=1,run=50000,scaled_trigger_bits=1<<30,scaled_bit30=1,vertex_z=1.,centrality=10.,event_weight=2.)
ph={**ev,"photon_encounter_ordinal":0,"photon_et":20.,"bdt_score":.6,"bdt_tight_threshold":.784375,"bdt_nontight_low_threshold":.4666667,"bdt_nontight_high_threshold":.715625,"bdt_is_tight":0,"bdt_is_nontight":1,"bdt_is_not_tight":1,"bdt_input_count":11,"iso_r03":0.,"iso_r04":0.,"iso_r03_threshold":1.23,"iso_r04_threshold":1.23,"iso_r03_nonisolated_threshold":2.03,"iso_r04_nonisolated_threshold":2.03,"iso_r03_pass":1,"iso_r04_pass":1,"truth_matched":1,"truth_barcode":42}
for i in range(14): ph["bdt_input_%02d"%i]=.5
ph["bdt_input_06"]=.7
ph["bdt_input_10"]=.9
ph2={**ph,"photon_encounter_ordinal":1,"bdt_score":.9,"bdt_is_tight":1,"bdt_is_nontight":0,"bdt_is_not_tight":0,"iso_r03":5.,"iso_r04":5.,"iso_r03_pass":0,"iso_r04_pass":0}
tr=dict(source_file_index=0,event_id_lo=1,truth_photon_pt=20.,truth_photon_eta=0.,truth_photon_phi=0.,prompt_class=1,source_role=1,generator_barcode=42,truth_isolation=0.)
def make_fixture(path, photons=None, event_types=None):
    f=ROOT.TFile(str(path),'RECREATE')
    write_tree(f,'events',event_types or event_schema,[ev])
    write_tree(f,'photons',photon_schema,[ph,ph2] if photons is None else photons)
    write_tree(f,'truthPhotons',truth_schema,[tr])
    f.Close()

MBD_NAMES=['h_truth_pT_vertexcut_0','h_truth_pT_vertexcut_mbd_cut_0',
 'h_truth_pT_vertexcut_mbd_north_cut_0','h_truth_pT_vertexcut_mbd_south_cut_0',
 'h_truth_pT_vertexcut_mbd_only_north_0','h_truth_pT_vertexcut_mbd_only_south_0',
 'h_truth_pT_vertexcut_mbd_neither_0']

def correction_file(path, edges=(8,22,45), names=MBD_NAMES, corrected=False):
    f=ROOT.TFile(str(path),'RECREATE')
    counts=[100,40,30,25,15,10,60] if corrected else [100]*7
    for name,content in zip(names,counts):
        h=ROOT.TH1D(name,name,len(edges)-1,array('d',edges))
        for i in range(1,len(edges)): h.SetBinContent(i,content)
        h.Write()
    marker=ROOT.TH1D('unrelated','preserve this object',1,0,1)
    marker.SetBinContent(1,123.)
    marker.Write()
    f.Close()
