# Photon-jet histogram makers for PhotonJetTrees v21

The histogram makers of `histMakers/` (`RecoEffCalculator_TTreeReader.C`, `EvtCharacter.C`)
rewritten for Bennett's flat PhotonJetTrees v21 exports: the same histograms, names and
selections, with the event loop reading the flat trees instead of the slimtrees. pp is run with
the PPG12 analysis block, Au+Au with this repository's `nom.yaml` analysis block. Nothing else.

| File in histmakers/ | Responsibility |
| --- | --- |
| PhotonJetObjects.h | Event, which owns its photons, truth photons, jets and truth jets. |
| PhotonJetReader.h | Read the flat trees of one part, each once, and rebuild the events by their key. |
| PhotonJetConfig.h | The Config struct: every value read from the YAML, its validation, the centrality x eta bin layout. |
| PhotonJetSelection.h | Event, photon and jet cuts, ABCD classification, calibration, weights, smearing draws, jet matching. |
| PhotonJetHistograms.h | Booking and writing primitives, and the event counter shared by both macros. |
| PhotonJetIO.h | Configure a job, read the input list, publish a ROOT file, the PPG12 event weight. |
| RecoEffCalculator_TTreeReader.C | Spectra, ABCD regions, x_Jgamma, truth-efficiency counts and the matched responses. |
| EvtCharacter.C | Vertex, centrality, calorimeter energy, photon and jet characterization. |
| FinalizeMerged.C | After hadd: check the event count, write the full manifest, build efficiencies and responses. |
| CrossSectionWeights.h, TruthVertexReweightLoader.h | PPG12's sample weights, truth-pT windows and truth-vertex reweight, copied from the PPG12 repository. |

Each macro books its own histograms at the top of the file, one per line, followed by the event
loop with explicit Fill calls. To add an observable, add the member and its booking line there
and fill it in the loop. Every booked histogram is written.

## Environment

Source the sPHENIX environment (`source /opt/sphenix/core/bin/sphenix_setup.sh -n new`). The
configuration is read with yaml-cpp, which is not in the release: the macros load
`/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so` and need its headers and library on the paths,

~~~bash
export ROOT_INCLUDE_PATH=/sphenix/u/shuhang98/install/include:$ROOT_INCLUDE_PATH
export LD_LIBRARY_PATH=/sphenix/u/shuhang98/install/lib64:$LD_LIBRARY_PATH
~~~

`condor/run_photonjet_job.sh` does both for the jobs. RooUnfold comes with the release.
Run everything from this directory: the configs use paths relative to it.

## Inputs

The promoted v21 exports under `/sphenix/tg/tg01/bulk/jbennett/PhotonJetTrees_v2/`, each a
`parts/` directory with a `files.txt`: data/pp, data/auau, simulation/pp_photonjet,
simulation/pp_inclusive, simulation/pp_di, simulation/auau_embedded_photonjet and
simulation/auau_embedded_inclusive.

~~~bash
./make_lists.sh        # lists/<system>_<product>.list, one path per line, from the files.txt
~~~

The pp double-interaction parts mix photon and jet samples, so they are appended to both the
signal and the inclusive pp list, and each macro keeps the events of its product's kind of
sample. The lists are not tracked by git.

## Run

A direct call processes a whole list in one process:

~~~bash
root -l -b -q 'histmakers/RecoEffCalculator_TTreeReader.C("configs/pp/config_pj_pp_nom.yaml","data")'
root -l -b -q 'histmakers/EvtCharacter.C("configs/pp/config_pj_pp_nom.yaml","sim_signal")'
~~~

The second argument is `data`, `sim_signal` or `sim_inclusive`. The optional third argument is
one ROOT file or a list of them instead of the configured list, the fourth an output tag.
An existing output file is never overwritten.

For production, split a product into Condor jobs and merge when they are done:

~~~bash
./submit_photonjet.sh configs/pp/config_pj_pp_nom.yaml spectrum sim_signal 50   # 50 parts per job
./submit_photonjet.sh configs/pp/config_pj_pp_nom.yaml spectrum data 50 --local  # same chunks, run here
./merge_photonjet.sh  configs/pp/config_pj_pp_nom.yaml spectrum sim_signal
~~~

The core is `spectrum` (RecoEffCalculator_TTreeReader.C) or `events` (EvtCharacter.C). Chunk
lists go to `chunks/<var_type>/<core>_<product>/`, job logs to `logs/`. A product can be
submitted once per `var_type`: remove its chunk directory and outputs, or use a new `var_type`,
to redo it. The full pp set (three products, both cores, 1440 + 3172 jobs) runs in about 20
minutes, 3 minutes and 200 MB per job.

`merge_photonjet.sh` refuses to merge unless every submitted chunk wrote its file, then runs
hadd and `histmakers/FinalizeMerged.C`, which checks that the merged cutflow holds exactly the
events of the chunks, writes the full input manifest and builds the efficiencies and the
RooUnfoldResponse objects from the merged counts. Outputs land in `output.directory`
(`output/` by default): `<core>_<product>_<var_type>[.chunkNNNNN].root`, each embedding its
configuration and input manifest.

## Configs

`configs/pp/config_pj_pp_nom.yaml` is the PPG12 nominal analysis block (vertex cut 60 cm,
R = 0.3 isolation for reconstructed and truth photons, PPG12 shower-shape windows and BDT
thresholds by formula, resolution smearing of the response) plus the x_Jgamma keys of this
repository (`jet_cone_size`, `jet_eta`, `b2bjet_pT_min`, `b2bjet_dphi`, `xjgamma_bins`).
`configs/auau/config_pj_auau_nom.yaml` is `histMakers/nom.yaml`'s analysis block with a
BDT-based photon identification placeholder (`bdt_on: 1`, `n_nt_fail: 0`, thresholds from the
tree), `trigger_used: [22]` (Photon 10 GeV + MBD), the run range for the per-run histograms and
the resolution parameters. The `photonjet` block holds what those blocks cannot express: the
system, the part lists, whether thresholds come from the tree or the formulas, the jet
collection (`towerinfo_calib` in pp, `towerinfo_sub1` in Au+Au), the sample weights and the run
periods of the PPG12 blend. Analysis keys the code does not read are listed once at start-up.

Simulation weights follow PPG12 and use nothing from the producer's weight columns: the sample
of each source file comes from the part's `source_files` metadata, `CrossSectionWeights.h`
gives its cross-section ratio and truth-pT window (leading prompt truth photon, or leading
R = 0.4 truth jet for jet samples), and

    weight = sample weight * sum over periods ( lumi_fraction * SI or DI fraction * truth-vertex reweight )

with the periods and `reweight/<period>/reweight.root` files in `photonjet.periods`. Au+Au
simulation has no PPG12 table yet and runs with weight 1 (`sample_weights: 0`).

## Histograms

Names follow `histMakers/RecoEffCalculator_TTreeReader.C`: `h_<selection>_cluster_<c>` spectra
and ABCD regions with `_signal`, `_background` and `_notmatch` splits in simulation,
`h_truth_pT*_<c>` efficiency counts, `h_pT_reco/truth_response_<c>` and `h_response_full_<c>`,
the x_Jgamma set `h_<selection>_xjgamma_cent<c>` (plus `_signal`, `_truthmatchreco_..._signal`,
`_truthjet_..._signal`, `_xjgamma_background`), `h_jet_pT_response_cent<c>`,
`h_dphi_clusterJets_tight_cent<c>_pt<i>`, `h_tight/nontight_isoET_cent<c>_pt<i>`,
`h_pT_reco_fake_<c>`, `h_all/tight_cluster_signal_<c>`, `h_all_cluster_Et_max_b2bjet_<c>`,
`h_vertexz`, and after the merge `eff_reco/iso/id/all/vertex_cent<c>` and
`response_matrix_full_<c>`, `response_matrix_xjgamma_2d_cent<c>`. x_Jgamma is an observable of
the leading accepted photon: one entry per back-to-back jet (`|eta| < jet_eta`,
`|dphi| > b2bjet_dphi`, `pT > b2bjet_pT_min`) at jet pT / photon ET and an x-underflow entry when
there is none. Jets are matched geometrically (dR < 0.2, truth jets within 0.2 of a truth photon
excluded). The jet pT is used as stored, without the 1/0.65 of the reference. `EvtCharacter.C`
follows `histMakers/EvtCharacter.C` with the same names and binning.

## Known properties of the v21 trees (September 2026)

- The jets carry no energy scale correction: `jet_pt` equals `jet_raw_pt` in both systems, and the
  pp R = 0.3 reco / truth pT is 0.59 at 5-8 GeV rising to 0.69 at 30-50 GeV (`h_jet_pT_response_cent0`).
- The photon rows of `recoTruthLinks` are a nearest-neighbour dR < 0.1 match at the reconstructed
  vertex. In double-interaction events that vertex is within 3 cm of the hard-scatter vertex only
  22 % of the time, the photon eta moves by up to 0.4, and two thirds of the DI photons lose their
  link. Until the producer rebuilds the links from the truth track id (`dominant_truth_track_id`
  on the photon rows, `native_track_id` on the truth rows), the blended pp reconstruction
  efficiency comes out at 0.80 where PPG12 has 0.92. Isolation and identification efficiencies,
  the data and the single-interaction samples are unaffected.
- The isolation columns are the producer's own cones (R = 0.3 bound definition, R = 0.4 witness)
  with PPG12's a = 0.49, b = 0.037, not PPG12's topo-cluster R = 0.4 column. With `cone_size: 3`
  the pp data ABCD yields are 1.08 (A), 0.87 (B), 1.17 (C), 0.89 (D) of PPG12 nominal, the
  isolation-independent `h_all_cluster` and `h_tight_cluster` agree to 1 %.
- The pp export holds photons down to 5 GeV and |eta| < 0.7 at the reconstructed vertex.
