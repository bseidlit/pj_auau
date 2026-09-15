# PhotonJetTrees_v1 → PPG12 histograms

The maker reads original ROOT parts and fills histograms in visible C++ event,
reco-photon and truth-photon loops. Histograms live in a typed
`histograms[centrality][eta]` array, booked from the configured bin edges.
Each worker uses one CPU; Condor runs independent jobs in parallel.

This directory is a self-contained copy of the direct PhotonJetTrees pipeline
from PPG12. It runs alongside the existing `histMakers/RecoEffCalculator_TTreeReader.C`
workflow, which reads different input trees. See [COPYING.md](COPYING.md) for the
source revision and relocation details.

## Run

From `histMakers/photonjet`, after sourcing the analysis environment:

```bash
export PJ_SETUP_SCRIPT=/sphenix/u/shuhang98/setup.sh
source "$PJ_SETUP_SCRIPT"

bash oneforall_photonjet.sh configs/pp/config_pj_pp_nom.yaml chunks
bash oneforall_photonjet.sh configs/pp/config_pj_pp_nom.yaml submit
bash oneforall_photonjet.sh configs/pp/config_pj_pp_nom.yaml status
bash oneforall_photonjet.sh configs/pp/config_pj_pp_nom.yaml merge
bash oneforall_photonjet.sh configs/pp/config_pj_pp_nom.yaml check
bash oneforall_photonjet.sh configs/pp/config_pj_pp_nom.yaml yield
```

Use `local` in place of `submit` to run the same assigned jobs sequentially.
Use Python 3.11 or newer with `ruamel.yaml`, ROOT with PyROOT, yaml-cpp and
RooUnfold on the library/include paths. The validated environment uses ROOT
6.32.06 and RooUnfold 3.1.0. Set `PJ_SETUP_SCRIPT` to the setup script available on
each Condor node; the submit wrapper passes it to workers. The SDCC default is
shown above. Python dependencies are listed in [requirements.txt](requirements.txt). Planning needs no ROOT reader;
processing and merging use the configured ROOT, yaml-cpp and RooUnfold environment.
Regression tests also use `uproot` for independent file inspection.

`chunks` writes indexed ROOT lists and seals their assignments in `jobs.json`
format 4. It records each original's canonical index/identity, path, size and
nanosecond modification time without opening the ROOT file. Configuration,
master-list order, code, histogram layout, job lists and scientific dependencies
are still hashed. Default chunks contain 200 data parts or 40 MC parts.

Each worker checks its assigned files, computes their full MD5 checksums before
and after processing, and records ROOT UUIDs and tree counts. Schema errors fail
the worker; duplicate UUIDs within a job fail there, and duplicates across jobs
fail merging. Workers publish a checked output generation with an exact binding
to the job plan.

Status, merge, check and yield verify output checksums, current dependencies,
exact assignments and coverage. They check original sizes/timestamps without
reading original contents. This workflow assumes a finalized, immutable input
release: the plan seals file metadata, while each worker records the bytes it
processed. A content edit that preserves both size and timestamp before a worker
starts is outside the planning guarantee. Regenerate chunks and outputs when
inputs change. No extra validation mode is needed.

Use a scratch config with separate output prefixes and `var_type` for bounded
work. For example:

```bash
bash oneforall_photonjet.sh /tmp/my-pj-config.yaml chunks \
  --products data,sim_signal --indices 0,4:7 \
  --chunk-size data=2,sim_signal=2
bash oneforall_photonjet.sh /tmp/my-pj-config.yaml local \
  --products data,sim_signal
```

Indices are zero-based positions in the original master lists; ranges exclude
their upper bound. Keep those lists and their ordering. A subset never defines
MC sample normalization. `check` and `yield` expect merged data, signal MC and
inclusive MC. They require full coverage unless `--allow-subset` explicitly
selects a bounded validation run.

`PJ_CHUNK_DIR` and `PJ_LOG_DIR` override job-list and log directories.
`PJ_RELEASE_ID` optionally overrides the release label, whose default is the
basename of `photonjet.release`. Keep that label stable: it is part of the response
RNG identity. Physics settings remain in YAML. `PJ_THREADS` defaults to 1; values
other than 1 are rejected. The former `prepare` step and prepared-store options are retired.

ROOT output names remain `{data_outfile,eff_outfile,response_outfile}_{var_type}`;
inclusive MC adds `_jet` to the efficiency prefix. Chunk ROOT files live in the
output directory's `chunks/` subdirectory. Input/job lists default to
`photonjet/chunks/`, stdout/stderr and yield logs to `photonjet/logs/`, and
Condor event logs to `/tmp`. AuAu histogram production uses the same commands;
the current yield/closure command supports pp only.

### Direct ROOT invocation

```bash
root -l -b -q \
  'histmakers/PhotonJetHistMaker.C("config.yaml","sim_signal","parts.list","check")'
```

The arguments are `(config, product, chunk_list = "", chunk_tag = "")`.
An empty input list uses the configured master list. A nonempty list contains
`original_index original_ROOT_path` rows, preserving canonical indices even for
subsets. When the list belongs to a sealed `jobs.json`, its assignment and tag
must match. Direct calls also check source identities, dependencies, output
collisions and complete processing before publication.

## Code structure

| File | Responsibility |
| --- | --- |
| [`histmakers/PhotonJetHistMaker.C`](histmakers/PhotonJetHistMaker.C) | Configuration, booking, visible event/reco/truth loops, direct ROOT fills and publication |
| [`histmakers/PhotonJetReader.h`](histmakers/PhotonJetReader.h) | Ordinary records, checked original tree/schema/identity loading and associations |
| [`histmakers/PhotonJetPhysics.h`](histmakers/PhotonJetPhysics.h) | Named cuts, matching, event/candidate weights and keyed response calculations |
| [`histmakers/PhotonJetHistograms.h`](histmakers/PhotonJetHistograms.h) | Typed cell ownership, booking/names/writing and per-cell RooUnfold merge correction |
| [`histmakers/PhotonJetConfig.h`](histmakers/PhotonJetConfig.h) | Configuration validation, bin layout and boundary helpers |
| [`histmakers/PhotonJetIO.h`](histmakers/PhotonJetIO.h) | Original inputs, model setup, dependencies, provenance and checked publication |
| [`histmakers/InjectMbdEff.C`](histmakers/InjectMbdEff.C) | Checked seven-histogram external MBD/vertex correction import |
| `scripts/{inputs,make_chunks,workflow}.py`, `scripts/CheckChunks.C` | Sealed original-input jobs, serial worker, status, checked merge and yield |
| [`tests/reference/`](tests/reference/) | Frozen older scalar-loop code used only for regression comparisons |
| `support/{CrossSectionWeights.h,CalculatePhotonYield.C}` | Bundled cross sections and pp yield/closure consumer |
| `plotting/`, `metadata/` | Plot consumers, luminosity/run lists and canonical sample maps |

A cut changes in its named physics helper. Adding a histogram changes its typed
member/booking and its visible fill in the maker. Adding a centrality class changes
the bin definitions; the existing booker creates the complete inventory for it.

## Centrality × eta bins

Eta edges come from `analysis.eta_bins`; AuAu centrality edges come from
`photonjet.centrality_bins`. Both accept finite, strictly increasing edge arrays.
For example, these edits to an AuAu config create six histogram cells:

```yaml
analysis:
  eta_bins: [-0.7, 0.0, 0.7]
photonjet:
  centrality_bins: [0, 20, 50, 80]
```

Centrality intervals use `[low, high)`. The outer eta acceptance remains strict:
`first_edge < eta < last_edge`; shared interior edges belong to the higher eta
bin. pp has one inclusive centrality slot and requires no centrality edges.
The shipped configs retain their existing single-cell selections.

One-cell outputs preserve the legacy names and titles. Multi-cell objects use
`_cent{centrality_index}_eta{eta_index}` suffixes, for example
`h_all_cluster_cent0_eta0` and `eff_reco_cent0_eta0`. Every cell owns its reco/ABCD
spectra, photon QA, truth spectra, efficiencies and signal response. Matching and
conditional efficiencies use reco and truth from the same cell; a cross-eta
association follows the unmatched policy in the reco cell. Cell indices never
enter the RNG key.

Seven job-level QA objects and source/event accounting stay outside the array.
Selected events are counted once per centrality row in metadata, independently
of the number of eta cells. Each part is loaded once; extra cells do not duplicate
input rows or rerun the event loop. AuAu tree thresholds retain their stored
centrality dependence.

The direct maker and checked merge support multiple cells. Yield, MBD injection
and local plot consumers reject unsupported multi-cell inputs. Explicit cell
selection and class-specific normalization are required before extending those
consumers to future AuAu production.

## Physics and completeness

- The reader keeps zero-reco events, original candidate order and producer sentinel
  values. It checks full event identities, branch types, repeated event fields,
  orphan/duplicate objects and association bounds before event selection.
- Stored pp MC weights equal to 1 are producer placeholders and are skipped.
  Recomputed/sample-map weights use canonical sample metadata and original indices.
  AuAu simulation requires sample-map weighting.
- Data run/trigger/vertex/centrality and skim-floor checks remain strict. pp formula
  thresholds are checked against stored flags; `flag_check: report` supports
  deliberate threshold variations. Requested missing link inputs are errors.
- MC scale/nonlinearity affect selection. Extra resolution affects only matched
  tight-and-isolated response entries, with width `SigmaExtraFrac(truth_pt)*truth_pt`.
- Each fiducial truth photon contributes once to an applicable efficiency.
  Conditional ID and combined efficiency require one candidate passing both tight
  ID and isolation. Each qualifying matched candidate fills the response; there
  are no Miss/Fake entries.
- The response RNG retains `fnv1a64-splitmix64-mixmax17-gaus-v1`, keyed by the
  configured seed and original release/product/part/event/candidate identity.
  File order, chunking and histogram cell indices do not change its draws.
- One-cell main outputs contain 40 physics/QA objects plus four metadata objects;
  signal response files contain four response objects plus four metadata objects.
  For `N` cells, the counts are `7 + 33*N + 4` and `4*N + 4`, respectively.
  Conversion efficiency is empty; `novtx` duplicates weighted truth;
  vertex/MBD histograms remain unit placeholders until correction import.

Every main output has a `.complete.yaml` marker binding its checksum, UUID and
generation to the corresponding response file when present. The backend is
`photonjet-cpp-v1`, with original schema `PhotonJetTrees_v1`, provenance format 4
and completion format 3, using `worker-md5-stat-v1`. Checkers enforce the configuration, current code,
layout/schema/RNG, dependencies and exact disjoint original-part coverage.
Old plans and outputs must be regenerated; verification contracts cannot be mixed.

Merging stages all products and MBD work before replacing a prior generation.
It merges efficiency passed/total counts, restores neutral global sample weight
and rebuilds RooUnfold from merged histograms independently for each cell, applying
constructor entry offsets once per complete sample. Yield stages both data and
MC closure outputs and checks finite final spectra. Processing failures preserve
prior files; an interrupted publication leaves no valid completion marker for
a mixed generation. There is no incomplete-merge override.

MBD import requires all seven histograms, matching truth axes, valid counts and
compatible archived selections. Empty `external_mbd_eff_file` keeps the unit
placeholders. The release lacks the truth/MBD information needed to derive that
correction directly. Trigger bits cannot replace luminosity scalers. AuAu
normalization still needs external MB event counts and T_AA; its `lumi: 1.0`
remains an exploratory placeholder.

Metadata archives the parsed YAML, resolved layout, original identities/counts,
ROOT version, worker-recorded content fingerprints and code/dependency identities.
MD5 detects content changes when a file is hashed; subsequent original freshness
checks use size/timestamp. Analysis checkpoints additionally freeze SHA-256 evidence.
These records describe the declared coverage, not release-wide normalization.

## External inputs and paths

Run commands from this package directory. The shell wrapper changes into it
before launching stages. Shipped output prefixes use `results/`, and shipped
sample maps use the included `metadata/` directory. The input master lists retain
their original SDCC release paths and canonical ordering.

Original ROOT parts, calibration ROOT files and luminosity scaler lists are
external data. Configure their locations in your YAML for your environment.
The pp `external_mbd_eff_file` still identifies the original PPG12 correction
input; provide that validated file or an equivalent explicitly reviewed
correction. The copy does not replace it with an efficiency of one. Configs
retain historical unused fields from their PPG12 source; physics values and
selection/weight settings are preserved.

Plot macros use `results/` and `figures/` relative to this directory. Create
`figures/` before plotting. Comparison macros also require the relevant PPG12
reference outputs in `results/`; those generated data files are not included.

Regenerate chunk plans and outputs in this checkout: code fingerprints include
the copied support files and must match the checkout that produces the results.

## Validation

```bash
python3 -m unittest discover -s tests -p 'test_cpp_*.py' -v
```

The suites cover raw schema/identity failures, selection and weight boundaries,
conditional efficiencies, response/RNG vectors, centrality × eta routing,
partitioned merging, stale coverage/dependencies, publication failures, MBD,
populated pp yield/closure and worker-side checksum verification. The copy's
validation results and its comparison to the source are recorded in
[COPYING.md](COPYING.md).
