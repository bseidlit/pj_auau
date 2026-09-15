# Pipeline source and validation

This is an additive copy of `efficiencytool/photonjet` and its required shared
helpers from PPG12 commit `638ddc8acb711a8f532a08567d8409755baa5912` on
`PPG12-2-simplify-histogram-maker` (2026-09-15). Source files and their hashes are
listed in [SOURCE.json](SOURCE.json) for the import snapshot at copy commit
`0f95772e2eccf31eeb5b8e844e8d21867d6ef7d3`; later changes are tracked in Git.
The original PPG12 checkout is retained.
The target baseline is pj_auau main `a9a9ea125ac54d5c2b609a40dad59cb12b755074`.

## Relocation changes

- Bundle `CrossSectionWeights.h` and `CalculatePhotonYield.C` in `support/`;
  adjust include paths and provenance inventories to use those copies.
- Keep output/log/figure paths within the new package and use included sample maps.
- Load yaml-cpp by library name from the configured environment; allow
  `PJ_SETUP_SCRIPT` for Condor workers, retaining the original SDCC default.
- Copy current scripts, configs, metadata, tests, frozen C++ reference and plot
  consumers. Keep private bookkeeping, generated data and retired RDF design and
  benchmark documentation in the source repository.

Histogram contents, selections, weights, centrality/eta layout and keyed response
RNG are unchanged. External scientific data retain their original identities;
no calibration or normalization is substituted during relocation. AuAu final
normalization remains outside the pp yield command's supported scope.

## Import validation

At clean copy commit `957b532d11191af062a717d46e626b919b6c8872`, using ROOT 6.32.06,
RooUnfold 3.1.0 and Python 3.13.0:

- All 39 reader/physics/workflow fixtures passed in 95.484 seconds, including
  centrality/eta grids, MBD, pp yield/closure, canonical RNG vectors and checksum
  safety checks.
- Six main/response histogram files matched saved source outputs exactly,
  including values, errors, moments, axes and titles.
- An actual Condor dry run generated six jobs with one CPU each. All six worker
  arguments matched the planned config/product/list/tag, copied directory and
  custom setup script. No jobs were submitted.
- A shipped AuAu config planned successfully when the wrapper was invoked from
  outside the package, using this copy's sample map.
- All pre-existing pj_auau files remained unchanged. The checkpoint closed with
  zero changed dependencies; private checkpoint files are excluded from this copy.

At the import revision, the event loop, reader, histogram booking, configuration
implementation and cross-section helper matched the source byte for byte. The
physics helper differed only in its include path to the bundled cross-section helper. The source had separately
passed broader real-input and exact response RNG comparisons before copying.
These checks establish bounded equivalence, not full-release throughput or AuAu
final-normalization validation.

## Readability follow-up

Production histogram code now spells out analysis records, histogram collections,
selection/response results and ROOT ownership types where they explain the local
code. Lambdas, iterators and declarations with obvious types retain `auto`. This
changes type spelling only, preserving const/reference qualifiers and ownership.
The hashes in SOURCE.json continue to describe the import snapshot above.
