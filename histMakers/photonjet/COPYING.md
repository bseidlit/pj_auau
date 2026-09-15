# Pipeline source and validation

This is an additive copy of `efficiencytool/photonjet` and its required shared
helpers from PPG12 commit `638ddc8acb711a8f532a08567d8409755baa5912` on
`PPG12-2-simplify-histogram-maker` (2026-09-15). Source files and their hashes are
listed in [SOURCE.json](SOURCE.json). The original PPG12 checkout is retained.
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

## Validation

Pending clean-checkout validation of this copy. The source passed 39 fixtures,
14 exact numerical comparisons and four exact response audits before copying.
