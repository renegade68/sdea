# sdea 5.8.0 for Stata

`sdea` estimates data envelopment analysis (DEA) models as observations arrive. It also computes radial technical efficiency for a cross-section using a reduced reference set.

## Run in Stata 16 or later

Download or clone sdea and the pinned [DEA-Stata 2.0.73 core](https://github.com/renegade68/DEA-Stata/tree/ad77287466b22c3edd1ae3cc074981225c95d424). Keep the core's complete `ado/` directory together. Open a fresh Stata session, change to the sdea repository root, and set `core` to that external directory. The loader adds the core followed by sdea to the front of the search path and prints both resolved commands.

```stata
local root "`c(pwd)'"
local core "../DEA-Stata/ado"
do "`root'/load.do" "`root'" "`core'"
help sdea
```

### Streaming arrivals and drift

```stata
do "`root'/examples/streaming.do" "`root'" "`core'"
```

This example orders observations by period and arrival, raises the simulated production frontier in period 7, and calls `method(buildhull)` with `driftalgo(frontier)`. Read `alpha` as hull distance, `frame` as a membership indicator, and positive `drift` as a signal from the selected diagnostic. The original DDM and EDDM detectors instead take an external binary error stream. The regression example explicitly supplies the simulated indicator with `driftmode(original) drifterror(err)`.

### Cross-sectional technical efficiency

```stata
do "`root'/examples/cross_section.do" "`root'" "`core'"
```

The four-firm example uses `method(standard) refset(buildhull)` to score all firms against a reduced VRS reference set. Firm 3 has input-oriented efficiency 0.75. Firms 1 and 2 are identical and efficient, and only one is needed as a reference. For output orientation, the reported score is the reciprocal of the radial expansion factor. Both examples accept an optional third argument naming an existing output directory for a CSV of their results.

## Reproduce the checks

Run these commands from the root with the same two locals. The scripts replace the dataset in memory.

```stata
do "tests/refset.do" "`root'" "`core'"
do "tests/benchmark_9040.do" "`root'" "`core'" 1
foreach check in validate_concept_drift_scenarios validate_concept_drift_heterogeneity validate_original_drift_detectors cert_frontier_gate cert_summarysave {
    do "tests/`check'.do" "`root'" "`core'"
}
```

The benchmark generates 9,040 firms with the mt64 generator and seed 119. It compares each firm's score under CRS and VRS, in both orientations, with `tests/fixtures/full_dea_9040.csv` at an absolute tolerance of 1e-7. The fixture contains the four full-sample DEA-Stata reference scores by firm. To reproduce those columns, run `examples/generate_9040.do` and use the pinned `dea x1 x2 x3 x4 x5 = y1 y2, dmu(firm) stage(1) engine(auto)` with the corresponding `rts()`, `orientation()` and `generate()` on a fresh copy for each model. The optional repeat count applies to the benchmark. `tests/validation_9040.csv` records two earlier runs on the development machine. Elapsed times depend on the host.

For the optional empirical comparison, supply the external data directory as a third argument to `examples/empirical.do`. File origins, variables and checksums are in the [data guide](data/README.md).

## Licensing

The external DEA-Stata core identifies its license as AFL-3.0. The sdea author's choice of a license for this independent package is pending.
