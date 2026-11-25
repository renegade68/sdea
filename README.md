# sdea 5.8.0 for Stata

`sdea` estimates data envelopment analysis (DEA) models as observations arrive. It also computes radial technical efficiency for a cross-section using a reduced reference set.

## Installation

You need **Stata 16 or later** and **DEA-Stata 2.0.73**. Stata includes Mata and its interior-point solver. DEA-Stata supplies the simplex solver, its helper files, and the `dea` command used for reference scores. These are the runtime requirements for sdea and its synthetic examples.

1. Download the sdea repository using **Code → Download ZIP** and extract it, or clone it with Git. Open a fresh Stata session and set the working directory to the folder containing `load.do`.
2. In Stata, run the following command **once** to download and extract the tested DEA-Stata release into that folder. For an offline installation, download the same ZIP in a browser and extract its top-level folder beside `load.do`, keeping the folder name and contents intact. [Download DEA-Stata 2.0.73 at the tested commit](https://github.com/renegade68/DEA-Stata/archive/ad77287466b22c3edd1ae3cc074981225c95d424.zip)

   ```stata
   unzipfile "https://github.com/renegade68/DEA-Stata/archive/ad77287466b22c3edd1ae3cc074981225c95d424.zip"
   ```

3. Run the following commands to load sdea and open the help. The loader checks the core version and solver files, adds both packages to the current session’s search path, and prints the locations of `dea` and `sdea`. In a new session, run `do load.do` again from the sdea folder.

   ```stata
   do load.do
   help sdea
   ```

If you already keep the pinned core elsewhere, supply the sdea folder as the first argument and the full path to the core’s `ado` folder as the second argument to `load.do` and to the example or test you run. With the installation above, omit both arguments.

### Streaming arrivals and drift

```stata
do examples/streaming.do
```

This example orders observations by period and arrival, raises the simulated production frontier in period 7, and calls `method(buildhull)` with `driftalgo(frontier)`. Read `alpha` as hull distance, `frame` as a membership indicator, and positive `drift` as a signal from the selected diagnostic. The original DDM and EDDM detectors instead take an external binary error stream. The regression example explicitly supplies the simulated indicator with `driftmode(original) drifterror(err)`.

### Cross-sectional technical efficiency

```stata
do examples/cross_section.do
```

The four-firm example uses `method(standard) refset(buildhull)` to score all firms against a reduced VRS reference set. Firm 3 has input-oriented efficiency 0.75. Firms 1 and 2 are identical and efficient, and only one is needed as a reference. For output orientation, the reported score is the reciprocal of the radial expansion factor. Both examples accept an optional third argument naming an existing output directory for a CSV of their results.

## Citing BuildHull

For the two-phase DEA algorithm that first identifies the extreme efficient units (the frame) and then scores units against that frame, cite J. H. Dulá (2011), “An Algorithm for Data Envelopment Analysis,” *INFORMS Journal on Computing* 23(2), 284–296. The `bhvariant(dula)` option follows the phase-one search using a separating hyperplane. The default `bhvariant(modified)` adds exterior arrivals to the current frame and supports optional pruning. [Dulá (2011)](https://doi.org/10.1287/ijoc.1100.0400)

For frontier analysis as observations arrive, cite J. H. Dulá and F. J. López (2013), “DEA with streaming data,” *Omega* 41(1), 41–47. When reporting sdea results, specify the method and variant so readers can identify the implementation used. [Dulá and López (2013)](https://doi.org/10.1016/j.omega.2011.07.010)

## Reproduce the checks

Run these commands from the folder containing `load.do`. The scripts replace the dataset in memory.

```stata
do tests/refset.do
do tests/benchmark_9040.do
foreach check in validate_concept_drift_scenarios validate_concept_drift_heterogeneity validate_original_drift_detectors cert_frontier_gate cert_summarysave {
    do "tests/`check'.do"
}
```

The benchmark generates 9,040 firms with the mt64 generator and seed 119. It compares each firm's score under CRS and VRS, in both orientations, with `tests/fixtures/full_dea_9040.csv` at an absolute tolerance of 1e-7. The fixture contains the four full-sample DEA-Stata reference scores by firm. To reproduce those columns, run `examples/generate_9040.do` and use the pinned `dea x1 x2 x3 x4 x5 = y1 y2, dmu(firm) stage(1) engine(auto)` with the corresponding `rts()`, `orientation()` and `generate()` on a fresh copy for each model. The benchmark accepts a repeat count as its third argument. `tests/validation_9040.csv` records two reference runs with Stata/MP 16.1. Elapsed times depend on the host.

For the optional empirical comparison, supply the external data directory as a third argument to `examples/empirical.do`. File origins, variables and checksums are in the [data guide](data/README.md).

## Licensing

Copyright © 2026 Choonjoo Lee. Licensed under the Academic Free License version 3.0 (AFL-3.0). The license covers the sdea code, documentation, examples, tests and project-generated benchmark files. [Full license text](LICENSE)

The separately downloaded DEA-Stata core carries its own AFL-3.0 notice. External empirical datasets retain their source-specific terms. [Example data and redistribution](data/README.md)

For bug reports and questions about sdea, contact Choonjoo Lee at bloom.rampike@gmail.com.
