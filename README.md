# TritonPart: A Constraints-Driven General Hypergraph Partitioner

TritonPart is a constraints-driven, multilevel hypergraph partitioner implemented within the [OpenROAD](https://github.com/The-OpenROAD-Project/OpenROAD) application.

The partitioner supports multiway partitioning, multidimensional vertex and hyperedge weights, fixed vertices, placement-aware partitioning, community constraints, and timing-driven optimization.

> **Development status:** TritonPart is maintained inside OpenROAD. Please refer to [par](https://github.com/The-OpenROAD-Project/OpenROAD/tree/master/src/par) for the latest implementation. All results reported in this README use OpenROAD v2.0-27290-g5e8f15c075.
> Experimental results should be reported together with the exact source commit, command-line configuration, random seeds, and balance constraint used.

---

## Table of Contents

* [Repository Organization](#repository-organization)
* [Features](#features)
* [Reproducibility of the Reported Titan23 Cuts](#reproducibility-of-the-reported-titan23-cuts)

  * [Pinned Tool Version](#pinned-tool-version)
  * [Experiment Summary](#experiment-summary)
  * [Exact TritonPart Command](#exact-tritonpart-command)
  * [Explicitly Set Parameters](#explicitly-set-parameters)
  * [Parameters Left at Their Defaults](#parameters-left-at-their-defaults)
  * [Important Balance-Constraint Semantics](#important-balance-constraint-semantics)
  * [Best-of-20 Titan23 Results](#best-of-20-titan23-results)
* [Building TritonPart](#building-tritonpart)
* [Basic Hypergraph Partitioning](#basic-hypergraph-partitioning)
* [Fixed-Vertex Constraints](#fixed-vertex-constraints)
* [Embedding and Placement Constraints](#embedding-and-placement-constraints)
* [Running the Regression Suite](#running-the-regression-suite)
* [Reporting TritonPart Results](#reporting-tritonpart-results)

---

## Repository Organization

The main repository components are:

* [`OpenROAD`](https://github.com/ABKGroup/TritonPart/tree/main/OpenROAD)
  Contains the OpenROAD source tree used to build TritonPart. The TritonPart implementation is located primarily in:

  ```text
  OpenROAD/src/par/
  ```

* [`titan23_benchmark`](https://github.com/ABKGroup/TritonPart/tree/main/titan23_benchmark)
  Contains the Titan23 hypergraph benchmark suite.

* [`regression`](https://github.com/ABKGroup/TritonPart/tree/main/regression)
  Contains scripts for running TritonPart regression experiments.

---

## Features

TritonPart provides a general-purpose hypergraph partitioning framework with support for:

* Multiway hypergraph partitioning
* Multilevel coarsening and refinement
* Multiple vertex-weight dimensions
* Multiple hyperedge-weight dimensions
* User-configurable partition-balance constraints
* Fixed-vertex constraints
* Community constraints
* Group constraints
* Placement-aware and embedding-aware partitioning
* Timing-driven partitioning
* V-cycle refinement
* FM-based refinement
* ILP-based optimization for sufficiently small subproblems
* User-configurable cost-function weights
* Reproducible randomized runs through explicit seeds
* Integration with the OpenROAD application
* A permissive open-source license

---

# Reproducibility of the Reported Titan23 Cuts

This section records the exact TritonPart source version, configuration, random seeds, benchmark set, execution environment, and result-selection procedure used to generate the Titan23 cuts reported later in the README.

---

## Pinned Tool Version

The experiments were generated using the following OpenROAD/TritonPart revision.

| Item                   | Value                                                                             |
| ---------------------- | --------------------------------------------------------------------------------- |
| OpenROAD version       | `OpenROAD v2.0-27290-g5e8f15c075`                                                 |
| Full source commit     | `5e8f15c075ddd9d45d750977ad9a856f5c1c2569`                                        |
| Short source commit    | `5e8f15c075`                                                                      |
| Branch                 | `master`                                                                          |
| Commit date            | `2025-12-22`                                                                      |
| TritonPart Tcl command | `triton_part_hypergraph`                                                          |
| Source directory       | `OpenROAD/src/par/`                                                               |

The source revision can be verified from a local checkout using:

```bash
git -C OpenROAD rev-parse HEAD
git -C OpenROAD describe --always --long --dirty
git -C OpenROAD show --no-patch --format=fuller \
    5e8f15c075ddd9d45d750977ad9a856f5c1c2569
```

To reproduce the reported cuts exactly, use:

```text
5e8f15c075ddd9d45d750977ad9a856f5c1c2569
```

Do not assume that results obtained from a newer TritonPart or OpenROAD revision are directly comparable. Changes to coarsening, refinement, random-number generation, balance checking, or default parameters may change the resulting cutsize.

---

## Tool and Execution Environment

| Item                       | Value                                                                        |
| -------------------------- | ---------------------------------------------------------------------------- |
| Hypergraph format          | hMETIS                                                                       |
| Hypergraph location        | `benchmark/Titan23_benchmark/<design>.hgr`                                   |
| Number of partitions       | 2                                                                            |
| Number of seeds            | 20 per design                                                                |
| Seeds                      | `0` through `19`                                                             |
| Number of designs          | 22                                                                           |
| Total number of runs       | 440                                                                          |
| Completed runs             | 440                                                                          |
| Cutsize metric             | `Cutcost`                                                                    |
| Cutsize reporting          | `set_debug_level PAR evaluation 1`                                           |
| Result selected per design | Minimum `Cutcost` over 20 seeds                                              |

The experiments did not use timing constraints, placement files, fixed vertices, community files, or group files.

---

## Experiment Summary

| Property                          | Configuration                 |
| --------------------------------- | ----------------------------- |
| Benchmark suite                   | Titan23                       |
| Designs                           | 22                            |
| Partitioning problem              | 2-way hypergraph partitioning |
| Requested balance constraint      | `0.1`                         |
| Effective balance interpretation  | `0.1%`                        |
| Effective target                  | Nearly exact 50/50 bisection  |
| Random seeds                      | `0–19`                        |
| Runs per design                   | 20                            |
| Total runs                        | `22 × 20 = 440`               |
| Initial solutions per run         | 50                            |
| Best initial solutions retained   | 10                            |
| Coarsening solutions              | 3                             |
| Refinement iterations             | 10                            |
| Maximum moves per refinement pass | 100                           |
| Maximum V-cycles                  | 3                             |
| Final reported result             | Best cutsize over 20 seeds    |

---

## Exact TritonPart Command

Every run used the following Tcl command. Only `-hypergraph_file` and `-seed` changed between runs.

```tcl
set_debug_level PAR evaluation 1

triton_part_hypergraph \
    -hypergraph_file <design>.hgr \
    -num_parts 2 \
    -balance_constraint 0.1 \
    -seed <0..19> \
    -num_initial_solutions 50 \
    -num_best_initial_solutions 10 \
    -num_coarsen_solutions 3 \
    -refiner_iters 10 \
    -max_moves 100 \
    -early_stop_ratio 0.5 \
    -v_cycle_flag true \
    -max_num_vcycle 3 \
    -thr_coarsen_hyperedge_size_skip 200 \
    -thr_coarsen_vertices 10 \
    -thr_coarsen_hyperedges 50 \
    -coarsening_ratio 1.6 \
    -max_coarsen_iters 30

exit
```

A representative shell invocation is:

```bash
singularity exec \
    -B /home \
    -B /tmp \
    /home/tool/singularity/images/ispd26.sif \
    /home/fetzfs_projects/TritonPart/OpenROAD/build/bin/openroad \
    -exit < run.tcl
```

---

## Explicitly Set Parameters

The following parameters were explicitly passed to `triton_part_hypergraph`.

The defaults listed below correspond to:

```text
OpenROAD/src/par/src/partitionmgr.tcl
```

at commit:

```text
5e8f15c075ddd9d45d750977ad9a856f5c1c2569
```

| Parameter                          | Value used | Default | Differs from default? | Description                                                              |
| ---------------------------------- | ---------: | ------: | :-------------------: | ------------------------------------------------------------------------ |
| `-num_parts`                       |        `2` |     `2` |           No          | Number of output partitions                                              |
| `-balance_constraint`              |      `0.1` |   `1.0` |        **Yes**        | Allowed imbalance, expressed in percent                                  |
| `-seed`                            |     `0–19` |     `0` |        **Yes**        | Random seed; 20 runs per design                                          |
| `-num_initial_solutions`           |       `50` |    `50` |           No          | Number of randomized initial solutions generated at the coarsest level   |
| `-num_best_initial_solutions`      |       `10` |    `10` |           No          | Number of best initial solutions retained for subsequent refinement      |
| `-num_coarsen_solutions`           |        `3` |     `3` |           No          | Number of independent coarsening solutions considered per V-cycle        |
| `-refiner_iters`                   |       `10` |    `10` |           No          | Number of FM refinement iterations or passes per level                   |
| `-max_moves`                       |      `100` |    `60` |        **Yes**        | Maximum number of moves allowed during a refinement pass                 |
| `-early_stop_ratio`                |      `0.5` |   `0.5` |           No          | Early-stopping threshold within refinement                               |
| `-v_cycle_flag`                    |     `true` |  `true` |           No          | Enables V-cycle refinement                                               |
| `-max_num_vcycle`                  |        `3` |     `1` |        **Yes**        | Maximum number of V-cycle iterations                                     |
| `-thr_coarsen_hyperedge_size_skip` |      `200` |   `200` |           No          | Hyperedges larger than this threshold are skipped during coarsening      |
| `-thr_coarsen_vertices`            |       `10` |    `10` |           No          | Stop coarsening when the number of vertices falls below this threshold   |
| `-thr_coarsen_hyperedges`          |       `50` |    `50` |           No          | Stop coarsening when the number of hyperedges falls below this threshold |
| `-coarsening_ratio`                |      `1.6` |   `1.6` |           No          | Target vertex-count reduction ratio between coarsening levels            |
| `-max_coarsen_iters`               |       `30` |    `30` |           No          | Maximum number of coarsening iterations or levels                        |

Only three non-seed configuration values differ from the tool defaults:

```text
-balance_constraint 0.1
-max_moves 100
-max_num_vcycle 3
```

The seed also differs across the 20 runs:

```text
-seed 0
-seed 1
...
-seed 19
```

All other explicitly passed options match the defaults at the pinned commit.

---

## Parameters Left at Their Defaults

The following parameters were not passed on the command line. They therefore used their default values from `partitionmgr.tcl` and the associated TritonPart implementation.

| Parameter                     | Effective value | Description                                                            |
| ----------------------------- | --------------- | ---------------------------------------------------------------------- |
| `-base_balance`               | `{0.5, 0.5}`    | Target fraction assigned to each of the two partitions                 |
| `-scale_factor`               | `{1.0, 1.0}`    | Per-partition scaling factors                                          |
| `-vertex_dimension`           | `1`             | Number of vertex-weight dimensions                                     |
| `-hyperedge_dimension`        | `1`             | Number of hyperedge-weight dimensions                                  |
| `-placement_dimension`        | `0`             | Placement optimization disabled                                        |
| `-e_wt_factors`               | `{1.0}`         | Hyperedge-weight objective factors                                     |
| `-v_wt_factors`               | `{1.0}`         | Vertex-weight factors                                                  |
| `-placement_wt_factors`       | `{1.0}`         | Placement-weight factors; inactive because placement dimension is zero |
| `-adj_diff_ratio`             | `0.0001`        | Adjacency-difference threshold                                         |
| `-min_num_vertices_each_part` | `4`             | Minimum number of vertices required in each partition                  |
| `-total_corking_passes`       | `25`            | Number of corking passes used by the refinement procedure              |
| `-num_vertices_threshold_ilp` | `50`            | Vertex-count threshold for invoking ILP-based optimization             |
| `-global_net_threshold`       | `1000`          | Threshold used to identify or handle global nets                       |
| `-fixed_file`                 | Empty           | No fixed-vertex constraints                                            |
| `-community_file`             | Empty           | No community constraints                                               |
| `-group_file`                 | Empty           | No group constraints                                                   |
| `-placement_file`             | Empty           | No placement or embedding constraints                                  |

For two-way partitioning, the Tcl default for `-base_balance` is normalized to:

```tcl
{0.5 0.5}
```

That is, the target partition weight is 50% for each block.

Timing-driven partitioning was disabled. Each run produced the corresponding TritonPart warning:

```text
WARNING PAR-0003
```

Therefore, the reported objective is the hypergraph cut cost rather than a timing-weighted partitioning objective.

---

## Important Balance-Constraint Semantics

### `-balance_constraint` is expressed in percent

The value passed to:

```tcl
-balance_constraint
```

is interpreted as a percentage, not as a fractional ratio.

In `Hypergraph::GetUpperVertexBalance`, the upper balance bound is computed using an expression equivalent to:

```text
(base_balance[i] + balance_constraint × 0.01)
    × total_vertex_weight
```

For these experiments:

```text
base_balance[i] = 0.5
balance_constraint = 0.1
```

Therefore:

```text
upper bound
= (0.5 + 0.1 × 0.01) × total weight
= (0.5 + 0.001) × total weight
= 0.501 × total weight
```

Each partition was consequently limited to approximately 50.1% of the total vertex weight.

For a two-way partition, this corresponds to partition weights near:

```text
49.9% / 50.1%
```

rather than:

```text
40% / 60%
```

### Effective imbalance used in the experiments

The command:

```tcl
-balance_constraint 0.1
```

means an allowed imbalance of:

```text
0.1%
```

It does **not** mean an allowed imbalance of 10%.

This interpretation is confirmed by the experiment logs:

* All 440 runs satisfy the balance constraint.

* Observed block-balance fractions fall within:

  ```text
  [0.49900, 0.50100]
  ```

* Every run reports:

  ```text
  Satisfy the balance constraint : true
  ```

The cuts reported in `summary_best_of_20.csv` therefore correspond to nearly exact bisections.

### Running with a true 10% imbalance

To allow an actual 10% imbalance, use:

```tcl
-balance_constraint 10
```

For a two-way partition with equal base balance, this permits an upper block size of approximately 60%:

```text
0.5 + 10 × 0.01 = 0.6
```

The complementary block may consequently contain approximately 40%, giving a feasible range near:

```text
40% / 60%
```

Allowing a larger imbalance generally enlarges the feasible solution space and may produce smaller cutsizes. Such results must not be compared directly with the nearly exact-balance results reported here without clearly identifying the different balance constraint.

---

## Best-of-20 Titan23 Results

For every Titan23 design, TritonPart was executed with seeds `0` through `19`.

The result reported in `summary_best_of_20.csv` is:

```text
minimum Cutcost over the 20 seeds
```

The corresponding partition and log files are:

```text
partitions/<design>_k2_ub0.1_seed<best_seed>.part
logs/<design>_k2_ub0.1_seed<best_seed>.log
```

| Design          | Best cutsize | Best seed | Average cutsize over 20 seeds |
| --------------- | -----------: | --------: | ----------------------------: |
| `LU230`         |         3637 |        16 |                       4115.70 |
| `LU_Network`    |          730 |        18 |                       1009.15 |
| `SLAM_spheric`  |         1134 |         0 |                       1304.00 |
| `bitcoin_miner` |         1491 |         1 |                       1561.30 |
| `bitonic_mesh`  |          592 |         4 |                        671.20 |
| `cholesky_bdti` |         1246 |        13 |                       1345.40 |
| `cholesky_mc`   |          333 |         1 |                        472.00 |
| `dart`          |          839 |        15 |                        910.60 |
| `denoise`       |          467 |        15 |                        570.05 |
| `des90`         |          421 |         7 |                        494.45 |
| `directrf`      |          683 |        11 |                        871.00 |
| `gsm_switch`    |         1783 |         8 |                       2100.60 |
| `mes_noc`       |          827 |        12 |                       1133.10 |
| `minres`        |          271 |         3 |                        313.10 |
| `neuron`        |          248 |        14 |                        305.30 |
| `openCV`        |          591 |         8 |                        654.55 |
| `segmentation`  |          194 |        12 |                        289.55 |
| `sparcT1_chip2` |         1094 |         1 |                       1667.45 |
| `sparcT1_core`  |         1112 |         1 |                       1265.60 |
| `sparcT2_core`  |         1248 |         7 |                       1595.90 |
| `stap_qrd`      |          450 |        13 |                        527.35 |
| `stereo_vision` |          182 |         2 |                        188.10 |

The best seed is part of the reproducibility record. Re-running only the benchmark name without the corresponding seed does not necessarily reproduce the reported cutsize.

---

# Building TritonPart

## Prerequisites

TritonPart is implemented within OpenROAD. Before building, install the dependencies required by the OpenROAD application.

Refer to the official OpenROAD build documentation:

* [OpenROAD repository](https://github.com/The-OpenROAD-Project/OpenROAD)
* [OpenROAD build documentation](https://openroad.readthedocs.io/)

TritonPart may also use the following optimization packages:

* [Google OR-Tools](https://developers.google.com/optimization/install)
* IBM CPLEX

Unload Anaconda or remove conflicting Anaconda paths from the environment before configuring the build.

For example:

```bash
conda deactivate
```

Repeat the command if multiple Conda environments are nested.

It may also be necessary to remove Conda-provided compilers and libraries from `PATH`, `LD_LIBRARY_PATH`, or `CMAKE_PREFIX_PATH`.

---

## Build Instructions

From the TritonPart repository root:

```bash
mkdir -p build
cd build
cmake ../OpenROAD/
make -j
```

The resulting OpenROAD executable is typically written to:

```text
build/bin/openroad
```

Verify the build before running experiments:

```bash
./bin/openroad -version
```

For the reported Titan23 results, the expected version string is:

```text
OpenROAD v2.0-27290-g5e8f15c075
```

---

# Basic Hypergraph Partitioning

TritonPart accepts hypergraphs in hMETIS format.

A basic three-way partitioning command is:

```tcl
triton_part_hypergraph \
    -hypergraph_file des90.hgr \
    -num_parts 3 \
    -balance_constraint 2 \
    -seed 2
```

Here:

* `des90.hgr` is the input hypergraph.
* `-num_parts 3` requests three partitions.
* `-balance_constraint 2` specifies a 2% balance tolerance.
* `-seed 2` sets the random seed.

To run the included example:

```bash
cd test
source run.sh
```

Alternatively, execute the Tcl script directly through OpenROAD:

```bash
../build/bin/openroad -exit run.tcl
```

---

# Fixed-Vertex Constraints

TritonPart supports fixed vertices through a fixed-vertex constraint file.

Example:

```tcl
triton_part_hypergraph \
    -hypergraph_file sample.hgr \
    -num_parts 2 \
    -balance_constraint 20 \
    -seed 2 \
    -fixed_file sample.fix
```

To run the example:

```bash
cd test
source run_fix.sh
```

The fixed-vertex file assigns selected vertices to predetermined partitions. These assignments are treated as hard constraints during partitioning.

---

# Embedding and Placement Constraints

TritonPart supports placement-aware or embedding-aware partitioning.

Example:

```tcl
triton_part_hypergraph \
    -hypergraph_file sample.hgr \
    -num_parts 2 \
    -balance_constraint 20 \
    -seed 2 \
    -placement_file sample.emb \
    -placement_wt_factors {1.0 1.0} \
    -placement_dimension 2
```

To run the example:

```bash
cd test
source run_emb.sh
```

In this example:

* `sample.emb` supplies a two-dimensional embedding.
* `-placement_dimension 2` enables two placement dimensions.
* `-placement_wt_factors {1.0 1.0}` assigns equal importance to both dimensions.

---

# Running the Regression Suite

Regression parameters can be configured in:

```text
regression/regression.py
```

The relevant parameter block begins near the section corresponding to lines 131–148 in the referenced regression script.

Run the regression suite using:

```bash
cd regression
python3 regression.py | tee summary.log
```

The `tee` command writes the complete output to:

```text
summary.log
```

When comparing regression results, record:

* TritonPart/OpenROAD commit
* Build configuration
* Solver versions
* Container or operating-system environment
* Hypergraph benchmark revision
* Number of partitions
* Balance constraint
* Random seeds
* Number of initial solutions
* Coarsening parameters
* Refinement parameters
* V-cycle configuration
* Result-selection procedure

---

# Reporting TritonPart Results

A TritonPart cutsize should not be reported without sufficient configuration information to reproduce it.

At minimum, report:

```text
Repository or source tree:
OpenROAD/TritonPart commit:
OpenROAD version string:
Benchmark suite and revision:
Hypergraph format:
Number of partitions:
Balance constraint:
Balance-constraint interpretation:
Random seed or seed range:
Number of runs:
Initial-solution parameters:
Coarsening parameters:
Refinement parameters:
V-cycle parameters:
Constraints enabled:
Objective or cut metric:
Best-result selection rule:
Container or system environment:
```

For example:

```text
Tool: TritonPart in OpenROAD
Commit: 5e8f15c075ddd9d45d750977ad9a856f5c1c2569
OpenROAD version: OpenROAD v2.0-27290-g5e8f15c075
Benchmark: Titan23
Problem: 2-way partitioning
Balance constraint: 0.1%
Seeds: 0–19
Runs per design: 20
Reported result: Minimum Cutcost over 20 seeds
Initial solutions: 50
Best initial solutions: 10
Coarsening solutions: 3
Refinement iterations: 10
Maximum moves: 100
V-cycle enabled: true
Maximum V-cycles: 3
Timing driven: false
Placement constraints: none
Fixed vertices: none
Community constraints: none
Group constraints: none
```

This information is necessary because cutsize can depend on:

* Random seed
* Balance tolerance
* Number of initial solutions
* Coarsening strategy
* Refinement effort
* Number of V-cycles
* Constraint files
* Source-code revision
* Default-parameter changes
* Benchmark preprocessing
* Result-selection methodology

For randomized partitioners, distinguish clearly between:

* A single-seed result
* The average over multiple seeds
* The median over multiple seeds
* The best result over multiple seeds

The Titan23 cuts documented in this README are specifically:

```text
best-of-20 results under an effective 0.1% balance constraint
```
