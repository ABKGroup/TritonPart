# TritonPart : A Constraints Driven General Hypergraph Partitioner

TritonPart is a constraints-driven general hypergraph partitioners.
Currently our repo is under active development.


## Table of Content
  - [OpenROAD](https://github.com/ABKGroup/TritonPart/tree/main/OpenROAD) includes the source for TritonPart. The codes are available in OpenROAD/src/par.
  - [titan23_benchmark](https://github.com/ABKGroup/TritonPart/tree/main/titan23_benchmark) provides the titan23 benchmark suites.
  - [regression](https://github.com/ABKGroup/TritonPart/tree/main/regression) provides the scripts for benchmarking TritonPart.
  
 
## Features
- Start of the art multiple-constraints driven partitioning “multi-tool”
- Optimizes cost function based on user requirement
- Permissive open-source license
- Solves multi-way partitioning with following features:
  - Multidimensional weights on vertices and hyperedges
  - Multilevel coarsening and refinement framework
  - Fixed vertices
  - Timing-driven partitioning framework (under development)
  - Community constraints: Groups of vertices need to be in same partition
  
 
## How to build TritonPart
The first step, independent of the build method, is to download the repository:
- We implement our TritonPart based on OpenROAD app.  Please check the dependencies of [OpenROAD app](https://github.com/The-OpenROAD-Project/OpenROAD.git) first.  
- We use Google OR-Tools as our ILP solver.  Please install Google OR-Tools following the [instructions](https://developers.google.com/optimization/install).
- Run following scripts
``` shell
mkdir build
cd build
cmake ../OpenROAD/
make
```

### How to partition a hypergraph like hMETIS
``` shell
triton_part_hypergraph -hypergraph_file des90.hgr -num_parts 5 -balance_constraint 2 -seed 2
```


## How to run regression test
We also users to specify their parameters by modifying the [regression.py](https://github.com/ABKGroup/TritonPart/blob/main/regression/regression.py) from [Line131](https://github.com/ABKGroup/TritonPart/blob/68e516145a7090c3b9bae7ac9bf2464e58758b69/regression/regression.py#L131) to 
[Line148](https://github.com/ABKGroup/TritonPart/blob/68e516145a7090c3b9bae7ac9bf2464e58758b69/regression/regression.py#L148).
You can run the regression scripts as following:
``` shell
python regression.py | tee summary.log
```

