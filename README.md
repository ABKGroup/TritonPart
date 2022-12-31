# TritonPart : A Constraints Driven General Hypergraph Partitioner

TritonPart is a constraints-driven general hypergraph partitioners.   

# Julia Requirement
Julia version 1.7.2 (must use this version)
You need following Julia packages:
- 



## Table of Content
  - [OpenROAD_TritonPart_src](https://github.com/ZhiangWang033/TritonPart_v1.0/tree/main/OpenROAD_TritonPart_src) includes the source for TritonPart. The codes are available in OpenROAD_TritonPart_src/src/par.
  - [hMETIS](https://github.com/ZhiangWang033/TritonPart_v1.0/tree/main/hMETIS) provides the hMETIS and khMETIS.
  - [titan23_benchmark](https://github.com/ZhiangWang033/TritonPart_v1.0/tree/main/titan23_benchmark) provides the titan23 benchmark suites.
  - [regression](https://github.com/ZhiangWang033/TritonPart_v1.0/tree/main/regression) provides the scripts for benchmarking TritonPart against hMETIS and khMETIS.
  
 
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
cmake ../OpenROAD_TritonPart_src/
make
```

## How to run regression test
We also users to specify their parameters by modifying the [regression.py](https://github.com/ZhiangWang033/TritonPart_v1.0/blob/main/regression/regression.py) from [Line131](https://github.com/ZhiangWang033/TritonPart_v1.0/blob/b5c37a90115ec6a856e77b37bdd495a9c5dc93df/regression/regression.py#L131) to [Line148](https://github.com/ZhiangWang033/TritonPart_v1.0/blob/b5c37a90115ec6a856e77b37bdd495a9c5dc93df/regression/regression.py#L148).
You can run the regression scripts as following:
``` shell
python regression.py | tee summary.log
```

## Flow
<p align="center">
<img src="./images/TritonPart_initial_flow.png" width= "1200"/>
</p>
<p align="center">
Figure 1.  TritonPart flow.  
</p>


## Results
<img src="./images/TritonPart_initial_result.png" width= "1200"/>
</p>
<p align="center">
Figure 2. Promising results (average of 10 runs) with current implementation.
</p>








