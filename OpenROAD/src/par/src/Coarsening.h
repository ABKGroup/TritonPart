///////////////////////////////////////////////////////////////////////////
//
// BSD 3-Clause License
//
// Copyright (c) 2020, The Regents of the University of California
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
///////////////////////////////////////////////////////////////////////////////
#pragma once
//
// This file contains the classes for Coarsening Phase
// Coarsening Phase is to generate a sequence of coarser hypergraph
// We define Coarsening as an operator class.
// It will accept a HGraph (std::shared_ptr<TPHypergraph>) as input
// and return a sequence of coarser hypergraphs 
//
#include "TPHypergraph.h"
#include "utl/Logger.h"

namespace par {

class Coarsening
{
 public:
  explicit Coarsening() {  }  // Use the default value for each parameter 

  // Constructor
  Coarsening(const std::vector<float> &e_wt_factors,
             const std::vector<float> &v_wt_factors,
             const std::vector<float> &p_wt_factors,
             float timing_factor,
             int path_traverse_step,
             const std::vector<float> &max_vertex_weights,
             int global_net_threshold, 
             int smallest_v_size_cgraph,
             int smallest_e_size_cgraph,
             float coarsening_ratio,
             int max_coarsen_iters, 
             int seed, 
             utl::Logger* logger)
  {
    e_wt_factors_ = e_wt_factors;
    v_wt_factors_ = v_wt_factors;
    p_wt_factors_ = p_wt_factors;
    timing_factor_ = timing_factor;
    path_traverse_step_ = path_traverse_step;
    max_vertex_weight_ = max_vertex_weights;
    global_net_threshold_ = global_net_threshold;
    smallest_v_size_cgraph_ = smallest_v_size_cgraph;
    smallest_e_size_cgraph_ = smallest_e_size_cgraph;
    coarsening_ratio_ = coarsening_ratio;
    max_coarsen_iters_ = max_coarsen_iters;
    seed_ = seed;
    logger_ = logger;
  }
   
  // Copy constructor
  Coarsening(const Coarsening& coarsening)
  {
    e_wt_factors_ = coarsening.e_wt_factors_;
    v_wt_factors_ = coarsening.v_wt_factors_;
    p_wt_factors_ = coarsening.p_wt_factors_;
    timing_factor_ = coarsening.timing_factor_;
    path_traverse_step_ = coarsening.path_traverse_step_;
    max_vertex_weight_ = coarsening.max_vertex_weight_;
    global_net_threshold_ = coarsening.global_net_threshold_;
    smallest_v_size_cgraph_ = coarsening.smallest_v_size_cgraph_;
    smallest_e_size_cgraph_ = coarsening.smallest_e_size_cgraph_;
    coarsening_ratio_ = coarsening.coarsening_ratio_;
    max_coarsen_iters_ = coarsening.max_coarsen_iters_;
    seed_ = coarsening.seed_;
    logger_ = coarsening.logger_;
  }

  // destructor
  ~Coarsening() = default;
  
  // Coarse a hypergraph into a sequence of coarser hypergraphs
  // using first-choice method
  // Note that here the hgraph CANNOT be CONST HGraph
  // because we need to update vertex_c_attr_ attribute of hgraph
  // during coarsening
  std::vector<HGraph> LazyFirstChoice(HGraph hgraph);
  
  // Coarse a hypergraph with first-choice method
  // Aggregate vertices into clusters.
  // Aggregate takes a hypergraph hgraph as input,
  // then return a coarser hypergrapg hgraph_c,
  // where each vertex in hgraph_c is a set of vertices in hgraph
  // Note that here the hgraph CANNOT be CONST HGraph
  // because we need to update vertex_c_attr_ attribute of hgraph
  // during coarsening
  HGraph Aggregate(HGraph hgraph);

 private:
  // e_wt_factors_: the factor for each dimension of hyperedge weights
  std::vector<float> e_wt_factors_;
  // v_wt_factors_ : the factor for each dimension of vertex weights
  std::vector<float> v_wt_factors_;
  // p_wt_factors_ : the factor for each dimension of placement attributes
  std::vector<float> p_wt_factors_;
  // timing factor : the factor for timing path
  float timing_factor_ = 0.0;
  int path_traverse_step_ = 2;  // start with a node in timing path and consider its 
                                // path_traverse_step_ neighbors

  // max_vertex_weight: the maximum allowed weight for vertex 
  std::vector<float> max_vertex_weight_;
  // global_net_threshold : the threshold for identifying a hyperedge is global
  int global_net_threshold_ = 5000;
  // match_global_net_threshold : the threshold for net when calculating neighbors
  int match_global_net_threshold_ = 500; 
  // smallest_v_size_cgraph : the minimum allowed number of vertices in coarsest hypergraph
  int smallest_v_size_cgraph_= 500;
  // smallest_e_size_cgraph : the minimum allowed number of hyperedges in coarsest hypergraph
  int smallest_e_size_cgraph_= 500;
  // coarsening_ratio : maxinum allowed ratio for number of vertices in adjacent hypergraph
  float coarsening_ratio_ = 1.7;
  // max_coarsen_iters : maximum coarsening iterations
  int max_coarsen_iters_ = 7;
  // if the adjacent graphs are similar (the difference of vertices), then stop coarsening
  float adj_diff_ratio_ = 0.0001; 
  

  int seed_ = 0;  // seed for random shuffling
  utl::Logger* logger_ = nullptr;

  // private functions

  // First-choice based node matching : map each vertex to its cluster
  // 1) : Fill vertex_c_attr which maps the vertex to its corresponding cluster
  //      vertex_c_attr has hgraph->num_vertices_ elements
  // 2) : Fill vertex_weights_c which stores the weight for each cluster
  // 3) : Fill placement_attr_c which stores the placement attribute for each vertex 
  //      (may be empty)
  // 4) : Fill community_attr_c and fixed_attr_
  // 5) : Remaining elements for hgraphc (coarser hypergraph) will be determined by Aggregate
  // In our implementation, fixed vertices does not involve in coarsening, which means 
  // fixed vertices will be not touched during coarsening. Before we start coarsening, 
  // the fixed vertices has been merged into clusters. 
  // Thus, each block can have only one fixed vertex (cluster)
  void VertexMatching(const HGraph       hgraph,
                      std::vector<int>   &vertex_c_attr,
                      std::vector<std::vector<float> > &vertex_weights_c,
                      std::vector<int>   &community_attr_c,
                      std::vector<int>   &fixed_attr_c,
                      std::vector<std::vector<float> > &placement_attr_c);

  // friend class
  friend class MultiLevelHierarchy;
};


// nickname for shared pointer of Coarsening
using CoarseningPtr = std::shared_ptr<Coarsening>;

}  // namespace par

