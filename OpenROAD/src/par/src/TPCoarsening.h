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
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <string>
#include <vector>

#include "TPHypergraph.h"
#include "db_sta/dbReadVerilog.hh"
#include "db_sta/dbSta.hh"
#include "odb/db.h"
#include "sta/Bfs.hh"
#include "sta/Graph.hh"
#include "sta/Liberty.hh"
#include "sta/Sta.hh"
#include "utl/Logger.h"

namespace par {
class Coarseners
{
 public:
  Coarseners(const Coarseners&) = delete;
  Coarseners(Coarseners&&) = delete;
  Coarseners& operator=(const Coarseners&) = delete;
  Coarseners& operator=(Coarseners&&) = delete;
  virtual ~Coarseners() = default;

 protected:
  // Constructor
  Coarsening(const std::vector<float>& e_wt_factors,
             const std::vector<float>& v_wt_factors,
             const std::vector<float>& p_wt_factors,
             float timing_factor,
             int path_traverse_step,
             const std::vector<float>& max_vertex_weights,
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
  Coarsening(const Coarsening& coarsening) {}
  Coarseners() = default;
  ~Coarseners() = default;
  virtual HGraph Aggregate(HGraph hgraph);
  void SetGlobalNetThreshold(int global_net_threshold)
  {
    global_net_threshold_ = global_net_threshold;
  }
  void SetMatchGlobalNetThreshold(int match_global_net_threshold)
  {
    match_global_net_threshold_ = match_global_net_threshold;
  }
  void SetSmallestVertexSizeCgraph(int smallest_v_size_cgraph)
  {
    smallest_v_size_cgraph_ = smallest_v_size_cgraph;
  }
  void SetSmallestEdgeSizeCgraph(int smallest_e_size_cgraph)
  {
    smallest_e_size_cgraph_ = smallest_e_size_cgraph;
  }
  void SetCoarseningRatio(float coarsening_ratio)
  {
    coarsening_ratio_ = coarsening_ratio;
  }
  void SetCoarsenIters(int max_iters) { max_coarsen_iters_ = max_iters; }
  void SetAdjDiffRatio(float adj_diff_ratio)
  {
    adj_diff_ratio_ = adj_diff_ratio;
  }
  int GetGlobalNetThreshold() const { return global_net_threshold_; }
  int GetMatchGlobalNetThreshold() const { return match_global_net_threshold_; }
  int GetSmallestVertexSizeCgraph() const { return smallest_v_size_cgraph_; }
  int GetSmallestEdgeSizeCgraph() const { return smallest_e_size_cgraph_; }
  void VertexMatching(const HGraph hgraph,
                      std::vector<int>& vertex_c_attr,
                      std::vector<std::vector<float>>& vertex_weights_c,
                      std::vector<int>& community_attr_c,
                      std::vector<int>& fixed_attr_c,
                      std::vector<std::vector<float>>& placement_attr_c);
  // Member variables
  std::vector<float> e_wt_factors_;
  std::vector<float> v_wt_factors_;
  std::vector<float> p_wt_factors_;
  float timing_factor_ = 0.0;
  int path_traverse_step_ = 2;
  std::vector<float> max_vertex_weight_;
  int global_net_threshold_
      = 5000;  // Threshold number of hyperedges to skip while clustering
  int match_global_net_threshold_
      = 500;  // Threshold number of hyperedges to skip while building coarse
              // hypergraph
  int smallest_v_size_cgraph_
      = 500;  // Threshold number of vertices at coarsest level
  int smallest_e_size_cgraph_
      = 500;  // Threshold number of hyperedges at coarsest level
  float coarsening_ratio_ = 1.7;  // Coarsening ratio of clustering
  int max_coarsen_iters_ = 7;     // Maximum number of coarsening iterations
  float adj_diff_ratio_ = 0.0001;
  int seed_ = 0;  // seed for random shuffling
  utl::Logger* logger_ = nullptr;
};

class FirstChoice : public Coarseners
{
 public:
  FirstChoice() = default;
  FirstChoice(FirstChoice&&) = default;
  FirstChoice& operator=(const FirstChoice&) = default;
  FirstChoice& operator=(FirstChoice&&) = default;
  FirstChoice(const std::vector<float>& e_wt_factors,
              const std::vector<float>& v_wt_factors,
              const std::vector<float>& p_wt_factors,
              float timing_factor,
              int path_traverse_step,
              const std::vector<float>& max_vertex_weights,
              int global_net_threshold,
              int smallest_v_size_cgraph,
              int smallest_e_size_cgraph,
              float coarsening_ratio,
              int max_coarsen_iters,
              int seed,
              utl::Logger* logger)
      : Coarseners(e_wt_factors,
                   v_wt_factors,
                   p_wt_factors,
                   timing_factor,
                   path_traverse_step,
                   max_vertex_weights,
                   global_net_threshold,
                   smallest_v_size_cgraph,
                   smallest_e_size_cgraph,
                   coarsening_ratio,
                   max_coarsen_iters,
                   seed,
                   logger)
  {
  }
  FirstChoice(const FirstChoice&)
  {
    e_wt_factors_ = FirstChoice.e_wt_factors_;
    v_wt_factors_ = FirstChoice.v_wt_factors_;
    p_wt_factors_ = FirstChoice.p_wt_factors_;
    timing_factor_ = FirstChoice.timing_factor_;
    path_traverse_step_ = FirstChoice.path_traverse_step_;
    max_vertex_weight_ = FirstChoice.max_vertex_weight_;
    global_net_threshold_ = FirstChoice.global_net_threshold_;
    smallest_v_size_cgraph_ = FirstChoice.smallest_v_size_cgraph_;
    smallest_e_size_cgraph_ = FirstChoice.smallest_e_size_cgraph_;
    coarsening_ratio_ = FirstChoice.coarsening_ratio_;
    max_coarsen_iters_ = FirstChoice.max_coarsen_iters_;
    seed_ = FirstChoice.seed_;
    logger_ = FirstChoice.logger_;
  }
  std::vector<HGraph> RunFirstChoice(HGraph hgraph);
  HGraph Aggregate(HGraph hgraph);
}
}  // namespace par
