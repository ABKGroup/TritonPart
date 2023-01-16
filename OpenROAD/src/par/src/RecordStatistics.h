///////////////////////////////////////////////////////////////////////////
//
// BSD 3-Clause License
//
// Copyright (c) 2022, The Regents of the University of California
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

namespace par {

class Stats
{
 public:
  Stats() = default;
  Stats(const float avg_fanin,
        const float avg_fanout,
        const float avg_hyperedge_size,
        const float avg_terminals)
      : avg_fanin_(avg_fanin),
        avg_fanout_(avg_fanout),
        avg_hyperedge_size_(avg_hyperedge_size),
        avg_terminals_(avg_terminals)
  {
  }
  Stats(const Stats&) = default;
  Stats& operator=(const Stats&) = default;
  Stats& operator=(Stats&&) = default;
  ~Stats() = default;
  float GetAvgFanIn() const { return avg_fanin_; }
  float GetAvgFanOut() const { return avg_fanout_; }
  float GetAvgHyperedgeSize() const { return avg_hyperedge_size_; }
  float GetAvgTerminals() const { return avg_terminals_; }

 private:
  float avg_fanin_;
  float avg_fanout_;
  float avg_hyperedge_size_;
  float avg_terminals_;
};

class Obfuscator
{
 public:
  Obfuscator()
      : global_net_threshold_(5000),
        contract_global_net_threshold_(500),
        contraction_factor_(0.5),
        seed_(0),
        avg_fanin_(0.0),
        avg_fanout_(0.0),
        avg_terminals_(0.0),
        timing_factor_(0.0),
        path_traverse_step_(2)
  {
  }
  Obfuscator(const Obfuscator&) = default;
  Obfuscator(Obfuscator&&) = default;
  Obfuscator& operator=(const Obfuscator&) = default;
  Obfuscator& operator=(Obfuscator&&) = default;
  Obfuscator(const std::vector<float>& e_wt_factors,
             const std::vector<float>& v_wt_factors,
             const std::vector<float>& p_wt_factors,
             float timing_factor,
             int path_traverse_step,
             const std::vector<float>& max_vertex_weights,
             int seed,
             utl::Logger* logger)
      : avg_fanin_(0.0),
        avg_fanout_(0.0),
        avg_terminals_(0.0),
        timing_factor_(0.0),
        path_traverse_step_(2),
        logger_(logger)
  {
    e_wt_factors_ = e_wt_factors;
    v_wt_factors_ = v_wt_factors;
    p_wt_factors_ = p_wt_factors;
    timing_factor_ = timing_factor;
    path_traverse_step_ = path_traverse_step;
    max_vertex_weight_ = max_vertex_weights;
    seed_ = seed;
  }
  ~Obfuscator() = default;
  std::shared_ptr<Stats> RecordStatistics(const HGraph hypergraph)
  {
    InitFanIns(hypergraph);
    InitFanOuts(hypergraph);
    InitTerminals(hypergraph);
    GenerateFanInsOuts(hypergraph);
    GenerateTerminals(hypergraph);
    FindAverageFanIn();
    FindAverageFanOut();
    FindAverageHyperedgeSize(hypergraph);
    FindAverageTerminals();
    return std::make_shared<Stats>(GetAvgFanIn(),
                                   GetAvgFanOut(),
                                   GetAvgHyperedgeSize(),
                                   GetAvgTerminals());
  }
  void SetContractionFactor(const float contraction_factor)
  {
    contraction_factor_ = contraction_factor;
  }
  void SetGlobalNetThreshold(const int global_net_threshold)
  {
    global_net_threshold_ = global_net_threshold;
  }
  void SetContractGlobalNetThreshold(const int contract_global_net_threshold)
  {
    contract_global_net_threshold_ = contract_global_net_threshold;
  }
  void InitFanIns(const HGraph hypergraph)
  {
    fanins_.resize(hypergraph->num_vertices_);
    std::fill(fanins_.begin(), fanins_.end(), 0);
  }
  void InitFanOuts(const HGraph hypergraph)
  {
    fanouts_.resize(hypergraph->num_vertices_);
    std::fill(fanouts_.begin(), fanouts_.end(), 0);
  }
  void InitTerminals(const HGraph hypergraph)
  {
    terminals_.resize(hypergraph->num_vertices_);
    std::fill(terminals_.begin(), terminals_.end(), 0);
  }
  void SetRandSeed(const int seed) { seed_ = seed; }
  float GetAvgFanIn() const { return avg_fanin_; }
  float GetAvgFanOut() const { return avg_fanout_; }
  float GetAvgTerminals() const { return avg_terminals_; }
  float GetAvgHyperedgeSize() const { return avg_hyperedge_size_; }
  void Obfuscate(HGraph hypergraph,
                 const float contraction_factor,
                 const int global_net_threshold,
                 const int contract_global_net_threshold,
                 const int seed);

 private:
  void GenerateFanInsOuts(const HGraph hypergraph);
  void GenerateTerminals(const HGraph hypergraph);
  void FindAverageFanIn()
  {
    avg_fanin_
        = std::accumulate(fanins_.begin(), fanins_.end(), 0.0) / fanins_.size();
  }
  void FindAverageFanOut()
  {
    avg_fanout_ = std::accumulate(fanouts_.begin(), fanouts_.end(), 0.0)
                  / fanouts_.size();
  }
  void FindAverageTerminals()
  {
    avg_terminals_ = std::accumulate(terminals_.begin(), terminals_.end(), 0.0)
                     / terminals_.size();
  }
  void FindAverageHyperedgeSize(const HGraph hypergraph);
  void QuickMatching(const HGraph hgraph,
                     std::vector<int>& vertex_c_attr,
                     std::vector<std::vector<float>>& vertex_weights_c,
                     std::vector<int>& community_attr_c,
                     std::vector<int>& fixed_attr_c,
                     std::vector<std::vector<float>>& placement_attr_c);
  HGraph OneLevelContraction(HGraph hypergraph);
  // member variables for recording hypergraph statistics
  int global_net_threshold_;
  int contract_global_net_threshold_;
  int seed_;
  float contraction_factor_;
  float avg_fanout_;
  float avg_fanin_;
  float avg_hyperedge_size_;
  float avg_terminals_;
  std::vector<int> fanins_;
  std::vector<int> fanouts_;
  std::vector<int> terminals_;
  // coarsening related member variables
  std::vector<float> e_wt_factors_;
  std::vector<float> v_wt_factors_;
  std::vector<float> p_wt_factors_;
  float timing_factor_;
  int path_traverse_step_;
  // max_vertex_weight: the maximum allowed weight for vertex
  std::vector<float> max_vertex_weight_;
  utl::Logger* logger_ = nullptr;
};

using ObfuscatorPtr = std::shared_ptr<Obfuscator>;
}  // namespace par