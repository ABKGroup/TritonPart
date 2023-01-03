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

#include <set>

#include "TPHypergraph.h"
#include "Utilities.h"
#include "utl/Logger.h"

namespace par {
template <typename T>
using matrix = std::vector<std::vector<T>>;
using partitiontoken = std::pair<float, std::vector<std::vector<float>>>;
class IlpGraph
{
 public:
  IlpGraph() = default;
  IlpGraph(const IlpGraph&) = default;
  IlpGraph(IlpGraph&&) = default;
  IlpGraph& operator=(const IlpGraph&) = default;
  IlpGraph& operator=(IlpGraph&&) = default;
  ~IlpGraph() = default;
  IlpGraph(const int vertex_dimensions,
           const int hyperedge_dimensions,
           const bool fixed_flag,
           std::vector<int> fixed,
           const matrix<int>& hyperedges,
           const matrix<float>& vertex_weights,
           const matrix<float>& hyperedge_weights);
  inline const int GetVertexDimensions() const { return vertex_dimensions_; }
  inline const int GetHyperedgeDimensions() const
  {
    return hyperedge_dimensions_;
  }
  inline const int GetNumVertices() const { return num_vertices_; }
  inline const int GetNumHyperedges() const { return num_hyperedges_; }
  inline std::pair<int, int> GetEdgeIndices(const int he) const
  {
    return std::make_pair(eptr_[he], eptr_[he + 1]);
  }
  inline std::vector<float> const& GetVertexWeight(const int& v) const
  {
    return vertex_weights_[v];
  }
  inline std::vector<float> const& GetHyperedgeWeight(const int& e) const
  {
    return hyperedge_weights_[e];
  }
  inline bool CheckFixedStatus(const int v) const { return fixed_[v] > -1; }
  inline int GetFixedPart(const int v) const { return fixed_[v]; }
  inline bool CheckFixedFlag() const { return fixed_flag_; }
  inline std::vector<float> GetTotalVertexWeights()
  {
    std::vector<float> total_wt(GetVertexDimensions(), 0.0);
    for (auto vWt : vertex_weights_) {
      total_wt = total_wt + vWt;
    }
    return total_wt;
  }
  std::vector<int> eind_;

 private:
  int hyperedge_dimensions_;
  int vertex_dimensions_;
  int num_vertices_;
  int num_hyperedges_;
  bool fixed_flag_;
  std::vector<int> fixed_;
  std::vector<int> eptr_;
  matrix<float> vertex_weights_;
  matrix<float> hyperedge_weights_;
};

class PartitionPair
{
 public:
  PartitionPair() = default;
  PartitionPair(const PartitionPair&) = default;
  PartitionPair(PartitionPair&&) = default;
  PartitionPair& operator=(const PartitionPair&) = default;
  PartitionPair& operator=(PartitionPair&&) = default;
  PartitionPair(int partition_x, int partition_y, float score)
      : partition_x_(partition_x), partition_y_(partition_y), score_(score)
  {
  }
  inline int GetPartitionX() const { return partition_x_; }
  inline int GetPartitionY() const { return partition_y_; }
  inline float GetScore() const { return score_; }
  inline void SetPartitionX(int partition_x) { partition_x_ = partition_x; }
  inline void SetPartitionY(int partition_y) { partition_y_ = partition_y; }
  inline void SetScore(float score) { score_ = score; }

 private:
  int partition_x_;
  int partition_y_;
  float score_;
};

auto comp = [](const PartitionPair& p1, const PartitionPair& p2) {
  return p1.GetScore() < p2.GetScore();
};

using partitionqueue = std::
    priority_queue<PartitionPair, std::vector<PartitionPair>, decltype(comp)>;

class IlpRefiner
{
 public:
  IlpRefiner() = default;
  IlpRefiner(int num_parts,
             int seed,
             int wavefront,
             std::vector<float> e_wt_factors,
             float path_wt_factor,
             float snaking_wt_factor,
             int max_passes)
  {
    num_parts_ = num_parts;
    seed_ = seed;
    wavefront_ = wavefront;
    e_wt_factors_ = e_wt_factors;
    path_wt_factor_ = path_wt_factor;
    snaking_wt_factor_ = snaking_wt_factor;
    max_passes_ = max_passes;
  }

  IlpRefiner(const IlpRefiner&) = default;
  IlpRefiner(IlpRefiner&&) = default;
  IlpRefiner& operator=(const IlpRefiner&) = default;
  IlpRefiner& operator=(IlpRefiner&&) = default;
  ~IlpRefiner() = default;
  int GetHyperedgeThreshold() const { return he_size_threshold_; }
  void Refine(HGraph hypergraph, std::vector<int>& partition);
  inline void SetHyperedgeThreshold(int he_size_threshold)
  {
    he_size_threshold_ = he_size_threshold;
  }
  inline void SetMaxBalance(matrix<float> max_balance)
  {
    max_block_balance_ = max_balance;
  }
  inline void SetCurrBalance(matrix<float> curr_balance)
  {
    curr_block_balance_ = curr_balance;
  }
  void SetClusterMap(const int total_ele) { cluster_map_.resize(total_ele); }
  inline void ResetClusterMap()
  {
    std::fill(cluster_map_.begin(), cluster_map_.end(), -1);
  }

 private:
  partitiontoken CutEvaluator(std::shared_ptr<IlpGraph> hypergraph,
                              std::vector<int>& solution,
                              bool print_flag);
  inline matrix<float> GetBlockBalance(std::shared_ptr<IlpGraph> hypergraph,
                                       std::vector<int>& partition);
  inline matrix<float> GetBlockBalance(const HGraph hgraph,
                                       std::vector<int>& partition);
  
  partitiontoken CutEvaluator(const HGraph hgraph,
                              std::vector<int>& solution,
                              bool print_flag);
  inline void Remap(std::vector<int>& partition,
                    std::vector<int>& refined_partition);
  inline void ResetWeights(const int& vertex_id,
                           const int& from_part,
                           const int& to_part)
  {
    curr_block_balance_[from_part]
        = curr_block_balance_[from_part] - vtx_weights_[vertex_id];
    curr_block_balance_[to_part]
        = curr_block_balance_[to_part] + vtx_weights_[vertex_id];
  }
  void OrderVertexSet(matrix<int> net_degs,
                      const std::vector<float>& cur_path_cost,
                      HGraph hypergraph,
                      std::vector<int>& boundary_vertices,
                      int from_pid,
                      int to_pid,
                      std::vector<int>& partition);
  matrix<float> GetMaxBlockBalance() { return max_block_balance_; }
  matrix<float> GetCurrBlockBalance() { return curr_block_balance_; }
  float FindPairScore(HGraph hypergraph,
                      std::vector<int>& partition,
                      int& pair_x,
                      int& pair_y);
  partitionqueue FindMaximalPairs(HGraph hypergraph,
                                  std::vector<int>& partition);
  void SolveIlpInstance(std::shared_ptr<IlpGraph> hypergraph,
                        std::vector<int>& partition,
                        const int& part_x,
                        const int& part_y);
  void DebugIlpInstance(const char* file_name);
  void FindNetDegs(HGraph hypergraph,
                   matrix<int>& net_degs,
                   std::vector<int>& partition);
  void UpdateNetDegs(HGraph hypergraph,
                     std::vector<int>& vertices,
                     matrix<int>& net_degs,
                     std::vector<int>& prev_partition,
                     std::vector<int>& new_partition);
  float CalculatePathCost(int path_id,
                          const HGraph hypergraph,
                          const std::vector<int>& partition,
                          int v = -1,
                          int to_pid = -1);
  float CalculateGain(int v,
                      int from_pid,
                      int to_pid,
                      const HGraph hypergraph,
                      const std::vector<int>& partition,
                      const std::vector<float>& cur_path_cost,
                      const std::vector<std::vector<int>>& net_degs);
  void GenerateCommunities(HGraph hypergraph,
                           const int& part_x,
                           const int& part_y,
                           std::vector<int>& ordered_set,
                           std::vector<int>& fixed_vertices,
                           matrix<float>& vertex_weights_coarsened,
                           std::vector<int>& partition);
  std::shared_ptr<IlpGraph> ContractHypergraph(
      HGraph hypergraph,
      matrix<float>& vertex_weights_c,
      std::vector<int>& fixed_vertices);
  bool CheckBalance(const int& vertex_id,
                    const int& from_part,
                    const int& to_part,
                    HGraph hypergraph);
  std::vector<int> FindBoundaryNodes(HGraph hypergraph,
                                     const int& part_x,
                                     const int& part_y,
                                     std::vector<int>& partition,
                                     matrix<int>& net_degs);
  std::priority_queue<PartitionPair> pairs_;
  std::vector<int> cluster_map_;
  matrix<float> curr_block_balance_;
  matrix<float> max_block_balance_;
  matrix<float> vtx_weights_;
  int num_parts_;
  int seed_;
  int wavefront_;
  int he_size_threshold_;
  std::vector<float> e_wt_factors_;
  float path_wt_factor_;
  float snaking_wt_factor_;
  int max_passes_;
  utl::Logger* logger_ = nullptr;
  std::vector<bool> boundary_;
};
using IlpRefinerPtr = std::shared_ptr<IlpRefiner>;
}  // namespace par
