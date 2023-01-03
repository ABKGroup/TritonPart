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

class Refiners
{
 public:
  Refiners() = default;
  Refiners(int refiner_choice,
           int num_parts,
           int seed,
           std::vector<float> e_wt_factor,
           float path_wt_factor,
           float snaking_wt_factor,
           int max_iters)
      : refiner_choice_(refiner_choice),
        curr_block_balance_(curr_block_balance),
        max_block_balance_(max_block_balance),
        num_parts_(num_parts),
        seed_(seed),
        he_size_threshold_(he_size_threshold),
        e_wt_factor_(e_wt_factor),
        path_wt_factor_(path_wt_factor),
        snaking_wt_factor_(snaking_wt_factor),
        max_iters_(max_iters)
  {
  }
  Refiners(const Refiners&) = default;
  Refiners(Refiners&&) = default;
  Refiners& operator=(const Refiners&) = default;
  Refiners& operator=(Refiners&&) = default;
  ~Refiners() = default;
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
  matrix<float> GetMaxBlockBalance() { return max_block_balance_; }
  matrix<float> GetCurrBlockBalance() { return curr_block_balance_; }
  inline void ResetWeights(const int& vertex_id,
                           const int& from_part,
                           const int& to_part)
  {
    curr_block_balance_[from_part]
        = curr_block_balance_[from_part] - vtx_weights_[vertex_id];
    curr_block_balance_[to_part]
        = curr_block_balance_[to_part] + vtx_weights_[vertex_id];
  }
  inline matrix<float> GetBlockBalance(std::shared_ptr<IlpGraph> hypergraph,
                                       std::vector<int>& partition);
  partitiontoken CutEvaluator(std::shared_ptr<IlpGraph> hypergraph,
                              std::vector<int>& solution,
                              bool print_flag);

 protected:
  bool CheckBalance(const int& vertex_id,
                    const int& from_part,
                    const int& to_part,
                    HGraph hypergraph);
  float CalculateGain(int v,
                      int from_pid,
                      int to_pid,
                      const HGraph hypergraph,
                      const std::vector<int>& partition,
                      const std::vector<float>& cur_path_cost,
                      const std::vector<std::vector<int>>& net_degs);
  float CalculatePathCost(int path_id,
                          const HGraph hypergraph,
                          const std::vector<int>& partition,
                          int v = -1,
                          int to_pid = -1);
  void FindNetDegs(HGraph hypergraph,
                   matrix<int>& net_degs,
                   std::vector<int>& partition);
  int refiner_choice_;
  matrix<float> curr_block_balance_;
  matrix<float> max_block_balance_;
  int num_parts_;
  int seed_;
  int he_size_threshold_;
  std::vector<float> e_wt_factor_;
  float path_wt_factor_;
  float snaking_wt_factor_;
  int max_iters_;
};

class IlpRefiner : public Refiners
{
 public:
  IlpRefiner() = default;
  IlpRefiner(int num_parts,
             int seed,
             int wavefront,
             std::vector<float> e_wt_factor,
             float path_wt_factor,
             float snaking_wt_factor,
             int max_iters)
      : Refiners(0,
                 num_parts,
                 seed,
                 e_wt_factor,
                 path_wt_factor,
                 snaking_wt_factor,
                 max_iters),
        wavefront_(wavefront)
  {
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
  void SetClusterMap(const int total_ele) { cluster_map_.resize(total_ele); }
  inline void ResetClusterMap()
  {
    std::fill(cluster_map_.begin(), cluster_map_.end(), -1);
  }

 private:
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
  void UpdateNetDegs(HGraph hypergraph,
                     std::vector<int>& vertices,
                     matrix<int>& net_degs,
                     std::vector<int>& prev_partition,
                     std::vector<int>& new_partition);
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
  std::vector<int> FindBoundaryNodes(HGraph hypergraph,
                                     const int& part_x,
                                     const int& part_y,
                                     std::vector<int>& partition,
                                     matrix<int>& net_degs);
  std::priority_queue<PartitionPair> pairs_;
  std::vector<int> cluster_map_;
  matrix<float> vtx_weights_;
  int wavefront_;
  std::vector<bool> boundary_;
};
using IlpRefinerPtr = std::shared_ptr<IlpRefiner>;

class PairwiseRefiner : public Refiners
{
 public:
  PairwiseRefiner() = default;
  PairwiseRefiner(int num_parts,
                  int seed,
                  const std::vector<float>& e_wt_factor,
                  float path_wt_factor,
                  float snaking_wt_factor,
                  float early_stop_ratio,
                  int max_iters)
      : Refiners(0,
                 num_parts,
                 seed,
                 e_wt_factor,
                 path_wt_factor,
                 snaking_wt_factor,
                 max_iters),
        early_stop_ratio_(early_stop_ratio)
  {
  }
  PairwiseRefiner(const PairwiseRefiner&) = default;
  PairwiseRefiner(PairwiseRefiner&&) = default;
  PairwiseRefiner& operator=(const PairwiseRefiner&) = default;
  PairwiseRefiner& operator=(PairwiseRefiner&&) = default;
  ~PairwiseRefiner() = default;
  float PairwisePass(const HGraph,
                std::vector<std::pair<int, int>>& partition_pairs,
                std::vector<int>& solution);
  void Refine(const HGraph hgraph,
                     std::vector<int>& solution);
  float CalculateSpan(const HGraph hgraph,
                         int& from_pid,
                         int& to_pid,
                         std::vector<int>& solution);
  static std::vector<int> FindBoundaryVertices(
      const HGraph hgraph,
      std::pair<int, int>& partition_pair,
      matrix<int>& net_degs);
  void SetThresholdMoves(const int moves) { max_num_moves_ = moves; }
  void InitVisitFlags(int total_ele)
  {
    visit_.resize(total_ele);
    std::fill(visit_.begin(), visit_.end(), false);
  }
  void ResetVisitFlags(int total_ele)
  {
    std::fill(visit_.begin(), visit_.end(), false);
  }
  void InitBoundaryFlags(int total_ele)
  {
    boundary_.resize(total_ele);
    std::fill(boundary_.begin(), boundary_.end(), false);
  }
  void ResetBoundaryFlags(int total_ele)
  {
    std::fill(boundary_.begin(), boundary_.end(), false);
  }
  void ResetVisited(int v) { visit_[v] = false; }
  void MarkVisited(int v) { visit_[v] = true; }
  bool GetVisitStatus(int v) const { return visit_[v]; }
  void MarkBoundary(int v) { boundary_[v] = true; }
  bool GetBoundaryStatus(int v) const { return boundary_[v]; }

 private:
  class VertexGain
  {
   public:
    VertexGain()
    {
      status_ = false;
      vertex_ = -1;
      source_part_ = -1;
      potential_move_ = -1;
    }
    VertexGain(int arg_vertex, float arg_gain)
        : vertex_(arg_vertex), gain_(arg_gain)
    {
      potential_move_ = -1;
      status_ = true;
    }
    VertexGain(int arg_vertex, float arg_gain, int part)
        : vertex_(arg_vertex), gain_(arg_gain), source_part_(part)
    {
      potential_move_ = -1;
      status_ = true;
    }
    VertexGain(int arg_vertex,
               float arg_gain,
               std::map<int, float> arg_path_cost)
        : vertex_(arg_vertex), gain_(arg_gain), path_cost_(arg_path_cost)
    {
      potential_move_ = -1;
      status_ = true;
    }
    VertexGain(int arg_vertex,
               float arg_gain,
               int part,
               std::map<int, float> arg_path_cost)
        : vertex_(arg_vertex),
          gain_(arg_gain),
          source_part_(part),
          path_cost_(arg_path_cost)
    {
      potential_move_ = -1;
      status_ = true;
    }
    VertexGain(const VertexGain&) = default;
    VertexGain& operator=(const VertexGain&) = default;
    VertexGain(VertexGain&&) = default;
    VertexGain& operator=(VertexGain&&) = default;
    bool operator<(const VertexGain& vertex_gain) const
    {
      if (gain_ > vertex_gain.gain_) {
        return true;
      } else if (gain_ == vertex_gain.gain_ && vertex_ < vertex_gain.vertex_) {
        return true;
      } else {
        return false;
      }
    }
    bool operator==(const VertexGain& vertex_gain) const
    {
      if (gain_ == vertex_gain.gain_ && path_cost_ == vertex_gain.path_cost_) {
        return true;
      } else {
        return false;
      }
    }
    int GetVertex() const { return vertex_; }
    int GetPotentialMove() const { return potential_move_; }
    float GetGain() const { return gain_; }
    bool GetStatus() const { return status_; }
    void SetVertex(int vertex) { vertex_ = vertex; }
    void SetGain(float gain) { gain_ = gain; }
    void SetActive() { status_ = true; }
    void SetDeactive() { status_ = false; }
    void SetPotentialMove(int move) { potential_move_ = move; }
    void SetPathCost(int path_id, float cost) { path_cost_[path_id] = cost; }
    float GetPathCost(int path_id) { return path_cost_[path_id]; }
    int GetTotalPaths() const { return path_cost_.size(); }
    int GetSourcePart() const { return source_part_; }

   private:
    int vertex_;
    int source_part_;
    int potential_move_;
    float gain_;
    bool status_;
    std::map<int, float>
        path_cost_;  // the updated path cost after moving vertex
  };

  // Priority queue implementation
  class PriorityQ
  {
   public:
    PriorityQ() = default;
    PriorityQ(int total_elements, HGraph hypergraph)
    {
      vertices_map_.resize(total_elements);
      std::fill(vertices_map_.begin(), vertices_map_.end(), -1);
      total_elements_ = 0;
      hypergraph_ = hypergraph;
      active_ = false;
    }
    PriorityQ(std::vector<std::shared_ptr<VertexGain>>& vertices,
              std::vector<int>& vertices_map)
        : vertices_(vertices), vertices_map_(vertices_map)
    {
    }
    PriorityQ(const PriorityQ&) = default;
    PriorityQ& operator=(const PriorityQ&) = default;
    PriorityQ(PriorityQ&&) = default;
    PriorityQ& operator=(PriorityQ&&) = default;
    ~PriorityQ() = default;
    inline void HeapifyUp(int index);
    void HeapifyDown(int index);
    void InsertIntoPQ(std::shared_ptr<VertexGain> element);
    std::shared_ptr<VertexGain> ExtractMax();
    void ChangePriority(int index, float priority);
    std::shared_ptr<VertexGain> GetMax() { return vertices_.front(); }
    void RemoveAt(int location);
    bool CheckIfEmpty() { return vertices_.empty(); }
    int GetTotalElements() const { return total_elements_; }
    int GetSizeOfPQ() const { return vertices_.size(); }
    int GetSizeOfMap() const { return vertices_map_.size(); }
    int GetLocationOfVertex(int v) const { return vertices_map_[v]; }
    std::shared_ptr<VertexGain> GetHeapVertex(int index)
    {
      return vertices_[index];
    }
    std::shared_ptr<VertexGain> GetHeapVertexFromId(int vertex_id)
    {
      assert(CheckIfVertexExists(vertex_id) == true);
      const int map_loc = GetLocationOfVertex(vertex_id);
      return vertices_[map_loc];
      // return &(vertices_[map_loc]);
    }
    void SetActive() { active_ = true; }
    void SetDeactive() { active_ = false; }
    bool GetStatus() const { return active_; }
    bool CheckIfVertexExists(int v)
    {
      return GetLocationOfVertex(v) > -1 ? true : false;
    }
    void Clear()
    {
      active_ = false;
      vertices_.clear();
      total_elements_ = 0;
      std::fill(vertices_map_.begin(), vertices_map_.end(), -1);
    }

   private:
    int Parent(int& element) { return floor((element - 1) / 2); }
    int LeftChild(int& element) { return ((2 * element) + 1); }
    int RightChild(int& element) { return ((2 * element) + 2); }

    bool active_;
    HGraph hypergraph_;
    std::vector<std::shared_ptr<VertexGain>> vertices_;
    std::vector<int> vertices_map_;
    int total_elements_;
  };

 public:
  // Alias for member classes PriorityQ, PQs and VertexGain
  using pq = PriorityQ;
  using vgain = VertexGain;
  using pqs = std::vector<std::shared_ptr<PriorityQ>>;

 private:
  void InitialGainsBetweenPairs(const HGraph hgraph,
                                const std::pair<int, int>& partition_pair,
                                const std::vector<int>& boundary_vertices,
                                const matrix<int>& net_degs,
                                const std::vector<int>& solution,
                                const std::vector<float>& cur_path_cost,
                                pqs& gain_buckets);
  std::vector<std::pair<int, int>> KPMfindPairs(
      const HGraph hgraph,
      std::vector<int>& solution,
      std::vector<float>& prev_scores);
  std::shared_ptr<VertexGain> PickVertexToMove(HGraph hgraph,
                                               pqs& gain_buckets);
  void AcceptMove(std::shared_ptr<VertexGain> vertex_to_move,
                  HGraph hgraph,
                  std::vector<VertexGain>& moves_trace,
                  float& total_gain,
                  float& total_delta_gain,
                  std::pair<int, int>& partition_pair,
                  std::vector<int>& solution,
                  std::vector<float>& paths_cost,
                  pqs& gain_buckets,
                  matrix<int>& net_degs);
  void RollBackMoves(std::vector<VertexGain>& trace,
                     matrix<int>& net_degs,
                     int& best_move,
                     float& total_delta_gain,
                     pqs& gain_buckets,
                     HGraph hgraph,
                     std::vector<float>& cur_path_cost,
                     std::vector<int>& solution);
  void UpdateNeighbors(const HGraph hgraph,
                       const std::pair<int, int>& partition_pair,
                       const std::set<int>& neighbors,
                       const std::vector<int>& solution,
                       const matrix<int>& net_degs,
                       pqs& gain_buckets);
  inline int GetConnectivity(const int& he, const matrix<int>& net_degs) const;
  inline bool CheckBoundaryVertex(const HGraph hgraph,
                                  const int& v,
                                  const std::pair<int, int> partition_pair,
                                  const matrix<int>& net_degs) const;
  int num_parts_;
  std::vector<float> e_wt_factor_;
  float path_wt_factor_;
  float snaking_wt_factor_;
  float early_stop_ratio_;
  float max_num_moves_;
  int max_num_fm_pass_;
  int max_stagnation_;
  int seed_;
  std::vector<bool> visit_;
  std::vector<bool> boundary_;
};
using kpm_heap = KPMRefinement::pq;
using vertex = KPMRefinement::vgain;
using kpm_heaps = KPMRefinement::pqs;
using PairwiseRefinerPtr = std::shared_ptr<PairwiseRefiner>;

class DirectRefiner : public Refiners
{
};
}  // namespace par