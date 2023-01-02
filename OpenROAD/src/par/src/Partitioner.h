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
//
// This file define the Paritioner class
// Paritioner is an operator class, which takes a hypergraph
// and perform partitioning based on multiple options.
// The Partitioner class implements following types of partitioning algorithm:
// 1) Greedy partitioning
// 2) Priority-queue FM-based partitioning
// 3) ILP-based partitioning
// 4) Iterative k-way partitioning
// 5) Direct k-way partitioning
// Currently we only consider the balance constraints on vertices
// TO DOs:
//    1) Handle the constraints on pins
//
//
#include <set>
#include "utl/Logger.h"
#include "TPHypergraph.h"
#include "Utilities.h"


namespace par {

// Define the partitioning algorithm
enum PartitionAlgorithm { Random, DirectGreedy, DirectFM, DirectILP,
                          RecursiveGreedy, RecursiveFM, RecursiveILP};

class Partitioners
{
  public:
    explicit Partitioners() {  }  // Use the default value for each parameter 
    Partitioners(int num_parts,
                 const std::vector<float> &e_wt_factors,
                 float path_wt_factor,
                 float snaking_wt_factor,
                 float early_stop_ratio,
                 int max_num_fm_pass,
                 int seed,
                 utl::Logger* logger)
    {
      num_parts_ = num_parts;
      e_wt_factors_ = e_wt_factors;
      path_wt_factor_ = path_wt_factor;
      snaking_wt_factor_ = snaking_wt_factor;
      early_stop_ratio_ = early_stop_ratio;
      max_num_fm_pass_ = max_num_fm_pass;
      seed_ = seed;
      logger_ = logger;
    }

    // destructor
    ~Partitioners() = default;

    // Golden Evaluator
    std::pair<float, std::vector<std::vector<float> > > 
          GoldenEvaluator(const HGraph hgraph, std::vector<int> &solution,
                          bool print_flag = true);
   
    // max_block_balance defines the balance of vertex weights in each block
    // for example, max_block_balance[0] is the maxinum allowed vertex weights for block 0
    // Random Partition
    void RandomPartition(const HGraph graph, 
                         const std::vector<std::vector<float> > &max_block_balance,
                         std::vector<int> &solution);

    // ILP based initial partitioning with Google OR-Tools
    void OptimalInitialPartition(const HGraph hgraph,
                                 const std::vector<std::vector<float> > &max_block_balance,
                                 std::vector<int> &solution);
    
    // ILP based speeded up initial partitioner with CPLEX
    void OptimalInitialPartitionCplexWarmStart(const HGraph hgraph,
                                 const std::vector<std::vector<float> > &max_block_balance,
                                 std::vector<int> &solution);

    // ILP based initial partitioning with CPLEX
    void OptimalInitialPartitionCplex(const HGraph hgraph,
                                 const std::vector<std::vector<float> > &max_block_balance,
                                 std::vector<int> &solution);

    // Balance a partition using FM 
    void BalancePartition(const HGraph hgraph, 
                          const std::vector<std::vector<float> > &max_block_balance,
                          std::vector<int> &solution);

    void DirectKWayFMWithImb(const HGraph hgraph,
                      const std::vector<std::vector<float> > &max_block_balance,
                      std::vector<int> &solution);

    // Direct Priority-queue FM-based partitioning
    void DirectKWayFM(const HGraph hgraph,
                      const std::vector<std::vector<float> > &max_block_balance,
                      std::vector<int> &solution);
    // Set seed
    void SetPartitionerSeed(int seed); 
    int GetPartitionerSeed() const; 

  private:
    int num_parts_ = 2;
    std::vector<float> e_wt_factors_; // the cost weight for hyperedge weights
    float path_wt_factor_ = 0.0; // the cost for cut of timing paths
    float snaking_wt_factor_ = 0.0; // the cost for snaking timing paths
    float early_stop_ratio_ = 0.01; // The most improvements are from the earlier moves
    float max_num_moves_ = 25; //50; // The maximum number of moves in each pass
    int max_num_fm_pass_ = 1; //10; // number of fm pass
    int seed_ = 0; // random seed
    utl::Logger* logger_ = nullptr;

    // VertexGain and CompareVertexGain are for priority-queue based max heap
    // define member class for Partitioners
    struct VertexGain {
      int vertex = -1;  // vertex
      float gain = 0.0; // the gain of the vertex
      std::map<int, float> path_cost; // the updated path cost after moving vertex

      // Constructor
      VertexGain(int arg_vertex, float arg_gain) 
        : vertex(arg_vertex), gain(arg_gain) {  }

      VertexGain(int arg_vertex, float arg_gain, std::map<int, float>  arg_path_cost) 
        : vertex(arg_vertex), gain(arg_gain), path_cost(arg_path_cost) {  }
      
      // Comparison (for set-based max heap)
      bool operator < (const VertexGain& vertex_gain) const {
        // first compare gain, the vertex
        if (gain > vertex_gain.gain)
          return true;
        else if (gain == vertex_gain.gain && vertex < vertex_gain.vertex)
          return true;
        return false;
      }

      bool operator == (const VertexGain& vertex_gain) const {
        if (gain == vertex_gain.gain && path_cost == vertex_gain.path_cost)
          return true;
        return false;
      }
    };
    
    
    // utility functions.
    // Get block balance 
    std::vector<std::vector<float> >  GetBlockBalance(const HGraph hgraph,
                                                      std::vector<int> &solution);
    // update the net degree for existing solution
    // for each hyperedge, calculate the number of vertices in each part
    std::vector<std::vector<int> > GetNetDegrees(const HGraph hgraph,
                                                 std::vector<int> &solution);
    
    // Direct Priority-queue FM-based partitioning
    float DirectKWayFMPass(const HGraph hgraph,
                           const std::vector<std::vector<float> > &max_block_balance,
                           std::vector<int> &solution);
    // We enable multithread for initializing gain bucket
    void InitializeGainBucket(std::set<VertexGain> *gain_bucket_ptr, 
                              int to_pid,
                              const std::vector<bool> *visited_ptr,
                              const HGraph hgraph,
                              const std::vector<int> *solution_ptr,
                              const std::vector<float> *cur_path_cost_ptr,
                              const std::vector<std::vector<int> > *net_degs_ptr);

    // We enable multithread for updating gain bucket after moving vertex v
    void UpdateGainBucket(int v, int from_pid, int to_pid,
                          std::set<VertexGain> *gain_bucket_ptr,
                          const HGraph hgraph,
                          const std::vector<bool> *unvisited_ptr,
                          const std::vector<int>  *solution_ptr,
                          const std::vector<float> *cur_path_cost_ptr,
                          const std::vector<std::vector<int> > *net_degs_ptr);
    
    // Calculate the cost for each timing path
    // In the default mode (v = -1, to_pid = -1), 
    // we just calculate the cost for the timing path path_id
    // In the replacement mode, we replace the block id of v to to_pid
    float CalculatePathCost(int path_id, const HGraph hgraph,
                            const std::vector<int> &solution,
                            int v = -1, int to_pid = -1);
    // Calculate the gain for a vertex v
    VertexGain CalculateGain(int v, int from_pid, int to_pid,
                             const HGraph hgraph,
                             const std::vector<int> &solution,
                             const std::vector<float> &cur_path_cost,
                             const std::vector<std::vector<int> > &net_degs);
};

// nickname for shared pointer of Partitioners
using PartitionersPtr = std::shared_ptr<Partitioners>;

}  // namespace par

