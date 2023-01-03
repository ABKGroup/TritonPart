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

#include "Coarsening.h"
#include "KPMRefinement.h"
#include "ILPbasedRefinement.h"
#include "Partitioner.h"
#include "TPHypergraph.h"
#include "Utilities.h"
#include "utl/Logger.h"

namespace par {

enum RefineType
{
  KFMREFINEMENT,  // direct k-way FM refinement
  KPMREFINEMENT,  // pair-wise k-way FM refinement
  KILPREFINEMENT  // pair-wise ILP-based refinement
};

class MultiLevelHierarchy
{
 public:
  explicit MultiLevelHierarchy() {}
  MultiLevelHierarchy(CoarseningPtr coarsening,
                      PartitionersPtr partitioners,
                      KPMrefinementPtr kpmrefiner,
                      IlpRefinerPtr ilprefiner,
                      int num_parts, // number of blocks
                      bool v_cycle_flag, // vcycle flag
                      int num_initial_solutions, // number of initial random solutions
                      int num_best_initial_solutions, // number of best initial solutions
                      int num_ubfactor_delta, // allowing marginal imbalance to improve QoR
                      int max_num_vcycle, // maximum number of vcycles
                      int seed, // random seed
                      float ub_factor, // ubfactor
                      RefineType refine_type, // refinement type
                      utl::Logger* logger)
  {
    coarsening_ = coarsening;
    partitioners_ = partitioners;
    kpmrefiner_ = kpmrefiner;
    ilprefiner_ = ilprefiner;
    num_parts_ = num_parts;
    v_cycle_flag_ = v_cycle_flag;          
    num_initial_solutions_ = num_initial_solutions; 
    num_best_initial_solutions_ = num_best_initial_solutions;
    num_ubfactor_delta_ = num_ubfactor_delta;
    max_num_vcycle_ = max_num_vcycle;
    seed_ = seed;                        
    ub_factor_ = ub_factor;
    refine_type_ = refine_type;
    logger_ = logger;
  }


  // New implementation
  std::vector<int> CallFlow(HGraph hgraph,
                            std::vector<std::vector<float>> max_block_balance);
 private:
  bool v_cycle_flag_ = true;                // enable v-cycle
  CoarseningPtr coarsening_ = nullptr;      // coarsening operator
  PartitionersPtr partitioners_ = nullptr;  // partitioners operator
  KPMrefinementPtr kpmrefiner_ = nullptr;   // kpmrefiner operator
  IlpRefinerPtr ilprefiner_ = nullptr;         // ilprefiner operator
  utl::Logger* logger_ = nullptr;
  int num_parts_ = 2;          // number of blocks
  int num_initial_solutions_ = 20; // number of initial random solutions
  int num_best_initial_solutions_ = 3; // number of best initial solutions
  int num_ubfactor_delta_ = 5; // allowing marginal imbalance to improve QoR
  int max_num_vcycle_ = 5; // maximum number of vcycles
  int seed_ = 0;                        // random seed
  float ub_factor_ = 1.0;               // ubfactor
  RefineType refine_type_ = KPMREFINEMENT; // refinement type
  
  // The main function
  void RunFlow(HGraph hgraph,
               std::vector<std::vector<float>> max_vertex_balance,
               std::vector<int>& solution);
 
  void SingleCycle(std::vector<HGraph> hgraph_vec,
                   std::vector<int>* solution,
                   std::vector<std::vector<float>> max_block_balance,
                   bool v_cycle);
};

using MultiLevelHierarchyPtr = std::shared_ptr<MultiLevelHierarchy>;

}  // namespace par
