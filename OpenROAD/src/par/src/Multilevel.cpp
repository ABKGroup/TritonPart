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

#include "Multilevel.h"

#include <julia.h>
#include <unistd.h>

#include <cstdlib>
#include <string>

#include "Partitioner.h"
#include "TPHypergraph.h"
#include "utl/Logger.h"

using utl::PAR;

namespace par {

std::vector<int> MultiLevelHierarchy::CallFlow(
    HGraph hgraph,  // input hypergraph
    std::vector<std::vector<float>> max_block_balance)
{
  std::vector<int> solution;
  // run the multilevel partitioner
  RunFlow(hgraph, max_block_balance, solution);
  return solution;
}

// The main function for solution
void MultiLevelHierarchy::RunFlow(
    HGraph hgraph,  // hypergraph
    std::vector<std::vector<float>> max_block_balance,
    // solution vector
    std::vector<int>& solution)
{
  // multilevel coarsening
  std::vector<HGraph> hgraphs = coarsening_->LazyFirstChoice(hgraph);
  // initial partitioning
  HGraph hgraph_c = hgraphs.back();  // get the coarsest hypergraph
  hgraphs.pop_back();
  // generate multiple random initial solutions
  // set the random seed
  std::mt19937 gen;
  gen.seed(seed_);
  std::uniform_real_distribution<> dist(0.0, 1.0);
  // set the solution set
  std::vector<std::vector<int>> solution_set;
  std::vector<int> cutsize_vec;
  int best_solution_id = -1;
  float best_cutsize = std::numeric_limits<float>::max();
  // generate num_initial_solutions_ random solutions
  for (int i = 0; i < num_initial_solutions_; ++i) {
    const int seed = std::numeric_limits<int>::max() * dist(gen);
    partitioners_->SetPartitionerSeed(seed);
    // random partitioning
    partitioners_->RandomPartition(hgraph_c, max_block_balance, solution);
    partitioners_->DirectKWayFM(hgraph_c, max_block_balance, solution);
    const float cutsize
        = partitioners_->GoldenEvaluator(hgraph_c, solution, false).first;
    // push solution
    cutsize_vec.push_back(cutsize);
    solution_set.push_back(solution);
    logger_->report(
        "[INIT PART]  solution id :  {}  Init cutsize : {}", i, cutsize);
    if (cutsize < best_cutsize) {
      best_cutsize = cutsize;
      best_solution_id = i;
    }
  }
  logger_->report("[INIT PART] Best cutsize of random solutions : {}",
                  best_cutsize);
  // run ILP-based initial partitioning
  // Running ILP with guidance from best random solution
  ilprefiner_->SetMaxBalance(max_block_balance);
  std::vector<int> ilp_warm_start = solution_set[best_solution_id];
  partitioners_->OptimalInitialPartitionCplexWarmStart(
      hgraph_c, max_block_balance, ilp_warm_start);
  const float ilp_cutsize
      = partitioners_->GoldenEvaluator(hgraph_c, ilp_warm_start, false).first;
  solution_set.push_back(ilp_warm_start);
  cutsize_vec.push_back(ilp_cutsize);
  logger_->report("[INIT PART]  ILP-based initial solution.  Init cutsize : {}",
                  ilp_cutsize);
  if (ilp_cutsize < best_cutsize) {
    best_cutsize = ilp_cutsize;
    best_solution_id = num_initial_solutions_;
  }
  // sort the solution based on cutsize
  std::vector<int> partition_ids(solution_set.size());
  std::iota(partition_ids.begin(), partition_ids.end(), 0);
  std::sort(partition_ids.begin(),
            partition_ids.end(),
            [&](const int x, const int y) {
              return cutsize_vec[x] < cutsize_vec[y];
            });
  logger_->report("[INIT PART] Best init cutsize : {}", best_cutsize);
  // running multilevel refinement
  // Refining top num_best_initial_solutions_ solutions across different threads
  logger_->report("[REFINE] Refining {} initial solutions in parallel",
                  num_best_initial_solutions_);
  std::vector<std::thread> threads;
  std::vector<std::vector<HGraph>> hgraphs_vec(num_best_initial_solutions_);
  for (int i = 0; i < num_best_initial_solutions_; i++) {
    for (auto& hgraph_new : hgraphs) {
      std::shared_ptr<TPHypergraph> hgraph_temp(new TPHypergraph(*hgraph_new));
      hgraphs_vec[i].push_back(hgraph_temp);
    }
  }
  for (int i = 0; i < num_best_initial_solutions_; ++i) {
    threads.push_back(std::thread(&par::MultiLevelHierarchy::SingleCycle,
                                  this,
                                  hgraphs,
                                  &solution_set[partition_ids[i]],
                                  max_block_balance,
                                  true));
  }
  for (auto& th : threads) {
    th.join();
  }
  threads.clear();

  best_cutsize = std::numeric_limits<float>::max();
  best_solution_id = -1;
  for (int i = 0; i < num_best_initial_solutions_; ++i) {
    const float refined_cutsize
        = partitioners_
              ->GoldenEvaluator(hgraph, solution_set[partition_ids[i]], false)
              .first;
    logger_->report(
        "[PARTITION] Refined solution id : {}, refined cutsize : {}",
        i,
        refined_cutsize);
    if (refined_cutsize < best_cutsize) {
      best_cutsize = refined_cutsize;
      best_solution_id = partition_ids[i];
    }
  }
  solution = solution_set[best_solution_id];
  logger_->report("[BEST CUT] Best partition has cutsize : {}", best_cutsize);
  if (v_cycle_flag_ == false)
    return;  // if there is no v_cycle, return

  // *************************************************************
  // Vcycle Refinement starts here
  // *************************************************************
  // allowing marginal imbalance to prevent FM from getting stuck
  std::vector<float> ubfactor_delta;
  for (int i = 1; i <= num_ubfactor_delta_; i++)
    ubfactor_delta.push_back(-1 * ub_factor_ / 2.0
                             + i * (ub_factor_ / 2.0 / num_ubfactor_delta_));
  // vcycle starts
  int v_cycle_iter = 0;
  float pre_cost
      = partitioners_->GoldenEvaluator(hgraph, solution, false).first;
  float delta_cost = std::numeric_limits<float>::max();
  bool last_run_flag = false;  // stops if there is no improvement
  while ((v_cycle_iter < max_num_vcycle_ && delta_cost > 0.0)
         || last_run_flag == true) {
    const float ub_factor
        = (last_run_flag == true)
              ? ub_factor_
              : ub_factor_ + ubfactor_delta[v_cycle_iter % num_ubfactor_delta_];
    const std::vector<std::vector<float>> temp_block_balance
        = hgraph->GetVertexBalance(num_parts_, ub_factor);
    logger_->info(PAR, 3940, "V-cycle Iteration {}\n", v_cycle_iter++);
    hgraph->community_attr_ = solution;
    hgraph->community_flag_ = true;
    // coarse the hgraph with initial_solution as community constraints
    std::vector<HGraph> hgraphs_c = coarsening_->LazyFirstChoice(hgraph);
    // update the initial solution as the coarsest hgraph
    solution.clear();
    for (auto value : hgraphs_c.back()->community_attr_)
      solution.push_back(value % num_parts_);
    // initial partitioning
    if (refine_type_ == KFMREFINEMENT) {
      partitioners_->DirectKWayFM(
          hgraphs_c.back(), max_block_balance, solution);
    } else if (refine_type_ == KPMREFINEMENT) {
      kpmrefiner_->KPMrefinement(hgraphs_c.back(), max_block_balance, solution);
      partitioners_->DirectKWayFM(
          hgraphs_c.back(), max_block_balance, solution);
    } else {
      ilprefiner_->SetCurrBalance(
          partitioners_->GoldenEvaluator(hgraphs_c.back(), solution, true)
              .second);
      ilprefiner_->Refine(hgraphs_c.back(), solution);
      partitioners_->DirectKWayFM(
          hgraphs_c.back(), max_block_balance, solution);
    }
    hgraphs_c.pop_back();
    // run single-cycle refinement
    SingleCycle(hgraphs_c, &solution, temp_block_balance, true);
    const float cur_cost
        = partitioners_->GoldenEvaluator(hgraph, solution, true).first;
    delta_cost = pre_cost - cur_cost;
    pre_cost = cur_cost;
    if (last_run_flag == true)
      break;
    if (delta_cost <= 0.0)
      last_run_flag = true;  // if there i no improvement
  }
}

// Single-cycle refinement
// project the current solution to previous level of hypergraph's hypergraph
// and call FM engine to refine
void MultiLevelHierarchy::SingleCycle(
    std::vector<HGraph> hgraph_vec,
    std::vector<int>* pre_solution_ptr,
    std::vector<std::vector<float>> max_block_balance,
    bool v_cycle)
{
  std::vector<int>& pre_solution = *pre_solution_ptr;
  std::vector<int> solution;
  int refine_iter = 0;
  for (auto it = hgraph_vec.crbegin(); it != hgraph_vec.crend(); ++it) {
    HGraph hgraph_c = *it;
    solution.clear();
    solution.resize(hgraph_c->num_vertices_);
    std::fill(solution.begin(), solution.end(), 0);
    for (int v = 0; v < hgraph_c->num_vertices_; v++) {
      solution[v] = pre_solution[hgraph_c->vertex_c_attr_[v]];
    }
    if (refine_type_ == KFMREFINEMENT) {
      partitioners_->DirectKWayFM(hgraph_c, max_block_balance, solution);
    } else if (refine_type_ == KPMREFINEMENT) {
      kpmrefiner_->KPMrefinement(hgraph_c, max_block_balance, solution);
      partitioners_->DirectKWayFM(hgraph_c, max_block_balance, solution);
    } else {
      std::cout << "ILP Prefiner" << std::endl;
      ilprefiner_->SetCurrBalance(
          partitioners_->GoldenEvaluator(hgraph_c, solution, false).second);
      ilprefiner_->Refine(hgraph_c, solution);
      partitioners_->DirectKWayFM(hgraph_c, max_block_balance, solution);
    }
    refine_iter++;
    pre_solution = solution;  // update the initial solution
  }
}

}  // namespace par
