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

#include "TPMultilevel.h"

#include <julia.h>
#include <unistd.h>

#include <cstdlib>
#include <string>

#include "TPHypergraph.h"
#include "TPPartitioner.h"
#include "utl/Logger.h"

using utl::PAR;

namespace par {

TP_partition TPmultilevelPartitioner::Partition(
    HGraph hgraph,
    HGraph hgraph_processed,
    matrix<float> max_vertex_balance,
    bool VCycle)
{
  TP_partition solution;
  MultilevelPart(
      hgraph, hgraph_processed, max_vertex_balance, solution, VCycle);
  std::cout << "[debug] cutsize "
            << partitioner_->GoldenEvaluator(hgraph, solution, false).first
            << std::endl;
  return solution;
}

std::pair<matrix<int>, std::vector<int>> TPmultilevelPartitioner::InitialPart(
    HGraph coarsest_hgraph,
    matrix<float>& max_vertex_balance,
    TP_partition& solution)
{
  std::mt19937 gen;
  gen.seed(seed_);
  std::uniform_real_distribution<> dist(0.0, 1.0);
  // set the solution set
  matrix<int> solution_set;
  std::vector<int> cutsize_vec;
  int best_solution_id = -1;
  float best_cutsize = std::numeric_limits<float>::max();
  partitioner_->SetPartitionerChoice(INIT_RANDOM);
  const auto start_timestamp = std::chrono::high_resolution_clock::now();
  std::cout << "[" << std::left << std::setw(5) << "#" << std::left
            << std::setw(10) << "Type" << std::left << std::setw(10)
            << "Cutsize"
            << "]" << std::endl;
  for (int i = 0; i < num_initial_solutions_; ++i) {
    const int seed = std::numeric_limits<int>::max() * dist(gen);
    partitioner_->SetPartitionerSeed(seed);
    partitioner_->Partition(coarsest_hgraph, max_vertex_balance, solution);
    two_way_refiner_->Refine(coarsest_hgraph, max_vertex_balance, solution);
    ilp_refiner_->Refine(coarsest_hgraph, max_vertex_balance, solution);
    greedy_refiner_->Refine(coarsest_hgraph, max_vertex_balance, solution);
    const float cutsize
        = partitioner_->GoldenEvaluator(coarsest_hgraph, solution, false).first;
    /*(std::cout << "[debug] cutsizes " << cutsize_pre_greedy << ", "
              << cutsize_post_ilp << ", " << cutsize_post_greedy << std::endl;*/
    /*std::cout << "[debug] greedy refine " << cutsize_pre_greedy << " and "
              << cutsize_post_greedy << std::endl;*/
    cutsize_vec.push_back(cutsize);
    solution_set.push_back(solution);
    std::cout << "[" << std::left << std::setw(5) << i << std::left
              << std::setw(10) << "Random" << std::left << std::setw(10)
              << cutsize << "]" << std::endl;
    if (cutsize < best_cutsize) {
      best_cutsize = cutsize;
      best_solution_id = i;
    }
  }
  // VILE initial partitioning
  partitioner_->SetPartitionerChoice(INIT_VILE);
  std::vector<int> vile_solution(coarsest_hgraph->num_vertices_);
  partitioner_->Partition(coarsest_hgraph, max_vertex_balance, vile_solution);
  ilp_refiner_->Refine(coarsest_hgraph, max_vertex_balance, vile_solution);
  greedy_refiner_->Refine(coarsest_hgraph, max_vertex_balance, vile_solution);
  float cutsize
      = partitioner_->GoldenEvaluator(coarsest_hgraph, vile_solution, false)
            .first;
  solution_set.push_back(vile_solution);
  cutsize_vec.push_back(cutsize);
  std::cout << "[" << std::left << std::setw(5) << solution_set.size() - 1
            << std::left << std::setw(10) << "vile" << std::left
            << std::setw(10) << cutsize << "]" << std::endl;
  if (cutsize < best_cutsize) {
    best_cutsize = cutsize;
    best_solution_id = solution_set.size() - 1;
  }
  // ILP initial partitioning
  std::vector<int> ilp_part = solution_set[best_solution_id];
  if (coarsest_hgraph->num_hyperedges_ > 1000) {
    partitioner_->SetPartitionerChoice(INIT_DIRECT_WARM_ILP);
    partitioner_->Partition(coarsest_hgraph, max_vertex_balance, ilp_part);

  } else {
    partitioner_->SetPartitionerChoice(INIT_DIRECT_ILP);
    partitioner_->Partition(coarsest_hgraph, max_vertex_balance, ilp_part);
  }
  greedy_refiner_->Refine(coarsest_hgraph, max_vertex_balance, ilp_part);
  cutsize
      = partitioner_->GoldenEvaluator(coarsest_hgraph, ilp_part, false).first;
  solution_set.push_back(ilp_part);
  cutsize_vec.push_back(cutsize);
  std::cout << "[" << std::left << std::setw(5) << solution_set.size() - 1
            << std::left << std::setw(10) << "ilp" << std::left << std::setw(10)
            << cutsize << "]" << std::endl;
  if (cutsize < best_cutsize) {
    best_cutsize = cutsize;
    best_solution_id = solution_set.size() - 1;
  }

  std::cout << "[Best cutsize from cut pool: " << best_cutsize << std::endl;
  auto end_timestamp = std::chrono::high_resolution_clock::now();
  double time_taken = std::chrono::duration_cast<std::chrono::nanoseconds>(
                          end_timestamp - start_timestamp)
                          .count()
                      * 1e-9;
  std::cout << "**Initial partitioning time** " << time_taken << "seconds"
            << std::endl;
  return std::make_pair(solution_set, cutsize_vec);
}

void TPmultilevelPartitioner::Vcycle(std::vector<HGraph> hgraph_vec,
                                     matrix<float>& max_vertex_balance,
                                     TP_partition& solution,
                                     TP_two_way_refining_ptr refiner,
                                     TP_ilp_refiner_ptr i_refiner,
                                     bool print)
{
  std::vector<int> refined_solution;
  if (print == true) {
    std::cout << "[" << std::left << std::setw(6) << "Iter" << std::left
              << std::setw(10) << "Vtxs" << std::left << std::setw(10) << "Hdgs"
              << std::left << std::setw(10) << "Cutsize"
              << "]" << std::endl;
  }
  int num_levels = hgraph_vec.size();
  int limit_levels = static_cast<int>(num_levels / 2);
  float ub_factor_delta = 2.0 / limit_levels;
  std::vector<float> ub_factors_tol(num_levels, 0.0);
  for (int i = 0; i < num_levels; ++i) {
    int tol = (i - limit_levels) * ub_factor_delta;
    if (tol < 0) {
      tol = 0.0;
    }
    ub_factors_tol[i] = ceil(ub_factor_ + tol);
  }
  int iter = 0;
  int hier_ptr = hgraph_vec.size();
  for (auto it = hgraph_vec.crbegin(); it != hgraph_vec.crend(); ++it) {
    --hier_ptr;
    HGraph coarse_hypergraph = *it;
    refined_solution.clear();
    refined_solution.resize(coarse_hypergraph->num_vertices_);
    std::fill(refined_solution.begin(), refined_solution.end(), 0);
    for (int v = 0; v < coarse_hypergraph->num_vertices_; ++v) {
      refined_solution[v] = solution[coarse_hypergraph->vertex_c_attr_[v]];
    }
    /*if (hier_ptr == 0) {
      refiner->SetHeSizeSkip(20);
      refiner->SetMaxMoves(25);
      refiner->SetMaxPasses(1);
    }*/
    /*refiner->BalancePartition(coarse_hypergraph, max_vertex_balance,
    refined_solution); const matrix<float> max_vertex_balance_tol =
    coarse_hypergraph->GetVertexBalance( num_parts_,
    ub_factors_tol[hgraph_vec.size() - 1]);*/
    refiner->Refine(coarse_hypergraph, max_vertex_balance, refined_solution);
    if (coarse_hypergraph->num_vertices_ < 1000) {
      i_refiner->Refine(
          coarse_hypergraph, max_vertex_balance, refined_solution);
    }
    greedy_refiner_->Refine(
        coarse_hypergraph, max_vertex_balance, refined_solution);
    solution = refined_solution;  // update the initial solution
    if (print == true) {
      std::cout << "[" << std::left << std::setw(6) << ++iter << std::left
                << std::setw(10) << coarse_hypergraph->num_vertices_

                << std::left << std::setw(10)
                << coarse_hypergraph->num_hyperedges_ << std::left
                << std::setw(10)
                << partitioner_
                       ->GoldenEvaluator(
                           coarse_hypergraph, refined_solution, false)
                       .first
                << "]" << std::endl;
    }
  }
}

void TPmultilevelPartitioner::GuidedVcycle(TP_partition& solution,
                                           HGraph hgraph,
                                           matrix<float>& max_vertex_balance,
                                           TP_two_way_refining_ptr refiner,
                                           TP_ilp_refiner_ptr i_refiner)
{
  int num_cycles = 0;
  float delta_cost = std::numeric_limits<float>::max();
  float cutsize_before_vcycle
      = partitioner_->GoldenEvaluator(hgraph, solution, false).first;
  int last_run = 0;
  const int stagnation_count = 2;
  int stagnation = 0;
  /*while ((num_cycles < max_num_vcycle_ && delta_cost > 0.0)
         || last_run == false) {*/
  while (true) {
    ++num_cycles;
    hgraph->community_attr_ = solution;
    hgraph->community_flag_ = true;
    // coarse the hgraph with initial_solution as community constraints
    auto hierarchy = coarsener_->LazyFirstChoice(hgraph);
    auto coarsest_hgraph = hierarchy.back();
    coarsest_hgraph->vertex_c_attr_.resize(coarsest_hgraph->num_vertices_);
    std::iota(coarsest_hgraph->vertex_c_attr_.begin(),
              coarsest_hgraph->vertex_c_attr_.end(),
              0);
    // update the initial solution as the coarsest hgraph
    solution.clear();
    for (auto value : hierarchy.back()->community_attr_)
      solution.push_back(value);
    Vcycle(hierarchy, max_vertex_balance, solution, refiner, i_refiner);
    float cutsize_after_vcycle
        = partitioner_->GoldenEvaluator(hgraph, solution, false).first;
    delta_cost = cutsize_before_vcycle - cutsize_after_vcycle;
    cutsize_before_vcycle = cutsize_after_vcycle;
    std::cout << "[v-cycle " << num_cycles << " delta cost " << delta_cost
              << "]" << std::endl;
    if (delta_cost == 0.0) {
      ++stagnation;
      if (stagnation == stagnation_count) {
        break;
      }
    }
    if (num_cycles >= max_num_vcycle_) {
      break;
    }
  }
}

void TPmultilevelPartitioner::MultilevelPart(HGraph hgraph,
                                             HGraph hgraph_processed,
                                             matrix<float>& max_vertex_balance,
                                             TP_partition& solution,
                                             bool VCycle)
{
  coarsener_->SetVertexOrderChoice(RANDOM);
  TP_coarse_graphs hierarchy = coarsener_->LazyFirstChoice(hgraph);
  HGraph coarsest_hgraph = hierarchy.back();
  hierarchy.pop_back();
  two_way_refiner_->SetHeSizeSkip(500);
  auto init_partitions
      = InitialPart(coarsest_hgraph, max_vertex_balance, solution);
  // sort the solution based on cutsize
  matrix<int> partitions_vec = init_partitions.first;
  std::vector<int> cutsize_vec = init_partitions.second;
  std::vector<int> partition_ids(partitions_vec.size());
  std::iota(partition_ids.begin(), partition_ids.end(), 0);
  std::sort(partition_ids.begin(),
            partition_ids.end(),
            [&](const int x, const int y) {
              return cutsize_vec[x] < cutsize_vec[y];
            });
  solution = partitions_vec[partition_ids.front()];
  two_way_refiner_->SetHeSizeSkip(50);
  greedy_refiner_->SetMaxMoves(std::min(coarsest_hgraph->num_hyperedges_, 10));
  ilp_refiner_->SetHeSizeSkip(50);
  std::vector<std::thread> threads;
  matrix<HGraph> hg_threads(GetBestInitSolns());
  std::vector<TP_two_way_refining_ptr> refiner_threads;
  std::vector<TP_ilp_refiner_ptr> i_refiner_threads;
  for (int i = 0; i < GetBestInitSolns(); ++i) {
    for (auto& hg : hierarchy) {
      std::shared_ptr<TPHypergraph> hg_thread(new TPHypergraph(*hg));
      hg_threads[i].push_back(hg_thread);
    }
    TP_two_way_refining_ptr refiner_thread
        = std::make_shared<TPtwoWayFM>(two_way_refiner_->GetNumParts(),
                                       two_way_refiner_->GetIters(),
                                       two_way_refiner_->GetMaxMoves(),
                                       two_way_refiner_->GetRefinerChoice(),
                                       two_way_refiner_->GetSeed(),
                                       two_way_refiner_->GetEdgeWtFactors(),
                                       two_way_refiner_->GetPathWtFactor(),
                                       two_way_refiner_->GetSnakingWtFactor(),
                                       two_way_refiner_->GetLogger());
    TP_ilp_refiner_ptr ilp_thread
        = std::make_shared<TPilpRefine>(ilp_refiner_->GetNumParts(),
                                        ilp_refiner_->GetIters(),
                                        two_way_refiner_->GetMaxMoves(),
                                        ilp_refiner_->GetRefinerChoice(),
                                        ilp_refiner_->GetSeed(),
                                        ilp_refiner_->GetEdgeWtFactors(),
                                        ilp_refiner_->GetPathWtFactor(),
                                        ilp_refiner_->GetSnakingWtFactor(),
                                        ilp_refiner_->GetLogger(),
                                        ilp_refiner_->GetWavefront());
    refiner_thread->SetHeSizeSkip(50);
    ilp_thread->SetHeSizeSkip(50);
    refiner_threads.push_back(refiner_thread);
    i_refiner_threads.push_back(ilp_thread);
  }
  for (int i = 0; i < GetBestInitSolns(); ++i) {
    threads.push_back(std::thread(&par::TPmultilevelPartitioner::Vcycle,
                                  this,
                                  std::ref(hg_threads[i]),
                                  std::ref(max_vertex_balance),
                                  std::ref(partitions_vec[partition_ids[i]]),
                                  refiner_threads[i],
                                  i_refiner_threads[i],
                                  false));
  }
  for (auto& th : threads) {
    th.join();
  }
  threads.clear();
  /*for (int i = 0; i < GetBestInitSolns(); ++i) {
    Vcycle(hierarchy,
           max_vertex_balance,
           partitions_vec[partition_ids[i]],
           two_way_refiner_,
           false);
  }*/
  float best_cut = std::numeric_limits<float>::max();
  int best_partition_id;
  for (int i = 0; i < GetBestInitSolns(); ++i) {
    const float refined_cut
        = two_way_refiner_
              ->CutEvaluator(hgraph, partitions_vec[partition_ids[i]])
              .first;
    std::cout << "[Partition #" << partition_ids[i] << " " << refined_cut << "]"
              << std::endl;
    if (refined_cut < best_cut) {
      best_cut = refined_cut;
      best_partition_id = partition_ids[i];
    }
  }
  solution = partitions_vec[best_partition_id];
  if (VCycle == false) {
    return;
  }
  GuidedVcycle(
      solution, hgraph, max_vertex_balance, two_way_refiner_, ilp_refiner_);
  /*
  std::vector<std::thread> guided_threads;
  std::vector<HGraph> guided_hgraphs;
  std::vector<TP_two_way_refining_ptr> guided_refiners;
  std::vector<TP_ilp_refiner_ptr> guided_irefiners;
  for (int i = 0; i < GetBestInitSolns(); ++i) {
    guided_hgraphs.push_back(hgraph);
    TP_two_way_refining_ptr refiner_thread
        = std::make_shared<TPtwoWayFM>(two_way_refiner_->GetNumParts(),
                                       two_way_refiner_->GetIters(),
                                       two_way_refiner_->GetMaxMoves(),
                                       two_way_refiner_->GetRefinerChoice(),
                                       two_way_refiner_->GetSeed(),
                                       two_way_refiner_->GetEdgeWtFactors(),
                                       two_way_refiner_->GetPathWtFactor(),
                                       two_way_refiner_->GetSnakingWtFactor(),
                                       two_way_refiner_->GetLogger());
    TP_ilp_refiner_ptr ilp_thread
        = std::make_shared<TPilpRefine>(ilp_refiner_->GetNumParts(),
                                        ilp_refiner_->GetIters(),
                                        two_way_refiner_->GetMaxMoves(),
                                        ilp_refiner_->GetRefinerChoice(),
                                        ilp_refiner_->GetSeed(),
                                        ilp_refiner_->GetEdgeWtFactors(),
                                        ilp_refiner_->GetPathWtFactor(),
                                        ilp_refiner_->GetSnakingWtFactor(),
                                        ilp_refiner_->GetLogger(),
                                        ilp_refiner_->GetWavefront());
    refiner_thread->SetHeSizeSkip(50);
    ilp_thread->SetHeSizeSkip(50);
    guided_refiners.push_back(refiner_thread);
    guided_irefiners.push_back(ilp_thread);
  }

  for (int i = 0; i < GetBestInitSolns(); ++i) {
    threads.push_back(std::thread(&par::TPmultilevelPartitioner::GuidedVcycle,
                                  this,
                                  std::ref(partitions_vec[partition_ids[i]]),
                                  std::ref(guided_hgraphs[i]),
                                  std::ref(max_vertex_balance),
                                  guided_refiners[i],
                                  guided_irefiners[i]));
  }
  for (auto& th : threads) {
    th.join();
  }
  threads.clear();

  for (int i = 0; i < GetBestInitSolns(); ++i) {
    std::cout << "[Partition " << i << " cutsize "
              << partitioner_
                     ->GoldenEvaluator(
                         hgraph, partitions_vec[partition_ids[i]], false)
                     .first
              << std::endl;
  }*/
  /*GuidedVcycle(
      solution, hgraph, max_vertex_balance, two_way_refiner_, ilp_refiner_);*/
}
}  // namespace par