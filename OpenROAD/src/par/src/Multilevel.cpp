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

#include "Multilevel.h"

#include <julia.h>
#include <unistd.h>

#include <cstdlib>
#include <string>

#include "Partitioner.h"
#include "TPHypergraph.h"
#include "utl/Logger.h"

JULIA_DEFINE_FAST_TLS

using utl::PAR;

namespace par {

void MultiLevelHierarchy::SpecifySpecParams(std::string hint_hgraph,
                                            int num_parts,
                                            int num_eigen_vectors,
                                            int num_solver_iters,
                                            int expander_cycles,
                                            int embed_placement_dimensions,
                                            int seed,
                                            float ub_factor)
{
  hint_hgraph_ = hint_hgraph;
  hint_solution_file_ = hint_hgraph + ".solution";
  num_parts_ = num_parts;
  num_eigen_vectors_ = num_eigen_vectors;
  num_solver_iters_ = num_solver_iters;
  expander_cycles_ = expander_cycles;
  embed_placement_dimensions_ = embed_placement_dimensions;
  seed_ = seed;
  ub_factor_ = ub_factor;
}

void MultiLevelHierarchy::CallJulia()
{
  jl_init();
  (void) jl_eval_string("import Zhiang");
  //(void)jl_eval_string("Zhiang.TestFunc()");
  jl_module_t* spec_module = (jl_module_t*) jl_eval_string("Zhiang");
  // jl_function_t *func = jl_get_function(spec_module,
  // "FindGlobalClustersConcatenated");
  jl_function_t* func = jl_get_function(spec_module, "FindGlobalClusters");
  jl_value_t* argument1 = jl_cstr_to_string(hint_hgraph_.c_str());
  jl_value_t* argument2 = jl_cstr_to_string(hint_solution_file_.c_str());
  jl_value_t* argument3 = jl_box_int64(num_parts_);
  jl_value_t* argument4 = jl_box_int64(num_eigen_vectors_);
  jl_value_t* argument5 = jl_box_int64(num_solver_iters_);
  jl_value_t* argument6 = jl_box_int64(seed_);
  jl_value_t* argument7 = jl_box_int64(expander_cycles_);
  jl_value_t* argument8 = jl_box_float64(ub_factor_);
  jl_value_t* arguments[8] = {argument1,
                              argument2,
                              argument3,
                              argument4,
                              argument5,
                              argument6,
                              argument7,
                              argument8};

  jl_call(func, arguments, 8);
  jl_atexit_hook(0);
}

void MultiLevelHierarchy::StandardSpec(std::string hypergraph_file)
{
  hint_hgraph_ = hypergraph_file;
  jl_init();
  (void) jl_eval_string("import StandardSpec");
  jl_module_t* spec_module = (jl_module_t*) jl_eval_string("StandardSpec");
  jl_function_t* func = jl_get_function(spec_module, "FindGlobalClusters");
  jl_value_t* argument1 = jl_cstr_to_string(hint_hgraph_.c_str());
  jl_value_t* argument2 = jl_box_int64(num_parts_);
  jl_value_t* argument3 = jl_box_int64(num_eigen_vectors_);
  jl_value_t* argument4 = jl_box_int64(num_solver_iters_);
  jl_value_t* argument5 = jl_box_int64(seed_);
  jl_value_t* argument6 = jl_box_int64(expander_cycles_);
  jl_value_t* arguments[6]
      = {argument1, argument2, argument3, argument4, argument5, argument6};

  jl_call(func, arguments, 6);
  jl_atexit_hook(0);
}

void MultiLevelHierarchy::CallSpecPart(HGraph hgraph,
                                       std::vector<int>& hint_solution)
{
  // write solution file
  WriteSolution(hint_solution_file_.c_str(), hint_solution);
  // int exit_code = CallJulia();
  CallJulia();
  // update community solution
  std::vector<int> community;
  std::ifstream file_input(cluster_file_);
  int part_id = 0;
  for (int i = 0; i < hgraph->num_vertices_; i++) {
    file_input >> part_id;
    community.push_back(part_id);
  }
  file_input.close();

  CutOverlay(hgraph, hint_solution);
  hgraph->community_flag_ = true;
  // hgraph->community_attr_ = hint_solution;
  hgraph->community_attr_ = CutOverlay(hgraph, community);

  /*
  std::cout << "Zhiang : Max community =  ";
  const auto [min, max] = std::minmax_element(begin(hgraph->community_attr_),
                                              end(hgraph->community_attr_));
  std::cout << *max << std::endl;
  std::set<int> num_comm;
  for (auto value : hgraph->community_attr_)
    num_comm.insert(value);
  std::cout << "Error:  Number of comm = " << num_comm.size() << std::endl;
  for (auto value : num_comm)
    std::cout << "value = " << value << std::endl;

  if (hgraph->community_flag_ == false) {
    hgraph->community_flag_ = true;
    hgraph->community_attr_ = community;
  } else {
    hgraph->community_flag_ = true;
    for (auto v = 0;  v < community.size(); v++)
      hgraph->community_attr_[v] =
          hgraph->community_attr_[v] * num_parts_ + community[v];
  }
  // check the number of solution
  for (int i = 0; i < hint_solution.size(); i++) {
    if (hint_solution[i] != community[i]) {
      std::cout << "The solution is different" << std::endl;
      break;
    }
  }
  */

  // hgraph->community_flag_ = true;
  // hgraph->community_attr_ = community;
  // hgraph->community_attr_ = hint_solution;
}

std::vector<int> MultiLevelHierarchy::CallFlow(
    std::string hypergraph_file,
    HGraph hgraph,
    HGraph hgraph_pr,
    std::vector<std::vector<float>> max_block_balance)
{
  std::vector<int> solution;
  edge_mask_.resize(hgraph->num_hyperedges_);
  std::fill(edge_mask_.begin(), edge_mask_.end(), false);
  // setting temp ub factor to 0.5 since it generates better results
  std::vector<std::vector<float>> temp_block_balance
      = hgraph->GetVertexBalance(num_parts_, 0.5);
  // StandardSpec(hypergraph_file);
  solution = RunFlow(hgraph, hgraph_pr, max_block_balance, solution, false);
  return solution;
}

std::vector<int> MultiLevelHierarchy::RunFlow(
    HGraph hgraph,
    HGraph hgraph_pr,
    std::vector<std::vector<float>> max_block_balance,
    std::vector<int> initial_solution,
    bool read_spec_embedding_flag)
{
  if (read_spec_embedding_flag == true) {
    // read embedding information from file
    hgraph->placement_flag_ = true;
    std::ifstream embed_file_input(embedding_file_);
    if (!embed_file_input.is_open()) {
      std::cerr << "Error! Can not open embedding file " << std::endl;
    }
    if (hgraph->placement_dimensions_ == embed_placement_dimensions_) {
      for (int i = 0; i < hgraph->num_vertices_; ++i) {
        std::vector<float> placement_info;
        float value;
        for (int j = 0; j < num_eigen_vectors_; ++j) {
          embed_file_input >> value;
          placement_info.push_back(value);
        }
        // update the placement info
        std::vector<float>::reverse_iterator old_iter
            = hgraph->placement_attr_[i].rbegin();
        std::vector<float>::reverse_iterator new_iter = placement_info.rbegin();
        while (new_iter != placement_info.rend()) {
          *old_iter++ = *new_iter++;
        }
      }
    } else if (hgraph->placement_dimensions_ > 0) {
      hgraph->placement_dimensions_ += num_eigen_vectors_;
      for (int i = 0; i < hgraph->num_vertices_; ++i) {
        std::vector<float> placement_info;
        float value;
        for (int j = 0; j < num_eigen_vectors_ * num_parts_; ++j) {
          embed_file_input >> value;
          placement_info.push_back(value);
        }
        hgraph->placement_attr_[i].insert(hgraph->placement_attr_[i].end(),
                                          placement_info.begin(),
                                          placement_info.end());
      }
    } else {
      hgraph->placement_dimensions_ += num_eigen_vectors_;
      for (int i = 0; i < hgraph->num_vertices_; ++i) {
        std::vector<float> placement_info;
        float value;
        for (int j = 0; j < num_eigen_vectors_; ++j) {
          embed_file_input >> value;
          placement_info.push_back(value);
        }
        hgraph->placement_attr_.push_back(placement_info);
      }
    }
    embed_file_input.close();
  }
  std::vector<HGraph> hgraphs = coarsening_->LazyFirstChoice(hgraph);
  // initial partitioning
  HGraph hgraph_c = hgraphs.back();
  hgraphs.pop_back();
  std::vector<std::vector<float>> temp_block_balance
      = hgraph->GetVertexBalance(num_parts_, 0.5);
  // partitioners_->OptimalInitialPartitionCplex(hgraph_c, temp_block_balance,
  // initial_solution); float cutsize = partitioners_->GoldenEvaluator(hgraph_c,
  // initial_solution, false).first; std::cout << "[INIT PART] ILP " << " Init
  // cutsize " << cutsize << std::endl;
  std::vector<std::vector<int>> solution_set;
  std::map<int, float> cutsize_vec;
  std::mt19937 gen;
  gen.seed(seed_);
  std::uniform_real_distribution<> dist(0.0, 1.0);
  std::vector<int> best_solution;
  float best_cutsize = std::numeric_limits<float>::max();
  for (int i = 0; i < 20; ++i) {
    int seed = 100000 * dist(gen);
    partitioners_->SetPartitionerSeed(seed);
    partitioners_->RandomPartition(
        hgraph_c, max_block_balance, initial_solution);
    partitioners_->DirectKWayFM(hgraph_c, max_block_balance, initial_solution);
    // kpmrefiner_->KPMrefinement(hgraph_c, temp_block_balance,
    // initial_solution);
    float cutsize
        = partitioners_->GoldenEvaluator(hgraph_c, initial_solution, false)
              .first;
    cutsize_vec.insert(std::make_pair(i, cutsize));
    solution_set.push_back(initial_solution);
    std::cout << "[INIT PART] " << i << " Init cutsize " << cutsize
              << std::endl;
    if (cutsize < best_cutsize) {
      best_cutsize = cutsize;
      best_solution = initial_solution;
    }
  }
  std::cout << "[INIT PART] Best init cutsize " << best_cutsize << std::endl;
  initial_solution = best_solution;
  ilprefiner_->SetMaxBalance(max_block_balance);
  /*ilprefiner_->SetCurrBalance(
      partitioners_->GoldenEvaluator(hgraph_c, initial_solution, false).second);*/
  // Running ILP with guidance from best random solution
  std::vector<int> ilp_warm_start = best_solution;
  partitioners_->OptimalInitialPartitionCplexWarmStart(
      hgraph_c, max_block_balance, ilp_warm_start);
  // partitioners_->OptimalInitialPartitionCplex(hgraph_c, max_block_balance,
  // ilp_warm_start);
  solution_set.push_back(ilp_warm_start);
  float ilp_cutsize
      = partitioners_->GoldenEvaluator(hgraph_c, ilp_warm_start, false).first;
  cutsize_vec.insert(std::make_pair(cutsize_vec.size(), ilp_cutsize));
  std::cout << "[INIT PART] 20 "
            << "Init cutsize " << ilp_cutsize << std::endl;
  if (ilp_cutsize < best_cutsize) {
    best_cutsize = ilp_cutsize;
    best_solution = ilp_warm_start;
  }
  std::vector<int> partition_ids(solution_set.size());
  std::iota(partition_ids.begin(), partition_ids.end(), 0);
  std::sort(partition_ids.begin(),
            partition_ids.end(),
            [&](const int x, const int y) {
              return cutsize_vec[x] < cutsize_vec[y];
            });
  std::cout << "[INIT PART] Best init cutsize " << best_cutsize << std::endl;
  initial_solution = best_solution;
  ilprefiner_->SetMaxBalance(max_block_balance);
  // SingleCycle(hgraphs, initial_solution, max_block_balance, false);
  // return initial_solution;

  // partitioners_->GoldenEvaluator(hgraph_c, best_solution, true);
  // running multilevel refinement

  // Refining all 20 solutions across different threads
  /*std::vector<std::thread> threads;
  std::vector<std::vector<HGraph>> hgraphs_vec;
  int num_best_partitions = 3;
  std::cout << "[REFINE] Refining " << num_best_partitions
            << " partitions in parallel" << std::endl;
  for (int i = 0; i < num_best_partitions; ++i) {
    threads.push_back(std::thread(&par::MultiLevelHierarchy::SingleCycle,
                                  this,
                                  std::ref(hgraphs),
                                  std::ref(solution_set[partition_ids[i]]),
                                  max_block_balance,
                                  false));
    // threads.push_back(std::thread(&par::MultiLevelHierarchy::VCycle, this,
    // std::ref(hgraphs), std::ref(solution_set[partition_ids[i]]),
    // max_block_balance));
  }
  for (auto& th : threads) {
    th.join();
  }
  threads.clear();
  best_cutsize = std::numeric_limits<float>::max();
  int best_partition = -1;
  for (int i = 0; i < num_best_partitions; ++i) {
    float refined_cutsize
        = partitioners_
              ->GoldenEvaluator(hgraph, solution_set[partition_ids[i]], false)
              .first;
    std::cout << "[PARTITION] " << i << " Refined cutsize " << refined_cutsize
              << std::endl;
    if (refined_cutsize < best_cutsize) {
      best_cutsize = refined_cutsize;
      best_partition = partition_ids[i];
    }
  }
  initial_solution = solution_set[best_partition];
  std::cout << "[BEST CUT] Best partition has cutsize " << best_cutsize
            << std::endl;
  return initial_solution;*/
  SingleCycle(hgraphs, initial_solution, max_block_balance, false);
  //return initial_solution;
  // allowing marginal imbalance to prevent FM from getting stuck
  std::vector<float> ubfactor_delta;
  int num_ubfactor_delta = 5;
  for (int i = 1; i <= num_ubfactor_delta; i++)
    ubfactor_delta.push_back(-1 * ub_factor_ / 2.0
                             + i * (ub_factor_ / 2.0 / num_ubfactor_delta));
  if (read_spec_embedding_flag == true) {
    ubfactor_delta = std::vector<float>(num_ubfactor_delta, 0.0);
  }
  const int max_num_vcycle = 5;
  if (v_cycle_flag_ == true) {
    int v_cycle_iter = 0;
    float pre_cost
        = partitioners_->GoldenEvaluator(hgraph, initial_solution, false).first;
    float delta_cost = std::numeric_limits<float>::max();
    bool last_run_flag = false;
    while ((v_cycle_iter < max_num_vcycle && delta_cost > 0.0)
           || last_run_flag == true) {
      float ub_factor
          = (last_run_flag == true)
                ? ub_factor_
                : ub_factor_
                      + ubfactor_delta[v_cycle_iter % num_ubfactor_delta];
      std::vector<std::vector<float>> temp_block_balance
          = hgraph->GetVertexBalance(num_parts_, ub_factor);
      logger_->info(PAR, 3939, "V-cycle Iteration {}\n", v_cycle_iter++);
      std::vector<int>& community = hgraph->community_attr_;  // for abbrev.
      // community = CutOverlay(hgraph, initial_solution);
      community = initial_solution;
      hgraph->community_flag_ = true;
      // community = initial_solution;
      for (auto i = 0; i < community.size(); ++i) {
        community[i] = (community[i] * num_parts_) + initial_solution[i];
      }
      // coarse the hgraph with initial_solution as community constraints
      std::vector<HGraph> hgraphs_c = coarsening_->LazyFirstChoice(hgraph);
      // update the initial solution as the coarsest hgraph
      initial_solution.clear();
      for (auto value : hgraphs_c.back()->community_attr_)
        initial_solution.push_back(value % num_parts_);
      /*partitioners_->DirectKWayFM(
          hgraphs_c.back(), max_block_balance, initial_solution);*/

      /*kpmrefiner_->KPMrefinement(
          hgraphs_c.back(), max_block_balance, initial_solution);
      partitioners_->DirectKWayFM(
          hgraphs_c.back(), max_block_balance, initial_solution);*/
      ilprefiner_->SetCurrBalance(
          partitioners_
              ->GoldenEvaluator(hgraphs_c.back(), initial_solution, true)
              .second);
      ilprefiner_->Refine(hgraphs_c.back(), initial_solution);
      partitioners_->DirectKWayFM(
          hgraphs_c.back(), max_block_balance, initial_solution);
      // initial partitioning
      /*partitioners_->DirectKWayFMWithImb(hgraphs_c.back(), max_block_balance,
                                  initial_solution);
      //exit(1);
      float cutsize = partitioners_->GoldenEvaluator(hgraphs_c.back(),
      initial_solution, false).first; */
      hgraphs_c.pop_back();
      // run single-cycle refinement
      SingleCycle(hgraphs_c, initial_solution, temp_block_balance, true);
      const float cur_cost
          = partitioners_->GoldenEvaluator(hgraph, initial_solution, true)
                .first;
      delta_cost = pre_cost - cur_cost;
      pre_cost = cur_cost;
      if (last_run_flag == true)
        break;
      if (delta_cost <= 0.0)
        last_run_flag = false;
    }
  }
  return initial_solution;
}

std::vector<int> MultiLevelHierarchy::SpecRun(
    HGraph hgraph,
    std::vector<std::vector<float>> max_block_balance)
{
  HGraph original_hgraph = hgraph;
  std::vector<int> solution;
  // run with random embedding
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dist(-0.5, 0.5);

  /*
  hgraph->placement_flag_ = true;
  hgraph->placement_dimensions_ = 3;
  hgraph->placement_attr_.clear();
  for (auto i = 0; i < hgraph->num_vertices_; i++) {
    std::vector<float> embed;
    for (auto j = 0; j < hgraph->placement_dimensions_; j++)
      embed.push_back(dist(gen));
    hgraph->placement_attr_.push_back(embed);
  }
  */

  // CallJulia();
  edge_mask_.resize(hgraph->num_hyperedges_);
  std::fill(edge_mask_.begin(), edge_mask_.end(), false);
  // std::vector<std::vector<float> > temp_block_balance =
  // hgraph->GetVertexBalance(num_parts_, 0.25);
  std::vector<std::vector<float>> temp_block_balance
      = hgraph->GetVertexBalance(num_parts_, 0.5);

  solution = RunZhiang(hgraph, temp_block_balance, solution, false);

  // return solution;

  partitioners_->GoldenEvaluator(hgraph, solution, true);
  CallSpecPart(hgraph, solution);
  // update the community
  for (auto v = 0; v < solution.size(); v++)
    hgraph->community_attr_[v]
        = hgraph->community_attr_[v] * num_parts_ + solution[v];

  std::cout << "Zhiang : Max community =  ";
  const auto [min, max] = std::minmax_element(begin(hgraph->community_attr_),
                                              end(hgraph->community_attr_));
  std::cout << *max << std::endl;
  std::set<int> num_comm;
  for (auto value : hgraph->community_attr_)
    num_comm.insert(value);
  std::cout << "Error:  Number of comm = " << num_comm.size() << std::endl;

  solution = RunZhiang(hgraph, max_block_balance, solution, true);
  // std::vector<int> solution = RunZhiang(hgraph, max_block_balance);
  /*
  std::vector<int> solution = InitialRun(hgraph, max_block_balance);
  if (spectral_flag_ == true) {
    CallSpecPart(hgraph, solution);
    coarsening_->max_coarsen_iters_ = coarsening_->max_coarsen_iters_ / 2;
    coarsening_->coarsening_ratio_ = 2.7;
    coarsening_->smallest_v_size_cgraph_ = 2000;
    coarsening_->seed_ = coarsening_->seed_ + 10001;
    solution = RunZhiang(hgraph, max_block_balance, solution,  true);
    //solution = RunZhiang(hgraph, max_block_balance, false);
  }
  */
  return solution;
}

std::vector<int> MultiLevelHierarchy::InitialRun(
    HGraph hgraph_origin,
    std::vector<std::vector<float>> max_block_balance)
{
  HGraph hgraph = std::make_shared<TPHypergraph>(*hgraph_origin);
  std::vector<std::vector<HGraph>> hgraphs_vec(Nruns_);
  std::vector<std::vector<int>> solutions(Nruns_);
  std::vector<int> seeds;
  for (int run_id = 0; run_id < Nruns_; run_id++)
    seeds.push_back(seed_ + run_id);

  std::vector<std::thread> threads;  // for parallel initial partitioning
  for (int run_id = 0; run_id < Nruns_; run_id++)
    threads.push_back(std::thread(&MultiLevelHierarchy::InitialPartitioning,
                                  this,
                                  hgraph,
                                  &max_block_balance,
                                  &hgraphs_vec[run_id],
                                  &solutions[run_id],
                                  seeds[run_id]));
  for (auto& t : threads)
    t.join();  // wait for all threads to finish
  threads.clear();

  std::cout << "Check all the initial solution" << std::endl;
  float best_cutsize = std::numeric_limits<float>::max();
  int best_id = -1;
  for (int i = 0; i < Nruns_; i++) {
    float cutsize
        = (partitioners_->GoldenEvaluator(hgraph, solutions[i], true)).first;
    if (cutsize < best_cutsize) {
      best_cutsize = cutsize;
      best_id = i;
    }
  }

  // SingleCycle(hgraphs_vec[0], solutions[0], max_block_balance);
  return solutions[best_id];
}

// Run Multiple Coarsening with different random seeds
// to improve the stability and reduce the variance
void MultiLevelHierarchy::InitialPartitioning(
    HGraph hgraph,
    const std::vector<std::vector<float>>* max_block_balance,
    std::vector<HGraph>* hgraphs,
    std::vector<int>* initial_solution,
    int seed)
{
  Coarsening coarsening(*coarsening_);
  coarsening.seed_ = seed;
  *hgraphs = coarsening.LazyFirstChoice(hgraph);
  // Initial partitioning
  HGraph hgraph_c = hgraphs->back();
  hgraphs->pop_back();  // get the coarsest hypergraph

  if (hgraph_c->num_hyperedges_ > 4000)
    partitioners_->DirectKWayFM(
        hgraph_c, *max_block_balance, *initial_solution);
  else
    partitioners_->OptimalInitialPartitionCplex(
        hgraph_c, *max_block_balance, *initial_solution);
  SingleCycle(*hgraphs, *initial_solution, *max_block_balance, true);
}

std::vector<int> MultiLevelHierarchy::RunZhiang(
    HGraph hgraph,
    std::vector<std::vector<float>> max_block_balance,
    std::vector<int> initial_solution,
    bool read_spec_embedding_flag)
{
  if (read_spec_embedding_flag == true) {
    // read embedding information
    hgraph->placement_flag_ = true;
    // embedding_file_ = "zhiang.txt";
    std::ifstream embed_file_input(embedding_file_);
    if (!embed_file_input.is_open()) {
      logger_->error(
          PAR, 2916, "Can not open the embedding file : {}", embedding_file_);
    }
    if (hgraph->placement_dimensions_
        == embed_placement_dimensions_) {  // update placement
      for (int i = 0; i < hgraph->num_vertices_; i++) {
        std::vector<float> placement_info;
        float value;
        for (int j = 0; j < num_eigen_vectors_; j++) {
          embed_file_input >> value;
          placement_info.push_back(value);
        }
        // update the placement_info
        std::vector<float>::reverse_iterator old_iter
            = hgraph->placement_attr_[i].rbegin();
        std::vector<float>::reverse_iterator new_iter = placement_info.rbegin();
        while (new_iter != placement_info.rend())
          *old_iter++ = *new_iter++;
      }
    } else if (hgraph->placement_dimensions_ > 0) {
      // add the placement information to default placement information
      hgraph->placement_dimensions_ += num_eigen_vectors_;
      for (int i = 0; i < hgraph->num_vertices_; i++) {
        std::vector<float> placement_info;
        float value;
        for (int j = 0; j < num_eigen_vectors_ * num_parts_; j++) {
          embed_file_input >> value;
          placement_info.push_back(value);
        }
        hgraph->placement_attr_[i].insert(hgraph->placement_attr_[i].end(),
                                          placement_info.begin(),
                                          placement_info.end());
      }
    } else {
      hgraph->placement_dimensions_ += num_eigen_vectors_;
      for (int i = 0; i < hgraph->num_vertices_; i++) {
        std::vector<float> placement_info;
        float value;
        for (int j = 0; j < num_eigen_vectors_ * num_parts_; j++) {
          embed_file_input >> value;
          placement_info.push_back(value);
        }
        hgraph->placement_attr_.push_back(placement_info);
      }
    }
    embed_file_input.close();
  }

  // Do the coarsening first
  std::vector<HGraph> hgraphs = coarsening_->LazyFirstChoice(hgraph);
  // std::vector<int> initial_solution;
  // Initial partitioning
  HGraph hgraph_c = hgraphs.back();
  hgraphs.pop_back();  // get the coarsest hypergraph
  // temp max block balance
  std::vector<std::vector<float>> temp_block_balance
      = hgraph->GetVertexBalance(num_parts_, 0.5);
  if (read_spec_embedding_flag == true) {
    temp_block_balance = hgraph->GetVertexBalance(num_parts_, ub_factor_);
  }

  if (hgraph_c->num_hyperedges_ > 4000 || read_spec_embedding_flag == true) {
    initial_solution.clear();
    for (auto value : hgraph_c->community_attr_)
      initial_solution.push_back(value % num_parts_);
    partitioners_->DirectKWayFM(hgraph_c, temp_block_balance, initial_solution);
  } else {
    partitioners_->OptimalInitialPartitionCplex(
        hgraph_c, temp_block_balance, initial_solution);
  }
  logger_->info(PAR,
                2912,
                "Initial Partitioning\n"
                "Hypergraph information : "
                "\tnum_vertices = {}"
                "\tnum_hyperedges = {}\n",
                hgraph_c->num_vertices_,
                hgraph_c->num_hyperedges_);
  // print basic statistics
  partitioners_->GoldenEvaluator(hgraph_c, initial_solution, true);

  // Multilevel Refinement
  SingleCycle(hgraphs, initial_solution, max_block_balance, true);
  std::cout << "Zhiang :  "
            << "Finish Original SingleCycle" << std::endl;
  // std::vector<float> ubfactor_delta {-1.0, -0.75, -0.5, -0.25, 0.0};
  std::vector<float> ubfactor_delta;
  int num_ubfactor_delta = 5;
  for (int i = 1; i <= num_ubfactor_delta; i++)
    ubfactor_delta.push_back(-1 * ub_factor_ / 2.0
                             + i * (ub_factor_ / 2.0 / num_ubfactor_delta));
  if (read_spec_embedding_flag == true) {
    ubfactor_delta = std::vector<float>(num_ubfactor_delta, 0.0);
  }

  const int max_num_vcycle = 2;

  // std::vector<float> ubfactor_delta {-0.5, -0.25, 0.0};
  // V-cycle improvement
  if (v_cycle_flag_ == true) {
    int v_cycle_iter = 0;
    float pre_cost
        = partitioners_->GoldenEvaluator(hgraph, initial_solution, true).first;
    float delta_cost = std::numeric_limits<float>::max();
    bool last_run_flag = false;

    while ((v_cycle_iter < max_num_vcycle && delta_cost > 0.0)
           || last_run_flag == true) {
      float ub_factor
          = (last_run_flag == true)
                ? ub_factor_
                : ub_factor_
                      + ubfactor_delta[v_cycle_iter % num_ubfactor_delta];
      std::cout << "Zhiang :  ubfactor = " << ub_factor << std::endl;
      max_block_balance = hgraph->GetVertexBalance(num_parts_, ub_factor);
      logger_->info(PAR, 2913, "V-cycle Iteration {}\n", v_cycle_iter++);

      std::vector<int>& community = hgraph->community_attr_;  // for abbrev.
      /*community = CutOverlay(hgraph, initial_solution);
      std::cout << "community.size() = " << community.size() << std::endl;
      hgraph->community_flag_ = true;
      for (auto v = 0;  v < community.size(); v++)
        community[v] = community[v] * num_parts_ + initial_solution[v];*/

      // coarse the hgraph with initial_solution as community constraints
      std::vector<HGraph> hgraphs_c = coarsening_->LazyFirstChoice(hgraph);
      // update the initial solution as the coarsest hgraph
      initial_solution.clear();
      for (auto value : hgraphs_c.back()->community_attr_)
        initial_solution.push_back(value % num_parts_);

      // initial_solution = hgraphs_c.back()->community_attr_ % num_parts_;
      // initial partitioning
      /*partitioners_->DirectKWayFM(hgraphs_c.back(), max_block_balance,
                                  initial_solution);*/

      hgraphs_c.pop_back();
      // run single-cycle refinement
      SingleCycle(hgraphs_c, initial_solution, max_block_balance, true);
      const float cur_cost
          = partitioners_->GoldenEvaluator(hgraph, initial_solution, true)
                .first;
      delta_cost = pre_cost - cur_cost;
      pre_cost = cur_cost;
      if (last_run_flag == true)
        break;
      if (delta_cost <= 0.0)
        last_run_flag = false;
    }
  }
  return initial_solution;
}

/*
// Single-cycle refinement
// project the current solution to previous level of hypergraph's hypergraph
// and call FM engine to refine
void MultiLevelHierarchy::SingleCycle(std::vector<HGraph> &hgraph_vec,
                  std::vector<int>& pre_solution,
                  std::vector<std::vector<float> > max_block_balance)
{
  std::vector<int> solution;
  int refine_iter = 0;
  while (hgraph_vec.empty() == false) {
    HGraph hgraph_c = hgraph_vec.back();
    hgraph_vec.pop_back();
    solution.clear();
    solution.resize(hgraph_c->num_vertices_);
    std::fill(solution.begin(), solution.end(), 0);
    for (int v = 0; v < hgraph_c->num_vertices_; v++) {
      solution[v] = pre_solution[hgraph_c->vertex_c_attr_[v]];
    }
    // FM-based refinement
    partitioners_->DirectKWayFM(hgraph_c, max_block_balance, solution);
    pre_solution = solution; // update the initial solution
    // print basic statistics
    logger_->info(PAR, 2914, "Refinement Iteration {}\n"
                             "Hypergraph information : "
                             "\tnum_vertices = {}"
                             "\tnum_hyperedges = {}\n",
                             refine_iter++,
                             hgraph_c->num_vertices_,
                             hgraph_c->num_hyperedges_);
    partitioners_->GoldenEvaluator(hgraph_c, solution, true);
  }
}
*/

// Single-cycle refinement
// project the current solution to previous level of hypergraph's hypergraph
// and call FM engine to refine
void MultiLevelHierarchy::SingleCycle(
    std::vector<HGraph>& hgraph_vec,
    std::vector<int>& pre_solution,
    std::vector<std::vector<float>> max_block_balance,
    bool v_cycle)
{
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
    ilprefiner_->SetCurrBalance(
        partitioners_->GoldenEvaluator(hgraph_c, solution, false).second);
    ilprefiner_->Refine(hgraph_c, solution);
    partitioners_->DirectKWayFM(hgraph_c, max_block_balance, solution);
    // FM-based refinement
    // partitioners_->DirectKWayFM(hgraph_c, max_block_balance, solution);
    // kpmrefiner_->KPMrefinement(hgraph_c, max_block_balance, solution);
    /*partitioners_->DirectKWayFM(
          hgraph_c, max_block_balance, solution);*/
    pre_solution = solution;  // update the initial solution
    /*float cutsize
        = partitioners_->GoldenEvaluator(hgraph_c, solution, false).first;
    std::cout << "[REFINE] ITER " << refine_iter++
              << " num_vertices = " << hgraph_c->num_vertices_
              << " num_hyperedges = " << hgraph_c->num_hyperedges_
              << " cutsize = " << cutsize << std::endl;*/
    /*if (v_cycle == true) {
      std::cout << "[REFINE] ITER " << refine_iter++
                << " num_vertices = " << hgraph_c->num_vertices_
                << " num_hyperedges = " << hgraph_c->num_hyperedges_
                << " cutsize = " << cutsize << std::endl;
    }*/
  }

  /*
  while (hgraph_vec.empty() == false) {
    HGraph hgraph_c = hgraph_vec.back();
    hgraph_vec.pop_back();
    solution.clear();
    solution.resize(hgraph_c->num_vertices_);
    std::fill(solution.begin(), solution.end(), 0);
    for (int v = 0; v < hgraph_c->num_vertices_; v++) {
      solution[v] = pre_solution[hgraph_c->vertex_c_attr_[v]];
    }
    // FM-based refinement
    partitioners_->DirectKWayFM(hgraph_c, max_block_balance, solution);
    pre_solution = solution; // update the initial solution
    float cutsize = partitioners_->GoldenEvaluator(hgraph_c, solution,
  false).first; std::cout << "[REFINE] ITER " << refine_iter++ << " num_vertices
  = " << hgraph_c->num_vertices_ << " num_hyperedges = " <<
  hgraph_c->num_hyperedges_ << " cutsize = " << cutsize << std::endl;
  }
  */
}

/*
// edge_mask_ : check if current hyperedge has been cut
std::vector<int> MultiLevelHierarchy::CutOverlay(HGraph hgraph,
                                const std::vector<int>& solution)
{
  //return solution;
  // check if the hyperedge is being cut
  for (int e = 0;  e < hgraph->num_hyperedges_; e++) {
    if (edge_mask_[e] == true)
      continue;
    // check the current solutions
    const int start_idx = hgraph->eptr_[e];
    const int end_idx = hgraph->eptr_[e+1];
    if (end_idx - start_idx <= 1) {
      edge_mask_[e] = true;
      continue;
    }
    int init_part = solution[hgraph->eind_[start_idx]];
    int idx = start_idx + 1;
    for (idx; idx < end_idx; idx++) {
      if (init_part != solution[hgraph->eind_[idx]])
        break;
    }
    // if current hyperedge is cut
    if (idx != end_idx)
      edge_mask_[e] = true; // mask the edge
  }
  // find the connected components using BFS
  std::vector<int> community(hgraph->num_vertices_, -1);
  std::queue<int> wave_front;
  int community_id = 0;
  while (true) {
    community_id++;
    int seed = -1;
    // find the first unvisited seed
    for (auto v = 0; v < community.size(); v++) {
      if (community[v] == -1) {
        seed = v;
        break;
      }
    }
    if (seed == -1)
      break;
    wave_front.push(seed);
    community[seed] = community_id;
    // traverse hyperedges
    while (wave_front.empty() == false) {
      const int v = wave_front.front();
      wave_front.pop();
      for (int e_idx = hgraph->vptr_[v]; e_idx < hgraph->vptr_[v+1]; e_idx++) {
        const int e = hgraph->vind_[e_idx];
        if (edge_mask_[e] == true) // remove this hyperedge
          continue;
        for (auto v_idx = hgraph->eptr_[e]; v_idx < hgraph->eptr_[e+1]; v_idx++)
{ if (community[hgraph->eind_[v_idx]] == -1) { community[hgraph->eind_[v_idx]] =
community_id; wave_front.push(hgraph->eind_[v_idx]);
          }
        }
      } // end traverse current hyperedge
    } // finish current hyperedge
  }

  // for debug
  // check the status of community
  std::set<int> comm_set;
  for (auto value : community)
    comm_set.insert(value);
  std::cout << "comm_set.size() = " << comm_set.size() << std::endl;

  return community;
}
*/

// edge_mask_ : check if current hyperedge has been cut
std::vector<int> MultiLevelHierarchy::CutOverlay(
    HGraph hgraph,
    const std::vector<int>& solution)
{
  for (int i = 0; i < hgraph->num_hyperedges_; ++i) {
    if (edge_mask_[i] == true)
      continue;
    const int first_valid_entry = hgraph->eptr_[i];
    const int first_invalid_entry = hgraph->eptr_[i + 1];
    if (first_invalid_entry - first_valid_entry <= 1) {
      edge_mask_[i] = true;
      continue;
    }
    int base_part = solution[hgraph->eind_[first_valid_entry]];
    for (int j = first_valid_entry + 1; j < first_invalid_entry; ++j) {
      const int v = hgraph->eind_[j];
      int ref_part = solution[v];
      if (base_part != ref_part) {
        edge_mask_[i] = true;
        break;
      }
    }
  }
  std::vector<int> community(hgraph->num_vertices_, -1);
  std::queue<int> wave_front;
  int community_id = -1;
  while (true) {
    int seed = -1;
    for (int i = 0; i < hgraph->num_vertices_; ++i) {
      if (community[i] == -1) {
        seed = i;
        break;
      }
    }
    if (seed == -1) {
      break;
    }
    wave_front.push(seed);
    community_id++;
    while (wave_front.empty() == false) {
      int v = wave_front.front();
      community[v] = community_id;
      wave_front.pop();
      const int first_valid_entry = hgraph->vptr_[v];
      const int first_invalid_entry = hgraph->vptr_[v + 1];
      for (int j = first_valid_entry; j < first_invalid_entry; ++j) {
        const int he = hgraph->vind_[j];
        if (edge_mask_[he] == true) {
          continue;
        }
        const int first_valid_entry_ = hgraph->eptr_[he];
        const int first_invalid_entry_ = hgraph->eptr_[he + 1];
        for (int k = first_valid_entry_; k < first_invalid_entry_; ++k) {
          const int v_ = hgraph->eind_[k];
          if (community[v_] == -1) {
            community[v_] = community_id;
            wave_front.push(v_);
          }
        }
      }
    }
  }
  return community;
}

}  // namespace par