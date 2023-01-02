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

#include "ILPbasedRefinement.h"

#include "TPHypergraph.h"
#include "Utilities.h"
#include "utl/Logger.h"
// for ILP solver in CPLEX
#include "ilcplex/cplex.h"
#include "ilcplex/ilocplex.h"

namespace par {

IlpGraph::IlpGraph(const int vertex_dimensions,
                   const int hyperedge_dimensions,
                   const bool fixed_flag,
                   std::vector<int> fixed,
                   const matrix<int>& hyperedges,
                   const matrix<float>& vertex_weights,
                   const matrix<float>& hyperedge_weights)
{
  num_vertices_ = static_cast<int>(vertex_weights.size());
  num_hyperedges_ = static_cast<int>(hyperedge_weights.size());
  vertex_dimensions_ = vertex_dimensions;
  hyperedge_dimensions_ = hyperedge_dimensions;
  vertex_weights_ = vertex_weights;
  hyperedge_weights_ = hyperedge_weights;
  fixed_flag_ = fixed_flag;
  fixed_ = fixed;
  // hyperedges
  eind_.clear();
  eptr_.clear();
  eptr_.push_back(static_cast<int>(eind_.size()));
  for (auto hyperedge : hyperedges) {
    eind_.insert(eind_.end(), hyperedge.begin(), hyperedge.end());
    eptr_.push_back(static_cast<int>(eind_.size()));
  }
}

void IlpRefiner::Refine(HGraph hypergraph, std::vector<int>& partition)
{
  auto ilp_refine_time_start = std::chrono::high_resolution_clock::now();
  matrix<int> net_degs(hypergraph->num_hyperedges_,
                       std::vector<int>(num_parts_, 0));
  SetClusterMap(hypergraph->num_vertices_);
  std::vector<float> paths_cost;
  paths_cost.resize(hypergraph->num_timing_paths_);
  for (int path_id = 0; path_id < hypergraph->num_timing_paths_; path_id++) {
    paths_cost[path_id] = CalculatePathCost(path_id, hypergraph, partition);
  }
  partitionqueue pairs = FindMaximalPairs(hypergraph, partition);
  // std::cout << "[DEBUG] Maximal pairs have been calculated " << std::endl;
  std::vector<int> boundary_nodes;
  std::vector<int> ordered_nodes;
  std::vector<int> old_partition = partition;
  matrix<float> vertex_weights_coarsened;
  while (pairs.empty() == false) {
    auto ele = pairs.top();
    const int part_x = ele.GetPartitionX();
    const int part_y = ele.GetPartitionY();
    /*std::cout << "[DEBUG] Pairs calculated " << part_x << " " << part_y
              << std::endl;*/
    FindNetDegs(hypergraph, net_degs, partition);
    boundary_nodes
        = FindBoundaryNodes(hypergraph, part_x, part_y, partition, net_degs);
    /*std::cout << "[DEBUG] Boundary nodes generated " << boundary_nodes.size()
              << std::endl;*/
    OrderVertexSet(net_degs,
                   paths_cost,
                   hypergraph,
                   boundary_nodes,
                   part_x,
                   part_y,
                   partition);
    // std::cout << "[DEBUG] Boundary nodes have been ordered " << std::endl;
    std::vector<int> fixed_vertices;
    GenerateCommunities(hypergraph,
                        part_x,
                        part_y,
                        boundary_nodes,
                        fixed_vertices,
                        vertex_weights_coarsened,
                        partition);
    auto ilpgraph = ContractHypergraph(
        hypergraph, vertex_weights_coarsened, fixed_vertices);
    std::cout << "[DEBUG] ILP hypergraph contains "
              << ilpgraph->GetNumVertices() << " and "
              << ilpgraph->GetNumHyperedges() << std::endl;
    float total_wt = 0.0;
    for (int i = 0; i < ilpgraph->GetNumVertices(); ++i) {
      auto wt = ilpgraph->GetVertexWeight(i);
      for (int j = 0; j < wt.size(); ++j) {
        total_wt += wt[j];
      }
    }
    // std::cout << "[DEBUG] Total weight " << total_wt << std::endl;
    std::vector<int> refined_partition(ilpgraph->GetNumVertices(), -1);
    SolveIlpInstance(ilpgraph, refined_partition, part_x, part_y);
    Remap(partition, refined_partition);
    UpdateNetDegs(
        hypergraph, boundary_nodes, net_degs, old_partition, partition);
    pairs.pop();
    old_partition = partition;
    /*std::cout << "[DEBUG] Partition post ilp-refinement "
              << CutEvaluator(ilpgraph, refined_partition, false).first
              << std::endl;
    std::cout << "[DEBUG] Partition post ilp-refinement "
              << CutEvaluator(hypergraph, partition, false).first <<
    std::endl;*/
    auto ilp_refine_time_end = std::chrono::high_resolution_clock::now();
    double total_ilp_refine_time
        = std::chrono::duration_cast<std::chrono::nanoseconds>(
              ilp_refine_time_end - ilp_refine_time_start)
              .count();
    total_ilp_refine_time *= 1e-9;
    // std::cout << "[REFINE] Total ilp refine runtime " <<
    // total_ilp_refine_time << std::endl;
  }
}

partitiontoken IlpRefiner::CutEvaluator(std::shared_ptr<IlpGraph> hypergraph,
                                        std::vector<int>& solution,
                                        bool print_flag)
{
  matrix<float> block_balance = GetBlockBalance(hypergraph, solution);
  float cost = 0.0;
  // check the cutsize
  for (int e = 0; e < hypergraph->GetNumHyperedges(); ++e) {
    auto indices = hypergraph->GetEdgeIndices(e);
    for (int idx = indices.first + 1; idx < indices.second; ++idx) {
      if (solution[hypergraph->eind_[idx]]
          != solution[hypergraph->eind_[idx - 1]]) {
        auto ewt = hypergraph->GetHyperedgeWeight(e);
        cost += std::inner_product(
            ewt.begin(), ewt.end(), e_wt_factors_.begin(), 0.0);
        break;  // this net has been cut
      }
    }  // finish hyperedge e
  }

  // check if the solution is valid
  for (auto value : solution)
    if (value < 0 || value >= num_parts_)
      std::cout << "[ERROR] The solution is invalid!! " << std::endl;

  if (print_flag == true) {
    // print cost
    std::cout << "[EVAL] Cutsize: " << cost << std::endl;
    // print block balance
    std::vector<float> tot_vertex_weights = hypergraph->GetTotalVertexWeights();
    for (auto block_id = 0; block_id < block_balance.size(); block_id++) {
      std::string line
          = "Vertex balance of block_" + std::to_string(block_id) + " : ";
      for (auto dim = 0; dim < tot_vertex_weights.size(); dim++) {
        std::stringstream ss;  // for converting float to string
        ss << std::fixed << std::setprecision(5)
           << block_balance[block_id][dim] / tot_vertex_weights[dim] << "  ( "
           << block_balance[block_id][dim] << " )  ";
        line += ss.str() + "  ";
      }
      std::cout << "cutsize : " << cost << std::endl;
      std::cout << "balance : " << line << std::endl;
    }  // finish block balance
  }
  return partitiontoken(cost, block_balance);
}

partitiontoken IlpRefiner::CutEvaluator(const HGraph hgraph,
                                        std::vector<int>& solution,
                                        bool print_flag)
{
  matrix<float> block_balance = GetBlockBalance(hgraph, solution);
  float cost = 0.0;
  // check the cutsize
  for (int e = 0; e < hgraph->num_hyperedges_; e++) {
    for (int idx = hgraph->eptr_[e] + 1; idx < hgraph->eptr_[e + 1]; idx++) {
      if (solution[hgraph->eind_[idx]] != solution[hgraph->eind_[idx - 1]]) {
        cost += std::inner_product(hgraph->hyperedge_weights_[e].begin(),
                                   hgraph->hyperedge_weights_[e].end(),
                                   e_wt_factors_.begin(),
                                   0.0);
        break;  // this net has been cut
      }
    }  // finish hyperedge e
  }
  // check timing paths
  for (int path_id = 0; path_id < hgraph->num_timing_paths_; path_id++)
    cost += CalculatePathCost(path_id, hgraph, solution);

  // check if the solution is valid
  for (auto value : solution)
    if (value < 0 || value >= num_parts_)
      std::cout << "[ERROR] The solution is invalid!! " << std::endl;

  if (print_flag == true) {
    // print cost
    std::cout << "[EVAL] Cutsize: " << cost << std::endl;
    // print block balance
    std::vector<float> tot_vertex_weights = hgraph->GetTotalVertexWeights();
    for (auto block_id = 0; block_id < block_balance.size(); block_id++) {
      std::string line
          = "Vertex balance of block_" + std::to_string(block_id) + " : ";
      for (auto dim = 0; dim < tot_vertex_weights.size(); dim++) {
        std::stringstream ss;  // for converting float to string
        ss << std::fixed << std::setprecision(5)
           << block_balance[block_id][dim] / tot_vertex_weights[dim] << "  ( "
           << block_balance[block_id][dim] << " )  ";
        line += ss.str() + "  ";
      }
      std::cout << "cutsize : " << cost << std::endl;
      std::cout << "balance : " << line << std::endl;
    }  // finish block balance
  }
  return partitiontoken(cost, block_balance);
}

inline matrix<float> IlpRefiner::GetBlockBalance(
    std::shared_ptr<IlpGraph> hypergraph,
    std::vector<int>& partition)
{
  matrix<float> block_balance(
      num_parts_, std::vector<float>(hypergraph->GetVertexDimensions(), 0.0));
  if (hypergraph->GetNumVertices() != static_cast<int>(partition.size()))
    return block_balance;
  // check if the solution is valid
  for (auto block_id : partition)
    if (block_id < 0 || block_id >= num_parts_)
      return block_balance;
  // update the block_balance
  for (int v = 0; v < hypergraph->GetNumVertices(); v++)
    block_balance[partition[v]]
        = block_balance[partition[v]] + hypergraph->GetVertexWeight(v);
  return block_balance;
}

inline matrix<float> IlpRefiner::GetBlockBalance(const HGraph hgraph,
                                                 std::vector<int>& partition)
{
  matrix<float> block_balance(
      num_parts_, std::vector<float>(hgraph->vertex_dimensions_, 0.0));
  if (hgraph->num_vertices_ != static_cast<int>(partition.size()))
    return block_balance;
  // check if the solution is valid
  for (auto block_id : partition)
    if (block_id < 0 || block_id >= num_parts_)
      return block_balance;
  // update the block_balance
  for (int v = 0; v < hgraph->num_vertices_; v++)
    block_balance[partition[v]]
        = block_balance[partition[v]] + hgraph->vertex_weights_[v];
  return block_balance;
}

inline void IlpRefiner::Remap(std::vector<int>& partition,
                              std::vector<int>& refined_partition)
{
  for (int i = 0; i < cluster_map_.size(); ++i) {
    partition[i] = refined_partition[cluster_map_[i]];
  }
}

void IlpRefiner::UpdateNetDegs(HGraph hypergraph,
                               std::vector<int>& vertices,
                               matrix<int>& net_degs,
                               std::vector<int>& prev_partition,
                               std::vector<int>& new_partition)
{
  int vertices_span
      = wavefront_ < vertices.size() ? wavefront_ : vertices.size();
  for (int i = 0; i < vertices_span; ++i) {
    const int v = vertices[i];
    const int first_valid_entry = hypergraph->vptr_[v];
    const int first_invalid_entry = hypergraph->vptr_[v + 1];
    const int old_part = prev_partition[v];
    const int new_part = new_partition[v];
    if (old_part == new_part) {
      continue;
    }
    for (int j = first_valid_entry; j < first_invalid_entry; ++j) {
      const int he = hypergraph->vind_[j];
      --net_degs[he][old_part];
      ++net_degs[he][new_part];
    }
  }
}

std::vector<int> IlpRefiner::FindBoundaryNodes(HGraph hypergraph,
                                               const int& part_x,
                                               const int& part_y,
                                               std::vector<int>& partition,
                                               matrix<int>& net_degs)
{
  std::vector<int> boundary_nodes;
  std::vector<bool> vertex_watermark(hypergraph->num_vertices_, false);
  for (int i = 0; i < hypergraph->num_vertices_; ++i) {
    if (partition[i] != part_x && partition[i] != part_y) {
      //|| hypergraph->fixed_attr_[i] > -1) {
      continue;
    }
    const int first_valid_entry = hypergraph->vptr_[i];
    const int first_invalid_entry = hypergraph->vptr_[i + 1];
    for (int j = first_valid_entry; j < first_invalid_entry; ++j) {
      const int he = hypergraph->vind_[j];
      if (net_degs[he][part_x] > 0 && net_degs[he][part_y] > 0) {
        boundary_nodes.push_back(i);
        vertex_watermark[i] = true;
        break;
      }
    }
  }
  std::set<int> one_hops;
  // Also add immediate neighbors of the boundary nodes
  for (auto v : boundary_nodes) {
    const int first_valid_entry = hypergraph->vptr_[v];
    const int first_invalid_entry = hypergraph->vptr_[v + 1];
    for (int i = first_valid_entry; i < first_invalid_entry; ++i) {
      const int he = hypergraph->vind_[i];
      const int first_valid_entry_e = hypergraph->eptr_[he];
      const int first_invalid_entry_e = hypergraph->eptr_[he + 1];
      for (int j = first_valid_entry_e; j < first_invalid_entry_e; ++j) {
        const int nbr = hypergraph->eind_[j];
        if (vertex_watermark[nbr] == true || nbr == v
            || (partition[nbr] != part_x && partition[nbr] != part_y)) {
          continue;
        }
        vertex_watermark[nbr] = true;
        one_hops.insert(nbr);
      }
    }
  }
  std::vector<int> one_hopes_vec(one_hops.begin(), one_hops.end());
  boundary_nodes.insert(
      boundary_nodes.end(), one_hopes_vec.begin(), one_hopes_vec.end());
  return boundary_nodes;
}

void IlpRefiner::SolveIlpInstance(std::shared_ptr<IlpGraph> hypergraph,
                                  std::vector<int>& refined_partition,
                                  const int& part_x,
                                  const int& part_y)
{
  // CPLEX-Based Implementation
  // define environment
  IloEnv myenv;
  IloModel mymodel(myenv);
  // Define constraints
  // For each vertex, define a variable x
  // For each hyperedge, define a variable y
  IloArray<IloNumVarArray> x(myenv, num_parts_);
  IloArray<IloNumVarArray> y(myenv, num_parts_);
  // only consider the top 100 hyperedges based on weight
  // order hyperedges based on decreasing order
  struct comp
  {
    // comparator function
    bool operator()(const std::pair<int, float>& l,
                    const std::pair<int, float>& r) const
    {
      if (l.second != r.second)
        return l.second > r.second;
      return l.first < r.first;
    }
  };

  // use set data structure to sort unvisited vertices
  std::set<std::pair<int, float>, comp> unvisited_hyperedges;
  for (auto e = 0; e < hypergraph->GetNumHyperedges(); ++e) {
    auto he_wt = hypergraph->GetHyperedgeWeight(e);
    const float score = std::inner_product(
        he_wt.begin(), he_wt.end(), e_wt_factors_.begin(), 0.0);
    unvisited_hyperedges.insert(std::pair<int, float>(e, score));
  }
  int max_num_hyperedges = 0;
  std::vector<int> edge_mask;
  float base_score = (*unvisited_hyperedges.begin()).second / 10.0;
  for (auto& value : unvisited_hyperedges) {
    /*if (value.second <= base_score || max_num_hyperedges >= 50)
      break;*/
    edge_mask.push_back(value.first);
    max_num_hyperedges++;
  }
  for (int i = 0; i < num_parts_; ++i) {
    x[i] = IloNumVarArray(myenv, hypergraph->GetNumVertices(), 0, 1, ILOINT);
    y[i] = IloNumVarArray(myenv, edge_mask.size(), 0, 1, ILOINT);
  }
  if (hypergraph->CheckFixedFlag() == true) {
    for (int i = 0; i < hypergraph->GetNumVertices(); ++i) {
      if (hypergraph->CheckFixedStatus(i) == true) {
        mymodel.add(x[hypergraph->GetFixedPart(i)][i] == 1);
      } else {
        for (int j = 0; j < num_parts_; ++j) {
          if (j != part_x && j != part_y) {
            mymodel.add(x[j][i] == 0);
          }
        }
      }
    }
  }
  // define constraints
  // balance constraint
  // check each dimension
  for (int i = 0; i < hypergraph->GetVertexDimensions(); ++i) {
    // allowed balance for each dimension
    for (int j = 0; j < num_parts_; ++j) {
      IloExpr balance_expr(myenv);
      for (int k = 0; k < hypergraph->GetNumVertices(); ++k) {
        auto vwt = hypergraph->GetVertexWeight(k);
        balance_expr += vwt[i] * x[j][k];
      }  // finish traversing vertices
      mymodel.add(IloRange(myenv, 0.0, balance_expr, max_block_balance_[j][i]));
      balance_expr.end();
    }  // finish traversing blocks
  }    // finish dimension check
  // each vertex can only belong to one part between part_x and part_y
  for (int i = 0; i < hypergraph->GetNumVertices(); ++i) {
    if (hypergraph->CheckFixedStatus(i) == true) {
      continue;
    } else {
      IloExpr vertex_expr(myenv);
      vertex_expr += x[part_x][i];
      vertex_expr += x[part_y][i];
      mymodel.add(vertex_expr == 1);
      vertex_expr.end();
    }
  }
  // Hyperedge constraint
  for (int i = 0; i < edge_mask.size(); ++i) {
    const int e = edge_mask[i];
    std::pair<int, int> indices = hypergraph->GetEdgeIndices(e);
    const int first_valid_entry = indices.first;
    const int first_invalid_entry = indices.second;
    for (int j = first_valid_entry; j < first_invalid_entry; ++j) {
      const int vertex_id = hypergraph->eind_[j];
      for (int k = 0; k < num_parts_; ++k) {
        mymodel.add(y[k][e] <= x[k][vertex_id]);
      }
    }
  }
  // Maximize cutsize objective
  IloExpr obj_expr(myenv);  // empty expression
  for (int i = 0; i < edge_mask.size(); ++i) {
    auto ewt = hypergraph->GetHyperedgeWeight(edge_mask[i]);
    const float cost_value = std::inner_product(
        ewt.begin(), ewt.end(), e_wt_factors_.begin(), 0.0);
    for (int j = 0; j < num_parts_; ++j) {
      obj_expr += cost_value * y[j][i];
    }
  }
  mymodel.add(IloMaximize(myenv, obj_expr));  // adding minimization objective
  obj_expr.end();                             // clear memory
  // Model Solution
  IloCplex mycplex(myenv);
  mycplex.extract(mymodel);
  mycplex.setOut(myenv.getNullStream());
  // mycplex.setParam(IloCplex::Param::MIP::Limits::Solutions, 5);
  // mycplex.setParam(IloCplex::Param::TimeLimit, 5);
  mycplex.solve();
  IloBool feasible = mycplex.solve();
  if (feasible == IloTrue) {
    for (int i = 0; i < hypergraph->GetNumVertices(); ++i) {
      for (int j = 0; j < num_parts_; ++j) {
        if (mycplex.getValue(x[j][i]) == 1.00) {
          refined_partition[i] = j;
        }
      }
    }
    // some solution may invalid due to the limitation of ILP solver
    for (auto& value : refined_partition)
      value = (value == -1) ? 0 : value;
  } else {
    /* for Ilp infeasibility debug
    mycplex.exportModel("model.mps");
    DebugIlpInstance("model.mps"); */
    std::cout << "[ILP-Refine] ILP could not be run into completion. Running "
                 "K-way-FM instead!"
              << std::endl;
  }
  // closing the model
  mycplex.clear();
  myenv.end();
}

void IlpRefiner::DebugIlpInstance(const char* file_name)
{
  IloEnv env;
  IloModel model(env);
  IloCplex cplex(env);
  IloObjective obj;
  IloNumVarArray var(env);
  IloRangeArray rng(env);
  IloSOS1Array sos1(env);
  IloSOS2Array sos2(env);
  IloRangeArray lazy(env);
  IloRangeArray cuts(env);
  cplex.importModel(model, file_name, obj, var, rng, sos1, sos2, lazy, cuts);
  cplex.extract(model);
  if (lazy.getSize() > 0)
    cplex.addLazyConstraints(lazy);
  if (cuts.getSize() > 0)
    cplex.addUserCuts(cuts);
  cplex.solve();
  if ((cplex.getStatus() == IloAlgorithm::Infeasible)
      || (cplex.getStatus() == IloAlgorithm::InfeasibleOrUnbounded)) {
    std::cout << std::endl
              << "No solution - starting Conflict refinement" << std::endl;

    IloConstraintArray infeas(env);
    IloNumArray preferences(env);

    infeas.add(rng);
    infeas.add(sos1);
    infeas.add(sos2);
    if (lazy.getSize() || cuts.getSize()) {
      std::cout << "Lazy Constraints and User Cuts ignored" << std::endl;
    }
    for (IloInt i = 0; i < var.getSize(); i++) {
      if (var[i].getType() != IloNumVar::Bool) {
        infeas.add(IloBound(var[i], IloBound::Lower));
        infeas.add(IloBound(var[i], IloBound::Upper));
      }
    }
    for (IloInt i = 0; i < infeas.getSize(); i++) {
      preferences.add(1.0);  // user may wish to assign unique preferences
    }
    if (cplex.refineConflict(infeas, preferences)) {
      IloCplex::ConflictStatusArray conflict = cplex.getConflict(infeas);
      env.getImpl()->useDetailedDisplay(IloTrue);
      std::cout << "Conflict :" << std::endl;
      for (IloInt i = 0; i < infeas.getSize(); i++) {
        if (conflict[i] == IloCplex::ConflictMember)
          std::cout << "Proved  : " << infeas[i] << std::endl;
        if (conflict[i] == IloCplex::ConflictPossibleMember)
          std::cout << "Possible: " << infeas[i] << std::endl;
      }
    } else
      std::cout << "Conflict could not be refined" << std::endl;
    std::cout << std::endl;
  }
  env.out() << "Solution status = " << cplex.getStatus() << std::endl;
  env.out() << "Solution value  = " << cplex.getObjValue() << std::endl;
  IloNumArray vals(env);
  cplex.getValues(vals, var);
  env.out() << "Values        = " << vals << std::endl;
}

partitionqueue IlpRefiner::FindMaximalPairs(HGraph hypergraph,
                                            std::vector<int>& partition)
{
  partitionqueue ppairs(comp);
  for (int i = 0; i < num_parts_; ++i) {
    for (int j = i + 1; j < num_parts_; ++j) {
      float score = FindPairScore(hypergraph, partition, i, j);
      ppairs.push(PartitionPair(i, j, score));
    }
  }
  return ppairs;
}

float IlpRefiner::FindPairScore(HGraph hypergraph,
                                std::vector<int>& partition,
                                int& pair_x,
                                int& pair_y)
{
  // Score between pairs are calculated based on the connectivity of vertices
  // belonging to the partitions
  float score = 0.0;
  for (int i = 0; i < hypergraph->num_hyperedges_; ++i) {
    const int first_valid_entry = hypergraph->eptr_[i];
    const int first_invalid_entry = hypergraph->eptr_[i + 1];
    bool flag_partition_x = false;
    bool flag_partition_y = false;
    for (int j = first_valid_entry; j < first_invalid_entry; ++j) {
      int v = hypergraph->eind_[j];
      if (partition[v] == pair_x) {
        flag_partition_x = true;
      }
      if (partition[v] == pair_y) {
        flag_partition_y = true;
      }
      if (flag_partition_x == true && flag_partition_y == true) {
        score += std::inner_product(hypergraph->hyperedge_weights_[i].begin(),
                                    hypergraph->hyperedge_weights_[i].end(),
                                    e_wt_factors_.begin(),
                                    0.0);
        break;
      }
    }
  }
  return score;
}

bool IlpRefiner::CheckBalance(const int& vertex_id,
                              const int& from_part,
                              const int& to_part,
                              HGraph hypergraph)
{
  if (curr_block_balance_[to_part] + hypergraph->vertex_weights_[vertex_id]
      > max_block_balance_[to_part]) {
    return false;
  } else {
    return true;
  }
}

void IlpRefiner::OrderVertexSet(matrix<int> net_degs,
                                const std::vector<float>& cur_path_cost,
                                HGraph hypergraph,
                                std::vector<int>& boundary_vertices,
                                int from_pid,
                                int to_pid,
                                std::vector<int>& partition)
{
  std::vector<float> gains_boundary_vertices(hypergraph->num_vertices_, 0.0);
  for (int i = 0; i < boundary_vertices.size(); ++i) {
    const int v = boundary_vertices[i];
    gains_boundary_vertices[v] = CalculateGain(
        v, from_pid, to_pid, hypergraph, partition, cur_path_cost, net_degs);
  }
  auto gain_comparison = [&](int x, int y) {
    if (gains_boundary_vertices[x] > gains_boundary_vertices[y]) {
      const int from = partition[x];
      const int to = partition[x] == from_pid ? to_pid : from_pid;
      if (CheckBalance(x, from, to, hypergraph) == true) {
        return true;
      }
      return false;
    }
    return false;
    // return gains_boundary_vertices[x] > gains_boundary_vertices[y]
    //       && CheckBalance(x, from_pid, to_pid, hypergraph) == true;
  };
  std::sort(
      boundary_vertices.begin(), boundary_vertices.end(), gain_comparison);
}

void IlpRefiner::FindNetDegs(HGraph hypergraph,
                             matrix<int>& net_degs,
                             std::vector<int>& partition)
{
  for (int e = 0; e < hypergraph->num_hyperedges_; ++e) {
    for (int idx = hypergraph->eptr_[e]; idx < hypergraph->eptr_[e + 1];
         ++idx) {
      net_degs[e][partition[hypergraph->eind_[idx]]]++;
    }
  }
}

float IlpRefiner::CalculatePathCost(int path_id,
                                    const HGraph hypergraph,
                                    const std::vector<int>& partition,
                                    int v,
                                    int to_pid)
{
  float cost = 0.0;  // cost for current path
  if (hypergraph->num_timing_paths_ == 0)
    return cost;                     // no timing paths
  std::vector<int> path;             // represent the path in terms of block_id
  std::map<int, int> block_counter;  // block_id counter
  for (auto idx = hypergraph->vptr_p_[path_id];
       idx < hypergraph->vptr_p_[path_id];
       idx++) {
    const int u = hypergraph->vind_p_[idx];  // current vertex
    const int block_id = (u == v) ? to_pid : partition[u];
    if (path.size() == 0 || path.back() != block_id) {
      path.push_back(block_id);
      if (block_counter.find(block_id) != block_counter.end())
        block_counter[block_id] += 1;
      else
        block_counter[block_id] = 1;
    }
  }
  if (path.size() <= 1)
    return cost;
  // num_cut = path.size() - 1
  cost = path_wt_factor_ * static_cast<float>(path.size() - 1);
  // get the snaking factor of the path (maximum repetition of block_id - 1)
  int snaking_factor = 0;
  for (auto [block_id, count] : block_counter)
    if (count > snaking_factor)
      snaking_factor = count;
  cost += snaking_wt_factor_ * static_cast<float>(snaking_factor - 1);
  return cost;
}

float IlpRefiner::CalculateGain(int v,
                                int from_pid,
                                int to_pid,
                                const HGraph hypergraph,
                                const std::vector<int>& partition,
                                const std::vector<float>& cur_path_cost,
                                const std::vector<std::vector<int>>& net_degs)
{
  float score = 0.0;
  std::map<int, float> path_cost;  // map path_id to latest score
  if (from_pid == to_pid)
    return score;
  // define two lamda function
  // 1) for checking connectivity
  // 2) for calculating the score of a hyperedge
  // function : check the connectivity for the hyperedge
  auto GetConnectivity = [&](int e) {
    int connectivity = 0;
    for (auto& num_v : net_degs[e])
      if (num_v > 0)
        connectivity++;
    return connectivity;
  };
  // function : calculate the score for the hyperedge
  auto GetHyperedgeScore = [&](int e) {
    return std::inner_product(hypergraph->hyperedge_weights_[e].begin(),
                              hypergraph->hyperedge_weights_[e].end(),
                              e_wt_factors_.begin(),
                              0.0);
  };
  // traverse all the hyperedges connected to v
  const int first_valid_entry = hypergraph->vptr_[v];
  const int first_invalid_entry = hypergraph->vptr_[v + 1];
  // std::cout << std::endl << "[DEBUG] Vertex " << v << std::endl;
  for (auto e_idx = first_valid_entry; e_idx < first_invalid_entry; e_idx++) {
    const int e = hypergraph->vind_[e_idx];  // hyperedge id
    const int connectivity = GetConnectivity(e);
    const float e_score = GetHyperedgeScore(e);
    if (connectivity == 1
        && net_degs[e][from_pid]
               > 1) {  // move from_pid to to_pid will have negative socre
      score -= e_score;
    } else if (connectivity == 2 && net_degs[e][from_pid] == 1
               && net_degs[e][to_pid] > 0) {
      score += e_score;  // after move, all the vertices in to_pid
    }
  }
  // check the timing path
  if (hypergraph->num_timing_paths_ > 0) {
    for (auto p_idx = hypergraph->pptr_v_[v];
         p_idx < hypergraph->pptr_v_[v + 1];
         p_idx++) {
      const int path_id = hypergraph->pind_v_[p_idx];
      const float cost
          = CalculatePathCost(path_id, hypergraph, partition, v, to_pid);
      path_cost[path_id] = cost;
      score += cur_path_cost[path_id] - cost;
    }
  }
  return score;
}

void IlpRefiner::GenerateCommunities(HGraph hypergraph,
                                     const int& part_x,
                                     const int& part_y,
                                     std::vector<int>& ordered_set,
                                     std::vector<int>& fixed_vertices,
                                     matrix<float>& vertex_weights_coarsened,
                                     std::vector<int>& partition)
{
  int wavefront = wavefront_;
  int cluster_id = 0;
  std::vector<int> partition_mask(num_parts_, -1);
  ResetClusterMap();
  std::vector<int> border_population(num_parts_, 0);
  int min_size = ordered_set.size() < wavefront_ ? ordered_set.size() : wavefront_;
  for (int i = 0; i < min_size; ++i) {
    const int p = partition[ordered_set[i]];
    border_population[p]++;
  }
  std::vector<int> global_population(num_parts_, 0);
  for (int i = 0; i < hypergraph->num_vertices_; ++i) {
    const int p = partition[i];
    global_population[p]++;
  }
  int partition_mask_idx = -1;
  for (int i = 0; i < num_parts_; ++i) {
    if (global_population[i] != border_population[i]) {
      partition_mask[i] = ++partition_mask_idx;
    }
  }
  std::cout << "[DEBUG] Wave front size " << wavefront_ << " and set "
            << ordered_set.size() << std::endl;
  if (wavefront_ > ordered_set.size()) {
    wavefront = ordered_set.size();
    //const int ilp_vertices = ordered_set.size() + num_parts_;
    const int ilp_vertices = ordered_set.size() + partition_mask_idx + 1;
    vertex_weights_coarsened.resize(ilp_vertices);
    fixed_vertices.resize(ilp_vertices);
    cluster_id = 0;
  } else {
    const int ilp_vertices = wavefront_ + num_parts_;
    vertex_weights_coarsened.resize(ilp_vertices);
    fixed_vertices.resize(ilp_vertices);
  }
  std::fill(vertex_weights_coarsened.begin(),
            vertex_weights_coarsened.end(),
            std::vector<float>(hypergraph->vertex_dimensions_, 0.0));
  for (int i = 0; i < wavefront; ++i) {
    const int vertex = ordered_set[i];
    cluster_map_[vertex] = cluster_id;
    vertex_weights_coarsened[cluster_id] = hypergraph->vertex_weights_[vertex];
    ++cluster_id;
  }
  std::fill(fixed_vertices.begin(), fixed_vertices.end(), -1);
  for (int i = wavefront; i < ordered_set.size(); ++i) {
    const int vertex = ordered_set[i];
    const int pid = partition[vertex];
    //const int cid = pid + cluster_id;
    const int cid = partition_mask[pid] + cluster_id;
    cluster_map_[vertex] = cid;
    vertex_weights_coarsened[cid]
        = vertex_weights_coarsened[cid] + hypergraph->vertex_weights_[vertex];
    if (fixed_vertices[cid] > -1) {
      continue;
    }
    fixed_vertices[cid] = pid;
  }
  for (int i = 0; i < hypergraph->num_vertices_; ++i) {
    const int pid = partition[i];
    //const int cid = pid + cluster_id;
    const int cid = partition_mask[pid] + cluster_id;
    if (cluster_map_[i] == -1) {
      cluster_map_[i] = cid;
      vertex_weights_coarsened[cid]
          = vertex_weights_coarsened[cid] + hypergraph->vertex_weights_[i];
    }
    if (fixed_vertices[cid] > -1) {
      continue;
    }
    fixed_vertices[cid] = pid;
  }
  std::cout << "[DEBUG] Max clusterid "
            << *std::max_element(cluster_map_.begin(), cluster_map_.end());
  // std::cout << "[DEBUG] Cluster map has been generated " << std::endl;
}

std::shared_ptr<IlpGraph> IlpRefiner::ContractHypergraph(
    HGraph hypergraph,
    matrix<float>& vertex_weights_c,
    std::vector<int>& fixed_vertices)
{
  matrix<int> hyperedges_c;
  matrix<float> hyperedge_weights_c;
  std::map<long long int, int> hash_map;
  for (int i = 0; i < hypergraph->num_hyperedges_; ++i) {
    const int first_valid_entry = hypergraph->eptr_[i];
    const int first_invalid_entry = hypergraph->eptr_[i + 1];
    const int he_size = first_invalid_entry - first_valid_entry;
    if (he_size <= 1 || he_size > GetHyperedgeThreshold()) {
      continue;
    }
    std::set<int> hyperedge_c;
    for (int j = first_valid_entry; j < first_invalid_entry; ++j) {
      const int vtx = hypergraph->eind_[j];
      hyperedge_c.insert(cluster_map_[vtx]);
    }
    if (hyperedge_c.size() <= 1) {
      continue;
    }
    const long long int hash_value
        = std::inner_product(hyperedge_c.begin(),
                             hyperedge_c.end(),
                             hyperedge_c.begin(),
                             static_cast<long long int>(0));
    if (hash_map.find(hash_value) == hash_map.end()) {
      hash_map[hash_value] = static_cast<int>(hyperedges_c.size());
      hyperedge_weights_c.push_back(hypergraph->hyperedge_weights_[i]);
      hyperedges_c.push_back(
          std::vector<int>(hyperedge_c.begin(), hyperedge_c.end()));
    } else {
      const int hash_id = hash_map[hash_value];
      std::vector<int> hyperedge_vec(hyperedge_c.begin(), hyperedge_c.end());
      if (hyperedges_c[hash_id] == hyperedge_vec) {
        hyperedge_weights_c[hash_id]
            = hyperedge_weights_c[hash_id] + hypergraph->hyperedge_weights_[i];
      } else {
        hyperedge_weights_c.push_back(hypergraph->hyperedge_weights_[i]);
        hyperedges_c.push_back(hyperedge_vec);
      }
    }
  }
  std::shared_ptr<IlpGraph> ilpgraph
      = std::make_shared<IlpGraph>(hypergraph->vertex_dimensions_,
                                   hypergraph->hyperedge_dimensions_,
                                   true,
                                   fixed_vertices,
                                   hyperedges_c,
                                   vertex_weights_c,
                                   hyperedge_weights_c);
  return ilpgraph;
}
}  // namespace par
