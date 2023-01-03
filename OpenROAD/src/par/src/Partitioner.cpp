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
#include "utl/Logger.h"
#include "TPHypergraph.h"
#include "Utilities.h"
#include "Partitioner.h"

// for ILP solver in google OR-Tools
#include <absl/flags/flag.h>
#include <absl/strings/match.h>
#include <absl/strings/string_view.h>
#include <ortools/base/commandlineflags.h>
//#include <ortools/base/init_google.h>
#include <ortools/base/logging.h>
#include <ortools/linear_solver/linear_solver.h>
#include <ortools/linear_solver/linear_solver.pb.h>

// for ILP solver in CPLEX
#include "ilcplex/cplex.h"
#include "ilcplex/ilocplex.h"

using utl::PAR;

// par 2912


namespace par {

using operations_research::MPSolver;
using operations_research::MPObjective;
using operations_research::MPConstraint;
using operations_research::MPVariable;


// Golden Evaluator
std::pair<float, std::vector<std::vector<float> > >
  Partitioners::GoldenEvaluator(const HGraph hgraph,
                                std::vector<int> &solution,
                                bool print_flag)
{
  std::vector<std::vector<float> > block_balance = GetBlockBalance(hgraph, solution);
  float cost = 0.0;
  // check the cutsize
  for (int e = 0; e < hgraph->num_hyperedges_; e++) {
    for (int idx = hgraph->eptr_[e] + 1;  idx < hgraph->eptr_[e+1]; idx++) {
      if (solution[hgraph->eind_[idx]] != solution[hgraph->eind_[idx - 1]]) {
        cost += std::inner_product(hgraph->hyperedge_weights_[e].begin(),
                                    hgraph->hyperedge_weights_[e].end(),
                                    e_wt_factors_.begin(), 0.0);
        break; // this net has been cut
      }
    } // finish hyperedge e
  } 
  // check timing paths
  for (int path_id = 0; path_id < hgraph->num_timing_paths_; path_id++)
    cost += CalculatePathCost(path_id, hgraph, solution);

  // check if the solution is valid
  for (auto value : solution)
    if (value < 0 || value >= num_parts_)
      logger_->error(PAR, 2960, "The solution is invalid!!!");


  if (print_flag == true) {
    // print cost
    logger_->info(PAR, 2910, "Cutsize (cost): {:0.2f}", cost);
    // print block balance
    std::vector<float> tot_vertex_weights = hgraph->GetTotalVertexWeights();
    for (auto block_id = 0; block_id < block_balance.size(); block_id++) {
      std::string line = "Vertex balance of block_" + std::to_string(block_id) + " : ";
      for (auto dim = 0; dim < tot_vertex_weights.size(); dim++) {
        std::stringstream ss; // for converting float to string
        ss << std::fixed << std::setprecision(5) << 
                            block_balance[block_id][dim] / tot_vertex_weights[dim]
                          << "  ( "
                          << block_balance[block_id][dim] 
                          << " )  ";
        line += ss.str() + "  ";
      }
      logger_->info(PAR, 2911, line);
      std::cout << "cutsize : " << cost << std::endl;
      std::cout << "balance : " << line << std::endl;
    } // finish block balance
  }
  return std::pair<float, std::vector<std::vector<float> > >(cost, block_balance);
}
    
// Get block balance 
std::vector<std::vector<float> > 
    Partitioners::GetBlockBalance(const HGraph hgraph,
                                  std::vector<int> &solution)
{
  std::vector<std::vector<float> > 
     block_balance(num_parts_, std::vector<float>(hgraph->vertex_dimensions_, 0.0));
  if (hgraph->num_vertices_ != static_cast<int>(solution.size()))
    return block_balance;
  // check if the solution is valid
  for (auto block_id : solution)
    if (block_id < 0 || block_id >= num_parts_)
      return block_balance;
  // update the block_balance
  for (int v = 0; v < hgraph->num_vertices_; v++) 
    block_balance[solution[v]] = block_balance[solution[v]]
                               + hgraph->vertex_weights_[v];
  return block_balance;
}

/*
// Find start vertices for BFS based region growth
int Partitioners::FindStartVertices(const HGraph hgraph, std::vector<int>& solution) {
  int min_degree_vtx = -1;
  int min_degree = numeric_limits<int>::max();
  for(int i = 0; i < hgraph->num_vertices_; ++i) {
    if (solution[i] > -1) {
      continue;
    }
    int first_valid_entry = hgraph->vptr[i];
    int first_invalid_entry = hgraph->vptr[i+1];
    int degree = first_invalid_entry - first_valid_entry; 
    if (degree < min_degree) {
      min_degree = degree; 
      min_degree_vtx = i;
    }
  }
  return min_degree_vtx;
}

// Random BFS based initial partitioning
void Partitioners::BFSInitPartition(const HGraph hgraph, 
                                   const std::vector<std::vector<float> > &max_block_balance,
                                   std::vector<int> &solution) {
  std::vector<int> mark_vertices();
  int part_iter = 0;
  while(true) {
    std::priority_queue<int> pq;
    int start_vtx = FindStartingVertex(hgraph, solution);
    solution[start_vtx] = part_iter; 
    pq.push();
    for (int i = 0; i < hgraph->num_vertices_; ++i) {
      std::vector<int> border_vertices = FindBorderVertices(hgraph, solution, part_iter);
    }
    
  }



  std::vector<std::vector<int> > mark_vertices(num_parts_, std::vector<int> (hgraph->num_vertices_, 0));
  std::vector<std::vector<int> > mark_hyperedges(num_parts_, std::vector<int> (hgraph->num_hyperedges_, 0));
  std::vector<int> start_vertices(num_parts_, -1);
  FindStartVertices(hgraph, start_vertices);
  std::vector<std::priority_queue<int> > pqs;
  for (int i = 0; i < num_parts_; ++i) {
    std::priority_queue<int> pq;
    pq.push(start_vertices[i]);
    pqs.push_back(pq);
  }

}
*/

// Random Partition: 
// Randomly shuffle the vertices and assign it to a block
void Partitioners::RandomPartition(const HGraph hgraph, 
                                   const std::vector<std::vector<float> > &max_block_balance,
                                   std::vector<int> &solution)
{
  // reset variable
  solution.clear();
  solution.resize(hgraph->num_vertices_);
  std::fill(solution.begin(), solution.end(), -1);
  // the summation of vertex weights for vertices in current block
  std::vector<std::vector<float> > block_balance(num_parts_,
      std::vector<float>(hgraph->vertex_dimensions_, static_cast<float>(0)));
  // determine all the free vertices
  std::vector<int> unvisited;
  unvisited.reserve(hgraph->num_vertices_);
  if (hgraph->fixed_vertex_flag_ == true) {
    for (auto v = 0; v < hgraph->num_vertices_; v++) {
      if (hgraph->fixed_attr_[v] > -1) {
        solution[v] = hgraph->fixed_attr_[v];
        block_balance[solution[v]] = block_balance[solution[v]]
                                   + hgraph->vertex_weights_[v];
      } else {
        unvisited.push_back(v);  
      }
    } // done traversal
  } else {
    unvisited.resize(hgraph->num_vertices_);
    std::iota(unvisited.begin(), unvisited.end(), 0); // Fill with 0, 1, ...
  }
  // random shuffle
  shuffle(unvisited.begin(), unvisited.end(), std::default_random_engine(seed_));
  // assign vertex to blocks
  int block_id = 0;
  for (auto v : unvisited) {
    solution[v] = block_id;
    block_balance[solution[v]] = block_balance[solution[v]]
                               + hgraph->vertex_weights_[v];
    if (block_balance[block_id] >= DivideFactor(max_block_balance[block_id], 10.0))
      block_id = (++block_id) % num_parts_; // adjust the block_id
  }
}

void Partitioners::SetPartitionerSeed(int seed) {
  seed_ = seed;
}

int Partitioners::GetPartitionerSeed() const {
  return seed_;
}

// CPLEX with warm start
void Partitioners::OptimalInitialPartitionCplexWarmStart(const HGraph hgraph,
          const std::vector<std::vector<float> > &max_block_balance,
                                        std::vector<int> &solution) {
  std::vector<std::vector<int> > x(num_parts_, std::vector<int> (hgraph->num_vertices_, 0));
  std::vector<std::vector<int> > y(num_parts_, std::vector<int> (hgraph->num_hyperedges_, 0));
  for (int i = 0; i < hgraph->num_vertices_; ++i) {
    int p = solution[i];
    x[p][i] = 1;
  }
  for (int i = 0; i < hgraph->num_hyperedges_; ++i) {
    int firstValidEntry = hgraph->eptr_[i];
    int firstInvalidEntry = hgraph->eptr_[i+1];
    std::set<int> unique_partitions;
    for (int j = firstValidEntry; j < firstInvalidEntry; ++j) {
      int v = hgraph->eind_[j];
      int p = solution[v];
      unique_partitions.insert(p);
    }
    for (const int &j : unique_partitions) {
      y[j][i] = 1;
    }
  }
  IloEnv myenv;
  IloModel mymodel(myenv);
  IloNumVarArray startVarX(myenv);
  IloNumArray startValX(myenv);
  IloArray<IloNumVarArray> _x(myenv, num_parts_);
  IloArray<IloNumVarArray> _y(myenv, num_parts_);
  for (int i = 0; i < num_parts_; ++i) {
    _x[i] = IloNumVarArray(myenv, hgraph->num_vertices_, 0, 1, ILOINT);
    _y[i] = IloNumVarArray(myenv, hgraph->num_hyperedges_, 0, 1, ILOINT);
  }
  for (int i = 0; i < num_parts_; ++i) {
    for (int j = 0; j < hgraph->num_vertices_; ++j) {
      startVarX.add(_x[i][j]);
      startValX.add(x[i][j]);
    }
    /*for (int j = 0; j < hgraph->num_hyperedges_; ++j) {
      startVarX.add(_y[i][j]);
      startValX.add(y[i][j]);
      //startVarY.add(y[i][j]);
    }*/
  }

  // define constraints
  // balance constraint 
  // check each dimension  
  for (int i = 0; i < hgraph->vertex_dimensions_; ++i) {
    // allowed balance for each dimension
    for (int j = 0; j < num_parts_; ++j) {
      IloExpr balance_expr(myenv);
      for (int k = 0; k < hgraph->num_vertices_; ++k) {
        balance_expr += hgraph->vertex_weights_[k][i] * _x[j][k];
      } // finish traversing vertices
      mymodel.add(IloRange(myenv, 0.0, balance_expr, max_block_balance[j][i]));
      balance_expr.end();   
    } // finish traversing blocks
  } // finish dimension check
  // Fixed vertices constraint
  if (hgraph->fixed_vertex_flag_ == true) {
    for (int i = 0; i < hgraph->num_vertices_; ++i) {
      if (hgraph->fixed_attr_[i] > -1) {
        mymodel.add(_x[hgraph->fixed_attr_[i]][i] == 1);
      } // fixed vertices should be placed at the specified block
    }
  }
  // each vertex can only belong to one part
  for (int i = 0; i < hgraph->num_vertices_; ++i) {
    IloExpr vertex_expr(myenv);
    for (int j = 0; j < num_parts_; ++j) {
      vertex_expr += _x[j][i];
    }
    mymodel.add(vertex_expr == 1);
    vertex_expr.end();
  }
  // Hyperedge constraint
  for (int i = 0; i < hgraph->num_hyperedges_; ++i) {
    const int start_idx = hgraph->eptr_[i];
    const int end_idx = hgraph->eptr_[i + 1];
    for (int j = start_idx; j < end_idx; j++) {
      const int vertex_id = hgraph->eind_[j];
      for (int k = 0; k < num_parts_; ++k) {
        mymodel.add(_y[k][i] <= _x[k][vertex_id]);
      }
    }
  }
  // Maximize cutsize objective    
  IloExpr obj_expr(myenv);  // empty expression
  for (int i = 0; i < hgraph->num_hyperedges_; ++i) {
    const float cost_value = std::inner_product(
        hgraph->hyperedge_weights_[i].begin(),
        hgraph->hyperedge_weights_[i].end(),
                            e_wt_factors_.begin(), 0.0);
    for (int j = 0; j < num_parts_; ++j) {
      obj_expr += cost_value * _y[j][i];
    }
  }  
  mymodel.add(IloMaximize(myenv, obj_expr));  // adding minimization objective
  obj_expr.end();                             // clear memory
  // Model Solution
  IloCplex mycplex(myenv);
  mycplex.extract(mymodel);
  //mycplex.setParam(IloCplex::Param::Preprocessing::Presolve, 0);
  //mycplex.setParam(IloCplex::Param::Read::Scale, 1);
  //mycplex.setParam(IloCplex::Param::Simplex::Tolerances::Markowitz, 0.99);
  mycplex.setOut(myenv.getNullStream());
  mycplex.addMIPStart(startVarX, startValX, IloCplex::MIPStartAuto,  "secondMIPStart");
  mycplex.setParam(IloCplex::Param::MIP::Limits::Solutions, 2);
  startVarX.end();
  startValX.end();
  mycplex.solve();
  IloBool feasible = mycplex.solve();
  if (feasible == IloTrue) {
    for (int i = 0; i < hgraph->num_vertices_; i++) {
      for (int j = 0; j < num_parts_; j++) {
        if (mycplex.getValue(_x[j][i]) == 1.00) {
          solution[i] = j;
        }
      }
    }
    // some solution may invalid due to the limitation of ILP solver
    for (auto& value : solution)
      value = (value == -1) ? 0 : value;
  } else {
    std::cout << "No feasible solution found with ILP --> Running K-way FM instead" << std::endl;
    DirectKWayFM(hgraph, max_block_balance, solution);
  }
  // closing the model
  mycplex.clear();
  myenv.end(); 
}

// CPLEX based ILP Solver
void Partitioners::OptimalInitialPartitionCplex(const HGraph hgraph,
          const std::vector<std::vector<float> > &max_block_balance,
                                         std::vector<int> &solution) 
{
  // reset variable
  solution.clear();
  solution.resize(hgraph->num_vertices_);
  std::fill(solution.begin(), solution.end(), -1); 
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
  struct comp {
    // comparator function
    bool operator()(const std::pair<int, float> &l,
                    const std::pair<int, float> &r) const
    {
      if (l.second != r.second)
        return l.second > r.second;
      return l.first < r.first;
    }
  };

  // use set data structure to sort unvisited vertices
  std::set<std::pair<int, float>, comp> unvisited_hyperedges; 
  for (auto e = 0; e < hgraph->num_hyperedges_; e++) {
    const float score = std::inner_product(hgraph->hyperedge_weights_[e].begin(),
                                           hgraph->hyperedge_weights_[e].end(),
                                           e_wt_factors_.begin(), 0.0);
    unvisited_hyperedges.insert(std::pair<int, float>(e, score));
  }

  /*
  const int max_num_hyperedges = 50;
  std::vector<int> edge_mask;
  int e  = 0;
  for (auto& value : unvisited_hyperedges) {
    std::cout << "weight = " << value.second << std::endl;
    edge_mask.push_back(value.first);
    e++;
    if (e > max_num_hyperedges)
      break;
  }
  */
  int max_num_hyperedges = 0;
  std::vector<int> edge_mask;
  float base_score = (*unvisited_hyperedges.begin()).second / 10.0;
  for (auto& value : unvisited_hyperedges) {
    /*if (value.second <= base_score || max_num_hyperedges >= 50)
      break;*/
      
    //std::cout << "weight = " << value.second << std::endl;
    edge_mask.push_back(value.first);
    max_num_hyperedges++;
  }
 

  for (int i = 0; i < num_parts_; i++) {
    x[i] = IloNumVarArray(myenv, hgraph->num_vertices_, 0, 1, ILOINT);
    y[i] = IloNumVarArray(myenv, edge_mask.size(), 0, 1, ILOINT);
  }

  // define constraints
  // balance constraint 
  // check each dimension  
  for (int i = 0; i < hgraph->vertex_dimensions_; i++) {
    // allowed balance for each dimension
    for (int j = 0; j < num_parts_; j++) {
      IloExpr balance_expr(myenv);
      for (int k = 0; k < hgraph->num_vertices_; k++) {
        balance_expr += hgraph->vertex_weights_[k][i] * x[j][k];
      } // finish traversing vertices
      mymodel.add(IloRange(myenv, 0.0, balance_expr, max_block_balance[j][i]));
      balance_expr.end();   
    } // finish traversing blocks
  } // finish dimension check
  // Fixed vertices constraint
  if (hgraph->fixed_vertex_flag_ == true) {
    for (int i = 0; i < hgraph->num_vertices_; i++) {
      if (hgraph->fixed_attr_[i] > -1) {
        mymodel.add(x[hgraph->fixed_attr_[i]][i] == 1);
      } // fixed vertices should be placed at the specified block
    }
  }
  // each vertex can only belong to one part
  for (int i = 0; i < hgraph->num_vertices_; i++) {
    IloExpr vertex_expr(myenv);
    for (int j = 0; j < num_parts_; j++) {
      vertex_expr += x[j][i];
    }
    mymodel.add(vertex_expr == 1);
    vertex_expr.end();
  }
  // Hyperedge constraint
  for (int i = 0; i < edge_mask.size(); i++) {
    const int e = edge_mask[i];
    const int start_idx = hgraph->eptr_[e];
    const int end_idx = hgraph->eptr_[e + 1];
    for (int j = start_idx; j < end_idx; j++) {
      const int vertex_id = hgraph->eind_[j];
      for (int k = 0; k < num_parts_; k++) {
        mymodel.add(y[k][i] <= x[k][vertex_id]);
      }
    }
  }
  // Maximize cutsize objective    
  IloExpr obj_expr(myenv);  // empty expression
  for (int i = 0; i < edge_mask.size(); i++) {
    const float cost_value = std::inner_product(
        hgraph->hyperedge_weights_[edge_mask[i]].begin(),
        hgraph->hyperedge_weights_[edge_mask[i]].end(),
                            e_wt_factors_.begin(), 0.0);
    for (int j = 0; j < num_parts_; j++) {
      obj_expr += cost_value * y[j][i];
    }
  }  
  mymodel.add(IloMaximize(myenv, obj_expr));  // adding minimization objective
  obj_expr.end();                             // clear memory
  // Model Solution
  IloCplex mycplex(myenv);
  mycplex.extract(mymodel);
  mycplex.solve();
  IloBool feasible = mycplex.solve();
  if (feasible == IloTrue) {
    for (int i = 0; i < hgraph->num_vertices_; i++) {
      for (int j = 0; j < num_parts_; j++) {
        if (mycplex.getValue(x[j][i]) == 1.00) {
          solution[i] = j;
        }
      }
    }
    // some solution may invalid due to the limitation of ILP solver
    for (auto& value : solution)
      value = (value == -1) ? 0 : value;
  } else {
    logger_->info(PAR, 2119, "There is no feasible solution in ILP. Call FM instead");
    DirectKWayFM(hgraph, max_block_balance, solution);
  }
  // closing the model
  mycplex.clear();
  myenv.end(); 
}


// Google OR-Tools based ILP Solver
void Partitioners::OptimalInitialPartition(const HGraph hgraph,
     const std::vector<std::vector<float> > &max_block_balance,
                                    std::vector<int> &solution) 
{  
  // reset variable
  solution.clear();
  solution.resize(hgraph->num_vertices_);
  std::fill(solution.begin(), solution.end(), -1); 
  // Google OR-Tools Implementation
  std::unique_ptr<MPSolver> solver(MPSolver::CreateSolver("SCIP"));
  if (!solver) {
    logger_->error(PAR, 2504, "SCIP solver unavailable!");
  }

  // Define constraints
  // For each vertex, define a variable x
  // For each hyperedge, define a variable y 
  std::vector<std::vector<const MPVariable*> > x(num_parts_,
              std::vector<const MPVariable*>(hgraph->num_vertices_));
  std::vector<std::vector<const MPVariable*> > y(num_parts_,
              std::vector<const MPVariable*>(hgraph->num_hyperedges_));
  // initialize variables
  for (auto& x_v_vector : x)
    for (auto& x_v : x_v_vector)
      x_v = solver->MakeIntVar(0.0, 1.0, ""); // represent whether the vertex is within block
  for (auto& y_e_vector : y)
    for (auto& y_e : y_e_vector)
      y_e = solver->MakeIntVar(0.0, 1.0, ""); // represent whether the hyperedge is within block
  // check number of variables
  logger_->info(PAR, 2505, "Number of variables = {}", solver->NumVariables());
  // define the inifity constant
  const double infinity = solver->infinity();
  // define constraints
  // balance constraint 
  // check each dimension
  for (int i = 0; i < hgraph->vertex_dimensions_; i++) {
    // allowed balance for each dimension
    for (int j = 0; j < num_parts_; j++) {        
      MPConstraint* constraint = solver->MakeRowConstraint(0.0, max_block_balance[j][i], "");
      for (int k = 0; k < hgraph->num_vertices_; k++) {
        constraint->SetCoefficient(x[j][k], hgraph->vertex_weights_[k][i]);
      } // finish travering vertices
    } // finish traversing blocks
  }
  // Fixed vertices constraint
  if (hgraph->fixed_vertex_flag_ == true) {
    for (int i = 0; i < hgraph->num_vertices_; i++) {
      if (hgraph->fixed_attr_[i] > -1) {
        MPConstraint* constraint = solver->MakeRowConstraint(1, 1, "");
        constraint->SetCoefficient(x[hgraph->fixed_attr_[i]][i], 1);
      } // fixed vertices should be placed at the specified block
    }
  }
  // each vertex can only belong to one part
  for (int i = 0; i < hgraph->num_vertices_; i++) {
    MPConstraint* constraint = solver->MakeRowConstraint(1, 1, "");
    for (int j = 0; j < num_parts_; j++) {
      constraint->SetCoefficient(x[j][i], 1);
    }
  }
  // Hyperedge constraint
  for (int i = 0; i < hgraph->num_hyperedges_; i++) {
    const int start_idx = hgraph->eptr_[i];
    const int end_idx = hgraph->eptr_[i + 1];
    for (int j = start_idx; j < end_idx; j++) {
      const int vertex_id = hgraph->eind_[j];
      for (int k = 0; k < num_parts_; k++) {
        MPConstraint* constraint = solver->MakeRowConstraint(0, infinity, "");
        constraint->SetCoefficient(x[k][vertex_id], 1);
        constraint->SetCoefficient(y[k][i], -1);
      }
    }
  }
  // check number of constraints
  logger_->info(PAR, 2506, "Number of constraints = {}",
                           solver->NumConstraints());
  // Maximize cutsize objective
  MPObjective* const obj_expr = solver->MutableObjective();
  for (int i = 0; i < hgraph->num_hyperedges_; i++) {
    const float cost_value = std::inner_product(hgraph->hyperedge_weights_[i].begin(),
                                                hgraph->hyperedge_weights_[i].end(),
                                                e_wt_factors_.begin(), 0.0);
    for (int j = 0; j < num_parts_; j++) {
      obj_expr->SetCoefficient(y[j][i], cost_value);
    }
  }
  obj_expr->SetMaximization();

  // Solve the ILP Problem
  const MPSolver::ResultStatus result_status = solver->Solve();
  // Check that the problem has an optimal solution.
  if (result_status == MPSolver::OPTIMAL) {
    for (int i = 0; i < hgraph->num_vertices_; i++) {
      for (int j = 0; j < num_parts_; j++) {
        if (x[j][i]->solution_value() == 1.0) {
          solution[i] = j;
        }
      }  
    }
    // debug, print solutions
    //for (auto value : solution)
    //  std::cout << value << std::endl;

  } else {
    logger_->info(PAR, 2507, "There is no feasible solution in ILP. Call FM instead");
    DirectKWayFM(hgraph, max_block_balance, solution);
  }
}

// Direct Priority-queue FM-based partitioning
// The input argument - solution can be used as the initial solution
// and it's also the improved the solution
void Partitioners::DirectKWayFM(const HGraph hgraph,
                                const std::vector<std::vector<float> > &max_block_balance,
                                std::vector<int> &solution)
{
  // If there is no initial solution, use random initial solution
  if (solution.size() < hgraph->vertex_weights_.size()) {
    std::cout <<" Trying random partition" << std::endl;
    RandomPartition(hgraph, max_block_balance, solution);
  }
  for (int num_pass = 0; num_pass < max_num_fm_pass_; num_pass++) {
    //std::cout << "\n\nnum_pass = " << num_pass << std::endl;
    GoldenEvaluator(hgraph, solution, false);
    //std::vector<int> pre_fm_solution = solution;
    const float tot_gain = DirectKWayFMPass(hgraph, max_block_balance, solution);
    //std::cout << "tot_gain = " << tot_gain << std::endl;
    //std::cout << std::endl;
    if (tot_gain < 0.0) {
      //solution = pre_fm_solution;
      break; // stop pass loop if there is no improvement
    }
  }
}

// Direct Priority-queue FM-based partitioning
// The input argument - solution can be used as the initial solution
// and it's also the improved the solution
void Partitioners::DirectKWayFMWithImb(const HGraph hgraph,
                                const std::vector<std::vector<float> > &max_block_balance,
                                std::vector<int> &solution)
{
  // If there is no initial solution, use random initial solution
  if (solution.size() < hgraph->vertex_weights_.size()) {
    std::cout <<" Trying random partition" << std::endl;
    RandomPartition(hgraph, max_block_balance, solution);
  }
  for (int num_pass = 0; num_pass < max_num_fm_pass_; num_pass++) {
    //std::cout << "\n\nnum_pass = " << num_pass << std::endl;
    GoldenEvaluator(hgraph, solution, false);
    //std::vector<int> pre_fm_solution = solution;
    std::vector<std::vector<float> > temp_block_balance = max_block_balance;
    if (num_pass == 0) {
      for (int i = 0; i < num_parts_; ++i) {
        temp_block_balance[i] = max_block_balance[i] * 1.5; // allowing a 5% imbalance during 1st pass
      }
    } else {
      std::vector<std::vector<float> > curr_balance = GetBlockBalance(hgraph, solution);
      bool balance_check = true; 
      for (int i = 0; i < num_parts_; ++i) {
        if (curr_balance[i] > max_block_balance[i]) {
          balance_check = false;
          break;
        }
      }
      if (balance_check == false) {
        std::cout << "[FM] Balancing the partition " << std::endl;
        BalancePartition(hgraph, max_block_balance, solution);
        GoldenEvaluator(hgraph, solution, true);
      }
    }
    const float tot_gain = DirectKWayFMPass(hgraph, temp_block_balance, solution);
    if (tot_gain < 0.0) {
      //solution = pre_fm_solution;
      break; // stop pass loop if there is no improvement
    }
  }
}

//Partition Balancing using FM
void Partitioners::BalancePartition(const HGraph hgraph, 
                                     const std::vector<std::vector<float> > &max_block_balance,
                                     std::vector<int> &solution) {
  std::vector<float> total_vertex_wts = hgraph->GetTotalVertexWeights();
  std::vector<std::vector<int> > net_degs = GetNetDegrees(hgraph, solution);
  std::vector<std::vector<float> > block_balance = GetBlockBalance(hgraph, solution);
  bool early_balance_check = true; 
  for (int i = 0; i < num_parts_; ++i) {
    if (block_balance[i] > max_block_balance[i]) {
      early_balance_check = false; 
      break;
    }
  }
  if (early_balance_check == true) {
    std::cout << "[BALANCE] Partition is already balanced! Exiting!" << std::endl;
    return; 
  }
  int num_visited_vertex = 0; 
  std::vector<bool> visited(hgraph->num_vertices_, false);  
  if (hgraph->fixed_vertex_flag_ == true) {
    for (auto v = 0; v < hgraph->num_vertices_; v++) {
      if (hgraph->fixed_attr_[v] > -1) {
        visited[v] = true;
        num_visited_vertex++;
      }      
    }
  }
  // Initialize cost for timing paths
  std::vector<float> paths_cost;
  paths_cost.reserve(hgraph->num_timing_paths_);
  for (int path_id = 0; path_id < hgraph->num_timing_paths_; path_id++) {
    paths_cost.push_back(CalculatePathCost(path_id, hgraph, solution));
  }
  std::vector<std::set<VertexGain> > gain_buckets(num_parts_);
  std::vector<std::thread> threads; // for parallel updating
  // parallel initialize the num_parts gain_buckets
  for (int to_pid = 0; to_pid < num_parts_; to_pid++) {
    threads.push_back(std::thread(&Partitioners::InitializeGainBucket, this, 
                                  &gain_buckets[to_pid], to_pid, &visited,
                                  hgraph, &solution, &paths_cost, &net_degs));
  }
  for (auto& t : threads) {
    t.join(); // wait for all threads to finish
  }
  threads.clear();
  size_t tot_bucket_size = 0;
  for (auto pid = 0; pid < num_parts_; pid++) {
    tot_bucket_size += gain_buckets[pid].size();
  }
  float tot_gain = 0.0;
  while (true) {
    bool stop_flag = true;
    int to_pid = -1;
    for (int i = 0; i < num_parts_; ++i) {
      if (block_balance[i] > max_block_balance[i]) {
        stop_flag = false;
        break;
      }
    }
    if (stop_flag == false) {
      std::vector<float> min_bal = total_vertex_wts;
      for (int i = 0; i < num_parts_; ++i) {
        if (block_balance[i] < min_bal) {
          min_bal = block_balance[i];
          to_pid = i;
        }
      }
    } else {
      break;
    }
    for (auto p_id = 0; p_id < num_parts_; ++p_id) {
      if (gain_buckets[p_id].size() == 0) {
        stop_flag = true;
        break;
      }
    }
    if (stop_flag == true) {
      break;
    }
    // Pick best vertex from side with least wt
    const VertexGain& candidate = *(gain_buckets[to_pid].cbegin());
    const int v = candidate.vertex; // candidate vertex
    tot_gain += candidate.gain;
    //std::cout << "[DEBUG] Total gain = " << tot_gain << std::endl;
    num_visited_vertex++;
    visited[v] = true;
    int from_pid = solution[v];
    solution[v] = to_pid;
    //GoldenEvaluator(hgraph, solution, true);
    // update path cost
    for (const auto& [path_id, cost] : candidate.path_cost)
      paths_cost[path_id] = cost;
    // update net_degs
    for (auto e_idx = hgraph->vptr_[v]; e_idx < hgraph->vptr_[v+1]; e_idx++) {
      const int e = hgraph->vind_[e_idx];
      net_degs[e][from_pid]--;
      net_degs[e][to_pid]++;
    }
    // update balance
    block_balance[from_pid]  = block_balance[from_pid]
                                - hgraph->vertex_weights_[v];
    block_balance[to_pid] = block_balance[to_pid]
                                      + hgraph->vertex_weights_[v];
    // update the neighbors of v for all gain buckets in parallel
    for (int j = 0; j < num_parts_; ++j)
      threads.push_back(std::thread(&Partitioners::UpdateGainBucket, this, v,
                            from_pid, j, &gain_buckets[j], hgraph, 
                            &visited, &solution,  &paths_cost, &net_degs));
    for (auto& t : threads)
      t.join(); // wait for all threads to finish
    threads.clear();
  }
}


// Direct Priority-queue FM-based partitioning
// The input argument - solution can be used as the initial solution
// and it's also the improved the solution
float Partitioners::DirectKWayFMPass(const HGraph hgraph, 
                                     const std::vector<std::vector<float> > &max_block_balance,
                                     std::vector<int> &solution)
{
  // If there is no initial solution, use random initial solution
  if (solution.size() < hgraph->vertex_weights_.size()) {
    std::cout << "Generating random partition because of no initial solution" << std::endl;
    RandomPartition(hgraph, max_block_balance, solution);
  }
  // record the moves in this pass
  // global_best_ver_gain stores the best gain after moving vertex (used as timestamp)
  // update_solution stores the current solution after moving vertex v
  VertexGain global_best_ver_gain(-1, -std::numeric_limits<float>::max());
  std::vector<int> update_solution = solution;
  std::vector<int> move_trace; // store the moved vertices in sequence
  float tot_gain = 0.0;

  // Initialize cost for timing paths
  std::vector<float> paths_cost;
  paths_cost.reserve(hgraph->num_timing_paths_);
  for (int path_id = 0; path_id < hgraph->num_timing_paths_; path_id++)
    paths_cost.push_back(CalculatePathCost(path_id, hgraph, solution));
  // update the net degree for existing solution
  // for each hyperedge, calculate the number of vertices in each part
  std::vector<std::vector<int> > net_degs = GetNetDegrees(hgraph, solution);
  // calculate the current balance for each block
  std::vector<std::vector<float> > block_balance = GetBlockBalance(hgraph, solution);
  
  // Checking fixed vertices
  // stop the FM pass if the number of visited equal to or more than num_early_stop_vertex
  int num_visited_vertex = 0; // record the number of vertices visited
  std::vector<bool> visited(hgraph->num_vertices_, false);  // if the vertex has been locked 

  if (hgraph->fixed_vertex_flag_ == true) {
    for (auto v = 0; v < hgraph->num_vertices_; v++) {
      if (hgraph->fixed_attr_[v] > -1) {
        visited[v] = true;
        num_visited_vertex++;
      }      
    }
  }

  // stop the FM pass if the number of visited equal to or more than num_early_stop_vertex
  early_stop_ratio_ = 1.0;
  const int num_early_stop_vertex = std::min(max_num_moves_,
     num_visited_vertex + (hgraph->num_vertices_ - num_visited_vertex) * early_stop_ratio_);

  // Initialize current gain in a multi-thread manner
  // set based on max heap (k set)
  // each block has its own max heap
  std::vector<std::set<VertexGain> > gain_buckets(num_parts_);
  std::vector<std::thread> threads; // for parallel updating
  // parallel initialize the num_parts gain_buckets
  for (int to_pid = 0; to_pid < num_parts_; to_pid++)
    threads.push_back(std::thread(&Partitioners::InitializeGainBucket, this, 
                                  &gain_buckets[to_pid], to_pid, &visited,
                                  hgraph, &solution, &paths_cost, &net_degs));
  for (auto& t : threads)
    t.join(); // wait for all threads to finish
  threads.clear();
  
  size_t tot_bucket_size = 0;
  for (auto pid = 0; pid < num_parts_; pid++)
    tot_bucket_size += gain_buckets[pid].size();

  // Main loop of FM pass
  while (num_visited_vertex < num_early_stop_vertex) {
    // Find the candidate vertex v to move
    // candidate
    int to_pid = -1;
    VertexGain candidate(-1, -std::numeric_limits<float>::max()); // best_candidate
    // best gain bucket for "corking effect"
    int best_to_pid = -1; // block id with best_gain
    float best_gain = -std::numeric_limits<float>::max();
    bool stop_flag = false;
    // check if empty
    for (auto p_id = 0; p_id < num_parts_; p_id++) {
      if (gain_buckets[p_id].size() == 0) {
        stop_flag = true;
        break;
      }
    }
    if (stop_flag == true)
      break;
    
    for (auto p_id = 0; p_id < num_parts_; p_id++) {
      const VertexGain& vertex_gain = *(gain_buckets[p_id].cbegin());
      // check if there is feasible move
      if ((vertex_gain.gain > candidate.gain) && 
          (block_balance[p_id] + hgraph->vertex_weights_[vertex_gain.vertex]
           < max_block_balance[p_id])) {
        to_pid = p_id;
        candidate = vertex_gain;
      }
      // for the case of "corking effect"
      if (vertex_gain.gain > best_gain) {
        best_gain = vertex_gain.gain;
        best_to_pid = p_id;
      }
    }
    // checking if "corking effect", i.e., no candidate
    if (to_pid == -1) {
      // traverse the gain_buckets[to_pid] to find a candidate
      for (auto iter = gain_buckets[best_to_pid].begin();
              iter != gain_buckets[best_to_pid].end(); iter++) {
        const std::vector<float> update_block_balance = block_balance[best_to_pid]
                                                      + hgraph->vertex_weights_[iter->vertex];
        if (update_block_balance < max_block_balance[best_to_pid]) { 
          to_pid = best_to_pid;
          candidate = *iter;
          break;
        }
      }
      // if there is no feasible solution, stop current path
      if (to_pid == -1)
        break; 
    }

    // update the state of the partitioning
    const int v = candidate.vertex; // candidate vertex
    tot_gain += candidate.gain;
    move_trace.push_back(v);
    num_visited_vertex++;
    visited[v] = true;
    update_solution[v] = to_pid;
    if (global_best_ver_gain.gain <= tot_gain) {
      global_best_ver_gain.vertex = candidate.vertex;
      global_best_ver_gain.gain = tot_gain;
    }
    // update path cost
    for (const auto& [path_id, cost] : candidate.path_cost)
      paths_cost[path_id] = cost;
    // update net_degs
    for (auto e_idx = hgraph->vptr_[v]; e_idx < hgraph->vptr_[v+1]; e_idx++) {
      const int e = hgraph->vind_[e_idx];
      net_degs[e][solution[v]]--;
      net_degs[e][update_solution[v]]++;
    }

    // update balance
    block_balance[solution[v]] = block_balance[solution[v]]
                                - hgraph->vertex_weights_[v];
    block_balance[update_solution[v]] = block_balance[update_solution[v]]
                                      + hgraph->vertex_weights_[v];
    // update the neighbors of v for all gain buckets in parallel
    for (int to_pid = 0; to_pid < num_parts_; to_pid++)
      threads.push_back(std::thread(&Partitioners::UpdateGainBucket, this, v,
                            solution[v], to_pid, &gain_buckets[to_pid], hgraph, 
                            &visited, &update_solution,  &paths_cost, &net_degs));
    for (auto& t : threads)
      t.join(); // wait for all threads to finish
    threads.clear();
  }

  if (global_best_ver_gain.gain < 0.0)
    return global_best_ver_gain.gain; // if there is no improvement,  return original solution

  // update the solution
  for (auto v : move_trace) {
    solution[v] = update_solution[v];
    if (v == global_best_ver_gain.vertex) 
      break; // stop updating, we are at the global optimum for current pass
  } // finish solution update
 
  return global_best_ver_gain.gain;
}

// We enable multithread for initializing gain bucket
void Partitioners::InitializeGainBucket(std::set<VertexGain> *gain_bucket_ptr, 
                                        int to_pid,
                                        const std::vector<bool> *visited_ptr,
                                        const HGraph hgraph,
                                        const std::vector<int> *solution_ptr,
                                        const std::vector<float> *cur_path_cost_ptr,
                                        const std::vector<std::vector<int> > *net_degs_ptr)
{
  std::set<VertexGain> &gain_bucket = *gain_bucket_ptr;
  const std::vector<bool> &visited = *visited_ptr;
  const std::vector<int> &solution = *solution_ptr;
  const std::vector<float> &cur_path_cost = *cur_path_cost_ptr;
  const std::vector<std::vector<int> > &net_degs = *net_degs_ptr;
  int idx = 0;
  for (auto v = 0; v < hgraph->num_vertices_; v++) {
    if (visited[v] == true || solution[v] == to_pid)
      continue;
    idx++;
    gain_bucket_ptr->insert(CalculateGain(v, solution[v], to_pid, hgraph, 
                                     solution, cur_path_cost, net_degs));
  }
}

// We enable multithread for updating gain bucket after moving vertex v
void Partitioners::UpdateGainBucket(int v, int from_pid, int to_pid,
                                    std::set<VertexGain> *gain_bucket_ptr, 
                                    const HGraph hgraph,
                                    const std::vector<bool> *unvisited_ptr,
                                    const std::vector<int>  *solution_ptr,
                                    const std::vector<float> *cur_path_cost_ptr,
                                    const std::vector<std::vector<int> > *net_degs_ptr)
{
  std::set<VertexGain> &gain_bucket = *gain_bucket_ptr;
  const std::vector<bool> &unvisited = *unvisited_ptr;
  const std::vector<int> &solution = *solution_ptr;
  const std::vector<float> &cur_path_cost = *cur_path_cost_ptr;
  const std::vector<std::vector<int> > &net_degs = *net_degs_ptr;
    
  // find all the neighbors of vertex v first
  std::set<int> neighbors;
  for (auto e_idx = hgraph->vptr_[v]; e_idx < hgraph->vptr_[v+1]; e_idx++) {
    const int e = hgraph->vind_[e_idx];
    for (auto v_idx = hgraph->eptr_[e]; v_idx < hgraph->eptr_[e+1]; v_idx++) {
      const int v = hgraph->eind_[v_idx];
      if (unvisited[v] == false)
        neighbors.insert(v);
    }
  }
  // If the current block_id is not equal to from_pid
  // we need to remove v out from the current gain_bucket
  if (from_pid != to_pid)
    neighbors.insert(v);  
  // find the position of each neighbor vertices
  std::vector<std::set<VertexGain>::iterator> iters;
  for (auto iter = gain_bucket.begin(); iter != gain_bucket.end(); iter++) {
    if (neighbors.find(iter->vertex) != neighbors.end())
      iters.push_back(iter);
    if (iters.size() >= neighbors.size())
      break; // find all the vertices
  }
  // Check the updated_gain for neighbors
  for (auto iter : iters) {
    // remove vertex v from the heap
    if (iter->vertex == v) {
      gain_bucket.erase(iter);
      continue;
    }
    // calculate the vertex_gain for the neighbor
    VertexGain vertex_gain = CalculateGain(iter->vertex,  
                                  solution[iter->vertex], to_pid,
                                  hgraph, solution, cur_path_cost, net_degs);
    // if the vertex_gain has been changed, then update the heap
    if (((*iter) == vertex_gain) == false) {
      gain_bucket.erase(iter);
      gain_bucket.insert(vertex_gain);
    }
  }
}



/*
// Improve the final cutsize based on greedy hyperedge
void Partitioners::GreedyRemoveHyperedge(const HGraph hgraph,
          std::vector<std::vector<float> > max_block_balance,
                                 std::vector<int> &solution)
{
  std::vector<std::vector<float> > block_balances = 
          (GoldenEvaluator(hgraph, solution, false)).second;
  // balance constraints
  std::vector<std::vector<float> > reduce_balance(num_parts_);
  for (int i = 0; i < num_parts_; i++)
    reduce_balance[i] = max_block_balance[i] - block_balances[i];
  int num_stop_hyperedge = 200;
  // check hyperedge
  for (int i = 0; i < hgraph->num_hyperedges_; i++) {
    std::vector<std::vector<int> > blocks(num_parts_);
    for (int idx = hgraph->eptr_[i]; idx < hgraph->eptr_[i+1]; idx++) {
      const int vertex_id = hgraph->eind_[idx];
      blocks[solution[vertex_id]].push_back(vertex_id);
    }
    int num_blocks = 0;
    for (auto& block : blocks)
      if (block.size() > 0)
        num_blocks++;
    if (num_blocks == 1)
      continue;
    // find the mininum balance
    int min_block_id = -1;
  }
}
*/


// update the net degree for existing solution
// for each hyperedge, calculate the number of vertices in each part
std::vector<std::vector<int> > Partitioners::GetNetDegrees(const HGraph hgraph,
                                                   std::vector<int> &solution) 
{
  std::vector<std::vector<int> > net_degs(hgraph->num_hyperedges_,
                                 std::vector<int>(num_parts_, 0));
  for (int e = 0; e < hgraph->num_hyperedges_; e++) {
    for (int idx = hgraph->eptr_[e]; idx < hgraph->eptr_[e+1]; idx++) {
      net_degs[e][solution[hgraph->eind_[idx]]]++;
    } 
  } 
  return net_degs;
}


// Calculate the cost for each timing path
// In the default mode (v = -1, to_pid = -1), 
// we just calculate the cost for the timing path path_id
// In the replacement mode, we replace the block id of v to to_pid
float Partitioners::CalculatePathCost(int path_id,
                                      const HGraph hgraph,
                                      const std::vector<int> &solution,
                                      int v,
                                      int to_pid)
{    
  float cost = 0.0; // cost for current path
  if (hgraph->num_timing_paths_ == 0)
    return cost; // no timing paths
  
  std::vector<int> path; // represent the path in terms of block_id
  std::map<int, int> block_counter; // block_id counter
  for (auto idx = hgraph->vptr_p_[path_id]; 
            idx < hgraph->vptr_p_[path_id]; idx++) {
    const int u = hgraph->vind_p_[idx]; // current vertex
    const int block_id = (u == v) ? to_pid : solution[u];
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

// Calculate the gain for a vertex v
// from_pid is the id of current block
// to_pid is the id of destination block
Partitioners::VertexGain Partitioners::CalculateGain(int v, int from_pid, int to_pid,
                                       const HGraph hgraph,
                                       const std::vector<int> &solution,
                                       const std::vector<float> &cur_path_cost,
                                       const std::vector<std::vector<int> > &net_degs)
{
  float score = 0.0;
  std::map<int, float> path_cost; // map path_id to latest score
  if (from_pid == to_pid)
    return VertexGain(v, score, path_cost); 

  // define two lamda function
  // 1) for checking connectivity
  // 2) for calculating the score of a hyperedge
  // function : check the connectivity for the hyperedge
  auto GetConnectivity = [&] (int e) {
    int connectivity = 0;
    for (auto& num_v : net_degs[e])
      if (num_v > 0)
        connectivity++;
    return connectivity;
  };
  
  // function : calculate the score for the hyperedge
  auto GetHyperedgeScore = [&] (int e)
  {
    return std::inner_product(hgraph->hyperedge_weights_[e].begin(),
                              hgraph->hyperedge_weights_[e].end(),
                              e_wt_factors_.begin(), 0.0);
  };

  // traverse all the hyperedges connected to v
  for (auto e_idx = hgraph->vptr_[v]; e_idx < hgraph->vptr_[v+1]; e_idx++) {
    const int e = hgraph->vind_[e_idx]; // hyperedge id
    const int connectivity = GetConnectivity(e);
    const float e_score = GetHyperedgeScore(e);
    if (connectivity == 1) { // move from_pid to to_pid will have negative socre
      score -= e_score;
    } else if (connectivity == 2 && net_degs[e][from_pid] == 1 && net_degs[e][to_pid] > 0) {
      score += e_score; // after move, all the vertices in to_pid
    }
  }

  // check the timing path
  if (hgraph->num_timing_paths_ > 0) {
    for (auto p_idx = hgraph->pptr_v_[v]; p_idx < hgraph->pptr_v_[v+1]; p_idx++) {
      const int path_id = hgraph->pind_v_[p_idx];
      const float cost = CalculatePathCost(path_id, hgraph, solution, v, to_pid);
      path_cost[path_id]  = cost;
      score += cur_path_cost[path_id] - cost;
    }
  }

  return VertexGain(v, score, path_cost);
}


/*
float Partitioners::KPMcalculateSpan(const HGraph hgraph,
                          int& from_pid, 
                          int& to_pid,
                          std::vector<int>& solution) {
  float span = 0.0;
  for (int i = 0; i < hgraph->num_hyperedges_; ++i) {
    const int first_valid_entry = hgraph->eptr_[i];
    const int first_invalid_entry = hgraph->eptr_[i+1];
    bool flag_partition_from = false; 
    bool flag_partition_to = false; 
    for (int j = first_valid_entry; j < first_invalid_entry; ++j) {
      const int v_id = hgraph->eind_[j];
      if (solution[v_id] == from_pid) {
        flag_partition_from = true; 
      } else if (solution[v_id] == to_pid) {
        flag_partition_to = true;
      }
    }
    if (flag_partition_from == true && flag_partition_to == true) {
      span += std::inner_product(hgraph->hyperedge_weights_[i].begin(),
                                hgraph->hyperedge_weights_[i].end(),
                                e_wt_factors_.begin(), 0.0);
    }
  }
  return span; 
}

std::vector<int> Partitioners::KPMfindBoundaryVertices(const HGraph hgraph,
                                            int& from_pid, 
                                            int& to_pid,
                                            std::vector<std::vector<int> > net_degs) {
  std::set<int> boundary_set; 
  std::vector<bool> boundary_net_flag(hgraph->num_hyperedges_, false); 
  for (int i = 0; i < hgraph->num_hyperedges_; ++i) {
    if (net_degs[i][from_pid] > 0 && net_degs[i][to_pid] > 0) {
      boundary_net_flag[i] = true;
    }
  }
  for (int i = 0; i < hgraph->num_vertices_; ++i) {
    const int first_valid_entry = hgraph->vptr_[i];
    const int first_invalid_entry = hgraph->vptr_[i+1];
    for (int j = first_valid_entry; j < first_invalid_entry; ++j) {
      const int h_id = hgraph->vind_[j];
      if (boundary_net_flag[h_id] == true) {
        boundary_set.insert(i);
        break;
      }
    }
  }
  std::vector<int> boundary_vertices (boundary_set.begin(), boundary_set.end());
  return boundary_vertices;
}
*/
}  // namespace par
