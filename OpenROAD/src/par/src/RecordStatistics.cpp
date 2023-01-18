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

#include "RecordStatistics.h"

#include <set>

#include "TPHypergraph.h"

namespace par {

void Obfuscator::QuickMatching(
    const HGraph hgraph,
    std::vector<int>& vertex_c_attr,
    std::vector<std::vector<float>>& vertex_weights_c,
    std::vector<int>& community_attr_c,
    std::vector<int>& fixed_attr_c,
    std::vector<std::vector<float>>& placement_attr_c)
{
  const int num_vtx_threshold = hgraph->num_vertices_ * contraction_factor_;
  // reset variables
  int num_clusters_est
      = hgraph->num_vertices_;  // the estimated number of clusters
                                // (existing clusters + remaining vertices)
                                // This variable is used for early stop
  int cluster_id = 0;           // the number of clusters in coarser hypergraph
  vertex_c_attr.clear();
  vertex_c_attr.resize(hgraph->num_vertices_);
  std::fill(
      vertex_c_attr.begin(), vertex_c_attr.end(), -1);  // reset vertex_c_attr
  vertex_weights_c.clear();  // reset vertex_weights_c
  placement_attr_c.clear();  // reset placement_attr_c
  community_attr_c.clear();  // reset community_attr_c
  fixed_attr_c.clear();      // reset fixed_attr_c

  // Compute the algebraic weights based on placement location
  std::vector<float> algebraic_weight = hgraph->ComputeAlgebraicWights();
  // Checking free vertices in hgraph
  std::vector<int> unvisited;
  // Step1 : Check if there exists fixed vertices
  if (hgraph->fixed_vertex_flag_ == false) {
    unvisited.resize(hgraph->num_vertices_);
    std::iota(unvisited.begin(), unvisited.end(), 0);  // Fill with 0, 1, ...
  } else {
    for (int v = 0; v < hgraph->num_vertices_;
         v++) {  // traverse fixed vertices
      if (hgraph->fixed_attr_[v] > -1) {
        vertex_c_attr[v] = cluster_id++;  // identify the fixed vertices first
        vertex_weights_c.push_back(hgraph->vertex_weights_[v]);
        fixed_attr_c.push_back(hgraph->fixed_attr_[v]);
        if (hgraph->community_flag_ == true)
          community_attr_c.push_back(hgraph->community_attr_[v]);
        if (hgraph->placement_flag_ == true)
          placement_attr_c.push_back(hgraph->placement_attr_[v]);
      } else {
        unvisited.push_back(v);
      }
    }  // done fixed vertices
  }
  // Step 2: Randomly shuffle remaining vertices and calculate node matching
  shuffle(
      unvisited.begin(), unvisited.end(), std::default_random_engine(seed_));
  for (auto v_iter = unvisited.begin(); v_iter != unvisited.end(); v_iter++) {
    const int v = *v_iter;  // current vertex being visited
    if (vertex_c_attr[v] != -1)
      continue;  // current vertex has been visited
    std::map<int, float>
        score_map;  // store the score for each candidate vertex
    const int v_start_idx = hgraph->vptr_[v];
    const int v_end_idx
        = hgraph->vptr_[v + 1];  // get all the hyperedges related to v
    for (int v_idx = v_start_idx; v_idx < v_end_idx;
         v_idx++) {                        // traverse hyperedge
      const int e = hgraph->vind_[v_idx];  // hyperedge_id
      const int e_start_idx = hgraph->eptr_[e];
      const int e_end_idx = hgraph->eptr_[e + 1];
      const int e_size = e_end_idx - e_start_idx;
      // ignore single-vertex net and large net
      if (e_size <= 1 || e_size > contract_global_net_threshold_)
        continue;
      // calculate the score from hyperedge
      const float e_score
          = std::inner_product(hgraph->hyperedge_weights_[e].begin(),
                               hgraph->hyperedge_weights_[e].end(),
                               e_wt_factors_.begin(),
                               0.0)
            / (e_size - 1) * algebraic_weight[e];
      // traverse current hyperedge
      for (int e_idx = e_start_idx; e_idx < e_end_idx; e_idx++) {
        const int u = hgraph->eind_[e_idx];  // candiate vertex
        if (u == v)
          continue;                                  // skip v itself
        if (score_map.find(u) != score_map.end()) {  // u is a known candidate
          score_map[u] += e_score;
        } else {  // u is a possible candidate
          // u is a good candidate if it can statisfy
          // fixed vertex, community and max_vertex_weights_ threshold
          const std::vector<float>& u_weight
              = vertex_c_attr[u] > -1 ? vertex_weights_c[vertex_c_attr[u]]
                                      : hgraph->vertex_weights_[u];
          if ((hgraph->fixed_vertex_flag_ == true
               && hgraph->fixed_attr_[u] > -1)
              || (hgraph->community_flag_ == true
                  && hgraph->community_attr_[u] != hgraph->community_attr_[v])
              || (hgraph->vertex_weights_[v] + u_weight > max_vertex_weight_)) {
            continue;
          }
          score_map[u] = e_score;  // u is a candidate
        }                          // end u is a possible candidate
      }                            // end hyperedge e
    }                              // end traverse hyperedges connected to v
    // If there is no candidate ...
    if (score_map.size() == 0) {  // no candidate
      vertex_c_attr[v] = cluster_id++;
      vertex_weights_c.push_back(hgraph->vertex_weights_[v]);
      // update placement, community and fixed attr
      if (hgraph->placement_flag_ == true)
        placement_attr_c.push_back(hgraph->placement_attr_[v]);
      if (hgraph->community_flag_ == true)
        community_attr_c.push_back(hgraph->community_attr_[v]);
      if (hgraph->fixed_vertex_flag_ == true)
        fixed_attr_c.push_back(-1);  // free cluster
      continue;                      // done for current vertex
    }
    // add score contriuted from timing path
    if (hgraph->num_timing_paths_ >= 1 && path_traverse_step_ > 0) {
      const int v_start_idx = hgraph->pptr_v_[v];
      const int v_end_idx = hgraph->pptr_v_[v + 1];
      std::map<int, float> timing_neighbors;  // neighbors in timing paths
      // traverse timing path related to v
      for (int v_idx = v_start_idx; v_idx < v_end_idx; v_idx++) {
        const int p = hgraph->pind_v_[v_idx];  // current timing path
        const int p_start_idx = hgraph->vptr_p_[p];
        const int p_end_idx = hgraph->vptr_p_[p + 1];
        // divide into weights based on distance
        const float timing_attr
            = hgraph->timing_attr_[p] * timing_factor_ / path_traverse_step_;
        int p_idx = hgraph->vptr_p_[p];  // loop variable
        while (p_idx < p_end_idx) {
          if (hgraph->vind_p_[p_idx] == v) {  // traverse around the vertex
            // use lamda function to traverse left and right
            auto TraversePath = [&](bool direction) {
              int i = 0;
              while (i++ < path_traverse_step_) {
                const int idx = (direction == true) ? p_idx + i : p_idx - i;
                if (idx < p_start_idx || idx >= p_end_idx)
                  return;  // exceed current path
                const int u = hgraph->vind_p_[idx];
                if (u == v)
                  return;  // stop when encounter v again
                if (timing_neighbors.find(u) == timing_neighbors.end())
                  timing_neighbors[u] = i * timing_attr;
                else
                  timing_neighbors[u] += i * timing_attr;
              }
            };
            TraversePath(true);   // right
            TraversePath(false);  // left
          }                       // done current node
          p_idx++;
        }  // done current path
      }    // done all the timing paths related to v

      // add timing contributions to the overall score
      for (auto& [u, score] : timing_neighbors)
        if (score_map.find(u) != score_map.end())
          score_map[u] += score;
    }  // done timing paths

    // add the score from placement and find the best score.
    // find the best candidate
    float best_score
        = -std::numeric_limits<float>::max();  // score will be always
                                               // nonnegative
    int best_u = -1;
    for (auto& [u, score] : score_map) {
      if (hgraph->placement_flag_ == true) {  // add score from placement
        const std::vector<float>& u_placement
            = vertex_c_attr[u] > -1 ? placement_attr_c[vertex_c_attr[u]]
                                    : hgraph->placement_attr_[u];
        score += norm2(u_placement - hgraph->placement_attr_[v], p_wt_factors_);
      }
      // find better best_u
      // break ties with unclustered vertex
      if ((score > best_score)
          || (score == best_score && vertex_c_attr[u] == -1)) {
        best_u = u;
        best_score = score;
      }
    }

    // merge v and best_u
    const int best_u_c_attr = vertex_c_attr[best_u];
    if (best_u_c_attr > -1) {            // best_u has been clustered
      vertex_c_attr[v] = best_u_c_attr;  // update cluster id for v
      // update placement_attr first
      if (hgraph->placement_flag_ == true) {
        const float w1 = norm2(vertex_weights_c[best_u_c_attr], v_wt_factors_);
        const float w2 = norm2(hgraph->vertex_weights_[v], v_wt_factors_);
        placement_attr_c[best_u_c_attr]
            = placement_attr_c[best_u_c_attr] * (w1 / (w1 + w2))
              + hgraph->placement_attr_[v] * (w2 / (w1 + w2));
      }
      vertex_weights_c[best_u_c_attr]
          = vertex_weights_c[best_u_c_attr]
            + hgraph->vertex_weights_[v];  // update weight
    } else {                               // best_u is unclustered
      vertex_c_attr[v] = cluster_id;
      vertex_c_attr[best_u] = cluster_id++;

      // update placement_attr first
      if (hgraph->placement_flag_ == true) {
        const float w1 = norm2(hgraph->vertex_weights_[best_u], v_wt_factors_);
        const float w2 = norm2(hgraph->vertex_weights_[v], v_wt_factors_);

        placement_attr_c.push_back(
            hgraph->placement_attr_[best_u] * (w1 / (w1 + w2))
            + hgraph->placement_attr_[v] * (w2 / (w1 + w2)));
      }
      // then update weight
      vertex_weights_c.push_back(hgraph->vertex_weights_[best_u]
                                 + hgraph->vertex_weights_[v]);
      // update community_attr and fixed_attr
      if (hgraph->community_flag_ == true)
        community_attr_c.push_back(hgraph->community_attr_[best_u]);
      if (hgraph->fixed_vertex_flag_ == true)
        fixed_attr_c.push_back(-1);  // free cluster
    }

    // check early stop condition
    num_clusters_est--;  // reduce one vertex in coarser hypergraph

    if (num_clusters_est < num_vtx_threshold) {
      while (v_iter != unvisited.end()) {
        const int& v = *v_iter++;
        if (vertex_c_attr[v] == -1) {  // if v is unclustered
          vertex_c_attr[v] = cluster_id++;
          vertex_weights_c.push_back(hgraph->vertex_weights_[v]);
          // update attr
          if (hgraph->fixed_vertex_flag_ == true)
            fixed_attr_c.push_back(-1);  // free cluster
          if (hgraph->placement_flag_ == true)
            placement_attr_c.push_back(hgraph->placement_attr_[v]);
          if (hgraph->community_flag_ == true)
            community_attr_c.push_back(hgraph->community_attr_[v]);
        }  // done for v
      }
      return;
    }
  }  // done for traversing vertices
}

HGraph Obfuscator::OneLevelContraction(HGraph hypergraph)
{
  std::vector<int> vertex_c_attr;
  std::vector<std::vector<float>> vertex_weights_c;
  std::vector<int> community_attr_c;
  std::vector<int> fixed_attr_c;
  std::vector<std::vector<float>> placement_attr_c;
  QuickMatching(hypergraph,
                vertex_c_attr,
                vertex_weights_c,
                community_attr_c,
                fixed_attr_c,
                placement_attr_c);
  // Step 2:
  // Update the hyperedges
  std::vector<std::vector<int>> hyperedges_c;  // hyperedge in hgraph_c
  std::vector<std::vector<float>> hyperedge_weights_c;
  // variables used to detect parallel nets and single nets
  std::map<long long int, int>
      hash_map;  // store hash value for hyperedge in hgraph_c
  // traverse all the hyperedges in hgraph
  for (int e = 0; e < hypergraph->num_hyperedges_; e++) {
    const int e_start_idx = hypergraph->eptr_[e];
    const int e_end_idx = hypergraph->eptr_[e + 1];
    const int e_size = e_end_idx - e_start_idx;
    // ignore single-vertex net and large net
    if (e_size <= 1 || e_size > global_net_threshold_)
      continue;
    std::set<int> hyperedge_c;
    for (int idx = e_start_idx; idx < e_end_idx; idx++)
      hyperedge_c.insert(vertex_c_attr[hypergraph->eind_[idx]]);
    if (hyperedge_c.size() <= 1)
      continue;  // ignore singl-vertex hyperedge
    const long long int hash_value
        = std::inner_product(hyperedge_c.begin(),
                             hyperedge_c.end(),
                             hyperedge_c.begin(),
                             static_cast<long long int>(0));
    if (hash_map.find(hash_value) == hash_map.end()) {
      hash_map[hash_value] = static_cast<int>(hyperedges_c.size());
      hyperedge_weights_c.push_back(hypergraph->hyperedge_weights_[e]);
      hyperedges_c.push_back(
          std::vector<int>(hyperedge_c.begin(), hyperedge_c.end()));
    } else {
      // check if the hyperedge is parallel to the first hyperedge with the same
      // hash_value there may be many hyperedges with the same hash_value but we
      // will only check the first one for simplicity
      const int hash_id = hash_map[hash_value];
      std::vector<int> hyperedge_vec(hyperedge_c.begin(), hyperedge_c.end());
      if (hyperedges_c[hash_id] == hyperedge_vec) {
        hyperedge_weights_c[hash_id]
            = hyperedge_weights_c[hash_id] + hypergraph->hyperedge_weights_[e];
      } else {
        hyperedge_weights_c.push_back(hypergraph->hyperedge_weights_[e]);
        hyperedges_c.push_back(hyperedge_vec);
      }
    }  // done for current hyperedge
  }    // done for hyperedge traversal
  // Step 3:
  // update the timing path
  std::vector<std::vector<int>> paths_c;  // timing paths in hgraph_c
  std::vector<float> timing_attr_c;       // weight for timing paths in hgraph_c
  hash_map.clear();  // use hash_map to detect the parallel timing paths
  if (hypergraph->num_timing_paths_ >= 1) {
    for (int p = 0; p < hypergraph->num_timing_paths_; p++) {
      const int p_start_idx = hypergraph->vptr_p_[p];
      const int p_end_idx = hypergraph->vptr_p_[p + 1];
      if (p_end_idx - p_start_idx <= 1)
        continue;
      // update path
      std::vector<int> path{hypergraph->vind_p_[p_start_idx]};
      for (int idx = p_start_idx + 1; idx < p_end_idx; idx++)
        if (path.back() != hypergraph->vind_p_[idx])
          path.push_back(hypergraph->vind_p_[idx]);
      // ignore single-vertex path
      if (path.size() <= 1)
        continue;
      // check possible parallel condition
      const long long int hash_value
          = std::inner_product(path.begin(),
                               path.end(),
                               path.begin(),
                               static_cast<long long int>(0));
      if (hash_map.find(hash_value) == hash_map.end()) {
        hash_map[hash_value] = static_cast<int>(paths_c.size());
        paths_c.push_back(path);
        timing_attr_c.push_back(hypergraph->timing_attr_[p]);
      } else {
        // check if the path is parallel to the first path with the same
        // hash_value there may be many paths with the same hash_value but we
        // will only check the first one for simplicity
        const int hash_id = hash_map[hash_value];
        if (paths_c[hash_id] == path) {
          timing_attr_c[hash_id]
              = timing_attr_c[hash_id] + hypergraph->timing_attr_[p];
        } else {
          paths_c.push_back(path);
          timing_attr_c.push_back(hypergraph->timing_attr_[p]);
        }
      }  // done for current path
    }    // done for all the timing path
  }      // done for timing cases

  // create the coarser hgraph
  HGraph hgraph_c
      = std::make_shared<TPHypergraph>(hypergraph->vertex_dimensions_,
                                       hypergraph->hyperedge_dimensions_,
                                       hyperedges_c,
                                       vertex_weights_c,
                                       hyperedge_weights_c,
                                       fixed_attr_c,
                                       community_attr_c,
                                       hypergraph->placement_dimensions_,
                                       placement_attr_c,
                                       paths_c,
                                       timing_attr_c,
                                       logger_);

  hypergraph->vertex_c_attr_ = vertex_c_attr;  // update clustering attr
  return hgraph_c;
}

void Obfuscator::FindAverageHyperedgeSize(const HGraph hypergraph)
{
  int total_size = 0;
  for (int i = 0; i < hypergraph->num_hyperedges_; ++i) {
    const int first_valid_entry = hypergraph->eptr_[i];
    const int first_invalid_entry = hypergraph->eptr_[i + 1];
    total_size += (first_invalid_entry - first_valid_entry);
  }
  avg_hyperedge_size_ = static_cast<float>(total_size) / static_cast<float>(hypergraph->num_hyperedges_);
}

void Obfuscator::GenerateFanInsOuts(const HGraph hypergraph)
{
  for (int i = 0; i < hypergraph->num_hyperedges_; ++i) {
    const int first_valid_entry = hypergraph->eptr_[i];
    const int first_invalid_entry = hypergraph->eptr_[i + 1];
    const int source = hypergraph->eind_[first_valid_entry];
    ++fanouts_[source];
    for (int j = first_valid_entry + 1; j < first_invalid_entry; ++j) {
      const int sink = hypergraph->eind_[j];
      ++fanins_[sink];
    }
  }
}

void Obfuscator::GenerateTerminals(const HGraph hypergraph)
{
  for (int i = 0; i < hypergraph->num_vertices_; ++i) {
    const int first_valid_entry = hypergraph->vptr_[i];
    const int first_invalid_entry = hypergraph->vptr_[i+1];
    terminals_[i] = (first_invalid_entry - first_valid_entry);
  }
}

void Obfuscator::Obfuscate(HGraph hypergraph,
                           const float contraction_factor,
                           const int global_net_threshold,
                           const int contract_global_net_threshold,
                           const int seed)
{
  // Do 1 level of coarsening here
  SetContractionFactor(contraction_factor);
  SetGlobalNetThreshold(global_net_threshold);
  SetContractGlobalNetThreshold(contract_global_net_threshold);
  SetRandSeed(seed);
  auto stats = RecordStatistics(hypergraph);
  std::cout << "[OBFUSCATE] Average fan in: " << stats->GetAvgFanIn() << std::endl;
  std::cout << "[OBFUSCATE] Average fan out: " << stats->GetAvgFanOut() << std::endl;
  std::cout << "[OBFUSCATE] Average terminals per vertex: " << stats->GetAvgTerminals() << std::endl;
  std::cout << "[OBFUSCATE] Average pins per hyperedge: " << stats->GetAvgHyperedgeSize() << std::endl;
  std::cout << "[OBFUSCATE] Finished setting coarsening parameters "
            << std::endl;
  auto hg = OneLevelContraction(hypergraph);
  std::cout << "[OBFUSCATE] Clustered netlist has " << hg->num_vertices_
            << " vertices and " << hg->num_hyperedges_ << " hyperedges" << std::endl;
  auto stats_c = RecordStatistics(hg);
  std::cout << "[OBFUSCATE] Average fan in: " << stats_c->GetAvgFanIn() << std::endl;
  std::cout << "[OBFUSCATE] Average fan out: " << stats_c->GetAvgFanOut() << std::endl;
  std::cout << "[OBFUSCATE] Average terminals per vertex: " << stats_c->GetAvgTerminals() << std::endl;
  std::cout << "[OBFUSCATE] Average pins per hyperedge: " << stats_c->GetAvgHyperedgeSize() << std::endl;
}

}  // namespace par
