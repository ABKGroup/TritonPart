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
// High-level description
// This is the interfaces for TritonPart
// It works in two ways:
// 1) default mode :  read a verilog netlist and extract the hypergraph based on
// netlist 2) classical mode : read the hypergraph file in hmetis format
///////////////////////////////////////////////////////////////////////////////
#include "TritonPart.h"

#include <string>

#include "Coarsening.h"
#include "ILPbasedRefinement.h"
#include "KPMRefinement.h"
#include "Multilevel.h"
#include "Partitioner.h"
#include "RecordStatistics.h"
#include "TPCoarsener.h"
#include "TPHypergraph.h"
#include "TPMultilevel.h"
#include "TPPartitioner.h"
#include "TPRefiner.h"
#include "Utilities.h"
#include "odb/db.h"
#include "sta/ArcDelayCalc.hh"
#include "sta/Bfs.hh"
#include "sta/Corner.hh"
#include "sta/DcalcAnalysisPt.hh"
#include "sta/ExceptionPath.hh"
#include "sta/FuncExpr.hh"
#include "sta/Graph.hh"
#include "sta/GraphDelayCalc.hh"
#include "sta/Liberty.hh"
#include "sta/Network.hh"
#include "sta/PathAnalysisPt.hh"
#include "sta/PathEnd.hh"
#include "sta/PathExpanded.hh"
#include "sta/PathRef.hh"
#include "sta/PatternMatch.hh"
#include "sta/PortDirection.hh"
#include "sta/Sdc.hh"
#include "sta/Search.hh"
#include "sta/SearchPred.hh"
#include "sta/Sequential.hh"
#include "sta/Sta.hh"
#include "sta/Units.hh"
#include "utl/Logger.h"
// julia interfaces
//#include "julia_init.h"
//#include "spectral.h"

// PAR 2502

using utl::PAR;

namespace par {

// Read hypergraph from input files
void TritonPart::ReadHypergraph(std::string hypergraph_file,
                                std::string fixed_file)
{
  std::ifstream hypergraph_file_input(hypergraph_file);
  if (!hypergraph_file_input.is_open()) {
    logger_->error(PAR,
                   2500,
                   "Can not open the input hypergraph file : {}",
                   hypergraph_file);
  }

  // Set the flag variables
  // vertex_dimensions_ = 1;
  // hyperedge_dimensions_ = 1;
  // placement_dimensions_ = 0;
  timing_aware_flag_ = false;

  // Check the number of vertices, number of hyperedges, weight flag
  std::string cur_line;
  std::getline(hypergraph_file_input, cur_line);
  std::istringstream cur_line_buf(cur_line);
  std::vector<int> stats{std::istream_iterator<int>(cur_line_buf),
                         std::istream_iterator<int>()};
  num_hyperedges_ = stats[0];
  num_vertices_ = stats[1];
  bool hyperedge_weight_flag = false;
  bool vertex_weight_flag = false;
  if (stats.size() == 3) {
    if ((stats[2] % 10) == 1)
      hyperedge_weight_flag = true;

    if (stats[2] >= 10)
      vertex_weight_flag = true;
  }

  // Read hyperedge information
  for (int i = 0; i < num_hyperedges_; i++) {
    std::getline(hypergraph_file_input, cur_line);
    if (hyperedge_weight_flag == true) {
      std::istringstream cur_line_buf(cur_line);
      std::vector<float> hvec{std::istream_iterator<float>(cur_line_buf),
                              std::istream_iterator<float>()};
      std::vector<float>::iterator breakpoint{hvec.begin()
                                              + hyperedge_dimensions_};
      // read first hyperedge_dimensions_ elements as hyperege weights
      // std::vector<float> hwts(hvec.begin(), std::advance(hvec.begin(),
      // hyperedge_dimensions_));
      std::vector<float> hwts(hvec.begin(), breakpoint);
      // read remaining elements as hyperedge
      // std::vector<int> hyperedge(std::advance(hvec.begin(),
      // hyperedge_dimensions_), hvec.end());
      std::vector<int> hyperedge(breakpoint, hvec.end());
      for (auto& value : hyperedge)
        value--;
      if (hyperedge.size() > global_net_threshold_)
        continue;

      hyperedge_weights_.push_back(hwts);
      hyperedges_.push_back(hyperedge);
    } else {
      std::istringstream cur_line_buf(cur_line);
      std::vector<int> hyperedge{std::istream_iterator<int>(cur_line_buf),
                                 std::istream_iterator<int>()};
      for (auto& value : hyperedge)
        value--;
      std::vector<float> hwts(hyperedge_dimensions_, 1.0);
      if (hyperedge.size() > global_net_threshold_)
        continue;
      hyperedge_weights_.push_back(hwts);
      hyperedges_.push_back(hyperedge);
    }
  }

  // Read weight for vertices
  for (int i = 0; i < num_vertices_; i++) {
    if (vertex_weight_flag == true) {
      std::istringstream cur_line_buf(cur_line);
      std::vector<float> vwts{std::istream_iterator<float>(cur_line_buf),
                              std::istream_iterator<float>()};
      vertex_weights_.push_back(vwts);
    } else {
      std::vector<float> vwts(vertex_dimensions_, 1.0);
      vertex_weights_.push_back(vwts);
    }
  }

  // Read fixed vertices
  if (fixed_file.size() > 0) {
    int part_id = -1;
    fixed_vertex_flag_ = true;
    std::ifstream fixed_file_input(fixed_file);
    if (!fixed_file_input.is_open()) {
      logger_->error(PAR, 2501, "Can not open the fixed file : {}", fixed_file);
    }
    for (int i = 0; i < num_vertices_; i++) {
      fixed_file_input >> part_id;
      fixed_attr_.push_back(part_id);
    }
    fixed_file_input.close();
  }
  num_vertices_ = vertex_weights_.size();
  num_hyperedges_ = hyperedge_weights_.size();
}

// Convert the netlist into hypergraphs
void TritonPart::ReadNetlist()
{
  // assign vertex_id property of each instance and each IO port
  int vertex_id = 0;
  for (auto term : block_->getBTerms()) {
    odb::dbIntProperty::create(term, "vertex_id", vertex_id++);
  }
  for (auto inst : block_->getInsts()) {
    odb::dbIntProperty::create(inst, "vertex_id", vertex_id++);
  }
  num_vertices_ = vertex_id;

  // Each net correponds to an hyperedge
  // Traverse the hyperedge and assign hyperedge id
  int hyperedge_id = 0;
  for (auto net : block_->getNets()) {
    odb::dbIntProperty::create(net, "hyperedge_id", -1);
    // ignore all the power net
    if (net->getSigType().isSupply())
      continue;
    // check the hyperedge
    int driver_id = -1;      // vertex id of the driver instance
    std::set<int> loads_id;  // vertex id of sink instances
    // check the connected instances
    for (odb::dbITerm* iterm : net->getITerms()) {
      odb::dbInst* inst = iterm->getInst();
      const int vertex_id
          = odb::dbIntProperty::find(inst, "vertex_id")->getValue();
      if (iterm->getIoType() == odb::dbIoType::OUTPUT)
        driver_id = vertex_id;
      else
        loads_id.insert(vertex_id);
    }
    // check the connected IO pins
    for (odb::dbBTerm* bterm : net->getBTerms()) {
      const int vertex_id
          = odb::dbIntProperty::find(bterm, "vertex_id")->getValue();
      if (bterm->getIoType() == odb::dbIoType::INPUT)
        driver_id = vertex_id;
      else
        loads_id.insert(vertex_id);
    }
    // check the hyperedges
    std::vector<int> hyperedge;
    if (driver_id != -1 && loads_id.size() > 0) {
      hyperedge.push_back(driver_id);
      for (auto& load_id : loads_id)
        if (load_id != driver_id)
          hyperedge.push_back(load_id);
    }
    // Ignore all the single-vertex hyperedge
    if (hyperedge.size() > 1 && hyperedge.size() <= global_net_threshold_) {
      hyperedges_.push_back(hyperedge);
      odb::dbIntProperty::find(net, "hyperedge_id")->setValue(hyperedge_id++);
    }
  }  // finish hyperedge

  num_hyperedges_ = static_cast<int>(hyperedges_.size());

  // initialize other parameters
  // Set the flag variables
  vertex_dimensions_ = 1;
  hyperedge_dimensions_ = 1;
  placement_dimensions_ = 0;
  timing_aware_flag_ = false;
  for (int i = 0; i < num_vertices_; i++) {
    std::vector<float> vwts(vertex_dimensions_, 1.0);
    vertex_weights_.push_back(vwts);
    std::vector<float> hwts(hyperedge_dimensions_, 1.0);
  }

  // add timing feature
  if (timing_aware_flag_ == true)
    BuildTimingPaths();  // create timing paths
}

// Find all the critical timing paths
// The codes below similar to gui/src/staGui.cpp
// Please refer to sta/Search/ReportPath.cc for how to check the timing path
void TritonPart::BuildTimingPaths()
{
  if (timing_aware_flag_ == false || top_n_ <= 0)
    return;

  sta_->ensureGraph();
  sta_->searchPreamble();

  sta::ExceptionFrom* e_from = nullptr;
  sta::ExceptionThruSeq* e_thrus = nullptr;
  sta::ExceptionTo* e_to = nullptr;
  bool include_unconstrained = false;
  bool get_max = true;  // max for setup check, min for hold check
  // Timing paths are grouped into path groups according to the clock
  // associated with the endpoint of the path, for example, path group for clk
  int group_count = top_n_;
  int endpoint_count = 1;  // The number of paths to report for each endpoint.
  bool unique_pins
      = true;  // Only the worst path through the set of pins is reported
  // Definition for findPathEnds function in Search.hh
  // PathEndSeq *findPathEnds(ExceptionFrom *from,
  //              ExceptionThruSeq *thrus,
  //              ExceptionTo *to,
  //              bool unconstrained,
  //              const Corner *corner,
  //              const MinMaxAll *min_max,
  //              int group_count,
  //              int endpoint_count,
  //              bool unique_pins,
  //              float slack_min,
  //              float slack_max,
  //              bool sort_by_slack,
  //              PathGroupNameSet *group_names,
  //              bool setup,
  //              bool hold,
  //              bool recovery,
  //              bool removal,
  //              bool clk_gating_setup,
  //              bool clk_gating_hold);
  sta::PathEndSeq* path_ends
      = sta_->search()->findPathEnds(  // from, thrus, to, unconstrained
          e_from,
          e_thrus,
          e_to,
          include_unconstrained,
          // corner, min_max,
          sta_->cmdCorner(),
          get_max ? sta::MinMaxAll::max() : sta::MinMaxAll::min(),
          // group_count, endpoint_count, unique_pins
          group_count,     // path_count used by GUI, not sure what happened
          endpoint_count,  // path_count used by GUI, not sure what happened
          true,
          -sta::INF,
          sta::INF,  // slack_min, slack_max,
          true,      // sort_by_slack
          nullptr,   // group_names
          // setup, hold, recovery, removal,
          get_max,
          !get_max,
          false,
          false,
          // clk_gating_setup, clk_gating_hold
          false,
          false);

  // check all the timing paths
  for (auto& path_end : *path_ends) {
    sta_->reportPathEnd(path_end);
    auto* path = path_end->path();
    TimingPath timing_path;               // create the timing path
    float slack = path_end->slack(sta_);  // slack information
    // slack = slack / sta_->search()->units()->timeUnit()->scale();
    timing_path.slack = slack;
    sta::PathExpanded expand(path, sta_);
    for (size_t i = 0; i < expand.size(); i++) {
      sta::PathRef* ref = expand.path(i);
      sta::Pin* pin = ref->vertex(sta_)->pin();
      auto net = network_->net(pin);  // sta::Net*
      if (net == nullptr)
        continue;  // check if the net exists
      if (network_->isTopLevelPort(pin) == true) {
        auto bterm = block_->findBTerm(network_->pathName(pin));
        int vertex_id
            = odb::dbIntProperty::find(bterm, "vertex_id")->getValue();
        if (std::find(
                timing_path.path.begin(), timing_path.path.end(), vertex_id)
            == timing_path.path.end())
          timing_path.path.push_back(vertex_id);
      } else {
        auto inst = network_->instance(pin);
        auto dbinst = block_->findInst(network_->pathName(inst));
        int vertex_id
            = odb::dbIntProperty::find(dbinst, "vertex_id")->getValue();
        if (std::find(
                timing_path.path.begin(), timing_path.path.end(), vertex_id)
            == timing_path.path.end())
          timing_path.path.push_back(vertex_id);
      }

      auto dbnet = block_->findNet(
          network_->pathName(net));  // convert sta::Net* to dbNet*
      int hyperedge_id
          = odb::dbIntProperty::find(dbnet, "hyperedge_id")->getValue();
      if (hyperedge_id == -1) {
        continue;
      }
      if (std::find(
              timing_path.arcs.begin(), timing_path.arcs.end(), hyperedge_id)
          == timing_path.arcs.end())
        timing_path.arcs.push_back(hyperedge_id);
    }
    // add timing path
    timing_paths_.push_back(timing_path);
  }

  // release memory
  delete path_ends;
}

// Create the hypergraph object hypergraph_
void TritonPart::BuildHypergraph()
{
  // add hyperedge
  std::vector<int> eind;
  std::vector<int> eptr;  // hyperedges
  eptr.push_back(static_cast<int>(eind.size()));
  for (auto hyperedge : hyperedges_) {
    eind.insert(eind.end(), hyperedge.begin(), hyperedge.end());
    eptr.push_back(static_cast<int>(eind.size()));
  }
  // add vertex
  // create vertices from hyperedges
  std::vector<std::vector<int>> vertices(num_vertices_);
  for (int i = 0; i < num_hyperedges_; i++)
    for (auto v : hyperedges_[i])
      vertices[v].push_back(i);  // i is the hyperedge id
  std::vector<int> vind;
  std::vector<int> vptr;  // vertices
  vptr.push_back(static_cast<int>(vind.size()));
  for (auto& vertex : vertices) {
    vind.insert(vind.end(), vertex.begin(), vertex.end());
    vptr.push_back(static_cast<int>(vind.size()));
  }

  // Convert the timing information
  std::vector<int> vind_p;  // each timing path is a sequences of vertices
  std::vector<int> vptr_p;
  std::vector<int>
      pind_v;  // store all the timing paths connected to the vertex
  std::vector<int> pptr_v;
  std::vector<float> timing_attr;

  // create TPHypergraph
  hypergraph_ = std::make_shared<TPHypergraph>(num_vertices_,
                                               num_hyperedges_,
                                               vertex_dimensions_,
                                               hyperedge_dimensions_,
                                               eind,
                                               eptr,
                                               vind,
                                               vptr,
                                               vertex_weights_,
                                               hyperedge_weights_,
                                               fixed_attr_,
                                               community_attr_,
                                               placement_dimensions_,
                                               placement_attr_,
                                               vind_p,
                                               vptr_p,
                                               pind_v,
                                               pptr_v,
                                               timing_attr,
                                               logger_);
}

// Partition the design
// The first step is to convert the netlist into a hypergraph
// The second step is to get all the features such as timing paths
void TritonPart::tritonPartDesign(unsigned int num_parts_arg,
                                  float balance_constraint_arg,
                                  unsigned int seed_arg)
{
  logger_->info(PAR, 2300, "Staring TritonPart (TritonPartDesign)");
  auto start_timestamp_global = std::chrono::high_resolution_clock::now();

  // Parameters
  num_parts_ = num_parts_arg;
  ub_factor_ = balance_constraint_arg;
  seed_ = seed_arg;
  srand(seed_);  // set the random seed
  logger_->info(PAR,
                2301,
                "num_parts = {}\n"
                "UBfactor = {}\n"
                "seed = {}\n",
                num_parts_,
                ub_factor_,
                seed_);
  // only for TritonPartDesign, we need the block_ information
  block_ = db_->getChip()->getBlock();

  // build hypergraph
  // for IO port and insts (std cells and macros),
  // there is an attribute for vertex_id
  ReadNetlist();
  BuildHypergraph();
  logger_->info(PAR,
                2302,
                "num_vertices = {}\n"
                "num_hyperedges = {}\n",
                num_vertices_,
                num_hyperedges_);

  auto end_timestamp_global = std::chrono::high_resolution_clock::now();
  double total_global_time
      = std::chrono::duration_cast<std::chrono::nanoseconds>(
            end_timestamp_global - start_timestamp_global)
            .count();
  total_global_time *= 1e-9;
  logger_->info(PAR, 2399, "Total Runtime: {:0.2f} sec", total_global_time);
  logger_->info(PAR, 2400, "Exiting TritonPart (TritonPartDesign)");
}

HGraph TritonPart::preProcessHypergraph()
{
  std::cout << "Pre-processing hypergraph by temporarily removing hyperedges "
               "of size**"
            << he_size_threshold_ << std::endl;
  std::vector<std::vector<int>> hyperedges_p;
  std::vector<std::vector<float>> hyperedge_weights_p;
  int num_hyperedges_p = 0;
  for (int i = 0; i < num_hyperedges_; ++i) {
    const int he_size = hyperedges_[i].size();
    if (he_size <= he_size_threshold_) {
      std::vector<int> he = hyperedges_[i];
      hyperedges_p.push_back(he);
      std::vector<float> hwt = hyperedge_weights_[i];
      hyperedge_weights_p.push_back(hwt);
      ++num_hyperedges_p;
    }
  }
  // add hyperedge
  std::vector<int> eind_p;
  std::vector<int> eptr_p;  // hyperedges
  eptr_p.push_back(static_cast<int>(eind_p.size()));
  for (auto hyperedge : hyperedges_p) {
    eind_p.insert(eind_p.end(), hyperedge.begin(), hyperedge.end());
    eptr_p.push_back(static_cast<int>(eind_p.size()));
  }
  // add vertex
  // create vertices from hyperedges
  std::vector<std::vector<int>> vertices_p(num_vertices_);
  for (int i = 0; i < num_hyperedges_p; i++)
    for (auto v : hyperedges_p[i])
      vertices_p[v].push_back(i);  // i is the hyperedge id
  std::vector<int> vind_p;
  std::vector<int> vptr_p;  // vertices
  vptr_p.push_back(static_cast<int>(vind_p.size()));
  for (auto& vertex : vertices_p) {
    vind_p.insert(vind_p.end(), vertex.begin(), vertex.end());
    vptr_p.push_back(static_cast<int>(vind_p.size()));
  }
  // Convert the timing information
  std::vector<int> vind_p_p;  // each timing path is a sequences of vertices
  std::vector<int> vptr_p_p;
  std::vector<int>
      pind_v_p;  // store all the timing paths connected to the vertex
  std::vector<int> pptr_v_p;
  std::vector<float> timing_attr_p;

  // create TPHypergraph
  HGraph hypergraph_p = std::make_shared<TPHypergraph>(num_vertices_,
                                                       num_hyperedges_p,
                                                       vertex_dimensions_,
                                                       hyperedge_dimensions_,
                                                       eind_p,
                                                       eptr_p,
                                                       vind_p,
                                                       vptr_p,
                                                       vertex_weights_,
                                                       hyperedge_weights_p,
                                                       fixed_attr_,
                                                       community_attr_,
                                                       placement_dimensions_,
                                                       placement_attr_,
                                                       vind_p_p,
                                                       vptr_p_p,
                                                       pind_v_p,
                                                       pptr_v_p,
                                                       timing_attr_p,
                                                       logger_);
  return hypergraph_p;
}

void TritonPart::tritonPartHypergraph(const char* hypergraph_file_arg,
                                      const char* fixed_file_arg,
                                      unsigned int num_parts_arg,
                                      float balance_constraint_arg,
                                      int vertex_dimension_arg,
                                      int hyperedge_dimension_arg,
                                      unsigned int seed_arg)
{
  TritonPart_PartTwoWay(hypergraph_file_arg,
                        fixed_file_arg,
                        num_parts_arg,
                        balance_constraint_arg,
                        vertex_dimension_arg,
                        hyperedge_dimension_arg,
                        seed_arg);
  logger_->report("Starting TritonPart Partitioner");
  auto start_time_stamp_global = std::chrono::high_resolution_clock::now();
  // Parameters
  num_parts_ = num_parts_arg;
  ub_factor_ = balance_constraint_arg;
  seed_ = seed_arg;
  vertex_dimensions_ = vertex_dimension_arg;
  hyperedge_dimensions_ = hyperedge_dimension_arg;
  // local parameters
  std::string hypergraph_file = hypergraph_file_arg;
  std::string fixed_file = fixed_file_arg;
  logger_->report("Partition Parameters**");
  logger_->report("Number of partitions = {}", num_parts_);
  logger_->report("UBfactor = {}", ub_factor_);
  logger_->report("Vertex dimensions = {}", vertex_dimensions_);
  logger_->report("Hyperedge dimensions = {}", hyperedge_dimensions_);
  // build hypergraph
  ReadHypergraph(hypergraph_file, fixed_file);
  BuildHypergraph();
  logger_->report("Hypergraph Information**");
  logger_->report("#Vertices = {}", num_vertices_);
  logger_->report("#Hyperedges = {}", num_hyperedges_);
  // create coarsening class
  std::vector<float> e_wt_factors(hyperedge_dimensions_, 1.0);
  std::vector<float> v_wt_factors(vertex_dimensions_, 1.0);
  std::vector<float> p_wt_factors(placement_dimensions_ + 100, 1.0);
  float timing_factor = 1.0;
  int path_traverse_step = 2;
  std::vector<float> tot_vertex_weights = hypergraph_->GetTotalVertexWeights();
  int alpha = 4;
  std::vector<float> max_vertex_weights
      = DivideFactor(hypergraph_->GetTotalVertexWeights(), alpha * num_parts_);
  int smallest_v_size_cgraph = 250;
  int smallest_e_size_cgraph = 50;
  float coarsening_ratio = 1.5;
  int max_coarsen_iters = 20;

  // Netlist obfuscation related
  /*ObfuscatorPtr obfuscation = std::make_shared<Obfuscator>(e_wt_factors,
                                                           v_wt_factors,
                                                           p_wt_factors,
                                                           timing_factor,
                                                           path_traverse_step,
                                                           max_vertex_weights,
                                                           seed_,
                                                           logger_);

  const float contraction_factor = 0.85;
  const int global_net_threshold = 5000;
  const int contract_global_net_threshold = 500;
  const int seed = 0;
  obfuscation->Obfuscate(hypergraph_,
                         contraction_factor,
                         global_net_threshold,
                         contract_global_net_threshold,
                         seed);

  exit(EXIT_SUCCESS);*/

  CoarseningPtr coarsening
      = std::make_shared<Coarsening>(e_wt_factors,
                                     v_wt_factors,
                                     p_wt_factors,
                                     timing_factor,
                                     path_traverse_step,
                                     max_vertex_weights,
                                     global_net_threshold_,
                                     smallest_v_size_cgraph,
                                     smallest_e_size_cgraph,
                                     coarsening_ratio,
                                     max_coarsen_iters,
                                     seed_,
                                     logger_);
  // create partitioner class
  std::vector<std::vector<float>> vertex_balance
      = hypergraph_->GetVertexBalance(num_parts_, ub_factor_);
  std::vector<std::vector<float>> hyperedge_balance;
  float path_wt_factor = 1.0;
  float snaking_wt_factor = 1.0;
  float early_stop_ratio = 0.5;
  int max_num_fm_pass = 10;
  PartitionersPtr partitioners
      = std::make_shared<Partitioners>(num_parts_,
                                       e_wt_factors,
                                       path_wt_factor,
                                       snaking_wt_factor,
                                       early_stop_ratio,
                                       max_num_fm_pass,
                                       seed_,
                                       logger_);
  // create the refiner class
  KPMrefinementPtr kpmrefiner
      = std::make_shared<KPMRefinement>(num_parts_,
                                        e_wt_factors,
                                        path_wt_factor,
                                        snaking_wt_factor,
                                        early_stop_ratio,
                                        max_num_fm_pass,
                                        seed_,
                                        logger_);

  // create the ilp refiner class
  int wavefront
      = 50;  // wavefront is the number of vertices the ILP will consider
  IlpRefinerPtr ilprefiner = std::make_shared<IlpRefiner>(num_parts_,
                                                          seed_,
                                                          wavefront,
                                                          e_wt_factors,
                                                          path_wt_factor,
                                                          snaking_wt_factor,
                                                          max_num_fm_pass);

  // create the multilevel class
  bool v_cycle_flag = true;
  RefineType refine_type = KPMREFINEMENT;
  int num_initial_solutions = 20;      // number of initial random solutions
  int num_best_initial_solutions = 3;  // number of best initial solutions
  int num_ubfactor_delta = 5;  // allowing marginal imbalance to improve QoR
  int max_num_vcycle = 5;      // maximum number of vcycles

  MultiLevelHierarchyPtr multilevel_hierarchy
      = std::make_shared<MultiLevelHierarchy>(coarsening,
                                              partitioners,
                                              kpmrefiner,
                                              ilprefiner,
                                              num_parts_,
                                              v_cycle_flag,
                                              num_initial_solutions,
                                              num_best_initial_solutions,
                                              num_ubfactor_delta,
                                              max_num_vcycle,
                                              seed_,
                                              ub_factor_,
                                              refine_type,
                                              logger_);

  HGraph hypergraph_processed = preProcessHypergraph();
  logger_->report("\nPost processing hypergraph information**");
  logger_->report("#Vertices = {}", hypergraph_processed->GetNumVertices());
  logger_->report("#Hyperedges = {}", hypergraph_processed->GetNumHyperedges());

  std::vector<int> solution
      = multilevel_hierarchy->CallFlow(hypergraph_processed, vertex_balance);
  // check the existing solution
  std::pair<float, std::vector<std::vector<float>>> cutsize_balance
      = partitioners->GoldenEvaluator(hypergraph_, solution, true);
  // write the solution
  std::string solution_file
      = hypergraph_file + std::string(".part.") + std::to_string(num_parts_);
  WriteSolution(solution_file.c_str(), solution);
  auto end_timestamp_global = std::chrono::high_resolution_clock::now();
  double total_global_time
      = std::chrono::duration_cast<std::chrono::nanoseconds>(
            end_timestamp_global - start_time_stamp_global)
            .count();
  total_global_time *= 1e-9;
  logger_->report("Total runtime {}", total_global_time);
  logger_->report("Exiting TritonPart");
}

void TritonPart::TritonPart_PartTwoWay(const char* hypergraph_file_arg,
                                       const char* fixed_file_arg,
                                       unsigned int num_parts_arg,
                                       float balance_constraint_arg,
                                       int vertex_dimension_arg,
                                       int hyperedge_dimension_arg,
                                       unsigned int seed_arg)
{
  logger_->report("Starting TritonPart Partitioner");
  auto start_time_stamp_global = std::chrono::high_resolution_clock::now();
  // Parameters
  num_parts_ = num_parts_arg;
  ub_factor_ = balance_constraint_arg;
  seed_ = seed_arg;
  vertex_dimensions_ = vertex_dimension_arg;
  hyperedge_dimensions_ = hyperedge_dimension_arg;
  // local parameters
  std::string hypergraph_file = hypergraph_file_arg;
  std::string fixed_file = fixed_file_arg;
  logger_->report("Partition Parameters**");
  logger_->report("Number of partitions = {}", num_parts_);
  logger_->report("UBfactor = {}", ub_factor_);
  logger_->report("Vertex dimensions = {}", vertex_dimensions_);
  logger_->report("Hyperedge dimensions = {}", hyperedge_dimensions_);
  // build hypergraph
  ReadHypergraph(hypergraph_file, fixed_file);
  BuildHypergraph();
  logger_->report("Hypergraph Information**");
  logger_->report("#Vertices = {}", num_vertices_);
  logger_->report("#Hyperedges = {}", num_hyperedges_);
  // process hypergraph
  HGraph hypergraph_processed = preProcessHypergraph();
  logger_->report("Post processing hypergraph information**");
  logger_->report("#Vertices = {}", hypergraph_processed->GetNumVertices());
  logger_->report("#Hyperedges = {}", hypergraph_processed->GetNumHyperedges());
  // create coarsening class
  const std::vector<float> e_wt_factors(hyperedge_dimensions_, 1.0);
  const std::vector<float> v_wt_factors(vertex_dimensions_, 1.0);
  const std::vector<float> p_wt_factors(placement_dimensions_ + 100, 1.0);
  const float timing_factor = 1.0;
  const int path_traverse_step = 2;
  const std::vector<float> tot_vertex_weights
      = hypergraph_->GetTotalVertexWeights();
  const int alpha = 4;
  const std::vector<float> max_vertex_weights
      = DivideFactor(hypergraph_->GetTotalVertexWeights(), alpha * num_parts_);
  const int thr_coarsen_hyperedge_size = 50;
  global_net_threshold_ = 500;
  const int thr_coarsen_vertices = 200;
  const int thr_coarsen_hyperedges = 50;
  const float coarsening_ratio = 1.5;
  const int max_coarsen_iters = 20;
  const float adj_diff_ratio = 0.0001;
  const int vertex_ordering_choice = RANDOM;
  TP_coarsening_ptr tritonpart_coarsener
      = std::make_shared<TPcoarsener>(e_wt_factors,
                                      v_wt_factors,
                                      p_wt_factors,
                                      timing_factor,
                                      path_traverse_step,
                                      max_vertex_weights,
                                      thr_coarsen_hyperedge_size,
                                      global_net_threshold_,
                                      thr_coarsen_vertices,
                                      thr_coarsen_hyperedges,
                                      coarsening_ratio,
                                      max_coarsen_iters,
                                      adj_diff_ratio,
                                      seed_,
                                      logger_);
  float path_wt_factor = 1.0;
  float snaking_wt_factor = 1.0;
  const int refiner_iters = 2;
  const int max_moves = 50;
  int refiner_choice = TWO_WAY_FM;
  TP_two_way_refining_ptr tritonpart_twoway_refiner
      = std::make_shared<TPtwoWayFM>(num_parts_,
                                     refiner_iters,
                                     max_moves,
                                     refiner_choice,
                                     seed_,
                                     e_wt_factors,
                                     path_wt_factor,
                                     snaking_wt_factor,
                                     logger_);
  const int greedy_refiner_iters = 2;
  const int greedy_max_moves = 10;
  refiner_choice = GREEDY;
  TP_greedy_refiner_ptr tritonpart_greedy_refiner
      = std::make_shared<TPgreedyRefine>(num_parts_,
                                         greedy_refiner_iters,
                                         greedy_max_moves,
                                         refiner_choice,
                                         seed_,
                                         e_wt_factors,
                                         path_wt_factor,
                                         snaking_wt_factor,
                                         logger_);
  int wavefront = 50;
  int he_thr = 50;
  TP_ilp_refiner_ptr tritonpart_ilp_refiner
      = std::make_shared<TPilpRefine>(num_parts_,
                                      greedy_refiner_iters,
                                      greedy_max_moves,
                                      refiner_choice,
                                      seed_,
                                      e_wt_factors,
                                      path_wt_factor,
                                      snaking_wt_factor,
                                      logger_,
                                      wavefront);
  float early_stop_ratio = 0.5;
  int max_num_fm_pass = 10;
  TP_partitioning_ptr tritonpart_partitioner
      = std::make_shared<TPpartitioner>(num_parts_,
                                        e_wt_factors,
                                        path_wt_factor,
                                        snaking_wt_factor,
                                        early_stop_ratio,
                                        max_num_fm_pass,
                                        seed_,
                                        tritonpart_twoway_refiner,
                                        logger_);
  bool v_cycle_flag = true;
  RefinerType refine_type = KPM_REFINEMENT;
  int num_initial_solutions = 50;       // number of initial random solutions
  int num_best_initial_solutions = 10;  // number of best initial solutions
  int num_ubfactor_delta = 5;  // allowing marginal imbalance to improve QoR
  int max_num_vcycle = 5;      // maximum number of vcycles
  TP_mlevel_partitioning_ptr tritonpart_mlevel_partitioner
      = std::make_shared<TPmultilevelPartitioner>(tritonpart_coarsener,
                                                  tritonpart_partitioner,
                                                  tritonpart_twoway_refiner,
                                                  tritonpart_greedy_refiner,
                                                  tritonpart_ilp_refiner,
                                                  num_parts_,
                                                  v_cycle_flag,
                                                  num_initial_solutions,
                                                  num_best_initial_solutions,
                                                  num_ubfactor_delta,
                                                  max_num_vcycle,
                                                  seed_,
                                                  ub_factor_,
                                                  refine_type,
                                                  logger_);
  bool vcycle = true;
  matrix<float> vertex_balance
      = hypergraph_->GetVertexBalance(num_parts_, ub_factor_);
  auto solution = tritonpart_mlevel_partitioner->Partition(
      hypergraph_, hypergraph_processed, vertex_balance, vcycle);
  std::string solution_file
      = hypergraph_file + std::string(".part.") + std::to_string(num_parts_);
  auto cut_pair
      = tritonpart_partitioner->GoldenEvaluator(hypergraph_, solution, true);
  WriteSolution(solution_file.c_str(), solution);
  auto end_timestamp_global = std::chrono::high_resolution_clock::now();
  double total_global_time
      = std::chrono::duration_cast<std::chrono::nanoseconds>(
            end_timestamp_global - start_time_stamp_global)
            .count();
  total_global_time *= 1e-9;
  logger_->report("Total runtime {}", total_global_time);
  logger_->report("Exiting TritonPart");
  exit(EXIT_SUCCESS);
}
}  // namespace par
