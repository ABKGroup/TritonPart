// OpenSTA, Static Timing Analyzer
// Copyright (c) 2023, Parallax Software, Inc.
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.

#include "WritePathSpice.hh"

#include <string>
#include <iostream>
#include <fstream>

#include "Debug.hh"
#include "Error.hh"
#include "Report.hh"
#include "StringUtil.hh"
#include "FuncExpr.hh"
#include "Units.hh"
#include "Sequential.hh"
#include "TableModel.hh"
#include "Liberty.hh"
#include "TimingArc.hh"
#include "PortDirection.hh"
#include "Network.hh"
#include "Graph.hh"
#include "Sdc.hh"
#include "DcalcAnalysisPt.hh"
#include "Parasitics.hh"
#include "PathAnalysisPt.hh"
#include "Path.hh"
#include "PathRef.hh"
#include "PathExpanded.hh"
#include "StaState.hh"
#include "Sim.hh"

namespace sta {

using std::ofstream;
using std::ifstream;
using std::max;

typedef Map<string, StringVector*> CellSpicePortNames;
typedef int Stage;
typedef Map<ParasiticNode*, int> ParasiticNodeMap;
typedef Map<LibertyPort*, LogicValue> LibertyPortLogicValues;

void
streamPrint(ofstream &stream,
	    const char *fmt,
	    ...) __attribute__((format (printf, 2, 3)));

////////////////////////////////////////////////////////////////

class WritePathSpice : public StaState
{
public:
  WritePathSpice(Path *path,
		 const char *spice_filename,
		 const char *subckt_filename,
		 const char *lib_subckt_filename,
		 const char *model_filename,
                 StdStringSet *off_path_pin_names,
		 const char *power_name,
		 const char *gnd_name,
		 const StaState *sta);
  ~WritePathSpice();
  void writeSpice();

private:
  void writeHeader();
  void writeStageInstances();
  void writeInputSource();
  void writeRampVoltSource(const Pin *pin,
			   const RiseFall *rf,
			   float slew,
			   float time,
			   int &volt_index);
  void writeStageSubckts();
  void writeInputStage(Stage stage);
  void writeMeasureStmts();
  void writeMeasureStmt(const Pin *pin);
  void writeGateStage(Stage stage);
  void writeSubcktInst(const Pin *input_pin);
  void writeSubcktInstVoltSrcs(Stage stage,
			       const Pin *input_pin,
			       int &volt_index,
			       LibertyPortLogicValues &port_values,
			       const Clock *clk,
			       DcalcAPIndex dcalc_ap_index);
  void writeStageParasitics(Stage stage);
  void writeStageParasiticNetwork(Pin *drvr_pin,
                                  Parasitic *parasitic,
                                  ParasiticAnalysisPt *parasitic_ap);
  void writeStagePiElmore(Pin *drvr_pin,
                          Parasitic *parasitic);
  void writeNullParasitics(Pin *drvr_pin);
  void writeSubckts();
  StdStringSet findPathCellnames();
  void findPathCellSubckts(StdStringSet &path_cell_names);
  void recordSpicePortNames(const char *cell_name,
			    StringVector &tokens);
  float maxTime();
  float pathMaxTime();
  const char *nodeName(ParasiticNode *node);
  void initNodeMap(const char *net_name);
  const char *spiceTrans(const RiseFall *rf);
  void writeMeasureDelayStmt(Stage stage,
			     Path *from_path,
			     Path *to_path);
  void writeMeasureSlewStmt(Stage stage,
			    Path *path);
  void gatePortValues(Stage stage,
		      // Return values.
		      LibertyPortLogicValues &port_values,
		      const Clock *&clk,
		      DcalcAPIndex &dcalc_ap_index);
  void regPortValues(Stage stage,
		     // Return values.
		     LibertyPortLogicValues &port_values,
		     const Clock *&clk,
		     DcalcAPIndex &dcalc_ap_index);
  void gatePortValues(const Instance *inst,
		      FuncExpr *expr,
		      LibertyPort *input_port,
		      // Return values.
		      LibertyPortLogicValues &port_values);
  void seqPortValues(Sequential *seq,
		     const RiseFall *rf,
		     // Return values.
		     LibertyPortLogicValues &port_values);
  void writeInputWaveform();
  void writeClkWaveform();
  void writeWaveformEdge(const RiseFall *rf,
			 float time,
			 float slew);
  void writeWaveformVoltSource(const Pin *pin,
                               DriverWaveform *drvr_waveform,
                               const RiseFall *rf,
                               float slew,
                               int &volt_index);
  void writeClkedStepSource(const Pin *pin,
			    const RiseFall *rf,
			    const Clock *clk,
			    DcalcAPIndex dcalc_ap_index,
			    int &volt_index);
  float clkWaveformTImeOffset(const Clock *clk);
  float findSlew(Path *path);
  float findSlew(Path *path,
		 const RiseFall *rf,
		 TimingArc *next_arc);
  float findSlew(Vertex *vertex,
		 const RiseFall *rf,
		 TimingArc *next_arc,
		 DcalcAPIndex dcalc_ap_index);
  LibertyPort *onePort(FuncExpr *expr);
  void writeVoltageSource(const char *inst_name,
			  const char *port_name,
			  float voltage,
			  int &volt_index);
  void writeVoltageSource(LibertyCell *cell,
			  const char *inst_name,
			  const char *subckt_port_name,
			  const char *pg_port_name,
			  float voltage,
			  int &volt_index);
  float slewAxisMinValue(TimingArc *arc);
  float pgPortVoltage(LibertyPgPort *pg_port);
  void writePrintStmt();
  float railToRailSlew(float slew,
                       const RiseFall *rf);

  // Stage "accessors".
  //
  //           stage
  //      |---------------|
  //        |\             |\   .
  // -------| >---/\/\/----| >---
  //  gate  |/ drvr    load|/
  //  input
  //
  // A path from an input port has no GateInputPath (the input port is the drvr).
  // Internally a stage index from stageFirst() to stageLast()
  // is turned into an index into path_expanded_.
  //
  Stage stageFirst();
  Stage stageLast();
  string stageName(Stage stage);
  int stageGateInputPathIndex(Stage stage);
  int stageDrvrPathIndex(Stage stage);
  int stageLoadPathIndex(Stage stage);
  PathRef *stageGateInputPath(Stage stage);
  PathRef *stageDrvrPath(Stage stage);
  PathRef *stageLoadPath(Stage stage);
  TimingArc *stageGateArc(Stage stage);
  TimingArc *stageWireArc(Stage stage);
  Edge *stageGateEdge(Stage stage);
  Edge *stageWireEdge(Stage stage);
  Pin *stageGateInputPin(Stage stage);
  Pin *stageDrvrPin(Stage stage);
  LibertyPort *stageGateInputPort(Stage stage);
  LibertyPort *stageDrvrPort(Stage stage);
  Pin *stageLoadPin(Stage stage);
  const char *stageGateInputPinName(Stage stage);
  const char *stageDrvrPinName(Stage stage);
  const char *stageLoadPinName(Stage stage);
  LibertyCell *stageLibertyCell(Stage stage);
  Instance *stageInstance(Stage stage);
  StdStringSet stageOffPathPinNames(Stage stage);

  Path *path_;
  const char *spice_filename_;
  const char *subckt_filename_;
  const char *lib_subckt_filename_;
  const char *model_filename_;
  StdStringSet *off_path_pin_names_;
  const char *power_name_;
  const char *gnd_name_;

  ofstream spice_stream_;
  PathExpanded path_expanded_;
  CellSpicePortNames cell_spice_port_names_;
  ParasiticNodeMap node_map_;
  int next_node_index_;
  const char *net_name_;
  float power_voltage_;
  float gnd_voltage_;
  LibertyLibrary *default_library_;
  // Resistance to use to simulate a short circuit between spice nodes.
  float short_ckt_resistance_;
  // Input clock waveform cycles.
  int clk_cycle_count_;
};

////////////////////////////////////////////////////////////////

class SubcktEndsMissing : public Exception
{
public:
  SubcktEndsMissing(const char *cell_name,
		    const char *subckt_filename);
  const char *what() const noexcept;

protected:
  string what_;
};

SubcktEndsMissing::SubcktEndsMissing(const char *cell_name,
				     const char *subckt_filename) :
  Exception()
{
  what_ = "spice subckt for cell ";
  what_ += cell_name;
  what_ += " missing .ends in ";
  what_ += subckt_filename;
}

const char *
SubcktEndsMissing::what() const noexcept
{
  return what_.c_str();
}

////////////////////////////////////////////////////////////////

void
writePathSpice(Path *path,
	       const char *spice_filename,
	       const char *subckt_filename,
	       const char *lib_subckt_filename,
	       const char *model_filename,
               StdStringSet *off_path_pin_names,
               const char *power_name,
	       const char *gnd_name,
	       StaState *sta)
{
  WritePathSpice writer(path, spice_filename, subckt_filename,
			lib_subckt_filename, model_filename,
			off_path_pin_names, power_name, gnd_name, sta);
  writer.writeSpice();
}

WritePathSpice::WritePathSpice(Path *path,
			       const char *spice_filename,
			       const char *subckt_filename,
			       const char *lib_subckt_filename,
			       const char *model_filename,
                               StdStringSet *off_path_pin_names,
			       const char *power_name,
			       const char *gnd_name,
			       const StaState *sta) :
  StaState(sta),
  path_(path),
  spice_filename_(spice_filename),
  subckt_filename_(subckt_filename),
  lib_subckt_filename_(lib_subckt_filename),
  model_filename_(model_filename),
  off_path_pin_names_(off_path_pin_names),
  power_name_(power_name),
  gnd_name_(gnd_name),
  path_expanded_(sta),
  net_name_(nullptr),
  default_library_(network_->defaultLibertyLibrary()),
  short_ckt_resistance_(.0001),
  clk_cycle_count_(3)
{
  bool exists;
  default_library_->supplyVoltage(power_name_, power_voltage_, exists);
  if (!exists) {
    DcalcAnalysisPt *dcalc_ap = path_->dcalcAnalysisPt(this);
    const OperatingConditions *op_cond = dcalc_ap->operatingConditions();
    if (op_cond == nullptr)
      op_cond = network_->defaultLibertyLibrary()->defaultOperatingConditions();
    power_voltage_ = op_cond->voltage();
  }
  default_library_->supplyVoltage(gnd_name_, gnd_voltage_, exists);
  if (!exists)
    gnd_voltage_ = 0.0;
}

WritePathSpice::~WritePathSpice()
{
  stringDelete(net_name_);
  cell_spice_port_names_.deleteContents();
}

void
WritePathSpice::writeSpice()
{
  spice_stream_.open(spice_filename_);
  if (spice_stream_.is_open()) {
    path_expanded_.expand(path_, true);
    // Find subckt port names as a side-effect of writeSubckts.
    writeSubckts();
    writeHeader();
    writePrintStmt();
    writeStageInstances();
    writeMeasureStmts();
    writeInputSource();
    writeStageSubckts();
    streamPrint(spice_stream_, ".end\n");
    spice_stream_.close();
  }
  else
    throw FileNotWritable(spice_filename_);
}

// Use c++17 fs::path(filename).stem()
static string
filenameStem(const char *filename)
{
  string filename1 = filename;
  const size_t last_slash_idx = filename1.find_last_of("\\/");
  if (last_slash_idx != std::string::npos)
    return filename1.substr(last_slash_idx + 1);
  else
    return filename1;
}

void
WritePathSpice::writeHeader()
{
  const MinMax *min_max = path_->minMax(this);
  Pvt *pvt = sdc_->operatingConditions(min_max);
  if (pvt == nullptr)
    pvt = default_library_->defaultOperatingConditions();
  Path *start_path = path_expanded_.startPath();
  streamPrint(spice_stream_, "* Path from %s %s to %s %s\n",
	      network_->pathName(start_path->pin(this)),
	      start_path->transition(this)->asString(),
	      network_->pathName(path_->pin(this)),
	      path_->transition(this)->asString());
  float temp = pvt->temperature();
  streamPrint(spice_stream_, ".temp %.1f\n", temp);
  streamPrint(spice_stream_, ".include \"%s\"\n", model_filename_);
  string subckt_filename_stem = filenameStem(subckt_filename_);
  streamPrint(spice_stream_, ".include \"%s\"\n", subckt_filename_stem.c_str());

  float max_time = maxTime();
  float time_step = 1e-13;
  streamPrint(spice_stream_, ".tran %.3g %.3g\n\n",
	      time_step, max_time);
  streamPrint(spice_stream_, ".options nomod\n");
}

void
WritePathSpice::writePrintStmt()
{
  streamPrint(spice_stream_, ".print tran");
  for (Stage stage = stageFirst(); stage <= stageLast(); stage++) {
    streamPrint(spice_stream_, " v(%s)", stageDrvrPinName(stage));
    streamPrint(spice_stream_, " v(%s)", stageLoadPinName(stage));
    StdStringSet off_path_names = stageOffPathPinNames(stage);
    for (const string &off_path_name : off_path_names)
      streamPrint(spice_stream_, " v(%s)", off_path_name.c_str());
  }
  streamPrint(spice_stream_, "\n\n");
}

float
WritePathSpice::maxTime()
{
  Stage input_stage = stageFirst();
  PathRef *input_path = stageDrvrPath(input_stage);
  if (input_path->isClock(this)) {
    const Clock *clk = input_path->clock(this);
    float period = clk->period();
    float first_edge_offset = period / 10;
    float max_time = period * clk_cycle_count_ + first_edge_offset;
    return max_time;
  }
  else
    return pathMaxTime();
}

// Make sure run time is long enough to see side load transitions along the path.
float
WritePathSpice::pathMaxTime()
{
  float max_time = 0.0;
  DcalcAPIndex dcalc_ap_index = path_->dcalcAnalysisPt(this)->index();
  for (size_t i = 0; i < path_expanded_.size(); i++) {
    PathRef *path = path_expanded_.path(i);
    const RiseFall *rf = path->transition(this);
    Vertex *vertex = path->vertex(this);
    float path_max_slew = railToRailSlew(findSlew(vertex,rf,nullptr,dcalc_ap_index),rf);
    if (vertex->isDriver(network_)) {
      VertexOutEdgeIterator edge_iter(vertex, graph_);
      while (edge_iter.hasNext()) {
        Edge *edge = edge_iter.next();
        Vertex *load = edge->to(graph_);
        float load_slew = railToRailSlew(findSlew(load, rf, nullptr, dcalc_ap_index),rf);
        if (load_slew > path_max_slew)
          path_max_slew = load_slew;
      }
    }
    float path_max_time = delayAsFloat(path->arrival(this)) + path_max_slew * 2.0;
    if (path_max_time > max_time)
      max_time = path_max_time;
  }
  return max_time;
}

float
WritePathSpice::railToRailSlew(float slew,
                               const RiseFall *rf)
{
  float lower = default_library_->slewLowerThreshold(rf);
  float upper = default_library_->slewUpperThreshold(rf);
  return slew / (upper - lower);
}

void
WritePathSpice::writeStageInstances()
{
  streamPrint(spice_stream_, "*****************\n");
  streamPrint(spice_stream_, "* Stage instances\n");
  streamPrint(spice_stream_, "*****************\n\n");

  for (Stage stage = stageFirst(); stage <= stageLast(); stage++) {
    string stage_name = stageName(stage);
    const char *stage_cname = stage_name.c_str();
    if (stage == stageFirst())
      streamPrint(spice_stream_, "x%s %s %s %s\n",
		  stage_cname,
		  stageDrvrPinName(stage),
		  stageLoadPinName(stage),
		  stage_cname);
    else {
      streamPrint(spice_stream_, "x%s %s %s %s",
		  stage_cname,
		  stageGateInputPinName(stage),
		  stageDrvrPinName(stage),
		  stageLoadPinName(stage));
      StdStringSet off_path_names = stageOffPathPinNames(stage);
      for (const string &off_path_name : off_path_names)
        streamPrint(spice_stream_, " %s", off_path_name.c_str());
      streamPrint(spice_stream_, " %s\n", stage_cname);
    }
  }
  streamPrint(spice_stream_, "\n");
}

float
WritePathSpice::pgPortVoltage(LibertyPgPort *pg_port)
{
  LibertyLibrary *liberty = pg_port->cell()->libertyLibrary();
  float voltage = 0.0;
  bool exists;
  const char *voltage_name = pg_port->voltageName();
  if (voltage_name) {
    liberty->supplyVoltage(voltage_name, voltage, exists);
    if (!exists) {
      if (stringEqual(voltage_name, power_name_))
	voltage = power_voltage_;
      else if (stringEqual(voltage_name, gnd_name_))
	voltage = gnd_voltage_;
      else
	report_->error(24, "pg_pin %s/%s voltage %s not found,",
		       pg_port->cell()->name(),
		       pg_port->name(),
		       voltage_name);
    }
  }
  else
    report_->error(25, "Liberty pg_port %s/%s missing voltage_name attribute,",
		   pg_port->cell()->name(),
		   pg_port->name());
  return voltage;
}

void
WritePathSpice::writeInputSource()
{
  streamPrint(spice_stream_, "**************\n");
  streamPrint(spice_stream_, "* Input source\n");
  streamPrint(spice_stream_, "**************\n\n");

  Stage input_stage = stageFirst();
  PathRef *input_path = stageDrvrPath(input_stage);
  if (input_path->isClock(this))
    writeClkWaveform();
  else
    writeInputWaveform();
  streamPrint(spice_stream_, "\n");
}

void
WritePathSpice::writeInputWaveform()
{
  Stage input_stage = stageFirst();
  PathRef *input_path = stageDrvrPath(input_stage);
  const RiseFall *rf = input_path->transition(this);
  TimingArc *next_arc = stageGateArc(input_stage + 1);
  float slew0 = findSlew(input_path, rf, next_arc);

  float threshold = default_library_->inputThreshold(rf);
  float dt = railToRailSlew(slew0, rf);
  float time0 = dt * threshold;

  int volt_index = 1;
  const Pin *drvr_pin = stageDrvrPin(input_stage);
  const Pin *load_pin = stageLoadPin(input_stage);
  const LibertyPort *load_port = network_->libertyPort(load_pin);
  DriverWaveform *drvr_waveform = nullptr;
  if (load_port)
    drvr_waveform = load_port->driverWaveform(rf);
  if (drvr_waveform)
    writeWaveformVoltSource(drvr_pin, drvr_waveform,
                            rf, slew0, volt_index);
  else
    writeRampVoltSource(drvr_pin, rf, slew0, time0, volt_index);
}

void
WritePathSpice::writeWaveformVoltSource(const Pin *pin,
                                        DriverWaveform *drvr_waveform,
                                        const RiseFall *rf,
                                        float slew,
                                        int &volt_index)
{
  float volt0, volt1, volt_factor;
  if (rf == RiseFall::rise()) {
    volt0 = gnd_voltage_;
    volt1 = power_voltage_;
    volt_factor = power_voltage_;
  }
  else {
    volt0 = power_voltage_;
    volt1 = gnd_voltage_;
    volt_factor = -power_voltage_;
  }
  streamPrint(spice_stream_, "v%d %s 0 pwl(\n",
	      volt_index,
	      network_->pathName(pin));
  streamPrint(spice_stream_, "+%.3e %.3e\n", 0.0, volt0);
  Table1 waveform = drvr_waveform->waveform(slew);
  const TableAxis *time_axis = waveform.axis1();
  for (size_t time_index = 0; time_index <  time_axis->size(); time_index++) {
    float time = time_axis->axisValue(time_index);
    float wave_volt = waveform.value(time_index);
    float volt = volt0 + wave_volt * volt_factor;
    streamPrint(spice_stream_, "+%.3e %.3e\n", time, volt);
  }
  streamPrint(spice_stream_, "+%.3e %.3e\n", maxTime(), volt1);
  streamPrint(spice_stream_, "+)\n");
  volt_index++;
}

void
WritePathSpice::writeRampVoltSource(const Pin *pin,
				    const RiseFall *rf,
				    float slew,
				    float time,
				    int &volt_index)
{
  float volt0, volt1;
  if (rf == RiseFall::rise()) {
    volt0 = gnd_voltage_;
    volt1 = power_voltage_;
  }
  else {
    volt0 = power_voltage_;
    volt1 = gnd_voltage_;
  }
  streamPrint(spice_stream_, "v%d %s 0 pwl(\n",
	      volt_index,
	      network_->pathName(pin));
  streamPrint(spice_stream_, "+%.3e %.3e\n", 0.0, volt0);
  writeWaveformEdge(rf, time, slew);
  streamPrint(spice_stream_, "+%.3e %.3e\n", maxTime(), volt1);
  streamPrint(spice_stream_, "+)\n");
  volt_index++;
}

void
WritePathSpice::writeClkWaveform()
{
  Stage input_stage = stageFirst();
  PathRef *input_path = stageDrvrPath(input_stage);
  TimingArc *next_arc = stageGateArc(input_stage + 1);
  const ClockEdge *clk_edge = input_path->clkEdge(this);
  const Clock *clk = clk_edge->clock();
  float period = clk->period();
  float time_offset = clkWaveformTImeOffset(clk);
  RiseFall *rf0, *rf1;
  float volt0;
  if (clk_edge->time() < period) {
    rf0 = RiseFall::rise();
    rf1 = RiseFall::fall();
    volt0 = gnd_voltage_;
  }
  else {
    rf0 = RiseFall::fall();
    rf1 = RiseFall::rise();
    volt0 = power_voltage_;
  }
  float slew0 = findSlew(input_path, rf0, next_arc);
  float slew1 = findSlew(input_path, rf1, next_arc);
  streamPrint(spice_stream_, "v1 %s 0 pwl(\n",
	      stageDrvrPinName(input_stage));
  streamPrint(spice_stream_, "+%.3e %.3e\n", 0.0, volt0);
  for (int cycle = 0; cycle < clk_cycle_count_; cycle++) {
    float time0 = time_offset + cycle * period;
    float time1 = time0 + period / 2.0;
    writeWaveformEdge(rf0, time0, slew0);
    writeWaveformEdge(rf1, time1, slew1);
  }
  streamPrint(spice_stream_, "+%.3e %.3e\n", maxTime(), volt0);
  streamPrint(spice_stream_, "+)\n");
}

float
WritePathSpice::clkWaveformTImeOffset(const Clock *clk)
{
  return clk->period() / 10;
}

float
WritePathSpice::findSlew(Path *path)
{
  Vertex *vertex = path->vertex(this);
  DcalcAPIndex dcalc_ap_index = path->dcalcAnalysisPt(this)->index();
  const RiseFall *rf = path->transition(this);
  return findSlew(vertex, rf, nullptr, dcalc_ap_index);
}

float
WritePathSpice::findSlew(Path *path,
			 const RiseFall *rf,
			 TimingArc *next_arc)
{
  Vertex *vertex = path->vertex(this);
  DcalcAPIndex dcalc_ap_index = path->dcalcAnalysisPt(this)->index();
  return findSlew(vertex, rf, next_arc, dcalc_ap_index);
}

float
WritePathSpice::findSlew(Vertex *vertex,
			 const RiseFall *rf,
			 TimingArc *next_arc,
			 DcalcAPIndex dcalc_ap_index)
{
  float slew = delayAsFloat(graph_->slew(vertex, rf, dcalc_ap_index));
  if (slew == 0.0 && next_arc)
    slew = slewAxisMinValue(next_arc);
  if (slew == 0.0)
    slew = units_->timeUnit()->scale();
  return slew;
}

// Look up the smallest slew axis value in the timing arc delay table.
float
WritePathSpice::slewAxisMinValue(TimingArc *arc)
{
  GateTableModel *gate_model = dynamic_cast<GateTableModel*>(arc->model());
  if (gate_model) {
    const TableModel *model = gate_model->delayModel();
    const TableAxis *axis1 = model->axis1();
    TableAxisVariable var1 = axis1->variable();
    if (var1 == TableAxisVariable::input_transition_time
	|| var1 == TableAxisVariable::input_net_transition)
      return axis1->axisValue(0);

    const TableAxis *axis2 = model->axis2();
    TableAxisVariable var2 = axis2->variable();
    if (var2 == TableAxisVariable::input_transition_time
	|| var2 == TableAxisVariable::input_net_transition)
      return axis2->axisValue(0);

    const TableAxis *axis3 = model->axis3();
    TableAxisVariable var3 = axis3->variable();
    if (var3 == TableAxisVariable::input_transition_time
	|| var3 == TableAxisVariable::input_net_transition)
      return axis3->axisValue(0);
  }
  return 0.0;
}

// Write PWL rise/fall edge that crosses threshold at time.
void
WritePathSpice::writeWaveformEdge(const RiseFall *rf,
				  float time,
				  float slew)
{
  float volt0, volt1;
  if (rf == RiseFall::rise()) {
    volt0 = gnd_voltage_;
    volt1 = power_voltage_;
  }
  else {
    volt0 = power_voltage_;
    volt1 = gnd_voltage_;
  }
  float threshold = default_library_->inputThreshold(rf);
  float dt = railToRailSlew(slew, rf);
  float time0 = time - dt * threshold;
  float time1 = time0 + dt;
  if (time0 > 0.0)
    streamPrint(spice_stream_, "+%.3e %.3e\n", time0, volt0);
  streamPrint(spice_stream_, "+%.3e %.3e\n", time1, volt1);
}

////////////////////////////////////////////////////////////////

void
WritePathSpice::writeMeasureStmts()
{
  streamPrint(spice_stream_, "********************\n");
  streamPrint(spice_stream_, "* Measure statements\n");
  streamPrint(spice_stream_, "********************\n\n");

  for (Stage stage = stageFirst(); stage <= stageLast(); stage++) {
    PathRef *gate_input_path = stageGateInputPath(stage);
    PathRef *drvr_path = stageDrvrPath(stage);
    PathRef *load_path = stageLoadPath(stage);
    if (gate_input_path) {
      // gate input -> gate output
      writeMeasureSlewStmt(stage, gate_input_path);
      writeMeasureDelayStmt(stage, gate_input_path, drvr_path);
    }
    writeMeasureSlewStmt(stage, drvr_path);
    // gate output | input port -> load
    writeMeasureDelayStmt(stage, drvr_path, load_path);
    if (stage == stageLast())
      writeMeasureSlewStmt(stage, load_path);
  }
  streamPrint(spice_stream_, "\n");
}

void
WritePathSpice::writeMeasureDelayStmt(Stage stage,
				      Path *from_path,
				      Path *to_path)
{
  const char *from_pin_name = network_->pathName(from_path->pin(this));
  const RiseFall *from_rf = from_path->transition(this);
  float from_threshold = power_voltage_ * default_library_->inputThreshold(from_rf);

  const char *to_pin_name = network_->pathName(to_path->pin(this));
  const RiseFall *to_rf = to_path->transition(this);
  float to_threshold = power_voltage_ * default_library_->inputThreshold(to_rf);
  string stage_name = stageName(stage);
  streamPrint(spice_stream_,
	      ".measure tran %s_%s_delay_%s\n",
	      stage_name.c_str(),
	      from_pin_name,
	      to_pin_name);
  streamPrint(spice_stream_,
	      "+trig v(%s) val=%.3f %s=last\n",
	      from_pin_name,
	      from_threshold,
	      spiceTrans(from_rf));
  streamPrint(spice_stream_,
	      "+targ v(%s) val=%.3f %s=last\n",
	      to_pin_name,
	      to_threshold,
	      spiceTrans(to_rf));
}

void
WritePathSpice::writeMeasureSlewStmt(Stage stage,
				     Path *path)
{
  const char *pin_name = network_->pathName(path->pin(this));
  const RiseFall *rf = path->transition(this);
  const char *spice_rf = spiceTrans(rf);
  float lower = power_voltage_ * default_library_->slewLowerThreshold(rf);
  float upper = power_voltage_ * default_library_->slewUpperThreshold(rf);
  float threshold1, threshold2;
  if (rf == RiseFall::rise()) {
    threshold1 = lower;
    threshold2 = upper;
  }
  else {
    threshold1 = upper;
    threshold2 = lower;
  }
  string stage_name = stageName(stage);
  streamPrint(spice_stream_,
	      ".measure tran %s_%s_slew\n",
	      stage_name.c_str(),
	      pin_name);
  streamPrint(spice_stream_,
	      "+trig v(%s) val=%.3f %s=last\n",
	      pin_name,
	      threshold1,
	      spice_rf);
  streamPrint(spice_stream_,
	      "+targ v(%s) val=%.3f %s=last\n",
	      pin_name,
	      threshold2,
	      spice_rf);
}

const char *
WritePathSpice::spiceTrans(const RiseFall *rf)
{
  if (rf == RiseFall::rise())
    return "RISE";
  else
    return "FALL";
}

void
WritePathSpice::writeStageSubckts()
{
  streamPrint(spice_stream_, "***************\n");
  streamPrint(spice_stream_, "* Stage subckts\n");
  streamPrint(spice_stream_, "***************\n\n");

  for (Stage stage = stageFirst(); stage <= stageLast(); stage++) {
    if (stage == stageFirst())
      writeInputStage(stage);
    else
      writeGateStage(stage);
  }
}

// Input port to first gate input.
void
WritePathSpice::writeInputStage(Stage stage)
{
  // Input arc.
  // External driver not handled.
  const char *drvr_pin_name = stageDrvrPinName(stage);
  const char *load_pin_name = stageLoadPinName(stage);
  string stage_name = stageName(stage);
  streamPrint(spice_stream_, ".subckt %s %s %s\n",
	      stage_name.c_str(),
	      drvr_pin_name,
	      load_pin_name);
  writeStageParasitics(stage);
  streamPrint(spice_stream_, ".ends\n\n");
}

// Gate and load parasitics.
void
WritePathSpice::writeGateStage(Stage stage)
{
  const Pin *input_pin = stageGateInputPin(stage);
  const char *input_pin_name = stageGateInputPinName(stage);
  const Pin *drvr_pin = stageDrvrPin(stage);
  const char *drvr_pin_name = stageDrvrPinName(stage);
  const Pin *load_pin = stageLoadPin(stage);
  const char *load_pin_name = stageLoadPinName(stage);
  streamPrint(spice_stream_, ".subckt stage%d %s %s %s",
	      stage,
	      input_pin_name,
	      drvr_pin_name,
	      load_pin_name);
  StdStringSet off_path_names = stageOffPathPinNames(stage);
  for (const string &off_path_name : off_path_names)
    streamPrint(spice_stream_, " %s", off_path_name.c_str());
  streamPrint(spice_stream_, "\n");

  // Driver subckt call.
  Instance *inst = stageInstance(stage);
  LibertyPort *input_port = stageGateInputPort(stage);
  LibertyPort *drvr_port = stageDrvrPort(stage);
  streamPrint(spice_stream_, "* Gate %s %s -> %s\n",
	      network_->pathName(inst),
	      input_port->name(),
	      drvr_port->name());
  writeSubcktInst(input_pin);
  LibertyPortLogicValues port_values;
  DcalcAPIndex dcalc_ap_index;
  const Clock *clk;
  int volt_index = 1;
  gatePortValues(stage, port_values, clk, dcalc_ap_index);
  writeSubcktInstVoltSrcs(stage, input_pin, volt_index,
			  port_values, clk, dcalc_ap_index);
  streamPrint(spice_stream_, "\n");

  port_values.clear();
  auto pin_iter = network_->connectedPinIterator(drvr_pin);
  while (pin_iter->hasNext()) {
    const Pin *pin = pin_iter->next();
    if (pin != drvr_pin
	&& pin != load_pin
	&& network_->direction(pin)->isAnyInput()
	&& !network_->isHierarchical(pin)
	&& !network_->isTopLevelPort(pin)) {
      streamPrint(spice_stream_, "* Side load %s\n", network_->pathName(pin));
      writeSubcktInst(pin);
      writeSubcktInstVoltSrcs(stage, pin, volt_index, port_values, nullptr, 0); 
      streamPrint(spice_stream_, "\n");
    }
  }
  delete pin_iter;

  writeStageParasitics(stage);
  streamPrint(spice_stream_, ".ends\n\n");
}

void
WritePathSpice::writeSubcktInst(const Pin *input_pin)
{
  const Instance *inst = network_->instance(input_pin);
  const char *inst_name = network_->pathName(inst);
  LibertyCell *cell = network_->libertyCell(inst);
  const char *cell_name = cell->name();
  StringVector *spice_port_names = cell_spice_port_names_[cell_name];
  streamPrint(spice_stream_, "x%s", inst_name);
  for (string subckt_port_name : *spice_port_names) {
    const char *subckt_port_cname = subckt_port_name.c_str();
    Pin *pin = network_->findPin(inst, subckt_port_cname);
    LibertyPgPort *pg_port = cell->findPgPort(subckt_port_cname);
    const char *pin_name;
    if (pin) {
      pin_name = network_->pathName(pin);
      streamPrint(spice_stream_, " %s", pin_name);
    }
    else if (pg_port)
      streamPrint(spice_stream_, " %s/%s", inst_name, subckt_port_cname);
    else if (stringEq(subckt_port_cname, power_name_)
	     || stringEq(subckt_port_cname, gnd_name_))
      streamPrint(spice_stream_, " %s/%s", inst_name, subckt_port_cname);
  }
  streamPrint(spice_stream_, " %s\n", cell_name);
}

// Power/ground and input voltage sources.
void
WritePathSpice::writeSubcktInstVoltSrcs(Stage stage,
					const Pin *input_pin,
					int &volt_index,
					LibertyPortLogicValues &port_values,
					const Clock *clk,
					DcalcAPIndex dcalc_ap_index)

{
  const Instance *inst = network_->instance(input_pin);
  LibertyCell *cell = network_->libertyCell(inst);
  const char *cell_name = cell->name();
  StringVector *spice_port_names = cell_spice_port_names_[cell_name];

  const Pin *drvr_pin = stageDrvrPin(stage);
  const LibertyPort *input_port = network_->libertyPort(input_pin);
  const LibertyPort *drvr_port = network_->libertyPort(drvr_pin);
  const char *input_port_name = input_port->name();
  const char *drvr_port_name = drvr_port->name();
  const char *inst_name = network_->pathName(inst);

  debugPrint(debug_, "write_spice", 2, "subckt %s", cell->name());
  for (string subckt_port_sname : *spice_port_names) {
    const char *subckt_port_name = subckt_port_sname.c_str();
    LibertyPgPort *pg_port = cell->findPgPort(subckt_port_name);
    debugPrint(debug_, "write_spice", 2, " port %s%s",
               subckt_port_name,
               pg_port ? " pwr/gnd" : "");
    if (pg_port)
      writeVoltageSource(inst_name, subckt_port_name,
			 pgPortVoltage(pg_port), volt_index);
    else if (stringEq(subckt_port_name, power_name_))
      writeVoltageSource(inst_name, subckt_port_name,
			 power_voltage_, volt_index);
    else if (stringEq(subckt_port_name, gnd_name_))
      writeVoltageSource(inst_name, subckt_port_name,
			 gnd_voltage_, volt_index);
    else if (!(stringEq(subckt_port_name, input_port_name)
		 || stringEq(subckt_port_name, drvr_port_name))) {
      // Input voltage to sensitize path from gate input to output.
      LibertyPort *port = cell->findLibertyPort(subckt_port_name);
      if (port
	  && port->direction()->isAnyInput()) {
	const Pin *pin = network_->findPin(inst, port);
	// Look for tie high/low or propagated constant values.
	LogicValue port_value = sim_->logicValue(pin);
	if (port_value == LogicValue::unknown) {
	  bool has_value;
	  LogicValue value;
	  port_values.findKey(port, value, has_value);
	  if (has_value)
	    port_value = value;
	}
	switch (port_value) {
	case LogicValue::zero:
	case LogicValue::unknown:
	  writeVoltageSource(cell, inst_name, subckt_port_name,
			     port->relatedGroundPin(),
			     gnd_voltage_,
			     volt_index);
	  break;
	case LogicValue::one:
	  writeVoltageSource(cell, inst_name, subckt_port_name,
			     port->relatedPowerPin(),
			     power_voltage_,
			     volt_index);
	  break;
	case LogicValue::rise:
	  writeClkedStepSource(pin, RiseFall::rise(), clk,
			       dcalc_ap_index, volt_index);
	  break;
	case LogicValue::fall:
	  writeClkedStepSource(pin, RiseFall::fall(), clk,
			       dcalc_ap_index, volt_index);
	  break;
	}
      }
    }
  }
}

// PWL voltage source that rises half way into the first clock cycle.
void
WritePathSpice::writeClkedStepSource(const Pin *pin,
				     const RiseFall *rf,
				     const Clock *clk,
				     DcalcAPIndex dcalc_ap_index,
				     int &volt_index)
{
  Vertex *vertex = graph_->pinLoadVertex(pin);
  float slew = findSlew(vertex, rf, nullptr, dcalc_ap_index);
  float time = clkWaveformTImeOffset(clk) + clk->period() / 2.0;
  writeRampVoltSource(pin, rf, slew, time, volt_index);
}

void
WritePathSpice::writeVoltageSource(const char *inst_name,
				   const char *port_name,
				   float voltage,
				   int &volt_index)
{
  streamPrint(spice_stream_, "v%d %s/%s 0 %.3f\n",
	      volt_index,
	      inst_name, port_name,
	      voltage);
  volt_index++;
}

void
WritePathSpice::writeVoltageSource(LibertyCell *cell,
				   const char *inst_name,
				   const char *subckt_port_name,
				   const char *pg_port_name,
				   float voltage,
				   int &volt_index)
{
  if (pg_port_name) {
    LibertyPgPort *pg_port = cell->findPgPort(pg_port_name);
    if (pg_port)
      voltage = pgPortVoltage(pg_port);
    else
      report_->error(26, "%s pg_port %s not found,",
		     cell->name(),
		     pg_port_name);

  }
  writeVoltageSource(inst_name, subckt_port_name, voltage, volt_index);
}

void
WritePathSpice::gatePortValues(Stage stage,
			       // Return values.
			       LibertyPortLogicValues &port_values,
			       const Clock *&clk,
			       DcalcAPIndex &dcalc_ap_index)
{
  clk = nullptr;
  dcalc_ap_index = 0;

  Edge *gate_edge = stageGateEdge(stage);
  LibertyPort *drvr_port = stageDrvrPort(stage);
  if (gate_edge
      && gate_edge->role()->genericRole() == TimingRole::regClkToQ()) 
    regPortValues(stage, port_values, clk, dcalc_ap_index);
  else if (drvr_port->function()) {
    Pin *input_pin = stageGateInputPin(stage);
    LibertyPort *input_port = network_->libertyPort(input_pin);
    Instance *inst = network_->instance(input_pin);
    gatePortValues(inst, drvr_port->function(), input_port, port_values);
  }
}

void
WritePathSpice::regPortValues(Stage stage,
			      // Return values.
			      LibertyPortLogicValues &port_values,
			      const Clock *&clk,
			      DcalcAPIndex &dcalc_ap_index)
{
  LibertyPort *drvr_port = stageDrvrPort(stage);
  FuncExpr *drvr_expr = drvr_port->function();
  if (drvr_expr) {
    LibertyPort *q_port = drvr_expr->port();
    if (q_port) {
      // Drvr (register/latch output) function should be a reference
      // to an internal port like IQ or IQN.
      LibertyCell *cell = stageLibertyCell(stage);
      Sequential *seq = cell->outputPortSequential(q_port);
      if (seq) {
	PathRef *drvr_path = stageDrvrPath(stage);
	const RiseFall *drvr_rf = drvr_path->transition(this);
	seqPortValues(seq, drvr_rf, port_values);
	clk = drvr_path->clock(this);
	dcalc_ap_index = drvr_path->dcalcAnalysisPt(this)->index();
      }
      else
	report_->error(27, "no register/latch found for path from %s to %s,",
		       stageGateInputPort(stage)->name(),
		       stageDrvrPort(stage)->name());
    }
  }
}

// Find the logic values for expression inputs to enable paths input_port.
void
WritePathSpice::gatePortValues(const Instance *inst,
			       FuncExpr *expr,
			       LibertyPort *input_port,
			       // Return values.
			       LibertyPortLogicValues &port_values)
{
  FuncExpr *left = expr->left();
  FuncExpr *right = expr->right();
  switch (expr->op()) {
  case FuncExpr::op_port:
    break;
  case FuncExpr::op_not:
    gatePortValues(inst, left, input_port, port_values);
    break;
  case FuncExpr::op_or:
    if (left->hasPort(input_port)
	&& right->op() == FuncExpr::op_port)
      port_values[right->port()] = LogicValue::zero;
    else if (left->hasPort(input_port)
	     && right->op() == FuncExpr::op_not
	     && right->left()->op() == FuncExpr::op_port)
	// input_port + !right_port
	port_values[right->left()->port()] = LogicValue::one;
    else if (right->hasPort(input_port)
	     && left->op() == FuncExpr::op_port)
	port_values[left->port()] = LogicValue::zero;
    else if (right->hasPort(input_port)
	     && left->op() == FuncExpr::op_not
	     && left->left()->op() == FuncExpr::op_port)
      // input_port + !left_port
      port_values[left->left()->port()] = LogicValue::one;
    else {
      gatePortValues(inst, left, input_port, port_values);
      gatePortValues(inst, right, input_port, port_values);
    }
    break;
  case FuncExpr::op_and:
    if (left->hasPort(input_port)
	&& right->op() == FuncExpr::op_port)
	port_values[right->port()] = LogicValue::one;
    else if (left->hasPort(input_port)
	     && right->op() == FuncExpr::op_not
	     && right->left()->op() == FuncExpr::op_port)
      // input_port * !right_port
      port_values[right->left()->port()] = LogicValue::zero;
    else if (right->hasPort(input_port)
	     && left->op() == FuncExpr::op_port)
      port_values[left->port()] = LogicValue::one;
    else if (right->hasPort(input_port)
	     && left->op() == FuncExpr::op_not
	     && left->left()->op() == FuncExpr::op_port)
      // input_port * !left_port
      port_values[left->left()->port()] = LogicValue::zero;
    else {
      gatePortValues(inst, left, input_port, port_values);
      gatePortValues(inst, right, input_port, port_values);
    }
    break;
  case FuncExpr::op_xor:
    // Need to know timing arc sense to get this right.
    if (left->port() == input_port
	&& right->op() == FuncExpr::op_port)
      port_values[right->port()] = LogicValue::zero;
    else if (right->port() == input_port
	     && left->op() == FuncExpr::op_port)
      port_values[left->port()] = LogicValue::zero;
    else {
      gatePortValues(inst, left, input_port, port_values);
      gatePortValues(inst, right, input_port, port_values);
    }
    break;
  case FuncExpr::op_one:
  case FuncExpr::op_zero:
    break;
  }
}

void
WritePathSpice::seqPortValues(Sequential *seq,
			      const RiseFall *rf,
			      // Return values.
			      LibertyPortLogicValues &port_values)
{
  FuncExpr *data = seq->data();
  LibertyPort *port = onePort(data);
  if (port) {
    TimingSense sense = data->portTimingSense(port);
    switch (sense) {
    case TimingSense::positive_unate:
      if (rf == RiseFall::rise())
	port_values[port] = LogicValue::rise;
      else
	port_values[port] = LogicValue::fall;
      break;
    case TimingSense::negative_unate:
      if (rf == RiseFall::rise())
	port_values[port] = LogicValue::fall;
      else
	port_values[port] = LogicValue::rise;
      break;
    case TimingSense::non_unate:
    case TimingSense::none:
    case TimingSense::unknown:
    default:
      break;
    }
  }
}

// Pick a port, any port...
LibertyPort *
WritePathSpice::onePort(FuncExpr *expr)
{
  FuncExpr *left = expr->left();
  FuncExpr *right = expr->right();
  LibertyPort *port;
  switch (expr->op()) {
  case FuncExpr::op_port:
    return expr->port();
  case FuncExpr::op_not:
    return onePort(left);
  case FuncExpr::op_or:
  case FuncExpr::op_and:
  case FuncExpr::op_xor:
    port = onePort(left);
    if (port == nullptr)
      port = onePort(right);
    return port;
  case FuncExpr::op_one:
  case FuncExpr::op_zero:
  default:
    return nullptr;
  }
}

// Sort predicate for ParasiticDevices.
class ParasiticDeviceLess
{
public:
  ParasiticDeviceLess(Parasitics *parasitics) :
    parasitics_(parasitics)
  {
  }
  bool operator()(const ParasiticDevice *device1,
		  const ParasiticDevice *device2) const
  {
    ParasiticNode *node1 = parasitics_->node1(device1);
    ParasiticNode *node2 = parasitics_->node1(device2);
    const char *name1 = parasitics_->name(node1);
    const char *name2 = parasitics_->name(node2);
    if (stringEq(name1, name2)) {
      ParasiticNode *node12 = parasitics_->node2(device1);
      ParasiticNode *node22 = parasitics_->node2(device2);
      const char *name12 = parasitics_->name(node12);
      const char *name22 = parasitics_->name(node22);
      return stringLess(name12, name22);
    }
    else 
      return stringLess(name1, name2);
  }
private:
  Parasitics *parasitics_;
};

// Sort predicate for ParasiticDevices.
class ParasiticNodeLess
{
public:
  ParasiticNodeLess(Parasitics *parasitics) :
    parasitics_(parasitics)
  {
  }
  bool operator()(const ParasiticNode *node1,
		  const ParasiticNode *node2) const
  {
    const char *name1 = parasitics_->name(node1);
    const char *name2 = parasitics_->name(node2);
    return stringLess(name1, name2);
  }
private:
  Parasitics *parasitics_;
};

void
WritePathSpice::writeStageParasitics(Stage stage)
{
  PathRef *drvr_path = stageDrvrPath(stage);
  Pin *drvr_pin = stageDrvrPin(stage);

  Net *net = network_->net(drvr_pin);
  const char *net_name = net ? network_->pathName(net) : network_->pathName(drvr_pin);
  initNodeMap(net_name);
  streamPrint(spice_stream_, "* Net %s\n", net_name);

  DcalcAnalysisPt *dcalc_ap = drvr_path->dcalcAnalysisPt(this);
  ParasiticAnalysisPt *parasitic_ap = dcalc_ap->parasiticAnalysisPt();
  Parasitic *parasitic = parasitics_->findParasiticNetwork(drvr_pin, parasitic_ap);
  if (parasitic)
    writeStageParasiticNetwork(drvr_pin, parasitic, parasitic_ap);
  else {
    const RiseFall *drvr_rf = drvr_path->transition(this);
    parasitic = parasitics_->findPiElmore(drvr_pin, drvr_rf, parasitic_ap);
    if (parasitic)
      writeStagePiElmore(drvr_pin, parasitic);
    else {
      streamPrint(spice_stream_, "* No parasitics found for this net.\n");
      writeNullParasitics(drvr_pin);
    }
  }
}

void
WritePathSpice::writeStageParasiticNetwork(Pin *drvr_pin,
                                           Parasitic *parasitic,
                                           ParasiticAnalysisPt *parasitic_ap)
{
  Set<const Pin*> reachable_pins;
  int res_index = 1;
  int cap_index = 1;

  // Sort devices for consistent regression results.
  Vector<ParasiticDevice*> devices;
  ParasiticDeviceIterator *device_iter1 = parasitics_->deviceIterator(parasitic);
  while (device_iter1->hasNext()) {
    ParasiticDevice *device = device_iter1->next();
    devices.push_back(device);
  }
  delete device_iter1;

  sort(devices, ParasiticDeviceLess(parasitics_));

  for (ParasiticDevice *device : devices) {
    float resistance = parasitics_->value(device, parasitic_ap);
    if (parasitics_->isResistor(device)) {
      ParasiticNode *node1 = parasitics_->node1(device);
      ParasiticNode *node2 = parasitics_->node2(device);
      streamPrint(spice_stream_, "R%d %s %s %.3e\n",
                  res_index,
                  nodeName(node1),
                  nodeName(node2),
                  resistance);
      res_index++;

      const Pin *pin1 = parasitics_->connectionPin(node1);
      reachable_pins.insert(pin1);
      const Pin *pin2 = parasitics_->connectionPin(node2);
      reachable_pins.insert(pin2);
    }
    else if (parasitics_->isCouplingCap(device)) {
      // Ground coupling caps for now.
      ParasiticNode *node1 = 	parasitics_->node1(device);
      float cap = parasitics_->value(device, parasitic_ap);
      streamPrint(spice_stream_, "C%d %s 0 %.3e\n",
                  cap_index,
                  nodeName(node1),
                  cap);
      cap_index++;
    }
  }

  // Add resistors from drvr to load for missing parasitic connections.
  auto pin_iter = network_->connectedPinIterator(drvr_pin);
  while (pin_iter->hasNext()) {
    const Pin *pin = pin_iter->next();
    if (pin != drvr_pin
	&& network_->isLoad(pin)
	&& !network_->isHierarchical(pin)
	&& !reachable_pins.hasKey(pin)) {
      streamPrint(spice_stream_, "R%d %s %s %.3e\n",
		  res_index,
		  network_->pathName(drvr_pin),
		  network_->pathName(pin),
		  short_ckt_resistance_);
      res_index++;
    }
  }
  delete pin_iter;

  // Sort node capacitors for consistent regression results.
  Vector<ParasiticNode*> nodes;
  ParasiticNodeIterator *node_iter = parasitics_->nodeIterator(parasitic);
  while (node_iter->hasNext()) {
    ParasiticNode *node = node_iter->next();
    nodes.push_back(node);
  }

  sort(nodes, ParasiticNodeLess(parasitics_));

  for (ParasiticNode *node : nodes) {
    float cap = parasitics_->nodeGndCap(node, parasitic_ap);
    // Spice has a cow over zero value caps.
    if (cap > 0.0) {
      streamPrint(spice_stream_, "C%d %s 0 %.3e\n",
                  cap_index,
                  nodeName(node),
                  cap);
      cap_index++;
    }
  }
  delete node_iter;
}

void
WritePathSpice::writeStagePiElmore(Pin *drvr_pin,
                                   Parasitic *parasitic)
{
  float c2, rpi, c1;
  parasitics_->piModel(parasitic, c2, rpi, c1);
  const char *c1_node = "n1";
  streamPrint(spice_stream_, "RPI %s %s %.3e\n",
              network_->pathName(drvr_pin),
              c1_node,
              rpi);
  if (c2 > 0.0)
    streamPrint(spice_stream_, "C2 %s 0 %.3e\n",
                network_->pathName(drvr_pin),
                c2);
  if (c1 > 0.0)
    streamPrint(spice_stream_, "C1 %s 0 %.3e\n",
                c1_node,
                c1);
  
  int load_index = 3;
  auto pin_iter = network_->connectedPinIterator(drvr_pin);
  while (pin_iter->hasNext()) {
    const Pin *load_pin = pin_iter->next();
    if (load_pin != drvr_pin
	&& network_->isLoad(load_pin)
	&& !network_->isHierarchical(load_pin)) {
      float elmore;
      bool exists;
      parasitics_->findElmore(parasitic, load_pin, elmore, exists);
      if (exists) {
        streamPrint(spice_stream_, "E%d el%d 0 %s 0 1.0\n",
                    load_index,
                    load_index,
                    network_->pathName(drvr_pin));
        streamPrint(spice_stream_, "R%d el%d %s 1.0\n",
                    load_index,
                    load_index,
                    network_->pathName(load_pin));
        streamPrint(spice_stream_, "C%d %s 0 %.3e\n",
                    load_index,
                    network_->pathName(load_pin),
                    elmore);
      }
      else
        // Add resistor from drvr to load for missing elmore.
        streamPrint(spice_stream_, "R%d %s %s %.3e\n",
                    load_index,
                    network_->pathName(drvr_pin),
                    network_->pathName(load_pin),
                    short_ckt_resistance_);
      load_index++;
    }
  }
  delete pin_iter;
}

void
WritePathSpice::writeNullParasitics(Pin *drvr_pin)
{
  int res_index = 1;
  // Add resistors from drvr to load for missing parasitic connections.
  auto pin_iter = network_->connectedPinIterator(drvr_pin);
  while (pin_iter->hasNext()) {
    const Pin *load_pin = pin_iter->next();
    if (load_pin != drvr_pin
	&& network_->isLoad(load_pin)
	&& !network_->isHierarchical(load_pin)) {
      streamPrint(spice_stream_, "R%d %s %s %.3e\n",
                  res_index,
                  network_->pathName(drvr_pin),
                  network_->pathName(load_pin),
                  short_ckt_resistance_);
      res_index++;
    }
  }
  delete pin_iter;
}

void
WritePathSpice::initNodeMap(const char *net_name)
{
  stringDelete(net_name_);
  node_map_.clear();
  next_node_index_ = 1;
  net_name_ = stringCopy(net_name);
}

const char *
WritePathSpice::nodeName(ParasiticNode *node)
{
  const Pin *pin = parasitics_->connectionPin(node);
  if (pin)
    return parasitics_->name(node);
  else {
    int node_index;
    bool node_index_exists;
    node_map_.findKey(node, node_index, node_index_exists);
    if (!node_index_exists) {
      node_index = next_node_index_++;
      node_map_[node] = node_index;
    }
    return stringPrintTmp("%s/%d", net_name_, node_index);
  }
}

////////////////////////////////////////////////////////////////

// Copy the subckt definition from lib_subckt_filename for
// each cell in path to path_subckt_filename.
void
WritePathSpice::writeSubckts()
{
  StdStringSet path_cell_names = findPathCellnames();
  findPathCellSubckts(path_cell_names);

  ifstream lib_subckts_stream(lib_subckt_filename_);
  if (lib_subckts_stream.is_open()) {
    ofstream subckts_stream(subckt_filename_);
    if (subckts_stream.is_open()) {
      string line;
      while (getline(lib_subckts_stream, line)) {
	// .subckt <cell_name> [args..]
	StringVector tokens;
	split(line, " \t", tokens);
	if (tokens.size() >= 2
	    && stringEqual(tokens[0].c_str(), ".subckt")) {
	  const char *cell_name = tokens[1].c_str();
	  if (path_cell_names.find(cell_name) != path_cell_names.end()) {
	    subckts_stream << line << "\n";
	    bool found_ends = false;
	    while (getline(lib_subckts_stream, line)) {
	      subckts_stream << line << "\n";
	      if (stringBeginEqual(line.c_str(), ".ends")) {
		subckts_stream << "\n";
		found_ends = true;
		break;
	      }
	    }
	    if (!found_ends)
	      throw SubcktEndsMissing(cell_name, lib_subckt_filename_);
	    path_cell_names.erase(cell_name);
	  }
	  recordSpicePortNames(cell_name, tokens);
	}
      }
      subckts_stream.close();
      lib_subckts_stream.close();

      if (!path_cell_names.empty()) {
	string missing_cells;
        for (const string &cell_name : path_cell_names) {
	  missing_cells += "\n";
	  missing_cells += cell_name;
        }
	report_->error(28, "The subkct file %s is missing definitions for %s",
		       lib_subckt_filename_,
                       missing_cells.c_str());
      }
    }
    else {
      lib_subckts_stream.close();
      throw FileNotWritable(subckt_filename_);
    }
  }
  else
    throw FileNotReadable(lib_subckt_filename_);
}

StdStringSet
WritePathSpice::findPathCellnames()
{
  StdStringSet path_cell_names;
  for (Stage stage = stageFirst(); stage <= stageLast(); stage++) {
    TimingArc *arc = stageGateArc(stage);
    if (arc) {
      LibertyCell *cell = arc->set()->libertyCell();
      if (cell) {
	debugPrint(debug_, "write_spice", 2, "cell %s", cell->name());
	path_cell_names.insert(cell->name());
      }
      // Include side receivers.
      Pin *drvr_pin = stageDrvrPin(stage);
      auto pin_iter = network_->connectedPinIterator(drvr_pin);
      while (pin_iter->hasNext()) {
	const Pin *pin = pin_iter->next();
	LibertyPort *port = network_->libertyPort(pin);
	if (port) {
	  LibertyCell *cell = port->libertyCell();
	  path_cell_names.insert(cell->name());
	}
      }
      delete pin_iter;
    }
  }
  return path_cell_names;
}

// Subckts can call subckts (asap7).
void
WritePathSpice::findPathCellSubckts(StdStringSet &path_cell_names)
{
  ifstream lib_subckts_stream(lib_subckt_filename_);
  if (lib_subckts_stream.is_open()) {
    string line;
    while (getline(lib_subckts_stream, line)) {
      // .subckt <cell_name> [args..]
      StringVector tokens;
      split(line, " \t", tokens);
      if (tokens.size() >= 2
          && stringEqual(tokens[0].c_str(), ".subckt")) {
        const char *cell_name = tokens[1].c_str();
        if (path_cell_names.find(cell_name) != path_cell_names.end()) {
          // Scan the subckt definition for subckt calls.
          string stmt;
          while (getline(lib_subckts_stream, line)) {
            if (line[0] == '+')
              stmt += line.substr(1);
            else {
              // Process previous statement.
              if (tolower(stmt[0]) == 'x') {
                split(stmt, " \t", tokens);
                string &subckt_cell = tokens[tokens.size() - 1];
                path_cell_names.insert(subckt_cell);
              }
              stmt = line;
            }
            if (stringBeginEqual(line.c_str(), ".ends"))
              break;
          }
        }
      }
    }
  }
  else
    throw FileNotReadable(lib_subckt_filename_);
}

void
WritePathSpice::recordSpicePortNames(const char *cell_name,
				     StringVector &tokens)
{
  LibertyCell *cell = network_->findLibertyCell(cell_name);
  if (cell) {
    StringVector *spice_port_names = new StringVector;
    for (size_t i = 2; i < tokens.size(); i++) {
      const char *port_name = tokens[i].c_str();
      LibertyPort *port = cell->findLibertyPort(port_name);
      LibertyPgPort *pg_port = cell->findPgPort(port_name);
      if (port == nullptr
	  && pg_port == nullptr
	  && !stringEqual(port_name, power_name_)
	  && !stringEqual(port_name, gnd_name_))
	report_->error(29, "subckt %s port %s has no corresponding liberty port, pg_port and is not power or ground.",
		       cell_name, port_name);
      spice_port_names->push_back(port_name);
    }
    cell_spice_port_names_[cell_name] = spice_port_names;
  }
}

////////////////////////////////////////////////////////////////

Stage
WritePathSpice::stageFirst()
{
  return 1;
}

Stage
WritePathSpice::stageLast()
{
  return (path_expanded_.size() + 1) / 2;
}

string
WritePathSpice::stageName(Stage stage)
{
  string name;
  stringPrint(name, "stage%d", stage);
  return name;
}

int
WritePathSpice::stageGateInputPathIndex(Stage stage)
{
  return stage * 2 - 3;
}

int
WritePathSpice::stageDrvrPathIndex(Stage stage)
{
  return stage * 2 - 2;
}

int
WritePathSpice::stageLoadPathIndex(Stage stage)
{
  return stage * 2 - 1;
}

PathRef *
WritePathSpice::stageGateInputPath(Stage stage)
{
  int path_index = stageGateInputPathIndex(stage);
  return path_expanded_.path(path_index);
}

PathRef *
WritePathSpice::stageDrvrPath(Stage stage)
{
  int path_index = stageDrvrPathIndex(stage);
  return path_expanded_.path(path_index);
}

PathRef *
WritePathSpice::stageLoadPath(Stage stage)
{
  int path_index = stageLoadPathIndex(stage);
  return path_expanded_.path(path_index);
}

TimingArc *
WritePathSpice::stageGateArc(Stage stage)
{
  int path_index = stageDrvrPathIndex(stage);
  if (path_index >= 0)
    return path_expanded_.prevArc(path_index);
  else
    return nullptr;
}

TimingArc *
WritePathSpice::stageWireArc(Stage stage)
{
  int path_index = stageLoadPathIndex(stage);
  return path_expanded_.prevArc(path_index);
}

Edge *
WritePathSpice::stageGateEdge(Stage stage)
{
  PathRef *path = stageDrvrPath(stage);
  TimingArc *arc = stageGateArc(stage);
  return path->prevEdge(arc, this);
}

Edge *
WritePathSpice::stageWireEdge(Stage stage)
{
  PathRef *path = stageLoadPath(stage);
  TimingArc *arc = stageWireArc(stage);
  return path->prevEdge(arc, this);
}

Pin *
WritePathSpice::stageGateInputPin(Stage stage)
{
  PathRef *path = stageGateInputPath(stage);
  return path->pin(this);
}

LibertyPort *
WritePathSpice::stageGateInputPort(Stage stage)
{
  Pin *pin = stageGateInputPin(stage);
  return network_->libertyPort(pin);
}

Pin *
WritePathSpice::stageDrvrPin(Stage stage)
{
  PathRef *path = stageDrvrPath(stage);
  return path->pin(this);
}

LibertyPort *
WritePathSpice::stageDrvrPort(Stage stage)
{
  Pin *pin = stageDrvrPin(stage);
  return network_->libertyPort(pin);
}

Pin *
WritePathSpice::stageLoadPin(Stage stage)
{
  PathRef *path = stageLoadPath(stage);
  return path->pin(this);
}

const char *
WritePathSpice::stageGateInputPinName(Stage stage)
{
  Pin *pin = stageGateInputPin(stage);
  return network_->pathName(pin);
}

const char *
WritePathSpice::stageDrvrPinName(Stage stage)
{
  Pin *pin = stageDrvrPin(stage);
  return network_->pathName(pin);
}

const char *
WritePathSpice::stageLoadPinName(Stage stage)
{
  Pin *pin = stageLoadPin(stage);
  return network_->pathName(pin);
}

StdStringSet
WritePathSpice::stageOffPathPinNames(Stage stage)
{
  StdStringSet pin_names;
  if (off_path_pin_names_) {
    const PathRef *path = stageDrvrPath(stage);
    Vertex *drvr = path->vertex(this);
    VertexOutEdgeIterator edge_iter(drvr, graph_);
    while (edge_iter.hasNext()) {
      Edge *edge = edge_iter.next();
      Vertex *load = edge->to(graph_);
      const Pin *load_pin = load->pin();
      string load_pin_name = network_->pathName(load_pin);
      if (off_path_pin_names_->find(load_pin_name) != off_path_pin_names_->end())
        pin_names.insert(load_pin_name);
    }
  }
  return pin_names;
}

Instance *
WritePathSpice::stageInstance(Stage stage)
{
  Pin *pin = stageDrvrPin(stage);
  return network_->instance(pin);
}

LibertyCell *
WritePathSpice::stageLibertyCell(Stage stage)
{
  Pin *pin = stageDrvrPin(stage);
  return network_->libertyPort(pin)->libertyCell();
}

////////////////////////////////////////////////////////////////

// fprintf for c++ streams.
// Yes, I hate formatted output to ostream THAT much.
void
streamPrint(ofstream &stream,
	    const char *fmt,
	    ...)
{
  va_list args;
  va_start(args, fmt);
  char *result = nullptr;
  if (vasprintf(&result, fmt, args) == -1)
    criticalError(267, "out of memory");
  stream << result;
  free(result);
  va_end(args);
}

} // namespace
