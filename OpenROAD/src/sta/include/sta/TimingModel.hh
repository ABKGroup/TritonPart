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

#pragma once

#include <string>

#include "Delay.hh"
#include "LibertyClass.hh"

namespace sta {

using std::string;

// Abstract base class for GateTimingModel and CheckTimingModel.
class TimingModel
{
public:
  TimingModel(LibertyCell *cell);
  virtual ~TimingModel() {}
  virtual void setIsScaled(bool is_scaled) = 0;

protected:
  LibertyCell *cell_;
};

// Abstract base class for LinearModel and TableModel.
class GateTimingModel : public TimingModel
{
public:
  GateTimingModel(LibertyCell *cell);
  // Gate delay calculation.
  virtual void gateDelay(const Pvt *pvt,
			 float in_slew,
			 float load_cap,
			 float related_out_cap,
			 bool pocv_enabled,
			 // Return values.
			 ArcDelay &gate_delay,
			 Slew &drvr_slew) const = 0;
  virtual string reportGateDelay(const Pvt *pvt,
                                 float in_slew,
                                 float load_cap,
                                 float related_out_cap,
                                 bool pocv_enabled,
                                 int digits) const = 0;
  virtual float driveResistance(const Pvt *pvt) const = 0;
};

// Abstract base class for timing check timing models.
class CheckTimingModel : public TimingModel
{
public:
  CheckTimingModel(LibertyCell *cell);
  // Timing check margin delay calculation.
  virtual void checkDelay(const Pvt *pvt,
			  float from_slew,
			  float to_slew,
			  float related_out_cap,
			  bool pocv_enabled,
			  // Return values.
			  ArcDelay &margin) const = 0;
  virtual string reportCheckDelay(const Pvt *pvt,
                                  float from_slew,
                                  const char *from_slew_annotation,
                                  float to_slew,
                                  float related_out_cap,
                                  bool pocv_enabled,
                                  int digits) const = 0;
};

} // namespace
