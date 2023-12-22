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

#include "TimingModel.hh"

namespace sta {

class GateLinearModel : public GateTimingModel
{
public:
  GateLinearModel(LibertyCell *cell,
                  float intrinsic,
                  float resistance);
  void gateDelay(const Pvt *pvt,
                 float in_slew,
                 float load_cap,
                 float related_out_cap,
                 bool pocv_enabled,
                 // Return values.
                 ArcDelay &gate_delay,
                 Slew &drvr_slew) const override;
  string reportGateDelay(const Pvt *pvt,
                         float in_slew,
                         float load_cap,
                         float related_out_cap,
                         bool pocv_enabled,
                         int digits) const override;
  float driveResistance(const Pvt *pvt) const override;

protected:
  void setIsScaled(bool is_scaled) override;

  float intrinsic_;
  float resistance_;
};

class CheckLinearModel : public CheckTimingModel
{
public:
  explicit CheckLinearModel(LibertyCell *cell,
                            float intrinsic);
  void checkDelay(const Pvt *pvt,
                  float from_slew,
                  float to_slew,
                  float related_out_cap,
                  bool pocv_enabled,
                  // Return values.
                  ArcDelay &margin) const override;
  string reportCheckDelay(const Pvt *pvt,
                          float from_slew,
                          const char *from_slew_annotation,
                          float to_slew,
                          float related_out_cap,
                          bool pocv_enabled,
                          int digits) const override;

protected:
  void setIsScaled(bool is_scaled) override;

  float intrinsic_;
};

} // namespace
