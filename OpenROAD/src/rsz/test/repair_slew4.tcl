# repair_design load max_slew << driver max_slew, no cap violation
source "helpers.tcl"
read_liberty Nangate45/Nangate45_typ.lib
read_liberty repair_slew4.lib
read_lef Nangate45/Nangate45.lef
read_lef repair_slew4.lef
read_def repair_slew4.def

source Nangate45/Nangate45.rc
set_wire_rc -layer metal3
estimate_parasitics -placement

report_check_types -max_slew -max_cap -max_fanout -violators
repair_design
report_check_types -max_slew -max_cap -max_fanout -violators
