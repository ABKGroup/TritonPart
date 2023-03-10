# repair_design long wire from input port
# in1-----------u1--u2--out1
#       1500u
source "helpers.tcl"
read_liberty Nangate45/Nangate45_typ.lib
read_lef Nangate45/Nangate45.lef
read_def repair_wire7.def

source Nangate45/Nangate45.rc
set_wire_rc -layer metal3
estimate_parasitics -placement

report_checks -unconstrained -fields {input slew cap} -digits 3 -rise_to out1
report_long_wires 3

# wire length = 1500u -> 2 buffers required
repair_design -max_wire_length 600

report_long_wires 3
report_checks -unconstrained -fields {input slew cap} -digits 3 -rise_to out1
