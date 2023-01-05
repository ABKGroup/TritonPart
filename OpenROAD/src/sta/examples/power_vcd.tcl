# read_vcd_activities gcd
read_liberty sky130hd_tt.lib
read_verilog gcd_sky130hd.v
link_design gcd

read_sdc gcd_sky130hd.sdc
read_spef gcd_sky130hd.spef
# Generate vcd file
#  iverilog -o gcd_tb gcd_tb.v
#  vvp gcd_tb
read_power_activities -scope gcd_tb/gcd1 -vcd gcd_sky130hd.vcd
report_power
