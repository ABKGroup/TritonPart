# Regression Tests

This folder containts scripts pertaining to running regression tests with TritonPart. The regression script is available in ```regression.py```. 

# How to run regression

We recommend users to set their input parameters:
- UBfactor_list: List of imbalance factors that the user want to experiment 
- Nparts_list: List of number of partitions that the user want to experiment 
- Nruns: The number of times each testcase will be run

We recommend users to set the partition executables by modifying line 153 to line 155 in ```regression.py```. 
We recommend users to set the benchmark directory by modifying lines 132 and 157 in ```regression.py```. 
