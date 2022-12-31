# Regression Tests

This folder containts scripts pertaining to running regression tests with TritonPart. The regression script is available in ```regression.py```. 

# How to run regression

We break the process into four primary steps

## Setting input parameters

We recommend users to set their input parameters:
- ```UBfactor_list```: List of imbalance factors that the user want to experiment 
- ```Nparts_list```: List of number of partitions that the user want to experiment 
- ```Nruns```: The number of times each testcase will be run

## Setting benchmark information

We recommend users to set the benchmark directory by modifying lines [132](https://github.com/ABKGroup/TritonPart/blob/68e516145a7090c3b9bae7ac9bf2464e58758b69/regression/regression.py#L132) and [157](https://github.com/ABKGroup/TritonPart/blob/68e516145a7090c3b9bae7ac9bf2464e58758b69/regression/regression.py#L157) in ```regression.py```. 


## Setting partitioner executable

We recommend users to set the partition executables by modifying line [153](https://github.com/ABKGroup/TritonPart/blob/68e516145a7090c3b9bae7ac9bf2464e58758b69/regression/regression.py#L153) to line [155](https://github.com/ABKGroup/TritonPart/blob/68e516145a7090c3b9bae7ac9bf2464e58758b69/regression/regression.py#L155) in ```regression.py```. 

## Running ```regression.py```
``` shell
python regression.py | tee summary.log
```



