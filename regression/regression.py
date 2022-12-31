############################################################################
### Regression Scripts, 2022/10/14
### All the files (hypergraph file, fixed vertex file and solutio file)
### follows the hMETIS definition. Please refer to hMETIS manual for details
############################################################################
import os
from os.path import exists
import json
import time
import argparse
import shutil

#########################################################################
# Golden Evaluator
def hMetisEvaluator(hypergraph_file, solution_file, Nparts, UBfactor):
    if(os.path.exists(solution_file)):
        pass
    else:
        return 1e9, None

    with open(hypergraph_file) as f:
        content = f.read().splitlines()
    f.close()

    if(len(content) == 0):
        return 1e9, None

    items = content[0].split()
    num_hyperedges = int(items[0])
    num_vertices = int(items[1])
    flag = 0
    if(len(items) == 3):
        flag = int(items[2])

    hyperedges = []
    vertices_weight = [1 for i in range(num_vertices)]
    hyperedges_weight = [1 for i in range(num_hyperedges)]

    if((flag % 10) == 1):
        for i in range(1, num_hyperedges + 1):
            items = [int(item) for item in content[i].split()]
            hyperedges_weight[i - 1] = items[0]
            hyperedge = [items[j] - 1 for j in range(1, len(items))]
            hyperedges.append(hyperedge)
    else:
        for i in range(1, num_hyperedges + 1):
            items = [int(item) - 1 for item in content[i].split()]
            hyperedges.append(items)

    if(flag >= 10):
        for i in range(num_vertices):
            idx = i + num_hyperedges + 1
            vertices_weight[i] = int(float(content[idx]))

    part = [-1 for i in range(num_vertices)]
    with open(solution_file) as f:
        content = f.read().splitlines()
    f.close()

    if(len(content) != num_vertices):
        return 1e9, None

    blocks_balance = [0.0 for i in range(Nparts)]
    total_weight = 0.0
    for i in range(len(content)):
        part_id = int(content[i])
        if(part_id == -1):
            return 1e9, None
        part[i] = part_id
        blocks_balance[part_id] += vertices_weight[i]
        total_weight += vertices_weight[i]

    max_balance = (100.0 / Nparts  + UBfactor) * 0.01;
    for i in range(len(blocks_balance)):
        blocks_balance[i] = blocks_balance[i] / total_weight

    flag = True
    for block_balance in blocks_balance:
        if (round(block_balance, 2) > round(max_balance, 2)):
            flag = False

    num_cut = 0
    for i in range(num_hyperedges):
        part_list = []
        for vertex in hyperedges[i]:
            if(part[vertex] not in part_list):
                part_list.append(part[vertex])

        if(len(part_list) > 1):
            num_cut += hyperedges_weight[i]

    return num_cut, blocks_balance


#########################################################################
# Run Partitioner (hmetis, khmetis, tritonpart)
def RunPartitioner(partitioner, exe, hypergraph_file, Nparts, UBfactor, seed):
    if (partitioner == "hmetis"):
        cmd = exe + " " + hypergraph_file + " " + str(Nparts) + " " + str(UBfactor) + "  "
        cmd += "10 1 1 1 0 24"
        os.system(cmd)
    elif (partitioner == "khmetis"):
        cmd = exe + " " + hypergraph_file + " " + str(Nparts) + " " + str(UBfactor) + "  "
        cmd += "10 1 1 2 24"
        os.system(cmd)
    elif (partitioner == "tritonpart"):
        cmd = "triton_part_hypergraph -hypergraph_file " + hypergraph_file + " "
        cmd += "-num_parts " + str(Nparts) + " "
        cmd += "-balance_constraint " + str(UBfactor) + " "
        cmd += "-seed " + str(seed)
        cmd += "\n"
        cmd += "exit\n"
        test_file = "run_openroad.tcl"
        f = open(test_file, "w")
        f.write(cmd)
        f.close()
        cmd = exe + " " + test_file
        os.system(cmd)
        cmd = "rm " + test_file
        os.system(cmd)
    else:
        print("[INFO] Error !  No such partitioner : ", partitioner)
        exit()

#########################################################################
### Main script starts !!!
#########################################################################
### Please specify user parameters here !
print("Running Regression Sweep!")

# Setup benckmark list
design_list = ["bitonic_mesh",
                "cholesky_mc",
                "dart",
                "des90",
                "neuron",
                "openCV",
                "segmentation",
                "SLAM_spheric",
                "sparcT1_core",
                "stereo_vision"]

# Setup sweep parameters
UBfactor_list = [5, 10, 15, 20]
Nparts_list = [2, 3, 4, 5, 6, 7]
Nruns = 10
seed_list = list(range(Nruns))

##########################################################################
### Do not touch the scripts from here
##########################################################################
# Setup exe files
tritonpart_exe = "../build/src/openroad"
khmetis_exe = "../hMETIS/khmetis"
hmetis_exe = "../hMETIS/hmetis"
# Setup benchmark directory
benchmark_dir = "../titan23_benchmark/"

# Setup work directory
cur_path = os.getcwd()
tritonpart_path = cur_path + "/tritonpart_rpt/"
if os.path.isdir(tritonpart_path):
    shutil.rmtree(tritonpart_path)
    os.mkdir(tritonpart_path)
else:
    os.mkdir(tritonpart_path)

khmetis_path = cur_path + "/khmetis_rpt/"
if os.path.isdir(khmetis_path):
    shutil.rmtree(khmetis_path)
    os.mkdir(khmetis_path)
else:
    os.mkdir(khmetis_path)

hmetis_path = cur_path + "/hmetis_rpt/"
if os.path.isdir(hmetis_path):
    shutil.rmtree(hmetis_path)
    os.mkdir(hmetis_path)
else:
    os.mkdir(hmetis_path)


partitioner_list = ["hmetis", "khmetis", "tritonpart"]
partitioner_exe = {"hmetis" : hmetis_exe, "khmetis" : khmetis_exe, "tritonpart" : tritonpart_exe}
work_dir = {"hmetis" : hmetis_path, "khmetis" : khmetis_path, "tritonpart" : tritonpart_path}

summary_report = "summary.csv"
f = open(summary_report, "w")
f.close()
line = "benchmark,UBfactor,Nparts,seed,"
for partitioner in partitioner_list:
    line += partitioner + "_runtime (s),"
    line += partitioner + "_cutsize,"
    line += partitioner + "_balance,"
line += "best_cutsize,best_partitioner"
f = open(summary_report, "w")
f.write(line + "\n")
f.close()

for Nparts in Nparts_list:
    for UBfactor in UBfactor_list:
        for design in design_list:
            print("[INFO] Running ", design, " with UBfactor ", UBfactor, " and Nparts ", Nparts)
            for seed in seed_list:
                line = str(design) + "," + str(UBfactor) + "," + str(Nparts) + ","
                line += str(seed) + ","
                cutsize_list = []
                for partitioner in partitioner_list:
                    print("[INFO] Running Partitioner : ", partitioner)
                    original_hypergraph_file = benchmark_dir + design + ".hgr"
                    hypergraph_file = work_dir[partitioner] + design + ".hgr"
                    cmd = "cp " + original_hypergraph_file + " " + hypergraph_file
                    os.system(cmd)
                    start_time = time.time()
                    RunPartitioner(partitioner, partitioner_exe[partitioner],
                                   hypergraph_file, Nparts, UBfactor, seed)
                    end_time = time.time()
                    runtime = round(end_time - start_time, 3)
                    original_solution_file = hypergraph_file + ".part." + str(Nparts)
                    solution_file = hypergraph_file +  ".k." + str(Nparts)
                    solution_file += ".UBfactor." + str(UBfactor)
                    solution_file += ".seed." + str(seed)
                    cmd = "mv  " + original_solution_file + " " + solution_file
                    os.system(cmd)
                    cutsize, balance = hMetisEvaluator(hypergraph_file, solution_file,
                                                       Nparts, UBfactor)
                    print("cutsize : ", cutsize)
                    print("balance : ", balance)
                    cutsize_list.append(cutsize)
                    line += str(runtime) + "," + str(cutsize) + ","
                    balance_str = "["
                    for value in balance:
                        balance_str += " " + str(value) + " "
                    balance_str += "]"
                    line += balance_str + ","
                min_cutsize = min(cutsize_list)
                min_idx = cutsize_list.index(min_cutsize)
                line += str(min_cutsize) + "," + partitioner_list[min_idx]
                f = open(summary_report, "a")
                f.write(line + "\n")
                f.close()

# Write the improvement compared to hMETIS and khMETIS
# win_ratio, average_improve
# fail_ratio, average_improve
# runtime comparison
with open(summary_report) as f:
    content = f.read().splitlines()
f.close()

hmetis_win = 0
hmetis_fail = 0
hmetis_improve = 0
hmetis_loss = 0
hmetis_runtime = 0

khmetis_win = 0
khmetis_fail = 0
khmetis_improve = 0
khmetis_loss = 0
khmetis_runtime = 0

for i in range(1, len(content)):
    items = content[i].split(",")
    hmetis_runtime = float(items[4])
    hmetis_cutsize = float(items[5])
    khmetis_runtime = float(items[7])
    khmetis_cutsize = float(items[8])
    tritonpart_runtime = float(items[10])
    tritonpart_cutsize = float(items[11])

    # check hmetis
    if (hmetis_cutsize < tritonpart_cutsize):
        hmetis_win += 1
        hmetis_improve += (tritonpart_cutsize - hmetis_cutsize) / hmetis_cutsize
    elif (hmetis_cutsize > tritonpart_cutsize):
        hmetis_fail += 1
        hmetis_loss += (hmetis_cutsize - tritonpart_cutsize) / hmetis_cutsize
    else:
        pass

    hmetis_runtime += tritonpart_runtime / hmetis_runtime

    # check khmetis
    if (khmetis_cutsize < tritonpart_cutsize):
        khmetis_win += 1
        khmetis_improve += (tritonpart_cutsize - khmetis_cutsize) / khmetis_cutsize
    elif (khmetis_cutsize > tritonpart_cutsize):
        khmetis_fail += 1
        khmetis_loss += (khmetis_cutsize - tritonpart_cutsize) / khmetis_cutsize
    else:
        pass

    khmetis_runtime += tritonpart_runtime / khmetis_runtime


num_runs = len(content) - 1
print("************************************************")
print("************* Result Summary *******************")
print('************************************************')
print("\n")
print("------------------------------------------------")
print("--------- hMETIS VS TritonPart -----------------")
print("Fail ratio of TritonPart : ", hmetis_win / num_runs, "% (", hmetis_win, "/", num_runs, ")")
print("Average QoR degradation of TritonPart : ", hmetis_improve / hmetis_win)
print("Win ratio of TritonPart : ", hmetis_fail / num_runs, "% (", hmetis_fail, "/", num_runs, ")")
print("Average QoR improvement of TritonPart : ", hmetis_loss / hmetis_fail)
print("Runtime of TritonPart : ", hmetis_runtime / num_runs)
print("\n\n")
print("------------------------------------------------")
print("--------- khMETIS VS TritonPart -----------------")
print("Fail ratio of TritonPart : ", khmetis_win / num_runs, "% (", khmetis_win, "/", num_runs, ")")
print("Average QoR degradation of TritonPart : ", khmetis_improve / khmetis_win)
print("Win ratio of TritonPart : ", khmetis_fail / num_runs, "% (", khmetis_fail, "/", num_runs, ")")
print("Average QoR improvement of TritonPart : ", khmetis_loss / khmetis_fail)
print("Runtime of TritonPart : ", khmetis_runtime / num_runs)



