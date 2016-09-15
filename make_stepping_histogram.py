#! /usr/bin/python2.7

import subprocess
import threading




import time # testing

#### Vars ####
num_threads = 4
num_runs = 10

T = 100
Ls = 15
Lt = 14
##############

class SimulationObject:
    def __init__(self, name, seed):
        self.name = name
        self.seed = seed

def run_simulation(sim_object):
    print("Running simulation " + sim_object.name)
    cmd_list = [
        "./generate_stepping_data",
        "--T", str(T),
        "--Ls", str(Ls),
        "--Lt", str(Lt),
        "--seed", str(sim_object.seed),
        "--name", str(sim_object.name)]
    ret = subprocess.call(cmd_list)
    assert(ret == 0)

def stepping_worker_thread(run_count, run_count_mutex, simulation_objects):
    while run_count_mutex.acquire():
        if (run_count[0] == 0):
            run_count_mutex.release()
            exit()
        else:
            run_num = run_count[0]
            run_count[0] -= 1
            run_count_mutex.release()
            run_simulation(simulation_objects[run_num-1])

run_count_mutex = threading.Lock()

threads = []
simulation_objects = []
run_count = [num_runs]
for i in range(run_count[0]):
    simulation_objects.append(SimulationObject("sim-" + str(i), i))

for i in range(num_threads):
    thread = threading.Thread(target=stepping_worker_thread,
                        args=(run_count, run_count_mutex, simulation_objects))
    threads.append(thread)
    thread.start()

for t in threads:
    t.join()

data_file_txts = []
for i in range(num_runs):
    data_file_txts.append(open("data/stepping_data_" + str(simulation_objects[i].name) + ".txt", 'r').read())
