#!/usr/bin/python2.7

import numpy as np
import subprocess

MT_BINDING_HEIGHT = 0.2

assert(subprocess.call("make histogram-stuff", shell=True) == 0)

k_b = 0 # s^-1
runtime = 1e6 * 1e-11

data_file_basename = "binding_time_fraction"
result_fname = 'simulations/simulation_results/binding_time_fraction.txt'

generate_cmd = [
    "./generate_stepping_data",
    "--k_b", str(k_b),
    "--name", str(data_file_basename),
    "--runtime", str(runtime),
    "--movie"
]
print ' '.join(generate_cmd)
subprocess.call(generate_cmd)

mv_cmd = ["mv", "data/stepping_movie_data_" + data_file_basename + ".txt",
          result_fname + ".data"]
print ' '.join(mv_cmd)
subprocess.call(mv_cmd)

rm_cmd = ["rm",
          "data/stepping_data_" + data_file_basename + ".txt",
          "data/stepping_config_" + data_file_basename + ".txt",
          "data/stepping_movie_config_" + data_file_basename + ".txt"]
print ' '.join(rm_cmd)
subprocess.call(rm_cmd)

movie_data = np.genfromtxt(result_fname + ".data", invalid_raise=False)
ubys = [datum[17] for datum in movie_data[1:]]

num_binding_timesteps = np.sum(np.where(ubys < MT_BINDING_HEIGHT))
num_timesteps = len(ubys)

print ubys

rm2_cmd = ["rm", result_fname + ".data"]
print ' '.join(rm2_cmd)
subprocess.call(rm2_cmd)

print "timesteps with uby < %g: %g/%g, fraction: %g" % (MT_BINDING_HEIGHT, num_binding_timesteps, num_timesteps, num_binding_timesteps / num_timesteps)

result_file = open(result_fname, 'w')
result_file.write("timesteps with uby  <%g: %g/%g, fraction: %g" % (MT_BINDING_HEIGHT, num_binding_timesteps, num_timesteps, num_binding_timesteps / num_timesteps))
