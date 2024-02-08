import matplotlib.pyplot as plt
import numpy as np
import os
import time
import sys
import configparser


#dir_path = "data1/"
dir_path = sys.argv[1]
config = configparser.ConfigParser()
config.read('input.ini')

ele = [file for file in os.listdir(dir_path) if file.startswith("e")]
ion = [file for file in os.listdir(dir_path) if file.startswith("i")]

negion = [file for file in os.listdir(dir_path) if file.startswith("n")]
beam = [file for file in os.listdir(dir_path) if file.startswith("b")]


ele = sorted(ele, key=lambda x: int(x[1:-4])) # excluding first letter and txt extension
ion = sorted(ion, key=lambda x: int(x[1:-4]))
negion = sorted(negion, key=lambda x: int(x[1:-4]))
beam = sorted(beam, key=lambda x: int(x[1:-4]))

k_b=1.38e-23
m=9.1e-31
T_i=2*11604.52
v_th=np.sqrt(k_b*T_i/m)

t = config.getint('diagnostics', 'write_interval_phase') 

fig, ax = plt.subplots()

#min_length = min(len(ele), len(ion))
#ele = ele[:min_length]
#ion = ion[:min_length]

for i in range(len(ele)):
    with open(os.path.join(dir_path, ele[i]), 'r') as f:
        data_ele = f.readlines()
    with open(os.path.join(dir_path, ion[i]), 'r') as f:
        data_ion = f.readlines()
    with open(os.path.join(dir_path, negion[i]), 'r') as f:
        data_negion = f.readlines()
    with open(os.path.join(dir_path,  beam[i]), 'r') as f:
       data_beam = f.readlines()




    x_ele = [float(line.split()[0]) for line in data_ele]
    v_ele = [float(line.split()[1]) for line in data_ele]

    x_ion = [float(line.split()[0]) for line in data_ion]
    v_ion = [float(line.split()[1]) for line in data_ion]

    x_negion = [float(line.split()[0]) for line in data_negion]
    v_negion = [float(line.split()[1]) for line in data_negion]

    x_beam = [float(line.split()[0]) for line in data_beam]
    v_beam = [float(line.split()[1]) for line in data_beam]

 
    ax.clear()

    ax.scatter(x_ele, v_ele, s=1, label=f"electron - Timestep: {i*t}") 
    #ax.scatter(x_ion, v_ion, s= 1, label=f"ion - Timestep: {i*t}")
    #ax.scatter(x_negion, v_negion, s=1, label=f"negion - Timestep: {i*t}")
    #ax.scatter(x_beam, v_beam, s=1, label=f"beam - Timestep: {i*t}")


    #ax.hist(v_ele, bins=200,label= "velocity distribution")
    #ax.set_ylim([0,2500])

    #ax.hist(v_ion,bins=200)

    #ax.hist(v_negion,bins=500)
    #ax.set_ylim([0,20000])

    #ax.set_xlabel("Velocity")
    #ax.set_ylabel("")


    ax.legend()
  

   
    #ax.set_xlim([0, 2000])
    #ax.set_ylim([-300, 300.0])

    #ax.set_title("Phase Space Plot")
    #ax.set_xlabel("Position")
    #ax.set_ylabel("Velocity/Momentum")
    #ax.legend()

    plt.pause(1e-9)


plt.show()
