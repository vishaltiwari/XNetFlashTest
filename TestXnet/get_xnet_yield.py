#!~/.conda/envs/myenv_vis/bin/python

import matplotlib
matplotlib.use('agg')

import pdb
import matplotlib.pyplot as plt

import pickle

def parse_extract(file_name):
  #pdb.set_trace()
  file_obj = open(file_name,"r")
  start_line = "0.000000000000000000E+00  0.000000000000000000E+00  0.000000000000000000E+00"
  flag = 0
  abundance_dict_time = {}
  all_nuclides = []
  count = 1
  total_mass = 0
  sim_time_arr = []
  for line in file_obj:
    #print(line)
    if 'time' in line:
      continue
    if start_line in line:
      # this is the first line.
      first_line = 1
      if flag == 1:
        #if count==50:
          #pdb.set_trace()
        abundance_dict_time[count] = all_nuclides
	count += 1
        all_nuclides = []
      flag = 1

    cols = [ele for ele in line.split(' ') if ele !='']
    if first_line == 1:
      total_mass = cols[1]
      sim_time = cols[0]
      sim_time_arr.append(sim_time)
      nuclies_abud = [ float(ele)/float(total_mass) for ele in cols[9:]]
      first_line = 0
    else:
      nuclies_abud = [ float(ele)/float(total_mass) for ele in cols[:]]
      #nuclies_abud = cols[:]
    
    all_nuclides= all_nuclides + nuclies_abud
  
  abundance_dict_time[count] = all_nuclides
  return abundance_dict_time , sim_time_arr

def parse_sym_file(file_name):
  #print("Opening file:" + file_name)
  file_obj = open(file_name,"r")
  sym_list = []
  #pdb.set_trace()
  for line in file_obj:
    sym = line.split(' ')[0].strip()
    print(sym)
    sym_list.append(sym)
  return sym_list

def save_yields(file_name , abundance_list , sym_list):
  out_filename = file_name.split('.')[0] + "_yields.txt"
  file_obj = open(out_filename,"w")
  for indx , nuclide_sym in enumerate(sym_list):
    print(nuclide_sym +" "+ str(abundance_list[indx]))
    file_obj.write(nuclide_sym +" "+ str(abundance_list[indx])+"\n")
  #pdb.set_trace()

def plot_abund(abundance_dict_time , sym_list,file_name,sim_time_arr):
  #pdb.set_trace()
  plot_name = file_name.split('.')[0]+"_mass_fraction.png"
  sym_abund = {}
  thres_abund = 1e-2
  all_legends = []
  all_handles = []
  for indx, sym in enumerate(sym_list):
    sym_abund[sym] = []
    mark = False
    for iter_indx in abundance_dict_time:
      #pdb.set_trace()
      abund_val = abundance_dict_time[iter_indx][indx]
      if abund_val > thres_abund:
        mark = True
      sym_abund[sym].append(abund_val)
    # plot
    #if sym == 'fe53':
    #  pdb.set_trace()
    # NOTE:This is to plot agaist the timestep
    #handle, = plt.plot(range(len(abundance_dict_time.keys())) , sym_abund[sym], label=sym)
    
    # NOTE: this is to plot agaist simulation time
    #handle, = plt.semilogx(sim_time_arr , sym_abund[sym], label=sym)
    handle, = plt.semilogx(sim_time_arr , sym_abund[sym], label=sym)
    if mark == True:
      all_handles.append(handle)
      all_legends.append(sym)

  plt.xlabel('Time step')
  plt.ylabel('Mass fraction')
  #pdb.set_trace()
  #plt.legend(all_handles,all_legends,bbox_to_anchor=(1.04,1), ncol=3)
  plt.legend(all_handles,all_legends, bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand", borderaxespad=0, ncol=3)
  #plt.legend(all_handles,all_legends, bbox_to_anchor=(1.04,0), loc="lower left", mode="expand", borderaxespad=0, ncol=3)
  plt.savefig(plot_name,bbox_inches="tight")
  plt.clf()

def main():
  print("Inside the main function")
  #filenames = ["xnet_SN150_rho2e3.dat","xnet_SN150_rho2e7.dat", "xnet_SN160_rho2e3.dat","xnet_SN160_rho2e7.dat", "xnet_SN231_rho2e3_long.dat","xnet_SN231_rho2e7.dat","xnet_SN150_rho1e9.dat","xnet_SN160_rho1e9.dat","xnet_SN231_rho1e9.dat"]
  #nuclide_list = ["Data_SN150_list.txt","Data_SN150_list.txt","Data_SN160_list.txt","Data_SN160_list.txt","Data_SN231_list.txt","Data_SN231_list.txt","Data_SN150_list.txt","Data_SN160_list.txt","Data_SN231_list.txt"]
  #filenames = ['xnet_SN150_rho1e9_T5e9/xnetFlash_SN150_rho1e9_T5e9_yields.dat']
  #filenames = ['xnet_SN150_rho1e9_T6e9/xnetFlash_SN150.dat']
  filenames = ['xnet_alpha_rho1e9_T6e9/xnetFlash_alph.dat']
  #filenames = ['xnet_SN150_GP.dat']
  nuclide_list = ['Data_alpha_list.txt']
  pickle_data = [] 
  for indx , file_name in enumerate(filenames):
    abundance_dict_time, sim_time_arr = parse_extract(file_name)
    total_iter = len(abundance_dict_time.keys())
    sym_list = parse_sym_file(nuclide_list[indx])
    save_yields(file_name, abundance_dict_time[total_iter] , sym_list)

    # save the data into a pickel file which can be used by plot_XnetFlash_vs_XnetStand.py
    pickle_data.append(abundance_dict_time)
    pickle_data.append(sym_list)
    pickle_data.append(sim_time_arr)
    pickle.dump(pickle_data , open(file_name.split('.')[0]+"_data.p","w"))

    plot_abund(abundance_dict_time , sym_list, file_name , sim_time_arr)

if __name__ == '__main__':
  main()
