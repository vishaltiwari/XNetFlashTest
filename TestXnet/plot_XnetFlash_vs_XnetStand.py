#!~/.conda/envs/myenv_vis/bin/python

import matplotlib
matplotlib.use('agg')

import pdb
import matplotlib.pyplot as plt

import pickle

def plot_abund(abundance_dict_time , sym_list,file_name,sim_time_arr):
  #pdb.set_trace()
  # change the name to compare
  plot_name = file_name.split('/')[0]+"XNETFLASH_VS_STANDALONE.png"
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
  print("This is Main function")
  # plot the pickel data and plot it.
  # XNETFLASH pickle file:
  Xnetflash_filename = 'xnet_SN150_rho1e9_T5e9/xnetFlash_SN150_rho1e9_T5e9_yields_data.p'
  pdb.set_trace()
  pickle_data = pickle.load(open(Xnetflash_filename,"r"))
  abundance_dict_time = pickle_data[0]
  sym_list = pickle_data[1]
  sim_time_arr = pickle_data[2]

  plot_abund(abundance_dict_time , sym_list , Xnetflash_filename , sim_time_arr)

if __name__ == "__main__":
  main()
