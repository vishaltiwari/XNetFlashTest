#!/usr/bin/python

import matplotlib
matplotlib.use('agg')

import pdb
import matplotlib.pyplot as plt

import numpy as np

plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True

def plot_diff(sym_diff_dict,skynet_yield , xnet_yield,base,target , err_plot_flag):
  if err_plot_flag == False:
    output = base + "_vs_"+target + ".png"
    fig, ax = plt.subplots()
  else:
    output = base + "_vs_"+target + "_error.png"
    ax = plt.gca()
    #fig, ax = plt.plot()
  #pdb.set_trace()
  width = 0.25
  for indx, sym in enumerate(sym_diff_dict):
    if err_plot_flag == False:
      xnet_bar, = ax.bar(indx - width/2, np.log(xnet_yield[sym]), width,color='blue',label=base)
      skynet_bar, = ax.bar(indx + width/2, np.log(skynet_yield[sym]), width ,color='red', label=target)
    else:
      err_bar, = ax.bar(indx , np.log(sym_diff_dict[sym]),width , color='orange' , label="Error")

  ax.set_xlabel('Nuclide')
  if err_plot_flag == False:
    ax.set_ylabel('log(Yield)')
  else:
    ax.set_ylabel('Log(Relative Error)')

  ax.set_xticks(xrange(len(sym_diff_dict.keys())))
  ax.set_xticklabels(sym_diff_dict.keys())
  if err_plot_flag == False:
    ax.legend([xnet_bar , skynet_bar],[base,target])
  else:
    ax.legend([err_bar],["Error"])
    
  plt.savefig(output,bbox_inches="tight")

def main():
  #skynet_results = ['skynet_SN150_yields.txt' , 'skynet_SN160_yields.txt' , 'skynet_SN231_yields.txt']
  #xnet_results = ['xnet_SN150_yields.txt' , 'xnet_SN160_yields.txt' , 'xnet_SN231_yields.txt']
  #skynet_results = [ 'skynet_SN160_T6e9_rho1e9yields.txt']
  #xnet_results = ['xnet_SN231_rho1e9_yields.txt']
  #skynet_results = [ 'xnet_standalone_alpha_yields.txt']
  #xnet_results = ['xnet_alpha_rho1e9_yields.txt']
  #skynet_results = ['/home/vtiwari/runs/XNetFlashTest/xnet_standalone/branches/public/my_code/xnetStandalone_SN150_yields.txt']
  #skynet_results = ['xnet_standalone_alpha_dtmax_1e2_yields.txt']
  err_plot_flag = False
  base_net = 'XnetStandalone'
  target_net = 'XnetFlash'
  skynet_results = ['xnet_standalone_alpha_rho1e9_T5e9_yields.txt']
  xnet_results = ['xnetFLASH_alpha_rho1e9_T5e9_yields.txt']

  for indx , file_name in enumerate(skynet_results):
    print("comparing:" + file_name + " " + xnet_results[indx])
    skynet_file_obj = open(file_name,"r")
    xnet_file_obj = open(xnet_results[indx],"r")
    skynet_abud_dic = {}
    xnet_abud_dic = {}
    for line in skynet_file_obj:
      sym = line.split(' ')[0]
      mass_frac = line.split(' ')[1]
      skynet_abud_dic[sym] = float(mass_frac)
    for line in xnet_file_obj:
      sym = line.split(' ')[0]
      mass_frac = line.split(' ')[1]
      xnet_abud_dic[sym] = float(mass_frac)

    # Actual comparision:
    eps = 1e-100
    diff_sum = 0
    sym_diff_dict = {}
    total_nuclide = len(skynet_abud_dic.keys())
    for sym in skynet_abud_dic:
      skynet_yield = skynet_abud_dic[sym]
      xnet_yield = xnet_abud_dic[sym]
      rel_diff = (skynet_yield - xnet_yield)**2/(skynet_yield+eps)
      print("sym:" + sym +" "+ str(skynet_yield) + " " +str(xnet_yield) + " "+str(rel_diff))
      sym_diff_dict[sym] = rel_diff
      diff_sum += rel_diff/total_nuclide

    plot_diff(sym_diff_dict,skynet_abud_dic , xnet_abud_dic,base_net , target_net,err_plot_flag)

    print("Relative Error:" + str(diff_sum))

if __name__ == '__main__':
  main()
