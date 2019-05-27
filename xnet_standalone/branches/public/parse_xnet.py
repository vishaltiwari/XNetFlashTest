#!/home/vtiwari/.conda/envs/myenv_vis/bin/python

import matplotlib
matplotlib.use('agg')

import pdb
import matplotlib.pyplot as plt

def plot_yields(abund_dict):
  #pdb.set_trace()
  all_handles = []
  all_labels = []
  for sym in abund_dict.keys():
    print(sym +" "+ str(abund_dict[sym][-1]))
    all_labels.append(sym)
    sym_vals = abund_dict[sym]
    handle, = plt.plot(range(len(sym_vals)) , sym_vals, label=sym)
    all_handles.append(handle)
  plt.xlabel('Time Step')
  plt.ylabel('Mass Fraction')
  plt.legend(all_handles,all_labels, bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand", borderaxespad=0, ncol=3)
  plt.savefig('xnet_standalone_alpha.png',bbox_inches="tight")

def main():
  #print('parsing the xnet results')
  filename = 'run1/ev1'
  file_obj = open(filename,"r")
  count = 0
  abund_dict = {}
  for line in file_obj:
    #print(line)
    cols = [ele for ele in line.split(' ') if ele !='']
    if count == 0:
      #pdb.set_trace()
      headers = cols
      nuclide_list = cols[6:-2]
      for nuclide_sym in nuclide_list:
        abund_dict[nuclide_sym] = []
      count += 1
      continue
    abunds = cols[6:-2]
    #pdb.set_trace()
    count += 1

    for indx , yield_amt in enumerate(abunds):
      abund_dict[nuclide_list[indx]].append(float(yield_amt))
  #pdb.set_trace()
  #plot the abundances:
  plot_yields(abund_dict)

if __name__ == '__main__':
  main()
