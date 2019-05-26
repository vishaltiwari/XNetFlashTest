#!/usr/bin/python

from SkyNet import *

import numpy as np
import matplotlib.pyplot as plt

import pdb

def parse_file(file_name):
  print("Opening file:" + file_name)
  file_obj = open(file_name,"r")
  dic_sym = {}
  for line in file_obj:
    sym = line.split(' ')[0].strip()
    print(sym)
    dic_sym[sym] = 0
  return dic_sym

def run_nse(nuclides_dic,file_name):
  #nucs = open("Data_alpha_list.txt", "r").readlines()
  nucs = open(file_name, "r").readlines()
  #nucs = open("nuclide_names", "r").readlines()
  nucs = [n.strip() for n in nucs]
  # cook everything at temp = 5e10 , rho = 2e7, ye=0.5
  temp = 5 # This is 5 billion kelvin
  rho = 2e7
  ye = 0.5

  print(SkyNetRoot)
  helm = HelmholtzEOS(SkyNetRoot + "/data/helm_table.dat")
  nuclib = NuclideLibrary.CreateFromWebnucleoXML(SkyNetRoot + "/data/webnucleo_nuc_v2.0.xml", nucs)
  
  nuc_species_index = {}
  nuc_species_As = {}
  for nuclide_sym in nuclides_dic:
    if nuclide_sym in nuclib.NuclideIdsVsNames():
      nuc_species_index[nuclide_sym] = nuclib.NuclideIdsVsNames()[nuclide_sym]
  
  nse_opts = NSEOptions()
  nse_opts.SetDoScreening(False)
  nse = NSE(nuclib, helm, None, nse_opts)

  nsc = nse.CalcFromTemperatureAndDensity(temp, rho, ye).Y()

  for nuclide_sym in nuclides_dic:
    if nuclide_sym in nuclib.NuclideIdsVsNames():
      print(nuclide_sym + " " + str(nsc[nuc_species_index[nuclide_sym]]))
      nuclides_dic[nuclide_sym] = nsc[nuc_species_index[nuclide_sym]] * nuclib.As()[nuc_species_index[nuclide_sym]]
  pdb.set_trace()
  return nuclides_dic

def write_to_file(nuclides_dic,filename):
  outfile = filename.split('.')[0] + "_out.txt"
  outfile_obj  = open(outfile,"w")
  for nuclide_sym in nuclides_dic:
    outfile_obj.write(nuclide_sym + " " + str(nuclides_dic[nuclide_sym])+"\n")
  outfile_obj.close()

def main():

  print("Calculating the nuclide yields")
  #filenames = ['Data_alpha_list.txt'  ,'Data_SN150_list.txt',  'Data_SN160_list.txt' , 'Data_SN231_list.txt']
  filenames = ['Data_SN150_list.txt',  'Data_SN160_list.txt' , 'Data_SN231_list.txt']
  for file_name in filenames:
    nuclides_dic = parse_file(file_name)
    # Run NSE for this nuclides_dic
    nuclides_dic = run_nse(nuclides_dic,file_name)
    write_to_file(nuclides_dic,file_name)

if __name__ == "__main__":
  main()
