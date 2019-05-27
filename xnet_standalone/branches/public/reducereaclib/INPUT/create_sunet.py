#!/home/vtiwari/runs/XNetFlashTest/xnet_standalone/branches/public/reducereaclib/INPUT

import pdb

def main():
  print("Main")
  # Read the sunet original file and match the nuclides from 
  filename1 = "sunet.original"
  filename2 = "Data_SN231_list.txt"
  sunet_filename = filename2.split('_')[1] + ".sunet"
  file_obj1 = open(filename1,"r")
  file_obj2 = open(filename2,"r")
  write_obj = open(sunet_filename , "w")

  #pdb.set_trace()
  nuclides1 = []
  for line in file_obj1:
    nuclides1.append(line.split('\n')[0])

  nuclides2 = []
  for line in file_obj2:
    nuclides2.append(line.split('\n')[0])
  
  #pdb.set_trace()
  out_nuclide = []
  for nuclide in nuclides1:
    sym = [ele for ele in nuclide.split(' ') if ele !=''][0]
    for nuc in nuclides2:
      if nuc == sym:
        write_obj.write(nuclide+"\n")
	break
    
  write_obj.close()
  print('Done creating the sunet file')

if __name__ == "__main__":
  main()
