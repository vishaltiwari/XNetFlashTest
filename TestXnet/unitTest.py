#!/usr/bin/python

def main():
  #skynet_results = ['skynet_SN150_yields.txt' , 'skynet_SN160_yields.txt' , 'skynet_SN231_yields.txt']
  #xnet_results = ['xnet_SN150_yields.txt' , 'xnet_SN160_yields.txt' , 'xnet_SN231_yields.txt']
  skynet_results = [ 'skynet_SN160_yields.txt']
  xnet_results = ['xnet_SN160_yields.txt']

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
    for sym in skynet_abud_dic:
      skynet_yield = skynet_abud_dic[sym]
      xnet_yield = xnet_abud_dic[sym]
      print("sym:" + sym +" "+ str(skynet_yield) + " " +str(xnet_yield))
      rel_diff = (skynet_yield - xnet_yield)**2/(skynet_yield+eps)
      diff_sum += rel_diff

    print("Relative Error:" + str(diff_sum))



if __name__ == '__main__':
  main()
