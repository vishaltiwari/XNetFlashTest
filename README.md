# Notes:

## Setups in FLASH :
(Taken from the doc file that was prepared)

### Test_burn(aprox19): 
./setup test_burn -objdir=object_testburn_aprox19 -3d +cube8 -auto -maxblocks=30 -unit=physics/sourceTerms/Burn/BurnMain/nuclearBurn/Aprox19 

### Aprox13 
./setup test_burn -objdir=object_testburn_aprox13 -3d +cube8 -auto -maxblocks=30 -unit=physics/sourceTerms/Burn/BurnMain/nuclearBurn/Aprox13 

### Xnet: 
./setup test_burn -objdir=object_testburn_xnet_SN231_GPU -3d +cube8 -auto -maxblocks=30 xnet=True xnetGPU=True xnetData=Data_SN231

#### 1D, 1 cell
./setup test_burn -objdir=object_testburn_xnet_alpha_GPU -1d -auto -maxblocks=30 xnet=True xnetData=Data_alpha xnetGPU=True +noio -nxb=1  

./setup test_burn -objdir=object_testburn_xnet_SN150_GPU -1d -auto -maxblocks=30 xnet=True xnetData=Data_SN150 xnetGPU=True +noio -nxb=1 -site=master.cl.umassd.edu-xnet

For now the temperature is set by changing the temperature variable in the Burn.F90 file.
`tmp(ii,jj,kk,thisBlock)  = 5.0e9` line 226, because these is some issue in setting the temperature in Simulation_initBlock.F90.

Other Xnet parameters are set in the Config file: `source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/XNet/Config`

Density is the average of minRho and maxRho in the flash.par file. Initial abundance is explicity set to 0.5 for c12 and 016.


### Changes in the Makefile: 

### TODO: commit the site folder in the 

+CUDA_PATH = /usr/local/cuda 

+LAPACK_PATH = /home/vtiwari/local_sw/lapack-3.8.0 

+FFLAGS_CUDA  = -I${CUDA_PATH}/include 

+FFLAGS_LAPACK = -I{LAPACK_PATH}/LAPACKE/include 

+CFLAGS_OPT   = -c -mcmodel=medium -O3 -g -D_LARGEFILE64_SOURCE -D_FORTIFY_SOURCE=2 

+CFLAGS_CUDA  = -I${CUDA_PATH}/include 

+CFLAGS_LAPACK = -I{LAPACK_PATH}/LAPACKE/include 

+LFLAGS_OPT   = -mcmodel=medium -fpic -O3 –o 

+LIB_CUDA  = -L${CUDA_PATH}/lib64 -lcublas -lcudart -lcuda 

+LIB_LAPACK = -L${LAPACK_PATH} -ltmglib -llapack -lrefblas –lcblas 

### Issues: 

/autofs/nas1/home/vtiwari/flash4_fishergroup/object_testburn_xnet/Burn.F90:192:(.text+0x6c0): relocation truncated to fit: R_X86_64_PC32 against symbol 'xnet_controls_mp_nzbatchmx_' defined in COMMON section in xnet_controls.o 

### Resolved: 

+FFLAGS_OPT   = -c -mcmodel=medium -g -r8 -i4 -O3 -real_size 64 -diag-disable 10120 

### Issues: 

When not using xnetGPU flag. (xnetGPU is set to false). 

flash4: error while loading shared libraries: libcublas.so.10.0: cannot open shared object file: No such file or directory 

(Not able to find cuda?). 

 
### Issues when using xnet nuclear reaction network: 
Primary job  terminated normally, but 1 process returned 

a non-zero exit code. Per user-direction, the job has been aborted. 

------------------------------------------------------- 

-------------------------------------------------------------------------- 

mpiexec noticed that process rank 9 with PID 0 on node node50 exited on signal 9 (Killed). 

-------------------------------------------------------------------------- 

forrtl: error (78): process killed (SIGTERM) 

Image              PC                Routine            Line        Source     

libifcoremt.so.5   00007F2D2BA36345  for__signal_handl     Unknown  Unknown 

libpthread-2.17.s  00007F2D2931A6D0  Unknown               Unknown  Unknown 

libintlc.so.5      00007F2D2956EE82  intel_avx_rep_m     Unknown  Unknown 

flash4             00000000004CCDB9  Unknown               Unknown  Unknown 

flash4             000000000044DF59  Unknown               Unknown  Unknown 

flash4             000000000042B2C2  Unknown               Unknown  Unknown 

flash4             000000000041039D  Unknown               Unknown  Unknown 

flash4             00000000004178A3  Unknown               Unknown  Unknown 

flash4             0000000000405642  Unknown               Unknown  Unknown 

libc-2.17.so       00007F2D28F60445  __libc_start_main     Unknown  Unknown 

flash4             0000000000405509  Unknown               Unknown  Unknown 
 
 
### Issues on running on multiple GPUs: 

Segmentation fault, all the processes try to fit on a single GPU. Using multiple GPUS need the code to be in MPI, not sure if xnetFlash is doing that, need to check!.

### Running Standalone xnet code:
Code: https://wikihost.nscl.msu.edu/talent/doku.php?id=rxnnetcode


#### Installation:
- Unzip the directory.
- Follow instructions on the path: $XNET/doc/Compiling.txt
- cd $XNET/source
- cp Makefile_local_ifort Makefile_local (carnie has ifort)
- make
- This will create the executable xnet in the source directory.

#### Running:
- $XNET/
- Current Directory should contain the "control" file. Should be named as "control"
- The control file points to the nuclear data directory.
   - eg: $XNET/test/Data_alpha
- The control files points to the Initial Abundance and Thermodynamic Trajectory Files
  - eg: $XNET/test/Data_alpha/ab_he (Specifies the initial abund)(Ye)
     and
        $XNET/th_constant (the temperature and density)
  - Other parameter for the numerical algorithm
  - List of nuclides to output.
- Run the executable in the source dorectory, where the control file is.

#### Results:
- Results will be in the file "ev".
- The log file is the file "tso".

#### Generating REACLIB for xnet:
- use the code: https://wikihost.nscl.msu.edu/talent/doku.php?id=rxnnetcode
- Files Inside the nuclear data directory
   - sunet: this file contains the list of nuclides.(Get the list of nuclides)
   - netwinv: This is generated by the above code.
   - netsu: This contains the information of different reaction rates.(Generated from the script)
- After extracting code; gunzip reducereaclib; cd reducereaclib/
- Add the list of nuclides in the INPUT directory. Name is as sunet.
- Complie the code: ifort reducereaclib.f, and execute: ./a.out
- netwinv and netsu are generated.

#### Analysizing the binary Time-series files(using matlab):
- tools/matlab/read_ts_file.m : Parse the extract and parse the file.
- tools/matlab/plot_time_mf.m : Get Z,A, mass fraction.
- Need to write scripts that will map the Z,A to a particular nuclide, and output the mass fraction in the file format used for testing.


### Changing the parameters for the XnetFlash.
Config file: source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/XNet/Config

| XNetStandalone  | XNetFlash  | Value | 
| ---------------- |:------------------:| :------------------:|
|Initial Zone  (Number)                                     |xnet_nzbatchmx      |  1 |
|Weak Reactions (yes=1,no=0,only=-1)                        |xnet_iweak          |  1 |
|Screening (yes=1)                                          |xnet_iscrn          |  1 |
|Process Nuclear Data (yes=1,no=0)                          |xnet_iprocess       |  1 |
|integration Scheme (1=Backward Euler, 2= Bader-Deufelhard) |xnet_isolv          |  1 |
|Max. number of timesteps (number)                          |xnet_kstmx          |  6000 |
|Max. iterations per step | xnet_kitmx | 5 |
| Convergence Condition (Mass Cons.=0, (dY/Y small)=1) | xnet_iconvc | 1 |
| Max. Abundance Change per timestep | xnet_changemx | 1.0e-1 |
| Smallest Abundance used in timestep calculation | xnet_yacc | 1.0e-7 | 
| Mass Conservation Limit | xnet_tolc | 1.00E-08 |
| Convergence Criterion | xnet_ymin | 1.0e-30 |
| Max. Factor to change dt in a timestep | xnet_tdel_maxmult | 2.0e+0 |


### There are clear differences in the Nuclear Network data from XNETFLash and XnetStandalone:

#### netwinv files
XNetFlash  
n14      14.000   7   7   1.0     2.86341671                                                                                               
1.00000E+00 1.00000E+00 1.00000E+00 1.00000E+00 1.00000E+00 1.00000E+00 1.00000E+00 1.00000E+00                                             
1.00000E+00 1.00000E+00 1.00000E+00 1.00000E+00 1.00000E+00 1.00000E+00 1.00000E+00 1.00000E+00                                             
1.00000E+00 1.00000E+00 1.00000E+00 1.00000E+00 1.00000E+00 1.00000E+00 1.00000E+00 1.00000E+00                                             

XNet Standalone  
n14      14.000   7   7   1.0     2.863                                                                                                     
1.00     1.00     1.00     1.00     1.00     1.00     1.00     1.00                                                                      
1.00     1.00     1.00     1.00     1.00     1.00     1.00     1.00                                                                      
1.00     1.00     1.00     1.00     1.00     1.00     1.00     1.00                                                                      





#### netsu:
The reaction with label *ffn* is missing in the XNetStandlone. 
Example:

XnetStandalone  
0     na22 ne22                            wc12w     2.84300E+0  
-0.185900E+02 0.000000E+00 0.000000E+00 0.000000E+00  
0.000000E+00 0.000000E+00 0.000000E+00  

XNetFlash  
0     na22 ne22                             ffn      2.33100E+00  
0.140000E+02 0.000000E+00 0.000000E+00 0.000000E+00  
0.000000E+00 0.000000E+00 0.000000E+00            

## Skynet on Carnie

Bitkucket: git clone https://bitbucket.org/jlippuner/skynet.git

#### Building GSL
- wget ftp://ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz
- gunzip and tar the file.
- ./configure --prefix=/home/vtiwari/lib/gsl
- make install -j 16

#### Building Boost
- Download the source from https://www.boost.org/.
- ./bootstrap.sh
- ./b2 -j 16
- ./b2 headers
- Libs in $BOOST/stage/lib
- Include files in $BOOST/boost
- <b>export BOOST_ROOT=/home/vtiwari/local_sw/boost/boost_1_70_0/</b>

#### Installing Pardiso
- Download the library from the website: https://pardiso-project.org/
- Put the lic in the file: /home/vtiwari/.local/lic/pardiso.lic

#### Building Swig:
- Download the source from :
- ./configure --prefix=/home/vtiwari/lib/swig
- make install -j 16

SWIG_DIR=
<b>export PATH=/home/vtiwari/lib/swig/bin:$PATH</b>

CMAKE PATH  
<b>export CMAKE_PREFIX_PATH=/home/vtiwari/lib/gsl/:/home/vtiwari/local_sw/lapack-3.8.0:/home/vtiwari/lib/pardiso:/home/vtiwari/local_sw/boost/boost_1_70_0:$CMAKE_PREFIX_PATH</b>

For HDF5:  
export LD_LIBRARY_PATH=/sw/lib:$LD_LIBRARY_PATH

INSTALL:
make -j 16 install

#### Issues:
/home/vtiwari/lib/pardiso/libpardiso600-GNU800-X86-64.so: undefined reference to `log2f@GLIBC_2.27`                                      
/home/vtiwari/lib/pardiso/libpardiso600-GNU800-X86-64.so: undefined reference to `logf@GLIBC_2.27`    

undefined reference to `H5::CommonFG::openDataSet`  

Not able to find SWIG, after changing the PATH variable, its not able to find libz files.


## torch nuclear reaction network.


## Comparing the error between Skynet(NSE) vs timmes Touch(NSE) vs XnetStandalone vs XNetFlash

