# Notes:

## Setups in FLASH :
(Taken from the doc file that was prepared)

### Test_burn(aprox19): 
./setup test_burn -objdir=object_testburn_aprox19 -3d +cube32 -auto -maxblocks=30 -unit=physics/sourceTerms/Burn/BurnMain/nuclearBurn/Aprox19 

### Aprox13 
./setup test_burn -objdir=object_testburn_aprox13 -3d +cube32 -auto -maxblocks=30 -unit=physics/sourceTerms/Burn/BurnMain/nuclearBurn/Aprox13 

### Xnet: 
./setup test_burn -objdir=object_testburn_xnet_SN231_GPU -3d +cube32 -auto -maxblocks=30 xnet=True xnetGPU=True xnetData=Data_SN231 

 

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

Segmentation fault, all the processes try to fit on a single GPU. 
