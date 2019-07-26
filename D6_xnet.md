### Setup:
./setup Super3D_hd -objdir=object_0610_mpollvl8_maxblock_30 -3d +cube32 -auto -maxblocks=30 +uhd +newMpole –with-unit=Particles +parallelio xnet=True xnetData=Data_SN150 xnetGPU=True -site=master.cl.umassd.edu-xnet

### Requesting GPUS and CPUS in XnetFlash:
When running the flash exe in the sbatch script: mpiexec -np <> ./flash4. 

`#SBATCH --nodes=2 --sockets-per-node=2 --cores-per-socket=12 --threads-per-core=1  #nodes / tasks distribution`  
This is leading to a segmentation fault. This configuration states that this node can "host" 24 mpi process, and the mpi code will start to give each node 24 mpi process, even though having 24 mpi process will exceed the memory for that node, which will lead to a segmentation fault.

Thus for carnie, this is the working configuration:

`#SBATCH --nodes=2 --sockets-per-node=2 --cores-per-socket=7 --threads-per-core=1  #nodes / tasks distribution`  

or  

`#SBATCH --nodes=6`  
`#SBATCH --ntasks-per-node=14`  

Also, its important to note that a total of 14 process can be fit onto a 2.4GB Tesla M2050 GPUs.
