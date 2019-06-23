## Building the Xnet Standalone

Link: https://github.com/starkiller-astro/XNet

### Download Timmes helmholz data:
- Download helmholtz.tbz from http://cococubed.asu.edu/code_pages/eos.shtml
- In the helmholtz.f90 file, get rid of the main function, as the compiler will start to cry about the of mains.
- In the $XNET_DIR/source/Makefile : Change HELMHOLTZ_PATH = ../../../EOS/Helmholtz --> HELMHOLTZ_PATH = ../helmholtz
- In Makefile.internal, remove/comment out the MKL_INC , MKL_LIB part, because there is a syntax error regarding space/tabs.

- After this do a make -j 16
- This will create the executable.

### Running xnet.
- Make a new directory.
- Copy the control file, and the Data_<Network> , Data_alpha, Data_SN150, etc. 
- Point the data files in the control files.
- Process Nuclear Data at Run Time should be set to one.
