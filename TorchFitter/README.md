## Kinematic Refit tool for hydra  
### libtorch implementation
You will need `cmake` to install this lib  
1. login to your virgo account.  
2. download and install cmake [cmake.org](https://cmake.org/install/)  
3. copy `libtorch-shared-with-deps-1.2.0.zip` in `lustre/hades/user/wesmail/` to your working directory.  
4. `unzip libtorch-shared-with-deps-1.2.0.zip`  
5. set the enviroment variable `TORCHLIB_PATH` to the `libtorch` location by typing `export TORCHLIB_PATH=/path/to/unzipped/libtorch`.  
6. clone the repo `git clone https://github.com/wesmail/HADES_Hydra.git`.  
7. `cd HADES_Hydra/TorchFitter`.  
8. create a build directory. `mkdir build`, then `cd build`.  
9. activate hydra. This lib is tested on `Hydra2.4.9v`.  
10. `cmake ../`
11. `make -j8`
12. load the lib in`rootlogon.C` by adding `gSystem->Load("libHFitter.so");`.  
13. enjoy!.  
