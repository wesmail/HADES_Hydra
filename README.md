## Deep Learig Based PID in HADES  
### How to run:  
⋅⋅* Open `Hades_Pid_NN.ipynb` in google colab and train and save your model.  
⋅⋅* `git clone https://github.com/wesmail/HADES_Hydra.git` and then `git checkout pid`.  
⋅⋅* `cd forHydra && mkdir build && cd build`.  
⋅⋅* Move the trained model in `build` directory. 
⋅⋅* `cmake ../ && make`.  
⋅⋅* `root -l`.  
⋅⋅* `gSystem->Load("libAnalysis.so`.  
⋅⋅* `Analysis *obj = new Analysis()`.  
⋅⋅* `l->Loop("DST.root",1000)`. 
⋅⋅* If you have any question feel free to ask `w.esmail@fz-juelich.de`. 