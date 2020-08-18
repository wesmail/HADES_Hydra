# Kinematic fitting tool for hydra

The tool is installed as custom module/library to the hydra.

## Requirements

1. Installed ROOT
2. Installed hydra2
3. exported `MYHADDIR`

## Installation

1. Clone the repository
```sh
git clone https://github.com/wesmail/HADES_Hydra
```

2. Enter the directory Fitter of the repository root directory
```sh
cd HADES_Hydra/Fitter
```

3. Compile and install the library
```sh
make
make install
```

4. Update the `rootlogon.C` script and add around line 88, after the last `common_libs +=` expression an additional:
```c++
    common_libs += "KineRefit";
```

The `rootlogon.C` script is either in the working directory, hardcoded in ~/.rootrc or exported as ROOTLOGON variable, depend on your system.


## Using

### The `HRefitCand` class

The `HRefitCand` object is the working object of the refitter. It derivates from `TLorentzVector` and is associated with `HVirtualCand` object. The only constructor takes `HVirtualCand` as an argument and the object is initialzied from the `HVirtualCand`.

The `HRefitCand` object has two members:

1. `reset()` - restores state (E, P, M) from the associated `HVirtualCand`
2. `update()` - updateds associated `HVirtualCand` with the refited values from the object

The `HFitter` constructor accepts single argument wihich is a vector of `HRefitCand` objects. After the `fit()` is executed, each object of `HRefitCand` is automaticaly updated with the new values of parameters. In order to copy the new values into original object, one has to call `updateDaughters()` which internally calls `HRefitCand::update()` on each object.

### I give example for Lambda decay (mass constraint)

1. Include the hitter header:
```c++
#include "hfitter.h"
```

2. Create `FParticleCand` as follows:
```c++
// set mass hypothesis for HParticleCand (e.g. called fcand1 and fcand2)
fcand1->calc4vectorProperties(938.272);
fcand2->calc4vectorProperties(139.570);

HRefitCand cand1(fcand1), cand2(fcand2);

/* define the covaraince here called e.g. cov
   cov should be 5x5
 */
TMatrixD cov1(5,5);
TMatrixD cov2(5,5);

cand1.setCovariance(cov1);
cand2.setCovariance(cov2);

std::vector<HRefitCand> list;
list.push_back(cand1);
list.push_back(cand2);
```

3. Instantiate `HFitter` object and fit:
```c++
HFitter fitter(list);
// add lambda mass (please check units I use here MeV)
fitter.addMassConstraint(1115.68); 
bool ok = fitter.fit();
// get fitted objects
HRefitCand fittedcand1 = fitter.getDaughter(0); // proton
HRefitCand fittedcand2 = fitter.getDaughter(1); // pion
// get more information
double chi2 = fitter.getChi2();
double prob = fitter.getProb();
// pull index can have values from 0 to 9 
// 0 = 1/p, 1 = theta, 2 = phi, 3 = r, 4 = z
double pull = fitter.getPull(index); 
```

4. Use the fitted objects:
```c++
fitter.updateDaughters();
TLorentzVector mother = *fcand1 + *fcand2;
// or do your other stuff
```

Enjoy!
If you have questions contact me `w.esmail@fz-juelich.de`