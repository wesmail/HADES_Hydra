## kinematic fitting tool for hydra  
### how to run (in your analysis code)  
### I give example for Lambda decay (mass constraint)
1. `#include "HFitter.h"`
2. You have to fill the `FParticleCand` as follows:  
```
// set mass hypo for HParticleCand  (e.g. called fcand1 and fcand2)
fcand1->calc4vectorProperties(938.272);
fcand2->calc4vectorProperties(139.570);
FParticleCand cand1, cand2;
cand1.SetXYZM(fcand1->Px(),fcand1->Py(),fcand1->Pz(),fcand1->M());
cand1.setR(fcand1->getR());
cand1.setZ(fcand1->getZ());
cand2.SetXYZM(fcand2->Px(),fcand2->Py(),fcand2->Pz(),fcand2->M());
cand2.setR(fcand2->getR());
cand2.setZ(fcand2->getZ());
/* define the covaraince here called e.g. cov
  cov should be 5x5
  TMatrixD cov(5,5);
  */
cand1.setCovariance(cov1);
cand2.setCovariance(cov2);
std::vector<FParticleCand> list;
list.push_back(cand1);
list.push_back(cand2);
```
3. instantiate `HFitter` object and fit   
```
HFitter fitter(2, list);
// add lambda mass (please check units I use here MeV)
fitter.addMassConstraint(1115.68); 
bool ok = fitter.fit();
// get fitted objects
FParticleCand fittedcand1 = fitter.getDaughter(0); // proton
FParticleCand fittedcand2 = fitter.getDaughter(1); // pion
// get more information
double chi2 = fitter.getChi2();
double prob = fitter.getProb();
// pull index can have values from 0 to 9 
// 0 = 1/p, 1 = theta, 2 = phi, 3 = r, 4 = z
double pull = fitter.getPull(index); 
```

Enjoy!
If you have questions contact me `w.esmail@fz-juelich.de`
