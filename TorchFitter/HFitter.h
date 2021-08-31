#ifndef HFitter_H_
#define HFitter_H_

#include "HRefitCand.h"
#include "TMath.h"
#include "TMatrix.h"
#include "TObject.h"

#include <cmath>
#include <iostream>
#include <vector>

class HFitter : public TObject {
  private:
    TMatrixD fA, fPull, fV;
    double fChi2, fProb;
    int fN, fDim, fSize, fNr; // event number

    // data members for constraints
    double fMass, fIMass, fMMass;
    int fNdf;
    std::vector<double> fM;
    TLorentzVector fInit;
    bool fMassConstraint, f3CConstraint, f2CConstraint, f4CConstraint;

    // additional for MC-Truth
    std::vector<int> fPid;
    std::vector<TVector3> fGMomentum;

  public:
    HFitter(int n);
    virtual ~HFitter(){};
    double getChi2() const { return fChi2; }
    double getProb() const { return fProb; }
    double getPull(int val = 0) { return fPull(val, val); }
    void addMassConstraint(double mass);
    void add2CConstraint(float inv_mass, float miss_mass, TLorentzVector init);
    void add3CConstraint(float inv_mass, float miss_mass, TLorentzVector init);
    void add4CConstraint(float inv_mass, float miss_mass, TLorentzVector init);
    bool fit(std::vector<HRefitCand> cands);
    HRefitCand getDaughter(int val);
    ClassDef(HFitter, 1);
};
#endif /* HFitter_H_ */
