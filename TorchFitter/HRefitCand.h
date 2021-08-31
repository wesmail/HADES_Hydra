#ifndef HRefitCand_H_
#define HRefitCand_H_

// ROOT includes
#include "TLorentzVector.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TVector3.h"

class HRefitCand : public TLorentzVector {
  private:
    Float_t fR, fZ;
    TVector3 fBase, fDir, fGP;
    TMatrixD fCov;
    Int_t fPid, fNr;

  public:
    HRefitCand();
    void setGMomentum(Float_t px, Float_t py, Float_t pz);
    void setPid(Int_t pid) { fPid = pid; }
    void setEventNr(Int_t number) { fNr = number; }
    void setCovariance(TMatrixD cov);
    void setR(Float_t val) { fR = val; }
    void setZ(Float_t val) { fZ = val; }
    void setMomentum(Float_t val);
    Int_t getPid() const { return fPid; }
    Int_t getEventNr() const { return fNr; }
    TVector3 getGMomentum() const { return fGP; }
    Float_t getR() const { return fR; }
    Float_t getZ() const { return fZ; }
    TMatrixD getCovariance() const { return fCov; }
    ClassDef(HRefitCand, 1);
};
#endif /* HRefitCand_H_ */
