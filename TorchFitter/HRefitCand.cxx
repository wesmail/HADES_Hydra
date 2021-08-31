#include "HRefitCand.h"

ClassImp(HRefitCand);
HRefitCand::HRefitCand() : TLorentzVector() {
    fCov.ResizeTo(5, 5);
    fCov.Zero();
}

void HRefitCand::setMomentum(Float_t val) {
    HRefitCand::SetXYZM(val * TMath::Sin(HRefitCand::Theta()) * TMath::Cos(HRefitCand::Phi()),
                        val * TMath::Sin(HRefitCand::Theta()) * TMath::Sin(HRefitCand::Phi()), val * TMath::Cos(HRefitCand::Theta()),
                        HRefitCand::M());
}

void HRefitCand::setGMomentum(Float_t px, Float_t py, Float_t pz) { fGP.SetXYZ(px, py, pz); }

void HRefitCand::setCovariance(TMatrixD cov) { fCov = cov; }
