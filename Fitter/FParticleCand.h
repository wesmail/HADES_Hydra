// ROOT includes
#include "TMath.h"
#include "TMatrixD.h"
#include "TLorentzVector.h"


class FParticleCand : public TLorentzVector{
    private:
        double fR, fZ;
        TMatrixD fCov;
    public:
        void setR(double val){ fR=val;}
        void setZ(double val){ fZ=val;}
        void setCovariance(TMatrixD cov);
        double getR() const { return fR;}
        double getZ() const { return fZ;}
        TMatrixD getCovariance() const {return fCov;}
};

void FParticleCand::setCovariance(TMatrixD cov){
    // 0 = 1/p
    // 1 = theta
    // 2 = phi
    // 3 = R
    // 4 = z
    fCov.ResizeTo(5,5);
    fCov = cov;
}
