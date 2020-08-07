/**
 * HFitter.h
 *
 * \brief Kinematic and Vertex fitting tool for Hydra
 * \author Waleed Esmail <w.esmail@fz-juelich.de>
 * \date Aug 07, 2020
 *
 */


// system includes
#include <iostream>
#include <vector>
#include <cmath>

// framework includes
#include "FParticleCand.h"

const double pi2 = TMath::PiOver2();

template<typename T>
void Print(T const &matrix)
{
    Int_t nrows = matrix.GetNrows();
    Int_t ncols = matrix.GetNcols();
    
    cout << "shape(" << nrows << "," << ncols << ")" << std::endl;

    for (Int_t i=0; i<nrows; i++){
        for (Int_t j=0; j<ncols; j++){
            Double_t element = matrix(i,j);
            if ( TMath::Abs(element) < 1e-10 ) element=0.;
            if ( element >= 0.)
              std::cout << " " << std::fixed << std::setw(8) << std::scientific << element << " ";
            else
              std::cout << std::fixed << std::setw(8) << std::scientific << element << " ";
        }
        std::cout <<endl;
    }

    std::cout << std::endl;
}

class HFitter{
    private:
        TMatrixD y, V, fPull;
	    double fChi2, fProb;
        bool   fConverged;
        int    fIteration, fN;

        // data members for constraints
        double fMass;
        int fNdf;
        std::vector<double> fM;
        TLorentzVector fInit;
        bool fMassConstraint, fMMConstraint, fMassVtxConstraint, fVtxConstraint;
    public:
      HFitter(int n, std::vector<FParticleCand> cands);
      ~HFitter(){};
      TMatrixD f_eval(const TMatrixD &m_iter);
      TMatrixD Feta_eval(const TMatrixD &m_iter);
      void   addMassConstraint(double mass);
      void   addMassVtxConstraint(double mass);
      void   addVtxConstraint();
      void   addMissingMassConstraint(double mass, TLorentzVector init);
      double getChi2() const {return fChi2;}
      double getProb() const {return fProb;}
      double getPull(int val=0){return fPull(val,val);}
      bool   isConverged() const {return fConverged;}
      int    getIteration() const {return fIteration;}
      void   setCovariance(TMatrixD &val){V=val;}
      void   setMeasurement(TMatrixD &val){y=val;}
      FParticleCand getDaughter(int val);
      bool fit();
};

HFitter::HFitter(int n, std::vector<FParticleCand> cands){
    // check correct number
    if (n!=cands.size()){
         Fatal("HFitter::HFitter", "Wrong number of particle candidates n MUST EQUAL cands!");
    }
    // n is the number of daughters e.g. (L->ppi-) n=2
	y.ResizeTo(n*5,1);
	V.ResizeTo(n*5,n*5);

    y.Zero();
    V.Zero();

    fConverged = false;
    fIteration = 0;
    fNdf = 0;
    fN = n;

    // set 'y=alpha' measurements
    // and the covariance
    for (int ix=0; ix<fN; ix++){
        FParticleCand cand = cands[ix];
        y(0+ix*5, 0) = 1./cand.P();
        y(1+ix*5, 0) = cand.Theta();
        y(2+ix*5, 0) = cand.Phi();
        y(3+ix*5, 0) = cand.getR();
        y(4+ix*5, 0) = cand.getZ();
        fM.push_back(cand.M());

        // FIX ME: only for diagonal elements
        TMatrixD covariance = cand.getCovariance();
        V(0+ix*5, 0+ix*5) = covariance(0,0);
        V(1+ix*5, 1+ix*5) = covariance(1,1);
        V(2+ix*5, 2+ix*5) = covariance(2,2);
        V(3+ix*5, 3+ix*5) = covariance(3,3);
        V(4+ix*5, 4+ix*5) = covariance(4,4);
    }    

    fMassConstraint    = false;
    fMMConstraint      = false;
    fMassVtxConstraint = false;
    fVtxConstraint     = false;

}

void HFitter::addMassConstraint(double mass){
    fMass = mass;
    fNdf += 1;
    fMassConstraint = true;
}

void HFitter::addMissingMassConstraint(double mass, TLorentzVector init){
    fMass = mass;
    fInit = init;
    fNdf += 1;
    fMMConstraint = true; 
}

void HFitter::addMassVtxConstraint(double mass){
    fMass = mass;
    fNdf += 2;
    fMassVtxConstraint = true;
}

void HFitter::addVtxConstraint(){
    fNdf += 1;
    fVtxConstraint = true;
}

TMatrixD HFitter::f_eval(const TMatrixD &m_iter){
    TMatrixD d;
    // invariant mass constraint
    if (fMassConstraint){
        d.ResizeTo(1,1);
        double Px=0., Py=0., Pz=0., E=0.;
        for (int q=0; q<fN; q++){
            E  += std::sqrt( (1./m_iter(0+q*5,0))*(1./m_iter(0+q*5,0)) + fM[q]*fM[q]);
            Px += (1./m_iter(0+q*5,0))*std::sin(m_iter(1+q*5,0))*std::cos(m_iter(2+q*5,0));
            Py += (1./m_iter(0+q*5,0))*std::sin(m_iter(1+q*5,0))*std::sin(m_iter(2+q*5,0));
            Pz += (1./m_iter(0+q*5,0))*std::cos(m_iter(1+q*5,0));
        }
        
        d(0,0) = std::pow(E,2) - std::pow(Px,2) - std::pow(Py,2) - std::pow(Pz,2) - fMass*fMass; 
    }


    // invariant mass + vertex constraint
    if (fMassVtxConstraint){
        d.ResizeTo(2,1);
        double Px=0., Py=0., Pz=0., E=0.;
        for (int q=0; q<fN; q++){
            E  += std::sqrt( (1./m_iter(0+q*5,0))*(1./m_iter(0+q*5,0)) + fM[q]*fM[q]);
            Px += (1./m_iter(0+q*5,0))*std::sin(m_iter(1+q*5,0))*std::cos(m_iter(2+q*5,0));
            Py += (1./m_iter(0+q*5,0))*std::sin(m_iter(1+q*5,0))*std::sin(m_iter(2+q*5,0));
            Pz += (1./m_iter(0+q*5,0))*std::cos(m_iter(1+q*5,0));
        }
 
        TVector3 base_1, base_2, dir_1, dir_2;
        base_1.SetXYZ( m_iter(3+0*5,0)*std::cos(m_iter(2+0*5,0)+TMath::PiOver2()),
                       m_iter(3+0*5,0)*std::sin(m_iter(2+0*5,0)+TMath::PiOver2()), m_iter(4+0*5,0) );
        base_2.SetXYZ( m_iter(3+1*5,0)*std::cos(m_iter(2+1*5,0)+TMath::PiOver2()),
                       m_iter(3+1*5,0)*std::sin(m_iter(2+1*5,0)+TMath::PiOver2()), m_iter(4+1*5,0) );

        dir_1.SetXYZ(std::sin(m_iter(1+0*5,0))*std::cos(m_iter(2+0*5,0)),
                     std::sin(m_iter(1+0*5,0))*std::sin(m_iter(2+0*5,0)), std::cos(m_iter(1+0*5,0)));
        dir_2.SetXYZ(std::sin(m_iter(1+1*5,0))*std::cos(m_iter(2+1*5,0)),
                     std::sin(m_iter(1+1*5,0))*std::sin(m_iter(2+1*5,0)), std::cos(m_iter(1+1*5,0)));

        d(0,0) = std::pow(E,2) - std::pow(Px,2) - std::pow(Py,2) - std::pow(Pz,2) - fMass*fMass;
        d(1,0) = std::fabs((dir_1.Cross(dir_2)).Dot((base_2-base_1)));
    }

    // missing mass constraint
    if (fMMConstraint){
        d.ResizeTo(1,1);
        double Px=0., Py=0., Pz=0., E=0.;
        Px += fInit.Px();
        Py += fInit.Py();
        Pz += fInit.Pz();
        E  += fInit.E();
        for (int q=0; q<fN; q++){
            E  -= std::sqrt( (1./m_iter(0+q*5,0))*(1./m_iter(0+q*5,0)) + fM[q]*fM[q] );
            Px -= (1./m_iter(0+q*5,0))*std::sin(m_iter(1+q*5,0))*std::cos(m_iter(2+q*5,0));
            Py -= (1./m_iter(0+q*5,0))*std::sin(m_iter(1+q*5,0))*std::sin(m_iter(2+q*5,0));
            Pz -= (1./m_iter(0+q*5,0))*std::cos(m_iter(1+q*5,0));
        }
        d(0,0) = std::pow(E,2) - std::pow(Px,2) - std::pow(Py,2) - std::pow(Pz,2) - fMass*fMass; 
    }

    // vertex constraint
    if (fVtxConstraint){
        d.ResizeTo(1,1); 
        TVector3 base_1, base_2, dir_1, dir_2;
        base_1.SetXYZ( m_iter(3+0*5,0)*std::cos(m_iter(2+0*5,0)+TMath::PiOver2()),
                       m_iter(3+0*5,0)*std::sin(m_iter(2+0*5,0)+TMath::PiOver2()), m_iter(4+0*5,0) );
        base_2.SetXYZ( m_iter(3+1*5,0)*std::cos(m_iter(2+1*5,0)+TMath::PiOver2()),
                       m_iter(3+1*5,0)*std::sin(m_iter(2+1*5,0)+TMath::PiOver2()), m_iter(4+1*5,0) );

        dir_1.SetXYZ(std::sin(m_iter(1+0*5,0))*std::cos(m_iter(2+0*5,0)),
                     std::sin(m_iter(1+0*5,0))*std::sin(m_iter(2+0*5,0)), std::cos(m_iter(1+0*5,0)));
        dir_2.SetXYZ(std::sin(m_iter(1+1*5,0))*std::cos(m_iter(2+1*5,0)),
                     std::sin(m_iter(1+1*5,0))*std::sin(m_iter(2+1*5,0)), std::cos(m_iter(1+1*5,0)));

        d(0,0) = std::fabs((dir_1.Cross(dir_2)).Dot((base_1-base_2)));
    }        

    return d;
}

TMatrixD HFitter::Feta_eval(const TMatrixD &m_iter){

    TMatrixD H;
    // invariant mass constraint
    if (fMassConstraint){
        H.ResizeTo(1,fN*5);
        H.Zero();
        double Px=0., Py=0., Pz=0., E=0.;
        for (int q=0; q<fN; q++){
            E  += std::sqrt( (1./m_iter(0+q*5,0))*(1./m_iter(0+q*5,0)) + fM[q]*fM[q]);
            Px += (1./m_iter(0+q*5,0))*std::sin(m_iter(1+q*5,0))*std::cos(m_iter(2+q*5,0));
            Py += (1./m_iter(0+q*5,0))*std::sin(m_iter(1+q*5,0))*std::sin(m_iter(2+q*5,0));
            Pz += (1./m_iter(0+q*5,0))*std::cos(m_iter(1+q*5,0));
        }
        for (int q=0; q<fN; q++){
            double Pi = 1./m_iter(0+q*5,0);
            double Ei = std::sqrt( Pi*Pi + fM[q]*fM[q] );
            H(0, 0+q*5) = -2*E*(std::pow(Pi,3)/Ei) + 2*std::pow(Pi,2)*std::sin(m_iter(1+q*5,0))*std::cos(m_iter(2+q*5,0))*Px
                                                   + 2*std::pow(Pi,2)*std::sin(m_iter(1+q*5,0))*std::sin(m_iter(2+q*5,0))*Py
                                                   + 2*std::pow(Pi,2)*std::cos(m_iter(1+q*5,0))*Pz;                                              

            H(0, 1+q*5) = -2*Pi*std::cos(m_iter(1+q*5,0))*std::cos(m_iter(2+q*5,0))*Px
                        -  2*Pi*std::cos(m_iter(1+q*5,0))*std::sin(m_iter(2+q*5,0))*Py
                        +  2*Pi*std::sin(m_iter(1+q*5,0))*Pz;

            H(0, 2+q*5) = 2*Pi*std::sin(m_iter(1+q*5,0))*std::sin(m_iter(2+q*5,0))*Px
                        - 2*Pi*std::sin(m_iter(1+q*5,0))*std::cos(m_iter(2+q*5,0))*Py;

            H(0, 3+q*5) = 0.;
            H(0, 4+q*4) = 0.;
        }
    }

    // invariant mass + vertex constraint
    if (fMassVtxConstraint){
        H.ResizeTo(2,fN*5);
        H.Zero();
        double Px=0., Py=0., Pz=0., E=0.;
        for (int q=0; q<fN; q++){
            E  += std::sqrt( (1./m_iter(0+q*5,0))*(1./m_iter(0+q*5,0)) + fM[q]*fM[q]);
            Px += (1./m_iter(0+q*5,0))*std::sin(m_iter(1+q*5,0))*std::cos(m_iter(2+q*5,0));
            Py += (1./m_iter(0+q*5,0))*std::sin(m_iter(1+q*5,0))*std::sin(m_iter(2+q*5,0));
            Pz += (1./m_iter(0+q*5,0))*std::cos(m_iter(1+q*5,0));
        }
        for (int q=0; q<fN; q++){
            double Pi = 1./m_iter(0+q*5,0);
            double Ei = std::sqrt( Pi*Pi + fM[q]*fM[q] );
            H(0, 0+q*5) = -2*E*(std::pow(Pi,3)/Ei) + 2*std::pow(Pi,2)*std::sin(m_iter(1+q*5,0))*std::cos(m_iter(2+q*5,0))*Px
                                                   + 2*std::pow(Pi,2)*std::sin(m_iter(1+q*5,0))*std::sin(m_iter(2+q*5,0))*Py
                                                   + 2*std::pow(Pi,2)*std::cos(m_iter(1+q*5,0))*Pz;                                                  

            H(0, 1+q*5) = -2*Pi*std::cos(m_iter(1+q*5,0))*std::cos(m_iter(2+q*5,0))*Px
                        -  2*Pi*std::cos(m_iter(1+q*5,0))*std::sin(m_iter(2+q*5,0))*Py
                        +  2*Pi*std::sin(m_iter(1+q*5,0))*Pz;

            H(0, 2+q*5) = 2*Pi*std::sin(m_iter(1+q*5,0))*std::sin(m_iter(2+q*5,0))*Px
                        - 2*Pi*std::sin(m_iter(1+q*5,0))*std::cos(m_iter(2+q*5,0))*Py;
            H(0, 3+q*5) = 0.;
            H(0, 4+q*5) = 0.;
            H(1, 0+q*5) = 0.;
        }

        H(1,1) = m_iter(3,0)*std::cos(m_iter(2,0)+pi2)*std::cos(m_iter(1,0))*std::sin(m_iter(2,0))*std::cos(m_iter(6,0)) +
                 m_iter(3,0)*std::cos(m_iter(2,0)+pi2)*std::sin(m_iter(6,0))*std::sin(m_iter(7,0))*std::sin(m_iter(2,0)) -
                 m_iter(8,0)*std::cos(m_iter(7,0)+pi2)*std::cos(m_iter(1,0))*std::sin(m_iter(2,0))*std::cos(m_iter(6,0)) -
                 m_iter(8,0)*std::cos(m_iter(7,0)+pi2)*std::sin(m_iter(6,0))*std::sin(m_iter(7,0))*std::sin(m_iter(1,0)) -
                 m_iter(3,0)*std::sin(m_iter(2,0)+pi2)*std::cos(m_iter(1,0))*std::cos(m_iter(2,0))*std::cos(m_iter(6,0)) -
                 m_iter(3,0)*std::sin(m_iter(2,0)+pi2)*std::sin(m_iter(6,0))*std::cos(m_iter(7,0))*std::sin(m_iter(1,0)) +
                 m_iter(8,0)*std::sin(m_iter(7,0)+pi2)*std::cos(m_iter(1,0))*std::cos(m_iter(2,0))*std::cos(m_iter(6,0)) +
                 m_iter(8,0)*std::sin(m_iter(7,0)+pi2)*std::sin(m_iter(6,0))*std::cos(m_iter(7,0))*std::sin(m_iter(1,0)) +
                 (m_iter(4,0)-m_iter(9,0))*( std::cos(m_iter(1,0))*std::cos(m_iter(2,0))*std::sin(m_iter(6,0))*std::sin(m_iter(7,0)) - 
                                             std::cos(m_iter(1,0))*std::sin(m_iter(2,0))*std::sin(m_iter(6,0))*std::cos(m_iter(7,0)));


        H(1,6) = -1*m_iter(3,0)*std::cos(m_iter(2,0)+pi2)*std::sin(m_iter(1,0))*std::sin(m_iter(2,0))*std::sin(m_iter(6,0)) -
                 m_iter(3,0)*std::cos(m_iter(2,0)+pi2)*std::cos(m_iter(6,0))*std::sin(m_iter(7,0))*std::cos(m_iter(1,0)) + //
                 m_iter(8,0)*std::cos(m_iter(7,0)+pi2)*std::sin(m_iter(1,0))*std::sin(m_iter(2,0))*std::sin(m_iter(6,0)) +
                 m_iter(8,0)*std::cos(m_iter(7,0)+pi2)*std::cos(m_iter(6,0))*std::sin(m_iter(7,0))*std::cos(m_iter(1,0)) +
                 m_iter(3,0)*std::sin(m_iter(2,0)+pi2)*std::sin(m_iter(1,0))*std::cos(m_iter(2,0))*std::sin(m_iter(6,0)) +
                 m_iter(3,0)*std::sin(m_iter(2,0)+pi2)*std::cos(m_iter(6,0))*std::cos(m_iter(7,0))*std::cos(m_iter(1,0)) -
                 m_iter(8,0)*std::sin(m_iter(7,0)+pi2)*std::sin(m_iter(1,0))*std::cos(m_iter(2,0))*std::sin(m_iter(6,0)) -
                 m_iter(8,0)*std::sin(m_iter(7,0)+pi2)*std::cos(m_iter(6,0))*std::cos(m_iter(7,0))*std::cos(m_iter(1,0)) +
                 (m_iter(4,0)-m_iter(9,0))*( std::sin(m_iter(1,0))*std::cos(m_iter(2,0))*std::cos(m_iter(6,0))*std::sin(m_iter(7,0)) - 
                                             std::sin(m_iter(1,0))*std::sin(m_iter(2,0))*std::cos(m_iter(6,0))*std::cos(m_iter(7,0)));

        H(1,2) = m_iter(3,0)*std::cos(m_iter(2,0)+pi2)*std::sin(m_iter(1,0))*std::cos(m_iter(2,0))*std::cos(m_iter(6,0)) -
                 m_iter(3,0)*std::sin(m_iter(2,0)+pi2)*std::sin(m_iter(1,0))*std::sin(m_iter(2,0))*std::cos(m_iter(6,0)) +
                 m_iter(3,0)*std::sin(m_iter(2,0)+pi2)*std::sin(m_iter(6,0))*std::sin(m_iter(7,0))*std::cos(m_iter(1,0)) -
                 m_iter(8,0)*std::cos(m_iter(7,0)+pi2)*std::sin(m_iter(1,0))*std::cos(m_iter(2,0))*std::cos(m_iter(6,0)) +
                 m_iter(3,0)*std::sin(m_iter(2,0)+pi2)*std::sin(m_iter(1,0))*std::sin(m_iter(2,0))*std::cos(m_iter(6,0)) -
                 m_iter(3,0)*std::cos(m_iter(2,0)+pi2)*std::sin(m_iter(1,0))*std::cos(m_iter(2,0))*std::cos(m_iter(6,0)) +
                 m_iter(3,0)*std::cos(m_iter(2,0)+pi2)*std::sin(m_iter(6,0))*std::cos(m_iter(7,0))*std::cos(m_iter(1,0)) - //
                 m_iter(8,0)*std::sin(m_iter(7,0)+pi2)*std::sin(m_iter(1,0))*std::sin(m_iter(2,0))*std::cos(m_iter(6,0)) +
                 (m_iter(9,0)-m_iter(4,0))*( std::sin(m_iter(1,0))*std::sin(m_iter(2,0))*std::sin(m_iter(6,0))*std::sin(m_iter(7,0)) - 
                                             std::sin(m_iter(1,0))*std::cos(m_iter(2,0))*std::sin(m_iter(6,0))*std::cos(m_iter(7,0)));
                                             
        H(1,7) = -1*m_iter(3,0)*std::cos(m_iter(2,0)+pi2)*std::sin(m_iter(6,0))*std::cos(m_iter(7,0))*std::cos(m_iter(1,0)) +
                 m_iter(8,0)*std::sin(m_iter(7,0)+pi2)*std::sin(m_iter(1,0))*std::sin(m_iter(2,0))*std::cos(m_iter(6,0)) +
                 m_iter(8,0)*std::cos(m_iter(7,0)+pi2)*std::sin(m_iter(6,0))*std::cos(m_iter(7,0))*std::cos(m_iter(1,0)) -
                 m_iter(8,0)*std::sin(m_iter(7,0)+pi2)*std::sin(m_iter(6,0))*std::sin(m_iter(7,0))*std::cos(m_iter(1,0)) -
                 m_iter(3,0)*std::sin(m_iter(2,0)+pi2)*std::sin(m_iter(6,0))*std::sin(m_iter(7,0))*std::cos(m_iter(1,0)) +
                 m_iter(8,0)*std::cos(m_iter(7,0)+pi2)*std::sin(m_iter(1,0))*std::cos(m_iter(2,0))*std::cos(m_iter(6,0)) +
                 m_iter(8,0)*std::sin(m_iter(7,0)+pi2)*std::sin(m_iter(6,0))*std::cos(m_iter(7,0))*std::cos(m_iter(1,0)) -
                 m_iter(8,0)*std::cos(m_iter(7,0)+pi2)*std::sin(m_iter(6,0))*std::cos(m_iter(7,0))*std::cos(m_iter(1,0)) +
                 (m_iter(4,0)-m_iter(9,0))*( std::sin(m_iter(1,0))*std::cos(m_iter(2,0))*std::sin(m_iter(6,0))*std::cos(m_iter(7,0)) - 
                                             std::sin(m_iter(1,0))*std::sin(m_iter(2,0))*std::sin(m_iter(6,0))*std::sin(m_iter(7,0)));

        H(1,3) = std::cos(m_iter(2,0)+pi2)*( std::sin(m_iter(1,0))*std::sin(m_iter(2,0))*std::cos(m_iter(6,0)) - 
                                             std::sin(m_iter(6,0))*std::sin(m_iter(7,0))*std::cos(m_iter(1,0))) -
                 std::sin(m_iter(2,0)+pi2)*( std::sin(m_iter(1,0))*std::cos(m_iter(2,0))*std::cos(m_iter(6,0)) - 
                                             std::sin(m_iter(6,0))*std::cos(m_iter(7,0))*std::cos(m_iter(1,0)));

        H(1,8) = -1*std::cos(m_iter(7,0)+pi2)*( std::sin(m_iter(1,0))*std::sin(m_iter(2,0))*std::cos(m_iter(6,0)) - 
                                                std::sin(m_iter(6,0))*std::sin(m_iter(7,0))*std::cos(m_iter(1,0))) +
                    std::sin(m_iter(7,0)+pi2)*( std::sin(m_iter(1,0))*std::cos(m_iter(2,0))*std::cos(m_iter(6,0)) - 
                                                std::sin(m_iter(6,0))*std::cos(m_iter(7,0))*std::cos(m_iter(1,0)));

        H(1,4) = std::sin(m_iter(1,0))*std::cos(m_iter(2,0))*std::sin(m_iter(6,0))*std::sin(m_iter(7,0)) - 
                 std::sin(m_iter(1,0))*std::sin(m_iter(2,0))*std::sin(m_iter(6,0))*std::cos(m_iter(7,0));

        H(1,9) = -1*std::sin(m_iter(1,0))*std::cos(m_iter(2,0))*std::sin(m_iter(6,0))*std::sin(m_iter(7,0)) +
                    std::sin(m_iter(1,0))*std::sin(m_iter(2,0))*std::sin(m_iter(6,0))*std::cos(m_iter(7,0));             

    }

    // missing mass constraint
    if (fMMConstraint){
        H.ResizeTo(1,fN*5);
        H.Zero();
        double Px=0., Py=0., Pz=0., E=0.;
        Px += fInit.Px();
        Py += fInit.Py();
        Pz += fInit.Pz();
        E  += fInit.E();
        for (int q=0; q<fN; q++){
            E  -= std::sqrt( (1./m_iter(0+q*5,0))*(1./m_iter(0+q*5,0)) + fM[q]*fM[q]);
            Px -= (1./m_iter(0+q*5,0))*std::sin(m_iter(1+q*5,0))*std::cos(m_iter(2+q*5,0));
            Py -= (1./m_iter(0+q*5,0))*std::sin(m_iter(1+q*5,0))*std::sin(m_iter(2+q*5,0));
            Pz -= (1./m_iter(0+q*5,0))*std::cos(m_iter(1+q*5,0));
        }

        for (int q=0; q<fN; q++){
            double Pi = 1./m_iter(0+q*5,0);
            double Ei = std::sqrt( Pi*Pi + fM[q]*fM[q] );
            H(0, 0+q*5) =  2*E*(std::pow(Pi,3)/Ei) - 2*std::pow(Pi,2)*std::sin(m_iter(1+q*5,0))*std::cos(m_iter(2+q*5,0))*Px
                                                   - 2*std::pow(Pi,2)*std::sin(m_iter(1+q*5,0))*std::sin(m_iter(2+q*5,0))*Py
                                                   - 2*std::pow(Pi,2)*std::cos(m_iter(1+q*5,0))*Pz;                                              

            H(0, 1+q*5) =  2*Pi*std::cos(m_iter(1+q*5,0))*std::cos(m_iter(2+q*5,0))*Px
                        +  2*Pi*std::cos(m_iter(1+q*5,0))*std::sin(m_iter(2+q*5,0))*Py
                        -  2*Pi*std::sin(m_iter(1+q*5,0))*Pz;

            H(0, 2+q*5) = -2*Pi*std::sin(m_iter(1+q*5,0))*std::sin(m_iter(2+q*5,0))*Px
                         + 2*Pi*std::sin(m_iter(1+q*5,0))*std::cos(m_iter(2+q*5,0))*Py;
                                                                                                           
        }
    }


    // vertex constraint
    if (fVtxConstraint){
        H.ResizeTo(1,fN*5);
        H.Zero();

        H(0,0) = 0.;
        H(0,5) = 0.;

        H(0,1) = m_iter(3,0)*std::cos(m_iter(2,0)+pi2)*std::cos(m_iter(1,0))*std::sin(m_iter(2,0))*std::cos(m_iter(6,0)) +
                 m_iter(3,0)*std::cos(m_iter(2,0)+pi2)*std::sin(m_iter(6,0))*std::sin(m_iter(7,0))*std::sin(m_iter(2,0)) -
                 m_iter(8,0)*std::cos(m_iter(7,0)+pi2)*std::cos(m_iter(1,0))*std::sin(m_iter(2,0))*std::cos(m_iter(6,0)) -
                 m_iter(8,0)*std::cos(m_iter(7,0)+pi2)*std::sin(m_iter(6,0))*std::sin(m_iter(7,0))*std::sin(m_iter(1,0)) -
                 m_iter(3,0)*std::sin(m_iter(2,0)+pi2)*std::cos(m_iter(1,0))*std::cos(m_iter(2,0))*std::cos(m_iter(6,0)) -
                 m_iter(3,0)*std::sin(m_iter(2,0)+pi2)*std::sin(m_iter(6,0))*std::cos(m_iter(7,0))*std::sin(m_iter(1,0)) +
                 m_iter(8,0)*std::sin(m_iter(7,0)+pi2)*std::cos(m_iter(1,0))*std::cos(m_iter(2,0))*std::cos(m_iter(6,0)) +
                 m_iter(8,0)*std::sin(m_iter(7,0)+pi2)*std::sin(m_iter(6,0))*std::cos(m_iter(7,0))*std::sin(m_iter(1,0)) +
                 (m_iter(4,0)-m_iter(9,0))*( std::cos(m_iter(1,0))*std::cos(m_iter(2,0))*std::sin(m_iter(6,0))*std::sin(m_iter(7,0)) - 
                                             std::cos(m_iter(1,0))*std::sin(m_iter(2,0))*std::sin(m_iter(6,0))*std::cos(m_iter(7,0)));


        H(0,6) = -1*m_iter(3,0)*std::cos(m_iter(2,0)+pi2)*std::sin(m_iter(1,0))*std::sin(m_iter(2,0))*std::sin(m_iter(6,0)) -
                 m_iter(3,0)*std::cos(m_iter(2,0)+pi2)*std::cos(m_iter(6,0))*std::sin(m_iter(7,0))*std::cos(m_iter(1,0)) + //
                 m_iter(8,0)*std::cos(m_iter(7,0)+pi2)*std::sin(m_iter(1,0))*std::sin(m_iter(2,0))*std::sin(m_iter(6,0)) +
                 m_iter(8,0)*std::cos(m_iter(7,0)+pi2)*std::cos(m_iter(6,0))*std::sin(m_iter(7,0))*std::cos(m_iter(1,0)) +
                 m_iter(3,0)*std::sin(m_iter(2,0)+pi2)*std::sin(m_iter(1,0))*std::cos(m_iter(2,0))*std::sin(m_iter(6,0)) +
                 m_iter(3,0)*std::sin(m_iter(2,0)+pi2)*std::cos(m_iter(6,0))*std::cos(m_iter(7,0))*std::cos(m_iter(1,0)) -
                 m_iter(8,0)*std::sin(m_iter(7,0)+pi2)*std::sin(m_iter(1,0))*std::cos(m_iter(2,0))*std::sin(m_iter(6,0)) -
                 m_iter(8,0)*std::sin(m_iter(7,0)+pi2)*std::cos(m_iter(6,0))*std::cos(m_iter(7,0))*std::cos(m_iter(1,0)) +
                 (m_iter(4,0)-m_iter(9,0))*( std::sin(m_iter(1,0))*std::cos(m_iter(2,0))*std::cos(m_iter(6,0))*std::sin(m_iter(7,0)) - 
                                             std::sin(m_iter(1,0))*std::sin(m_iter(2,0))*std::cos(m_iter(6,0))*std::cos(m_iter(7,0)));

        H(0,2) = m_iter(3,0)*std::cos(m_iter(2,0)+pi2)*std::sin(m_iter(1,0))*std::cos(m_iter(2,0))*std::cos(m_iter(6,0)) -
                 m_iter(3,0)*std::sin(m_iter(2,0)+pi2)*std::sin(m_iter(1,0))*std::sin(m_iter(2,0))*std::cos(m_iter(6,0)) +
                 m_iter(3,0)*std::sin(m_iter(2,0)+pi2)*std::sin(m_iter(6,0))*std::sin(m_iter(7,0))*std::cos(m_iter(1,0)) -
                 m_iter(8,0)*std::cos(m_iter(7,0)+pi2)*std::sin(m_iter(1,0))*std::cos(m_iter(2,0))*std::cos(m_iter(6,0)) +
                 m_iter(3,0)*std::sin(m_iter(2,0)+pi2)*std::sin(m_iter(1,0))*std::sin(m_iter(2,0))*std::cos(m_iter(6,0)) -
                 m_iter(3,0)*std::cos(m_iter(2,0)+pi2)*std::sin(m_iter(1,0))*std::cos(m_iter(2,0))*std::cos(m_iter(6,0)) +
                 m_iter(3,0)*std::cos(m_iter(2,0)+pi2)*std::sin(m_iter(6,0))*std::cos(m_iter(7,0))*std::cos(m_iter(1,0)) - //
                 m_iter(8,0)*std::sin(m_iter(7,0)+pi2)*std::sin(m_iter(1,0))*std::sin(m_iter(2,0))*std::cos(m_iter(6,0)) +
                 (m_iter(9,0)-m_iter(4,0))*( std::sin(m_iter(1,0))*std::sin(m_iter(2,0))*std::sin(m_iter(6,0))*std::sin(m_iter(7,0)) - 
                                             std::sin(m_iter(1,0))*std::cos(m_iter(2,0))*std::sin(m_iter(6,0))*std::cos(m_iter(7,0)));
                                             
        H(0,7) = -1*m_iter(3,0)*std::cos(m_iter(2,0)+pi2)*std::sin(m_iter(6,0))*std::cos(m_iter(7,0))*std::cos(m_iter(1,0)) +
                 m_iter(8,0)*std::sin(m_iter(7,0)+pi2)*std::sin(m_iter(1,0))*std::sin(m_iter(2,0))*std::cos(m_iter(6,0)) +
                 m_iter(8,0)*std::cos(m_iter(7,0)+pi2)*std::sin(m_iter(6,0))*std::cos(m_iter(7,0))*std::cos(m_iter(1,0)) -
                 m_iter(8,0)*std::sin(m_iter(7,0)+pi2)*std::sin(m_iter(6,0))*std::sin(m_iter(7,0))*std::cos(m_iter(1,0)) -
                 m_iter(3,0)*std::sin(m_iter(2,0)+pi2)*std::sin(m_iter(6,0))*std::sin(m_iter(7,0))*std::cos(m_iter(1,0)) +
                 m_iter(8,0)*std::cos(m_iter(7,0)+pi2)*std::sin(m_iter(1,0))*std::cos(m_iter(2,0))*std::cos(m_iter(6,0)) +
                 m_iter(8,0)*std::sin(m_iter(7,0)+pi2)*std::sin(m_iter(6,0))*std::cos(m_iter(7,0))*std::cos(m_iter(1,0)) -
                 m_iter(8,0)*std::cos(m_iter(7,0)+pi2)*std::sin(m_iter(6,0))*std::cos(m_iter(7,0))*std::cos(m_iter(1,0)) +
                 (m_iter(4,0)-m_iter(9,0))*( std::sin(m_iter(1,0))*std::cos(m_iter(2,0))*std::sin(m_iter(6,0))*std::cos(m_iter(7,0)) - 
                                             std::sin(m_iter(1,0))*std::sin(m_iter(2,0))*std::sin(m_iter(6,0))*std::sin(m_iter(7,0)));

        H(0,3) = std::cos(m_iter(2,0)+pi2)*( std::sin(m_iter(1,0))*std::sin(m_iter(2,0))*std::cos(m_iter(6,0)) - 
                                             std::sin(m_iter(6,0))*std::sin(m_iter(7,0))*std::cos(m_iter(1,0))) -
                 std::sin(m_iter(2,0)+pi2)*( std::sin(m_iter(1,0))*std::cos(m_iter(2,0))*std::cos(m_iter(6,0)) - 
                                             std::sin(m_iter(6,0))*std::cos(m_iter(7,0))*std::cos(m_iter(1,0)));

        H(0,8) = -1*std::cos(m_iter(7,0)+pi2)*( std::sin(m_iter(1,0))*std::sin(m_iter(2,0))*std::cos(m_iter(6,0)) - 
                                                std::sin(m_iter(6,0))*std::sin(m_iter(7,0))*std::cos(m_iter(1,0))) +
                    std::sin(m_iter(7,0)+pi2)*( std::sin(m_iter(1,0))*std::cos(m_iter(2,0))*std::cos(m_iter(6,0)) - 
                                                std::sin(m_iter(6,0))*std::cos(m_iter(7,0))*std::cos(m_iter(1,0)));

        H(0,4) = std::sin(m_iter(1,0))*std::cos(m_iter(2,0))*std::sin(m_iter(6,0))*std::sin(m_iter(7,0)) - 
                 std::sin(m_iter(1,0))*std::sin(m_iter(2,0))*std::sin(m_iter(6,0))*std::cos(m_iter(7,0));

        H(0,9) = -1*std::sin(m_iter(1,0))*std::cos(m_iter(2,0))*std::sin(m_iter(6,0))*std::sin(m_iter(7,0)) +
                    std::sin(m_iter(1,0))*std::sin(m_iter(2,0))*std::sin(m_iter(6,0))*std::cos(m_iter(7,0));             

    }    

    return H;
}

bool HFitter::fit(){
	double lr = 0.5;
	TMatrixD alpha0(fN*5,1), alpha(fN*5,1);
    TMatrixD A0(y), V0(V);
	alpha0 = y;
	alpha  = alpha0;
	double chi2=1e6;
    TMatrixD D = Feta_eval(alpha);
    TMatrixD d = f_eval(alpha);

	for (int q=0; q<5; q++){
        TMatrixD DT(D.GetNcols(), D.GetNrows());
		DT.Transpose(D);
		TMatrixD VD = D*V*DT;
		VD.Invert();

		TMatrixD delta_alpha = alpha - alpha0;
		TMatrixD lambda = VD*D*delta_alpha + VD*d;
		TMatrixD lambdaT(lambda.GetNcols(), lambda.GetNrows());
		lambdaT.Transpose(lambda);
		TMatrixD neu_alpha(fN*5,1);
		neu_alpha = alpha - lr*V*DT*lambda;

        double chisqrd = 0.;
        
        for (int p=0; p<lambda.GetNrows(); p++){
            chisqrd = lambdaT(0,p)*d(p,0);
        }

        /* for checking convergence
        // three parameters are checked
        // 1. difference between measurements (for successive iterations) y
        // 2. difference between constraints (for successive iterations)  d
        // 3. difference between chi2 (for successive iterations)  chisqrd
        // check converge for 'y' measurements
        double sum0 = 0;
        for(int p=0; p<(fN*5); p++){
            sum0 += (neu_alpha(p,0)-alpha(p,0))*(neu_alpha(p,0)-alpha(p,0));
        }

        double d_const = fabs(d(0,0));
		if(fabs(chi2-chisqrd)<1e-3 && d_const<10 && sqrt(sum0)<1e-3){
            fIteration = q;
            fConverged = true;
            break;
        }
        */
		chi2 = chisqrd;		
		alpha0 = alpha;
		alpha = neu_alpha;
		V     = V - lr*V*DT*VD*D*V;
		D     = Feta_eval(alpha);
		d     = f_eval(alpha);
	}


	y = alpha;
	fChi2 = chi2;
	fProb = TMath::Prob(chi2, fNdf);

    // -----------------------------------------
    // Pull
    // -----------------------------------------
    fPull.ResizeTo(fN*5,fN*5);
    for (int b=0; b<(fN*5); b++) fPull(b,b) = -10000;

    if (true){
        for (int b=0; b<(fN*5); b++){
            double num = A0(b,0) - alpha(b,0);
            double dem = V0(b,b) - V(b,b);
            if (dem>0){
                fPull(b,b) = num/std::sqrt(dem);
            }
        }
    }

    //return fConverged; // for number of iterations greater than 1
    return true; // for number of iterations equal to 1
}

FParticleCand HFitter::getDaughter(int val=0){
    FParticleCand cand;
    double Px = (1./y(0+val*5,0))*std::sin(y(1+val*5,0))*std::cos(y(2+val*5,0));
    double Py = (1./y(0+val*5,0))*std::sin(y(1+val*5,0))*std::sin(y(2+val*5,0));
    double Pz = (1./y(0+val*5,0))*std::cos(y(1+val*5,0));
    double M  = fM[val];
    cand.SetXYZM(Px, Py, Pz, M);
    cand.setR(y(3+val*5,0));
    cand.setZ(y(4+val*5,0));

    // ---------------------------------------------------------------------------
    // set covariance
    // ---------------------------------------------------------------------------
    TMatrixD cov(5,5);
    cov(0,0) = V(0+val*5,0+val*5);
    cov(1,1) = V(1+val*5,1+val*5);
    cov(2,2) = V(2+val*5,2+val*5);
    cov(3,3) = V(3+val*5,3+val*5);
    cov(4,4) = V(4+val*5,4+val*5);
    cand.setCovariance(cov);
    // ---------------------------------------------------------------------------    

    return cand;
}
