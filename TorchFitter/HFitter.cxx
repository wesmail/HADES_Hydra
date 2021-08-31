#include "HFitter.h"
#include <ATen/ATen.h>
#include <torch/torch.h>

ClassImp(HFitter);

void eval(at::Tensor y_, at::Tensor& M, at::Tensor& c, Int_t n, std::vector<double> m, double im, double mm, TLorentzVector init, TString key) {
    Int_t size = 5;
    // set measurements to be differentiable
    auto y = y_.detach().clone();
    y.set_requires_grad(true);

    if (key == "Mass") {
        auto tE = torch::zeros({1, 1});
        auto tPx = torch::zeros({1, 1});
        auto tPy = torch::zeros({1, 1});
        auto tPz = torch::zeros({1, 1});
        // ----------------------------------------------------------------------
        // constraint equations
        // ----------------------------------------------------------------------
        for (int q = 0; q < n; q++) {
            tE -= torch::sqrt((1. / y[0 + q * size][0]) * (1. / y[0 + q * size][0]) + m[q] * m[q]);
            tPx -= (1. / y[0 + q * size][0]) * torch::sin(y[1 + q * size][0]) * torch::cos(y[2 + q * size][0]);
            tPy -= (1. / y[0 + q * size][0]) * torch::sin(y[1 + q * size][0]) * torch::sin(y[2 + q * size][0]);
            tPz -= (1. / y[0 + q * size][0]) * torch::cos(y[1 + q * size][0]);
        }

        auto f1 = at::pow(tE, 2) - at::pow(tPx, 2) - at::pow(tPy, 2) - at::pow(tPz, 2) - (im * im);
        c.narrow(0, 0, 1) = f1;
        f1.backward();                                    // perform the partial derivate of f() w.r.t measured variables
        M.narrow(0, 0, 1) = y.grad().view({1, n * size}); // store the partial derivates w.r.t measured variables
        y.grad().zero_();                                 // zero out the gradient
    }

    else if (key == "2C") {
        auto t1E = torch::zeros({1, 1});
        auto t1Px = torch::zeros({1, 1});
        auto t1Py = torch::zeros({1, 1});
        auto t1Pz = torch::zeros({1, 1});
        auto t2E = torch::zeros({1, 1});
        auto t2Px = torch::zeros({1, 1});
        auto t2Py = torch::zeros({1, 1});
        auto t2Pz = torch::zeros({1, 1});
        t2Px += init.Px();
        t2Py += init.Py();
        t2Pz += init.Pz();
        t2E += init.E();
        // ----------------------------------------------------------------------
        // constraint equations
        // ----------------------------------------------------------------------
        for (int q = 0; q < n; q++) {
            if (q < 2) {
                t1E += torch::sqrt((1. / y[0 + q * size][0]) * (1. / y[0 + q * size][0]) + m[q] * m[q]);
                t1Px += (1. / y[0 + q * size][0]) * torch::sin(y[1 + q * size][0]) * torch::cos(y[2 + q * size][0]);
                t1Py += (1. / y[0 + q * size][0]) * torch::sin(y[1 + q * size][0]) * torch::sin(y[2 + q * size][0]);
                t1Pz += (1. / y[0 + q * size][0]) * torch::cos(y[1 + q * size][0]);
            }
            t2E -= torch::sqrt((1. / y[0 + q * size][0]) * (1. / y[0 + q * size][0]) + m[q] * m[q]);
            t2Px -= (1. / y[0 + q * size][0]) * torch::sin(y[1 + q * size][0]) * torch::cos(y[2 + q * size][0]);
            t2Py -= (1. / y[0 + q * size][0]) * torch::sin(y[1 + q * size][0]) * torch::sin(y[2 + q * size][0]);
            t2Pz -= (1. / y[0 + q * size][0]) * torch::cos(y[1 + q * size][0]);
        }

        auto f1 = at::pow(t1E, 2) - at::pow(t1Px, 2) - at::pow(t1Py, 2) - at::pow(t1Pz, 2) - (im * im);
        auto f2 = at::pow(t2E, 2) - at::pow(t2Px, 2) - at::pow(t2Py, 2) - at::pow(t2Pz, 2) - (mm * mm);

        c.narrow(0, 0, 1) = f1;
        c.narrow(0, 1, 1) = f2;

        f1.backward();                                    // perform the partial derivate of f() w.r.t measured variables
        M.narrow(0, 0, 1) = y.grad().view({1, n * size}); // store the partial derivates w.r.t measured variables
        y.grad().zero_();                                 // zero out the gradient
        f2.backward();                                    // perform the partial derivate of f() w.r.t measured variables
        M.narrow(0, 1, 1) = y.grad().view({1, n * size}); // store the partial derivates w.r.t measured variables
        y.grad().zero_();                                 // zero out the gradient
    }

    else if (key == "3C") {
        auto t1E = torch::zeros({1, 1});
        auto t1Px = torch::zeros({1, 1});
        auto t1Py = torch::zeros({1, 1});
        auto t1Pz = torch::zeros({1, 1});
        auto t2E = torch::zeros({1, 1});
        auto t2Px = torch::zeros({1, 1});
        auto t2Py = torch::zeros({1, 1});
        auto t2Pz = torch::zeros({1, 1});
        t2Px += init.Px();
        t2Py += init.Py();
        t2Pz += init.Pz();
        t2E += init.E();
        // ----------------------------------------------------------------------
        // constraint equations
        // ----------------------------------------------------------------------
        for (int q = 0; q < n; q++) {
            if (q < 2) {
                t1E += torch::sqrt((1. / y[0 + q * size][0]) * (1. / y[0 + q * size][0]) + m[q] * m[q]);
                t1Px += (1. / y[0 + q * size][0]) * torch::sin(y[1 + q * size][0]) * torch::cos(y[2 + q * size][0]);
                t1Py += (1. / y[0 + q * size][0]) * torch::sin(y[1 + q * size][0]) * torch::sin(y[2 + q * size][0]);
                t1Pz += (1. / y[0 + q * size][0]) * torch::cos(y[1 + q * size][0]);
            }
            t2E -= torch::sqrt((1. / y[0 + q * size][0]) * (1. / y[0 + q * size][0]) + m[q] * m[q]);
            t2Px -= (1. / y[0 + q * size][0]) * torch::sin(y[1 + q * size][0]) * torch::cos(y[2 + q * size][0]);
            t2Py -= (1. / y[0 + q * size][0]) * torch::sin(y[1 + q * size][0]) * torch::sin(y[2 + q * size][0]);
            t2Pz -= (1. / y[0 + q * size][0]) * torch::cos(y[1 + q * size][0]);
        }

        auto B1 = torch::zeros({1, 3});
        auto B2 = torch::zeros({1, 3});
        auto D1 = torch::zeros({1, 3});
        auto D2 = torch::zeros({3, 1});
        B1[0][0] = y[3 + 2 * size][0] * (torch::cos(y[2 + 2 * size][0] + TMath::PiOver2()));
        B1[0][1] = y[3 + 2 * size][0] * (torch::sin(y[2 + 2 * size][0] + TMath::PiOver2()));
        B1[0][2] = y[4 + 2 * size][0];
        B2[0][0] = y[3 + 3 * size][0] * (torch::cos(y[2 + 3 * size][0] + TMath::PiOver2()));
        B2[0][1] = y[3 + 3 * size][0] * (torch::sin(y[2 + 3 * size][0] + TMath::PiOver2()));
        B2[0][2] = y[4 + 3 * size][0];
        D1[0][0] = torch::sin(y[1 + 2 * size][0]) * torch::cos(y[2 + 2 * size][0]); // v1x
        D1[0][1] = torch::sin(y[1 + 2 * size][0]) * torch::sin(y[2 + 2 * size][0]); // v1y
        D1[0][2] = torch::cos(y[1 + 2 * size][0]);                                  // v1z
        D2[0][0] = torch::sin(y[1 + 3 * size][0]) * torch::cos(y[2 + 3 * size][0]); // v2x
        D2[1][0] = torch::sin(y[1 + 3 * size][0]) * torch::sin(y[2 + 3 * size][0]); // v2y
        D2[2][0] = torch::cos(y[1 + 3 * size][0]);                                  // v2z

        auto f1 = at::pow(t1E, 2) - at::pow(t1Px, 2) - at::pow(t1Py, 2) - at::pow(t1Pz, 2) - (im * im);
        auto f2 = torch::abs(torch::mm(B2 - B1, torch::cross(D1.t(), D2))) / torch::norm(torch::cross(D1.t(), D2));
        auto f3 = at::pow(t2E, 2) - at::pow(t2Px, 2) - at::pow(t2Py, 2) - at::pow(t2Pz, 2) - (mm * mm);

        c.narrow(0, 0, 1) = f1;
        c.narrow(0, 1, 1) = f2;
        c.narrow(0, 2, 1) = f3;

        f1.backward();                                    // perform the partial derivate of f() w.r.t measured variables
        M.narrow(0, 0, 1) = y.grad().view({1, n * size}); // store the partial derivates w.r.t measured variables
        y.grad().zero_();                                 // zero out the gradient
        f2.backward();                                    // perform the partial derivate of f() w.r.t measured variables
        M.narrow(0, 1, 1) = y.grad().view({1, n * size}); // store the partial derivates w.r.t measured variables
        y.grad().zero_();                                 // zero out the gradient
        f3.backward();                                    // perform the partial derivate of f() w.r.t measured variables
        M.narrow(0, 2, 1) = y.grad().view({1, n * size}); // store the partial derivates w.r.t measured variables
        y.grad().zero_();                                 // zero out the gradient
    } else if (key = "4C") {
        auto t1E = torch::zeros({1, 1});
        auto t1Px = torch::zeros({1, 1});
        auto t1Py = torch::zeros({1, 1});
        auto t1Pz = torch::zeros({1, 1});
        auto t2E = torch::zeros({1, 1});
        auto t2Px = torch::zeros({1, 1});
        auto t2Py = torch::zeros({1, 1});
        auto t2Pz = torch::zeros({1, 1});
        auto t3E = torch::zeros({1, 1});
        auto t3Px = torch::zeros({1, 1});
        auto t3Py = torch::zeros({1, 1});
        auto t3Pz = torch::zeros({1, 1});
        t2Px += init.Px();
        t2Py += init.Py();
        t2Pz += init.Pz();
        t2E += init.E();
        t3Px += init.Px();
        t3Py += init.Py();
        t3Pz += init.Pz();
        t3E += init.E();
        // ----------------------------------------------------------------------
        // constraint equations
        // ----------------------------------------------------------------------
        for (int q = 0; q < n; q++) {
            if (q < 2) {
                t1E += torch::sqrt((1. / y[0 + q * size][0]) * (1. / y[0 + q * size][0]) + m[q] * m[q]);
                t1Px += (1. / y[0 + q * size][0]) * torch::sin(y[1 + q * size][0]) * torch::cos(y[2 + q * size][0]);
                t1Py += (1. / y[0 + q * size][0]) * torch::sin(y[1 + q * size][0]) * torch::sin(y[2 + q * size][0]);
                t1Pz += (1. / y[0 + q * size][0]) * torch::cos(y[1 + q * size][0]);
            }

            if (q >= 2) {
                t3E -= torch::sqrt((1. / y[0 + q * size][0]) * (1. / y[0 + q * size][0]) + m[q] * m[q]);
                t3Px -= (1. / y[0 + q * size][0]) * torch::sin(y[1 + q * size][0]) * torch::cos(y[2 + q * size][0]);
                t3Py -= (1. / y[0 + q * size][0]) * torch::sin(y[1 + q * size][0]) * torch::sin(y[2 + q * size][0]);
                t3Pz -= (1. / y[0 + q * size][0]) * torch::cos(y[1 + q * size][0]);
            }
            t2E -= torch::sqrt((1. / y[0 + q * size][0]) * (1. / y[0 + q * size][0]) + m[q] * m[q]);
            t2Px -= (1. / y[0 + q * size][0]) * torch::sin(y[1 + q * size][0]) * torch::cos(y[2 + q * size][0]);
            t2Py -= (1. / y[0 + q * size][0]) * torch::sin(y[1 + q * size][0]) * torch::sin(y[2 + q * size][0]);
            t2Pz -= (1. / y[0 + q * size][0]) * torch::cos(y[1 + q * size][0]);
        }

        auto b1 = torch::zeros({1, 3});
        auto b2 = torch::zeros({1, 3});
        auto d1 = torch::zeros({1, 3});
        auto d2 = torch::zeros({3, 1});
        b1[0][0] = y[3 + 0 * size][0] * (torch::cos(y[2 + 0 * size][0] + TMath::PiOver2()));
        b1[0][1] = y[3 + 0 * size][0] * (torch::sin(y[2 + 0 * size][0] + TMath::PiOver2()));
        b1[0][2] = y[4 + 0 * size][0];
        b2[0][0] = y[3 + 1 * size][0] * (torch::cos(y[2 + 1 * size][0] + TMath::PiOver2()));
        b2[0][1] = y[3 + 1 * size][0] * (torch::sin(y[2 + 1 * size][0] + TMath::PiOver2()));
        b2[0][2] = y[4 + 1 * size][0];
        d1[0][0] = torch::sin(y[1 + 0 * size][0]) * torch::cos(y[2 + 0 * size][0]); // v1x
        d1[0][1] = torch::sin(y[1 + 0 * size][0]) * torch::sin(y[2 + 0 * size][0]); // v1y
        d1[0][2] = torch::cos(y[1 + 0 * size][0]);                                  // v1z
        d2[0][0] = torch::sin(y[1 + 1 * size][0]) * torch::cos(y[2 + 1 * size][0]); // v2x
        d2[1][0] = torch::sin(y[1 + 1 * size][0]) * torch::sin(y[2 + 1 * size][0]); // v2y
        d2[2][0] = torch::cos(y[1 + 1 * size][0]);                                  // v2z

        auto B1 = torch::zeros({1, 3});
        auto B2 = torch::zeros({1, 3});
        auto D1 = torch::zeros({1, 3});
        auto D2 = torch::zeros({3, 1});
        B1[0][0] = y[3 + 2 * size][0] * (torch::cos(y[2 + 2 * size][0] + TMath::PiOver2()));
        B1[0][1] = y[3 + 2 * size][0] * (torch::sin(y[2 + 2 * size][0] + TMath::PiOver2()));
        B1[0][2] = y[4 + 2 * size][0];
        B2[0][0] = y[3 + 3 * size][0] * (torch::cos(y[2 + 3 * size][0] + TMath::PiOver2()));
        B2[0][1] = y[3 + 3 * size][0] * (torch::sin(y[2 + 3 * size][0] + TMath::PiOver2()));
        B2[0][2] = y[4 + 3 * size][0];
        D1[0][0] = torch::sin(y[1 + 2 * size][0]) * torch::cos(y[2 + 2 * size][0]); // v1x
        D1[0][1] = torch::sin(y[1 + 2 * size][0]) * torch::sin(y[2 + 2 * size][0]); // v1y
        D1[0][2] = torch::cos(y[1 + 2 * size][0]);                                  // v1z
        D2[0][0] = torch::sin(y[1 + 3 * size][0]) * torch::cos(y[2 + 3 * size][0]); // v2x
        D2[1][0] = torch::sin(y[1 + 3 * size][0]) * torch::sin(y[2 + 3 * size][0]); // v2y
        D2[2][0] = torch::cos(y[1 + 3 * size][0]);                                  // v2z

        auto f1 = at::pow(t1E, 2) - at::pow(t1Px, 2) - at::pow(t1Py, 2) - at::pow(t1Pz, 2) - (im * im);
        auto f2 = torch::abs(torch::mm(b2 - b1, torch::cross(d1.t(), d2))) / torch::norm(torch::cross(d1.t(), d2));
        auto f3 = torch::abs(torch::mm(B2 - B1, torch::cross(D1.t(), D2))) / torch::norm(torch::cross(D1.t(), D2));
        auto f4 = at::pow(t2E, 2) - at::pow(t2Px, 2) - at::pow(t2Py, 2) - at::pow(t2Pz, 2) - (mm * mm);

        c.narrow(0, 0, 1) = f1;
        c.narrow(0, 1, 1) = f2;
        c.narrow(0, 2, 1) = f3;
        c.narrow(0, 3, 1) = f4;

        f1.backward();                                    // perform the partial derivate of f() w.r.t measured variables
        M.narrow(0, 0, 1) = y.grad().view({1, n * size}); // store the partial derivates w.r.t measured variables
        y.grad().zero_();                                 // zero out the gradient
        f2.backward();                                    // perform the partial derivate of f() w.r.t measured variables
        M.narrow(0, 1, 1) = y.grad().view({1, n * size}); // store the partial derivates w.r.t measured variables
        y.grad().zero_();                                 // zero out the gradient
        f3.backward();                                    // perform the partial derivate of f() w.r.t measured variables
        M.narrow(0, 2, 1) = y.grad().view({1, n * size}); // store the partial derivates w.r.t measured variables
        y.grad().zero_();                                 // zero out the gradient
        f4.backward();                                    // perform the partial derivate of f() w.r.t measured variables
        M.narrow(0, 3, 1) = y.grad().view({1, n * size}); // store the partial derivates w.r.t measured variables
        y.grad().zero_();                                 // zero out the gradient
    }
}

HFitter::HFitter(int n) {
    fSize = 5;
    fNdf = 0;
    fN = n; // n is the number of daughters e.g. (L->ppi-) n=2

    fMassConstraint = false;
    f2CConstraint = false;
    f3CConstraint = false;
    f4CConstraint = false;

    fA.ResizeTo(fN * fSize, 1);
    fV.ResizeTo(fN * fSize, fN * fSize);
    fA.Zero();
    fV.Zero();
    // -----------------------------------------
    // Pull
    // -----------------------------------------
    fPull.ResizeTo(fN * fSize, fN * fSize);
    for (int b = 0; b < (fN * fSize); b++)
        fPull(b, b) = -10000;

    // -----------------------------------------
    // MC-Truth
    // -----------------------------------------
    // for (int b = 0; b < fN; b++)
    //    fPid[b] = -1;
}

void HFitter::addMassConstraint(double mass) {
    fMass = mass;
    fNdf += 1;
    fMassConstraint = true;
    fDim = 1;
}

void HFitter::add2CConstraint(float inv_mass, float miss_mass, TLorentzVector init) {
    fIMass = inv_mass;
    fMMass = miss_mass;
    fInit = init;
    fNdf += 2;
    f2CConstraint = true;
    fDim = 2;
}

void HFitter::add3CConstraint(float inv_mass, float miss_mass, TLorentzVector init) {
    fIMass = inv_mass;
    fMMass = miss_mass;
    fInit = init;
    fNdf += 3;
    f3CConstraint = true;
    fDim = 3;
}

void HFitter::add4CConstraint(float inv_mass, float miss_mass, TLorentzVector init) {
    fIMass = inv_mass;
    fMMass = miss_mass;
    fInit = init;
    fNdf += 4;
    f4CConstraint = true;
    fDim = 4;
}

bool HFitter::fit(std::vector<HRefitCand> cands) {
    // ----------------------------------------------------------------------
    // set values
    // ----------------------------------------------------------------------
    // set 'y=alpha' measurements
    // and the covariance
    auto y = torch::zeros({fN * fSize, 1}); // measurement vector
    auto V = torch::zeros(fN * fSize);      // covariance matrix
    for (int ix = 0; ix < fN; ix++) {
        HRefitCand cand = cands[ix];
        y.narrow(0, 0 + ix * fSize, 1) = 1. / cand.P(); // set values
        y.narrow(0, 1 + ix * fSize, 1) = cand.Theta();  // set values
        y.narrow(0, 2 + ix * fSize, 1) = cand.Phi();    // set values
        y.narrow(0, 3 + ix * fSize, 1) = cand.getR();   // set values
        y.narrow(0, 4 + ix * fSize, 1) = cand.getZ();   // set values
        fM.push_back(cand.M());

        // FIX ME: only for diagonal elements
        TMatrixD covariance = cand.getCovariance();
        V.narrow(0, 0 + ix * fSize, 1) = covariance(0, 0);
        V.narrow(0, 1 + ix * fSize, 1) = covariance(1, 1);
        V.narrow(0, 2 + ix * fSize, 1) = covariance(2, 2);
        V.narrow(0, 3 + ix * fSize, 1) = covariance(3, 3);
        V.narrow(0, 4 + ix * fSize, 1) = covariance(4, 4);

        // MC-Truth
        fNr = cand.getEventNr();
        fPid.push_back(cand.getPid());
        TVector3 tP(cand.getGMomentum().X(), cand.getGMomentum().Y(), cand.getGMomentum().Z());
        fGMomentum.push_back(tP);
    }
    // square covariance matrix
    V = torch::diag(V);
    // set measurements to be differentiable
    auto alpha0 = y.detach().clone(); // measurement vector
    auto alpha = y.detach().clone();  // measurement vector
    // y.set_requires_grad(true);
    // ----------------------------------------------------------------------
    // constraint equation as a matrix
    auto f = torch::zeros({fDim, 1});
    auto D = torch::zeros({fDim, fN * fSize}); // Jacobian matrix w.r.t measured variables
    // ----------------------------------------------------------------------
    // constraint equations
    // ----------------------------------------------------------------------
    if (fMassConstraint) {
        eval(y, D, f, fN, fM, fIMass, fMMass, fInit, "Mass");
    }
    if (f2CConstraint) {
        eval(y, D, f, fN, fM, fIMass, fMMass, fInit, "2C");
    }
    if (f3CConstraint) {
        eval(y, D, f, fN, fM, fIMass, fMMass, fInit, "3C");
    }
    if (f4CConstraint) {
        eval(y, D, f, fN, fM, fIMass, fMMass, fInit, "4C");
    }

    // ----------------------------------------------------------------------
    // fit the data here
    // ----------------------------------------------------------------------
    auto V0 = V.clone(); // covariance matrix
    auto d = f;          // constraint equation
    auto DT = D.transpose(1, 0);
    float lr = 1., chisqrd = 1e6;
    for (Int_t r = 0; r < 1; r++) {
        auto VD = torch::mm(D, torch::mm(V, DT));
        VD = torch::abs(torch::inverse(VD));
        auto delta_alpha = alpha - alpha0;
        auto lambda = torch::mm(VD, torch::mm(D, delta_alpha)) + torch::mm(VD, d);
        auto lambdaT = lambda.transpose(1, 0);
        if (r > 0)
            alpha0 = alpha;
        alpha = alpha0 - lr * torch::mm(V, torch::mm(DT, lambda));
        auto chi2 = torch::mm(lambdaT, d);
        // update covariance matrix
        V = V - lr * torch::mm(V, torch::mm(DT, torch::mm(VD, torch::mm(D, V))));
        for (int p = 0; p < 1; p++) {
            chisqrd = chi2.narrow(0, p, 1).item<float>();
        }

        if (fMassConstraint) {
            eval(y, D, f, fN, fM, fIMass, fMMass, fInit, "Mass");
        }
        if (f2CConstraint) {
            eval(y, D, f, fN, fM, fIMass, fMMass, fInit, "2C");
        }
        if (f3CConstraint) {
            eval(y, D, f, fN, fM, fIMass, fMMass, fInit, "3C");
        }
        if (f4CConstraint) {
            eval(y, D, f, fN, fM, fIMass, fMMass, fInit, "4C");
        }
    }

    fChi2 = chisqrd;
    fProb = TMath::Prob(fChi2, fNdf);
    // ----------------------------------------------------------------------

    // ----------------------------------------------------------------------
    // Pull distributions
    // ----------------------------------------------------------------------
    if (true) {
        for (int b = 0; b < (fN * fSize); b++) {
            fA(b, 0) = alpha[b][0].item<double>();
            double num = (y[b][0] - alpha[b][0]).item<double>();
            double dem = V0[b][b].item<double>() - V[b][b].item<double>();
            if (dem > 0) {
                fPull(b, b) = num / std::sqrt(dem);
            }

            for (int c = 0; c < (fN * fSize); c++) {
                fV(b, c) = V[b][c].item<double>();
            }
        }
    }
    // ----------------------------------------------------------------------
    return true;
}

HRefitCand HFitter::getDaughter(int val) {
    HRefitCand cand;
    double Px = (1. / fA(0 + val * fSize, 0)) * std::sin(fA(1 + val * fSize, 0)) * std::cos(fA(2 + val * fSize, 0));
    double Py = (1. / fA(0 + val * fSize, 0)) * std::sin(fA(1 + val * fSize, 0)) * std::sin(fA(2 + val * fSize, 0));
    double Pz = (1. / fA(0 + val * fSize, 0)) * std::cos(fA(1 + val * fSize, 0));
    double M = fM[val];
    cand.SetXYZM(Px, Py, Pz, M);
    cand.setR(fA(3 + val * fSize, 0));
    cand.setZ(fA(4 + val * fSize, 0));

    // MC-Truth
    int pid = fPid[val];
    TVector3 tp = fGMomentum[val]; // true momentum
    cand.setPid(pid);
    cand.setGMomentum(tp.X(), tp.Y(), tp.Z());
    cand.setEventNr(fNr);
    // ---------------------------------------------------------------------------
    // set covariance
    // ---------------------------------------------------------------------------
    TMatrixD cov(fSize, fSize);
    cov.Zero();
    for (int e = 0; e < fSize; e++) {
        for (int r = 0; r < fSize; r++) {
            cov(e, r) = fV(e + val * fSize, r + val * fSize);
        }
    }
    cand.setCovariance(cov);
    // ---------------------------------------------------------------------------
    return cand;
}
