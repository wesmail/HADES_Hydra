/**
 * KinFitTutorial.h
 *
 * \brief Kinematic fitting tutorial
 * \author Waleed Esmail <w.esmail@fz-juelich.de>
 * \date Aug 07, 2020
 *
 */

// A print function for matrices
template <typename T>
void Print(T const &matrix)
{
    Int_t nrows = matrix.GetNrows();
    Int_t ncols = matrix.GetNcols();

    cout << "shape(" << nrows << "," << ncols << ")" << std::endl;

    for (Int_t i = 0; i < nrows; i++)
    {
        for (Int_t j = 0; j < ncols; j++)
        {
            Double_t element = matrix(i, j);
            if (TMath::Abs(element) < 1e-10)
                element = 0.;
            if (element >= 0.)
                std::cout << " " << std::fixed << std::setw(8) << std::scientific << element << " ";
            else
                std::cout << std::fixed << std::setw(8) << std::scientific << element << " ";
        }
        std::cout << endl;
    }

    std::cout << std::endl;
}

// ----------------------------------------------------------------------------------
// Kinematic Fitter class
// ----------------------------------------------------------------------------------
class HKinFitter
{
private:
    TMatrixD y, V, fPull;
    double fMass, fTheta;
    double fChi2, fProb;
    bool fConverged;
    int fIteration;

public:
    HKinFitter(double inputs[]);
    ~HKinFitter(){};
    TMatrixD d_eval(const TMatrixD &m_iter);
    TMatrixD D_eval(const TMatrixD &m_iter);
    double getDaughter(int val = 0) { return y(val, 0); }
    double getChi2() const { return fChi2; }
    double getProb() const { return fProb; }
    double getPull(int val = 0) { return fPull(val, val); }
    bool isConverged() const { return fConverged; }
    int getIteration() const { return fIteration; }
    void setCovariance(TMatrixD &val) { V = val; }
    void setMeasurement(TMatrixD &val) { y = val; }
    bool fit();
};

HKinFitter::HKinFitter(double inputs[])
{
    y.ResizeTo(2, 1);
    V.ResizeTo(2, 2);

    y.Zero();
    V.Zero();

    fTheta = inputs[0];
    fMass = inputs[1];

    fConverged = false;
    fIteration = 0.;
}

TMatrixD HKinFitter::d_eval(const TMatrixD &m_iter)
{
    TMatrixD d(1, 1);

    d(0, 0) = 2 * m_iter(0, 0) * m_iter(1, 0) * (1 - cos(fTheta)) - fMass * fMass;
    return d;
}

TMatrixD HKinFitter::D_eval(const TMatrixD &m_iter)
{
    TMatrixD H(1, 2);

    H(0, 0) = 2 * m_iter(1, 0) * (1 - cos(fTheta));
    H(0, 1) = 2 * m_iter(0, 0) * (1 - cos(fTheta));

    return H;
}

bool HKinFitter::fit()
{
    double lr = 1;
    TMatrixD alpha0(2, 1), alpha(2, 1);
    TMatrixD A0(y), V0(V);
    alpha0 = y;
    alpha = alpha0;
    double chi2 = 1e6;
    double convFalg = 1e9;
    TMatrixD D = D_eval(alpha);
    TMatrixD d = d_eval(alpha);

    for (int q = 0; q < 1; q++)
    {
        TMatrixD DT(2, 1);
        DT.Transpose(D);
        TMatrixD VD = D * V * DT;
        VD.Invert();
        TMatrixD delta_alpha = alpha - alpha0;
        TMatrixD lambda = VD * D * delta_alpha + VD * d;
        TMatrixD lambdaT(1, 1);
        lambdaT.Transpose(lambda);
        TMatrixD neu_alpha(2, 1);
        // fitted variables
        neu_alpha = alpha - lr * V * DT * lambda;

        double chisqrd = lambdaT(0, 0) * d(0, 0);
        if (fabs(chi2 - chisqrd) < 1e-3)
        {
            // for number of iterations > 1
            fIteration = q;
            fConverged = true;
            break;
        }

        chi2 = chisqrd;

        alpha0 = alpha;
        alpha = neu_alpha;
        // improved covariance matrix
        V = V - lr * V * DT * VD * D * V;
        D = D_eval(alpha);
        d = d_eval(alpha);
    }
    y = alpha;
    fChi2 = chi2;
    fProb = TMath::Prob(chi2, 1);
    //Print (y);
    // -----------------------------------------
    // Pull
    // -----------------------------------------
    fPull.ResizeTo(2, 2);
    for (int b = 0; b < 2; b++)
        fPull(b, b) = -10000;

    if (true)
    {
        for (int b = 0; b < 2; b++)
        {
            double num = A0(b, 0) - alpha(b, 0);
            double dem = V0(b, b) - V(b, b);
            if (dem > 0)
            {
                fPull(b, b) = num / std::sqrt(dem);
            }
        }
    }

    return true;
}

// ----------------------------------------------------------------------------------

void ToyMCTutorial_1(Int_t nevts = 100000)
{
    TLorentzVector pi0;
    pi0.SetXYZM(0., 0., 500.0, 134.9766);

    TH1F *h1 = new TH1F("MassPrefit", ";M [MeV];counts [a.u.]", 100, 55, 255);
    TH1F *h2 = new TH1F("MassRefit", ";M [MeV];counts [a.u.]", 100, 55, 255);
    TH1F *h3 = new TH1F("chi2", ";#chi^{2};counts [a.u.]", 100, 0, 5);
    TH1F *h4 = new TH1F("chi2_prob", ";P(#chi^{2});counts [a.u.]", 100, 0, 1);
    TH1F *h5 = new TH1F("pull", ";pull;counts [a.u.]", 100, -5, 5);



    // (Momentum, Energy units are Gev/C, GeV)
    const double masses[] = {0.0, 0.0};
    TGenPhaseSpace event;
    event.SetDecay(pi0, 2, masses);
    // add noise to momentum
    TRandom3 *noise = new TRandom3();

    for (Int_t n = 0; n < nevts; n++)
    {
        if (n % 10000 == 0)
            cout << "Event ------------------------------------- : " << n << endl;
        Double_t weight = event.Generate();
        TLorentzVector *p1 = event.GetDecay(0);
        TLorentzVector *p2 = event.GetDecay(1);
        TLorentzVector g1, g2;

        // smear the momentum
        double mom1 = p1->E() + noise->Gaus(0, 12);
        double mom2 = p2->E() + noise->Gaus(0, 12);

        g1.SetXYZM(mom1 * sin(p1->Theta()) * cos(p1->Phi()),
                   mom1 * sin(p1->Theta()) * sin(p1->Phi()),
                   mom1 * cos(p1->Theta()), 0.0);
        g2.SetXYZM(mom2 * sin(p2->Theta()) * cos(p2->Phi()),
                   mom2 * sin(p2->Theta()) * sin(p2->Phi()),
                   mom2 * cos(p2->Theta()), 0.0);
        h1->Fill((g1 + g2).M());

        double angle = g1.Angle(g2.Vect());
        // --------------------------------------------------------------------------------------
        // Calculate Covariance Matrix
        // --------------------------------------------------------------------------------------
        TMatrixD y(2, 1), V(2, 2);

        y(0, 0) = g1.E();
        y(1, 0) = g2.E();

        V(0, 0) = pow(12, 2);
        V(1, 1) = pow(12, 2);
        // --------------------------------------------------------------------------------------
        double Inputs[] = {angle, 134.9766};

        HKinFitter fitter(Inputs);
        fitter.setCovariance(V);
        fitter.setMeasurement(y);
        bool ok = fitter.fit();
        if (!ok)
            continue;
        h3->Fill(fitter.getChi2());
        h4->Fill(fitter.getProb());
        if (fitter.getProb() > 0.01)
            h5->Fill(fitter.getPull(0));

        TLorentzVector c1, c2;

        c1.SetXYZM(fitter.getDaughter(0) * sin(g1.Theta()) * cos(g1.Phi()),
                   fitter.getDaughter(0) * sin(g1.Theta()) * sin(g1.Phi()),
                   fitter.getDaughter(0) * cos(g1.Theta()), 0.0);
        c2.SetXYZM(fitter.getDaughter(1) * sin(g2.Theta()) * cos(g2.Phi()),
                   fitter.getDaughter(1) * sin(g2.Theta()) * sin(g2.Phi()),
                   fitter.getDaughter(1) * cos(g2.Theta()), 0.0);

        if (fitter.getProb() > 0.01)
            h2->Fill((c1 + c2).M());

    } // end event loop
    TCanvas *c0 = new TCanvas();
    c0->cd();
    h2->SetLineColor(kRed);
    h1->Draw("e1");
    h2->Draw("same&e1");

    TCanvas *c1 = new TCanvas();
    c1->cd();
    h5->Draw("e1");

    TCanvas *c2 = new TCanvas();
    c2->cd();
    h3->Draw("e1");

    TCanvas *c3 = new TCanvas();
    c3->cd();
    h4->Draw("e1");
}
