#include "hades.h"
#include "hcategorymanager.h"
#include "henergylosscorrpar.h"
#include "hhistmap.h"
#include "hloop.h"
#include "hparticleanglecor.h"
#include "hparticlepairmaker.h"
#include "hparticletool.h"
#include "hparticletracksorter.h"
#include "hphysicsconstants.h"
#include "htool.h"

#include "hcategory.h"
#include "hlinearcategory.h"
#include "hparticlecand.h"
#include "hparticlecandsim.h"
#include "hparticlegeantpair.h"
#include "hparticlepair.h"
#include "hrichhit.h"
#include "hrichhitsim.h"

#include "hgeantkine.h"
#include "hparticledef.h"
#include "hstartdef.h"
#include "richdef.h"

#include "hparticlecutrange.h"
#include "hparticlegeant.h"
#include "hparticlegeantdecay.h"
#include "hparticlegeantevent.h"

#include "TTree.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>

#include "hfitter.h"

using namespace std;
using namespace Particle;

void FillData(HParticleCand* cand, HRefitCand& outcand, double arr[],
              double mass)
{
    double deg2rad = TMath::DegToRad();

    TMatrixD cov(5, 5);
    cov(0, 0) = std::pow(arr[0], 2);
    cov(1, 1) = std::pow(arr[1], 2);
    cov(2, 2) = std::pow(arr[2], 2);
    cov(3, 3) = std::pow(arr[3], 2);
    cov(4, 4) = std::pow(arr[4], 2);

    outcand.SetXYZM(cand->getMomentum() * std::sin(cand->getTheta() * deg2rad) *
                        std::cos(cand->getPhi() * deg2rad),
                    cand->getMomentum() * std::sin(cand->getTheta() * deg2rad) *
                        std::sin(cand->getPhi() * deg2rad),
                    cand->getMomentum() * std::cos(cand->getTheta() * deg2rad),
                    mass);
    outcand.setR(cand->getR());
    outcand.setZ(cand->getZ());
    outcand.setCovariance(cov);
}

Bool_t selectHadrons(HParticleCand* pcand)
{
    // build in selection function for hadron candidates.
    // Requires besides an RK + META and fitted
    // inner+outer segment.

    Bool_t test = kFALSE;
    if (pcand->isFlagAND(4, Particle::kIsAcceptedHitInnerMDC,
                         Particle::kIsAcceptedHitOuterMDC,
                         Particle::kIsAcceptedHitMETA,
                         Particle::kIsAcceptedRK) &&
        pcand->getInnerSegmentChi2() > 0 && pcand->getChi2() < 10000 // RK
    )
        test = kTRUE;

    if (!test) return kFALSE;

    if (test) test = pcand->getMetaMatchQuality() < 3 ? kTRUE : kFALSE;

    return test;
}

Int_t analysis(TString infileList = "inputDST.root", Int_t nEvents = 10000)
{

    TStopwatch timer;
    timer.Start();

    // -----------------------------------------------------------------------
    // define output file and some histograms
    // -----------------------------------------------------------------------
    // set ouput file
    TFile* outfile = new TFile("kinfit_example.root", "recreate");

    TH1F* h01 = new TH1F("hLambdaMassPreFit", "", 100, 1070, 1170);
    h01->SetXTitle(" M_{p#pi^{-}} [MeV/c^{2}]");
    h01->SetYTitle(" events ");

    TH1F* h02 = new TH1F("hChi2", "", 100, 0, 10);
    h02->SetXTitle("#chi^{2}");
    h02->SetYTitle(" counts ");

    TH1F* h03 = new TH1F("hPChi2", "", 100, 0, 1);
    h03->SetXTitle("P(#chi^{2})");
    h03->SetYTitle(" counts ");

    TH1F* h04 = new TH1F("hLambdaMassPostFit", "", 100, 1070, 1170);
    h04->SetXTitle(" M_{p#pi^{-}} [MeV/c^{2}]");
    h04->SetYTitle(" events ");

    TH1F* h05 = new TH1F("hPull", "", 100, -5, 5);
    h05->SetXTitle("Pull(1/P_{p})");
    h05->SetYTitle(" counts ");
    // -----------------------------------------------------------------------

    HLoop loop(kTRUE);
    Bool_t ret = loop.addFiles(infileList);
    if (ret == 0)
    {
        cout << "READBACK: ERROR : cannot find inputfiles : "
             << infileList.Data() << endl;
        return 1;
    }

    // select categories here
    if (!loop.setInput("-*,+HParticleCandSim,+HGeantKine"))
    {
        cout << "READBACK: ERROR : cannot read input !" << endl;
        exit(1);
    } // read all categories

    loop.printCategories();
    loop.printChain();

    HCategory* catParticle = loop.getCategory("HParticleCandSim");
    if (!catParticle)
    {
        std::cout << "No particleCat in input!" << std::endl;
        exit(1);
    }
    HCategory* catGeant = loop.getCategory("HGeantKine");
    if (!catGeant)
    {
        std::cout << "No kineCat in input!" << std::endl;
        exit(1);
    }

    Int_t entries = loop.getEntries();
    if (nEvents < entries && nEvents >= 0) entries = nEvents;

    // start of the event loop
    for (Int_t i = 1; i < nEvents; i++)
    {
        //----------break if last event is reached-------------
        if (loop.nextEvent(i) <= 0)
        {
            cout << " end recieved " << endl;
            break;
        } // last event reached
        HTool::printProgress(i, nEvents, 1, "Analysing evt# :");

        // for each event there are a number of tracks
        Int_t ntracks = catParticle->getEntries();

        std::vector<HRefitCand> protons, pions;
        for (Int_t j = 0; j < ntracks; j++)
        {
            HParticleCandSim* cand =
                HCategoryManager::getObject(cand, catParticle, j);
            // skip ghost tracks (only avalible for MC events)
            if (cand->isGhostTrack()) continue;
            // select "good" tracks
            if (!cand->isFlagBit(Particle::kIsUsed)) continue;

            HRefitCand candidate(cand);
            // select particles based on MC info
            // proton pdg==14, pion pdg==9
            // error values obtained from resoultion plots
            if (cand->getGeantPID() == 14)
            {
                double errors[] = {1.469 * 1e-5, 2.410 * 1e-3, 5.895 * 1e-3,
                                   1.188, 2.652};
                FillData(cand, candidate, errors, 938.272);
                protons.push_back(candidate);
            }
            else if (cand->getGeantPID() == 9)
            {
                double errors[] = {5.959 * 1e-5, 9.316 * 1e-3, 1.991 * 1e-2,
                                   4.006, 7.629};
                FillData(cand, candidate, errors, 139.570);
                pions.push_back(candidate);
            }
            else
                continue;
        } // end track loop

        // -----------------------------------------------------------------------
        // looking at Lambda invariant mass here
        // -----------------------------------------------------------------------
        for (size_t n = 0; n < protons.size(); n++)
        {
            HRefitCand cand1 = protons[n];
            for (size_t m = 0; m < pions.size(); m++)
            {
                HRefitCand cand2 = pions[m];
                // mass prefit
                h01->Fill((cand1 + cand2).M());

                // ---------------------------------------------------------------------------------
                // begin kinfit here
                // ---------------------------------------------------------------------------------
                std::vector<HRefitCand> cands;
                cands.push_back(cand1);
                cands.push_back(cand2);

                HFitter fitter(cands);
                fitter.addMassConstraint(1115.68);
                fitter.fit();

                // get fitted objects fittedcand1 and fittedcand2
                HRefitCand fcand1 = fitter.getDaughter(0); // proton
                HRefitCand fcand2 = fitter.getDaughter(1); // pion

                h02->Fill(fitter.getChi2());
                h03->Fill(fitter.getProb());
                h04->Fill((fcand1 + fcand2).M());

                // get Pull example (1/P for the fitted proton)
                h05->Fill(fitter.getPull(0));
                // ---------------------------------------------------------------------------------
            }
        }
        // -----------------------------------------------------------------------

    } // end of the events loop

    // write histograms to the output file
    outfile->cd();
    h01->Write();
    h02->Write();
    h03->Write();
    h04->Write();
    h05->Write();
    outfile->Close();

    return 0;
} // end of the macro
