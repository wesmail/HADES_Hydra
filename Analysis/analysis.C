#include "hades.h"
#include "hloop.h"
#include "htool.h"
#include "hcategorymanager.h"
#include "hparticleanglecor.h"
#include "hparticlepairmaker.h"
#include "hparticletool.h"
#include "hphysicsconstants.h"
#include "hhistmap.h"
#include "hparticletracksorter.h"
#include "henergylosscorrpar.h"


#include "hcategory.h"
#include "hlinearcategory.h"
#include "hrichhit.h"
#include "hrichhitsim.h"
#include "hparticlecand.h"
#include "hparticlecandsim.h"
#include "hparticlepair.h"
#include "hparticlegeantpair.h"

#include "hgeantkine.h"
#include "hparticledef.h"
#include "hstartdef.h"
#include "richdef.h"


#include "hparticlegeant.h"
#include "hparticlegeantdecay.h"
#include "hparticlegeantevent.h"
#include "hparticlecutrange.h"

#include "TTree.h"

#include <vector>
#include <algorithm>
#include <map>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace Particle;

Bool_t selectHadrons(HParticleCand* pcand)
{
    // build in selection function for hadron candidates.
    // Requires besides an RK + META and fitted
    // inner+outer segment.

    Bool_t test = kFALSE;
    if( pcand->isFlagAND(4,
			 Particle::kIsAcceptedHitInnerMDC,
			 Particle::kIsAcceptedHitOuterMDC,
			 Particle::kIsAcceptedHitMETA,
			 Particle::kIsAcceptedRK
			)
       &&
       pcand->getInnerSegmentChi2() > 0
       &&
       pcand->getChi2()             < 10000      // RK
      ) test = kTRUE;

    if(!test) return kFALSE;

    if(test) test = pcand->getMetaMatchQuality() < 3 ? kTRUE : kFALSE ;

    return test;
}

Int_t analysis(TString infileList="inputDST.root", Int_t nEvents=10000){
    
    TStopwatch timer;
    timer.Start();

    // -----------------------------------------------------------------------
    // define output file and some histograms
    // -----------------------------------------------------------------------
    // set ouput file
    TFile *outfile = new TFile("example.root","recreate");

    TH2F *h01 = new TH2F("h01","",200, -2000, 2000, 200, 0, 20);
    h01->SetXTitle(" pxq [MeV/c]");
    h01->SetYTitle(" MdcdEdx [a.u.]");

    TH2F *h02 = new TH2F("h02","",200, -2000, 2000, 200, 0, 20);
    h02->SetXTitle(" pxq [MeV/c]");
    h02->SetYTitle(" TofdEdx [a.u.]");    

    TH2F *h03 = new TH2F("h03","",200, -2000, 2000, 120, 0, 1.2);
    h03->SetXTitle(" pxq [MeV/c]");
    h03->SetYTitle(" #beta ");

    TH2F *h04 = new TH2F("h04","",100, 0, 360, 100, 0, 90);
    h04->SetXTitle(" #phi");
    h04->SetYTitle(" #theta "); 

    TH1F *h05 = new TH1F("h05","",100, 1070, 1170);
    h05->SetXTitle(" M_{p#pi^{-}} [MeV/c^{2}]");
    h05->SetYTitle(" events ");          
    // -----------------------------------------------------------------------

    HLoop loop(kTRUE);
    Bool_t ret = loop.addFiles(infileList);
    if(ret == 0) {
        cout<<"READBACK: ERROR : cannot find inputfiles : "<<infileList.Data()<<endl;
        return 1;
    }

    // select categories here
    if(!loop.setInput("-*,+HParticleCandSim,+HGeantKine")) {
        cout<<"READBACK: ERROR : cannot read input !"<<endl;
	    exit(1);
    } // read all categories

    loop.printCategories();
    loop.printChain();

    HCategory *catParticle = loop.getCategory("HParticleCandSim"); 
    if (!catParticle) { std::cout<<"No particleCat in input!"<<std::endl; exit(1);}
    HCategory *catGeant    = loop.getCategory("HGeantKine");
    if (!catGeant) { std::cout<<"No kineCat in input!"<<std::endl; exit(1);}

    Int_t entries = loop.getEntries();
    if(nEvents < entries && nEvents >= 0 ) entries = nEvents;

    // start of the event loop
    for(Int_t i=1; i<nEvents; i++){
    //----------break if last event is reached-------------
    if(loop.nextEvent(i) <= 0) { cout<<" end recieved "<<endl; break; } // last event reached
    HTool::printProgress(i,nEvents,1,"Analysing evt# :");

    // for each event there are a number of tracks
    Int_t ntracks = catParticle->getEntries();

    std::vector<HParticleCandSim*> protons, pions;
    for(Int_t j=0; j<ntracks; j++){
      HParticleCandSim *cand = HCategoryManager::getObject(cand,catParticle,j);
      // skip ghost tracks (only avalible for MC events)
      if ( cand->isGhostTrack() ) continue;
      // select "good" tracks
      if (!fCand->isFlagBit(Particle::kIsUsed)) continue;

        // looking at some observables
        Float_t mom        = cand->getMomentum();
        Float_t charge     = cand->getCharge();
        Float_t beta       = cand->getBeta();
        Float_t theta      = cand->getTheta();
        Float_t phi        = cand->getPhi();
        Float_t mdcdEdx    = cand->getMdcdEdx();
        Float_t tofdEdx    = cand->getTofdEdx();

        h01->Fill(mom*charge, mdcdEdx);
        h02->Fill(mom*charge, tofdEdx);
        h03->Fill(mom*charge, beta);
        h04->Fill(phi, theta);

        // select particles based on MC info
        // proton pdg==14, pion pdg==9
        if (cand->getGeantPID()==14) protons.push_back(cand);
        else if (cand->getGeantPID()==9) pions.push_back(cand);
    } // end track loop

    // -----------------------------------------------------------------------
    // looking at Lambda invariant mass here
    // -----------------------------------------------------------------------
    for(size_t m=0; m<protons.size(); m++){
        HParticleCandSim *cand1 = protons[m];
        for(size_t n=0; n<pions.size(); n++){
            HParticleCandSim *cand2 = pions[n];

            // set mass hypothesis for protons and pions
            cand1->calc4vectorProperties(938.272);
            cand2->calc4vectorProperties(139.570);

            h05->Fill( (*cand1 + *cand2).M() );
            
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
