#include "Analysis.h"
#include "hgeomvector.h"
#include "hparticletool.h"
#include "hphysicsconstants.h"

#include <torch/script.h> // One-stop header.
#include <torch/torch.h>

// -----------------------------------------------------------------------
// PyTorch
// -----------------------------------------------------------------------
torch::Tensor predict(const HParticleCand* cand, torch::jit::script::Module* model) {
    std::vector<torch::jit::IValue> inputs;

    Float_t var1 = cand->getMomentum() * 1e-3;
    Float_t var2 = cand->getTheta() * TMath::DegToRad();
    Float_t var3 = cand->getPhi() * TMath::DegToRad();
    Float_t var4 = cand->getMdcdEdx();
    Float_t var5 = cand->getTofdEdx();
    Float_t var6 = cand->getDistanceToMetaHit() * 1e-3;
    Float_t var7 = cand->getTofRec();

    at::Tensor input = torch::tensor({var1, var2, var3, var4, var5, var6, var7}, torch::kFloat32).view({1, 7});

    inputs.push_back(input);
    inputs.push_back(input);

    auto outputs = model->forward(inputs).toTuple();
    torch::Tensor out1 = outputs->elements()[1].toTensor();
    // torch::Tensor out2 = outputs->elements()[0].toTensor(); // not needed
    return torch::exp(out1).view(-1);
}
// -----------------------------------------------------------------------

// ------------------------------------------------------------------------
void Analysis::Loop(TString infileList, Int_t nEvents) {
    TStopwatch timer;
    timer.Start();

    // ******************************************************
    torch::jit::script::Module fModule = torch::jit::load("model_cpp.pt");
    // ******************************************************

    HLoop loop(kTRUE);
    Bool_t ret = loop.addFiles(infileList);
    if (ret == 0) {
        cout << "READBACK: ERROR : cannot find inputfiles : " << infileList.Data() << endl;
        exit(1);
    }

    // select categories here
    if (!loop.setInput("-*,+HParticleCandSim,+HGeantKine")) {
        cout << "READBACK: ERROR : cannot read input !" << endl;
        exit(1);
    } // read all categories

    loop.printCategories();
    loop.printChain();

    HCategory* catParticle = loop.getCategory("HParticleCandSim");
    if (!catParticle) {
        std::cout << "No particleCat in input!" << std::endl;
        exit(1);
    }
    HCategory* catGeant = loop.getCategory("HGeantKine");
    if (!catGeant) {
        std::cout << "No kineCat in input!" << std::endl;
        exit(1);
    }

    Int_t entries = loop.getEntries();
    if (nEvents < entries && nEvents >= 0)
        entries = nEvents;

    // start of the event loop
    for (Int_t i = 1; i < nEvents; i++) {
        //----------break if last event is reached-------------
        if (loop.nextEvent(i) <= 0) {
            cout << " end recieved " << endl;
            break;
        } // last event reached
        // HTool::printProgress(i, nEvents, 1, "Analysing evt# :");

        std::cout << " *********************** Event Number " << i << " **************************************** " << std::endl;

        // for each event there are a number of tracks
        Int_t ntracks = catParticle->getEntries();
        for (Int_t j = 0; j < ntracks; j++) {
            HParticleCandSim* cand = HCategoryManager::getObject(cand, catParticle, j);
            // skip ghost tracks (only avalible for MC events)
            if (cand->isGhostTrack())
                continue;
            // select "good" tracks
            if (!cand->isFlagBit(Particle::kIsUsed))
                continue;

            // make predictions using the trained network
            auto output = predict(cand, &fModule);
            auto Pion_node = output[0].item<float>();   // pi+
            auto Proton_node = output[1].item<float>(); // proton
            auto Kaon_node = output[2].item<float>();   // K+

            std::cout << "Pion:   " << std::setw(10) << Pion_node << std::endl;
            std::cout << "Proton: " << std::setw(10) << Proton_node << std::endl;
            std::cout << "Kaon:   " << std::setw(10) << Kaon_node << std::endl;

        } // end track loop

        std::cout << " ********************************************************************************* " << std::endl;

    } // end event loop

    timer.Stop();
    timer.Print();
} // end loop method
ClassImp(Analysis)