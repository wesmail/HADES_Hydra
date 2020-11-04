void Generate(UInt_t nEvents=10000, TString outdir="./", TString outfile="output")
{
    Int_t  rootOut      = 0;    // write pluto root file
    Int_t  calcVertex   = 1;
    Int_t  asciiOut     = 1;    // write pluto ascci output for HGeant (==0 if we use HGeantOutput)

    //################################################################################
    //########################## REACTION ############################################
    //################################################################################
    makeDistributionManager();
    //=================================================================

    // Define the reaction: pp -> Sigma1385+ @ 4.5 kinetic beam energy
    PReaction *my_reaction = new PReaction("_T1=4.5", "p", "p", "n K+ Sigma1385+ [ Lambda pi+ ]", Form("%s/%s",outdir.Data(),outfile.Data()), rootOut, 0,calcVertex,asciiOut);

    my_reaction->setHGeant(0);   // set to 1, if PLUTO run from HGeant prompt
    my_reaction->Preheating(1000);
    cout << my_reaction->Loop(nEvents) << " events recorded" << endl;
    my_reaction->Print();
}
