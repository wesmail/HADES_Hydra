#define FWDET 1

//
#include "hades.h"
#include "hspectrometer.h"
#include "hdetector.h"
#include "hruntimedb.h"
#include "htask.h"
#include "hevent.h"
#include "hcategory.h"
#include "hdst.h"
#include "htime.h"
#include "hsrckeeper.h"


#include "hstarttaskset.h"
#include "hrichtaskset.h"
#include "hrich700taskset.h"
#include "hmdctaskset.h"
#include "htoftaskset.h"
#include "hrpctaskset.h"
#include "hshowertaskset.h"
#include "hemctaskset.h"
#include "hwalltaskset.h"
#include "hsplinetaskset.h"
#ifdef FWDET
#include "hfwdettaskset.h"
#include "hfwdetdetector.h"
#include "hfwdetvectorfinder.h"
#endif

#include "hparticlevertexfind.h"
#include "hparticlecandfiller.h"
#include "hparticletrackcleaner.h"
#include "hparticleevtinfofiller.h"
#include "hparticlebt.h"
#include "hparticlet0reco.h"
#include "hqamaker.h"

#include "haddef.h"
#include "richdef.h"
#include "hmdcdef.h"
#include "hmdctrackddef.h"
#include "hmdctrackgdef.h"
#include "showerdef.h"
#include "emcdef.h"
#include "rpcdef.h"
#include "tofdef.h"
#include "walldef.h"



#include "hrich700digimappar.h"
#include "hrich700ringfinderpar.h"

#include "hmdcsetup.h"
#include "hmagnetpar.h"
#include "hmdclayercorrpar.h"
#include "hmdcdedx2maker.h"
#include "hmdcdigitizer.h"
#include "hrich700digitizer.h"      // rich700
#include "hmdctrackdset.h"
#include "hmdc12fit.h"
#include "hmetamatchF2.h"
#include "hstart2hitfsim.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TString.h"
#include "TStopwatch.h"
#include "TDatime.h"


#include <iostream>
#include <cstdio>

using namespace std;





Int_t analysisDST(TString inFile, TString outdir,Int_t nEvents=1, Int_t startEvt=0)
{
    new Hades;
    TStopwatch timer;
    gHades->setTreeBufferSize(8000);
    gHades->makeCounter(10);
    gHades->setBeamTimeID(Particle::kApr12);
    gHades->getSrcKeeper()->addSourceFile("analysisDST.cc");
    gHades->getSrcKeeper()->addSourceFile("run_job.py");


    //####################################################################
    //######################## CONFIGURATION #############################
    printf("Setting configuration...+++\n");

    HRuntimeDb *rtdb = gHades -> getRuntimeDb();
    //Int_t refId = 12001; // apr12sim_mediumfieldalign_auau;
    //Int_t refId = 15000; // prob2017 ag1650ag;
    Int_t refId = 14100; // aug14;
    Bool_t HighMult = kFALSE;  // kTRUE : ag1650ag,ag4500ag , kFALSE : pp,pAg


    TString beamtime     ="apr12"; // used for spectrometer setup and track params
    TString beamtimeTrack="apr12"; // highmult = apr12, lowmult = aug14 used for  track params

    if(HighMult) {
	beamtimeTrack="apr12";
    } else {
        beamtimeTrack="aug14";
    }
    beamtime = "aug14";


    //-------------- Default Settings for the File names -----------------

    TString baseDir = outdir;
    if(!baseDir.EndsWith("/")) baseDir+="/";
    TString outDir   = baseDir;
    TString outDirQA = outDir+"qa/";

    TString outFileSuffix = "_BG.root";

    TString asciiParFile = "./params/all_pars_ag1650ag_26032018_fwdet_4.txt";
    TString rootParFile = "./allParam_AUG14_gen2_26022016_170632_refID_14100.root";
    TString paramSource = "ascii,root"; // root, ascii, oracle
    TString paramRelease = "APR12SIM_dst_gen8a"; // 27012016

    Float_t scaleDeltaElectrons=0.35; // Au+Au to Ag+Ag (47*47 / 79*79)
    Bool_t kParamFile        = kFALSE;
    Bool_t doExtendedFit     = kTRUE; // switch on/off fit for initial params of segment fitter (10 x slower!)
    Bool_t doMetaMatch       = kFALSE;  // default : kTRUE, kFALSE switch off metamatch in clusterfinder
    Bool_t useOffVertex      = kTRUE;  // default : kTRUE,  kTRUE=use off vertex procedure  (apr12 gen8:kFALSE, gen9:kTRUE)
    Bool_t doMetaMatchScale  = kTRUE;
    Bool_t useWireStat       = kTRUE;
    Float_t metaScale        = 1.5;
    Bool_t doTree            = kTRUE;
    Bool_t doMDCDeltaElectron = kTRUE;
    Bool_t doRICHDeltaElectron= kFALSE; // kFALSE : suppress inserted delta electrons
    Bool_t doReverseField     = kFALSE; // kTRUE : switch magnet polarity

    //--------------------------------------------------------------------
    // switch off unwanted categories in the poutput
    Cat_t notPersistentCat[] =
    {
	//catRich, catMdc, catShower, catTof, catTracks,
	catRichRaw,
	//catRichHitHdr,  // rich700
        catRichTrack,
	//catRichDirClus,
	//catRichHit,
	//catRichCal,
	catMdcRaw,
	catMdcCal1,   // changed
	catMdcCal2, catMdcHit,
	catMdcClusInf,
	catMdcSeg,
	catMdcTrkCand, catMdcRawEventHeader,

	catShowerCal, catShowerPID, catShowerHitHdr, catShowerRaw,
	catShowerHit,  // changed
        catShowerGeantWire,catShowerRawMatr,catShowerTrack,

        catEmcRaw,catEmcCal,
        //catEmcCluster

	catTofRaw,
	catTofHit,     // changed
	catTofCluster, // changed

	catRpcRaw,
	catRpcHit,      // changed
	catRpcCluster,  // changed
	catRpcCal,

	catRKTrackB, catSplineTrack,
	catMetaMatch,
	//catParticleCandidate, catParticleEvtInfo,
        catParticleMdc,
	catWallRaw, catWallOneHit, catWallCal,

#ifdef FWDET
//        catFwDetGeantRaw,
#endif
	catMdcGeantRaw,catTofGeantRaw,catRpcGeantRaw,catShowerGeantRaw,catEmcGeantRaw,catWallGeantRaw,
	catStartGeantRaw,catRichGeantRaw,catRichGeantRaw+1,catRichGeantRaw+2

    };
    //--------------------------------------------------------------------

    //####################################################################
    //####################################################################


    //------------- Operations on the filenames --------------------------
    TString rootSuffix =".root";
    TString nFile;    // needed to build outputfile
    TString dirname;  // needed for datasource
    TString filename; // needed for datasource
    TString outFile;

    Int_t sourcetype = 3; // root source

    if(inFile.Contains(",")){ // comma seperated list for geant merge source
	sourcetype = 4;
	inFile.ReplaceAll(" ","");
	TObjArray* ar = inFile.Tokenize(",");
	TString firstfile;
	if(ar){
	    if(ar->GetEntries()>0) {
		firstfile = ((TObjString*)ar->At(0))->GetString();
	    }
	    delete ar;
	}
	nFile     = gSystem->BaseName(firstfile.Data());
	filename  = inFile;
	dirname   = "";

    }  else {  // root source

	nFile     = gSystem->BaseName(inFile.Data());
	filename  = gSystem->BaseName(inFile.Data());
	dirname   = gSystem->DirName(inFile.Data());
    }
    if (nFile.EndsWith(rootSuffix)) nFile.ReplaceAll(rootSuffix,"");
    outFile  = outDir+nFile+outFileSuffix;
    outFile.ReplaceAll("//", "/");

    if(gSystem->AccessPathName(outDir.Data()) != 0){
	cout<<"Creating output dir :"<<outDir.Data()<<endl;
	gSystem->Exec(Form("mkdir -p %s",outDir.Data()));
    }
    if(gSystem->AccessPathName(outDirQA.Data()) != 0){
	cout<<"Creating output qadir :"<<outDirQA.Data()<<endl;
	gSystem->Exec(Form("mkdir -p %s",outDirQA.Data()));
    }
    //--------------------------------------------------------------------




    Int_t mdcMods[6][4]=
    { {1,1,1,1},
    {1,1,1,1},
    {1,1,1,1},
    {1,1,1,1},
    {1,1,1,1},
    {1,1,1,1} };

    // recommendations from Vladimir+Olga
    // according to params from 28.04.2011
    Int_t nLayers[6][4] = {
	{6,6,5,6},
	{6,6,5,6},
	{6,6,5,6},
	{6,6,5,6},
	{6,6,5,6},
	{6,6,5,6} };
    Int_t nLevel[4] = {10,50000,10,5000};

    HDst::setupSpectrometer(beamtime,mdcMods,"start,rich,mdc,tof,rpc,emc,wall,fwdet");
    // beamtime mdcMods_apr12, mdcMods_full
    // Int_t mdcset[6][4] setup mdc. If not used put NULL (default).
    // if not NULL it will overwrite settings given by beamtime
    // detectors (default)= rich,mdc,tof,rpc,shower,wall,tbox,start
#ifdef FWDET
    HSpectrometer* spec = gHades->getSetup();
    Int_t fwdetMods[]   = {1,1,0,0,1,0,1,0,0};
    spec->addDetector(new HFwDetDetector);
    spec->getDetector("FwDet")->setModules(-1,fwdetMods);
#endif


    HDst::setupParameterSources(paramSource,asciiParFile,rootParFile,paramRelease);
    //HDst::setupParameterSources("oracle",asciiParFile,rootParFile,"now"); // use to create param file
    // parsource = oracle,ascii,root (order matters)
    // if source is "ascii" a ascii param file has to provided
    // if source is "root" a root param file has to provided
    // The histDate paramter (default "now") is used wit the oracle source

    cout<<"dir="<<dirname<<", file="<<filename<<", refid="<<refId<<endl;

    HDst::setDataSource(sourcetype,dirname,filename,refId); // Int_t sourceType,TString inDir,TString inFile,Int_t refId, TString eventbuilder"
    // datasource 0 = hld, 1 = hldgrep 2 = hldremote, 3 root, 4 geantmerge
    // like "lxhadeb02.gsi.de"  needed by dataosoure = 2
    // inputDir needed by datasoure = 1,2
    // inputFile needed by datasoure = 1,3
    // for datasource 4 inputFile is a comma seprated list
    // "file1_with_path,file2_with_path,file3_with_path"


    if(kParamFile) {

	TDatime time;
        TString paramfilename= Form("allParam_APR12_sim_run_%i_gen8_%02i%02i%i",refId,time.GetDay(),time.GetMonth(),time.GetYear());  // without .root

	if(gSystem->AccessPathName(Form("%s.root",paramfilename.Data())) == 0){
	    gSystem->Exec(Form("rm -f %s.root",paramfilename.Data()));
	}
	if(gSystem->AccessPathName(Form("%s.log",paramfilename.Data())) == 0){
	    gSystem->Exec(Form("rm -f %s.log" ,paramfilename.Data()));
	}

	if (!rtdb->makeParamFile(Form("%s.root",paramfilename.Data()),"apr12sim","20-MAR-2012 00:00:00","20-MAR-2012 00:00:01")) {
	    delete gHades;
	    exit(1);
	}
    }

    //--------------------------------------------------------------------
    // ----------- Build TASK SETS (using H***TaskSet::make) -------------
    HStartTaskSet        *startTaskSet        = new HStartTaskSet();
    HRich700TaskSet      *richTaskSet         = new HRich700TaskSet();
    HRpcTaskSet          *rpcTaskSet          = new HRpcTaskSet();
    HEmcTaskSet          *emcTaskSet          = new HEmcTaskSet();
    HTofTaskSet          *tofTaskSet          = new HTofTaskSet();
    HWallTaskSet         *wallTaskSet         = new HWallTaskSet();
    HMdcTaskSet          *mdcTaskSet          = new HMdcTaskSet();
    //    mdcTaskSet->setVersionDeDx(1); // 0 = no dEdx, 1 = HMdcDeDx2
#ifdef FWDET
    HFwDetTaskSet* fwdetTaskSet = new HFwDetTaskSet("FwDetStrawTaskSet","FwDetStrawTaskSet");
#endif

    HMagnetPar* magnet = (HMagnetPar*)rtdb->getContainer("MagnetPar");
    rtdb->initContainers(refId);
    magnet->setStatic();
    magnet->setCurrent(3200); //  <=========

    if(doReverseField){
        Int_t current = magnet->getCurrent();
	magnet->setCurrent(-1*current);
    }
    magnet->printParams();

    HMdcSetup* mysetup = (HMdcSetup*)rtdb->getContainer("MdcSetup");
    rtdb->initContainers(refId);
    mysetup->setStatic();

    //mysetup->getMdcCommonSet()->setIsSimulation(1);                 // fit
    //mysetup->getMdcCommonSet()->setAnalysisLevel(1);                // fit
    //mysetup->getMdcTrackFinderSet()->setIsCoilOff(kFALSE);          // field is on
    //mysetup->getMdcTrackFinderSet()->setNLayers(nLayers[0]);
    //mysetup->getMdcTrackFinderSet()->setNLevel(nLevel);
    //mysetup->getMdc12FitSet()->setMdc12FitSet(2,1,0,kFALSE,kFALSE); // tuned fitter, seg

    HTask *startTasks         = startTaskSet       ->make("simulation");
    HTask *richTasks          = richTaskSet        ->make("simulation"); //"NORINGFINDER"
    HTask *tofTasks           = tofTaskSet         ->make("simulation");
    HTask *wallTasks          = wallTaskSet        ->make("simulation");
    HTask *rpcTasks           = rpcTaskSet         ->make("simulation");
    HTask *emcTasks           = emcTaskSet         ->make("simulation");
    HTask *mdcTasks           = mdcTaskSet         ->make("rtdb","");

#ifdef FWDET
    HTask* fwdetTasks = fwdetTaskSet->make("simulation","strawcal,rpccal,rpchitf");
#endif
    HMdcDigitizer* digi = mdcTaskSet->getDigitizer();
    if(digi){
	digi->setDeltaElectronUse(doMDCDeltaElectron,kFALSE,109,-750.,600.,20.,scaleDeltaElectrons);
	digi->setDeltaElectronMinMomCut(2.,2.,4.5,2.,2.,4.5);  // take care of glass mirrors in sec 2+5
	digi->setTimeCutUse(kTRUE);
    }
    HRich700Digitizer* richdigi = HRich700Digitizer::getDigitizer();
    if(richdigi){
	richdigi->setDeltaElectronUse(doRICHDeltaElectron,kFALSE,109,20.,0.66*scaleDeltaElectrons); // 1 - prob  0.5 Mus RICH / 1.5 mus MDC
	richdigi->setDeltaElectronMinMomCut(0.,0.,0.,0.,0.,0.);
    }


    //----------------SPLINE and RUNGE TACKING----------------------------------------
    HSplineTaskSet         *splineTaskSet       = new HSplineTaskSet("","");
    HTask *splineTasks     = splineTaskSet      ->make("","spline,runge");

    HParticleCandFiller    *pParticleCandFiller = new HParticleCandFiller   ("particlecandfiller","particlecandfiller","NORICHVERTEXCORR");
    pParticleCandFiller->setFillMdc(kFALSE); // new : testing close pair
    HParticleTrackCleaner  *pParticleCleaner    = new HParticleTrackCleaner ("particlecleaner"   ,"particlecleaner");
    HParticleVertexFind    *pParticleVertexFind = new HParticleVertexFind   ("particlevertexfind","particlevertexfind",kTRUE);
    HParticleEvtInfoFiller *pParticleEvtInfo    = new HParticleEvtInfoFiller("particleevtinfo"   ,"particleevtinfo",beamtime);
    HParticleBt            *pParticleBt         = new HParticleBt("RichBackTracking","RichBackTracking",beamtime);
#ifdef FWDET
    HFwDetVectorFinder     *pFwDetVectFind      = new HFwDetVectorFinder("FwDetVectorFinder","FwDetVectorFinder");
#endif


    //----------------------- Quality Assessment -------------------------
    HQAMaker *qaMaker =0;
    if (!outDirQA.IsNull())
    {
	qaMaker = new HQAMaker("qamaker","qamaker");
	//qaMaker->setUseSlowPar(kFALSE);
	qaMaker->setOutputDir((Text_t *)outDirQA.Data());
	//qaMaker->setPSFileName((Text_t *)hldFile.Data());
	qaMaker->setSamplingRate(1);
	qaMaker->setIntervalSize(50);
    }



    //------------------------ Master task set ---------------------------
    HTaskSet *masterTaskSet = gHades->getTaskSet("all");

    masterTaskSet->add(startTasks);
    masterTaskSet->add(tofTasks);
    masterTaskSet->add(wallTasks);
    masterTaskSet->add(rpcTasks);

    masterTaskSet->add(richTasks);

    masterTaskSet->add(emcTasks);
    masterTaskSet->add(mdcTasks);
    masterTaskSet->add(splineTasks);
#ifdef FWDET
    masterTaskSet->add(fwdetTasks);
#endif
    masterTaskSet->add(pParticleCandFiller);
    masterTaskSet->add(pParticleCleaner);
    masterTaskSet->add(pParticleVertexFind); // run after track cleaning
    masterTaskSet->add(pParticleEvtInfo);
    //masterTaskSet->add(pParticleBt);  // rich700
    masterTaskSet->add(new HParticleT0Reco("T0","T0",beamtime));
#ifdef FWDET
    masterTaskSet->add(pFwDetVectFind);
#endif

    //if (qaMaker) masterTaskSet->add(qaMaker);

    HMdcTrackDSet::setTrackParam(beamtimeTrack);
    HMdcTrackDSet::setFindOffVertTrkFlag(useOffVertex);
    if(!doMetaMatch)HMdcTrackDSet::setMetaMatchFlag(kFALSE,kFALSE);  //do not user meta match in clusterfinder
    if(doMetaMatchScale)HMetaMatchF2::setScaleCut(metaScale,metaScale,metaScale); // (tof,rpc,shower) increase matching window, but do not change normalization of MetaQA

    HStart2HitFSim* starthitf = HStart2HitFSim::getHitFinder() ;
    if(starthitf) starthitf->setResolution(0.06);    // 60 ps start res

    //--------------------------------------------------------------------
    // find best initial params for segment fit (takes long!)
    if(doExtendedFit) {
 	HMdcTrackDSet::setCalcInitialValue(1);  // -  1 : for all clusters 2 : for not fitted clusters
    }
    //--------------------------------------------------------------------

    if (!gHades->init()){
	Error("init()","Hades did not initialize ... once again");
	exit(1);
    }




    //--------------------------------------------------------------------
    //----------------- Set not persistent categories --------------------
    HEvent *event = gHades->getCurrentEvent();

    for(UInt_t i=0;i<sizeof(notPersistentCat)/sizeof(Cat_t);i++){
	HCategory *cat = ((HCategory *)event->getCategory(notPersistentCat[i]));
	if(cat)cat->setPersistency(kFALSE);
    }
    //--------------------------------------------------------------------
    if(doTree){
	// output file
	gHades->setOutputFile((Text_t*)outFile.Data(),"RECREATE","Test",2);
	gHades->makeTree();
    }
    Int_t nProcessed = gHades->eventLoop(nEvents,startEvt);
    printf("Events processed: %i\n",nProcessed);

    cout<<"--Input file      : "<<inFile  <<endl;
    cout<<"--QA directory is : "<<outDirQA<<endl;
    cout<<"--Output file is  : "<<outFile <<endl;

    printf("Real time: %f\n",timer.RealTime());
    printf("Cpu time: %f\n",timer.CpuTime());
    if (nProcessed) printf("Performance: %f s/ev\n",timer.CpuTime()/nProcessed);

    if(kParamFile) rtdb->saveOutput();

    delete gHades;
    timer.Stop();

    return 0;

}

#ifndef __CINT__
int main(int argc, char **argv)
{
    TROOT AnalysisDST("AnalysisDST","compiled analysisDST macros");


    TString nevents,startevent;
    switch (argc)
    {
    case 3:
	return analysisDST(TString(argv[1]),TString(argv[2])); // inputfile + outdir
	break;
    case 4:  // inputfile + outdir + nevents
	nevents=argv[3];

	return analysisDST(TString(argv[1]),TString(argv[2]),nevents.Atoi());
	break;
	// inputfile + nevents + startevent
    case 5: // inputfile + outdir + nevents + startevent
	nevents   =argv[3];
	startevent=argv[4];
	return analysisDST(TString(argv[1]),TString(argv[2]),nevents.Atoi(),startevent.Atoi());
	break;
    default:
	cout<<"usage : "<<argv[0]<<" inputfile outputdir [nevents] [startevent]"<<endl;
	return 0;
    }
}
#endif
