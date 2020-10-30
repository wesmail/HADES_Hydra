#include "analysisSimDst.h"
#include <math.h>
using namespace std;

Bool_t analysisSimDst(TString inputDir, const char *inpFile, Int_t nEvents = 0, Int_t startEvt = 0)
{
  TString inputFile(inpFile);
  Hades *myHades = NULL;
  TStopwatch timer;
  Int_t evN = 0;

  myHades = new Hades;
  gHades->setTreeBufferSize(2000000000);
  printf("Setting configuration...\n");
  //---------------  Set batch (needed for TCanvas's) ------------------
  gROOT->SetBatch();
  //-----------------------------------------------------------------

  TString outputDir = "./out/";
  TString outputDirNtuples = "./out/";
  TString outDirQA = "./out/";
  TString paramSource = "ORACLE, ASCII";

  TString inFile = inputDir + inputFile;
  TString outputFile = outputDir + inputFile + "_dst.root";
  TString ntupleFile = outputDirNtuples + inputFile + "_ntuple.root";
  TString outPidNtFile = outputDirNtuples + inputFile + "_PIDntuple.root";
  Bool_t kDumpParam = kFALSE;

  TString asciiParFile1 = "/lustre/hebe/hades/user/iciepal/pNb/dst/ascii_par.txt";

  cout << "----------------------------------------------------" << endl;
  cout << "input file: " << inFile.Data() << endl;
  cout << "----------------------------------------------------" << endl;

  //---------------  Reference run id -----------------------------

  Int_t refId = 9502;

  if (refId == 0)
    Warning("reference run for parameters not set properly", "");

  Cat_t notPersistentCat[] =
      {
          catRichCal, catShowerCal, catMdcCal1, catMdcCal2, catTofCal,
          catTofinoCal, catShowerHit, catShowerPID,
          catMdcGeantCell, catShowerGeantWire, catRichGeantRaw,
          catRichGeantRaw + 1, catRichGeantRaw + 2, catMatchUDiLeptons,
          catMatchULeptons, catRichPID, catShowerPID,
          catMdcTrkCand, catMdcHit, catStartHit, catRichHit, catRichHitHdr,
          catRichDirClus, catRichTrack, catHardRichHit, catMUEMUDiLeptons, catMUEMULeptons,
          catTofCluster, catTofRaw, catTofCluster, catTofHit, catShowerHitHdr,
          catShowerHitTof, catShowerHitTrack, catShowerRawMatr, catShowerHitTofTrack,
          catShowerTrack, catMdcHit, catMdcSeg, catRKTrackB, catWallRaw, catWallOneHit, catWallCal, catWallGeantRaw,
          catSplineTrack, catKickTrack, catKickTrackB, catMetaMatch, catRichHitFit, catMdcClusInf};

  //----------------- Detector setup configuration ---------------------

  Int_t richMods[] = {1};

  Int_t mdcMods[6][4] = {{1, 1, 1, 1},
                         {1, 1, 1, 1},
                         {1, 1, 0, 1},
                         {1, 1, 1, 1},
                         {1, 1, 1, 1},
                         {1, 1, 1, 1}};

  Int_t tofMods[22] = {1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  Int_t showerMods[3] = {1, 1, 1};
  Int_t tofinoMods[6][1] = {{1}, {1}, {1}, {1}, {1}, {1}};
  Int_t startMods[6] = {1, 1, 1, 0, 0, 0};
  Int_t trigMods[] = {1};
  Int_t wallMods[] = {1}; //w.esmail

  // ------------ Set input data file: NO NEED TO CHANGE  --------------
  HRootSource *source = new HRootSource;
  //source->replaceHeaderVersion(0,kTRUE);
  source->setDirectory(((Text_t *)inputDir.Data()));
  source->addFile((Text_t *)inFile.Data());
  source->setGlobalRefId(refId);
  gHades->setDataSource(source);

  // ----------Add detectors to the setup: NO NEED TO CHANGE -----------

  HSpectrometer *spec = gHades->getSetup();
  //spec->addDetector(new HStartDetector); // w.esmail
  spec->addDetector(new HRichDetector);
  spec->addDetector(new HMdcDetector);
  spec->addDetector(new HTofDetector);
  spec->addDetector(new HTofinoDetector);
  spec->addDetector(new HShowerDetector);
  spec->addDetector(new HTBoxDetector);
  //spec->addDetector(new HTriggerDetector);
  spec->addDetector(new HWallDetector); //w.esmail

  // ----- Set active modules for each detector: NO NEED TO CHANGE -----

  spec->getDetector("TBox")->setModules(0, trigMods);
  spec->getDetector("Wall")->setModules(0, wallMods); //w.esmail
  for (Int_t is = 0; is < 6; is++)
  {
    spec->getDetector("Shower")->setModules(is, showerMods);
    spec->getDetector("Tof")->setModules(is, tofMods);
    spec->getDetector("Tofino")->setModules(is, tofinoMods[is]);
    spec->getDetector("Mdc")->setModules(is, mdcMods[is]);
    spec->getDetector("Rich")->setModules(is, richMods);
  }
  // -------------  RuntimeDb input: CHANGE PARAM CONTEXT HERE ----------

  HRuntimeDb *rtdb = gHades->getRuntimeDb();
  if (paramSource.Contains("ASCII"))
  {
    HParAsciiFileIo *input2 = new HParAsciiFileIo();
    input2->open((Text_t *)asciiParFile1.Data(), "in");
    rtdb->setSecondInput(input2);
  }

  // ----- SECOND PARAM INPUT FOR PID
  if (paramSource.Contains("ORACLE"))
  {
    HParOraIo *ora = new HParOraIo;
    ora->open();
    ora->setHistoryDate("now");
    rtdb->setFirstInput(ora);
  }

  // -------------  RuntimeDb output: save the containers init from Oracle ----
  if (paramSource.Contains("ORACLE") && kDumpParam)
  {
    TString outputName("SIMDST_par_");
    TDatime dt;
    dt.Set();
    outputName = outputName + "refId.";
    outputName += refId;
    outputName += "+";
    outputName += dt.GetDate();
    outputName += dt.GetTime();
    outputName += ".root";

    HParRootFileIo *output = new HParRootFileIo;
    TString outputDir = "./";
    output->open(((Text_t *)(outputDir + outputName).Data()), "RECREATE");
    rtdb->setOutput(output);
  }

  // ----------- Build TASK SETS (using H***TaskSet::make) -------------
  //HStartTaskSet *startTaskSet = new HStartTaskSet(); //w.esmail
  //HTriggerTaskSet      *triggerTaskSet      = new HTriggerTaskSet();
  HRichTaskSet *richTaskSet = new HRichTaskSet();
  //HRichIPUTaskSet      *richIPUTaskSet      = new HRichIPUTaskSet();
  HShowerTaskSet *showerTaskSet = new HShowerTaskSet();
  HTofTaskSet *tofTaskSet = new HTofTaskSet();
  HTofinoTaskSet *tofinoTaskSet = new HTofinoTaskSet();
  HShowerTofinoTaskSet *showerTofinoTaskSet = new HShowerTofinoTaskSet();
  HMdcTaskSet *mdcTaskSet = new HMdcTaskSet();
  //HMdcTrackDSet::setMixCuts(1,2);// <------------ required to compensate for missing motherboards
  //HKickTaskSet         *kickTaskSet         = new HKickTaskSet();
  HMdcLookUpTb::setUseFloatLevel();
  HMdc34ClFinder::setUseFloatLevel();
  //--------------------------------------------------------------------
  //
  HHodoTaskSet *hodoTaskSet = new HHodoTaskSet();
  HWallTaskSet *wallTaskSet = new HWallTaskSet(); //w.esmail

  //----------- TRIGGER --------------------------------------------------
  //HTask *triggerTasks       = triggerTaskSet     ->make("simulation");
  //----------- START --------------------------------------------------
  //HTask *startTasks = startTaskSet->make("simulation"); // w.esmail
  //  ----------- RICH -----------------------------------------------------
  HTask *richTasks = richTaskSet->make("simulation", "noiseon");
  //HTask *richIPUTasks       = richIPUTaskSet     ->make("simulation");
  //  ----------- SHOWER ---------------------------------------------------
  // HTask *showerTasks        = showerTaskSet      ->make("","simulation,leprecognition");
  HTask *showerTasks = showerTaskSet->make("", "simulation,lowshowerefficiency");
  //HTask *showerTasks        =showerTaskSet->make("","simulation");
  //  ----------- TOF ------------------------------------------------------
  HTask *tofTasks = tofTaskSet->make("simulation");
  HTask *wallTasks = wallTaskSet->make("simulation"); // w.esmail
  //  ----------- TOFINO ---------------------------------------------------
  HTask *tofinoTasks = tofinoTaskSet->make("", "simulation");
  //  ----------- SHOWERTOFINO ---------------------------------------------
  //HTask *showerTofinoTasks  = showerTofinoTaskSet->make("","simulation,leprecognition");
  HTask *showerTofinoTasks = showerTofinoTaskSet->make("", "simulation");

  HTask *mdcTasks = mdcTaskSet->make("rtdb", "");
  // ***** make MDC time resolution 4 x times worse *****
  Float_t scaleErr = 4;
  mdcTaskSet->getDigitizer()->setScalerTime1Err(scaleErr, scaleErr, scaleErr, scaleErr);

  // LVL1 filter
  // HTaskSet *evtfilter = new HTaskSet("evtfiltertasks","evtfiltertasks");
  //----------------SPLINETACKING----------------------------------------
  HSplineTaskSet *splineTaskSet = new HSplineTaskSet("", "");
  HTask *splineTasks = splineTaskSet->make("", "spline,reslowb,&old,tofclust,simulation,runge");

  //---------------------------- PID --------------------------------------
  TString opt_pidrec = "pdf,ALG_RUNGEKUTTA";
  HPidReconstructor *pPidRec = new HPidReconstructor(opt_pidrec.Data());

  Short_t nParticles[1] = {14}; // only lepton ID
  pPidRec->setParticleIds(nParticles, sizeof(nParticles) / sizeof(Short_t));

  TString cuts = "RICHRINGCUTS,TOFPVSBETA,TOFINOPVSBETA,SHOWERSUMDIFFVSP";
  HPidAlgStandCuts *pPidAlgStandCuts = new HPidAlgStandCuts(outPidNtFile.Data(), cuts.Data());

  pPidRec->addAlgorithm(pPidAlgStandCuts);

  HPidParticleFiller *pPartFiller = new HPidParticleFiller("RUNGEKUTTA");
  pPartFiller->setAlgorithm(7); // 7-check values from HPidAlgStandCuts algorithm
  pPartFiller->print();

  HTask *pTrackFiller = new HPidTrackFiller("makesimcategory,NOCHI2SEG1");

  //---------------------------- Project PidTrackCand to Ntuple ------------------------------------
  //HTask  *pPdfMaker = new HPidPdfMaker((Char_t *)ntupleFile.Data(),kTRUE);

  //------------------------ PAIRS -----------------------------

  HPairFiller *pPairFiller = new HPairFiller("HPairFiller", "HPairFiller");
  HPairFilter *pPairFilter = new HPairFilter("HPairFilter", "HPairFilter");
  //pPairFilter->setQAFileName((Text_t*)outPairFiltNtFile.Data());

  //-----------------------------------------------------------------------

  HQAMaker *qaMaker = new HQAMaker("qa.maker", "qa.maker");
  qaMaker->setOutputDir((Char_t *)outDirQA.Data());
  qaMaker->setPSFileName((Char_t *)inputFile.Data()); // very important
  qaMaker->setSamplingRate(1);
  qaMaker->setIntervalSize(5000);

  // WARNING - HERE ************** LVL1 FILTER ********************
  //HLvl1EvtFilter *evtflt = new HLvl1EvtFilter("eventfilter","eventfilter","sim,metamult,tofinomult,opsec",2,1);
  HLvl1EvtFilter *evtflt = new HLvl1EvtFilter("eventfilter", "eventfilter", "sim,metamult", 3);
  // evtfilter->connect(evtflt);

  //------------------------ Master task set --------------------------

  HTaskSet *masterTaskSet = gHades->getTaskSet("simulation");
  //masterTaskSet->add(startTasks); //W.Esmail
  masterTaskSet->add(tofTasks);
  masterTaskSet->add(tofinoTasks);
  masterTaskSet->add(richTasks);
  //masterTaskSet->add(richIPUTasks);
  masterTaskSet->add(mdcTasks);
  masterTaskSet->add(showerTasks);
  //  masterTaskSet->add(triggerTasks);
  //
  // WARNING - HERE ************** LVL1 FILTER ********************
  masterTaskSet->add(evtflt);

  masterTaskSet->add(splineTasks);
  masterTaskSet->add(pTrackFiller);
  //masterTaskSet->add(pidtofrectask); ////start time recalibration W.Esmail
  masterTaskSet->add(pPidRec);
  masterTaskSet->add(wallTasks); //after PID to get reconstructed START TIME //w.esmail
  masterTaskSet->add(pPartFiller);
  ////masterTaskSet->add(pPairFiller);
  ////masterTaskSet->add(pPairFilter);
  //masterTaskSet->add(qaMaker);
  //	masterTaskSet->add(pPdfMaker);
  masterTaskSet->isTimed(kTRUE);
  gHades->makeCounter(1000); //TU! bylo 100

  //------------------------ Initialization ----------------------------
  /*
  ((HMUEmulationSim*)(triggerTasks->getTask("trigger.emu")))->setParThresholds(riputhr2);
  ((HMUEmulationExp*)(triggerTasks->getTask("trigger.emu")))->setParLocMax(ripuloc2);
*/

  if (!gHades->init())
  {
    printf("Error in initialization, exiting\n");
    delete myHades;
    return EXIT_FAILURE;
  }

  //----------------- Set non-persistent categories --------------------

  HEvent *event = gHades->getCurrentEvent();
  for (UInt_t i = 0; i < sizeof(notPersistentCat) / sizeof(Cat_t); i++)
  {
    HCategory *cat = event->getCategory(notPersistentCat[i]);
    if (cat)
      cat->setPersistency(kFALSE);
  }

  //-------------------------- Set output ------------------------------

  gHades->setOutputFile((Char_t *)outputFile.Data(), "RECREATE", "Test", 2);
  gHades->makeTree();
  //------------------------------------------------------------------
  //                          EVENT LOOP
  printf("Processing events...\n");
  timer.Reset();
  timer.Start();
  if (nEvents < 1)
  {
    evN = gHades->eventLoop();
  }
  else
  {
    evN = gHades->eventLoop(nEvents, startEvt);
  }

  timer.Stop();
  //------------------------------------------------------------------
  if (paramSource.Contains("ORACLE") && kDumpParam)
  {
    if (!rtdb->initContainers(refId))
      Error("analysisSimDST.C", "RTDB not initialized for ref run");
    if (!rtdb->initContainers(refId))
      Error("analysisSimDST.C", "RTDB not initialized for ref run");

    rtdb->writeContainers();
    rtdb->setContainersStatic();
    rtdb->saveOutput();
    rtdb->print();
  }
  //------------------------------------------------------------------

  rtdb->saveOutput();
  rtdb->print();
  myHades->getTaskSet("simulation")->printTimer();

  //------------------------------------------------------------------
  delete myHades;
  //------------------------------------------------------------------
  printf("------------------------------------------------------\n");
  printf("Events processed\t: %i\n", evN);
  printf("Real time\t\t: %f\n", timer.RealTime());
  printf("Cpu time\t\t: %f\n", timer.CpuTime());
  if (evN)
    printf("Performance\t\t: %f cpu-s/ev\n", timer.CpuTime() / evN);
  if (evN)
    printf("Performance (real)\t: %f ev/s\n", evN / timer.RealTime());
  return EXIT_SUCCESS;
}

void usage(char *name)
{
  printf("Usage:\n");
  printf("\t%s %s %s\n", "", gSystem->BaseName(name), "");
}

int main(int argc, char **argv)
{

  TROOT AnalysisSimDst("AnalysisSimDst", "compiled analysis macros");

  if (argc < 2)
    std::cerr << "- Usage "
              << "./analysisSimDst inDir inFile nEvents " << std::endl;
  else
    return analysisSimDst(TString(argv[1]), TString(argv[2]), atoi(argv[3]));
}
