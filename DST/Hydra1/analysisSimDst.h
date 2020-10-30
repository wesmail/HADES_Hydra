#ifndef HANALYSISSIMDST_H
#define HANALYSISSIMDST_H

#ifdef __CINT__
#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1
#endif

#ifndef __CINT__
#include <cstdio>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <TCanvas.h>
#include <TCutG.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLatex.h>
#include <TMath.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TStopwatch.h>
#include <TString.h>
#include <TStyle.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TUnixSystem.h>
#include "hades.h"
#include "haddef.h"
#include "heventheader.h"
#include "hiterator.h"
#include "hmatrixcategory.h"
#include "hrecevent.h"
#include "hpartialevent.h"
#include "htree.h"
#include "hrootsource.h"
#include "hrichdetector.h"
#include "hmdcdetector.h"
#include "htofdetector.h"
#include "htofrec.h"
#include "hpidtofrec.h"
#include "hpidevtinfofiller.h"
#include "hpidparticlefiller.h"
#include "hpidalgmomvsbeta.h"
#include "hpidalgstandcuts.h"
#include "hpidtrackcleaner.h"
#include "htofinodetector.h"
#include "hshowerdetector.h"
#include "hstartdetector.h"
#include "hlatchunpacker.h"
#include "htboxdetector.h"
#include "htriggerdetector.h"
#include "htriggertaskset.h"
#include "hrichtaskset.h"
#include "hrichIPUtaskset.h"
#include "hrichchernovringfitter.h"
#include "hkicktaskset.h"
#include "hmdctaskset.h"
#include "hmdctrackdset.h"
#include "hshowertaskset.h"
#include "hshowertofinotaskset.h"
#include "hstarttaskset.h"
#include "htoftaskset.h"
#include "htofinotaskset.h"
#include "hqamaker.h"
#include "hrichunpackerraw99.h"
#include "hrichunpackercal99.h"
#include "hmdcunpacker.h"
#include "htofunpacker.h"
#include "htofinounpacker.h"
#include "hshowerunpacker.h"
#include "hstartunpacker.h"
#include "hmatchuunpacker.h"
#include "muEmulationSim.h"
#include "muEmulationExp.h"
#include "hrichIPUremakeSim.h"
#include "hspectrometer.h"
#include "showerdef.h"
#include "hstartdef.h"
#include "tofinodef.h"
#include "tofdef.h"
#include "richdef.h"
#include "hmdcdef.h"
#include "hmdctrackddef.h"
#include "triggerinfodef.h"
#include "hmudata.h"
#include "muDilepEmulation.h"
#include "hmuleptons.h"
#include "hrootsource.h"
#include "hldfilesource.h"
#include "hmdcraw.h"
#include "hmdccal1.h"
#include "hmdccal1sim.h"
#include "hmdccal2.h"
#include "hmdccal2sim.h"
#include "hmdcclus.h"
#include "hmdchit.h"
#include "hmdchitsim.h"
#include "hmdcseg.h"
#include "hmdcsegsim.h"
#include "hmdcsetup.h"
#include "hmdcvertexfind.h"
#include "hmdcbitflipcor.h"
#include "hmdclookuptb.h"   //-newly-added--
#include "hmdc34clfinder.h" //-newly-added--
#include "hsplinetrackF.h"
#include "hsplinetaskset.h"
#include "hkicktrackbaseF.h"
#include "hmetamatchF.h"
#include "hrktrackBF.h"
#include "hparoraio.h"
#include "hparrootfileio.h"
#include "hparasciifileio.h"
#include "hruntimedb.h"
#include "htofraw.h"
#include "htofhit.h"
#include "htofhitsim.h"
#include "htofinocal.h"
#include "htofinocalsim.h"
#include "hkicktrack.h"
#include "hrtmdctrk.h"
#include "hmagnetpar.h"
#include "hmdctrackgdef.h"
#include "hsplinetaskset.h"
#include "hrichIPUparthresholds.h"
#include "hrichIPUparlocmax.h"
#include "hshowerhitfinder.h"
#include "hshowerhittoftrackmatcher.h"
#include "hshowerhittrackmatcher.h"
#include "hshowertofinocorrelator.h"
#include "hlvl1evtfilter.h"
#include "hkicktaskset.h"
#include "kickdef.h"
#include "hgeantkine.h"
#include "hgeantmdc.h"
#include "showertofinodef.h"
#include "hpidtrackfiller.h"
#include "hpidreconstructor.h"
#include "hpidpdfmaker.h"
#include "hqamaker.h"

#include "hpair.h"
#include "hpaircontfact.h"
#include "hpaircutpar.h"
#include "hpairdata.h"
#include "hpaireffpar.h"
#include "hpairevtfilter.h"
#include "hpairevtmixer.h"
#include "hpairext.h"
#include "hpairfiller.h"
#include "hpairfilter.h"
#include "hpairfl.h"
#include "hpairgeantdata.h"
#include "hpairhisto.h"
#include "hpairqa.h"
#include "hpairsim.h"
#include "hpairsimext.h"
#include "hwallcalibrater.h"
#include "hwallcalpar.h"
#include "hwallcontfact.h"
#include "hwalldetector.h"
#include "hwallhit.h"
#include "hwallhitf.h"
#include "hwalllookup.h"
#include "hwallonehit.h"
#include "hwallonehitf.h"
#include "hwallparasciifileio.h"
#include "hwallparoraio.h"
#include "hwallparrootfileio.h"
#include "hwallraw.h"
#include "hwallrefwinpar.h"
#include "hwalltaskset.h"
#include "hwalltrbunpacker.h"
#include "hwallunpacker.h"
#include "hhodotaskset.h"
#include "walldef.h"

#endif

#endif