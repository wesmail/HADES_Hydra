//
// start up a button bar with standart HGeant commands
//
// R.Holzmann  21/05/99
// modified by I. Koenig 03/05/2016 (button "Draw Event" added)
//
void menu(void) {
   TControlBar *menu = new TControlBar("vertical","HGeant Menu",20,650);
   menu->AddButton("KUIP","goKuip()","Go to KUIP command interpreter");
   menu->AddButton("Draw Event","DrawEvent()","Simulate and draw one event");
   menu->AddButton("Draw Cave","DrawCave()","Draw cave into HIGZ window");
   menu->AddButton("Do Event","TriggerEv(1)","Simulate one event");
   menu->AddButton("Do 10 Events","TriggerEv(10)","Simulate 10 events");
   menu->AddButton("Do 100 Events","TriggerEv(100)","Simulate 100 events");
   menu->AddButton(" Do 1000 Events ","TriggerEv(1000)","Simulate 1000 events");
   menu->AddButton("Clear Screen","ClearHIGZ()","Clear HIGZ window");
   menu->AddButton("Open New File","OpenNewFile()","Open new ROOT file");
   menu->AddButton("EXIT",".q","Exit HGeant");

   gROOT->SaveContext();
   menu->Show();
}

void TriggerEv(Int_t nEv) {
   char command[20];
   sprintf(command,"trigger %d \0",nEv);
   doGeant(command);
}

void DrawEvent() {
   doGeant("next");
   doGeant("dcut cave 1 0 2 10 .025 .025");
   doGeant("swit 2 6");
   doGeant("trigger 1");
}

void DrawCave() {
   doGeant("draw cave 90 0 0 10 12 0.03 0.03");
   doGeant("swit 2 3");
}

void ClearHIGZ() {
   doGeant("next");
}

void OpenNewFile() {
   doGeant("hades/newfile");
}





