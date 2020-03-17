void plot(histo_name="hHyperonMass", Int_t binning=0) 
{

  TFile* fdata1 = new TFile("file:s1.root");
  TFile* fdata2 = new TFile("file:s2.root");
  TFile* fdata3 = new TFile("file:s3.root");
  TFile* fdata4  = new TFile("file:chan01.root");
  TFile* fdata5  = new TFile("file:chan05.root");
  TFile* fdata6  = new TFile("file:chan12.root");
  TFile* fdata7  = new TFile("file:chan06.root");
  TFile* fdata8  = new TFile("file:chan03.root");
  TFile* fdata9  = new TFile("file:chan14.root");
  TFile* fdata10 = new TFile("file:chan11.root");
  TFile* fdata42 = new TFile("file:chan42.root");

  // hDocapK, hPvtxZ, hDocappi, hPocappi, hHDecayLength, hAngleHades, hLambdaMassHades, hDocafdpi, hAngleFD
  // hFDDecayLength, hLambdaMassFD, hLambdaMass, hPhotonEvsBeta, hMissMass1, hMissMass0, hPi0Mass, hMissEnergy0
  // hOAPhotonLambda, hHyperonMass0, hHyperonMass, hMMvsME

  TH1D *hd1 = (TH1D*)fdata1->Get(histo_name);
  TH1D *hd2 = (TH1D*)fdata2->Get(histo_name);
  TH1D *hd3 = (TH1D*)fdata3->Get(histo_name);
  TH1D *hd4  = (TH1D*)fdata4->Get(histo_name);
  TH1D *hd5  = (TH1D*)fdata5->Get(histo_name);
  TH1D *hd6  = (TH1D*)fdata6->Get(histo_name);
  TH1D *hd7  = (TH1D*)fdata7->Get(histo_name);
  TH1D *hd8  = (TH1D*)fdata8->Get(histo_name);
  TH1D *hd9  = (TH1D*)fdata9->Get(histo_name);
  TH1D *hd10 = (TH1D*)fdata10->Get(histo_name);
  TH1D *hd42 = (TH1D*)fdata42->Get(histo_name);

  //Add different backgrounds
  
  hd1->SetLineColor(38);              // Set fill color
  hd1->SetFillColor(38);
  hd1->SetLineWidth(3);              // Set fill color
  hd2->SetLineColor(30);             // Set fill color
  hd2->SetFillColor(30);
  hd2->SetLineWidth(3);              // Set fill color
  hd3->SetLineColor(46);             // Set fill color
  hd3->SetFillColor(46);
  hd3->SetLineWidth(3);              // Set fill color

  hd1->Scale(0.1640);    //0.1856
  hd2->Scale(0.0045);    //0.0239
  hd3->Scale(0.1550);   //0.2245

  //hd1->Scale(0.1856);    // old
  //hd2->Scale(0.0239);    // old
  //hd3->Scale(0.2245);    // old  
  
  hd4->Scale(21.755);    //35.6
  hd5->Scale(0.4627);    //1.935
  hd6->Scale(0.4576);    //1.935
  hd7->Scale(0.2255);    //0.2
  hd8->Scale(0.7231);    //0.5
  hd9->Scale(0.0225);    //0.19
  hd10->Scale(0.305);   //0.23
  hd42->Scale(0.0968);

  hd1->Sumw2(false);
  hd2->Sumw2(false);
  hd3->Sumw2(false);
  hd4->Sumw2(false);
  hd5->Sumw2(false);
  hd6->Sumw2(false);
  hd7->Sumw2(false);
  hd8->Sumw2(false);
  hd9->Sumw2(false);
  hd10->Sumw2(false);
  hd42->Sumw2(false);
  
  TH1D *histo = (TH1D*) hd5->Clone();
  histo->Add(hd6);
  histo->Add(hd7);
  histo->Add(hd8);
  histo->Add(hd9);
  histo->Add(hd10);
  histo->Add(hd42);

  TH1D *sig  =(TH1D*) hd1->Clone();
  sig->Add(hd2);
  sig->Add(hd3);
  sig->Add(histo);
  //histo->Add(sig);


  TCanvas *c1 = new TCanvas("c1", "c1", 0, 0, 2000, 2200);
  c1->cd();
  c1->SetGrid();
  c1->SetTickx();
  c1->SetTicky();
  gStyle->SetOptStat(0);

 if (binning){
  histo->Rebin(binning);
  hd1->Rebin(binning);
  hd2->Rebin(binning);
  hd3->Rebin(binning);
 }

  histo->Draw("e1"); histo->SetLineWidth(2); histo->SetMarkerSize(3); histo->SetMarkerStyle(34); histo->SetMarkerColor(kBlack);//histo->SetFillStyle(3003); histo->SetFillColor(kRed); histo->SetLineColor(kRed);
  hd1->Draw("same");
  hd2->Draw("same");
  hd3->Draw("same");
  //sig->Draw("colz");

  //THStack *hs =  new THStack("hs", "");
  //hs->SetMinimum(0);hs->SetMaximum(100);
  //hs->GetXaxis()->SetTitle("M_{#Lambda#gamma} [GeV/c^{2}]");//hs->GetXaxis()->SetTitleSize(0.05);
  //hs->GetYaxis()->SetTitle("dN/dM");//hs->GetYaxis()->SetTitleSize(0.05);    
  //hs->Add(histo); histo->SetLineWidth(2); histo->SetMarkerSize(3); histo->SetMarkerStyle(34); histo->SetMarkerColor(kBlack);//histo->SetFillStyle(3003); histo->SetFillColor(kRed); histo->SetLineColor(kRed);
  //hs->Add(hd1);
  //hs->Add(hd2);
  //hs->Add(hd3);
  //hs->Draw("nostack&e1");

  std::cout << "Sigma(1385)  Significance S/Sqrt(S+B) = " << hd1->Integral()/TMath::Sqrt(histo->Integral()+hd1->Integral()) << std::endl;
  std::cout << "Lambda(1405) Significance S/Sqrt(S+B) = " << hd2->Integral()/TMath::Sqrt(histo->Integral()+hd2->Integral()) << std::endl;
  std::cout << "Lambda(1520) Significance S/Sqrt(S+B) = " << hd3->Integral()/TMath::Sqrt(histo->Integral()+hd3->Integral()) << std::endl;

  std::cout << std::endl;

  std::cout << "Sigma(1385)  Purity S/S+B = " << hd1->Integral()/(histo->Integral()+hd1->Integral()) << std::endl;
  std::cout << "Lambda(1405) Purity S/S+B = " << hd2->Integral()/(histo->Integral()+hd2->Integral()) << std::endl;
  std::cout << "Lambda(1520) Purity S/S+B = " << hd3->Integral()/(histo->Integral()+hd3->Integral()) << std::endl;  

  std::cout << std::endl;

  std::cout << "Sigma(1385)   Reconstruction Efficiency = " << (hd1->GetEntries()/2.5e6)*100 << std::endl;
  std::cout << "Lambda(1405)  Reconstruction Efficiency = " << (hd2->GetEntries()/2.5e6)*100 << std::endl;
  std::cout << "Lambda(1520)  Reconstruction Efficiency = " << (hd3->GetEntries()/2.5e6)*100 << std::endl;

  std::cout << std::endl;

  std::cout << "Signal Yield: " << hd1->Integral() << "  " << hd2->Integral() << "  " << hd3->Integral() << std::endl;

  std::cout << std::endl;

  std::cout << "Background Yield: " << histo->Integral() << std::endl;

  // Here is a legend so that multiple histograms can be plotted together. 
  TLegend *legend = new TLegend(0.7, 0.7, .88, .88);
  legend->AddEntry(hd3, "#Lambda(1520) #rightarrow #Lambda #gamma", "f");
  legend->AddEntry(hd2, "#Lambda(1405) #rightarrow #Lambda #gamma", "f");
  legend->AddEntry(hd1, "#Sigma(1385) #rightarrow #Lambda #gamma", "f");
  legend->AddEntry(histo, "Background", "lep");
  //legend->AddEntry(hd5, "pp #rightarrow pK^{+}#Lambda#pi^{0}", "l");
  //legend->AddEntry(hd6, "pp #rightarrow pK^{+}#Sigma#pi^{0}", "l");
  //legend->AddEntry(hd7, "pp #rightarrow pK^{+}#Lambda#pi^{+}#pi^{-}", "l");
  //legend->AddEntry(hd8, "pp #rightarrow pK^{+}#Lambda", "l");
  //legend->AddEntry(hd9, "pp #rightarrow pK^{+}#Sigma^{0}#pi^{+}#pi^{-}", "l");
  //legend->AddEntry(hd10, "pp #rightarrow pK^{+}#Sigma^{0}","l");

  legend->Draw("same");

    TPaveText *Pt = new TPaveText(1.15, 37, 1.3, 37);
    Pt->SetTextSize(0.035);
    Pt->SetFillColor(19);
    Pt->SetTextColor(15);
    pt->SetTextAlign(12);
    Pt->AddText("HADES Simulations");  
    //Pt->Draw("same");

}
