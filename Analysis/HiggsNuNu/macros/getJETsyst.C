#include "TH1F.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TStyle.h"
#include "TTree.h"
#include "TCanvas.h"

#include <iostream>
#include <sstream>


int getJETsyst(){//main

  const unsigned nS = 5;
  float sys[nS] = {0.2,0.5,1,2,2.5};
  float syserr[nS] = {0.,0,0,0,0};
  std::string syslist[nS] = {"0d2","0d5","","2","2d5"};


  TCanvas *myc = new TCanvas("myc","myc",1);


  std::string type[4] = {"JESUP","JESDOWN","JERBETTER","JERWORSE"};

  TFile *finref = TFile::Open("/vols/cms/magnan/Hinvisible/RunIILT/output_lighttree_170227/MC_Powheg-VBFHtoinv-mH125.root");
  finref->cd();
  TTree *treeref = (TTree*)gDirectory->Get("LightTree");
  TH1F *j2ptref = new TH1F("j2ptref",";jet2 p_{T} (GeV)",100,40,400);

  std::string lcut = "(jet1_eta*jet2_eta<0 && abs(jet1_eta)<4.7 && abs(jet2_eta)<4.7 && dijet_M>1300 && jet1_pt>80 && dijet_deta>4.0 && jet2_pt>40 && metnomuons>200 && nloosephotons==0 && dijet_dphi<1.5 && alljetsmetnomu_mindphi>0.5 && nvetomuons==0&&nvetoelectrons==0&&nvetotaus==0)*total_weight_lepveto*weight_0b_alljets";
  myc->cd();
  gStyle->SetOptStat(1111111);
  treeref->Draw("jet2_pt>>j2ptref",lcut.c_str());
  myc->Update();
  double valref = j2ptref->Integral(0,-1);
  std::cout << " Reference val: " << j2ptref->GetEntries() << " " << j2ptref->Integral(0,-1) << std::endl;
  //return 1;

  TFile *fin[4][nS];
  TH1F *j2pt[4][nS];
  float val[4][nS];
  float err[4][nS];
  TGraphErrors *gr[4];
  for (unsigned iT(0); iT<4; ++iT){//loop on type
    for (unsigned iS(0); iS<nS; ++iS){//loop on syst
    std::ostringstream lname;
    lname << "/vols/cms/magnan/Hinvisible/RunIILT/output_lighttree_170227/" << type[iT] << syslist[iS];
    lname << "/MC_Powheg-VBFHtoinv-mH125.root";
    fin[iT][iS] = TFile::Open(lname.str().c_str());
    fin[iT][iS]->cd();
    TTree *tree = (TTree*)gDirectory->Get("LightTree");
    lname.str("");
    lname << "j2pt_" << type[iT] << "_" << syslist[iS];
    j2pt[iT][iS] = (TH1F*)j2ptref->Clone(lname.str().c_str());
    tree->Draw(("jet2_pt>>"+lname.str()).c_str(),lcut.c_str());
    double tmp = 0;
    val[iT][iS] = j2pt[iT][iS]->IntegralAndError(0,-1,tmp);
    val[iT][iS] = (val[iT][iS]-valref)/valref;
    err[iT][iS] = tmp/valref;
    std::cout << type[iT] << " " << syslist[iS] << " " << val[iT][iS] << "+/-" << err[iT][iS] << std::endl;

    }//loop on syst
  }
  
  myc->Clear();
  for (unsigned iT(0); iT<4; ++iT){//loop on type
    gr[iT] = new TGraphErrors(nS,sys,val[iT],syserr,err[iT]);
    gr[iT]->SetTitle(";nSigma;Variation");
    gr[iT]->SetMinimum(-0.5);
    gr[iT]->SetMaximum(0.5);
    gr[iT]->SetMarkerStyle(20+iT);
    gr[iT]->SetMarkerColor(1+iT);
    gr[iT]->SetLineColor(1+iT);
    gr[iT]->Draw(iT==0?"APL":"PL");

  }//loop on type
  myc->Update();
  myc->Print("effectOfJESJER.png");


  return 0;
}//main
