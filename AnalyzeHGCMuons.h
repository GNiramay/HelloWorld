/* This part of the comment was edited from cloned repository */




#ifndef ANALYZEHGCMuons_H
#define ANALYZEHGCMuons_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "HGCNtupleVariables.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TGraph.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TDirectory.h"
#include<TF1.h>
class AnalyzeHGCMuons : public HGCNtupleVariables{

 public:
  AnalyzeHGCMuons(const TString &inputFileList="foo.txt", const char *outFileName="histo.root",const char *dataset="data");
  ~AnalyzeHGCMuons();
  Bool_t   FillChain(TChain *chain, TChain *chain2, const TString &inputFileList);
  Long64_t LoadTree(Long64_t entry);
  void     EventLoop(const char *);
  void     BookHistogram(const char *);
  
  TFile *oFile;
  /* TDirectory *d_Layer1; */
  /* TH1F *h_ADChg[128]; */
  TH1F* h_nrechits;
  TH1F* h_nrechit_1evt[5];	// Tracking hits in a single event
  void Assignment1(int jentry);		// This will plot everything related to Layer 15 problem
  //////// Assignment 1
  /* TH1F* h_layered_nrechits[28];	// To plot nrechits for each layer */
  /* TH1F* h_Debug15;		// To seen problem with layer no 15 */
  /* TH2F* h_Debug15geom;		// To see geometrical distribution of rechits for layer 15 */
  /* TH2F* h_Chip3_rechits;	// Layer 15 has dip here */
  ////////////// Shower Developement/////////////
  TH1F* h_PerEventShower[3];	// Shower profile of 3 events
  TH3F* h_ShowerDev[3];		// Distinguishing early and late showers
  TH1F* h_ShowerX[3];
  TH1F* h_ShowerY[3];
  TH3F* h_Electron;		// Shower developement of electron candidates
  TH3F* h_Pion;			// Shower developement of pion candidates
  TH1F* h_ratio;
  TH2F* h_RatioCheck;		// histogram to check e-pi distinction for different n1's
  TH2F* h_CumuNRechits;		// histogram of cumulative rechits
  TH2F* h_CumuEnergy;
  TH2F* h_NRvsE;		// relation between nrechits and actual energy

  /* Reproducing Shubham's Analysis */
  float* Reprod(int jentry); /* This will return array of avg rechits energy in each layer */
  bool My_P_Cut(int option);   /* cut on Pion. int option=0 -> no cut */
  /* TF1 *edc; */
  TH1F* h_Tot_En;		/* Total rechit energy of electrons */
  TH2F* h_ratio_En;		/* normalised rechit energy per layer */
  TH2F* h_NRvsER[28];		/* Layerwise NRechits vs energy ratio */
  TH1F* h_EnR[28];		/* Layerwise energy ratio histograms */
  TH2F* h_NRvsMolE;		// NRechits vs energy contained in 1 Moller radius
  TH3F* h_TestMR;		/* Checking Moller radius of data */
  TH2F* h_XvsNR;		/* X=alpha/beta bs NRechits */
  /////////// Chi Square Fit
  TH2F* h_ChivsNR;		// Chisquare vs NRechits
  Double_t EnDecay(Double_t*,Double_t*); /* Energy Decay function */
  float MyChiFit(float*,float*); /* For Chisquare fit of avg layerwise energy distribution */
  float layers[28];
  /////////// My E-Pi differentiation
  void My_E_Pi_Diff();		/* This function will perform all the plotting tasks required */
  TH1F* h_TMaxE;			/* Matching TMax with ln(E/Ec)-1*/
  TH1F* h_TMaxN;			/* rechits wise shower maxima */
  TH2F* h_NECorr;		/* To check correlation between above two qtys */
  TH2F* h_TransProf;		/* Trnsverse Profile */
  TH2F* h_TEvsE;		/* Energy maximum vs total energy detected */
  TH2F* h_TNvsE;		/* nrechit maximum vs total energy detected */
  TH1F* h_ChiFit;		/* Chi2/NDF fit of the shower profile */
  TH2F* h_radial;		/* radial shower profile */
  TH1F* h_COG_E;		/* Center of Gravity w.r.t. Energy */
  TGraph* g_shower;		/* shower profile of single event */
  TH2F* h_CG_NR;		/* 2D Histogram of COG vs nrechits */
  /* TCanvas* tc_temp;		/\* temp tc for above tgraph *\/ */
  /* const float Ec=0.5;			/\* Critical energy *\/ */
};
#endif

#ifdef ANALYZEHGCMuons_cxx

void AnalyzeHGCMuons::BookHistogram(const char *outFileName) {

//  char hname[200], htit[200];
//  double xlow = 0.0,  xhigh = 2000.0;
//  int nbins = 2000;

  oFile = new TFile(outFileName, "recreate");
  /* d_Layer1 = oFile->mkdir("Layer1"); */
  /* h_nrechits = new TH1F("h_nrechits","NRechits",500,0,1000); */
  /* char* hname= new char[200]; */
  /* d_Layer1-> cd(); */
  /* for(int ii=0; ii<128; ii++) { */
  /*   sprintf(hname, "ADC_HG_%d", ii+1); */
  /*   h_ADChg[ii] = new TH1F(hname, hname, 100, 0, 400); */
  /* } */
  ////////////// Assignment 1/////////////
  /* h_Debug15 = new TH1F("Layer15","Rechits at Layer 15 per event",10057,0,10057); */
  /* h_Debug15geom=new TH2F("Layer15_Geom","Chipwise rechits distribution",4,0,4,64,0,64); */
  /* h_Chip3_rechits=new TH2F("Chip3_rechits","Chip3 channelwise distribution",30,0,30,64,0,64); */
  /* char* HName=new char[50]; */
  /* for(int i=0;i<28;i++){ */
  /*   sprintf(HName,"RC_layer_%d",i+1); */
  /*   h_layered_nrechits[i]=new TH1F(HName,"Rechits distribution for one layer",75,0,75); */
  /* } */

  ////////////// Shower Developement/////////////
  /* h_PerEventShower[0] = new TH1F("LayerEvt10","10th event shower profile",30,0,30); */
  /* h_PerEventShower[1] = new TH1F("LayerEvt20","20th event shower profile",30,0,30); */
  /* h_PerEventShower[2] = new TH1F("LayerEvt30","30th event shower profile",30,0,30); */
  /* h_ShowerDev[0] = new TH3F("LocationEvt10","Event 10 Shower",30,0,30,20,-10,10,20,-10,10); */
  /* h_ShowerDev[1] = new TH3F("LocationEvt20","Event 20 Shower",30,0,30,20,-10,10,20,-10,10); */
  /* h_ShowerDev[2] = new TH3F("LocationEvt30","Event 30 Shower",30,0,30,20,-10,10,20,-10,10); */
  /* h_ShowerX[0] = new TH1F("ShowerXEvt10","Event 10 X",20,-10,10); */
  /* h_ShowerX[1] = new TH1F("ShowerXEvt20","Event 20 X",20,-10,10); */
  /* h_ShowerX[2] = new TH1F("ShowerXEvt30","Event 30 X",20,-10,10); */

  /* h_ShowerY[0] = new TH1F("ShowerYEvt10","Event 10 Y",20,-10,10); */
  /* h_ShowerY[1] = new TH1F("ShowerYEvt20","Event 20 Y",20,-10,10); */
  /* h_ShowerY[2] = new TH1F("ShowerYEvt30","Event 30 Y",20,-10,10); */

  /* h_Electron = new TH3F("EShower","Electron candidate shower",30,0,30,20,-10,10,20,-10,10); */
  /* h_Pion = new TH3F("PShower","Pion candidate shower",30,0,30,20,-10,10,20,-10,10); */
  /* h_ratio=new TH1F("LayerRatio","nrechits ratio first 15 layers/last 15 layers",210,-1,20); */
  /* h_RatioCheck=new TH2F("RatioCheck","initial to final layer rechits ratio",26,1,27,600,-1,5); /\* for 300,-1,2 a distinct feature was seen *\/ */
  /* h_CumuNRechits=new TH2F("CumuNrechits","ratio of cumulative n-rechits",27,1,28,100,0,1); */
  /* h_CumuEnergy=new TH2F("CumuEnergy","ratio of cumulative energy",27,1,28,100,0,1); */
  /* h_NRvsE=new TH2F("NRvsE","Energy vs NRechits relation",15000,0,1500,1500,0,1500); */

  ///////////////////////Reprod_Shubham
  /* h_Tot_En=new TH1F("EnTotal","rechits energy distribution",200,0,20000); */
  /* h_ratio_En=new TH2F("EnRatio","Layerwise Rechits energy distribution",28,0,28,500,0,1); */
  /* for(int i=0;i<28;i++){ */
  /*   char a[50],b[50],c[50],d[50]; */
  /*   /\* sprintf(a,"EnR%d",i+1); *\/ */
  /*   sprintf(a,"En%d",i+1); */
  /*   sprintf(b,"Energy ratio distribution of Layer %i",i+1); */
  /*   sprintf(c,"NRvsER%d",i+1); */
  /*   sprintf(d,"NRechits vs energy fraction in Layer %i",i+1); */
  /*   /\* h_EnR[i]=new TH1F(a,b,200,0,1); *\/ */
  /*   h_EnR[i]=new TH1F(a,b,1000,0,1000); */
  /*   h_NRvsER[i]=new TH2F(c,d,200,0,1,1000,0,1000); */
  /* } */
  /* h_NRvsMolE=new TH2F("MolEvsNR","energy fraction in 1 moller radius vs NRechits",1000,0,1000,100,0,1); */
  /* h_ChivsNR= new TH2F("ChivsNR","ChiSquare/NDF vs NRechits",2000,0,2000,100,0,1000); */
  /* h_TestMR=new TH3F("MRTesting","Choosing suitable Moller radius",1000,0,1000,100,0,1,10,1,2); */
  /* h_XvsNR=new TH2F("XvsNR","NRechits vs Alpha/Beta",1000,0,100,1000,0,1000); */

  /* h_TMaxE=new TH1F("TMaxE","Layer with maximum deposited energy",28,1,29); */
  /* h_TMaxN=new TH1F("TMaxN","Layer with maximum no. of rechits",28,1,29); */
  h_TMaxE=new TH1F("TMaxE","Layer with maximum deposited energy",28,1,29);
  h_TMaxN=new TH1F("TMaxN","Layer with maximum no. of rechits",28,1,29);
  /* h_NECorr=new TH2F("NECorr","energy maxima vs nrechits maxima",28,1,29,28,1,29); */
  h_NECorr=new TH2F("NECorr","energy maxima vs nrechits maxima",29,0,29,29,0,29);
  h_TEvsE=new TH2F("TEvsE","Energy maximum vs Total detected energy",400,2,6,50,0,50);
  h_TNvsE=new TH2F("TNvsE","NRechits maximum vs Total detected energy",400,2,6,50,0,50);
  h_TransProf=new TH2F("TransProf","rechit x-y coordinates",14,-7,7,14,-7,7);
  h_ChiFit=new TH1F("ChiFit","Chi^{2}/NDf",1000,0,20000);
  h_radial=new TH2F("radial","Radial shower profile",70,0,7,150,0,1.5);
  h_COG_E=new TH1F("COG_E","Center of Gravity of shower",150,0,60);
  h_CG_NR = new TH2F("CG_NR","Center of gravity vs nrechits",500,0,1000,150,0,60);
  /* g_shower=new TGraph(28,); */
  std::cout<<"HahaLol\n":
}


AnalyzeHGCMuons::AnalyzeHGCMuons(const TString &inputFileList, const char *outFileName, const char* dataset) {

  TChain *tree = new TChain("rechitntupler/hits");
  TChain *tree2 = new TChain("trackimpactntupler/impactPoints");

  if( ! FillChain(tree, tree2, inputFileList) ) {
    std::cerr << "Cannot get the tree " << std::endl;
  } else {
    std::cout << "Initiating analysis of dataset " << dataset << std::endl;
  }

  HGCNtupleVariables::Init(tree, tree2);

  BookHistogram(outFileName);
  
}

Bool_t AnalyzeHGCMuons::FillChain(TChain *chain, TChain *chain2, const TString &inputFileList) {

  ifstream infile(inputFileList, ifstream::in);
  std::string buffer;

  if(!infile.is_open()) {
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
    return kFALSE;
  }

  std::cout << "TreeUtilities : FillChain " << std::endl;
  while(1) {
    infile >> buffer;
    if(!infile.good()) break;
    //std::cout << "Adding tree from " << buffer.c_str() << std::endl;                                                   
    chain->Add(buffer.c_str());
    chain2->Add(buffer.c_str());
  }
  std::cout << "No. of Entries in chain  : " << chain->GetEntries() << std::endl;
  std::cout << "No. of Entries in chain2 : " << chain2->GetEntries() << std::endl;
  return kTRUE;
}

Long64_t AnalyzeHGCMuons::LoadTree(Long64_t entry) {
  // Set the environment to read one entry                                                                                          
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (!fChain->InheritsFrom(TChain::Class()))  return centry;
  TChain *chain = (TChain*)fChain;
  if (chain->GetTreeNumber() != fCurrent) {
    fCurrent = chain->GetTreeNumber();
    //    Notify();
  }

  if (!fChain2) return -5;
  Long64_t centry2 = fChain2->LoadTree(entry);
  if (centry2 < 0) return centry2;
  if (!fChain2->InheritsFrom(TChain::Class()))  return centry2;
  TChain *chain2 = (TChain*)fChain2;
  if (chain2->GetTreeNumber() != fCurrent) {
    fCurrent = chain->GetTreeNumber();
    //    Notify();
  }

  if (centry==centry2)
    return centry;
  else return -1;
}

AnalyzeHGCMuons::~AnalyzeHGCMuons() { 

  if (!fChain || !fChain2) return;
  delete fChain->GetCurrentFile();
  delete fChain2->GetCurrentFile();
  oFile->cd();
  oFile->Write();
  oFile->Close();

}

#endif
