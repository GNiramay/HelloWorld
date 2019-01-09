#define ANALYZEHGCMuons_cxx

#include <iostream>
#include <vector>
#include <cstring>
#include "AnalyzeHGCMuons.h"
#include <TProfile.h>
#include<TCanvas.h>
#include<TF1.h>

using namespace std;



// chip 3022,44,3028




int main(int argc, char* argv[])
{

  if (argc < 2) {
    cerr << "Please give 3 arguments " << "runList " << " " << "outputFileName" << " " << "dataset" << endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *data          = argv[3];

  AnalyzeHGCMuons hgcmuons(inputFileList, outFileName, data);
  cout << "dataset " << data << " " << endl;

  hgcmuons.EventLoop(data);
  return 0;
}

void AnalyzeHGCMuons::EventLoop(const char *data) {
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  cout << "nentries " << nentries << endl;
  cout << "Analyzing dataset " << data << " " << endl;

  Long64_t nbytes = 0, nb = 0;
  Long64_t nbytes2 = 0, nb2 = 0;
  int decade = 0;

  // ======== For Chisquare fit
  // float* layer=new float[28];
  // for(float i=0;i<28;i++){layer[(int)i]=i;}
  // AnalyzeHGCMuons::edc=new TF1("ENDECAY",[](double*x,double*p){return(p[0]*pow(x[0],p[1])*TMath::Exp(-p[2]*x[0]));},0,30,3);
  // AnalyzeHGCMuons::edc->SetParameters(10,2,0.1);
  // ========

  // float* AvgEn = new float[28];	// For reproducing Shubham's data
  // for(int i=0;i<28;i++){AvgEn[i]=0;}
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    // if(jentry != 10044 ) continue;
    // ==============print number of events done == == == == == == == =
    double progress = 10.0 * jentry / (1.0 * nentries);
    int k = int (progress);
    if (k > decade)
      cout << 10 * k << " %" << endl;
    decade = k;
    
    // ===============read this entry == == == == == == == == == == == 
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    nb2 = fChain2->GetEntry(jentry);   nbytes2 += nb2;

    // h_nrechits->Fill(NRechits);
    // std::cout << jentry << " Run " << run << "  Event " << event << " beamEnergy " << beamEnergy 
    // 	      << " NRechits " << NRechits << std::endl;

    // if(jentry<1000){
    //   cout<<"Rechit_channel= "<<rechit_channel->size();
    //   cout<<"\tRechit_chip= "<<rechit_chip->size();
    //   cout<<"\t\tRechit_layer= "<<rechit_layer->size()<<endl;
    // }
    // Assignment1(jentry);

  // float* AvgEn=Reprod(jentry);
  // h_ChivsNR->Fill(MyChiFit(layer,AvgEn),NRechits);
  // delete(AvgEn);

    // if(jentry==1)break;
    My_E_Pi_Diff();
    // if(jentry==1000){for(int i=0;i<NRechits;i++){printf("%e\n",rechit_z->at(i));}}
  } // loop over entries
}

void AnalyzeHGCMuons::Assignment1(int jentry)
{
  // int n1=rechit_layer->size();
  // for(int i=0;i<rechit_channel->size();i++){
  //   if(rechit_chip->at(i)==2){
  //     h_Chip3_rechits->Fill(rechit_layer->at(i),rechit_channel->at(i));
  //   }
  // }
  // for(int i=0;i<n1;i++){
  //   int* q=&i;
  //   int temp=rechit_layer->at(i);
  //   if(temp==28){
  //     h_layered_nrechits[27]->Fill(n1-i);
  //     *q=n1;
  //     break;
  //   }
  //   for(int j=0;j<n1-i;j++){
  //     if(temp==15){		// Loop for debugging Layer 15
  // 	h_Debug15->Fill(jentry);
  // 	h_Debug15geom->Fill(rechit_chip->at(i+j),rechit_channel->at(i+j));
  //     }
  //     if(rechit_layer->at(i+j)!=temp){
  // 	h_layered_nrechits[temp-1]->Fill(j);
  // 	*q=*q+j-1;
  // 	break;
  //     }
  //   }
  // }
  if(NRechits<100){
    for(int i=0;i<rechit_x->size();i++){
      h_Electron->Fill(rechit_z->at(i),rechit_y->at(i),rechit_x->at(i));
	}
  }
  if(NRechits>500){
    for(int i=0;i<rechit_x->size();i++){
      h_Pion->Fill(rechit_z->at(i),rechit_y->at(i),rechit_x->at(i));
	}
  }
  if(jentry==10){
    for(int i=0;i<rechit_layer->size();i++){
      h_PerEventShower[0]->Fill(rechit_layer->at(i));
      h_ShowerDev[0]->Fill(rechit_z->at(i),rechit_y->at(i),rechit_x->at(i));
      h_ShowerX[0]->Fill(rechit_x->at(i));
      h_ShowerY[0]->Fill(rechit_y->at(i));
    }
  }
  if(jentry==20){
    for(int i=0;i<rechit_layer->size();i++){
      h_PerEventShower[1]->Fill(rechit_layer->at(i));
      h_ShowerDev[1]->Fill(rechit_z->at(i),rechit_y->at(i),rechit_x->at(i));
      h_ShowerX[1]->Fill(rechit_x->at(i));
      h_ShowerY[1]->Fill(rechit_y->at(i));
      // cout<<"";
    }
  }
  if(jentry==30){
    for(int i=0;i<rechit_layer->size();i++){
      h_PerEventShower[2]->Fill(rechit_layer->at(i));
      h_ShowerDev[2]->Fill(rechit_z->at(i),rechit_y->at(i),rechit_x->at(i));
      h_ShowerX[2]->Fill(rechit_x->at(i));
      h_ShowerY[2]->Fill(rechit_y->at(i));
    }
  }
  // float in=0;
  // float fin=0;
  // int n1=15;			// number of layers under consideration
  // for(int i=0;i<rechit_layer->size();i++){
  //   if(rechit_layer->at(i)<=n1){in++;}
  //   if(rechit_layer->at(i)>28-n1){fin++;}
  // }
  // if(fin==0){h_ratio->Fill(-1);}
  // else{h_ratio->Fill(in/fin);}
  float in_NR=0;
  float in_En=0;

  float fin_NR=rechit_layer->size();
  float fin_En=0;
  for(int i=0;i<rechit_energy->size();i++){
    fin_En+=rechit_energy->at(i);
  }
  h_NRvsE->Fill(fin_En,NRechits);
  for(int n=1;n<=28;n++){
    float in_ratio=0;
    float fin_ratio=0;
    
    int n1=n;			// number of layers under consideration
    for(int i=0;i<rechit_layer->size();i++){
      // if(rechit_layer->at(i)<=n1){in_ratio++;}
      // if(rechit_layer->at(i)>28-n1){fin_ratio++;}
      if(rechit_layer->at(i)>n){break;}
      if(rechit_layer->at(i)==n){in_NR++;in_En+=rechit_energy->at(i);}
    }
    // if(fin_ratio==0){h_RatioCheck->Fill(n1,-1);}
    // else{h_RatioCheck->Fill(n1,in_ratio/fin_ratio);}
    h_CumuNRechits->Fill(n,in_NR/fin_NR);
    h_CumuEnergy->Fill(n,in_En/fin_En);
  }
}
float* AnalyzeHGCMuons::Reprod(int jentry)
{
  float* en =new float[28];		// This will store layer-wise energy
  float* LyrN =new float[28];		// This will store no. of rechits in  layer
  en[0]=0;
  LyrN[0]=0;
  float TE=0;				// Total energy detected in the event
  // float* MolEn1=new float[11];
  float MolEn=0;			// Energy fraction in 1 Moller radius
  float MR=1.5;				// Moller radius
  for(int ii=0,j=0;ii<NRechits;ii++){
    en[j]+=rechit_energy->at(ii);
    LyrN[j]++;
    TE+=rechit_energy->at(ii);
    // for(float MR=1.0;MR<2.1;MR+=0.1){
      // if(rechit_x->at(ii)*rechit_x->at(ii)+rechit_y->at(ii)*rechit_y->at(ii)<=MR*MR){MolEn1[(int)(10*MR-10)]+=rechit_energy->at(ii);}
      if(rechit_x->at(ii)*rechit_x->at(ii)+rechit_y->at(ii)*rechit_y->at(ii)<=MR*MR){MolEn+=rechit_energy->at(ii);}
    // }
    if(ii<NRechits-1&&rechit_layer->at(ii)<rechit_layer->at(ii+1)){
      j++;en[j]=0;LyrN[j]=0;
    }
  }
  for(int ii=0;ii<28;ii++){
    float k=en[ii]/TE;
    // if(pdgID==211&&(k<0.01||k>0.04)){}
    // else{

    // h_EnR[ii]->Fill(k);
    h_EnR[ii]->Fill(en[ii]);
    h_ratio_En->Fill(ii,k);
    h_NRvsER[ii]->Fill(k,NRechits);
    // }
    // For the Chi squared fit :
    // en[ii]=en[ii]/LyrN[ii];
  }
  // for(float j=0;j<11;j++){h_TestMR->Fill(NRechits,MolEn1[(int)j]/TE,1.0+j/10.0);}
  h_Tot_En->Fill(TE);
  h_NRvsMolE->Fill(NRechits,MolEn/TE);
  return(en);
}

// bool AnalyzeHGCMuons::My_P_Cut(int opt)
// {
//   if(opt==1){return(true);}
//   return(true);
// }

// float AnalyzeHGCMuons::MyChiFit(float* x, float* y)
// {
//   // The following line was supposed to make use of EnDecay() function
//   // But since that did not work , I had to write the function in the following fasion

//   TF1 *edc = new TF1("ENDECAY",[](double*x,double*p){return(p[0]*pow(x[0],p[1])*TMath::Exp(-p[2]*x[0]));},0,30,3);
//   // TF1 *edc = new TF1("ENDECAY",AnalyzeHGCMuons::EnDecay,0,30,3);
//   edc->SetParameters(10,2,0.1);
//   // edc->SetParName(0,"E");
//   // edc->SetParName(1,"a");
//   // edc->SetParName(2,"b");
//   TGraph* tg=new TGraph(28,x,y);
//   tg->Fit("ENDECAY","q");
//   // Double_t CS=edc->GetChisquare();
//   // Double_t ND=edc->GetNDF();

//   h_XvsNR->Fill(edc->GetParameter(1)/edc->GetParameter(2),NRechits);
//   // return((float)(CS/ND));
//   return(edc->GetChisquare()/edc->GetNDF());
// }

Double_t AnalyzeHGCMuons::EnDecay(Double_t* X,Double_t* par)
{
  Double_t E0=par[0];
  Double_t A=par[1];
  Double_t B=par[2];
  Double_t x=X[0];

  Double_t Ed=E0*pow(x,A)*TMath::Exp(-B*x);
  return(Ed);
}
Double_t TempEnDecay(Double_t* X,Double_t* par)
{
  Double_t E0=par[0];
  Double_t A=par[1];
  Double_t B=par[2];
  Double_t x=X[0];

  Double_t Ed=E0*pow(x,A)*TMath::Exp(-B*x);
  return(Ed);
}

void AnalyzeHGCMuons::My_E_Pi_Diff()
{
  ////////////////////1. TMax fit/////////////////////
  int EM=0;			// energy maximmum
  int NM=0;			// rechits maximmum
  float TE=0.0;			// Total detected energy
  float COG_E=0.0;		// Will store center of gravity
  vector<float>enloc,nrloc,zloc,rad,enrad,rechit_r;
  enloc.clear();    // stores total energy at corresponding z location
  nrloc.clear();    // strores total nrechits for corresponding z location
  zloc.clear();	    // stores all the available z locations (in cm)
  rad.clear();	    // will store radial dist of rechit from origin
  enrad.clear();    // will store energy at given raduis
  rechit_r.clear();
  float radbins=70; // No. of bins for h_radial
  float radend=7;   // End point of radius
  float* RAD=new float[radbins];
  for(float i=0;i<radbins;i++){RAD[(int)i]=(i+1)*radend/radbins;}
  vector<float>ENRAD;
  ENRAD.clear();
  if(NRechits>0){
    enloc.push_back(rechit_energy->at(0));
    nrloc.push_back(1.0);
    zloc.push_back(rechit_z->at(0));
    rad.push_back(sqrt(rechit_x->at(0)*rechit_x->at(0)+rechit_y->at(0)*rechit_y->at(0)));
    enrad.push_back(rechit_energy->at(0));
    rechit_r.push_back(sqrt(rechit_x->at(0)*rechit_x->at(0)+rechit_y->at(0)*rechit_y->at(0)));
    for(int i=1;i<NRechits;i++){
      h_TransProf->Fill(rechit_x->at(i),rechit_y->at(i));
      TE+=rechit_energy->at(i);
      if(rechit_z->at(i)!=zloc[zloc.size()-1]){
	enloc.push_back(rechit_energy->at(i));
	nrloc.push_back(1.0);
	zloc.push_back(rechit_z->at(i));
      }
      else{
	int n=zloc.size()-1;
	enloc[n]+=rechit_energy->at(i);
	nrloc[n]++;
      }
      float rtemp=sqrt(rechit_x->at(i)*rechit_x->at(i)+rechit_y->at(i)*rechit_y->at(i));
      rechit_r.push_back(rtemp);
      COG_E+=rechit_energy->at(i)*rechit_z->at(i);
    }
    for(int i=1;i<zloc.size();i++)
      {
	if(enloc[EM]<enloc[i]){EM=i;}
	if(nrloc[NM]<nrloc[i]){NM=i;}
      }
    h_TMaxE->Fill(EM);
    h_TMaxN->Fill(NM);
    h_TEvsE->Fill(log(trueBeamEnergy),zloc[EM]);
    h_TNvsE->Fill(log(trueBeamEnergy),zloc[NM]);
    h_NECorr->Fill(NM,EM);
    // for(int i=0;i<rad.size();i++){h_radial->Fill(rad[i],enrad[i]/TE);
    // }
    h_COG_E->Fill(COG_E/TE);
    h_CG_NR->Fill(NRechits,COG_E/TE);
  }
  ////////////////////2. Chi square fit///////////////////// 

  // TProfile* tp = new TProfile("tp","dad",30,0,60);
  // int nl=zloc.size();
  // if(nl>1){
  // float* x=new float[nl-1];
  // float* y=new float[nl-1];
  // for(int i=0;i<nl-1;i++){
  //   x[i]=zloc[i];
  //   // x[i]=i;
  //   // y[i]=()((enloc[i])/(zloc[i+1]-zloc[i]));
  //   y[i]=enloc[i]/TE;
  //   tp->Fill(x[i],y[i],1);
  // }
  // AnalyzeHGCMuons::edc->SetParameters(10,2,0.1);
  // TF1* edc = new TF1("ENDECAY",[](double*x,double*p){return(p[0]*pow(x[0],p[1])*TMath::Exp(-p[2]*x[0]));},0,50,3);
  // TF1* edc = new TF1("EnDecay",EnDecay,0,60,3);
  // Double_t TempEnDecay(Double_t*,Double_t*);
  // TF1* edc = new TF1("ENDECAY",TempEnDecay,0,60,3);
  // edc->SetParameters(1,2,0.1);

  // TCanvas* tc_temp=new TCanvas("asdfg","asdfg",800,600);
  // g_shower=new TGraph(28,x,y);
  // g_shower->Fit("ENDECAY","r",0,60);
  // g_shower->Draw("ap");
  // TFile* ff = new TFile("out.root","recreate");
  // char* r=new char[50];
  // sprintf(r,"OneEvent%devents.root",NRechits);
  // g_shower->Write();
  // tp->Write();
  // ff->Close();
  // tc_temp->SaveAs(r);
  // float E0=(float)edc->GetParameter(0);
  // float A=(float)edc->GetParameter(1);
  // float B=(float)edc->GetParameter(2);
  // h_ChiFit->Fill(edc->GetChisquare()/edc->GetNDF());
  // delete(tc_temp);
  // delete tp;

  // TGraph* tg=new TGraph(28,x,y);
  // tg->Fit("ENDECAY","QN0R");
  // // tg->Draw();
  // float E0=(float)edc->GetParameter(0);
  // float A=(float)edc->GetParameter(1);
  // float B=(float)edc->GetParameter(2);
  // h_ChiFit->Fill(edc->GetChisquare()/edc->GetNDF());
  // delete(edc);
  // delete(tg);
  // }
  // h_XvsNR->Fill(edc->GetParameter(1)/edc->GetParameter(2),NRechits);
  // return((float)(CS/ND));
  // return(edc->GetChisquare()/edc->GetNDF());
 
}
