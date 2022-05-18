#include <iostream>
#include <fstream>
#include <string>

#include "TMath.h"
#include "TMinuit.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TF1.h"

std::vector<double> mResult = {0., 0., 0., 0.};
std::vector<double> mError = {0., 0., 0., 0.};

TTree *tt;
//Float_t mResidual, mTrackX, mTrackY, mR;
Float_t mResidual, mLocalPoint[3], mGlobalPoint[3];
Long64_t mEvents;

double MuonResidualsFitter_logPureGaussian(double residual, double center, double sigma) {
  sigma = fabs(sigma);
  static const double cgaus = 0.5 * log( 2.*M_PI );
  return (-pow(residual - center, 2) *0.5 / sigma / sigma) - cgaus - log(sigma);
} 

double getResidual(double delta_x, double delta_y, double delta_phiz, double track_x, double track_y, double R) {
  return delta_x - (track_x/R - 3.*pow(track_x/R, 3)) * delta_y - track_y * delta_phiz; 
}

void MuonResiduals3DOFFitter_FCN(int &npar, double *gin, double &fval, double *par, int iflag) {
  const double dx = par[0];
  const double dy = par[1];
  const double dphiz = par[2];
  const double sig = par[3];
  
  fval = 0.;
  for (Long64_t i=0;i<mEvents;i++) {
    tt->GetEntry(i);
    double residual = mResidual; double trackX = mLocalPoint[0]; double trackY = mLocalPoint[1]; double R = pow(pow(mGlobalPoint[0],2) + pow(mGlobalPoint[1],2),0.5);
    double residpeak = getResidual(dx, dy, dphiz, trackX, trackY, R);
    fval += -1.*MuonResidualsFitter_logPureGaussian(residual, residpeak, sig); 
  }
}

void doFit(bool doDx, bool doDy, bool doDphiz) {
  TMinuit mfit(4);
  mfit.SetFCN(MuonResiduals3DOFFitter_FCN);
  double par[4] = {0., 0., 0., 0.5};
  mfit.DefineParameter(0, "dx", par[0], 0.1, 0, 0); 
  mfit.DefineParameter(1, "dy", par[1], 0.1, 0, 0); 
  mfit.DefineParameter(2, "dphiz", par[2], 0.001, 0, 0); 
  mfit.DefineParameter(3, "sig", par[3], 0.01, 0, 0); 
  mfit.FixParameter(3);
  if (!doDx) mfit.FixParameter(0);
  if (!doDy) mfit.FixParameter(1);
  if (!doDphiz) mfit.FixParameter(2);

  double arglist[10];
  int ierflg;
  int smierflg;
  
  for (int i = 0;  i < 10;  i++) arglist[i] = 0.;
  arglist[0] = 0.5;
  ierflg = 0;
  smierflg = 0;
  mfit.mnexcm("SET ERR", arglist, 1, ierflg);
  for (int i = 0;  i < 10;  i++) arglist[i] = 0.;
  arglist[0] = 2;
  ierflg = 0;
  mfit.mnexcm("SET STR", arglist, 1, ierflg);

  bool try_again = false;
  for (int i = 0;  i < 10;  i++) arglist[i] = 0.;
  arglist[0] = 50000;
  ierflg = 0;
  mfit.mnexcm("MIGRAD", arglist, 1, ierflg);
  if (ierflg != 0) try_again = true;
  if (try_again){
    std::cout << "try again" << std::endl;
    for (int i = 0;  i < 10;  i++) arglist[i] = 0.;
    arglist[0] = 50000;
    mfit.mnexcm("MIGRAD", arglist, 1, smierflg);
  }

  Double_t fmin, fedm, errdef;
  Int_t npari, nparx, istat;
  mfit.mnstat(fmin, fedm, errdef, npari, nparx, istat);
  if (istat != 3) {
    for (int i = 0;  i < 10;  i++) arglist[i] = 0.;
    ierflg = 0;
    mfit.mnexcm("HESSE", arglist, 0, ierflg);
  }
  for (int i = 0;  i < 3;  i++){
    double v,e;
    mfit.GetParameter(i,v,e);
    mResult[i] = v;
    mError[i] = e;
  }
}

int main() {
  TFile *tf = new TFile("input.root");							//Name of input root file
  TTree *tmpTr = (TTree*)tf->Get("analyzer/ME11Seg_Prop");				//TTree directory in the ntuple MAKE SURE TO USE CORRECT PROPAGATION TREE
  //TTree *tmpTr = (TTree*)tf->Get("ME11ana/Inner_Prop");					//Tree name for ME11 analyzer ****

  TFile* tmpTF = new TFile("tmp1.root","recreate");					//tmp files created to stop memory errors
  std::cout << "Copying Tree" << std::endl;
  TTree *cutEn = tmpTr->CopyTree(Form("n_ME11_segment == 1 && muon_pt > 5 && abs(RdPhi) < 100 && has_fidcut"));                 //Basic cut on full tree
  std::cout << "Copied" << std::endl;
  std::cout << "Closing input" << std::endl;
  tf->Close();
  int nCuts = 2;							//Number of cuts to fit on, each cut reduces the total sample remaining by 2 (n/2, 3n/4, 7n/8, 15n/16, ...)
  for (int nCut = 1; nCut < nCuts; nCut++) {				//0 is the initial misalignment, start fitter counter at 1
    std::cout << "Starting cut number " << nCut << std::endl;
    std::cout << "Current number of entries = " << cutEn->GetEntries() << std::endl;
    if(nCut == 1){
      std::cout << "nCut 1, using full file" << std::endl;
    }
    else{
      std::cout << "Slimming tmp file" << std::endl;
      int total_entries = cutEn->GetEntries();
      int nCut_entries = int(total_entries/2.0);
      std::cout << "nCut = " << nCut << std::endl;
      std::cout << "total_entries = " << total_entries << std::endl;
      std::cout << "nCut_entries = " << nCut_entries << std::endl;
      std::cout << "Taking first half" << std::endl;
      cutEn = cutEn->CloneTree(nCut_entries);
    }
    std::cout << "New number of entries = " << cutEn->GetEntries() << std::endl;
    std::ofstream myfile;
    myfile.open (Form("out_CutNumber%d.csv",nCut));				//Name of output .csv file with suggested alignments
    double dx, dy, dz, dphix, dphiy, dphiz;
    int detNum;
    dz = 0.0; dphix = 0.0; dphiy = 0.0;
    bool doDx = true;								//Option to turn on or off 3 dof alignments
    bool doDy = true;
    bool doDphiz = true;
    std::cout << "Starting Chamber loop" << std::endl;
    for (int j = -1; j < 2; j = j + 2){
      for (int i = 0; i<36;i++){
        detNum = j*(i+101);
        std::cout << "at chamber " << detNum << std::endl;

        TFile* tmpTF = new TFile("tmp2.root","recreate");
        std::cout <<"About to copy tree" << std::endl;
        TTree* tt_tmp = cutEn->CopyTree(Form("rechit_detId==%d",detNum));		//Only fits 1 chamber at a time (det_id)
        std::cout << "Entries are on chamber are " << tt_tmp->GetEntries() << std::endl;

        TH1F *h1 = new TH1F("h1", "h1 title", 100, -20, 20);
        tt_tmp->Project("h1", "RdPhi", "");

        TF1 f1 = TF1("f1", "gaus", -2, 2);
        f1.SetParLimits(1, -2, 2);
        f1.SetParLimits(2, 0, 2);
        h1->Fit("f1", "R");
        float fitMean = f1.GetParameter(1);
        float fitStd = f1.GetParameter(2);

        std::cout << "Starting small copy" << std::endl;
        tt = tt_tmp->CopyTree(Form("RdPhi <= (%f + (1.6*%f)) && RdPhi >= (%f - (1.6*%f))", fitMean, fitStd, fitMean, fitStd));

        if (tt->GetEntries() == 0){						//If there are no events on the chamber it is skipped
          myfile << detNum << ", " << 0 << ", " << 0 << ", " << 0 << ", " << 0 << ", " << 0 << ", " << 0 << ", " << 0 << "\n";
          delete tmpTF;
          delete tt_tmp;
          continue;
        }
        tt->SetBranchAddress("RdPhi", &mResidual);			//Variables used, set for CSC propagations
        tt->SetBranchAddress("prop_LP", &mLocalPoint);
        tt->SetBranchAddress("prop_GP", &mGlobalPoint);
        mEvents = tt->GetEntries();
        doFit(doDx, doDy, doDphiz);
        dx = mResult[0];
        dy = mResult[1];
        dphiz = mResult[2];
        myfile << detNum << ", " <<dx << ", " << dy << ", " << dz << ", " << dphix << ", " << dphiy << ", " << dphiz << ", " << mEvents << "\n";
        delete tt_tmp;
        delete tmpTF;
      }
    }
    myfile.close();
  }
  //tf->Close(); 
}
