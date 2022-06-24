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

using namespace std;

// vector for recording the calculated misalignment parameters (delta_smth) and their errors
// delta_x, delta_y, delta_z, delta_phi_x, delta_phi_y, delta_phi_z, sigma_x, sigma_y, sigma_dxdz, sigma_dydz
vector<double> mResult = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
vector<double> mError = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

TTree *tt;
Float_t mResidual_dx, mResidual_dy, mResidual_dxdz, mResidual_dydz, mLocalPoint[3], mGlobalPoint[3], mLocation[3];
Long64_t mEvents;

double MuonResidualsFitter_logPureGaussian(double residual, double center, double sigma) {
  sigma = fabs(sigma);
  static const double cgaus = 0.5 * log( 2.*M_PI );
  return (-pow(residual - center, 2) *0.5 / sigma / sigma) - cgaus - log(sigma);
}

// residuals from track coordinates
double getResidual_dx(double delta_x, double delta_y, double delta_z,
                    double delta_phi_x, double delta_phi_y, double delta_phi_z,
                    double track_x, double track_y, double track_dxdz, double track_dydz)
  {
    return delta_x
        - track_dxdz * delta_z
        - track_y * track_dxdz * delta_phi_x
        + track_x * track_dxdz * delta_phi_y
        - track_y * delta_phi_z;
  }

double getResidual_dy(double delta_x, double delta_y, double delta_z,
                    double delta_phi_x, double delta_phi_y, double delta_phi_z,
                    double track_x, double track_y, double track_dxdz, double track_dydz)
  {
    return delta_y
        - track_dydz * delta_z
        - track_y * track_dydz * delta_phi_x
        + track_x * track_dydz * delta_phi_y
        + track_x * delta_phi_z;
  }

  double getResidual_dxdz(double delta_x, double delta_y, double delta_z,
                       double delta_phi_x, double delta_phi_y, double delta_phi_z,
                       double track_x, double track_y, double track_dxdz, double track_dydz)
  {
      return -track_dxdz * track_dydz * delta_phi_x
        + (1.0 + track_dxdz * track_dxdz) * delta_phi_y
        - track_dydz * delta_phi_z;
  }

  double getResidual_dydz(double delta_x, double delta_y, double delta_z,
                       double delta_phi_x, double delta_phi_y, double delta_phi_z,
                       double track_x, double track_y, double track_dxdz, double track_dydz)
  {
    return -(1.0 + track_dydz * track_dydz) * delta_phi_x
        + track_dxdz * track_dydz * delta_phi_y
        + track_dxdz * delta_phi_z;
  }

// wtf is *gin??
void MuonResiduals6DOFFitter_FCN(int &npar, double *gin, double &fval, double *par, int iflag) {
  const double delta_x = par[0];
  const double delta_y = par[1];
  const double delta_z = par[2];
  const double delta_phi_x = par[3];
  const double delta_phi_y = par[4];
  const double delta_phi_z = par[5];
  const double sig_x = par[6];
  const double sig_y = par[7]; // need separate sigma parameter for y?
  const double sig_dxdz = par[8];
  const double sig_dydz = par[9];

  fval = 0.0;
  for (Long64_t i=0; i<mEvents; i++) {
    tt->GetEntry(i);
    // int station = mLocation[1]; this parameter is not passed to this function in contrast to e.g. mLocalPoint. Why?
    double residual_dx;
    double residual_dy;
    double residual_dxdz;
    double residual_dydz;
    double track_x = mLocalPoint[0];
    double track_y = mLocalPoint[1];
    double track_dxdz = mResidual_dxdz;
    double track_dydz = mResidual_dydz;
    double residpeak_x;
    double residpeak_y;
    double slopepeak_x;
    double slopepeak_y;
    residual_dx = mResidual_dx;
    residual_dy = mResidual_dy;
    residual_dxdz = mResidual_dxdz;
    residual_dydz = mResidual_dydz;
    residpeak_x = getResidual_dx(delta_x, delta_y, delta_z, delta_phi_x, delta_phi_y, delta_phi_z, track_x, track_y, track_dxdz, track_dydz);
    residpeak_y = getResidual_dy(delta_x, delta_y, delta_z, delta_phi_x, delta_phi_y, delta_phi_z, track_x, track_y, track_dxdz, track_dydz);
    slopepeak_x = getResidual_dxdz(delta_x, delta_y, delta_z, delta_phi_x, delta_phi_y, delta_phi_z, track_x, track_y, track_dxdz, track_dydz);
    slopepeak_y = getResidual_dydz(delta_x, delta_y, delta_z, delta_phi_x, delta_phi_y, delta_phi_z, track_x, track_y, track_dxdz, track_dydz);
    fval += -1.0*MuonResidualsFitter_logPureGaussian(residual_dx, residpeak_x, sig_x);
    fval += -1.0*MuonResidualsFitter_logPureGaussian(residual_dy, residpeak_y, sig_y);
    fval += -1.0*MuonResidualsFitter_logPureGaussian(residual_dxdz, slopepeak_x, sig_dxdz);
    fval += -1.0*MuonResidualsFitter_logPureGaussian(residual_dydz, slopepeak_y, sig_dydz);
  }
}

void MuonResiduals6DOFFitter_FCN_station_4(int &npar, double *gin, double &fval, double *par, int iflag) {
  const double delta_x = par[0];
  const double delta_y = par[1];
  const double delta_z = par[2];
  const double delta_phi_x = par[3];
  const double delta_phi_y = par[4];
  const double delta_phi_z = par[5];
  const double sig_x = par[6];
  const double sig_y = par[7];
  const double sig_dxdz = par[8];
  const double sig_dydz = par[9];

  fval = 0.0;
  for (Long64_t i=0; i<mEvents; i++) {
    tt->GetEntry(i);
    double residual_dx;
    double residual_dxdz;
    double track_x = mLocalPoint[0];
    double track_y = mLocalPoint[1];
    double track_dxdz = mResidual_dxdz;
    double track_dydz = 0.0; // need to initialze y slope by zero in station 4 to calculate slopepeak_x correctly, since all getResidual_ functions have the same signature
    double residpeak_x;
    double slopepeak_x;
    residual_dx = mResidual_dx;
    residual_dxdz = mResidual_dxdz;
    residpeak_x = getResidual_dx(delta_x, delta_y, delta_z, delta_phi_x, delta_phi_y, delta_phi_z, track_x, track_y, track_dxdz, track_dydz);
    slopepeak_x = getResidual_dxdz(delta_x, delta_y, delta_z, delta_phi_x, delta_phi_y, delta_phi_z, track_x, track_y, track_dxdz, track_dydz);
    fval += -1.0*MuonResidualsFitter_logPureGaussian(residual_dx, residpeak_x, sig_x);
    fval += -1.0*MuonResidualsFitter_logPureGaussian(residual_dxdz, slopepeak_x, sig_dxdz);
  }
}

void doFit(bool do_delta_x, bool do_delta_y, bool do_delta_z,
    bool do_delta_phi_x, bool do_delta_phi_y, bool do_delta_phi_z, int station)
  {
      TMinuit mfit(10);
      if (station == 4) {
          mfit.SetFCN(MuonResiduals6DOFFitter_FCN_station_4);
      }
      else {
          mfit.SetFCN(MuonResiduals6DOFFitter_FCN);
      }
      // par[10] = {delta_x, delta_y, delta_z, delta_phi_x, delta_phi_y, delta_phi_z, sig_x, sig_y, sig_dxdz, sig_dydz}
      double par[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5};
      mfit.DefineParameter(0, "delta_x", par[0], 0.1, 0, 0);
      mfit.DefineParameter(1, "delta_y", par[1], 0.1, 0, 0);
      mfit.DefineParameter(2, "delta_z", par[2], 0.1, 0, 0);
      mfit.DefineParameter(3, "delta_phi_x", par[3], 0.001, 0, 0);
      mfit.DefineParameter(4, "delta_phi_y", par[4], 0.001, 0, 0);
      mfit.DefineParameter(5, "delta_phi_z", par[5], 0.001, 0, 0);
      mfit.DefineParameter(6, "sig_x", par[6], 0.01, 0, 0);
      mfit.DefineParameter(7, "sig_y", par[7], 0.01, 0, 0);
      mfit.DefineParameter(8, "sig_dxdz", par[8], 0.01, 0, 0);
      mfit.DefineParameter(9, "sig_dydz", par[9], 0.01, 0, 0);

      // fix sigmas
      mfit.FixParameter(6);
      mfit.FixParameter(7);
      mfit.FixParameter(8);
      mfit.FixParameter(9);

      if (do_delta_x == false) mfit.FixParameter(0);
      if (do_delta_y == false) mfit.FixParameter(1);
      if (do_delta_z == false) mfit.FixParameter(2);
      if (do_delta_phi_x == false) mfit.FixParameter(3);
      if (do_delta_phi_y == false) mfit.FixParameter(4);
      if (do_delta_phi_z == false) mfit.FixParameter(5);

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
      // i<6 for each of 6 misalignments we're fitting
      for (int i = 0;  i < 6;  i++){
        double v,e;
        mfit.GetParameter(i,v,e);
        mResult[i] = v;
        mError[i] = e;
      }
  }

int main(){
    // input file
    TFile *tf = new TFile("/afs/cern.ch/user/a/abolshov/private/CMSSW_12_3_0/src/GEMCSCBendingAnalyzer/GEM_Alignment/test/DT_collisionMC_062222.root");
    // directory in the .root files; here the directory for DT_tbma_collisionMC_allRes.root
    TTree *tmpTr = (TTree*)tf->Get("DT_tbma/Inner_Prop_ChamberLevel");
    cout << "Entries in tmpTr: " << tmpTr->GetEntries() << endl;

    //tmp files created to stop memory errors
    TFile* tmpTF = new TFile("tmp1.root", "recreate");
    cout << "Copying Tree" << endl;
    // Basic cut on full tree. Which cut should I impose for DTs?
    TTree *cutEn = tmpTr->CopyTree(Form("pt > 5 && abs(res_dx) < 100 && (abs(res_dy) < 100 && abs(res_dxdz) < 100 && abs(res_dydz) < 100 || location[1] == 4)"));
    cout << "Copied" << endl;
    cout << "Closing input" << endl;
    tf->Close();
    cout << "entries " << cutEn->GetEntries() << endl;


    // i'll be using the full file, so i'll skip bulshit with nCuts and stuff
    // cloning the tree
    cutEn = cutEn->CloneTree();
    cout << "Entries in cutEn: " << cutEn->GetEntries() << endl;

    ofstream myfile;
    myfile.open(Form("out.csv"));
    double delta_x, delta_y, delta_z, delta_phi_x, delta_phi_y, delta_phi_z;

    // choose which dof we align
    bool do_delta_x = true;
    bool do_delta_y = true;
    bool do_delta_z = true;
    bool do_delta_phi_x = true;
    bool do_delta_phi_y = true;
    bool do_delta_phi_z = true;

    // looping over all sectors
    cout << "Starting Sector loop" << endl;
    for (int wheel = -2; wheel < 3; wheel++)
    {

        cout << "Wheel " << wheel << ":" << endl;
        // do_delta_y = true;

        for (int station = 1; station < 5; station++)
        {

            cout << "Station " << station << ":" << endl;
            // y is not aligned in station 4, turn of fit of delta_y
            if (station == 4) {
                do_delta_y = false;
                delta_y = 0.0;
            }

            for (int sector = 1; sector < 13; sector++)
            {

                //copying tree corresponding to the current sector
                TFile* tmpTF1 = new TFile("tmp2.root", "recreate");
                TString cut = Form("location[0] == %i && location[1] == %i && location[2] == %i", wheel, station, sector);
                TTree* tt_tmp = cutEn->CopyTree(Form(cut));
                cout << "Entries in sector " << sector << "  are " << tt_tmp->GetEntries() << endl;
                TString sectorInfo = Form("%i/%i/%i", wheel, station, sector);

                // define TH1F here and extract some parameters from the histogram;
                // extrat sigma and mean value needed for minimization

                TH1F *h_x = new TH1F("h_x", "h_x title", 100, -20, 20);
                tt_tmp->Project("h_x", "res_dx", "");

                TF1 f_x = TF1("f_x", "gaus", -2, 2);
                f_x.SetParLimits(1, -2, 2);
                f_x.SetParLimits(2, 0, 2);
                h_x->Fit("f_x", "R");
                float fitMean_x = f_x.GetParameter(1);
                float fitStd_x = f_x.GetParameter(2);

                TH1F *h_slope_x = new TH1F("h_slope_x", "h_slope_x title", 100, -20, 20);
                tt_tmp->Project("h_slope_x", "res_dxdz", "");

                TF1 f_slope_x = TF1("f_slope_x", "gaus", -2, 2);
                f_slope_x.SetParLimits(1, -2, 2);
                f_slope_x.SetParLimits(2, 0, 2);
                h_slope_x->Fit("f_slope_x", "R");
                float fitMean_slope_x = f_slope_x.GetParameter(1);
                float fitStd_slope_x = f_slope_x.GetParameter(2);

                // conditiion for making smaller copy
                TString cond = Form("res_dx <= (%f + (1.6*%f)) && res_dx >= (%f - (1.6*%f)) && res_dxdz <= (%f + (1.6*%f)) && res_dxdz >= (%f - (1.6*%f))", fitMean_x, fitStd_x, fitMean_x, fitStd_x, fitMean_slope_x, fitStd_slope_x, fitMean_slope_x, fitStd_slope_x);

                if (station != 4) {
                    TH1F *h_y = new TH1F("h_y", "h_y title", 100, -20, 20);
                    tt_tmp->Project("h_y", "res_dy", "");

                    TF1 f_y = TF1("f_y", "gaus", -2, 2);
                    f_y.SetParLimits(1, -2, 2);
                    f_y.SetParLimits(2, 0, 2);
                    h_y->Fit("f_y", "R");
                    float fitMean_y = f_y.GetParameter(1);
                    float fitStd_y = f_y.GetParameter(2);

                    TH1F *h_slope_y = new TH1F("h_slope_y", "h_slope_y title", 100, -20, 20);
                    tt_tmp->Project("h_slope_y", "res_dydz", "");

                    TF1 f_slope_y = TF1("f_slope_y", "gaus", -2, 2);
                    f_slope_y.SetParLimits(1, -2, 2);
                    f_slope_y.SetParLimits(2, 0, 2);
                    h_slope_y->Fit("f_slope_y", "R");
                    float fitMean_slope_y = f_slope_y.GetParameter(1);
                    float fitStd_slope_y = f_slope_y.GetParameter(2);

                    cond.Append(Form(" && res_dy <= (%f + (1.6*%f)) && res_dy >= (%f - (1.6*%f)) && res_dydz <= (%f + (1.6*%f)) && res_dydz >= (%f - (1.6*%f))", fitMean_y, fitStd_y, fitMean_y, fitStd_y, fitMean_slope_y, fitStd_slope_y, fitMean_slope_y, fitStd_slope_y));
                }

                cout << "Starting small copy" << endl;
                tt = tt_tmp->CopyTree(Form(cond));
                cout << "Number of entries in the small copy is " << tt->GetEntries() << endl;

                if (tt->GetEntries() == 0) {						//If there are no events on the chamber it is skipped
                  myfile << sectorInfo << ", " << 0 << ", " << 0 << ", " << 0 << ", " << 0 << ", " << 0 << ", " << 0 << ", " << 0 << "\n";
                  delete tmpTF1;
                  delete tt_tmp;
                  continue;
                }

                tt->SetBranchAddress("res_dx", &mResidual_dx);
                tt->SetBranchAddress("res_dy", &mResidual_dy);
                tt->SetBranchAddress("res_dxdz", &mResidual_dxdz);
                tt->SetBranchAddress("res_dydz", &mResidual_dydz);
                tt->SetBranchAddress("sectorLevel_LP", &mLocalPoint);
                tt->SetBranchAddress("sectorLevel_GP", &mGlobalPoint);
                tt->SetBranchAddress("location", &mLocation);
                mEvents = tt->GetEntries();
                doFit(do_delta_x, do_delta_y, do_delta_z, do_delta_phi_x, do_delta_phi_y, do_delta_phi_z, station);
                delta_x = mResult[0];
                delta_y = mResult[1];
                delta_z = mResult[2];
                delta_phi_x = mResult[3];
                delta_phi_y = mResult[4];
                delta_phi_z = mResult[5];
                myfile << sectorInfo << ", " << delta_x << ", " << delta_y << ", " << delta_z << ", " << delta_phi_x << ", " << delta_phi_y << ", " << delta_phi_z << ", " << mEvents << "\n";
                delete tt_tmp;
                delete tmpTF1;

            }

            //turn delta_y fit back on for all other stations
            do_delta_y = true;

        }
    }
    myfile.close();

}
