#include <iostream>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <vector>
#include <sstream>
#include <fstream>

#include <TFile.h>
#include <TF1.h>
#include <TF2.h>
#include <TTree.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TSpectrum.h>
#include <TStyle.h>
#include <TChain.h>
#include <TLegend.h>
#include <TCut.h>
#include <TROOT.h>

#include "TCanvas.h"
#include "TMath.h"
#include "TRandom.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"

//Maths function for factorial
//this is used in the Vino Fit Function
int f(int n) { return TMath::Factorial(n); }

//Defining B which has different values depending
//on the indices
//this is also included in the Vino function
int B(int i, int k) {
    if(i == 0 && k == 0) return 1;
    else if(i == 0 && k > 0) return 0;
    else return f(k)*f(k - 1)/f(i)/f(i - 1)/f(k - i);
}
//This class contains the Vino Function that will fit the data
class LL {
private:
    vector<double> n_; //number of events in peak
    double N_; //Total number of events in spectrum
public:
    LL (vector<double> n, double N) : n_(n), N_(N) {}

    double fk (double *x, double *p, int k) {
        double lambda = x[0];
        double P = p[0];
        double fk_val = 0;
        double nk = n_[k]/N_;
        double Sum = 0;
        for(int i = 0; i <= k; i++) Sum += B(i, k)*pow(lambda*(1 - P), i)*pow(P, k - i)/f(k);
        fk_val = exp(-lambda)*Sum;
        return fk_val;
    }

    double Evaluate(double *x, double *p) {
        double lambda = x[0];
        double P = x[1];
        
        double ll = 0;
        double fk_val;
        for (size_t k = 0; k < n_.size(); k++) {
            //kth normalised value/ number of normalised value
            double nk = n_[k]/N_;
            fk_val = fk(&lambda, &P, k);            
            
            if (fk_val > 0) ll += 2*(fk_val - nk + nk*log(nk/fk_val)); //ask ashley about this
            /////
        }
        return ll;
    }
};
//Creating the Histogram 'Finger Plot'
void ManVino(int nbin=300, double nmin= -1, int nmax=20, int npeaks=15){
    
    //this sets up for looping over multiple data sets and needs to be changed depending on the data
    int numrun = 5;     //number of runs
    int run[numrun];
    int start = 65;     //starting voltage
    int dif = 2;        //difference in voltage between each run
    double snr[numrun];
    double K[numrun];
    double Prob[numrun];

    for (int jk = 0; jk <5; jk++)
    {
        run[jk] = start + (jk*dif);
    }
    double gainplot[numrun];
    double gainerrplot[numrun];
    ////////////////////////////////////////////////////////////////

    //loading the .ROOT  file
    TCanvas *c1[numrun*2];
    for (int lk = 0; lk < numrun; lk++)
    {
    printf("Analysing run voltage %dV\n,",run[lk]);
    TFile *f = new TFile(Form("data_11Nov_%dV.root",run[lk]), "read");
    //TFile *f = new TFile(Form("data_10Oct_%dV_laserOn.root",run[lk]), "read");
    //Selecting the tree from the .Root file
    TTree *data = (TTree*)f->Get("dstree");
        
        //creating a histgram hist to load the branch in
        c1[lk] = new TCanvas(Form("c%d", lk),"Charge Integral Finger Plot", 200,10,600,400);
        TH1* hist = new TH1F(Form("hist"),Form("ch_roi"), nbin, nmin, nmax);
        TF1 *sum;
        data->Draw(Form("ch_roi>>hist")); //loading the branch and adding it to hist
        data->Draw("ch_roi>>hist");
        hist->SetTitle(Form("Charge Integral %dV",run[lk]));
        hist->GetXaxis()->SetRangeUser(nmin, nmax);
        hist->GetXaxis()->SetTitle("Charge (ADC*Sample)");
        hist->GetXaxis()->SetTitleOffset(1.16);
        hist->GetYaxis()->SetTitle("Counts");
        hist->GetYaxis()->SetTitleOffset(1.16);
        hist->SetLineWidth(1);
        gStyle->SetOptFit(1111);
        hist->Sumw2();
    
        //Use TSpectrum to find the peak candidates
        int npeak;
        TSpectrum *s = new TSpectrum(2*npeaks);
        int nfound = s->Search(hist,4,"",0.005); //(name, sigma, option t, threshold)
        printf("Found %d candidate peaks to fit\n",nfound);
        
        //if there is less peaks found then what is expected use the expected
        if (nfound<=npeaks) 
        {
        npeak=nfound;
        } else 
        {
        npeak=npeaks;
        }
        printf("Found %d useful peaks to fit\n",npeak);

        //array with X-positions of the centroids found by TSpectrum
        double *xus = s->GetPositionX(); 
        //printing the peak x-positions
        /*for (int i = 0; i < npeak; i++)
        {
            printf("The %i PE peak x-position is %f\n", i, xus[i]);
        };*/
        
        //Using x-positions of the centroids to find the average distance
        vector<double> x(xus, xus + npeak);
        sort(x.begin(), x.end());
        float dmu = 0;
        //average distance between the peaks
        for (int j = 0; j < npeak - 1; j++)
        {
            dmu += x[j + 1] - x[j];
        }
        dmu = dmu/npeak-1;
        //printf("The average distance across %d peaks is %f\n",npeak,dmu);
        //printf("Now fitting\n");
        

        //Fitting
        //Fitting one gaussian to the histogram
        TF1 *g[npeak];
      for (int p=0;p<npeak;p++) 
      {
            g[p] = new TF1("gaus","gaus",x[p] -dmu, x[p] + dmu);
	
            g[p]->SetLineWidth(2);
            g[p]->SetLineColor(kRed);
            hist->Fit(g[p],"R+Q+0");
      }
      //Fitting multiple gaussians to the histogram
      string sgaus = "gaus(0) ";
      for (int ss = 1; ss < npeak; ss++)
      {
            sgaus += Form("+ gaus(%d) ", 3*ss);
      } 
      
      sum = new TF1("mysum",sgaus.c_str(),x[0] - dmu, x[npeak - 1] + dmu);
      sum->SetNpx(1000);
      //This sections sets the paramters of the gaussian
      //fixes everything but the x paramter to fit the x for k = 0,1,2,3... false, true, false, false, true, false
      //fits only x
      for (int k=0;k<3*npeak;k++)
      {
	        sum->SetParameter(k,g[(k-k%3)/3]->GetParameter(k%3));
	        if(!(k-1)%3) sum->FixParameter(k,g[(k-k%3)/3]->GetParameter(k%3));
            //look in the range +/- dmu/3
	        if(!(k-1)%3) sum->SetParLimits(k, sum->GetParameter(k) - dmu/3,sum->GetParameter(k) + dmu/3);    
      }
      
      hist->Fit(sum,"R+Q+0");
        // this refines the fit, fixes the other paramters and refines the x
       for (int k=0;k<3*npeak;k++)
      {
	        sum->SetParameter(k,g[(k-k%3)/3]->GetParameter(k%3));
	        if(!(k-1)%3) sum->ReleaseParameter(k);
            //look in the range =/- dmu/3
            if(!(k-1)%3) sum->SetParLimits(k, sum->GetParameter(k) - dmu/3,sum->GetParameter(k) + dmu/3);
      }
      
       sum->SetLineWidth(2);
       sum->SetLineColor(kBlack);
       hist->Fit(sum,"R+Q");
       vector<double>xpfit;
       vector<double>xperr;
       //extracting fit paramters and their associated errors
       for (int w=0;w<npeak;w++)
       {
           xpfit.push_back(sum->GetParameter((3*w)+1));
           xperr.push_back(sum->GetParError((3*w)+1));
           printf("The %i PE peak x-position is %f +/- %f\n", w, xpfit[w], xperr[w]);
       }
       float dmuf = 0;
       float dmuferr = 0;
        //average fitted distance between the peaks
        for (int t = 0; t < npeak -1 ; t++)
        {
            dmuf += xpfit[t + 1] - xpfit[t];
            //printf("%f\n",dmuf);
        }
        for (int v =0; v < npeak; v++)
        {
            dmuferr += pow(xperr[v],2);
        }
        dmuf = dmuf/(npeak-1);
        double xerror = sqrt(dmuferr)/(npeak-1);
        float gain = ((dmuf*10e-10)/10e3)/1.6e-19;
        float gainerr = ((xerror*10e-10)/10e3)/1.6e-19;
        gainplot[lk] = gain;
        gainerrplot[lk] = gainerr;
        printf("The average fitted distance across %d peaks is %f +/- %f\n",npeak,dmuf,xerror);
        printf("The Gain is %f +/-%f\n",gain, gainerr);
        printf("Now fitting\n");
      
      //Now log-likelihood fit to Vinogradov function
      vector<double>N; //This will be filled with the number of events in each peak
      double N_tot = hist->GetEntries(); //read the number of entries in hist
       
      // 
      for (int l=0; l<npeak; l++)
      {
	
            N.push_back(sum->GetParameter(3*l)*sqrt(2*TMath::Pi())*sum->GetParameter(2+3*l)/hist->GetBinWidth(1));
      }
	 
      LL* oll = new LL(N, N_tot);
      
      TF2* f2 = new TF2("f2", oll, &LL::Evaluate, 0.01, 10, 0.01, 0.99, 0, "LL", "Evaluate"); //Meanmin, Meanmax, Pmin, Pmax
      double lambda, P;
      //minimise log likelyhood
      f2->GetMinimumXY(lambda, P);
      
      double k_dup_ll = P/(1-P);
      P = k_dup_ll/(1+k_dup_ll);
      double ff_ll = 1.+(2.*k_dup_ll);
      cout << " " << "p: " << P << " " << " k_dup: " << k_dup_ll << " Var: " << ff_ll << endl;
      K[lk] = k_dup_ll;
      Prob[lk] = P;

      //SNR Extenstion Analysis
      
      double Varroot[numrun];
      Varroot[lk] = sqrt(ff_ll);
      double mean[numrun];
      mean[lk] = lambda/(1-P);
      //snr[lk] = mean[lk]/Varroot[lk];
      snr[lk] = sqrt(lambda)/sqrt(1+P);
      printf("The SNR is %f for %dV\n",snr[lk],run[lk]);


      //fitting the Vinofunction to the normalised events
      
      int n = (npeak);
      double xv[n], y[n];
      double VinoDist[n];
      for ( int h = 0; h < n; h++)
      {
          xv[h] = h;
          y[h] = N[h]/N_tot;
          VinoDist[h] =oll->fk(&lambda, &P, h);
          //cout << oll->fk(&lambda, &P, h) << endl;
      }

      c1[lk+numrun] = new TCanvas(Form("c%d", lk+numrun),"Normalised events", 200,10,600,400);
      TGraph *gr1 = new TGraph (n, xv, y);
      gr1->SetMarkerColor(6);
      TGraph *gr2 = new TGraph (n, xv, VinoDist);

        // Create a TMultiGraph and draw it:
        TMultiGraph *mg = new TMultiGraph();
        mg->Add(gr1);
        mg->Add(gr2);
        mg->Draw("A*");
        
        
        mg->SetTitle(Form("Number of Peak Events Normalised Compared to Model at %dV", run[lk]));
        mg->GetXaxis()->SetRangeUser(0, (npeak -1));
        mg->GetYaxis()->SetRangeUser(0, 1);
        mg->GetXaxis()->SetTitle("Peak Number");
        mg->GetYaxis()->SetTitle("Relative %");
        mg->GetXaxis()->SetTitleOffset(1.16);
        mg->GetYaxis()->SetTitleOffset(1.16);

        delete oll, f2;
    }

   TCanvas *c = new TCanvas("c","SNR", 200,10,600,400);
    Int_t p = numrun;
    Double_t x[p], y[p];
    Double_t ex[p];
    Double_t yKdup[p];
    Double_t yP[p];

    for (Int_t i=0; i<p; i++) 
    {
      x[i] = 1.0*run[i];
      y[i] = 1.0*snr[i];
      yKdup[i] = 1.0*K[i];
      yP[i] = 1.0*Prob[i];
      ex[i] = 0;
   }
    TGraph *gr = new TGraph (p, x, y);
    gr->SetTitle("SNR as a Function of Bias, Calculated from the Vinogradov Function");
    gr->GetXaxis()->SetTitle("Voltage (V)");
    gr->GetYaxis()->SetTitle("SNR");
    gr->Draw("A*");

    TCanvas *ce = new TCanvas("ce","Gain Plot", 200,10,600,400);

    TGraphErrors *gre = new TGraphErrors (p, x, gainplot, ex, gainerrplot);
    gre->SetTitle("Gain as a function of voltage");
    gre->GetXaxis()->SetTitle("Voltage (V)");
    gre->GetYaxis()->SetTitle("Gain");

   TF1 *f = new TF1("f","[0]*x + [1]", 0, run[numrun]);
   f->SetNpx(1000);
   gre->Fit(f,"Q");
   gre->Draw("A*");


   double bdv;
   double bdverr;
   float grad = f->GetParameter(0);
   float cept = f->GetParameter(1);
   float graderr = f->GetParError(0);
   float cepterr = f->GetParError(1);

   bdv = cept/(-grad);
   bdverr = bdv*(sqrt(pow(cepterr/cept,2)+pow(graderr/grad,2)));
   printf("The breakdown voltage for this data set is %f +/-%f\n", bdv,bdverr);

    
   TCanvas *cKP = new TCanvas("KP","K_dup and P as a function of bias", 200,10,600,400);
   TGraph *grK = new TGraph (p, x, yKdup);
   grK->SetMarkerColor(6);
   TGraph *grP = new TGraph (p, x, yP);

   // Create a TMultiGraph and draw it:
        TMultiGraph *mg5 = new TMultiGraph();
        mg5->Add(grK);
        mg5->Add(grP);
        mg5->SetTitle("K_dup and P as a function of Bias");
        mg5->GetXaxis()->SetTitle("Voltage (V)");
        mg5->GetYaxis()->SetTitle("Probability");
        mg5->Draw("A*");



}





