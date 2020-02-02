/*
* LGADFits.cxx
*
*
*      Author: Gkougkousis Evangelos - Leonidas
*              egkougko@cern.ch
*               IFAE-BARCELONA
*/

#include "LGADUtils/LGADBase.h"
#include "TVirtualFFT.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooFFTConvPdf.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooPolynomial.h"
#include "RooPlot.h"
#include "TFitResult.h"
#include "TVirtualFitter.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#ifndef ROO_MSG_SERVICE
#include "RooMsgService.h"
#endif

using namespace RooFit;

// ROOT FIT RESULTS STATUS CODES
/*
* Fit status parameter for MINUIT optimizer using MINOS
* 0 - Sucssesfull fit
* <0 - Error not associated with minimization proccess (eg incorrect funciton)
* 1 - MIGRAD command is blank, ignored
* 2 - MIGRAD command line unreadable, ignored
* 3 - unknown MIGRAD command, ignored
* 4 - MIGRAD not converged
* 10 - MIGRAD END command
* 11 - MIGRAD EXIT or STOP command
* 12 - MIGRAD RETURN command
*/

// ROOFIT FIT RESULTS STATUS CODES
/*
* Fit quality parameter  MINUIT and NEWMINUIT optimizers
* (3 - Full accurate covariafornce matrix,
* 2 - Full matrix, but forced positive-definite (i.e. not accurate),
* 1 - Diagonal approximation only, not accurate,
* 0 - Error matrix not calculated at all)
* */

int LGADBase::IterativeFit(std::vector<double> *w, std::pair<double, double> &gmean, std::pair<double, double> &gsigma, TH1D* &FitHist,
                           double &minchi2, std::string methode, std::pair<int, int> points)
{
    gmean = std::make_pair(-1., -1.);
    gsigma = std::make_pair(-1., -1.);
    minchi2 = -1.;
    int last, first;
    if (points.second < 0) last = w->size();
    else last = points.second;
    if (points.first < 0) first = 0;
    else first = points.first;

    if ((methode == "Gauss" && abs(last - first) < 5) || (methode == "LandauXGauss" && abs(last - first) < 10))
       {
        if (m_verbose > 0) std::cout << __FUNCTION__ << " WARNING: Inadequate number of points for " << methode << " calculation["
                                     << abs(last - first) << "] -> will not cintimue!" << std::endl;
        return -1;
       }

    // Set and populate the 1D histogram on which the fit will be perfirmed
    double mean = Mean(w, first, last);
    double strdv = Stdev(w, first, last);
    double wmax = *std::max_element(w->begin() + first, w->begin() + last);
    double wmin = *std::min_element(w->begin() + first, w->begin() + last);

    if (mean == -99 || strdv == -99 || strdv == 0)
       {
        if (m_verbose > 0) std::cout << __FUNCTION__ << " WARNING: Unabale to determine mean and stdv ["
                                     << mean << ", " << strdv << "] -> will not cintimue!" << std::endl;
        return -2;
       }

    // Histogtam bins and limits
    m_bins = 0;
    double limit1 = 0.0;
    double limit2 = 0.0;
    double chi2 = 0.0;
    bool stop = false;
    // Variables for fit range
    double rmin = 0.0;
    double rmax = 0.0;
    int quality = -99;
    std::vector<double> mag;
    mag.reserve(35);
    std::vector<double> magErr;
    magErr.reserve(35);
    std::vector<double> Sigma;
    Sigma.reserve(35);
    std::vector<double> SigmaErr;
    SigmaErr.reserve(35);
    std::vector<double> ChiSq;
    ChiSq.reserve(35);
    std::vector<int> globindx;
    globindx.reserve(35);
    std::vector<TH1D* > Fits;
    Fits.resize(35, NULL);
    std::vector<unsigned int> itr;
    itr.reserve(35);

    TAxis *xaxis;
    Int_t bin1 = 0;
    Int_t bin2 = 0;

    if (m_verbose == 2) std::cout << __FUNCTION__ << " INFO: Calculating before fit| Average: " << mean << " , Standard Deviation: " 
                                  << strdv << ", Min Value: " << wmin << ", Max Value: "  << wmax << std::endl;

    double fact = pow(10, fabs(floor(log10(fabs(mean)))));
    // iterative re-fitting 1, affects the limits of the histogram only
    for (double gk = 7; gk > 2; gk--) // 5 cases
        {
         if (stop) break;
         if (methode == "Gauss")
            {
             limit1 = floor((mean-(gk/2)*strdv)*fact) / fact; // Symentric limit variaton
             limit2 = ceil((mean+(gk/2)*strdv)*fact) / fact;
            }
         else if (methode == "LandauXGauss")
                 {
                  limit1 = floor((wmin + ((gk-5)/10)*strdv) * fact) / fact; // Asymentric limit variaton
                  limit2 = ceil((wmax - ((gk-5)/2)*strdv) * fact) / fact;
                 }
         // Iterative re-fitting 2, affects only the width of the bin but not the range
         for (double tmb = 14; tmb > 7; tmb--) // 7 cases
             {
              int iter = ((7-gk)*7)+(14-tmb);
              if (stop) break;
              if (methode == "Gauss") m_bins = ceil(sqrt(last - (first+1))*(fabs(limit2 - limit1)/((tmb/4)*strdv)));
              else if (methode == "LandauXGauss") m_bins = ceil(sqrt(last - (first+1))*(fabs(limit2 - limit1)/(0.5*strdv*(tmb/4))));
              Fits.at(iter) = new TH1D(Form("MagFit%02u", iter), Form("MagFit%02u", iter), m_bins, limit1, limit2);
              itr.push_back(iter);
              for (int i = first; i < last; i++) Fits.at(iter)->Fill(w->at(i));
              if (Fits.at(iter)->Integral("width") == 0)
                 {
                  if (m_verbose > 0) std::cout << __FUNCTION__ << " WARNING: Integral of " << methode 
                                               << " distribution is 0 -> will not continue!" << std::endl;
                  continue;
                 }
              double mx = Fits.at(iter)->GetMean();
              double rms = Fits.at(iter)->GetRMS();
              if (rms == 0) 
                 {
                  if (m_verbose > 0) std::cout << __FUNCTION__ << " WARNING: Standard deviation of " << methode 
                                               << " distribution is 0 -> will not continue!" << std::endl;
                  continue;
                 }
              if (m_verbose == 2) std::cout << __FUNCTION__ << " INFO: Magnitude histo mean: " << mx << " , stdev: " << rms << ", integral: " 
                                            << Fits.at(iter)->Integral("width") << ", bins: " << m_bins << ", from " << limit1 << " to " << limit2 << std::endl;
              TFitResultPtr fitResult(0);
              xaxis = Fits.at(iter)->GetXaxis();
              if (methode == "Gauss" || methode == "GaussInt")
                 {
                  rmin = mx - (gk/2) * rms; // symetric rangre
                  rmax = mx + (gk/2) * rms;
                  bin1 = xaxis->FindBin(rmin);
                  bin2 = xaxis->FindBin(rmax);
                  if (((bin2-bin1) < 2) || (Fits.at(iter)->Integral(bin1 + 1, bin2 - 1) == 0)) continue;
                  if (methode == "GaussInt") fitResult = LGADBase::Gauss(rmin, rmax, strdv, mean, Fits.at(iter), "I");
                  else fitResult = LGADBase::Gauss(rmin, rmax, strdv, mean, Fits.at(iter), "N");
                 }
              else if (methode == "LandauXGauss" || methode == "LandauXGaussInt")
                      {
                       rmin = wmin + ((gk-5)/10)*rms; // asymentric range
                       rmax = wmax - ((gk-5)/2)*rms;
                       bin1 = xaxis->FindBin(rmin);
                       bin2 = xaxis->FindBin(rmax);
                       if (((bin2 - bin1) < 2) || (Fits.at(iter)->Integral(bin1 + 1, bin2 - 1) == 0)) continue;
                       if (methode == "LandauXGaussInt") fitResult = LGADBase::GauXLandau(rmin, rmax, strdv, Fits.at(iter), "I");
                       else fitResult = LGADBase::GauXLandau(rmin, rmax, strdv, Fits.at(iter), "N");
                      }
              else {
                    std::cout << __FUNCTION__ << " ERROR: Uknown fit function!" << std::endl;
                    // Cleanup Time
                    for (unsigned int k = 0; k < itr.size(); k++) 
                        {
                         delete Fits.at(itr.at(k));
                         Fits.at(itr.at(k)) = NULL;
                        }
                    return -3;
                   }
              if (fitResult->IsEmpty() || !(fitResult->IsValid()))
                 {
                  if (m_verbose == 2) std::cout << __FUNCTION__ << " WARNING: Fit result not there -> will not calculate magnitude!" << std::endl;
                  quality = -4;
                 }
              else {
                    quality = fitResult;
                    if (quality == 0 && (fitResult->Ndf()) != 0)
                       {
                        chi2 = fabs(1 - (fitResult->MinFcnValue() / fitResult->Ndf()));
                        if (methode == "Gauss")
                           {
                            mag.push_back(fitResult->Parameter(1));
                            magErr.push_back(fitResult->ParError(1));
                            Sigma.push_back(fitResult->Parameter(2));
                            SigmaErr.push_back(fitResult->ParError(2));
                            ChiSq.push_back(chi2);
                            globindx.push_back(iter);
                            if (chi2 < 0.03) stop = true;
                           }
                        else if (methode == "LandauXGauss")
                                {
                                 double max = -99;
                                 double FWHM = -99;
                                 if (langaupro(fitResult, max, FWHM) == 0)
                                    {
                                     mag.push_back(max);
                                     magErr.push_back(fitResult->ParError(1));
                                     Sigma.push_back(FWHM/2);
                                     SigmaErr.push_back(sqrt(pow(2 * fitResult->ParError(0), 2) + pow(2 * fitResult->ParError(3), 2)));
                                     ChiSq.push_back(chi2);
                                     globindx.push_back(iter);
                                     if (chi2 < 0.03) stop = true;
                                    }
                                }
                       }
                    }
              if (m_verbose == 2) std::cout << __FUNCTION__ << " INFO: Quality of performed fit: " << quality << " at iteration: " 
                                            << gk << " " << tmb << std::endl;
             }
        }

    if (quality < 0 && ChiSq.size() == 0)
       {
        if (m_verbose > 0) std::cout << __FUNCTION__ << " WARNING: Fit result not there -> will not calculate mean value!" << std::endl;
        // Cleanup Time
        for (unsigned int k = 0; k < itr.size(); k++)
            {
             delete Fits.at(itr.at(k));
             Fits.at(itr.at(k)) = NULL;
            }
        return -4;
       }
    else if (quality > 0 && ChiSq.size() == 0)
            {
             if (m_verbose > 0) std::cout << __FUNCTION__ << " WARNING: Convolution quality is bad [" << quality 
                                          << "] -> will not calculate magnitutde!" << std::endl;
             // Cleanup Time
             for (unsigned int k = 0; k < itr.size(); k++) 
                 {
                  delete Fits.at(itr.at(k));
                  Fits.at(itr.at(k)) = NULL;
                 }
             return -5;
            }
    else if (quality == 0 && ChiSq.size() == 0)
            {
             if (m_verbose > 0) std::cout << __FUNCTION__ << " WARNING: Fit calculated values deviate from Standard Deviation and Average" << std::endl;
             // Cleanup Time
             for (unsigned int k = 0; k < itr.size(); k++)
                 {
                  delete Fits.at(itr.at(k));
                  Fits.at(itr.at(k)) = NULL;
                 }
             return -6;
            }
    else {
          int indx = std::distance(ChiSq.begin(), std::min_element(ChiSq.begin(), ChiSq.end()));
          gmean = std::make_pair(mag.at(indx), magErr.at(indx));
          gsigma = std::make_pair(Sigma.at(indx), SigmaErr.at(indx));
          FitHist = (TH1D*)(Fits.at(globindx.at(indx)))->Clone("OutHist");
          // Fits.at(globindx.at(indx))->SetDirectory(0);
          // Fits.at(globindx.at(indx))->Draw();
          minchi2 = ChiSq.at(indx);
          if (m_verbose == 2) std::cout << __FUNCTION__ << " INFO: Fit result " << ChiSq.at(indx) << " --> " << mag.at(indx) << " +/- " 
                                        << magErr.at(indx) << " , " << Sigma.at(indx) << " +/- " << SigmaErr.at(indx) << std::endl;
          // Cleanup Time
          for (unsigned int k = 0; k < itr.size(); k++)
              {
               delete Fits.at(itr.at(k));
               Fits.at(itr.at(k)) = NULL;
              }
          return 0;
         }

}
// --------------------------------------------------------------------------------------------------------------
TFitResultPtr LGADBase::Gauss(double rmin, double rmax, double strdv, double mean, TH1D* &magHist, std::string integ)
{
    TF1* mygau = new TF1("GaussFit", "gaus(0)", rmin, rmax);
    double limitdown = mean - (5 * strdv);
    double limitup = mean + (5 * strdv);
    if (limitdown < rmin) limitdown = rmin;
    if (limitup > rmax) limitup = rmax;
    mygau->SetParameter(0, 0.3989422804014*(1/strdv));
    if (mean < limitdown || mean > limitup)
       {
        double a = fabs(mean - limitdown);
        double b = fabs(limitup - mean);
        if (b < a) mygau->SetParameter(1, (b*(limitup-limitdown)/((limitup-limitdown)+2*b)) + (limitdown + ((limitup - limitdown)/2)));
        else mygau->SetParameter(1, (limitdown + ((limitup - limitdown) / 2)) - (a*(limitup - limitdown) / ((limitup - limitdown) + 2 * a)));
       }
    else mygau->SetParameter(1, mean);
    mygau->SetParameter(2, strdv);
    //mygau->SetParLimits(0, 0.3989422804014*(1/(0.01*strdv)), 0.3989422804014*(1/(10*strdv)));
    mygau->SetParLimits(1, limitdown, limitup);
    mygau->SetParLimits(2, 0.3*strdv, 3*strdv);
    TFitResultPtr fitResult;
    // Fit within specified range, use ParLimits, do not plot
    if (integ == "I") fitResult = magHist->Fit(mygau, "S Q N R I");
    else fitResult = magHist->Fit(mygau, "S Q N R");
    if (m_verbose == 2) fitResult->Print("V");
    delete mygau;
    return fitResult;
}
// --------------------------------------------------------------------------------------------------------------
TFitResultPtr LGADBase::GauXLandau(double rmin, double rmax, double strdv, TH1D* &magHist, std::string integ)
{
    TF1 *mygauland = new TF1("LandXGauFun", this, &LGADBase::LandXGauFun, rmin, rmax, 4,"LGADUtils", "LandXGauFun");
    mygauland->SetParNames("Width", "MPV", "Area", "GaussSigma");
    TAxis *axis = magHist->GetXaxis();
    double wmin = axis->GetBinCenter(magHist->FindFirstBinAbove(0, 1));
    double wmax = axis->GetBinCenter(magHist->FindLastBinAbove(0, 1));
    double MPV = axis->GetBinCenter(magHist->GetMaximumBin());
    double integral = magHist->Integral("width");
    double limitdown = MPV - 3 * strdv;
    double limitup = MPV + 3 * strdv;
    if (limitdown < rmin) limitdown = rmin;
    if (limitup > rmax) limitup = rmax;
    if (MPV < limitdown || MPV > limitup)
       {
        double a = fabs(MPV - limitdown);
        double b = fabs(limitup - MPV);
        if (b < a) mygauland->SetParameter(1, (b*(limitup - limitdown) / ((limitup - limitdown) + 2 * b)) + (limitdown + ((limitup - limitdown) / 2)));
        else mygauland->SetParameter(1, (limitdown + ((limitup - limitdown) / 2)) - (a*(limitup - limitdown) / ((limitup - limitdown) + 2 * a)));
       }
    else mygauland->SetParameter(1, MPV);
    mygauland->SetParameter(0, (wmax-wmin)/10);    
    mygauland->SetParameter(2, integral);
    mygauland->SetParameter(3, strdv/4);
    //mygauland->SetParLimits(0, (wmax - wmin)/100, 100*(wmax - wmin));
    mygauland->SetParLimits(1, rmin, rmax);
    mygauland->SetParLimits(2, 0.1*integral, 10*integral);
    //mygauland->SetParLimits(3, 0.01*(strdv/2), 100*strdv/2);

    // Fit within specified range, use ParLimits, do not plot
    TFitResultPtr fitResult;
    if (integ == "I") fitResult = magHist->Fit(mygauland, "S Q N R I");
    else fitResult = magHist->Fit(mygauland, "S Q N R");
    if (m_verbose == 2) fitResult->Print("V");
    delete mygauland;
    return fitResult;
}
// --------------------------------------------------------------------------------------------------------------
Double_t LGADBase::LandXGauFun(Double_t *x, Double_t *par)
{
    // Fit parameters:
    // par[0] = Width (scale) parameter of Landau density
    // par[1] = Most Probable (MPV, location) parameter of Landau density
    // par[2] = Total area (integral -inf to inf, normalization constant)
    // par[3] = Width (sigma) of convoluted Gaussian function
    // In the Landau distribution (represented by the CERNLIB approximation),
    // the maximum is located at x=-0.22278298 with respect to the 
    // location parameter 0. This shift is corrected within this function,
    // so that the actual maximum is identical to the MPV parameter.

    // Numeric constants
    Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)

    // Control constants
    double sc = 5.0;        // convolution extends to +-sc Gaussian sigmas

    // MP shift correction
    double mpc = par[1] - 0.22278298 * par[0];

    // Range of convolution integral
    double xlow = x[0] - sc * par[3];
    double xupp = x[0] + sc * par[3];
    double step = (xupp - xlow) / m_bins;

    double xx;
    double fland;
    double sum = 0.0;
    // Convolution integral of Landau and Gaussian by sum
    for (double i = 1.0; i <= m_bins/2; i++)
        {
         xx = xlow + (i - .5) * step;
         fland = TMath::Landau(xx, mpc, par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0], xx, par[3]);
         xx = xupp - (i - .5) * step;
         fland = TMath::Landau(xx, mpc, par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0], xx, par[3]);
        }

    return (par[2] * step * sum * invsq2pi / par[3]);
}
// --------------------------------------------------------------------------------------------------------------
int LGADBase::langaupro(TFitResultPtr fitResult, double &maxx, double &FWHM)
{
    // Seaches for the location (x value) at the maximum of the
    // Landau-Gaussian convolute and its full width at half-maximum.
    // The search is probably not very efficient, but it's a first try.

    double params[4];
    for (unsigned int k = 0; k < 4; k++) params[k] = fitResult->Parameter(k);

    double p, x, half, fxr, fxl;
    double step;
    double l, lold;
    int i = 0;
    int MAXCALLS = 10000;

    // Search for maximum
    p = params[1] - 0.1 * params[0];
    step = 0.05 * params[0];
    lold = -2.0;
    l = -1.0;
    while ((l != lold) && (i < MAXCALLS)) 
          {
           i++;
           lold = l;
           x = p + step;
           l = LandXGauFun(&x, params);
           if (l < lold) step = -step / 10;
           p += step;
          }
    if (i == MAXCALLS) return (-1);
    maxx = x;
    half = l / 2;

    // Search for right x location of half
    p = maxx + params[0];
    step = params[0];
    lold = -2.0;
    l = -1e300;
    i = 0;
    while ((l != lold) && (i < MAXCALLS)) 
          {
           i++;
           lold = l;
           x = p + step;
           l = TMath::Abs(LandXGauFun(&x, params) - half);
           if (l > lold) step = -step / 10;
           p += step;
          }
    if (i == MAXCALLS) return (-2);
    fxr = x;

    // Search for left x location of half
    p = maxx - 0.5 * params[0];
    step = -params[0];
    lold = -2.0;
    l = -1e300;
    i = 0;
    while ((l != lold) && (i < MAXCALLS)) 
          {
           i++;
           lold = l;
           x = p + step;
           l = TMath::Abs(LandXGauFun(&x, params) - half);
           if (l > lold) step = -step / 10;
           p += step;
          }
    if (i == MAXCALLS) return (-3);
    fxl = x;

    // Calculate FWHM
    FWHM = fxr - fxl;

    return (0);
}
// --------------------------------------------------------------------------------------------------------------
int LGADBase::RooConvFit(std::vector<double>* vec, std::pair<double, double> &magMPV, std::pair<double, double> &magSigma, std::string conv)
{
    unsigned int npoints = vec->size();
    if (npoints < 10)
       {
        std::cout << __FUNCTION__ << " WARNING: Less than 10 points for this calculation["
                  << npoints << "] -> will not cintimue!" << std::endl;
        return -1;
       }

    // Set and populate the 1D histogram on which the fit will be perfirmed
    double mean = Mean(vec);
    double strd = Stdev(vec);

    if (mean == -99 || strd == -99 || strd <= 0)
       {
        std::cout << __FUNCTION__ << " WARNING: Unabale to determine meand and stdv ["
                  << mean << "," << strd << "] -> will not cintimue!" << std::endl;
        return -2;
       }

    double fact = pow(10, fabs(floor(log10(fabs(mean)))));
    double limit1 = floor((mean - 5 * strd) * fact) / fact;
    double limit2 = ceil((mean + 5 * strd) * fact) / fact;
    int bins = ceil(fabs(limit2 - limit1) / (strd / 3));
    // Bin increase in high statistics mode
    if (npoints / ((3 / 2) * bins) > 1) bins = ceil((npoints / ((3 / 2)*bins))*bins);
    TH1F* mag = new TH1F("mag", "mag", bins, limit1, limit2);
    for (unsigned int a = 0; a < vec->size(); a++) mag->Fill(vec->at(a));
    if (mag->Integral() <= 0)
       {
        std::cout << __FUNCTION__ << " ERROR: Integral of magnitude distribution is <=0 -> will not continue!" << std::endl;
        delete mag;
        return -3;
       }
    if (m_verbose == 2) std::cout << __FUNCTION__ << " INFO: Calculating mean and std before fit: " << mean << " , " << strd << std::endl;

    // Start RooFit
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    RooMsgService::instance().setSilentMode(true);

    // Construct observable
    RooRealVar t("t", "t", limit1, limit2);

    // Define fit range
    double rmin = mean - 5 * strd;
    double rmax = mean + 5 * strd;
    t.setRange("range", rmin, rmax);
    t.setBins(5000, "cache");

    // Construct gauss(t,mg,sg)
    RooRealVar mg("mg", "mg", 0.0, -2*strd, 2*strd);
    RooRealVar sg("sg", "sg", strd, 0.001*strd, 2*strd);
    RooGaussian gauss("gauss", "gauss", t, mg, sg);

    RooRealVar var1;
    RooRealVar var2;
    RooFFTConvPdf* Conv;
    RooLandau landau;
    RooPolynomial poly1;

    if (conv == "LanXGau")
       { 
        // Construct landau(t,ml,sl) ;
        var1 = RooRealVar("ml", "mean landau", mean, -4 * fabs(mean), 2 * fabs(mean));
        var2 = RooRealVar("sl", "sigma landau", strd, 0.001*strd, 2 * strd);
        landau = RooLandau("lx", "lx", t, var1, var2);
        // Construct landau (x) gauss
        Conv = new RooFFTConvPdf("lxg", "landau (X) gauss", t, landau, gauss, 2);
       }
    else if (conv == "LinXGau")
            {
             // Construct linear(t,slope,interc)
             var1 = RooRealVar("ilin", "intersect linear", mean, mean - 3 * strd, mean + 3 * strd);
             var2 = RooRealVar("slin", "slope linear", strd, 0.01*strd, 20.*strd);
             poly1 = RooPolynomial("linx", "linx", t, RooArgList(var1, var2), 1);
             // Construct linear (x) gauss
             Conv = new RooFFTConvPdf("linXg", "linear (X) gauss", t, poly1, gauss, 2);
            }
    else {
          std::cout << __FUNCTION__ << " ERROR: Unknown convolution type requested!" << std::endl;
          delete mag;
          return -4;
         }

    // S a m p l e ,   f i t   a n d   p l o t   c o n v o l u t e d   p d f 
    // ----------------------------------------------------------------------
    RooDataHist* data = new RooDataHist("dh", "dh", t, mag);
    // Fit gxlx to data
    RooFitResult* fitResult = Conv->fitTo(*data, Save(kTRUE), Range("range"));

    double par1 = -99.;
    double par1Err = -99.;
    double par2 = -99.;
    double par2Err = -99.;
    if (conv == "LinXGau")
       {
        par1 = var1.getVal();
        par1Err = var1.getError();
        par2 = var2.getVal();
        par2Err = var2.getError();
       }
    else if (conv == "LanXGau")
            { 
             par1 = var1.getVal();
             par1Err = var1.getError();
             par2 = var2.getVal();
             par2Err = var2.getError();
            }
    double meanG = mg.getVal();
    double meanGerror = mg.getError();
    double sigmaG = sg.getVal();
    double sigmaGerror = sg.getError();
    int quality = 0;

    if (m_verbose == 2) 
       {
        std::cout << __FUNCTION__ << " INFO: Results from the";
        if (conv == "LanXGau") std::cout << " Landau x Gauss fit" << std::endl << "MPV Landau: " << par1 << "+/-"
                                         << par1Err << " - sigma: " << par2 << "+/-" << par2Err << std::endl;
        else if (conv == "LinXGau")  std::cout << " poly x Gauss fit" << std::endl << "Linear slope: " << par2 << "+/-" << par2Err
                                               << " - intercept: " << par1 << "+/-" << par1Err << std::endl;
        std::cout << "mean Gauss: " << meanG << "+/-" << meanGerror << " - sigma: " << sigmaG << "+/-" << sigmaGerror << std::endl;
       }

    if (!fitResult)
       {
        std::cout << __FUNCTION__ << " WARNING: Fit result not there -> will not calculate magnitude!" << std::endl;
        quality = -5;
       }
    else {
          quality = (int)fitResult->covQual();
          if (m_verbose == 2) std::cout << __FUNCTION__ << " INFO: Quality of Landau x Gauss fit: " << quality << std::endl;
          if (quality < 1)
             {
              std::cout << __FUNCTION__ << " WARNING: Convolution quality is bad [" << quality << "] -> will not calculate magnitutde!" << std::endl;
              quality = -6;
             }
          else {
                if (conv == "LanXGau")
                   { 
                    magMPV = std::make_pair(par1 + meanG, sqrt(pow(par1Err, 2) + pow(meanGerror, 2)));
                    magSigma = std::make_pair(par2 + sigmaG, sqrt(pow(par2Err, 2) + pow(sigmaGerror, 2)));
                   }
                else if (conv == "LinXGau") 
                        {
                         magMPV = std::make_pair(meanG, meanGerror);
                         magSigma = std::make_pair(sigmaG, sigmaGerror);
                        }
                quality = 0;
               }
         }

    delete fitResult;
    delete data;
    delete Conv;
    delete mag;
    return quality;
}
// --------------------------------------------------------------------------------------------------------------
int LGADBase::LinearFit(std::vector<double>* vec, std::pair<double, double> &slope, std::pair<double, double> &intersept, std::vector<double>* vecErr)
{
    unsigned int pulses = vec->size();
    if (pulses < 10)
       {
        std::cout << __FUNCTION__ << " WARNING: Less than 60 points for this falculation["
                  << pulses << "] -> will not cintimue!" << std::endl;
        return -1;
       }

    double mean = Mean(vec);
    double strdv = Stdev(vec);
    double wmax = *std::max_element(vec->begin(), vec->end());
    double wmin = *std::min_element(vec->begin(), vec->end());
    double fact = pow(10, fabs(floor(log10(fabs(mean)))));
    double limit1 = floor(wmin * fact) / fact;
    double limit2 = ceil(wmax * fact) / fact;

    TGraphErrors *ped = new TGraphErrors();
    for (unsigned int p = 0; p < pulses; p++)
        { 
         ped->SetPoint(p, p, vec->at(p));
         if (vecErr != NULL) ped->SetPointError(p, 1/2, vecErr->at(p));
        }
    ped->SetMinimum(limit1);
    ped->SetMaximum(limit2);

    if (ped->Integral() == 0)
       {
        std::cout << __FUNCTION__ << " WARNING: Integral of histogram is 0 -> will not continue with linear fitting!" << std::endl;
        delete ped;
        return -2;
       }

    // if (m_verbose == 2) std::cout << __FUNCTION__ << " INFO: Calculating before fit| Min: " << wmin << " , Max: " << wmax 
    //                                  << ", Integral.: " << ped->Integral() << ", Linits: " << limit1 << " to " << limit2 << std::endl;
     
    TF1 *line = new TF1("line", "pol1", 0, pulses);
    line->SetParameter(0, mean);
    line->SetParameter(1, strdv);
    TFitResultPtr myfit = ped->Fit(line, "S Q N R", "+rob=0.75"); // 'q' to suppress messages  
    // if (m_verbose == 2) myfit->Print("V");

    int qual = 0;
    double chi2 = 0.0;
    if (!(myfit->IsValid()) || myfit->IsEmpty())
       {
        std::cout << __FUNCTION__ << ": WARNING : No fit result!" << std::endl;
        qual = -3;
       }
    else if (myfit->Ndf() == 0)
            {
             std::cout << __FUNCTION__ << ": WARNING : Zero degrees of freedom, fit unreliable!" << std::endl;
             qual = -4;
            }
    else {
          chi2 = fabs(1 - (myfit->MinFcnValue() / myfit->Ndf()));
          slope.first = myfit->Parameter(1);
          slope.second = myfit->ParError(1);
          intersept.first = myfit->Parameter(0);
          intersept.second = myfit->ParError(0);              
          }

    // if (m_verbose == 2) std::cout << __FUNCTION__ << " INFO: Fit quality " << qual << " and fit chi2/NDF: " << chi2 << ", Slope: "
    //                                 << myfit->Parameter(1) << " +/- " << myfit->ParError(1) << " , Intercept: " << myfit->Parameter(0) 
    //                                 << " +/- " << myfit->ParError(0) << std::endl;

    delete ped;
    delete line;
    return qual;
}
// --------------------------------------------------------------------------------------------------------------
double LGADBase::LinearInter(double x1, double y1, double x2, double y2, double y3)
{
    double slope, x3;
    if (x2 != x1) 
       {
        slope = (y2 - y1) / (x2 - x1);
        if (slope != 0) x3 = ((y3 - y1)/slope) + x1;
        else x3 = x1;
       }
    else x3 = x1;
    return x3;
}
// --------------------------------------------------------------------------------------------------------------
double LGADBase::FFT(vector<double> *w, Long64_t snrate, int start, int stop)
{
    if (start < 0) start = 0;
    if (stop <= 0) stop = w->size() - 1;
    if (start >= stop) { start = 0; stop = w->size()-1; }
    double freq = -1.;
    int npoints = (stop - start) + 1;
    double* kk = new double[npoints];
    for (int i = start; i <= stop; i++) kk[i-start] = w->at(i);
            
    if (m_verbose == 2) std::cout << __FUNCTION__ << " INFO : start indx: " << start << ", stop indx: " << stop << ", no. of points: " << npoints << " " << std::endl;
    TVirtualFFT *fftWav = TVirtualFFT::FFT(1, &npoints, "R2C");
    fftWav->SetPoints(kk);
    fftWav->Transform();
    TH1 *hm = 0;
    hm = TH1::TransformHisto(fftWav, hm, "MAG");
    hm->GetXaxis()->SetRange(0, npoints / 2);
    freq = (snrate/(double)(npoints-1))*hm->GetBinCenter(hm->GetMaximumBin());
    if (m_verbose == 2) std::cout << __FUNCTION__ << " INFO : estimated frequecy: " << freq  << ", total frequency: " << snrate/(double)(npoints-1) 
                                  << ", maxbin: " << hm->GetMaximumBin() << ", max bin center: " << hm->GetBinCenter(hm->GetMaximumBin()) << ", Max bin content: "
                                  << hm->GetBinContent(hm->GetMaximumBin()) << std::endl;
    delete fftWav;
    delete hm;
    delete[] kk; 

    if (freq != -1) return freq;
    else {
          std::cout << __FUNCTION__ << " ERROR: Failed to perform FFT transform between "
                    << start << " and " << stop << " points of the voltage vector!" << std::endl;
          return -99.;
         }
}
