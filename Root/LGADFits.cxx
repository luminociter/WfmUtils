/*
* LGADFits.cxx
*
*
*      Author: Gkougkousis Evangelos - Leonidas
*              egkougko@cern.ch
*               IFAE-BARCELONA
*/

#include "../LGADUtils/LGADBase.h"
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
//#include "TCanvas.h"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#ifdef __ROOTCLING__
#pragma "RooGlobalFunc.h"
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

int LGADBase::IterativeFit(std::vector<double> *w, std::pair <double, double> &gmean, std::pair <double, double> &gsigma, TH1D* &FitHist,
                           double &minchi2, std::string methode, std::pair <int, int> points, bool discr)
{
    std::vector<double> w2; // Temporary container for used element subset per itereation
    w2.reserve(w->size());
    gmean = std::make_pair(-1., -1.);
    gsigma = std::make_pair(-1., -1.);
    minchi2 = -1.;
    int last, first;
    if (points.second < 0) last = w->size();
    else last = points.second;
    if (points.first < 0) first = 0;
    else first = points.first;

    if (((methode == "Gauss" ||methode == "GaussVarBin" || methode == "GaussInt" || methode == "GaussIntVarBin" ) && abs(last - first) < 5) 
        || ((methode == "LandauXGauss" || methode == "LandauXGaussVarBin" || methode == "LandauXGaussInt" || methode == "LandauXGaussIntVarBin") && abs(last - first) < 10))
       {
        if (LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Inadequate number of points for " << methode << " calculation["
                                     << abs(last - first) << "] -> will not cuntinue!" << std::endl;
        return -1;
       }

    // Outlier rejection
    w2 = LGADBase::OutlierReject(w, 3, 0.2, first, last);
    
    // Calculate dataset mean, std, min and max value for limits and binning variations
    double mean = LGADBase::Mean(&w2);
    double strdv = LGADBase::Stdev(&w2);
    double wmax = *std::max_element(w2.begin(), w2.end());
    double wmin = *std::min_element(w2.begin(), w2.end());
    if (mean == -99 || strdv == -99 || strdv == 0)
       {
        if (m_verbose > 0) std::cout << __FUNCTION__ << " WARNING: Unabale to determine mean and stdv ["
                                     << mean << ", " << strdv << "] -> will not cintimue!" << std::endl;
        return -2;
       }

    double median = -99.;
    // Calculate minimum element distance to define minimum bin width
    if (methode == "GaussVarBin" || methode == "LandauXGaussVarBin" || methode == "GaussIntVarBin" || methode == "LandauXGaussIntVarBin") discr = true;
    double acc = -1;
    double FWHM = -1;
    if (discr) 
       {
        acc = LGADBase::CalResolution<double>(&w2, 5);
        if (methode == "LandauXGaussVarBin" || methode == "LandauXGaussIntVarBin") median = LGADBase::MaxDensity<double>(&w2, acc);
       }
    else if (methode == "LandauXGauss" || methode == "LandauXGaussVarBin" || methode == "LandauXGaussInt" || methode == "LandauXGaussIntVarBin")
            { 
             median = LGADBase::MaxDensity<double>(&w2);
             FWHM = LGADBase::CalcFWHM(&w2, median);
            }

    // Histogtam bins and limits
    double limit1 = 0.0;
    double limit2 = 0.0;
    double chi2 = 0.0;
    double gdnes = 0.0;
    bool stop = false;
    // Variables for fit range
    double rmin = 0.0;
    double rmax = 0.0;
    int quality = -99;
    std::vector<double> mag;
    mag.reserve(42);
    std::vector<double> magErr;
    magErr.reserve(42);
    std::vector<double> Sigma;
    Sigma.reserve(42);
    std::vector<double> SigmaErr;
    SigmaErr.reserve(42);
    std::vector<double> ChiSq;
    ChiSq.reserve(42);
    std::vector<double> GdNess;
    GdNess.reserve(42);
    std::vector<int> globindx;
    globindx.reserve(42);
    std::vector<TH1D* > Fits;
    Fits.resize(42, NULL);
    std::vector<unsigned int> itr;
    itr.reserve(42);
    TAxis *xaxis;
    Int_t bin1 = 0;
    Int_t bin2 = 0;
    int n_bins[7] = {0};

    std::vector<double> VarBinX; // vector only necessry if rebining 

    if (m_verbose >= 2) 
       {
        std::cout << __FUNCTION__ << " INFO: Calculating before fit | Average: " << mean << ", Standard Deviation: " 
                                  << strdv << ", Median: " << median << ", Full width half maximum: " << FWHM << std::endl;
        std::cout << __FUNCTION__ << " INFO: Fitting options | Methode: " << methode << ", Discrete: " << discr << ", Outlieres: " 
                                  << w->size() - w2.size() << ", fit elements vector size: " << w->size() << ", first - last element: "
                                  << first << "-" << last << std::endl;
       }

    double fact = pow(10, fabs(floor(log10(fabs(strdv))))); // Rounding factor for limits
    int bins_max = 0; // Absolute maximum number of bins extrapolated from minimum resolution and limtis
    int n_elements = 0; // Number of points included inth e fit
    double strdv2 = -99; // standard deviation of modified dataset
    double mean2 = -99; // mean value of modified dataset
    double median2 = -99; // median of modified dataset
    double FWHM2 = -99; // full width half maximum of modified dataset
    std::vector<double> wb;
    std::vector<double> w2b;
    // Iterative re-fitting 1, affects the limits of the histogram only
    for (double gk = 12; gk > 6; gk--) // 6 cases
        {
         if (stop) break;
         if (methode == "Gauss" || methode == "GaussVarBin" || methode == "GaussInt" || methode == "GaussIntVarBin")
            {
             limit1 = floor((mean - (gk/3)*strdv)*fact) / fact; // Symentric limit variaton
             limit2 = ceil((mean + (gk/3)*strdv)*fact) / fact;
            }
         else if (methode == "LandauXGauss" || methode == "LandauXGaussVarBin" || methode == "LandauXGaussInt" || methode == "LandauXGaussIntVarBin")
                 {
                  limit1 = floor((median - ((gk-5)/7)*FWHM)*fact) / fact; // Asymentric limit variaton
                  limit2 = ceil((median + ((gk-2)/4)*FWHM)*fact) / fact;
                 }
         wb.clear();
         w2b.clear();
         for (int gb = first; gb < last; gb++) if (w->at(gb) >= limit1 && w->at(gb) <= limit2) wb.push_back(w->at(gb));
         for (unsigned int gb = 0; gb < w2.size(); gb++) if (w2.at(gb) >= limit1 && w2.at(gb) <= limit2) w2b.push_back(w2.at(gb));
         n_elements = wb.size();
         if (n_elements < 5) continue; // minimum number of elements check
         // Compute numbers for modified datasheet
         strdv2 = LGADBase::Stdev(&w2b);
         mean2 = LGADBase::Mean(&w2b);
         // Calculate maximum number of bins that makes sense for discreete histograms
         if (discr) 
            {
             bins_max = (int)(fabs(limit2 - limit1) / acc);
             if (methode == "LandauXGaussVarBin" || methode == "LandauXGaussIntVarBin") median2 = LGADBase::MaxDensity<double>(&w2b, acc);
            }
         else if (methode == "LandauXGauss" || methode == "LandauXGaussVarBin" || methode == "LandauXGaussInt" || methode == "LandauXGaussIntVarBin")
                 {
                  median2 = LGADBase::MaxDensity<double>(&w2b);
                  FWHM2 = LGADBase::CalcFWHM(&w2b, median2);
                 }
         // Calculate bin vector for variable bin histograms
         if (methode == "GaussVarBin" || methode == "LandauXGaussVarBin" || methode == "GaussIntVarBin" || methode == "LandauXGaussIntVarBin")
            {
             VarBinX = LGADBase::ConrtVarBinX(&wb, limit1, limit2, bins_max);
            }
         LGADBase::CalcuRebin(discr, n_elements, bins_max, limit1, limit2, strdv2, n_bins);
         // Iterative re-fitting 2, affects only the width of the bins but not the fit range
         bool varbinbase = false;
         for (unsigned int tmb = 0; tmb < 7; tmb++)
             {
              if (stop) break;
              if (n_bins[tmb] == -1 || n_bins[tmb] == 0) continue;
              unsigned int iter = (12 - gk) * 7 + tmb;
              // Variable bin algorithm
              if (methode == "GaussVarBin" || methode == "LandauXGaussVarBin" || methode == "GaussIntVarBin" || methode == "LandauXGaussIntVarBin")
                 {
                  std::vector<double> ReBinXbins;
                  ReBinXbins.clear();
                  ReBinXbins.reserve(n_bins[tmb]);
                  ReBinXbins.push_back(VarBinX.at(0));
                  for (int pr = 1; pr < bins_max; pr++)
                      {
                       double diff = fabs(VarBinX.at(pr) - VarBinX.at(pr-1));
                       //if (m_verbose >= 2) std::cout << pr << "/" << bins_max << ", bins: " << n_bins[tmb] << ", factor " << (float)bins_max/(float)n_bins[tmb] << ", diff: " << diff << ", mod.diff: " << diff*((float)bins_max / (float)n_bins[tmb]) 
                       //                              << ", previous bin: " << ReBinXbins.at(pr - 1) << ", next mod.bin: " << ReBinXbins.at(pr - 1) + diff << ", bin upper lim:" << VarBinX.back() << std::endl;
                       diff = diff*((float)bins_max/(float)n_bins[tmb]);
                       if ((ReBinXbins.at(pr-1) + diff) <= VarBinX.back()) ReBinXbins.push_back(ReBinXbins.at(pr-1) + diff);
                       else {
                             ReBinXbins.push_back(VarBinX.back());
                             break;
                            }
                      }
                  // if (m_verbose >= 2) std::cout << tmb << " " << ReBinXbins.size() << " " << bins_max << std::endl;
                  if (varbinbase) continue;
                  if ((int)ReBinXbins.size() == bins_max) varbinbase = true;
                  m_bins = ReBinXbins.size();
                  Fits.at(iter) = new TH1D(Form("MagFit%02u", iter), Form("MagFit%02u", iter), m_bins - 1, &(ReBinXbins[0]));
                 }
              // Normal binning algorithm
              else {
                    m_bins = n_bins[tmb];
                    Fits.at(iter) = new TH1D(Form("MagFit%02u", iter), Form("MagFit%02u", iter), m_bins, limit1, limit2);
                   }
              // Create the re-binned histogram and calculate quanitites
              for (int i = 0; i < n_elements; i++) Fits.at(iter)->Fill(wb.at(i));
              itr.push_back(iter);
              if (Fits.at(iter)->Integral("width") == 0)
                 {
                  if (m_verbose > 0) std::cout << __FUNCTION__ << " WARNING: Integral of " << methode 
                                               << " distribution is 0 -> will not continue!" << std::endl;
                  continue;
                 }
              double mean3 = Fits.at(iter)->GetMean();
              double strdv3 = Fits.at(iter)->GetRMS();
              double median3 = (Fits.at(iter))->GetXaxis()->GetBinCenter((Fits.at(iter))->GetMaximumBin());
              double FWHM3 = (Fits.at(iter))->GetXaxis()->GetBinCenter((Fits.at(iter))->FindLastBinAbove(0.5*(Fits.at(iter))->GetBinContent((Fits.at(iter))->GetMaximumBin())));
              FWHM3 = FWHM3 - (Fits.at(iter))->GetXaxis()->GetBinCenter((Fits.at(iter))->FindFirstBinAbove(0.5 * (Fits.at(iter))->GetBinContent((Fits.at(iter))->GetMaximumBin())));
              if (strdv3 == 0)
                 {
                  if (m_verbose > 0) std::cout << __FUNCTION__ << " WARNING: Standard deviation of " << methode 
                                               << " distribution is 0 -> will not continue!" << std::endl;
                  continue;
                 }
              if (m_verbose >= 3) std::cout << __FUNCTION__ << " INFO: Magnitude histo mean: " << mean3 << " , stdev: " << strdv3 << ", integral: "
                                            << Fits.at(iter)->Integral("width") << ", bins: " << m_bins << ", from " << limit1 << " to " << limit2 << std::endl;
              TFitResultPtr fitResult(0);
              xaxis = Fits.at(iter)->GetXaxis();
              // Perform the fit and get the parameters
              if (methode == "Gauss" || methode == "GaussVarBin" || methode == "GaussInt" || methode == "GaussIntVarBin")
                 {
                  rmin = mean3 - (gk/3) * strdv3; // symetric range
                  rmax = mean3 + (gk/3) * strdv3;
                  bin1 = xaxis->FindBin(rmin);
                  bin2 = xaxis->FindBin(rmax);
                  if (((bin2-bin1) < 2) || (Fits.at(iter)->Integral(bin1 + 1, bin2 - 1) == 0)) continue;
                  if (methode == "GaussInt" || methode == "GaussIntVarBin") fitResult = LGADBase::Gauss(rmin, rmax, strdv2, mean2, Fits.at(iter), "I");
                  else fitResult = LGADBase::Gauss(rmin, rmax, strdv2, mean2, Fits.at(iter));
                 }
              else if (methode == "LandauXGauss" || methode == "LandauXGaussInt" || methode == "LandauXGaussIntVarBin" ||  methode == "LandauXGaussVarBin")
                      {
                       rmin = median3 - ((gk-5)/7)*FWHM3; // asymentric range
                       rmax = median3 + ((gk-2)/4)*FWHM3;
                       // std::cout << "Median: " << (fabs(median3-median2)/median2)*100 << "%, FWHM: " << (fabs(FWHM3-FWHM2)/FWHM2)*100  << "%, Lower Limit: " << (fabs(rmin-limit1)/limit1)*100  << "%,  Upper Limit: " << (fabs(rmax-limit2)/limit2)*100 << "%" << std::endl;
                       bin1 = xaxis->FindBin(rmin);
                       bin2 = xaxis->FindBin(rmax);
                       if (((bin2 - bin1) < 2) || (Fits.at(iter)->Integral(bin1 + 1, bin2 - 1) == 0)) continue;
                       if (methode == "LandauXGaussInt" || methode == "LandauXGaussIntVarBin") fitResult = LGADBase::GauXLandau(rmin, rmax, median2, FWHM2, Fits.at(iter), "I");
                       else fitResult = LGADBase::GauXLandau(rmin, rmax, median2, FWHM2, Fits.at(iter));
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
              // Fill the fit parameters to the vectors
              if (fitResult->IsEmpty() || !(fitResult->IsValid()))
                 {
                  if (m_verbose > 1) std::cout << __FUNCTION__ << " WARNING: Fit result not there -> will not calculate magnitude!" << std::endl;
                  quality = -4;
                 }
              else {
                    quality = fitResult;
                    if (quality == 0 && (fitResult->Ndf()) != 0)
                       {
                        chi2 = fabs((fitResult->MinFcnValue() - (float)fitResult->Ndf())/sqrt(2*(float)fitResult->Ndf()));
                        gdnes = fabs(1 - (fitResult->MinFcnValue()/(float)fitResult->Ndf()));
                        if (methode == "Gauss" || methode == "GaussVarBin" || methode == "GaussInt" || methode == "GaussIntVarBin")
                           {
                            mag.push_back(fitResult->Parameter(1));
                            magErr.push_back(fitResult->ParError(1));
                            Sigma.push_back(fitResult->Parameter(2));
                            SigmaErr.push_back(fitResult->ParError(2));
                            ChiSq.push_back(chi2);
                            GdNess.push_back(gdnes);
                            globindx.push_back(iter);
                            // if (gdnes < 0.03) stop = true;
                            if (m_verbose >= 2) std::cout << __FUNCTION__ << " INFO: Fit goodness: " << gdnes << ", magnitude: " << fitResult->Parameter(1)  << ", sigma: " 
                                                          << fitResult->Parameter(2) << " at iteration: " << iter << " from combinations " << gk << " & " << tmb << std::endl;
                           }
                        else if (methode == "LandauXGauss" || methode == "LandauXGaussInt" || methode == "LandauXGaussIntVarBin" ||  methode == "LandauXGaussVarBin")
                                {
                                 double max = -99;
                                 double FWHM4 = -99;
                                 if (LGADBase::langaupro(fitResult, max, FWHM4) == 0)
                                    {
                                     mag.push_back(max);
                                     magErr.push_back(fitResult->ParError(1));
                                     Sigma.push_back(FWHM4/2);
                                     SigmaErr.push_back(sqrt(pow(2 * fitResult->ParError(0), 2) + pow(2 * fitResult->ParError(3), 2)));
                                     ChiSq.push_back(chi2);
                                     GdNess.push_back(gdnes);
                                     globindx.push_back(iter);
                                     // if (gdnes < 0.03) stop = true;
                                     if (m_verbose >= 2) std::cout << __FUNCTION__ << " INFO: Fit goodness: " << gdnes << ", magnitude: " << max << ", sigma: "
                                                                   << FWHM4/2 << " at iteration: " << iter << " from combinations " << gk << " & " << tmb << std::endl;
                                    }
                                 else { if (m_verbose > 1) std::cout << __FUNCTION__ << " WARNING: LandauXGauss renormalization failed!" << std::endl; }
                                }
                       }
                     else { if (m_verbose > 1) std::cout << __FUNCTION__ << " WARNING: Bad fit quality (" << quality << ") or problematic Ndf (Ndf = " << fitResult->Ndf() << ")" << std::endl; }
                    }
              if (m_verbose >= 2) std::cout << __FUNCTION__ << " INFO: Quality of performed fit: " << quality << " at iteration: " 
                                            << iter << " from combinations " << gk << " & " << tmb << std::endl;
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
          minchi2 = GdNess.at(indx);
          if (m_verbose >= 2) std::cout << __FUNCTION__ << " INFO: Fit result goodness: " << GdNess.at(indx) << " --> " << mag.at(indx) << " +/- "
                                       << magErr.at(indx) << " , " << Sigma.at(indx) << " +/- " << SigmaErr.at(indx) << " at iteration " << indx << std::endl;
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
    // function to fit: [0]*exp(-0.5*((x-[1])/[2])**2)
    mygau->SetParNames("Normalization", "Mean", "Sigma");
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
    // Fit within specified range (R), save (S), do not print fitting messages (Q)
    if (integ == "I") fitResult = magHist->Fit(mygau, "S Q R I"); // Integral (I) option
    else fitResult = magHist->Fit(mygau, "S Q R");
    if (m_verbose >= 3) fitResult->Print("V");
    delete mygau;
    return fitResult;
}
// --------------------------------------------------------------------------------------------------------------
TFitResultPtr LGADBase::GauXLandau(double rmin, double rmax, double median, double FWHM, TH1D* &magHist, std::string integ)
{
    TF1 *mygauland = new TF1("LandXGauFun", this, &LGADBase::LandXGauFun, rmin, rmax, 4,"LGADUtils", "LandXGauFun");
    mygauland->SetParNames("Width", "MPV", "Area", "GaussSigma");
    TAxis *axis = magHist->GetXaxis();
    double wmin = axis->GetBinCenter(magHist->FindFirstBinAbove(0, 1));
    if (m_verbose >= 3) std::cout << __FUNCTION__ << " INFO: First populated bin " << magHist->FindFirstBinAbove(0, 1) << ", last populated bin: " << magHist->FindLastBinAbove(0, 1)
                                                  << ", maximum in bin " << magHist->GetMaximumBin() << " of " << magHist->GetSize() << " total bins" << std::endl;
    double wmax = axis->GetBinCenter(magHist->FindLastBinAbove(0, 1)); //abs
    double integral = magHist->Integral("width");
    double limitdown = median - 0.5 * FWHM;
    double limitup = median + 1.4 * FWHM;
    if (limitdown < rmin) limitdown = rmin;
    if (limitup > rmax) limitup = rmax;
    if (median <= limitdown || median >= limitup)
       {
        double a = fabs(median - limitdown);
        double b = fabs(limitup - median);
        if (b < a) mygauland->SetParameter(1, (b*(limitup - limitdown) / ((limitup - limitdown) + 2 * b)) + (limitdown + ((limitup - limitdown) / 2)));
        else mygauland->SetParameter(1, (limitdown + ((limitup - limitdown) / 2)) - (a*(limitup - limitdown) / ((limitup - limitdown) + 2 * a)));
       }
    else mygauland->SetParameter(1, median);
    mygauland->SetParameter(0, (wmax - wmin)/10);
    mygauland->SetParameter(2, integral);
    mygauland->SetParameter(3, FWHM /10);
    mygauland->SetParLimits(0, (wmax-wmin)/10, 10*(wmax-wmin));
    mygauland->SetParLimits(1, rmin, rmax);
    mygauland->SetParLimits(2, 0.1*integral, 10*integral);
    mygauland->SetParLimits(3, 0.01*(FWHM / 2), 100 * (FWHM / 2));
    //mygauland->FixParameter(3, strdv);

    // Fit within specified range (R), save (S), quiet mode (Q)
    TFitResultPtr fitResult;
    if (integ == "I") fitResult = magHist->Fit(mygauland, "S Q R I");
    else fitResult = magHist->Fit(mygauland, "S Q R");
    if (m_verbose >= 3) fitResult->Print("V");
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
    Double_t invsq2pi = 0.3989422804014; // (2 pi)^(-1/2)

    // Control constants
    double sc = 5.0; // convolution extends to +-sc Gaussian sigmas

    // MP shift correction
    double mpc = par[1] + 0.22278298 * par[0];

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
int LGADBase::RooConvFit(std::vector<double>* vec, std::pair <double, double> &magMPV, std::pair <double, double> &magSigma, std::string conv)
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
    if (m_verbose >= 2) std::cout << __FUNCTION__ << " INFO: Calculating mean and std before fit: " << mean << " , " << strd << std::endl;

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
    RooLandau landau("lx", "lx", t, var1, var2);
    RooPolynomial poly1("linx", "linx", t, RooArgList(var1, var2), 1);

    if (conv == "LanXGau")
       { 
        // Construct landau(t,ml,sl)
        var1.SetNameTitle("ml", "mean landau"); var1.setMin(-4*fabs(mean)); var1.setMax(2*fabs(mean)); var1.setVal(mean);
        var2.SetNameTitle("sl", "sigma landau"); var2.setMin(0.001*strd); var2.setMax(2*strd); var2.setVal(strd);
        // Construct landau (x) gauss
        Conv = new RooFFTConvPdf("lxg", "landau (X) gauss", t, landau, gauss, 2);
       }
    else if (conv == "LinXGau")
            {
             // Construct linear(t,slope,interc)
             var1.SetNameTitle("ilin", "intersect linear"); var1.setMin(mean-3*strd); var1.setMax(mean+3*strd); var1.setVal(mean);
             var2.SetNameTitle("slin", "slope linear"); var2.setMin(0.01*strd); var2.setMax(20.*strd); var2.setVal(strd);
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

    if (m_verbose >= 2) 
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
          if (m_verbose >= 2) std::cout << __FUNCTION__ << " INFO: Quality of Landau x Gauss fit: " << quality << std::endl;
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
int LGADBase::LinearFit(std::vector<double>* vec, std::pair <double, double> &slope, std::pair <double, double> &intersept, std::vector<double>* vecErr)
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

    // if (m_verbose >= 2) std::cout << __FUNCTION__ << " INFO: Calculating before fit| Min: " << wmin << " , Max: " << wmax 
    //                                  << ", Integral.: " << ped->Integral() << ", Linits: " << limit1 << " to " << limit2 << std::endl;
     
    TF1 *line = new TF1("line", "pol1", 0, pulses);
    line->SetParameter(0, mean);
    line->SetParameter(1, 0);
    TFitResultPtr myfit = ped->Fit(line, "S Q N E", "+rob=0.75"); // 'q' to suppress messages  
    // if (m_verbose >= 2) myfit->Print("V");

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

    // if (m_verbose >= 2) std::cout << __FUNCTION__ << " INFO: Fit quality " << qual << " and fit chi2/NDF: " << chi2 << ", Slope: "
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
double LGADBase::FFT(std::vector<double> *w, Long64_t snrate, int start, int stop)
{
    if (start < 0) start = 0;
    if (stop <= 0) stop = w->size() - 1;
    if (start >= stop) { start = 0; stop = w->size()-1; }
    double freq = -1.;
    int npoints = (stop - start) + 1;
    double* kk = new double[npoints];
    for (int i = start; i <= stop; i++) kk[i-start] = w->at(i);
            
    if (m_verbose >= 2) std::cout << __FUNCTION__ << " INFO : start indx: " << start << ", stop indx: " << stop << ", no. of points: " << npoints << " " << std::endl;
    TVirtualFFT *fftWav = TVirtualFFT::FFT(1, &npoints, "R2C");
    if (strcmp (fftWav->GetDefaultFFT(), "") == 0)
       {
        if (m_verbose >= 2) std::cout << __FUNCTION__ << " WARNING : FFT libraryu not installed, no Furrier analys performed!" << std::endl;
        return  -99.;
       }
    fftWav->SetPoints(kk);
    fftWav->Transform();
    TH1 *hm = 0;
    hm = TH1::TransformHisto(fftWav, hm, "MAG");
    hm->GetXaxis()->SetRange(0, npoints / 2);
    freq = (snrate/(double)(npoints-1))*hm->GetBinCenter(hm->GetMaximumBin());
    if (m_verbose >= 2) std::cout << __FUNCTION__ << " INFO : estimated frequecy: " << freq  << ", total frequency: " << snrate/(double)(npoints-1) 
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
// --------------------------------------------------------------------------------------------------------------
std::vector<double> LGADBase::ConrtVarBinX(std::vector<double> *wmod, double limUp, double limDown, int &nbins)
{
    std::vector<double> ReBinX;
    ReBinX.reserve(nbins);
    // Find minimum bin width and create temporary histo to rebin
    TH1D* FTempReBin = new TH1D("FTempReBin", "FTempReBin", nbins + 1, limUp, limDown);
    for (unsigned int ga = 0; ga < wmod->size(); ga++) FTempReBin->Fill(wmod->at(ga));
    TAxis *TmpRbinXAxis = FTempReBin->GetXaxis();
    for (int i = 0; i < nbins; i++)
        {
         for (int j = i + 1; j < nbins + 1; j++)
             {
              if (FTempReBin->GetBinContent(j) != 0)
                 {
                  ReBinX.push_back(TmpRbinXAxis->GetBinLowEdge(j));
                  i = j - 1;
                  break;
                 }
              else {
                    if (j < nbins) continue;
                    else {
                          ReBinX.push_back(TmpRbinXAxis->GetBinLowEdge(j + 1));
                          i = nbins + 1;
                         }
                   }
             }
        }
    delete FTempReBin;
    nbins = ReBinX.size(); // Modify maximum amount of bins for variable binning hisots
    return ReBinX;
}
// --------------------------------------------------------------------------------------------------------------
bool LGADBase::CalcuRebin(bool discr, int n_elements, int nbins, double limDown, double limUp, double stdev, int (&Nofbins)[7])
{
    // Rebinning Algorithm
    for (unsigned int u = 0; u < 7; u++) Nofbins[u] = -1;
    int varlow = 0;
    int varhigh = 0;
    if ((discr && sqrt(n_elements) < nbins) || !discr)
       {
        if (fabs((limDown - limUp)/stdev) < sqrt(n_elements) || !discr)
           {
            unsigned int dfact = 1;
            if (fabs((limDown - limUp)/stdev) < sqrt(n_elements))
               {
                if (discr) // Define the upper part of the bin arrray for descrete datasets
                   {
                    varhigh = floor((float)(fabs(nbins - sqrt(n_elements))) / 3);
                    // Reduction of number of bins at high statistics mod
                    if (nbins > 10*ceil(sqrt(n_elements))) dfact = pow(10, floor(fabs(log10(nbins/ceil(sqrt(n_elements))))))/2;
                    else if (nbins > 5*ceil(sqrt(n_elements))) dfact = 5/2;
                    if (m_verbose >= 2) std::cout << "Case 1 discr " << dfact << std::endl;
                   }
                else {
                      if (m_verbose >= 2) std::cout << "Case 2, !discr" << std::endl;
                      varhigh = floor((float)fabs((limDown - limUp) / stdev) / 3);
                     }
                varlow = ceil((float)(sqrt(n_elements) - (fabs(limDown - limUp)/stdev))/3);
                Nofbins[3] = ceil(sqrt(n_elements));
               }
            else if (fabs((limDown - limUp)/stdev) >= sqrt(n_elements) && !discr)
                    {
                     Nofbins[3] = ceil(fabs(limDown - limUp)/stdev);
                     varlow = ceil((float)((fabs(limDown - limUp)/stdev)-sqrt(n_elements))/3);
                     varhigh = floor(sqrt(n_elements)/3);
                     if (m_verbose >= 2) std::cout << "Case 3, !discr" << std::endl;
                    }
            if (m_verbose >= 2)  std::cout << "Numbers: " << varlow << " " << varhigh << " " << sqrt(n_elements) << " " << n_elements << " " << fabs((limDown - limUp) / stdev) << " " << stdev << " " << limDown << " " << limUp << std::endl;
            for (unsigned int gkl = 0; gkl < 3; gkl++)
                {
                 int bin1 = ceil(Nofbins[3] - (gkl+1)*varlow);
                 int bin2 = floor((float)(Nofbins[3] + (gkl + 1)*varhigh) / (float)dfact);
                 // Checking the array for existing values
                 bool exist1 = false;
                 bool exist2 = false;
                 for (unsigned amr = 0; amr < gkl+1; amr++) 
                     {
                      if (dfact == 1)
                         {
                          if (!exist1) { if (Nofbins[3 - amr] == bin1) exist1 = true; }
                          if (!exist2) { if (Nofbins[3 + amr] == bin2) exist2 = true; } 
                         }
                      else {
                            if (bin2 == Nofbins[3+amr] || bin2 == Nofbins[3-amr])
                               {
                                bool acrt = true;
                                do { 
                                    acrt = false;
                                    bin2++;
                                    for (unsigned amb = 0; amb < gkl + 1; amb++) 
                                        {
                                         if (Nofbins[3 - amb] == bin2 || Nofbins[3 + amb] == bin2) 
                                            {
                                             acrt = true; 
                                             break; 
                                            }
                                        }
                                   }
                                while (acrt && bin2 < nbins);
                                if (bin2 > nbins) exist2 = true;
                               }
                            if (bin1 == Nofbins[3+amr] || bin2 == Nofbins[3-amr])
                               {
                                bool acrt = true;
                                do { 
                                    acrt = false;
                                    bin1--;
                                    for (unsigned amb = 0; amb < gkl + 1; amb++) 
                                        {
                                         if (Nofbins[3 - amb] == bin1 || Nofbins[3 + amb] == bin1 || (!exist2 && bin1 == bin2)) 
                                            {
                                             acrt = true; 
                                             break; 
                                            }
                                        }
                                  }
                                while (acrt && bin1 > 0);
                                if (bin1 <= 0) exist1 = true;
                               }
                           }
                      if (exist1 && exist2) break;
                     }
                 if (!exist1) Nofbins[2 - gkl] = bin1;
                 else Nofbins[2 - gkl] = -1;
                 if (!exist2) Nofbins[4 + gkl] = bin2;
                 else Nofbins[4 + gkl] = -1;
                }
           }
        else { 
              if (m_verbose >= 2) std::cout << "Case 4, dicr" << std::endl;
              varhigh = floor((float)(fabs(nbins-sqrt(n_elements)))/6);
              if (varhigh == 0) varhigh = ceil((float)(fabs(nbins - sqrt(n_elements)))/6);
              Nofbins[0] = floor(sqrt(n_elements));
              for (unsigned int gkl = 1; gkl < 7; gkl++) 
                  {
                   float dfact = 1;
                   int bin = floor(sqrt(n_elements) + gkl*varhigh);
                   if (bin > n_elements) dfact = (float)n_elements/(float)bin;
                   bin = floor(bin*dfact);
                   bool exist = false;
                   for (unsigned int kh = 0; kh < gkl; kh++) if (Nofbins[kh] == bin) { exist = true; break; }
                   if (!exist && bin <= nbins) Nofbins[gkl] = bin;
                   else Nofbins[gkl] = -1;
                  }
             }
      }
    else {
          if (m_verbose >= 2) std::cout << "Case 5, dicr" << std::endl;
          varhigh = floor((float)nbins/7);
          if (varhigh == 0) varhigh = ceil((float)nbins/7);
          for (unsigned int gkl = 1; gkl < 8; gkl++) 
              {
               int bin = varhigh*gkl;
               bool exist = false;
               for (unsigned int kh = 0; kh < gkl-1; kh++) if (Nofbins[kh] == bin) { exist = true; break; }
               if (!exist && bin <= nbins) Nofbins[gkl-1] = bin;
               else Nofbins[gkl-1] = -1;
              }
         }

    // Printout calculated bins for each case
    if (m_verbose >= 2)
       {
        std::cout << "Bin array: ";
        for (unsigned int kh = 0; kh < 7; kh++) 
            {
             std::cout << Nofbins[kh];
             if (kh < 6) std::cout << ":";
            }
        std::cout << std::endl;
       }

    return true;
}