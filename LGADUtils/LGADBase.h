/*
* LGADBase.h
*
*
*      Author: Gkougkousis Evangelos - Leonidas
*              egkougko@cern.ch
*               IFAE-BARCELONA
*/

#ifndef __LGADBase__
#define __LGADBase__

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm> 
#include <numeric>
#include <sys/stat.h>
#include <sys/types.h> 

#include "TSystemFile.h"
#include "TSystemDirectory.h"
#include "TFitResultPtr.h"
#include "TObject.h"
#include "TClass.h"
#include "TROOT.h"
#include "TFile.h"
#include "TKey.h"
#include "TSelector.h"
#include "TTree.h"
#include "TBranch.h"
#include "TIterator.h"
#include "TBits.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TString.h"
#include "RooFitResult.h"

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <Windows4Root.h>
#include <windows.h>
#else
#include <dirent.h>
#include <bits/stdc++.h> 
#endif

class LGADBase : public TSelector {

 public:

  enum AqInstrument { Sampic, LabTXT, InfiniiumScope, TektronixScope, LeCroyWRBin, Unasigned };
  
  bool m_convert;
  bool m_TrackComb;
  bool m_fei4Eff;
  bool m_WaveShape;
  int m_verbose;
  TString m_ofname;
  TString m_ofdir;
  TString m_TransFileName;  // Transimpedence Simulation File
  TFile* m_TransFile;
  TFile* m_ofile;
  TTree* m_tree;
  unsigned int m_event;
  unsigned int m_nchan;
  AqInstrument m_instrument;
  
  // nTuple Variables
  std::vector<double> m_scale;
  std::vector<unsigned int> m_npoints;
  std::vector<Long64_t> m_srate; 
  std::vector<double> m_triggTime;
  std::vector<double> m_physt;
  std::vector<double> m_ordrt;
  std::vector<unsigned int> m_channels;
  double  m_trigtime;
  std::string m_fitopt;

  // Constructor and destructor
  LGADBase();
  virtual ~LGADBase();
  
  // Set Methodes
  void SetInstrument(unsigned instr = 2);
  bool SetSRate(double rate, unsigned int ch);
  bool SetNPoints(unsigned int points, unsigned int ch);
  void SetInDataNames (TString DataDir = "", TString DataName = "", TString ext = "");
  void SetOutDataNames (TString DataDir = "", TString DataName = "");
  void SetTrackComb (bool comb = false);
  void SetFEi4Eff(bool FEi4Eff = false);
  void SetTrackInDataNames (TString DataDir = "", TString DataName = "");
  void SetVerbose(int verbose = 1) { m_verbose = verbose; };
  void SetFitMethode(std::string method);
  void SetConvertSucess(bool convert) { m_convert = convert;};
  bool SetRootFile(const char* file) { m_ofile = TFile::Open(file); return m_ofile->IsOpen(); };
  bool SetRootTree(TFile* f, std::string name = "");
  void SetWaveShape(bool doShape = false) { m_WaveShape = doShape; };

  // Get Methodes
  AqInstrument GetInstrument();
  double GetSRate(unsigned int ch);
  unsigned int GetNPoints(unsigned int ch);
  bool GetTrackComb();
  bool GetFEi4Eff();
  int GetVerbosity() { return m_verbose; };
  std::string GetFitMethode() { return m_fitopt; };
  unsigned int GetChNo() { return m_nchan; };
  std::vector<unsigned int> GetChannels() { return m_channels; };
  unsigned int GetEvntNo() { return m_event; };
  TString GetTransFileName() { return m_TransFileName; };
  TFile* GetTransFile() { return m_TransFile; };
  TFile* GetRootFile() { return m_ofile; };
  TTree* GetRootTree() { return m_tree; };
  bool GetConvertSucess() { return m_convert; };
  bool GetWaveShape() { return m_WaveShape; };
  TString GetExtention() { return m_ext; };
  TString GetDataDir() { return m_datadir; };
  TString GetDataName() { return m_dataname; };
  TString GetOutFileaDir() { return m_ofdir; };
  TString GetOutFileName() { return m_ofname; };

  // Class funcitons
  void SetStartStopEvnt(unsigned int Evnt1 = 0, unsigned int Evnt2 = 0);
  bool ConvertData();
  void Initialize();
  void SetVectorSize(unsigned int ch);

  // General helper functions
  std::string reduce(const std::string str, const std::string fill, const std::string whitespace);
  std::string trim(const std::string str, const std::string whitespace);
  unsigned int CountFiles(const char* dir, const char* ext);
  int DirExists(const char* path);
  int CreateDir(const char* path);
  int RecursMkDir(const char* dirname);
  std::vector<std::string> ListFileNames(const char* path, const char* ext);
  template <typename T> T Derivate(T *w, int start = 1);
  template <typename U> double Mean(U *w, int start = -1, int stop = -1);
  template <typename V> double Stdev(V *w, int start = -1, int stop = -1, double mean = 0.0);
  template <typename V, typename T> double BayesianErr(V *w, T value = 1);
  template <typename V> double CalcMeadian(V *vec);
  void ProgressBar(Long64_t evt, Long64_t total);

  // LGAD Fits
  int IterativeFit(std::vector<double> *w, std::pair<double, double> &gmean, std::pair<double, double> &gsigma, TH1D* &FitHist,
                   double &minchi2, std::string methode = "Gauss", std::pair<int, int> points = std::make_pair(-1, -1));
  int LinearFit(std::vector<double>* vec, std::pair<double, double> &slope, std::pair<double, double> &intersept, std::vector<double>* vecErr = NULL);
  int RooConvFit(std::vector<double>* vec, std::pair<double, double> &magMPV, std::pair<double, double> &magSigma, std::string conv);
  double LinearInter(double x1, double y1, double x2, double y2, double y3);
  double FFT(std::vector<double> *w, Long64_t snrate, int start, int stop);

 private:

  // Member variables
  TString m_dataname;
  TString m_datadir;
  TString m_Trackdatadir;
  TString m_Trackdataname;
  TString m_ext;
  unsigned int m_evnt1;
  unsigned int m_evnt2;
  unsigned int m_bins; // number of bins for fitting histos

  // Ntuple Variables
  std::vector<std::vector<double> > m_t;
  std::vector<std::vector<double> > m_w;
  TH1F* m_trigDt;

  bool WriteSampic(const char* dir, const char* name, const char* ext);
  bool WriteLecroyTXT(const char* dir, const char* name, const char* ext, unsigned int evt1 = 0, unsigned int evt2 = 0);
  bool WriteTectronixTXT(const char* dir, const char* name, const char* ext);
  bool WriteAgilentBinary(const char* dir, const char* name, const char* ext, unsigned int evt1 = 0, unsigned int evt2 = 0);
  bool WriteLecroyBinary(const char* dir, const char* name, const char* ext, unsigned int evt1 = 0, unsigned int evt2 = 0);
  bool CombineTrack(const char* dir, const char* name);
  bool CreateOutputFile(const char* dir, const char* ofname, std::vector<unsigned int> nchan);
  bool SetScale(std::vector<unsigned int> channel, unsigned int nchan, std::vector<double>* scale);

  // Internal functions for iterative re-fitting
  TFitResultPtr Gauss(double rmin, double rmax, double strdv, double mean, TH1D* &Magnitude, std::string integ = "N");
  TFitResultPtr GauXLandau(double rmin, double rmax, double strdv, TH1D* &Magnitude, std::string integ = "N");
  Double_t LandXGauFun(Double_t *x, Double_t *par);
  int langaupro(TFitResultPtr fitResult, double &maxx, double &FWHM);

  // ClassDef(LGADBase,0);

};

#endif