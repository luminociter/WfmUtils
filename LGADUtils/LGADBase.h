/*
* LGADBase.h
*
*
*      Author: Gkougkousis Evangelos - Leonidas
*              egkougko@cern.ch
*                CERN-GENEVA
*/

#ifndef __LGADBase__
#define __LGADBase__

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm> 
#include <numeric>

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
#include "TCanvas.h"
#include "TFunction.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TPad.h"
#include "TText.h"
#include "TGaxis.h"

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <Windows4Root.h>
#include <windows.h>
#else
#include <dirent.h>
#include <bits/stdc++.h> 
#include <sys/stat.h>
#include <sys/types.h> 
#endif

enum AqInstrument { Sampic, LabTXT, TestBeamBin, TektronixScope, LeCroyWRBin, Unasigned, TestBeamBin1, TestBeamBin2,};
enum AqBoard { SingleCh, FourCh, IN2P3, KU, SiPM };
enum SecStage { MinCircuits, Particulars, CIVIDEC, none };
struct TrCrHist { AqBoard Board; float Capacitance; TH2D* TransHist; };

class LGADBase : public TSelector {

 public:

  bool m_convert;
  bool m_TrackComb;
  bool m_fei4Eff;
  bool m_WaveShape;
  bool m_TrnsCorr;
  TString m_ofname;
  TString m_ofdir;
  TString m_TransFileName;  // Transimpedence Simulation File
  TFile* m_TransFile;
  TFile* m_ofile;
  TTree* m_tree;
  std::vector<TrCrHist> m_TrsHists;
  unsigned int m_event;
  unsigned int m_nchan;
  unsigned int m_evnt1;
  unsigned int m_evnt2;
  AqInstrument m_instrument;

  // nTuple Variables
  std::vector<float> m_scale;
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
  void SetInstrument(AqInstrument instr = LabTXT);
  bool SetSRate(Long64_t rate, unsigned int ch);
  bool SetNPoints(unsigned int points, unsigned int ch);
  void SetInDataNames (TString DataDir = "", TString DataName = "", TString ext = "");
  void SetOutDataNames (TString DataDir = "", TString DataName = "");
  void SetTransFileName (TString filename = "");
  void SetDoTrnsCorr(bool TrnsCorr = false) { m_TrnsCorr = TrnsCorr; };
  void SetTrackComb (bool comb = false);
  void SetFEi4Eff(bool FEi4Eff = false);
  void SetTrackInDataNames (TString DataDir = "", TString DataName = "");
  void SetVerbose(int verbose = 1) { m_verbose = verbose; };
  void SetFitMethode(std::string method);
  void SetConvertSucess(bool convert) { m_convert = convert;};
  bool SetRootFile(const char* file) { m_ofile = TFile::Open(file); return m_ofile->IsOpen(); };
  bool SetRootTree(TFile* f, std::string name = "");
  void SetWaveShape(bool doShape = false) { m_WaveShape = doShape; };
  void SetTreeName(std::string treename) { m_treename = treename; };
  void SetDUTName(int ChId, std::string Name);
  void SetDUTNames(std::vector<std::string> DUTNames);
  void SetDUTBoard(int ChId, AqBoard Brd);
  void SetDUTBoards(std::vector<AqBoard> Brds);
  void SetDUTransImp(int ChId, float Transimp);
  void SetDUTransImps(std::vector<float> Transimps);
  void SetDUTSecStage(int ChId, SecStage Amp);
  void SetDUTSecStages(std::vector<SecStage> Amps);
  void SetDUTAmpGain(int ChId, int gain);
  void SetDUTAmpGains(std::vector<int> Gains);
  void SetDUTrigg(int ChId, float trigg);
  void SetDUTriggs(std::vector<float> Triggs);
  void SetDUTFract(int ChId, float frac);
  void SetDUTFracts(std::vector<float> Fracs);
  void SetDUTCap(int ChId, double cap);
  void SetDUTCaps(std::vector<double> Caps);
  void SetPlaneDMgt(int PlaneId1, int PlaneId2, double DMgt, int Qt);
  void SetPlaneDMgts(std::vector<double> DMgts, int Qt);
  void SetChMag(int ChId, double ChMags, int Qt);
  void SetChMags(std::vector<double> ChMags, int Qt);

  // Get Methodes
  AqInstrument GetInstrument();
  Long64_t GetSRate(unsigned int ch);
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
  bool GetDoTrnsCorr() { return m_TrnsCorr; };
  std::vector<TrCrHist> GetTransHistos() { return m_TrsHists; };
  TFile* GetRootFile() { return m_ofile; };
  TTree* GetRootTree() { return m_tree; };
  bool GetConvertSucess() { return m_convert; };
  bool GetWaveShape() { return m_WaveShape; };
  TString GetExtention() { return m_ext; };
  TString GetDataDir() { return m_datadir; };
  TString GetDataName() { return m_dataname; };
  TString GetOutFileaDir() { return m_ofdir; };
  TString GetOutFileName() { return m_ofname; };
  std::string GetTreeName() { return m_treename; };
  std::vector<std::pair<int, std::string>> GetDUTNames() { return m_DUTChsNames; };
  std::vector<std::pair<int, AqBoard>> GetDUTBoards() { return m_DUTChsBrd; };
  std::vector<std::pair<int, float>> GetDUTransImps() { return m_DUTChsTrns; };
  std::vector<std::pair<int, SecStage>> GetDUTSecStages() { return m_DUTChsAmp; };
  std::vector<std::pair<int, int>> GetDUTAmpGains() { return m_DUTChsAmpGn; };
  std::vector<std::pair<int, float>> GetDUTriggs() { return m_DUTChsTrigg; };
  std::vector<std::pair<int, float>> GetDUTFracts() { return m_DUTChsFrc; };
  std::vector<std::pair<int, double>> GetDUTCaps() { return m_DUTChsCaps; };
  std::vector<std::pair<std::pair<int, int>, double>> GetPlaneDTs() { return m_PlaneDTs; };
  std::vector<std::pair<std::pair<int, int>, double>> GetPlaneDCs() { return m_PlaneDCs; };
  std::vector<std::pair<int, double>> GetChVoltCuts() { return m_ChVoltCuts; };
  std::vector<std::pair<int, double>> GetChNoiseCuts() { return m_ChNoiseCuts; };
  std::vector<std::pair<int, double>> GetChJitterCuts() { return m_ChJitterCuts; };
  std::vector<std::pair<int, double>> GetChChargeCuts() { return m_ChChargeCuts; };

  // Class funcitons
  void SetStartStopEvnt(int Evnt1 = 0, int Evnt2 = 0);
  std::pair <unsigned int, unsigned int> GetStartStopEvnt();
  bool ConvertData();
  void Initialize();
  void SetVectorSize (unsigned int ch);

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
  template <typename V> double CalcMeadian(V *vec, int start = 0, int stop = 0);
  template <typename V> V OutlierReject(V *w, unsigned int order, float elem, int start = 0, int stop = 0);
  template <typename T, typename V> T CalResolution(V *w, unsigned int order, int start = 0, int stop = 0);
  template <typename V> double CalcFWHM(V* vec, double median = -99, int start = 0, int stop = 0);
  template <typename T, typename V> T MaxDensity(V* w, T res = -99, int start = 0, int stop = 0);

  bool ProgressBar(Long64_t evt, Long64_t total);
  int Addoriel(int val);
  bool PrintFitInfo(TH1D* histo, TCanvas** ca, std::string funcName = "none");

  // LGAD Fits
  int IterativeFit(std::vector<double> *w, std::pair <double, double> &gmean, std::pair <double, double> &gsigma, TH1D* &FitHist,
                   double &minchi2, std::string methode = "Gauss", std::pair <int, int> points = std::make_pair(-1, -1), bool discr = true);
  int LinearFit(std::vector<double>* vec, std::pair <double, double> &slope, std::pair <double, double> &intersept, std::vector<double>* vecErr = NULL);
  int RooConvFit(std::vector<double>* vec, std::pair <double, double> &magMPV, std::pair <double, double> &magSigma, std::string conv);
  double LinearInter(double x1, double y1, double x2, double y2, double y3);
  double FFT(std::vector<double> *w, Long64_t snrate, int start, int stop);

 private:

  // Member variables
  TString m_dataname;
  TString m_datadir;
  TString m_Trackdatadir;
  TString m_Trackdataname;
  TString m_ext;
  unsigned int m_bins; // number of bins for fitting histos
  int m_verbose;
  std::string m_treename;

  // Variables to be passed to Channels
  std::vector<std::pair<int, AqBoard>> m_DUTChsBrd;
  std::vector<std::pair<int, SecStage>> m_DUTChsAmp;
  std::vector<std::pair<int, float>> m_DUTChsFrc;
  std::vector<std::pair<int, float>> m_DUTChsTrigg;
  std::vector<std::pair<int, float>> m_DUTChsTrns;
  std::vector<std::pair<int, int>> m_DUTChsAmpGn;
  std::vector<std::pair<int, std::string>> m_DUTChsNames;
  std::vector<std::pair<int, double>> m_DUTChsCaps;
  std::vector<std::pair<std::pair<int, int>, double>> m_PlaneDTs;
  std::vector<std::pair<std::pair<int, int>, double>> m_PlaneDCs;
  std::vector<std::pair<int, double>> m_ChNoiseCuts;
  std::vector<std::pair<int, double>> m_ChVoltCuts;
  std::vector<std::pair<int, double>> m_ChChargeCuts;
  std::vector<std::pair<int, double>> m_ChJitterCuts;

  // Ntuple Variables
  std::vector<std::vector<double> > m_t;
  std::vector<std::vector<double> > m_w;
  TH1F* m_trigDt;
  TH1F* m_trigFr;

  bool WriteSampic(const char* dir, const char* name, const char* ext, int evt1, int evt2);
  bool WriteLabTXT(const char* dir, const char* name, const char* ext, int evt1, int evt2);
  bool WriteTectronixTXT(const char* dir, const char* name, const char* ext);
  bool WriteTestBeamBinary(const char* dir, const char* name, const char* ext, unsigned int evt1 = 0, unsigned int evt2 = 0);
  bool WriteLecroyBinary(const char* dir, const char* name, const char* ext, unsigned int evt1 = 0, unsigned int evt2 = 0);
  bool CombineTrack(const char* dir, const char* name);
  bool CreateOutputFile(const char* dir, const char* ofname, std::vector<unsigned int> nchan);
  bool SetScale(std::vector<unsigned int> channel, unsigned int nchan, std::vector<float>* scale);
  std::vector<double> ConrtVarBinX(std::vector<double> *wmod, double limUp, double limDown, int &nbins);
  bool CalcuRebin(bool discr, int n_elements, int nbins, double limDown, double limUp, double stdev, int (&Nofbins)[7]);
  void CalcTrigFr(std::vector<double>  EvTrigTime, TH1F* TrigPer, TH1F* TrigFrq, unsigned int entriesNo);
  bool OpenTransFile(TString filename = "");

  // Internal functions for iterative re-fitting
  TFitResultPtr Gauss(double rmin, double rmax, double strdv, double mean, TH1D* &Magnitude, std::string integ = "");
  TFitResultPtr GauXLandau(double rmin, double rmax, double median, double FWHM, TH1D* &Magnitude, std::string integ = "");
  Double_t LandXGauFun(Double_t *x, Double_t *par);
  int langaupro(TFitResultPtr fitResult, double &maxx, double &FWHM);

  // ClassDef(LGADBase, 1);
};

#endif