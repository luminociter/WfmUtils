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
#include <unordered_set>

#include "TSystem.h"
#include "TSystemFile.h"
#include "TSystemDirectory.h"
#include "TFitResultPtr.h"
#include "TObject.h"
#include "TClass.h"
#include "TROOT.h"
#include "TFile.h"
#include "TKey.h"
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
#include "TDirectory.h"

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

enum AqInstrument { Sampic, LabTXT, TestBeamBin, TektronixScope, LeCroyWRBin, Unasigned, TestBeamBin1, TestBeamBin2};
enum AqBoard { SingleCh, FourCh, IN2P3, KU, SiPM };
enum SecStage { MinCircuits, Particulars, CIVIDEC, none };
struct TrCrHist { AqBoard Board; float Capacitance; TH2D* TransHist; };

class LGADBase
{
 public:

  std::vector<TrCrHist> m_TrsHists;
  unsigned int m_event;
  unsigned int m_nchan;
  unsigned int m_evnt1;
  unsigned int m_evnt2;
  std::vector<double> m_oscdel;
  double m_trigclk;
  AqInstrument m_instrument;

  // nTuple Variables
  std::vector<std::vector<double>> m_t; 
  std::vector<std::vector<double>> m_w;
  std::vector<float> m_scale;
  std::vector<unsigned int> m_npoints;
  std::vector<unsigned int> m_scope;
  std::vector<Long64_t> m_srate;
  std::vector<double> m_triggTime;
  std::vector<double> m_physt;
  std::vector<double> m_ordrt;
  std::vector<unsigned int> m_channels;
  double m_trigtime;

  // Usefull common analysis varyables
  std::vector<unsigned int> m_nplanes;
  std::vector<unsigned int> m_nRef;
  std::vector<unsigned int> m_nDUT;
  std::vector<unsigned int> m_DUTs; // flattened 2 dimentional vector of DUT numbers
  std::vector<unsigned int> m_Refs; // flattened 2 dimentional vector of Reference numbers
  std::vector<unsigned int> m_planes; // flattened 2 dimentional vector of Plane numbers
  std::vector<unsigned int> m_ntp_planes;
  std::vector<unsigned int> m_ntp_Ref;
  std::vector<unsigned int> m_ntp_DUT;
  
  // ROI branches
  std::vector<TBranch*> b_m_ROI;
  std::vector<std::vector<double>*> m_ROI_Dbl;
  std::vector<std::vector<int>*> m_ROI_Int;
  // Reference branches
  std::vector<TBranch*> b_m_Ref;
  std::vector<std::vector<double>*> m_Ref_Dbl;
  std::vector<std::vector<int>*> m_Ref_Int;
  // DUT branches
  std::vector<TBranch*> b_m_DUT;
  std::vector<std::vector<double>*> m_DUT_Dbl;
  std::vector<std::vector<int>*> m_DUT_Int;
  // Telescope Plane Branches
  std::vector<TBranch*> b_m_TelPlane;
  std::vector<std::vector<double>*> m_TelPlane_Dbl;
  std::vector<std::vector<int>*> m_TelPlane_Int;

  // Setup Hits tree
  std::vector<TBranch*> b_m_Hits;
  std::vector<std::vector<int>*> m_HitsInt;
  // Setup Tracks tree
  std::vector<TBranch*> b_m_Tracks;
  std::vector<std::vector<int>*> m_TracksInt;
  std::vector<std::vector<double>*> m_TracksDbl;
  std::vector<Int_t> m_TrsInt;

  // Tracking leaf types for Cory
  // ROI Brnaches
  std::vector<double>* m_ClusterSizeX_ROI;
  std::vector<double>* m_ClusterSizeY_ROI;
  std::vector<double>* m_PixelX_ROI;
  std::vector<double>* m_PixelY_ROI;
  std::vector<int>* m_ClusterNumPixels_ROI;
  std::vector<double>* m_InterceptX_ROI;
  std::vector<double>* m_InterceptY_ROI;
  std::vector<double>* m_InterceptResX_ROI;
  std::vector<double>* m_InterceptResY_ROI;
  std::vector<int>* m_HasAsso;
  // Reference Branches
  std::vector<std::vector<double>*> m_InterceptX_Ref;
  std::vector<std::vector<double>*> m_InterceptY_Ref;
  std::vector<std::vector<double>*> m_InterceptResX_Ref;
  std::vector<std::vector<double>*> m_InterceptResY_Ref;
  // DUT Branches
  std::vector<std::vector<double>*> m_InterceptX_DUT;
  std::vector<std::vector<double>*> m_InterceptY_DUT;
  std::vector<std::vector<double>*> m_InterceptResX_DUT;
  std::vector<std::vector<double>*> m_InterceptResY_DUT;
  // Telescope Plane Branches
  std::vector<std::vector<double>*> m_InterceptX_TelPlane;
  std::vector<std::vector<double>*> m_InterceptY_TelPlane;
  std::vector<std::vector<double>*> m_InterceptResX_TelPlane;
  std::vector<std::vector<double>*> m_InterceptResY_TelPlane;
  std::vector<std::vector<double>*> m_ClusterSizeX_TelPlane;
  std::vector<std::vector<double>*> m_ClusterSizeY_TelPlane;
  std::vector<std::vector<int>*> m_ClusterNumPixels_TelPlane;

  // Tracking leaf types for EU Telescope
  std::vector<double>* m_xTr;
  std::vector<double>* m_yTr;
  std::vector<double>* m_omegaTr;
  std::vector<double>* m_phiTr;
  std::vector<double>* m_kinkxTr;
  std::vector<double>* m_kinkyTr;
  std::vector<double>* m_chi2Tr;
  std::vector<int>* m_plane_id;
  std::vector<int>* m_track_id;
  std::vector<int>* m_ndofTr;
  std::vector<int>* m_ID;
  std::vector<double>* m_xHit;
  std::vector<double>* m_yHit;
  std::vector<double>* m_zHit;
  Int_t       m_trigger_id;
  Int_t       m_timestampTr;

  std::vector<TTree*> m_trackTrees;

  // Constructor and destructor
  LGADBase();
  LGADBase(LGADBase& base);
  virtual ~LGADBase();
  
  // Set Methodes
  void SetChannels(std::vector<unsigned int> chan) { m_channels = chan; m_nchan = m_channels.size(); };
  void SetInstrument(AqInstrument instr = LabTXT);
  bool SetSRate(Long64_t rate, unsigned int ch);
  bool SetNPoints(unsigned int points, unsigned int ch);
  // Set in/out filenames and locations
  void SetInDataNames(TString DataDir = "", TString DataName = "", TString ext = "");
  void SetInFileNames(std::vector<TString> datanames) { m_datanames = datanames; };
  void SetInDataDir(TString dir) { m_datadir = dir; };
  void SetExtention(TString ext = "") { m_ext = ext; };
  int SetInRootFiles(std::vector<TString>& files, std::vector<unsigned int>& indx);
  void SetInRootFiles(std::vector<TFile*> files) { m_infiles = files; };
  int SetInRootFile(TFile* file);
  void SetOutDataNames(TString DataDir = "", TString DataName = "");
  void SetOutFileName(TString filename) { m_ofname = filename; };
  void SetOutFileDir(TString dir) { m_ofdir = dir; };
  bool SetOutRootFile(TFile* file);
  bool SetOutRootFile(const char* file);
  bool SetRootTrees(std::vector<TFile*> files, std::vector<unsigned int>& indx, std::string name = "");
  void SetRootTrees(std::vector<TTree*> trees) { m_trees = trees; };
  void SetTreeName(std::string treename) { m_treename = treename; };
  void SetAnaStage(unsigned int stage = 1) { m_stage = stage; };
  void SetFileLVL(unsigned int lvl = 0) { m_filelvl = lvl; };
  void SetTransFileName (TString filename = "");
  void SetTransFile(TFile* file) { m_TransFile = file; };
  void SetTrackInDataNames(TString DataDir = "", TString DataName = "");
  // Set all the boolean variables
  void SetDoTrnsCorr(bool TrnsCorr = false) { m_TrnsCorr = TrnsCorr; };
  void SetTrackComb(bool comb = false) { m_TrackComb = comb; };
  void SetHasFFT(bool hasfft = false) { m_hasfft = hasfft; };
  void SetHasTrck(bool hastrck = false) { m_hastrck = hastrck; };
  void SetDoFFT(bool dofft = false) { m_dofft = dofft; };
  void SetFEi4Eff(bool FEi4Eff = false) { m_fei4Eff = FEi4Eff; };
  void SetWaveShape(bool doShape = false) { m_WaveShape = doShape; };
  void SetHasWaveShape(bool hasShape = false) { m_hasWaveShape = hasShape; };
  void SetConvertSucess(bool convert = false) { m_convert = convert;};
  void SetVerbose(int verbose = 1) { m_verbose = verbose; };
  void SetTrigClk (double clk = 25e-9) { m_trigclk = clk; };
  void SetScopeDelay(std::vector<double> delay = {50e-9, 162e-9}) { m_oscdel = delay; };
  void SetFitMethode(std::string method);
  void SetTrackPackage(std::string package) { m_trkpck = package; };
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
  void SetTestEvn(int Evn = -1) { m_testEvn = Evn; };
  void SetEvntNo(unsigned int evnt) { m_event = evnt; };
  void SetTransHistos(std::vector<TrCrHist> histos) { m_TrsHists = histos; };
  void SetExcludeTrackFiles(std::vector<unsigned int> trackExclude) { m_trackExclude = trackExclude; };

  // Get Methodes
  AqInstrument GetInstrument();
  Long64_t GetSRate(unsigned int ch);
  unsigned int GetNPoints(unsigned int ch);
  unsigned int GetChNo() { return m_nchan; };
  std::vector<unsigned int> GetChannels() { return m_channels; };
  // Input file names and directory
  std::vector<TString> GetInFileNames() { return m_datanames; };
  TString GetInDataDir() { return m_datadir; };
  TString GetExtention() { return m_ext; };
  std::vector<TFile*> GetInRootFiles() { return m_infiles; };
  TFile* GetInRootFile(unsigned int indx) { return m_infiles.at(indx); };
  // Output file name and directory
  TString GetOutFileName() { return m_ofname; };
  TString GetOutFileDir() { return m_ofdir; };
  TFile* GetOutRootFile() { return m_ofile; };
  // Data Tree
  std::vector<TTree*> GetRootTrees() { return m_trees; };
  TTree* GetRootTree(unsigned int indx) { return m_trees.at(indx); };
  std::string GetTreeName() { return m_treename; };
  // Transimpedence file name
  TString GetTransFileName() { return m_TransFileName; };
  TFile* GetTransFile() { return m_TransFile; };
  // Boolean Switches
  bool GetTrackComb() { return m_TrackComb; };
  bool GetDoFFT() { return m_dofft; };
  bool GetHasFFT() { return m_hasfft; };
  bool GetHasTrck() { return m_hastrck; };
  bool GetFEi4Eff() { return m_fei4Eff; };
  bool GetDoTrnsCorr() { return m_TrnsCorr; };
  bool GetConvertSucess() { return m_convert; };
  bool GetWaveShape() { return m_WaveShape; };
  bool GetHasWaveShape() { return m_hasWaveShape; };
  // Debug levele
  int GetVerbosity() { return m_verbose; };
  // Analysis options
  unsigned int GetAnaStage() { return m_stage; };
  unsigned int GetFileLVL() { return m_filelvl; };
  std::string GetFitMethode() { return m_fitopt; };
  std::string GetTrackPackage() { return m_trkpck; };
  unsigned int GetEvntNo() { return m_event; };
  std::vector<TrCrHist> GetTransHistos() { return m_TrsHists; };
  std::vector<double> GetScopeDelay() { return m_oscdel; };
  double GetTrigClk() { return m_trigclk; };
  // Chanel by channel properties and interchanel cuts
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
  bool GetTestEvn(UInt_t Evn);
  UInt_t GetTestEvtNm() { return m_testEvn; };
  std::vector<unsigned int> GetExcludeTrackFiles() { return m_trackExclude; };

  // Class funcitons
  void SetStartStopEvnt(int Evnt1 = 0, int Evnt2 = 0);
  std::pair <unsigned int, unsigned int> GetStartStopEvnt();
  bool ConvertData();
  void Initialize();
  void SetVectorSize (unsigned int ch);

  // General helper functions
  std::string reduce(const std::string str, const std::string fill, const std::string whitespace);
  std::string trim(const std::string str, const std::string whitespace);
  unsigned int CountFiles(const char* dir, const char* ext = nullptr);
  int DirExists(const char* path);
  int CreateDir(const char* path);
  int RecursMkDir(const char* dirname);
  std::vector<std::string> ListFileNames(const char* path, const char* ext = nullptr);
  template <typename T> T Derivate(T *w, int start = 1);
  template <typename A> bool IsVecEqual(std::vector<A>& first, std::vector<A>& second);
  template <typename U> double Mean(U *w, int start = -1, int stop = -1);
  template <typename V> double Stdev(V *w, int start = -1, int stop = -1, double mean = 0.0);
  template <typename V, typename T> double BayesianErr(V *w, T value = 1);
  template <typename V> double CalcMeadian(V *vec, int start = 0, int stop = 0);
  template <typename V> V OutlierReject(V *w, unsigned int order, float elem, int start = 0, int stop = 0);
  template <typename T, typename V> T CalResolution(V *w, unsigned int order, int start = 0, int stop = 0);
  template <typename V> double CalcFWHM(V* vec, double median = -99, int start = 0, int stop = 0);
  template <typename T, typename V> T MaxDensity(V* w, T res = -99, int start = 0, int stop = 0);

  void PrintBranches(TTree* tree, std::string name);
  void PrintObjects(TFile* rootfile, std::string ignore);
  bool ProgressBar(Long64_t evt, Long64_t total);
  int Addoriel(int val);
  bool PrintFitInfo(TH1D* histo, TCanvas** ca, std::string funcName = "none");
  bool OpenTrackFile(TString trackdirdir, std::vector<TString> trackdataname, std::vector<unsigned int>& filesIndx);
  bool AddTrckBranches(TTree* c_tree, std::string trckpkg);
  void InitEUTelPointers();
  void InitTBAPointers();
  void InitCoryPointers(unsigned int nRef, unsigned int nDUT, unsigned int nplanes);
  bool CombineTrackCory(int EvntNmbr, std::vector<unsigned int>& evntlast, unsigned int file = 0);
  bool CombineTrackEUTel(int EvntNmbr, std::vector<unsigned int>& evntlast, unsigned int file = 0);
  bool CombineTrackTBA(int EvntNmbr, std::vector<unsigned int>& evntlast, unsigned int file = 0);
  void CloseTrackFiles();
  void InportProperties(LGADBase& base);
  void ExportProperties(LGADBase& base);
  int ArrangeCahnnels(std::vector<unsigned int>& channels, std::vector<unsigned int>& ChScope);
  bool FixBaseQuant(std::vector<std::vector<unsigned int>>& filevect, std::vector<unsigned int>& basevect, unsigned int& basecounter);

  // LGAD Fits
  int IterativeFit(std::vector<double> *w, std::pair <double, double> &gmean, std::pair <double, double> &gsigma, TH1D* &FitHist,
                   double &minchi2, std::string methode = "Gauss", std::pair <int, int> points = std::make_pair(-1, -1), bool discr = true);
  int LinearFit(std::vector<double>* vec, std::pair <double, double> &slope, std::pair <double, double> &intersept, std::vector<double>* vecErr = NULL);
  int RooConvFit(std::vector<double>* vec, std::pair <double, double> &magMPV, std::pair <double, double> &magSigma, std::string conv);
  double LinearInter(double x1, double y1, double x2, double y2, double y3);
  double FFT(std::vector<double> *w, Long64_t snrate, int start, int stop);

 private:

  // Member variables
  std::vector<TString> m_datanames;
  std::vector<TFile*> m_trackfiles;
  TString m_datadir;
  TString m_ext;
  std::vector<TFile*> m_infiles;
  TString m_ofname;
  TString m_ofdir;
  TFile* m_ofile;
  std::vector<TTree*> m_trees;
  TString m_Trackdatadir;
  std::vector<TString> m_Trackdataname;
  TString m_TransFileName;  // Transimpedence Simulation File
  TFile* m_TransFile;
  unsigned int m_bins; // number of bins for fitting histos
  int m_verbose;
  unsigned int m_stage;
  unsigned int m_filelvl;
  std::string m_treename;
  bool m_convert;
  bool m_fei4Eff;
  bool m_WaveShape;
  bool m_TrnsCorr;
  bool m_TrackComb;
  bool m_dofft;
  bool m_hasfft;
  bool m_hastrck;
  bool m_hasWaveShape;
  std::string m_fitopt;
  std::string m_trkpck;
  int m_testEvn;

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
  std::vector<unsigned int> m_trackExclude;

  // Ntuple Variables
  TH1F* m_trigDt;
  TH1F* m_trigFr;

  bool WriteSampic(const char* dir, std::vector<TString> names, const char* ext, unsigned int evt1 = 0, unsigned int evt2 = 0);
  bool WriteLabTXT(const char* dir, std::vector<TString> names, const char* ext, unsigned int evt1 = 0, unsigned int evt2 = 0);
  bool WriteTectronixTXT(const char* dir, std::vector<TString> names, const char* ext, unsigned int evt1 = 0, unsigned int evt2 = 0);
  bool WriteTestBeamBinary(const char* dir, std::vector<TString> names, const char* ext, unsigned int evt1 = 0, unsigned int evt2 = 0);
  bool WriteLecroyBinary(const char* dir, std::vector<TString> names, const char* ext, unsigned int evt1 = 0, unsigned int evt2 = 0);
  bool CreateOutputFile(const char* dir, const char* ofname, std::vector<unsigned int> nchan);
  bool SetScale(std::vector<unsigned int> channel, unsigned int nchan, std::vector<float>* scale);
  bool CheckStartStopEvnt(unsigned int* evt1, unsigned int* evt2, unsigned int min_evnt, unsigned int max_evnt);
  int RunCount(std::vector<TString>& names, std::vector<unsigned int>& nfiles, std::vector<TString>& Runs, std::vector<TString>& ext, 
               std::vector<unsigned int>& runIndx, const char* inext, std::vector<std::vector<unsigned int>>* nchannels = NULL);
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