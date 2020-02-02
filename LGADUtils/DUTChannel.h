/*
* DUTChannel.h
*
*
*      Author: Gkougkousis Evangelos - Leonidas
*              egkougko@cern.ch
*               IFAE-BARCELONA
*/

#ifndef __DUTChannel__
#define __DUTChannel__

#include "LGADUtils/WaveForm.h"
#include "LGADUtils/LGADBase.h"

class DUTChannel : public LGADBase {

public:

  enum AqBorad { SingleCh, FourCh, IN2P3, KU, SiPM};
  enum SecStage { MinCircuits, Particulars, CIVIDEC, none};

  // Constructors
  DUTChannel();
  DUTChannel(int channelID, TString name);
  DUTChannel(int channelID);
  DUTChannel(TString name);
  DUTChannel(int channelID, TString name, unsigned int board, unsigned int amp, float fraction, float trigger);
  DUTChannel(LGADBase *tBase);
  DUTChannel(int channelID, TString name, LGADBase *tBase);
  DUTChannel(int channelID, LGADBase *tBase);
  DUTChannel(TString name, LGADBase *tBase);
  DUTChannel(int channelID, TString name, unsigned int board, unsigned int amp, float fraction, float trigger, LGADBase *tBase);
  virtual ~DUTChannel();
  
  // Set Functions
  void SetBoard(unsigned int board);
  void SetSecStage(unsigned int ampl);
  bool SetChName(TString name);
  bool SetChID(int channelID);
  bool SetChRate(Long64_t rate);
  bool SetChTransImp(float transimp);
  bool SetChAmpGain(int gain);
  bool SetChFract(float fraction);
  bool SetChTrigg(double trigger);
  int ChDump();

  // Get Methods
  AqBorad GetBoard() { return m_board; };
  SecStage GetAmp() { return m_amp; };
  TString GetChName() { return m_channelName; };
  int GetChID() { return m_channelID; };
  Long64_t GetChRate() { return m_ChRate; };
  float GetChTransImp() { return m_ChTransimp; };
  int GetChAmpGain() { return m_ChAmpgain; };
  float GetChFract() { return m_ChFraction; };
  double GetChFTrigg() { return m_ChTrigg; };
  std::pair<int, int> GetChMaxIndx(int methode = 1);
  std::pair<double, double> GetChMaxVolt(int methode = 1);
  std::pair<int, int> GetChMinIndx(int methode = 1);
  std::pair<double, double> GetChMinVolt(int methode = 1);
  std::pair<double, double> GetChIsSignal();
  std::pair<double, double> GetChPolarity();
  std::pair<double, double> GetChInWindow();
  std::pair<double, double> GetChVoltStr();
  std::pair<int, int> GetChStrIndx(int methode = 1);
  std::pair<int, int> GetChEndIndx(int methode = 1);
  std::pair<double, double> GetChNoise(int methode = 1);
  std::pair<double, double> GetChPedestal();
  std::pair<double, double> GetChNoiseSlope();
  std::pair<double, double> GetChPedestalSlope();
  std::pair<double, double> GetChNoiseErr(int methode = 1);
  std::pair<double, double> GetChPedestalErr(int methode = 1);
  std::pair<double, double> GetChMaxTime(int methode = 1);
  std::pair<double, double> GetChMinTime(int methode = 1);
  std::pair<double, double> GetChCharge(int methode = 1);
  std::pair<double, double> GetChRiseTime(int methode = 1);
  std::pair<double, double> GetChJitNdVdT(int methode = 1);
  std::pair<double, double> GetChJitRiseSNR(int methode = 1);
  std::vector<std::pair<double, double>> GetChCFDTime(int methode = 1);
  std::pair<double, double> GetChTriggTime(int methode = 1);
  std::pair<double, double> GetChDVDTMax(int methode = 1);
  std::pair<double, double> GetChDVDTCFD(int methode = 1);
  std::pair<double, double> GetChCFDToT(int methode = 1);
  std::pair<double, double> GetChTriggToT(int methode = 1);
  std::pair<double, double> GetChSignalFFT(int methode = 1);
  std::pair<double, double> GetChNoiseFFT(int methode = 1);
  std::pair<double, double> GetChComplete();
  std::pair<double, double> GetChSoNR(int methode = 1);
  std::vector<std::pair<double, double>> GetMeanChPulse() { return m_MeanChPulse; };
  double GetFitChi2(std::string fit, int i = 3);
  TH1D* GetFit(std::string fit, int i = 3);

  // Set of Inilined functions to get single event variables
  // They return current value if one assigned
  // mirrored form waveform class
  int GetEvntPolarity() { return GetEvntQuantInt(&m_ChPol); };
  int GetEvntIsSignal() { return GetEvntQuantBool(&m_ChIsSignal); };
  int GetEvntInWindow() { return GetEvntQuantBool(&m_ChIsInWindow); };
  int GetEvntVoltStr() { return GetEvntQuantBool(&m_ChVoltStr); };
  int GetEvntMaxindx() { return GetEvntQuantInt(&m_ChMaxIndx); };
  double GetEvntMaxTime() { return GetEvntQuantDouble(&m_ChMaxTime); };
  double GetEvntMaxVolt() { return GetEvntQuantDouble(&m_ChMaxVolt); };
  int GetEvntMinIndx() { return GetEvntQuantInt(&m_ChMinIndx); };
  double GetEvntMinTime() { return GetEvntQuantDouble(&m_ChMinTime); };
  double GetEvntMinVolt() { return GetEvntQuantDouble(&m_ChMinVolt); };
  int GetEvntStrIndx() { return GetEvntQuantInt(&m_ChStrIndx); };
  int GetEvntEndIndx() { return GetEvntQuantInt(&m_ChEndIndx); };
  double GetEvntNoise() { return GetEvntQuantDouble(&m_ChNoise); };
  double GetEvntPedest() { return GetEvntQuantDouble(&m_ChPedestal); };
  double GetEvntNoiseErr() { return GetEvntQuantDouble(&m_ChNoiseErr); };
  double GetEvntPedestErr() { return GetEvntQuantDouble(&m_ChPedestalErr); };
  double GetEvntCharge() { return GetEvntQuantDouble(&m_ChCharge); };
  double GetEvntRiseT() { return GetEvntQuantDouble(&m_ChRiseTime); };
  double GetEvntJitNdVdT() { return GetEvntQuantDouble(&m_ChJitNdVdT); };
  double GetEvntJitRiseSNR() { return GetEvntQuantDouble(&m_ChJitRiseSNR); };
  std::vector<int> GetEvntHalfPts() { return m_ChHalfPts; };
  std::vector<double> GetEvntVAdjCFD() { return m_ChVAdjCFD; };
  double GetEvntTriggTime() { return GetEvntQuantDouble(&m_ChTriggTime); };
  double GetEvntDVDTMax() { return GetEvntQuantDouble(&m_ChDVDTMax); };
  double GetEvntDVDTCFD() { return GetEvntQuantDouble(&m_ChDVDTCFD); };
  double GetEvntCFDToT() { return GetEvntQuantDouble(&m_ChCFDToT); };
  double GetEvntTriggToT() { return GetEvntQuantDouble(&m_ChTriggToT); };
  double GetEvntSignalFFT() { return GetEvntQuantDouble(&m_ChSignalFFT); };
  double GetEvntNoiseFFT() { return GetEvntQuantDouble(&m_ChNoiseFFT); };
  int GetEvntComplete() { return GetEvntQuantInt(&m_ChComplete); };
  double GetEvntCFDTime(int bin) { return GetEvntQuantDouble(&(m_ChCFDTime[bin])); };
  WaveForm* GetWaveform() { return m_waveform; };

  // Publicly accessible methods
  void ChannelInit(/*unsigned int fevent, unsigned int points*/);
  bool AppendEvent(std::vector<double> *volt, std::vector<double>* time, bool waveshape = false);
  int updateChProperties(bool waveshape = false, TTree* wavetree = NULL);

private:

  WaveForm* m_waveform;
  int m_channelID;
  TString m_channelName;
  Long64_t m_ChRate;
  float m_ChTransimp;
  int m_ChAmpgain;
  float m_ChFraction;
  double m_ChTrigg;
  int m_ChNoPulses;
  AqBorad m_board;
  SecStage m_amp;

  // Outword variables, Mean pulse propertirs
  int m_ChMeanPolarity;
  int m_ChMeanMaxIndx;
  double m_ChMeanMax;
  int m_ChMeanMinIndx;
  double m_ChMeanMin;
  int m_ChMeanStrIndx;
  int m_ChMeanEndIndx;
  double m_ChMeanMaxTime;
  double m_ChMeanMinTime;
  double m_ChMeanNoise;
  double m_ChMeanNoiseErr;
  double m_chMeanMaxErr;
  double m_ChMeanCharge;
  double m_ChMeanRiseTime;
  double m_ChMeanJitterNdVdT;
  double m_ChMeanJitterRiseSNR;
  double m_ChMeanTriggTime;
  double m_ChMeanDVDTMax;
  double m_ChMeanDVDTCFD;
  double m_ChMeanCFDToT;
  double m_ChMeanTriggToT;
  double m_ChMeanSignalFFT;
  double m_ChMeanNoiseFFT;
  std::pair<double, double> m_ChMeanSoN;
  std::vector<double> m_ChMeanCFDTime;

  // Outword variables, propertirs from pulse fit
  std::pair<int, int> m_ChFitMaxIndx;
  std::pair<double, double> m_ChFitMaxVolt;
  std::pair<int, int> m_ChFitMinIndx;
  std::pair<double, double> m_ChFitMinVolt;
  std::pair<double, double> m_ChFitIsSignal;
  std::pair<double, double> m_ChFitIsInWindow;
  std::pair<double, double> m_ChFitVoltStr;
  std::pair<int, int> m_ChFitStrIndx;
  std::pair<int, int> m_ChFitEndIndx;
  std::pair<double, double> m_ChFitNoise;
  std::pair<double, double> m_ChFitNoiseSlope;
  std::pair<double, double> m_ChFitPedestalSlope;
  std::pair<double, double> m_ChFitPedestal;
  std::pair<double, double> m_ChFitNoiseErr;
  std::pair<double, double> m_ChFitPedestalErr;
  std::pair<double, double> m_ChFitMaxTime;
  std::pair<double, double> m_ChFitMinTime;
  std::pair<double, double> m_ChFitCharge;
  std::pair<double, double> m_ChFitRiseTime;
  std::pair<double, double> m_ChFitJitNdVdT;
  std::pair<double, double> m_ChFitJitRiseSNR;
  std::vector<std::pair<double, double>> m_FitChCFDTime;
  std::pair<double, double> m_ChFitTriggTime;
  std::pair<double, double> m_ChFitDVDTMax;
  std::pair<double, double> m_ChFitDVDTCFD;
  std::pair<double, double> m_ChFitCFDToT;
  std::pair<double, double> m_ChFitTriggToT;
  std::pair<double, double> m_ChFitSignalFFT;
  std::pair<double, double> m_ChFitNoiseFFT;
  std::pair<double, double> m_ChFitComplete;
  std::pair<double, double> m_ChFitSoN;

  TH1D* m_ChFitMaxIndxFt;
  TH1D* m_ChFitMaxVoltFt;
  TH1D* m_ChFitMinIndxFt;
  TH1D* m_ChFitMinVoltFt;
  TH1D* m_ChFitStrIndxFt;
  TH1D* m_ChFitEndIndxFt;
  TH1D* m_ChNoiseErrFt;
  TH1D* m_ChNoiseFt;
  TH1D* m_ChPedestalFt;
  TH1D* m_ChPedestalErFt;
  TH1D* m_ChMaxTimeFt;
  TH1D* m_ChMinTimeFt;
  TH1D* m_ChChargeFt;
  TH1D* m_ChFitRiseTimeFt;
  std::vector<TH1D* > m_ChCFDTimeFt;
  TH1D* m_ChTriggTimeFt;
  TH1D* m_ChDVDTMaxFt;
  TH1D* m_ChDVDTCFDFt;
  TH1D* m_ChCFDToTFt;
  TH1D* m_ChTriggToTFt;
  TH1D* m_ChSignalFFTFt;
  TH1D* m_ChNoiseFFTFt;

  double m_ChFitMaxIndxChi2;
  double m_ChFitMaxVoltChi2;
  double m_ChFitMinIndxChi2;
  double m_ChFitMinVoltChi2;
  double m_ChFitStrIndxChi2;
  double m_ChFitEndIndxChi2;
  double m_ChNoiseErrChi2;
  double m_ChNoiseChi2;
  double m_ChPedestalChi2;
  double m_ChPedestalErChi2;
  double m_ChMaxTimeChi2;
  double m_ChMinTimeChi2;
  double m_ChChargeChi2;
  double m_ChFitRiseTimeChi2;
  std::vector<double> m_ChCFDTimeChi2;
  double m_ChTriggTimeChi2;
  double m_ChDVDTMaxChi2;
  double m_ChDVDTCFDChi2;
  double m_ChCFDToTChi2;
  double m_ChTriggToTChi2;
  double m_ChSignalFFTChi2;
  double m_ChNoiseFFTChi2;

  // Internal vectors
  std::vector<double> m_ChVAdjCFD;
  std::vector<bool> m_ChIsSignal;
  std::vector<bool> m_ChIsInWindow;
  std::vector<bool> m_ChVoltStr;
  std::vector<int> m_ChComplete;
  std::vector<int> m_ChPol;
  std::vector<int> m_ChMaxIndx;
  std::vector<double> m_ChMaxTime;
  std::vector<double> m_ChMaxVolt;
  std::vector<int> m_ChMinIndx;
  std::vector<double> m_ChMinTime;
  std::vector<double> m_ChMinVolt;
  std::vector<int> m_ChStrIndx;
  std::vector<int> m_ChEndIndx;
  std::vector<double> m_ChNoise;
  std::vector<double> m_ChPedestal;
  std::vector<double> m_ChNoiseErr;
  std::vector<double> m_ChPedestalErr;
  std::vector<double> m_ChCharge;
  std::vector<double> m_ChRiseTime;
  std::vector<double> m_ChJitNdVdT;
  std::vector<double> m_ChJitRiseSNR;
  std::vector<double> m_ChCFDTime[19];
  std::vector<double> m_ChTriggTime;
  std::vector<double> m_ChDVDTMax;
  std::vector<double> m_ChDVDTCFD;
  std::vector<double> m_ChCFDToT;
  std::vector<double> m_ChTriggToT;
  std::vector<double> m_ChSignalFFT;
  std::vector<double> m_ChNoiseFFT;
  std::vector<int> m_ChHalfPts;

  std::vector<std::pair<double, double> > m_MeanChPulse;
  std::vector<TH1D*> m_MeanChPlFt;
  std::vector<double> m_MeanChPlFtChi2;

  bool FillDefault(int k);
  bool AppendWaveform(WaveForm* k);
  int CalculateShape(TTree* wavetree);
  std::pair<double, double> CalculateMeanG(vector<double>* vec, std::vector<int>* signal, TH1D* &hitHist,
                                           double &fitChi2, int level, std::string fitopt = "root");
  std::pair<double, double> CalculateMeanGLandau(std::vector<double>* vec, std::vector<int>* signal, TH1D* &hitHist,
                                                double &fitChi2, int level, std::string fitopt = "root");
  std::pair<double, double> CalculateMeanLinear(vector<double>* vec, std::vector<bool>* signal, vector<double>* vecErr = NULL);
  std::pair<double, double> CalculateMeanGLinear(vector<double>* vec, std::vector<int>* signal, TH1D* &hitHist,
                                                 double &fitChi2, int level, std::string fitopt = "root");
  std::pair<double, double> calculateSoN(double noise, double noiseErr, double VMax, double VMaxErr);
  std::pair<double, double> GetChValue(int med, double ChMean, double ChMeanErr, std::pair<double, double> ChFit);
  int GetEvntQuantInt(std::vector<int>* ChQuant);
  int GetEvntQuantBool(std::vector<bool>* ChQuant);
  double GetEvntQuantDouble(std::vector<double>* ChQuant);

  // ClassDef(DUTChannel, 0);
};

#endif
