/*
* Waveform.h
*
*
*      Author: Gkougkousis Evangelos - Leonidas
*              egkougko@cern.ch
*               CERN-GENEVA
*/

#ifndef __WaveForm__
#define __WaveForm__

#include "LGADBase.h"

class WaveForm : public LGADBase 
{
public:

    LGADBase* m_WvBase;
    enum polarity {pos, neg, undef};
    enum basepos {poside, negside, zerocross, undet};

    WaveForm();
    WaveForm(std::vector<double>* voltage, std::vector<double>* time);
    WaveForm(std::vector<double>* voltage, Long64_t snrate);
    WaveForm(std::vector<double>* voltage, std::vector<double>* time, Long64_t snrate);
    WaveForm(LGADBase *Base);
    WaveForm(std::vector<double>* voltage, std::vector<double>* time, LGADBase *Base);
    WaveForm(std::vector<double>* voltage, Long64_t snrate, LGADBase *Base);
    WaveForm(std::vector<double>* voltage, std::vector<double>* time, Long64_t snrate, LGADBase *Base);
    virtual ~WaveForm();
    void InitializeWaveForm(LGADBase* tBase, int level = 0);
    bool Calculate();

    void SetVoltage(std::vector<double>* volt);
    void SetTime(std::vector<double>* time);
    void SetPolarity(unsigned int pol);
    void SetSnRate(Long64_t snrate);
    void SetCFDfraction(double fraction);
    void SetTrigg(double tri);
    void SetBasePos(unsigned int pos);
    void SetTransimp(float transimp);
    void SetAmpGain(int gain);
    void SetWaveId(int id);
    void SetTrnsHist(TH2D* hist);

    int GetPolarity();
    Long64_t GetSnRate();
    double GetTrigg();
    float GetCFDfraction();
    std::vector<double>* GetVoltage();
    std::vector<double>* GetTime();
    std::vector<double> GetAdjVoltage();
    std::vector<double> GetTimeAdjVolt(float fraction = -99);
    std::vector<int> GetHalfPoints();
    int GetBasePos();
    bool GetIsSignal();
    bool GetIsInWidow();
    bool GetSaturation();
    float GetTransimp();
    float GetAmpGain();
    int GetMaxIndx();
    double GetMaxTime();
    double GetMax();
    double GetMaxErr();
    int GetMinIndx();
    double GetMinTime();
    double GetMin();
    int GetStrIndx();
    int GetEndIndx();
    double GetNoise();
    double GetNoiseErr();
    double GetPedestal();
    double GetPedestalErr();
    double GetCharge();
    double GetRiseTime(float bottom = -99, float top = -99);
    double GetCFDTime(float fraction = -99);
    double GetTriggTime(float trigg = -99);
    double GetdVdTMax();
    double GetdVdTCFD(float fraction = -99, int ndif = 0);
    double GetCFDToT(float fraction = -99);
    double GetTriggToT(float trigg = -99);
    double GetFrequency(int start = -1, int stop = -1);
    int GetFitQual() { return m_NoiseFtQl; };
    double GetNoiseFitChi2() { return m_fitchi2; };
    TH1D* GetNoiseFit() { return (TH1D*)m_noiseFit->Clone(); };
    int GetWaveId() { return m_waveId; };
    double GetJitterNdVdT(float fraction = -99, int ndif = 0);
    double GetJitterRiseSNR();
    double GetSignalFFT();
    double GetNoiseFFT();
    TH2D* GetTrnsHist();
    LGADBase* GetBase() { return m_WvBase; };
    bool dump();

private:

    polarity m_pol;
    basepos m_basepos;
    Long64_t m_SnRate;

    int m_waveId;
    int m_maxIndx;
    double m_maxTime;
    double m_maxVolt;
    double m_maxVoltErr;
    int m_minIndx;
    double m_minTime;
    double m_minVolt;
    int m_StrIndx;
    int m_EndIndx;

    //Non NTuple variables usefull though the class
    double m_maxavg;
    double m_minavg;
    double m_wfmavg;
    double m_wfmstdv;

    double m_noise;
    double m_noiseErr;
    double m_pedestal;
    double m_pedestalErr;
    double m_fitchi2;
    TH1D* m_noiseFit;
    TH2D* m_trnsHist;
    int m_NoiseFtQl;

    double m_charge;
    double m_RiseTime;
    double m_jitter1;
    double m_jitter2;
    double m_CFDTime;
    double m_dVdTMax;
    double m_dVdTCFD;
    double m_TriggTime;

    double m_CFDToT;
    double m_TriggToT;
    double m_signalFFT;
    double m_noiseFFT;

    float m_transimp;
    int m_ampgain;
    float m_fraction;
    float m_trigg;

    bool m_IsSignal;
    bool m_IsInWindow;
    bool m_VoltStr;

    std::vector<double> *m_voltage;
    std::vector<double> *m_time;
    std::vector<double> m_intime;
    std::vector<double> m_voltageAdj;
    std::vector<double> m_vAdjCFD;
    std::vector<int> m_halfPtsIndx;

    basepos FindPosition(std::vector<double> *w);
    polarity FindPolarity(std::vector<double> *w);

    int VoltMaxIndx(std::vector<double> *w, bool adj, bool poldef = true);
    double VoltMax(std::vector<double> *w, bool adj, bool poldef = true);
    double VoltMaxErr(std::vector<double> *w, bool adj);
    double TimeMax(std::vector<double> *w, Long64_t snrate, bool adj, bool poldef = true);
    int VoltMinIndx(std::vector<double> *w, bool adj, bool poldef = true);
    double VoltMin(std::vector<double> *w, bool adj, bool poldef = true);
    double TimeMin(std::vector<double> *w, Long64_t snrate, bool adj, bool poldef = true);

    bool IsSignal(std::vector<int> *w);
    bool InWindow(std::vector<double> *w);
    bool VoltSatur(std::vector<double> *w, bool poldef = true);
    std::pair <unsigned int, unsigned int> FindNoisePoints(std::vector<double> *w);
    std::vector<double> PulseAdj(std::vector<double> *w, double baseline, int factor = 1);
    std::vector<double> PulseTimeVoltAdj(std::vector<double> *t, std::vector<double> *w, float fraction);
    double CollectedCharge(std::vector<double> *w, Long64_t snrate, float transimp, float ampgain, int start = -1, int stop = -1);

    double RiseTimeLinear(std::vector<double> *w, Long64_t snrate, float top = 0.9, float bottom = 0.1);
    double CFDTimeLinear(std::vector<double> *w, std::vector<double> *t, float fraction);
    double CFDToTLinear(std::vector<double> *w, std::vector<double> *t, float fraction);
    double TriggToTLinear(std::vector<double> *w, std::vector<double> *t, double trig);
    double dVdTMaxLinear(std::vector<double> *w, Long64_t snrate);
    double dVdTCFDLinear(std::vector<double> *w, Long64_t snrate, float fraction, int ndif = 0);
    float GetTrsFromHisto(double signalFFT, TH2D* transhist);
        
    int StartIndx(std::vector<double> *w, bool adj = false, bool poldef = true);
    int EndIndx(std::vector<double> *w, bool adj = false, bool poldef = true);
    std::vector<int> FindHalfPoints(std::vector<double> *w, bool adj = false);
    double FirstTimeForVoltage(std::vector<double> *w, std::vector<double> *t, double volt);
    double SecondTimeForVoltage(std::vector<double> *w, std::vector<double> *t, double volt);
    std::vector<double>* FillTime(std::vector<double>* voltage, Long64_t snrate);

    // ClassDef(WaveForm, 1);
};

#endif