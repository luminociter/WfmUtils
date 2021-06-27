/*
* LGADURun.h
*
*
*      Author: Gkougkousis Evangelos - Leonidas
*              egkougko@cern.ch
*               IFAE-BARCELONA
*/

#ifndef __LGADRun__
#define __LGADRun__

#include "DUTChannel.h"
#include "LGADSel.h"
#include <TParameter.h>

class LGADRun : public LGADSel {

public:

    /// constructors
    LGADRun (TTree * /*tree*/ = 0);
    LGADRun(LGADBase* rBase, TTree* /*tree*/ = 0);

    /// Destructor
    virtual ~LGADRun();

    /// Overall initialization
    void Begin(TTree *tree);

    /// Initialization for thread
    void SlaveBegin(TTree *tree);
    
    /// The Notify function called everytime a file is accessed
    Bool_t Notify();
     
    /// Process an event
    Bool_t Process(Long64_t entry);

    /// End of thread
    void SlaveTerminate();

    /// End of processing
    void Terminate();

    bool SetPlaneDQuants(std::vector<double> dQunat, std::vector<std::pair<std::pair<int, int>, double>> DTs);

private:

    void VarInit();
    void FirstInit();

    // Need to asign Ntuple variables here
    std::vector<int> m_EvPol;
    std::vector<double> m_EvCharge;
    std::vector<int> m_EvIsSignal;
    std::vector<int> m_EvInWindow;
    std::vector<int> m_EvVoltStr;
    std::vector<int> m_EvMaxindx;
    std::vector<double> m_EvMaxTime;
    std::vector<double> m_EvMaxVolt;
    std::vector<int> m_EvMinIndx;
    std::vector<double> m_EvMinTime;
    std::vector<double> m_EvMinVolt;
    std::vector<int> m_EvStrIndx;
    std::vector<int> m_EvEndIndx;
    std::vector<double> m_EvNoise;
    std::vector<double> m_EvPedest;
    std::vector<double> m_EvNoiseErr;
    std::vector<double> m_EvPedestErr;
    std::vector<double> m_EvRiseT;
    std::vector<double> m_EvJitRiseSNR;
    std::vector<double> m_EvTriggTime;
    std::vector<double> m_EvDVDTMax;
    std::vector<double> m_EvTriggToT;
    std::vector<double> m_EvSignalFFT;
    std::vector<double> m_EvNoiseFFT;
    std::vector<int> m_EvComplete;
    std::vector< std::vector<double>> m_EvVAdjCFD;
    std::vector< std::vector<double>> m_EvCFDTime;
    std::vector< std::vector<double>> m_EvDVDTCFD;
    std::vector< std::vector<double>> m_EvCFDToT;
    std::vector< std::vector<double>> m_EvJitNdVdT;
    std::vector<std::vector<std::vector<double>>> m_EvTmDiff;
    std::vector<std::vector<TH1D* >> h_TmDiffCFD;
    std::vector<std::vector<TCanvas* >> h_TmDiffCFDCanv;
    std::vector<TH2D* > h_TmSigmas;

    // InterChannel Cut Criteria
    std::vector<double> m_dTMax;
    std::vector<double> m_dCMax;

    // Branches pointers for Ntuple
    std::vector<TBranch*> b_EvPol;
    std::vector<TBranch*> b_EvCharge;
    std::vector<TBranch*> b_EvIsSignal;
    std::vector<TBranch*> b_EvInWindow;
    std::vector<TBranch*> b_EvVoltStr;
    std::vector<TBranch*> b_EvMaxindx;
    std::vector<TBranch*> b_EvMaxTime;
    std::vector<TBranch*> b_EvMaxVolt;
    std::vector<TBranch*> b_EvMinIndx;
    std::vector<TBranch*> b_EvMinTime;
    std::vector<TBranch*> b_EvMinVolt;
    std::vector<TBranch*> b_EvStrIndx;
    std::vector<TBranch*> b_EvEndIndx;
    std::vector<TBranch*> b_EvNoise;
    std::vector<TBranch*> b_EvPedest;
    std::vector<TBranch*> b_EvNoiseErr;
    std::vector<TBranch*> b_EvPedestErr;
    std::vector<TBranch*> b_EvRiseT;
    std::vector<TBranch*> b_EvJitNdVdT;
    std::vector<TBranch*> b_EvJitRiseSNR;
    std::vector<TBranch*> b_EvTriggTime;
    std::vector<TBranch*> b_EvDVDTMax;
    std::vector<TBranch*> b_EvDVDTCFD;
    std::vector<TBranch*> b_EvCFDToT;
    std::vector<TBranch*> b_EvTriggToT;
    std::vector<TBranch*> b_EvSignalFFT;
    std::vector<TBranch*> b_EvNoiseFFT;
    std::vector<TBranch*> b_EvComplete;
    std::vector<TBranch*> b_EvVAdjCFD;
    std::vector<TBranch*> b_EvCFDTime;

    std::vector<TCanvas* > h_MaxIndxFt;
    std::vector<TCanvas* > h_MinIndxFt;
    std::vector<TCanvas* > h_MaxVoltFt;
    std::vector<TCanvas* > h_MinVoltFt;
    std::vector<TCanvas* > h_StrIndxFt;
    std::vector<TCanvas* > h_EndIndxFt;
    std::vector<TCanvas* > h_NoiseFt;
    std::vector<TCanvas* > h_NoiseErrFt;
    std::vector<TCanvas* > h_PedestFt;
    std::vector<TCanvas* > h_NPedestErrFt;
    std::vector<TCanvas* > h_MaxTimeFt;
    std::vector<TCanvas* > h_MinTimeFt;
    std::vector<TCanvas* > h_ChargeFt;
    std::vector<TCanvas* > h_RiseTimeFt;
    std::vector<TCanvas* > h_TriggTimeFt;
    std::vector<TCanvas* > h_DVDTMaxFt;
    std::vector<TCanvas* > h_TriggToTFt;
    std::vector<TCanvas* > h_SignalFFTFt;
    std::vector<TCanvas* > h_NoiseFFTFt;
    std::vector< std::vector<TCanvas* > > h_CFDTimeFt;
    std::vector< std::vector<TCanvas* > > h_DVDTCFDFt;
    std::vector< std::vector<TCanvas* > > h_CFDToTFt;

    std::vector<DUTChannel*> m_RunDUTCh;
    // std::vector<WaveForm*> m_RunWaveform;
    LGADBase* m_RunBase;
    unsigned int m_init;

    // ClassDef(LGADRun, 1);
};

#endif
