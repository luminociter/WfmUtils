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

#include "LGADUtils/DUTChannel.h"
#include "LGADUtils/LGADSel.h"
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

private:

    void VarInit();
    void FirstInit();

    // Nedd to asign Ntuple variables here
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
    std::vector<double> m_EvJitNdVdT;
    std::vector<double> m_EvJitRiseSNR;
    std::vector<double> m_EvTriggTime;
    std::vector<double> m_EvDVDTMax;
    std::vector<double> m_EvDVDTCFD;
    std::vector<double> m_EvCFDToT;
    std::vector<double> m_EvTriggToT;
    std::vector<double> m_EvSignalFFT;
    std::vector<double> m_EvNoiseFFT;
    std::vector<int> m_EvComplete;
    std::vector< std::vector<double>> m_EvVAdjCFD;
    std::vector< std::vector<double>> m_EvCFDTime;

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

    std::vector<TH1D* > h_MaxIndxFt;
    std::vector<TH1D* > h_MinIndxFt;
    std::vector<TH1D* > h_MaxVoltFt;
    std::vector<TH1D* > h_MinVoltFt;
    std::vector<TH1D* > h_StrIndxFt;
    std::vector<TH1D* > h_EndIndxFt;
    std::vector<TH1D* > h_NoiseFt;
    std::vector<TH1D* > h_NoiseErrFt;
    std::vector<TH1D* > h_PedestFt;
    std::vector<TH1D* > h_NPedestErrFt;
    std::vector<TH1D* > h_MaxTimeFt;
    std::vector<TH1D* > h_MinTimeFt;
    std::vector<TH1D* > h_ChargeFt;
    std::vector<TH1D* > h_RiseTimeFt;
    std::vector<TH1D* > h_TriggTimeFt;
    std::vector<TH1D* > h_DVDTMaxFt;
    std::vector<TH1D* > h_DVDTCFDFt;
    std::vector<TH1D* > h_CFDToTFt;
    std::vector<TH1D* > h_TriggToTFt;
    std::vector<TH1D* > h_SignalFFTFt;
    std::vector<TH1D* > h_NoiseFFTFt;
    std::vector< std::vector<TH1D* > > h_CFDTimeFt;

    std::vector<DUTChannel*> m_RunDUTCh;
    // std::vector<WaveForm*> m_RunWaveform;
    // ClassDef(LGADRun, 0);
    LGADBase* m_RunBase;

    unsigned int m_init;
};

#endif
