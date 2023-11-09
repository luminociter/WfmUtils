/*
* WfmURun.h
*
*
*      Author: Gkougkousis Evangelos - Leonidas
*              egkougko@cern.ch
*               CERN-GENEVA
*/

#ifndef __WfmRun__
#define __WfmRun__

#include "DUTChannel.h"
#include "TParameter.h"

class WfmRun : public WfmBase {

public:

    WfmRun();                // Default constructor
    WfmRun(WfmBase* rBase); // Constructor from WfmBase object
    ~WfmRun();               // Destructor
    bool InitAnalysis(std::vector<unsigned int>& fileindx);

private:

    bool InitNtuple(const std::vector<TTree*>& inputrees, const std::vector<TFile*>& inputfiles);
    bool SetBranchAccess(unsigned int file);
    bool GetEvntBranches(Long64_t entry, std::vector<bool> FFTShape, unsigned int filelvl, unsigned int trck, unsigned int fileindx);
    bool GetChannelBranches(Long64_t channel, unsigned int fileindx, std::vector<bool> FFTShape, unsigned int filelvl);
    bool InitOutputBranches();
    bool SetOutputBranches();
    // Function to add the 2nd and 3rd stage analysis histograms to the output file
    bool AddOutputHistos();
    // Function to only write to the output tree the branches of added variable in case the framework is adding first stage analysis resutlts on existing file
    bool WriteBranches(unsigned int fileindx);

    void VarInit(WfmBase* lbase, unsigned int file);
    bool SetPlaneDQuants(std::vector<double> dQunat, std::vector<std::pair<std::pair<int, int>, double>> DQs);
    void FirstInit(std::pair<unsigned int, unsigned int> filelvl, unsigned int anastage);
    void FillFistLvl(WfmBase* lbase, WaveForm* wvmf, unsigned int ch);
    void WvfFillDefault(WfmBase* lbase, unsigned int completeness, unsigned int ch);
    void FillSecondLvl(WfmBase* lbase, DUTChannel* DUTCh, unsigned int ch);
    bool FillSecFirst(WfmBase* lbase, DUTChannel* DUTCh, unsigned int ch);

    std::vector<std::vector<unsigned int>> m_filechannel;
    std::vector<std::vector<unsigned int>> m_fileplanes;
    std::vector<std::vector<unsigned int>> m_fileRefs;
    std::vector<std::vector<unsigned int>> m_fileDUTs;
    unsigned int m_nRecDUTs;
    unsigned int m_nRecRefs;
    unsigned int m_nRecPlanes;

    long long int m_totalevnts;

    std::vector<unsigned int> m_filelevels;
    std::vector<std::vector<bool>> m_FFTShDef;
    std::vector<std::pair<bool, bool>> m_waveshape;
    std::vector<unsigned int> m_filetrack;
    std::vector<TTree*> m_ChTrees;
    std::vector<TTree*> m_TimeResTrees;
    std::vector<TTree*> m_ShapeTrees;

    // Inter-Channel Cut Criteria
    std::vector<double> m_dTMax;
    std::vector<double> m_dCMax;

    TTree* m_OutTree;
    TTree* m_ChFitParam;
    TTree* m_ChShapeParam;
    TTree* m_TimeResParam;

    int m_cp;

    //==================Branches for the I/O filess==================
 
    // List of 0th stage branches
    std::vector<TBranch*> b_t;   //!
    std::vector<TBranch*> b_w;   //!
    std::vector<TBranch*> b_vScale;   //!
    std::vector<TBranch*> b_nPoints;   //!
    std::vector<TBranch*> b_SnRate;   //!
    std::vector<TBranch*> b_ordrt;   //!   
    std::vector<TBranch*> b_physt;   //!
    std::vector<TBranch*> b_EvnTimeCh;   //!
    std::vector<TBranch*> b_Scope;   //!
    TBranch* b_EvnTime;   //!
    TBranch* b_EvnNo;   //!

    // Output Ntuple Branches in case of copying converted file to new one
    std::vector<TBranch*> b_t_out;   //!
    std::vector<TBranch*> b_w_out;   //!
    std::vector<TBranch*> b_vScale_out;   //!
    std::vector<TBranch*> b_nPoints_out;   //!
    std::vector<TBranch*> b_SnRate_out;   //!
    std::vector<TBranch*> b_ordrt_out;   //!   
    std::vector<TBranch*> b_physt_out;   //!
    std::vector<TBranch*> b_EvnTimeCh_out;   //!
    std::vector<TBranch*> b_Scope_out;   //!
    TBranch* b_EvnTime_out;   //!
    TBranch* b_EvnNo_out;   //!

    // List of 1st stage analysis branches
    std::vector<TBranch*> b_EvPol;   //!
    std::vector<TBranch*> b_EvCharge;   //!
    std::vector<TBranch*> b_EvIsSignal;   //!
    std::vector<TBranch*> b_EvInWindow;   //!
    std::vector<TBranch*> b_EvVoltStr;   //!
    std::vector<TBranch*> b_EvMaxindx;   //!
    std::vector<TBranch*> b_EvMaxTime;   //!
    std::vector<TBranch*> b_EvMaxVolt;   //!
    std::vector<TBranch*> b_EvMinIndx;   //!
    std::vector<TBranch*> b_EvMinTime;   //!
    std::vector<TBranch*> b_EvMinVolt;   //!
    std::vector<TBranch*> b_EvStrIndx;   //!
    std::vector<TBranch*> b_EvEndIndx;   //!
    std::vector<TBranch*> b_EvNoise;   //!
    std::vector<TBranch*> b_EvPedest;   //!
    std::vector<TBranch*> b_EvNoiseErr;   //!
    std::vector<TBranch*> b_EvPedestErr;   //!
    std::vector<TBranch*> b_EvRiseT;   //!
    std::vector<TBranch*> b_EvJitNdVdT;   //!
    std::vector<TBranch*> b_EvJitRiseSNR;   //!
    std::vector<TBranch*> b_EvTriggTime;   //!
    std::vector<TBranch*> b_EvDVDTMax;   //!
    std::vector<TBranch*> b_EvDVDTCFD;   //!
    std::vector<TBranch*> b_EvCFDToT;   //!
    std::vector<TBranch*> b_EvTriggToT;   //!
    std::vector<TBranch*> b_EvSignalFFT;   //!
    std::vector<TBranch*> b_EvNoiseFFT;   //!
    std::vector<TBranch*> b_EvComplete;   //!
    std::vector<TBranch*> b_EvVAdjCFD;   //!
    std::vector<TBranch*> b_EvCFDTime;   //!
    TBranch* b_FitDiffChi2;   //!
    TBranch* b_FitDiffMag;   //!
    TBranch* b_FitDiffMagErr;   //!
    TBranch* b_FitDiffSigma;   //!
    TBranch* b_FitDiffSigmaErr;   //!
    TBranch* b_FitDiffTime;   //!

    // Output Ntuple branches in case of copying 1st stage analysis to new file
    std::vector<TBranch*> b_EvPol_out;   //!
    std::vector<TBranch*> b_EvCharge_out;   //!
    std::vector<TBranch*> b_EvIsSignal_out;   //!
    std::vector<TBranch*> b_EvInWindow_out;   //!
    std::vector<TBranch*> b_EvVoltStr_out;   //!
    std::vector<TBranch*> b_EvMaxindx_out;   //!
    std::vector<TBranch*> b_EvMaxTime_out;   //!
    std::vector<TBranch*> b_EvMaxVolt_out;   //!
    std::vector<TBranch*> b_EvMinIndx_out;   //!
    std::vector<TBranch*> b_EvMinTime_out;   //!
    std::vector<TBranch*> b_EvMinVolt_out;   //!
    std::vector<TBranch*> b_EvStrIndx_out;   //!
    std::vector<TBranch*> b_EvEndIndx_out;   //!
    std::vector<TBranch*> b_EvNoise_out;   //!
    std::vector<TBranch*> b_EvPedest_out;   //!
    std::vector<TBranch*> b_EvNoiseErr_out;   //!
    std::vector<TBranch*> b_EvPedestErr_out;   //!
    std::vector<TBranch*> b_EvRiseT_out;   //!
    std::vector<TBranch*> b_EvJitNdVdT_out;   //!
    std::vector<TBranch*> b_EvJitRiseSNR_out;   //!
    std::vector<TBranch*> b_EvTriggTime_out;   //!
    std::vector<TBranch*> b_EvDVDTMax_out;   //!
    std::vector<TBranch*> b_EvDVDTCFD_out;   //!
    std::vector<TBranch*> b_EvCFDToT_out;   //!
    std::vector<TBranch*> b_EvTriggToT_out;   //!
    std::vector<TBranch*> b_EvSignalFFT_out;   //!
    std::vector<TBranch*> b_EvNoiseFFT_out;   //!
    std::vector<TBranch*> b_EvComplete_out;   //!
    std::vector<TBranch*> b_EvVAdjCFD_out;   //!
    std::vector<TBranch*> b_EvCFDTime_out;   //!
    TBranch* b_FitDiffChi2_out;   //!
    TBranch* b_FitDiffMag_out;   //!
    TBranch* b_FitDiffMagErr_out;   //!
    TBranch* b_FitDiffSigma_out;   //!
    TBranch* b_FitDiffSigmaErr_out;   //!
    TBranch* b_FitDiffTime_out;   //!

    // List of 2nd stage analysis branches
    TBranch* b_chName;   //!
    TBranch* b_SubChName;   //!
    TBranch* b_ChFitMaxVolt;   //!
    TBranch* b_ChFitMaxVoltErr;   //!
    TBranch* b_ChFitMinVolt;   //!
    TBranch* b_ChFitMinVoltErr;   //!
    TBranch* b_ChFitMxIndx;   //!
    TBranch* b_ChFitMxIndxErr;   //!
    TBranch* b_ChFitMinIndx;   //!
    TBranch* b_ChFitMinIndxErr;   //!
    TBranch* b_ChFitStrIndx;   //!
    TBranch* b_ChFitStrIndxErr;   //!
    TBranch* b_ChFitEndIndx;   //!
    TBranch* b_ChFitEndIndxErr;   //!
    TBranch* b_ChFitIsSig;   //!
    TBranch* b_ChFitIsSigErr;   //!
    TBranch* b_ChFitInWin;   //!
    TBranch* b_ChFitInWinErr;   //!
    TBranch* b_ChFitVStr;   //!
    TBranch* b_ChFitVStrErr;   //!
    TBranch* b_ChFitNoise;   //!
    TBranch* b_ChFitNoiseUnc;   //!
    TBranch* b_ChFitPedestal;   //!
    TBranch* b_ChFitPedestalUnc;   //!
    TBranch* b_ChFitNoiseSlp;   //!
    TBranch* b_ChFitNoiseSlpErr;   //!
    TBranch* b_ChFitPedestalSlp;   //!
    TBranch* b_ChFitPedestalSlpErr;   //!
    TBranch* b_ChFitNoiseErr;   //!
    TBranch* b_ChFitNoiseErrUnc;   //!
    TBranch* b_ChFitPedestalErr;   //!
    TBranch* b_ChFitPedestalErrUnc;   //!
    TBranch* b_ChFitMaxTime;   //!
    TBranch* b_ChFitMaxTimeErr;   //!
    TBranch* b_ChFitMinTime;   //!
    TBranch* b_ChFitMinTimeErr;   //!
    TBranch* b_ChFitCharge;   //!
    TBranch* b_ChFitChargeErr;   //!
    TBranch* b_ChFitRiseT;   //!
    TBranch* b_ChFitRiseTErr;   //!
    TBranch* b_ChFitJitRiseSNR;   //!
    TBranch* b_ChFitJitRiseSNRErr;   //!
    TBranch* b_ChFitCFDTime;   //!
    TBranch* b_CFDTimeFtChi2;   //!
    TBranch* b_ChFitDVDTCFD;   //!
    TBranch* b_DVDTCFDFtChi2;   //!
    TBranch* b_ChFitCFDToT;   //!
    TBranch* b_CFDToTFtChi2;   //!
    TBranch* b_ChFitJitNdVdT;   //!
    TBranch* b_ChFitTriggTime;   //!
    TBranch* b_ChFitTriggTimeErr;   //!
    TBranch* b_ChFitDVDTMax;   //!
    TBranch* b_ChFitDVDTMaxErr;   //!
    TBranch* b_ChFitTriggToT;   //!
    TBranch* b_ChFitTriggToTErr;   //!
    TBranch* b_ChFitComplete;   //!
    TBranch* b_ChFitCompleteErr;   //!
    TBranch* b_ChFitSoNR;   //!
    TBranch* b_ChFitSoNRErr;   //!
    TBranch* b_ChFitSignalFFT;   //!
    TBranch* b_ChFitSignalFFTErr;   //!
    TBranch* b_ChFitNoiseFFT;   //!
    TBranch* b_ChFitNoiseFFTErr;   //!

    TBranch* b_SignalFFTFtChi2;   //!
    TBranch* b_NoiseFFTFtChi2;   //!
    TBranch* b_MaxIndxFtChi2;   //!
    TBranch* b_MinIndxFtChi2;   //!
    TBranch* b_MaxVoltFtChi2;   //!
    TBranch* b_MinVoltFtChi2;   //!
    TBranch* b_StrIndxFtChi2;   //!
    TBranch* b_EndIndxFtChi2;   //!
    TBranch* b_NoiseFtChi2;   //!
    TBranch* b_NoiseErrFtChi2;   //!
    TBranch* b_PedestFtChi2;   //!
    TBranch* b_PedestErrFtChi2;   //!
    TBranch* b_MaxTimeFtChi2;   //!
    TBranch* b_MinTimeFtChi2;   //!
    TBranch* b_ChargeFtChi2;   //!
    TBranch* b_RiseTimeFtChi2;   //!
    TBranch* b_TriggTimeFtChi2;   //!
    TBranch* b_DVDTMaxFtChi2;   //!
    TBranch* b_TriggToTFtChi2;   //!

    // List of WaveShape 2nd Stage analysis branches
    TBranch* b_ChShapeName;   //!
    TBranch* b_SubChShName;   //!
    TBranch* b_ChShapePolarity;   //!
    TBranch* b_ChShapeMaxIndx;   //!
    TBranch* b_ChShapeMax;   //!
    TBranch* b_ChShapeMinIndx;   //!
    TBranch* b_ChShapeMin;   //!
    TBranch* b_ChShapeStrIndx;   //!
    TBranch* b_ChShapeEndIndx;   //!
    TBranch* b_ChShapeMaxTime;   //!
    TBranch* b_ChShapeMinTime;   //!
    TBranch* b_ChShapeNoise;   //!
    TBranch* b_ChShapeNoiseErr;   //!
    TBranch* b_ChShapeCharge;   //!;
    TBranch* b_ChShapeRiseTime;   //!
    TBranch* b_ChShapeJitRiseSNR;   //!
    TBranch* b_ChShapeTriggTime;   //!
    TBranch* b_ChShapeDVDTMax;   //!
    TBranch* b_ChShapeTriggToT;   //!
    TBranch* b_ChShapeSoN;   //!
    TBranch* b_ChShapeSoNUnc;   //!
    TBranch* b_ChShapeCFDTime;   //!
    TBranch* b_ChShapeDVDTCFD;   //!
    TBranch* b_ChShapeCFDToT;   //!
    TBranch* b_ChShapeJitNdVdT;   //!
    TBranch* b_ChMeanVolt;   //!;
    TBranch* b_ChShapeSignalFFT;   //!
    TBranch* b_ChShapeNoiseFFT;   //!

    // Output Ntuple branches in case of copying 2nd stage analysis to new file
    TBranch* b_chName_out;   //!
    TBranch* b_SubChName_out;   //!
    TBranch* b_ChFitMaxVolt_out;   //!
    TBranch* b_ChFitMaxVoltErr_out;   //!
    TBranch* b_ChFitMinVolt_out;   //!
    TBranch* b_ChFitMinVoltErr_out;   //!
    TBranch* b_ChFitMxIndx_out;   //!
    TBranch* b_ChFitMxIndxErr_out;   //!
    TBranch* b_ChFitMinIndx_out;   //!
    TBranch* b_ChFitMinIndxErr_out;   //!
    TBranch* b_ChFitStrIndx_out;   //!
    TBranch* b_ChFitStrIndxErr_out;   //!
    TBranch* b_ChFitEndIndx_out;   //!
    TBranch* b_ChFitEndIndxErr_out;   //!
    TBranch* b_ChFitIsSig_out;   //!
    TBranch* b_ChFitIsSigErr_out;   //!
    TBranch* b_ChFitInWin_out;   //!
    TBranch* b_ChFitInWinErr_out;   //!
    TBranch* b_ChFitVStr_out;   //!
    TBranch* b_ChFitVStrErr_out;   //!
    TBranch* b_ChFitNoise_out;   //!
    TBranch* b_ChFitNoiseUnc_out;   //!
    TBranch* b_ChFitPedestal_out;   //!
    TBranch* b_ChFitPedestalUnc_out;   //!
    TBranch* b_ChFitNoiseSlp_out;   //!
    TBranch* b_ChFitNoiseSlpErr_out;   //!
    TBranch* b_ChFitPedestalSlp_out;   //!
    TBranch* b_ChFitPedestalSlpErr_out;   //!
    TBranch* b_ChFitNoiseErr_out;   //!
    TBranch* b_ChFitNoiseErrUnc_out;   //!
    TBranch* b_ChFitPedestalErr_out;   //!
    TBranch* b_ChFitPedestalErrUnc_out;   //!
    TBranch* b_ChFitMaxTime_out;   //!
    TBranch* b_ChFitMaxTimeErr_out;   //!
    TBranch* b_ChFitMinTime_out;   //!
    TBranch* b_ChFitMinTimeErr_out;   //!
    TBranch* b_ChFitCharge_out;   //!
    TBranch* b_ChFitChargeErr_out;   //!
    TBranch* b_ChFitRiseT_out;   //!
    TBranch* b_ChFitRiseTErr_out;   //!
    TBranch* b_ChFitJitRiseSNR_out;   //!
    TBranch* b_ChFitJitRiseSNRErr_out;   //!
    TBranch* b_ChFitCFDTime_out;   //!
    TBranch* b_CFDTimeFtChi2_out;   //!
    TBranch* b_ChFitDVDTCFD_out;   //!
    TBranch* b_DVDTCFDFtChi2_out;   //!
    TBranch* b_ChFitCFDToT_out;   //!
    TBranch* b_CFDToTFtChi2_out;   //!
    TBranch* b_ChFitJitNdVdT_out;   //!
    TBranch* b_ChFitTriggTime_out;   //!
    TBranch* b_ChFitTriggTimeErr_out;   //!
    TBranch* b_ChFitDVDTMax_out;   //!
    TBranch* b_ChFitDVDTMaxErr_out;   //!
    TBranch* b_ChFitTriggToT_out;   //!
    TBranch* b_ChFitTriggToTErr_out;   //!
    TBranch* b_ChFitComplete_out;   //!
    TBranch* b_ChFitCompleteErr_out;   //!
    TBranch* b_ChFitSoNR_out;   //!
    TBranch* b_ChFitSoNRErr_out;   //!
    TBranch* b_ChFitSignalFFT_out;   //!
    TBranch* b_ChFitSignalFFTErr_out;   //!
    TBranch* b_ChFitNoiseFFT_out;   //!
    TBranch* b_ChFitNoiseFFTErr_out;   //!
    TBranch* b_SignalFFTFtChi2_out;   //!
    TBranch* b_NoiseFFTFtChi2_out;   //!
    TBranch* b_MaxIndxFtChi2_out;   //!
    TBranch* b_MinIndxFtChi2_out;   //!
    TBranch* b_MaxVoltFtChi2_out;   //!
    TBranch* b_MinVoltFtChi2_out;   //!
    TBranch* b_StrIndxFtChi2_out;   //!
    TBranch* b_EndIndxFtChi2_out;   //!
    TBranch* b_NoiseFtChi2_out;   //!
    TBranch* b_NoiseErrFtChi2_out;   //!
    TBranch* b_PedestFtChi2_out;   //!
    TBranch* b_PedestErrFtChi2_out;   //!
    TBranch* b_MaxTimeFtChi2_out;   //!
    TBranch* b_MinTimeFtChi2_out;   //!
    TBranch* b_ChargeFtChi2_out;   //!
    TBranch* b_RiseTimeFtChi2_out;   //!
    TBranch* b_TriggTimeFtChi2_out;   //!
    TBranch* b_DVDTMaxFtChi2_out;   //!
    TBranch* b_TriggToTFtChi2_out;   //!

    // Output Ntuple branches in case of copying WaveShape 2nd Stage analysis to a new file
    TBranch* b_ChShapeName_out;   //!
    TBranch* b_SubChShName_out;   //!
    TBranch* b_ChShapePolarity_out;   //!
    TBranch* b_ChShapeMaxIndx_out;   //!
    TBranch* b_ChShapeMax_out;   //!
    TBranch* b_ChShapeMinIndx_out;   //!
    TBranch* b_ChShapeMin_out;   //!
    TBranch* b_ChShapeStrIndx_out;   //!
    TBranch* b_ChShapeEndIndx_out;   //!
    TBranch* b_ChShapeMaxTime_out;   //!
    TBranch* b_ChShapeMinTime_out;   //!
    TBranch* b_ChShapeNoise_out;   //!
    TBranch* b_ChShapeNoiseErr_out;   //!
    TBranch* b_ChShapeCharge_out;   //!;
    TBranch* b_ChShapeRiseTime_out;   //!
    TBranch* b_ChShapeJitRiseSNR_out;   //!
    TBranch* b_ChShapeTriggTime_out;   //!
    TBranch* b_ChShapeDVDTMax_out;   //!
    TBranch* b_ChShapeTriggToT_out;   //!
    TBranch* b_ChShapeSoN_out;   //!
    TBranch* b_ChShapeSoNUnc_out;   //!
    TBranch* b_ChShapeCFDTime_out;   //!
    TBranch* b_ChShapeDVDTCFD_out;   //!
    TBranch* b_ChShapeCFDToT_out;   //!
    TBranch* b_ChShapeJitNdVdT_out;   //!
    TBranch* b_ChMeanVolt_out;   //!;
    TBranch* b_ChShapeSignalFFT_out;   //!
    TBranch* b_ChShapeNoiseFFT_out;   //!

    // List of 3rd stage analysis branches
    TBranch* b_TmDiffChi2;   //!
    TBranch* b_TmDiffFitMean;   //!
    TBranch* b_TmDiffFitErr;   //!
    TBranch* b_CFD_ChA;   //!
    TBranch* b_ChA;   //!
    TBranch* b_ChB;   //!

    // Output Ntuple branches in case of copying 3rd stage analysis to new file
    TBranch* b_TmDiffChi2_out;   //!
    TBranch* b_TmDiffFitMean_out;   //!
    TBranch* b_TmDiffFitErr_out;   //!
    TBranch* b_CFD_ChA_out;   //!
    TBranch* b_ChA_out;   //!
    TBranch* b_ChB_out;   //!

    // Tracking Branches Pointers for EU Telescope
    TBranch* b_plane_id;   //!
    TBranch* b_track_id;   //!
    TBranch* b_trigger_id;   //!
    TBranch* b_timestampTr;   //!
    TBranch* b_xTr;   //!
    TBranch* b_yTr;   //!
    TBranch* b_omegaTr;   //!
    TBranch* b_phiTr;   //!
    TBranch* b_kinkxTr;   //!
    TBranch* b_kinkyTr;   //!
    TBranch* b_chi2Tr;   //!
    TBranch* b_ndofTr;   //!
    TBranch* b_xHit;   //!
    TBranch* b_yHit;   //!
    TBranch* b_zHit;   //!
    TBranch* b_ID;   //!

    // Copies of the tracking braches in case of copying file to new one
    TBranch* b_plane_id_out;   //!
    TBranch* b_track_id_out;   //!
    TBranch* b_trigger_id_out;   //!
    TBranch* b_timestampTr_out;   //!
    TBranch* b_xTr_out;   //!
    TBranch* b_yTr_out;   //!
    TBranch* b_omegaTr_out;   //!
    TBranch* b_phiTr_out;   //!
    TBranch* b_kinkxTr_out;   //!
    TBranch* b_kinkyTr_out;   //!
    TBranch* b_chi2Tr_out;   //!
    TBranch* b_ndofTr_out;   //!
    TBranch* b_xHit_out;   //!
    TBranch* b_yHit_out;   //!
    TBranch* b_zHit_out;   //!
    TBranch* b_ID_out;   //!

    // Tracking Branches Pointers for Cory

    // ROI Branches
    TBranch* b_ClusterSizeX_ROI;   //!
    TBranch* b_ClusterSizeY_ROI;   //!
    TBranch* b_PixelX_ROI;   //!
    TBranch* b_PixelY_ROI;   //!
    TBranch* b_ClusterNumPixels_ROI;   //!
    TBranch* b_InterceptX_ROI;   //!
    TBranch* b_InterceptY_ROI;   //!
    TBranch* b_InterceptResX_ROI;   //!
    TBranch* b_InterceptResY_ROI;   //!
    TBranch* b_HasAsso;   //!
    // Reference Branches
    std::vector<TBranch*> b_InterceptX_Ref;   //!
    std::vector<TBranch*> b_InterceptY_Ref;   //!
    std::vector<TBranch*> b_InterceptResX_Ref;   //!
    std::vector<TBranch*> b_InterceptResY_Ref;   //!
    // DUT Branches
    std::vector<TBranch*> b_InterceptX_DUT;   //!
    std::vector<TBranch*> b_InterceptY_DUT;   //!
    std::vector<TBranch*> b_InterceptResX_DUT;   //!
    std::vector<TBranch*> b_InterceptResY_DUT;   //!
    // Telescope Plane Branches
    std::vector<TBranch*> b_InterceptX_TelPlane;   //!
    std::vector<TBranch*> b_InterceptY_TelPlane;   //!
    std::vector<TBranch*> b_InterceptResX_TelPlane;   //!
    std::vector<TBranch*> b_InterceptResY_TelPlane;   //!
    std::vector<TBranch*> b_ClusterSizeX_TelPlane;   //!
    std::vector<TBranch*> b_ClusterSizeY_TelPlane;   //!
    std::vector<TBranch*> b_ClusterNumPixels_TelPlane;   //!

    // Copies of the tracking braches in case of copying file to new one
    TBranch* b_ClusterSizeX_ROI_out;   //!
    TBranch* b_ClusterSizeY_ROI_out;   //!
    TBranch* b_PixelX_ROI_out;   //!
    TBranch* b_PixelY_ROI_out;   //!
    TBranch* b_ClusterNumPixels_ROI_out;   //!
    TBranch* b_InterceptX_ROI_out;   //!
    TBranch* b_InterceptY_ROI_out;   //!
    TBranch* b_InterceptResX_ROI_out;   //!
    TBranch* b_InterceptResY_ROI_out;   //!
    TBranch* b_HasAsso_out;   //!
    // Reference Branches
    std::vector<TBranch*> b_InterceptX_Ref_out;   //!
    std::vector<TBranch*> b_InterceptY_Ref_out;   //!
    std::vector<TBranch*> b_InterceptResX_Ref_out;   //!
    std::vector<TBranch*> b_InterceptResY_Ref_out;   //!
    // DUT Branches
    std::vector<TBranch*> b_InterceptX_DUT_out;   //!
    std::vector<TBranch*> b_InterceptY_DUT_out;   //!
    std::vector<TBranch*> b_InterceptResX_DUT_out;   //!
    std::vector<TBranch*> b_InterceptResY_DUT_out;   //!
    // Telescope Plane Branches
    std::vector<TBranch*> b_InterceptX_TelPlane_out;   //!
    std::vector<TBranch*> b_InterceptY_TelPlane_out;   //!
    std::vector<TBranch*> b_InterceptResX_TelPlane_out;   //!
    std::vector<TBranch*> b_InterceptResY_TelPlane_out;   //!
    std::vector<TBranch*> b_ClusterSizeX_TelPlane_out;   //!
    std::vector<TBranch*> b_ClusterSizeY_TelPlane_out;   //!
    std::vector<TBranch*> b_ClusterNumPixels_TelPlane_out;  //!

    //==================Vactors for assgning the branches to==================

    UInt_t  EvnNo;

    // 1st stage analysis variables and pointers
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
    std::vector<std::vector<double>> m_EvVAdjCFD;
    std::vector<std::vector<double>> m_EvCFDTime;
    std::vector<std::vector<double>> m_EvJitNdVdT;
    std::vector<std::vector<double>> m_EvDVDTCFD;
    std::vector<std::vector<double>> m_EvCFDToT;
    double m_FitDiffChi2;
    double m_FitDiffMag;
    double m_FitDiffMagErr;
    double m_FitDiffSigma;
    double m_FitDiffSigmaErr;
    double m_FitDiffTime;

    // 2nd stage analysis variables and pointers
    int m_chName;
    std::pair<int, int> m_SubChName;
    double m_MaxIndxFtChi2;
    double m_MinIndxFtChi2;
    double m_MaxVoltFtChi2;
    double m_MinVoltFtChi2;
    double m_StrIndxFtChi2;
    double m_EndIndxFtChi2;
    double m_NoiseFtChi2;
    double m_NoiseErrFtChi2;
    double m_PedestFtChi2;
    double m_PedestErrFtChi2;
    double m_MaxTimeFtChi2;
    double m_MinTimeFtChi2;
    double m_ChargeFtChi2;
    double m_RiseTimeFtChi2;
    double m_TriggTimeFtChi2;
    double m_DVDTMaxFtChi2;
    double m_TriggToTFtChi2;
    double m_SignalFFTFtChi2;
    double m_NoiseFFTFtChi2;
    std::pair <double, double> m_ChFitMxVlt;
    std::pair <double, double> m_ChFitMinVlt;
    std::pair <int, int> m_ChFitMxIndx;
    std::pair <int, int> m_ChFitMinIndx;
    std::pair <int, int> m_ChFitStrIndx;
    std::pair <int, int> m_ChFitEndIndx;
    std::pair <double, double> m_ChFitIsSig;
    std::pair <double, double> m_ChFitInWin;
    std::pair <double, double> m_ChFitVStr;
    std::pair <double, double> m_ChFitNoise;
    std::pair <double, double> m_ChFitPedestal;
    std::pair <double, double> m_ChFitNoiseSlp;
    std::pair <double, double> m_ChFitPedestalSlp;
    std::pair <double, double> m_ChFitNoiseErr;
    std::pair <double, double> m_ChFitPedestalErr;
    std::pair <double, double> m_ChFitMaxTime;
    std::pair <double, double> m_ChFitMinTime;
    std::pair <double, double> m_ChFitCharge;
    std::pair <double, double> m_ChFitRiseT;
    std::pair <double, double> m_ChFitJitRiseSNR;
    std::pair <double, double> m_ChFitTriggTime;
    std::pair <double, double> m_ChFitDVDTMax;
    std::pair <double, double> m_ChFitTriggToT;
    std::pair <double, double> m_ChFitSignalFFT;
    std::pair <double, double> m_ChFitNoiseFFT;
    std::pair <double, double> m_ChFitComplete;
    std::pair <double, double> m_ChFitSoNR;
    std::vector<std::pair<double, double>>* m_ChFitCFDTime;
    std::vector<std::pair<double, double>>* m_ChFitDVDTCFD;
    std::vector<std::pair<double, double>>* m_ChFitCFDToT;
    std::vector<std::pair<double, double>>* m_ChFitJitNdVdT;
    std::vector<double>* m_CFDTimeFtChi2;
    std::vector<double>* m_DVDTCFDFtChi2;
    std::vector<double>* m_CFDToTFtChi2;

    // 2nd stage waveshape analysis variables and pointers
    int m_ChShapeName;
    int m_SubChShName;
    int m_ChShapePolarity;
    int m_ChShapeMaxIndx;
    int m_ChShapeMinIndx;
    int m_ChShapeStrIndx;
    int m_ChShapeEndIndx;
    double m_ChShapeMax;
    double m_ChShapeMin;
    double m_ChShapeMaxTime;
    double m_ChShapeMinTime;
    double m_ChShapeCharge;
    double m_ChShapeRiseTime;
    double m_ChShapeJitRiseSNR;
    double m_ChShapeTriggTime;
    double m_ChShapeDVDTMax;
    double m_ChShapeTriggToT;
    double m_ChShapeSignalFFT;
    double m_ChShapeNoiseFFT;
    std::pair <double, double> m_ChShapeNoise;
    std::pair <double, double> m_ChShapeSoN;
    std::vector<double>* m_ChShapeCFDTime;
    std::vector<double>* m_ChShapeDVDTCFD;
    std::vector<double>* m_ChShapeCFDToT;
    std::vector<double>* m_ChShapeJitNdVdT;
    std::vector<std::pair <double, double>>* m_ChMeanVolt;

    // Third leevl analysis vectors
    std::vector<double>* m_TmDiffChi2;
    std::vector<double>* m_TmDiffFitMean;
    std::vector<double>* m_TmDiffFitErr;
    unsigned int m_CFD_ChA;
    unsigned int m_ChA;
    unsigned int m_ChB;

    // Histograms
    std::vector<TH1F*> m_trigDt_In;
    std::vector<TH1F*> m_trigFr_In;

    TH2D* m_FitMagVar;
    TH2D* m_FitMagErrVar;
    TH2D* m_FitSigmaVar;
    TH2D* m_FitSigmaErrVar;
    TH2D* m_FitChi2Var;
    std::vector<TH2D*>[5] m_FitVarIn;
    unsigned int m_NormCount[5];

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
    std::vector<std::vector<TH1D* >> h_TmDiffCFD;
    std::vector<std::vector<TCanvas* >> h_TmDiffCFDCanv;
    std::vector<TH2D* > h_TmSigmas;

    std::vector<DUTChannel*> m_RunDUTCh;
    std::vector<std::pair<unsigned int, std::vector<DUTChannel*>>> m_RunDUTSubCh;
    std::vector<WaveForm*> m_RunWaveform;
    WfmBase* m_RunBase;
    unsigned int m_init;

    // ClassDef(WfmRun, 1);
};

#endif