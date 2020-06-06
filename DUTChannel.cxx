/*
* DUTChannel.cxx
*
*
*      Author: Gkougkousis Evangelos - Leonidas
*              egkougko@cern.ch
*               IFAE-BARCELONA
*/

#include "LGADUtils/DUTChannel.h"

//#if !defined(__CINT__)
//ClassImp(DUTChannel);
//#endif

// Empty constructor
DUTChannel::DUTChannel()
{
    m_waveform = new WaveForm();
    ChannelInit();
}
// --------------------------------------------------------------------------------------------------------------
DUTChannel::DUTChannel(int channelID, TString name)
{
  m_waveform = new WaveForm();
  ChannelInit();
  SetChName(name);
  SetChID(channelID);
}
// --------------------------------------------------------------------------------------------------------------
DUTChannel::DUTChannel(int channelID)
{
    m_waveform = new WaveForm();
    ChannelInit();
    SetChID(channelID);
}
// --------------------------------------------------------------------------------------------------------------
DUTChannel::DUTChannel(TString name)
{
    m_waveform = new WaveForm();
    ChannelInit();
    SetChName(name);
}
// --------------------------------------------------------------------------------------------------------------
DUTChannel::DUTChannel(int channelID, TString name, unsigned int board, unsigned int ampl, float fraction, float trigger)
{
    m_waveform = new WaveForm();
    ChannelInit();
    SetChName(name);
    SetChID(channelID);
    SetBoard(board);
    if (board == 1 || board == 2 || board == 3) SetSecStage(none);
    else SetSecStage(ampl);
    SetChFract(fraction);
    SetChTrigg(trigger);
}
// --------------------------------------------------------------------------------------------------------------
DUTChannel::DUTChannel(LGADBase *tBase)
{
    m_waveform = new WaveForm(tBase);
    ChannelInit();
    m_instrument = tBase->LGADBase::GetInstrument();
    m_TrackComb = tBase->LGADBase::GetTrackComb();
    m_verbose = tBase->LGADBase::GetVerbosity();
    m_TransFile = tBase->LGADBase::GetTransFile();
    m_TransFileName = tBase->LGADBase::GetTransFileName();
    m_fitopt = tBase->LGADBase::GetFitMethode();
    m_WaveShape = tBase->LGADBase::GetWaveShape();
}
// --------------------------------------------------------------------------------------------------------------
DUTChannel::DUTChannel(int channelID, TString name, LGADBase *tBase)
{
    m_waveform = new WaveForm(tBase);
    ChannelInit();
    SetChName(name);
    SetChID(channelID);
    m_instrument = tBase->LGADBase::GetInstrument();
    m_TrackComb = tBase->LGADBase::GetTrackComb();
    m_verbose = tBase->LGADBase::GetVerbosity();
    m_TransFile = tBase->LGADBase::GetTransFile();
    m_TransFileName = tBase->LGADBase::GetTransFileName();
    m_fitopt = tBase->LGADBase::GetFitMethode();
    m_WaveShape = tBase->LGADBase::GetWaveShape();
}
// --------------------------------------------------------------------------------------------------------------
DUTChannel::DUTChannel(int channelID, LGADBase *tBase)
{
    m_waveform = new WaveForm(tBase);
    ChannelInit();
    SetChID(channelID);
    m_instrument = tBase->LGADBase::GetInstrument();
    m_TrackComb = tBase->LGADBase::GetTrackComb();
    m_verbose = tBase->LGADBase::GetVerbosity();
    m_TransFile = tBase->LGADBase::GetTransFile();
    m_TransFileName = tBase->LGADBase::GetTransFileName();
    m_fitopt = tBase->LGADBase::GetFitMethode();
    m_WaveShape = tBase->LGADBase::GetWaveShape();
}
// --------------------------------------------------------------------------------------------------------------
DUTChannel::DUTChannel(TString name, LGADBase *tBase)
{
    m_waveform = new WaveForm(tBase);
    ChannelInit();
    SetChName(name);
    m_instrument = tBase->LGADBase::GetInstrument();
    m_TrackComb = tBase->LGADBase::GetTrackComb();
    m_verbose = tBase->LGADBase::GetVerbosity();
    m_TransFile = tBase->LGADBase::GetTransFile();
    m_TransFileName = tBase->LGADBase::GetTransFileName();
    m_fitopt = tBase->LGADBase::GetFitMethode();
    m_WaveShape = tBase->LGADBase::GetWaveShape();
}
// --------------------------------------------------------------------------------------------------------------
DUTChannel::DUTChannel(int channelID, TString name, unsigned int board, unsigned int ampl, float fraction, float trigger, LGADBase *tBase)
{
    m_waveform = new WaveForm(tBase);
    ChannelInit();
    SetChName(name);
    SetChID(channelID);
    SetBoard(board);
    if (board == 1 || board == 2 || board == 3) SetSecStage(none);
    else SetSecStage(ampl);
    SetChFract(fraction);
    SetChTrigg(trigger);
    m_instrument = tBase->LGADBase::GetInstrument();
    m_TrackComb = tBase->LGADBase::GetTrackComb();
    m_verbose = tBase->LGADBase::GetVerbosity();
    m_TransFile = tBase->LGADBase::GetTransFile();
    m_TransFileName = tBase->LGADBase::GetTransFileName();
    m_fitopt = tBase->LGADBase::GetFitMethode();
    m_WaveShape = tBase->LGADBase::GetWaveShape();
}
// --------------------------------------------------------------------------------------------------------------
bool DUTChannel::SetChName(TString name)
{
    if (!(name.IsNull())) { m_channelName = name; return true; }
    else {
          std::cout << __FUNCTION__ << " WARNING: Channel name cannot be empty" << std::endl;
          return false;
         }
}
// --------------------------------------------------------------------------------------------------------------
bool DUTChannel::SetChID(int channelID)
{
    if (channelID >= 0) 
       {
        m_channelID = channelID; 
        m_waveform->WaveForm::SetWaveId(channelID);
        return true; 
       }
    else {
          std::cout << __FUNCTION__ << " WARNING: Channel ID not valid, reventing to 0" << std::endl;
          return false;
         }
}
// --------------------------------------------------------------------------------------------------------------
bool DUTChannel::SetChRate(Long64_t rate)
{
    if (rate > 0)
       { 
        m_ChRate = rate; 
        m_waveform->WaveForm::SetSnRate(m_ChRate);
        return true; 
       }
    else {
          std::cout << __FUNCTION__ << " WARNING: Invalid sampling rate value!" << std::endl;
          return false;
         }
}
// --------------------------------------------------------------------------------------------------------------
bool DUTChannel::SetChTransImp(float transimp)
{
    if (transimp > 0) 
       {
        m_ChTransimp = transimp; 
        m_waveform->WaveForm::SetTransimp(transimp);
        return true; 
       }
    else {
          std::cout << __FUNCTION__ << " WARNING: Transimpedence value not valid, reventing to default: " << m_ChTransimp << std::endl;
          return false;
         }
}
// --------------------------------------------------------------------------------------------------------------
bool DUTChannel::SetChAmpGain(int gain)
{
    if (gain >= 1) 
       { 
        m_ChAmpgain = gain; 
        m_waveform->WaveForm::SetAmpGain(gain);
        return true; 
       }
    else {
          std::cout << __FUNCTION__ << " WARNING: Amplifier gain not valid, reventing to default value: " << m_ChAmpgain << std::endl;
          return false;
         }
}
// --------------------------------------------------------------------------------------------------------------
bool DUTChannel::SetChFract(float fraction)
{
    if (fraction < 1 && fraction > 0) 
       { 
        m_ChFraction = fraction; 
        m_waveform->WaveForm::SetCFDfraction(fraction);
        return true; 
       }
    else {
          std::cout << __FUNCTION__ << " WARNING: CFD value not valid, reventing to default: " << m_ChFraction << std::endl;
          return false;
         }
}
// --------------------------------------------------------------------------------------------------------------
bool DUTChannel::SetChTrigg(double trigger)
{
    if (trigger > 0 && trigger < 1.5) 
       { 
        m_ChTrigg = trigger; 
        m_waveform->WaveForm::SetTrigg(trigger);
        return true; 
       }
    else {
          std::cout << __FUNCTION__ << " WARNING: Trigger value out of acceptable range, reventing to default: " << m_ChTrigg << std::endl;
          return false;
         }
}
// --------------------------------------------------------------------------------------------------------------
void DUTChannel::SetBoard(unsigned int board)
{   
    if (board == 0) { m_board = SingleCh; SetChTransImp(470); }
    else if (board == 1) { m_board = FourCh; SetChTransImp(2700); }
    else if (board == 2) { m_board = IN2P3; SetChTransImp(500); }
    else if (board == 3) { m_board = KU; SetChTransImp(400); }
    else if (board == 4) { m_board = SiPM; SetChTransImp(50); }
    else {
          m_board = SingleCh;
          std::cout << "ERROR: board type " << board 
                    << " not allowed. Allowed values are SingleCh, FourCh, IN2P3, SiPM and KU. Falling back to " 
                    << m_board << std::endl;
         } 
}
// --------------------------------------------------------------------------------------------------------------
void DUTChannel::SetSecStage(unsigned int ampl)
{
    if (ampl == 0) { m_amp = MinCircuits; SetChAmpGain(10); }
    else if (ampl == 1) { m_amp = Particulars; SetChAmpGain(60); }
    else if (ampl == 2) { m_amp = CIVIDEC; SetChAmpGain(100); }
    else if (ampl == 3) { m_amp = none; SetChAmpGain(1); }
    else {
          m_amp = MinCircuits;
          std::cout << "ERROR: second stagne mplifier " << ampl
                    << " not allowed. Allowed values are MinCircuits, Particulars, CIVIDEC and none. Falling back to "
                    << m_amp << std::endl;
        }
}
// --------------------------------------------------------------------------------------------------------------
// Destructor
DUTChannel::~DUTChannel()
{
    delete m_waveform;
    if (m_ChFitMaxIndxChi2 != -99 && m_ChFitMaxIndxChi2 != -1) delete m_ChFitMaxIndxFt;
    if (m_ChFitMaxVoltChi2 != -99 && m_ChFitMaxVoltChi2 != -1) delete m_ChFitMaxVoltFt;
    if (m_ChFitMinIndxChi2 != -99 && m_ChFitMinIndxChi2 != -1) delete m_ChFitMinIndxFt;
    if (m_ChFitMinVoltChi2 != -99 && m_ChFitMinVoltChi2 != -1) delete m_ChFitMinVoltFt;
    if (m_ChFitStrIndxChi2 != -99 && m_ChFitStrIndxChi2 != -1) delete m_ChFitStrIndxFt;
    if (m_ChFitEndIndxChi2 != -99 && m_ChFitEndIndxChi2 != -1) delete m_ChFitEndIndxFt;
    if (m_ChNoiseErrChi2 != -99 && m_ChNoiseErrChi2 != -1) delete m_ChNoiseErrFt;
    if (m_ChNoiseChi2 != -99 && m_ChNoiseChi2 != -1) delete m_ChNoiseFt;
    if (m_ChPedestalChi2 != -99 && m_ChPedestalChi2 != -1) delete m_ChPedestalFt;
    if (m_ChPedestalErChi2 != -99 && m_ChPedestalErChi2 != -1) delete m_ChPedestalErFt;
    if (m_ChMaxTimeChi2 != -99 && m_ChMaxTimeChi2 != -1) delete m_ChMaxTimeFt;
    if (m_ChMinTimeChi2 != -99 && m_ChMinTimeChi2 != -1) delete m_ChMinTimeFt;
    if (m_ChChargeChi2 != -99 && m_ChChargeChi2 != -1) delete m_ChChargeFt;
    if (m_ChFitRiseTimeChi2 != -99 && m_ChFitRiseTimeChi2 != -1) delete m_ChFitRiseTimeFt;
    for (unsigned int i = 0; i < 19; i++) { if (m_ChCFDTimeChi2[i] != -99 && m_ChCFDTimeChi2[i] != -1) delete m_ChCFDTimeFt[i]; }
    if (m_ChTriggTimeChi2 != -99 && m_ChTriggTimeChi2 != -1) delete m_ChTriggTimeFt;
    if (m_ChDVDTMaxChi2 != -99 && m_ChDVDTMaxChi2 != -1) delete m_ChDVDTMaxFt;
    if (m_ChDVDTCFDChi2 != -99 && m_ChDVDTCFDChi2 != -1) delete m_ChDVDTCFDFt;
    if (m_ChCFDToTChi2 != -99 && m_ChCFDToTChi2 != -1) delete m_ChCFDToTFt;
    if (m_ChTriggToTChi2 != -99 && m_ChTriggToTChi2 != -1) delete m_ChTriggToTFt;
    if (m_ChSignalFFTChi2 != -99 && m_ChSignalFFTChi2 != -1) delete m_ChSignalFFTFt;
    if (m_ChNoiseFFTChi2 != -99 && m_ChNoiseFFTChi2 != -1) delete m_ChNoiseFFTFt;

}
// --------------------------------------------------------------------------------------------------------------
std::pair<double, double> DUTChannel::GetChValue(int med, double ChMean, double ChMeanErr, std::pair<double, double> ChFit)
{
    std::pair<double, double> mag = std::make_pair(-99., -99.);
    if (med == 1)
       {
        if (ChMean == -99 || ChMeanErr == -99)
           {
            int err = updateChProperties();
            if (err != 0) std::cout << __FUNCTION__ << " ERROR: Failed to calculate channel shape properties!" << std::endl;
            else { mag.first = ChMean; mag.second = ChMeanErr;}
           }
        else {
              mag.first = ChMean;
              mag.second = ChMeanErr;
             }
       }
    else if (med == 2)
            {
             if (ChFit.first == -99 || ChFit.second == -99)
                {
                 int err = updateChProperties();
                 if (err != 0 || err!= -3) std::cout << __FUNCTION__ << " ERROR: Failed to calculate channel properties!" << std::endl;
                 else mag = ChFit;
                }
             else mag = ChFit;
            }
    else std::cout << __FUNCTION__ << " ERROR: Invalid methode for properties estimation, returning -99" << std::endl;
    return mag;
}
// --------------------------------------------------------------------------------------------------------------
double DUTChannel::GetFitChi2(std::string fit, int i)
{

    if (fit == "MaxIndx") return m_ChFitMaxIndxChi2;
    else if (fit == "MaxVolt") return m_ChFitMaxVoltChi2;
    else if (fit == "MinIndx") return m_ChFitMinIndxChi2;
    else if (fit == "MinVotl") return m_ChFitMinVoltChi2;
    else if (fit == "Strindx") return m_ChFitStrIndxChi2;
    else if (fit == "EndIndx") return m_ChFitEndIndxChi2;
    else if (fit == "NoisErr") return m_ChNoiseErrChi2;
    else if (fit == "Noise") return m_ChNoiseChi2;
    else if (fit == "Pedestl") return m_ChPedestalChi2;
    else if (fit == "PdslErr") return m_ChPedestalErChi2;
    else if (fit == "MaxTime") return m_ChMaxTimeChi2;
    else if (fit == "MinTime") return m_ChMinTimeChi2;
    else if (fit == "Charge") return m_ChChargeChi2;
    else if (fit == "RiseTime") return m_ChFitRiseTimeChi2;
    else if (fit == "TriggTime") return m_ChTriggTimeChi2;
    else if (fit == "DVDTMax") return m_ChDVDTMaxChi2;
    else if (fit == "DVDTCFD") return m_ChDVDTCFDChi2;
    else if (fit == "CFDToT") return m_ChCFDToTChi2;
    else if (fit == "TrigToT") return m_ChTriggToTChi2;
    else if (fit == "SigFFT") return m_ChSignalFFTChi2;
    else if (fit == "NoiseFFT") return m_ChNoiseFFTChi2;
    else if (fit == "CFDTime") return m_ChCFDTimeChi2.at(i);
    else {
          std::cout << __FUNCTION__ << " ERROR : Please use one of the appropriate options: " << std::endl
                    << " \t \t " << " MaxVolt \tMinIndx \tMinVotl  \tStrindx\tEndIndx\tMaxIndx" << std::endl
                    << " \t \t " << " NoisErrt\tNoise   \tPedestl  \tPdslErr\tMaxTime\tMinTime" << std::endl
                    << " \t \t " << " Charget \tRiseTime\tTriggTime\tDVDTMax\tDVDTCFD\tCFDToT" << std::endl
                    << " \t \t " << " TrigToTt\tSigFFT  \tNoiseFFT \tCFDTime" << std::endl;
          return -99.;
         }

}
// --------------------------------------------------------------------------------------------------------------
TH1D* DUTChannel::GetFit(std::string fit, int i)
{

    if (fit == "MaxIndx") return (TH1D*)m_ChFitMaxIndxFt->Clone();
    else if (fit == "MaxVolt") return (TH1D*)m_ChFitMaxVoltFt->Clone();
    else if (fit == "MinIndx") return (TH1D*)m_ChFitMinIndxFt->Clone();
    else if (fit == "MinVotl") return (TH1D*)m_ChFitMinVoltFt->Clone();
    else if (fit == "Strindx") return (TH1D*)m_ChFitStrIndxFt->Clone();
    else if (fit == "EndIndx") return (TH1D*)m_ChFitEndIndxFt->Clone();
    else if (fit == "NoisErr") return (TH1D*)m_ChNoiseErrFt->Clone();
    else if (fit == "Noise") return (TH1D*)m_ChNoiseFt->Clone();
    else if (fit == "Pedestl") return (TH1D*)m_ChPedestalFt->Clone();
    else if (fit == "PdslErr") return (TH1D*)m_ChPedestalErFt->Clone();
    else if (fit == "MaxTime") return (TH1D*)m_ChMaxTimeFt->Clone();
    else if (fit == "MinTime") return (TH1D*)m_ChMinTimeFt->Clone();
    else if (fit == "Charge") return (TH1D*)m_ChChargeFt->Clone();
    else if (fit == "RiseTime") return (TH1D*)m_ChFitRiseTimeFt->Clone();
    else if (fit == "TriggTime") return (TH1D*)m_ChTriggTimeFt->Clone();
    else if (fit == "DVDTMax") return (TH1D*)m_ChDVDTMaxFt->Clone();
    else if (fit == "DVDTCFD") return (TH1D*)m_ChDVDTCFDFt->Clone();
    else if (fit == "CFDToT") return (TH1D*)m_ChCFDToTFt->Clone();
    else if (fit == "TrigToT") return (TH1D*)m_ChTriggToTFt->Clone();
    else if (fit == "SigFFT") return (TH1D*)m_ChSignalFFTFt->Clone();
    else if (fit == "NoiseFFT") return (TH1D*)m_ChNoiseFFTFt->Clone();
    else if (fit == "CFDTime") return (TH1D*)(m_ChCFDTimeFt.at(i))->Clone();
    else {
          std::cout << __FUNCTION__ << " ERROR : Please use one of the appropriate options: " << std::endl
                    << " \t " << " MaxVolt\tMinIndx\tMinVotl\tStrindx\tEndIndx\tEndIndx" << std::endl
                    << " \t " << " NoisErrt\tNoise\tPedestl\tPdslErr\tMaxTime\tMinTime" << std::endl
                    << " \t " << " Charget\tRiseTime\tTriggTime\tDVDTMax\tDVDTCFD\tCFDToT" << std::endl
                    << " \t " << " TrigToTt\tSigFFT\tNoiseFFT\tCFDTime" << std::endl;
          return NULL;
         }

}
// --------------------------------------------------------------------------------------------------------------
std::pair<int, int> DUTChannel::GetChMaxIndx(int methode)
{
  return GetChValue(methode, m_ChMeanMaxIndx, 0.0, m_ChFitMaxIndx);
}
// --------------------------------------------------------------------------------------------------------------
std::pair<double, double> DUTChannel::GetChMaxVolt(int methode)
{
  return GetChValue(methode, m_ChMeanMax, m_chMeanMaxErr, m_ChFitMaxVolt);
}
// --------------------------------------------------------------------------------------------------------------
std::pair<int, int> DUTChannel::GetChMinIndx(int methode)
{
 return GetChValue(methode, m_ChMeanMinIndx, 0.0, m_ChFitMinIndx);
}
// --------------------------------------------------------------------------------------------------------------
std::pair<double, double> DUTChannel::GetChMinVolt(int methode)
{
 return GetChValue(methode, m_ChMeanMin, 0.0, m_ChFitMinVolt);
}
// --------------------------------------------------------------------------------------------------------------
std::pair<double, double> DUTChannel::GetChIsSignal()
{
 return GetChValue(2, 0.0, 0.0, m_ChFitIsSignal);
}
// --------------------------------------------------------------------------------------------------------------
std::pair<double, double> DUTChannel::GetChInWindow()
{
 return GetChValue(2, 0.0, 0.0, m_ChFitIsInWindow);
}
// --------------------------------------------------------------------------------------------------------------
std::pair<double, double> DUTChannel::GetChVoltStr()
{
 return GetChValue(2, 0.0, 0.0, m_ChFitVoltStr);
}
// --------------------------------------------------------------------------------------------------------------
std::pair<int, int> DUTChannel::GetChStrIndx(int methode)
{
 return GetChValue(methode, (double)m_ChMeanStrIndx, 0.0, (pair<double, double>)m_ChFitStrIndx);
}
// --------------------------------------------------------------------------------------------------------------
std::pair<int, int> DUTChannel::GetChEndIndx(int methode)
{
 return GetChValue(methode, (double)m_ChMeanEndIndx, 0.0, (pair<double, double>)m_ChFitEndIndx);
}
// --------------------------------------------------------------------------------------------------------------
std::pair<double, double> DUTChannel::GetChNoise(int methode)
{
 return GetChValue(methode, m_ChMeanNoise, m_ChMeanNoiseErr, m_ChFitNoise);
}
// --------------------------------------------------------------------------------------------------------------
std::pair<double, double> DUTChannel::GetChPedestal()
{
 return GetChValue(2, 0.0, 0.0, m_ChFitPedestal);
}
// --------------------------------------------------------------------------------------------------------------
std::pair<double, double> DUTChannel::GetChNoiseSlope()
{
 return GetChValue(2, 0.0, 0.0, m_ChFitNoiseSlope);
}
// --------------------------------------------------------------------------------------------------------------
std::pair<double, double> DUTChannel::GetChPedestalSlope()
{
 return GetChValue(2, 0.0, 0.0, m_ChFitPedestalSlope);
}
// --------------------------------------------------------------------------------------------------------------
std::pair<double, double> DUTChannel::GetChNoiseErr(int methode)
{
 return GetChValue(2, m_ChMeanNoiseErr, 0.0, m_ChFitNoiseErr);
}
// --------------------------------------------------------------------------------------------------------------
std::pair<double, double> DUTChannel::GetChPedestalErr(int methode)
{
 return GetChValue(2, 0.0, 0.0, m_ChFitPedestalErr);
}
// --------------------------------------------------------------------------------------------------------------
std::pair<double, double> DUTChannel::GetChMaxTime(int methode)
{
 return GetChValue(methode, m_ChMeanMaxTime, 0.0, m_ChFitMaxTime);
}
// --------------------------------------------------------------------------------------------------------------
std::pair<double, double> DUTChannel::GetChMinTime(int methode)
{
 return GetChValue(methode, m_ChMeanMinTime, 0.0, m_ChFitMinTime);
}
// --------------------------------------------------------------------------------------------------------------
std::pair<double, double> DUTChannel::GetChCharge(int methode)
{
 return GetChValue(methode, m_ChMeanCharge, 0.0, m_ChFitCharge);
}
// --------------------------------------------------------------------------------------------------------------
std::pair<double, double> DUTChannel::GetChRiseTime(int methode)
{
 return GetChValue(methode, m_ChMeanRiseTime, 0.0, m_ChFitRiseTime);
}
// --------------------------------------------------------------------------------------------------------------
std::pair<double, double> DUTChannel::GetChJitNdVdT(int methode)
{
 return GetChValue(methode, m_ChMeanJitterNdVdT, 0.0, m_ChFitJitNdVdT);
}
// --------------------------------------------------------------------------------------------------------------
std::pair<double, double> DUTChannel::GetChJitRiseSNR(int methode)
{
 return GetChValue(methode, m_ChMeanJitterRiseSNR, 0.0, m_ChFitJitRiseSNR);
}
// --------------------------------------------------------------------------------------------------------------
std::vector<std::pair<double, double>> DUTChannel::GetChCFDTime(int methode)
{
  std::vector<std::pair<double, double>> func;
  func.clear();
  if (methode == 1)
     {
      if (m_ChMeanCFDTime.size() == 0)
         {
          int err = updateChProperties();
          if (err != 0) std::cout << __FUNCTION__ << " ERROR: Failed to calculate channel CFD!" << std::endl;
          else for (unsigned int h = 0; h < m_ChMeanCFDTime.size(); h++) func.push_back(std::make_pair(m_ChMeanCFDTime.at(h), -1));
         }
      else for (unsigned int h = 0; h < m_ChMeanCFDTime.size(); h++) func.push_back(std::make_pair(m_ChMeanCFDTime.at(h), -1));
     }
  else if (methode == 2)
          {
           if (m_FitChCFDTime.size() == 0)
              {
               int err = updateChProperties();
               if (err != 0 || err!= -3) std::cout << __FUNCTION__ << " ERROR: Failed to calculate channel CFD!" << std::endl;
               else func = m_FitChCFDTime;
              }
           else func = m_FitChCFDTime;
          }
  else std::cout << __FUNCTION__ << " ERROR: Invalid methode for CFD computation, returning empty vector" << std::endl;
  return func;
}
// --------------------------------------------------------------------------------------------------------------
std::pair<double, double> DUTChannel::GetChTriggTime(int methode)
{
 return GetChValue(methode, m_ChMeanTriggTime, 0.0, m_ChFitTriggTime);
}
// --------------------------------------------------------------------------------------------------------------
std::pair<double, double> DUTChannel::GetChDVDTMax(int methode)
{
 return GetChValue(methode, m_ChMeanDVDTMax, 0.0, m_ChFitDVDTMax);
}
// --------------------------------------------------------------------------------------------------------------
std::pair<double, double> DUTChannel::GetChDVDTCFD(int methode)
{
 return GetChValue(methode, m_ChMeanDVDTCFD, 0.0, m_ChFitDVDTCFD);
}
// --------------------------------------------------------------------------------------------------------------
std::pair<double, double> DUTChannel::GetChCFDToT(int methode)
{
 return GetChValue(methode, m_ChMeanCFDToT, 0.0, m_ChFitCFDToT);
}
// --------------------------------------------------------------------------------------------------------------
std::pair<double, double> DUTChannel::GetChTriggToT(int methode)
{
 return GetChValue(methode, m_ChMeanTriggToT, 0.0, m_ChFitTriggToT);
}
// --------------------------------------------------------------------------------------------------------------
std::pair<double, double> DUTChannel::GetChSignalFFT(int methode)
{
 return GetChValue(methode, m_ChMeanSignalFFT, 0.0, m_ChFitSignalFFT);
}
// --------------------------------------------------------------------------------------------------------------
std::pair<double, double> DUTChannel::GetChNoiseFFT(int methode)
{
 return GetChValue(methode, m_ChMeanNoiseFFT, 0.0, m_ChFitNoiseFFT);
}
// --------------------------------------------------------------------------------------------------------------
std::pair<double, double> DUTChannel::GetChComplete()
{
 return GetChValue(2, 0.0, 0.0, m_ChFitComplete);
}
// --------------------------------------------------------------------------------------------------------------
std::pair<double, double> DUTChannel::GetChPolarity()
{
return GetChValue(1, m_ChMeanPolarity, 0.0, std::make_pair(0.0, 0.0));
}
// --------------------------------------------------------------------------------------------------------------
std::pair<double, double> DUTChannel::GetChSoNR(int methode)
{
 return GetChValue(methode, m_ChMeanSoN.first, m_ChMeanSoN.second, m_ChFitSoN);
}
// --------------------------------------------------------------------------------------------------------------
int DUTChannel::GetEvntQuantInt(std::vector<int>* ChQuant)
{
    if (ChQuant->size() != 0) return ChQuant->back();
    else return -99;

}
// --------------------------------------------------------------------------------------------------------------
int DUTChannel::GetEvntQuantBool(std::vector<bool>* ChQuant)
{
    if (ChQuant->size() != 0 && ChQuant->back()) return 1;
    else if (ChQuant->size() != 0 && !ChQuant->back()) return 0;
    else return -99;
}
// --------------------------------------------------------------------------------------------------------------
double DUTChannel::GetEvntQuantDouble(std::vector<double>* ChQuant)
{
    if (ChQuant->size() != 0) return ChQuant->back();
    else return -99.;
}
// --------------------------------------------------------------------------------------------------------------
void DUTChannel::ChannelInit(/*unsigned int fevent, unsigned int points*/)
{
    // Constants
    SetChAmpGain(10);
    SetChTrigg(0.010);
    m_ChRate = -99.;
    SetChTransImp(470);
    SetChFract(0.2);
    m_channelID = -99;
    m_channelName = "n/a";
    m_ChNoPulses = 0;
    m_amp = MinCircuits;
    m_board = SingleCh;

    // Outwards variables for Mean pulse
    m_ChMeanPolarity = -99;
    m_ChMeanMaxIndx = -99;
    m_ChMeanMax = -99;
    m_ChMeanMinIndx = -99;
    m_ChMeanMin = -99;
    m_ChMeanStrIndx = -99;
    m_ChMeanEndIndx = -99;
    m_ChMeanMaxTime = -99;
    m_ChMeanMinTime = -99;
    m_ChMeanNoise = -99;
    m_ChMeanNoiseErr = -99;
    m_chMeanMaxErr = -99;
    m_ChMeanCharge = -99;
    m_ChMeanRiseTime = -99;
    m_ChMeanJitterNdVdT = -99;
    m_ChMeanJitterRiseSNR = -99;
    m_ChMeanTriggTime = -99;
    m_ChMeanDVDTMax = -99;
    m_ChMeanDVDTCFD = -99;
    m_ChMeanCFDToT = -99;
    m_ChMeanTriggToT = -99;
    m_ChMeanSignalFFT = -99;
    m_ChMeanNoiseFFT = -99;
    m_ChMeanSoN = std::make_pair(-99, -99);
    m_ChMeanCFDTime.clear();

    m_ChFitPedestalSlope = std::make_pair(-99, -99);
    m_ChFitNoiseSlope = std::make_pair(-99, -99);
    m_ChFitMaxIndx = std::make_pair(-99, -99);
    m_ChFitMaxVolt = std::make_pair(-99, -99);
    m_ChFitMinIndx = std::make_pair(-99, -99);
    m_ChFitMinVolt = std::make_pair(-99, -99);
    m_ChFitIsSignal = std::make_pair(-99, -99);
    m_ChFitIsInWindow = std::make_pair(-99, -99);
    m_ChFitVoltStr = std::make_pair(-99, -99);
    m_ChFitStrIndx = std::make_pair(-99, -99);
    m_ChFitEndIndx = std::make_pair(-99, -99);
    m_ChFitNoise = std::make_pair(-99, -99);
    m_ChFitPedestal = std::make_pair(-99, -99);
    m_ChFitNoiseErr = std::make_pair(-99, -99);
    m_ChFitPedestalErr = std::make_pair(-99, -99);
    m_ChFitMaxTime = std::make_pair(-99, -99);
    m_ChFitMinTime = std::make_pair(-99, -99);
    m_ChFitCharge = std::make_pair(-99, -99);
    m_ChFitRiseTime = std::make_pair(-99, -99);
    m_ChFitJitNdVdT = std::make_pair(-99, -99);
    m_ChFitJitRiseSNR = std::make_pair(-99, -99);
    m_FitChCFDTime.clear();
    m_ChFitTriggTime = std::make_pair(-99, -99);
    m_ChFitDVDTMax = std::make_pair(-99, -99);
    m_ChFitDVDTCFD = std::make_pair(-99, -99);
    m_ChFitCFDToT = std::make_pair(-99, -99);
    m_ChFitTriggToT = std::make_pair(-99, -99);
    m_ChFitSignalFFT = std::make_pair(-99, -99);
    m_ChFitNoiseFFT = std::make_pair(-99, -99);
    m_ChFitSoN = std::make_pair(-99, -99);
    m_ChFitComplete = std::make_pair(-99, -99);

    // Vectors
    m_ChVAdjCFD.clear();
    m_ChIsSignal.clear();
    m_ChIsInWindow.clear();
    m_ChPol.clear();
    m_ChVoltStr.clear();
    m_ChMaxIndx.clear();
    m_ChMaxTime.clear();
    m_ChMaxVolt.clear();
    m_ChMinIndx.clear();
    m_ChMinTime.clear();
    m_ChMinVolt.clear();
    m_ChStrIndx.clear();
    m_ChEndIndx.clear();
    m_ChNoise.clear();
    m_ChPedestal.clear();
    m_ChNoiseErr.clear();
    m_ChPedestalErr.clear();
    m_ChCharge.clear();
    m_ChRiseTime.clear();
    m_ChJitNdVdT.clear();
    m_ChJitRiseSNR.clear();
    m_ChHalfPts.clear();

    m_MeanChPlFt.clear();
    m_MeanChPlFtChi2.clear();
    m_MeanChPulse.clear();

    // Fits and Chi2 values
    m_ChCFDTimeChi2.clear();
    m_ChCFDTimeChi2.resize(19, -99);
    m_ChFitMaxIndxChi2 = -99;
    m_ChFitMaxVoltChi2 = -99;
    m_ChFitMinIndxChi2 = -99;
    m_ChFitMinVoltChi2 = -99;
    m_ChFitStrIndxChi2 = -99;
    m_ChFitEndIndxChi2 = -99;
    m_ChNoiseErrChi2 = -99;
    m_ChNoiseChi2 = -99;
    m_ChPedestalChi2 = -99;
    m_ChPedestalErChi2 = -99;
    m_ChMaxTimeChi2 = -99;
    m_ChMinTimeChi2 = -99;
    m_ChChargeChi2 = -99;
    m_ChFitRiseTimeChi2 = -99;
    m_ChDVDTMaxChi2 = -99;
    m_ChDVDTCFDChi2 = -99;
    m_ChCFDToTChi2 = -99;
    m_ChTriggToTChi2 = -99;
    m_ChSignalFFTChi2 = -99;
    m_ChNoiseFFTChi2 = -99;
    m_ChTriggTimeChi2 = -99;

    m_ChCFDTimeFt.clear();
    m_ChCFDTimeFt.resize(19, NULL);
    m_ChTriggTimeFt = NULL;
    m_ChDVDTMaxFt = NULL;
    m_ChDVDTCFDFt = NULL;
    m_ChCFDToTFt = NULL;
    m_ChTriggToTFt = NULL;
    m_ChSignalFFTFt = NULL;
    m_ChNoiseFFTFt = NULL;
    m_ChFitMaxIndxFt = NULL;
    m_ChFitMaxVoltFt = NULL;
    m_ChFitMinIndxFt = NULL;
    m_ChFitMinVoltFt = NULL;
    m_ChFitStrIndxFt = NULL;
    m_ChFitEndIndxFt = NULL;
    m_ChNoiseErrFt = NULL;
    m_ChNoiseFt = NULL;
    m_ChPedestalFt = NULL;
    m_ChPedestalErFt = NULL;
    m_ChMaxTimeFt = NULL;
    m_ChMinTimeFt = NULL;
    m_ChChargeFt = NULL;
    m_ChFitRiseTimeFt = NULL;

    for (unsigned int i = 0; i < 19; i++) (m_ChCFDTime[i]).clear(); 
    m_ChTriggTime.clear();
    m_ChDVDTMax.clear();
    m_ChDVDTCFD.clear();
    m_ChCFDToT.clear();
    m_ChTriggToT.clear();
    m_ChSignalFFT.clear();
    m_ChNoiseFFT.clear();
    m_ChComplete.clear();

}
// --------------------------------------------------------------------------------------------------------------
bool DUTChannel::AppendEvent(std::vector<double> *volt, std::vector<double>* time, bool waveshape)
{
    m_waveform->WaveForm::InitializeWaveForm(1);
    m_waveform->LGADBase::SetVerbose(m_verbose);
    m_waveform->WaveForm::SetVoltage(volt);
    m_waveform->WaveForm::SetTime(time);
    bool append = AppendWaveform(m_waveform);
    // Set waveform rate to channel rate
    if (m_waveform->WaveForm::GetSnRate() != m_ChRate) m_ChRate = m_waveform->WaveForm::GetSnRate();
    return append;
}
// --------------------------------------------------------------------------------------------------------------
bool DUTChannel::AppendWaveform(WaveForm* k)
{
    // Calculate Waveform properties
    k->WaveForm::Calculate();
    // Fill vectors
    if ((k->WaveForm::GetVoltage())->size() != 0 && (k->WaveForm::GetTime())->size() != 0)
       {           
        bool wfmvsat = k->WaveForm::GetSaturation();
        m_ChVoltStr.push_back(wfmvsat);
        m_ChPol.push_back(k->WaveForm::GetPolarity());
        if (!wfmvsat) 
           {
            m_ChMaxIndx.push_back(k->WaveForm::GetMaxIndx());
            m_ChMaxVolt.push_back(k->WaveForm::GetMax());
           }
        else {
              m_ChMaxIndx.push_back(-99);
              m_ChMaxVolt.push_back(-99);
             }
        m_ChMinIndx.push_back(k->WaveForm::GetMinIndx());
        m_ChMinVolt.push_back(k->WaveForm::GetMin());
        m_ChHalfPts = k->WaveForm::GetHalfPoints();
        m_ChIsSignal.push_back(k->WaveForm::GetIsSignal());
        if (k->WaveForm::GetIsSignal())
           {
            m_ChIsInWindow.push_back(k->WaveForm::GetIsInWidow());
            m_ChStrIndx.push_back(k->WaveForm::GetStrIndx());
            m_ChEndIndx.push_back(k->WaveForm::GetEndIndx());
            if (k->WaveForm::GetIsInWidow())
               {
                m_ChNoise.push_back(k->WaveForm::GetNoise());
                m_ChPedestal.push_back(k->WaveForm::GetPedestal());
                m_ChNoiseErr.push_back(k->WaveForm::GetNoiseErr());
                m_ChPedestalErr.push_back(k->WaveForm::GetPedestalErr());
                if (!wfmvsat) m_ChMaxTime.push_back(k->WaveForm::GetMaxTime());
                else m_ChMaxTime.push_back(-99);
                m_ChMinTime.push_back(k->WaveForm::GetMinTime());
                if (k->WaveForm::GetPolarity() != 2 && !wfmvsat)
                   {
                    m_ChCharge.push_back(k->WaveForm::GetCharge());
                    m_ChRiseTime.push_back(k->WaveForm::GetRiseTime());
                    m_ChSignalFFT.push_back(k->WaveForm::GetSignalFFT());
                    m_ChNoiseFFT.push_back(k->WaveForm::GetNoiseFFT());
                    m_ChJitRiseSNR.push_back(k->WaveForm::GetJitterRiseSNR());
                    m_ChDVDTMax.push_back(k->WaveForm::GetdVdTMax());
                    for (unsigned int i = 0; i < 19; i++) (m_ChCFDTime[i]).push_back(k->WaveForm::GetCFDTime(0.05*(i + 1)));
                    if (k->WaveForm::GetCFDTime() != -1)
                       {
                        m_ChVAdjCFD = k->WaveForm::GetTimeAdjVolt();
                        m_ChJitNdVdT.push_back(k->WaveForm::GetJitterNdVdT());
                        m_ChTriggTime.push_back(k->WaveForm::GetTriggTime());
                        m_ChDVDTCFD.push_back(k->WaveForm::GetdVdTCFD());
                        if (k->WaveForm::GetCFDToT() > 0.0)
                           {
                            m_ChCFDToT.push_back(k->WaveForm::GetCFDToT());
                            if (k->WaveForm::GetTriggToT() > 0)
                               {
                                m_ChTriggToT.push_back(k->WaveForm::GetTriggToT());
                                m_ChComplete.push_back(6);
                                m_ChNoPulses++;
                               }
                            else {
                                  FillDefault(6);
                                  m_ChComplete.push_back(5);
                                  if (m_verbose > 0) std::cout << __FUNCTION__ << " WARNING: Failed to calculate trigger CFD time: " 
                                                               << k->WaveForm::GetTriggToT() << " at " << k->WaveForm::GetTrigg() << std::endl;
                                  return true;
                                 }
                           }
                        else {
                              FillDefault(5);
                              m_ChComplete.push_back(4);
                              if (m_verbose > 0) std::cout << __FUNCTION__ << " WARNING: Failed to calculate ToT CFD time: " 
                                                           << k->WaveForm::GetCFDToT() << " at " <<  k->WaveForm::GetCFDfraction() << std::endl;
                              return true;
                             }
                       }
                    else {
                          FillDefault(4);
                          m_ChComplete.push_back(3);
                          if (m_verbose > 0) std::cout << __FUNCTION__ << " WARNING: Failed to calculate CFD time" << k->WaveForm::GetCFDTime()
                                                       << " for standard set value!" << std::endl;
                          return true;
                         }
                   }
                else {
                      FillDefault(3);
                      m_ChComplete.push_back(2);
                      if (m_verbose > 0) std::cout << __FUNCTION__ << " WARNING: Failed to calculate polarity [" << k->WaveForm::GetPolarity() 
                                                   << "] or signal saturated [" << wfmvsat << "] for cahnnel [" << m_channelID 
                                                   << "] -> Stopping at Min/Max estimation!" << std::endl;
                       return true;
                     }
               }
            else {
                  FillDefault(2);
                  m_ChComplete.push_back(1);
                  if (m_verbose > 0) std::cout << __FUNCTION__ << " WARNING: Signal out of window for channel [" << m_channelID 
                                               << "] -> Stopping at Min/Max estimation!" << std::endl;
                  return true;
                 }
           }
        else {
              FillDefault(1);
              m_ChComplete.push_back(0);
              if (m_verbose > 0) std::cout << __FUNCTION__ << " WARNING: No signal found in channel [" << m_channelID 
                                           << "] -> Stopping at Min/Max estimation!" << std::endl;
              return true;
             }
       }
    else {
          if (m_verbose > 0) std::cout << __FUNCTION__ << " WARNING: No pulses for this channel[" 
                                       << m_channelID << "] -> will not append properties!" << std::endl;
          return false;
         }
  return true;
}
// --------------------------------------------------------------------------------------------------------------
bool DUTChannel::FillDefault(int k)
{
    if (k < 2) 
       {
        m_ChIsInWindow.push_back(1);
        m_ChStrIndx.push_back(-99);
        m_ChEndIndx.push_back(-99);
       }
    if (k < 3)
       {
        m_ChNoise.push_back(-99);
        m_ChPedestal.push_back(-99);
        m_ChNoiseErr.push_back(-99);
        m_ChPedestalErr.push_back(-99);
        m_ChMaxTime.push_back(-99);
        m_ChMinTime.push_back(-99);
       }
    if (k < 4)
       {
        m_ChCharge.push_back(-99);
        m_ChRiseTime.push_back(-99);
        m_ChSignalFFT.push_back(-99);
        m_ChNoiseFFT.push_back(-99);
        m_ChJitRiseSNR.push_back(-99);
        m_ChDVDTMax.push_back(-99);
        for (unsigned int l = 0; l < 19; l++) (m_ChCFDTime[l]).push_back(-99);
       }
    if (k < 5)
       {
        m_ChVAdjCFD = vector<double>();
        m_ChJitNdVdT.push_back(-99);
        m_ChTriggTime.push_back(-99);
        m_ChDVDTCFD.push_back(-99);
       }
    if (k < 6) m_ChCFDToT.push_back(-99);
    m_ChTriggToT.push_back(-99);

    return true;
}
// --------------------------------------------------------------------------------------------------------------
int DUTChannel::CalculateShape(TTree* wavetree)
{

    // Find events are available for full calculation
    unsigned int events = 0;
    for (unsigned int a = 0; a < m_ChComplete.size(); a++) { if (m_ChComplete.at(a) > 3) events++; }

    // Protection if you could not find enougph complete events
    if (events < 10)
       {
        std::cout << __FUNCTION__ << " ERROR: Too few cpoplete pulses for this channel["
                  << events << "] -> will not perform waveform shape calculation!" << std::endl;
        return -1;
       }

    // Check if the tree exists
    if (wavetree == NULL)
       {
        std::cout << __FUNCTION__ << " ERROR: Pointer to event tree is empty, no shape calculation possible!" << std::endl;
        return -2;
       }

    unsigned int apoints = 2;
    std::vector<double> pointvolt;
    pointvolt.reserve(events);

    std::vector<double> *VAdjCFDEv;
    TBranch             *b_VAdjCFDEv;    //!
    VAdjCFDEv = 0;                       // Don't remove this or the previous line!!!!
    wavetree->SetBranchAddress(Form("VAdjCFD%02u", m_channelID), &VAdjCFDEv, &b_VAdjCFDEv);
    std::cout << __FUNCTION__ << " INFO : Mean shape calculation for channel " << m_channelID << " using " << events << " events..." << std::endl;
    for (unsigned int k = 0; k < (apoints-1); k++)
        {
         for (Long64_t ievt = 0; ievt < wavetree->GetEntries(); ievt++)
             {
              b_VAdjCFDEv->GetEntry(ievt);
              if (apoints == 2)
                 {
                  // Aassume all events have same amount of points
                  if (m_ChComplete.at(ievt) > 3) apoints = VAdjCFDEv->size();
                  else continue;
                  // ignore first point, set to zero from afjustment
                  m_MeanChPulse.resize(apoints - 1, std::make_pair(-99, -99));
                  m_MeanChPlFt.resize(apoints - 1, NULL);
                  m_MeanChPlFtChi2.resize(apoints-1, -99);
                 }
              if (m_ChComplete.at(ievt) > 3) pointvolt.push_back(VAdjCFDEv->at(k+1));
             }
         if (m_fitopt == "rootInt") m_MeanChPulse.at(k) = CalculateMeanGLandau(&pointvolt, &m_ChComplete, m_MeanChPlFt[k], m_MeanChPlFtChi2[k], 3, "root");
         else m_MeanChPulse.at(k) = CalculateMeanGLandau(&pointvolt, &m_ChComplete, m_MeanChPlFt[k], m_MeanChPlFtChi2[k], 3, m_fitopt);
         LGADBase::ProgressBar(k, (apoints-1));
         pointvolt.clear();
        }

    // Check if the average waveform is computed correctly
    if (m_MeanChPulse.size() == (apoints-1))
       {
        // Fill the time and voltage vectors
        std::vector<double> meanvolt;
        std::vector<double> meanvoltErr;
        meanvolt.reserve(apoints);
        meanvoltErr.reserve(apoints);
        meanvolt.push_back(0.0);
        meanvoltErr.push_back(0.0);
        for (unsigned int p = 0; p < m_MeanChPulse.size(); p++)
            {
             meanvolt.push_back((m_MeanChPulse.at(p)).first);
             meanvoltErr.push_back((m_MeanChPulse.at(p)).second);
            }
        // Define a waveformn and calculate all caracteristics for the mean pulse
        m_waveform->WaveForm::InitializeWaveForm(1);
        m_waveform->LGADBase::SetVerbose(m_verbose);
        m_waveform->WaveForm::SetVoltage(&meanvolt);
        m_waveform->WaveForm::SetSnRate(m_ChRate);
        m_waveform->WaveForm::Calculate();
        // Save mean pulse quantites, Individual pulses already cleaned and noise substracted
        // no point in saving adjasted voltages and pedestal from here!!
        m_ChMeanPolarity = m_waveform->WaveForm::GetPolarity();
        m_ChMeanMaxIndx = m_waveform->WaveForm::GetMaxIndx();
        m_ChMeanMax = m_waveform->WaveForm::GetMax();
        m_ChMeanMinIndx = m_waveform->WaveForm::GetMinIndx();
        m_ChMeanMin = m_waveform->WaveForm::GetMin();
        if (m_waveform->WaveForm::GetIsSignal())
           {
            m_ChMeanStrIndx = m_waveform->WaveForm::GetStrIndx();
            m_ChMeanEndIndx = m_waveform->WaveForm::GetEndIndx();
            if (m_waveform->WaveForm::GetPolarity() != 2)
               {
                m_ChMeanNoise = m_waveform->WaveForm::GetNoise();
                m_ChMeanNoiseErr = m_waveform->WaveForm::GetNoiseErr();
                m_ChMeanMaxTime = m_waveform->WaveForm::GetMaxTime();
                m_chMeanMaxErr = m_waveform->WaveForm::GetMaxErr();
                m_ChMeanMinTime = m_waveform->WaveForm::GetMinTime();
                m_ChMeanCharge = m_waveform->WaveForm::GetCharge();
                m_ChMeanRiseTime = m_waveform->WaveForm::GetRiseTime();
                m_ChMeanJitterNdVdT = m_waveform->WaveForm::GetJitterNdVdT();
                m_ChMeanJitterRiseSNR = m_waveform->WaveForm::GetJitterRiseSNR();
                m_ChMeanSoN = calculateSoN(m_ChMeanNoise, m_ChMeanNoiseErr, m_ChMeanMax, m_chMeanMaxErr);
                for (unsigned int i = 1; i < 20; i++) m_ChMeanCFDTime.push_back(m_waveform->WaveForm::GetCFDTime(0.05 * i));
                m_ChMeanTriggTime = m_waveform->WaveForm::GetTriggTime();
                m_ChMeanDVDTMax = m_waveform->WaveForm::GetdVdTMax();
                m_ChMeanDVDTCFD = m_waveform->WaveForm::GetdVdTCFD();
                m_ChMeanCFDToT = m_waveform->WaveForm::GetCFDToT();
                m_ChMeanTriggToT = m_waveform->WaveForm::GetTriggToT();
                m_ChMeanSignalFFT = m_waveform->WaveForm::GetSignalFFT();
                m_ChMeanNoiseFFT = m_waveform->WaveForm::GetNoiseFFT();
               }
           }
        return 0;
       }
    else {
          std::cout << __FUNCTION__ << " ERROR : Average pulse shape calculation for channel["
                                    << m_channelID << "] failed! -> don't use mean properties!" << std::endl;
          return -3;
         }

}
// --------------------------------------------------------------------------------------------------------------
int DUTChannel::updateChProperties(bool waveshape, TTree* wavetree)
{
  // Value to be returned 
  int back = 0;

  if (m_ChNoPulses < 10)
     {
      std::cout << __FUNCTION__ << " ERROR : Less than 10 pulses for this channel ["
                << m_ChNoPulses << " valid] -> will not update properties!" << std::endl;
      return -1;
     }

  if (waveshape == true) CalculateShape(wavetree);
  else {
        std::cout << __FUNCTION__ << " WARNING : Did not calculate mean pulse shape for this channel["
                                  << m_channelID << "] -> don't use mean properties!" << std::endl;
        back = 1;
       }

  // Mean of pulse properties form Gaussian fits
  std::vector<double> ChMaxIndx(m_ChMaxIndx.begin(), m_ChMaxIndx.end());
  std::vector<double> ChMinIndx(m_ChMinIndx.begin(), m_ChMinIndx.end());
  std::vector<double> ChStrIndx(m_ChStrIndx.begin(), m_ChStrIndx.end());
  std::vector<double> ChEndIndx(m_ChEndIndx.begin(), m_ChEndIndx.end());
  m_ChFitMaxIndx = CalculateMeanGLandau(&ChMaxIndx, &m_ChComplete, m_ChFitMaxIndxFt, m_ChFitMaxIndxChi2, 2, m_fitopt);
  if (m_verbose >= 1) std::cout << __FUNCTION__ << " INFO : Channel fit Max Index: " << m_ChFitMaxIndx.first << " +/- "
                                << m_ChFitMaxIndx.second << ", Chi2: " << m_ChFitMaxIndxChi2 << std::endl;
  if (m_ChFitMaxIndxChi2 != -1) 
     {
      m_ChFitMaxIndxFt->SetNameTitle(Form("Or_MaxIndx%02u", m_channelID), Form("Or_MaxIndx%02u", m_channelID));
      m_ChFitMaxIndxFt->SetDirectory(nullptr); // detach histo from open directory
     }
  m_ChFitMaxVolt = CalculateMeanGLandau(&m_ChMaxVolt, &m_ChComplete, m_ChFitMaxVoltFt, m_ChFitMaxVoltChi2, 2, m_fitopt);
  if (m_verbose >= 1) std::cout << __FUNCTION__ << " INFO : Channel fit Max Voltage: " << m_ChFitMaxVolt.first << " +/- "
                                << m_ChFitMaxVolt.second << ", Chi2: " << m_ChFitMaxVoltChi2 << std::endl;
  if (m_ChFitMaxVoltChi2 != -1) 
     {
      m_ChFitMaxVoltFt->SetNameTitle(Form("Or_MaxVolt%02u", m_channelID), Form("Or_MaxVolt%02u", m_channelID));
      m_ChFitMaxVoltFt->SetDirectory(nullptr); // detach histo from open directory
     }
  m_ChFitMinIndx = CalculateMeanG(&ChMinIndx, &m_ChComplete, m_ChFitMinIndxFt, m_ChFitMinIndxChi2, 2, m_fitopt);
  if (m_verbose >= 1) std::cout << __FUNCTION__ << " INFO : Channel fit Min Index: " << m_ChFitMinIndx.first << " +/- "
                                << m_ChFitMinIndx.second << ", Chi2: " << m_ChFitMinIndxChi2 << std::endl;
  if (m_ChFitMinIndxChi2 != -1) 
     {
      m_ChFitMinIndxFt->SetNameTitle(Form("Or_MinIndx%02u", m_channelID), Form("Or_MinIndx%02u", m_channelID));
      m_ChFitMinIndxFt->SetDirectory(nullptr); // detach histo from open directory
      }
  m_ChFitMinVolt = CalculateMeanG(&m_ChMinVolt, &m_ChComplete, m_ChFitMinVoltFt, m_ChFitMinVoltChi2, 2, m_fitopt);
  if (m_verbose >= 1) std::cout << __FUNCTION__ << " INFO : Channel fit Min Voltage: " << m_ChFitMinVolt.first << " +/- "
                                << m_ChFitMinVolt.second << ", Chi2: " << m_ChFitMinVoltChi2 << std::endl;
  if (m_ChFitMinVoltChi2 != -1) 
     {
      m_ChFitMinVoltFt->SetNameTitle(Form("Or_MinVolt%02u", m_channelID), Form("Or_MinVolt%02u", m_channelID));
      m_ChFitMinVoltFt->SetDirectory(nullptr); // detach histo from open directory
     }
  m_ChFitStrIndx = CalculateMeanGLandau(&ChStrIndx, &m_ChComplete, m_ChFitStrIndxFt, m_ChFitStrIndxChi2, 2, m_fitopt);
  if (m_verbose >= 1) std::cout << __FUNCTION__ << " INFO : Channel fit Start Index: " << m_ChFitStrIndx.first << " +/- "
                                << m_ChFitStrIndx.second << ", Chi2: " << m_ChFitStrIndxChi2 << std::endl;
  if (m_ChFitStrIndxChi2 != -1) 
     {
      m_ChFitStrIndxFt->SetNameTitle(Form("Or_StrIndx%02u", m_channelID), Form("Or_StrIndx%02u", m_channelID));
      m_ChFitStrIndxFt->SetDirectory(nullptr); // detach histo from open directory
     }
  m_ChFitEndIndx = CalculateMeanGLandau(&ChEndIndx, &m_ChComplete, m_ChFitEndIndxFt, m_ChFitEndIndxChi2, 2, m_fitopt);
  if (m_verbose >= 1) std::cout << __FUNCTION__ << " INFO : Channel fit End Index: " << m_ChFitEndIndx.first << " +/- "
                                << m_ChFitEndIndx.second << ", Chi2: " << m_ChFitEndIndxChi2 << std::endl;
  if (m_ChFitEndIndxChi2 != -1) 
     {
      m_ChFitEndIndxFt->SetNameTitle(Form("Or_EndIndx%02u", m_channelID), Form("Or_EndIndx%02u", m_channelID));
      m_ChFitEndIndxFt->SetDirectory(nullptr); // detach histo from open directory
     }

  // Noise and Pedestal Calculations
  // Use linear fits to offset baseline increase and substrat th effect
  // Obly apply to sifnals within the aquisition window
  m_ChFitNoiseSlope = CalculateMeanLinear(&m_ChNoise, &m_ChIsInWindow, &m_ChNoiseErr);
  if (fabs(m_ChFitNoiseSlope.first) < 0.1)
     {
      std::cout << __FUNCTION__ << " INFO : Noise slope negligable ("<< fabs(m_ChFitNoiseSlope.first) <<"), applying Gaussian fit" << std::endl;
      m_ChFitNoiseErr = CalculateMeanG(&m_ChNoiseErr, &m_ChComplete, m_ChNoiseErrFt, m_ChNoiseErrChi2, 1, m_fitopt);
      if (m_ChNoiseErrChi2 != -1) 
         {
          m_ChNoiseErrFt->SetNameTitle(Form("Or_NoiseErr%02u", m_channelID), Form("Or_NoiseErr%02u", m_channelID));
          m_ChNoiseErrFt->SetDirectory(nullptr); // detach histo from open directory
         }
      m_ChFitNoise = CalculateMeanG(&m_ChNoise, &m_ChComplete, m_ChNoiseFt, m_ChNoiseChi2, 1, m_fitopt);
      if (m_ChNoiseChi2 != -1) 
         {
          m_ChNoiseFt->SetNameTitle(Form("Or_Noise%02u", m_channelID), Form("Or_Noise%02u", m_channelID));
          m_ChNoiseFt->SetDirectory(nullptr); // detach histo from open directory
         }
     }
  else {
        std::cout << __FUNCTION__ << " INFO : Significant noise slope detected ("<< fabs(m_ChFitNoiseSlope.first) 
                  <<"), applying Linear x Gaussian fit" << std::endl;
        m_ChFitNoiseErr = CalculateMeanGLinear(&m_ChNoiseErr, &m_ChComplete, m_ChNoiseErrFt, m_ChNoiseErrChi2, 1, m_fitopt);
        if (m_ChNoiseErrChi2 != -1) 
           {
            m_ChNoiseErrFt->SetNameTitle(Form("Or_NoiseErr%02u", m_channelID), Form("Or_NoiseErr%02u", m_channelID));
            m_ChNoiseErrFt->SetDirectory(nullptr); // detach histo from open directory
           }
        m_ChFitNoise = CalculateMeanGLinear(&m_ChNoise, &m_ChComplete, m_ChNoiseFt, m_ChNoiseChi2, 1, m_fitopt);
        if (m_ChNoiseChi2 != -1) 
           {
            m_ChNoiseFt->SetNameTitle(Form("Or_Noise%02u", m_channelID), Form("Or_Noise%02u", m_channelID));
            m_ChNoiseFt->SetDirectory(nullptr); // detach histo from open directory
           }
       }
  
  m_ChFitPedestalSlope = CalculateMeanLinear(&m_ChPedestal, &m_ChIsInWindow, &m_ChPedestalErr);
  if (fabs(m_ChFitPedestalSlope.first) < 0.1)
     {
      std::cout << __FUNCTION__ << " INFO : Pedestal slope negligable (" << fabs(m_ChFitPedestalSlope.first) << "), applying Gaussian fit" << std::endl;
      m_ChFitPedestal = CalculateMeanG(&m_ChPedestal, &m_ChComplete, m_ChPedestalFt, m_ChPedestalChi2, 1, m_fitopt);
      if (m_ChPedestalChi2 != -1) 
         {
          m_ChPedestalFt->SetNameTitle(Form("Or_Pedest%02u", m_channelID), Form("Or_Pedest%02u", m_channelID));
          m_ChPedestalFt->SetDirectory(nullptr); // detach histo from open directory
         }
      m_ChFitPedestalErr = CalculateMeanG(&m_ChPedestalErr, &m_ChComplete, m_ChPedestalErFt, m_ChPedestalErChi2, 1, m_fitopt);
      if (m_ChPedestalErChi2 != -1) 
         {
          m_ChPedestalErFt->SetNameTitle(Form("Or_PedestErr%02u", m_channelID), Form("Or_PedestErr%02u", m_channelID));
          m_ChPedestalErFt->SetDirectory(nullptr); // detach histo from open directory
         }
     }
  else {
        std::cout << __FUNCTION__ << " INFO : Significant pedestal slope detected ("<< fabs(m_ChFitNoiseSlope.first)
                                  <<"), applying Linear x Gaussian fit" << std::endl;
        m_ChFitPedestal = CalculateMeanGLinear(&m_ChPedestal, &m_ChComplete, m_ChPedestalFt, m_ChPedestalChi2, 1, m_fitopt);
        if (m_ChPedestalChi2 != -1) 
           {
            m_ChPedestalFt->SetNameTitle(Form("Or_Pedest%02u", m_channelID), Form("Or_Pedest%02u", m_channelID));
            m_ChPedestalFt->SetDirectory(nullptr); // detach histo from open directory
           }
        m_ChFitPedestalErr = CalculateMeanGLinear(&m_ChPedestalErr, &m_ChComplete, m_ChPedestalErFt, m_ChPedestalErChi2, 1, m_fitopt);
        if (m_ChPedestalErChi2 != -1) 
           {
            m_ChPedestalErFt->SetNameTitle(Form("Or_PedestErr%02u", m_channelID), Form("Or_PedestErr%02u", m_channelID));
            m_ChPedestalErFt->SetDirectory(nullptr); // detach histo from open directory
           }
       }

  m_ChFitMaxTime = CalculateMeanG(&m_ChMaxTime, &m_ChComplete, m_ChMaxTimeFt, m_ChMaxTimeChi2, 2, m_fitopt);
  if (m_verbose >= 1) std::cout << __FUNCTION__ << " INFO : Channel fit Max time: " << m_ChFitMaxTime.first << " +/- "
                                << m_ChFitMaxTime.second << ", Chi2: " << m_ChMaxTimeChi2 << std::endl;
  if (m_ChMaxTimeChi2 != -1) 
     {
      m_ChMaxTimeFt->SetNameTitle(Form("Or_MaxTime%02u", m_channelID), Form("Or_MaxTime%02u", m_channelID));
      m_ChMaxTimeFt->SetDirectory(nullptr); // detach histo from open directory
     }
  m_ChFitMinTime = CalculateMeanG(&m_ChMinTime, &m_ChComplete, m_ChMinTimeFt, m_ChMinTimeChi2, 2, m_fitopt);
  if (m_verbose >= 1) std::cout << __FUNCTION__ << " INFO : Channel fit Min time: " << m_ChFitMinTime.first << " +/- "
                                << m_ChFitMinTime.second << ", Chi2: " << m_ChMinTimeChi2 << std::endl;
  if (m_ChMinTimeChi2 != -1) 
     {
      m_ChMinTimeFt->SetNameTitle(Form("Or_MinTime%02u", m_channelID), Form("Or_MinTime%02u", m_channelID));
      m_ChMinTimeFt->SetDirectory(nullptr); // detach histo from open directory
     }
  m_ChFitCharge = CalculateMeanGLandau(&m_ChCharge, &m_ChComplete, m_ChChargeFt, m_ChChargeChi2, 2, m_fitopt);
  if (m_verbose >= 1) std::cout << __FUNCTION__ << " INFO : Channel fit Charge: " << m_ChFitCharge.first << " +/- "
                                << m_ChFitCharge.second << ", Chi2: " << m_ChChargeChi2 << std::endl;
  if (m_ChChargeChi2 != -1) 
     {
      m_ChChargeFt->SetNameTitle(Form("Or_Charge%02u", m_channelID), Form("Or_Charge%02u", m_channelID));
      m_ChChargeFt->SetDirectory(nullptr); // detach histo from open directory
     }
  m_ChFitRiseTime = CalculateMeanG(&m_ChRiseTime, &m_ChComplete, m_ChFitRiseTimeFt, m_ChFitRiseTimeChi2, 2, m_fitopt);
  if (m_verbose >= 1) std::cout << __FUNCTION__ << " INFO : Channel fit Rise time: " << m_ChFitRiseTime.first << " +/- "
                                << m_ChFitRiseTime.second << ", Chi2: " << m_ChFitRiseTimeChi2 << std::endl;
  if (m_ChFitRiseTimeChi2 != -1) 
     {
      m_ChFitRiseTimeFt->SetNameTitle(Form("Or_RiseT%02u", m_channelID), Form("Or_RiseT%02u", m_channelID));
      m_ChFitRiseTimeFt->SetDirectory(nullptr); // detach histo from open directory
     }

  std::cout << __FUNCTION__ << " INFO : Performing CFD fits for all bins.." << std::endl;
  for (unsigned int i = 0; i < 19; i++) 
      {
       LGADBase::ProgressBar(i, 19);
       m_FitChCFDTime.push_back(CalculateMeanG(&m_ChCFDTime[i], &m_ChComplete, m_ChCFDTimeFt[i], m_ChCFDTimeChi2[i], 2, m_fitopt));
       if (m_verbose >= 1) std::cout << __FUNCTION__ << " INFO : Channel fit CFD time at " << (i + 1) * 5 << "%: "
                                     << m_FitChCFDTime[i].first << " +/- " << m_FitChCFDTime[i].second
                                     << ", Chi2: " << m_ChCFDTimeChi2[i] << std::endl;
       if (m_ChCFDTimeChi2[i] != -1) 
          {
           m_ChCFDTimeFt[i]->SetNameTitle(Form("Or_CFDT%02u%02u", m_channelID, i), Form("Or_CFDT%02u%02u", m_channelID, i));
           m_ChCFDTimeFt[i]->SetDirectory(nullptr); // detach histo from open directory
          }
      }

  m_ChFitTriggTime = CalculateMeanG(&m_ChTriggTime, &m_ChComplete, m_ChTriggTimeFt, m_ChTriggTimeChi2, 3, m_fitopt);
  if (m_verbose >= 1) std::cout << __FUNCTION__ << " INFO : Channel fit Trigger time: " << m_ChFitTriggTime.first << " +/- "
                                << m_ChFitTriggTime.second << ", Chi2: " << m_ChTriggTimeChi2 << std::endl;
  if (m_ChTriggTimeChi2 != -1) 
     {
      m_ChTriggTimeFt->SetNameTitle(Form("Or_TriggT%02u", m_channelID), Form("Or_TriggT%02u", m_channelID));
      m_ChTriggTimeFt->SetDirectory(nullptr); // detach histo from open directory
     }
  m_ChFitDVDTMax = CalculateMeanGLandau(&m_ChDVDTMax, &m_ChComplete, m_ChDVDTMaxFt, m_ChDVDTMaxChi2, 2, m_fitopt);
  if (m_verbose >= 1) std::cout << __FUNCTION__ << " INFO : Channel fit Max dV/dT: " << m_ChFitDVDTMax.first << " +/- "
                                << m_ChFitDVDTMax.second << ", Chi2: " << m_ChDVDTMaxChi2 << std::endl;
  if (m_ChDVDTMaxChi2 != -1) 
     {
      m_ChDVDTMaxFt->SetNameTitle(Form("Or_DVDTMax%02u", m_channelID), Form("Or_DVDTMax%02u", m_channelID));
      m_ChDVDTMaxFt->SetDirectory(nullptr); // detach histo from open directory
     }
  m_ChFitDVDTCFD = CalculateMeanGLandau(&m_ChDVDTCFD, &m_ChComplete, m_ChDVDTCFDFt, m_ChDVDTCFDChi2, 3, m_fitopt);
  if (m_verbose >= 1) std::cout << __FUNCTION__ << " INFO : Channel fit dV/dT at CFD value: " << m_ChFitDVDTCFD.first << " +/- "
                                << m_ChFitDVDTCFD.second << ", Chi2: " << m_ChDVDTCFDChi2 << std::endl;
  if (m_ChDVDTCFDChi2 != -1) 
     {
      m_ChDVDTCFDFt->SetNameTitle(Form("Or_DVDTCFD%02u", m_channelID), Form("Or_DVDTCFD%02u", m_channelID));
      m_ChDVDTCFDFt->SetDirectory(nullptr); // detach histo from open directory
     }
  m_ChFitCFDToT = CalculateMeanG(&m_ChCFDToT, &m_ChComplete, m_ChCFDToTFt, m_ChCFDToTChi2, 4, m_fitopt);
  if (m_verbose >= 1) std::cout << __FUNCTION__ << " INFO : Channel fit CFD Time over Threshold: " << m_ChFitCFDToT.first << " +/- "
                                << m_ChFitCFDToT.second << ", Chi2: " << m_ChCFDToTChi2 << std::endl;
  if (m_ChCFDToTChi2 != -1) 
     {
      m_ChCFDToTFt->SetNameTitle(Form("Or_CFDToT%02u", m_channelID), Form("Or_CFDToT%02u", m_channelID));
      m_ChCFDToTFt->SetDirectory(nullptr); // detach histo from open directory
     }
  m_ChFitTriggToT = CalculateMeanG(&m_ChTriggToT, &m_ChComplete, m_ChTriggToTFt, m_ChTriggToTChi2, 5, m_fitopt);
  if (m_verbose >= 1) std::cout << __FUNCTION__ << " INFO : Channel fit Trigger Time over Threshold: " << m_ChFitTriggToT.first << " +/- "
                                << m_ChFitTriggToT.second << ", Chi2: " << m_ChTriggToTChi2 << std::endl;
  if (m_ChTriggToTChi2 != -1) 
     {
      m_ChTriggToTFt->SetNameTitle(Form("Or_TriggToT%02u", m_channelID), Form("Or_TriggToT%02u", m_channelID));
      m_ChTriggToTFt->SetDirectory(nullptr); // detach histo from open directory
     }
  m_ChFitSignalFFT = CalculateMeanG(&m_ChSignalFFT, &m_ChComplete, m_ChSignalFFTFt, m_ChSignalFFTChi2, 2, m_fitopt);
  if (m_verbose >= 1) std::cout << __FUNCTION__ << " INFO : Channel fit Signal FFT: " << m_ChFitSignalFFT.first << " +/- "
                                << m_ChFitSignalFFT.second << ", Chi2: " << m_ChSignalFFTChi2 << std::endl;
  if (m_ChSignalFFTChi2 != -1) 
     {
      m_ChSignalFFTFt->SetNameTitle(Form("Or_SignalFFT%02u", m_channelID), Form("Or_SignalFFT%02u", m_channelID));
      m_ChSignalFFTFt->SetDirectory(nullptr); // detach histo from open directory
     }
  m_ChFitNoiseFFT = CalculateMeanG(&m_ChNoiseFFT, &m_ChComplete, m_ChNoiseFFTFt, m_ChNoiseFFTChi2, 2, m_fitopt);
  if (m_verbose >= 1) std::cout << __FUNCTION__ << " INFO : Channel fit Noise FFT: " << m_ChFitNoiseFFT.first << " +/- "
                                << m_ChFitNoiseFFT.second << ", Chi2: " << m_ChNoiseFFTChi2 << std::endl;
  if (m_ChNoiseFFTChi2 != -1) 
     {
      m_ChNoiseFFTFt->SetNameTitle(Form("Or_NoiseFFT%02u", m_channelID), Form("Or_NoiseFFT%02u", m_channelID));
      m_ChNoiseFFTFt->SetDirectory(nullptr); // detach histo from open directory
     }
 
  m_ChFitIsSignal.first = LGADBase::Mean(&m_ChIsSignal);
  m_ChFitIsSignal.second = LGADBase::BayesianErr(&m_ChIsSignal, true);
  if (m_verbose >= 1) std::cout << __FUNCTION__ << " INFO : Channel fit isSignal: " << m_ChFitIsSignal.first << " +/- "
                                << m_ChFitIsSignal.second << std::endl;
  m_ChFitIsInWindow.first = LGADBase::Mean(&m_ChIsInWindow);
  m_ChFitIsInWindow.second = LGADBase::BayesianErr(&m_ChIsInWindow, true);
  if (m_verbose >= 1) std::cout << __FUNCTION__ << " INFO : Channel fit inWindow: " << m_ChFitIsInWindow.first << " +/- "
                                << m_ChFitIsInWindow.second << std::endl;
  m_ChFitVoltStr.first = LGADBase::Mean(&m_ChVoltStr);
  m_ChFitVoltStr.second = LGADBase::BayesianErr(&m_ChVoltStr, true);
  if (m_verbose >= 1) std::cout << __FUNCTION__ << " INFO : Channel fit isSaturated: " << m_ChFitVoltStr.first << " +/- "
                                <<  m_ChFitVoltStr.second << std::endl;
  m_ChFitComplete.first = (LGADBase::Mean(&m_ChComplete)) / 6;
  m_ChFitComplete.second = LGADBase::BayesianErr(&m_ChComplete, 6);
  if (m_verbose >= 1) std::cout << __FUNCTION__ << " INFO : Channel fit is Complete: " << m_ChFitComplete.first << " +/- "
                                << m_ChFitComplete.second << std::endl;

  // Signal over Noise ratio
  m_ChFitSoN = calculateSoN(m_ChFitNoise.first, m_ChFitNoiseErr.first, m_ChFitMaxVolt.first, m_ChFitMaxVolt.second);
  if (m_verbose >= 1) std::cout << __FUNCTION__ << " INFO : Channel fit SoN: " << m_ChFitSoN.first << " +/- "
                                << m_ChFitSoN.second << ", Chi2: " << m_ChNoiseFFTChi2 << std::endl;

  // Jitter calculation
  m_ChFitJitNdVdT.first = (m_ChFitNoise.first)/(m_ChFitDVDTCFD.first);
  m_ChFitJitNdVdT.second = sqrt(pow(m_ChFitNoise.second/m_ChFitDVDTCFD.first, 2) + pow((m_ChFitNoise.first/pow(m_ChFitDVDTCFD.first, 2))*m_ChFitDVDTCFD.second, 2));
  if (m_verbose >= 1) std::cout << __FUNCTION__ << " INFO : Channel fit Jitter (from dV/dt): " << m_ChFitJitNdVdT.first << " +/- "
                                << m_ChFitJitNdVdT.second << std::endl;
  m_ChFitJitRiseSNR.first = (m_ChFitRiseTime.first)*(m_ChFitMaxVolt.first) / (m_ChFitNoise.first);
  m_ChFitJitRiseSNR.second = sqrt(pow(m_ChFitRiseTime.second*m_ChFitNoise.first / m_ChFitMaxVolt.first, 2) +
                                  pow(m_ChFitRiseTime.first*m_ChFitNoise.second / m_ChFitMaxVolt.first, 2) + 
                                  pow((m_ChFitMaxVolt.second*m_ChFitRiseTime.first*m_ChFitNoise.first)/pow(m_ChFitNoise.first, 2), 2));
  if (m_verbose >= 1) std::cout << __FUNCTION__ << " INFO : Channel fit Jitter (from SNR): " << m_ChFitJitRiseSNR.first << " +/- "
                                << m_ChFitJitRiseSNR.second << std::endl;

  if (m_verbose >= 1) ChDump();

  return back;
}
// --------------------------------------------------------------------------------------------------------------
std::pair<double, double> DUTChannel::CalculateMeanGLandau(std::vector<double>* vec, std::vector<int>* signal, TH1D* &hitHist, double &fitChi2, int level, std::string fitopt)
{
    std::pair<double, double> VoltMPV = std::make_pair( -99, -99 );
    std::pair<double, double>  VoltSigma = std::make_pair( -99, -99 );
    std::vector<double> volt;
    unsigned int as = vec->size();
    int qual = -99;
    if (as == 0) std::cout << __FUNCTION__ << " EROOR : Cannot perform Gaussian x Landau fit on empty vector!!" << std::endl;
    else {
          volt.reserve(as);
          for (unsigned int k = 0; k < as; k++) if (signal->at(k) > level) volt.push_back(vec->at(k));
          if (fitopt == "RooFit") qual = LGADBase::RooConvFit(&volt, VoltMPV, VoltSigma, "LanXGau");
          else if (fitopt == "rootInt") qual = LGADBase::IterativeFit(&volt, VoltMPV, VoltSigma, hitHist, fitChi2, "LandauXGaussInt");
          else qual = LGADBase::IterativeFit(&volt, VoltMPV, VoltSigma, hitHist, fitChi2, "LandauXGauss");
          if (qual == 0) VoltMPV.second = VoltSigma.first;
          else {
                std::cout << __FUNCTION__ << " WARNING : Cannot perform Gaussian x Landau fit, reventing to Gaussian!!" << std::endl;
                VoltMPV = CalculateMeanG(vec, signal, hitHist, fitChi2, level, fitopt);
               }
         }
    return VoltMPV;
}    
// --------------------------------------------------------------------------------------------------------------
std::pair<double, double> DUTChannel::CalculateMeanG(vector<double>* vec, std::vector<int>* signal, TH1D* &hitHist, double &fitChi2, int level, std::string fitopt)
{
    std::pair<double, double> mean = std::make_pair(-99, -99 );
    std::pair<double, double> sigma = std::make_pair( -99, -99 );
    std::vector<double> real;
    unsigned int evnt = vec->size();
    int qual = -99;
    if (evnt == 0) std::cout << __FUNCTION__ << " EROOR : Cannot perform Gaussian fit on empty vector!!" << std::endl;
    else {
          real.reserve(evnt);
          for (unsigned int k = 0; k < evnt; k++) if (signal->at(k) > level && vec->at(k) != -1.) real.push_back(vec->at(k));
          if (fitopt == "rootInt") qual = LGADBase::IterativeFit(&real, mean, sigma, hitHist, fitChi2, "GaussInt");
          else qual = LGADBase::IterativeFit(&real, mean, sigma, hitHist, fitChi2, "Gauss");
          if (qual == 0) mean.second = sigma.first;
          else {
                mean.first = LGADBase::Mean(&real);
                mean.second = LGADBase::Stdev(&real);
                std::cout << __FUNCTION__ << " WARNING: Gaussian fit failed, reverting to mean and std!" << std::endl;
               }
         }
    return mean;
}
// --------------------------------------------------------------------------------------------------------------
std::pair<double, double> DUTChannel::CalculateMeanLinear(vector<double>* vec, std::vector<bool>* signal, vector<double>* vecErr)
{
    std::pair<double, double> slope = std::make_pair(-99, -99 );
    std::pair<double, double> intercept = std::make_pair( -99, -99 );
    std::vector<double> real;
    std::vector<double> realErr;
    unsigned int evnt = vec->size();
    if (evnt == 0) std::cout << __FUNCTION__ << " EROOR : Cannot perform Linear fit on empty vector!!" << std::endl;
    else {
          real.reserve(evnt);
          if (vecErr != NULL) realErr.reserve(evnt);
          for (unsigned int k = 0; k < evnt; k++) 
              {
               if (signal->at(k) == true) 
                  { 
                   if (vecErr != NULL) { real.push_back(vec->at(k)); realErr.push_back(vecErr->at(k)); }
                   else real.push_back(vec->at(k));
                  }
               }
          int res = -99;
          if (vecErr != NULL) res = LGADBase::LinearFit(&real, slope, intercept, &realErr);
          else res = LGADBase::LinearFit(&real, slope, intercept);
          if (res != 0) 
             {
              std::vector <double> deriv = LGADBase::Derivate(&real);
              slope.first = LGADBase::Mean(&deriv);
              slope.second = LGADBase::Stdev(&deriv);
              std::cout << __FUNCTION__ << " WARNING: Linear fit failed, reverting to derivative mean!" << std::endl;
             }
        }
    return slope;
}
// --------------------------------------------------------------------------------------------------------------
std::pair<double, double> DUTChannel::CalculateMeanGLinear(vector<double>* vec, std::vector<int>* signal, TH1D* &hitHist, double &fitChi2, int level, std::string fitopt)
{
    std::pair<double, double> mean = std::make_pair(-99, -99);
    std::pair<double, double> sigma = std::make_pair(-99, -99);
    std::vector<double> real;
    unsigned int evnt = vec->size();
    if (evnt == 0) std::cout << __FUNCTION__ << " EROOR : Cannot perform Linear fit on empty vector!!" << std::endl;
    else {
          real.reserve(evnt);
          for (unsigned int k = 0; k < evnt; k++) if (signal->at(k) > level) real.push_back(vec->at(k));
          int res = LGADBase::RooConvFit(&real, mean, sigma, "LinXGau") == 0;
          if (res != 0) 
             {
              std::cout << __FUNCTION__ << " WARNING : Cannot perform Linear x Gaussian fit reventing to Gaussian!!" << std::endl;
              mean = CalculateMeanG(vec, signal, hitHist, fitChi2, level, fitopt);
             }
        }
    return mean;
}
// --------------------------------------------------------------------------------------------------------------
std::pair<double, double>  DUTChannel::calculateSoN(double noise, double noiseErr, double VMax, double VMaxErr)
{
  std::pair<double, double> SNR = std::make_pair(-99, -99 );
  if (noise != -99 && VMax != -99)
     {
      SNR.first = VMax / noise;
      if (noiseErr != -99 && VMaxErr != -99) SNR.second = sqrt((pow((VMaxErr/noise), 2) + pow(((VMaxErr*noiseErr)/pow(noise, 2)), 2)));
      else std::cout << __FUNCTION__ << " WARNIRNG : Will not calculate SNR error, values are not provided!" << std::endl;
     }
  else std::cout << __FUNCTION__ << " ERROR : Cannot calculate SNR, values are not provided!" << std::endl;
  return SNR;
}
// --------------------------------------------------------------------------------------------------------------
int DUTChannel::ChDump()
{
  std::cout << "-----------------------------------------------------------------------------------------" << std::endl;
  std::cout << "==> Dumping channel info for channel " << m_channelName << " [" << m_channelID << "]" << std::endl;
  std::cout << "\t Amplifier gain:  " << m_ChAmpgain << " - set transimpedance: " << m_ChTransimp << std::endl
            << "\t CFD value:       " << m_ChFraction << " - trigger thershold:  " << m_ChTrigg << std::endl
            << "\t Sampling rate:   " << m_ChRate << std::endl << std::endl;

  if (m_ChNoPulses < 1)
     {
      std::cout << __FUNCTION__ << " ERROR : No pulses for this channel -> will not dump info!" << std::endl;
      return 1;
     }

  std::cout << "\t Events: " << m_ChPol.size() << " - events with pulses: " << m_ChNoPulses << " and polarity: " << m_ChMeanPolarity << std::endl 
            << "\t Pulse Shape max index            : " << m_ChMeanMaxIndx << ", max index from gausisan fit: " << m_ChFitMaxIndx.first << " +/- " << m_ChFitMaxIndx.second << std::endl 
            << "\t Pulse Shape amplitude            : " << m_ChMeanMax << " +/- " << m_chMeanMaxErr << ", MPV amplitude from gauss x Landau fit: " << m_ChFitMaxVolt.first << " +/- " << m_ChFitMaxVolt.second << std::endl 
            << "\t Pulse Shape min index            : " << m_ChMeanMinIndx << ", min index from gausisan fit: " << m_ChFitMinIndx.first << " +/- " << m_ChFitMinIndx.second << std::endl 
            << "\t Pulse Shape min value            : " << m_ChMeanMin << ", MPV min from gauss x Landau fit: " << m_ChFitMinVolt.first << " +/- " << m_ChFitMinVolt.second << std::endl 
            << "\t Pulse Shape start index          : " << m_ChMeanStrIndx << ", start index from gausisan fit: " << m_ChFitStrIndx.first << " +/- " << m_ChFitStrIndx.second << std::endl 
            << "\t Pulse Shape end index            : " << m_ChMeanEndIndx << ", end index from gausisan fit: " << m_ChFitEndIndx.first << " +/- " << m_ChFitEndIndx.second << std::endl 
            << "\t Pulse Shape noise                : " << m_ChMeanNoise << ", noise from Gausisan x Linear fit: " << m_ChFitNoise.first << " +/- " << m_ChFitNoise.second << std::endl 
            << "\t Pulse Shape noise uncertenty     : " << m_ChMeanNoiseErr << ", noise uncertenty from gausisan fit: " << m_ChFitNoiseErr.first << " +/- " << m_ChFitNoiseErr.second << std::endl 
            << "\t Pulse Shape max time             : " << m_ChMeanMaxTime << ", max time from gausisan fit: " << m_ChFitMaxTime.first << " +/- " << m_ChFitMaxTime.second << std::endl 
            << "\t Pulse Shape min time             : " << m_ChMeanMinTime << ", min time from gausisan fit: " << m_ChFitMinTime.first << " +/- " << m_ChFitMinTime.second << std::endl 
            << "\t Pulse Shape collected charge     : " << m_ChMeanCharge << ", collected charge MPV from Lamdau x Gaussian fit: " << m_ChFitCharge.first << " +/- " << m_ChFitCharge.second << std::endl;
  std::cout << "\t Pulse Shape rise time            : " << m_ChMeanRiseTime << ", rise time from Gaussian fit: " << m_ChFitRiseTime.first << " +/- " << m_ChFitRiseTime.second << std::endl 
            << "\t Pulse Shape jitter [RMS/(dV/dT)] : " << m_ChMeanJitterNdVdT << ", jitter [RMS/(dV/dT)] derived from Gaussian fit values: " << m_ChFitJitNdVdT.first << " +/- " << m_ChFitJitNdVdT.second << std::endl 
            << "\t Pulse Shape jitter [Rise/SNR]    : " << m_ChMeanJitterRiseSNR << ", jitter [Rise/SNR] derived from Gaussian fit values: " << m_ChFitJitRiseSNR.first << " +/- " << m_ChFitJitRiseSNR.second << std::endl 
            << "\t Pulse Shape Trigger ToA          : " << m_ChMeanTriggTime << ", Trigger ToA from Gaussian fit: " << m_ChFitTriggTime.first << " +/- " << m_ChFitTriggTime.second << std::endl 
            << "\t Pulse Shape max dV/dT            : " << m_ChMeanDVDTMax << ", max dV/dT from Gaussian fit: " << m_ChFitDVDTMax.first << " +/- " << m_ChFitDVDTMax.second << std::endl 
            << "\t Pulse Shape CFD dV/dT            : " << m_ChMeanDVDTCFD << ", CFD dV/dT from Gaussian fit: " << m_ChFitDVDTCFD.first << " +/- " << m_ChFitDVDTCFD.second << std::endl 
            << "\t Pulse Shape CFD ToT              : " << m_ChMeanCFDToT << ", CFD ToT from Gaussian fit: " << m_ChFitCFDToT.first << " +/- " << m_ChFitCFDToT.second << std::endl 
            << "\t Pulse Shape Trigger ToT          : " << m_ChMeanTriggToT << ", Trigger ToT from Gaussian fit: " << m_ChFitTriggToT.first << " +/- " << m_ChFitTriggToT.second << std::endl 
            << "\t Pulse Shape signal FFT           : " << m_ChMeanSignalFFT << ", signal FFT from Gaussian fit: " << m_ChFitSignalFFT.first << " +/- " << m_ChFitSignalFFT.second << std::endl 
            << "\t Pulse Shape noise FFT            : " << m_ChMeanNoiseFFT << ", noise FFT from Gaussian fit: " << m_ChFitNoiseFFT.first << " +/- " << m_ChFitNoiseFFT.second << std::endl
            << "\t Pulse Shape SNR                  : " << m_ChMeanSoN.first << " +/- " << m_ChMeanSoN.second << ", SNR derived from Gaussian fit values: " << m_ChFitSoN.first << " +/- " << m_ChFitSoN.second << std::endl;

  std::cout << "\t Channel is isgnal from gaussian fit              : " << m_ChFitIsSignal.first << " +/- " << m_ChFitIsSignal.second << std::endl 
            << "\t Channel is in window from gaussian fit           : " << m_ChFitIsInWindow.first << " +/- " << m_ChFitIsInWindow.second << std::endl
            << "\t Channel saturation from gaussian fit             : " << m_ChFitVoltStr.first << " +/- " << m_ChFitVoltStr.second << std::endl
            << "\t Channel noise slope from linear fit              : " << m_ChFitNoiseSlope.first << " +/- " << m_ChFitNoiseSlope.second << std::endl 
            << "\t Channel pedestal slope from linear fit           : " << m_ChFitPedestalSlope.first << " +/- " << m_ChFitPedestalSlope.second << std::endl 
            << "\t Channel pedestal mean from Linear x Gaussian fit : " << m_ChFitPedestal.first << " +/- " << m_ChFitPedestal.second << std::endl 
            << "\t Channel pedestal error from linear fit           : " << m_ChFitPedestalErr.first << " +/- " << m_ChFitPedestalErr.second << std::endl 
            << "\t Channel is complete from gaussian fit            : " << m_ChFitComplete.first << " +/- " << m_ChFitComplete.second << std::endl;
  std::cout << "-----------------------------------------------------------------------------------------" << std::endl;

  return 0;
}
// --------------------------------------------------------------------------------------------------------------
