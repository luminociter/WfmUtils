/*
* Waveform.cxx
*
*
*      Author: Gkougkousis Evangelos - Leonidas
*              egkougko@cern.ch
*               IFAE-BARCELONA
*/

#include "LGADUtils/WaveForm.h"

//#if !defined(__CINT__)
//ClassImp(WaveForm);
//#endif

// Default Empty constructor
WaveForm::WaveForm()
{
    m_fitchi2 = -99;
    m_WvBase = new LGADBase();
    InitializeWaveForm(m_WvBase);
}
// --------------------------------------------------------------------------------------------------------------
// Constructor with time and voltage input
WaveForm::WaveForm(std::vector<double>* voltage, std::vector<double>* time)
{
    m_fitchi2 = -99;
    m_WvBase = new LGADBase();
    InitializeWaveForm(m_WvBase);
    SetVoltage(voltage);
    SetTime(time);
    SetSnRate((long long int)ceil((Long64_t)(m_time->size() / abs(m_time->back() - m_time->front()))));
}
// --------------------------------------------------------------------------------------------------------------
// Constructor with rate and voltage input (time vector is computed internaly)
WaveForm::WaveForm(std::vector<double>* voltage, Long64_t snrate)
{
    m_fitchi2 = -99;
    m_WvBase = new LGADBase();
    InitializeWaveForm(m_WvBase);
    SetVoltage(voltage);
    SetSnRate(snrate);
}
// --------------------------------------------------------------------------------------------------------------
// Constructor with time, voltage and sampling rate input
WaveForm::WaveForm(std::vector<double>* voltage, std::vector<double>* time, Long64_t snrate)
{
    m_fitchi2 = -99;
    m_WvBase = new LGADBase();
    InitializeWaveForm(m_WvBase);
    SetVoltage(voltage);
    SetTime(time);
    SetSnRate(snrate);
}
// --------------------------------------------------------------------------------------------------------------
WaveForm::WaveForm(LGADBase *Base)
{
    m_fitchi2 = -99;
    m_WvBase = Base;
    InitializeWaveForm(m_WvBase);
}
// --------------------------------------------------------------------------------------------------------------
// Constructor from LGADBase with time and voltage input
WaveForm::WaveForm(std::vector<double>* voltage, std::vector<double>* time, LGADBase *Base)
{
    m_fitchi2 = -99;
    m_WvBase = Base;
    InitializeWaveForm(m_WvBase);
    SetVoltage(voltage);
    SetTime(time);
    SetSnRate((long long int)ceil((Long64_t)(m_time->size() / abs(m_time->back() - m_time->front()))));
}
// --------------------------------------------------------------------------------------------------------------
// Constructor from LGADBase with rate and voltage input (time vector is computed internaly)
WaveForm::WaveForm(std::vector<double>* voltage, Long64_t snrate, LGADBase *Base)
{
    m_fitchi2 = -99;
    m_WvBase = Base;
    InitializeWaveForm(m_WvBase);
    SetVoltage(voltage);
    SetSnRate(snrate);
}
// --------------------------------------------------------------------------------------------------------------
// Constructor from LGADBase with time, voltage and sampling rate input
WaveForm::WaveForm(std::vector<double>* voltage, std::vector<double>* time, Long64_t snrate, LGADBase *Base)
{
    m_fitchi2 = -99;
    m_WvBase = Base;
    InitializeWaveForm(m_WvBase);
    SetVoltage(voltage);
    SetTime(time);
    SetSnRate(snrate);
}
// --------------------------------------------------------------------------------------------------------------
// Default destructor
WaveForm::~WaveForm()
{
    delete m_WvBase;
    if (m_fitchi2 != -99 && m_fitchi2 != -1) delete m_noiseFit;
}
// --------------------------------------------------------------------------------------------------------------
void WaveForm::InitializeWaveForm(LGADBase* tBase, int level)
{
     m_instrument = tBase->LGADBase::GetInstrument();
     m_TrackComb = tBase->LGADBase::GetTrackComb();
     m_TransFile = tBase->LGADBase::GetTransFile();
     m_TransFileName = tBase->LGADBase::GetTransFileName();
     m_fitopt = tBase->LGADBase::GetFitMethode();
     m_WaveShape = tBase->LGADBase::GetWaveShape();

     // Pointers
     m_voltage = NULL;
     m_time = NULL;

     // Vectors
     m_voltageAdj.clear();
     m_vAdjCFD.clear();
     m_halfPtsIndx.clear();
     m_intime.clear();

     //Boleans
     m_IsSignal = false;
     m_IsInWindow = false;
     m_VoltStr = false;

     // Enums
     m_basepos = undet;
     m_pol = undef;

     // Re-Initiallizable variables
     m_maxIndx = -99;
     m_maxTime = -99;
     m_maxVolt = -99;
     m_maxVoltErr = -99;
     m_minIndx = -99;
     m_minTime = -99;
     m_minVolt = -99;
     m_StrIndx = -99;
     m_EndIndx = -99;
     m_maxavg = -99;
     m_minavg = -99;

     // Single initialization varibales
     m_noise = -99;
     m_pedestal = -99;
     m_noiseErr = -99;
     m_pedestalErr = -99;

     m_charge = -99;
     m_RiseTime = -99;
     m_jitter1 = -99;
     m_jitter2 = -99;
     m_CFDTime = -99;
     m_TriggTime = -99;
     m_dVdTMax = -99;
     m_dVdTCFD = -99;
     m_CFDToT = -99;
     m_TriggToT = -99;
     m_signalFFT = -99;
     m_noiseFFT = -99;
     m_NoiseFtQl = -99;
    if (m_fitchi2 != -99 && m_fitchi2 != -1) delete m_noiseFit;
     m_noiseFit = NULL;

     // Used values in the code
     if (level == 0)
        {
         m_fraction = 0.2;
         m_trigg = 0.010;
         m_transimp = 470;
         m_ampgain = 10;
         m_SnRate = -99;
         m_waveId = -99;
        }
}
// --------------------------------------------------------------------------------------------------------------
bool WaveForm::Calculate()
{

    // for (unsigned int srtr = 0; srtr < m_voltage->size(); srtr++) std::cout << m_time->at(srtr) << "\t " << m_voltage->at(srtr) << std::endl;
    if (m_voltage == NULL || m_time == NULL)
       {
        std::cout << __FUNCTION__ << " ERROR: Time and Voltage not set for waveform, cannot perform any operation!" << std::endl;
        return false;
       }
    m_pol = FindPolarity(m_voltage);
    if (m_WvBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Polarity: " << m_pol << std::endl;
    m_maxIndx = VoltMaxIndx(m_voltage, false, false);
    if (m_WvBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Max Voltage Index: " << m_maxIndx << std::endl;
    m_maxVolt = VoltMax(m_voltage, false, false);
    if (m_WvBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Max Voltage: " << m_maxVolt << std::endl;
    m_minIndx = VoltMinIndx(m_voltage, false, false);
    if (m_WvBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Min Voltage Index: " << m_minIndx << std::endl;
    m_minVolt = VoltMin(m_voltage, false, false);
    if (m_WvBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Min Voltage: " << m_minVolt << std::endl;
    m_halfPtsIndx = FindHalfPoints(m_voltage, false);
    if (m_halfPtsIndx.size() < 1)
       {
        if (m_WvBase->LGADBase::GetVerbosity() >= 1) std::cout << __FUNCTION__ << " WARNING: Could not determine half points, stopping..." << std::endl;
        return false;
       }
    m_IsSignal = IsSignal(&m_halfPtsIndx);    
    if (!m_IsSignal)
       {
        if (m_WvBase->LGADBase::GetVerbosity() >= 1) std::cout << __FUNCTION__<< " WARNING: Found several reocuring signals, classified as noise, stopping...." << std::endl;
        m_IsInWindow = true;
        return false;
       }
    if (m_WvBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: No. of Half points: " << m_halfPtsIndx.size() << " and is signal: " << m_IsSignal << std::endl;
    m_StrIndx = StartIndx(m_voltage, false, false);
    if (m_WvBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Start index: " << m_StrIndx << std::endl;
    m_EndIndx = EndIndx(m_voltage, false, false);
    if (m_WvBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: End index : " << m_EndIndx << std::endl;
    m_IsInWindow = InWindow(m_voltage);
    if (!m_IsInWindow)
       {
        if (m_WvBase->LGADBase::GetVerbosity() >= 1) std::cout << __FUNCTION__ << " WARNING: Signal ouside of DAQ window, stopping...." << std::endl;
        return false;
       }
    if (m_WvBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Pulse is within window: " << m_IsInWindow << std::endl;
    std::pair <unsigned int, unsigned int> noisepts = FindNoisePoints(m_voltage);
    if ((int)(noisepts.second - noisepts.first) < 4)
       {    
        if (m_WvBase->LGADBase::GetVerbosity() >= 1) std::cout << __FUNCTION__ << " WARNING: Too few noise points determined, stopping...." << std::endl;
        return false;
       }
    if (m_WvBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Using range for noise from " << noisepts.first << " to " << noisepts.second << std::endl;
    std::pair <double, double> baseline (-1., -1.);
    std::pair <double, double> noise_rms (-1., -1.);
    if (m_fitopt == "rootInt") m_NoiseFtQl = m_WvBase->LGADBase::IterativeFit(m_voltage, baseline, noise_rms, m_noiseFit, m_fitchi2, "GaussInt", noisepts, true);
    else m_NoiseFtQl = m_WvBase->LGADBase::IterativeFit(m_voltage, baseline, noise_rms, m_noiseFit, m_fitchi2, "Gauss", noisepts, true);
    if (m_fitchi2 != -1) 
       {
        m_noiseFit->SetNameTitle(Form("Or_WavNoise%02u", m_waveId), Form("Or_WavNoise%02u", m_waveId));
        m_noiseFit->SetDirectory(nullptr); // detach histo from open directory
       }
    if (m_NoiseFtQl == 0)
       {
        m_pedestal = baseline.first;
        m_pedestalErr = baseline.second;
        m_noise = noise_rms.first;
        m_noiseErr = noise_rms.second;
       }
    else {
          m_pedestal = m_WvBase->LGADBase::Mean(m_voltage, noisepts.first, noisepts.second);
          m_noise = m_WvBase->LGADBase::Stdev(m_voltage, noisepts.first, noisepts.second);
          if (m_WvBase->LGADBase::GetVerbosity() >= 1) std::cout << __FUNCTION__ << " WARNING: Could not calculate GAUSSIAN noise, reventing to mean and stdv...." << std::endl;
         }
    if (m_WvBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Noise: " << m_noise << " +/- " << m_noiseErr << " and pedestal: " 
                                                           << m_pedestal << " +/- " << m_pedestalErr << std::endl;
    m_VoltStr = VoltSatur(m_voltage, false);
    if (m_WvBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Saturation test: " << m_VoltStr << std::endl;
    // Adjust the pulse to positive, remove pedestal and recalculate max, min values and start, stop points
    if (m_pol == pos && !m_VoltStr) m_voltageAdj = PulseAdj(m_voltage, m_pedestal, 1);
    else if (m_pol == neg && !m_VoltStr) m_voltageAdj = PulseAdj(m_voltage, m_pedestal, -1);
    else if (m_pol == undef && !m_VoltStr)
            {
             if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Undetermined polarity for adjustment, stopping...." << std::endl;
             m_maxVoltErr = VoltMaxErr(m_voltage, false);
             m_minTime = TimeMin(m_voltage, m_SnRate, false, false);
             m_maxTime = TimeMax(m_voltage, m_SnRate, false, false);
             return false;
            }
    else if (m_VoltStr) 
            {
             if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Saturated signal, stopping...." << std::endl;
             m_minTime = TimeMin(m_voltage, m_SnRate, false, false);
             return false;
            }
    if (m_WvBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Pulse has been adjasted " << m_voltageAdj.size() << std::endl;
    m_RiseTime = RiseTimeLinear(&m_voltageAdj, m_SnRate);
    if (m_WvBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Rise time calculation " << m_RiseTime << std::endl;
    m_CFDTime = CFDTimeLinear(&m_voltageAdj, m_time, m_fraction);
    if (m_WvBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: CFD time calculation " << m_CFDTime << std::endl;
    m_vAdjCFD = PulseTimeVoltAdj(m_time, &m_voltageAdj, m_fraction);
    if (m_WvBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: TimeWalk corrected pulse: " << m_vAdjCFD.size() << std::endl;
    m_dVdTMax = dVdTMaxLinear(&m_voltageAdj, m_SnRate);
    if (m_WvBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: dV/dT max linear: " << m_dVdTMax << std::endl;
    m_dVdTCFD = dVdTCFDLinear(&m_voltageAdj, m_SnRate, m_fraction);
    if (m_WvBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: dV/dT CFD linear: " << m_dVdTCFD << std::endl;
    m_TriggTime = FirstTimeForVoltage(&m_voltageAdj, m_time, m_trigg);
    if (m_WvBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Trigger time: " << m_TriggTime << std::endl;
    m_CFDToT = CFDToTLinear(&m_voltageAdj, m_time, m_fraction);
    if (m_WvBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: CFD ToT: " << m_CFDToT << std::endl;
    m_TriggToT = TriggToTLinear(&m_voltageAdj, m_time, m_trigg);
    if (m_WvBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Triger ToT: " << m_TriggToT << std::endl;
    m_signalFFT = LGADBase::FFT(&m_voltageAdj, m_SnRate, m_StrIndx, m_EndIndx);
    if (m_WvBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Signal FFT: " << m_signalFFT << std::endl;
    m_noiseFFT = LGADBase::FFT(&m_voltageAdj, m_SnRate, noisepts.first, noisepts.second);
    if (m_WvBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Noise FFT: " << m_noiseFFT << std::endl;
    if (m_signalFFT != -99 && m_WvBase->LGADBase::GetDoTrnsCorr()) m_transimp = GetTransimpFromFile(m_signalFFT);
    m_charge = CollectedCharge(&m_voltageAdj, m_SnRate, m_transimp, m_ampgain, m_StrIndx, m_EndIndx);
    if (m_WvBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Collected Charge: " << m_charge << std::endl;
    m_jitter1 = m_noise / m_dVdTCFD;
    if (m_WvBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Jitter from Noise/dVdT: " << m_charge << std::endl;
    m_jitter2 = m_RiseTime*m_noise / m_maxVolt;
    if (m_WvBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Jitter from RieTime/pMax: " << m_charge << std::endl;
    if (m_WvBase->LGADBase::GetVerbosity() >= 3) dump();

    return true;
}
// --------------------------------------------------------------------------------------------------------------
void WaveForm::SetVoltage(std::vector<double> *volt)
{
    if (volt->size() != 0 && volt != NULL)
       {
        m_voltage = volt;
        if (m_time == NULL && m_SnRate > 0) m_time = FillTime(m_voltage, m_SnRate);
       }
    else std::cout << __FUNCTION__ << " ERROR: Empty voltage vector provided!!" << std::endl;
}
// --------------------------------------------------------------------------------------------------------------
void WaveForm::SetTime(std::vector<double>* time)
{
    if (time->size() != 0) 
       {
        m_time = time;
        // Calculate rate and compare it to set rate
        Long64_t rate = ceil((long long int)((time->size()) / abs(m_time->front() - m_time->back())));
        if (abs(rate - m_SnRate) / rate > 0.01) m_SnRate = rate; 
       }
    else std::cout << __FUNCTION__ << " ERROR: Empty time vector provided!!" << std::endl;
}
// --------------------------------------------------------------------------------------------------------------
void WaveForm::SetPolarity(unsigned int pol)
{
    if (pol == 0) m_pol = pos;
    else if (pol == 1) m_pol = neg;
    else if (pol == 2) m_pol = undef;
    else {
          std::cout << __FUNCTION__ << " WARNING: Polarity has none of the accepted values, setting it to undefined!" << std::endl;
          m_pol = undef;
         }
}
// --------------------------------------------------------------------------------------------------------------
void WaveForm::SetSnRate(Long64_t snrate)
{
    if (snrate <= 0) std::cout << __FUNCTION__ << " ERROR: Unacceptable sampling rate value " << snrate << "!!" << std::endl;
    if (m_time!= NULL)
       {
        int npoints = m_time->size() - 1;
        if (npoints < 0) 
           {
            std::cout << __FUNCTION__ << " WARNING: Rate test failed !!" << std::endl;
            m_SnRate = snrate;
           }
        else {
              Long64_t rate = ceil((long long int)(npoints / abs(m_time->front() - m_time->back())));
              if (abs(rate - snrate)/ rate > 0.01)
                 {
                  std::cout << __FUNCTION__ << " WARNING: Trying to set rate " << snrate << " is in disagreement with vector size and aq. time! " ;
                  std::cout << "Correct rate " << rate << ", will be set" << std::endl;
                  std::cout << npoints << " " << m_time->front() << " " << m_time->back() << " " << abs(m_time->front() - m_time->back()) << " " 
                            << npoints / abs(m_time->front() - m_time->back()) << std::endl;
                  m_SnRate = rate;
                 }
              else m_SnRate = snrate;
             }
       }
    else {
          m_SnRate = snrate;
          m_time = FillTime(m_voltage, m_SnRate);
         }
}
// --------------------------------------------------------------------------------------------------------------
std::vector<double>* WaveForm::FillTime(std::vector<double>* voltage, Long64_t snrate)
{
    std::vector<double>* time = NULL;
    if (voltage != NULL)
      { 
       m_intime.clear();
       m_intime.reserve(voltage->size());
       for (unsigned int b = 0; b < voltage->size(); b++) m_intime.push_back((double)b / (double)snrate);
       time = &m_intime;
      }
   else { 
         if (m_WvBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " WARNING: Rate set but no voltage vector found. Cannot establsh time vector, remember to do it before calculations!" 
                                                               << std::endl;
        }
   return time;
}
// --------------------------------------------------------------------------------------------------------------
void WaveForm::SetCFDfraction(double fraction)
{
    if (fraction <= 0 || fraction > 1)
       {
        std::cout << __FUNCTION__ << " ERROR: Unacceptable CFD value!!" << std::endl;
       }
    else m_fraction = fraction;
}
// --------------------------------------------------------------------------------------------------------------
void WaveForm::SetTrigg(double tri)
{
    m_trigg = tri;
}
// --------------------------------------------------------------------------------------------------------------
void WaveForm::SetBasePos(unsigned int pos)
{
    if (pos == 0) m_basepos = poside;
    else if (pos == 1) m_basepos = negside;
    else if (pos == 2) m_basepos = zerocross;
    else if (pos == 3) m_basepos = undet; 
    else {
          std::cout << __FUNCTION__ << " WARNING: Base position has none of the accepted values, seting it undefined!" << std::endl;
          m_basepos = undet;
         }
}
// --------------------------------------------------------------------------------------------------------------
void WaveForm::SetTransimp(float transimp)
{
    if (transimp > 0) m_transimp = transimp;
    else std::cout << __FUNCTION__ << " WARNING: Transimpedence value not correct!!" << std::endl;
}
// --------------------------------------------------------------------------------------------------------------
void WaveForm::SetAmpGain(int gain)
{
    if (gain > 1) m_ampgain = gain;
    else {
          std::cout << __FUNCTION__ << " WARNING: Gain value not correct, setting gain to 1!!" << std::endl;
          m_ampgain = 1;
         }
}
// --------------------------------------------------------------------------------------------------------------
void WaveForm::SetWaveId(int id)
{
    m_waveId = id;
}
// --------------------------------------------------------------------------------------------------------------
std::vector<double>* WaveForm::GetVoltage()
{
    return m_voltage;
}
// --------------------------------------------------------------------------------------------------------------
std::vector<double>* WaveForm::GetTime()
{
    return m_time;
}
// --------------------------------------------------------------------------------------------------------------
int WaveForm::GetPolarity()
{
    return m_pol;
}
// --------------------------------------------------------------------------------------------------------------
Long64_t WaveForm::GetSnRate()
{
    if (m_SnRate != -99) return m_SnRate;
    else if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Rate not set!!" << std::endl;
    return -1;
}
// --------------------------------------------------------------------------------------------------------------
double WaveForm::GetTrigg()
{
    return m_trigg;
}
// --------------------------------------------------------------------------------------------------------------
float WaveForm::GetCFDfraction()
{
   return m_fraction;
}
// --------------------------------------------------------------------------------------------------------------
int WaveForm::GetBasePos()
{
    if (m_basepos == undet && m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Base Position not established!!" << std::endl;
    return m_basepos;
}
// --------------------------------------------------------------------------------------------------------------
bool WaveForm::GetIsSignal()
{
    return m_IsSignal;
}
// --------------------------------------------------------------------------------------------------------------
bool WaveForm::GetIsInWidow()
{
    return m_IsInWindow;
}// --------------------------------------------------------------------------------------------------------------
float WaveForm::GetTransimp()
{
    return m_transimp;
}
// --------------------------------------------------------------------------------------------------------------
float WaveForm::GetTransimpFromFile(double signalFFT)
{
    m_TransFile->cd();
	/*
	if (m_debug) std::cout << "Name of the jet Energy Rescale File: " << m_rescalefile.c_str() << std::endl;
	TFile *file = TFile::Open(m_rescalefile.c_str());
	TH2F *h2d = (TH2F*)file->Get("h_dpt_bj");
	TH1F *hcorX = (TH1F*)h2d->ProjectionX();
	int hcor_xbins = hcorX->GetNbinsX() + 1;
	hcorX->Reset();
	for (int i = 1; i<hcor_xbins; i++)
	{
	TH1D* htp = h2d->ProjectionY("_pfy", i, i);
	if (htp->Integral("width")>0)
	{
	double m, dm;
	double mx = htp->GetMean();
	double rms = htp->GetRMS();
	double rmin = mx - 3 * rms;
	double rmax = mx + 3 * rms;
	TF1* mygau = new TF1("mygau", "gaus", rmin, rmax);
	htp->Fit(mygau, "S R Q N");
	m = mygau->GetParameter(1);
	dm = mygau->GetParError(1);
	hcorX->SetBinContent(i, m);
	hcorX->SetBinError(i, dm);
	}
	}
	m_pol_JES = new TF1("mypol", "pol2", 20, 250);
	hcorX->Fit(m_pol_JES, "S R Q N");
	hcorX->SetDirectory(0);
	file->Close();
	*/
    // Code to get the histogram from the file and the transimpedence value
    return m_transimp;
}
// --------------------------------------------------------------------------------------------------------------
float WaveForm::GetAmpGain()
{
    return m_ampgain;
}
// --------------------------------------------------------------------------------------------------------------
std::vector<double> WaveForm::GetAdjVoltage()
{
    if (m_voltageAdj.size() == 0)
       {
        if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Voltage adjustment for pedestal not yet performed!" << std::endl;
        return std::vector<double>();
       }
    else return m_voltageAdj;
}
// --------------------------------------------------------------------------------------------------------------
std::vector<double> WaveForm::GetTimeAdjVolt(float fraction)
{
    if (m_vAdjCFD.size() == 0)
       { 
        if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Time adjustment for CFD not yet performed!" << std::endl;
        return std::vector<double>();
       }
    else if (fraction == -99) return m_vAdjCFD;
    else if (fraction > 0 && fraction < 1) return PulseTimeVoltAdj(m_time, &m_voltageAdj, fraction);
    else {
          if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Ubacceptable CFD Value, using 20%..." << std::endl;
          return m_vAdjCFD;
         }
}
// --------------------------------------------------------------------------------------------------------------
int WaveForm::GetMaxIndx()
{
    if (m_maxIndx != -99) return m_maxIndx;
    else if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << "WARNING: Maximum not defined!!" << std::endl;
    return -1.;
}
// --------------------------------------------------------------------------------------------------------------
double WaveForm::GetMax()
{
    if(m_maxVolt != -99) return m_maxVolt;
    else if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Maximum not defined!!" << std::endl;
    return -1.;
}
// --------------------------------------------------------------------------------------------------------------
double WaveForm::GetMaxErr()
{
    if (m_maxVoltErr != -99) return m_maxVoltErr;
    else if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Maximum time not defined!!" << std::endl;
    return -1.;
}
// --------------------------------------------------------------------------------------------------------------
std::vector<int> WaveForm::GetHalfPoints()
{
    if (m_halfPtsIndx.size() > 0) return m_halfPtsIndx;
    else {
          if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << "WARNING: Half Points not calculated!!" << std::endl;
          return std::vector<int>();
         }
}
// --------------------------------------------------------------------------------------------------------------
double WaveForm::GetMaxTime()
{
    if (m_maxTime != -99) return m_maxTime;
    else if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Maximum time not defined!!" << std::endl;
    return -1.;
}
// --------------------------------------------------------------------------------------------------------------
int WaveForm::GetMinIndx()
{
    if (m_minIndx != - 99) return m_minIndx;
    else if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Minimum not defined!!" << std::endl;
    return -1.;
}
// --------------------------------------------------------------------------------------------------------------
double WaveForm::GetMin()
{
    if (m_minVolt != -99) return m_minVolt;
    else if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Minimum not defined!!" << std::endl;
    return -1.;
}
// --------------------------------------------------------------------------------------------------------------
double WaveForm::GetMinTime()
{
    if (m_minTime != -99) return m_minTime;
    else if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Minimum not defined!!" << std::endl;
    return -1.;
}
// --------------------------------------------------------------------------------------------------------------
int WaveForm::GetStrIndx()
{
    if (m_StrIndx != -99) return m_StrIndx;
    else if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Signal start not defined!!" << std::endl;
    return -1.;
}
// --------------------------------------------------------------------------------------------------------------
int WaveForm::GetEndIndx()
{
    if (m_EndIndx != -99) return m_EndIndx;
    else std::cout << __FUNCTION__ << " WARNING: Signal end not defined!!" << std::endl;
    return -1.;
}
// --------------------------------------------------------------------------------------------------------------
double WaveForm::GetNoise()
{
    if (m_noise != -99) return m_noise;
    else if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Noise not calculated!!" << std::endl;
    return -1.;
}
// --------------------------------------------------------------------------------------------------------------
double WaveForm::GetNoiseErr()
{
    if (m_noiseErr != -99) return m_noiseErr;
    else if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Noise uncertenty not calculated!!" << std::endl;
    return -1.;
}
// --------------------------------------------------------------------------------------------------------------
double WaveForm::GetPedestal()
{
    if (m_pedestal != -99) return m_pedestal;
    else if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Pedestal not determined!!" << std::endl;
    return -1.;
}
// --------------------------------------------------------------------------------------------------------------
double WaveForm::GetPedestalErr()
{
    if (m_pedestalErr != -99) return m_pedestalErr;
    else if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Pedestal uncertenty not determined!!" << std::endl;
    return -1.;
}
// --------------------------------------------------------------------------------------------------------------
double WaveForm::GetCharge()
{
    if (m_charge != -99) return m_charge;
    else if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Charge has not been calculated!!" << std::endl;
    return -1.;
}
// --------------------------------------------------------------------------------------------------------------
double WaveForm::GetRiseTime(float bottom, float top)
{
    if (m_RiseTime != -99)
       {
        if (bottom == -99 || top == -99) return m_RiseTime;
        else if (bottom > 0 && bottom < 1 && top > 0 && top < 1) return RiseTimeLinear(&m_voltageAdj, m_SnRate, top, bottom);
        else {
              std::cout << __FUNCTION__ << " WARNING: Unacceptable percentage values, using 20 - 80%..." << std::endl;
              return m_RiseTime;
             }
       }
    else if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Rise time has not been estimated!!" << std::endl;
    return -1.;
}
// --------------------------------------------------------------------------------------------------------------
double WaveForm::GetCFDTime(float fraction)
{
    if (m_CFDTime != -99) 
       {
        if (fraction == -99) return m_CFDTime;
        else if (fraction > 0 && fraction < 1 ) return CFDTimeLinear(&m_voltageAdj, m_time, fraction);
        else {
              std::cout << __FUNCTION__ << " WARNING: Unacceptable CFD value, using 20%..." << std::endl;
              return m_CFDTime;
             }
       }
    else if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: CFD time has not been estimated!!" << std::endl;
    return -1.;
}
// --------------------------------------------------------------------------------------------------------------
double WaveForm::GetTriggTime(float trigg)
{
    if (m_TriggTime != -99)
       {
        if (trigg == -99) return m_TriggTime;
        else if (trigg <= 0) return FirstTimeForVoltage(&m_voltageAdj, m_time, trigg);
        else {
              std::cout << __FUNCTION__ << " WARNING: Unacceptable values for trigger, using standard 50 mV" << std::endl;
              return m_TriggTime;
            }
       }
    else if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Trigger time has not been estimated!!" << std::endl;
    return -1;
}
// --------------------------------------------------------------------------------------------------------------
double WaveForm::GetdVdTMax()
{
    if (m_dVdTMax != -99) return m_dVdTMax;
    else if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Maximum dV/dT has not been computed!!" << std::endl;
    return -1.;
}
// --------------------------------------------------------------------------------------------------------------
double WaveForm::GetdVdTCFD(float fraction, int ndif)
{
    if (m_dVdTCFD != -99) 
       {
        if (fraction == -99 && ndif == 0) return m_dVdTCFD;
        else if (fraction < 1 && fraction < 0 && ndif >= 0) return dVdTCFDLinear(&m_voltageAdj, m_SnRate, fraction, ndif);
        else {
              std::cout << __FUNCTION__ << " WARNING: Unacceptable values for CFD and number of integration points"
                                        << ", using standard 20% and single point integration" << std::endl;
              return m_dVdTCFD;
            }
       }
    else if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: dV/dT at CFD value has not been computed!!" << std::endl;
    return -1.;
}
// --------------------------------------------------------------------------------------------------------------
double WaveForm::GetCFDToT(float fraction)
{
    if (m_CFDToT != -99)
       {
        if (fraction == -99) return m_CFDToT;
        else if (fraction < 1 && fraction < 0) return CFDToTLinear(&m_voltageAdj, m_time, fraction);
        else {
              std::cout << __FUNCTION__ << " WARNING: Unacceptable value for CFD, using standard 20%" << std::endl;
              return m_CFDToT;
            }
       }
    else if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: ToT for CFD value has not been computed!!" << std::endl;
    return -1.;
}
// --------------------------------------------------------------------------------------------------------------
double WaveForm::GetTriggToT(float trigg)
{
    if (m_TriggToT != -99)
       {
        if (trigg == -99) return m_TriggToT;
        else if (trigg <= 0) return TriggToTLinear(&m_voltageAdj, m_time, trigg);
        else {
              std::cout << __FUNCTION__ << " WARNING: Unacceptable trigger value, using standard 50mV" << std::endl;
              return m_TriggToT;
            }
       }
    else if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: ToT for trigger value has not been computed!!" << std::endl;
    return -1.;
}
// --------------------------------------------------------------------------------------------------------------
double WaveForm::GetSignalFFT()
{
    if (m_signalFFT != -99) return m_signalFFT;
    else if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Fast Fourier Transform for signal not performed!!" << std::endl;
    return -1.;
}
// --------------------------------------------------------------------------------------------------------------
double WaveForm::GetNoiseFFT()
{
    if (m_noiseFFT != -99) return m_noiseFFT;
    else if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Fast Fourier Transform for noise not performed!!" << std::endl;
    return -1.;
}
// --------------------------------------------------------------------------------------------------------------
double WaveForm::GetJitterNdVdT()
{
    if (m_jitter1 != -99) return m_jitter1;
    else if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Jitter [noise/dVdt] calculation not performed!!" << std::endl;
    return -1.;
}
// --------------------------------------------------------------------------------------------------------------
double WaveForm::GetJitterRiseSNR()
{
    if (m_jitter2 != -99) return m_jitter2;
    else if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Jitter [Rise Time/SNR] calculation not performed!!" << std::endl;
    return -1.;
}
// --------------------------------------------------------------------------------------------------------------
double WaveForm::GetFrequency(int start, int stop)
{
    if (m_signalFFT != -99)
       {
        if (start > 0 && stop > 0 && start < stop) return m_WvBase->LGADBase::FFT(&m_voltageAdj, m_SnRate, start, stop);
        else {
              std::cout << __FUNCTION__ << " WARNING: Unacceptable values for start and stop, abording..." << std::endl;
              return -1.;
            }
       }
    else if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: FFT not yet performed!!" << std::endl;
    return -1.;
}
// --------------------------------------------------------------------------------------------------------------
bool WaveForm::GetSaturation()
{
    if (m_maxIndx != -99 && m_StrIndx != -99 && m_EndIndx != -99) return m_VoltStr;
    else if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Saturation test not performed!!" << std::endl;
    return false;
}
// --------------------------------------------------------------------------------------------------------------
WaveForm::basepos WaveForm::FindPosition(std::vector<double> *w)
{
    basepos position;
    bool zerocr = false;
    for (int v = 0; v < (int)w->size(); v++)
        {
        if ((w->at(v) < 0 && w->at(2) > 0) || (w->at(v) > 0 && w->at(2) < 0))
           {
            zerocr = true;
            break;
           }
        }
    if (zerocr == false && w->at(2) > 0) position = poside;
    else if (zerocr == false && w->at(2) < 0) position = negside;
    else if (zerocr == true) position = zerocross;
    else position = undet;
    return position;
}
// --------------------------------------------------------------------------------------------------------------
WaveForm::polarity WaveForm::FindPolarity(std::vector<double> *w)
{
    polarity pol = undef;
    std::vector<double> firstder = m_WvBase->LGADBase::Derivate(w, 3);

    if (firstder.size() != 0)
       { 
        std::vector<double> secordr = m_WvBase->LGADBase::Derivate(&firstder, 1);
        if (secordr.size() != 0)
           { 
            std::vector<int> maxindx;
            std::vector<double> maxvolt;
            std::vector<int> minindx;
            std::vector<double> minvolt;
            for (unsigned int h = 0; h < secordr.size(); h++)
                { 
                 // We find two consecutive points with different sign, examine the sign of second derivative
                 // No linear interpolation needed, it is already incorporated in derivation
                 if ((firstder.at(h) >= 0 && firstder.at(h+1) < 0) || (firstder.at(h) < 0 && firstder.at(h+1) >= 0) ||
                     (firstder.at(h) > 0 && firstder.at(h+1) <= 0) || (firstder.at(h) <= 0 && firstder.at(h + 1) > 0))
                    {
                     if (secordr.at(h) >= 0) { minvolt.push_back(w->at(h+3)); minindx.push_back(h+3); }
                     else if (secordr.at(h) < 0) { maxvolt.push_back(w->at(h+3)); maxindx.push_back(h+3); }
                    }
                }
            if (maxvolt.size() !=0 && minvolt.size()!=0)
               {
                m_maxavg = m_WvBase->LGADBase::Mean(&maxvolt);
                m_minavg = m_WvBase->LGADBase::Mean(&minvolt);
                double maxstdev = m_WvBase->LGADBase::Stdev(&maxvolt, -1, -1, m_maxavg);
                double minstdev = m_WvBase->LGADBase::Stdev(&minvolt, -1, -1, m_minavg);
                // caculate extremum distances
                double minextr = fabs(*std::min_element(minvolt.begin(), minvolt.end()) - m_minavg);
                double maxextr = fabs(*std::max_element(maxvolt.begin(), maxvolt.end()) - maxstdev);
                // Combine extremum std, extremum limit values and base position to determine polarity
                if (maxstdev >= minstdev && maxextr > minextr) pol = pos;
                else if (maxstdev <= minstdev && maxextr < minextr) pol = neg;
                else if (maxextr > minextr && minstdev < (1.05*maxstdev)) pol = pos;
                else if (maxextr < minextr && maxstdev < (1.05*minstdev)) pol = neg;
                else if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Cannot determine polarity with available info!!" << std::endl;
               }
            else if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " ERROR: Cannot determine polarity, no maxima or minima found!" << std::endl;
           }
        else if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " ERROR: Cannot determine polarity, second derivative size too small!" << std::endl;
       }
    else if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " ERROR: Cannot determine polarity, first derivative size too small!" << std::endl;
    return pol;
}
// --------------------------------------------------------------------------------------------------------------
// Function to find the RELATIVE to polarity maxima of the pulse and return the inddex of that value
int WaveForm::VoltMaxIndx(std::vector<double> *w, bool adj, bool poldef)
{
    int ind_max = -1;
    int npoints = w->size();
    double pmax;
    bool posedge = false;
    if (poldef && m_pol == undef) m_pol = FindPolarity(w);
    if (adj == true || m_pol == pos || m_pol == undef) posedge = true;
    if (posedge) pmax = -990;
    else pmax = 990;
    // Works for positive and negative polarity pulses in any baseline position
    for (int j = 0; j < npoints; j++)
        { 
         if ((w->at(j) > pmax && posedge) ||
             (w->at(j) < pmax && !posedge))
            {
             pmax = w->at(j); // update pulse max
             ind_max = j;
            }
        }
    return ind_max;
}
// --------------------------------------------------------------------------------------------------------------
// Returna the voltage value of the maxima
double WaveForm::VoltMax(std::vector<double> *w, bool adj, bool poldef)
{
    if (m_maxIndx == -99) m_maxIndx = VoltMaxIndx(w, adj, false);
    if (m_maxIndx != -1) return w->at(m_maxIndx);
    else return 0;
}
// --------------------------------------------------------------------------------------------------------------
// Returna the uncertenty of the maximum voltage
double WaveForm::VoltMaxErr(std::vector<double> *w, bool adj)
{
    if (m_maxVolt == -99) m_maxVolt = VoltMax(w, adj);
    double VoltUp = 0.0;
    double VoltDo = 0.0;
    if ((w->size() - m_maxIndx) > 1) VoltUp = LinearInter(m_maxVolt, m_maxIndx, w->at(m_maxIndx + 1), m_maxIndx + 1, m_maxIndx + 0.5);
    else {
          VoltDo = LinearInter(w->at(m_maxIndx - 1), m_maxIndx - 1, m_maxVolt, m_maxIndx, m_maxIndx - 0.5);
          VoltUp = VoltDo;
          return 0.5*(fabs(m_maxVolt - VoltDo) + fabs(VoltUp - m_maxVolt));
         }
    if (m_maxIndx > 0) VoltDo = LinearInter(w->at(m_maxIndx - 1), m_maxIndx - 1, m_maxVolt, m_maxIndx, m_maxIndx - 0.5);
    else {
          VoltUp = LinearInter(m_maxVolt, m_maxIndx, w->at(m_maxIndx + 1), m_maxIndx + 1, m_maxIndx + 0.5);
          VoltDo = VoltUp;
          return 0.5*(fabs(m_maxVolt - VoltDo) + fabs(VoltUp - m_maxVolt));
         }
    return 0.5*(fabs(m_maxVolt - VoltDo) + fabs(VoltUp - m_maxVolt));
}
// --------------------------------------------------------------------------------------------------------------
// Returns the time when the pulse mnaxima occured
double WaveForm::TimeMax(std::vector<double> *w, Long64_t snrate, bool adj, bool poldef)
{
    if (m_maxIndx == -99) m_maxIndx = VoltMaxIndx(w, adj, false);
    return m_maxIndx*(1 / (double)snrate);
}
// --------------------------------------------------------------------------------------------------------------
// Function to find the RELATIVE to polarity minima of the pulse. It uses the output of the maximum functyion to define the signal polarity
int WaveForm::VoltMinIndx(std::vector<double> *w, bool adj, bool poldef)
{
    int ind_min = -1;
    int npoints = w->size();
    if (m_maxVolt == -99) m_maxVolt = VoltMax(w, adj, poldef);
    double pmin = m_maxVolt;
    bool posedge = false;
    if (m_pol == pos || adj == true || m_pol == undef) posedge = true;

    for (int j = 0; j < npoints; j++)
        {
         if ((posedge && w->at(j) < pmin) || (!posedge && w->at(j) > pmin))
            {
             pmin = w->at(j); // update pulse min
             ind_min = j;
            }
        }
    return ind_min;
}
// --------------------------------------------------------------------------------------------------------------
// Returns the voltge of the minimum vlaue
double WaveForm::VoltMin(std::vector<double> *w, bool adj, bool poldef)
{
    if (m_minIndx == -99) m_minIndx = VoltMinIndx(w, adj, poldef);
    if (m_minIndx !=-1) return w->at(m_minIndx);
    else return 0;
}
// --------------------------------------------------------------------------------------------------------------
// Returns the time when the pulse minima is observed
double WaveForm::TimeMin(std::vector<double> *w, Long64_t snrate, bool adj, bool poldef)
{
    if (m_minIndx == -99) m_minIndx = VoltMinIndx(w, adj, poldef);
    return m_minIndx*(1/ (double)snrate);
}
// --------------------------------------------------------------------------------------------------------------
// Find the points where the voltage value is 90%, 80%, 70, or 60% of the difference between max and min values consecutively
std::vector<int> WaveForm::FindHalfPoints(std::vector<double> *w, bool adj)
{
    if (m_maxVolt == -99) m_maxVolt = VoltMax(w, adj);
    if (m_minVolt == -99) m_minVolt = VoltMin(w, adj);
    double halfvalue = 0.0;
    std::vector<int> halfPtsIndx;

    // Calculate 90% of the distance between max and min voiltages
    if (m_pol == pos || m_pol == undef || adj) halfvalue = m_minVolt + abs(m_maxVolt - m_minVolt)*0.9;
    else halfvalue = m_minVolt - abs(m_maxVolt - m_minVolt)*0.8;
    for (unsigned int d = 0; d < (w->size()-1); d++)
        {
         if ((w->at(d) < halfvalue && w->at(d + 1) > halfvalue) ||
             (w->at(d) > halfvalue && w->at(d + 1) < halfvalue) || 
             (w->at(d) == halfvalue && w->at(d+1) != halfvalue)) halfPtsIndx.push_back(d);
        }

    if (halfPtsIndx.size() > 2 && (halfPtsIndx.at(halfPtsIndx.size()-1)-halfPtsIndx.at(0)) >= 8)
       {
        halfPtsIndx.clear();
        if (m_pol == pos || m_pol == undef) halfvalue = m_minVolt + abs(m_maxVolt - m_minVolt)*0.8;
        else halfvalue = m_minVolt - abs(m_maxVolt - m_minVolt)*0.7;
        for (unsigned int d = 0; d < (w->size()-1); d++)
            {
             if ((w->at(d) < halfvalue && w->at(d + 1) > halfvalue) ||
                (w->at(d) > halfvalue && w->at(d + 1) < halfvalue) || 
                (w->at(d) == halfvalue && w->at(d+1) != halfvalue)) halfPtsIndx.push_back(d);
            }
        if (halfPtsIndx.size() > 2 && (halfPtsIndx.at(halfPtsIndx.size()-1)-halfPtsIndx.at(0)) >= 12)
           {
            halfPtsIndx.clear();
            if (m_pol == pos || m_pol == undef) halfvalue = m_minVolt + abs(m_maxVolt - m_minVolt)*0.7;
            else halfvalue = m_minVolt - abs(m_maxVolt - m_minVolt)*0.6;
            for (unsigned int d = 0; d < (w->size()-1); d++)
                {
                 if ((w->at(d) < halfvalue && w->at(d + 1) > halfvalue) ||
                    (w->at(d) > halfvalue && w->at(d + 1) < halfvalue) || 
                    (w->at(d) == halfvalue && w->at(d + 1) != halfvalue)) halfPtsIndx.push_back(d);
                }
            if (halfPtsIndx.size() > 2 && (halfPtsIndx.at(halfPtsIndx.size()-1)-halfPtsIndx.at(0)) >= 16)
               { 
                halfPtsIndx.clear();
                if (m_pol == pos || m_pol == undef) halfvalue = m_minVolt + abs(m_maxVolt - m_minVolt)*0.6;
                else halfvalue = m_minVolt - abs(m_maxVolt - m_minVolt)*0.5;
                for (unsigned int d = 0; d < (w->size()-1); d++)
                    {
                     if ((w->at(d) < halfvalue && w->at(d + 1) > halfvalue) ||
                         (w->at(d) > halfvalue && w->at(d + 1) < halfvalue) || 
                         (w->at(d) == halfvalue && w->at(d + 1) != halfvalue)) halfPtsIndx.push_back(d);
                    }
                if (halfPtsIndx.size() > 2 && (halfPtsIndx.at(halfPtsIndx.size()-1)-halfPtsIndx.at(0)) >= 20) halfPtsIndx.push_back(0);
                else halfPtsIndx.push_back(1);              
               }
            else halfPtsIndx.push_back(1);
           }
        else halfPtsIndx.push_back(1);
       }
      else halfPtsIndx.push_back(1);

    return halfPtsIndx;
}
// --------------------------------------------------------------------------------------------------------------
bool WaveForm::IsSignal(std::vector<int>* halfPtsIndx)
{
    bool signal = false;
    unsigned int hpoint = halfPtsIndx->size();
    if (hpoint > 2 && halfPtsIndx->at(hpoint -1) == 0) signal = false;
    else signal = true;
    return signal; 
}
// --------------------------------------------------------------------------------------------------------------
// Signal saturation test
bool WaveForm::VoltSatur(std::vector<double> *w, bool poldef)
{
    bool satur = false;
    if (m_maxIndx == -99) m_maxIndx = VoltMaxIndx(w, poldef);
    if (m_StrIndx == -99) m_StrIndx = StartIndx(w, false, false);
    if (m_EndIndx == -99) m_EndIndx = EndIndx(w, false);

    // Find number of points with sanve value as max
    std::vector<unsigned int> indexvect;
    if (m_maxIndx >= 0)
       {
        for (int l = 0; l < (int)(w->size()); l++) 
            {
             if (w->at(l) == w->at(m_maxIndx) && l >= m_StrIndx && l <= m_EndIndx) indexvect.push_back(l);
            }
       }

   // Allow for exactly two consecutive max voltages
   if (indexvect.size() == 2)
      {
       unsigned int v = 0;
       for (unsigned int k = 0; k < (indexvect.size()-1); k++) if (indexvect.at(k + 1) - indexvect.at(k) == 1) v++;
       if (v == 1) satur = false;
       else satur = true;
      }
   else if (indexvect.size() > 2) satur = true;

    return satur;
}
// --------------------------------------------------------------------------------------------------------------
// Function to find if pulse is within window. 
bool WaveForm::InWindow(std::vector<double>* w)
{
    bool window = false;
    // Populate the half points indices vector if not populated already
    if (m_halfPtsIndx.size() == 0) m_halfPtsIndx = FindHalfPoints(w, false);
    if (m_maxIndx == -99) m_maxIndx = VoltMaxIndx(w, false);
    // Points after the end of the signal if any
    if (m_EndIndx == -99) m_EndIndx = EndIndx(w, false);
    int endindx = -1;
    if (m_EndIndx != -1) endindx = m_voltage->size() - endindx;
    // points before the begining of the signal if any
    if (m_StrIndx == -99) m_StrIndx = StartIndx(w, false);

    if (m_halfPtsIndx.size() >= 2)
       {
        if (endindx == -1 && (m_StrIndx > (int)w->size()/2))
           {
            if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Pulse within window but located at the higher extremity!" << std::endl;
            window = true;
           }
        if (m_StrIndx == -1 && (endindx > (int)w->size()/2))
           {
            if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Pulse outside window, located at the lower extremity, rise time calculation will fail!" << std::endl;
            window = false;
           }
        else window = true;
       }
    else if (m_halfPtsIndx.size() == 1)
            {
             if (endindx == -1 && (m_StrIndx > (int)w->size() / 2))
                {
                 if (m_maxIndx == (int)w->size()) 
                    { 
                      window = false;
                      if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Pulse rising edge ouside window at high side!" << std::endl;
                    }
                 else {
                       window = true;
                       if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Pulse rising edge within window but located at the higher extremity!" << std::endl;
                      }
                }
             else if (m_StrIndx == -1 && (endindx > (int)w->size() / 2))
                     {
                      if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Pulse outside window, located at the lower extremity, rise time calculation will fail!" << std::endl;
                      window = false;
                     }
             else window = false;
            }
    else  window = false;
    return window;
}
// --------------------------------------------------------------------------------------------------------------
// Function to determin the start of the signal. If called without the baseline value, it gives a rougph estimate of 
// the zero crossing value with respect to the maximum. Configured for both positive and negative signals.
int WaveForm::StartIndx(std::vector<double> *w, bool adj, bool poldef)
{
    int istart = -1;
    double basevolt = 0.0;
    if (poldef && m_pol == undef) m_pol = FindPolarity(w);
    bool posedge = true;
    if (!adj)
       {
        if (m_maxIndx == -99) m_maxIndx = VoltMaxIndx(w, adj);
        if (m_pol == neg) { basevolt = m_minavg; posedge = false;}
        else basevolt = m_maxavg; 
       }
    else {
          posedge = true;
          if (m_noise != -99) basevolt += fabs(m_noise);
          if (m_pedestalErr != -99) basevolt += m_pedestalErr;
         }

    // find index of start of pulse
    for (int j = m_maxIndx; j > 0; j--)
        {
         if ((posedge && (w->at(j) > basevolt && w->at(j - 1) < basevolt)) ||
             (!posedge && (w->at(j) < basevolt && w->at(j - 1) > basevolt)))
            {
             istart = j;
             break;
            }
         else if (w->at(j) == basevolt) { istart = j; break; }
        }
    // Sfeguard agains signals at the begining of DAQ window
    if (adj && istart == -1) istart = 0;

    if (istart == -1 && m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Signal start index not found!! " << std::endl;
    return istart;
}
// --------------------------------------------------------------------------------------------------------------
// Without adjustment, find the point that is 20% higher than the lowest detectable voltage of the pulse
int WaveForm::EndIndx(std::vector<double> *w, bool adj, bool poldef)
{
    int iend = -1;
    double basevolt = 0.0;
    if (poldef && m_pol == undef) m_pol = FindPolarity(w);
    bool posedge = false;
    if (!adj)
       {
        if (m_maxIndx == -99) m_maxIndx = VoltMaxIndx(w, adj);
        if (m_pol == neg) { basevolt = m_minavg; posedge = false; }
        else basevolt = m_maxavg; 
       }
    else {
          posedge = true;
          if (m_noise != -99) basevolt += fabs(m_noise);
          if (m_pedestalErr != -99) basevolt += m_pedestalErr;
         }

     // find index of end of pulse
    for (int j = m_maxIndx; j < (int)w->size() - 1; j++)
        {
         if ((posedge && (w->at(j) >= basevolt && w->at(j + 1) <= basevolt)) ||
             (!posedge && (w->at(j) <= basevolt && w->at(j + 1) >= basevolt)))
            {
             iend = j;
             break;
            }
        }

    if (iend == -1 && m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Signal end index not found!! " << std::endl;
    return iend;
}
// --------------------------------------------------------------------------------------------------------------
// Functon to determine the indices of the points used for noise calculation. Defined as the 80% of the point preceeding the
// begining of the signal. It ignors the first and last 10% of points before and after the signal respectively
std::pair <unsigned int, unsigned int> WaveForm::FindNoisePoints(std::vector<double> *w)
{
    int npoints = w->size();
    if (m_StrIndx == -99) m_StrIndx = StartIndx(w, false, false);
    std::pair <unsigned int, unsigned int> noisepoints (0, 0);
    if (ceil(0.2*m_StrIndx) > 1) noisepoints.first = ceil(0.2*m_StrIndx);
    else noisepoints.first = 1;
    // Enforce a 4 points noise calculation
    if (ceil(0.8*m_StrIndx) > ceil(0.2*m_StrIndx) && ceil(0.8*m_StrIndx) >= 10) noisepoints.second = ceil(0.8*m_StrIndx);
    else {
          if (m_EndIndx == -99) m_EndIndx = EndIndx(w, false, false);
          if (m_EndIndx <= ((int)(w->size()) - 5)) 
             {
              noisepoints.second = m_EndIndx + ceil(0.8*(w->size() - 3 - m_EndIndx));
              noisepoints.first = m_EndIndx + ceil(0.2*(w->size() - 3 - m_EndIndx));
             }
          else if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Noise Determination Failed, starting point"
                                                                     << noisepoints.first << ", ending point: " << ceil(0.9*m_StrIndx) << "!!!!" << std::endl;
         }
    return noisepoints;
}
// --------------------------------------------------------------------------------------------------------------
std::vector<double> WaveForm::PulseAdj(std::vector<double> *w, double baseline, int factor)
{
    int npoints = w->size();
    std::vector<double> wadj;
    wadj.reserve(npoints);
    for (int i = 0; i < npoints; i++) wadj.push_back(factor*(w->at(i)-baseline));
    m_maxIndx = VoltMaxIndx(&wadj, true);
    m_maxVolt = VoltMax(&wadj, true);
    m_maxVoltErr = VoltMaxErr(&wadj, true);
    m_maxTime = TimeMax(&wadj, m_SnRate, true);
    m_minIndx = VoltMinIndx(&wadj, true);
    m_minVolt = VoltMin(&wadj, true);
    m_minTime = TimeMin(&wadj, m_SnRate, true);
    m_halfPtsIndx = FindHalfPoints(&wadj, true);
    m_StrIndx = StartIndx(&wadj, true);
    m_EndIndx = EndIndx(&wadj, true);
    if (m_WvBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: After pedestal subbtraction and inversion : baseline " << baseline 
                                                           << ", corrected vector size: " << wadj .size() << ", factor " << factor << std::endl;

    return wadj;
}
// --------------------------------------------------------------------------------------------------------------
double WaveForm::CollectedCharge(std::vector<double> *w, Long64_t snrate, float transimp, float ampgain, int start, int stop)
{
    double charge, time_difference;
    charge = 0.;
    time_difference = 1 / (double)snrate;
    if (stop < 0 || stop > (int)w->size() || stop < start) stop = (int)w->size();
    if (start < 0 || start > (int)w->size() || start > stop) start = 0;
    if (stop < (int)(w->size())) for (int j = start; j < stop + 1; j++) { charge += abs(w->at(j)); }
    else  for (int j = start; j < stop; j++) { charge += abs(w->at(j)); }
    charge = charge* time_difference / (transimp * (float)ampgain);
    if (m_WvBase->LGADBase::GetVerbosity() >= 2)  std::cout << __FUNCTION__ << " INFO: Charge Calculation-> charge " << charge << ", gain: " << ampgain 
                                                            << ", transimpedence: "<< transimp  << ", timebin size: " << time_difference << std::endl;
    return charge; 
}
// --------------------------------------------------------------------------------------------------------------
double WaveForm::FirstTimeForVoltage(std::vector<double> *w, std::vector<double> *t, double volt)
{
    double tar = -1.;
    int ind_tar = -1;
    if (m_StrIndx == -99) m_StrIndx = StartIndx(w, true);
    if (m_maxIndx == -99) m_maxIndx = VoltMaxIndx(w, true);

    int dx = 0;
    if (m_StrIndx > 0) dx = m_StrIndx - 1;
    else dx = m_StrIndx;
    for (int j = m_maxIndx; j >= dx; j--)
        {
         if (w->at(j) <= volt)
            {
             ind_tar = j;
             break;
            }
        } // index of first trigger

    if (ind_tar != -1 && ind_tar != m_maxIndx)
       {
        if (w->at(ind_tar) == volt) tar = t->at(ind_tar); // if the lower value is exacly at one of the recorded points
        else  tar = LinearInter(t->at(ind_tar), w->at(ind_tar), t->at(ind_tar + 1), w->at(ind_tar + 1), volt);
        if (tar < t->back() && tar > t->front()) return tar;
        else {
              if (m_WvBase->LGADBase::GetVerbosity() > 0 ) std::cout << __FUNCTION__ << " ERROR: Failed to find voltage time for " << volt << " trigger!" << std::endl;
              return -1.;
             }
       }
    else if (ind_tar == m_maxIndx)
            {
             if (m_WvBase->LGADBase::GetVerbosity() > 0 ) std::cout << __FUNCTION__ << " ERROR: Voltage value " << volt << " is higher than pulse maximum!" << std::endl;
             return -1.;
            }
    else {
          if (m_WvBase->LGADBase::GetVerbosity() > 0 ) std::cout << __FUNCTION__ << " ERROR: Failed to find position where voltage " << volt << " first occurrs!" << std::endl;
          return -1.;
         }
}
// --------------------------------------------------------------------------------------------------------------
double WaveForm::CFDTimeLinear(std::vector<double> *w, std::vector<double> *t, float fraction)
{
    if (m_maxVolt == -99) m_maxVolt = VoltMax(w, true);
    if (m_minVolt == -99) m_minVolt = VoltMin(w, true);
    if (fraction*m_maxVolt > m_minVolt) return FirstTimeForVoltage(w, t, fraction*m_maxVolt);
    else return -1.;
}
// --------------------------------------------------------------------------------------------------------------
std::vector<double> WaveForm::PulseTimeVoltAdj(std::vector<double> *t, std::vector<double> *w, float fraction)
{
    unsigned int npoints = t->size();
    std::vector<double> vadjCFD;
    vadjCFD.resize(npoints, -1.);
    double time = 0.0;
    if (m_CFDTime == -99 || fraction != m_fraction) time = CFDTimeLinear(w, t, fraction);
    else time = m_CFDTime;

    // We need to recalculate the voltage vector wwith linear interpolation for the relative time shift
    // This is to ensure adding with different pulses in a bined histogram and avoid creating gost effects
    // Counter intuitive but needs to be applied.
    // Recalculate time vector with 0 set at CFD value
    int indxCFD = -99;
    for (unsigned int a = 0; a < (npoints-1); a++) if (time >= t->at(a) && time < t->at(a + 1)) indxCFD = a;
    if (indxCFD != -99)
       {
        double shift = fabs(time - t->at(indxCFD));
        for (unsigned int g = 0; g < npoints; g++)
            {
             if (g == 0) vadjCFD.at(g) = 0.0; // maintain number of points
             else vadjCFD.at(g) = LinearInter(w->at(g-1), t->at(g-1), w->at(g), t->at(g), t->at(g-1) + shift);
            }
        if (m_WvBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: After timewalk correction : CFD fraction " 
                                                               << fraction << ", corrected vector size: " << vadjCFD.size() << std::endl;
       }
    else if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " ERROR: Failed to find time index of CFD time: " << time << std::endl;

    return vadjCFD;
}
// --------------------------------------------------------------------------------------------------------------
double WaveForm::dVdTMaxLinear(std::vector<double> *w, Long64_t snrate)
{
    // function to calculate max dv/dt of the pulse 
    double dvdt = 0.; 
    if (m_maxIndx == -99) m_maxIndx = VoltMaxIndx(w, true);
    if (m_StrIndx == -99) m_StrIndx = StartIndx(w, true);

    int dx = 0;
    if (m_StrIndx > 0) dx = m_StrIndx;
    else dx = m_StrIndx + 1;
    for (int j = m_maxIndx; j >= dx; j--)
        {
         if (fabs(w->at(j) - w->at(j-1)) > fabs(dvdt)) dvdt = w->at(j) - w->at(j - 1);
        } 
    if (dvdt != 0) return dvdt*snrate;
    else {
          if (m_WvBase->LGADBase::GetVerbosity() > 0 ) std::cout << __FUNCTION__ << " ERROR: Failed to calculate maximum dV/dt value!" << std::endl;
          return -1.;
         }
}
// --------------------------------------------------------------------------------------------------------------
double WaveForm::dVdTCFDLinear(std::vector<double> *w, Long64_t snrate, float fraction, int ndif)
{
    // Function to calculate dv/dt at a given constant fraction value
    if (m_maxIndx == -99) m_maxIndx = VoltMaxIndx(w, true);
    if (m_StrIndx == -99) m_StrIndx = StartIndx(w, true);
    if (m_maxVolt == -99) m_maxVolt = VoltMax(w, true);

    int ifraction = -99;
    double dvdt = -0.;

    // find index of first point before constant fraction of pulse
    int dx = 0;
    if (m_StrIndx > 0) dx = m_StrIndx-1;
    else dx = m_StrIndx;
    for (int j = m_maxIndx; j >= dx; j--)
        {
          if (w->at(j) <= m_maxVolt*fraction)
             { 
              ifraction = j; 
              break; 
             }
        }
    if (ifraction != -99)
       {
        if (ndif == 0) { dvdt = (w->at(ifraction + 1) - w->at(ifraction))*snrate; }
        else { dvdt = (w->at(ifraction + ndif - 1) - w->at(ifraction + ndif - 1))*snrate / (2 * ndif); }
        if (dvdt != 0) return dvdt;
        else {
              if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " ERROR: Failed to calculate dV/dT at a CFD of " << fraction << std::endl;
              return -1.;
             }
       }
    else {
          if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " ERROR: Failed to determine first point at a CFD of " << fraction << std::endl;
          return -1.;
         }
}
// --------------------------------------------------------------------------------------------------------------
double WaveForm::RiseTimeLinear(std::vector<double> *w, Long64_t snrate, float top, float bottom)
{
    double tbottom = -1., ttop = -1.;
    bool top_found = false, bottom_found = false;
    int top_indx = -1, bottom_indx = -1;
    if (m_maxVolt == -99) m_maxVolt = VoltMax(w, true);
    if (m_StrIndx == -99) m_StrIndx = StartIndx(w, true);

    if (m_maxIndx > m_StrIndx && m_StrIndx >= 0 && m_maxIndx >= 0)
       {
        int ndx = 0;
        if (m_StrIndx > 0) ndx = (m_StrIndx - 1);
        else ndx = m_StrIndx;
        for (int j = m_maxIndx; j >= ndx; j--)
            { 
             if (w->at(j) >= m_maxVolt*top)
                {
                 top_indx = j;
                 top_found = true;
                }
             if ((w->at(j) <= m_maxVolt*bottom) && !bottom_found)
                {
                 bottom_found = true;
                 bottom_indx = j;
                }
             if (bottom_found && top_found) break;
            }
        // Safeguard for very small pulses where top and bottom indexes coincide with start and max indexes
        if (top_found && !bottom_found) bottom_indx = m_StrIndx;
       }

    if (bottom_indx != -1 && top_indx != -1 && top_indx > bottom_indx)
       { 
        if (fabs(w->at(bottom_indx)) == fabs(m_maxVolt*bottom)) tbottom = (1 / (double)snrate)*bottom_indx; // if the lower value is exacly at one of the recorded points
        else tbottom = LinearInter((1/(double)snrate)*bottom_indx, w->at(bottom_indx), (1/(double)snrate)*(bottom_indx+1), w->at(bottom_indx+1), m_maxVolt*bottom);
        if (fabs(w->at(top_indx)) == fabs(m_maxVolt*top)) ttop = (1 / (double)snrate)*top_indx; // if the top value is exacly at one of the recorded points
        else ttop = LinearInter((1 / (double)snrate)*(top_indx - 1), w->at(top_indx-1), (1 / (double)snrate)*top_indx, w->at(top_indx), m_maxVolt*top);
        if (fabs(ttop - tbottom) > 0) return fabs(ttop - tbottom); // rise time
        else {
              if (m_WvBase->LGADBase::GetVerbosity() > 0 ) std::cout << __FUNCTION__ << " ERROR: Linear interpolation failed!" << std::endl;
              return -1.;
             }
       }
    else {
          if (m_WvBase->LGADBase::GetVerbosity() > 0 ) std::cout << __FUNCTION__ << " ERROR: Failed to determine upper and lower values" << std::endl;
          return -1.;
         }
}
// --------------------------------------------------------------------------------------------------------------
double WaveForm::SecondTimeForVoltage(std::vector<double> *w, std::vector<double> *t, double volt)
{
    if (m_EndIndx == -99) m_EndIndx = EndIndx(w, true);
    if (m_maxIndx == -99) m_maxIndx = VoltMaxIndx(w, true);
    double tar = -1.;
    int ind_t = -1;
    for (int j = m_maxIndx; j <= m_EndIndx+1; j++)
        { 
         if (w->at(j) <= volt)
            {  
             ind_t = j;
             break;
            } 
        } // index of second voltage
    if (ind_t != -1 && ind_t != m_maxIndx)
       { 
        if (w->at(ind_t) == volt) tar = t->at(ind_t); // if the top value is exacly at one of the recorded points
        else tar = LinearInter(t->at(ind_t - 1), w->at(ind_t - 1), t->at(ind_t), w->at(ind_t), volt);
        if (tar != -1 && tar < t->back() && tar > t->front()) return tar; // Time over threshold
        else {
              if (m_WvBase->LGADBase::GetVerbosity() > 0 ) std::cout << __FUNCTION__ << " ERROR: Failed to calculate second time for voltage: " << volt << std::endl;
              return -1.;
             }
       }
    else if (ind_t == m_maxIndx)
            {
             if (m_WvBase->LGADBase::GetVerbosity() > 0 ) std::cout << __FUNCTION__ << " ERROR: Voltage value " << volt << " is larger than pulse maximum!" << std::endl;
             return -1.;
            }
    else {
          if (m_WvBase->LGADBase::GetVerbosity() > 0 ) std::cout << __FUNCTION__ << " ERROR: Failed to determine second index for voltage " << volt << std::endl;
          return -1.;
         }
}
// --------------------------------------------------------------------------------------------------------------
double WaveForm::CFDToTLinear(std::vector<double> *w, std::vector<double> *t, float fraction)
{
    if (m_EndIndx == -99) m_EndIndx = EndIndx(w, true);
    if (m_maxIndx == -99) m_maxIndx = VoltMaxIndx(w, true);
    if (m_maxVolt == -99) m_maxVolt = VoltMax(w, true);
    double st_time = -99.;
    if (fraction != m_fraction || m_CFDTime == -99) st_time = CFDTimeLinear(w, t, fraction);
    else  st_time = m_CFDTime;
    double nd_time = SecondTimeForVoltage(w, t, fabs(fraction*m_maxVolt));
    if (nd_time != -99 && st_time != -99 && nd_time != st_time) return fabs(nd_time - st_time);
    else {
          if (m_WvBase->LGADBase::GetVerbosity() > 0 ) std::cout << __FUNCTION__ << " ERROR: Failed to calculate CFD ToT!" << std::endl;
          return -1.;
         }
}
// --------------------------------------------------------------------------------------------------------------
double WaveForm::TriggToTLinear(std::vector<double> *w, std::vector<double> *t, double trig)
{
    if (m_EndIndx == -99) m_EndIndx = EndIndx(w, true);
    if (m_maxIndx == -99) m_maxIndx = VoltMaxIndx(w, true);
    double t1 = -99.;
    if (m_TriggTime == -99 || m_trigg != trig) t1 = FirstTimeForVoltage(w, t, trig);
    else t1 = m_TriggTime;
    double t2 = SecondTimeForVoltage(w, t, trig);
    if (t1 != -99 && t2 !=-99 && t1 != t2) return fabs(t2 - t1); // Time over threshold
    else {
          if (m_WvBase->LGADBase::GetVerbosity() > 0) std::cout << __FUNCTION__ << " ERROR: Failed to calculate Trigger ToT!" << std::endl;
          return -1.;
         }
}
// --------------------------------------------------------------------------------------------------------------
bool WaveForm::dump()
{
    if (m_voltage == NULL || m_time == NULL)
       {
        std::cout << __FUNCTION__ << " ERROR: Time and Voltage not set for waveform, cannot perform any operation!" << std::endl;
        return false;
       }
    std::cout << "==> Dumping Waveform info" << std::endl;
    std::cout << "\t Rate: " << m_SnRate << " - check against rate from size of w/t vectors: "
              << (long long int)ceil(((Long64_t)(m_time->size() / abs(m_time->back() - m_time->front())))) << std::endl;
    std::cout <<
        "\t Polarity                  :" << m_pol << std::endl <<
        "\t Is Signal                 :" << m_IsSignal << std::endl <<
        "\t In Window                 :" << m_IsInWindow << std::endl <<
        "\t Average of Maxima         :" << m_maxavg << std::endl <<
        "\t Average of minima         :" << m_minavg << std::endl <<
        "\t Pedestal                  :" << m_pedestal << "+/-" << m_pedestalErr << std::endl <<
        "\t Noise                     :" << m_noise << "+/-" << m_noiseErr << std::endl <<
        "\t Max positipon             :" << m_maxIndx << std::endl <<
        "\t Max Time                  :" << m_maxTime << std::endl <<
        "\t Max Voltage               :" << m_maxVolt << "+/-" << m_maxVoltErr << std::endl <<
        "\t Min position              :" << m_minIndx << std::endl <<
        "\t Min Time                  :" << m_minTime << std::endl <<
        "\t Min Voltage               :" << m_minVolt << std::endl <<
        "\t Signal Start position     :" << m_StrIndx << std::endl <<
        "\t Signal Stop position      :" << m_EndIndx << std::endl <<
        "\t Collected Charge          :" << m_charge << std::endl <<
        "\t Risetime [10-90%]         :" << m_RiseTime << std::endl <<
        "\t CFD time at " << m_fraction << "           :" << m_CFDTime << std::endl <<
        "\t Trigger time at " << m_trigg << "      :" << m_TriggTime << std::endl <<
        "\t Maximum dV/dT             :" << m_dVdTMax << std::endl <<
        "\t dV/dT at CFD calue of " << m_fraction << " :" << m_dVdTCFD << std::endl <<
        "\t CFD ToT for " << m_fraction << "           :" << m_CFDToT << std::endl <<
        "\t Trigger ToT at " << m_trigg << " mV    :" << m_TriggToT << std::endl <<
        "\t Signal FFT                :" << m_signalFFT << std::endl <<
        "\t Noise FFT                 :" << m_noiseFFT << std::endl <<
        "\t Amplifier transimpedence  :" << m_transimp << " and gain: " << m_ampgain << std::endl <<
        "\t jitter [noise/dVdt]       :" << m_jitter1 << std::endl <<
        "\t jitter [Rise time/SNR]    :" << m_jitter2 << std::endl;

    return true;
}