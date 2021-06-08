/*
* LGADRun.cxx
*
*
*      Author: Gkougkousis Evangelos - Leonidas
*              egkougko@cern.ch
*               IFAE-BARCELONA
*/

#include "../LGADUtils/LGADRun.h"

//#if !defined(__CINT__)
//ClassImp(LGADRun);
//#endif
//
//#if !defined(__CLING__)
//ClassImp(LGADRun);
//#endif
//
//#ifdef __CINT__
//#pragma link C++ class LGADRun;
//#endif
//
//#ifdef __ROOTCLING__
//#pragma link C++ class LGADRun;
//#endif

LGADRun::LGADRun(TTree*)
{
    std::cout << setfill('-') << setw(89) << "-" << setfill(' ') << std::endl;
}
//------------------------------------------------------------------------------
LGADRun::LGADRun(LGADBase* rBase, TTree*)
{
    m_instrument = rBase->LGADBase::GetInstrument();
    m_TrackComb = rBase->LGADBase::GetTrackComb();
    m_channels = rBase->LGADBase::GetChannels();
    m_nchan = rBase->LGADBase::GetChNo();
    m_event = rBase->LGADBase::GetEvntNo();
    m_TransFile = rBase->LGADBase::GetTransFile();
    m_TransFileName = rBase->LGADBase::GetTransFileName();
    m_ofile = rBase->LGADBase::GetRootFile();
    m_tree = rBase->LGADBase::GetRootTree();
    m_WaveShape = rBase->LGADBase::GetWaveShape();
    m_RunBase = rBase;

    std::cout << setfill('-') << setw(89) << "-" << setfill(' ') << std::endl;
}
//------------------------------------------------------------------------------
LGADRun::~LGADRun()
{
    delete m_RunBase;
    for (unsigned int gh = 0; gh < m_RunDUTCh.size(); gh++) delete m_RunDUTCh.at(gh);
}
//------------------------------------------------------------------------------
bool LGADRun::SetChDQunats(std::vector<double> dQunat, std::vector<std::pair<std::pair<int, int>, double>> DQs)
{
    if (m_nchan == 0)
       {
        std::cout << __FUNCTION__ << " WARNING: No channels found!" << std::endl;
        return false;
       }
    else { 
          dQunat.clear();
          bool found = false;
          for (unsigned int k = 0; k < m_nchan; k++)
              {
               for (unsigned g = k + 1; g < m_nchan; g++)
                   {
                    found = false;
                    for (unsigned gh = 0; gh < DQs.size(); gh++)
                         {
                          if (((DQs.at(gh)).first).first == m_channels.at(k) && ((DQs.at(gh)).first).second == m_channels.at(g))
                             {
                              found = true;
                              dQunat.push_back((DQs.at(gh)).second);
                              break;
                             }
                         }
                    if (!found) dQunat.push_back(-99);
                   }
              }
          for (unsigned int hrt = 0; hrt < DQs.size(); hrt++)
              {
               if (((DQs.at(hrt)).first).first == -99 && ((DQs.at(hrt)).first).second == -99)
                  {
                   for (unsigned int kr = 0; kr < dQunat.size(); kr++)
                       {
                        if (dQunat.at(kr) == -99) dQunat.at(kr) = (DQs.at(hrt)).second;
                       }
                  }
              }
          for (unsigned int kp = 0; kp < dQunat.size(); kp++) if (dQunat.at(kp) == -99) dQunat.at(kp) = -1;
          return true;
         }
}
//------------------------------------------------------------------------------
void LGADRun::Begin(TTree*)
{
    printf("Client master begin!\n");
    TString option = GetOption();

    // Variable initialisation 

    std::cout << setfill('-') << setw(89) << "-" << std::endl << setfill(' ');
}
//------------------------------------------------------------------------------
void LGADRun::SlaveBegin(TTree*)
{
    printf("Client slave begin!\n");
    TString option = GetOption();
    m_init = 0;
    std::cout << setfill('-') << setw(89) << "-" << std::endl << setfill(' ');
}
//------------------------------------------------------------------------------
// The Notify() function is called when a new file is opened. This
// can be either for a new TTree in a TChain or when a new TTree is
// started when using PROOF. It is normally not necessary to make 
// changes to the generated code, but the routine can be extended  
// by the user if needed. The return value is currently not used.
// This is called twice in the begining, not sure yet why, possibly once for
// the cain and once for the tree
Bool_t LGADRun::Notify()
 {
  Info("Notify", "Opening a new file ...");

  if (fChain)                             // We are running locally... 
     {
      std::cout << "--> Running Localy" << std::endl;
      m_ofname = m_ofile->GetName(); 
      if (m_event == 0) m_event = fChain->GetTree()->GetEntriesFast();
     }
  else {                                   // We are running on PROOF:
        std::cout << "--> Running on PROOF" << std::endl;
        m_ofile = fChain->GetCurrentFile();
        m_ofname = fChain->GetName();
        if (m_event == 0) m_event = fChain->GetTree()->GetEntriesFast();
       }
  std::cout << "--> File Accessed: " << m_ofname << std::endl;
  std::cout << "--> Total number of entries in file: " << m_event << std::endl;

  m_RunBase->SetInstrument(m_instrument); // need to set the instrument for the prototype object
                                          // after selector automaticaly scans it 

  // should only be called once, not every time a new file is opened
  // Needs to be called AFTER selector has been initialised (not before Notify)
  if (m_init < 2) 
     {
      FirstInit();
      if ((m_RunBase->GetChDTs()).size() > 0) SetChDQunats(m_dTMax, m_RunBase->GetChDTs());
      if ((m_RunBase->GetChDCs()).size() > 0) SetChDQunats(m_dCMax, m_RunBase->GetChDCs());
     }

  std::cout << setfill('-') << setw(89) << "-" << std::endl << setfill(' ');
  return kTRUE;
 }
//------------------------------------------------------------------------------
Bool_t LGADRun::Process(Long64_t entry)
{
    if (fChain == 0) return kFALSE;

    b_EvnNo->GetEntry(entry);
    if (entry == 0)
       { 
        m_evnt1 = m_RunBase->LGADBase::GetStartStopEvnt().first;
        m_evnt2 = m_RunBase->LGADBase::GetStartStopEvnt().second;
        // outside of paramterts
        if (m_evnt2 == 0 || m_evnt2 <= m_evnt1 || m_evnt2 > (m_event + EvnNo)) 
           {
            if (m_RunBase->GetInstrument() == Sampic) m_evnt2 = (m_event + EvnNo)*m_nchan;
            else m_evnt2 = m_event + EvnNo;
           }
        else if (m_evnt2 <= EvnNo && (m_evnt2 - m_evnt1) <= m_event && m_evnt1 < m_evnt2) { m_evnt2 += EvnNo; m_evnt1 += EvnNo; }
        else if (m_evnt2 > EvnNo && m_evnt1 < EvnNo && (m_evnt2 - EvnNo) < m_event && m_evnt1 < m_evnt2) m_evnt1 += EvnNo;
       }
    else VarInit();

    if (m_evnt2 != 0 && m_evnt2 > m_evnt1 && EvnNo > m_evnt2) return false;
    m_RunBase->LGADBase::ProgressBar(entry, m_event);

    if (m_RunBase->LGADBase::GetVerbosity() >= 2) 
       {
        std::cout << "-------------------------------------------------------------------" << std::endl;
        std::cout << __FUNCTION__ << " INFO: Number of active channels: " << m_nchan << ", instrument: " << m_instrument << ", event: " << entry << std::endl;
        std::cout << __FUNCTION__ << " INFO: Channels: ";
        for (unsigned int k = 0; k < m_nchan; k++) std::cout << m_channels.at(k) << " ";
        std::cout << std::endl;
       }

    // Get general branches
    for (unsigned int ich = 0; ich < m_nchan; ich++)
        {
         b_t.at(ich)->GetEntry(entry);
         b_w.at(ich)->GetEntry(entry);
        }
    // Get instrument specific branches
    if (m_instrument == Sampic) 
       { 
        b_nPoints.at(0)->GetEntry(entry);
        b_SnRate.at(0)->GetEntry(entry);
        for (unsigned int ich = 0; ich < m_nchan; ich++)
            {
             m_RunDUTCh.at(ich)->DUTChannel::SetChRate(m_srate.at(0));
             b_ordrt.at(ich)->GetEntry(entry);
             b_physt.at(ich)->GetEntry(entry);
            }
       }
    else {
          for (unsigned int ich = 0; ich < m_nchan; ich++)
              {
               b_nPoints.at(ich)->GetEntry(entry);
               b_SnRate.at(ich)->GetEntry(entry);
               m_RunDUTCh.at(ich)->DUTChannel::SetChRate(m_srate.at(ich));
               b_vScale.at(ich)->GetEntry(entry);
               if (m_instrument == TektronixScope) b_TriggTime.at(ich)->GetEntry(entry);
              }
          if (m_instrument != InfiniiumScope) b_trigtime->GetEntry(entry);
         }
    if (m_RunBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Branches retrieved!" << std::endl;

    for (unsigned int ich = 0; ich < m_nchan; ich++)
        {
         // m_RunDUTCh.at(ich)->SetVerbose(2);
         //m_RunWaveform.at(ich)->WaveForm::InitializeWaveForm(1);
         //m_RunWaveform.at(ich)->WaveForm::SetVoltage(w.at(ich));
         //m_RunWaveform.at(ich)->WaveForm::SetTime(t.at(ich));
         //m_RunWaveform.at(ich)->WaveForm::Calculate();
         m_RunDUTCh.at(ich)->DUTChannel::AppendEvent(w.at(ich), t.at(ich));
         if (m_RunBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Event Added to channel: " << ich << "/" << m_nchan << "!" << std::endl;

         // Get Event Variables
         m_EvPol.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetEvntPolarity();
         m_EvCharge.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetEvntCharge();
         m_EvIsSignal.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetEvntIsSignal();
         m_EvInWindow.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetEvntInWindow();
         m_EvVoltStr.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetEvntVoltStr();
         m_EvMaxindx.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetEvntMaxindx();
         m_EvMaxTime.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetEvntMaxTime();
         m_EvMaxVolt.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetEvntMaxVolt();
         m_EvMinIndx.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetEvntMinIndx();
         m_EvMinTime.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetEvntMinTime();
         m_EvMinVolt.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetEvntMinVolt();
         m_EvStrIndx.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetEvntStrIndx();
         m_EvEndIndx.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetEvntEndIndx();
         m_EvNoise.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetEvntNoise();
         m_EvPedest.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetEvntPedest();
         m_EvNoiseErr.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetEvntNoiseErr();
         m_EvPedestErr.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetEvntPedestErr();
         m_EvRiseT.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetEvntRiseT();
         m_EvJitRiseSNR.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetEvntJitRiseSNR();
         m_EvTriggTime.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetEvntTriggTime();
         m_EvDVDTMax.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetEvntDVDTMax();
         m_EvTriggToT.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetEvntTriggToT();
         m_EvSignalFFT.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetEvntSignalFFT();
         m_EvNoiseFFT.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetEvntNoiseFFT();
         m_EvComplete.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetEvntComplete();
         m_EvVAdjCFD.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetEvntVAdjCFD();
         if (m_RunBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Got event variables: " << ich+1 << "/" << m_nchan << "!" << std::endl;

         for (unsigned int bin = 0; bin < 19; bin++)
             { 
              (m_EvCFDTime.at(ich)).at(bin) = m_RunDUTCh.at(ich)->DUTChannel::GetEvntCFDTime(bin);
              (m_EvDVDTCFD.at(ich)).at(bin) = m_RunDUTCh.at(ich)->DUTChannel::GetEvntDVDTCFD(bin);
              (m_EvCFDToT.at(ich)).at(bin) = m_RunDUTCh.at(ich)->DUTChannel::GetEvntCFDToT(bin);
              (m_EvJitNdVdT.at(ich)).at(bin) = m_RunDUTCh.at(ich)->DUTChannel::GetEvntJitNdVdT(bin);
             }
         if (m_RunBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Got CFD points: " << ich+1 << "/" << m_nchan << "!" << std::endl;

         // Write ntuple branches
         b_EvPol.at(ich)->Fill();
         b_EvCharge.at(ich)->Fill();
         b_EvIsSignal.at(ich)->Fill();
         b_EvInWindow.at(ich)->Fill();
         b_EvVoltStr.at(ich)->Fill();
         b_EvMaxindx.at(ich)->Fill();
         b_EvMaxTime.at(ich)->Fill();
         b_EvMaxVolt.at(ich)->Fill();
         b_EvMinIndx.at(ich)->Fill();
         b_EvMinTime.at(ich)->Fill();
         b_EvMinVolt.at(ich)->Fill();
         b_EvStrIndx.at(ich)->Fill();
         b_EvEndIndx.at(ich)->Fill();
         b_EvNoise.at(ich)->Fill();
         b_EvPedest.at(ich)->Fill();
         b_EvNoiseErr.at(ich)->Fill();
         b_EvPedestErr.at(ich)->Fill();
         b_EvRiseT.at(ich)->Fill();
         b_EvJitNdVdT.at(ich)->Fill();
         b_EvJitRiseSNR.at(ich)->Fill();
         b_EvTriggTime.at(ich)->Fill();
         b_EvDVDTMax.at(ich)->Fill();
         b_EvDVDTCFD.at(ich)->Fill();
         b_EvCFDToT.at(ich)->Fill();
         b_EvTriggToT.at(ich)->Fill();
         b_EvSignalFFT.at(ich)->Fill();
         b_EvNoiseFFT.at(ich)->Fill();
         b_EvComplete.at(ich)->Fill();
         b_EvCFDTime.at(ich)->Fill();
         b_EvVAdjCFD.at(ich)->Fill();
         if (m_RunBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Written all brancehs: " << ich << "/" << m_nchan << "!" << std::endl;
        }

    bool dTCutPass = false;
    bool dChCutPass = false;
    for (unsigned int ic1 = 0; ic1 < m_nchan; ic1++)
        {
         for (unsigned int ic2 = ic1+1; ic2 < m_nchan; ic2++)
             {
              double dCh = fabs(m_EvCharge.at(ic1) - m_EvCharge.at(ic2));
              for (unsigned int bin1 = 0; bin1 < 19; bin1++)
                  {
                   for (unsigned int bin2 = 0; bin2 < 19; bin2++)
                      {
                       dTCutPass = false;
                       dChCutPass = false;
                       // populate the time difference vectors
                       if ((m_EvCFDTime.at(ic1)).at(bin1) != -1 && (m_EvCFDTime.at(ic2)).at(bin2) != -1
                            && (m_EvCFDTime.at(ic1)).at(bin1) != -99 && (m_EvCFDTime.at(ic2)).at(bin2) != -99)
                          {
                           double dT = fabs((m_EvCFDTime.at(ic1)).at(bin1) - (m_EvCFDTime.at(ic2)).at(bin2));
                           if (m_dTMax.size() == 0) dTCutPass = true;
                           else if (m_dTMax.at((ic1*(ic2-1))/2) != -1 && fabs(dT) < m_dTMax.at((ic1*(ic2-1))/2)) dTCutPass = true;
                           else if (m_dTMax.at((ic1*(ic2-1))/2) == -1) dTCutPass = true;
                           if (m_dCMax.size() == 0) dChCutPass = true;
                           else if (m_dCMax.at((ic1*(ic2-1))/2) != -1 && fabs(dCh) < m_dCMax.at((ic1*(ic2-1))/2)) dChCutPass = true;
                           else if (m_dCMax.at((ic1*(ic2-1))/2) == -1)  dChCutPass = true;
                           if (dTCutPass && dChCutPass) ((m_EvTmDiff.at((ic1*(ic2-1))/2)).at(bin1*19+bin2)).push_back(dT);
                          }
                      }
                  }
             }
        }

    return kTRUE;
}
//------------------------------------------------------------------------------
void LGADRun::SlaveTerminate()
{
    std::cout << setfill('-') << setw(89) << "-" << std::endl << setfill(' '); 

    m_ofile->cd();
    m_tree->Write("", TObject::kOverwrite);

    // Tree for Channel propertied form Fits
    TTree *ChFitParam = new TTree("ChFitParam", "Channel parameters form fits");
    TTree *ChShapeParam = nullptr;
    ChFitParam->SetAutoFlush(10000); // Set autoflush to 30 MB

    TList *histoslist = new TList();
    std::vector<TList *> ChHist;
    std::vector<TList *> CFDTimeHist;
    std::vector<TList *> CFDToTHist;
    std::vector<TList *> CFDdVdTHist;
    std::vector<TList *> TimeDiff;

    int chName = 0;
    // Channel properties from respective fits
    std::pair <double, double> ChFitMxVlt;
    std::pair <double, double> ChFitMinVlt;
    std::pair <int, int> ChFitMxIndx;
    std::pair <int, int> ChFitMinIndx;
    std::pair <int, int> ChFitStrIndx;
    std::pair <int, int> ChFitEndIndx;
    std::pair <double, double> ChFitIsSig;
    std::pair <double, double> ChFitInWin;
    std::pair <double, double> ChFitVStr;
    std::pair <double, double> ChFitNoise;
    std::pair <double, double> ChFitPedestal;
    std::pair <double, double> ChFitNoiseSlp;
    std::pair <double, double> ChFitPedestalSlp;
    std::pair <double, double> ChFitNoiseErr;
    std::pair <double, double> ChFitPedestalErr;
    std::pair <double, double> ChFitMaxTime;
    std::pair <double, double> ChFitMinTime;
    std::pair <double, double> ChFitCharge;
    std::pair <double, double> ChFitRiseT;
    std::pair <double, double> ChFitJitRiseSNR;
    std::vector<std::pair <double, double>> ChFitCFDTime;
    std::vector<std::pair <double, double>> ChFitDVDTCFD;
    std::vector<std::pair <double, double>> ChFitCFDToT;
    std::vector<std::pair <double, double>> ChFitJitNdVdT;
    std::pair <double, double> ChFitTriggTime;
    std::pair <double, double> ChFitDVDTMax;
    std::pair <double, double> ChFitTriggToT;
    std::pair <double, double> ChFitSignalFFT;
    std::pair <double, double> ChFitNoiseFFT;
    std::pair <double, double> ChFitComplete;
    std::pair <double, double> ChFitSoNR;
    std::vector<std::vector<std::pair<double, double>>> TmDiffFitQts;

    // Chi2 of the Fits for channel properties
    double MaxIndxFtChi2;
    double MinIndxFtChi2;
    double MaxVoltFtChi2;
    double MinVoltFtChi2;
    double StrIndxFtChi2;
    double EndIndxFtChi2;
    double NoiseFtChi2;
    double NoiseErrFtChi2;
    double PedestFtChi2;
    double PedestErrFtChi2;
    double MaxTimeFtChi2;
    double MinTimeFtChi2;
    double ChargeFtChi2;
    double RiseTimeFtChi2;
    double TriggTimeFtChi2;
    double DVDTMaxFtChi2;
    double TriggToTFtChi2;
    double SignalFFTFtChi2;
    double NoiseFFTFtChi2;
    std::vector<double> CFDTimeFtChi2;
    std::vector<double> DVDTCFDFtChi2;
    std::vector<double> CFDToTFtChi2;
    std::vector<std::vector<double>> TmDiffChi2;

    ChFitParam->Branch("chName", &chName, "chName/i");

    ChFitParam->Branch("ChFitMaxVolt", &ChFitMxVlt.first, "ChFitMxVlt.first/D");
    ChFitParam->Branch("ChFitMaxVoltErr", &ChFitMxVlt.second, "ChFitMxVlt.second/D");
    ChFitParam->Branch("ChFitMinVolt", &ChFitMinVlt.first, "ChFitMinVlt.first/D");
    ChFitParam->Branch("ChFitMinVoltErr", &ChFitMinVlt.second, "ChFitMinVlt.second/D");
    ChFitParam->Branch("ChFitMxIndx", &ChFitMxIndx.first, "ChFitMxIndx.first/I");
    ChFitParam->Branch("ChFitMxIndxErr", &ChFitMxIndx.second, "ChFitMxIndx.second/I");
    ChFitParam->Branch("ChFitMinIndx", &ChFitMinIndx.first, "ChFitMinIndx.first/I");
    ChFitParam->Branch("ChFitMinIndxErr", &ChFitMinIndx.second, "ChFitMinIndx.second/I");
    ChFitParam->Branch("ChFitStrIndx", &ChFitStrIndx.first, "ChFitStrIndx.first/D");
    ChFitParam->Branch("ChFitStrIndxErr", &ChFitEndIndx.second, "ChFitEndIndx.second/D");
    ChFitParam->Branch("ChFitEndIndx", &ChFitEndIndx.first, "ChFitEndIndx.first/D");
    ChFitParam->Branch("ChFitEndIndxErr", &ChFitEndIndx.second, "ChFitEndIndx.second/D");
    ChFitParam->Branch("ChFitIsSig", &ChFitIsSig.first, "ChFitIsSig.first/D");
    ChFitParam->Branch("ChFitIsSigErr", &ChFitIsSig.second, "ChFitIsSig.second/D");
    ChFitParam->Branch("ChFitInWin", &ChFitInWin.first, "ChFitInWin.first/D");
    ChFitParam->Branch("ChFitInWinErr", &ChFitInWin.second, "ChFitInWin.second/D");
    ChFitParam->Branch("ChFitVStr", &ChFitVStr.first, "ChFitVStr.first/D");
    ChFitParam->Branch("ChFitVStrErr", &ChFitVStr.second, "ChFitVStr.second/D");
    ChFitParam->Branch("ChFitNoise", &ChFitNoise.first, "ChFitNoise.first/D");
    ChFitParam->Branch("ChFitNoiseErr", &ChFitNoise.second, "ChFitNoise.second/D");
    ChFitParam->Branch("ChFitPedestal", &ChFitPedestal.first, "ChFitPedestal.first/D");
    ChFitParam->Branch("ChFitPedestalErr", &ChFitPedestal.second, "ChFitPedestal.second/D");
    ChFitParam->Branch("ChFitNoiseSlp", &ChFitNoiseSlp.first, "ChFitNoiseSlp.first/D");
    ChFitParam->Branch("ChFitNoiseSlpErr", &ChFitNoiseSlp.second, "ChFitNoiseSlp.second/D");
    ChFitParam->Branch("ChFitPedestalSlp", &ChFitPedestalSlp.first, "ChFitPedestalSlp.first/D");
    ChFitParam->Branch("ChFitPedestalSlpErr", &ChFitPedestalSlp.second, "ChFitPedestalSlp.second/D");
    ChFitParam->Branch("ChFitNoiseErr", &ChFitNoiseErr.first, "ChFitNoiseErr.first/D");
    ChFitParam->Branch("ChFitNoiseErrUnc", &ChFitNoiseErr.second, "ChFitNoiseErr.second/D");
    ChFitParam->Branch("ChFitPedestalErr", &ChFitPedestalErr.first, "ChFitPedestalErr.first/D");
    ChFitParam->Branch("ChFitPedestalErrUnc", &ChFitPedestalErr.second, "ChFitPedestalErr.second/D");
    ChFitParam->Branch("ChFitMaxTime", &ChFitMaxTime.first, "ChFitMaxTime.first/D");
    ChFitParam->Branch("ChFitMaxTimeErr", &ChFitMaxTime.second, "ChFitMaxTime.second/D");
    ChFitParam->Branch("ChFitMinTime", &ChFitMinTime.first, "ChFitMinTime.first/D");
    ChFitParam->Branch("ChFitMinTimeErr", &ChFitMinTime.second, "ChFitMinTime.second/D");
    ChFitParam->Branch("ChFitCharge", &ChFitCharge.first, "ChFitCharge.first/D");
    ChFitParam->Branch("ChFitChargeErr", &ChFitCharge.second, "ChFitCharge.second/D");
    ChFitParam->Branch("ChFitRiseT", &ChFitRiseT.first, "ChFitRiseT.first/D");
    ChFitParam->Branch("ChFitRiseTErr", &ChFitRiseT.second, "ChFitRiseT.second/D");
    ChFitParam->Branch("ChFitJitRiseSNR", &ChFitJitRiseSNR.first, "ChFitJitRiseSNR.first/D");
    ChFitParam->Branch("ChFitJitRiseSNRErr", &ChFitJitRiseSNR.second, "ChFitJitRiseSNR.second/D");
    ChFitParam->Branch("ChFitCFDTime", &ChFitCFDTime);
    ChFitParam->Branch("ChFitDVDTCFD", &ChFitDVDTCFD);
    ChFitParam->Branch("ChFitCFDToT", &ChFitCFDToT);
    ChFitParam->Branch("ChFitJitNdVdT", &ChFitJitNdVdT);
    ChFitParam->Branch("TmDiffFitQts", &TmDiffFitQts);
    ChFitParam->Branch("ChFitTriggTime", &ChFitTriggTime.first, "ChFitTriggTime.first/D");
    ChFitParam->Branch("ChFitTriggTimeErr", &ChFitTriggTime.second, "ChFitTriggTime.second/D");
    ChFitParam->Branch("ChFitDVDTMax", &ChFitDVDTMax.first, "ChFitDVDTMax.first/D");
    ChFitParam->Branch("ChFitDVDTMaxErr", &ChFitDVDTMax.second, "ChFitDVDTMax.second/D");
    ChFitParam->Branch("ChFitTriggToT", &ChFitTriggToT.first, "ChFitTriggToT.first/D");
    ChFitParam->Branch("ChFitTriggToTErr", &ChFitTriggToT.second, "ChFitTriggToT.second/D");
    ChFitParam->Branch("ChFitSignalFFT", &ChFitSignalFFT.first, "ChFitSignalFFT.first/D");
    ChFitParam->Branch("ChFitSignalFFTErr", &ChFitSignalFFT.second, "ChFitSignalFFT.second/D");
    ChFitParam->Branch("ChFitNoiseFFT", &ChFitNoiseFFT.first, "ChFitNoiseFFT.first/D");
    ChFitParam->Branch("ChFitNoiseFFTErr", &ChFitNoiseFFT.second, "ChFitNoiseFFT.second/D");
    ChFitParam->Branch("ChFitComplete", &ChFitComplete.first, "ChFitComplete.first/D");
    ChFitParam->Branch("ChFitCompleteErr", &ChFitComplete.second, "ChFitComplete.second/D");
    ChFitParam->Branch("ChFitSoNR", &ChFitSoNR.first, "ChFitSoNR.first/D");
    ChFitParam->Branch("ChFitSoNRErr", &ChFitSoNR.second, "ChFitSoNR.second/D");

    ChFitParam->Branch("MaxIndxFtChi2", &MaxIndxFtChi2, "MaxIndxFtChi2/D");
    ChFitParam->Branch("MinIndxFtChi2", &MinIndxFtChi2, "MinIndxFtChi2/D");
    ChFitParam->Branch("MaxVoltFtChi2", &MaxVoltFtChi2, "MaxVoltFtChi2/D");
    ChFitParam->Branch("MinVoltFtChi2", &MinVoltFtChi2, "MinVoltFtChi2/D");
    ChFitParam->Branch("StrIndxFtChi2", &StrIndxFtChi2, "StrIndxFtChi2/D");
    ChFitParam->Branch("EndIndxFtChi2", &EndIndxFtChi2, "EndIndxFtChi2/D");
    ChFitParam->Branch("NoiseFtChi2", &NoiseFtChi2, "NoiseFtChi2/D");
    ChFitParam->Branch("NoiseErrFtChi2", &NoiseErrFtChi2, "NoiseErrFtChi2/D");
    ChFitParam->Branch("PedestFtChi2", &PedestFtChi2, "PedestFtChi2/D");
    ChFitParam->Branch("PedestErrFtChi2", &PedestErrFtChi2, "PedestErrFtChi2/D");
    ChFitParam->Branch("MaxTimeFtChi2", &MaxTimeFtChi2, "MaxTimeFtChi2/D");
    ChFitParam->Branch("MinTimeFtChi2", &MinTimeFtChi2, "MinTimeFtChi2/D");
    ChFitParam->Branch("ChargeFtChi2", &ChargeFtChi2, "ChargeFtChi2/D");
    ChFitParam->Branch("RiseTimeFtChi2", &RiseTimeFtChi2, "RiseTimeFtChi2/D");
    ChFitParam->Branch("TriggTimeFtChi2", &TriggTimeFtChi2, "TriggTimeFtChi2/D");
    ChFitParam->Branch("DVDTMaxFtChi2", &DVDTMaxFtChi2, "DVDTMaxFtChi2/D");
    ChFitParam->Branch("TriggToTFtChi2", &TriggToTFtChi2, "TriggToTFtChi2/D");
    ChFitParam->Branch("SignalFFTFtChi2", &SignalFFTFtChi2, "SignalFFTFtChi2/D");
    ChFitParam->Branch("NoiseFFTFtChi2", &NoiseFFTFtChi2, "NoiseFFTFtChi2/D");
    ChFitParam->Branch("CFDTimeFtChi2", &CFDTimeFtChi2);
    ChFitParam->Branch("DVDTCFDFtChi2", &DVDTCFDFtChi2);
    ChFitParam->Branch("CFDToTFtChi2", &CFDToTFtChi2);
    ChFitParam->Branch("TmDiffChi2", &TmDiffChi2);

    // Channel Shape properties
    int ChShapePolarity;
    int ChShapeMaxIndx;
    double ChShapeMax;
    int ChShapeMinIndx;
    double ChShapeMin;
    int ChShapeStrIndx;
    int ChShapeEndIndx;
    double ChShapeMaxTime;
    double ChShapeMinTime;
    std::pair <double, double> ChShapeNoise;
    double ChShapeCharge;
    double ChShapeRiseTime;
    double ChShapeJitRiseSNR;
    double ChShapeTriggTime;
    double ChShapeDVDTMax;
    double ChShapeTriggToT;
    double ChShapeSignalFFT;
    double ChShapeNoiseFFT;
    std::pair <double, double> ChShapeSoN;
    std::vector<double> ChShapeCFDTime;
    std::vector<double> ChShapeDVDTCFD;
    std::vector<double> ChShapeCFDToT;
    std::vector<double> ChShapeJitNdVdT;
    std::vector<std::pair <double, double>> ChMeanVolt;

    if (m_WaveShape)
       {
        ChShapeParam = new TTree("ChShapeParam", "Channel parameters form wwaveform shape");
        ChShapeParam->SetAutoFlush(10000); // Set autoflush to 30 MB

        ChShapeParam->Branch("chName", &chName, "chName/i");
        ChShapeParam->Branch("ChShapePolarity", &ChShapePolarity, "ChShapePolarity/I");
        ChShapeParam->Branch("ChShapeMaxIndx", &ChShapeMaxIndx, "ChShapeMaxIndx/I");
        ChShapeParam->Branch("ChShapeMax", &ChShapeMax, "ChShapeMax/D");
        ChShapeParam->Branch("ChShapeMinIndx", &ChShapeMinIndx, "ChShapeMinIndx/I");
        ChShapeParam->Branch("ChShapeMin", &ChShapeMin, "ChShapeMin/D");
        ChShapeParam->Branch("ChShapeStrIndx", &ChShapeStrIndx, "ChShapeStrIndx/I");
        ChShapeParam->Branch("ChShapeEndIndx", &ChShapeEndIndx, "ChShapeEndIndx/I");
        ChShapeParam->Branch("ChShapeMaxTime", &ChShapeMaxTime, "ChShapeMaxTime/D");
        ChShapeParam->Branch("ChShapeMinTime", &ChShapeMinTime, "ChShapeMinTime/D");
        ChShapeParam->Branch("ChShapeNoise", &ChShapeNoise.first, "ChShapeNoise.first/D");
        ChShapeParam->Branch("ChShapeNoiseErr", &ChShapeNoise.second, "ChShapeNoise.second/D");
        ChShapeParam->Branch("ChShapeCharge", &ChShapeCharge, "ChShapeCharge/D");
        ChShapeParam->Branch("ChShapeRiseTime", &ChShapeRiseTime, "ChShapeRiseTime/D");
        ChShapeParam->Branch("ChShapeJitRiseSNR", &ChShapeJitRiseSNR, "ChShapeJitRiseSNR/D");
        ChShapeParam->Branch("ChShapeTriggTime", &ChShapeTriggTime, "ChShapeTriggTime/D");
        ChShapeParam->Branch("ChShapeDVDTMax", &ChShapeDVDTMax, "ChShapeDVDTMax/D");
        ChShapeParam->Branch("ChShapeTriggToT", &ChShapeTriggToT, "ChShapeTriggToT/D");
        ChShapeParam->Branch("ChShapeSignalFFT", &ChShapeSignalFFT, "ChShapeSignalFFT/D");
        ChShapeParam->Branch("ChShapeNoiseFFT", &ChShapeNoiseFFT, "ChShapeNoiseFFT/D");
        ChShapeParam->Branch("ChShapeSoN", &ChShapeSoN.first, "ChShapeSoN.first/D");
        ChShapeParam->Branch("ChShapeSoNUnc", &ChShapeSoN.second, "ChShapeSoN.second/D");
        ChShapeParam->Branch("ChShapeCFDTime", &ChShapeCFDTime);
        ChShapeParam->Branch("ChShapeDVDTCFD", &ChShapeDVDTCFD);
        ChShapeParam->Branch("ChShapeCFDToT", &ChShapeCFDToT);
        ChShapeParam->Branch("ChShapeJitNdVdT", &ChShapeJitNdVdT);
        ChShapeParam->Branch("ChMeanVolt", &ChMeanVolt);
       }

    for (unsigned int ich = 0; ich < m_nchan; ich++)
        {
         std::cout << __FUNCTION__ << " INFO: Updating Channel " << ich + 1 << " of " << m_nchan << " designated: " << m_channels.at(ich) << std::endl;
         m_RunDUTCh.at(ich)->DUTChannel::updateChProperties(m_WaveShape, m_tree);
         if (m_RunBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: After channel update: " << ich + 1 << "/" << m_nchan << std::endl;

         ChHist.push_back(new TList());
         CFDTimeHist.push_back(new TList());
         CFDToTHist.push_back(new TList());
         CFDdVdTHist.push_back(new TList());
         for (unsigned int d = ich + 1; d < m_nchan; d++) TimeDiff.push_back(new TList());

         chName = m_channels.at(ich);

         // Channel Fit variables initialization part
         ChFitMxVlt = std::make_pair(-99, -99);
         ChFitMinVlt = std::make_pair(-99, -99);
         ChFitMxIndx = std::make_pair(-99, -99);
         ChFitMinIndx = std::make_pair(-99, -99);
         ChFitStrIndx = std::make_pair(-99, -99);
         ChFitEndIndx = std::make_pair(-99, -99);
         ChFitIsSig = std::make_pair(-99, -99);
         ChFitInWin = std::make_pair(-99, -99);
         ChFitVStr = std::make_pair(-99, -99);
         ChFitNoise = std::make_pair(-99, -99);
         ChFitPedestal = std::make_pair(-99, -99);
         ChFitNoiseSlp = std::make_pair(-99, -99);
         ChFitPedestalSlp = std::make_pair(-99, -99);
         ChFitNoiseErr = std::make_pair(-99, -99);
         ChFitPedestalErr = std::make_pair(-99, -99);
         ChFitMaxTime = std::make_pair(-99, -99);
         ChFitMinTime = std::make_pair(-99, -99);
         ChFitCharge = std::make_pair(-99, -99);
         ChFitRiseT = std::make_pair(-99, -99);
         ChFitJitRiseSNR = std::make_pair(-99, -99);
         ChFitCFDTime.clear();
         ChFitCFDTime.reserve(19);
         ChFitDVDTCFD.clear();
         ChFitDVDTCFD.reserve(19);
         ChFitCFDToT.clear();
         ChFitCFDToT.reserve(19);
         ChFitJitNdVdT.clear();
         ChFitJitNdVdT.reserve(19);
         ChFitTriggTime = std::make_pair(-99, -99);
         ChFitDVDTMax = std::make_pair(-99, -99);
         ChFitTriggToT = std::make_pair(-99, -99);
         ChFitSignalFFT = std::make_pair(-99, -99);
         ChFitNoiseFFT = std::make_pair(-99, -99);
         ChFitComplete = std::make_pair(-99, -99);
         ChFitSoNR = std::make_pair(-99, -99);

         int iter = (m_nchan-(ich+1));
         TmDiffChi2.clear();
         TmDiffFitQts.clear();
         TmDiffChi2.resize(iter);
         TmDiffFitQts.resize(iter);
         for (unsigned int k = 0; k < iter; k++)
             {
              (TmDiffChi2.at(k)).resize(19*19, -99.0);
              (TmDiffFitQts.at(k)).resize(19*19, std::make_pair(-99, -99));
             }

         MaxIndxFtChi2 = -99.;
         MinIndxFtChi2 = -99.;
         MaxVoltFtChi2 = -99.;
         MinVoltFtChi2 = -99.;
         StrIndxFtChi2 = -99.;
         EndIndxFtChi2 = -99.;
         NoiseFtChi2 = -99.;
         NoiseErrFtChi2 = -99.;
         PedestFtChi2 = -99.;
         PedestErrFtChi2 = -99.;
         MaxTimeFtChi2 = -99.;
         MinTimeFtChi2 = -99.;
         ChargeFtChi2 = -99.;
         RiseTimeFtChi2 = -99.;
         TriggTimeFtChi2 = -99.;
         DVDTMaxFtChi2 = -99.;
         TriggToTFtChi2 = -99.;
         SignalFFTFtChi2 = -99.;
         NoiseFFTFtChi2 = -99.;
         CFDTimeFtChi2.resize(19, -99);
         DVDTCFDFtChi2.resize(19, -99);
         CFDToTFtChi2.resize(19, -99);

         if (m_RunBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Iinitialized fit variables for channel: " << ich +1 << "/" << m_nchan << std::endl;

         ChFitMxVlt = m_RunDUTCh.at(ich)->DUTChannel::GetChMaxVolt(2);
         ChFitMinVlt = m_RunDUTCh.at(ich)->DUTChannel::GetChMinVolt(2);
         ChFitMxIndx = m_RunDUTCh.at(ich)->DUTChannel::GetChMaxIndx(2);
         ChFitMinIndx = m_RunDUTCh.at(ich)->DUTChannel::GetChMinIndx(2);
         ChFitStrIndx = m_RunDUTCh.at(ich)->DUTChannel::GetChStrIndx(2);
         ChFitEndIndx = m_RunDUTCh.at(ich)->DUTChannel::GetChEndIndx(2);
         ChFitIsSig = m_RunDUTCh.at(ich)->DUTChannel::GetChIsSignal();
         ChFitInWin = m_RunDUTCh.at(ich)->DUTChannel::GetChInWindow();
         ChFitVStr = m_RunDUTCh.at(ich)->DUTChannel::GetChVoltStr();
         ChFitNoise = m_RunDUTCh.at(ich)->DUTChannel::GetChNoise(2);
         ChFitPedestal = m_RunDUTCh.at(ich)->DUTChannel::GetChPedestal();
         ChFitNoiseSlp = m_RunDUTCh.at(ich)->DUTChannel::GetChNoiseSlope();
         ChFitPedestalSlp = m_RunDUTCh.at(ich)->DUTChannel::GetChPedestalSlope();
         ChFitNoiseErr = m_RunDUTCh.at(ich)->DUTChannel::GetChNoiseErr(2);
         ChFitPedestalErr = m_RunDUTCh.at(ich)->DUTChannel::GetChPedestalErr(2);
         ChFitMaxTime = m_RunDUTCh.at(ich)->DUTChannel::GetChMaxTime(2);
         ChFitMinTime = m_RunDUTCh.at(ich)->DUTChannel::GetChMinTime(2);
         ChFitCharge = m_RunDUTCh.at(ich)->DUTChannel::GetChCharge(2);
         ChFitRiseT = m_RunDUTCh.at(ich)->DUTChannel::GetChRiseTime(2);
         ChFitJitRiseSNR = m_RunDUTCh.at(ich)->DUTChannel::GetChJitRiseSNR(2);
         ChFitCFDTime = m_RunDUTCh.at(ich)->DUTChannel::GetChCFDMag(2, "CFDTime");
         ChFitDVDTCFD = m_RunDUTCh.at(ich)->DUTChannel::GetChCFDMag(2, "DVDTCFD");
         ChFitCFDToT = m_RunDUTCh.at(ich)->DUTChannel::GetChCFDMag(2, "CFDToT");
         ChFitJitNdVdT = m_RunDUTCh.at(ich)->DUTChannel::GetChCFDMag(2, "JitterNdVdT");
         ChFitTriggTime = m_RunDUTCh.at(ich)->DUTChannel::GetChTriggTime(2);
         ChFitDVDTMax = m_RunDUTCh.at(ich)->DUTChannel::GetChDVDTMax(2);
         ChFitTriggToT = m_RunDUTCh.at(ich)->DUTChannel::GetChTriggToT(2);
         ChFitSignalFFT = m_RunDUTCh.at(ich)->DUTChannel::GetChSignalFFT(2);
         ChFitNoiseFFT = m_RunDUTCh.at(ich)->DUTChannel::GetChNoiseFFT(2);
         ChFitComplete = m_RunDUTCh.at(ich)->DUTChannel::GetChComplete();
         ChFitSoNR = m_RunDUTCh.at(ich)->DUTChannel::GetChSoNR(2);
         if (m_RunBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Recovered fit variables for channel: " << ich +1 << "/" << m_nchan << std::endl;

         // Get the fits  and relevent chi2 to the Ntuple
         MaxIndxFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("MaxIndx");
         if (MaxIndxFtChi2 != -99 && MaxIndxFtChi2 != -1)
            {
             h_MaxIndxFt.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetFit("MaxIndx");
             h_MaxIndxFt.at(ich)->SetName(Form("MaxIndx%02u", m_channels.at(ich)));
             ChHist.at(ich)->Add(h_MaxIndxFt.at(ich));
            }
         MinIndxFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("MinIndx");
         if (MinIndxFtChi2 != -99 && MinIndxFtChi2 != -1)
            {
             h_MinIndxFt.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetFit("MinIndx");
             h_MinIndxFt.at(ich)->SetName(Form("MinIndx%02u", m_channels.at(ich)));
             ChHist.at(ich)->Add(h_MinIndxFt.at(ich));
            }
         MaxVoltFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("MaxVolt");
         if (MaxVoltFtChi2 != -99 && MaxVoltFtChi2 != -1)
            {
             h_MaxVoltFt.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetFit("MaxVolt");
             h_MaxVoltFt.at(ich)->SetName(Form("MaxVolt%02u", m_channels.at(ich)));
             ChHist.at(ich)->Add(h_MaxVoltFt.at(ich));
            }
         MinVoltFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("MinVotl");
         if (MinVoltFtChi2 != -99 && MinVoltFtChi2 != -1)
            {
             h_MinVoltFt.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetFit("MinVotl");
             h_MinVoltFt.at(ich)->SetName(Form("MinVotl%02u", m_channels.at(ich)));
             ChHist.at(ich)->Add(h_MinVoltFt.at(ich));
            }
         StrIndxFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("Strindx");
         if (StrIndxFtChi2 != -99 && StrIndxFtChi2 != -1)
            {
             h_StrIndxFt.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetFit("Strindx");
             h_StrIndxFt.at(ich)->SetName(Form("Strindx%02u", m_channels.at(ich)));
             ChHist.at(ich)->Add(h_StrIndxFt.at(ich));
            }
         EndIndxFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("EndIndx");
         if (EndIndxFtChi2 != -99 && EndIndxFtChi2 != -1)
            {
             h_EndIndxFt.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetFit("EndIndx");
             h_EndIndxFt.at(ich)->SetName(Form("EndIndx%02u", m_channels.at(ich)));
             ChHist.at(ich)->Add(h_EndIndxFt.at(ich));
            }
         NoiseFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("Noise");
         if (NoiseFtChi2 != -99 && NoiseFtChi2 != -1)
            {
             h_NoiseFt.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetFit("Noise");
             h_NoiseFt.at(ich)->SetName(Form("Noise%02u", m_channels.at(ich)));
             ChHist.at(ich)->Add(h_NoiseFt.at(ich));
            }
         NoiseErrFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("NoisErr");
         if (NoiseErrFtChi2 != -99 && NoiseErrFtChi2 != -1)
            {
             h_NoiseErrFt.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetFit("NoisErr");
             h_NoiseErrFt.at(ich)->SetName(Form("NoisErr%02u", m_channels.at(ich)));
             ChHist.at(ich)->Add(h_NoiseErrFt.at(ich));
            }
         PedestFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("Pedestl");
         if (PedestFtChi2 != -99 && PedestFtChi2 != -1)
            {
             h_PedestFt.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetFit("Pedestl");
             h_PedestFt.at(ich)->SetName(Form("Pedestl%02u", m_channels.at(ich)));
             ChHist.at(ich)->Add(h_PedestFt.at(ich));
            }
         PedestErrFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("PdslErr");
         if (PedestErrFtChi2 != -99 && PedestErrFtChi2 != -1)
            {
             h_NPedestErrFt.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetFit("PdslErr");
             h_NPedestErrFt.at(ich)->SetName(Form("PdslErr%02u", m_channels.at(ich)));
             ChHist.at(ich)->Add(h_NPedestErrFt.at(ich));
            }
         MaxTimeFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("MaxTime");
         if (MaxTimeFtChi2 != -99 && MaxTimeFtChi2 != -1)
            {
             h_MaxTimeFt.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetFit("MaxTime");
             h_MaxTimeFt.at(ich)->SetName(Form("MaxTime%02u", m_channels.at(ich)));
             ChHist.at(ich)->Add(h_MaxTimeFt.at(ich));
            }
         MinTimeFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("MinTime");
         if (MinTimeFtChi2 != -99 && MinTimeFtChi2 != -1)
            {
             h_MinTimeFt.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetFit("MinTime");
             h_MinTimeFt.at(ich)->SetName(Form("MinTime%02u", m_channels.at(ich)));
             ChHist.at(ich)->Add(h_MinTimeFt.at(ich));
            }
         ChargeFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("Charge");
         if (ChargeFtChi2 != -99 && ChargeFtChi2 != -1)
            {
             h_ChargeFt.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetFit("Charge");
             h_ChargeFt.at(ich)->SetName(Form("Charge%02u", m_channels.at(ich)));
             ChHist.at(ich)->Add(h_ChargeFt.at(ich));
            }
         RiseTimeFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("RiseTime");
         if (RiseTimeFtChi2 != -99 && RiseTimeFtChi2 != -1)
            {
             h_RiseTimeFt.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetFit("RiseTime");
             h_RiseTimeFt.at(ich)->SetName(Form("RiseTime%02u", m_channels.at(ich)));
             ChHist.at(ich)->Add(h_RiseTimeFt.at(ich));
            }
         TriggTimeFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("TriggTime");
         if (TriggTimeFtChi2 != -99 && TriggTimeFtChi2 != -1)
            {
             h_TriggTimeFt.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetFit("TriggTime");
             h_TriggTimeFt.at(ich)->SetName(Form("TriggTime%02u", m_channels.at(ich)));
             ChHist.at(ich)->Add(h_TriggTimeFt.at(ich));
            }
         DVDTMaxFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("DVDTMax");
         if (DVDTMaxFtChi2 != -99 && DVDTMaxFtChi2 != -1)
            {
             h_DVDTMaxFt.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetFit("DVDTMax");
             h_DVDTMaxFt.at(ich)->SetName(Form("DVDTMax%02u", m_channels.at(ich)));
             ChHist.at(ich)->Add(h_DVDTMaxFt.at(ich));
            }
         TriggToTFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("TrigToT");
         if (TriggToTFtChi2 != -99 && TriggToTFtChi2 != -1)
            {
             h_TriggToTFt.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetFit("TrigToT");
             h_TriggToTFt.at(ich)->SetName(Form("TrigToT%02u", m_channels.at(ich)));
             ChHist.at(ich)->Add(h_TriggToTFt.at(ich));
            }
         SignalFFTFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("SigFFT");
         if (SignalFFTFtChi2 != -99 && SignalFFTFtChi2 != -1)
            {
             h_SignalFFTFt.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetFit("SigFFT");
             h_SignalFFTFt.at(ich)->SetName(Form("SigFFT%02u", m_channels.at(ich)));
             ChHist.at(ich)->Add(h_SignalFFTFt.at(ich));
            }
         NoiseFFTFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("NoiseFFT");
         if (NoiseFFTFtChi2 != -99 && NoiseFFTFtChi2 != -1)
            {
             h_NoiseFFTFt.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetFit("NoiseFFT");
             h_NoiseFFTFt.at(ich)->SetName(Form("NoiseFFT%02u", m_channels.at(ich)));
             ChHist.at(ich)->Add(h_NoiseFFTFt.at(ich));
            }

         for (unsigned int h = 0; h < 19; h++)
             {
              char title[20];
              CFDTimeFtChi2.at(h) = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("CFDTime", h);
              if (CFDTimeFtChi2.at(h) != -99 && CFDTimeFtChi2.at(h) != -1)
                 {
                  sprintf(title, "CFDTime%02u-%02u%%", m_channels.at(ich), (h*5+5));
                  (h_CFDTimeFt.at(ich)).at(h) = m_RunDUTCh.at(ich)->DUTChannel::GetFit("CFDTime", h);
                  (h_CFDTimeFt.at(ich)).at(h)->SetName((const char*)title);
                  CFDTimeHist.at(ich)->Add((h_CFDTimeFt.at(ich)).at(h));
                 }
              DVDTCFDFtChi2.at(h) = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("DVDTCFD", h);
              if (DVDTCFDFtChi2.at(h) != -99 && DVDTCFDFtChi2.at(h) != -1)
                 {
                  memset(title, '0', sizeof(title));
                  sprintf(title, "DVDTCFD%02u-%02u%%", m_channels.at(ich), (h * 5 + 5));
                  (h_DVDTCFDFt.at(ich)).at(h) = m_RunDUTCh.at(ich)->DUTChannel::GetFit("DVDTCFD", h);
                  (h_DVDTCFDFt.at(ich)).at(h)->SetName((const char*)title);
                  CFDdVdTHist.at(ich)->Add((h_DVDTCFDFt.at(ich)).at(h));
                }
              CFDToTFtChi2.at(h) = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("CFDToT", h);
              if (CFDToTFtChi2.at(h) != -99 && CFDToTFtChi2.at(h) != -1)
                 {
                  memset(title, '0', sizeof(title));
                  sprintf(title, "CFDToT%02u-%02u%%", m_channels.at(ich), (h * 5 + 5));
                  (h_CFDToTFt.at(ich)).at(h) = m_RunDUTCh.at(ich)->DUTChannel::GetFit("CFDToT", h);
                  (h_CFDToTFt.at(ich)).at(h)->SetName((const char*)title);
                  CFDToTHist.at(ich)->Add((h_CFDToTFt.at(ich)).at(h));
                 }
             }
         if (m_RunBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Recovered fits for channel: " << ich + 1 << "/" << m_nchan << std::endl;
         
         unsigned int a = 0;
         unsigned int b = 0;
         // Lets create the Time difference histos and polupate the time difference sigmas histo
         for (int ich2 = ich + 1; ich2 < m_nchan; ich2++)
             {
              a = ((float)(ich*(ich2 - 1)))/2;
              std::pair <double, double> mean;
              std::pair <double, double> sigma;
              int qual;
              char title[20];
              h_TmSigmas.at(a) = new TH2D(Form("h_TimeMap%02u-%02u", m_channels.at(ich), m_channels.at(ich2)), 
                                          Form("CFD Map, Channels %02u-%02u", m_channels.at(ich), m_channels.at(ich2)),
                                          19, 0.025, 0.975, 19, 0.025, 0.975); // initialize the 2D sigmas histos
              for (int bin1 = 0; bin1 < 19; bin1++)
                  {
                   for (int bin2 = 0; bin2 < 19; bin2++)
                       {
                        qual = -99;
                        mean = std::make_pair(-99, -99);
                        sigma = std::make_pair(-99, -99);
                        b = bin1*19+bin2;
                        if (m_RunBase->LGADBase::GetFitMethode() == "rootInt")
                           {
                            qual = m_RunBase->LGADBase::IterativeFit(&(m_EvTmDiff.at(a)).at(b), mean, sigma, 
                                                                    (h_TmDiffCFD.at(a)).at(b), 
                                                                    (TmDiffChi2.at(a)).at(b), 
                                                                    "GaussInt", std::make_pair(-1, -1), false);
                            }
                        else if (m_RunBase->LGADBase::GetFitMethode() == "root" || m_RunBase->LGADBase::GetFitMethode() == "RooFit")
                                {
                                 qual = m_RunBase->LGADBase::IterativeFit(&(m_EvTmDiff.at(a)).at(b), mean, sigma, 
                                                                          (h_TmDiffCFD.at(a)).at(b), 
                                                                          (TmDiffChi2.at(a)).at(b), 
                                                                          "Gauss", std::make_pair(-1, -1), false);
                                }
                        if (qual == 0) mean.second = sigma.first;
                        else {
                              mean.first = m_RunBase->LGADBase::Mean(&(m_EvTmDiff.at(a)).at(b));
                              mean.second = m_RunBase->LGADBase::Stdev(&(m_EvTmDiff.at(a)).at(b));
                              std::cout << __FUNCTION__ << " WARNING: Gaussian fit failed, reverting to mean and std for CFD map between channels " 
                                                        << m_channels.at(ich) << "-" << m_channels.at(ich2) << " and bins " 
                                                        << bin1 * 5 + 5 <<"% - "<< bin2 * 5 + 5 << "%!" << std::endl;
                             }
                        (TmDiffFitQts.at(a)).at(b) = mean;
                        if ((TmDiffChi2.at(a)).at(b) != -99 && (TmDiffChi2.at(a)).at(b) != -1)
                           {
                            memset(title, '0', sizeof(title));
                            sprintf(title, "TmDiffCFD%02u%%-%02u%%CH%02u-%02u", bin1*5+5, bin2*5+5, m_channels.at(ich), m_channels.at(ich2));
                            (h_TmDiffCFD.at(a)).at(b)->SetName((const char*)title);
                            (h_TmDiffCFD.at(a)).at(b)->SetTitle(Form("#Deltat Channel %02u (%02u%% CFD) - Channel %02u (%02u%% CFD)", m_channels.at(ich), bin1*5+5, m_channels.at(ich2), bin2*5+5));
                            TimeDiff.at(a)->Add((h_TmDiffCFD.at(a)).at(b));
                            h_TmSigmas.at(a)->SetBinContent(bin1 + 1, bin2 + 1, ((TmDiffFitQts.at(a)).at(b)).second);
                           }
                        else h_TmSigmas.at(a)->SetBinContent(bin1 + 1, bin2 + 1, 0.0);
                       }
                  }
              histoslist->Add(h_TmSigmas.at(a));
             }

         ChFitParam->Fill();

         // Mean pulseShape properties if shape is calculated
         if (m_WaveShape)
            {
             ChShapePolarity = -99;
             ChShapeMaxIndx = -99;
             ChShapeMax = -99.;
             ChShapeMinIndx = -99;
             ChShapeMin = -99.;
             ChShapeStrIndx = -99;
             ChShapeEndIndx = -99;
             ChShapeMaxTime = -99.;
             ChShapeMinTime = -99.;
             ChShapeNoise = std::make_pair(-99, -99);
             ChShapeCharge = -99.;
             ChShapeRiseTime = -99.;
             ChShapeJitRiseSNR = -99.;
             ChShapeTriggTime = -99.;
             ChShapeDVDTMax = -99.;
             ChShapeTriggToT = -99.;
             ChShapeSignalFFT = -99.;
             ChShapeNoiseFFT = -99.;
             ChShapeSoN = std::make_pair(-99, -99);
             ChShapeCFDTime.clear();
             ChShapeCFDTime.reserve(19);
             ChShapeDVDTCFD.clear();
             ChShapeDVDTCFD.reserve(19);
             ChShapeCFDToT.clear();
             ChShapeCFDToT.reserve(19);
             ChShapeJitNdVdT.clear();
             ChShapeJitNdVdT.reserve(19);
             ChMeanVolt.clear();
             if (m_RunBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Initialized waveshape variables for channel: " << ich + 1 << "/" << m_nchan << std::endl;

             ChShapePolarity = (m_RunDUTCh.at(ich)->DUTChannel::GetChPolarity()).first;
             ChShapeMaxIndx = (m_RunDUTCh.at(ich)->DUTChannel::GetChMaxIndx(1)).first;
             ChShapeMax = (m_RunDUTCh.at(ich)->DUTChannel::GetChMaxVolt(1)).first;
             ChShapeMinIndx = (m_RunDUTCh.at(ich)->DUTChannel::GetChMinIndx(1)).first;
             ChShapeMin = (m_RunDUTCh.at(ich)->DUTChannel::GetChMinVolt(1)).first;
             ChShapeStrIndx = (m_RunDUTCh.at(ich)->DUTChannel::GetChStrIndx(1)).first;
             ChShapeEndIndx = (m_RunDUTCh.at(ich)->DUTChannel::GetChEndIndx(1)).first;
             ChShapeMaxTime = (m_RunDUTCh.at(ich)->DUTChannel::GetChMaxTime(1)).first;
             ChShapeMinTime = (m_RunDUTCh.at(ich)->DUTChannel::GetChMinTime(1)).first;
             ChShapeNoise = m_RunDUTCh.at(ich)->DUTChannel::GetChNoise(1);
             ChShapeCharge = (m_RunDUTCh.at(ich)->DUTChannel::GetChCharge(1)).first;
             ChShapeRiseTime = (m_RunDUTCh.at(ich)->DUTChannel::GetChRiseTime(1)).first;
             ChShapeJitRiseSNR = (m_RunDUTCh.at(ich)->DUTChannel::GetChJitRiseSNR(1)).first;
             ChShapeTriggTime = (m_RunDUTCh.at(ich)->DUTChannel::GetChTriggTime(1)).first;
             ChShapeDVDTMax = (m_RunDUTCh.at(ich)->DUTChannel::GetChDVDTMax(1)).first;
             ChShapeTriggToT = (m_RunDUTCh.at(ich)->DUTChannel::GetChTriggToT(1)).first;
             ChShapeSignalFFT = (m_RunDUTCh.at(ich)->DUTChannel::GetChSignalFFT(1)).first;
             ChShapeNoiseFFT = (m_RunDUTCh.at(ich)->DUTChannel::GetChNoiseFFT(1)).first;
             ChShapeSoN = m_RunDUTCh.at(ich)->DUTChannel::GetChSoNR(1);
             std::vector<std::pair <double, double>> CFDTimeTmp = m_RunDUTCh.at(ich)->DUTChannel::GetChCFDMag(1, "CFDTime");
             std::vector<std::pair <double, double>> DVDTCFDTmp = m_RunDUTCh.at(ich)->DUTChannel::GetChCFDMag(1, "DVDTCFD");
             std::vector<std::pair <double, double>> CFDToTTmp = m_RunDUTCh.at(ich)->DUTChannel::GetChCFDMag(1, "CFDToT");
             std::vector<std::pair <double, double>> JitterNdVdTTmp = m_RunDUTCh.at(ich)->DUTChannel::GetChCFDMag(1, "JitterNdVdT");
             for (unsigned int k = 0; k < CFDTimeTmp.size(); k++) 
                 {
                  ChShapeCFDTime.push_back((CFDTimeTmp.at(k)).first);
                  ChShapeDVDTCFD.push_back((DVDTCFDTmp.at(k)).first);
                  ChShapeCFDToT.push_back((CFDToTTmp.at(k)).first);
                  ChShapeJitNdVdT.push_back((JitterNdVdTTmp.at(k)).first);
                 }
             ChMeanVolt = m_RunDUTCh.at(ich)->DUTChannel::GetMeanChPulse();
             if (m_RunBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Recovered waveshape variables for channel: " 
                                                                      << ich +1 << "/" << m_nchan << std::endl;
             ChShapeParam->Fill();
            }

         if (m_RunBase->LGADBase::GetVerbosity() == 0) m_RunDUTCh.at(ich)->DUTChannel::ChDump();

         // Create folder for each channel and write corresponding histos
         m_ofile->mkdir(Form("ChHistos%02u", m_channels.at(ich)));
         m_ofile->cd(Form("ChHistos%02u", m_channels.at(ich)));
         ChHist.at(ich)->Write();

         m_ofile->mkdir(Form("ChHistos%02u/CFDTimeHists%02u", m_channels.at(ich), m_channels.at(ich)));
         m_ofile->cd(Form("ChHistos%02u/CFDTimeHists%02u", m_channels.at(ich), m_channels.at(ich)));
         CFDTimeHist.at(ich)->Write();
         m_ofile->mkdir(Form("ChHistos%02u/CFDToTHists%02u", m_channels.at(ich), m_channels.at(ich)));
         m_ofile->cd(Form("ChHistos%02u/CFDToTHists%02u", m_channels.at(ich), m_channels.at(ich)));
         CFDToTHist.at(ich)->Write();
         m_ofile->mkdir(Form("ChHistos%02u/CFDdVdTHists%02u", m_channels.at(ich), m_channels.at(ich)));
         m_ofile->cd(Form("ChHistos%02u/CFDdVdTHists%02u", m_channels.at(ich), m_channels.at(ich)));
         CFDdVdTHist.at(ich)->Write();
         for (unsigned int gp = ich + 1; gp < m_nchan; gp++) 
             {
              m_ofile->mkdir(Form("ChHistos%02u/TimeDiff%02u-%02u", m_channels.at(ich), m_channels.at(ich), m_channels.at(gp)));
              m_ofile->cd(Form("ChHistos%02u/TimeDiff%02u-%02u", m_channels.at(ich), m_channels.at(ich), m_channels.at(gp)));
              TimeDiff.at(((float)(ich*(gp-1)))/2)->Write();
             }
         m_ofile->cd();
        }

    ChFitParam->Write();
    histoslist->Write();
    // deleting the histogram lists, deletion happens from end to start
    for (unsigned int a = m_nchan-1; a == 0; a--)
        {
         delete CFDdVdTHist.at(a);
         delete CFDToTHist.at(a);
         delete CFDTimeHist.at(a);
         delete ChHist.at(a);
         for (unsigned int gp = m_nchan - 1; gp == a+1; gp--) delete TimeDiff.at(((float)(a*(gp-1)))/2);
        }
    delete histoslist;
    delete ChFitParam;
    if (m_WaveShape) 
       {
        ChShapeParam->Write();
        delete ChShapeParam;
       }
    delete m_tree;

    m_ofile->Close();
    std::cout << __FUNCTION__ << " INFO: Analysis Finished successfully!" << std::endl;
}
//------------------------------------------------------------------------------
void LGADRun::Terminate()
{
    // for (unsigned int a = 0; a < m_nchan; a++) delete m_RunWaveform.at(a);
    std::cout << setfill('-') << setw(89) << "-" << std::endl << setfill(' ');
    //std::cout << __FUNCTION__ << " INFO: Cleaning up objeects..." << std::endl;
    //for (unsigned int a = 0; a < m_nchan; a++) delete m_RunDUTCh.at(a);
}
//------------------------------------------------------------------------------
void LGADRun::FirstInit()
{
  // Create new DUT channels
  m_RunDUTCh.resize(m_nchan, nullptr);
  std::vector<bool> foundPr[8];
  for (unsigned int h = 0; h < 8; h++) foundPr[h].resize(m_nchan, false);
  int sizes[8];
  sizes[0] = (m_RunBase->LGADBase::GetDUTNames()).size();
  sizes[1] = (m_RunBase->LGADBase::GetDUTBoards()).size();
  sizes[2] = (m_RunBase->LGADBase::GetDUTransImps()).size();
  sizes[3] = (m_RunBase->LGADBase::GetDUTSecStages()).size();
  sizes[4] = (m_RunBase->LGADBase::GetDUTAmpGains()).size();
  sizes[5] = (m_RunBase->LGADBase::GetDUTriggs()).size();
  sizes[6] = (m_RunBase->LGADBase::GetDUTFracts()).size();
  sizes[7] = (m_RunBase->LGADBase::GetDUTCaps()).size();
  for (unsigned int a = 0; a < m_nchan; a++)
      {
       m_RunDUTCh.at(a) = new DUTChannel(m_channels.at(a), m_RunBase);
       for (unsigned int u = 0; u < 8; u++)
           {
            for (unsigned int w = 0; w < sizes[u]; w++) 
                {
                 if (u == 0)
                    {
                     if (m_channels.at(a) == ((m_RunBase->LGADBase::GetDUTNames()).at(w)).first) 
                        {
                         m_RunDUTCh.at(a)->DUTChannel::SetChName(((m_RunBase->LGADBase::GetDUTNames()).at(w)).second);
                         (foundPr[u]).at(a) = true;
                        }
                    }
                 if (u == 1)
                    {
                     if (m_channels.at(a) == ((m_RunBase->LGADBase::GetDUTBoards()).at(w)).first)
                        {
                         m_RunDUTCh.at(a)->DUTChannel::SetChBoard(((m_RunBase->LGADBase::GetDUTBoards()).at(w)).second);
                         (foundPr[u]).at(a) = true;
                        }
                    }
                 if (u == 2)
                    {
                     if (m_channels.at(a) == ((m_RunBase->LGADBase::GetDUTransImps()).at(w)).first)
                        {
                         m_RunDUTCh.at(a)->DUTChannel::SetChTransImp(((m_RunBase->LGADBase::GetDUTransImps()).at(w)).second);
                         (foundPr[u]).at(a) = true;
                        }
                    }
                 if (u == 3)
                    {
                     if (m_channels.at(a) == ((m_RunBase->LGADBase::GetDUTSecStages()).at(w)).first)
                        {
                         m_RunDUTCh.at(a)->DUTChannel::SetSecStage(((m_RunBase->LGADBase::GetDUTSecStages()).at(w)).second);
                         (foundPr[u]).at(a) = true;
                        }
                    }
                 if (u == 4)
                    {
                     if (m_channels.at(a) == ((m_RunBase->LGADBase::GetDUTAmpGains()).at(w)).first)
                        {
                         m_RunDUTCh.at(a)->DUTChannel::SetChAmpGain(((m_RunBase->LGADBase::GetDUTAmpGains()).at(w)).second);
                         (foundPr[u]).at(a) = true;
                        }
                    }
                 if (u == 5)
                    {
                     if (m_channels.at(a) == ((m_RunBase->LGADBase::GetDUTriggs()).at(w)).first)
                        {
                         m_RunDUTCh.at(a)->DUTChannel::SetChTrigg(((m_RunBase->LGADBase::GetDUTriggs()).at(w)).second);
                         (foundPr[u]).at(a) = true;
                        }
                    }
                 if (u == 6)
                    {
                     if (m_channels.at(a) == ((m_RunBase->LGADBase::GetDUTFracts()).at(w)).first)
                        {
                         m_RunDUTCh.at(a)->DUTChannel::SetChFract(((m_RunBase->LGADBase::GetDUTFracts()).at(w)).second);
                         (foundPr[u]).at(a) = true;
                        }
                    }
                 if (u == 7)
                    {
                     if (m_channels.at(a) == ((m_RunBase->LGADBase::GetDUTCaps()).at(w)).first)
                        {
                         m_RunDUTCh.at(a)->DUTChannel::SetChCap(((m_RunBase->LGADBase::GetDUTCaps()).at(w)).second);
                         (foundPr[u]).at(a) = true;
                        }
                    }
                } 
            }
      }

  if (m_RunBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Intialized Channel ID properties for " << m_nchan << " channels" << std::endl;

  std::vector<unsigned int> skip[8];
  // Fill the properties of the channel without index
  for (unsigned int a = 0; a < m_nchan; a++)
      {
       for (unsigned int u = 0; u < 8; u++)
           {
            if (!(foundPr[u]).at(a))
               {
                for (unsigned int w = 0; w < sizes[u]; w++)
                    {
                     for (unsigned int tc = 0; tc < (skip[u]).size(); tc++) if (w == (skip[u]).at(tc)) continue;
                     if (u == 0 && ((m_RunBase->LGADBase::GetDUTNames()).at(w)).first == -99)
                        {
                         m_RunDUTCh.at(a)->DUTChannel::SetChName(((m_RunBase->LGADBase::GetDUTNames()).at(w)).second);
                         (skip[u]).push_back(w);
                        }
                     if (u == 1 && ((m_RunBase->LGADBase::GetDUTBoards()).at(w)).first == -99)
                        {
                         m_RunDUTCh.at(a)->DUTChannel::SetChBoard(((m_RunBase->LGADBase::GetDUTBoards()).at(w)).second);
                         (skip[u]).push_back(w);
                        }
                     if (u == 2 && ((m_RunBase->LGADBase::GetDUTransImps()).at(w)).first == -99)
                        {
                         m_RunDUTCh.at(a)->DUTChannel::SetChTransImp(((m_RunBase->LGADBase::GetDUTransImps()).at(w)).second);
                         (skip[u]).push_back(w);
                        }
                     if (u == 3 && ((m_RunBase->LGADBase::GetDUTSecStages()).at(w)).first == -99)
                        {
                         m_RunDUTCh.at(a)->DUTChannel::SetSecStage(((m_RunBase->LGADBase::GetDUTSecStages()).at(w)).second);
                         (skip[u]).push_back(w);
                        }
                     if (u == 4 && ((m_RunBase->LGADBase::GetDUTAmpGains()).at(w)).first == -99)
                        {
                         m_RunDUTCh.at(a)->DUTChannel::SetChAmpGain(((m_RunBase->LGADBase::GetDUTAmpGains()).at(w)).second);
                         (skip[u]).push_back(w);
                        }
                     if (u == 5 && ((m_RunBase->LGADBase::GetDUTriggs()).at(w)).first == -99)
                        {
                         m_RunDUTCh.at(a)->DUTChannel::SetChTrigg(((m_RunBase->LGADBase::GetDUTriggs()).at(w)).second);
                         (skip[u]).push_back(w);
                        }
                     if (u == 6 && ((m_RunBase->LGADBase::GetDUTFracts()).at(w)).first == -99)
                        {
                         m_RunDUTCh.at(a)->DUTChannel::SetChFract(((m_RunBase->LGADBase::GetDUTFracts()).at(w)).second);
                         (skip[u]).push_back(w);
                        }
                     if (u == 7 && ((m_RunBase->LGADBase::GetDUTCaps()).at(w)).first == -99)
                        {
                         m_RunDUTCh.at(a)->DUTChannel::SetChCap(((m_RunBase->LGADBase::GetDUTCaps()).at(w)).second);
                         (skip[u]).push_back(w);
                        }
                    }
               }
           }
      }

  VarInit();

  // Create new Branches
  b_EvPol.clear();
  b_EvCharge.clear();
  b_EvIsSignal.clear();
  b_EvInWindow.clear();
  b_EvVoltStr.clear();
  b_EvMaxindx.clear();
  b_EvMaxTime.clear();
  b_EvMaxVolt.clear();
  b_EvMinIndx.clear();
  b_EvMinTime.clear();
  b_EvMinVolt.clear();
  b_EvStrIndx.clear();
  b_EvEndIndx.clear();
  b_EvNoise.clear();
  b_EvPedest.clear();
  b_EvNoiseErr.clear();
  b_EvPedestErr.clear();
  b_EvRiseT.clear();
  b_EvJitNdVdT.clear();
  b_EvJitRiseSNR.clear();
  b_EvTriggTime.clear();
  b_EvDVDTMax.clear();
  b_EvDVDTCFD.clear();
  b_EvCFDToT.clear();
  b_EvTriggToT.clear();
  b_EvSignalFFT.clear();
  b_EvNoiseFFT.clear();
  b_EvComplete.clear();
  b_EvCFDTime.clear();
  b_EvVAdjCFD.clear();

  b_EvPol.resize(m_nchan);
  b_EvCharge.resize(m_nchan);
  b_EvIsSignal.resize(m_nchan);
  b_EvInWindow.resize(m_nchan);
  b_EvVoltStr.resize(m_nchan);
  b_EvMaxindx.resize(m_nchan);
  b_EvMaxTime.resize(m_nchan);
  b_EvMaxVolt.resize(m_nchan);
  b_EvMinIndx.resize(m_nchan);
  b_EvMinTime.resize(m_nchan);
  b_EvMinVolt.resize(m_nchan);
  b_EvStrIndx.resize(m_nchan);
  b_EvEndIndx.resize(m_nchan);
  b_EvNoise.resize(m_nchan);
  b_EvPedest.resize(m_nchan);
  b_EvNoiseErr.resize(m_nchan);
  b_EvPedestErr.resize(m_nchan);
  b_EvRiseT.resize(m_nchan);
  b_EvJitNdVdT.resize(m_nchan);
  b_EvJitRiseSNR.resize(m_nchan);
  b_EvTriggTime.resize(m_nchan);
  b_EvDVDTMax.resize(m_nchan);
  b_EvDVDTCFD.resize(m_nchan);
  b_EvCFDToT.resize(m_nchan);
  b_EvTriggToT.resize(m_nchan);
  b_EvSignalFFT.resize(m_nchan);
  b_EvNoiseFFT.resize(m_nchan);
  b_EvComplete.resize(m_nchan);
  b_EvCFDTime.resize(m_nchan);
  b_EvVAdjCFD.resize(m_nchan);

  // Fits Initialization
  h_MaxIndxFt.clear();
  h_MinIndxFt.clear();
  h_MaxVoltFt.clear();
  h_MinVoltFt.clear();
  h_StrIndxFt.clear();
  h_EndIndxFt.clear();
  h_NoiseFt.clear();
  h_NoiseErrFt.clear();
  h_PedestFt.clear();
  h_NPedestErrFt.clear();
  h_MaxTimeFt.clear();
  h_MinTimeFt.clear();
  h_ChargeFt.clear();
  h_RiseTimeFt.clear();
  h_TriggTimeFt.clear();
  h_DVDTMaxFt.clear();
  h_DVDTCFDFt.clear();
  h_CFDToTFt.clear();
  h_TriggToTFt.clear();
  h_SignalFFTFt.clear();
  h_NoiseFFTFt.clear();
  h_CFDTimeFt.clear();

  h_MaxIndxFt.resize(m_nchan, nullptr);
  h_MinIndxFt.resize(m_nchan, nullptr);
  h_MaxVoltFt.resize(m_nchan, nullptr);
  h_MinVoltFt.resize(m_nchan, nullptr);
  h_StrIndxFt.resize(m_nchan, nullptr);
  h_EndIndxFt.resize(m_nchan, nullptr);
  h_NoiseFt.resize(m_nchan, nullptr);
  h_NoiseErrFt.resize(m_nchan, nullptr);
  h_PedestFt.resize(m_nchan, nullptr);
  h_NPedestErrFt.resize(m_nchan, nullptr);
  h_MaxTimeFt.resize(m_nchan, nullptr);
  h_MinTimeFt.resize(m_nchan, nullptr);
  h_ChargeFt.resize(m_nchan, nullptr);
  h_RiseTimeFt.resize(m_nchan, nullptr);
  h_TriggTimeFt.resize(m_nchan, nullptr);
  h_DVDTMaxFt.resize(m_nchan, nullptr);
  h_TriggToTFt.resize(m_nchan, nullptr);
  h_SignalFFTFt.resize(m_nchan, nullptr);
  h_NoiseFFTFt.resize(m_nchan, nullptr);
  h_CFDTimeFt.resize(m_nchan);
  h_DVDTCFDFt.resize(m_nchan);
  h_CFDToTFt.resize(m_nchan);
  for (unsigned int s = 0; s < m_nchan; s++)
      {
       (h_CFDTimeFt.at(s)).clear();
       (h_CFDTimeFt.at(s)).resize(19, nullptr);
       (h_DVDTCFDFt.at(s)).clear();
       (h_DVDTCFDFt.at(s)).resize(19, nullptr);
       (h_CFDToTFt.at(s)).clear();
       (h_CFDToTFt.at(s)).resize(19, nullptr);
      }

  m_ofile->ReOpen("UPDATE");
  m_ofile->SetCompressionLevel(6);
  // m_ofile->GetObject(m_tree->GetName(), m_tree); // Should not need that but nont sure I am refering to the same object

  // Connect the branches to the vectors
  for (unsigned int a = 0; a < m_nchan; a++)
      {
       b_EvPol.at(a) = m_tree->Branch(Form("Polarity%02u", m_channels.at(a)), &(m_EvPol.at(a)), "m_EvPol.at(a)/I");
       b_EvCharge.at(a) = m_tree->Branch(Form("Charge%02u", m_channels.at(a)), &(m_EvCharge.at(a)), "m_EvCharge.at(a)/D");
       b_EvIsSignal.at(a) = m_tree->Branch(Form("IsSignal%02u", m_channels.at(a)), &(m_EvIsSignal.at(a)), "m_EvIsSignal.at(a)/I");
       b_EvInWindow.at(a) = m_tree->Branch(Form("InWindow%02u", m_channels.at(a)), &(m_EvInWindow.at(a)), "m_EvInWindow.at(a)/I");
       b_EvVoltStr.at(a) = m_tree->Branch(Form("VoltStr%02u", m_channels.at(a)), &(m_EvVoltStr.at(a)), "m_EvVoltStr.at(a)/I");
       b_EvMaxindx.at(a) = m_tree->Branch(Form("MaxIndx%02u", m_channels.at(a)), &(m_EvMaxindx.at(a)), "m_EvMaxindx.at(a)/I");
       b_EvMaxTime.at(a) = m_tree->Branch(Form("MaxTime%02u", m_channels.at(a)), &(m_EvMaxTime.at(a)), "m_EvMaxTime.at(a)/D");
       b_EvMaxVolt.at(a) = m_tree->Branch(Form("MaxVolt%02u", m_channels.at(a)), &(m_EvMaxVolt.at(a)), "m_EvMaxVolt.at(a)/D");
       b_EvMinIndx.at(a) = m_tree->Branch(Form("MinIndx%02u", m_channels.at(a)), &(m_EvMinIndx.at(a)), "m_EvMinIndx.at(a)/I");
       b_EvMinTime.at(a) = m_tree->Branch(Form("MinTime%02u", m_channels.at(a)), &(m_EvMinTime.at(a)), "m_MinTime.at(a)/D");
       b_EvMinVolt.at(a) = m_tree->Branch(Form("MinVolt%02u", m_channels.at(a)), &(m_EvMinVolt.at(a)), "m_EvMinVolt.at(a)/D");
       b_EvStrIndx.at(a) = m_tree->Branch(Form("Strindx%02u", m_channels.at(a)), &(m_EvStrIndx.at(a)), "m_EvStrIndx.at(a)/I");
       b_EvEndIndx.at(a) = m_tree->Branch(Form("EndIndx%02u", m_channels.at(a)), &(m_EvEndIndx.at(a)), "m_EvEndIndx.at(a)/I");
       b_EvNoise.at(a) = m_tree->Branch(Form("Noise%02u", m_channels.at(a)), &(m_EvNoise.at(a)), "m_EvNoise.at(a)/D");
       b_EvPedest.at(a) = m_tree->Branch(Form("Pedest%02u", m_channels.at(a)), &(m_EvPedest.at(a)), "m_EvPedest.at(a)/D");
       b_EvNoiseErr.at(a) = m_tree->Branch(Form("NoiseErr%02u", m_channels.at(a)), &(m_EvNoiseErr.at(a)), "m_EvNoiseErr.at(a)/D");
       b_EvPedestErr.at(a) = m_tree->Branch(Form("PedestErr%02u", m_channels.at(a)), &(m_EvPedestErr.at(a)), "m_EvPedestEr.at(a)/D");
       b_EvRiseT.at(a) = m_tree->Branch(Form("RiseT%02u", m_channels.at(a)), &(m_EvRiseT.at(a)), "m_EvRiseTe.at(a)/D");
       b_EvJitRiseSNR.at(a) = m_tree->Branch(Form("JitRiseSNR%02u", m_channels.at(a)), &(m_EvJitRiseSNR.at(a)), "m_EvJitRiseSNR.at(a)/D");
       b_EvTriggTime.at(a) = m_tree->Branch(Form("TriggTime%02u", m_channels.at(a)), &(m_EvTriggTime.at(a)), "m_EvTriggTime.at(a)/D");
       b_EvDVDTMax.at(a) = m_tree->Branch(Form("DVDTMax%02u", m_channels.at(a)), &(m_EvDVDTMax.at(a)), "m_EvDVDTMax.at(a)/D");
       b_EvTriggToT.at(a) = m_tree->Branch(Form("TriggToT%02u", m_channels.at(a)), &(m_EvTriggToT.at(a)), "m_EvTriggToT.at(a)/D");
       b_EvSignalFFT.at(a) = m_tree->Branch(Form("SignalFFT%02u", m_channels.at(a)), &(m_EvSignalFFT.at(a)), "m_EvSignalFFT.at(a)/D");
       b_EvNoiseFFT.at(a) = m_tree->Branch(Form("NoiseFFT%02u", m_channels.at(a)), &(m_EvNoiseFFT.at(a)), "m_EvNoiseFFT.at(a)/D");
       b_EvComplete.at(a) = m_tree->Branch(Form("Complete%02u", m_channels.at(a)), &(m_EvComplete.at(a)), "m_EvComplete.at(a)/I");
       b_EvVAdjCFD.at(a) = m_tree->Branch(Form("VAdjCFD%02u", m_channels.at(a)), &(m_EvVAdjCFD.at(a)));
       b_EvCFDTime.at(a) = m_tree->Branch(Form("CFDTime%02u", m_channels.at(a)), &(m_EvCFDTime.at(a)));
       b_EvDVDTCFD.at(a) = m_tree->Branch(Form("DVDTCFD%02u", m_channels.at(a)), &(m_EvDVDTCFD.at(a)));
       b_EvCFDToT.at(a) = m_tree->Branch(Form("CFDToT%02u", m_channels.at(a)), &(m_EvCFDToT.at(a)));
       b_EvJitNdVdT.at(a) = m_tree->Branch(Form("JitNdVdT%02u", m_channels.at(a)), &(m_EvJitNdVdT.at(a)));
      }

  m_init++;
}
//------------------------------------------------------------------------------
void LGADRun::VarInit()
{
  // Initiakize event variables
  m_EvPol.clear();
  m_EvCharge.clear();
  m_EvIsSignal.clear();
  m_EvInWindow.clear();
  m_EvVoltStr.clear();
  m_EvMaxindx.clear();
  m_EvMaxTime.clear();
  m_EvMaxVolt.clear();
  m_EvMinIndx.clear();
  m_EvMinTime.clear();
  m_EvMinVolt.clear();
  m_EvStrIndx.clear();
  m_EvEndIndx.clear();
  m_EvNoise.clear();
  m_EvPedest.clear();
  m_EvNoiseErr.clear();
  m_EvPedestErr.clear();
  m_EvRiseT.clear();
  m_EvJitNdVdT.clear();
  m_EvJitRiseSNR.clear();
  m_EvTriggTime.clear();
  m_EvDVDTMax.clear();
  m_EvTriggToT.clear();
  m_EvSignalFFT.clear();
  m_EvNoiseFFT.clear();
  m_EvComplete.clear();
  m_EvCFDTime.clear();
  m_EvDVDTCFD.clear();
  m_EvCFDToT.clear();
  m_EvVAdjCFD.clear();
  m_dTMax.clear();
  m_dCMax.clear();

  // resize vectors to proper size to avoid memory fragmentation
  m_EvPol.resize(m_nchan, -99);
  m_EvCharge.resize(m_nchan, -99);
  m_EvIsSignal.resize(m_nchan, -99);
  m_EvInWindow.resize(m_nchan, -99);
  m_EvVoltStr.resize(m_nchan, -99);
  m_EvMaxindx.resize(m_nchan, -99);
  m_EvMaxTime.resize(m_nchan, -99);
  m_EvMaxVolt.resize(m_nchan, -99);
  m_EvMinIndx.resize(m_nchan, -99);
  m_EvMinTime.resize(m_nchan, -99);
  m_EvMinVolt.resize(m_nchan, -99);
  m_EvStrIndx.resize(m_nchan, -99);
  m_EvEndIndx.resize(m_nchan, -99);
  m_EvNoise.resize(m_nchan, -99);
  m_EvPedest.resize(m_nchan, -99);
  m_EvNoiseErr.resize(m_nchan, -99);
  m_EvPedestErr.resize(m_nchan, -99);
  m_EvRiseT.resize(m_nchan, -99);
  m_EvJitRiseSNR.resize(m_nchan, -99);
  m_EvTriggTime.resize(m_nchan, -99);
  m_EvDVDTMax.resize(m_nchan, -99);
  m_EvTriggToT.resize(m_nchan, -99);
  m_EvSignalFFT.resize(m_nchan, -99);
  m_EvNoiseFFT.resize(m_nchan, -99);
  m_EvComplete.resize(m_nchan, -99);
  m_EvCFDTime.resize(m_nchan);
  m_EvJitNdVdT.resize(m_nchan);
  m_EvDVDTCFD.resize(m_nchan);
  m_EvCFDToT.resize(m_nchan);
  for (unsigned int s = 0; s < m_nchan; s++)
      {
       (m_EvCFDTime.at(s)).clear();
       (m_EvCFDTime.at(s)).resize(19, -99);
       (m_EvJitNdVdT.at(s)).clear();
       (m_EvJitNdVdT.at(s)).resize(19, -99);
       (m_EvDVDTCFD.at(s)).clear();
       (m_EvDVDTCFD.at(s)).resize(19, -99);
       (m_EvCFDToT.at(s)).clear();
       (m_EvCFDToT.at(s)).resize(19, -99);
      }
  m_EvVAdjCFD.resize(m_nchan);

  int combi = (m_nchan*(m_nchan-1))/2;
  m_EvTmDiff.resize(combi);
  h_TmDiffCFD.resize(combi);
  h_TmSigmas.resize(combi);
  for (unsigned int k = 0; k < combi; k++)
      {
       (m_EvTmDiff.at(k)).resize(19*19);
       (h_TmDiffCFD.at(k)).resize(19*19, nullptr);
       h_TmSigmas.push_back(nullptr);
      }
}
//------------------------------------------------------------------------------
