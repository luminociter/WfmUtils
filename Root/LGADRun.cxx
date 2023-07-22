/*
*                 WFMRun.cxx
*
*
*      Author: Gkougkousis Evangelos - Leonidas
*             egkougko@cern.ch - 2023
*                    CERN-GENEVA
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

LGADRun::LGADRun()
{
    m_RunBase = new LGADBase();
    std::cout << setfill('-') << setw(89) << "-" << setfill(' ') << std::endl;
}
// --------------------------------------------------------------------------------------------------------------
LGADRun::LGADRun(LGADBase* rBase)
{
    LGADBase::CopyProperties(*rBase); // setting the basic properties LGADRun object
    m_RunBase = rBase; // createing a class-wide pointer of LGADBase connected to the prototype
    if (m_RunBase->LGADBase::GetVerbosity() >= 2) std::cout << "LGAD Run constructor called and completed!" << std::endl;
    std::cout << setfill('-') << setw(89) << "-" << setfill(' ') << std::endl;
}
// --------------------------------------------------------------------------------------------------------------
LGADRun::~LGADRun()
{
    delete m_RunBase;
    for (unsigned int gh = 0; gh < m_RunDUTCh.size(); gh++) delete m_RunDUTCh.at(gh);
}
// --------------------------------------------------------------------------------------------------------------
bool LGADRun::InitAnalysis(std::vector<unsigned int>& fileindx)
{
    // Function to be run at the beginning to iniialise all the pointers and estalish analysis stage and options of input files
    LGADRun::InitNtuple(m_RunBase->LGADBase::GetRootTrees(), m_RunBase->LGADBase::GetInRootFiles());
    m_event = 0;
    // here we need to create the output file and tree
    LGADRun::InitOutputBranches();

    if (m_RunBase->LGADBase::GetTrackComb())
       {
        // Exclude the files that are set to not be processed with tracking
        std::vector<unsigned int> l_fileindx_trck = fileindx;
        if (m_RunBase->LGADBase::GetExcludeTrackFiles()).size() > 0)
           {
            for (unsigned int q = 0; q < l_fileindx_trck.size(); q++
                {
                 for (unsigned int p = 0; p < m_RunBase->LGADBase::GetExcludeTrackFiles()).size(); p++)
                     {
                      if (l_fileindx_trck.at(q) == m_RunBase->LGADBase::GetExcludeTrackFiles().at(p))
                         {
                          l_fileindx_trck.erase(l_fileindx_trck.begin() + q)
                          q--;
                          break;
                         }
                     }
                }
           }
        // Loading tracking fuctions, just openening at this stage
        if (!OpenTrackFile(m_Trackdatadir, m_Trackdataname, l_fileindx_trck))
           {
            std::cout << __FUNCTION__ << " ERROR: Failed oppening tracking files in directory: " << m_Trackdatadir << std::endl;
            return false;
           }
        else {
              // here we need to initialise the tracking branches using the ones from both the inpus files and the tracking files
              if (!AddTrckBranches(, m_RunBase->LGADBase::GetTrackPackage()))
                 {
                  std::cout << __FUNCTION__ << " Warning: Failed adding the tracking part on the output file: " << m_ofname << std::endl;
                  std::cout << "                                       Switching to non-trackimg mode..." << std::endl;
                  LGADBase::SetTrackComb(false);
                 }
             }
       }



    // Start looping on all input files. The vetor of input files has the ssame size as the vector of fileindx. That does ont necesseraly mean that 
    // all root files in the directory were succesfully opeened!
    for (unsigned int a = 0; a < fileindx.size(); a++)
        {
         std::cout << "--> File Accessed: " << m_RunBase->LGADBase::GetInRootFile(a)->GetName() << std::endl;
         unsigned int evnt = m_RunBase->LGADBase::GetRootTree(indx)->GetEntriesFast();
         std::cout << "--> Total number of entries in file: " << evnt << std::endl;
         m_event += evnt;

         m_RunBase->LGADBase::GetInRootFile(a);
        }




    bool SetBranchAccess(TTree * inputree, unsigned int FFTShape, unsigned int filelvl, std::vector<unsigned int> filechannel,
                         std::vector<unsigned int> fileplanes, std::vector<unsigned int> fileRefs, std::vector<unsigned int> fileDUTs,
                         TTree * ChTree, TTree * ShapeTree, TTree * TimeResTree);
    bool GetEvntBranches(Long64_t entry, unsigned int FFTShape, unsigned int filelvl, std::vector<unsigned int> filechannel, unsigned int fileindx,
                         std::vector<unsigned int> fileplanes, std::vector<unsigned int> fileRefs, std::vector<unsigned int> fileDUTs);
    bool GetChanelBranches(Long64_t channel, unsigned int fileindx, unsigned int stage, unsigned int FFTShape, unsigned int filelvl);
    bool SetOutputBranchAdresses(unsigned int FFTShape, unsigned int filelvl, std::vector<unsigned int> filechannel, unsigned int fileindx,
                                 std::vector<unsigned int> fileplanes, std::vector<unsigned int> fileRefs, std::vector<unsigned int> fileDUTs);
    bool AddOutputHistos();
    bool WriteBranches(unsigned int fileindx);







    printf("Client slave begin!\n");
    m_init = 0;
    m_nplanes.clear();
    m_nRef.clear();
    m_nDUT.clear();
    m_DUTs.clear();
    m_Refs.clear();
    m_planes.clear();
    std::cout << setfill('-') << setw(89) << "-" << std::endl << setfill(' ');
}
// --------------------------------------------------------------------------------------------------------------
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
      std::cout << "--> File Accessed: " << m_RunBase->LGADBase::GetInRootFile()->GetName() << std::endl;
      if (m_event == 0) m_event = fChain->GetTree()->GetEntriesFast();
     }
  else {                                   // We are running on PROOF:
        std::cout << "--> Running on PROOF" << std::endl;
        m_RunBase->LGADBase::SetInRootFile(fChain->GetCurrentFile());
        std::cout << "--> File Accessed: " << fChain->GetName() << std::endl;
        if (m_event == 0) m_event = fChain->GetTree()->GetEntriesFast();
       }
  std::cout << "--> Total number of entries in file: " << m_event << std::endl;

  // should only be called once, not every time a new file is opened
  // Needs to be called AFTER selector has been initialised (not before Notify)
  if (m_init > 0) 
     {
      m_RunBase->LGADBase::SetInstrument(m_instrument);
      m_RunBase->LGADBase::SetChannels(LGADBase::GetChannels());
      m_RunBase->LGADBase::SetFileLVL(LGADBase::GetFileLVL());
      m_RunBase->LGADBase::SetTrackPackage(LGADBase::GetTrackPackage());
      m_RunBase->LGADBase::SetTrackComb(LGADBase::GetTrackComb());
      m_RunBase->LGADBase::SetFEi4Eff(LGADBase::GetFEi4Eff());
      FirstInit(m_RunBase->LGADBase::GetFileLVL(), m_RunBase->LGADBase::GetAnaStage());
      if ((m_RunBase->LGADBase::GetPlaneDTs()).size() > 0) SetPlaneDQuants(m_dTMax, m_RunBase->GetPlaneDTs());
      if ((m_RunBase->LGADBase::GetPlaneDCs()).size() > 0) SetPlaneDQuants(m_dCMax, m_RunBase->GetPlaneDCs());
     }
  m_init++;
  std::cout << setfill('-') << setw(89) << "-" << std::endl << setfill(' ');
  return kTRUE;
 }
// --------------------------------------------------------------------------------------------------------------
Bool_t LGADRun::Process(Long64_t entry)
{
    if (fChain == 0) return kFALSE;
    b_EvnNo->GetEntry(entry);
    // Do Event tests for debugging
    if (m_RunBase->LGADBase::GetTestEvn(EvnNo)) return kTRUE;

    if (entry == 0)
       { 
        m_evnt1 = m_RunBase->LGADBase::GetStartStopEvnt().first;
        m_evnt2 = m_RunBase->LGADBase::GetStartStopEvnt().second;
        if (m_evnt2 == 0 || m_evnt2 <= m_evnt1 || m_evnt2 > (m_event + EvnNo)) 
           {
            if (m_RunBase->LGADBase::GetInstrument() == Sampic) m_evnt2 = (m_event + EvnNo)*m_nchan;
            else m_evnt2 = m_event + EvnNo;
           }
        else if (m_evnt2 <= EvnNo && (m_evnt2 - m_evnt1) <= m_event && m_evnt1 < m_evnt2) { m_evnt2 += EvnNo; m_evnt1 += EvnNo; }
        else if (m_evnt2 > EvnNo && m_evnt1 < EvnNo && (m_evnt2 - EvnNo) < m_event && m_evnt1 < m_evnt2) m_evnt1 += EvnNo;
       }
    else {
          // Stage one, always has to happen if we run this code anyway
          if (m_RunBase->LGADBase::GetFileLVL() < 1 || m_RunBase->LGADBase::GetFileLVL() > 4) VarInit();
          // Stage three analysis only
          if (m_RunBase->LGADBase::GetAnaStage() > 2 && m_RunBase->LGADBase::GetFileLVL() < 3) ThirdStageInit();
         }

    if (m_evnt2 != 0 && m_evnt2 > m_evnt1 && EvnNo > m_evnt2) return false;
    LGADBase::ProgressBar(entry, m_evnt2-m_evnt1+1);

    if (m_RunBase->LGADBase::GetVerbosity() >= 2) 
       {
        std::cout << "-------------------------------------------------------------------" << std::endl;
        std::cout << __FUNCTION__ << " INFO: Number of active channels: " << m_nchan << ", instrument: " << m_instrument << ", event: " << entry << std::endl;
        std::cout << __FUNCTION__ << " INFO: Channels: ";
        for (unsigned int k = 0; k < m_nchan; k++) std::cout << m_channels.at(k) << " ";
        std::cout << std::endl;
       }

    GetBranches((Long64_t entry, m_base);
        
    for (unsigned int ich = 0; ich < m_nchan; ich++)
        {
         //m_RunWaveform.at(ich)->WaveForm::InitializeWaveForm(1);
         //m_RunWaveform.at(ich)->WaveForm::SetVoltage(w.at(ich));
         //m_RunWaveform.at(ich)->WaveForm::SetTime(t.at(ich));
         //m_RunWaveform.at(ich)->WaveForm::Calculate();
         if (m_RunBase->LGADBase::GetInstrument() == Sampic) m_RunDUTCh.at(ich)->DUTChannel::SetChRate(m_srate.at(0));
         else m_RunDUTCh.at(ich)->DUTChannel::SetChRate(m_srate.at(ich));
         m_RunDUTCh.at(ich)->DUTChannel::AppendEvent(&m_w.at(ich), &m_t.at(ich));
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
         if (m_RunBase->LGADBase::GetDoFFT())
            {
             m_EvSignalFFT.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetEvntSignalFFT();
             m_EvNoiseFFT.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetEvntNoiseFFT();
            }
         m_EvComplete.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetEvntComplete();
         m_EvVAdjCFD.at(ich) = m_RunDUTCh.at(ich)->DUTChannel::GetEvntVAdjCFD();
         if (m_RunBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Got event variables: " << ich + 1 << "/" << m_nchan << "!" << std::endl;

         for (unsigned int bin = 0; bin < 19; bin++)
             {
              (m_EvCFDTime.at(ich)).at(bin) = m_RunDUTCh.at(ich)->DUTChannel::GetEvntCFDTime(bin);
              (m_EvDVDTCFD.at(ich)).at(bin) = m_RunDUTCh.at(ich)->DUTChannel::GetEvntDVDTCFD(bin);
              (m_EvCFDToT.at(ich)).at(bin) = m_RunDUTCh.at(ich)->DUTChannel::GetEvntCFDToT(bin);
              (m_EvJitNdVdT.at(ich)).at(bin) = m_RunDUTCh.at(ich)->DUTChannel::GetEvntJitNdVdT(bin);
             }
         if (m_RunBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Got CFD points: " << ich + 1 << "/" << m_nchan << "!" << std::endl;
        }
    if (m_RunBase->LGADBase::GetAnaStage() > 1) ThirdStageProcess();
    WriteBranches(m_base);

    return kTRUE;
}
// --------------------------------------------------------------------------------------------------------------
void LGADRun::SlaveTerminate()
{
    std::cout << setfill('-') << setw(89) << "-" << std::endl << setfill(' '); 
    m_OutTree->Write("", TObject::kOverwrite);
    TTree* ChFitParam = new TTree("ChFitParam", "Channel parameters form fits");
    if (m_RunBase->LGADBase::GetAnaStage() > 1) SecondStageAnalysis(ChFitParam); // Channel Properties
    if (m_RunBase->LGADBase::GetAnaStage() > 2) ThirdStageAnalysis(ChFitParam);  // CFD Maps and chi2 minimization
    if (m_RunBase->LGADBase::GetAnaStage() > 3) TrckEffAnalaysis(); // Tracking efficiency
    delete ChFitParam;
    delete m_OutTree;

    if (m_RunBase->LGADBase::GetInRootFile() != m_RunBase->LGADBase::GetOutRootFile()) 
       {
        m_RunBase->LGADBase::GetOutRootFile()->Close();
        m_RunBase->LGADBase::GetInRootFile()->cd();
       }
    m_RunBase->LGADBase::GetInRootFile()->Close();
    std::cout << __FUNCTION__ << " INFO: Analysis Finished successfully!" << std::endl;
}
// --------------------------------------------------------------------------------------------------------------
void LGADRun::Terminate()
{
    // for (unsigned int a = 0; a < m_nchan; a++) delete m_RunWaveform.at(a);
    std::cout << setfill('-') << setw(89) << "-" << std::endl << setfill(' ');
    //std::cout << __FUNCTION__ << " INFO: Cleaning up objeects..." << std::endl;
    //for (unsigned int a = 0; a < m_nchan; a++) delete m_RunDUTCh.at(a);
}
// --------------------------------------------------------------------------------------------------------------
bool LGADRun::SetPlaneDQuants(std::vector<double> dQunat, std::vector<std::pair<std::pair<int, int>, double>> DQs)
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
                          if (((DQs.at(gh)).first).first == (int)(m_channels.at(k)) && ((DQs.at(gh)).first).second == (int)(m_channels.at(g)))
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
// --------------------------------------------------------------------------------------------------------------
void LGADRun::FirstInit()
{
  // Create new DUT channels
  m_RunDUTCh.resize(m_nchan, nullptr);
  std::vector<bool> foundPr[12];
  for (unsigned int h = 0; h < 12; h++) foundPr[h].resize(m_nchan, false);
  unsigned int sizes[12];
  sizes[0] = (m_RunBase->LGADBase::GetDUTNames()).size();
  sizes[1] = (m_RunBase->LGADBase::GetDUTBoards()).size();
  sizes[2] = (m_RunBase->LGADBase::GetDUTransImps()).size();
  sizes[3] = (m_RunBase->LGADBase::GetDUTSecStages()).size();
  sizes[4] = (m_RunBase->LGADBase::GetDUTAmpGains()).size();
  sizes[5] = (m_RunBase->LGADBase::GetDUTriggs()).size();
  sizes[6] = (m_RunBase->LGADBase::GetDUTFracts()).size();
  sizes[7] = (m_RunBase->LGADBase::GetDUTCaps()).size();
  sizes[8] = (m_RunBase->LGADBase::GetChVoltCuts()).size();
  sizes[9] = (m_RunBase->LGADBase::GetChNoiseCuts()).size();
  sizes[10] = (m_RunBase->LGADBase::GetChJitterCuts()).size();
  sizes[11] = (m_RunBase->LGADBase::GetChChargeCuts()).size();

  for (unsigned int a = 0; a < m_nchan; a++)
      {
       m_RunDUTCh.at(a) = new DUTChannel(m_channels.at(a), m_RunBase);
       for (unsigned int u = 0; u < 12; u++)
           {
            for (unsigned int w = 0; w < sizes[u]; w++) 
                {
                 if (u == 0)
                    {
                     if ((int)(m_channels.at(a)) == ((m_RunBase->LGADBase::GetDUTNames()).at(w)).first) 
                        {
                         m_RunDUTCh.at(a)->DUTChannel::SetChName(((m_RunBase->LGADBase::GetDUTNames()).at(w)).second);
                         (foundPr[u]).at(a) = true;
                        }
                    }
                 if (u == 1)
                    {
                     if ((int)(m_channels.at(a)) == ((m_RunBase->LGADBase::GetDUTBoards()).at(w)).first)
                        {
                         m_RunDUTCh.at(a)->DUTChannel::SetChBoard(((m_RunBase->LGADBase::GetDUTBoards()).at(w)).second);
                         (foundPr[u]).at(a) = true;
                        }
                    }
                 if (u == 2)
                    {
                     if ((int)(m_channels.at(a)) == ((m_RunBase->LGADBase::GetDUTransImps()).at(w)).first)
                        {
                         m_RunDUTCh.at(a)->DUTChannel::SetChTransImp(((m_RunBase->LGADBase::GetDUTransImps()).at(w)).second);
                         (foundPr[u]).at(a) = true;
                        }
                    }
                 if (u == 3)
                    {
                     if ((int)(m_channels.at(a)) == ((m_RunBase->LGADBase::GetDUTSecStages()).at(w)).first)
                        {
                         m_RunDUTCh.at(a)->DUTChannel::SetSecStage(((m_RunBase->LGADBase::GetDUTSecStages()).at(w)).second);
                         (foundPr[u]).at(a) = true;
                        }
                    }
                 if (u == 4)
                    {
                     if ((int)(m_channels.at(a)) == ((m_RunBase->LGADBase::GetDUTAmpGains()).at(w)).first)
                        {
                         m_RunDUTCh.at(a)->DUTChannel::SetChAmpGain(((m_RunBase->LGADBase::GetDUTAmpGains()).at(w)).second);
                         (foundPr[u]).at(a) = true;
                        }
                    }
                 if (u == 5)
                    {
                     if ((int)(m_channels.at(a)) == ((m_RunBase->LGADBase::GetDUTriggs()).at(w)).first)
                        {
                         m_RunDUTCh.at(a)->DUTChannel::SetChTrigg(((m_RunBase->LGADBase::GetDUTriggs()).at(w)).second);
                         (foundPr[u]).at(a) = true;
                        }
                    }
                 if (u == 6)
                    {
                     if ((int)(m_channels.at(a)) == ((m_RunBase->LGADBase::GetDUTFracts()).at(w)).first)
                        {
                         m_RunDUTCh.at(a)->DUTChannel::SetChFract(((m_RunBase->LGADBase::GetDUTFracts()).at(w)).second);
                         (foundPr[u]).at(a) = true;
                        }
                    }
                 if (u == 7)
                    {
                     if ((int)(m_channels.at(a)) == ((m_RunBase->LGADBase::GetDUTCaps()).at(w)).first)
                        {
                         m_RunDUTCh.at(a)->DUTChannel::SetChCap(((m_RunBase->LGADBase::GetDUTCaps()).at(w)).second);
                         (foundPr[u]).at(a) = true;
                        }
                    }
                 if (u == 8)
                    {
                     if ((int)(m_channels.at(a)) == ((m_RunBase->LGADBase::GetChVoltCuts()).at(w)).first)
                        {
                         m_RunDUTCh.at(a)->DUTChannel::SetVoltCut(((m_RunBase->LGADBase::GetChVoltCuts()).at(w)).second);
                         (foundPr[u]).at(a) = true;
                        }
                    }
                 if (u == 9)
                    {
                     if ((int)(m_channels.at(a)) == ((m_RunBase->LGADBase::GetChNoiseCuts()).at(w)).first)
                        {
                         m_RunDUTCh.at(a)->DUTChannel::SetNoiseCut(((m_RunBase->LGADBase::GetChNoiseCuts()).at(w)).second);
                         (foundPr[u]).at(a) = true;
                        }
                    }
                 if (u == 10)
                    {
                     if ((int)(m_channels.at(a)) == ((m_RunBase->LGADBase::GetChJitterCuts()).at(w)).first)
                        {
                         m_RunDUTCh.at(a)->DUTChannel::SetJitterCut(((m_RunBase->LGADBase::GetChJitterCuts()).at(w)).second);
                         (foundPr[u]).at(a) = true;
                        }
                    }
                 if (u == 11)
                    {
                     if ((int)(m_channels.at(a)) == ((m_RunBase->LGADBase::GetChChargeCuts()).at(w)).first)
                        {
                         m_RunDUTCh.at(a)->DUTChannel::SetChargeCut(((m_RunBase->LGADBase::GetChChargeCuts()).at(w)).second);
                         (foundPr[u]).at(a) = true;
                        }
                    }
                } 
            }
      }

  if (m_RunBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Intialized Channel ID properties for " << m_nchan << " channels" << std::endl;

  std::vector<unsigned int> skip[12];
  // Fill the properties of the channel without index
  for (unsigned int a = 0; a < m_nchan; a++)
      {
       for (unsigned int u = 0; u < 12; u++)
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
                     if (u == 8 && ((m_RunBase->LGADBase::GetChVoltCuts()).at(w)).first == -99)
                        {
                         m_RunDUTCh.at(a)->DUTChannel::SetVoltCut(((m_RunBase->LGADBase::GetChVoltCuts()).at(w)).second);
                         (skip[u]).push_back(w);
                        }
                     if (u == 9 && ((m_RunBase->LGADBase::GetChNoiseCuts()).at(w)).first == -99)
                        {
                         m_RunDUTCh.at(a)->DUTChannel::SetNoiseCut(((m_RunBase->LGADBase::GetChNoiseCuts()).at(w)).second);
                         (skip[u]).push_back(w);
                        }
                     if (u == 10 && ((m_RunBase->LGADBase::GetChJitterCuts()).at(w)).first == -99)
                        {
                         m_RunDUTCh.at(a)->DUTChannel::SetJitterCut(((m_RunBase->LGADBase::GetChJitterCuts()).at(w)).second);
                         (skip[u]).push_back(w);
                        }
                     if (u == 11 && ((m_RunBase->LGADBase::GetChChargeCuts()).at(w)).first == -99)
                        {
                         m_RunDUTCh.at(a)->DUTChannel::SetChargeCut(((m_RunBase->LGADBase::GetChChargeCuts()).at(w)).second);
                         (skip[u]).push_back(w);
                        }
                    }
               }
           }
      }

  // If we want to open several fiels we need to scan ahead and identify the files. Then reserve the apropriate vectors and branches and then access then one by one 
  // in  a non-selector loop, combining each one with the apropriate tracking file. A selector loop is not convinient in this case unless we open all of them at the
  // initialise function, declare and create branches there and then close the files and allow the notify fucntion to oproceed. The problem is what

  VarInit(m_base, 0);

  if (m_RunBase->LGADBase::GetInRootFile() != m_RunBase->LGADBase::GetOutRootFile()) 
     {
      if (!(m_RunBase->LGADBase::GetOutRootFile()->IsOpen())) m_RunBase->LGADBase::GetOutRootFile()->ReOpen("UPDATE");
      m_RunBase->LGADBase::GetOutRootFile()->SetCompressionLevel(6);
      m_OutTree = new TTree("wfm", "Per-event quantities of charge collection measurements");
      std::cout << m_RunBase->LGADBase::GetOutRootFile()->IsOpen() << " " << m_RunBase->LGADBase::GetOutRootFile()->GetName() << " " << m_OutTree->GetName() << std::endl;
      InitCopyFile(m_base);
     }
  else {
        m_RunBase->LGADBase::GetInRootFile()->ReOpen("UPDATE");
        m_RunBase->LGADBase::GetInRootFile()->SetCompressionLevel(6);
        m_OutTree = m_RunBase->LGADBase::GetRootTree();
       }
  m_OutTree->SetAutoFlush(-30000000); // Set autoflush to 30 MB
  AddOutputBranches(m_base);
}
// --------------------------------------------------------------------------------------------------------------
void LGADRun::VarInit(LGADBase* lbase, int level)
{
  // Initialise event variables
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
}
// --------------------------------------------------------------------------------------------------------------
void LGADRun::SecondStageAnalysis(TTree *ChFitParam)
{
  // Tree for Channel properties form Fits
  TTree *ChShapeParam = nullptr;
  ChFitParam->SetAutoFlush(10000); // Set autoflush to 30 MB

  std::vector<TList* > ChCanv;
  std::vector<TList* > CFDTimeCanv;
  std::vector<TList* > CFDToTCanv;
  std::vector<TList* > CFDdVdTCanv;

  // Channel Shape properties
  if (m_RunBase->LGADBase::GetWaveShape())
     {
      ChShapeParam = new TTree("ChShapeParam", "Channel parameters form wwaveform shape");
      ChShapeParam->SetAutoFlush(10000); // Set autoflush to 30 MB
     }

  for (unsigned int ich = 0; ich < m_nchan; ich++)
      {
       std::cout << __FUNCTION__ << " INFO: Updating Channel " << ich + 1 << " of " << m_nchan << " designated: " << m_channels.at(ich) << std::endl;
       m_RunDUTCh.at(ich)->DUTChannel::updateChProperties(m_RunBase->LGADBase::GetWaveShape(), m_OutTree);
       if (m_RunBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: After channel update: " << ich + 1 << "/" << m_nchan << std::endl;
       std::cout << "---------------------------------------------------------------------------------------------------------------------------" << std::endl;
       ChCanv.push_back(new TList());
       CFDTimeCanv.push_back(new TList());
       CFDToTCanv.push_back(new TList());
       CFDdVdTCanv.push_back(new TList());

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
       if (m_RunBase->LGADBase::GetDoFFT())
          {
           ChFitSignalFFT = m_RunDUTCh.at(ich)->DUTChannel::GetChSignalFFT(2);
           ChFitNoiseFFT = m_RunDUTCh.at(ich)->DUTChannel::GetChNoiseFFT(2);
          }
       ChFitComplete = m_RunDUTCh.at(ich)->DUTChannel::GetChComplete();
       ChFitSoNR = m_RunDUTCh.at(ich)->DUTChannel::GetChSoNR(2);
       if (m_RunBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Recovered fit variables for channel: " << ich +1 << "/" << m_nchan << std::endl;

       // Get the fits and relevent chi2 to the Ntuple
       MaxIndxFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("MaxIndx");
       if (MaxIndxFtChi2 != -99 && MaxIndxFtChi2 != -1)
          {
           if (LGADBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("MaxIndx"), &h_MaxIndxFt.at(ich))) 
              {
               ChCanv.at(ich)->Add(h_MaxIndxFt.at(ich));
              }
          }
       MinIndxFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("MinIndx");
       if (MinIndxFtChi2 != -99 && MinIndxFtChi2 != -1)
          {
           if (LGADBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("MinIndx"), &h_MinIndxFt.at(ich)))
              {
               ChCanv.at(ich)->Add(h_MinIndxFt.at(ich));
              }
         }
       MaxVoltFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("MaxVolt");
       if (MaxVoltFtChi2 != -99 && MaxVoltFtChi2 != -1)
          {
           if (LGADBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("MaxVolt"), &h_MaxVoltFt.at(ich)))
              {
               ChCanv.at(ich)->Add(h_MaxVoltFt.at(ich));
              }
          }
       MinVoltFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("MinVotl");
       if (MinVoltFtChi2 != -99 && MinVoltFtChi2 != -1)
          {
           if (LGADBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("MinVotl"), &h_MinVoltFt.at(ich)))
              {
               ChCanv.at(ich)->Add(h_MinVoltFt.at(ich));
              }
          }
       StrIndxFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("Strindx");
       if (StrIndxFtChi2 != -99 && StrIndxFtChi2 != -1)
          {
           if (LGADBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("Strindx"), &h_StrIndxFt.at(ich)))
              {
               ChCanv.at(ich)->Add(h_StrIndxFt.at(ich));
              }
          }
       EndIndxFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("EndIndx");
       if (EndIndxFtChi2 != -99 && EndIndxFtChi2 != -1)
          {
           if (LGADBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("EndIndx"), &h_EndIndxFt.at(ich)))
              {
               ChCanv.at(ich)->Add(h_EndIndxFt.at(ich));
              }
          }
       NoiseFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("Noise");
       if (NoiseFtChi2 != -99 && NoiseFtChi2 != -1)
          {
           if (LGADBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("Noise"), &h_NoiseFt.at(ich)))
              {
               ChCanv.at(ich)->Add(h_NoiseFt.at(ich));
              }
          }
       NoiseErrFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("NoisErr");
       if (NoiseErrFtChi2 != -99 && NoiseErrFtChi2 != -1)
          {
           if (LGADBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("NoisErr"), &h_NoiseErrFt.at(ich)))
              {
               ChCanv.at(ich)->Add(h_NoiseErrFt.at(ich));
              }
          }
       PedestFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("Pedestl");
       if (PedestFtChi2 != -99 && PedestFtChi2 != -1)
          {
           if (LGADBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("Pedestl"), &h_PedestFt.at(ich)))
              {
               ChCanv.at(ich)->Add(h_PedestFt.at(ich));
              }
          }
       PedestErrFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("PdslErr");
       if (PedestErrFtChi2 != -99 && PedestErrFtChi2 != -1)
          {
           if (LGADBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("PdslErr"), &h_NPedestErrFt.at(ich)))
              {
               ChCanv.at(ich)->Add(h_NPedestErrFt.at(ich));
              }
          }
       MaxTimeFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("MaxTime");
       if (MaxTimeFtChi2 != -99 && MaxTimeFtChi2 != -1)
          {
           if (LGADBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("MaxTime"), &h_MaxTimeFt.at(ich)))
              {
               ChCanv.at(ich)->Add(h_MaxTimeFt.at(ich));
              }
          }
       MinTimeFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("MinTime");
       if (MinTimeFtChi2 != -99 && MinTimeFtChi2 != -1)
          {
           if (LGADBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("MinTime"), &h_MinTimeFt.at(ich)))
              {
               ChCanv.at(ich)->Add(h_MinTimeFt.at(ich));
              }
          }
       ChargeFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("Charge");
       if (ChargeFtChi2 != -99 && ChargeFtChi2 != -1)
          {
           if (LGADBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("Charge"), &h_ChargeFt.at(ich)))
              {
               ChCanv.at(ich)->Add(h_ChargeFt.at(ich));
              }
          }
       RiseTimeFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("RiseTime");
       if (RiseTimeFtChi2 != -99 && RiseTimeFtChi2 != -1)
          {
           if (LGADBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("RiseTime"), &h_RiseTimeFt.at(ich)))
              {
               ChCanv.at(ich)->Add(h_RiseTimeFt.at(ich));
              }
          }
       TriggTimeFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("TriggTime");
       if (TriggTimeFtChi2 != -99 && TriggTimeFtChi2 != -1)
          {
           if (LGADBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("TriggTime"), &h_TriggTimeFt.at(ich)))
              {
               ChCanv.at(ich)->Add(h_TriggTimeFt.at(ich));
              }
          }
       DVDTMaxFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("DVDTMax");
       if (DVDTMaxFtChi2 != -99 && DVDTMaxFtChi2 != -1)
          {
           if (LGADBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("DVDTMax"), &h_DVDTMaxFt.at(ich)))
              {
               ChCanv.at(ich)->Add(h_DVDTMaxFt.at(ich));
              }
          }
       TriggToTFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("TrigToT");
       if (TriggToTFtChi2 != -99 && TriggToTFtChi2 != -1)
          {
           if (LGADBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("TrigToT"), &h_TriggToTFt.at(ich)))
              {
               ChCanv.at(ich)->Add(h_TriggToTFt.at(ich));
              }
          }

       if (m_RunBase->LGADBase::GetDoFFT())
          {
           SignalFFTFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("SigFFT");
           if (SignalFFTFtChi2 != -99 && SignalFFTFtChi2 != -1)
              {
               if (LGADBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("SigFFT"), &h_SignalFFTFt.at(ich)))
                  {
                   ChCanv.at(ich)->Add(h_SignalFFTFt.at(ich));
                  }
              }
           NoiseFFTFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("NoiseFFT");
           if (NoiseFFTFtChi2 != -99 && NoiseFFTFtChi2 != -1)
              {
               if (LGADBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("NoiseFFT"), &h_NoiseFFTFt.at(ich)))
                  {
                   ChCanv.at(ich)->Add(h_NoiseFFTFt.at(ich));
                  }
              }
          }

       for (unsigned int h = 0; h < 19; h++)
           {
            char title[20];
            CFDTimeFtChi2.at(h) = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("CFDTime", h);
            if (CFDTimeFtChi2.at(h) != -99 && CFDTimeFtChi2.at(h) != -1)
               {
                if (LGADBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("CFDTime", h), &(h_CFDTimeFt.at(ich)).at(h)))
                   {
                    sprintf(title, "CFDTime%02u-%02u%%", m_channels.at(ich), (h*5+5));
                    (h_CFDTimeFt.at(ich)).at(h)->SetName((const char*)title);
                    CFDTimeCanv.at(ich)->Add((h_CFDTimeFt.at(ich)).at(h));
                   }
               }
            DVDTCFDFtChi2.at(h) = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("DVDTCFD", h);
            if (DVDTCFDFtChi2.at(h) != -99 && DVDTCFDFtChi2.at(h) != -1)
               {
                if (LGADBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("DVDTCFD", h), &(h_DVDTCFDFt.at(ich)).at(h)))
                   {
                    memset(title, '0', sizeof(title));
                    sprintf(title, "DVDTCFD%02u-%02u%%", m_channels.at(ich), (h * 5 + 5));
                    (h_DVDTCFDFt.at(ich)).at(h)->SetName((const char*)title);
                    CFDdVdTCanv.at(ich)->Add((h_DVDTCFDFt.at(ich)).at(h));
                   }
              }
            CFDToTFtChi2.at(h) = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("CFDToT", h);
            if (CFDToTFtChi2.at(h) != -99 && CFDToTFtChi2.at(h) != -1)
               {
                if (LGADBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("CFDToT", h), &(h_CFDToTFt.at(ich)).at(h)))
                   {
                    memset(title, '0', sizeof(title));
                    sprintf(title, "CFDToT%02u-%02u%%", m_channels.at(ich), (h * 5 + 5));
                    (h_CFDToTFt.at(ich)).at(h)->SetName((const char*)title);
                    CFDToTCanv.at(ich)->Add((h_CFDToTFt.at(ich)).at(h));
                    if (h == 18) (h_CFDToTFt.at(ich)).at(h)->Draw();
                   }
               }
           }
       if (m_RunBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Recovered fits for channel: " << ich + 1 << "/" << m_nchan << std::endl;
         
       ChFitParam->Fill();

       // Mean pulseShape properties if shape is calculated
       if (m_RunBase->LGADBase::GetWaveShape())
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
           if (m_RunBase->LGADBase::GetDoFFT())
              {
               ChShapeSignalFFT = (m_RunDUTCh.at(ich)->DUTChannel::GetChSignalFFT(1)).first;
               ChShapeNoiseFFT = (m_RunDUTCh.at(ich)->DUTChannel::GetChNoiseFFT(1)).first;
              }
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
       m_RunBase->LGADBase::GetOutRootFile()->mkdir(Form("ChHistos%02u", m_channels.at(ich)));
       m_RunBase->LGADBase::GetOutRootFile()->cd(Form("ChHistos%02u", m_channels.at(ich)));
       ChCanv.at(ich)->Write();

       m_RunBase->LGADBase::GetOutRootFile()->mkdir(Form("ChHistos%02u/CFDTimeHists%02u", m_channels.at(ich), m_channels.at(ich)));
       m_RunBase->LGADBase::GetOutRootFile()->cd(Form("ChHistos%02u/CFDTimeHists%02u", m_channels.at(ich), m_channels.at(ich)));
       CFDTimeCanv.at(ich)->Write();
       m_RunBase->LGADBase::GetOutRootFile()->mkdir(Form("ChHistos%02u/CFDdVdTHists%02u", m_channels.at(ich), m_channels.at(ich)));
       m_RunBase->LGADBase::GetOutRootFile()->cd(Form("ChHistos%02u/CFDdVdTHists%02u", m_channels.at(ich), m_channels.at(ich)));
       CFDdVdTCanv.at(ich)->Write();
       m_RunBase->LGADBase::GetOutRootFile()->mkdir(Form("ChHistos%02u/CFDToTHists%02u", m_channels.at(ich), m_channels.at(ich)));
       m_RunBase->LGADBase::GetOutRootFile()->cd(Form("ChHistos%02u/CFDToTHists%02u", m_channels.at(ich), m_channels.at(ich)));
       CFDToTCanv.at(ich)->Write();
       m_RunBase->LGADBase::GetOutRootFile()->cd();
      }

  ChFitParam->Write();
  // deleting the histogram lists, deletion happens from end to start
  for (unsigned int a = m_nchan-1; a == 0; a--)
      {
       delete CFDdVdTCanv.at(a);
       delete CFDToTCanv.at(a);
       delete CFDTimeCanv.at(a);
       delete ChCanv.at(a);
      }

  if (m_RunBase->LGADBase::GetWaveShape())
     {
      ChShapeParam->Write();
      delete ChShapeParam;
     }
}
// --------------------------------------------------------------------------------------------------------------
bool LGADRun::ThirdStageProcess()
{
   // Carefull this is third level analyis that is here and obviously won't work as it referes to channel combinaitons
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

   return true;
}
// --------------------------------------------------------------------------------------------------------------
void LGADRun::ThirdStageAnalysis(TTree* tree)
{
    std::vector<TList* > TimeDiffCanv;
    TList* histoslist = new TList();
    std::vector<std::vector<double>> TmDiffChi2;
    std::vector<std::vector<std::pair<double, double>>> TmDiffFitQts;

    // Estimate the number of unique combinations between different channels
    TmDiffChi2.clear();
    TmDiffFitQts.clear();
    for (unsigned int ich = 0; ich < m_nchan; ich++)
        {
         for (unsigned int w = ich + 1; w < m_nchan; w++) 
             {
              if (m_RunBase->LGADBase::GetInstrument() == TestBeamBin1 || m_RunBase->LGADBase::GetInstrument() == TestBeamBin2 || m_RunBase->LGADBase::GetInstrument() == LeCroyWRBin)
                 { 
                  if (m_scope(ich) != m_scope.at(w)) continue;
                 }
              TmDiffChi2.push_back(std::vector<double>());
              TmDiffFitQts.push_back(std::vector<std::pair<double, double>>());
              tree->Branch(Form("TmDiffChi2-%02u/%02u", m_channels.at(ich), m_channels.at(w)), &TmDiffChi2.back());
              tree->Branch(Form("TmDiffFitQts-%02u/%02u", m_channels.at(ich), m_channels.at(w)), &TmDiffFitQts.back());
              (TmDiffChi2.back()).resize(19*19, -99.0);
              (TmDiffFitQts.back()).resize(19*19, std::make_pair(-99, -99));
             }
        }

    for (unsigned int ich = 0; ich < m_nchan; ich++)
        {
         unsigned int a = 0;
         unsigned int b = 0;
         // Lets create the Time difference histos and polupate the time difference sigmas histo
         std::cout << "Probing channel: " << ich+1 << " from " << ich+1 << " to " << m_nchan << std::endl;
         for (unsigned int ich2 = ich + 1; ich2 < m_nchan; ich2++)
             {
              TimeDiffCanv.push_back(new TList());
              a = ((float)(ich*(ich2-1)))/2;
              std::cout << ich2 << " " << m_nchan  << " alpha: " << a << std::endl;
              std::pair <double, double> mean;
              std::pair <double, double> sigma;
              int qual;
              char title[24];
              h_TmSigmas.at(a) = new TH2D(Form("h_TimeMap%02u-%02u", m_channels.at(ich), m_channels.at(ich2)), 
                                          Form("CFD Map, Channels %02u-%02u;CFD %% Channel %02u;CFD %% Channel %02u;#sigma_{t(%02u)-t(%02u)}", 
                                               m_channels.at(ich), m_channels.at(ich2), m_channels.at(ich), m_channels.at(ich2), m_channels.at(ich), m_channels.at(ich2)),
                                               19, 0.025, 0.975, 19, 0.025, 0.975); // initialize the 2D sigmas histos
              h_TmSigmas.at(a)->GetXaxis()->SetLabelSize(13);
              h_TmSigmas.at(a)->GetYaxis()->SetLabelSize(13);
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
                                 std::cout << a << " " << b << " " << m_EvTmDiff.size() << " " << (m_EvTmDiff.at(a)).size() << " " << h_TmDiffCFD.size() << " " << (h_TmDiffCFD.at(a)).size() << " " << TmDiffChi2.size() << " " << (TmDiffChi2.at(a)).size() << std::endl;
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
                            snprintf(title, sizeof(title), "TmDiffCFD%02u%%-%02u%%CH%02u-%02u", bin1*5+5, bin2*5+5, m_channels.at(ich), m_channels.at(ich2));
                            (h_TmDiffCFD.at(a)).at(b)->SetName((const char*)title);
                            (h_TmDiffCFD.at(a)).at(b)->SetTitle(Form("#Deltat Channel %02u (%02u%% CFD) - Channel %02u (%02u%% CFD);#DeltaT_{%02u-%02u} [sec];Entries", 
                                                                     m_channels.at(ich), bin1*5+5, m_channels.at(ich2), bin2*5+5, m_channels.at(ich), m_channels.at(ich2)));
                            if (LGADBase::PrintFitInfo((h_TmDiffCFD.at(a)).at(b), &(h_TmDiffCFDCanv.at(a)).at(b)))
                               {
                                TimeDiffCanv.at(ich)->Add((h_TmDiffCFDCanv.at(a)).at(b));
                               }
                            h_TmSigmas.at(a)->SetBinContent(bin1 + 1, bin2 + 1, ((TmDiffFitQts.at(a)).at(b)).second);
                           }
                        else h_TmSigmas.at(a)->SetBinContent(bin1 + 1, bin2 + 1, 0.0);
                       }
                  }
              histoslist->Add(h_TmSigmas.at(a));
             }
         tree->Fill();
         for (unsigned int gp = ich + 1; gp < m_nchan; gp++) 
             {
              m_RunBase->LGADBase::GetOutRootFile()->mkdir(Form("ChHistos%02u/TimeDiff%02u-%02u", m_channels.at(ich), m_channels.at(ich), m_channels.at(gp)));
              m_RunBase->LGADBase::GetOutRootFile()->cd(Form("ChHistos%02u/TimeDiff%02u-%02u", m_channels.at(ich), m_channels.at(ich), m_channels.at(gp)));
              TimeDiffCanv.at(((float)(ich*(gp-1)))/2)->Write();
             }
        }

    // Chi2 system level minimization algorithm

    for (unsigned int a = m_nchan - 1; a == 0; a--)
        {
         for (unsigned int gp = m_nchan - 1; gp == a + 1; gp--) delete TimeDiffCanv.at(((float)(a * (gp - 1))) / 2);
        }
    histoslist->Write();
    delete histoslist;
}
// --------------------------------------------------------------------------------------------------------------
void LGADRun::TrckEffAnalaysis()
{

}
// --------------------------------------------------------------------------------------------------------------