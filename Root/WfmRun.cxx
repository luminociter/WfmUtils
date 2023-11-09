/*
*                 WFMRun.cxx
*
*
*      Author: Gkougkousis Evangelos - Leonidas
*             egkougko@cern.ch - 2023
*                  CERN-GENEVA
*/

#include "../WfmUtils/WfmRun.h"

//#if !defined(__CINT__)
//ClassImp(WfmRun);
//#endif
//
//#if !defined(__CLING__)
//ClassImp(WfmRun);
//#endif
//
//#ifdef __CINT__
//#pragma link C++ class WfmRun;
//#endif
//
//#ifdef __ROOTCLING__
//#pragma link C++ class WfmRun;
//#endif

WfmRun::WfmRun()
{
    m_RunBase = new WfmBase();
    m_dTMax.clear();
    m_dCMax.clear();
    std::cout << setfill('-') << setw(89) << "-" << setfill(' ') << std::endl;
}
// --------------------------------------------------------------------------------------------------------------
WfmRun::WfmRun(WfmBase* rBase)
{
    WfmBase::InportProperties(*rBase); // setting the basic properties WfmRun object
    m_RunBase = rBase; // createing a class-wide pointer of WfmBase connected to the prototype
    m_dTMax.clear();
    m_dCMax.clear();
    if (m_RunBase->WfmBase::GetVerbosity() >= 2) std::cout << "Wfm Run constructor called and completed!" << std::endl;
    std::cout << setfill('-') << setw(89) << "-" << setfill(' ') << std::endl;
}
// --------------------------------------------------------------------------------------------------------------
WfmRun::~WfmRun()
{
    delete m_RunBase;
    for (unsigned int gh = 0; gh < m_RunDUTCh.size(); gh++) delete m_RunDUTCh.at(gh);
}
// --------------------------------------------------------------------------------------------------------------
// A vector with the index of xorretly opened files is provided for matching with tracking
bool WfmRun::InitAnalysis(std::vector<unsigned int>& fileindx)
{
    // Iniialise all the pointers and estalish analysis stage and options of input files
    if (!WfmRun::InitNtuple(m_RunBase->WfmBase::GetRootTrees(), m_RunBase->WfmBase::GetInRootFiles()))
       {
        std::cout << __FUNCTION__ << " ERROR: Failed at ntuple list processing, aborting..." << std::endl;
        return false;
       }
    m_event = 0;

    std::vector<unsigned int> l_fileindx_trck;
    if (m_RunBase->WfmBase::GetTrackComb())
       {
        // Exclude files set to not be combined with tracking ntpules
        l_fileindx_trck = fileindx;
        if ((m_RunBase->WfmBase::GetExcludeTrackFiles()).size() > 0)
           {
            for (unsigned int q = 0; q < l_fileindx_trck.size(); q++)
                {
                 for (unsigned int p = 0; p < (m_RunBase->WfmBase::GetExcludeTrackFiles()).size(); p++)
                     {
                      if (l_fileindx_trck.at(q) == m_RunBase->WfmBase::GetExcludeTrackFiles().at(p))
                         {
                          if (m_filetrack.at(l_fileindx_trck.at(q)) != 4)
                             {
                              l_fileindx_trck.erase(l_fileindx_trck.begin() + q);
                              q--;
                             }
                          else {
                                std::cout << __FUNCTION__ << " WARNING: File " << (m_RunBase->WfmBase::GetInRootFile(l_fileindx_trck.at(q)))->GetName()
                                          << " has not tracking branches but has been specified to not be included in tacking combination. Ignoring..." << std::endl;
                               }
                          break;
                         }
                     }
                }
           }
        std::vector<unsigned int> l_fileindx_trck_old = l_fileindx_trck;
        // Loading tracking functions, opening relevent tracking files, initializing pointers and vectors for reading track data
        if (!WfmBase::OpenTrackFile(m_RunBase->WfmBase::GetTrackInDir(), m_RunBase->WfmBase::GetTrackDataName(), l_fileindx_trck))
           {
            std::cout << __FUNCTION__ << " ERROR: Failed oppening tracking files in directory: " 
                      << m_RunBase->WfmBase::GetTrackInDir() << std::endl;
            return false;
           }
        // Handle files that were not oppened by the tracking function and exclude them from loop
        if (l_fileindx_trck_old.size() != l_fileindx_trck.size())
           {
            std::vector<unsigned int> missing_indxs;
            bool fnd = false;
            for (unsigned int we = 0; we < l_fileindx_trck_old.size(); we++)
                {
                 fnd = false;
                 for (unsigned int yu = 0; yu < l_fileindx_trck.size(); yu++) if (l_fileindx_trck.at(yu) == l_fileindx_trck_old.at(we)) { fnd = true; break; }
                 if (!fnd) missing_indxs.push_back(we);
                }
            // Need to remove the relevent root files from the input file vector and trees and then call the InitNtuple() function again
            for (unsigned int rt = 0; rt < missing_indxs.size(); rt++)
                {
                 (m_RunBase->WfmBase::GetInRootFiles()).erase((m_RunBase->WfmBase::GetInRootFiles()).begin() + missing_indxs.at(rt));
                 (m_RunBase->WfmBase::GetRootTrees()).erase((m_RunBase->WfmBase::GetRootTrees()).begin() + missing_indxs.at(rt));
                 fileindx.erase(fileindx.begin() + missing_indxs.at(rt));
                }
            if (!WfmRun::InitNtuple(m_RunBase->WfmBase::GetRootTrees(), m_RunBase->WfmBase::GetInRootFiles()))
               {
                std::cout << __FUNCTION__ << " ERROR: Failed at ntuple list processing, aborting..." << std::endl;
                return false;
               }
           }
        // Checking for tracking packages consistency and asigning the new tracking package tot he common base
        if (m_RunBase->WfmBase::GetHasTrck() && (WfmBase::GetTrackPackage() != m_RunBase->WfmBase::GetTrackPackage()))
           {
            std::cout << __FUNCTION__ << " ERROR: Trying to mix files from different tracking algorithms: " << m_RunBase->WfmBase::GetTrackPackage()
                << " and " << m_RunBase->WfmBase::GetTrackPackage() << "! Abording..." << std::endl;
            return false;
           }
        else m_RunBase->WfmBase::SetTrackPackage(WfmBase::GetTrackPackage());
        if (m_RunBase->WfmBase::GetTrackPackage() == "Cory")
           { 
            bool fd1 = false;
            for (unsigned int dp = 0; dp < m_planes.size(); dp++)
                {
                 fd1 = false;
                 for (unsigned int dp1 = 0; dp1 < m_ntp_planes.size(); dp1++) { if (m_planes.at(dp) == m_ntp_planes.at(dp1)) { fd1 = true; break; } }
                 if (!fd1) m_ntp_planes.push_back(m_planes.at(dp));
                }
            for (unsigned int dp = 0; dp < m_Refs.size(); dp++)
                {
                 fd1 = false;
                 for (unsigned int dp1 = 0; dp1 < m_ntp_Ref.size(); dp1++) { if (m_Refs.at(dp) == m_ntp_Ref.at(dp1)) { fd1 = true; break; } }
                 if (!fd1) m_ntp_Ref.push_back(m_Refs.at(dp));
                }
            for (unsigned int dp = 0; dp < m_DUTs.size(); dp++)
                {
                 fd1 = false;
                 for (unsigned int dp1 = 0; dp1 < m_ntp_DUT.size(); dp1++) { if (m_DUTs.at(dp) == m_ntp_DUT.at(dp1)) { fd1 = true; break; } }
                 if (!fd1) m_ntp_DUT.push_back(m_DUTs.at(dp));
                }
            WfmBase::InitCoryPointers(m_ntp_Ref.size(), m_ntp_DUT.size(), m_ntp_planes.size());
           }
       }

    // Create waveforms or DUTChannels for the analysis, setup the output tree and propagete any cut values 
    WfmRun::FirstInit(m_RunBase->WfmBase::GetFileLVL(), m_RunBase->WfmBase::GetAnaStage());

    if ((m_RunBase->WfmBase::GetPlaneDTs()).size() > 0) SetPlaneDQuants(m_dTMax, m_RunBase->GetPlaneDTs());
    if ((m_RunBase->WfmBase::GetPlaneDCs()).size() > 0) SetPlaneDQuants(m_dCMax, m_RunBase->GetPlaneDCs());

    // Create Sub-Channel vectors and expand properties
    if (((m_RunBase->WfmBase::GetFileLVL().first < 3  || m_RunBase->WfmBase::GetFileLVL().first == 4) && (m_RunBase->WfmBase::GetAnaStage() == 3 || m_RunBase->WfmBase::GetAnaStage() == 5)) ||
        (m_RunBase->WfmBase::GetFileLVL().first == 3 && m_RunBase->WfmBase::GetAnaStage() == 3) || ((m_RunBase->WfmBase::GetFileLVL().first == 5 && m_RunBase->WfmBase::GetAnaStage() == 5)))
       {
        m_RunDUTSubCh.clear();
        // if there is no tracking and there will be no combination, create the appropriate number of channels
        if (m_RunBase->WfmBase::GetTrackComb() || m_RunBase->WfmBase::GetHasTrck())
           {
            for (unsigned int q = 0; q < (m_RunBase->WfmBase::GetChannels()).size(); q++)
                {
                 // Create Square Channels
                 if (m_RunDUTCh.at(q)->DUTChannel::GetChBinSz() != -99 &&
                     (m_RunDUTCh.at(q)->DUTChannel::GetChLowLim()).first != -99 && (m_RunDUTCh.at(q)->DUTChannel::GetChLowLim()).second != -99 &&
                     (m_RunDUTCh.at(q)->DUTChannel::GetChUpLim()).first != -99 && (m_RunDUTCh.at(q)->DUTChannel::GetChUpLim()).second != -99)
                    {
                     m_RunDUTSubCh.push_back(std::make_pair((m_RunBase->WfmBase::GetChannels()).at(q), std::vector<DUTChannel*>()));
                     if (m_RunDUTCh.at(q)->DUTChannel::GetChBinShape() == "square")
                        { 
                         double area = fabs((m_RunDUTCh.at(q)->DUTChannel::GetChUpLim().first-(m_RunDUTCh.at(q)->DUTChannel::GetChLowLim()).first))*
                                       fabs((m_RunDUTCh.at(q)->DUTChannel::GetChUpLim().second-(m_RunDUTCh.at(q)->DUTChannel::GetChLowLim()).second));
                         unsigned int nbins = ceil((area*1e6)/m_RunDUTCh.at(q)->DUTChannel::GetChBinSz());
                         double dx = fabs((m_RunDUTCh.at(q)->DUTChannel::GetChUpLim().first - (m_RunDUTCh.at(q)->DUTChannel::GetChLowLim()).first)) / (double)nbins;
                         double dy = fabs((m_RunDUTCh.at(q)->DUTChannel::GetChUpLim().second - (m_RunDUTCh.at(q)->DUTChannel::GetChLowLim()).second)) / (double)nbins;
                         ((m_RunDUTSubCh.at(q)).second).resize(nbins, nullptr);
                         for (unsigned int k = 0; k < nbins; k++) 
                             {
                              ((m_RunDUTSubCh.at(q)).second).at(k) = new DUTChannel(m_RunDUTCh.at(q));
                              ((m_RunDUTSubCh.at(q)).second).at(k)->DUTChannel::SetChLowLim(std::make_pair((m_RunDUTCh.at(q)->DUTChannel::GetChLowLim()).first+k*dx,
                                                                                                           (m_RunDUTCh.at(q)->DUTChannel::GetChLowLim()).first+k*dy));
                              ((m_RunDUTSubCh.at(q)).second).at(k)->DUTChannel::SetChUpLim(std::make_pair((m_RunDUTCh.at(q)->DUTChannel::GetChLowLim()).first+(k+1)*dx, 
                                                                                                          (m_RunDUTCh.at(q)->DUTChannel::GetChLowLim()).second+(k+1)*dy));
                             }
                        }
                     else if (m_RunDUTCh.at(q)->DUTChannel::GetChBinShape() == "circle")
                             {
                              std::pair<double, double> center  = std::make_pair(((m_RunDUTCh.at(q)->DUTChannel::GetChLowLim()).first + (m_RunDUTCh.at(q)->DUTChannel::GetChUpLim()).first)/2,
                                                                                 ((m_RunDUTCh.at(q)->DUTChannel::GetChLowLim()).second+(m_RunDUTCh.at(q)->DUTChannel::GetChUpLim()).second)/2);
                              double distance = sqrt(pow(((m_RunDUTCh.at(q)->DUTChannel::GetChUpLim()).first-(m_RunDUTCh.at(q)->DUTChannel::GetChLowLim()).first), 2) + 
                                                     pow(((m_RunDUTCh.at(q)->DUTChannel::GetChUpLim()).second-(m_RunDUTCh.at(q)->DUTChannel::GetChLowLim()).second), 2));
                              unsigned int nbins = ceil((distance*1e3)/m_RunDUTCh.at(q)->DUTChannel::GetChBinSz());
                              ((m_RunDUTSubCh.at(q)).second).resize(nbins, nullptr);
                              for (unsigned int k = 0; k < nbins; k++)
                                  {
                                   ((m_RunDUTSubCh.at(q)).second).at(k) = new DUTChannel(m_RunDUTCh.at(q));
                                   ((m_RunDUTSubCh.at(q)).second).at(k)->DUTChannel::SetChLowLim(std::make_pair(center.first, center.second));
                                   ((m_RunDUTSubCh.at(q)).second).at(k)->DUTChannel::SetChUpLim(std::make_pair(nbins*m_RunDUTCh.at(q)->DUTChannel::GetChBinSz(),
                                                                                                               (nbins+1)* m_RunDUTCh.at(q)->DUTChannel::GetChBinSz()));
                                  }
                             }
                    }
                }
           }
       }

    // Initialise output branch pointers, output root file is already available and created within the Analyse function of general wrapper WfmUtils
    // Accessing output file via: m_RunBase->WfmBase::GetOutRootFile() or using m_ofile as the WfmRun base object is a copy of the m_base.
    if (!WfmRun::InitOutputBranches())
       {
        std::cout << __FUNCTION__ << " ERROR: Failed initializing output branch pointers, quitting.... " << std::endl;
        return false;
       }
    else if (!WfmRun::SetOutputBranches()) // Create the output Ntuple
            {
             std::cout << __FUNCTION__ << " ERROR: Failed creating the output Ntuple " << m_RunBase->WfmBase::GetOutRootFile()->GetName()
                       << "! Aborting ..." << std::endl;
             return false;
            }
    if (m_RunBase->WfmBase::GetTrackComb() && m_RunBase->WfmBase::GetTrackPackage() != "")
       {
        if (!WfmBase::AddTrckBranches(m_OutTree,m_RunBase->WfmBase::GetTrackPackage()))
           {
            std::cout << __FUNCTION__ << " Warning: Failed adding the tracking part on the output file: " << m_RunBase->WfmBase::GetOutRootFile()->GetName() << std::endl;
            std::cout << "                                    Switching to non-trackimg mode..." << std::endl;
            m_RunBase->WfmBase::SetTrackComb(false);
           }
       }

    std::cout << setfill('-') << setw(89) << "-" << std::endl << setfill(' ');

    // Correct event estimation thoughout all input files
    std::vector<int> fdTrackFile;
    fdTrackFile.resize(fileindx.size(), -99);
    std::vector<std::pair<unsigned int, unsigned int>> limEvnt;
    limEvnt.resize(fileindx.size(), std::make_pair(99, 99));
    m_totalevnts = 0;
    // Start - Stop event correction, loop over all files and all trees and get first and last event berfore closing each file
    for (unsigned int sp = 0; sp < fileindx.size(); sp++)
        {
         (m_RunBase->WfmBase::GetRootTree(sp))->SetBranchAddress("EvnNo", &EvnNo, &b_EvnNo);
         // Find if there is a tracking file to combine with the oscilloscope file we are processing
         if (m_RunBase->WfmBase::GetTrackComb())
            {
             for (unsigned int w = 0; w < l_fileindx_trck.size(); w++) 
                 { 
                  if (l_fileindx_trck.at(w) == fileindx.at(sp)) 
                     { 
                      fdTrackFile.at(sp) = w;
                      break; 
                     } 
                 }
            }
         unsigned int min_evnt, max_evnt;
         b_EvnNo->GetEntry(m_RunBase->WfmBase::GetRootTree(sp)->GetEntriesFast()-1);
         max_evnt = EvnNo;
         if (m_RunBase->WfmBase::GetTrackComb() && fdTrackFile.at(sp) != -99 && m_filetrack.at(sp) == 4) 
            {
             b_EvnNo->GetEntry(1);
             min_evnt = EvnNo - 1;
             max_evnt = max_evnt - 1;
            }
         else {
               b_EvnNo->GetEntry(0);
               min_evnt = EvnNo;
              }
         limEvnt.at(sp).first = m_RunBase->WfmBase::GetStartStopEvnt().first;
         limEvnt.at(sp).second = m_RunBase->WfmBase::GetStartStopEvnt().second;
         WfmBase::CheckStartStopEvnt((limEvnt.at(sp)).first, (limEvnt.at(sp)).second, min_evnt, max_evnt);
         m_totalevnts += (limEvnt.at(sp)).second - (limEvnt.at(sp)).first;
        }

    m_cp != -99
    // So these are relatedto min and max. If max has them, then recover them, if min does not though, then create them
    // We are doing first level in at least one of the input fiels (so check minimum with respect to analysis level) or we are copying to a new file
          
    if ((m_cp == -99 && (m_RunBase->WfmBase::GetAnaStage() > 0 || m_RunBase->WfmBase::GetFileLVL().second > 0)) ||
        (m_cp != -99 && m_RunBase->WfmBase::GetAnaStage() > 0 && (m_RunBase->WfmBase::GetFileLVL().first == 0 ||
            (m_RunBase->WfmBase::GetDoFFT() && (m_RunBase->WfmBase::GetHasFFT().at(0)).first == 0))))
       {
        if (m_cp == -99)
           {

           }
        else {
              
             }
        m_NormCount[5] = { 0, 0, 0, 0, 0 };
       }

    // Relevent vector necessery for track combination
    std::vector<unsigned int> trckevntlast; 
    bool frqhist = (m_RunBase->WfmBase::GetInstrument() != TestBeamBin1 && m_RunBase->WfmBase::GetInstrument() != LeCroyWRBin);
    std::vector<double> l_triggtimeVec;
    long long int l_evnt = 0;
    long long unsigned int l_prosEvnt = 0;
    m_event = 0;
    // Start looping on all input files. The vetor of input files has the ssame size as the vector of fileindx. That does ont necesseraly mean that 
    // all root files in the directory were succesfully opeened!
    for (unsigned int a = 0; a < fileindx.size(); a++)
        {
         if (a > 0) l_prosEvnt += (limEvnt.at(a)).second - (limEvnt.at(a)).first;
         std::cout << "--> File Accessed: " << m_RunBase->WfmBase::GetInRootFile(a)->GetName() << std::endl;
         l_evnt = m_RunBase->WfmBase::GetRootTree(a)->GetEntriesFast();
         std::cout << "--> Total number of entries in file: " << l_evnt << std::endl;

         // Set branch adresses for oppened Ntuple
         if (!WfmRun::SetBranchAccess(a))
            {
             std::cout << __FUNCTION__ << " ERROR: Failed while setting branch adresses for Ntuple " << m_RunBase->WfmBase::GetOutRootFile()->GetName()
                       << "! Skipping ..." << std::endl;
             continue;
            }

         // Fix vectors and variables if tracking file to combine with the oscilloscope file we are processing
         if (m_RunBase->WfmBase::GetTrackComb() && fdTrackFile.at(a) != -99)
            {
             trckevntlast.clear();
             if (m_RunBase->WfmBase::GetTrackPackage() == "Cory") trckevntlast.resize(1 + m_nplanes.at(a) + m_nDUT.at(a) + m_nRef.at(a), 0);
             else trckevntlast.push_back(0);
            }

         // If combining complete files or writing at an existing output file, no need to recalculate trigger rates, clone histograms
         if (a > 0 && frqhist && (((m_RunBase->WfmBase::GetStartStopEvnt().first == 0 && m_RunBase->WfmBase::GetStartStopEvnt().second == 0) || m_cp != -99) || m_RunBase->WfmBase::GetInstrument() == Sampic))
            {
             m_trigDt->Add(m_trigDt_In.at(a));
             m_trigFr->Add(m_trigFr_In.at(a));
            }
         else if (frqhist && (((m_RunBase->WfmBase::GetStartStopEvnt().first != 0 || m_RunBase->WfmBase::GetStartStopEvnt().second != 0) && m_cp == -99) || m_RunBase->WfmBase::GetInstrument() != Sampic))
                 {
                  l_triggtimeVec.clear();
                  l_triggtimeVec.resize(((limEvnt.at(a)).second-(limEvnt.at(a)).first), 0);
                 }

         unsigned int entry = 0;
         if (m_RunBase->WfmBase::GetTrackComb() && fdTrackFile.at(a) != -99 && m_filetrack.at(a) == 4) entry = 1;

         // Start the event loop for openned Ntuple
         for (; entry < l_evnt; entry++)
             {
              b_EvnNo->GetEntry(entry);
              // Check if event within readable scope
              if ((m_RunBase->WfmBase::GetTrackComb() && fdTrackFile.at(a) != -99 && m_filetrack.at(a) == 4 && (long int)(EvnNo-1) < (limEvnt.at(a)).first) ||
                  ((!m_RunBase->WfmBase::GetTrackComb() || (m_RunBase->WfmBase::GetTrackComb() && (fdTrackFile.at(a) == -99 || m_filetrack.at(a) != 4))) 
                      && (long int)EvnNo < (limEvnt.at(a)).first)) continue;
              if ((m_RunBase->WfmBase::GetTrackComb() && fdTrackFile.at(a) != -99 && m_filetrack.at(a) == 4 && (long int)(EvnNo-1) >= (limEvnt.at(a)).second) ||
                  ((!m_RunBase->WfmBase::GetTrackComb() || (m_RunBase->WfmBase::GetTrackComb() && (fdTrackFile.at(a) == -99 || m_filetrack.at(a) != 4))) 
                      && (long int)EvnNo >= (limEvnt.at(a)).second)) break;

              m_event++;
              WfmBase::ProgressBar(m_event, m_totalevnts);

              // Do Event tests for debugging
              if (m_RunBase->WfmBase::GetTestEvn(std::make_pair(EvnNo, a))) continue;

              if (m_RunBase->WfmBase::GetVerbosity() >= 2) 
                 {
                  std::cout << "-------------------------------------------------------------------" << std::endl;
                  std::cout << __FUNCTION__ << " INFO: Number of total channels: " << (m_RunBase->WfmBase::GetChannels()).size() << ", instrument: " << m_RunBase->WfmBase::GetInstrument() << std::endl;
                  std::cout << __FUNCTION__ << " INFO: Active channels in current file: " << (m_filechannel.at(a)).size() << ", event: " << entry << std::endl;
                  std::cout << __FUNCTION__ << " INFO: Channels: ";
                  for (unsigned int k = 0; k < (m_filechannel.at(a)).size(); k++) std::cout << (m_filechannel.at(a)).at(k) << " ";
                  std::cout << std::endl;
                 }

              // Initialise existing variables and vectors to -88 (-99 is defualt initialization for non-existan channels)
              WfmRun::VarInit(m_RunBase, a);

              // Get the branches for this event with respect to the file level and analysis stage
              if (!WfmRun::GetEvntBranches(entry, m_FFTShDef.at(a), m_filelevels.at(a), m_filetrack.at(a), a))
                 {
                  std::cout << __FUNCTION__ << " ERROR: Failed recovering branches for event: " << EvnNo  
                            << ", of file: " << m_RunBase->WfmBase::GetInRootFile(a)->GetName() << std::endl;
                  std::cout << __FUNCTION__ << " ERROR: Skipping...: ";
                  continue;
                 } 
              // Do track Combination here
              if (m_RunBase->WfmBase::GetTrackComb() && fdTrackFile.at(a) != -99)
                 {
                  if (m_RunBase->WfmBase::GetTrackPackage() == "Cory") WfmBase::CombineTrackCory(EvnNo-1, trckevntlast, fdTrackFile.at(a));
                  else if (m_RunBase->WfmBase::GetTrackPackage() == "EUTelescope") WfmBase::CombineTrackEUTel(EvnNo-1, trckevntlast, fdTrackFile.at(a));
                  else if (m_RunBase->WfmBase::GetTrackPackage() == "TBA") WfmBase::CombineTrackTBA(EvnNo-1, trckevntlast, fdTrackFile.at(a));
                  else {
                        std::cout << __FUNCTION__ << "        ERROR: Unsopported tracking package " << m_RunBase->WfmBase::GetTrackPackage() << std::endl;
                        return false;
                       }
                 }

              // Channel loop for current file
              for (unsigned int ich = 0; ich < (m_RunBase->WfmBase::GetChannels()).size(); ich++)
                  {
                   // Skip non existant channels in this file
                   if ((m_t.at(ich))->size() == 1 && (m_t.at(ich))->at(0) == -99) continue;
                   if (m_RunBase->WfmBase::GetAnaStage() == 1 && m_filelevels.at(a) < 2) FillFistLvl(m_RunBase, m_RunWaveform.at(ich), ich);
                   else if (m_RunBase->WfmBase::GetAnaStage() > 1 && m_filelevels.at(a) == 0) FillSecondLvl(m_RunBase, m_RunDUTCh.at(ich), ich);
                   else if (m_RunBase->WfmBase::GetAnaStage() > 1 && m_filelevels.at(a) == 1 ||
                            m_RunBase->WfmBase::GetAnaStage() == 2 && m_filelevels.at(a) == 2) FillSecFirst(m_RunBase, m_RunDUTCh.at(ich), ich);
                  }

              // Populate the vector for trigger frequency calculations when apropriate instrument, not processecing all events and outputing to a new file
              if (frqhist && (((m_RunBase->WfmBase::GetStartStopEvnt().first != 0 || m_RunBase->WfmBase::GetStartStopEvnt().second != 0) && m_cp == -99) || m_RunBase->WfmBase::GetInstrument() != Sampic))
                 {
                  if (m_RunBase->WfmBase::GetInstrument() == LabTXT || m_RunBase->WfmBase::GetInstrument() == TestBeamBin2) l_triggtimeVec.at(m_event - 1 - l_prosEvnt) = m_EvnTime;
                  else if (m_RunBase->WfmBase::GetInstrument() == TektronixScope) l_triggtimeVec.at(m_event - 1 - l_prosEvnt) = m_EvnTimeCh.at(0);
                 }

              if (!WriteBranches(a)) 
                 {
                  std::cout << __FUNCTION__ << " ERROR: Failed writing to file for event: " << EvnNo << ", aboding..." << std::endl;
                  return false;
                 }

             } // End of event loop for current file
          
          // Correctly populate trigger frequency histograms when apropriate instrument, not processecing all events and outputing to a new file
         if (frqhist && (((m_RunBase->WfmBase::GetStartStopEvnt().first != 0 || m_RunBase->WfmBase::GetStartStopEvnt().second != 0) && m_cp == -99) || m_RunBase->WfmBase::GetInstrument() != Sampic))
            {
             WfmBase::CalcTrigFr(l_triggtimeVec, m_trigDt, m_trigFr, (limEvnt.at(a)).second-(limEvnt.at(a)).first);
            }

        } // End of file loop
    
    // Closing and saving files
    m_RunBase->WfmBase::GetOutRootFile()->Write("", TObject::kOverwrite);
    if (WfmBase::GetTrackComb()) WfmBase::CloseTrackFiles();
    for (unsigned int k = 0; k < m_RunBase->WfmBase::GetInRootFiles().size(); k++) m_RunBase->WfmBase::GetInRootFile(k)->Close();
    if (m_cp != -99) m_RunBase->WfmBase::GetOutRootFile()->Close();

    return true;
}
// --------------------------------------------------------------------------------------------------------------
void WfmRun::FirstInit(std::pair<unsigned int, unsigned int> filelvl, unsigned int anastage)
{
  // Create new DUT channels
  if (anastage == 1) m_RunWaveform.resize((m_RunBase->WfmBase::GetChannels()).size(), nullptr);
  else if ((filelvl.first < 2 && anastage > 1) || (filelvl.first == 2 && (anastage == 2 || anastage == 3 || anastage == 5)) ||
           (filelvl.first == 3 && anastage == 3)) m_RunDUTCh.resize((m_RunBase->WfmBase::GetChannels()).size(), nullptr);

  std::vector<bool> foundPr[16];
  for (unsigned int h = 0; h < 16; h++) foundPr[h].resize((m_RunBase->WfmBase::GetChannels()).size(), false);
  unsigned int sizes[16];
  sizes[0] = (m_RunBase->WfmBase::GetDUTNames()).size();
  sizes[1] = (m_RunBase->WfmBase::GetDUTBoards()).size();
  sizes[2] = (m_RunBase->WfmBase::GetDUTransImps()).size();
  sizes[3] = (m_RunBase->WfmBase::GetDUTSecStages()).size();
  sizes[4] = (m_RunBase->WfmBase::GetDUTAmpGains()).size();
  sizes[5] = (m_RunBase->WfmBase::GetDUTriggs()).size();
  sizes[6] = (m_RunBase->WfmBase::GetDUTFracts()).size();
  sizes[7] = (m_RunBase->WfmBase::GetDUTCaps()).size();
  sizes[8] = (m_RunBase->WfmBase::GetDUTUpLims()).size();
  sizes[9] = (m_RunBase->WfmBase::GetDUTLowLims()).size();
  sizes[10] = (m_RunBase->WfmBase::GetDUTBinSizes()).size();
  sizes[11] = (m_RunBase->WfmBase::GetChVoltCuts()).size();
  sizes[12] = (m_RunBase->WfmBase::GetChNoiseCuts()).size();
  sizes[13] = (m_RunBase->WfmBase::GetChJitterCuts()).size();
  sizes[14] = (m_RunBase->WfmBase::GetChChargeCuts()).size();
  sizes[15] = (m_RunBase->WfmBase::GetDUTBinShapes()).size();

  for (unsigned int a = 0; a < (m_RunBase->WfmBase::GetChannels()).size(); a++)
      {
       if (m_RunWaveform.size() != 0) m_RunWaveform.at(a) = new WaveForm(m_RunBase);
       else {
             m_RunDUTCh.at(a) = new DUTChannel(m_channels.at(a), m_RunBase);
             for (unsigned int u = 0; u < 16; u++)
                 {
                  for (unsigned int w = 0; w < sizes[u]; w++) 
                      {
                       if (u == 0 && (int)(m_channels.at(a)) == ((m_RunBase->WfmBase::GetDUTNames()).at(w)).first) 
                          {
                           m_RunDUTCh.at(a)->DUTChannel::SetChName(((m_RunBase->WfmBase::GetDUTNames()).at(w)).second);
                           (foundPr[u]).at(a) = true;
                          }
                       if (u == 1 && (int)(m_channels.at(a)) == ((m_RunBase->WfmBase::GetDUTBoards()).at(w)).first)
                          {
                           m_RunDUTCh.at(a)->DUTChannel::SetChBoard(((m_RunBase->WfmBase::GetDUTBoards()).at(w)).second);
                           (foundPr[u]).at(a) = true;
                          }
                       if (u == 2 && (int)(m_channels.at(a)) == ((m_RunBase->WfmBase::GetDUTransImps()).at(w)).first)
                          {
                           m_RunDUTCh.at(a)->DUTChannel::SetChTransImp(((m_RunBase->WfmBase::GetDUTransImps()).at(w)).second);
                           (foundPr[u]).at(a) = true;
                          }
                       if (u == 3 && (int)(m_channels.at(a)) == ((m_RunBase->WfmBase::GetDUTSecStages()).at(w)).first)
                          {
                           m_RunDUTCh.at(a)->DUTChannel::SetSecStage(((m_RunBase->WfmBase::GetDUTSecStages()).at(w)).second);
                           (foundPr[u]).at(a) = true;
                          }                 
                       if (u == 4 && (int)(m_channels.at(a)) == ((m_RunBase->WfmBase::GetDUTAmpGains()).at(w)).first)
                          {
                           m_RunDUTCh.at(a)->DUTChannel::SetChAmpGain(((m_RunBase->WfmBase::GetDUTAmpGains()).at(w)).second);
                           (foundPr[u]).at(a) = true;
                          }
                       if (u == 5 && (int)(m_channels.at(a)) == ((m_RunBase->WfmBase::GetDUTriggs()).at(w)).first)
                          {
                           m_RunDUTCh.at(a)->DUTChannel::SetChTrigg(((m_RunBase->WfmBase::GetDUTriggs()).at(w)).second);
                           (foundPr[u]).at(a) = true;
                          }
                       if (u == 6 && (int)(m_channels.at(a)) == ((m_RunBase->WfmBase::GetDUTFracts()).at(w)).first)
                          {
                           m_RunDUTCh.at(a)->DUTChannel::SetChFract(((m_RunBase->WfmBase::GetDUTFracts()).at(w)).second);
                           (foundPr[u]).at(a) = true;
                          }
                       if (u == 7 && (int)(m_channels.at(a)) == ((m_RunBase->WfmBase::GetDUTCaps()).at(w)).first)
                          {
                           m_RunDUTCh.at(a)->DUTChannel::SetChCap(((m_RunBase->WfmBase::GetDUTCaps()).at(w)).second);
                           (foundPr[u]).at(a) = true;
                          }
                       if (u == 8 && (int)(m_channels.at(a)) == ((m_RunBase->WfmBase::GetDUTUpLims()).at(w)).first)
                          {
                           m_RunDUTCh.at(a)->DUTChannel::SetChUpLim(((m_RunBase->WfmBase::GetDUTUpLims()).at(w)).second);
                           (foundPr[u]).at(a) = true;
                          }
                       if (u == 9 && (int)(m_channels.at(a)) == ((m_RunBase->WfmBase::GetDUTLowLims()).at(w)).first)
                          {
                           m_RunDUTCh.at(a)->DUTChannel::SetChLowLim(((m_RunBase->WfmBase::GetDUTLowLims()).at(w)).second);
                           (foundPr[u]).at(a) = true;
                          }
                       if (u == 10 && (int)(m_channels.at(a)) == ((m_RunBase->WfmBase::GetDUTBinSizes()).at(w)).first)
                          {
                           m_RunDUTCh.at(a)->DUTChannel::SetChBinSz(((m_RunBase->WfmBase::GetDUTBinSizes()).at(w)).second);
                           (foundPr[u]).at(a) = true;
                          }
                       if (u == 11 && (int)(m_channels.at(a)) == ((m_RunBase->WfmBase::GetChVoltCuts()).at(w)).first)
                          {
                           m_RunDUTCh.at(a)->DUTChannel::SetVoltCut(((m_RunBase->WfmBase::GetChVoltCuts()).at(w)).second);
                           (foundPr[u]).at(a) = true;
                          }
                       if (u == 12 && (int)(m_channels.at(a)) == ((m_RunBase->WfmBase::GetChNoiseCuts()).at(w)).first)
                          {
                           m_RunDUTCh.at(a)->DUTChannel::SetNoiseCut(((m_RunBase->WfmBase::GetChNoiseCuts()).at(w)).second);
                           (foundPr[u]).at(a) = true;
                          }
                       if (u == 13 && (int)(m_channels.at(a)) == ((m_RunBase->WfmBase::GetChJitterCuts()).at(w)).first)
                          {
                           m_RunDUTCh.at(a)->DUTChannel::SetJitterCut(((m_RunBase->WfmBase::GetChJitterCuts()).at(w)).second);
                           (foundPr[u]).at(a) = true;
                          }
                       if (u == 14 && (int)(m_channels.at(a)) == ((m_RunBase->WfmBase::GetChChargeCuts()).at(w)).first)
                          {
                           m_RunDUTCh.at(a)->DUTChannel::SetChargeCut(((m_RunBase->WfmBase::GetChChargeCuts()).at(w)).second);
                           (foundPr[u]).at(a) = true;
                          }
                        if (u == 15 && (int)(m_channels.at(a)) == ((m_RunBase->WfmBase::GetDUTBinShapes()).at(w)).first)
                          {
                           m_RunDUTCh.at(a)->DUTChannel::SetChBinShape(((m_RunBase->WfmBase::GetDUTBinShapes()).at(w)).second);
                           (foundPr[u]).at(a) = true;
                          }
                      } 
                 }
            }
      }

  if (m_RunBase->WfmBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Intialized Channel ID properties for " 
                                                         << (m_RunBase->WfmBase::GetChannels()).size() << " channels" << std::endl;
  
  // Fill the properties of the channel without index
  if (m_RunWaveform.size() == 0)
     {
      std::vector<unsigned int> skip[16];
      for (unsigned int a = 0; a < (m_RunBase->WfmBase::GetChannels()).size(); a++)
          {
           for (unsigned int u = 0; u < 16; u++)
               {
                if (!(foundPr[u]).at(a))
                   {
                    for (unsigned int w = 0; w < sizes[u]; w++)
                        {
                         for (unsigned int tc = 0; tc < (skip[u]).size(); tc++) if (w == (skip[u]).at(tc)) continue;
                         if (u == 0 && ((m_RunBase->WfmBase::GetDUTNames()).at(w)).first == -99)
                            {
                             m_RunDUTCh.at(a)->DUTChannel::SetChName(((m_RunBase->WfmBase::GetDUTNames()).at(w)).second);
                             (skip[u]).push_back(w);
                            }
                         if (u == 1 && ((m_RunBase->WfmBase::GetDUTBoards()).at(w)).first == -99)
                            {
                             m_RunDUTCh.at(a)->DUTChannel::SetChBoard(((m_RunBase->WfmBase::GetDUTBoards()).at(w)).second);
                             (skip[u]).push_back(w);
                            }
                         if (u == 2 && ((m_RunBase->WfmBase::GetDUTransImps()).at(w)).first == -99)
                            {
                             m_RunDUTCh.at(a)->DUTChannel::SetChTransImp(((m_RunBase->WfmBase::GetDUTransImps()).at(w)).second);
                             (skip[u]).push_back(w);
                            }
                         if (u == 3 && ((m_RunBase->WfmBase::GetDUTSecStages()).at(w)).first == -99)
                            {
                             m_RunDUTCh.at(a)->DUTChannel::SetSecStage(((m_RunBase->WfmBase::GetDUTSecStages()).at(w)).second);
                             (skip[u]).push_back(w);
                            }
                         if (u == 4 && ((m_RunBase->WfmBase::GetDUTAmpGains()).at(w)).first == -99)
                            {
                             m_RunDUTCh.at(a)->DUTChannel::SetChAmpGain(((m_RunBase->WfmBase::GetDUTAmpGains()).at(w)).second);
                             (skip[u]).push_back(w);
                            }
                         if (u == 5 && ((m_RunBase->WfmBase::GetDUTriggs()).at(w)).first == -99)
                            {
                             m_RunDUTCh.at(a)->DUTChannel::SetChTrigg(((m_RunBase->WfmBase::GetDUTriggs()).at(w)).second);
                             (skip[u]).push_back(w);
                            }
                         if (u == 6 && ((m_RunBase->WfmBase::GetDUTFracts()).at(w)).first == -99)
                            {
                             m_RunDUTCh.at(a)->DUTChannel::SetChFract(((m_RunBase->WfmBase::GetDUTFracts()).at(w)).second);
                             (skip[u]).push_back(w);
                            }
                         if (u == 7 && ((m_RunBase->WfmBase::GetDUTCaps()).at(w)).first == -99)
                            {
                             m_RunDUTCh.at(a)->DUTChannel::SetChCap(((m_RunBase->WfmBase::GetDUTCaps()).at(w)).second);
                             (skip[u]).push_back(w);
                            }
                         if (u == 8 && ((m_RunBase->WfmBase::GetDUTUpLims()).at(w)).first == -99)
                            {
                             m_RunDUTCh.at(a)->DUTChannel::SetChUpLim(((m_RunBase->WfmBase::GetDUTUpLims()).at(w)).second);
                             (skip[u]).push_back(w);
                            }
                         if (u == 9 && ((m_RunBase->WfmBase::GetDUTLowLims()).at(w)).first == -99)
                            {
                             m_RunDUTCh.at(a)->DUTChannel::SetChLowLim(((m_RunBase->WfmBase::GetDUTLowLims()).at(w)).second);
                             (skip[u]).push_back(w);
                            }
                         if (u == 10 && ((m_RunBase->WfmBase::GetDUTBinSizes()).at(w)).first == -99)
                           {
                             m_RunDUTCh.at(a)->DUTChannel::SetChBinSz(((m_RunBase->WfmBase::GetDUTBinSizes()).at(w)).second);
                             (skip[u]).push_back(w);
                            }
                         if (u == 11 && ((m_RunBase->WfmBase::GetChVoltCuts()).at(w)).first == -99)
                            {
                             m_RunDUTCh.at(a)->DUTChannel::SetVoltCut(((m_RunBase->WfmBase::GetChVoltCuts()).at(w)).second);
                             (skip[u]).push_back(w);
                            }
                         if (u == 12 && ((m_RunBase->WfmBase::GetChNoiseCuts()).at(w)).first == -99)
                            {
                             m_RunDUTCh.at(a)->DUTChannel::SetNoiseCut(((m_RunBase->WfmBase::GetChNoiseCuts()).at(w)).second);
                             (skip[u]).push_back(w);
                            }
                         if (u == 13 && ((m_RunBase->WfmBase::GetChJitterCuts()).at(w)).first == -99)
                            {
                             m_RunDUTCh.at(a)->DUTChannel::SetJitterCut(((m_RunBase->WfmBase::GetChJitterCuts()).at(w)).second);
                             (skip[u]).push_back(w);
                            }
                         if (u == 14 && ((m_RunBase->WfmBase::GetChChargeCuts()).at(w)).first == -99)
                            {
                             m_RunDUTCh.at(a)->DUTChannel::SetChargeCut(((m_RunBase->WfmBase::GetChChargeCuts()).at(w)).second);
                             (skip[u]).push_back(w);
                            }
                         if (u == 14 && ((m_RunBase->WfmBase::GetDUTBinShapes()).at(w)).first == -99)
                            {
                             m_RunDUTCh.at(a)->DUTChannel::SetChBinShape(((m_RunBase->WfmBase::GetDUTBinShapes()).at(w)).second);
                             (skip[u]).push_back(w);
                            }
                        }
                   }
               }
          }
     }

  m_cp = -99;
  for (unsigned int k = 0; k < m_RunBase->WfmBase::GetInRootFiles().size(); k++)
      {
       if (m_RunBase->WfmBase::GetInRootFile(k) == m_RunBase->WfmBase::GetOutRootFile())
          {
           m_cp = k;
           break;
          }
      }

  // Create the output trees with respect to the analysis level
  if (m_cp == -99)
     {
      if (!(m_RunBase->WfmBase::GetOutRootFile()->IsOpen())) m_RunBase->WfmBase::GetOutRootFile()->ReOpen("UPDATE");
      m_RunBase->WfmBase::GetOutRootFile()->cd();
      // Create the output tree for the event by event information
      m_OutTree = new TTree("wfm", "Per-event quantities of charge collection measurements");
      if ((m_RunBase->WfmBase::GetFileLVL()).second > 1 || m_RunBase->WfmBase::GetAnaStage() > 1)
         {
          m_ChFitParam = new TTree("ChFitParam", "Channel parameters form fits");
          m_ChFitParam->SetAutoFlush(-10000000); // Set autoflush to 30 MB
          if (m_RunBase->WfmBase::GetHasWaveShape().second || m_RunBase->WfmBase::GetDoWaveShape())
             {
              m_ChShapeParam = new TTree("ChShapeParam", "Channel parameters form wwaveform shape");
              m_ChShapeParam->SetAutoFlush(-10000000); // Set autoflush to 30 MB
             }
          if ((m_RunBase->WfmBase::GetFileLVL()).second > 2 || m_RunBase->WfmBase::GetAnaStage() > 2)
             {
              m_TimeResParam = new TTree("TimeResParam", "Time Resolution quantities and vectors for all planes");
              m_TimeResParam->SetAutoFlush(-10000000); // Set autoflush to 30 MB
             }
         }
     }
  else {
        m_RunBase->WfmBase::GetInRootFile(m_cp)->ReOpen("UPDATE");
        m_OutTree = m_RunBase->WfmBase::GetRootTree(m_cp);
        if ((m_RunBase->WfmBase::GetFileLVL()).second > 1 && m_RunBase->WfmBase::GetAnaStage() > 1)
           {
            (m_RunBase->WfmBase::GetOutRootFile())->GetObject("ChFitParam", m_ChFitParam);
            m_ChFitParam->SetAutoFlush(-10000000); // Set autoflush to 10 MB
            if (m_RunBase->WfmBase::GetHasWaveShape().second && m_RunBase->WfmBase::GetDoWaveShape())
               {
                (m_RunBase->WfmBase::GetOutRootFile())->GetObject("ChShapeParam", m_ChShapeParam);
                m_ChShapeParam->SetAutoFlush(-10000000); // Set autoflush to 10 MB
               }
            if ((m_RunBase->WfmBase::GetFileLVL()).second > 2 && m_RunBase->WfmBase::GetAnaStage() > 2)
               {
                m_TimeResParam->SetAutoFlush(-10000000); // Set autoflush to 10 MB
                (m_RunBase->WfmBase::GetOutRootFile())->GetObject("TimeResParam", m_TimeResParam);
               }
           }
        else if ((m_RunBase->WfmBase::GetFileLVL()).second < 2 && m_RunBase->WfmBase::GetAnaStage() > 1)
                {
                 m_ChFitParam = new TTree("ChFitParam", "Channel parameters form fits");
                 m_ChFitParam->SetAutoFlush(-10000000); // Set autoflush to 10 MB      
                 if (!(m_RunBase->WfmBase::GetHasWaveShape().first) && m_RunBase->WfmBase::GetDoWaveShape())
                    {
                     m_ChShapeParam = new TTree("ChShapeParam", "Channel parameters form wwaveform shape");
                     m_ChShapeParam->SetAutoFlush(-10000000); // Set autoflush to 10 MB
                    }
                 if ((m_RunBase->WfmBase::GetFileLVL()).second < 3 && m_RunBase->WfmBase::GetAnaStage() > 2)
                    {
                     m_TimeResParam = new TTree("TimeResParam", "Time Resolution quantities and vectors for all planes");
                     m_TimeResParam->SetAutoFlush(-10000000); // Set autoflush to 10 MB
                    }
                }
       }
  m_RunBase->WfmBase::GetOutRootFile()->SetCompressionLevel(6);
  m_OutTree->SetAutoFlush(-30000000); // Set autoflush to 30 MB
}
// --------------------------------------------------------------------------------------------------------------
bool WfmRun::SetPlaneDQuants(std::vector<double> dQunat, std::vector<std::pair<std::pair<int, int>, double>> DQs)
{
    if ((m_RunBase->WfmBase::GetChannels()).size() == 0)
       {
        std::cout << __FUNCTION__ << " WARNING: No channels found!" << std::endl;
        return false;
       }
    else { 
          dQunat.clear();
          bool found = false;
          for (unsigned int k = 0; k < (m_RunBase->WfmBase::GetChannels()).size(); k++)
              {
               for (unsigned g = k + 1; g < (m_RunBase->WfmBase::GetChannels()).size(); g++)
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
void WfmRun::VarInit(WfmBase* lbase, unsigned int file)
{

  // Fit test variables to be used for all chaennels
  m_FitDiffChi2 = 0;
  m_FitDiffMag = 0;
  m_FitDiffMagErr = 0;
  m_FitDiffSigma = 0;
  m_FitDiffSigmaErr = 0;
  m_FitDiffTime = 0;

  // Only initialize the channels present in this file, the rest are already set to -99 by the function that sets branch addresses
  for (unsigned int k = 0; k < (lbase->WfmBase::GetChannels()).size(); k++)
      {
       if (b_t.at(k) != nullptr)
          {
           m_EvPol.at(k) = -88;
           m_EvCharge.at(k) = -88;
           m_EvIsSignal.at(k) = -88;
           m_EvInWindow.at(k) = -88;
           m_EvVoltStr.at(k) = -88;
           m_EvMaxindx.at(k) = -88;
           m_EvMaxTime.at(k) = -88;
           m_EvMaxVolt.at(k) = -88;
           m_EvMinIndx.at(k) = -88;
           m_EvMinTime.at(k) = -88;
           m_EvMinVolt.at(k) = -88;
           m_EvStrIndx.at(k) = -88;
           m_EvEndIndx.at(k) = -88;
           m_EvNoise.at(k) = -88;
           m_EvPedest.at(k) = -88;
           m_EvNoiseErr.at(k) = -88;
           m_EvPedestErr.at(k) = -88;
           m_EvRiseT.at(k) = -88;
           m_EvJitRiseSNR.at(k) = -88;
           m_EvTriggTime.at(k) = -88;
           m_EvDVDTMax.at(k) = -88;
           m_EvTriggToT.at(k) = -88;
           if ((m_FFTShDef.at(file)).at(0) || m_RunBase->WfmBase::GetDoFFT())
              {
               m_EvSignalFFT.at(k) = -88;
               m_EvNoiseFFT.at(k) = -88;
              }
           m_EvComplete.at(k) = -88;
           // Vector initialization
           (m_EvCFDTime.at(k)).clear();
           (m_EvCFDTime.at(k)).resize(19, -88);
           (m_EvJitNdVdT.at(k)).clear();
           (m_EvJitNdVdT.at(k)).resize(19, -88);
           (m_EvDVDTCFD.at(k)).clear();
           (m_EvDVDTCFD.at(k)).resize(19, -88);
           (m_EvCFDToT.at(k)).clear();
           (m_EvCFDToT.at(k)).resize(19, -88);
           // Vector for the CFD adjasted waveform values of each channel
           (m_EvVAdjCFD.at(k)).clear();
           (m_EvVAdjCFD.at(k)).push_back(-88);
          }
      }
}
// --------------------------------------------------------------------------------------------------------------
void WfmRun::FillFistLvl(WfmBase* lbase, WaveForm* wvmf, unsigned int ch)
{
  bool hlplogic[3] = { false, false, false };
  int hlpint = -99;

  if (lbase->WfmBase::GetInstrument() == Sampic) wvmf->WaveForm::SetSnRate(lbase->WfmBase::GetSRate(0));
  else wvmf->WaveForm::SetSnRate(WfmBase::GetSRate(ch));

  wvmf->WaveForm::InitializeWaveForm(1);
  wvmf->WaveForm::SetVoltage(m_w.at(ch));
  wvmf->WaveForm::SetTime(m_t.at(ch));
  wvmf->WaveForm::Calculate();
  if (lbase->WfmBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Event Added to channel: " << ch+1 << "/" 
                                                     << (lbase->WfmBase::GetChannels()).size() << "!" << std::endl;

  // Recover the main control variables
  hlplogic[0] = wvmf->WaveForm::GetSaturation();
  hlplogic[1] = wvmf->WaveForm::GetIsSignal();
  hlplogic[2] = wvmf->WaveForm::GetIsInWidow();
  hlpint = wvmf->WaveForm::GetPolarity();
  // Start populating the event variables
  m_EvVoltStr.at(ch) = hlplogic[0];
  m_EvPol.at(ch) = hlpint;
  // Populate max and min voltages
  if (!hlplogic[0])
     {
      m_EvMaxindx.at(ch) = wvmf->WaveForm::GetMaxIndx();
      m_EvMaxVolt.at(ch) = wvmf->WaveForm::GetMax();
     }
  else {
        m_EvMaxindx.at(ch) = -99;
        m_EvMaxVolt.at(ch) = -99;
       }
  m_EvMinIndx.at(ch) = wvmf->WaveForm::GetMinIndx();
  m_EvMinVolt.at(ch) = wvmf->WaveForm::GetMin();
  m_EvIsSignal.at(ch) = hlplogic[1];
  if (hlplogic[1])
     {
      m_EvInWindow.at(ch) = hlplogic[2];
      m_EvStrIndx.at(ch) = wvmf->WaveForm::GetStrIndx();
      m_EvEndIndx.at(ch) = wvmf->WaveForm::GetEndIndx();
      if (hlplogic[2])
         {
          m_FitDiffChi2 += wvmf->WaveForm::GetFitComp(0);
          m_FitDiffMag += wvmf->WaveForm::GetFitComp(1);
          m_FitDiffMagErr += wvmf->WaveForm::GetFitComp(2);
          m_FitDiffSigma += wvmf->WaveForm::GetFitComp(3);
          m_FitDiffSigmaErr += wvmf->WaveForm::GetFitComp(4);
          m_FitDiffTime += wvmf->WaveForm::GetFitComp(5);
          for (unsigned int tmb = 0; tmb < 6; tmb++)
              {
               for (unsigned int rest = 0; rest < 7; rest++)
                   {
                    if (wvmf->WaveForm::GetFitComp(0, true, tmb, rest) != std::numeric_limits<double>::max()) 
                       {
                        m_FitChi2Var->SetBinContent(tmb, rest, m_FitChi2Var->GetBinContent(tmb, rest) + wvmf->WaveForm::GetFitComp(0, true, tmb, rest));
                        m_NormCount[0]++;
                       }
                    if (wvmf->WaveForm::GetFitComp(1, true, tmb, rest) != -99)
                       { 
                        m_FitMagVar->SetBinContent(tmb, rest, m_FitMagVar->GetBinContent(tmb, rest) + wvmf->WaveForm::GetFitComp(1, true, tmb, rest));
                        m_NormCount[1]++;
                       }
                    if (wvmf->WaveForm::GetFitComp(2, true, tmb, rest) != -99) 
                       {
                        m_FitMagErrVar->SetBinContent(tmb, rest, m_FitMagErrVar->GetBinContent(tmb, rest) + wvmf->WaveForm::GetFitComp(2, true, tmb, rest));
                        m_NormCount[2]++;
                       }
                    if (wvmf->WaveForm::GetFitComp(3, true, tmb, rest) != -99) 
                       {
                        m_FitSigmaVar->SetBinContent(tmb, rest, m_FitSigmaVar->GetBinContent(tmb, rest) + wvmf->WaveForm::GetFitComp(3, true, tmb, rest));
                        m_NormCount[3]++;
                       }
                    if (wvmf->WaveForm::GetFitComp(4, true, tmb, rest) != -99) 
                       {
                        m_FitSigmaErrVar->SetBinContent(tmb, rest, m_FitSigmaErrVar->GetBinContent(tmb, rest) + wvmf->WaveForm::GetFitComp(1, true, tmb, rest));
                        m_NormCount[4]++;
                       }
                   }
              }
          m_EvNoise.at(ch) = wvmf->WaveForm::GetNoise();
          m_EvPedest.at(ch) = wvmf->WaveForm::GetPedestal();
          m_EvNoiseErr.at(ch) = wvmf->WaveForm::GetNoiseErr();
          m_EvPedestErr.at(ch) = wvmf->WaveForm::GetPedestalErr();
          if (!hlplogic[0]) m_EvMaxTime.at(ch) = wvmf->WaveForm::GetMaxTime();
          else m_EvMaxTime.at(ch) = - 99;
          m_EvMinTime.at(ch) = wvmf->WaveForm::GetMinTime();
          if (hlpint != 2 && !hlplogic[0])
             {
              m_EvCharge.at(ch) = wvmf->WaveForm::GetCharge();
              m_EvRiseT.at(ch) = wvmf->WaveForm::GetRiseTime();
              if (lbase->WfmBase::GetDoFFT())
                 {
                  m_EvSignalFFT.at(ch) = wvmf->WaveForm::GetSignalFFT();
                  m_EvNoiseFFT.at(ch) = wvmf->WaveForm::GetNoiseFFT();
                 }
              m_EvJitRiseSNR.at(ch) = wvmf->WaveForm::GetJitterRiseSNR();
              m_EvDVDTMax.at(ch) = wvmf->WaveForm::GetdVdTMax();
              for (unsigned int i = 0; i < 19; i++) (m_EvCFDTime.at(ch)).at(i) = wvmf->WaveForm::GetCFDTime(0.05*(i+1));
              if (wvmf->WaveForm::GetCFDTime() != -1)
                 {
                  m_EvVAdjCFD.at(ch) = wvmf->WaveForm::GetTimeAdjVolt();
                  for (unsigned int i = 0; i < 19; i++)
                      {
                       (m_EvDVDTCFD.at(ch)).at(i) = wvmf->WaveForm::GetdVdTCFD(0.05*(i+1));
                       (m_EvJitNdVdT.at(ch)).at(i) = wvmf->WaveForm::GetJitterNdVdT(0.05*(i+1));
                       }
                  m_EvTriggTime.at(ch) = wvmf->WaveForm::GetTriggTime();
                  if (wvmf->WaveForm::GetCFDToT() > 0.0)
                     {
                      for (unsigned int i = 0; i < 19; i++) (m_EvCFDToT.at(ch)).at(i) = wvmf->WaveForm::GetCFDToT(0.05*(i+1));
                      if (wvmf->WaveForm::GetTriggToT() > 0)
                         {
                          m_EvTriggToT.at(ch) = wvmf->WaveForm::GetTriggToT();
                          m_EvComplete.at(ch) = 6;
                         }
                      else {
                            WvfFillDefault(lbase, 6, ch);
                            m_EvComplete.at(ch) = 5;
                            if (lbase->GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Failed to calculate trigger CFD time: "
                                                                        << wvmf->WaveForm::GetTriggToT() << " at " << wvmf->WaveForm::GetTrigg() << std::endl;
                           }
                     }
                  else {
                        WvfFillDefault(lbase, 5, ch);
                        m_EvComplete.at(ch) = 4;
                        if (lbase->GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Failed to calculate ToT CFD time: "
                                                                 << wvmf->WaveForm::GetCFDToT() << " at " <<  wvmf->WaveForm::GetCFDfraction() << std::endl;
                       }
                 }
              else {
                    WvfFillDefault(lbase, 4, ch);
                    m_EvComplete.at(ch) = 3;
                    if (lbase->GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Failed to calculate CFD time"
                                                             << wvmf->WaveForm::GetCFDTime() << " for standard set value!" << std::endl;
                   }
             }
          else {
                WvfFillDefault(lbase, 3, ch);
                m_EvComplete.at(ch) = 2;
                if (lbase->GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Failed to calculate polarity [" << wvmf->WaveForm::GetPolarity()
                                                         << "] or signal saturated [" << hlplogic[0] << "] for cahnnel [" << (lbase->WfmBase::GetChannels()).at(ch)
                                                         << "] -> Stopping at Min/Max estimation!" << std::endl;
               }
         }
      else {
            WvfFillDefault(lbase, 2, ch);
            m_EvComplete.at(ch) = 1;
            if (lbase->GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: Signal out of window for channel ["
                                                     << (lbase->WfmBase::GetChannels()).at(ch) << "] -> Stopping at Min/Max estimation!" << std::endl;
           }
     }
  else {
        WvfFillDefault(lbase, 1, ch);
        m_EvComplete.at(ch) = 0;
        if (lbase->GetVerbosity() > 0) std::cout << __FUNCTION__ << " WARNING: No signal found in channel ["
                                                 << (lbase->WfmBase::GetChannels()).at(ch) << "] -> Stopping at Min/Max estimation!" << std::endl;
       }

  if (lbase->WfmBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Got event variables: " << ch+1 << "/"
                                                     << (lbase->WfmBase::GetChannels()).size() << "!" << std::endl;

}
// --------------------------------------------------------------------------------------------------------------
void WfmRun::WvfFillDefault(WfmBase* lbase, unsigned int completeness, unsigned int ch)
{
  if (completeness < 2)
     {
      m_EvInWindow.at(ch) = 1;
      m_EvStrIndx.at(ch) = -99;
      m_EvEndIndx.at(ch) = -99;
     }
  if (completeness < 3)
     {
      m_EvNoise.at(ch) = -99;
      m_EvPedest.at(ch) = -99;
      m_EvNoiseErr.at(ch) = -99;
      m_EvPedestErr.at(ch) = -99;
      m_EvMaxTime.at(ch) = -99;
      m_EvMinTime.at(ch) = -99;
     }
  if (completeness < 4)
     {
      m_EvCharge.at(ch) = -99;
      m_EvRiseT.at(ch) = -99;
      if (lbase->WfmBase::GetDoFFT())
         {
          m_EvSignalFFT.at(ch) = -99;
          m_EvNoiseFFT.at(ch) = -99;
         }
      m_EvJitRiseSNR.at(ch) = -99;
      m_EvDVDTMax.at(ch) = -99;
      for (unsigned int l = 0; l < 19; l++) (m_EvCFDTime.at(ch)).at(l) = -99;
     }
  if (completeness < 5)
     {
      m_EvVAdjCFD.at(ch) = std::vector<double>();
      m_EvTriggTime.at(ch) = - 99;
      for (unsigned int l = 0; l < 19; l++) (m_EvDVDTCFD.at(ch)).at(l) = -99;
      for (unsigned int l = 0; l < 19; l++) (m_EvJitNdVdT.at(ch)).at(l) = -99;
     }
  if (completeness < 6) for (unsigned int l = 0; l < 19; l++) (m_EvCFDToT.at(ch)).at(l) = -99;
  m_EvTriggToT.at(ch) = -99;
}
// --------------------------------------------------------------------------------------------------------------
void WfmRun::FillSecondLvl(WfmBase* lbase, DUTChannel* DUTCh, unsigned int ch)
{
  if (lbase->WfmBase::GetInstrument() == Sampic) DUTCh->DUTChannel::SetChRate(lbase->WfmBase::GetSRate(0));
  else DUTCh->DUTChannel::SetChRate(lbase->WfmBase::GetSRate(ch));
  DUTCh->DUTChannel::AppendEvent(m_w.at(ch), m_t.at(ch));
  if (lbase->WfmBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Event Added to channel: " << ch+1 << "/" 
                                                     << (lbase->WfmBase::GetChannels()).size() << "!" << std::endl;

  // Get Event Variables
  m_EvPol.at(ch) = DUTCh->DUTChannel::GetEvntPolarity();
  m_EvCharge.at(ch) = DUTCh->DUTChannel::GetEvntCharge();
  m_EvIsSignal.at(ch) = DUTCh->DUTChannel::GetEvntIsSignal();
  m_EvInWindow.at(ch) = DUTCh->DUTChannel::GetEvntInWindow();
  m_EvVoltStr.at(ch) = DUTCh->DUTChannel::GetEvntVoltStr();
  m_EvMaxindx.at(ch) = DUTCh->DUTChannel::GetEvntMaxindx();
  m_EvMaxTime.at(ch) = DUTCh->DUTChannel::GetEvntMaxTime();
  m_EvMaxVolt.at(ch) = DUTCh->DUTChannel::GetEvntMaxVolt();
  m_EvMinIndx.at(ch) = DUTCh->DUTChannel::GetEvntMinIndx();
  m_EvMinTime.at(ch) = DUTCh->DUTChannel::GetEvntMinTime();
  m_EvMinVolt.at(ch) = DUTCh->DUTChannel::GetEvntMinVolt();
  m_EvStrIndx.at(ch) = DUTCh->DUTChannel::GetEvntStrIndx();
  m_EvEndIndx.at(ch) = DUTCh->DUTChannel::GetEvntEndIndx();
  m_EvNoise.at(ch) = DUTCh->DUTChannel::GetEvntNoise();
  m_EvPedest.at(ch) = DUTCh->DUTChannel::GetEvntPedest();
  m_EvNoiseErr.at(ch) = DUTCh->DUTChannel::GetEvntNoiseErr();
  m_EvPedestErr.at(ch) = DUTCh->DUTChannel::GetEvntPedestErr();
  m_EvRiseT.at(ch) = DUTCh->DUTChannel::GetEvntRiseT();
  m_EvJitRiseSNR.at(ch) = DUTCh->DUTChannel::GetEvntJitRiseSNR();
  m_EvTriggTime.at(ch) = DUTCh->DUTChannel::GetEvntTriggTime();
  m_EvDVDTMax.at(ch) = DUTCh->DUTChannel::GetEvntDVDTMax();
  m_EvTriggToT.at(ch) = DUTCh->DUTChannel::GetEvntTriggToT();
  if (lbase->WfmBase::GetDoFFT())
     {
      m_EvSignalFFT.at(ch) = DUTCh->DUTChannel::GetEvntSignalFFT();
      m_EvNoiseFFT.at(ch) = DUTCh->DUTChannel::GetEvntNoiseFFT();
     }
  m_EvComplete.at(ch) = DUTCh->DUTChannel::GetEvntComplete();
  m_EvVAdjCFD.at(ch) = DUTCh->DUTChannel::GetEvntVAdjCFD();
  if (lbase->WfmBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Got event variables: " << ch+1 << "/" 
                                                     << (lbase->WfmBase::GetChannels()).size() << "!" << std::endl;
  // Get CFD values and populate the evctors
  for (unsigned int bin = 0; bin < 19; bin++)
      {
       (m_EvCFDTime.at(ch)).at(bin) = DUTCh->DUTChannel::GetEvntCFDTime(bin);
       (m_EvDVDTCFD.at(ch)).at(bin) = DUTCh->DUTChannel::GetEvntDVDTCFD(bin);
       (m_EvCFDToT.at(ch)).at(bin) = DUTCh->DUTChannel::GetEvntCFDToT(bin);
       (m_EvJitNdVdT.at(ch)).at(bin) = DUTCh->DUTChannel::GetEvntJitNdVdT(bin);
      }
  if (lbase->WfmBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Got CFD points: " << ch+1 << "/" 
                                                     << (lbase->WfmBase::GetChannels()).size() << "!" << std::endl;
 }
// --------------------------------------------------------------------------------------------------------------
bool WfmRun::FillSecFirst(WfmBase* lbase, DUTChannel* DUTCh, unsigned int ch)
{

  if (lbase->WfmBase::GetInstrument() == Sampic) DUTCh->DUTChannel::SetChRate(lbase->WfmBase::GetSRate(0));
  else DUTCh->DUTChannel::SetChRate(lbase->WfmBase::GetSRate(ch));
  if (lbase->WfmBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Event Added to channel: " << ch+1 << "/" 
                                                     << (lbase->WfmBase::GetChannels()).size() << "!" << std::endl;

  // Cuts implementation
  if (DUTCh->DUTChannel::GetNoiseCut() != -99 || DUTCh->DUTChannel::GetVoltCut() != -99 || 
      DUTCh->DUTChannel::GetChargetCut() != -99 || DUTCh->DUTChannel::GetJitterCut() != -99)
     {
      if (m_EvIsSignal.at(ch) && m_EvInWindow.at(ch))
         {
          if (DUTCh->DUTChannel::GetNoiseCut() != -99 && m_EvNoise.at(ch) > DUTCh->DUTChannel::GetNoiseCut()) return true;
          if (m_EvVoltStr.at(ch) != 1)
             {
              if (DUTCh->DUTChannel::GetVoltCut() != -99 && m_EvMaxVolt.at(ch) < DUTCh->DUTChannel::GetVoltCut()) return true;
              if (m_EvPol.at(ch) != 2)
                 {
                  if (DUTCh->DUTChannel::GetChargetCut() != -99 && m_EvCharge.at(ch) < DUTCh->DUTChannel::GetChargetCut()) return true;
                  if (DUTCh->DUTChannel::GetJitterCut() != -99 && m_EvJitRiseSNR.at(ch) > DUTCh->DUTChannel::GetJitterCut()) return true;
                 }
              else if (DUTCh->DUTChannel::GetChargetCut() != -99 || DUTCh->DUTChannel::GetJitterCut() != -99) return true;
             }
          else if (DUTCh->DUTChannel::GetVoltCut() != -99 || DUTCh->DUTChannel::GetChargetCut() != -99 || DUTCh->DUTChannel::GetJitterCut() != -99) return true;
         }
      else return true;
     }

  // Get Event Variables
  DUTCh->DUTChannel::SetEvntPolarity(m_EvPol.at(ch));
  DUTCh->DUTChannel::SetEvntCharge(m_EvCharge.at(ch));
  DUTCh->DUTChannel::SetEvntIsSignal(m_EvIsSignal.at(ch));
  DUTCh->DUTChannel::SetEvntInWindow(m_EvInWindow.at(ch));
  DUTCh->DUTChannel::SetEvntVoltStr(m_EvVoltStr.at(ch));
  DUTCh->DUTChannel::SetEvntMaxindx(m_EvMaxindx.at(ch));
  DUTCh->DUTChannel::SetEvntMaxTime(m_EvMaxTime.at(ch));
  DUTCh->DUTChannel::SetEvntMaxVolt(m_EvMaxVolt.at(ch));
  DUTCh->DUTChannel::SetEvntMinIndx(m_EvMinIndx.at(ch));
  DUTCh->DUTChannel::SetEvntMinTime(m_EvMinTime.at(ch));
  DUTCh->DUTChannel::SetEvntMinVolt(m_EvMinVolt.at(ch));
  DUTCh->DUTChannel::SetEvntStrIndx(m_EvStrIndx.at(ch));
  DUTCh->DUTChannel::SetEvntEndIndx(m_EvEndIndx.at(ch));
  DUTCh->DUTChannel::SetEvntNoise(m_EvNoise.at(ch));
  DUTCh->DUTChannel::SetEvntPedest(m_EvPedest.at(ch));
  DUTCh->DUTChannel::SetEvntNoiseErr(m_EvNoiseErr.at(ch));
  DUTCh->DUTChannel::SetEvntPedestErr(m_EvPedestErr.at(ch));
  DUTCh->DUTChannel::SetEvntRiseT(m_EvRiseT.at(ch));
  DUTCh->DUTChannel::SetEvntJitRiseSNR(m_EvJitRiseSNR.at(ch));
  DUTCh->DUTChannel::SetEvntTriggTime(m_EvTriggTime.at(ch));
  DUTCh->DUTChannel::SetEvntDVDTMax(m_EvDVDTMax.at(ch));
  DUTCh->DUTChannel::SetEvntTriggToT(m_EvTriggToT.at(ch));
  DUTCh->DUTChannel::SetEvntComplete(m_EvComplete.at(ch));
  DUTCh->DUTChannel::SetEvntVAdjCFD(m_EvVAdjCFD.at(ch));
  if (lbase->WfmBase::GetDoFFT())
     {
      DUTCh->DUTChannel::SetEvntSignalFFT(m_EvSignalFFT.at(ch));
      DUTCh->DUTChannel::SetEvntNoiseFFT(m_EvNoiseFFT.at(ch));
     }
    if (lbase->WfmBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Populated event variables: " << ch+1 << "/" 
                                                       << (lbase->WfmBase::GetChannels()).size() << "!" << std::endl;
  // Get CFD values and populate the evctors
  for (unsigned int bin = 0; bin < 19; bin++)
      {
       DUTCh->DUTChannel::SetEvntCFDTime((m_EvCFDTime.at(ch)).at(bin), bin);
       DUTCh->DUTChannel::SetEvntDVDTCFD((m_EvDVDTCFD.at(ch)).at(bin), bin);
       DUTCh->DUTChannel::SetEvntCFDToT((m_EvCFDToT.at(ch)).at(bin), bin);
       DUTCh->DUTChannel::SetEvntJitNdVdT((m_EvJitNdVdT.at(ch)).at(bin), bin);
      }
  if (lbase->WfmBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Populated CFD points: " << ch+1 << "/" 
                                                     << (lbase->WfmBase::GetChannels()).size() << "!" << std::endl;
  return true;
 }
// --------------------------------------------------------------------------------------------------------------
















    // Some general documentation
    // We have two vectors here, one that has all the files for processing and combining:       fileindx
    // one that has only the files on which tracking is to be combined with:                    l_fileindx_trck
    // A vector for files where tracking combination is not to be done exists:                  m_RunBase->WfmBase::GetExcludeTrackFiles()  
    // Planes for tracking, initialization already done:                                        std::vector<std::vector<unsigned int>> m_fileplanes;
    //                                                                                          std::vector<std::vector<unsigned int>> m_fileRefs;
    //                                                                                          std::vector<std::vector<unsigned int>> m_fileDUTs;
    //                                                                                          unsigned int m_nRecDUTs;
    //                                                                                          unsigned int m_nRecRefs;
    //                                                                                          unsigned int m_nRecPlanes;
    // std::vector<unsigned int> m_nplanes; // number of planes per file
    // std::vector<unsigned int> m_nRef; // number of References per file
    // std::vector<unsigned int> m_nDUT; // Number of DUTs per file
    // std::vector<unsigned int> m_DUTs; // flattened 2 dimentional vector of DUT numbers
    // std::vector<unsigned int> m_Refs; // flattened 2 dimentional vector of Reference numbers
    // std::vector<unsigned int> m_planes; // flattened 2 dimentional vector of Plane numbers
 /*
void WfmRun::SlaveTerminate()
{
    std::cout << setfill('-') << setw(89) << "-" << std::endl << setfill(' '); 
    m_OutTree->Write("", TObject::kOverwrite);
    TTree* ChFitParam = new TTree("ChFitParam", "Channel parameters form fits");
    if (m_RunBase->WfmBase::GetAnaStage() > 1) SecondStageAnalysis(ChFitParam); // Channel Properties
    if (m_RunBase->WfmBase::GetAnaStage() > 2) ThirdStageAnalysis(ChFitParam);  // CFD Maps and chi2 minimization
    if (m_RunBase->WfmBase::GetAnaStage() > 3) TrckEffAnalaysis(); // Tracking efficiency
    delete ChFitParam;
    delete m_OutTree;

    if (m_RunBase->WfmBase::GetInRootFile() != m_RunBase->WfmBase::GetOutRootFile()) 
       {
        m_RunBase->WfmBase::GetOutRootFile()->Close();
        m_RunBase->WfmBase::GetInRootFile()->cd();
       }
    m_RunBase->WfmBase::GetInRootFile()->Close();
    std::cout << __FUNCTION__ << " INFO: Analysis Finished successfully!" << std::endl;
}
// --------------------------------------------------------------------------------------------------------------
void WfmRun::Terminate()
{
    // for (unsigned int a = 0; a < m_nchan; a++) delete m_RunWaveform.at(a);
    std::cout << setfill('-') << setw(89) << "-" << std::endl << setfill(' ');
    //std::cout << __FUNCTION__ << " INFO: Cleaning up objeects..." << std::endl;
    //for (unsigned int a = 0; a < m_nchan; a++) delete m_RunDUTCh.at(a);
}
// --------------------------------------------------------------------------------------------------------------
void WfmRun::SecondStageAnalysis(TTree *ChFitParam)
{
  // Tree for Channel properties form Fits
  TTree *ChShapeParam = nullptr;
  ChFitParam->SetAutoFlush(10000); // Set autoflush to 30 MB

  std::vector<TList* > ChCanv;
  std::vector<TList* > CFDTimeCanv;
  std::vector<TList* > CFDToTCanv;
  std::vector<TList* > CFDdVdTCanv;

  // Channel Shape properties
  if (m_RunBase->WfmBase::GetDoWaveShape())
     {
      ChShapeParam = new TTree("ChShapeParam", "Channel parameters form wwaveform shape");
      ChShapeParam->SetAutoFlush(10000); // Set autoflush to 30 MB
     }

  for (unsigned int ich = 0; ich < m_nchan; ich++)
      {
       std::cout << __FUNCTION__ << " INFO: Updating Channel " << ich + 1 << " of " << m_nchan << " designated: " << m_channels.at(ich) << std::endl;
       m_RunDUTCh.at(ich)->DUTChannel::updateChProperties(m_RunBase->WfmBase::GetDoWaveShape(), m_OutTree);
       if (m_RunBase->WfmBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: After channel update: " << ich + 1 << "/" << m_nchan << std::endl;
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

       if (m_RunBase->WfmBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Iinitialized fit variables for channel: " << ich +1 << "/" << m_nchan << std::endl;

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
       if (m_RunBase->WfmBase::GetDoFFT())
          {
           ChFitSignalFFT = m_RunDUTCh.at(ich)->DUTChannel::GetChSignalFFT(2);
           ChFitNoiseFFT = m_RunDUTCh.at(ich)->DUTChannel::GetChNoiseFFT(2);
          }
       ChFitComplete = m_RunDUTCh.at(ich)->DUTChannel::GetChComplete();
       ChFitSoNR = m_RunDUTCh.at(ich)->DUTChannel::GetChSoNR(2);
       if (m_RunBase->WfmBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Recovered fit variables for channel: " << ich +1 << "/" << m_nchan << std::endl;

       // Get the fits and relevent chi2 to the Ntuple
       MaxIndxFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("MaxIndx");
       if (MaxIndxFtChi2 != -99 && MaxIndxFtChi2 != -1)
          {
           if (WfmBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("MaxIndx"), &h_MaxIndxFt.at(ich))) 
              {
               ChCanv.at(ich)->Add(h_MaxIndxFt.at(ich));
              }
          }
       MinIndxFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("MinIndx");
       if (MinIndxFtChi2 != -99 && MinIndxFtChi2 != -1)
          {
           if (WfmBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("MinIndx"), &h_MinIndxFt.at(ich)))
              {
               ChCanv.at(ich)->Add(h_MinIndxFt.at(ich));
              }
         }
       MaxVoltFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("MaxVolt");
       if (MaxVoltFtChi2 != -99 && MaxVoltFtChi2 != -1)
          {
           if (WfmBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("MaxVolt"), &h_MaxVoltFt.at(ich)))
              {
               ChCanv.at(ich)->Add(h_MaxVoltFt.at(ich));
              }
          }
       MinVoltFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("MinVotl");
       if (MinVoltFtChi2 != -99 && MinVoltFtChi2 != -1)
          {
           if (WfmBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("MinVotl"), &h_MinVoltFt.at(ich)))
              {
               ChCanv.at(ich)->Add(h_MinVoltFt.at(ich));
              }
          }
       StrIndxFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("Strindx");
       if (StrIndxFtChi2 != -99 && StrIndxFtChi2 != -1)
          {
           if (WfmBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("Strindx"), &h_StrIndxFt.at(ich)))
              {
               ChCanv.at(ich)->Add(h_StrIndxFt.at(ich));
              }
          }
       EndIndxFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("EndIndx");
       if (EndIndxFtChi2 != -99 && EndIndxFtChi2 != -1)
          {
           if (WfmBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("EndIndx"), &h_EndIndxFt.at(ich)))
              {
               ChCanv.at(ich)->Add(h_EndIndxFt.at(ich));
              }
          }
       NoiseFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("Noise");
       if (NoiseFtChi2 != -99 && NoiseFtChi2 != -1)
          {
           if (WfmBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("Noise"), &h_NoiseFt.at(ich)))
              {
               ChCanv.at(ich)->Add(h_NoiseFt.at(ich));
              }
          }
       NoiseErrFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("NoisErr");
       if (NoiseErrFtChi2 != -99 && NoiseErrFtChi2 != -1)
          {
           if (WfmBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("NoisErr"), &h_NoiseErrFt.at(ich)))
              {
               ChCanv.at(ich)->Add(h_NoiseErrFt.at(ich));
              }
          }
       PedestFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("Pedestl");
       if (PedestFtChi2 != -99 && PedestFtChi2 != -1)
          {
           if (WfmBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("Pedestl"), &h_PedestFt.at(ich)))
              {
               ChCanv.at(ich)->Add(h_PedestFt.at(ich));
              }
          }
       PedestErrFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("PdslErr");
       if (PedestErrFtChi2 != -99 && PedestErrFtChi2 != -1)
          {
           if (WfmBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("PdslErr"), &h_NPedestErrFt.at(ich)))
              {
               ChCanv.at(ich)->Add(h_NPedestErrFt.at(ich));
              }
          }
       MaxTimeFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("MaxTime");
       if (MaxTimeFtChi2 != -99 && MaxTimeFtChi2 != -1)
          {
           if (WfmBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("MaxTime"), &h_MaxTimeFt.at(ich)))
              {
               ChCanv.at(ich)->Add(h_MaxTimeFt.at(ich));
              }
          }
       MinTimeFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("MinTime");
       if (MinTimeFtChi2 != -99 && MinTimeFtChi2 != -1)
          {
           if (WfmBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("MinTime"), &h_MinTimeFt.at(ich)))
              {
               ChCanv.at(ich)->Add(h_MinTimeFt.at(ich));
              }
          }
       ChargeFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("Charge");
       if (ChargeFtChi2 != -99 && ChargeFtChi2 != -1)
          {
           if (WfmBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("Charge"), &h_ChargeFt.at(ich)))
              {
               ChCanv.at(ich)->Add(h_ChargeFt.at(ich));
              }
          }
       RiseTimeFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("RiseTime");
       if (RiseTimeFtChi2 != -99 && RiseTimeFtChi2 != -1)
          {
           if (WfmBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("RiseTime"), &h_RiseTimeFt.at(ich)))
              {
               ChCanv.at(ich)->Add(h_RiseTimeFt.at(ich));
              }
          }
       TriggTimeFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("TriggTime");
       if (TriggTimeFtChi2 != -99 && TriggTimeFtChi2 != -1)
          {
           if (WfmBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("TriggTime"), &h_TriggTimeFt.at(ich)))
              {
               ChCanv.at(ich)->Add(h_TriggTimeFt.at(ich));
              }
          }
       DVDTMaxFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("DVDTMax");
       if (DVDTMaxFtChi2 != -99 && DVDTMaxFtChi2 != -1)
          {
           if (WfmBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("DVDTMax"), &h_DVDTMaxFt.at(ich)))
              {
               ChCanv.at(ich)->Add(h_DVDTMaxFt.at(ich));
              }
          }
       TriggToTFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("TrigToT");
       if (TriggToTFtChi2 != -99 && TriggToTFtChi2 != -1)
          {
           if (WfmBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("TrigToT"), &h_TriggToTFt.at(ich)))
              {
               ChCanv.at(ich)->Add(h_TriggToTFt.at(ich));
              }
          }

       if (m_RunBase->WfmBase::GetDoFFT())
          {
           SignalFFTFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("SigFFT");
           if (SignalFFTFtChi2 != -99 && SignalFFTFtChi2 != -1)
              {
               if (WfmBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("SigFFT"), &h_SignalFFTFt.at(ich)))
                  {
                   ChCanv.at(ich)->Add(h_SignalFFTFt.at(ich));
                  }
              }
           NoiseFFTFtChi2 = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("NoiseFFT");
           if (NoiseFFTFtChi2 != -99 && NoiseFFTFtChi2 != -1)
              {
               if (WfmBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("NoiseFFT"), &h_NoiseFFTFt.at(ich)))
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
                if (WfmBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("CFDTime", h), &(h_CFDTimeFt.at(ich)).at(h)))
                   {
                    sprintf(title, "CFDTime%02u-%02u%%", m_channels.at(ich), (h*5+5));
                    (h_CFDTimeFt.at(ich)).at(h)->SetName((const char*)title);
                    CFDTimeCanv.at(ich)->Add((h_CFDTimeFt.at(ich)).at(h));
                   }
               }
            DVDTCFDFtChi2.at(h) = m_RunDUTCh.at(ich)->DUTChannel::GetFitChi2("DVDTCFD", h);
            if (DVDTCFDFtChi2.at(h) != -99 && DVDTCFDFtChi2.at(h) != -1)
               {
                if (WfmBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("DVDTCFD", h), &(h_DVDTCFDFt.at(ich)).at(h)))
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
                if (WfmBase::PrintFitInfo(m_RunDUTCh.at(ich)->DUTChannel::GetFit("CFDToT", h), &(h_CFDToTFt.at(ich)).at(h)))
                   {
                    memset(title, '0', sizeof(title));
                    sprintf(title, "CFDToT%02u-%02u%%", m_channels.at(ich), (h * 5 + 5));
                    (h_CFDToTFt.at(ich)).at(h)->SetName((const char*)title);
                    CFDToTCanv.at(ich)->Add((h_CFDToTFt.at(ich)).at(h));
                    if (h == 18) (h_CFDToTFt.at(ich)).at(h)->Draw();
                   }
               }
           }
       if (m_RunBase->WfmBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Recovered fits for channel: " << ich + 1 << "/" << m_nchan << std::endl;
         
       ChFitParam->Fill();

       // Mean pulseShape properties if shape is calculated
       if (m_RunBase->WfmBase::GetDoWaveShape())
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
           if (m_RunBase->WfmBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Initialized waveshape variables for channel: " << ich + 1 << "/" << m_nchan << std::endl;

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
           if (m_RunBase->WfmBase::GetDoFFT())
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
           if (m_RunBase->WfmBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Recovered waveshape variables for channel: " 
                                                                      << ich +1 << "/" << m_nchan << std::endl;
           ChShapeParam->Fill();
          }

       if (m_RunBase->WfmBase::GetVerbosity() == 0) m_RunDUTCh.at(ich)->DUTChannel::ChDump();

       // Create folder for each channel and write corresponding histos
       m_RunBase->WfmBase::GetOutRootFile()->mkdir(Form("ChHistos%02u", m_channels.at(ich)));
       m_RunBase->WfmBase::GetOutRootFile()->cd(Form("ChHistos%02u", m_channels.at(ich)));
       ChCanv.at(ich)->Write();

       m_RunBase->WfmBase::GetOutRootFile()->mkdir(Form("ChHistos%02u/CFDTimeHists%02u", m_channels.at(ich), m_channels.at(ich)));
       m_RunBase->WfmBase::GetOutRootFile()->cd(Form("ChHistos%02u/CFDTimeHists%02u", m_channels.at(ich), m_channels.at(ich)));
       CFDTimeCanv.at(ich)->Write();
       m_RunBase->WfmBase::GetOutRootFile()->mkdir(Form("ChHistos%02u/CFDdVdTHists%02u", m_channels.at(ich), m_channels.at(ich)));
       m_RunBase->WfmBase::GetOutRootFile()->cd(Form("ChHistos%02u/CFDdVdTHists%02u", m_channels.at(ich), m_channels.at(ich)));
       CFDdVdTCanv.at(ich)->Write();
       m_RunBase->WfmBase::GetOutRootFile()->mkdir(Form("ChHistos%02u/CFDToTHists%02u", m_channels.at(ich), m_channels.at(ich)));
       m_RunBase->WfmBase::GetOutRootFile()->cd(Form("ChHistos%02u/CFDToTHists%02u", m_channels.at(ich), m_channels.at(ich)));
       CFDToTCanv.at(ich)->Write();
       m_RunBase->WfmBase::GetOutRootFile()->cd();
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

  if (m_RunBase->WfmBase::GetDoWaveShape())
     {
      ChShapeParam->Write();
      delete ChShapeParam;
     }
}
// --------------------------------------------------------------------------------------------------------------
bool WfmRun::ThirdStageProcess()
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
void WfmRun::ThirdStageAnalysis(TTree* tree)
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
              if (m_RunBase->WfmBase::GetInstrument() == TestBeamBin1 || m_RunBase->WfmBase::GetInstrument() == TestBeamBin2 || m_RunBase->WfmBase::GetInstrument() == LeCroyWRBin)
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
                        if (m_RunBase->WfmBase::GetFitMethode() == "rootInt")
                           {
                            qual = m_RunBase->WfmBase::IterativeFit(&(m_EvTmDiff.at(a)).at(b), mean, sigma, 
                                                                    (h_TmDiffCFD.at(a)).at(b), 
                                                                    (TmDiffChi2.at(a)).at(b), 
                                                                    "GaussInt", std::make_pair(-1, -1), false);
                           }
                        else if (m_RunBase->WfmBase::GetFitMethode() == "root" || m_RunBase->WfmBase::GetFitMethode() == "RooFit")
                                {
                                 std::cout << a << " " << b << " " << m_EvTmDiff.size() << " " << (m_EvTmDiff.at(a)).size() << " " << h_TmDiffCFD.size() << " " << (h_TmDiffCFD.at(a)).size() << " " << TmDiffChi2.size() << " " << (TmDiffChi2.at(a)).size() << std::endl;
                                 qual = m_RunBase->WfmBase::IterativeFit(&(m_EvTmDiff.at(a)).at(b), mean, sigma, 
                                                                          (h_TmDiffCFD.at(a)).at(b), 
                                                                          (TmDiffChi2.at(a)).at(b), 
                                                                          "Gauss", std::make_pair(-1, -1), false);
                                }
                        if (qual == 0) mean.second = sigma.first;
                        else {
                              mean.first = m_RunBase->WfmBase::Mean(&(m_EvTmDiff.at(a)).at(b));
                              mean.second = m_RunBase->WfmBase::Stdev(&(m_EvTmDiff.at(a)).at(b));
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
                            if (WfmBase::PrintFitInfo((h_TmDiffCFD.at(a)).at(b), &(h_TmDiffCFDCanv.at(a)).at(b)))
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
              m_RunBase->WfmBase::GetOutRootFile()->mkdir(Form("ChHistos%02u/TimeDiff%02u-%02u", m_channels.at(ich), m_channels.at(ich), m_channels.at(gp)));
              m_RunBase->WfmBase::GetOutRootFile()->cd(Form("ChHistos%02u/TimeDiff%02u-%02u", m_channels.at(ich), m_channels.at(ich), m_channels.at(gp)));
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
*/