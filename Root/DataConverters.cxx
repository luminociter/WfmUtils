/*
* DataConverter.cxx
*
*
*      Author: Gkougkousis Evangelos - Leonidas
*              egkougko@cern.ch
*               IFAE-BARCELONA
*/
 
#include "../LGADUtils/LGADBase.h"
#include <time.h>
#include <stdio.h>

void LGADBase::SetInDataNames(TString DataDir, TString DataName, TString ext)
{
    m_datadir.Clear();
    m_ext.Clear();
    m_datanames.clear();

    // Check if the directory name has extention and if yes decompose accordingly
    std::size_t found = DataDir.Last('.');
    // If we found an extention (dot not at the begining, not at the end and extetniton at leat three charaters)
    if (found != std::string::npos && found != 0 && DataDir(found+1, DataDir.Length()-(found+1)).Length() > 2)
       {
        m_ext = DataDir(found, DataDir.Length()-found);
        found = DataDir.Last('/');
        m_datadir = DataDir(0, found+1);
        m_datanames.push_back(DataDir(found+1, DataDir.Length()-(m_ext.Length()+found+1)));
       }
    else {
          TString sLast = DataDir[DataDir.Length() - 1];
          if (!sLast.Contains("/")) DataDir = DataDir + "/";
          m_datadir = DataDir;
         }

    // Check the data file name if populated and remove extention if it exists only if extention not found in 
    // directory name. If extention was found in directory name, everything below is ignored.
    if (m_ext.Length() == 0 && DataName.Length() != 0)
       {
        found = DataName.Last('.'); 
        if (found != std::string::npos && DataName(found + 1, DataName.Length() - (found + 1)).Length() > 2)
           {
            m_ext = DataName(found, DataName.Length() - found);
            m_datanames.push_back(DataName(0, DataName.Length() - found + m_ext.Length()));
           }
        else {
              m_datanames.push_back(DataName);
              m_ext = ext;
              std::cout << __FUNCTION__ << " Scanning for input files in folder " << m_datadir << " ...";
              std::vector<std::string> temp;
              if (m_ext != "") // Add all files with that extention and filename in the directory and remove their extention
                 {
                  temp = ListFileNames(m_datadir + m_datanames.back(), m_ext);
                  m_datanames.pop_back();
                  std::cout << temp.size() << " found!" << std::endl;
                  for (unsigned int i = 0; i < temp.size(); i++)
                      {
                       found = (temp.at(i)).find_last_of(".");
                       m_datanames.at(i) = (TString)((temp.at(i)).erase(found, string::npos));
                      }
                 }
              else { // Add all files with filename and any extention without removeing the extention
                    temp = ListFileNames(m_datadir + m_datanames.back(), nullptr);
                    m_datanames.pop_back();
                    std::cout << temp.size() << " found!" << std::endl;
                    for (unsigned int i = 0; i < temp.size(); i++) m_datanames.push_back((TString)(temp.at(i)));
                   }  
             }
       }
    else if (m_ext.Length() == 0 && DataName.Length() == 0) // case were no filename is provided
            { 
             std::cout << __FUNCTION__ << " Scanning for input files in folder " << m_datadir << " ...";
             std::vector<std::string> temp;
             if (ext != "") // Add all files with this extention in directory and remove their extention
                {
                 char arr[5];
                 strcpy(arr, ext.Data());
                 m_datanames.resize(CountFiles(m_datadir, arr));
                 std::cout << m_datanames.size() << " found!" << std::endl;
                 temp = ListFileNames(m_datadir, arr);
                 for (unsigned int i = 0; i < temp.size(); i++) 
                     {
                      found = (temp.at(i)).find_last_of(".");
                      m_datanames.at(i) = (TString)((temp.at(i)).erase(found, std::string::npos));
                     }
                }
             else { // Add all files present in the directoy without removing their extentions
                   m_datanames.resize(CountFiles(m_datadir, nullptr));
                   std::cout << m_datanames.size() << " found!" << std::endl;
                   temp = ListFileNames(m_datadir, nullptr);
                   for (unsigned int i = 0; i < temp.size(); i++) m_datanames.at(i) = (TString)(temp.at(i));
                  }
            }
}
// --------------------------------------------------------------------------------------------------------------
void LGADBase::SetOutDataNames(TString DataDir, TString DataName)
{
    TString sLast = DataDir[DataDir.Length() - 1];
    if (!sLast.Contains("/")) DataDir = DataDir + "/";
    if (DataName.Contains(".root")) m_ofname = DataName(0, DataName.Length() - 5);
    else m_ofname = DataName;
    m_ofdir = DataDir;
}
// --------------------------------------------------------------------------------------------------------------
bool LGADBase::ConvertData()
{
    if (m_datadir.IsNull()) { std::cout << __FUNCTION__ << " ERROR: Input data directory not set! Abording..." << std::endl; return false; }
    if (m_datanames.size() == 0) { std::cout << __FUNCTION__ << " ERROR: No inpput data files found! Abording..." << std::endl; return false; }

    if (m_ofname.IsNull())
       {
        std::cout << "--> Warning: Output Ntuple name not set! Using input data filename: ";
        if (m_instrument == LabTXT || m_instrument == TektronixScope || m_instrument == Sampic)
           {
            m_ofname = m_datadir;
            TString sLast2 = m_ofname[m_ofname.Length() - 1];
            if (sLast2.Contains("/")) m_ofname = m_ofname(0, m_ofname.Length() - 1);
            m_ofname = m_ofname(m_ofname.Last('/') + 1, m_ofname.Length() - m_ofname.Last('/') - 1);
           }
        else m_ofname = m_datanames.at(0);
        std::cout << m_ofname << std::endl;
       }
    if (m_ofdir.IsNull())
       {
        std::cout << "--> Warning: Output Ntuple location not set! Using input data diretory" << std::endl;
        m_ofdir = m_datadir;
       }
    if (m_TrackComb && (m_Trackdatadir.IsNull() || m_Trackdataname.size()==0))
       { 
        std::cout << "--> Warning: Tracking input file not defined, switching to shadow test beam mode" << std::endl;
        m_TrackComb = false;
       }

    m_nchan = 0;
    m_channels.clear();
    m_convert = false;
    m_event = 0;
    if (m_instrument == LeCroyWRBin) m_convert = WriteLecroyBinary(m_datadir.Data(), m_datanames, m_ext.Data(), m_evnt1, m_evnt2);
    else if (m_instrument == TestBeamBin1 || m_instrument == TestBeamBin2 || m_instrument == TestBeamBin) m_convert = WriteTestBeamBinary(m_datadir.Data(), m_datanames, m_ext.Data(), m_evnt1, m_evnt2);
    else if (m_instrument == Sampic) m_convert = WriteSampic(m_datadir.Data(), m_datanames, m_ext.Data(), m_evnt1, m_evnt2);
    else if (m_instrument == LabTXT) m_convert = WriteLabTXT(m_datadir.Data(), m_datanames, m_ext.Data(), m_evnt1, m_evnt2);
    else if (m_instrument == TektronixScope) m_convert = WriteTectronixTXT(m_datadir.Data(), m_datanames, m_ext.Data(), m_evnt1, m_evnt2);
    else {
          std::cout << __FUNCTION__ << " ERROR: Unrecongnised instrument!" << std::endl;
          return false;
         }
    m_datanames.clear();
    if (m_TrackComb) m_TrackComb = false;
    return true;
}
// --------------------------------------------------------------------------------------------------------------
// Functon to read Sampic data. Because of the way sampic outputs events in a single consecutive
// .dat file, start and stop event numbers are of no conseuence here. Functon defines automatically
// the number of active channels and their respective position. Only accepts evetns that are within
// one buffer window. Time corrections for bin shifting is applied tking into account the ordered cell0
// time. All information are included in the output root file.
bool LGADBase::WriteSampic(const char* dir, std::vector<TString> names, const char* ext, unsigned int evt1, unsigned int evt2)
{
  // Temporary variables to be used in the function
  TBits* b = new TBits;
  std::istringstream iss(std::ios_base::app | std::ios_base::in | std::ios_base::out);
  char ifname[2048];
  unsigned int prv = 0;
  unsigned int indx = 0;
  std::string line;
  std::ifstream file;
  bool fnd = false;
  double tri = -99.;
  int channel = -1;
  double phystime;
  double ordrtime;

  std::vector<unsigned int> l_nfiles;  // Files per run
  std::vector<unsigned int> l_runIndx; // indexes of each run for synchronization with tracking files
  std::vector<TString> l_Runs;         // Base name of the run file
  std::vector<TString> l_ext;          // Extentions of the files in the directory
  unsigned int nRuns = RunCount(names, l_nfiles, l_Runs, l_ext, l_runIndx, ext);

  // Check that the number of extention in the extention vector is equal to the maximum number of
  // files in a run, as each run creates multiple files iterating on the extention prefix
  if (l_ext.size() != *std::max_element(l_nfiles.begin(), l_nfiles.end()))
     {
      std::cout << __FUNCTION__ << " ERROR: Inconsistent number of files with respect to expected number of runs in: " << dir << std::endl;
      std::cout << __FUNCTION__ << "        Abording......." << std::endl;
      return false;
     }

  std::vector<int> l_TotHits;                       // Total number of events prer run
  std::vector<Long64_t> l_srate;                    // Sampling rate per run
  std::vector<unsigned int> l_nChans;               // numer of channels per Run 
  std::vector<std::vector<unsigned int>> l_chName;  // channels names for each run
  unsigned long int l_totalevnts = 0;                 // Total available events from all runs (non reconstituted)
  prv = 0;
  for (unsigned int f = 0; f < nRuns; f++)
      {
       // Locate last file of the run    
       if (f > 0) prv += l_nfiles.at(f-1);
       indx = prv + l_nfiles.at(f)-1;
fileprob1:
       memset(&ifname[0], 0, sizeof(ifname));
       if (strcmp(ext, "") == 0) sprintf(ifname, "%s%s", dir, (names.at(indx)).Data());
       else {
             if (indx == prv) sprintf(ifname, "%s%s.%s", dir, (names.at(indx)).Data(), ext);
             else sprintf(ifname, "%s%s.%s_%04d", dir, (names.at(indx)).Data(), ext, l_nfiles.at(f)-1);
            }
       file.open(ifname);
       if (!file)
          {
fileprob2:
           std::cout << __FUNCTION__ << " ERROR: Unable to find input file: " << ifname;
           names.erase(names.begin() + indx);
           if (indx > prv) 
              {
               std::cout << ", skipping..." << std::endl;               
               l_nfiles.at(f) = l_nfiles.at(f) - 1;
               if (l_ext.size() != *std::max_element(l_nfiles.begin(), l_nfiles.end())) 
                   l_ext.resize(*std::max_element(l_nfiles.begin(), l_nfiles.end()));
               indx--;
               goto fileprob1;
              }
           else {
                 std::cout << ". This was the first file, run will be ignored!" << std::endl;
                 l_nfiles.erase(l_nfiles.begin() + f);
                 l_Runs.erase(l_Runs.begin() + f);
                 l_runIndx.erase(l_runIndx.begin() + f);
                 if (l_ext.size() != *std::max_element(l_nfiles.begin(), l_nfiles.end()))
                     l_ext.resize(*std::max_element(l_nfiles.begin(), l_nfiles.end()));
                 nRuns--; f--;
                 continue;
                }
          }
       if (!(file.is_open())) goto fileprob2; // Start reading from the back until you get the 3rd line from last
       else {
             l_TotHits.push_back(0);
             l_srate.push_back(0);
             l_nChans.push_back(0);
             l_chName.push_back(std::vector<unsigned int>());
             file.seekg(-2, file.end);
             prv = 0; // reuse integer
#ifdef _WIN32
             unsigned int ps = 4;
#else
             unsigned int ps = 3;
#endif
             while (prv < ps)
                   {
                    char cdfcv;
                    file.get(cdfcv);
                    if ((int)file.tellg() <= 1)
                       {
                        file.seekg(0);
                        prv = 4;
                       }
                    else if (cdfcv == '\n') 
                            {
                             file.seekg(-3, ios_base::cur); 
                             prv++;
                            }
                    else file.seekg(-2, ios_base::cur);
                   }
             file.seekg(3, ios_base::cur);
             getline(file, line);
             line = reduce(line, "", "=");
             sscanf(line.c_str(), "%*s %u %*s %*lf", &l_TotHits.back());
             // Read form the start of the file to get channels and properties
             file.seekg(0, std::ios::beg);
             indx = 0;  // reuse integer
             // Find the total number of channels in the system
             for (unsigned int h = 0; h < 2; h++) { line.clear(); getline(file, line); }
             line = reduce(line, "", "===");
             line = LGADBase::trim(line, " ");
             sscanf(line.c_str(), "%*s %*s %*f %*s %*s %*s %*s %*s %d", &indx);
             // Find the number of active channels by converting the hex register value to a decimal
             line.clear();
             getline(file, line);
             line = reduce(line, "", "===");
             line = LGADBase::trim(line, " ");
             sscanf(line.c_str(), "%*s %*s %llu %*s %*s %*s %*s %x", &l_srate.back(), &l_nChans.back());
             b->Set(indx, &l_nChans.back());
             l_nChans.back() = 0;
             for (unsigned int i = 0; i < indx; i++) if (b->TestBitNumber(i)) { (l_chName.back()).push_back(i + 1); l_nChans.back()++; }
             file.close();
             l_totalevnts += l_TotHits.back();
            }
      }

  // Fix the number of channels
  if (!FixBaseQuant(l_chName, m_channels, m_nchan)) return false;

  // Create the Ntuple taking into accoutn the number of channels in all runs to be combined, also setting the vector sizes for later on
  if (!CreateOutputFile(m_ofdir, m_ofname, m_channels))
     {
      std::cout << __FUNCTION__ << " ERROR: Output Ntuple creation failed: " << m_ofname << std::endl;
      return false;
     }

  // Tracking Mode related actions:
  // - Reduce the number of events by one to account for the initial TLU isseud trigger at the start of each run
  // - Check the oscilloscope delay vector and populate it correctly if needed
  // - Open the required data file and get the collection of trees
  // - Add the required parts to the Ntuple
  if (LGADBase::GetTrackComb())
     {
      std::vector<unsigned int> l_runIndx_old = l_runIndx;
      if (!OpenTrackFile(m_Trackdatadir, m_Trackdataname, l_runIndx))
         {
          std::cout << __FUNCTION__ << " ERROR: Failed oppening tacking files in directory: " << m_Trackdatadir << std::endl;
          return false;
         }
      else {
            if (!AddTrckBranches(m_trees.at(0), m_trkpck))
               {
                std::cout << __FUNCTION__ << " Warning: Failed adding the tracking part on the output file: " << m_ofname << std::endl;
                std::cout << "                                       Switching to non-trackimg mode..." << std::endl;
                LGADBase::SetTrackComb(false);
               }
           }
      if (l_runIndx.size() != l_runIndx_old.size())
         {
          std::vector<unsigned int> missing_indxs;
          for (unsigned int we = 0; we < l_runIndx_old.size(); we++)
              {
               fnd = false;
               for (unsigned int yu = 0; yu < l_runIndx.size(); yu++) if (l_runIndx.at(yu) == l_runIndx_old.at(we)) { fnd = true; break; }
               if (!fnd) missing_indxs.push_back(we);
              }
          for (unsigned int rt = 0; rt < missing_indxs.size(); rt++)
              {
               prv = 0;
               for (unsigned int a= 0; a < missing_indxs.at(rt); a++) prv += l_nfiles.at(a);
               names.erase(names.begin() + prv, names.begin() + l_nfiles.at(missing_indxs.at(rt)));
               l_totalevnts -= l_TotHits.at(missing_indxs.at(rt));
               l_nfiles.erase(l_nfiles.begin() + missing_indxs.at(rt));
               l_Runs.erase(l_Runs.begin() + missing_indxs.at(rt));
               l_runIndx.erase(l_runIndx.begin() + missing_indxs.at(rt));
               if (l_ext.size() != *std::max_element(l_nfiles.begin(), l_nfiles.end()))
                   l_ext.resize(*std::max_element(l_nfiles.begin(), l_nfiles.end()));
               l_TotHits.erase(l_TotHits.begin() + missing_indxs.at(rt));
               l_srate.erase(l_srate.begin() + missing_indxs.at(rt));
               l_nChans.erase(l_nChans.begin() + missing_indxs.at(rt));
               l_chName.erase(l_chName.begin() + missing_indxs.at(rt));
               nRuns--;
              }
         }
     }

  std::vector<double> l_triggtimeVec;
  std::vector<unsigned int> trckevntlast;
  std::vector<std::vector<double>> l_voltage;
  std::vector<long double> l_phystime;
  std::vector<long double> l_ordrtime;
  std::vector<unsigned int> l_channel;  // temporary vector to hold channel number for current event
  std::vector<long int> l_nEvntN;
  unsigned int lastEvnt;                // Last event processed from this run
  unsigned int firstEvnt;               // First event processed from this run
  unsigned int totEvnt_old;             // Last event from the l_Hits vector before it was modified
  unsigned long int ievent = 0;
  unsigned long int tevent = 0;        // Events writen out that actually have concidences (tevent = l_totalevnts/m_nchan for 100% efficiency)
  unsigned long int trackevnt = 0;    // Event number for tracking, to be reset at each run start
  bool skip = false;  // Used when a files is skipped to trigger the recalculation at the next available file
  prv = 0;            // Do not reuse, keeping track of processed files
  for (unsigned int d = 0; d < nRuns; d++)
      {
       trackevnt = 0;
       if (d > 0) prv += l_nfiles.at(d - 1);
       fnd = false; // flag for the first file of the run (not necesserely the first in the directory but the first accessible)
       l_triggtimeVec.clear();
       lastEvnt = 0;
       firstEvnt = 0;
       totEvnt_old = l_TotHits.at(d);
       for (unsigned int ifile = 0; ifile < l_nfiles.at(d); ifile++)
           {
            // Create filename ot open
            memset(&ifname[0], 0, sizeof(ifname));
            indx = prv + ifile;
retry:
            if (strcmp(ext, "") == 0) sprintf(ifname, "%s%s", dir, (names.at(indx)).Data());
            else {
                  if (ifile == 0) sprintf(ifname, "%s%s.%s", dir, (names.at(indx)).Data(), ext);
                  else sprintf(ifname, "%s%s.%s_%04d", dir, (names.at(indx)).Data(), ext, ifile);
                 }
            file.open(ifname);
            // File Open checks and skip conditions
            if (!file)
               {
                std::cout << __FUNCTION__ << " ERROR: Unable to find input file: " << ifname << ", skipping..." << std::endl;
                names.erase(names.begin() + indx);
                l_nfiles.at(d) = l_nfiles.at(d) - 1;
                if (l_ext.size() != *std::max_element(l_nfiles.begin(), l_nfiles.end()))
                    l_ext.resize(*std::max_element(l_nfiles.begin(), l_nfiles.end()));
                ifile--;
                skip = true; // variable to correct the number of events later on
                if (indx < (prv + l_nfiles.at(d) - 1)) goto retry;
                else {
                      if (l_nfiles.at(d) == 0)
                         {
                          l_nfiles.erase(l_nfiles.begin() + d);
                          if (l_ext.size() != *std::max_element(l_nfiles.begin(), l_nfiles.end())) l_ext.resize(*std::max_element(l_nfiles.begin(), l_nfiles.end()));
                          l_totalevnts -= l_TotHits.at(d);
                          l_Runs.erase(l_Runs.begin() + d);
                          l_runIndx.erase(l_runIndx.begin() + d);
                          l_TotHits.erase(l_TotHits.begin() + d);
                          l_srate.erase(l_srate.begin() + d);
                          l_nChans.erase(l_nChans.begin() + d);
                          l_chName.erase(l_chName.begin() + d);
                          nRuns--; d--;
                         }
                      else {
                            l_TotHits.at(d) -= (totEvnt_old - lastEvnt);
                            l_totalevnts -= (totEvnt_old - lastEvnt);
                           }
                      break;
                     }
               }
            std::cout << "--> Openning file: " << ifname << std::endl;
            indx = 0; // Reuse indx for  number of points
            for (unsigned int i = 0; i < 7; i++) {line.clear(); getline(file, line);} // skip the first 7 lines, may include information later
            // Read out time and voltage of all events in file
            l_voltage.clear();
            l_phystime.clear();
            l_ordrtime.clear();
            l_channel.clear();
            l_nEvntN.clear();
            while (file.good())
                  {
                   phystime = 0.0;
                   ordrtime = 0.0;
                   channel = -1;
                   line.clear();
                   getline(file, line);
                   sscanf(line.c_str(), "%*s %lu", &ievent);
                   // operations at the begining of the file
                   if (l_voltage.size() == 0 && skip)
                      {
                       l_TotHits.at(d) -= (ievent - lastEvnt);
                       l_totalevnts -= (ievent - lastEvnt);
                       CheckStartStopEvnt(&evt1, &evt2, (long int)firstEvnt, l_TotHits.at(d));
                       skip = false;
                      }
                   // Find out if start - stop events correspond to expectations (first event of first file) and fix it
                   if (!fnd) 
                      {
                       fnd = CheckStartStopEvnt(&evt1, &evt2, (long int)ievent, l_TotHits.at(d));
                       firstEvnt = ievent;
                      }
                   if ((long int)ievent < evt1) continue;
                   if ((long int)ievent >= evt2) break;
                   LGADBase::ProgressBar(floor((float)(ievent-evt1)/(float)(2*m_nchan)+tevent*m_nchan), ceil((float)l_totalevnts/(float)m_nchan));
                   line.clear();
                   getline(file, line);
                   // if (l_voltage.size() == 0) sscanf(line.c_str(), "%*s %u %*s %lf %*s %*u %*s %*lf %*s %lf %*s %u", &channel, &phystime, &ordrtime, &indx);
                   // else sscanf(line.c_str(), "%*s %u %*s %lf %*s %*u %*s %*lf %*s %lf %*s %*u", &channel, &phystime, &ordrtime);
                   if (l_voltage.size() == 0) sscanf(line.c_str(), "%*s %u %*s %lf %*s %*d %*s %*lf %*s %lf %*s %*lf %*s %*lf %*s %*lf %*s %u", &channel, &phystime, &ordrtime, &indx);
                   else sscanf(line.c_str(), "%*s %u %*s %lf %*s %*d %*s %*lf %*s %lf %*s %*lf %*s %*lf %*s %*lf %*s %*u", &channel, &phystime, &ordrtime);
                   if (channel == -1) continue;
                   l_voltage.push_back(std::vector<double>());
                   (l_voltage.back()).reserve(indx);
                   line.clear();
                   getline(file, line);
                   line = reduce(line, "", "DataSamples");
                   line = trim(line, " ");
                   iss.str(std::string());
                   iss.clear();
                   iss.str(line);
                   for (unsigned int i = 0; i < indx; i++)
                       {
                        iss >> tri; 
                        (l_voltage.back()).push_back(tri);
                       }
                   l_phystime.push_back(phystime);
                   l_ordrtime.push_back(ordrtime);
                   l_channel.push_back(channel+1);
                   l_nEvntN.push_back(ievent);
                  }
            lastEvnt = ievent;
            file.close();
            trckevntlast.clear();
            if (LGADBase::GetTrackComb()) trckevntlast.resize(1 + m_nplanes.at(d) + m_nDUT.at(d) + m_nRef.at(d), 0);
            // Trying to find conicidences within the file we proces
            ievent = 0;
            std::vector<unsigned int> ChNo;
            fnd = false;
            unsigned int temp = 0;
            for (unsigned int i = 0; i < l_voltage.size(); i++)
                {
                 temp = i;
                 ChNo.clear();
                 ChNo.push_back(l_channel.at(i));
                 if ((l_chName.at(d)).size() > 1)
                    {
                     // Loop though all events in the file, find if there is any in coincidence with the one we have and populate vectors
                     for (unsigned int j = i + 1; j < l_voltage.size(); j++)
                         {
                          if (ChNo.at(0) == l_channel.at(j)) continue; // Ignore events in the same chanel as the one under test
                          else if (fabs(l_phystime.at(i) - l_phystime.at(j)) <= indx*(1e3/l_srate.at(d))) // cut for coicidence within acuisition time
                                  {  
                                   int wch = -1;
                                   for (unsigned int h = 0; h < m_channels.size(); h++) if (m_channels.at(h) == l_channel.at(j)) wch = h;
                                   if (wch == -1)
                                      {
                                       std::cout << __FUNCTION__ << " ERROR: Cannot determine which channel this event corrsponds to" << std::endl;
                                       continue;
                                      }
                                   m_physt.at(wch) = l_phystime.at(j);
                                   m_ordrt.at(wch) = l_ordrtime.at(j);
                                   m_w.at(wch) = l_voltage.at(j);
                                   m_t.at(wch).clear();
                                   m_t.at(wch).reserve(indx);
                                   // This is channel jitter realignement, it assumes time resolutions < sampling interval
                                   int bin = nearbyint((l_ordrtime.at(i) - l_ordrtime.at(j))*l_srate.at(d)*0.001);
                                   for (unsigned int k = 0; k < indx; k++) m_t.at(wch).push_back((1 / (double)l_srate.at(d))*((int)k - bin)*(1e-6));
                                   ChNo.push_back(l_channel.at(j));
                                   i = j; // Aassumes time-ordered aquisition
                                  }
                          else continue;
                         }
                     }
                 // We found events in coincidednce (at least one)
                 if (((temp-i) > 0 && (l_chName.at(d)).size() > 1) || ((temp-i) == 0 && (l_chName.at(d)).size() == 1))
                    {
                     // If you have multiple active channels, check that you have found events in conincidence 
                     // in all channels, if not fill the ones that are not found with -1.
                     for (unsigned int f = 0; f < m_nchan; f++) 
                         {
                          fnd = false;
                          for (unsigned int g = 0; g < (temp-i)+1; g++) if (m_channels.at(f) == ChNo.at(g)) fnd = true;
                          if (!fnd)
                             {
                              m_t.at(f).clear();
                              m_physt.at(f) = -1;
                              m_ordrt.at(f) = -1;
                              m_w.at(f).clear();
                              m_t.at(f).reserve(indx);
                              m_w.at(f).reserve(indx);
                              for (unsigned int k = 0; k < indx; k++)
                                  {
                                   m_w.at(f).push_back(-1);
                                   m_t.at(f).push_back((1 / (double)l_srate.at(d))*k*(1e-6));
                                  }
                             }     
                         }
                     // Save the event we were checking for conicidences with at the previous setp
                     int wch1 = -1;
                     for (unsigned int w = 0; w < m_channels.size(); w++) if (m_channels.at(w) == ChNo.at(0)) wch1 = w;
                     if (wch1 == -1)
                        {
                         std::cout << __FUNCTION__ << " ERROR: Cannot determine which channel this event corrsponds to" << std::endl;
                         continue;
                        }
                     m_w.at(wch1) = l_voltage.at(temp);
                     (m_t.at(wch1)).clear();
                     (m_t.at(wch1)).reserve(indx);
                     for (unsigned int k = 0; k < indx; k++) (m_t.at(wch1)).push_back((1/(double)l_srate.at(d))*k*(1e-6));
                     m_physt.at(wch1) = l_phystime.at(temp);
                     m_ordrt.at(wch1) = l_ordrtime.at(temp);
                     m_npoints.at(0) = indx;
                     m_srate.at(0) = ceil(l_srate.at(d)*1e6);
                     m_event = trackevnt + temp;
                     l_triggtimeVec.push_back(l_phystime.at(temp)/1e9);
                     if (LGADBase::GetTrackComb())
                        {
                         if (LGADBase::GetTrackPackage() == "Cory") LGADBase::CombineTrackCory(m_event, trckevntlast, d);
                         else if (LGADBase::GetTrackPackage() == "EUTelescope") LGADBase::CombineTrackEUTel(m_event, trckevntlast, d);
                         else if (LGADBase::GetTrackPackage() == "TBA") LGADBase::CombineTrackTBA(m_event, trckevntlast, d);
                         else {
                               std::cout << __FUNCTION__ << "          ERROR: Unsopported tracking package " << LGADBase::GetTrackPackage() << std::endl;
                               return false;
                              }
                        }
                     m_trees.at(0)->Fill();
                     ievent++; // events in this file that actually have conicdences
                    }
                 LGADBase::ProgressBar(floor((float)(l_voltage.size()+temp)/(float)(2*m_nchan)+tevent), ceil((float)l_totalevnts/(float)m_nchan));
                }
            trackevnt += temp; // event umber used for trakcing, reset per run
            tevent += ievent; // events actually writen out which have conicidences, not rest per run
           }
       LGADBase::CalcTrigFr(l_triggtimeVec, m_trigDt, m_trigFr, ievent);
      }

  std::cout << "--> Wrote " << tevent << " events into output file " << m_ofile->GetName() << ", Efficiency = " << (float)(tevent*m_nchan) / (float)l_totalevnts << std::endl;
  LGADBase::CloseTrackFiles();
  delete b;
  m_event = tevent;

  TListIter l_next_object(gDirectory->GetList());
  l_next_object.Reset();
  TObject* l_obj;
  while ((l_obj=l_next_object()))
        {
         l_obj->Write("", 6);
         delete l_obj;
        }
  m_ofile->Close();

  return true;
}
// --------------------------------------------------------------------------------------------------------------
// Readout function for the Labview TiCAs program, Lecroy scope and old tektronix style scope on a txt file saving mode 
// A seperate txt file created per channel and per event.
bool LGADBase::WriteLabTXT(const char* dir, std::vector<TString> names, const char* ext, unsigned int evt1, unsigned int evt2)
{
    // File name format string for scanf with placeholders for the dir, channel, the event number and extension
    std::istringstream iss(std::ios_base::app | std::ios_base::in | std::ios_base::out);
    double time;
    double voltage;
    std::string line;
    std::ifstream file;
    std::vector <unsigned int> l_npoints;
    std::vector<double> l_triggtimeVec;
    l_npoints.clear();
    bool l_labview = false;
    bool end = false;

    // Test input file number
    unsigned int nfiles = names.size();
    if (nfiles == 0)
       {
        std::cout << __FUNCTION__ << " ERROR: Unable to find any input files in directory: " << dir << std::endl;
        return false;
       }

    // Format the filename vector with respect to extention
    std::string fname; // temporary
    if (strcmp(ext, "") == 0)
       {
        std::size_t found;
        for (unsigned int k = 0; k < nfiles; k++)
            {
             if (!(names.at(k).Contains(".")))
                {
                 names.erase(names.begin() + k);
                 nfiles--;
                 k--;
                 continue;
                }
             else {
                   found = (names.at(k)).Last('.');
                   if (k == 0) ext = (names.at(k)(found + 1, names.at(k).Length())).Data();
                   else { 
                         fname = names.at(k)(found + 1, names.at(k).Length());
                         if (fname != ext)
                            {
                             std::cout << __FUNCTION__ << " ERROR: Different file extentions found: " << fname.c_str()
                                       << " and " << ext << "! Abording..." << std::endl;
                             return false;
                            }
                        } 
                   names.at(k) = names.at(k)(0, found);
                  }
            }
        }

    // Scan event numbers, channels and do necessery tests
    char ifname[2048];          // temporary
    unsigned int ch = 0;        // temporary
    bool got = false;           // temporary
    unsigned int evtfd = 0;     // temporary
    unsigned int max_ent = 0;
    unsigned int min_ent = 0;
    fname.clear();
    for (unsigned r = 0; r < nfiles; r++)
        {
         memset(&ifname[0], 0, sizeof(ifname));
         evtfd = 0; ch = 0;
         sscanf(names.at(r), "C%u%s", &ch, ifname);
         evtfd = stoi(((std::string)ifname).substr(((std::string)ifname).length() - 5, 5));
         got = false;
         for (unsigned int a = 0; a < m_channels.size(); a++) 
             {
              if (m_channels.at(a) == ch) { got = true; break; }
             }
         if (!got) { m_channels.push_back(ch); m_nchan++; }
         if (r == 0) 
            { 
             max_ent = evtfd;  
             min_ent = evtfd;
             fname = ((std::string)ifname).substr(0, ((std::string)ifname).length() - 5);
            }
         else {
               if (evtfd > max_ent) max_ent = evtfd;
               if (evtfd < min_ent) min_ent = evtfd;
               line = ((std::string)ifname).substr(0, ((std::string)ifname).length() - 5);
               if (!fname.compare(line))
                  {
                   std::cout << __FUNCTION__ << " ERROR: Several file names found: " << line.c_str()
                             << " and " << fname << "! Abording..." << std::endl;
                   return false;
                  }
              }
        }

    // Check available events for processing
    CheckStartStopEvnt(&evt1, &evt2, min_ent, max_ent);

    // Using frist event number find the no. of points per cahannel and oscilloscope type
    char oscillos[64];
    unsigned int npoints = 0;     // temporary
    int count = 0;                // temporary
    for (unsigned int i = 1; i < m_nchan; i++)
        { 
         memset(&ifname[0], 0, sizeof(ifname));
         sprintf(ifname, "%sC%d%s%05d.%s", dir, m_channels.at(i), fname.c_str(), evt1, ext);  // build the filename for the first event in the first channel
         file.open(ifname);
         if (file)
            {
             line.clear();
             memset(&oscillos[0], 0, sizeof(oscillos));
             getline(file, line);
             std::replace(line.begin(), line.end(), ',', ' ');
             count = sscanf(line.c_str(), "%*s %*lu %s", oscillos);
             if (count && std::string(oscillos) == "Waveform")
                {
                 line.clear();
                 getline(file, line);
                 std::replace(line.begin(), line.end(), ',', ' ');  // replace all commas (if any) to spaces
                 sscanf(line.c_str(), "%*s %*lu %*s %u ", &npoints);
                 l_labview = false;
                }
             else {
                   l_labview = true;
                   memset(&oscillos[0], 0, sizeof(oscillos));
                   if (sscanf(line.c_str(), "%*u %*u %u %*u %*f %*f %*f %*f %*f %*f %*u %*f %*f %*f %*f \"%*2d %*3s %*4d\" \"%*2d:%*2d:%*2d:%*2d\" \"%*s\" %*u %*u %*u %*u %*lu %*lu ", &npoints) > 0)
                      {
                       strcpy(oscillos, "Agilent");
                      }
                   else if (sscanf(line.c_str(), "%*u;%*u;%*3s;%*2s;%*3s;\"%*3s  %*2s %*8s  %*f%*6s  %*f%*6s  %*u %*6s  %*4s %*3s %*4s\";%u;%*s;\"%*s\";%*f;%*f;%*u;\"%*s\";%*f;%*f;%*f;%*u", &npoints) > 0)
                           {
                            strcpy(oscillos, "Tektronix");
                           }
                   else if (sscanf(line.c_str(), "%*100s %*100s Vertical Offset %*f; Vertical Gain %*f; Horizontal Offset %*f; Horizontal Interval %*f; %u points  from %*d to %*d; Trigger time %*2d:%*2d:%*2d.%*4d  %*2d/%*2d/%*4d", &npoints) > 0)
                           {
                            strcpy(oscillos, "LabViewCom");
                           }
                   else {
                         std::cout << __FUNCTION__ << " ERROR: Failed recognizing waveform file header: " << line.c_str() << std::endl;
                         return false;
                        }
                  }
             l_npoints.push_back(npoints);
             file.close();
            }
         else {
               std::cout << __FUNCTION__ << " ERROR: Could not open file: " << ifname << ", first event file for channel " << m_channels.at(i) << ". Channel will be ignored!" << std::endl;
               m_channels.erase(m_channels.begin() + i);
               nfiles = nfiles - ((float)nfiles / (float)m_nchan);
               m_nchan--;             
               i--;
              }
        }

    // Time of the standard time_t type: the number of seconds since Jan 1, 1970 
    double trigtime_first = 0.0; // time of the first event
    m_trigtime = 0.0;

    // create the ntuple
    if (!LGADBase::CreateOutputFile(m_ofdir, m_ofname, m_channels))
       {
        std::cout << __FUNCTION__ << " ERROR: Output Ntuple creation failed: " << m_ofname << std::endl;
        return false;
       }
    if (!l_labview) LGADBase::SetScale(m_channels, m_nchan, &m_scale);
    npoints = *min_element(l_npoints.begin(), l_npoints.end());

    bool comma = false;
    unsigned int oi = 0;
    unsigned int zr = 0;
    float oi1 = 0;
    float a = 1.;
    // Loop over all event files
    unsigned int n = 0; // Writen events counter
    char VUnit[6]; // Voltage Units placeholder
    for (m_event = evt1; true; m_event++)  // NB: an infinite loop -- for (int m_event=0; m_event<10; ++m_event)
        {
         if (((int)m_event >= evt2) || (m_event*m_nchan == nfiles)) break;
         LGADBase::ProgressBar(m_event - evt1, evt2 - evt1);
         if (end) break;
         bool first = false;
         bool bad = false; // Set per event, if one channel is bad, ignore the event
         for (unsigned int ich = 0; ich < m_nchan; ich++) // loop over all cahnnels
             {
              if (!bad) // event is bad, we ignore all additional cahnnels
                 {
                  m_w[ich].clear();
                  m_t[ich].clear();
                  m_w[ich].resize(npoints);
                  m_t[ich].resize(npoints);
                  m_npoints.at(ich) = l_npoints.at(ich);
                  memset(&ifname[0], 0, sizeof(ifname));
                  sprintf(ifname, "%sC%d%s%05d.%s", dir, m_channels.at(ich), fname.c_str(), m_event, ext); // build the filename
                  file.open(ifname);
                  if (!file)
                     {
                      if (m_trees.at(0)->GetEntries() == 0)
                         {
                          std::cout << __FUNCTION__ << " ERROR: Could not open any input file" << ifname << std::endl;
                          return false;
                         }
                      else if ((m_event*m_nchan) < nfiles) continue; // If a file is missing just go to the next
                      else {end = true; break;}
                     }
                  if (!l_labview) for (unsigned int h = 0; h < 4; h++) { line.clear(); getline(file, line); } // read the line #1,07-Mar-2007 07:39:21,0
                  else {
                        line.clear();
                        getline(file, line);
                        std::replace(line.begin(), line.end(), ',', ' ');
                        if (std::string(oscillos) == "Agilent" || std::string(oscillos) == "LabViewCom" ) 
                           {
                            // get the scale for each channel event
                            if (std::string(oscillos) == "Agilent")
                               {
                                sscanf(line.c_str(), "%*u %*u %*u %*u %*f %*f %*f %f %*f %*f %*u %*f %*f %*f %*f \"%*2d %*3s %*4d\" \"%*2d:%*2d:%*2d:%*2d\" \"%*s\" %*u %*u %*u %*u %*lu %*lu ", &m_scale.at(ich));
                                m_scale.at(ich) = m_scale.at(ich) / 0.000131732;
                               }
                            else sscanf(line.c_str(), "%*100s %*100s Vertical Offset %*f; Vertical Gain %f; Horizontal Offset %*f; Horizontal Interval %*f; %*u points  from %*d to %*d; Trigger time %*2d:%*2d:%*2d.%*4d  %*2d/%*2d/%*4d", &m_scale.at(ich));
                           }
                        else if ( std::string(oscillos) == "Tektronix") 
                                {
                                 memset(&VUnit[0], 0, sizeof(VUnit));
                                 sscanf(line.c_str(), "%*u;%*u;%*3s;%*2s;%*3s;\"%*3s  %*2s %*8s  %f%6s  %*f%*6s  %*u %*6s  %*4s %*3s %*4s\";%*u;%*s;\"%*s\";%*f;%*f;%*u;\"%*s\";%*f;%*f;%*f;%*u", &m_scale.at(ich), VUnit);
                                 if (std::string(VUnit) == "mV/div") m_scale.at(ich) = m_scale.at(ich)*0.001;
                                 line.clear();
                                 getline(file, line);
                                }
                       }
                  if (!first)  // use the time stamp from the first channel only
                     {
                      char month_str[3];
                      struct tm timeinfo;
                      int year = 0;
                      if (!l_labview)
                         {
                          std::replace(line.begin(), line.end(), ',', ' ');
                          sscanf(line.c_str(), "#%*d %2d-%3s-%4d %2d:%2d:%2d %*d ", &timeinfo.tm_mday, month_str, &year, &timeinfo.tm_hour, &timeinfo.tm_min, &timeinfo.tm_sec);
                         }
                      else {
                            if (std::string(oscillos) == "Agilent")
                               {
                                sscanf(line.c_str(), "%*u %*u %*u %*u %*lf %*lf %*lf %*lf %*lf %*lf %*u %*lf %*lf %*lf %*lf \"%2d %3s %4d\" \"%2d:%2d:%2d:%2u\" \"%*20s\" %*u %*u %*u %*u %*lf %*lu ", 
                                                    &timeinfo.tm_mday, month_str, &year, &timeinfo.tm_hour, &timeinfo.tm_min, &timeinfo.tm_sec, &oi);
                               }
                            else if (std::string(oscillos) == "Tektronix")
                                    {
                                     std::replace(line.begin(), line.end(), ',', ' ');
                                     sscanf(line.c_str(), "%2d/%2d/%4d;%2d:%2d:%2d %3s", &timeinfo.tm_mon, &timeinfo.tm_mday, &year, &timeinfo.tm_hour, &timeinfo.tm_min, &timeinfo.tm_sec, month_str);
                                     if (std::string(month_str) == "PM") timeinfo.tm_hour = timeinfo.tm_hour + 12;
                                    }
                            else if (std::string(oscillos) == "LabViewCom")
                                    {
                                     sscanf(line.c_str(), "%*100s %*100s Vertical Offset %*f; Vertical Gain %*f; Horizontal Offset %*f; Horizontal Interval %*f; %*u points  from %*d to %*d; Trigger time %2d:%2d:%2d.%4u  %2d/%2d/%4d",
                                            &timeinfo.tm_hour, &timeinfo.tm_min, &timeinfo.tm_sec, &oi, &timeinfo.tm_mday, &timeinfo.tm_mon, &year);
                                    }
                            else {
                                  std::cout << __FUNCTION__ << " ERROR: Failed recognizing oscilloscope type: " << oscillos << std::endl;
                                  return false;
                                 }
                           }
                      timeinfo.tm_year = year - 1900;
                      if (!l_labview) timeinfo.tm_mon = std::string("Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec").find(month_str) / 4;   // convert the three-letter month name into a number 0..11
                      else if (std::string(oscillos) == "Agilent") timeinfo.tm_mon = std::string("JAN FEB MAR APR MAY JUN JUL AUG SEP OCT NOV DEC").find(month_str) / 4;
                      m_trigtime = mktime(&timeinfo);
                      if (l_labview && std::string(oscillos) != "Tektronix") 
                         {
                          if (std::string(oscillos) == "LabViewCom") a = 10000;
                          else  a = 100;
                         }
                      if ((int)m_event == evt1) // time of the first event
                         {
                          trigtime_first = m_trigtime;
                          if (l_labview && std::string(oscillos) != "Tektronix") oi1 = (float)oi/a;
                         }
                      m_trigtime = (long long int)m_trigtime - (long long int)trigtime_first;
                      if (m_trigtime < 0) { while(m_trigtime < 0) m_trigtime = m_trigtime + 3600; }
                      else if (m_trigtime >= 3600 && n < 12) {trigtime_first = trigtime_first + 3600; m_trigtime = m_trigtime-3600;}
                      if (l_labview && std::string(oscillos) != "Tektronix") m_trigtime = m_trigtime + ((float)oi/a) - oi1;
                      l_triggtimeVec.push_back(m_trigtime);
                      if (n != 0) first = true;
                     }
                  line.clear();
                  if ((l_labview && std::string(oscillos) != "Tektronix") || !l_labview) getline(file, line);    // skip the line Time,Ampl
                  // Read the data
                  Int_t nlines = 0;
                  for (unsigned int i = 0; i < m_npoints.at(ich); i++)
                      {
                       line.clear();
                       if (getline(file,line)) // get line from the file
                          {
                           if (!first && n == 0) 
                              {
                               first = true;
                               if ((int)m_event == evt1) { if (std::count(line.begin(), line.end(), ',') > 1) comma = true; }
                              }
                           if (comma)
                              {
                               if (std::count(line.begin(), line.end(), ',') == 3)
                                  {
                                   line.replace(line.find(","), 1, ".");
                                   size_t found = line.find(",");
                                   line.replace(line.find(",", found+1), 1, ".");
                                  }
                               else if (std::count(line.begin(), line.end(), ',') == 2)
                                       {
                                        size_t found = line.find(",");
                                        if (found <= line.size()/2) line.replace(line.find(",", found + 1), 1, ".");
                                        else line.replace(line.find(","), 1, ".");
                                       }
                              }
                           std::replace(line.begin(), line.end(), ',', ' ');
                           iss.str(std::string());
                           iss.clear();
                           iss.str(line);
                           iss >> time >> voltage; // read the space separated data from the line
                           m_w[ich].at(i) = voltage;
                           m_t[ich].at(i) = time; // convert the time to ns
                           ++nlines;
                          }
                       else {
                             std::cout << __FUNCTION__ << " ERROR: file " << ifname << " has " << nlines << " points instead of " << m_npoints.at(ich) << " npoints." << std::endl;
                             bad = true;
                             l_triggtimeVec.pop_back();
                             break;
                            }
                      }
                  if (!bad) // if event is not bad we take it into acount
                     {
                      m_srate[ich] = (Long64_t)ceil((double)(m_npoints[ich] - 1)/abs(m_t[ich].front() - m_t[ich].back()));
                      n++;
                     }
                  file.close();
                 }
             }
         if (n > 0 && !bad) m_trees.at(0)->Fill();
        }

    // Fix the numbers of events
    m_event = n / m_nchan;
    std::cout << "Wrote " << m_event << " events into output file " << m_ofile->GetName() << std::endl;
    LGADBase::CalcTrigFr(l_triggtimeVec, m_trigDt, m_trigFr, m_trees.at(0)->GetEntries());

    TListIter l_next_object(gDirectory->GetList());
    l_next_object.Reset();
    TObject* l_obj;
    while ((l_obj=l_next_object()))
          {
           l_obj->Write("", 6);
           delete l_obj;
          }
    m_ofile->Close();

    return true;
}
// --------------------------------------------------------------------------------------------------------------
// Tektronix ASCI file format. The foramt is not structured as one event per channel and per aquistion but as a single ASCI file
// containing all events for a specified channel in a continues stream. The functions expects to find as
// many files at the run directorry as channels and all with the same name with the suffix _Ch. THis is a cross vetween the Sampic 
// and the  LabTXT data functions, since one can have several runs in the same directory which would inevitably have different names.
bool LGADBase::WriteTectronixTXT(const char* dir, std::vector<TString> names, const char* ext, unsigned int evt1, unsigned int evt2)
{
  // Temporary variables to be uses internaly
  char ifname[2048];
  std::string fname;
  long int lCurPos;
  static FILE* fd_file = NULL;
  unsigned int ch = 0;
  unsigned int npts = 0;
  double volt = 0;
  double time = 0;
  unsigned int frame = 0;
  double interval = 0;
  bool fnd = false;

  // Determine the number of runs contained in the list of files 
  // to be analysed and the number of channels per run
  std::vector<TString> l_Runs;          // Base name of the run file
  std::vector<TString> l_ext;           // File extention
  std::vector<unsigned int> l_nfiles;   // Files per run
  std::vector<unsigned int> l_runIndx;  // Run indexes for reconstruciton synchronization
  std::vector<std::vector<unsigned int>> l_nchannels;
  unsigned int nRuns = 0;
  nRuns = RunCount(names, l_nfiles, l_Runs, l_ext, l_runIndx, ext, &l_nchannels);

  // Checking that we only have run files in the input folder
  if (l_ext.size() > 1)
     {
      std::cout << __FUNCTION__ << " ERROR: Multiple file extentions found in input directory: " << dir << std::endl;
      std::cout << __FUNCTION__ << "        Abording......." << std::endl;
      return false;
     }

  // Checking that each run only has as many files as channels
  for (unsigned int hg = 0; hg < l_nfiles.size(); hg++)
      {
       if (l_nfiles.at(hg) != (l_nchannels.at(hg)).size())
          {
           std::cout << __FUNCTION__ << " ERROR: Run with either incomplete or additional files found: " << l_Runs.at(hg) << ", excluding..." << std::endl;
           l_nfiles.erase(l_nfiles.begin() + hg);
           l_Runs.erase(l_Runs.begin() + hg);
           l_runIndx.erase(l_runIndx.begin() + hg);
           l_nchannels.erase(l_nchannels.begin() + hg);
           nRuns--; hg--;
          }
      }

  // Open the files and test the number of points and sizes to veryfy everything
  std::vector<std::vector<int>>      l_npoints;
  std::vector<unsigned int>          l_nRunEvnts;
  std::vector<std::vector<Long64_t>> l_snrate;
  for (unsigned int r = 0; r < nRuns; r++)
      {
       for (unsigned int i = 0; i < (l_nchannels.at(r)).size(); i++)
           {
            frame = 0;
            memset(&ifname[0], 0, sizeof(ifname));
            sprintf(ifname, "%s%s_Ch%i.%s", dir, (l_Runs.at(r)).Data(), (l_nchannels.at(r)).at(i), (l_ext.at(0)).Data());
            fd_file = fopen(ifname, "r");
            if (fd_file == NULL)
               {
fileprob1:
                std::cout << __FUNCTION__ << " ERROR: Failed to open txt file: " << ifname << ", skipping channel..." << std::endl;
fileprob2:
                l_nfiles.at(r)--; i--;
                (l_nchannels.at(r)).erase((l_nchannels.at(r)).begin() + i);                
                if (l_nfiles.at(r) == 0)
                   {
                    std::cout << __FUNCTION__ << " ERROR: No channel files for this run processed, skipping run " << l_Runs.at(r) << "..." << std::endl;
                    l_nchannels.erase(l_nchannels.begin() + r);
                    l_nfiles.erase(l_nfiles.begin() + r);
                    l_Runs.erase(l_Runs.begin() + r);
                    l_runIndx.erase(l_runIndx.begin() + r);
                    nRuns--; r--;
                   }
                fclose(fd_file);
                continue;
               }
            // GetFile size
#ifdef _WIN32
            lCurPos = _ftelli64(fd_file);
            _fseeki64(fd_file, 0, 2);
            size_t size = _ftelli64(fd_file);
            _fseeki64(fd_file, lCurPos, 0);
#else
            lCurPos = ftell(fd_file);
            fseek(fd_file, 0, 2);
            size_t size = ftell(fd_file);
            fseek(fd_file, lCurPos, 0);
#endif
            if (size <= 0) goto fileprob1;
            // Resize the vectors for the properties of the channel
            if (i == 0)
               {
                l_npoints.push_back(std::vector<int>());
                l_snrate.push_back(std::vector<Long64_t>());
               }
            // Get the number of points and the rate for each channel by scannig the first 2 lines
            for (unsigned int ap = 0; ap < 2; ap++)
                {
                 memset(&ifname[0], 0, sizeof(ifname));
                 if (fgets(ifname, sizeof(ifname), fd_file) == NULL) goto fileprob1;
                 if (ap == 0) ch = sscanf(ifname, "%*s %*s %u %*s %lf %lf", &npts, &time, &volt);
                 else ch = sscanf(ifname, "%*s %*s %lf %*s %lf %lf", &interval, &time, &volt);
                 if (ch != 3 || npts == 0)
                    {
                     std::cout << __FUNCTION__ << " ERROR: Invalid parameters from ASCI file, skipping channel " << ifname << std::endl;
                     goto fileprob2;
                    }
                }
            // Ignore the next four lines
            for (unsigned int ap = 0; ap < 4; ap++)
                {
                 memset(&ifname[0], 0, sizeof(ifname));
                 if (fgets(ifname, sizeof(ifname), fd_file) == NULL) goto fileprob1;
                }    
            // Find the events for the run but only for the first channel, line 7
            if (i == 0)
               {
                // Get the number of events in the run
                memset(&ifname[0], 0, sizeof(ifname));
                if (fgets(ifname, sizeof(ifname), fd_file) == NULL) goto fileprob1;
                ch = sscanf(ifname, "%*s %*s %u %*s %lf %lf", &frame, &time, &volt);
                if (ch != 3 || frame <= 0)
                   {
                    std::cout << __FUNCTION__ << " ERROR: Invalid parameters from ASCI file, skipping channel " << ifname << std::endl;
                    goto fileprob2;
                   }
                l_nRunEvnts.push_back(frame);
               }
             (l_npoints.at(r)).push_back(npts);
             (l_snrate.at(r)).push_back((Long64_t)(ceil((double)1/interval)));
             fclose(fd_file);
            }
        }

  // Fix the number of channels
  if (!FixBaseQuant(l_nchannels, m_channels, m_nchan)) return false;

  // Create the ntuple
  if (!CreateOutputFile(m_ofdir, m_ofname, m_channels))
     {
      std::cout << __FUNCTION__ << " ERROR: Output Ntuple creation failed: " << m_ofname << std::endl;
      return false;
     }

  // Tracking Mode related actions:
  // - Check the oscilloscope delay vector and populate it correctly if needed
  // - Open the required data file and get the collection of trees
  // - Add the required parts to the Ntuple
  if (LGADBase::GetTrackComb())
     {
      std::vector<unsigned int> l_runIndx_old = l_runIndx;
      if (!OpenTrackFile(m_Trackdatadir, m_Trackdataname, l_runIndx))
         {
          std::cout << __FUNCTION__ << " ERROR: Failed oppening tracking files in directory: " << m_Trackdatadir << std::endl;
          return false;
         }
      else {
            if (!AddTrckBranches(m_trees.at(0), m_trkpck))
               {
                std::cout << __FUNCTION__ << " Warning: Failed adding the tracking part on the output file: " << m_ofname << std::endl;
                std::cout << "                                       Switching to non-trackimg mode..." << std::endl;
                LGADBase::SetTrackComb(false);
               }
           }
      if (l_runIndx.size() != l_runIndx_old.size())
         {
          std::vector<unsigned int> missing_indxs;
          for (unsigned int we = 0; we < l_runIndx_old.size(); we++)
              {
               fnd = false;
               for (unsigned int yu = 0; yu < l_runIndx.size(); yu++) if (l_runIndx.at(yu) == l_runIndx_old.at(we)) { fnd = true; break; }
               if (!fnd) missing_indxs.push_back(we);
              }
          for (unsigned int rt = 0; rt < missing_indxs.size(); rt++)
              {
               l_nRunEvnts.erase(l_nRunEvnts.begin() + rt);
               l_npoints.erase(l_npoints.begin() + rt);
               l_snrate.erase(l_snrate.begin() + rt);
               l_nfiles.erase(l_nfiles.begin() + rt);
               l_Runs.erase(l_Runs.begin() + rt);
               l_runIndx.erase(l_runIndx.begin() + rt);
               l_nchannels.erase(l_nchannels.begin() + rt);
               nRuns--;
              }
         }
     }

  unsigned int np = 0;  // index for number of saved points to the vector
  unsigned int trigpoint = 0;
  double trigtime = 0;
  std::vector<unsigned int> trckevntlast;
  std::vector<std::vector<std::vector<double>>> l_w;  // channel, event, actuall vector
  std::vector<std::vector<std::vector<double>>> l_t;  // channel, event, actuall vector
  std::vector<std::vector<double>> l_trigtime;        // channel, event
  std::vector<std::vector<unsigned int>> l_trigpoint; // channel, event
  std::vector<unsigned int> l_envtNo;
  npts = 0; // index for number or information red in scanf
  for (unsigned int r = 0; r < nRuns; r++)
      {
       LGADBase::SetScale(m_channels, m_nchan, &m_scale);
       l_w.clear();
       l_t.clear();
       l_trigtime.clear();
       l_trigpoint.clear();
       trckevntlast.clear();
       l_w.resize((l_nchannels.at(r)).size());
       l_t.resize((l_nchannels.at(r)).size());
       l_trigtime.resize((l_nchannels.at(r)).size());
       l_trigpoint.resize((l_nchannels.at(r)).size());
       // Reduce the number of events by one to account for the initial TLU isseud trigger at the start of each run
       if (LGADBase::GetTrackComb()) 
          {
           l_nRunEvnts.at(r) = l_nRunEvnts.at(r) - 1;
           trckevntlast.resize(1 + m_nplanes.at(r) + m_nDUT.at(r) + m_nRef.at(r), 0);
          }
       for (unsigned int i = 0; i < (l_nchannels.at(r)).size(); i++)
           { 
            memset(&ifname[0], 0, sizeof(ifname));
            sprintf(ifname, "%s%s_Ch%i.%s", dir, (l_Runs.at(r)).Data(), (l_nchannels.at(r)).at(i), (l_ext.at(0)).Data());
            fd_file = fopen(ifname, "rb");
            if (fd_file == NULL)
               {
                std::cout << __FUNCTION__ << " ERROR: Could not open file " << std::string(ifname) << ". Skipping channel..." << std::endl;
 fileprob3:
                l_nfiles.at(r)--; i--;
                (l_nchannels.at(r)).erase((l_nchannels.at(r)).begin() + i);
                l_w.erase(l_w.begin() + i);
                l_t.erase(l_t.begin() + i);
                l_trigtime.erase(l_trigtime.begin() + i);
                l_trigpoint.erase(l_trigpoint.begin() + i);
                if (l_nfiles.at(r) == 0)
                   {
                    std::cout << __FUNCTION__ << " ERROR: No channel files for this run processed, skipping run " << l_Runs.at(r) << "..." << std::endl;
                    l_nchannels.erase(l_nchannels.begin() + r);
                    l_nfiles.erase(l_nfiles.begin() + r);
                    l_Runs.erase(l_Runs.begin() + r);
                    l_runIndx.erase(l_runIndx.begin() + r);
                    nRuns--; r--;
                   }
                fclose(fd_file);
                continue;
               }
            else if (i == 0) std::cout << "    Opening run: " << l_Runs.at(r) << std::endl;
            fname = ifname;
            (l_w.at(i)).clear();
            (l_t.at(i)).clear();
            (l_trigtime.at(i)).clear();
            (l_trigpoint.at(i)).clear();
            (l_w.at(i)).reserve(l_nRunEvnts.at(r));
            (l_t.at(i)).reserve(l_nRunEvnts.at(r));
            (l_trigtime.at(i)).reserve(l_nRunEvnts.at(r));
            (l_trigpoint.at(i)).reserve(l_nRunEvnts.at(r));
            if (i == 0) l_envtNo.reserve(l_nRunEvnts.at(r));
            // Bring the number of events to its original value
            if (LGADBase::GetTrackComb()) frame = l_nRunEvnts.at(r) + 1;
            else frame = l_nRunEvnts.at(r);
            for (unsigned int icount = 0; icount < frame; icount++)
                {
                 // Fix the start - stop event with respect to events available in the run
                 if (LGADBase::GetTrackComb())
                    {
                     if (icount == 0 && i == 0) CheckStartStopEvnt(&evt1, &evt2, 0, l_nRunEvnts.at(r)-1);
                     if ((long int)icount < evt1+1)
                        {
                         for (unsigned int wes = 0; wes < (l_npoints.at(r)).at(i); wes++) fgets(ifname, sizeof(ifname), fd_file);
                         icount++;
                         continue;
                        }
                     LGADBase::ProgressBar(icount-1+((evt2-evt1)*r)*(l_nchannels.at(r).size())+(evt2-evt1)*i, ((evt2-evt1)*nRuns)*(l_nchannels.at(r).size()));
                     if ((long int)icount > evt2+1) break;
                    }
                 else {
                       if (icount == 0 && i == 0) CheckStartStopEvnt(&evt1, &evt2, 0, l_nRunEvnts.at(r));
                       if ((long int)icount < evt1)
                          {
                           for (unsigned int wes = 0; wes < (l_npoints.at(r)).at(i); wes++) fgets(ifname, sizeof(ifname), fd_file);
                           icount++;
                           continue;
                          }
                       LGADBase::ProgressBar(icount+((evt2-evt1)*r)*(l_nchannels.at(r).size())+(evt2-evt1)*i, ((evt2-evt1)*nRuns)*(l_nchannels.at(r).size()));
                       if ((long int)icount > evt2) break;
                      }
                 np = 0;
                 (l_w.at(i)).push_back(std::vector<double>());
                 ((l_w.at(i)).back()).reserve((l_npoints.at(r)).at(i));
                 (l_t.at(i)).push_back(std::vector<double>());
                 ((l_t.at(i)).back()).reserve((l_npoints.at(r)).at(i));
                 (l_trigtime.at(i)).push_back(0.0);
                 (l_trigpoint.at(i)).push_back(0);
                 fnd = true;
                 // Read the first 7 lines in. Q: Is this done for each event or just for the first event?
                 for (unsigned int tr = 0; tr < 7; tr++)
                     {
                      memset(&ifname[0], 0, sizeof(ifname));
                      if (fgets(ifname, sizeof(ifname), fd_file) == NULL) 
                         {
                          std::cout << __FUNCTION__ << " ERROR: Error getting data from file " << std::string(ifname) << ". Skipping channel..." << std::endl;
                          goto fileprob3;
                         }
                      if (tr == 0) npts = sscanf(ifname, "%*s %*s %*u %*s %lf %lf", &time, &volt);
                      else if (tr == 1) npts = sscanf(ifname, "%*s %*s %*lf %*s %lf %lf", &time, &volt);
                      else if (tr == 2) npts = sscanf(ifname, "%*s %*s %u %*s %lf %lf", &trigpoint, &time, &volt);
                      else if (tr == 3) npts = sscanf(ifname, "%*s %*s %lf %*s %lf %lf", &trigtime, &time, &volt);
                      else if (tr == 4) npts = sscanf(ifname, " %*s %lf %lf", &time, &volt);
                      else if (tr == 5) npts = sscanf(ifname, "%*s %*s %*lf %*s %lf %lf", &time, &volt);
                      else {
                            if (std::string(ifname).find("\"FastFrame Count\"") != std::string::npos) npts = sscanf(ifname, "%*s %*s %*u %*s %lf %lf", &time, &volt);
                            else npts = sscanf(ifname, "%lf %lf", &time, &volt);
                           }
                      if (((tr < 2 || (tr > 3 && tr < 7)) && npts == 2) || ((tr > 1 && tr < 4 && npts == 3)))
                         {
                          ((l_t.at(i)).back()).push_back(time); ((l_w.at(i)).back()).push_back(volt);
                          if (tr == 2) l_trigpoint.at(i).back() = trigpoint;
                          else if (tr == 3) (l_trigtime.at(i)).back() = trigtime;
                          np++;
                         }
                      else { 
                            std::cout << __FUNCTION__ << " ERROR: Invalid parameters from ASCI file, skipping event ";
                            if (LGADBase::GetTrackComb()) std::cout << icount-1 << " of file " << fname << std::endl;
                            else std::cout << icount << " of file " << fname << std::endl;
                            ((l_t.at(i)).back()).clear(); ((l_w.at(i)).back()).clear();
                            ((l_t.at(i)).back()).resize((l_npoints.at(r)).at(i), -1); ((l_w.at(i)).back()).resize((l_npoints.at(r)).at(i), -1);
                            (l_trigpoint.at(i)).back() = 0; (l_trigtime.at(i)).back() = -1;
                            while (np < (l_npoints.at(r)).at(i)) { fgets(ifname, sizeof(ifname), fd_file); np++; }
                            icount++; fnd = false;
                            break;
                           }
                     }
                 while (np < (l_npoints.at(r)).at(i))
                       {
                        memset(&ifname[0], 0, sizeof(ifname));
                        if (fgets(ifname, sizeof(ifname), fd_file) == NULL) return false;
                        ch = sscanf(ifname, "%lf %lf", &time, &volt);
                        if (ch == 2)
                           {
                            ((l_t.at(i)).back()).push_back(time);
                            ((l_w.at(i)).back()).push_back(volt);
                            np++;
                           }
                        else {
                              std::cout << __FUNCTION__ << " ERROR: Invalid parameters from ASCI file, skipping event ";
                              if (LGADBase::GetTrackComb()) std::cout << icount - 1 << " of file " << fname << std::endl;
                              else std::cout << icount << " of file " << fname << std::endl;
                              ((l_t.at(i)).back()).clear(); ((l_w.at(i)).back()).clear();
                              ((l_t.at(i)).back()).resize((l_npoints.at(r)).at(i), -1); ((l_w.at(i)).back()).resize((l_npoints.at(r)).at(i), -1);
                              (l_trigpoint.at(i)).back() = 0; (l_trigtime.at(i)).back() = -1;
                              while (np < (l_npoints.at(r)).at(i)) { fgets(ifname, sizeof(ifname), fd_file); np++; }
                              icount++; fnd = false;
                              break;
                             }
                       }
                 if (fnd) icount++;
                 np = 0;
                 if (i == 0 && !LGADBase::GetTrackComb()) l_envtNo.push_back(icount);
                 else if (i == 0 && LGADBase::GetTrackComb()) l_envtNo.push_back(icount-1);
                }  // End of the event loop for this channel of the run
            fclose(fd_file);
           } // End of the loop on all events and chanels for this run
       // loop through all recorded events for all channels of this run and write them to Ntuple
       for (unsigned int gf = 0; gf < (l_w.at(0)).size(); gf++)
           {
            if (LGADBase::GetTrackComb())
               {
                if (LGADBase::GetTrackPackage() == "Cory") LGADBase::CombineTrackCory(l_envtNo.at(gf), trckevntlast, r);
                else if (LGADBase::GetTrackPackage() == "EUTelescope") LGADBase::CombineTrackEUTel(l_envtNo.at(gf), trckevntlast, r);
                else if (LGADBase::GetTrackPackage() == "TBA") LGADBase::CombineTrackTBA(l_envtNo.at(gf), trckevntlast, r);
                else {
                      std::cout << __FUNCTION__ << "          ERROR: Unsopported tracking package " << LGADBase::GetTrackPackage() << std::endl;
                      return false;
                     }
               }
            // Per channel re-organization and event composition
            for (unsigned int dr = 0; dr < (l_nchannels.at(r)).size(); dr++) // loop thoguh available channels
                { 
                 // Find the correct index for this channel at the Ntuple Vectors
                 ch = -1;  // channel index for the ntuple vectors
                 for (unsigned int prt = 0; prt < m_channels.size(); prt++)
                     {
                      if ((l_nchannels.at(r)).at(dr) == m_channels.at(prt))
                         {
                          ch = prt; 
                          break;
                         }
                     }
                 // Re-initalise event vectors and copy event over
                 (m_w.at(ch)).clear(); (m_t.at(ch)).clear();
                 (m_w.at(ch)) = (l_w.at(gf)).at(dr); (m_w.at(ch)) = (l_t.at(gf)).at(dr);
                 m_triggTime.at(ch) = (l_trigtime.at(gf)).at(dr);
                 m_npoints.at(ch) = (l_npoints.at(r)).at(dr);
                 m_srate.at(ch) = (l_snrate.at(r)).at(dr);
                }
            // Per event cross-checks before writing to the file
            m_trigtime = (l_trigpoint.at(0)).at(gf);
            m_event = l_envtNo.at(gf);
            // Before filling the tree one has to look if all channels are populated
            fnd = false;
            for (unsigned int prt = 0; prt < m_channels.size(); prt++)
                {
                 fnd = false;
                 for (unsigned int hj = 0; hj < (l_nchannels.at(r)).size(); hj++)
                     {
                      if ((l_nchannels.at(r)).at(hj) == m_channels.at(prt))
                         {
                          fnd = true;
                          break;
                         }
                     }
                 if (!fnd)
                    {
                     (m_w.at(prt)).clear(); (m_t.at(prt)).clear();
                     (m_w.at(prt)).resize((l_npoints.at(r)).at(0), -1); 
                     (m_w.at(prt)).resize((l_npoints.at(r)).at(0), -1);
                     m_triggTime.at(prt) = -1;
                     m_npoints.at(prt) = -1;
                     m_srate.at(prt) = -1;
                    }
                }
            m_trees.at(0)->Fill();
           }
       std::cout << "Wrote " << l_envtNo.size() << " events from run " << std::string(fname) << std::endl;
       LGADBase::CalcTrigFr(l_trigtime.at(0), m_trigDt, m_trigFr, (l_w.at(0)).size());
      } 
  LGADBase::CloseTrackFiles();
  TListIter l_next_object(gDirectory->GetList());
  l_next_object.Reset();
  TObject* l_obj;
  while ((l_obj=l_next_object()))
        {
         l_obj->Write("", 6);
         delete l_obj;
        }
  m_ofile->Close();

  return true;
}
// --------------------------------------------------------------------------------------------------------------
// Function to read Lecroy oscilloscope binary data from DESY test beams
// LeCroy Binary format is different than the Agilent and Tektronix one that is process by the WriteTestBeamBinary function
// These fiels are writen in an idiotic way tough, storring the same info over and over for each event
// To be merged with the Testbeam binary function that follows, NOT TESTED
// LeCroy has a non-dynamic allocated memory with a fixed number of eents per cycle due to the way the control program is writeen
// To be better writen and incorporated in the TestBeamBinary function later on
bool LGADBase::WriteLecroyBinary(const char* dir, std::vector<TString> names, const char* ext, unsigned int evt1, unsigned int evt2)
{
  // Temporary variables that can be reused
  std::size_t found;
  std::size_t size;
  bool fnd = false;
  char ifname[2048];
  char temp[20480];
  char read[2048];
  char* fileBuf;  // used to store chunks of the binary file when read
  static FILE* fd_file = NULL;
  unsigned int nel = 0;
  int ch = -1;

  std::vector<unsigned int> l_nfiles;  // Files per run
  std::vector<unsigned int> l_runIndx; // indexes of each run for synchronization with tracking files
  std::vector<TString> l_Runs;         // Base name of the run file
  std::vector<TString> l_ext;          // Extentions of the files in the directory
  int nRuns = RunCount(names, l_nfiles, l_Runs, l_ext, l_runIndx, ext);

  // Checking that we only have run files in the input folder
  if (l_ext.size() > 2)
     {
      std::cout << __FUNCTION__ << " ERROR: Multiple file extentions found in input directory: " << dir << std::endl;
      std::cout << __FUNCTION__ << "        Abording......." << std::endl;
      return false;
     }

  // Checking that each run only has 2 files
  for (unsigned int hj = 0; hj < l_nfiles.size(); hj++)
      {
       if (l_nfiles.at(hj) != 2)
          {
           std::cout << __FUNCTION__ << " ERROR: Run with either incomplete or additional files found: " << l_Runs.at(hj) << ", excluding..." << std::endl;
           l_nfiles.erase(l_nfiles.begin() + hj);
           l_Runs.erase(l_Runs.begin() + hj);
           l_runIndx.erase(l_runIndx.begin() + hj);
           nRuns--; hj--;
          }
      }

  // Open all the txt files and get the properties for all channels and runs
  std::vector<std::vector<unsigned int>> l_channels;
  std::vector<int> l_scope;
  std::vector<std::vector< unsigned int>> l_ChScope;
  std::vector < std::vector<unsigned int>> l_EvCount;
  std::vector<std::vector<double>> l_xoffset;
  std::vector<std::vector<double>> l_yoffset;
  std::vector<std::vector<double>> l_yscale;
  std::vector<std::vector<double>> l_snrate;
  std::vector<std::vector<unsigned int>> l_points;
  std::vector<int> l_nEvents;
  // Open the .txt file of each run and recover the information about the number of events
  for (unsigned int a = 0; a < nRuns; a++)
      {
       memset(&ifname[0], 0, sizeof(ifname));
       strncpy(ifname, dir, sizeof(ifname));
       strcat(ifname, l_Runs.at(a));
       strcat(ifname, ".txt");
       fd_file = fopen(ifname, "r");
       // Check if file openned
       if (fd_file == NULL) 
          { 
fileprob1:
           std::cout << __FUNCTION__ << " ERROR: Failed to open txt file: " << ifname << ", skipping run..." << std::endl;
fileprob2:
           l_nfiles.erase(l_nfiles.begin()+a);
           l_Runs.erase(l_Runs.begin()+a);
           l_runIndx.erase(l_runIndx.begin()+a);
           nRuns--; a--;
           fclose(fd_file);
           continue;
          }
       else std::cout << "    Opening text file: " << std::string(ifname) << std::endl;
#ifdef _WIN32
       found = _ftelli64(fd_file);
       _fseeki64(fd_file, 0, 2);
       size = _ftelli64(fd_file);
       _fseeki64(fd_file, found, 0);
#else
       found = ftell(fd_file);
       fseek(fd_file, 0, 2);
       size = ftell(fd_file);
       fseek(fd_file, found, 0);
#endif
       if (size <= 0) goto fileprob1;
       l_channels.push_back(std::vector<unsigned int>());
       l_scope.push_back(0);
       l_ChScope.push_back(std::vector<unsigned int>());
       // Figure out number of recorded channels and their names per run
       for (int i = 0; true; i++)
           {
            // This is problematic since we don't know how the double oscillsocope is handled
            memset(&temp[0], 0, sizeof(temp));
            if (fgets(temp, sizeof(temp), fd_file) == NULL)
               { 
                l_channels.pop_back();
                l_scope.pop_back();
                goto fileprob1;
               }
            memset(&read[0], 0, sizeof(read));
            if (sscanf(temp, "LECROY,WAVERUNNER%u,%14s,%*u.%*u.%*u", &ch, read) == 2) 
               {
                l_scope.back()++;
                memset(&temp[0], 0, sizeof(temp));
                fgets(temp, sizeof(temp), fd_file);
               }
            else if (sscanf(temp, "Channel %d will be saved", &ch) == 1) 
                    {
                     (l_channels.back()).push_back(ch);
                     memset(&temp[0], 0, sizeof(temp));
                     fgets(temp, sizeof(temp), fd_file);
                    }
            else break;
           }
       // Before moving on need to disentangle number of channels and fix channel names in vector
       if (l_scope.at(a) > 1) ArrangeCahnnels(l_channels.at(a), l_ChScope.at(a));
       l_xoffset.push_back(std::vector<double>());
       l_yoffset.push_back(std::vector<double>());
       l_snrate.push_back(std::vector<double>());
       l_yscale.push_back(std::vector<double>());
       l_points.push_back(std::vector<unsigned int>());
       l_EvCount.push_back(std::vector<unsigned int>());
       (l_xoffset.at(a)).resize((l_channels.at(a)).size(), 0.0);
       (l_yoffset.at(a)).resize((l_channels.at(a)).size(), 0.0);
       (l_yscale.at(a)).resize((l_channels.at(a)).size(), 0.0);
       (l_snrate.at(a)).resize((l_channels.at(a)).size(), 0.0);
       (l_points.at(a)).resize((l_channels.at(a)).size(), 0);
       // We are at apoint that we already have a line but it did not meet any of criteria for channel of scope subtrction
       unsigned int iline = 0; // number of lines read
       while (fgets(temp, sizeof(temp), fd_file))
             {
              iline++;
              if (iline < 10) continue; // igonore the first 10 lines
              unsigned int limit = 0;
              for (unsigned int ichan = 0; ichan < (l_channels.at(a)).size(); ichan++)
                  {
                   iline = 0;
                   if (ichan == (l_channels.at(a)).size() - 1) limit = 12;
                   else limit = 49;
                   for (; iline < limit; iline++)
                       {
                        memset(&temp[0], 0, sizeof(temp));
                        fgets(temp, sizeof(temp), fd_file);
                        if (iline == 0 || (iline > 20 && iline < 30))
                           {
                            if (strncmp("NOM_SUBARRAY_COUNT", temp, 18) == 0)                               
                               {
                                if (ichan == 0) (l_EvCount.at(a)).push_back(std::stoi(((std::string)temp).substr(21)));
                                else {
                                      if ((l_EvCount.at(a)).back() != std::stoi(((std::string)temp).substr(21)))
                                         {
                                          std::cout << __FUNCTION__ << " Unequal events recorded in different cahnnels, skipping run " << ifname << std::endl;
                                          l_channels.pop_back(); l_scope.pop_back(); l_EvCount.pop_back(); l_yscale.pop_back();
                                          l_points.pop_back(); l_yoffset.pop_back(); l_xoffset.pop_back(); l_snrate.pop_back();
                                          goto fileprob2;
                                         }
                                    }
                               }
                            else if (strncmp("WAVE_ARRAY_1", temp, 12) == 0)
                                    {
                                     if ((l_EvCount.at(a)).size() == 1) (l_points.at(a)).push_back(std::stoi(((std::string)temp).substr(21))/(l_EvCount.at(a)).back());
                                     else {
                                           if ((l_points.at(a)).at(ichan) != std::stoi(((std::string)temp).substr(21))/(l_EvCount.at(a)).back())
                                              {
                                               std::cout << __FUNCTION__ << " Unequal points for same channel between cycless, skipping run " << ifname << std::endl;
                                               l_channels.pop_back(); l_scope.pop_back(); l_EvCount.pop_back(); l_yscale.pop_back();
                                               l_points.pop_back(); l_yoffset.pop_back(); l_xoffset.pop_back(); l_snrate.pop_back();
                                               goto fileprob2;
                                              }
                                          }
                                    }
                            else if (strncmp ("VERTICAL_OFFSET", temp, 15) == 0) 
                                    {
                                     if ((l_EvCount.at(a)).size() == 1) (l_yoffset.at(a)).push_back(std::stof(((std::string)temp).substr(21)));
                                     else {
                                           if ((l_yoffset.at(a)).at(ichan) != std::stof(((std::string)temp).substr(21)))
                                              {
                                               std::cout << __FUNCTION__ << " Different voltage offsets for same channel between cycless, skipping run " << ifname << std::endl;
                                               l_channels.pop_back(); l_scope.pop_back(); l_EvCount.pop_back(); l_yscale.pop_back();
                                               l_points.pop_back(); l_yoffset.pop_back(); l_xoffset.pop_back(); l_snrate.pop_back();
                                               goto fileprob2;
                                              }
                                          }
                                    }
                            else if (strncmp ("HORIZ_OFFSET", temp, 12) == 0) 
                                             {
                                              if ((l_EvCount.at(a)).size() == 1) (l_xoffset.at(a)).push_back(std::stof(((std::string)temp).substr(21)));
                                              else {
                                                    if ((l_xoffset.at(a)).at(ichan) != std::stof(((std::string)temp).substr(21)))
                                                       {
                                                        std::cout << __FUNCTION__ << " Unequal time offset for same channel between cycless, skipping run " << ifname << std::endl;
                                                        l_channels.pop_back(); l_scope.pop_back(); l_EvCount.pop_back(); l_yscale.pop_back();
                                                        l_points.pop_back(); l_yoffset.pop_back(); l_xoffset.pop_back(); l_snrate.pop_back();
                                                        goto fileprob2;
                                                       }
                                                   }
                                             }
                            else if (strncmp ("HORIZ_INTERVAL", temp, 14) == 0)
                                    {
                                     if ((l_EvCount.at(a)).size() == 1) (l_snrate.at(a)).push_back((double)1/std::stod(((std::string)temp).substr(21)));
                                     else {
                                           if ((l_snrate.at(a)).at(ichan) != ((double)1/std::stod(((std::string)temp).substr(21))))
                                              {
                                               std::cout << __FUNCTION__ << " Different smapling rate for same channel between cycless, skipping run " << ifname << std::endl;
                                               l_channels.pop_back(); l_scope.pop_back(); l_EvCount.pop_back(); l_yscale.pop_back();
                                               l_points.pop_back(); l_yoffset.pop_back(); l_xoffset.pop_back(); l_snrate.pop_back();
                                               goto fileprob2;
                                              }
                                          }
                                    }
                            else if (strncmp ("VERTICAL_GAIN", temp,13) == 0) 
                                    {
                                     if ((l_EvCount.at(a)).size() == 1) (l_yscale.at(a)).push_back(std::stof(((std::string)temp).substr(21)));
                                     else {
                                           if ((l_yscale.at(a)).at(ichan) != std::stof(((std::string)temp).substr(21)))
                                              {
                                               std::cout << __FUNCTION__ << " Different voltage scale same channel between cycless, skipping run " << ifname << std::endl;
                                               l_channels.pop_back(); l_scope.pop_back(); l_EvCount.pop_back(); l_yscale.pop_back();
                                               l_points.pop_back(); l_yoffset.pop_back(); l_xoffset.pop_back(); l_snrate.pop_back();
                                               goto fileprob2;
                                              }
                                          }
                                    }
                           }
                       }
                   }
              iline = 0;
              memset(&temp[0], 0, sizeof(temp));
             }
       l_nEvents.push_back(0);
       for (unsigned int lk = 0; lk < (l_EvCount.back()).size(); lk++) l_nEvents.back() += (l_EvCount.back()).at(lk);
       fclose(fd_file);
      }

  if (m_verbose == 2) 
     {
      for (unsigned int a = 0; a < nRuns; a++)
          {
           std::cout << __FUNCTION__ << " INFO: Total number of events in run " << a << " is " << l_nEvents.at(a) << std::endl;
           std::cout << __FUNCTION__ << " INFO: NUmber of cycles in runs " << a << " is " << l_EvCount.size() << std::endl;
          }
     }

  // Binary file integrity and consistency check
  for (unsigned int a = 0; a < nRuns; a++)
      {
       memset(&ifname[0], 0, sizeof(ifname));
       strncpy(ifname, dir, sizeof(ifname));
       strcat(ifname, l_Runs.at(a));
       if (sizeof(ifname) == 0) 
          { 
           std::cout << __FUNCTION__ << " ERROR: Could not properly form ASCII file name, skipping run...." << std::endl; 
fileprob3:
           l_channels.erase(l_channels.begin() + a); l_scope.erase(l_scope.begin() + a); l_EvCount.erase(l_EvCount.begin() + a);
           l_yscale.erase(l_yscale.begin() + a); l_points.erase(l_points.begin() + a); l_yoffset.erase(l_yoffset.begin() + a);
           l_xoffset.erase(l_xoffset.begin() + a); l_snrate.erase(l_snrate.begin() + a); l_nfiles.erase(l_nfiles.begin() + a);
           l_nEvents.erase(l_nEvents.begin() + a); l_Runs.erase(l_Runs.begin() + a); l_runIndx.erase(l_runIndx.begin() + a);
           nRuns--; a--;
           fclose(fd_file);
           continue;
          }
       nel = 0; // reuse integer to find the fird non txt extention in the extention vector
       for (unsigned int ert = 0; ert < l_ext.size(); ert++) if (l_ext.at(ert) != "txt") { nel = ert;  break; }
       strcat(ifname, ".");
       strcat(ifname, l_ext.at(nel));
       // Check if file openned
       fd_file = fopen(std::string(ifname).c_str(), "rb");
       if (fd_file == NULL)
          { 
           std::cout << __FUNCTION__ << " ERROR: Could not open file " << std::string(ifname) << ". Skipping run..." << std::endl;
           goto fileprob3;
          }
       // Check the file size of the dat file
#ifdef _WIN32
       found = _ftelli64(fd_file);
       _fseeki64(fd_file, 0, 2);
       size = _ftelli64(fd_file);
       _fseeki64(fd_file, found, 0);
#else
       found = ftell(fd_file);
       fseek(fd_file, 0, 2);
       size = ftell(fd_file);
       fseek(fd_file, found, 0);
#endif
       if (size <= 0)
          {
           std::cout << __FUNCTION__ << " ERROR: Failed getting data file size for run " << std::string(ifname) << ". Skipping..." << std::endl;
           goto fileprob3;
          }
       // Get the expected size of the binary file
       size_t exp_size = 0;
       for (unsigned int d = 0; d < l_EvCount.size(); d++)
           {
            // 16 bits as a header, 3 bits at the end of each sequece, number of points, number of event in a cycle
            for (unsigned f = 0; f < (l_channels.at(a)).size(); f++) exp_size += (16 + 3 + (l_points.at(a).at(f))*(l_EvCount.at(a)).at(d));
           }
       if (size != exp_size)
          { 
           std::cout << __FUNCTION__ << " ERROR: Data file size for run " << std::string(ifname) 
                     << " not corresponding to expected value. Skipping..." << std::endl;
           goto fileprob3;
          }
       fclose(fd_file);
      }

  // Fix the number of channels
  if (!FixBaseQuant(l_channels, m_channels, m_nchan)) return false;

  // Create the ntuple
  if (!CreateOutputFile(m_ofdir, m_ofname, m_channels))
     {
      std::cout << __FUNCTION__ << " ERROR: Output Ntuple creation failed: " << m_ofname << std::endl;
      return false;
     }

  // Here we need all the Tracking preparatory actions
  // Tracking Mode related actions:
  // - Check the oscilloscope delay vector and populate it correctly if needed
  // - Open the required data file and get the collection of trees
  // - Add the required parts to the Ntuple
  if (LGADBase::GetTrackComb())
     {
      std::vector<unsigned int> l_runIndx_old = l_runIndx;
      if (!OpenTrackFile(m_Trackdatadir, m_Trackdataname, l_runIndx))
         {
          std::cout << __FUNCTION__ << " ERROR: Failed oppening tracking files in directory: " << m_Trackdatadir << std::endl;
          return false;
         }
      else {
            if (!AddTrckBranches(m_trees.at(0), m_trkpck))
               {
                std::cout << __FUNCTION__ << " Warning: Failed adding the tracking part on the output file: " << m_ofname << std::endl;
                std::cout << "                                       Switching to non-trackimg mode..." << std::endl;
                LGADBase::SetTrackComb(false);
               }
           }
      if (l_runIndx.size() != l_runIndx_old.size())
         {
          std::vector<unsigned int> missing_indxs;
          for (unsigned int we = 0; we < l_runIndx_old.size(); we++)
              {
               fnd = false;
               for (unsigned int yu = 0; yu < l_runIndx.size(); yu++) if (l_runIndx.at(yu) == l_runIndx_old.at(we)) { fnd = true; break; }
               if (!fnd) missing_indxs.push_back(we);
              }
          for (unsigned int rt = 0; rt < missing_indxs.size(); rt++)
              {
               l_channels.erase(l_channels.begin() + rt);
               l_scope.erase(l_scope.begin() + rt);
               l_EvCount.erase(l_EvCount.begin() + rt);
               l_yscale.erase(l_yscale.begin() + rt);
               l_points.erase(l_points.begin() + rt);
               l_yoffset.erase(l_yoffset.begin() + rt);
               l_xoffset.erase(l_xoffset.begin() + rt);
               l_snrate.erase(l_snrate.begin() + rt);
               l_nfiles.erase(l_nfiles.begin() + rt);
               l_nEvents.erase(l_nEvents.begin() + rt);
               l_Runs.erase(l_Runs.begin() + rt);
               l_runIndx.erase(l_runIndx.begin() + rt);
               nRuns--;
              }
         }
     }

  // Start Processing the binary file
  unsigned int totalEvnt = 0;
  std::vector<unsigned int> trckevntlast;
  for (unsigned int a = 0; a < nRuns; a++)
      {
       totalEvnt = 0;
       memset(&ifname[0], 0, sizeof(ifname));
       strncpy(ifname, dir, sizeof(ifname));
       strcat(ifname, l_Runs.at(a));
       nel = 0;
       for (unsigned int ert = 0; ert < l_ext.size(); ert++) if (l_ext.at(ert) != "txt") { nel = ert;  break; }
       strcat(ifname, ".");
       strcat(ifname, l_ext.at(nel));
       // Check if file openned
       fd_file = fopen(std::string(ifname).c_str(), "rb");
       if (fd_file == NULL)
          { 
           std::cout << __FUNCTION__ << " ERROR: Could not open file " << std::string(ifname) << ". Skipping run..." << std::endl;
fileprob4:
           l_channels.erase(l_channels.begin() + a); l_scope.erase(l_scope.begin() + a); l_EvCount.erase(l_EvCount.begin() + a);
           l_yscale.erase(l_yscale.begin() + a); l_points.erase(l_points.begin() + a); l_yoffset.erase(l_yoffset.begin() + a);
           l_xoffset.erase(l_xoffset.begin() + a); l_snrate.erase(l_snrate.begin() + a); l_nfiles.erase(l_nfiles.begin() + a);
           l_nEvents.erase(l_nEvents.begin() + a); l_Runs.erase(l_Runs.begin() + a); l_runIndx.erase(l_runIndx.begin() + a);
           nRuns--; a--;
           fclose(fd_file);
           continue;
          }
       else std::cout << "    Opening binary file: " << std::string(ifname) << std::endl;  
       // Fix the start - stop event with respect to events available in the run
       if (LGADBase::GetTrackComb()) 
          {
           trckevntlast.clear();
           l_nEvents.at(a) = l_nEvents.at(a) - 1;
           trckevntlast.resize(1 + m_nplanes.at(a) + m_nDUT.at(a) + m_nRef.at(a), 0);
          }
       CheckStartStopEvnt(&evt1, &evt2, 0, l_nEvents.at(a));
       nel = 0; // already processed events
       unsigned int offset = 0; // Offset for reading the waveform in the buffer
       unsigned int npts = 0;   // Number of points from all other channels aprat the one read
       // Start copying binary file to buffer and converting events
       for (unsigned int d = 0; d < (l_EvCount.at(a)).size(); d++)
           {
            if (d > 0) nel += (l_EvCount.at(a)).at(d-1);
            // First test of events, to not even perform size estimation and buffer reservation
            if ((LGADBase::GetTrackComb() && (long int)(nel+(l_EvCount.at(a)).at(d)-1) < evt1) ||
                (!LGADBase::GetTrackComb() && (long int)(nel+(l_EvCount.at(a)).at(d)) < evt1)) continue;
            if ((LGADBase::GetTrackComb() && (long int)(nel+(l_EvCount.at(a)).at(d)-1) > evt2) ||
                (!LGADBase::GetTrackComb() && (long int)(nel+(l_EvCount.at(a)).at(d)) > evt2))
               {
jumpnextrun:
                break;
               }
            size = 0;
            for (unsigned f = 0; f < (l_channels.at(a)).size(); f++) size += (16+3+(l_points.at(a)).at(f)*(l_EvCount.at(a)).at(d));
            fileBuf = (char*)malloc(size);
            size_t found = fread(fileBuf, size, 1, fd_file);
            if (fileBuf == NULL || found != 1) 
               { 
                std::cout << __FUNCTION__ << " ERROR: Failed reading buffer for run " << std::string(ifname) << " skipping...." << std::endl;
                if (d == 0) goto fileprob4;
                else break;
               }
            for (unsigned int etr = 0; etr < (l_EvCount.at(a)).at(d); etr++)
                { 
                 nel += etr;
                 if (LGADBase::GetTrackComb())
                    { 
                     if ((long int)nel < evt1+1) continue;
                     if ((long int)nel > evt2+1) goto jumpnextrun;
                     m_event = nel - 1;
                     LGADBase::ProgressBar((evt2 - evt1)* a + nel-1, (evt2 - evt1)* nRuns);
                    }
                 else {
                       if ((long int)nel < evt1) continue;
                       if ((long int)nel > evt2) goto jumpnextrun;
                       m_event = nel;
                       LGADBase::ProgressBar((evt2 - evt1)* a + nel, (evt2 - evt1)* nRuns);
                      }
                 for (unsigned int gtr = 0; gtr < (l_channels.at(a)).size(); gtr++)
                     {
                      // Find the correct index for the vector of the entuple, use ch variable
                      for (unsigned int prt = 0; prt < m_channels.size(); prt++)
                          {
                           if ((l_channels.at(a)).at(gtr) == m_channels.at(prt)) { ch = prt; break; }
                          }
                      m_scale.at(ch) = (l_yscale.at(a)).at(gtr);
                      m_npoints.at(ch) = (l_points.at(a)).at(gtr);
                      m_srate.at(ch) = (l_snrate.at(a)).at(gtr);
                      m_w.at(ch).clear();
                      m_w.at(ch).reserve((l_points.at(a)).at(gtr));
                      m_t.at(ch).clear();
                      m_t.at(ch).reserve((l_points.at(a)).at(gtr));
                      m_scope.at(ch) = (l_ChScope.at(a)).at(gtr);
                      // The output format of the data file seems to probably be following the structure:
                      // channels x {header + all waveforms + footer} with header = 16 and footer = 3 and no separator between waveforms
                      // The file is not organized in an event by event but in a channel by channel approach and one has to jump around
                      offset = 0; npts = 0;
                      for (unsigned int yr = 0; yr < gtr; yr++) npts += (l_points.at(a)).at(yr);
                      offset = (gtr+1)*16 + gtr*3 + npts*(l_EvCount.at(a)).at(d) + etr*(l_points.at(a)).at(gtr);
                      for (unsigned int w = 0; w < (l_points.at(a)).at(gtr); w++)
                          {
                           if ((w > ((l_points.at(a)).at(gtr) - 4)) && (((fileBuf[w+offset]*(l_yscale.at(a)).at(gtr)) - (l_yoffset.at(a)).at(gtr)) > 0.04)) break;
                           else if ((w < 10) && (((fileBuf[w+offset]*(l_yscale.at(a)).at(gtr) - (l_yoffset.at(a)).at(gtr)) > 0.04) ||
                                                 ((fileBuf[w+offset]*(l_yscale.at(a)).at(gtr) - (l_yoffset.at(a)).at(gtr)) < -0.04))) continue;
                           else { // amplitude of signal is given in volts                            
                                 m_w.at(ch).push_back(fileBuf[w+offset]*(l_yscale.at(a)).at(gtr) - (l_yoffset.at(a)).at(gtr));
                                 m_t.at(ch).push_back(w*((double)1/(l_snrate.at(a)).at(gtr))+(l_xoffset.at(a)).at(gtr));
                                }
                          }    
                     }
                 // Before filling the tree have to populate any not included channel
                 for (unsigned int prt = 0; prt < m_channels.size(); prt++)
                     {
                      fnd = false;
                      for (unsigned int hj = 0; hj < (l_channels.at(a)).size(); hj++)
                          {
                           if (l_channels.at(a).at(hj) == m_channels.at(prt)) { fnd = true; break; }
                          }
                      if (!fnd)
                         {
                          (m_w.at(prt)).clear(); (m_t.at(prt)).clear();
                          (m_w.at(prt)).resize((l_points.at(a)).at(0), -1); 
                          (m_w.at(prt)).resize((l_points.at(a)).at(0), -1);
                          m_scale.at(prt) = -1;
                          m_npoints.at(prt) = -1;
                          m_srate.at(prt) = -1;
                         }
                     }
                 if (LGADBase::GetTrackComb())
                    {
                     if (LGADBase::GetTrackPackage() == "Cory") LGADBase::CombineTrackCory(m_event, trckevntlast, a);
                     else if (LGADBase::GetTrackPackage() == "EUTelescope") LGADBase::CombineTrackEUTel(m_event, trckevntlast, a);
                     else if (LGADBase::GetTrackPackage() == "TBA") LGADBase::CombineTrackTBA(m_event, trckevntlast, a);
                     else {
                           std::cout << __FUNCTION__ << "          ERROR: Unsopported tracking package " << LGADBase::GetTrackPackage() << std::endl;
                           return false;
                          }
                    }
                 totalEvnt++;
                 m_trees.at(0)->Fill(); // fill every event
                }
            free(fileBuf);
           }
       std::cout << "Wrote " << totalEvnt << " events from run " << std::string(ifname) << std::endl;
       fclose(fd_file);
      }

  LGADBase::CloseTrackFiles();
  TListIter l_next_object(gDirectory->GetList());
  l_next_object.Reset();
  TObject* l_obj;
  while ((l_obj=l_next_object()))
        {
         l_obj->Write("", 6);
         delete l_obj;
        }
  m_ofile->Close();

  return true;
}
// --------------------------------------------------------------------------------------------------------------
// Function to read binary files from multiple oscilloscopes 
// Works with Agilent and Tektronix models and supports automatic
// scope detection and multi-file processing for combined runs.
bool LGADBase::WriteTestBeamBinary(const char* dir, std::vector<TString> names, const char* ext, unsigned int evt1, unsigned int evt2)
{
  // Temporary variables that can be reused
  std::size_t found;
  std::size_t size;
  std::size_t lPrevPos;
  bool fnd = false;
  static FILE* fd_file = NULL;  // internal
  char* fileBuf;  // used to store chunks of the binary file when read
  char ifname[2048], temp[20480];
  char time0[1024], time1[1024], time2[1024];
  char sn0[1024], sn1[1024], sn2[1024];
  unsigned int nel = 0;
  unsigned int nev = 0;
  unsigned int nrec = 0;
  unsigned int scope = 0;

  std::vector<unsigned int> l_nfiles;  // Files per run
  std::vector<unsigned int> l_runIndx; // indexes of each run for synchronization with tracking files
  std::vector<TString> l_Runs;         // Base name of the run file
  std::vector<TString> l_ext;          // Extentions of the files in the directory
  int nRuns = RunCount(names, l_nfiles, l_Runs, l_ext, l_runIndx, ext);

  // Checking that we only have run files in the input folder
  if (l_ext.size() > 2)
     {
      std::cout << __FUNCTION__ << " ERROR: Multiple file extentions found in input directory: " << dir << std::endl;
      std::cout << __FUNCTION__ << "        Abording......." << std::endl;
      return false;
     }

  // Checking that each run only has 2 files
  for (unsigned int hg = 0; hg < l_nfiles.size(); hg++)
      {
       if (l_nfiles.at(hg) != 2)
          {
           std::cout << __FUNCTION__ << " ERROR: Run with either incomplete or additional files found: " << l_Runs.at(hg) << ", excluding..." << std::endl;
           l_nfiles.erase(l_nfiles.begin() + hg);
           l_Runs.erase(l_Runs.begin() + hg);
           l_runIndx.erase(l_runIndx.begin() + hg);
           nRuns--; hg--;
          }
      }

  // Open all the header files and get the properties for all channels and runs
  // Define oscilloscope type TestBeamBin1 = Agilent Scope
  //                          TestBeamBin2 = Tektronix Scope
  m_instrument = TestBeamBin;
  std::vector<std::vector<unsigned int>> l_channels;   // Channels within a run
  std::vector<std::vector<unsigned int>> l_ChScope;    // Oscilloscope at which a channel belongs to
  std::vector<unsigned int> l_nScopes;                 // Oscilloscopes per run, important for delays
  std::vector<std::vector<unsigned int>> l_EvCount;  // Numer of events per sequence and per run
  std::vector<int> l_nEvents;                          // Total number of events in the run
  std::vector<std::vector<double>> l_xoffset;
  std::vector<std::vector<double>> l_yoffset;
  std::vector<std::vector<double>> l_yscale;
  std::vector<std::vector<double>> l_snrate;
  std::vector<std::vector<double>> l_triggtimeVec;
  std::vector<std::vector<unsigned int>> l_points;
  std::vector<std::vector<unsigned int>> l_dattype;
  for (unsigned int a = 0; a < nRuns; a++)
      {
       memset(&ifname[0], 0, sizeof(ifname));
       strncpy(ifname, dir, sizeof(ifname));
       strcat(ifname, l_Runs.at(a));
       strcat(ifname, ".txt");
       fd_file = fopen(ifname, "r");
       // Check if file openned
       if (fd_file == NULL)
          {
           memset(&ifname[0], 0, sizeof(ifname));
           strncpy(ifname, dir, sizeof(ifname));
           strcat(ifname, l_Runs.at(a));
           strcat(ifname, "_head.txt");
           fd_file = fopen(ifname, "r");
           if (fd_file == NULL) 
              { 
fileprob1:
               std::cout << __FUNCTION__ << " ERROR: Failed to open txt file: " << std::string(ifname) << ", skipping run..." << std::endl;
fileprob2:
               l_nfiles.erase(l_nfiles.begin()+a);
               l_Runs.erase(l_Runs.begin() + a);
               l_runIndx.erase(l_runIndx.begin() + a);
               nRuns--; a--;
               fclose(fd_file);
               continue;
              }
           else { std::cout << "    Opening text file: " << std::string(ifname) << std::endl; m_instrument = TestBeamBin2; }
          }
       else { std::cout << "    Opening text file: " << std::string(ifname) << std::endl; m_instrument = TestBeamBin1; }
#ifdef _WIN32
       found = _ftelli64(fd_file); 
       _fseeki64(fd_file, 0, 2);     // carefull the fseek third argument is as follows: 0 -> start of file, 1-> current position, 2-> end of file
       size = _ftelli64(fd_file);
       _fseeki64(fd_file, found, 0); // carefull the fseek third argument is as follows: 0 -> start of file, 1-> current position, 2-> end of file
#else
       found = ftell(fd_file); 
       fseek(fd_file, 0, 2);         // carefull the fseek third argument is as follows: 0 -> start of file, 1-> current position, 2-> end of file
       size = ftell(fd_file);  
       fseek(fd_file, found, 0);     // carefull the fseek third argument is as follows: 0 -> start of file, 1-> current position, 2-> end of file
#endif
       if (size <= 0) goto fileprob1;
       scope = 0;
       l_channels.push_back(std::vector<unsigned int>());
       l_ChScope.push_back(std::vector<unsigned int>());
       l_yscale.push_back(std::vector<double>());
       l_yoffset.push_back(std::vector<double>());
       l_points.push_back(std::vector<unsigned int>());
       l_snrate.push_back(std::vector<double>());
       l_EvCount.push_back(std::vector<unsigned int>());
       if (m_instrument == TestBeamBin2) 
          {
           l_xoffset.push_back(std::vector<double>());
           l_dattype.push_back(std::vector<unsigned int>());
           l_triggtimeVec.push_back(std::vector<double>());
          }
       memset(&temp[0], 0, sizeof(temp));
       fnd = false;
       bool start = false;
       // Read first lines to determine number of channels
       while (!fnd && fgets(temp, sizeof(temp), fd_file))
             {
              nrec = 0;
              unsigned int npts, xorg;
              float yinc, yorg;
              double sr;
              memset(&sn2[0], 0, sizeof(sn2));
              if (m_instrument == TestBeamBin1)
                 {
                  memset(&time2[0], 0, sizeof(time2));
                  nel = sscanf(temp, "%*i,%*i,%i,%*i,%lf,%*f,%*f,%f,%*f,%*f,%*d,%*f,%*f,%*f,%f,%*[^,],%[^,],%[^,],%*d,%*d,%*d,%*d,%*f,%*f,%d", &npts, &sr, &yinc, &yorg, time2, sn2, &nev);
                  if (!(nel == 7 && npts > 0 && nev > 0 && sr > 0))
                     { 
                      std::cout << __FUNCTION__ << " ERROR: Invalid parameters from ASCI file, skippimg run " << l_Runs.at(a) << "..." << std::endl;
                      l_channels.pop_back(); l_ChScope.pop_back(); l_yscale.pop_back(); 
                      l_yoffset.pop_back(); l_points.pop_back(); l_snrate.pop_back(); l_EvCount.pop_back();                  
                      goto fileprob2;
                     }
                  if (!start)
                     {
                      (l_yscale.at(a)).push_back(yinc);
                      (l_yoffset.at(a)).push_back(yorg);
                      (l_points.at(a)).push_back(npts);
                      (l_snrate.at(a)).push_back((double)1/sr);
                      if (strncmp(sn1, sn2, 1024) != 0 && strncmp(time1, time2, 1024) != 0 && strncmp(sn0, sn2, 1024) != 0) // every time we change oscilloscope
                         {
                          memcpy(sn1, sn2, 1024);
                          memcpy(time1, time2, 1024);
                          if (scope == 0) // first oscilloscope
                             {
                              (l_EvCount.at(a)).push_back(nev);
                              memcpy(sn0, sn2, 1024);
                              memcpy(time0, time2, 1024);
                             }
                          (l_channels.at(a)).push_back((scope*4) + 1);
                          scope++;
                         }
                      else if (strncmp(sn1, sn2, 1024) == 0 && strncmp(time1, time2, 1024) == 0 && nev == (l_EvCount.at(a)).at(0))
                              {
                               (l_channels.at(a)).push_back((l_channels.at(a)).back() + 1); // Same scope
                              }
                      else if (strncmp(sn0, sn2, 1024) == 0 && strncmp(time0, time2, 1024) != 0 && npts == (l_points.at(a)).at(0) && // back to first scope
                               ceil(1/sr) == (l_snrate.at(a)).at(0) && yinc == (l_yscale.at(a)).at(0) && yorg == (l_yoffset.at(a)).at(0)) start = true;
                      else {
                            std::cout << __FUNCTION__ << " ERROR: Not able to define number of used channels, skipping run "  << l_Runs.at(a) << "..." << std::endl;
                            l_channels.pop_back(); l_ChScope.pop_back(); l_yscale.pop_back();
                            l_yoffset.pop_back(); l_points.pop_back(); l_snrate.pop_back(); l_EvCount.pop_back();
                            goto fileprob2;
                           }
                      (l_ChScope.at(a)).push_back(scope);
                     }
                  if ((l_channels.at(a)).size() != 0 && (nrec % (l_channels.at(a)).size() == 0)) (l_EvCount.at(a)).push_back(nev);
                  nrec++;  
                 }
              else if (m_instrument == TestBeamBin2)
                      {
                       memset(&sn1[0], 0, sizeof(sn1));
                       memset(&sn0[0], 0, sizeof(sn0));
                       nel = sscanf(temp, "%i;%*i;%*[^;];%[^;];%*[^;];%[^;];%[^;];%i;%*[^;];%*[^;];%*[^;];%lf;%*lf;%i;%*[^;];%f;%*f;%f;%*[^;];%*[^;];%*f;%*f;%*f;%*i", &nev, sn0, sn2, sn1, &npts, &sr, &xorg, &yinc, &yorg);
                       // Either it fails becasue the fromat in incorect (no channels read) or because it progressed to the date/time part of the file
                       if (!(nel == 9 && npts > 0 && nev > 0))
                          {
                           if ((l_channels.at(a)).size() == 0) // failed right from the start, no channels were read
                              {
                               std::cout << __FUNCTION__ << " ERROR: Invalid parameters from header file, skipping run " << l_Runs.at(a) << "..." << std::endl;
                               l_channels.pop_back(); l_ChScope.pop_back(); l_yscale.pop_back(); l_dattype.pop_back(); l_triggtimeVec.pop_back();
                               l_yoffset.pop_back(); l_points.pop_back(); l_snrate.pop_back(); l_EvCount.pop_back(); l_xoffset.pop_back(); 
                               goto fileprob2;
                              }
                           else {
                                 fnd = true; // we finhsed reading the part containing channel info, move along to timestamps
                                 nel = 0;
                                 for (unsigned int k = 0; k < sizeof(temp); k++) { if (temp[k] != '\0') nel++; else break; }
#ifdef _WIN32
                                 if ((nel +1) == (int)(sizeof(temp))) _fseeki64(fd_file, _ftelli64(fd_file)- nel -4, 0);
                                 else if ((nel +1) < (int)(sizeof(temp))) _fseeki64(fd_file, _ftelli64(fd_file)- nel -5, 0);
#else
                                 if ((nel +1) == (int)(sizeof(temp))) fseek(fd_file, ftell(fd_file)- nel, 0);
                                 else if ((nel +1) < (int)(sizeof(temp))) fseek(fd_file, ftell(fd_file)- nel, 0);
#endif  
                                 nel = 0;
                                }
                          } 
                       (l_yscale.at(a)).push_back(yinc);
                       (l_yoffset.at(a)).push_back(yorg);
                       (l_points.at(a)).push_back(npts);
                       (l_snrate.at(a)).push_back((double)1 / sr);
                       (l_xoffset.at(a)).push_back(xorg);
                       // Unfold waveform data type following table:
                       unsigned int type = 0;
                       if (nev == 1) type = 100;
                       else if (nev == 2) type = 200;
                       else {
                             std::cout << __FUNCTION__ << " ERROR: COuld not determine waveform byte size for run " << l_Runs.at(a) << ". Skipping..." << std::endl;
                             l_channels.pop_back(); l_ChScope.pop_back(); l_yscale.pop_back(); l_dattype.pop_back(); l_triggtimeVec.pop_back();
                             l_yoffset.pop_back(); l_points.pop_back(); l_snrate.pop_back(); l_EvCount.pop_back(); l_xoffset.pop_back();
                             goto fileprob2;
                            }
                       if (strcmp(sn2, "MSB") == 0) type += 10;
                       else if (strcmp(sn2, "LSB") == 0) type += 20;
                       else {
                             std::cout << __FUNCTION__ << " ERROR: Could not determine waveform endianess for run " << l_Runs.at(a) << ". Skipping..." << std::endl;
                             l_channels.pop_back(); l_ChScope.pop_back(); l_yscale.pop_back(); l_dattype.pop_back(); l_triggtimeVec.pop_back();
                             l_yoffset.pop_back(); l_points.pop_back(); l_snrate.pop_back(); l_EvCount.pop_back(); l_xoffset.pop_back();
                             goto fileprob2;
                            }
                       if (strcmp(sn0, "RI") == 0) type += 1;
                       else if (strcmp(sn0, "RP") == 0) type += 2;
                       else {
                             std::cout << __FUNCTION__ << " ERROR: Could not determine wavwform integer representation for run " << l_Runs.at(a) << ". Skipping..." << std::endl;
                             l_channels.pop_back(); l_ChScope.pop_back(); l_yscale.pop_back(); l_dattype.pop_back(); l_triggtimeVec.pop_back();
                             l_yoffset.pop_back(); l_points.pop_back(); l_snrate.pop_back(); l_EvCount.pop_back(); l_xoffset.pop_back();
                             goto fileprob2;
                            }
                       (l_dattype.at(a)).push_back(type);
                       sscanf(sn1, "\"Ch%u, %*[^\"]\"", &type); // Reusing the integer here
                       (l_channels.at(a)).push_back(type);
                      }
              memset(&temp[0], 0, sizeof(temp));
             }
       // Part only relevent  for the tektronix scope
       if (m_instrument == TestBeamBin2)
          {
           scope = ArrangeCahnnels(l_channels.at(a), l_ChScope.at(a));
           double trigtime_first = 0.0; // time of the first event
           m_trigtime = 0.0;  // reuse of class variable, no meaning in this sence, used as a temporary value
           nel = 0; nev = 0;  // reusing integers from before, decalred at the begining and still in scope
           unsigned int year, month;
           struct tm timeinfo;
           long double milisec = 0;
           long double milisec1 = 0;
#ifdef _WIN32
           _fseeki64(fd_file, _ftelli64(fd_file) - 2*(scope-1), 0);
           found = _ftelli64(fd_file);
           if (scope > 1) lPrevPos = found;
           while (_ftelli64(fd_file) < (int)size)
#else
           fseek(fd_file, ftell(fd_file), 0); // LXplus mod
           found = ftell(fd_file);
           if (scope > 1) lPrevPos = found;
           while (ftell(fd_file) < (int)size)
#endif
                 {
                  if (!(fread(temp, 1, 64, fd_file) == 64)) // do not read more than 64 charcters, you should have found the next event if it exists within
                     {
                      std::cout << __FUNCTION__ << " ERROR: Read Error while pharsing date/time file for run " << l_Runs.at(a) << ". Skipping..." << std::endl;
                      l_channels.pop_back(); l_ChScope.pop_back(); l_yscale.pop_back(); l_dattype.pop_back(); l_triggtimeVec.pop_back();
                      l_yoffset.pop_back(); l_points.pop_back(); l_snrate.pop_back(); l_EvCount.pop_back(); l_xoffset.pop_back();
                      goto fileprob2;
                     } 
                  for (int i = 1; i < 64; i++) 
                      {         
                       if (temp[i] == ',' || temp[i] == '\n' || temp[i] == '"') // search for the comma separator or end of line
                          {
                           nrec = 99; // for creating the string to read in later, reusing the integer
                           if (temp[i] == '"' && temp[i + 1] == '"') nrec = 0;
                           else if (temp[i] == ',' || temp[i] == '\n') nrec = 1;
                           if (nrec != 99)
                              {
                               // Move ittertor to process the next string
                               found += i + 1;
                               if (scope > 1 && temp[i] == '\n')
                                  {
                                   found += scope*(found - lPrevPos) + 1;
                                   lPrevPos = found;
                                  }
                               // Restructure read string to contain exactly one event
#ifdef _WIN32                        
                               memset(&temp[i - nrec + 1], 0, sizeof(temp) - i + nrec - 1);
                               _fseeki64(fd_file, found, 0);
#else
                               memset(&temp[i - nrec + 2], 0, sizeof(temp) - i + nrec - 1);  // LXplus mod
                               fseek(fd_file, found, 0);
#endif                     
                               break; // break to spare some time
                              }
                          }
                      } 
               strcpy(temp, (LGADBase::trim(std::string(temp), "\"")).c_str());
               nel = sscanf(temp, "%i: %2u.%2u.%4u.%2u:%2u:%2u.%Lf", &nrec, &timeinfo.tm_mday, &month, &year, &timeinfo.tm_hour, &timeinfo.tm_min, &timeinfo.tm_sec, &milisec);
               if (!(nel == 8 && nrec > 0 && milisec > 0 && year > 0))
                  { 
                   std::cout << __FUNCTION__ << " ERROR: Invalid parameters from header file for event date/time for run " << l_Runs.at(a) << ". Skipping..." << std::endl;
                   l_channels.pop_back(); l_ChScope.pop_back(); l_yscale.pop_back(); l_dattype.pop_back(); l_triggtimeVec.pop_back();
                   l_yoffset.pop_back(); l_points.pop_back(); l_snrate.pop_back(); l_EvCount.pop_back(); l_xoffset.pop_back();
                   goto fileprob2;
                  }               
               timeinfo.tm_year = year - 1900;
               timeinfo.tm_mon = month - 1;
               m_trigtime = mktime(&timeinfo);
               if (nrec == 1) // trigger time from first event in the sequence
                  {
                   nev = 0;
                   trigtime_first = m_trigtime;
                   while (milisec > 1) { milisec = milisec / ((double)10); nev++; }
                   milisec1 = milisec;
                  }
               m_trigtime = (long long int)m_trigtime - (long long int)trigtime_first;
               if (m_trigtime < 0) { while (m_trigtime < 0) m_trigtime = m_trigtime + 3600; }
               else if (m_trigtime >= 3600 && nrec < 12) {trigtime_first = trigtime_first + 3600; m_trigtime = m_trigtime-3600;}
               if (nrec == 1) m_trigtime = m_trigtime + ((long double)milisec - (long double)milisec1);
               else m_trigtime = m_trigtime + ((long double)(milisec / pow(10, nev)) - (long double)milisec1);
               (l_triggtimeVec.at(a)).push_back(m_trigtime);
               if (nrec != 1) (l_EvCount.at(a)).back() = nrec; // update the number of events for the sequence
               else (l_EvCount.at(a)).push_back(nrec); // add an entry to the vector, start a new sequence
               memset(&temp[0], 0, sizeof(temp)); // clear the char buffer before reusing it
              }
          }
       if (m_instrument == TestBeamBin1)
          {
           if (l_EvCount.size() != (nrec / (l_channels.at(a)).size()))
              {
               std::cout << __FUNCTION__ << " ERROR: Channel and sequences numbers not correctly decoded, skipping run " << l_Runs.at(a) << "..." << std::endl;
               l_channels.pop_back(); l_ChScope.pop_back(); l_yscale.pop_back();
               l_yoffset.pop_back(); l_points.pop_back(); l_snrate.pop_back(); l_EvCount.pop_back();
               goto fileprob2;
              }
          }
       fclose(fd_file);
      }

   // Fix the number of oscilloscopes per run
   for (unsigned int jh = 0; jh < nRuns; jh++)
       {
       l_nScopes.push_back(1);
        unsigned int tre = (l_ChScope.at(jh)).at(0);
        for (unsigned int sd = 1; sd < (l_ChScope.at(jh)).size(); jh++)
            {
             if ((l_ChScope.at(jh)).at(sd) != tre) l_nScopes.back()++; tre = (l_ChScope.at(jh)).at(sd);
            }
       }
   
  // File size check and calculations of the number of bytes and digits of data decoding
  std::vector<std::vector<unsigned int>> l_bytes;
  std::vector<std::vector<unsigned int>> l_digits;
  // Start events processing from the binary files
  for (unsigned int a = 0; a < l_runIndx.size(); a++)
      {
       memset(&ifname[0], 0, sizeof(ifname));
       strncpy(ifname, dir, sizeof(ifname));
       strcat(ifname, l_Runs.at(a));
       if (sizeof(ifname) == 0) 
          { 
           std::cout << __FUNCTION__ << " ERROR: Could not properly form ASCII file name, skipping run...." << std::endl; 
fileprob3:
           l_channels.erase(l_channels.begin() + a); l_ChScope.erase(l_ChScope.begin() + a); l_yscale.erase(l_yscale.begin() + a);
           l_yoffset.erase(l_yoffset.begin() + a); l_points.erase(l_points.begin() + a); l_snrate.erase(l_snrate.begin() + a); 
           l_EvCount.erase(l_EvCount.begin() + a); l_nScopes.erase(l_nScopes.begin() + a);
           if (m_instrument == TestBeamBin2) 
              { 
               l_dattype.erase(l_dattype.begin() + a); 
               l_triggtimeVec.erase(l_triggtimeVec.begin() + a); 
               l_xoffset.erase(l_xoffset.begin() + a);
              }
           l_nfiles.erase(l_nfiles.begin() + a); l_Runs.erase(l_Runs.begin() + a); l_runIndx.erase(l_runIndx.begin() + a);
           nRuns--; a--;
           fclose(fd_file);
           continue;
          }
       nel = 0;
       for (unsigned int ert = 0; ert < l_ext.size(); ert++) if (l_ext.at(ert) != "txt") { nel = ert;  break; }
       strcat(ifname, ".");
       strcat(ifname, l_ext.at(nel));
       // Check if file openned
       fd_file = fopen(std::string(ifname).c_str(), "rb");
       if (fd_file == NULL)
          { 
           std::cout << __FUNCTION__ << " ERROR: Could not open file " << std::string(ifname) << ". Skipping run..." << std::endl;
           goto fileprob3;
          }
       else std::cout << "    Opening binary file: " << std::string(ifname) << std::endl;   
       // Get the file size of the dat file
#ifdef _WIN32
       found = _ftelli64(fd_file);
       _fseeki64(fd_file, 0, 2); 
       size = _ftelli64(fd_file);
       _fseeki64(fd_file, found, 0);
#else
       found = ftell(fd_file);
       fseek(fd_file, 0, 2);
       size = ftell(fd_file);
       fseek(fd_file, found, 0);
#endif
       if (size <= 0)
          {
           std::cout << __FUNCTION__ << " ERROR: Failed getting data file size for run " << std::string(ifname) << ". Skipping..." << std::endl;
           goto fileprob3;
          }
       // Check file size of dat file
       found = 0; // re-use the integer
       // Local vectors to store number of bytes and number of points digit
       if (m_instrument == TestBeamBin2)
          {
           l_bytes.push_back(std::vector<unsigned int>());
           l_digits.push_back(std::vector<unsigned int>());
           (l_bytes.at(a)).reserve((l_channels.at(a)).size());
           (l_digits.at(a)).reserve((l_channels.at(a)).size());
          }
       for (unsigned f = 0; f < (l_channels.at(a)).size(); f++)
           { 
            if (m_instrument == TestBeamBin1) found += ((l_points.at(a)).at(f) + 1) * 2 * l_nEvents.at(a);
            // Structure of data file for tektronix oscilloscope: 
            // Per sequence: #<number_of_digits><number_of_bytes><waveforme>; ..... ; next event ; last event
            else if (m_instrument == TestBeamBin2)
                    {
                     nel = 0; // reuse the integer
                     if ((l_dattype.at(a)).at(f) < 200) (l_bytes.at(a)).push_back(1);
                     else if ((l_dattype.at(a)).at(f) > 200) (l_bytes.at(a)).push_back(2);
                     else {
                           std::cout << __FUNCTION__ << " ERROR: Incorect data size for channel " << (l_channels.at(a)).at(f) << " of run " 
                                     << std::string(ifname) << ". Skipping...." << std::endl;
                           l_bytes.erase(l_bytes.begin() + a);
                           goto fileprob3;
                          }
                     nrec = (l_points.at(a)).at(f)*(l_bytes.at(a)).at(f); // reuse the integer
                     do { nrec /= 10; nel++; } while (nrec != 0);
                     (l_digits.at(a)).push_back(nel);
                     found += (((l_points.at(a)).at(f) * (l_bytes.at(a)).at(f)) + 3 + (l_digits.at(a)).at(f)) * l_nEvents.at(a) - (l_EvCount.at(a)).size();
                    }
           }
       if (size != found)
          { 
           std::cout << __FUNCTION__ << " ERROR: Data file size for run " << std::string(ifname) << " not corresponding to expected value!" << std::endl; 
           if (size > found) 
              {
               std::cout << "                                     Size bigger than estimated, will try to recover expected events from header but CHECK DATA!!!!" << std::endl;
               std::cout << "                                     Clculated size: "<< found << " Real siez: " << size << " Differencce: " << size - found << std::endl;
              }
           else { 
                 std::cout << "                                     Skipping..." << std::endl; 
                 if (m_instrument == TestBeamBin2)
                    {
                     l_bytes.erase(l_bytes.begin() + a);
                     l_digits.erase(l_digits.begin() + a);
                    }
                 goto fileprob3; 
                }
          }
       fclose(fd_file);
      }

  // Fix the number of channels
  if (!FixBaseQuant(l_channels, m_channels, m_nchan)) return false;

  // Create the ntuple
  if (!CreateOutputFile(m_ofdir, m_ofname, m_channels))
     {
      std::cout << __FUNCTION__ << " ERROR: Output Ntuple creation failed: " << m_ofname << std::endl;
      return false;
     }

  // Tracking Mode related actions:
  // - Reduce the number of events by one to account for the initial TLU isseud trigger at the start of each run
  // - Check the oscilloscope delay vector and populate it correctly if needed
  // - Open the required data file and get the collection of trees
  // - Add the required parts to the Ntuple
  if (LGADBase::GetTrackComb())
     {
      if (l_nScopes.size() > m_oscdel.size()) while (m_oscdel.size() < l_nScopes.size()) m_oscdel.push_back(-99.);
      std::vector<unsigned int> l_runIndx_old = l_runIndx;
      if (!OpenTrackFile(m_Trackdatadir, m_Trackdataname, l_runIndx))
         {
          std::cout << __FUNCTION__ << " ERROR: Failed oppening tracking files in directory: " << m_Trackdatadir << std::endl;
          return false;
         }
      else {
            if (!AddTrckBranches(m_trees.at(0), m_trkpck))
               {
                std::cout << __FUNCTION__ << " Warning: Failed adding the tracking part on the output file: " << m_ofname << std::endl;
                std::cout << "                                       Switching to non-trackimg mode..." << std::endl;
                LGADBase::SetTrackComb(false);
               }
           }
      if (l_runIndx.size() != l_runIndx_old.size())
         {
          std::vector<unsigned int> missing_indxs;
          for (unsigned int we = 0; we < l_runIndx_old.size(); we++)
              {
               fnd = false;
               for (unsigned int yu = 0; yu < l_runIndx.size(); yu++) if (l_runIndx.at(yu) == l_runIndx_old.at(we)) { fnd = true; break; }
               if (!fnd) missing_indxs.push_back(we);
              }
          for (unsigned int rt = 0; rt < missing_indxs.size(); rt++)
              {
               l_channels.erase(l_channels.begin() + missing_indxs.at(rt)); l_ChScope.erase(l_ChScope.begin() + missing_indxs.at(rt));
               l_yscale.erase(l_yscale.begin() + missing_indxs.at(rt)); l_yoffset.erase(l_yoffset.begin() + missing_indxs.at(rt));
               l_points.erase(l_points.begin() + missing_indxs.at(rt)); l_snrate.erase(l_snrate.begin() + missing_indxs.at(rt));
               l_EvCount.erase(l_EvCount.begin() + missing_indxs.at(rt)); l_nScopes.erase(l_nScopes.begin() + missing_indxs.at(rt));
               l_nfiles.erase(l_nfiles.begin() + missing_indxs.at(rt)); l_Runs.erase(l_Runs.begin() + missing_indxs.at(rt));
               l_runIndx.erase(l_runIndx.begin() + missing_indxs.at(rt));
               if (m_instrument == TestBeamBin2)
                  { 
                   l_dattype.erase(l_dattype.begin() + missing_indxs.at(rt));
                   l_triggtimeVec.erase(l_triggtimeVec.begin() + missing_indxs.at(rt));
                   l_xoffset.erase(l_xoffset.begin() + missing_indxs.at(rt));
                   l_bytes.erase(l_bytes.begin() + missing_indxs.at(rt)); 
                   l_digits.erase(l_digits.begin() + missing_indxs.at(rt));
                  }       
               nRuns--;
              }
         }
     }

  // nev = 0; // Total number of events on all runs. Do not reuse
  // for (unsigned int qw = 0; qw < nRuns; < qw++) l_nEvents.at(qw) = l_nEvents.at(qw) - 1;
  // for (unsigned int i = 0; i < l_EvCount.size(); i++) nev += l_EvCount.at(i);

  // Binary file processing and event deoding 
  std::vector<unsigned int> trckevntlast;
  unsigned int totalEvnt = 0;
  for (unsigned int a = 0; a < l_runIndx.size(); a++)
      {
       totalEvnt = 0;
       memset(&ifname[0], 0, sizeof(ifname));
       strncpy(ifname, dir, sizeof(ifname));
       strcat(ifname, l_Runs.at(a));
       fd_file = fopen(std::string(ifname).c_str(), "rb");
       if (fd_file == NULL)
          { 
           std::cout << __FUNCTION__ << " ERROR: Could not open file " << std::string(ifname) << ". Skipping run..." << std::endl;
fileprob4:
           l_channels.erase(l_channels.begin() + a); l_ChScope.erase(l_ChScope.begin() + a);
           l_yscale.erase(l_yscale.begin() + a); l_yoffset.erase(l_yoffset.begin() + a);
           l_points.erase(l_points.begin() + a); l_snrate.erase(l_snrate.begin() + a);
           l_EvCount.erase(l_EvCount.begin() + a); l_nScopes.erase(l_nScopes.begin() + a);
           l_nfiles.erase(l_nfiles.begin() + a); l_Runs.erase(l_Runs.begin() + a);
           l_runIndx.erase(l_runIndx.begin() + a);
           if (m_instrument == TestBeamBin2)
              {
               l_dattype.erase(l_dattype.begin() + a); l_triggtimeVec.erase(l_triggtimeVec.begin() + a);
               l_xoffset.erase(l_xoffset.begin() + a); l_bytes.erase(l_bytes.begin() + a);
               l_digits.erase(l_digits.begin() + a);
              }
           nRuns--; a--;
           fclose(fd_file);
           continue;
          }
       printf("Starting loop over %d events on %llu channels, number of samples per channel and event: %d", l_nEvents.at(a), (l_channels.at(a)).size(), (l_points.at(a)).at(0));
       // Fix the start - stop event with respect to events available in the run
       if (LGADBase::GetTrackComb()) 
          {
           trckevntlast.clear();
           l_nEvents.at(a) = l_nEvents.at(a) - 1;
           trckevntlast.resize(1 + m_nplanes.at(a) + m_nDUT.at(a) + m_nRef.at(a), 0);
          }
       CheckStartStopEvnt(&evt1, &evt2, 0, l_nEvents.at(a));
       nel = 0; // Parsed events per run, do not reuse, used to skip events if necessery
       for (unsigned int d = 0; d < (l_EvCount.at(a)).size(); d++)
           {
            if (d > 0) nel += (l_EvCount.at(a)).at(d-1);
            // First test of events, to not even perform size estimation and buffer reservation
            if ((LGADBase::GetTrackComb() && (long int)(nel + (l_EvCount.at(a)).at(d) - 1) < evt1) ||
                (!LGADBase::GetTrackComb() && (long int)(nel + (l_EvCount.at(a)).at(d)) < evt1)) continue;
            if ((LGADBase::GetTrackComb() && (long int)(nel + (l_EvCount.at(a)).at(d) - 1) > evt2) ||
                (!LGADBase::GetTrackComb() && (long int)(nel + (l_EvCount.at(a)).at(d)) > evt2))
               {
jumpnextrun:
                break;
               }
            // Calculate the buffer size for the current event segment
            found = 0; // reuse the integer
            for (unsigned f = 0; f < (l_channels.at(a)).size(); f++)
                {
                 if (m_instrument == TestBeamBin1) found += ((l_points.at(a)).at(f) + 1) * 2 * (l_EvCount.at(a)).at(d);
                 else if (m_instrument == TestBeamBin2) found += (((l_points.at(a)).at(f) * (l_bytes.at(a)).at(f)) + 3 + (l_digits.at(a)).at(f)) * (l_EvCount.at(a)).at(d) - 1;
                }
            // Copy file to buffer, Reading sequence by sequence
            fileBuf = (char*)malloc(found);
            size = fread(fileBuf, found, 1, fd_file);
            if (fileBuf == NULL || size != 1) 
               { 
                std::cout << __FUNCTION__ << " ERROR: Buffer is empty, check data file " << std::string(ifname) << " skipping...." << std::endl;
                if (d == 0) goto fileprob4;
                else break;
               }
            // Start Converting events in buffer
            for (unsigned int i = 0; i < (l_EvCount.at(a)).at(d); i++)
                {
                 nel += i;
                 if (LGADBase::GetTrackComb())
                    { 
                     if ((long int)nel < evt1+1) continue;
                     if ((long int)nel > evt2+1) goto jumpnextrun;
                     m_event = nel-1;
                     LGADBase::ProgressBar((evt2-evt1)*a + nel-1, (evt2 - evt1)* nRuns);
                    }
                 else {
                       if ((long int)nel < evt1) continue;
                       if ((long int)nel > evt2) goto jumpnextrun;
                       m_event = nel;
                       LGADBase::ProgressBar((evt2-evt1)*a + nel, (evt2 - evt1)* nRuns);
                      }
                 fnd = false; // reuse the boolean for skipping events
                 for (unsigned int j = 0; j < (l_channels.at(a)).size(); j++)
                     {  
                      if (m_instrument == TestBeamBin2 && j == 0) m_trigtime = (l_triggtimeVec.at(a)).at(0);
                      nrec = 0; // used for index corrsponddance, do not reuse
                      for (unsigned int prt = 0; prt < m_channels.size(); prt++)
                          {
                           if ((l_channels.at(a)).at(j) == m_channels.at(prt)) { nrec = prt; break; }
                          }
                      if (m_event == evt1)
                         { 
                          if (m_instrument == TestBeamBin1) m_scale.at(nrec) = (l_yscale.at(a)).at(j) / 2.64674e-4;
                          else {
                                if ((l_dattype.at(a)).at(j) < 200) m_scale.at(nrec) = (l_yscale.at(a)).at(j) / 1.5625e-4; // better ADC than Agilent
                                else if ((l_dattype.at(a)).at(j) > 200) m_scale.at(nrec) = (l_yscale.at(a)).at(j) / 4e-2;
                               }
                          m_npoints.at(nrec) = (l_points.at(a)).at(j);
                          m_srate.at(nrec) = (l_snrate.at(a)).at(j);
                          m_scope.at(nrec) = (l_ChScope.at(a)).at(j);
                         }
                      (m_w.at(nrec)).clear();
                      (m_t.at(nrec)).clear();
                      // Fix the number of reserved points to reject unecessety data
                      if (LGADBase::GetTrackComb() && m_oscdel.at((l_ChScope.at(a)).at(j)) != -99 && m_trigclk != -99)
                         {
                          // cacluclate the number of points required for +/- 2 clock cycles
                          if (m_verbose >=2 ) std::cout <<"Channel: " << j << ", Scope: " << (l_ChScope.at(a)).at(j) << " delay points: " 
                                                        << 2*ceil(m_oscdel.at((l_ChScope.at(a)).at(j))*(l_snrate.at(a)).at(j))
                                                        << " trigger clock points: " << 2*ceil(2*m_trigclk*(l_snrate.at(a)).at(j)) << " " 
                                                        << m_oscdel.at((l_ChScope.at(a)).at(j)) << " " << m_trigclk << " " << (l_snrate.at(a)).at(j) << std::endl;
                          if (ceil(m_oscdel.at((l_ChScope.at(a)).at(j))*(l_snrate.at(a)).at(j)) < ceil(2*m_trigclk*(l_snrate.at(a)).at(j)))
                             {
                              m_w.at(nrec).reserve(2*ceil(m_oscdel.at((l_ChScope.at(a)).at(j))*(l_snrate.at(a)).at(j)));
                              m_t.at(nrec).reserve(2*ceil(m_oscdel.at((l_ChScope.at(a)).at(j))*(l_snrate.at(a)).at(j)));
                             }
                          else {
                                m_w.at(nrec).reserve(ceil(4*m_trigclk*(l_snrate.at(a)).at(j)));
                                m_t.at(nrec).reserve(ceil(4*m_trigclk*(l_snrate.at(a)).at(j)));
                               }
                         }
                      else { 
                            m_w.at(nrec).reserve((l_points.at(a)).at(j));
                            m_t.at(nrec).reserve((l_points.at(a)).at(j));
                           }
                      // Oscilloscope data readout is performed channel by channel, reconstruct each envent by pointing the buffer
                      // to the correct position for each chunnel block. Buffer positions are related to buffer type
                      found = 0; // re-use the integer
                      for (unsigned int h = 0; h < j; h++) 
                          {
                           if (m_instrument == TestBeamBin1) found += ((l_points.at(a)).at(h)+1)*2*(l_EvCount.at(a)).at(d);
                           else if (m_instrument == TestBeamBin2) found += (((l_points.at(a)).at(h)*(l_bytes.at(a)).at(h))+3+(l_digits.at(a)).at(h))*(l_EvCount.at(a)).at(d)-1;
                          }
                      // previous loop gets the itterator to skip the channels we already read, we still need to get to the 
                      // correct position within the channel block currently being read for the specific event
                      if (m_instrument == TestBeamBin1) found += i*2*((l_points.at(a)).at(j)+1)-1;
                      if (m_instrument == TestBeamBin2) found += i*(((l_points.at(a)).at(j)*(l_bytes.at(a)).at(j))+3+(l_digits.at(a)).at(j))+2+(l_digits.at(a)).at(j);
                      for (unsigned int w = 0; w < (l_points.at(a)).at(j); w++)
                          {                       
                           if (m_instrument == TestBeamBin2) 
                              {
                               // lets cut off unnecessery points, set to accomodate +/- 2 x 25 nsec FEi4 cycles
                               if (LGADBase::GetTrackComb() && m_oscdel.at((l_ChScope.at(a)).at(j)) != -99)
                                  {
                                   // skip dpoint until you get top the point of interest
                                   if (m_verbose >= 2) std::cout << w << " : " << m_oscdel.at((l_ChScope.at(a)).at(j))*(l_snrate.at(a)).at(j) << " " << m_w.at(nrec).capacity()/2 
                                                                 << " : " << ceil(m_oscdel.at((l_ChScope.at(a)).at(j))*(l_snrate.at(a)).at(j)-(m_w.at(nrec).capacity()/2)) << std::endl;
                                   if (w < ceil(m_oscdel.at((l_ChScope.at(a)).at(j))*(l_snrate.at(a)).at(j)-(m_w.at(nrec).capacity())/2)) continue;
                                   // Don't process points after the point of interest, you need to adjust the position here to accoutn for the skipped poitns
                                   if (w >= ceil(m_oscdel.at((l_ChScope.at(a)).at(j))*(l_snrate.at(a)).at(j)+((m_w.at(nrec)).capacity())/2)) break;
                                  }
                               // 1 byte, MSB, signed int RI (SHOULD ADD +128?)
                               if ((l_dattype.at(a)).at(j) == 111 || (l_dattype.at(a)).at(j) == 121) m_w.at(nrec).push_back((double)(((int8_t)fileBuf[found+w]))*(l_yscale.at(a)).at(j) + (l_yoffset.at(a)).at(j));
                               // 1 byte, MSB, positive int RP
                               else if ((l_dattype.at(a)).at(j) == 112 || (l_dattype.at(a)).at(j) == 122) m_w.at(nrec).push_back((double)((uint8_t)fileBuf[found+w])*(l_yscale.at(a)).at(j) + (l_yoffset.at(a)).at(j));
                               // 2 byte, MSB, signed int RI (SHOULD ADD +32,768?)
                               else if ((l_dattype.at(a)).at(j) == 211) m_w.at(nrec).push_back((double)(((int8_t)fileBuf[found+2*w]<<8) | ((uint8_t)fileBuf[found+2*w+1]))*(l_yscale.at(a)).at(j) + (l_yoffset.at(a)).at(j));
                               // 2 byte, MSB, positive int RP 
                               else if ((l_dattype.at(a)).at(j) == 212) m_w.at(nrec).push_back((double)(((uint8_t)fileBuf[found+2*w]<<8) | ((uint8_t)fileBuf[found+2*w+1]))*(l_yscale.at(a)).at(j)+(l_yoffset.at(a)).at(j));
                               // 2 byte, LSB, signed int RI (SHOULD ADD +32,768?)
                               else if ((l_dattype.at(a)).at(j) == 221) m_w.at(nrec).push_back((double)((uint8_t)fileBuf[found+2*w] | ((int8_t)fileBuf[found+2*w+1]<<8))*(l_yscale.at(a)).at(j) + (l_yoffset.at(a)).at(j));
                               // 2 byte, LSB, positive int RP 
                               else if ((l_dattype.at(a)).at(j) == 222) m_w.at(nrec).push_back((double)((uint8_t)fileBuf[found+2*w] | ((uint8_t)fileBuf[found+2*w+1]<<8))*(l_yscale.at(a)).at(j) + (l_yoffset.at(a)).at(j));
                               m_t.at(nrec).push_back((double)(w+(l_xoffset.at(a)).at(j))/(l_snrate.at(a)).at(j));
                              }
                           else { // Assumes LSB representation (SMALL ENDIANESS)
                                 m_w.at(nrec).push_back((short int)(fileBuf[found+2*w]+fileBuf[found+2*w+1])*256*(l_yscale.at(a)).at(j) + (l_yoffset.at(a)).at(j));
                                 m_w.at(nrec).push_back((double)w/(l_snrate.at(a)).at(j));
                                }
                          }
                      // Consistency check that oscilloscope data changed from channel to channel
                      if (j > 0) 
                         { 
                          scope = 0;
                          for (unsigned int prt = 0; prt < m_channels.size(); prt++)
                              {
                               if ((l_channels.at(a)).at(j-1) == m_channels.at(prt)) { scope = prt; break; }
                              }
                          if (LGADBase::IsVecEqual(m_w.at(nrec), m_w.at(scope)))
                             {
                              std::cout << __FUNCTION__ << "          ERROR: Removing event " << nel+1 << " of run " << std::string(ifname) 
                                        << ", detected no channel change at oscilloscope readout between channels "
                                        << (l_channels.at(a)).at(j) << " and " << (l_channels.at(a)).at(j - 1) << "!!" << std::endl;
                              fnd = true;
                              break;
                             }
                         }
                     }
                 if (!fnd)
                    {
                     if (LGADBase::GetTrackComb())
                        {
                         if (LGADBase::GetTrackPackage() == "Cory") LGADBase::CombineTrackCory(m_event, trckevntlast, a);
                         else if (LGADBase::GetTrackPackage() == "EUTelescope") LGADBase::CombineTrackEUTel(m_event, trckevntlast, a);
                         else if (LGADBase::GetTrackPackage() == "TBA") LGADBase::CombineTrackTBA(m_event, trckevntlast, a);
                         else {
                               std::cout << __FUNCTION__ << "          ERROR: Unsopported tracking package " << LGADBase::GetTrackPackage() << std::endl;
                               return false;
                              }
                        }
                     // Before filling the tree have to populate any not included channel
                     for (unsigned int prt = 0; prt < m_channels.size(); prt++)
                         {
                          fnd = false;
                          for (unsigned int hj = 0; hj < (l_channels.at(a)).size(); hj++)
                              {
                               if ((l_channels.at(a)).at(hj) == m_channels.at(prt)) { fnd = true; break; }
                              }
                          if (!fnd)
                             {
                              (m_w.at(prt)).clear(); (m_t.at(prt)).clear();
                              (m_w.at(prt)).resize((l_points.at(a)).at(0), -1); 
                              (m_w.at(prt)).resize((l_points.at(a)).at(0), -1);
                              m_scope.at(prt) = -1;
                              m_scale.at(prt) = -1;
                              m_npoints.at(prt) = -1;
                              m_srate.at(prt) = -1;
                              m_trigtime = -1;
                             }
                         }
                     totalEvnt++;
                     m_trees.at(0)->Fill(); // fill every event
                    }
                }
            free(fileBuf);
           }
       std::cout << "Wrote " << totalEvnt << " events from run " << std::string(ifname) << std::endl;
       if (m_instrument == TestBeamBin2) LGADBase::CalcTrigFr(l_triggtimeVec.at(a), m_trigDt, m_trigFr, totalEvnt);
       fclose(fd_file);
      }

   LGADBase::CloseTrackFiles();
   TListIter l_next_object(gDirectory->GetList());
   l_next_object.Reset();  
   TObject* l_obj;
   while ((l_obj=l_next_object()))
         {
          l_obj->Write("", 6);
          delete l_obj;
         }
   m_ofile->Close();

   return true;
}
// --------------------------------------------------------------------------------------------------------------
bool LGADBase::CreateOutputFile(const char* dir, const char* ofname, std::vector<unsigned int> nchan)
{
    // Define Output and reserve branches
#ifdef _WIN32
	char ifname[MAX_PATH];
#else
	char ifname[PATH_MAX];
#endif
    strncpy(ifname, dir, sizeof(ifname));

    std::string print = ifname;
    print.erase(print.size() - 1);
    int a = LGADBase::DirExists(print.c_str());
    if (a == 0)
       {
        std::cout << __FUNCTION__ << " WARNING: Selected output directory does not exist" << std::endl;
        a = LGADBase::RecursMkDir(print.c_str());
        if (a != 1)
           {
            std::cout << __FUNCTION__ << " ERROR: Failed to create " << print.c_str()<< "!" << std::endl;
            return false;
           }
        else std::cout << __FUNCTION__ << " WARNING: Succesfully created directory: " << print.c_str() << std::endl;
       }
    else if (a == -1)
            {
             std::cout << __FUNCTION__ << " ERROR: Selected path " << print.c_str()<< " is not a directory" << std::endl;
             return false;
            }

    strcat(ifname, ofname);
    m_ofile = new TFile(Form("%s.root", ifname), "recreate");
    m_ofile->SetCompressionLevel(6);
    if (!GetTrackComb()) std::cout << "    Creating output ROOT file: " << m_ofile->GetName() << std::endl;
    m_trees.push_back(new TTree("wfm", "Timing part of charge collection measurements"));
    m_trees.at(0)->SetAutoFlush(-300000000); // Set autoflush to 300 MB
    LGADBase::SetVectorSize(nchan.size());

    for (unsigned int ich = 0; ich < nchan.size(); ich++)
        {
         if (m_instrument == Sampic)
            {
             m_trees.at(0)->Branch(Form("physt%02u", nchan.at(ich)), &(m_physt.at(ich)), "m_physt.at(ich)/D");
             m_trees.at(0)->Branch(Form("ordrt%02u", nchan.at(ich)), &(m_ordrt.at(ich)), "m_ordrt.at(ich)/D");
            }
         m_trees.at(0)->Branch(Form("t%02u", nchan.at(ich)), &(m_t.at(ich)));
         m_trees.at(0)->Branch(Form("w%02u", nchan.at(ich)), &(m_w.at(ich)));
         if (m_instrument == TestBeamBin1 || m_instrument == TestBeamBin2 || m_instrument == LabTXT || m_instrument == TektronixScope || m_instrument == LeCroyWRBin)
            {
             m_trees.at(0)->Branch(Form("vScale%02u", nchan.at(ich)), &(m_scale.at(ich)), "m_scale.at(ich)/F");
             m_trees.at(0)->Branch(Form("nPoints%02u", nchan.at(ich)), &(m_npoints.at(ich)), "m_npoints.at(ich)/i");
             m_trees.at(0)->Branch(Form("SnRate%02u", nchan.at(ich)), &(m_srate.at(ich)), "m_srate.at(ich)/L");
             if (m_instrument == TektronixScope) m_trees.at(0)->Branch(Form("triggTime%02u", nchan.at(ich)), &(m_triggTime.at(ich)), "m_triggTime.at(ich)/D");
             if (m_instrument == TestBeamBin1 || m_instrument == TestBeamBin2 || m_instrument == LeCroyWRBin) m_trees.at(0)->Branch(Form("Scope%02u", nchan.at(ich)), &(m_scope.at(ich)), "m_scope.at(ich)/i");
            }
        }
    if (m_instrument == LabTXT || m_instrument == TestBeamBin2) m_trees.at(0)->Branch("trigtime", &m_trigtime, "m_trigtime/D");
    m_trees.at(0)->Branch("EvnNo", &m_event, "m_event/i");
    if (m_instrument == Sampic)
       {
        m_trees.at(0)->Branch("nPoints", &(m_npoints.at(0)), "m_npoints.at(0)/i");
        m_trees.at(0)->Branch("SnRate", &(m_srate.at(0)), "m_srate.at(0)/L");
       }

    if (m_instrument != TestBeamBin1 || m_instrument != LeCroyWRBin || m_instrument != TektronixScope)
       {
        m_trigDt = new TH1F("trigDt", "DT time", 6000, 0.0, 60.0);
        m_trigFr = new TH1F("trigFr", "Trigger Frequency", 50000, 0.0, 5000);
       }

    // We are printing this after the addition of the tracking branch planes
    if (!GetTrackComb())  
       {
        PrintBranches(m_trees.at(0), "Waveform");
        PrintObjects(m_ofile, "wfm");
       }

   return true;
}
// --------------------------------------------------------------------------------------------------------------
void LGADBase::CalcTrigFr(std::vector<double> EvTrigTime, TH1F* TrigPer, TH1F* TrigFrq, unsigned int entriesNo)
{
  float dT1 = 0;
  std::vector<double> DeltaTrigTime;
  for (unsigned int hp = 0; hp < EvTrigTime.size() - 1; hp++)
      {
       if (hp > 0 && EvTrigTime.at(hp) == 0) continue;
       dT1 = fabs(EvTrigTime.at(hp + 1) - EvTrigTime.at(hp));
       if (dT1 > 30) continue;
       else DeltaTrigTime.push_back(dT1);
       TrigPer->Fill(dT1);
      }

  float tbin = 0;
  int count = 0;
  for (unsigned int ka = 0; ka < EvTrigTime.size() - 1; ka++)
      {
       count = 1;
       dT1 = -1;
       for (unsigned int ma = ka + 1; ma < EvTrigTime.size(); ma++)
           {
            if ((EvTrigTime.at(ma) - EvTrigTime.at(ka)) == 0)
               {
                if (ma < (EvTrigTime.size()-1) && (EvTrigTime.at(ma) != 0))
                   {
                    count++;
                    continue;
                   }
                else {
                      if (tbin == 0) 
                         { 
                          unsigned int scope = 0;
                          tbin = EvTrigTime.at(ma);
                          tbin = fabs(tbin);
                          tbin = tbin - (int)tbin;
                          while (tbin != 0)
                               {
                                tbin = tbin * 10;
                                scope++;
                                tbin = tbin - (int)tbin;
                               }
                          tbin = 1/pow(10, scope);
                         } // THIS HAS TO BE CHANGED
                      dT1 = (float)count/tbin;
                      TrigFrq->Fill(dT1);
                      ka = ma;
                     }
               }
            else {
                  if (EvTrigTime.at(ma) != 0)
                     {
                      tbin = EvTrigTime.at(ma) - EvTrigTime.at(ka);
                      dT1 = (float)count/tbin;
                      TrigFrq->Fill(dT1);
                      ka = ma - 1;
                     }
                  else ka = ma;
                  break;
                 }
           }
      }

  if (EvTrigTime.back() != 0) dT1 = EvTrigTime.back() - EvTrigTime.at(0);
  else dT1 = EvTrigTime.at(EvTrigTime.size()-2) - EvTrigTime.at(0);
  std::cout << "Event rate = " << entriesNo << " / " << dT1 << " = " << (double)(entriesNo / dT1) << " Hz" << std::endl;
  dT1 = LGADBase::Mean(&DeltaTrigTime);
  std::cout << "Mean Trigger: " << dT1;
  if (dT1 != 0 && dT1 != -1) std::cout << ", mean rate: " << 1 / dT1 << std::endl;
  else std::cout << ", no mean rate calculation possible!" << std::endl;
  dT1 = LGADBase::CalcMeadian(&DeltaTrigTime);
  std::cout << "Median Trigger: " << dT1;
  if (dT1 != 0 && dT1 != -1) std::cout << ", median rate: " << 1 / dT1 << std::endl;
  else std::cout << ", no median rate calculation possible!" << std::endl;
}
// --------------------------------------------------------------------------------------------------------------
int LGADBase::RunCount(std::vector<TString> &names, std::vector<unsigned int> &nfiles, std::vector<TString> &Runs, std::vector<TString>& ext,
                       std::vector<unsigned int> &runIndx, const char* inext, std::vector<std::vector<unsigned int>>* nchannels)
{
  int nRuns = -1;
  std::size_t found = 0;
  unsigned int ch = 0;
  std::string fname;
  std::string fname_old;
  bool fnd = false;
  // Determine the number of runs contained in the list of files to analyse
  for (unsigned int hwr = 0; hwr < names.size(); hwr++)
      {
       if (strcmp(inext, "") == 0)
          {
           if (!(names.at(hwr).Contains("."))) { names.erase(names.begin()+hwr); hwr--; continue; }
           else found = (names.at(hwr)).Last('.');
          }
       else found = (names.at(hwr)).Length();
       if (hwr == 0) 
          { 
           Runs.push_back(names.at(0)(0, found));
           if (m_instrument == TektronixScope)
              {
               nchannels->push_back(std::vector<unsigned int>());
               found = ((std::string)Runs.at(0)).rfind("_Ch");
               fname_old = (std::string)(Runs.at(0)(0, (Runs.at(0)).Length() - found));
               ch = stoi((std::string)(Runs.at(0)(found + 3, names.at(hwr).Length() - (found + 3))));
               Runs.at(0) = (TString)fname_old;
               (nchannels->at(0)).push_back(ch);
              }
           runIndx.push_back(0);
           nfiles.push_back(1);
           if (strcmp(inext, "") == 0) ext.push_back(names.at(0)(found + 1, names.at(0).Length()));
           else {
                 ext.push_back(inext);
                 if (m_instrument == TestBeamBin || m_instrument == LeCroyWRBin || 
                     m_instrument == TestBeamBin1 || m_instrument == TestBeamBin2) ext.push_back("txt");
                }
          }
       else {
             // Extention check that has anywway to be done
             if (strcmp(inext, "") == 0)
                {
                 fnd = false;
                 for (unsigned int a = 0; a < ext.size(); a++)
                     {
                      if (ext.at(a) == names.at(hwr)(found + 1, names.at(0).Length())) 
                         { 
                          fnd = true; 
                          break; 
                         }
                     }
                 if (!fnd) ext.push_back(names.at(hwr)(found + 1, names.at(0).Length()));
                }
             fname = (names.at(hwr)(0, found));
             if (m_instrument == TektronixScope)
                {
                 found = fname.rfind("_Ch");
                 fname_old = fname.substr(0, fname.length() - found);
                 ch = stoi((std::string)(fname.substr(found + 3, fname.length() - (found + 3))));
                 fname = fname_old;
                }
             if (fname != Runs.back())
                {
                 Runs.push_back((TString)fname);
                 runIndx.push_back(Runs.size());
                 nfiles.push_back(1);
                 if (m_instrument == TektronixScope) nchannels->push_back(std::vector<unsigned int>());
                }
             else nfiles.back()++;
             if (m_instrument == TektronixScope) (nchannels->back()).push_back(ch);
            }
      }
  if (Runs.size() > 0) nRuns = Runs.size();
  return nRuns;
}
// --------------------------------------------------------------------------------------------------------------
int LGADBase::ArrangeCahnnels(std::vector<unsigned int> &channels, std::vector<unsigned int> &ChScope)
{
   // Find the last channel of each scope
   std::vector<unsigned int> lclch;
   lclch.push_back(0);
   unsigned int we = 0;
   for (unsigned int d = 0; d < channels.size(); d++)
       {
        if (channels.at(d) > lclch.at(we)) lclch.at(we) = channels.at(d);
        else { lclch.push_back(channels.at(d)); we = we + 1; }
       }
   // Determine the maximum amount of channels that each socpe can have (multitudes of 4)
   std::vector<unsigned int> maxch;
   maxch.resize(lclch.size(), 0);
   for (unsigned int m = 0; m < lclch.size(); m++)
       {
        int remainder = lclch.at(m) % 4;
        if (remainder == 0) return maxch.at(m) = lclch.at(m);
        else maxch.at(m) = lclch.at(m) + 4 - remainder;
       }
   unsigned int rw = 0;   // index of the max channel and scope channel vectors
   unsigned int chft = 0; // commulative number of channels up to the current scope
   ChScope.resize(channels.size(), 0);
   for (unsigned int d = 0; d < channels.size(); d++)
       {
        channels.at(d) += chft;
        ChScope.at(d) = rw + 1;
        if (channels.at(d) == lclch.at(rw)) 
           {  
            rw++;
            chft += maxch.at(rw);
           }
       }
   return lclch.size();
}
// --------------------------------------------------------------------------------------------------------------
bool LGADBase::FixBaseQuant(std::vector<std::vector<unsigned int>> &filevect, std::vector<unsigned int> &basevect, unsigned int &basecounter)
{
  bool fnd = false;
  unsigned int nRuns = filevect.size();
  // Populate the channel vector from all runs
  basevect.clear();
  for (unsigned int a = 0; a < nRuns; a++)
      {
       for (unsigned int g = 0; g < (filevect.at(a)).size(); g++)
           {
            fnd = false;
            for (unsigned int b = 0; b < basevect.size(); b++)
               {
                if ((filevect.at(a)).at(g) == basevect.at(b)) { fnd = true; break; }
               }
            if (!fnd) basevect.push_back((filevect.at(a)).at(g));
           }
      }
  basecounter = basevect.size();
  if (basecounter == 0)
     {
      std::cout << __FUNCTION__ << " ERROR: Failed to find correctly determine the numbeer of unique elements!" << std::endl;
      return false;
     }
  return true;
}
// ---------------------------------------------------------------------------------------------------------------
bool LGADBase::CheckStartStopEvnt(unsigned int *evt1, unsigned int *evt2, unsigned int min_evnt, unsigned int max_evnt)
{
  if (*evt1 == 0) evt1 = &min_evnt;
  if (*evt2 == 0) evt2 = &max_evnt;
  if (*evt1 == min_evnt || *evt1 < min_evnt)
     {
      if ((max_evnt - min_evnt) <= (*evt2 - *evt1)) evt2 = &min_evnt + (*evt2 - *evt1);
      else evt2 = &max_evnt;
      if (*evt1 < min_evnt) evt1 = &min_evnt;
      std::cout << __FUNCTION__ << "          Adjusting start/stop event to available from files to " << *evt1 << " and " << evt2 << std::endl;
     }
  else { 
        if ((max_evnt - *evt1) > (*evt2 - *evt1)) evt2 = &max_evnt;
        std::cout << __FUNCTION__ << "          Adjusting stop event to available number of events from files to " << *evt2 << std::endl;
       }
  return true;
}
// --------------------------------------------------------------------------------------------------------------