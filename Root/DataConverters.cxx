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
    TString sLast = DataDir[DataDir.Length() - 1];
    if (!sLast.Contains("/")) DataDir = DataDir + "/";
    m_datadir = DataDir; 
    if (DataName.Length() > 4 && ext =="")
       {
        if (DataName.Contains(".txt")) { m_dataname = DataName(0, DataName.Length() - 4); m_ext = ".txt"; }
        else if (DataName.Contains(".dat")) { m_dataname = DataName(0, DataName.Length() - 4); m_ext = ".dat"; }
        else if (DataName.Contains(".root")) { m_dataname = DataName(0, DataName.Length() - 5); m_ext = ".root"; }
        else { m_dataname = DataName; m_ext = ext; }
       }
    else if (DataName != "" && ext != "") { m_dataname = DataName; m_ext = ext; }
    else { m_dataname = ""; m_ext = ""; }
}
//---------------------------------------------------------------------------------------------------------------
void LGADBase::SetTrackInDataNames(TString DataDir, TString DataName)
{
    TString sLast = DataDir[DataDir.Length() - 1];
    if (!sLast.Contains("/")) DataDir = DataDir + "/";
    m_Trackdatadir = DataDir;
    TString sLast2 = DataName[DataDir.Length() - 5];
    if (!sLast2.Contains(".root")) m_Trackdatadir = DataName + ".root";
    else m_Trackdatadir = DataName;
}
// --------------------------------------------------------------------------------------------------------------
void LGADBase::SetOutDataNames(TString DataDir, TString DataName)
{
    TString sLast = DataDir[DataDir.Length() - 1];
    if (!sLast.Contains("/")) DataDir = DataDir + "/";
    m_ofname = DataName;
    m_ofdir = DataDir;
}
// --------------------------------------------------------------------------------------------------------------
bool LGADBase::ConvertData()
{
    if (m_datadir.IsNull()) { std::cout << __FUNCTION__ << " ERROR: Input data directory not set! Abording..." << std::endl; return false; }

    if (m_dataname.IsNull() && (m_instrument == LabTXT || m_instrument == TektronixScope))
       { 
        std::cout << "--> Warning: Input data file name not set! Trying with default for instrument: --Trace-- " << std::endl;
        m_dataname = "--Trace--";
       }
    else if (m_dataname.IsNull() && m_instrument == Sampic)
            {
             m_dataname = m_datadir;
             TString sLast1 = m_dataname[m_dataname.Length() - 1];
             if (sLast1.Contains("/")) m_dataname = m_dataname(0, m_dataname.Length() - 1);
             m_dataname = m_dataname(m_dataname.Last('/') + 1, m_dataname.Length() - m_dataname.Last('/') - 1);
             std::cout << "--> Warning: Input data file name not set! Trying with same name as containing folder: " << m_dataname << std::endl;           
            }
    else if (m_dataname.IsNull()) { std::cout << __FUNCTION__ << " ERROR: Input data file name not set! Abording..." << std::endl; return false; }

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
        else m_ofname = m_dataname;
        std::cout << m_ofname << std::endl;
       }
    if (m_ofdir.IsNull())
       {
        std::cout << "--> Warning: Output Ntuple location not set! Using input data diretory" << std::endl;
        m_ofdir = m_datadir;
       }
    if (m_TrackComb && (m_Trackdatadir.IsNull() || m_Trackdataname.IsNull()))
       { 
        std::cout << "--> Warning: Tracking input file not defined, switching to non-test beam mode" << std::endl;
        m_TrackComb = false;
       }

    m_nchan = 0;
    m_channels.clear();
    m_convert = false;
    if (m_instrument == LeCroyWRBin || m_instrument == TektronixBinary || m_instrument == Sampic || m_instrument == InfiniiumScope)
       {
        if (m_ext.IsNull()) 
           {
            if (m_instrument == Sampic) m_convert = WriteSampic(m_datadir.Data(), m_dataname.Data(), "dat", m_evnt1, m_evnt2);
            else if (m_instrument == LeCroyWRBin) m_convert = WriteLecroyBinary(m_datadir.Data(), m_dataname.Data(), "dat", m_evnt1, m_evnt2);
            else if (m_instrument == InfiniiumScope) m_convert = WriteAgilentBinary(m_datadir.Data(), m_dataname.Data(), "dat", m_evnt1, m_evnt2);
            else (m_instrument == TektronixBinary) m_convert = WriteTekBinary(m_datadir.Data(), m_dataname.Data(), "dat", m_evnt1, m_evnt2);
           }
        else {
              if (m_instrument == LeCroyWRBin) m_convert = WriteLecroyBinary(m_datadir.Data(), m_dataname.Data(), m_ext.Data(), m_evnt1, m_evnt2);
              else if (m_instrument == InfiniiumScope) m_convert = WriteAgilentBinary(m_datadir.Data(), m_dataname.Data(), m_ext.Data(), m_evnt1, m_evnt2);
              else if (m_instrument == Sampic) m_convert = WriteSampic(m_datadir.Data(), m_dataname.Data(), m_ext.Data(), m_evnt1, m_evnt2);
              else (m_instrument == TektronixBinary) m_convert = WriteTekBinary(m_datadir.Data(), m_dataname.Data(), m_ext.Data(), m_evnt1, m_evnt2);
              }
       }
    else if (m_instrument == LabTXT)
            {
             if (m_ext.IsNull()) m_convert = WriteLabTXT(m_datadir.Data(), m_dataname.Data(), "txt", m_evnt1, m_evnt2);
             else m_convert = WriteLabTXT(m_datadir.Data(), m_dataname.Data(), m_ext.Data(), m_evnt1, m_evnt2);
            }
    else if (m_instrument == TektronixScope)
            {
             if (m_ext.IsNull()) m_convert = WriteTectronixTXT(m_datadir.Data(), m_dataname.Data(), "txt");
             else m_convert = WriteTectronixTXT(m_datadir.Data(), m_dataname.Data(), m_ext.Data());
            }

  // this has to go into each of the testbeam supported instruments (WriteAgilentBinary, WriteLecroyBinary)
  if (m_TrackComb) CombineTrack(m_Trackdatadir.Data(), m_Trackdataname.Data());

  return true;
}
// --------------------------------------------------------------------------------------------------------------
// Functon to read Sampic data. Because of the way sampic outputs events in a single consecutive
// .dat file, start and stop event numbers are of no conseuence here. Functon defines automatically
// the number of active channels and their respective position. Only accepts evetns that are within
// one buffer window. Time corrections for bin shifting is applied tking into account the ordered cell0
// time. All information are included in the output root file.
bool LGADBase::WriteSampic(const char* dir, const char* name, const char* ext, int evt1, int evt2)
{
  std::ifstream file;
  TBits *b = new TBits;
  char ifname[2048];
  std::string line;
  std::istringstream iss(std::ios_base::app | std::ios_base::in | std::ios_base::out);
  double tri = -99.;
  int TotHits = 0.0;
  std::vector<double> voltage;
  std::vector<double> l_triggtimeVec;
  double phystime;
  double ordrtime;
  unsigned int npoints;
  int channel;
  std::vector<std::vector<double> > l_voltage;
  std::vector<long double> l_phystime;
  std::vector<long double> l_ordrtime;
  std::vector<long int> l_nEvntN;
  std::vector<unsigned int> l_channel;
  unsigned long int ievent = 0;
  unsigned long int tevent = 0;
  Long64_t srate = -1;

  unsigned int nfiles = CountFiles(dir, ext);
  if (nfiles == 0)
     {
      std::cout << __FUNCTION__ << " ERROR: No data files found in " << dir << ", aborting ..." << std::endl;
      return false;
     }

  // Find the last event in the recorded series
  // Open the last file in the folder
  if (nfiles == 1) sprintf(ifname, "%s%s.%s", dir, name, ext);  // build the filename for the first event in the first channel
  else sprintf(ifname, "%s%s.%s_%04d", dir, name, ext, nfiles - 1);
  file.open(ifname);
  if (file.is_open())  // Start reading from the back until you get the 3rd line from last
     {
      file.seekg(-2, file.end);
      unsigned int rdl = 0;
#ifdef _WIN32
      unsigned int ps = 4;
#else
      unsigned int ps = 3;
#endif
      while (rdl < ps)
            {
             char cdfcv;
             file.get(cdfcv);
             if ((int)file.tellg() <= 1)
                {
                 file.seekg(0);
                 rdl = 4;
                }
             else if (cdfcv == '\n') 
                     {
                      file.seekg(-3, ios_base::cur); 
                      rdl++;
                     }
             else file.seekg(-2, ios_base::cur);
            }
      file.seekg(3, ios_base::cur);
      getline(file, line);
      line = reduce(line, "", "=");
      sscanf(line.c_str(), "%*s %u %*s %*lf", &TotHits);
      file.close();
     }

  l_triggtimeVec.clear();
  for (unsigned int ifile = 0; true; ifile++)
      {
       memset(ifname, '0', sizeof(ifname));
       if (ifile == 0)
          { 
           const char* ifname_format = "%s%s.%s";
           sprintf(ifname, ifname_format, dir, name, ext);  // build the filename for the first event in the first channel
          }
       else {
             const char* ifname_format = "%s%s.%s_%04d";
             sprintf(ifname, ifname_format, dir, name, ext, ifile);
            }
       file.open(ifname);
       if (!file)
          {
           if (ifile == 0) 
              {
               std::cout << __FUNCTION__ << " ERROR: Unable to find input file: " << dir << name << "." << ext << std::endl;
               return false;
              }
           else if (ifile < nfiles) continue;
           else break;
          }
       std::cout << "--> Openning file: " << ifname << std::endl;
       npoints = 0;

       // If first file then read all info for rate and number of active channels
       if (ifile == 0)
          {
           unsigned int sysch = 0;
           line.clear();
           getline(file, line); // skip the first line, may include information later
           // Find the total number of channels in the system
           line.clear();
           getline(file, line);
           line = reduce(line, "", "===");
           line = LGADBase::trim(line, " ");
           sscanf(line.c_str(), "%*s %*s %*f %*s %*s %*s %*s %*s %d", &sysch);
           // Find the number of active channels by converting the hex register value to a decimal
           line.clear();
           getline(file, line);
           line = reduce(line, "", "===");
           line = LGADBase::trim(line, " ");
           sscanf(line.c_str(), "%*s %*s %llu %*s %*s %*s %*s %x", &srate, &m_nchan);
           b->Set(sysch, &m_nchan);
           m_nchan = 0;
           for (unsigned int i = 0; i < sysch; i++) if (b->TestBitNumber(i)) { m_channels.push_back(i + 1); m_nchan++; }
           for (unsigned int i = 0; i < 4; i++) {line.clear(); getline(file, line);} // skip the following 4 lines, may include information later
           // create the ntuple
           if (!CreateOutputFile(m_ofdir, m_ofname, m_channels))
              {
               std::cout << __FUNCTION__ << " ERROR: Output Ntuple creation failed: " << m_ofname << std::endl;
               return false;
              }
          }
       else for (unsigned int i = 0; i < 7; i++) {line.clear(); getline(file, line);} // skip the first 7 lines, may include information later

       // Read out time and voltage of all events in file
       l_voltage.clear();
       l_phystime.clear();
       l_ordrtime.clear();
       l_channel.clear();
       l_nEvntN.clear();
       while (file.good())
             {
              voltage.clear();
              phystime = 0.0;
              ordrtime = 0.0;
              channel = -1;
              line.clear();
              getline(file, line);
              sscanf(line.c_str(), "%*s %lu", &ievent);
              // Find out if start - stop events correspond to expectations (first event of first file) and fix it
              if (ifile == 0 && l_nEvntN.size() == 0)
                 {
                  // Start event doeskys not correspond to expectations, find first event
                  if (evt1 != 0 && evt1 < (long int)ievent)
                     {
                      evt1 = ievent;
                      std::cout << __FUNCTION__ << " WARNING: Given start event number " << evt1 << " in disagreement with first event in file: " << ievent << ", attempting to recover... " << std::endl;
                      if (evt2 != 0 && (TotHits - (long int)ievent) <= (evt2 - evt1)) evt2 = (long int)ievent + (evt2 - evt1);
                      else evt2 = TotHits;
                      std::cout << __FUNCTION__ << "          Adjusting start/stop event to available files from " << evt1 << " to " << evt2 << std::endl;
                     }
                  else {
                        evt1 = ievent;
                        if (evt2 != 0 && evt2 > TotHits) 
                           {
                            std::cout << __FUNCTION__ << " WARNING: Given stop event number " << evt2 << " in disagreement with total events in file: " << TotHits << ", attempting to recover... " << std::endl;
                            evt2 = TotHits;
                            std::cout << __FUNCTION__ << "          Adjusting stop event to available files " << evt2 << std::endl;
                           }
                        else if (evt2 == 0) evt2 = TotHits;
                       }
                 }
              if ((long int)ievent >= evt2) break;
              LGADBase::ProgressBar(floor((float)(ievent-evt1)/(float)(2*m_nchan)), ceil((float)evt2/(float)m_nchan));
              line.clear();
              getline(file, line);
              // if (l_voltage.size() == 0) sscanf(line.c_str(), "%*s %u %*s %lf %*s %*u %*s %*lf %*s %lf %*s %u", &channel, &phystime, &ordrtime, &npoints);
              // else sscanf(line.c_str(), "%*s %u %*s %lf %*s %*u %*s %*lf %*s %lf %*s %*u", &channel, &phystime, &ordrtime);
              if (l_voltage.size() == 0) sscanf(line.c_str(), "%*s %u %*s %lf %*s %*d %*s %*lf %*s %lf %*s %*lf %*s %*lf %*s %*lf %*s %u", &channel, &phystime, &ordrtime, &npoints);
              else sscanf(line.c_str(), "%*s %u %*s %lf %*s %*d %*s %*lf %*s %lf %*s %*lf %*s %*lf %*s %*lf %*s %*u", &channel, &phystime, &ordrtime);
              if (channel == -1) continue;
              voltage.reserve(npoints);
              line.clear();
              getline(file, line);
              line = reduce(line, "", "DataSamples");
              line = trim(line, " ");
              iss.str(std::string());
              iss.clear();
              iss.str(line);
              for (unsigned int i = 0; i < npoints; i++)
                  {
                   iss >> tri; 
                   voltage.push_back(tri);
                  }
              l_voltage.push_back(voltage);
              l_phystime.push_back(phystime);
              l_ordrtime.push_back(ordrtime);
              l_channel.push_back(channel+1);
              l_nEvntN.push_back(ievent);
             }
       file.close();

       // Trying to find conicidences
       ievent = 0;
       unsigned int n = 0;
       std::vector<unsigned int> ChNo;
       bool done = false;
       unsigned int temp = 0;
       for (unsigned int i = 0; i < l_voltage.size(); i++)
           {
            temp = i;
            ChNo.clear();
            ChNo.push_back(l_channel.at(i));
            if (m_channels.size() > 1)
               {
                n = 0;
                // Loop though all events in the file, find if there is any in coincidence with the one we have and populate vectors
                for (unsigned int j = i + 1; j < l_voltage.size(); j++)
                    {
                     if (ChNo.at(0) == l_channel.at(j)) continue; // Ignore events in the same chanel as the one we have
                     else if (fabs(l_phystime.at(i) - l_phystime.at(j)) <= npoints*(1e3/srate))
                             {
                              int wch = -1;
                              for (unsigned int h = 0; h < m_channels.size(); h++) if (m_channels.at(h) == l_channel.at(j)) wch = h;
                              if (wch == -1)
                                 {
                                  std::cout << __FUNCTION__ << " ERROR: Cannot determine which channel this event corrsponds to" << std::endl;
                                  return false;
                                 }
                              m_physt[wch] = l_phystime.at(j);
                              m_ordrt[wch] = l_ordrtime.at(j);
                              m_w[wch] = l_voltage.at(j);
                              m_t[wch].clear();
                              m_t[wch].reserve(npoints);
                              // This is channel jitter realignement, it assumes time resolutions < sampling interval
                              int bin = nearbyint((l_ordrtime.at(i) - l_ordrtime.at(j))*srate*0.001);
                              for (unsigned int k = 0; k < npoints; k++) m_t[wch].push_back((1 / (double)srate)*((int)k - bin)*(1e-6));
                              ChNo.push_back(l_channel.at(j));
                              n++;
                              i = j; // Aassumes time-ordered aquisition
                             }
                     else continue;
                    }
               }
            // We found events in coincidednce (at least one)
            if ((n > 0 && m_channels.size() > 1) || (n == 0 && m_channels.size() == 1))
               {
                // If you have multiple active channels, check that you have found events in conincidence in all channels, if not
                // fill the ones that are not found with -1.
                for (unsigned int f = 0; f < m_nchan; f++) 
                    {
                     bool found = false;
                     for (unsigned int g = 0; g < n+1; g++) if (m_channels.at(f) == ChNo.at(g)) found = true;
                     if (!found)
                        {
                         m_t[f].clear();
                         m_physt[f] = 0;
                         m_ordrt[f] = 0;
                         m_w[f].clear();
                         m_t[f].reserve(npoints);
                         m_w[f].reserve(npoints);
                         for (unsigned int k = 0; k < npoints; k++)
                             {
                              m_w[f].push_back(-1);
                              m_t[f].push_back((1 / (double)srate)*k*(1e-6));
                             }
                        }     
                    }
                // Save the event we were checking for conicidences with in the previous setp
                int wch1 = -1;
                for (unsigned int d = 0; d < m_channels.size(); d++) if (m_channels.at(d) == ChNo.at(0)) wch1 = d;
                if (wch1 == -1)
                   {
                    std::cout << __FUNCTION__ << " ERROR: Cannot determine which channel this event corrsponds to" << std::endl;
                    return false;
                   }
                 m_w[wch1] = l_voltage.at(temp);
                 m_t[wch1].clear();
                 m_t[wch1].reserve(npoints);
                 for (unsigned int k = 0; k < npoints; k++) m_t[wch1].push_back((1/(double)srate)*k*(1e-6));
                 m_physt[wch1] = l_phystime.at(temp);
                 m_ordrt[wch1] = l_ordrtime.at(temp);
                 m_npoints.at(0) = npoints;
                 m_srate.at(0) = ceil(srate*1e6);
                 m_event = l_nEvntN.at(temp);
                 l_triggtimeVec.push_back(l_phystime.at(temp)/1e9);
                 m_tree->Fill();
                 ievent++;
                }
            if (!done) done = LGADBase::ProgressBar(floor((float)(l_voltage.size()+temp)/(float)(2*m_nchan)),ceil((float)l_voltage.size()/(float)m_nchan));
           }
       tevent += ievent;
      }
  std::cout << "--> Wrote " << tevent << " events into output file " << m_ofile->GetName() << ", Efficiency = " << (float)(ievent*m_nchan) / (float)TotHits << std::endl;
  LGADBase::CalcTrigFr(l_triggtimeVec, m_trigDt, m_trigFr, m_tree->GetEntries());

  delete b;
  m_event = tevent;
  m_ofile->Write();
  delete m_trigDt;
  delete m_trigFr;
  delete m_tree;
  m_ofile->Close();
  return true;
}
// --------------------------------------------------------------------------------------------------------------
// Readout function for the LeCroy scope on a txt file saving mode where each channel of each event is writen to a seperate txt file.
// Multiple scopes are supported with a non-limited number of channels as long as the run names are the same
bool LGADBase::WriteLabTXT(const char* dir, const char* name, const char* ext, int evt1, int evt2)
{
    // File name format string for scanf with placeholders for the dir, channel, the event number and extension
    const char* ifname_format = "%sC%d%s%05d.%s";
    std::istringstream iss(std::ios_base::app | std::ios_base::in | std::ios_base::out);
    double time;
    double voltage;
    char ifname[2048];
    std::string line;
    std::ifstream file;
    unsigned int npoints = 0; 
    std::vector <unsigned int> l_npoints;
    std::vector<double> l_triggtimeVec;
    l_npoints.clear();
    bool l_labview = false;
    int evtfd = -1;

    unsigned int nfiles = LGADBase::CountFiles(dir, ext);
    bool end = false;
   
    if (nfiles == 0)
       {
        std::cout << __FUNCTION__ << " ERROR: Unable to find any input files in directory: " << dir << std::endl;
        return false;
       }

    int count = 0;
    char oscillos[64];
    // Test the number of chanels
    findchannel:
    for (unsigned int i = 1; i < 65; i++) 
        { 
         sprintf(ifname, ifname_format, dir, i, name, evt1, ext);  // build the filename for the first event in the first channel
         file.open(ifname);
         if (file)
            {
             line.clear();
             memset(oscillos, 0, sizeof(oscillos));
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
                   memset(oscillos, 0, sizeof(oscillos));
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
                         std::cout << __FUNCTION__ << " ERROR: Failed recognizing waveform file header: " << line.c_str()<< std::endl;
                         return false;
                        }
                  }
             l_npoints.push_back(npoints);
             m_channels.push_back(i);
             m_nchan++;
             file.close();
            }
         else continue;
        }

    // Start event does not correspond to expectations, find first event
    if (m_nchan == 0)
       {
        std::cout << __FUNCTION__ << " WARNING: Given start event number " << evt1 << " in disagreement with files in: " << dir << ", attempting to recover... " << std::endl;
        std::vector<std::string> filenames = LGADBase::ListFileNames(dir, ext);
        std::string tr = name;
        tr = "C%*u" + tr + "%05u.%*s";
        sscanf((filenames.at(0)).c_str(), tr.c_str(), &evtfd);
        if (evtfd > 0)
           {
            std::cout << __FUNCTION__ << "          Found files starting from event " << evtfd << std::endl;
            evt1 = evtfd;
            evt2 += evtfd;
            goto findchannel;
           }
        else {
              std::cout << __FUNCTION__ << " ERROR: Unsaccsesful in finding first event in folder: " << dir << std::endl;
              return false;
             }
       }
    // Verify and adjust stop event
    if (evtfd < 0 && (evt2 > (int)(nfiles/m_nchan)-evt1))
       {
        evt2 = (nfiles/m_nchan)-evt1;
        std::cout << __FUNCTION__ << "          Adjusting stop event to available files " << evt2 << std::endl;
       }
    else if (evtfd > 0 && evt2 != 0 && (evt2 - evt1) >(int)(nfiles/m_nchan))
            { 
             evt2 = int(nfiles / m_nchan);
             std::cout << __FUNCTION__ << "          Adjusting stop event to available files " << evt2 << std::endl;
            }

    // time of the standard time_t type: the number of seconds since Jan 1, 1970 
    double trigtime_first = 0.0; // time of the first event
    m_trigtime = 0.0;

    // create the ntuple
    if (!LGADBase::CreateOutputFile(m_ofdir, m_ofname, m_channels))
       {
        std::cout << __FUNCTION__ << " ERROR: Output Ntuple creation failed: " << m_ofname << std::endl;
        return false;
       }
     
    if (!l_labview) LGADBase::SetScale(m_channels, m_nchan, &m_scale);

    int points = *min_element(l_npoints.begin(), l_npoints.end());

    bool comma = false;
    unsigned int oi = 0;
    unsigned int zr = 0;
    float oi1 = 0;
    float a = 1.;
    // Loop over all event files
    count = 0; // Display coiunter for user refernce
    unsigned int n = 0; // Writen events counter
    char VUnit[6]; // Voltage Units placeholder
    for (m_event = evt1; true; m_event++)  // NB: an infinite loop -- for (int m_event=0; m_event<10; ++m_event)
        {
         if ((evt2 != 0 && evt2 > evt1 && (int)m_event >= evt2) || (m_event*m_nchan == nfiles)) break;
         if (evt2 != 0 && evt2 != evtfd) LGADBase::ProgressBar(m_event-evt1, evt2 - evt1);
         else if (evtfd == -1 && evt2 == 0) LGADBase::ProgressBar(m_event-evt1, (nfiles/m_nchan)-evt1);
         else LGADBase::ProgressBar(m_event - evt1, ((nfiles / m_nchan) + evtfd) - evt1);
         if (end) break;
         bool first = false;
         bool bad = false; // Set per event, if one channel is bad, ignore the event
         for (unsigned int ich = 0; ich < m_nchan; ich++) // loop over all cahnnels
             {
              if (!bad) // event is bad, we ignore all additional cahnnels
                 {
                  m_w[ich].clear();
                  m_t[ich].clear();
                  m_w[ich].resize(points);
                  m_t[ich].resize(points);
                  m_npoints[ich] = points;
                  memset(ifname, '0', sizeof(ifname));
                  sprintf(ifname, ifname_format, dir, m_channels[ich], name, m_event, ext); // build the filename
                  file.open(ifname);
                  if (!file)
                     {
                      if (m_tree->GetEntries() == 0) 
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
                                 memset(VUnit, 0, sizeof(VUnit));
                                 sscanf(line.c_str(), "%*u;%*u;%*3s;%*2s;%*3s;\"%*3s  %*2s %*8s  %f%6s  %*f%*6s  %*u %*6s  %*4s %*3s %*4s\";%*u;%*s;\"%*s\";%*f;%*f;%*u;\"%*s\";%*f;%*f;%*f;%*u", &m_scale.at(ich), VUnit);
                                 if (std::string(VUnit) == "mV/div") m_scale.at(ich) = m_scale.at(ich)*0.001;
                                 line.clear();
                                 getline(file, line);
                                }
                       }
                  if (!first)                                                  // use the time stamp from the first channel only
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
                      if (m_trigtime < 0)
                         {
                          do m_trigtime = m_trigtime + 3600;
                          while (m_trigtime < 0);
                         }
                      else if (m_trigtime >= 3600 && n < 12) {trigtime_first = trigtime_first + 3600; m_trigtime = m_trigtime-3600;}
                      if (l_labview && std::string(oscillos) != "Tektronix") 
                         {
                          float rd =  (float)oi/a;
                          m_trigtime = m_trigtime + ((float)oi/a) - oi1;
                         }
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
         if (n > 0 && !bad) m_tree->Fill();
         count++; // Display counter, not taken into account
        }

    // Fix the numbers of events
    m_event = n / m_nchan;
    std::cout << "Wrote " << m_event << " events into output file " << m_ofile->GetName() << std::endl;
    LGADBase::CalcTrigFr(l_triggtimeVec, m_trigDt, m_trigFr, m_tree->GetEntries());

    m_ofile->Write();
    delete m_trigDt;
    delete m_trigFr;
    delete m_tree;
    m_ofile->Close();
    return true;
}
// --------------------------------------------------------------------------------------------------------------
bool LGADBase::WriteTectronixTXT(const char* dir, const char* name, const char* ext)
{
    static FILE* fd_file = NULL;  // internal
    std::string ifname(dir);
    long int lCurPos;
    char line[2048];
    std::vector<unsigned int> l_npoints;
    l_npoints.clear();

    // Open file
    ifname += name;
    ifname += "_Ch";
    unsigned int nel = 0;
    unsigned int npoints = 0;
    double volt = 0;
    double time = 0;
    for (unsigned int i = 1; i < 65; i++) 
        { 
         char fn[2];
         sprintf(fn, "%i.", i);
         ifname += fn;
         ifname += ext;
         fd_file = fopen(ifname.c_str(), "r");
         if (fd_file != NULL)
            { 
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
             if (size <= 0) { std::cout << __FUNCTION__ << " ERROR: txt file is zero size, skipping..." << std::endl; continue; }
             if(fgets(line, sizeof(line), fd_file) == NULL) return false;
             nel = sscanf(line, "%*s %*s %u %*s %lf %lf", &npoints, &time, &volt);
             if (nel !=3 && npoints == 0)
                 {
                  std::cout << __FUNCTION__ << " ERROR: Invalid parameters from ASCI file, skipping file" << std::endl;
                  continue;
                 }
             l_npoints.push_back(npoints);
             m_nchan++;
             m_channels.push_back(i);
             fclose(fd_file);
             ifname = ifname.substr(0, ifname.find_last_of("_Ch")+1);
            }
         else ifname = ifname.substr(0, ifname.find_last_of("_Ch")+1);
        }
    if (m_nchan == 0)
       {
        std::cout << __FUNCTION__ << " ERROR: Failed to find ny active channels within files: " << ifname << std::endl;
        return false;
       }
    
    // create the ntuple
    if (!CreateOutputFile(m_ofdir, m_ofname, m_channels))
       {
        std::cout << __FUNCTION__ << " ERROR: Output Ntuple creation failed: " << m_ofname << std::endl;
        return false;
       }

    SetScale(m_channels, m_nchan, &m_scale);

    for (unsigned int i = 0; i < m_nchan; i++)
        { 
         char fn[2];
         sprintf(fn, "%i.", m_channels.at(i));
         ifname += fn;
         ifname += ext;
         fd_file = fopen(ifname.c_str(), "r");
         unsigned int frame = 0;
         double interval = 0;
         double trigtime = 0;
         unsigned int trigpoint = 0;
         unsigned int np = 0;
         m_w[i].clear();
         m_npoints[i] = l_npoints[i];
         m_w[i].reserve(m_npoints.at(i));
         m_t[i].clear();
         m_t[i].reserve(m_npoints.at(i));
         if(fgets(line, sizeof(line), fd_file) == NULL) return false;
         nel = sscanf(line, "%*s %*s %*lf %*s %lf %lf", &time, &volt);
         if (nel == 2) 
            {
             m_t[i].push_back(time); m_w[i].push_back(volt);
             np++;
            }
         else { std::cout << __FUNCTION__ << " ERROR: Invalid parameters from ASCI file, skipping file" << ifname << std::endl; continue; }
         if(fgets(line, sizeof(line), fd_file) == NULL) return false;
         nel = sscanf(line, "%*s %*s %lf %*s %lf %lf", &interval, &time, &volt);
         if (nel == 3) 
            {
             m_srate.at(i) = (Long64_t)(ceil((double)1/interval));
             m_t[i].push_back(time); m_w[i].push_back(volt);
             np++;
            }
         else { std::cout << __FUNCTION__ << " ERROR: Invalid parameters from ASCI file, skipping file" << ifname << std::endl; continue; }
         if(fgets(line, sizeof(line), fd_file) == NULL) return false;
         nel = sscanf(line, "%*s %*s %u %*s %lf %lf", &trigpoint, &time, &volt);
         if (nel == 3) 
            {
             m_trigtime = trigpoint;
             m_t[i].push_back(time); m_w[i].push_back(volt);
             np++;
            }
         else { std::cout << "ERROR: Invalid parameters from ASCI file, skipping file" << ifname << std::endl; continue; }
         if(fgets(line, sizeof(line), fd_file) == NULL) return false;
         nel = sscanf(line, "%*s %*s %lf %*s %lf %lf", &trigtime, &time, &volt);
         if (nel == 3) 
            {
             m_triggTime.at(i) = trigtime;
             m_t[i].push_back(time); m_w[i].push_back(volt);
             np++;
            }
         else { std::cout << __FUNCTION__ << " ERROR: Invalid parameters from ASCI file, skipping file" << ifname << std::endl; continue; }
         if(fgets(line, sizeof(line), fd_file) == NULL) return false;
         nel = sscanf(line, " %*s %lf %lf", &time, &volt);
         if (nel == 2) 
            {
             m_t[i].push_back(time); m_w[i].push_back(volt);
             np++;
            }
         else { std::cout << __FUNCTION__ << " ERROR: Invalid parameters from ASCI file, skipping file" << ifname << std::endl; continue; }
         if(fgets(line, sizeof(line), fd_file) == NULL) return false;
         nel = sscanf(line, "%*s %*s %*lf %*s %lf %lf", &time, &volt);
         if (nel == 2) 
            {
             m_t[i].push_back(time); m_w[i].push_back(volt);
             np++;
            }
         else { std::cout << __FUNCTION__ << " ERROR: Invalid parameters from ASCI file, skipping file" << ifname << std::endl; continue; }
         if(fgets(line, sizeof(line), fd_file) == NULL) return false;
         if (std::string(line).find("\"FastFrame Count\"") != std::string::npos) 
            {
             nel = sscanf(line, "%*s %*s %u %*s %lf %lf", &frame, &time, &volt);
             if (nel == 3)
                {
                 m_t[i].push_back(time); 
                 m_w[i].push_back(volt);
                 np++;
                }
             else { std::cout << __FUNCTION__ << " ERROR: Invalid parameters from ASCI file, skipping file" << ifname << std::endl; continue; }
            }
         else {
               nel = sscanf(line, "%lf %lf", &time, &volt);
               if (nel == 2)
                  {
                   m_t[i].push_back(time);
                   m_w[i].push_back(volt);
                   np++;
                  }
               else { std::cout << __FUNCTION__ << " ERROR: Invalid parameters from ASCI file, skipping file" << ifname << std::endl; continue; }
              }
         for (unsigned int icount = 0; icount < frame; icount++)
             {
             if (icount %1000 == 0) printf("Converting event number %d in channel %d \n", icount, m_channels.at(i));
             if (np == 0)
                 { 
                  m_w[i].clear();
                  m_w[i].reserve(m_npoints.at(i));
                  m_t[i].clear();
                  m_t[i].reserve(m_npoints.at(i));
                 }
              while (np < m_npoints.at(i))
                    {
                     if(fgets(line, sizeof(line), fd_file) == NULL) return false;
                     nel = sscanf(line, "%lf %lf", &time, &volt);
                     if (nel == 2)
                        {
                         m_t[i].push_back(time); 
                         m_w[i].push_back(volt);
                         np++;
                        }
                    else { std::cout << __FUNCTION__ << " ERROR: Invalid parameters from ASCI file, skipping file" << ifname << std::endl; continue; }
                   }
             np = 0;
             m_event++;
            }
         m_tree->Fill();
         ifname = ifname.substr(0, ifname.find_last_of("_Ch") + 1);
        }

    m_ofile->Write();
    delete m_trigDt;
    delete m_trigFr;
    delete m_tree;
    m_ofile->Close();
    return true;
}
// --------------------------------------------------------------------------------------------------------------
// Function to merge tracking information from EUTelescope reconstruction
//
// Author: Luci�a Castillo Garci�a - lucia.castillo.garcia@cern.ch - IFAE-BARCELONA
//
bool LGADBase::CombineTrack(const char* dir, const char* name)
{
   // // Check the directory where the tracking file is and set the output combined file 
   // if (m_Trackdatadir.IsNull()) { std::cout << __FUNCTION__ << " ERROR: Input tracking data directory not set! Abording..." << std::endl; return false; }

   // if (m_Trackdataname.IsNull() && (m_instrument == LabTXT || m_instrument == TektronixScope))
   //    { 
   //     std::cout << "--> Warning: Input tracking data file name not set! Trying with default for instrument: --Trace-- " << std::endl;
   //     m_Trackdataname = "--Trace--";
   //    }
   // else if (m_Trackdataname.IsNull() && m_instrument == Sampic)
   //         {
   //          m_Trackdataname = m_Trackdatadir;
   //          TString sLast1 = m_Trackdataname[m_Trackdataname.Length() - 1];
   //          if (sLast1.Contains("/")) m_Trackdataname = m_Trackdataname(0, m_Trackdataname.Length() - 1);
   //          m_Trackdataname = m_Trackdataname(m_Trackdataname.Last('/') + 1, m_Trackdataname.Length() - m_Trackdataname.Last('/') - 1);
   //          std::cout << "--> Warning: Input tracking data file name not set! Trying with same name as containing folder: " << m_Trackdataname << std::endl;           
   //         }
   // else if (m_Trackdataname.IsNull()) { std::cout << __FUNCTION__ << " ERROR: Input data file name not set! Abording..." << std::endl; return false; }

   // if (!m_ofname.IsNull())
   //    {
   //     std::cout << "--> Output Ntuple name is combined using output data filename and tracking filename: ";
   //     if (m_instrument == LabTXT || m_instrument == TektronixScope || m_instrument == Sampic)
   //        {
   //         m_ofname = m_Trackdatadir;
   //         TString sLast2 = m_ofname[m_ofname.Length() - 1];
   //         if (sLast2.Contains("/")) m_ofname = m_ofname(0, m_ofname.Length() - 1);
   //         m_ofname = m_ofname(m_ofname.Last('/') + 1, m_ofname.Length() - m_ofname.Last('/') - 1);
   //        }
   //     else m_ofname = m_ofname + "_" + m_Trackdataname;
   //     std::cout << m_ofname << std::endl;
   //    }
   // if (m_ofdir.IsNull())
   //    {
   //     std::cout << "--> Warning: Output combined Ntuple location not set! Using input data diretory" << std::endl;
   //     m_ofdir = m_Trackdatadir;
   //    }
   // if (m_TrackComb && (m_Trackdatadir.IsNull() || m_Trackdataname.IsNull()))
   //    { 
   //     std::cout << "--> Warning: Tracking input file not defined, switching to non-test beam mode" << std::endl;
   //     m_TrackComb = false;
   //    }
   // 
   // // Open oscilloscope file
   // char scopename[2048];
   // strncpy(scopename, dir, sizeof(scopename));
   // strcat(scopename, m_dataname);
   // strcat(scopename, ".root");
   // TFile *scopefile = new TFile(scopename,"READ");
   // if (scopefile == NULL) { std::cout << __FUNCTION__ << " ERROR: Failed to open oscilloscope root file: " << scopename << std::endl; return false; }
   // else std::cout << "    Opening oscilloscope root file: " << std::string(scopename) << endl;
   // TTree *wfm = (TTree*)scopefile->Get("wfm");
   // if (!wfm) 
   //    {
   //     std::cerr << "Could not read oscilloscope tree" << std::endl;
   //     return -1;
   //    }
   // 
   // unsigned int nentries_scope = wfm->GetEntries();
   // std::cout << "Entries oscilloscope data file: " << nentries_scope << std::endl;    
   // 
   // // Declaration of leaf types
   // vector<double>  *t01, *w01, *t02, *w02, *t03, *w03, *t04, *w04;
   // Double_t        vScale01, vScale02, vScale03, vScale04;
   // UInt_t          nPoints01, nPoints02, nPoints03, nPoints04;
   // Long64_t        SnRate01, SnRate02, SnRate03, SnRate04;
   // UInt_t          EvnNo;

   //TBranch        *b_m_event;   //!
   // //Tree branches declaration
   // wfm->SetBranchStatus("*",0);
   // wfm->SetBranchStatus("t01",1);
   // wfm->SetBranchStatus("w01",1);
   // wfm->SetBranchStatus("vScale01",1);
   // wfm->SetBranchStatus("nPoints01",1);
   // wfm->SetBranchStatus("SnRate01",1);
   // wfm->SetBranchStatus("t02",1);
   // wfm->SetBranchStatus("w02",1);
   // wfm->SetBranchStatus("vScale02",1);
   // wfm->SetBranchStatus("nPoints02",1);
   // wfm->SetBranchStatus("SnRate02",1);
   // wfm->SetBranchStatus("t03",1);
   // wfm->SetBranchStatus("w03",1);
   // wfm->SetBranchStatus("vScale03",1);
   // wfm->SetBranchStatus("nPoints03",1);
   // wfm->SetBranchStatus("SnRate03",1);
   // wfm->SetBranchStatus("t04",1);
   // wfm->SetBranchStatus("w04",1);
   // wfm->SetBranchStatus("vScale04",1);
   // wfm->SetBranchStatus("nPoints04",1);
   // wfm->SetBranchStatus("SnRate04",1);
   // wfm->SetBranchStatus("EvnNo",1);
   // 
   // EvnNo = 0;
   // wfm->SetMakeClass(1);
   // wfm->SetBranchAddress("EvnNo", &EvnNo, &b_m_event);
   // 
   // // Open the tracking file        
   // char telescopename[2048];    
   // strncpy(telescopename, dir, sizeof(telescopename));
   // strcat(telescopename, name);
   // strcat(telescopename, ".root");
   // TFile *telescopefile = new TFile(telescopename,"READ");
   // if (telescopefile == NULL) { std::cout << __FUNCTION__ << " ERROR: Failed to open tracking root file: " << telescopename << std::endl; return false; }
   // else std::cout << "    Opening tracking root file: " << std::string(telescopename) << endl;
   // TTree *Tracks = (TTree*)telescopefile->Get("Tracks");
   // if (!Tracks) 
   //    {
   //     std::cerr << "Could not read telescope tree" << std::endl;
   //     return -1;
   //    }
   // TTree *ZeroSuppressed = (TTree*)telescopefile->Get("ZeroSuppressed");
   // if (!ZeroSuppressed) 
   //    {
   //     std::cerr << "Could not read telescope tree" << std::endl;
   //     return -1;
   //    }
   //    
   // unsigned int nentries_telescope = Tracks->GetEntries();
   // std::cout << "Entries tracking data file: " << nentries_telescope << std::endl;
   // 
   // unsigned int nentries_telescope_hits = ZeroSuppressed->GetEntries();
   // std::cout << "Entries tracking hits data file: " << nentries_telescope_hits << std::endl;
   // 
   // unsigned int nentries_diff = nentries_scope - nentries_telescope;
   // std::cout << "Difference in entries between oscilloscope and tracking data files: " << nentries_diff << std::endl;
   //    
   // // Definition of variables Tracks tree
   // Int_t eventNumber, triggerID, timestamp;
   // std::vector<double> *xPos, *yPos, *omega, *phi, *kinkx, *kinky, *chi2;
   // std::vector<int> *planeID, *trackID, *ndof;
   // // Definition of variables Hits tree
   // std::vector<int> *ID;
   // std::vector<int> *xPosH, *yPosH;//, *zPosH;
   // std::vector<double> *signalH;
   // std::vector<int> *timeH;
   // 
   // // List of branches Tracks tree
   // TBranch        *b_eventNumber;   //!
   // TBranch        *b_planeID;   //!
   // TBranch        *b_trackID;   //!
   // TBranch        *b_triggerID;   //!
   // TBranch        *b_timestamp;   //!
   // TBranch        *b_xPos;   //!
   // TBranch        *b_yPos;   //!
   // TBranch        *b_omega;   //!
   // TBranch        *b_phi;   //!
   // TBranch        *b_kinkx;   //!
   // TBranch        *b_kinky;   //!
   // TBranch        *b_chi2;   //!
   // TBranch        *b_ndof;   //!
   // // List of branches Hits tree
   // TBranch        *b_ID;   //!
   // TBranch        *b_xPosH;   //!
   // TBranch        *b_yPosH;   //!
   // //TBranch        *b_zPosH;   //!
   // TBranch        *b_signalH;   //!
   // TBranch        *b_timeH;   //!
   // 
   // // Set object pointer Tracks tree
   // planeID = 0;
   // trackID = 0;
   // xPos = 0;
   // yPos = 0;
   // omega = 0;
   // phi = 0;
   // kinkx = 0;
   // kinky = 0;
   // chi2 = 0;
   // ndof = 0;
   // // Set object pointer Hits tree
   // ID = 0;
   // xPosH = 0;
   // yPosH = 0;
   // //zPosH = 0;
   // signalH = 0;
   // timeH = 0;
   // 
   // // Tracks tree branches declaration
   // Tracks->SetMakeClass(1);
   // Tracks->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   // Tracks->SetBranchAddress("planeID", &planeID, &b_planeID);
   // Tracks->SetBranchAddress("trackID", &trackID, &b_trackID);
   // Tracks->SetBranchAddress("triggerID", &triggerID, &b_triggerID);
   // Tracks->SetBranchAddress("timestamp", &timestamp, &b_timestamp);
   // Tracks->SetBranchAddress("xPos", &xPos, &b_xPos);
   // Tracks->SetBranchAddress("yPos", &yPos, &b_yPos);
   // Tracks->SetBranchAddress("omega", &omega, &b_omega);
   // Tracks->SetBranchAddress("phi", &phi, &b_phi);
   // Tracks->SetBranchAddress("kinkx", &kinkx, &b_kinkx);
   // Tracks->SetBranchAddress("kinky", &kinky, &b_kinky);
   // Tracks->SetBranchAddress("chi2", &chi2, &b_chi2);
   // Tracks->SetBranchAddress("ndof", &ndof, &b_ndof);
   // // Tracks tree branches declaration
   // ZeroSuppressed->SetMakeClass(1);
   // ZeroSuppressed->SetBranchAddress("ID", &ID, &b_ID);
   // ZeroSuppressed->SetBranchAddress("xPos", &xPosH, &b_xPosH);
   // ZeroSuppressed->SetBranchAddress("yPos", &yPosH, &b_yPosH);
   // //ZeroSuppressed->SetBranchAddress("zPos", &zPosH, &b_zPosH);
   // ZeroSuppressed->SetBranchAddress("Signal", &signalH, &b_signalH);
   // ZeroSuppressed->SetBranchAddress("Time", &timeH, &b_timeH);

   // Int_t eventNum, evNum, trigger_id, timestampTr;
   // std::vector<double> xTr, yTr, omegaTr, phiTr, kinkxTr, kinkyTr, chi2Tr;
   // std::vector<int> plane_id, track_id, ndofTr;
   // 
   // std::vector<int> dutID;
   // std::vector<int> xHit, yHit, timeHit;//, zHit;
   // std::vector<double> signalHit;
   //     
   // // Printout registered Branches
   // TIter next(Tracks->GetListOfBranches());
   // TBranch* b2;
   // std::cout << "--> Cloning oscilloscope ROOT file and adding tracking Branches: " << std::endl << "\t";
   // int l_space = 0;
   // while ((b2 = (TBranch*)next()))
   //       {
   //        std::cout << std::setw(10) << std::left << b2->GetName();
   //        l_space++;
   //        if (l_space % 4 == 0) std::cout << std::endl << "\t";
   //       }
   // std::cout << std::endl;
   // 
   // TIter next2(ZeroSuppressed->GetListOfBranches());
   // TBranch* b3;
   // std::cout << "--> Cloning oscilloscope ROOT file and adding tracking Branches: " << std::endl << "\t";
   // int l_space2 = 0;
   // while ((b3 = (TBranch*)next2()))
   //       {
   //        std::cout << std::setw(10) << std::left << b3->GetName();
   //        l_space2++;
   //        if (l_space2 % 4 == 0) std::cout << std::endl << "\t";
   //       }
   // std::cout << std::endl;
   // 
   // // Create a new file
   // char newname[2048];
   // strncpy(newname, dir, sizeof(newname));
   // strcat(newname, m_ofname);
   // strcat(newname, ".root");
   // TFile *combinedfile = new TFile(newname,"RECREATE");
   // TTree *comb = wfm->CloneTree(0);
   // //comb->Branch("evNum",&evNum);
   // comb->Branch("plane_id","std::vector<int>",&plane_id);
   // comb->Branch("track_id","std::vector<int>",&track_id);
   // comb->Branch("trigger_id",&trigger_id);
   // comb->Branch("timestampTr",&timestampTr);
   // comb->Branch("xTr","std::vector<double>",&xTr);
   // comb->Branch("yTr","std::vector<double>",&yTr);
   // comb->Branch("omegaTr",&omegaTr);
   // comb->Branch("phiTr",&phiTr);
   // comb->Branch("kinkxTr",&kinkxTr);
   // comb->Branch("kinkyTr",&kinkyTr);
   // comb->Branch("chi2Tr",&chi2Tr);
   // comb->Branch("ndofTr",&ndofTr);
   // comb->Branch("dutID","std::vector<int>",&dutID);
   // comb->Branch("xHit","std::vector<int>",&xHit);
   // comb->Branch("yHit","std::vector<int>",&yHit);
   // comb->Branch("signalHit","std::vector<double>",&signalHit);
   // comb->Branch("timeHit","std::vector<int>",&timeHit);

   // // Loop over all events and fill the combined tree
   // for(unsigned int i=0; i<nentries_telescope; i++) 
   //     {
   //     if(i%1000 == 0)
   //         {
   //         //std::cout << "------------------------------------------------------------" << std::endl;
   //         std::cout << "Converting event number: " << i << std::endl;
   //         //std::cout << "------------------------------------------------------------" << std::endl;            
   //         }
   //     wfm->GetEntry(i);  
   //     eventNum = EvnNo;
   //     b_eventNumber->GetEntry(i);
   //     evNum = eventNumber;
   //     b_planeID->GetEntry(i);
   //     b_trackID->GetEntry(i);
   //     b_triggerID->GetEntry(i);
   //     trigger_id = triggerID;
   //     b_timestamp->GetEntry(i);
   //     timestampTr = timestamp;
   //     b_xPos->GetEntry(i);
   //     b_yPos->GetEntry(i);
   //     b_omega->GetEntry(i);
   //     b_phi->GetEntry(i);
   //     b_kinkx->GetEntry(i);
   //     b_kinky->GetEntry(i);
   //     b_chi2->GetEntry(i);
   //     b_ndof->GetEntry(i);
   //     b_ID->GetEntry(evNum-1);
   //     b_xPosH->GetEntry(evNum-1);
   //     b_yPosH->GetEntry(evNum-1);
   //     //b_zPosH->GetEntry(evNum-1);
   //     b_signalH->GetEntry(evNum-1);
   //     b_timeH->GetEntry(evNum-1);
   //     //comb->Fill();
   //     plane_id.clear();
   //     track_id.clear();
   //     xTr.clear();
   //     yTr.clear();
   //     omegaTr.clear();
   //     phiTr.clear();
   //     kinkxTr.clear();
   //     kinkyTr.clear();
   //     chi2Tr.clear();
   //     ndofTr.clear();
   //     dutID.clear();
   //     xHit.clear();
   //     yHit.clear();
   //     //zHit.clear();
   //     signalHit.clear();
   //     timeHit.clear();
   //     
   //     if(i%1000 == 0)
   //         {
   //         std::cout << "---> First check of events: " << std::endl;
   //         std::cout << "Oscilloscope event: " << eventNum << " EUTelescope event: " << evNum << std::endl;
   //         }
   //     if((eventNum - evNum) != 1)
   //         {
   //         for(unsigned int d=i+1; d<nentries_scope; d++) 
   //             {
   //             b_m_event->GetEntry(d);
   //             eventNum = EvnNo;
   //             if((eventNum - evNum) == 1)
   //                 {                
   //                 if(i%1000 == 0)
   //                     {
   //                     std::cout << "---> Checking events: " << std::endl;
   //                     std::cout << "Oscilloscope event: " << eventNum << " EUTelescope event: " << evNum << std::endl;
   //                     }
   //                 break;
   //                 }
   //             }
   //         }
   //     if(i%1000 == 0)
   //         {
   //          //Start building histograms
   //         std::cout << "---> Current event is: " << eventNum << " " << evNum << std::endl;
   //         std::cout << "Number of tracks: " << trackID->size()/4 << std::endl;                
   //         std::cout << "Number of hits in FE-I4: " << ID->size() << std::endl;                
   //         }
   //     for(unsigned int j=0; j<planeID->size(); j++)
   //         {
   //         plane_id.push_back(planeID->at(j));                
   //         track_id.push_back(trackID->at(j));                
   //         xTr.push_back(xPos->at(j));                
   //         yTr.push_back(yPos->at(j));                
   //         omegaTr.push_back(omega->at(j));                
   //         phiTr.push_back(phi->at(j));                
   //         kinkxTr.push_back(kinkx->at(j));                
   //         kinkyTr.push_back(kinky->at(j));                
   //         chi2Tr.push_back(chi2->at(j));                
   //         ndofTr.push_back(ndof->at(j));
   //         }
   // 
   //     if(ID->size() != 0){    
   //         for(unsigned int k=0; k<ID->size(); k++)
   //             {
   //             dutID.resize(3*(k+1)+k-1);
   //             dutID.resize(3*(k+1)+k,ID->at(k));
   //             dutID.resize(3*(k+1)+k+1);
   //             xHit.resize(3*(k+1)+k-1);
   //             xHit.resize(3*(k+1)+k,xPosH->at(k));
   //             xHit.resize(3*(k+1)+k+1);
   //             yHit.resize(3*(k+1)+k-1);
   //             yHit.resize(3*(k+1)+k,yPosH->at(k));
   //             yHit.resize(3*(k+1)+k+1);
   //             signalHit.resize(3*(k+1)+k-1);
   //             signalHit.resize(3*(k+1)+k,signalH->at(k));
   //             signalHit.resize(3*(k+1)+k+1);
   //             timeHit.resize(3*(k+1)+k-1);
   //             timeHit.resize(3*(k+1)+k,timeH->at(k));
   //             timeHit.resize(3*(k+1)+k+1);
   //             }
   //     }else{
   //         dutID.resize(4);
   //         xHit.resize(4);
   //         yHit.resize(4);
   //         //zHit.resize(4);
   //         signalHit.resize(4);
   //         timeHit.resize(4);
   //     }
   //     comb->Fill();
   //     }
   // comb->Print();
   // unsigned int nentries_comb = comb->GetEntries();
   // std::cout << "Entries combined data file: " << nentries_comb << std::endl;
   // combinedfile->Write();
   // combinedfile->Close();

   // delete scopefile;
   // delete telescopefile;
   // delete combinedfile;
   return true;
}
// --------------------------------------------------------------------------------------------------------------
// Function to read Lecroy oscilloscope data from DESY test beams
// Author: Luc�a Castillo Garc�a - lucia.castillo.garcia@cern.ch - IFAE-BARCELONA
bool LGADBase::WriteLecroyBinary(const char* dir, const char* name, const char* ext, unsigned int evt1, unsigned int evt2)
{
    std::vector<unsigned int> nevt_seq;  // internal, number of events per sequence, identical for all channels
    static FILE* fd_file = NULL;  // internal
    char* fileBuf;  // internal
    size_t lCurPos;
    int nrec = 0;

    std::vector<float> yscale; //VERTICAL_GAIN
    std::vector<float> yoffset; //VERTICAL_OFFSET
    std::vector<float> xoffset; //HORIZ_OFFSET
    std::vector<float> xinterval; //HORIZ_INTERVAL
    std::vector<unsigned int> ncount; //NOM_SUBARRAY_COUNT
    std::vector<unsigned int> warray; //WAVE_ARRAY_1
    std::vector<unsigned int> l_npoints;
    std::vector<Long64_t> l_srate;

    char ifname[2048];    
    char* buscador;
    unsigned int ievent = 0;
    int scope = 0;
    char oscillos[1024];
    m_nchan = 0;
    std::string line;
    int ncycle = 0;
    
    const char* search_ch1 = "Channel 1 will be saved";
    const char* search_ch2 = "Channel 2 will be saved";
    const char* search_ch3 = "Channel 3 will be saved";
    const char* search_ch4 = "Channel 4 will be saved";
    char* search;
    const char* buscador1 = "NOM_SUBARRAY_COUNT";
    const char* buscador2 = "WAVE_ARRAY_1";
    const char* buscador3 = "VERTICAL_OFFSET";
    const char* buscador4 = "HORIZ_OFFSET";
    const char* buscador5 = "HORIZ_INTERVAL";
    const char* buscador6 = "VERTICAL_GAIN";

    yscale.clear();
    yoffset.clear();
    nevt_seq.clear();
    l_npoints.clear();
    l_srate.clear();
    xinterval.clear();
    xoffset.clear();
    ncount.clear();
    warray.clear();

    // Open file
    strncpy(ifname, dir, sizeof(ifname));
    strcat(ifname, name);
    strcat(ifname, ".txt");
    fd_file = fopen(ifname, "r");
    if (fd_file == NULL) { std::cout << __FUNCTION__ << " ERROR: Failed to open txt file: " << ifname << std::endl; return false; }
    else std::cout << "    Opening text file: " << std::string(ifname) << std::endl;
    
    // GetFile size
    lCurPos = ftell(fd_file);
    fseek(fd_file, 0, 2);
    size_t size = ftell(fd_file);
    fseek(fd_file, lCurPos, 0);
    if (size <= 0) { std::cout << __FUNCTION__ << " ERROR: txt file is zero size, skipping..." << std::endl; return false; }
    
    // Figure out number of recorded channels  
    for(int i=0;i<9;i++)
        {    
        search = fgets(ifname, sizeof(ifname), fd_file);
        if ((strncmp (search_ch1,search,23) == 0) || (strncmp (search_ch2,search,23) == 0) || (strncmp (search_ch3,search,23) == 0) || (strncmp (search_ch4,search,23) == 0))
            {
            std::cout << search << std::endl;
            if (m_nchan == 0) m_channels.push_back((scope * 4)+ 1);
            else m_channels.push_back(m_channels.back() + 1); // Same scope    
            m_nchan++;
            }        
        }
    std::cout << __FUNCTION__ << " INFO: Number of recorded channels: " << m_nchan << "\n";
        
    //Get information from txt file
    unsigned int iline;
    std::string number;
    std::string number2;
    std::string number3;
    std::string number4;
    std::string number5;
    std::string number6;
    unsigned int num;
    unsigned int num2;
    float num3;
    float num4;
    float num5;
    float num6;
    unsigned int nev;            
    while (fgets(ifname, sizeof(ifname), fd_file))
        {
        unsigned int nel; // npts;
        nel = sscanf(ifname, "%s,%*s,%*s,%*s", oscillos);        
        iline = 11;
        for(unsigned ichan = 1; ichan < m_nchan+1; ichan++)
           {
            for(; iline < (58*ichan+10); iline++)
               {
                buscador = fgets(ifname, sizeof(ifname), fd_file);
                if (strncmp (buscador1,buscador,18) == 0)
                   {
                    number = buscador;
                    num = std::stoi(number.substr(21));
                    ncount.push_back(num);
                   }
                if (strncmp (buscador2,buscador,12) == 0)
                   {
                    number2 = buscador;
                    num2 = std::stoi(number2.substr(21));
                    warray.push_back(num2);
                   }
                if (strncmp (buscador3,buscador,15) == 0)
                   {
                    number3 = buscador;
                    num3 = std::stof(number3.substr(21));
                    yoffset.push_back(num3);
                   }    
                if (strncmp (buscador4,buscador,12) == 0)
                   {
                    number4 = buscador;
                    num4 = std::stof(number4.substr(21));
                    xoffset.push_back(num4);
                   }        
                if (strncmp (buscador5,buscador,14) == 0)
                   {
                    number5 = buscador;
                    num5 = std::stod(number5.substr(21));
                    xinterval.push_back(num5);
                    l_srate.push_back(1/num5);
                   }
                if (strncmp (buscador6,buscador,13) == 0)
                   {
                    number6 = buscador;
                    num6 = std::stof(number6.substr(21));
                    yscale.push_back(num6);
                   }
                }
            l_npoints.push_back(num2/num);
            nrec++;
            }        
        ncycle++;                   
        nev = num;
        nevt_seq.push_back(nev); // Write the number of events per sequence
        memset(&ifname[0], 0, sizeof(ifname));
        }
        if (m_verbose == 2) 
           {
            std::cout << __FUNCTION__ << " INFO: Number of events in each cycle: " << nev << std::endl;
            std::cout << __FUNCTION__ << " INFO: Size of vector number of events per cycle: " << nevt_seq.size() << std::endl;
           }
        
        std::cout << " INFO: Number of cycles in units of 4096 events: " << ncycle << std::endl;

    if (nevt_seq.size() != (nrec / m_nchan)) {std::cout << __FUNCTION__ << " ERROR: Unable to correctly determine number of cycles and total number of events" << std::endl; return false; }
    for (unsigned int i = 0; i < nevt_seq.size(); i++) ievent += nevt_seq.at(i);
    fclose(fd_file);
    
    // Create the data filename and open
    strncpy(ifname, dir, sizeof(ifname));
    strcat(ifname, name);
    if (sizeof(ifname) == 0) {std::cout << __FUNCTION__ << " ERROR: Could not properly form ASCII file name" << std::endl; return false;};
    strcat(ifname, ".");
    strcat(ifname, ext);
    if ((fd_file = fopen(std::string(ifname).c_str(), "rb")) == NULL) { std::cout << __FUNCTION__ << " ERROR: Could not open file " << std::string(ifname) << std::endl; return false; }
    else std::cout << "    Opening binary file: " << std::string(ifname) << std::endl;
     
    // Check the file size of the dat file
#ifdef _WIN32
    lCurPos = _ftelli64(fd_file);
    _fseeki64(fd_file, 0, 2); 
    size = _ftelli64(fd_file);
    _fseeki64(fd_file, lCurPos, 0);
#else
    lCurPos = ftell(fd_file);
    fseek(fd_file, 0, 2);
    size = ftell(fd_file);
    fseek(fd_file, lCurPos, 0);
#endif

    if (size <= 0) 
       { 
        std::cout << __FUNCTION__ << " ERROR: Failed getting data file size" << std::endl;
        return false;
       }
    size_t exp_size = 0; 

    // Get the expected size of the binary file
    for (unsigned int d = 0; d <nevt_seq.size(); d++)
        {
         // 16 bits as a header, 3 bits at the end of each sequece, number of points, number of event in a cycle
         for (unsigned f = 0; f < m_nchan; f++) exp_size += (16 + 3 + l_npoints.at(m_nchan*d+f)*nev);
        }
    if (size != exp_size)
       { 
        std::cout << __FUNCTION__ << " ERROR: Data file size not corresponding to expected value" << std::endl; 
        return false; 
       }

    // create the ntuple
    if (!CreateOutputFile(m_ofdir, m_ofname, m_channels))
       {
        std::cout << __FUNCTION__ << " ERROR: Output Ntuple creation failed: " << m_ofname << std::endl;
        return false;
       }

    // Copy  file to buffer
    // need to cut the file in chenks of events if it is big, cannot copy it in memory at once!
    // Reading sequence by sequence
    printf("Starting loop over %d events on %d channels, number of samples per channel and event: %d \n", ievent, m_nchan, l_npoints.at(0));
    lCurPos = 0;
    int offst = 0;
    int firsthead = 16;
    int trail = 0;
    size_t calc_size = 0;
    size_t prev_size = 0;
    for (unsigned int d = 0; d <nevt_seq.size(); d++)
        {
         exp_size = 0;
         calc_size = 0;
         for (unsigned f = 0; f < m_nchan; f++) exp_size += (16 + 3 + l_npoints.at(m_nchan*d+f)*nevt_seq.at(d));
         std::cout << "Expected size of sequence in binary file " << exp_size << std::endl;
         fileBuf = (char*)malloc(exp_size);
         size_t result = fread(fileBuf, exp_size, 1, fd_file);
         if (fileBuf == NULL || result != 1) { std::cout << __FUNCTION__ << " ERROR: Buffer is empty, check data file!!!" << std::endl; return false; }
         // Start Converting events in buffer 
         if(offst == 0) offst += firsthead * 4;
         else offst += firsthead*4 + trail*4;    
         for (unsigned int i = 0; i < nevt_seq.at(d); i++)
             {        
              if (m_event %1000 == 0) printf("Converting event number: %d \n", m_event);
              if (m_event < evt1) continue;
              if (m_event > evt2 && evt2 > 0) break;         
              for (unsigned int j = 0; j < m_nchan; j++)
                  {
                   m_scale[j] = yscale[m_nchan*d+j];
                   m_npoints[j] = l_npoints[m_nchan*d+j];
                   m_srate[j] = l_srate[m_nchan*d+j];
                  int offset = 0;
                  if (j == 0) offset += firsthead + l_npoints[m_nchan*d+j]*i;
                  else offset += firsthead*(j+1) + trail*j + l_npoints[m_nchan*d+j]*i;
                  for (unsigned int h = 0; h < j; h++) offset += l_npoints.at(m_nchan*d+h)*nevt_seq.at(d);
                  m_w[j].clear();
                  m_w[j].reserve(l_npoints.at(m_nchan*d+j));
                  m_t[j].clear();
                  m_t[j].reserve(l_npoints.at(m_nchan*d+j));
                  for (unsigned int w = 0; w < (l_npoints[m_nchan*d+j]); w++)
                      {
                       if ((w > (l_npoints[m_nchan*d + j]-4)) && (((fileBuf[w+offset] * yscale[m_nchan*d+j]) - yoffset[m_nchan*d+j]) > 0.04)) 
                          {
                           break;
                          } 
                       else if ((w < 10) && ((((fileBuf[w+offset] * yscale[m_nchan*d+j]) - yoffset[m_nchan*d+j]) > 0.04) || (((fileBuf[w+offset] * yscale[m_nchan*d+j]) - yoffset[m_nchan*d+j]) < -0.04)))
                              {
                              }
                       else {
                             // amplitude of signal is given in volts
                             m_w[j].push_back( ((fileBuf[w+offset] * yscale[m_nchan*d+j]) - yoffset[m_nchan*d+j])); 
                             m_t[j].push_back(((w) * xinterval[m_nchan*d+j]) + xoffset[m_nchan*d+j]);
                            }
                      }                
                  offst += l_npoints[m_nchan*d+j];            
                 }                
              m_event++;
              m_tree->Fill(); // fill every event
             }
         trail = 3;
         calc_size = offst + trail*4 - prev_size;
         std::cout << "Calculated size of sequence in binary file " << calc_size << std::endl;
         prev_size += calc_size;
         free(fileBuf);
         if (m_event > evt2 && evt2 > 0) break;
        }

    fclose(fd_file);
    m_event = m_event - 1;
    m_ofile->Write();
    delete m_trigDt;
    delete m_trigFr;
    delete m_tree;
    m_ofile->Close();
    return true;
}
// --------------------------------------------------------------------------------------------------------------
bool LGADBase::WriteTekBinary(const char* dir, const char* name, const char* ext, unsigned int evt1, unsigned int evt2)
{
}
// --------------------------------------------------------------------------------------------------------------
// Function to read binary files from multiple oscilloscopes supporting up to 64 channels, or a total of 16 non-connected scopes or 
// 31 dasy chaned scopes with common clock and syncronistion with Agilent readout code option for multimanaging.
bool LGADBase::WriteAgilentBinary(const char* dir, const char* name, const char* ext, unsigned int evt1, unsigned int evt2)
{
    std::vector<unsigned int> nevt_seq;  // internal, number of events per sequence, identical for all channels
    static FILE* fd_file = NULL;  // internal
    short int* fileBuf;  // internal
    size_t lCurPos;
    std::vector<float> yscale;
    std::vector<float> yoffset;
    std::vector<unsigned int> l_npoints;
    std::vector<Long64_t> l_srate;
    char ifname[2048];
    char sn0[1024], sn1[1024], sn2[1024];
    char time0[1024], time1[1024], time2[1024];
    unsigned int ievent = 0;
    int nrec = 0;
    int scope = 0;
    m_nchan = 0;

    yscale.clear();
    yoffset.clear();
    nevt_seq.clear();
    l_npoints.clear();
    l_srate.clear();

    // Open file
    strncpy(ifname, dir, sizeof(ifname));
    strcat(ifname, name);
    strcat(ifname, ".txt");
    fd_file = fopen(ifname, "r");
    if (fd_file == NULL) { std::cout << __FUNCTION__ << " ERROR: Failed to open txt file: " << ifname << std::endl; return false; }
    // GetFile size
    lCurPos = ftell(fd_file);
    fseek(fd_file, 0, 2);
    size_t size = ftell(fd_file);
    fseek(fd_file, lCurPos, 0);
    if (size <= 0) { std::cout << __FUNCTION__ << " ERROR: txt file is zero size, skipping..." << std::endl; return false; }
    while (fgets(ifname, sizeof(ifname), fd_file))
          {
           unsigned int nel, npts, nev;
           float yinc, yorg;
           double sr;
           nel = sscanf(ifname, "%*i,%*i,%i,%*i,%lf,%*f,%*f,%f,%*f,%*f,%*d,%*f,%*f,%*f,%f,%*[^,],%[^,],%[^,],%*d,%*d,%*d,%*d,%*f,%*f,%d", &npts, &sr, &yinc, &yorg, time2, sn2, &nev);
           if (!(nel == 7 && npts > 0 && nev > 0 && sr > 0))
              { 
               std::cout << __FUNCTION__ << " ERROR: Invalid parameters from ASCI file" << std::endl; 
               return false; 
              }
           // Figure out number of channels if not known
           if (m_nchan == 0)
              {
               yscale.push_back(yinc);
               yoffset.push_back(yorg);
               l_npoints.push_back(npts);
               l_srate.push_back(1/sr);
               if (strncmp(sn1, sn2, 1024) != 0 && strncmp(time1, time2, 1024) != 0 && strncmp(sn0, sn2, 1024) != 0) // every time we change oscilloscope
                  {
                   memcpy(sn1, sn2, 1024);
                   memcpy(time1, time2, 1024);
                   if (scope == 0) // first oscilloscope
                      {
                       nevt_seq.push_back(nev);
                       memcpy(sn0, sn2, 1024);
                       memcpy(time0, time2, 1024);
                      }
                   m_channels.push_back((scope * 4)+ 1);
                   scope++;
                  }
               else if (strncmp(sn1, sn2, 1024) == 0 && strncmp(time1, time2, 1024) == 0 && nev == nevt_seq.at(0))
                       {
                        m_channels.push_back(m_channels.back() + 1); // Same scope 
                       }
               else if (strncmp(sn0, sn2, 1024) == 0 && strncmp(time0, time2, 1024) != 0 && npts == l_npoints.at(0) && // back to first scope
                        ceil(1/sr) == l_srate.at(0) && yinc == yscale.at(0) && yorg == yoffset.at(0))
                      {
                       m_nchan = m_channels.size();
                      }
               else { 
                     std::cout << __FUNCTION__ << " ERROR: Not able to define number of used channels!! " << std::endl; 
                     return false; 
                    }
              }           
           if (m_nchan !=0 && (nrec % m_channels.size() == 0)) nevt_seq.push_back(nev); // Wrtite the number of events per sequence
           nrec++;
           memset(&ifname[0], 0, sizeof(ifname));
          }
    if (m_nchan == 0) m_nchan = nrec;
    if (nevt_seq.size() != (nrec / m_nchan)) {std::cout << __FUNCTION__ << " ERROR: Unable to correctly determine number of cycles and total number of events" << std::endl; return false; }
    for (unsigned int i = 0; i < nevt_seq.size(); i++) ievent += nevt_seq.at(i);
    fclose(fd_file);

    // Create the data filename and open
    strncpy(ifname, dir, sizeof(ifname));
    strcat(ifname, name);
    if (sizeof(ifname) == 0) {std::cout << __FUNCTION__ << " EROOR: Could not properly form ASCI file name" << std::endl; return false;};
    strcat(ifname, ".");
    strcat(ifname, ext);
    if ((fd_file = fopen(std::string(ifname).c_str(), "rb")) == NULL) { std::cout << "ERROR: Could not open file " << std::string(ifname) << std::endl; return false; }
    else std::cout << "    Opening binary file: " << std::string(ifname) << std::endl;
     
    // Check the file size of the dat file
#ifdef _WIN32
    lCurPos = _ftelli64(fd_file);
    _fseeki64(fd_file, 0, 2); 
    size = _ftelli64(fd_file);
    _fseeki64(fd_file, lCurPos, 0);
#else
    lCurPos = ftell(fd_file);
    fseek(fd_file, 0, 2);
    size = ftell(fd_file);
    fseek(fd_file, lCurPos, 0);
#endif
    if (size <= 0) 
       { 
        std::cout << __FUNCTION__ << " ERROR: Faild getting data file size" << std::endl;
        return false;
       }
    size_t exp_size = 0; 
    for (unsigned f = 0; f < m_nchan; f++) exp_size += (l_npoints.at(f) + 1) * 2 * ievent;
    if (size != exp_size)
       { 
        std::cout << __FUNCTION__ << " ERROR: Data file size not corresponding to expected value" << std::endl; 
        return false; 
       }

    // create the ntuple
    if (!CreateOutputFile(m_ofdir, m_ofname, m_channels))
       {
        std::cout << __FUNCTION__ << " ERROR: Output Ntuple creation failed: " << m_ofname << std::endl;
        return false;
       }

    // Copy  file to buffer
    // need to cut the file in chenks of events if it is big, cannot copy it in memory at once!
    // Reading sequence by sequence
    printf("Starting loop over %d events on %d channels, number of samples per channel and event: %d \n", ievent, m_nchan, l_npoints.at(0));
    lCurPos = 0;
    
    for (unsigned int d = 0; d <nevt_seq.size(); d++)
        {
         exp_size = 0;
         for (unsigned f = 0; f < m_nchan; f++) exp_size += (l_npoints.at(f) + 1) * 2 * nevt_seq.at(d);
         fileBuf = (short int*)malloc(exp_size);
         size_t result = fread(fileBuf, exp_size, 1, fd_file);
         if (fileBuf == NULL || result != 1) { std::cout << __FUNCTION__ << " ERROR: Buffer is empty, check data file!!!" << std::endl; return false; }
          // Start Converting events in buffer
         for (unsigned int i = 0; i < nevt_seq.at(d); i++)
             { 
              if (m_event %1000 == 0) printf("Converting event number: %d \n", m_event);
              if (m_event < evt1) continue;
              if (m_event > evt2 && evt2 > 0) break;
              for (unsigned int j = 0; j < m_nchan; j++)
                  {  
                   if (m_event == evt1)
                      { 
                       m_scale[j] = yscale[j] / 2.64674e-4;
                       m_npoints[j] = l_npoints[j];
                       m_srate[j] = l_srate[j];
                      }
                   int offst = 0;
                   for (unsigned int h = 0; h < j; h++) offst += l_npoints.at(h)*nevt_seq.at(d) + i;
                   m_w[j].clear();
                   m_w[j].reserve(l_npoints.at(j));
                   m_t[j].clear();
                   m_t[j].reserve(l_npoints.at(j));
                   for (unsigned int w = 0; w < l_npoints.at(j); w++)
                       {
                        m_w[j].push_back(fileBuf[offst + w] * yscale[j] + yoffset[j]);
                        m_t[j].push_back((double)w/l_srate.at(j));
                       }
                  }   
              m_event++;
              m_tree->Fill(); // fill every event
             }
         free(fileBuf);
         if (m_event > evt2 && evt2 > 0) break;
        }

    fclose(fd_file);
    m_event = m_event - 1;
    m_ofile->Write();
    delete m_trigDt;
    delete m_trigFr;
    delete m_tree;
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
        a = RecursMkDir(print.c_str());
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
    std::cout << "    Creating output ROOT file: " << m_ofile->GetName() << std::endl;
    m_tree = new TTree("wfm", "Timing part of the LGAD measurements");
    m_tree->SetAutoFlush(10000); // Set autoflush to 30 MB
    LGADBase::SetVectorSize(nchan.size());

    for (unsigned int ich = 0; ich < nchan.size(); ich++)
        {
         if (m_instrument == Sampic)
            {
             m_tree->Branch(Form("physt%02u", nchan.at(ich)), &(m_physt.at(ich)), "m_physt.at(ich)/D");
             m_tree->Branch(Form("ordrt%02u", nchan.at(ich)), &(m_ordrt.at(ich)), "m_ordrt.at(ich)/D");
            }
         m_tree->Branch(Form("t%02u", nchan.at(ich)), &(m_t.at(ich)));
         m_tree->Branch(Form("w%02u", nchan.at(ich)), &(m_w.at(ich)));
         if (m_instrument == InfiniiumScope || m_instrument == TektronixBinary || m_instrument == LabTXT || m_instrument == TektronixScope || m_instrument == LeCroyWRBin)
            {
             m_tree->Branch(Form("vScale%02u", nchan.at(ich)), &(m_scale.at(ich)), "m_scale.at(ich)/F");
             m_tree->Branch(Form("nPoints%02u", nchan.at(ich)), &(m_npoints.at(ich)), "m_npoints.at(ich)/i");
             m_tree->Branch(Form("SnRate%02u", nchan.at(ich)), &(m_srate.at(ich)), "m_srate.at(ich)/L");
             if (m_instrument == TektronixScope) m_tree->Branch(Form("triggTime%02u", nchan.at(ich)), &(m_triggTime.at(ich)), "m_triggTime.at(ich)/D");
            }
        }
    if (m_instrument == LabTXT || m_instrument == TektronixScope) m_tree->Branch("trigtime", &m_trigtime, "m_trigtime/D");
    m_tree->Branch("EvnNo", &m_event, "m_event/i");
    if (m_instrument == Sampic)
       {
        m_tree->Branch("nPoints", &(m_npoints.at(0)), "m_npoints.at(0)/i");
        m_tree->Branch("SnRate", &(m_srate.at(0)), "m_srate.at(0)/L");
       }

    m_trigDt = new TH1F("trigDt", "DT time", 6000, 0.0, 60.0);
    m_trigFr = new TH1F("trigFr", "Trigger Frequency", 50000, 0.0, 5000);

    // Printout registered Branches
    TIter next(m_tree->GetListOfBranches());
    TBranch* b;
    std::cout << "--> Booking Branches: " << std::endl << "\t";
    int l_space = 0;
    while ((b = (TBranch*)next()))
          {
           std::cout << std::setw(10) << std::left << b->GetName();
           l_space++;
           if (l_space % 4 == 0) std::cout << std::endl << "\t";
          }
  if (l_space % 4 != 0) std::cout << std::endl;

  TListIter l_next_object(gDirectory->GetList());
  l_next_object.Reset();
  TObject* l_obj;
  std::cout << "\r--> Booking Histos: " << std::endl << "\t";
  l_space = 0;
  while ((l_obj=l_next_object()))
        {
         TString l_objname = l_obj->GetName();
         if (l_objname.Contains("wfm")) continue;
         std::cout << std::setw(21) << std::left << l_objname ;
         l_space++;
         if (l_space %5 == 0) cout << std::endl << "\t" ;
        }
  if (l_space % 5 != 0) std::cout << std::endl;
  std::cout << std::endl;

   return true;
}
// --------------------------------------------------------------------------------------------------------------
void LGADBase::CalcTrigFr ( std::vector<double>  EvTrigTime, TH1F* TrigPer, TH1F* TrigFrq, unsigned int entriesNo)
{
  float dT1 = 0;
  std::vector<double> DeltaTrigTime;
  for (unsigned int hp = 0; hp < EvTrigTime.size() - 1; hp++)
      {
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
                if (ma < (EvTrigTime.size() - 1))
                   {
                    count++;
                    continue;
                   }
                else {
                      if (tbin == 0) tbin = 1;
                      dT1 = (float)count/tbin;
                      TrigFrq->Fill(dT1);
                      ka = ma;
                     }
               }
            else {
                  tbin = EvTrigTime.at(ma) - EvTrigTime.at(ka);
                  dT1 = (float)count/tbin;
                  TrigFrq->Fill(dT1);
                  ka = ma - 1;
                  break;
                 }
           }
      }

  dT1 = EvTrigTime.back() - EvTrigTime.at(0);
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
