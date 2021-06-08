/*
* LGADBase.cxx
*
*
*      Author: Gkougkousis Evangelos - Leonidas
*              egkougko@cern.ch
*               IFAE-BARCELONA
*/

// Modify code to be able to rerun processed Ntuple and produce waveform shape
// Write Canvas with the selected fit
// Modify mean shape and channel properties only for selected channel
// Add board, transimpedence, amp, gain, rate, ect per channel at final Ntuple
// Add CFD ToT as a vector for all CFD values
// Add a CFD - CFD for all bins and all channels
// Add Transimpedence consideration on charge
// Fix small signals charge issue

#include "../LGADUtils/LGADBase.h"

//#if !defined(__CINT__)
//ClassImp(LGADBase);
//#endif
//
//#if !defined(__CLING__)
//ClassImp(LGADBase);
//#endif
//
//#ifdef __CINT__
//#pragma link C++ class LGADBase;
//#endif
//
//#ifdef __ROOTCLING__
//#pragma link C++ class LGADBase;
//#endif

LGADBase::LGADBase()
{
    // This is the default constructor, so we use default settings
    m_instrument = Unasigned;
    Initialize();
    m_nchan = 0;
    m_event = 0;
    m_trigtime = (double)std::time(0);
    m_fitopt = "root";
}
// --------------------------------------------------------------------------------------------------------------
LGADBase::~LGADBase()
{
}
// --------------------------------------------------------------------------------------------------------------
void LGADBase::SetFitMethode(std::string method)
{
    if (method == "RooFit" || method == "root" || method == "rootInt") m_fitopt = method;
    else {
          m_fitopt = "root";
          std::cout << "ERROR: methode type " << method
                    << " not allowed. Available fit options are root and RooFit. Falling back to "
                    << m_fitopt << std::endl;
    }
}
// --------------------------------------------------------------------------------------------------------------
void LGADBase::SetInstrument(AqInstrument instr)
{

    if (instr == LeCroyWRBin || instr == Unasigned || instr == Sampic || instr == LabTXT || instr == InfiniiumScope || instr == TektronixScope) m_instrument = instr;
    else {
          std::cout << __FUNCTION__ << " ERROR: instrument value " << instr
                    << " not allowed. Allowed values are in ascending order: Sampic, LecroyWP725Zi, TektronixScope, InfiniumScope or Undefined. Falling back to "
                    << m_instrument << std::endl;
         }
}
// --------------------------------------------------------------------------------------------------------------
AqInstrument LGADBase::GetInstrument()
{
    if (m_instrument == Unasigned) std::cout << __FUNCTION__ << " WARNING: Instrument is not asigned!!" << std::endl;
    return m_instrument;
}
// --------------------------------------------------------------------------------------------------------------
void LGADBase::SetTrackComb(bool comb)
{
    m_TrackComb = comb;
}
// --------------------------------------------------------------------------------------------------------------
void LGADBase::SetTransFileName(TString filename) 
{
     m_TransFileName = filename;
}
//---------------------------------------------------------------------------------------------------------------
void LGADBase::SetDUTName(int ChId, std::string Name)
{
    if (0 < ChId && ChId > 64)
       {
        if (Name.length() != 0 && !(trim(Name, " ")).empty()) m_DUTChsNames.push_back(std::make_pair(ChId, Name));
        else std::cout << __FUNCTION__ << " WARNING: Cannot set empty name for channel " << ChId << "!" << std::endl;
       }
    else std::cout << __FUNCTION__ << " WARNING: Cahnnel ID not supported, 1 - 64 channels allowed." << std::endl;
}
//---------------------------------------------------------------------------------------------------------------
void LGADBase::SetDUTNames(std::vector<std::string> DUTNames)
{
    for (unsigned int k = 0; k < DUTNames.size(); k++)
        {
         if ((DUTNames.at(k)).length() != 0 && !(trim((DUTNames.at(k)), " ")).empty()) m_DUTChsNames.push_back(std::make_pair(-99, (DUTNames.at(k))));
         else std::cout << __FUNCTION__ << " WARNING: Cannot set empty name on channel name vector entry " << k << ", ignoring..." << std::endl;
        }
}
//---------------------------------------------------------------------------------------------------------------
void LGADBase::SetDUTBoard(int ChId, AqBoard Brd)
{
    if (0 < ChId && ChId > 64)
       {
        if (Brd == SingleCh || Brd == FourCh || Brd == IN2P3 || Brd == KU || Brd == SiPM) m_DUTChsBrd.push_back(std::make_pair(ChId, Brd));
        else std::cout << __FUNCTION__ << " WARNING: Board type for channel " << ChId << " not supported!" << std::endl;
       }
    else std::cout << __FUNCTION__ << " WARNING: Cahnnel ID not supported, 1 - 64 channels allowed." << std::endl;
}
//---------------------------------------------------------------------------------------------------------------
void LGADBase::SetDUTBoards(std::vector<AqBoard> Brds)
{
    for (unsigned int k = 0; k < Brds.size(); k++)
        {
         if (Brds.at(k) == SingleCh || Brds.at(k) == FourCh || Brds.at(k) == IN2P3 || Brds.at(k) == KU || Brds.at(k) == SiPM) m_DUTChsBrd.push_back(std::make_pair(-99, Brds.at(k)));
         else std::cout << __FUNCTION__ << " WARNING: Board type on channel board type vector entry " << k << " not supported, ignoring..." << std::endl;
        }
}
//---------------------------------------------------------------------------------------------------------------
void LGADBase::SetDUTransImp(int ChId, float Transimp)
{
    if (0 < ChId && ChId > 64)
       {
        if (Transimp > 0) m_DUTChsTrns.push_back(std::make_pair(ChId, Transimp));
        else std::cout << __FUNCTION__ << " WARNING: Invalid transimpedence value for channel " << ChId << "!" << std::endl;
       }
    else std::cout << __FUNCTION__ << " WARNING: Cahnnel ID not supported, 1 - 64 channels allowed." << std::endl;
}
//---------------------------------------------------------------------------------------------------------------
void LGADBase::SetDUTransImps(std::vector<float> Transimps)
{
    for (unsigned k = 0; k < Transimps.size(); k++)
        {
         if (Transimps.at(k) > 0) m_DUTChsTrns.push_back(std::make_pair(-99, Transimps.at(k)));
         else std::cout << __FUNCTION__ << " WARNING: Invalid transimpedence value on channel transimpedemce vector entry " << k << ", ignoring..." << std::endl;
        }
}
//---------------------------------------------------------------------------------------------------------------
void LGADBase::SetDUTSecStage(int ChId, SecStage Amp)
{
    if (0 < ChId && ChId > 64)
       {
        if (Amp == MinCircuits || Amp == Particulars || Amp == CIVIDEC || Amp == none) m_DUTChsAmp.push_back(std::make_pair(ChId, Amp));
        else std::cout << __FUNCTION__ << " WARNING: Amplifier type for channel " << ChId << " not supported!" << std::endl;
       }
    else std::cout << __FUNCTION__ << " WARNING: Cahnnel ID not supported, 1 - 64 channels allowed." << std::endl;
}
//---------------------------------------------------------------------------------------------------------------
void LGADBase::SetDUTSecStages(std::vector<SecStage> Amps)
{
    for (unsigned int k = 0; k < Amps.size(); k++)
        {
         if (Amps.at(k) == MinCircuits || Amps.at(k) == Particulars || Amps.at(k) == CIVIDEC || Amps.at(k) == none) m_DUTChsAmp.push_back(std::make_pair(-99, Amps.at(k)));
         else std::cout << __FUNCTION__ << " WARNING: Invalid amplifier type on channel amplifier type vector entry " << k << ", ignoring..." << std::endl;
        }
}
//---------------------------------------------------------------------------------------------------------------
void LGADBase::SetDUTAmpGain(int ChId, int gain)
{
    if (0 < ChId && ChId > 64)
       {
        if (gain > 1) m_DUTChsAmpGn.push_back(std::make_pair(ChId, gain));
        else std::cout << __FUNCTION__ << " WARNING: Invalid amplifier gain value for channel " << ChId << "!" << std::endl;
       }
    else std::cout << __FUNCTION__ << " WARNING: Cahnnel ID not supported, 1 - 64 channels allowed." << std::endl;
}
//---------------------------------------------------------------------------------------------------------------
void LGADBase::SetDUTAmpGains(std::vector<int> Gains)
{
    for (unsigned int k = 0; k < Gains.size(); k++)
        {
         if (Gains.at(k) > 1) m_DUTChsAmpGn.push_back(std::make_pair(-99, Gains.at(k)));
         else std::cout << __FUNCTION__ << " WARNING: Invalid amplifier gain value on channel gain vector entry " << k << ", ignoring..." << std::endl;
        }
}
//---------------------------------------------------------------------------------------------------------------
void LGADBase::SetDUTrigg(int ChId, float trigg)
{
    if (0 < ChId && ChId > 64)
       {
        if (trigg != 0) m_DUTChsTrigg.push_back(std::make_pair(ChId, trigg));
        else std::cout << __FUNCTION__ << " WARNING: Invalid trigger voltage value for channel " << ChId << "!" << std::endl;
       }
    else std::cout << __FUNCTION__ << " WARNING: Cahnnel ID not supported, 1 - 64 channels allowed." << std::endl;
}
//---------------------------------------------------------------------------------------------------------------
void LGADBase::SetDUTriggs(std::vector<float> Triggs)
{
    for (unsigned int k = 0; k < Triggs.size(); k++)
        {
         if (Triggs.at(k) != 0) m_DUTChsTrigg.push_back(std::make_pair(-99, Triggs.at(k)));
         else std::cout << __FUNCTION__ << " WARNING: Invalid trigger voltage value on channel trigger vector entry " << k << ", ignoring..." << std::endl;
        }
}
//---------------------------------------------------------------------------------------------------------------
void LGADBase::SetDUTFract(int ChId, float frac)
{
    if (0 < ChId && ChId > 64)
       {
        if (frac >= 0.5 && frac <=0.95) m_DUTChsFrc.push_back(std::make_pair(ChId, frac));
        else std::cout << __FUNCTION__ << " WARNING: Invalid CFD fraction value for channel " << ChId << "!" << std::endl;
       }
    else std::cout << __FUNCTION__ << " WARNING: Cahnnel ID not supported, 1 - 64 channels allowed." << std::endl;
}
//---------------------------------------------------------------------------------------------------------------
void LGADBase::SetDUTFracts(std::vector<float> Fracs)
{
    for (unsigned int k = 0; k < Fracs.size(); k++)
        {
         if (Fracs.at(k) >= 0.5 && Fracs.at(k) <=0.95) m_DUTChsFrc.push_back(std::make_pair(-99, Fracs.at(k)));
         else std::cout << __FUNCTION__ << " WARNING: Invalid CFD fraction value on channel CFD vector entry " << k << ", ignoring..." << std::endl;
        }
}
//---------------------------------------------------------------------------------------------------------------
void LGADBase::SetDUTCap(int ChId, double cap)
{
    if (0 < ChId && ChId > 64)
       {
        if (cap > 0) m_DUTChsCaps.push_back(std::make_pair(ChId, cap));
        else std::cout << __FUNCTION__ << " WARNING: Invalid capacitance value for channel " << ChId << "!" << std::endl;
       }
    else std::cout << __FUNCTION__ << " WARNING: Cahnnel ID not supported, 1 - 64 channels allowed." << std::endl;
}
//---------------------------------------------------------------------------------------------------------------
void LGADBase::SetDUTCaps(std::vector<double> Caps)
{
    for (unsigned int k = 0; k < Caps.size(); k++)
        {
         if (Caps.at(k) > 0) m_DUTChsCaps.push_back(std::make_pair(-99, Caps.at(k)));
         else std::cout << __FUNCTION__ << " WARNING: Invalid capacitance value on channel capacitance vector entry " << k << ", ignoring..." << std::endl;
        }
}
//---------------------------------------------------------------------------------------------------------------
void LGADBase::SetChDMgt(int ChId1, int ChId2, double DMgt, int Qt)
{
    if (0 < ChId1 && ChId1 > 64)
       {
        if (ChId2 > ChId1 && ChId2 > 64)
           {
            if (DMgt > 0) 
               {
                if (Qt == 0) m_ChsDTs.push_back(std::make_pair(std::make_pair(ChId1, ChId2), DMgt));
                else if (Qt == 1) m_ChsDCs.push_back(std::make_pair(std::make_pair(ChId1, ChId2), DMgt));
               }
            else {
                  std::cout << __FUNCTION__ << " WARNING: Invalid"; 
                  if (Qt == 0) std::cout << " time difference";
                  else if (Qt ==1) std::cout << " charge difference";
                  std::cout << " value between channels " << ChId1 << " and " << ChId2 << "!" << std::endl;
                 }
           }
        else std::cout << __FUNCTION__ << " WARNING: Cahnnel ID2 needs to be grater thand ChId1 and between 2-64!" << std::endl;
       }
    else std::cout << __FUNCTION__ << " WARNING: Cahnnel ID not supported, 1 - 64 channels allowed." << std::endl;
}
// --------------------------------------------------------------------------------------------------------------
void LGADBase::SetChDMgts(std::vector<double> DMgt, int Qt)
{
    for (unsigned int k = 0; k < DMgt.size(); k++)
        {
         if (DMgt.at(k) > 0) 
            {
             if (Qt == 0) m_ChsDTs.push_back(std::make_pair(std::make_pair(-99, -99), DMgt.at(k)));
             if (Qt == 1) m_ChsDCs.push_back(std::make_pair(std::make_pair(-99, -99), DMgt.at(k)));
            }
         else {
               if (Qt == 0) std::cout << __FUNCTION__ << " WARNING: Invalid time difference value on channel dT";
               if (Qt == 1) std::cout << __FUNCTION__ << " WARNING: Invalid cahrge difference value on channel dC";
               std::cout <<" vector entry " << k << ", ignoring..." << std::endl;
              }
        }
}
//---------------------------------------------------------------------------------------------------------------
bool LGADBase::OpenTransFile(TString filename)
{
    if (filename == "" || ((trim((std::string)filename, " ")).empty()))
       { 
        std::cout << __FUNCTION__ << " WARNING: No transimpedence file set, ignoring transimpedence correction..." << std::endl;
        m_TrnsCorr = false; 
        return false; 
       }
    else {
          TFile *m_TransFile = TFile::Open(filename);
          if (m_TransFile->IsOpen()) 
             {
              // Clone histogeams
              unsigned int n = 0;
              TIter nextkey(m_TransFile->GetListOfKeys());
              TKey *key = 0;
              TrCrHist info;
              TObject *obj = key->ReadObj();
              while ((key = (TKey*)nextkey()))
                    {
                     if (obj->IsA()->InheritsFrom(TH2D::Class()))
                        {
                         if (!(strcmp(obj->GetName(), "hs1_3p")) || !(strcmp(obj->GetName(), "hs1_4p")) || !(strcmp(obj->GetName(), "hs2_3p")) || !(strcmp(obj->GetName(), "hs2_4p"))) 
                            {
                             info.TransHist = (TH2D*)obj->Clone();
                             info.TransHist->SetDirectory(0);
                             if (!(strcmp(obj->GetName(), "hs1_3p")) || !(strcmp(obj->GetName(), "hs1_4p"))) info.Board = SingleCh;
                             else if (!(strcmp(obj->GetName(), "hs2_3p")) || !(strcmp(obj->GetName(), "hs2_4p"))) info.Board = FourCh;
                             if (!(strcmp(obj->GetName(), "hs1_3p")) || !(strcmp(obj->GetName(), "hs2_3p"))) info.Capacitance = 3;
                             else if (!(strcmp(obj->GetName(), "hs1_4p")) || !(strcmp(obj->GetName(), "hs2_4p"))) info.Capacitance = 4;
                             m_TrsHists.push_back(info);
                             n++; 
                            }
                        }
                   }
              if (n > 0) return true;
              else {
                    std::cout << __FUNCTION__ << " WARNING: No transimpedence histograms found in " << filename << std::endl;
                    m_TrnsCorr = false;
                    return false;
                   }
             }
          else {
                std::cout << __FUNCTION__ << " WARNING: Failed openning Tansimpedence file: " << filename << std::endl;
                m_TrnsCorr = false;
                return false;
               }
         }
}
// --------------------------------------------------------------------------------------------------------------
bool LGADBase::GetTrackComb()
{
    return m_TrackComb;
}
// --------------------------------------------------------------------------------------------------------------
void LGADBase::SetFEi4Eff(bool FEi4Eff)
{
    m_fei4Eff = FEi4Eff;
}
// --------------------------------------------------------------------------------------------------------------
bool LGADBase::GetFEi4Eff()
{
    return m_fei4Eff;
}
// --------------------------------------------------------------------------------------------------------------
bool LGADBase::SetSRate(Long64_t rate, unsigned int ch)
{
    if (ch > m_srate.size())
       { 
        std::cout << __FUNCTION__ << " ERROR: Invalid chanel number" << std::endl; 
        return false;
       }
    m_srate[ch] = rate;
    return true;
}
// --------------------------------------------------------------------------------------------------------------
Long64_t LGADBase::GetSRate(unsigned int ch)
{
    return m_srate.at(ch);
}
// --------------------------------------------------------------------------------------------------------------
bool LGADBase::SetNPoints(unsigned int points, unsigned int ch)
{
  if (ch > m_npoints.size())
     { 
      std::cout << __FUNCTION__ << " ERROR: Invalid chanel number" << std::endl;
      return false;
     }
  m_npoints[ch] = points;
  return true;
}
// --------------------------------------------------------------------------------------------------------------
unsigned int LGADBase::GetNPoints(unsigned int ch)
{
    return m_npoints.at(ch);
}
// --------------------------------------------------------------------------------------------------------------
void LGADBase::Initialize()
{
    m_verbose = 1;
    m_channels.clear();
    m_channels.push_back(0);
    m_scale.clear();
    m_scale.push_back(0);
    m_physt.clear();
    m_physt.push_back(0);
    m_ordrt.clear();
    m_ordrt.push_back(0);
    m_npoints.clear();
    m_srate.clear();
    m_DUTChsBrd.clear();
    m_DUTChsAmp.clear();
    m_DUTChsFrc.clear();
    m_DUTChsTrigg.clear();
    m_DUTChsTrns.clear();
    m_DUTChsAmpGn.clear();
    m_DUTChsNames.clear();
    m_DUTChsCaps.clear();
    m_ChsDTs.clear();
    m_ChsDCs.clear();
    if (m_instrument == TektronixScope)
       {
        m_npoints.push_back(1024);
        m_srate.push_back(5e9);
       }
    else if (m_instrument == Sampic)
            {
             m_npoints.push_back(64);
             m_srate.push_back(8512e6);
            }
    else if (m_instrument == LabTXT || m_instrument == Unasigned)
            {
             m_npoints.push_back(1024);
             m_srate.push_back(2e10);
            }
    else if (m_instrument == InfiniiumScope || m_instrument == LeCroyWRBin)
            {
             m_npoints.push_back(1024);
             m_srate.push_back(4e10);
            }

    m_evnt1 = 0;
    m_evnt2 = 0;
    m_TrackComb = false;
    m_fei4Eff = false;
    m_convert = false;
    m_WaveShape = false;
    m_TrnsCorr = false;
    m_ext = "";
    m_datadir = "";
    m_dataname = "";
    m_treename = "";
    m_ofdir = "";
    m_ofname = "";
    if (m_TrnsCorr) OpenTransFile(m_TransFileName);
    if (m_verbose == 2) std::cout << __FUNCTION__ << " INFO: Initializing sampling rate to default value : " << m_srate.at(0) / 1e6 << "MS/sec and no. of points to : " << m_npoints.at(0) << std::endl;
}
// --------------------------------------------------------------------------------------------------------------
void LGADBase::SetStartStopEvnt(int Evnt1, int Evnt2)
{
    if (Evnt1 >= 0) m_evnt1 = Evnt1;
    else std::cout << __FUNCTION__ << " WARNING: Start event value not coherent, setting to zero!" << std::endl;
    if (Evnt2 > Evnt1 && Evnt2 > 0) m_evnt2 = Evnt2;
    else if (Evnt2 != 0) { m_evnt2 = 0; std::cout << __FUNCTION__ << " WARNING: Stop event value not coherent, setting to zero!" << std::endl;}
    else m_evnt2 = 0;
}
// --------------------------------------------------------------------------------------------------------------
std::pair <unsigned int, unsigned int> LGADBase::GetStartStopEvnt()
{
    return  std::make_pair(m_evnt1, m_evnt2);
}
// --------------------------------------------------------------------------------------------------------------
void LGADBase::SetVectorSize(unsigned int nch)
{
    m_physt.clear();
    m_ordrt.clear();
    m_t.clear();
    m_w.clear();
    m_scale.clear();
    m_triggTime.clear();
    m_npoints.clear();
    m_srate.clear();

    m_t.resize(nch);
    m_w.resize(nch);

    if (m_instrument == TektronixScope || m_instrument == Unasigned) m_triggTime.resize(nch);
    if (m_instrument == Sampic)
       {
        m_physt.resize(nch);
        m_ordrt.resize(nch);
        m_npoints.resize(1);
        m_srate.resize(1);
       }
    else {
          m_npoints.resize(nch);
          m_srate.resize(nch);
          m_scale.resize(nch);
         }

}
// --------------------------------------------------------------------------------------------------------------
int LGADBase::Addoriel(int val)
{
    int fract = 0;
    if (val > 0) for (int i = 0; i <= val; i++) fract += i;
    else {
          fract = -1;
          std::cout << __FUNCTION__ << " ERROR: Trying to calculate factorial of negative number!" << std::endl;
         }
    return fract;
}
// --------------------------------------------------------------------------------------------------------------
bool LGADBase::ProgressBar(Long64_t evt, Long64_t total)
{
    float prd = (float)total / (float)100;
    if (((remainder((float)evt, prd) < 0 && ceil(remainder((float)evt, prd)) == 0)) && (evt + 1) < total)
       {
        std::cout << "<" << std::setfill('=') << std::setw(floor((0.2*((float)evt/prd)))) << "";
        std::cout << std::setfill(' ') << std::setw(20 - floor((0.2*((float)evt / prd))) + 2) << std::right << "> :";
        std::cout << std::left << round((float)evt / prd) << "%" << ", Processed entries: " << evt + 1 << " / " << total << "\r" << std::setfill(' ');
        return false;
       }
    else if ((evt + 1) == total) 
            {
             std::cout << "<" << std::setfill('=') << std::setw(19) << "=" << std::setfill(' ') 
                       << "> :" << 100 << "%" << ", Processed entries: "  << evt+1 << " / " 
                       << total << "\r" << std::setfill(' ') << std::endl;
             return true;
            }
    return false;
}
// --------------------------------------------------------------------------------------------------------------
bool LGADBase::SetRootTree(TFile* f, std::string name)
{
int n = 0;
TTree* data_tree = NULL;
TIter nextkey(f->GetListOfKeys());
TKey *key = 0;

while ((key = (TKey*)nextkey())) 
      {
       TObject *obj = key->ReadObj();
       if (obj->IsA()->InheritsFrom(TTree::Class())) 
          {
           if (name != "" && obj->GetName() == name)
              { 
               data_tree = (TTree*)obj;
               n = 1;
               break;
              }
           else {
                 data_tree = (TTree*)obj;
                 n++;
                }
          }
      }

if (n == 1) { m_tree = data_tree; return true; }
else if (n == 0) return false;
else return false;
}
// --------------------------------------------------------------------------------------------------------------
// Helper Function integrated to trim strings by reoving spaces or other defined characters from it's begining and 
// end. 
std::string LGADBase::trim(const std::string str, const std::string whitespace)
{
    int strBegin = str.find_first_not_of(whitespace);
    if ((int)std::string::npos == strBegin) return ""; // no content
    unsigned int strEnd = str.find_last_not_of(whitespace);
    unsigned int strRange = strEnd - strBegin + 1;
    return str.substr(strBegin, strRange);
}
// --------------------------------------------------------------------------------------------------------------
// Function to remove characters and/or spaces from the totality, beginning and end of the string.
std::string LGADBase::reduce(const std::string str, const std::string fill, const std::string whitespace)
{
    std::string result = LGADBase::trim(str, whitespace); // trim first
                                                // Replace sub ranges
    int beginSpace = result.find_first_of(whitespace);
    while ((int)std::string::npos != beginSpace)
          {  
           unsigned int endSpace = result.find_first_not_of(whitespace, beginSpace);
           unsigned int range = endSpace - beginSpace;
           result.replace(beginSpace, range, fill);
           unsigned int newStart = beginSpace + fill.length();
           beginSpace = result.find_first_of(whitespace, newStart);
          }

    return result;
}
// --------------------------------------------------------------------------------------------------------------
// Function to create directories in linux and windows
int LGADBase::RecursMkDir(const char* dirname)
{
    int check = 0;
    const size_t len = strlen(dirname);
#ifdef _WIN32
    char _path[MAX_PATH];
#else
    char _path[PATH_MAX];
#endif
    char *p;

    /* Copy string so its mutable */
    if (len > sizeof(_path)-1) 
       {
        std::cout << __FUNCTION__ << " ERROR: Path too long!" << std::endl;
        return -1;
       }
    strcpy(_path, dirname);

    /* Iterate the string */
    unsigned int t = 0;
    for (p = _path + 1; *p; p++) 
        {
         if (*p == '/') 
            {
             t++;
             if (t == 1) continue;
             *p = '\0'; // Temporarily truncate 
             if (LGADBase::CreateDir(_path) != 1) return -1;
             *p = '/';
            }
        }

    if (LGADBase::CreateDir(_path) != 1) return -1;

    return 1;
}
// --------------------------------------------------------------------------------------------------------------
int LGADBase::CreateDir(const char* path)
{
    int check = 0;

    if (DirExists(path) != 1)
       {
#ifdef _WIN32
       check = 0;
       check = CreateDirectory(path, NULL);
#else
       check = mkdir(path, 0777);
       if (check == -1) check = 0;
       else check = 1;
#endif
       if (check == 1 && DirExists(path) != 1)
          {
           std::cout << __FUNCTION__ << " ERROR: Folder creation failed: " << path << "!" << std::endl;
           return -1;
          }
       }
    else check = 1;

    return check;
}
// --------------------------------------------------------------------------------------------------------------
// Function to check if directory exists
int LGADBase::DirExists(const char* path)
{
    struct stat info;
    if (stat(path, &info) != 0) return 0;
    else if (info.st_mode & S_IFDIR) return 1;
    else return -1;
}
// --------------------------------------------------------------------------------------------------------------
// Function to return vector of filenames in a directory                                                                                                                        
std::vector<std::string> LGADBase::ListFileNames(const char* path, const char* ext)
{
    std::vector<std::string> filenames;
    std::string search_path = path;
    std::string search_ext = ext;
#ifdef _WIN32
    search_path += "*.";
    search_path += search_ext;
    WIN32_FIND_DATA fd;
    HANDLE hFind = ::FindFirstFile(search_path.c_str(), &fd);
    if (hFind != INVALID_HANDLE_VALUE)
       {
        do {
            if (!(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) filenames.push_back(fd.cFileName);
           } 
        while (::FindNextFile(hFind, &fd));
        ::FindClose(hFind);
       }
#else
    struct dirent *pdir;
    DIR *path_dir = opendir(path);
    while ((pdir = readdir(path_dir)))
          {
           search_path = pdir->d_name;
           if (search_path.find(search_ext) != std::string::npos) filenames.push_back(pdir->d_name);
          }
#endif

    return filenames;
}
// --------------------------------------------------------------------------------------------------------------
bool LGADBase::SetScale(std::vector<unsigned int> channel, unsigned int nchan, std::vector<float>* scale)
{
    // Input scales for the diferent channels
    for (unsigned int c = 0; c < nchan; c++)
        { 
         scale->at(c) = -1;
         scaleset: std::cout << "Please enter vertical scale for channel " << channel.at(c) << " in mV/DiV : ";
         std::cin >> scale->at(c);
         if (scale->at(c) != 5 && scale->at(c) != 10 && scale->at(c) != 20 && scale->at(c) != 50 && scale->at(c) != 100
             && scale->at(c) != 200 && scale->at(c) != 300 && scale->at(c) != 500 && scale->at(c) != 1000)
            {
             std::cout << __FUNCTION__ << " ERROR: Incorect scale, please try again!" << std::endl;
             goto scaleset;
            }
        }
    return true;
}
// --------------------------------------------------------------------------------------------------------------
unsigned int LGADBase::CountFiles(const char* dir, const char* ext)
{
    // Count the number of files in the idrectory
    unsigned int nfiles = 0;
    TSystemDirectory *directory = new TSystemDirectory("", dir);
    TList *files = directory->GetListOfFiles();
    if (files)
       {
        TSystemFile *sfile;
        TIter next(files);
        TString fname;
        while ((sfile=(TSystemFile*)next()))
              {
               fname = sfile->GetName();
               if (!sfile->IsDirectory())
                  {
                   if (fname.EndsWith(ext)) nfiles++;
                   else if (m_instrument == Sampic && (TString(fname(0, fname.Length()-5))).EndsWith(ext)) nfiles++;
                  }
              }
       }
    return nfiles;
}
// --------------------------------------------------------------------------------------------------------------
template <typename T> T LGADBase::Derivate(T *w, int start)
{
    int sz = w->size();
    T deriv;
    if (start <= 0) 
       {
        std::cout << __FUNCTION__ << " WARNING: Start index invalid, starting from element 0!" << std::endl;
        start = 1;
       }
    if (sz > start)
       {
        deriv.reserve(sz - start);
        for (int k = start-1; k < (sz-1); k++) deriv.push_back(w->at(k+1) - w->at(k));
       }
    else std::cout << __FUNCTION__ << " ERROR: Vector has too few elements to derivate!" << std::endl;
    return deriv;
}
// --------------------------------------------------------------------------------------------------------------
template <typename U> double LGADBase::Mean(U *w, int start, int stop)
{
    double mean = -99.;
    int sz = 0;
    if (start <= 0 && stop <= 0) sz = w->size();
    else if (start <= 0 && stop > 0) sz = stop;
    else if (stop <= 0 && start > 0) sz = w->size() - start;
    else sz = stop - start;
    if (sz > 0)
       { 
        if (start <= 0 && stop <= 0) mean = std::accumulate(w->begin(), w->end(), 0.0) / sz;
        else if (start <= 0 && stop > 0) mean = std::accumulate(w->begin(), w->begin() + stop, 0.0) / sz;
        else if (stop <= 0 && start > 0) mean = std::accumulate(w->begin() + start, w->end(), 0.0) / sz;
        else mean = std::accumulate(w->begin() + start, w->begin() + stop, 0.0) / sz;
       }
    else std::cout << __FUNCTION__ << " ERROR: Trying to calculate mean of empty vector or inproper range!" << std::endl;
    return mean;
}
// --------------------------------------------------------------------------------------------------------------
template <typename V> double LGADBase::Stdev(V *w, int start, int stop, double mean)
{
    double stdev = -99.;
    double E = 0;
    int sz = 0;
    if (start <= 0 && stop <= 0) sz = w->size();
    else if (start <= 0 && stop > 0) sz = stop;
    else if (stop <= 0 && start > 0) sz = w->size() - start;
    else sz = stop - start;
    if (sz > 0)
       {
        double ave = 0.0;
        if (mean == 0.0) ave = Mean(w, start, stop);
        else ave = mean;
        if ((start <= 0 && stop <= 0) || (start <= 0 && stop > 0)) for (int i = 0; i < sz; i++) E += pow((ave - w->at(i)), 2);
        else if (stop <= 0 && start > 0) for (unsigned int i = start; i < w->size(); i++) E += pow((ave - w->at(i)), 2);
        else for (int i = start; i < stop; i++) E += pow((ave - w->at(i)), 2);
        stdev = sqrt(E / (sz-1));
       }
    else std::cout << __FUNCTION__ << " ERROR: Trying to calculate standard deviation of empty vector or improper range!" << std::endl;
    return stdev;
}
// --------------------------------------------------------------------------------------------------------------
template <typename V, typename T> double LGADBase::BayesianErr(V *w, T value)
{
    double events = w->size();
    double pass = 0;
    for (unsigned int k = 0; k < events; k++) if ((w->at(k)) >= value) pass++;
    if (events > 0) return sqrt((((pass + 1)*(pass + 2)) / ((events + 2)*(events + 3))) - (pow(pass + 1, 2) / pow(events + 2, 2)));
    else {
          std::cout << __FUNCTION__ << " ERROR: Trying to calculate baysian uncertenty of empty distribution!" << std::endl;
          return 99.;
         }
}
// --------------------------------------------------------------------------------------------------------------
template <typename V> double LGADBase::CalcMeadian(V *vec, int start, int stop)
{
    if (stop <= 0 || stop > (int)(vec->size()) || stop <= start) stop = vec->size();
    if (start < 0 || start >= (int)(vec->size())) start = 0;

    V wmod;
    for (int ga = start; ga < stop; ga++) wmod.push_back(vec->at(ga));
    size_t size = wmod.size();
    if (size == 0)
       {
        std::cout << __FUNCTION__ << " ERROR: Trying to calculate median of empty vector!" << std::endl;
        return -1;  // Undefined, realy.
       }
    else {
          sort(wmod.begin(), wmod.end());
          if (size % 2 == 0) return (wmod.at((size / 2) - 1) + wmod.at(size / 2)) / 2;
          else return wmod.at(ceil((float)size / 2));
         }
}
// --------------------------------------------------------------------------------------------------------------
template <typename V> V LGADBase::OutlierReject(V *w, unsigned int order, float elem, int start, int stop)
{
    if (stop <= 0 || stop > (int)(w->size()) || stop <= start) stop = w->size();
    if (start < 0 || start >= (int)(w->size())) start = 0;

    V wmod;
    for (int ga = start; ga < stop; ga++) wmod.push_back(w->at(ga));
    sort(wmod.begin(), wmod.end());
    double mnt = 0.0;
    int rd = ceil(((float)(wmod.size()*elem))/2);
    // if (m_verbose == 2) std::cout << __FUNCTION__ << ": " << wmod.size() << " " << elem << " " << rd << std::endl; 
    for (unsigned int i = 0; i < rd; i++) mnt += (wmod.at(ceil(((float)wmod.size()/2)+i)) + wmod.at(floor(((float)wmod.size()/2)-i)))/2;
    mnt = mnt/((float)rd);
    // if (m_verbose == 2) std::cout << __FUNCTION__ << ": " << mnt << std::endl;
    for (unsigned int de = 0; de < wmod.size(); de++)
        {
         if (wmod.at(de) !=0)
            {
             if (fabs(log10(wmod.at(de)/mnt)) >= order)
                {
                 wmod.erase(wmod.begin() + de);
                 de--;
                }
            }
        }

    return wmod;
}
// --------------------------------------------------------------------------------------------------------------
template <typename T, typename V> T LGADBase::CalResolution(V *w, unsigned int order, int start, int stop)
{ 
    if (stop <= 0 || stop > (int)(w->size()) || stop <= start) stop = w->size();
    if (start < 0 || start >= (int)(w->size())) start = 0;

    T res = -1;
    V mag;
    mag.reserve(LGADBase::Addoriel(stop-start-1));
    double l = 0;
    for (int ga = start; ga < stop - 1; ga++)
        {
         for (int ma = ga + 1; ma < stop; ma++)
             {
              l = fabs((double)(w->at(ma) - w->at(ga)));
              if (l > 0) mag.push_back(l);
             }
        }
    // Remove extreme elements from vectors to account for accuracy issues
    l = LGADBase::Mean(&mag);
    for (unsigned int kr = 0; kr < mag.size(); kr++)
            {
             if (fabs(log10(mag.at(kr)/l)) >= order)
                {
                 mag.erase(mag.begin() + kr);
                 kr--;
                }
            }
    res = *std::min_element(mag.begin(), mag.end());

    return res;
}
// --------------------------------------------------------------------------------------------------------------
// Explicit template instantanization
template vector<int> LGADBase::Derivate<vector<int> >(std::vector<int> *, int);
template vector<float> LGADBase::Derivate<vector<float> >(std::vector<float> *, int);
template vector<double> LGADBase::Derivate<vector<double> >(std::vector<double> *, int);
template double LGADBase::Mean<vector<int> >(std::vector<int> *, int, int);
template double LGADBase::Mean<vector<bool> >(std::vector<bool> *, int, int);
template double LGADBase::Mean<vector<float> >(std::vector<float> *, int, int);
template double LGADBase::Mean<vector<double> >(std::vector<double> *, int, int);
template double LGADBase::Stdev<vector<int> >(std::vector<int> *, int, int, double);
template double LGADBase::Stdev<vector<bool> >(std::vector<bool> *, int, int, double);
template double LGADBase::Stdev<vector<float> >(std::vector<float> *, int, int, double);
template double LGADBase::Stdev<vector<double> >(std::vector<double> *, int, int, double);
template double LGADBase::BayesianErr<vector<int>, int>(std::vector<int> *, int);
template double LGADBase::BayesianErr<vector<bool>, bool>(std::vector<bool> *, bool);
template double LGADBase::CalcMeadian<vector<int> >(std::vector<int> *, int, int);
template double LGADBase::CalcMeadian<vector<float> >(std::vector<float> *, int, int);
template double LGADBase::CalcMeadian<vector<double> >(std::vector<double> *, int, int);
template vector<int> LGADBase::OutlierReject<vector<int> >(std::vector<int> *, unsigned int, float, int, int);
template vector<float> LGADBase::OutlierReject<vector<float> >(std::vector<float> *, unsigned int, float, int, int);
template vector<double> LGADBase::OutlierReject<vector<double> >(std::vector<double> *, unsigned int, float, int, int);
template int LGADBase::CalResolution<int, vector<int> >(std::vector<int> *, unsigned int, int, int);
template float LGADBase::CalResolution<float, vector<float> >(std::vector<float> *, unsigned int, int, int);
template double LGADBase::CalResolution<double, vector<double> >(std::vector<double> *, unsigned int, int, int);
