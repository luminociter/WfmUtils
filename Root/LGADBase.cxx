/*
* LGADBase.cxx
*
*
*      Author: Gkougkousis Evangelos - Leonidas
*              egkougko@cern.ch
*                CERN-GENEVA
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
    m_TrnsCorr = false;
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
    if (instr == Sampic || instr == LabTXT || instr == TestBeamBin || instr == TektronixScope || instr == LeCroyWRBin || instr == Unasigned || instr == TestBeamBin1 || instr == TestBeamBin2) m_instrument = instr;
    else {
          std::cout << __FUNCTION__ << " ERROR: instrument value " << instr
                    << " not allowed. Allowed values are in ascending order: Sampic, LabTXT, TestBeamBin, TektronixScope, LeCroyWRBin, or Unasigned. Falling back to "
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
void LGADBase::SetChMag(int ChId, double ChMag, int Qt)
{
    if (0 < ChId && ChId < 65)
       {
        if (ChMag > 0)
           {
            if (Qt == 0) m_ChVoltCuts.push_back(std::make_pair(ChId, ChMag));
            else if (Qt == 1) m_ChNoiseCuts.push_back(std::make_pair(ChId, ChMag));
            else if (Qt == 2) m_ChChargeCuts.push_back(std::make_pair(ChId, ChMag));
            else if (Qt == 3) m_ChJitterCuts.push_back(std::make_pair(ChId, ChMag));
           }
        else {
              std::cout << __FUNCTION__ << " WARNING: Invalid"; 
              if (Qt == 0) std::cout << " voltage cut";
              else if (Qt == 1) std::cout << " noise cut";
              else if (Qt == 2) std::cout << " charge cut";
              else if (Qt == 3) std::cout << " jitter cut";
              std::cout << " value for channel " << ChId << "!" << std::endl;
             }
       }
    else std::cout << __FUNCTION__ << " WARNING: Plane id not supported." << std::endl; 
}
//---------------------------------------------------------------------------------------------------------------
void LGADBase::SetChMags(std::vector<double> ChMags, int Qt)
{
    for (unsigned int k = 0; k < ChMags.size(); k++)
        {
         if (ChMags.at(k) > 0)
            {
             if (Qt == 0) m_ChVoltCuts.push_back(std::make_pair(-99, ChMags.at(k)));
             else if (Qt == 1) m_ChNoiseCuts.push_back(std::make_pair( -99, ChMags.at(k)));
             else if (Qt == 2) m_ChChargeCuts.push_back(std::make_pair(-99, ChMags.at(k)));
             else if (Qt == 3) m_ChJitterCuts.push_back(std::make_pair(-99, ChMags.at(k)));
            }
         else {
               if (Qt == 0) std::cout << __FUNCTION__ << " WARNING: Invalid value on chanel voltage cut";
               else if (Qt == 1) std::cout << __FUNCTION__ << " WARNING: Invalid value on chanel noise cut";
               else if (Qt == 2) std::cout << __FUNCTION__ << " WARNING: Invalid value on chanel charge cut";
               else if (Qt == 3) std::cout << __FUNCTION__ << " WARNING: Invalid value on chanel jitter cut";
               std::cout <<" vector entry " << k << ", ignoring..." << std::endl;
              }
        }
}
//---------------------------------------------------------------------------------------------------------------
void LGADBase::SetPlaneDMgt(int PlaneId1, int PlaneId2, double DMgt, int Qt)
{
    if (0 < PlaneId1)
       {
        if (PlaneId2 > PlaneId1)
           {
            if (DMgt > 0) 
               {
                if (Qt == 0) m_PlaneDTs.push_back(std::make_pair(std::make_pair(PlaneId1, PlaneId2), DMgt));
                else if (Qt == 1) m_PlaneDCs.push_back(std::make_pair(std::make_pair(PlaneId1, PlaneId2), DMgt));
               }
            else {
                  std::cout << __FUNCTION__ << " WARNING: Invalid"; 
                  if (Qt == 0) std::cout << " time difference";
                  else if (Qt ==1) std::cout << " charge difference";
                  std::cout << " value between channels " << PlaneId1 << " and " << PlaneId2 << "!" << std::endl;
                 }
           }
        else std::cout << __FUNCTION__ << " WARNING: Plane id 2 needs to be grater than Plane Id 1!" << std::endl;
       }
    else std::cout << __FUNCTION__ << " WARNING: Plane id not supported." << std::endl;
}
// --------------------------------------------------------------------------------------------------------------
void LGADBase::SetPlaneDMgts(std::vector<double> DMgts, int Qt)
{
    for (unsigned int k = 0; k < DMgts.size(); k++)
        {
         if (DMgts.at(k) > 0) 
            {
             if (Qt == 0) m_PlaneDTs.push_back(std::make_pair(std::make_pair(-99, -99), DMgts.at(k)));
             else if (Qt == 1) m_PlaneDCs.push_back(std::make_pair(std::make_pair(-99, -99), DMgts.at(k)));
            }
         else {
               if (Qt == 0) std::cout << __FUNCTION__ << " WARNING: Invalid time difference value on plane time";
               else if (Qt == 1) std::cout << __FUNCTION__ << " WARNING: Invalid cahrge difference value on plane charge";
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
    m_PlaneDTs.clear();
    m_PlaneDCs.clear();
    m_ChNoiseCuts.clear();
    m_ChVoltCuts.clear();
    m_ChChargeCuts.clear();
    m_ChJitterCuts.clear();
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
    else if (m_instrument == TestBeamBin1 || TestBeamBin2 || TestBeamBin || m_instrument == LeCroyWRBin)
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
    float prd = (float)total/(float)100;
    if ((remainder((float)evt, prd) <= 0  && ceil(remainder((float)evt, prd)) == 0) && (evt + 1) < total)
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
bool LGADBase::PrintFitInfo(TH1D* histo, TCanvas** ca, std::string funcName)
{
  TIter nextfunc(histo->GetListOfFunctions());
  TFunction* func = 0;
  std::vector <std::string> FitNames;
  while ((func = (TFunction*)nextfunc())) FitNames.push_back(func->GetName());
  bool found = false;
  if (m_verbose > 2)
     {
      std::cout << __FUNCTION__ << " INFO: Found functions in histogram " << histo->GetName() << ": ";
      for (unsigned int k = 0; k < FitNames.size(); k++) std::cout << FitNames.at(k) << " ";
      std::cout << std::endl;
     }
  if (funcName == "none" && FitNames.size() == 1) { found = true; funcName = FitNames.at(0); }
  else if (funcName != "none" && FitNames.size() > 1)
          {
           for (unsigned int k = 0; k < FitNames.size(); k++)
               {
                if (FitNames.at(k) == funcName) { found = true; break; }
               }
          }
  if (found)
     {
      std::string fitname;
      std::string par0;
      std::string par1;
      if (funcName == "LandXGauFun") 
         {
          fitname = "Landau X Guass Fit";
          par0 = Form("MPV: %e +/- %e", histo->GetFunction(funcName.c_str())->GetParameter(1), histo->GetFunction(funcName.c_str())->GetParError(1));
          par1 = Form("#lambda : %e +/- %e", histo->GetFunction(funcName.c_str())->GetParameter(3), histo->GetFunction(funcName.c_str())->GetParError(3));
         }
      else if (funcName == "GaussFit") 
              {
               fitname = "Gauss Fit";
               par0 = Form("Mean: %e +/- %e", histo->GetFunction(funcName.c_str())->GetParameter(1), histo->GetFunction(funcName.c_str())->GetParError(1));
               par1 = Form("#sigma %e +/- %e", histo->GetFunction(funcName.c_str())->GetParameter(2), histo->GetFunction(funcName.c_str())->GetParError(2));
              }
      else if (funcName == "line")
              {
               fitname = "Linear Fit";
               par0 = Form("Intercept: %e +/- %e", histo->GetFunction(funcName.c_str())->GetParameter(0), histo->GetFunction(funcName.c_str())->GetParError(0));
               par1 = Form("Slope: %e +/- %e", histo->GetFunction(funcName.c_str())->GetParameter(1), histo->GetFunction(funcName.c_str())->GetParError(1));
              }
      else if (funcName == "LanXGau")
              {
               fitname = "RooFit Landau X Gauss Fit";
              }
      else if (funcName == "LinXGau")
              {
               fitname = "RooFit Linear X Gauss Fit";
              }
      else {
            std::cout << __FUNCTION__ << " WARNING: Unsupported fit type for histogram " << histo->GetName() << ": " << funcName << std::endl;
            return false;
           }    
      gROOT->SetBatch(kTRUE);
      TGaxis::SetMaxDigits(3);
      gStyle->SetOptStat(kFALSE);
      (*ca) = new TCanvas(Form("%s_Ca", histo->GetName()), histo->GetTitle(), 500, 500);
      //(*ca)->cd(1);
      gPad->SetGrid(1, 1);
      gPad->SetTicky(1);
      gPad->SetTickx(1);
      //gROOT->ForceStyle();
      //(*ca)->Draw();
      histo->Draw(/*"PE1"*/);
      histo->GetYaxis()->SetTitleOffset(1.3);
      TPaveText* paveA = new TPaveText(0.50, 0.60, 0.90, 0.75, "NDCNB");
      paveA->SetTextAlign(11);
      paveA->SetFillStyle(0);
      paveA->SetBorderSize(0);
      paveA->AddText(fitname.c_str());
      ((TText*)paveA->GetListOfLines()->Last())->SetTextColor(1);
      ((TText*)paveA->GetListOfLines()->Last())->SetTextFont(32);
      ((TText*)paveA->GetListOfLines()->Last())->SetTextSize(0.035);
      paveA->AddText(Form("Entires: %lu", (long unsigned int)(histo->GetEntries())));
      ((TText*)paveA->GetListOfLines()->Last())->SetTextColor(histo->GetLineColor());
      paveA->AddText(par0.c_str());
      ((TText*)paveA->GetListOfLines()->Last())->SetTextColor(histo->GetLineColor());
      paveA->AddText(par1.c_str());
      ((TText*)paveA->GetListOfLines()->Last())->SetTextColor(histo->GetLineColor());
      par0 = std::to_string(fabs(1 - (histo->GetFunction(funcName.c_str())->GetChisquare()/histo->GetFunction(funcName.c_str())->GetNDF()))*100);
      par0 = "Goodness : " + par0 + "%";
      paveA->AddText(par0.c_str());
      ((TText*)paveA->GetListOfLines()->Last())->SetTextColor(histo->GetLineColor());
      paveA->Draw();
      (*ca)->Update();
      gROOT->SetBatch(kFALSE);
      return true;
     }
  else {
        std::cout << __FUNCTION__ << " WARNING: Requested function not found in histogram " << histo->GetName() << std::endl;
        return false;
       }
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
    if (start != 0 && stop != 0)
       {
        for (int ga = start; ga < stop; ga++) wmod.push_back(vec->at(ga));
       }
    else wmod = *vec;
    unsigned int Wsize = wmod.size();
    if (Wsize == 0)
       {
        std::cout << __FUNCTION__ << " ERROR: Trying to calculate median of empty vector!" << std::endl;
        return -1;  // Undefined, realy.
       }
    else {
          sort(wmod.begin(), wmod.end());
          if (Wsize % 2 == 0) return (wmod.at((Wsize/2)-1) + wmod.at(Wsize/2))/2;
          else return wmod.at(ceil((float)Wsize / 2));
         }
}
// --------------------------------------------------------------------------------------------------------------
template <typename V> double LGADBase::CalcFWHM(V* vec, double median, int start, int stop)
{
    if (stop <= 0 || stop > (int)(vec->size()) || stop <= start) stop = vec->size();
    if (start < 0 || start >= (int)(vec->size())) start = 0;

    V wmod;
    if (start != 0 && stop != 0)
       {
        for (int ga = start; ga < stop; ga++) wmod.push_back(vec->at(ga));
       }
    else wmod = *vec;
    unsigned int Wsize = wmod.size();
    if (Wsize == 0)
       {
        std::cout << __FUNCTION__ << " ERROR: Trying to calculate FWHM of empty vector!" << std::endl;
        return -1;  // Undefined, realy.
       }
    else {
          if (median == -99) median = LGADBase::CalcMeadian(&wmod);
          if (median == -1)
             {
              std::cout << __FUNCTION__ << " ERROR: Median calcuulation failed, will not continue!" << std::endl;
              return -1;
             }
          sort(wmod.begin(), wmod.end());
          V wb;
          unsigned int a = 0;
          do {
              wb.push_back(wmod.at(a));
              a++;
             } 
          while (wmod.at(a) < median);
          float b = 0;
          for (; a < wmod.size(); a++) if (wmod.at(a) == median) b++;
          b = ceil(b/2);
          if (b > 0) for (int k = 0; k < b; k++) wb.push_back(median);
          V wb2;
          a = 1;
          do {
              wb2.push_back(wmod.at(wmod.size()-a));
              a++;
             } 
          while (a <= (wmod.size() - wb.size()));
          if (m_verbose >= 3) std::cout << __FUNCTION__  << " INFO: Original size: " << wmod.size() << " first half size: " << wb.size() << " second half size: " << wb2.size() << std::endl;
          double FWHM1 = LGADBase::CalcMeadian(&wb);
          double FWHM2 = LGADBase::CalcMeadian(&wb2);
          if (m_verbose >=3 ) std::wcout << __FUNCTION__ << " INFO: First median: " << FWHM1 << " second median: " << FWHM2 << " Full width half max: " << FWHM2 - FWHM1 << std::endl;
          if (FWHM2 != -1 && FWHM1 != -1) return (FWHM2 - FWHM1);
          else {
                std::cout << __FUNCTION__ << " ERROR: Failed determingnf edge of peak!" << std::endl;
                return -1;
               }
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
    for (int i = 0; i < rd; i++) mnt += (wmod.at(ceil(((float)wmod.size()/2)+i)) + wmod.at(floor(((float)wmod.size()/2)-i)))/2;
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
template <typename T, typename V> T LGADBase::MaxDensity(V* w, T res, int start, int stop)
{
    if (stop <= 0 || stop > (int)(w->size()) || stop <= start) stop = w->size();
    if (start < 0 || start >= (int)(w->size())) start = 0;

    V wmod;
    for (int ga = start; ga < stop; ga++) wmod.push_back(w->at(ga));
    sort(wmod.begin(), wmod.end());

   // if (res == -99)  // Non discreate case
   //    {

        // Clcualte the density of observations
        std::vector<double> density;
        for (unsigned int k = 0; k < (wmod.size() - 1); k++) 
            {
             if (wmod.at(k + 1) != wmod.at(k)) density.push_back(1.0/(fabs(wmod.at(k + 1) - wmod.at(k))));
             else density.push_back(0.0);
             // std::cout << __FUNCTION__ << " Density at " << k << " " << density.back() << " " << wmod.at(k) << " " << wmod.at(k + 1) << std::endl;
            }
        // Find points of zero density and correct them with respect to the frequency
        for (unsigned int gi = 0; gi < density.size(); gi++)
            {
             if (density.at(gi) == 0)
                {
                 T wmod_value = wmod.at(gi);
                 T wmod_valueaft = -99;
                 unsigned int frq = 2;
                 for (unsigned int gt = gi + 1; gt < density.size(); gt++)
                     {
                      if (density.at(gt) == 0) frq++;
                      else { wmod_valueaft = wmod.at(gt); break; }
                     }
                 density.at(gi) = frq/(wmod_valueaft-wmod_value);
                 for (unsigned int gt = gi + 1; gt < density.size(); gt++)
                     {
                      if (density.at(gt) == 0) density.at(gt) = frq / (wmod_valueaft - wmod_value);
                      else break;
                     }
                }
            }
        double maxdens = *std::max_element(density.begin(), density.end());
        // std::cout << __FUNCTION__ << " Max density: " << maxdens << " " << density.size() << std::endl;

        // Find if points of increased densitty exist 
        std::vector<unsigned int> dens_points[5];
        for (unsigned int m = 0; m < density.size(); m++)
            {
             for (unsigned int d = 0; d < 5; d++)
                 {
                  if (density.at(m) >= (0.9-(d*0.1))*maxdens) dens_points[d].push_back(m);
                  break; 
                 }
             }
        // Minimize diatances of the mean between the elements with respect to the median
        double med = CalcMeadian(&wmod);
        // std::cout << __FUNCTION__ << " Median of non-descrete vector: " << med << std::endl;
        int maxdensindx[5];
        for (unsigned int h = 0; h < 5; h++)
            {
             maxdensindx[h] = -1;
             // std::cout << __FUNCTION__ << " No of points at " << (0.9-h*0.1)*100 << "% of max density: "  << dens_points[h].size() << std::endl;
             if (dens_points[h].size() > 2)
                { 
                 double mimdist = fabs(med); 
                 double cal = 2*mimdist;
                 for (unsigned int t = 0; t < (dens_points[h].size() - 1); t++)
                     {
                      int indx = ceil((dens_points[h].at(t) + dens_points[h].at(t+1))/2);
                      cal = fabs(wmod.at(indx) - med);
                      // std::cout << __FUNCTION__ << " inside density points " << h << ":" << t << "/" << (dens_points[h].size()-1) << " " << indx << " " << cal << " " << mimdist << " " << dens_points[h].at(t) << " " << dens_points[h].at(t + 1) << std::endl;
                      // std::cout << cal - mimdist << std::endl;
                      if ((cal-mimdist)<(1e-10*mimdist))
                         {
                          mimdist = cal;                         
                          maxdensindx[h] = std::max_element(std::next(density.begin(), dens_points[h].at(t)), std::next(density.begin(), dens_points[h].at(t + 1))) - density.begin();
                          // std::cout << __FUNCTION__ << " Replacing minimization values " << mimdist << " " << maxdensindx[h] << std::endl;
                         }
                     }
                }
             else if (dens_points[h].size() == 2) maxdensindx[h] = std::max_element(std::next(density.begin(), dens_points[h].at(0)), std::next(density.begin(), dens_points[h].at(1))) - density.begin();
             // std::cout << __FUNCTION__ << " At " << (0.9 - h * 0.1) * 100 << "%, max density point: " << maxdensindx[h] << std::endl;
            }
        // Find mean of element closest to median and return it
        T masdensval = 0.0;
        float crt = 0.;
        for (unsigned int o = 0; o < 5; o++) if (maxdensindx[o] != -1) { masdensval += wmod.at(maxdensindx[o]); crt += 1.;}
        // std::cout << __FUNCTION__ << " we are at the end, " << masdensval << " " << crt << " " << maxdens << std::endl;
        if (crt > 0) return masdensval/crt;
        else return wmod.at(std::max_element(density.begin(), density.end()) - density.begin());

      // }
   /* else {
          // will work only if all bins of the binned magnitude are occupied, if not it won't work.....
          std::vector<unsigned int> wmod_count;
          V wmod_value;
          std::cout << __FUNCTION__ << " we are on the discrete case with resolution" << res << std::endl;
          for (unsigned int i = 0; i < wmod.size(); i++)
              {
               wmod_value.push_back(wmod.at(i));
               unsigned int j = i;
               do j++;
               while (wmod.at(j) == wmod_value.at(i) || j < wmod.size());
               wmod_count.push_back(j - i);
              }
          // Find out if we have contiues elements or if we are missing some discrete elements
          float miselem = 0.0;
          for (unsigned int a = ceil(0.2*wmod_value.size()); a < ceil(0.8*wmod_value.size()); a++)
              {
               if (fabs(wmod_value.at(a + 1) - wmod_value.at(a)) > res) miselem+=1;
              }
          miselem = miselem/(ceil(0.8*wmod_value.size())-ceil(0.2*wmod_value.size()));
          int maxfrq = std::max_element(wmod_count.begin(), wmod_count.end()) - wmod_count.begin();
          return wmod_value.at(maxfrq);
         }*/
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
template double LGADBase::CalcFWHM<vector<int> >(std::vector<int>*, double, int, int);
template double LGADBase::CalcFWHM<vector<float> >(std::vector<float>*, double, int, int);
template double LGADBase::CalcFWHM<vector<double> >(std::vector<double>*, double, int, int);
template double LGADBase::CalcMeadian<vector<int> >(std::vector<int> *, int, int);
template double LGADBase::CalcMeadian<vector<float> >(std::vector<float> *, int, int);
template double LGADBase::CalcMeadian<vector<double> >(std::vector<double> *, int, int);
template int LGADBase::MaxDensity<int, vector<int> >(std::vector<int> *, int, int, int);
template float LGADBase::MaxDensity<float, vector<float> >(std::vector<float> *, float, int, int);
template double LGADBase::MaxDensity<double, vector<double> >(std::vector<double> *, double, int, int);
template vector<int> LGADBase::OutlierReject<vector<int> >(std::vector<int> *, unsigned int, float, int, int);
template vector<float> LGADBase::OutlierReject<vector<float> >(std::vector<float> *, unsigned int, float, int, int);
template vector<double> LGADBase::OutlierReject<vector<double> >(std::vector<double> *, unsigned int, float, int, int);
template int LGADBase::CalResolution<int, vector<int> >(std::vector<int> *, unsigned int, int, int);
template float LGADBase::CalResolution<float, vector<float> >(std::vector<float> *, unsigned int, int, int);
template double LGADBase::CalResolution<double, vector<double> >(std::vector<double> *, unsigned int, int, int);
