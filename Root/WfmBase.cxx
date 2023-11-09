/*
* WfmBase.cxx
*
*
*      Author: Gkougkousis Evangelos - Leonidas
*              egkougko@cern.ch
*                CERN-GENEVA
*/

#include "../WfmUtils/WfmBase.h"

//#if !defined(__CINT__)
//ClassImp(WfmBase);
//#endif
//
//#if !defined(__CLING__)
//ClassImp(WfmBase);
//#endif
//
//#ifdef __CINT__
//#pragma link C++ class WfmBase;
//#endif
//
//#ifdef __ROOTCLING__
//#pragma link C++ class WfmBase;
//#endif

WfmBase::WfmBase()
{
    // This is the default constructor, so we use default settings
    m_instrument = Unasigned;
    m_TrnsCorr = false;
    m_TrackComb = false;
    Initialize();
    m_nchan = 0;
    m_event = 0;
    m_EvnTime = (double)std::time(0);
    m_fitopt = "root";
}
// --------------------------------------------------------------------------------------------------------------
WfmBase::~WfmBase()
{
}
// --------------------------------------------------------------------------------------------------------------
WfmBase::WfmBase(WfmBase& base)
{
    InportProperties(base);
}
// --------------------------------------------------------------------------------------------------------------
void WfmBase::InportProperties(WfmBase& base)
{
    // Basic object properties
    m_instrument = base.GetInstrument();
    m_channels = base.GetChannels();
    m_DUTChsNames = base.GetDUTNames();
    m_DUTChsBrd = base.GetDUTBoards();
    m_DUTChsTrns = base.GetDUTransImps();
    m_DUTChsAmp = base.GetDUTSecStages();
    m_DUTChsAmpGn = base.GetDUTAmpGains();
    m_DUTChsTrigg = base.GetDUTriggs();
    m_DUTChsFrc = base.GetDUTFracts();
    m_DUTChsCaps = base.GetDUTCaps();
    m_DUTChsUpLims = base.GetDUTUpLims();
    m_DUTChsLowLims = base.GetDUTLowLims();
    m_DUTChsBinSzs = base.GetDUTBinSizes();
    m_DUTChsBinShps = base.GetDUTBinShapes();
    m_nchan = base.GetChNo();
    for (unsigned int a = 0; a < m_nchan; a++)
        { 
         m_srate.push_back(base.GetSRate(a));
         m_npoints.push_back(base.GetNPoints(a));
        } 
    // Input file names and directory
    m_datanames = base.GetInFileNames();
    m_datadir = base.GetInDataDir();
    m_infiles = base.GetInRootFiles();
    m_ext = base.GetExtention();
    m_ofname = base.GetOutFileName();
    m_ofdir = base.GetOutFileDir();
    m_ofile = base.GetOutRootFile();
    // Data Tree
    m_trees = base.GetRootTrees();
    m_treename = base.GetTreeName();
    // Transimpedence file name
    m_TransFileName = base.GetTransFileName();
    m_TransFile =  base.GetTransFile();
    // Boolean Switches
    m_TrackComb = base.GetTrackComb();
    m_dofft = base.GetDoFFT();
    m_hasfft = base.GetHasFFT();
    m_hastrck = base.GetHasTrck();
    m_TrckEff = base.GetTrckEff();
    m_TrnsCorr = base.GetDoTrnsCorr();
    m_convert = base.GetConvertSucess();
    m_WaveShape = base.GetDoWaveShape();
    m_hasWaveShape = base.GetHasWaveShape();
    // Debug level
    m_verbose = base.GetVerbosity();
    // Basic analyis options
    m_stage = base.GetAnaStage();
    m_filelvl = base.GetFileLVL();
    m_fitopt = base.GetFitMethode();
    m_trkpck = base.GetTrackPackage();
    m_event = base.GetEvntNo();
    m_TrsHists = base.GetTransHistos();
    m_oscdel = base.GetScopeDelay();
    m_trigclk = base.GetTrigClk();
    // Channel level analysis cuts
    m_ChVoltCuts = base.GetChVoltCuts();
    m_ChNoiseCuts = base.GetChNoiseCuts();
    m_ChJitterCuts = base.GetChJitterCuts();
    m_ChChargeCuts = base.GetChChargeCuts();
    // Interchannel level analysis cuts
    m_PlaneDTs = base.GetPlaneDTs();
    m_PlaneDCs = base.GetPlaneDCs();
    m_testEvn = base.GetTestEvtNm();
    m_trackExclude = base.GetExcludeTrackFiles();
}
// --------------------------------------------------------------------------------------------------------------
void WfmBase::ExportProperties(WfmBase& base)
{
    // Basic object properties
    base.SetInstrument(m_instrument);
    base.SetChannels(m_channels);
    for (unsigned int a = 0; a < m_nchan; a++)
        { 
         if (m_DUTChsNames.size() > a) base.SetDUTName((m_DUTChsNames.at(a)).first, (m_DUTChsNames.at(a)).second);
         if (m_DUTChsBrd.size() > a) base.SetDUTBoard((m_DUTChsBrd.at(a)).first, (m_DUTChsBrd.at(a)).second);
         if (m_DUTChsTrns.size() > a) base.SetDUTransImp((m_DUTChsTrns.at(a)).first, (m_DUTChsTrns.at(a)).second);
         if (m_DUTChsAmp.size() > a) base.SetDUTSecStage((m_DUTChsAmp.at(a)).first, (m_DUTChsAmp.at(a)).second);
         if (m_DUTChsAmpGn.size() > a) base.SetDUTAmpGain((m_DUTChsAmpGn.at(a)).first, (m_DUTChsAmpGn.at(a)).second);
         if (m_DUTChsTrigg.size() > a) base.SetDUTrigg((m_DUTChsTrigg.at(a)).first, (m_DUTChsTrigg.at(a)).second);
         if (m_DUTChsFrc.size() > a) base.SetDUTFract((m_DUTChsFrc.at(a)).first, (m_DUTChsFrc.at(a)).second);
         if (m_DUTChsCaps.size() > a) base.SetDUTCap((m_DUTChsCaps.at(a)).first, (m_DUTChsCaps.at(a)).second);
         if (m_DUTChsUpLims.size() > a) base.SetDUTUpLim((m_DUTChsUpLims.at(a)).first, (m_DUTChsUpLims.at(a)).second);
         if (m_DUTChsLowLims.size() > a) base.SetDUTLowLim((m_DUTChsLowLims.at(a)).first, (m_DUTChsLowLims.at(a)).second);
         if (m_DUTChsBinSzs.size() > a) base.SetDUTBinSize((m_DUTChsBinSzs.at(a)).first, (m_DUTChsBinSzs.at(a)).second);
         if (m_DUTChsBinShps.size() > a) base.SetDUTBinShape((m_DUTChsBinShps.at(a)).first, (m_DUTChsBinShps.at(a)).second);
         if (m_srate.size() > a) base.SetSRate(m_srate.at(a), a);
         if (m_npoints.size() > a) base.SetNPoints(m_npoints.at(a), a);
        } 
    // Input file names and directory
    base.SetInFileNames(m_datanames);
    base.SetInDataDir(m_datadir);
    base.SetInRootFiles(m_infiles);
    base.SetExtention(m_ext);
    base.SetOutFileName(m_ofname);
    base.SetOutFileDir(m_ofdir);
    base.SetOutRootFile(m_ofile);
    // Data Tree
    base.SetRootTrees(m_trees);
    base.SetTreeName(m_treename);
    // Transimpedence file name
    base.SetTransFileName(m_TransFileName);
    base.SetTransFile(m_TransFile);
    // Boolean Switches
    base.SetTrackComb(m_TrackComb);
    base.SetDoFFT(m_dofft);
    base.SetHasFFT(m_hasfft);
    base.SetHasTrck(m_hastrck);
    base.SetTrckEff(m_TrckEff);
    base.SetDoTrnsCorr(m_TrnsCorr);
    base.SetConvertSucess(m_convert);
    base.SetWaveShape(m_WaveShape);
    base.SetHasWaveShape(m_hasWaveShape);
    // Debug level
    base.SetVerbose(m_verbose);
    // Basic analyis options
    base.SetAnaStage(m_stage);
    base.SetFileLVL(m_filelvl);
    base.SetFitMethode(m_fitopt);
    base.SetTrackPackage(m_trkpck);
    base.SetEvntNo(m_event);
    base.SetTransHistos(m_TrsHists);
    base.SetScopeDelay(m_oscdel);
    base.SetTrigClk(m_trigclk);
    base.SetExcludeTrackFiles(m_trackExclude);
    // Channel level analysis cuts
    for (unsigned int a = 0; a < m_ChVoltCuts.size(); a++) base.SetChMag((m_ChVoltCuts.at(a)).first, (m_ChVoltCuts.at(a)).second, 0);
    for (unsigned int a = 0; a < m_ChNoiseCuts.size(); a++) base.SetChMag((m_ChNoiseCuts.at(a)).first, (m_ChNoiseCuts.at(a)).second, 1);
    for (unsigned int a = 0; a < m_ChJitterCuts.size(); a++) base.SetChMag((m_ChJitterCuts.at(a)).first, (m_ChJitterCuts.at(a)).second, 2);
    for (unsigned int a = 0; a < m_ChChargeCuts.size(); a++) base.SetChMag((m_ChChargeCuts.at(a)).first, (m_ChChargeCuts.at(a)).second, 3);
    // Interchannel level analysis cuts
    for (unsigned int a = 0; a < m_PlaneDTs.size(); a++) 
        base.SetPlaneDMgt(((m_PlaneDTs.at(a)).first).first, ((m_PlaneDTs.at(a)).first).second, (m_PlaneDTs.at(a)).second, 0);
    for (unsigned int a = 0; a < m_PlaneDCs.size(); a++) 
        base.SetPlaneDMgt(((m_PlaneDCs.at(a)).first).first, ((m_PlaneDCs.at(a)).first).second, (m_PlaneDCs.at(a)).second, 1);
    base.SetTestEvn(m_testEvn);
}
// --------------------------------------------------------------------------------------------------------------
void WfmBase::SetFitMethode(std::string method)
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
void WfmBase::SetInstrument(AqInstrument instr)
{    
    if (instr == Sampic || instr == LabTXT || instr == TestBeamBin || instr == TektronixScope || instr == LeCroyWRBin || instr == Unasigned || instr == TestBeamBin1 || instr == TestBeamBin2)
       {
        m_instrument = instr;
       }
    else {
          std::cout << __FUNCTION__ << " ERROR: instrument value " << instr
                    << " not allowed. Allowed values are in ascending order: Sampic, LabTXT, TestBeamBin, TektronixScope, LeCroyWRBin, or Unasigned. Falling back to "
                    << m_instrument << std::endl;
         }
}
// --------------------------------------------------------------------------------------------------------------
AqInstrument WfmBase::GetInstrument()
{
    if (m_instrument == Unasigned) std::cout << __FUNCTION__ << " WARNING: Instrument is not asigned!!" << std::endl;
    return m_instrument;
}
// --------------------------------------------------------------------------------------------------------------
void WfmBase::SetTransFileName(TString filename) 
{
     m_TransFileName = filename;
}
//---------------------------------------------------------------------------------------------------------------
void WfmBase::SetDUTName(int ChId, std::string Name)
{
    if (1 < ChId && ChId > 128)
       {
        if (Name.length() != 0 && !(trim(Name, " ")).empty()) m_DUTChsNames.push_back(std::make_pair(ChId, Name));
        else std::cout << __FUNCTION__ << " WARNING: Cannot set empty name for channel " << ChId << "!" << std::endl;
       }
    else std::cout << __FUNCTION__ << " WARNING: Cahnnel ID not supported, 1 - 128 channels allowed." << std::endl;
}
//---------------------------------------------------------------------------------------------------------------
void WfmBase::SetDUTNames(std::vector<std::string> DUTNames)
{
    for (unsigned int k = 0; k < DUTNames.size(); k++)
        {
         if ((DUTNames.at(k)).length() != 0 && !(trim((DUTNames.at(k)), " ")).empty()) m_DUTChsNames.push_back(std::make_pair(-99, (DUTNames.at(k))));
         else {
               std::cout << __FUNCTION__ << " WARNING: Cannot set empty name on channel name vector entry " << k << ", ignoring..." << std::endl;
               m_DUTChsNames.push_back(std::make_pair(-99, "invalid"));
              }
        }
}
//---------------------------------------------------------------------------------------------------------------
void WfmBase::SetDUTBoard(int ChId, AqBoard Brd)
{
    if (1 < ChId && ChId > 128)
       {
        if (Brd == SingleCh || Brd == FourCh || Brd == IN2P3 || Brd == KU || Brd == SiPM) m_DUTChsBrd.push_back(std::make_pair(ChId, Brd));
        else std::cout << __FUNCTION__ << " WARNING: Board type for channel " << ChId << " not supported!" << std::endl;
       }
    else std::cout << __FUNCTION__ << " WARNING: Cahnnel ID not supported, 1 - 128 channels allowed." << std::endl;
}
//---------------------------------------------------------------------------------------------------------------
void WfmBase::SetDUTBoards(std::vector<AqBoard> Brds)
{
    for (unsigned int k = 0; k < Brds.size(); k++)
        {
         if (Brds.at(k) == SingleCh || Brds.at(k) == FourCh || Brds.at(k) == IN2P3 || Brds.at(k) == KU || Brds.at(k) == SiPM) m_DUTChsBrd.push_back(std::make_pair(-99, Brds.at(k)));
         else {
               std::cout << __FUNCTION__ << " WARNING: Board type on channel board type vector entry " << k << " not supported, ignoring..." << std::endl;
               m_DUTChsBrd.push_back(std::make_pair(-99, undef));
              }
        }
}
//---------------------------------------------------------------------------------------------------------------
void WfmBase::SetDUTransImp(int ChId, float Transimp)
{
    if (1 < ChId && ChId > 128)
       {
        if (Transimp > 0) m_DUTChsTrns.push_back(std::make_pair(ChId, Transimp));
        else std::cout << __FUNCTION__ << " WARNING: Invalid transimpedence value for channel " << ChId << "!" << std::endl;
       }
    else std::cout << __FUNCTION__ << " WARNING: Cahnnel ID not supported, 1 - 128 channels allowed." << std::endl;
}
//---------------------------------------------------------------------------------------------------------------
void WfmBase::SetDUTransImps(std::vector<float> Transimps)
{
    for (unsigned k = 0; k < Transimps.size(); k++)
        {
         if (Transimps.at(k) > 0) m_DUTChsTrns.push_back(std::make_pair(-99, Transimps.at(k)));
         else {
               std::cout << __FUNCTION__ << " WARNING: Invalid transimpedence value on channel transimpedemce vector entry " << k << ", ignoring..." << std::endl;
               m_DUTChsTrns.push_back(std::make_pair(-99, 99));
              }
        }
}
//---------------------------------------------------------------------------------------------------------------
void WfmBase::SetDUTSecStage(int ChId, SecStage Amp)
{
    if (1 < ChId && ChId > 128)
       {
        if (Amp == MinCircuits || Amp == Particulars || Amp == CIVIDEC || Amp == none) m_DUTChsAmp.push_back(std::make_pair(ChId, Amp));
        else std::cout << __FUNCTION__ << " WARNING: Amplifier type for channel " << ChId << " not supported!" << std::endl;
       }
    else std::cout << __FUNCTION__ << " WARNING: Cahnnel ID not supported, 1 - 128 channels allowed." << std::endl;
}
//---------------------------------------------------------------------------------------------------------------
void WfmBase::SetDUTSecStages(std::vector<SecStage> Amps)
{
    for (unsigned int k = 0; k < Amps.size(); k++)
        {
         if (Amps.at(k) == MinCircuits || Amps.at(k) == Particulars || Amps.at(k) == CIVIDEC || Amps.at(k) == none) m_DUTChsAmp.push_back(std::make_pair(-99, Amps.at(k)));
         else {
               std::cout << __FUNCTION__ << " WARNING: Invalid amplifier type on channel amplifier type vector entry " << k << ", ignoring..." << std::endl;
               m_DUTChsAmpGn.push_back(std::make_pair(-99, und));
              }
        }
}
//---------------------------------------------------------------------------------------------------------------
void WfmBase::SetDUTAmpGain(int ChId, int gain)
{
    if (1 < ChId && ChId > 128)
       {
        if (gain > 1) m_DUTChsAmpGn.push_back(std::make_pair(ChId, gain));
        else std::cout << __FUNCTION__ << " WARNING: Invalid amplifier gain value for channel " << ChId << "!" << std::endl;
       }
    else std::cout << __FUNCTION__ << " WARNING: Cahnnel ID not supported, 1 - 128 channels allowed." << std::endl;
}
//---------------------------------------------------------------------------------------------------------------
void WfmBase::SetDUTAmpGains(std::vector<int> Gains)
{
    for (unsigned int k = 0; k < Gains.size(); k++)
        {
         if (Gains.at(k) > 1) m_DUTChsAmpGn.push_back(std::make_pair(-99, Gains.at(k)));
         else {
               std::cout << __FUNCTION__ << " WARNING: Invalid amplifier gain value on channel gain vector entry " << k << ", ignoring..." << std::endl;
               m_DUTChsAmpGn.push_back(std::make_pair(-99, 99));
              }
        }
}
//---------------------------------------------------------------------------------------------------------------
void WfmBase::SetDUTrigg(int ChId, float trigg)
{
    if (1 < ChId && ChId > 128)
       {
        if (trigg != 0) m_DUTChsTrigg.push_back(std::make_pair(ChId, trigg));
        else std::cout << __FUNCTION__ << " WARNING: Invalid trigger voltage value for channel " << ChId << "!" << std::endl;
       }
    else std::cout << __FUNCTION__ << " WARNING: Cahnnel ID not supported, 1 - 128 channels allowed." << std::endl;
}
//---------------------------------------------------------------------------------------------------------------
void WfmBase::SetDUTriggs(std::vector<float> Triggs)
{
    for (unsigned int k = 0; k < Triggs.size(); k++)
        {
         if (Triggs.at(k) != 0) m_DUTChsTrigg.push_back(std::make_pair(-99, Triggs.at(k)));
         else {
               std::cout << __FUNCTION__ << " WARNING: Invalid trigger voltage value on channel trigger vector entry " << k << ", ignoring..." << std::endl;
               m_DUTChsTrigg.push_back(std::make_pair(-99, 99));
              }
        }
}
//---------------------------------------------------------------------------------------------------------------
void WfmBase::SetDUTFract(int ChId, float frac)
{
    if (1 < ChId && ChId > 128)
       {
        if (frac >= 0.5 && frac <=0.95) m_DUTChsFrc.push_back(std::make_pair(ChId, frac));
        else std::cout << __FUNCTION__ << " WARNING: Invalid CFD fraction value for channel " << ChId << "!" << std::endl;
       }
    else std::cout << __FUNCTION__ << " WARNING: Cahnnel ID not supported, 1 - 128 channels allowed." << std::endl;
}
//---------------------------------------------------------------------------------------------------------------
void WfmBase::SetDUTFracts(std::vector<float> Fracs)
{
    for (unsigned int k = 0; k < Fracs.size(); k++)
        {
         if (Fracs.at(k) >= 0.5 && Fracs.at(k) <=0.95) m_DUTChsFrc.push_back(std::make_pair(-99, Fracs.at(k)));
         else {
               std::cout << __FUNCTION__ << " WARNING: Invalid CFD fraction value on channel CFD vector entry " << k << ", ignoring..." << std::endl;
               m_DUTChsFrc.push_back(std::make_pair(-99, 99));
              }
        }
}
//---------------------------------------------------------------------------------------------------------------
void WfmBase::SetDUTCap(int ChId, double cap)
{
    if (1 < ChId && ChId > 128)
       {
        if (cap > 0) m_DUTChsCaps.push_back(std::make_pair(ChId, cap));
        else std::cout << __FUNCTION__ << " WARNING: Invalid capacitance value for channel " << ChId << "!" << std::endl;
       }
    else std::cout << __FUNCTION__ << " WARNING: Cahnnel ID not supported, 1 - 128 channels allowed." << std::endl;
}
//---------------------------------------------------------------------------------------------------------------
void WfmBase::SetDUTCaps(std::vector<double> Caps)
{
    for (unsigned int k = 0; k < Caps.size(); k++)
        {
         if (Caps.at(k) > 0) m_DUTChsCaps.push_back(std::make_pair(-99, Caps.at(k)));
         else {
               std::cout << __FUNCTION__ << " WARNING: Invalid capacitance value on channel capacitance vector entry " << k << ", ignoring..." << std::endl;
               m_DUTChsCaps.push_back(std::make_pair(-99, 99));
              }
        }
}
//---------------------------------------------------------------------------------------------------------------
void WfmBase::SetDUTUpLim(int ChId, std::pair<double, double> UpLim)
{
    if (1 < ChId && ChId > 128)
       {
        if (UpLim.first >= -10 && UpLim.first <= 10 && UpLim.second >= -10 && UpLim.second <= 10)
           { 
            bool fnd = false;
            for (unsigned int k = 0; k < m_DUTChsLowLims.size(); k++)
                {
                 if ((m_DUTChsLowLims.at(k)).first == ChId)
                    {
                     fnd = true;
                     if (((m_DUTChsLowLims.at(k)).second).first > UpLim.first && UpLim.second > ((m_DUTChsLowLims.at(k)).second).second)
                        {
                         m_DUTChsUpLims.push_back(std::make_pair(ChId, UpLim));
                        }
                     else {
                           std::cout << __FUNCTION__ << " WARNING: Channel " << ChId;
                           if (((m_DUTChsLowLims.at(k)).second).first <= UpLim.first) 
                              {
                               std::cout << " x high limit (" << UpLim.first << ") set higher than the low limit (" 
                                         << ((m_DUTChsLowLims.at(k)).second).first << "), setting to -10!"<< std::endl;
                               m_DUTChsUpLims.push_back(std::make_pair(ChId, std::make_pair(-10, UpLim.second)));
                              }
                           else {
                                 std::cout << " y high limit (" << UpLim.second << ") set lower than the low limit ("
                                           << ((m_DUTChsLowLims.at(k)).second).second << "), setting to 10!" << std::endl;
                                 m_DUTChsUpLims.push_back(std::make_pair(ChId, std::make_pair(UpLim.first, 10)));
                                }
                          }
                     break;
                    }
                }
            if (!fnd) m_DUTChsUpLims.push_back(std::make_pair(ChId, UpLim));
           }
        else {
              if (UpLim.first < -10 || UpLim.first > 10) std::cout << __FUNCTION__ << " WARNING: Channel " << ChId 
                                                                   << " x upper limit out of bounds (-10, 10)!" << std::endl;
              else  std::cout << __FUNCTION__ << " WARNING: Channel " << ChId << " y upper limit out of bounds (-10, 10)!" << std::endl;
             }
       }
    else std::cout << __FUNCTION__ << " WARNING: Cahnnel ID not supported, 1 - 128 channels allowed." << std::endl;
}
//---------------------------------------------------------------------------------------------------------------
void WfmBase::SetDUTUpLims(std::vector<std::pair<double, double>> upLims)
{
    for (unsigned int k = 0; k < upLims.size(); k++)
        {
         if ((upLims.at(k)).first >= -10 && (upLims.at(k)).first <= 10 && (upLims.at(k)).second >= -10 && (upLims.at(k)).second <= 10)
            {
             if (m_DUTChsLowLims.size() > k)
                {
                 if (((m_DUTChsLowLims.at(k)).second).first > (upLims.at(k)).first && (upLims.at(k)).second > ((m_DUTChsLowLims.at(k)).second).second)
                    {
                     m_DUTChsUpLims.push_back(std::make_pair(-99, upLims.at(k)));
                    }
                 else {
                       std::cout << __FUNCTION__ << " WARNING: DUT vector entry no " << k;
                       if (((m_DUTChsLowLims.at(k)).second).first <= (upLims.at(k)).first)
                          {
                           std::cout << " x high limit (" << (upLims.at(k)).first << ") set higher than the low limit ("
                                     << ((m_DUTChsLowLims.at(k)).second).first << "), setting to -10!"<< std::endl;
                           m_DUTChsUpLims.push_back(std::make_pair(-99, std::make_pair(-10, (upLims.at(k)).second)));
                          }
                       else {
                             std::cout << " y high limit (" << (upLims.at(k)).second << ") set lower than the low limit ("
                                       << ((m_DUTChsLowLims.at(k)).second).second << "), setting to 10!" << std::endl;
                             m_DUTChsUpLims.push_back(std::make_pair(- 99, std::make_pair((upLims.at(k)).first, 10)));
                            }
                      }
                }
             else m_DUTChsUpLims.push_back(std::make_pair(-99, upLims.at(k)));
            }
         else {
               std::cout << __FUNCTION__ << " WARNING: DUT vector entry no " << k;
               if ((upLims.at(k)).first < -10 || (upLims.at(k)).first > 10)
                  {
                   std::cout << " x upper limit out of bounds (-10, 10), ignoring...." << std::endl;
                   m_DUTChsUpLims.push_back(std::make_pair(-99, std::make_pair(-10, 10)));
                  }
               else {
                     std::cout << " y upper limit out of bounds (-10, 10), ignoring..." << std::endl;   
                     m_DUTChsUpLims.push_back(std::make_pair(-99, std::make_pair(-10, 10)));
                    }
              }
        }
}
//---------------------------------------------------------------------------------------------------------------
void WfmBase::SetDUTLowLim(int ChId, std::pair<double, double> LowLim)
{
    if (1 < ChId && ChId > 128)
       {
        if (LowLim.first >= -10 && LowLim.first <= 10 && LowLim.second >= -10 && LowLim.second <= 10)
           { 
            bool fnd = false;
            for (unsigned int k = 0; k < m_DUTChsUpLims.size(); k++)
                {
                 if ((m_DUTChsUpLims.at(k)).first == ChId)
                    {
                     fnd = true;
                     if (((m_DUTChsUpLims.at(k)).second).first < LowLim.first && LowLim.second < ((m_DUTChsUpLims.at(k)).second).second)
                        {
                         m_DUTChsLowLims.push_back(std::make_pair(ChId, LowLim));
                        }
                     else {
                           std::cout << __FUNCTION__ << " WARNING: Channel " << ChId;
                           if (((m_DUTChsUpLims.at(k)).second).first >= LowLim.first)
                              {
                               std::cout << " x low limit (" << LowLim.first << ") set lower than the high limit ("
                                         << ((m_DUTChsUpLims.at(k)).second).first << "), setting to 10!"<< std::endl;
                               m_DUTChsLowLims.push_back(std::make_pair(ChId, std::make_pair(10, LowLim.second)));
                              }
                           else {
                                 std::cout << " y low limit (" << LowLim.second << ") set higher than the high limit ("
                                           << ((m_DUTChsUpLims.at(k)).second).second << "), setting to -10!" << std::endl;
                                 m_DUTChsLowLims.push_back(std::make_pair(ChId, std::make_pair(LowLim.first, -10)));
                                }
                          }
                     break;
                    }
                }
            if (!fnd) m_DUTChsLowLims.push_back(std::make_pair(ChId, LowLim));
           }
        else {
              if (LowLim.first < -10 || LowLim.first > 10) std::cout << __FUNCTION__ << " WARNING: Channel " << ChId
                                                                     << " x upper limit out of bounds (-10, 10)!" << std::endl;
              else  std::cout << __FUNCTION__ << " WARNING: Channel " << ChId << " y upper limit out of bounds (-10, 10)!" << std::endl;
             }
       }
    else std::cout << __FUNCTION__ << " WARNING: Cahnnel ID not supported, 1 - 128 channels allowed." << std::endl;
}
//---------------------------------------------------------------------------------------------------------------
void WfmBase::SetDUTLowLims(std::vector<std::pair<double, double>> lowLims)
{
    for (unsigned int k = 0; k < lowLims.size(); k++)
        {
         if ((lowLims.at(k)).first >= -10 && (lowLims.at(k)).first <= 10 && (lowLims.at(k)).second >= -10 && (lowLims.at(k)).second <= 10)
            {
             if (m_DUTChsUpLims.size() > k)
                {
                 if (((m_DUTChsUpLims.at(k)).second).first < (lowLims.at(k)).first && (lowLims.at(k)).second <((m_DUTChsUpLims.at(k)).second).second)
                    {
                     m_DUTChsLowLims.push_back(std::make_pair(-99, lowLims.at(k)));
                    }
                 else {
                       std::cout << __FUNCTION__ << " WARNING: DUT vector entry no " << k;
                       if (((m_DUTChsUpLims.at(k)).second).first >= (lowLims.at(k)).first)
                          {
                           std::cout << " x low limit (" << (lowLims.at(k)).first << ") set lower than the high limit ("
                                     << ((m_DUTChsUpLims.at(k)).second).first << "), setting to 10!"<< std::endl;
                           m_DUTChsLowLims.push_back(std::make_pair(-99, std::make_pair(10, (lowLims.at(k)).second)));
                          }
                       else {
                             std::cout << " y low limit (" << (lowLims.at(k)).second << ") set higher than the high limit ("
                                       << ((m_DUTChsUpLims.at(k)).second).second << "), setting to -10!" << std::endl;
                             m_DUTChsLowLims.push_back(std::make_pair(- 99, std::make_pair((lowLims.at(k)).first, -10)));
                            }
                      }
                }
             else m_DUTChsLowLims.push_back(std::make_pair(-99, lowLims.at(k)));
            }
         else {
               std::cout << __FUNCTION__ << " WARNING: DUT vector entry no " << k;
               if ((lowLims.at(k)).first < -10 || (lowLims.at(k)).first > 10)
                  {
                   std::cout << " x upper limit out of bounds (-10, 10), ignoring...." << std::endl;
                   m_DUTChsLowLims.push_back(std::make_pair(-99, std::make_pair(-10, 10))); 
                  }
               else {
                     std::cout << " y upper limit out of bounds (-10, 10), ignoring..." << std::endl;   
                     m_DUTChsLowLims.push_back(std::make_pair(-99, std::make_pair(-10, 10)));
                    }
              }
        }
}
//---------------------------------------------------------------------------------------------------------------
void WfmBase::SetDUTBinSize(int ChId, double size)
{
    if (1 < ChId && ChId > 128)
       {
        if (size > 0 && size < 10000) m_DUTChsBinSzs.push_back(std::make_pair(ChId, size));
        else std::cout << __FUNCTION__ << " WARNING: Invalid track bin size value for channel " << ChId << "!" << std::endl;
       }
    else std::cout << __FUNCTION__ << " WARNING: Cahnnel ID not supported, 1 - 128 channels allowed." << std::endl;
}
//---------------------------------------------------------------------------------------------------------------
void WfmBase::SetDUTBinSizes(std::vector<double> sizes)
{
    for (unsigned int k = 0; k < sizes.size(); k++)
        {
         if (sizes.at(k) > 0 && sizes.at(k) < 10000) m_DUTChsBinSzs.push_back(std::make_pair(-99, sizes.at(k)));
         else {
               std::cout << __FUNCTION__ << " WARNING: Invalid track bin size value on vector entry " << k << ", ignoring..." << std::endl;
               m_DUTChsBinSzs.push_back(std::make_pair(-99, 10000));
              }
        }
}
//---------------------------------------------------------------------------------------------------------------
void WfmBase::SetDUTBinShape(int ChId, std::string shape)
{
    if (1 < ChId && ChId > 128)
       {
        if (shape == "square" || shape == "circle") m_DUTChsBinShps.push_back(std::make_pair(ChId, shape));
        else std::cout << __FUNCTION__ << " WARNING: Invalid bin shape for channel " << ChId << ", please specify either square or circle!" << std::endl;
       }
    else std::cout << __FUNCTION__ << " WARNING: Cahnnel ID not supported, 1 - 128 channels allowed." << std::endl;
}
//---------------------------------------------------------------------------------------------------------------
void WfmBase::SetDUTBinShapes(std::vector<std::string> shapes)
{
    for (unsigned int k = 0; k < shapes.size(); k++)
        {
         if (shapes.at(k) == "square" || shapes.at(k) == "circle") m_DUTChsBinShps.push_back(std::make_pair(-99, shapes.at(k)));
         else {
               std::cout << __FUNCTION__ << " WARNING: Invalid bin shape on vector entry " << k << ", reverting to square bins..." << std::endl;
               m_DUTChsBinShps.push_back(std::make_pair(-99, "square"));
              }
        }
}
//---------------------------------------------------------------------------------------------------------------
bool WfmBase::GetTestEvn(std::pair<int, int> Evn)
{
    if (m_testEvn.second == -1)
       {
        if (m_testEvn.first == -1) return false;
        else {
              if (int(Evn.first) != m_testEvn.first) return true;
              else return false;
             }
       }
    else {
          if (m_testEvn.second != Evn.second) return true;
          else {
                if (int(Evn.first) != m_testEvn.first) return true;
                else return false;
               }
         }
}
//---------------------------------------------------------------------------------------------------------------
void WfmBase::SetChMag(int ChId, double ChMag, int Qt)
{
    if (0 < ChId && ChId < 128)
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
void WfmBase::SetChMags(std::vector<double> ChMags, int Qt)
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
void WfmBase::SetPlaneDMgt(int PlaneId1, int PlaneId2, double DMgt, int Qt)
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
void WfmBase::SetPlaneDMgts(std::vector<double> DMgts, int Qt)
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
bool WfmBase::OpenTransFile(const TString &filename)
{
    if (filename == "" || ((trim((std::string)filename, " ")).empty()))
       { 
        std::cout << __FUNCTION__ << " WARNING: No transimpedence file set, ignoring transimpedence correction..." << std::endl;
        m_TrnsCorr = false; 
        return false; 
       }
    else {
          TFile *m_TransFile = TFile::Open(filename, "READ");
          if (m_TransFile && !(m_TransFile->IsZombie()))
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
bool WfmBase::SetSRate(Long64_t rate, unsigned int ch)
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
Long64_t WfmBase::GetSRate(unsigned int ch)
{
    return m_srate.at(ch);
}
// --------------------------------------------------------------------------------------------------------------
bool WfmBase::SetNPoints(unsigned int points, unsigned int ch)
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
unsigned int WfmBase::GetNPoints(unsigned int ch)
{
    return m_npoints.at(ch);
}
// --------------------------------------------------------------------------------------------------------------
void WfmBase::Initialize()
{
    m_trees.clear();
    m_verbose = 0;
    m_trigclk = -99;
    m_oscdel.clear();
    m_stage = 1;
    m_filelvl = std::make_pair(0, 0);
    m_channels.clear();
    m_channels.push_back(0);
    m_scale.clear();
    m_scale.push_back(0);
    m_physt.clear();
    m_physt.push_back(0);
    m_ordrt.clear();
    m_ordrt.push_back(0);
    m_npoints.clear();
    m_scope.clear();
    m_srate.clear();
    m_DUTChsBrd.clear();
    m_DUTChsAmp.clear();
    m_DUTChsFrc.clear();
    m_DUTChsTrigg.clear();
    m_DUTChsTrns.clear();
    m_DUTChsAmpGn.clear();
    m_DUTChsNames.clear();
    m_DUTChsCaps.clear();
    m_DUTChsUpLims.clear();
    m_DUTChsLowLims.clear();
    m_DUTChsBinSzs.clear();
    m_DUTChsBinShps.clear();
    m_PlaneDTs.clear();
    m_PlaneDCs.clear();
    m_ChNoiseCuts.clear();
    m_ChVoltCuts.clear();
    m_ChChargeCuts.clear();
    m_ChJitterCuts.clear();
    m_trackExclude.clear();
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
    if (m_TrackComb)
       {
        m_Trackdatadir.Clear();
        m_Trackdataname.clear();
        m_trackTrees.clear();
        m_ntp_planes.clear();
        m_ntp_Ref.clear();
        m_ntp_DUT.clear();
       }
    m_evnt1 = 0;
    m_evnt2 = 0;
    m_testEvn = std::make_pair(-1, -1);
    m_TrckEff = false;
    m_convert = false;
    m_WaveShape = false;
    m_hasWaveShape = std::make_pair(99, 99);
    m_dofft = false;
    m_hasfft.resize(5, std::make_pair(true, false));
    m_hastrck = false;
    m_ext = "";
    m_datadir = "";
    m_datanames.clear();
    m_treename = "wfm";
    m_ofdir = "";
    m_ofname = "";
    if (m_TrnsCorr) OpenTransFile(m_TransFileName);
    if (m_verbose >= 2) std::cout << __FUNCTION__ << " INFO: Initializing sampling rate to default value : " << m_srate.at(0) / 1e6 << " MSamples/sec and no. of points to : " << m_npoints.at(0) << std::endl;
}
// --------------------------------------------------------------------------------------------------------------
void WfmBase::SetStartStopEvnt(int Evnt1, int Evnt2)
{
    if (Evnt1 >= 0) m_evnt1 = Evnt1;
    else std::cout << __FUNCTION__ << " WARNING: Start event value not coherent, setting to zero!" << std::endl;
    if (Evnt2 > Evnt1 && Evnt2 > 0) m_evnt2 = Evnt2;
    else if (Evnt2 != 0) { m_evnt2 = 0; std::cout << __FUNCTION__ << " WARNING: Stop event value not coherent, setting to zero!" << std::endl;}
    else m_evnt2 = 0;
}
// --------------------------------------------------------------------------------------------------------------
std::pair <unsigned int, unsigned int> WfmBase::GetStartStopEvnt()
{
    return std::make_pair(m_evnt1, m_evnt2);
}
// --------------------------------------------------------------------------------------------------------------
void WfmBase::SetVectorSize(unsigned int nch, unsigned int stage)
{
    m_physt.clear();
    m_ordrt.clear();
    m_t.clear();
    m_w.clear();
    m_scale.clear();
    m_EvnTimeCh.clear();
    m_npoints.clear();
    m_srate.clear();
    m_scope.clear();
    m_t.resize(nch, nullptr);
    m_w.resize(nch, nullptr);
    // This opperation cannot be done be resize, the adress atrributed is idential for all elements.
    if (stage == 1)
       {
        for (unsigned int qrt = 0; qrt < nch; qrt++)
            {
             m_w.at(qrt) = new std::vector<double>();
             m_t.at(qrt) = new std::vector<double>();
            }
       }

    if (m_instrument == TektronixScope || m_instrument == Unasigned) m_EvnTimeCh.resize(nch);
    if (m_instrument == Sampic)
       {
        m_physt.resize(nch, -99.);
        m_ordrt.resize(nch, -99.);
        m_npoints.resize(1, 99);
        m_srate.resize(1, -99.);
       }
    else {
          m_npoints.resize(nch, 99);
          m_srate.resize(nch, -99.);
          m_scale.resize(nch, -99.);
          m_scope.resize(nch, 99);
         }
}
// --------------------------------------------------------------------------------------------------------------
int WfmBase::Addoriel(int val)
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
bool WfmBase::ProgressBar(Long64_t evt, Long64_t total)
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
bool WfmBase::SetRootTrees(std::vector<TFile*> files, std::vector<unsigned int> &indx, std::string name)
{
    int n;
    TTree* data_tree = NULL;
    TKey* key = NULL;
    
    for (unsigned int r = 0; r < files.size(); r++)
        {
         TIter nextkey(files.at(r)->GetListOfKeys());
         data_tree = NULL;
         key = NULL;
         n = 0;
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
         if (n == 1) m_trees.push_back(data_tree); 
         else {
               std::cout << __FUNCTION__ << " WARNIGIN: Data file " << files.at(r)->GetName() << " containts " << n << " trees with name " << name << ". Skipping..." << std::endl;
               indx.erase(indx.begin() + r);
               continue;
              }
        }

    if (indx.size() > 0 && m_trees.size() > 0) return true;
    else return false;
}
// --------------------------------------------------------------------------------------------------------------
int WfmBase::SetInRootFiles(std::vector<TString> &files, std::vector<unsigned int> &indx)
{
    int opfiles = 0;
    for (unsigned int k = 0; k < files.size(); k++)
        {
         m_infiles.push_back(TFile::Open(GetInDataDir() + files.at(k) + "." + GetExtention()));
         if (m_infiles.back() && !(m_infiles.back()->IsZombie())) { opfiles++; indx.push_back(k); }
         else {
               std::cout << __FUNCTION__ << " ERROR: Data file " << GetInDataDir() + (GetInFileNames()).at(k) << "." << GetExtention() << " cannot be openned!" << std::endl;
               opfiles--;
               m_infiles.pop_back();
               files.erase(files.begin()+k);
               k--;
              }
        }
    return opfiles;
}
// --------------------------------------------------------------------------------------------------------------
int WfmBase::SetInRootFile(TFile* file) 
{ 
    m_infiles.push_back(file); 
    if (m_infiles.back()) return 1;
    else {
          std::cout << __FUNCTION__ << " ERROR: Data file " << file->GetName() << " cannot be openned! Probably zombie???" << std::endl;
          return -1;
         }
}
// --------------------------------------------------------------------------------------------------------------
// One of the methodes to set the output root file. Another exists that only takes as iunput a TFile
bool WfmBase::SetOutRootFile(const char* file)
{
    m_ofile = new TFile(file, "recreate");
    if (m_ofile) return true;
    else {
          std::cout << __FUNCTION__ << " ERROR: Could not create output file (" << m_ofile << ")!" << std::endl;
          return false;
         }
}
// --------------------------------------------------------------------------------------------------------------
bool WfmBase::SetOutRootFile(TFile* file)
{
    m_ofile = file;
    if (m_ofile) return true;
    else {
          std::cout << __FUNCTION__ << " ERROR: Could not create output file (" << m_ofile << ")!" << std::endl;
          return false;
         }
}
// --------------------------------------------------------------------------------------------------------------
// Helper Function integrated to trim strings by reoving spaces or other defined characters from it's begining and end. 
std::string WfmBase::trim(const std::string str, const std::string whitespace)
{
    int strBegin = str.find_first_not_of(whitespace);
    if ((int)std::string::npos == strBegin) return ""; // no content
    unsigned int strEnd = str.find_last_not_of(whitespace);
    unsigned int strRange = strEnd - strBegin + 1;
    return str.substr(strBegin, strRange);
}
// --------------------------------------------------------------------------------------------------------------
// Function to remove characters and/or spaces from the totality, beginning and end of the string.
std::string WfmBase::reduce(const std::string str, const std::string fill, const std::string whitespace)
{
    std::string result = trim(str, whitespace); // trim first
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
int WfmBase::RecursMkDir(const char* dirname)
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
             if (CreateDir(_path) != 1) return -1;
             *p = '/';
            }
        }

    if (CreateDir(_path) != 1) return -1;

    return 1;
}
// --------------------------------------------------------------------------------------------------------------
int WfmBase::CreateDir(const char* path)
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
int WfmBase::DirExists(const char* path)
{
    struct stat info;
    if (stat(path, &info) != 0) return 0;
    else if (info.st_mode & S_IFDIR) return 1;
    else return -1;
}
// --------------------------------------------------------------------------------------------------------------
// Function insensitive to the extention having or not a dot included at astart
unsigned int WfmBase::CountFiles(const char* dir, const char* ext)
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
               // Ignore directories and windows temporary backup files
               if (!sfile->IsDirectory() && !fname.BeginsWith("~$"))
                  {
                   if (ext != nullptr && fname.EndsWith(ext)) nfiles++;
                   else if (m_instrument == Sampic && (TString(fname(0, fname.Length() - 5))).EndsWith(ext)) nfiles++;
                   else if (ext == nullptr) nfiles++;
                  }
              }
       }
    return nfiles;
}
// --------------------------------------------------------------------------------------------------------------
// Function to return vector of filenames in a directory, insensitive to extention having or not a dot as first character                                                                                                                    
std::vector<std::string> WfmBase::ListFileNames(const char* path, const char* ext)
{
    std::vector<std::string> filenames;
    std::string search_path = path;
    std::string search_ext;
    if (ext != nullptr) 
       {
        search_ext = ext;
        if (search_ext.rfind(".", 0) == 0) search_ext.erase(0, 1);
        if (m_instrument == Sampic) search_ext += "*";
       }
#ifdef _WIN32
    search_path += "*.";
    if (ext != nullptr) search_path += search_ext;
    else search_path.pop_back();
    WIN32_FIND_DATA fd;
    HANDLE hFind = ::FindFirstFile(search_path.c_str(), &fd);
    if (hFind != INVALID_HANDLE_VALUE)
       {
        do {
            // It will ingore directories and temporary hidden windws backup files
            if (!(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) && ((std::string)(fd.cFileName)).rfind("~$", 0) != 0) filenames.push_back(fd.cFileName);
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
           if (ext != nullptr && search_path.find(search_ext) != std::string::npos) filenames.push_back(pdir->d_name);
           else if (strcmp(pdir->d_name, ".") != 0 && strcmp(pdir->d_name, "..") != 0) filenames.push_back(pdir->d_name);
          }
#endif

    return filenames;
}
// --------------------------------------------------------------------------------------------------------------
bool WfmBase::SetScale(std::vector<unsigned int> channel, unsigned int nchan, std::vector<float>* scale)
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
bool WfmBase::PrintFitInfo(TH1D* histo, TCanvas** ca, std::string funcName)
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
void WfmBase::PrintBranches(TTree* tree, std::string name)
{
    // Printout registered Branches
    TIter next(tree->GetListOfBranches());
    TBranch* b2;
    std::cout << "--> Booking " << name.c_str() << " tree branches: " << std::endl << "\t";
    // Get the names of all the branches
    std::vector<string> brcname;
    while ((b2 = (TBranch*)next())) brcname.push_back(b2->GetName());
    // Find the mazimumm branch length
    unsigned int max_size = 0;
    for (unsigned int k = 0; k < brcname.size(); k++)
        { 
         if ((brcname.at(k)).size() > max_size) max_size = (brcname.at(k)).size();
        }
    // Print out the list of object names
    int l_space = 0;
    for (unsigned int k = 0; k < brcname.size(); k++)
        {
         std::cout << std::setw(max_size+2) << std::left << brcname.at(k);
         l_space++;
         if (l_space % 4 == 0) std::cout << std::endl << "\t";
        }
    if (l_space % 4 != 0) std::cout << std::endl;
}
// --------------------------------------------------------------------------------------------------------------
void WfmBase::PrintObjects(TFile* rootfile, std::string ignore)
{
    rootfile->cd();
    TListIter l_next_object(gDirectory->GetList());
    l_next_object.Reset();
    TObject* l_obj;
    std::cout << "\r--> Booking Objects: " << std::endl << "\t";
    int l_space = 0;
    while ((l_obj=l_next_object()))
          {
           TString l_objname = l_obj->GetName();
           if (l_objname.Contains(ignore)) continue;
           std::cout << std::setw(21) << std::left << l_objname ;
           l_space++;
           if (l_space %5 == 0) cout << std::endl << "\t" ;
          }
    if (l_space % 5 != 0) std::cout << std::endl;
    std::cout << std::endl;
}
// --------------------------------------------------------------------------------------------------------------
template <typename A> bool WfmBase::IsVecEqual(const std::vector<A> &first, const std::vector<A> &second)
{
    if (first.size() != second.size()) return false;
    std::unordered_multiset<double> s1(first.begin(), first.end());
    std::unordered_multiset<double> s2(second.begin(), second.end());
    if (s1 == s2) return true;
    else return false;
}
// --------------------------------------------------------------------------------------------------------------
template <typename T> T WfmBase::Derivate(const T &w, int start)
{
    int sz = w.size();
    T deriv;
    if (start <= 0) 
       {
        std::cout << __FUNCTION__ << " WARNING: Start index invalid, starting from element 0!" << std::endl;
        start = 1;
       }
    if (sz > start)
       {
        deriv.reserve(sz - start);
        for (int k = start-1; k < (sz-1); k++) deriv.push_back(w.at(k+1) - w.at(k));
       }
    else std::cout << __FUNCTION__ << " ERROR: Vector has too few elements to derivate!" << std::endl;
    return deriv;
}
// --------------------------------------------------------------------------------------------------------------
template <typename U> double WfmBase::Mean(const U &w, int start, int stop)
{
    double mean = -99.;
    int sz = 0;
    if (start <= 0 && stop <= 0) sz = w.size();
    else if (start <= 0 && stop > 0) sz = stop;
    else if (stop <= 0 && start > 0) sz = w.size() - start;
    else sz = stop - start;
    if (sz > 0)
       { 
        if (start <= 0 && stop <= 0) mean = std::accumulate(w.begin(), w.end(), 0.0) / sz;
        else if (start <= 0 && stop > 0) mean = std::accumulate(w.begin(), w.begin() + stop, 0.0) / sz;
        else if (stop <= 0 && start > 0) mean = std::accumulate(w.begin() + start, w.end(), 0.0) / sz;
        else mean = std::accumulate(w.begin() + start, w.begin() + stop, 0.0) / sz;
       }
    else std::cout << __FUNCTION__ << " ERROR: Trying to calculate mean of empty vector or inproper range!" << std::endl;
    return mean;
}
// --------------------------------------------------------------------------------------------------------------
template <typename V> double WfmBase::Stdev(const V &w, int start, int stop, double mean)
{
    double stdev = -99.;
    double E = 0;
    int sz = 0;
    if (start <= 0 && stop <= 0) sz = w.size();
    else if (start <= 0 && stop > 0) sz = stop;
    else if (stop <= 0 && start > 0) sz = w.size() - start;
    else sz = stop - start;
    if (sz > 0)
       {
        double ave = 0.0;
        if (mean == 0.0) ave = Mean(w, start, stop);
        else ave = mean;
        if ((start <= 0 && stop <= 0) || (start <= 0 && stop > 0)) for (int i = 0; i < sz; i++) E += pow((ave - w.at(i)), 2);
        else if (stop <= 0 && start > 0) for (unsigned int i = start; i < w.size(); i++) E += pow((ave - w.at(i)), 2);
        else for (int i = start; i < stop; i++) E += pow((ave - w.at(i)), 2);
        stdev = sqrt(E / (sz-1));
       }
    else std::cout << __FUNCTION__ << " ERROR: Trying to calculate standard deviation of empty vector or improper range!" << std::endl;
    return stdev;
}
// --------------------------------------------------------------------------------------------------------------
template <typename V, typename T> double WfmBase::BayesianErr(const V &w, T value)
{
    double events = w.size();
    double pass = 0;
    for (unsigned int k = 0; k < events; k++) if ((w.at(k)) >= value) pass++;
    if (events > 0) return sqrt((((pass + 1)*(pass + 2)) / ((events + 2)*(events + 3))) - (pow(pass + 1, 2) / pow(events + 2, 2)));
    else {
          std::cout << __FUNCTION__ << " ERROR: Trying to calculate baysian uncertenty of empty distribution!" << std::endl;
          return 99.;
         }
}
// --------------------------------------------------------------------------------------------------------------
template <typename V> double WfmBase::CalcMeadian(const V &vec, int start, int stop)
{
    if (stop <= 0 || stop > static_cast<int>(vec.size()) || stop <= start) stop = vec.size();
    if (start < 0 || start >= static_cast<int>(vec.size())) start = 0;

    // Temporary container for used element subset per itereation
    V wmod(vec.begin() + start, vec.begin() + stop);
    if (wmod.empty()) return -99.;

    unsigned int Wsize = wmod.size();
    if (Wsize == 0)
       {
        std::cout << __FUNCTION__ << " ERROR: Trying to calculate median of empty vector!" << std::endl;
        return -1;  // Undefined, realy.
       }
    else if (Wsize == 1)
            {
             std::cout << __FUNCTION__ << " Warning: Trying to calculate median of a single element vector!" << std::endl;
             return wmod.at(0);
            }
    else {
          sort(wmod.begin(), wmod.end());
          if (Wsize % 2 == 0) return (wmod.at((Wsize/2)-1) + wmod.at(Wsize/2))/2;
          else return wmod.at(ceil((float)Wsize / 2));
         }
}
// --------------------------------------------------------------------------------------------------------------
template <typename V> double WfmBase::CalcFWHM(const V &vec, double median, int start, int stop)
{
    if (stop <= 0 || stop > static_cast<int>(vec.size()) || stop <= start) stop = vec.size();
    if (start < 0 || start >= static_cast<int>(vec.size())) start = 0;

    // Temporary container for used element subset per itereation
    V wmod(vec.begin() + start, vec.begin() + stop);
    if (wmod.empty()) return -99.;

    unsigned int Wsize = wmod.size();
    if (Wsize == 0)
       {
        std::cout << __FUNCTION__ << " ERROR: Trying to calculate FWHM of empty vector!" << std::endl;
        return -1;  // Undefined, realy.
       }
    else if (Wsize == 1)
       {
        std::cout << __FUNCTION__ << " Warning: FWHM of a single element vector is 0!" << std::endl;
        return 0;
       }
    else {
          if (median == -99) median = CalcMeadian(wmod);
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
          double FWHM1 = CalcMeadian(wb);
          double FWHM2 = CalcMeadian(wb2);
          if (m_verbose >=3 ) std::wcout << __FUNCTION__ << " INFO: First median: " << FWHM1 << " second median: " << FWHM2 << " Full width half max: " << FWHM2 - FWHM1 << std::endl;
          if (FWHM2 != -1 && FWHM1 != -1) return (FWHM2 - FWHM1);
          else {
                std::cout << __FUNCTION__ << " ERROR: Failed determingnf edge of peak!" << std::endl;
                return -1;
               }
         }
}
// --------------------------------------------------------------------------------------------------------------
template <typename V> V WfmBase::OutlierReject(const V &w, unsigned int order, float elem, int start, int stop)
{
    if (stop <= 0 || stop > static_cast<int>(w.size()) || stop <= start) stop = w.size();
    if (start < 0 || start >= static_cast<int>(w.size())) start = 0;

    // Temporary container for used element subset per itereation
    V wmod(w.begin() + start, w.begin() + stop);
    if (wmod.empty()) return static_cast<V>(-99);
    sort(wmod.begin(), wmod.end());

    double mnt = 0.0;
    int rd = ceil((static_cast<float>(wmod.size()*elem))/2);
    for (int i = 0; i < rd; i++) mnt += (wmod.at(ceil(((float)wmod.size()/2)+i)) + wmod.at(floor(((float)wmod.size()/2)-i)))/2;
    mnt = mnt/(static_cast<float>(rd));
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
template <typename T, typename V> T WfmBase::CalResolution(const V &w, unsigned int order, int start, int stop)
{ 
    // Consistency checks
    if (stop <= 0 || stop > static_cast<int>(w.size()) || stop <= start) stop = w.size();
    if (start < 0 || start >= static_cast<int>(w.size())) start = 0;

    // temporary vector and value
    V mag;
    mag.reserve(WfmBase::Addoriel(stop-start-1));
    double l = 0;

    // Calculate distances from next element for all vector elements
    for (int ga = start; ga < stop - 1; ga++)
        {
         for (int ma = ga + 1; ma < stop; ma++)
             {
              l = fabs((double)(w.at(ma) - w.at(ga)));
              if (l > 0) mag.push_back(l);
             }
        }

    // Remove extreme elements from vectors to account for accuracy issues
    l = WfmBase::Mean(mag);
    for (unsigned int kr = 0; kr < mag.size(); kr++)
        {
         if (fabs(std::log10(mag.at(kr)/l)) >= order)
            {
             mag.erase(mag.begin() + kr);
             kr--;
            }
        }

    if (mag.empty()) return static_cast<T>(-99);
    else return *std::min_element(mag.begin(), mag.end());
}
// --------------------------------------------------------------------------------------------------------------
template <typename T, typename V> T WfmBase::MaxDensity(const V &w, int start, int stop)
{
    if (stop <= 0 || stop > static_cast<int>(w.size()) || stop <= start) stop = w.size();
    if (start < 0 || start >= static_cast<int>(w.size())) start = 0;

    // Temporary container for used element subset per itereation
    V wmod(w.begin() + start, w.begin() + stop);
    if (wmod.empty()) return static_cast<T>(-99);

    // Estimate the sum of distances of every value to the rest of the values in the vector and get the lowest
    double maxDensity = std::numeric_limits<double>::max();
    T mostProbableValue = wmod[0];
    double density = 0.0;
    for (unsigned int k = 0; k < wmod.size(); k++)
        {
         density = 0.0;
         for (unsigned le = 0; le < wmod.size(); le++) density += fabs(wmod.at(k) - wmod.at(le));
         if (density < maxDensity)
            {
             maxDensity = density;
             mostProbableValue = wmod[k];
            }
        }
    return mostProbableValue;
}
// --------------------------------------------------------------------------------------------------------------
// Explicit template instantanization
template vector<int> WfmBase::Derivate<vector<int> >(const std::vector<int>&, int);
template vector<float> WfmBase::Derivate<vector<float> >(const std::vector<float>&, int);
template vector<double> WfmBase::Derivate<vector<double> >(const std::vector<double>&, int);
template bool WfmBase::IsVecEqual(const std::vector<int>&, const std::vector<int>& );
template bool WfmBase::IsVecEqual(const std::vector<bool>&, const std::vector<bool>& );
template bool WfmBase::IsVecEqual(const std::vector<float>&, const std::vector<float>& );
template bool WfmBase::IsVecEqual(const std::vector<double>&, const std::vector<double>& );
template double WfmBase::Mean<vector<int> >(const std::vector<int>&, int, int);
template double WfmBase::Mean<vector<bool> >(const std::vector<bool>&, int, int);
template double WfmBase::Mean<vector<float> >(const std::vector<float>&, int, int);
template double WfmBase::Mean<vector<double> >(const std::vector<double>&, int, int);
template double WfmBase::Stdev<vector<int> >(const std::vector<int>&, int, int, double);
template double WfmBase::Stdev<vector<bool> >(const std::vector<bool>&, int, int, double);
template double WfmBase::Stdev<vector<float> >(const std::vector<float>&, int, int, double);
template double WfmBase::Stdev<vector<double> >(const std::vector<double>&, int, int, double);
template double WfmBase::BayesianErr<vector<int>, int>(const std::vector<int>&, int);
template double WfmBase::BayesianErr<vector<bool>, bool>(const std::vector<bool>&, bool);
template double WfmBase::CalcFWHM<vector<int> >(const std::vector<int>&, double, int, int);
template double WfmBase::CalcFWHM<vector<float> >(const std::vector<float>&, double, int, int);
template double WfmBase::CalcFWHM<vector<double> >(const std::vector<double>&, double, int, int);
template double WfmBase::CalcMeadian<vector<int> >(const std::vector<int>&, int, int);
template double WfmBase::CalcMeadian<vector<float> >(const std::vector<float>&, int, int);
template double WfmBase::CalcMeadian<vector<double> >(const std::vector<double>&, int, int);
template int WfmBase::MaxDensity<int, vector<int> >(const std::vector<int>&, int, int);
template float WfmBase::MaxDensity<float, vector<float> >(const std::vector<float>&, int, int);
template double WfmBase::MaxDensity<double, vector<double> >(const std::vector<double>&, int, int);
template vector<int> WfmBase::OutlierReject<vector<int> >(const std::vector<int>&, unsigned int, float, int, int);
template vector<float> WfmBase::OutlierReject<vector<float> >(const std::vector<float>&, unsigned int, float, int, int);
template vector<double> WfmBase::OutlierReject<vector<double> >(const std::vector<double>&, unsigned int, float, int, int);
template int WfmBase::CalResolution<int, vector<int> >(const std::vector<int>&, unsigned int, int, int);
template float WfmBase::CalResolution<float, vector<float> >(const std::vector<float>&, unsigned int, int, int);
template double WfmBase::CalResolution<double, vector<double> >(const std::vector<double>&, unsigned int, int, int);