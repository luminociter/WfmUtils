    /*
* LGADUtils.cxx
*
*
*      Author: Gkougkousis Evangelos - Leonidas
*              egkougko@cern.ch
*                CERN-GENEVA
*/

#include "../LGADUtils/LGADUtils.h"
#include "LinkDef.h"

//#if !defined(__CINT__)
//ClassImp(LGADUtils);
//#endif
//
//#if !defined(__CLING__)
//ClassImp(LGADUtils);
//#endif
//
//#ifdef __CINT__
//#pragma link C++ class LGADUtils;
//#endif
//
//#ifdef __ROOTCLING__
//#pragma link C++ class LGADUtils;
//#endif

// Default constructor creating a general base class object
LGADUtils::LGADUtils()
{
    m_Base = new LGADBase();
}
// --------------------------------------------------------------------------------------------------------------
// Default destructur deleting the general base class object
LGADUtils::~LGADUtils()
{
    delete m_Base;
    if (m_Run) delete m_Run;
    if (m_chain) delete m_chain;
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetTrackComb(bool comb)
{
    m_Base->LGADBase::SetTrackComb(comb);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetFEi4Eff(bool FEi4Eff)
{
    m_Base->LGADBase::SetFEi4Eff(FEi4Eff);
}
// --------------------------------------------------------------------------------------------------------------
bool LGADUtils::SetSRate(Long64_t rate, unsigned int ch)
{
    return m_Base->LGADBase::SetSRate(rate, ch);
}
// --------------------------------------------------------------------------------------------------------------
bool LGADUtils::SetNPoints(unsigned int points, unsigned int ch)
{
    return m_Base->LGADBase::SetNPoints(points, ch);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetFitMethode(std::string method)
{
    m_Base->LGADBase::SetFitMethode(method);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetInstrument(AqInstrument instr)
{
    if (instr == TestBeamBin1 || instr == TestBeamBin2) instr = TestBeamBin;
    m_Base->LGADBase::SetInstrument(instr);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetStartStopEvnt(int Evnt1, int Evnt2)
{
    m_Base->LGADBase::SetStartStopEvnt(Evnt1, Evnt2);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::Initialize()
{
    m_Base->LGADBase::Initialize();
}
// --------------------------------------------------------------------------------------------------------------
bool LGADUtils::ConvertData()
{
    return m_Base->LGADBase::ConvertData();
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetInDataNames(TString DataDir, TString DataName, TString ext)
{
    m_Base->LGADBase::SetInDataNames(DataDir, DataName, ext);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetOutDataNames(TString DataDir, TString DataName)
{
    m_Base->LGADBase::SetOutDataNames(DataDir, DataName);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetTrackInDataNames(TString DataDir, TString DataName)
{
    m_Base->LGADBase::SetTrackInDataNames(DataDir, DataName);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetVerbose(int verbose)
{
    m_Base->LGADBase::SetVerbose(verbose);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetWaveShape(bool shape)
{
    m_Base->LGADBase::SetWaveShape(shape);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetTransFileName(TString filename)
{
    m_Base->LGADBase::SetTransFileName(filename);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetDoTrnsCorr(bool TrnsCorr)
{
    m_Base->LGADBase::SetDoTrnsCorr(TrnsCorr);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetTreeName(std::string treename)
{
    m_Base->LGADBase::SetTreeName(treename);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetDUTName(int ChId, std::string Name)
{
    m_Base->LGADBase::SetDUTName(ChId, Name);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetDUTNames(std::vector<std::string> DUTNames)
{
    m_Base->LGADBase::SetDUTNames(DUTNames);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetDUTBoard(int ChId, AqBoard Brd)
{
    m_Base->LGADBase::SetDUTBoard(ChId, Brd);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetDUTBoards(std::vector<AqBoard> Brds)
{
    m_Base->LGADBase::SetDUTBoards(Brds);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetDUTransImp(int ChId, float Transimp)
{
    m_Base->LGADBase::SetDUTransImp(ChId, Transimp);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetDUTransImps(std::vector<float> Transimps)
{
    m_Base->LGADBase::SetDUTransImps(Transimps);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetDUTSecStage(int ChId, SecStage Amp)
{
    m_Base->LGADBase::SetDUTSecStage(ChId, Amp);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetDUTSecStages(std::vector<SecStage> Amps)
{
    m_Base->LGADBase::SetDUTSecStages(Amps);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetDUTAmpGain(int ChId, int gain)
{
    m_Base->LGADBase::SetDUTAmpGain(ChId, gain);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetDUTAmpGains(std::vector<int> Gains)
{
    m_Base->LGADBase::SetDUTAmpGains(Gains);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetDUTrigg(int ChId, float trigg)
{
    m_Base->LGADBase::SetDUTrigg(ChId, trigg);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetDUTriggs(std::vector<float> Triggs)
{
    m_Base->LGADBase::SetDUTriggs(Triggs);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetDUTFract(int ChId, float frac)
{
    m_Base->LGADBase::SetDUTFract(ChId, frac);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetDUTFracts(std::vector<float> Fracs)
{
    m_Base->LGADBase::SetDUTFracts(Fracs);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetDUTCap(int ChId, double cap)
{
    m_Base->LGADBase::SetDUTCap(ChId, cap);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetDUTCaps(std::vector<double> Caps)
{
    m_Base->LGADBase::SetDUTCaps(Caps);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetChVoltCut(int ChId, double VoltCut)
{
    m_Base->LGADBase::SetChMag(ChId, VoltCut, 0);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetChVoltCuts(std::vector<double> VoltCuts)
{
    m_Base->LGADBase::SetChMags(VoltCuts, 0);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetChNoiseCut(int ChId, double NoiseCut)
{
    m_Base->LGADBase::SetChMag(ChId, NoiseCut, 1);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetChNoiseCuts(std::vector<double> NoiseCuts)
{
    m_Base->LGADBase::SetChMags(NoiseCuts, 1);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetChChargeCut(int ChId, double ChargeCut)
{
    m_Base->LGADBase::SetChMag(ChId, ChargeCut, 2);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetChChargeCuts(std::vector<double> ChargeCuts)
{
    m_Base->LGADBase::SetChMags(ChargeCuts, 2);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetChJitterCut(int ChId, double JitterCut)
{
    m_Base->LGADBase::SetChMag(ChId, JitterCut, 3);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetChJitterCuts(std::vector<double> JitterCuts)
{
    m_Base->LGADBase::SetChMags(JitterCuts, 3);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetPlaneDT(int PlaneId1, int PlaneId2, double DT)
{
    m_Base->LGADBase::SetPlaneDMgt(PlaneId1, PlaneId2, DT, 0);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetPlaneDTs(std::vector<double> DTs)
{
    m_Base->LGADBase::SetPlaneDMgts(DTs, 0);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetPlaneDCrg(int PlaneId1, int PlaneId2, double DCrg)
{
    m_Base->LGADBase::SetPlaneDMgt(PlaneId1, PlaneId2, DCrg, 1);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetPlaneDCrgs(std::vector<double> DCrgs)
{
    m_Base->LGADBase::SetPlaneDMgts(DCrgs, 1);
}
// --------------------------------------------------------------------------------------------------------------
bool LGADUtils::Analyse(Long64_t nprocess, std::string tree)
{
#ifdef _WIN32
    char print[MAX_PATH];
#else
    char print[PATH_MAX];
#endif

    int setup = 0;

    if (m_Base->LGADBase::GetExtention() == ".root" || m_Base->LGADBase::GetConvertSucess())
       {
        bool file = false;
        if (!m_Base->LGADBase::GetConvertSucess()) file = m_Base->SetRootFile(m_Base->GetDataDir() + m_Base->GetDataName() + m_Base->GetExtention());
        else file = m_Base->SetRootFile(m_Base->GetOutFileaDir() + m_Base->GetOutFileName() + ".root");
        if (file)
           {
            if (tree == "")
               {
                if (m_Base->SetRootTree(m_Base->GetRootFile(), m_Base->GetTreeName()))
                   {
                    m_chain = new TChain(m_Base->GetRootTree()->GetName(), "");
                    setup = 1;
                   }
                else {
                      std::cout << __FUNCTION__ << " ERROR: Failed to find tree or more than one trees present!!" << std::endl;
                      return false;
                     }
               }
            else {
                  strcpy(print, tree.c_str());
                  if (m_Base->SetRootTree(m_Base->GetRootFile(), std::string(print)))
                     {
                      m_chain = new TChain(print, "");
                      setup = 2;
                     }
                  else {
                        std::cout << __FUNCTION__ << " ERROR: Failed to find tree " << std::string(print) << std::endl;
                        return false;
                       }
                 }
           }
        else {
              std::cout << __FUNCTION__ << " ERROR: Failed to open file: ";
              if (!m_Base->LGADBase::GetConvertSucess()) std::cout << (m_Base->GetDataDir() + m_Base->GetDataName() + m_Base->GetExtention()).Data() << "!" << std::endl;
              else std::cout << m_Base->GetRootFile()->GetName() << "!" << std::endl;
              return false;
             }
        }
    else {
          std::cout << __FUNCTION__ << " ERROR: No input root file set!!"<< std::endl;
          return false;
         }
   
    if (setup == 1 || setup == 2)
      {
       m_chain->SetCacheSize(500 * 1024 * 1024);
       m_chain->Add((const char*)(m_Base->GetRootFile()->GetName()));
       m_Run = new LGADRun(m_Base);
       if (nprocess != 0 && nprocess < m_chain->GetEntries()) m_chain->Process(m_Run, "", nprocess);
       else m_chain->Process(m_Run);
       if (setup == 1) m_Base->LGADBase::SetConvertSucess(false);
      }

    return true;
}
