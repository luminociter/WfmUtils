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
void LGADUtils::SetDoFFT(bool dofft)
{
    m_Base->LGADBase::SetDoFFT(dofft);
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
void LGADUtils::SetScopeDelay(std::vector <double> delay)
{
    m_Base->LGADBase::SetScopeDelay(delay);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetTrigClk(double clk)
{
    m_Base->LGADBase::SetTrigClk(clk);
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
void LGADUtils::SetTestEvn(int Evn)
{
    m_Base->LGADBase::SetTestEvn(Evn);
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
void LGADUtils::SetExcludeTrackFiles(std::vector<unsigned int> trackExclude)
{
    m_Base->LGADBase::SetExcludeTrackFiles(trackExclude);
}
// --------------------------------------------------------------------------------------------------------------
bool LGADUtils::Analyse(unsigned int stage, std::string tree)
{
#ifdef _WIN32
    char print[MAX_PATH];
#else
    char print[PATH_MAX];
#endif

    int file1 = -1;
    bool file2 = false;
    std::size_t found = 0;
    std::vector<TString> rootfiles;
    srd::vector<unsigned int> fileindx;
    // m_convert takes the return value of the data cnersion function and m_datanames is cleanes at the end of the data conereter
    if (m_Base->LGADBase::GetConvertSucess() && (m_Base->GetInFileNames()).size() == 0) // finsihed conversion but did not reset the input files        
       { 
        file1 = m_Base->LGADBase::SetInRootFile(m_Base->GetOutRootFile()->GetName());
        file2 = m_Base->LGADBase::SetOutRootFile(m_Base->LGADBase::GetOutFileDir() + m_Base->LGADBase::GetOutFileName() + ".root");
       }
    else if ((m_Base->GetInFileNames()).size() != 0)
            {
             rootfiles.clear();
             for (unsigned int k = 0; k < (m_Base->GetInFileNames()).size(); k++)
                 { 
                  rootfiles.push_back((m_Base->GetInFileNames()).at(k));
                  if (m_Base->LGADBase::GetExtention() == "")
                     { 
                      if (!(rootfiles.back().Contains("."))) { rootfiles.pop_back(); continue; }
                      else found = (rootfiles.back()).Last('.');
                      if (rootfiles.back()(found + 1, (rootfiles.back()).Length()) != "root")
                         {
                          rootfiles.pop_back(); 
                          continue;
                         }
                      else rootfiles.back() = rootfiles.back()(0, found);
                     }
                 }
             if (m_Base->LGADBase::GetExtention() == "") m_Base->LGADBase::SetExtention("root");
             // Fix the filenames and extentions in a single vector
             for (unsigned int k = 0; k < rootfiles.size(); k++) rootfiles.at() = m_Base->GetInDataDir() + rootfiles.at() "." + m_Base->GetExtention();
             file1 = m_Base->LGADBase::SetInRootFiles(rootfiles, fileindx);
            // Append on the input file
            if ((m_Base->LGADBase::GetOutFileName() == "" || (m_Base->LGADBase::GetOutFileName() == (m_Base->LGADBase::GetInFileNames()).at(0)) && (m_Base->LGADBase::GetOutFileDir() == m_Base->LGADBase::GetInDataDir())))
               {
                file2 = true;
                m_Base->LGADBase::SetOutRootFile((m_Base->GetInRootFiles()).at(0));
               }
            // Create a new file for the analysis output
            else file2 = m_Base->LGADBase::SetOutRootFile(m_Base->LGADBase::GetOutFileDir() + m_Base->LGADBase::GetOutFileName() + ".root");
            }
    else {
          std::cout << __FUNCTION__ << " ERROR: No input root file set!!" << std::endl;
          return false;
         }
        
    // Block to determine the name of the data tree to analyze for the input files. Normal is wfm
    if (file1 != -1 && file2)
       {
        if (tree == "")
           {
            if (!(m_Base->SetRootTree(m_Base->GetInRootFiles(), fileindx, m_Base->GetTreeName())))
               {
                std::cout << __FUNCTION__ << " ERROR: Failed to find tree or more than one trees present!!" << std::endl;
                return false;
               }
           }
        else {
              if (!(m_Base->SetRootTree(m_Base->GetInRootFile(), fileindx, tree)))
                 {
                  std::cout << __FUNCTION__ << " ERROR: Failed to find tree " << tree << std::endl;
                  return false;
                 }
             }
       }
    else {
          if (!file1) std::cout << __FUNCTION__ << " ERROR: Failed to open file: " << m_Base->GetInDataDir() + (m_Base->GetInFileNames()).at(0) + "." + m_Base->GetExtention() << "!" << std::endl;
          if (!file2) std::cout << __FUNCTION__ << " ERROR: Output root file not correctly set!!" << std::endl;
          return false;
         }

    // Setting analysis level
    // Stage 1: Only waveform class run in the event loop 
    // Stage 2: The Channel class is run after the event loop to update channel histograms
    // Stage 3: The CFD maps and chi2 minimizagion is run alognside the inter-channel purity cuts if they are defined
    // Stage 4: The efficiency tracking class is run to produce the 3D timing maps with the CFD selectedcuts
    if (stage > 0 && stage < 5) m_Base->LGADBase::SetAnaStage(stage);
    else {
          std::cout << __FUNCTION__ << " WARNING: Analysis level invalid, (" << stage << "), setting to first level..." << std::endl;
          m_Base->LGADBase::SetAnaStage(1);
         }

    if (tree == "") m_Base->LGADBase::SetConvertSucess(false);
    m_Base->LGADBase::SetAnaStage(stage);
    m_Run = new LGADRun(m_Base);
    // Input root files are already open at this point, the m_infiles and m_trees vectors should be populated with the appropriate pointers which can be immediately used
    m_Run->LGADRun::InitAnalysis(fileindx);

    return true;
}
// --------------------------------------------------------------------------------------------------------------