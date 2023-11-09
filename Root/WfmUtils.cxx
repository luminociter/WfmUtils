/*
* WfmUtils.cxx
*
*
*      Author: Gkougkousis Evangelos - Leonidas
*                egkougko@cern.ch
*            Imperial Collge - London
*/

#include "../WfmUtils/WfmUtils.h"
#include "LinkDef.h"

//#if !defined(__CINT__)
//ClassImp(WfmUtils);
//#endif
//
//#if !defined(__CLING__)
//ClassImp(WfmUtils);
//#endif
//
//#ifdef __CINT__
//#pragma link C++ class WfmUtils;
//#endif
//
//#ifdef __ROOTCLING__
//#pragma link C++ class WfmUtils;
//#endif

// Default constructor creating a general base class object
WfmUtils::WfmUtils()
{
    m_Base = new WfmBase();
}
// --------------------------------------------------------------------------------------------------------------
// Default destructur deleting the general base class object
WfmUtils::~WfmUtils()
{
    delete m_Base;
    if (m_Run) delete m_Run;
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetTrackComb(bool comb)
{
    m_Base->WfmBase::SetTrackComb(comb);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetDoFFT(bool dofft)
{
    m_Base->WfmBase::SetDoFFT(dofft);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetTrckEff(bool TrckEff)
{
    m_Base->WfmBase::SetTrckEff(TrckEff);
}
// --------------------------------------------------------------------------------------------------------------
bool WfmUtils::SetSRate(Long64_t rate, unsigned int ch)
{
    return m_Base->WfmBase::SetSRate(rate, ch);
}
// --------------------------------------------------------------------------------------------------------------
bool WfmUtils::SetNPoints(unsigned int points, unsigned int ch)
{
    return m_Base->WfmBase::SetNPoints(points, ch);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetScopeDelay(std::vector <double> delay)
{
    m_Base->WfmBase::SetScopeDelay(delay);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetTrigClk(double clk)
{
    m_Base->WfmBase::SetTrigClk(clk);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetFitMethode(std::string method)
{
    m_Base->WfmBase::SetFitMethode(method);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetInstrument(AqInstrument instr)
{
    if (instr == TestBeamBin1 || instr == TestBeamBin2) instr = TestBeamBin;
    m_Base->WfmBase::SetInstrument(instr);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetStartStopEvnt(int Evnt1, int Evnt2)
{
    m_Base->WfmBase::SetStartStopEvnt(Evnt1, Evnt2);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetTestEvn(int Evn, int file)
{
    std::pair<int, int> testEvn = std::make_pair(Evn, file);
    m_Base->WfmBase::SetTestEvn(testEvn);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::Initialize()
{
    m_Base->WfmBase::Initialize();
}
// --------------------------------------------------------------------------------------------------------------
bool WfmUtils::ConvertData()
{
    return m_Base->WfmBase::ConvertData();
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetInDataNames(TString DataDir, TString DataName, TString ext)
{
    m_Base->WfmBase::SetInDataNames(DataDir, DataName, ext);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetOutDataNames(TString DataDir, TString DataName)
{
    m_Base->WfmBase::SetOutDataNames(DataDir, DataName);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetTrackInDataNames(TString DataDir, TString DataName)
{
    m_Base->WfmBase::SetTrackInDataNames(DataDir, DataName);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetVerbose(int verbose)
{
    m_Base->WfmBase::SetVerbose(verbose);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetWaveShape(bool shape)
{
    m_Base->WfmBase::SetWaveShape(shape);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetTransFileName(TString filename)
{
    m_Base->WfmBase::SetTransFileName(filename);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetDoTrnsCorr(bool TrnsCorr)
{
    m_Base->WfmBase::SetDoTrnsCorr(TrnsCorr);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetTreeName(std::string treename)
{
    m_Base->WfmBase::SetTreeName(treename);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetExcludeTrackFiles(std::vector<unsigned int> trackExclude)
{
    m_Base->WfmBase::SetExcludeTrackFiles(trackExclude);
}
// --------------------------------------------------------------------------------------------------------------
// #################################### Channel Properties setting functions ####################################
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetDUTName(int ChId, std::string Name)
{
    m_Base->WfmBase::SetDUTName(ChId, Name);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetDUTNames(std::vector<std::string> DUTNames)
{
    m_Base->WfmBase::SetDUTNames(DUTNames);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetDUTBoard(int ChId, AqBoard Brd)
{
    m_Base->WfmBase::SetDUTBoard(ChId, Brd);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetDUTBoards(std::vector<AqBoard> Brds)
{
    m_Base->WfmBase::SetDUTBoards(Brds);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetDUTransImp(int ChId, float Transimp)
{
    m_Base->WfmBase::SetDUTransImp(ChId, Transimp);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetDUTransImps(std::vector<float> Transimps)
{
    m_Base->WfmBase::SetDUTransImps(Transimps);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetDUTSecStage(int ChId, SecStage Amp)
{
    m_Base->WfmBase::SetDUTSecStage(ChId, Amp);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetDUTSecStages(std::vector<SecStage> Amps)
{
    m_Base->WfmBase::SetDUTSecStages(Amps);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetDUTAmpGain(int ChId, int gain)
{
    m_Base->WfmBase::SetDUTAmpGain(ChId, gain);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetDUTAmpGains(std::vector<int> Gains)
{
    m_Base->WfmBase::SetDUTAmpGains(Gains);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetDUTrigg(int ChId, float trigg)
{
    m_Base->WfmBase::SetDUTrigg(ChId, trigg);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetDUTriggs(std::vector<float> Triggs)
{
    m_Base->WfmBase::SetDUTriggs(Triggs);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetDUTFract(int ChId, float frac)
{
    m_Base->WfmBase::SetDUTFract(ChId, frac);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetDUTFracts(std::vector<float> Fracs)
{
    m_Base->WfmBase::SetDUTFracts(Fracs);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetDUTCap(int ChId, double cap)
{
    m_Base->WfmBase::SetDUTCap(ChId, cap);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetDUTCaps(std::vector<double> Caps)
{
    m_Base->WfmBase::SetDUTCaps(Caps);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetDUTUpLim(int ChId, double UpLimX, double UpLimY) 
{
    m_Base->WfmBase::SetDUTUpLim(ChId, std::make_pair(UpLimX, UpLimY));
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetDUTUpLims(std::vector<std::pair<double, double>> upLims)
{
    m_Base->WfmBase::SetDUTUpLims(upLims);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetDUTLowLim(int ChId, double LowLimX, double LowLimY)
{
    m_Base->WfmBase::SetDUTLowLim(ChId, std::make_pair(LowLimX, LowLimY));
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetDUTLowLims(std::vector<std::pair<double, double>> lowLims)
{
    m_Base->WfmBase::SetDUTLowLims(lowLims);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetDUTBinSize(int ChId, double size)
{
    m_Base->WfmBase::SetDUTBinSize(ChId, size);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetDUTBinSizes(std::vector<double> sizes)
{
    m_Base->WfmBase::SetDUTBinSizes(sizes);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetDUTBinShape(int ChId, std::string shape)
{
    m_Base->WfmBase::SetDUTBinShape(ChId, shape);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetDUTBinShapes(std::vector<std::string> shapes)
{
    m_Base->WfmBase::SetDUTBinShapes(shapes);
}
// --------------------------------------------------------------------------------------------------------------
// ################################### Channel related cuts setting functions ###################################
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetChVoltCut(int ChId, double VoltCut)
{
    m_Base->WfmBase::SetChMag(ChId, VoltCut, 0);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetChVoltCuts(std::vector<double> VoltCuts)
{
    m_Base->WfmBase::SetChMags(VoltCuts, 0);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetChNoiseCut(int ChId, double NoiseCut)
{
    m_Base->WfmBase::SetChMag(ChId, NoiseCut, 1);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetChNoiseCuts(std::vector<double> NoiseCuts)
{
    m_Base->WfmBase::SetChMags(NoiseCuts, 1);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetChChargeCut(int ChId, double ChargeCut)
{
    m_Base->WfmBase::SetChMag(ChId, ChargeCut, 2);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetChChargeCuts(std::vector<double> ChargeCuts)
{
    m_Base->WfmBase::SetChMags(ChargeCuts, 2);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetChJitterCut(int ChId, double JitterCut)
{
    m_Base->WfmBase::SetChMag(ChId, JitterCut, 3);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetChJitterCuts(std::vector<double> JitterCuts)
{
    m_Base->WfmBase::SetChMags(JitterCuts, 3);
}
// --------------------------------------------------------------------------------------------------------------
// ################################## Inter-Channel cuts fpr timing and charge ##################################
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetPlaneDT(int PlaneId1, int PlaneId2, double DT)
{
    m_Base->WfmBase::SetPlaneDMgt(PlaneId1, PlaneId2, DT, 0);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetPlaneDTs(std::vector<double> DTs)
{
    m_Base->WfmBase::SetPlaneDMgts(DTs, 0);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetPlaneDCrg(int PlaneId1, int PlaneId2, double DCrg)
{
    m_Base->WfmBase::SetPlaneDMgt(PlaneId1, PlaneId2, DCrg, 1);
}
// --------------------------------------------------------------------------------------------------------------
void WfmUtils::SetPlaneDCrgs(std::vector<double> DCrgs)
{
    m_Base->WfmBase::SetPlaneDMgts(DCrgs, 1);
}
// --------------------------------------------------------------------------------------------------------------
bool WfmUtils::Analyse(unsigned int stage, std::string tree)
{

    int file1 = -1;
    bool file2 = false;
    std::size_t found = 0;
    std::vector<TString> rootfiles;
    std::vector<unsigned int> fileindx;
    // m_convert takes the return value of the data conversion function and m_datanames is cleaned at the end of the data conereter
    if (m_Base->WfmBase::GetConvertSucess() && (m_Base->WfmBase::GetInFileNames()).size() == 0) // finsihed conversion but did not reset the input files        
       { 
        file1 = m_Base->WfmBase::SetInRootFile(m_Base->WfmBase::GetOutRootFile());
        file2 = m_Base->WfmBase::SetOutRootFile(m_Base->WfmBase::GetOutFileDir() + m_Base->WfmBase::GetOutFileName() + ".root");
       }
    else if ((m_Base->WfmBase::GetInFileNames()).size() != 0)
            {
             rootfiles.clear();
             for (unsigned int k = 0; k < (m_Base->WfmBase::GetInFileNames()).size(); k++)
                 { 
                  rootfiles.push_back((m_Base->GetInFileNames()).at(k));
                  if (m_Base->WfmBase::GetExtention() == "")
                     { 
                      if (!((rootfiles.back()).Contains("."))) { rootfiles.pop_back(); continue; }
                      else found = (rootfiles.back()).Last('.');
                      if (rootfiles.back()(found + 1, (rootfiles.back()).Length()) != "root")
                         {
                          rootfiles.pop_back(); 
                          continue;
                         }
                      else rootfiles.back() = rootfiles.back()(0, found);
                     }
                 }
             if (m_Base->WfmBase::GetExtention() == "") m_Base->WfmBase::SetExtention("root");
             file1 = m_Base->WfmBase::SetInRootFiles(rootfiles, fileindx);
             // Append on the input file
             if ((m_Base->WfmBase::GetOutFileName() == "" || (m_Base->WfmBase::GetOutFileName() == (m_Base->WfmBase::GetInFileNames()).at(0)) && (m_Base->WfmBase::GetOutFileDir() == m_Base->WfmBase::GetInDataDir())))
                {
                 file2 = true;
                 m_Base->WfmBase::SetOutRootFile((m_Base->WfmBase::GetInRootFiles()).at(0));
                }
             // Create a new file for the analysis output
             else file2 = m_Base->WfmBase::SetOutRootFile(m_Base->WfmBase::GetOutFileDir() + m_Base->WfmBase::GetOutFileName() + ".root");
            }
    else {
          std::cout << __FUNCTION__ << " ERROR: No input root file set!!" << std::endl;
          return false;
         }
        
    // Block to determine the name of the data tree to analyze for the input files. Normal is wfm
    if (file1 > 0 && file2)
       {
        if (tree == "")
           {
            if (!(m_Base->WfmBase::SetRootTrees(m_Base->WfmBase::GetInRootFiles(), fileindx, m_Base->WfmBase::GetTreeName())))
               {
                std::cout << __FUNCTION__ << " ERROR: Failed to find tree or more than one trees present!!" << std::endl;
                return false;
               }
           }
        else {
              if (!(m_Base->WfmBase::SetRootTrees(m_Base->WfmBase::GetInRootFiles(), fileindx, tree)))
                 {
                  std::cout << __FUNCTION__ << " ERROR: Failed to find tree " << tree << std::endl;
                  return false;
                 }
             }
       }
    else {
          if (file1 <= 0) std::cout << __FUNCTION__ << " ERROR: Failed to open file: " << m_Base->GetInDataDir() + (m_Base->GetInFileNames()).at(0) + "." + m_Base->GetExtention() << "!" << std::endl;
          if (!file2) std::cout << __FUNCTION__ << " ERROR: Output root file not correctly set!!" << std::endl;
          return false;
         }

    if (stage > 0 && stage < 6) m_Base->WfmBase::SetAnaStage(stage);
    else {
          std::cout << __FUNCTION__ << " WARNING: Analysis level invalid, (" << stage << "), setting to first level..." << std::endl;
          m_Base->WfmBase::SetAnaStage(1);
         }

    if (tree == "") m_Base->WfmBase::SetConvertSucess(false);
    m_Run = new WfmRun(m_Base);
    // Input root files are already open at this point, the m_infiles and m_trees vectors should be populated with the appropriate pointers which can be immediately used
    m_Run->WfmRun::InitAnalysis(fileindx);

    return true;
}
// --------------------------------------------------------------------------------------------------------------