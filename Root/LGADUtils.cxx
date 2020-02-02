/*
* LGADUtils.cxx
*
*
*      Author: Gkougkousis Evangelos - Leonidas
*              egkougko@cern.ch
*               IFAE-BARCELONA
*/

#include "LGADUtils/LGADUtils.h"
#include "LinkDef.h"

//#if !defined(__CINT__)
//ClassImp(LGADUtils);
//#endif

// Default constructor creating a general base class object
LGADUtils::LGADUtils()
{
    m_Base = new LGADBase();
}
//------------------------------------------------------------------------------
// Default destructur deleting the general base class object
LGADUtils::~LGADUtils()
{
    delete m_Base;
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
bool LGADUtils::SetSRate(double rate, unsigned int ch)
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
void LGADUtils::SetInstrument(unsigned int instr)
{
    m_Base->LGADBase::SetInstrument(instr);
}
// --------------------------------------------------------------------------------------------------------------
void LGADUtils::SetStartStopEvnt(unsigned int Evnt1, unsigned int Evnt2)
{
    m_Base->LGADBase::SetStartStopEvnt(Evnt1, Evnt2);
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
bool LGADUtils::Analyse(Long64_t nprocess, std::string tree)
{
#ifdef _WIN32
	char print[MAX_PATH];
#else
	char print[PATH_MAX];
#endif

    if (m_Base->LGADBase::GetExtention() == ".root" || m_Base->LGADBase::GetConvertSucess())
       {
        bool file = false;
        if (!m_Base->LGADBase::GetConvertSucess()) file = m_Base->SetRootFile(m_Base->GetDataDir() + m_Base->GetDataName() + m_Base->GetExtention());
        else file = m_Base->SetRootFile(m_Base->GetOutFileaDir() + m_Base->GetOutFileName() + ".root");
        if (file)
           {
            if (tree == "")
               {
                if (m_Base->SetRootTree(m_Base->GetRootFile()))
                   {
                    m_chain = new TChain(m_Base->GetRootTree()->GetName(), "");
                    m_chain->SetCacheSize(500 * 1024 * 1024);
                    m_chain->Add((const char* )(m_Base->GetRootFile()->GetName()));
                    m_Run = new LGADRun(m_Base);
                    if (nprocess != 0 && nprocess < (m_chain->GetTree()->GetEntriesFast())) m_chain->Process(m_Run,"", nprocess);
                    else m_chain->Process(m_Run);
                    m_Base->LGADBase::SetConvertSucess(false);
                   }
                else {
                      std::cout << __FUNCTION__ << " ERROR: Failed to find tree or more than one trees present!!"<< std::endl;
                      return false;
                     }
               }
            else {
                  strcpy(print, tree.c_str());
                  m_chain = new TChain(print, "");
                  m_chain->SetCacheSize(500 * 1024 * 1024);
                  m_chain->Add((const char* )(m_Base->GetRootFile()->GetName()));
                  m_Run = new LGADRun(m_Base);
                  if (nprocess != 0 && nprocess < (m_chain->GetTree()->GetEntriesFast())) m_chain->Process(m_Run,"", nprocess);
                  else m_chain->Process(m_Run);
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

    return true;
}