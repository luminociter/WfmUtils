/*
* LGADUtils.h
*
*
*      Author: Gkougkousis Evangelos - Leonidas
*              egkougko@cern.ch
*               IFAE-BARCELONA
*/

#ifndef __LGADUtils__
#define __LGADUtils__

#include "LGADUtils/LGADRun.h"
#include "TChain.h"

class LGADUtils : public LGADRun {

public:

   LGADUtils();
   ~LGADUtils();

   void SetTrackComb(bool comb);
   void SetFEi4Eff(bool FEi4Eff);
   bool SetSRate(double rate, unsigned int ch);
   bool SetNPoints(unsigned int points, unsigned int ch);
   void Initialize();
   void SetFitMethode(std::string method);
   void SetInstrument(unsigned int instr);
   void SetStartStopEvnt(int Evnt1 = 0, int Evnt2 = 0);
   bool ConvertData();
   void SetTreeName(std::string treename);
   void SetInDataNames(TString DataDir = "", TString DataName = "", TString ext = "");
   void SetOutDataNames(TString DataDir = "", TString DataName = "");
   void SetTrackInDataNames(TString DataDir = "", TString DataName = "");
   void SetVerbose(int verbose);
   void SetWaveShape(bool shape = false);
   void SetDoTrnsCorr(bool TrnsCorr = false);
   void SetTransFileName(TString filename = "");
   bool Analyse(Long64_t nprocess = 0, std::string tree = "");
   std::string GetTreeName() { return LGADBase::GetTreeName(); };

private:

   LGADRun *m_Run;
   TChain *m_chain;
   LGADBase *m_Base;

};

#endif