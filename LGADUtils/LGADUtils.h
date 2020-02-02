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
#include "LGADUtils/LGADBase.h"
#include "TChain.h"

class LGADUtils : public LGADRun {

public:

   LGADUtils();
   ~LGADUtils();

   void SetTrackComb(bool comb);
   void SetFEi4Eff(bool FEi4Eff);
   bool SetSRate(double rate, unsigned int ch);
   bool SetNPoints(unsigned int points, unsigned int ch);
   void SetFitMethode(std::string method);
   void SetInstrument(unsigned int instr);
   void SetStartStopEvnt(unsigned int Evnt1 = 0, unsigned int Evnt2 = 0);
   bool ConvertData();
   void SetInDataNames(TString DataDir = "", TString DataName = "", TString ext = "");
   void SetOutDataNames(TString DataDir = "", TString DataName = "");
   void SetTrackInDataNames(TString DataDir = "", TString DataName = "");
   void SetVerbose(int verbose);
   void SetWaveShape(bool shape = false);
   bool Analyse(Long64_t nprocess = 0, std::string tree = "");

private:

   LGADBase* m_Base;
   LGADRun *m_Run;
   TChain *m_chain;

};

#endif