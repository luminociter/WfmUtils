/*
* WfmUtils.h
*
*
*      Author: Gkougkousis Evangelos - Leonidas
*                egkougko@cern.ch
*            Imperial Collge - London
*/

#ifndef __WfmUtils__
#define __WfmUtils__

#include "WfmRun.h"

class WfmUtils : public WfmRun {

public:

   WfmUtils();
   ~WfmUtils();

   void SetTrackComb(bool comb);
   void SetExcludeTrackFiles(std::vector<unsigned int> trackExclude);
   void SetDoFFT(bool dofft);
   void SetTrckEff(bool TrckEff);
   bool SetSRate(Long64_t rate, unsigned int ch);
   bool SetNPoints(unsigned int points, unsigned int ch);
   void SetScopeDelay(std::vector<double> delay);
   void SetTrigClk(double clk);
   void Initialize();
   void SetFitMethode(std::string method);
   void SetInstrument(AqInstrument instr);
   void SetStartStopEvnt(int Evnt1 = 0, int Evnt2 = 0);
   void SetTestEvn(int Evn = 1, int file = -1);
   bool ConvertData();
   void SetTreeName(std::string treename);
   void SetInDataNames(TString DataDir = "", TString DataName = "", TString ext = "");
   void SetOutDataNames(TString DataDir = "", TString DataName = "");
   void SetTrackInDataNames(TString DataDir = "", TString DataName = "");
   void SetVerbose(int verbose);
   void SetWaveShape(bool shape = false);
   void SetDoTrnsCorr(bool TrnsCorr = false);
   void SetTransFileName(TString filename = "");

   // DUT level settings (have to be subdivided to cahnnel and plane level)
   void SetDUTName(int ChId, std::string Name);
   void SetDUTNames(std::vector<std::string> DUTNames);
   void SetDUTBoard(int ChId, AqBoard Brd);
   void SetDUTBoards(std::vector<AqBoard> Brds);
   void SetDUTransImp(int ChId, float Transimp);
   void SetDUTransImps(std::vector<float> Transimps);
   void SetDUTSecStage(int ChId, SecStage Amp);
   void SetDUTSecStages(std::vector<SecStage> Amps);
   void SetDUTAmpGain(int ChId, int gain);
   void SetDUTAmpGains(std::vector<int> Gains);
   void SetDUTrigg(int ChId, float trigg);
   void SetDUTriggs(std::vector<float> Triggs);
   void SetDUTFract(int ChId, float frac);
   void SetDUTFracts(std::vector<float> Fracs);
   void SetDUTCap(int ChId, double cap);
   void SetDUTCaps(std::vector<double> Caps);
   void SetDUTUpLim(int ChId, double UpLimX, double UpLimY);
   void SetDUTUpLims(std::vector<std::pair<double, double>> upLims);
   void SetDUTLowLim(int ChId, double LowLimX, double LowLimY);
   void SetDUTLowLims(std::vector<std::pair<double, double>> lowLims);
   void SetDUTBinSize(int ChId, double size);
   void SetDUTBinSizes(std::vector<double> sizes);
   void SetDUTBinShape(int ChId, std::string shape);
   void SetDUTBinShapes(std::vector<std::string> shapes);

   // Channel level cuts
   void SetChVoltCut(int ChId, double VoltCut);
   void SetChVoltCuts(std::vector<double> VoltCuts);
   void SetChNoiseCut(int ChId, double NoiseCut);
   void SetChNoiseCuts(std::vector<double> NoiseCuts);
   void SetChChargeCut(int ChId, double ChargeCut);
   void SetChChargeCuts(std::vector<double> ChargeCuts);
   void SetChJitterCut(int ChId, double JitterCut);
   void SetChJitterCuts(std::vector<double> JitterCuts);

   // Plane level cuts
   void SetPlaneDT(int PlaneId1, int PlaneId2, double DT);
   void SetPlaneDTs(std::vector<double> DTs);
   void SetPlaneDCrg(int PlaneId1, int PlaneId2, double DCrg);
   void SetPlaneDCrgs(std::vector<double> DCrgs);

   bool Analyse(unsigned int stage = 1, std::string tree = "");
   std::string GetTreeName() { return WfmBase::GetTreeName(); };

private:

   WfmRun *m_Run;
   WfmBase *m_Base;

   // ClassDef(WfmUtls, 1);
};

#endif