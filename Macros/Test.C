 {
    gInterpreter->Reset(); // Reset Root

    gInterpreter->LoadMacro("C:/Users/Vagelis/Desktop/LGADUtils/LGADUtils/Root/LGADBase.cxx++");
    gInterpreter->LoadMacro("C:/Users/Vagelis/Desktop/LGADUtils/LGADUtils/Root/LGADFits.cxx++");
    gInterpreter->LoadMacro("C:/Users/Vagelis/Desktop/LGADUtils/LGADUtils/Root/TrackCombine.cxx++");
    gInterpreter->LoadMacro("C:/Users/Vagelis/Desktop/LGADUtils/LGADUtils/Root/DataConverters.cxx++");
    gInterpreter->LoadMacro("C:/Users/Vagelis/Desktop/LGADUtils/LGADUtils/Root/WaveForm.cxx++");
    gInterpreter->LoadMacro("C:/Users/Vagelis/Desktop/LGADUtils/LGADUtils/Root/DUTChannel.cxx++");
    gInterpreter->LoadMacro("C:/Users/Vagelis/Desktop/LGADUtils/LGADUtils/Root/IOHandler.cxx++");
    gInterpreter->LoadMacro("C:/Users/Vagelis/Desktop/LGADUtils/LGADUtils/Root/LGADRun.cxx++");
    gInterpreter->LoadMacro("C:/Users/Vagelis/Desktop/LGADUtils/LGADUtils/Root/LGADUtils.cxx++");

    gROOT->ProcessLine("LGADUtils* Utils = new LGADUtils()");
    gROOT->ProcessLine("Utils->SetVerbose(0)");
    gEnv->SetValue("TFile.AsyncPrefetching", 1);
/*
   // Basic analysis options
   // For the instrument type, use the AqInstrument type definition values as described at the forward declaration section of LGADBase
   // AqInstrument: { Sampic, LabTXT, TestBeamBin, TektronixScope, LeCroyWRBin, Unasigned }
    gROOT->ProcessLine("Utils->Initialize()"); // Reset the analysis object if you want to run a second analysis within the same session
    gROOT->ProcessLine("Utils->SetFitMethode(std::string method)"); // Prefered methode for the itterative reffiter, allowed values: RooFit, root, rootInt
    gROOT->ProcessLine("Utils->SetInstrument(AqInstrument instr)"); // Set the instrument used for data taking. See type definition
    gROOT->ProcessLine("Utils->ConvertData();");
    gROOT->ProcessLine("Utils->Analyse()");
    gROOT->ProcessLine("Utils->SetStartStopEvnt(int Evnt1, int Evnt2)");
    gROOT->ProcessLine("Utils->SetTreeName(std::string treename)"); // Root tree name where converted data are store if its not the default wfm
    gROOT->ProcessLine("Utils->SetInDataNames(TString DataDir, TString DataName, TString ext)");
    gROOT->ProcessLine("Utils->SetOutDataNames(TString DataDir, TString DataName)");
    gROOT->ProcessLine("Utils->SetVerbose(int verbose)"); // Code verbosity level, rnaging from 0 - 3
    gROOT->ProcessLine("Utils->SetWaveShape(bool shape)"); // Enable the wavefornm shape sit methode for more accurate resutls (computationaly intensive)
    gROOT->ProcessLine("Utils->SeDoFFT(bool dofft)"); // Enable or disable the per evetn FFT and overall FFT fits for the run
    gROOT->ProcessLine("Utils->SetTestEvn(int EVN)"); // Test only a specific event in the analysis (1st stage) 

    // Transimpedence correction with respect to signal frquency when available simulation data
    gROOT->ProcessLine("Utils->SetDoTrnsCorr(bool TrnsCorr)");
    gROOT->ProcessLine("Utils->SetTransFileName(TString filename)");

    //TestBeam Related Options
    gROOT->ProcessLine("Utils->SetTrackComb(bool comb)");
    gROOT->ProcessLine("Utils->SetFEi4Eff(bool FEi4Eff)");
    gROOT->ProcessLine("Utils->SetTrackPackage(std::string package)");
    gROOT->ProcessLine("Utils->SetScopeDelay(std::vector<double> delay)"); // Delay between the ROI trigger and the oscilloscope trigger. THis is intially an empty vector, initialized to 
                                                                           // default valus 50e-9 and 162e-9. In case of more oscilloscopes, additional instruments are intialized to -99.
    gROOT->ProcessLine("Utils->SetTrigClk(double clk"); // Trigger clock of the ROI. Dfault value 25 nsec initialized only when processing testbeam instruments and only
                                                        // has an effect when used in conjuction with TrackComb variable. 
    gROOT->ProcessLine("Utils->SetTrackInDataNames(TString DataDir, TString DataName)");

    // When seeting channel lelvel propertires use the correct definition of the enumerators for the secind stage amplifier and the 
    // aquisition board as the y are defined in LGADBase.h in the form of forward declarations. Implement values:
    // AqBoard: { SingleCh, FourCh, IN2P3, KU, SiPM }
    // SecStage: { MinCircuits, Particulars, CIVIDEC, none }
    // Options for setting channel level properties
    gROOT->ProcessLine("Utils->SetSRate(Long64_t rate, unsigned int ch)");
    gROOT->ProcessLine("Utils->SetNPoints(unsigned int points, unsigned int ch)");
    gROOT->ProcessLine("Utils->SetDUTName(int ChId, std::string Name)");
    gROOT->ProcessLine("Utils->SetDUTNames(std::vector<std::string> DUTNames)");
    gROOT->ProcessLine("Utils->SetDUTBoard(int ChId, AqBoard Brd)");
    gROOT->ProcessLine("Utils->SetDUTBoards(std::vector<AqBoard> Brds)");
    gROOT->ProcessLine("Utils->SetDUTransImp(int ChId, float Transimp)");
    gROOT->ProcessLine("Utils->SetDUTransImps(std::vector<float> Transimps)");
    gROOT->ProcessLine("Utils->SetDUTSecStage(int ChId, SecStage Amp)");
    gROOT->ProcessLine("Utils->SetDUTSecStages(std::vector<SecStage> Amps)");
    gROOT->ProcessLine("Utils->SetDUTAmpGain(int ChId, int gain)");
    gROOT->ProcessLine("Utils->SetDUTAmpGains(std::vector<int> Gains)");
    gROOT->ProcessLine("Utils->SetDUTrigg(int ChId, float trigg)");
    gROOT->ProcessLine("Utils->SetDUTriggs(std::vector<float> Triggs)");
    gROOT->ProcessLine("Utils->SetDUTFract(int ChId, float frac)");
    gROOT->ProcessLine("Utils->SetDUTFracts(std::vector<float> Fracs)");
    gROOT->ProcessLine("SetDUTCap(int ChId, double cap)"); // Set sensor capacitance in pf for selecting correct transimpedence simulation histogram when appplying fft corrections
    gROOT->ProcessLine("SetDUTCaps(std::vector<double> Caps)"); // Set vector of sensors capacitance for a multitudde of channels. Values adjasted sequentially to all non-set channels

    // Options for setting channel level analyisis cuts 
    gROOT->ProcessLine("Utils->SetChVoltCut(int ChId, double VoltCut)");
    gROOT->ProcessLine("Utils->SetChVoltCuts(std::vector<double> VoltCuts)");
    gROOT->ProcessLine("Utils->SetChNoiseCut(int ChId, double NoiseCut)");
    gROOT->ProcessLine("Utils->SetChNoiseCuts(std::vector<double> NoiseCuts)");
    gROOT->ProcessLine("Utils->SetChChargeCut(int ChId, double ChargeCut)");
    gROOT->ProcessLine("Utils->SetChChargeCuts(std::vector<double> ChargeCuts)");
    gROOT->ProcessLine("Utils->SetChJitterCut(int ChId, double JitterCut)");
    gROOT->ProcessLine("Utils->SetChJitterCuts(std::vector<double> JitterCuts)");

    // Options for setting inter-channel level cuts (soon to be defined as plane level cuts)
    gROOT->ProcessLine("SetPlaneDT(int PlaneId1, int PlaneId2, double DT)");
    gROOT->ProcessLine("SetPlaneDTs(std::vector<double> DTs)");
    gROOT->ProcessLine("SetPlaneDCrg(int PlaneId1, int PlaneId2, double DCrg)");
    gROOT->ProcessLine("SetPlaneDCrgs(std::vector<double> DCrgs)");
*/

    // TestBeam Example
    gROOT->ProcessLine("Utils->SetTrackComb(true)");
    gROOT->ProcessLine("Utils->SetTrigClk(25e-9)");
    gROOT->ProcessLine("Utils->SetScopeDelay({50e-9, 162e-9})");
    gROOT->ProcessLine("Utils->SetTrackInDataNames(\"C:/Users/Vagelis/Desktop/LGADUtils/LGADUtils/SampleDatasets/TestbeamTek/run000426.root\")");
    gROOT->ProcessLine("Utils->SetInstrument(TestBeamBin)");
    gROOT->ProcessLine("Utils->SetStartStopEvnt(0, 500)");
    gROOT->ProcessLine("Utils->SetInDataNames(\"C:/Users/Vagelis/Desktop/LGADUtils/LGADUtils/SampleDatasets/TestbeamTek\", \"data_1657215023\")");
    gROOT->ProcessLine("Utils->SetOutDataNames(\"C:/Users/Vagelis/Desktop/LGADUtils/LGADUtils/SampleDatasets/TestbeamTek\", \"Combined.root\")");
    gROOT->ProcessLine("Utils->ConvertData();");
    // gROOT->ProcessLine("Utils->SetInDataNames(\"C:/Users/Vagelis/Desktop/LGADUtils/LGADUtils/SampleDatasets/TestbeamTek\", \"Converted_Reducted.root\")");
    // gROOT->ProcessLine("Utils->SetOutDataNames(\"C:/Users/Vagelis/Desktop/LGADUtils/LGADUtils/SampleDatasets/TestbeamTek\", \"Analyzsed1stStageTest.root\")");
    // gROOT->ProcessLine("Utils->SetWaveShape(false)");
    // gROOT->ProcessLine("Utils->SetTestEvn(496)");
    // gROOT->ProcessLine("Utils->Analyse(1)");
 
    // Sampic Example
    // gROOT->ProcessLine("Utils->SetInstrument(Sampic)");
    // gROOT->ProcessLine("Utils->SetInDataNames(\"C:/Users/Vagelis/Desktop/HPK_-10C_-90V_Sampic\")");
    // gROOT->ProcessLine("Utils->SetOutDataNames(\"C:/Users/Vagelis/Desktop/HPK_-10C_-90V_Sampic\")");
    // gROOT->ProcessLine("Utils->ConvertData();");
}