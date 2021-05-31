{
    gInterpreter->Reset(); // Reset Root

    gInterpreter->AddIncludePath("<path>/LGADUtils/");

    gInterpreter->LoadMacro("<path>/LGADUtils/Root/LGADBase.cxx+");
    gInterpreter->LoadMacro("<path>/LGADUtils/Root/LGADFits.cxx+");
    gInterpreter->LoadMacro("<path>/LGADUtils/Root/DataConverters.cxx+");
    gInterpreter->LoadMacro("<path>/LGADUtils/Root/LGADSel.cxx+");
    gInterpreter->LoadMacro("<path>/LGADUtils/Root/WaveForm.cxx+");
    gInterpreter->LoadMacro("<path>/LGADUtils/Root/DUTChannel.cxx+");
    gInterpreter->LoadMacro("<path>/LGADUtils/Root/LGADRun.cxx+");
    gInterpreter->LoadMacro("<path>/LGADUtils/Root/LGADUtils.cxx+");

    gROOT->ProcessLine("LGADUtils* Utils = new LGADUtils()");
    gROOT->ProcessLine("Utils->SetVerbose(0)");
    gEnv->SetValue("TFile.AsyncPrefetching", 1);

    // Sampic Example
    //gROOT->ProcessLine("Utils->SetInstrument(0)");
    // gROOT->ProcessLine("Utils->SetStartStopEvnt(0, 200)");
    //gROOT->ProcessLine("Utils->SetInDataNames(\"<path>/LGADUtils/Sample Datasets/Run_SAMPIC_Vagelis3_Data_2_26_2018_Ascii\")");
    //gROOT->ProcessLine("Utils->SetOutDataNames(\"<path>/LGADUtils/Sample Datasets\")");

    // Lab example
    gROOT->ProcessLine("Utils->SetInstrument(1)");
    gROOT->ProcessLine("Utils->SetStartStopEvnt(0, 5000)");
    gROOT->ProcessLine("Utils->SetInDataNames(\"<path>/LGADUtils/Sample Datasets/LecroyTXT\")");
    gROOT->ProcessLine("Utils->SetOutDataNames(\"<path>/LGADUtils/Sample Datasets\")");

    // Agilent Binary example
    //gROOT->ProcessLine("Utils->SetInstrument(2)");
    //gROOT->ProcessLine("Utils->SetStartStopEvnt(0, 0)");
    //gROOT->ProcessLine("Utils->SetInDataNames(\"<path>/LGADUtils/Sample Datasets/AgilentBinary\", \"data_1536561052\")");
    //gROOT->ProcessLine("Utils->SetOutDataNames(\"<path>/LGADUtils/Sample Datasets\")");

    // LeCroy Binary example
    //gROOT->ProcessLine("Utils->SetInstrument(4)");
    //gROOT->ProcessLine("Utils->SetStartStopEvnt(0, 0)");
    //gROOT->ProcessLine("Utils->SetInDataNames(\"<path>/LGADUtils/Sample Datasets/LecroyBinary\", \"data_1553874457\")");
    //gROOT->ProcessLine("Utils->SetOutDataNames(\"<path>/LGADUtils/Sample Datasets\")");

    // Options part to set for the run
/*
    gROOT->ProcessLine("Utils->SetTrackComb(true)"); // incluides EUDAT track file
    gROOT->ProcessLine("Utils->SetFEi4Eff(true)");   // Include FEi4 exclusion regions
    gROOT->ProcessLine("Utils->SetSRate(Long64_t rate, unsigned int ch)");   // Set sampling rate for individual channels IF not using autodetect (default)
    gROOT->ProcessLine("Utils->SetNPoints(unsigned int points, unsigned int ch)"); // Set number of points per channel IF not using autodetect (default)
    gROOT->ProcessLine("Utils->Initialize()");       // Reset all the LGADBase Object properties torestart a new run after finishing one 
    gROOT->ProcessLine("Utils->SetFitMethode()");    // Set the methode for the iterative refiiter (Standar Minuit, RooFit, Integral optimization)
    gROOT->ProcessLine("Utils->SetInstrument()");    // Set the instrument data were taken for the conversion
    gROOT->ProcessLine("Utils->ConvertData();");     // Convert waveform data to root format
    gROOT->ProcessLine("Utils->SetTreeName()");      // Set the name of the input tree for waveform data IF not using autodetect (defualt)
    gROOT->ProcessLine("Utils->SetInDataNames()");   // Set the folder and extention for input data
    gROOT->ProcessLine("Utils->SetOutDataNames()");  // Set folder and name where to wrtie output
    gROOT->ProcessLine("Utils->SetTrackInDataNames()"); // Set directory for tracking data file
    gROOT->ProcessLine("Utils->SetVerbose(0)");   // Set Verbosity level for debug
    gROOT->ProcessLine("Utils->SetWaveShape(true)");  // Do waveform Shape calculations 
    gROOT->ProcessLine("Utils->SetDoTrnsCorr(true)");  // Do transimpedence correction vs signal FFT in an event by event basis
    gROOT->ProcessLine("Utils->SetTransFileName()");   // Set the input file containing the transimpedence vs frequency simulation data

    gROOT->ProcessLine("Utils->SetDUTName()");  // Set name for a channel
    gROOT->ProcessLine("Utils->SetDUTNames()"); // Set a vector of names for several channels
    gROOT->ProcessLine("Utils->SetDUTBoard()");  // Set board type for a specific channel (supported UCSC single channel, UCSC four Channel, IN2P3 four channel , Kansas U board, SiPM Infineon board)
    gROOT->ProcessLine("Utils->SetDUTBoards()"); // Set board types as a vector for multiple channels
    gROOT->ProcessLine("Utils->SetDUTransImp()");  // Set transimpedence value IF not using event by event transimpedence correction
    gROOT->ProcessLine("Utils->SetDUTransImps()");  // Set transimpedence values as vector for multiple channels
    gROOT->ProcessLine("Utils->SetDUTSecStage()");  // Set Second stage amplifier, supported are MiniCircuits, Particulars, Cividec, none
    gROOT->ProcessLine("Utils->SetDUTSecStages()");  // Set Second stage amplifiers as a vector for multiple channels
    gROOT->ProcessLine("Utils->SetDUTAmpGain()");    // Set amplifier gain if not using one of th estandard suported amplifieres
    gROOT->ProcessLine("Utils->SetDUTAmpGains()");   // Set amplifiers gains as a vector for multiple channels
    gROOT->ProcessLine("Utils->SetDUTrigg()");   // Set trigger threshold in V for a specific channel, no cut applied
    gROOT->ProcessLine("Utils->SetDUTriggs()");   // Set trigger threshold as vector for multiple channels
    gROOT->ProcessLine("Utils->SetDUTFract()");   // Set standard CFD percentrage for a channel
    gROOT->ProcessLine("Utils->SetDUTFracts()");  // Set CFD percentages as a vector for multiple channels
    gROOT->ProcessLine("Utils->SetChDT()");  // This is a cut (!) specifying the time difference between two chanhels for an event to be considered
    gROOT->ProcessLine("Utils->SetChDTs()");  // Same as above but as a vector including mu;ltiple DT cuts across multiple channel combinations
    gROOT->ProcessLine("Utils->SetChDCh()"); // This is a cut (!) for a minimum collected charge for an event to be considered
    gROOT->ProcessLine("Utils->SetChDChs()");  // Same as above but as a vector for mulitple channels
*/
    gROOT->ProcessLine("Utils->ConvertData()");
    gROOT->ProcessLine("Utils->SetWaveShape(false)");
    gROOT->ProcessLine("Utils->Analyse()");
}
