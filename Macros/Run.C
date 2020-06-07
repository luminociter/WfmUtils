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
    // Start stop event not supported (no effect)
    //gROOT->ProcessLine("Utils->SetInstrument(0)");
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

    gROOT->ProcessLine("Utils->ConvertData()");
    gROOT->ProcessLine("Utils->SetWaveShape(false)");
    gROOT->ProcessLine("Utils->Analyse()");
}