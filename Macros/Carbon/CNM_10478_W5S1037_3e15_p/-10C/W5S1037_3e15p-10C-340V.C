{
    gInterpreter->Reset(); // Reset Root

    gInterpreter->AddIncludePath("/afs/cern.ch/work/e/egkougko/private/LGADUtils/");

    gInterpreter->LoadMacro("/afs/cern.ch/work/e/egkougko/private/LGADUtils/Root/LGADBase.cxx+");
    gInterpreter->LoadMacro("/afs/cern.ch/work/e/egkougko/private/LGADUtils/Root/LGADFits.cxx+");
    gInterpreter->LoadMacro("/afs/cern.ch/work/e/egkougko/private/LGADUtils/Root/DataConverters.cxx+");
    gInterpreter->LoadMacro("/afs/cern.ch/work/e/egkougko/private/LGADUtils/Root/LGADSel.cxx+");
    gInterpreter->LoadMacro("/afs/cern.ch/work/e/egkougko/private/LGADUtils/Root/WaveForm.cxx+");
    gInterpreter->LoadMacro("/afs/cern.ch/work/e/egkougko/private/LGADUtils/Root/DUTChannel.cxx+");
    gInterpreter->LoadMacro("/afs/cern.ch/work/e/egkougko/private/LGADUtils/Root/LGADRun.cxx+");
    gInterpreter->LoadMacro("/afs/cern.ch/work/e/egkougko/private/LGADUtils/Root/LGADUtils.cxx+");

    gROOT->ProcessLine("LGADUtils* Utils = new LGADUtils()");
    gROOT->ProcessLine("Utils->SetVerbose(0)");
    gEnv->SetValue("TFile.AsyncPrefetching", 1);

    // Lab example
    gROOT->ProcessLine("Utils->SetInstrument(1)");
    gROOT->ProcessLine("Utils->SetStartStopEvnt(0, 0)");
    gROOT->ProcessLine("Utils->SetInDataNames(\"/eos/home-e/egkougko/Root/Carbon/CNM_10478_W5S1037_3e15_p/-10C\",\"340V.root\")");
    gROOT->ProcessLine("Utils->SetOutDataNames(\"/eos/home-e/egkougko/Root/Carbon/CNM_10478_W5S1037_3e15_p/-10C\")");

    // gROOT->ProcessLine("Utils->ConvertData()");
    gROOT->ProcessLine("Utils->SetWaveShape(true)");
    gROOT->ProcessLine("Utils->Analyse()");
}