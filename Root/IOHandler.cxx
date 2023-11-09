/*
*                 WFMRun1Stage.cxx
*
*
*      Author: Gkougkousis Evangelos - Leonidas
*              egkougko@cern.ch - 2023
*                    CERN-GENEVA
*/

#include "../WfmUtils/WfmRun.h"

bool WfmRun::InitNtuple(const std::vector<TTree*>& inputrees, const std::vector<TFile*>& inputfiles)
{
   // For each file figure out the number of channels and the instrument
   m_filechannel.clear();
   m_filelevels.clear();
   m_ChTrees.clear();
   m_ShapeTrees.clear();
   m_TimeResTrees.clear();
   m_FFTShDef.clear();
   m_fileplanes.clear();
   m_fileRefs.clear();
   m_fileDUTs.clear();
   m_filetrack.clear();
   m_RunDUTSubCh.clear();
   bool cp = true;
   unsigned int nfiles = inputrees.size();

   std::vector<AqInstrument> fileInstr;
   fileInstr.resize(nfiles, Unasigned);
   std::vector<bool> fftShDef;
   bool temp = false;
   for (unsigned int a = 0; a < nfiles; a++)
       {
        fftShDef.clear();
        fftShDef = { false, false, false, false, false, false, false };
        m_ShapeTrees.push_back(NULL);
        m_ChTrees.push_back(NULL);
        m_TimeResTrees.push_back(NULL);
        // Determine channels present in each file
        m_filechannel.push_back(std::vector<unsigned int>());
        for (unsigned int ich = 1; ich < 128; ich++) if (inputrees.at(a)->GetBranch(Form("t%02u", ich))) (m_filechannel.back()).push_back(ich);
        // Define instrument
        if (inputrees.at(a)->GetBranch(Form("ordrt%02u", (m_filechannel.back()).at(0)))) fileInstr.at(a) = Sampic;
        else if (inputrees.at(a)->GetBranch(Form("EvnTimeCh%02u", (m_filechannel.back()).at(0)))) fileInstr.at(a) = TektronixScope;
        else if (inputrees.at(a)->GetBranch(Form("Scope%02u", (m_filechannel.back()).at(0))))
                {
                 if (inputrees.at(a)->GetBranch("EvnTime")) fileInstr.at(a) = TestBeamBin2;
                 else fileInstr.at(a) = TestBeamBin1;
                }
        else if (inputrees.at(a)->GetBranch("EvnTime") && !inputrees.at(a)->GetBranch(Form("Scope%02u", (m_filechannel.back()).at(0)))) fileInstr.at(a) = LabTXT;
        else fileInstr.at(a) = LeCroyWRBin;
        if (m_RunBase->WfmBase::GetVerbosity() >= 1) std::cout << __FUNCTION__ << " INFO: Aquisition instrument for file " << inputfiles.at(a)->GetName()
                                                               << " identified as " << fileInstr.at(a) << "." << std::endl;
        // Define file Analysis, FFT and WaveShape levels, if time resolution branches exist
        if ((TTree*)(inputfiles.at(a)->Get("ChFitParam")))
           {
            temp = false;
            if ((TTree*)(inputfiles.at(a)->Get("TimeResParam"))) 
               {
                inputfiles.at(a)->GetObject("TimeResParam", m_TimeResTrees.back());
                temp = true;
               }
            inputfiles.at(a)->GetObject("ChFitParam", m_ChTrees.back());
            if (m_ChTrees.back()->GetBranch("SubChName"))
               {
                if (temp) m_filelevels.push_back(5);
                else m_filelevels.push_back(3);
                m_ChTrees.back()->SetBranchAddress("SubChName", &m_SubChName, &b_SubChName);
                if (m_ChTrees.back()->GetBranch("ChFitSignalFFT")) 
                   {
                    m_ChTrees.back()->SetBranchAddress("ChFitSignalFFT", &m_ChFitSignalFFT.first, &b_ChFitSignalFFT);
                    for (unsigned int a = 0; a < m_ChTrees.back()->GetEntriesFast(); a++)
                        {
                         b_SubChName->GetEntry(a);
                         b_ChFitSignalFFT->GetEntry(a);
                         if (m_SubChName.second != 0 && m_ChFitSignalFFT.first != -99) fftShDef.at(2) = true;
                         if (m_SubChName.second == 0 && m_ChFitSignalFFT.first != -99) fftShDef.at(1) = true;
                         if (fftShDef.at(1) && fftShDef.at(2)) break;
                        }                  
                   }
               }
            else {
                  if (temp) m_filelevels.push_back(4);
                  else m_filelevels.push_back(2);
                  if (m_ChTrees.back()->GetBranch("ChFitSignalFFT")) fftShDef.at(1) = true;
                 }
            if (fftShDef.at(1) || fftShDef.at(2)) fftShDef.at(0) = true;
            if ((TTree*)(inputfiles.at(a)->Get("ChShapeParam")))
               {
                temp = false;
                if (m_ShapeTrees.back()->GetBranch("ChShapeSignalFFT"))
                   {
                    m_ShapeTrees.back()->SetBranchAddress("ChShapeSignalFFT", &m_ChShapeSignalFFT, &b_ChShapeSignalFFT);
                    temp = true;
                   }
                inputfiles.at(a)->GetObject("ChShapeParam", m_ShapeTrees.back());
                if (m_ShapeTrees.back()->GetBranch("SubChShName"))
                   {   
                    fftShDef.at(6) = true;
                    m_ShapeTrees.back()->SetBranchAddress("SubChShName", &m_SubChShName, &b_SubChShName);
                    for (unsigned int a = 0; a < m_ShapeTrees.back()->GetEntriesFast(); a++)
                        {
                         b_SubChShName->GetEntry(a);
                         if (temp) b_ChShapeSignalFFT->GetEntry(a);
                         if (m_SubChShName == 0)
                            {
                             fftShDef.at(5) = true;
                             if (temp && m_ChShapeSignalFFT != -99) fftShDef.at(3) = true;
                            }
                         else if (temp && m_ChShapeSignalFFT != -99) fftShDef.at(4) = true;
                         if (fftShDef.at(3) && fftShDef.at(4) && fftShDef.at(5)) break;
                        }
                   }
                else {
                      fftShDef.at(5) = true;
                      if (temp) fftShDef.at(3) = true;
                     }
               }
           }
        else if (inputrees.at(a)->GetBranch(Form("Polarity%02u", (m_filechannel.back()).at(0)))) 
                {
                 m_filelevels.push_back(1);
                 if (inputrees.at(a)->GetBranch(Form("SignalFFT%02u", (m_filechannel.back()).at(0)))) fftShDef.at(0) = true;
                }
        else m_filelevels.push_back(0);
        if (m_RunBase->WfmBase::GetVerbosity() >= 1) std::cout << __FUNCTION__ << " INFO: Analysis level of file " << inputfiles.at(a)->GetName() << " :" << m_filelevels.back() << std::endl;
        // Presence of tracking information and trakcing package    
        if (inputrees.at(a)->GetBranch("trigger_id") && inputrees.at(a)->GetBranch("timestampTr")) m_filetrack.push_back(1); // EUTelescpe
        else if (inputrees.at(a)->GetBranch("HasAsso")) // Cory
                {
                 m_fileplanes.push_back(std::vector<unsigned int>());
                 m_fileRefs.push_back(std::vector<unsigned int>());
                 m_fileDUTs.push_back(std::vector<unsigned int>());
                 for (unsigned int pl = 0; pl < 8; pl++)
                     {
                      if (inputrees.at(a)->GetBranch(Form("InterceptX_TelPlane%u", pl))) (m_fileplanes.back()).push_back(pl);
                      if (inputrees.at(a)->GetBranch(Form("InterceptX_Ref%u", pl))) (m_fileRefs.back()).push_back(pl);
                      if (inputrees.at(a)->GetBranch(Form("InterceptX_DUT%u", pl))) (m_fileDUTs.back()).push_back(pl);
                     }
                 m_filetrack.push_back(2);
                }
        else if (inputrees.at(a)->GetBranch("blabalblabla")) m_filetrack.push_back(3); // TBA
        else m_filetrack.push_back(4); // No tracking
        m_FFTShDef.push_back(fftShDef);
        if (m_RunBase->WfmBase::GetInRootFile(a) == m_RunBase->WfmBase::GetOutRootFile()) cp = false;
       }

   // Fix the number of channels for reserving the correct amount of pointers
   m_RunBase->WfmBase::FixBaseQuant(m_filechannel, m_channels, m_nchan);
   m_RunBase->WfmBase::SetChannels(m_channels);
   if (m_RunBase->WfmBase::GetVerbosity() >= 1) std::cout << __FUNCTION__ << " INFO: Total number of channels found " << m_nchan << " in " << nfiles << " files." << std::endl;  // Debug

   // Verify that a single tracking frameworks and same qaquisition isntrument is used for all files
   unsigned int trcq = 0;
   if (m_filetrack.at(0) != 4) trcq = m_filetrack.at(0);
   for (unsigned int a = 1; a < nfiles; a++)
       {
        if (trcq == 0 && m_filetrack.at(a) != 4) trcq = m_filetrack.at(a);
        if (m_filetrack.at(a) != 4 && trcq != 0 && trcq != m_filetrack.at(a))
           { 
            std::cout << __FUNCTION__ << " ERROR: Trying to mix different trakcing frameworks on input files, aborting...." << std::endl;
            return false;
           }
        if (fileInstr.at(0) != fileInstr.at(a))
           { 
            std::cout << __FUNCTION__ << " ERROR: Trying to mix data from different measuring instruments, aborting...." << std::endl;
            return false;
           }
       }

   // Fix vector sizes for Cory files
   if (trcq == 2)
      {
       for (unsigned int a = 0; nfiles; a++)
           {
            if (m_filetrack.at(a) == 2) continue;
            else {
                  m_fileplanes.insert(m_fileplanes.begin() + a, std::vector<unsigned int>());
                  m_fileRefs.insert(m_fileRefs.begin() + a, std::vector<unsigned int>());
                  m_fileDUTs.insert(m_fileDUTs.begin() + a, std::vector<unsigned int>());
                 }
           }
      }
   if (trcq != 0)
      { 
       m_RunBase->WfmBase::SetHasTrck(true);
       if (trcq == 1) m_RunBase->WfmBase::SetTrackPackage("EUTelescope ");
       else if (trcq == 2) m_RunBase->WfmBase::SetTrackPackage("Cory");
       else if (trcq == 3) m_RunBase->WfmBase::SetTrackPackage("TBA");
      }
   else {
         m_RunBase->WfmBase::SetHasTrck(false);
         m_RunBase->WfmBase::SetTrackPackage("");
        }
   m_RunBase->WfmBase::SetInstrument(fileInstr.at(0));

   // Fix the file level for initializing the right pointers
   std::pair<unsigned int, unsigned int> lvl = std::make_pair(*std::min_element(m_filelevels.begin(), m_filelevels.end()), 
                                                              *std::max_element(m_filelevels.begin(), m_filelevels.end()));
   m_RunBase->WfmBase::SetFileLVL(lvl);

   // Establish if waveshape and FFT pointers are to be reserved
   std::vector<std::pair<bool, bool>> fft;
   fft.resize(5, std::make_pair(true, false));
   std::pair<bool, bool> shape_min = std::make_pair(true, true);
   std::pair<bool, bool> shape_max = std::make_pair(false, false);
   for (unsigned int p = 0; p < m_filelevels.size(); p++)
       { 
        // This has to be modified now
        for (unsigned int k = 0; k < 5; k++) 
            {
             if ((m_FFTShDef.at(p)).at(k) && !(fft.at(k)).second) (fft.at(k)).second = true;
             else if (!(m_FFTShDef.at(p)).at(k) && (fft.at(k)).first) (fft.at(k)).first = false;
            }
        if (!(m_FFTShDef.at(p)).at(5)) { if (shape_min.first) shape_min.first = false; }
        else { if (!shape_max.first) shape_max.first = true; }
        if (!(m_FFTShDef.at(p)).at(6)) { if (shape_min.second) shape_min.second = false; }
        else { if (!shape_max.second) shape_max.second = true; }
        if ((fft.at(0)).second && !(fft.at(0)).first && (fft.at(1)).second && !(fft.at(1)).first &&
            (fft.at(2)).second && !(fft.at(2)).first && (fft.at(3)).second && !(fft.at(3)).first &&
            (fft.at(4)).second && !(fft.at(4)).first && !shape_min.first && !shape_min.second &&
            shape_max.first && shape_max.second) break;
       }
   m_RunBase->WfmBase::SetHasFFT(fft);
   // Set WaveShape variable
   std::pair<unsigned int, unsigned int> WaveShape = std::make_pair(99, 99);
   if (shape_min.first && shape_min.second) WaveShape.first = 3;
   else if (!shape_min.first && shape_min.second) WaveShape.first = 2;
   else if (shape_min.first && !shape_min.second) WaveShape.first = 1;
   else WaveShape.first = 0;
   if (shape_max.first && shape_max.second) WaveShape.second = 3;
   else if (!shape_max.first && shape_max.second) WaveShape.second = 2;
   else if (shape_max.first && !shape_max.second) WaveShape.second = 1;
   else WaveShape.second = 0;
   m_RunBase->WfmBase::SetHasWaveShape(WaveShape);

   // Initialize zero level Pointers and vectors
   b_t.resize(m_nchan, nullptr);
   b_w.resize(m_nchan, nullptr);
   b_ordrt.resize(m_nchan, nullptr);
   b_physt.resize(m_nchan, nullptr);
   b_vScale.resize(m_nchan, nullptr);
   b_nPoints.resize(m_nchan, nullptr);
   b_SnRate.resize(m_nchan, nullptr);
   b_Scope.resize(m_nchan, nullptr);
   b_EvnTimeCh.resize(m_nchan, nullptr);
   WfmBase::SetVectorSize(m_nchan, 0);
   if (m_RunBase->WfmBase::GetInstrument() != TestBeamBin1 && m_RunBase->WfmBase::GetInstrument() != LeCroyWRBin)
      { 
       m_trigDt_In.resize(nfiles, nullptr); 
       m_trigFr_In.resize(nfiles, nullptr);
      }
   // THis is the Fit study algorithm
   unsigned int ftlvl1 = 0;
   for (unsigned int p = 0; p < m_filelevels.size(); p++) if (m_filelevels.at(p) > 0) ftlvl1++;
   for (unsigned int h = 0; h < 5; h++) m_FitVarIn[h].reserve(ftlvl1);
      
   // Initialize Cory tracking pointers size
   if (m_RunBase->WfmBase::GetHasTrck() && !m_RunBase->WfmBase::GetTrackComb())
      {
       if (m_RunBase->WfmBase::GetTrackPackage() == "Cory")
          {
           m_ntp_planes.clear();
           m_ntp_Ref.clear();
           m_ntp_DUT.clear();
           m_nRecDUTs = 0;
           m_nRecRefs = 0;
           m_nRecPlanes = 0;
           m_RunBase->WfmBase::FixBaseQuant(m_fileDUTs, m_ntp_DUT, m_nRecDUTs);
           m_RunBase->WfmBase::FixBaseQuant(m_fileRefs, m_ntp_Ref, m_nRecRefs);
           m_RunBase->WfmBase::FixBaseQuant(m_fileplanes, m_ntp_planes, m_nRecPlanes);

           // Reference Branches
           b_InterceptX_Ref.resize(m_nRecRefs, nullptr);
           b_InterceptY_Ref.resize(m_nRecRefs, nullptr);
           b_InterceptResX_Ref.resize(m_nRecRefs, nullptr);
           b_InterceptResY_Ref.resize(m_nRecRefs, nullptr);
           // DUT Branches
           b_InterceptX_DUT.resize(m_nRecDUTs, nullptr);
           b_InterceptY_DUT.resize(m_nRecDUTs, nullptr);
           b_InterceptResX_DUT.resize(m_nRecDUTs, nullptr);
           b_InterceptResY_DUT.resize(m_nRecDUTs, nullptr);
           // Telescope Plane Branches
           b_InterceptX_TelPlane.resize(m_nRecPlanes, nullptr);
           b_InterceptY_TelPlane.resize(m_nRecPlanes, nullptr);
           b_InterceptResX_TelPlane.resize(m_nRecPlanes, nullptr);
           b_InterceptResY_TelPlane.resize(m_nRecPlanes, nullptr);
           b_ClusterSizeX_TelPlane.resize(m_nRecPlanes, nullptr);
           b_ClusterSizeY_TelPlane.resize(m_nRecPlanes, nullptr);
           b_ClusterNumPixels_TelPlane.resize(m_nRecPlanes, nullptr);
           m_RunBase->WfmBase::InitCoryPointers(m_nRecRefs, m_nRecDUTs, m_nRecPlanes);
          }
       if (m_RunBase->WfmBase::GetTrackPackage() == "EUTelescope") m_RunBase->WfmBase::InitEUTelPointers();
       if (m_RunBase->WfmBase::GetTrackPackage() == "TBA") m_RunBase->WfmBase::InitTBAPointers();
      }

   // Intialize 1st level pointers and vectors
   if (lvl.second >=0 || m_RunBase->WfmBase::GetAnaStage() > 0)
      {
       if (lvl.second >= 0)
          {
           b_EvPol.resize(m_nchan, nullptr);
           b_EvCharge.resize(m_nchan, nullptr);
           b_EvIsSignal.resize(m_nchan, nullptr);
           b_EvInWindow.resize(m_nchan, nullptr);
           b_EvVoltStr.resize(m_nchan, nullptr);
           b_EvMaxindx.resize(m_nchan, nullptr);
           b_EvMaxTime.resize(m_nchan, nullptr);
           b_EvMaxVolt.resize(m_nchan, nullptr);
           b_EvMinIndx.resize(m_nchan, nullptr);
           b_EvMinTime.resize(m_nchan, nullptr);
           b_EvMinVolt.resize(m_nchan, nullptr);
           b_EvStrIndx.resize(m_nchan, nullptr);
           b_EvEndIndx.resize(m_nchan, nullptr);
           b_EvNoise.resize(m_nchan, nullptr);
           b_EvPedest.resize(m_nchan, nullptr);
           b_EvNoiseErr.resize(m_nchan, nullptr);
           b_EvPedestErr.resize(m_nchan, nullptr);
           b_EvRiseT.resize(m_nchan, nullptr);
           b_EvJitNdVdT.resize(m_nchan, nullptr);
           b_EvJitRiseSNR.resize(m_nchan, nullptr);
           b_EvTriggTime.resize(m_nchan, nullptr);
           b_EvDVDTMax.resize(m_nchan, nullptr);
           b_EvDVDTCFD.resize(m_nchan, nullptr);
           b_EvCFDToT.resize(m_nchan, nullptr);
           b_EvTriggToT.resize(m_nchan, nullptr);
           b_EvComplete.resize(m_nchan, nullptr);
           b_EvVAdjCFD.resize(m_nchan, nullptr);
           b_EvCFDTime.resize(m_nchan, nullptr);
           if ((fft.at(1)).second)
              {
               b_EvSignalFFT.resize(m_nchan, nullptr);
               b_EvNoiseFFT.resize(m_nchan, nullptr);
              }
          }
       m_EvPol.resize(m_nchan, -99);
       m_EvIsSignal.resize(m_nchan, -99);
       m_EvInWindow.resize(m_nchan, -99);
       m_EvVoltStr.resize(m_nchan, -99);
       m_EvMaxindx.resize(m_nchan, -99);
       m_EvStrIndx.resize(m_nchan, -99);
       m_EvEndIndx.resize(m_nchan, -99);
       m_EvMinIndx.resize(m_nchan, -99);
       m_EvComplete.resize(m_nchan, -99);
       m_EvCharge.resize(m_nchan, -99.0);
       m_EvMaxTime.resize(m_nchan, -99.0);
       m_EvMaxVolt.resize(m_nchan, -99.0);
       m_EvMinTime.resize(m_nchan, -99.0);
       m_EvMinVolt.resize(m_nchan, -99.0);
       m_EvNoise.resize(m_nchan, -99.0);
       m_EvPedest.resize(m_nchan, -99.0);
       m_EvNoiseErr.resize(m_nchan, -99.0);
       m_EvPedestErr.resize(m_nchan, -99.0);
       m_EvRiseT.resize(m_nchan, -99.0);
       m_EvJitRiseSNR.resize(m_nchan, -99.0);
       m_EvTriggTime.resize(m_nchan, -99.0);
       m_EvDVDTMax.resize(m_nchan, -99.0);
       m_EvTriggToT.resize(m_nchan, -99.0);
       // Initalization of vector pointrers
       m_EvVAdjCFD.clear();
       m_EvCFDTime.clear();
       m_EvJitNdVdT.clear();
       m_EvDVDTCFD.clear();
       m_EvCFDToT.clear();
       m_EvVAdjCFD.resize(m_nchan, std::vector<double>());
       m_EvCFDTime.resize(m_nchan, std::vector<double>());
       m_EvJitNdVdT.resize(m_nchan, std::vector<double>());
       m_EvDVDTCFD.resize(m_nchan, std::vector<double>());
       m_EvCFDToT.resize(m_nchan, std::vector<double>());
       if ((fft.at(1)).second || m_RunBase->WfmBase::GetDoFFT())
          {
           m_EvSignalFFT.resize(m_nchan, -99.0);
           m_EvNoiseFFT.resize(m_nchan, -99.0);
          }
      }

    if (lvl.second > 1 || m_RunBase->WfmBase::GetAnaStage() > 1)
       {
        m_CFDTimeFtChi2 = 0;
        m_DVDTCFDFtChi2 = 0;
        m_CFDToTFtChi2 = 0;
        m_ChFitCFDTime = 0;
        m_ChFitDVDTCFD = 0;
        m_ChFitCFDToT = 0;
        m_ChFitJitNdVdT = 0;
        if (shape_max.first || shape_max.second)
           {
            m_ChShapeCFDTime = 0;
            m_ChShapeDVDTCFD = 0;
            m_ChShapeCFDToT = 0;
            m_ChShapeJitNdVdT = 0;
            m_ChMeanVolt = 0;
           }
       }

    // Intialize 4th level pointers and vectors
    if (lvl.second > 3 || m_RunBase->WfmBase::GetAnaStage() > 3)
       {
        m_TmDiffChi2 = 0;
        m_TmDiffFitMean = 0;
        m_TmDiffFitErr = 0;
        m_CFD_ChA = 0;
       }

   return true;
}
// --------------------------------------------------------------------------------------------------------------
// Function to be called every time a new NTuple is opened to correctly set the branch addresses for the varibales present in the dataset
bool WfmRun::SetBranchAccess(unsigned int file)
{
   // Variables to assign for non-existent branches need to be populated
   double no_ch_double = -99;
   int no_ch_int = -99;
   unsigned int no_ch_uint = 99;
   float no_ch_flot = -99.;
   Long64_t no_ch_long = -99;
   std::vector<double> no_ch_doub_vec;
   no_ch_doub_vec.resize(1, -99);
   std::vector<int> no_ch_int_vec;
   no_ch_int_vec.resize(1,-99);

   // Add the pointers to the histograms from the input files either whwn processing the entire file or when merging to one of the files
   if ((m_RunBase->WfmBase::GetInstrument() != TestBeamBin1 && m_RunBase->WfmBase::GetInstrument() != LeCroyWRBin &&
       ((m_RunBase->WfmBase::GetStartStopEvnt().first == 0 && m_RunBase->WfmBase::GetStartStopEvnt().second == 0) || m_cp != -99)) || m_RunBase->WfmBase::GetInstrument() == Sampic)
      {
       m_trigDt_In.at(file) = dynamic_cast<TH1F*>(m_RunBase->WfmBase::GetInRootFile(file)->Get("trigDt"));
       m_trigFr_In.at(file) = dynamic_cast<TH1F*>(m_RunBase->WfmBase::GetInRootFile(file)->Get("trigFr"));
      }

   // If branches are available assign them independantly if used or not. This will be decided in the analysis stage
   (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress("EvnNo", &EvnNo, &b_EvnNo);
   for (unsigned int ich = 0; ich < (m_filechannel.at(file)).size(); ich++)
       {
        unsigned int chindx = 999;
        for (unsigned int s = 0; s < (m_RunBase->WfmBase::GetChannels()).size(); s++) { if ((m_filechannel.at(file)).at(ich) == (m_RunBase->WfmBase::GetChannels()).at(s)) { chindx = s; break; } }
        if (chindx == 999)
           {
            std::cout << __FUNCTION__ << " ERROR: Could not match channel " << (m_filechannel.at(file)).at(ich) << " of this file, skipping...." << std::endl;
            return false;
           }
        if (m_RunBase->WfmBase::GetVerbosity() >= 1) std::cout << __FUNCTION__ << " INFO: Channel: " << ich+1  << " indexed as: " << (m_filechannel.at(file)).at(ich) << " with index @ main channel vector: " 
                                                               << chindx << " Vector sizes: " << m_t.size() << ", " << m_w.size() << ", " << b_t.size() << ", " << b_w.size() << std::endl;
        (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("t%02u", (m_filechannel.at(file)).at(ich)), &m_t.at(chindx), &b_t.at(chindx));
        (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("w%02u", (m_filechannel.at(file)).at(ich)), &m_w.at(chindx), &b_w.at(chindx));
        if (m_RunBase->WfmBase::GetInstrument() == Sampic)
           {
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("ordrt%02u", (m_filechannel.at(file)).at(ich)), &m_ordrt.at(chindx), &b_ordrt.at(chindx));
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("physt%02u", (m_filechannel.at(file)).at(ich)), &m_physt.at(chindx), &b_physt.at(chindx));
           }
        if (m_RunBase->WfmBase::GetInstrument() == TestBeamBin1 || m_RunBase->WfmBase::GetInstrument() == TestBeamBin2 || m_RunBase->WfmBase::GetInstrument() == LabTXT || 
            m_RunBase->WfmBase::GetInstrument() == TektronixScope || m_RunBase->WfmBase::GetInstrument() == LeCroyWRBin)
           {
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("nPoints%02u", (m_filechannel.at(file)).at(ich)), &m_npoints.at(chindx), &b_nPoints.at(chindx));
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("SnRate%02u", (m_filechannel.at(file)).at(ich)), &m_srate.at(chindx), &b_SnRate.at(chindx));
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("vScale%02u", (m_filechannel.at(file)).at(ich)), &m_scale.at(chindx), &b_vScale.at(chindx));
            if (m_RunBase->WfmBase::GetInstrument() == TektronixScope) 
                (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("EvnTimeCh%02u", (m_filechannel.at(file)).at(ich)), &m_EvnTimeCh.at(chindx), &b_EvnTimeCh.at(chindx));
            if (m_RunBase->WfmBase::GetInstrument() == TestBeamBin1 || m_RunBase->WfmBase::GetInstrument() == TestBeamBin2 || m_RunBase->WfmBase::GetInstrument() == LeCroyWRBin)
                (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("Scope%02u", (m_filechannel.at(file)).at(ich)), &m_scope.at(chindx), &b_Scope.at(chindx));
           }
        // First stage analysis branches if they exist
        if (m_filelevels.at(file) > 0)
           {
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("Polarity%02u", (m_filechannel.at(file)).at(ich)), &m_EvPol.at(chindx), &b_EvPol.at(chindx));
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("Charge%02u", (m_filechannel.at(file)).at(ich)), &m_EvCharge.at(chindx), &b_EvCharge.at(chindx));
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("IsSignal%02u", (m_filechannel.at(file)).at(ich)), &m_EvIsSignal.at(chindx), &b_EvIsSignal.at(chindx));
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("InWindow%02u", (m_filechannel.at(file)).at(ich)), &m_EvInWindow.at(chindx), &b_EvInWindow.at(chindx));
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("VoltStr%02u", (m_filechannel.at(file)).at(ich)), &m_EvVoltStr.at(chindx), &b_EvVoltStr.at(chindx));
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("MaxIndx%02u", (m_filechannel.at(file)).at(ich)), &m_EvMaxindx.at(chindx), &b_EvMaxindx.at(chindx));
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("MaxTime%02u", (m_filechannel.at(file)).at(ich)), &m_EvMaxTime.at(chindx), &b_EvMaxTime.at(chindx));
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("MaxVolt%02u", (m_filechannel.at(file)).at(ich)), &m_EvMaxVolt.at(chindx), &b_EvMaxVolt.at(chindx));
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("MinIndx%02u", (m_filechannel.at(file)).at(ich)), &m_EvMinIndx.at(chindx), &b_EvMinIndx.at(chindx));
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("MinTime%02u", (m_filechannel.at(file)).at(ich)), &m_EvMinTime.at(chindx), &b_EvMinTime.at(chindx));
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("MinVolt%02u", (m_filechannel.at(file)).at(ich)), &m_EvMinVolt.at(chindx), &b_EvMinVolt.at(chindx));
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("Strindx%02u", (m_filechannel.at(file)).at(ich)), &m_EvStrIndx.at(chindx), &b_EvStrIndx.at(chindx));
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("EndIndx%02u", (m_filechannel.at(file)).at(ich)), &m_EvEndIndx.at(chindx), &b_EvEndIndx.at(chindx));
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("Noise%02u", (m_filechannel.at(file)).at(ich)), &m_EvNoise.at(chindx), &b_EvNoise.at(chindx));
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("Pedest%02u", (m_filechannel.at(file)).at(ich)), &m_EvPedest.at(chindx), &b_EvPedest.at(chindx));
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("NoiseErr%02u", (m_filechannel.at(file)).at(ich)), &m_EvNoiseErr.at(chindx), &b_EvNoiseErr.at(chindx));
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("PedestErr%02u", (m_filechannel.at(file)).at(ich)), &m_EvPedestErr.at(chindx), &b_EvPedestErr.at(chindx));
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("RiseT%02u", (m_filechannel.at(file)).at(ich)), &m_EvRiseT.at(chindx), &b_EvRiseT.at(chindx));
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("JitRiseSNR%02u", (m_filechannel.at(file)).at(ich)), &m_EvJitRiseSNR.at(chindx), &b_EvJitRiseSNR.at(chindx));
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("TriggTime%02u", (m_filechannel.at(file)).at(ich)), &m_EvTriggTime.at(chindx), &b_EvTriggTime.at(chindx));
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("DVDTMax%02u", (m_filechannel.at(file)).at(ich)), &m_EvDVDTMax.at(chindx), &b_EvDVDTMax.at(chindx));
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("TriggToT%02u", (m_filechannel.at(file)).at(ich)), &m_EvTriggToT.at(chindx), &b_EvTriggToT.at(chindx));
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("Complete%02u", (m_filechannel.at(file)).at(ich)), &m_EvComplete.at(chindx), &b_EvComplete.at(chindx));
            if ((m_FFTShDef.at(file)).at(0))
               {
                (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("SignalFFT%02u", (m_filechannel.at(file)).at(ich)), &m_EvSignalFFT.at(chindx), &b_EvSignalFFT.at(chindx));
                (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("NoiseFFT%02u", (m_filechannel.at(file)).at(ich)), &m_EvNoiseFFT.at(chindx), &b_EvNoiseFFT.at(chindx));
               }
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("VAdjCFD%02u", (m_filechannel.at(file)).at(ich)), &m_EvVAdjCFD.at(chindx), &b_EvVAdjCFD.at(chindx));
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("CFDTime%02u", (m_filechannel.at(file)).at(ich)), &m_EvCFDTime.at(chindx), &b_EvCFDTime.at(chindx));
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("DVDTCFD%02u", (m_filechannel.at(file)).at(ich)), &m_EvDVDTCFD.at(chindx), &b_EvDVDTCFD.at(chindx));
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("CFDToT%02u", (m_filechannel.at(file)).at(ich)), &m_EvCFDToT.at(chindx), &b_EvCFDToT.at(chindx));
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("JitNdVdT%02u", (m_filechannel.at(file)).at(ich)), &m_EvJitNdVdT.at(chindx), &b_EvJitNdVdT.at(chindx));
           }
       }
   if (m_filelevels.at(file) > 0)
      {
       m_FitVarIn[0].push_back(dynamic_cast<TH1F*>(m_RunBase->WfmBase::GetInRootFile(file)->Get("VarFit_Mag")));
       m_FitVarIn[1].push_back(dynamic_cast<TH1F*>(m_RunBase->WfmBase::GetInRootFile(file)->Get("VarFit_MagErr")));
       m_FitVarIn[2].push_back(dynamic_cast<TH1F*>(m_RunBase->WfmBase::GetInRootFile(file)->Get("VarFit_Sigma")));
       m_FitVarIn[3].push_back(dynamic_cast<TH1F*>(m_RunBase->WfmBase::GetInRootFile(file)->Get("VarFit_SigmaErr")));
       m_FitVarIn[4].push_back(dynamic_cast<TH1F*>(m_RunBase->WfmBase::GetInRootFile(file)->Get("VarFit_Chi2")));
       (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress("FitDiffChi2", &m_FitDiffChi2, &b_FitDiffChi2);
       (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress("FitDiffMag", &m_FitDiffMag, &b_FitDiffMag);
       (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress("FitDiffMagErr", &m_FitDiffMagErr, &b_FitDiffMagErr);
       (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress("FitDiffSigma", &m_FitDiffSigma, &b_FitDiffSigma);
       (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress("FitDiffSigmaErr", &m_FitDiffSigmaErr, &b_FitDiffSigmaErr);
       (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress("FitDiffTime", &m_FitDiffTime, &b_FitDiffTime;
      }
   // Non-existent channels should have reset branches with zero pointers vectors set to -99
   std::vector<unsigned int> not_channels;
   bool found = false;
   for (unsigned int p = 0; p < (m_RunBase->WfmBase::GetChannels()).size(); p++)
       { 
        found = false;
        for (unsigned int a = 0; a < (m_filechannel.at(file)).size(); a++)
            {
             if ((m_RunBase->WfmBase::GetChannels()).at(p) == (m_filechannel.at(file)).at(a))
                {
                 found = true;
                 break;
                }
            }
        if (!found) not_channels.push_back(p);
       }
   if (not_channels.size() > 0)
      { 
       for (unsigned int s = 0; s < not_channels.size(); s++)
           {
            b_t.at(not_channels.at(s)) = nullptr;
            b_w.at(not_channels.at(s)) = nullptr;
            *(m_t.at(not_channels.at(s))) = no_ch_doub_vec;
            *(m_w.at(not_channels.at(s))) = no_ch_doub_vec;
            if (m_RunBase->WfmBase::GetInstrument() == Sampic)
               {
                b_ordrt.at(not_channels.at(s)) = nullptr;
                b_physt.at(not_channels.at(s)) = nullptr;
                m_ordrt.at(not_channels.at(s)) = no_ch_double;
                m_physt.at(not_channels.at(s)) = no_ch_double;
               }
            if (m_RunBase->WfmBase::GetInstrument() == TestBeamBin1 || m_RunBase->WfmBase::GetInstrument() == TestBeamBin2 || m_RunBase->WfmBase::GetInstrument() == LabTXT || 
                m_RunBase->WfmBase::GetInstrument() == TektronixScope || m_RunBase->WfmBase::GetInstrument() == LeCroyWRBin)
               {
                b_nPoints.at(not_channels.at(s)) = nullptr;
                b_SnRate.at(not_channels.at(s)) = nullptr;
                b_vScale.at(not_channels.at(s)) = nullptr;
                m_npoints.at(not_channels.at(s)) = no_ch_uint;
                m_srate.at(not_channels.at(s)) = no_ch_long;
                m_scale.at(not_channels.at(s)) = no_ch_flot;
                if (m_RunBase->WfmBase::GetInstrument() == TektronixScope) 
                   {
                    b_EvnTimeCh.at(not_channels.at(s)) = nullptr;
                    m_EvnTimeCh.at(not_channels.at(s)) = no_ch_double;
                   }
                if (m_RunBase->WfmBase::GetInstrument() == TestBeamBin1 || m_RunBase->WfmBase::GetInstrument() == TestBeamBin2 || m_RunBase->WfmBase::GetInstrument() == LeCroyWRBin)
                   {
                    b_Scope.at(not_channels.at(s)) = nullptr;
                    m_scope.at(not_channels.at(s)) = no_ch_uint;
                   }
               }
            // First level analysis branches
            if (m_filelevels.at(file) > 0)
               {
                b_EvPol.at(not_channels.at(s)) = nullptr;
                b_EvCharge.at(not_channels.at(s)) = nullptr;
                b_EvIsSignal.at(not_channels.at(s)) = nullptr;
                b_EvInWindow.at(not_channels.at(s)) = nullptr;
                b_EvVoltStr.at(not_channels.at(s)) = nullptr;
                b_EvMaxindx.at(not_channels.at(s)) = nullptr;
                b_EvMaxTime.at(not_channels.at(s)) = nullptr;
                b_EvMaxVolt.at(not_channels.at(s)) = nullptr;
                b_EvMinIndx.at(not_channels.at(s)) = nullptr;
                b_EvMinTime.at(not_channels.at(s)) = nullptr;
                b_EvMinVolt.at(not_channels.at(s)) = nullptr;
                b_EvStrIndx.at(not_channels.at(s)) = nullptr;
                b_EvEndIndx.at(not_channels.at(s)) = nullptr;
                b_EvNoise.at(not_channels.at(s)) = nullptr;
                b_EvPedest.at(not_channels.at(s)) = nullptr;
                b_EvNoiseErr.at(not_channels.at(s)) = nullptr;
                b_EvPedestErr.at(not_channels.at(s)) = nullptr;
                b_EvRiseT.at(not_channels.at(s)) = nullptr;
                b_EvJitRiseSNR.at(not_channels.at(s)) = nullptr;
                b_EvTriggTime.at(not_channels.at(s)) = nullptr;
                b_EvDVDTMax.at(not_channels.at(s)) = nullptr;
                b_EvTriggToT.at(not_channels.at(s)) = nullptr;
                b_EvComplete.at(not_channels.at(s)) = nullptr;
                m_EvPol.at(not_channels.at(s)) = no_ch_int;
                m_EvCharge.at(not_channels.at(s)) = no_ch_double;
                m_EvIsSignal.at(not_channels.at(s)) = no_ch_int;
                m_EvInWindow.at(not_channels.at(s)) = no_ch_int;
                m_EvVoltStr.at(not_channels.at(s)) = no_ch_int;
                m_EvMaxindx.at(not_channels.at(s)) = no_ch_int;
                m_EvMaxTime.at(not_channels.at(s)) = no_ch_double;
                m_EvMaxVolt.at(not_channels.at(s)) = no_ch_double;
                m_EvMinIndx.at(not_channels.at(s)) = no_ch_int;
                m_EvMinTime.at(not_channels.at(s)) = no_ch_double;
                m_EvMinVolt.at(not_channels.at(s)) = no_ch_double;
                m_EvStrIndx.at(not_channels.at(s)) = no_ch_int;
                m_EvEndIndx.at(not_channels.at(s)) = no_ch_int;
                m_EvNoise.at(not_channels.at(s)) = no_ch_double;
                m_EvPedest.at(not_channels.at(s)) = no_ch_double;
                m_EvNoiseErr.at(not_channels.at(s)) = no_ch_double;
                m_EvPedestErr.at(not_channels.at(s)) = no_ch_double;
                m_EvRiseT.at(not_channels.at(s)) = no_ch_double;
                m_EvJitRiseSNR.at(not_channels.at(s)) = no_ch_double;
                m_EvTriggTime.at(not_channels.at(s)) = no_ch_double;
                m_EvDVDTMax.at(not_channels.at(s)) = no_ch_double;
                m_EvTriggToT.at(not_channels.at(s)) = no_ch_double;
                m_EvComplete.at(not_channels.at(s)) = no_ch_int;
                if ((m_FFTShDef.at(file)).at(0))
                   {
                    b_EvSignalFFT.at(not_channels.at(s)) = nullptr;
                    b_EvNoiseFFT.at(not_channels.at(s)) = nullptr;
                    m_EvSignalFFT.at(not_channels.at(s)) = no_ch_double;
                    m_EvNoiseFFT.at(not_channels.at(s)) = no_ch_double;
                   }
                b_EvVAdjCFD.at(not_channels.at(s)) = nullptr;
                b_EvCFDTime.at(not_channels.at(s)) = nullptr;
                b_EvDVDTCFD.at(not_channels.at(s)) = nullptr;
                b_EvCFDToT.at(not_channels.at(s)) = nullptr;
                b_EvJitNdVdT.at(not_channels.at(s)) = nullptr;
                m_EvVAdjCFD.at(not_channels.at(s)) = no_ch_doub_vec;
                m_EvCFDTime.at(not_channels.at(s)) = no_ch_doub_vec;
                m_EvDVDTCFD.at(not_channels.at(s)) = no_ch_doub_vec;
                m_EvCFDToT.at(not_channels.at(s)) = no_ch_doub_vec;
                m_EvJitNdVdT.at(not_channels.at(s)) = no_ch_doub_vec;
               }
           }
      }

   if (m_RunBase->WfmBase::GetInstrument() == Sampic)
      {
       (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress("nPoints", &m_npoints.at(0), &b_nPoints.at(0));
       (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress("SnRate", &m_srate.at(0), &b_SnRate.at(0));
      }
   if (m_RunBase->WfmBase::GetInstrument() == LabTXT || m_RunBase->WfmBase::GetInstrument() == TestBeamBin2) 
        (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress("EvnTime", &m_EvnTime, &b_EvnTime);

   // EU Telescope Branches
   if (m_RunBase->WfmBase::GetTrackPackage() == "EUTelescope" && m_filetrack.at(file) == 1)
      {
       (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress("xTr", &m_xTr, &b_xTr);
       (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress("yTr", &m_yTr, &b_yTr);
       (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress("omegaTr", &m_omegaTr, &b_omegaTr);
       (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress("phiTr", &m_phiTr, &b_phiTr);
       (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress("kinkxTr", &m_kinkxTr, &b_kinkxTr);
       (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress("kinkyTr", &m_kinkyTr, &b_kinkyTr);
       (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress("chi2Tr", &m_chi2Tr, &b_chi2Tr);
       (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress("plane_id", &m_plane_id, &b_plane_id);
       (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress("track_id", &m_track_id, &b_track_id);
       (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress("ndofTr", &m_ndofTr, &b_ndofTr);
       (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress("trigger_id", &m_trigger_id, &b_trigger_id);
       (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress("timestampTr", &m_timestampTr, &b_timestampTr);
       (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress("xHit", &m_xHit, &b_xHit);
       (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress("yHit", &m_yHit, &b_yHit);
       (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress("zHit", &m_zHit, &b_zHit);
       (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress("ID", &m_ID, &b_ID);
      }
   else if (m_RunBase->WfmBase::GetTrackPackage() == "Cory" && m_filetrack.at(file) == 2)
           {
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress("ClusterSizeX_ROI", &m_ClusterSizeX_ROI, &b_ClusterSizeX_ROI);
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress("ClusterSizeY_ROI", &m_ClusterSizeY_ROI, &b_ClusterSizeY_ROI);
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress("PixelX_ROI", &m_PixelX_ROI, &b_PixelX_ROI);
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress("PixelY_ROI", &m_PixelY_ROI, &b_PixelY_ROI);
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress("ClusterNumPixels_ROI", &m_ClusterNumPixels_ROI, &b_ClusterNumPixels_ROI);
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress("InterceptX_ROI", &m_InterceptX_ROI, &b_InterceptX_ROI);
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress("InterceptY_ROI", &m_InterceptY_ROI, &b_InterceptY_ROI);
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress("InterceptResX_ROI", &m_InterceptResX_ROI, &b_InterceptResX_ROI);
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress("InterceptResY_ROI", &m_InterceptResY_ROI, &b_InterceptResY_ROI);
            (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress("HasAsso", &m_HasAsso, &b_HasAsso);
            for (unsigned int ref = 0; ref < (m_fileRefs.at(file)).size(); ref++)
                {
                 unsigned int refindx = 999;
                 for (unsigned int s = 0; s < m_ntp_Ref.size(); s++) { if (m_ntp_Ref.at(s) == (m_fileRefs.at(file)).at(ref)) { refindx = s; break; } }
                 if (refindx == 999)
                    {
                     std::cout << __FUNCTION__ << " ERROR: Could not match refernce plane " << (m_fileRefs.at(file)).at(ref) << " for this run, skipping...." << std::endl;
                     return false;
                    }
                 (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("InterceptX_Ref%u", (m_fileRefs.at(file)).at(ref)), m_InterceptX_Ref.at(refindx), &b_InterceptX_Ref.at(refindx));
                 (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("InterceptY_Ref%u", (m_fileRefs.at(file)).at(ref)), m_InterceptY_Ref.at(refindx), &b_InterceptY_Ref.at(refindx));
                 (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("InterceptResX_Ref%u", (m_fileRefs.at(file)).at(ref)), m_InterceptResX_Ref.at(refindx), &b_InterceptResX_Ref.at(refindx));
                 (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("InterceptResY_Ref%u", (m_fileRefs.at(file)).at(ref)), m_InterceptResY_Ref.at(refindx), &b_InterceptResY_Ref.at(refindx));
                }
            // Non-existent References need to be populated correctly
            not_channels.clear();
            for (unsigned int p = 0; p < m_ntp_Ref.size(); p++)
                {
                 found = false;
                 for (unsigned int a = 0; a < (m_fileRefs.at(file)).size(); a++)
                     {
                      if (m_ntp_Ref.at(p) == (m_fileRefs.at(file)).at(a))
                         {
                          found = true;
                          break;
                         }
                     }
                 if (!found) not_channels.push_back(p);
                }
            if (not_channels.size() > 0)
               { 
                for (unsigned int s = 0; s < not_channels.size(); s++)
                    {
                     b_InterceptX_Ref.at(not_channels.at(s)) = nullptr;
                     b_InterceptY_Ref.at(not_channels.at(s)) = nullptr;
                     b_InterceptResX_Ref.at(not_channels.at(s)) = nullptr;
                     b_InterceptResY_Ref.at(not_channels.at(s)) = nullptr;
                     m_InterceptX_Ref.at(not_channels.at(s)) = &no_ch_doub_vec;
                     m_InterceptY_Ref.at(not_channels.at(s)) = &no_ch_doub_vec;
                     m_InterceptResX_Ref.at(not_channels.at(s)) = &no_ch_doub_vec;
                     m_InterceptResY_Ref.at(not_channels.at(s)) = &no_ch_doub_vec;
                    }
               }
            for (unsigned int dut = 0; dut < (m_fileDUTs.at(file)).size(); dut++)
                {
                 unsigned int dutindx = 999;
                 for (unsigned int s = 0; s < m_ntp_DUT.size(); s++) { if (m_ntp_DUT.at(s) == (m_fileDUTs.at(file)).at(dut)) { dutindx = s; break; } }
                 if (dutindx == 999)
                    {
                     std::cout << __FUNCTION__ << " ERROR: Could not match DUT " << (m_fileDUTs.at(file)).at(dut) << " for this run, skipping...." << std::endl;
                     return false;
                    }
                 (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("InterceptX_DUT%u", (m_fileDUTs.at(file)).at(dut)), m_InterceptX_DUT.at(dutindx), &b_InterceptX_DUT.at(dutindx));
                 (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("InterceptY_DUT%u", (m_fileDUTs.at(file)).at(dut)), m_InterceptY_DUT.at(dutindx), &b_InterceptY_DUT.at(dutindx));
                 (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("InterceptResX_DUT%u", (m_fileDUTs.at(file)).at(dut)), m_InterceptResX_DUT.at(dutindx), &b_InterceptResX_DUT.at(dutindx));
                 (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("InterceptResY_DUT%u", (m_fileDUTs.at(file)).at(dut)), m_InterceptResY_DUT.at(dutindx), &b_InterceptResY_DUT.at(dutindx));
                }
            // Non existant DUTs need to be populated correctly
            not_channels.clear();
            for (unsigned int p = 0; p < m_ntp_DUT.size(); p++)
                { 
                 found = false;
                 for (unsigned int a = 0; a < (m_fileDUTs.at(file)).size(); a++)
                     {
                      if (m_ntp_DUT.at(p) == (m_fileDUTs.at(file)).at(a))
                         {
                          found = true;
                          break;
                         }
                     }
                 if (!found) not_channels.push_back(p);
                }
            if (not_channels.size() > 0)
               { 
                for (unsigned int s = 0; s < not_channels.size(); s++)
                    {
                     b_InterceptX_DUT.at(not_channels.at(s)) = nullptr;
                     b_InterceptY_DUT.at(not_channels.at(s)) = nullptr;
                     b_InterceptResX_DUT.at(not_channels.at(s)) = nullptr;
                     b_InterceptResY_DUT.at(not_channels.at(s)) = nullptr;
                     m_InterceptX_DUT.at(not_channels.at(s)) = &no_ch_doub_vec;
                     m_InterceptY_DUT.at(not_channels.at(s)) = &no_ch_doub_vec;
                     m_InterceptResX_DUT.at(not_channels.at(s)) = &no_ch_doub_vec;
                     m_InterceptResY_DUT.at(not_channels.at(s)) = &no_ch_doub_vec;
                    }
               }
            for (unsigned int pln = 0; pln < (m_fileplanes.at(file)).size(); pln++)
                {
                 unsigned int plnindx = 999;
                 for (unsigned int s = 0; s < m_ntp_planes.size(); s++) { if (m_ntp_planes.at(s) == (m_fileplanes.at(file)).at(pln)) { plnindx = s; break; } }
                 if (plnindx == 999)
                    {
                     std::cout << __FUNCTION__ << " ERROR: Could not match telescope plane " << (m_fileplanes.at(file)).at(pln) << " on this run, skipping...." << std::endl;
                     return false;
                    }
                 (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("InterceptX_TelPlane%u", (m_fileplanes.at(file)).at(pln)), m_InterceptX_TelPlane.at(plnindx), &b_InterceptX_TelPlane.at(plnindx));
                 (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("InterceptY_TelPlane%u", (m_fileplanes.at(file)).at(pln)), m_InterceptY_TelPlane.at(plnindx), &b_InterceptY_TelPlane.at(plnindx));
                 (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("InterceptResX_TelPlane%u", (m_fileplanes.at(file)).at(pln)), m_InterceptResX_TelPlane.at(plnindx), &b_InterceptResX_TelPlane.at(plnindx));
                 (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("InterceptResY_TelPlane%u", (m_fileplanes.at(file)).at(pln)), m_InterceptResY_TelPlane.at(plnindx), &b_InterceptResY_TelPlane.at(plnindx));
                 (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("ClusterSizeX_TelPlane%u", (m_fileplanes.at(file)).at(pln)), m_ClusterSizeX_TelPlane.at(plnindx), &b_ClusterSizeX_TelPlane.at(plnindx));
                 (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("ClusterSizeY_TelPlane%u", (m_fileplanes.at(file)).at(pln)), m_ClusterSizeY_TelPlane.at(plnindx), &b_ClusterSizeY_TelPlane.at(plnindx));
                 (m_RunBase->WfmBase::GetRootTree(file))->SetBranchAddress(Form("ClusterNumPixels_TelPlane%u", (m_fileplanes.at(file)).at(pln)), m_ClusterNumPixels_TelPlane.at(plnindx), &b_ClusterNumPixels_TelPlane.at(plnindx));
                }
            // Non existant telescope planes need to be populated correctly
            not_channels.clear();
            for (unsigned int p = 0; p < m_ntp_planes.size(); p++)
                { 
                 found = false;
                 for (unsigned int a = 0; a < (m_fileplanes.at(file)).size(); a++)
                     {
                      if (m_ntp_planes.at(p) == (m_fileplanes.at(file)).at(a))
                         {
                          found = true;
                          break;
                         }
                     }
                 if (!found) not_channels.push_back(p);
                }
            if (not_channels.size() > 0)
               { 
                for (unsigned int s = 0; s < not_channels.size(); s++)
                    {
                     b_InterceptX_TelPlane.at(not_channels.at(s)) = nullptr;
                     b_InterceptY_TelPlane.at(not_channels.at(s)) = nullptr;
                     b_InterceptResX_TelPlane.at(not_channels.at(s)) = nullptr;
                     b_InterceptResY_TelPlane.at(not_channels.at(s)) = nullptr;
                     b_ClusterSizeX_TelPlane.at(not_channels.at(s)) = nullptr;
                     b_ClusterSizeY_TelPlane.at(not_channels.at(s)) = nullptr;
                     b_ClusterNumPixels_TelPlane.at(not_channels.at(s)) = nullptr;
                     m_InterceptX_TelPlane.at(not_channels.at(s)) = &no_ch_doub_vec;
                     m_InterceptY_TelPlane.at(not_channels.at(s)) = &no_ch_doub_vec;
                     m_InterceptResX_TelPlane.at(not_channels.at(s)) = &no_ch_doub_vec;
                     m_InterceptResY_TelPlane.at(not_channels.at(s)) = &no_ch_doub_vec;
                     m_ClusterSizeX_TelPlane.at(not_channels.at(s)) = &no_ch_doub_vec;
                     m_ClusterSizeY_TelPlane.at(not_channels.at(s)) = &no_ch_doub_vec;
                     m_ClusterNumPixels_TelPlane.at(not_channels.at(s)) = &no_ch_int_vec;
                     }
               }
           }
   else if (m_RunBase->WfmBase::GetTrackPackage() == "TBA" && m_filetrack.at(file) == 3)
           {
            // Add TBA relevent code
           }

   // Second stage analysis branches
   if (m_filelevels.at(file) > 1)
      {
       // Chi2 for all fits
       (m_ChTrees.at(file))->SetBranchAddress("chName", &m_chName, &b_chName);
       if (m_filelevels.at(file) == 3 || m_filelevels.at(file) == 5) 
           (m_ChTrees.at(file))->SetBranchAddress("SubChName", &m_SubChName, &b_SubChName);
       (m_ChTrees.at(file))->SetBranchAddress("MaxIndxFtChi2", &m_MaxIndxFtChi2, &b_MaxIndxFtChi2);
       (m_ChTrees.at(file))->SetBranchAddress("MinIndxFtChi2", &m_MinIndxFtChi2, &b_MinIndxFtChi2);
       (m_ChTrees.at(file))->SetBranchAddress("MaxVoltFtChi2", &m_MaxVoltFtChi2, &b_MaxVoltFtChi2);
       (m_ChTrees.at(file))->SetBranchAddress("MinVoltFtChi2", &m_MinVoltFtChi2, &b_MinVoltFtChi2);
       (m_ChTrees.at(file))->SetBranchAddress("StrIndxFtChi2", &m_StrIndxFtChi2, &b_StrIndxFtChi2);
       (m_ChTrees.at(file))->SetBranchAddress("EndIndxFtChi2", &m_EndIndxFtChi2, &b_EndIndxFtChi2);
       (m_ChTrees.at(file))->SetBranchAddress("NoiseFtChi2", &m_NoiseFtChi2, &b_NoiseFtChi2);
       (m_ChTrees.at(file))->SetBranchAddress("NoiseErrFtChi2", &m_NoiseErrFtChi2, &b_NoiseErrFtChi2);
       (m_ChTrees.at(file))->SetBranchAddress("PedestFtChi2", &m_PedestFtChi2, &b_PedestFtChi2);
       (m_ChTrees.at(file))->SetBranchAddress("PedestErrFtChi", &m_PedestErrFtChi2, &b_PedestErrFtChi2);
       (m_ChTrees.at(file))->SetBranchAddress("MaxTimeFtChi2", &m_MaxTimeFtChi2, &b_MaxTimeFtChi2);
       (m_ChTrees.at(file))->SetBranchAddress("MinTimeFtChi2", &m_MinTimeFtChi2, &b_MinTimeFtChi2);
       (m_ChTrees.at(file))->SetBranchAddress("ChargeFtChi2", &m_ChargeFtChi2, &b_ChargeFtChi2);
       (m_ChTrees.at(file))->SetBranchAddress("RiseTimeFtChi2", &m_RiseTimeFtChi2, &b_RiseTimeFtChi2);
       (m_ChTrees.at(file))->SetBranchAddress("TriggTimeFtChi2", &m_TriggTimeFtChi2, &b_TriggTimeFtChi2);
       (m_ChTrees.at(file))->SetBranchAddress("DVDTMaxFtChi2", &m_DVDTMaxFtChi2, &b_DVDTMaxFtChi2);
       (m_ChTrees.at(file))->SetBranchAddress("TriggToTFtChi2", &m_TriggToTFtChi2, &b_TriggToTFtChi2);
       (m_ChTrees.at(file))->SetBranchAddress("SignalFFTFtChi2", &m_SignalFFTFtChi2, &b_SignalFFTFtChi2);
       (m_ChTrees.at(file))->SetBranchAddress("NoiseFFTFtChi2", &m_NoiseFFTFtChi2, &b_NoiseFFTFtChi2);
       // Fit parameters for each channel
       (m_ChTrees.at(file))->SetBranchAddress("ChFitMaxVolt", &m_ChFitMxVlt.first, &b_ChFitMaxVolt);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitMaxVoltErr", &m_ChFitMxVlt.second, &b_ChFitMaxVoltErr);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitMinVolt", &m_ChFitMinVlt.first, &b_ChFitMinVolt);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitMinVoltErr", &m_ChFitMinVlt.second, &b_ChFitMinVoltErr);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitMxIndx", &m_ChFitMxIndx.first, &b_ChFitMxIndx);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitMxIndxErr", &m_ChFitMxIndx.second, &b_ChFitMxIndxErr);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitMinIndx", &m_ChFitMinIndx.first, &b_ChFitMinIndx);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitMinIndxErr", &m_ChFitMinIndx.second, &b_ChFitMinIndxErr);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitStrIndx", &m_ChFitStrIndx.first, &b_ChFitStrIndx);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitStrIndxErr", &m_ChFitStrIndx.second, &b_ChFitStrIndxErr);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitEndIndx", &m_ChFitEndIndx.first, &b_ChFitEndIndx);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitEndIndxErr", &m_ChFitEndIndx.second, &b_ChFitEndIndxErr);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitIsSig", &m_ChFitIsSig.first, &b_ChFitIsSig);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitIsSigErr", &m_ChFitIsSig.second, &b_ChFitIsSigErr);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitInWin", &m_ChFitInWin.first, &b_ChFitInWin);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitInWinErr", &m_ChFitInWin.second, &b_ChFitInWinErr);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitVStr", &m_ChFitVStr.first, &b_ChFitVStr);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitVStrErr", &m_ChFitVStr.second, &b_ChFitVStrErr);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitNoise", &m_ChFitNoise.first, &b_ChFitNoise);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitNoiseUnc", &m_ChFitNoise.second, &b_ChFitNoiseUnc);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitPedestal", &m_ChFitPedestal.first, &b_ChFitPedestal);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitPedestalUnc", &m_ChFitPedestal.second, &b_ChFitPedestalUnc);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitNoiseSlp", &m_ChFitNoiseSlp.first, &b_ChFitNoiseSlp);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitNoiseSlpErr", &m_ChFitNoiseSlp.second, &b_ChFitNoiseSlpErr);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitPedestalSlp", &m_ChFitPedestalSlp.first, &b_ChFitPedestalSlp);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitPedestalSlpErr", &m_ChFitPedestalSlp.second, &b_ChFitPedestalSlpErr);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitNoiseErr", &m_ChFitNoiseErr.first, &b_ChFitNoiseErr);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitNoiseErrUnc", &m_ChFitNoiseErr.second, &b_ChFitNoiseErrUnc);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitPedestalErr", &m_ChFitPedestalErr.first, &b_ChFitPedestalErr);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitPedestalErrUnc", &m_ChFitPedestalErr.second, &b_ChFitPedestalErrUnc);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitMaxTime", &m_ChFitMaxTime.first, &b_ChFitMaxTime);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitMaxTimeErr", &m_ChFitMaxTime.second, &b_ChFitMaxTimeErr);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitMinTime", &m_ChFitMinTime.first, &b_ChFitMinTime);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitMinTimeErr", &m_ChFitMinTime.second, &b_ChFitMinTimeErr);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitCharge", &m_ChFitCharge.first, &b_ChFitCharge);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitChargeErr", &m_ChFitCharge.second, &b_ChFitChargeErr);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitRiseT", &m_ChFitRiseT.first, &b_ChFitRiseT);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitRiseTErr", &m_ChFitRiseT.second, &b_ChFitRiseTErr);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitJitRiseSNR", &m_ChFitJitRiseSNR.first, &b_ChFitJitRiseSNR);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitJitRiseSNRErr", &m_ChFitJitRiseSNR.second, &b_ChFitJitRiseSNRErr);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitTriggTime", &m_ChFitTriggTime.first, &b_ChFitTriggTime);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitTriggTimeErr", &m_ChFitTriggTime.second, &b_ChFitTriggTimeErr);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitDVDTMax", &m_ChFitDVDTMax.first, &b_ChFitDVDTMax);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitDVDTMaxErr", &m_ChFitDVDTMax.second, &b_ChFitDVDTMaxErr);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitTriggToT", &m_ChFitTriggToT.first, &b_ChFitTriggToT);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitTriggToTErr", &m_ChFitTriggToT.second, &b_ChFitTriggToTErr);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitComplete", &m_ChFitComplete.first, &b_ChFitComplete);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitCompleteErr", &m_ChFitComplete.second, &b_ChFitCompleteErr);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitSoNR", &m_ChFitSoNR.first, &b_ChFitSoNR);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitSoNRErr", &m_ChFitSoNR.second, &b_ChFitSoNRErr);
       // FFT Second stage
       if ((m_FFTShDef.at(file)).at(1) || (m_FFTShDef.at(file)).at(2))
          {
           (m_ChTrees.at(file))->SetBranchAddress("ChFitSignalFFT", &m_ChFitSignalFFT.first, &b_ChFitSignalFFT);
           (m_ChTrees.at(file))->SetBranchAddress("ChFitSignalFFTErr", &m_ChFitSignalFFT.second, &b_ChFitSignalFFTErr);
           (m_ChTrees.at(file))->SetBranchAddress("ChFitNoiseFFT", &m_ChFitNoiseFFT.first, &b_ChFitNoiseFFT);
           (m_ChTrees.at(file))->SetBranchAddress("ChFitNoiseFFTErr", &m_ChFitNoiseFFT.second, &b_ChFitNoiseFFT);
          }
       // CFD anf DV/DT vectors
       (m_ChTrees.at(file))->SetBranchAddress("ChFitCFDTime", &m_ChFitCFDTime, &b_ChFitCFDTime);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitDVDTCFD", &m_ChFitDVDTCFD, &b_ChFitDVDTCFD);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitCFDToT", &m_ChFitCFDToT, &b_ChFitCFDToT);
       (m_ChTrees.at(file))->SetBranchAddress("ChFitJitNdVdT", &m_ChFitJitNdVdT, &b_ChFitJitNdVdT);
       (m_ChTrees.at(file))->SetBranchAddress("CFDTimeFtChi2", &m_CFDTimeFtChi2, &b_CFDTimeFtChi2);
       (m_ChTrees.at(file))->SetBranchAddress("DVDTCFDFtChi2", &m_DVDTCFDFtChi2, &b_DVDTCFDFtChi2);
       (m_ChTrees.at(file))->SetBranchAddress("CFDToTFtChi2", &m_CFDToTFtChi2, &b_CFDToTFtChi2);
       // Shape paramters for lvl2
       if ((m_FFTShDef.at(file)).at(5) || (m_FFTShDef.at(file)).at(6))
          {
           (m_ShapeTrees.at(file))->SetBranchAddress("ChShapeName", &m_ChShapeName, &b_ChShapeName);
           if ((m_FFTShDef.at(file)).at(6)) (m_ShapeTrees.at(file))->SetBranchAddress("SubChShName", &m_SubChShName, &b_SubChShName);
           (m_ShapeTrees.at(file))->SetBranchAddress("ChShapePolarity", &m_ChShapePolarity, &b_ChShapePolarity);
           (m_ShapeTrees.at(file))->SetBranchAddress("ChShapeMaxIndx", &m_ChShapeMaxIndx, &b_ChShapeMaxIndx);
           (m_ShapeTrees.at(file))->SetBranchAddress("ChShapeMinIndx", &m_ChShapeMinIndx, &b_ChShapeMinIndx);
           (m_ShapeTrees.at(file))->SetBranchAddress("ChShapeStrIndx", &m_ChShapeStrIndx, &b_ChShapeStrIndx);
           (m_ShapeTrees.at(file))->SetBranchAddress("ChShapeEndIndx", &m_ChShapeEndIndx, &b_ChShapeEndIndx);
           (m_ShapeTrees.at(file))->SetBranchAddress("ChShapeMax", &m_ChShapeMax, &b_ChShapeMax);
           (m_ShapeTrees.at(file))->SetBranchAddress("ChShapeMin", &m_ChShapeMin, &b_ChShapeMin);
           (m_ShapeTrees.at(file))->SetBranchAddress("ChShapeMaxTime", &m_ChShapeMaxTime, &b_ChShapeMaxTime);
           (m_ShapeTrees.at(file))->SetBranchAddress("ChShapeMinTime", &m_ChShapeMinTime, &b_ChShapeMinTime);
           (m_ShapeTrees.at(file))->SetBranchAddress("ChShapeCharge", &m_ChShapeCharge, &b_ChShapeCharge);
           (m_ShapeTrees.at(file))->SetBranchAddress("ChShapeRiseTime", &m_ChShapeRiseTime, &b_ChShapeRiseTime);
           (m_ShapeTrees.at(file))->SetBranchAddress("ChShapeJitRiseSNR", &m_ChShapeJitRiseSNR, &b_ChShapeJitRiseSNR);
           (m_ShapeTrees.at(file))->SetBranchAddress("ChShapeTriggTime", &m_ChShapeTriggTime, &b_ChShapeTriggTime);
           (m_ShapeTrees.at(file))->SetBranchAddress("ChShapeDVDTMax", &m_ChShapeDVDTMax, &b_ChShapeDVDTMax);
           (m_ShapeTrees.at(file))->SetBranchAddress("ChShapeTriggToT", &m_ChShapeTriggToT, &b_ChShapeTriggToT);
           if ((m_FFTShDef.at(file)).at(3) || (m_FFTShDef.at(file)).at(4))
              {
               (m_ShapeTrees.at(file))->SetBranchAddress("ChShapeSignalFFT", &m_ChShapeSignalFFT, &b_ChShapeSignalFFT);
               (m_ShapeTrees.at(file))->SetBranchAddress("ChShapeNoiseFFT", &m_ChShapeNoiseFFT, &b_ChShapeNoiseFFT);
              }
           (m_ShapeTrees.at(file))->SetBranchAddress("ChShapeNoise", &m_ChShapeNoise.first, &b_ChShapeNoise);
           (m_ShapeTrees.at(file))->SetBranchAddress("ChShapeNoiseErr", &m_ChShapeNoise.second, &b_ChShapeNoiseErr);
           (m_ShapeTrees.at(file))->SetBranchAddress("ChShapeSoN", &m_ChShapeSoN.first, &b_ChShapeSoN);
           (m_ShapeTrees.at(file))->SetBranchAddress("ChShapeSoNUnc", &m_ChShapeSoN.second, &b_ChShapeSoNUnc);
           (m_ShapeTrees.at(file))->SetBranchAddress("ChShapeCFDTime", &m_ChShapeCFDTime, &b_ChShapeCFDTime);
           (m_ShapeTrees.at(file))->SetBranchAddress("ChShapeDVDTCFD", &m_ChShapeDVDTCFD, &b_ChShapeDVDTCFD);
           (m_ShapeTrees.at(file))->SetBranchAddress("ChShapeCFDToT", &m_ChShapeCFDToT, &b_ChShapeCFDToT);
           (m_ShapeTrees.at(file))->SetBranchAddress("ChShapeJitNdVdT", &m_ChShapeJitNdVdT, &b_ChShapeJitNdVdT);
           (m_ShapeTrees.at(file))->SetBranchAddress("ChMeanVolt", &m_ChMeanVolt, &b_ChMeanVolt);
          }

       // Initialize 3rd stage analysis pointers
       if (m_filelevels.at(file) > 3)
          {
           (m_TimeResTrees.at(file))->SetBranchAddress("TmDiffChi2", &m_TmDiffChi2, &b_TmDiffChi2);
           (m_TimeResTrees.at(file))->SetBranchAddress("TmDiffFitMean", &m_TmDiffFitMean, &b_TmDiffFitMean);
           (m_TimeResTrees.at(file))->SetBranchAddress("TmDiffFitErr", &m_TmDiffFitErr, &b_TmDiffFitErr);
           (m_TimeResTrees.at(file))->SetBranchAddress("CFD_ChA", &m_CFD_ChA, &b_CFD_ChA);
           (m_TimeResTrees.at(file))->SetBranchAddress("ChA", &m_ChA, &b_ChA);
           (m_TimeResTrees.at(file))->SetBranchAddress("ChB", &m_ChB, &b_ChB);
          }
      }

   return true;
}
// --------------------------------------------------------------------------------------------------------------
// Function to be called in every event to fetch the per event parameteres for processing - copying.
bool WfmRun::GetEvntBranches(Long64_t entry, std::vector<bool> FFTShape, unsigned int filelvl, unsigned int trck, unsigned int fileindx)
{
    // Get these branches if FileLVL = 0 and analysis stage > 0 (always), copying or not: filelvl == 0 && m_RunBase->WfmBase::GetAnaStage() > 0
    //                    if FileLVL > 0 and analysis stage = 1 (redo first stage analysis), copying or not:  filelvl > 0 && m_RunBase->WfmBase::GetAnaStage() == 1
    //                    if copying to a new file independent of analysis stage and FileLVL: m_RunBase->WfmBase::GetInRootFile(fileindx) != m_RunBase->WfmBase::GetOutRootFile()
    //                    if FFT is requested but no FFT is present at first stage independent of analysis stage and file level, copying or not: m_RunBase->WfmBase::GetDoFFT() && !FFTShape.at(0)
    //                    if WaveShape is requeested but no WaveeShape is presentm indipendent of the level of the waveshape requestt
    if ((filelvl == 0 && m_RunBase->WfmBase::GetAnaStage() > 0) || (filelvl > 0 && m_RunBase->WfmBase::GetAnaStage() == 1) ||
        (m_RunBase->WfmBase::GetInRootFile(fileindx) != m_RunBase->WfmBase::GetOutRootFile()) || (m_RunBase->WfmBase::GetDoFFT() && !FFTShape.at(0)) ||
        (m_RunBase->WfmBase::GetDoWaveShape() && (!FFTShape.at(5) || !FFTShape.at(6))))
       {
        for (unsigned int ich = 0; ich < (m_RunBase->WfmBase::GetChannels()).size(); ich++)
            {
             if (b_t.at(ich) != nullptr)
                {
                 b_t.at(ich)->GetEntry(entry);
                 b_w.at(ich)->GetEntry(entry);
                 if (m_RunBase->WfmBase::GetInstrument() == Sampic)
                    {
                     b_ordrt.at(ich)->GetEntry(entry);
                     b_physt.at(ich)->GetEntry(entry);
                    }
                 if (m_RunBase->WfmBase::GetInstrument() == TestBeamBin1 || m_RunBase->WfmBase::GetInstrument() == TestBeamBin2 || m_RunBase->WfmBase::GetInstrument() == LabTXT ||
                     m_RunBase->WfmBase::GetInstrument() == TektronixScope || m_RunBase->WfmBase::GetInstrument() == LeCroyWRBin)
                    {
                     b_nPoints.at(ich)->GetEntry(entry);
                     b_SnRate.at(ich)->GetEntry(entry);
                     b_vScale.at(ich)->GetEntry(entry);
                     if (m_RunBase->WfmBase::GetInstrument() == TektronixScope) b_EvnTimeCh.at(ich)->GetEntry(entry);
                     if (m_RunBase->WfmBase::GetInstrument() == TestBeamBin1 || m_RunBase->WfmBase::GetInstrument() == TestBeamBin2 ||
                         m_RunBase->WfmBase::GetInstrument() == LeCroyWRBin) b_Scope.at(ich)->GetEntry(entry);
                    }
                }
            }
        if (m_RunBase->WfmBase::GetInstrument() == Sampic)
           {
            b_nPoints.at(0)->GetEntry(entry);
            b_SnRate.at(0)->GetEntry(entry);
           }
        if (m_RunBase->WfmBase::GetInstrument() == LabTXT || m_RunBase->WfmBase::GetInstrument() == TestBeamBin2) b_EvnTime->GetEntry(entry);
        if (m_RunBase->WfmBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Basic branches retrieved!" << std::endl;
       }

    // Get these branches if FileLVL < 2 && Analysis Stage > 1, copying or not
    //                    if FileLVL > 1 && Analysis Stage == 2 (redo second stage), copying or not
    //                    if copying and fileLVL > 1
    //                    if m_RunBase->WfmBase::GetDoFFT() and !FFTShape.at(0) and analysis LVL > 1 indipendent of file level, copying or not
    if ((m_RunBase->WfmBase::GetAnaStage() > 1 && (filelvl == 1 || (m_RunBase->WfmBase::GetDoFFT() && !FFTShape.at(1)))) ||
        (filelvl > 1 && m_RunBase->WfmBase::GetAnaStage() == 2) || (m_RunBase->WfmBase::GetInRootFile(fileindx) != m_RunBase->WfmBase::GetOutRootFile() && filelvl > 1))
       { 
        for (unsigned int ich = 0; ich < (m_RunBase->WfmBase::GetChannels()).size(); ich++)
            {
             if (b_EvPol.at(ich) != nullptr)
                   {
                    b_EvPol.at(ich)->GetEntry(entry);
                    b_EvCharge.at(ich)->GetEntry(entry);
                    b_EvIsSignal.at(ich)->GetEntry(entry);
                    b_EvInWindow.at(ich)->GetEntry(entry);
                    b_EvVoltStr.at(ich)->GetEntry(entry);
                    b_EvMaxindx.at(ich)->GetEntry(entry);
                    b_EvMaxTime.at(ich)->GetEntry(entry);
                    b_EvMaxVolt.at(ich)->GetEntry(entry);
                    b_EvMinIndx.at(ich)->GetEntry(entry);
                    b_EvMinTime.at(ich)->GetEntry(entry);
                    b_EvMinVolt.at(ich)->GetEntry(entry);
                    b_EvStrIndx.at(ich)->GetEntry(entry);
                    b_EvEndIndx.at(ich)->GetEntry(entry);
                    b_EvNoise.at(ich)->GetEntry(entry);
                    b_EvPedest.at(ich)->GetEntry(entry);
                    b_EvNoiseErr.at(ich)->GetEntry(entry);
                    b_EvPedestErr.at(ich)->GetEntry(entry);
                    b_EvRiseT.at(ich)->GetEntry(entry);
                    b_EvJitRiseSNR.at(ich)->GetEntry(entry);
                    b_EvTriggTime.at(ich)->GetEntry(entry);
                    b_EvDVDTMax.at(ich)->GetEntry(entry);
                    b_EvTriggToT.at(ich)->GetEntry(entry);
                    b_EvComplete.at(ich)->GetEntry(entry);
                    if (FFTShape.at(0))
                       {
                        b_EvSignalFFT.at(ich)->GetEntry(entry);
                        b_EvNoiseFFT.at(ich)->GetEntry(entry);
                       }
                    b_EvVAdjCFD.at(ich)->GetEntry(entry);
                    b_EvCFDTime.at(ich)->GetEntry(entry);
                    b_EvDVDTCFD.at(ich)->GetEntry(entry);
                    b_EvCFDToT.at(ich)->GetEntry(entry);
                    b_EvJitNdVdT.at(ich)->GetEntry(entry);
                   }
            }
        b_FitDiffChi2->GetEntry(entry);
        b_FitDiffMag->GetEntry(entry);
        b_FitDiffMagErr->GetEntry(entry);
        b_FitDiffSigma->GetEntry(entry);
        b_FitDiffSigmaErr->GetEntry(entry);
        b_FitDiffTime->GetEntry(entry);
       }

    // Tracking branches
    // Get these branches if tracking exists for this file and the file is copied: m_RunBase->WfmBase::GetInRootFile(fileindx) != m_RunBase->WfmBase::GetOutRootFile() && trck != 4
    //                    if analysis stage > 1 and tracking efficiency requested: m_RunBase->WfmBase::GetAnaStage() > 1 && m_RunBase->WfmBase::GetTrckEff()  
   if ((m_RunBase->WfmBase::GetInRootFile(fileindx) != m_RunBase->WfmBase::GetOutRootFile() && trck != 4) || (m_RunBase->WfmBase::GetAnaStage() > 1 && m_RunBase->WfmBase::GetTrckEff()))
      {
       // EU Telescope branches
       if (m_RunBase->WfmBase::GetTrackPackage() == "EUTelescope")
          {
           b_plane_id->GetEntry(entry);
           b_track_id->GetEntry(entry);
           b_trigger_id->GetEntry(entry);
           b_timestampTr->GetEntry(entry);
           b_xTr->GetEntry(entry);
           b_yTr->GetEntry(entry);
           b_omegaTr->GetEntry(entry);
           b_phiTr->GetEntry(entry);
           b_kinkxTr->GetEntry(entry);
           b_kinkyTr->GetEntry(entry);
           b_chi2Tr->GetEntry(entry);
           b_ndofTr->GetEntry(entry);
           b_xHit->GetEntry(entry);
           b_yHit->GetEntry(entry);
           b_zHit->GetEntry(entry);
           b_ID->GetEntry(entry);
          }
       // Bonn TBA Tracking barnches
       if (m_RunBase->WfmBase::GetTrackPackage() == "TBA")
          {
           //  Add TBA relevent code
          }
       // Cory Branches
       if (m_RunBase->WfmBase::GetTrackPackage() == "Cory")
          {
           b_ClusterSizeX_ROI->GetEntry(entry);
           b_ClusterSizeY_ROI->GetEntry(entry);
           b_PixelX_ROI->GetEntry(entry);
           b_PixelY_ROI->GetEntry(entry);
           b_ClusterNumPixels_ROI->GetEntry(entry);
           b_InterceptX_ROI->GetEntry(entry);
           b_InterceptY_ROI->GetEntry(entry);
           b_InterceptResX_ROI->GetEntry(entry);
           b_InterceptResY_ROI->GetEntry(entry);
           b_HasAsso->GetEntry(entry);
           for (unsigned int ref = 0; ref < m_ntp_Ref.size(); ref++)
               {
                if (b_InterceptX_Ref.at(ref) != nullptr)
                   {
                    b_InterceptX_Ref.at(ref)->GetEntry(entry);
                    b_InterceptY_Ref.at(ref)->GetEntry(entry);
                    b_InterceptResX_Ref.at(ref)->GetEntry(entry);
                    b_InterceptResY_Ref.at(ref)->GetEntry(entry);
                   }
               }
           for (unsigned int dut = 0; dut < m_ntp_DUT.size(); dut++)
               {
                if (b_InterceptX_DUT.at(dut) != nullptr)
                   {
                    b_InterceptX_DUT.at(dut)->GetEntry(entry);
                    b_InterceptY_DUT.at(dut)->GetEntry(entry);
                    b_InterceptResX_DUT.at(dut)->GetEntry(entry);
                    b_InterceptResY_DUT.at(dut)->GetEntry(entry);
                   }
               }
           for (unsigned int pln = 0; pln < m_ntp_planes.size(); pln++)
               {
                if (b_InterceptX_TelPlane.at(pln) != nullptr)
                   {
                    b_InterceptX_TelPlane.at(pln)->GetEntry(entry);
                    b_InterceptY_TelPlane.at(pln)->GetEntry(entry);
                    b_InterceptResX_TelPlane.at(pln)->GetEntry(entry);
                    b_InterceptResY_TelPlane.at(pln)->GetEntry(entry);
                    b_ClusterSizeX_TelPlane.at(pln)->GetEntry(entry);
                    b_ClusterSizeY_TelPlane.at(pln)->GetEntry(entry);
                    b_ClusterNumPixels_TelPlane.at(pln)->GetEntry(entry);
                   }
               }
          }
       if (m_RunBase->WfmBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Tracking branches retrieved!" << std::endl;
      }

   return true;
}
// --------------------------------------------------------------------------------------------------------------
// Function to be called at the second or third stage analysis to fetch the corresponding parameters from the Ntuple.
// THis function need the stage as an argument and should be called for every channel at the second stage or for 
// every CFD percentage at the third stage.
bool WfmRun::GetChannelBranches(Long64_t channel, unsigned int fileindx, std::vector<bool> FFTShape, unsigned int filelvl)
{
   // Second stage analysis branches
   if (filelvl > 1)
      {         
       if ((m_RunBase->WfmBase::GetAnaStage() == 2 || m_RunBase->WfmBase::GetAnaStage() > 3) || 
           (m_RunBase->WfmBase::GetInRootFile(fileindx) != m_RunBase->WfmBase::GetOutRootFile()) ||
           (m_RunBase->WfmBase::GetDoFFT() && (FFTShape.at(5) && !FFTShape.at(3)) || (FFTShape.at(6) && !FFTShape.at(5))))
          { 
           b_chName->GetEntry(channel);
           if (filelvl == 3 || filelvl == 5) b_SubChName->GetEntry(channel);
           b_MaxIndxFtChi2->GetEntry(channel);
           b_MinIndxFtChi2->GetEntry(channel);
           b_MaxVoltFtChi2->GetEntry(channel);
           b_MinVoltFtChi2->GetEntry(channel);
           b_StrIndxFtChi2->GetEntry(channel);
           b_EndIndxFtChi2->GetEntry(channel);
           b_NoiseFtChi2->GetEntry(channel);
           b_NoiseErrFtChi2->GetEntry(channel);
           b_PedestFtChi2->GetEntry(channel);
           b_PedestErrFtChi2->GetEntry(channel);
           b_MaxTimeFtChi2->GetEntry(channel);
           b_MinTimeFtChi2->GetEntry(channel);
           b_ChargeFtChi2->GetEntry(channel);
           b_RiseTimeFtChi2->GetEntry(channel);
           b_TriggTimeFtChi2->GetEntry(channel);
           b_DVDTMaxFtChi2->GetEntry(channel);
           b_TriggToTFtChi2->GetEntry(channel);
           b_SignalFFTFtChi2->GetEntry(channel);
           b_NoiseFFTFtChi2->GetEntry(channel);
           // Fit parameters for each channel
           b_ChFitMaxVolt->GetEntry(channel);
           b_ChFitMaxVoltErr->GetEntry(channel);
           b_ChFitMinVolt->GetEntry(channel);
           b_ChFitMinVoltErr->GetEntry(channel);
           b_ChFitMxIndx->GetEntry(channel);
           b_ChFitMxIndxErr->GetEntry(channel);
           b_ChFitMinIndx->GetEntry(channel);
           b_ChFitMinIndxErr->GetEntry(channel);
           b_ChFitStrIndx->GetEntry(channel);
           b_ChFitStrIndxErr->GetEntry(channel);
           b_ChFitEndIndx->GetEntry(channel);
           b_ChFitEndIndxErr->GetEntry(channel);
           b_ChFitIsSig->GetEntry(channel);
           b_ChFitIsSigErr->GetEntry(channel);
           b_ChFitInWin->GetEntry(channel);
           b_ChFitInWinErr->GetEntry(channel);
           b_ChFitVStr->GetEntry(channel);
           b_ChFitVStrErr->GetEntry(channel);
           b_ChFitNoise->GetEntry(channel);
           b_ChFitNoiseUnc->GetEntry(channel);
           b_ChFitPedestal->GetEntry(channel);
           b_ChFitPedestalUnc->GetEntry(channel);
           b_ChFitNoiseSlp->GetEntry(channel);
           b_ChFitNoiseSlpErr->GetEntry(channel);
           b_ChFitPedestalSlp->GetEntry(channel);
           b_ChFitPedestalSlpErr->GetEntry(channel);
           b_ChFitNoiseErr->GetEntry(channel);
           b_ChFitNoiseErrUnc->GetEntry(channel);
           b_ChFitPedestalErr->GetEntry(channel);
           b_ChFitPedestalErrUnc->GetEntry(channel);
           b_ChFitMaxTime->GetEntry(channel);
           b_ChFitMaxTimeErr->GetEntry(channel);
           b_ChFitMinTime->GetEntry(channel);
           b_ChFitMinTimeErr->GetEntry(channel);
           b_ChFitCharge->GetEntry(channel);
           b_ChFitChargeErr->GetEntry(channel);
           b_ChFitRiseT->GetEntry(channel);
           b_ChFitRiseTErr->GetEntry(channel);
           b_ChFitJitRiseSNR->GetEntry(channel);
           b_ChFitJitRiseSNRErr->GetEntry(channel);
           b_ChFitTriggTime->GetEntry(channel);
           b_ChFitTriggTimeErr->GetEntry(channel);
           b_ChFitDVDTMax->GetEntry(channel);
           b_ChFitDVDTMaxErr->GetEntry(channel);
           b_ChFitTriggToT->GetEntry(channel);
           b_ChFitTriggToTErr->GetEntry(channel);
           b_ChFitComplete->GetEntry(channel);
           b_ChFitCompleteErr->GetEntry(channel);
           b_ChFitSoNR->GetEntry(channel);
           b_ChFitSoNRErr->GetEntry(channel);
           // FFT First stage
           if (FFTShape.at(1) || FFTShape.at(2))
              {
               b_ChFitSignalFFT->GetEntry(channel);
               b_ChFitSignalFFTErr->GetEntry(channel);
               b_ChFitNoiseFFT->GetEntry(channel);
               b_ChFitNoiseFFT->GetEntry(channel);
              }
           // CFD anf DV/DT vectors
           b_ChFitCFDTime->GetEntry(channel);
           b_ChFitDVDTCFD->GetEntry(channel);
           b_ChFitCFDToT->GetEntry(channel);
           b_ChFitJitNdVdT->GetEntry(channel);
           b_CFDTimeFtChi2->GetEntry(channel);
           b_DVDTCFDFtChi2->GetEntry(channel);
           b_CFDToTFtChi2->GetEntry(channel);
           // Shape paramters for lvl2
           if (FFTShape.at(5) || FFTShape.at(6))
              {
               b_ChShapeName->GetEntry(channel);
               if (FFTShape.at(5)) b_SubChShName->GetEntry(channel);
               b_ChShapePolarity->GetEntry(channel);
               b_ChShapeMaxIndx->GetEntry(channel);
               b_ChShapeMinIndx->GetEntry(channel);
               b_ChShapeStrIndx->GetEntry(channel);
               b_ChShapeEndIndx->GetEntry(channel);
               b_ChShapeMax->GetEntry(channel);
               b_ChShapeMin->GetEntry(channel);
               b_ChShapeMaxTime->GetEntry(channel);
               b_ChShapeMinTime->GetEntry(channel);
               b_ChShapeCharge->GetEntry(channel);
               b_ChShapeRiseTime->GetEntry(channel);
               b_ChShapeJitRiseSNR->GetEntry(channel);
               b_ChShapeTriggTime->GetEntry(channel);
               b_ChShapeDVDTMax->GetEntry(channel);
               b_ChShapeTriggToT->GetEntry(channel);
               if (FFTShape.at(3) || FFTShape.at(4))
                  {
                   b_ChShapeSignalFFT->GetEntry(channel);
                   b_ChShapeNoiseFFT->GetEntry(channel);
                  }
               b_ChShapeNoise->GetEntry(channel);
               b_ChShapeNoiseErr->GetEntry(channel);
               b_ChShapeSoN->GetEntry(channel);
               b_ChShapeSoNUnc->GetEntry(channel);
               b_ChShapeCFDTime->GetEntry(channel);
               b_ChShapeDVDTCFD->GetEntry(channel);
               b_ChShapeCFDToT->GetEntry(channel);
               b_ChShapeJitNdVdT->GetEntry(channel);
               b_ChMeanVolt->GetEntry(channel);
              }
           return true;
          }
       // Third stage analysis branches, get per CFD percentage parameters
       else if ((filelvl > 3 && m_RunBase->WfmBase::GetAnaStage() > 3) || 
                (m_RunBase->WfmBase::GetInRootFile(fileindx) != m_RunBase->WfmBase::GetOutRootFile()))
               {
                b_TmDiffChi2->GetEntry(channel);
                b_TmDiffFitMean->GetEntry(channel);
                b_TmDiffFitErr->GetEntry(channel);
                b_CFD_ChA->GetEntry(channel);
                b_ChA->GetEntry(channel);
                b_ChB->GetEntry(channel);
                return true;
               }
       else return false;
      }
   else return false;
}
// --------------------------------------------------------------------------------------------------------------
// Function to be called at the begining after InitNtuple funvtion, to setup the required amount of pointers
// for the output file. Weather this will be just an addition on an existing file or a complete copy is handled internaly.
bool WfmRun::InitOutputBranches()
{
    bool newfile = (m_cp == -99 || (m_cp != -99 && m_RunBase->WfmBase::GetInRootFiles().size() > 1));
    // Typically  following branches always exist in input filse. Only reason to re-intialise them for output is either
    // creating a new file or appending to an existing file a series of input files
    if (newfile)
       {
        // Initialize the branch pointers
        b_t_out.clear();
        b_w_out.clear();
        b_t_out.resize((m_RunBase->WfmBase::GetChannels()).size());
        b_w_out.resize((m_RunBase->WfmBase::GetChannels()).size());
        if (m_RunBase->WfmBase::GetInstrument() == Sampic)
           {
            b_ordrt_out.clear();
            b_physt_out.clear();
            b_nPoints_out.clear();
            b_SnRate_out.clear();
            b_ordrt_out.resize((m_RunBase->WfmBase::GetChannels()).size());
            b_physt_out.resize((m_RunBase->WfmBase::GetChannels()).size());
            b_nPoints_out.resize(1, 0);
            b_SnRate_out.resize(1, 0);
           }
        if (m_RunBase->WfmBase::GetInstrument() == TestBeamBin1 || m_RunBase->WfmBase::GetInstrument() == TestBeamBin2 || m_RunBase->WfmBase::GetInstrument() == LabTXT ||
            m_RunBase->WfmBase::GetInstrument() == TektronixScope || m_RunBase->WfmBase::GetInstrument() == LeCroyWRBin)
           {
            b_nPoints_out.clear();
            b_SnRate_out.clear();
            b_vScale_out.clear();
            b_nPoints_out.resize((m_RunBase->WfmBase::GetChannels()).size());
            b_SnRate_out.resize((m_RunBase->WfmBase::GetChannels()).size());
            b_vScale_out.resize((m_RunBase->WfmBase::GetChannels()).size());
            if (m_RunBase->WfmBase::GetInstrument() == TektronixScope)
               {
                b_EvnTimeCh_out.clear();
                b_EvnTimeCh_out.resize((m_RunBase->WfmBase::GetChannels()).size());
               }
            if (m_RunBase->WfmBase::GetInstrument() == TestBeamBin1 || m_RunBase->WfmBase::GetInstrument() == TestBeamBin2 || m_RunBase->WfmBase::GetInstrument() == LeCroyWRBin)
               {
                b_Scope_out.clear();
                b_Scope_out.resize((m_RunBase->WfmBase::GetChannels()).size());
               }
           }
        if (WfmBase::GetVerbosity() >= 1) std::cout << __FUNCTION__ << " INFO: Initialization of level 0 branches complete, no of channels " << (m_RunBase->WfmBase::GetChannels()).size()
                                                    << " for instuemnet" << m_RunBase->WfmBase::GetInstrument() << " while newfile: " << newfile << std::endl;
       }

    // We need to split into cases:
    // - Creating new file, intialize if analysis level > 0 or high file level > 0 : newfile && (m_RunBase->WfmBase::GetAnaStage() > 0 || 
    //                                                                               m_RunBase->WfmBase::GetFileLVL().second > 0)
    // - Apending to existing file (!newfile) that either:
    //             1/ Has no first level branches while analysis level > 0: !newfile && m_RunBase->WfmBase::GetFileLVL().first == 0 && m_RunBase->WfmBase::GetAnaStage() > 0
    //             2/ Has no FFT while FFT is requested and analysis level > 0: !newfile && m_RunBase->WfmBase::GetDoFFT() && (m_RunBase->WfmBase::GetHasFFT().at(0)).first == 0 
    //                                                                           && m_RunBase->WfmBase::GetAnaStage() > 0
    if ((newfile && (m_RunBase->WfmBase::GetAnaStage() > 0 || m_RunBase->WfmBase::GetFileLVL().second > 0)) ||
        (!newfile && m_RunBase->WfmBase::GetAnaStage() > 0 && (m_RunBase->WfmBase::GetFileLVL().first == 0 || 
         (m_RunBase->WfmBase::GetDoFFT() && (m_RunBase->WfmBase::GetHasFFT().at(0)).first == 0))))
       {
        // If newfile OR if appending to a file AND have not first level branches while analysis level > 0
        if ((newfile && (m_RunBase->WfmBase::GetAnaStage() > 0 || m_RunBase->WfmBase::GetFileLVL().second > 0)) || 
            (!newfile && m_RunBase->WfmBase::GetFileLVL().first == 0 && m_RunBase->WfmBase::GetAnaStage() > 0))
           {
            b_EvPol_out.clear();
            b_EvCharge_out.clear();
            b_EvIsSignal_out.clear();
            b_EvInWindow_out.clear();
            b_EvVoltStr_out.clear();
            b_EvMaxindx_out.clear();
            b_EvMaxTime_out.clear();
            b_EvMaxVolt_out.clear();
            b_EvMinIndx_out.clear();
            b_EvMinTime_out.clear();
            b_EvMinVolt_out.clear();
            b_EvStrIndx_out.clear();
            b_EvEndIndx_out.clear();
            b_EvNoise_out.clear();
            b_EvPedest_out.clear();
            b_EvNoiseErr_out.clear();
            b_EvPedestErr_out.clear();
            b_EvRiseT_out.clear();
            b_EvJitNdVdT_out.clear();
            b_EvJitRiseSNR_out.clear();
            b_EvTriggTime_out.clear();
            b_EvDVDTMax_out.clear();
            b_EvDVDTCFD_out.clear();
            b_EvCFDToT_out.clear();
            b_EvTriggToT_out.clear();
            b_EvComplete_out.clear();
            b_EvVAdjCFD_out.clear();
            b_EvCFDTime_out.clear();
            b_EvPol_out.resize((m_RunBase->WfmBase::GetChannels()).size(), nullptr);
            b_EvCharge_out.resize((m_RunBase->WfmBase::GetChannels()).size(), nullptr);
            b_EvIsSignal_out.resize((m_RunBase->WfmBase::GetChannels()).size(), nullptr);
            b_EvInWindow_out.resize((m_RunBase->WfmBase::GetChannels()).size(), nullptr);
            b_EvVoltStr_out.resize((m_RunBase->WfmBase::GetChannels()).size(), nullptr);
            b_EvMaxindx_out.resize((m_RunBase->WfmBase::GetChannels()).size(), nullptr);
            b_EvMaxTime_out.resize((m_RunBase->WfmBase::GetChannels()).size(), nullptr);
            b_EvMaxVolt_out.resize((m_RunBase->WfmBase::GetChannels()).size(), nullptr);
            b_EvMinIndx_out.resize((m_RunBase->WfmBase::GetChannels()).size(), nullptr);
            b_EvMinTime_out.resize((m_RunBase->WfmBase::GetChannels()).size(), nullptr);
            b_EvMinVolt_out.resize((m_RunBase->WfmBase::GetChannels()).size(), nullptr);
            b_EvStrIndx_out.resize((m_RunBase->WfmBase::GetChannels()).size(), nullptr);
            b_EvEndIndx_out.resize((m_RunBase->WfmBase::GetChannels()).size(), nullptr);
            b_EvNoise_out.resize((m_RunBase->WfmBase::GetChannels()).size(), nullptr);
            b_EvPedest_out.resize((m_RunBase->WfmBase::GetChannels()).size(), nullptr);
            b_EvNoiseErr_out.resize((m_RunBase->WfmBase::GetChannels()).size(), nullptr);
            b_EvPedestErr_out.resize((m_RunBase->WfmBase::GetChannels()).size(), nullptr);
            b_EvRiseT_out.resize((m_RunBase->WfmBase::GetChannels()).size(), nullptr);
            b_EvJitNdVdT_out.resize((m_RunBase->WfmBase::GetChannels()).size(), nullptr);
            b_EvJitRiseSNR_out.resize((m_RunBase->WfmBase::GetChannels()).size(), nullptr);
            b_EvTriggTime_out.resize((m_RunBase->WfmBase::GetChannels()).size(), nullptr);
            b_EvDVDTMax_out.resize((m_RunBase->WfmBase::GetChannels()).size(), nullptr);
            b_EvDVDTCFD_out.resize((m_RunBase->WfmBase::GetChannels()).size(), nullptr);
            b_EvCFDToT_out.resize((m_RunBase->WfmBase::GetChannels()).size(), nullptr);
            b_EvTriggToT_out.resize((m_RunBase->WfmBase::GetChannels()).size(), nullptr);
            b_EvComplete_out.resize((m_RunBase->WfmBase::GetChannels()).size(), nullptr);
            b_EvVAdjCFD_out.resize((m_RunBase->WfmBase::GetChannels()).size(), nullptr);
            b_EvCFDTime_out.resize((m_RunBase->WfmBase::GetChannels()).size(), nullptr);
            if (WfmBase::GetVerbosity() >= 1) std::cout << __FUNCTION__ << " INFO: Initialization of level 1 branches complete, no of channels " << (m_RunBase->WfmBase::GetChannels()).size() 
                                                        << " DoFFT " << m_RunBase->WfmBase::GetDoFFT() << " HasFFT " << (m_RunBase->WfmBase::GetHasFFT().at(0)).first << " Analysis stage " 
                                                        << m_RunBase->WfmBase::GetAnaStage() << " File Level " << (m_RunBase->WfmBase::GetFileLVL()).first << std::endl;
            if (m_RunBase->WfmBase::GetDoFFT() || (((m_RunBase->WfmBase::GetHasFFT()).at(0)).second && newfile))
               {
                b_EvSignalFFT_out.clear();
                b_EvNoiseFFT_out.clear();
                b_EvSignalFFT_out.resize((m_RunBase->WfmBase::GetChannels()).size(), nullptr);
                b_EvNoiseFFT_out.resize((m_RunBase->WfmBase::GetChannels()).size(), nullptr);
                if (WfmBase::GetVerbosity() >= 1) std::cout << __FUNCTION__ << " INFO: FFT brancehs initialized! doFFT: " << m_RunBase->WfmBase::GetDoFFT() << std::endl;
               }
           }
        else {
              b_EvSignalFFT_out.clear();
              b_EvNoiseFFT_out.clear();
              b_EvSignalFFT_out.resize((m_RunBase->WfmBase::GetChannels()).size(), nullptr);
              b_EvNoiseFFT_out.resize((m_RunBase->WfmBase::GetChannels()).size(), nullptr);
              if (WfmBase::GetVerbosity() >= 1) std::cout << __FUNCTION__ << " INFO: Only FFT level 1 brancehs initialized! doFFT: " << m_RunBase->WfmBase::GetDoFFT() << std::endl;
             }
       }
    
   // Creatiing a new file, tracking exisit in input files and correctly identified tracking algorithm
   // Appending to existing file that has tracking, trakcing is cory and does not have all the planes, DUTs, Reference or ROI
   if (newfile && m_RunBase->WfmBase::GetTrackPackage() != "" && m_RunBase->WfmBase::GetHasTrck())
      { 
       if (WfmBase::GetVerbosity() >= 1) std::cout << __FUNCTION__ << " INFO: Tracking level branches intialization: " << m_RunBase->WfmBase::GetTrackPackage() 
                                                   << " Has tracking: " << m_RunBase->WfmBase::GetHasTrck() << " copying: " << newfile << std::endl;
       if (m_RunBase->WfmBase::GetTrackPackage() == "Cory")
          {
           // Reference Branches
           b_InterceptX_Ref_out.clear();
           b_InterceptY_Ref_out.clear();
           b_InterceptResX_Ref_out.clear();
           b_InterceptResY_Ref_out.clear();
           b_InterceptX_Ref_out.resize(m_nRecRefs, nullptr);
           b_InterceptY_Ref_out.resize(m_nRecRefs, nullptr);
           b_InterceptResX_Ref_out.resize(m_nRecRefs, nullptr);
           b_InterceptResY_Ref_out.resize(m_nRecRefs, nullptr);
           // DUT Branches
           b_InterceptX_DUT_out.clear();
           b_InterceptY_DUT_out.clear();
           b_InterceptResX_DUT_out.clear();
           b_InterceptResY_DUT_out.clear();
           b_InterceptX_DUT_out.resize(m_nRecDUTs, nullptr);
           b_InterceptY_DUT_out.resize(m_nRecDUTs, nullptr);
           b_InterceptResX_DUT_out.resize(m_nRecDUTs, nullptr);
           b_InterceptResY_DUT_out.resize(m_nRecDUTs, nullptr);
           // Telescope Plane Branches
           b_InterceptX_TelPlane_out.clear();
           b_InterceptY_TelPlane_out.clear();
           b_InterceptResX_TelPlane_out.clear();
           b_InterceptResY_TelPlane_out.clear();
           b_ClusterSizeX_TelPlane_out.clear();
           b_ClusterSizeY_TelPlane_out.clear();
           b_ClusterNumPixels_TelPlane_out.clear();
           b_InterceptX_TelPlane_out.resize(m_nRecPlanes, nullptr);
           b_InterceptY_TelPlane_out.resize(m_nRecPlanes, nullptr);
           b_InterceptResX_TelPlane_out.resize(m_nRecPlanes, nullptr);
           b_InterceptResY_TelPlane_out.resize(m_nRecPlanes, nullptr);
           b_ClusterSizeX_TelPlane_out.resize(m_nRecPlanes, nullptr);
           b_ClusterSizeY_TelPlane_out.resize(m_nRecPlanes, nullptr);
           b_ClusterNumPixels_TelPlane_out.resize(m_nRecPlanes, nullptr);
           if (WfmBase::GetVerbosity() >= 1) std::cout << __FUNCTION__ << " INFO: Finished cory branches initialization!" << std::endl;
          }
       else if (m_RunBase->WfmBase::GetTrackPackage() == "TBA")
               {
                //  TBA Branches should go here
               }
      }
   // No need for initialization of second and third stage pointers here, the structure of these trees with a channel as an event
   // automatically acounts for the varied number of channels per file, since one branch is globaly needed for each quntity.
   
   return true;
}
// --------------------------------------------------------------------------------------------------------------
bool WfmRun::SetOutputBranches()
{
   // Intial test used later on in the funciton
    bool newfile = (m_cp == -99 || (m_cp != -99 && m_RunBase->WfmBase::GetInRootFiles().size() > 1));

   // Only propagate event number if creating new file
    if (newfile)
       {
        b_EvnNo_out = m_OutTree->Branch("EvnNo", &EvnNo, "EvnNo/i");
        // If Fit Histograms exist clone first in the series, else create it
        if (m_RunBase->WfmBase::GetFileLVL().second > 0)
           {
            double m_FitDiffChi2;
            double m_FitDiffMag;
            double m_FitDiffMagErr;
            double m_FitDiffSigma;
            double m_FitDiffSigmaErr;
            double m_FitDiffTime;




           }
        else {

             }
       }
    else {
          // Copying on same file, recover esisting histogram to allow for modificaitons
          if (m_RunBase->WfmBase::GetFileLVL().first > 0)
             {
 
             }
          else {

               }
         }

   // Assigned branches in every case
   for (unsigned int chindx = 0; chindx < (m_RunBase->WfmBase::GetChannels()).size(); chindx++)
       {
        if (newfile)
           {
            if (m_RunBase->WfmBase::GetVerbosity() > 1)
               {
                std::cout << __FUNCTION__ << " INFO: Assigning lelvel 0 branches, copying file: " << newfile << std::endl;
                std::cout << "                                Processing channel " << chindx << " from " << (m_RunBase->WfmBase::GetChannels()).size() << std::endl;
               }
            if (m_RunBase->WfmBase::GetInstrument() == Sampic)
               {
                b_physt_out.at(chindx) = m_OutTree->Branch(Form("physt%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_physt.at(chindx), "m_physt.at(chindx)/D");
                b_ordrt_out.at(chindx) = m_OutTree->Branch(Form("ordrt%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_ordrt.at(chindx), "m_ordrt.at(chindx)/D");
               }
            b_t_out.at(chindx) = m_OutTree->Branch(Form("t%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &(m_t.at(chindx)));
            b_w_out.at(chindx) = m_OutTree->Branch(Form("w%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &(m_w.at(chindx)));
            if (m_RunBase->WfmBase::GetInstrument() == TestBeamBin1 || m_RunBase->WfmBase::GetInstrument() == TestBeamBin2 ||
                m_RunBase->WfmBase::GetInstrument() == LabTXT || m_RunBase->WfmBase::GetInstrument() == TektronixScope ||
                m_RunBase->WfmBase::GetInstrument() == LeCroyWRBin)
               {
                b_vScale_out.at(chindx) = m_OutTree->Branch(Form("vScale%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_scale.at(chindx), "m_scale.at(chindx)/F");
                b_nPoints_out.at(chindx) = m_OutTree->Branch(Form("nPoints%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_npoints.at(chindx), "m_npoints.at(chindx)/i");
                b_SnRate_out.at(chindx) = m_OutTree->Branch(Form("SnRate%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_srate.at(chindx), "m_srate.at(chindx)/L");
                if (m_RunBase->WfmBase::GetInstrument() == TektronixScope)
                    b_EvnTimeCh_out.at(chindx) = m_OutTree->Branch(Form("EvnTimeCh%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_EvnTimeCh.at(chindx), "m_EvnTimeCh.at(chindx)/D");
                if (m_RunBase->WfmBase::GetInstrument() == TestBeamBin1 || m_RunBase->WfmBase::GetInstrument() == TestBeamBin2)
                    b_Scope_out.at(chindx) = m_OutTree->Branch(Form("Scope%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_scope.at(chindx), "m_scope.at(chindx)/i");
               }
           }
        // Level 1 analysis requirements
        // Creating new file while file level > 0 or while analysis stage > 0
        // Not creating new file, file level 0 with analysis stage > 0 
        // Not creating new file but reuqested FFT while FFT not available and analysis stage > 0
        if ((newfile && (m_RunBase->WfmBase::GetAnaStage() > 0 || m_RunBase->WfmBase::GetFileLVL().second > 0)) ||
            (!newfile && m_RunBase->WfmBase::GetAnaStage() > 0 && (m_RunBase->WfmBase::GetFileLVL().first == 0 ||
                (m_RunBase->WfmBase::GetDoFFT() && (m_RunBase->WfmBase::GetHasFFT().at(0)).first == 0))))
           {
            if (!(!newfile && m_RunBase->WfmBase::GetAnaStage() > 0 && m_RunBase->WfmBase::GetFileLVL().first > 0 && 
                 m_RunBase->WfmBase::GetDoFFT() && (m_RunBase->WfmBase::GetHasFFT().at(0)).first == 0))
               {
                if (m_RunBase->WfmBase::GetVerbosity() > 1)
                  {
                    std::cout << __FUNCTION__ << " INFO: Assigning lelvel 1 branches, copying file: " << newfile << "\t Do FFT: " << m_RunBase->WfmBase::GetDoFFT() << std::endl;
                    std::cout << "                                                               file level: " << (m_RunBase->WfmBase::GetFileLVL()).second 
                              << "\t FFT presence in event level " << (m_RunBase->WfmBase::GetHasFFT().at(0)).first << std::endl;
                  }
                b_EvPol_out.at(chindx) = m_OutTree->Branch(Form("Polarity%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_EvPol.at(chindx), "m_EvPol.at(chindx)/I");
                b_EvCharge_out.at(chindx) = m_OutTree->Branch(Form("Charge%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_EvCharge.at(chindx), "m_EvCharge.at(chindx)/D");
                b_EvIsSignal_out.at(chindx) = m_OutTree->Branch(Form("IsSignal%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_EvIsSignal.at(chindx), "m_EvIsSignal.at(chindx)/I");
                b_EvInWindow_out.at(chindx) = m_OutTree->Branch(Form("InWindow%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_EvInWindow.at(chindx), "m_EvInWindow.at(chindx)/I");
                b_EvVoltStr_out.at(chindx) = m_OutTree->Branch(Form("VoltStr%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_EvVoltStr.at(chindx), "m_EvVoltStr.at(chindx)/I");
                b_EvMaxindx_out.at(chindx) = m_OutTree->Branch(Form("MaxIndx%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_EvMaxindx.at(chindx), "m_EvMaxindx.at(chindx)/I");
                b_EvMaxTime_out.at(chindx) = m_OutTree->Branch(Form("MaxTime%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_EvMaxTime.at(chindx), "m_EvMaxTime.at(chindx)/D");
                b_EvMaxVolt_out.at(chindx) = m_OutTree->Branch(Form("MaxVolt%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_EvMaxVolt.at(chindx), "m_EvMaxVolt.at(chindx)/D");
                b_EvMinIndx_out.at(chindx) = m_OutTree->Branch(Form("MinIndx%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_EvMinIndx.at(chindx), "m_EvMinIndx.at(chindx)/I");
                b_EvMinTime_out.at(chindx) = m_OutTree->Branch(Form("MinTime%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_EvMinTime.at(chindx), "m_MinTime.at(chindx)/D");
                b_EvMinVolt_out.at(chindx) = m_OutTree->Branch(Form("MinVolt%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_EvMinVolt.at(chindx), "m_EvMinVolt.at(chindx)/D");
                b_EvStrIndx_out.at(chindx) = m_OutTree->Branch(Form("Strindx%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_EvStrIndx.at(chindx), "m_EvStrIndx.at(chindx)/I");
                b_EvEndIndx_out.at(chindx) = m_OutTree->Branch(Form("EndIndx%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_EvEndIndx.at(chindx), "m_EvEndIndx.at(chindx)/I");
                b_EvNoise_out.at(chindx) = m_OutTree->Branch(Form("Noise%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_EvNoise.at(chindx), "m_EvNoise.at(chindx)/D");
                b_EvPedest_out.at(chindx) = m_OutTree->Branch(Form("Pedest%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_EvPedest.at(chindx), "m_EvPedest.at(chindx)/D");
                b_EvNoiseErr_out.at(chindx) = m_OutTree->Branch(Form("NoiseErr%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_EvNoiseErr.at(chindx), "m_EvNoiseErr.at(chindx)/D");
                b_EvPedestErr_out.at(chindx) = m_OutTree->Branch(Form("PedestErr%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_EvPedestErr.at(chindx), "m_EvPedestEr.at(chindx)/D");
                b_EvRiseT_out.at(chindx) = m_OutTree->Branch(Form("RiseT%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_EvRiseT.at(chindx), "m_EvRiseTe.at(chindx)/D");
                b_EvJitRiseSNR_out.at(chindx) = m_OutTree->Branch(Form("JitRiseSNR%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_EvJitRiseSNR.at(chindx), "m_EvJitRiseSNR.at(chindx)/D");
                b_EvTriggTime_out.at(chindx) = m_OutTree->Branch(Form("TrggTime%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_EvTriggTime.at(chindx), "m_EvTriggTime.at(chindx)/D");
                b_EvDVDTMax_out.at(chindx) = m_OutTree->Branch(Form("DVDTMax%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_EvDVDTMax.at(chindx), "m_EvDVDTMax.at(chindx)/D");
                b_EvTriggToT_out.at(chindx) = m_OutTree->Branch(Form("TriggToT%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_EvTriggToT.at(chindx), "m_EvTriggToT.at(chindx)/D");
                b_EvComplete_out.at(chindx) = m_OutTree->Branch(Form("Complete%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_EvComplete.at(chindx), "m_EvComplete.at(chindx)/I");
                b_EvVAdjCFD_out.at(chindx) = m_OutTree->Branch(Form("VAdjCFD%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_EvVAdjCFD.at(chindx));
                b_EvCFDTime_out.at(chindx) = m_OutTree->Branch(Form("CFDTime%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_EvCFDTime.at(chindx));
                b_EvDVDTCFD_out.at(chindx) = m_OutTree->Branch(Form("DVDTCFD%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_EvDVDTCFD.at(chindx));
                b_EvCFDToT_out.at(chindx) = m_OutTree->Branch(Form("CFDToT%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_EvCFDToT.at(chindx));
                b_EvJitNdVdT_out.at(chindx) = m_OutTree->Branch(Form("JitNdVdT%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_EvJitNdVdT.at(chindx));
                if (m_RunBase->WfmBase::GetDoFFT() || (((m_RunBase->WfmBase::GetHasFFT()).at(0)).second && newfile))
                   {
                    if (m_RunBase->WfmBase::GetVerbosity() > 1)
                       {
                        std::cout << __FUNCTION__ << " INFO: Assigning lelvel 1 FFT branches, copying file: " << newfile << "\t Do FFT: " << m_RunBase->WfmBase::GetDoFFT() << std::endl;
                        std::cout << "                                                                     file level: " << (m_RunBase->WfmBase::GetFileLVL()).second
                                  << "\t FFT presence at event level: " << (m_RunBase->WfmBase::GetHasFFT().at(0)).first << std::endl;
                       }
                    b_EvSignalFFT_out.at(chindx) = m_OutTree->Branch(Form("SignalFFT%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_EvSignalFFT.at(chindx), "m_EvSignalFFT.at(chindx)/D");
                    b_EvNoiseFFT_out.at(chindx) = m_OutTree->Branch(Form("NoiseFFT%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_EvNoiseFFT.at(chindx), "m_EvNoiseFFT.at(chindx)/D");
                   }
               }
            else {
                  if (m_RunBase->WfmBase::GetVerbosity() > 1)
                     {
                      std::cout << __FUNCTION__ << " INFO: Assigning lelvel 1 FFT branches, copying file: " << newfile << "\t Do FFT: " << m_RunBase->WfmBase::GetDoFFT() << std::endl;
                      std::cout << "                                                                     file level: " << (m_RunBase->WfmBase::GetFileLVL()).second
                                << "\t FFT presence at event level: " << (m_RunBase->WfmBase::GetHasFFT().at(0)).first << std::endl;
                     }
                  b_EvSignalFFT_out.at(chindx) = m_OutTree->Branch(Form("SignalFFT%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_EvSignalFFT.at(chindx), "m_EvSignalFFT.at(chindx)/D");
                  b_EvNoiseFFT_out.at(chindx) = m_OutTree->Branch(Form("NoiseFFT%02u", (m_RunBase->WfmBase::GetChannels()).at(chindx)), &m_EvNoiseFFT.at(chindx), "m_EvNoiseFFT.at(chindx)/D");
                 }
           }
       }

   // Assigned brancehs in all cases but non channel related
   if (m_RunBase->WfmBase::GetInstrument() == Sampic && newfile)
      {
       b_nPoints_out.at(0) = m_OutTree->Branch("nPoints", &m_npoints.at(0), "m_npoints.at(0)/i");
       b_SnRate_out.at(0) = m_OutTree->Branch("SnRate", &m_srate.at(0), "m_srate.at(0)/L");
      }
   if ((m_RunBase->WfmBase::GetInstrument() == LabTXT || m_RunBase->WfmBase::GetInstrument() == TestBeamBin2) && newfile)
       b_EvnTime_out = m_OutTree->Branch("EvnTime", &m_EvnTime, "m_EvnTime/D");

   // Tracking part of the file, only prppagated when copying. If tracking combination is happening at this stage, the TrackCombine functinos will create the necessery branches.
   if (newfile && m_RunBase->WfmBase::GetTrackPackage() != "" && m_RunBase->WfmBase::GetHasTrck())  // if we copy and either there is tracking in the original
      { 
       if (m_RunBase->WfmBase::GetTrackPackage() == "Cory")
          {
           b_ClusterSizeX_ROI_out = m_OutTree->Branch("ClusterSizeX_ROI", &m_ClusterSizeX_ROI);
           b_ClusterSizeY_ROI_out = m_OutTree->Branch("ClusterSizeY_ROI", &m_ClusterSizeY_ROI);
           b_PixelX_ROI_out = m_OutTree->Branch("PixelX_ROI", &m_PixelX_ROI);
           b_PixelY_ROI_out = m_OutTree->Branch("PixelY_ROI", &m_PixelY_ROI);
           b_ClusterNumPixels_ROI_out = m_OutTree->Branch("ClusterNumPixels_ROI", &m_ClusterNumPixels_ROI);
           b_InterceptX_ROI_out = m_OutTree->Branch("InterceptX_ROI", &m_InterceptX_ROI);
           b_InterceptY_ROI_out = m_OutTree->Branch("InterceptY_ROI", &m_InterceptY_ROI);
           b_InterceptResX_ROI_out = m_OutTree->Branch("InterceptResX_ROI", &m_InterceptResX_ROI);
           b_InterceptResY_ROI_out = m_OutTree->Branch("InterceptResY_ROI", &m_InterceptResY_ROI);
           b_HasAsso_out = m_OutTree->Branch("HasAsso", &m_HasAsso);
           for (unsigned int refsindx = 0; refsindx < m_ntp_Ref.size(); refsindx++)
               {
                b_InterceptX_Ref_out.at(refsindx) = m_OutTree->Branch(Form("InterceptX_Ref%u", m_ntp_Ref.at(refsindx)), m_InterceptX_Ref.at(refsindx));
                b_InterceptY_Ref_out.at(refsindx) = m_OutTree->Branch(Form("InterceptY_Ref%u", m_ntp_Ref.at(refsindx)), m_InterceptY_Ref.at(refsindx));
                b_InterceptResX_Ref_out.at(refsindx) = m_OutTree->Branch(Form("InterceptResX_Ref%u", m_ntp_Ref.at(refsindx)), m_InterceptResX_Ref.at(refsindx));
                b_InterceptResY_Ref_out.at(refsindx) = m_OutTree->Branch(Form("InterceptResY_Ref%u", m_ntp_Ref.at(refsindx)), m_InterceptResY_Ref.at(refsindx));
               }
           for (unsigned int dutindx = 0; dutindx < m_ntp_DUT.size(); dutindx++)
               {
                b_InterceptX_DUT_out.at(dutindx) = m_OutTree->Branch(Form("InterceptX_DUT%u", m_ntp_DUT.at(dutindx)), m_InterceptX_DUT.at(dutindx));
                b_InterceptY_DUT_out.at(dutindx) = m_OutTree->Branch(Form("InterceptY_DUT%u", m_ntp_DUT.at(dutindx)), m_InterceptY_DUT.at(dutindx));
                b_InterceptResX_DUT_out.at(dutindx) = m_OutTree->Branch(Form("InterceptResX_DUT%u", m_ntp_DUT.at(dutindx)), m_InterceptResX_DUT.at(dutindx));
                b_InterceptResY_DUT_out.at(dutindx) = m_OutTree->Branch(Form("InterceptResY_DUT%u", m_ntp_DUT.at(dutindx)), m_InterceptResY_DUT.at(dutindx));
               }
           for (unsigned int plindx = 0; plindx < m_ntp_planes.size(); plindx++)
               {
                b_InterceptX_TelPlane_out.at(plindx) = m_OutTree->Branch(Form("InterceptX_TelPlane%u", m_ntp_planes.at(plindx)), m_InterceptX_TelPlane.at(plindx));
                b_InterceptY_TelPlane_out.at(plindx) = m_OutTree->Branch(Form("InterceptY_TelPlane%u", m_ntp_planes.at(plindx)), m_InterceptY_TelPlane.at(plindx));
                b_InterceptResX_TelPlane_out.at(plindx) = m_OutTree->Branch(Form("InterceptResX_TelPlane%u", m_ntp_planes.at(plindx)), m_InterceptResX_TelPlane.at(plindx));
                b_InterceptResY_TelPlane_out.at(plindx) = m_OutTree->Branch(Form("InterceptResY_TelPlane%u", m_ntp_planes.at(plindx)), m_InterceptResY_TelPlane.at(plindx));
                b_ClusterSizeX_TelPlane_out.at(plindx) = m_OutTree->Branch(Form("ClusterSizeX_TelPlane%u", m_ntp_planes.at(plindx)), m_ClusterSizeX_TelPlane.at(plindx));
                b_ClusterSizeY_TelPlane_out.at(plindx) = m_OutTree->Branch(Form("ClusterSizeY_TelPlane%u", m_ntp_planes.at(plindx)), m_ClusterSizeY_TelPlane.at(plindx));
                b_ClusterNumPixels_TelPlane_out.at(plindx) = m_OutTree->Branch(Form("ClusterNumPixels_TelPlane%u", m_ntp_planes.at(plindx)), m_ClusterNumPixels_TelPlane.at(plindx));
               }
          }
       else if (m_RunBase->WfmBase::GetTrackPackage() == "EUTelescope")
               { 
                b_plane_id_out = m_OutTree->Branch("plane_id", &m_plane_id);
                b_track_id_out = m_OutTree->Branch("track_id", &m_track_id);
                b_trigger_id_out = m_OutTree->Branch("trigger_id", &m_trigger_id, "m_trigger_id/i");
                b_timestampTr_out = m_OutTree->Branch("timestampTr", &m_timestampTr, "timestampTr/i");
                b_xTr_out = m_OutTree->Branch("xTr", &m_xTr);
                b_yTr_out = m_OutTree->Branch("yTr", &m_yTr);
                b_omegaTr_out = m_OutTree->Branch("omegaTr", &m_omegaTr);
                b_phiTr_out = m_OutTree->Branch("phiTr", &m_phiTr);
                b_kinkxTr_out = m_OutTree->Branch("kinkxTr", &m_kinkxTr);
                b_kinkyTr_out = m_OutTree->Branch("kinkyTr", &m_kinkyTr);
                b_chi2Tr_out = m_OutTree->Branch("chi2Tr", &m_chi2Tr);
                b_ndofTr_out = m_OutTree->Branch("ndofTr", &m_ndofTr);
                b_xHit_out = m_OutTree->Branch("xHit", &m_xHit);
                b_yHit_out = m_OutTree->Branch("yHit", &m_yHit);
                b_zHit_out = m_OutTree->Branch("zHit", &m_zHit);
                b_ID_out = m_OutTree->Branch("ID", &m_ID);
               }
       else if (m_RunBase->WfmBase::GetTrackPackage() == "TBA")
               { 
                //  TBA Branches should go here
               }
       else {
             std::cout << "    " << __FUNCTION__ << " ERROR: Unknown reconstruction framework!" << std::endl;
             return false;
            }
      }

   // Creating new file while file level max > 1 irrelevant of analysis stage				
   // Creating new file, file level min < 2 while analysis stage > 1 					
   // Creating new file, file level min < 3 while analysis stage > 2 and analysis stage != 4	
   // Not creating new file but file level min < 2 while analysis stage > 1				
   // Not creating new file but file level min < 3 while analysis stage > 2 and analysis stage != 4	
   bool ChAna = (m_RunBase->WfmBase::GetAnaStage() == 2 || m_RunBase->WfmBase::GetAnaStage() == 4);
   bool SubChAna = (m_RunBase->WfmBase::GetAnaStage() == 3 || m_RunBase->WfmBase::GetAnaStage() == 5);
   bool lv2lv3 = ((newfile && m_RunBase->WfmBase::GetFileLVL().second > 1) || 
                  (m_RunBase->WfmBase::GetFileLVL().first < 2 && m_RunBase->WfmBase::GetAnaStage() > 1) ||
                  ((m_RunBase->WfmBase::GetFileLVL().first < 3 || m_RunBase->WfmBase::GetFileLVL().first == 4) &&
                  (m_RunBase->WfmBase::GetAnaStage() > 2 && m_RunBase->WfmBase::GetAnaStage() != 4)));
   // WaveShape provision for channels
   bool lv2Wsh = (ChAna && (m_RunBase->WfmBase::GetHasWaveShape().first != 1 || m_RunBase->WfmBase::GetHasWaveShape().first != 4) &&
                  m_RunBase->WfmBase::GetDoWaveShape());
   // WaveShape provision for sub-channels                                     
   bool lv3Wsh = (SubChAna && m_RunBase->WfmBase::GetDoWaveShape() && m_RunBase->WfmBase::GetHasWaveShape().first < 2);
   // FFT provision for channels
   bool lv2FFT = ((m_RunBase->WfmBase::GetAnaStage() > 1 && (m_RunBase->WfmBase::GetAnaStage() != 3 || m_RunBase->WfmBase::GetAnaStage() != 5)) && 
                  !((m_RunBase->WfmBase::GetHasFFT()).at(1)).first && m_RunBase->WfmBase::GetDoFFT());
   // FFT provision for subchannels
   bool lv3FFT = ((m_RunBase->WfmBase::GetAnaStage() > 2 && m_RunBase->WfmBase::GetAnaStage() != 4) && 
                  !((m_RunBase->WfmBase::GetHasFFT()).at(2)).first && m_RunBase->WfmBase::GetDoFFT());
   // FFT in WaveShape provision for channels
   bool lv2Wsh_FFT = ((((m_RunBase->WfmBase::GetHasWaveShape().first == 1 || m_RunBase->WfmBase::GetHasWaveShape().first == 3) && 
                        !(m_RunBase->WfmBase::GetHasFFT().at(3)).first && m_RunBase->WfmBase::GetDoFFT()) ||
                      ((m_RunBase->WfmBase::GetHasWaveShape().first == 0 || m_RunBase->WfmBase::GetHasWaveShape().first == 2) && 
                        m_RunBase->WfmBase::GetDoWaveShape() && m_RunBase->WfmBase::GetDoFFT())) && ChAna);
   // FFT in WaveShape provision for sub-channels
   bool lv3Wsh_FFT = (((m_RunBase->WfmBase::GetHasWaveShape().first > 1 && !(m_RunBase->WfmBase::GetHasFFT().at(3)).first && m_RunBase->WfmBase::GetDoFFT()) ||
                      (m_RunBase->WfmBase::GetHasWaveShape().first < 2 && m_RunBase->WfmBase::GetDoWaveShape() && m_RunBase->WfmBase::GetDoFFT())) && SubChAna);
   if (lv2lv3 || lv2Wsh || lv3Wsh || lv2FFT || lv3FFT || lv2Wsh_FFT || lv2Wsh_FFT)
      {
       if (lv2lv3)
          {
           // Chi2 for all fits
           b_chName_out = m_ChFitParam->Branch("chName", &m_chName, "m_chName/I");
           if (SubChAna || (newfile && (m_RunBase->WfmBase::GetFileLVL().second == 3 || m_RunBase->WfmBase::GetFileLVL().second == 5)))
               b_SubChName_out = m_ChFitParam->Branch("SubChName", &m_SubChName);
           b_MaxIndxFtChi2_out = m_ChFitParam->Branch("MaxIndxFtChi2", &m_MaxIndxFtChi2, "m_MaxIndxFtChi2/D");
           b_MinIndxFtChi2_out = m_ChFitParam->Branch("MinIndxFtChi2", &m_MinIndxFtChi2, "m_MinIndxFtChi2/D");
           b_MaxVoltFtChi2_out = m_ChFitParam->Branch("MaxVoltFtChi2", &m_MaxVoltFtChi2, "m_MaxVoltFtChi2/D");
           b_MinVoltFtChi2_out = m_ChFitParam->Branch("MinVoltFtChi2", &m_MinVoltFtChi2, "m_MinVoltFtChi2/D");
           b_StrIndxFtChi2_out = m_ChFitParam->Branch("StrIndxFtChi2", &m_StrIndxFtChi2, "m_StrIndxFtChi2/D");
           b_EndIndxFtChi2_out = m_ChFitParam->Branch("EndIndxFtChi2", &m_EndIndxFtChi2, "m_EndIndxFtChi2/D");
           b_NoiseFtChi2_out = m_ChFitParam->Branch("NoiseFtChi2", &m_NoiseFtChi2, "m_NoiseFtChi2/D");
           b_NoiseErrFtChi2_out = m_ChFitParam->Branch("NoiseErrFtChi2", &m_NoiseErrFtChi2, "m_NoiseErrFtChi2/D");
           b_PedestFtChi2_out = m_ChFitParam->Branch("PedestFtChi2", &m_PedestFtChi2, "m_PedestFtChi2/D");
           b_PedestErrFtChi2_out = m_ChFitParam->Branch("PedestErrFtChi", &m_PedestErrFtChi2, "m_PedestErrFtChi2/D");
           b_MaxTimeFtChi2_out = m_ChFitParam->Branch("MaxTimeFtChi2", &m_MaxTimeFtChi2, "m_MaxTimeFtChi2/D");
           b_MinTimeFtChi2_out = m_ChFitParam->Branch("MinTimeFtChi2", &m_MinTimeFtChi2, "m_MinTimeFtChi2/D");
           b_ChargeFtChi2_out = m_ChFitParam->Branch("ChargeFtChi2", &m_ChargeFtChi2, "m_ChargeFtChi2/D");
           b_RiseTimeFtChi2_out = m_ChFitParam->Branch("RiseTimeFtChi2", &m_RiseTimeFtChi2, "m_RiseTimeFtChi2/D");
           b_TriggTimeFtChi2_out = m_ChFitParam->Branch("TriggTimeFtChi2", &m_TriggTimeFtChi2, "m_TriggTimeFtChi2/D");
           b_DVDTMaxFtChi2_out = m_ChFitParam->Branch("DVDTMaxFtChi2", &m_DVDTMaxFtChi2, "m_DVDTMaxFtChi2/D");
           b_TriggToTFtChi2_out = m_ChFitParam->Branch("TriggToTFtChi2", &m_TriggToTFtChi2, "m_TriggToTFtChi2/D");
           b_SignalFFTFtChi2_out = m_ChFitParam->Branch("SignalFFTFtChi2", &m_SignalFFTFtChi2, "m_SignalFFTFtChi2/D");
           b_NoiseFFTFtChi2_out = m_ChFitParam->Branch("NoiseFFTFtChi2", &m_NoiseFFTFtChi2, "m_NoiseFFTFtChi2/D");
           // Fit parameters for each channel
           b_ChFitMaxVolt_out = m_ChFitParam->Branch("ChFitMaxVolt", &m_ChFitMxVlt.first, "m_ChFitMxVlt.first/D");
           b_ChFitMaxVoltErr_out = m_ChFitParam->Branch("ChFitMaxVoltErr", &m_ChFitMxVlt.second, "m_ChFitMxVlt.second/D");
           b_ChFitMinVolt_out = m_ChFitParam->Branch("ChFitMinVolt", &m_ChFitMinVlt.first, "m_ChFitMinVlt.first/D");
           b_ChFitMinVoltErr_out = m_ChFitParam->Branch("ChFitMinVoltErr", &m_ChFitMinVlt.second, "m_ChFitMinVlt.second/D");
           b_ChFitMxIndx_out = m_ChFitParam->Branch("ChFitMxIndx", &m_ChFitMxIndx.first, "m_ChFitMxIndx.first/I");
           b_ChFitMxIndxErr_out = m_ChFitParam->Branch("ChFitMxIndxErr", &m_ChFitMxIndx.second, "m_ChFitMxIndx.second/I");
           b_ChFitMinIndx_out = m_ChFitParam->Branch("ChFitMinIndx", &m_ChFitMinIndx.first, "m_ChFitMinIndx.first/I");
           b_ChFitMinIndxErr_out = m_ChFitParam->Branch("ChFitMinIndxErr", &m_ChFitMinIndx.second, "m_ChFitMinIndx.second/I");
           b_ChFitStrIndx_out = m_ChFitParam->Branch("ChFitStrIndx", &m_ChFitStrIndx.first, "m_ChFitStrIndx.first/I");
           b_ChFitStrIndxErr_out = m_ChFitParam->Branch("ChFitStrIndxErr", &m_ChFitStrIndx.second, "m_ChFitStrIndx.second/I");
           b_ChFitEndIndx_out = m_ChFitParam->Branch("ChFitEndIndx", &m_ChFitEndIndx.first, "m_ChFitEndIndx.first/I");
           b_ChFitEndIndxErr_out = m_ChFitParam->Branch("ChFitEndIndxErr", &m_ChFitEndIndx.second, "m_ChFitEndIndx.second/I");
           b_ChFitIsSig_out = m_ChFitParam->Branch("ChFitIsSig", &m_ChFitIsSig.first, "m_ChFitIsSig.first/D");
           b_ChFitIsSigErr_out = m_ChFitParam->Branch("ChFitIsSigErr", &m_ChFitIsSig.second, "m_ChFitIsSig.second/D");
           b_ChFitInWin_out = m_ChFitParam->Branch("ChFitInWin", &m_ChFitInWin.first, "m_ChFitInWin.first/D");
           b_ChFitInWinErr_out = m_ChFitParam->Branch("ChFitInWinErr", &m_ChFitInWin.second, "m_ChFitInWin.second/D");
           b_ChFitVStr_out = m_ChFitParam->Branch("ChFitVStr", &m_ChFitVStr.first, "m_ChFitVStr.first/D");
           b_ChFitVStrErr_out = m_ChFitParam->Branch("ChFitVStrErr", &m_ChFitVStr.second, "m_ChFitVStr.second/D");
           b_ChFitNoise_out = m_ChFitParam->Branch("ChFitNoise", &m_ChFitNoise.first, "m_ChFitNoise.first/D");
           b_ChFitNoiseUnc_out = m_ChFitParam->Branch("ChFitNoiseUnc", &m_ChFitNoise.second, "m_ChFitNoise.second/D");
           b_ChFitPedestal_out = m_ChFitParam->Branch("ChFitPedestal", &m_ChFitPedestal.first, "m_ChFitPedestal.first/D");
           b_ChFitPedestalUnc_out = m_ChFitParam->Branch("ChFitPedestalUnc", &m_ChFitPedestal.second, "m_ChFitPedestal.second/D");
           b_ChFitNoiseSlp_out = m_ChFitParam->Branch("ChFitNoiseSlp", &m_ChFitNoiseSlp.first, "m_ChFitNoiseSlp.first/D");
           b_ChFitNoiseSlpErr_out = m_ChFitParam->Branch("ChFitNoiseSlpErr", &m_ChFitNoiseSlp.second, "m_ChFitNoiseSlp.second/D");
           b_ChFitPedestalSlp_out = m_ChFitParam->Branch("ChFitPedestalSlp", &m_ChFitPedestalSlp.first, "m_ChFitPedestalSlp.first/D");
           b_ChFitPedestalSlpErr_out = m_ChFitParam->Branch("ChFitPedestalSlpErr", &m_ChFitPedestalSlp.second, "m_ChFitPedestalSlp.second/D");
           b_ChFitNoiseErr_out = m_ChFitParam->Branch("ChFitNoiseErr", &m_ChFitNoiseErr.first, "m_ChFitNoiseErr.first/D");
           b_ChFitNoiseErrUnc_out = m_ChFitParam->Branch("ChFitNoiseErrUnc", &m_ChFitNoiseErr.second, "m_ChFitNoiseErr.second/D");
           b_ChFitPedestalErr_out = m_ChFitParam->Branch("ChFitPedestalErr", &m_ChFitPedestalErr.first, "m_ChFitPedestalErr.first/D");
           b_ChFitPedestalErrUnc_out = m_ChFitParam->Branch("ChFitPedestalErrUnc", &m_ChFitPedestalErr.second, "m_ChFitPedestalErr.second/D");
           b_ChFitMaxTime_out = m_ChFitParam->Branch("ChFitMaxTime", &m_ChFitMaxTime.first, "m_ChFitMaxTime.first/D");
           b_ChFitMaxTimeErr_out = m_ChFitParam->Branch("ChFitMaxTimeErr", &m_ChFitMaxTime.second, "m_ChFitMaxTime.second/D");
           b_ChFitMinTime_out = m_ChFitParam->Branch("ChFitMinTime", &m_ChFitMinTime.first, "m_ChFitMinTime.first/D");
           b_ChFitMinTimeErr_out = m_ChFitParam->Branch("ChFitMinTimeErr", &m_ChFitMinTime.second, "m_ChFitMinTime.second/D");
           b_ChFitCharge_out = m_ChFitParam->Branch("ChFitCharge", &m_ChFitCharge.first, "m_ChFitCharge.first/D");
           b_ChFitChargeErr_out = m_ChFitParam->Branch("ChFitChargeErr", &m_ChFitCharge.second, "m_ChFitCharge.second/D");
           b_ChFitRiseT_out = m_ChFitParam->Branch("ChFitRiseT", &m_ChFitRiseT.first, "m_ChFitRiseT.first/D");
           b_ChFitRiseTErr_out = m_ChFitParam->Branch("ChFitRiseTErr", &m_ChFitRiseT.second, "m_ChFitRiseT.second/D");
           b_ChFitJitRiseSNR_out = m_ChFitParam->Branch("ChFitJitRiseSNR", &m_ChFitJitRiseSNR.first, "m_ChFitJitRiseSNR.first/D");
           b_ChFitJitRiseSNRErr_out = m_ChFitParam->Branch("ChFitJitRiseSNRErr", &m_ChFitJitRiseSNR.second, "m_ChFitJitRiseSNR.second/D");
           b_ChFitTriggTime_out = m_ChFitParam->Branch("ChFitTriggTime", &m_ChFitTriggTime.first, "m_ChFitTriggTime.first/D");
           b_ChFitTriggTimeErr_out = m_ChFitParam->Branch("ChFitTriggTimeErr", &m_ChFitTriggTime.second, "m_ChFitTriggTime.second/D");
           b_ChFitDVDTMax_out = m_ChFitParam->Branch("ChFitDVDTMax", &m_ChFitDVDTMax.first, "m_ChFitDVDTMax.first/D");
           b_ChFitDVDTMaxErr_out = m_ChFitParam->Branch("ChFitDVDTMaxErr", &m_ChFitDVDTMax.second, "/m_ChFitDVDTMax.secondD");
           b_ChFitTriggToT_out = m_ChFitParam->Branch("ChFitTriggToT", &m_ChFitTriggToT.first, "m_ChFitTriggToT.first/D");
           b_ChFitTriggToTErr_out = m_ChFitParam->Branch("ChFitTriggToTErr", &m_ChFitTriggToT.second, "m_ChFitTriggToT.second/D");
           b_ChFitComplete_out = m_ChFitParam->Branch("ChFitComplete", &m_ChFitComplete.first, "m_ChFitComplete.first/D");
           b_ChFitCompleteErr_out = m_ChFitParam->Branch("ChFitCompleteErr", &m_ChFitComplete.second, "m_ChFitComplete.second/D");
           b_ChFitSoNR_out = m_ChFitParam->Branch("ChFitSoNR", &m_ChFitSoNR.first, "m_ChFitSoNR.first/D");
           b_ChFitSoNRErr_out = m_ChFitParam->Branch("ChFitSoNRErr", &m_ChFitSoNR.second, "m_ChFitSoNR.second/D");
           // CFD and DV/DT vectors
           b_ChFitCFDTime_out = m_ChFitParam->Branch("ChFitCFDTime", &m_ChFitCFDTime);
           b_ChFitDVDTCFD_out = m_ChFitParam->Branch("ChFitDVDTCFD", &m_ChFitDVDTCFD);
           b_ChFitCFDToT_out = m_ChFitParam->Branch("ChFitCFDToT", &m_ChFitCFDToT);
           b_ChFitJitNdVdT_out = m_ChFitParam->Branch("ChFitJitNdVdT", &m_ChFitJitNdVdT);
           b_CFDTimeFtChi2_out = m_ChFitParam->Branch("CFDTimeFtChi2", &m_CFDTimeFtChi2);
           b_DVDTCFDFtChi2_out = m_ChFitParam->Branch("DVDTCFDFtChi2", &m_DVDTCFDFtChi2);
           b_CFDToTFtChi2_out = m_ChFitParam->Branch("CFDToTFtChi2", &m_CFDToTFtChi2);
           if (m_RunBase->WfmBase::GetDoFFT() || (newfile && ((m_RunBase->WfmBase::GetHasFFT().at(1)).second || (m_RunBase->WfmBase::GetHasFFT().at(2)).second)))
              {
               b_ChFitSignalFFT_out = m_ChFitParam->Branch("ChFitSignalFFT", &m_ChFitSignalFFT.first, "m_ChFitSignalFFT.first/D");
               b_ChFitSignalFFTErr_out = m_ChFitParam->Branch("ChFitSignalFFTErr", &m_ChFitSignalFFT.second, "m_ChFitSignalFFT.second/D");
               b_ChFitNoiseFFT_out = m_ChFitParam->Branch("ChFitNoiseFFT", &m_ChFitNoiseFFT.first, "m_ChFitNoiseFFT.first/D");
               b_ChFitNoiseFFT_out = m_ChFitParam->Branch("ChFitNoiseFFTErr", &m_ChFitNoiseFFT.second, "m_ChFitNoiseFFT.second/D");
              }
           if (m_RunBase->WfmBase::GetDoWaveShape() || (newfile && (m_RunBase->WfmBase::GetHasWaveShape()).second > 0))
              {
               b_ChShapeName_out = m_ChShapeParam->Branch("ChShapeName", &m_ChShapeName, "m_ChShapeName/I");
               // mot enougph to fo SbChannel, waveshape has to be on. The file level, replace with subchannel waveshape
               if ((SubChAna && m_RunBase->WfmBase::GetDoWaveShape()) || (newfile && (m_RunBase->WfmBase::GetHasWaveShape()).second > 1))
                   b_SubChShName_out = m_ChShapeParam->Branch("SubChShName", &m_SubChShName, "m_SubChShNam/I");
               b_ChShapePolarity_out = m_ChShapeParam->Branch("ChShapePolarity", &m_ChShapePolarity, "m_ChShapePolarity/I");
               b_ChShapeMaxIndx_out = m_ChShapeParam->Branch("ChShapeMaxIndx", &m_ChShapeMaxIndx, "m_ChShapeMaxIndx/I");
               b_ChShapeMinIndx_out = m_ChShapeParam->Branch("ChShapeMinIndx", &m_ChShapeMinIndx, "m_ChShapeMinIndx/I");
               b_ChShapeStrIndx_out = m_ChShapeParam->Branch("ChShapeStrIndx", &m_ChShapeStrIndx, "m_ChShapeStrIndx/I");
               b_ChShapeEndIndx_out = m_ChShapeParam->Branch("ChShapeEndIndx", &m_ChShapeEndIndx, "m_ChShapeEndIndx/I");
               b_ChShapeMax_out = m_ChShapeParam->Branch("ChShapeMax", &m_ChShapeMax, "m_ChShapeMax/D");
               b_ChShapeMin_out = m_ChShapeParam->Branch("ChShapeMin", &m_ChShapeMin, "m_ChShapeMin/D");
               b_ChShapeMaxTime_out = m_ChShapeParam->Branch("ChShapeMaxTime", &m_ChShapeMaxTime, "m_ChShapeMaxTime/D");
               b_ChShapeMinTime_out = m_ChShapeParam->Branch("ChShapeMinTime", &m_ChShapeMinTime, "m_ChShapeMinTime/D");
               b_ChShapeCharge_out = m_ChShapeParam->Branch("ChShapeCharge", &m_ChShapeCharge, "m_ChShapeCharge/D");
               b_ChShapeRiseTime_out = m_ChShapeParam->Branch("ChShapeRiseTime", &m_ChShapeRiseTime, "m_ChShapeRiseTime/D");
               b_ChShapeJitRiseSNR_out = m_ChShapeParam->Branch("ChShapeJitRiseSNR", &m_ChShapeJitRiseSNR, "m_ChShapeJitRiseSNR/D");
               b_ChShapeTriggTime_out = m_ChShapeParam->Branch("ChShapeTriggTime", &m_ChShapeTriggTime, "m_ChShapeTriggTime/D");
               b_ChShapeDVDTMax_out = m_ChShapeParam->Branch("ChShapeDVDTMax", &m_ChShapeDVDTMax, "m_ChShapeDVDTMax/D");
               b_ChShapeTriggToT_out = m_ChShapeParam->Branch("ChShapeTriggToT", &m_ChShapeTriggToT, "m_ChShapeTriggToT/D");
               b_ChShapeNoise_out = m_ChShapeParam->Branch("ChShapeNoise", &m_ChShapeNoise.first, "m_ChShapeNoise.first/D");
               b_ChShapeNoiseErr_out = m_ChShapeParam->Branch("ChShapeNoiseErr", &m_ChShapeNoise.second, "m_ChShapeNoise.second/D");
               b_ChShapeSoN_out = m_ChShapeParam->Branch("ChShapeSoN", &m_ChShapeSoN.first, "m_ChShapeSoN.first/D");
               b_ChShapeSoNUnc_out = m_ChShapeParam->Branch("ChShapeSoNUnc", &m_ChShapeSoN.second, "m_ChShapeSoN.second/D");
               // CFD and DV/DT Waveshape vectors
               b_ChShapeCFDTime_out = m_ChShapeParam->Branch("ChShapeCFDTime", &m_ChShapeCFDTime);
               b_ChShapeDVDTCFD_out = m_ChShapeParam->Branch("ChShapeDVDTCFD", &m_ChShapeDVDTCFD);
               b_ChShapeCFDToT_out = m_ChShapeParam->Branch("ChShapeCFDToT", &m_ChShapeCFDToT);
               b_ChShapeJitNdVdT_out = m_ChShapeParam->Branch("ChShapeJitNdVdT", &m_ChShapeJitNdVdT);
               b_ChMeanVolt_out = m_ChShapeParam->Branch("ChMeanVolt", &m_ChMeanVolt);
               if ((m_RunBase->WfmBase::GetDoFFT() && m_RunBase->WfmBase::GetDoWaveShape()) ||
                   (newfile && ((m_RunBase->WfmBase::GetHasFFT().at(3)).second || (m_RunBase->WfmBase::GetHasFFT().at(4)).second)))
                  {
                   b_ChShapeSignalFFT_out = m_ChShapeParam->Branch("ChShapeSignalFFT", &m_ChShapeSignalFFT, "m_ChShapeSignalFFT/D");
                   b_ChShapeNoiseFFT_out = m_ChShapeParam->Branch("ChShapeNoiseFFT", &m_ChShapeNoiseFFT, "m_ChShapeNoiseFFT/D");
                  }
              }
          }
       else {
             if (lv2FFT || lv3FFT)
                {
                 b_ChFitSignalFFT_out = m_ChFitParam->Branch("ChFitSignalFFT", &m_ChFitSignalFFT.first, "m_ChFitSignalFFT.first/D");
                 b_ChFitSignalFFTErr_out = m_ChFitParam->Branch("ChFitSignalFFTErr", &m_ChFitSignalFFT.second, "m_ChFitSignalFFT.second/D");
                 b_ChFitNoiseFFT_out = m_ChFitParam->Branch("ChFitNoiseFFT", &m_ChFitNoiseFFT.first, "m_ChFitNoiseFFT.first/D");
                 b_ChFitNoiseFFT_out = m_ChFitParam->Branch("ChFitNoiseFFTErr", &m_ChFitNoiseFFT.second, "m_ChFitNoiseFFT.second/D");
                }
             if (lv2Wsh || lv3Wsh)
                {
                 b_ChShapeName_out = m_ChShapeParam->Branch("ChShapeName", &m_ChShapeName, "m_ChShapeName/I");
                 // mot enougph to fo SbChannel, waveshape has to be on. The file level, replace with subchannel waveshape
                 if (SubChAna) b_SubChShName_out = m_ChShapeParam->Branch("SubChShName", &m_SubChShName, "m_SubChShNam/I");
                 b_ChShapePolarity_out = m_ChShapeParam->Branch("ChShapePolarity", &m_ChShapePolarity, "m_ChShapePolarity/I");
                 b_ChShapeMaxIndx_out = m_ChShapeParam->Branch("ChShapeMaxIndx", &m_ChShapeMaxIndx, "m_ChShapeMaxIndx/I");
                 b_ChShapeMinIndx_out = m_ChShapeParam->Branch("ChShapeMinIndx", &m_ChShapeMinIndx, "m_ChShapeMinIndx/I");
                 b_ChShapeStrIndx_out = m_ChShapeParam->Branch("ChShapeStrIndx", &m_ChShapeStrIndx, "m_ChShapeStrIndx/I");
                 b_ChShapeEndIndx_out = m_ChShapeParam->Branch("ChShapeEndIndx", &m_ChShapeEndIndx, "m_ChShapeEndIndx/I");
                 b_ChShapeMax_out = m_ChShapeParam->Branch("ChShapeMax", &m_ChShapeMax, "m_ChShapeMax/D");
                 b_ChShapeMin_out = m_ChShapeParam->Branch("ChShapeMin", &m_ChShapeMin, "m_ChShapeMin/D");
                 b_ChShapeMaxTime_out = m_ChShapeParam->Branch("ChShapeMaxTime", &m_ChShapeMaxTime, "m_ChShapeMaxTime/D");
                 b_ChShapeMinTime_out = m_ChShapeParam->Branch("ChShapeMinTime", &m_ChShapeMinTime, "m_ChShapeMinTime/D");
                 b_ChShapeCharge_out = m_ChShapeParam->Branch("ChShapeCharge", &m_ChShapeCharge, "m_ChShapeCharge/D");
                 b_ChShapeRiseTime_out = m_ChShapeParam->Branch("ChShapeRiseTime", &m_ChShapeRiseTime, "m_ChShapeRiseTime/D");
                 b_ChShapeJitRiseSNR_out = m_ChShapeParam->Branch("ChShapeJitRiseSNR", &m_ChShapeJitRiseSNR, "m_ChShapeJitRiseSNR/D");
                 b_ChShapeTriggTime_out = m_ChShapeParam->Branch("ChShapeTriggTime", &m_ChShapeTriggTime, "m_ChShapeTriggTime/D");
                 b_ChShapeDVDTMax_out = m_ChShapeParam->Branch("ChShapeDVDTMax", &m_ChShapeDVDTMax, "m_ChShapeDVDTMax/D");
                 b_ChShapeTriggToT_out = m_ChShapeParam->Branch("ChShapeTriggToT", &m_ChShapeTriggToT, "m_ChShapeTriggToT/D");
                 b_ChShapeNoise_out = m_ChShapeParam->Branch("ChShapeNoise", &m_ChShapeNoise.first, "m_ChShapeNoise.first/D");
                 b_ChShapeNoiseErr_out = m_ChShapeParam->Branch("ChShapeNoiseErr", &m_ChShapeNoise.second, "m_ChShapeNoise.second/D");
                 b_ChShapeSoN_out = m_ChShapeParam->Branch("ChShapeSoN", &m_ChShapeSoN.first, "m_ChShapeSoN.first/D");
                 b_ChShapeSoNUnc_out = m_ChShapeParam->Branch("ChShapeSoNUnc", &m_ChShapeSoN.second, "m_ChShapeSoN.second/D");
                 // CFD and DV/DT Waveshape vectors
                 b_ChShapeCFDTime_out = m_ChShapeParam->Branch("ChShapeCFDTime", &m_ChShapeCFDTime);
                 b_ChShapeDVDTCFD_out = m_ChShapeParam->Branch("ChShapeDVDTCFD", &m_ChShapeDVDTCFD);
                 b_ChShapeCFDToT_out = m_ChShapeParam->Branch("ChShapeCFDToT", &m_ChShapeCFDToT);
                 b_ChShapeJitNdVdT_out = m_ChShapeParam->Branch("ChShapeJitNdVdT", &m_ChShapeJitNdVdT);
                 b_ChMeanVolt_out = m_ChShapeParam->Branch("ChMeanVolt", &m_ChMeanVolt);
                 if (m_RunBase->WfmBase::GetDoFFT())
                    {
                     b_ChShapeSignalFFT_out = m_ChShapeParam->Branch("ChShapeSignalFFT", &m_ChShapeSignalFFT, "m_ChShapeSignalFFT/D");
                     b_ChShapeNoiseFFT_out = m_ChShapeParam->Branch("ChShapeNoiseFFT", &m_ChShapeNoiseFFT, "m_ChShapeNoiseFFT/D");
                    }
                }
             if (lv2Wsh_FFT || lv3Wsh_FFT)
                {
                 b_ChShapeSignalFFT_out = m_ChShapeParam->Branch("ChShapeSignalFFT", &m_ChShapeSignalFFT, "m_ChShapeSignalFFT/D");
                 b_ChShapeNoiseFFT_out = m_ChShapeParam->Branch("ChShapeNoiseFFT", &m_ChShapeNoiseFFT, "m_ChShapeNoiseFFT/D");
                }
            }
      }

   // Level 4 - 5 analysis requirements 			
   // Creating new file while file level max > 3 irrelevant of analysis stage
   // Creating new file, file level min < 4 while analysis stage > 3 
   // Not creating new file but file level min < 4 while analysis stage > 3
   if ((newfile && m_RunBase->WfmBase::GetFileLVL().second > 3) || (m_RunBase->WfmBase::GetFileLVL().first < 4 && m_RunBase->WfmBase::GetAnaStage() > 3))
      {
       b_TmDiffChi2_out = m_TimeResParam->Branch("TmDiffChi2", &m_TmDiffChi2);
       b_TmDiffFitMean_out = m_TimeResParam->Branch("TmDiffFitMean", &m_TmDiffFitMean);
       b_TmDiffFitErr_out = m_TimeResParam->Branch("TmDiffFitErr", &m_TmDiffFitErr);
       b_CFD_ChA_out = m_TimeResParam->Branch("CFD_ChA", &m_CFD_ChA);
       b_ChA_out = m_TimeResParam->Branch("ChA", &m_ChA);
       b_ChB_out = m_TimeResParam->Branch("ChB", &m_ChB);
      }

   return true;
}
// --------------------------------------------------------------------------------------------------------------
bool WfmRun::AddOutputHistos()
{

    // Intial test used later on in the funciton
    bool newfile = (m_cp == -99 || (m_cp != -99 && m_RunBase->WfmBase::GetInRootFiles().size() > 1));

   // Only propagate event number if creating new file
    if (newfile)
       {
        // If processing all events clone first histogram in the series and add other on top of it else re-create it
        if (m_RunBase->WfmBase::GetInstrument() != TestBeamBin1 || m_RunBase->WfmBase::GetInstrument() != LeCroyWRBin)
           {
            if ((m_RunBase->WfmBase::GetStartStopEvnt().first == 0 && m_RunBase->WfmBase::GetStartStopEvnt().second == 0) || m_RunBase->WfmBase::GetInstrument() == Sampic)
               {
                for (unsigned int a = 0; a < m_RunBase->WfmBase::GetInRootFiles().size(); a++)
                    {
                     if (m_trigDt_In.at(a) && m_trigFr_In.at(a)) 
                        { 
                         m_trigDt = static_cast<TH1F*>(m_trigDt_In.at(a)->Clone("trigDt"));
                         m_trigFr = static_cast<TH1F*>(m_trigFr_In.at(a)->Clone("trigFr"));
                         break;
                        }
                     else continue;
                    }
               }
            else {
                  m_trigDt = new TH1F("trigDt", "DT time", 6000, 0.0, 60.0);
                  m_trigFr = new TH1F("trigFr", "Trigger Frequency", 50000, 0.0, 5000);
                 }
           }
        // If Fit Histograms exist clone first in the series, else re-create it
        if (m_RunBase->WfmBase::GetFileLVL().second > 0)
           {
            unsigned int lvl1ft = 0;
            for (unsigned int a = 0; a < m_RunBase->WfmBase::GetInRootFiles().size(); a++)
                {
                 if (m_filelevels.at(a) > 0)
                    {
                     if (m_FitVarIn[0].at(lvl1ft) && m_FitVarIn[1].at(lvl1ft) && m_FitVarIn[2].at(lvl1ft) && m_FitVarIn[3].at(lvl1ft) && m_FitVarIn[4].at(lvl1ft))
                        {
                         m_FitMagVar = static_cast<TH2D*>(m_FitVarIn[0].at(lvl1ft)->Clone("VarFit_Mag"));
                         m_FitMagErrVar = static_cast<TH2D*>(m_FitVarIn[1].at(lvl1ft)->Clone("VarFit_MagErr"));
                         m_FitSigmaVar = static_cast<TH2D*>(m_FitVarIn[2].at(lvl1ft)->Clone("VarFit_Sigma"));
                         m_FitSigmaErrVar = static_cast<TH2D*>(m_FitVarIn[3].at(lvl1ft)->Clone("VarFit_SigmaErr"));
                         m_FitChi2Var = static_cast<TH2D*>(m_FitVarIn[4].at(lvl1ft)->Clone("VarFit_Chi2"));
                         break;
                        }
                     else continue;
                     lvl1ft++;
                    }
                }
           }
        else {
              m_FitMagVar = new TH2D("VarFit_Mag", "VarFit_Mag", 6, 0, 6, 7, 0, 7);
              m_FitMagErrVar = new TH2D(Form("VarFit_MagErr", "VarFit_MagErr", 6, 0, 6, 7, 0, 7);
              m_FitSigmaVar = new TH2D(Form("VarFit_Sigma", "VarFit_Sigma", 6, 0, 6, 7, 0, 7);
              m_FitSigmaErrVar = new TH2D("VarFit_SigmaErr", "VarFit_SigmaErr", 6, 0, 6, 7, 0, 7);
              m_FitChi2Var = new TH2D(Form("VarFit_Chi2", "VarFit_Chi2", 6, 0, 6, 7, 0, 7);
             }
       }
    else {
          // Copying on same file, recover esisting histogram to allow for modificaitons
          if (m_RunBase->WfmBase::GetInstrument() != TestBeamBin1 || m_RunBase->WfmBase::GetInstrument() != LeCroyWRBin)
             {
              m_trigDt = dynamic_cast<TH1F*>((m_RunBase->WfmBase::GetOutRootFile())->Get("trigDt"));
              m_trigFr = dynamic_cast<TH1F*>((m_RunBase->WfmBase::GetOutRootFile())->Get("trigFr"));
             }
          if (m_RunBase->WfmBase::GetFileLVL().first > 0)
             {
              m_FitMagVar = dynamic_cast<TH2D*>((m_RunBase->WfmBase::GetOutRootFile())->Get("VarFit_Mag"));
              m_FitMagErrVar = dynamic_cast<TH2D*>((m_RunBase->WfmBase::GetOutRootFile())->Get("VarFit_MagErr"));
              m_FitSigmaVar = dynamic_cast<TH2D*>((m_RunBase->WfmBase::GetOutRootFile())->Get("VarFit_Sigma"));
              m_FitSigmaErrVar = dynamic_cast<TH2D*>((m_RunBase->WfmBase::GetOutRootFile())->Get("VarFit_SigmaErr"));
              m_FitChi2Var = dynamic_cast<TH2D*>((m_RunBase->WfmBase::GetOutRootFile())->Get("VarFit_Chi2"));
             }
          else {
                m_FitMagVar = new TH2D("VarFit_Mag", "VarFit_Mag", 6, 0, 6, 7, 0, 7);
                m_FitMagErrVar = new TH2D(Form("VarFit_MagErr", "VarFit_MagErr", 6, 0, 6, 7, 0, 7);
                m_FitSigmaVar = new TH2D(Form("VarFit_Sigma", "VarFit_Sigma", 6, 0, 6, 7, 0, 7);
                m_FitSigmaErrVar = new TH2D("VarFit_SigmaErr", "VarFit_SigmaErr", 6, 0, 6, 7, 0, 7);
                m_FitChi2Var = new TH2D(Form("VarFit_Chi2", "VarFit_Chi2", 6, 0, 6, 7, 0, 7);
               }
         }

   // Second stage analysis object initialization 
   if (((m_RunBase->WfmBase::GetFileLVL()).second < 2 && m_RunBase->WfmBase::GetAnaStage() > 1) || m_RunBase->WfmBase::GetAnaStage() == 2 || m_RunBase->WfmBase::GetAnaStage() == 3)
      {
       // Fits Initialization
       h_MaxIndxFt.clear();
       h_MinIndxFt.clear();
       h_MaxVoltFt.clear();
       h_MinVoltFt.clear();
       h_StrIndxFt.clear();
       h_EndIndxFt.clear();
       h_NoiseFt.clear();
       h_NoiseErrFt.clear();
       h_PedestFt.clear();
       h_NPedestErrFt.clear();
       h_MaxTimeFt.clear();
       h_MinTimeFt.clear();
       h_ChargeFt.clear();
       h_RiseTimeFt.clear();
       h_TriggTimeFt.clear();
       h_DVDTMaxFt.clear();
       h_DVDTCFDFt.clear();
       h_CFDToTFt.clear();
       h_TriggToTFt.clear();
       h_SignalFFTFt.clear();
       h_NoiseFFTFt.clear();
       h_CFDTimeFt.clear();
 
       h_MaxIndxFt.resize(m_nchan, nullptr);
       h_MinIndxFt.resize(m_nchan, nullptr);
       h_MaxVoltFt.resize(m_nchan, nullptr);
       h_MinVoltFt.resize(m_nchan, nullptr);
       h_StrIndxFt.resize(m_nchan, nullptr);
       h_EndIndxFt.resize(m_nchan, nullptr);
       h_NoiseFt.resize(m_nchan, nullptr);
       h_NoiseErrFt.resize(m_nchan, nullptr);
       h_PedestFt.resize(m_nchan, nullptr);
       h_NPedestErrFt.resize(m_nchan, nullptr);
       h_MaxTimeFt.resize(m_nchan, nullptr);
       h_MinTimeFt.resize(m_nchan, nullptr);
       h_ChargeFt.resize(m_nchan, nullptr);
       h_RiseTimeFt.resize(m_nchan, nullptr);
       h_TriggTimeFt.resize(m_nchan, nullptr);
       h_DVDTMaxFt.resize(m_nchan, nullptr);
       h_TriggToTFt.resize(m_nchan, nullptr);
       h_SignalFFTFt.resize(m_nchan, nullptr);
       h_NoiseFFTFt.resize(m_nchan, nullptr);
       h_CFDTimeFt.resize(m_nchan);
       h_DVDTCFDFt.resize(m_nchan);
       h_CFDToTFt.resize(m_nchan);
       for (unsigned int s = 0; s < m_nchan; s++)
           {
            (h_CFDTimeFt.at(s)).clear();
            (h_CFDTimeFt.at(s)).resize(19, nullptr);
            (h_DVDTCFDFt.at(s)).clear();
            (h_DVDTCFDFt.at(s)).resize(19, nullptr);
            (h_CFDToTFt.at(s)).clear();
            (h_CFDToTFt.at(s)).resize(19, nullptr);
           }
      }

   // Third stage analysis histogram intialization
   if (((m_RunBase->WfmBase::GetFileLVL()).second < 3 && m_RunBase->WfmBase::GetAnaStage() > 2) || m_RunBase->WfmBase::GetAnaStage() == 3)
      {

      }

   return true;
}
// --------------------------------------------------------------------------------------------------------------
bool WfmRun::WriteBranches(unsigned int fileindx)
{
    // Only if I analyzed, not if I am coyinhg or doing second level!!!!
    if ((m_RunBase->WfmBase::GetAnaStage() == 1 && m_filelevels.at(fileindx) < 2) || (m_filelevels.at(fileindx) == 0 && &&m_RunBase->WfmBase::GetAnaStage() > 1))
       {
        m_FitDiffChi2 = m_FitDiffChi2/(m_filechannel.at(fileindx)).size();
        m_FitDiffMag = m_FitDiffMag/(m_filechannel.at(fileindx).size();
        m_FitDiffMagErr = m_FitDiffMagErr/(m_filechannel.at(fileindx)).size();
        m_FitDiffSigma = m_FitDiffSigma/(m_filechannel.at(fileindx)).size();
        m_FitDiffSigmaErr = m_FitDiffSigmaErr/(m_filechannel.at(fileindx)).size();
        m_FitDiffTime = m_FitDiffTime/(m_filechannel.at(fileindx)).size();
       }

   // If combining to one of the input files and current file not one where merged output would be or if creating new file
   if (m_RunBase->WfmBase::GetInRootFile(fileindx) != m_RunBase->WfmBase::GetOutRootFile()) m_OutTree->Fill();
   else {
         // Current file one on which combined output will go or not creating new file
         if ((m_filelevels.at(fileindx) < 2 && m_RunBase->WfmBase::GetAnaStage() == 1) || (m_filelevels.at(fileindx) < 1 && &&m_RunBase->WfmBase::GetAnaStage() > 1))
            {
             b_FitDiffChi2_out->Fill();
             b_FitDiffMag_out->Fill();
             b_FitDiffMagErr_out->Fill();
             b_FitDiffSigma_out->Fill();
             b_FitDiffSigmaErr_out->Fill();
             b_FitDiffTime_out->Fill();
             for (unsigned int ich = 0; ich < m_nchan; ich++)
                 {
                  b_EvPol_out.at(ich)->Fill();
                  b_EvCharge_out.at(ich)->Fill();
                  b_EvIsSignal_out.at(ich)->Fill();
                  b_EvInWindow_out.at(ich)->Fill();
                  b_EvVoltStr_out.at(ich)->Fill();
                  b_EvMaxindx_out.at(ich)->Fill();
                  b_EvMaxTime_out.at(ich)->Fill();
                  b_EvMaxVolt_out.at(ich)->Fill();
                  b_EvMinIndx_out.at(ich)->Fill();
                  b_EvMinTime_out.at(ich)->Fill();
                  b_EvMinVolt_out.at(ich)->Fill();
                  b_EvStrIndx_out.at(ich)->Fill();
                  b_EvEndIndx_out.at(ich)->Fill();
                  b_EvNoise_out.at(ich)->Fill();
                  b_EvPedest_out.at(ich)->Fill();
                  b_EvNoiseErr_out.at(ich)->Fill();
                  b_EvPedestErr_out.at(ich)->Fill();
                  b_EvRiseT_out.at(ich)->Fill();
                  b_EvJitNdVdT_out.at(ich)->Fill();
                  b_EvJitRiseSNR_out.at(ich)->Fill();
                  b_EvTriggTime_out.at(ich)->Fill();
                  b_EvDVDTMax_out.at(ich)->Fill();
                  b_EvDVDTCFD_out.at(ich)->Fill();
                  b_EvCFDToT_out.at(ich)->Fill();
                  b_EvTriggToT_out.at(ich)->Fill();
                  if (m_RunBase->WfmBase::GetDoFFT())
                     {
                      b_EvSignalFFT_out.at(ich)->Fill();
                      b_EvNoiseFFT_out.at(ich)->Fill();
                     }
                  b_EvComplete_out.at(ich)->Fill();
                  b_EvCFDTime_out.at(ich)->Fill();
                  b_EvVAdjCFD_out.at(ich)->Fill();
                  if (m_RunBase->WfmBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Written all brancehs: " << ich << "/" << m_nchan << "!" << std::endl;
                 }
            }
        }
   return true;
}