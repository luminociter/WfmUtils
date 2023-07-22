/*
*                 WFMRun1Stage.cxx
*
*
*      Author: Gkougkousis Evangelos - Leonidas
*              egkougko@cern.ch - 2023
*                    CERN-GENEVA
*/

#include "../LGADUtils/LGADRun.h"

bool LGADRun::InitNtuple(std::vector<TTree*>& inputrees, std::vector<TFile*>& inputfiles)
{
   // For each file figure out the number of channels and the instrument
   m_filechannel.clear();
   m_filelvl.clear();
   m_ChTrees.clear();
   m_ShapeTrees.clear();
   m_TimeResTrees.clear();
   m_FFTShape.clear();
   m_fileplanes.clear();
   m_fileRefs.clear();
   m_fileDUTs.clear();

   std::vector<AqInstrument> fileInstr;
   std::vector<unsigned int> filetrack;
   filetrack.clear(); fileInstr.clear();
   for (unsigned int a = 0; a < inputrees.size(); a++)
       {
        // Determine present channels in file
       m_filechannel.push_back(std::vector<unsigned int>());
        for (unsigned int ich = 1; ich < 65; ich++)
            {
             if (inputrees.at(a)->GetBranch(Form("t%02u", ich))) (m_filechannel.back()).push_back(ich);
            }
        // Define instrument
        if (inputrees.at(a)->GetBranch(Form("ordrt%02u", (m_filechannel.back()).at(0)))) fileInstr.push_back(Sampic);
        else if (inputrees.at(a)->GetBranch(Form("triggTime%02u", (m_filechannel.back()).at(0)))) fileInstr.push_back(TektronixScope);
        else if (inputrees.at(a)->GetBranch(Form("Scope%02u", (m_filechannel.back()).at(0))))
                {
                 if (inputrees.at(a)->GetBranch("trigtime")) fileInstr.push_back(TestBeamBin2);
                 else fileInstr.push_back(TestBeamBin1);
                }
        else if (inputrees.at(a)->GetBranch("trigtime") && !inputrees.at(a)->GetBranch(Form("Scope%02u", (m_filechannel.back()).at(0))))
            fileInstr.push_back(LabTXT);
        else fileInstr.push_back(Unasigned);
        if (m_RunBase->LGADBase::GetVerbosity() >= 1) std::cout << __FUNCTION__ << " INFO: Aquisition instrument for file " << inputfiles.at(a)->GetName()
                                      << " identified as " << fileInstr.back() << "." << std::endl;
        // Define file Analysis level
        if ((TTree*)(inputfiles.at(a)->Get("TimeResParam")))
           {
            m_TimeResTrees.push_back(NULL);
            inputfiles.at(a)->GetObject("TimeResParam", m_TimeResTrees.back());
            m_filelvl.push_back(3);
            inputfiles.at(a)->GetObject("ChFitParam", m_ChTrees.back());
            if (m_ChTrees.back()->GetBranch("ChFitSignalFFT")) m_FFTShape.push_back(20);
            else m_FFTShape.push_back(10);
            if ((TTree*)(inputfiles.at(a)->Get("ChShapeParam"))) 
               {
                m_ShapeTrees.push_back(NULL);
                inputfiles.at(a)->GetObject("ChShapeParam", m_ShapeTrees.back());
                m_FFTShape.back() += 2;
               }
            else m_FFTShape.back() += 1;
           }
        else if ((TTree*)(inputfiles.at(a)->Get("ChFitParam")))
                {
                 m_ChTrees.push_back(NULL);
                 inputfiles.at(a)->GetObject("ChFitParam", m_ChTrees.back());
                 m_filelvl.push_back(2);
                 if (m_ChTrees.back()->GetBranch("ChFitSignalFFT")) m_FFTShape.push_back(20);
                 else m_FFTShape.push_back(10);
                 if ((TTree*)(inputfiles.at(a)->Get("ChShapeParam")))
                    {
                     m_ShapeTrees.push_back(NULL);
                     inputfiles.at(a)->GetObject("ChShapeParam", m_ShapeTrees.back());
                     m_FFTShape.back() += 2;
                    }
                 else m_FFTShape.back() += 1;
                }
        else if (inputrees.at(a)->GetBranch(Form("Polarity%02u", (m_filechannel.back()).at(0)))) 
                {
                 m_filelvl.push_back(1);
                 if (inputrees.at(a)->GetBranch(Form("SignalFFT%02u", (m_filechannel.back()).at(0)))) m_FFTShape.push_back(21);
                 else m_FFTShape.push_back(11);
                }
        else m_filelvl.push_back(0);
        if (m_RunBase->LGADBase::GetVerbosity() >= 1) std::cout << __FUNCTION__ << " INFO: Analysis level of file " << inputfiles.at(a)->GetName() << " :" << m_filelvl.back() << std::endl;
        // Presence of tracking information and trakcing package    
        if (inputrees.at(a)->GetBranch("trigger_id") && inputrees.at(a)->GetBranch("timestampTr")) filetrack.push_back(1);
        else if (inputrees.at(a)->GetBranch("HasAsso"))
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
                 filetrack.push_back(2);
                }
        else if (inputrees.at(a)->GetBranch("blabalblabla")) filetrack.push_back(3);
        else filetrack.push_back(4);
       }

   // Fix the number of channels for reserving the correct amount of pointers
   m_RunBase->LGADBase::FixBaseQuant(m_filechannel, m_channels, m_nchan);
   if (m_RunBase->LGADBase::GetVerbosity() >= 1) std::cout << __FUNCTION__ << " INFO: Total number of channels found " << m_nchan << " in " << inputrees.size() << " files." << std::endl;  // Debug

   // Verify that a single tracking frameworks and same qaquisition isntrument is used for all files
   for (unsigned int a = 1; a < inputrees.size(); a++)
       {
        if (filetrack.at(a) != filetrack.at(a-1)) 
           { 
            std::cout << __FUNCTION__ << " ERROR: Trying to mix different trakcing frameworks on input files, aborting...." << std::endl;
            return true;
           }
        if (fileInstr.at(a) != fileInstr.at(a - 1))
           { 
            std::cout << __FUNCTION__ << " ERROR: Trying to mix data from different measuring instruments, aborting...." << std::endl;
            return true;
           }
       }
   if (filetrack.size() > 0)
      { 
       m_RunBase->LGADBase::SetHasTrck(true);
       if (filetrack.at(0) == 1) m_RunBase->LGADBase::SetTrackPackage("EUTelescope ");
       else if (filetrack.at(0) == 2) m_RunBase->LGADBase::SetTrackPackage("Cory");
       else if (filetrack.at(0) == 3) m_RunBase->LGADBase::SetTrackPackage("TBA");
       else if (filetrack.at(0) == 4) m_RunBase->LGADBase::SetTrackPackage("Udef");
      }
   m_RunBase->LGADBase::SetInstrument(fileInstr.at(0));

   // Fix the file level for initializing the right pointers
   unsigned int lvl = *std::max_element(m_filelvl.begin(), m_filelvl.end());
   m_RunBase->LGADBase::SetFileLVL(lvl);
   // Establish if waveshape and FFT pointers are to be reserved
   bool fft = false;
   bool shape = false;
   for (unsigned int p = 0; p < m_filelvl.size(); p++) 
       { 
        if (m_FFTShape.at(p) == 21 || m_FFTShape.at(p) == 22) fft = true;
        if (m_FFTShape.at(p) == 22 || m_FFTShape.at(p) == 12) shape = true;
       }
   m_RunBase->LGADBase::SetHasFFT(fft);
   m_RunBase->LGADBase::SetHasWaveShape(shape);

   // Initialize zero level Pointers and vectors
   b_t.resize(m_nchan, nullptr);
   b_w.resize(m_nchan, nullptr);
   b_ordrt.resize(m_nchan, nullptr);
   b_physt.resize(m_nchan, nullptr);
   b_vScale.resize(m_nchan, nullptr);
   b_nPoints.resize(m_nchan, nullptr);
   b_SnRate.resize(m_nchan, nullptr);
   b_Scope.resize(m_nchan, nullptr);
   b_triggTime.resize(m_nchan, nullptr);
   m_RunBase->LGADBase::SetVectorSize(m_nchan);

   // initialize Cory tracking pointers size
   if (!m_RunBase->LGADBase::GetTrackComb())
      {
       if (m_RunBase->LGADBase::GetTrackPackage() == "Cory")
          {
           m_ntp_planes.clear();
           m_ntp_Ref.clear();
           m_ntp_DUT.clear();
           m_nRecDUTs = 0;
           m_nRecRefs = 0;
           m_nRecPlanes = 0;
           m_RunBase->LGADBase::FixBaseQuant(m_fileDUTs, m_ntp_DUT, m_nRecDUTs);
           m_RunBase->LGADBase::FixBaseQuant(m_fileDUTs, m_ntp_Ref, m_nRecRefs);
           m_RunBase->LGADBase::FixBaseQuant(m_fileplanes, m_ntp_planes, m_nRecPlanes);

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
           m_RunBase->LGADBase::InitCoryPointers(m_nRecRefs, m_nRecDUTs, m_nRecPlanes);
          }
       if (m_RunBase->LGADBase::GetTrackPackage() == "EUTelescope") m_RunBase->LGADBase::InitEUTelPointers();
       if (m_RunBase->LGADBase::GetTrackPackage() == "TBA") m_RunBase->LGADBase::InitTBAPointers();
      }

   // Intialize 1st level pointers and vectors
   if (lvl > 0)
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

       m_EvPol.resize(m_nchan, 0);
       m_EvIsSignal.resize(m_nchan, 0);
       m_EvInWindow.resize(m_nchan, 0);
       m_EvVoltStr.resize(m_nchan, 0);
       m_EvMaxindx.resize(m_nchan, 0);
       m_EvStrIndx.resize(m_nchan, 0);
       m_EvEndIndx.resize(m_nchan, 0);
       m_EvMinIndx.resize(m_nchan, 0);
       m_EvComplete.resize(m_nchan, 0);
       m_EvCharge.resize(m_nchan, 0.0);
       m_EvMaxTime.resize(m_nchan, 0.0);
       m_EvMaxVolt.resize(m_nchan, 0.0);
       m_EvMinTime.resize(m_nchan, 0.0);
       m_EvMinVolt.resize(m_nchan, 0.0);
       m_EvNoise.resize(m_nchan, 0.0);
       m_EvPedest.resize(m_nchan, 0.0);
       m_EvNoiseErr.resize(m_nchan, 0.0);
       m_EvPedestErr.resize(m_nchan, 0.0);
       m_EvRiseT.resize(m_nchan, 0.0);
       m_EvJitRiseSNR.resize(m_nchan, 0.0);
       m_EvTriggTime.resize(m_nchan, 0.0);
       m_EvDVDTMax.resize(m_nchan, 0.0);
       m_EvTriggToT.resize(m_nchan, 0.0);

       // Initalization of vector pointrers
       m_EvVAdjCFD.clear();
       m_EvCFDTime.clear();
       m_EvJitNdVdT.clear();
       m_EvDVDTCFD.clear();
       m_EvCFDToT.clear();
       m_EvVAdjCFD.resize(m_nchan, nullptr);
       m_EvCFDTime.resize(m_nchan, nullptr);
       m_EvJitNdVdT.resize(m_nchan, nullptr);
       m_EvDVDTCFD.resize(m_nchan, nullptr);
       m_EvCFDToT.resize(m_nchan, nullptr);

       if (fft)
          {
           b_EvSignalFFT.resize(m_nchan, nullptr);
           b_EvNoiseFFT.resize(m_nchan, nullptr);
           m_EvSignalFFT.resize(m_nchan, 0.0);
           m_EvNoiseFFT.resize(m_nchan, 0.0);
          }

       // Intialize 2nd level pointers and vectors
       if (lvl > 1)
         {
          m_CFDTimeFtChi2 = 0;
          m_DVDTCFDFtChi2 = 0;
          m_CFDToTFtChi2 = 0;
          m_ChFitCFDTime = 0;
          m_ChFitDVDTCFD = 0;
          m_ChFitCFDToT = 0;
          m_ChFitJitNdVdT = 0;
          if (shape)
             {
              m_ChShapeCFDTime = 0;
              m_ChShapeDVDTCFD = 0;
              m_ChShapeCFDToT = 0;
              m_ChShapeJitNdVdT = 0;
              m_ChMeanVolt = 0;
             }
          // Intialize 3rd level pointers and vectors
          if (lvl > 2)
             {
              m_TmDiffChi2 = 0;
              m_TmDiffFitMean = 0;
              m_TmDiffFitErr = 0;
              m_CFD_ChA = 0;
             }
         }
      }

   return true;
}
// --------------------------------------------------------------------------------------------------------------
// Function to be called every time a new Ntuple is oppened to corectly set the branch adresses for the devices present in the dataset.
bool LGADRun::SetBranchAccess(TTree* inputree, unsigned int FFTShape, unsigned int filelvl, std::vector<unsigned int> filechannel,
                              std::vector<unsigned int> fileplanes, std::vector<unsigned int> fileRefs, std::vector<unsigned int> fileDUTs, 
                              TTree* ChTree, TTree* ShapeTree, TTree* TimeResTree)
{
   inputree->SetBranchAddress("EvnNo", &EvnNo, &b_EvnNo);
   for (unsigned int ich = 0; ich < filechannel.size(); ich++)
       {
        unsigned int chindx = 999;
        for (unsigned int s = 0; s < m_channels.size(); s++) if (m_channels.at(s) == filechannel.at(ich)) { chindx = s; break; }
        if (chindx == 999)
           {
            std::cout << __FUNCTION__ << " ERROR: Could not match channel " << filechannel.at(ich) << " of this file, skipping...." << std::endl;
            return false;
           }
        inputree->SetBranchAddress(Form("t%02u", filechannel.at(ich)), &m_t.at(chindx), &b_t.at(chindx));
        inputree->SetBranchAddress(Form("w%02u", filechannel.at(ich)), &m_w.at(chindx), &b_w.at(chindx));
        if (m_RunBase->LGADBase::GetInstrument() == Sampic)
           {
            inputree->SetBranchAddress(Form("ordrt%02u", filechannel.at(ich)), &m_ordrt.at(chindx), &b_ordrt.at(chindx));
            inputree->SetBranchAddress(Form("physt%02u", filechannel.at(ich)), &m_physt.at(chindx), &b_physt.at(chindx));
           }
        if (m_RunBase->LGADBase::GetInstrument() == TestBeamBin1 || m_RunBase->LGADBase::GetInstrument() == TestBeamBin2 || m_RunBase->LGADBase::GetInstrument() == LabTXT || 
            m_RunBase->LGADBase::GetInstrument() == TektronixScope || m_RunBase->LGADBase::GetInstrument() == LeCroyWRBin)
           {
            inputree->SetBranchAddress(Form("nPoints%02u", filechannel.at(ich)), &m_npoints.at(chindx), &b_nPoints.at(chindx));
            inputree->SetBranchAddress(Form("SnRate%02u", filechannel.at(ich)), &m_srate.at(chindx), &b_SnRate.at(chindx));
            inputree->SetBranchAddress(Form("vScale%02u", filechannel.at(ich)), &m_scale.at(chindx), &b_vScale.at(chindx));
            if (m_RunBase->LGADBase::GetInstrument() == TektronixScope) inputree->SetBranchAddress(Form("triggTime%02u", filechannel.at(ich)), &m_triggTime.at(chindx), &b_triggTime.at(chindx));
            if (m_RunBase->LGADBase::GetInstrument() == TestBeamBin1 || m_RunBase->LGADBase::GetInstrument() == TestBeamBin2 || m_RunBase->LGADBase::GetInstrument() == LeCroyWRBin)
                inputree->SetBranchAddress(Form("Scope%02u", filechannel.at(ich)), &m_scope.at(chindx), &b_Scope.at(chindx));
           }
        // First stage analysis branches if they exist
        if (filelvl > 0)
           {
            inputree->SetBranchAddress(Form("Polarity%02u", filechannel.at(ich)), &m_EvPol.at(chindx), &b_EvPol.at(chindx));
            inputree->SetBranchAddress(Form("Charge%02u", filechannel.at(ich)), &m_EvCharge.at(chindx), &b_EvCharge.at(chindx));
            inputree->SetBranchAddress(Form("IsSignal%02u", filechannel.at(ich)), &m_EvIsSignal.at(chindx), &b_EvIsSignal.at(chindx));
            inputree->SetBranchAddress(Form("InWindow%02u", filechannel.at(ich)), &m_EvInWindow.at(chindx), &b_EvInWindow.at(chindx));
            inputree->SetBranchAddress(Form("VoltStr%02u", filechannel.at(ich)), &m_EvVoltStr.at(chindx), &b_EvVoltStr.at(chindx));
            inputree->SetBranchAddress(Form("MaxIndx%02u", filechannel.at(ich)), &m_EvMaxindx.at(chindx), &b_EvMaxindx.at(chindx));
            inputree->SetBranchAddress(Form("MaxTime%02u", filechannel.at(ich)), &m_EvMaxTime.at(chindx), &b_EvMaxTime.at(chindx));
            inputree->SetBranchAddress(Form("MaxVolt%02u", filechannel.at(ich)), &m_EvMaxVolt.at(chindx), &b_EvMaxVolt.at(chindx));
            inputree->SetBranchAddress(Form("MinIndx%02u", filechannel.at(ich)), &m_EvMinIndx.at(chindx), &b_EvMinIndx.at(chindx));
            inputree->SetBranchAddress(Form("MinTime%02u", filechannel.at(ich)), &m_EvMinTime.at(chindx), &b_EvMinTime.at(chindx));
            inputree->SetBranchAddress(Form("MinVolt%02u", filechannel.at(ich)), &m_EvMinVolt.at(chindx), &b_EvMinVolt.at(chindx));
            inputree->SetBranchAddress(Form("Strindx%02u", filechannel.at(ich)), &m_EvStrIndx.at(chindx), &b_EvStrIndx.at(chindx));
            inputree->SetBranchAddress(Form("EndIndx%02u", filechannel.at(ich)), &m_EvEndIndx.at(chindx), &b_EvEndIndx.at(chindx));
            inputree->SetBranchAddress(Form("Noise%02u", filechannel.at(ich)), &m_EvNoise.at(chindx), &b_EvNoise.at(chindx));
            inputree->SetBranchAddress(Form("Pedest%02u", filechannel.at(ich)), &m_EvPedest.at(chindx), &b_EvPedest.at(chindx));
            inputree->SetBranchAddress(Form("NoiseErr%02u", filechannel.at(ich)), &m_EvNoiseErr.at(chindx), &b_EvNoiseErr.at(chindx));
            inputree->SetBranchAddress(Form("PedestErr%02u", filechannel.at(ich)), &m_EvPedestErr.at(chindx), &b_EvPedestErr.at(chindx));
            inputree->SetBranchAddress(Form("RiseT%02u", filechannel.at(ich)), &m_EvRiseT.at(chindx), &b_EvRiseT.at(chindx));
            inputree->SetBranchAddress(Form("JitRiseSNR%02u", filechannel.at(ich)), &m_EvJitRiseSNR.at(chindx), &b_EvJitRiseSNR.at(chindx));
            inputree->SetBranchAddress(Form("TriggTime%02u", filechannel.at(ich)), &m_EvTriggTime.at(chindx), &b_EvTriggTime.at(chindx));
            inputree->SetBranchAddress(Form("DVDTMax%02u", filechannel.at(ich)), &m_EvDVDTMax.at(chindx), &b_EvDVDTMax.at(chindx));
            inputree->SetBranchAddress(Form("TriggToT%02u", filechannel.at(ich)), &m_EvTriggToT.at(chindx), &b_EvTriggToT.at(chindx));
            inputree->SetBranchAddress(Form("Complete%02u", filechannel.at(ich)), &m_EvComplete.at(chindx), &b_EvComplete.at(chindx));
            if (FFTShape == 21 || FFTShape == 22)
               {
                inputree->SetBranchAddress(Form("SignalFFT%02u", filechannel.at(ich)), &m_EvSignalFFT.at(chindx), &b_EvSignalFFT.at(chindx));
                inputree->SetBranchAddress(Form("NoiseFFT%02u", filechannel.at(ich)), &m_EvNoiseFFT.at(chindx), &b_EvNoiseFFT.at(chindx));
               }
            inputree->SetBranchAddress(Form("VAdjCFD%02u", filechannel.at(ich)), &m_EvVAdjCFD.at(chindx), &b_EvVAdjCFD.at(chindx));
            inputree->SetBranchAddress(Form("CFDTime%02u", filechannel.at(ich)), &m_EvCFDTime.at(chindx), &b_EvCFDTime.at(chindx));
            inputree->SetBranchAddress(Form("DVDTCFD%02u", filechannel.at(ich)), &m_EvDVDTCFD.at(chindx), &b_EvDVDTCFD.at(chindx));
            inputree->SetBranchAddress(Form("CFDToT%02u", filechannel.at(ich)), &m_EvCFDToT.at(chindx), &b_EvCFDToT.at(chindx));
            inputree->SetBranchAddress(Form("JitNdVdT%02u", filechannel.at(ich)), &m_EvJitNdVdT.at(chindx), &b_EvJitNdVdT.at(chindx));
           }
       }
   if (m_RunBase->LGADBase::GetInstrument() == Sampic)
      {
       inputree->SetBranchAddress("nPoints", &m_npoints.at(0), &b_nPoints.at(0));
       inputree->SetBranchAddress("SnRate", &m_srate.at(0), &b_SnRate.at(0));
      }
   if (m_RunBase->LGADBase::GetInstrument() == LabTXT || m_RunBase->LGADBase::GetInstrument() == TestBeamBin2) inputree->SetBranchAddress("trigtime", &m_trigtime, &b_trigtime);

   // EU Telescope Branches
   if (m_RunBase->LGADBase::GetTrackPackage() == "EUTelescope")
      {
       inputree->SetBranchAddress("xTr", &m_xTr, &b_xTr);
       inputree->SetBranchAddress("yTr", &m_yTr, &b_yTr);
       inputree->SetBranchAddress("omegaTr", &m_omegaTr, &b_omegaTr);
       inputree->SetBranchAddress("phiTr", &m_phiTr, &b_phiTr);
       inputree->SetBranchAddress("kinkxTr", &m_kinkxTr, &b_kinkxTr);
       inputree->SetBranchAddress("kinkyTr", &m_kinkyTr, &b_kinkyTr);
       inputree->SetBranchAddress("chi2Tr", &m_chi2Tr, &b_chi2Tr);
       inputree->SetBranchAddress("plane_id", &m_plane_id, &b_plane_id);
       inputree->SetBranchAddress("track_id", &m_track_id, &b_track_id);
       inputree->SetBranchAddress("ndofTr", &m_ndofTr, &b_ndofTr);
       inputree->SetBranchAddress("trigger_id", &m_trigger_id, &b_trigger_id);
       inputree->SetBranchAddress("timestampTr", &m_timestampTr, &b_timestampTr);
       inputree->SetBranchAddress("xHit", &m_xHit, &b_xHit);
       inputree->SetBranchAddress("yHit", &m_yHit, &b_yHit);
       inputree->SetBranchAddress("zHit", &m_zHit, &b_zHit);
       inputree->SetBranchAddress("ID", &m_ID, &b_ID);
      }

   // Cory Branches
   if (m_RunBase->LGADBase::GetTrackPackage() == "Cory")
      {
       inputree->SetBranchAddress("ClusterSizeX_ROI", &m_ClusterSizeX_ROI, &b_ClusterSizeX_ROI);
       inputree->SetBranchAddress("ClusterSizeY_ROI", &m_ClusterSizeY_ROI, &b_ClusterSizeY_ROI);
       inputree->SetBranchAddress("PixelX_ROI", &m_PixelX_ROI, &b_PixelX_ROI);
       inputree->SetBranchAddress("PixelY_ROI", &m_PixelY_ROI, &b_PixelY_ROI);
       inputree->SetBranchAddress("ClusterNumPixels_ROI", &m_ClusterNumPixels_ROI, &b_ClusterNumPixels_ROI);
       inputree->SetBranchAddress("InterceptX_ROI", &m_InterceptX_ROI, &b_InterceptX_ROI);
       inputree->SetBranchAddress("InterceptY_ROI", &m_InterceptY_ROI, &b_InterceptY_ROI);
       inputree->SetBranchAddress("InterceptResX_ROI", &m_InterceptResX_ROI, &b_InterceptResX_ROI);
       inputree->SetBranchAddress("InterceptResY_ROI", &m_InterceptResY_ROI, &b_InterceptResY_ROI);
       inputree->SetBranchAddress("HasAsso", &m_HasAsso, &b_HasAsso);
       for (unsigned int ref = 0; ref < fileRefs.size(); ref++)
           {
            unsigned int refindx = 999;
            for (unsigned int s = 0; s < m_ntp_Ref.size(); s++) { if (m_ntp_Ref.at(s) == fileRefs.at(ref)) { refindx = s; break; } }
            if (refindx == 999)
               {
                std::cout << __FUNCTION__ << " ERROR: Could not match refernce plane " << fileRefs.at(ref) << " for this run, skipping...." << std::endl;
                return false;
               }
            inputree->SetBranchAddress(Form("InterceptX_Ref%u", fileRefs.at(ref)), &m_InterceptX_Ref.at(refindx), &b_InterceptX_Ref.at(refindx));
            inputree->SetBranchAddress(Form("InterceptY_Ref%u", fileRefs.at(ref)), &m_InterceptY_Ref.at(refindx), &b_InterceptY_Ref.at(refindx));
            inputree->SetBranchAddress(Form("InterceptResX_Ref%u", fileRefs.at(ref)), &m_InterceptResX_Ref.at(refindx), &b_InterceptResX_Ref.at(refindx));
            inputree->SetBranchAddress(Form("InterceptResY_Ref%u", fileRefs.at(ref)), &m_InterceptResY_Ref.at(refindx), &b_InterceptResY_Ref.at(refindx));
           }
       for (unsigned int dut = 0; dut < fileDUTs.size(); dut++)
           {
            unsigned int dutindx = 999;
            for (unsigned int s = 0; s < m_ntp_DUT.size(); s++) { if (m_ntp_DUT.at(s) == fileDUTs.at(dut)) { dutindx = s; break; } }
            if (dutindx == 999)
               {
                std::cout << __FUNCTION__ << " ERROR: Could not match DUT " << fileDUTs.at(dut) << " for this run, skipping...." << std::endl;
                return false;
               }
            inputree->SetBranchAddress(Form("InterceptX_DUT%u", fileDUTs.at(dut)), &m_InterceptX_DUT.at(dutindx), &b_InterceptX_DUT.at(dutindx));
            inputree->SetBranchAddress(Form("InterceptY_DUT%u", fileDUTs.at(dut)), &m_InterceptY_DUT.at(dutindx), &b_InterceptY_DUT.at(dutindx));
            inputree->SetBranchAddress(Form("InterceptResX_DUT%u", fileDUTs.at(dut)), &m_InterceptResX_DUT.at(dutindx), &b_InterceptResX_DUT.at(dutindx));
            inputree->SetBranchAddress(Form("InterceptResY_DUT%u", fileDUTs.at(dut)), &m_InterceptResY_DUT.at(dutindx), &b_InterceptResY_DUT.at(dutindx));
           }
       for (unsigned int pln = 0; pln < fileplanes.size(); pln++)
           {
            unsigned int plnindx = 999;
            for (unsigned int s = 0; s < m_ntp_planes.size(); s++) { if (m_ntp_planes.at(s) == fileplanes.at(pln)) { plnindx = s; break; } }
            if (plnindx == 999)
               {
                std::cout << __FUNCTION__ << " ERROR: Could not match telescope plane " << fileplanes.at(pln) << " on this run, skipping...." << std::endl;
                return false;
               }
            inputree->SetBranchAddress(Form("InterceptX_TelPlane%u", fileplanes.at(pln)), &m_InterceptX_TelPlane.at(plnindx), &b_InterceptX_TelPlane.at(plnindx));
            inputree->SetBranchAddress(Form("InterceptY_TelPlane%u", fileplanes.at(pln)), &m_InterceptY_TelPlane.at(plnindx), &b_InterceptY_TelPlane.at(plnindx));
            inputree->SetBranchAddress(Form("InterceptResX_TelPlane%u", fileplanes.at(pln)), &m_InterceptResX_TelPlane.at(plnindx), &b_InterceptResX_TelPlane.at(plnindx));
            inputree->SetBranchAddress(Form("InterceptResY_TelPlane%u", fileplanes.at(pln)), &m_InterceptResY_TelPlane.at(plnindx), &b_InterceptResY_TelPlane.at(plnindx));
            inputree->SetBranchAddress(Form("ClusterSizeX_TelPlane%u", fileplanes.at(pln)), &m_ClusterSizeX_TelPlane.at(plnindx), &b_ClusterSizeX_TelPlane.at(plnindx));
            inputree->SetBranchAddress(Form("ClusterSizeY_TelPlane%u", fileplanes.at(pln)), &m_ClusterSizeY_TelPlane.at(plnindx), &b_ClusterSizeY_TelPlane.at(plnindx));
            inputree->SetBranchAddress(Form("ClusterNumPixels_TelPlane%u", fileplanes.at(pln)), &m_ClusterNumPixels_TelPlane.at(plnindx), &b_ClusterNumPixels_TelPlane.at(plnindx));
           }
      }

   // Bonn TBA Tracking barnches
   if (m_RunBase->LGADBase::GetTrackPackage() == "TBA")
      {
       // Add TBA relevent code
      }

   // Second stage analysis branches
   if (filelvl > 1)
      {
       // Chi2 for all fits
       ChTree->SetBranchAddress("chName", &m_chName, &b_chName);
       ChTree->SetBranchAddress("MaxIndxFtChi2", &m_MaxIndxFtChi2, &b_MaxIndxFtChi2);
       ChTree->SetBranchAddress("MinIndxFtChi2", &m_MinIndxFtChi2, &b_MinIndxFtChi2);
       ChTree->SetBranchAddress("MaxVoltFtChi2", &m_MaxVoltFtChi2, &b_MaxVoltFtChi2);
       ChTree->SetBranchAddress("MinVoltFtChi2", &m_MinVoltFtChi2, &b_MinVoltFtChi2);
       ChTree->SetBranchAddress("StrIndxFtChi2", &m_StrIndxFtChi2, &b_StrIndxFtChi2);
       ChTree->SetBranchAddress("EndIndxFtChi2", &m_EndIndxFtChi2, &b_EndIndxFtChi2);
       ChTree->SetBranchAddress("NoiseFtChi2", &m_NoiseFtChi2, &b_NoiseFtChi2);
       ChTree->SetBranchAddress("NoiseErrFtChi2", &m_NoiseErrFtChi2, &b_NoiseErrFtChi2);
       ChTree->SetBranchAddress("PedestFtChi2", &m_PedestFtChi2, &b_PedestFtChi2);
       ChTree->SetBranchAddress("PedestErrFtChi", &m_PedestErrFtChi2, &b_PedestErrFtChi2);
       ChTree->SetBranchAddress("MaxTimeFtChi2", &m_MaxTimeFtChi2, &b_MaxTimeFtChi2);
       ChTree->SetBranchAddress("MinTimeFtChi2", &m_MinTimeFtChi2, &b_MinTimeFtChi2);
       ChTree->SetBranchAddress("ChargeFtChi2", &m_ChargeFtChi2, &b_ChargeFtChi2);
       ChTree->SetBranchAddress("RiseTimeFtChi2", &m_RiseTimeFtChi2, &b_RiseTimeFtChi2);
       ChTree->SetBranchAddress("TriggTimeFtChi2", &m_TriggTimeFtChi2, &b_TriggTimeFtChi2);
       ChTree->SetBranchAddress("DVDTMaxFtChi2", &m_DVDTMaxFtChi2, &b_DVDTMaxFtChi2);
       ChTree->SetBranchAddress("TriggToTFtChi2", &m_TriggToTFtChi2, &b_TriggToTFtChi2);
       ChTree->SetBranchAddress("SignalFFTFtChi2", &m_SignalFFTFtChi2, &b_SignalFFTFtChi2);
       ChTree->SetBranchAddress("NoiseFFTFtChi2", &m_NoiseFFTFtChi2, &b_NoiseFFTFtChi2);
       // Fit parameters for each channel
       ChTree->SetBranchAddress("ChFitMaxVolt", &m_ChFitMxVlt.first, &b_ChFitMaxVolt);
       ChTree->SetBranchAddress("ChFitMaxVoltErr", &m_ChFitMxVlt.second, &b_ChFitMaxVoltErr);
       ChTree->SetBranchAddress("ChFitMinVolt", &m_ChFitMinVlt.first, &b_ChFitMinVolt);
       ChTree->SetBranchAddress("ChFitMinVoltErr", &m_ChFitMinVlt.second, &b_ChFitMinVoltErr);
       ChTree->SetBranchAddress("ChFitMxIndx", &m_ChFitMxIndx.first, &b_ChFitMxIndx);
       ChTree->SetBranchAddress("ChFitMxIndxErr", &m_ChFitMxIndx.second, &b_ChFitMxIndxErr);
       ChTree->SetBranchAddress("ChFitMinIndx", &m_ChFitMinIndx.first, &b_ChFitMinIndx);
       ChTree->SetBranchAddress("ChFitMinIndxErr", &m_ChFitMinIndx.second, &b_ChFitMinIndxErr);
       ChTree->SetBranchAddress("ChFitStrIndx", &m_ChFitStrIndx.first, &b_ChFitStrIndx);
       ChTree->SetBranchAddress("ChFitStrIndxErr", &m_ChFitStrIndx.second, &b_ChFitStrIndxErr);
       ChTree->SetBranchAddress("ChFitEndIndx", &m_ChFitEndIndx.first, &b_ChFitEndIndx);
       ChTree->SetBranchAddress("ChFitEndIndxErr", &m_ChFitEndIndx.second, &b_ChFitEndIndxErr);
       ChTree->SetBranchAddress("ChFitIsSig", &m_ChFitIsSig.first, &b_ChFitIsSig);
       ChTree->SetBranchAddress("ChFitIsSigErr", &m_ChFitIsSig.second, &b_ChFitIsSigErr);
       ChTree->SetBranchAddress("ChFitInWin", &m_ChFitInWin.first, &b_ChFitInWin);
       ChTree->SetBranchAddress("ChFitInWinErr", &m_ChFitInWin.second, &b_ChFitInWinErr);
       ChTree->SetBranchAddress("ChFitVStr", &m_ChFitVStr.first, &b_ChFitVStr);
       ChTree->SetBranchAddress("ChFitVStrErr", &m_ChFitVStr.second, &b_ChFitVStrErr);
       ChTree->SetBranchAddress("ChFitNoise", &m_ChFitNoise.first, &b_ChFitNoise);
       ChTree->SetBranchAddress("ChFitNoiseUnc", &m_ChFitNoise.second, &b_ChFitNoiseUnc);
       ChTree->SetBranchAddress("ChFitPedestal", &m_ChFitPedestal.first, &b_ChFitPedestal);
       ChTree->SetBranchAddress("ChFitPedestalUnc", &m_ChFitPedestal.second, &b_ChFitPedestalUnc);
       ChTree->SetBranchAddress("ChFitNoiseSlp", &m_ChFitNoiseSlp.first, &b_ChFitNoiseSlp);
       ChTree->SetBranchAddress("ChFitNoiseSlpErr", &m_ChFitNoiseSlp.second, &b_ChFitNoiseSlpErr);
       ChTree->SetBranchAddress("ChFitPedestalSlp", &m_ChFitPedestalSlp.first, &b_ChFitPedestalSlp);
       ChTree->SetBranchAddress("ChFitPedestalSlpErr", &m_ChFitPedestalSlp.second, &b_ChFitPedestalSlpErr);
       ChTree->SetBranchAddress("ChFitNoiseErr", &m_ChFitNoiseErr.first, &b_ChFitNoiseErr);
       ChTree->SetBranchAddress("ChFitNoiseErrUnc", &m_ChFitNoiseErr.second, &b_ChFitNoiseErrUnc);
       ChTree->SetBranchAddress("ChFitPedestalErr", &m_ChFitPedestalErr.first, &b_ChFitPedestalErr);
       ChTree->SetBranchAddress("ChFitPedestalErrUnc", &m_ChFitPedestalErr.second, &b_ChFitPedestalErrUnc);
       ChTree->SetBranchAddress("ChFitMaxTime", &m_ChFitMaxTime.first, &b_ChFitMaxTime);
       ChTree->SetBranchAddress("ChFitMaxTimeErr", &m_ChFitMaxTime.second, &b_ChFitMaxTimeErr);
       ChTree->SetBranchAddress("ChFitMinTime", &m_ChFitMinTime.first, &b_ChFitMinTime);
       ChTree->SetBranchAddress("ChFitMinTimeErr", &m_ChFitMinTime.second, &b_ChFitMinTimeErr);
       ChTree->SetBranchAddress("ChFitCharge", &m_ChFitCharge.first, &b_ChFitCharge);
       ChTree->SetBranchAddress("ChFitChargeErr", &m_ChFitCharge.second, &b_ChFitChargeErr);
       ChTree->SetBranchAddress("ChFitRiseT", &m_ChFitRiseT.first, &b_ChFitRiseT);
       ChTree->SetBranchAddress("ChFitRiseTErr", &m_ChFitRiseT.second, &b_ChFitRiseTErr);
       ChTree->SetBranchAddress("ChFitJitRiseSNR", &m_ChFitJitRiseSNR.first, &b_ChFitJitRiseSNR);
       ChTree->SetBranchAddress("ChFitJitRiseSNRErr", &m_ChFitJitRiseSNR.second, &b_ChFitJitRiseSNRErr);
       ChTree->SetBranchAddress("ChFitTriggTime", &m_ChFitTriggTime.first, &b_ChFitTriggTime);
       ChTree->SetBranchAddress("ChFitTriggTimeErr", &m_ChFitTriggTime.second, &b_ChFitTriggTimeErr);
       ChTree->SetBranchAddress("ChFitDVDTMax", &m_ChFitDVDTMax.first, &b_ChFitDVDTMax);
       ChTree->SetBranchAddress("ChFitDVDTMaxErr", &m_ChFitDVDTMax.second, &b_ChFitDVDTMaxErr);
       ChTree->SetBranchAddress("ChFitTriggToT", &m_ChFitTriggToT.first, &b_ChFitTriggToT);
       ChTree->SetBranchAddress("ChFitTriggToTErr", &m_ChFitTriggToT.second, &b_ChFitTriggToTErr);
       ChTree->SetBranchAddress("ChFitComplete", &m_ChFitComplete.first, &b_ChFitComplete);
       ChTree->SetBranchAddress("ChFitCompleteErr", &m_ChFitComplete.second, &b_ChFitCompleteErr);
       ChTree->SetBranchAddress("ChFitSoNR", &m_ChFitSoNR.first, &b_ChFitSoNR);
       ChTree->SetBranchAddress("ChFitSoNRErr", &m_ChFitSoNR.second, &b_ChFitSoNRErr);
       // FFT First stage
       if (FFTShape == 21 || FFTShape == 22)
          {
           ChTree->SetBranchAddress("ChFitSignalFFT", &m_ChFitSignalFFT.first, &b_ChFitSignalFFT);
           ChTree->SetBranchAddress("ChFitSignalFFTErr", &m_ChFitSignalFFT.second, &b_ChFitSignalFFTErr);
           ChTree->SetBranchAddress("ChFitNoiseFFT", &m_ChFitNoiseFFT.first, &b_ChFitNoiseFFT);
           ChTree->SetBranchAddress("ChFitNoiseFFTErr", &m_ChFitNoiseFFT.second, &b_ChFitNoiseFFT);
          }
       // CFD anf DV/DT vectors
       ChTree->SetBranchAddress("ChFitCFDTime", &m_ChFitCFDTime, &b_ChFitCFDTime);
       ChTree->SetBranchAddress("ChFitDVDTCFD", &m_ChFitDVDTCFD, &b_ChFitDVDTCFD);
       ChTree->SetBranchAddress("ChFitCFDToT", &m_ChFitCFDToT, &b_ChFitCFDToT);
       ChTree->SetBranchAddress("ChFitJitNdVdT", &m_ChFitJitNdVdT, &b_ChFitJitNdVdT);
       ChTree->SetBranchAddress("CFDTimeFtChi2", &m_CFDTimeFtChi2, &b_CFDTimeFtChi2);
       ChTree->SetBranchAddress("DVDTCFDFtChi2", &m_DVDTCFDFtChi2, &b_DVDTCFDFtChi2);
       ChTree->SetBranchAddress("CFDToTFtChi2", &m_CFDToTFtChi2, &b_CFDToTFtChi2);
       // Shape paramters for lvl2
       if (FFTShape == 22 || FFTShape == 12)
          {
           ShapeTree->SetBranchAddress("ChShapeName", &m_WvSh_chName, &b_WvSh_chName);
           ShapeTree->SetBranchAddress("ChShapePolarity", &m_ChShapePolarity, &b_ChShapePolarity);
           ShapeTree->SetBranchAddress("ChShapeMaxIndx", &m_ChShapeMaxIndx, &b_ChShapeMaxIndx);
           ShapeTree->SetBranchAddress("ChShapeMinIndx", &m_ChShapeMinIndx, &b_ChShapeMinIndx);
           ShapeTree->SetBranchAddress("ChShapeStrIndx", &m_ChShapeStrIndx, &b_ChShapeStrIndx);
           ShapeTree->SetBranchAddress("ChShapeEndIndx", &m_ChShapeEndIndx, &b_ChShapeEndIndx);
           ShapeTree->SetBranchAddress("ChShapeMax", &m_ChShapeMax, &b_ChShapeMax);
           ShapeTree->SetBranchAddress("ChShapeMin", &m_ChShapeMin, &b_ChShapeMin);
           ShapeTree->SetBranchAddress("ChShapeMaxTime", &m_ChShapeMaxTime, &b_ChShapeMaxTime);
           ShapeTree->SetBranchAddress("ChShapeMinTime", &m_ChShapeMinTime, &b_ChShapeMinTime);
           ShapeTree->SetBranchAddress("ChShapeCharge", &m_ChShapeCharge, &b_ChShapeCharge);
           ShapeTree->SetBranchAddress("ChShapeRiseTime", &m_ChShapeRiseTime, &b_ChShapeRiseTime);
           ShapeTree->SetBranchAddress("ChShapeJitRiseSNR", &m_ChShapeJitRiseSNR, &b_ChShapeJitRiseSNR);
           ShapeTree->SetBranchAddress("ChShapeTriggTime", &m_ChShapeTriggTime, &b_ChShapeTriggTime);
           ShapeTree->SetBranchAddress("ChShapeDVDTMax", &m_ChShapeDVDTMax, &b_ChShapeDVDTMax);
           ShapeTree->SetBranchAddress("ChShapeTriggToT", &m_ChShapeTriggToT, &b_ChShapeTriggToT);
           if (FFTShape == 22)
              {
               ShapeTree->SetBranchAddress("ChShapeSignalFFT", &m_ChShapeSignalFFT, &b_ChShapeSignalFFT);
               ShapeTree->SetBranchAddress("ChShapeNoiseFFT", &m_ChShapeNoiseFFT, &b_ChShapeNoiseFFT);
              }
           ShapeTree->SetBranchAddress("ChShapeNoise", &m_ChShapeNoise.first, &b_ChShapeNoise);
           ShapeTree->SetBranchAddress("ChShapeNoiseErr", &m_ChShapeNoise.second, &b_ChShapeNoiseErr);
           ShapeTree->SetBranchAddress("ChShapeSoN", &m_ChShapeSoN.first, &b_ChShapeSoN);
           ShapeTree->SetBranchAddress("ChShapeSoNUnc", &m_ChShapeSoN.second, &b_ChShapeSoNUnc);
           ShapeTree->SetBranchAddress("ChShapeCFDTime", &m_ChShapeCFDTime, &b_ChShapeCFDTime);
           ShapeTree->SetBranchAddress("ChShapeDVDTCFD", &m_ChShapeDVDTCFD, &b_ChShapeDVDTCFD);
           ShapeTree->SetBranchAddress("ChShapeCFDToT", &m_ChShapeCFDToT, &b_ChShapeCFDToT);
           ShapeTree->SetBranchAddress("ChShapeJitNdVdT", &m_ChShapeJitNdVdT, &b_ChShapeJitNdVdT);
           ShapeTree->SetBranchAddress("ChMeanVolt", &m_ChMeanVolt, &b_ChMeanVolt);
          }

       // Initialize 3rd stage analysis pointers
       if (filelvl > 2)
          {
           TimeResTree->SetBranchAddress("TmDiffChi2", &m_TmDiffChi2, &b_TmDiffChi2);
           TimeResTree->SetBranchAddress("TmDiffFitMean", &m_TmDiffFitMean, &b_TmDiffFitMean);
           TimeResTree->SetBranchAddress("TmDiffFitErr", &m_TmDiffFitErr, &b_TmDiffFitErr);
           TimeResTree->SetBranchAddress("CFD_ChA", &m_CFD_ChA, &b_CFD_ChA);
           TimeResTree->SetBranchAddress("ChA", &m_ChA, &b_ChA);
           TimeResTree->SetBranchAddress("ChB", &m_ChB, &b_ChB);
          }
      }

   return true;
}
// --------------------------------------------------------------------------------------------------------------
// Function to be called in every event to fetch the per event parameteres for processing - copying.
bool LGADRun::GetEvntBranches(Long64_t entry, unsigned int FFTShape, unsigned int filelvl, std::vector<unsigned int> filechannel, unsigned int fileindx, 
                              std::vector<unsigned int> fileplanes, std::vector<unsigned int> fileRefs, std::vector<unsigned int> fileDUTs)
{
   // Get general branches, needed in any config
   if (m_RunBase->LGADBase::GetInRootFile(fileindx) != m_RunBase->LGADBase::GetOutRootFile() || (m_RunBase->LGADBase::GetAnaStage() > 0 && filelvl == 0))
      {
       b_EvnNo->GetEntry(entry);
       for (unsigned int ich = 0; ich < filechannel.size(); ich++)
           {
            unsigned int chindx = 999;
            for (unsigned int s = 0; s < m_channels.size(); s++) if (m_channels.at(s) == filechannel.at(ich)) { chindx = s; break; }
            if (chindx == 999)
               {
                std::cout << __FUNCTION__ << " ERROR: Could not match channel " << filechannel.at(ich) << " for this file, skipping...." << std::endl;
                return false;
               }
            b_t.at(chindx)->GetEntry(entry);
            b_w.at(chindx)->GetEntry(entry);
            if (m_RunBase->LGADBase::GetInstrument() == Sampic)
               {
                b_ordrt.at(chindx)->GetEntry(entry);
                b_physt.at(chindx)->GetEntry(entry);
               }
            if (m_RunBase->LGADBase::GetInstrument() == TestBeamBin1 || m_RunBase->LGADBase::GetInstrument() == TestBeamBin2 || m_RunBase->LGADBase::GetInstrument() == LabTXT || 
                m_RunBase->LGADBase::GetInstrument() == TektronixScope || m_RunBase->LGADBase::GetInstrument() == LeCroyWRBin)
               {
                b_nPoints.at(chindx)->GetEntry(entry);
                b_SnRate.at(chindx)->GetEntry(entry);
                b_vScale.at(chindx)->GetEntry(entry);
                if (m_RunBase->LGADBase::GetInstrument() == TektronixScope) b_triggTime.at(chindx)->GetEntry(entry);
                if (m_RunBase->LGADBase::GetInstrument() == TestBeamBin1 || m_RunBase->LGADBase::GetInstrument() == TestBeamBin2 || m_RunBase->LGADBase::GetInstrument() == LeCroyWRBin) b_Scope.at(chindx)->GetEntry(entry);
               }
            if ((m_RunBase->LGADBase::GetInRootFile(fileindx) != m_RunBase->LGADBase::GetOutRootFile() || m_RunBase->LGADBase::GetAnaStage() > 1) && filelvl > 0)
               {
                b_EvPol.at(chindx)->GetEntry(entry);
                b_EvCharge.at(chindx)->GetEntry(entry);
                b_EvIsSignal.at(chindx)->GetEntry(entry);
                b_EvInWindow.at(chindx)->GetEntry(entry);
                b_EvVoltStr.at(chindx)->GetEntry(entry);
                b_EvMaxindx.at(chindx)->GetEntry(entry);
                b_EvMaxTime.at(chindx)->GetEntry(entry);
                b_EvMaxVolt.at(chindx)->GetEntry(entry);
                b_EvMinIndx.at(chindx)->GetEntry(entry);
                b_EvMinTime.at(chindx)->GetEntry(entry);
                b_EvMinVolt.at(chindx)->GetEntry(entry);
                b_EvStrIndx.at(chindx)->GetEntry(entry);
                b_EvEndIndx.at(chindx)->GetEntry(entry);
                b_EvNoise.at(chindx)->GetEntry(entry);
                b_EvPedest.at(chindx)->GetEntry(entry);
                b_EvNoiseErr.at(chindx)->GetEntry(entry);
                b_EvPedestErr.at(chindx)->GetEntry(entry);
                b_EvRiseT.at(chindx)->GetEntry(entry);
                b_EvJitRiseSNR.at(chindx)->GetEntry(entry);
                b_EvTriggTime.at(chindx)->GetEntry(entry);
                b_EvDVDTMax.at(chindx)->GetEntry(entry);
                b_EvTriggToT.at(chindx)->GetEntry(entry);
                b_EvComplete.at(chindx)->GetEntry(entry);
                if (FFTShape == 21 || FFTShape == 22)
                   {
                    b_EvSignalFFT.at(chindx)->GetEntry(entry);
                    b_EvNoiseFFT.at(chindx)->GetEntry(entry);
                   }
                b_EvVAdjCFD.at(chindx)->GetEntry(entry);
                b_EvCFDTime.at(chindx)->GetEntry(entry);
                b_EvDVDTCFD.at(chindx)->GetEntry(entry);
                b_EvCFDToT.at(chindx)->GetEntry(entry);
                b_EvJitNdVdT.at(chindx)->GetEntry(entry);
               }
           }
       if (m_RunBase->LGADBase::GetInstrument() == Sampic)
          {
           b_nPoints.at(0)->GetEntry(entry);
           b_SnRate.at(0)->GetEntry(entry);
          }
       if (m_RunBase->LGADBase::GetInstrument() == LabTXT || m_RunBase->LGADBase::GetInstrument() == TestBeamBin2) b_trigtime->GetEntry(entry);
       if (m_RunBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Basic branches retrieved!" << std::endl;
      }

   // Recover traking branches either if file is cpoied and they are already contained inside or if file lvl is < 4 and the analysis level is 4th stage
   if (m_RunBase->LGADBase::GetInRootFile(fileindx) != m_RunBase->LGADBase::GetOutRootFile() || (filelvl < 4 && m_RunBase->LGADBase::GetAnaStage() > 3))
      {
       // EU Telescope branches
       if (m_RunBase->LGADBase::GetTrackPackage() == "EUTelescope")
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
       if (m_RunBase->LGADBase::GetTrackPackage() == "TBA")
          {
           //  Add TBA relevent code
          }
       // Cory Branches
       if (m_RunBase->LGADBase::GetTrackPackage() == "Cory")
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
           for (unsigned int ref = 0; ref < fileRefs.size(); ref++)
               {
                unsigned int refindx = 999;
                for (unsigned int s = 0; s < m_ntp_Ref.size(); s++) { if (m_ntp_Ref.at(s) == fileRefs.at(ref)) { refindx = s; break; } }
                if (refindx == 999)
                   {
                    std::cout << __FUNCTION__ << " ERROR: Could not match refernce plane " << fileRefs.at(ref) << " for this file, skipping...." << std::endl;
                    return false;
                   }
                b_InterceptX_Ref.at(refindx)->GetEntry(entry);
                b_InterceptY_Ref.at(refindx)->GetEntry(entry);
                b_InterceptResX_Ref.at(refindx)->GetEntry(entry);
                b_InterceptResY_Ref.at(refindx)->GetEntry(entry);
               }
           for (unsigned int dut = 0; dut < fileDUTs.size(); dut++)
               {
                unsigned int dutindx = 999;
                for (unsigned int s = 0; s < m_ntp_DUT.size(); s++) { if (m_ntp_DUT.at(s) == fileDUTs.at(dut)) { dutindx = s; break; } }
                if (dutindx == 999)
                   {
                    std::cout << __FUNCTION__ << " ERROR: Could not match DUT " << fileDUTs.at(dut) << " for this file, skipping...." << std::endl;
                    return false;
                   }
                b_InterceptX_DUT.at(dutindx)->GetEntry(entry);
                b_InterceptY_DUT.at(dutindx)->GetEntry(entry);
                b_InterceptResX_DUT.at(dutindx)->GetEntry(entry);
                b_InterceptResY_DUT.at(dutindx)->GetEntry(entry);
               }
           for (unsigned int pln = 0; pln < fileplanes.size(); pln++)
               {
                unsigned int plnindx = 999;
                for (unsigned int s = 0; s < m_ntp_planes.size(); s++) { if (m_ntp_planes.at(s) == fileplanes.at(pln)) { plnindx = s; break; } }
                if (plnindx == 999)
                   {
                    std::cout << __FUNCTION__ << " ERROR: Could not match telescope plane " << fileplanes.at(pln) << " for this file, skipping...." << std::endl;
                    return false;
                   }
                b_InterceptX_TelPlane.at(plnindx)->GetEntry(entry);
                b_InterceptY_TelPlane.at(plnindx)->GetEntry(entry);
                b_InterceptResX_TelPlane.at(plnindx)->GetEntry(entry);
                b_InterceptResY_TelPlane.at(plnindx)->GetEntry(entry);
                b_ClusterSizeX_TelPlane.at(plnindx)->GetEntry(entry);
                b_ClusterSizeY_TelPlane.at(plnindx)->GetEntry(entry);
                b_ClusterNumPixels_TelPlane.at(plnindx)->GetEntry(entry);
               }
          }
       if (m_RunBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Tracking branches retrieved!" << std::endl;
      }

   return true;
}
// --------------------------------------------------------------------------------------------------------------
// Function to be called at the second or third stage analysis to fetch the corresponding parameters from the Ntuple.
// THis function need the stage as an argument and should be called for every channel at the second stage or for 
// every CFD percentage at the third stage.
bool LGADRun::GetChanelBranches(Long64_t channel, unsigned int fileindx, unsigned int stage, unsigned int FFTShape, unsigned int filelvl)
{
   // Second stage analysis branches
   if (filelvl > 1)
      {
       // Get per channel parameters
       if (stage == 1 || m_RunBase->LGADBase::GetInRootFile(fileindx) != m_RunBase->LGADBase::GetOutRootFile())
          { 
           b_chName->GetEntry(channel);
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
           if (FFTShape == 21 || FFTShape == 22)
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
           if (FFTShape == 22 || FFTShape == 12)
              {
               b_WvSh_chName->GetEntry(channel);
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
               if (FFTShape == 22)
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
       else if (filelvl > 2 && (stage == 2 || m_RunBase->LGADBase::GetInRootFile(fileindx) != m_RunBase->LGADBase::GetOutRootFile()))
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
// Function to be called at the begining after the InitNtuple funvtion, to setup the required amount of pointers
// for the output file. Weather this will be just an addition on the existing or a complete copy is handled internaly.
bool LGADRun::InitOutputBranches()
{
    // Typically the following branches always exist in the input file. THe only reason to re-intialise them for the output is copying
    if (m_RunBase->LGADBase::GetInRootFile(0) != m_RunBase->LGADBase::GetOutRootFile() || (m_RunBase->LGADBase::GetInRootFiles()).size() > 1)
       {
        // Initialize the branch pointers
        b_t_out.clear();
        b_w_out.clear();
        b_t_out.resize((m_RunBase->LGADBase::GetChannels()).size());
        b_w_out.resize((m_RunBase->LGADBase::GetChannels()).size());
        if (m_RunBase->LGADBase::GetInstrument() == Sampic)
           {
            b_ordrt_out.clear();
            b_physt_out.clear();
            b_nPoints_out.clear();
            b_SnRate_out.clear();
            b_ordrt_out.resize((m_RunBase->LGADBase::GetChannels()).size());
            b_physt_out.resize((m_RunBase->LGADBase::GetChannels()).size());
            b_nPoints_out.resize(1, 0);
            b_SnRate_out.resize(1, 0);
           }
        if (m_RunBase->LGADBase::GetInstrument() == TestBeamBin1 || m_RunBase->LGADBase::GetInstrument() == TestBeamBin2 || m_RunBase->LGADBase::GetInstrument() == LabTXT ||
            m_RunBase->LGADBase::GetInstrument() == TektronixScope || m_RunBase->LGADBase::GetInstrument() == LeCroyWRBin)
           {
            b_nPoints_out.clear();
            b_SnRate_out.clear();
            b_vScale_out.clear();
            b_nPoints_out.resize((m_RunBase->LGADBase::GetChannels()).size());
            b_SnRate_out.resize((m_RunBase->LGADBase::GetChannels()).size());
            b_vScale_out.resize((m_RunBase->LGADBase::GetChannels()).size());
            if (m_RunBase->LGADBase::GetInstrument() == TektronixScope)
               {
                b_triggTime_out.clear();
                b_triggTime_out.resize((m_RunBase->LGADBase::GetChannels()).size());
               }
            if (m_RunBase->LGADBase::GetInstrument() == TestBeamBin1 || m_RunBase->LGADBase::GetInstrument() == TestBeamBin2 || m_RunBase->LGADBase::GetInstrument() == LeCroyWRBin)
               {
                b_Scope_out.clear();
                b_Scope_out.resize((m_RunBase->LGADBase::GetChannels()).size());
               }
           }
       }

   if ((m_RunBase->LGADBase::GetAnaStage() > 0 && m_RunBase->LGADBase::GetFileLVL() == 0) ||
       (m_RunBase->LGADBase::GetFileLVL() > 0 && (m_RunBase->LGADBase::GetInRootFile(0) != m_RunBase->LGADBase::GetOutRootFile() || (m_RunBase->LGADBase::GetInRootFiles()).size() > 1)))
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
       b_EvPol_out.resize((m_RunBase->LGADBase::GetChannels()).size(), nullptr);
       b_EvCharge_out.resize((m_RunBase->LGADBase::GetChannels()).size(), nullptr);
       b_EvIsSignal_out.resize((m_RunBase->LGADBase::GetChannels()).size(), nullptr);
       b_EvInWindow_out.resize((m_RunBase->LGADBase::GetChannels()).size(), nullptr);
       b_EvVoltStr_out.resize((m_RunBase->LGADBase::GetChannels()).size(), nullptr);
       b_EvMaxindx_out.resize((m_RunBase->LGADBase::GetChannels()).size(), nullptr);
       b_EvMaxTime_out.resize((m_RunBase->LGADBase::GetChannels()).size(), nullptr);
       b_EvMaxVolt_out.resize((m_RunBase->LGADBase::GetChannels()).size(), nullptr);
       b_EvMinIndx_out.resize((m_RunBase->LGADBase::GetChannels()).size(), nullptr);
       b_EvMinTime_out.resize((m_RunBase->LGADBase::GetChannels()).size(), nullptr);
       b_EvMinVolt_out.resize((m_RunBase->LGADBase::GetChannels()).size(), nullptr);
       b_EvStrIndx_out.resize((m_RunBase->LGADBase::GetChannels()).size(), nullptr);
       b_EvEndIndx_out.resize((m_RunBase->LGADBase::GetChannels()).size(), nullptr);
       b_EvNoise_out.resize((m_RunBase->LGADBase::GetChannels()).size(), nullptr);
       b_EvPedest_out.resize((m_RunBase->LGADBase::GetChannels()).size(), nullptr);
       b_EvNoiseErr_out.resize((m_RunBase->LGADBase::GetChannels()).size(), nullptr);
       b_EvPedestErr_out.resize((m_RunBase->LGADBase::GetChannels()).size(), nullptr);
       b_EvRiseT_out.resize((m_RunBase->LGADBase::GetChannels()).size(), nullptr);
       b_EvJitNdVdT_out.resize((m_RunBase->LGADBase::GetChannels()).size(), nullptr);
       b_EvJitRiseSNR_out.resize((m_RunBase->LGADBase::GetChannels()).size(), nullptr);
       b_EvTriggTime_out.resize((m_RunBase->LGADBase::GetChannels()).size(), nullptr);
       b_EvDVDTMax_out.resize((m_RunBase->LGADBase::GetChannels()).size(), nullptr);
       b_EvDVDTCFD_out.resize((m_RunBase->LGADBase::GetChannels()).size(), nullptr);
       b_EvCFDToT_out.resize((m_RunBase->LGADBase::GetChannels()).size(), nullptr);
       b_EvTriggToT_out.resize((m_RunBase->LGADBase::GetChannels()).size(), nullptr);
       if ((m_RunBase->LGADBase::GetAnaStage() > 0 && m_RunBase->LGADBase::GetFileLVL() == 0 && m_RunBase->LGADBase::GetDoFFT()) ||
           ((m_RunBase->LGADBase::GetInRootFile(0) != m_RunBase->LGADBase::GetOutRootFile() || (m_RunBase->LGADBase::GetInRootFiles()).size() > 1) && m_RunBase->LGADBase::GetHasFFT()))
          {
           b_EvSignalFFT_out.clear();
           b_EvNoiseFFT_out.clear();
           b_EvSignalFFT_out.resize((m_RunBase->LGADBase::GetChannels()).size(), nullptr);
           b_EvNoiseFFT_out.resize((m_RunBase->LGADBase::GetChannels()).size(), nullptr);
          }
       b_EvComplete_out.resize((m_RunBase->LGADBase::GetChannels()).size(), nullptr);
       b_EvVAdjCFD_out.resize((m_RunBase->LGADBase::GetChannels()).size(), nullptr);
       b_EvCFDTime_out.resize((m_RunBase->LGADBase::GetChannels()).size(), nullptr);
      }

   // we only intialise these pointers if we need to copy, if the 
  if ((m_RunBase->LGADBase::GetInRootFile(0) != m_RunBase->LGADBase::GetOutRootFile() || (m_RunBase->LGADBase::GetInRootFiles()).size() > 1) && 
      m_RunBase->LGADBase::GetHasTrck() && m_RunBase->LGADBase::GetTrackPackage() != "undef" && m_RunBase->LGADBase::GetTrackComb())
      { 
       if (m_RunBase->LGADBase::GetTrackPackage() == "Cory")
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
          }
       else if (m_RunBase->LGADBase::GetTrackPackage() == "TBA")
               {
                //  TBA Branches should go here
               }
      }
   // No need for initialization of second and third stage pointers here, the structure of these trees with a channel as an event
   // automatically acounts for the varied number of channels per file, since one branch is globaly needed for each quntity.

   return true;
}
// --------------------------------------------------------------------------------------------------------------
// Function to be called for every Ntuple if the framework is writing the analysis in a different file than  the input data.
// THis function sets the correct pointers adress for the file that is opened to the output Ntuple. Function to be called for 
// every new file. It has to be called with filelvl <= 1 at the analysis stage and again with filelvl > 1 at the second/third analysis stage.
bool LGADRun::SetOutputBranchAdresses(unsigned int FFTShape, unsigned int filelvl, std::vector<unsigned int> filechannel, unsigned int fileindx,
                                      std::vector<unsigned int> fileplanes, std::vector<unsigned int> fileRefs, std::vector<unsigned int> fileDUTs)
{
   // Variables to reference for non-existant branches that need to be populated
   double no_ch_double = -99;
   int no_ch_int = -99;
   unsigned int no_ch_uint = 99;
   float no_ch_flot = -99.;
   Long64_t no_ch_long = -99;
   std::vector<double> no_ch_doub_vec;
   no_ch_doub_vec.push_back(-99);
   std::vector<int> no_ch_int_vec;
   no_ch_int_vec.push_back(-99);

   // Only propagate the event number if you are copying to a new file
   if (m_RunBase->LGADBase::GetInRootFile(fileindx) != m_RunBase->LGADBase::GetOutRootFile()) b_EvnNo_out = m_OutTree->Branch("EvnNo", &EvnNo, "EvnNo/i");
   // Do channel matching and assign branches in each case
   for (unsigned int a = 0; a < filechannel.size(); a++)
       {
        unsigned int chindx = 999;
        for (unsigned int k = 0; k < m_RunBase->LGADBase::GetChannels().size(); k++) { if (m_RunBase->LGADBase::GetChannels().at(k) == filechannel.at(a)) { chindx = k; break; } }
        if (chindx == 999)
           {
            std::cout << __FUNCTION__ << " ERROR: Could not match channel " << filechannel.at(a) << " of this file, skipping...." << std::endl;
            return false;
           }
        if (m_RunBase->LGADBase::GetInRootFile(fileindx) != m_RunBase->LGADBase::GetOutRootFile())
           {
            if (m_RunBase->LGADBase::GetInstrument() == Sampic)
               {
                b_physt_out.at(chindx) = m_OutTree->Branch(Form("physt%02u", filechannel.at(a)), &m_physt.at(chindx), "m_physt.at(chindx)/D");
                b_ordrt_out.at(chindx) = m_OutTree->Branch(Form("ordrt%02u", filechannel.at(a)), &m_ordrt.at(chindx), "m_ordrt.at(chindx)/D");
               }
            b_t_out.at(chindx) = m_OutTree->Branch(Form("t%02u", filechannel.at(a)), &m_t.at(chindx));
            b_w_out.at(chindx) = m_OutTree->Branch(Form("w%02u", filechannel.at(a)), &m_w.at(chindx));
            if (m_RunBase->LGADBase::GetInstrument() == TestBeamBin1 || m_RunBase->LGADBase::GetInstrument() == TestBeamBin2 ||
                m_RunBase->LGADBase::GetInstrument() == LabTXT || m_RunBase->LGADBase::GetInstrument() == TektronixScope ||
                m_RunBase->LGADBase::GetInstrument() == LeCroyWRBin)
               {
                b_vScale_out.at(chindx) = m_OutTree->Branch(Form("vScale%02u", filechannel.at(a)), &m_scale.at(chindx), "m_scale.at(chindx)/F");
                b_nPoints_out.at(chindx) = m_OutTree->Branch(Form("nPoints%02u", filechannel.at(a)), &m_npoints.at(chindx), "m_npoints.at(chindx)/i");
                b_SnRate_out.at(chindx) = m_OutTree->Branch(Form("SnRate%02u", filechannel.at(a)), &m_srate.at(chindx), "m_srate.at(chindx)/L");
                if (m_RunBase->LGADBase::GetInstrument() == TektronixScope)
                    b_triggTime_out.at(chindx) = m_OutTree->Branch(Form("triggTime%02u", filechannel.at(a)), &m_triggTime.at(chindx), "m_triggTime.at(chindx)/D");
                if (m_RunBase->LGADBase::GetInstrument() == TestBeamBin1 || m_RunBase->LGADBase::GetInstrument() == TestBeamBin2)
                    b_Scope_out.at(chindx) = m_OutTree->Branch(Form("Scope%02u", filechannel.at(a)), &m_scope.at(chindx), "m_scope.at(chindx)/i");
               }
           }
        // Level 1 analysis requirements
        if ((m_RunBase->LGADBase::GetAnaStage() > 0 && filelvl == 0) || (filelvl > 0 && m_RunBase->LGADBase::GetInRootFile(fileindx) != m_RunBase->LGADBase::GetOutRootFile()))
           {
            b_EvPol_out.at(chindx) = m_OutTree->Branch(Form("Polarity%02u", filechannel.at(a)), &m_EvPol.at(chindx), "m_EvPol.at(chindx)/I");
            b_EvCharge_out.at(chindx) = m_OutTree->Branch(Form("Charge%02u", filechannel.at(a)), &m_EvCharge.at(chindx), "m_EvCharge.at(chindx)/D");
            b_EvIsSignal_out.at(chindx) = m_OutTree->Branch(Form("IsSignal%02u", filechannel.at(a)), &m_EvIsSignal.at(chindx), "m_EvIsSignal.at(chindx)/I");
            b_EvInWindow_out.at(chindx) = m_OutTree->Branch(Form("InWindow%02u", filechannel.at(a)), &m_EvInWindow.at(chindx), "m_EvInWindow.at(chindx)/I");
            b_EvVoltStr_out.at(chindx) = m_OutTree->Branch(Form("VoltStr%02u", filechannel.at(a)), &m_EvVoltStr.at(chindx), "m_EvVoltStr.at(chindx)/I");
            b_EvMaxindx_out.at(chindx) = m_OutTree->Branch(Form("MaxIndx%02u", filechannel.at(a)), &m_EvMaxindx.at(chindx), "m_EvMaxindx.at(chindx)/I");
            b_EvMaxTime_out.at(chindx) = m_OutTree->Branch(Form("MaxTime%02u", filechannel.at(a)), &m_EvMaxTime.at(chindx), "m_EvMaxTime.at(chindx)/D");
            b_EvMaxVolt_out.at(chindx) = m_OutTree->Branch(Form("MaxVolt%02u", filechannel.at(a)), &m_EvMaxVolt.at(chindx), "m_EvMaxVolt.at(chindx)/D");
            b_EvMinIndx_out.at(chindx) = m_OutTree->Branch(Form("MinIndx%02u", filechannel.at(a)), &m_EvMinIndx.at(chindx), "m_EvMinIndx.at(chindx)/I");
            b_EvMinTime_out.at(chindx) = m_OutTree->Branch(Form("MinTime%02u", filechannel.at(a)), &m_EvMinTime.at(chindx), "m_MinTime.at(chindx)/D");
            b_EvMinVolt_out.at(chindx) = m_OutTree->Branch(Form("MinVolt%02u", filechannel.at(a)), &m_EvMinVolt.at(chindx), "m_EvMinVolt.at(chindx)/D");
            b_EvStrIndx_out.at(chindx) = m_OutTree->Branch(Form("Strindx%02u", filechannel.at(a)), &m_EvStrIndx.at(chindx), "m_EvStrIndx.at(chindx)/I");
            b_EvEndIndx_out.at(chindx) = m_OutTree->Branch(Form("EndIndx%02u", filechannel.at(a)), &m_EvEndIndx.at(chindx), "m_EvEndIndx.at(chindx)/I");
            b_EvNoise_out.at(chindx) = m_OutTree->Branch(Form("Noise%02u", filechannel.at(a)), &m_EvNoise.at(chindx), "m_EvNoise.at(chindx)/D");
            b_EvPedest_out.at(chindx) = m_OutTree->Branch(Form("Pedest%02u", filechannel.at(a)), &m_EvPedest.at(chindx), "m_EvPedest.at(chindx)/D");
            b_EvNoiseErr_out.at(chindx) = m_OutTree->Branch(Form("NoiseErr%02u", filechannel.at(a)), &m_EvNoiseErr.at(chindx), "m_EvNoiseErr.at(chindx)/D");
            b_EvPedestErr_out.at(chindx) = m_OutTree->Branch(Form("PedestErr%02u", filechannel.at(a)), &m_EvPedestErr.at(chindx), "m_EvPedestEr.at(chindx)/D");
            b_EvRiseT_out.at(chindx) = m_OutTree->Branch(Form("RiseT%02u", filechannel.at(a)), &m_EvRiseT.at(chindx), "m_EvRiseTe.at(chindx)/D");
            b_EvJitRiseSNR_out.at(chindx) = m_OutTree->Branch(Form("JitRiseSNR%02u", filechannel.at(a)), &m_EvJitRiseSNR.at(chindx), "m_EvJitRiseSNR.at(chindx)/D");
            b_EvTriggTime_out.at(chindx) = m_OutTree->Branch(Form("TriggTime%02u", filechannel.at(a)), &m_EvTriggTime.at(chindx), "m_EvTriggTime.at(chindx)/D");
            b_EvDVDTMax_out.at(chindx) = m_OutTree->Branch(Form("DVDTMax%02u", filechannel.at(a)), &m_EvDVDTMax.at(chindx), "m_EvDVDTMax.at(chindx)/D");
            b_EvTriggToT_out.at(chindx) = m_OutTree->Branch(Form("TriggToT%02u", filechannel.at(a)), &m_EvTriggToT.at(chindx), "m_EvTriggToT.at(chindx)/D");
            if (((FFTShape == 21 || FFTShape == 22) && m_RunBase->LGADBase::GetInRootFile(fileindx) != m_RunBase->LGADBase::GetOutRootFile()) ||
                (m_RunBase->LGADBase::GetAnaStage() > 0 && filelvl == 0 && m_RunBase->LGADBase::GetDoFFT()))
               {
                b_EvSignalFFT_out.at(chindx) = m_OutTree->Branch(Form("SignalFFT%02u", filechannel.at(a)), &m_EvSignalFFT.at(chindx), "m_EvSignalFFT.at(chindx)/D");
                b_EvNoiseFFT_out.at(chindx) = m_OutTree->Branch(Form("NoiseFFT%02u", filechannel.at(a)), &m_EvNoiseFFT.at(chindx), "m_EvNoiseFFT.at(chindx)/D");
               }
            b_EvComplete_out.at(chindx) = m_OutTree->Branch(Form("Complete%02u", filechannel.at(a)), &m_EvComplete.at(chindx), "m_EvComplete.at(chindx)/I");
            b_EvVAdjCFD_out.at(chindx) = m_OutTree->Branch(Form("VAdjCFD%02u", filechannel.at(a)), &m_EvVAdjCFD.at(chindx));
            b_EvCFDTime_out.at(chindx) = m_OutTree->Branch(Form("CFDTime%02u", filechannel.at(a)), &m_EvCFDTime.at(chindx));
            b_EvDVDTCFD_out.at(chindx) = m_OutTree->Branch(Form("DVDTCFD%02u", filechannel.at(a)), &m_EvDVDTCFD.at(chindx));
            b_EvCFDToT_out.at(chindx) = m_OutTree->Branch(Form("CFDToT%02u", filechannel.at(a)), &m_EvCFDToT.at(chindx));
            b_EvJitNdVdT_out.at(chindx) = m_OutTree->Branch(Form("JitNdVdT%02u", filechannel.at(a)), &m_EvJitNdVdT.at(chindx));
           }
       }
   // Non existant channels should be coupled to values that make no sens for filling the tree corectly at the ntuple
   std::vector<unsigned int> nochan;
   bool found; // variable to be re-used in this function for all othervectors to be populated
   for (unsigned int k = 0; k < m_RunBase->LGADBase::GetChannels().size(); k++)
       {
        found = false;
        for (unsigned int hgp = 0; hgp < filechannel.size(); hgp++) { if (m_RunBase->LGADBase::GetChannels().at(k) == filechannel.at(hgp)) { found = true; break; } }
        if (!found) nochan.push_back(k);
       }
   if (nochan.size() > 0)
      {
       for (unsigned int gt = 0; gt < nochan.size(); gt++)
           {
            // These branches are to be filled only if saving output to a new file, if not don't migrate them 
            if (m_RunBase->LGADBase::GetInRootFile(fileindx) != m_RunBase->LGADBase::GetOutRootFile())
               {
                if (m_RunBase->LGADBase::GetInstrument() == Sampic)
                   {
                    b_physt_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("physt%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_double, "no_ch_double/D");
                    b_ordrt_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("ordrt%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_double, "no_ch_double/D");
                   }
                b_t_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("t%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_doub_vec);
                b_w_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("w%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_doub_vec);
                if (m_RunBase->LGADBase::GetInstrument() == TestBeamBin1 || m_RunBase->LGADBase::GetInstrument() == TestBeamBin2 ||
                    m_RunBase->LGADBase::GetInstrument() == LabTXT || m_RunBase->LGADBase::GetInstrument() == TektronixScope ||
                    m_RunBase->LGADBase::GetInstrument() == LeCroyWRBin)
                   {
                    b_vScale_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("vScale%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_flot, "no_ch_flot/F");
                    b_nPoints_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("nPoints%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_uint, "no_ch_uint/i");
                    b_SnRate_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("SnRate%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_long, "no_ch_long/L");
                    if (m_RunBase->LGADBase::GetInstrument() == TektronixScope) b_triggTime_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("triggTime%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_double, "no_ch_double/D");
                    if (m_RunBase->LGADBase::GetInstrument() == TestBeamBin1 || m_RunBase->LGADBase::GetInstrument() == TestBeamBin2)
                        b_Scope_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("Scope%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_uint, "no_ch_uint/i");
                   }
               }
            if ((m_RunBase->LGADBase::GetAnaStage() > 0 && filelvl == 0) || (filelvl > 0 && m_RunBase->LGADBase::GetInRootFile(fileindx) != m_RunBase->LGADBase::GetOutRootFile()))
               {
                b_EvPol_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("Polarity%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_int, "no_ch_int/I");
                b_EvCharge_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("Charge%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_double, "no_ch_double/D");
                b_EvIsSignal_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("IsSignal%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_int, "no_ch_int/I");
                b_EvInWindow_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("InWindow%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_int, "no_ch_int/I");
                b_EvVoltStr_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("VoltStr%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_int, "no_ch_int/I");
                b_EvMaxindx_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("MaxIndx%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_int, "no_ch_int/I");
                b_EvMaxTime_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("MaxTime%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_double, "no_ch_double/D");
                b_EvMaxVolt_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("MaxVolt%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_double, "no_ch_double/D");
                b_EvMinIndx_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("MinIndx%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_int, "no_ch_int/I");
                b_EvMinTime_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("MinTime%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_double, "no_ch_double/D");
                b_EvMinVolt_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("MinVolt%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_double, "no_ch_double/D");
                b_EvStrIndx_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("Strindx%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_int, "no_ch_int/I");
                b_EvEndIndx_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("EndIndx%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_int, "no_ch_int/I");
                b_EvNoise_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("Noise%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_double, "/D");
                b_EvPedest_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("Pedest%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_double, "no_ch_double/D");
                b_EvNoiseErr_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("NoiseErr%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_double, "no_ch_double/D");
                b_EvPedestErr_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("PedestErr%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_double, "no_ch_double/D");
                b_EvRiseT_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("RiseT%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_double, "no_ch_double/D");
                b_EvJitRiseSNR_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("JitRiseSNR%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_double, "no_ch_double/D");
                b_EvTriggTime_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("TriggTime%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_double, "no_ch_double/D");
                b_EvDVDTMax_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("DVDTMax%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_double, "no_ch_double/D");
                b_EvTriggToT_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("TriggToT%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_double, "no_ch_double/D");
                if (((FFTShape == 21 || FFTShape == 22) && m_RunBase->LGADBase::GetInRootFile(fileindx) != m_RunBase->LGADBase::GetOutRootFile()) ||
                    (m_RunBase->LGADBase::GetAnaStage() > 0 && filelvl == 0 && m_RunBase->LGADBase::GetDoFFT()))
                   {
                    b_EvSignalFFT_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("SignalFFT%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_double, "no_ch_double/D");
                    b_EvNoiseFFT_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("NoiseFFT%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_double, "no_ch_double/D");
                   }
                b_EvComplete_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("Complete%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_int, "no_ch_int/I");
                b_EvVAdjCFD_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("VAdjCFD%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_doub_vec);
                b_EvCFDTime_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("CFDTime%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_doub_vec);
                b_EvDVDTCFD_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("DVDTCFD%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_doub_vec);
                b_EvCFDToT_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("CFDToT%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_doub_vec);
                b_EvJitNdVdT_out.at(nochan.at(gt)) = m_OutTree->Branch(Form("JitNdVdT%02u", (m_RunBase->LGADBase::GetChannels()).at(nochan.at(gt))), &no_ch_doub_vec);
               }
           }
      }

   // Tracking part of the file, only prppagated when copying. If tracking combination is happening at this stage, the TrackCombine functinos will create the necessery brancehs.
   if (m_RunBase->LGADBase::GetInRootFile(fileindx) != m_RunBase->LGADBase::GetOutRootFile() && m_RunBase->LGADBase::GetTrackComb() && m_RunBase->LGADBase::GetTrackPackage() != "undef")
      { 
       if (m_RunBase->LGADBase::GetTrackPackage() == "Cory")
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
           for (unsigned int refs = 0; refs < fileRefs.size(); refs++)
               {
                unsigned int refsindx = 999;
                for (unsigned int k = 0; k < m_ntp_Ref.size(); k++) { if (m_ntp_Ref.at(k) == fileRefs.at(refs)) { refsindx = k; break; } }
                if (refsindx == 999)
                   {
                    std::cout << __FUNCTION__ << " ERROR: Could not match reference " << fileRefs.at(refs) << " of this file, skipping...." << std::endl;
                    return false;
                   }
                b_InterceptX_Ref_out.at(refsindx) = m_OutTree->Branch(Form("InterceptX_Ref%u", fileRefs.at(refs)), &m_InterceptX_Ref.at(refsindx));
                b_InterceptY_Ref_out.at(refsindx) = m_OutTree->Branch(Form("InterceptY_Ref%u", fileRefs.at(refs)), &m_InterceptY_Ref.at(refsindx));
                b_InterceptResX_Ref_out.at(refsindx) = m_OutTree->Branch(Form("InterceptResX_Ref%u", fileRefs.at(refs)), &m_InterceptResX_Ref.at(refsindx));
                b_InterceptResY_Ref_out.at(refsindx) = m_OutTree->Branch(Form("InterceptResY_Ref%u", fileRefs.at(refs)), &m_InterceptResY_Ref.at(refsindx));
               }
           for (unsigned int duts = 0; duts < fileDUTs.size(); duts++)
               {
                unsigned int dutindx = 999;
                for (unsigned int k = 0; k < m_ntp_DUT.size(); k++) { if (m_ntp_DUT.at(k) == fileDUTs.at(duts)) { dutindx = k; break; } }
                if (dutindx == 999)
                   {
                    std::cout << __FUNCTION__ << " ERROR: Could not match DUT " << fileDUTs.at(duts) << " of this file, skipping...." << std::endl;
                    return false;
                   }
                b_InterceptX_DUT_out.at(dutindx) = m_OutTree->Branch(Form("InterceptX_DUT%u", fileDUTs.at(duts)), &m_InterceptX_DUT.at(dutindx));
                b_InterceptY_DUT_out.at(dutindx) = m_OutTree->Branch(Form("InterceptY_DUT%u", fileDUTs.at(duts)), &m_InterceptY_DUT.at(dutindx));
                b_InterceptResX_DUT_out.at(dutindx) = m_OutTree->Branch(Form("InterceptResX_DUT%u", fileDUTs.at(duts)), &m_InterceptResX_DUT.at(dutindx));
                b_InterceptResY_DUT_out.at(dutindx) = m_OutTree->Branch(Form("InterceptResY_DUT%u", fileDUTs.at(duts)), &m_InterceptResY_DUT.at(dutindx));
               }
           for (unsigned int pls = 0; pls < fileplanes.size(); pls++)
               {
                unsigned int plindx = 999;
                for (unsigned int k = 0; k < m_ntp_planes.size(); k++) { if (m_ntp_planes.at(k) == fileplanes.at(pls)) { plindx = k; break; } }
                if (plindx == 999)
                   {
                    std::cout << __FUNCTION__ << " ERROR: Could not match DUT " << fileplanes.at(pls) << " of this file, skipping...." << std::endl;
                    return false;
                   }   
                b_InterceptX_TelPlane_out.at(plindx) = m_OutTree->Branch(Form("InterceptX_TelPlane%u", fileplanes.at(pls)), &m_InterceptX_TelPlane.at(plindx));
                b_InterceptY_TelPlane_out.at(plindx) = m_OutTree->Branch(Form("InterceptY_TelPlane%u", fileplanes.at(pls)), &m_InterceptY_TelPlane.at(plindx));
                b_InterceptResX_TelPlane_out.at(plindx) = m_OutTree->Branch(Form("InterceptResX_TelPlane%u", fileplanes.at(pls)), &m_InterceptResX_TelPlane.at(plindx));
                b_InterceptResY_TelPlane_out.at(plindx) = m_OutTree->Branch(Form("InterceptResY_TelPlane%u", fileplanes.at(pls)), &m_InterceptResY_TelPlane.at(plindx));
                b_ClusterSizeX_TelPlane_out.at(plindx) = m_OutTree->Branch(Form("ClusterSizeX_TelPlane%u", fileplanes.at(pls)), &m_ClusterSizeX_TelPlane.at(plindx));
                b_ClusterSizeY_TelPlane_out.at(plindx) = m_OutTree->Branch(Form("ClusterSizeY_TelPlane%u", fileplanes.at(pls)), &m_ClusterSizeY_TelPlane.at(plindx));
                b_ClusterNumPixels_TelPlane_out.at(plindx) = m_OutTree->Branch(Form("ClusterNumPixels_TelPlane%u", fileplanes.at(pls)), &m_ClusterNumPixels_TelPlane.at(plindx));
               }
           // Populate missing telescope references for this run
           nochan.clear();
           for (unsigned int k = 0; k < m_ntp_Ref.size(); k++)
               {
                found = false;
                for (unsigned int noref = 0; noref < fileRefs.size(); noref++) { if (m_ntp_Ref.at(k) == fileRefs.at(noref)) { found = true; break; } }
                if (!found) nochan.push_back(k);
               }
           if (nochan.size() > 0)
              {
               for (unsigned int w = 0; w < nochan.size(); w++)
                   {
                    b_InterceptX_Ref_out.at(nochan.at(w)) = m_OutTree->Branch(Form("InterceptX_Ref%u", m_ntp_Ref.at(nochan.at(w))), &no_ch_doub_vec);
                    b_InterceptY_Ref_out.at(nochan.at(w)) = m_OutTree->Branch(Form("InterceptY_Ref%u", m_ntp_Ref.at(nochan.at(w))), &no_ch_doub_vec);
                    b_InterceptResX_Ref_out.at(nochan.at(w)) = m_OutTree->Branch(Form("InterceptResX_Ref%u", m_ntp_Ref.at(nochan.at(w))), &no_ch_doub_vec);
                    b_InterceptResY_Ref_out.at(nochan.at(w)) = m_OutTree->Branch(Form("InterceptResY_Ref%u", m_ntp_Ref.at(nochan.at(w))), &no_ch_doub_vec);
                   }
              }
           // Populate missing planes for this run
           nochan.clear();
           for (unsigned int k = 0; k < m_ntp_planes.size(); k++)
               {
                found = false;
                for (unsigned int nopls = 0; nopls < fileplanes.size(); nopls++) { if (m_ntp_planes.at(k) == fileplanes.at(nopls)) { found = true; break; } }
                if (!found) nochan.push_back(k);
               }
           if (nochan.size() > 0)
              {
               for (unsigned int w = 0; w < nochan.size(); w++)
                   {
                    b_InterceptX_TelPlane_out.at(nochan.at(w)) = m_OutTree->Branch(Form("InterceptX_TelPlane%u", m_ntp_planes.at(nochan.at(w))), &no_ch_doub_vec);
                    b_InterceptY_TelPlane_out.at(nochan.at(w)) = m_OutTree->Branch(Form("InterceptY_TelPlane%u", m_ntp_planes.at(nochan.at(w))), &no_ch_doub_vec);
                    b_InterceptResX_TelPlane_out.at(nochan.at(w)) = m_OutTree->Branch(Form("InterceptResX_TelPlane%u", m_ntp_planes.at(nochan.at(w))), &no_ch_doub_vec);
                    b_InterceptResY_TelPlane_out.at(nochan.at(w)) = m_OutTree->Branch(Form("InterceptResY_TelPlane%u", m_ntp_planes.at(nochan.at(w))), &no_ch_doub_vec);
                    b_ClusterSizeX_TelPlane_out.at(nochan.at(w)) = m_OutTree->Branch(Form("ClusterSizeX_TelPlane%u", m_ntp_planes.at(nochan.at(w))), &no_ch_doub_vec);
                    b_ClusterSizeY_TelPlane_out.at(nochan.at(w)) = m_OutTree->Branch(Form("ClusterSizeY_TelPlane%u", m_ntp_planes.at(nochan.at(w))), &no_ch_doub_vec);
                    b_ClusterNumPixels_TelPlane_out.at(nochan.at(w)) = m_OutTree->Branch(Form("ClusterNumPixels_TelPlane%u", m_ntp_planes.at(nochan.at(w))), &no_ch_int_vec);
                   }
              }
           // Populate missing DUTs for this run
           nochan.clear();
           for (unsigned int k = 0; k < m_ntp_DUT.size(); k++)
               {
                found = false;
                for (unsigned int noDUTs = 0; noDUTs < fileDUTs.size(); noDUTs++) { if (m_ntp_DUT.at(k) == fileDUTs.at(noDUTs)) { found = true; break; } }
                if (!found) nochan.push_back(k);
               }
           if (nochan.size() > 0)
              {
               for (unsigned int w = 0; w < nochan.size(); w++)
                   {
                    b_InterceptX_DUT_out.at(nochan.at(w)) = m_OutTree->Branch(Form("InterceptX_DUT%u", m_ntp_DUT.at(nochan.at(w))), &no_ch_doub_vec);
                    b_InterceptY_DUT_out.at(nochan.at(w)) = m_OutTree->Branch(Form("InterceptY_DUT%u", m_ntp_DUT.at(nochan.at(w))), &no_ch_doub_vec);
                    b_InterceptResX_DUT_out.at(nochan.at(w)) = m_OutTree->Branch(Form("InterceptResX_DUT%u", m_ntp_DUT.at(nochan.at(w))), &no_ch_doub_vec);
                    b_InterceptResY_DUT_out.at(nochan.at(w)) = m_OutTree->Branch(Form("InterceptResY_DUT%u", m_ntp_DUT.at(nochan.at(w))), &no_ch_doub_vec);
                   }
              }
          }
       else if (m_RunBase->LGADBase::GetTrackPackage() == "EUTelescope")
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
       else if (m_RunBase->LGADBase::GetTrackPackage() == "TBA")
               { 
                //  TBA Branches should go here
               }
       else {
             std::cout << "    " << __FUNCTION__ << " ERROR: Unknown reconstruction framework!" << std::endl;
             return false;
            }
      }

   // Level 2 analysis requirements
   if ((m_RunBase->LGADBase::GetInRootFile(fileindx) != m_RunBase->LGADBase::GetOutRootFile() && filelvl > 1) || (m_RunBase->LGADBase::GetAnaStage() > 1 && filelvl <= 1))
      {
       // Chi2 for all fits
       b_chName_out = m_OutChTree->Branch("chName", &m_chName, "m_chName/I");
       b_MaxIndxFtChi2_out = m_OutChTree->Branch("MaxIndxFtChi2", &m_MaxIndxFtChi2, "m_MaxIndxFtChi2/D");
       b_MinIndxFtChi2_out = m_OutChTree->Branch("MinIndxFtChi2", &m_MinIndxFtChi2, "m_MinIndxFtChi2/D");
       b_MaxVoltFtChi2_out = m_OutChTree->Branch("MaxVoltFtChi2", &m_MaxVoltFtChi2, "m_MaxVoltFtChi2/D");
       b_MinVoltFtChi2_out = m_OutChTree->Branch("MinVoltFtChi2", &m_MinVoltFtChi2, "m_MinVoltFtChi2/D");
       b_StrIndxFtChi2_out = m_OutChTree->Branch("StrIndxFtChi2", &m_StrIndxFtChi2, "m_StrIndxFtChi2/D");
       b_EndIndxFtChi2_out = m_OutChTree->Branch("EndIndxFtChi2", &m_EndIndxFtChi2, "m_EndIndxFtChi2/D");
       b_NoiseFtChi2_out = m_OutChTree->Branch("NoiseFtChi2", &m_NoiseFtChi2, "m_NoiseFtChi2/D");
       b_NoiseErrFtChi2_out = m_OutChTree->Branch("NoiseErrFtChi2", &m_NoiseErrFtChi2, "m_NoiseErrFtChi2/D");
       b_PedestFtChi2_out = m_OutChTree->Branch("PedestFtChi2", &m_PedestFtChi2, "m_PedestFtChi2/D");
       b_PedestErrFtChi2_out = m_OutChTree->Branch("PedestErrFtChi", &m_PedestErrFtChi2, "m_PedestErrFtChi2/D");
       b_MaxTimeFtChi2_out = m_OutChTree->Branch("MaxTimeFtChi2", &m_MaxTimeFtChi2, "m_MaxTimeFtChi2/D");
       b_MinTimeFtChi2_out = m_OutChTree->Branch("MinTimeFtChi2", &m_MinTimeFtChi2, "m_MinTimeFtChi2/D");
       b_ChargeFtChi2_out = m_OutChTree->Branch("ChargeFtChi2", &m_ChargeFtChi2, "m_ChargeFtChi2/D");
       b_RiseTimeFtChi2_out = m_OutChTree->Branch("RiseTimeFtChi2", &m_RiseTimeFtChi2, "m_RiseTimeFtChi2/D");
       b_TriggTimeFtChi2_out = m_OutChTree->Branch("TriggTimeFtChi2", &m_TriggTimeFtChi2, "m_TriggTimeFtChi2/D");
       b_DVDTMaxFtChi2_out = m_OutChTree->Branch("DVDTMaxFtChi2", &m_DVDTMaxFtChi2, "m_DVDTMaxFtChi2/D");
       b_TriggToTFtChi2_out = m_OutChTree->Branch("TriggToTFtChi2", &m_TriggToTFtChi2, "m_TriggToTFtChi2/D");
       b_SignalFFTFtChi2_out = m_OutChTree->Branch("SignalFFTFtChi2", &m_SignalFFTFtChi2, "m_SignalFFTFtChi2/D");
       b_NoiseFFTFtChi2_out = m_OutChTree->Branch("NoiseFFTFtChi2", &m_NoiseFFTFtChi2, "m_NoiseFFTFtChi2/D");
       // Fit parameters for each channel
       b_ChFitMaxVolt_out = m_OutChTree->Branch("ChFitMaxVolt", &m_ChFitMxVlt.first, "m_ChFitMxVlt.first/D");
       b_ChFitMaxVoltErr_out = m_OutChTree->Branch("ChFitMaxVoltErr", &m_ChFitMxVlt.second, "m_ChFitMxVlt.second/D");
       b_ChFitMinVolt_out = m_OutChTree->Branch("ChFitMinVolt", &m_ChFitMinVlt.first, "m_ChFitMinVlt.first/D");
       b_ChFitMinVoltErr_out = m_OutChTree->Branch("ChFitMinVoltErr", &m_ChFitMinVlt.second, "m_ChFitMinVlt.second/D");
       b_ChFitMxIndx_out = m_OutChTree->Branch("ChFitMxIndx", &m_ChFitMxIndx.first, "m_ChFitMxIndx.first/I");
       b_ChFitMxIndxErr_out = m_OutChTree->Branch("ChFitMxIndxErr", &m_ChFitMxIndx.second, "m_ChFitMxIndx.second/I");
       b_ChFitMinIndx_out = m_OutChTree->Branch("ChFitMinIndx", &m_ChFitMinIndx.first, "m_ChFitMinIndx.first/I");
       b_ChFitMinIndxErr_out = m_OutChTree->Branch("ChFitMinIndxErr", &m_ChFitMinIndx.second, "m_ChFitMinIndx.second/I");
       b_ChFitStrIndx_out = m_OutChTree->Branch("ChFitStrIndx", &m_ChFitStrIndx.first, "m_ChFitStrIndx.first/I");
       b_ChFitStrIndxErr_out = m_OutChTree->Branch("ChFitStrIndxErr", &m_ChFitStrIndx.second, "m_ChFitStrIndx.second/I");
       b_ChFitEndIndx_out = m_OutChTree->Branch("ChFitEndIndx", &m_ChFitEndIndx.first, "m_ChFitEndIndx.first/I");
       b_ChFitEndIndxErr_out = m_OutChTree->Branch("ChFitEndIndxErr", &m_ChFitEndIndx.second, "m_ChFitEndIndx.second/I");
       b_ChFitIsSig_out = m_OutChTree->Branch("ChFitIsSig", &m_ChFitIsSig.first, "m_ChFitIsSig.first/D");
       b_ChFitIsSigErr_out = m_OutChTree->Branch("ChFitIsSigErr", &m_ChFitIsSig.second, "m_ChFitIsSig.second/D");
       b_ChFitInWin_out = m_OutChTree->Branch("ChFitInWin", &m_ChFitInWin.first, "m_ChFitInWin.first/D");
       b_ChFitInWinErr_out = m_OutChTree->Branch("ChFitInWinErr", &m_ChFitInWin.second, "m_ChFitInWin.second/D");
       b_ChFitVStr_out = m_OutChTree->Branch("ChFitVStr", &m_ChFitVStr.first, "m_ChFitVStr.first/D");
       b_ChFitVStrErr_out = m_OutChTree->Branch("ChFitVStrErr", &m_ChFitVStr.second, "m_ChFitVStr.second/D");
       b_ChFitNoise_out = m_OutChTree->Branch("ChFitNoise", &m_ChFitNoise.first, "m_ChFitNoise.first/D");
       b_ChFitNoiseUnc_out = m_OutChTree->Branch("ChFitNoiseUnc", &m_ChFitNoise.second, "m_ChFitNoise.second/D");
       b_ChFitPedestal_out = m_OutChTree->Branch("ChFitPedestal", &m_ChFitPedestal.first, "m_ChFitPedestal.first/D");
       b_ChFitPedestalUnc_out = m_OutChTree->Branch("ChFitPedestalUnc", &m_ChFitPedestal.second, "m_ChFitPedestal.second/D");
       b_ChFitNoiseSlp_out = m_OutChTree->Branch("ChFitNoiseSlp", &m_ChFitNoiseSlp.first, "m_ChFitNoiseSlp.first/D");
       b_ChFitNoiseSlpErr_out = m_OutChTree->Branch("ChFitNoiseSlpErr", &m_ChFitNoiseSlp.second, "m_ChFitNoiseSlp.second/D");
       b_ChFitPedestalSlp_out = m_OutChTree->Branch("ChFitPedestalSlp", &m_ChFitPedestalSlp.first, "m_ChFitPedestalSlp.first/D");
       b_ChFitPedestalSlpErr_out = m_OutChTree->Branch("ChFitPedestalSlpErr", &m_ChFitPedestalSlp.second, "m_ChFitPedestalSlp.second/D");
       b_ChFitNoiseErr_out = m_OutChTree->Branch("ChFitNoiseErr", &m_ChFitNoiseErr.first, "m_ChFitNoiseErr.first/D");
       b_ChFitNoiseErrUnc_out = m_OutChTree->Branch("ChFitNoiseErrUnc", &m_ChFitNoiseErr.second, "m_ChFitNoiseErr.second/D");
       b_ChFitPedestalErr_out = m_OutChTree->Branch("ChFitPedestalErr", &m_ChFitPedestalErr.first, "m_ChFitPedestalErr.first/D");
       b_ChFitPedestalErrUnc_out = m_OutChTree->Branch("ChFitPedestalErrUnc", &m_ChFitPedestalErr.second, "m_ChFitPedestalErr.second/D");
       b_ChFitMaxTime_out = m_OutChTree->Branch("ChFitMaxTime", &m_ChFitMaxTime.first, "m_ChFitMaxTime.first/D");
       b_ChFitMaxTimeErr_out = m_OutChTree->Branch("ChFitMaxTimeErr", &m_ChFitMaxTime.second, "m_ChFitMaxTime.second/D");
       b_ChFitMinTime_out = m_OutChTree->Branch("ChFitMinTime", &m_ChFitMinTime.first, "m_ChFitMinTime.first/D");
       b_ChFitMinTimeErr_out = m_OutChTree->Branch("ChFitMinTimeErr", &m_ChFitMinTime.second, "m_ChFitMinTime.second/D");
       b_ChFitCharge_out = m_OutChTree->Branch("ChFitCharge", &m_ChFitCharge.first, "m_ChFitCharge.first/D");
       b_ChFitChargeErr_out = m_OutChTree->Branch("ChFitChargeErr", &m_ChFitCharge.second, "m_ChFitCharge.second/D");
       b_ChFitRiseT_out = m_OutChTree->Branch("ChFitRiseT", &m_ChFitRiseT.first, "m_ChFitRiseT.first/D");
       b_ChFitRiseTErr_out = m_OutChTree->Branch("ChFitRiseTErr", &m_ChFitRiseT.second, "m_ChFitRiseT.second/D");
       b_ChFitJitRiseSNR_out = m_OutChTree->Branch("ChFitJitRiseSNR", &m_ChFitJitRiseSNR.first, "m_ChFitJitRiseSNR.first/D");
       b_ChFitJitRiseSNRErr_out = m_OutChTree->Branch("ChFitJitRiseSNRErr", &m_ChFitJitRiseSNR.second, "m_ChFitJitRiseSNR.second/D");
       b_ChFitTriggTime_out = m_OutChTree->Branch("ChFitTriggTime", &m_ChFitTriggTime.first, "m_ChFitTriggTime.first/D");
       b_ChFitTriggTimeErr_out = m_OutChTree->Branch("ChFitTriggTimeErr", &m_ChFitTriggTime.second, "m_ChFitTriggTime.second/D");
       b_ChFitDVDTMax_out = m_OutChTree->Branch("ChFitDVDTMax", &m_ChFitDVDTMax.first, "m_ChFitDVDTMax.first/D");
       b_ChFitDVDTMaxErr_out = m_OutChTree->Branch("ChFitDVDTMaxErr", &m_ChFitDVDTMax.second, "/m_ChFitDVDTMax.secondD");
       b_ChFitTriggToT_out = m_OutChTree->Branch("ChFitTriggToT", &m_ChFitTriggToT.first, "m_ChFitTriggToT.first/D");
       b_ChFitTriggToTErr_out = m_OutChTree->Branch("ChFitTriggToTErr", &m_ChFitTriggToT.second, "m_ChFitTriggToT.second/D");
       b_ChFitComplete_out = m_OutChTree->Branch("ChFitComplete", &m_ChFitComplete.first, "m_ChFitComplete.first/D");
       b_ChFitCompleteErr_out = m_OutChTree->Branch("ChFitCompleteErr", &m_ChFitComplete.second, "m_ChFitComplete.second/D");
       b_ChFitSoNR_out = m_OutChTree->Branch("ChFitSoNR", &m_ChFitSoNR.first, "m_ChFitSoNR.first/D");
       b_ChFitSoNRErr_out = m_OutChTree->Branch("ChFitSoNRErr", &m_ChFitSoNR.second, "m_ChFitSoNR.second/D");
       // FFT Second stage
       if (((FFTShape == 21 || FFTShape == 22) && m_RunBase->LGADBase::GetInRootFile(fileindx) != m_RunBase->LGADBase::GetOutRootFile()) ||
           (m_RunBase->LGADBase::GetAnaStage() > 1 && filelvl <= 1 && m_RunBase->LGADBase::GetDoFFT()))
          {
           b_ChFitSignalFFT_out = m_OutChTree->Branch("ChFitSignalFFT", &m_ChFitSignalFFT.first, "m_ChFitSignalFFT.first/D");
           b_ChFitSignalFFTErr_out = m_OutChTree->Branch("ChFitSignalFFTErr", &m_ChFitSignalFFT.second, "m_ChFitSignalFFT.second/D");
           b_ChFitNoiseFFT_out = m_OutChTree->Branch("ChFitNoiseFFT", &m_ChFitNoiseFFT.first, "m_ChFitNoiseFFT.first/D");
           b_ChFitNoiseFFT_out = m_OutChTree->Branch("ChFitNoiseFFTErr", &m_ChFitNoiseFFT.second, "m_ChFitNoiseFFT.second/D");
          }
       // CFD and DV/DT vectors
       b_ChFitCFDTime_out = m_OutChTree->Branch("ChFitCFDTime", &m_ChFitCFDTime);
       b_ChFitDVDTCFD_out = m_OutChTree->Branch("ChFitDVDTCFD", &m_ChFitDVDTCFD);
       b_ChFitCFDToT_out = m_OutChTree->Branch("ChFitCFDToT", &m_ChFitCFDToT);
       b_ChFitJitNdVdT_out = m_OutChTree->Branch("ChFitJitNdVdT", &m_ChFitJitNdVdT);
       b_CFDTimeFtChi2_out = m_OutChTree->Branch("CFDTimeFtChi2", &m_CFDTimeFtChi2);
       b_DVDTCFDFtChi2_out = m_OutChTree->Branch("DVDTCFDFtChi2", &m_DVDTCFDFtChi2);
       b_CFDToTFtChi2_out = m_OutChTree->Branch("CFDToTFtChi2", &m_CFDToTFtChi2);
       // Shape paramters for lvl2
       if (((FFTShape == 22 || FFTShape == 12) && m_RunBase->LGADBase::GetInRootFile(fileindx) != m_RunBase->LGADBase::GetOutRootFile()) ||
            (m_RunBase->LGADBase::GetAnaStage() > 1 && filelvl <= 1 && m_RunBase->LGADBase::GetWaveShape()))
          {
           b_WvSh_chName_out = m_OutShapeTree->Branch("ChShapeName", &m_WvSh_chName, "m_WvSh_chName/I");
           b_ChShapePolarity_out = m_OutShapeTree->Branch("ChShapePolarity", &m_ChShapePolarity, "m_ChShapePolarity/I");
           b_ChShapeMaxIndx_out = m_OutShapeTree->Branch("ChShapeMaxIndx", &m_ChShapeMaxIndx, "m_ChShapeMaxIndx/I");
           b_ChShapeMinIndx_out = m_OutShapeTree->Branch("ChShapeMinIndx", &m_ChShapeMinIndx, "m_ChShapeMinIndx/I");
           b_ChShapeStrIndx_out = m_OutShapeTree->Branch("ChShapeStrIndx", &m_ChShapeStrIndx, "m_ChShapeStrIndx/I");
           b_ChShapeEndIndx_out = m_OutShapeTree->Branch("ChShapeEndIndx", &m_ChShapeEndIndx, "m_ChShapeEndIndx/I");
           b_ChShapeMax_out = m_OutShapeTree->Branch("ChShapeMax", &m_ChShapeMax, "m_ChShapeMax/D");
           b_ChShapeMin_out = m_OutShapeTree->Branch("ChShapeMin", &m_ChShapeMin, "m_ChShapeMin/D");
           b_ChShapeMaxTime_out = m_OutShapeTree->Branch("ChShapeMaxTime", &m_ChShapeMaxTime, "m_ChShapeMaxTime/D");
           b_ChShapeMinTime_out = m_OutShapeTree->Branch("ChShapeMinTime", &m_ChShapeMinTime, "m_ChShapeMinTime/D");
           b_ChShapeCharge_out = m_OutShapeTree->Branch("ChShapeCharge", &m_ChShapeCharge, "m_ChShapeCharge/D");
           b_ChShapeRiseTime_out = m_OutShapeTree->Branch("ChShapeRiseTime", &m_ChShapeRiseTime, "m_ChShapeRiseTime/D");
           b_ChShapeJitRiseSNR_out = m_OutShapeTree->Branch("ChShapeJitRiseSNR", &m_ChShapeJitRiseSNR, "m_ChShapeJitRiseSNR/D");
           b_ChShapeTriggTime_out = m_OutShapeTree->Branch("ChShapeTriggTime", &m_ChShapeTriggTime, "m_ChShapeTriggTime/D");
           b_ChShapeDVDTMax_out = m_OutShapeTree->Branch("ChShapeDVDTMax", &m_ChShapeDVDTMax, "m_ChShapeDVDTMax/D");
           b_ChShapeTriggToT_out = m_OutShapeTree->Branch("ChShapeTriggToT", &m_ChShapeTriggToT, "m_ChShapeTriggToT/D");
           if ((FFTShape == 22 && m_RunBase->LGADBase::GetInRootFile(fileindx) != m_RunBase->LGADBase::GetOutRootFile()) ||
               (m_RunBase->LGADBase::GetAnaStage() > 1 && filelvl <= 1 && m_RunBase->LGADBase::GetWaveShape() && m_RunBase->LGADBase::GetDoFFT()))
              {
               b_ChShapeSignalFFT_out = m_OutShapeTree->Branch("ChShapeSignalFFTT", &m_ChShapeSignalFFT, "m_ChShapeSignalFFT/D");
               b_ChShapeNoiseFFT_out = m_OutShapeTree->Branch("ChShapeNoiseFFTT", &m_ChShapeNoiseFFT, "m_ChShapeNoiseFFT/D");
              }
           b_ChShapeNoise_out = m_OutShapeTree->Branch("ChShapeNoise", &m_ChShapeNoise.first, "m_ChShapeNoise.first/D");
           b_ChShapeNoiseErr_out = m_OutShapeTree->Branch("ChShapeNoiseErr", &m_ChShapeNoise.second, "m_ChShapeNoise.second/D");
           b_ChShapeSoN_out = m_OutShapeTree->Branch("ChShapeSoN", &m_ChShapeSoN.first, "m_ChShapeSoN.first/D");
           b_ChShapeSoNUnc_out = m_OutShapeTree->Branch("ChShapeSoNUnc", &m_ChShapeSoN.second, "m_ChShapeSoN.second/D");
           b_ChShapeCFDTime_out = m_OutShapeTree->Branch("ChShapeCFDTime", &m_ChShapeCFDTime);
           b_ChShapeDVDTCFD_out = m_OutShapeTree->Branch("ChShapeDVDTCFD", &m_ChShapeDVDTCFD);
           b_ChShapeCFDToT_out = m_OutShapeTree->Branch("ChShapeCFDToT", &m_ChShapeCFDToT);
           b_ChShapeJitNdVdT_out = m_OutShapeTree->Branch("ChShapeJitNdVdT", &m_ChShapeJitNdVdT);
           b_ChMeanVolt_out = m_OutShapeTree->Branch("ChMeanVolt", &m_ChMeanVolt);
          }
       // Third stage analysis rReuquirements
       if ((filelvl > 2 && m_RunBase->LGADBase::GetInRootFile(fileindx) != m_RunBase->LGADBase::GetOutRootFile()) ||
           (m_RunBase->LGADBase::GetAnaStage() > 2 && filelvl <= 2))
          {
           b_TmDiffChi2_out = m_OutTimeResTree->Branch("TmDiffChi2", &m_TmDiffChi2);
           b_TmDiffFitMean_out = m_OutTimeResTree->Branch("TmDiffFitMean", &m_TmDiffFitMean);
           b_TmDiffFitErr_out = m_OutTimeResTree->Branch("TmDiffFitErr", &m_TmDiffFitErr);
           b_CFD_ChA_out = m_OutTimeResTree->Branch("CFD_ChA", &m_CFD_ChA);
           b_ChA_out = m_OutTimeResTree->Branch("ChA", &m_ChA);
           b_ChB_out = m_OutTimeResTree->Branch("ChB", &m_ChB);
          }
       // 4th stage analysis reuqirements
       if (m_RunBase->LGADBase::GetAnaStage() > 3 && filelvl <= 3)
          {
          }
      }

   return true;
}
// --------------------------------------------------------------------------------------------------------------
bool LGADRun::AddOutputHistos()
{
   // Second stage analysis object initialization 
   if ((m_RunBase->LGADBase::GetFileLVL() < 2 && m_RunBase->LGADBase::GetAnaStage() > 1) || m_RunBase->LGADBase::GetAnaStage() == 2)
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
   if ((m_RunBase->LGADBase::GetFileLVL() < 3 && m_RunBase->LGADBase::GetAnaStage() > 2) || m_RunBase->LGADBase::GetAnaStage() == 3)
      {

      }

   return true;
}
// --------------------------------------------------------------------------------------------------------------
bool LGADRun::WriteBranches(unsigned int fileindx)
{
   // If we are copying to a new file just fill all the recovered and calculated info to the tree
   if (m_RunBase->LGADBase::GetInRootFile(fileindx) != m_RunBase->LGADBase::GetOutRootFile()) m_OutTree->Fill();
   else {
         if (m_RunBase->LGADBase::GetFileLVL() < 1 || m_RunBase->LGADBase::GetAnaStage() == 1)
            {
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
                  if (m_RunBase->LGADBase::GetDoFFT())
                     {
                      b_EvSignalFFT_out.at(ich)->Fill();
                      b_EvNoiseFFT_out.at(ich)->Fill();
                     }
                  b_EvComplete_out.at(ich)->Fill();
                  b_EvCFDTime_out.at(ich)->Fill();
                  b_EvVAdjCFD_out.at(ich)->Fill();
                  if (m_RunBase->LGADBase::GetVerbosity() >= 2) std::cout << __FUNCTION__ << " INFO: Written all brancehs: " << ich << "/" << m_nchan << "!" << std::endl;
                 }
            }
        }
   return true;
}