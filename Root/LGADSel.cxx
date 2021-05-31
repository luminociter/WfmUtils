// The class definition in LGADSel.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("MySelector.C")
// Root > T->Process("MySelector.C","some options")
// Root > T->Process("MySelector.C+")
//

#include "LGADUtils/LGADSel.h"
#include <vector> 

//#if !defined(__CINT__)
//ClassImp(LGADSel);
//#endif
//
//#if !defined(__CLING__)
//ClassImp(LGADSel);
//#endif
//
//#ifdef __CINT__
//#pragma link C++ class LGADSel;
//#endif
//
//#ifdef __ROOTCLING__
//#pragma link C++ class LGADSel;
//#endif

void LGADSel::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

void LGADSel::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t LGADSel::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}

Bool_t LGADSel::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either MySelector::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.


   return kTRUE;
}

void LGADSel::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void LGADSel::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

}

void LGADSel::Init(TTree *tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;
    fChain->SetMakeClass(1);

    // Count number of present channels
    m_nchan = 0;
    m_channels.clear();
    for (unsigned int ich = 1; ich < 65; ich++)
        {
         if (fChain->GetBranch(Form("t%02u", ich))) m_channels.push_back(ich);
         if (fChain->GetBranch(Form("ordrt%02u", ich))) m_instrument = Sampic;
        }
    m_nchan = m_channels.size();

    // Initialize Object Pointers
    t.reserve(m_nchan);
    w.reserve(m_nchan);
    b_t.resize(m_nchan);
    b_w.resize(m_nchan);
    b_ordrt.resize(m_nchan);
    b_physt.resize(m_nchan);
    b_vScale.resize(m_nchan);
    b_nPoints.resize(m_nchan);
    b_SnRate.resize(m_nchan);
    b_TriggTime.resize(m_nchan);
    for (unsigned int i = 0; i < m_nchan; i++)
        {
         t.push_back(new std::vector<double>(i));
         w.push_back(new std::vector<double>(i));
         t.at(i) = 0;
         w.at(i) = 0;
        }

    xTr = 0;
    yTr = 0;
    omegaTr = 0;
    phiTr = 0;
    kinkxTr = 0;
    kinkyTr = 0;
    chi2Tr = 0;
    plane_id = 0;
    track_id = 0;
    ndofTr = 0;
    dutID = 0;
    xHit = 0;
    yHit = 0;
    zHit = 0;

   fChain->SetBranchAddress("EvnNo", &EvnNo, &b_EvnNo);

   if (fChain->GetBranch("SnRate") && fChain->GetBranch("nPoints")) m_instrument = Sampic;
   else if (fChain->GetBranch("trigtime"))
           {
            m_instrument = Unasigned;
            fChain->SetBranchAddress("trigtime", &m_trigtime, &b_trigtime);
           }
   else m_instrument = InfiniiumScope;

   for (unsigned int ich = 0; ich < m_nchan; ich++)
       {
        if (fChain->GetBranch(Form("t%02u", m_channels.at(ich))))
           {
            fChain->SetBranchAddress(Form("t%02u", m_channels.at(ich)), &t.at(ich), &b_t.at(ich));
            fChain->SetBranchAddress(Form("w%02u", m_channels.at(ich)), &w.at(ich), &b_w.at(ich));
           }
       }

   LGADBase::SetVectorSize(m_nchan);

   if (m_instrument == Sampic)
      {
       fChain->SetBranchAddress("nPoints", &m_npoints.at(0), &b_nPoints.at(0));
       fChain->SetBranchAddress("SnRate", &m_srate.at(0), &b_SnRate.at(0));
       for (unsigned int ich = 0; ich < m_nchan; ich++)
           {
            if (fChain->GetBranch(Form("ordrt%02u", m_channels.at(ich))))
               {
                fChain->SetBranchAddress(Form("ordrt%02u", m_channels.at(ich)), &m_ordrt.at(ich), &b_ordrt.at(ich));
                fChain->SetBranchAddress(Form("physt%02u", m_channels.at(ich)), &m_physt.at(ich), &b_physt.at(ich));
               }
           }
      }
   else {
         for (unsigned int ich = 0; ich < m_nchan; ich++)
             {
              if (fChain->GetBranch(Form("nPoints%02u", m_channels.at(ich))))
                 { 
                  fChain->SetBranchAddress(Form("nPoints%02u", m_channels.at(ich)), &m_npoints.at(ich), &b_nPoints.at(ich));
                  fChain->SetBranchAddress(Form("SnRate%02u", m_channels.at(ich)), &m_srate.at(ich), &b_SnRate.at(ich));
                  fChain->SetBranchAddress(Form("vScale%02u", m_channels.at(ich)), &m_scale.at(ich), &b_vScale.at(ich));
                  if (fChain->GetBranch(Form("TriggTime%02u", m_channels.at(ich))))
                     { 
                      m_instrument = TektronixScope;
                      fChain->SetBranchAddress(Form("TriggTime%02u", m_channels.at(ich)), &m_triggTime.at(ich), &b_TriggTime.at(ich));
                     }
                  else if (m_instrument != InfiniiumScope) m_instrument = LabTXT;
                 }
             }
        }

   if (fChain->GetBranch("plane_id"))
      {
       LGADBase::SetTrackComb(true);
       fChain->SetBranchAddress("xTr", &xTr, &b_xTr);
       fChain->SetBranchAddress("yTr", &yTr, &b_yTr);
       fChain->SetBranchAddress("omegaTr", &omegaTr, &b_omegaTr);
       fChain->SetBranchAddress("phiTr", &phiTr, &b_phiTr);
       fChain->SetBranchAddress("kinkxTr", &kinkxTr, &b_kinkxTr);
       fChain->SetBranchAddress("kinkyTr", &kinkyTr, &b_kinkyTr);
       fChain->SetBranchAddress("chi2Tr", &chi2Tr, &b_chi2Tr);
       fChain->SetBranchAddress("plane_id", &plane_id, &b_plane_id);
       fChain->SetBranchAddress("track_id", &track_id, &b_track_id);
       fChain->SetBranchAddress("ndofTr", &ndofTr, &b_ndofTr);
       fChain->SetBranchAddress("trigger_id", &trigger_id, &b_trigger_id);
       fChain->SetBranchAddress("timestampTr", &timestampTr, &b_timestampTr);
      }	  
   if (fChain->GetBranch("dutID"))
      {
       LGADBase::SetFEi4Eff(true);
       fChain->SetBranchAddress("xHit", &xHit, &b_xHit);
       fChain->SetBranchAddress("yHit", &yHit, &b_yHit);
       fChain->SetBranchAddress("zHit", &zHit, &b_zHit);
       fChain->SetBranchAddress("dutID", &dutID, &b_dutID);
      }

}
