/*
*               TrackCombine.cxx
*
*
*      Author: Gkougkousis Evangelos - Leonidas
*              egkougko@cern.ch - 2022
*                    CERN-GENEVA
*/

#include "../LGADUtils/LGADBase.h"

void LGADBase::SetTrackInDataNames(TString DataDir, TString DataName)
{
    // Test if the directory contains extention
    m_Trackdatadir.Clear();
    m_Trackdataname.clear();
    if (DataDir.Length() > 4 && DataDir.Contains(".root"))
       {
        m_Trackdatadir = DataDir(0, DataDir.Length() - 5);  // remove extention
        m_Trackdatadir = DataDir(0, DataDir.Last('/') + 1); // get directory
        m_Trackdataname.push_back(DataDir(DataDir.Last('/') + 1, DataDir.Length() - DataDir.Last('/') - 6)); // get file
       }
    else {
          TString sLast = DataDir[DataDir.Length() - 1];
          if (!(sLast.Contains("/"))) DataDir = DataDir + "/"; // get directory
          m_Trackdatadir = DataDir;
         }

    // Have directory, need to get the filename or list of filenames if not already done
    if (DataName != "")
       {
        if (DataName.Contains(".root")) m_Trackdataname.push_back(DataName(0, DataName.Length() - 5));
        else m_Trackdataname.push_back(DataName);
       }
    else { // dataname is not set but the file was not in the name of the folder, so searching for list
          if (m_Trackdataname.size() == 0)
             {
              std::cout << __FUNCTION__ << " Scanning for tracking files in folder " << m_Trackdatadir << " ...";
              m_Trackdataname.resize(CountFiles(m_Trackdatadir, "root"));
              std::cout << m_Trackdataname.size() << " found!" << std::endl;
              std::vector <string> temp = ListFileNames(m_Trackdatadir, "root");
              for (unsigned int i = 0; i < temp.size(); i++) m_Trackdataname.at(i) = (TString)((temp.at(i)).erase((temp.at(i)).size() - 5));
             }
         }
}
// --------------------------------------------------------------------------------------------------------------
bool LGADBase::OpenTrackFile(TString trackdirdir, std::vector<TString> trackdataname, std::vector<unsigned int>& filesIndx)
{
    if (trackdirdir.IsNull())
       { 
        std::cout << "    " << __FUNCTION__ << " ERROR: Input tracking data directory not set! Abording..." << std::endl;
        return false; 
       }
    if (trackdataname.size() == 0)
       {
        std::cout << "    " << __FUNCTION__ << " ERROR: No tracking files set in directory " << trackdirdir << "!Abording..." << std::endl;
        return false;
       }

    m_trackfiles.clear();
    m_nplanes.clear();
    m_planes.clear();
    m_nRef.clear();
    m_Refs.clear();
    m_nDUT.clear();
    m_DUTs.clear();
    unsigned int files = 0;
    bool found = false;
    unsigned int binaryindx = 0; 

    // Figure out the number of files and the track package to use
    for (unsigned int i = 0; i < trackdataname.size(); i++)
        {
         // Check if the index is within the vector of indexes proided by the binary data conversion funtion or if the run is not to be consideed
         found = false;
         for (unsigned int as = 0; as < filesIndx.size(); as++) { if (filesIndx.at(as) == i) { found = true; binaryindx = as; break; } }
         if (!found) continue;
         m_trackfiles.push_back(TFile::Open(trackdirdir+trackdataname.at(i)+".root"));
         std::cout << "    " << __FUNCTION__ << " INFO: Opening tracking file " << (TString)(trackdirdir + trackdataname.at(i) + ".root") << "..." << std::endl;
         files++;
         if (m_trackfiles.at(files-1))
            { 
             // ##################################### CoryRecon integration part #####################################
             if ((TTree*)(m_trackfiles.at(files-1))->Get("TreeWriterPassive/USBPIX_GEN2_BOARD_131_20/tracks"))
                {     
                 m_trackTrees.push_back(NULL);
                 m_trackfiles.at(files - 1)->GetObject("TreeWriterPassive/USBPIX_GEN2_BOARD_131_20/tracks", m_trackTrees.back());            
                 if (GetTrackPackage() == "") SetTrackPackage("Cory");
                 else if (GetTrackPackage() != "Cory") 
                         {
                          std::cout << "    " << __FUNCTION__ << " WARNING: Trying to mix different reco framework files! File: " << trackdirdir + trackdataname.at(i) + ".root" << " skipped!" << std::endl;
                          std::cout << "    " << "               Current file is from Cory while previous runs are from: " << GetTrackPackage() << "!" << std::endl;
                          files--;
                          m_trackTrees.pop_back();
                          m_trackfiles.back()->Close();
                          m_trackfiles.pop_back();
                          filesIndx.erase(filesIndx.begin() + binaryindx);
                          continue;
                         }
                 // Get Telescope planes trees
                 m_nplanes.push_back(0);
                 for (unsigned int k = 0; k < 8; k++)
                     { 
                      if ((TTree*)(m_trackfiles.at(files-1))->Get(Form("TreeWriterPassive/MIMOSA26_%01u/tracks", k)))
                         {
                          m_trackTrees.push_back(NULL);
                          m_trackfiles.at(files-1)->GetObject(Form("TreeWriterPassive/MIMOSA26_%01u/tracks", k), m_trackTrees.back());
                          m_nplanes.back() += 1;
                          m_planes.push_back(k);
                         }
                     }
                 // Reference trees
                 m_nRef.push_back(0);
                 for (unsigned int k = 0; k < 3; k++)
                     {
                      std::string treedir;
                      if (k == 0) treedir = "TreeWriterPassive/LGAD_F/tracks";
                      else if (k == 1) treedir = "TreeWriterPassive/LGAD_B/tracks";
                      else if (k == 2) treedir = "TreeWriterPassive/LGAD_M/tracks";
                      if ((TTree*)(m_trackfiles.at(files - 1))->Get(treedir.c_str()))
                         {
                          m_trackTrees.push_back(NULL);
                          m_trackfiles.at(files - 1)->GetObject(treedir.c_str(), m_trackTrees.back());
                          m_nRef.back() += 1;
                          m_Refs.push_back(k);
                         }
                     }
                 // Get DUT trees
                 m_nDUT.push_back(0);
                 for (unsigned int k = 0; k < 8; k++)
                     {                      
                      if ((TTree*)(m_trackfiles.at(files - 1))->Get(Form("TreeWriterPassive/Plane_%01u/tracks", k)))
                         { 
                          m_trackTrees.push_back(NULL);
                          m_trackfiles.at(files - 1)->GetObject(Form("TreeWriterPassive/Plane_%01u/tracks", k), m_trackTrees.back());
                          m_nDUT.back() += 1;
                          m_DUTs.push_back(k);
                         }
                     }
                }
             // ##################################### EUTelescope integration part #####################################
             else if ((TTree*)(m_trackfiles.at(files - 1))->Get("ZeroSupressed"))
                     {
                      m_trackTrees.push_back(NULL);
                      m_trackfiles.at(files - 1)->GetObject("Hits", m_trackTrees.back());
                      if (GetTrackPackage() == "") SetTrackPackage("EUTelescope");
                      else if (GetTrackPackage() != "EUTelescope") 
                              {
                               std::cout << "    " << __FUNCTION__ << " WARNING: Trying to mix different reco framework files! File: " << trackdirdir + trackdataname.at(i) + ".root" << " skipped!" << std::endl;
                               std::cout << "    " << "               Current file is from EUTelescope while previous runs are from: " << GetTrackPackage() << "!" << std::endl;
                               files--;
                               filesIndx.erase(filesIndx.begin() + binaryindx);
                               m_trackTrees.pop_back();
                               m_trackfiles.back()->Close();
                               m_trackfiles.pop_back();
                               continue;
                              }
                      // Get the Tracks tree
                      if (!((TTree*)(m_trackfiles.at(files - 1))->Get("Tracks")))
                         {
                          std::cout << "    " << __FUNCTION__  << " WARNING: Could not read tracks tree for EUTelescope file " << trackdirdir + trackdataname.at(i) + ".root" << ". Skipping..." << std::endl;
                          files--;
                          filesIndx.erase(filesIndx.begin() + binaryindx);
                          m_trackTrees.pop_back();
                          m_trackfiles.back()->Close();
                          m_trackfiles.pop_back();
                          continue;
                         }
                      m_trackTrees.push_back(NULL);
                      m_trackfiles.at(files - 1)->GetObject("Tracks", m_trackTrees.back());
                     }
             // ##################################### TBA integration part #####################################
             else if ((TTree*)(m_trackfiles.at(files-1))->Get(""))
                     {
                      m_trackTrees.push_back(NULL);
                      m_trackfiles.at(files-1)->GetObject("", m_trackTrees.back());
                      if (GetTrackPackage() == "") SetTrackPackage("TBA");
                      else if (GetTrackPackage() != "TBA") 
                              {
                               std::cout << "    " << __FUNCTION__ << " WARNING: Trying to mix different reco framework files! File: " << trackdirdir + trackdataname.at(i) + ".root" << " skipped!" << std::endl;
                               std::cout << "    " << "               Current file is from TBA framework while previous runs are from: " << GetTrackPackage() << "!" << std::endl;
                               files--;
                               filesIndx.erase(filesIndx.begin() + binaryindx);
                               m_trackTrees.pop_back();
                               m_trackfiles.back()->Close();
                               m_trackfiles.pop_back();
                               continue;
                              }
                      if (files == 1) // Pointer vector initialization
                         {
                            //
                            // Missingn code for TBA from Bonn, use example fiel to implement
                            // 
                         }
                     }
             // ##################################### Error part #####################################
             else {                  
                   std::cout << "    " << __FUNCTION__ << " WARNING: Unrecognised file structure for " << trackdirdir + trackdataname.at(i) + ".root"
                             << "! Unknown reconstruction framework." << std::endl;
                   files--;
                   filesIndx.erase(filesIndx.begin() + binaryindx);
                   m_trackfiles.back()->Close();
                   m_trackfiles.pop_back();
                   continue;
                  }
            }
         else {
               std::cout << "    " << __FUNCTION__ << " WARNING: Failed oppening tracking file " << trackdirdir + trackdataname.at(i) + ".root" << "! Probably Zombie?" << std::endl;
               files--;
               filesIndx.erase(filesIndx.begin() + binaryindx);
               m_trackfiles.back()->Close();
               m_trackfiles.pop_back();
               continue;
              }
        }

    // Clean and asign the appropriate number of pointers for the branches
    int vectorindx = 0;
    unsigned int NofPlanes = 0;
    unsigned int NofRefs = 0;
    unsigned int NofDUTs = 0;
    if (GetTrackPackage() == "Cory")
       {
        // Initialize pointer vectors
        b_m_ROI.clear();
        m_ROI_Dbl.clear();
        m_ROI_Int.clear();
        b_m_Ref.clear();
        m_Ref_Dbl.clear();
        m_Ref_Int.clear();
        b_m_DUT.clear();
        m_DUT_Dbl.clear();
        m_DUT_Int.clear();
        b_m_TelPlane.clear();
        m_TelPlane_Dbl.clear();
        m_TelPlane_Int.clear();
        m_ROI_Dbl.reserve(8 * files);
        m_ROI_Int.reserve(3 * files);
        m_TelPlane_Dbl.reserve(6*m_planes.size());
        m_TelPlane_Int.reserve(2 * m_planes.size());
        m_Ref_Dbl.reserve(4 * m_Refs.size());
        m_Ref_Int.reserve(m_Refs.size());
        m_DUT_Dbl.reserve(4 * m_Refs.size());
        m_DUT_Int.reserve(m_DUTs.size());
        for (unsigned int f = 0; f < files; f++)
            {
             if (f > 0) 
                {
                 NofPlanes += m_nplanes.at(f-1);
                 NofRefs += m_nRef.at(f-1);
                 NofDUTs+= m_nDUT.at(f-1);
                }
             // ROI Pointers
             for (unsigned w = 0; w < 11; w++) b_m_ROI.push_back(NULL);
             for (unsigned l = 0; l < 8; l++) m_ROI_Dbl.push_back(NULL);
             for (unsigned m = 0; m < 3; m++) m_ROI_Int.push_back(NULL);
             // Telescope plane pointers
             for (unsigned int k = 0; k < m_nplanes.at(f); k++)
                 {
                  for (unsigned w = 0; w < 8; w++) b_m_TelPlane.push_back(NULL);   //!
                  for (unsigned l = 0; l < 6; l++) m_TelPlane_Dbl.push_back(NULL);
                  for (unsigned m = 0; m < 2; m++) m_TelPlane_Int.push_back(NULL);
                 }
             // Reference Pointers
             for (unsigned int k = 0; k < m_nRef.at(f); k++)
                 {
                  for (unsigned w = 0; w < 5; w++) b_m_Ref.push_back(NULL);   //!
                  for (unsigned l = 0; l < 4; l++) m_Ref_Dbl.push_back(NULL);
                  m_Ref_Int.push_back(NULL);
                 }
             // Initialise the DUT Pointers
             for (unsigned int k = 0; k < m_nDUT.at(f); k++)
                 {
                  for (unsigned w = 0; w < 5; w++) b_m_DUT.push_back(NULL);   //!  
                  for (unsigned l = 0; l < 4; l++) m_DUT_Dbl.push_back(NULL);
                  m_DUT_Int.push_back(NULL);
                 }
            }
       }
    else if (GetTrackPackage() == "EUTelescope")
            {
             // Initialize pointer vectors
             b_m_Hits.clear();
             m_HitsInt.clear();
             b_m_Tracks.clear();
             m_TracksInt.clear();
             m_TracksDbl.clear();
             m_TrsInt.clear();
             m_HitsInt.reserve(4*files);
             m_TracksInt.reserve(3*files);
             m_TracksDbl.reserve(7*files);
             m_TrsInt.reserve(3*files);
             for (unsigned int f = 0; f < files; f++)
                 {
                  for (unsigned int w = 0; w < 4; w++) 
                      {
                       b_m_Hits.push_back(NULL);   //!
                       m_HitsInt.push_back(NULL);
                      }
                  for (unsigned int l = 0; l < 13; l++) b_m_Tracks.push_back(NULL);   //!
                  for (unsigned int m = 0; m < 3; m++) 
                      {
                       m_TracksInt.push_back(NULL);
                       m_TrsInt.push_back(0);
                      }
                  for (unsigned int n = 0; n < 7; n++) m_TracksDbl.push_back(NULL);
                 }
            }
    else if (GetTrackPackage() == "TBA")
            {
               //
               // Missingn code for TBA from Bonn, use example fiel to implement
               // 
            }

    int treeindx = 0;
    int brchindx = 0;
    NofPlanes = 0;
    NofRefs = 0;
    NofDUTs = 0;
    for (unsigned int f = 0; f < files; f++)
        {
         if (GetTrackPackage() == "Cory")
            {
             if (f > 0) 
                {
                 NofPlanes += m_nplanes.at(f-1);
                 NofRefs += m_nRef.at(f-1);
                 NofDUTs+= m_nDUT.at(f-1);
                }
             // Get ROI tree
             treeindx = NofPlanes+NofRefs+NofDUTs+f;
             brchindx = 11*f;
             vectorindx = 8*f;
             m_trackTrees.at(treeindx)->SetBranchAddress("InterceptX", &m_ROI_Dbl.at(vectorindx), &b_m_ROI.at(brchindx));
             m_trackTrees.at(treeindx)->SetBranchAddress("InterceptY", &m_ROI_Dbl.at(vectorindx + 1), &b_m_ROI.at(brchindx + 1));
             m_trackTrees.at(treeindx)->SetBranchAddress("InterceptResolutionX", &m_ROI_Dbl.at(vectorindx + 2), &b_m_ROI.at(brchindx + 2));
             m_trackTrees.at(treeindx)->SetBranchAddress("InterceptResolutionY", &m_ROI_Dbl.at(vectorindx + 3), &b_m_ROI.at(brchindx + 3));
             m_trackTrees.at(treeindx)->SetBranchAddress("clusterSizeX", &m_ROI_Dbl.at(vectorindx + 4), &b_m_ROI.at(brchindx + 4));
             m_trackTrees.at(treeindx)->SetBranchAddress("clusterSizeY", &m_ROI_Dbl.at(vectorindx + 5), &b_m_ROI.at(brchindx + 5));
             m_trackTrees.at(treeindx)->SetBranchAddress("pixelX", &m_ROI_Dbl.at(vectorindx + 6), &b_m_ROI.at(brchindx + 6));
             m_trackTrees.at(treeindx)->SetBranchAddress("pixelY", &m_ROI_Dbl.at(vectorindx + 7), &b_m_ROI.at(brchindx + 7));
             vectorindx = 3*f;
             m_trackTrees.at(treeindx)->SetBranchAddress("clusterNumPixels", &m_ROI_Int.at(vectorindx), &b_m_ROI.at(brchindx + 8));
             m_trackTrees.at(treeindx)->SetBranchAddress("HasAsso", &m_ROI_Int.at(vectorindx + 1), &b_m_ROI.at(brchindx + 9));
             m_trackTrees.at(treeindx)->SetBranchAddress("EventID", &m_ROI_Int.at(vectorindx + 2), &b_m_ROI.at(brchindx + 10));
             PrintBranches(m_trackTrees.at(treeindx), "ROI");
             // Get Telescope planes trees
             for (unsigned int k = 0; k < m_nplanes.at(f); k++)
                 {
                  treeindx = NofPlanes+NofRefs+NofDUTs+f+1+k;
                  brchindx = 8*(NofPlanes+k);
                  vectorindx = 6*(NofPlanes+k);
                  m_trackTrees.at(treeindx)->SetBranchAddress("InterceptX", &m_TelPlane_Dbl.at(vectorindx), &b_m_TelPlane.at(brchindx));
                  m_trackTrees.at(treeindx)->SetBranchAddress("InterceptY", &m_TelPlane_Dbl.at(vectorindx + 1), &b_m_TelPlane.at(brchindx + 1));
                  m_trackTrees.at(treeindx)->SetBranchAddress("InterceptResolutionX", &m_TelPlane_Dbl.at(vectorindx + 2), &b_m_TelPlane.at(brchindx + 2));
                  m_trackTrees.at(treeindx)->SetBranchAddress("InterceptResolutionY", &m_TelPlane_Dbl.at(vectorindx + 3), &b_m_TelPlane.at(brchindx + 3));
                  m_trackTrees.at(treeindx)->SetBranchAddress("clusterSizeX", &m_TelPlane_Dbl.at(vectorindx + 4), &b_m_TelPlane.at(brchindx + 4));
                  m_trackTrees.at(treeindx)->SetBranchAddress("clusterSizeY", &m_TelPlane_Dbl.at(vectorindx + 5), &b_m_TelPlane.at(brchindx + 5));
                  vectorindx = 2*(NofPlanes+k);
                  m_trackTrees.at(treeindx)->SetBranchAddress("clusterNumPixels", &m_TelPlane_Int.at(vectorindx), &b_m_TelPlane.at(brchindx+6));
                  m_trackTrees.at(treeindx)->SetBranchAddress("EventID", &m_TelPlane_Int.at(vectorindx+1), &b_m_TelPlane.at(brchindx+7));
                  PrintBranches(m_trackTrees.at(treeindx), Form("TelescopePlane_%01u", k));
                 }
             // Reference trees
             m_nRef.push_back(0);
             for (unsigned int k = 0; k < m_nRef.at(f); k++)
                 {
                  treeindx = NofPlanes+NofRefs+NofDUTs+f+1+m_nplanes.at(f)+k;
                  brchindx = 5*(NofRefs+k);
                  vectorindx = 4*(NofRefs+k);
                  m_trackTrees.at(treeindx)->SetBranchAddress("InterceptX", &m_Ref_Dbl.at(vectorindx), &b_m_Ref.at(brchindx));
                  m_trackTrees.at(treeindx)->SetBranchAddress("InterceptY", &m_Ref_Dbl.at(vectorindx + 1), &b_m_Ref.at(brchindx + 1));
                  m_trackTrees.at(treeindx)->SetBranchAddress("InterceptResolutionX", &m_Ref_Dbl.at(vectorindx + 2), &b_m_Ref.at(brchindx + 2));
                  m_trackTrees.at(treeindx)->SetBranchAddress("InterceptResolutionY", &m_Ref_Dbl.at(vectorindx + 3), &b_m_Ref.at(brchindx + 3));
                  m_trackTrees.at(treeindx)->SetBranchAddress("EventID", &m_Ref_Int.at(NofRefs+k), &b_m_Ref.at(brchindx + 4));
                  PrintBranches(m_trackTrees.at(treeindx), Form("Reference_%01u", k));
                 }
             // Get DUT trees
             m_nDUT.push_back(0);
             for (unsigned int k = 0; k < m_nDUT.at(f); k++)
                 {
                  treeindx = NofPlanes+NofRefs+NofDUTs+f+1+m_nplanes.at(f)+m_nRef.at(f)+k;
                  brchindx = 5*(NofDUTs+k);
                  vectorindx = 4*(NofDUTs+k);
                  m_trackTrees.at(treeindx)->SetBranchAddress("InterceptX", &m_DUT_Dbl.at(vectorindx), &b_m_DUT.at(brchindx));
                  m_trackTrees.at(treeindx)->SetBranchAddress("InterceptY", &m_DUT_Dbl.at(vectorindx + 1), &b_m_DUT.at(brchindx + 1));
                  m_trackTrees.at(treeindx)->SetBranchAddress("InterceptResolutionX", &m_DUT_Dbl.at(vectorindx + 2), &b_m_DUT.at(brchindx + 2));
                  m_trackTrees.at(treeindx)->SetBranchAddress("InterceptResolutionY", &m_DUT_Dbl.at(vectorindx + 3), &b_m_DUT.at(brchindx + 3));
                  m_trackTrees.at(treeindx)->SetBranchAddress("EventID", &m_DUT_Int.at(NofDUTs+k), &b_m_DUT.at(brchindx+4));
                  PrintBranches(m_trackTrees.at(treeindx), Form("DUT_%01u", k));
                 }
            }
         else if (GetTrackPackage() == "EUTelescope")
                 {
                  // Assign pointers to tree branches for Hits tree
                  vectorindx = 4*f;
                  treeindx = 2*f;
                  brchindx = 4*f;
                  m_trackTrees.at(treeindx)->SetBranchAddress("ID", &m_HitsInt.at(vectorindx), &b_m_Hits.at(brchindx));
                  m_trackTrees.at(treeindx)->SetBranchAddress("xPosH", &m_HitsInt.at(vectorindx+1), &b_m_Hits.at(brchindx+1));
                  m_trackTrees.at(treeindx)->SetBranchAddress("yPosH", &m_HitsInt.at(vectorindx+2), &b_m_Hits.at(brchindx+2));
                  m_trackTrees.at(treeindx)->SetBranchAddress("zPosH", &m_HitsInt.at(vectorindx+3), &b_m_Hits.at(brchindx+3));
                  PrintBranches(m_trackTrees.at(treeindx), "Hits");
                  vectorindx = 7*f;
                  treeindx = 2*f+1;
                  brchindx = 13*f;
                  // Assign pointers to tree branches for the Tracks tree
                  m_trackTrees.at(treeindx)->SetBranchAddress("xPos", &m_TracksDbl.at(vectorindx), &b_m_Tracks.at(brchindx));
                  m_trackTrees.at(treeindx)->SetBranchAddress("yPos", &m_TracksDbl.at(vectorindx+1), &b_m_Tracks.at(brchindx+1));
                  m_trackTrees.at(treeindx)->SetBranchAddress("omega", &m_TracksDbl.at(vectorindx+2), &b_m_Tracks.at(brchindx+2));
                  m_trackTrees.at(treeindx)->SetBranchAddress("phi", &m_TracksDbl.at(vectorindx+3), &b_m_Tracks.at(brchindx+3));
                  m_trackTrees.at(treeindx)->SetBranchAddress("kinkx", &m_TracksDbl.at(vectorindx+4), &b_m_Tracks.at(brchindx+4));
                  m_trackTrees.at(treeindx)->SetBranchAddress("kinky", &m_TracksDbl.at(vectorindx+5), &b_m_Tracks.at(brchindx+5));
                  m_trackTrees.at(treeindx)->SetBranchAddress("chi2", &m_TracksDbl.at(vectorindx+6), &b_m_Tracks.at(brchindx+6));
                  vectorindx = 3*f;
                  m_trackTrees.at(treeindx)->SetBranchAddress("planeID", &m_TracksInt.at(vectorindx), &b_m_Tracks.at(brchindx+7));
                  m_trackTrees.at(treeindx)->SetBranchAddress("trackID", &m_TracksInt.at(vectorindx+1), &b_m_Tracks.at(brchindx+8));
                  m_trackTrees.at(treeindx)->SetBranchAddress("ndof", &m_TracksInt.at(vectorindx+2), &b_m_Tracks.at(brchindx+9));
                  m_trackTrees.at(treeindx)->SetBranchAddress("eventNumber", &m_TrsInt.at(vectorindx), &b_m_Tracks.at(brchindx+10));
                  m_trackTrees.at(treeindx)->SetBranchAddress("triggerID", &m_TrsInt.at(vectorindx+1), &b_m_Tracks.at(brchindx+11));
                  m_trackTrees.at(treeindx)->SetBranchAddress("timestamp", &m_TrsInt.at(vectorindx+2), &b_m_Tracks.at(brchindx+12));
                  PrintBranches(m_trackTrees.at(treeindx), "Tracks");
                 }
         else if (GetTrackPackage() == "TBA")
                 {
                   //
                   // Missingn code for TBA from Bonn, use example fiel to implement
                   //    
                 }
        }

    std::cout << std::endl;
    if (m_trackTrees.size() == 0) return false;
    else return true;
}
// --------------------------------------------------------------------------------------------------------------
void LGADBase::InitCoryPointers(unsigned int nRef, unsigned int nDUT, unsigned int nplanes)
{
    // Pointer Initialization
    m_ClusterSizeX_ROI = NULL;
    m_ClusterSizeY_ROI = NULL;
    m_PixelX_ROI = NULL;
    m_PixelY_ROI = NULL;
    m_ClusterNumPixels_ROI = NULL;
    m_InterceptX_ROI = NULL;
    m_InterceptY_ROI = NULL;
    m_InterceptResX_ROI = NULL;
    m_InterceptResY_ROI = NULL;
    m_HasAsso = NULL;

    // Initialization of reference pointers
    m_InterceptX_Ref.clear();
    m_InterceptY_Ref.clear();
    m_InterceptResX_Ref.clear();
    m_InterceptResY_Ref.clear();
    m_InterceptX_Ref.resize(nRef, nullptr);
    m_InterceptY_Ref.resize(nRef, nullptr);
    m_InterceptResX_Ref.resize(nRef, nullptr);
    m_InterceptResY_Ref.resize(nRef, nullptr);

    // Initialization of DUT pointers
    m_InterceptX_DUT.clear();
    m_InterceptY_DUT.clear();
    m_InterceptResX_DUT.clear();
    m_InterceptResY_DUT.clear();
    m_InterceptX_DUT.resize(nDUT, nullptr);
    m_InterceptY_DUT.resize(nDUT, nullptr);
    m_InterceptResX_DUT.resize(nDUT, nullptr);
    m_InterceptResY_DUT.resize(nDUT, nullptr);

    // Initialization of Telescope plane pointers
    m_InterceptX_TelPlane.clear();
    m_InterceptY_TelPlane.clear();
    m_InterceptResX_TelPlane.clear();
    m_InterceptResY_TelPlane.clear();
    m_ClusterSizeX_TelPlane.clear();
    m_ClusterSizeY_TelPlane.clear();
    m_ClusterNumPixels_TelPlane.clear();
    m_InterceptX_TelPlane.resize(nplanes, nullptr);
    m_InterceptY_TelPlane.resize(nplanes, nullptr);
    m_InterceptResX_TelPlane.resize(nplanes, nullptr);
    m_InterceptResY_TelPlane.resize(nplanes, nullptr);
    m_ClusterSizeX_TelPlane.resize(nplanes, nullptr);
    m_ClusterSizeY_TelPlane.resize(nplanes, nullptr);
    m_ClusterNumPixels_TelPlane.resize(nplanes, nullptr);
}
// --------------------------------------------------------------------------------------------------------------
void LGADBase::InitEUTelPointers()
{
    // Pointer Initialization
    m_xTr = NULL;
    m_yTr = NULL;
    m_omegaTr = NULL;
    m_phiTr = NULL;
    m_kinkxTr = NULL;
    m_kinkyTr = NULL;
    m_chi2Tr = NULL;
    m_plane_id = NULL;
    m_track_id = NULL;
    m_ndofTr = NULL;
    m_ID = NULL;
    m_xHit = NULL;
    m_yHit = NULL;
    m_zHit = NULL;
}
// --------------------------------------------------------------------------------------------------------------
void LGADBase::InitTBAPointers()
{
     //
     // Missingn code for TBA from Bonn, use example fiel to implement
     //    
}
// --------------------------------------------------------------------------------------------------------------
bool LGADBase::AddTrckBranches(TTree* c_tree, std::string trckpkg)
{
    std::cout << "    Creating output ROOT file: " << m_ofile->GetName() << std::endl;

    if (trckpkg == "Cory")
       {
        // Find out the number of unique elements for the references, DUTs and planes vectors
        m_ntp_planes.clear();
        m_ntp_Ref.clear();
        m_ntp_DUT.clear();
        bool exist = false;
        for (unsigned int y = 0; y < 3; y++)
            {
             std::vector<unsigned int> inputvec;
             std::vector<unsigned int> outputvec;
             inputvec.clear();
             outputvec.clear();
             if (y == 0) inputvec = m_planes; 
             else if (y == 1) inputvec = m_DUTs;
             else inputvec = m_Refs;
             for (unsigned int k = 0; k < inputvec.size(); k++)
                 {
                  exist = false;
                  for (unsigned int a = 0; a < outputvec.size(); a++)
                      {
                       if (inputvec.at(k) == outputvec.at(a))
                          {
                           exist = true;
                           break;
                          }
                      }
                  if (!exist) outputvec.push_back(inputvec.at(k));
                 }
             if (y == 0) m_ntp_planes = outputvec;
             else if (y == 1) m_ntp_DUT = outputvec;
             else m_ntp_Ref = outputvec;
            }
        std::sort(m_ntp_planes.begin(), m_ntp_planes.end());
        std::sort(m_ntp_DUT.begin(), m_ntp_DUT.end());
        std::sort(m_ntp_Ref.begin(), m_ntp_Ref.end());

        // Initialise enougph pointers to contain the maxium number of unique elements
        LGADBase::InitCoryPointers(m_ntp_Ref.size(), m_ntp_DUT.size(), m_ntp_planes.size());

        // Create the ROI branches at the output Ntuple 
        c_tree->Branch("ClusterSizeX_ROI", &m_ClusterSizeX_ROI);
        c_tree->Branch("ClusterSizeY_ROI", &m_ClusterSizeY_ROI);
        c_tree->Branch("PixelX_ROI", &m_PixelX_ROI);
        c_tree->Branch("PixelY_ROI", &m_PixelY_ROI);
        c_tree->Branch("ClusterNumPixels_ROI", &m_ClusterNumPixels_ROI);
        c_tree->Branch("InterceptX_ROI", &m_InterceptX_ROI);
        c_tree->Branch("InterceptY_ROI", &m_InterceptY_ROI);
        c_tree->Branch("InterceptResX_ROI", &m_InterceptResX_ROI);
        c_tree->Branch("InterceptResY_ROI", &m_InterceptResY_ROI);
        c_tree->Branch("HasAsso", &m_HasAsso);
        // Create the Reference branches at the output Ntuple      
        for (unsigned int ich = 0; ich < m_ntp_Ref.size(); ich++)
            {
             c_tree->Branch(Form("InterceptX_Ref%u", m_ntp_Ref.at(ich)), &m_InterceptX_Ref.at(ich));
             c_tree->Branch(Form("InterceptY_Ref%u", m_ntp_Ref.at(ich)), &m_InterceptY_Ref.at(ich));
             c_tree->Branch(Form("InterceptResX_Ref%u", m_ntp_Ref.at(ich)), &m_InterceptResX_Ref.at(ich));
             c_tree->Branch(Form("InterceptResY_Ref%u", m_ntp_Ref.at(ich)), &m_InterceptResY_Ref.at(ich));
            }
        // Create the DUT branches at the output Ntuple  
        for (unsigned int ich = 0; ich < m_ntp_DUT.size(); ich++)
            {
             c_tree->Branch(Form("InterceptX_DUT%u", m_ntp_DUT.at(ich)), &m_InterceptX_DUT.at(ich));
             c_tree->Branch(Form("InterceptY_DUT%u", m_ntp_DUT.at(ich)), &m_InterceptY_DUT.at(ich));
             c_tree->Branch(Form("InterceptResX_DUT%u", m_ntp_DUT.at(ich)), &m_InterceptResX_DUT.at(ich));
             c_tree->Branch(Form("InterceptResY_DUT%u", m_ntp_DUT.at(ich)), &m_InterceptResY_DUT.at(ich));
            }
        // Create the planes branches at the output Ntuple 
        for (unsigned int ich = 0; ich < m_ntp_planes.size(); ich++)
            {
             c_tree->Branch(Form("InterceptX_TelPlane%u", m_ntp_planes.at(ich)), &m_InterceptX_TelPlane.at(ich));
             c_tree->Branch(Form("InterceptY_TelPlane%u", m_ntp_planes.at(ich)), &m_InterceptY_TelPlane.at(ich));
             c_tree->Branch(Form("InterceptResX_TelPlane%u", m_ntp_planes.at(ich)), &m_InterceptResX_TelPlane.at(ich));
             c_tree->Branch(Form("InterceptResY_TelPlane%u", m_ntp_planes.at(ich)), &m_InterceptResY_TelPlane.at(ich));
             c_tree->Branch(Form("ClusterSizeX_TelPlane%u", m_ntp_planes.at(ich)), &m_ClusterSizeX_TelPlane.at(ich));
             c_tree->Branch(Form("ClusterSizeY_TelPlane%u", m_ntp_planes.at(ich)), &m_ClusterSizeY_TelPlane.at(ich));
             c_tree->Branch(Form("ClusterNumPixels_TelPlane%u", m_ntp_planes.at(ich)), &m_ClusterNumPixels_TelPlane.at(ich));
            }
        // Printout the created branches
        PrintBranches(c_tree, "Waveform Cory combined");
        PrintObjects(m_ofile, "wfm");
        return true;
       }
    else if (trckpkg == "EUTelescope")
            {
             LGADBase::InitEUTelPointers();
             c_tree->Branch("xTr", &m_xTr);
             c_tree->Branch("yTr", &m_yTr);
             c_tree->Branch("omegaTr", &m_omegaTr);
             c_tree->Branch("phiTr", &m_phiTr);
             c_tree->Branch("kinkxTr", &m_kinkxTr);
             c_tree->Branch("kinkyTr", &m_kinkyTr);
             c_tree->Branch("chi2Tr", &m_chi2Tr);
             c_tree->Branch("plane_id", &m_plane_id);
             c_tree->Branch("track_id", &m_track_id);
             c_tree->Branch("ndofTr", &m_ndofTr);
             c_tree->SetBranchAddress("xHit", &m_xHit);
             c_tree->SetBranchAddress("yHit", &m_yHit);
             c_tree->SetBranchAddress("zHit", &m_zHit);
             c_tree->SetBranchAddress("ID", &m_ID);
             c_tree->Branch("trigger_id", &m_trigger_id, "m_trigger_id/i");
             c_tree->Branch("timestampTr", &m_timestampTr, "timestampTr/i");

             // Printout the created branches
             PrintBranches(c_tree, "Waveform EUTelescpe combined");
             PrintObjects(m_ofile, "wfm");
             return true;
            }
    else if (trckpkg == "TBA")
            {
             LGADBase::InitTBAPointers();
                //
                // Missingn code for TBA from Bonn, use example fiel to implement
                //               
             // Printout the created branches
             PrintBranches(c_tree, "Waveform BTA combined");
             PrintObjects(m_ofile, "wfm");
             return true;
            }
    else {
          std::cout << __FUNCTION__ << " WARNING: Unrecognised tracking framework: " << trckpkg << std::endl;
          return false;
         }
}
// --------------------------------------------------------------------------------------------------------------
bool LGADBase::CombineTrackCory(int EvntNmbr, std::vector<unsigned int> &evntlast, unsigned int file)
{
    if (EvntNmbr < 0)
       {
        std::cout << __FUNCTION__ << " EROOR: Negative event number passed to the tracking function, abording.... " << std::endl;
        return false;
       }
    // Vector for controlling that all registered branches are getting populated
    std::vector<unsigned int> control;

    // Vector pointer cleanup from previous event, once cannot destroy and reacreate the pointers as the InitCory pointers function does
    // because the assignement t branches will be lost. Vectors of ointers need to be maintained ans individual allements set to NULL
    m_ClusterSizeX_ROI = NULL;
    m_ClusterSizeY_ROI = NULL;
    m_PixelX_ROI = NULL;
    m_PixelY_ROI = NULL;
    m_InterceptX_ROI = NULL;
    m_InterceptY_ROI = NULL;
    m_InterceptResX_ROI = NULL;
    m_InterceptResY_ROI = NULL;
    m_HasAsso = NULL;
    m_ClusterNumPixels_ROI = NULL;
    for (unsigned int k = 0; k < m_ntp_planes.size(); k++)
        {
         (m_InterceptX_TelPlane.at(k)) = NULL;
         (m_InterceptY_TelPlane.at(k)) = NULL;
         (m_InterceptResX_TelPlane.at(k)) = NULL;
         (m_InterceptResY_TelPlane.at(k)) = NULL;
         (m_ClusterSizeX_TelPlane.at(k)) = NULL;
         (m_ClusterSizeY_TelPlane.at(k)) = NULL;
         (m_ClusterNumPixels_TelPlane.at(k)) = NULL;
        }
    for (unsigned int k = 0; k < m_ntp_Ref.size(); k++)
        {
         (m_InterceptX_Ref.at(k)) = NULL;
         (m_InterceptY_Ref.at(k)) = NULL;
         (m_InterceptResX_Ref.at(k)) = NULL;
         (m_InterceptResY_Ref.at(k)) = NULL;
        }
    for (unsigned int k = 0; k < m_ntp_DUT.size(); k++)
        {
         (m_InterceptX_DUT.at(k)) = NULL;
         (m_InterceptY_DUT.at(k)) = NULL;
         (m_InterceptResX_DUT.at(k)) = NULL;
         (m_InterceptResY_DUT.at(k)) = NULL;
        }

    // Beaceause of the flattened vectors for trees and branches, calculate the correct indexes to recover data
    unsigned int tmpindx = 0; // Temporary index for several later calculations
    unsigned int tmptreeindx = 0;  // Temporary index for trees
    unsigned int plnprfx = 0;
    unsigned int dutprfx = 0;
    unsigned int refprfx = 0;
    for (unsigned int i = 0; i < file; i++)
        {
         plnprfx += m_nplanes.at(i);
         dutprfx += m_nDUT.at(i);
         refprfx += m_nRef.at(i);
        }

    bool found = false;
    int multindx = -1;
  
    // ============================ Process the ROI tree ============================
    tmptreeindx = plnprfx + dutprfx + refprfx + file;
    for (unsigned int k = evntlast.at(0); k < m_trackTrees.at(tmptreeindx)->GetEntries(); k++)
        {
         tmpindx = 11*file;
         // Get the EventId
         b_m_ROI.at(tmpindx+10)->GetEntry(k);
         if ((m_ROI_Int.at(3*file+2))->at(0) == EvntNmbr)
            { 
             found = true;
             for (unsigned int m = 0; m < 10; m++) b_m_ROI.at(tmpindx + m)->GetEntry(k);
             tmpindx = 8*file;
             m_InterceptX_ROI = m_ROI_Dbl.at(tmpindx);
             m_InterceptY_ROI = m_ROI_Dbl.at(tmpindx+1);
             m_InterceptResX_ROI = m_ROI_Dbl.at(tmpindx+2);
             m_InterceptResY_ROI = m_ROI_Dbl.at(tmpindx+3);
             m_ClusterSizeX_ROI = m_ROI_Dbl.at(tmpindx+4);
             m_ClusterSizeY_ROI = m_ROI_Dbl.at(tmpindx+5);
             m_PixelX_ROI = m_ROI_Dbl.at(tmpindx+6);
             m_PixelY_ROI = m_ROI_Dbl.at(tmpindx+7);
             tmpindx = 3*file;
             m_ClusterNumPixels_ROI = m_ROI_Int.at(tmpindx);
             m_HasAsso = m_ROI_Int.at(tmpindx+1);
             evntlast.at(0)++;
            }
         else break;
        }
    if (!found) // If there is no coresponding event in the cory file
       { 
        m_InterceptX_ROI->push_back(-99);
        m_InterceptY_ROI->push_back(-99);
        m_InterceptResX_ROI->push_back(-99);
        m_InterceptResY_ROI->push_back(-99);
        m_ClusterSizeX_ROI->push_back(-99);
        m_ClusterSizeY_ROI->push_back(-99);
        m_PixelX_ROI->push_back(-99);
        m_PixelY_ROI->push_back(-99);
        m_ClusterNumPixels_ROI->push_back(-99);
        m_HasAsso->push_back(-99);
       }

    // ============================ Process telescope plane trees ============================
    control = m_ntp_planes;
    for (unsigned int y = 0; y < m_nplanes.at(file); y++)
        {
         // Step 1: Subtract the current plane from the vector of planes to make sure that all planes were present at the end
         for (unsigned int g = 0; g < control.size(); g++)  
             {
              if (control.at(g) == m_planes.at(plnprfx+y))
                 {
                  control.erase(control.begin()+g);
                  break;
                 }
              }
         // Step 2: Find out which index of the output Ntuple vectors this plane's quntities need to bbe assigned to
         multindx = -1;
         for (unsigned int a = 0; a < m_ntp_planes.size(); a++)
             {
              if (m_ntp_planes.at(a) == m_planes.at(plnprfx + y))
                 {
                  multindx = a;
                  break;
                 }
             }
         if (multindx == -1)
            {
             std::cout << __FUNCTION__ << " ERROR: Cannot determine which Ntuple branch to store current plane: " << m_planes.at(plnprfx+y) << std::endl;
             return false;
            }
         // Step 3: Find out if the event we are looking for has an entry and recover branches
         found = false;
         tmptreeindx = plnprfx+dutprfx+refprfx+file+1+y;
         for (unsigned int k = evntlast.at(y+1); k < m_trackTrees.at(tmptreeindx)->GetEntries(); k++)
             {
              tmpindx = 8*(plnprfx+y);
              b_m_TelPlane.at(tmpindx+7)->GetEntry(k);
              if (m_TelPlane_Int.at(2*(plnprfx+y)+1)->at(0) == EvntNmbr)
                 {
                  found = true;
                  for (unsigned int m = 0; m < 7; m++) b_m_TelPlane.at(tmpindx + m)->GetEntry(k);
                  tmpindx = 6*(plnprfx+y);
                  m_InterceptX_TelPlane.at(multindx) = m_TelPlane_Dbl.at(tmpindx);
                  m_InterceptY_TelPlane.at(multindx) = m_TelPlane_Dbl.at(tmpindx+1);
                  m_InterceptResX_TelPlane.at(multindx) = m_TelPlane_Dbl.at(tmpindx+2);
                  m_InterceptResY_TelPlane.at(multindx) = m_TelPlane_Dbl.at(tmpindx+3);
                  m_ClusterSizeX_TelPlane.at(multindx) = m_TelPlane_Dbl.at(tmpindx+4);
                  m_ClusterSizeY_TelPlane.at(multindx) = m_TelPlane_Dbl.at(tmpindx+5);
                  m_ClusterNumPixels_TelPlane.at(multindx) = m_TelPlane_Int.at(2*(plnprfx+y));
                  evntlast.at(y+1)++;
                 }
              else break;
             }
         // Step 4: If the event was not found in the previous step, then populate the vectors with -99
         if (!found) // If there is no coresponding event in the cory file
            {
             (m_InterceptX_TelPlane.at(multindx))->push_back(-99);
             (m_InterceptY_TelPlane.at(multindx))->push_back(-99);
             (m_InterceptResX_TelPlane.at(multindx))->push_back(-99);
             (m_InterceptResY_TelPlane.at(multindx))->push_back(-99);
             (m_ClusterSizeX_TelPlane.at(multindx))->push_back(-99);
             (m_ClusterSizeY_TelPlane.at(multindx))->push_back(-99);
             (m_ClusterNumPixels_TelPlane.at(multindx))->push_back(-99);
            }
         if (LGADBase::GetVerbosity() >= 1) std::cout << __FUNCTION__ << " INFO: Finished planes tree " << tmptreeindx << std::endl;
        }
    // Step 5: If the plane was not present in this file, populate the vectors only with one entry of -100
    for (unsigned int h = 0; h < control.size(); h++)
        {
         (m_InterceptX_TelPlane.at(control.at(h)))->push_back(-100);
         (m_InterceptY_TelPlane.at(control.at(h)))->push_back(-100);
         (m_InterceptResX_TelPlane.at(control.at(h)))->push_back(-100);
         (m_InterceptResY_TelPlane.at(control.at(h)))->push_back(-100);
         (m_ClusterSizeX_TelPlane.at(control.at(h)))->push_back(-100);
         (m_ClusterSizeY_TelPlane.at(control.at(h)))->push_back(-100);
         (m_ClusterNumPixels_TelPlane.at(control.at(h)))->push_back(-100);
        }
    if (LGADBase::GetVerbosity() >= 1) std::cout << __FUNCTION__ << " INFO: We finished looking at all telescope plane trees" << std::endl;

    // ============================ Process References trees ============================
    control = m_ntp_Ref;
    for (unsigned int y = 0; y < m_nRef.at(file); y++)
        {
         // Step 1: Subtract the current reference from the vector of references to make sure that all references were present at the end
         for (unsigned int g = 0; g < control.size(); g++)  
             {
              if (control.at(g) == m_Refs.at(refprfx+y))
                 {
                  control.erase(control.begin()+g);
                  break;
                 }
              } 
         // Step 2: Find out which index of the output Ntuple vectors this reference's quntities need to bbe assigned to
         multindx = -1;
         for (unsigned int a = 0; a < m_ntp_Ref.size(); a++)
             {
              if (m_ntp_Ref.at(a) == m_Refs.at(refprfx + y))
                 {
                  multindx = a;
                  break;
                 }
             }
         if (multindx == -1)
            {
             std::cout << __FUNCTION__ << " ERROR: Cannot determine which Ntuple branch to store current reference: " << m_Refs.at(refprfx+y) << std::endl;
             return false;
            }
         // Step 3: Find out if the event we are looking for has an entry and recover branches
         found = false;
         tmptreeindx = plnprfx+dutprfx+refprfx+file+1+m_nplanes.at(file)+y;
         for (unsigned int k = evntlast.at(m_nplanes.at(file)+1+y); k < m_trackTrees.at(tmptreeindx)->GetEntries(); k++)
             {
              tmpindx = 5*(refprfx+y);
              b_m_Ref.at(tmpindx+4)->GetEntry(k);
              if ((m_Ref_Int.at(refprfx+y))->at(0) == EvntNmbr)
                 {
                  found = true;
                  for (unsigned int m = 0; m < 4; m++) b_m_Ref.at(tmpindx + m)->GetEntry(k);
                  tmpindx = 4*(refprfx+y);
                  m_InterceptX_Ref.at(multindx) = m_Ref_Dbl.at(tmpindx);
                  m_InterceptY_Ref.at(multindx) = m_Ref_Dbl.at(tmpindx+1);
                  m_InterceptResX_Ref.at(multindx) = m_Ref_Dbl.at(tmpindx+2);
                  m_InterceptResY_Ref.at(multindx) = m_Ref_Dbl.at(tmpindx+3);
                  evntlast.at(m_nplanes.at(file)+1+y)++;
                 }
              else break;
             }
         // Step 4: If the event was not found in the previous step, then populate the vectors with -99
         if (!found) // If there is no coresponding event in the cory file
            {
             (m_InterceptX_Ref.at(multindx))->push_back(-99);
             (m_InterceptY_Ref.at(multindx))->push_back(-99);
             (m_InterceptResX_Ref.at(multindx))->push_back(-99);
             (m_InterceptResY_Ref.at(multindx))->push_back(-99);
            }
         if (LGADBase::GetVerbosity() >= 1) std::cout << __FUNCTION__ << " INFO: Finished reference tree " << tmptreeindx << std::endl;
        }
    // Step 5: If the reference was not present in this file, populate the vectors only with one entry of -100
    for (unsigned int h = 0; h < control.size(); h++)
        {
         (m_InterceptX_Ref.at(control.at(h)))->push_back(-100);
         (m_InterceptY_Ref.at(control.at(h)))->push_back(-100);
         (m_InterceptResX_Ref.at(control.at(h)))->push_back(-100);
         (m_InterceptResY_Ref.at(control.at(h)))->push_back(-100);
        }
    if (LGADBase::GetVerbosity() >= 1) std::cout << __FUNCTION__ << " INFO: We finished looking at all reference plane trees" << std::endl;

    // ============================ Process DUT trees ============================
    control = m_ntp_DUT;
    for (unsigned int y = 0; y < m_nDUT.at(file); y++)
        {
        // Step 1: Subtract the current DUT from the vector of DUTs to make sure that all devices were present at the end
         for (unsigned int g = 0; g < control.size(); g++)  
             {
              if (control.at(g) == m_DUTs.at(dutprfx+y))
                 {
                  control.erase(control.begin()+g);
                  break;
                 }
              } 
         // Step 2: Find out which index of the output Ntuple vectors this DUT's quntities need to bbe assigned to
         multindx = -1;
         for (unsigned int a = 0; a < m_ntp_DUT.size(); a++)
             {
              if (m_ntp_DUT.at(a) == m_DUTs.at(dutprfx + y))
                 {
                  multindx = a;
                  break;
                 }
             }
         if (multindx == -1)
            {
             std::cout << __FUNCTION__ << " ERROR: Cannot determine which Ntuple branch to store current DUT: " << m_DUTs.at(dutprfx+y) << std::endl;
             return false;
            }
         // Step 3: Find out if the event we are looking for has an entry and recover branches
         found = false;
         tmptreeindx = plnprfx+dutprfx+refprfx+file+1+m_nplanes.at(file)+m_nRef.at(file)+y;
         for (unsigned int k = evntlast.at(m_nplanes.at(file)+1+m_nRef.at(file)+y); k < m_trackTrees.at(tmptreeindx)->GetEntries(); k++)
             {
              tmpindx = 5*(dutprfx+y);
              b_m_DUT.at(tmpindx+4)->GetEntry(k);
              if ((m_DUT_Int.at(dutprfx+y))->at(0) == EvntNmbr)
                 {
                  found = true;
                  for (unsigned int m = 0; m < 4; m++) b_m_DUT.at(tmpindx + m)->GetEntry(k);
                  tmpindx = 4*(dutprfx+y);
                  m_InterceptX_DUT.at(multindx) = m_DUT_Dbl.at(tmpindx);
                  m_InterceptY_DUT.at(multindx) = m_DUT_Dbl.at(tmpindx+1);
                  m_InterceptResX_DUT.at(multindx) = m_DUT_Dbl.at(tmpindx+2);
                  m_InterceptResY_DUT.at(multindx) = m_DUT_Dbl.at(tmpindx+3);
                  evntlast.at(m_nplanes.at(file)+1+m_nRef.at(file)+y)++;
                 }
              else break;
             }
         // Step 4: If the event was not found in the previous step, then populate the vectors with -99
         if (!found) // If there is no coresponding event in the cory file
            {
             (m_InterceptX_DUT.at(multindx))->push_back(-99);
             (m_InterceptY_DUT.at(multindx))->push_back(-99);
             (m_InterceptResX_DUT.at(multindx))->push_back(-99);
             (m_InterceptResY_DUT.at(multindx))->push_back(-99);
            }
         if (LGADBase::GetVerbosity() >= 1) std::cout << __FUNCTION__ << " INFO: Finished DUT tree " << tmptreeindx << std::endl;
        }
    // Step 5: If the DUT was not present in this file, populate the vectors only with one entry of -100
    for (unsigned int h = 0; h < control.size(); h++)
        {
         (m_InterceptX_DUT.at(control.at(h)))->push_back(-100);
         (m_InterceptY_DUT.at(control.at(h)))->push_back(-100);
         (m_InterceptResX_DUT.at(control.at(h)))->push_back(-100);
         (m_InterceptResY_DUT.at(control.at(h)))->push_back(-100);
        }
    if (LGADBase::GetVerbosity() >= 1) std::cout << __FUNCTION__ << " INFO: We finished looking at all DUT trees" << std::endl;

   return true;
}
// --------------------------------------------------------------------------------------------------------------
bool LGADBase::CombineTrackEUTel(int EvntNmbr, std::vector<unsigned int>& evntlast, unsigned int file)
{
    // Clear pointers from previous event
    LGADBase::InitEUTelPointers();
    m_trigger_id = 0;
    m_timestampTr = 0;

    //
    // Missing code for EUTelescope, use example fiel to implement
    //    
    return true;
}
// --------------------------------------------------------------------------------------------------------------
bool LGADBase::CombineTrackTBA(int EvntNmbr, std::vector<unsigned int>& evntlast, unsigned int file)
{
    //
    // Missing code for TBA from Bonn, use example fiel to implement
    //   
    return true;
}
// --------------------------------------------------------------------------------------------------------------
void LGADBase::CloseTrackFiles()
{
    for (unsigned int k = 0; k < m_trackfiles.size(); k++) m_trackfiles.at(k)->Close();
    SetTrackComb(false);
} 
// --------------------------------------------------------------------------------------------------------------