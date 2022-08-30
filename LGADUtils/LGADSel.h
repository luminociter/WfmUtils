/*
* LGADSel.h
*
*
*      Author: Gkougkousis Evangelos - Leonidas
*              egkougko@cern.ch
*               IFAE-BARCELONA
*
* File derived from a staandard root selector with additions on 
* automatic instrument asignement and branch matching
*/

#ifndef __LGADSel__
#define __LGADSel__

#include <TChain.h>
#include <TH2.h>
#include <TStyle.h>

#include "WaveForm.h"

// Fixed size dimensions of array or collections stored in the TTree if any.

class LGADSel : public LGADBase {

public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   std::vector<std::vector<double>*> t;
   std::vector<std::vector<double>*> w;

   // Tracking leaf types 
   std::vector<double> *xTr;
   std::vector<double> *yTr;
   std::vector<double> *omegaTr;
   std::vector<double> *phiTr;
   std::vector<double> *kinkxTr;
   std::vector<double> *kinkyTr;
   std::vector<double> *chi2Tr;
   std::vector<int> *plane_id;
   std::vector<int> *track_id;
   std::vector<int> *ndofTr;

   std::vector<int> *dutID;
   std::vector<double> *xHit;
   std::vector<double> *yHit;
   std::vector<double> *zHit;

   UInt_t          EvnNo;
   Int_t       trigger_id;
   Int_t       timestampTr;

   // List of branches
   std::vector<TBranch*> b_t;   //!
   std::vector<TBranch*> b_w;   //!
   std::vector<TBranch*> b_vScale;   //!
   std::vector<TBranch*> b_nPoints;   //!
   std::vector<TBranch*> b_SnRate;   //!
   std::vector<TBranch*> b_ordrt;   //!   
   std::vector<TBranch*> b_physt;   //!
   std::vector<TBranch*> b_triggTime;   //!
   std::vector<TBranch*> b_Scope;   //!

   TBranch        *b_trigtime;   //!
   TBranch        *b_EvnNo;   //!

   // Tracking Branches Pointers
   TBranch        *b_plane_id;   //!
   TBranch        *b_track_id;   //!
   TBranch        *b_trigger_id;   //!
   TBranch        *b_timestampTr;   //!
   TBranch        *b_xTr;   //!
   TBranch        *b_yTr;   //!
   TBranch        *b_omegaTr;   //!
   TBranch        *b_phiTr;   //!
   TBranch        *b_kinkxTr;   //!
   TBranch        *b_kinkyTr;   //!
   TBranch        *b_chi2Tr;   //!
   TBranch        *b_ndofTr;   //!
   
   TBranch        *b_dutID;   //!
   TBranch        *b_xHit;   //!
   TBranch        *b_yHit;   //!
   TBranch        *b_zHit;   //!

   LGADSel(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~LGADSel() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList   *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   // ClassDef(LGADSel, 1);
};

#endif
