//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Apr  1 21:20:30 2019 by ROOT version 6.12/07
// from TTree ntp1/ntp1
// found on file: /eos/cms/store/group/phys_pps/ProtonRecoValidation/Dileptons/2018/ProtonReconstruction_DoubleMu2018B_merge.root
//////////////////////////////////////////////////////////

#ifndef ntp1_h
#define ntp1_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class ntp1 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          Run;
   UInt_t          LumiSection;
   UInt_t          BX;
   UInt_t          EventNum;
   Float_t         CrossingAngle;
   UInt_t          nHLT;
   Int_t           HLT_Accept[3];   //[nHLT]
   Int_t           HLT_Prescl[3];   //[nHLT]
   vector<string>  *HLT_Name;
   UInt_t          nMuonCand;
   Double_t        MuonCand_pt[5];   //[nMuonCand]
   Double_t        MuonCand_eta[5];   //[nMuonCand]
   Double_t        MuonCand_phi[5];   //[nMuonCand]
   Double_t        MuonCand_e[5];   //[nMuonCand]
   Int_t           MuonCand_charge[5];   //[nMuonCand]
   Double_t        MuonCand_vtxx[5];   //[nMuonCand]
   Double_t        MuonCand_vtxy[5];   //[nMuonCand]
   Double_t        MuonCand_vtxz[5];   //[nMuonCand]
   Double_t        MuonCand_dxy[5];   //[nMuonCand]
   Int_t           MuonCand_nstatseg[5];   //[nMuonCand]
   Int_t           MuonCand_ntrklayers[5];   //[nMuonCand]
   Int_t           MuonCand_npxlhits[5];   //[nMuonCand]
   Int_t           MuonCand_isglobal[5];   //[nMuonCand]
   Int_t           MuonCand_istracker[5];   //[nMuonCand]
   Int_t           MuonCand_isstandalone[5];   //[nMuonCand]
   Int_t           MuonCand_ispfmuon[5];   //[nMuonCand]
   Int_t           MuonCand_istight[5];   //[nMuonCand]
   Int_t           MuonCandTrack_nmuchits[5];   //[nMuonCand]
   Double_t        MuonCandTrack_chisq[5];   //[nMuonCand]
   Double_t        MuonCand_innerTrackPt[5];   //[nMuonCand]
   Double_t        MuonCand_innerTrackEta[5];   //[nMuonCand]
   Double_t        MuonCand_innerTrackPhi[5];   //[nMuonCand]
   Double_t        MuonCand_innerTrackVtxz[5];   //[nMuonCand]
   UInt_t          nPrimVertexCand;
   Double_t        PrimVertexCand_x[133];   //[nPrimVertexCand]
   Double_t        PrimVertexCand_y[133];   //[nPrimVertexCand]
   Double_t        PrimVertexCand_z[133];   //[nPrimVertexCand]
   Double_t        PrimVertexCand_chi2[133];   //[nPrimVertexCand]
   UInt_t          PrimVertexCand_ndof[133];   //[nPrimVertexCand]
   UInt_t          PrimVertexCand_tracks[133];   //[nPrimVertexCand]
   UInt_t          nPair;
   Int_t           Pair_lepton1[6];   //[nPair]
   Int_t           Pair_lepton2[6];   //[nPair]
   Double_t        Pair_mass[6];   //[nPair]
   Double_t        Pair_pt[6];   //[nPair]
   Double_t        Pair_eta[6];   //[nPair]
   Double_t        Pair_phi[6];   //[nPair]
   Double_t        Pair_dpt[6];   //[nPair]
   Double_t        Pair_dphi[6];   //[nPair]
   Double_t        Pair_3Dangle[6];   //[nPair]
   UInt_t          Pair_extratracks0p5mm[6];   //[nPair]
   UInt_t          Pair_extratracks1mm[6];   //[nPair]
   UInt_t          Pair_extratracks2mm[6];   //[nPair]
   UInt_t          Pair_extratracks3mm[6];   //[nPair]
   UInt_t          Pair_extratracks4mm[6];   //[nPair]
   UInt_t          Pair_extratracks5mm[6];   //[nPair]
   UInt_t          Pair_extratracks1cm[6];   //[nPair]
   UInt_t          Pair_extratracks2cm[6];   //[nPair]
   UInt_t          Pair_extratracks3cm[6];   //[nPair]
   UInt_t          Pair_extratracks4cm[6];   //[nPair]
   UInt_t          Pair_extratracks5cm[6];   //[nPair]
   UInt_t          Pair_extratracks10cm[6];   //[nPair]
   Double_t        KalmanVertexCand_x[6];   //[nPair]
   Double_t        KalmanVertexCand_y[6];   //[nPair]
   Double_t        KalmanVertexCand_z[6];   //[nPair]
   UInt_t          nLocalProtCand;
   Double_t        LocalProtCand_x[34];   //[nLocalProtCand]
   Double_t        LocalProtCand_y[34];   //[nLocalProtCand]
   Double_t        LocalProtCand_t[34];   //[nLocalProtCand]
   Double_t        LocalProtCand_xSigma[34];   //[nLocalProtCand]
   Double_t        LocalProtCand_ySigma[34];   //[nLocalProtCand]
   Double_t        LocalProtCand_tSigma[34];   //[nLocalProtCand]
   Int_t           LocalProtCand_arm[34];   //[nLocalProtCand]
   Int_t           LocalProtCand_station[34];   //[nLocalProtCand]
   Int_t           LocalProtCand_pot[34];   //[nLocalProtCand]
   Int_t           LocalProtCand_rpid[34];   //[nLocalProtCand]
   UInt_t          nRecoProtCand;
   Double_t        ProtCand_xi[28];   //[nRecoProtCand]
   Double_t        ProtCand_t[28];   //[nRecoProtCand]
   Double_t        ProtCand_ThX[28];   //[nRecoProtCand]
   Double_t        ProtCand_ThY[28];   //[nRecoProtCand]
   Int_t           ProtCand_rpid[28];   //[nRecoProtCand]
   Int_t           ProtCand_arm[28];   //[nRecoProtCand]
   Int_t           ProtCand_ismultirp[28];   //[nRecoProtCand]
   UInt_t          nQualityExtraTrack;
   Double_t        ClosestExtraTrack_vtxdxyz[6];   //[nPair]
   Int_t           ClosestExtraTrack_id[6];   //[nPair]
   Double_t        ClosestHighPurityExtraTrack_vtxdxyz[6];   //[nPair]
   Int_t           ClosestHighPurityExtraTrack_id[6];   //[nPair]
   Double_t        Weight;
   Double_t        PUWeightTrue;

   // List of branches
   TBranch        *b_Run;   //!
   TBranch        *b_LumiSection;   //!
   TBranch        *b_BX;   //!
   TBranch        *b_EventNum;   //!
   TBranch        *b_CrossingAngle;   //!
   TBranch        *b_nHLT;   //!
   TBranch        *b_HLT_Accept;   //!
   TBranch        *b_HLT_Prescl;   //!
   TBranch        *b_HLT_Name;   //!
   TBranch        *b_nMuonCand;   //!
   TBranch        *b_MuonCand_pt;   //!
   TBranch        *b_MuonCand_eta;   //!
   TBranch        *b_MuonCand_phi;   //!
   TBranch        *b_MuonCand_e;   //!
   TBranch        *b_MuonCand_charge;   //!
   TBranch        *b_MuonCand_vtxx;   //!
   TBranch        *b_MuonCand_vtxy;   //!
   TBranch        *b_MuonCand_vtxz;   //!
   TBranch        *b_MuonCand_dxy;   //!
   TBranch        *b_MuonCand_nstatseg;   //!
   TBranch        *b_MuonCand_ntrklayers;   //!
   TBranch        *b_MuonCand_npxlhits;   //!
   TBranch        *b_MuonCand_isglobal;   //!
   TBranch        *b_MuonCand_istracker;   //!
   TBranch        *b_MuonCand_isstandalone;   //!
   TBranch        *b_MuonCand_ispfmuon;   //!
   TBranch        *b_MuonCand_istight;   //!
   TBranch        *b_MuonCandTrack_nmuchits;   //!
   TBranch        *b_MuonCandTrack_chisq;   //!
   TBranch        *b_MuonCand_innerTrackPt;   //!
   TBranch        *b_MuonCand_innerTrackEta;   //!
   TBranch        *b_MuonCand_innerTrackPhi;   //!
   TBranch        *b_MuonCand_innerTrackVtxz;   //!
   TBranch        *b_nPrimVertexCand;   //!
   TBranch        *b_PrimVertexCand_x;   //!
   TBranch        *b_PrimVertexCand_y;   //!
   TBranch        *b_PrimVertexCand_z;   //!
   TBranch        *b_PrimVertexCand_chi2;   //!
   TBranch        *b_PrimVertexCand_ndof;   //!
   TBranch        *b_PrimVertexCand_tracks;   //!
   TBranch        *b_nPair;   //!
   TBranch        *b_Pair_lepton1;   //!
   TBranch        *b_Pair_lepton2;   //!
   TBranch        *b_Pair_mass;   //!
   TBranch        *b_Pair_pt;   //!
   TBranch        *b_Pair_eta;   //!
   TBranch        *b_Pair_phi;   //!
   TBranch        *b_Pair_dpt;   //!
   TBranch        *b_Pair_dphi;   //!
   TBranch        *b_Pair_3Dangle;   //!
   TBranch        *b_Pair_extratracks0p5mm;   //!
   TBranch        *b_Pair_extratracks1mm;   //!
   TBranch        *b_Pair_extratracks2mm;   //!
   TBranch        *b_Pair_extratracks3mm;   //!
   TBranch        *b_Pair_extratracks4mm;   //!
   TBranch        *b_Pair_extratracks5mm;   //!
   TBranch        *b_Pair_extratracks1cm;   //!
   TBranch        *b_Pair_extratracks2cm;   //!
   TBranch        *b_Pair_extratracks3cm;   //!
   TBranch        *b_Pair_extratracks4cm;   //!
   TBranch        *b_Pair_extratracks5cm;   //!
   TBranch        *b_Pair_extratracks10cm;   //!
   TBranch        *b_KalmanVertexCand_x;   //!
   TBranch        *b_KalmanVertexCand_y;   //!
   TBranch        *b_KalmanVertexCand_z;   //!
   TBranch        *b_nLocalProtCand;   //!
   TBranch        *b_LocalProtCand_x;   //!
   TBranch        *b_LocalProtCand_y;   //!
   TBranch        *b_LocalProtCand_t;   //!
   TBranch        *b_LocalProtCand_xSigma;   //!
   TBranch        *b_LocalProtCand_ySigma;   //!
   TBranch        *b_LocalProtCand_tSigma;   //!
   TBranch        *b_LocalProtCand_arm;   //!
   TBranch        *b_LocalProtCand_station;   //!
   TBranch        *b_LocalProtCand_pot;   //!
   TBranch        *b_LocalProtCand_rpid;   //!
   TBranch        *b_nRecoProtCand;   //!
   TBranch        *b_ProtCand_xi;   //!
   TBranch        *b_ProtCand_t;   //!
   TBranch        *b_ProtCand_ThX;   //!
   TBranch        *b_ProtCand_ThY;   //!
   TBranch        *b_ProtCand_rpid;   //!
   TBranch        *b_ProtCand_arm;   //!
   TBranch        *b_ProtCand_ismultirp;   //!
   TBranch        *b_nQualityExtraTrack;   //!
   TBranch        *b_ClosestExtraTrack_vtxdxyz;   //!
   TBranch        *b_ClosestExtraTrack_id;   //!
   TBranch        *b_ClosestHighPurityExtraTrack_vtxdxyz;   //!
   TBranch        *b_ClosestHighPurityExtraTrack_id;   //!
   TBranch        *b_Weight;   //!
   TBranch        *b_PUWeightTrue;   //!

   ntp1(TTree *tree=0);
   virtual ~ntp1();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ntp1_cxx
ntp1::ntp1(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/eos/cms/store/group/phys_pps/ProtonRecoValidation/Dileptons/2017/ProtonReconstruction_DoubleMu2017F_merge.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/eos/cms/store/group/phys_pps/ProtonRecoValidation/Dileptons/2017/ProtonReconstruction_DoubleMu2017F_merge.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/eos/cms/store/group/phys_pps/ProtonRecoValidation/Dileptons/2017/ProtonReconstruction_DoubleMu2017F_merge.root:/ggll_aod");
      dir->GetObject("ntp1",tree);

   }
   Init(tree);
}

ntp1::~ntp1()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ntp1::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ntp1::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void ntp1::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   HLT_Name = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("LumiSection", &LumiSection, &b_LumiSection);
   fChain->SetBranchAddress("BX", &BX, &b_BX);
   fChain->SetBranchAddress("EventNum", &EventNum, &b_EventNum);
   fChain->SetBranchAddress("CrossingAngle", &CrossingAngle, &b_CrossingAngle);
   fChain->SetBranchAddress("nHLT", &nHLT, &b_nHLT);
   fChain->SetBranchAddress("HLT_Accept", HLT_Accept, &b_HLT_Accept);
   fChain->SetBranchAddress("HLT_Prescl", HLT_Prescl, &b_HLT_Prescl);
   fChain->SetBranchAddress("HLT_Name", &HLT_Name, &b_HLT_Name);
   fChain->SetBranchAddress("nMuonCand", &nMuonCand, &b_nMuonCand);
   fChain->SetBranchAddress("MuonCand_pt", MuonCand_pt, &b_MuonCand_pt);
   fChain->SetBranchAddress("MuonCand_eta", MuonCand_eta, &b_MuonCand_eta);
   fChain->SetBranchAddress("MuonCand_phi", MuonCand_phi, &b_MuonCand_phi);
   fChain->SetBranchAddress("MuonCand_e", MuonCand_e, &b_MuonCand_e);
   fChain->SetBranchAddress("MuonCand_charge", MuonCand_charge, &b_MuonCand_charge);
   fChain->SetBranchAddress("MuonCand_vtxx", MuonCand_vtxx, &b_MuonCand_vtxx);
   fChain->SetBranchAddress("MuonCand_vtxy", MuonCand_vtxy, &b_MuonCand_vtxy);
   fChain->SetBranchAddress("MuonCand_vtxz", MuonCand_vtxz, &b_MuonCand_vtxz);
   fChain->SetBranchAddress("MuonCand_dxy", MuonCand_dxy, &b_MuonCand_dxy);
   fChain->SetBranchAddress("MuonCand_nstatseg", MuonCand_nstatseg, &b_MuonCand_nstatseg);
   fChain->SetBranchAddress("MuonCand_ntrklayers", MuonCand_ntrklayers, &b_MuonCand_ntrklayers);
   fChain->SetBranchAddress("MuonCand_npxlhits", MuonCand_npxlhits, &b_MuonCand_npxlhits);
   fChain->SetBranchAddress("MuonCand_isglobal", MuonCand_isglobal, &b_MuonCand_isglobal);
   fChain->SetBranchAddress("MuonCand_istracker", MuonCand_istracker, &b_MuonCand_istracker);
   fChain->SetBranchAddress("MuonCand_isstandalone", MuonCand_isstandalone, &b_MuonCand_isstandalone);
   fChain->SetBranchAddress("MuonCand_ispfmuon", MuonCand_ispfmuon, &b_MuonCand_ispfmuon);
   fChain->SetBranchAddress("MuonCand_istight", MuonCand_istight, &b_MuonCand_istight);
   fChain->SetBranchAddress("MuonCandTrack_nmuchits", MuonCandTrack_nmuchits, &b_MuonCandTrack_nmuchits);
   fChain->SetBranchAddress("MuonCandTrack_chisq", MuonCandTrack_chisq, &b_MuonCandTrack_chisq);
   fChain->SetBranchAddress("MuonCand_innerTrackPt", MuonCand_innerTrackPt, &b_MuonCand_innerTrackPt);
   fChain->SetBranchAddress("MuonCand_innerTrackEta", MuonCand_innerTrackEta, &b_MuonCand_innerTrackEta);
   fChain->SetBranchAddress("MuonCand_innerTrackPhi", MuonCand_innerTrackPhi, &b_MuonCand_innerTrackPhi);
   fChain->SetBranchAddress("MuonCand_innerTrackVtxz", MuonCand_innerTrackVtxz, &b_MuonCand_innerTrackVtxz);
   fChain->SetBranchAddress("nPrimVertexCand", &nPrimVertexCand, &b_nPrimVertexCand);
   fChain->SetBranchAddress("PrimVertexCand_x", PrimVertexCand_x, &b_PrimVertexCand_x);
   fChain->SetBranchAddress("PrimVertexCand_y", PrimVertexCand_y, &b_PrimVertexCand_y);
   fChain->SetBranchAddress("PrimVertexCand_z", PrimVertexCand_z, &b_PrimVertexCand_z);
   fChain->SetBranchAddress("PrimVertexCand_chi2", PrimVertexCand_chi2, &b_PrimVertexCand_chi2);
   fChain->SetBranchAddress("PrimVertexCand_ndof", PrimVertexCand_ndof, &b_PrimVertexCand_ndof);
   fChain->SetBranchAddress("PrimVertexCand_tracks", PrimVertexCand_tracks, &b_PrimVertexCand_tracks);
   fChain->SetBranchAddress("nPair", &nPair, &b_nPair);
   fChain->SetBranchAddress("Pair_lepton1", Pair_lepton1, &b_Pair_lepton1);
   fChain->SetBranchAddress("Pair_lepton2", Pair_lepton2, &b_Pair_lepton2);
   fChain->SetBranchAddress("Pair_mass", Pair_mass, &b_Pair_mass);
   fChain->SetBranchAddress("Pair_pt", Pair_pt, &b_Pair_pt);
   fChain->SetBranchAddress("Pair_eta", Pair_eta, &b_Pair_eta);
   fChain->SetBranchAddress("Pair_phi", Pair_phi, &b_Pair_phi);
   fChain->SetBranchAddress("Pair_dpt", Pair_dpt, &b_Pair_dpt);
   fChain->SetBranchAddress("Pair_dphi", Pair_dphi, &b_Pair_dphi);
   fChain->SetBranchAddress("Pair_3Dangle", Pair_3Dangle, &b_Pair_3Dangle);
   fChain->SetBranchAddress("Pair_extratracks0p5mm", Pair_extratracks0p5mm, &b_Pair_extratracks0p5mm);
   fChain->SetBranchAddress("Pair_extratracks1mm", Pair_extratracks1mm, &b_Pair_extratracks1mm);
   fChain->SetBranchAddress("Pair_extratracks2mm", Pair_extratracks2mm, &b_Pair_extratracks2mm);
   fChain->SetBranchAddress("Pair_extratracks3mm", Pair_extratracks3mm, &b_Pair_extratracks3mm);
   fChain->SetBranchAddress("Pair_extratracks4mm", Pair_extratracks4mm, &b_Pair_extratracks4mm);
   fChain->SetBranchAddress("Pair_extratracks5mm", Pair_extratracks5mm, &b_Pair_extratracks5mm);
   fChain->SetBranchAddress("Pair_extratracks1cm", Pair_extratracks1cm, &b_Pair_extratracks1cm);
   fChain->SetBranchAddress("Pair_extratracks2cm", Pair_extratracks2cm, &b_Pair_extratracks2cm);
   fChain->SetBranchAddress("Pair_extratracks3cm", Pair_extratracks3cm, &b_Pair_extratracks3cm);
   fChain->SetBranchAddress("Pair_extratracks4cm", Pair_extratracks4cm, &b_Pair_extratracks4cm);
   fChain->SetBranchAddress("Pair_extratracks5cm", Pair_extratracks5cm, &b_Pair_extratracks5cm);
   fChain->SetBranchAddress("Pair_extratracks10cm", Pair_extratracks10cm, &b_Pair_extratracks10cm);
   fChain->SetBranchAddress("KalmanVertexCand_x", KalmanVertexCand_x, &b_KalmanVertexCand_x);
   fChain->SetBranchAddress("KalmanVertexCand_y", KalmanVertexCand_y, &b_KalmanVertexCand_y);
   fChain->SetBranchAddress("KalmanVertexCand_z", KalmanVertexCand_z, &b_KalmanVertexCand_z);
   fChain->SetBranchAddress("nLocalProtCand", &nLocalProtCand, &b_nLocalProtCand);
   fChain->SetBranchAddress("LocalProtCand_x", LocalProtCand_x, &b_LocalProtCand_x);
   fChain->SetBranchAddress("LocalProtCand_y", LocalProtCand_y, &b_LocalProtCand_y);
   fChain->SetBranchAddress("LocalProtCand_t", LocalProtCand_t, &b_LocalProtCand_t);
   fChain->SetBranchAddress("LocalProtCand_xSigma", LocalProtCand_xSigma, &b_LocalProtCand_xSigma);
   fChain->SetBranchAddress("LocalProtCand_ySigma", LocalProtCand_ySigma, &b_LocalProtCand_ySigma);
   fChain->SetBranchAddress("LocalProtCand_tSigma", LocalProtCand_tSigma, &b_LocalProtCand_tSigma);
   fChain->SetBranchAddress("LocalProtCand_arm", LocalProtCand_arm, &b_LocalProtCand_arm);
   fChain->SetBranchAddress("LocalProtCand_station", LocalProtCand_station, &b_LocalProtCand_station);
   fChain->SetBranchAddress("LocalProtCand_pot", LocalProtCand_pot, &b_LocalProtCand_pot);
   fChain->SetBranchAddress("LocalProtCand_rpid", LocalProtCand_rpid, &b_LocalProtCand_rpid);
   fChain->SetBranchAddress("nRecoProtCand", &nRecoProtCand, &b_nRecoProtCand);
   fChain->SetBranchAddress("ProtCand_xi", ProtCand_xi, &b_ProtCand_xi);
   fChain->SetBranchAddress("ProtCand_t", ProtCand_t, &b_ProtCand_t);
   fChain->SetBranchAddress("ProtCand_ThX", ProtCand_ThX, &b_ProtCand_ThX);
   fChain->SetBranchAddress("ProtCand_ThY", ProtCand_ThY, &b_ProtCand_ThY);
   fChain->SetBranchAddress("ProtCand_rpid", ProtCand_rpid, &b_ProtCand_rpid);
   fChain->SetBranchAddress("ProtCand_arm", ProtCand_arm, &b_ProtCand_arm);
   fChain->SetBranchAddress("ProtCand_ismultirp", ProtCand_ismultirp, &b_ProtCand_ismultirp);
   fChain->SetBranchAddress("nQualityExtraTrack", &nQualityExtraTrack, &b_nQualityExtraTrack);
   fChain->SetBranchAddress("ClosestExtraTrack_vtxdxyz", ClosestExtraTrack_vtxdxyz, &b_ClosestExtraTrack_vtxdxyz);
   fChain->SetBranchAddress("ClosestExtraTrack_id", ClosestExtraTrack_id, &b_ClosestExtraTrack_id);
   fChain->SetBranchAddress("ClosestHighPurityExtraTrack_vtxdxyz", ClosestHighPurityExtraTrack_vtxdxyz, &b_ClosestHighPurityExtraTrack_vtxdxyz);
   fChain->SetBranchAddress("ClosestHighPurityExtraTrack_id", ClosestHighPurityExtraTrack_id, &b_ClosestHighPurityExtraTrack_id);
   fChain->SetBranchAddress("Weight", &Weight, &b_Weight);
   fChain->SetBranchAddress("PUWeightTrue", &PUWeightTrue, &b_PUWeightTrue);
   Notify();
}

Bool_t ntp1::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ntp1::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ntp1::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ntp1_cxx
