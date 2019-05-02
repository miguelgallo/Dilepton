#define ntp1_MC_cxx
#include "ntp1_MC.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TSpline.h>
#include <math.h>

#define PI 3.14159265358979323846

#include <iostream>
#include <vector>
#include <set>
using namespace std;

void ntp1_MC::Loop()
{
//   In a ROOT session, you can do:
//      root> .L ntp1_MC.C
//      root> ntp1_MC t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

	TH1D *mumu_mass = new TH1D("mumu_mass", "#mu^{+}#mu^{-} Mass", 100, 0, 1000);
	TH1D *mumu_pt = new TH1D("mumu_pt", "#mu^{+}#mu^{-} p_{T}", 100, 0, 200);
	TH1D *mumu_y = new TH1D("mumu_y", "#mu^{+}#mu^{-} y", 100, -3, 3);
	TH2D *xi_left = new TH2D("xi_left", "#xi Left Correlation", 100, 0, 0.5, 100, 0, 0.5); 
	TH2D *xi_left_wrong = new TH2D("xi_left_wrong", "#xi Left Correlation", 100, 0, 0.5, 100, 0, 0.5); 
	TH2D *xi_right = new TH2D("x_right", "#xi Right Correlation", 100, 0, 0.5, 100, 0, 0.5);
	TH2D *xi_right_wrong = new TH2D("x_right_wrong", "#xi Right Correlation", 100, 0, 0.5, 100, 0, 0.5);

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
  
		double a = 1 - fabs(Pair_dphi[0])/PI;
		double xi_pair_left = (1./13000) * ( MuonCand_pt[0]*exp( MuonCand_eta[0] ) + MuonCand_pt[1]*exp( MuonCand_eta[1] ) );
      double xi_pair_right = (1./13000) * ( MuonCand_pt[0]*exp( -MuonCand_eta[0] ) + MuonCand_pt[1]*exp( -MuonCand_eta[1] ) );
		double xi_diff_left = (xi_pair_left - ProtCand_xi[0])/(xi_pair_left);
      double xi_diff_right = (xi_pair_right - ProtCand_xi[0])/(xi_pair_right);

		TLorentzVector v;

		if(MuonCand_pt[0] > 50 && MuonCand_pt[1] > 50){
			if(MuonCand_charge[0] * MuonCand_charge[1] < 0){
				if(MuonCand_istight[0] && MuonCand_istight[1]){
        			if(Pair_mass[0] > 110.){
              		if(fabs(KalmanVertexCand_z[0] < 15.)){
                  	if(ClosestExtraTrack_vtxdxyz[0] > 0.05){
								if(a < 0.009){
					
									//cout << *ProtCand_arm << endl;
									
									v.SetPtEtaPhiM(Pair_pt[0],Pair_eta[0],Pair_phi[0],Pair_mass[0]);		

									mumu_mass->Fill(Pair_mass[0]);
									mumu_pt->Fill(Pair_pt[0]);
									mumu_y->Fill(v.Rapidity());
									for ( int idx = 0; idx < nRecoProtCand; ++idx ) {
										if ((ProtCand_rpid[idx] == 3) || (ProtCand_rpid[idx] == 23)) {
											if (ProtCand_ismultirp[idx] == 0) {
												xi_left->Fill(ProtCand_xi[idx], xi_pair_left); 
												xi_left_wrong->Fill(ProtCand_xi[idx], xi_pair_right); 
											}		
										}
										if ((ProtCand_rpid[idx] == 103) || (ProtCand_rpid[idx] == 123)) {                           
            	                  if (ProtCand_ismultirp[idx] == 0) {
												xi_right->Fill(ProtCand_xi[idx], xi_pair_right);
												xi_right_wrong->Fill(ProtCand_xi[idx], xi_pair_left);
											}	
               	            }
                           }
								}
							}
						}
					}
				}
			}
		}	
	}

	TFile* f = new TFile("out_MC_Pt50_xangle150_divergence30.root", "RECREATE");
   mumu_mass->DrawCopy();
   mumu_pt->DrawCopy();
   mumu_y->DrawCopy();
   xi_left->DrawCopy();
   xi_right->DrawCopy();
   mumu_mass->Write();
   mumu_pt->Write();
   mumu_y->Write();
   xi_left->Write();
   xi_right->Write();
	
	f->Close();
}
	
int run() {
	ntp1_MC m;
 	m.Loop();

 	return 0;
 	}

