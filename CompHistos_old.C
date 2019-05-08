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
#include <TH1F.h>
#include <THStack.h>

#define PI 3.14159265358979323846

#include <iostream>
#include <vector>
#include <set>
using namespace std;

int CompHistos_old()
{
		  gStyle->SetOptStat(0);

		  //char* const hist[] = {"mumu_mass", "mumu_pt", "mumu_y", "xi_left", "xi_right"}; 

		  //size_t Nhist = sizeof(hist)/sizeof(hist[0]);

		  TFile *f_elastic_x130 = new TFile("MC_2017_LPAIRmumuElastic/out_MC_Pt50_xangle130.root");
		  TFile *f_elastic_x130_div30 = new TFile("MC_2017_LPAIRmumuElastic/out_MC_Pt50_xangle130_divergence30.root");
		  TFile *f_elastic_x140_div30 = new TFile("MC_2017_LPAIRmumuElastic/out_MC_Pt50_xangle140_divergence30.root");
		  TFile *f_elastic_x150_div30 = new TFile("MC_2017_LPAIRmumuElastic/out_MC_Pt50_xangle150_divergence30.root");
		  TFile *f_singlediss = new TFile("MC_2017_LPAIRmumuSingleDiss/out_MC_Pt50_xangle130.root");
		  TFile *f_data = new TFile("2017/out_data_run2017.root");

		  //for (int i = 0; i < Nhist; i++) {

		  double scale_factor_elastic_x130, scale_factor_elastic_x130_div30, scale_factor_elastic_x140_div30, scale_factor_elastic_x150_div30, scale_factor_singlediss;
		  double eff_trigger_elastic_x130 = 1;
		  double eff_trigger_elastic_x130_div30 = 1;
		  double eff_trigger_elastic_x140_div30 = 1;
		  double eff_trigger_elastic_x150_div30 = 1;
		  double eff_trigger_singlediss = 1;

		  TH1* h_elastic_x130 = 0;
		  TH1* h_elastic_x130_div30 = 0;
		  TH1* h_elastic_x140_div30 = 0;
		  TH1* h_elastic_x150_div30 = 0;
		  TH1* h_singlediss = 0;
		  TH1* h_data = 0;


		  f_elastic_x130->GetObject("mumu_mass", h_elastic_x130);
		  f_elastic_x130_div30->GetObject("mumu_mass", h_elastic_x130_div30);
		  f_elastic_x140_div30->GetObject("mumu_mass", h_elastic_x140_div30);
		  f_elastic_x150_div30->GetObject("mumu_mass", h_elastic_x150_div30);
		  f_singlediss->GetObject("mumu_mass", h_singlediss);
		  f_data->GetObject("mumu_mass", h_data);

        double n_events_h_elastic_x130 = 2472;
        double n_events_h_elastic_x130_div30 = 2459;
        double n_events_h_elastic_x140_div30 = 2299;
        double n_events_h_elastic_x150_div30 = 2139;
        double n_events_h_singlediss = 1532;

		  //xsec1 = eff_trigger1*37.022*0.017254036*1000*(h1->GetEntries()/2000) ; //
		  //xsec2 = eff_trigger2*37.022*0.025643100*1000*(h2->GetEntries()/2000) ;//
		  scale_factor_elastic_x130 = eff_trigger_elastic_x130*37.022*0.017254036*1000 / n_events_h_elastic_x130 ; //
		  scale_factor_elastic_x130_div30 = eff_trigger_elastic_x130_div30*37.022*0.017254036*1000 / n_events_h_elastic_x130_div30 ; //
		  scale_factor_elastic_x140_div30 = eff_trigger_elastic_x140_div30*37.022*0.017254036*1000 / n_events_h_elastic_x140_div30 ; //
		  scale_factor_elastic_x150_div30 = eff_trigger_elastic_x150_div30*37.022*0.017254036*1000 / n_events_h_elastic_x150_div30 ; //
		  scale_factor_singlediss = eff_trigger_singlediss*37.022*0.025643100*1000 / n_events_h_singlediss ; //

		  //	  h_WJ_PT_1->Add(h_WJ_PT_2);
		  //	  h_WJ_PT_1->Add(h_WJ_PT_3);
		  //	  h_WJ_PT_1->Add(h_WJ_PT_4);
		  //	  h_WJ_PT_1->Add(h_WJ_PT_5);

		  h_elastic_x130->SetLineColor(kRed);
		  h_elastic_x130_div30->SetLineColor(kRed);
		  h_elastic_x140_div30->SetLineColor(kRed);
		  h_elastic_x150_div30->SetLineColor(kRed);
		  h_singlediss->SetLineColor(kBlue);
		  h_elastic_x130->SetFillColor(kRed);
		  h_elastic_x130_div30->SetFillColor(kRed);
		  h_elastic_x140_div30->SetFillColor(kRed);
		  h_elastic_x150_div30->SetFillColor(kRed);
		  h_singlediss->SetFillColor(kBlue);
		  h_data->SetMarkerStyle(8);
		  h_data->Sumw2();

		  //h1->Scale(xsec1/h1->GetEntries());
		  //h2->Scale(xsec2/h2->GetEntries());
		  h_elastic_x130->Scale(scale_factor_elastic_x130);
		  h_elastic_x130_div30->Scale(scale_factor_elastic_x130_div30);
		  h_elastic_x140_div30->Scale(scale_factor_elastic_x140_div30);
		  h_elastic_x150_div30->Scale(scale_factor_elastic_x150_div30);
		  h_singlediss->Scale(scale_factor_singlediss);

        double integral_h_elastic_x130 = h_elastic_x130->Integral( 0, ( h_elastic_x130->GetNbinsX() + 1 ) ); // With under and overflow
        double integral_h_elastic_x130_div30 = h_elastic_x130_div30->Integral( 0, ( h_elastic_x130_div30->GetNbinsX() + 1 ) ); // With under and overflow
        double integral_h_elastic_x140_div30 = h_elastic_x140_div30->Integral( 0, ( h_elastic_x140_div30->GetNbinsX() + 1 ) ); // With under and overflow
        double integral_h_elastic_x150_div30 = h_elastic_x150_div30->Integral( 0, ( h_elastic_x150_div30->GetNbinsX() + 1 ) ); // With under and overflow
        double integral_h_singlediss = h_singlediss->Integral( 0, ( h_singlediss->GetNbinsX() + 1 ) ); // With under and overflow
        double sum_integral = integral_h_elastic_x130 + integral_h_elastic_x130_div30 + integral_h_elastic_x140_div30 + integral_h_elastic_x150_div30  + integral_h_singlediss;
        double scale_data = h_data->GetEntries() / sum_integral ;
        cout << "Scale factor data normalization: " << scale_data << endl;        

        // Clone MC histograms and scale
        
        TH1F *h_elastic_x130_clone = (TH1F*)h_elastic_x130->Clone("h_elastic_x130_clone");
        TH1F *h_elastic_x130_div30_clone = (TH1F*)h_elastic_x130_div30->Clone("h_elastic_x130_div30_clone");
        TH1F *h_elastic_x140_div30_clone = (TH1F*)h_elastic_x140_div30->Clone("h_elastic_x140_div30_clone");
        TH1F *h_elastic_x150_div30_clone = (TH1F*)h_elastic_x150_div30->Clone("h_elastic_x150_div30_clone");
        TH1F *h_singlediss_clone = (TH1F*)h_singlediss->Clone("h_singlediss_clone");

		  h_elastic_x130_clone->Scale(scale_data);
        h_elastic_x130_div30_clone->Scale(scale_data);
        h_elastic_x140_div30_clone->Scale(scale_data);
        h_elastic_x150_div30_clone->Scale(scale_data);
        h_singlediss_clone->Scale(scale_data);
		  
	     h_elastic_x130_clone->Add(h_elastic_x130_div30_clone);
		  h_elastic_x130_clone->Add(h_elastic_x140_div30_clone);
		  h_elastic_x130_clone->Add(h_elastic_x150_div30_clone);
	
		  THStack hs("hs","");
		  hs.Add(h_elastic_x130_clone);
		  hs.Add(h_singlediss_clone);

		  TFile* f4 = new TFile("2017/comp_data_MC_2017.root", "RECREATE");
		  hs.Draw("HIST");
		  //h_elastic_x130->Draw("HIST");
		  //h_elastic_x130_div30->Draw("HIST SAME");
		  //h_elastic_x140_div30->Draw("HIST SAME");
		  //h_elastic_x150_div30->Draw("HIST SAME");
		  //h_singlediss->Draw("HIST");
		  h_data->Draw("SAME");
		  hs.Write();
	     h_data->Write();
		  f4->Close();
		 // }
		  return 0;
}
