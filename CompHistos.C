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

int CompHistos()
{
		  gStyle->SetOptStat(0);

		  //char* const hist[] = {"mumu_mass", "mumu_pt", "mumu_y", "xi_left", "xi_right"}; 

		  //size_t Nhist = sizeof(hist)/sizeof(hist[0]);

		  //TFile *f_elastic_x130 = new TFile("MC_2017_LPAIRmumuElastic/out_MC_Pt50_xangle130.root");
		  TFile *f_elastic_x130_div30 = new TFile("MC_2017_LPAIRmumuElastic/out_MC_Pt50_xangle130_divergence30.root");
		  TFile *f_elastic_x140_div30 = new TFile("MC_2017_LPAIRmumuElastic/out_MC_Pt50_xangle140_divergence30.root");
		  TFile *f_elastic_x150_div30 = new TFile("MC_2017_LPAIRmumuElastic/out_MC_Pt50_xangle150_divergence30.root");
		  TFile *f_singlediss = new TFile("MC_2017_LPAIRmumuSingleDiss/out_MC_Pt50_xangle130.root");
		  TFile *f_data = new TFile("2017/out_data_run2017.root");

		  //for (int i = 0; i < Nhist; i++) {

		  double scale_factor_elastic_x130, scale_factor_elastic_x130_div30, scale_factor_elastic_x140_div30, scale_factor_elastic_x150_div30, scale_factor_singlediss;
		  //double eff_trigger_elastic_x130 = 1;
		  double eff_trigger_elastic_x130_div30 = 1;
		  double eff_trigger_elastic_x140_div30 = 1;
		  double eff_trigger_elastic_x150_div30 = 1;
		  double eff_trigger_singlediss = 1;

		  //TH1* h_elastic_x130_mass = 0;					
		  TH1* h_elastic_x130_div30_mass = 0;
		  TH1* h_elastic_x140_div30_mass = 0;
		  TH1* h_elastic_x150_div30_mass = 0;
		  TH1* h_singlediss_mass = 0;
		  TH1* h_data_mass = 0;
		  //TH1* h_elastic_x130_pt = 0;					
		  TH1* h_elastic_x130_div30_pt = 0;
		  TH1* h_elastic_x140_div30_pt = 0;
		  TH1* h_elastic_x150_div30_pt = 0;
		  TH1* h_singlediss_pt = 0;
		  TH1* h_data_pt = 0;
		  //TH1* h_elastic_x130_y = 0;					
		  TH1* h_elastic_x130_div30_y = 0;
		  TH1* h_elastic_x140_div30_y = 0;
		  TH1* h_elastic_x150_div30_y = 0;
		  TH1* h_singlediss_y = 0;
		  TH1* h_data_y = 0;
		  //TH1* h_elastic_x130_xil = 0;					
		  TH1* h_elastic_x130_div30_xil = 0;
		  TH1* h_elastic_x140_div30_xil = 0;
		  TH1* h_elastic_x150_div30_xil = 0;
		  TH1* h_singlediss_xil = 0;
		  TH1* h_data_xil = 0;
		  //TH1* h_elastic_x130_xir = 0;					
		  TH1* h_elastic_x130_div30_xir = 0;
		  TH1* h_elastic_x140_div30_xir = 0;
		  TH1* h_elastic_x150_div30_xir = 0;
		  TH1* h_singlediss_xir = 0;
		  TH1* h_data_xir = 0;


		  //f_elastic_x130->GetObject("mumu_mass", h_elastic_x130_mass);
		  f_elastic_x130_div30->GetObject("mumu_mass", h_elastic_x130_div30_mass);
		  f_elastic_x140_div30->GetObject("mumu_mass", h_elastic_x140_div30_mass);
		  f_elastic_x150_div30->GetObject("mumu_mass", h_elastic_x150_div30_mass);
		  f_singlediss->GetObject("mumu_mass", h_singlediss_mass);
		  f_data->GetObject("mumu_mass", h_data_mass);
		  //f_elastic_x130->GetObject("mumu_pt", h_elastic_x130_pt);
		  f_elastic_x130_div30->GetObject("mumu_pt", h_elastic_x130_div30_pt);
		  f_elastic_x140_div30->GetObject("mumu_pt", h_elastic_x140_div30_pt);
		  f_elastic_x150_div30->GetObject("mumu_pt", h_elastic_x150_div30_pt);
		  f_singlediss->GetObject("mumu_pt", h_singlediss_pt);
		  f_data->GetObject("mumu_pt", h_data_pt);
		  //f_elastic_x130->GetObject("mumu_y", h_elastic_x130_y);
		  f_elastic_x130_div30->GetObject("mumu_y", h_elastic_x130_div30_y);
		  f_elastic_x140_div30->GetObject("mumu_y", h_elastic_x140_div30_y);
		  f_elastic_x150_div30->GetObject("mumu_y", h_elastic_x150_div30_y);
		  f_singlediss->GetObject("mumu_y", h_singlediss_y);
		  f_data->GetObject("mumu_y", h_data_y);
		  //f_elastic_x130->GetObject("xi_left", h_elastic_x130_xil);
		  f_elastic_x130_div30->GetObject("xi_left", h_elastic_x130_div30_xil);
		  f_elastic_x140_div30->GetObject("xi_left", h_elastic_x140_div30_xil);
		  f_elastic_x150_div30->GetObject("xi_left", h_elastic_x150_div30_xil);
		  f_singlediss->GetObject("xi_left", h_singlediss_xil);
		  f_data->GetObject("xi_left", h_data_xil);
		  //f_elastic_x130->GetObject("x_right", h_elastic_x130_xir);
		  f_elastic_x130_div30->GetObject("x_right", h_elastic_x130_div30_xir);
		  f_elastic_x140_div30->GetObject("x_right", h_elastic_x140_div30_xir);
		  f_elastic_x150_div30->GetObject("x_right", h_elastic_x150_div30_xir);
		  f_singlediss->GetObject("x_right", h_singlediss_xir);
		  f_data->GetObject("x_right", h_data_xir);

        //double n_events_h_elastic_x130 = 2472;
        double n_events_h_elastic_x130_div30 = 2459;
        double n_events_h_elastic_x140_div30 = 2299;
        double n_events_h_elastic_x150_div30 = 2139;
		  double n_events_h_elastic = n_events_h_elastic_x130_div30 + n_events_h_elastic_x140_div30 + n_events_h_elastic_x150_div30;
        double n_events_h_singlediss = 1532;

		  //xsec1 = eff_trigger1*37.022*0.017254036*1000*(h1->GetEntries()/2000) ; //
		  //xsec2 = eff_trigger2*37.022*0.025643100*1000*(h2->GetEntries()/2000) ;//
		  //scale_factor_elastic_x130 = eff_trigger_elastic_x130*37.022*0.017254036*1000 / n_events_h_elastic_x130 ; //
		  scale_factor_elastic_x130_div30 = eff_trigger_elastic_x130_div30*37.022*0.017254036*1000 / n_events_h_elastic ; //
		  scale_factor_elastic_x140_div30 = eff_trigger_elastic_x140_div30*37.022*0.017254036*1000 / n_events_h_elastic ; //
		  scale_factor_elastic_x150_div30 = eff_trigger_elastic_x150_div30*37.022*0.017254036*1000 / n_events_h_elastic ; //
		  scale_factor_singlediss = eff_trigger_singlediss*37.022*0.025643100*1000 / n_events_h_singlediss ; //

		  //	  h_WJ_PT_1->Add(h_WJ_PT_2);
		  //	  h_WJ_PT_1->Add(h_WJ_PT_3);
		  //	  h_WJ_PT_1->Add(h_WJ_PT_4);
		  //	  h_WJ_PT_1->Add(h_WJ_PT_5);

		  //h_elastic_x130_mass->SetLineColor(kRed);
		  h_elastic_x130_div30_mass->SetLineColor(kRed);
		  h_elastic_x140_div30_mass->SetLineColor(kRed);
		  h_elastic_x150_div30_mass->SetLineColor(kRed);
		  h_singlediss_mass->SetLineColor(kBlue);
		  //h_elastic_x130_mass->SetFillColor(kRed);
		  h_elastic_x130_div30_mass->SetFillColor(kRed);
		  h_elastic_x140_div30_mass->SetFillColor(kRed);
		  h_elastic_x150_div30_mass->SetFillColor(kRed);
		  h_singlediss_mass->SetFillColor(kBlue);
		  h_data_mass->SetMarkerStyle(8);
		  h_data_mass->Sumw2();
		  //h_elastic_x130_pt->SetLineColor(kRed);
		  h_elastic_x130_div30_pt->SetLineColor(kRed);
		  h_elastic_x140_div30_pt->SetLineColor(kRed);
		  h_elastic_x150_div30_pt->SetLineColor(kRed);
		  h_singlediss_pt->SetLineColor(kBlue);
		  //h_elastic_x130_pt->SetFillColor(kRed);
		  h_elastic_x130_div30_pt->SetFillColor(kRed);
		  h_elastic_x140_div30_pt->SetFillColor(kRed);
		  h_elastic_x150_div30_pt->SetFillColor(kRed);
		  h_singlediss_pt->SetFillColor(kBlue);
		  h_data_pt->SetMarkerStyle(8);
		  h_data_pt->Sumw2();
		  //h_elastic_x130_y->SetLineColor(kRed);
		  h_elastic_x130_div30_y->SetLineColor(kRed);
		  h_elastic_x140_div30_y->SetLineColor(kRed);
		  h_elastic_x150_div30_y->SetLineColor(kRed);
		  h_singlediss_y->SetLineColor(kBlue);
		  //h_elastic_x130_y->SetFillColor(kRed);
		  h_elastic_x130_div30_y->SetFillColor(kRed);
		  h_elastic_x140_div30_y->SetFillColor(kRed);
		  h_elastic_x150_div30_y->SetFillColor(kRed);
		  h_singlediss_y->SetFillColor(kBlue);
		  h_data_y->SetMarkerStyle(8);
		  h_data_y->Sumw2();
		  //h_elastic_x130_xil->SetLineColor(kRed);
		  h_elastic_x130_div30_xil->SetLineColor(kRed);
		  h_elastic_x140_div30_xil->SetLineColor(kRed);
		  h_elastic_x150_div30_xil->SetLineColor(kRed);
		  h_singlediss_xil->SetLineColor(kBlue);
		  //h_elastic_x130_xil->SetFillColor(kRed);
		  h_elastic_x130_div30_xil->SetFillColor(kRed);
		  h_elastic_x140_div30_xil->SetFillColor(kRed);
		  h_elastic_x150_div30_xil->SetFillColor(kRed);
		  h_singlediss_xil->SetFillColor(kBlue);
		  h_data_xil->SetMarkerStyle(8);
		  h_data_xil->Sumw2();
		  //h_elastic_x130_xir->SetLineColor(kRed);
		  h_elastic_x130_div30_xir->SetLineColor(kRed);
		  h_elastic_x140_div30_xir->SetLineColor(kRed);
		  h_elastic_x150_div30_xir->SetLineColor(kRed);
		  h_singlediss_xir->SetLineColor(kBlue);
		  //h_elastic_x130_xir->SetFillColor(kRed);
		  h_elastic_x130_div30_xir->SetFillColor(kRed);
		  h_elastic_x140_div30_xir->SetFillColor(kRed);
		  h_elastic_x150_div30_xir->SetFillColor(kRed);
		  h_singlediss_xir->SetFillColor(kBlue);
		  h_data_xir->SetMarkerStyle(8);
		  h_data_xir->Sumw2();

		  //h1->Scale(xsec1/h1->GetEntries());
		  //h2->Scale(xsec2/h2->GetEntries());
		  //h_elastic_x130_mass->Scale(scale_factor_elastic_x130);
		  h_elastic_x130_div30_mass->Scale(scale_factor_elastic_x130_div30);
		  h_elastic_x140_div30_mass->Scale(scale_factor_elastic_x140_div30);
		  h_elastic_x150_div30_mass->Scale(scale_factor_elastic_x150_div30);
		  h_singlediss_mass->Scale(scale_factor_singlediss);
		  //h_elastic_x130_pt->Scale(scale_factor_elastic_x130);
		  h_elastic_x130_div30_pt->Scale(scale_factor_elastic_x130_div30);
		  h_elastic_x140_div30_pt->Scale(scale_factor_elastic_x140_div30);
		  h_elastic_x150_div30_pt->Scale(scale_factor_elastic_x150_div30);
		  h_singlediss_pt->Scale(scale_factor_singlediss);
		  //h_elastic_x130_y->Scale(scale_factor_elastic_x130);
		  h_elastic_x130_div30_y->Scale(scale_factor_elastic_x130_div30);
		  h_elastic_x140_div30_y->Scale(scale_factor_elastic_x140_div30);
		  h_elastic_x150_div30_y->Scale(scale_factor_elastic_x150_div30);
		  h_singlediss_y->Scale(scale_factor_singlediss);
		  //h_elastic_x130_xil->Scale(scale_factor_elastic_x130);
		  h_elastic_x130_div30_xil->Scale(scale_factor_elastic_x130_div30);
		  h_elastic_x140_div30_xil->Scale(scale_factor_elastic_x140_div30);
		  h_elastic_x150_div30_xil->Scale(scale_factor_elastic_x150_div30);
		  h_singlediss_xil->Scale(scale_factor_singlediss);
		  //h_elastic_x130_xir->Scale(scale_factor_elastic_x130);
		  h_elastic_x130_div30_xir->Scale(scale_factor_elastic_x130_div30);
		  h_elastic_x140_div30_xir->Scale(scale_factor_elastic_x140_div30);
		  h_elastic_x150_div30_xir->Scale(scale_factor_elastic_x150_div30);
		  h_singlediss_xir->Scale(scale_factor_singlediss);

        //double integral_h_elastic_x130_mass = h_elastic_x130_mass->Integral( 0, ( h_elastic_x130_mass->GetNbinsX() + 1 ) ); // With under and overflow
        double integral_h_elastic_x130_div30_mass = h_elastic_x130_div30_mass->Integral( 0, ( h_elastic_x130_div30_mass->GetNbinsX() + 1 ) ); // With under and overflow
        double integral_h_elastic_x140_div30_mass = h_elastic_x140_div30_mass->Integral( 0, ( h_elastic_x140_div30_mass->GetNbinsX() + 1 ) ); // With under and overflow
        double integral_h_elastic_x150_div30_mass = h_elastic_x150_div30_mass->Integral( 0, ( h_elastic_x150_div30_mass->GetNbinsX() + 1 ) ); // With under and overflow
        double integral_h_singlediss_mass = h_singlediss_mass->Integral( 0, ( h_singlediss_mass->GetNbinsX() + 1 ) ); // With under and overflow
        double sum_integral_mass = integral_h_elastic_x130_div30_mass + integral_h_elastic_x140_div30_mass + integral_h_elastic_x150_div30_mass  + integral_h_singlediss_mass;
        double scale_data_mass = h_data_mass->GetEntries() / sum_integral_mass ;
        cout << "Scale factor mass data normalization: " << scale_data_mass << endl;        
        //double integral_h_elastic_x130_pt = h_elastic_x130_pt->Integral( 0, ( h_elastic_x130_pt->GetNbinsX() + 1 ) ); // With under and overflow
        double integral_h_elastic_x130_div30_pt = h_elastic_x130_div30_pt->Integral( 0, ( h_elastic_x130_div30_pt->GetNbinsX() + 1 ) ); // With under and overflow
        double integral_h_elastic_x140_div30_pt = h_elastic_x140_div30_pt->Integral( 0, ( h_elastic_x140_div30_pt->GetNbinsX() + 1 ) ); // With under and overflow
        double integral_h_elastic_x150_div30_pt = h_elastic_x150_div30_pt->Integral( 0, ( h_elastic_x150_div30_pt->GetNbinsX() + 1 ) ); // With under and overflow
        double integral_h_singlediss_pt = h_singlediss_pt->Integral( 0, ( h_singlediss_pt->GetNbinsX() + 1 ) ); // With under and overflow
        double sum_integral_pt = integral_h_elastic_x130_div30_pt + integral_h_elastic_x140_div30_pt + integral_h_elastic_x150_div30_pt  + integral_h_singlediss_pt;
        double scale_data_pt = h_data_pt->GetEntries() / sum_integral_pt ;
        cout << "Scale factor pt data normalization: " << scale_data_pt << endl;        
        //double integral_h_elastic_x130_y = h_elastic_x130_y->Integral( 0, ( h_elastic_x130_y->GetNbinsX() + 1 ) ); // With under and overflow
        double integral_h_elastic_x130_div30_y = h_elastic_x130_div30_y->Integral( 0, ( h_elastic_x130_div30_y->GetNbinsX() + 1 ) ); // With under and overflow
        double integral_h_elastic_x140_div30_y = h_elastic_x140_div30_y->Integral( 0, ( h_elastic_x140_div30_y->GetNbinsX() + 1 ) ); // With under and overflow
        double integral_h_elastic_x150_div30_y = h_elastic_x150_div30_y->Integral( 0, ( h_elastic_x150_div30_y->GetNbinsX() + 1 ) ); // With under and overflow
        double integral_h_singlediss_y = h_singlediss_y->Integral( 0, ( h_singlediss_y->GetNbinsX() + 1 ) ); // With under and overflow
        double sum_integral_y = integral_h_elastic_x130_div30_y + integral_h_elastic_x140_div30_y + integral_h_elastic_x150_div30_y  + integral_h_singlediss_y;
        double scale_data_y = h_data_y->GetEntries() / sum_integral_y ;
        cout << "Scale factor y data normalization: " << scale_data_y << endl;        
        //double integral_h_elastic_x130_xil = h_elastic_x130_xil->Integral( 0, ( h_elastic_x130_xil->GetNbinsX() + 1 ) ); // With under and overflow
        double integral_h_elastic_x130_div30_xil = h_elastic_x130_div30_xil->Integral( 0, ( h_elastic_x130_div30_xil->GetNbinsX() + 1 ) ); // With under and overflow
        double integral_h_elastic_x140_div30_xil = h_elastic_x140_div30_xil->Integral( 0, ( h_elastic_x140_div30_xil->GetNbinsX() + 1 ) ); // With under and overflow
        double integral_h_elastic_x150_div30_xil = h_elastic_x150_div30_xil->Integral( 0, ( h_elastic_x150_div30_xil->GetNbinsX() + 1 ) ); // With under and overflow
        double integral_h_singlediss_xil = h_singlediss_xil->Integral( 0, ( h_singlediss_xil->GetNbinsX() + 1 ) ); // With under and overflow
        double sum_integral_xil = integral_h_elastic_x130_div30_xil + integral_h_elastic_x140_div30_xil + integral_h_elastic_x150_div30_xil  + integral_h_singlediss_xil;
        double scale_data_xil = h_data_xil->GetEntries() / sum_integral_xil ;
        cout << "Scale factor xil data normalization: " << scale_data_xil << endl;        
        //double integral_h_elastic_x130_xir = h_elastic_x130_xir->Integral( 0, ( h_elastic_x130_xir->GetNbinsX() + 1 ) ); // With under and overflow
        double integral_h_elastic_x130_div30_xir = h_elastic_x130_div30_xir->Integral( 0, ( h_elastic_x130_div30_xir->GetNbinsX() + 1 ) ); // With under and overflow
        double integral_h_elastic_x140_div30_xir = h_elastic_x140_div30_xir->Integral( 0, ( h_elastic_x140_div30_xir->GetNbinsX() + 1 ) ); // With under and overflow
        double integral_h_elastic_x150_div30_xir = h_elastic_x150_div30_xir->Integral( 0, ( h_elastic_x150_div30_xir->GetNbinsX() + 1 ) ); // With under and overflow
        double integral_h_singlediss_xir = h_singlediss_xir->Integral( 0, ( h_singlediss_xir->GetNbinsX() + 1 ) ); // With under and overflow
        double sum_integral_xir = integral_h_elastic_x130_div30_xir + integral_h_elastic_x140_div30_xir + integral_h_elastic_x150_div30_xir  + integral_h_singlediss_xir;
        double scale_data_xir = h_data_xir->GetEntries() / sum_integral_xir ;
        cout << "Scale factor xir data normalization: " << scale_data_xir << endl;        

        // Clone MC histograms and scale
        
        //TH1F *h_elastic_x130_clone_mass = (TH1F*)h_elastic_x130_mass->Clone("h_elastic_x130_clone_mass");
        TH1F *h_elastic_x130_div30_clone_mass = (TH1F*)h_elastic_x130_div30_mass->Clone("h_elastic_x130_div30_clone_mass");
        TH1F *h_elastic_x140_div30_clone_mass = (TH1F*)h_elastic_x140_div30_mass->Clone("h_elastic_x140_div30_clone_mass");
        TH1F *h_elastic_x150_div30_clone_mass = (TH1F*)h_elastic_x150_div30_mass->Clone("h_elastic_x150_div30_clone_mass");
        TH1F *h_singlediss_clone_mass = (TH1F*)h_singlediss_mass->Clone("h_singlediss_clone_mass");
        //TH1F *h_elastic_x130_clone_pt = (TH1F*)h_elastic_x130_pt->Clone("h_elastic_x130_clone_pt");
        TH1F *h_elastic_x130_div30_clone_pt = (TH1F*)h_elastic_x130_div30_pt->Clone("h_elastic_x130_div30_clone_pt");
        TH1F *h_elastic_x140_div30_clone_pt = (TH1F*)h_elastic_x140_div30_pt->Clone("h_elastic_x140_div30_clone_pt");
        TH1F *h_elastic_x150_div30_clone_pt = (TH1F*)h_elastic_x150_div30_pt->Clone("h_elastic_x150_div30_clone_pt");
        TH1F *h_singlediss_clone_pt = (TH1F*)h_singlediss_pt->Clone("h_singlediss_clone_pt");
        //TH1F *h_elastic_x130_clone_y = (TH1F*)h_elastic_x130_y->Clone("h_elastic_x130_clone_y");
        TH1F *h_elastic_x130_div30_clone_y = (TH1F*)h_elastic_x130_div30_y->Clone("h_elastic_x130_div30_clone_y");
        TH1F *h_elastic_x140_div30_clone_y = (TH1F*)h_elastic_x140_div30_y->Clone("h_elastic_x140_div30_clone_y");
        TH1F *h_elastic_x150_div30_clone_y = (TH1F*)h_elastic_x150_div30_y->Clone("h_elastic_x150_div30_clone_y");
        TH1F *h_singlediss_clone_y = (TH1F*)h_singlediss_y->Clone("h_singlediss_clone_y");
        //TH1F *h_elastic_x130_clone_xil = (TH1F*)h_elastic_x130_xil->Clone("h_elastic_x130_clone_xil");
        TH1F *h_elastic_x130_div30_clone_xil = (TH1F*)h_elastic_x130_div30_xil->Clone("h_elastic_x130_div30_clone_xil");
        TH1F *h_elastic_x140_div30_clone_xil = (TH1F*)h_elastic_x140_div30_xil->Clone("h_elastic_x140_div30_clone_xil");
        TH1F *h_elastic_x150_div30_clone_xil = (TH1F*)h_elastic_x150_div30_xil->Clone("h_elastic_x150_div30_clone_xil");
        TH1F *h_singlediss_clone_xil = (TH1F*)h_singlediss_xil->Clone("h_singlediss_clone_xil");
        //TH1F *h_elastic_x130_clone_xir = (TH1F*)h_elastic_x130_xir->Clone("h_elastic_x130_clone_xir");
        TH1F *h_elastic_x130_div30_clone_xir = (TH1F*)h_elastic_x130_div30_xir->Clone("h_elastic_x130_div30_clone_xir");
        TH1F *h_elastic_x140_div30_clone_xir = (TH1F*)h_elastic_x140_div30_xir->Clone("h_elastic_x140_div30_clone_xir");
        TH1F *h_elastic_x150_div30_clone_xir = (TH1F*)h_elastic_x150_div30_xir->Clone("h_elastic_x150_div30_clone_xir");
        TH1F *h_singlediss_clone_xir = (TH1F*)h_singlediss_xir->Clone("h_singlediss_clone_xir");

		  //h_elastic_x130_clone_mass->Scale(scale_data_mass);
        h_elastic_x130_div30_clone_mass->Scale(scale_data_mass);
        h_elastic_x140_div30_clone_mass->Scale(scale_data_mass);
        h_elastic_x150_div30_clone_mass->Scale(scale_data_mass);
        h_singlediss_clone_mass->Scale(scale_data_mass);
		  //h_elastic_x130_clone_pt->Scale(scale_data_pt);
        h_elastic_x130_div30_clone_pt->Scale(scale_data_pt);
        h_elastic_x140_div30_clone_pt->Scale(scale_data_pt);
        h_elastic_x150_div30_clone_pt->Scale(scale_data_pt);
        h_singlediss_clone_pt->Scale(scale_data_pt);
		  //h_elastic_x130_clone_y->Scale(scale_data_y);
        h_elastic_x130_div30_clone_y->Scale(scale_data_y);
        h_elastic_x140_div30_clone_y->Scale(scale_data_y);
        h_elastic_x150_div30_clone_y->Scale(scale_data_y);
        h_singlediss_clone_y->Scale(scale_data_y);
		  //h_elastic_x130_clone_xil->Scale(scale_data_xil);
        h_elastic_x130_div30_clone_xil->Scale(scale_data_xil);
        h_elastic_x140_div30_clone_xil->Scale(scale_data_xil);
        h_elastic_x150_div30_clone_xil->Scale(scale_data_xil);
        h_singlediss_clone_xil->Scale(scale_data_xil);
		  //h_elastic_x130_clone_xir->Scale(scale_data_xir);
        h_elastic_x130_div30_clone_xir->Scale(scale_data_xir);
        h_elastic_x140_div30_clone_xir->Scale(scale_data_xir);
        h_elastic_x150_div30_clone_xir->Scale(scale_data_xir);
        h_singlediss_clone_xir->Scale(scale_data_xir);
		  
	     h_elastic_x130_div30_clone_mass->Add(h_elastic_x140_div30_clone_mass);
	     h_elastic_x130_div30_clone_mass->Add(h_elastic_x150_div30_clone_mass);
	     h_elastic_x130_div30_clone_pt->Add(h_elastic_x140_div30_clone_pt);
	     h_elastic_x130_div30_clone_pt->Add(h_elastic_x150_div30_clone_pt);
	     h_elastic_x130_div30_clone_y->Add(h_elastic_x140_div30_clone_y);
	     h_elastic_x130_div30_clone_y->Add(h_elastic_x150_div30_clone_y);
	     h_elastic_x130_div30_clone_xil->Add(h_elastic_x140_div30_clone_xil);
	     h_elastic_x130_div30_clone_xil->Add(h_elastic_x150_div30_clone_xil);
	     h_elastic_x130_div30_clone_xir->Add(h_elastic_x140_div30_clone_xir);
	     h_elastic_x130_div30_clone_xir->Add(h_elastic_x150_div30_clone_xir);
	
		  THStack hs_mass("hs_mass","");
		  hs_mass.Add(h_elastic_x130_div30_clone_mass);
		  hs_mass.Add(h_singlediss_clone_mass);
		  THStack hs_pt("hs_pt","");
		  hs_pt.Add(h_elastic_x130_div30_clone_pt);
		  hs_pt.Add(h_singlediss_clone_pt);
		  THStack hs_y("hs_y","");
		  hs_y.Add(h_elastic_x130_div30_clone_y);
		  hs_y.Add(h_singlediss_clone_y);
		  THStack hs_xil("hs_xil","");
		  hs_xil.Add(h_elastic_x130_div30_clone_xil);
		  hs_xil.Add(h_singlediss_clone_xil);
		  THStack hs_xir("hs_xir","");
		  hs_xir.Add(h_elastic_x130_div30_clone_xir);
		  hs_xir.Add(h_singlediss_clone_xir);

		  TFile* f4 = new TFile("2017/comp_data_MC_2017.root", "RECREATE");
		  hs_mass.Draw("HIST");
		  hs_pt.Draw("HIST");
		  hs_y.Draw("HIST");
		  hs_xil.Draw("HIST");
		  hs_xir.Draw("HIST");
		  //h_elastic_x130->Draw("HIST");
		  //h_elastic_x130_div30->Draw("HIST SAME");
		  //h_elastic_x140_div30->Draw("HIST SAME");
		  //h_elastic_x150_div30->Draw("HIST SAME");
		  //h_singlediss->Draw("HIST");
		  h_data_mass->Draw("SAME");
		  h_data_pt->Draw("SAME");
		  h_data_y->Draw("SAME");
		  h_data_xil->Draw("SAME");
		  h_data_xir->Draw("SAME");
		  hs_mass.Write();
		  hs_pt.Write();
		  hs_y.Write();
		  hs_xil.Write();
		  hs_xir.Write();
	     h_data_mass->Write();
	     h_data_pt->Write();
	     h_data_y->Write();
	     h_data_xil->Write();
	     h_data_xir->Write();
		  f4->Close();
		 // }
		  return 0;
}
