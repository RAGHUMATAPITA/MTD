#include <iostream>
#include "TROOT.h"
#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TVector.h"
#include <cstdlib>
#include <string>
#include "TString.h"
#include <chrono>
#include <fstream>

void setPalette()
{
   gStyle->SetPalette(55);
   const Int_t NRGBs = 5;
   const Int_t NCont = 255;
   Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
   Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
   Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
   Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
   TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
   gStyle->SetNumberContours(NCont);
}

void dedx_copy()
{
  setPalette();
  //gStyle->SetPalette(1,0);
  gStyle->SetTitleOffset(1.35, "y");
  gStyle->SetTitleOffset(1.2, "x");
  //TFile *file = new TFile("HiMTDTree_pp_fulldata_MB.root");
  //TFile *file = new TFile("HiMTDTree_ppRERECO_fulldata_MB.root");
  //TFile *file = new TFile("HiMTDTree_PbPbRERECO_fulldata_MB.root");
  TFile *file = new TFile("HiMTDTree_PbPbRERECO_B2p0T_fulldata_MB.root");
  TDirectory *dir = (TDirectory*)file->Get("timeAna");
  
  TTree *MTDtree =(TTree*)dir->Get("HiMTDTree");

  vector<int> *vertex_ntrk_eta3=0;
  
  vector<float> *vertex_t0=0;
  vector<float> *vertex_t0Err=0;

  vector<int> *reco_track_vtxidx=0;
  vector<float> *reco_track_pathlength=0;
  vector<float> *reco_track_p=0;
  vector<float> *reco_track_tMTD=0;
  vector<float> *reco_track_tMTDerr=0;
  vector<float> *reco_track_eta=0;
  vector<float> *reco_track_pt=0;
  vector<float> *reco_track_phi=0;

  vector<float> *gen_track_eta=0;
  vector<float> *gen_track_pt=0;
  vector<float> *gen_track_phi=0;
  vector<float> *gen_track_mass=0;
  vector<int> *gen_track_pdgid=0;

  TBranch *b_vertex_ntrk_eta3;
  TBranch *b_vertex_t0;
  TBranch *b_vertex_t0Err;

  TBranch *b_reco_track_vtxidx;
  TBranch *b_reco_track_pathlength;
  TBranch *b_reco_track_p;
  TBranch *b_reco_track_tMTD;
  TBranch *b_reco_track_tMTDerr;
  TBranch *b_reco_track_eta;
  TBranch *b_reco_track_pt;
  TBranch *b_reco_track_phi;
  
  TBranch *b_gen_track_eta;
  TBranch *b_gen_track_pt;
  TBranch *b_gen_track_mass;
  TBranch *b_gen_track_phi;
  TBranch *b_gen_track_pdgid;


  MTDtree->SetBranchAddress("Reco_Vertex_nTrk_AEta3", &vertex_ntrk_eta3, &b_vertex_ntrk_eta3);
  
  MTDtree->SetBranchAddress("Reco_Vertex_t0", &vertex_t0, &b_vertex_t0);
  MTDtree->SetBranchAddress("Reco_Vertex_t0Err", &vertex_t0Err, &b_vertex_t0Err);

  MTDtree->SetBranchAddress("Reco_Track_vtxIdx", &reco_track_vtxidx, &b_reco_track_vtxidx);
  MTDtree->SetBranchAddress("Reco_Track_pathLength", &reco_track_pathlength, &b_reco_track_pathlength);
  MTDtree->SetBranchAddress("Reco_Track_p", &reco_track_p, &b_reco_track_p);
  MTDtree->SetBranchAddress("Reco_Track_tMTD", &reco_track_tMTD, &b_reco_track_tMTD);
  MTDtree->SetBranchAddress("Reco_Track_tMTDErr", &reco_track_tMTDerr, &b_reco_track_tMTDerr);
  MTDtree->SetBranchAddress("Reco_Track_eta", &reco_track_eta, &b_reco_track_eta);
  MTDtree->SetBranchAddress("Reco_Track_pt", &reco_track_pt, &b_reco_track_pt);
  MTDtree->SetBranchAddress("Reco_Track_phi", &reco_track_phi, &b_reco_track_phi);

  MTDtree->SetBranchAddress("Gen_Track_eta", &gen_track_eta, &b_gen_track_eta);
  MTDtree->SetBranchAddress("Gen_Track_pt", &gen_track_pt, &b_gen_track_pt);
  MTDtree->SetBranchAddress("Gen_Track_phi", &gen_track_phi, &b_gen_track_phi);
  MTDtree->SetBranchAddress("Gen_Track_mass", &gen_track_mass, &b_gen_track_mass);
  MTDtree->SetBranchAddress("Gen_Track_pdgId", &gen_track_pdgid, &b_gen_track_pdgid);


  TH1F* hDeltaR = new TH1F("hDeltaR", "", 1000, 0., 0.1);

  TH1F* hMTDErr = new TH1F("hMTDErr", "", 160, 0., 0.08);

  TH2F* hpT_Eta_before = new TH2F("hpT_Eta_before", "", 600, -3.0, 3.0, 500, 0., 5.);
  TH2F* hpT_Eta_after = new TH2F("hpT_Eta_after", "", 600, -3.0, 3.0, 500, 0., 5.);

  TH2F* hp_Eta_before = new TH2F("hp_Eta_before", "", 600, -3.0, 3.0, 1000, 0., 10.);
  TH2F* hp_Eta_after = new TH2F("hp_Eta_after", "", 600, -3.0, 3.0, 1000, 0., 10.);
  
  TH2F* hinvbeta_p_CAL_UMatch= new TH2F("hinvbeta_p_CAL_UMatch", "", 1000, 0, 5, 1000, 0.9, 1.7);
  TH2F* hinvbeta_p_CAL_Match= new TH2F("hinvbeta_p_CAL_Match", "", 1000, 0, 5, 1000, 0.9, 1.7);
  TH2F* hinvbeta_p_ETL= new TH2F("hinvbeta_p_ETL", "", 1000, 0, 5, 1000, 0.9, 1.7);
  TH2F* hinvbeta_p_BTL= new TH2F("hinvbeta_p_BTL", "", 1000, 0, 5, 1000, 0.9, 1.7);
  
  TH2F* hinvbeta_pion_reco = new TH2F("hinvbeta_pion_reco", "", 500, 0., 10., 1000, 0.9, 1.7);
  TH2F* hinvbeta_pion_gen = new TH2F("hinvbeta_pion_gen", "", 500, 0., 10., 1000, 0.9, 1.7);
  TH2F* hinvbeta_pion_reco_recohypo= new TH2F("hinvbeta_pion_reco_recohypo", "", 500, 0., 10., 250, -0.1, 0.1);

  TH2F* hinvbeta_pion_ETL_reco = new TH2F("hinvbeta_pion_ETL_reco", "", 500, 0., 10., 1000, 0.9, 1.7);
  TH2F* hinvbeta_pion_ETL_gen = new TH2F("hinvbeta_pion_ETL_gen", "", 500, 0., 10., 1000, 0.9, 1.7);
  TH2F* hinvbeta_pion_ETL_reco_recohypo= new TH2F("hinvbeta_pion_ETL_reco_recohypo", "", 500, 0., 10., 250, -0.1, 0.1);

  TH2F* hinvbeta_pion_BTL_reco = new TH2F("hinvbeta_pion_BTL_reco", "", 500, 0., 10., 1000, 0.9, 1.7);
  TH2F* hinvbeta_pion_BTL_gen = new TH2F("hinvbeta_pion_BTL_gen", "", 500, 0., 10., 1000, 0.9, 1.7);
  TH2F* hinvbeta_pion_BTL_reco_recohypo= new TH2F("hinvbeta_pion_BTL_reco_recohypo", "", 500, 0., 10., 250, -0.1, 0.1);

  TH2F* hinvbeta_kaon_reco = new TH2F("hinvbeta_kaon_reco", "", 500, 0., 10., 1000, 0.9, 1.7);
  TH2F* hinvbeta_kaon_gen = new TH2F("hinvbeta_kaon_gen", "", 500, 0., 10., 1000, 0.9, 1.7);
  TH2F* hinvbeta_kaon_reco_recohypo= new TH2F("hinvbeta_kaon_reco_recohypo", "", 500, 0., 10., 250, -0.1, 0.1);

  TH2F* hinvbeta_kaon_ETL_reco = new TH2F("hinvbeta_kaon_ETL_reco", "", 500, 0., 10., 1000, 0.9, 1.7);
  TH2F* hinvbeta_kaon_ETL_gen = new TH2F("hinvbeta_kaon_ETL_gen", "", 500, 0., 10., 1000, 0.9, 1.7);
  TH2F* hinvbeta_kaon_ETL_reco_recohypo= new TH2F("hinvbeta_kaon_ETL_reco_recohypo", "", 500, 0., 10., 250, -0.1, 0.1);

  TH2F* hinvbeta_kaon_BTL_reco = new TH2F("hinvbeta_kaon_BTL_reco", "", 500, 0., 10., 1000, 0.9, 1.7);
  TH2F* hinvbeta_kaon_BTL_gen = new TH2F("hinvbeta_kaon_BTL_gen", "", 500, 0., 10., 1000, 0.9, 1.7);
  TH2F* hinvbeta_kaon_BTL_reco_recohypo= new TH2F("hinvbeta_kaon_BTL_reco_recohypo", "", 500, 0., 10., 250, -0.1, 0.1);

  TH2F* hinvbeta_proton_reco = new TH2F("hinvbeta_proton_reco", "", 500, 0., 10., 1000, 0.9, 1.7);
  TH2F* hinvbeta_proton_gen = new TH2F("hinvbeta_proton_gen", "", 500, 0., 10., 1000, 0.9, 1.7);
  TH2F* hinvbeta_proton_reco_recohypo= new TH2F("hinvbeta_proton_reco_recohypo", "", 500, 0., 10., 250, -0.1, 0.1);

  TH2F* hinvbeta_proton_ETL_reco = new TH2F("hinvbeta_proton_ETL_reco", "", 500, 0., 10., 1000, 0.9, 1.7);
  TH2F* hinvbeta_proton_ETL_gen = new TH2F("hinvbeta_proton_ETL_gen", "", 500, 0., 10., 1000, 0.9, 1.7);
  TH2F* hinvbeta_proton_ETL_reco_recohypo= new TH2F("hinvbeta_proton_ETL_reco_recohypo", "", 500, 0., 10., 250, -0.1, 0.1);
  
  TH2F* hinvbeta_proton_BTL_reco = new TH2F("hinvbeta_proton_BTL_reco", "", 500, 0., 10., 1000, 0.9, 1.7);
  TH2F* hinvbeta_proton_BTL_gen = new TH2F("hinvbeta_proton_BTL_gen", "", 500, 0., 10., 1000, 0.9, 1.7);
  TH2F* hinvbeta_proton_BTL_reco_recohypo= new TH2F("hinvbeta_proton_BTL_reco_recohypo", "", 500, 0., 10., 250, -0.1, 0.1);
  
  const float c_cm_ns   = 2.99792458e1; //[cm/ns]

  Int_t nevent = MTDtree->GetEntries();

  std::cout<<"nevent is: "<<nevent<<std::endl;
  for(Int_t ievt = 0; ievt < nevent; ievt++)
    {
      Long64_t tentry = MTDtree->LoadTree(ievt);
      if(tentry < 0) break;

      b_vertex_ntrk_eta3->GetEntry(tentry);

      b_vertex_t0->GetEntry(tentry);
      b_vertex_t0Err->GetEntry(tentry);

      b_reco_track_vtxidx->GetEntry(tentry);
      b_reco_track_pathlength->GetEntry(tentry);
      b_reco_track_p->GetEntry(tentry);
      b_reco_track_tMTD->GetEntry(tentry);
      b_reco_track_tMTDerr->GetEntry(tentry);
      b_reco_track_eta->GetEntry(tentry);
      b_reco_track_pt->GetEntry(tentry);
      b_reco_track_phi->GetEntry(tentry);

      b_gen_track_eta->GetEntry(tentry);
      b_gen_track_pt->GetEntry(tentry);
      b_gen_track_phi->GetEntry(tentry);
      b_gen_track_mass->GetEntry(tentry);
      b_gen_track_pdgid->GetEntry(tentry);

      const std::vector<float> MASS = {0.13957018, 0.493677, 0.9382720813};

      std::vector<float>reco_Time_vtx;
      
      for(int j = 0; j < vertex_t0->size(); j++)
	{
	  reco_Time_vtx.push_back((*vertex_t0)[j]);
	  
	  for (unsigned int i = 0; i< reco_track_vtxidx->size(); i++)
	    {
	      bool match_track = kFALSE;
	      float vtxidx = (*reco_track_vtxidx)[i];
	      float reco_eta = (*reco_track_eta)[i];
	      float reco_totmom = (*reco_track_p)[i];
	      float reco_Time_MTD = (*reco_track_tMTD)[i];
	      float reco_Time_MTDerr = (*reco_track_tMTDerr)[i];
	      float reco_Path_length = (*reco_track_pathlength)[i];
	      float reco_Time_vertex = reco_Time_vtx[j];
	      float reco_pt = (*reco_track_pt)[i];
	      float reco_phi = (*reco_track_phi)[i];
	      
	      float reco_tdiff = (reco_Time_MTD - reco_Time_vertex);
	      float beta_cal = reco_Path_length/(c_cm_ns*reco_tdiff);

	      hpT_Eta_before->Fill(reco_eta, reco_pt);
	      hp_Eta_before->Fill(reco_eta, reco_totmom);
	      
	      if(reco_Time_MTDerr <= 0.) continue;

	      hpT_Eta_after->Fill(reco_eta, reco_pt);
	      hp_Eta_after->Fill(reco_eta, reco_totmom);
	      
	      if(TMath::Abs(reco_eta) < 1.4)
		{
		  if(reco_pt <= 0.8) continue;
		}
	      else
		{
		  if(reco_totmom <= 0.7) continue;
		}
	      
	      if(reco_Time_MTDerr > 0.)
		{
		  hMTDErr->Fill(reco_Time_MTDerr);

		  hinvbeta_p_CAL_UMatch->Fill(reco_totmom, 1./beta_cal);
		  
		  if(TMath::Abs(reco_eta) < 1.5)
		    {
		      hinvbeta_p_BTL->Fill(reco_totmom, 1./beta_cal);
		    }
		  if(TMath::Abs(reco_eta) > 1.6)
		    {
		      hinvbeta_p_ETL->Fill(reco_totmom, 1./beta_cal);
		    }
		}
	      	      
	      float match_mass = -99.;
	      int match_pdgid = -99;
	      
	      for(unsigned int j = 0; j< gen_track_phi->size(); j++)
		{
		  float gen_eta = (*gen_track_eta)[j];
		  float gen_phi = (*gen_track_phi)[j];
		  float gen_mass = (*gen_track_mass)[j];
		  int gen_pdgid = (*gen_track_pdgid)[j];
		  
		  float dEta = (reco_eta - gen_eta);
		  float dPhi = (reco_phi - gen_phi);
		  
		  float Delta_R = TMath::Sqrt(pow(dEta,2) + pow(dPhi,2));
		  hDeltaR->Fill(TMath::Abs(Delta_R));
		  
		  if(TMath::Abs(Delta_R) < 0.01)
		    {
		      match_track = kTRUE;
		      match_mass = gen_mass;
		      match_pdgid = gen_pdgid; 
		      break;
		    }
		  else match_track = kFALSE;
		  
		}//~~~~~~gen track loop end
	      
	      if(match_track)
		{
		  if(reco_Time_MTDerr <= 0.) continue;
		  
		  if(TMath::Abs(reco_eta) < 1.4)
		    {
		      if(reco_pt <= 0.8) continue;
		    }
		  else
		    {
		      if(reco_totmom <= 0.7) continue;
		    }
		  
		  hinvbeta_p_CAL_Match->Fill(reco_totmom, 1./beta_cal);

		  /*
		  if(TMath::Abs(reco_eta) < 1.5)
		    {
                      hinvbeta_p_BTL->Fill(reco_totmom, 1./beta_cal);
                    }
                  if(TMath::Abs(reco_eta) > 1.6)
	            {
                      hinvbeta_p_ETL->Fill(reco_totmom, 1./beta_cal);
                    }
		  */
		  
		  if(match_pdgid == -99. || match_mass == -99.) continue;
		  
		  if(TMath::Abs(match_pdgid) == 211) //~~~~~~for pion
		    {
		      float Invbeta_pion_reco = 1./beta_cal;
		      
		      //float sterm_pion = (0.1395701*pow(1.,2))/reco_totmom;
		      //float Invbeta_pion_gen = TMath::Sqrt(1.+ pow(sterm_pion,2));
		      
		      float Invbeta_pion_gen = TMath::Sqrt(pow(0.1395701,2) + pow(reco_totmom,2))/reco_totmom;
		      float Invbeta_diff_pion = Invbeta_pion_reco - Invbeta_pion_gen;
		      
		      hinvbeta_pion_reco->Fill(reco_totmom, Invbeta_pion_reco);
		      hinvbeta_pion_gen->Fill(reco_totmom, Invbeta_pion_gen);
		      hinvbeta_pion_reco_recohypo->Fill(reco_totmom, Invbeta_diff_pion);

		      if(TMath::Abs(reco_eta) > 1.6)
			{
			  hinvbeta_pion_ETL_reco->Fill(reco_totmom, Invbeta_pion_reco);
			  hinvbeta_pion_ETL_gen->Fill(reco_totmom, Invbeta_pion_gen);
			  hinvbeta_pion_ETL_reco_recohypo->Fill(reco_totmom, Invbeta_diff_pion);
			}
		      if(TMath::Abs(reco_eta) < 1.5)
			{
			  hinvbeta_pion_BTL_reco->Fill(reco_totmom, Invbeta_pion_reco);
                          hinvbeta_pion_BTL_gen->Fill(reco_totmom, Invbeta_pion_gen);
                          hinvbeta_pion_BTL_reco_recohypo->Fill(reco_totmom, Invbeta_diff_pion);
			}
		    }
		  
		  else if(TMath::Abs(match_pdgid) == 321) //~~~~~~for kaon
		    {
		      float Invbeta_kaon_reco = 1./beta_cal;
		      float Invbeta_kaon_gen = TMath::Sqrt(pow(0.493677,2) + pow(reco_totmom,2))/reco_totmom;
		      float Invbeta_diff_kaon = Invbeta_kaon_reco - Invbeta_kaon_gen;
		      
		      hinvbeta_kaon_reco->Fill(reco_totmom, Invbeta_kaon_reco);
		      hinvbeta_kaon_gen->Fill(reco_totmom, Invbeta_kaon_gen);
		      hinvbeta_kaon_reco_recohypo->Fill(reco_totmom, Invbeta_diff_kaon);

		      if(TMath::Abs(reco_eta) > 1.6)
                        {
                          hinvbeta_kaon_ETL_reco->Fill(reco_totmom, Invbeta_kaon_reco);
                          hinvbeta_kaon_ETL_gen->Fill(reco_totmom, Invbeta_kaon_gen);
                          hinvbeta_kaon_ETL_reco_recohypo->Fill(reco_totmom, Invbeta_diff_kaon);
                        }
                      if(TMath::Abs(reco_eta) < 1.5)
                        {
                          hinvbeta_kaon_BTL_reco->Fill(reco_totmom, Invbeta_kaon_reco);
                          hinvbeta_kaon_BTL_gen->Fill(reco_totmom, Invbeta_kaon_gen);
                          hinvbeta_kaon_BTL_reco_recohypo->Fill(reco_totmom, Invbeta_diff_kaon);
                        }

		    }
		  else if(TMath::Abs(match_pdgid) == 2212) //~~~~~~for proton
		    {
		      float Invbeta_proton_reco = 1./beta_cal;
		      float Invbeta_proton_gen = TMath::Sqrt(pow(0.9382720813,2) + pow(reco_totmom,2))/reco_totmom;
		      float Invbeta_diff_proton = Invbeta_proton_reco - Invbeta_proton_gen;
		      
		      hinvbeta_proton_reco->Fill(reco_totmom, Invbeta_proton_reco);
		      hinvbeta_proton_gen->Fill(reco_totmom, Invbeta_proton_gen);
		      hinvbeta_proton_reco_recohypo->Fill(reco_totmom, Invbeta_diff_proton);

		      if(TMath::Abs(reco_eta) > 1.6)
                        {
                          hinvbeta_proton_ETL_reco->Fill(reco_totmom, Invbeta_proton_reco);
                          hinvbeta_proton_ETL_gen->Fill(reco_totmom, Invbeta_proton_gen);
                          hinvbeta_proton_ETL_reco_recohypo->Fill(reco_totmom, Invbeta_diff_proton);
                        }
                      if(TMath::Abs(reco_eta) < 1.5)
                        {
                          hinvbeta_proton_BTL_reco->Fill(reco_totmom, Invbeta_proton_reco);
                          hinvbeta_proton_BTL_gen->Fill(reco_totmom, Invbeta_proton_gen);
                          hinvbeta_proton_BTL_reco_recohypo->Fill(reco_totmom, Invbeta_diff_proton);
                        }
		    }
		  else continue;
		} //~~~~~ if match_track condition end
	    }//~~~~~reco track loop end
	}
    } //~~~~~~event loop end

  TFile* fout = new TFile("MTD_out_B2p0T_newEta.root", "recreate");
  hpT_Eta_before->Write();
  hpT_Eta_after->Write();
  hp_Eta_before->Write();
  hp_Eta_after->Write();
  fout->Write();
  fout->Close();
  /*
  TCanvas* c_DeltaR = new TCanvas("c_DeltaR", "", 400, 380);
  c_DeltaR->cd();
  hDeltaR->Draw("EP");
  c_DeltaR->SaveAs("pdf_files/PbPbRERECO_B2p0T_DeltaR.png");
    
  TCanvas* c_MTDErr = new TCanvas("c_MTDErr", "", 400, 380);
  c_MTDErr->cd();
  hMTDErr->SetMarkerStyle(kOpenCircle);
  hMTDErr->SetMarkerSize(0.5);
  hMTDErr->GetYaxis()->SetTitle("counts");
  hMTDErr->GetYaxis()->CenterTitle(1);
  hMTDErr->GetXaxis()->SetTitle("tMTDError");
  hMTDErr->GetXaxis()->CenterTitle(1);
  hMTDErr->Draw("COLZ");
  c_MTDErr->SaveAs("pdf_files/PbPbRERECO_B2p0T_MTDErr.png");
 
  TCanvas* c_CAL_UMatch = new TCanvas("c_CAL_UMatch", "", 400, 380);
  c_CAL_UMatch->cd();
  c_CAL_UMatch->SetLogz();
  hinvbeta_p_CAL_UMatch->Draw("COLZ");
  hinvbeta_p_CAL_UMatch->GetYaxis()->SetTitle("1/#beta");
  hinvbeta_p_CAL_UMatch->GetYaxis()->CenterTitle(1);
  hinvbeta_p_CAL_UMatch->GetXaxis()->SetTitle("p [GeV]");
  hinvbeta_p_CAL_UMatch->GetXaxis()->CenterTitle(1);
  c_CAL_UMatch->SaveAs("pdf_files/PbPbRERECO_B2p0T_IntBeta_p_UMatch.png");
  
  TCanvas* c_CAL_Match = new TCanvas("c_CAL_Match", "", 400, 380);
  c_CAL_Match->cd();
  c_CAL_Match->SetLogz();
  hinvbeta_p_CAL_Match->Draw("COLZ");
  hinvbeta_p_CAL_Match->GetYaxis()->SetTitle("1/#beta");
  hinvbeta_p_CAL_Match->GetYaxis()->CenterTitle(1);
  hinvbeta_p_CAL_Match->GetXaxis()->SetTitle("p [GeV]");
  hinvbeta_p_CAL_Match->GetXaxis()->CenterTitle(1);
  c_CAL_Match->SaveAs("pdf_files/PbPbRERECO_B2p0T_IntBeta_p_Match.png");

  TCanvas* c_ETL = new TCanvas("c_ETL", "", 400, 380);
  c_ETL->cd();
  c_ETL->SetLogz();
  hinvbeta_p_ETL->Draw("COLZ");
  hinvbeta_p_ETL->GetYaxis()->SetTitle("1/#beta");
  hinvbeta_p_ETL->GetYaxis()->CenterTitle(1);
  hinvbeta_p_ETL->GetXaxis()->SetTitle("p [GeV]");
  hinvbeta_p_ETL->GetXaxis()->CenterTitle(1);
  c_ETL->SaveAs("pdf_files/PbPbRERECO_B2p0T_IntBeta_p_ETL.png");
  
  TCanvas* c_BTL = new TCanvas("c_BTL", "", 400, 380);
  c_BTL->cd();
  c_BTL->SetLogz();
  hinvbeta_p_BTL->Draw("COLZ");
  hinvbeta_p_BTL->GetYaxis()->SetTitle("1/#beta");
  hinvbeta_p_BTL->GetYaxis()->CenterTitle(1);
  hinvbeta_p_BTL->GetXaxis()->SetTitle("p [GeV]");
  hinvbeta_p_BTL->GetXaxis()->CenterTitle(1);
  c_BTL->SaveAs("pdf_files/PbPbRERECO_B2p0T_IntBeta_p_BTL.png");
  
  TCanvas* c_pion = new TCanvas("c_pion", "", 600, 300);
  c_pion->Divide(2,1);
  c_pion->cd(1);
  c_pion->SetLogz();
  hinvbeta_pion_reco->Draw("COLZ");
  hinvbeta_pion_reco->GetYaxis()->SetTitle("1/#beta");
  hinvbeta_pion_reco->GetYaxis()->CenterTitle(1);
  hinvbeta_pion_reco->GetXaxis()->SetTitle("p [GeV]");
  hinvbeta_pion_reco->GetXaxis()->CenterTitle(1);

  
  c_pion->cd(2);
  c_pion->SetLogz();
  hinvbeta_pion_gen->Draw("COLZ");
  hinvbeta_pion_gen->GetYaxis()->SetTitle("1/#beta_{#pi}");
  hinvbeta_pion_gen->GetYaxis()->CenterTitle(1);
  hinvbeta_pion_gen->GetXaxis()->SetTitle("p [GeV]");
  hinvbeta_pion_gen->GetXaxis()->CenterTitle(1);

  c_pion->SaveAs("pdf_files/PbPbRERECO_B2p0T_IntBeta_p_Gen_Reco_pion.png");
  
  TCanvas* c_pion_reco_hypo = new TCanvas("c_pion_reco_hypo", "", 400, 380);
  c_pion_reco_hypo->SetLogz();
  c_pion_reco_hypo->cd();
  hinvbeta_pion_reco_recohypo->Draw("COLZ");
  hinvbeta_pion_reco_recohypo->GetYaxis()->SetTitle("1/#beta - 1/#beta_{#pi}");
  hinvbeta_pion_reco_recohypo->GetYaxis()->CenterTitle(1);
  hinvbeta_pion_reco_recohypo->GetXaxis()->SetTitle("p [GeV]");
  hinvbeta_pion_reco_recohypo->GetXaxis()->CenterTitle(1);
  c_pion_reco_hypo->SaveAs("pdf_files/PbPbRERECO_B2p0T_IntBeta_p_Recohypo_pion.png");
  
  TCanvas* c_pion_ETL = new TCanvas("c_pion_ETL", "", 600, 300);
  c_pion_ETL->Divide(2,1);
  c_pion_ETL->cd(1);
  c_pion_ETL->SetLogz();
  hinvbeta_pion_ETL_reco->Draw("COLZ");
  hinvbeta_pion_ETL_reco->GetYaxis()->SetTitle("1/#beta");
  hinvbeta_pion_ETL_reco->GetYaxis()->CenterTitle(1);
  hinvbeta_pion_ETL_reco->GetXaxis()->SetTitle("p [GeV]");
  hinvbeta_pion_ETL_reco->GetXaxis()->CenterTitle(1);
  
  c_pion_ETL->cd(2);
  c_pion_ETL->SetLogz();
  hinvbeta_pion_ETL_gen->Draw("COLZ");
  hinvbeta_pion_ETL_gen->GetYaxis()->SetTitle("1/#beta_{#pi}");
  hinvbeta_pion_ETL_gen->GetYaxis()->CenterTitle(1);
  hinvbeta_pion_ETL_gen->GetXaxis()->SetTitle("p [GeV]");
  hinvbeta_pion_ETL_gen->GetXaxis()->CenterTitle(1);
  c_pion_ETL->SaveAs("pdf_files/PbPbRERECO_B2p0T_IntBeta_p_Gen_Reco_pion_ETL.png");
  
  TCanvas* c_pion_ETL_reco_hypo = new TCanvas("c_pion_ETL_reco_hypo", "", 400, 380);
  c_pion_ETL_reco_hypo->SetLogz();
  c_pion_ETL_reco_hypo->cd();
  hinvbeta_pion_ETL_reco_recohypo->Draw("COLZ");
  hinvbeta_pion_ETL_reco_recohypo->GetYaxis()->SetTitle("1/#beta - 1/#beta_{#pi}");
  hinvbeta_pion_ETL_reco_recohypo->GetYaxis()->CenterTitle(1);
  hinvbeta_pion_ETL_reco_recohypo->GetXaxis()->SetTitle("p [GeV]");
  hinvbeta_pion_ETL_reco_recohypo->GetXaxis()->CenterTitle(1);
  c_pion_ETL_reco_hypo->SaveAs("pdf_files/PbPbRERECO_B2p0T_IntBeta_p_Recohypo_pion_ETL.png");

  TCanvas* c_pion_BTL = new TCanvas("c_pion_BTL", "", 600, 300);
  c_pion_BTL->Divide(2,1);
  c_pion_BTL->cd(1);
  c_pion_BTL->SetLogz();
  hinvbeta_pion_BTL_reco->Draw("COLZ");
  hinvbeta_pion_BTL_reco->GetYaxis()->SetTitle("1/#beta");
  hinvbeta_pion_BTL_reco->GetYaxis()->CenterTitle(1);
  hinvbeta_pion_BTL_reco->GetXaxis()->SetTitle("p [GeV]");
  hinvbeta_pion_BTL_reco->GetXaxis()->CenterTitle(1);
  
  c_pion_BTL->cd(2);
  c_pion_BTL->SetLogz();
  hinvbeta_pion_BTL_gen->Draw("COLZ");
  hinvbeta_pion_BTL_gen->GetYaxis()->SetTitle("1/#beta_{#pi}");
  hinvbeta_pion_BTL_gen->GetYaxis()->CenterTitle(1);
  hinvbeta_pion_BTL_gen->GetXaxis()->SetTitle("p [GeV]");
  hinvbeta_pion_BTL_gen->GetXaxis()->CenterTitle(1);
  c_pion_BTL->SaveAs("pdf_files/PbPbRERECO_B2p0T_IntBeta_p_Gen_Reco_pion_BTL.png");
    
  TCanvas* c_pion_BTL_reco_hypo = new TCanvas("c_pion_BTL_reco_hypo", "", 400, 380);
  c_pion_BTL_reco_hypo->SetLogz();
  c_pion_BTL_reco_hypo->cd();
  hinvbeta_pion_BTL_reco_recohypo->Draw("COLZ");
  hinvbeta_pion_BTL_reco_recohypo->GetYaxis()->SetTitle("1/#beta - 1/#beta_{#pi}");
  hinvbeta_pion_BTL_reco_recohypo->GetYaxis()->CenterTitle(1);
  hinvbeta_pion_BTL_reco_recohypo->GetXaxis()->SetTitle("p [GeV]");
  hinvbeta_pion_BTL_reco_recohypo->GetXaxis()->CenterTitle(1);
  c_pion_BTL_reco_hypo->SaveAs("pdf_files/PbPbRERECO_B2p0T_IntBeta_p_Recohypo_pion_BTL.png");

  
  //TF1* f_pion = new TF1("f_pion","[0]+([1]*pow((x-2),2))", 0.5, 5.0);
  //f_pion->SetLineColor(kBlack);
  //f_pion->SetLineWidth(3);
  
  //f_pion->SetParameter(0, 0.);
  //f_pion->SetParameter(0, 0.5);

  //hinvbeta_pion_reco_recohypo->Fit(f_pion, "R");

  //return;
  
  
  TCanvas* c_kaon = new TCanvas("c_kaon", "", 600, 300);
  c_kaon->Divide(2,1);
  c_kaon->cd(1);
  c_kaon->SetLogz();
  hinvbeta_kaon_reco->Draw("COLZ");
  hinvbeta_kaon_reco->GetYaxis()->SetTitle("1/#beta");
  hinvbeta_kaon_reco->GetYaxis()->CenterTitle(1);
  hinvbeta_kaon_reco->GetXaxis()->SetTitle("p [GeV]");
  hinvbeta_kaon_reco->GetXaxis()->CenterTitle(1);

  c_kaon->cd(2);
  c_kaon->SetLogz();
  hinvbeta_kaon_gen->Draw("COLZ");
  hinvbeta_kaon_gen->GetYaxis()->SetTitle("1/#beta_{#pi}");
  hinvbeta_kaon_gen->GetYaxis()->CenterTitle(1);
  hinvbeta_kaon_gen->GetXaxis()->SetTitle("p [GeV]");
  hinvbeta_kaon_gen->GetXaxis()->CenterTitle(1);
  c_kaon->SaveAs("pdf_files/PbPbRERECO_B2p0T_IntBeta_p_Gen_Reco_kaon.png");
  
  TCanvas* c_kaon_reco_hypo = new TCanvas("c_kaon_reco_hypo", "", 400, 380);
  c_kaon_reco_hypo->SetLogz();
  c_kaon_reco_hypo->cd();
  hinvbeta_kaon_reco_recohypo->Draw("COLZ");
  hinvbeta_kaon_reco_recohypo->GetYaxis()->SetTitle("1/#beta - 1/#beta_{K}");
  hinvbeta_kaon_reco_recohypo->GetYaxis()->CenterTitle(1);
  hinvbeta_kaon_reco_recohypo->GetXaxis()->SetTitle("p [GeV]");
  hinvbeta_kaon_reco_recohypo->GetXaxis()->CenterTitle(1);
  c_kaon_reco_hypo->SaveAs("pdf_files/PbPbRERECO_B2p0T_IntBeta_p_Recohypo_kaon.png");
 
  TCanvas* c_kaon_ETL = new TCanvas("c_kaon_ETL", "", 600, 300);
  c_kaon_ETL->Divide(2,1);
  c_kaon_ETL->cd(1);
  c_kaon_ETL->SetLogz();
  hinvbeta_kaon_ETL_reco->Draw("COLZ");
  hinvbeta_kaon_ETL_reco->GetYaxis()->SetTitle("1/#beta");
  hinvbeta_kaon_ETL_reco->GetYaxis()->CenterTitle(1);
  hinvbeta_kaon_ETL_reco->GetXaxis()->SetTitle("p [GeV]");
  hinvbeta_kaon_ETL_reco->GetXaxis()->CenterTitle(1);

  c_kaon_ETL->cd(2);
  c_kaon_ETL->SetLogz();
  hinvbeta_kaon_ETL_gen->Draw("COLZ");
  hinvbeta_kaon_ETL_gen->GetYaxis()->SetTitle("1/#beta_{#pi}");
  hinvbeta_kaon_ETL_gen->GetYaxis()->CenterTitle(1);
  hinvbeta_kaon_ETL_gen->GetXaxis()->SetTitle("p [GeV]");
  hinvbeta_kaon_ETL_gen->GetXaxis()->CenterTitle(1);
  c_kaon_ETL->SaveAs("pdf_files/PbPbRERECO_B2p0T_IntBeta_p_Gen_Reco_kaon_ETL.png");
 
  TCanvas* c_kaon_ETL_reco_hypo = new TCanvas("c_kaon_ETL_reco_hypo", "", 400, 380);
  c_kaon_ETL_reco_hypo->SetLogz();
  c_kaon_ETL_reco_hypo->cd();
  hinvbeta_kaon_ETL_reco_recohypo->Draw("COLZ");
  hinvbeta_kaon_ETL_reco_recohypo->GetYaxis()->SetTitle("1/#beta - 1/#beta_{K}");
  hinvbeta_kaon_ETL_reco_recohypo->GetYaxis()->CenterTitle(1);
  hinvbeta_kaon_ETL_reco_recohypo->GetXaxis()->SetTitle("p [GeV]");
  hinvbeta_kaon_ETL_reco_recohypo->GetXaxis()->CenterTitle(1);
  c_kaon_ETL_reco_hypo->SaveAs("pdf_files/PbPbRERECO_B2p0T_IntBeta_p_Recohypo_kaon_ETL.png");
  
  TCanvas* c_kaon_BTL = new TCanvas("c_kaon_BTL", "", 600, 300);
  c_kaon_BTL->Divide(2,1);
  c_kaon_BTL->cd(1);
  c_kaon_BTL->SetLogz();
  hinvbeta_kaon_BTL_reco->Draw("COLZ");
  hinvbeta_kaon_BTL_reco->GetYaxis()->SetTitle("1/#beta");
  hinvbeta_kaon_BTL_reco->GetYaxis()->CenterTitle(1);
  hinvbeta_kaon_BTL_reco->GetXaxis()->SetTitle("p [GeV]");
  hinvbeta_kaon_BTL_reco->GetXaxis()->CenterTitle(1);

  c_kaon_BTL->cd(2);
  c_kaon_BTL->SetLogz();
  hinvbeta_kaon_BTL_gen->Draw("COLZ");
  hinvbeta_kaon_BTL_gen->GetYaxis()->SetTitle("1/#beta_{#pi}");
  hinvbeta_kaon_BTL_gen->GetYaxis()->CenterTitle(1);
  hinvbeta_kaon_BTL_gen->GetXaxis()->SetTitle("p [GeV]");
  hinvbeta_kaon_BTL_gen->GetXaxis()->CenterTitle(1);
  c_kaon_BTL->SaveAs("pdf_files/PbPbRERECO_B2p0T_IntBeta_p_Gen_Reco_kaon_BTL.png");
  
  TCanvas* c_kaon_BTL_reco_hypo = new TCanvas("c_kaon_BTL_reco_hypo", "", 400, 380);
  c_kaon_BTL_reco_hypo->SetLogz();
  c_kaon_BTL_reco_hypo->cd();
  hinvbeta_kaon_BTL_reco_recohypo->Draw("COLZ");
  hinvbeta_kaon_BTL_reco_recohypo->GetYaxis()->SetTitle("1/#beta - 1/#beta_{K}");
  hinvbeta_kaon_BTL_reco_recohypo->GetYaxis()->CenterTitle(1);
  hinvbeta_kaon_BTL_reco_recohypo->GetXaxis()->SetTitle("p [GeV]");
  hinvbeta_kaon_BTL_reco_recohypo->GetXaxis()->CenterTitle(1);
  c_kaon_BTL_reco_hypo->SaveAs("pdf_files/PbPbRERECO_B2p0T_IntBeta_p_Recohypo_kaon_BTL.png");


  TCanvas* c_proton = new TCanvas("c_proton", "", 600, 300);
  c_proton->Divide(2,1);
  c_proton->cd(1);
  c_proton->SetLogz();
  hinvbeta_proton_reco->Draw("COLZ");
  hinvbeta_proton_reco->GetYaxis()->SetTitle("1/#beta");
  hinvbeta_proton_reco->GetYaxis()->CenterTitle(1);
  hinvbeta_proton_reco->GetXaxis()->SetTitle("p [GeV]");
  hinvbeta_proton_reco->GetXaxis()->CenterTitle(1);
  
  c_proton->cd(2);
  c_proton->SetLogz();
  hinvbeta_proton_gen->Draw("COLZ");
  hinvbeta_proton_gen->GetYaxis()->SetTitle("1/#beta_{#pi}");
  hinvbeta_proton_gen->GetYaxis()->CenterTitle(1);
  hinvbeta_proton_gen->GetXaxis()->SetTitle("p [GeV]");
  hinvbeta_proton_gen->GetXaxis()->CenterTitle(1);
  c_proton->SaveAs("pdf_files/PbPbRERECO_B2p0T_IntBeta_p_Gen_Reco_proton.png");
  
  TCanvas* c_proton_reco_hypo = new TCanvas("c_proton_reco_hypo", "", 400, 380);
  c_proton_reco_hypo->SetLogz();
  c_proton_reco_hypo->cd();
  hinvbeta_proton_reco_recohypo->Draw("COLZ");
  hinvbeta_proton_reco_recohypo->GetYaxis()->SetTitle("1/#beta - 1/#beta_{P}");
  hinvbeta_proton_reco_recohypo->GetYaxis()->CenterTitle(1);
  hinvbeta_proton_reco_recohypo->GetXaxis()->SetTitle("p [GeV]");
  hinvbeta_proton_reco_recohypo->GetXaxis()->CenterTitle(1);
  c_proton_reco_hypo->SaveAs("pdf_files/PbPbRERECO_B2p0T_IntBeta_p_Recohypo_proton.png");
  
  TCanvas* c_proton_ETL = new TCanvas("c_proton_ETL", "", 600, 300);
  c_proton_ETL->Divide(2,1);
  c_proton_ETL->cd(1);
  c_proton_ETL->SetLogz();
  hinvbeta_proton_ETL_reco->Draw("COLZ");
  hinvbeta_proton_ETL_reco->GetYaxis()->SetTitle("1/#beta");
  hinvbeta_proton_ETL_reco->GetYaxis()->CenterTitle(1);
  hinvbeta_proton_ETL_reco->GetXaxis()->SetTitle("p [GeV]");
  hinvbeta_proton_ETL_reco->GetXaxis()->CenterTitle(1);
  
  c_proton_ETL->cd(2);
  c_proton_ETL->SetLogz();
  hinvbeta_proton_ETL_gen->Draw("COLZ");
  hinvbeta_proton_ETL_gen->GetYaxis()->SetTitle("1/#beta_{#pi}");
  hinvbeta_proton_ETL_gen->GetYaxis()->CenterTitle(1);
  hinvbeta_proton_ETL_gen->GetXaxis()->SetTitle("p [GeV]");
  hinvbeta_proton_ETL_gen->GetXaxis()->CenterTitle(1);
  c_proton_ETL->SaveAs("pdf_files/PbPbRERECO_B2p0T_IntBeta_p_Gen_Reco_proton_ETL.png");
  
  TCanvas* c_proton_ETL_reco_hypo = new TCanvas("c_proton_ETL_reco_hypo", "", 400, 380);
  c_proton_ETL_reco_hypo->SetLogz();
  c_proton_ETL_reco_hypo->cd();
  hinvbeta_proton_ETL_reco_recohypo->Draw("COLZ");
  hinvbeta_proton_ETL_reco_recohypo->GetYaxis()->SetTitle("1/#beta - 1/#beta_{P}");
  hinvbeta_proton_ETL_reco_recohypo->GetYaxis()->CenterTitle(1);
  hinvbeta_proton_ETL_reco_recohypo->GetXaxis()->SetTitle("p [GeV]");
  hinvbeta_proton_ETL_reco_recohypo->GetXaxis()->CenterTitle(1);
  c_proton_ETL_reco_hypo->SaveAs("pdf_files/PbPbRERECO_B2p0T_IntBeta_p_Recohypo_proton_ETL.png");
  
  TCanvas* c_proton_BTL = new TCanvas("c_proton_BTL", "", 600, 300);
  c_proton_BTL->Divide(2,1);
  c_proton_BTL->cd(1);
  c_proton_BTL->SetLogz();
  hinvbeta_proton_BTL_reco->Draw("COLZ");
  hinvbeta_proton_BTL_reco->GetYaxis()->SetTitle("1/#beta");
  hinvbeta_proton_BTL_reco->GetYaxis()->CenterTitle(1);
  hinvbeta_proton_BTL_reco->GetXaxis()->SetTitle("p [GeV]");
  hinvbeta_proton_BTL_reco->GetXaxis()->CenterTitle(1);
  
  c_proton_BTL->cd(2);
  c_proton_BTL->SetLogz();
  hinvbeta_proton_BTL_gen->Draw("COLZ");
  hinvbeta_proton_BTL_gen->GetYaxis()->SetTitle("1/#beta_{#pi}");
  hinvbeta_proton_BTL_gen->GetYaxis()->CenterTitle(1);
  hinvbeta_proton_BTL_gen->GetXaxis()->SetTitle("p [GeV]");
  hinvbeta_proton_BTL_gen->GetXaxis()->CenterTitle(1);
  c_proton_BTL->SaveAs("pdf_files/PbPbRERECO_B2p0T_IntBeta_p_Gen_Reco_proton_BTL.png");
  
  TCanvas* c_proton_BTL_reco_hypo = new TCanvas("c_proton_BTL_reco_hypo", "", 400, 380);
  c_proton_BTL_reco_hypo->SetLogz();
  c_proton_BTL_reco_hypo->cd();
  hinvbeta_proton_BTL_reco_recohypo->Draw("COLZ");
  hinvbeta_proton_BTL_reco_recohypo->GetYaxis()->SetTitle("1/#beta - 1/#beta_{P}");
  hinvbeta_proton_BTL_reco_recohypo->GetYaxis()->CenterTitle(1);
  hinvbeta_proton_BTL_reco_recohypo->GetXaxis()->SetTitle("p [GeV]");
  hinvbeta_proton_BTL_reco_recohypo->GetXaxis()->CenterTitle(1);
  c_proton_BTL_reco_hypo->SaveAs("pdf_files/PbPbRERECO_B2p0T_IntBeta_p_Recohypo_proton_BTL.png");

  */
}

