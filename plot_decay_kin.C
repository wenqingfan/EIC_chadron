#include "bins.h"
#include "diff_det.h"
#include "Pad2x1.h"
#include "Pad3x1.h"
#include "TStyle.h"
#include "TGraphErrors.h"
using namespace std;
unsigned int verbosity = 0;

const float R_DIRC = 0.76; // used to be 90cm in ATHENA but now at 76cm in ePIC

TH2D* h2d_K_D0_p_vs_eta[smearbin][bfieldbin][pidbin][etabin][pptbin] = {0};
TH2D* h2d_pi_D0_p_vs_eta[smearbin][bfieldbin][pidbin][etabin][pptbin] = {0};

TH2D* h2d_K_Lc_p_vs_eta[smearbin][bfieldbin][pidbin][etabin][pptbin] = {0};
TH2D* h2d_pi_Lc_p_vs_eta[smearbin][bfieldbin][pidbin][etabin][pptbin] = {0};
TH2D* h2d_p_Lc_p_vs_eta[smearbin][bfieldbin][pidbin][etabin][pptbin] = {0};

const float EPIC_FIELD_14 = 1.4;
const float EPIC_FIELD_17 = 1.7;

TGraph* kin_tracking_ePIC_14 = NULL;
TGraph* kin_tracking_ePIC_17 = NULL;

TGraph* kin_DIRC_ePIC_14 = NULL;
TGraph* kin_DIRC_ePIC_17 = NULL;

// DIRC firing (detectable) threshold
const float PI_DIRC[2] = {0.13, 0.25};
const float K_DIRC[2] = {0.47, 1.1};
const float P_DIRC[2] = {0.88, 2.15};

// dRICH firing (detectable) threshold
const float PI_DRICH[2] = {0.69, 1};
const float K_DRICH[2] = {2.46, 3};
const float P_DRICH[2] = {4.67, 5};

TLine* line_pi_dRICH_h = NULL;
TLine* line_pi_dRICH_e = NULL;
TLine* line_pi_DIRC = NULL;

TLine* line_K_dRICH_h = NULL;
TLine* line_K_dRICH_e = NULL;
TLine* line_K_DIRC = NULL;

double eta_to_theta(const double eta)
{
  return 2*atan(exp(-eta))/TMath::Pi()*180;
}

void set_tracking_limits()
{ // this looks like the tracking threshold values from the YR
  double kin_eta[52] = {0};
  double kin_pt_ePIC_14[52] = {0};
  double kin_pt_ePIC_17[52] = {0};
  double kin_p_ePIC_14[52] = {0};
  double kin_p_ePIC_17[52] = {0};

  kin_eta[0] = 3.0; kin_pt_ePIC_14[0] = 0.15; kin_pt_ePIC_17[0] = 0.3;
  kin_eta[1] = 2.5; kin_pt_ePIC_14[1] = 0.15; kin_pt_ePIC_17[1] = 0.3;
  kin_eta[2] = 2.5; kin_pt_ePIC_14[2] = 0.13; kin_pt_ePIC_17[2] = 0.22;
  kin_eta[3] = 2.0; kin_pt_ePIC_14[3] = 0.13; kin_pt_ePIC_17[3] = 0.22;
  kin_eta[3] = 2.0; kin_pt_ePIC_14[3] = 0.07; kin_pt_ePIC_17[3] = 0.16;
  kin_eta[3] = 1.5; kin_pt_ePIC_14[3] = 0.07; kin_pt_ePIC_17[3] = 0.16;
  kin_eta[4] = 1.5; kin_pt_ePIC_14[4] = 0.15; kin_pt_ePIC_17[4] = 0.3;
  kin_eta[5] = 1.0; kin_pt_ePIC_14[5] = 0.15; kin_pt_ePIC_17[5] = 0.3;
  for (int ibin = 6; ibin < 46; ++ibin)
  {
    if (ibin<26)
    {
      kin_eta[ibin] = 1-1./20*(ibin-6);
      kin_pt_ePIC_14[ibin] = 0.2; // All-Si 0.09
      kin_pt_ePIC_17[ibin] = 0.4; // All-Si 0.195
    }
    else
    {
      kin_eta[ibin] = 0+(-1.)/20*(ibin-6-20);
      kin_pt_ePIC_14[ibin] = 0.2; // All-Si 0.09
      kin_pt_ePIC_17[ibin] = 0.4; // All-Si 0.195
    }
  }
  kin_eta[51-0] = -3.0; kin_pt_ePIC_14[51-0] = 0.15; kin_pt_ePIC_17[51-0] = 0.3;
  kin_eta[51-1] = -2.5; kin_pt_ePIC_14[51-1] = 0.15; kin_pt_ePIC_17[51-1] = 0.3;
  kin_eta[51-2] = -2.5; kin_pt_ePIC_14[51-2] = 0.13; kin_pt_ePIC_17[51-2] = 0.22;
  kin_eta[51-3] = -2.0; kin_pt_ePIC_14[51-3] = 0.13; kin_pt_ePIC_17[51-3] = 0.22;
  kin_eta[51-3] = -2.0; kin_pt_ePIC_14[51-3] = 0.07; kin_pt_ePIC_17[51-3] = 0.16;
  kin_eta[51-3] = -1.5; kin_pt_ePIC_14[51-3] = 0.07; kin_pt_ePIC_17[51-3] = 0.16;
  kin_eta[51-4] = -1.5; kin_pt_ePIC_14[51-4] = 0.15; kin_pt_ePIC_17[51-4] = 0.3;
  kin_eta[51-5] = -1.0; kin_pt_ePIC_14[51-5] = 0.15; kin_pt_ePIC_17[51-5] = 0.3;

  for (int ibin = 0; ibin < 52; ++ibin)
  {
    if (ibin<52/2)
    {
      double kin_theta = eta_to_theta(kin_eta[ibin]);
      kin_p_ePIC_14[ibin] = kin_pt_ePIC_14[ibin]*1/cos((90-kin_theta)/180*TMath::Pi());
      kin_p_ePIC_17[ibin] = kin_pt_ePIC_17[ibin]*1/cos((90-kin_theta)/180*TMath::Pi());
    }
    else
    {
      double kin_theta = eta_to_theta(kin_eta[ibin]);
      kin_p_ePIC_14[ibin] = kin_pt_ePIC_14[ibin]*1/cos((kin_theta-90)/180*TMath::Pi());
      kin_p_ePIC_17[ibin] = kin_pt_ePIC_17[ibin]*1/cos((kin_theta-90)/180*TMath::Pi()); 
    }
  }

  kin_tracking_ePIC_14 = new TGraph(52,kin_p_ePIC_14,kin_eta);
  kin_tracking_ePIC_14->SetLineWidth(2);
  kin_tracking_ePIC_14->SetLineColor(kGray+2);
  kin_tracking_ePIC_14->SetMarkerColor(kGray+2);

  kin_tracking_ePIC_17 = new TGraph(52,kin_p_ePIC_17,kin_eta);
  kin_tracking_ePIC_17->SetLineWidth(2);
  kin_tracking_ePIC_17->SetLineStyle(2);
  kin_tracking_ePIC_17->SetLineColor(kGray+2);
  kin_tracking_ePIC_17->SetMarkerColor(kGray+2);
}

void set_PID_limits(const int option = 1)
{
  line_pi_dRICH_h = new TLine(PI_DRICH[option],1,PI_DRICH[option],3);
  line_pi_dRICH_h->SetLineColor(kRed);
  line_pi_dRICH_h->SetLineWidth(2);

  line_pi_dRICH_e = new TLine(PI_DRICH[option],-1,PI_DRICH[option],-3);
  line_pi_dRICH_e->SetLineColor(kRed);
  line_pi_dRICH_e->SetLineWidth(2);

  line_pi_DIRC = new TLine(PI_DIRC[option],-1,PI_DIRC[option],1);
  line_pi_DIRC->SetLineColor(kRed);
  line_pi_DIRC->SetLineWidth(2);

  line_K_dRICH_h = new TLine(K_DRICH[option],1,K_DRICH[option],3);
  line_K_dRICH_h->SetLineColor(kBlue);
  line_K_dRICH_h->SetLineWidth(2);

  line_K_dRICH_e = new TLine(K_DRICH[option],-1,K_DRICH[option],-3);
  line_K_dRICH_e->SetLineColor(kBlue);
  line_K_dRICH_e->SetLineWidth(2);

  line_K_DIRC = new TLine(K_DIRC[option],-1,K_DIRC[option],1);
  line_K_DIRC->SetLineColor(kBlue);
  line_K_DIRC->SetLineWidth(2);

  const float R_in = R_DIRC;
  float thr_kin_ePIC_14 = 0.3*EPIC_FIELD_14*R_in/2; 
  float thr_kin_ePIC_17 = 0.3*EPIC_FIELD_17*R_in/2; 

  double kin_eta[40] = {0};
  double kin_p_ePIC_14[40] = {0};
  double kin_p_ePIC_17[40] = {0};

  for (int ibin = 0; ibin < 40; ++ibin)
  {
    if (ibin<20)
    {
      kin_eta[ibin] = 1-1./20*(ibin);
      double kin_theta = eta_to_theta(kin_eta[ibin]);
      kin_p_ePIC_14[ibin] = thr_kin_ePIC_14*1/cos((90-kin_theta)/180*TMath::Pi());
      kin_p_ePIC_17[ibin] = thr_kin_ePIC_17*1/cos((90-kin_theta)/180*TMath::Pi());
    }
    else
    {
      kin_eta[ibin] = 0+(-1.)/20*(ibin-20);
      double kin_theta = eta_to_theta(kin_eta[ibin]);
      kin_p_ePIC_14[ibin] = thr_kin_ePIC_14*1/cos((kin_theta-90)/180*TMath::Pi());
      kin_p_ePIC_17[ibin] = thr_kin_ePIC_17*1/cos((kin_theta-90)/180*TMath::Pi());
    }
  }

  kin_DIRC_ePIC_14 = new TGraph(40,kin_p_ePIC_14,kin_eta);
  kin_DIRC_ePIC_14->SetLineWidth(2);
  kin_DIRC_ePIC_14->SetLineStyle(2);
  kin_DIRC_ePIC_14->SetLineColor(kViolet+1);
  kin_DIRC_ePIC_14->SetMarkerColor(kViolet+1);

  kin_DIRC_ePIC_17 = new TGraph(40,kin_p_ePIC_17,kin_eta);
  kin_DIRC_ePIC_17->SetLineWidth(2);
  kin_DIRC_ePIC_17->SetLineColor(kViolet+1);
  kin_DIRC_ePIC_17->SetMarkerColor(kViolet+1);
}

void plot_D0_decay(const int smear_option = 2, const int bfield_option = 0, const int pid_option = 2)
{
  {
    gStyle->SetPalette(52);
    gStyle->SetTitleFont(63,"X");
    gStyle->SetTitleFont(63,"Y");
    gStyle->SetLabelFont(63,"X");
    gStyle->SetLabelFont(63,"Y");
    gStyle->SetTextFont(63);

    gStyle->SetStatStyle(0);
    gStyle->SetTitleStyle(0);
    gStyle->SetCanvasBorderSize(0);
    gStyle->SetFrameBorderSize(0);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetStatBorderSize(0);
    gStyle->SetTitleBorderSize(0);

    gStyle->SetStatStyle(0);
    gStyle->SetTitleStyle(0);
    gStyle->SetCanvasBorderSize(0);
    gStyle->SetFrameBorderSize(0);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetStatBorderSize(0);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetOptStat(0);
    gStyle->SetFillStyle(4000);
    gStyle->SetStatStyle(0);
    gStyle->SetFrameFillColor(0);
    gStyle->SetLegendBorderSize(0);

    int ismear = smear_option; // DM smearing
    int ibfield = bfield_option; // 1.4T
    int ipid = pid_option;

    for (int ieta = 0; ieta < etabin; ++ieta)
    {
      for (int ipt = 0; ipt < pptbin; ++ipt)
      {
        {
            Pad2x1 *pad2x1 = new Pad2x1("c1",500,500,100,45,70,100);
            pad2x1->Draw();

            // first pad
            TPad* mypad = pad2x1->GetPad(1);
            mypad->SetLogz();
            mypad->SetLogx();

            // h2d_pi_D0_p_vs_eta[ismear][ibfield][ipid][ieta][ipt]->SetMaximum(1E3);
            h2d_pi_D0_p_vs_eta[ismear][ibfield][ipid][ieta][ipt]->SetMinimum(1E0);

            TH2F* htemp = new TH2F("htemp","",10,0.08,9.999,10,-4,6);
            htemp->Draw();
            htemp->GetXaxis()->SetTitle("p [GeV/c]");
            htemp->GetYaxis()->SetTitle("#eta");
            htemp->GetXaxis()->SetLabelFont(63);
            htemp->GetYaxis()->SetLabelFont(63);
            htemp->GetXaxis()->SetLabelSize(30);
            htemp->GetYaxis()->SetLabelSize(30);
            htemp->GetXaxis()->SetLabelOffset(0.01);
            htemp->GetYaxis()->SetLabelOffset(0.01);
            htemp->GetXaxis()->CenterTitle(1);
            htemp->GetYaxis()->CenterTitle(1);
            htemp->GetXaxis()->SetTitleSize(0);
            htemp->GetXaxis()->SetTitleOffset(1.5); 
            htemp->GetYaxis()->SetTitleSize(45);
            htemp->GetYaxis()->SetTitleOffset(1.0);

            h2d_pi_D0_p_vs_eta[ismear][ibfield][ipid][ieta][ipt]->Draw("colzsame");

            TLine* line1 = new TLine(0,-1,10,-1);
            line1->SetLineStyle(2);
            TLine* line2 = new TLine(0,1,10,1);
            line2->SetLineStyle(2);

            line1->Draw("same");
            line2->Draw("same");

            // kin_tracking_ePIC_17->Draw("same");
            // kin_tracking_ePIC_14->Draw("same");

            kin_DIRC_ePIC_17->Draw("same");
            kin_DIRC_ePIC_14->Draw("same");

            line_pi_DIRC->Draw("same");
            line_pi_dRICH_h->Draw("same");
            line_pi_dRICH_e->Draw("same");

            // second pad
            TPad* mypad = pad2x1->GetPad(2);
            mypad->SetLogz();
            mypad->SetLogx();

            // h2d_K_D0_p_vs_eta[ismear][ibfield][ipid][ieta][ipt]->SetMaximum(1E3);
            h2d_K_D0_p_vs_eta[ismear][ibfield][ipid][ieta][ipt]->SetMinimum(1E0);

            htemp->Draw();

            h2d_K_D0_p_vs_eta[ismear][ibfield][ipid][ieta][ipt]->Draw("colzsame");

            line1->Draw("same");
            line2->Draw("same");

            // kin_tracking_ePIC_17->Draw("same");
            // kin_tracking_ePIC_14->Draw("same");

            kin_DIRC_ePIC_17->Draw("same");
            kin_DIRC_ePIC_14->Draw("same");

            line_K_DIRC->Draw("same");
            line_K_dRICH_h->Draw("same");
            line_K_dRICH_e->Draw("same");

            // get the text pad
            pad2x1->GetPad(3);

            TLegend* leg1 = new TLegend(0.14,0.73,0.37,0.83);
            leg1->SetBorderSize(0);
            leg1->SetTextSize(0.04);
            leg1->SetFillStyle(0);
            leg1->SetMargin(0.2);

            leg1->AddEntry(kin_DIRC_ePIC_17,"low p limit to reach DIRC (1.7T)","l");
            leg1->AddEntry(kin_DIRC_ePIC_14,"low p limit to reach DIRC (1.4T)","l");
            leg1->Draw("same");

            TLegend* leg2 = new TLegend(0.63,0.73,0.86,0.83);
            leg2->SetBorderSize(0);
            leg2->SetTextSize(0.04);
            leg2->SetFillStyle(0);
            leg2->SetMargin(0.2);

            leg2->AddEntry(line_pi_DIRC,"#pi threshold","l");
            leg2->AddEntry(line_K_DIRC,"K threshold","l");
            leg2->Draw("same");

            pave1 = new TPaveText(0.30,0.88,0.28,0.88,"NDCNB");
            pave1->AddText("D^{0} #rightarrow #pi");
            pave1->SetTextFont(63);
            pave1->SetTextSize(25);
            pave1->Draw();

            pave2 = new TPaveText(0.30+0.43,0.88,0.28+0.43,0.88,"NDCNB");
            pave2->AddText("D^{0} #rightarrow K");
            pave2->SetTextFont(63);
            pave2->SetTextSize(25);
            pave2->Draw();

            pavelabelX = new TPaveText(0.51,0.04,0.54,0.04,"NDCNB");
            pavelabelX->AddText("p [GeV/c]");
            pavelabelX->SetTextFont(63);
            pavelabelX->SetTextSize(32);
            pavelabelX->SetFillColorAlpha(kWhite,0);
            pavelabelX->Draw();

            TPaveText* pave = new TPaveText(0.10,0.93,0.90,0.95,"NDCNB");
            pave->AddText(Form("Pythia e+p @ 10+100 GeV, Min Bias (Q^{2} > 10 GeV^{2}), D^{0} in #eta [%.1f, %.1f], p_{T} [%.1f, %.1f] GeV",eta_lo[ieta],eta_hi[ieta],ppt_lo[ipt],ppt_hi[ipt]));
            pave->SetTextAlign(21);
            pave->SetTextFont(63);
            pave->SetTextSize(25);
            pave->Draw();

            TPaveText* pave4 = new TPaveText(0.07,0.16,0.70,0.18,"NDCNB");
            pave4->AddText(Form("%s",pid_full_name[pid_option]));
            pave4->SetTextAlign(11);
            pave4->SetTextFont(63);
            pave4->SetTextSize(23);
            pave4->Draw();

            gROOT->ProcessLine( Form("c1->Print(\"figs/D0_decay_eta%d_pt%d.pdf\")",ieta,ipt) );
        }
      }
    }   
  }
}

void plot_Lc_decay(const int smear_option = 2, const int bfield_option = 0, const int pid_option = 2)
{
  {
    gStyle->SetTitleFont(63,"X");
    gStyle->SetTitleFont(63,"Y");
    gStyle->SetLabelFont(63,"X");
    gStyle->SetLabelFont(63,"Y");
    gStyle->SetTextFont(63);

    gStyle->SetStatStyle(0);
    gStyle->SetTitleStyle(0);
    gStyle->SetCanvasBorderSize(0);
    gStyle->SetFrameBorderSize(0);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetStatBorderSize(0);
    gStyle->SetTitleBorderSize(0);

    gStyle->SetStatStyle(0);
    gStyle->SetTitleStyle(0);
    gStyle->SetCanvasBorderSize(0);
    gStyle->SetFrameBorderSize(0);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetStatBorderSize(0);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetOptStat(0);
    gStyle->SetFillStyle(4000);
    gStyle->SetStatStyle(0);
    gStyle->SetFrameFillColor(0);
    gStyle->SetLegendBorderSize(0);

    int ismear = smear_option; // DM smearing
    int ibfield = bfield_option; // 1.4T
    int ipid = pid_option;

    for (int ieta = 0; ieta < etabin; ++ieta)
    {
      for (int ipt = 0; ipt < pptbin; ++ipt)
      {
        {
            Pad3x1 *pad3x1 = new Pad3x1("c1",500,500,100,100,70,120);
            pad3x1->Draw();

            // first pad
            TPad* mypad = pad3x1->GetPad(1);
            mypad->SetLogz();
            mypad->SetLogx();

            // h2d_pi_Lc_p_vs_eta[ismear][ibfield][ipid][ieta][ipt]->SetMaximum(1E3);
            h2d_pi_Lc_p_vs_eta[ismear][ibfield][ipid][ieta][ipt]->SetMinimum(1E0);

            TH2F* htemp = new TH2F("htemp","",10,0.08,9.999,10,-4,6);
            htemp->Draw();
            htemp->GetXaxis()->SetTitle("p [GeV/c]");
            htemp->GetYaxis()->SetTitle("#eta");
            htemp->GetXaxis()->SetLabelFont(63);
            htemp->GetYaxis()->SetLabelFont(63);
            htemp->GetXaxis()->SetLabelSize(30);
            htemp->GetYaxis()->SetLabelSize(30);
            htemp->GetXaxis()->SetLabelOffset(0.01);
            htemp->GetYaxis()->SetLabelOffset(0.01);
            htemp->GetXaxis()->CenterTitle(1);
            htemp->GetYaxis()->CenterTitle(1);
            htemp->GetXaxis()->SetTitleSize(0);
            htemp->GetXaxis()->SetTitleOffset(4.0); 
            htemp->GetYaxis()->SetTitleSize(45);
            htemp->GetYaxis()->SetTitleOffset(1.0);

            h2d_pi_Lc_p_vs_eta[ismear][ibfield][ipid][ieta][ipt]->Draw("colzsame");

            TLine* line1 = new TLine(0,-1,10,-1);
            line1->SetLineStyle(2);
            TLine* line2 = new TLine(0,1,10,1);
            line2->SetLineStyle(2);

            line1->Draw("same");
            line2->Draw("same");

            // kin_tracking_ePIC_17->Draw("same");
            // kin_tracking_ePIC_14->Draw("same");

            kin_DIRC_ePIC_17->Draw("same");
            kin_DIRC_ePIC_14->Draw("same");

            line_pi_DIRC->Draw("same");
            line_pi_dRICH_h->Draw("same");
            line_pi_dRICH_e->Draw("same");

            // second pad
            TPad* mypad = pad3x1->GetPad(2);
            mypad->SetLogz();
            mypad->SetLogx();

            // h2d_K_Lc_p_vs_eta[ismear][ibfield][ipid][ieta][ipt]->SetMaximum(1E3);
            h2d_K_Lc_p_vs_eta[ismear][ibfield][ipid][ieta][ipt]->SetMinimum(1E0);

            htemp->Draw();

            h2d_K_Lc_p_vs_eta[ismear][ibfield][ipid][ieta][ipt]->Draw("colzsame");

            line1->Draw("same");
            line2->Draw("same");

            kin_DIRC_ePIC_17->Draw("same");
            kin_DIRC_ePIC_14->Draw("same");

            line_K_DIRC->Draw("same");
            line_K_dRICH_h->Draw("same");
            line_K_dRICH_e->Draw("same");

            // third pad
            TPad* mypad = pad3x1->GetPad(3);
            mypad->SetLogz();
            mypad->SetLogx();

            // h2d_p_Lc_p_vs_eta[ismear][ibfield][ipid][ieta][ipt]->SetMaximum(1E3);
            h2d_p_Lc_p_vs_eta[ismear][ibfield][ipid][ieta][ipt]->SetMinimum(1E0);

            htemp->Draw();

            h2d_p_Lc_p_vs_eta[ismear][ibfield][ipid][ieta][ipt]->Draw("colzsame");

            line1->Draw("same");
            line2->Draw("same");

            // kin_tracking_ePIC_17->Draw("same");
            // kin_tracking_ePIC_14->Draw("same");

            kin_DIRC_ePIC_17->Draw("same");
            kin_DIRC_ePIC_14->Draw("same");

            line_K_DIRC->Draw("same");
            line_K_dRICH_h->Draw("same");
            line_K_dRICH_e->Draw("same");

            // get the text pad
            pad3x1->GetPad(4);

            TLegend* leg1 = new TLegend(0.08,0.72,0.31,0.82);
            leg1->SetBorderSize(0);
            leg1->SetTextSize(0.04);
            leg1->SetFillStyle(0);
            leg1->SetMargin(0.2);

            leg1->AddEntry(kin_DIRC_ePIC_17,"low p limit to reach DIRC (1.7T)","l");
            leg1->AddEntry(kin_DIRC_ePIC_14,"low p limit to reach DIRC (1.4T)","l");
            leg1->Draw("same");

            TLegend* leg2 = new TLegend(0.42,0.72,0.65,0.82);
            leg2->SetBorderSize(0);
            leg2->SetTextSize(0.04);
            leg2->SetFillStyle(0);
            leg2->SetMargin(0.2);

            leg2->AddEntry(line_pi_DIRC,"#pi threshold","l");
            leg2->AddEntry(line_K_DIRC,"K threshold","l");
            leg2->Draw("same");

            pave1 = new TPaveText(0.19,0.87,0.23,0.87,"NDCNB");
            pave1->AddText("#Lambda_{c} #rightarrow #pi");
            pave1->SetTextFont(63);
            pave1->SetTextSize(30);
            pave1->Draw();

            pave2 = new TPaveText(0.19+0.3,0.87,0.23+0.3,0.87,"NDCNB");
            pave2->AddText("#Lambda_{c} #rightarrow K");
            pave2->SetTextFont(63);
            pave2->SetTextSize(30);
            pave2->Draw();

            pave3 = new TPaveText(0.19+0.3*2,0.87,0.23+0.3*2,0.87,"NDCNB");
            pave3->AddText("#Lambda_{c} #rightarrow p");
            pave3->SetTextFont(63);
            pave3->SetTextSize(30);
            pave3->Draw();

            pavelabelX = new TPaveText(0.51,0.04,0.54,0.04,"NDCNB");
            pavelabelX->AddText("p [GeV/c]");
            pavelabelX->SetTextFont(63);
            pavelabelX->SetTextSize(32);
            pavelabelX->SetFillColorAlpha(kWhite,0);
            pavelabelX->Draw();

            TPaveText* pave = new TPaveText(0.10,0.93,0.90,0.95,"NDCNB");
            pave->AddText(Form("Pythia e+p @ 10+100 GeV, Min Bias (Q^{2} > 10 GeV^{2}), #Lambda_{c} in #eta [%.1f, %.1f], p_{T} [%.1f, %.1f] GeV",eta_lo[ieta],eta_hi[ieta],ppt_lo[ipt],ppt_hi[ipt]));
            pave->SetTextAlign(21);
            pave->SetTextFont(63);
            pave->SetTextSize(30);
            pave->Draw();

            TPaveText* pave4 = new TPaveText(0.04,0.16,0.70,0.18,"NDCNB");
            pave4->AddText(Form("%s",pid_full_name[pid_option]));
            pave4->SetTextAlign(11);
            pave4->SetTextFont(63);
            pave4->SetTextSize(23);
            pave4->Draw();

            gROOT->ProcessLine( Form("c1->Print(\"figs/Lc_decay_eta%d_pt%d.pdf\")",ieta,ipt) );
        }
      }
    }   
  }
}

void plot_decay_kin(const int smear_option = 2, const int bfield_option = 0, const int pid_option = 2)
{
  mcs(-1);

  TFile* fin[smearbin][bfieldbin][pidbin] = {0};
  
  for (int ismear = 0; ismear < smearbin; ++ismear)
  {
    for (int ibfield = 0; ibfield < bfieldbin; ++ibfield)
    {
      for (int ipid = 0; ipid < pidbin; ++ipid)
      {
        if (ismear!=smear_option) continue;
        if (ibfield!=bfield_option) continue;
        if (ipid!=pid_option) continue;

        fin[ismear][ibfield][ipid] = new TFile(Form("hists_highQ2_%s_%s_%s.root",smear_abbr[ismear],bfield_abbr[ibfield],pid_abbr[ipid]),"READ"); 
        fin[ismear][ibfield][ipid]->ls();
        
        for (int ieta = 0; ieta < etabin; ++ieta)
        {
          for (int ipt = 0; ipt < pptbin; ++ipt)
          {
            cout << "ieta " << ieta << " ipt " << ipt << endl; 
            h2d_K_D0_p_vs_eta[ismear][ibfield][ipid][ieta][ipt] = (TH2D*)fin[ismear][ibfield][ipid]->Get(Form("h2d_K_D0_p_vs_eta_%d_%d",ieta,ipt));
            h2d_pi_D0_p_vs_eta[ismear][ibfield][ipid][ieta][ipt] = (TH2D*)fin[ismear][ibfield][ipid]->Get(Form("h2d_pi_D0_p_vs_eta_%d_%d",ieta,ipt));

            h2d_K_Lc_p_vs_eta[ismear][ibfield][ipid][ieta][ipt] = (TH2D*)fin[ismear][ibfield][ipid]->Get(Form("h2d_K_Lc_p_vs_eta_%d_%d",ieta,ipt));
            h2d_pi_Lc_p_vs_eta[ismear][ibfield][ipid][ieta][ipt] = (TH2D*)fin[ismear][ibfield][ipid]->Get(Form("h2d_pi_Lc_p_vs_eta_%d_%d",ieta,ipt));    
            h2d_p_Lc_p_vs_eta[ismear][ibfield][ipid][ieta][ipt] = (TH2D*)fin[ismear][ibfield][ipid]->Get(Form("h2d_p_Lc_p_vs_eta_%d_%d",ieta,ipt));    
          }
        }
      }
    }
  }

  set_tracking_limits();

  set_PID_limits(1);

  plot_D0_decay(smear_option, bfield_option, pid_option);

  plot_Lc_decay(smear_option, bfield_option, pid_option);
}