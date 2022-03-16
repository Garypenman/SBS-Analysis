#include "GMnTree.C"
#include <TROOT.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <iostream>
#include <fstream>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TSystem.h>

using namespace std;

// options
const Bool_t   ApplyFidu  = true;
const Bool_t   ApplyElec  = true;
const Bool_t   ApplyElas  = false;
const Bool_t   ApplyPion  = false;
const Bool_t   Apply1Bar  = false;
const Bool_t   DoFit      = true;

const Bool_t   PlotHodo   = false;
const Bool_t   PlotBBCal  = true;
const Bool_t   PlotHCal   = true;
const Bool_t   PlotKine   = true;
const Bool_t   PlotOptics = false;
const Bool_t   PlotDiff   = false;
const Bool_t   PlotfADC   = false;

TF1 *f1, *f2;

double finter(double* x, double* par){
  return fabs(f1->EvalPar(x,par) - f2->EvalPar(x,par));
}

void Quasielastics(const Int_t kin_no = 8) { 

  //-----------------------------------------------------------------------------------------------------------------------------

  TChain* C = new TChain("T");

  Double_t Eb, th_sbs, th_bb, pcent, pres, runtime, avI;
  Double_t pdiff_off, hcal_dist;
  
  //Elastic Cuts Based on Kinematic Setting
  Double_t sh_min, sh_max, sh_e, ps_min, W_min, W_max, W_minc, W_maxc;
  Double_t hcal_ysig, hcal_ymean, hcal_ycut;
  Double_t hcal_xsigp, hcal_xsign,  hcal_xmeanp, hcal_xmeann, hcal_xcutp, hcal_xcutn;
  Double_t pdiffcut;

  sh_e = 0.70;
  
  W_min   = 0.0; 
  W_max   = 4.0; 
  W_minc = 0.25;
  W_maxc = 1.5;
  
  if( kin_no == 4) { //SBS-4
    C->Add("$OUT_DIR/LD2/*11593*.root");
    
    Eb        = 3.7278;  
    th_bb     = 36.0; 
    th_sbs    = 31.9; 
    hcal_dist = 11.0;
    
    pcent     = 2.122;  
    pres      = 0.02;   
    
    runtime   = 59. * 60.;    
    avI       = 1.75;     

    pdiff_off = 0.;

    sh_min  = 0.55;
    sh_max  = 0.85;
    ps_min  = 0.1;

    hcal_xmeanp = -0.425;
    hcal_xsigp = 0.193;
    
    hcal_ymean = -0.40;
    hcal_ysig = 0.15;
    
    hcal_xmeann = -1.575;
    hcal_xsign = 0.175;
    
    pdiffcut = 0.1;

    W_minc = 0.0;
    W_maxc = 4.0;

  }  
  else if( kin_no == 7) { //SBS-7
    //Need full LD2 Runs
    //C->Add("");
    
    Eb      = 7.906;   
    th_bb   = 40.0;   
    th_sbs  = 16.1; 
    hcal_dist = 14.0;
    
    pcent   = 2.670;  
    pres    = 0.02;   
    
    runtime = 0. * 60.;
    avI = 0.;
      
    pdiff_off = 0.23;

    sh_min  = 0.60;
    sh_max  = 0.95;
    ps_min  = 0.10;
      
    hcal_xmeanp = -1.0;
    hcal_xsigp = 0.1;
    
    hcal_ymean = -0.2;
    hcal_ysig = 0.1;
      
    hcal_xmeann = -1.0;
    hcal_xsign = 0.1;
    
    pdiffcut = 0.05;

    //W_minc = 0.0;
    //W_maxc = 4.0;
  }
  else if( kin_no == 11) { //SBS-11
    //Need full LD2 Runs
    //C->Add("");

    Eb      = 9.91;   
    th_bb   = 42.0; 
    th_sbs  = 13.3; 
    hcal_dist = 14.5;

    pcent   = 2.670;  
    pres    = 0.02;   
    
    runtime = 0. * 60.;
    avI     = 0.;
    
    pdiff_off = 0.23;

    sh_min  = 0.75;
    sh_max  = 1.05;
    ps_min  = 0.07;

    hcal_xmeanp = -1.0;
    hcal_xsigp = 0.1;
   
    hcal_ymean = -0.2;
    hcal_ysig = 0.1;
    
    hcal_xmeann = -1.0;
    hcal_xsign = 0.1;
    
    pdiffcut = 0.1;
  }
  else if( kin_no == 14) { //SBS-14
    //Need full LD2 Runs
    //C->Add("");
    
    Eb      = 5.9648;   
    th_bb   = 46.5; 
    th_sbs  = 17.3; 
    hcal_dist = 14.;
    
    pcent   = 2.0;  
    pres    = 0.02;   
    
    runtime = 0. * 60.;
    avI = 0.;
    
    pdiff_off = 0.23;

    sh_min  = 0.75;
    sh_max  = 1.05;
    ps_min  = 0.07;
    
    hcal_xmeanp = -0.75;
    hcal_xsigp = 0.125;
    
    hcal_ymean = -0.4;
    hcal_ysig = 0.125;
    
    hcal_xmeann = -1.0;
    hcal_xsign = 0.1;
    
    pdiffcut = 0.1;
  }
  else if( kin_no == 8) { //SBS-8
    //C->Add("$OUT_DIR/LD2/e1209019_fullreplay_13545_stream0_seg8*.root");
    C->Add("$OUT_DIR/LD2/e1209019_fullreplay_13545*.root"); 
    
    Eb      = 6.0;   
    th_bb   = 26.5; 
    th_sbs  = 29.9; 
    hcal_dist = 11.0;
    
    pcent   = 3.59;  
    pres    = 0.02;   
    
    runtime = 63. * (52./146.) * 60.;
    avI     = 5.0;
    
    pdiff_off = -0.03;

    sh_min  = 0.70;
    sh_max  = 1.05;
    ps_min  = 0.07;
      
    hcal_xmeanp = -0.33;
    hcal_xsigp = 0.2;
    
    hcal_ymean = -0.171;
    hcal_ysig = 0.155;
    
    hcal_xmeann = -1.48;
    hcal_xsign = 0.2;
    
    pdiffcut = 0.1;
   }
  else if( kin_no == 9) { //SBS-9
    //Need full LD2 Runs
    //C->Add("");
    
    Eb      = 4.014;   
    th_bb   = 49.0; 
    th_sbs  = 22.5; 
    hcal_dist = 11.0;
    
    pcent   = 1.63;  
    pres    = 0.02;   
    
    runtime = 0. * 60.; 
    avI     = 0.;
    
    pdiff_off = 0.29;

    sh_min  = 0.70;
    sh_max  = 0.95;
    ps_min  = 0.12;
    
    hcal_xmeanp = -0.634;
    hcal_xsigp = 0.094;
			
    hcal_ymean = -0.460;
    hcal_ysig = 0.155;
    
    hcal_xmeann = -1.0;
    hcal_xsign = 0.1;
    
    pdiffcut = 0.2;
  }
  
  GMnTree* T = new GMnTree(C);
  Long64_t nentries = C->GetEntries();
  cout << "Processing " << nentries << endl;

  //-----------------------------------------------------------------------------------------------------------------------------

  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  gStyle->SetCanvasColor(0);
  gStyle->SetFrameFillColor(0);

  gStyle->SetPadTopMargin(.05);
  gStyle->SetPadLeftMargin(.18);
  gStyle->SetPadRightMargin(.18);
  gStyle->SetPadBottomMargin(.15);

  gStyle->SetTitleOffset(1.1, "X");
  gStyle->SetTitleOffset(1.5, "Y");
  gStyle->SetTitleFont(42,"X");
  gStyle->SetTitleFont(42,"Y");
  gStyle->SetTitleSize(0.055,"X");
  gStyle->SetTitleSize(0.055,"Y");

  gStyle->SetLabelOffset(0.01, "X");
  gStyle->SetLabelOffset(0.01, "Y");
  gStyle->SetLabelFont(42,"X");
  gStyle->SetLabelFont(42,"Y");
  gStyle->SetLabelSize(0.045,"X");
  gStyle->SetLabelSize(0.045,"Y");

  gStyle->SetNdivisions(105,"X");
  gStyle->SetNdivisions(105,"Y");

  gStyle->SetStripDecimals(kFALSE);

  TH1D* hth_tmult   = new TH1D("hth_tmult","",6,0,6);
  TH1D* hth_hmult   = new TH1D("hth_cmult","",6,0,6);
  TH1D* hth_csize   = new TH1D("hth_csize","",6,0,6);

  TH1D* hth_tx      = new TH1D("hth_tx","",90,-1.125,1.125);
  TH1D* hth_ty      = new TH1D("hth_ty","",90,-0.30,0.30);

  TH1D* hth_xmean   = new TH1D("hth_xmean","",90,-1.125,1.125);
  TH1D* hth_xdiff   = new TH1D("hth_xdiff","",100,-0.2,0.2);
  TH1D* hth_ymean   = new TH1D("hth_ymean","",90,-0.30,0.30);
  TH1D* hth_ydiff   = new TH1D("hth_ydiff","",100,-0.2,0.2);

  TH2D* hth2d_txy    = new TH2D("hth2d_txy","", 1,-0.30,0.30, 90,-1.125,1.125);
  TH2D* hth2d_xymean = new TH2D("hth2d_xymean","", 1,-0.30,0.30, 90,-1.125,1.125);
  TH2D* hth2d_xdiff  = new TH2D("hth2d_xdiff","",20,-0.08,0.08,90,0,90);
  TH2D* hth2d_ydiff  = new TH2D("hth2d_ydiff","",20,-0.16,0.16,90,0,90);
  TH2D* hth2d_tdiff  = new TH2D("hth2d_tdiff","",100,-100.,100.,90,0,90);
  TH2D* hth2d_tmean  = new TH2D("hth2d_tmean","",100,-20.,20.,90,0,90);
  TH2D* hth2d_Diff   = new TH2D("hth2d_Diff", "", 50, -0.35, 0.35, 50, -5., 5. );

  TH2D* hth2d_fadcr1  = new TH2D("hth2d_fadcr1","",125,0.,250.,32,0,32);
  TH2D* hth2d_fadcr2  = new TH2D("hth2d_fadcr2","",125,0.,250.,32,0,32);
  TH2D* hth2d_fadcl1  = new TH2D("hth2d_fadcl1","",125,0.,250.,32,0,32);
  TH2D* hth2d_fadcl2  = new TH2D("hth2d_fadcl2","",125,0.,250.,32,0,32);
  TH1D* hth_fadcl1    = new TH1D("hth_fadcl1","",125,0.,250.);
  TH1D* hth_fadcl2    = new TH1D("hth_fadcl2","",125,0.,250.);
  TH1D* hth_fadcr1    = new TH1D("hth_fadcr1","",125,0.,250.);
  TH1D* hth_fadcr2    = new TH1D("hth_fadcr2","",125,0.,250.);

  TH1D* hkin_p        = new TH1D("hkin_p","",100,0.25*pcent,1.25*pcent);
  TH1D* hkin_th       = new TH1D("hkin_th","",100,-0.3,0.3);
  TH1D* hkin_ph       = new TH1D("hkin_ph","",100,-0.1,0.1);
  TH1D* hkin_x        = new TH1D("hkin_x","",100,-1.0,1.0);
  TH1D* hkin_y        = new TH1D("hkin_y","",100,-0.4,0.4);
  TH1D* hkin_yt       = new TH1D("hkin_yt","",100,-0.15,0.15);
  TH1D* hkin_pdiff    = new TH1D("hkin_pdiff","",50,-0.5,0.5);
  TH1D* hkin_W        = new TH1D("hkin_W","",50,W_min,W_max);
  TH2D* hkin2d_thp    = new TH2D("hkin2d_thp","",100, th_bb-6,th_bb+6.,100,0.25*pcent,1.25*pcent);
  
  TH1D* hkin_pc       = new TH1D("hkin_pc","",100,0.25*pcent,1.25*pcent);
  TH1D* hkin_thc      = new TH1D("hkin_thc","",100,-0.3,0.3);
  TH1D* hkin_phc      = new TH1D("hkin_phc","",100,-0.1,0.1);
  TH1D* hkin_xc       = new TH1D("hkin_xc","",100,-1.0,1.0);
  TH1D* hkin_yc       = new TH1D("hkin_yc","",100,-0.4,0.4);
  TH1D* hkin_ytc      = new TH1D("hkin_ytc","",100,-0.15,0.15);
  TH1D* hkin_pdiffc   = new TH1D("hkin_pdiffc","",50,-0.5,0.9);
  TH1D* hkin_Wc       = new TH1D("hkin_Wc","",50,W_minc,W_maxc);
  TH2D* hkin2d_thpc    = new TH2D("hkin2d_thpc","",100, th_bb-6,th_bb+6.,100,0.25*pcent,1.25*pcent);
  
  TH1D *hbbcal_psE    = new TH1D("hbbcal_psE","",100,0,1.25*pcent/2.); 
  TH1D *hbbcal_shE    = new TH1D("hbbcal_shE","",100,0,1.25*pcent); 
  TH1D* hbbcal_cale   = new TH1D("hbbcal_cale","",100,0.25*pcent,1.25*pcent);
  TH1D* hbbcal_edivp  = new TH1D("hbbcal_edivp","",100,-1.0,1.0);
  TH1D* hbbcal_xdiff  = new TH1D("hbbcal_xdiff","",100,-0.2,0.2);
  TH1D* hbbcal_ydiff  = new TH1D("hbbcal_ydiff","",100,-0.2,0.2);
  TH2D *hbbcal2d_pss  = new TH2D("hbbcal2d_pssh","",100,0.,1.2, 100,0.,1.2); 
  TH1D *hbbcal_psEc   = new TH1D("hbbcal_psEc","",100,0,1.25*pcent/2.); 
  TH1D *hbbcal_shEc   = new TH1D("hbbcal_shEc","",100,0,1.25*pcent); 
  TH1D* hbbcal_calec  = new TH1D("hbbcal_calec","",100,0.25*pcent,1.25*pcent);
  TH1D* hbbcal_edivpc = new TH1D("hbbcal_edivpc","",100,-1.0,1.0);
  TH1D* hbbcal_xdiffc = new TH1D("hbbcal_xdiffc","",100,-0.2,0.2);
  TH1D* hbbcal_ydiffc = new TH1D("hbbcal_ydiffc","",100,-0.2,0.2);
  TH2D *hbbcal2d_pssc  = new TH2D("hbbcal2d_psshc","",100,0.,1.2, 100,0.,1.2); 

  TH2D* hhcal_xbb   = new TH2D("hhcal_xbb","",100,-1.2,1.2,100,-1.2,1.2);
  TH2D* hhcal_ybb   = new TH2D("hhcal_ybb","",100,-0.3,0.3,100,-0.3,0.3);
  TH1D* hhcal_xdiff = new TH1D("hhcal_xdiff","",100,-1.2,1.2);
  TH1D* hhcal_ydiff = new TH1D("hhcal_ydiff","",100,-1.2,1.2);
  TH1D* hhcal_xdiffc = new TH1D("hhcal_xdiffc","",100,-1.2,1.2);
  TH1D* hhcal_xdiffcW = new TH1D("hhcal_xdiffcW","",100,-1.2,1.2);
  TH1D* hhcal_ydiffc = new TH1D("hhcal_ydiffc","",100,-1.2,1.2);

  TH1D* hhcal_x = new TH1D("hhcal_x","",100,-2.5,2.5);
  TH1D* hhcal_xc = new TH1D("hhcal_xc","",100,-2.5,2.5);
  TH1D* hhcal_y = new TH1D("hhcal_y","",100,-2.,2.);
  TH1D* hhcal_yc = new TH1D("hhcal_yc","",100,-2.,2.);
  TH2D* hhcal_xy = new TH2D("hhcal_xy","",50,-2.,2.,100,-2.5,1.);
  TH2D* hhcal_xyc = new TH2D("hhcal_xyc","",50,-2.,2.,100,-2.5,1.);
  
  TH1D* hhcal_predx = new TH1D("hhcal_predx","",100,-2.5,2.5);
  TH1D* hhcal_predxc = new TH1D("hhcal_predxc","",100,-2.5,2.5);
  TH1D* hhcal_predy = new TH1D("hhcal_predy","",100,-2.,2.);
  TH1D* hhcal_predyc = new TH1D("hhcal_predyc","",100,-2.,2.);
  TH2D* hhcal_predxy = new TH2D("hhcal_predxy","",50,-2.,2.,50,-2.,2.);
  TH2D* hhcal_predxyc = new TH2D("hhcal_predxyc","",50,-2.,2.,50,-2.,2.);
  
  TH1D* hhcal_deltax = new TH1D("hhcal_deltax","",100,-2.5,2.5);
  TH1D* hhcal_deltaxc = new TH1D("hhcal_deltaxc","",100,-2.5,2.5);
  TH1D* hhcal_deltay = new TH1D("hhcal_deltay","",100,-2,2.);
  TH1D* hhcal_deltayc = new TH1D("hhcal_deltayc","",100,-2.,2.);
  TH2D* hhcal_deltaxy = new TH2D("hhcal_deltaxy","",50,-2.,2.,50,-2.5,2.5);
  TH2D* hhcal_deltaxyc = new TH2D("hhcal_deltaxyc","",50,-2.5,2.5,50,-2.5,2.5);
  TH2D* hkin_deltaxyc = new TH2D("hkin_deltaxyc","",50,-2.5,2.5,50,-2.5,2.5);
  
  TH1D* hbbcal_pxdiff = new TH1D("hbbcal_pxdiff","",100,-0.5,0.5);
  TH1D* hbbcal_pydiff = new TH1D("hbbcal_pydiff","",100,-0.5,0.5);
  TH1D* hbbcal_pzdiff = new TH1D("hbbcal_pzdiff","",100,-0.5,0.5);
  
  TH2D* hhcalsh_blkdiff = new TH2D("hhcalbbsh_blkdiff","",5,0,5,5,0,5);
  TH2D* hhcalps_blkdiff = new TH2D("hhcalbbps_blkdiff","",5,0,5,5,0,5);
  
  TH2D* hopt_x_pd     = new TH2D("hopt_x_pd","",  20,-0.8,0.8,   20,-0.1, +0.1);
  TH2D* hopt_y_pd     = new TH2D("hopt_y_pd","",  20,-0.3,0.3,   20,-0.1, +0.1);
  TH2D* hopt_xp_pd     = new TH2D("hopt_xp_pd","",20,-0.25,0.25, 20,-0.1, +0.1);
  TH2D* hopt_yp_pd     = new TH2D("hopt_yp_pd","",20,-0.1,0.1,   20,-0.1, +0.1);
  TH2D* hopt_p_pd      = new TH2D("hopt_p_pd","",20,pcent-0.2,pcent+0.2,    20,-0.1, +0.1);
  TH2D* hopt_th_pd     = new TH2D("hopt_th_pd","",20,-0.3,0.3,   20,-0.1, +0.1);
  TH2D* hopt_ph_pd     = new TH2D("hopt_ph_pd","",20,-0.1,0.1,   20,-0.1, +0.1);
  TH2D* hopt_yt_pd     = new TH2D("hopt_yt_pd","",20,-0.1,0.1,   20,-0.1, +0.1);

  TH2D* hopt_x_pc     = new TH2D("hopt_x_pc","",  20,-0.8,0.8,   20,pcent-0.2,pcent+0.2);
  TH2D* hopt_y_pc     = new TH2D("hopt_y_pc","",  20,-0.3,0.3,   20,pcent-0.2,pcent+0.2);
  TH2D* hopt_xp_pc     = new TH2D("hopt_xp_pc","",20,-0.25,0.25, 20,pcent-0.2,pcent+0.2);
  TH2D* hopt_yp_pc     = new TH2D("hopt_yp_pc","",20,-0.1,0.1,   20,pcent-0.2,pcent+0.2);
  TH2D* hopt_p_pc      = new TH2D("hopt_p_pc","",20,pcent-0.2,pcent+0.2,   20,pcent-0.2,pcent+0.2);
  TH2D* hopt_th_pc     = new TH2D("hopt_th_pc","",20,-0.3,0.3,   20,pcent-0.2,pcent+0.2);
  TH2D* hopt_ph_pc     = new TH2D("hopt_ph_pc","",20,-0.1,0.1,   20,pcent-0.2,pcent+0.2);
  TH2D* hopt_yt_pc     = new TH2D("hopt_yt_pc","",20,-0.1,0.1,   20,pcent-0.2,pcent+0.2);

  const Int_t    nBarsTDC   = 90;
  TH1D* hResx[nBarsTDC];
  TH1D* hResy[nBarsTDC];
  TH2D* hDiff[nBarsTDC];
  
  for(Int_t i=0; i<nBarsTDC; i++) {
    hDiff[i] = new TH2D(Form("hDiff_%d",i), "", 50, -0.35, 0.35, 50, -7., 7. );
  }
  
  Double_t Mp = 0.93827;
  TLorentzVector Tp4(0,0,0,Mp); //target 4vec
  TLorentzVector kp4(0,0,Eb,Eb); //beam 4vec
  TLorentzVector Qp4, kpp4, Rp4; //q, recoil electron, recoil nucleon

  th_bb = th_bb * M_PI/180.;
  th_sbs = th_sbs * M_PI/180.;
  avI = avI * 0.85;
  
  //-----------------------------------------------------------------------------------------------------------------------------

  for(Long64_t ev=0; ev<nentries;ev++) {

    T->GetEntry(ev);
    
    if( ev%10000 == 0 )
      cout << ev << endl;

    //-----------------------------------------------------------------------------------------------------------------------------
    // Pre-cuts -- these should be in the replay cdef
    //-----------------------------------------------------------------------------------------------------------------------------

    if( T->bb_tr_n <= 0 ) continue; 
    if( T->bb_ps_e < 0.05 ) continue;

    //-----------------------------------------------------------------------------------------------------------------------------
    // Hodo fADCs (no cuts)
    //-----------------------------------------------------------------------------------------------------------------------------

//     for(int i = 0 ; i < T->Ndata_bb_hodoadc_bar_adc_L_ap ; i++ ) {
//       hth2d_fadcl1->Fill(T->bb_hodoadc_bar_adc_L_ap[i], T->bb_hodoadc_bar_adc_id[i]);
// 	hth_fadcl1->Fill(T->bb_hodoadc_bar_adc_L_ap[i]);
//     }
//     for(int i = 0 ; i < T->Ndata_bb_hodoadc_bar_adc_R_ap ; i++ )  {
//       hth2d_fadcr1->Fill(T->bb_hodoadc_bar_adc_R_ap[i], T->bb_hodoadc_bar_adc_id[i]);
// 	hth_fadcr1->Fill(T->bb_hodoadc_bar_adc_R_ap[i]);
//     }

    //-----------------------------------------------------------------------------------------------------------------------------
    // Kinematic cuts
    //-----------------------------------------------------------------------------------------------------------------------------

    const Double_t hodo_dist  = 1.8545;
    const Double_t show_dist  = 1.902;

    Double_t Mp      = 0.93827;
    Double_t th      = acos(T->bb_tr_pz[0]/T->bb_tr_p[0]);
    Double_t pexp_th = 2.*Mp*Eb*(Mp+Eb)*cos(th) / (Mp*Mp + 2.*Mp*Eb + (Eb*sin(th)*Eb*sin(th))); // e mom from angle
    Double_t p       = T->bb_tr_p[0];
    Double_t pdiff   = p - pexp_th; 
      
    pdiff = pdiff - pdiff_off;
    pdiff = pdiff-0.07*T->bb_tr_x[0];
    pdiff = pdiff-0.32*T->bb_tr_y[0];
    
    Double_t pc = p;
    
    Double_t px1 = pc * TMath::Sin( th_bb+T->bb_tr_tg_ph[0] ) * TMath::Cos( T->bb_tr_tg_th[0] );
    Double_t py1 = pc * TMath::Sin( th_bb+T->bb_tr_tg_ph[0] ) * TMath::Sin( T->bb_tr_tg_th[0] );
    Double_t pz1 = pc * TMath::Cos( th_bb+T->bb_tr_tg_ph[0] );

    // _tg = target plane,  tr_ = local transport coord system (focal plane variables) 
    Double_t px = T->bb_tr_px[0];
    Double_t py = T->bb_tr_py[0];
    Double_t pz = T->bb_tr_pz[0];
    
    hbbcal_pxdiff->Fill(px-px1);
    hbbcal_pydiff->Fill(py+py1);
    hbbcal_pzdiff->Fill(pz-pz1);
    
    hhcalsh_blkdiff->Fill(T->Ndata_bb_sh_clus_blk_id,T->Ndata_sbs_hcal_clus_blk_id);
    hhcalps_blkdiff->Fill(T->Ndata_bb_ps_clus_blk_id,T->Ndata_sbs_hcal_clus_blk_id);
    
    kpp4.SetPxPyPzE(px,py,pz,pc);
    Qp4 = kp4 - kpp4;
    Rp4 = Tp4 + Qp4;
    Double_t W2 = Rp4.M2(); 

    Double_t trx_sh  = (T->bb_tr_x[0] + (show_dist) * T->bb_tr_th[0] );
    Double_t try_sh  = (T->bb_tr_y[0] + (show_dist) * T->bb_tr_ph[0] );
    
    Rp4.RotateY(th_sbs);
    
    Double_t hcal_th = TMath::ATan(Rp4.Px()/Rp4.Pz());
    Double_t hcal_ph = TMath::ATan(Rp4.Py()/Rp4.Pz());

    //offsets gained from adjusting delta peaks to 0 with magnet off.
    Double_t hcal_xoff = 0.;
    Double_t hcal_yoff = 0.; 
    
    Double_t hcal_x = T->sbs_hcal_x;
    Double_t pred_x = (hcal_dist * TMath::Sin( hcal_ph )) + hcal_xoff;
    
    Double_t hcal_y = T->sbs_hcal_y;
    Double_t pred_y   = -hcal_dist * TMath::Sin(hcal_th) + hcal_yoff;
    
    Double_t delta_x = hcal_x - pred_x;
    Double_t delta_y = hcal_y - pred_y;
    
    if( ApplyFidu ) {
      Double_t thtg_cut   = 0.16;
      Double_t phtg_cut   = 0.07;
      Double_t ytg_cut    = 0.05;
      if( fabs(T->bb_tr_th[0])>0.1  ) continue;
      if( fabs(T->bb_tr_tg_th[0]) > thtg_cut ) continue;
      if( fabs(T->bb_tr_tg_ph[0]) > phtg_cut ) continue;
      if( fabs(T->bb_tr_tg_y[0]) > ytg_cut ) continue;
    }
   
    hhcal_xdiff->Fill( trx_sh - T->sbs_hcal_y);
    hhcal_ydiff->Fill( try_sh - T->sbs_hcal_x);
    
    hkin_p->Fill(p);
    hkin_x->Fill(T->bb_tr_x[0]);
    hkin_y->Fill(T->bb_tr_y[0]);
    hkin_th->Fill(T->bb_tr_tg_th[0]);
    hkin_ph->Fill(T->bb_tr_tg_ph[0]);
    hkin_yt->Fill(T->bb_tr_tg_y[0]);
    hkin_pdiff->Fill( pdiff );
    hkin_W->Fill(T->e_kine_W2);
    hkin2d_thp->Fill(57.3*th, p);

    hbbcal_psE->Fill( T->bb_ps_e );
    hbbcal_shE->Fill( T->bb_sh_e );
    hbbcal_cale->Fill((T->bb_ps_e + T->bb_sh_e));
    hbbcal_edivp->Fill(((T->bb_ps_e + T->bb_sh_e) - p));
    hbbcal_xdiff->Fill( T->bb_sh_x - trx_sh );
    hbbcal_ydiff->Fill( T->bb_sh_y - try_sh );
    hbbcal2d_pss->Fill( T->bb_sh_e/T->bb_tr_p[0], T->bb_ps_e/T->bb_tr_p[0] );
    
    hhcal_xbb->Fill( trx_sh,  T->sbs_hcal_x);
    hhcal_ybb->Fill( try_sh, T->sbs_hcal_y);
    hhcal_xdiffc->Fill( trx_sh - T->sbs_hcal_x);
    hhcal_ydiffc->Fill( try_sh - T->sbs_hcal_y);
    
    hhcal_x->Fill(hcal_x);
    hhcal_y->Fill(hcal_y);
    hhcal_xy->Fill(hcal_y,hcal_x);
    
    hhcal_predx->Fill(pred_x);
    hhcal_predy->Fill(pred_y);
    hhcal_predxy->Fill(pred_y,pred_x);
    
    hhcal_deltax->Fill(delta_x);
    hhcal_deltay->Fill(delta_y);
    hhcal_deltaxy->Fill(delta_y,delta_x);
    
    hcal_xcutp = 2. * hcal_xsigp;
    hcal_ycut = 3. * hcal_ysig;
    hcal_xcutn = 2. * hcal_xsign;
    
    if( ApplyElec ) { 
      if( T->bb_ps_e/T->bb_tr_p[0] < ps_min ) continue;
      if( (T->bb_ps_e/T->bb_tr_p[0] + sh_e*T->bb_sh_e/T->bb_tr_p[0])< sh_min ) continue;
      if( (T->bb_ps_e/T->bb_tr_p[0] + sh_e*T->bb_sh_e/T->bb_tr_p[0])> sh_max) continue;
      if( T->bb_tr_p[0] < 0.25*pcent ) continue;
      if( T->e_kine_W2 < W_min ) continue; 
      if( T->e_kine_W2 > W_max ) continue;   

    }
    else if (ApplyPion) { 
      if( T->bb_ps_e/T->bb_tr_p[0] > ps_min ) continue;
      if( (T->bb_ps_e/T->bb_tr_p[0] + sh_e*T->bb_sh_e/T->bb_tr_p[0]) > sh_min ) continue;
    }
    
    if( fabs(delta_x - hcal_xmeanp) < hcal_xcutp ){
      if( fabs(delta_y - hcal_ymean) < hcal_ycut ){
	
	hkin_pdiffc->Fill( pdiff );
	
	//if( T->e_kine_W2 < W_minc ) continue; 
	//if( T->e_kine_W2 > W_maxc ) continue;   
	hkin_pc->Fill(pc);
	hkin_xc->Fill(T->bb_tr_x[0]);
	hkin_yc->Fill(T->bb_tr_y[0]);
	hkin_thc->Fill(T->bb_tr_tg_th[0]);
	hkin_phc->Fill(T->bb_tr_tg_ph[0]);
	hkin_ytc->Fill(T->bb_tr_tg_y[0]);
	hkin2d_thpc->Fill(57.3*th, pc);
	hkin_Wc->Fill(T->e_kine_W2);
	hkin_deltaxyc->Fill(delta_y,delta_x);
	
      }	
    }
    
    if( fabs(delta_x - hcal_xmeann) < hcal_xcutn){
      if( fabs(delta_y - hcal_ymean) < hcal_ycut ){
	
	//hkin_pdiffc->Fill( pdiff );
	
	//if( T->e_kine_W2 < W_minc ) continue; 
	//if( T->e_kine_W2 > W_maxc ) continue;   
	//hkin_pc->Fill(pc);
	//hkin_xc->Fill(T->bb_tr_x[0]);
	//hkin_yc->Fill(T->bb_tr_y[0]);
	//hkin_thc->Fill(T->bb_tr_tg_th[0]);
	//hkin_phc->Fill(T->bb_tr_tg_ph[0]);
	//hkin_ytc->Fill(T->bb_tr_tg_y[0]);
	//hkin2d_thpc->Fill(57.3*th, pc);
	//hkin_Wc->Fill(T->e_kine_W2);
	hkin_deltaxyc->Fill(delta_y,delta_x);
	
      }	
    }
    
    
    //if( T->e_kine_W2 < W_minc ) continue; 
    //if( T->e_kine_W2 > W_maxc ) continue;   
    if( fabs(pdiff) < pdiffcut ) {
    //if( true ){
      hhcal_xc->Fill(hcal_x);
      hhcal_yc->Fill(hcal_y);
      hhcal_xyc->Fill(hcal_y,hcal_x);
    
      hhcal_predxc->Fill(pred_x);
      hhcal_predyc->Fill(pred_y);
      hhcal_predxyc->Fill(pred_y,pred_x);
    
      hhcal_deltaxc->Fill(delta_x);
      hhcal_deltayc->Fill(delta_y);
      hhcal_deltaxyc->Fill(delta_y,delta_x);
      
      hbbcal_psEc->Fill( T->bb_ps_e );
      hbbcal_shEc->Fill( T->bb_sh_e );
      hbbcal_calec->Fill((T->bb_ps_e + T->bb_sh_e));
      hbbcal_edivpc->Fill(((T->bb_ps_e + T->bb_sh_e) - pc));
      hbbcal_xdiffc->Fill( T->bb_sh_x - (T->bb_tr_x[0] + (show_dist) * T->bb_tr_th[0] ) ); 
      hbbcal_ydiffc->Fill( T->bb_sh_y - (T->bb_tr_y[0] + (show_dist) * T->bb_tr_ph[0] ) );
      hbbcal2d_pssc->Fill( T->bb_sh_e/T->bb_tr_p[0], T->bb_ps_e/T->bb_tr_p[0] );
      
      if ( PlotOptics ){
	hopt_x_pd->Fill(T->bb_tr_x[0] , pdiff);
	hopt_y_pd->Fill(T->bb_tr_y[0] , pdiff);
	hopt_xp_pd->Fill(T->bb_tr_th[0] , pdiff);
	hopt_yp_pd->Fill(T->bb_tr_ph[0] , pdiff);
	hopt_p_pd->Fill(p , pdiff);
	hopt_th_pd->Fill(T->bb_tr_tg_th[0] , pdiff);
	hopt_ph_pd->Fill(T->bb_tr_tg_ph[0] , pdiff);
	hopt_yt_pd->Fill(T->bb_tr_tg_y[0] , pdiff);
      
	hopt_x_pc->Fill(T->bb_tr_x[0] , p);
	hopt_y_pc->Fill(T->bb_tr_y[0] , p);
	hopt_xp_pc->Fill(T->bb_tr_th[0] , p);
	hopt_yp_pc->Fill(T->bb_tr_ph[0] , p);
	hopt_p_pc->Fill(p , p);
	hopt_th_pc->Fill(T->bb_tr_tg_th[0] , p);
	hopt_ph_pc->Fill(T->bb_tr_tg_ph[0] , p);
	hopt_yt_pc->Fill(T->bb_tr_tg_y[0] , p);
      }
    }
      
      //-----------------------------------------------------------------------------------------------------------------------------
      // GEM-Hodoscope track matching
      //-----------------------------------------------------------------------------------------------------------------------------
      if ( PlotHodo ){
	hth_tx->Fill( T->bb_tr_x[0] + hodo_dist * T->bb_tr_th[0] ); // track x at hodo (dispersive)
	hth_ty->Fill( T->bb_tr_y[0] + hodo_dist * T->bb_tr_ph[0] ); // track y at hodo (non-dispersive)
    
	hth2d_txy->Fill( T->bb_tr_y[0] + hodo_dist * T->bb_tr_ph[0], T->bb_tr_x[0] + hodo_dist * T->bb_tr_th[0] );
    
	hth_tmult->Fill(T->bb_tr_n-1); // BB track "id" 
    
	Bool_t Singlebarhit = kFALSE;
	if( T->bb_hodotdc_clus_size[0] == 1 && T->bb_hodotdc_clus_trackindex[0] == 0 )
	  Singlebarhit = kTRUE;
	
	if( Apply1Bar && !Singlebarhit ) continue;
    
	//-----------------------------------------------------------------------------------------------------------------------------
	
	if( T->bb_hodotdc_clus_trackindex[0] == 0 ) {  // only hodo clusters that match bb track id = 0 (I guess redundndat for single track events)	
      
	  hth_csize->Fill( T->bb_hodotdc_clus_size[0] );
	  hth_hmult->Fill( T->bb_hodotdc_clus_trackindex[0] ); // track id that is matched
      
	  hth_xmean->Fill( T->bb_hodotdc_clus_xmean[0] ); // mean x position of hodo cluster
	  hth_xdiff->Fill( (T->bb_hodotdc_clus_xmean[0] - (T->bb_tr_x[0] + hodo_dist * T->bb_tr_th[0])) );
      
	  hth_ymean->Fill( T->bb_hodotdc_clus_ymean[0] ); // mean y position of hodo cluster
	  hth_ydiff->Fill( (T->bb_hodotdc_clus_ymean[0] - (T->bb_tr_y[0] + hodo_dist * T->bb_tr_ph[0])) );
      
	  Int_t maxbar = (Int_t)T->bb_hodotdc_clus_id[0]; 
      
	  hth2d_xdiff->Fill( (T->bb_hodotdc_clus_xmean[0] - (T->bb_tr_x[0] + hodo_dist * T->bb_tr_th[0])), maxbar );
	  hth2d_ydiff->Fill( (T->bb_hodotdc_clus_ymean[0] - (T->bb_tr_y[0] + hodo_dist * T->bb_tr_ph[0])), maxbar );
	  hth2d_tmean->Fill( T->bb_hodotdc_clus_tmean[0], maxbar );

	  hth2d_xymean->Fill( T->bb_hodotdc_clus_ymean[0], T->bb_hodotdc_clus_xmean[0] );
      
	  hth2d_Diff->Fill((T->bb_tr_y[0] + hodo_dist * T->bb_tr_ph[0]), 0.5*T->bb_hodotdc_clus_tdiff[0] ); // note th 1/2
	  
	  hDiff[maxbar]->Fill((T->bb_tr_y[0] + hodo_dist * T->bb_tr_ph[0]), 0.5*T->bb_hodotdc_clus_tdiff[0] ); // note th 1/2
	}
      }
  }
  
  //-----------------------------------------------------------------------------------------------------------------------------
  // GEM-Hodoscope track matching plots
  //-----------------------------------------------------------------------------------------------------------------------------

  TLatex* tex;
  TLatex* tex2;
  
  if( PlotHodo ) {
    TCanvas* c1 = new TCanvas("c1","",1200,800);
    c1->Divide(4,2);
    
    c1->cd(1);
    hth_tx->Draw("");
    hth_tx->SetLineColor(2);
    hth_xmean->Draw("same");
    hth_tx->GetXaxis()->SetTitle("x_{FP} [m]");
    
    tex = new TLatex( 0.35, 0.3, "GEM");
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(2);
    tex->SetTextSize(0.055);
    tex->Draw();
    
    tex = new TLatex( 0.35, 0.22, "BBHodo");
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(4);
    tex->SetTextSize(0.055);
    tex->Draw();
    
    c1->cd(2);
    hth_ty->Draw("");
    hth_ty->SetLineColor(2);
    hth_ymean->Draw("same");
    hth_ty->GetXaxis()->SetTitle("y_{FP} [m]");
    
    tex = new TLatex( 0.35, 0.3, "GEM");
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(2);
    tex->SetTextSize(0.055);
    tex->Draw();
    
    tex = new TLatex( 0.35, 0.22, "BBHodo");
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(4);
    tex->SetTextSize(0.055);
    tex->Draw();
    
    c1->cd(3);
    
    hth_hmult->Divide(hth_tmult);
    hth_hmult->Draw("");
    hth_hmult->GetXaxis()->SetTitle("BB track ID");
    hth_csize->Draw("");
    hth_csize->SetLineColor(2);
    hth_csize->GetXaxis()->SetTitle("BBHodo Cluster Size");
    
    cout << "Mean cluster size " << hth_csize->GetMean() << endl;
    
    tex = new TLatex( 0.25, 0.8, Form("Mean Size = %3.2f bars", hth_csize->GetMean()) );
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.055);
    tex->Draw();
    
    tex = new TLatex( 0.25, 0.72, Form("Efficiency = %3.2f %%", 100*hth_hmult->GetBinContent(1)) );
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.055);
    tex->Draw();
        
    c1->cd(4);
    hth2d_xymean->Divide(hth2d_txy);
    hth2d_xymean->Draw("colz");
    hth2d_xymean->GetYaxis()->SetTitle("x_{FP} [m]");
    hth2d_xymean->GetXaxis()->SetTitle("y_{FP} [m]");
    hth2d_xymean->SetMaximum(1.3);
        
    c1->cd(5);
    TH1D* hth_xmean1 = (TH1D*)hth_xmean->Clone("hth_xmean1");
    hth_xmean->Sumw2();
    hth_xmean1->Sumw2();
    hth_tx->Sumw2();
    hth_xmean1->Divide(hth_tx);
    hth_xmean1->Draw("");
    hth_xmean1->GetXaxis()->SetTitle("x_{FP} [m]");
    hth_xmean1->GetYaxis()->SetTitle("Track Match Efficiency");
    hth_xmean1->SetMaximum(1.5);
    
    c1->cd(6);
    TH1D* hth_ymean1 = (TH1D*)hth_ymean->Clone("hth_ymean1");
    hth_ymean->Sumw2();
    hth_ymean1->Sumw2();
    hth_ty->Sumw2();
    hth_ymean1->Divide(hth_ty);
    hth_ymean1->Draw("");
    hth_ymean1->GetXaxis()->SetTitle("y_{FP} [m]");
    hth_ymean1->GetYaxis()->SetTitle("Track Match Efficiency");
    hth_ymean1->SetMaximum(1.5);
    
    c1->cd(7);
    hth2d_xdiff->Draw("colz");
    hth2d_xdiff->GetXaxis()->SetTitle("#delta x_{FP} (GEM-BBHodo) [m]");
    hth2d_xdiff->GetYaxis()->SetTitle("Hodoscope Bar ID");
    
    c1->cd(8);
    hth2d_ydiff->Draw("colz");
    hth2d_ydiff->GetXaxis()->SetTitle("#delta y_{FP} (GEM-BBHodo) [m]");
    hth2d_ydiff->GetYaxis()->SetTitle("Hodoscope Bar ID");
    
    c1->Print("bbhodo_temp1.pdf");
    c1->Close();
    //-----------------------------------------------------------------------------------------------------------------------------
    // Hodoscope bar efficiencies and resolution plots
    //-----------------------------------------------------------------------------------------------------------------------------

    TCanvas* c2 = new TCanvas("c2","",1200,800);
    c2->Divide(3,1);
    c2->cd(1)->SetRightMargin(.05);
    
    Double_t Eff[nBarsTDC], eEff[nBarsTDC];
    Double_t Bar[nBarsTDC], eBar[nBarsTDC];
    
    for(Int_t i=0; i<nBarsTDC; i++){
      Bar[i]  = (Double_t)(nBarsTDC -1 - i);
      eBar[i] = 0.0;
      Eff[i]  = hth_xmean1->GetBinContent(i);
      eEff[i]  = hth_xmean1->GetBinError(i);
    }
    
    TGraphErrors* gEff = new TGraphErrors( nBarsTDC, Bar, Eff, eBar, eEff );
    gEff->SetMarkerColor( 4 );
    gEff->SetLineColor( 4 );
    gEff->SetMarkerStyle( 20 );
    gEff->SetMarkerSize( 1.5 );
    gEff->Draw("AP");
    gEff->GetYaxis()->SetRangeUser( 0.0, 1.8 );
    gEff->GetXaxis()->SetTitle( "Hodoscope Bar ID");
    gEff->GetYaxis()->SetTitle( "Track Match Efficiency");
    gEff->RemovePoint(59);
    
    TF1* meaneff = new TF1("meaneff","pol0", 12., 78.);  
    meaneff->SetLineColor( 4 );
    gEff->Fit(meaneff,"Q","",12, 78);
    
    tex = new TLatex( 0.35, 0.9, Form("Mean = %3.2f %%",100.* meaneff->GetParameter(0)));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(4);
    tex->SetTextSize(0.055);
    tex->Draw();
      
    //========================================================= Resolutions
    
    
    TF1* gausf = new TF1("gausf","gaus(0)",-0.20, 0.20);
    
    Double_t ulim, dlim;
    Double_t sig1, sig2, esig1, esig2;
    Double_t sigx[nBarsTDC], esigx[nBarsTDC];
    Double_t sigy[nBarsTDC], esigy[nBarsTDC];
    
    for(Int_t i=0; i<nBarsTDC; i++){ 
      
      Bar[i]  = (Double_t)i;
      
      hResx[i] = (TH1D*)hth2d_xdiff->ProjectionX(Form("htempx%d",i),i,i+1);
      
      ulim = hResx[i]->GetMean() + (5*hResx[i]->GetRMS());
      dlim = hResx[i]->GetMean() - (5*hResx[i]->GetRMS()); 
      gausf->SetParameter(1, hResx[i]->GetMean() ); 
      gausf->SetParameter(2,(3*hResx[i]->GetRMS()) ); 
      
      hResx[i]->Fit(gausf,"Q","",dlim, ulim);

      sigx[i]  = gausf->GetParameter(2); 
      esigx[i] = gausf->GetParError(2);

      hResy[i] = (TH1D*)hth2d_ydiff->ProjectionX(Form("htempy%d",i),i,i+1);
      
      ulim = hResy[i]->GetMean() + (5*hResy[i]->GetRMS());
      dlim = hResy[i]->GetMean() - (5*hResy[i]->GetRMS()); 
      gausf->SetParameter(1, hResy[i]->GetMean() ); 
      gausf->SetParameter(2,(3*hResy[i]->GetRMS()) ); 
      
      hResy[i]->Fit(gausf,"Q","",dlim, ulim);

      sigy[i]  = gausf->GetParameter(2); 
      esigy[i] = gausf->GetParError(2);
    }    

    c2->cd(2)->SetRightMargin(.05);
  
    TGraphErrors* gResx = new TGraphErrors( nBarsTDC, Bar, sigx, eBar, esigx );
    gResx->SetMarkerColor( 2 );
    gResx->SetLineColor( 2 );
    gResx->SetMarkerStyle( 20 );
    gResx->SetMarkerSize( 1.5 );
    gResx->Draw("AP");
    gResx->GetXaxis()->SetTitle( "Hodoscope Bar ID");
    gResx->GetYaxis()->SetTitle( "Vertical Position Resolution [m]");
    gResx->GetYaxis()->SetRangeUser( 0.0, 0.03 );
    gResx->RemovePoint(59);
    gResx->RemovePoint(29);

    TF1* meanresx = new TF1("meanresx","pol0", 12., 78.);  
    meanresx->SetLineColor( 2 );
    gResx->Fit(meanresx,"Q","",12., 78.);
    
    tex = new TLatex( 0.34, 0.9, Form("Mean = %4.3f m", meanresx->GetParameter(0)));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(2);
    tex->SetTextSize(0.055);
    tex->Draw();

    c2->cd(3)->SetRightMargin(.05);
  
    TGraphErrors* gResy = new TGraphErrors( nBarsTDC, Bar, sigy, eBar, esigy );
    gResy->SetMarkerColor( 1 );
    gResy->SetLineColor( 1 );
    gResy->SetMarkerStyle( 20 );
    gResy->SetMarkerSize( 1.5 );
    gResy->Draw("AP");
    gResy->GetXaxis()->SetTitle( "Hodoscope Bar ID");
    gResy->GetYaxis()->SetTitle( "Horizontal Position Resolution [m]");
    gResy->GetYaxis()->SetRangeUser( 0.0, 0.06 );
    gResy->RemovePoint(59);
    gResy->RemovePoint(29);

    TF1* meanresy = new TF1("meanresy","pol0", 12., 78.);  
    meanresy->SetLineColor( 1 );
    gResy->Fit(meanresy,"Q","",12., 78.);
    
    tex = new TLatex( 0.34, 0.9, Form("Mean = %4.3f m", meanresy->GetParameter(0)));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.055);
    tex->Draw();
    
    c2->Print("bbhodo_temp2.pdf");
    c2->Close();
  }
  
  //-----------------------------------------------------------------------------------------------------------------------------
  // Pre-shower and Shower correlation plots
  //-----------------------------------------------------------------------------------------------------------------------------

  if( PlotBBCal ) {

    TCanvas* c3 = new TCanvas("c3","",1200,800);
    c3->Divide(4,2);
    c3->cd(1)->SetLogy(1);
    hbbcal_psE->Draw("");
    if( kin_no != 4 ) 
      hbbcal_psE->GetYaxis()->SetRangeUser(1,6e5);
    hbbcal_psEc->SetLineColor(2);
    hbbcal_psEc->Draw("same");
    hbbcal_psE->GetXaxis()->SetTitle("BBCal PS Energy [GeV]");

    tex = new TLatex( 0.48, 0.8, "All tracks");
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(4);
    tex->SetTextSize(0.055);
    tex->Draw();
    
    tex = new TLatex( 0.48, 0.72, "After cuts");
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(2);
    tex->SetTextSize(0.055);
    tex->Draw();

    c3->cd(2)->SetLogy(1);
    hbbcal_shE->Draw("");
    hbbcal_shEc->SetLineColor(2);
    hbbcal_shEc->Draw("same");
    hbbcal_shE->GetXaxis()->SetTitle("BBCal Shower Energy [GeV]");

    c3->cd(3)->SetLogy(1);
    hbbcal_cale->Draw();
    hbbcal_calec->SetLineColor(2);
    hbbcal_calec->Draw("same");
    hbbcal_cale->GetXaxis()->SetTitle("E_{bbcal} [GeV]");

    c3->cd(4);
    hbbcal2d_pss->Draw("colz");
    hbbcal2d_pss->GetXaxis()->SetTitle("E_{SH}/p_{bbtrack}");
    hbbcal2d_pss->GetYaxis()->SetTitle("E_{PS}/p_{bbtrack}");
    TLine *line1 = new TLine(0,sh_min,sh_min/sh_e,0); 
    line1->SetLineColor(2); 
    line1->SetLineWidth(2); 
    line1->Draw(); 
    TLine *line3 = new TLine(0,sh_max,sh_max/sh_e,0); 
    line3->SetLineColor(2); 
    line3->SetLineWidth(2); 
    line3->Draw(); 
    TLine *line2 = new TLine(0,ps_min,1.2,ps_min); 
    line2->SetLineColor(2); 
    line2->SetLineWidth(2); 
    line2->Draw(); 

    c3->cd(8);
    hbbcal2d_pssc->Draw("colz");
    hbbcal2d_pssc->GetXaxis()->SetTitle("E_{SH}/p_{bbtrack}");
    hbbcal2d_pssc->GetYaxis()->SetTitle("E_{PS}/p_{bbtrack}");

    c3->cd(5)->SetLogy(1);
    hbbcal_edivp->Draw();
    if( kin_no != 4 ) 
      hbbcal_edivp->GetYaxis()->SetRangeUser(1,2e5);
    hbbcal_edivpc->SetLineColor(2);
    hbbcal_edivpc->Draw("same");
    TF1* gausfedp = new TF1("gausfedp","gaus(0)",-1.00, 1.00);
    gausfedp->SetLineColor(1);
    hbbcal_edivpc->Fit(gausfedp,"Q");
    hbbcal_edivp->GetXaxis()->SetTitle("(E_{bbcal} - p_{bbtrack}) [GeV]");

    tex = new TLatex( 0.54, 0.9, Form("#sigma = %3.2f %%", 
				      100.*TMath::Sqrt( (gausfedp->GetParameter(2)/pcent) * 
							(gausfedp->GetParameter(2)/pcent) - pres*pres ) ) );
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(2);
    tex->SetTextSize(0.055);
    tex->Draw();

    c3->cd(6)->SetLogy(1);
    
    TF1* gausf = new TF1("gausf","gaus(0)",-0.20, 0.20);
    hbbcal_xdiff->Draw("");
    gausf->SetLineColor(1);
    hbbcal_xdiffc->SetLineColor(2);
    hbbcal_xdiffc->Draw("same");
    hbbcal_xdiff->GetXaxis()->SetTitle("#delta x_{FP} (GEM-BBSH) [m]");

    TF1* gausfc = new TF1("gausfc","gaus(0)",-0.20, 0.20);
    gausfc->SetLineColor(1);
    hbbcal_xdiffc->Fit(gausfc,"Q");

    tex = new TLatex( 0.54, 0.9, Form("#sigma = %2.1f cm", 100.*gausfc->GetParameter(2)));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(2);
    tex->SetTextSize(0.055);
    tex->Draw();


    c3->cd(7)->SetLogy(1);

    hbbcal_ydiff->Draw("");
    //    hbbcal_xdiff->Fit(gausf,"Q");
    hbbcal_ydiffc->SetLineColor(2);
    hbbcal_ydiffc->Draw("same");
    gausfc->SetLineColor(1);
    hbbcal_ydiffc->Fit(gausfc,"Q");
    hbbcal_ydiff->GetXaxis()->SetTitle("#delta y_{FP} (GEM-BBSH) [m]");

    tex = new TLatex( 0.54, 0.9, Form("#sigma = %2.1f cm", 100.*gausfc->GetParameter(2)));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(2);
    tex->SetTextSize(0.055);
    tex->Draw();

    c3->Print("kinematics_1.pdf");
    c3->Print("bbhodo_temp3.pdf");
    //c3->Close();
  }

  //-----------------------------------------------------------------------------------------------------------------------------
  // Kinematic plots
  //-----------------------------------------------------------------------------------------------------------------------------

  if( PlotKine ) {

    TCanvas* c4 = new TCanvas("c4","",1200,800);
    c4->Divide(4,2);

    c4->cd(1);
    hkin_p->Draw();
    hkin_p->GetXaxis()->SetTitle("p_{bbtrack} [GeV/c]");

    tex = new TLatex( 0.45, 0.8, "All tracks");
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(4);
    tex->SetTextSize(0.05);
    tex->Draw();

    tex = new TLatex( 0.45, 0.75, "(E_{ps} > 0.05 GeV)");
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(4);
    tex->SetTextSize(0.05);
    tex->Draw();

    c4->cd(2);
    hkin_th->Draw();
    hkin_th->GetXaxis()->SetTitle("#theta_tgt_{bbtrack} [rad]");

    c4->cd(3);
    hkin_ph->Draw();
    hkin_ph->GetXaxis()->SetTitle("#phi_tgt_{bbtrack} [rad]");

    c4->cd(4);
    hkin_yt->Draw();
    hkin_yt->GetXaxis()->SetTitle("y_tgt_{bbtrack} [m]");

    c4->cd(5);
    hkin2d_thp->Draw("colz");
    hkin2d_thp->GetXaxis()->SetTitle("#theta_{bblab} [degrees]");
    hkin2d_thp->GetYaxis()->SetTitle("p_{bbtrack} [GeV/c]");

    c4->cd(6);
    hhcal_deltaxy->Draw("colz");
    hhcal_deltaxy->GetXaxis()->SetTitle("HCal (Meas - Pred) y[m]");
    hhcal_deltaxy->GetYaxis()->SetTitle("HCal (Meas - Pred) x[m]");
    
    c4->cd(7);
    hkin_W->Draw();
    hkin_W->GetXaxis()->SetTitle("W^{2} [GeV^{2}]");
    
    c4->cd(8);
    hkin_pdiff->Draw();
    hkin_pdiff->GetXaxis()->SetTitle("p_{bbtrack} - p_{#theta exp} [GeV/c]");

    c4->Print("bbhodo_temp4.pdf");
    c4->Print("kinematics_2.pdf");
    //c4->Close();
    
    TCanvas* c5 = new TCanvas("c5","",1200,800);
    c5->Divide(4,2);

    c5->cd(1);
    hkin_pc->Draw();
    hkin_pc->SetLineColor(2);
    hkin_pc->GetXaxis()->SetTitle("p_{bbtrack} [GeV/c]");

    tex = new TLatex( 0.58, 0.8, "After cuts");
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(2);
    tex->SetTextSize(0.05);
    tex->Draw();

    c5->cd(2);
    hkin_thc->Draw();
    hkin_thc->SetLineColor(2);
    hkin_thc->GetXaxis()->SetTitle("#theta_tgt_{bbtrack} [rad]");

    c5->cd(3);
    hkin_phc->Draw();
    hkin_phc->SetLineColor(2);
    hkin_phc->GetXaxis()->SetTitle("#phi_tgt_{bbtrack} [rad]");

    c5->cd(4);
    hkin_ytc->Draw();
    hkin_ytc->SetLineColor(2);
    hkin_ytc->GetXaxis()->SetTitle("y_tgt_{bbtrack} [m]");

    c5->cd(5);
    hkin2d_thpc->Draw("colz");
    hkin2d_thpc->GetXaxis()->SetTitle("#theta_{bblab} [degrees]");
    hkin2d_thpc->GetYaxis()->SetTitle("p_{bbtrack} [GeV/c]");

    cout << hkin2d_thpc->GetMean(1) << "\t" << hkin2d_thpc->GetMean(2) << endl;
    
    c5->cd(6);
    hkin_deltaxyc->Draw("colz");
    hkin_deltaxyc->GetXaxis()->SetTitle("Hcal (Pred - Meas) y[m]");
    hkin_deltaxyc->GetYaxis()->SetTitle("Hcal (Pred - Meas) x[m]");
    
    c5->cd(7);
    hkin_Wc->Draw();
    hkin_Wc->SetLineColor(2);
    hkin_Wc->GetXaxis()->SetTitle("W^{2} [GeV^{2}]");

    
    c5->cd(8);
    hkin_pdiffc->Draw();
    hkin_pdiffc->SetLineColor(2);
    hkin_pdiffc->GetXaxis()->SetTitle("p_{bbtrack} - p_{#theta exp} [GeV/c]");

    if( DoFit ) {

      Float_t binwidth = hkin_pdiffc->GetXaxis()->GetBinWidth(1); 
      
      //this method works when the peak is clear and higher than the bg
      int maxbin = hkin_pdiffc->GetMaximumBin();
      float a = hkin_pdiffc->GetMaximum();
      float b = hkin_pdiffc->GetBinCenter(maxbin);

      //need this method when the peak is less clear (high Q2 kinematics)
      //float b = 0;
      //float maxbin = hkin_pdiffc->FindBin(0);
      //float a = hkin_pdiffc->GetBinContent(maxbin);
      
      float pmean    = b;
      float pwidth   = 0.06;

      float pmin   = pmean - 3*pwidth;    // peak minimum 
      float pmax   = pmean + 3*pwidth;    // peak maximum 
      
      float cmin   = -0.3;    // fit minimum 
      float cmax   = 0.3;    // fit maximum 

      int pbmin   = hkin_pdiffc->FindBin( pmin ); 
      int pbmax   = hkin_pdiffc->FindBin( pmax ); 
      int pbrange = pbmax - pbmin; 
      float pbins[pbrange], perr[pbrange]; 
      for( int i = pbmin; i < pbmax; i++) { 
	pbins[i - pbmin] = hkin_pdiffc->GetBinContent( i ); 
	perr[i - pbmin]  = hkin_pdiffc->GetBinError( i ); 
	hkin_pdiffc->SetBinContent( i, 0 ); 
	hkin_pdiffc->SetBinError( i, 0 ); 
      } 
      
      TF1* back = new TF1("back", "pol1(0)", cmin, cmax);     // 1-D root function 
      //TF1* back = new TF1("back", "gaus(0)", -0.5, 0.5);     // 1-D root function 
      hkin_pdiffc->Fit("back","Q","",cmin,cmax); 
      
      float abg = back->Eval(b);
      float par0 = back->GetParameter(0); 
      float par1 = back->GetParameter(1); 
      //float par2 = back->GetParameter(2); 
      //float par3 = back->GetParameter(3); 
      
      back->SetLineWidth(2); 
      back->SetLineColor(4); 
      
      
      for( int i = pbmin; i < pbmax; i++) { 
	hkin_pdiffc->SetBinContent( i, pbins[i - pbmin] ); 
	hkin_pdiffc->SetBinError( i, perr[i - pbmin] ); 
      }
      
      
      TF1* peakbg = new TF1("peakbg", "pol1(0)+gaus(4)", cmin, cmax); 
      peakbg->FixParameter( 0, par0   ); 
      peakbg->FixParameter( 1, par1   ); 
      //peakbg->FixParameter( 2, par2   ); 
      //peakbg->FixParameter( 3, par3   ); 
      peakbg->FixParameter( 4, a-abg);
      peakbg->FixParameter( 5, b); 
      peakbg->SetParameter( 6, pwidth );    
      peakbg->SetLineWidth(2); 
      peakbg->SetLineColor(1); 
      hkin_pdiffc->Fit("peakbg","Q","",cmin,cmax); 
      peakbg->Draw("same"); 
      back->Draw("same"); 
      
      float chi2              = peakbg->GetChisquare(); 
      float NDF               = peakbg->GetNDF(); 
      float fitted_mean       = peakbg->GetParameter(5); 
      float fitted_mean_err   = peakbg->GetParError(5); 
      float fitted_sigma      = peakbg->GetParameter(6); 
      float fitted_sigma_err  = peakbg->GetParError(6); 
      
      float integral     = (peakbg->Integral(cmin, cmax) - back->Integral(cmin,cmax)) / binwidth; 
      float integral_err = TMath::Sqrt( integral ); 
      
      if( kin_no == 4)
	tex = new TLatex( 0.43, 0.9, Form("R_{elastic} = %3.1fk /uAh", (0.001*integral/(runtime*avI/3600))));
      else 
	tex = new TLatex( 0.43, 0.9, Form("R_{elastic} = %3.1f /uAh", (integral/(runtime*avI/3600))));
     
      tex->SetNDC(1);
      tex->SetTextFont(42);
      tex->SetTextColor(1);
      tex->SetTextSize(0.05);
      tex->Draw();
      
      tex = new TLatex( 0.43, 0.83, Form("#sigma = %3.2f %%", 100*(fitted_sigma/pcent)));
      tex->SetNDC(1);
      tex->SetTextFont(42);
      tex->SetTextColor(1);
      tex->SetTextSize(0.05);
      tex->Draw();

    }

    c5->Print("bbhodo_temp5.pdf");
    c5->Print("kinematics_3.pdf");
    //c5->Print("kinematics_aftercuts.C");
    //c5->Close();
  }

  //-----------------------------------------------------------------------------------------------------------------------------
  // Optics plots
  //-----------------------------------------------------------------------------------------------------------------------------

  if( PlotOptics ) {

    TCanvas* c6 = new TCanvas("c6","",1200,800);
    c6->Divide(4,2);

    TF1 *f11 = new TF1("f11","[0]+[1]*x",-0.4,0.4);
    f11->SetLineColor(kRed);    


    TF1 *f12 = new TF1("f12","[0]+[1]*x",-0.15,0.15);
    f12->SetLineColor(kRed);    

    c6->cd(1);
    hopt_x_pd->Draw("colz");
    hopt_x_pd->GetXaxis()->SetTitle("x_FP [m]");
    hopt_x_pd->GetYaxis()->SetTitle("p_{bbtrack} - p_{#theta exp} [GeV/c]");
    hopt_x_pd->Fit(f11,"Q");
    
    c6->cd(2);
    hopt_xp_pd->Draw("colz");
    hopt_xp_pd->GetXaxis()->SetTitle("xp_FP [rad]");
    hopt_xp_pd->GetYaxis()->SetTitle("p_{bbtrack} - p_{#theta exp} [GeV/c]");

    c6->cd(3);
    hopt_y_pd->Draw("colz");
    hopt_y_pd->GetXaxis()->SetTitle("y_FP [m]");
    hopt_y_pd->GetYaxis()->SetTitle("p_{bbtrack} - p_{#theta exp} [GeV/c]");
    hopt_y_pd->Fit(f11,"Q");

    c6->cd(4);
    hopt_yp_pd->Draw("colz");
    hopt_yp_pd->GetXaxis()->SetTitle("yp_FP [rad]");
    hopt_yp_pd->GetYaxis()->SetTitle("p_{bbtrack} - p_{#theta exp} [GeV/c]");

    c6->cd(5);
    hopt_p_pd->Draw("colz");
    hopt_p_pd->GetXaxis()->SetTitle("p_{bbtrack} [GeV/c]");
    hopt_p_pd->GetYaxis()->SetTitle("p_{bbtrack} - p_{#theta exp} [GeV/c]");

    c6->cd(6);
    hopt_th_pd->Draw("colz");
    hopt_th_pd->GetXaxis()->SetTitle("#theta_tgt [rad]");
    hopt_th_pd->GetYaxis()->SetTitle("p_{bbtrack} - p_{#theta exp} [GeV/c]");

    c6->cd(7);
    hopt_yt_pd->Draw("colz");
    hopt_yt_pd->GetXaxis()->SetTitle("y_tgt [m]");
    hopt_yt_pd->GetYaxis()->SetTitle("p_{bbtrack} - p_{#theta exp} [GeV/c]");

    c6->cd(8);
    hopt_ph_pd->Draw("colz");
    hopt_ph_pd->GetXaxis()->SetTitle("#phi_tgt [rad]");
    hopt_ph_pd->GetYaxis()->SetTitle("p_{bbtrack} - p_{#theta exp} [GeV/c]");

    c6->Print("bbhodo_temp6.pdf");
    //c6->Print("kinematics_4.pdf");
    c6->Close();
    
    TCanvas* c7 = new TCanvas("c7","",1200,800);
    c7->Divide(4,2);

    c7->cd(1);
    hopt_x_pc->Draw("colz");
    hopt_x_pc->GetXaxis()->SetTitle("x_FP [m]");
    hopt_x_pc->GetYaxis()->SetTitle("p_{bbtrack} [GeV/c]");

    c7->cd(2);
    hopt_xp_pc->Draw("colz");
    hopt_xp_pc->GetXaxis()->SetTitle("xp_FP [rad]");
    hopt_xp_pc->GetYaxis()->SetTitle("p_{bbtrack} [GeV/c]");

    c7->cd(3);
    hopt_y_pc->Draw("colz");
    hopt_y_pc->GetXaxis()->SetTitle("y_FP [m]");
    hopt_y_pc->GetYaxis()->SetTitle("p_{bbtrack} [GeV/c]");

    c7->cd(4);
    hopt_yp_pc->Draw("colz");
    hopt_yp_pc->GetXaxis()->SetTitle("yp_FP [rad]");
    hopt_yp_pc->GetYaxis()->SetTitle("p_{bbtrack} [GeV/c]");

    c7->cd(5);
    hopt_p_pc->Draw("colz");
    hopt_p_pc->GetXaxis()->SetTitle("p_{bbtrack} [GeV/c]");
    hopt_p_pc->GetYaxis()->SetTitle("p_{bbtrack} [GeV/c]");

    c7->cd(6);
    hopt_th_pc->Draw("colz");
    hopt_th_pc->GetXaxis()->SetTitle("#theta_tgt [rad]");
    hopt_th_pc->GetYaxis()->SetTitle("p_{bbtrack} [GeV/c]");

    c7->cd(7);
    hopt_yt_pc->Draw("colz");
    hopt_yt_pc->GetXaxis()->SetTitle("y_tgt [m]");
    hopt_yt_pc->GetYaxis()->SetTitle("p_{bbtrack} [GeV/c]");

    c7->cd(8);
    hopt_ph_pc->Draw("colz");
    hopt_ph_pc->GetXaxis()->SetTitle("#phi_tgt [rad]");
    hopt_ph_pc->GetYaxis()->SetTitle("p_{bbtrack} [GeV/c]");

    c7->Print("bbhodo_temp7.pdf");
    c7->Close();
  }

  
  //-----------------------------------------------------------------------------------------------------------------------------
  // HCal plots
  //-----------------------------------------------------------------------------------------------------------------------------

  if( PlotHCal ) {

    TCanvas* c8 = new TCanvas("c8","",1200,800);
    c8->Divide(4,3);
    c8->cd(1)->SetLogy(1);
    hhcal_x->Draw("");
    hhcal_x->GetXaxis()->SetTitle("HCal x[m]");
    hhcal_xc->SetLineColor(2);
    hhcal_xc->Draw("same");
    hhcal_x->GetYaxis()->SetRangeUser(1,1.2*hhcal_x->GetMaximum());
    
    tex = new TLatex( 0.48, 0.8, "All BB tracks");
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(4);
    tex->SetTextSize(0.055);
    tex->Draw();
    
    tex = new TLatex( 0.48, 0.72, "After elastics cuts");
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(2);
    tex->SetTextSize(0.055);
    tex->Draw();

    c8->cd(2)->SetLogy(1);
    hhcal_y->Draw("");
    hhcal_y->GetXaxis()->SetTitle("HCal y[m]");
    hhcal_yc->SetLineColor(2);
    hhcal_yc->Draw("same");
    hhcal_y->GetYaxis()->SetRangeUser(1,1.2*hhcal_y->GetMaximum());
    
    c8->cd(3);
    hhcal_xy->Draw("colz");
    hhcal_xy->GetXaxis()->SetTitle("HCal y[m]");
    hhcal_xy->GetYaxis()->SetTitle("HCal x[m]");
    
    c8->cd(4);
    hhcal_xyc->Draw("colz");
    hhcal_xyc->GetXaxis()->SetTitle("HCal y[m]");
    hhcal_xyc->GetYaxis()->SetTitle("HCal x[m]");
    
    c8->cd(5)->SetLogy(1);
    hhcal_predx->Draw("");
    hhcal_predx->GetXaxis()->SetTitle("HCal predicted x[m]");
    hhcal_predxc->SetLineColor(2);
    hhcal_predxc->Draw("same");
    
    c8->cd(6)->SetLogy(1);
    hhcal_predy->Draw("");
    hhcal_predy->GetXaxis()->SetTitle("HCal predicted y[m]");
    hhcal_predyc->SetLineColor(2);
    hhcal_predyc->Draw("same");
    
    c8->cd(7);
    hhcal_predxy->Draw("colz");
    hhcal_predxy->GetXaxis()->SetTitle("HCal Pred y[m]");
    hhcal_predxy->GetYaxis()->SetTitle("HCal Pred x[m]");

    c8->cd(8);
    hhcal_predxyc->Draw("colz");
    hhcal_predxyc->GetXaxis()->SetTitle("HCal Pred y[m]");
    hhcal_predxyc->GetYaxis()->SetTitle("HCal Pred x[m]");
   
    c8->cd(9)->SetLogy(1);
    hhcal_deltax->Draw("");
    hhcal_deltax->GetXaxis()->SetTitle("HCal (Meas - Predicted) x[m]");
    hhcal_deltaxc->SetLineColor(2);
    hhcal_deltaxc->Draw("same");
    
    c8->cd(10)->SetLogy(1);
    hhcal_deltay->Draw("");
    hhcal_deltay->GetXaxis()->SetTitle("HCal (Meas - Pred) y[m]");
    hhcal_deltayc->SetLineColor(2);
    hhcal_deltayc->Draw("same");
    
    c8->cd(11);
    hhcal_deltaxy->Draw("colz");
    hhcal_deltaxy->GetXaxis()->SetTitle("HCal (Meas - Pred) y[m]");
    hhcal_deltaxy->GetYaxis()->SetTitle("HCal (Meas - Pred) x[m]");
    
    c8->cd(12);
    hhcal_deltaxyc->Draw("colz");
    hhcal_deltaxyc->GetXaxis()->SetTitle("HCal (Meas - Pred) y[m]");
    hhcal_deltaxyc->GetYaxis()->SetTitle("HCal (Meas - Pred) x[m]");
    
    c8->Print("kinematics_5.pdf");
    c8->Print("bbhodo_temp8.pdf");
    //c8->Close();
  
    
    TCanvas *c9 = new TCanvas("c9","",1200,800);
    c9->Divide(2,1);
    c9->cd(1);
    hhcal_deltax->Draw("");
    hhcal_deltax->GetXaxis()->SetTitle("HCal (Meas - Pred) x[m]");
    hhcal_deltaxc->SetLineColor(2);
    hhcal_deltaxc->Draw("same");
    
    c9->cd(2);
    hhcal_deltaxc->GetXaxis()->SetTitle("HCal (Meas - Pred) x[m]");
    hhcal_deltaxc->SetLineColor(2);
    hhcal_deltaxc->Draw("");
    

    
    //delta_x fit
    float binwidth = hhcal_deltaxc->GetXaxis()->GetBinWidth(1); 
    float binsig = 3.0;
    //proton peak
    int maxbin = hhcal_deltaxc->GetMaximumBin();

    float ap = hhcal_deltaxc->GetMaximum();
    float bp = hhcal_deltaxc->GetBinCenter(maxbin);
    
    float pmean    = bp;
    float pwidth   = 0.15;
    
    float pmin   = pmean - binsig*pwidth;    // peak minimum 
    float pmax   = pmean + binsig*pwidth;    // peak maximum 
    
    float cmin   = -0.95;    // fit minimum 
    float cmax   = 0.2;    // fit maximum 
    //float cmin = -2.5;
    //float cmax = 2.5;
    
    int pbmin   = hhcal_deltaxc->FindBin( pmin ); 
    int pbmax   = hhcal_deltaxc->FindBin( pmax ); 
    
    int pbrange = pbmax - pbmin; 
    float pbins[pbrange], perr[pbrange];
    
    for( int i = pbmin; i < pbmax; i++) { 
      pbins[i - pbmin] = hhcal_deltaxc->GetBinContent( i ); 
      perr[i - pbmin]  = hhcal_deltaxc->GetBinError( i ); 
      hhcal_deltaxc->SetBinContent( i, 0 ); 
      hhcal_deltaxc->SetBinError( i, 0 ); 
    } 
    
    //neutron peak
    int nmaxbin = hhcal_deltaxc->GetMaximumBin();

    float an = hhcal_deltaxc->GetMaximum();
    float bn = hhcal_deltaxc->GetBinCenter(nmaxbin);
    
    float nmean    = bn;
    float nwidth   = 0.15;

    float nmin   = nmean - binsig*nwidth;    // peak minimum 
    float nmax   = nmean + binsig*nwidth;    // peak maximum 
    
    float ncmin = -2.05; // fit minimum
    float ncmax = -1.05; // fit maximum
    
    int nbmin   = hhcal_deltaxc->FindBin( nmin ); 
    int nbmax   = hhcal_deltaxc->FindBin( nmax ); 
    int nbrange = nbmax - nbmin; 
    float nbins[nbrange], nerr[nbrange];
    
    for( int i = nbmin; i < nbmax; i++) { 
      nbins[i - nbmin] = hhcal_deltaxc->GetBinContent( i ); 
      nerr[i - nbmin]  = hhcal_deltaxc->GetBinError( i ); 
      hhcal_deltaxc->SetBinContent( i, 0 ); 
      hhcal_deltaxc->SetBinError( i, 0 ); 
    } 
    
    //fit background without the two peaks
    //TF1* bg = new TF1("bg","gaus(0)",-2.5,2.5);
    TF1* bg = new TF1("bg","pol3(0)",-2.5,2.5);
    bg->SetLineColor(kBlue);
    hhcal_deltaxc->Fit("bg","Q","",-2.5,2.5);
    
    double par0 = bg->GetParameter(0);
    double par1 = bg->GetParameter(1);
    double par2 = bg->GetParameter(2);
    double par3 = bg->GetParameter(3);
    
    TFormula* bgform = bg->GetFormula();
    double apbg = bgform->Eval(bp);
    double anbg = bgform->Eval(bn);
    
    
    //add back in the peaks
    for( int i = pbmin; i < pbmax; i++) { 
      hhcal_deltaxc->SetBinContent( i, pbins[i - pbmin] ); 
      hhcal_deltaxc->SetBinError( i, perr[i - pbmin] ); 
    } 
    for( int i = nbmin; i < nbmax; i++) { 
      hhcal_deltaxc->SetBinContent( i, nbins[i - nbmin] ); 
      hhcal_deltaxc->SetBinError( i, nerr[i - nbmin] ); 
    }
    
    
    //TF1* ppeakbg = new TF1("ppeakbg","gaus(0)+gaus(4)",cmin,cmax);
    //TF1* npeakbg = new TF1("npeakbg","gaus(0)+gaus(4)",ncmin,ncmax);
    TF1* ppeakbg = new TF1("ppeakbg","pol3(0)+gaus(4)",cmin,cmax);
    TF1* npeakbg = new TF1("npeakbg","pol3(0)+gaus(4)",ncmin,ncmax);
    ppeakbg->SetLineColor(kRed);
    npeakbg->SetLineColor(kRed);
    
    ppeakbg->FixParameter(0,par0);
    ppeakbg->FixParameter(1,par1);
    ppeakbg->FixParameter(2,par2);
    ppeakbg->FixParameter(3,par3);
    ppeakbg->FixParameter(4,ap-apbg);
    ppeakbg->FixParameter(5,bp);
    ppeakbg->SetParameter(6,pwidth);
    
    npeakbg->FixParameter(0,par0);
    npeakbg->FixParameter(1,par1);
    npeakbg->FixParameter(2,par2);
    ppeakbg->FixParameter(3,par3);
    npeakbg->FixParameter(4,an-anbg);
    npeakbg->FixParameter(5,bn);
    npeakbg->SetParameter(6,nwidth);
    
    hhcal_deltaxc->Fit("ppeakbg","Q","",cmin,cmax);
    hhcal_deltaxc->Fit("npeakbg","Q","",ncmin,ncmax);
    bg->Draw("same");
    ppeakbg->Draw("same");
    npeakbg->Draw("same");
    
    float pfitted_sigma      = ppeakbg->GetParameter(6); 
    float pfitted_sigma_err  = ppeakbg->GetParError(6); 
    
    float nfitted_sigma      = npeakbg->GetParameter(6); 
    float nfitted_sigma_err  = npeakbg->GetParError(6); 
    
    float pintegral = (ppeakbg->Integral(cmin,cmax) - bg->Integral(cmin,cmax)) / binwidth;
    float nintegral = (npeakbg->Integral(ncmin,ncmax) - bg->Integral(ncmin,ncmax)) / binwidth;

    if( kin_no == 4){
      tex = new TLatex( 0.5, 0.9, Form("R_{proton} = %3.1fk /uAh", (0.001*pintegral/(runtime*avI/3600))));
    tex2 = new TLatex( 0.5, 0.7, Form("R_{neutron} = %3.1fk /uAh", (0.001*nintegral/(runtime*avI/3600))));
    }else{ 
      tex = new TLatex( 0.5, 0.9, Form("R_{proton} = %3.1f /uAh", (pintegral/(runtime*avI/3600))));
      tex2 = new TLatex( 0.5, 0.7, Form("R_{neutron} = %3.1f /uAh", (nintegral/(runtime*avI/3600))));
    }
    
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.05);
    tex->Draw();
      
    tex2->SetNDC(1);
    tex2->SetTextFont(42);
    tex2->SetTextColor(1);
    tex2->SetTextSize(0.05);
    tex2->Draw();
    
    tex = new TLatex( 0.5, 0.83, Form("#sigma = %3.2f m",pfitted_sigma));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.05);
    tex->Draw();

    tex2 = new TLatex( 0.5, 0.63, Form("#sigma = %3.2f m",nfitted_sigma));
    tex2->SetNDC(1);
    tex2->SetTextFont(42);
    tex2->SetTextColor(1);
    tex2->SetTextSize(0.05);
    tex2->Draw();
    
    
    c9->Print("kinematics_6.pdf");
    //c9->Close();
  }

  //gSystem->Exec(Form("pdfunite  bbhodo_temp*.pdf fullanalysis_SBS-%d.pdf", kin_no));  
  gSystem->Exec("rm  bbhodo_temp*.pdf");  

  gSystem->Exec(Form("pdfunite  kinem*.pdf Deep_SBS-%d.pdf", kin_no));  
  gSystem->Exec("rm  kinema*.pdf");  

}

//-----------------------------------------------------------------------------------------------------------------------------
