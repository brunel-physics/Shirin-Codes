//clang++ -Isrc -std=c++17 -march=native -pipe -O3 -Wall -Wextra -Wpedantic -o build/ttreehisto src/ttreehisto.cxx `root-config --libs` -lm

// TODO:: extend cmet to include MET and CMS (passing the same value), it is needed for them to be on the same  stack
// TODO:: Shape uncertainties -> implement lhepdfWeight for up and down based on the formula
// TODO:: change filter to include deltaphi cuts

#include <ROOT/RDataFrame.hxx>//#include <ROOT/RCsvDS.hxx>
#include <Math/Vector4D.h>
#include <TRandom3.h>// used Gaussian, uniform each once
//#include <execution>// need to link -ltbb in Makefile
#include <TChain.h>
#include <TF1.h>
#include "csv.h"
#include "json.hpp"

enum      channel      {elnu,munu};
constexpr channel
          channelAll[]={elnu,munu};

enum      dataSource	  {tzq,zz,tz1,tz2,ww,wz,met,wjt,wjx,st,stb,stw,stbw,wzll,wjqq,zjt1,zjt2,zjt3,zjqq,cms,ttb,ttl,ttj};
constexpr dataSource
          dataSourceAll[]={tzq,zz,tz1,tz2,ww,wz,met,wjt,wjx,st,stb,stw,stbw,wzll,wjqq,zjt1,zjt2,zjt3,zjqq,cms,ttb,ttl,ttj};


enum      PtEtaPhiM	 {pt,eta,phi,m};
constexpr PtEtaPhiM
          PtEtaPhiMall[]={pt,eta,phi,m};

using doubles = ROOT::VecOps::RVec<double>;
using  floats = ROOT::VecOps::RVec<float>;
using    ints = ROOT::VecOps::RVec<int>;
using   bools = ROOT::VecOps::RVec<bool>;
using strings = ROOT::VecOps::RVec<std::string>;

  constexpr  float   MU__PT_MIN   = 29.f;
  constexpr  float   MU_ETA_MAX   = 2.4f;
  constexpr  float   MU_LOOSE_ISO = .25f;
  constexpr  float   MU_TIGHT_ISO = .15f;

//constexpr  float    MET__PT_MIN = 40.f;
  constexpr  float    MET_EL_PT   = 30.f;//80.f;// TODO: Need new values
  constexpr  float    MET_MU_PT   = 45.f;//40.f;

  constexpr double     Z_MASS     =  91.1876;
  constexpr double     Z_MASS_CUT =  20.    ;
  constexpr double     W_MASS     =  80.385 ;
  constexpr double     W_MASS_CUT =  20.    ;
  constexpr double   TOP_MASS     = 172.5   ;
//constexpr double   TOP_MASS_CUT =  20.    ;
  constexpr double DELTA___R_ZL   = 1.6;
  constexpr double DELTA_PHI_ZW   = 2. ;
  constexpr double DELTA_PHI_ZMET = 2. ;


inline auto met_pt_cut(const channel ch){
	float lower_bound;
	switch(ch){
		case elnu:{lower_bound = MET_EL_PT;break;}
		case munu:{lower_bound = MET_MU_PT;break;}
		default  :throw std::invalid_argument(
			"Unimplemented ch (met_pt_cut)");
	}
	return [=](const float   met_pt)->bool
	   {return lower_bound < met_pt;};
}

inline auto easy_mass_cut(const double theo,const double cut){
	return [=](const double ours){return std::abs(ours-theo)<cut;};
}
/*auto deltaphi_cut(const double    x){
      return  [=](const doubles& dps){
		return  std::any(dps >= x);
	};
}*/


void ttreehisto(const channel ch,const dataSource ds){
        ROOT::EnableImplicitMT(4);// SYNC WITH CONDOR JOBS!
	// Open data files even if unused
	// then automatically choose which one to read from
	// No penalty for opening and leaving unused
	// Can even open multiple times at once in parallel
	// Open MC data source EVEN IF UNUSED
std::string chN,dsN,
	    temp_header="BDTInput/Histoffile_",
            temp_opener,temp_footer=".root";

switch(ds){
	case  tzq:{dsN = "tzq";break;}
	case   ww:{dsN =  "ww";break;}
	case   wz:{dsN =  "wz";break;}
	case   zz:{dsN =  "zz";break;}
	case   st:{dsN =  "st";break;}
	case  stb:{dsN = "stb";break;}
	case  stw:{dsN = "stw";break;}
	case stbw:{dsN ="stbw";break;}
	case wjqq:{dsN ="wjqq";break;}
	case wzll:{dsN ="wzll";break;}
	case zjt1:{dsN ="zjt1";break;}
	case zjt2:{dsN ="zjt2";break;}
	case zjt3:{dsN ="zjt3";break;}
	case zjqq:{dsN ="zjqq";break;}
	case  ttb:{dsN = "ttb";break;}
	case  ttl:{dsN = "ttl";break;}
	case  ttj:{dsN = "ttj";break;}
	case  wjt:{dsN = "wjt";break;}
	case  wjx:{dsN = "wjx";break;}
	case  tz1:{dsN = "tz1";break;}
	case  tz2:{dsN = "tz2";break;}
	case  cms:{dsN = "cms";break;}
} // end switch

if(ch == elnu)chN = "elnu";
else{	      chN = "munu";
}// end if
temp_opener= temp_header + chN + dsN + temp_footer;
std::cout<< "chN and dsN are "<<chN<<" "<<dsN<<std::endl;
ROOT::RDataFrame finalDF("Events",temp_opener);// Monte Carlo
std::cout<<"opened file is "<< temp_opener<<std::endl;

//auto finalDF
//	= DF
//	.Filter(met_pt_cut(ch),{"met__pt"},"MET Pt cut")
//	.Filter(easy_mass_cut(W_MASS,W_MASS_CUT),{"tw_lep_mas"},"W mass cut")
//	.Filter( easy_mass_cut(Z_MASS,Z_MASS_CUT),{"z_mas"},"z mass cut")
//	.Filter(      deltaphi_cut(DELTA_PHI_ZW),
//	       {   "zw_Dph"},"delta phi ZW cut")
//	.Filter(      deltaphi_cut(DELTA_PHI_ZMET),
//	       { "zmet_Dph"},"Z met cut ")
	;
// now we make the histogram names and titles
switch(ch){
	case elnu:{temp_header = "elnu_";
		           temp_footer = "electron-neutrino";break;}
	case munu:{temp_header = "munu_";
		           temp_footer = "muon"  "-neutrino";break;}
//		default  :throw std::invalid_argument(
//			"Unimplemented ch (hist titles)");
} // end switch
temp_footer = "pt vs eta in " + temp_footer + " channel for ";
switch(ds){
	case  tzq:{temp_header+= "tzq";temp_footer+="tZq" ;break;}
	case   ww:{temp_header+= "_ww";temp_footer+=" WW" ;break;}
	case   wz:{temp_header+= "_wz";temp_footer+=" WZ" ;break;}
	case   zz:{temp_header+= "_zz";temp_footer+=" ZZ" ;break;}
	case   st:{temp_header+= "_st";temp_footer+=" ST" ;break;}
	case  stb:{temp_header+= "stb";temp_footer+="STB" ;break;}
	case  stw:{temp_header+= "stw";temp_footer+="STW" ;break;}
	case stbw:{temp_header+="stbw";temp_footer+="STBW";break;}
	case wjqq:{temp_header+="wjqq";temp_footer+="wjqq";break;}
	case wzll:{temp_header+="wzll";temp_footer+="WZLL";break;}
	case zjt1:{temp_header+="zjt1";temp_footer+="ZJT1";break;}
	case zjt2:{temp_header+="zjt2";temp_footer+="ZJT2";break;}
	case zjt3:{temp_header+="zjt3";temp_footer+="ZJT3";break;}
	case zjqq:{temp_header+="zjqq";temp_footer+="ZJQQ";break;}
	case  wjt:{temp_header+= "wjt";temp_footer+="Wjt" ;break;}
        case  wjx:{temp_header+= "wjx";temp_footer+="Wjx" ;break;}
	case  ttb:{temp_header+= "ttb";temp_footer+="ttb" ;break;}
        case  ttl:{temp_header+= "ttl";temp_footer+="ttl" ;break;}
        case  ttj:{temp_header+= "ttj";temp_footer+="ttj" ;break;}
        case  tz1:{temp_header+= "tz1";temp_footer+="tz1" ;break;}
        case  tz2:{temp_header+= "tz2";temp_footer+="tz2" ;break;}
	case  cms:{temp_header+= "cms";temp_footer+="CMS" ;break;}
//		default :throw std::invalid_argument(
//			"Unimplemented ds (hist titles)");
} // end switch
	// Histogram names sorted, now branch into MC vs exptData
if(ds != cms){


	const std::vector<double> ptBins
	= {0.,20.,30.,40.,45.,50.,55.,60.,65.,75.,85.,100.,150.,500.};
	const int ptBinsSize = ptBins.size() - 1;


	// Copied to below, skip MC-only, ADD MET_sumEt!
	// Assuming temp_header and footer and all are set per (hist titles)!
// MC only


	auto h_sfj = finalDF.Histo1D({
	("sfj_"+temp_header).c_str(),
	("sfj "+temp_header).c_str(),
	5000,0,1.3},"sfj");
	h_sfj->GetXaxis()->SetTitle("sf_{j}");
	h_sfj->GetYaxis()->SetTitle("Event");
	h_sfj->SetLineStyle(kSolid);
/*
	auto h_p_ei = finalDF.Histo1D({
	("p_ei_"+temp_header).c_str(),
	("p_ei "+temp_header).c_str(),
	5000,0,1.3},"P___ei");
	h_p_ei->GetXaxis()->SetTitle("\\prod_{i} e_{i}");
	h_p_ei->GetYaxis()->SetTitle("Event");
	h_p_ei->SetLineStyle(kSolid);
*/
	auto h_p_ej = finalDF.Histo1D({
	("p_ej_"+temp_header).c_str(),
	("p_ej "+temp_header).c_str(),
	5000,0,1.3},"P___ej");
	h_p_ej->GetXaxis()->SetTitle("\\prod_{j} 1 - e_{j}");
	h_p_ej->GetYaxis()->SetTitle("Event");
	h_p_ej->SetLineStyle(kSolid);

	auto h_p_sf_i = finalDF.Histo1D({
	("p_sf_i_"+temp_header).c_str(),
	("p_sf_i "+temp_header).c_str(),
	5000,0,1.3},"P_sf_i");
	h_p_sf_i->GetXaxis()->SetTitle("\\prod_{i} \\text{sf}_{i}");// e_{i}");
	h_p_sf_i->GetYaxis()->SetTitle("Event");
	h_p_sf_i->SetLineStyle(kSolid);

	auto h_p_sfej = finalDF.Histo1D({
	("p_sfej_"+temp_header).c_str(),
	("p_sfej "+temp_header).c_str(),
	5000,0,1.3},"P_sfej");
	h_p_sfej->GetXaxis()->SetTitle("\\prod_{j} 1 - \\text{sf}_{j} e_{j}");
	h_p_sfej->GetYaxis()->SetTitle("Event");
	h_p_sfej->SetLineStyle(kSolid);

	auto h_btag_w = finalDF.Histo1D({
	("btag_w_"+temp_header).c_str(),
	("btag_W "+temp_header).c_str(),
	5000,-1,1.3},"btag_w");
	h_btag_w->GetXaxis()->SetTitle("btag weight");
	h_btag_w->GetYaxis()->SetTitle("Event");
	h_btag_w->SetLineStyle(kSolid);

        auto h_mostSF = finalDF.Histo1D({
        ("mostSF_"+temp_header).c_str(),
        ("mostSF "+temp_header).c_str(),
        5000,-1,1.3},"btag_w");
        h_mostSF->GetXaxis()->SetTitle("most SF");
        h_mostSF->GetYaxis()->SetTitle("Event");
        h_mostSF->SetLineStyle(kSolid);

        auto h_ttbSF = finalDF.Histo1D({
        ("ttbSF_"+temp_header).c_str(),
        ("ttbSF "+temp_header).c_str(),
        5000,-1,1.3},"ttbSF");
        h_btag_w->GetXaxis()->SetTitle("ttb SF");
        h_btag_w->GetYaxis()->SetTitle("Event");
        h_btag_w->SetLineStyle(kSolid);



	auto h_cmet_sEt = finalDF.Histo1D({
	("cmet_sEt_"+temp_header).c_str(),
	("cmet_sEt "+temp_header).c_str(),
	100,0,600},"cmet_sEt","nw_cmet_sEt");
	h_cmet_sEt->GetXaxis()->SetTitle("corrected MET Sum Et (GeV)");
	h_cmet_sEt->GetYaxis()->SetTitle("Event");
	h_cmet_sEt->SetLineStyle(kSolid);

	auto h_cmet__pt = finalDF.Histo1D({
	("cmet__pt_"+temp_header).c_str(),
	("cmet__pt "+temp_header).c_str(),
	50,0,300},"cmet__pt","nw_cmet__pt");
	h_cmet__pt->GetXaxis()->SetTitle("corrected MET p_{T} (GeV/c)");
	h_cmet__pt->GetYaxis()->SetTitle("Event");
	h_cmet__pt->SetLineStyle(kSolid);

        auto h_cmet_phi = finalDF.Histo1D({
        ("cmet_phi_"+temp_header).c_str(),
        ("cmet_phi "+temp_header).c_str(),
        50,0,8},"cmet_phi","nw_cmet_phi");
        h_cmet_phi->GetXaxis()->SetTitle("corrected MET #phi /rad");
	h_cmet_phi->GetYaxis()->SetTitle("Event");
        h_cmet_phi->SetLineStyle(kSolid);

	auto h_cmet_dpx = finalDF.Histo1D({
	("cmet_dpx_"+temp_header).c_str(),
	("cmet_dpx "+temp_header).c_str(),
	100,-300,300},"cmet_dpx");
	h_cmet_dpx->GetXaxis()->SetTitle("MET correction p_{x} (GeV/c)");
	h_cmet_dpx->GetYaxis()->SetTitle("Event");
	h_cmet_dpx->SetLineStyle(kSolid);

	auto h_cmet_dpy = finalDF.Histo1D({
	("cmet_dpy_"+temp_header).c_str(),
	("cmet_dpy "+temp_header).c_str(),
	100,-300,300},"cmet_dpy");
	h_cmet_dpy->GetXaxis()->SetTitle("MET correction p_{y} (GeV/c)");
	h_cmet_dpy->GetYaxis()->SetTitle("Event");
	h_cmet_dpy->SetLineStyle(kSolid);
// end MC only

	auto h_met_sEt = finalDF.Histo1D({
	("met_sEt_"+temp_header).c_str(),
	("met_sEt "+temp_header).c_str(),
	100,0,600},"MET_sumEt");
	h_met_sEt->GetXaxis()->SetTitle("MET Sum Et (GeV)");
	h_met_sEt->GetYaxis()->SetTitle("Event");
	h_met_sEt->SetLineStyle(kSolid);

	auto h_met__pt = finalDF.Histo1D({
	("met__pt_"+temp_header).c_str(),
	("met__pt "+temp_header).c_str(),
	50,0,300},"MET_pt");
	h_met__pt->GetXaxis()->SetTitle("MET p_{T} (GeV/c)");
	h_met__pt->GetYaxis()->SetTitle("Event");
	h_met__pt->SetLineStyle(kSolid);

	auto h_trans_T = finalDF.Histo1D({
	(          "tTm_"     + temp_header).c_str(),
	("Transverse T mass " + temp_header).c_str(),
	50,0,250},
	"ttop_mas","nw_ttop_mas");
	//h_trans_T->Fit("Gaus"); Didn't work!
	TF1 *f1 = new TF1("f1","gaus",0,250);
	f1->SetParameters(h_trans_T->GetMaximum()
			, h_trans_T->GetMean()
			, h_trans_T->GetRMS() );
	h_trans_T->Fit("f1");
	h_trans_T->GetXaxis()->SetTitle("\\text{mass GeV/}c^{2}");
	h_trans_T->GetYaxis()->SetTitle("Event");
	h_trans_T->SetLineStyle(kSolid);

	auto h_ttop_pt = finalDF.Histo1D({
        ("ttop_pt_"+temp_header).c_str(),
        ("Transverse top pt "+temp_header).c_str(),
        50,0,300},
	"ttop__pt","nw_ttop__pt");
        h_ttop_pt->GetXaxis()->SetTitle("transverse top p_{T} (GeV/c)");
        h_ttop_pt->GetYaxis()->SetTitle("Event");
        h_ttop_pt->SetLineStyle(kSolid);

	auto h_trans_w = finalDF.Histo1D({
	(          "tWm_"     + temp_header).c_str(),
	("Transverse W mass " + temp_header).c_str(),
	50,0,160},
	"tw_lep_mas","nw_tw_lep_mas");
	//h_trans_w->Fit("Gaus"); Didn't work!
	TF1 *f2 = new TF1("f2","gaus",0,160);
	f2->SetParameters(h_trans_w->GetMaximum()
	                 ,h_trans_w->GetMean()
	                 ,h_trans_w->GetRMS ()  );
	h_trans_w->Fit("f2");
	h_trans_w->GetXaxis()->SetTitle("\\text{mass GeV/}c^{2}");
	h_trans_w->GetYaxis()->SetTitle("Event");
	h_trans_w->SetLineStyle(kSolid);

	auto h_Winvmas = finalDF.Histo1D({
	("W_invariant_mass_" + temp_header).c_str(),
	("W invariant mass " + temp_header).c_str(),
	50,0,160},
	"lep_nu_invmass","nw_lep_nu_invmass");
	h_Winvmas->GetXaxis()->SetTitle("\\text{mass GeV/}c^{2}");
	h_Winvmas->GetYaxis()->SetTitle("Event");
	h_Winvmas->SetLineStyle(kSolid);

	auto h_ev_w = finalDF.Histo1D({
	(   "ev_w_"     +temp_header).c_str(),
	("Event weight "+temp_header).c_str(),
	50,-1,3},"sf");
	h_ev_w->GetXaxis()->SetTitle("weight");
	h_ev_w->GetYaxis()->SetTitle("Event" );
	h_ev_w->SetLineStyle(kSolid);

	auto h_z_mas = finalDF.Histo1D({
	(        "zmas_"  + temp_header).c_str(),
	("Recon. Z mass " + temp_header).c_str(),
	50,0,150},
	"z_mas","nw_z_mas");
	h_z_mas->GetXaxis()->SetTitle("\\text{mass GeV/}c^{2}");
	h_z_mas->GetYaxis()->SetTitle("Event");
	h_z_mas->SetLineStyle(kSolid);

	auto h_zw_Dph = finalDF.Histo1D({
	("Z_W_Delta_Phi_" + temp_header).c_str(),
	("Z W Delta Phi " + temp_header).c_str(),
	50,0,3.2},
	"zw_Dph","nw_zw_Dph");
	h_zw_Dph->GetXaxis()->SetTitle("Z & W #Delta#phi/rad");
	h_zw_Dph->GetYaxis()->SetTitle("Event");
	h_zw_Dph->SetLineStyle(kSolid);

	auto h_zmet_Dph = finalDF.Histo1D({
	("Z_MET_Delta_Phi_" + temp_header).c_str(),
	("Z MET Delta Phi " + temp_header).c_str(),
	50,0,3.2},
	"zmet_Dph","nw_zmet_Dph");
	h_zmet_Dph->GetXaxis()->SetTitle("Z & MET #Delta#phi/rad");
	h_zmet_Dph->GetYaxis()->SetTitle("Event");
	h_zmet_Dph->SetLineStyle(kSolid);

        auto h_WZ_deltaR = finalDF.Histo1D({
        ("WZ_DeltaR_"  + temp_header).c_str(),
        ("W Z deltaR"  + temp_header).c_str(),
        50,0,7},
        "WZ_deltaR","nw_WZ_deltaR");
        h_WZ_deltaR->GetXaxis()->SetTitle("Z & W #DeltaR");
        h_WZ_deltaR->GetYaxis()->SetTitle("Event");
        h_WZ_deltaR->SetLineStyle(kSolid);


	auto h_z_daughters_Dph = finalDF.Histo1D({
	("Z_pair_jets_Delta_Phi_" + temp_header).c_str(),
	("Z pair jets Delta Phi " + temp_header).c_str(),
	50,0,7},
	"z_pair_Dph","nw_z_pair_Dph");
	h_z_daughters_Dph->GetXaxis()->SetTitle("Z pair jets #Delta#phi/rad");
	h_z_daughters_Dph->GetYaxis()->SetTitle("Event");
	h_z_daughters_Dph->SetLineStyle(kSolid);

/*	auto h_npl = finalDF.Histo1D({
	("npl_" + temp_header).c_str(),
	("npl " + temp_header).c_str(),
	500,-0.1,0.1},
	"npl_est");
	h_npl->GetXaxis()->SetTitle("NPL");
	h_npl->GetYaxis()->SetTitle("Event");
	h_npl->SetLineStyle(kSolid);
*/
	auto h_tWmVsZmass = finalDF.Histo2D({
	("tWmVsZmass_" + temp_header).c_str(),
	("tWmVsZmass " + temp_header).c_str(),
	50,0,180,50,0,150},
	"tw_lep_mas","z_mas");
	h_tWmVsZmass->GetXaxis()->SetTitle("\\text{  tWm  GeV/}c^{2}");
	h_tWmVsZmass->GetYaxis()->SetTitle("\\text{Z mass GeV/}c^{2}");

	// cut histograms


	// write histograms to a root file
	// ASSUMES temp_header is correct!
	TFile hf(("histo/BDT_"+temp_header+".root").c_str(),"RECREATE");
// MC only
//		hf.WriteTObject(h_sfi                  .GetPtr());hf.Flush();sync();
		hf.WriteTObject(h_sfj                  .GetPtr());hf.Flush();sync();
//		hf.WriteTObject(h_p_ei                 .GetPtr());hf.Flush();sync();
		hf.WriteTObject(h_p_ej                 .GetPtr());hf.Flush();sync();
		hf.WriteTObject(h_p_sf_i               .GetPtr());hf.Flush();sync();
		hf.WriteTObject(h_p_sfej               .GetPtr());hf.Flush();sync();
		hf.WriteTObject(h_btag_w               .GetPtr());hf.Flush();sync();
                hf.WriteTObject(h_mostSF               .GetPtr());hf.Flush();sync();
                hf.WriteTObject(h_ttbSF                .GetPtr());hf.Flush();sync();
		hf.WriteTObject(h_cmet_sEt             .GetPtr());hf.Flush();sync();
		hf.WriteTObject(h_cmet__pt             .GetPtr());hf.Flush();sync();
                hf.WriteTObject(h_cmet_phi             .GetPtr());hf.Flush();sync();
		hf.WriteTObject(h_cmet_dpx             .GetPtr());hf.Flush();sync();
		hf.WriteTObject(h_cmet_dpy             .GetPtr());hf.Flush();sync();

		//hf.WriteTObject(h_is_btag_numer_PtVsEta.GetPtr());hf.Flush();sync();
		//hf.WriteTObject(h_no_btag_numer_PtVsEta.GetPtr());hf.Flush();sync();
		//hf.WriteTObject(h_is_btag_denom_PtVsEta.GetPtr());hf.Flush();sync();
		//hf.WriteTObject(h_no_btag_denom_PtVsEta.GetPtr());hf.Flush();sync();
		//hf.WriteTObject(is_btag_ratio)                   ;hf.Flush();sync();
		//hf.WriteTObject(no_btag_ratio)                   ;hf.Flush();sync();
// end MC only;
	hf.WriteTObject(h_met_sEt        .GetPtr());hf.Flush();sync();
	hf.WriteTObject(h_met__pt        .GetPtr());hf.Flush();sync();
	hf.WriteTObject(h_ttop_pt        .GetPtr());hf.Flush();sync();
	hf.WriteTObject(h_trans_T        .GetPtr());hf.Flush();sync();
	hf.WriteTObject(h_trans_w        .GetPtr());hf.Flush();sync();
	hf.WriteTObject(h_Winvmas        .GetPtr());hf.Flush();sync();

	hf.WriteTObject(h_ev_w           .GetPtr());hf.Flush();sync();
	hf.WriteTObject(h_z_mas          .GetPtr());hf.Flush();sync();
	hf.WriteTObject(  h_zw_Dph       .GetPtr());hf.Flush();sync();
	hf.WriteTObject(h_zmet_Dph       .GetPtr());hf.Flush();sync();
        hf.WriteTObject(h_WZ_deltaR      .GetPtr());hf.Flush();sync();
	hf.WriteTObject(h_z_daughters_Dph.GetPtr());hf.Flush();sync();

	//hf.WriteTObject(h_npl            .GetPtr());hf.Flush();sync();
	hf.WriteTObject(h_tWmVsZmass     .GetPtr());hf.Flush();sync();
        //hf.WriteTObject(h_ttop_mas_w     .GetPtr());hf.Flush();sync();
        //hf.WriteTObject(h_ttop_mas_z     .GetPtr());hf.Flush();sync();
        //hf.WriteTObject(h_ttop_mas_wz    .GetPtr());hf.Flush();sync();
        //hf.WriteTObject(h_ttop_pt_w      .GetPtr());hf.Flush();sync();
        //hf.WriteTObject(h_ttop_pt_z      .GetPtr());hf.Flush();sync();
        //hf.WriteTObject(h_ttop_pt_wz     .GetPtr());hf.Flush();sync();

	// the following two for loops stack correctly
	for(std::string particle:{"fin_jets","lep","bjet"})
	for(PtEtaPhiM k:PtEtaPhiMall){
//		if( e == k ) continue;
		std::string  kstring  = "_";
		std::string  xAxisStr = " ";
		double xmin,xmax;
		switch(k){
			case pt :{kstring += "_pt";xmin =  0;xmax = 200;
			          xAxisStr = "p_{T} (GeV/c)"           ;break;}
			case eta:{kstring += "eta";xmin = -3;xmax =  3 ;
			          xAxisStr = "PseudoRapidity  #eta"    ;break;}
			case phi:{kstring += "phi";xmin = -7;xmax =  7 ;
			          xAxisStr = "Azimuthal angle #phi/rad";break;}
			case  m :{kstring += "mas";xmin =  0;xmax = 200;
			          xAxisStr = "\\text{mass GeV/}c^{2}"  ;break;}
//			case  e :
			default :throw std::invalid_argument(
				"Unimplemented component (histo)");
		}
		temp_footer = particle    + kstring;
		temp_opener = temp_header + "_" + temp_footer;
		auto
		h = finalDF.Histo1D(
		{temp_opener.c_str()
		,temp_opener.c_str()
		,50,xmin,xmax}
		,       temp_footer .c_str()
		,("nw_"+temp_footer).c_str()
		);
		h->SetLineStyle(kSolid);
		h->GetXaxis()->SetTitle(xAxisStr.c_str());
		h->GetYaxis()->SetTitle("Event");
		hf.WriteTObject(h.GetPtr());hf.Flush();sync();
}}else{


	auto h_met_sEt = finalDF.Histo1D({
	("met_sEt_"+temp_header).c_str(),
	("met_sEt "+temp_header).c_str(),
	100,0,600},"MET_sumEt");
	h_met_sEt->GetXaxis()->SetTitle("MET Sum Et (GeV)");
	h_met_sEt->GetYaxis()->SetTitle("Event");
	h_met_sEt->SetLineStyle(kSolid);

	auto h_met__pt = finalDF.Histo1D({
	("met__pt_"+temp_header).c_str(),
	("met__pt "+temp_header).c_str(),
	50,0,300},"MET_pt");
	h_met__pt->GetXaxis()->SetTitle("MET p_{T} (GeV/c)");
	h_met__pt->GetYaxis()->SetTitle("Event");
	h_met__pt->SetLineStyle(kSolid);

	auto h_trans_T = finalDF.Histo1D({
	(          "tTm_"     + temp_header).c_str(),
	("Transverse T mass " + temp_header).c_str(),
	50,0,250},
	"ttop_mas","nw_ttop_mas");
	h_trans_T->GetXaxis()->SetTitle("\\text{mass GeV/}c^{2}");
	h_trans_T->GetYaxis()->SetTitle("Event");
	h_trans_T->SetLineStyle(kSolid);

	auto h_trans_w = finalDF.Histo1D({
	(          "tWm_"     + temp_header).c_str(),
	("Transverse W mass " + temp_header).c_str(),
	50,0,160},
	"tw_lep_mas","nw_tw_lep_mas");
	h_trans_w->GetXaxis()->SetTitle("\\text{mass GeV/}c^{2}");
	h_trans_w->GetYaxis()->SetTitle("Event");
	h_trans_w->SetLineStyle(kSolid);

	auto h_Winvmas = finalDF.Histo1D({
	("W_invariant_mass_" + temp_header).c_str(),
	("W invariant mass " + temp_header).c_str(),
	50,0,160},
	"lep_nu_invmass","nw_lep_nu_invmass");
	h_Winvmas->GetXaxis()->SetTitle("\\text{mass GeV/}c^{2}");
	h_Winvmas->GetYaxis()->SetTitle("Event");
	h_Winvmas->SetLineStyle(kSolid);

	auto h_ev_w = finalDF.Histo1D({
	(   "ev_w_"     +temp_header).c_str(),
	("Event weight "+temp_header).c_str(),
	50,-1,3},"sf");
	h_ev_w->GetXaxis()->SetTitle("weight");
	h_ev_w->GetYaxis()->SetTitle("Event" );
	h_ev_w->SetLineStyle(kSolid);

	auto h_z_mas = finalDF.Histo1D({
	(        "zmas_"  + temp_header).c_str(),
	("Recon. Z mass " + temp_header).c_str(),
	50,0,150},
	"z_mas","nw_z_mas");
	h_z_mas->GetXaxis()->SetTitle("\\text{mass GeV/}c^{2}");
	h_z_mas->GetYaxis()->SetTitle("Event");
	h_z_mas->SetLineStyle(kSolid);

	auto h_zw_Dph = finalDF.Histo1D({
	("Z_W_Delta_Phi_" + temp_header).c_str(),
	("Z W Delta Phi " + temp_header).c_str(),
	50,0,3.2},
	"zw_Dph","nw_zw_Dph");
	h_zw_Dph->GetXaxis()->SetTitle("Z & W #Delta#phi/rad");
	h_zw_Dph->GetYaxis()->SetTitle("Event");
	h_zw_Dph->SetLineStyle(kSolid);

	auto h_zmet_Dph = finalDF.Histo1D({
	("Z_MET_Delta_Phi_" + temp_header).c_str(),
	("Z MET Delta Phi " + temp_header).c_str(),
	50,0,3.2},
	"zmet_Dph","nw_zmet_Dph");
	h_zmet_Dph->GetXaxis()->SetTitle("Z & MET #Delta#phi/rad");
	h_zmet_Dph->GetYaxis()->SetTitle("Event");
	h_zmet_Dph->SetLineStyle(kSolid);

	auto h_z_daughters_Dph = finalDF.Histo1D({
	("Z_pair_jets_Delta_Phi_" + temp_header).c_str(),
	("Z pair jets Delta Phi " + temp_header).c_str(),
	50,0,7},
	"z_pair_Dph","nw_z_pair_Dph");
	h_z_daughters_Dph->GetXaxis()->SetTitle("Z pair jets #Delta#phi/rad");
	h_z_daughters_Dph->GetYaxis()->SetTitle("Event");
	h_z_daughters_Dph->SetLineStyle(kSolid);

	auto h_tWmVsZmass = finalDF.Histo2D({
	("tWmVsZmass_" + temp_header).c_str(),
	("tWmVsZmass " + temp_header).c_str(),
	50,0,160,50,0,160},
	"tw_lep_mas","z_mas");
	h_tWmVsZmass->GetXaxis()->SetTitle("\\text{  tWm  GeV/}c^{2}");
	h_tWmVsZmass->GetYaxis()->SetTitle("\\text{Z mass GeV/}c^{2}");
	// No SetLineStyle here

	//Reporting the filter
	finalDF.Report() ->Print();
	// write histograms to a root file
	// ASSUMES temp_header is correct!
	TFile hf(("histo/BDT_"+temp_header+".root").c_str(),"RECREATE");
	hf.WriteTObject(h_met_sEt        .GetPtr());hf.Flush();sync();
	hf.WriteTObject(h_met__pt        .GetPtr());hf.Flush();sync();
	hf.WriteTObject(h_trans_T        .GetPtr());hf.Flush();sync();
	hf.WriteTObject(h_trans_w        .GetPtr());hf.Flush();sync();
	hf.WriteTObject(h_Winvmas        .GetPtr());hf.Flush();sync();

	hf.WriteTObject(h_ev_w           .GetPtr());hf.Flush();sync();

	hf.WriteTObject(h_z_mas          .GetPtr());hf.Flush();sync();
	hf.WriteTObject(  h_zw_Dph       .GetPtr());hf.Flush();sync();
	hf.WriteTObject(h_zmet_Dph       .GetPtr());hf.Flush();sync();
	hf.WriteTObject(h_z_daughters_Dph.GetPtr());hf.Flush();sync();
	hf.WriteTObject(h_tWmVsZmass     .GetPtr());hf.Flush();sync();

	// the following two for loops stack correctly
	for(std::string particle:{"fin_jets","lep","bjet"})
	for(PtEtaPhiM k:PtEtaPhiMall){
//		if( e == k ) continue;
		std::string  kstring  = "_";
		std::string  xAxisStr = " ";
		double xmin,xmax;
		switch(k){
			case pt :{kstring += "_pt";xmin =  0;xmax = 200;
			          xAxisStr = "p_{T} (GeV/c)"           ;break;}
			case eta:{kstring += "eta";xmin = -3;xmax =  3 ;
			          xAxisStr = "PseudoRapidity  #eta"    ;break;}
			case phi:{kstring += "phi";xmin = -7;xmax =  7 ;
			          xAxisStr = "Azimuthal angle #phi/rad";break;}
			case  m :{kstring += "mas";xmin =  0;xmax = 200;
			          xAxisStr = "\\text{mass GeV/}c^{2}"  ;break;}
//			case  e :
			default :throw std::invalid_argument(
				"Unimplemented component (histo)");
		}
		temp_footer = particle    + kstring;
		temp_opener = temp_header + "_" + temp_footer;
		auto
		h = finalDF.Histo1D(
		{temp_opener.c_str()
		,temp_opener.c_str()
		,50,xmin,xmax}
		,       temp_footer .c_str()
		,("nw_"+temp_footer).c_str()
		);
		h->SetLineStyle(kSolid);
		h->GetXaxis()->SetTitle(xAxisStr.c_str());
		h->GetYaxis()->SetTitle("Event");
		hf.WriteTObject(h.GetPtr());hf.Flush();sync();
	}}
	std::cout<<"ttreehisto successfully completed"<<std::endl;
}
int main ( int argc , char *argv[] ){
        if ( argc < 2 ) {
                std::cout << "Error: no command provided" << std::endl ;
                return 1 ;
        }
        if ( argc < 3 ) {
                   std::cout
                << "Error: NPL_run needs channel and data source"
                << std::endl
                << "e.g.   NPL_run elnu DY"
                << std::endl
                ;
                return 2 ;
        }
	channel c ; dataSource d ;
                if ( const auto chN = std::string_view( argv[1] ) ; false ) ;

        else if ( "elnu"  == chN ) c = elnu ;
        else if ( "munu"  == chN ) c = munu ;
        else { std::cout << "Error: channel " << chN
                << " not recognised" << std::endl ;
                return 3 ;
        }
             if ( const auto dsN = std::string_view( argv[2] ) ; false ) ;
        else if ( "tzq"  ==  dsN ){ d = tzq   ;}
        else if ( "ttb"  ==  dsN ){ d = ttb   ;}
        else if ( "ttl"  ==  dsN ){ d = ttl   ;}
        else if ( "ttj"  ==  dsN ){ d = ttj   ;}
        else if ( "tz1"  ==  dsN ){ d = tz1   ;}
        else if ( "tz2"  ==  dsN ){ d = tz2   ;}
        else if ( "wjt"  ==  dsN ){ d = wjt   ;}
        else if ( "wjx"  ==  dsN ){ d = wjx   ;}
        else if ( "ww"   ==  dsN ){ d = ww    ;}
        else if ( "wz"   ==  dsN ){ d = wz    ;}
        else if ( "zz"   ==  dsN ){ d = zz    ;}
        else if ( "st"   ==  dsN ){ d = st    ;}
        else if ( "stb"  ==  dsN ){ d = stb   ;}
        else if ( "stw"  ==  dsN ){ d = stw   ;}
        else if ( "stbw" ==  dsN ){ d = stbw  ;}
        else if ( "wzll" ==  dsN ){ d = wzll  ;}
        else if ( "wjqq" ==  dsN ){ d = wjqq  ;}
        else if ( "zjt1" ==  dsN ){ d = zjt1  ;}
        else if ( "zjt2" ==  dsN ){ d = zjt2  ;}
        else if ( "zjt3" ==  dsN ){ d = zjt3  ;}
        else if ( "zjqq" ==  dsN ){ d = zjqq  ;}
        else if ( "cms"  ==  dsN ){ d = cms   ;}
        else { std::cout << "Error: data source " << dsN
                << " not recognised" << std::endl ;
                return 4 ;
        }
                ttreehisto(c,d) ;
                return 0 ;
}

