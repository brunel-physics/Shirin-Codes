//clang++ -Isrc -std=c++17 -march=native -pipe -O3 -Wall -Wextra -Wpedantic -o build/addhists src/AddhistsCalc.cxx `root-config --libs` -lm
#include <ROOT/RDataFrame.hxx>//#include <ROOT/RCsvDS.hxx>

#include "calchisto.hpp"

int debug = 1;

namespace{
	std::string
	NPLc = "elnu",
	NPLds = "tzq",
	temp_header = "histo/NPL_run_",
	temp_footer = ".root",
	temp_opener
	;
	TH2D *hbtagw, *hlp_sf, *htpt_w,
 	     *hcmtet, *hcmtpt, *hcmtph,
	     *hmteta, *hmt_pt, /**hmetph*/,
	     *ht__pt, *ht_mas, *hWinvm,
	     *hev_sf, *hz_mas, *hzwdph,
	     *hzmdph, *hwz_dr, *hzjdph,
	     *h__npl, *hjt_pt, *hjteta,
	     *hjtphi, *hjtmas, *hlp_pt,
	     *hlpeta, *hlpphi, *hlpmas,
	     *hb__pt, *hb_eta, *hb_phi,
	     *hb_mas;// for taking from files
	TH2D *fbtagw, *flp_sf, *ftpt_w,
             *fcmtet, *fcmtpt, *fcmtph,
             *fmteta, *fmt_pt, /**hmetph*/,
             *ft__pt, *ft_mas, *fWinvm,
             *fev_sf, *fz_mas, *fzwdph,
             *fzmdph, *fwz_dr, *fzjdph,
             *f__npl, *fjt_pt, *fjteta,
             *fjtphi, *fjtmas, *flp_pt,
             *flpeta, *flpphi, *flpmas,
             *fb__pt, *fb_eta, *fb_phi,
             *fb_mas;      // finals go here
}

void addhistsCalc(const channel ch){
	// do tzq first so that tpr & tln are not null
	if(munu == ch) NPLc = "munu";
		temp_opener = temp_header + NPLc + "_" + NPLds + temp_footer;
		std::cout << "Opening file " << temp_opener << std::endl;
		TFile zq(temp_opener.c_str());
		if( ! zq.IsOpen()) throw std::runtime_error("File is not opened");

		zq ->GetObject(("btag_w_" + NPLc + "_" + NPLds).c_str(),hbtagw);
		hbtagw->SetDirectory(nullptr);// make it stay even if file closed

		zq ->GetObject(("mostSF_" + NPLc + "_" + NPLds).c_str(),hlp_sf);
		hlp_sf->SetDirectory(nullptr);// make it stay even if file closed

		zq ->GetObject(( "ttbSF_" + NPLc + "_" + NPLds).c_str(),htpt_w);
		htpt_w->SetDirectory(nullptr);// make it stay even if file closed

		zq ->GetObject(("cmet_sEt_" + NPLc + "_" + NPLds).c_str(),hcmtet);
		hcmtet->SetDirectory(nullptr);// make it stay even if file close

		zq ->GetObject(("cmet__pt_" + NPLc + "_" + NPLds).c_str(),hcmtpt);
		hcmtpt->SetDirectory(nullptr);// make it stay even if file closed

		zq ->GetObject(("cmet_phi_" + NPLc + "_" + NPLds).c_str(),hcmtph);
		hcmtph->SetDirectory(nullptr);// make it stay even if file closed

		zq ->GetObject(("met_sEt_" + NPLc + "_" + NPLds).c_str(),hmteta);
		hmteta->SetDirectory(nullptr);// make it stay even if file closed

		zq ->GetObject(("met__pt_" + NPLc + "_" + NPLds).c_str(),hmt_pt);
		hmt_pt->SetDirectory(nullptr);// make it stay even if file closed

		zq ->GetObject(("ttop_pt_" + NPLc + "_" + NPLds).c_str(),ht__pt);
		ht__pt->SetDirectory(nullptr);// make it stay even if file closed

		zq ->GetObject(("tTm_" + NPLc + "_" + NPLds).c_str(),ht_mas);
		ht_mas->SetDirectory(nullptr);// make it stay even if file closed

		zq ->GetObject(("W_invariant_mass_" + NPLc + "_" + NPLds).c_str(),hWinvm);
		hWinvm->SetDirectory(nullptr);// make it stay even if file closed

		zq ->GetObject(("ev_w_" + NPLc + "_" + NPLds).c_str(),hev_sf);
		hev_sf->SetDirectory(nullptr);// make it stay even if file closed

		zq ->GetObject(("zmas_" + NPLc + "_" + NPLds).c_str(),hz_mas);
		hz_mas->SetDirectory(nullptr);// make it stay even if file closed

		zq ->GetObject(("Z_W_Delta_Phi_" + NPLc + "_" + NPLds).c_str(),hzwdph);
		hzwdph->SetDirectory(nullptr);// make it stay even if file closed

		zq ->GetObject(("Z_MET_Delta_Phi_" + NPLc + "_" + NPLds).c_str(),hzmdph);
		hzmdph->SetDirectory(nullptr);// make it stay even if file closed

		zq ->GetObject(("WZ_DeltaR_" + NPLc + "_" + NPLds).c_str(),hwz_dr);
		hwz_dr->SetDirectory(nullptr);// make it stay even if file closed

		zq ->GetObject(("Z_pair_jets_Delta_Phi_" + NPLc + "_" + NPLds).c_str(),hzjdph);
		hzjdph->SetDirectory(nullptr);// make it stay even if file closed

		zq ->GetObject(("npl_" + NPLc + "_" + NPLds).c_str(),h__npl);
		h__npl->SetDirectory(nullptr);// make it stay even if file closed

		zq ->GetObject((NPLc + "_" + NPLds + "_fin_jets__pt").c_str(),hjt_pt);
		hjt_pt->SetDirectory(nullptr);// make it stay even if file closed

		zq ->GetObject((NPLc + "_" + NPLds + "_fin_jets_eta").c_str(),hjteta);
		hjteta->SetDirectory(nullptr);// make it stay even if file close

		zq ->GetObject((NPLc + "_" + NPLds + "_fin_jets_phi").c_str(),hjtphi);
		hjtphi->SetDirectory(nullptr);// make it stay even if file closed

		zq ->GetObject((NPLc + "_" + NPLds + "_fin_jets_mas").c_str(),hjtmas);
		hjtmas->SetDirectory(nullptr);// make it stay even if file closed

		zq ->GetObject((NPLc + "_" + NPLds + "_lep__pt").c_str(),hlp_pt);
		hlp_pt->SetDirectory(nullptr);// make it stay even if file closed

		zq ->GetObject((NPLc + "_" + NPLds + "_lep_eta").c_str(),hlpeta);
		hlpeta->SetDirectory(nullptr);// make it stay even if file closed

		zq ->GetObject((NPLc + "_" + NPLds + "_lep_phi").c_str(),hlpphi);
		hlpphi->SetDirectory(nullptr);// make it stay even if file closed

		zq ->GetObject((NPLc + "_" + NPLds + "_lep_mas").c_str(),hlpmas);
		hlpmas->SetDirectory(nullptr);// make it stay even if file closed

		zq ->GetObject((NPLc + "_" + NPLds + "_bjet__pt").c_str(),hb__pt);
		hb__pt->SetDirectory(nullptr);// make it stay even if file closed

		zq ->GetObject((NPLc + "_" + NPLds + "_bjet_eta").c_str(),hb_eta);
		hb_eta->SetDirectory(nullptr);// make it stay even if file closed

		zq ->GetObject((NPLc + "_" + NPLds + "_bjet_phi").c_str(),hb_phi);
		hb_phi->SetDirectory(nullptr);// make it stay even if file closed

		zq ->GetObject((NPLc + "_" + NPLds + "_bjet_mas").c_str(),hb_mas);
		hb_mas->SetDirectory(nullptr);// make it stay even if file closed

		zq.Close();
	/*NPLds = "cms"; // now do the odd one out
		temp_opener = temp_header + NPLc + "_" + NPLds + temp_footer;
		std::cout << "Opening file " << temp_opener << std::endl;
		TFile dc(temp_opener.c_str());
		if( ! dc.IsOpen()) throw std::runtime_error("File is not opened");
		dc.GetObject(("N_data_LnT_" + NPLc + "_" + NPLds).c_str(),dcm);
		if(!dcm) throw std::runtime_error("N_data_LnT not found");
		dcm->SetDirectory(nullptr);// make it stay even if file closed
		dc.Close();*/
	for(dataSource ds:dataSourceAll){// to go through all ds get the files open
		if(tzq == ds || met == ds || cms == ds)continue;
		switch (ds){// only needs to be done for tz1 and tz2 for calchisto.cpp
		case  tzq:{NPLds =  "tzq";break;}
		case   ww:{NPLds =   "ww";break;}
		case   wz:{NPLds =   "wz";break;}
		case   zz:{NPLds =   "zz";break;}
		case   st:{NPLds =   "st";break;}
		case  stb:{NPLds =  "stb";break;}
		case  stw:{NPLds =  "stw";break;}
		case stbw:{NPLds = "stbw";break;}
		case wjqq:{NPLds = "wjqq";break;}
		case wzll:{NPLds = "wzll";break;}
		case  ttb:{NPLds =  "ttb";break;}
		case  ttl:{NPLds =  "ttl";break;}
		case  ttj:{NPLds =  "ttj";break;}
		case  wjt:{NPLds =  "wjt";break;}
		case  tz1:{NPLds =  "tz1";break;}
		case  tz2:{NPLds =  "tz2";break;}
		case  met:{NPLds =  "met";break;}
		case  cms:{NPLds =  "cms";break;}
		default  :throw std::invalid_argument(
			"Unimplemented ds (NPL file reading)");
		}
		temp_opener = temp_header + NPLc + "_" + NPLds + temp_footer;
		if(0<debug) std::cout << "Opening file " << temp_opener << std::endl;
		TFile tf(temp_opener.c_str());
		if( ! tf.IsOpen()) throw std::runtime_error("File is not opened");

		tF ->GetObject(("btag_w_" + NPLc + "_" + NPLds).c_str(),hbtagw);
		fbtagw->Add(hbtagw);

		tF ->GetObject(("mostSF_" + NPLc + "_" + NPLds).c_str(),hlp_sf);
		flp_sf->Add(hlp_sf);

		tF ->GetObject(( "ttbSF_" + NPLc + "_" + NPLds).c_str(),htpt_w);
		ftpt_w->Add(htbt_w);

		tF ->GetObject(("cmet_sEt_" + NPLc + "_" + NPLds).c_str(),hcmtet);
		fcmtet->Add(hcmtet);

		tF ->GetObject(("cmet__pt_" + NPLc + "_" + NPLds).c_str(),hcmtpt);
		fcmtpt->Add(hcmtpt);

		tF ->GetObject(("cmet_phi_" + NPLc + "_" + NPLds).c_str(),hcmtph);
		fcmtph->Add(hcmtph);

		tF ->GetObject(("met_sEt_" + NPLc + "_" + NPLds).c_str(),hmteta);
		fmteta->Add(hmteta);

		tF ->GetObject(("met__pt_" + NPLc + "_" + NPLds).c_str(),hmt_pt);
		fmt_pt->Add(hmt_pt);

		tF ->GetObject(("ttop_pt_" + NPLc + "_" + NPLds).c_str(),ht__pt);
		ft__pt->Add(ht__pt);

		tF ->GetObject(("tTm_" + NPLc + "_" + NPLds).c_str(),ht_mas);
		ft_mas->Add(ht_mas);

		tF ->GetObject(("W_invariant_mass_" + NPLc + "_" + NPLds).c_str(),hWinvm);
		fWinvm->Add(hWinvm);

		tF ->GetObject(("ev_w_" + NPLc + "_" + NPLds).c_str(),hev_sf);
		fev_sf->Add(hev_sf);

		tF ->GetObject(("zmas_" + NPLc + "_" + NPLds).c_str(),hz_mas);
		fz_mas->Add(hz_mas);

		tF ->GetObject(("Z_W_Delta_Phi_" + NPLc + "_" + NPLds).c_str(),hzwdph);
		fzwdph->Add(hzwdph);

		tF ->GetObject(("Z_MET_Delta_Phi_" + NPLc + "_" + NPLds).c_str(),hzmdph);
		fzmdph->Add(hzmdph);

		tF ->GetObject(("WZ_DeltaR_" + NPLc + "_" + NPLds).c_str(),hwz_dr);
		fwz_dr->Add(hwz_dr);

		tF ->GetObject(("Z_pair_jets_Delta_Phi_" + NPLc + "_" + NPLds).c_str(),hzjdph);
		fzjdph->Add(hzjdph);

		tF ->GetObject(("npl_" + NPLc + "_" + NPLds).c_str(),h__npl);
		f__npl->Add(h__npl);

		tF ->GetObject((NPLc + "_" + NPLds + "_fin_jets__pt").c_str(),hjt_pt);
		fjt_pt->Add(hjt_pt);

		tF ->GetObject((NPLc + "_" + NPLds + "_fin_jets_eta").c_str(),hjteta);
		fjteta->Add(hjteta);

		tF ->GetObject((NPLc + "_" + NPLds + "_fin_jets_phi").c_str(),hjtphi);
		fjtphi->Add(hjtphi);

		tF ->GetObject((NPLc + "_" + NPLds + "_fin_jets_mas").c_str(),hjtmas);
		fjtmas->Add(hjtmas);

		tF ->GetObject((NPLc + "_" + NPLds + "_lep__pt").c_str(),hlp_pt);
		flp_pt->Add(hlp_pt);

		tF ->GetObject((NPLc + "_" + NPLds + "_lep_eta").c_str(),hlpeta);
		flpeta->Add(hlpeta);

		tF ->GetObject((NPLc + "_" + NPLds + "_lep_phi").c_str(),hlpphi);
		flpphi->Add(hlpphi);

		tF ->GetObject((NPLc + "_" + NPLds + "_lep_mas").c_str(),hlpmas);
		flpmas->Add(hlpmas);

		tF ->GetObject((NPLc + "_" + NPLds + "_bjet__pt").c_str(),hb__pt);
		fb__pt->Add(hb__pt);

		tF ->GetObject((NPLc + "_" + NPLds + "_bjet_eta").c_str(),hb_eta);
		fb_eta->Add(hb_eta);

		tF ->GetObject((NPLc + "_" + NPLds + "_bjet_phi").c_str(),hb_phi);
		fb_phi->Add(hb_phi);

		tF ->GetObject((NPLc + "_" + NPLds + "_bjet_mas").c_str(),hb_mas);
		fb_mas->Add(hb_mas);

		tf.Close();
	}// for
	// try to associate pointers correctly and store them
	if(0<debug) std::cout<<"all objects added"<<std::endl;
	TFile hf(("histo/NPL_run_" + NPLc + ".root").c_str(),"RECREATE");
	if(0<debug) std::cout<<"file created"<<std::endl;
	fbtagw->SetName(("btag_w_"   + NPLc +"_NPL".root").c_str());
	flp_sf->SetName(("mostSF_"   + NPLc +"_NPL".root").c_str());
	ftpt_w->SetName(( "ttbSF_"   + NPLc +"_NPL".root").c_str());
	fcmtet->SetName(("cmet_sEt_" + NPLc +"_NPL".root").c_str());
	fcmtpt->SetName(("cmet__pt_" + NPLc +"_NPL".root").c_str());
	fcmtph->SetName(("cmet_phi_" + NPLc +"_NPL".root").c_str());
	fmteta->SetName(( "met_sEt_" + NPLc +"_NPL".root").c_str());
	fmt_pt->SetName(( "met__pt_" + NPLc +"_NPL".root").c_str());
	ft__pt->SetName(( "ttop_pt_" + NPLc +"_NPL".root").c_str());
	ft_mas->SetName((     "tTm_" + NPLc +"_NPL".root").c_str());
	fWinvm->SetName(("W_invariant_mass_" + NPLc +"_NPL".root").c_str());
	fev_sf->SetName((     "ev_w" + NPLc +"_NPL".root").c_str());
	fz_mas->SetName((    "zmas_" + NPLc +"_NPL".root").c_str());
	fzwdph->SetName((  "Z_W_Delta_Phi" + NPLc +"_NPL".root").c_str());
	fzmdph->SetName(("Z_MET_Delta_Phi" + NPLc +"_NPL".root").c_str());
	fwz_dr->SetName(("WZ_DeltaR" + NPLc +"_NPL".root").c_str());
	fzjdph->SetName(("Z_pair_jets_Delta_Phi" + NPLc +"_NPL".root").c_str());
	f__npl->SetName(("npl_" + NPLc +"_NPL".root").c_str());
	fjt_pt->SetName((NPLc +"_NPL_" + "fin_jets__pt".root").c_str());
	fjteta->SetName((NPLc +"_NPL_" + "fin_jets_eta".root").c_str());
	fjtphi->SetName((NPLc +"_NPL_" + "fin_jets_phi".root").c_str());
	fjtmas->SetName((NPLc +"_NPL_" + "fin_jets_mas".root").c_str());
	flp_pt->SetName((NPLc +"_NPL_" + "lep__pt".root").c_str());
	flpeta->SetName((NPLc +"_NPL_" + "lep_eta".root").c_str());
	flpphi->SetName((NPLc +"_NPL_" + "lep_phi".root").c_str());
	flpmas->SetName((NPLc +"_NPL_" + "lep_mas".root").c_str());
	fb__pt->SetName((NPLc +"_NPL_" + "bjet__pt".root").c_str());
	fb_eta->SetName((NPLc +"_NPL_" + "bjet_eta".root").c_str());
	fb_phi->SetName((NPLc +"_NPL_" + "bjet_phi".root").c_str());
	fb_mas->SetName((NPLc +"_NPL_" + "bjet_mas".root").c_str());

	hf.WriteTObject(fbtagw);hf.Flush();sync();
	hf.WriteTObject(flp_sf);hf.Flush();sync();
	hf.WriteTObject(ftpt_w);hf.Flush();sync();
	hf.WriteTObject(fcmtet);hf.Flush();sync();
	hf.WriteTObject(fcmtpt);hf.Flush();sync();
	hf.WriteTObject(fcmtph);hf.Flush();sync();
	hf.WriteTObject(fmteta);hf.Flush();sync();
	hf.WriteTObject(fmt_pt);hf.Flush();sync();
	hf.WriteTObject(ft__pt);hf.Flush();sync();
	hf.WriteTObject(ft_mas);hf.Flush();sync();
	hf.WriteTObject(fWinvm);hf.Flush();sync();
	hf.WriteTObject(fev_sf);hf.Flush();sync();
	hf.WriteTObject(fz_mas);hf.Flush();sync();
	hf.WriteTObject(fzwdph);hf.Flush();sync();
	hf.WriteTObject(fzmdph);hf.Flush();sync();
	hf.WriteTObject(fwz_dr);hf.Flush();sync();
	hf.WriteTObject(fzjdph);hf.Flush();sync();
	hf.WriteTObject(f__npl);hf.Flush();sync();
	hf.WriteTObject(fjt_pt);hf.Flush();sync();
	hf.WriteTObject(fjteta);hf.Flush();sync();
	hf.WriteTObject(fjtphi);hf.Flush();sync();
	hf.WriteTObject(fjtmas);hf.Flush();sync();
	hf.WriteTObject(flp_pt);hf.Flush();sync();
	hf.WriteTObject(flpeta);hf.Flush();sync();
	hf.WriteTObject(flpphi);hf.Flush();sync();
	hf.WriteTObject(flpmas);hf.Flush();sync();
	hf.WriteTObject(fb__pt);hf.Flush();sync();
	hf.WriteTObject(fb_eta);hf.Flush();sync();
	hf.WriteTObject(fb_phi);hf.Flush();sync();
	hf.WriteTObject(fb_mas);hf.Flush();sync();

	if(0<debug) std::cout<<"all objects stored on file"<<std::endl;
}// void

int main ( int argc , char *argv[] ){
	addhists(elnu);
	addhists(munu);
	return 0;
}
