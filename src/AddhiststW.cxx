//clang++ -Isrc -std=c++17 -march=native -pipe -O3 -Wall -Wextra -Wpedantic -o build/addhiststW src/AddhiststW.cxx `root-config --libs` -lm
#include <ROOT/RDataFrame.hxx>//#include <ROOT/RCsvDS.hxx>

#include "calchisto.hpp"

int debug = 1;

namespace{
	std::string
	NPLc = "elnu",
	NPLds = "stw",
	temp_header = "histo/",
	temp_footer = ".root",
	temp_opener
	;
	TH1D *hbtagw, *hlp_sf, *htpt_w,
 	     *hcmtet, *hcmtpt, *hcmtph,
	     *hmteta, *hmt_pt, *htWinm,
	     *ht__pt, *ht_mas, *hWinvm,
	     *hev_sf, *hz_mas, *hzwdph,
	     *hzmdph, *hwz_dr, *hzjdph,
	     *h__npl, *hjt_pt, *hjteta,
	     *hjtphi, *hjtmas, *hlp_pt,
	     *hlpeta, *hlpphi, *hlpmas,
	     *hb__pt, *hb_eta, *hb_phi,
	     *hb_mas;// for taking from files
	TH1D *fbtagw, *flp_sf, *ftpt_w,
             *fcmtet, *fcmtpt, *fcmtph,
             *fmteta, *fmt_pt, *ftWinm,
             *ft__pt, *ft_mas, *fWinvm,
             *fev_sf, *fz_mas, *fzwdph,
             *fzmdph, *fwz_dr, *fzjdph,
             *f__npl, *fjt_pt, *fjteta,
             *fjtphi, *fjtmas, *flp_pt,
             *flpeta, *flpphi, *flpmas,
             *fb__pt, *fb_eta, *fb_phi,
             *fb_mas;      // finals go here
}

void addhiststW(const channel ch){
	// do tzq first so that tpr & tln are not null
	if(munu == ch) {NPLc = "munu";NPLds = "stw";}
		temp_opener = temp_header + NPLc + "_" + NPLds + temp_footer;
		std::cout << "Opening file " << temp_opener << std::endl;
		TFile zq(temp_opener.c_str());
		if( ! zq.IsOpen()) throw std::runtime_error("File is not opened");
		std::cout<<"file is open"<<std::endl;
		zq.GetObject(("btag_w_" + NPLc + "_" + NPLds).c_str(),fbtagw);
		fbtagw->SetDirectory(nullptr);// make it stay even if file closed

		zq.GetObject(("mostSF_" + NPLc + "_" + NPLds).c_str(),flp_sf);
		flp_sf->SetDirectory(nullptr);// make it stay even if file closed

		zq.GetObject(( "ttbSF_" + NPLc + "_" + NPLds).c_str(),ftpt_w);
		ftpt_w->SetDirectory(nullptr);// make it stay even if file closed

		zq.GetObject(("cmet_sEt_" + NPLc + "_" + NPLds).c_str(),fcmtet);
		fcmtet->SetDirectory(nullptr);// make it stay even if file close

		zq.GetObject(("cmet__pt_" + NPLc + "_" + NPLds).c_str(),fcmtpt);
		fcmtpt->SetDirectory(nullptr);// make it stay even if file closed

		zq.GetObject(("cmet_phi_" + NPLc + "_" + NPLds).c_str(),fcmtph);
		fcmtph->SetDirectory(nullptr);// make it stay even if file closed

		zq.GetObject(("met_sEt_" + NPLc + "_" + NPLds).c_str(),fmteta);
		fmteta->SetDirectory(nullptr);// make it stay even if file closed

		zq.GetObject(("met__pt_" + NPLc + "_" + NPLds).c_str(),fmt_pt);
		fmt_pt->SetDirectory(nullptr);// make it stay even if file closed

		zq.GetObject(("ttop_pt_" + NPLc + "_" + NPLds).c_str(),ft__pt);
		ft__pt->SetDirectory(nullptr);// make it stay even if file closed

		zq.GetObject(("tTm_" + NPLc + "_" + NPLds).c_str(),ft_mas);
		ft_mas->SetDirectory(nullptr);// make it stay even if file closed

		zq.GetObject(("W_invariant_mass_" + NPLc + "_" + NPLds).c_str(),fWinvm);
		fWinvm->SetDirectory(nullptr);// make it stay even if file closed

                zq.GetObject(("tWm_" + NPLc + "_" + NPLds).c_str(),ftWinm);
                ftWinm->SetDirectory(nullptr);

		zq.GetObject(("ev_w_" + NPLc + "_" + NPLds).c_str(),fev_sf);
		fev_sf->SetDirectory(nullptr);// make it stay even if file closed

		zq.GetObject(("zmas_" + NPLc + "_" + NPLds).c_str(),fz_mas);
		fz_mas->SetDirectory(nullptr);// make it stay even if file closed

		zq.GetObject(("Z_W_Delta_Phi_" + NPLc + "_" + NPLds).c_str(),fzwdph);
		fzwdph->SetDirectory(nullptr);// make it stay even if file closed

		zq.GetObject(("Z_MET_Delta_Phi_" + NPLc + "_" + NPLds).c_str(),fzmdph);
		fzmdph->SetDirectory(nullptr);// make it stay even if file closed

		zq.GetObject(("WZ_DeltaR_" + NPLc + "_" + NPLds).c_str(),fwz_dr);
		fwz_dr->SetDirectory(nullptr);// make it stay even if file closed

		zq.GetObject(("Z_pair_jets_Delta_Phi_" + NPLc + "_" + NPLds).c_str(),fzjdph);
		fzjdph->SetDirectory(nullptr);// make it stay even if file closed

		zq.GetObject(("npl_" + NPLc + "_" + NPLds).c_str(),f__npl);
		f__npl->SetDirectory(nullptr);// make it stay even if file closed

		zq.GetObject((NPLc + "_" + NPLds + "_fin_jets__pt").c_str(),fjt_pt);
		fjt_pt->SetDirectory(nullptr);// make it stay even if file closed

		zq.GetObject((NPLc + "_" + NPLds + "_fin_jets_eta").c_str(),fjteta);
		fjteta->SetDirectory(nullptr);// make it stay even if file close

		zq.GetObject((NPLc + "_" + NPLds + "_fin_jets_phi").c_str(),fjtphi);
		fjtphi->SetDirectory(nullptr);// make it stay even if file closed

		zq.GetObject((NPLc + "_" + NPLds + "_fin_jets_mas").c_str(),fjtmas);
		fjtmas->SetDirectory(nullptr);// make it stay even if file closed

		zq.GetObject((NPLc + "_" + NPLds + "_lep__pt").c_str(),flp_pt);
		flp_pt->SetDirectory(nullptr);// make it stay even if file closed

		zq.GetObject((NPLc + "_" + NPLds + "_lep_eta").c_str(),flpeta);
		flpeta->SetDirectory(nullptr);// make it stay even if file closed

		zq.GetObject((NPLc + "_" + NPLds + "_lep_phi").c_str(),flpphi);
		flpphi->SetDirectory(nullptr);// make it stay even if file closed

		zq.GetObject((NPLc + "_" + NPLds + "_lep_mas").c_str(),flpmas);
		flpmas->SetDirectory(nullptr);// make it stay even if file closed

		zq.GetObject((NPLc + "_" + NPLds + "_bjet__pt").c_str(),fb__pt);
		fb__pt->SetDirectory(nullptr);// make it stay even if file closed

		zq.GetObject((NPLc + "_" + NPLds + "_bjet_eta").c_str(),fb_eta);
		fb_eta->SetDirectory(nullptr);// make it stay even if file closed

		zq.GetObject((NPLc + "_" + NPLds + "_bjet_phi").c_str(),fb_phi);
		fb_phi->SetDirectory(nullptr);// make it stay even if file closed

		zq.GetObject((NPLc + "_" + NPLds + "_bjet_mas").c_str(),fb_mas);
		fb_mas->SetDirectory(nullptr);// make it stay even if file closed

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
	for(dataSource ds:{stw,stbw}){// to go through all ds get the files open
		if(stw == ds)continue;
		switch (ds){// only needs to be done for tz1 and tz2 for calchisto.cpp
		case  tzq:{NPLds =  "tzq";break;}
		case   ww:{NPLds =   "_ww";break;}
		case   wz:{NPLds =   "_wz";break;}
		case   zz:{NPLds =   "_zz";break;}
		case   st:{NPLds =   "_st";break;}
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
		std::cout<<"0"<<std::endl;
		tf.GetObject(("btag_w_" + NPLc + "_" + NPLds).c_str(),hbtagw);
                if(!hbtagw) throw std::runtime_error("hbtagw not found");
		std::cout<<"0.5"<< ("btag_w_" + NPLc + "_" + NPLds).c_str() <<std::endl;
		fbtagw->Add(hbtagw);
                std::cout<<"1"<<std::endl;
		tf.GetObject(("mostSF_" + NPLc + "_" + NPLds).c_str(),hlp_sf);
		flp_sf->Add(hlp_sf);

		tf.GetObject(( "ttbSF_" + NPLc + "_" + NPLds).c_str(),htpt_w);
		ftpt_w->Add(htpt_w);

		tf.GetObject(("cmet_sEt_" + NPLc + "_" + NPLds).c_str(),hcmtet);
		fcmtet->Add(hcmtet);

		tf.GetObject(("cmet__pt_" + NPLc + "_" + NPLds).c_str(),hcmtpt);
		fcmtpt->Add(hcmtpt);

		tf.GetObject(("cmet_phi_" + NPLc + "_" + NPLds).c_str(),hcmtph);
		fcmtph->Add(hcmtph);

		tf.GetObject(("met_sEt_" + NPLc + "_" + NPLds).c_str(),hmteta);
		fmteta->Add(hmteta);

		tf.GetObject(("met__pt_" + NPLc + "_" + NPLds).c_str(),hmt_pt);
		fmt_pt->Add(hmt_pt);

		tf.GetObject(("ttop_pt_" + NPLc + "_" + NPLds).c_str(),ht__pt);
		ft__pt->Add(ht__pt);

		tf.GetObject(("tTm_" + NPLc + "_" + NPLds).c_str(),ht_mas);
		ft_mas->Add(ht_mas);

		tf.GetObject(("W_invariant_mass_" + NPLc + "_" + NPLds).c_str(),hWinvm);
		fWinvm->Add(hWinvm);

                tf.GetObject(("tWm_" + NPLc + "_" + NPLds).c_str(),htWinm);
                ftWinm->Add(htWinm);

		tf.GetObject(("ev_w_" + NPLc + "_" + NPLds).c_str(),hev_sf);
		fev_sf->Add(hev_sf);

		tf.GetObject(("zmas_" + NPLc + "_" + NPLds).c_str(),hz_mas);
		fz_mas->Add(hz_mas);

		tf.GetObject(("Z_W_Delta_Phi_" + NPLc + "_" + NPLds).c_str(),hzwdph);
		fzwdph->Add(hzwdph);

		tf.GetObject(("Z_MET_Delta_Phi_" + NPLc + "_" + NPLds).c_str(),hzmdph);
		fzmdph->Add(hzmdph);

		tf.GetObject(("WZ_DeltaR_" + NPLc + "_" + NPLds).c_str(),hwz_dr);
		fwz_dr->Add(hwz_dr);

		tf.GetObject(("Z_pair_jets_Delta_Phi_" + NPLc + "_" + NPLds).c_str(),hzjdph);
		fzjdph->Add(hzjdph);

		tf.GetObject(("npl_" + NPLc + "_" + NPLds).c_str(),h__npl);
		f__npl->Add(h__npl);

		tf.GetObject((NPLc + "_" + NPLds + "_fin_jets__pt").c_str(),hjt_pt);
		fjt_pt->Add(hjt_pt);

		tf.GetObject((NPLc + "_" + NPLds + "_fin_jets_eta").c_str(),hjteta);
		fjteta->Add(hjteta);

		tf.GetObject((NPLc + "_" + NPLds + "_fin_jets_phi").c_str(),hjtphi);
		fjtphi->Add(hjtphi);

		tf.GetObject((NPLc + "_" + NPLds + "_fin_jets_mas").c_str(),hjtmas);
		fjtmas->Add(hjtmas);

		tf.GetObject((NPLc + "_" + NPLds + "_lep__pt").c_str(),hlp_pt);
		flp_pt->Add(hlp_pt);

		tf.GetObject((NPLc + "_" + NPLds + "_lep_eta").c_str(),hlpeta);
		flpeta->Add(hlpeta);

		tf.GetObject((NPLc + "_" + NPLds + "_lep_phi").c_str(),hlpphi);
		flpphi->Add(hlpphi);

		tf.GetObject((NPLc + "_" + NPLds + "_lep_mas").c_str(),hlpmas);
		flpmas->Add(hlpmas);

		tf.GetObject((NPLc + "_" + NPLds + "_bjet__pt").c_str(),hb__pt);
		fb__pt->Add(hb__pt);

		tf.GetObject((NPLc + "_" + NPLds + "_bjet_eta").c_str(),hb_eta);
		fb_eta->Add(hb_eta);

		tf.GetObject((NPLc + "_" + NPLds + "_bjet_phi").c_str(),hb_phi);
		fb_phi->Add(hb_phi);

		tf.GetObject((NPLc + "_" + NPLds + "_bjet_mas").c_str(),hb_mas);
		fb_mas->Add(hb_mas);

		tf.Close();
	}// for
	// try to associate pointers correctly and store them
	if(0<debug) std::cout<<"all objects added"<<std::endl;
	TFile hf(("histo/" + NPLc + "__tW.root").c_str(),"RECREATE");
	if(0<debug) std::cout<<"file created"<<std::endl;
	fbtagw->SetName(("btag_w_"   + NPLc +"__tW").c_str());
	flp_sf->SetName(("mostSF_"   + NPLc +"__tW").c_str());
	ftpt_w->SetName(( "ttbSF_"   + NPLc +"__tW").c_str());
	fcmtet->SetName(("cmet_sEt_" + NPLc +"__tW").c_str());
	fcmtpt->SetName(("cmet__pt_" + NPLc +"__tW").c_str());
	fcmtph->SetName(("cmet_phi_" + NPLc +"__tW").c_str());
	fmteta->SetName(( "met_sEt_" + NPLc +"__tW").c_str());
	fmt_pt->SetName(( "met__pt_" + NPLc +"__tW").c_str());
	ft__pt->SetName(( "ttop_pt_" + NPLc +"__tW").c_str());
	ft_mas->SetName((     "tTm_" + NPLc +"__tW").c_str());
	fWinvm->SetName(("W_invariant_mass_" + NPLc +"__tW").c_str());
        ftWinm->SetName((     "tWm_" + NPLc +"__tW").c_str());
	fev_sf->SetName((     "ev_w_"+ NPLc +"__tW").c_str());
	fz_mas->SetName((    "zmas_" + NPLc +"__tW").c_str());
	fzwdph->SetName((  "Z_W_Delta_Phi_" + NPLc +"__tW").c_str());
	fzmdph->SetName(("Z_MET_Delta_Phi_" + NPLc +"__tW").c_str());
	fwz_dr->SetName(("WZ_DeltaR_"+ NPLc +"__tW").c_str());
	fzjdph->SetName(("Z_pair_jets_Delta_Phi_" + NPLc +"__tW").c_str());
	f__npl->SetName(("npl_" + NPLc +"__tW").c_str());
	fjt_pt->SetName((NPLc +"__tW_" + "fin_jets__pt").c_str());
	fjteta->SetName((NPLc +"__tW_" + "fin_jets_eta").c_str());
	fjtphi->SetName((NPLc +"__tW_" + "fin_jets_phi").c_str());
	fjtmas->SetName((NPLc +"__tW_" + "fin_jets_mas").c_str());
	flp_pt->SetName((NPLc +"__tW_" + "lep__pt").c_str());
	flpeta->SetName((NPLc +"__tW_" + "lep_eta").c_str());
	flpphi->SetName((NPLc +"__tW_" + "lep_phi").c_str());
	flpmas->SetName((NPLc +"__tW_" + "lep_mas").c_str());
	fb__pt->SetName((NPLc +"__tW_" + "bjet__pt").c_str());
	fb_eta->SetName((NPLc +"__tW_" + "bjet_eta").c_str());
	fb_phi->SetName((NPLc +"__tW_" + "bjet_phi").c_str());
	fb_mas->SetName((NPLc +"__tW_" + "bjet_mas").c_str());

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
        hf.WriteTObject(ftWinm);hf.Flush();sync();
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
	addhiststW(elnu);
	addhiststW(munu);
	return 0;
}
