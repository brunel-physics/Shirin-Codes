#include <ROOT/RDataFrame.hxx>// big guns
#include <TCanvas.h>
#include <TLegend.h>
#include <THStack.h>
//#include <TStyle.h>

#include "tdrstyle.C"
#include "src/calchisto.hpp"
//#include "calchisto.cpp" // only for compiling reasons

int plotstacks(){
	setTDRStyle();
	
	std::string opener="elnu_tzq";// become 2 for loops
	TFile hf((opener+".histo").c_str());
	//TH1D *htransT, *htransW, *hzmetdph, *hzwdph, *hzjetdphi, *hWinvmass;
	for(std::string particle:{"fin_jets","lep","bjet"}){
	
	for(PtEtaPhiM k:PtEtaPhiMall){
		std::string kstring = "_" ;
		std::string xAxisStr;
		switch (k){
			case pt :{kstring += "_pt";
				  xAxisStr = "pT/GeV"		       ;break;}
			case eta:{kstring += "eta";
				  xAxisStr = "PseudoRapidity eta"      ;break;}
			case phi:{kstring += "phi";
				  xAxisStr = "Azimuthal angle, phi/rad";break;}
			case m  :{kstring += "mas";
				  xAxisStr = "mass GeV/C^2"	       ;break;}
		}
	std::string hobjname = (opener+"_"+particle+kstring).c_str();
	std::string   stname = (particle + kstring).c_str();
        THStack *hstack = new THStack(stname.c_str(),stname.c_str());
	TH1D *hobj;
	hf.GetObject(hobjname.c_str(),hobj);
	// clone hobj by DrawClone or add to stack here
	//THStack *hstack = new THStack(stname.c_str(),stname.c_str());
	hobj->SetLineColor(kBlack);// TBC when other ds are inlcuded.
	hstack->Add(static_cast<TH1D*> (hobj->Clone()));
	auto hcanvas =
	new TCanvas(stname.c_str(), stname.c_str(),10,10,900,900);
	hobj->GetXaxis()->SetTitle(xAxisStr.c_str());
        hobj->GetYaxis()->SetTitle("Event");
	hobj->DrawClone("SAME");
        hcanvas->cd(2);// From here downward
        hstack->Draw("HIST");// should be done
        hcanvas->BuildLegend();// once all data sources
        hcanvas->SaveAs((stname + ".root").c_str());// are included
        hcanvas->SaveAs((stname + ".pdf" ).c_str());
	}}
	
        // Adding other Plots
	
        std::string hTransTm = (opener + " Transverse Top mass").c_str();
        std::string   stname = (opener+  "_transTm").c_str();
        THStack *h_transT_stack = new THStack(stname.c_str(),stname.c_str());
        TH1D *htransT;
        hf.GetObject("h_trans_T", htransT);
        // clone hobj by DrawClone or add to stack here
        //THStack *hstack = new THStack(stname.c_str(),stname.c_str());
        htransT->SetLineColor(kBlack);// TBC when other ds are inlcuded.
        htransT->Add(static_cast<TH1D*> (htransT->Clone()));
        auto h_transT_canvas =
        new TCanvas(stname.c_str(), stname.c_str(),10,10,900,900);
        htransT->GetXaxis()->SetTitle("tranverse top mass GeV/C^2");
        htransT->GetYaxis()->SetTitle("Event");
        htransT->DrawClone("SAME");
        h_transT_canvas->cd(2);// From here downward
        h_transT_stack->Draw("HIST");// should be done
        h_transT_canvas->BuildLegend();// once all data sources
        h_transT_canvas->SaveAs((stname + ".root").c_str());// are included
        h_transT_canvas->SaveAs((stname + ".pdf" ).c_str());

// *htransW, *hzmetdph, *hzwdph, *hzjetdphi, *hWinvmass;
	return 0; // end of file
/*        TH2D *h2numer,*h2denom;
	hf.GetObject("is_numer_elnu_tzq",h2numer);
	hf.GetObject("is_denom_elnu_tzq",h2denom);
	auto h_events_is_btag_PtVsEta_canvas
	= new TCanvas("is b tag pt Vs eta" , "is b tag pt Vs eta",10,10,900,900);
	TH2D *is_btag_ratio = new TH2D("ei", "is b tag ei",50,0,400,50,-3,3);
	is_btag_ratio = static_cast<TH2D*>(h2numer->Clone());
	is_btag_ratio->GetXaxis()->SetTitle( "is b tag pt");
	is_btag_ratio->GetYaxis()->SetTitle( "is b tag eta");
	is_btag_ratio->Divide(             h2denom->Clone());//.GetPtr());
	h_events_is_btag_PtVsEta_canvas->BuildLegend();
	is_btag_ratio->Draw("COLZ");
	//h_events_is_btag_PtVsEta_canvas->SaveAs("h_events_is_btag_PtVsEta_canvas.root");
	h_events_is_btag_PtVsEta_canvas->SaveAs("h_events_is_btag_PtVsEta_canvas.pdf");
	
	hf.GetObject("no_numer_elnu_tzq",h2numer);
	hf.GetObject("no_denom_elnu_tzq",h2denom);
	auto h_events_no_btag_PtVsEta_canvas
	= new TCanvas("no b tag pt Vs eta" , "no b tag pt Vs eta",10,10,900,900);
	TH2D *no_btag_ratio = new TH2D("ej", "no b tag ej",50,0,400,50,-3,3);
	no_btag_ratio = static_cast<TH2D*>(h2numer->Clone());
	no_btag_ratio->GetXaxis()->SetTitle( "no b tag pt");
	no_btag_ratio->GetYaxis()->SetTitle( "no b tag eta");
	no_btag_ratio->Divide(             h2denom->Clone());//.GetPtr());
	h_events_no_btag_PtVsEta_canvas->BuildLegend();
	no_btag_ratio->Draw("COLZ");
	//h_events_no_btag_PtVsEta_canvas->SaveAs("h_events_no_btag_PtVsEta_canvas.root");
	h_events_no_btag_PtVsEta_canvas->SaveAs("h_events_no_btag_PtVsEta_canvas.pdf");
*/	
	hf.Close();
}

/*

///////////////////////////////////////////////////////////////////// THSTACKS ///////////////////////////////////////////////////////////////////////////////

TLegend legend_ed= TLegend(0.7,0.7,0.94,0.94);
legend_ed.SetFillStyle(1001);
legend_ed.SetBorderSize(1);
legend_ed.SetFillColor(kWhite);
//gStyle->SetOptStat(1111111);

TLegend legend_ebg= TLegend(0.7,0.7,0.94,0.94);
legend_ebg.SetFillStyle(1001);
legend_ebg.SetBorderSize(1);
legend_ebg.SetFillColor(kWhite);
//gStyle->SetOptStat(1111111);

TLegend legend_md= TLegend(0.7,0.7,0.94,0.94);
legend_md.SetFillStyle(1001);
legend_md.SetBorderSize(1);
legend_md.SetFillColor(kWhite);
//gStyle->SetOptStat(1111111);

TLegend legend_mbg= TLegend(0.7,0.7,0.94,0.94);
legend_mbg.SetFillStyle(1001);
legend_mbg.SetBorderSize(1);
legend_mbg.SetFillColor(kWhite);
//gStyle->SetOptStat(1111111);

TLegend legend_sed= TLegend(0.7,0.7,0.94,0.94);
legend_sed.SetFillStyle(1001);
legend_sed.SetBorderSize(1);
legend_sed.SetFillColor(kWhite);
//gStyle->SetOptStat(1111111);

TLegend legend_emet= TLegend(0.7,0.7,0.94,0.94);
legend_emet.SetFillStyle(1001);
legend_emet.SetBorderSize(1);
legend_emet.SetFillColor(kWhite);
//gStyle->SetOptStat(1111111);

TLegend legend_smd= TLegend(0.7,0.7,0.94,0.94);
legend_smd.SetFillStyle(1001);
legend_smd.SetBorderSize(1);
legend_smd.SetFillColor(kWhite);
gStyle->SetOptStat(1111111);

TLegend legend_mmet= TLegend(0.7,0.7,0.94,0.94);
legend_mmet.SetFillStyle(1001);
legend_mmet.SetBorderSize(1);
legend_mmet.SetFillColor(kWhite);
//gStyle->SetOptStat(1111111);

auto h_d_enu_events_ept = d_enu_P_btag.Histo1D({"MC electron_pt_enu_Channel","MC electron pt in electron-neutrino channel",50,0,250},"tight_ele_pt", "nw_tight_ele_pt");
auto h_d_munu_events_mupt = d_munu_P_btag.Histo1D({"MC muon_pt_munu_Channel","MC muon pt in muon-neutrino channel",50,0,250},"tight_mu_pt","nw_tight_mu_pt");
auto h_se_enu_events_ept = se_enu_top_selection.Histo1D({"Single Electron electron_pt_enu_Channel","Single Electron electron pt in electron-neutrino channel",50,0,250}, "tight_ele_pt");
//auto h_met_enu_events_ept = met_enu_P_btag.Histo1D({"MET electron_pt_enu_Channel","MET electron pt in electron-neutrino channel",50,0,250},"tight_ele_pt","nw_tight_ele_pt");
auto h_smu_munu_events_mupt = sm_munu_top_selection.Histo1D({"Single Muon muon_pt_munu_Channel","Single muon muon pt in muon-neutrino channel",50,0,250},"tight_mu_pt");
//auto h_met_munu_events_mupt = met_munu_P_btag.Histo1D({"MET muon_pt_munu_Channel","MET muon pt in muon-neutrino channel",50,0,250},"tight_mu_pt", "nw_tight_mu_pt");
auto h_ww_enu_events_ept = ww_enu_P_btag.Histo1D({"WW electron_pt_enu_Channel","WW electron pt in electron-neutrino channel",50,0,250},"tight_ele_pt","nw_tight_ele_pt");
auto h_wz_enu_events_ept = wz_enu_P_btag.Histo1D({"WZ electron_pt_enu_Channel","WZ electron pt in electron-neutrino channel",50,0,250},"tight_ele_pt","nw_tight_ele_pt");
auto h_zz_enu_events_ept = zz_enu_P_btag.Histo1D({"ZZ electron_pt_enu_Channel","ZZ electron pt in electron-neutrino channel",50,0,250},"tight_ele_pt","nw_tight_ele_pt");
auto h_ttZ_enu_events_ept = ttZ_enu_P_btag.Histo1D({"ttZ electron_pt_enu_Channel","ttZ electron pt in electron-neutrino channel",50,0,250},"tight_ele_pt","nw_tight_ele_pt");
auto h_ww_munu_events_mupt = ww_munu_P_btag.Histo1D({"WW muon_pt_munu_Channel","WW muon pt in muon-neutrino channel",50,0,250},"tight_mu_pt", "nw_tight_mu_pt");
auto h_wz_munu_events_mupt = wz_munu_P_btag.Histo1D({"WZ muon_pt_munu_Channel","WZ muon pt in muon-neutrino channel",50,0,250},"tight_mu_pt", "nw_tight_mu_pt");
auto h_zz_munu_events_mupt = zz_munu_P_btag.Histo1D({"ZZ muon_pt_munu_Channel","ZZ muon pt in muon-neutrino channel",50,0,250},"tight_mu_pt", "nw_tight_mu_pt");
auto h_ttZ_munu_events_mupt = ttZ_munu_P_btag.Histo1D({"ttZ muon_pt_munu_Channel","ttZ muon pt in muon-neutrino channel",50,0,250},"tight_mu_pt", "nw_tight_mu_pt");

THStack *lep_pt_Stack = new THStack("MC_Stack","Lepton Transverse Momentum ");
h_d_enu_events_ept->SetLineColor(kBlack);
h_d_munu_events_mupt->SetLineColor(kGreen);
h_se_enu_events_ept->SetLineColor(kPink);
//h_met_enu_events_ept->SetLineColor(kCherry);
h_smu_munu_events_mupt->SetLineColor(kViolet);
//h_met_munu_events_mupt->SetLineColor(kRose);
h_ww_enu_events_ept->SetLineColor(kGray);
h_wz_enu_events_ept->SetLineColor(kRed);
h_zz_enu_events_ept->SetLineColor(kBlue);
h_ttZ_enu_events_ept->SetLineColor(kCyan);
h_ww_munu_events_mupt->SetLineColor(kOrange);
h_wz_munu_events_mupt->SetLineColor(kSpring);
h_zz_munu_events_mupt->SetLineColor(kAzure);
h_ttZ_munu_events_mupt->SetLineColor(kTeal);

lep_pt_Stack->Add((TH1*)&h_d_enu_events_ept.GetValue());
lep_pt_Stack->Add((TH1*)&h_d_munu_events_mupt.GetValue());
lep_pt_Stack->Add((TH1*)&h_se_enu_events_ept.GetValue());
//lep_pt_Stack->Add((TH1*)&h_met_enu_events_ept.GetValue());
lep_pt_Stack->Add((TH1*)&h_smu_munu_events_mupt.GetValue());
//lep_pt_Stack->Add((TH1*)&h_met_munu_events_mupt.GetValue());
lep_pt_Stack->Add((TH1*)&h_ww_enu_events_ept.GetValue());
lep_pt_Stack->Add((TH1*)&h_wz_enu_events_ept.GetValue());
lep_pt_Stack->Add((TH1*)&h_zz_enu_events_ept.GetValue());
lep_pt_Stack->Add((TH1*)&h_ttZ_enu_events_ept.GetValue());
lep_pt_Stack->Add((TH1*)&h_ww_munu_events_mupt.GetValue());
lep_pt_Stack->Add((TH1*)&h_wz_munu_events_mupt.GetValue());
lep_pt_Stack->Add((TH1*)&h_zz_munu_events_mupt.GetValue());
lep_pt_Stack->Add((TH1*)&h_ttZ_munu_events_mupt.GetValue());

auto h_events_lep_pt_canvas = new TCanvas("electron pt", "electron pt",10,10,900,900);
h_d_enu_events_ept->GetXaxis()->SetTitle("Pt/GeV");
h_d_enu_events_ept->GetYaxis()->SetTitle("Events");
legend_ed.AddEntry(h_d_enu_events_ept.GetPtr(),"tZq MC,electron pt","l");
legend_md.AddEntry(h_d_munu_events_mupt.GetPtr(),"tZq MC,muon pt","l");
legend_sed.AddEntry(h_se_enu_events_ept.GetPtr(),"Single electron pt","l");
//legend_emet.AddEntry(h_se_enu_events_ept.GetPtr(),"MET data,electron pt","l");
legend_smd.AddEntry(h_smu_munu_events_mupt.GetPtr(),"Signal muon pt","l");
//legend_mmet.AddEntry(h_met_munu_events_mupt.GetPtr(),"MET data,muon pt","l");

//h_d_enu_events_ept->Fit("jacob_fit");
h_d_enu_events_ept->Draw();
h_d_munu_events_mupt->Draw("SAME");
h_se_enu_events_ept->Draw("SAME");
//h_met_enu_events_ept->Draw("SAME");
h_smu_munu_events_mupt->Draw("SAME");
//h_met_munu_events_mupt->Draw("SAME");
h_ww_enu_events_ept->Draw("SAME");
h_wz_enu_events_ept->Draw("SAME");
h_zz_enu_events_ept->Draw("SAME");
h_ttZ_enu_events_ept->Draw("SAME");
h_ww_munu_events_mupt->Draw("SAME");
h_wz_munu_events_mupt->Draw("SAME");
h_zz_munu_events_mupt->Draw("SAME");
h_ttZ_munu_events_mupt->Draw("SAME");

h_events_lep_pt_canvas->cd(2);
lep_pt_Stack->Draw("HIST");

h_events_lep_pt_canvas->BuildLegend();
h_events_lep_pt_canvas->SaveAs("hist_lep_pt.root");
h_events_lep_pt_canvas->SaveAs("hist_lep_pt.pdf");
legend_ed.Clear();
legend_md.Clear();
legend_sed.Clear();
//legend_emet.Clear();
legend_smd.Clear();
//legend_mmet.Clear();
*/
