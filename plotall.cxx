#include <ROOT/RDataFrame.hxx>// big guns
#include <TCanvas.h>
#include <TLegend.h>
#include <THStack.h>
#include <iterator>// just for std::size

#include "src/tdrstyle.C"
#include "src/calchisto.hpp"

std::string allNamesArray[][3] = {// histogram id, histogram title, x axis string
	 //{"cmet_sEt_",   "Corrected MET " "sum E_{T} ","Sum E_{T} (GeV)"}
	 { "met_sEt_","Un-corrected MET " "sum E_{T} ","Sum E_{T} (GeV)"}
	//,{"cmet__pt_",   "Corrected MET "     "p_{T} ","p_{T} (GeV/c)"}
	//,{"cmet_dpx_",   "Corrected MET #Delta p_{x} ","p_{x} (GeV/c)"}
	//,{"cmet_dpy_",   "Corrected MET #Delta p_{y} ","p_{y} (GeV/c)"}
	,{ "met__pt_","Un-corrected MET "     "p_{T} ","p_{T} (GeV/c)"}
	,{"W_invariant_mass_","W invariant mass ","\\text{W m_{ } (GeV/}c^{2})"}
	,{     "tWm_","Transverse W mass ","\\text{W m_{T} (GeV/}c^{2})"}
	,{     "tTm_","Transverse T mass ","\\text{T m_{T} (GeV/}c^{2})"}
	,{    "zmas_",    "Recon. Z mass ","\\text{Z m_{ } (GeV/}c^{2})"}
	,{        "Z_W_Delta_Phi_",      "Z  W  #Delta#phi ",    "Z &  W  #Delta#phi (rad)"}
	,{      "Z_MET_Delta_Phi_",      "Z MET #Delta#phi ",    "Z & MET #Delta#phi (rad)"}
	,{"Z_pair_jets_Delta_Phi_","Z pair jets #Delta#phi ","Z pair jets #Delta#phi (rad)"}
	,{"ev_w_", "Event Weight ","Weight"}
};

int plotall(){
	gROOT->SetBatch(kTRUE);// no open canvas window
	setTDRStyle();
	TFile     cf("plots/plotall.root","RECREATE");
	TCanvas canv("name to reset","title to reset",10,10,900,900);
	TH1D   *hobj;
	// now we open ALL the files
	std::map<std::pair<channel,dataSource>,TFile*> hFd;
	for(channel ch:{elnu}){//channelAll){//for(channel ch:channelAll){
	std::string chN;
	switch     (ch){
		case elnu:  {chN ="elnu_";break;}
		case munu:  {chN ="munu_";break;}
	}
	for(dataSource ds:dataSourceAll){
	std::string  opener  =  chN ;
	switch  (ds){
		case tzq:{opener += "tzq";break;}
		case  ww:{opener += "_ww";break;}
		case  wz:{opener += "_wz";break;}
		case  zz:{opener += "_zz";break;}
		case ttz:{opener += "ttz";break;}
		case ttb:{opener += "ttb";break;}
		case met:{opener += "met";break;}
		case cms:{opener += "cms";break;}
	}
	hFd[std::make_pair(ch,ds)]
		= new TFile(("histo/" + opener + ".root").c_str());
	}}// now we have a histogram file dictionary of all the files miahahaha

	for(size_t i=0; i < std::size(allNamesArray) ;++i){
	for(channel ch:{elnu}){//channelAll){
	std::string chN;
	switch     (ch){
		case elnu:  {chN ="elnu";break;}
		case munu:  {chN ="munu";break;}
	}// switch
	std::string   title = allNamesArray[i][1] + chN;
	std::string  stname = allNamesArray[i][0] + chN;
	canv.SetName(stname.c_str());canv.SetTitle(stname.c_str());
	THStack stac(stname.c_str(),title.c_str());
	for(dataSource ds:dataSourceAll){
	std::string  opener  = chN + "_";
	int colour;
	switch  (ds){
		case tzq:{opener += "tzq";colour = 6;break;}// magenta
		case  ww:{opener += "_ww";colour = 2;break;}// red
		case  wz:{opener += "_wz";colour = 3;break;}// green
		case  zz:{opener += "_zz";colour = 4;break;}// blue
		case ttz:{opener += "ttz";colour = 5;break;}// yellow
		case ttb:{opener += "ttb";colour = 7;break;}// cyan
		case met:{opener += "met";colour = 9;break;}// violet
		case cms:{opener += "cms";colour = 1;break;}// black
	}
	std::string    hobjN = allNamesArray[i][0] + opener;
	std::cout<<"hobjN is "<< hobjN<<std::endl;
	hFd[std::make_pair(ch,ds)]->GetObject(hobjN.c_str(),hobj);
	                              hobj->SetDirectory(nullptr);
	if( cms == ds || met == ds )  hobj->SetLineColor( colour);
	else                          hobj->SetFillColor( colour);
	stac .Add(                    hobj);
	}// dataSource
	canv .cd();// pick me to draw?
	stac .Draw("HIST");// must draw before set axes
	stac .GetXaxis()->SetTitle(allNamesArray[i][2].c_str());
	stac .GetYaxis()->SetTitle("Event");
	canv .BuildLegend();
	canv .Update();
	cf   .WriteTObject(&canv);
//	canv .SaveAs(("plots/" + stname + ".root").c_str());// slow
//	canv .SaveAs(("plots/" + stname + ".pdf" ).c_str());
	canv .Clear();// clear rather than delete, reusable!
	// stac will auto clean up since it is not new-ed
	}}// channel, p
	for(auto &x : hFd){
		x.second->Close(); x.second = nullptr;
	}
	hobj = nullptr;
	// file cf will automatically close
	gROOT->SetBatch(kFALSE);// re-enable TBrowser
	return 0;
}
