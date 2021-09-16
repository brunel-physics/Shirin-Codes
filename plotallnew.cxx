// TODO:: Change title name to be like e.g ev Jet pt -> v , greek nu
// TODO:: legend should be only ds
// TODO:: change the mass unit \\text{W m_{T} (GeV/}c^{2}) doesn't work!
// TODO:: add more histograms for filtred mass and angles.

#include <ROOT/RDataFrame.hxx>// big guns
#include <TCanvas.h>
#include <TLegend.h>
#include <THStack.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TF1.h>
#include <iterator>// just for std::size

#include "src/tdrstyle.C"
//#include "src/calchisto.hpp"


enum      channel      {elnu,munu};
constexpr channel
          channelAll[]={elnu,munu};

enum      dataSource	  {tzq,ttz,wj,vv,st,tw,dy,ttbar,npl/*,wjt,stb,stw,stbw,ttl,ttj,ttb*/,cms};
constexpr dataSource
          dataSourceAll[]={tzq,ttz,wj,vv,st,tw,dy,ttbar,npl/*,wjt,stb,stw,stbw,ttl,ttj,ttb*/,cms};
int debug = 1;

std::string allNamesArray[][3] = {// histogram id, histogram title, x axis string
//	 {"cmet_sEt_",   "Corrected MET " "sum E_{T} ","Sum E_{T} (GeV)"}
	 { "met_sEt_","Un-corrected MET " "sum E_{T} ","Sum E_{T} (GeV)"}
//	,{"cmet__pt_",   "Corrected MET "     "p_{T} ","p_{T} (GeV/c)"}
//	,{"cmet_dpx_",   "Corrected MET #Delta p_{x} ","p_{x} (GeV/c)"}
//	,{"cmet_dpy_",   "Corrected MET #Delta p_{y} ","p_{y} (GeV/c)"}
	,{ "met__pt_","Un-corrected MET "     "p_{T} ","p_{T} (GeV/c)"}
	,{"W_invariant_mass_","W invariant mass ","W\\ m_{ }\\ (\\text{GeV/}c^{2})"}
	,{     "tWm_","Transverse W mass ","W\\ m_{T}\\ (\\text{GeV/}c^{2})"}
	,{     "tTm_","Transverse T mass ","T\\ m_{T}\\ (\\text{GeV/}c^{2})"}
//	,{ "ttop_pt_","Transverse T p_{t}","T\\ p_{T}\\ (Gev)"}
	,{    "zmas_",    "Recon. Z mass ","Z\\ m_{ }\\ (\\text{GeV/}c^{2})"}
	,{        "Z_W_Delta_Phi_",      "Z  W  #Delta#phi ",    "Z &  W  #Delta#phi (rad)"}
	,{      "Z_MET_Delta_Phi_",      "Z MET #Delta#phi ",    "Z & MET #Delta#phi (rad)"}
	,{"Z_pair_jets_Delta_Phi_","Z pair jets #Delta#phi ","Z pair jets #Delta#phi (rad)"}
	,{"ev_w_", "Event Weight ","Weight"}
};

int plotallnew(){
	gROOT->SetBatch(kTRUE);// no open canvas window
//	setTDRStyle();
	TFile     cf("plots/plotall.root","RECREATE");
	TCanvas canv("name to reset","title to reset",10,10,900,900);
	TH1D   *hobj;
	// now we open ALL the files
	std::map<std::pair<channel,dataSource>,TFile*> hFd;
	for(channel ch:channelAll){
	std::string chN;
	switch     (ch){
		case elnu:  {chN ="elnu_";break;}
		case munu:  {chN ="munu_";break;}
	}
	for(dataSource ds:dataSourceAll){
	std::string  opener  =  chN ;
	switch  (ds){
		case    tzq:{opener +=  "tzq";break;}
		case     vv:{opener +=  "_vv";break;}
		case     st:{opener +=  "_ST";break;}
		case    ttz:{opener +=  "ttz";break;}
		case  ttbar:{opener +="ttbar";break;}
		case     wj:{opener +=  "_Wj";break;}
		case     dy:{opener +=  "_DY";break;}
		//case    met:{opener +=  "met";break;}
		case    cms:{opener +=  "cms";break;}
		case     tw:{opener +=  "_tW";break;}
		case    npl:{opener  = "NPL_run_" + chN;break;}
//		case    npl:{opener  ="NPL_run_" + chN;break;}
	}
	if(debug > 0) std::cout<<ds<<std::endl;
	hFd[std::make_pair(ch,ds)]
		= new TFile(("histo/" + opener + ".root").c_str());
	}}// now we have a histogram file dictionary of all the files miahahaha

	for(size_t i=0; i < std::size(allNamesArray) ;++i){
	for(channel ch:channelAll){
	std::string chN,chF,lgN; // chf : channel title
	switch     (ch){
		case elnu:  {chN ="elnu";chF =  "e#nu"; break;}
		case munu:  {chN ="munu";chF ="#mu#nu"; break;}
	}// switch
	double max = -1;
	std::string   title = allNamesArray[i][1] + chF;
	std::string  stname = allNamesArray[i][0] + chN;
	canv.SetName(stname.c_str());canv.SetTitle(stname.c_str());
	THStack stac(stname.c_str(),title.c_str());
//	THStack stad(stname.c_str(),title.c_str());//testing stack for data
	TLegend legS = TLegend(0.8,0.6,0.95,0.9);
/*	legS.SetFillStyle(1001);
	legS.SetBorderSize(1);
	legS.SetFillColor(0);
	legS.SetLineColor(0);
	legS.SetShadowColor(0);
	legS.SetFillColor(kWhite);
	legS.SetTextSize(0.02);
*/	gStyle->SetOptStat(0);
	for(dataSource ds:dataSourceAll){
	//if(met == ds) continue;
	std::string  opener  = chN + "_";
	int colour;
	switch  (ds){
		case   tzq:{opener+="tzq";lgN="tZq"          ;colour= 6;break;}// magenta
		case    vv:{opener+="_vv";lgN="VV "          ;colour= 2;break;}// red
		case    st:{opener+="_ST";lgN="Single t"     ;colour=95;break;}//
		case   ttz:{opener+="ttz";lgN="t#bar{t}Z"    ;colour= 5;break;}// yellow
		case ttbar:{opener+="ttb";lgN="t#bar{t}"     ;colour= 7;break;}// cyan
		case    wj:{opener+="_Wj";lgN="W+Jets"       ;colour=55;break;}//
		case    dy:{opener+="_DY";lgN="Z/#gamma+Jets";colour=79;break;}//
		//case   met:{opener+="met";lgN="MET"          ;colour= 9;break;}// violet
		case   cms:{opener+="cms";lgN="data"         ;colour= 1;break;}// black
		case    tw:{opener+="_tW";lgN="tW"           ;colour=75;break;}//
		case   npl:{opener+="NPL";lgN="NPL"          ;colour=40;break;}//
	}
	std::string hobjN = allNamesArray[i][0] + opener;
	hFd[std::make_pair(ch,ds)]->GetObject(hobjN.c_str(),        hobj);
//	std::cout<< "ds and ch "<<ds<< " "<<ch<<std::endl;
	if( 0  < debug)std::cout<<"passed make pair "<<hobjN<<std::endl;
	                                      hobj->SetDirectory(nullptr);
	if( 0  < debug)std::cout<<"passed setDir"<<std::endl;
	if(max < hobj->GetMaximum()) max =    hobj->GetMaximum(         );
	if( 0  < debug)std::cout<<"passed get max"<<std::endl;
	// note that MET is already skipped above
	// WARNING: We require CMS to be the last thing in dataSourceAll !
	if( cms != ds ){
		hobj->Scale(1./hobj->Integral("width"));
		hobj->SetFillColor(colour);
		stac .Add(hobj);
		legS .AddEntry(hobj,lgN.c_str(),"f");
	}else{// CMS is last, so we plot at this time!
//		canv .cd();// pick me to draw?
		TH1D *lsht = static_cast<TH1D*>(stac.GetStack()->Last());
//		lsht->Scale(1./ls->Integral("width"));
		stac .Draw("hist");// TODO: histe
		stac .GetXaxis()->SetTitle(allNamesArray[i][2].c_str());
		stac .GetYaxis()->SetTitle("Event");
		stac .SetMaximum(max*1.2);
		// now plot is set, plot CMS on it
		hobj->SetLineColor(  colour);
		hobj->SetMarkerColor(colour);
		hobj->SetMarkerStyle(20);
		hobj->SetMarkerSize(1.0);
		legS .AddEntry(hobj,lgN.c_str(),"lep");
		hobj->Draw("e0 same");
		lsht->Draw("same");
		legS .Draw();
	}}// else & dataSource
//	canv.cd();
//	canv.BuildLegend();
	canv.Update();
	cf  .WriteTObject(&canv);
	legS.Clear();
	canv.Clear();// clear rather than delete, reusable!
	// stac will auto clean up since it is not new-ed
	}}// channel, i
	for(auto &x : hFd){
		x.second->Close(); x.second = nullptr;
	}
	hobj = nullptr;
	// file cf will automatically close
	gROOT->SetBatch(kFALSE);// re-enable TBrowser
	return 0;
}
