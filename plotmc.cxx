#include <ROOT/RDataFrame.hxx>// big guns
#include <TCanvas.h>
#include <TLegend.h>
#include <THStack.h>
#include <TROOT.h>
#include <TStyle.h>
#include "src/tdrstyle.C"
#include "src/calchisto.hpp"

int plotmc(){
	gROOT->SetBatch(kTRUE);// no open canvas window
	setTDRStyle();
	TFile     cf("plots/MC.root","RECREATE");
	TCanvas canv("name to reset","title to reset",10,10,900,900);
	TH1D   *hobj;
	// now we open ALL the files
	std::map<std::pair<channel,dataSource>,TFile*> hFd;
	// Channel and dataSource taken, to map to a TFile pointer
	for(channel ch:channelAll){
	std::string chN;
	switch     (ch){
		case elnu:  {chN ="elnu_";break;}
		case munu:  {chN ="munu_";break;}
	}
	for(dataSource ds:dataSourceAll){
	if(cms == ds || met == ds)continue;
	std::string  opener  =  chN ;
	switch  (ds){
		case tzq:{opener += "tzq";break;}
		case  ww:{opener += "_ww";break;}
		case  wz:{opener += "_wz";break;}
		case  zz:{opener += "_zz";break;}
		case ttb:{opener += "ttb";break;}
		case ttz:{opener += "ttz";break;}
		case met:{opener += "met";break;}
		case cms:{opener += "cms";break;}
	}
	hFd[std::make_pair(ch,ds)]
		= new TFile(("histo/Result/" + opener + ".root").c_str());
	}}// now we have a histogram file dictionary of all the files miahahaha
	for(std::string sf:
	{"sfi","sfj","p_ej","p_sfej","p_sf_i","btag_w","cmet_sEt","cmet__pt","cmet_dpx","cmet_dpy"}){
	std::string xAxisStr, staT; //staT = stack title
	if(false);
	else if(  "sfi"    == sf) {staT = "btag SF_{i}";xAxisStr = "sf_{i}";}
	else if(  "sfj"    == sf) {staT = "btag SF_{j}";xAxisStr = "sf_{j}";}
	else if("p_ej"     == sf) {staT = "\\prod_{j} - e_{j}";xAxisStr = "\\prod_{j} 1 -                e_{j}";}
	else if("p_sfej"   == sf) {staT = "\\prod_{j} 1 - \\text{sf}_{j} e_{j}";xAxisStr = "\\prod_{j} 1 - \\text{sf}_{j} e_{j}";}
	else if("p_sf_i"   == sf) {staT = "\\prod_{i}     \\text{sf}_{i}";xAxisStr = "\\prod_{i}     \\text{sf}_{i}      ";}
	else if("btag_w"   == sf) {staT = "btag weight";xAxisStr = "btag weight";}
	else if("cmet_sEt" == sf) {staT = "corrected MET Sum E_{T}";xAxisStr = "Sum E_{T} (GeV)";}
	else if("cmet__pt" == sf) {staT = "corrected MET p_{T}";xAxisStr = "p_{T} (GeV/c)";}
	else if("cmet_dpx" == sf) {staT = "corrected dp_{x}";xAxisStr = "p_{x} (GeV/c)";}
	else if("cmet_dpy" == sf) {staT = "corrected dp_{y}";xAxisStr = "p_{y} (GeV/c)";}

	for(channel ch:channelAll){
	std::string chN,chF,lgN; // chF channel title,lgN legend name
	switch     (ch){
		case elnu:  {chN ="elnu";chF =  "#e#nu";break;}
		case munu:  {chN ="munu";chF = "#mu#nu";break;}
	}
	std::string                title = staT + " " + chF;
	std::string  stname =(sf + "_" + chN).c_str() ;
	canv.SetName(stname.c_str());canv.SetTitle(stname.c_str());
	THStack stac(stname.c_str(),title.c_str());
	TLegend legS = TLegend(0.8,0.6,0.95,0.9);
	gStyle->SetOptStat(0);
	for(dataSource ds:dataSourceAll){
	std::string  opener  = chN + "_";
	if(cms == ds || met == ds)continue;// NOTE: Skip since NOT Monte Carlo
	int colour;
	switch  (ds){
		case tzq:{opener += "tzq";lgN ="tzq";colour = 6;break;}// magenta
		case  ww:{opener += "_ww";lgN = "WW";colour = 2;break;}// red
		case  wz:{opener += "_wz";lgN = "WZ";colour = 3;break;}// green
		case  zz:{opener += "_zz";lgN = "ZZ";colour = 4;break;}// blue
		case ttb:{opener += "ttb";lgN ="ttb";colour = 7;break;}// cyan
		case ttz:{opener += "ttz";lgN ="ttZ";colour = 5;break;}// yellow
		case met:{opener += "met";lgN ="MET";colour = 9;break;}// violet
		case cms:{opener += "cms";lgN ="CMS";colour = 1;break;}// black
	}
	std::string    hobjN = sf + "_" + opener ;
	hFd[std::make_pair(ch,ds)]->GetObject(hobjN.c_str(),hobj);
	                              hobj->SetDirectory(nullptr);
	if( cms == ds || met == ds )  hobj->SetLineColor( colour);
	else                          hobj->SetFillColor( colour);
	stac .Add(                    hobj);
	if( cms == ds || met == ds )legS .AddEntry(hobj,lgN.c_str(),"l");
	else			    legS .AddEntry(hobj,lgN.c_str(),"f");
	}// dataSource
	canv .cd();// pick me to draw?
	stac .Draw("HIST");// must draw before set axes
	stac .GetXaxis()->SetTitle(xAxisStr.c_str());
	stac .GetYaxis()->SetTitle("Event");
	legS.Draw();
	//canv .BuildLegend();
	canv .Update();
	cf   .WriteTObject(&canv);
//	canv .SaveAs(("plots/" + stname + ".root").c_str());// slow
//	canv .SaveAs(("plots/" + stname + ".pdf" ).c_str());
	legS.Clear();
	canv .Clear();// clear rather than delete, reusable!
	// stac will auto clean up since it is not new-ed
	}}// channel, particle
	for(auto &x : hFd){
		x.second->Close(); x.second = nullptr;
	}
	hobj = nullptr;
	// file cf will automatically close
	gROOT->SetBatch(kFALSE);// re-enable TBrowser
	return 0;
}
