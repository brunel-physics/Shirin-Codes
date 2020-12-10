#include <ROOT/RDataFrame.hxx>// big guns
#include <TCanvas.h>
#include <TLegend.h>
#include <THStack.h>

#include "src/tdrstyle.C"
#include "src/calchisto.hpp"

int plotcomp(){
	gROOT->SetBatch(kTRUE);// no open canvas window
	setTDRStyle();
	TFile     cf("plots/components.root","RECREATE");
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
//	if (ttb == ds || met == ds || cms == ds)continue;
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
		= new TFile(("histo/" + opener + ".root").c_str());
	}}// now we have a histogram file dictionary of all the files miahahaha
	for(std::string particle:{"fin_jets","lep","bjet"}){
	for(PtEtaPhiM k:PtEtaPhiMall){
//	if ( e == k ) continue;
	std::string   kstring = "_", tkstr = " ", xAxisStr;
	switch(k){
	case pt :{kstring += "_pt";
	          tkstr    = " p_{T}";
	          xAxisStr = " p_{T}/GeV";
	          break;}
	case eta:{kstring += "eta";
	          tkstr   += "eta";
	          xAxisStr = "PseudoRapidity #eta";
	          break;}
	case phi:{kstring += "phi";
	          tkstr   += "phi";
	          xAxisStr = "Azimuthal angle #varphi/rad";
	          break;}
	case  m :{kstring += "mas";
	          tkstr   += "mass";
	          xAxisStr = "\\text{mass GeV/}c^{2}";
	          break;}
//	case  e :throw std::logic_error("can't plot energy");
	}
	for(channel ch:channelAll){
	std::string chN;
	switch     (ch){
		case elnu:  {chN ="elnu";break;}
		case munu:  {chN ="munu";break;}
	}
	std::string                title = chN + " " + particle;
	if("fin_jets" == particle) title = chN + " jets";
	std::string  stname =(chN+"_"+particle + kstring).c_str() ;
	canv.SetName(stname.c_str());canv.SetTitle(stname.c_str());
	THStack stac(stname.c_str(),(title + tkstr).c_str());
	for(dataSource ds:dataSourceAll){
//	if (ttb == ds || met == ds || cms == ds)continue;
	std::string  opener  = chN + "_";
	int colour;
	switch  (ds){
		case tzq:{opener += "tzq";colour = 6;break;}// magenta
		case  ww:{opener += "_ww";colour = 2;break;}// red
		case  wz:{opener += "_wz";colour = 3;break;}// green
		case  zz:{opener += "_zz";colour = 4;break;}// blue
		case ttb:{opener += "ttb";colour = 7;break;}// cyan
		case ttz:{opener += "ttz";colour = 5;break;}// yellow
		case met:{opener += "met";colour = 9;break;}// violet
		case cms:{opener += "cms";colour = 1;break;}// black
	}
	std::string    hobjN = opener + "_" + particle + kstring ;
	hFd[std::make_pair(ch,ds)]->GetObject(hobjN.c_str(),hobj);
	                              hobj->SetDirectory(nullptr);
	if( cms == ds || met == ds )  hobj->SetLineColor( colour);
	else                          hobj->SetFillColor( colour);
	stac .Add(                    hobj);
	}// dataSource
	canv .cd();// pick me to draw?
	stac .Draw("HIST");// must draw before set axes
	stac .GetXaxis()->SetTitle(xAxisStr.c_str());
	stac .GetYaxis()->SetTitle("Event");
	canv .BuildLegend();
	canv .Update();
	cf   .WriteTObject(&canv);
//	canv .SaveAs(("plots/" + stname + ".root").c_str());// slow
//	canv .SaveAs(("plots/" + stname + ".pdf" ).c_str());
	canv .Clear();// clear rather than delete, reusable!
	// stac will auto clean up since it is not new-ed
	}}}// channel, particle, component
	for(auto &x : hFd){
		x.second->Close(); x.second = nullptr;
	}
	hobj = nullptr;
	// file cf will automatically close
	gROOT->SetBatch(kFALSE);// re-enable TBrowser
	return 0;
}
