#include <ROOT/RDataFrame.hxx>// big guns
#include <TCanvas.h>
#include <TLegend.h>
#include <THStack.h>

#include "src/tdrstyle.C"
#include "src/calchisto.hpp"

int plotstacks(){
	setTDRStyle();
	TCanvas canv("name to reset","title to reset",10,10,900,900);
	for(std::string particle:{"fin_jets","lep","bjet"}){
	for(PtEtaPhiM k:PtEtaPhiMall){
	std::string   kstring = "_", tkstr = " ", xAxisStr;
	switch(k){
	case pt :{kstring += "_pt";
	          tkstr    = " pT";
	          xAxisStr = " p_{T}/GeV";
	          break;}
	case eta:{kstring += "eta";
	          tkstr   += "eta";
	          xAxisStr = "PseudoRapidity #eta";
	          break;}
	case phi:{kstring += "phi";
	          tkstr   += "phi";
	          xAxisStr = "Azimuthal angle, #phi/rad";
	          break;}
	case  m :{kstring += "mas";
	          tkstr   += "mass";
	          xAxisStr = "mass GeV/#c^{2}";
	          break;}
	}
	for(channel ch:channelAll){
	std::string chN;
	switch   (ch){
		case elnu:{chN ="elnu";break;}
		case munu:{chN ="munu";break;}
	}
	std::string                title = chN + " " + particle;
	if("fin_jets" == particle) title = chN + " jets";
	std::string  stname =(chN+"_"+particle + kstring).c_str() ;
	canv.SetName(stname.c_str());canv.SetTitle(stname.c_str());
	THStack stac(stname.c_str(),(title + tkstr).c_str());
	TH1D   *hobj;
	for(dataSource ds:{ww,ttz}){// TODO: dataSourceAll){
	std::string  opener  = chN + "_";
	int colour;
	switch  (ds){
		case tzq:{opener += "tzq";colour = 6;break;}// magneta
		case  ww:{opener += "_ww";colour = 2;break;}// red
		case  wz:{opener += "_wz";colour = 3;break;}// green
		case  zz:{opener += "_zz";colour = 4;break;}// blue
		case ttz:{opener += "ttz";colour = 5;break;}// yellow
		case met:{opener += "met";colour = 9;break;}// violet
		case cms:{opener += "cms";colour = 1;break;}// cyan
	}
	TFile hf(("histo/" + opener + ".histo").c_str());// source of slow
	std::string  hobjname = (opener+"_"+particle+kstring).c_str();
	hf.GetObject(hobjname.c_str(),hobj);
	hobj->SetDirectory(nullptr);
	//if(ds =="cms"|| ds=="met")hobj->SetLineColor(colour);
	/*else{*/hobj->SetFillColor(colour);/*}*/
	stac .Add(static_cast<TH1D*>(hobj));
	hf   .Close();
	}// dataSource
	canv .cd();// pick me to draw?
	stac .Draw("HIST");// must draw before set axes
	stac .GetXaxis()->SetTitle(xAxisStr.c_str());
	stac .GetYaxis()->SetTitle("Event");
	canv .BuildLegend();
	canv .Update();
	canv .SaveAs(("plots/" + stname + ".root").c_str());// another slow
//	canv .SaveAs(("plots/" + stname + ".pdf" ).c_str());
	canv .Clear();// clear rather than delete, reusable!
	// stac will auto clean up since it is not new-ed
	}}}// channel, particle, component
	return 0;
}
