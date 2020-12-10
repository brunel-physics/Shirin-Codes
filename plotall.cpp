#include <ROOT/RDataFrame.hxx>// big guns
#include <TCanvas.h>
#include <TLegend.h>
#include <THStack.h>

#include "src/tdrstyle.C"
#include "src/calchisto.hpp"

enum      Plots	 	 {"cmet_sEt","cmt__pt","cmet_dpx","cmet_dpy"
                         ,"met_sEt" ,"met__pt","tTm","tWm","ev_w"
                         ,"zmas","Z_W_Delta_Phi","Z_MET_Delta_Phi"
                         ,"Z_pair_jets_Delta_Phi"};
constexpr Plots
          PlotsAll[]	={"cmet_sEt","cmt__pt","cmet_dpx","cmet_dpy"
                         ,"met_sEt" ,"met__pt","tTm","tWm","ev_w"
                         ,"zmas","Z_W_Delta_Phi","Z_MET_Delta_Phi"
                         ,"Z_pair_jets_Delta_Phi"};
int plotall(){
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
	//for(std::string particle:{"fin_jets","lep","bjet"}){
	/*for(std::string particle:{"cmet_sEt","cmt__pt","cmet_dpx","cmet_dpy"
				 ,"met_sEt" ,"met__pt","tTm","tWm","ev_w"
				 ,"zmas","Z_W_Delta_Phi","Z_MET_Delta_Phi"
				 ,"Z_pair_jets_Delta_Phi"}){*/
	for(Plots p: PlotsAll){
	std::string   kstring = "_", tkstr = " ", xAxisStr;
	switch (p){case		     cmet_sEt:{kstring += "cmet_sEt";
                			       tkstr    = " Sum E_{t}";
                                  	       xAxisStr = " Sum E_{t}/GeV"	       ;break;}
		   case 	     cmet__pt:{kstring += "cmet__pt";
               			  	       tkstr    = " p_{t}";
                		  	       xAxisStr = " p_{t}/GeV"    	       ;break;}
		   case 	     cmet_dpx:{kstring += "cmet_dpx";
                		               tkstr    = " p_{x}";
                		  	       xAxisStr = " p_{x}/GeV"    	       ;break;}
		   case 	     cmet_dpy:{kstring += "cmet_dpy";
                		  	       tkstr    = " p_{y}";
                		  	       xAxisStr = " p_{y}/GeV"    	       ;break;}
		   case 	      met_sEt:{kstring += "met_sEt";
                		  	       tkstr    = " Sum E_{T}";
                                  	       xAxisStr = " Sum E_{T}/GeV"	       ;break;}
		   case 	      met__pt:{kstring += "met__pt";
                		 	       tkstr    = " MET p_{T}";
                		   	       xAxisStr = " MET p_{T}/GeV"	       ;break;}
		   case     	   	  tTm:{kstring += "tTm";
                                               tkstr    = " top m_{T}";
                                               xAxisStr = " top \\text{m_{T} GeV/}c^{2}";break;}
		   case     	          tWm:{kstring += "tWm";
                                               tkstr    = " W m_{t}";
                                               xAxisStr = " W \\text{m_{T} GeV/}c^{2}"  ;break;}
		   case     	         ev_w:{kstring += "ev_w";
                                               tkstr    = " event weight";
                                               xAxisStr = " SFs"		       ;break;}
		   case    	  	 zmas:{kstring += "zmas";
                                               tkstr    = " Z mass";
                                               xAxisStr = " Z \\text{mass GeV/}c^{2}"    ;break;}
		   case   	 Z_W_Deltaphi:{kstring += "Z_W_Deltaphi";
                                               tkstr    = " Z and W \delta #varphi";
                                               xAxisStr = " Z and W \delta #varphi"    ;break;}
		   case        Z_MET_DeltaPhi:{kstring += "Z_MET_Deltaphi";
                                               tkstr    = " Z and MET \delta #varphi/rad";
                                               xAxisStr = " Z and MET \delta #varphi/rad"    ;break;}
		   case Z_pair_jets_Delta_Phi:{kstring += "Z_pair_jets_Delta_Phi";
                                               tkstr    = " Z jet pairs \delta #varphi/rad";
                                               xAxisStr = " Z jet pairs \delta #varphi/rad"    ;break;}
	}// switch

	//}// for loop
	/*for(PtEtaPhiM k:PtEtaPhiMall){
	//if ( e == k ) continue;
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
	}*/
	for(channel ch:channelAll){
	std::string chN;
	switch     (ch){
		case elnu:  {chN ="elnu";break;}
		case munu:  {chN ="munu";break;}
	}
	std::string                title = chN + " " + p;
	//if("fin_jets" == particle) title = chN + " jets";
	std::string  stname =(chN+"_"+p).c_str() ;
	canv.SetName(stname.c_str());canv.SetTitle(stname.c_str());
	THStack stac(stname.c_str(),(chN +" "+ tkstr).c_str());
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
	}}// channel, p
	for(auto &x : hFd){
		x.second->Close(); x.second = nullptr;
	}
	hobj = nullptr;
	// file cf will automatically close
	gROOT->SetBatch(kFALSE);// re-enable TBrowser
	return 0;
}
