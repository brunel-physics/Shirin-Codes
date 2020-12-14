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
#include <iterator>// just for std::size

#include "src/tdrstyle.C"
#include "src/calchisto.hpp"

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
	,{    "zmas_",    "Recon. Z mass ","Z\\ m_{ }\\ (\\text{GeV/}c^{2})"}
	,{        "Z_W_Delta_Phi_",      "Z  W  #Delta#phi ",    "Z &  W  #Delta#phi (rad)"}
	,{      "Z_MET_Delta_Phi_",      "Z MET #Delta#phi ",    "Z & MET #Delta#phi (rad)"}
	,{"Z_pair_jets_Delta_Phi_","Z pair jets #Delta#phi ","Z pair jets #Delta#phi (rad)"}
	,{"ev_w_", "Event Weight ","Weight"}
};

int plotall(){
	gROOT->SetBatch(kTRUE);// no open canvas window
	//setTDRStyle();
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
	std::string chN,chF,lgN; // chf : channel title
	switch     (ch){
		case elnu:  {chN ="elnu";chF =  "e#nu"; break;}
		case munu:  {chN ="munu";chF ="#mu#nu"; break;}
	}// switch
	std::string   title = allNamesArray[i][1] + chF;
	std::string  stname = allNamesArray[i][0] + chN;
	canv.SetName(stname.c_str());canv.SetTitle(stname.c_str());
	THStack stac(stname.c_str(),title.c_str());
	TLegend legS = TLegend(0.8,0.6,0.95,0.9);
	//legS.SetFillStyle(1001);
	//legS.SetBorderSize(1);
	//legS.SetFillColor(0);
	//legS.SetLineColor(0);
	//legS.SetShadowColor(0);
	//legS.SetFillColor(kWhite);
	//legS.SetTextSize(0.02);
	gStyle->SetOptStat(0);
	for(dataSource ds:dataSourceAll){
	std::string  opener  = chN + "_";
	int colour;
	switch  (ds){
		case tzq:{opener += "tzq";lgN = "tZq";colour = 6;break;}// magenta
		case  ww:{opener += "_ww";lgN = "WW ";colour = 2;break;}// red
		case  wz:{opener += "_wz";lgN = "WZ ";colour = 3;break;}// green
		case  zz:{opener += "_zz";lgN = "ZZ ";colour = 4;break;}// blue
		case ttz:{opener += "ttz";lgN = "ttZ";colour = 5;break;}// yellow
		case ttb:{opener += "ttb";lgN = "ttb";colour = 7;break;}// cyan
		case met:{opener += "met";lgN = "MET";colour = 9;break;}// violet
		case cms:{opener += "cms";lgN = "CMS";colour = 1;break;}// black
	}
	std::string    hobjN = allNamesArray[i][0] + opener;
	hFd[std::make_pair(ch,ds)]->GetObject(	     hobjN.c_str(),hobj);
	                              hobj->SetDirectory(	nullptr);
	if( cms == ds || met == ds )  hobj->SetLineColor(	 colour);
	else                          hobj->SetFillColor( 	 colour);
	stac .Add(                    hobj);
	if( cms == ds || met == ds )legS .AddEntry(hobj,lgN.c_str(),"l");
	else			    legS .AddEntry(hobj,lgN.c_str(),"f");
	}// dataSource
	canv .cd();// pick me to draw?
	stac .Draw("HIST");// must draw before set axes
	stac .GetXaxis()->SetTitle(allNamesArray[i][2].c_str());
	stac .GetYaxis()->SetTitle("Event");
	legS .Draw();
	//canv .BuildLegend();
	canv .Update();
	cf   .WriteTObject(&canv);
//	canv .SaveAs(("plots/" + stname + ".root").c_str());// slow
//	canv .SaveAs(("plots/" + stname + ".pdf" ).c_str());
	legS .Clear();
	canv .Clear();// clear rather than delete, reusable!
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

/*
Working legend


 TLegend legend_= TLegend(0.8,0.6,0.95,0.9);
 legend_.SetFillStyle(1001);
 legend_.SetBorderSize(1);
 legend_.SetFillColor(0);
 legend_.SetLineColor(0);
 legend_.SetShadowColor(0);
 legend_.SetFillColor(kWhite);
 legend_.SetTextSize(0.02);
 gStyle->SetOptStat(0);


 std :: map<string,string> IdLegend;

 IdLegend["BNumHist"]="b quark";
 IdLegend["BBarNumHist"]="b bar quark";
for (const auto histId : HistPerCanvas){
   TCanvas *tempcanvas = new TCanvas(histId.first.c_str(),histId.first.c_str(),200,10,700,500);
   tempcanvas->cd();
   int colour = 2;
   float max = -1;
   for (auto histo : HistPerCanvas[histId.first.c_str()]){
     legend_.AddEntry(&IdHist[histo.c_str()],IdLegend[histo.c_str()].c_str(),"l");
     if(colour == 10)
       {// to get rid of the white colour graph in the plots
         colour++;
       }// to get rid of the white colour graph in the plots
     IdHist[histo.c_str()].SetLineColor(colour);
     colour++;
     IdHist[histo.c_str()].Draw("same");
     IdHist[histo.c_str()].GetXaxis()->SetTitle(histId.first.c_str());
     IdHist[histo.c_str()].GetYaxis()->SetTitle("Number of events");
     if (IdHist[histo.c_str()].GetMaximum() > max) max = IdHist[histo.c_str()].GetMaximum();
   }
   IdHist[HistPerCanvas[histId.first.c_str()][0].c_str()].SetMaximum(max*1.2);
   legend_.Draw();
   tempcanvas->SaveAs(("/scratch/eepgssg/plots/plots2017/tZq_lnu_Z_qq_4f/PtEtaECut/"+histId.first+".root").c_str());
   tempcanvas->SaveAs(("/scratch/eepgssg/plots/plots2017/tZq_lnu_Z_qq_4f/PtEtaECut/"+histId.first+".png").c_str());
   legend_.Clear();
}

*/

