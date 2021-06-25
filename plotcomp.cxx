#include <ROOT/RDataFrame.hxx>// big guns
#include <TCanvas.h>
#include <TLegend.h>
#include <THStack.h>

#include "src/tdrstyle.C"
//#include "src/calchisto.hpp"

enum      channel      {elnu,munu};
constexpr channel
          channelAll[]={elnu,munu};

enum      dataSource	  {tzq,vv,ttz,met,wj,st,tw,cms,ttbar,npl};//,wjt,met,st,stb,stw,stbw,ttl,ttj,ttb,cms};
constexpr dataSource
          dataSourceAll[]={tzq,vv,ttz,met,wj,st,tw,cms,ttbar,npl};//,wjt,met,st,stb,stw,stbw,ttl,ttj,ttb,cms};

enum      PtEtaPhiM	 {pt,eta,phi,m};
constexpr PtEtaPhiM
          PtEtaPhiMall[]={pt,eta,phi,m};


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
	//if (ds == ttz)continue;
	std::string  opener  =  chN ;
	switch  (ds){
		case     tzq:{opener +=  "tzq" ;break;}
		case     vv:{opener +=   "_vv" ;break;}
		case     st:{opener +=   "_ST" ;break;}
		case    ttz:{opener +=   "ttz" ;break;}
		case  ttbar:{opener += "ttbar" ;break;}
                case     wj:{opener +=   "_Wj" ;break;}
		case    met:{opener +=   "met" ;break;}
		case    cms:{opener +=   "cms" ;break;}
                case     tw:{opener +=   "_tW" ;break;}
		case    npl:{opener  ="NPL_run_" + chN;break;}
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
	std::string chN, chF; //chF channel title
	switch     (ch){
		case elnu:  {chN ="elnu";chF =  "#e#nu";break;}
		case munu:  {chN ="munu";chF = "#mu#nu";break;}
	}
	std::string                title = chN + " " + particle;
	if("fin_jets" == particle) title = chN + " jets";
	std::string  stname =(chN+"_"+particle + kstring).c_str() ;
	canv.SetName(stname.c_str());canv.SetTitle(stname.c_str());
	THStack stac(stname.c_str(),(title + tkstr).c_str());
	TLegend legS = TLegend(0.8,0.6,0.95,0.9);
	gStyle->SetOptStat(0);
	int W = 800, H = 600, H_ref = 600, W_ref = 800;
	// references for T, B, L, R
	float T = 0.08 * H_ref, B = 0.12 * H_ref,
	      L = 0.12 * W_ref, R = 0.04 * W_ref;
	double max = -1;
	TPad *pads = new TPad("pad","pad",0.01, 0.315, 0.99, 0.99);//pads = pad stack
	pads->SetTopMargin(0);
	pads->SetFillColor(0);
	pads->SetBorderMode(0);
	pads->SetFrameFillStyle(0);
	pads->SetFrameBorderMode(0);
	pads->SetLeftMargin(L / W);
	pads->SetRightMargin(R / W);
	pads->SetTopMargin(T / H);
	pads->SetBottomMargin(B / H * 0.3);
	pads->SetTickx(0);
	pads->SetTicky(0);
	pads->Draw();
	pads->cd();
	TH1D * rp;
	for(dataSource ds:dataSourceAll){
	//if (ttz == ds)continue;
	std::string  opener  = chN + "_";
	std::string  lgN;
	int colour;
	switch  (ds){
		case    tzq:{opener += "tzq" ;lgN = "tZq"           ;colour =  6 ;break;}// magenta
		case     vv:{opener += "_vv" ;lgN = "VV "           ;colour =  2 ;break;}// red
                case     st:{opener += "_ST" ;lgN = "Single t"      ;colour =  95;break;}//
		case    ttz:{opener += "ttz" ;lgN = "t#bar{t}Z"     ;colour =  5 ;break;}// yellow
		case  ttbar:{opener += "ttb" ;lgN = "t#bar{t}"      ;colour =  7 ;break;}// cyan
                case     wj:{opener += "_Wj" ;lgN = "W/#gamma+Jets" ;colour =  55;break;}//
		case    met:{opener += "met" ;lgN = "MET"           ;colour =  9 ;break;}// violet
		case    cms:{opener += "cms" ;lgN = "data"          ;colour =  1 ;break;}// black
                case     tw:{opener += "_tW" ;lgN = "tW"            ;colour =  75;break;}//
		case    npl:{opener +=  "NPL";lgN = "NPL"           ;colour =  40;break;}
	}
	std::string    hobjN = opener + "_" + particle + kstring ;
	hFd[std::make_pair(ch,ds)]->GetObject(hobjN.c_str(),hobj);
	if(hobj->GetMaximum() > max)  max =   hobj->GetMaximum( );
	                              hobj->SetDirectory(nullptr);
	if( cms == ds || met == ds ){ if(met == ds) 	 continue;
				      stac .Draw("hist");
		                      //stac .SetMaximum(max*1.2);// now plot is set, plot CMS on it
		                      hobj->SetLineColor(  colour);
		                      hobj->SetMarkerColor(colour);
		                      hobj->SetMarkerStyle(20);
		                      hobj->SetMarkerSize(1.0);
		                      legS .AddEntry(hobj,lgN.c_str(),"lep");
		                      hobj->Draw("E0 SAME");
		                      legS .Draw();
		                      rp = (TH1D*)(hobj->Clone());
	}else{                        hobj->SetFillColor( colour);
		                      stac .Add(hobj            );
				      stac .SetMaximum(max*1.2);// now plot is set, plot CMS on it
		                      legS .AddEntry(hobj,lgN.c_str(),"f");
	}// else
	}// dataSource
	canv .cd();// pick me to draw?
	//stac .Draw("HIST");// must draw before set axes
	stac .GetXaxis()->SetTitle(xAxisStr.c_str());
	stac .GetYaxis()->SetTitle("Event");
        rp->Divide((TH1D*)stac.GetStack()->Last());
        TPad *padr = new TPad("pad", "pad", 0.01, 0.01, 0.99, 0.2800);// padr = pad ratio
        padr->SetTopMargin(0);
        padr->SetFillColor(0);
        padr->SetBorderMode(0);
        padr->SetFrameFillStyle(0);
        padr->SetFrameBorderMode(0);
        padr->SetLeftMargin(L / W);
        padr->SetRightMargin(R / W);
        padr->SetTopMargin(T / H);
        padr->SetBottomMargin(B / H * 2.1);
        padr->SetTickx(0);
        padr->SetTicky(0);
        padr->SetGridy(1);
	canv.cd();
        padr->Draw();
        padr->cd();
	rp->Draw();
	//canv .BuildLegend();
	canv .Update();
	cf   .WriteTObject(&canv);
//	canv .SaveAs(("plots/" + stname + ".root").c_str());// slow
//	canv .SaveAs(("plots/" + stname + ".pdf" ).c_str());
	legS .Clear();
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
