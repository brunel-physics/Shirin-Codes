//clang++ -Isrc -std=c++17 -march=native -pipe -O3 -Wall -Wextra -Wpedantic -o build/addhists src/Addhists.cxx `root-config --libs` -lm
#include <ROOT/RDataFrame.hxx>//#include <ROOT/RCsvDS.hxx>
#include <Math/Vector4D.h>
#include <TRandom3.h>// used Gaussian, uniform each once
//#include <execution>// need to link -ltbb in Makefile
#include <TChain.h>
#include <TF1.h>

#include "csv.h"
#include "json.hpp"
#include "calchisto.hpp"
#include "eval_complex.hpp"
#include "roccor.Run2.v3/RoccoR.cc"

int debug = 1;

void Addhists(const channel ch){
std::string NPLc, NPLds;
switch (ch){
case elnu:{NPLc = "elnu";break;}
case munu:{NPLc = "munu";break;}
default  :throw std::invalid_argument(
	"Unimplemented ch (NPL file reading)");
}
TFile *tF;
TH2D  *tHdpr, *tHdln, *tHdcm;// for taking from files
TH2D    *Tpr,           *Tln;// totals go here
std::string temp_header = "histo/NPL_",
	    temp_footer =      ".root",
	    temp_opener ;

for(dataSource ds:dataSourceAll){// to go through all ds get the files open

switch (ds){
case  tzq:{NPLds =  "tzq";break;}
case   ww:{NPLds =   "ww";break;}
case   wz:{NPLds =   "wz";break;}
case   zz:{NPLds =   "zz";break;}
case   st:{NPLds =   "st";break;}
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
temp_opener = temp_header+NPLc+"_"+ NPLds+temp_footer;
if(debug >0)std::cout<< "NPLc and NPLds are "<<NPLc<<" "<<NPLds<<std::endl;
if(debug >0)std::cout<<"temp_opener: "<<temp_opener<<std::endl;
//tF = TFile::Open(("histo/NPL_"+ NPLc + "_" + NPLds +".root").c_str());
tF = TFile::Open(temp_opener.c_str());
if(debug > 0) std::cout<<"opening file failed"<<std::endl;
if(!(ds == cms || ds == met)){
	tF ->GetObject(("prompt_LnT_"  + NPLc + "_" + NPLds).c_str(),tHdpr);
	if(debug > 0) std::cout<<"Get Objct 1st failed"<<std::endl;
	tHdpr->SetDirectory(nullptr);// make it stay even if file closed
	tF ->GetObject(("TL_eff_"      + NPLc + "_" + NPLds).c_str(),tHdln);
	if(debug > 0) std::cout<<"Get objct 2nd failed"<<std::endl;
	tHdln->SetDirectory(nullptr);// make it stay even if file closed
	Tpr->Add(tHdpr);
	if(debug > 0) std::cout<<"add 1st objct failed"<<std::endl;
	Tpr->Scale(1/Tpr->Integral());// freqeuncy probability in each bin
	if(debug > 0) std::cout<<"1st renormalization failed"<<std::endl;
	Tln->Add(tHdln);
        if(debug > 0) std::cout<<"add 2nd objct failed"<<std::endl;
	Tln->Scale(1/Tln->Integral());// frequency probability in each bin
        if(debug > 0) std::cout<<"2nd renormalization failed"<<std::endl;
	tF ->Close();
}else{
	tF ->GetObject(("N_data_LnT_"  + NPLc + "_" + NPLds).c_str(),tHdcm);
        if(debug > 0) std::cout<<"Get Objct 3rd failed"<<std::endl;
	tHdcm->SetDirectory(nullptr);// make it stay even if file closed
	tF ->Close();
}}// for and if
// try to associate pointeres correctly and store them
if(debug > 0) std::cout<<"all objct added"<<std::endl;
const TH2D* const prompt_LnT = static_cast<TH2D*>(Tpr  );
if(debug > 0) std::cout<<"moving pointer 1"<<std::endl;
const TH2D* const TL_eff_QCD = static_cast<TH2D*>(Tln  );
if(debug > 0) std::cout<<"moving pointer 2"<<std::endl;
const TH2D* const dt_LnT_cms = static_cast<TH2D*>(tHdcm);
if(debug > 0) std::cout<<"moving pointer 3"<<std::endl;
TFile hf(("aux/NPL/NPL_" + NPLc + ".root").c_str(),"RECREATE");
if(debug > 0) std::cout<<"file created"<<std::endl;
hf.WriteTObject(prompt_LnT  );hf.Flush();sync();
hf.WriteTObject(TL_eff_QCD  );hf.Flush();sync();
hf.WriteTObject(dt_LnT_cms  );hf.Flush();sync();
if(debug > 0) std::cout<<"all objcts stored on file"<<std::endl;
}// void

int main ( int argc , char *argv[] ){
	Addhists(elnu);
	Addhists(munu);
	return 0;
}
