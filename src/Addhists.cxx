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
TH2D     *Tpr,    *Tln;// totals go here

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
tF = TFile::Open(("aux/NPL/NPL"+ NPLc + NPLds +".root").c_str());

if(!(ds == cms || ds == met)){
	tF ->GetObject(("prompt_LnT_"  + NPLc + "_" + NPLds).c_str(),tHdpr);
	tHdpr->SetDirectory(nullptr);// make it stay even if file closed
	tF ->GetObject(("TL_eff_"      + NPLc + "_" + NPLds).c_str(),tHdln);
	tHdln->SetDirectory(nullptr);// make it stay even if file closed
	Tpr->Add(tHdpr);
	Tpr->Scale(1/Tpr->Integral());// freqeuncy probability in each bin
	Tln->Add(tHdln);
	Tln->Scale(1/Tln->Integral());// frequency probability in each bin
	tF ->Close();
}else{
	tF ->GetObject(("N_data_LnT_"  + NPLc + "_" + NPLds).c_str(),tHdcm);
	tHdcm->SetDirectory(nullptr);// make it stay even if file closed
	tF ->Close();
}}// for and if
if(!(ds == cms || ds == met)){// now that it has gone through each file,
// try to associate pointeres correctly and store them
const TH2D* const prompt_LnT = static_cast<TH2D*>(Tpr  );
const TH2D* const TL_eff_QCD = static_cast<TH2D*>(Tln  );
TFile hf(("histo/NPL_" + NPLc + ".root").c_str(),"RECREATE");
hf.WriteTObject(prompt_LnT  );hf.Flush();sync();
hf.WriteTObject(TL_eff_QCD  );hf.Flush();sync();
}else{
const TH2D* const dt_LnT_cms = static_cast<TH2D*>(tHdcm);
TFile hf(("histo/NPL_" + NPLc + ".root").c_str(),"RECREATE");
hf.WriteTObject(dt_LnT_cms  );hf.Flush();sync();
}// if
}// void

int main ( int argc , char *argv[] ){
	Addhists(elnu);
	Addhists(munu);
	return 0;
}
