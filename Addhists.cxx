clang++ -Isrc -std=c++17 -march=native -pipe -O3 -Wall -Wextra -Wpedantic -o build/Adhis src/Addhists.cxx `root-config --libs` -lm
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

void Addhists(const channel ch,const dataSource ds){
// NEED To Split them accordinly:
std:: string NPLc, NPLds;
switch (ch){
         case elnu:{NPLc = "_elnu";break;}
         case munu:{NPLc = "_munu";break;}
         default  :throw std::invalid_argument(
         "Unimplemented ch (NPL file reading)");
           }
switch (ds){
         case  tzq:{NPLds =  "_tzq";break;}
         case   ww:{NPLds =   "_ww";break;}
         case   wz:{NPLds =   "_wz";break;}
         case   zz:{NPLds =   "_zz";break;}
         case   st:{NPLds =   "_st";break;}
         case  stb:{NPLds =  "_stb";break;}
         case  stw:{NPLds =  "_stw";break;}
         case stbw:{NPLds = "_stbw";break;}
         case wjqq:{NPLds = "_wjqq";break;}
         case wzll:{NPLds = "_wzll";break;}
         case  ttb:{NPLds =  "_ttb";break;}
         case  ttl:{NPLds =  "_ttl";break;}
         case  ttj:{NPLds =  "_ttj";break;}
         case  wjt:{NPLds =  "_wjt";break;}
         case  tz1:{NPLds =  "_tz1";break;}
         case  tz2:{NPLds =  "_tz2";break;}
         case  met:{NPLds =  "_met";break;}
         case  cms:{NPLds =  "_cms";break;}
         default :throw std::invalid_argument(
         "Unimplemented ds (NPL file reading)");
           }
TFile *tF; // We need as many Tfs as ds so maybe also temporary file?
	   // All files needs to be open so they can be added to each other?
	   // can we open one file then next add to each other
	   // save as new close old ones, open next one to add the already saved one
	   // itirate till finish?
TH2D *tHdpr, *tHdln, *tHdcm; // similar to tF
TH2D Tpr, Tln;
for(dataSource ds:dataSourceAll){// to go through all ds get the files open
	tF = TFile::Open(("aux/NPL/NPL"+ NPLc + NPLds +".root").c_str());
if(!(ds == cms || ds == met)){

	tF ->GetObject(("prompt_LnT"   + NPLc + NPLds).c_str(),THdpr);
	tHdpr->SetDirectory(nullptr);// make it stay even if file closed
	tF ->GetObject(("TL_eff"       + NPLc + NPLds).c_str(),THdln);
	tHdln->SetDirectory(nullptr);// make it stay even if file closed
	// need to declare the th and also keep the files open so hists can be added
	// save as new file and close.
	Tpr.Add(THdpr);
	Tln.Add(THdln);
	tF ->Close();
}else{
	tF ->GetObject(("N_data_LnT"   + NPLc + NPLds).c_str(),THdcm);
	tHd->SetDirectory(nullptr);// make it stay even if file closed
	tF ->Close();
}}// for and if
	TFile hf(("histo/NPL_"+NPLc+".root").c_str(),"RECREATE");
	hf.WriteTObject(Tpr            .GetPtr());hf.Flush();sync();
	hf.WriteTObject(Tln            .GetPtr());hf.Flush();sync();
	hf.WriteTObject(THdcm			);hf.Flush();sync();
} // void
// need a main to write
