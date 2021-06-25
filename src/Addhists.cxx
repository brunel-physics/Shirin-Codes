//clang++ -Isrc -std=c++17 -march=native -pipe -O3 -Wall -Wextra -Wpedantic -o build/addhists src/Addhists.cxx `root-config --libs` -lm
#include <ROOT/RDataFrame.hxx>//#include <ROOT/RCsvDS.hxx>

#include "calchisto.hpp"

int debug = 1;

namespace{
	std::string
	NPLc = "elnu",
	NPLds = "tzq",
	temp_header = "histo/NPL_",
	temp_footer = ".root",
	temp_opener
	;
	TH2D *hpr, *hln, *dcm;// for taking from files
//	TH2D *cpr, *cln;      // for clones
	TH2D *tpr, *tln;      // totals go here
}

void addhists(const channel ch){
	// do tzq first so that tpr & tln are not null
	if(munu == ch) {NPLc = "munu"; NPLds = "tzq";}
		temp_opener = temp_header + NPLc + "_" + NPLds + temp_footer;
		std::cout << "Opening file " << temp_opener << std::endl;
		TFile zq(temp_opener.c_str());
		if( ! zq.IsOpen()) throw std::runtime_error("File is not opened");
		zq.GetObject(("prompt_LnT_" + NPLc + "_" + NPLds).c_str(),tpr);
		if(!tpr) throw std::runtime_error("prompt not found");
		tpr->SetDirectory(nullptr);// make it stay even if file closed
		zq.GetObject((    "TL_eff_" + NPLc + "_" + NPLds).c_str(),tln);
		if(!tln) throw std::runtime_error("TL eff not found");
		tln->SetDirectory(nullptr);// make it stay even if file closed
		zq.Close();
	NPLds = "cms"; // now do the odd one out
		temp_opener = temp_header + NPLc + "_" + NPLds + temp_footer;
		std::cout << "Opening file " << temp_opener << std::endl;
		TFile dc(temp_opener.c_str());
		if( ! dc.IsOpen()) throw std::runtime_error("File is not opened");
		dc.GetObject(("N_data_LnT_" + NPLc + "_" + NPLds).c_str(),dcm);
		if(!dcm) throw std::runtime_error("N_data_LnT not found");
		dcm->SetDirectory(nullptr);// make it stay even if file closed
		dc.Close();
	for(dataSource ds:dataSourceAll){// to go through all ds get the files open
		if(tzq == ds || met == ds || cms == ds)continue;
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
		temp_opener = temp_header + NPLc + "_" + NPLds + temp_footer;
		if(0<debug) std::cout << "Opening file " << temp_opener << std::endl;
		TFile tf(temp_opener.c_str());
		if( ! tf.IsOpen()) throw std::runtime_error("File is not opened");
		tf.GetObject(("prompt_LnT_" + NPLc + "_" + NPLds).c_str(),hpr);
		if(!hpr) throw std::runtime_error("prompt not found");
		tpr->Add(hpr);
		if(0<debug) std::cout << "added prompt" << std::endl;
		//tpr->Scale(1/tpr->Integral());// frequency probability in each bin
		//if(0<debug) std::cout << "normalised prompt" << std::endl;
//		hpr->SetDirectory(nullptr);// make it stay even if file closed
		tf.GetObject((    "TL_eff_" + NPLc + "_" + NPLds).c_str(),hln);
		if(!hln) throw std::runtime_error("TL eff not found");
		tln->Add(hln);
		if(0<debug) std::cout << "added TL eff" << std::endl;
		tln->Scale(1/tln->Integral());// frequency probability in each bin
		if(0<debug) std::cout << "normalised TL eff" << std::endl;
//		hln->SetDirectory(nullptr);// make it stay even if file closed
		tf.Close();
	}// for
	// try to associate pointers correctly and store them
	if(0<debug) std::cout<<"all objects added"<<std::endl;
	TFile hf(("aux/NPL/NPL_" + NPLc + ".root").c_str(),"RECREATE");
	if(0<debug) std::cout<<"file created"<<std::endl;
	tpr->SetName("prompt_LnT");
	tln->SetName("TL_eff_QCD");
	dcm->SetName("dt_LnT_cms");
	hf.WriteTObject(tpr);hf.Flush();sync();
	hf.WriteTObject(tln);hf.Flush();sync();
	hf.WriteTObject(dcm);hf.Flush();sync();
	if(0<debug) std::cout<<"all objects stored on file"<<std::endl;
}// void

int main ( int argc , char *argv[] ){
	addhists(elnu);
	addhists(munu);
	return 0;
}
