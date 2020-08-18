// Compilor:
//clang++ -Isrc -std=c++17 -march=native -pipe -O3 -Wall -Wextra -Wpedantic -o build/dptest src/TestDeep.cxx `root-config --libs` -lm


#include <ROOT/RDataFrame.hxx>//#include <ROOT/RCsvDS.hxx>
#include <TChain.h>
#include "csv.h"
#include "calchisto.hpp"
#include "eval_complex.hpp"

using doubles = ROOT::VecOps::RVec<double>;
using  floats = ROOT::VecOps::RVec<float>;
using    ints = ROOT::VecOps::RVec<int>;
using   bools = ROOT::VecOps::RVec<bool>;
using strings = ROOT::VecOps::RVec<std::string>;


namespace{
  constexpr    int debug = 0;
  constexpr double DBTAG_DISC_MIN =  .4941; // DeepBTagCSV

template<typename T,typename U>
[[gnu::const]] bool all_equal(const T& t, const U& u){return t == u;}
template<typename T,typename U,typename... Types>
[[gnu::const]] bool all_equal(const T& t, const U& u, Types const&... args)
	{return t == u && all_equal(u, args...);}

inline auto btagCSVv2(const bool check_CSVv2){
	return [=](
		 const  floats& btag
		,const doubles& pt
		,const doubles& eta
		,const    ints& flav
	){
	if(0<debug)std::cout<<"btagCSVv2 entered"<<std::endl;
	const size_t     size = pt  .size();
	strings formulae(size  ,"1");//check_CSVv2 ? "1" : "0");
	doubles  results(size );
	if(!all_equal(   size  ,flav.size(),
	             eta.size(),btag.size())) throw std::logic_error(
		"Collections must be the same size in btagCSVv2");
	if(0 == size) throw std::logic_error(
		"Collections must not be empty in btagCSVv2");
	std::string  measureType,sysType,rawFormula;
	ints   aflv   = abs(flav);
	bool   b;
	int    jetFlav;// TODO: jet flav 5 or 0? BTag Website, B = 0;
	double CSVv2  ,
	       pt_Min , pt_Max,
	       etaMin , etaMax,
	       CSVmin , CSVmax;
	//io::CSVReader<11> thisCSVfile("aux/CSVv2_94XSF_V2_B_F.csv");
	io::CSVReader<11> thisCSVfile("aux/DeepCSVv2_94XSF_V5_B_F.csv");
	thisCSVfile.next_line();// we happen to not need the header line
	// The following nests too much, so we do not indent
	// Each blank line means nesting deeper
	while(thisCSVfile.read_row(CSVv2,measureType,sysType,jetFlav,
	      etaMin,etaMax,pt_Min,pt_Max,CSVmin,CSVmax,rawFormula)){
	// CSVv2 column = Operating point
	if(check_CSVv2){
		b= DBTAG_DISC_MIN <= CSVv2
		&& "mujets" == measureType && 0 == jetFlav;
	}else{
		b=   "incl" == measureType && 0 != jetFlav;
	}
	if(b && "central" ==  sysType ){

	for(size_t i=0; i < pt.size() ;++i){
	if(check_CSVv2){
		b= CSVmin < btag[i] && btag[i] < CSVmax
		&& aflv[i] == 5;
//		&& aflv.at(i) == 5;
	}else{
		b= aflv[i] <  5 || 21 == aflv[i] ;
//		b= aflv.at(i) <  5 || 21 == aflv.at(i) ;
	}
	std::string tempFormula = rawFormula;
	if(b
	&& etaMin < eta [i] && eta [i] < etaMax
	&& pt_Min < pt  [i] && pt  [i] < pt_Max){
	std::cout<< "pt and eta are " <<
		 <<pt  [i] << eta [i] <<std::endl;

	if(1 == formulae[i].length()){// only 1st found wins

	if(std::string::npos != tempFormula.find("x")){

	// now guaranteed one good formula with x in it
	// Time to replace all x with pt ourselves~~ No need boost::replace_all
	std::string  ptstr = std::to_string(pt[i]);
	size_t pos  = tempFormula.find("x");
	while( pos != std::string::npos){
	              tempFormula.replace( pos,1,ptstr);// 1 = "x".size()
	       pos  = tempFormula.find("x",pos + ptstr.size());
	}
	std::cout<<"the formula is "<<tempFormula<<std::endl;
	formulae[i] = tempFormula;
	}}}}}
	}// No need to close file after this while loop.
	// resume indentation
	for(size_t j=0; j < formulae.size() ;++j){
		Eval ev;// numbers calculator
		results[j] = ev.eval(const_cast<char*>(
		             formulae[j].c_str())).real();
	}
	std::cout<<"result "<<results<<std::end;
	if(0<debug)std::cout<<"btagCSVv2 exiting"<<std::endl;
	return results;};
}

}// namespace
void TestDeep(const channel ch,const dataSource ds){

	std::string temp_header="/data/disk0/nanoAOD_2017/",
	temp_opener,temp_footer="/*.root";/**/
	switch(ds){// CMS and MET MUST do some OPENABLE file ; reject later
	case tzq:{temp_opener="/data/disk3/nanoAOD_2017/tZqlvqq/*.root"  ;break;}/**/
	}
	ROOT::RDataFrame mc__df("Events",temp_opener);// Monte Carlo
	auto df = [&,ch,ds](){// Get correct data frame
		switch(ds){
			case tzq:{     return mc__df;break;}
		}();
	switch(ch){
		case elnu:{temp_header = "Electron_";break;}
		case munu:{temp_header =     "Muon_";break;}
//		default  :throw std::invalid_argument(
//			"Unimplemented ch (init)");
	}
	// make test runs faster by restriction. Real run should not
	//auto dfr = df.Range(1000);// remember to enable MT when NOT range
	auto init_selection = df// remove one letter to do all
	.Define("sfi",btagCSVv2( true),//,btagDF),// checks btag
	       { "Jet_btagDeepB","Jet_pt","Jet_eta","Jet_partonFlavour"})
	;
	switch(ch){// laugh at muon-neutrino below
		case elnu:{temp_header = "elnu_";
		           temp_footer = "electron-neutrino";break;}
		case munu:{temp_header = "munu_";
		           temp_footer = "muon"  "-neutrino";break;}
//		default  :throw std::invalid_argument(
//			"Unimplemented ch (hist titles)");
	}
	temp_footer = "pt vs eta in " + temp_footer + " channel for ";
	switch(ds){
		case tzq:{temp_header+="tzq";temp_footer+="tZq";break;}
	}
	auto h_sfi = init_selection.Histo1D
			({("sfi_"+temp_header).c_str(),
			  ("sfi "+temp_header).c_str(),50,0,2},"sfi");
	h_sfi->GetXaxis()->SetTitle("sf_{i}");
	h_sfi->GetYaxis()->SetTitle("Event");
	h_sfi->SetLineStyle(kSolid);

	TFile hf(("histo/"+temp_header+".root").c_str(),"RECREATE");
		hf.WriteTObject(h_sfi                  .GetPtr());hf.Flush();sync();
	std::cout<<"TestDeep is finished"<<std::endl;
}// void
int main ( int argc , char *argv[] ){
	if ( argc < 2 ) {
		std::cout << "Error: no command provided" << std::endl ;
		return 1 ;
	}
	if ( argc < 3 ) {
		   std::cout
		<< "Error: tsf needs channel and data source"
		<< std::endl
		<< "e.g.   tsf elnu DY"
		<< std::endl
		;
		return 2 ;
	}
	channel c ; dataSource d ;
	     if ( const auto chN = std::string_view( argv[1] ) ; false ) ;
	else if ( "elnu"  == chN ) c = elnu ;
	else if ( "munu"  == chN ) c = munu ;
	else { std::cout << "Error: channel " << chN
		<< " not recognised" << std::endl ;
		return 3 ;
	}
	     if ( const auto dsN = std::string_view( argv[2] ) ; false ) ;
	else if ( "tzq"  ==  dsN ){ d = tzq  ;}
		else { std::cout << "Error: data source " << dsN
		<< " not recognised" << std::endl ;
		return 4 ;
	}
		TestDeep(c,d) ;
		return 0 ;
}
