
// TODO:: PILE UP NEEDS TO BE IMPLEMENTED WITHOUT ANY OTHER SF ONLY FOR MC (ttz)
// TODO:: APPLY LEP JET min Dr as a filter aswell


//clang++ -Isrc -std=c++17 -march=native -pipe -Ofast -Wall -Wextra -Wpedantic -o build/cht src/cht.cxx build/eval_complex.o `root-config --libs` -lm
#include <ROOT/RDataFrame.hxx>//#include <ROOT/RCsvDS.hxx>
#include <TRandom3.h>// used Gaussian, uniform each once

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
  constexpr  float ENDCAP_ETA_MIN = 1.566f;
  constexpr  float BARREL_ETA_MAX = 1.4442f;
//constexpr    int EL_MAX_NUM   = 1;
  constexpr  float EL__PT_MIN   = 30.f;//{15}//min 12, AP 45,
//constexpr  float EL_LPT_MIN   = 35.f;// Leading
  constexpr  float EL_ETA_MAX   = 2.5f;
  constexpr    int EL_LOOSE_ID  = 1;
  constexpr    int EL_TIGHT_ID  = 4;

//constexpr    int MU_MAX_NUM   = 1;
  constexpr  float MU__PT_MIN   = 29.f;//min 33, AP 40,
//constexpr  float MU_LPT_MIN   = 26.f;// Leading
  constexpr  float MU_ETA_MAX   = 2.4f;
  constexpr  float MU_LOOSE_ISO = .15f;
  constexpr  float MU_TIGHT_ISO = .25f;

//constexpr  float MET__PT_MIN  = 40.f;
  constexpr  float MET_EL_PT    = 20.f;//80.f;
  constexpr  float MET_MU_PT    = 25.f;//40.f;

constexpr float    JET_ETA_MAX =  4.7f;
constexpr float    JET__PT_MIN = 30.f;
constexpr double       JET_ISO =   .4;
constexpr unsigned    JETS_MIN = 4;
constexpr unsigned    JETS_MAX = 6;

constexpr double  BJET_ETA_MAX = 2.4;
constexpr double BTAG_DISC_MIN =  .8838;
constexpr unsigned   BJETS_MIN = 1;
constexpr unsigned   BJETS_MAX = 3;

//constexpr double DELTA___R_ZL   = 1.6;
//constexpr double DELTA_PHI_ZW   = 2.;
//constexpr double DELTA_PHI_ZMET = 2.;

constexpr double ak4RconeBy2 =  .2;
constexpr double ak8RconeBy2 =  .4;

// This Pi is more accurate than binary256; good for eternity
template <typename T> constexpr T  PI = T(3.14159265358979323846264338327950288419716939937510582097494459230781640628620899);
template <typename T> constexpr T TPI = PI<T> * 2;

enum      puSf      {puW,upW,dnW};

inline auto triggers(channel ch){
	return [=](
	 const bool el  // HLT_Ele32_WPTight_Gsf_L1DoubleEG
	,const bool enu// HLT_Ele28_eta2p1_WPTight_Gsf_HT150
	,const bool jt  // HLT_PFHT180
	,const bool mu  // HLT_IsoMu27
){
	switch(ch){
//	case elnu:return  el;
//	case elnu:return  enu;
//	case munu:return (mu && jt);
//	case munu:return  mu;
	case elnu:return  jt;
	case munu:return  jt;
	}};
}
auto lep_sel(    const channel ch){
      return [=](const  bools& isPFs,
                 const floats& pts,
                 const floats& etas,
                 const   ints& elids,
                 const  bools& muids,
                 const floats& isos){
		const auto abs_etas = abs(etas);
		switch(ch){
			case elnu:return (isPFs && pts >  EL__PT_MIN
			                 && ((abs_etas <  EL_ETA_MAX
			                 &&   abs_etas >  ENDCAP_ETA_MIN)
			                 ||  (abs_etas <  BARREL_ETA_MAX))
			                 &&      elids >= EL_LOOSE_ID);
			case munu:return (    muids    && isPFs
			                 &&   pts      >  MU__PT_MIN
			                 &&   abs_etas <  MU_ETA_MAX
			                 &&   isos     <= MU_LOOSE_ISO);// TODO: WEIRD
//			default  :throw std::invalid_argument(
//				"Unimplemented ch (lep_sel)");
		}
	};
}
auto lep_tight_cut(const channel ch){
        return [=](const   ints& mask,
                   const   ints& elids,
                   const floats& isos){
// TODO: In the case we want 1 loose lepton, comment the 2 NOTE lines below
		bool result;
		      if(ch==elnu){
			ints   temp = elids[mask];
			result = temp.size() == 1;// Choosing 1 Tight Lepton
			result = result && temp[0] >= EL_TIGHT_ID;// NOTE
		}else if(ch==munu){
			floats temp = isos[mask];
			result = temp.size() == 1;// TODO: WEIRD ISO
			result = result && temp[0] <= MU_TIGHT_ISO;// NOTE
		}else{throw std::invalid_argument(
			"Unimplemented ch (lep_tight_cut)");}
		return result;
	};
}
template <typename T>
[[gnu::const]] auto fastPrincipalRangeReductor(const T diff_phi){
	// This function just reduces input from [-2pi,+2pi] to [-pi,+pi]
	// Domain correctness is the user's responsibility(eg.subtract phis)
	// A more general function is easier: just std::remainder
	if(    diff_phi >  PI<T>) return diff_phi - TPI<T>;
	if(    diff_phi < -PI<T>) return diff_phi + TPI<T>;
	return diff_phi;
}
template <typename T> inline constexpr auto delta_phi(const T dp)
	{return fastPrincipalRangeReductor(dp);}

[[gnu::const]] auto deltaR(
	const double eta1,const double phi1,
	const double eta2,const double phi2)
	{return std::hypot(eta1-eta2,delta_phi(phi1-phi2));}//hypot from geometry

template<typename T,typename U>
[[gnu::const]] bool all_equal(const T& t, const U& u){return t == u;}
template<typename T,typename U,typename... Types>
[[gnu::const]] bool all_equal(const T& t, const U& u, Types const&... args)
	{return t == u && all_equal(u, args...);}
template <typename T>
auto jet_lep_min_deltaR(const    T& jet_etas,
                        const    T& jet_phis,
                        const float lep_eta ,
                        const float lep_phi){
	if(!all_equal(jet_etas.size(),jet_phis.size()))
		throw std::logic_error(
		      "Collections must be the same size (jet-lep dR)");
	if(jet_phis.empty())
		throw std::logic_error(
		      "Collections must not be empty for (jet-lep dR)");
	//if(0<debug) std::cout<<"jet_lep_min_dR"<<std::endl;
	doubles min_dRs; min_dRs.reserve(jet_etas.size());
	std::transform(
		jet_etas.cbegin(),
		jet_etas.  cend(),
		jet_phis.cbegin(),
		std::back_inserter(min_dRs),
		[=](double jet_eta, double jet_phi){
			return deltaR(jet_eta,jet_phi,lep_eta,lep_phi);});
	return min_dRs;
}
inline auto pile(
	 const TH1D* const &PuWd
	,const TH1D* const &PuUd
	,const TH1D* const &PuDd
){
	return[=](const int npv){
	std::map<puSf,double> dict = {{puW,1.},{upW,1.},{dnW,1.}};
	dict[puW] = PuWd->GetBinContent(PuWd->GetXaxis()->FindBin(npv));
	dict[upW] = PuUd->GetBinContent(PuUd->GetXaxis()->FindBin(npv));
	dict[dnW] = PuDd->GetBinContent(PuDd->GetXaxis()->FindBin(npv));
	if(0<debug)std::cout<<"pile "<<npv<<" "
		<<dict[puW]<<" "<<dict[upW]<<" "<<dict[dnW]<<std::endl;
	return dict;};
}
auto event_cleaning(
const bool Flag_goodVertices,
const bool Flag_globalSuperTightHalo2016Filter,
const bool Flag_HBHENoiseFilter,
const bool Flag_HBHENoiseIsoFilter,
const bool Flag_EcalDeadCellTriggerPrimitiveFilter,
const bool Flag_BadPFMuonFilter,
const bool Flag_BadChargedCandidateFilter,
const bool Flag_ecalBadCalibFilter,
const bool Flag_eeBadScFilter
	){
	return
		   Flag_goodVertices
		|| Flag_globalSuperTightHalo2016Filter
		|| Flag_HBHENoiseFilter
		|| Flag_HBHENoiseIsoFilter
		|| Flag_EcalDeadCellTriggerPrimitiveFilter
		|| Flag_BadPFMuonFilter
		|| Flag_BadChargedCandidateFilter
		|| Flag_ecalBadCalibFilter
		|| Flag_eeBadScFilter
	;
}
}// namespace
void cht ( const channel ch , const dataSource ds ){
	std::string temp_header="/data/disk0/nanoAOD_2017/",
	temp_opener,temp_footer="/*.root";/**/
	switch(ds){// tzq and exptData use disk3!
	case ttz:{temp_opener=temp_header+ "ttZToQQ"       +temp_footer;break;}
	case cms:{temp_opener=temp_header+ "ttZToQQ"       +temp_footer;break;}
//	default :throw std::invalid_argument("Unimplemented ds (rdfopen)");
	}// CMS and MET MUST do some OPENABLE file ; reject later
	switch(ch){
		case elnu:{temp_header = "Electron_";break;}
		case munu:{temp_header =     "Muon_";break;}
//		default  :throw std::invalid_argument(
//			"Unimplemented ch (init)");
	}
	ROOT::EnableImplicitMT();
	ROOT::RDataFrame mc__df("Events",temp_opener);// Monte Carlo
	auto df  = mc__df ;
	// make test runs faster by restriction. Real run should not
//	auto dfr = df.Range(10000);// remember to enable MT when NOT range
	auto reco = df// remove one letter to do all
	// lepton selection first
/*	.Filter(triggers(ch),
	{ "HLT_Ele32_WPTight_Gsf_L1DoubleEG"//"HLT_Ele28_eta2p1_WPTight_Gsf_HT150"
	 ,"HLT_Ele28_eta2p1_WPTight_Gsf_HT150"//"HLT_PFMET120_PFMHT120_IDTight"
	 ,"HLT_PFJet15"
	 ,"HLT_IsoMu27"
	},"Triggers Filter")*/
	// lepton selection first
	.Define("loose_leps",lep_sel(ch),
	       {temp_header+"isPFcand",
	        temp_header+"pt" ,
	        temp_header+"eta",
	        "Electron_cutBased",
	        "Muon_tightId",
	        "Muon_pfRelIso04_all"})
	.Filter(lep_tight_cut(ch),{"loose_leps",
	        "Electron_cutBased",// left with 1 tight lepton, no loose
	        "Muon_pfRelIso04_all"},"lepton cut")
	// jets selection follows; tW done and lepton selected
	.Define("jet_lep_min_dR"   , jet_lep_min_deltaR<floats>,
	       {"Jet_eta","Jet_phi",    "lep_eta" ,"lep_phi"})
	.Filter(event_cleaning,{
	        "Flag_goodVertices",
	        "Flag_globalSuperTightHalo2016Filter",
	        "Flag_HBHENoiseFilter",
	        "Flag_HBHENoiseIsoFilter",
	        "Flag_EcalDeadCellTriggerPrimitiveFilter",
	        "Flag_BadPFMuonFilter",
	        "Flag_BadChargedCandidateFilter",
	        "Flag_ecalBadCalibFilter",// TODO: v2?
	        "Flag_eeBadScFilter"
	       },
	        "Event Cleaning filter")
	;
	// now we make the histogram names and titles
	switch(ch){// laugh at muon-neutrino below
		case elnu:{temp_header = "elnu_";
		           temp_footer = "electron-neutrino";break;}
		case munu:{temp_header = "munu_";
		           temp_footer = "muon"  "-neutrino";break;}
//		default  :throw std::invalid_argument(
//			"Unimplemented ch (hist titles)");
	}
	}
}
int main ( int argc , char *argv[] ){
	if ( argc < 2 ) {
		std::cout << "Error: no command provided" << std::endl ;
		return 1 ;
	}
	if ( argc < 3 ) {
		   std::cout
		<< "Error: cht needs channel and data source"
		<< std::endl
		<< "e.g.   cht elnu tzq"
		<< std::endl
		;
		return 2 ;
	}
	channel c ; dataSource d ;
	     if ( const auto chN = std::string_view( argv[2] ) ;
	          "elnu"  == chN ) c = elnu ;
	else if ( "munu"  == chN ) c = munu ;
	else { std::cout << "Error: channel " << chN
		<< " not recognised" << std::endl ;
		return 3 ;
	}
	     if ( const auto dsN = std::string_view( argv[3] ) ;
	else if ( "ttz"  ==  dsN ) d = ttz ;
	else if ( "cms"  ==  dsN ) d = cms ;
	else { std::cout << "Error: data source " << dsN
		<< " not recognised" << std::endl ;
		return 4 ;
	}
	if (cms == d ) {
		   std::cout << "Error: currently only possible to test MC"
		<< std::endl ;
		return 5 ;
	}
		cht(c,d) ;
		return 0 ;
}
