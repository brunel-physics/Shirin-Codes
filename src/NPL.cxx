// This script is designed to compute the NPL estimation of the data. The main NPL calculation is in calchisto.cpp

// Compiler:
//clang++ -Isrc -std=c++17 -march=native -pipe -O3 -Wall -Wextra -Wpedantic -o build/npl src/NPL.cxx `root-config --libs` -lm

/* Methodology:
// Non-prompt-lepton estimation
inline auto Npl(
         const int Nqcd_data // Number of fake lepton in QCD region in data
        ,const int Nqcd_rmid // Number of real lepton + mis-id in MC
        ,const int N_sig__MC // Number of lepton in signal region MC
        ,const int N_qcd__MC // Number of lepton in  QCD region MC
){// Need to use the Count in RDF to get these constants after blinding is done
  // Blinding will allow us to know the sidebands which basically are responsible
  // for the QCD regions so by using the report for the filter and the count
  // all the Ns will be found..
  // all Nqcd need to come from the loosend isolation
  // loose isolation 20 GeV for absolute isolation
  // 0.8 for relative isolation
  // impact parameters |dxy| and |dz| loosened to 0.1, 0.5
  // Nqcd_rmid = the qcd with GenPart_Statusflag = 0 (prompt)
        return    (Nqcd_data - Nqcd_rmid)
                * (N_sig__MC / N_qcd__MC);
}

Therefore we need to compute the Nqcd_data which are the loose (excluding the tight criteria) leptons which are the remanant of the blinding procedure.
*/

#include <ROOT/RDataFrame.hxx>//#include <ROOT/RCsvDS.hxx>
#include <Math/Vector4D.h>
#include <TRandom3.h>// used Gaussian, uniform each once
//#include <execution>// need to link -ltbb in Makefile
#include <TChain.h>
#include <TF1.h>
#include "calchisto.hpp"
#include "csv.h"
#include "json.hpp"
#include "eval_complex.hpp"
#include "roccor.Run2.v3/RoccoR.cc"
/*
#if !defined(__FMA__) && defined(__AVX2__)
    #define __FMA__ 1
#endif
*/
using doubles = ROOT::VecOps::RVec<double>;
using  floats = ROOT::VecOps::RVec<float>;
using    ints = ROOT::VecOps::RVec<int>;
using   bools = ROOT::VecOps::RVec<bool>;
using strings = ROOT::VecOps::RVec<std::string>;
namespace{
  constexpr    int debug = 0;
//constexpr    int EL_MAX_NUM      = 1       ;
  constexpr  float EL__PT_MIN      = 35.f    ;
  constexpr  float EL_ETA_MAX      = 2.5f    ;
  constexpr    int EL_LOOSE_ID     = 1       ;
  constexpr    int EL_TIGHT_ID     = 4       ;
  constexpr  float ENDCAP_DZ__TIGHT=  .10f   ;
  constexpr  float ENDCAP_DXY_TIGHT=  .02f   ;
  constexpr  float ENDCAP_DZ__LOOSE=  .20f   ;
  constexpr  float ENDCAP_DXY_LOOSE=  .02f   ;
  constexpr  float ENDCAP_ETA_MIN  = 1.5660f ;
  constexpr  float BARREL_ETA_MAX  = 1.4442f ;
  constexpr  float BARREL_DXY_TIGHT=  .02f   ;
  constexpr  float BARREL_DZ__TIGHT=  .10f   ;
  constexpr  float BARREL_DXY_LOOSE=  .02f   ;
  constexpr  float BARREL_DZ__LOOSE=  .20f   ;

//constexpr    int   MU_MAX_NUM   = 1   ;
  constexpr  float   MU__PT_MIN   = 29.f;
  constexpr  float   MU_ETA_MAX   = 2.4f;
  constexpr  float   MU_LOOSE_ISO = .25f;
  constexpr  float   MU_TIGHT_ISO = .15f;

//constexpr  float    MET__PT_MIN = 40.f;
  constexpr  float    MET_EL_PT   = 20.f;//80.f;// TODO: Need new values
  constexpr  float    MET_MU_PT   = 25.f;//40.f;

  constexpr double     Z_MASS     =  91.1876;
  constexpr double     Z_MASS_CUT =  20.    ;
  constexpr double     W_MASS     =  80.385 ;
  constexpr double     W_MASS_CUT =  80.385 ;
  constexpr double   TOP_MASS     = 172.5   ;
//constexpr double   TOP_MASS_CUT =  20.    ;

  constexpr float     JET_ETA_MAX =  4.7f;
  constexpr float     JET_PT__MIN = 30.0f;
  constexpr double        JET_ISO =   .4 ;
  constexpr unsigned     JETS_MIN =  4   ;
  constexpr unsigned     JETS_MAX =  6   ;

  constexpr double   BJET_ETA_MAX = 2.4   ;
  constexpr double  BTAG_DISC_MIN =  .8838;
//constexpr double DBTAG_DISC_MIN =  .4941; // DeepBTagCSV
  constexpr unsigned    BJETS_MIN = 1     ;
  constexpr unsigned    BJETS_MAX = 3     ;

//constexpr double DELTA___R_ZL   = 1.6;
//constexpr double DELTA_PHI_ZW   = 2. ;
//constexpr double DELTA_PHI_ZMET = 2. ;

  constexpr double    ak4RconeBy2 =  .2;
//constexpr double    ak8RconeBy2 =  .4;

template <typename T> constexpr T  PI = T(3.14159265358979323846264338327950288419716939937510582097494459230781640628620899);
template <typename T> constexpr T TPI = PI<T> * 2;


inline auto lep_sel(const channel ch){
	return [=](
		 const  bools& isPFs
		,const floats& pts
		,const floats& etas
		,const   ints& elids
		,const  bools& muids
		,const floats& isos
		,const floats& dxy
		,const floats& dz
	){
		const auto abs_etas = abs(etas);
		switch(ch){
			case elnu:return ( true
				&&     isPFs
				&&       pts  >   EL__PT_MIN
				&& ((abs_etas <   EL_ETA_MAX
				&&        dz  <=  ENDCAP_DZ__TIGHT
				&&        dxy <=  ENDCAP_DXY_TIGHT
				&&   abs_etas >   ENDCAP_ETA_MIN)
				||  (abs_etas <   BARREL_ETA_MAX
				&&        dxy <=  BARREL_DXY_TIGHT
				&&        dz  <=  BARREL_DZ__TIGHT))
				&&      elids >=  EL_LOOSE_ID
			);
			case munu:return ( true
				&&   muids
				&&   isPFs
				&&     pts  >     MU__PT_MIN
				&& abs_etas <     MU_ETA_MAX
				&&     isos <=    MU_LOOSE_ISO
			);
			default  :throw std::invalid_argument(
				"Unimplemented ch (lep_sel)");
		}
	};
}
inline auto not_tight(const channel ch){
        return [=](
                 const  bools& isPFs
                ,const floats& pts
                ,const floats& etas
                ,const   ints& elids
                ,const  bools& muids
                ,const floats& isos
                ,const floats& dxy
                ,const floats& dz
        ){
                const auto abs_etas = abs(etas);
                switch(ch){
                        case elnu:return ( true
                                &&     isPFs
                                &&	 pts  >   EL__PT_MIN
                                && ((abs_etas <   EL_ETA_MAX
                                &&        dz  <=  ENDCAP_DZ__LOOSE
                                &&        dxy <=  ENDCAP_DXY_LOOSE
                                &&   abs_etas >   ENDCAP_ETA_MIN)
                                ||  (abs_etas <   BARREL_ETA_MAX
                                &&        dxy <=  BARREL_DXY_LOOSE
                                &&        dz  <=  BARREL_DZ__LOOSE))
                                &&	elids <=  EL_LOOSE_ID // Make it not tightable
                        );
                        case munu:return ( true
                                &&   muids
                                &&   isPFs
                                &&     pts  >  MU__PT_MIN
                                && abs_etas <  MU_ETA_MAX
                                &&     isos >= MU_LOOSE_ISO // Make it not tightable
                        );
                        default  :throw std::invalid_argument(
                                "Unimplemented ch (lep_sel)");
                }
        };
}
inline auto lep_tight_cut(const channel ch){
	return [=](const   ints& mask
		  ,const   ints& elids
		  ,const floats& isos
	){
// TODO: In the case we want 1 loose lepton, comment the 4 NOTE lines below
		bool result;
		 if(false) ;
		 else if(ch==elnu){
			ints   temp = elids[mask];
			result = temp.size() == 1;// Choosing 1 Tight Lepton
			result = result &&    temp [0] >= EL_TIGHT_ID ;// NOTE
		}else if(ch==munu){
			floats temp = isos[mask];
			result = temp.size() == 1;
			result = result &&    temp [0] <= MU_TIGHT_ISO;// NOTE
		}else{throw std::invalid_argument(
			"Unimplemented ch (lep_tight_cut)");}
		return result;
	};
}
inline auto not_tight_cut(const channel ch){
        return [=](const   ints& mask
                  ,const   ints& elids
                  ,const floats& isos
        ){
// TODO: In the case we want 1 loose lepton, comment the 4 NOTE lines below
                bool result;
                 if(false) ;
                 else if(ch==elnu){
                        ints   temp = elids[mask];
                        result = temp.size() == 1;// Choosing 1 Tight Lepton
                }else if(ch==munu){
                        floats temp = isos[mask];
                        result = temp.size() == 1;
                }else{throw std::invalid_argument(
                        "Unimplemented ch (not_tight_cut)");}
                return result;
        };
}
inline auto Triggers(const channel ch){
        return [=](
	 const bool el// HLT_Ele32_WPTight_Gsf_L1DoubleEG
        ,const bool mu// HLT_IsoMu27
	,const bool ht// HLT_PFHT250
){
        switch(ch){
        case elnu:return el && ht;
        case munu:return mu && ht;
        }};
}

template <typename T>
constexpr auto fastPrincipalRangeReductor(const T diff_phi){
	// This function just reduces input from [-2pi,+2pi] to [-pi,+pi]
	// Domain correctness is the user's responsibility(eg.subtract phis)
	// A more general function is easier: just std::remainder
	if(    diff_phi >  PI<T>) return diff_phi - TPI<T>;
	if(    diff_phi < -PI<T>) return diff_phi + TPI<T>;
	return diff_phi;
}
template <typename T> constexpr auto delta_phi(const T dp)
	{return fastPrincipalRangeReductor(dp);}

inline auto deltaR(
	const double eta1,const double phi1,
	const double eta2,const double phi2)
	{return std::hypot(eta1-eta2,delta_phi(phi1-phi2));}//hypot from geometry

template<typename T,typename U>
[[gnu::const]] bool all_equal(const T& t, const U& u){return t == u;}
template<typename T,typename U,typename... Types>
[[gnu::const]] bool all_equal(const T& t, const U& u, Types const&... args)
	{return t == u && all_equal(u, args...);}
auto jet_lep_min_deltaR(
	 const doubles &jet_etas
	,const doubles &jet_phis
	,const double   lep_eta
	,const double   lep_phi
){
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
inline auto tight_jet_id(
	 const doubles& jet_lep_min_dRs
	,const  floats& pts
	,const  floats& etas
	,const    ints& ids
){
	//if(0<debug) std::cout<<"tight_jet_id"<<std::endl;
	return ( JET_PT__MIN <  pts             )
	&&     (  abs(etas)  <  JET_ETA_MAX     )
	&&     ( JET_ISO     <  jet_lep_min_dRs )
	&&     (      2      <= ids             );
}
inline auto jetCutter(const unsigned jmin,const unsigned jmax){
	if(3<debug) std::cout<<"jet cut "<< jmin <<" "<< jmax <<std::endl;
	return[=](const ints& jetmask){
		const auto nj = std:: count_if(
			jetmask.cbegin(),
			jetmask.  cend(),
			[](int i){return i;});// 0 is false
		return jmin <= nj && nj <= jmax;
	};
}
inline auto lep_gpt(const channel ch){
	return [=](
		 const floats &pt
		,const ints &mask
		,const ints &eidx
		,const ints &midx
	){
	int i;
	switch(ch){
		case elnu:{i = eidx[mask][0];break;}
		case munu:{i = midx[mask][0];break;}
	}
	if(-1!=i) return static_cast<double>(pt[i]); else return -1.;
	};
}
inline auto find_lead_mask(const doubles& vals,const ints& mask){
	//if(0<debug) std::cout<<"lead mask"<<std::endl;
	if(!all_equal(mask.size(),vals.size())) throw std::logic_error(
		"Collections must be the same size in lead_mask");
	if(mask.empty()) throw std::logic_error(
		"Collections must not be empty in lead_mask");
	const auto      masked_vals = mask * vals;
	ints  lead_mask(masked_vals.size());// vector of zeroes
	const auto max_idx = static_cast<size_t>(std::distance(
		masked_vals.cbegin(),
		max_element( masked_vals.cbegin(),
		             masked_vals.  cend())));// Leading bjet == highest pt.
	lead_mask[max_idx] = 1;
	return lead_mask;
}
inline double transverse_w_mass(
	 const double lep__pt
	,const double lep_phi
	,const double met__pt
	,const double met_phi
){
	return 2.
		* std::abs (std::sin( delta_phi(lep_phi
		          - static_cast<double>(met_phi))*.5))
		* std::sqrt(static_cast<double>(met__pt)
		                              * lep__pt);
}
auto runLBfilter(
	const std::map<size_t,std::vector<std::pair<size_t,size_t>>>
	&runLBdict
	,const bool MC
){
	return [&,MC](const unsigned int run,const unsigned int LB){
		if(MC) return true;
		auto search =  runLBdict.find(run);
		if(  search == runLBdict.cend()) return false;// Not Found TODO: true?
		if(std::any_of(search->second.cbegin(),
		               search->second.  cend(),
		               [LB](std::pair<size_t,size_t> p)
		               {return p.first <= LB && LB <= p.second;}))
		     return  true;
		else return false;
	};
}
inline auto blinding(
          const doubles bjets
         ,const doubles tjets){// This is the jet multiplicity technique,
        // suggested by Dr Duncan Leggat.
        // We take the number of bjets and remaining tight jets based on
        // the number of existing final jets. Therefore, we expect either
        // 3 bjets at final state or , 1 bjet and two other jets.
        return (bjets.size() != 3 && tjets.size() <= JETS_MAX) || // the fourth jet
               (bjets.size() != 1 && tjets.size() <= JETS_MAX)  ; // is the recoiled jet
        // from t-channel process
}
inline auto btagP(const doubles &eta){return abs(eta) < BJET_ETA_MAX;}
inline auto btagB(const ints  &btagP,const floats &btags){
	return   btagP && ( BTAG_DISC_MIN < btags );// all tJ length
}
} // namespace

void NPL(const channel ch,const dataSource ds){
	nlohmann::json JSONdict;
	std::ifstream(// open this JSON file once as a stream
	"aux/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt")
	>> JSONdict;// and read into this json object, then fix type of key
	std::map<size_t,std::vector<std::pair<size_t,size_t>>> runLBdict;
	for(const auto& [key,value] : JSONdict.items())
		runLBdict.emplace(std::stoi(key),value);// key:string->size_t
	// we now have one single copy of a wonderfully clean dictionary

	std::string temp_header="/data/disk3/nanoAOD_2017/",
	temp_opener,temp_footer="/*.root";/**/
	switch(ds){// CMS and MET MUST do some OPENABLE file ; reject later
	case tzq:{temp_opener="/data/disk3/nanoAOD_2017/tZqlvqq/*.root"  ;break;}/**/
	case  ww:{temp_opener=temp_header+ "WWToLNuQQ"       +temp_footer;break;}
	case  wz:{temp_opener=temp_header+ "WZTo1L1Nu2Q"     +temp_footer;break;}
	case  zz:{temp_opener=temp_header+ "ZZTo2L2Q"        +temp_footer;break;}
	case ttb:{temp_opener=temp_header+ "TTToSemileptonic"+temp_footer;break;}
	case ttz:{temp_opener=            "ttz_dir"          +temp_footer;break;}
	case cms:{temp_opener=temp_header+"ttZToQQ"          +temp_footer;break;}
	default :throw std::invalid_argument("Unimplemented ds (rdfopen)");
	}
	ROOT::RDataFrame mc__df("Events",temp_opener);// Monte Carlo
	// Open chains of exptData EVEN IF UNUSED
	TChain elnuCMS("Events");
	TChain munuCMS("Events");
	temp_footer = "/*.root" ;/* safety redefinition now saving us */
	temp_header =
		"/data/disk3/nanoAOD_2017/SingleElectron_NanoAOD25Oct2019_Run";
	for(std::string c:{"B","C","D","E","F"}){// guaranteed sequential
		temp_opener=temp_header+ c +temp_footer;
		elnuCMS.Add(temp_opener. c_str());
	}
	temp_header="/data/disk3/nanoAOD_2017/SingleMuon_NanoAOD25Oct2019_Run";
	for(std::string c:{"B","C","D","E","F"}){// guaranteed sequential
		temp_opener=temp_header+ c +temp_footer;
		munuCMS.Add(temp_opener. c_str());
	}

	ROOT::RDataFrame  elnudf(elnuCMS);
	ROOT::RDataFrame  munudf(munuCMS);
	const bool MC = !(met == ds ||cms == ds);
	auto df = [&,ch,ds](){// Get correct data frame
		switch(ds){
			case tzq:
			case  ww:// fall through!
			case  wz:
			case  zz:
			case ttb:
			case ttz:{           return mc__df;break;}
			case cms :{switch(ch){// MC is already false
			          case elnu:{return elnudf;break;}
			          case munu:{return munudf;break;}
			          default  :throw std::invalid_argument(
				"Unimplemented ch (rdf set)");
				}break;}
			default :throw std::invalid_argument(
				"Unimplemented ds (rdf set)");
		}
	}();
	switch(ch){
		case elnu:{temp_header = "Electron_";break;}
		case munu:{temp_header =     "Muon_";break;}
		default  :throw std::invalid_argument(
			"Unimplemented ch (init)");
	}
	// make test runs faster by restriction. Real run should not
//	auto dfr = df.Range(10000);// remember to enable MT when NOT range
	auto origi = df// toggle one letter to do all
	.Define("lep_pts","static_cast<ROOT::RVec<double>>("+temp_header+"pt)")
	;
	auto lumclean = origi
	.Filter(runLBfilter(runLBdict,MC),{"run","luminosityBlock"},
			   "LuminosityBlock filter")
	;
	auto offlep = lumclean
	// lepton selection first
	.Filter(Triggers(ch),
	{ "HLT_Ele32_WPTight_Gsf_L1DoubleEG"
	 ,"HLT_IsoMu27"
	 ,"HLT_PFHT250"
	},"Triggers Filter")
	.Define("loose_leps",lep_sel(ch),{
	        temp_header+"isPFcand"
	       ,temp_header+"pt"
	       ,temp_header+"eta"
	       ,"Electron_cutBased"
	       ,"Muon_tightId"
	       ,"Muon_pfRelIso04_all"
	       ,temp_header+"dxy"
	       ,temp_header+"dz"
	       })
	.Define("not_tight",not_tight(ch),{
		temp_header+"isPFcand"
               ,temp_header+"pt"
               ,temp_header+"eta"
               ,"Electron_cutBased"
               ,"Muon_tightId"
               ,"Muon_pfRelIso04_all"
               ,temp_header+"dxy"
               ,temp_header+"dz"
	       })
	.Define("rawJet_eta","static_cast<ROOT::RVec<double>>(Jet_eta )")
	.Define("rawJet_phi","static_cast<ROOT::RVec<double>>(Jet_phi )")
	.Define("tight_jets"       ,tight_jet_id,
	{"jet_lep_min_dR"   ,"Jet_pt","Jet_eta","Jet_jetId"})
	.Define("tJ_btagCSVv2"  ,"Jet_btagCSVV2[tight_jets]")
        .Define("lep__pt","static_cast<ROOT::RVec<double>>("+temp_header+     "pt[loose_leps])")
        .Define("lep_eta","static_cast<ROOT::RVec<double>>("+temp_header+    "eta[loose_leps])")
        .Define("lep_phi","static_cast<ROOT::RVec<double>>("+temp_header+    "phi[loose_leps])")
	.Define("lep_mas","static_cast<ROOT::RVec<double>>("+temp_header+   "mass[loose_leps])")
        .Define("jet_lep_min_dR"   ,jet_lep_min_deltaR,// later reused with doubles
               {"rawJet_eta","rawJet_phi","lep_eta","lep_phi"})// gcc fail template
	.Define("fin_jets_eta","static_cast<ROOT::RVec<double>>(Jet_eta  [tight_jets])")
	.Define("fin_jets__pt","static_cast<ROOT::RVec<double>>(Jet_pt   [tight_jets])")
	.Define("btagP"            ,btagP  ,{"fin_jets_eta"})// suPer vs suBset
	.Define("btagB"            ,btagB  ,{"btagP","tJ_btagCSVv2"})
	.Define( "lead_bjet"      , find_lead_mask,{"fin_jets__pt","btagB"})
        .Filter("All(fin_jets__pt > 35)", "after cjer jet pt cut")
        .Filter(blinding,{"bjet__pt","fin_jets__pt"},"blinding:jet multiplicity")
	// now for transverse W; lepton done
	. Alias("tw_lep__pt","lep__pt")
	. Alias("tw_lep_eta","lep_eta")
	. Alias("tw_lep_phi","lep_phi")
	.Define("tw_lep_mas",transverse_w_mass,
	       {"tw_lep__pt",
	        "tw_lep_phi","MET_pt","MET_phi"})
	.Filter(lep_tight_cut(ch),{"loose_leps",
	        "Electron_cutBased",// edit function for  tight -> loose
	        "Muon_pfRelIso04_all"},"lepton cut")// left with 1 tight lepton
        .Filter(not_tight_cut(ch),{"not_tight",
                "Electron_cutBased",// edit function for  tight -> loose
                "Muon_pfRelIso04_all"},"lepton cut")// left with 1 not_tight lepton
	;
	auto QCD_lep = offlep
	.Filter("All(tw_lep_mas < 30 && tw_lep_mass > 130)","Transverse W mass cut for QCD region,N_data")
	//.Filter("All(MET_sumEt  < 40)","Sum ET < 40 GeV   cut for QCD region")
	;
	auto sig_lep = offlep
	.Filter("All(tw_lep_mas > 30 && tw_lep_mass < 130)","Transverse W mass cut for Signal, real-misid")
	;
	if(MC){
	auto prompt_QCD_lep = QCD_lep //prompt for MC in QCD
	.Define("lep_gsf",lep_gpt(ch),{"GenPart_StatusFlags","not_tight"
                                      ,"Electron_genPartIdx"
                                      ,    "Muon_genPartIdx"})
	.Filter("All(lep_gsf == 0)","QCD PROMPT MC")
	;
        auto prompt_sig_lep = sig_lep //prompt for MC in sig
        .Define("lep_gsf",lep_gpt(ch),{"GenPart_StatusFlags","not_tight"
                                      ,"Electron_genPartIdx"
                                      ,    "Muon_genPartIdx"})
        .Filter("All(lep_gsf == 0)","Signal prompt MC")
        ;
	prompt_QCD_lep.Report() ->Print();
	prompt_sig_lep.Report() ->Print();
	sig_lep.Report() ->Print();// N_Real mis-ID for MC in signal region
	}
	else{
	// N_QCD_data
	QCD_lep.Report() ->Print();
	;
	}

} //void

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
	else if ( "tzq"  ==  dsN ){ d = tzq;}
        else if ( "ttz"  ==  dsN ){ d = ttz;}
        else if ( "ttb"  ==  dsN ){ d = ttb;}
        else if (  "ww"  ==  dsN ){ d =  ww;}
        else if (  "wz"  ==  dsN ){ d =  wz;}
        else if (  "zz"  ==  dsN ){ d =  zz;}
	else if ( "cms" ==  dsN  ){ d = cms;}
        else if ( "met"  ==  dsN ){ d = met;}
	else { std::cout << "Error: data source " << dsN
		<< " not recognised" << std::endl ;
		return 4 ;
	}
		NPL(c,d) ;
		return 0 ;
}
