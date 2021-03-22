// This script is designed to compute the NPL estimation of the data. The main NPL calculation is in calchisto.cpp

// Compiler:
//clang++ -Isrc -std=c++17 -march=native -pipe -O3 -Wall -Wextra -Wpedantic -o build/npl src/NPL.cxx `root-config --libs` -lm

/* Methodology:
// Non-prompt-lepton estimation
inline auto Npl
	const int N_L!T_data  // Number of loose not tight lepton in QCD region (control region) in data
	,const int N_L!T_prompt// Number of loose not tight prompt lepton in QCD region of enriched QCD sample (ttb,ttl,ttj,tz1,tz2,wz,zz,st,stb,stw,stbw)
	,const int eff_TL      // Number of tight to loose lepton ratio in QCD enriched sample in QCD region
){// produce three p_T vs eta histograms for the parameters above and save them as root files to be used in calchisto.cpp
  // The QCD regions are the regions outside 40 GeV of the nominal W mass
  // all Nqcd need to come from the loosend isolation
  // loose and tight criteria based on the isolation and identification cut,
  // 0.8 for relative isolation
  // impact parameters |dxy| and |dz| loosened to 0.1, 0.5
  // N SR_fake = (eff_TL / 1 -eff_TL) . (N L!T data - N L!T prompt)
  // this formula is implemented in the calchisto in the lepeffgiver
  // where the associated lepton pt and eta will be chosen from the three histograms
  // the by applying the formula the result for N_signalRegion_fake is found.
}
*/

#include <ROOT/RDataFrame.hxx>//#include <ROOT/RCsvDS.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <Math/Vector4D.h>
#include <TRandom3.h>// used Gaussian, uniform each once
//#include <execution>// need to link -ltbb in Makefile
#include <TChain.h>
#include <TF1.h>
//#include "calchisto.hpp"
#include "csv.h"
#include "json.hpp"
#include "eval_complex.hpp"
#include "roccor.Run2.v3/RoccoR.cc"

enum dataSource  {tzq,ww,wz,zz,ttb,ttj,ttl,tz1,tz2,st,stb,wjt,stw,stbw,wjqq,wzll,met,mc,cms};// mc : QCD enriched samples: tzq,tz1,tz2,ttb,ttl,ttj,st,stb,stw,stbw,zz,wz,ww
enum channel     {elnu,munu};

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

  constexpr double          TZQ_W =   .00128;
  constexpr double       WWLNQQ_W =   .21740;
  constexpr double       WZLNQQ_W =   .02335;
  constexpr double        TTBLV_W =   .13791;
  constexpr double        TZQQ1_W =   .02826;// ttz
  constexpr double        TZQQ2_W =   .00237;// ttz
  constexpr double       ZZLLQQ_W =   .00485;
  constexpr double	    WJT_W = 73.26469;
  constexpr double           ST_W =   .03837;
  constexpr double          STB_W =   .04433;
  constexpr double        TTBLL_W =   .05303;
  constexpr double        TTBJJ_W =   .12066;
  constexpr double          STW_W =   .18247;
  constexpr double         STBW_W =   .18750;
  constexpr double         WJQQ_W =   .17827;
  constexpr double         WZLL_W =   .00844;

  constexpr float    TRIG_SF_ELNU = 1.070511816990938379014537410816376133957408282644307576982;
  constexpr float    TRIG_SF_MUNU = 1.089437304676686344318303492342796025659913489719426417037;

template <typename T> constexpr T  PI = T(3.14159265358979323846264338327950288419716939937510582097494459230781640628620899);
template <typename T> constexpr T TPI = PI<T> * 2;

enum      elSf      {Eff,Smr};
enum      muSf      {Id_N,Id_Y,Id_A,Id_T,//  Id, Idsys, IdsysStat, IdsysSyst,
                     IsoN,IsoY,IsoA,IsoT};//Iso,Isosys,IsosysStat,IsosysSyst};

inline auto lep_sel(const channel ch){
	if(0<debug) std::cout<<"lep sel start"<<std::endl;
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
	if(0<debug) std::cout<<"not tight"<<std::endl;
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
				&&	 pts  >  EL__PT_MIN
				&& ((abs_etas <  EL_ETA_MAX
				&&        dz  <  ENDCAP_DZ__TIGHT
				&&        dxy <  ENDCAP_DXY_TIGHT
				&&   abs_etas >  ENDCAP_ETA_MIN)
				||  (abs_etas <  BARREL_ETA_MAX
				&&        dxy <  BARREL_DXY_TIGHT
				&&        dz  <  BARREL_DZ__TIGHT))
				&&	elids <  EL_TIGHT_ID // Make it not tightable
			);
			case munu:return ( true
				&&   muids
				&&   isPFs
				&&     pts  >  MU__PT_MIN
				&& abs_etas <  MU_ETA_MAX
				&&     isos >  MU_TIGHT_ISO // Make it not tightable
			);
			default  :throw std::invalid_argument(
				"Unimplemented ch (lep_sel)");
		}
	};
}
inline auto lep_tight_cut(const channel ch){
	if(0<debug) std::cout<<"lep tight cut"<<std::endl;
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
	if(0<debug) std::cout<<"not tight cut"<<std::endl;
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
	if(0<debug) std::cout<<"Triggers"<<std::endl;
	return [=](
	 const bool el// HLT_Ele32_WPTight_Gsf_L1DoubleEG
	,const bool mu// HLT_IsoMu27
){
	switch(ch){
	case elnu:return el;
	case munu:return mu;
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
	if(0<debug) std::cout<<"jet_lep_min_dR"<<std::endl;
	doubles min_dRs; min_dRs.reserve(jet_etas.size());
	std::transform(
		jet_etas.cbegin(),
		jet_etas.  cend(),
		jet_phis.cbegin(),
		std::back_inserter(min_dRs),
		[=](double jet_eta, double jet_phi){
			return deltaR(jet_eta,jet_phi,lep_eta,lep_phi);});
	if(0<debug) std::cout<<"jet lep min finish"<<std::endl;
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
	if(3<debug) std::cout<<"jet cut "<<std::endl;
	return[=](const ints& jetmask){
		const auto nj = std:: count_if(
			jetmask.cbegin(),
			jetmask.  cend(),
			[](int i){return i;});// 0 is false
		if(0<debug) std::cout<<"jet cut finish"<<std::endl;
		return jmin <= nj && nj <= jmax;
	};
}
inline auto lep_gpsf(const channel ch){
	if(0<debug) std::cout<<"lep gpsf "<<std::endl;
	return [=](
		 const ints &gpsf
		,const ints &mask
		,const ints &eidx
		,const ints &midx
	){
	ints i;
	switch(ch){
		if(0<debug) std::cout<<"lep gpsf in switch"<<std::endl;
		case elnu:{i = eidx[mask];break;}
		case munu:{i = midx[mask];break;}
	}
	// GenPart_statusFlafs is stored as bits. it first should be sorted
	// by The idx and then masked correctly (to the tight lepton or loose ones)
	// The zero bit shows whether the particle is prompt or not. Thus we check
	// whether the remainer of the gpsf % 2 is 1 or not.
	// * the -1!=i is responsible for the indexcies which are not related to a
	// lepton or jet level therefore could be taken out.
	if(0<debug) std::cout<<"lep gpsf before Take"<<std::endl;
	ints g = Take(gpsf,i[-1!=i]) & (1 << 0);
	if(0<debug) std::cout<<"lep gpsf finish"<<std::endl;
	return Any(g == 1);
	};
}
inline auto find_lead_mask(const doubles& vals,const ints& mask){
	//if(0<debug) std::cout<<"lead mask"<<std::endl;
	if(!all_equal(mask.size(),vals.size())) throw std::logic_error(
		"Collections must be the same size in lead_mask");
	if(mask.empty()) throw std::logic_error(
		"Collections must not be empty in lead_mask");
	if(0<debug) std::cout<<"lead mask"<<std::endl;
	const auto      masked_vals = mask * vals;
	ints  lead_mask(masked_vals.size());// vector of zeroes
	const auto max_idx = static_cast<size_t>(std::distance(
		masked_vals.cbegin(),
		max_element( masked_vals.cbegin(),
		             masked_vals.  cend())));// Leading bjet == highest pt.
	lead_mask[max_idx] = 1;
	if(0<debug) std::cout<<"lead mask"<<std::endl;
	return lead_mask;
}
inline double transverse_w_mass(
	 const double lep__pt
	,const double lep_phi
	,const  float met__pt
	,const  float met_phi
){
	if(0<debug) std::cout<<"transverse w mass "<<std::endl;
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
	if(0<debug) std::cout<<"runLBfilter"<<std::endl;
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
inline auto btagP(const doubles &eta){return abs(eta) < BJET_ETA_MAX;}
inline auto btagB(const ints  &btagP,const floats &btags){
	return   btagP && ( BTAG_DISC_MIN < btags );// all tJ length
}
inline auto sf(const  dataSource ds,
               const  channel    ch)
{
	if(0<debug) std::cout<<"scale factor "<<std::endl;
	return [=](
		 const double mostSF
	){
		// TODO: trigger efficiency
		double result;
		bool MC = true;
		switch(ds){
			case  tzq:{result =    TZQ_W;break;}
			case   ww:{result = WWLNQQ_W;break;}
			case   wz:{result = WZLNQQ_W;break;}
			case   zz:{result = ZZLLQQ_W;break;}
			case  wjt:{result =    WJT_W;break;}
			case   st:{result =     ST_W;break;}
			case  stb:{result =    STB_W;break;}
			case  stw:{result =    STW_W;break;}
			case stbw:{result =   STBW_W;break;}
			case wjqq:{result =   WJQQ_W;break;}
			case wzll:{result =   WZLL_W;break;}
			case  ttb:{result =  TTBLV_W;break;}
			case  ttl:{result =  TTBLL_W;break;}
			case  ttj:{result =  TTBJJ_W;break;}
			case  tz1:{result =  TZQQ1_W;break;}
			case  tz2:{result =  TZQQ2_W;break;}
			case  mc:{result = WZLNQQ_W+
			                   ZZLLQQ_W+
			                    TTBLV_W+
					    TTBLL_W+
					    TTBJJ_W+
			                    TZQQ1_W+
					    TZQQ2_W+
					      WJT_W+
					       ST_W+
					      STB_W+
					      STW_W+
					     STBW_W+
					     WJQQ_W ;break;}
			case met:// fall through to cms and mc
			case cms:{result = 1.;MC=false;break;}// ignore btag wt
			default :throw std::invalid_argument(
				"Unimplemented ds (sf)");
		}
		switch(ch){
			case elnu:{result *= TRIG_SF_ELNU;break;}
			case munu:{result *= TRIG_SF_MUNU;break;}
		}
		result *=  mostSF;
		if(0<debug) std::cout<<"sf finish"<<std::endl;
		return result;
	};
}
auto roccorSF(
	 RoccoR          &rc
	,const channel    ch
	,const bool       MC
){
	return [=](
		 const double     pt,const double eta
		,const double    phi,const   int  Q
		,const double gen_pt,const   int  nl
	){
	if(0 < debug)std::cout<< "roccor SF"<<std::endl;
		double roc = 1.;
		switch(ch){
		case elnu:{break;}
		case munu:{if(MC){
			if(0. < gen_pt){
				roc = rc.kSpreadMC(Q,pt,eta,phi,gen_pt,0,0);
			}else{
				auto u = gRandom->Rndm();
				roc = rc. kSmearMC(Q,pt,eta,phi,nl,u,0,0);
			}
		}else{
				roc = rc. kScaleDT(Q,pt,eta,phi,0,0);
		}break;}}// not MC, case munu, switch
		if(0 < debug) std::cout << "roccor finish" << roc << std::endl;
		return roc;
	};
}
auto elEffGiver(
	 const double             pt
	,const double             eta
	,const TH2F* const &recoLowEt
	,const TH2F* const &reco_pass
	,const TH2F* const &tight_94x
){
	if(0<debug) std::cout<<"elEffGiver"<<std::endl;
	std::map<elSf,double> dict = {{Eff,1.},{Smr,1.}};
	// eff == electron    regression    corrections
	// smr == energy scale and smearing corrections
	if(2.5 < std::abs(eta)) return dict;
	int PtBin,EtaBin;
	if(2<debug)std::cout<<"el eff giver"<<std::endl;
	if( pt      < 20.){
		EtaBin   = recoLowEt->GetXaxis()-> FindBin(eta);
		 PtBin   = recoLowEt->GetYaxis()-> FindBin(pt );
		dict[Eff]= recoLowEt->GetBinContent(EtaBin,PtBin);
	}else{
		EtaBin   = reco_pass->GetXaxis()-> FindBin(eta);
		 PtBin   = reco_pass->GetYaxis()-> FindBin(pt );
		dict[Eff]= reco_pass->GetBinContent(EtaBin,PtBin);
	}//else{// TODO: this need clarification
		EtaBin   = tight_94x->GetXaxis()-> FindBin(eta);
		 PtBin   = tight_94x->GetYaxis()-> FindBin(pt );
		dict[Smr]= tight_94x->GetBinContent(EtaBin,PtBin);
//	}
	// TODO: GetBinError
	if(5<debug)std::cout<<dict[Eff]<<" eff , smr "<<dict[Smr]<<std::endl;
	// TODO: Eff == 0 cross check
	if(FP_NORMAL != std::fpclassify(dict[Eff])) dict[Eff] = 1.;
	if(0<debug) std::cout<<"elleffgiver finish"<<std::endl;
	return dict;
}
auto muEffGiver(
	 const double        pt
	,const double        eta
	,const TH2D* const &id_N
	,const TH2D* const &id_Y
	,const TH2D* const &id_A
	,const TH2D* const &id_T
	,const TH2D* const &isoN
	,const TH2D* const &isoY
	,const TH2D* const &isoA
	,const TH2D* const &isoT
){
	if(0<debug) std::cout<<"muEffGiver"<<std::endl;
	const double ata = std::abs(eta);
	std::map<muSf,double> dict
	={{Id_N,1.},{Id_Y,1.},{Id_A,1.},{Id_T,1. },
	  {IsoN,1.},{IsoY,1.},{IsoA,1.},{IsoT,1.}};
	if(2<debug)std::cout<<"mu eff giver"<<std::endl;
	if(pt < 20 || pt > 120 || ata > 2.4) return dict;
	int PtBin,EtaBin;
	 PtBin     = id_N->GetXaxis()->FindBin(pt );
	EtaBin     = id_N->GetYaxis()->FindBin(ata);
	dict[Id_N] = id_N->GetBinContent(PtBin,EtaBin);
	 PtBin     = id_Y->GetXaxis()->FindBin(pt );
	EtaBin     = id_Y->GetYaxis()->FindBin(ata);
	dict[Id_Y] = id_Y->GetBinContent(PtBin,EtaBin);
	 PtBin     = id_A->GetXaxis()->FindBin(pt );
	EtaBin     = id_A->GetYaxis()->FindBin(ata);
	dict[Id_A] = id_A->GetBinError  (PtBin,EtaBin);
	 PtBin     = id_T->GetXaxis()->FindBin(pt );
	EtaBin     = id_T->GetYaxis()->FindBin(ata);
	dict[Id_T] = id_T->GetBinError  (PtBin,EtaBin);

	 PtBin     = isoN->GetXaxis()->FindBin(pt );
	EtaBin     = isoN->GetYaxis()->FindBin(ata);
	dict[IsoN] = isoN->GetBinContent(PtBin,EtaBin);
	 PtBin     = isoY->GetXaxis()->FindBin(pt );
	EtaBin     = isoY->GetYaxis()->FindBin(ata);
	dict[IsoY] = isoY->GetBinContent(PtBin,EtaBin);
	 PtBin     = isoA->GetXaxis()->FindBin(pt );
	EtaBin     = isoA->GetYaxis()->FindBin(ata);
	dict[IsoA] = isoA->GetBinError  (PtBin,EtaBin);
	 PtBin     = isoT->GetXaxis()->FindBin(pt );
	EtaBin     = isoT->GetYaxis()->FindBin(ata);
	dict[IsoT] = isoT->GetBinError  (PtBin,EtaBin);

	if(5<debug)std::cout<<dict[IsoA]
	                    <<"\t:\t"<<dict[IsoT]<<std::endl;
	if(0<debug)std::cout<<"mueffgiver finsih"<<std::endl;
	return dict;
}
auto lepEffGiver(
	 const channel      ch
	,const bool         MC
	,const TH2F* const &recoLowEt
	,const TH2F* const &reco_pass
	,const TH2F* const &tight_94x
	,const TH2D* const &id_N
	,const TH2D* const &id_Y
	,const TH2D* const &id_A
	,const TH2D* const &id_T
	,const TH2D* const &isoN
	,const TH2D* const &isoY
	,const TH2D* const &isoA
	,const TH2D* const &isoT
){
	if(0<debug) std::cout<<"lepEffGiver"<<std::endl;
	return [=](const double pt,const double eta){
	if(!MC) return 1.;
//	if(0 < debug)std::cout<< "lep eff giver"<<std::endl;
	double id = 1., iso = 1., eff = 1., smr = 1.;
	switch(ch){
	case elnu:{
		auto  dict = elEffGiver(pt,eta,recoLowEt,reco_pass,tight_94x);
		eff = dict[Eff];
		smr = dict[Smr];
		if(0<debug)std::cout<<"eff "<<eff<<" smr "<<smr<<std::endl;
/*		if(eff < 0)std::cout<<"eff lt 0 ";
		if(smr < 0)std::cout<<"smr lt 0 ";
		if(eff < 0 || smr < 0)std::cout<<std::endl;
*/		break;}
	case munu:{
		auto  dict = muEffGiver(pt,eta
		                       ,id_N,id_Y,id_A,id_T
		                       ,isoN,isoY,isoA,isoT
		);
		id  = dict[Id_N];
		iso = dict[IsoN];
		if(0<debug)std::cout<<"id  "<<id <<" iso "<<iso<<std::endl;
/*		if(id  < 0)std::cout<<"id  lt 0 ";
		if(iso < 0)std::cout<<"iso lt 0 ";
		if(id  < 0 || iso < 0)std::cout<<std::endl;
*/		break;}
	}
	if(0<debug) std::cout<<"lepeffgiver finish"<<std::endl;
	return id * iso * eff * smr;};
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
	// now we repeat for some other files
	TFile *tF ;
	TH2F  *tHf; TH2D *tHd; TH1D* t1d;
	// electron efficiencies
	tF = TFile::Open(
		"aux/elEff/egammaEffi.txt_EGM2D_runBCDEF_passingRECO_lowEt.root");
	tF ->GetObject("EGamma_SF2D",tHf);tHf->SetDirectory(nullptr);
	const TH2F* const recoLowEt = static_cast<TH2F*>(tHf);
	tF ->Close();
	tF = TFile::Open(
		"aux/elEff/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root");
	tF ->GetObject("EGamma_SF2D",tHf);tHf->SetDirectory(nullptr);
	const TH2F* const reco_pass = static_cast<TH2F*>(tHf);
	tF ->Close();
	tF = TFile::Open(
		"aux/elEff/egammaEffi.txt_EGM2D_runBCDEF_passingTight94X.root");
	tF ->GetObject("EGamma_SF2D",tHf);tHf->SetDirectory(nullptr);
	const TH2F* const tight_94x = static_cast<TH2F*>(tHf);
	tF ->Close();
	// muon efficiencies
	tF = TFile::Open("aux/muEff/Muon_RunBCDEF_SF_ID.root");
	tF ->GetObject("NUM_TightID_DEN_genTracks_pt_abseta",tHd);
	tHd->SetDirectory(nullptr);// make it stay even if file closed
	const TH2D* const id_N = static_cast<TH2D*>(tHd);
	tF ->Close();
	tF = TFile::Open("aux/muEff/Muon_RunBCDEF_SF_ID_syst.root");
	tF ->GetObject("NUM_TightID_DEN_genTracks_pt_abseta",tHd);
	tHd->SetDirectory(nullptr);
	const TH2D* const id_Y = static_cast<TH2D*>(tHd);
	tF ->GetObject("NUM_TightID_DEN_genTracks_pt_abseta_stat",tHd);
	tHd->SetDirectory(nullptr);
	const TH2D* const id_A = static_cast<TH2D*>(tHd);
	tF ->GetObject("NUM_TightID_DEN_genTracks_pt_abseta_syst",tHd);
	tHd->SetDirectory(nullptr);
	const TH2D* const id_T = static_cast<TH2D*>(tHd);
	tF ->Close();
	tF = TFile::Open("aux/muEff/Muon_RunBCDEF_SF_ISO.root");
	tF ->GetObject("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta",tHd);
	tHd->SetDirectory(nullptr);
	const TH2D* const isoN = static_cast<TH2D*>(tHd);
	tF ->Close();
	tF = TFile::Open("aux/muEff/Muon_RunBCDEF_SF_ISO_syst.root");
	tF ->GetObject("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta",tHd);
	tHd->SetDirectory(nullptr);
	const TH2D* const isoY = static_cast<TH2D*>(tHd);
	tF ->GetObject("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta_stat",tHd);
	tHd->SetDirectory(nullptr);
	const TH2D* const isoA = static_cast<TH2D*>(tHd);
	tF ->GetObject("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta_syst",tHd);
	tHd->SetDirectory(nullptr);
	const TH2D* const isoT = static_cast<TH2D*>(tHd);
	tF ->Close();
	tF	= nullptr; tHf = nullptr; tHd = nullptr; t1d = nullptr;
	RoccoR rc("src/roccor.Run2.v3/RoccoR2017.txt");
	// All require files for lepton SF are read

	std::string temp_header="/data/disk0/nanoAOD_2017/",
	temp_opener,temp_footer="/*.root";/**/
	switch(ds){// CMS and MET MUST do some OPENABLE file ; reject later
	case cms:{temp_opener=temp_header+ "ttZToQQ"         +temp_footer;break;}
	case  mc:{temp_opener=temp_header+ "ttZToQQ"         +temp_footer;break;}
	default :throw std::invalid_argument("Unimplemented ds (rdfopen)");
	}
	// Open chains of exptData EVEN IF UNUSED
	std::string temp_header0 = "/data/disk0/nanoAOD_2017/",
		    temp_header1 = "/data/disk1/nanoAOD_2017_new/",
		    temp_header3 = "/data/disk3/nanoAOD_2017/";
	// QCD Enriched MCs file openning in chain
	TChain QCDofMC("Events"); // including the enriched QCD datasets of MC
	QCDofMC.Add((temp_header3+ "tZqlvqq"                      +temp_footer).c_str());
	QCDofMC.Add((temp_header0+ "WWToLNuQQ"                    +temp_footer).c_str());
	QCDofMC.Add((temp_header0+ "WZTo1L1Nu2Q"                  +temp_footer).c_str());
	QCDofMC.Add((temp_header0+ "ZZTo2L2Q"                     +temp_footer).c_str());
	QCDofMC.Add((temp_header3+ "WPlusJets_NanoAODv5"          +temp_footer).c_str());
	QCDofMC.Add((temp_header0+ "TTToSemileptonic"             +temp_footer).c_str());
	QCDofMC.Add((temp_header1+ "TT_2l2nu_nanoAODv5"           +temp_footer).c_str());
	QCDofMC.Add((temp_header0+ "TTToHadronic"                 +temp_footer).c_str());
	QCDofMC.Add((temp_header1+ "ST_tchannel_top_nanoAODv5"    +temp_footer).c_str());
	QCDofMC.Add((temp_header1+ "ST_tchannel_antitop_nanoAODv5"+temp_footer).c_str());
	QCDofMC.Add((temp_header0+ "ST_tW"                        +temp_footer).c_str());
	QCDofMC.Add((temp_header0+ "ST_tbarW"                     +temp_footer).c_str());
	QCDofMC.Add((temp_header0+ "WZTo2L2Q"                     +temp_footer).c_str());
	QCDofMC.Add((temp_header0+ "WPlusJetsToQQ"                +temp_footer).c_str());
	QCDofMC.Add((temp_header0+ "ttZToQQ"                      +temp_footer).c_str());
	QCDofMC.Add((temp_header0+ "ttZToQQ_ext"                  +temp_footer).c_str());
	ROOT::RDataFrame mc__df(QCDofMC); // QCD enriched MCs

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
			case tz1:
			case tz2:
			case  mc:{           return mc__df;break;}
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
	auto dfr = df.Range(10000);// remember to enable MT when NOT range
	auto origi = df// toggle one letter to do all
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
	;
	// QCD Region loose
	auto loose_lep = offlep // QCD region selection and filter
	.Filter(not_tight_cut(ch),{"not_tight",
	        "Electron_cutBased",// edit function for  tight -> loose
	        "Muon_pfRelIso04_all"},"lepton cut")// left with 1 not_tight lepton
	.Define("lepB_pt","static_cast<double>("+temp_header+     "pt[not_tight][0])")
	.Define("lep_eta","static_cast<double>("+temp_header+    "eta[not_tight][0])")
	.Define("lep_phi","static_cast<double>("+temp_header+    "phi[not_tight][0])")
	.Define("lep_mas","static_cast<double>("+temp_header+   "mass[not_tight][0])")
	.Define("lep___q",temp_header+"charge[not_tight][0]")// int, only for RoccoR
	.Define("roccorSF", roccorSF(rc,ch,MC)// used in allReconstruction
	                  ,{"lepB_pt","lep_eta",
	                    "lep_phi","lep___q", // last 4 unused
	                    "lep_phi","lep___q"})// last 2 repeat is fine
	.Define("lep__pt","lepB_pt * roccorSF")// only pt and mass scales
	.Define("lepSF"  , lepEffGiver(ch,MC// not in reco because files
	                 , recoLowEt,reco_pass,tight_94x
	                 , id_N,id_Y,id_A,id_T
	                 , isoN,isoY,isoA,isoT
	                ),{"lep__pt","lep_eta"})
	. Alias("mostSF" , "lepSF")
	.Define("sf",sf(ds,ch),{"mostSF"})
	.Define("jet_lep_min_dR"   ,jet_lep_min_deltaR,// later reused with doubles
	       {"rawJet_eta","rawJet_phi","lep_eta","lep_phi"})// gcc fail template
	.Define("tight_jets" 	   ,tight_jet_id,
	{"jet_lep_min_dR"   ,"Jet_pt","Jet_eta","Jet_jetId"})
	.Filter( jetCutter(JETS_MIN,JETS_MAX),{"tight_jets" },"Jet cut")
	.Define("tJ_btagCSVv2"  ,"Jet_btagCSVV2[tight_jets]")
	.Define("fin_jets_eta","static_cast<ROOT::RVec<double>>(Jet_eta  [tight_jets])")
	.Define("fin_jets__pt","static_cast<ROOT::RVec<double>>(Jet_pt   [tight_jets])")
	.Define("btagP"            ,btagP  ,{"fin_jets_eta"})// suPer vs suBset
	.Define("btagB"            ,btagB  ,{"btagP","tJ_btagCSVv2"})
	.Define( "lead_bjet"      , find_lead_mask,{"fin_jets__pt","btagB"})
	.Define(      "bjet__pt"  ,"fin_jets__pt[lead_bjet]")// Leading bj
	// now for transverse W; lepton done
	. Alias("tw_lep__pt","lep__pt")
	. Alias("tw_lep_eta","lep_eta")
	. Alias("tw_lep_phi","lep_phi")
	.Define("tw_lep_mas",transverse_w_mass,
	       {"tw_lep__pt",
	        "tw_lep_phi","MET_pt","MET_phi"})
	.Filter("tw_lep_mas < 30 || tw_lep_mas > 130","Transverse W mass cut for loose QCD region,loose QCD")
	;
	// QCD region tight
	auto tight_lep = offlep // Signal region selection and filter
	.Filter(lep_tight_cut(ch),{"loose_leps",
	        "Electron_cutBased",// edit function for  tight -> loose
	        "Muon_pfRelIso04_all"},"lepton cut")// left with 1 not_tight lepton
	.Define("lepB_pt","static_cast<double>("+temp_header+     "pt[loose_leps][0])")
	.Define("lep_eta","static_cast<double>("+temp_header+    "eta[loose_leps][0])")
	.Define("lep_phi","static_cast<double>("+temp_header+    "phi[loose_leps][0])")
	.Define("lep_mas","static_cast<double>("+temp_header+   "mass[loose_leps][0])")
	.Define("lep___q",temp_header+"charge[loose_leps][0]")// int, only for RoccoR
	.Define("roccorSF", roccorSF(rc,ch,MC)// used in allReconstruction
	                  ,{"lepB_pt","lep_eta",
	                    "lep_phi","lep___q", // last 4 unused
	                    "lep_phi","lep___q"})// last 2 repeat is fine
	.Define("lep__pt","lepB_pt * roccorSF")// only pt and mass scales
	.Define("lepSF"  , lepEffGiver(ch,MC// not in reco because files
	                 , recoLowEt,reco_pass,tight_94x
	                 , id_N,id_Y,id_A,id_T
	                 , isoN,isoY,isoA,isoT
	                ),{"lep__pt","lep_eta"})
	. Alias("mostSF" , "lepSF")
	.Define("sf",sf(ds,ch),{"mostSF"})
	.Define("jet_lep_min_dR"   ,jet_lep_min_deltaR,// later reused with doubles
	       {"rawJet_eta","rawJet_phi","lep_eta","lep_phi"})// gcc fail template
	.Define("tight_jets" 	   ,tight_jet_id,
	{"jet_lep_min_dR"   ,"Jet_pt","Jet_eta","Jet_jetId"})
	.Filter( jetCutter(JETS_MIN,JETS_MAX),{"tight_jets" },"Jet cut")
	.Define("tJ_btagCSVv2"  ,"Jet_btagCSVV2[tight_jets]")
	.Define("fin_jets_eta","static_cast<ROOT::RVec<double>>(Jet_eta  [tight_jets])")
	.Define("fin_jets__pt","static_cast<ROOT::RVec<double>>(Jet_pt   [tight_jets])")
	.Define("btagP"            ,btagP  ,{"fin_jets_eta"})// suPer vs suBset
	.Define("btagB"            ,btagB  ,{"btagP","tJ_btagCSVv2"})
	.Define( "lead_bjet"      , find_lead_mask,{"fin_jets__pt","btagB"})
	.Define(      "bjet__pt"  ,"fin_jets__pt[lead_bjet]")// Leading bj
	// now for transverse W; lepton done
	. Alias("tw_lep__pt","lep__pt")
	. Alias("tw_lep_eta","lep_eta")
	. Alias("tw_lep_phi","lep_phi")
	.Define("tw_lep_mas",transverse_w_mass,
	       {"tw_lep__pt",
	        "tw_lep_phi","MET_pt","MET_phi"})
	.Filter("All(fin_jets__pt > 35)", "after cjer jet pt cut")
	.Filter("tw_lep_mas > 130 || tw_lep_mas < 30","Transverse W mass cut for tight QCD region, tight in QCD")
	;// making sure signal region
	switch(ch){
		case elnu:{temp_header="elnu_";temp_footer="elnu_";break;}
		case munu:{temp_header="munu_";temp_footer="munu_";break;}
	}
	switch(ds){
		case  mc:{temp_header+="QCD";temp_footer+="QCD";break;}// QCD enriched sample
		case cms:{temp_header+="cms";temp_footer+="CMS";break;}
//		default :throw std::invalid_argument(
//			"Unimplemented ds (hist titles)");
	}
	// Bins definition -> this method lead to splitting bins which are more dense
	// also groups the less dense bin together
	const std::vector<double> ptBins
	= {0.,20.,30.,40.,45.,50.,55.,60.,65.,75.,85.,100.,150.,500.};
	const int ptBinsSize = ptBins.size() - 1;

	if(MC){
	auto h_tight = tight_lep.Histo2D({ // histo for eff tight to loose , tight
	("tight_Lepton_" + temp_header).c_str(),
	("tight Lepton " + temp_footer).c_str(),
	ptBinsSize,ptBins.data(),20,-2.5,2.5},
	"lep__pt","lep_eta")
	;
	auto h_loose = loose_lep.Histo2D({ // histo for eff tight to loose , loose
	("loose_Lepton_" + temp_header).c_str(),
	("loose Lepton " + temp_footer).c_str(),
	ptBinsSize,ptBins.data(),20,-2.5,2.5},
	"lep__pt","lep_eta")
	;
	auto
	h_eff_TL_ratio = static_cast<TH2D*>(h_tight->Clone());
	h_eff_TL_ratio->Divide(             h_loose.GetPtr());
	h_eff_TL_ratio->Draw("COLZ");
	h_eff_TL_ratio->SetNameTitle(
	("TL_eff_" + temp_header).c_str(),
	("TL_eff " + temp_header).c_str())
	;
	auto prompt_loose_lep = loose_lep //prompt for MC in QCD
	.Filter(lep_gpsf(ch),{"GenPart_statusFlags","not_tight"
	                              ,"Electron_genPartIdx"
	                              ,    "Muon_genPartIdx"},"loose not tight prompt MC")
	;
	auto h_prompt_loose_lep = prompt_loose_lep.Histo2D({// N_L!T prompt MC
	("prompt_LnT_"    + temp_header).c_str(),
	("ptompt LnT MC " + temp_footer).c_str(),
	ptBinsSize,ptBins.data(),20,-2.5,2.5},
	"lep__pt","lep_eta")
	;

	TFile hf(("histo/NPL_"+temp_header+".root").c_str(),"RECREATE");
	// MC only
	hf.WriteTObject(h_prompt_loose_lep  .GetPtr());hf.Flush();sync();
	hf.WriteTObject(h_eff_TL_ratio               );hf.Flush();sync();
	}
	else{
	// N_L!T_data
	auto h_LT_data = loose_lep.Histo2D({
	("N_data_LnT_" + temp_header).c_str(),
	("N data LnT " + temp_footer).c_str(),
	ptBinsSize,ptBins.data(),20,-2.5,2.5},
	"lep__pt","lep_eta")
	; // L!T data
	TFile hf(("histo/NPL_"+temp_header+".root").c_str(),"RECREATE");
	// CMS only
	hf.WriteTObject(h_LT_data           .GetPtr());hf.Flush();sync();

	}
	std::cout<<"NPL successfully finished"<<std::endl;
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
	else if ( "mc"  ==   dsN  ){d  = mc;}
	else if ( "cms" ==   dsN  ){d = cms;}
	else { std::cout << "Error: data source " << dsN
		<< " not recognised" << std::endl ;
		return 4 ;
	}
		NPL(c,d) ;
		return 0 ;
}
