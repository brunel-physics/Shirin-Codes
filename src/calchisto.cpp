// TODO: Adding Triggers
// TODO: lepton trigger efficiency
// TODO: Shape uncertainties
// TODO: top pt reweighing -> fetch pt min and max from ttb
// TODO: non prompt lepton corrections
// TODO: Pile Up uncertainties (done: weight)
// TODO: MET unclustering correction

#include <ROOT/RDataFrame.hxx>//#include <ROOT/RCsvDS.hxx>
#include <Math/Vector4D.h>
#include <TRandom3.h>// used Gaussian, uniform each once
//#include <execution>// need to link -ltbb in Makefile
#include <TChain.h>

#include "csv.h"
#include "json.hpp"
#include "calchisto.hpp"
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
//constexpr    int EL_MAX_NUM     = 1      ;
  constexpr  float EL__PT_MIN     = 35.f   ;
  constexpr  float EL_ETA_MAX     = 2.5f   ;
  constexpr    int EL_LOOSE_ID    = 1      ;
  constexpr    int EL_TIGHT_ID    = 4      ;
  constexpr  float ENDCAP_DZ      =  .2f   ;
  constexpr  float ENDCAP_DXY     =  .1f   ;
  constexpr  float ENDCAP_ETA_MIN = 1.5660f;
  constexpr  float BARREL_ETA_MAX = 1.4442f;
  constexpr  float BARREL_DXY     =  .05f  ;
  constexpr  float BARREL_DZ      =  .10f  ;

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
  constexpr double     W_MASS_CUT =  20.    ;
  constexpr double   TOP_MASS     = 172.5   ;
//constexpr double   TOP_MASS_CUT =  20.    ;

  constexpr float     JET_ETA_MAX =  4.7f;
  constexpr float     JET_PT__MIN = 30.0f;
  constexpr double        JET_ISO =   .4 ;
  constexpr unsigned     JETS_MIN =  4   ;
  constexpr unsigned     JETS_MAX =  6   ;

  constexpr double   BJET_ETA_MAX = 2.4   ;
  constexpr double  BTAG_DISC_MIN =  .8838;
  constexpr unsigned    BJETS_MIN = 1     ;
  constexpr unsigned    BJETS_MAX = 3     ;

//constexpr double DELTA___R_ZL   = 1.6;
//constexpr double DELTA_PHI_ZW   = 2. ;
//constexpr double DELTA_PHI_ZMET = 2. ;

  constexpr double    ak4RconeBy2 =  .2;
//constexpr double    ak8RconeBy2 =  .4;

  constexpr double          TZQ_W =  .0128;
  constexpr double       WWLNQQ_W = 2.1740;
  constexpr double       WZLNQQ_W =  .2335;
  constexpr double        TTBLV_W = 1.3791;
  constexpr double        TTZQQ_W =  .0237;
  constexpr double       ZZLLQQ_W =  .0485;

// This Pi is more accurate than binary256; good for eternity
template <typename T> constexpr T  PI = T(3.14159265358979323846264338327950288419716939937510582097494459230781640628620899);
template <typename T> constexpr T TPI = PI<T> * 2;

enum      puSf      {puW,upW,dnW};
enum      elSf      {Eff,Smr};
enum      muSf      {Id_N,Id_Y,Id_A,Id_T,//  Id, Idsys, IdsysStat, IdsysSyst,
                     IsoN,IsoY,IsoA,IsoT};//Iso,Isosys,IsosysStat,IsosysSyst};
/*
constexpr muSf
          muSfAll[]={Id_N,Id_Y,Id_A,Id_T,//  Id, Idsys, IdsysStat, IdsysSyst,
                     IsoN,IsoY,IsoA,IsoT};//Iso,Isosys,IsosysStat,IsosysSyst};
constexpr elSf
          elSfAll[]={Eff,Smr};
*/

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
				&&       pts  >  EL__PT_MIN
				&& ((abs_etas <  EL_ETA_MAX
				&&        dz  <  ENDCAP_DZ
				&&        dxy <  ENDCAP_DXY
				&&   abs_etas >  ENDCAP_ETA_MIN)
				||  (abs_etas <  BARREL_ETA_MAX
				&&        dxy <  BARREL_DXY
				&&        dz  <  BARREL_DZ))
				&&      elids >= EL_LOOSE_ID
			);
			case munu:return ( true
				&&   muids
				&&   isPFs
				&&     pts  >  MU__PT_MIN
				&& abs_etas <  MU_ETA_MAX
				&&     isos <= MU_LOOSE_ISO
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
//	if(0 < debug)std::cout<< "roccor SF"<<std::endl;
		double roc = 1.;
		switch(ch){
		case elnu:{break;}
		case munu:{if(MC){
			if(0. < gen_pt){
				std::cout<<"rocco 1"<<std::endl;
				roc = rc.kSpreadMC(Q,pt,eta,phi,gen_pt,0,0);
			}else{
				std::cout<<"rocco 2"<<std::endl;
				auto u = gRandom->Rndm();
				roc = rc. kSmearMC(Q,pt,eta,phi,nl,u,0,0);
			}
		}else{
				std::cout<<"rocco 3"<<std::endl;
				roc = rc. kScaleDT(Q,pt,eta,phi,0,0);
		}break;}}// not MC, case munu, switch
		if(0 < debug) std::cout << "roc " << roc << std::endl;
		return roc;
	};
}
/*
template<typename T>// allow us to return w/o knowing data type
auto retVar(const T& v){return[&](){return v;};}
*/
// Jet Energy Resolution and Jet Energy Smearing
auto jet_smear_pt_resol(
	 const floats &pt,const floats &eta
	,const float  rho
){
	if(0<debug) std::cout<<"jet smear pt resol"<<std::endl;
	doubles resol(pt.size());
	if(!all_equal(pt.size(),eta.size())) throw std::logic_error(
		"Collections must be the same size (jet_smear_pt_resol)");
	if(pt.empty()) throw std::logic_error(
		"Collections must not be empty in  (jet_smear_pt_resol)");
	const  doubles ptD = static_cast<doubles>(pt);
	double etaMin,etaMax,rhoMin,rhoMax;
	double pt_Min,pt_Max;
	double a,b,c,d;
	char    filepath[  ] = "aux/Fall17_V3_MC_PtResolution_AK4PFchs.txt";
	io::CSVReader<10> thisCSVfile(filepath);
	thisCSVfile.read_header(io::ignore_extra_column,
	"etaMin","etaMax","rhoMin","rhoMax","ptMin","ptMax","a","b","c","d");
	while(thisCSVfile.read_row(
	 etaMin , etaMax , rhoMin , rhoMax ,pt_Min ,pt_Max , a , b , c , d)){
		// The following if for if stacks correctly
		if(rhoMin < rho && rho < rhoMax)
		for(size_t i=0; i < pt.size() ;++i)
		if(etaMin < eta[i] && eta[i] < etaMax
		&& pt_Min < ptD[i] && ptD[i] < pt_Max){
			resol[i] += std::sqrt(  a * std::abs(a)/(ptD[i]*ptD[i])
			                    + b*b * std::pow(ptD[i],d) + c*c);
		}
	}// No need to close file after this while loop
	return resol;//};
}
auto jet_smear_Sjer(const floats &etas){
	//if(0<debug) std::cout<<"jet smear sjer"<<std::endl;
	doubles Sjers(etas.size());
	double  etaMin,etaMax;
	double  centralSF,dnSF,upSF;
	char    filepath[  ] = "aux/Fall17_V3_MC_SF_AK4PFchs.txt";
	io::CSVReader<5> thisCSVfile(filepath);
	thisCSVfile.read_header(io::ignore_extra_column,
	                          "etaMin","etaMax","centralSF","dnSF","upSF");
	while(thisCSVfile.read_row(etaMin , etaMax , centralSF , dnSF , upSF)){
		for(size_t i=0; i  <  etas.size() ;++i)
		   if(      etaMin <  etas[i] && etas[i] < etaMax)
		          Sjers[i] += centralSF;
	}
	return Sjers;//};
}
constexpr auto ramp(const double Sjer){
	if(0. < Sjer) return Sjer;
	else          return   0.;
}
auto delta_R_jet_smear(
		 const floats &jpt ,const floats &jeta,const floats &jphi
		,const   ints &jma  // jet gen match
		,const floats &gpt ,const floats &geta,const floats &gphi
//		,const   ints &mask,// mask is tight jets
		,const float   rho
	){
	if(0<debug) std::cout<<"delta r jet smear"<<std::endl;
	const size_t size = jpt.size();
	doubles cjers; cjers.reserve(/*tJ_final_*/size);
	if(!all_equal(    size  ,jeta.size(),jphi.size(),jma.size())//,mask.size())
	|| !all_equal(gpt.size(),geta.size(),gphi.size()))
		throw std::logic_error(
			"Collections must be the same size (deltaR_Jsmear)");
	if(!size) throw std::logic_error(
			"Collections must not be empty for (deltaR_Jsmear)");
	// the method used in here is the Jets Smearing Hybrid Method
	double temp,RconeBy2 = ak4RconeBy2;
	auto resol = jet_smear_pt_resol/*(ptRcsv)*/(jpt,jeta,rho);
	auto  Sjer = jet_smear_Sjer   /*(sjerCsv)*/(    jeta    );
	for(size_t j=0; j < size ;++j){// j goes with jet, i goes with gen
	    int      i ;// we will store the index from jma here(see below)
//	if(    0 !=      mask[j]){// only do stuff if tight jet
	   if(-1 != (i = jma [j])// -1 == jma when gen level info is absent
	   && std::abs(  jpt [j] - gpt [i]) <  3.*resol[j]*jpt [j]
	   && deltaR(    geta[i],  gphi[i],jeta[j],jphi[j]) < RconeBy2){
		   temp = (1.+(1.+Sjer[j])// Scaling method
		             *(1.-gpt [i]/jpt[j]));
	   }else{// Stochastic smearing
		   double Normdist = gRandom->Gaus(0.,Sjer[j]);
		   double  max_val =        Sjer[j] * Sjer[j] - 1.;
		   temp = (1.+Normdist*std::sqrt(ramp(max_val)));
	   }
	   if(temp < 0.) temp = 0.;
	   cjers.emplace_back(temp);
	}
	return cjers;//};
}
auto metCjer(
	 const  floats &jpt ,const floats &jeta
	,const  floats &jphi,const floats &jmas
	,const doubles &cjer
){
	if(0<debug) std::cout<<"MET Correction"<<std::endl;
	if(!all_equal(jpt .size(),jeta.size(),jphi.size(),
	              jmas.size(),cjer.size()))throw std::logic_error(
		"Collections must be the same size (MET correction)");
	if(jpt.empty()) throw std::logic_error(
		"Collections must not be empty for (MET correction)");
	ROOT::Math::PxPyPzMVector acc;// initial corrections all zero
	for(size_t i=0; i  < jpt.size() ;++i){
		double omc = 1. - cjer[i];
		ROOT::Math::PtEtaPhiMVector
		jet(omc*jpt[i],jeta[i],jphi[i],omc*jmas[i]);
		acc +=  jet;// same as unsmeared jp - smeared jp
	}
	return acc;
}
auto metCorrection(
	 ROOT::Math::PxPyPzMVector  &cTJer
//	,ROOT::Math::PxPyPzMVector  &cFJer
	,const float mpt ,const float mphi,const float mSEt
){
	ROOT::Math::PxPyPzEVector
//	corr((cTJer+cFJer).Px(),(cTJer+cFJer).Py(),0,(cTJer+cFJer).Et());
	corr( cTJer       .Px(), cTJer       .Py(),0, cTJer       .Et());
	ROOT::Math::PtEtaPhiEVector cmet(mpt,0.,mphi,mSEt);
	cmet += corr;
	return  cmet;
}
// btag suPer set and btag suB set
inline auto btagP(const doubles &eta){return abs(eta) < BJET_ETA_MAX;}
inline auto btagB(const ints  &btagP,const floats &btags){
	return   btagP && ( BTAG_DISC_MIN < btags );// all tJ length
}
inline auto	isBquark(   const ints &id , const ints &mask ){
	return mask &&   (   5 ==  abs(  id ) );
}
inline auto	isUDSCgluon(const ints &id , const ints &mask ){
	ints aid = abs ( id );
	return mask && ( 21 == aid || /*0 < aid &&*/ aid < 5 );// EXCLUDE 5
}
auto btagCSVv2(const bool check_CSVv2){
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
	io::CSVReader<11> thisCSVfile("aux/CSVv2_94XSF_V2_B_F.csv");
	thisCSVfile.next_line();// we happen to not need the header line
	// The following nests too much, so we do not indent
	// Each blank line means nesting deeper
	while(thisCSVfile.read_row(CSVv2,measureType,sysType,jetFlav,
	      etaMin,etaMax,pt_Min,pt_Max,CSVmin,CSVmax,rawFormula)){

	if(check_CSVv2){
		b= BTAG_DISC_MIN <= CSVv2
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
	formulae[i] = tempFormula;
	}}}}}
	}// No need to close file after this while loop.
	// resume indentation
	for(size_t j=0; j < formulae.size() ;++j){
		Eval ev;// numbers calculator
		results[j] = ev.eval(const_cast<char*>(
		             formulae[j].c_str())).real();
	}
	if(0<debug)std::cout<<"btagCSVv2 exiting"<<std::endl;
	return results;};
}
inline auto met_pt_cut(const channel ch){
	float lower_bound;
	switch(ch){
		case elnu:{lower_bound = MET_EL_PT;break;}
		case munu:{lower_bound = MET_MU_PT;break;}
		default  :throw std::invalid_argument(
			"Unimplemented ch (met_pt_cut)");
	}
	return [=](const float   met_pt)->bool
	   {return lower_bound < met_pt;};
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
auto lep_nu_invmass(
	 const double lep_pt
	,const double lep_eta
	,const double lep_phi
	,const double lep_mass
	,const double cal_metpt// TODO: If NOT cal, rename
	,const double cal_metphi
//	,const float  cal_metEt// cal_metEt is unused
){
	// this function computes the invariant mass of charged lepton
	// and neutrino system, in order to calculate the W mass later on.
	ROOT::Math::PtEtaPhiMVector
	lep(   lep_pt,lep_eta,   lep_phi,lep_mass),
	neu(cal_metpt,lep_eta,cal_metphi,   0.   );
	return (lep+neu).M();
}
inline auto easy_mass_cut(const double theo,const double cut){
               return [=](const double ours){return std::abs(ours-theo)<cut;};
}
constexpr auto abs_deltaphi(const double Zphi,const double Wphi)
	{return std::abs(delta_phi(Zphi-Wphi));}

inline auto jet_deltaphi(const doubles& phis){
	if(0<debug) std::cout<<"jet deltaphi"<<std::endl;
	doubles deltaphis;// half of anti-symmetric matrix has n(n-1)/2 elements
	deltaphis.reserve((phis.size()*(phis.size()-1))/2);// reserving size
	// The following two for loops stack correctly     // yet leave
	for(size_t   i=0; i < phis.size()-1 ;++i)          // empty.
	for(size_t j=i+1; j < phis.size()   ;++j)
		deltaphis.emplace_back(abs_deltaphi(phis[i],phis[j]));
	return deltaphis;
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
auto find_z_pair(
	 const doubles& pts
	,const doubles& etas
	,const doubles& phis
	,const doubles& ms
	,const    ints& lead_bjet
){
	if(0<debug) std::cout<<"find z pair"<<std::endl;
	// This function finds the pair nearest to z mass
	double  z_reco_mass = std::numeric_limits<double>::infinity();
	size_t  jet_index_1 = std::numeric_limits<size_t>::max();
	size_t  jet_index_2 = std::numeric_limits<size_t>::max();
	const size_t  njets = pts.size();
	if(!all_equal(njets,etas.size(),phis.size(),ms.size()))
		throw std::logic_error(
		"Collections must be the same size in Z-pair");
	if(pts.size()==0)	throw std::logic_error(
		"Collections must not be empty in Z-pair");
	ints z_pair(njets, 0);// vector of zeroes
	// The next two for loops stack correctly with the if
	for(size_t   i=0; i < njets-1 ;++i)// only try to z pair if
	for(size_t j=i+1; j < njets   ;++j)// both are NOT lead bjet
	if(lead_bjet[i] == 0 && lead_bjet[j] == 0){
		ROOT::Math::PtEtaPhiMVector
		jet1(pts[i],etas[i],phis[i],ms[i]),
		jet2(pts[j],etas[j],phis[j],ms[j]);
		if (const double reco_mass = (jet1+jet2).M();
		std::abs(Z_MASS- reco_mass) < std::abs(Z_MASS-z_reco_mass)){
		   z_reco_mass = reco_mass;// found nearer pair to z mass
		   jet_index_1 = i;
		   jet_index_2 = j;
		}
	}
	z_pair[jet_index_1] = 1;
	z_pair[jet_index_2] = 1;
	//if(1<debug) std::cout<<"z pair"<<z_pair<<std::endl;
	return z_pair;
}
auto LVpairAdd(
	 const doubles& pt__pair// Create LorentzV from 2 jets 4-mom
	,const doubles& eta_pair
	,const doubles& phi_pair
	,const doubles& mas_pair
){
	//if(0<debug) std::cout<<"LVpairAdd"<<std::endl;
	if(2 != pt__pair.size()) throw std::logic_error(
		"Not pair of Z (LVpairAdd)");
	ROOT::Math::PtEtaPhiMVector
	v(pt__pair[0],eta_pair[0],phi_pair[0],mas_pair[0]),
	p(pt__pair[1],eta_pair[1],phi_pair[1],mas_pair[1]);
	return v+p;
}
/*
auto deltaphi_cut(const double    x){
      return  [=](const doubles& dps){
		return  any(dps >= x);
	};
}
*/
auto top_reconst(
	 const doubles& bjets__pt
	,const doubles& bjets_eta
	,const doubles& bjets_phi
	,const doubles& bjets_mas
	,const double   wpair__pt// TODO: rename
	,const double   wpair_eta
	,const double   wpair_phi
	,const double   wpair_mas
){
	if(0<debug) std::cout<<"top recnst"<<std::endl;
	// This function finds the closest to top mass
	double     cur_mass = std::numeric_limits<double>::infinity();
	const size_t  nbjets=bjets__pt.size();
	if(!all_equal(nbjets,bjets_eta.size(),bjets_phi.size(),bjets_mas.size()))
		throw std::logic_error(
		"Collections must be the same size (top_recnst)");
	if(nbjets == 0) throw std::logic_error(
		"Collections must not be empty in  (top_recnst)");
	ROOT::Math::PtEtaPhiMVector reco_top,
	   RecoW(wpair__pt   ,wpair_eta   ,wpair_phi   ,wpair_mas   );
	// The following if for stacks correctly
	if(std::abs(RecoW.M() - W_MASS) < W_MASS_CUT)
	for(size_t i=0; i < nbjets ;++i){
	   ROOT::Math::PtEtaPhiMVector
	   temp (bjets__pt[i],bjets_eta[i],bjets_phi[i],bjets_mas[i]);
	   temp += RecoW;
	   if(double reco_mass = temp.M();
	     std::abs(TOP_MASS-reco_mass) < std::abs(TOP_MASS-cur_mass)){
		    cur_mass = reco_mass;// found closer to top mass
		   reco_top  = temp;
	   }
	}
	return reco_top;
}
template <typename T>
auto allReconstruction(T &rdf){
	return rdf
//	.Filter(met_pt_cut(ch),{"cmet__pt"},"MET Pt cut")// TODO: Re-enable!
	.Define("lep__pt","lepB_pt * roccorSF")// only pt and mass scales
	.Define("lep_mas","lepBmas * roccorSF")// eta and phi untouched
	// now for transverse W; lepton done
	. Alias("tw_lep__pt","lep__pt")
	. Alias("tw_lep_eta","lep_eta")
	. Alias("tw_lep_phi","lep_phi")
	.Define("tw_lep_mas",transverse_w_mass,
	       {"tw_lep__pt",
	        "tw_lep_phi","cmet__pt","cmet_phi"})
	.Define(   "lep_nu_invmass",lep_nu_invmass,
	       {   "lep__pt",
	           "lep_eta",
	           "lep_phi",
	           "lep_mas",
	          "cmet__pt" ,
	          "cmet_phi"})
//	       "CaloMET_sumEt"})// TODO: add this back
//	.Filter(easy_mass_cut(W_MASS,W_MASS_CUT),{"tw_lep_mas"},"W mass cut")
	// reconstruct z avoiding bjets first
	.Define(   "fin_jets_Dph" , jet_deltaphi  ,{"fin_jets_phi"})
	.Define( "lead_bjet"      , find_lead_mask,{"fin_jets__pt","btagB"})
	.Define(      "bjet__pt"  ,"fin_jets__pt[lead_bjet]")// Leading bjet
	.Define(      "bjet_eta"  ,"fin_jets_eta[lead_bjet]")// 4-momentum
	.Define(      "bjet_phi"  ,"fin_jets_phi[lead_bjet]")// used for top
	.Define(      "bjet_mas"  ,"fin_jets_mas[lead_bjet]")// reconstruction
	.Define("z_reco_jets"     , find_z_pair   ,
	       { "fin_jets__pt"   ,
	         "fin_jets_eta"   ,
	         "fin_jets_phi"   ,//easier to push back
	         "fin_jets_mas"   , "lead_bjet" })
	.Define(   "z_pair__pt"   ,   "fin_jets__pt[z_reco_jets]")
	.Define(   "z_pair_eta"   ,   "fin_jets_eta[z_reco_jets]")
	.Define(   "z_pair_phi"   ,   "fin_jets_phi[z_reco_jets]")
	.Define(   "z_pair_mas"   ,   "fin_jets_mas[z_reco_jets]")
	.Define(   "z_LV"         , LVpairAdd    ,
	       {   "z_pair__pt"   ,
	           "z_pair_eta"   ,
	           "z_pair_phi"   ,
	           "z_pair_mas"  })
	.Define(   "z_mas"        ,"z_LV. M ()")
//	.Filter( easy_mass_cut(Z_MASS,Z_MASS_CUT),{"z_mas"},"z mass cut")
	.Define(   "z__pt"        ,"z_LV.Pt ()")
	.Define(   "z_eta"        ,"z_LV.Eta()")
	.Define(   "z_phi"        ,"z_LV.Phi()")
	.Define(   "z_lep_min_dR" , jet_lep_min_deltaR,// TODO: Check input?
	       {   "z_pair_eta"   ,
	           "z_pair_phi"   ,
	              "lep_eta"   ,
	              "lep_phi"  })// TODO: wrong input below
//	.Filter(      deltaR_z_l,  {  "z_lep_min_dR"},"delta R ZL")
	.Define(   "zw_Dph" , abs_deltaphi,{"z_phi","tw_lep_phi"})
//	.Filter(      deltaphi_cut(DELTA_PHI_ZW),
//	       {   "zw_Dph"},"delta phi ZW cut")
//	.Define("rawMET_phi","static_cast<double>(MET_phi)")
	.Define( "zmet_Dph" , abs_deltaphi,{"z_phi","cmet_phi"})
//	.Filter(      deltaphi_cut(DELTA_PHI_ZMET),
//	       { "zmet_Dph"},"Z met cut ")
	.Define("z_pair_Dph",[](doubles p){return
	                      abs_deltaphi(p[0],p[1]);},{"z_pair_phi"})
	// reconstruct transverse top follows; z reconstructed
	.Define("recoTtop",top_reconst,
	       {"bjet__pt",
	        "bjet_eta",
	        "bjet_phi",
	        "bjet_mas",
	      "tw_lep__pt",
	      "tw_lep_eta",
	      "tw_lep_phi",
	      "tw_lep_mas"})
	.Define("ttop__pt","recoTtop.Pt ()")
	.Define("ttop_eta","recoTtop.Eta()")
	.Define("ttop_phi","recoTtop.Phi()")
	.Define("ttop_mas","recoTtop. M ()")
	.Filter("ttop_mas > 1.","tTm tiny mass filter")
	;
}
// Btagging for eff i and eff j
auto BTaggedEffGiver(TH2D*& ratio,bool b){
return [&,b](const doubles& pt,const doubles& eta){
	if(0<debug) std::cout<<"bt eff giver in  "<< ratio <<std::endl;
	if(!all_equal(pt.size(),eta.size())) throw std::logic_error(
		"Collections must be the same size (effGiver)");
	if(pt.empty()) throw std::logic_error(
		"Collections must not be empty in  (effGiver)");
	doubles BTaggedEff(pt.size(),b ? 1. : 0.);// TODO: testing min bad
	for(size_t   i=0; i <  pt.size() ;++i){
		int  PtBin = ratio->GetXaxis()->FindBin(         pt [i] );
		int EtaBin = ratio->GetYaxis()->FindBin(std::abs(eta[i]));
		double eff = ratio->GetBinContent(PtBin,EtaBin);
		// NOTE: The below comments are outdated. Now 0 <= e <= 1
//		if(FP_NORMAL == std::fpclassify(eff))// if eff non-zero/inf/NaN
			BTaggedEff[i] = eff;
		// above only pushed back nonzero nice eff
		// what do we do with eff==0? check with kathryn
		// below, if ej, and near 1., we put 0. instead
//		if(!b && std::abs(eff-1.) <= 2*std::numeric_limits<double>::epsilon())
//			BTaggedEff[i] = 0.;
	}
	if(5<debug) std::cout<<"bt eff giver "<< BTaggedEff <<" for b "<<b<<std::endl;
	return BTaggedEff;};
}
inline auto product(const doubles &sfi){
	return std::accumulate(sfi.cbegin(),sfi.cend(),1.,std::multiplies<>());
}
/*
inline auto IsEffBTaggedProduct(const doubles& IsEffBTagged){
//	double result = std::reduce(//std::execution::par_unseq,
//		IsEffBTagged. cbegin(),// parallel multiply
//		IsEffBTagged.   cend(),// un-sequentially
//		1.,std::multiplies<>()
//	);
	double result = product(IsEffBTagged);
	if(20<debug) std::cout<<"    IsEffBTaggedProduct "<<result<<std::endl;
	return result;
}
*/
inline auto NoEffBTaggedProduct(const doubles& NoEffBTagged){
	doubles values = 1. - NoEffBTagged;
//	double  result = std::reduce(//std::execution::par_unseq,
//	        values.cbegin(),// parallel multiply
//	        values.  cend(),// un-sequentially
//	        1.,std::multiplies<>()
//	);
	double  result = product(values);
	if(5<debug)  std::cout<<"    NoEffBTaggedProduct "<<result<<std::endl;
	return  result;
}
/*
inline auto Sfi_IsEffBTaggedProduct(const doubles& IsEffBTagged,const doubles& sfi){
	if(IsEffBTagged.size() != sfi.size()) throw std::logic_error(
		"Sfi_IsEffBTaggedProduct got diff sizes");// both tJ len
	doubles values = sfi * IsEffBTagged;
	double  result = std::accumulate(
	        values.cbegin(),
	        values.  cend(),
	        1.,std::multiplies<>()
	);
//	double  result = std::transform_reduce(//std::execution::par_unseq,
//		         sfi. cbegin(),// parallel multiply
//		         sfi.   cend(),// un-sequentially
//		IsEffBTagged. cbegin(),// stops correctly
//		1.,std::multiplies<>(),std::multiplies<>()
//	);
	if(20<debug) std::cout<<"Sfi_IsEffBTaggedProduct "<<result<<std::endl;
	return result;
}
*/
inline auto Sfj_NoEffBTaggedProduct(const doubles &NoEffBTagged,const doubles &sfj){
	if(NoEffBTagged.size() != sfj.size()) throw std::logic_error(
		"Sfj_NoEffBTaggedProduct got diff sizes");// both tJ len
	doubles values = 1. - sfj * NoEffBTagged;
	double  result = product(values);
/*	double  result = std::transform_reduce(//std::execution::par_unseq,
		         sfj. cbegin(),// parallel multiply
		         sfj.   cend(),// un-sequentially
		NoEffBTagged. cbegin(),// stops correctly
		1.,std::multiplies<>(),
//		[](double x,double y){return -std::fma(x,y,-1.);}
		[](double x,double y){return 1. - x*y;}
	);
*/
	if(5<debug)  std::cout<<"Sfj_EffNoBTaggedProduct "<<result<<std::endl;
	return result;
}
/*
constexpr auto btag_weight(const double p_data,const double p_MC){
	double weight = p_data / p_MC;
//	if(FP_NORMAL != std::fpclassify(weight))weight=1.;// Rids non-zero/inf/NaN
//	if(0 < debug)   std::cout<<"btag_w is   "<<weight<<std::endl;
//	if(weight < 0)  std::cout<<"btag_w lt 0 "        <<std::endl;
//	if(9999<weight){std::cout<<"btag_w huge "<<weight
//	<<" huge;replaced"<<std::endl; weight = 1.;}
//	<<std::endl;}
	return  weight;
}
*/
auto top_pt_sf(const dataSource ds){
    return [=](const ints& gId, const ints& flag, const floats& pt){
	// gId = GenPart_pdgId, flag = GenPart_statusFlags, pt = GenPart_pt
	// This reweighing is only computed for the regions where ttbar dominates
	// the signal, so by finding those regions, i.e (pt_min--pt_max)
	// only for the ttbar sample.
	// Top pt reweighing is computed using MC truth and applied as weight across
	// the analysis for the ttbar sample.
	if(ttb != ds) return 1.;
	if(0<debug)std::cout<<"top pt sf"<<std::endl;
	if(!all_equal(gId.size(),pt.size())) throw std::logic_error(
		"Collections must be the same size (top_pt_reweigh)");
	if(gId.empty()) throw std::logic_error(
		"Collections must not be empty for (top_pt_reweigh)");
	ints idxs = (6 == abs(gId) && 13 == flag);// 13 == lastCopy, 6 == top
	auto  pts =  pt[idxs];
	pts = pts[0.<pts && pts < 500.];// TODO: new constexpr ptMin ptMax
	double sf = 1.;
	if(0 < pts.size()) sf = std::sqrt(product(Map(pts,[](float p){
		return exp(.0615 - .0005*static_cast<double>(p));})));
	if(0<debug) std::cout<<"top pt sf is "<<sf<<std::endl;
	return sf;};
}
// Lepton efficiencies
auto elEffGiver(
	 const double             pt
	,const double             eta
	,const TH2F* const &recoLowEt
	,const TH2F* const &reco_pass
	,const TH2F* const &tight_94x
){
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
	return id * iso * eff * smr;};
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
// Simulation correction Scale Factors
inline auto sf(
	 const  dataSource    ds
	,const TH1D* const &PuWd
	,const TH1D* const &PuUd
	,const TH1D* const &PuDd
){
	if(0<debug) std::cout<<"scale factor "<<std::endl;
	return [=](
		 const double b
		,const double mostSF
		,const    int npv
	){
		// TODO: trigger efficiency
		double result;
		bool MC = true;
		switch(ds){
			case tzq:{result =    TZQ_W;break;}
			case  ww:{result = WWLNQQ_W;break;}
			case  wz:{result = WZLNQQ_W;break;}
			case  zz:{result = ZZLLQQ_W;break;}
			case ttb:{result =  TTBLV_W;break;}
			case ttz:{result =  TTZQQ_W;break;}
			case met:// fall through to cms
			case cms:{result = 1.;MC=false;break;}// ignore btag wt
//			default :throw std::invalid_argument(
//				"Unimplemented ds (sf)");
		}
		if(MC) result *= pile(PuWd,PuUd,PuDd)(npv)[puW];
		if(result < 0)std::cout<<"pile lt 0"<<std::endl;
		if(5<debug)   std::cout<<"b_w "<< b
		<<" few_SF " << result//<<std::endl;
		<<" mostSF " << mostSF  <<std::endl;
		result *=   b * mostSF;
		if(result < 0)std::cout<<"final sf lt 0 "<<result<<std::endl;
		return result;
	};
}
inline auto rep_const(const double sf,const doubles& iRVec){
	// this function just repeats sf, for the size of iRVec
	if(0<debug) std::cout<<"repeat const "<<iRVec.size()<<std::endl;
	doubles weight(iRVec.size(),sf);
	return  weight;
}
template <typename T>
auto finalScaling(
	 const dataSource     ds
	,const TH1D* const &PuWd
	,const TH1D* const &PuUd
	,const TH1D* const &PuDd
	,T &rdf
){
	return rdf
	.Define("sf",sf(ds,PuWd,PuUd,PuDd),{"btag_w","mostSF","PV_npvs"})
	. Alias("nw_lep__pt"       ,"sf")// is just one value, == sf
	. Alias("nw_lep_eta"       ,"sf")
	. Alias("nw_lep_phi"       ,"sf")
	. Alias("nw_lep_mas"       ,"sf")
	.Define("nw_fin_jets__pt"  ,rep_const,{"sf","fin_jets__pt"  })
	.Define("nw_fin_jets_eta"  ,rep_const,{"sf","fin_jets_eta"  })
	.Define("nw_fin_jets_phi"  ,rep_const,{"sf","fin_jets_phi"  })
	.Define("nw_fin_jets_mas"  ,rep_const,{"sf","fin_jets_mas"  })
	.Define("nw_bjet__pt"	   ,rep_const,{"sf","bjet__pt"      })
	.Define("nw_bjet_eta"      ,rep_const,{"sf","bjet_eta"      })
	.Define("nw_bjet_phi"      ,rep_const,{"sf","bjet_phi"      })
	.Define("nw_bjet_mas"      ,rep_const,{"sf","bjet_mas"      })
	.Define("nw_jet_lep_min_dR",rep_const,{"sf","jet_lep_min_dR"})
	.Define(  "nw_z_lep_min_dR",rep_const,{"sf",  "z_lep_min_dR"})
	. Alias( "nw_tw_lep_mas"   ,"sf")
	. Alias(  "nw_z_mas"       ,"sf")
	.Define("nw_fin_jets_Dph"  ,rep_const,{"sf","fin_jets_Dph"})
	.Define(  "nw_z_pair_Dph"	,"sf")
	. Alias(    "nw_zmet_Dph"  ,"sf")
	. Alias(      "nw_zw_Dph"  ,"sf")
	. Alias("nw_lep_nu_invmass","sf")
	. Alias("nw_ttop__pt"      ,"sf")
	. Alias("nw_ttop_mas"      ,"sf")
	;
}
auto runLBfilter(
	const std::map<size_t,std::vector<std::pair<size_t,size_t>>>
	&runLBdict
){
	return [&](const unsigned int run,const unsigned int LB){
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
}// namespace
void calchisto(const channel ch,const dataSource ds){
	ROOT::EnableImplicitMT(4);// SYNC WITH CONDOR JOBS!
	// Open LB file even if Monte Carlo will NOT use it
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
	// pile up
	tF = TFile::Open("aux/pileupMC.root");// denom because wrong
	tF ->GetObject("pileup",t1d);t1d->SetDirectory(nullptr);
	t1d->Scale(1.0/t1d->Integral());
	const TH1D* const PuDe = static_cast<TH1D*>(t1d);
	tF ->Close();
	tF = TFile::Open("aux/truePileupTest.root");
	tF ->GetObject("pileup",t1d);t1d->SetDirectory(nullptr);
	t1d->Scale(1.0/t1d->Integral());
	t1d->Divide(PuDe);
	const TH1D* const PuWd = static_cast<TH1D*>(t1d);
	tF ->Close();
	tF = TFile::Open("aux/truePileupUp.root");
	tF ->GetObject("pileup",t1d);t1d->SetDirectory(nullptr);
	t1d->Scale(1.0/t1d->Integral());
	t1d->Divide(PuDe);
	const TH1D* const PuUd = static_cast<TH1D*>(t1d);
	tF ->Close();
	tF = TFile::Open("aux/truePileupDown.root");
	tF ->GetObject("pileup",t1d);t1d->SetDirectory(nullptr);
	t1d->Scale(1.0/t1d->Integral());
	t1d->Divide(PuDe);
	const TH1D* const PuDd = static_cast<TH1D*>(t1d);
	tF ->Close();
	tF	= nullptr; tHf = nullptr; tHd = nullptr; t1d = nullptr;
	RoccoR rc("src/roccor.Run2.v3/RoccoR2017.txt");
//	std::cout<<"Auxiliary files processed"       <<std::endl;

	// Open data files even if unused
	// then automatically choose which one to read from
	// No penalty for opening and leaving unused
	// Can even open multiple times at once in parallel
	// Open MC data source EVEN IF UNUSED
	std::string temp_header="/data/disk0/nanoAOD_2017/",
	temp_opener,temp_footer="/*.root";/**/
	switch(ds){// CMS and MET MUST do some OPENABLE file ; reject later
	case tzq:{temp_opener="/data/disk3/nanoAOD_2017/tZqlvqq/*.root"  ;break;}/**/
	case  ww:{temp_opener=temp_header+ "WWToLNuQQ"       +temp_footer;break;}
	case  wz:{temp_opener=temp_header+ "WZTo1L1Nu2Q"     +temp_footer;break;}
	case  zz:{temp_opener=temp_header+ "ZZTo2L2Q"        +temp_footer;break;}
	case ttb:{temp_opener=temp_header+ "TTToSemileptonic"+temp_footer;break;}
	case ttz:{temp_opener=            "ttz_dir"          +temp_footer;break;}
	case met:{temp_opener=temp_header+"ttZToQQ_ext"      +temp_footer;break;}
	case cms:{temp_opener=temp_header+"ttZToQQ"          +temp_footer;break;}
//	default :throw std::invalid_argument("Unimplemented ds (rdfopen)");
	}
	ROOT::RDataFrame mc__df("Events",temp_opener);// Monte Carlo
	// Open chains of exptData EVEN IF UNUSED
	TChain elnuCMS("Events");
	TChain munuCMS("Events");
	TChain bothMET("Events");
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
	temp_header="/data/disk0/nanoAOD_2017/METRun2017";
	for(std::string c:{"B","C","D","E","F"}){// guaranteed sequential
		temp_opener=temp_header+ c +temp_footer;
		bothMET.Add(temp_opener.c_str());
	}
	ROOT::RDataFrame  elnudf(elnuCMS);
	ROOT::RDataFrame  munudf(munuCMS);
	ROOT::RDataFrame  bothdf(bothMET);
	const bool MC = !(met == ds || cms == ds);
	auto df = [&,ch,ds](){// Get correct data frame
		switch(ds){
			case tzq:
			case  ww:// fall through!
			case  wz:
			case  zz:
			case ttb:
			case ttz:{           return mc__df;break;}
			case met:{           return bothdf;break;}
			case cms:{switch(ch){// MC is already false
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
//		default  :throw std::invalid_argument(
//			"Unimplemented ch (init)");
	}
	// make test runs faster by restriction. Real run should not
//	auto dfr = df.Range(1000000);// remember to enable MT when NOT range
	auto init_selection = df// remove one letter to do all
/*	.Filter(triggers(ch),
		{ "HLT_Ele32_WPTight_Gsf_L1DoubleEG"//"HLT_Ele28_eta2p1_WPTight_Gsf_HT150"
		 ,"HLT_Ele28_eta2p1_WPTight_Gsf_HT150"//"HLT_PFMET120_PFMHT120_IDTight"
		 ,"HLT_PFJet15"
		 ,"HLT_IsoMu27"
		},"Triggers Filter")*/
	// lepton selection first
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
	.Filter(lep_tight_cut(ch),{"loose_leps",
	        "Electron_cutBased",// edit function for  tight -> loose
	        "Muon_pfRelIso04_all"},"lepton cut")// left with 1 tight lepton
	// B for Bare, Before muon RoccoR; only pt mass scales, eta phi untouched
	.Define("lepB_pt","static_cast<double>("+temp_header+ "pt [loose_leps][0])")
	.Define("lep_eta","static_cast<double>("+temp_header+ "eta[loose_leps][0])")
	.Define("lep_phi","static_cast<double>("+temp_header+ "phi[loose_leps][0])")
	.Define("lepBmas","static_cast<double>("+temp_header+"mass[loose_leps][0])")
	.Define("lep___q",temp_header+"charge[loose_leps][0]")// int, only for RoccoR
	// jets selection follows; lepton selected
	.Define("rawJet_eta","static_cast<ROOT::RVec<double>>(Jet_eta )")
	.Define("rawJet_phi","static_cast<ROOT::RVec<double>>(Jet_phi )")
	.Define("jet_lep_min_dR"   ,jet_lep_min_deltaR,// later reused with doubles
	       {"rawJet_eta","rawJet_phi","lep_eta","lep_phi"})// gcc fail template
	.Define("tight_jets"       ,tight_jet_id,
	       {"jet_lep_min_dR"   ,"Jet_pt","Jet_eta","Jet_jetId"})
	.Filter( jetCutter(JETS_MIN,JETS_MAX),{"tight_jets" },"Jet cut")
/*	.Define("tight_jets__pt"   ,"rawJet__pt[tight_jets]")
	.Define("tight_jets_eta"   ,"rawJet_eta[tight_jets]")
	.Define("tight_jets_phi"   ,"rawJet_phi[tight_jets]")
	.Define("tight_jets_mas"   ,"rawJet_mas[tight_jets]")
*/	.Define("tJ_btagCSVv2"  ,"Jet_btagCSVV2[tight_jets]")// leave as floats
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
	temp_footer = "pt vs eta in " + temp_footer + " channel for ";
	switch(ds){
		case tzq:{temp_header+="tzq";temp_footer+="tZq";break;}
		case  ww:{temp_header+="_ww";temp_footer+=" WW";break;}
		case  wz:{temp_header+="_wz";temp_footer+=" WZ";break;}
		case  zz:{temp_header+="_zz";temp_footer+=" ZZ";break;}
		case ttb:{temp_header+="ttb";temp_footer+="ttb";break;}
		case ttz:{temp_header+="ttz";temp_footer+="ttZ";break;}
		case met:{temp_header+="met";temp_footer+="MET";break;}
		case cms:{temp_header+="cms";temp_footer+="CMS";break;}
//		default :throw std::invalid_argument(
//			"Unimplemented ds (hist titles)");
	}
	// Histogram names sorted, now branch into MC vs exptData
	if(MC){
	auto jecs_bjets// JEC == Jet Energy Correction, only for MC
	   = init_selection
	.Filter("Flag_goodVertices"
	    " || Flag_globalSuperTightHalo2016Filter"
	    " || Flag_HBHENoiseFilter"
	    " || Flag_HBHENoiseIsoFilter"
	    " || Flag_EcalDeadCellTriggerPrimitiveFilter"
	    " || Flag_BadPFMuonFilter"
	    " || Flag_BadChargedCandidateFilter"
	    " || Flag_ecalBadCalibFilter"// TODO: v2?
	    " || Flag_eeBadScFilter"
	       ,"Event Cleaning filter")
	// next 3 defines are just for muons, which need RoccoR
	.Define("lep_gpt",lep_gpt(ch),{"GenPart_pt","loose_leps"
	                              ,"Electron_genPartIdx"
	                              ,    "Muon_genPartIdx"})
	.Define("lep__nl",[=](ints L,ints m){if(munu==ch)return L[m][0];
	                                     else        return      0 ;},
	       {"Muon_nTrackerLayers","loose_leps"})
	.Define("roccorSF", roccorSF(rc,ch,MC)
	                  ,{"lepB_pt","lep_eta",
	                    "lep_phi","lep___q",
	                    "lep_gpt","lep__nl"})
	// make lep__pt and lep_mas using roccorSF in allReconstruction later
	.Define("cjer"     ,delta_R_jet_smear,
	       {   "Jet_pt",   "Jet_eta",   "Jet_phi","Jet_genJetIdx",
	        "GenJet_pt","GenJet_eta","GenJet_phi",//"tight_jets",
	        "fixedGridRhoFastjetAll"})
	.Define("cTJer" ,metCjer,
	       {"Jet_pt","Jet_eta","Jet_phi","Jet_mass","cjer"})
	.Define("CmetLV",metCorrection,{"cTJer"
	       ,"MET_pt"          ,"MET_phi","MET_sumEt"})
	.Define("cmet_dpx","cTJer .Px ()")
	.Define("cmet_dpy","cTJer .Py ()")
	.Define("cmet__pt","CmetLV.Pt ()")
	.Define("cmet_phi","CmetLV.Phi()")
	.Define("cmet_sEt","CmetLV. E ()")
	// fin = final; Monte Carlo needs JEC
	.Define("fin_jets__pt","                        (cjer * Jet_pt  )[tight_jets] ")
	.Define("fin_jets_eta","static_cast<ROOT::RVec<double>>(Jet_eta  [tight_jets])")
	.Define("fin_jets_phi","static_cast<ROOT::RVec<double>>(Jet_phi  [tight_jets])")
	.Define("fin_jets_mas","                        (cjer * Jet_mass)[tight_jets] ")
	// jets selected, now bjets and btagging preliminaries
	.Define("btagP"            ,btagP  ,{"fin_jets_eta"})// suPer vs suBset
	.Define("btagB"            ,btagB  ,{"btagP","tJ_btagCSVv2"})
	.Filter(jetCutter(BJETS_MIN,BJETS_MAX),{"btagB"},"b jet cut")
	.Define("tJpF"             ,"Jet_partonFlavour [tight_jets]")
	.Define("is_btag_numer"    ,isBquark   ,{"tJpF","btagB"})
	.Define("is_btag_denom"    ,isBquark   ,{"tJpF","btagP"})
	.Define("no_btag_numer"    ,isUDSCgluon,{"tJpF","btagB"})
	.Define("no_btag_denom"    ,isUDSCgluon,{"tJpF","btagP"})
	.Define("is_btag_numer__pt","    fin_jets__pt[is_btag_numer] ")
	.Define("is_btag_denom__pt","    fin_jets__pt[is_btag_denom] ")
	.Define("no_btag_numer__pt","    fin_jets__pt[no_btag_numer] ")
	.Define("no_btag_denom__pt","    fin_jets__pt[no_btag_denom] ")
	.Define("is_btag_numer_eta","abs(fin_jets_eta[is_btag_numer])")
	.Define("is_btag_denom_eta","abs(fin_jets_eta[is_btag_denom])")
	.Define("no_btag_numer_eta","abs(fin_jets_eta[no_btag_numer])")
	.Define("no_btag_denom_eta","abs(fin_jets_eta[no_btag_denom])")
	.Define("sfi",btagCSVv2( true),//,btagDF),// checks btag
	       {  "tJ_btagCSVv2","fin_jets__pt","fin_jets_eta","tJpF"})
	.Define("sfj",btagCSVv2(false),//,btagDF),// ignore btag
	       {  "tJ_btagCSVv2","fin_jets__pt","fin_jets_eta","tJpF"})
	;
	auto reco = allReconstruction(
	     jecs_bjets )
	;

	const std::vector<double> ptBins
	= {0.,20.,30.,40.,45.,50.,55.,60.,65.,75.,85.,100.,150.,500.};
	const int ptBinsSize = ptBins.size() - 1;
	auto
	h_no_btag_numer_PtVsEta =
	     reco.Histo2D({
	(        "no_numer_" + temp_header).c_str(),
	("MC no btag numer " + temp_footer).c_str(),
//	12,0,200,12,0,3},
	ptBinsSize,ptBins.data(),12,0,3},
	"no_btag_numer__pt",
	"no_btag_numer_eta")
	;
	h_no_btag_numer_PtVsEta->Draw("COLZ");
	h_no_btag_numer_PtVsEta->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	h_no_btag_numer_PtVsEta->GetYaxis()->SetTitle("PseudoRapidity #eta");

	auto
	h_no_btag_denom_PtVsEta =
	     reco.Histo2D({
	(        "no_denom_" + temp_header).c_str(),
	("MC no btag denom " + temp_footer).c_str(),
	ptBinsSize,ptBins.data(),12,0,3},
	"no_btag_denom__pt",
	"no_btag_denom_eta")
	;
	h_no_btag_denom_PtVsEta->Draw("COLZ");
	h_no_btag_denom_PtVsEta->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	h_no_btag_denom_PtVsEta->GetYaxis()->SetTitle("PseudoRapidity #eta");

	TH2D *
	no_btag_ratio = static_cast<TH2D*>(h_no_btag_numer_PtVsEta->Clone());
	no_btag_ratio->Divide(             h_no_btag_denom_PtVsEta.GetPtr());
	no_btag_ratio->Draw("COLZ");
	no_btag_ratio->SetNameTitle("ej", "no b tag ej");

	auto has_btag_eff
	   = reco
//	.Define("IsEffBTagged",BTaggedEffGiver(is_btag_ratio, true),
//	       {"fin_jets__pt","fin_jets_eta"})// TODO: check sensibility
	.Define("NoEffBTagged",BTaggedEffGiver(no_btag_ratio,false),
	       {"fin_jets__pt","fin_jets_eta"})// of this eff formula
//	.Define("P___ei" ,    IsEffBTaggedProduct,{"IsEffBTagged"})
//	.Define("P_sfei" ,Sfi_IsEffBTaggedProduct,{"IsEffBTagged","sfi"})
	.Define("P___ej" ,    NoEffBTaggedProduct,{"NoEffBTagged"})
	.Define("P_sfej" ,Sfj_NoEffBTaggedProduct,{"NoEffBTagged","sfj"})
//	.Define("P_MC"   ,"P___ei * P___ej")
//	.Define("P_Data" ,"P_sfei * P_sfej")
//	.Define("btag_w" ,btag_weight,{"P_Data","P_MC"})
	.Define("P_sf_i",product,{"sfi"})
	.Define("btag_w","P_sf_i * P_sfej / P___ej")
	.Define("ttbSF"  , top_pt_sf(ds),{"GenPart_pdgId",
	                                  "GenPart_statusFlags",
	                                  "GenPart_pt"})
	.Define("lepSF"  , lepEffGiver(ch,MC// not in reco because files
	                 , recoLowEt,reco_pass,tight_94x
	                 , id_N,id_Y,id_A,id_T
	                 , isoN,isoY,isoA,isoT
	                ),{"lep__pt","lep_eta"})
	.Define("mostSF" , "lepSF" /* ttbSF"*/)
	;
	auto finalDF = finalScaling(ds,PuWd,PuUd,PuDd,
	     has_btag_eff )
	;
	// Copied to below, skip MC-only, ADD MET_sumEt!
	// Assuming temp_header and footer and all are set per (hist titles)!
// MC only
		auto
	h_is_btag_numer_PtVsEta =
	     reco.Histo2D({
	(        "is_numer_" + temp_header).c_str(),
	("MC is btag numer " + temp_footer).c_str(),
	ptBinsSize,ptBins.data(),12,0,3},
	"is_btag_numer__pt",
	"is_btag_numer_eta")
	;
	h_is_btag_numer_PtVsEta->Draw("COLZ");
	h_is_btag_numer_PtVsEta->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	h_is_btag_numer_PtVsEta->GetYaxis()->SetTitle("PseudoRapidity #eta");

	auto
	h_is_btag_denom_PtVsEta =
	     reco.Histo2D({
	(        "is_denom_" + temp_header).c_str(),
	("MC is btag denom " + temp_footer).c_str(),
	ptBinsSize,ptBins.data(),12,0,3},
	"is_btag_denom__pt",
	"is_btag_denom_eta")
	;
	h_is_btag_denom_PtVsEta->Draw("COLZ");
	h_is_btag_denom_PtVsEta->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	h_is_btag_denom_PtVsEta->GetYaxis()->SetTitle("PseudoRapidity #eta");

	TH2D *
	is_btag_ratio = static_cast<TH2D*>(h_is_btag_numer_PtVsEta->Clone());
	is_btag_ratio->Divide(             h_is_btag_denom_PtVsEta.GetPtr());
	is_btag_ratio->Draw("COLZ");// trigger getting everything done
	is_btag_ratio->SetNameTitle("ei", "is b tag ei");

	auto h_sfi = finalDF.Histo1D({
	("sfi_"+temp_header).c_str(),
	("sfi "+temp_header).c_str(),
	5000,0,1.3},"sfi");
	h_sfi->GetXaxis()->SetTitle("sf_{i}");
	h_sfi->GetYaxis()->SetTitle("Event");
	h_sfi->SetLineStyle(kSolid);

	auto h_sfj = finalDF.Histo1D({
	("sfj_"+temp_header).c_str(),
	("sfj "+temp_header).c_str(),
	5000,0,1.3},"sfj");
	h_sfj->GetXaxis()->SetTitle("sf_{j}");
	h_sfj->GetYaxis()->SetTitle("Event");
	h_sfj->SetLineStyle(kSolid);
/*
	auto h_p_ei = finalDF.Histo1D({
	("p_ei_"+temp_header).c_str(),
	("p_ei "+temp_header).c_str(),
	5000,0,1.3},"P___ei");
	h_p_ei->GetXaxis()->SetTitle("\\prod_{i} e_{i}");
	h_p_ei->GetYaxis()->SetTitle("Event");
	h_p_ei->SetLineStyle(kSolid);
*/
	auto h_p_ej = finalDF.Histo1D({
	("p_ej_"+temp_header).c_str(),
	("p_ej "+temp_header).c_str(),
	5000,0,1.3},"P___ej");
	h_p_ej->GetXaxis()->SetTitle("\\prod_{j} 1 - e_{j}");
	h_p_ej->GetYaxis()->SetTitle("Event");
	h_p_ej->SetLineStyle(kSolid);

	auto h_p_sf_i = finalDF.Histo1D({
	("p_sf_i_"+temp_header).c_str(),
	("p_sf_i "+temp_header).c_str(),
	5000,0,1.3},"P_sf_i");
	h_p_sf_i->GetXaxis()->SetTitle("\\prod_{i} \\text{sf}_{i}");// e_{i}");
	h_p_sf_i->GetYaxis()->SetTitle("Event");
	h_p_sf_i->SetLineStyle(kSolid);

	auto h_p_sfej = finalDF.Histo1D({
	("p_sfej_"+temp_header).c_str(),
	("p_sfej "+temp_header).c_str(),
	5000,0,1.3},"P_sfej");
	h_p_sfej->GetXaxis()->SetTitle("\\prod_{j} 1 - \\text{sf}_{j} e_{j}");
	h_p_sfej->GetYaxis()->SetTitle("Event");
	h_p_sfej->SetLineStyle(kSolid);

	auto h_btag_w = finalDF.Histo1D({
	("btag_w_"+temp_header).c_str(),
	("btag_W "+temp_header).c_str(),
	5000,-1,1.3},"btag_w");
	h_btag_w->GetXaxis()->SetTitle("btag weight");
	h_btag_w->GetYaxis()->SetTitle("Event");
	h_btag_w->SetLineStyle(kSolid);

	auto h_cmet_sEt = finalDF.Histo1D({
	("cmet_sEt_"+temp_header).c_str(),
	("cmet_sEt "+temp_header).c_str(),
	100,0,600},"cmet_sEt");
	h_cmet_sEt->GetXaxis()->SetTitle("corrected MET Sum Et (GeV)");
	h_cmet_sEt->GetYaxis()->SetTitle("Event");
	h_cmet_sEt->SetLineStyle(kSolid);

	auto h_cmet__pt = finalDF.Histo1D({
	("cmet__pt_"+temp_header).c_str(),
	("cmet__pt "+temp_header).c_str(),
	50,0,300},"cmet__pt");
	h_cmet__pt->GetXaxis()->SetTitle("corrected MET p_{T} (GeV/c)");
	h_cmet__pt->GetYaxis()->SetTitle("Event");
	h_cmet__pt->SetLineStyle(kSolid);

	auto h_cmet_dpx = finalDF.Histo1D({
	("cmet_dpx_"+temp_header).c_str(),
	("cmet_dpx "+temp_header).c_str(),
	100,-300,300},"cmet_dpx");
	h_cmet_dpx->GetXaxis()->SetTitle("MET correction p_{x} (GeV/c)");
	h_cmet_dpx->GetYaxis()->SetTitle("Event");
	h_cmet_dpx->SetLineStyle(kSolid);

	auto h_cmet_dpy = finalDF.Histo1D({
	("cmet_dpy_"+temp_header).c_str(),
	("cmet_dpy "+temp_header).c_str(),
	100,-300,300},"cmet_dpy");
	h_cmet_dpy->GetXaxis()->SetTitle("MET correction p_{y} (GeV/c)");
	h_cmet_dpy->GetYaxis()->SetTitle("Event");
	h_cmet_dpy->SetLineStyle(kSolid);
// end MC only

	auto h_met_sEt = finalDF.Histo1D({
	("met_sEt_"+temp_header).c_str(),
	("met_sEt "+temp_header).c_str(),
	100,0,600},"MET_sumEt");
	h_met_sEt->GetXaxis()->SetTitle("MET Sum Et (GeV)");
	h_met_sEt->GetYaxis()->SetTitle("Event");
	h_met_sEt->SetLineStyle(kSolid);

	auto h_met__pt = finalDF.Histo1D({
	("met__pt_"+temp_header).c_str(),
	("met__pt "+temp_header).c_str(),
	50,0,300},"MET_pt");
	h_met__pt->GetXaxis()->SetTitle("MET p_{T} (GeV/c)");
	h_met__pt->GetYaxis()->SetTitle("Event");
	h_met__pt->SetLineStyle(kSolid);

	auto h_trans_T = finalDF.Histo1D({
	(          "tTm_"     + temp_header).c_str(),
	("Transverse T mass " + temp_header).c_str(),
	50,0,250},
	"ttop_mas","nw_ttop_mas");
	h_trans_T->GetXaxis()->SetTitle("\\text{mass GeV/}c^{2}");
	h_trans_T->GetYaxis()->SetTitle("Event");
	h_trans_T->SetLineStyle(kSolid);

	auto h_trans_w = finalDF.Histo1D({
	(          "tWm_"     + temp_header).c_str(),
	("Transverse W mass " + temp_header).c_str(),
	50,0,180},
	"tw_lep_mas","nw_tw_lep_mas");
	h_trans_w->GetXaxis()->SetTitle("\\text{mass GeV/}c^{2}");
	h_trans_w->GetYaxis()->SetTitle("Event");
	h_trans_w->SetLineStyle(kSolid);

	auto h_Winvmas = finalDF.Histo1D({
	("W_invariant_mass_" + temp_header).c_str(),
	("W invariant mass " + temp_header).c_str(),
	50,0,180},
	"lep_nu_invmass","nw_lep_nu_invmass");
	h_Winvmas->GetXaxis()->SetTitle("\\text{mass GeV/}c^{2}");
	h_Winvmas->GetYaxis()->SetTitle("Event");
	h_Winvmas->SetLineStyle(kSolid);

	auto h_ev_w = finalDF.Histo1D({
	(   "ev_w_"     +temp_header).c_str(),
	("Event weight "+temp_header).c_str(),
	50,-1,3},"sf");
	h_ev_w->GetXaxis()->SetTitle("weight");
	h_ev_w->GetYaxis()->SetTitle("Event" );
	h_ev_w->SetLineStyle(kSolid);

	auto h_z_mas = finalDF.Histo1D({
	(        "zmas_"  + temp_header).c_str(),
	("Recon. Z mass " + temp_header).c_str(),
	50,0,150},
	"z_mas","nw_z_mas");
	h_z_mas->GetXaxis()->SetTitle("\\text{mass GeV/}c^{2}");
	h_z_mas->GetYaxis()->SetTitle("Event");
	h_z_mas->SetLineStyle(kSolid);

	auto h_zw_Dph = finalDF.Histo1D({
	("Z_W_Delta_Phi_" + temp_header).c_str(),
	("Z W Delta Phi " + temp_header).c_str(),
	50,0,3.2},
	"zw_Dph","nw_zw_Dph");
	h_zw_Dph->GetXaxis()->SetTitle("Z & W #Delta#phi/rad");
	h_zw_Dph->GetYaxis()->SetTitle("Event");
	h_zw_Dph->SetLineStyle(kSolid);

	auto h_zmet_Dph = finalDF.Histo1D({
	("Z_MET_Delta_Phi_" + temp_header).c_str(),
	("Z MET Delta Phi " + temp_header).c_str(),
	50,0,3.2},
	"zmet_Dph","nw_zmet_Dph");
	h_zmet_Dph->GetXaxis()->SetTitle("Z & MET #Delta#phi/rad");
	h_zmet_Dph->GetYaxis()->SetTitle("Event");
	h_zmet_Dph->SetLineStyle(kSolid);

	auto h_z_daughters_Dph = finalDF.Histo1D({
	("Z_pair_jets_Delta_Phi_" + temp_header).c_str(),
	("Z pair jets Delta Phi " + temp_header).c_str(),
	50,0,7},
	"z_pair_Dph","nw_z_pair_Dph");
	h_z_daughters_Dph->GetXaxis()->SetTitle("Z pair jets #Delta#phi/rad");
	h_z_daughters_Dph->GetYaxis()->SetTitle("Event");
	h_z_daughters_Dph->SetLineStyle(kSolid);

	auto h_tWmVsZmass = finalDF.Histo2D({
	("tWmVsZmass_" + temp_header).c_str(),
	("tWmVsZmass " + temp_header).c_str(),
	50,0,180,50,0,150},
	"tw_lep_mas","z_mas");
	h_tWmVsZmass->GetXaxis()->SetTitle("\\text{  tWm  GeV/}c^{2}");
	h_tWmVsZmass->GetYaxis()->SetTitle("\\text{Z mass GeV/}c^{2}");
	// No SetLineStyle here

	// write histograms to a root file
	// ASSUMES temp_header is correct!
	TFile hf(("histo/"+temp_header+".root").c_str(),"RECREATE");
// MC only
		hf.WriteTObject(h_sfi                  .GetPtr());hf.Flush();sync();
		hf.WriteTObject(h_sfj                  .GetPtr());hf.Flush();sync();
//		hf.WriteTObject(h_p_ei                 .GetPtr());hf.Flush();sync();
		hf.WriteTObject(h_p_ej                 .GetPtr());hf.Flush();sync();
		hf.WriteTObject(h_p_sf_i               .GetPtr());hf.Flush();sync();
		hf.WriteTObject(h_p_sfej               .GetPtr());hf.Flush();sync();
		hf.WriteTObject(h_btag_w               .GetPtr());hf.Flush();sync();
		hf.WriteTObject(h_cmet_sEt             .GetPtr());hf.Flush();sync();
		hf.WriteTObject(h_cmet__pt             .GetPtr());hf.Flush();sync();
		hf.WriteTObject(h_cmet_dpx             .GetPtr());hf.Flush();sync();
		hf.WriteTObject(h_cmet_dpy             .GetPtr());hf.Flush();sync();
		hf.WriteTObject(h_is_btag_numer_PtVsEta.GetPtr());hf.Flush();sync();
		hf.WriteTObject(h_no_btag_numer_PtVsEta.GetPtr());hf.Flush();sync();
		hf.WriteTObject(h_is_btag_denom_PtVsEta.GetPtr());hf.Flush();sync();
		hf.WriteTObject(h_no_btag_denom_PtVsEta.GetPtr());hf.Flush();sync();
		hf.WriteTObject(is_btag_ratio)                   ;hf.Flush();sync();
		hf.WriteTObject(no_btag_ratio)                   ;hf.Flush();sync();
// end MC only;
	hf.WriteTObject(h_met_sEt        .GetPtr());hf.Flush();sync();
	hf.WriteTObject(h_met__pt        .GetPtr());hf.Flush();sync();
	hf.WriteTObject(h_trans_T        .GetPtr());hf.Flush();sync();
	hf.WriteTObject(h_trans_w        .GetPtr());hf.Flush();sync();
	hf.WriteTObject(h_Winvmas        .GetPtr());hf.Flush();sync();
	hf.WriteTObject(h_ev_w           .GetPtr());hf.Flush();sync();
	hf.WriteTObject(h_z_mas          .GetPtr());hf.Flush();sync();
	hf.WriteTObject(  h_zw_Dph       .GetPtr());hf.Flush();sync();
	hf.WriteTObject(h_zmet_Dph       .GetPtr());hf.Flush();sync();
	hf.WriteTObject(h_z_daughters_Dph.GetPtr());hf.Flush();sync();
	hf.WriteTObject(h_tWmVsZmass     .GetPtr());hf.Flush();sync();
	// the following two for loops stack correctly
	for(std::string particle:{"fin_jets","lep","bjet"})
	for(PtEtaPhiM k:PtEtaPhiMall){
//		if( e == k ) continue;
		std::string  kstring  = "_";
		std::string  xAxisStr;
		double xmin,xmax;
		switch(k){
			case pt :{kstring += "_pt";xmin =  0;xmax = 200;
			          xAxisStr = "p_{T} (GeV/c)"           ;break;}
			case eta:{kstring += "eta";xmin = -3;xmax =  3 ;
			          xAxisStr = "PseudoRapidity  #eta"    ;break;}
			case phi:{kstring += "phi";xmin = -7;xmax =  7 ;
			          xAxisStr = "Azimuthal angle #phi/rad";break;}
			case  m :{kstring += "mas";xmin =  0;xmax = 200;
			          xAxisStr = "\\text{mass GeV/}c^{2}"  ;break;}
//			case  e :
			default :throw std::invalid_argument(
				"Unimplemented component (histo)");
		}
		temp_footer = particle    + kstring;
		temp_opener = temp_header + "_" + temp_footer;
		auto
		h = finalDF.Histo1D(
		{temp_opener.c_str()
		,temp_opener.c_str()
		,50,xmin,xmax}
		,       temp_footer .c_str()
		,("nw_"+temp_footer).c_str()
		);
		h->SetLineStyle(kSolid);
		h->GetXaxis()->SetTitle(xAxisStr.c_str());
		h->GetYaxis()->SetTitle("Event");
		hf.WriteTObject(h.GetPtr());hf.Flush();sync();
	}} else {
	auto expt_bjets
	   = init_selection
	.Filter(runLBfilter(runLBdict),{"run","luminosityBlock"},
	        "LuminosityBlock filter")
	.Define("roccorSF", roccorSF(rc,ch,MC)// used in allReconstruction
	                  ,{"lepB_pt","lep_eta",
	                    "lep_phi","lep___q", // last 4 unused
	                    "lep_phi","lep___q"})// last 2 repeat is fine
	.Define("cmet__pt"    ,"static_cast<double>(MET_pt )")// no need sumEt
	.Define("cmet_phi"    ,"static_cast<double>(MET_phi)")
	.Define("fin_jets__pt","static_cast<ROOT::RVec<double>>(Jet_pt  [tight_jets])")
	.Define("fin_jets_eta","static_cast<ROOT::RVec<double>>(Jet_eta [tight_jets])")
	.Define("fin_jets_phi","static_cast<ROOT::RVec<double>>(Jet_phi [tight_jets])")
	.Define("fin_jets_mas","static_cast<ROOT::RVec<double>>(Jet_mass[tight_jets])")
	.Define("btagP"            ,btagP  ,{"fin_jets_eta"})// suPer vs suBset
	.Define("btagB"            ,btagB  ,{"btagP","tJ_btagCSVv2"})
	.Filter(jetCutter(BJETS_MIN,BJETS_MAX),{"btagB"},"b jet cut")
	// TODO: Always check that the previous 3 lines are copies of earlier
	;
	auto reco = allReconstruction(
	     expt_bjets )
	;
	auto not_btag_eff
	   = reco
	.Define("btag_w" , [](){return 1.;})
	.Define("lepSF"  , lepEffGiver(ch,MC// not in reco because files
	                 , recoLowEt,reco_pass,tight_94x
	                 , id_N,id_Y,id_A,id_T
	                 , isoN,isoY,isoA,isoT
	                ),{"lep__pt","lep_eta"})
	. Alias("mostSF" , "lepSF")
	;
	auto finalDF = finalScaling(ds,PuWd,PuUd,PuDd,// unused but send pile
	     not_btag_eff )
	;
	// Copied from earlier, delete MC-only!
	// Assuming temp_header and footer and all are set per (hist titles)!
	auto h_met_sEt = finalDF.Histo1D({
	("met_sEt_"+temp_header).c_str(),
	("met_sEt "+temp_header).c_str(),
	100,0,600},"MET_sumEt");
	h_met_sEt->GetXaxis()->SetTitle("MET Sum Et (GeV)");
	h_met_sEt->GetYaxis()->SetTitle("Event");
	h_met_sEt->SetLineStyle(kSolid);

	auto h_met__pt = finalDF.Histo1D({
	("met__pt_"+temp_header).c_str(),
	("met__pt "+temp_header).c_str(),
	50,0,300},"MET_pt");
	h_met__pt->GetXaxis()->SetTitle("MET p_{T} (GeV/c)");
	h_met__pt->GetYaxis()->SetTitle("Event");
	h_met__pt->SetLineStyle(kSolid);

	auto h_trans_T = finalDF.Histo1D({
	(          "tTm_"     + temp_header).c_str(),
	("Transverse T mass " + temp_header).c_str(),
	50,0,250},
	"ttop_mas","nw_ttop_mas");
	h_trans_T->GetXaxis()->SetTitle("\\text{mass GeV/}c^{2}");
	h_trans_T->GetYaxis()->SetTitle("Event");
	h_trans_T->SetLineStyle(kSolid);

	auto h_trans_w = finalDF.Histo1D({
	(          "tWm_"     + temp_header).c_str(),
	("Transverse W mass " + temp_header).c_str(),
	50,0,180},
	"tw_lep_mas","nw_tw_lep_mas");
	h_trans_w->GetXaxis()->SetTitle("\\text{mass GeV/}c^{2}");
	h_trans_w->GetYaxis()->SetTitle("Event");
	h_trans_w->SetLineStyle(kSolid);

	auto h_Winvmas = finalDF.Histo1D({
	("W_invariant_mass_" + temp_header).c_str(),
	("W invariant mass " + temp_header).c_str(),
	50,0,180},
	"lep_nu_invmass","nw_lep_nu_invmass");
	h_Winvmas->GetXaxis()->SetTitle("\\text{mass GeV/}c^{2}");
	h_Winvmas->GetYaxis()->SetTitle("Event");
	h_Winvmas->SetLineStyle(kSolid);

	auto h_ev_w = finalDF.Histo1D({
	(   "ev_w_"     +temp_header).c_str(),
	("Event weight "+temp_header).c_str(),
	50,-1,3},"sf");
	h_ev_w->GetXaxis()->SetTitle("weight");
	h_ev_w->GetYaxis()->SetTitle("Event" );
	h_ev_w->SetLineStyle(kSolid);

	auto h_z_mas = finalDF.Histo1D({
	(        "zmas_"  + temp_header).c_str(),
	("Recon. Z mass " + temp_header).c_str(),
	50,0,150},
	"z_mas","nw_z_mas");
	h_z_mas->GetXaxis()->SetTitle("\\text{mass GeV/}c^{2}");
	h_z_mas->GetYaxis()->SetTitle("Event");
	h_z_mas->SetLineStyle(kSolid);

	auto h_zw_Dph = finalDF.Histo1D({
	("Z_W_Delta_Phi_" + temp_header).c_str(),
	("Z W Delta Phi " + temp_header).c_str(),
	50,0,3.2},
	"zw_Dph","nw_zw_Dph");
	h_zw_Dph->GetXaxis()->SetTitle("Z & W #Delta#phi/rad");
	h_zw_Dph->GetYaxis()->SetTitle("Event");
	h_zw_Dph->SetLineStyle(kSolid);

	auto h_zmet_Dph = finalDF.Histo1D({
	("Z_MET_Delta_Phi_" + temp_header).c_str(),
	("Z MET Delta Phi " + temp_header).c_str(),
	50,0,3.2},
	"zmet_Dph","nw_zmet_Dph");
	h_zmet_Dph->GetXaxis()->SetTitle("Z & MET #Delta#phi/rad");
	h_zmet_Dph->GetYaxis()->SetTitle("Event");
	h_zmet_Dph->SetLineStyle(kSolid);

	auto h_z_daughters_Dph = finalDF.Histo1D({
	("Z_pair_jets_Delta_Phi_" + temp_header).c_str(),
	("Z pair jets Delta Phi " + temp_header).c_str(),
	50,0,7},
	"z_pair_Dph","nw_z_pair_Dph");
	h_z_daughters_Dph->GetXaxis()->SetTitle("Z pair jets #Delta#phi/rad");
	h_z_daughters_Dph->GetYaxis()->SetTitle("Event");
	h_z_daughters_Dph->SetLineStyle(kSolid);

	auto h_tWmVsZmass = finalDF.Histo2D({
	("tWmVsZmass_" + temp_header).c_str(),
	("tWmVsZmass " + temp_header).c_str(),
	50,0,180,50,0,150},
	"tw_lep_mas","z_mas");
	h_tWmVsZmass->GetXaxis()->SetTitle("\\text{  tWm  GeV/}c^{2}");
	h_tWmVsZmass->GetYaxis()->SetTitle("\\text{Z mass GeV/}c^{2}");
	// No SetLineStyle here

	// write histograms to a root file
	// ASSUMES temp_header is correct!
	TFile hf(("histo/"+temp_header+".root").c_str(),"RECREATE");
	hf.WriteTObject(h_met_sEt        .GetPtr());hf.Flush();sync();
	hf.WriteTObject(h_met__pt        .GetPtr());hf.Flush();sync();
	hf.WriteTObject(h_trans_T        .GetPtr());hf.Flush();sync();
	hf.WriteTObject(h_trans_w        .GetPtr());hf.Flush();sync();
	hf.WriteTObject(h_Winvmas        .GetPtr());hf.Flush();sync();
	hf.WriteTObject(h_ev_w           .GetPtr());hf.Flush();sync();
	hf.WriteTObject(h_z_mas          .GetPtr());hf.Flush();sync();
	hf.WriteTObject(  h_zw_Dph       .GetPtr());hf.Flush();sync();
	hf.WriteTObject(h_zmet_Dph       .GetPtr());hf.Flush();sync();
	hf.WriteTObject(h_z_daughters_Dph.GetPtr());hf.Flush();sync();
	hf.WriteTObject(h_tWmVsZmass     .GetPtr());hf.Flush();sync();
	// the following two for loops stack correctly
	for(std::string particle:{"fin_jets","lep","bjet"})
	for(PtEtaPhiM k:PtEtaPhiMall){
//		if( e == k ) continue;
		std::string  kstring  = "_";
		std::string  xAxisStr;
		double xmin,xmax;
		switch(k){
			case pt :{kstring += "_pt";xmin =  0;xmax = 200;
			          xAxisStr = "p_{T} (GeV/c)"           ;break;}
			case eta:{kstring += "eta";xmin = -3;xmax =  3 ;
			          xAxisStr = "PseudoRapidity  #eta"    ;break;}
			case phi:{kstring += "phi";xmin = -7;xmax =  7 ;
			          xAxisStr = "Azimuthal angle #phi/rad";break;}
			case  m :{kstring += "mas";xmin =  0;xmax = 200;
			          xAxisStr = "\\text{mass GeV/}c^{2}"  ;break;}
//			case  e :
			default :throw std::invalid_argument(
				"Unimplemented component (histo)");
		}
		temp_footer = particle    + kstring;
		temp_opener = temp_header + "_" + temp_footer;
		auto
		h = finalDF.Histo1D(
		{temp_opener.c_str()
		,temp_opener.c_str()
		,50,xmin,xmax}
		,       temp_footer .c_str()
		,("nw_"+temp_footer).c_str()
		);
		h->SetLineStyle(kSolid);
		h->GetXaxis()->SetTitle(xAxisStr.c_str());
		h->GetYaxis()->SetTitle("Event");
		hf.WriteTObject(h.GetPtr());hf.Flush();sync();
	}}
	std::cout<<"calchisto successfully completed"<<std::endl;
}
