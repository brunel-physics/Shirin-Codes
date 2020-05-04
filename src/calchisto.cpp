// clang++ -std=c++17 src/calchisto.cpp src/eval_complex.cpp -o calchisto.o ` root-config --libs` 
#include <ROOT/RDataFrame.hxx>
#include <TLorentzVector.h>
#include <TChain.h>
#include <TRandom3.h>

#include <boost/algorithm/string.hpp>

#include "csv.h"
#include "calchisto.hpp"
#include "eval_complex.hpp"

using doubles = ROOT::VecOps::RVec<double>;
using  floats = ROOT::VecOps::RVec<float>;
using    ints = ROOT::VecOps::RVec<int>;
using   bools = ROOT::VecOps::RVec<bool>;
using   chars = ROOT::VecOps::RVec<UChar_t>; // aka 1 byte ints
using strings = ROOT::VecOps::RVec<std::string>;

namespace{
//constexpr double EL_MAX_NUM  = 1;
constexpr double EL__PT_MIN  = 45;//{15}//min 12, AP 45, 
//constexpr float  EL_LPT_MIN  = 35.f;// Leading
constexpr double EL_ETA_MAX  = 2.5;
constexpr int    EL_LOOSE_ID = 1;
constexpr int    EL_TIGHT_ID = 4;

constexpr double ENDCAP_ETA_MIN = 1.566;
constexpr double BARREL_ETA_MAX = 1.4442;

//constexpr double MU_MAX_NUM   = 1;
constexpr double MU__PT_MIN   = 40;//min 33, AP 40, 
//constexpr float  MU_LPT_MIN   = 26.f; // Leading
constexpr double MU_ETA_MAX   = 2.4;
//constexpr float  MU_LOOSE_ISO = .15f;
//constexpr float  MU_TIGHT_ISO = .25f;

//constexpr double MET__PT_MIN = 40;
constexpr double MET_EL_PT   = 80;
//constexpr double MET_MU_PT   = 40;

constexpr double   Z_MASS     = 91.1876;
//constexpr double   Z_MASS_CUT = 20.;
constexpr double   W_MASS     = 80.385;
constexpr double   W_MASS_CUT = 20.;
constexpr float  TOP_MASS     = 172.5f;
//constexpr float  TOP_MASS_CUT = 20.f;

constexpr float JET_ETA_MAX = 4.7f;
constexpr float JET__PT_MIN = 30.f;
constexpr float     JET_ISO = .4f;
constexpr unsigned JETS_MIN = 4;
constexpr unsigned JETS_MAX = 6;

constexpr float  BJET_ETA_MAX = 2.4f;
constexpr float BTAG_DISC_MIN =  .8838f;
constexpr unsigned  BJETS_MIN = 1;
constexpr unsigned  BJETS_MAX = 3;

constexpr float RconeBy2 = .2f;

//constexpr float DELTA___R_ZL   = 1.6f;
//constexpr float DELTA_PHI_ZW   = 2;
//constexpr float DELTA_PHI_ZMET = 2;

constexpr float    TZQ_W =  .0128f;
constexpr float WWLNQQ_W = 2.1740f;
constexpr float WZLNQQ_W =  .2335f;
constexpr float  TTZQQ_W =  .0237f;
constexpr float ZZLLQQ_W =  .0485f;

// This Pi is more accurate than binary256; good for eternity
template <typename T> constexpr T  PI = T(3.14159265358979323846264338327950288419716939937510582097494459230781640628620899);
template <typename T> constexpr T TPI = PI<T> * 2;

enum PtEtaPhiM {pt,eta,phi,m};

template <typename T>
[[gnu::const]] T select(const T& a, const ints& mask){return a[mask];}
auto   el_sel(const    int  target_id){
	return [=](const  bools& isPFs,
	           const floats& pts,
	           const floats& etas,
	           const   ints& ids){
		const auto  abs_etas =  abs(etas);
		return (isPFs && pts >  EL__PT_MIN
		       && ((abs_etas <  EL_ETA_MAX
		       &&   abs_etas >  ENDCAP_ETA_MIN)
		       ||  (abs_etas <  BARREL_ETA_MAX))
		       &&        ids >= target_id);
	};
}
auto   mu_sel(const float   target_iso){
	return [=](const  bools& isPFs,
	           const floats& pts,
	           const floats& etas,
	           const  bools& ids,
	           const floats& isos){
		const auto abs_etas = abs(etas);
		return (   ids      && isPFs
		       &&  pts      >  MU__PT_MIN
		       &&  abs_etas <  MU_ETA_MAX 
		       &&  isos     <= target_iso);
	};
}
/*
auto is_good_el 
	=[](const int     target_id,
		 const bools&  isPFs,
		 const floats& pts,
		 const floats& etas,
		 const ints&   ids){
	const auto  abs_etas = abs(etas);
	return (isPFs && pts > EL__PT_MIN
		    && ((abs_etas < EL_ETA_MAX 
		    &&   abs_etas > ENDCAP_ETA_MIN) 
		    ||  (abs_etas < BARREL_ETA_MAX))
		    && ids >= target_id);
};
auto is_tight_el
	=[]
	 (const bools&  isPFs,
	  const floats& pts,
	  const floats& etas,
	  const ints&   ids){return is_good_el(EL_TIGHT_ID,isPFs,pts,etas,ids);};
auto is_loose_el
	=[]
	 (const bools&  isPFs,
	  const floats& pts,
	  const floats& etas,
	  const ints&   ids){return is_good_el(EL_LOOSE_ID,isPFs,pts,etas,ids);};
auto is_good_mu
	=[](const float   target_iso,
		 const bools&  isPFs,
		 const floats& pts,
		 const floats& etas,
		 const bools&  ids,
		 const floats& isos){
	const auto abs_etas = abs(etas);
	return (isPFs && ids
		    &&  pts      >  MU__PT_MIN
		    &&  abs_etas <  MU_ETA_MAX 
		    &&  isos     <= target_iso);
};
auto is_tight_mu
	=[]
	 (const bools&  isPFs,
	  const floats& pts,
	  const floats& etas,
	  const bools&  ids,
	  const floats& isos)
	{return is_good_mu(MU_TIGHT_ISO,isPFs,pts,etas,ids,isos);};
auto is_loose_mu
	=[]
	 (const bools&  isPFs, 
	  const floats& pts, 
	  const floats& etas, 
	  const bools&  ids, 
	  const floats& isos)
	{return is_good_mu(MU_LOOSE_ISO,isPFs,pts,etas,ids,isos);};
auto el_cut
	=[](const floats& tight_el_pts,
	    const floats& loose_el_pts){
	const bool el_cut
		=tight_el_pts.size()==1&&tight_el_pts.size()==loose_el_pts.size();
	//bool lead_pt_cut{false};
	//lead_pt_cut = tight_el_pts.empty() ? false : *std::max_element(tight_el_pts.begin(), tight_el_pts.end()) > EL__PT_MIN;
	return el_cut;
	//return lead_pt_cut && ele_cut;
};
auto mu_cut
	=[](const floats& tight_mu_pts,
	    const floats& loose_mu_pts){
	const bool mu_cut
		=tight_mu_pts.size()==1&&tight_mu_pts.size()==loose_mu_pts.size();
	//bool lead_pt_cut{false};
	//lead_pt_cut = tight_mu_pts.empty() ? false : *std::max_element(tight_mu_pts.begin(), tight_mu_pts.end()) > MU__PT_MIN;
	//return lead_pt_cut && mu_cut;
	return mu_cut;
};
*/
auto lep_cut = [](const floats& tight_pts,const floats& loose_pts){
	// Do we need to make this a function of EL__PT_MIN? 
	const bool result 
		= tight_pts.size()==1&&tight_pts.size()==loose_pts.size();
	//const bool lead_pt_cut 
	/*	= tight_pts.empty() ? false : 
		*std::max_element(tight_pts.begin(),
		                  tight_pts  .end()) > MU__PT_MIN;
	 */
	//return lead_pt_cut && result;
	return result;
};
/*
auto get_w_el_string
	=[](const std::string& s){return "Electron_" + s + "[tight_els]";};
auto get_w_mu_string
	=[](const std::string& s){return     "Muon_" + s + "[tight_mus]";};
*/
auto get_w_lep_string(const channels ch,const PtEtaPhiM ss){
	std::string result;
	switch (ss){
		case  pt:{result="pt"  ;break;}
		case eta:{result="eta" ;break;}
		case phi:{result="phi" ;break;}
		case   m:{result="mass";break;} // Expect to never use this? 
		default :{throw std::invalid_argument("Unimplemented PtEtaPhiM");}
	}
	switch (ch){
		case elnu:{result="Electron_"+ result +"[tight_els]";break;}
		case munu:{result=    "Muon_"+ result +"[tight_mus]";break;}
		default  :{throw std::invalid_argument("Unimplemented ch");}
	}
	return result;
}
auto lep_met_selector(const float lower_bound){
	return [=](const float&  met_lep_pt_sel)->bool
	   {return lower_bound < met_lep_pt_sel;};
}
/*
auto el_met_selector
	=[](const float&    MET_el_pt_sel)->bool
	{return MET_EL_PT < MET_el_pt_sel;};
auto mu_met_selector
	=[](const float&    MET_mu_pt_sel)->bool
	{return MET_MU_PT < MET_mu_pt_sel;};
*/
template <typename T>
[[gnu::const]] auto delta_phi(const T diff_phi){
	// This function just reduces input from [-2pi,+2pi] to [-pi,+pi]
	// Domain correctness is the user's responsibility
	     if(diff_phi >  PI<T>)  return diff_phi - TPI<T>;
	else if(diff_phi < -PI<T>)  return diff_phi + TPI<T>;
	else                        return diff_phi;
}
[[gnu::const]] auto deltaR(
	const float eta1,const float phi1,const float eta2,const float phi2)
{return std::hypot(eta1-eta2,delta_phi(phi1-phi2));}
auto transverse_w_mass
	=[](const floats& lep__pt,
	    const floats& lep_phi,
	    const float & met_pt , // not typo
	    const float & met_phi){
	//float w_reco_mass{std::numeric_limits<float>::infinity()};
	//size_t  l_index_1{std::numeric_limits<size_t>::max()};
	floats   w_mass_vec(lep__pt.size());
	for(size_t i=0; i < lep__pt.size() ;++i){
		w_mass_vec[i] = static_cast<float>(
		                std::sqrt(2*lep__pt[i]*met_pt
		          *(1-cos(delta_phi(lep_phi[i]-met_phi)))));
		/*if(std::abs(W_MASS-reco_mass)<std::abs(W_MASS-w_reco_mass)){
			w_reco_mass = reco_mass;
			l_index_1 = i;
		}*/
	}
	return w_mass_vec;
};
/*auto w_mass_cut	= [](const floats& w_mass){
	return std::abs(w_mass-W_MASS)<W_MASS_CUT;
};*/
/*auto z_mass_cut = [](const float& z_mass){
	return std::abs(z_mass-Z_MASS)<Z_MASS_CUT;
};*/
// TODO: top_mass_cut belongs here
template<typename T,typename U>
[[gnu::const]] bool all_equal(const T& t, const U& u){return t == u;}
template<typename T,typename U,typename... Types>
[[gnu::const]] bool all_equal(const T& t, const U& u, Types const&... args)
	{return t == u && all_equal(u, args...);}
auto lep_nu_invmass
	=[](const floats& lep_pt    ,
	    const floats& lep_eta   ,
	    const floats& lep_phi   ,
	    const floats& lep_mass  ,
	    const float & cal_metphi,
	    const float & cal_metpt , 
	    const float & cal_metEt ){ // TODO: cal_metEt unused
	// this function computes the invariant mass of charged lepton 
	// and neutrino system, in order to calculate the W mass later on.
	auto lep=TLorentzVector{};
	auto neu=TLorentzVector{};
	if(!all_equal(lep_pt.size(),lep_eta.size(),lep_phi.size(),lep_mass.size())) 
		throw std::logic_error("Collections must be the same size");
	else if(lep_pt.size() != 1)
		throw std::logic_error("Should always only have one lepton");
	float lepnu_mass = 0.f;
	for(size_t i=0; i < lep_pt.size() ;++i){
		lep.SetPtEtaPhiM(lep_pt[i],lep_eta[i],lep_phi[i],lep_mass[i]);
		neu.SetPtEtaPhiM(cal_metpt,lep_eta[i],cal_metphi,0.);
		lepnu_mass=static_cast<float>((lep+neu).M());
	}
	return lepnu_mass;
};
auto jet_lep_min_deltaR
	=[](const floats& jet_etas,
	    const floats& jet_phis,
	    const floats& lep_etas,
	    const floats& lep_phis) {
	floats min_dRs;
	std::transform(
		jet_etas.begin(),
		jet_etas  .end(),
		jet_phis.begin(),
		std::back_inserter(min_dRs),
		[&](float jet_eta, float jet_phi){
			return deltaR(jet_eta,jet_phi,lep_etas.at(0),lep_phis.at(0));});
	return min_dRs;
};
auto tight_jet_id
	=[](const floats& jet_lep_min_dRs, 
	    const floats& pts, 
	    const floats& etas, 
	    const   ints& ids){
	return pts>JET__PT_MIN&&etas<JET_ETA_MAX&&jet_lep_min_dRs>JET_ISO&&ids>=2;
};
auto jet_deltaphi = [](const floats& phis){
	floats deltaphis;
	// The following two for loops stack correctly
	for(size_t   i=0; i < phis.size()-1 ;++i) 
	for(size_t j=i+1; j < phis.size()   ;++j)
		deltaphis.push_back(  static_cast<float>(
			std::abs(delta_phi(phis[i]-phis[j]))));
	return deltaphis;
};
auto  jet_cut = [](const ints& tight_jets){
	auto njet = std::count_if(
		tight_jets.begin(),
		tight_jets  .end(),
		[](int i){return i;});
	return (njet >= JETS_MIN) && (njet <= JETS_MAX);
};
auto bjet_cut = [](const ints& bjets){
	const auto nbjet = std::count_if(
		bjets.begin(),
		bjets  .end(),
		[](int i){return i;});
	return nbjet >= BJETS_MIN && nbjet <= BJETS_MAX;
};
auto is_bjet_id
	=[](const ints& tight_jets,const floats& etas,const floats& btags)
	{return tight_jets && (    btags > BTAG_DISC_MIN) 
	                   && (abs(etas) < BJET_ETA_MAX );
};
auto no_bjet_id    = [](const ints& tight_jets,const floats& etas)
	{return tight_jets && (abs(etas) < BJET_ETA_MAX);};
auto is_bjet_numer = [](const ints& id,const ints& is_bjet)
	{return abs(id)==5 && is_bjet;};
auto is_bjet_denom = [](const ints& id,const ints& no_bjet)
// using non_bjet_id particles not matching btag criteria
	{return abs(id)==5 && no_bjet;};
auto no_bjet_numer = [](const ints& id,const ints& is_bjet){
// using bjets which has satisfied btag conditions
	const auto aid = abs(id);
	return ((aid > 0 && aid <= 4) || aid == 21 || aid != 5) && is_bjet;
};
auto no_bjet_denom = [](const ints& id,const ints& no_bjet){
// using bjets which has satisfied non btag condition
	const auto aid = abs(id);
	return ((aid > 0 && aid <= 4) || aid == 21 || aid != 5) && no_bjet;
};
auto is_btag_CSVv2
	=[](const floats& btag,
	    const floats& pt,
	    const floats& eta){
	strings formulae(pt.size(),"0"); // vector of "0"
	floats   results(pt.size());
	if(!all_equal(pt.size(),eta.size(),btag.size()))
		throw std::logic_error("Collections must be the same size");
	else if(pt.empty())
		throw std::logic_error("Collections must not be empty");
	std::string  measure_type,sys_type,rawFormula;
	float CSVv2  ,jet_flav;
	float  pt_min,  pt_max;
	float eta_min, eta_max;
	float CSV_min, CSV_max;
	io::CSVReader<11> thisCSVfile("CSVv2_94XSF_V2_B_F.csv");
	// The following nests too much, so we do not indent
	// Each blank line means nesting deeper
	while(thisCSVfile.read_row(CSVv2,measure_type,sys_type,jet_flav,
		eta_min,eta_max,pt_min,pt_max,CSV_min,CSV_max,rawFormula)){

	if(measure_type == "comb" && CSVv2 >= BTAG_DISC_MIN
	&& sys_type  == "central" && jet_flav == 0.f){

	for(size_t i=0; i < pt.size() ;++i){

	std::string tempFormula = rawFormula;
	if( eta[i] > eta_min &&  eta[i] < eta_max
	&&   pt[i] > pt_min  &&   pt[i] <  pt_max
	&& btag[i] > CSV_min && btag[i] < CSV_max){

	if(formulae[i] == "0"){ // only 1st found wins

	if(tempFormula.find("x") != std::string::npos){

	boost::replace_all(tempFormula,"x",std::to_string(pt[i]));
	formulae[i] = tempFormula;
	}}}}}
	} // No need to close file after this while loop.
	// resume indentation
	for(size_t j=0; j < formulae.size() ;++j){
		Eval ev;
		results[j] = static_cast<float>(
		ev.eval(      const_cast<char*>(
		        formulae[j].c_str())).real());
	}
	return results;
};
auto no_btag_CSVv2
	=[](const floats& btag,
	    const floats& pt,
	    const floats& eta){
	strings formulae(pt.size(),"0"); // vector of "0" 
	floats   results(pt.size());
	if(!all_equal(   pt.size(),eta.size(),btag.size()))
		throw std::logic_error("Collections must be the same size");
	else if(pt.empty())
		throw std::logic_error("Collections must not be empty");
	std::string  measure_type,sys_type,rawFormula;
	float CSVv2  ,jet_flav;
	float  pt_min,  pt_max;
	float eta_min, eta_max;
	float CSV_min, CSV_max;
	io::CSVReader<11> thisCSVfile("CSVv2_94XSF_V2_B_F.csv");
	// The following nests too much, so we do not indent
	// Each blank line means nesting deeper
	while(thisCSVfile.read_row(CSVv2,measure_type,sys_type,jet_flav,
		eta_min,eta_max,pt_min,pt_max,CSV_min,CSV_max,rawFormula)){

	if(measure_type == "comb" // No CSVv2 here
	&& sys_type ==  "central" && jet_flav == 0.f){

	for(size_t i=0; i < pt.size() ;++i){

	std::string tempFormula = rawFormula;
	if( eta[i] > eta_min &&  eta[i] < eta_max 
	&&   pt[i] > pt_min  &&   pt[i] <  pt_max){
	// No btag checks

	if(formulae[i] == "0"){ // only 1st found wins 

	if(tempFormula.find("x") != std::string::npos){

	boost::replace_all(tempFormula,"x",std::to_string(pt[i]));
	formulae[i] = tempFormula; 
	}}}}}
	} // No need to close file after this while loop.
	// resume indentation
	for(size_t j=0; j < formulae.size() ;++j){
		Eval ev;
		results[j] = static_cast<float>(
		ev.eval(      const_cast<char*>(
		        formulae[j].c_str())).real());
	}
	return results;
};
auto jet_smear_pt_resol
	=[](const floats& pt,const floats& eta,const float& rho){
	float min_eta,max_eta,min_rho,max_rho;
	int   z; // unwanted
	float min_pt,max_pt;
	float a,b,c,d;
	floats  resol(pt.size(),0.f); // vec of 0.
	if(!all_equal(pt.size(),eta.size()))
		throw std::logic_error("Collections must be the same size");
	else if(pt.empty())
		throw std::logic_error("Collections must not be empty");
	io::CSVReader<11> thisCSVfile("Fall17_V3_MC_PtResolution_AK4PFchs.txt");
	while(thisCSVfile.read_row(
		min_eta,max_eta,min_rho,max_rho,z,min_pt,max_pt,a,b,c,d)){
		// The following if for if stacks correctly
		if(rho > min_rho && rho < max_rho)
		for(size_t i=0; i < pt.size() ;++i)
		if(eta[i] > min_eta && eta[i] < max_eta
		&&  pt[i] > min_pt  &&  pt[i] < max_pt){
			resol[i] += std::sqrt(  a * std::abs(a)/(pt[i]*pt[i])
			                    + b*b * std::pow(pt[i],d) + c*c);
		}
	} // No need to close file after this while loop
	return resol;
};
auto jet_smear_Sjer=[](const floats& eta){
	float min_eta,max_eta;
	int   z; // unwanted
	float central_SF,SF_dn,SF_up;
	floats     Sjer(eta.size(),0.f); // vec of 0.
	io::CSVReader<6> thisCSVfile("Fall17_V3_MC_SF_AK4PF.txt");
	while(thisCSVfile.read_row(min_eta,max_eta,z,central_SF,SF_dn,SF_up)){
		for(size_t i=0; i  <  eta.size() ;++i)
		   if(      eta[i] >  min_eta && eta[i] < max_eta)
		           Sjer[i] += central_SF;
	}
	return Sjer;
};
[[gnu::const]] auto gapRamp(const float Sjer,const float x){
	if(Sjer > x) return Sjer;
	else return 0.f;
}
auto delta_R_jet_smear
	=[](const floats& pt,
	    const floats& gen_pt,
	    const floats& resol,
	    const floats& Sjer,
	    const floats& deltaR){
	if(!all_equal(pt.size(),resol.size(),Sjer.size(),deltaR.size()))
		throw std::logic_error("Collections must be the same size");
	else if(pt.empty())
		throw std::logic_error("Collections must not be empty");
	// TODO: Finalise size equality guarantees
	long size_diff = static_cast<long>(gen_pt.size()) - static_cast<long>(pt.size());
	if(  size_diff < 0 ) throw std::logic_error("Insufficient gen_pt");
	const size_t  size = pt.size();
//	         = size_diff > 0 ?
//	           gen_pt.size() : pt.size();
//	size_diff = std::abs(size_diff);
	floats cjers(    size,0.f); //  correction factor // TODO: wrong?
	for(size_t i=0;i<size;++i){
		if(deltaR.at(i) < RconeBy2
		&& std::abs(pt.at(i)-gen_pt.at(i)) < 3*resol.at(i)*pt.at(i)){
			cjers[i] += (1+(1+Sjer.at(i))
			     * ((pt.at(i)-gen_pt.at(i))/pt.at(i)));
		}
		else{ // needs TRandom3 library
			float Normdist = static_cast<float>(gRandom->Gaus(0,Sjer.at(i)));
			float  max_val = Sjer.at(i) * Sjer.at(i) - 1;
			cjers[i] += (1+Normdist*std::sqrt(gapRamp(max_val,0)));
		}
	}
	return cjers;
};
auto cjer = [](const floats& jet,const floats& cjer){
	floats  weighted(jet.size());
	if(!all_equal(   jet.size(),cjer.size()))
		throw std::logic_error("Collections must be the same size");
	else if(jet.empty())
		throw std::logic_error("Collections must not be empty");
	for(size_t i=0; i < jet.size() ;++i) weighted[i] = jet[i]*cjer[i];
	return weighted;
};
auto find_lead_mask = [](const ints& mask,const floats& vals){
	const auto masked_vals = mask * vals;
	const auto max_idx = static_cast<size_t>(std::distance(
		masked_vals.begin(),
		max_element(masked_vals.begin(),
		            masked_vals  .end())));
	ints lead_mask(masked_vals .size(),0);
	lead_mask[max_idx] = 1;
	return lead_mask;
};
auto find_z_pair
	=[](const floats& pts,
	    const floats& etas,
	    const floats& phis,
	    const floats& ms,
	    const   ints& tight_jets,
	    const   ints& lead_bjet){
	// This function finds the pair nearest to z mass
	double z_reco_mass = std::numeric_limits<double>::infinity();
	size_t jet_index_1 = std::numeric_limits<size_t>::max();
	size_t jet_index_2 = std::numeric_limits<size_t>::max();
	const size_t njets = tight_jets.size(); 
	if(!all_equal(pts.size(),etas.size(),phis.size(),ms.size(),
		njets,lead_bjet.size()))
		throw std::logic_error("Collections must be the same size");
	else if(njets==0)
		throw std::logic_error("Collections must not be empty");
	// The next two for loops stack correctly with the if 
	for(size_t   i=0; i < njets-1 ;++i){
	for(size_t j=i+1; j < njets   ;++j){
	if(tight_jets[i] == 0 || tight_jets[j] == 0 || lead_bjet[i] == 1 || lead_bjet[j] == 1){
		auto jet1 = TLorentzVector{};
		auto jet2 = TLorentzVector{};
		jet1.SetPtEtaPhiM(pts[i],etas[i],phis[i],ms[i]);
		jet2.SetPtEtaPhiM(pts[j],etas[j],phis[j],ms[j]);
		if (const double reco_mass=(jet1+jet2).M();
		std::abs(Z_MASS-reco_mass) < std::abs(Z_MASS-z_reco_mass)){ 
			z_reco_mass = reco_mass; // found nearer pair to z mass
			jet_index_1 = i; 
			jet_index_2 = j;
		}
	}}}
	ints z_pair(njets,0); // all zeroed 
	z_pair[jet_index_1] = 1;
	z_pair[jet_index_2] = 1;
	return z_pair;
};
[[gnu::const]] auto inv_mass(
	const floats& pts,const floats& etas,const floats& phis,const floats& ms){
	if(!all_equal(pts.size(),etas.size(),phis.size(),ms.size())) 
		throw std::logic_error("Collections must be the same size");
	else if(pts.empty())
		throw std::logic_error("Collections must not be empty");
	TLorentzVector vec;
	for(size_t i=0; i < pts.size() ;++i){
		TLorentzVector p;
		p.SetPtEtaPhiM(pts[i],etas[i],phis[i],ms[i]);
		vec += p;
	}
	return static_cast<float>(vec.M());
}
auto zw_deltaphi=[](const floats& phis1,const floats& phis2){
	const size_t   p1s = phis1.size();
	const size_t   p2s = phis2.size();
	floats results(p2s * p1s);
	// The following two for loops stack correctly
	for(size_t i=0; i < p1s ;++i)
	for(size_t j=0; j < p2s ;++j)
		results[p2s*i+j] = static_cast<float>(
		std::abs(delta_phi(phis1[i]-phis2[j])));
	// This is a non-symmetric!! matrix of differences stored as RVec
	return results; 
};
//auto zw_deltaphi_cut=[](const floats& deltaphi){
//	return any_of(deltaphi.cbegin(),
//	              deltaphi  .cend(),
//	              [](float delta){return delta >= DELTA_PHI_ZW;});
//};
auto zmet_deltaphi=[](const floats& z_phi,const float& met_pt){
	floats      results(z_phi.size());
	for(size_t i=0; i < z_phi.size() ;++i)
		results[i] = std::abs(delta_phi(z_phi[i]-met_pt));
	return results;
};
//auto zmet_deltaphi_cut=[](const floats& deltaPhi){
//	return any_of(deltaPhi.cbegin(),
//	              deltaPhi  .cend(),
//	              [](float delta){return delta >= DELTA_PHI_ZMET;});
//};
auto bjet_variable
	=[](const        floats&  Jet_variable,
	    const unsigned int & nJet,
	    const          ints& lead_bjet){
	floats vec;
	for(size_t i=0; i < nJet ;++i)
		if(lead_bjet.at(i) == 1)
			vec.push_back(Jet_variable.at(i));
	return vec;
};
auto numberofbjets = [](const ints& bjets){
	const auto nbjet = std::count_if(
		bjets.begin(),
		bjets  .end(),
		[](int i){return i;});
	return nbjet;
};
auto top_reconst
	=[](const floats& bjets_pt,
	    const floats& bjets_eta,
	    const floats& bjets_phi,
	    const floats& bjets_mass,
	    const floats& wpair_pt,
	    const floats& wpair_eta,
	    const floats& wpair_phi,
	    const floats& wpair_mass){
	// This function finds the closest to top mass
	float   t_reco_mass = std::numeric_limits<float>::infinity();
	const size_t nbjets = bjets_pt.size();
	const size_t    nWs = wpair_pt.size();
	if(!all_equal(nbjets,bjets_eta.size(),bjets_phi.size(),bjets_mass.size())
	&& !all_equal(nWs   ,wpair_eta.size(),wpair_phi.size(),wpair_mass.size()))
		throw std::logic_error("Collections must be the same size");
	else if(nbjets == 0 || nWs == 0)
		throw std::logic_error("Collections must not be empty");
//	size_t bjet_index = std::numeric_limits<size_t>::max();
//	size_t    W_index = std::numeric_limits<size_t>::max();
	auto    BJets = TLorentzVector{};
	auto    RecoW = TLorentzVector{};
	auto reco_top = TLorentzVector{};
	// The following two for loops stack correctly
	for(size_t i=0; i < nbjets ;++i)
	for(size_t j=0; j < nWs    ;++j){
		BJets.SetPtEtaPhiM(bjets_pt[i],bjets_eta[i],bjets_phi[i],bjets_mass[i]);
		RecoW.SetPtEtaPhiM(wpair_pt[j],wpair_eta[j],wpair_phi[j],wpair_mass[j]);
		// The following two IF stacks correctly
		if(std::abs(RecoW.M() - W_MASS) < W_MASS_CUT)
		if (float reco_mass = static_cast<float>((RecoW + BJets).M());
		    std::abs(TOP_MASS-reco_mass) < std::abs(TOP_MASS-t_reco_mass)){
			t_reco_mass = reco_mass; // found closer to top mass 
			reco_top = RecoW + BJets;
		}
	}
	return reco_top;
};
auto TLVex(const PtEtaPhiM what){
	return [=](const TLorentzVector& object){
		float result;
		switch(what){
			case  pt:{result = static_cast<float>(object .Pt());break;}
			case eta:{result = static_cast<float>(object.Eta());break;}
			case phi:{result = static_cast<float>(object.Phi());break;}
			case   m:{result = static_cast<float>(object  .M());break;}
			default :{throw std::invalid_argument(
				"TLorentzVector extraction not recognised");}
		}
		return result;
	};
}
auto BTaggedEffGiver(TH2D* ratio){
	return [=](const floats& pts,const floats& etas){
		if(!all_equal(pts.size(),etas.size()))
			throw std::logic_error("Collections must be the same size");
		else if(pts.empty())
			throw std::logic_error("Collections must not be empty");
		floats BTaggedEff;
		for(size_t   i=0; i < pts.size() ;++i){
			int  PtBin = ratio->GetXaxis()->FindBin( pts[i]);
			int EtaBin = ratio->GetYaxis()->FindBin(etas[i]);
			float  eff = static_cast<float>(ratio->GetBinContent(PtBin,EtaBin));
			if(eff != 0) BTaggedEff.push_back(eff);//what if eff is zero? check with kathryn.
		}
		return BTaggedEff;
	};
}
auto EffIsBTaggedProduct=[](const floats& EffIsBTagged){
	float      result  = 1;
	for(size_t i=0;  i < EffIsBTagged.size() ;++i)
	           result *= EffIsBTagged[i];
	return result;
};
auto EffNoBTaggedProduct=[](const floats& EffNoBTagged){
	float  result  = 1;
	for(size_t i=0;  i < EffNoBTagged.size();++i)
	       result *= 1 - EffNoBTagged[i];
	return result;
};
auto Sfi_EffIsBTaggedProduct=[](const floats& EffIsBTagged,const floats& sfi){
	float  result = 1;
	size_t b = EffIsBTagged.size(), s = sfi.size();
	if(b!=s)std::cout<<"Sfi_EffIsBTaggedProduct got diff sizes"<<std::endl;
	size_t   size = b < s ? b : s;
	for(size_t i=0; i < size ;++i)
	       result    *= sfi[i] * EffIsBTagged[i];
	return result;
};
auto Sfj_EffNoBTaggedProduct=[](const floats& EffNoBTagged,const floats& sfj){
	float  result = 1;
	size_t b = EffNoBTagged.size(), s = sfj.size();
	if(b!=s)std::cout<<"Sfj_EffNoBTaggedProduct got diff sizes"<<std::endl;
	size_t   size = b < s ? b : s;
	for(size_t i=0; i < size ;++i)
	     result  *= 1 - EffNoBTagged[i] * sfj[i];
	return result;
};
auto btag_weight = [](const float& p_data,const float& p_MC){
	float weight;
	if(p_data != 0 && p_MC != 0)
		weight  = p_data / p_MC;
	else
		weight  = 0;
	return weight;
};
////////////// SCALE FACTORS /////////////
// Normalization * btag weights
//auto tzq_sf = [](const float& b){ return    TZQ_W*b; };
//auto  ww_sf = [](const float& b){ return WWLNQQ_W*b; };
//auto  wz_sf = [](const float& b){ return WZLNQQ_W*b; };
//auto  zz_sf = [](const float& b){ return ZZLLQQ_W*b; };
//auto ttz_sf = [](const float& b){ return  TTZQQ_W*b; };
auto sf(const dataSource ds){
	return [=](const float& b){
		float result;
		switch(ds){
			case tzq:{result =    TZQ_W;break;}
			case  ww:{result = WWLNQQ_W;break;}
			case  wz:{result = WZLNQQ_W;break;}
			case  zz:{result = ZZLLQQ_W;break;}
			case ttz:{result =  TTZQQ_W;break;}
			default :{throw std::invalid_argument("Unimplemented ds");}
		}
		return result * b;
	};
}
// Event Weight, incl. btag & Normalization
auto rvec_rep_const=[](const float& sf,const floats& iRVec){
// this function just repeats sf, for the size of iRVec
	floats weight(iRVec.size(),sf); 
	return weight;
};
} // namespace
///////////////////////////////////////////////////////////////////////////////
/*
int main(){
	//for(auto ds:dataSourceAll)
	//	calchisto(ds);
	calchisto(tzq);
	return 0;
};
*/
void calchisto(const dataSource ds){
	//ROOT::EnableImplicitMT();// parallel functioning

	// Read MC data source
	std::string temp_header="/data/disk0/nanoAOD_2017/",
	temp_opener,temp_footer="/*.root"; /**/
	switch(ds){ // tzq has a different disk !
	case tzq:{temp_opener="/data/disk3/nanoAOD_2017/tZqlvqq/*.root";break;} /**/
	case  ww:{temp_opener=temp_header+  "WWToLNuQQ"    +temp_footer;break;}
	case  wz:{temp_opener=temp_header+  "WZTo1L1Nu2Q"  +temp_footer;break;}
	case  zz:{temp_opener=temp_header+  "ZZTo2L2Q"     +temp_footer;break;}
	case ttz:{temp_opener=temp_header+ "ttZToQQ"       +temp_footer;break;}
	default :{throw std::invalid_argument("Unimplemented ds");}
	}
	ROOT::RDataFrame dssbdfc{"Events",temp_opener};

	// Read a chain of exptData
	TChain elnuEvents("Events");
	TChain munuEvents("Events");
	temp_footer = "/*.root"; // just to be sure
	temp_header =
		"/data/disk3/nanoAOD_2017/SingleElectron_NanoAOD25Oct2019_Run";
	for(std::string c:{"B","C","D","E","F"}){ // guaranteed sequential
		temp_opener=temp_header+ c +temp_footer;
		elnuEvents.Add(temp_opener.c_str());
	}
	temp_header="/data/disk3/nanoAOD_2017/SingleMuon_NanoAOD25Oct2019_Run";
	for(std::string c:{"B","C","D","E","F"}){ // guaranteed sequential
		temp_opener=temp_header+ c +temp_footer;
		munuEvents.Add(temp_opener.c_str());
	}
	ROOT::RDataFrame metEvents{"Events" , // channels unified?
		"/data/disk0/nanoAOD_2017/MET*/*.root"};
	ROOT::RDataFrame elnudfc(elnuEvents); // TODO: CMS and MET
	ROOT::RDataFrame munudfc(munuEvents); // if channels unified still make two
	auto dssbdf=dssbdfc.Range(0,100); // make test runs faster by restriction
	auto elnudf=elnudfc.Range(0,100); // real run should not Range
	auto munudf=munudfc.Range(0,100);

	auto elnu_event_selection = dssbdf //data source signal background data frame
	.Define("tight_els"    ,el_sel(EL_TIGHT_ID),
	   {"Electron_isPFcand","Electron_pt","Electron_eta","Electron_cutBased"})
	.Define("tight_el_pt"  ,select<floats>, {"Electron_pt" , "tight_els"})
	.Define("tight_el_eta" ,select<floats>, {"Electron_eta", "tight_els"})
	.Define("tight_el_phi" ,select<floats>, {"Electron_phi", "tight_els"})
	.Define("loose_els"    ,el_sel(EL_LOOSE_ID),
	   {"Electron_isPFcand","Electron_pt","Electron_eta","Electron_cutBased"})
	.Define("loose_el_pt"  ,select<floats>, {"Electron_pt" , "loose_els"})
	.Define("MET___phi_sel",{"MET_phi"})
	.Define("MET_el_pt_sel",{"MET_pt" })
	.Filter(lep_met_selector(MET_EL_PT),{"MET_el_pt_sel"}, "MET PT CUT")
	.Filter(lep_cut, {"tight_el_pt", "loose_el_pt"}, "lepton cut");

	auto elnu_w_selection
	   = elnu_event_selection
	.Define("w_el_pt"  ,get_w_lep_string(elnu,pt )) // TODO: Check if [tight_eles] came from nanoAOD
	.Define("w_el_eta" ,get_w_lep_string(elnu,eta))
	.Define("w_el_phi" ,get_w_lep_string(elnu,phi))
	.Define("w_el_mass",transverse_w_mass, {"w_el_pt","w_el_phi",
	        "MET_el_pt_sel","MET___phi_sel"})
	.Define("elnu_invmass" ,lep_nu_invmass, {"tight_el_pt","tight_el_eta",
	        "tight_el_phi" ,"Electron_mass",
	         "CaloMET_phi" ,"CaloMET_pt", "CaloMET_sumEt"});
	//.Filter(w_mass_cut, {"w_el_mass"}, "W mass cut");

	// There is a Histogram1D done on w_selection

	auto elnu_jets_selection
	   = elnu_w_selection
	.Define("jet_el_min_dR"  , jet_lep_min_deltaR,
	       {"Jet_eta","Jet_phi","w_el_eta","w_el_phi"})
	.Define("tight_jets"     , tight_jet_id   ,
	       {"jet_el_min_dR"  , "Jet_pt","Jet_eta","Jet_jetId"})
	.Define("tight_jets_pt"  , select<floats> , {"Jet_pt"  ,"tight_jets"})
	.Define("tight_jets_eta" , select<floats> , {"Jet_eta" ,"tight_jets"})
	.Define("tight_jets_phi" , select<floats> , {"Jet_phi" ,"tight_jets"})
	.Define("tight_jets_mass", select<floats> , {"Jet_mass","tight_jets"})
	.Define("tight_jets_deltaphi",jet_deltaphi, {"tight_jets_phi"})
	.Filter(jet_cut, {"tight_jets"}, "Jet cut");
	
	auto elnu_jets_bjets_selection
	   = elnu_jets_selection
	.Define("is_bjets"         ,is_bjet_id    ,{"tight_jets","Jet_eta","Jet_btagCSVV2"})
	.Define("no_bjets"         ,no_bjet_id    ,{"tight_jets","Jet_eta"})
	.Define("is_btag_numer"    ,is_bjet_numer ,{"Jet_partonFlavour","is_bjets"})
	.Define("is_btag_denom"    ,is_bjet_denom ,{"Jet_partonFlavour","no_bjets"})
	.Define("no_btag_numer"    ,no_bjet_numer ,{"Jet_partonFlavour","is_bjets"})
	.Define("no_btag_denom"    ,no_bjet_denom ,{"Jet_partonFlavour","no_bjets"})
	.Define("is_btag_numer_pt" ,select<floats>,{"Jet_pt" ,"is_btag_numer"})
	.Define("is_btag_numer_eta",select<floats>,{"Jet_eta","is_btag_numer"})
	.Define("is_btag_denom_pt" ,select<floats>,{"Jet_pt" ,"is_btag_denom"})
	.Define("is_btag_denom_eta",select<floats>,{"Jet_eta","is_btag_denom"})
	.Define("no_btag_numer_pt" ,select<floats>,{"Jet_pt" ,"no_btag_numer"})
	.Define("no_btag_numer_eta",select<floats>,{"Jet_eta","no_btag_numer"})
	.Define("no_btag_denom_pt" ,select<floats>,{"Jet_pt" ,"no_btag_denom"})
	.Define("no_btag_denom_eta",select<floats>,{"Jet_eta","no_btag_denom"})
	.Define("sfi",is_btag_CSVv2,{"Jet_btagCSVV2","tight_jets_pt","tight_jets_eta"})
	.Define("sfj",no_btag_CSVv2,{"Jet_btagCSVV2","tight_jets_pt","tight_jets_eta"})
	.Filter(bjet_cut, {"is_bjets"}, "b jet cut");
	
	auto elnu_jec_selection
	   = elnu_jets_bjets_selection
	.Define("pt_resol",jet_smear_pt_resol, {"tight_jets_pt",
	            "tight_jets_eta","fixedGridRhoFastjetAll"})
	.Define("Sjer",      jet_smear_Sjer  , {"tight_jets_eta"})
	.Define("cjer",    delta_R_jet_smear , {"tight_jets_pt","GenJet_pt",
	    "pt_resol","Sjer","jet_el_min_dR"})
	.Define("jec_tight_jets_pt"   , cjer , {"tight_jets_pt"   , "cjer"})
	.Define("jec_tight_jets_eta"  , cjer , {"tight_jets_eta"  , "cjer"})
	.Define("jec_tight_jets_phi"  , cjer , {"tight_jets_phi"  , "cjer"})
	.Define("jec_tight_jets_mass" , cjer , {"tight_jets_mass" , "cjer"});
	
	auto elnu_z_rec_selection
	   = elnu_jec_selection
	.Define(    "lead_bjet", find_lead_mask,{"is_bjets", "Jet_pt"})
	.Define(  "z_reco_jets", find_z_pair   ,{"jec_tight_jets_pt",
	   "jec_tight_jets_phi","jec_tight_jets_eta","jec_tight_jets_mass",
	           "tight_jets","lead_bjet"})
	.Define(  "z_pair_pt"  , select<floats>,{"jec_tight_jets_pt"  ,"z_reco_jets"})
	.Define(  "z_pair_eta" , select<floats>,{"jec_tight_jets_eta" ,"z_reco_jets"})
	.Define(  "z_pair_phi" , select<floats>,{"jec_tight_jets_phi" ,"z_reco_jets"})
	.Define(  "z_pair_mass", select<floats>,{"jec_tight_jets_mass","z_reco_jets"})
	.Define(  "z_mass"     , inv_mass,
	       {  "z_pair_pt"  ,"z_pair_eta","z_pair_phi","z_pair_mass"})
	.Define(  "z_el_min_dR", jet_lep_min_deltaR,
	       {  "z_pair_eta" ,"z_pair_phi","tight_el_eta","tight_el_phi"})
	.Define(  "zw_deltaphi",   zw_deltaphi, {"z_pair_phi","w_el_phi"})
	.Define("zmet_deltaphi", zmet_deltaphi,
	       {"z_pair_phi", "MET_el_pt_sel"});
	//.Filter(deltaR_z_l,{"z_e_min_dR"}, "delta R ZL")
	//.Filter(ZW_deltaphi_cut, {"ZW_deltaphi"}, "delta phi ZW cut")
	//.Filter(ZMet_deltaphi_cut, {"ZMet_deltaphi"}, "Z met cut ");
	//.Filter(z_mass_cut, {"z_mass"}, "z mass cut");
	
	auto elnu_brec_selection
	   = elnu_z_rec_selection
	.Define("bjetpt"  ,bjet_variable ,{"Jet_pt"  ,"nJet","lead_bjet"})
	.Define("bjeteta" ,bjet_variable ,{"Jet_eta" ,"nJet","lead_bjet"})
	.Define("bjetphi" ,bjet_variable ,{"Jet_phi" ,"nJet","lead_bjet"})
	.Define("bjetmass",bjet_variable ,{"Jet_mass","nJet","lead_bjet"})
	.Define("nbjets"  ,numberofbjets ,{"is_bjets"});
	
	auto elnu_top_selection
	   = elnu_brec_selection
	.Define("recoTop" ,top_reconst,
	       {"bjetpt"  ,"bjeteta" , "bjetphi" ,  "bjetmass",
	        "w_el_pt" ,"w_el_eta", "w_el_phi", "w_el_mass"})
	.Define("Top_Pt"  ,TLVex(pt ),{"recoTop"})
	.Define("Top_Eta" ,TLVex(eta),{"recoTop"})
	.Define("Top_Phi" ,TLVex(phi),{"recoTop"})
	.Define("Top_Mass",TLVex( m ),{"recoTop"});
	
	temp_header = "_elnu_tzq"; //this string will be used for all the enu channel therefore added.
	temp_footer = " pt vs eta in electron-neutrino channel for tZq"; //used in the histograms
	auto h_elnu_mass
	   = elnu_w_selection
	.Histo1D({"elnu_invmass", "elnu_invmass",50,0,200},"elnu_invmass");// lepton and neutrino system invariant mass histogram
	auto h_elnu_is_btag_numer_PtVsEta
	   = elnu_top_selection
	.Histo2D({static_cast<const char*>((        "is_numer" + temp_header).c_str()),
	          static_cast<const char*>(("MC is btag numer" + temp_header).c_str()),
	          50,0,400,50,-3,3},
	          "is_btag_numer_pt",
	          "is_btag_numer_eta");
	auto h_elnu_no_btag_numer_PtVsEta
	   = elnu_top_selection
	.Histo2D({static_cast<const char*>((        "no_numer" + temp_header).c_str()),
	          static_cast<const char*>(("MC no btag numer" + temp_header).c_str()),
	          50,0,400,50,-3,3},
	          "no_btag_numer_pt",
	          "no_btag_numer_eta");
	auto h_elnu_is_btag_denom_PtVsEta
	   = elnu_top_selection
	.Histo2D({static_cast<const char*>((        "is_denom" + temp_header).c_str()),
	          static_cast<const char*>(("MC is btag denom" + temp_header).c_str()),
	          50,0,400,50,-3,3},
	          "is_btag_denom_pt",
	          "is_btag_denom_eta");
	auto h_elnu_no_btag_denom_PtVsEta
	   = elnu_top_selection
	.Histo2D({static_cast<const char*>((        "no_denom" + temp_header).c_str()),
	          static_cast<const char*>(("MC no btag denom" + temp_header).c_str()),
	          50,0,400,50,-3,3},
	          "no_btag_denom_pt",
	          "no_btag_denom_eta");
//	auto h_events_is_btag_PtVsEta_canvas
//	= new TCanvas("is b tag pt Vs eta" , "is b tag pt Vs eta",10,10,900,900);
	TH2D *is_btag_ratio = new TH2D("ei", "is b tag ei",50,0,400,50,-3,3);
	is_btag_ratio = static_cast<TH2D*>(h_elnu_is_btag_numer_PtVsEta->Clone());
//	is_btag_ratio->GetXaxis()->SetTitle( "is b tag pt");
//	is_btag_ratio->GetYaxis()->SetTitle( "is b tag eta");
	is_btag_ratio->Divide(  h_elnu_is_btag_denom_PtVsEta.GetPtr());
//	h_events_is_btag_PtVsEta_canvas->BuildLegend();
//	is_btag_ratio->Draw("COLZ");
//	h_events_is_btag_PtVsEta_canvas->SaveAs("h_events_is_btag_PtVsEta_canvas.root");
//	h_events_is_btag_PtVsEta_canvas->SaveAs("h_events_is_btag_PtVsEta_canvas.pdf");

//	auto h_events_no_btag_PtVsEta_canvas
//	= new TCanvas("no b tag pt Vs eta" , "no b tag pt Vs eta",10,10,900,900);
	TH2D *no_btag_ratio = new TH2D("ej", "no b tag ei",50,0,400,50,-3,3);
	no_btag_ratio = static_cast<TH2D*>(h_elnu_no_btag_numer_PtVsEta->Clone());
//	no_btag_ratio->GetXaxis()->SetTitle( "no b tag pt");
//	no_btag_ratio->GetYaxis()->SetTitle( "no b tag eta");
	no_btag_ratio->Divide(  h_elnu_no_btag_denom_PtVsEta.GetPtr());
//	h_events_no_btag_PtVsEta_canvas->BuildLegend();
//	no_btag_ratio->Draw("COLZ");
//	h_events_no_btag_PtVsEta_canvas->SaveAs("h_events_no_btag_PtVsEta_canvas.root");
//	h_events_no_btag_PtVsEta_canvas->SaveAs("h_events_no_btag_PtVsEta_canvas.pdf");
/*
	auto IsBTaggedBinFunction
		=[&is_btag_ratio](const floats& pts,const floats& etas){
		if(!all_equal(pts.size(),etas.size()))
			throw std::logic_error("Collections must be the same size");
		else if(pts.empty())
			throw std::logic_error("Collections must not be empty");
		floats IsBTaggedEff;
		for(int i=0; i < pts.size() ;++i){
			int  PtBin = is_btag_ratio->GetXaxis()->FindBin( pts[i]);
			int EtaBin = is_btag_ratio->GetYaxis()->FindBin(etas[i]);
			float  eff = is_btag_ratio->GetBinContent(PtBin,EtaBin);
			if(eff != 0) IsBTaggedEff.push_back(eff);
		}
		return IsBTaggedEff;
	};

	auto NoBTaggedBinFunction
		=[&no_btag_ratio](const floats& pts,const floats& etas){
		if(!all_equal(pts.size(),etas.size()))
			throw std::logic_error("Collections must be the same size");
		else if(pts.empty())
			throw std::logic_error("Collections must not be empty");
		floats NoBTaggedEff;
		for(int i=0; i < pts.size() ;++i){
			int  PtBin = no_btag_ratio->GetXaxis()->FindBin( pts[i]);
			int EtaBin = no_btag_ratio->GetYaxis()->FindBin(etas[i]);
			float  eff = no_btag_ratio->GetBinContent(PtBin,EtaBin);
			if(eff != 0) NoBTaggedEff.push_back(eff);
		}
		return NoBTaggedEff;
	};
*/
	auto elnu_btag_eff
	   = elnu_top_selection
	.Define("IsEffBTagged",BTaggedEffGiver(is_btag_ratio),{"tight_jets_pt","tight_jets_eta"})// remember to use the jec version
	.Define("NoEffBTagged",BTaggedEffGiver(no_btag_ratio),{"tight_jets_pt","tight_jets_eta"});

	//auto h_d_enu_events_btag_eff
	//   = d_enu_btag_eff
	// .Histo1D({    "MC btag EFF",    "MC btag EFF electro and neutrino channel",50,0,400},   "EffBTagged");
	//auto h_d_enu_events_non_btag_eff
	//   = d_enu_btag_eff
	// .Histo1D({"MC non btag EFF","MC non btag EFF electro and neutrino channel",50,0,400},"NonEffBTagged");

	auto elnu_P_btag
	   = elnu_btag_eff
	.Define("Pi_ei"  ,    EffIsBTaggedProduct,{"IsEffBTagged"})
	.Define("Pi_ej"  ,    EffNoBTaggedProduct,{"NoEffBTagged"})
	.Define("Pi_sfei",Sfi_EffIsBTaggedProduct,{"IsEffBTagged","sfi"})
	.Define("Pi_sfej",Sfj_EffNoBTaggedProduct,{"NoEffBTagged","sfj"})
	.Define("P_MC"   ,  "Pi_ei * Pi_ej"  )
	.Define("P_Data" ,"Pi_sfei * Pi_sfej")
	.Define("btag_w" ,btag_weight,{"P_Data","P_MC"})
	.Define("sf"     ,sf(ds)     ,{"btag_w"})
	.Define("nw_tight_el_pt"     ,rvec_rep_const,{"sf","tight_el_pt"    })
	.Define("nw_tight_el_eta"    ,rvec_rep_const,{"sf","tight_el_eta"   })
	.Define("nw_tight_jets_pt"   ,rvec_rep_const,{"sf","tight_jets_pt"  })
	.Define("nw_tight_jets_eta"  ,rvec_rep_const,{"sf","tight_jets_eta" })
	.Define("nw_tight_jets_phi"  ,rvec_rep_const,{"sf","tight_jets_phi" })
	.Define("nw_tight_jets_mass" ,rvec_rep_const,{"sf","tight_jets_mass"})
	.Define("nw_jet_el_min_dR"   ,rvec_rep_const,{"sf",  "jet_el_min_dR"})
	.Define(  "nw_z_el_min_dR"   ,rvec_rep_const,{"sf",    "z_el_min_dR"})
	.Define(  "nw_w_el_mass"     ,rvec_rep_const,{"sf",    "w_el_mass"  })
	.Define(    "nw_z_mass"      ,"sf") // nw_z_mass is just one value, = sf
	.Define("nw_tight_jets_deltaphi",rvec_rep_const,{"sf","tight_jets_deltaphi"})
	.Define(      "nw_zmet_deltaphi",rvec_rep_const,{"sf",      "zmet_deltaphi"})
	.Define(        "nw_zw_deltaphi",rvec_rep_const,{"sf",        "zw_deltaphi"});

	//Write histogram to a root file:
	temp_opener  = "elnu_";
	switch(ds){
		case tzq:{temp_opener+="tzq";break;}
		case  ww:{temp_opener+="ww_";break;}
		case  wz:{temp_opener+="wz_";break;}
		case  zz:{temp_opener+="zz_";break;}
		case ttz:{temp_opener+="ttz";break;}
		default :{throw std::invalid_argument("Unimplemented ds");}
	}
	temp_opener += ".histo";
	TFile outfile(static_cast<const char*>(temp_opener.c_str()),"RECREATE");

	h_elnu_mass->Write();
	h_elnu_is_btag_numer_PtVsEta->Write();
	h_elnu_no_btag_numer_PtVsEta->Write();
	h_elnu_is_btag_denom_PtVsEta->Write();
	h_elnu_no_btag_denom_PtVsEta->Write();

	outfile.Close();
/*
        auto h_events_is_btag_numer_PtVsEta_canvas = new TCanvas("is b tag pt Vs eta", "is b tag pt Vs eta",10,10,900,900);

        h_d_enu_events_is_btag_numer_PtVsEta->GetXaxis()->SetTitle("is b tag numer pt");
        h_d_enu_events_is_btag_numer_PtVsEta->GetYaxis()->SetTitle("is b tag numer eta");

        h_events_is_btag_numer_PtVsEta_canvas->BuildLegend();
        h_d_enu_events_is_btag_numer_PtVsEta->Draw("COLZ");

        h_events_is_btag_numer_PtVsEta_canvas->SaveAs("hist_is_btag_numer_PtVsEta.root");
        h_events_is_btag_numer_PtVsEta_canvas->SaveAs("hist_is_btag_numer_PtVsEta.pdf");


        auto h_events_no_btag_numer_PtVsEta_canvas = new TCanvas("no b tag pt Vs eta", "no b tag pt Vs eta",10,10,900,900);

        h_d_enu_events_no_btag_numer_PtVsEta->GetXaxis()->SetTitle("no b tag numer pt");
        h_d_enu_events_no_btag_numer_PtVsEta->GetYaxis()->SetTitle("no b tag numer eta");

        h_events_no_btag_numer_PtVsEta_canvas->BuildLegend();
        h_d_enu_events_no_btag_numer_PtVsEta->Draw("COLZ");

        h_events_no_btag_numer_PtVsEta_canvas->SaveAs("hist_no_btag_numer_PtVsEta.root");
        h_events_no_btag_numer_PtVsEta_canvas->SaveAs("hist_no_btag_numer_PtVsEta.pdf");


        auto h_events_btag_denom_PtVsEta_canvas = new TCanvas("b tag pt Vs eta", "b tag pt Vs eta",10,10,900,900);

        h_d_enu_events_btag_denom_PtVsEta->GetXaxis()->SetTitle("b tag denom pt");
        h_d_enu_events_btag_denom_PtVsEta->GetYaxis()->SetTitle("b tag denom eta");

        h_events_btag_denom_PtVsEta_canvas->BuildLegend();
        h_d_enu_events_btag_denom_PtVsEta->Draw("COLZ");

        h_events_btag_denom_PtVsEta_canvas->SaveAs("hist_btag_denom_PtVsEta.root");
        h_events_btag_denom_PtVsEta_canvas->SaveAs("hist_btag_denom_PtVsEta.pdf");


        auto h_events_non_btag_denom_PtVsEta_canvas = new TCanvas("non b tag pt Vs eta", "non b tag pt Vs eta",10,10,900,900);

        h_d_enu_events_non_btag_denom_PtVsEta->GetXaxis()->SetTitle("non b tag denom pt");
        h_d_enu_events_non_btag_denom_PtVsEta->GetYaxis()->SetTitle("non b tag denom eta");

        h_events_non_btag_denom_PtVsEta_canvas->BuildLegend();
        h_d_enu_events_non_btag_denom_PtVsEta->Draw("COLZ");

        h_events_non_btag_denom_PtVsEta_canvas->SaveAs("hist_non_btag_denom_PtVsEta.root");
        h_events_non_btag_denom_PtVsEta_canvas->SaveAs("hist_non_btag_denom_PtVsEta.pdf");
*/
/*	auto h_d_enu_Pi_ei = d_enu_P_btag.Histo1D({"Pi ei histogram","Pi ei histogram",50,0,400},"Pi_ei"); h_d_enu_Pi_ei->Write();
	auto h_d_enu_Pi_ej = d_enu_P_btag.Histo1D({"Pi ej histogram","Pi ej histogram",50,0,400},"Pi_ej"); h_d_enu_Pi_ej->Write();
	auto h_d_enu_sfi = d_enu_P_btag.Histo1D({"sfi histogram","sfi histogram",50,0,400},"sfi"); h_d_enu_sfi->Write();
	auto h_d_enu_sfj = d_enu_P_btag.Histo1D({"sfj histogram","sfj histogram",50,0,400},"sfj"); h_d_enu_sfj->Write();
	auto h_d_enu_Pi_sfei = d_enu_P_btag.Histo1D({"Pi sfei histogram","Pi sfei histogram",50,0,400},"Pi_sfei"); h_d_enu_Pi_sfei->Write();
	auto h_d_enu_Pi_sfej = d_enu_P_btag.Histo1D({"Pi sfej histogram","Pi sfej histogram",50,0,400},"Pi_sfej"); h_d_enu_Pi_sfej->Write();
	auto h_d_enu_P_MC = d_enu_P_btag.Histo1D({"P(MC) histogram","P(MC) histogram",50,0,400},"P_MC"); h_d_enu_P_MC->Write();
	auto h_d_enu_P_Data = d_enu_P_btag.Histo1D({"P(Data) histogram","P(Data) histogram",50,0,400},"P_Data"); h_d_enu_P_Data->Write();
	auto btag_w = d_enu_P_btag.Histo1D({"btag w","btag w",50,0,300},"btag_w"); btag_w->Write();
        auto histo_jetsmearing_pt_resol = d_enu_P_btag.Histo1D({"i am test","i am test",50,0,10}, "pt_resol");histo_jetsmearing_pt_resol->Write();
        auto histo_jetsmearing_Sjer = d_enu_P_btag.Histo1D({"i am test","i am test", 50,0,10}, "Sjer");histo_jetsmearing_Sjer->Write();
        auto histo_jetsmearing_deltaR = d_enu_P_btag.Histo1D({"i am test", "i am test", 50,0,10}, "cjer");histo_jetsmearing_deltaR->Write();
*/
	
}
