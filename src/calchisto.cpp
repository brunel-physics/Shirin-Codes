// TODO: fix the break segmentation what comes after inv_mass?
#include <ROOT/RDataFrame.hxx>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TChain.h>

#include <boost/algorithm/string.hpp>

#include "csv.h"
#include "calchisto.hpp"
#include "eval_complex.hpp"

using doubles = ROOT::VecOps::RVec<double>;
using  floats = ROOT::VecOps::RVec<float>;
using    ints = ROOT::VecOps::RVec<int>;
using   bools = ROOT::VecOps::RVec<bool>;
using   chars = ROOT::VecOps::RVec<UChar_t>;// aka 1 byte ints
using strings = ROOT::VecOps::RVec<std::string>;

namespace{
constexpr int    debug = 1;
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
//constexpr float  MU_LPT_MIN   = 26.f;// Leading
constexpr double MU_ETA_MAX   = 2.4;
constexpr float  MU_LOOSE_ISO = .15f;
constexpr float  MU_TIGHT_ISO = .25f;

//constexpr double MET__PT_MIN = 40;
constexpr double MET_EL_PT   = 80;
constexpr double MET_MU_PT   = 40;

constexpr double   Z_MASS     = 91.1876;
constexpr float    Z_MASS_CUT = 20.f;
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

auto met_pt_cut(const channels ch){
	float lower_bound;
	switch(ch){
		case elnu:{lower_bound = MET_EL_PT;break;}
		case munu:{lower_bound = MET_MU_PT;break;}
		default  :throw std::invalid_argument(
		"Unimplemented ch (met_pt_cut)");
	}
	return [=](const float&  met_lep_pt_sel)->bool
	   {return lower_bound < met_lep_pt_sel;};
}
auto lep_sel(    const channels ch){
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
			default  :throw std::invalid_argument(
			                "Unimplemented ch (lep_sel)");
		}
	};
}
auto lep_tight_cut(const channels ch){
        return [=](const   ints& mask,
                   const   ints& elids,
                   const floats& isos){
		bool result;
		      if(ch==elnu){
			ints   temp = elids[mask];
			result = temp.size() == 1;// Choosing 1 Tight Lepton
			result = result && temp[0] >= EL_TIGHT_ID;
		}else if(ch==munu){
			floats temp = isos[mask];
			result = temp.size() == 1;// TODO: WEIRD ISO
			result = result && temp[0] <= MU_TIGHT_ISO;
		}else{throw std::invalid_argument(
			"Unimplemented ch (lep_tight_cut)"
		);}
		return result;
	};
}
template<typename T,typename U>
[[gnu::const]] bool all_equal(const T& t, const U& u){return t == u;}
template<typename T,typename U,typename... Types>
[[gnu::const]] bool all_equal(const T& t, const U& u, Types const&... args)
	{return t == u && all_equal(u, args...);}
template <typename T>
[[gnu::const]] auto delta_phi(const T diff_phi){
	// This function just reduces input from [-2pi,+2pi] to [-pi,+pi]
	// Domain correctness is the user's responsibility
	// A more general function is just one fmod function away
	     if(diff_phi >  PI<T>)  return diff_phi - TPI<T>;
	else if(diff_phi < -PI<T>)  return diff_phi + TPI<T>;
	else                        return diff_phi;
}
[[gnu::const]] auto deltaR(
	const float eta1,const float phi1,const float eta2,const float phi2)
{return std::hypot(eta1-eta2,delta_phi(phi1-phi2));}
auto transverse_w_mass(const floats& lep__pt,
                       const floats& lep_phi,
                       const float & met__pt,
                       const float & met_phi){
	if(!all_equal(lep__pt.size(),lep_phi.size()))
		throw std::logic_error("Collections must be the same size in tWm");
	else if(lep__pt.empty())
		throw std::logic_error("Collections must not be empty in tWm");
	floats   w_mass_vec(lep__pt.size());
	for(size_t i=0; i < lep__pt.size() ;++i){
		w_mass_vec[i] = static_cast<float>(
		                std::sqrt(2*lep__pt[i]*met__pt
		          *(1-cos(delta_phi(lep_phi[i]-met_phi)))));
		/*
		if(std::abs(W_MASS-reco_mass)<std::abs(W_MASS-w_reco_mass)){
			w_reco_mass = reco_mass;
			l_index_1 = i;
		}*/
	}
	return w_mass_vec;
}
auto lep_nu_invmass(const floats& lep_pt    ,
                    const floats& lep_eta   ,
                    const floats& lep_phi   ,
                    const floats& lep_mass  ,
                    const float & cal_metpt ,
                    const float & cal_metphi){//,
                  //const float & cal_metEt ){// cal_metEt is unused
	// this function computes the invariant mass of charged lepton
	// and neutrino system, in order to calculate the W mass later on.
	TLorentzVector lep,neu;
	if(!all_equal(lep_pt.size(),lep_eta.size(),lep_phi.size(),lep_mass.size()))
		throw std::logic_error("Collections must be the same size in lep-nu invmass");
	else if(lep_pt.size() != 1)
		throw std::logic_error("Should always only have one lepton in lep-nu invmass");
	float lepnu_mass = 0.f;
	for(size_t i=0; i < lep_pt.size()  ;++i){
		lep.SetPtEtaPhiM(lep_pt[i],lep_eta[i],lep_phi[i],lep_mass[i]);
		neu.SetPtEtaPhiM(cal_metpt,lep_eta[i],cal_metphi,0.);
		lepnu_mass=static_cast<float>((lep+neu).M());
	}
	return lepnu_mass;
}
/*
auto w_mass_cut = [](const floats& w_mass){
	return std::abs(w_mass-W_MASS)<W_MASS_CUT;
};*/
auto z_mass_cut = [](const float& z_mass){
	return std::abs(z_mass-Z_MASS)<Z_MASS_CUT;
};
// TODO: top_mass_cut belongs here
auto jet_lep_min_deltaR(const floats& jet_etas,
                        const floats& jet_phis,
                        const floats& lep_etas,
                        const floats& lep_phis){
	if(!all_equal(lep_etas.size(),lep_phis.size()))
		throw std::logic_error(
		      "Collections must be the same size (lep)(jet-lep dR)");
	else if(lep_phis.empty())
		throw std::logic_error(
		      "Collections must not be empty for (lep)(jet-lep dR)");
	if(!all_equal(jet_etas.size(),jet_phis.size()))
		throw std::logic_error(
		      "Collections must be the same size (jet)(jet-lep dR)");
	else if(lep_phis.empty())
		throw std::logic_error(
		      "Collections must not be empty for (jet)(jet-lep dR)");
	floats min_dRs;
	std::transform(
		jet_etas.begin(),
		jet_etas  .end(),
		jet_phis.begin(),
		std::back_inserter(min_dRs),
		[&](float jet_eta, float jet_phi){// TODO: unwrap single lep
			return deltaR(jet_eta,jet_phi,lep_etas[0],lep_phis[0]);});
	return min_dRs;
}
auto tight_jet_id(const floats& jet_lep_min_dRs,
                  const floats& pts,
                  const floats& etas,
                  const   ints& ids){
	if(debug>0) std::cout<<"tight_jet_id"<<std::endl;
	return pts>JET__PT_MIN&&abs(etas)<JET_ETA_MAX&&jet_lep_min_dRs>JET_ISO&&ids>=2;
}
auto jet_deltaphi(const floats& phis){
	if(debug>0) std::cout<<"jet deltaphi"<<std::endl;
	floats deltaphis;// half of anti-symmetric matrix has n(n-1)/2 elements
	deltaphis.reserve((phis.size()*(phis.size()-1))/2);
	// The following two for loops stack correctly
	for(size_t   i=0; i < phis.size()-1 ;++i)
	for(size_t j=i+1; j < phis.size()   ;++j)
		deltaphis.push_back(  static_cast<float>(
			std::abs(delta_phi(phis[i]-phis[j]))));
	return deltaphis;
}
auto jets_gen_select(const floats& gen, const floats& jet){
	if(debug>0) std::cout<<"gen select"<<std::endl;
// select jets that have generated level info; used @ JEC
// use this function ONLY for Monte Carlo, not CMS / MET
	if(gen.size()  < jet.size()){
		floats good_jets = jet;//       deep copy
		good_jets.resize(gen.size());// discard tail
		return good_jets;
	}else	return   jet;
}
auto jetCutter(const unsigned jmin, const unsigned jmax){
	if(debug>0) std::cout<<"jet cut"<< jmin << jmax <<std::endl;
	return[=](const ints& jets){
		const auto nj = std:: count_if(
			jets.begin(),
			jets  .end(),
			[](int i){return i;});// 0 is false
		return jmin <= nj && nj <= jmax;
	};
}
/*
auto  jet_cut(const ints& tight_jets){
	if(debug>0) std::cout<<"jet cut"<<std::endl;
	const auto njet = std::count_if(
		tight_jets.begin(),
		tight_jets  .end(),
		[](int i){return i;});// 0 is false
	return  njet >=  JETS_MIN &&  njet <=  JETS_MAX;
}
auto bjet_cut(const ints& bjets){
	if(debug>0) std::cout<<"bjet cut"<<std::endl;
	const auto nbjet = std::count_if(
		bjets.begin(),
		bjets  .end(),
		[](int i){return i;});// 0 is false
	return nbjet >= BJETS_MIN && nbjet <= BJETS_MAX;
}*/
auto jet_smear_pt_resol
	(const floats& pt,const floats& eta,const float& rho){
	if(debug>0) std::cout<<"jet smear pt resol"<<std::endl;
	float min_eta,max_eta,min_rho,max_rho;
	float min_pt,max_pt;
	int   z;// CSV provided parameter we are not using
	float a,b,c,d;
	floats  resol(pt.size(),0.f);// vector of 0.
	if(!all_equal(pt.size(),eta.size()))
		throw std::logic_error("Collections must be the same size in jet_smear");
	else if(pt.empty())
		throw std::logic_error("Collections must not be empty in jet_smear");
	io::CSVReader<11> thisCSVfile("Fall17_V3_MC_PtResolution_AK4PFchs.txt");
	while(thisCSVfile.read_row(
	      min_eta,max_eta,min_rho,max_rho,z,min_pt,max_pt,a,b,c,d)){
		// the following if for if stacks correctly
		if(rho > min_rho && rho < max_rho)
		for(size_t i=0; i < pt.size() ;++i)
		if(eta[i] > min_eta && eta[i] < max_eta
		&&  pt[i] > min_pt  &&  pt[i] < max_pt){
			resol[i] += std::sqrt(  a * std::abs(a)/(pt[i]*pt[i])
			                    + b*b * std::pow(pt[i],d) + c*c);
		}
	}// no need to close file after this while loop
	return resol;
}
auto jet_smear_Sjer(const floats& eta){
	if(debug>0) std::cout<<"jet smear sjer"<<std::endl;
	float min_eta,max_eta;
	int   z;// unwanted
	float central_SF,SF_dn,SF_up;
	floats     Sjer(eta.size(),0.f);// vector of zeroes
	io::CSVReader<6> thisCSVfile("Fall17_V3_MC_SF_AK4PF.txt");
	while(thisCSVfile.read_row(min_eta,max_eta,z,central_SF,SF_dn,SF_up)){
		for(size_t i=0; i  <  eta.size() ;++i)
		   if(      eta[i] >  min_eta && eta[i] < max_eta)
		           Sjer[i] += central_SF;
	}
	return Sjer;
}
[[gnu::const]] auto ramp(const float Sjer){
	if(Sjer > 0.f) return Sjer;
	else           return  0.f;
}
auto delta_R_jet_smear(const floats& pt,
                       const floats& gen_pt,
                       const floats& resol,
                       const floats& Sjer,
                       const floats& deltaR){
	if(debug>0) std::cout<<"delta r jet smear"<<std::endl;
	if(!all_equal(pt.size(),resol.size(),Sjer.size()))
		throw std::logic_error("Collections must be the same size in deltaR_Jsmear");
	else if(gen_pt.size() < pt.size())
		throw std::logic_error("gen_pt shorter than pt");
	else if(pt.empty())
		throw std::logic_error("Collections must not be empty in deltaR_Jsmear");
	else if(pt.size() > deltaR.size())
		throw std::logic_error("deltaR in Delta_R_jet_smear lacks sufficient data");
	const size_t  size = pt.size();
	floats cjers(       size ,0.f);// correction factor
	for(size_t i=0; i < size ;++i){
		if(deltaR[i] < RconeBy2
		&& std::abs( pt[i]-gen_pt[i]) < 3*resol[i]*pt[i]){
			cjers[i] += (1+(1+Sjer[i])
			     * (( pt[i]-gen_pt[i])/pt[i]));
		}
		else{// needs TRandom3 library// TODO: why -ve?
			float Normdist = static_cast<float>(gRandom->Gaus(0,Sjer[i]));
			float  max_val = Sjer[i] * Sjer[i] - 1;
			cjers[i] += (1+Normdist*std::sqrt(ramp(max_val)));
		}
	}
	return cjers;
}
/*
auto cjer(const floats& jet,const floats& cjer){
	if(debug>0) std::cout<<"cjer"<<std::endl;
	floats  weighted(    jet.size());// to use on jets 4-momenta, PtEtaPhiM
	if(!all_equal(       jet.size(),cjer.size()))
		throw std::logic_error("Collections must be the same size in Cjer");
	else if(jet.empty())
		throw std::logic_error("Collections must not be empty in Cjer");
	for(size_t i=0; i <  jet.size() ;++i) weighted[i] = jet[i]*cjer[i];
	return weighted;
}*/
auto   	is_bjet_id(const floats& etas,const floats& btags){// added jec_eta
		ints is_bjets(etas.size(),0);// vector of zero
		for(size_t i=0;i<etas.size();++i)// etas size <= 6
		if((btags[i]>BTAG_DISC_MIN)&&(abs(etas[i])<BJET_ETA_MAX))is_bjets[i]+=1;
	return is_bjets;
}
auto      no_bjet_id(const floats& etas)
	{return abs(etas) < BJET_ETA_MAX;}

auto      is_bjet_numer(const ints& id,const ints& is_bjet){
	ints   is_bjet_numer = id == 5;
	       is_bjet_numer.resize(is_bjet.size());// discard tail
	return is_bjet_numer;
}
auto      no_bjet_numer(const ints& id,const ints& is_bjet){
// using bjets which has satisfied is btag conditions
	auto aid = abs(id);
	aid = (( aid >  0  && aid <= 4)
	      || aid == 21 || aid != 5);
	aid.resize(is_bjet.size());// discard tail
	return aid;
}
auto      is_bjet_denom(const ints& id,const ints& no_bjet){
// using no_bjet_id particles not matching btag criteria
	ints   is_bjet_denom = id == 5;
	       is_bjet_denom.resize(no_bjet.size());// discard tail
	return is_bjet_denom;
}
auto      no_bjet_denom(const ints& id,const ints& no_bjet){
// using bjets which has satisfied no btag condition
	auto aid = abs(id);
	aid = (( aid >  0  && aid <= 4)
	      || aid == 21 || aid != 5);
	aid.resize(no_bjet.size());// discard tail
	return aid;
}
auto btagCSVv2(const bool    check_CSV){
   return  [=](const floats& btag,
               const floats& pt,
               const floats& eta){
		bool b;// magic btag checker; heavily reused
		strings formulae(pt.size() ,"0");// vector of "0"
		floats   results(pt.size());
		if(!all_equal(pt.size(),eta.size()))
			throw std::logic_error("Collections must be the same size in btagCSVv2");
		else if(pt.empty())
			throw std::logic_error("Collections must not be empty in btagCSVv2");
		else if(btag.size() < pt.size())
			throw std::logic_error("insufficient btagCSVv2");
		std::string  measure_type,sys_type,rawFormula;
		float CSVv2  ,jet_flav,
		       pt_min,  pt_max,
		      eta_min, eta_max,
		      CSV_min, CSV_max;
		io::CSVReader<11> thisCSVfile("CSVv2_94XSF_V2_B_F.csv");
		// The following nests too much, so we do not indent
		// Each blank line means nesting deeper
		while(thisCSVfile.read_row(CSVv2,measure_type,sys_type,jet_flav,
		      eta_min,eta_max,pt_min,pt_max,CSV_min,CSV_max,rawFormula)){
		
		// always true if dont check CSV
		b = (!check_CSV) || BTAG_DISC_MIN <= CSVv2;
		if(measure_type == "comb" && b
		&& sys_type  == "central" && jet_flav == 0.f){
		
		for(size_t i=0; i < pt.size() ;++i){
		
		// always true if dont check CSV
		b = (!check_CSV)
		||(btag[i] > CSV_min && btag[i] < CSV_max);
		std::string tempFormula = rawFormula;
		if( eta[i] > eta_min &&  eta[i] < eta_max
		&&   pt[i] > pt_min  &&   pt[i] <  pt_max
		&&   b){
		
		if(formulae[i] == "0"){// only 1st found wins
		
		if(tempFormula.find("x") != std::string::npos){
		
		boost::replace_all(tempFormula,"x",std::to_string(pt[i]));
		formulae[i] = tempFormula;
		}}}}}
		}// No need to close file after this while loop.
		// resume indentation
		for(size_t j=0; j < formulae.size() ;++j){
			Eval ev;
			results[j] = static_cast<float>(
			ev.eval(      const_cast<char*>(
			        formulae[j].c_str())).real());
		}
		return results;
	};
}
auto find_lead_mask(const ints& mask,const floats& vals){
	if(debug>0) std::cout<<"lead mask"<<std::endl;
	if(!all_equal(mask.size(),vals.size()))
		throw std::logic_error("Collections must be the same size in lead_mask");
	else if(mask.empty())
		throw std::logic_error("Collections must not be empty in lead_mask");
	const auto      masked_vals = mask * vals;
	ints  lead_mask(masked_vals .size(),0);// vector of zeroes
	const auto max_idx = static_cast<size_t>(std::distance(
		masked_vals.begin(),
		max_element( masked_vals.begin(),
		             masked_vals  .end())));// Leading bjet == bjet with highest pt.
	lead_mask[max_idx] = 1;
	return lead_mask;
}
/*
auto bjet_variable(const floats& Jet_variable,
                   const   ints& lead_bjet){
	// this function needs to be modified to include jec
	if(debug>0) std::cout<<"b jet variable"<<std::endl;
	floats vec;// TODO: input length equality checks
	for(size_t i=0; i < Jet_variable.size() ;++i)
		if(lead_bjet[i] != 0)
			vec.push_back(Jet_variable[i]);
	return vec;
}
auto numberofbjets(const ints& bjets){
	if(debug>0) std::cout<<"number of bjets"<<std::endl;
	const auto nbjet = std::count_if(
		bjets.begin(),
		bjets  .end(),
		[](int i){return i;});
	return nbjet;
}*/
auto find_z_pair(const floats& pts,
                 const floats& etas,
                 const floats& phis,
                 const floats& ms,
                 const   ints& tight_jets,
                 const   ints& lead_bjet){
	if(debug>0) std::cout<<"find z pair"<<std::endl;
	// This function finds the pair nearest to z mass
	double z_reco_mass = std::numeric_limits<double>::infinity();
	size_t jet_index_1 = std::numeric_limits<size_t>::max();
	size_t jet_index_2 = std::numeric_limits<size_t>::max();
	const size_t njets = pts.size();
	if(!all_equal(pts.size(),etas.size(),phis.size(),ms.size()))
		throw std::logic_error("Collections must be the same size in Z-pair");
	else if(    pts.size()==0)
		throw std::logic_error("Collections must not be empty in Z-pair");
	ints z_pair(njets, 0);// vector of zeroes
	// The next two for loops stack correctly with the if
	for(size_t   i=0; i < njets-1 ;++i){
	for(size_t j=i+1; j < njets   ;++j){
	if(tight_jets[i] != 0 || tight_jets[j] != 0// TODO: logic??
	||  lead_bjet[i] != 1 ||  lead_bjet[j] != 1){
		TLorentzVector jet1,jet2;
		jet1.SetPtEtaPhiM(pts[i],etas[i],phis[i],ms[i]);
		jet2.SetPtEtaPhiM(pts[j],etas[j],phis[j],ms[j]);
		if (const double reco_mass=(jet1+jet2).M();
		std::abs(Z_MASS-reco_mass) < std::abs(Z_MASS-z_reco_mass)){
			z_reco_mass = reco_mass;// found nearer pair to z mass
			jet_index_1 = i;
			jet_index_2 = j;
		}
	}}}
	z_pair[jet_index_1] = 1;
	z_pair[jet_index_2] = 1;
	if(debug>1) std::cout<<"z pair"<<z_pair<<std::endl;
	return z_pair;
}
[[gnu::const]] auto inv_mass(
	const floats& pts,const floats& etas,const floats& phis,const floats& ms){
	if(debug>0) std::cout<<"inv mass"<<std::endl;
	if(!all_equal(pts.size(),etas.size(),phis.size(),ms.size()))
		throw std::logic_error("Collections must be the same size in inv_mass");
	else if(pts.empty())
		throw std::logic_error("Collections must not be empty in inv_mass");
	TLorentzVector vec,p;
	for(size_t i=0; i < pts.size() ;++i){
		p.SetPtEtaPhiM(pts[i],etas[i],phis[i],ms[i]);
		vec += p;
	}
	if(debug>0) std::cout<<"inv mass "<< static_cast<float>(vec.M()) <<std::endl;
	return static_cast<float>(vec.M());
}
auto zw_deltaphi(const floats& phis1,const floats& phis2){
	if(debug>0) std::cout<<"zw deltaphi"<<std::endl;
	const size_t   p1s = phis1.size();
	const size_t   p2s = phis2.size();
	floats results(p2s * p1s);
	// The following two for loops stack correctly
	for(size_t i=0; i < p1s ;++i)
	for(size_t j=0; j < p2s ;++j)
		results[p2s*i+j] = static_cast<float>(
		std::abs(delta_phi(phis1[i]-phis2[j])));
	// This is a non-symmetric!! matrix of differences stored as 1D RVec
	return results;
}
auto zmet_deltaphi(const floats& z_phi,const float& met_pt){
	if(debug>0) std::cout<<"zmet deltaphi"<<std::endl;
	floats      results(z_phi.size());
	for(size_t i=0; i < z_phi.size() ;++i)
		results[i] = std::abs(delta_phi(z_phi[i]-met_pt));
	return results;
}
/*
auto zw_deltaphi_cut(const floats& deltaphi){
	return any_of(deltaphi.cbegin(),
	              deltaphi  .cend(),
	              [](float delta){return delta >= DELTA_PHI_ZW;});
}*/
/*
auto zmet_deltaphi_cut(const floats& deltaPhi){
	return any_of(deltaPhi.cbegin(),
	              deltaPhi  .cend(),
	              [](float delta){return delta >= DELTA_PHI_ZMET;});
}*/
auto top_reconst(const floats& bjets_pt,
                 const floats& bjets_eta,
                 const floats& bjets_phi,
                 const floats& bjets_mass,
                 const floats& wpair_pt,// TODO: rename
                 const floats& wpair_eta,
                 const floats& wpair_phi,
                 const floats& wpair_mass){
	if(debug>0) std::cout<<"top recnst"<<std::endl;
	// This function finds the closest to top mass
	float   t_reco_mass = std::numeric_limits<float>::infinity();
	const size_t nbjets = bjets_pt.size();
	const size_t    nWs = wpair_pt.size();
	if(!all_equal(nbjets,bjets_eta.size(),bjets_phi.size(),bjets_mass.size())
	&& !all_equal(nWs   ,wpair_eta.size(),wpair_phi.size(),wpair_mass.size()))
		throw std::logic_error("Collections must be the same size in top_reconst");
	else if(nbjets == 0 || nWs == 0)
		throw std::logic_error("Collections must not be empty in top_reconst");
//	size_t bjet_index = std::numeric_limits<size_t>::max();
//	size_t    W_index = std::numeric_limits<size_t>::max();
	TLorentzVector BJets,RecoW,reco_top;
	// The following two for loops stack correctly
	for(size_t i=0; i < nbjets ;++i)
	for(size_t j=0; j < nWs    ;++j){
		BJets.SetPtEtaPhiM(bjets_pt[i],bjets_eta[i],bjets_phi[i],bjets_mass[i]);
		RecoW.SetPtEtaPhiM(wpair_pt[j],wpair_eta[j],wpair_phi[j],wpair_mass[j]);
		// The following two IF stacks correctly
		if(std::abs(RecoW.M() - W_MASS) < W_MASS_CUT)
		if (float reco_mass = static_cast<float>((RecoW + BJets).M());
		    std::abs(TOP_MASS-reco_mass) < std::abs(TOP_MASS-t_reco_mass)){
			t_reco_mass = reco_mass;// found closer to top mass
			reco_top = RecoW + BJets;
		}
	}
	return reco_top;
}
auto TLVex(   const PtEtaPhiM         what){
	return [=](const TLorentzVector& object){
		if(debug>0) std::cout<<"TLVex"<<std::endl;
		float result;
		switch(what){
			case  pt:{result = static_cast<float>(object .Pt());break;}
			case eta:{result = static_cast<float>(object.Eta());break;}
			case phi:{result = static_cast<float>(object.Phi());break;}
			case   m:{result = static_cast<float>(object  .M());break;}
			default :throw std::invalid_argument(
				"TLorentzVector extraction not recognised");
		}
		return result;
	};
}
auto BTaggedEffGiver(TH2D* ratio){
	return [=](const floats& pts,const floats& etas){
		if(debug>0) std::cout<<"bt eff giver"<<std::endl;
		if(!all_equal(pts.size(),etas.size()))
			throw std::logic_error(
			      "Collections must be the same size (EffGiver)");
		else if(pts.empty())
			throw std::logic_error(
			      "Collections must not be empty in  (EffGiver)");
		floats BTaggedEff;
		for(size_t   i=0; i < pts.size() ;++i){
			int  PtBin = ratio->GetXaxis()->FindBin( pts[i]);
			int EtaBin = ratio->GetYaxis()->FindBin(etas[i]);
			float  eff = static_cast<float>(ratio->GetBinContent(PtBin,EtaBin));
			if(eff != 0) BTaggedEff.push_back(eff);
			// what if eff==0? check with kathryn
		}
		return BTaggedEff;
	};
}
auto EffIsBTaggedProduct(const floats& EffIsBTagged){
	float      result  = 1;
	for(size_t i=0;  i < EffIsBTagged.size() ;++i)
	           result *= EffIsBTagged[i];
	return result;
}
auto EffNoBTaggedProduct(const floats& EffNoBTagged){
	float  result  = 1;
	for(size_t i=0;  i < EffNoBTagged.size();++i)
	       result *= 1 - EffNoBTagged[i];
	return result;
}
auto Sfi_EffIsBTaggedProduct(const floats& EffIsBTagged,const floats& sfi){
	float  result = 1;
	size_t b = EffIsBTagged.size(), s = sfi.size();
	if(b!=s)std::cout<<"Sfi_EffIsBTaggedProduct got diff sizes"<<std::endl;
	size_t   size = b < s ? b : s;
	for(size_t i=0; i < size ;++i)
	       result    *= sfi[i] * EffIsBTagged[i];
	return result;
}
auto Sfj_EffNoBTaggedProduct(const floats& EffNoBTagged,const floats& sfj){
	float  result = 1;
	size_t b = EffNoBTagged.size(), s = sfj.size();
	if(b!=s)std::cout<<"Sfj_EffNoBTaggedProduct got diff sizes"<<std::endl;
	size_t   size = b < s ? b : s;
	for(size_t i=0; i < size ;++i)
	     result  *= 1 - EffNoBTagged[i] * sfj[i];
	return result;
}
	////////////// SCALE FACTORS /////////////
auto btag_weight(const float& p_data,const float& p_MC){
	float  weight = p_data / p_MC;
/*	if(p_MC != 0.f)// TODO: cannot compare float to exact
	       weight = p_data / p_MC;
	else
	       weight = 0.f;*/
	if(std::isinf(weight) || std::isnan(weight)) weight = 0.f;
	return weight;
}
	// Normalization * btag weights
auto sf(const dataSource ds){
	return [=](const float& b){
		float result;
		switch(ds){
			case tzq:{result =    TZQ_W;break;}
			case  ww:{result = WWLNQQ_W;break;}
			case  wz:{result = WZLNQQ_W;break;}
			case  zz:{result = ZZLLQQ_W;break;}
			case ttz:{result =  TTZQQ_W;break;}
			default :throw std::invalid_argument("Unimplemented ds (infile)");
		}
		return result * b;
	};
}
	// Event Weight, incl. btag & Normalization
auto rep_const(const float& sf,const floats& iRVec){
	// this function just repeats sf, for the size of iRVec
	floats weight(iRVec.size(),sf);
	return weight;
}
}// namespace
void calchisto(const dataSource ds){
//	ROOT::EnableImplicitMT();// parallel functioning
	
	// Read MC data source
	std::string temp_header="/data/disk0/nanoAOD_2017/",
	temp_opener,temp_footer="/*.root";/**/
	switch(ds){// tzq and exptData use disk3!
	case tzq:{temp_opener="/data/disk3/nanoAOD_2017/tZqlvqq/*.root";break;}/**/
	case  ww:{temp_opener=temp_header+  "WWToLNuQQ"    +temp_footer;break;}
	case  wz:{temp_opener=temp_header+  "WZTo1L1Nu2Q"  +temp_footer;break;}
	case  zz:{temp_opener=temp_header+  "ZZTo2L2Q"     +temp_footer;break;}
	case ttz:{temp_opener=temp_header+ "ttZToQQ"       +temp_footer;break;}
	default :throw std::invalid_argument("Unimplemented ds (rdfopen)");
	}
	ROOT::RDataFrame dssbdfc{"Events",temp_opener};
	
	// Read a chain of exptData
	TChain elnuEvents("Events");
	TChain munuEvents("Events");
	temp_footer = "/*.root";/* just to be sure */
	temp_header =
		"/data/disk3/nanoAOD_2017/SingleElectron_NanoAOD25Oct2019_Run";
	for(std::string c:{"B","C","D","E","F"}){// guaranteed sequential
		temp_opener=temp_header+ c +temp_footer;
		elnuEvents.Add(temp_opener.c_str());
	}
	temp_header="/data/disk3/nanoAOD_2017/SingleMuon_NanoAOD25Oct2019_Run";
	for(std::string c:{"B","C","D","E","F"}){// guaranteed sequential
		temp_opener=temp_header+ c +temp_footer;
		munuEvents.Add(temp_opener.c_str());
	}
	ROOT::RDataFrame metEvents{"Events" ,// channels unified?
		"/data/disk0/nanoAOD_2017/MET*/*.root"};/**/
	ROOT::RDataFrame elnudfc(elnuEvents);// TODO: CMS and MET
	ROOT::RDataFrame munudfc(munuEvents);// if channels unified still make two
	auto dssbdf=dssbdfc.Range(0,1000000);// make test runs faster by restriction
	auto elnudf=elnudfc.Range(0,1000000);// real run should not Range
	auto munudf=munudfc.Range(0,1000000);
	
	channels ch = elnu;
	switch(ch){
		case elnu:{temp_header = "Electron_";break;}
		case munu:{temp_header =     "Muon_";break;}
		default  :throw std::invalid_argument("Unimplemented ch (rdfopen)");
	}
	auto w_selection = dssbdf
	.Filter(met_pt_cut(ch),{"MET_pt"},"MET Pt cut")
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
	.Define(   "lep__pt",temp_header+  "pt[loose_leps]")
	.Define(   "lep_eta",temp_header+ "eta[loose_leps]")
	.Define(   "lep_phi",temp_header+ "phi[loose_leps]")
	.Define(   "lep_mas",temp_header+"mass[loose_leps]")
	.Define("tw_lep__pt" ,"lep__pt")
	.Define("tw_lep_eta" ,"lep_eta")
	.Define("tw_lep_phi" ,"lep_phi")
	.Define("tw_lep_mas",transverse_w_mass,
	       {"tw_lep__pt",
	        "tw_lep_phi","MET_pt","MET_phi"})
	.Define(   "lep_nu_invmass",lep_nu_invmass,
	       {   "lep__pt",
	           "lep_eta",
	           "lep_phi",
	           "lep_mas",
	       "CaloMET_pt" ,//    TODO: unwrap RVec of 1 element
	       "CaloMET_phi"})
//	       "CaloMET_sumEt"})// TODO: add this back
//	.Filter(w_mass_cut, {"w_el_mass"}, "W mass cut")
	;
	
	// There is a Histogram1D done on w_selection
	
	auto jets_selection
	   = w_selection
	.Define("jet_lep_min_dR" , jet_lep_min_deltaR,
	       {"Jet_eta","Jet_phi","lep_eta","lep_phi"})
	.Define("tight_jets"     , tight_jet_id   ,
	       {"jet_lep_min_dR" ,  "Jet_pt","Jet_eta","Jet_jetId"})
	.Define("tight_jets__pt" ,  "Jet_pt[tight_jets]")
	.Define("tight_jets_eta" , "Jet_eta[tight_jets]")
	.Define("tight_jets_phi" , "Jet_phi[tight_jets]")
	.Define("tight_jets_mas" ,"Jet_mass[tight_jets]")
	.Define("tight_jets_deltaphi",jet_deltaphi, {"tight_jets_phi"})
	// TODO: jet_deltaphi from jec
	.Define( "good_jets__pt" ,jets_gen_select , {"GenJet_pt"  ,"tight_jets__pt"})
	.Define( "good_jets_eta" ,jets_gen_select , {"GenJet_eta" ,"tight_jets_eta"})
	.Define( "good_jets_phi" ,jets_gen_select , {"GenJet_phi" ,"tight_jets_phi"})
	.Define( "good_jets_mas" ,jets_gen_select , {"GenJet_mass","tight_jets_mas"})
	.Filter(jetCutter(JETS_MIN,JETS_MAX),{"tight_jets"},"Jet cut")
	;
	// JEC == tight_jets inc. Jet Energy Correction
	// good_jets (defined above) are tight jets which have generated level info
	// by neglecting the ones which exclude generated level info,
	// we make sure the jec is computing correctly
	// & jets w/o JEC will be excluded too.
	auto jec_selection
	   = jets_selection
	.Define("pt_resol",jet_smear_pt_resol, {"good_jets__pt",
	        "good_jets_eta","fixedGridRhoFastjetAll"})
	.Define("Sjer",      jet_smear_Sjer  , {"good_jets_eta"})
	.Define("cjer",    delta_R_jet_smear , {"good_jets__pt","GenJet_pt",
	        "pt_resol","Sjer","jet_lep_min_dR"})
	.Define("jec_jets__pt" , "good_jets__pt * cjer")
	.Define("jec_jets_eta" , "good_jets_eta * cjer")
	.Define("jec_jets_phi" , "good_jets_phi * cjer")
	.Define("jec_jets_mas" , "good_jets_mas * cjer")
//	.Define("jec_jets__pt" , cjer , {"good_jets__pt" , "cjer"})
//	.Define("jec_jets_eta" , cjer , {"good_jets_eta" , "cjer"})
//	.Define("jec_jets_phi" , cjer , {"good_jets_phi" , "cjer"})
//	.Define("jec_jets_mas" , cjer , {"good_jets_mas" , "cjer"})
	;
	auto jets_bjets_selection
	   = jec_selection
	// TODO: Fix these
	.Define("is_bjets"         ,is_bjet_id   ,{"jec_jets_eta","Jet_btagCSVV2"})
	.Define("no_bjets"         ,no_bjet_id   ,{"jec_jets_eta"})
	.Define("is_btag_numer"    ,is_bjet_numer,{"Jet_partonFlavour","is_bjets"})
	.Define("no_btag_numer"    ,no_bjet_numer,{"Jet_partonFlavour","is_bjets"})
	.Define("is_btag_denom"    ,is_bjet_denom,{"Jet_partonFlavour","no_bjets"})
	.Define("no_btag_denom"    ,no_bjet_denom,{"Jet_partonFlavour","no_bjets"})
	.Define("is_btag_numer__pt","jec_jets__pt[is_btag_numer]")
	.Define("no_btag_numer__pt","jec_jets__pt[no_btag_numer]")
	.Define("is_btag_denom__pt","jec_jets__pt[is_btag_denom]")
	.Define("no_btag_denom__pt","jec_jets__pt[no_btag_denom]")
	.Define("is_btag_numer_eta","jec_jets_eta[is_btag_numer]")
	.Define("no_btag_numer_eta","jec_jets_eta[no_btag_numer]")
	.Define("is_btag_denom_eta","jec_jets_eta[is_btag_denom]")
	.Define("no_btag_denom_eta","jec_jets_eta[no_btag_denom]")
	.Define("sfi",btagCSVv2( true),// checks CSV
	        {"Jet_btagCSVV2","jec_jets__pt","jec_jets_eta"})
	.Define("sfj",btagCSVv2(false),// ignore CSV
	        {"Jet_btagCSVV2","jec_jets__pt","jec_jets_eta"})
	.Define("lead_bjet"     ,find_lead_mask,{"is_bjets","jec_jets__pt"})
	.Define(     "bjet__pt" ,"jec_jets__pt[lead_bjet]")// Leading bjet
	.Define(     "bjet_eta" ,"jec_jets_eta[lead_bjet]")// 4-momentum
	.Define(     "bjet_phi" ,"jec_jets_phi[lead_bjet]")// used for top
	.Define(     "bjet_mas" ,"jec_jets_mas[lead_bjet]")// reconstruction
//	.Define(     "bjet__pt" ,bjet_variable,{"jec_jets__pt","lead_bjet"})
//	.Define(     "bjet_eta" ,bjet_variable,{"jec_jets_eta","lead_bjet"})
//	.Define(     "bjet_phi" ,bjet_variable,{"jec_jets_phi","lead_bjet"})
//	.Define(     "bjet_mas" ,bjet_variable,{"jec_jets_mas","lead_bjet"})
//	.Define(    "nbjets"    ,numberofbjets,{"is_bjets"})
	.Filter(jetCutter(BJETS_MIN,BJETS_MAX), {"is_bjets"}, "b jet cut")
	// TODO: The three highest pts should be kept
	;
	auto z_rec_selection
	   = jets_bjets_selection
	.Define( "z_reco_jets"  , find_z_pair  ,
	       {"jec_jets__pt"  ,
	        "jec_jets_eta"  ,
	        "jec_jets_phi"  ,
	        "jec_jets_mas"  , "tight_jets" ,"lead_bjet" })
	.Define(  "z_pair__pt"  , "jec_jets__pt[z_reco_jets]")
	.Define(  "z_pair_eta"  , "jec_jets_eta[z_reco_jets]")
	.Define(  "z_pair_phi"  , "jec_jets_phi[z_reco_jets]")
	.Define(  "z_pair_mas"  , "jec_jets_mas[z_reco_jets]")
	.Define(  "z_mass"      , inv_mass     ,
	       {  "z_pair__pt"  ,
	          "z_pair_eta"  ,
	          "z_pair_phi"  ,
	          "z_pair_mas" })
	.Define(  "z_lep_min_dR", jet_lep_min_deltaR,
	       {  "z_pair_eta"  ,
	          "z_pair_phi"  ,
	             "lep_eta"  ,
	             "lep_phi" })
	.Define(  "zw_deltaphi" ,   zw_deltaphi,{"z_pair_phi","tw_lep_phi"})
	.Define("zmet_deltaphi" , zmet_deltaphi,{"z_pair_phi","MET_pt"})
//	.Filter(deltaR_z_l,{"z_e_min_dR"}, "delta R ZL")
//	.Filter(zw_deltaphi_cut, {"zw_deltaphi"}, "delta phi ZW cut")
//	.Filter(zmet_deltaphi_cut, {"zmet_deltaphi"}, "Z met cut ");
	.Filter(z_mass_cut, {"z_mass"}, "z mass cut")
	;
	auto top_selection
	   = z_rec_selection
	.Define(  "recoTop" ,top_reconst,
	       {  "bjet__pt",
	          "bjet_eta",
	          "bjet_phi",
	          "bjet_mas",
	        "tw_lep__pt",
	        "tw_lep_eta",
	        "tw_lep_phi",
	        "tw_lep_mas"})
	.Define("ttop__pt"  ,TLVex(pt ),{"recoTop"})
	.Define("ttop_eta"  ,TLVex(eta),{"recoTop"})
	.Define("ttop_phi"  ,TLVex(phi),{"recoTop"})
	.Define("ttop_mas"  ,TLVex( m ),{"recoTop"})
	;
	// now we make the histogram names and titles
	switch(ch){
		case elnu:{temp_header = "_elnu";temp_opener="elnu_";
		           temp_footer = " electron-neutrino";break;}
		case munu:{temp_header = "_munu";temp_opener="munu_";
		           temp_footer =     " muon-neutrino";break;}
		default  :throw std::invalid_argument(
		                "Unimplemented ch (hist titles)");
	}
	temp_footer = " pt vs eta in" + temp_footer + " channel for ";
	switch(ds){
		case tzq:{temp_header+="_tzq";temp_footer+="tZq";break;}
		case  ww:{temp_header+="__ww";temp_footer+=" WW";break;}
		case  wz:{temp_header+="__wz";temp_footer+=" WZ";break;}
		case  zz:{temp_header+="__zz";temp_footer+=" ZZ";break;}
		case ttz:{temp_header+="_ttz";temp_footer+="ttZ";break;}
		default :throw std::invalid_argument(
		                "Unimplemented ds (hist titles)");
	}
	
	std::string stempName,  stmpTitle;
	char       *ctempName, *ctmpTitle;
	stempName =          "invmass" + temp_header;
	stmpTitle =          "invmass" + temp_header;// not mistake
	ctempName = new char[stempName.length()+1];
	ctmpTitle = new char[stmpTitle.length()+1];
	std::strcpy(ctempName,stempName.c_str());
	std::strcpy(ctmpTitle,stmpTitle.c_str());
	auto h_invmass// lepton-neutrino invariant mass histogram
	   = w_selection
	.Histo1D({ctempName,ctmpTitle,50,0,200},"lep_nu_invmass");
	delete[] ctempName; delete[] ctmpTitle;
	
	stempName =         "is_numer" + temp_header;
	stmpTitle = "MC is btag numer" + temp_footer;
	ctempName = new char[stempName.length()+1];
	ctmpTitle = new char[stmpTitle.length()+1];
	std::strcpy(ctempName,stempName.c_str());
	std::strcpy(ctmpTitle,stmpTitle.c_str());
	auto h_is_btag_numer_PtVsEta
	   = top_selection
	.Histo2D({ctempName,ctmpTitle,50,0,400,50,-3,3},
	          "is_btag_numer__pt",
	          "is_btag_numer_eta");
	delete[] ctempName; delete[] ctmpTitle;
	
	stempName =         "no_numer" + temp_header;
	stmpTitle = "MC no btag numer" + temp_footer;
	ctempName = new char[stempName.length()+1];
	ctmpTitle = new char[stmpTitle.length()+1];
	std::strcpy(ctempName,stempName.c_str());
	std::strcpy(ctmpTitle,stmpTitle.c_str());
	auto h_no_btag_numer_PtVsEta
	   = top_selection
	.Histo2D({ctempName,ctmpTitle,50,0,400,50,-3,3},
	          "no_btag_numer__pt",
	          "no_btag_numer_eta");
	delete[] ctempName; delete[] ctmpTitle;
	
	stempName =         "is_denom" + temp_header;
	stmpTitle = "MC is btag denom" + temp_footer;
	ctempName = new char[stempName.length()+1];
	ctmpTitle = new char[stmpTitle.length()+1];
	std::strcpy(ctempName,stempName.c_str());
	std::strcpy(ctmpTitle,stmpTitle.c_str());
	auto h_is_btag_denom_PtVsEta
	   = top_selection
	.Histo2D({ctempName,ctmpTitle,50,0,400,50,-3,3},
	          "is_btag_denom__pt",
	          "is_btag_denom_eta");
	delete[] ctempName; delete[] ctmpTitle;
	
	stempName =         "no_denom" + temp_header;
	stmpTitle = "MC no btag denom" + temp_footer;
	ctempName = new char[stempName.length()+1];
	ctmpTitle = new char[stmpTitle.length()+1];
	std::strcpy(ctempName,stempName.c_str());
	std::strcpy(ctmpTitle,stmpTitle.c_str());
	auto h_no_btag_denom_PtVsEta
	   = top_selection
	.Histo2D({ctempName,ctmpTitle,50,0,400,50,-3,3},
	          "no_btag_denom__pt",
	          "no_btag_denom_eta");
	delete[] ctempName; delete[] ctmpTitle;
	
	TH2D *
	is_btag_ratio = new TH2D("ei", "is b tag ei",50,0,400,50,-3,3);
	is_btag_ratio = static_cast<TH2D*>(h_is_btag_numer_PtVsEta->Clone());
	is_btag_ratio->Divide(             h_is_btag_denom_PtVsEta.GetPtr());
	
	TH2D *
	no_btag_ratio = new TH2D("ej", "no b tag ei",50,0,400,50,-3,3);
	no_btag_ratio = static_cast<TH2D*>(h_no_btag_numer_PtVsEta->Clone());
	no_btag_ratio->Divide(             h_no_btag_denom_PtVsEta.GetPtr());
	
	auto btag_eff
	   = top_selection
	.Define("IsEffBTagged",BTaggedEffGiver(is_btag_ratio),
	       {"jec_jets__pt","jec_jets_eta"})// remember to use the jec version
	.Define("NoEffBTagged",BTaggedEffGiver(no_btag_ratio),
	       {"jec_jets__pt","jec_jets_eta"})
	;
/*	auto h_is_btag_eff
	   =      btag_eff
	 .Histo1D({     "MC btag EFF",// TODO: switch cases
	                "MC btag EFF electron-neutrino channel",
	            50,0,400},"IsEffBTagged");
	auto h_no_btag_eff
	   =      btag_eff
	 .Histo1D({ "MC no btag EFF",
	            "MC no btag EFF electron-neutrino channel",
	           50,0,400},"NoEffBTagged");
*/
	auto P_btag
	   = btag_eff
	.Define("Pi___ei",     EffIsBTaggedProduct,{"IsEffBTagged"})
	.Define("Pi___ej",     EffNoBTaggedProduct,{"NoEffBTagged"})
	.Define("Pi_sfei", Sfi_EffIsBTaggedProduct,{"IsEffBTagged","sfi"})
	.Define("Pi_sfej", Sfj_EffNoBTaggedProduct,{"NoEffBTagged","sfj"})
	.Define("P_MC"   , "Pi___ei * Pi___ej")
	.Define("P_Data" , "Pi_sfei * Pi_sfej")
	.Define("btag_w" , btag_weight  ,{"P_Data","P_MC"})
	.Define("sf"     , sf(ds)       ,{"btag_w"})
	.Define("nw_lep__pt"            ,rep_const,{"sf",     "lep__pt"   })
	.Define("nw_lep_eta"            ,rep_const,{"sf",     "lep_eta"   })
	.Define("nw_tight_jets__pt"     ,rep_const,{"sf","jec_jets__pt"   })
	.Define("nw_tight_jets_eta"     ,rep_const,{"sf","jec_jets_eta"   })
	.Define("nw_tight_jets_phi"     ,rep_const,{"sf","jec_jets_phi"   })
	.Define("nw_tight_jets_mas"     ,rep_const,{"sf","jec_jets_mas"   })
	.Define("nw_jet_lep_min_dR"     ,rep_const,{"sf", "jet_lep_min_dR"})
	.Define( "nw__z_lep_min_dR"     ,rep_const,{"sf",   "z_lep_min_dR"})
	.Define( "nw_tw_lep_mas"        ,rep_const,{"sf",  "tw_lep_mas"   })
	.Define( "nw__z_mass"           ,"sf")// is just one value, == sf
	.Define("nw_tight_jets_deltaphi",rep_const,{"sf","tight_jets_deltaphi"})
	.Define(      "nw_zmet_deltaphi",rep_const,{"sf",      "zmet_deltaphi"})
	.Define(        "nw_zw_deltaphi",rep_const,{"sf",        "zw_deltaphi"})
//	.Define("nw_ttop__pt"  ,"sf")
//	.Define("nw_ttop_mas","sf")
	;
	// Testing top mass histo
//	auto h_transTopmass = P_btag.Histo1D({
//		static_cast<const char*>((       "ttopmass" + temp_header).c_str()),
//		static_cast<const char*>(("trans. Top mass" + temp_header).c_str()),
//		50,0,200},
//		"ttop__pt","nw_ttop__pt");
	auto 
	h_tWmVsZmass_calc = P_btag.Histo2D({
		static_cast<const char*>(("tWmVsZmass" + temp_header).c_str()),
		static_cast<const char*>(("tWmVsZmass" + temp_header).c_str()),
		50,0,150,50,0,150},
		"tw_lep_mas","z_mass");
	h_tWmVsZmass_calc->GetXaxis()->SetTitle("tWm   GeV/C^2");
	h_tWmVsZmass_calc->GetYaxis()->SetTitle("Zmass GeV/C^2");

	// write histograms to a root file
	switch(ch){
		case elnu:{temp_opener ="elnu_";break;}
		case munu:{temp_opener ="munu_";break;}
		default  :throw std::invalid_argument(
		                "Unimplemented ch (outfile)");
	}
	switch(ds){
		case tzq: {temp_opener+=  "tzq";break;}
		case  ww: {temp_opener+=  "_ww";break;}
		case  wz: {temp_opener+=  "_wz";break;}
		case  zz: {temp_opener+=  "_zz";break;}
		case ttz: {temp_opener+=  "ttz";break;}
		default :throw std::invalid_argument(
		               "Unimplemented ds (outfile)");
	}
	temp_opener += ".histo";
	TFile hf(static_cast<const char*>(temp_opener.c_str()),"RECREATE");
	
	h_invmass->Write();
	h_is_btag_numer_PtVsEta->Write();
	h_no_btag_numer_PtVsEta->Write();
	h_is_btag_denom_PtVsEta->Write();
	h_no_btag_denom_PtVsEta->Write();
//	h_transTopmass->Write();
	h_tWmVsZmass_calc->Write();
	
	hf.Close();
//	std::cout << "btag weight is " << std::endl;
//	<< *btag_w << std::endl;
//	<< P_btag.Take<float>("btag_w").GetValue()[0] << std::endl;
}
