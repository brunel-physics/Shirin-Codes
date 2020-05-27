// TODO: Adding Triggers
// TODO: lepton trigger efficiency to be integrated inside calchisto
// TODO: PILE UP, Shape uncertainties, top pt reweighting, non promp lepton corrections
// TODO: MET unclustering correction
// TODO: Matching genjet to jets via using:Muon_genPartIdx,Electron_genPartIdx,Jet_genJetIdx
// TODO: plotsstack for all dataSources and channels

#include <ROOT/RDataFrame.hxx>//#include <ROOT/RCsvDS.hxx>
#include <TLorentzVector.h>
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
using   chars = ROOT::VecOps::RVec<UChar_t>;// aka 1 byte ints
using strings = ROOT::VecOps::RVec<std::string>;

namespace{
  constexpr    int debug = 0;
  constexpr  float ENDCAP_ETA_MIN = 1.566f;
  constexpr  float BARREL_ETA_MAX = 1.4442f;
//constexpr    int EL_MAX_NUM   = 1;
  constexpr  float EL__PT_MIN   = 45.f;//{15}//min 12, AP 45,
//constexpr  float EL_LPT_MIN   = 35.f;// Leading
  constexpr  float EL_ETA_MAX   = 2.5f;
  constexpr    int EL_LOOSE_ID  = 1;
  constexpr    int EL_TIGHT_ID  = 4;

//constexpr    int MU_MAX_NUM   = 1;
  constexpr  float MU__PT_MIN   = 40.f;//min 33, AP 40,
//constexpr  float MU_LPT_MIN   = 26.f;// Leading
  constexpr  float MU_ETA_MAX   = 2.4f;
  constexpr  float MU_LOOSE_ISO = .15f;
  constexpr  float MU_TIGHT_ISO = .25f;

//constexpr  float MET__PT_MIN  = 40.f;
  constexpr  float MET_EL_PT    = 80.f;
  constexpr  float MET_MU_PT    = 40.f;

  constexpr double   Z_MASS     =  91.1876;
  constexpr double   Z_MASS_CUT =  20.;
  constexpr double   W_MASS     =  80.385;
  constexpr double   W_MASS_CUT =  20.;
  constexpr double TOP_MASS     = 172.5;
//constexpr double TOP_MASS_CUT =  20.;

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

//constexpr double RconeBy2 =  .2;

constexpr double    TZQ_W =  .0128;
constexpr double WWLNQQ_W = 2.1740;
constexpr double WZLNQQ_W =  .2335;
constexpr double  TTZQQ_W =  .0237;
constexpr double ZZLLQQ_W =  .0485;

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

auto met_pt_cut(const channel ch){
	float lower_bound;
	switch(ch){
		case elnu:{lower_bound = MET_EL_PT;break;}
		case munu:{lower_bound = MET_MU_PT;break;}
//		default  :throw std::invalid_argument(
//			"Unimplemented ch (met_pt_cut)");
	}
	return [=](const float   met_lep_pt_sel)->bool
	   {return lower_bound < met_lep_pt_sel;};
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
	// A more general function is just one fmod function away
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

double transverse_w_mass(const double lep__pt,
                         const double lep_phi,
                         const  float met__pt,
                         const  float met_phi){
	return std::sqrt(                   2*(lep__pt)
	                                     *(met__pt)
	*(1.-cos(delta_phi(static_cast<double>(lep_phi)
	                  -static_cast<double>(met_phi)))));
}
auto lep_nu_invmass(const float lep_pt    ,
                    const float lep_eta   ,
                    const float lep_phi   ,
                    const float lep_mass  ,
                    const float cal_metpt ,
                    const float cal_metphi){//,
                  //const float cal_metEt ){// cal_metEt is unused
	// this function computes the invariant mass of charged lepton
	// and neutrino system, in order to calculate the W mass later on.
	TLorentzVector lep,neu;
	lep.SetPtEtaPhiM(   lep_pt,lep_eta,   lep_phi,lep_mass);
	neu.SetPtEtaPhiM(cal_metpt,lep_eta,cal_metphi,0.);
	return (lep+neu).M();
}
auto easy_mass_cut(const double theo,const double cut){
	return [=](const double ours){return std::abs(ours-theo)<cut;};
}
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
auto tight_jet_id(const doubles& jet_lep_min_dRs,
                  const  floats& pts,
                  const  floats& etas,
                  const    ints& ids){
	//if(0<debug) std::cout<<"tight_jet_id"<<std::endl;
	return pts>JET__PT_MIN&&abs(etas)<JET_ETA_MAX&&jet_lep_min_dRs>JET_ISO&&ids>=2;
}
auto jetCutter(const unsigned jmin, const unsigned jmax){
	if(3<debug) std::cout<<"jet cut "<< jmin <<" "<< jmax <<std::endl;
	return[=](const ints& jetmask){
		const auto nj = std:: count_if(
			jetmask.cbegin(),
			jetmask.  cend(),
			[](int i){return i;});// 0 is false
		return jmin <= nj && nj <= jmax;
	};
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
/*
template<typename T>// allow us to return w/o knowing data type
auto retVar(const T& v){return[&](){return v;};}
*/
// Jet Energy Resolution and Jet Energy Smearing
auto jet_smear_pt_resol(const floats& pt,const floats& eta,const float rho,const bool fatjet){
	//if(0<debug) std::cout<<"jet smear pt resol"<<std::endl;
	doubles resol(pt.size());
	if(!all_equal(pt.size(),eta.size())) throw std::logic_error(
		"Collections must be the same size (jet_smear_pt_resol)");
	if(pt.empty()) throw std::logic_error(
		"Collections must not be empty in  (jet_smear_pt_resol)");
	double etaMin,etaMax,rhoMin,rhoMax;
	double pt_Min,pt_Max;
	double a,b,c,d;
	std:string temp_header="aux/Fall17_V3_MC_PtResolution_AK",
		   temp_footer = "PFchs.txt",temp_opener;
	if(!fatjet)
	     temp_opener = temp_header + "4" + temp_footer;
	else{temp_opener = temp_header + "8" + temp_footer;}
	io::CSVReader<10> thisCSVfile(temp_opener.c_str());
	thisCSVfile.read_header(io::ignore_extra_column,
	"etaMin","etaMax","rhoMin","rhoMax","ptMin","ptMax","a","b","c","d");
	while(thisCSVfile.read_row(
	 etaMin , etaMax , rhoMin , rhoMax ,pt_Min ,pt_Max , a , b , c , d)){
		// The following if for if stacks correctly
		if(rhoMin < rho && rho < rhoMax)
		for(size_t i=0; i < pt.size() ;++i)
		if(etaMin < eta[i] && eta[i] < etaMax
		&& pt_Min < pt [i] && pt [i] < pt_Max){
			resol[i] += std::sqrt(  a * std::abs(a)/(pt[i]*pt[i])
			                    + b*b * std::pow(pt[i],d) + c*c);
		}
	}// No need to close file after this while loop
	return resol;//};
}
auto jet_smear_Sjer(const floats& etas,const bool fatjet){
	//if(0<debug) std::cout<<"jet smear sjer"<<std::endl;
   doubles Sjers(            etas.size());
	double  etaMin,etaMax;
	double  centralSF,dnSF,upSF;
        std:string temp_header="aux/Fall17_V3_MC_SF_AK",
                   temp_footer = "PFchs.txt",temp_opener;
        if(!fatjet)
             temp_opener = temp_header + "4" + temp_footer;
        else{temp_opener = temp_header + "8" + temp_footer;}
        io::CSVReader<10> thisCSVfile(temp_opener.c_str());
	thisCSVfile.read_header(io::ignore_extra_column,
	                          "etaMin","etaMax","centralSF","dnSF","upSF");
	while(thisCSVfile.read_row(etaMin , etaMax , centralSF , dnSF , upSF)){
		for(size_t i=0; i  <  etas.size() ;++i)
		   if(      etaMin <  etas[i] && etas[i] < etaMax)
		          Sjers[i] += centralSF;
	}
	return Sjers;//};
}
[[gnu::const]] auto ramp(const double Sjer){
	if(Sjer > 0.) return Sjer;
	else          return   0.;
}
auto delta_R_jet_smear(bool fatjet){
	return[&,fatjet](const floats& jpt,const floats& jeta,const floats& jphi,
                       const float   rho,
                       const floats& gpt,const floats& geta,const floats& gphi){
	if(0<debug) std::cout<<"delta r jet smear"<<std::endl;
	if(!all_equal(jpt.size(),jeta.size(),jphi.size())
	|| !all_equal(gpt.size(),geta.size(),gphi.size()))
		throw std::logic_error(
			"Collections must be the same size (deltaR_Jsmear)");
	if(jpt.empty()) throw std::logic_error(
			"Collections must not be empty for (deltaR_Jsmear)");
	// the method used in here is the Jets Smearing Hybrid Method
	double RconeBy2;
	if(fatjet)RconeBy2 = 0.8/2;
	else{RconeBy2 = 0.4/2;}
	const size_t size = jpt.size();
	double temp;
	doubles cjers; cjers.reserve(size);
	ints gen(gpt.size(),1); gen.resize(size);
	auto resol = jet_smear_pt_resol/*(ptRcsv)*/(jpt,jeta,rho,fatjet);
	auto  Sjer = jet_smear_Sjer   /*(sjerCsv)*/(    jeta	,fatjet);
	for(size_t i=0; i < size ;++i){
		if(0!=gen [i] && std::abs(jpt [i]-gpt [i]) < 3*resol[i]*jpt[i]
		&& deltaR(geta[i],gphi[i],jeta[i],jphi[i]) < RconeBy2){
			temp = (1+(1+Sjer[i])// Scaling method
			    *((jpt[i]-gpt[i])/jpt[i]));
		}else{// Stochastic smearing
			double Normdist = gRandom->Gaus(0,Sjer[i]);
			double  max_val = Sjer[i] * Sjer[i] - 1;
			temp = (1+Normdist*std::sqrt(ramp(max_val)));
		}
		if(temp < 0) temp = 0;
		cjers.emplace_back(temp);
	}
	return cjers;
};}
auto      is_bjet_id(const doubles& etas,const floats& btags){// added jec_eta
   ints   is_bjet_id = BTAG_DISC_MIN < btags;// A MASK!
          is_bjet_id.resize(etas.size(),0);// discard tail or pad zeroes
   return abs(etas) < BJET_ETA_MAX && is_bjet_id;// now same size,no more error
}
auto      no_bjet_id(const doubles& etas){
   return abs(etas) < BJET_ETA_MAX;}

auto      is_bjet_numer(const ints& id,const ints& is_bjet){
   ints   is_bjet_numer = abs(id) == 5;
          is_bjet_numer.resize(is_bjet.size(),0);// discard tail or pad zeroes
   return is_bjet_numer;
}
auto      no_bjet_numer(const ints& id,const ints& is_bjet){
// using bjets which has satisfied is btag conditions
	ints     aid =  abs(id);
	aid = (( aid >  0  && aid <= 4)
	      || aid == 21 || aid != 5);
	aid.resize(is_bjet.size(),0);// discard tail or pad zeroes
	return aid;
}
constexpr auto& is_bjet_denom = is_bjet_numer;// same functions, just
constexpr auto& no_bjet_denom = no_bjet_numer;// different input mask
//template<typename T>
auto btagCSVv2(const bool check_CSVv2){//,
//ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void>& btagDF){
    return [=](const  floats& btag,
               const doubles& pt  ,
               const doubles& eta ,
               const    ints& flav){
	if(0<debug)std::cout<<"btagCSVv2 entered"<<std::endl;
	const size_t     size = pt  .size();
	strings formulae(size  ,"1");//,check_CSVv2 ? "1" : "0");
	doubles  results(size );
	if(!all_equal(   size  ,flav.size(),
	             eta.size(),btag.size())) throw std::logic_error(
		"Collections must be the same size in btagCSVv2");
	if(0 == size) throw std::logic_error(
		"Collections must not be empty in btagCSVv2");
/*	if(btag.size() < pt.size()) throw std::logic_error(
		"insufficient btagCSVv2");
*/
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

	if("1" == formulae[i]){// only 1st found wins

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
	return results; };
}
template <typename T>
auto abs_deltaphi(const double Zphi,const T Wphi)
	{return std::abs(delta_phi(Zphi-Wphi));}

auto jet_deltaphi(const doubles& phis){
	if(0<debug) std::cout<<"jet deltaphi"<<std::endl;
	doubles deltaphis;// half of anti-symmetric matrix has n(n-1)/2 elements
	deltaphis.reserve((phis.size()*(phis.size()-1))/2);// reserving size
	// The following two for loops stack correctly     // yet leave
	for(size_t   i=0; i < phis.size()-1 ;++i)	         // empty.
	for(size_t j=i+1; j < phis.size()   ;++j)
		deltaphis.emplace_back(abs_deltaphi(phis[i],phis[j]));
	return deltaphis;
}
auto find_lead_mask(const doubles& vals,const ints& mask){
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
auto find_z_pair(const doubles& pts,
                 const doubles& etas,
                 const doubles& phis,
                 const doubles& ms,
                 const    ints& lead_bjet){
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
		TLorentzVector jet1,jet2;
		jet1.SetPtEtaPhiM(pts[i],etas[i],phis[i],ms[i]);
		jet2.SetPtEtaPhiM(pts[j],etas[j],phis[j],ms[j]);
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
auto TLVpairAdd(const doubles& pt__pair,// Create TLorentzV from 2 jets 4-mom
                const doubles& eta_pair,
                const doubles& phi_pair,
                const doubles& mas_pair){
	//if(0<debug) std::cout<<"TLVpairAdd"<<std::endl;
	if(2 != pt__pair.size()) throw std::logic_error(
		"Not pair of Z (TLVpairAdd)");
	TLorentzVector v,p;
	v.SetPtEtaPhiM(pt__pair[0],eta_pair[0],phi_pair[0],mas_pair[0]);
	p.SetPtEtaPhiM(pt__pair[1],eta_pair[1],phi_pair[1],mas_pair[1]);
	return v+p;
}
auto TLVex(   const PtEtaPhiM         what){
   return [=](const TLorentzVector& object){
		//if(0<debug) std::cout<<"TLVex"<<std::endl;
		double result;
		switch(what){
			case  pt:{result = object.Pt ();break;}
			case eta:{result = object.Eta();break;}
			case phi:{result = object.Phi();break;}
			case   m:{result = object. M ();break;}
//			default :throw std::invalid_argument(
//				"TLorentzVector extraction not recognised");
		}
		return result;
	};
}
/*
auto deltaphi_cut(const double    x){
      return  [=](const doubles& dps){
		return  any(dps >= x);
	};
}
*/
auto top_reconst(const doubles& bjets_pt,
                 const doubles& bjets_eta,
                 const doubles& bjets_phi,
                 const doubles& bjets_mass,
                 const double   wpair_pt,// TODO: rename
                 const double   wpair_eta,
                 const double   wpair_phi,
                 const double   wpair_mass){
	if(0<debug) std::cout<<"top recnst"<<std::endl;
	// This function finds the closest to top mass
	double     cur_mass = std::numeric_limits<double>::infinity();
	const size_t nbjets = bjets_pt.size();
	if(!all_equal(nbjets,bjets_eta.size(),bjets_phi.size(),bjets_mass.size()))
		throw std::logic_error(
		"Collections must be the same size (top_recnst)");
	if(nbjets == 0) throw std::logic_error(
		"Collections must not be empty in  (top_recnst)");
	TLorentzVector BJets,RecoW,reco_top,temp;
		RecoW.SetPtEtaPhiM(wpair_pt   ,wpair_eta   ,wpair_phi   ,wpair_mass);
	// The following if for stacks correctly
	if(std::abs(RecoW.M() - W_MASS) < W_MASS_CUT)
	for(size_t i=0; i < nbjets ;++i){
		BJets.SetPtEtaPhiM(bjets_pt[i],bjets_eta[i],bjets_phi[i],bjets_mass[i]);
		temp = RecoW + BJets;
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
	// reconstruct z avoiding bjets first
	.Define(   "fin_jets_Dph" , jet_deltaphi  ,{"fin_jets_phi"})
	.Define( "lead_bjet"      , find_lead_mask,{"fin_jets__pt","is_bjets"})
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
	.Define(   "z_TLV"	     , TLVpairAdd,
	       {   "z_pair__pt"   ,
	           "z_pair_eta"   ,
	           "z_pair_phi"   ,
	           "z_pair_mas"  })
	.Define(   "z_mas"        ,    TLVex( m ),{"z_TLV"})
	.Filter( easy_mass_cut(Z_MASS,Z_MASS_CUT),{"z_mas"},"z mass cut")
	.Define(   "z__pt"        ,    TLVex(pt ),{"z_TLV"})
	.Define(   "z_eta"        ,    TLVex(eta),{"z_TLV"})
	.Define(   "z_phi"        ,    TLVex(phi),{"z_TLV"})
	.Define(   "z_lep_min_dR" , jet_lep_min_deltaR<doubles>,// TODO: Check input?
	       {   "z_pair_eta"   ,
	           "z_pair_phi"   ,
	              "lep_eta"   ,
	              "lep_phi"  })// TODO: wrong input below
//	.Filter(      deltaR_z_l,  {  "z_lep_min_dR"},"delta R ZL")
	.Define(   "zw_Dph" , abs_deltaphi<double>,{"z_phi","tw_lep_phi"})
//	.Filter(      deltaphi_cut(DELTA_PHI_ZW),
//	       {   "zw_Dph"},"delta phi ZW cut")
	.Define( "zmet_Dph" , abs_deltaphi<float >,{"z_phi","MET_phi"})
//	.Filter(      deltaphi_cut(DELTA_PHI_ZMET),
//	       { "zmet_Dph"},"Z met cut ")
	// Adding a mask to pick the z_pair_phi from the fin_jets_phi
	.Define("z_jets_Dph", jet_deltaphi   ,{"z_pair_phi"})
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
	.Define("ttop__pt",TLVex(pt ),{"recoTtop"})
	.Define("ttop_eta",TLVex(eta),{"recoTtop"})
	.Define("ttop_phi",TLVex(phi),{"recoTtop"})
	.Define("ttop_mas",TLVex( m ),{"recoTtop"})
	;
}
// Btagging for eff i and eff j
auto BTaggedEffGiver(TH2D*& ratio,bool b){
return [&,b](const doubles& pts,const doubles& etas){
	if(0<debug) std::cout<<"bt eff giver in  "<< ratio <<std::endl;
	if(!all_equal(pts.size(),etas.size())) throw std::logic_error(
		"Collections must be the same size (effGiver)");
	if(pts.empty()) throw std::logic_error(
		"Collections must not be empty in  (effGiver)");
	doubles BTaggedEff(pts.size(),1);//pts.size(), b ? 1. : 0.);
	for(size_t   i=0; i <  pts.size() ;++i){
		int  PtBin = ratio->GetXaxis()->FindBin(         pts [i] );
		int EtaBin = ratio->GetYaxis()->FindBin(std::abs(etas[i]));
		double eff = ratio->GetBinContent(PtBin,EtaBin);
		if(FP_NORMAL == std::fpclassify(eff))// if eff non-zero/inf/NaN
			BTaggedEff[i] = eff;
		// above only pushed back nonzero nice eff
		// what do we do with eff==0? check with kathryn
		// below, if ej, and near 1., we put 0. instead
		//if(!b && std::abs(eff-1.) <= 2*std::numeric_limits<double>::epsilon())
		//	BTaggedEff[i] = 0.;
	}
	if(5<debug) std::cout<<"bt eff giver "<< BTaggedEff <<" for b "<<b<<std::endl;
	return BTaggedEff;};
}
auto EffIsBTaggedProduct(const doubles& EffIsBTagged){
	double     result  = 1.;
	for(size_t i=0; i  < EffIsBTagged.size() ;++i)
	           result *= EffIsBTagged[i];
	if(5<debug) std::cout<<"EffIsBTaggedProduct "<<result<<std::endl;
	return     result;
}
auto EffNoBTaggedProduct(const doubles& EffNoBTagged){
	double result  = 1.;
	for(size_t i=0;  i < EffNoBTagged.size();++i)
	       result *= 1.- EffNoBTagged[i];
	if(5<debug) std::cout<<"EffNoBTaggedProduct "<<result<<std::endl;
	return result;
}
auto Sfi_EffIsBTaggedProduct(const doubles& EffIsBTagged,const doubles& sfi){
	double result = 1.;
	size_t b = EffIsBTagged.size(), s = sfi.size();
	if(b!=s)std::cout<<"Sfi_EffIsBTaggedProduct got diff sizes"<<std::endl;
	size_t   size = b < s ? b : s;
	for(size_t i=0; i < size ;++i)
	       result    *= sfi[i] * EffIsBTagged[i];
	if(5<debug)std::cout<<"sfi = "<<sfi<<" product = "<<result<<std::endl;
	return result;
}
auto Sfj_EffNoBTaggedProduct(const doubles& EffNoBTagged,const doubles& sfj){
	double result = 1.;
	size_t b = EffNoBTagged.size(), s = sfj.size();
	if(b!=s)std::cout<<"Sfj_EffNoBTaggedProduct got diff sizes"<<std::endl;
//	if(__FMA__)std::cout<<"FMA not accelerated"<<std::endl;
	size_t   size = b < s ? b : s;
	for(size_t i=0; i < size ;++i)
//	       result*= -std::fma(EffNoBTagged[i],sfj[i],-1.);
	       result*= 1. - EffNoBTagged[i]*sfj[i];
	if(5<debug)std::cout<<"sfj = "<<sfj<<" product = "<<result<<std::endl;
	return result;
}
/*
auto EffIsBTaggedProduct(const doubles& IsEffBTagged){
	double result = std::reduce(//std::execution::par_unseq,
		IsEffBTagged. cbegin(),// parallel multiply
		IsEffBTagged.   cend(),// un-sequentially
		1.,std::multiplies<>()
	);
	if(1<debug) std::cout<<"EffIsBTaggedProduct "<<result<<std::endl;
	return result;
}
auto EffNoBTaggedProduct(const doubles& NoEffBTagged){
	doubles values = 1. - NoEffBTagged;
	double  result = std::reduce(//std::execution::par_unseq,
		     values.cbegin(),// parallel multiply
		     values.  cend(),// un-sequentially
		     1.,std::multiplies<>()
	);
	if(1<debug) std::cout<<"EffNoBTaggedProduct "<<result<<std::endl;
	return  result;
}
auto Sfi_EffIsBTaggedProduct(const doubles& IsEffBTagged,const doubles& sfi){
	if(IsEffBTagged.size() != sfi.size()) std::cout
		<< "Sfi_EffIsBTaggedProduct got diff sizes"  <<std::endl;
	double result;
	if(IsEffBTagged.size() <= sfi.size()){
		result = std::transform_reduce(//std::execution::par_unseq,
			IsEffBTagged. cbegin(),// parallel multiply
			IsEffBTagged.   cend(),// un-sequentially
			         sfi. cbegin(),// stops correctly
			1.,std::multiplies<>(),std::multiplies<>()
		);
	} else{
		result = std::transform_reduce(//std::execution::par_unseq,
			         sfi. cbegin(),// parallel multiply
			         sfi.   cend(),// un-sequentially
			IsEffBTagged. cbegin(),// stops correctly
			1.,std::multiplies<>(),std::multiplies<>()
		);
	}
	if(1<debug)std::cout<<"Sfi_EffIsBTaggedProduct "<<result<<std::endl;
	return result;
}
auto Sfj_EffNoBTaggedProduct(const doubles& NoEffBTagged,const doubles& sfj){
	if(NoEffBTagged.size() != sfj.size()) std::cout
		<< "Sfj_EffNoBTaggedProduct got diff sizes"  <<std::endl;
	double result;
	if(NoEffBTagged.size() <= sfj.size()){
		result = std::transform_reduce(//std::execution::par_unseq,
			NoEffBTagged. cbegin(),// parallel multiply
			NoEffBTagged.   cend(),// un-sequentially
			         sfj. cbegin(),// stops correctly
			1.,std::multiplies<>(),
			[](double x,double y){return -std::fma(x,y,-1.);}
		);// this multiplies together a lot of 1 - ej * sfj
	} else{
		result = std::transform_reduce(//std::execution::par_unseq,
			         sfj. cbegin(),// parallel multiply
			         sfj.   cend(),// un-sequentially
			NoEffBTagged. cbegin(),// stops correctly
			1.,std::multiplies<>(),
			[](double x,double y){return -std::fma(x,y,-1.);}
		);
	}
	if(1<debug)std::cout<<"Sfj_EffNoBTaggedProduct "<<result<<std::endl;
	return result;
}
*/
auto btag_weight(const double p_data,const double p_MC){
	double weight =p_data/p_MC;
	if(FP_NORMAL !=std::fpclassify(weight))weight=1.;// Rids non-zero/inf/NaN
//	if(0 < debug)  std::cout<<"btag_w is "<<weight<<std::endl;
	return  weight;
}
auto  lep_gpt(const channel ch){
   return [=](const ints &id, const floats &pt){
		int idV;
		switch(ch){
			case elnu:{idV = 11;break;}
			case munu:{idV = 13;break;}
		}
		auto gAll = pt[ abs(id) == idV ];
		if(5<debug)std::cout<<"gen pt picking from "<<gAll<<std::endl;
		return *max_element(gAll.cbegin(),gAll.cend());
   };
}
// Lepton efficiencies
auto elEffGiver(const float              pt ,
                const float              eta,
                const TH2F* const &recoLowEt,
                const TH2F* const &reco_pass,
                const TH2F* const &tight_94x){
	std::map<elSf,double> dict = {{Eff,1.},{Smr,1.}};
	// eff == electron    regression    corrections
	// smr == energy scale and smearing corrections
	if(2.5 < std::abs(eta)) return dict;
	int PtBin,EtaBin;
	if(2<debug)std::cout<<"el eff giver"<<std::endl;
	if( pt      < 20.f){
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
auto muEffGiver(const float         pt ,
                const float         eta,
                const TH2D* const &id_N,
                const TH2D* const &id_Y,
                const TH2D* const &id_A,
                const TH2D* const &id_T,
                const TH2D* const &isoN,
                const TH2D* const &isoY,
                const TH2D* const &isoA,
                const TH2D* const &isoT){
	const float ata = abs(eta);
	std::map<muSf,double> dict
	={{Id_N,1.},{Id_Y,1.},{Id_A,1.},{Id_T,1. },
	  {IsoN,1.},{IsoY,1.},{IsoA,1.},{IsoT,1.}};
	if(2<debug)std::cout<<"mu eff giver"<<std::endl;
	if(20.f <= pt && pt <= 120.f && ata <= 2.4f) return dict;
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
auto lepEffGiver(const channel      ch,
                 const bool         MC,
                 const TH2F* const &recoLowEt,
                 const TH2F* const &reco_pass,
                 const TH2F* const &tight_94x,
                 const TH2D* const &id_N,
                 const TH2D* const &id_Y,
                 const TH2D* const &id_A,
                 const TH2D* const &id_T,
                 const TH2D* const &isoN,
                 const TH2D* const &isoY,
                 const TH2D* const &isoA,
                 const TH2D* const &isoT){
      return [=](const float     pt,const float eta,
                 const float    phi,const   int Q,
                 const float gen_pt,const   int nl){
//	if(0 < debug)std::cout<< "lep eff giver"<<std::endl;
	double sf = 1., id = 1., iso = 1., eff = 1., smr = 1.;
	RoccoR rc("src/roccor.Run2.v3/RoccoR2017.txt");
	if(MC){switch(ch){
	case elnu:{
		auto  dict = elEffGiver(pt,eta,recoLowEt,reco_pass,tight_94x);
		eff = dict[Eff];
		smr = dict[Smr];
		break;}
	case munu:{
		auto  dict = muEffGiver(pt,eta
		                       ,id_N,id_Y,id_A,id_T
		                       ,isoN,isoY,isoA,isoT
		);
		id  = dict[Id_N];
		iso = dict[IsoN];// muEffGiver done; rocco follows
		if(gen_pt != 0){
			std::cout<<"rocco 1"<<std::endl;
			sf = rc.kSpreadMC(Q,pt,eta,phi,gen_pt,0,0);
		}else{
			std::cout<<"rocco 2"<<std::endl;
			auto u = gRandom->Rndm();
			std::cout<<// TODO: not sure if gRandom works!
			"Warning, u must be between 0 and 1, u is "
			<<u<<std::endl;
			sf = rc. kSmearMC(Q,pt,eta,phi,nl,u,0,0);
/* Rocco scale factor desc.
scale factors for momentum of each muon:
// data
double dtSF = rc.kScaleDT(Q, pt, eta, phi, s=0, m=0);
// (recommended),MC scale and resolution correction when gen muon exists
double mcSF = rc.kSpreadMC(Q, pt, eta, phi, genPt, s=0, m=0);
// MC scale and extra smearing when matched gen muon does not exist
double mcSF = rc.kSmearMC(Q, pt, eta, phi, nl, u, s=0, m=0);
----------------------------------------------------------
Here:
Q is charge
nl is trackerLayersWithMeasurement
u is a random number distributed uniformly between 0 and 1
(gRandom->Rndm());
s is error set    (default is 0)
m is error member (default is 0, ranges from 0 to nmembers-1)
For MC, when switching to different error sets/members for
a given muon, random number (u) should remain unchanged.
*/
		}
		break;}
	}}else{switch(ch){
		case elnu:break;
		case munu:{
			std::cout<<"rocco 3"<<std::endl;
			sf = rc. kScaleDT(Q,pt,eta,phi,0,0);
			break;}
		}
	}
	if(0 < debug) std::cout
		<< "most " << sf
		<< " id  " << id
		<< " iso " << iso
		<< " eff " << eff
		<< " smr " << smr
		<< " for " << ch << std::endl;
	return sf * id * iso * eff * smr;};
}
auto pile(const TH1D* const &PuWd,
          const TH1D* const &PuUd,
          const TH1D* const &PuDd){
return[=](const int   npv){
	std::map<puSf,double> dict = {{puW,1.},{upW,1.},{dnW,1.}};
	dict[puW] = PuWd->GetBinContent(PuWd->GetXaxis()->FindBin(npv));
	dict[upW] = PuUd->GetBinContent(PuUd->GetXaxis()->FindBin(npv));
	dict[dnW] = PuDd->GetBinContent(PuDd->GetXaxis()->FindBin(npv));
	if(-1<debug)std::cout<<"pile "<<npv<<" "
		<<dict[puW]<<" "<<dict[upW]<<" "<<dict[dnW]<<std::endl;
	return dict;
};
}
// Simulation correction Scale Factors
auto sf(const  dataSource ds){
	if(0<debug) std::cout<<"scale factor "<<std::endl;
	return [=](const double b,
	           const double mostSF){
		// TODO: lepton smearing and trigger efficiency
		double result;
		switch(ds){
			case tzq:{result =    TZQ_W;break;}
			case  ww:{result = WWLNQQ_W;break;}
			case  wz:{result = WZLNQQ_W;break;}
			case  zz:{result = ZZLLQQ_W;break;}
			case ttz:{result =  TTZQQ_W;break;}
			case met:// fall through to cms
			case cms:{result = 1.;break;}// ignore btag wt
//			default :throw std::invalid_argument(
//				"Unimplemented ds (sf)");
		}
		if(5<debug)std::cout<<"b_w "<<b
		<<" sf "<<result*b//<<std::endl;
		<<" mostSF " << mostSF<<std::endl;
		return result*b*mostSF;
	};
}
auto rep_const(const double sf,const doubles& iRVec){
	// this function just repeats sf, for the size of iRVec
	if(0<debug) std::cout<<"repeat const "<<iRVec.size()<<std::endl;
	doubles weight(iRVec.size(),sf);
	return  weight;
}
template <typename T>
auto finalScaling(const dataSource ds,T &rdf){
	return rdf
	.Define("sf"     ,  sf(ds) ,{"btag_w","mostSF"})
	. Alias("nw_lep__pt"       ,"sf")// is just one value, == sf
	. Alias("nw_lep_eta"       ,"sf")// LOL WHY SO DUMB, weigh the
	. Alias("nw_lep_phi"       ,"sf")// hist BY "sf" then!
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
	.Define("nw_z_jets_Dph"	   ,rep_const,{"sf",  "z_jets_Dph"})
	. Alias(    "nw_zmet_Dph"  ,"sf")
	. Alias(      "nw_zw_Dph"  ,"sf")
	. Alias("nw_lep_nu_invmass","sf")
	. Alias("nw_ttop__pt"      ,"sf")
	. Alias("nw_ttop_mas"      ,"sf")
	;
}
auto runLBfilter(
	const std::map<size_t,std::vector<std::pair<size_t,size_t>>>
	&runLBdict){
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
	const TH2F* const recoLowEt = static_cast<const TH2F* const>(tHf);
	tF->Close();
	tF = TFile::Open(
		"aux/elEff/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root");
	tF ->GetObject("EGamma_SF2D",tHf);tHf->SetDirectory(nullptr);
	const TH2F* const reco_pass = static_cast<const TH2F* const>(tHf);
	tF->Close();
	tF = TFile::Open(
		"aux/elEff/egammaEffi.txt_EGM2D_runBCDEF_passingTight94X.root");
	tF ->GetObject("EGamma_SF2D",tHf);tHf->SetDirectory(nullptr);
	const TH2F* const tight_94x = static_cast<const TH2F* const>(tHf);
	tF->Close();
	// muon efficiencies
	tF = TFile::Open("aux/muEff/Muon_RunBCDEF_SF_ID.root");
	tF ->GetObject("NUM_TightID_DEN_genTracks_pt_abseta",tHd);
	tHd->SetDirectory(nullptr);// make it stay even if file closed
	const TH2D* const id_N = static_cast<const TH2D* const>(tHd);
	tF ->Close();
	tF = TFile::Open("aux/muEff/Muon_RunBCDEF_SF_ID_syst.root");
	tF ->GetObject("NUM_TightID_DEN_genTracks_pt_abseta",tHd);
	tHd->SetDirectory(nullptr);
	const TH2D* const id_Y = static_cast<const TH2D* const>(tHd);
	tF ->GetObject("NUM_TightID_DEN_genTracks_pt_abseta_stat",tHd);
	tHd->SetDirectory(nullptr);
	const TH2D* const id_A = static_cast<const TH2D* const>(tHd);
	tF ->GetObject("NUM_TightID_DEN_genTracks_pt_abseta_syst",tHd);
	tHd->SetDirectory(nullptr);
	const TH2D* const id_T = static_cast<const TH2D* const>(tHd);
	tF ->Close();
	tF = TFile::Open("aux/muEff/Muon_RunBCDEF_SF_ISO.root");
	tF ->GetObject("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta",tHd);
	tHd->SetDirectory(nullptr);
	const TH2D* const isoN = static_cast<const TH2D* const>(tHd);
	tF ->Close();
	tF = TFile::Open("aux/muEff/Muon_RunBCDEF_SF_ISO_syst.root");
	tF ->GetObject("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta",tHd);
	tHd->SetDirectory(nullptr);
	const TH2D* const isoY = static_cast<const TH2D* const>(tHd);
	tF ->GetObject("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta_stat",tHd);
	tHd->SetDirectory(nullptr);
	const TH2D* const isoA = static_cast<const TH2D* const>(tHd);
	tF ->GetObject("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta_syst",tHd);
	tHd->SetDirectory(nullptr);
	const TH2D* const isoT = static_cast<const TH2D* const>(tHd);
	tF ->Close();
	// pile up
	tF = TFile::Open("aux/pileupMC.root");// denom because wrong
	tF ->GetObject("pileup",t1d);t1d->SetDirectory(nullptr);
	t1d->Scale(1.0/t1d->Integral());
	const TH1D* const PuDe = static_cast<const TH1D* const>(t1d);
	tF ->Close();
	tF = TFile::Open("aux/truePileupTest.root");
	tF ->GetObject("pileup",t1d);t1d->SetDirectory(nullptr);
	t1d->Scale(1.0/t1d->Integral());
	t1d->Divide(PuDe);
	const TH1D* const PuWd = static_cast<const TH1D* const>(t1d);
	tF ->Close();
	tF = TFile::Open("aux/truePileupUp.root");
	tF ->GetObject("pileup",t1d);t1d->SetDirectory(nullptr);
	t1d->Scale(1.0/t1d->Integral());
	t1d->Divide(PuDe);
	const TH1D* const PuUd = static_cast<const TH1D* const>(t1d);
	tF ->Close();
	tF = TFile::Open("aux/truePileupDown.root");
	tF ->GetObject("pileup",t1d);t1d->SetDirectory(nullptr);
	t1d->Scale(1.0/t1d->Integral());
	t1d->Divide(PuDe);
	const TH1D* const PuDd = static_cast<const TH1D* const>(t1d);
	tF ->Close();
	tF	= nullptr; tHf = nullptr; tHd = nullptr; t1d = nullptr;
//	std::cout<<"Auxiliary files processed"       <<std::endl;

	// Open data files even if unused
	// then automatically choose which one to read from
	// No penalty for opening and leaving unused
	// Can even open multiple times at once in parallel
	// Open MC data source EVEN IF UNUSED
	std::string temp_header="/data/disk0/nanoAOD_2017/",
	temp_opener,temp_footer="/*.root";/**/
	switch(ds){// tzq and exptData use disk3!
	case tzq:{temp_opener="/data/disk3/nanoAOD_2017/tZqlvqq/*.root";break;}/**/
	case  ww:{temp_opener=temp_header+  "WWToLNuQQ"    +temp_footer;break;}
	case  wz:{temp_opener=temp_header+  "WZTo1L1Nu2Q"  +temp_footer;break;}
	case  zz:{temp_opener=temp_header+  "ZZTo2L2Q"     +temp_footer;break;}
	case ttz:{temp_opener=temp_header+ "ttZToQQ"       +temp_footer;break;}
	case met:{temp_opener=temp_header+ "ttZToQQ"       +temp_footer;break;}
	case cms:{temp_opener=temp_header+ "ttZToQQ"       +temp_footer;break;}
//	default :throw std::invalid_argument("Unimplemented ds (rdfopen)");
	}// CMS and MET MUST do some OPENABLE file ; reject later
	ROOT::RDataFrame mc__df("Events",temp_opener);// Monte Carlo
	// Open chains of exptData EVEN IF UNUSED
	TChain elnuCMS("Events");
	TChain munuCMS("Events");
	TChain bothMET("Events");
	temp_footer = "/*.root" ;/* just to be sure */
	temp_header =
		"/data/disk3/nanoAOD_2017/SingleElectron_NanoAOD25Oct2019_Run";
	for(std::string c:{"B","C","D","E","F"}){// guaranteed sequential
		temp_opener=temp_header+ c +temp_footer;
		elnuCMS.Add(temp_opener.c_str());
	}
	temp_header="/data/disk3/nanoAOD_2017/SingleMuon_NanoAOD25Oct2019_Run";
	for(std::string c:{"B","C","D","E","F"}){// guaranteed sequential
		temp_opener=temp_header+ c +temp_footer;
		munuCMS.Add(temp_opener.c_str());
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
			case ttz:{           return mc__df;break;}
			case met:{           return bothdf;break;}
			case cms:{switch(ch){// MC is already false
			          case elnu:{return elnudf;break;}
			          case munu:{return munudf;break;}
//			          default  :throw std::invalid_argument(
//				"Unimplemented ch (rdf set)");
				}break;}
//			default :throw std::invalid_argument(
//				"Unimplemented ds (rdf set)");
		}
	}();
	switch(ch){
		case elnu:{temp_header = "Electron_";break;}
		case munu:{temp_header =     "Muon_";break;}
//		default  :throw std::invalid_argument(
//			"Unimplemented ch (init)");
	}
//	ROOT::EnableImplicitMT();
	// make test runs faster by restriction. Real run should not
	auto dfr = df.Range(10000);// remember to enable MT when NOT range
	auto init_selection = dfr// remove one letter to do all
	// lepton selection first
//	.Filter(met_pt_cut(ch),{"MET_pt"},"MET Pt cut")// TODO: Re-enable!
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
	.Define(   "lep__pt",temp_header+   "pt [loose_leps][0]")
	.Define(   "lep_eta",temp_header+   "eta[loose_leps][0]")
	.Define(   "lep_phi",temp_header+   "phi[loose_leps][0]")
	.Define(   "lep_mas",temp_header+  "mass[loose_leps][0]")
	.Define(   "lep___q",temp_header+"charge[loose_leps][0]")// only for rocco
	// now for transverse W; lepton selected
	.Define("tw_lep__pt",[](float x){return static_cast<double>(x);},{"lep__pt"})
	.Define("tw_lep_eta",[](float x){return static_cast<double>(x);},{"lep_eta"})
	.Define("tw_lep_phi",[](float x){return static_cast<double>(x);},{"lep_phi"})
	.Define("tw_lep_mas",transverse_w_mass,
	       {"tw_lep__pt",
	        "tw_lep_phi","MET_pt","MET_phi"})
	.Define(   "lep_nu_invmass",lep_nu_invmass,
	       {   "lep__pt",
	           "lep_eta",
	           "lep_phi",
	           "lep_mas",
	           "MET_pt" ,
	           "MET_phi"})
//	       "CaloMET_sumEt"})// TODO: add this back
//	.Filter(easy_mass_cut(W_MASS,W_MASS_CUT),{"tw_lep_mas"},"W mass cut")
	// jets selection follows; tW done and lepton selected
	.Define("jet_lep_min_dR" , jet_lep_min_deltaR<floats>,
	       {"Jet_eta","Jet_phi",  "lep_eta","lep_phi"})
	.Define("tight_jets"     , tight_jet_id   ,
	       {"jet_lep_min_dR" ,     "Jet_pt","Jet_eta","Jet_jetId"})
	.Filter( jetCutter(JETS_MIN,JETS_MAX),{"tight_jets"},"Jet cut")
	.Define("tight_jets__pt" ,     "Jet_pt  [tight_jets]")
	.Define("tight_jets_eta" ,     "Jet_eta [tight_jets]")
	.Define("tight_jets_phi" ,     "Jet_phi [tight_jets]")
	.Define("tight_jets_mas" ,     "Jet_mass[tight_jets]")
	.Define("tJ_btagCSVv2"   ,"Jet_btagCSVV2[tight_jets]")
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
	temp_footer = "pt vs eta in" + temp_footer + " channel for ";
	switch(ds){
		case tzq:{temp_header+="tzq";temp_footer+="tZq";break;}
		case  ww:{temp_header+="_ww";temp_footer+=" WW";break;}
		case  wz:{temp_header+="_wz";temp_footer+=" WZ";break;}
		case  zz:{temp_header+="_zz";temp_footer+=" ZZ";break;}
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
	.Define("cjer",    delta_R_jet_smear(false)  ,
	       {   "Jet_pt",   "Jet_eta",   "Jet_phi","fixedGridRhoFastjetAll",
	        "GenJet_pt","GenJet_eta","GenJet_phi"})
	.Define("fat_cjer",    delta_R_jet_smear(true)  ,
               {   "FatJet_pt",   "FatJet_eta",   "FatJet_phi","fixedGridRhoFastjetAll",
                "GenJetAK8_pt","GenJetAK8_eta","GenJetAK8_phi"})
	.Define("tight_cjer"  ,"cjer [tight_jets]")
	.Define("fin_jets__pt" ,"tight_jets__pt * tight_cjer")// these are JEC
	.Define("fin_jets_eta" ,"tight_jets_eta * tight_cjer")// Monte Carlo
	.Define("fin_jets_phi" ,"tight_jets_phi * tight_cjer")// needs JEC
	.Define("fin_jets_mas" ,"tight_jets_mas * tight_cjer")// fin = final
	// jets selected, now bjets and btagging preliminaries
	.Define("is_bjets"         ,is_bjet_id   ,{"fin_jets_eta","tJ_btagCSVv2"})
	.Filter(jetCutter(BJETS_MIN,BJETS_MAX)   ,{"is_bjets"},"b jet cut")
	.Define("no_bjets"         ,no_bjet_id   ,{"fin_jets_eta"   })
	.Define("tJpF"             ,"Jet_partonFlavour  [tight_jets]")
	.Define("is_btag_numer"    ,is_bjet_numer,{"tJpF","is_bjets"})
	.Define("no_btag_numer"    ,no_bjet_numer,{"tJpF","is_bjets"})
	.Define("is_btag_denom"    ,is_bjet_denom,{"tJpF","no_bjets"})
	.Define("no_btag_denom"    ,no_bjet_denom,{"tJpF","no_bjets"})
	.Define("is_btag_numer__pt","fin_jets__pt[is_btag_numer]")
	.Define("no_btag_numer__pt","fin_jets__pt[no_btag_numer]")
	.Define("is_btag_denom__pt","fin_jets__pt[is_btag_denom]")
	.Define("no_btag_denom__pt","fin_jets__pt[no_btag_denom]")
	.Define("is_btag_numer_eta",[](ints x,doubles y){return abs(y[x]);},
	       {"is_btag_numer"    ,"fin_jets_eta"})
	.Define("no_btag_numer_eta",[](ints x,doubles y){return abs(y[x]);},
	       {"no_btag_numer"    ,"fin_jets_eta"})
	.Define("is_btag_denom_eta",[](ints x,doubles y){return abs(y[x]);},
	       {"is_btag_denom"    ,"fin_jets_eta"})
	.Define("no_btag_denom_eta",[](ints x,doubles y){return abs(y[x]);},
	       {"no_btag_denom"    ,"fin_jets_eta"})
	.Define("sfi",btagCSVv2( true),//,btagDF),// checks btag
	       { "tJ_btagCSVv2","fin_jets__pt","fin_jets_eta", "tJpF"})
	.Define("sfj",btagCSVv2(false),//,btagDF),// ignore btag
	       { "tJ_btagCSVv2","fin_jets__pt","fin_jets_eta", "tJpF"})
	;
	auto reco = allReconstruction(
	     jecs_bjets )
	;

	auto h_is_btag_numer_PtVsEta
	   = reco.Histo2D({
	(        "is_numer_" + temp_header).c_str(),
	("MC is btag numer " + temp_footer).c_str(),
	12,0,200,12,0,3},
	"is_btag_numer__pt",
	"is_btag_numer_eta");
	h_is_btag_numer_PtVsEta->GetXaxis()->SetTitle("pT/GeV");
	h_is_btag_numer_PtVsEta->GetYaxis()->SetTitle("PseudoRapidity eta");

	auto h_no_btag_numer_PtVsEta
	   = reco.Histo2D({
	(        "no_numer_" + temp_header).c_str(),
	("MC no btag numer " + temp_footer).c_str(),
	12,0,200,12,0,3},
	"no_btag_numer__pt",
	"no_btag_numer_eta");
	h_no_btag_numer_PtVsEta->GetXaxis()->SetTitle("pT/GeV");
	h_no_btag_numer_PtVsEta->GetYaxis()->SetTitle("PseudoRapidity eta");

	auto h_is_btag_denom_PtVsEta
	   = reco.Histo2D({
	(        "is_denom_" + temp_header).c_str(),
	("MC is btag denom " + temp_footer).c_str(),
	12,0,200,12,0,3},
	"is_btag_denom__pt",
	"is_btag_denom_eta");
	h_is_btag_denom_PtVsEta->GetXaxis()->SetTitle("pT/GeV");
	h_is_btag_denom_PtVsEta->GetYaxis()->SetTitle("PseudoRapidity eta");

	auto h_no_btag_denom_PtVsEta
	   = reco.Histo2D({
	(        "no_denom_" + temp_header).c_str(),
	("MC no btag denom " + temp_footer).c_str(),
	12,0,200,12,0,3},
	"no_btag_denom__pt",
	"no_btag_denom_eta");
	h_no_btag_denom_PtVsEta->GetXaxis()->SetTitle("pT/GeV");
	h_no_btag_denom_PtVsEta->GetYaxis()->SetTitle("PseudoRapidity eta");

	TH2D *
	is_btag_ratio = static_cast<TH2D*>(h_is_btag_numer_PtVsEta->Clone());
	is_btag_ratio->Divide(             h_is_btag_denom_PtVsEta.GetPtr());
	is_btag_ratio->SetNameTitle("ei", "is b tag ei");
	is_btag_ratio->Draw("COLZ");// trigger getting everything done
	TH2D *
	no_btag_ratio = static_cast<TH2D*>(h_no_btag_numer_PtVsEta->Clone());
	no_btag_ratio->Divide(             h_no_btag_denom_PtVsEta.GetPtr());
	no_btag_ratio->SetNameTitle("ej", "no b tag ej");
	no_btag_ratio->Draw("COLZ");// DELETED @ EoF

	auto has_btag_eff
	   = reco
	.Define("IsEffBTagged",BTaggedEffGiver(is_btag_ratio, true),
	       {"fin_jets__pt","fin_jets_eta"})// TODO: check sensibility
	.Define("NoEffBTagged",BTaggedEffGiver(no_btag_ratio,false),
	       {"fin_jets__pt","fin_jets_eta"})// of this eff formula
	.Define("P___ei",    EffIsBTaggedProduct,{"IsEffBTagged"})
	.Define("P___ej",    EffNoBTaggedProduct,{"NoEffBTagged"})
	.Define("P_sfei",Sfi_EffIsBTaggedProduct,{"IsEffBTagged","sfi"})
	.Define("P_sfej",Sfj_EffNoBTaggedProduct,{"NoEffBTagged","sfj"})
	.Define("P_MC"  ,"P___ei * P___ej")
	.Define("P_Data","P_sfei * P_sfej")
	.Define("btag_w",btag_weight,{"P_Data","P_MC"})
	// next few lines are just for muons
	.Define("lep_gpt"    ,lep_gpt(ch),{"GenPart_pdgId","GenPart_pt"})
	.Define("lep__nl"    ,[=](ints L,ints m){if(munu==ch)return L[m][0];
	                                         else        return      0 ;},
	       {"Muon_nTrackerLayers","loose_leps"})
	.Define("mostSF"     ,lepEffGiver(ch,MC
	                       , recoLowEt,reco_pass,tight_94x
	                       , id_N,id_Y,id_A,id_T
	                       , isoN,isoY,isoA,isoT
	                      ),{"lep__pt","lep_eta",
	                         "lep_phi","lep___q",
	                         "lep_gpt","lep__nl"})
//	.Define("puSf",pile(PuWd,PuUd,PuDd),{"PV_npvs"})
	;
	auto finalDF = finalScaling(ds,
	     has_btag_eff )
	;
	// Copied to below, delete MC-only
	// Assuming temp_header and footer and all are set per (hist titles)!
	auto h_trans_w = finalDF.Histo1D({
	(          "tWm_"     + temp_header).c_str(),
	("Transverse W mass " + temp_header).c_str(),
	50,0,180},
	"tw_lep_mas","nw_tw_lep_mas");
	h_trans_w->GetXaxis()->SetTitle("mass GeV/C^2");
	h_trans_w->GetYaxis()->SetTitle("Event");
	h_trans_w->SetLineStyle(kSolid);

	auto h_Winvmas = finalDF.Histo1D({
	("W_invariant_mass_" + temp_header).c_str(),
	("W invariant mass " + temp_header).c_str(),
	50,0,180},
	"lep_nu_invmass","nw_lep_nu_invmass");
	h_Winvmas->GetXaxis()->SetTitle("mass GeV/C^2");
	h_Winvmas->GetYaxis()->SetTitle("Event");
	h_Winvmas->SetLineStyle(kSolid);

	auto h_z_mas = finalDF.Histo1D({
	(        "zmas_" + temp_header).c_str(),
	("Recon. Z mass" + temp_header).c_str(),
	50,0,200},
	"z_mas","nw_z_mas");

	auto h_trans_T = finalDF.Histo1D({
	(          "tTm_"     + temp_header).c_str(),
	("Transverse T mass " + temp_header).c_str(),
	50,0,180},
	"ttop_mas","nw_ttop_mas");
	h_trans_T->GetXaxis()->SetTitle("mass GeV/C^2");
	h_trans_T->GetYaxis()->SetTitle("Event");
	h_trans_T->SetLineStyle(kSolid);

	auto h_zmet_Dph = finalDF.Histo1D({
	("Z_MET_Delta_Phi_" + temp_header).c_str(),
	("Z MET Delta Phi " + temp_header).c_str(),
	50,-7,7},
	"zmet_Dph","nw_zmet_Dph");
	h_zmet_Dph->GetXaxis()->SetTitle("Z & MET delta phi/rad");
	h_zmet_Dph->GetYaxis()->SetTitle("Event");
	h_zmet_Dph->SetLineStyle(kSolid);

	auto h_zw_Dph = finalDF.Histo1D({
	("Z_W_Delta_Phi_" + temp_header).c_str(),
	("Z W Delta Phi " + temp_header).c_str(),
	50,-7,7},
	"zw_Dph","nw_zw_Dph");
	h_zw_Dph->GetXaxis()->SetTitle("Z & W delta phi/rad");
	h_zw_Dph->GetYaxis()->SetTitle("Event");
	h_zw_Dph->SetLineStyle(kSolid);

	auto h_z_daughters_Dph = finalDF.Histo1D({
	("Z_pair_jets_Delta_Phi_" + temp_header).c_str(),
	("Z pair jets Delta Phi " + temp_header).c_str(),
	50,-7,7},
	"z_jets_Dph","nw_z_jets_Dph");
	h_z_daughters_Dph->GetXaxis()->SetTitle("Z pair jets Delta phi/rad");
	h_z_daughters_Dph->GetYaxis()->SetTitle("Event");
	h_z_daughters_Dph->SetLineStyle(kSolid);

	auto h_tWmVsZmass = finalDF.Histo2D({
	("tWmVsZmass_" + temp_header).c_str(),
	("tWmVsZmass " + temp_header).c_str(),
	50,0,200,50,0,200},
	"tw_lep_mas","z_mas");
	h_tWmVsZmass->GetXaxis()->SetTitle("tWm   GeV/C^2");
	h_tWmVsZmass->GetYaxis()->SetTitle("Zmass GeV/C^2");

	auto h_sfi = finalDF.Histo1D({
	("sfi_"+temp_header).c_str(),
	("sfi "+temp_header).c_str(),
	50,-10,10},"sfi");

	auto h_sfj = finalDF.Histo1D({
	("sfj_"+temp_header).c_str(),
	("sfj "+temp_header).c_str(),
	50,-10,10},"sfj");

	auto h_p_ei = finalDF.Histo1D({
	("p_ei_"+temp_header).c_str(),
	("p_ei "+temp_header).c_str(),
	50,-10,10},"P___ei");

	auto h_p_ej = finalDF.Histo1D({
	("p_ej_"+temp_header).c_str(),
	("p_ej "+temp_header).c_str(),
	50,-10,10},"P___ej");

	auto h_p_sfei= finalDF.Histo1D({
	("p_sfei_"+temp_header).c_str(),
	("p_sfei "+temp_header).c_str(),
	50,-10,10},"P_sfei");

	auto h_p_sfej= finalDF.Histo1D({
	("p_sfej_"+temp_header).c_str(),
	("p_sfej "+temp_header).c_str(),
	50,-10,10},"P_sfej");

	auto h_btag_w= finalDF.Histo1D({
	("btag_w_"+temp_header).c_str(),
	("btag_W "+temp_header).c_str(),
	50,-100,100},"btag_w");

	auto h_ev_w = finalDF.Histo1D({
	(   "ev_w_"    +temp_header).c_str(),
	("Event weight"+temp_header).c_str(),
	50,-50,50},"sf");

	// write histograms to a root file
	// ASSUMES temp_header is correct!
	TFile hf(("histo/"+temp_header+".histo").c_str(),"RECREATE");
		hf.WriteTObject(h_sfi   .GetPtr());
		hf.WriteTObject(h_sfj   .GetPtr());
		hf.WriteTObject(h_p_ei  .GetPtr());
		hf.WriteTObject(h_p_ej  .GetPtr());
		hf.WriteTObject(h_p_sfei.GetPtr());
		hf.WriteTObject(h_p_sfej.GetPtr());
		hf.WriteTObject(h_btag_w.GetPtr());
		hf.WriteTObject(h_is_btag_numer_PtVsEta.GetPtr());
		hf.WriteTObject(h_no_btag_numer_PtVsEta.GetPtr());
		hf.WriteTObject(h_is_btag_denom_PtVsEta.GetPtr());
		hf.WriteTObject(h_no_btag_denom_PtVsEta.GetPtr());
		hf.WriteTObject(is_btag_ratio);
		hf.WriteTObject(no_btag_ratio);
	hf.WriteTObject(h_trans_T .GetPtr());
	hf.WriteTObject(h_trans_w .GetPtr());
	hf.WriteTObject(h_Winvmas .GetPtr());
	hf.WriteTObject(h_ev_w    .GetPtr());
	hf.WriteTObject(h_z_mas   .GetPtr());
	hf.WriteTObject(  h_zw_Dph.GetPtr());
	hf.WriteTObject(h_zmet_Dph.GetPtr());
	hf.WriteTObject(h_z_daughters_Dph.GetPtr());
	hf.WriteTObject(h_tWmVsZmass.GetPtr());
	// the following two for loops stack correctly
	for(std::string particle:{"fin_jets","lep","bjet"})
	for(PtEtaPhiM k:PtEtaPhiMall){
		std::string  kstring  = "_";
		std::string  xAxisStr;
		double xmin,xmax;
		switch(k){
			case pt :{kstring += "_pt";xmin =  0;xmax = 200;
			          xAxisStr = "pT/GeV"                  ;break;}
			case eta:{kstring += "eta";xmin = -3;xmax =  3 ;
			          xAxisStr = "PseudoRapidity eta"      ;break;}
			case phi:{kstring += "phi";xmin = -7;xmax =  7 ;
			          xAxisStr = "Azimuthal angle, phi/rad";break;}
			case  m :{kstring += "mas";xmin =  0;xmax = 200;
			          xAxisStr = "mass GeV/C^2"            ;break;}
//			default :throw std::invalid_argument(
//				"Unimplemented component (histo)");
		}
		temp_footer = particle + kstring;
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
		hf.WriteTObject(h.GetPtr());
	}
	hf.Flush();
	} else {
	auto expt_bjets
	   = init_selection
	.Filter(runLBfilter(runLBdict),{"run","luminosityBlock"},
	        "LuminosityBlock filter")
	.Define("fin_jets__pt",[](floats& x){return static_cast<doubles>(x);},
	     {"tight_jets__pt"})
	.Define("fin_jets_eta",[](floats& x){return static_cast<doubles>(x);},
	     {"tight_jets_eta"})
	.Define("fin_jets_phi",[](floats& x){return static_cast<doubles>(x);},
	     {"tight_jets_phi"})
	.Define("fin_jets_mas",[](floats& x){return static_cast<doubles>(x);},
	     {"tight_jets_mas"})
	.Define("is_bjets"         ,is_bjet_id   ,{"fin_jets_eta","tJ_btagCSVv2"})
	.Filter(jetCutter(BJETS_MIN,BJETS_MAX)   ,{"is_bjets"},"b jet cut")
	// TODO: Always check that the previous 2 lines are copies of earlier
	;
	auto reco = allReconstruction(
	     expt_bjets )
	;
	auto not_btag_eff
	   = reco
	.Define("btag_w"     ,[](){return 1.;})
	.Define("mostSF"     ,lepEffGiver(ch,MC
	                       , recoLowEt,reco_pass,tight_94x
	                       , id_N,id_Y,id_A,id_T
	                       , isoN,isoY,isoA,isoT
	                      ),{"lep__pt","lep_eta",
	                         "lep_phi","lep___q", // last 4 unused
	                         "lep_phi","lep___q"})// last 2 repeat is fine
//	.Define("puSf",[](){return std::map<puSf,double>
//	                   {{puW,1.},{upW,1.},{dnW,1.}};})
	;
	auto finalDF = finalScaling(ds,
	     not_btag_eff )
	;
	// Copied from earlier, delete MC-only
	// Assuming temp_header and footer and all are set per (hist titles)!
	auto h_trans_w = finalDF.Histo1D({
	(          "tWm_"     + temp_header).c_str(),
	("Transverse W mass " + temp_header).c_str(),
	50,0,180},
	"tw_lep_mas","nw_tw_lep_mas");
	h_trans_w->GetXaxis()->SetTitle("mass GeV/C^2");
	h_trans_w->GetYaxis()->SetTitle("Event");
	h_trans_w->SetLineStyle(kSolid);

	auto h_Winvmas = finalDF.Histo1D({
	("W_invariant_mass_" + temp_header).c_str(),
	("W invariant mass " + temp_header).c_str(),
	50,0,180},
	"lep_nu_invmass","nw_lep_nu_invmass");
	h_Winvmas->GetXaxis()->SetTitle("mass GeV/C^2");
	h_Winvmas->GetYaxis()->SetTitle("Event");
	h_Winvmas->SetLineStyle(kSolid);

	auto h_z_mas = finalDF.Histo1D({
	(        "zmas_" + temp_header).c_str(),
	("Recon. Z mass" + temp_header).c_str(),
	50,0,200},
	"z_mas","nw_z_mas");

	auto h_trans_T = finalDF.Histo1D({
	(          "tTm_"     + temp_header).c_str(),
	("Transverse T mass " + temp_header).c_str(),
	50,0,180},
	"ttop_mas","nw_ttop_mas");
	h_trans_T->GetXaxis()->SetTitle("mass GeV/C^2");
	h_trans_T->GetYaxis()->SetTitle("Event");
	h_trans_T->SetLineStyle(kSolid);

	auto h_zmet_Dph = finalDF.Histo1D({
	("Z_MET_Delta_Phi_" + temp_header).c_str(),
	("Z MET Delta Phi " + temp_header).c_str(),
	50,-7,7},
	"zmet_Dph","nw_zmet_Dph");
	h_zmet_Dph->GetXaxis()->SetTitle("Z & MET delta phi/rad");
	h_zmet_Dph->GetYaxis()->SetTitle("Event");
	h_zmet_Dph->SetLineStyle(kSolid);

	auto h_zw_Dph = finalDF.Histo1D({
	("Z_W_Delta_Phi_" + temp_header).c_str(),
	("Z W Delta Phi " + temp_header).c_str(),
	50,-7,7},
	"zw_Dph","nw_zw_Dph");
	h_zw_Dph->GetXaxis()->SetTitle("Z & W delta phi/rad");
	h_zw_Dph->GetYaxis()->SetTitle("Event");
	h_zw_Dph->SetLineStyle(kSolid);

	auto h_z_daughters_Dph = finalDF.Histo1D({
	("Z_pair_jets_Delta_Phi_" + temp_header).c_str(),
	("Z pair jets Delta Phi " + temp_header).c_str(),
	50,-7,7},
	"z_jets_Dph","nw_z_jets_Dph");
	h_z_daughters_Dph->GetXaxis()->SetTitle("Z pair jets Delta phi/rad");
	h_z_daughters_Dph->GetYaxis()->SetTitle("Event");
	h_z_daughters_Dph->SetLineStyle(kSolid);

	auto h_tWmVsZmass = finalDF.Histo2D({
	("tWmVsZmass_" + temp_header).c_str(),
	("tWmVsZmass " + temp_header).c_str(),
	50,0,200,50,0,200},
	"tw_lep_mas","z_mas");
	h_tWmVsZmass->GetXaxis()->SetTitle("tWm   GeV/C^2");
	h_tWmVsZmass->GetYaxis()->SetTitle("Zmass GeV/C^2");

	auto h_ev_w = finalDF.Histo1D({
	(   "ev_w_"    +temp_header).c_str(),
	("Event weight"+temp_header).c_str(),
	50,-50,50},"sf");

	// write histograms to a root file
	// ASSUMES temp_header is correct!
	TFile hf(("histo/"+temp_header+".histo").c_str(),"RECREATE");
/*	if(MC){
		hf.WriteTObject(h_sfi   .GetPtr());
		hf.WriteTObject(h_sfj   .GetPtr());
		hf.WriteTObject(h_p_ei  .GetPtr());
		hf.WriteTObject(h_p_ej  .GetPtr());
		hf.WriteTObject(h_p_sfei.GetPtr());
		hf.WriteTObject(h_p_sfej.GetPtr());
		hf.WriteTObject(h_btag_w.GetPtr());
		hf.WriteTObject(h_is_btag_numer_PtVsEta.GetPtr());
		hf.WriteTObject(h_no_btag_numer_PtVsEta.GetPtr());
		hf.WriteTObject(h_is_btag_denom_PtVsEta.GetPtr());
		hf.WriteTObject(h_no_btag_denom_PtVsEta.GetPtr());
		hf.WriteTObject(is_btag_ratio);
		hf.WriteTObject(no_btag_ratio);
	}
*/
	hf.WriteTObject(h_trans_T .GetPtr());
	hf.WriteTObject(h_trans_w .GetPtr());
	hf.WriteTObject(h_Winvmas .GetPtr());
	hf.WriteTObject(h_ev_w    .GetPtr());
	hf.WriteTObject(h_z_mas   .GetPtr());
	hf.WriteTObject(  h_zw_Dph.GetPtr());
	hf.WriteTObject(h_zmet_Dph.GetPtr());
	hf.WriteTObject(h_z_daughters_Dph.GetPtr());
	hf.WriteTObject(h_tWmVsZmass.GetPtr());
	// the following two for loops stack correctly
	for(std::string particle:{"fin_jets","lep","bjet"})
	for(PtEtaPhiM k:PtEtaPhiMall){
		std::string  kstring  = "_";
		std::string  xAxisStr;
		double xmin,xmax;
		switch(k){
			case pt :{kstring += "_pt";xmin =  0;xmax = 200;
			          xAxisStr = "pT/GeV"                  ;break;}
			case eta:{kstring += "eta";xmin = -3;xmax =  3 ;
			          xAxisStr = "PseudoRapidity eta"      ;break;}
			case phi:{kstring += "phi";xmin = -7;xmax =  7 ;
			          xAxisStr = "Azimuthal angle, phi/rad";break;}
			case  m :{kstring += "mas";xmin =  0;xmax = 200;
			          xAxisStr = "mass GeV/C^2"            ;break;}
//			default :throw std::invalid_argument(
//				"Unimplemented component (histo)");
		}
		temp_footer = particle + kstring;
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
		hf.WriteTObject(h.GetPtr());
	}
	hf.Flush();
	}
	std::cout<<"calchisto successfully completed"<<std::endl;
}
