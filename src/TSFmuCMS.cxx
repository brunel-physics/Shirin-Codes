//clang++ -Isrc -std=c++17 -march=native -pipe -Ofast -Wall -Wextra -Wpedantic -o build/tsfmuCMS src/TSFmuCMS.cxx `root-config --libs` -lm


#include <ROOT/RDataFrame.hxx>//#include <ROOT/RCsvDS.hxx>
#include <TRandom3.h>// used Gaussian, uniform each once
#include <TChain.h>

#include "csv.h"
#include "json.hpp"

using doubles = ROOT::VecOps::RVec<double>;
using  floats = ROOT::VecOps::RVec<float>;
using    ints = ROOT::VecOps::RVec<int>;
using   bools = ROOT::VecOps::RVec<bool>;
using strings = ROOT::VecOps::RVec<std::string>;

namespace{
enum dataSource  { muB,muC,muD,muE,muF };
enum channel     {munu};
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
/*
//constexpr  float    MET__PT_MIN = 40.f;
  constexpr  float    MET_EL_PT   = 20.f;//80.f;// TODO: Need new values
  constexpr  float    MET_MU_PT   = 25.f;//40.f;

  constexpr double     Z_MASS     =  91.1876;
  constexpr double     Z_MASS_CUT =  20.    ;*/
  constexpr double     W_MASS     =  80.385 ;
  constexpr double     W_MASS_CUT =  35.    ;/*
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
//constexpr double    ak8RconeBy2 =  .4;*/

  constexpr double          TZQ_W =  .0128;/*
  constexpr double       WWLNQQ_W = 2.1740;
  constexpr double       WZLNQQ_W =  .2335;
  constexpr double        TTBLV_W = 1.3791;
  constexpr double        TTZQQ_W =  .0237;
  constexpr double       ZZLLQQ_W =  .0485;
*/
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

inline auto Ptriggers(dataSource ds){
	return [=](
	 const bool el1// HLT_Ele35_WPTight_Gsf
	,const bool mu1// HLT_IsoMu27
	,const bool mu2// HLT_IsoMu24_eta2p1
	//const bool mu3// HLT_L1SingleMu25 keeping just incase
){
	switch(ds){
	case muB:
	case muC:
	case muD:return mu1 || mu2;
	case muE:
	case muF:return mu1;
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
	return [=](
		 const   ints& mask
		,const   ints& elids
		,const floats& isos
	){
// TODO: In the case we want 1 loose lepton, comment the 4 NOTE lines below
		bool result;
		if(false);
		else if(munu==ch){
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
inline double transverse_w_mass(
	 const double lep__pt
	,const double lep_phi
	,const float  met__pt
	,const float  met_phi
){
	return 2.
		* std::abs (std::sin( delta_phi(lep_phi
		          - static_cast<double>(met_phi))*.5))
		* std::sqrt(static_cast<double>(met__pt)
		                              * lep__pt);
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

inline auto sf(
	 const bool MC
	,const TH1D* const &PuWd
	,const TH1D* const &PuUd
	,const TH1D* const &PuDd
){
	if(0<debug) std::cout<<"scale factor "<<std::endl;
	return [=](const    int npv
	){
		// TODO: trigger efficiency
		double result  = 1;
		if(MC) result *= TZQ_W * pile(PuWd,PuUd,PuDd)(npv)[puW];
		return result;
	};
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
void TSFmuCMS ( const channel ch , const dataSource ds , const char b ){
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
	// pile up
	TFile *tF ;TH1D *t1d;
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
	tF = nullptr;t1d = nullptr;

	std::string temp_header,temp_opener,temp_footer="/*.root";/**/
	temp_header="/data/disk3/nanoAOD_2017/SingleMuon_NanoAOD25Oct2019_Run";
	switch(ds){
	case muB:{temp_opener = temp_header+"B"+temp_footer;break;}
        case muC:{temp_opener = temp_header+"C"+temp_footer;break;}
        case muD:{temp_opener = temp_header+"D"+temp_footer;break;}
        case muE:{temp_opener = temp_header+"E"+temp_footer;break;}
        case muF:{temp_opener = temp_header+"F"+temp_footer;break;}
	}
	ROOT::RDataFrame  munuCMS_df("Events",temp_opener);

	const bool MC = false;
	auto df = [&,ds](){// Get correct data frame
		switch(ds){
			case muB:// fall through :)
			case muC:
			case muD:
			case muE:
			case muF:{return munuCMS_df;break;}
			default :throw std::invalid_argument(
				"Unimplemented ds (rdf set)");
		}
	}();
	switch(ch){
		case munu:{temp_header =     "Muon_";break;}
		default  :throw std::invalid_argument(
			"Unimplemented ch (init)");
	}
	// make test runs faster by restriction. Real run should not
//	auto dfr = df.Range(10000);// remember to enable MT when NOT range
	auto origi = df// toggle one letter to do all
	.Define("lep_pts","static_cast<ROOT::RVec<double>>("+temp_header+"pt)")
	;
	auto clean = origi
        .Filter(runLBfilter(runLBdict),{"run","luminosityBlock"},
           "LuminosityBlock filter")
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
	;
	auto tight = clean
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
	.Define("lep__pt","static_cast<double>("+temp_header+  "pt[loose_leps][0])")
	.Define("lep_eta","static_cast<double>("+temp_header+ "eta[loose_leps][0])")
	.Define("lep_phi","static_cast<double>("+temp_header+ "phi[loose_leps][0])")
	.Define("tw_lep_mas",transverse_w_mass,
	       {"lep__pt",
	        "lep_phi","MET_pt","MET_phi"})
	// jets selection follows; lepton selected
	.Define("rawJet_eta","static_cast<ROOT::RVec<double>>(Jet_eta )")
	.Define("rawJet_phi","static_cast<ROOT::RVec<double>>(Jet_phi )")
	.Define("jet_lep_min_dR"   ,jet_lep_min_deltaR,// later reused with doubles
	       {"rawJet_eta","rawJet_phi","lep_eta","lep_phi"})// gcc fail template
	.Define("sf",sf(MC,PuWd,PuUd,PuDd),{"PV_npvs"})
	;
	auto Ptrig = tight
	.Filter(Ptriggers(ds),
		{ "HLT_Ele35_WPTight_Gsf"
		 ,"HLT_IsoMu27"
		 ,"HLT_IsoMu24_eta2p1"
		},"Tag Triggers Filter")
	;

	std::string opener(1,b); opener.reserve(10);
	switch(ch){
	case munu:{opener += "_munu_";break;}
	}
	switch(ds){
	case  muB:{opener += "muB"   ;break;}
	case  muC:{opener += "muC"   ;break;}
	case  muD:{opener += "muD"   ;break;}
	case  muE:{opener += "muE"   ;break;}
	case  muF:{opener += "muF"   ;break;}
	}
	TFile tsf(("histo/tsf"+opener+".root").c_str(),"RECREATE");
	auto origiPt = origi.Histo1D({
	    "origiPt"       ,
	    "Original P_{T}",
	50,0,400},"lep_pts");
	auto cleanPt = clean.Histo1D({
	    "cleanPt"       ,
	    "Cleaned P_{T}" ,
	50,0,400},"lep_pts");
	auto tightPt = tight.Histo1D({
	    "tightPt"       ,
	    "Tight P_{T} "  ,
	50,0,400},"lep__pt" ,"sf");
	tsf.WriteTObject(origiPt.GetPtr());tsf.Flush();sync();
	tsf.WriteTObject(cleanPt.GetPtr());tsf.Flush();sync();
	tsf.WriteTObject(tightPt.GetPtr());tsf.Flush();sync();
	switch(b)
	{case 'l':{
	Ptrig.Report() ->Print();
	auto  trigPt = Ptrig.Histo1D({
	    "PtrigPt"           ,
	"Probe trigger P_{T} " ,
	50,0,400},"lep__pt"     ,"sf");
	tsf.WriteTObject( trigPt.GetPtr());tsf.Flush();sync();
	break;
	}case 'x':{
	Ptrig.Report() ->Print();
	auto  trigPt = Ptrig.Histo1D({
	    "TtrigPt"           ,
	" Tag trigger P_{T} " ,
	50,0,400},"lep__pt"     ,"sf");
	tsf.WriteTObject( trigPt.GetPtr());tsf.Flush();sync();
	break;
	}default:throw std::logic_error("Only L or C triggers please");}
	std::cout << "TriggerSF completed successfully" << std::endl;
}
int main ( int argc , char *argv[] ){
	if ( argc < 2 ) {
		std::cout << "Error: no command provided" << std::endl ;
		return 1 ;
	}
	if ( argc < 3 ) {
		   std::cout
		<< "Error: tsf needs channel and data source+L or C (for cross trigger)"
		<< std::endl
		<< "e.g.   tsf elnu tzqL"
		<< std::endl
		;
		return 2 ;
	}
	channel c ; dataSource d ; char b ;
	     if ( const auto chN = std::string_view( argv[1] ) ; false ) ;
	else if ( "munu"  == chN ) c = munu ;
	else { std::cout << "Error: channel " << chN
		<< " not recognised" << std::endl ;
		return 3 ;
	}
	     if ( const auto dsN = std::string_view( argv[2] ) ; false ) ;
	else if ( "muBC" ==  dsN ){ d = muB ; b = 'x' ; }
	else if ( "muCC" ==  dsN ){ d = muC ; b = 'x' ; }
	else if ( "muDC" ==  dsN ){ d = muD ; b = 'x' ; }
	else if ( "muEC" ==  dsN ){ d = muE ; b = 'x' ; }
	else if ( "muFC" ==  dsN ){ d = muF ; b = 'x' ; }
	else { std::cout << "Error: data source " << dsN
		<< " not recognised" << std::endl ;
		return 4 ;
	}
		TSFmuCMS(c,d,b) ;
		return 0 ;
}
