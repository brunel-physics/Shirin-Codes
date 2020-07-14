// TODO:: APPLY LEP JET min Dr as a filter as well

//clang++ -Isrc -std=c++17 -march=native -pipe -Ofast -Wall -Wextra -Wpedantic -o build/tsf src/TriggerSF.cxx `root-config --libs` -lm
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
enum dataSource  { ttb,cms };
enum channel     {elnu,munu};
  constexpr    int debug = 0;
//constexpr    int EL_MAX_NUM     = 1      ;
  constexpr  float EL__PT_MIN     = 15.f   ;// TODO: plot -> pick
  constexpr  float EL_LPT_MIN     = 15.f   ;// TODO: plot -> pick
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
  constexpr  float   MU__PT_MIN   = 15.f;// TODO: plot -> pick
  constexpr  float   MU_LPT_MIN   = 15.f;// TODO: plot -> pick
  constexpr  float   MU_ETA_MAX   = 2.4f;
  constexpr  float   MU_LOOSE_ISO = .25f;
  constexpr  float   MU_TIGHT_ISO = .15f;
/*
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
  constexpr double       WZLNQQ_W =  .2335;*/
  constexpr double        TTBLV_W = 1.3791;
/*constexpr double        TTZQQ_W =  .0237;
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

inline auto ltriggers(channel ch){
	return [=](
	 const bool el// HLT_Ele32_WPTight_Gsf_L1DoubleEG
	,const bool mu// HLT_IsoMu27
){
	switch(ch){
	case elnu:return el;
	case munu:return mu;
	}};
}
inline auto xtriggers(channel ch){
	return [=](
	 const bool en// HLT_Ele28_eta2p1_WPTight_Gsf_HT150
	,const bool jt// HLT_PFHT180
	,const bool mu// HLT_IsoMu27
){
	switch(ch){
	case elnu:return en;// el && jt
	case munu:return mu && jt;
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
	return [=](
		 const floats& pt
		,const   ints& mask
		,const   ints& elids
		,const floats& isos
	){
// TODO: In the case we want 1 loose lepton, comment the 4 NOTE lines below
		bool result;
		 if(false) ;
		 else if(elnu==ch){
			ints   temp = elids[mask];
			result = temp.size() == 1;// Choosing 1 Tight Lepton
			result = result && pt[mask][0] >= EL_LPT_MIN  ;// NOTE
			result = result &&    temp [0] >= EL_TIGHT_ID ;// NOTE
		}else if(munu==ch){
			floats temp = isos[mask];
			result = temp.size() == 1;
			result = result && pt[mask][0] >= MU_LPT_MIN  ;// NOTE
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
		if(MC) result *= TTBLV_W * pile(PuWd,PuUd,PuDd)(npv)[puW];
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
void TriggerSF ( const channel ch , const dataSource ds , const char b ){
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

	std::string temp_header="/data/disk0/nanoAOD_2017/",
	temp_opener,temp_footer="/*.root";/**/
	switch(ds){// CMS and MET MUST do some OPENABLE file ; reject later
	case ttb:{temp_opener=temp_header+"TTToSemileptonic"+temp_footer;break;}
	case cms:{temp_opener=temp_header+"TTToSemileptonic"+temp_footer;break;}
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
	const bool MC = ttb == ds;
	auto df = [&,ch,ds](){// Get correct data frame
		switch(ds){
			case ttb:{           return mc__df;break;}
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
		default  :throw std::invalid_argument(
			"Unimplemented ch (init)");
	}
	// make test runs faster by restriction. Real run should not
//	auto dfr = df.Range(10000);// remember to enable MT when NOT range
	auto origi = df// toggle one letter to do all
	.Define("lep_pts","static_cast<ROOT::RVec<double>>("+temp_header+"pt)")
	;
	auto clean = origi
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
	.Filter(lep_tight_cut(ch),{temp_header+"pt","loose_leps",
	        "Electron_cutBased",// edit function for  tight -> loose
	        "Muon_pfRelIso04_all"},"lepton cut")// left with 1 tight lepton
	.Define("lep__pt","static_cast<double>("+temp_header+  "pt[loose_leps][0])")
	.Define("lep_eta","static_cast<double>("+temp_header+ "eta[loose_leps][0])")
	.Define("lep_phi","static_cast<double>("+temp_header+ "phi[loose_leps][0])")
	// jets selection follows; lepton selected
	.Define("rawJet_eta","static_cast<ROOT::RVec<double>>(Jet_eta )")
	.Define("rawJet_phi","static_cast<ROOT::RVec<double>>(Jet_phi )")
	.Define("jet_lep_min_dR"   ,jet_lep_min_deltaR,// later reused with doubles
	       {"rawJet_eta","rawJet_phi","lep_eta","lep_phi"})// gcc fail template
	.Define("sf",sf(MC,PuWd,PuUd,PuDd),{"PV_npvs"})
	;
	auto ltrig = tight
	.Filter(ltriggers(ch),
		{ "HLT_Ele32_WPTight_Gsf_L1DoubleEG"//"HLT_Ele28_eta2p1_WPTight_Gsf_HT150"
		 ,"HLT_IsoMu27"
		},"Lepton Triggers Filter")
	;
	auto xtrig = tight
	.Filter(xtriggers(ch),
		{ "HLT_Ele28_eta2p1_WPTight_Gsf_HT150"//"HLT_PFMET120_PFMHT120_IDTight"
		 ,"HLT_PFHT180"
		 ,"HLT_IsoMu27"
		},"Cross Triggers Filter")
	;
	if(MC){
	switch(b)
	{case 'l':{
	auto ttb_df = ltrig.Report();
	ttb_df->Print();
	}case 'x':{
	auto ttb_df = xtrig.Report();
	ttb_df->Print();
	}default:throw std::logic_error("Only L or C triggers please");}
	}else{
	switch(b)
	{case 'l':{
	auto CMS_df = ltrig
	.Filter(runLBfilter(runLBdict),{"run","luminosityBlock"},
	        "LuminosityBlock filter")
	.Report()
	;
	CMS_df->Print();
	}case 'x':{
	auto CMS_df = xtrig
	.Filter(runLBfilter(runLBdict),{"run","luminosityBlock"},
	        "LuminosityBlock filter")
	.Report()
	;
	CMS_df->Print();
	}default:throw std::logic_error("Only L or C triggers please");}
	}
	std::string opener(1,b); opener.reserve(10);
	switch(ch){
	case elnu:{opener += "_elnu_";break;}
	case munu:{opener += "_munu_";break;}
	}
	switch(ds){
	case  ttb:{opener += "ttb";break;}
	case  cms:{opener += "cms";break;}
	}
	TFile tsf(("histo/tsf"+opener+".root").c_str(),"RECREATE");
	auto origiPt =    df.Histo1D({
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
	auto triggPt = ltrig.Histo1D({
	   "ltrigPt"            ,
	"Lepton trigger P_{T} " ,
	50,0,400},"lep__pt"     ,"sf");
	tsf.WriteTObject(triggPt.GetPtr());tsf.Flush();sync();
	}case 'x':{
	auto triggPt = xtrig.Histo1D({
	   "xtrigPt"            ,
	 "Cross trigger P_{T} " ,
	50,0,400},"lep__pt"     ,"sf");
	tsf.WriteTObject(triggPt.GetPtr());tsf.Flush();sync();
	}default:throw std::logic_error("Only L or C triggers please");}
	tsf.Close();
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
		<< "e.g.   tsf elnu ttbL"
		<< std::endl
		;
		return 2 ;
	}
	channel c ; dataSource d ; char b ;
	     if ( const auto chN = std::string_view( argv[1] ) ; false ) ;
	else if ( "elnu"  == chN ) c = elnu ;
	else if ( "munu"  == chN ) c = munu ;
	else { std::cout << "Error: channel " << chN
		<< " not recognised" << std::endl ;
		return 3 ;
	}
	     if ( const auto dsN = std::string_view( argv[2] ) ; false ) ;
	else if ( "ttbL" ==  dsN ){ d = ttb ; b = 'l' ; }
	else if ( "ttbC" ==  dsN ){ d = ttb ; b = 'x' ; }
	else if ( "cmsL" ==  dsN ){ d = cms ; b = 'l' ; }
	else if ( "cmsC" ==  dsN ){ d = cms ; b = 'x' ; }
	else { std::cout << "Error: data source " << dsN
		<< " not recognised" << std::endl ;
		return 4 ;
	}
		TriggerSF(c,d,b) ;
		return 0 ;
}
