//clang++ -Isrc -std=c++17 -march=native -pipe -O3 -Wall -Wextra -Wpedantic -o build/tsf src/TriggerSF.cxx `root-config --libs` -lm

// TODO :: 1. How to make sure 1 lepton is tight and 1 lepton is loose for the tight_lep_cut
// TODO :: 2.Compute the Ns from the reports
// TODO :: 3.divide data dist. / MC dist.


// Tag and Probe Methodology:
/*
The idea of the tag and probe is to use Z->ll enriched samples
(for MC, while for data I guess you get what you get) and try
to reconstruct the Z mass through the double lepton invariant mass.
One lepton is the tag, the other one is the probe. So you will
have events where the reconstruction is successful (what I called
“pass” before) and events where the probe didn’t pass your selection
criteria, so this will go into a distribution of “fail”. Eventually,
your distribution of “fail” will contain both background events,
that is events where the probe wasn’t actually coming from a Z,
and events where the probe was from a Z but just didn’t pass the selection.
With a proper fit you can distinguish the former from the latter and obtain
the number of good Z events which ended up in the fail distribution.
Then the efficiency is simply:

N(pass)/N(total), where N(total) = N(pass)+N(fail which were actually good Z events)

You do this for data and for MC, you divide the efficiencies obtained and you get the SFs.
*/

#include <ROOT/RDataFrame.hxx>//#include <ROOT/RCsvDS.hxx>
#include <TRandom3.h>// used Gaussian, uniform each once
#include <TChain.h>
#include <TF1.h>
#include <Math/Vector4D.h>

#include "csv.h"
#include "json.hpp"

using doubles = ROOT::VecOps::RVec<double>;
using  floats = ROOT::VecOps::RVec<float>;
using    ints = ROOT::VecOps::RVec<int>;
using   bools = ROOT::VecOps::RVec<bool>;
using strings = ROOT::VecOps::RVec<std::string>;

namespace{
enum dataSource  { dy,cms };// DY : DYJetsToLL
enum channel     {elnu,munu};

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
  constexpr  float    MET_MU_PT   = 25.f;//40.f;*/

  constexpr double     Z_MASS     =  91.1876;
  constexpr double     Z_MASS_CUT =  40.    ;/*
  constexpr double     W_MASS     =  80.385 ;
  constexpr double     W_MASS_CUT =  35.    ;
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
*/
  constexpr double      DYTOLL_W = 22.8750;
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
	 const bool el// HLT_Ele32_WPTight_Gsf_L1DoubleEG // TBC to: HLT_Ele32_WPTight_Gsf
	,const bool mu// HLT_IsoMu27
	//const bool mu3// HLT_L1SingleMu25 keeping just incase
){
	switch(ch){
	case elnu:return el;
	case munu:return mu;
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
		 const   ints& mask
		,const   ints& elids
		,const floats& isos
		,const   ints& chg // electric charge
	){
// TODO: In the case we want 1 loose lepton, comment the 4 NOTE lines below
		bool result;
		// For MC (DYJetsToLL) we need only one tight lepton
                // to use the Single lepton trigger on, however,
                // for data (cms single lepton) we choose a tight
                // lepton and apply a veto lepton (same channel)
		// i.e. single electron has , 1 tight e and 1 veto e
                ints   chmk   = chg[mask];
		if(false) ;
		 else if(elnu==ch){
			ints   temp = elids[mask];
			result = temp.size() == 2;// Choosing 2 Leptons
			result = (result && temp [0] >= EL_TIGHT_ID // NOTE
                               		 && temp [1] >= EL_LOOSE_ID // NOTE
					 && temp [1] <  EL_TIGHT_ID// NOTE
					 && chmk [0] != chmk [1])
			      ||
				 (result && temp [1] >= EL_TIGHT_ID // NOTE
                                         && temp [0] >= EL_LOOSE_ID // NOTE
                                         && temp [0] <  EL_TIGHT_ID// NOTE
                                         && chmk [1] != chmk [0])

					;
		}else if(munu==ch){
			floats temp = isos[mask];
			result = temp.size() == 2;
			result = (result && temp [0] <= MU_TIGHT_ISO // NOTE
			       	 	 && temp [1] <= MU_LOOSE_ISO // NOTE
					 && temp [1] >  MU_TIGHT_ISO // NOTE
					 && chmk [0] != chmk [1])
			      ||
				 (result && temp [1] <= MU_TIGHT_ISO // NOTE
                                         && temp [0] <= MU_LOOSE_ISO // NOTE
                                         && temp [0] >  MU_TIGHT_ISO // NOTE
                                         && chmk [1] != chmk [0])
					;
		}else{throw std::invalid_argument(
			"Unimplemented ch (lep_tight_cut)");}
		return result;
	};
}
inline auto lep_num(const doubles pt){
	// Accepting only 2 or more leptons
	return pt.size() >= 2;
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
auto find_z_pair(
	 const doubles& pts
	,const doubles& etas
	,const doubles& phis
	,const doubles& ms
	,const    ints& ch //charge
){
	if(0<debug) std::cout<<"find z pair"<<std::endl;
	// This function finds the pair nearest to z mass
	double  z_reco_mass = std::numeric_limits<double>::infinity();
	size_t  lep_index_1 = std::numeric_limits<size_t>::max();
	size_t  lep_index_2 = std::numeric_limits<size_t>::max();
	const size_t  nleps = pts.size();
	if(!all_equal(nleps,etas.size(),phis.size(),ms.size(),ch.size()))
		throw std::logic_error(
		"Collections must be the same size in Z-pair");
	if(pts.size()==0)	throw std::logic_error(
		"Collections must not be empty in Z-pair");
	ints z_pair(nleps, 0);// vector of zeroes
	for(size_t   i=0; i < nleps-1 ;++i)
	for(size_t j=i+1; j < nleps   ;++j)
	if(ch[i] != ch[j]){// Must be opposite sign lepton
		ROOT::Math::PtEtaPhiMVector
		lep1(pts[i],etas[i],phis[i],ms[i]),
		lep2(pts[j],etas[j],phis[j],ms[j]);
		if (const double reco_mass = (lep1+lep2).M();
		std::abs(Z_MASS- reco_mass) < std::abs(Z_MASS-z_reco_mass)){
		   z_reco_mass = reco_mass;// found nearer pair to z mass
		   lep_index_1 = i;
		   lep_index_2 = j;
		}
	}
	z_pair[lep_index_1] = 1;
	z_pair[lep_index_2] = 1;
	//if(1<debug) std::cout<<"z pair"<<z_pair<<std::endl;
	return z_pair;
}
//template <typename T>
inline auto z_num(const ints z_pair){
	std::cout<<"in z_num"<<std::endl;
	bool z_size = false;
	if(z_pair.size() >= 2) z_size = true;
	return z_size; //should never be zero
}
inline auto pt_pair(const doubles pt_pair){
	return pt_pair.size() >=2;
}
inline auto easy_mass_cut(const double theo,const double cut){
               return [=](const double ours){return std::abs(ours-theo)<cut;};
}
auto LVpairAdd(
	 const doubles& pt__pair// Create LorentzV from 2 jets 4-mom
	,const doubles& eta_pair
	,const doubles& phi_pair
	,const doubles& mas_pair
){
	std::cout<< pt__pair.size()<<std::endl;
	//if(0<debug) std::cout<<"LVpairAdd"<<std::endl;
	if(2 != pt__pair.size())throw std::logic_error(
		"Not pair of Z (LVpairAdd)");
	ROOT::Math::PtEtaPhiMVector
	v(pt__pair[0],eta_pair[0],phi_pair[0],mas_pair[0]),
	p(pt__pair[1],eta_pair[1],phi_pair[1],mas_pair[1]);
	return v+p;
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
		if(MC) result *= DYTOLL_W * pile(PuWd,PuUd,PuDd)(npv)[puW];
		return result;
	};
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

}// namespace
void TriggerSF ( const channel ch , const dataSource ds){
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

	std::string temp_header="/data/disk3/nanoAOD_2017/",
	temp_opener,temp_footer="/*.root";/**/
	switch(ds){// CMS and MET MUST do some OPENABLE file ; reject later
	case dy :{temp_opener=temp_header+"DYJetsToLL_M-50"+temp_footer;break;}
	case cms:{temp_opener=temp_header+"DYJetsToLL_M-50"+temp_footer;break;}
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
	const bool MC = dy == ds;
	auto df = [&,ch,ds](){// Get correct data frame
		switch(ds){
			case dy :{           return mc__df;break;}
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
        .Filter(triggers(ch),
                { "HLT_Ele32_WPTight_Gsf_L1DoubleEG"
                 ,"HLT_IsoMu27"
                },"Triggers Filter")
	;
	auto lumclean = clean
	.Filter(runLBfilter(runLBdict,MC),{"run","luminosityBlock"},
			   "LuminosityBlock filter")
	;
	auto offlep = lumclean
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
        .Define("lep__pt","static_cast<ROOT::RVec<double>>("+temp_header+     "pt[loose_leps])")
        .Define("lep_eta","static_cast<ROOT::RVec<double>>("+temp_header+    "eta[loose_leps])")
        .Define("lep_phi","static_cast<ROOT::RVec<double>>("+temp_header+    "phi[loose_leps])")
	.Define("lep_mas","static_cast<ROOT::RVec<double>>("+temp_header+   "mass[loose_leps])")
	.Define("lep_chg","static_cast<ROOT::RVec<int>>   ("+temp_header+ "charge[loose_leps])")
	.Filter(lep_num,{"lep__pt"},"lepton number cut, accepting 2 or more leptons")
	.Define("z_reco_leps"       , find_z_pair,
	       { "lep__pt",
	         "lep_eta",
	         "lep_phi",//easier to push back
	         "lep_mas",
		 "lep_chg"})
        .Define(   "z_pair__pt"   ,   "lep__pt[z_reco_leps]")
        .Define(   "z_pair_eta"   ,   "lep_eta[z_reco_leps]")
        .Define(   "z_pair_phi"   ,   "lep_phi[z_reco_leps]")
        .Define(   "z_pair_mas"   ,   "lep_mas[z_reco_leps]")
	.Filter(z_num  ,{"z_reco_leps"},"z pairs should exist")
	.Filter(pt_pair,{"z_pair__pt" },"z pairs should exist")
	;
	//auto z_lep = offlep
	auto probe = offlep
	/*.Define(   "z_pair__pt"   ,   "lep__pt[z_reco_leps]")
	.Define(   "z_pair_eta"   ,   "lep_eta[z_reco_leps]")
	.Define(   "z_pair_phi"   ,   "lep_phi[z_reco_leps]")
	.Define(   "z_pair_mas"   ,   "lep_mas[z_reco_leps]")*/
	.Define(   "z_LV"         , LVpairAdd    ,
	       {   "z_pair__pt"   ,
	           "z_pair_eta"   ,
	           "z_pair_phi"   ,
	           "z_pair_mas"  })
	.Define(   "z_mas"        ,"z_LV. M ()")
	.Define(   "z__pt"        ,"z_LV.Pt ()")
	.Define(   "z_eta"        ,"z_LV.Eta()")
	.Filter( easy_mass_cut(Z_MASS,Z_MASS_CUT),{"z_mas"},"z mass cut")
	// jets selection follows; lepton selected
	.Define("sf",sf(MC,PuWd,PuUd,PuDd),{"PV_npvs"})
	;
	/*auto probe = z_lep
	.Filter(triggers(ch),
		{ "HLT_Ele32_WPTight_Gsf_L1DoubleEG"
		 ,"HLT_IsoMu27"
		},"Triggers Filter")
	;*/
	auto  tag = probe
        .Filter(lep_tight_cut(ch),{"loose_leps",
        "Electron_cutBased",// edit function for  tight -> loose
        "Muon_pfRelIso04_all",
         temp_header+"charge"},"lepton cut")// left with 1 tight lepton
	;
	std::string opener;
	switch(ch){
	case elnu:{opener += "_elnu_"       ;break;}
	case munu:{opener += "_munu_"       ;break;}
	}
	switch(ds){
	case  dy :{opener += "DYJetsToLL"   ;break;}
	case  cms:{opener += "cms"          ;break;}
	}
	TFile tsf(("histo/tsf"+opener+".root").c_str(),"RECREATE");
	auto ProbeZmass = probe.Histo1D({
	    "probeZmass"     ,
	    "probe Z\\text{mass GeV/}c^{2}" ,
	50,60,150},"z_mas","sf");
	TF1 *f1 = new TF1("f1","gaus",60,150);
	f1->SetParameters(ProbeZmass->GetMaximum()
			, ProbeZmass->GetMean()
			, ProbeZmass->GetRMS() );
	ProbeZmass->Fit("f1");

	auto tagZmass  =  tag.Histo1D({
	    "tagZmass"       ,
	    "tag Z\\text{mass GeV/}c^{2}"   ,
	50,60,150},"z_mas","sf");
        TF1 *f2 = new TF1("f2","gaus",60,150);
        f2->SetParameters(tagZmass->GetMaximum()
                        , tagZmass->GetMean()
                        , tagZmass->GetRMS() );
        tagZmass->Fit("f2");

	tsf.WriteTObject(ProbeZmass.GetPtr());tsf.Flush();sync();
	tsf.WriteTObject(  tagZmass.GetPtr());tsf.Flush();sync();
	tag.Report() ->Print();
	std::cout << "TriggerSF completed successfully" << std::endl;
}
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
	else if ( "dy"  ==  dsN ){ d = dy  ;}
	else if ( "cms" ==  dsN ){ d = cms ;}
	else { std::cout << "Error: data source " << dsN
		<< " not recognised" << std::endl ;
		return 4 ;
	}
		TriggerSF(c,d) ;
		return 0 ;
}
