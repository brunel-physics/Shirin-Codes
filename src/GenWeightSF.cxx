//clang++ -Isrc -std=c++17 -march=native -pipe -O3 -Wall -Wextra -Wpedantic -o build/genW src/GenWeightSF.cxx `root-config --libs` -lm

#include <ROOT/RDataFrame.hxx>//#include <ROOT/RCsvDS.hxx>
#include <TRandom3.h>// used Gaussian, uniform each once
#include <TChain.h>
#include <TF1.h>
#include <Math/Vector4D.h>
#include "calchisto.hpp"
#include "csv.h"
#include "json.hpp"

using doubles = ROOT::VecOps::RVec<double>;
using  floats = ROOT::VecOps::RVec<float>;
using    ints = ROOT::VecOps::RVec<int>;
using   bools = ROOT::VecOps::RVec<bool>;
using strings = ROOT::VecOps::RVec<std::string>;

void GenWeightSF(const channel ch,const dataSource ds){

	// Open data files even if unused
	// then automatically choose which one to read from
	// No penalty for opening and leaving unused
	// Can even open multiple times at once in parallel
	// Open MC data source EVEN IF UNUSED
	std::string temp_header="/data/disk0/nanoAOD_2017/",
	temp_opener,temp_footer="/*.root";/**/
	switch(ds){// CMS and MET MUST do some OPENABLE file ; reject later
	case  tzq:{temp_opener="/data/disk3/nanoAOD_2017/tZqlvqq/*.root"  ;break;}/**/
	case   ww:{temp_opener=temp_header+ "WWToLNuQQ"       +temp_footer;break;}
	case   wz:{temp_opener=temp_header+ "WZTo1L1Nu2Q"     +temp_footer;break;}
	case   zz:{temp_opener=temp_header+ "ZZTo2L2Q"        +temp_footer;break;}
	case  wjt:{temp_opener="/data/disk3/nanoAOD_2017/WPlusJets_NanoAODv5/*.root";break;}/**/
        case  wjx:{temp_opener="/nfs/data/eepgssg/WJetsToLNu_ext_NanoAODv5"+temp_footer;break;}
	case  ttb:{temp_opener=temp_header+"TTToSemileptonic" +temp_footer;break;}
	case  tz1:{temp_opener=temp_header+"ttZToQQ"          +temp_footer;break;}
	case  tz2:{temp_opener=temp_header+"ttZToQQ_ext"      +temp_footer;break;}
        case  ttl:{temp_opener="/data/disk1/nanoAOD_2017_new/TT_2l2nu_nanoAODv5"+temp_footer;break;}
        case  ttj:{temp_opener=temp_header+"TTToHadronic"     +temp_footer;break;}
        case   st:{temp_opener="/data/disk1/nanoAOD_2017_new/ST_tchannel_top_nanoAODv5"+temp_footer;break;}
        case  stb:{temp_opener="/data/disk1/nanoAOD_2017_new/ST_tchannel_antitop_nanoAODv5"+temp_footer;break;}
        case  stw:{temp_opener=temp_header+"ST_tW"            +temp_footer;break;}
        case stbw:{temp_opener=temp_header+"ST_tbarW"         +temp_footer;break;}
	case wzll:{temp_opener=temp_header+"WZTo2L2Q"         +temp_footer;break;}
	case wjqq:{temp_opener=temp_header+"WPlusJetsToQQ"    +temp_footer;break;}
	case zjt1:{temp_opener="/data/disk3/nanoAOD_2017/ZPlusJets_M10To50_NanoAODv5"+temp_footer;break;}
	case zjt2:{temp_opener="/nfs/data/eepgssg/ZPlusJets_M50_NanoAODv5"+temp_footer;break;}
	case zjt3:{temp_opener="/nfs/data/eepgssg/ZPlusJets_M50_ext_NanoAODv5"+temp_footer;break;}
	case zjqq:{temp_opener="/data/disk3/nanoAOD_2017/DYToQQ"+temp_footer;break;}
	case  met:{temp_opener=temp_header+"ttZToQQ_ext"      +temp_footer;break;}
	case  cms:{temp_opener=temp_header+"ttZToQQ"          +temp_footer;break;}

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
                        case  tzq:
                        case   ww:// fall through!
                        case   wz:
                        case   zz:
                        case  ttb:
			case  ttj:
			case  ttl:
                        case  wjt:
			case  wjx:
                        case  tz1:
                        case  tz2:
			case   st:
			case  stb:
			case  stw:
			case stbw:
			case wzll:
			case zjt1:
			case zjt2:
			case zjt3:
			case zjqq:
			case wjqq:{           return mc__df;break;}
                        case  met:{           return bothdf;break;}
                        case  cms:{switch(ch){// MC is already false
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
//              default  :throw std::invalid_argument(
//                      "Unimplemented ch (init)");
        }
	// make test runs faster by restriction. Real run should not
        auto dfr = df.Range(100000);// remember to enable MT when NOT range
        auto genWpos = df// remove one letter to do all
	.Filter("genWeight >= 0", "genPosCount")
	;
	genWpos.Report() ->Print();
}
int main ( int argc , char *argv[] ){
        if ( argc < 2 ) {
                std::cout << "Error: no command provided" << std::endl ;
                return 1 ;
        }
	if ( argc < 3 ) {
                   std::cout
                << "Error: GenWeightSF needs channel and data source"
                << std::endl
                << "e.g.   GenWeight elnu DY"
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
        else if ( "tzq"  ==  dsN ){ d = tzq   ;}
        else if ( "ttb"  ==  dsN ){ d = ttb   ;}
        else if ( "ttl"  ==  dsN ){ d = ttl   ;}
        else if ( "ttj"  ==  dsN ){ d = ttj   ;}
        else if ( "tz1"  ==  dsN ){ d = tz1   ;}
        else if ( "tz2"  ==  dsN ){ d = tz2   ;}
        else if ( "wjt"  ==  dsN ){ d = wjt   ;}
        else if ( "wjx"  ==  dsN ){ d = wjx   ;}
        else if ( "ww"   ==  dsN ){ d = ww    ;}
        else if ( "wz"   ==  dsN ){ d = wz    ;}
        else if ( "zz"   ==  dsN ){ d = zz    ;}
        else if ( "st"   ==  dsN ){ d = st    ;}
        else if ( "stb"  ==  dsN ){ d = stb   ;}
        else if ( "stw"  ==  dsN ){ d = stw   ;}
        else if ( "stbw" ==  dsN ){ d = stbw  ;}
        else if ( "wzll" ==  dsN ){ d = wzll  ;}
        else if ( "wjqq" ==  dsN ){ d = wjqq  ;}
        else if ( "zjt1" ==  dsN ){ d = zjt1  ;}
        else if ( "zjt2" ==  dsN ){ d = zjt2  ;}
        else if ( "zjt3" ==  dsN ){ d = zjt3  ;}
        else if ( "zjqq" ==  dsN ){ d = zjqq  ;}
        else if ( "cms"  ==  dsN ){ d = cms   ;}
        else { std::cout << "Error: data source " << dsN
                << " not recognised" << std::endl ;
                return 4 ;
        }
                GenWeightSF(c,d) ;
                return 0 ;
}
