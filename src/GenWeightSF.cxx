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
	case wjt:{temp_opener="/data/disk3/nanoAOD_2017/WPlusJets_NanoAODv5/*.root";break;}/**/
	case zjt:{temp_opener=temp_header+"DYJetsToQQ"       +temp_footer;break;}// not downloaded yet
	case ttb:{temp_opener=temp_header+"TTToSemileptonic" +temp_footer;break;}
	case tz1:{temp_opener=temp_header+"ttZToQQ"          +temp_footer;break;}
	case tz2:{temp_opener=temp_header+"ttZToQQ_ext"      +temp_footer;break;}
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
                        case wjt:
                        case zjt:
                        case tz1:
                        case tz2:{           return mc__df;break;}
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
//              default  :throw std::invalid_argument(
//                      "Unimplemented ch (init)");
        }
	// make test runs faster by restriction. Real run should not
        auto dfr = df.Range(100000);// remember to enable MT when NOT range
        auto genWpos = df// remove one letter to do all
	.Filter("genWeight >= 0", "genPosCount")
	;
	genWpos.Report() ->Print();
