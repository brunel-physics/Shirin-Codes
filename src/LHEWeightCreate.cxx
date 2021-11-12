// This script is designed to take WW WZ ZZ datasets
// which lack the branch LHEPDFWeight and create a dummy branch with value 1 .
// Compiler:
//clang++ -Isrc -std=c++17 -march=native -pipe -O3 -Wall -Wextra -Wpedantic -o build/lhe src/LHEWeightCreate.cxx `root-config --libs` -lm

#include <ROOT/RDataFrame.hxx>//#include <ROOT/RCsvDS.hxx>
#include <Math/Vector4D.h>
#include <TRandom3.h>// used Gaussian, uniform each once
//#include <execution>// need to link -ltbb in Makefile
#include <TChain.h>
#include <TF1.h>

enum      dataSource      {ww,wz,zz};
constexpr dataSource
          dataSourceAll[]={ww,wz,zz};

enum      channel      {elnu,munu};
constexpr channel
          channelAll[]={elnu,munu};


using doubles = ROOT::VecOps::RVec<double>;
using  floats = ROOT::VecOps::RVec<float>;
using    ints = ROOT::VecOps::RVec<int>;
using   bools = ROOT::VecOps::RVec<bool>;
using strings = ROOT::VecOps::RVec<std::string>;

std:: string temp_header1="/nfs/data/eepgssg/",
             temp_footer ="/*.root",
             temp_opener;

auto dummy(floats jet_pt){
floats dummys(jet_pt.size(),1.);
return dummys;
}

void LHEWeightCreate(const channel ch, const dataSource ds){
// open ww, wz ,zz
ROOT::EnableImplicitMT(4);// SYNC WITH CONDOR JOBS!
        // Open LB file even if Monte Carlo will NOT use it
switch (ds){
        case   ww:{temp_opener=temp_header1+ "WW"       +temp_footer;break;}
        case   wz:{temp_opener=temp_header1+ "WZ"       +temp_footer;break;}
        case   zz:{temp_opener=temp_header1+ "ZZ"       +temp_footer;break;}
        default :throw std::invalid_argument("Unimplemented ds (rdfopen)");
}
ROOT::RDataFrame mc__df("Events",temp_opener);// Monte Carlo
// define dummy variable
auto dfr = mc__df
	 .Define("LHEPdfWeight",dummy,{"Jet_pt"})
	 ;
switch (ds){
	case   ww:{temp_opener="/nfs/data/eepgssg/WW/WW_v7/file.root";break;}
        case   wz:{temp_opener="/nfs/data/eepgssg/WZ/WZ_v7/file.root";break;}
        case   zz:{temp_opener="/nfs/data/eepgssg/ZZ/ZZ_v7/file.root";break;}
        default :throw std::invalid_argument("Unimplemented ds(temp_open)");
}
std::cout<<"temp_opener is "<<temp_opener<<std::endl;
//snap shot in the nfs/data/eepgssg/WW/WW_ dir WZ and ZZ respectively too
auto SnapRDF = dfr.Snapshot("Events",temp_opener);// SNAPPED!
}// void
int main ( int argc , char *argv[] ){
        if ( argc < 2 ) {
                std::cout << "Error: no command provided" << std::endl ;
                return 1 ;
        }
        if ( argc < 3 ) {
                   std::cout
                << "Error: NPL_run needs channel and data source"
                << std::endl
                << "e.g.   NPL_run elnu DY"
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
        else if ( "ww"  ==  dsN ){ d = ww   ;}
        else if ( "wz"  ==  dsN ){ d = wz   ;}
        else if ( "zz"  ==  dsN ){ d = zz   ;}
        else { std::cout << "Error: data source " << dsN
                << " not recognised" << std::endl ;
                return 4 ;
        }
                LHEWeightCreate(c,d) ;
                return 0 ;
}
