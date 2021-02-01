//#include "analyse.hpp"
#include "badbranches.hpp"
#include "dedupe.hpp"
#include "calchisto.hpp"

#include <iostream>
#include <string>

int main(int argc,char* argv[]){
	if(argc < 2){
		std::cout << "Error: no command provided" << std::endl;
		return 1;
	}
	     if(const auto command = std::string_view(argv[1]);
//	        command ==     "analyse"){  }//analyse(    argc,argv);}
//	else if(command == "badbranches"){badbranches();}//argc,argv);}
	        command == "badbranches"){badbranches();}//argc,argv);}
	else if(command ==      "dedupe"){dedupe();}//     argc,argv);}
	else if(command ==   "calchisto"){
		if(argc < 4){
			std::cout
			<< "Error: calchisto now needs channel and data source"
			<< std::endl
			<< "e.g. eta calchisto elnu tzq"
			<< std::endl
			;
			return 2;
		}
		channel ch; dataSource ds;
		     if(const auto chN = std::string_view(argv[2]);
		        "elnu"  == chN){ch = elnu;}
		else if("munu"  == chN){ch = munu;}
		else{std::cout << "Error: channel " << chN
		               << " not recognised" << std::endl;
		     return 3;
		}
		     if(const auto dsN = std::string_view(argv[3]);
		        "tzq"   == dsN){ds = tzq;}
		else if( "ww"   == dsN){ds =  ww;}
		else if( "wz"   == dsN){ds =  wz;}
		else if( "zz"   == dsN){ds =  zz;}
		else if("ttb"   == dsN){ds = ttb;}
		else if("tz1"   == dsN){ds = tz1;}
		else if("tz2"   == dsN){ds = tz2;}
		else if("met"   == dsN){ds = met;}
		else if("cms"   == dsN){ds = cms;}
		else{std::cout << "Error: data source " << dsN
		               << " not recognised    " << std::endl;
		     return 4;}
		calchisto(ch,ds);
	}
	else{
		std::cout
		<< "Error:  command "
		<<          command
		<< " not recognised "
		<< std::endl
		;
		return 5;
	}
	return 0;
}
