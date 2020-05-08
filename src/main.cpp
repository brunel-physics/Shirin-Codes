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
	     if(const auto command  = std::string_view(argv[1]);
	        command ==     "analyse"){  }//analyse(    argc,argv);}
	else if(command == "badbranches"){badbranches();}//argc,argv);}
	else if(command ==      "dedupe"){dedupe();}//     argc,argv);}
	else if(command ==   "calchisto"){calchisto(elnu,tzq);}
	else{
		std::cout << "Error: command " 
		          <<         command 
			  << " not recognised" 
			  << std::endl;
		return 2;
	}
	return 0;
}
