#include <iostream>
#include <iterator>
#include <vector>
#include <string>
#include <boost/program_options.hpp>

namespace po = boost::program_options;
using std::cout;
using std::vector;
using std::string;

// Returns true, if parsing went ok
bool parse(int ac, char* av[], po::variables_map& vm);


int main(int ac, char* av[]){
	
	po::variables_map args;
	if (parse(ac, av, args)){
		vector<string>::const_iterator iter = 
			args["input-file"].as< vector<string> >().begin();
//	cout << args["input-file"].as< vector<string> >();
	}
	return 0;
}


bool parse(int ac, char* av[], po::variables_map& vm){
	
	// Declare regular options
	po::options_description desc("Allowed options");
	desc.add_options()
    		("help,h", "produce help message")
    		("version,v", "print version information")
		("optimization,O", po::value<int>(), "set optimization level")
		("input-file", po::value< vector<string> >(),
		 "input file")
	;
	
	// Declare positional options
	po::positional_options_description p;
	p.add("input-file", -1);

	// Parse
	po::store(po::command_line_parser(ac,av).
			options(desc).positional(p).run(), vm);
	po::notify(vm);    

	if (vm.count("help")) {
		cout << "Usage: <input-file> [options]\n";
		cout << desc << "\n";
	} else if (vm.count("version")) {
		cout << "V1.0, Nov 22nd 2011\n";
	} else if (!vm.count("input-file")) {
		cout << "Please specify input file\n";
	} else {
		return true;
	}

	return false;
}

