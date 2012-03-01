#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <iostream>
#include <fstream>
#include <iterator>


// Reading from configuration file (default: po.inp)
void configTest(int ac, char* av[]){

        std::string config_file;
        
        po::options_description cmdline_options("Generic options");
        cmdline_options.add_options()
            ("input,i",
             po::value<std::string>(&config_file)->default_value("po.inp"),
             "name of a file of a configuration.")
            ("test-string", po::value< std::string >(), "test option")
            ("test-double", po::value< double >(), "test option")
            ("test-uint", po::value< unsigned int >(), "test option")
            ;
    
        po::variables_map vm;
        store(po::command_line_parser(ac, av).
              options(cmdline_options).run(), vm);
        notify(vm);
        
        std::ifstream ifs(config_file.c_str());
        if (!ifs) std::cout << "can not open config file: " 
            << config_file << "\n";
        else  {
            store(parse_config_file(ifs, cmdline_options), vm);
            notify(vm);
        }
    
        if (vm.count("test-string")) std::cout << "Test option " << 
            vm["test-string"].as< std::string >() << " specified\n";
        if (vm.count("test-uint")) std::cout << "Test option " << 
            vm["test-uint"].as< unsigned int >() << " specified\n";
        if (vm.count("test-double")) std::cout << "Test option " << 
            vm["test-double"].as< double>() << " specified\n";
}


int main(int ac, char* av[]){
    configTest(ac, av);

    return 0;
}


