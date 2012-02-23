#include <fstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <string>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>


#include "io.hpp"
#include "types.hpp"

namespace io {

using namespace types;
namespace fs = boost::filesystem;

bool readFile(String filename, String& content) {
    std::ifstream file;
    file.open(filename.c_str());
    if (file.is_open()) {
        String line;

        while (file.good()) {
            std::getline(file,line);
            content += line;
        }
        file.close();
        return true;

    }
    else throw types::fileAccessError() << boost::errinfo_file_name(filename);

    return false;
}

bool readBinary(String filename, Binary& content) {
    // Open as binary
    using namespace std;
    ifstream file(filename.c_str(), ios_base::binary);
    if (file.is_open()) {
        // Set position to end of file
        file.seekg(0, ios_base::end);
        // Tell position
        ifstream::pos_type length = file.tellg();
        // go back to beginning
        file.seekg(0, ios_base::beg);

        // I tested reserving, but it seems for ifstream::read we must initialize
        content.resize(length);
        file.read( &content.front(), length );

        file.close();

        return true;
    }
    else throw types::fileAccessError() << boost::errinfo_file_name(filename);

    return false;
}

bool writeBinary(String filename, const Binary& data) {
    // Open as binary
    using namespace std;
    ofstream file(filename.c_str(), ios_base::binary);
    if (file.is_open()) {
        // Get length of data
        ifstream::pos_type length = data.size();
        file.write( &data.front(), length );

        file.close();

        return true;
    }
    else throw types::fileAccessError() << boost::errinfo_file_name(filename);

    return false;
}

bool writeStream(String filename, const Stream& data) {
    // Open as binary
    using namespace std;
    ofstream file(filename.c_str(), ios_base::binary);
    if (file.is_open()) {
        // Get length of data
        ifstream::pos_type length = data.size();
    
        file.write( data.c_str(), length );

        file.close();

        return true;
    }
    else throw types::fileAccessError() << boost::errinfo_file_name(filename);

    return false;
}

types::String getAbsolutePath(types::String relativePath){
    fs::path fullPath( fs::initial_path<fs::path>());
    fullPath = fs::system_complete( fs::path( relativePath ));
    if( !fs::exists( fullPath ) )
        throw types::fileAccessError() << boost::errinfo_file_name(fullPath.string());
    else return fullPath.string();
}

types::String getFileName(types::String path){
    fs::path p(path);
    return p.filename().string();
}

types::String getFileExtension(types::String path){
    fs::path p(path);
    return p.extension().string();
}

}
