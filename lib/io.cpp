#include <fstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <string>
#include <boost/spirit/include/qi.hpp>

#include "io.h"
#include "types.h"

namespace io {

using namespace types;

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

template<typename T>
bool parse(const Binary& content, std::vector<T> parsed) {

    // Note: phrase_parse does back-insertion, so we do *not*
    // want to fill our vector with zeros. We might reserve some space however
    //buf.reserve(content.size());

    using namespace boost::spirit;
    Binary::iterator begin = content.begin();
    Binary::iterator end = content.end();

    if(qi::phrase_parse(begin, end, qi::float_ , ascii::space, parsed))
        return true;

    return false;
}


}
