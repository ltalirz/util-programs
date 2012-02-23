#include <fstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <string>

#include <boost/spirit/include/qi.hpp>
#include <ctime>

#include "types.hpp"
#include "io.hpp"

using namespace std;

const string filename = "2MOL-18M.cube";
const int NFLOATS = 3000000;

void boostException(){
    std::string fileName = "notexisting", content;
    io::readFile(fileName, content);
}

void spiritVersion() {
    cout << "Using boost::spirit::qi::phrase_parse\n";

    time_t t = clock();

    //If we want to play with positions, need to open as binary
    ifstream file(filename.c_str(), ios_base::binary);

    // Set position to end of file
    file.seekg(0, ios_base::end);
    // Tell position
    ifstream::pos_type length = file.tellg();
    // go back to beginning
    file.seekg(0, ios_base::beg);

    // I tested reserving, but it seems for ifstream::read we must initialize
    vector<char> content(length);
    file.read( &content.front(), length );

    cout << "Time to read : " << (clock() -t)/1000.0 << " ms\n";
    t = clock();

    // Note: phrase_parse does back-insertion, so we do *not*
    // want to fill our vector with zeros
    vector<float> buf;
    buf.reserve(NFLOATS);
    length = buf.size();


    namespace sp = boost::spirit;
    vector<char>:: iterator begin = content.begin();
    vector<char>:: iterator end = content.end();

    sp::qi::phrase_parse(begin, end, *sp::qi::float_ , sp::ascii::space, buf);
    cout << "Time to parse : " << (clock() -t)/1000.0 << " ms\n";

    cout << buf[0] << endl;

}


void ifstreamVersion() {

    cout << "Using ifstream<float> iterators\n";
    time_t t = clock();
    ifstream file(filename.c_str());


    vector<float> numbers(3e6);
    cout << "Time to reserve vector : " << (clock() -t)/1000.0 << " ms\n";
    copy(istream_iterator<float>(file), istream_iterator<float>(), numbers.begin());
    //copy(istream_iterator<float>(file), istream_iterator<float>(), back_inserter(numbers));
    cout << "Time to parse : " << (clock() -t)/1000.0 << " ms\n";

    cout << numbers[0] << endl;
}


int main() {

//    spiritVersion();
//
//    ifstreamVersion();

   boostException();
    return 0;
}

