//#include <algorithm>

#include <boost/spirit/include/karma_generate_attr.hpp>
#include <boost/spirit/include/karma_char.hpp>
#include <boost/spirit/include/karma_operator.hpp>
#include <boost/spirit/include/karma_numeric.hpp>


#include "karma.hpp"
#include <iostream>
#include <string>
#include <iterator>

void testColumnGenerator() {
using boost::spirit::karma::int_;
using boost::spirit::karma::space;
using boost::spirit::karma::generate;
using boost::spirit::karma::generate_delimited;

    std::vector<int> numbers;
    for(int i=0; i<100; ++i){
        numbers.push_back(i);
    }
    std::string s;
    std::back_insert_iterator<std::string> sink(s);
    generate_delimited(sink, custom_generator::columns[*int_], space, numbers);
    std::cout << s; 
}

int main() {
    testColumnGenerator();
    return 0;
}
