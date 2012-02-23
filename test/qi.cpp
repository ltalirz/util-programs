//#include <algorithm>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_bind.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>


#include <boost/spirit/include/qi_core.hpp>
#include <boost/spirit/include/qi_eol.hpp>

#include <iostream>
#include <string>

#include "atomistic.hpp"
#include "io.hpp"
#include "types.hpp"

using boost::spirit::ascii::space;
using boost::spirit::ascii::print;
using boost::spirit::lit;
using boost::spirit::qi::eol;
using boost::spirit::qi::parse;
using boost::spirit::qi::phrase_parse;
using boost::spirit::qi::rule;
using boost::spirit::qi::double_;
using boost::spirit::qi::int_;
using boost::spirit::qi::repeat;
using boost::spirit::ascii::space_type;
using boost::spirit::_val;
using boost::spirit::_1;

using boost::phoenix::bind;
using boost::phoenix::ref;


// The space skipper eats line endings
// before they are given to the actual parser.
void parseNumbers(){
    // Read list of bias voltages

    std::string biasListFileName = "list";
    types::Binary biasBin;
    io::readBinary(biasListFileName, biasBin);
    typedef types::Binary::const_iterator binIt;
    binIt it = biasBin.begin(), end = biasBin.end();
    std::vector<types::Real> biasList;
    if (! phrase_parse(
                it,
                end,
                *double_,
                space,
                biasList
                )) throw types::parseError() << types::errinfo_parse("list of bias voltages");
    std::cout << biasList[0] << biasList[1];
}



// Parse into string instead of std::vector<char>
void toString() {
    std::string s = "Hi there";
    std::string parsed;

    phrase_parse(
        s.begin(),
        s.end(),
        *print,
        space,
        parsed
    );
    std::cout << parsed;
}



// Pass attribute to class member using phoenix::bind
void readDoubles() {
	using namespace atomistic;
    std::string s = "123.0 456 789";
    std::vector<char> v(s.begin(), s.end());
    std::vector<char>::const_iterator it = v.begin(), end = v.end();

    std::vector<Atom> mouden;
    rule<std::vector<char>::const_iterator, Atom(), space_type> vectorRule =  (*double_)[bind(&Atom::coordinates, _val) =  _1];
    phrase_parse(
        it,
        end,
        vectorRule,
        space,
        mouden
    );
    std::cout << mouden[0].coordinates[0] << std::endl;
}



// Read line from char-vector into string using phoenix
// and separately declared rule
void readPhoenixRule() {
    std::string s = "myline\n2ndline";
    std::vector<char> v(s.begin(), s.end());
    std::vector<char>::const_iterator it = v.begin(), end = v.end();
    std::vector<char> result;

    // Although we want to stuff the whole match into a ref, we don't
    // need brackets here, if we declare the rule separately
    rule<std::vector<char>::const_iterator, std::vector<char>()> lineRule = *print >> eol;
    parse(it, end, lineRule[ref(result) = _1]);
    std::cout << std::string(result.begin(), result.end()) << std::endl;
}


// Read line from char-vector into string using phoenix
void readPhoenixVector() {
    std::string s = "myline\n2ndline";
    std::vector<char> v(s.begin(), s.end());
    std::vector<char>::const_iterator it = v.begin(), end = v.end();
    std::vector<char> result;

    // One could be tempted to do [result = _1] but this doesn't work here.
    rule<std::vector<char>::const_iterator> lineRule = (*print >> eol)[ref(result) = _1];
    parse(it, end, lineRule);
    std::cout << std::string(result.begin(), result.end()) << std::endl;
}

// Read line from string into string using phoenix
void readPhoenix() {

    std::string s = "myline\n2ndline";
    std::string::const_iterator it = s.begin(), end = s.end();
    std::vector<char> result;

    // One could be tempted to do [result = _1] but this doesn't work here.
    rule<std::string::const_iterator> lineRule = (*print >> eol)[ref(result) = _1];
    parse(it, end, lineRule);
    std::cout << std::string(result.begin(), result.end()) << std::endl;
}

// Read line from string into string using attribute
void readAttribute() {

    std::string s = "myline\n2ndline";
    std::string::const_iterator it = s.begin(), end = s.end();
    std::string result;
    // parse is doing something like push_back(result, match).
    // Apparently we can use a string directly here.
    parse(it, end, *print >> eol, result);
    std::cout << result << std::endl;
}

int main() {
//    readPhoenix();
//    readPhoenixVector();
//    readAttribute();
//    readPhoenixRule();
//    readDoubles();
//    toString();
    parseNumbers();
    return 0;
}
