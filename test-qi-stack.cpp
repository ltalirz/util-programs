#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <algorithm>
#include <boost/spirit/include/qi.hpp>
#include <iostream>
#include <string>
using boost::spirit::ascii::space;
using boost::spirit::lit;
using boost::spirit::qi::eol;
using boost::spirit::qi::phrase_parse;

struct fix : std::unary_function<char, void> {
    fix(std::string &result) : result(result) {}
    void operator() (char c) {
        if      (c == '\n') result += "\\n";
        else if (c == '\r') result += "\\r";
        else                result += c;
    }
    std::string &result;
};

template <typename Parser>
void parse(const std::string &s, const Parser &p) {
    std::string::const_iterator it = s.begin(), end = s.end();
    bool r = phrase_parse(it, end, p, space - eol);
    std::string label;
    fix f(label);
    std::for_each(s.begin(), s.end(), f);
    std::cout << '"' << label << "\":\n" << "  - ";
    if (r && it == end) std::cout << "success!\n";
    else std::cout << "parse failed; r=" << r << '\n';
}

void hoden(){

	using boost::phoenix::ref;
using boost::spirit::qi::parse;
using boost::spirit::_1;
using boost::spirit::ascii::print;
	std::string s = "myline\n2ndline";
    std::string::const_iterator it = s.begin(), end = s.end();
  	std::string result;
	std::vector<char> test;
    parse(it, end, (*print >> eol)[ref(test) = _1]);
    std::cout << test[0] << std::endl;
    std::cout << std::string(test.begin(),test.end()) << std::endl;
}

// Read line into string
void readLine(){
using boost::spirit::qi::parse;
using boost::spirit::ascii::print;
	std::string s = "myline\n2ndline";
    std::string::const_iterator it = s.begin(), end = s.end();
  	std::string result;
    parse(it, end, *print >> eol, result);
    std::cout << result << std::endl;
    parse(it, end, *print, result);
    std::cout << result << std::endl;
}
// YEEEEEEESSSSS! It works!
void hoden4(){
using boost::spirit::qi::parse;
using boost::spirit::ascii::print;
	std::string s = "myline\n2ndline";
    std::string::const_iterator it = s.begin(), end = s.end();
  	std::string result;
    parse(it, end, *print, result);
    std::cout << result << std::endl;
}
// This version works also!
void hoden3(){
using boost::spirit::qi::parse;
	std::string s = "myline\n2ndline";
    std::string::const_iterator it = s.begin(), end = s.end();
  	std::string result;
    parse(it, end, lit("myline"), result);
    std::cout << result << std::endl;
}
// This version works!
void hoden2(){
	std::string s = "myline\n2ndline";
    std::string::const_iterator it = s.begin(), end = s.end();
  	std::string result;
    phrase_parse(it, end, lit("myline"), space-eol,result);
    std::cout << result;
}
int main() {
    parse("foo",     lit("foo"));
    parse("foo\n",   lit("foo") >> eol);
    parse("foo\r\n", lit("foo") >> eol);
    hoden();
}
