#include <boost/spirit/include/qi_core.hpp>
#include <boost/spirit/include/qi_eol.hpp>

#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_bind.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>

#include <iostream>
#include <string>
#include "types.hpp"

int main() {

    using boost::spirit::_1;
    using boost::spirit::_2;
    using boost::spirit::_3;
    using boost::spirit::_a;
    using boost::spirit::_b;
    using boost::spirit::_c;
    using boost::spirit::_val;
    using boost::spirit::double_;
    using boost::spirit::int_;
    using boost::spirit::uint_;

    using boost::spirit::qi::eol;
    using boost::spirit::qi::parse;
    using boost::spirit::qi::phrase_parse;
    using boost::spirit::qi::rule;
    using boost::spirit::qi::repeat;

    using boost::spirit::ascii::print;
    using boost::spirit::ascii::char_;
    using boost::spirit::ascii::space_type;
    using boost::spirit::ascii::space;

    using boost::phoenix::push_back;
    using boost::phoenix::val;
    using boost::phoenix::ref;
    using boost::phoenix::bind;

    std::string content = "HOOOOI";
    std::string title, description;

    typedef std::string::const_iterator binIt;
    binIt it = content.begin(), end = content.end();

    // title and description
    rule<binIt, types::String()> lineRule = *(char_ - eol) >> eol;
    if (! phrase_parse(
        it,
        end,
        lineRule[ref(title) = _1] >>
        lineRule[ref(description)  = _1],
        space
        )) throw types::parseError() << types::errinfo_parse("title or description");




    std::cout << "Hey\n";
    return 0;
}
