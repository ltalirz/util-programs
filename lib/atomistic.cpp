#include <fstream>

#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_bind.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/home/phoenix/bind/bind_function.hpp>

#include <eigen3/Eigen/Dense>

#include "atomistic.h"
#include "io.h"

namespace atomistic {

typedef boost::sregex_iterator RegIt;
typedef boost::smatch Match;
typedef const boost::regex Regex;

using boost::lexical_cast;

EnergyLevels::EnergyLevels(std::vector<Real> levels, Real fermi): levels(levels), fermi(fermi) {
    this->sort();
};

void EnergyLevels::sort() {
    std::sort(this->levels.begin(), this->levels.end());
}

Uint EnergyLevels::count() {
    return this->levels.size();
}

void EnergyLevels::print() {
    std::cout << "Fermi energy: " << this->fermi << std::endl;
    for(std::vector<Real>::iterator it = this->levels.begin(); it != this->levels.end(); it++) {
        std::cout << *it << std::endl;
    }
}



void Spectrum::print() {
    for(std::vector<EnergyLevels>::iterator it = this->spins.begin(); it != this->spins.end(); it++) {
        it->print();
    }
}


void Spectrum::readFromCp2k(String filename) {
    String content;
    io::readFile(filename, content);

    RegIt occIt(
        content.begin(),
        content.end(),
        boost::regex("Eigenvalues of the occupied subspace spin([\\.\\-\\d\\s]*)"));
    RegIt fermiIt(
        content.begin(),
        content.end(),
        boost::regex("Fermi Energy \\[eV\\] :([\\.\\-\\d\\s]*)"));
    RegIt unoccIt(
        content.begin(),
        content.end(),
        boost::regex("Eigenvalues of the unoccupied subspace spin.*?iterations([\\.\\-\\d\\s]*)"));
    RegIt dummyIt;


    Regex fermiRegex("\\-?\\d\\.\\d{6}");
    Regex levelRegex("\\-?\\d\\.\\d{8}");

    // Iterate over spins
    while(occIt != dummyIt) {
        String levelData = occIt->str().append(unoccIt->str());
        std::vector<Real> levels;
        RegIt levelIt(levelData.begin(), levelData.end(), levelRegex);
        while(levelIt != dummyIt) {
            levels.push_back(lexical_cast<Real>(levelIt->str()));
            ++levelIt;
        }

        Match fermiMatch;
        regex_search(fermiIt->str(), fermiMatch, fermiRegex);
        Real fermi = lexical_cast<Real>(fermiMatch);

        EnergyLevels energyLevels = EnergyLevels(levels, fermi);
        this->spins.push_back(energyLevels);

        ++occIt;
        ++unoccIt;
        ++fermiIt;
    }

}


void Cube::readCubeFile(String filename) {
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
    using boost::spirit::qi::lexeme;
    using boost::spirit::qi::locals;

    using boost::spirit::ascii::print;
    using boost::spirit::ascii::char_;
    using boost::spirit::ascii::space_type;
    using boost::spirit::ascii::space;

    using boost::phoenix::push_back;
    using boost::phoenix::val;
    using boost::phoenix::ref;
    using boost::phoenix::bind;

    using Eigen::MatrixXd;
    using Eigen::Vector3d;
    using Eigen::Map;

    io::Binary content;
    io::readBinary(filename, content);
    typedef io::Binary::const_iterator binIt;
    binIt it = content.begin(), end = content.end();

    /**
    *     THE FORMAT IS AS FOLLOWS (LAST CHECKED AGAINST GAUSSIAN 98):
    *
    *     LINE   FORMAT      CONTENTS
    *     ===============================================================
    *      1     A           TITLE
    *      2     A           DESCRIPTION OF PROPERTY STORED IN CUBEFILE
    *      3     I5,3F12.6   #ATOMS, X-,Y-,Z-COORDINATES OF ORIGIN
    *      4-6   I5,3F12.6   #GRIDPOINTS, INCREMENT VECTOR
    *      #ATOMS LINES OF ATOM COORDINATES:
    *      ...   I5,4F12.6   ATOM NUMBER, CHARGE, X-,Y-,Z-COORDINATE
    *      REST: 6E13.5      CUBE DATA (WITH Z INCREMENT MOVING FASTEST, THEN
    *                        Y AND THEN X)
    *
    *     FOR ORBITAL CUBE FILES, #ATOMS WILL BE < 0 AND THERE WILL BE ONE
    *     ADDITIONAL LINE AFTER THE FINAL ATOM GIVING THE NUMBER OF ORBITALS
    *     AND THEIR RESPECTIVE NUMBERS. ALSO THE ORBITAL NUMBER WILL BE
    *     THE FASTEST MOVING INCREMENT.
    *
    *     ALL COORDINATES ARE GIVEN IN ATOMIC UNITS.
    */

    // title and description
    rule<binIt, io::Binary()> lineRule = *(char_ - eol) >> eol;
    phrase_parse(
        it,
        end,
        lineRule[ref(this->title) = _1] >>
        lineRule[ref(this->description)  = _1],
        space
    );

    // nat, grid
    Uint nat;
    rule<binIt, std::vector<Real>(), space_type> vectorRule =
        repeat(3)[double_];
    phrase_parse(
        it,
        end,
        (uint_)[ref(nat)=_1]  >> vectorRule[ref(this->grid.originVector)=_1] >>
        repeat(3)[
            (uint_)[push_back(ref(this->grid.incrementCounts), _1)]  >>
            vectorRule[push_back(ref(this->grid.incrementVectors), _1)]
        ],
        space
    );
    
    // atoms
    rule<binIt, Atom(), space_type> atomRule =
        uint_[bind(&Atom::number, _val) = _1] >>
        double_[bind(&Atom::charge, _val) = _1] >>
        vectorRule[bind(&Atom::coordinates, _val) = _1];
    phrase_parse(
        it,
        end,
        repeat(val(nat))[
            atomRule[push_back(ref(this->atoms), _1)]
        ],
        space
    );

    // Can map it as an Eigen::Vector3d if one likes
    //Map<Vector3d> originVector (&origin[0]);

    this->print();

    // Todo:
    // 1. Write parser for cube files using qi
    // 2. Convert cube file to Eigen Matrix (or some Eigen object that can be used
    // with phrase_parse)

}


void Grid::print() {
    std::cout << this->originVector[0] << std::endl;
}

void Atom::print() {
    std::cout << this->number << coordinates[0] << coordinates[1] << std::endl;
}

void Cube::print() {
    // Todo: Learn boost format and print out cube file as it was
    std::cout << String(this->title.begin(), this->title.end()) << std::endl;
    std::cout << String(this->description.begin(), this->description.end()) << std::endl;

    this->grid.print();
    for(std::vector<Atom>::iterator it = this->atoms.begin(); it != this->atoms.end(); it++) {
        it->print();
    }
}
}
