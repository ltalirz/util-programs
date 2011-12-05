#include <boost/format.hpp>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_bind.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/home/phoenix/bind/bind_function.hpp>

#include <blitz/array.h>

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

Uint EnergyLevels::count() const {
    return this->levels.size();
}

void EnergyLevels::print() const {
    std::cout << "Fermi energy: " << this->fermi << std::endl;
    for(std::vector<Real>::const_iterator it = this->levels.begin(); it != this->levels.end(); it++) {
        std::cout << *it << std::endl;
    }
}



void Spectrum::print() const {
    for(std::vector<EnergyLevels>::const_iterator it = this->spins.begin(); it != this->spins.end(); it++) {
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
    rule<binIt, Direction(), space_type> directionRule =
        uint_[bind(&Direction::incrementCount, _val) = _1] >>
        vectorRule[bind(&Direction::incrementVector, _val) = _1];
    phrase_parse(
        it,
        end,
        (uint_)[ref(nat)=_1]  >> vectorRule[ref(this->grid.originVector)=_1] >>
        repeat(3)[
            directionRule[push_back(ref(this->grid.directions), _1)]
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

    // cube data
    Uint npoints = this->grid.countPoints();
    this->grid.data.reserve(npoints);
    bool finddata = phrase_parse(
        it,
        end,
        repeat(npoints)[double_],
        space,
        this->grid.data
    );


    // Can map it as an Eigen::Vector3d if one likes
    //Map<Vector3d> originVector (&origin[0]);

}


void Grid::printHeader() const {
    using boost::format;
    format formatter("%5i %12.6f %12.6f %12.6f \n");

    std::vector<Direction>::const_iterator it;
    for(it = directions.begin(); it!= directions.end(); ++it) {
        std::cout << formatter % it->incrementCount % it->incrementVector[0] % it->incrementVector[1] % it->incrementVector[2];
    }

}

void Grid::printData() const {
    using boost::format;
    format formatter("%13.5e");

    std::vector<Real>::const_iterator it = data.begin();
    // Fastest direction is z, stored in directions[2]
    Uint i = 1, nZ = directions[2].incrementCount;
    while(it!= data.end()) {
        std::cout << formatter % *it;
	if(i % 6 == 0) std::cout << std::endl; 
	else if(i % nZ == 0){
		std::cout << std::endl;
		i=0;
	}
	++it;
	++i;
    }

}

void Grid::squareValues(){
    std::vector<Real>::iterator it = data.begin(), end = data.end();
    while(it != end){
        (*it) *= (*it);
	++it;
    }
}


void Grid::sumXY(std::vector<Real>& reduced) const {
    /** The Blitz++ way
    using namespace blitz;
    namespace t = tensor;
    
    Array<Real,3> dataArray(&data[0], shape(directions[0].incrementCount, directions[1].incrementCount, directions[2].incrementCount));
    // reduce second dimension
    Array<Real,2> reducedY(sum(dataArray(t::i,t::k,t::j), t::k));
    //std::cout << reduceY;
    Array<Real,1> reducedXY(sum(reducedY(t::j,t::i), t::j));
    std::cout << reduceXY;
    
     **/

    reduced = std::vector<Real> (directions[2].incrementCount, 0.0);
    std::vector<Real>::const_iterator itData=data.begin(), endData=data.end();
    std::vector<Real>::iterator itReduced=reduced.begin(), endReduced=reduced.end();
    // z is the fast index of the cube file, so we just need to sum
    // all of the z-compartments together
    while(itData != endData){
	if(itReduced == endReduced){
            itReduced = reduced.begin();
	}
        *itReduced += *itData;
	++itData;
	++itReduced;
    }
	    
    
    
}

Uint Grid::countPoints() const {
    std::vector<Direction>::const_iterator it;
    Uint points = 1;
    for(it = directions.begin(); it!= directions.end(); ++it) {
        points *= it->incrementCount;
    }
    return points;
}


void Atom::print() const {
    using boost::format;
    format formatter("%5i %12.6f %12.6f %12.6f %12.6f \n");
    std::cout << formatter % number % charge % coordinates[0] % coordinates[1] % coordinates[2];
}



void Cube::print() const {
    std::cout << String(title.begin(), title.end()) << std::endl;
    std::cout << String(description.begin(), description.end()) << std::endl;

    using boost::format;
    format formatter("%5i %12.6f %12.6f %12.6f \n");
    std::cout << formatter % countAtoms() % grid.originVector[0] % grid.originVector[1] % grid.originVector[2];
    grid.printHeader();
    for(std::vector<Atom>::const_iterator it = atoms.begin(); it != atoms.end(); ++it) {
        it->print();
    }
    // Sending the whole data over std::out is not very efficient
    //grid.printData();
}


Uint Cube::countAtoms() const {
    return atoms.size();
}


}

