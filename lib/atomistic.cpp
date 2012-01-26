#include <boost/format.hpp>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_bind.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/home/phoenix/bind/bind_function.hpp>

#include "atomistic.h"
#include "io.h"
#include "la.h"


namespace atomistic {

typedef boost::sregex_iterator RegIt;
typedef boost::smatch Match;
typedef const boost::regex Regex;

using boost::lexical_cast;
using namespace types;


EnergyLevels::EnergyLevels(std::vector<Real> levels, Real fermi): levels(levels), fermi(fermi) {
    this->sort();
};

void EnergyLevels::sort() {
    std::sort(this->levels.begin(), this->levels.end());
}

Uint EnergyLevels::count() const {
    return this->levels.size();
}

Real EnergyLevels::getLevel(Uint i) const {
    if( i > 0 && i <= this->levels.size() ){
        return this->levels[i-1];
    } else {
        throw std::range_error("Level index out of bounds.");
        return 0;
    }
}

void EnergyLevels::print() const {
    std::cout << "Fermi energy: " << this->fermi << std::endl;
    for(std::vector<Real>::const_iterator it = this->levels.begin(); it != this->levels.end(); it++) {
        std::cout << *it << std::endl;
    }
}

void EnergyLevels::shift(Real deltaE) {
    for(std::vector<Real>::iterator it = this->levels.begin(); it != this->levels.end(); it++) {
        *it += deltaE;
    }
    fermi += deltaE;
}


void Spectrum::shift(Real deltaE) {
    for(std::vector<EnergyLevels>::iterator it = this->spins.begin(); it != this->spins.end(); it++) {
        it->shift(deltaE);
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

Uint Cube::countPoints() const {
    return grid.countPoints();
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

    types::Binary content;
    io::readBinary(filename, content);
    typedef types::Binary::const_iterator binIt;
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
    rule<binIt, types::Binary()> lineRule = *(char_ - eol) >> eol;
    phrase_parse(
        it,
        end,
        lineRule[ref(this->title) = _1] >>
        lineRule[ref(this->description)  = _1],
        space
    );

    // nat, grid
    using la::Direction;
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


void Cube::print() const {
    Stream text;
    addHeader(text);
    //std::copy(text.begin(), text.end(), std::ostream_iterator<char>(std::cout));
    std::cout << text;
} 


void Cube::addHeader(Stream &stream) const {
    stream.append(title.begin(), title.end()); stream += '\n';
    stream.append(description.begin(), description.end()); stream += '\n';
    using boost::spirit::karma::right_align;
    using boost::spirit::karma::repeat;
    using boost::spirit::karma::int_;
    using boost::spirit::karma::columns;
    using boost::spirit::karma::eol;
    using boost::spirit::karma::generate;

    std::back_insert_iterator<Stream> sink(stream);
    // Origin vector
    generate(sink, 
            right_align(5)[int_] << repeat(3)[right_align(12)[types::real6]] << eol,
            countAtoms(),
            grid.getOriginVector());
// Grid vectors    
        std::vector<la::Direction>::const_iterator it;
        const std::vector<la::Direction> &directions = grid.getDirections();
    for(it = directions.begin(); it!= directions.end(); ++it) {
        generate(sink, 
                right_align(5)[int_] << repeat(3)[right_align(12)[types::real6]] << eol,
                it->incrementCount, 
                it->incrementVector);
    }
    // Atoms
    for(std::vector<Atom>::const_iterator it = atoms.begin(); it != atoms.end(); ++it) {
    generate(sink, 
            right_align(5)[int_] << right_align(12)[types::real6] << repeat(3)[right_align(12)[types::real6]] << eol,
            it->getNumber(),
            it->getCharge(),
            it->getCoordinates());
    }
}


void Cube::addData(Stream &stream) const {
    using boost::spirit::karma::right_align;
    using boost::spirit::karma::repeat;
    using boost::spirit::karma::int_;
    using boost::spirit::karma::columns;
    using boost::spirit::karma::eol;
    using boost::spirit::karma::generate;

    std::back_insert_iterator<Stream> sink(stream);
    const std::vector<Real> &data = grid.getData();
    std::vector<Real>::const_iterator dataIt = data.begin();
    // Fastest direction is z, stored in directions[2]
    const std::vector<la::Direction> &directions = grid.getDirections();
    Uint i = 1, nZ = directions[2].incrementCount;
    while(dataIt != data.end()) {
        generate(sink, right_align(13)[types::sci5], *dataIt);
        if(i % 6 == 0) generate(sink, eol);
        else if(i % nZ == 0) {
            generate(sink, eol);
            i=0;
        }
        ++dataIt;
        ++i;
    }

//    // One could do the column logic completely within boost::spirit::karma
//    // However I have checked that the performance gain by the version below
//    // (which does not put line endings after nZ) is marginal (< 10%)
//    generate(sink, columns(6)[*right_align(13)[types::sci5]], data);
}


void Cube::writeCubeFile(String fileName) const {
    Stream data;
    this->addHeader(data);
    this->addData(data);
    io::writeStream(fileName, data);
}


Uint Cube::countAtoms() const {
    return atoms.size();
}

void WfnCube::readCubeFile(String filename) {
    Cube::readCubeFile(filename);

    using boost::spirit::uint_;
    using boost::spirit::_1;

    using boost::spirit::qi::phrase_parse;

    using boost::spirit::ascii::space;

    using boost::phoenix::ref;

    typedef types::Binary::const_iterator binIt;
    binIt it = description.begin(), end = description.end();
    bool success = phrase_parse(
                       it,
                       end,
                       "WAVEFUNCTION" >>
                       uint_[ref(this->wfn) = _1] >>
                       "spin" >>
                       uint_[ref(this->spin) = _1],
                       space);
    if (!success) {
        throw std::runtime_error("The description of this cubefile does not contain information on the wave function.");
//}
    }
}

}

