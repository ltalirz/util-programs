#include <boost/spirit/include/qi_core.hpp>
#include <boost/spirit/include/qi_eol.hpp>

#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_bind.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>

#include <fstream>

#include "formats/cube.h"
#include "atomistic/fundamental.h"
#include "io.h"
#include "la.h"
#include "types.hpp"


namespace formats {

using namespace types;

Uint Cube::countPoints() const {
    return grid.countPoints();
}

bool Cube::readCubeFile() {
    if(! this->fileName.empty() ) return this->readCubeFile(this->fileName);
    else throw types::runtimeError() << types::errinfo_runtime("Cannot read cube file - no file name given.");
}

// Simply adds the data of two cubefiles with the same grid
const Cube & Cube::operator+=(const Cube &c){
    this->grid += c.grid;
    return *this;
}

bool Cube::readCubeFile(String filename) {
    using boost::spirit::_1;
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

    using atomistic::Atom;
    types::Binary content;
    io::readBinary(filename, content);
    this->fileName = filename;
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
    rule<binIt, types::String()> lineRule = *(char_ - eol) >> eol;
    if (! phrase_parse(
        it,
        end,
        lineRule[ref(this->title) = _1] >>
        lineRule[ref(this->description)  = _1],
        space
        )) throw types::parseError() << types::errinfo_parse("title or description");

    // nat, grid
    using la::Direction;
    Uint nat;
    rule<binIt, std::vector<Real>(), space_type> vectorRule =
        repeat(3)[double_];
    rule<binIt, Direction(), space_type> directionRule =
        uint_[bind(&Direction::incrementCount, _val) = _1] >>
        vectorRule[bind(&Direction::incrementVector, _val) = _1];
    if(! phrase_parse(
        it,
        end,
        (uint_)[ref(nat)=_1]  >> vectorRule[ref(this->grid.originVector)=_1] >>
        repeat(3)[
            directionRule[push_back(ref(this->grid.directions), _1)]
        ],
        space
        )) throw types::parseError() << types::errinfo_parse("origin vector or grid");

    // atoms
    rule<binIt, Atom(), space_type> atomRule =
        uint_[bind(&Atom::number, _val) = _1] >>
        double_[bind(&Atom::charge, _val) = _1] >>
        vectorRule[bind(&Atom::coordinates, _val) = _1];
    if(! phrase_parse(
        it,
        end,
        repeat(val(nat))[
            atomRule[push_back(ref(this->atoms), _1)]
        ],
        space
        )) throw types::parseError() << types::errinfo_parse("list of atoms");

    // cube data
    Uint npoints = this->grid.countPoints();
    this->grid.data.reserve(npoints);
    if(! phrase_parse(
                        it,
                        end,
                        repeat(npoints)[double_],
                        space,
                        this->grid.data
        )) throw types::parseError() << types::errinfo_parse("data points");


    // Can map it as an Eigen::Vector3d if one likes
    //Map<Vector3d> originVector (&origin[0]);

    return true;
}


void Cube::print() const {
    Stream text;
    addHeader(text);
    //std::copy(text.begin(), text.end(), std::ostream_iterator<char>(std::cout));
    std::cout << text;
} 


void Cube::addHeader(Stream &stream) const {
    stream.append(title); stream += '\n';
    stream.append(description); stream += '\n';
    using boost::spirit::karma::right_align;
    using boost::spirit::karma::repeat;
    using boost::spirit::karma::int_;
    using boost::spirit::karma::columns;
    using boost::spirit::karma::eol;
    using boost::spirit::karma::generate;

    using atomistic::Atom;
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


bool Cube::writeCubeFile(String fileName) const {
    Stream data;
    this->addHeader(data);
    this->addData(data);
    
    return io::writeStream(fileName, data);
}

bool Cube::writeZProfile(String fileName) const {
        Stream data;
        this->addZProfile(data, "Z profile of cube file\n");
        return io::writeStream(fileName, data);
}

bool Cube::writeZProfile(String fileName, String header) const {
        Stream data;
        this->addZProfile(data, header);
        return io::writeStream(fileName, data);
}


/**
 * So far implemented only for cartesian grids with vectors
 * along x,y,z
 */
void Cube::addZProfile(Stream &stream, String header) const {
        using boost::spirit::karma::right_align;
        using boost::spirit::karma::double_;
        using boost::spirit::karma::eol;
        using boost::spirit::karma::generate;

        stream.append(header);
        stream.append( "z [a0]\t data\n");

        std::back_insert_iterator<Stream> sink(stream);
        std::vector<Real> data(grid.directions[2].incrementCount, 0.0);
        grid.averageXY(data);
        types::Real dZ = grid.directions[2].incrementVector[2];
        types::Real z = 0;

        std::vector<Real>::const_iterator dataIt = data.begin();
        // Fastest direction is z, stored in directions[2]
        while(dataIt != data.end()) {
                generate(sink, right_align(5)[double_] << '\t' << right_align(11)[types::sci5] << eol, z, *dataIt);
                ++dataIt;
                z += dZ;
        }
}

Uint Cube::countAtoms() const {
    return atoms.size();
}


bool Cube::readDescription(String filename) {
    std::ifstream file;
    file.open(filename.c_str());
    if (file.is_open()) {
        if(file.good()) std::getline(file, title);
        if(file.good()) std::getline(file, description);
    }
    else throw types::fileAccessError() << boost::errinfo_file_name(filename);

    this->fileName = filename;
   
    return true;
} 


bool WfnCube::readDescription(String filename) {
    Cube::readDescription(filename);
    
    using boost::spirit::uint_;
    using boost::spirit::_1;
    using boost::spirit::qi::phrase_parse;
    using boost::spirit::ascii::space;
    using boost::phoenix::ref;

    if( !phrase_parse(
                description.begin(),
                description.end(),
                "WAVEFUNCTION" >>
                uint_[ref(this->wfn) = _1] >>
                "spin" >>
                uint_[ref(this->spin) = _1],
                space))
        throw types::parseError() << types::errinfo_parse("The description of this cubefile does not contain information on the wave function.");

    return true;
}


bool WfnCube::readCubeFile(){
    return Cube::readCubeFile();
}


bool WfnCube::readCubeFile(String filename) {
    Cube::readCubeFile(filename);

    using boost::spirit::uint_;
    using boost::spirit::_1;

    using boost::spirit::qi::phrase_parse;

    using boost::spirit::ascii::space;

    using boost::phoenix::ref;

    if( !phrase_parse(
                description.begin(),
                description.end(),
                "WAVEFUNCTION" >>
                uint_[ref(this->wfn) = _1] >>
                "spin" >>
                uint_[ref(this->spin) = _1],
                space))
        throw types::parseError() << types::errinfo_parse("The description of this cubefile does not contain information on the wave function.");
    return true;
}

}

