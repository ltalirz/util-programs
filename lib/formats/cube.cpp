#include <boost/spirit/include/qi_core.hpp>
#include <boost/spirit/include/qi_eol.hpp>

#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_bind.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>

#include <blitz/array.h>

#include <fstream>

#include "formats/cube.hpp"
#include "atomistic/fundamental.hpp"
#include "io.hpp"
#include "la.hpp"
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
    this->atoms.clear();
    this->grid = CubeGrid();

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
        uint_[bind(&Direction::nElements, _val) = _1] >>
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
                it->nElements, 
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
    Uint i = 1, nZ = directions[2].nElements;
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
        std::vector<Real> data;
        this->averageXY(data);
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


/**
 * 3d wrapper for general stride function
 */
void CubeGrid::stride(types::Uint sX, types::Uint sY, types::Uint sZ) {
    std::vector<Uint> tempStrides;
    tempStrides.push_back(sX);
    tempStrides.push_back(sY);
    tempStrides.push_back(sZ);
    stride(tempStrides);
}

/**
 * 3d wrapper for general resize function
 */
void CubeGrid::resize(types::Uint nX, types::Uint nY, types::Uint nZ) {
    std::vector<Uint> tempCounts;
    tempCounts.push_back(nX);
    tempCounts.push_back(nY);
    tempCounts.push_back(nZ);
    resize(tempCounts);
}

// 3d wrapper for nd getNearestDataPoint
Real CubeGrid::getNearestDataPoint(Real x, Real y, Real z) const {
    std::vector<Real> coordinates;
    coordinates.push_back(x);
    coordinates.push_back(y);
    coordinates.push_back(z);

    return getNearestDataPoint(coordinates);
}


// 3d wrapper for nd getDataPoint
Real CubeGrid::getDataPoint(Uint x, Uint y, Uint z) const {
    std::vector<Uint> indices;
    indices.push_back(x);
    indices.push_back(y);
    indices.push_back(z);

    return getDataPoint(indices);
}



//types::Real CubeGrid::interpolateDataPoint(types::Real x, types::Real y, types::Real z) const {
//    std::vector<Real> coordinates, i;
//    coordinates.push_back(x);
//    coordinates.push_back(y);
//    coordinates.push_back(z);
//    getFractionalCoordinates(coordinates, i);
//
//
//    return  getDataPoint( Uint(i[0]), Uint(i[1]), Uint(i[2])) * (1-i[0]) * (1-i[1]) * (1-i[2])
//};
//
//void CubeGrid::interpolatedZPlane(const std::vector<Real> &zProfile,
//            std::vector<types::Real> &plane){
//    std::vector<Uint> indices;
//    std::vector<Real>::const_iterator it = zProfile.begin(),
//          end = zProfile.end();
//
//   while(it != end){
//       getNearest
//

void CubeGrid::zPlane(Uint zIndex, std::vector<types::Real> &plane) {
    using namespace blitz;
    std::vector<Real> test;
    const Array<Real,3> dataArray(
        &data[0],
        shape(directions[0].getNElements(),
              directions[1].getNElements(),
              directions[2].getNElements()),
        neverDeleteData);
    Array<types::Real,2> blitzPlane =
        dataArray(Range::all(), Range::all(), zIndex);
    plane.assign(blitzPlane.begin(), blitzPlane.end());
}


/**
 * Fill vector reduced with the sum over XY
 */
void CubeGrid::sumXY(std::vector<Real>& reduced) const {
    /** The Blitz++ way
    using namespace blitz;
    namespace t = tensor;

    Array<Real,3> dataArray(&data[0], shape(directions[0].nElements, directions[1].nElements, directions[2].nElements));
    // reduce second dimension
    Array<Real,2> reducedY(sum(dataArray(t::i,t::k,t::j), t::k));
    //std::cout << reduceY;
    Array<Real,1> reducedXY(sum(reducedY(t::j,t::i), t::j));
    std::cout << reduceXY;

     **/

    std::vector<Real>::const_iterator itData=data.begin(), endData=data.end();
    reduced = std::vector<Real>(directions[2].getNElements(), 0.0);
    std::vector<Real>::iterator itReduced=reduced.begin(), endReduced=reduced.end();
    // z is the fast index of the cube file, so we just need to sum
    // all of the z-compartments together
    while(itData != endData) {
        if(itReduced == endReduced) {
            itReduced = reduced.begin();
        }
        *itReduced += *itData;
        ++itData;
        ++itReduced;
    }


}


void CubeGrid::averageXY(std::vector<Real>& reduced) const {
    sumXY(reduced);
    Uint points = countPoints() / directions[2].getNElements();

    std::vector<Real>::iterator it;
    for(it = reduced.begin(); it!= reduced.end(); ++it) {
        *it /= points;
    }

}

Real Cube::topZCoordinate() const{
    
    // Get highest z coordinate
    std::vector< atomistic::Atom >::const_iterator it = atoms.begin(),
        end = atoms.end();
    types::Real zTop = it->coordinates[2];
    ++it;
    while(it != end) {
        if( it->coordinates[2] > zTop ) {
            zTop = it->coordinates[2];
        }
        ++it;
    }

    return zTop;
}

void CubeGrid::zIsoSurface(
        types::Real isoValue, 
        std::vector<types::Real> &surface) const {

    Uint zPoints = directions[2].getNElements();
    Real dZ = directions[2].getIncrementVector()[2];
    Uint planePoints = directions[0].getNElements()
        * directions[1].getNElements();
    surface.reserve(planePoints);
    
    
    std::vector<Real>::const_iterator dataIt, dataStop;
    for(Uint p = 0; p < planePoints; ++p){
        // Start from high z, going down
        dataIt = data.begin();
        dataIt += (p+1) * zPoints - 1;
        Uint z = zPoints -1;
        dataStop = data.begin();
        dataStop += p * zPoints;
        while(dataIt != dataStop){
            // As soon as we enter the isosurface
            if(*dataIt > isoValue){
                // If we are still at the top of the grid
                if(z == zPoints-1) surface.push_back(z * dZ);
                // else we need to extrapolate
                else surface.push_back( 
                            z * dZ +
                            (*dataIt - isoValue)/(*dataIt - *(dataIt+1)) * dZ
                            );
                break;
            }
            --dataIt;
            --z;
        }
        if(dataIt == dataStop) throw types::runtimeError() 
                    << types::errinfo_runtime("Missed an isosurface value");
    }
}

}

