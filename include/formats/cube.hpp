/*
 * cube.h
 */
#ifndef FORMATS_CUBE_H
#define FORMATS_CUBE_H

#include "atomistic/fundamental.hpp"
#include "types.hpp"


namespace formats {

struct Cube {
	std::vector<atomistic::Atom> atoms;
	la::Grid grid;
    types::String title;
    types::String description;
    types::String fileName;

    Cube(types::String filename) { readCubeFile(filename); }
    Cube() {};

    const Cube & operator+=(const Cube &c);

	bool readCubeFile(types::String filename);
	bool readCubeFile();
    bool readDescription(types::String filename);
	bool writeCubeFile(types::String filename) const;
	bool writeZProfile(types::String filename, types::String header) const;
	bool writeZProfile(types::String filename) const;

    void squareValues() { grid.squareValues(); }
    void stride(std::vector<types::Uint> s) { grid.stride(s);}
    void averageXY(std::vector<types::Real> &data) const {grid.averageXY(data);}
    std::vector<types::Real> getZProfile() const;
    void addZProfile(types::Stream &stream, types::String header) const;
    void getZPlane(types::Uint zIndex, std::vector<types::Real> &data){ 
        grid.zPlane(zIndex, data);}

    void print() const;
	void addHeader(types::Stream &stream) const;
	void addData(types::Stream &stream) const;

    types::Uint countAtoms() const;
	types::Uint countPoints() const;
};

/**
 * Cube file of a wave function.
 * CP2K will write the level and the spin in the title
 */
struct WfnCube : public Cube {
    types::Uint spin;       /**< as read from CP2K cube file */
    types::Uint wfn;        /**< as read from CP2K cube file */
    types::Real energy;     /**< [a.u.], read from CP2K output */
   
    bool readCubeFile(types::String filename);
    bool readCubeFile();
    bool readDescription(types::String filename);
};

/**
 * Specialized grid with several member functions for 3d
 */
struct CubeGrid : la::Grid {

    void sumXY(std::vector<types::Real>& reduced) const;
    void zPlane(types::Uint index, std::vector<types::Real> &plane);
    void averageXY(std::vector<types::Real>& reduced) const;
    /**
     * In lack of a proper resampling method. No new values are calculated
     */
    void stride(types::Uint, types::Uint, types::Uint);
    void resize(types::Uint nX, types::Uint nY, types::Uint nZ);

    types::Real getNearestDataPoint(types::Real x, types::Real y, types::Real z) const;
    types::Real getDataPoint(types::Uint x, types::Uint y, types::Uint z) const;

private:
};




}

#endif
