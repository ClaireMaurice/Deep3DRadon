#pragma once
#include "unit_cell.h"
#include "reflector.h"
#include "orientation.h" // pour 'Euler'
#include "microscope.h"  // pour 'Microscope'

class Crystal
{
public:
    Crystal();
    ~Crystal(); 

    void setUnitCell(const UnitCell& unitCell) { m_unitCell = unitCell; }
    const UnitCell& getUnitCell() const { return m_unitCell; }  

    void buildUnitCell(const std::string& element, const std::string& structureType, double latticeParameter);
    void buildReflectors(); // this function will calculate the diffraction pattern based on the unit cell parameters and the atomic positions, and populate the list of reflectors accordingly

    void getReflectors(std::vector<Reflector>& reflectors) const {
        reflectors = m_reflectors;
    }
    std::vector<Reflector> getReflectors() const {
        return m_reflectors;
    }


    // Cherche la section "public:" et ajoute cette ligne :
    Eigen::Vector3d getNormalForVP(int reflectorIndex, const Euler& orient, const Eigen::Vector3d& VP, const Microscope& microscope) const;
    


    void dump() const; // for debugging purposes, print the crystal parameters

    private:
        // TODO : add lattice parameters, unit cell, etc.
        UnitCell m_unitCell;
        std::vector<Reflector> m_reflectors; 
};