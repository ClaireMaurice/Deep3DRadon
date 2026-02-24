#include <iostream>
#include "crystal.h"

Crystal::Crystal() : m_unitCell(std::vector<double>{1,1,1,90,90,90}, std::vector<Atom>{}) {
    // default constructor corresponding to a simple cubic unit cell with one atom at the origin
}

Crystal::~Crystal() {
    // nothing to do for now, but if we had allocated resources we should release them here
}

void Crystal::buildUnitCell(const std::string& element, const std::string& structureType, double latticeParameter) {
    // for now we will just handle a few simple cases, but this function can be extended to handle more complex structures
    if (element == "Copper") {
        // for simplicity, we will just use atomic number 29 for all elements, but this can be extended to use the actual atomic number of the element (e.g. 29 for Copper when we have a more complete database of elements and want to do dynamic simulations)
        // we will also ignore the actual atomic positions for now, and just use the ideal positions for the given structure type
        // TODO : we should also handle the case where the lattice parameters are not all equal, and where the angles are not all 90 degrees, but for now we will just assume a cubic structure for simplicity
        
        if (structureType == "FCC") {
            // Face-Centered Cubic structure
            std::vector<double> latticeParameters = {latticeParameter, latticeParameter, latticeParameter, 90, 90, 90};
            std::vector<Atom> atoms = {
                Atom(29, 0.0, 0.0, 0.0), // atom at the origin
                Atom(29, 0.5, 0.5, 0.0), // atom at the face center
                Atom(29, 0.5, 0.0, 0.5), // atom at the face center
                Atom(29, 0.0, 0.5, 0.5)  // atom at the face center
            };
            m_unitCell = UnitCell(latticeParameters, atoms);
        } else if (structureType == "BCC") {
            // Body-Centered Cubic structure
            std::vector<double> latticeParameters = {latticeParameter, latticeParameter, latticeParameter, 90, 90, 90};
            std::vector<Atom> atoms = {
                Atom(1, 0.0, 0.0, 0.0), // atom at the origin
                Atom(1, 0.5, 0.5, 0.5)  // atom at the body center
            };
            m_unitCell = UnitCell(latticeParameters, atoms);
        } else if (structureType == "SimpleCubic") {
            // Simple Cubic structure
            std::vector<double> latticeParameters = {latticeParameter, latticeParameter, latticeParameter, 90, 90, 90};
            std::vector<Atom> atoms = {
                Atom(1, 0.0, 0.0, 0.0) // atom at the origin
            };
            m_unitCell = UnitCell(latticeParameters, atoms);
        } else {
            std::cerr << "Unknown structure type: " << structureType << std::endl;
        }
    } else {
        std::cerr << "Unknown element: " << element << std::endl;
    }
}

void Crystal::buildReflectors() {
    // this function populates the list of reflectors based on the unit cell parameters 
    // for now we will just create a few dummy reflectors for testing purposes, but this function can be extended to calculate the actual diffraction pattern based on the unit cell parameters and the atomic positions
    m_reflectors.clear(); // clear any existing reflectors
    m_reflectors.push_back(Reflector(Eigen::Vector3i(1,1,1))); 
    m_reflectors.push_back(Reflector(Eigen::Vector3i(-1,1,1))); 
    m_reflectors.push_back(Reflector(Eigen::Vector3i(1,-1,1))); 
    m_reflectors.push_back(Reflector(Eigen::Vector3i(1,1,-1)));   
    m_reflectors.push_back(Reflector(Eigen::Vector3i(2,0,0))); 
    m_reflectors.push_back(Reflector(Eigen::Vector3i(0,2,0))); 
    m_reflectors.push_back(Reflector(Eigen::Vector3i(0,0,2)));  
    m_reflectors.push_back(Reflector(Eigen::Vector3i(2,2,0)));
    m_reflectors.push_back(Reflector(Eigen::Vector3i(-2,2,0)));
    m_reflectors.push_back(Reflector(Eigen::Vector3i(2,0,2)));
    m_reflectors.push_back(Reflector(Eigen::Vector3i(-2,0,2)));
    m_reflectors.push_back(Reflector(Eigen::Vector3i(0,2,2)));      
    m_reflectors.push_back(Reflector(Eigen::Vector3i(0,-2,2)));  
    m_reflectors.push_back(Reflector(Eigen::Vector3i(1,1,3)));
    m_reflectors.push_back(Reflector(Eigen::Vector3i(-1,1,3)));
    m_reflectors.push_back(Reflector(Eigen::Vector3i(1,-1,3)));
    m_reflectors.push_back(Reflector(Eigen::Vector3i(1,1,-3)));
    m_reflectors.push_back(Reflector(Eigen::Vector3i(3,1,1)));
    m_reflectors.push_back(Reflector(Eigen::Vector3i(-3,1,1)));
    m_reflectors.push_back(Reflector(Eigen::Vector3i(3,-1,1)));
    m_reflectors.push_back(Reflector(Eigen::Vector3i(3,1,-1))); 
    m_reflectors.push_back(Reflector(Eigen::Vector3i(1,3,1)));
    m_reflectors.push_back(Reflector(Eigen::Vector3i(-1,3,1)));
    m_reflectors.push_back(Reflector(Eigen::Vector3i(1,-3,1)));
    m_reflectors.push_back(Reflector(Eigen::Vector3i(1,3,-1)));
}

void Crystal::dump() const {
    // for debugging purposes, print the crystal parameters
    std::cout << "Unit Cell Parameters: ";
    for (const auto& param : m_unitCell.getLatticeParameters()) {
        std::cout << param << " ";
    }
    std::cout << std::endl;

    std::cout << "Atoms in Unit Cell: " << std::endl;
    for (const auto& atom : m_unitCell.getAtoms()) {
        std::cout << "Atomic Number: " << atom.getAtomicNumber() 
                  << ", Position: (" << atom.getX() << ", " << atom.getY() << ", " << atom.getZ() << ")" 
                  << std::endl;
    }

    std::cout << "Reflectors: " << std::endl;
    for (const auto& reflector : m_reflectors) {
        Eigen::Vector3i hkl;
        reflector.getHKL(hkl);
        std::cout << "HKL: (" << hkl(0) << ", " << hkl(1) << ", " << hkl(2) << ")" << std::endl;
    }
}




