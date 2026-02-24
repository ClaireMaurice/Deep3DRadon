#pragma once

#include <type_traits>
#include <vector>   
#include <string>
#include <iostream>
#include <eigen3/Eigen/Dense>

class Reflector {
public:
    Reflector() : m_hkl(0,0,0) {
        // default constructor
    }
    ~Reflector() {
        // nothing to do for now, but if we had allocated resources we should release them here
    }
    
    Reflector(const Eigen::Vector3i& hkl) : m_hkl(hkl) {        // constructor with parameters
    }

    void setHKL(const Eigen::Vector3i& hkl) {
        m_hkl = hkl;
    }
    void getHKL(Eigen::Vector3i& hkl) const {
        hkl = m_hkl;
    }

    Eigen::Vector3i getHKL() const {
        return m_hkl;
    }
    

    private:
    // we will need to store the unit cell parameters and the atomic positions in order to calculate the structure factor and the diffraction pattern, but for now we will just store a reference to the unit cell and the list of atoms in the unit cell, and we can extend this class later to include
    // more complex functionality such as handling multiple unit cells, handling different types of atoms, etc.
    //UnitCell m_unitCell;
    Eigen::Vector3i m_hkl; // this will store the calculated diffraction pattern, which we can then use to visualize the results or to compare with experimental data    



};