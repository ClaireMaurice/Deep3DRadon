#pragma once
#include <iostream>

#include "eigen3/Eigen/Core"
#include "orientation.h"
#include "unit_cell.h"

class SourcePoint
{
public:
    SourcePoint() {}

    SourcePoint(const Euler& orientation, const Eigen::Vector3d& position):
        m_orientation(orientation),
        m_position(position)
    {}

    SourcePoint(const Euler& orientation, const Eigen::Vector3d& position, const UnitCell& unitCell):
        m_orientation(orientation),
        m_position(position),
        m_unitCell(unitCell)
    {}    
    
    ~SourcePoint() {
        // nothing to do for now, but if we had allocated resources we should release them here
    }
     
    Quat getOrientation() const {
        return Quat(m_orientation);
    }

    Eigen::Vector3d getPosition() const {
        return m_position;
    }

    Eigen::Vector3d getPositionInPixels(int img_width, int img_height) const {
        // convert the position from normalized coordinates (0 to 1) to pixel coordinates based on the image width and height
        return Eigen::Vector3d(m_position(0) * img_width, (1 - m_position(1)) * img_height, m_position(2) * img_width);
    }


    UnitCell getUnitCell() const {
        return m_unitCell;
    }
    

    void dump() const {
        std::cout << "SourcePoint parameters:" << std::endl;
        std::cout << "Orientation (Euler angles): " << m_orientation << std::endl;
        std::cout << "Position: " << m_position.transpose() << std::endl;
    }; // for debugging purposes, to print the source point parameters

private:
    Euler m_orientation; // the orientation of the source point, which can be used to calculate the diffraction pattern based on the orientation of the crystal and the microscope
    Eigen::Vector3d m_position; // the position of the source point in 3D space, which can be used to calculate the diffraction pattern based on the position of the crystal and the microscope
    UnitCell m_unitCell; // the unit cell of the crystal, which can be used to calculate the diffraction pattern based on the lattice parameters and atomic positions of the crystal    

};
