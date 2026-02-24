#pragma once
#include <iostream>

#include "eigen3/Eigen/Core"
#include "orientation.h"

class SourcePoint
{
public:
    SourcePoint();
    SourcePoint(const Euler& orientation, const Eigen::Vector3d& position):
        m_orientation(orientation),
        m_position(position)
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
    

    void dump() const {
        std::cout << "SourcePoint parameters:" << std::endl;
        std::cout << "Orientation (Euler angles): " << m_orientation << std::endl;
        std::cout << "Position: " << m_position.transpose() << std::endl;
    }; // for debugging purposes, to print the source point parameters



private:
    Euler m_orientation; // the orientation of the source point, which can be used to calculate the diffraction pattern based on the orientation of the crystal and the microscope
    Eigen::Vector3d m_position; // the position of the source point in 3D space, which can be used to calculate the diffraction pattern based on the position of the crystal and the microscope

};
