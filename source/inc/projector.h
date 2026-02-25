#pragma once

#include "eigen3/Eigen/Core"

#include "microscope.h"
#include "crystal.h"
#include "source_point.h"

class Projector
{
public:
    Projector() {};
    ~Projector() {};

    void buildProjectorList(Microscope &microscope, Crystal &crystal, SourcePoint &source_point); // this function will build the list of projection normals based on the microscope, crystal and source point parameters, and the resulting projection normals will be stored in the projection_normals vector
    const std::vector<Eigen::Vector4d>& getProjectorList() const { return projector_list; } // this function will return the list of projection normals, which will be used to simulate the diffraction pattern on the detector based on the microscope and source point parameters

private:
    std::vector<Eigen::Vector4d> projector_list; // this vector will store the projection normals for each reflector of the crystal, which will be used to simulate the diffraction pattern on the detector based on the microscope and source point parameters
};