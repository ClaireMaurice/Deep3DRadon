#pragma once

#include "pattern.h"
#include "microscope.h"
#include "crystal.h"
#include "source_point.h"
#include "local_radon_box.h"


#include "eigen3/Eigen/Core"



class RadonTransform
{
public:
    RadonTransform();
    ~RadonTransform();  

    void dump() const; // for debugging purposes, to print the parameters of the Radon transform, such as the angles and the distances, etc.

    void save(const std::string &basename, Crystal &crystal);

    void compute(Microscope& microscope, Crystal& crystal, SourcePoint& viewPoint, Pattern& pattern);

    void save(const std::string& filename);
    
    private:
        std::vector<LocalRadonBox> m_radonBoxes; // the list of local Radon boxes, which will be used to compute the Radon transform of the pattern based on the microscope, crystal and source point parameters, and which will be used as input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern
        
};
