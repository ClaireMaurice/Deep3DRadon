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
    RadonTransform(bool precompute = true, int box_size_precompute = 64, int box_size_final = 16, double increment_precompute = 0.001, double increment_final = 0.0005) : m_radon_precompute(precompute), m_radon_box_size_precompute(box_size_precompute), m_radon_box_size_final(box_size_final), m_radon_increment_precompute(increment_precompute), m_radon_increment_final(increment_final) {} // constructor to initialize the Radon transform with the given parameters, which will be used to compute the Radon transform of the pattern based on the microscope, crystal and source point parameters, and which will be used as input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern;
    ~RadonTransform();  

    void dump() const; // for debugging purposes, to print the parameters of the Radon transform, such as the angles and the distances, etc.

    void save(const std::string &basename, Crystal &crystal);

    void compute(Microscope& microscope, Crystal& crystal, SourcePoint& viewPoint, Pattern& pattern, const std::string& output_folder = ".");

    void save(const std::string& filename);
    
    private:
        std::vector<LocalRadonBox> m_radonBoxes; // the list of local Radon boxes, which will be used to compute the Radon transform of the pattern based on the microscope, crystal and source point parameters, and which will be used as input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern
        std::vector<int> activeReflectors; // the list of indices of the active reflectors for this view point, which will be used to compute the Radon transform of the pattern based on the microscope, crystal and source point parameters, and which can be used as additional input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern      
        bool m_radon_precompute; // whether to precompute the local Radon boxes for this view point, which will be used to compute the Radon transform of the pattern based on the microscope, crystal and source point parameters, and which will be used as input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern
        int m_radon_box_size_precompute; // the size of the local Radon boxes for precomputation, which will be used to compute the Radon transform of the pattern based on the microscope, crystal and source point parameters, and which will be used as input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern
        int m_radon_box_size_final; // the size of the local Radon boxes for final computation, which will be used to compute the Radon transform of the pattern based on the microscope, crystal and source point parameters, and which will be used as input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern
        double m_radon_increment_precompute; // the increment of the local Radon boxes for precomputation, which will be used to compute the Radon transform of the pattern based on the microscope, crystal and source point parameters, and which will be used as input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern
        double m_radon_increment_final; // the increment of the local
    };
