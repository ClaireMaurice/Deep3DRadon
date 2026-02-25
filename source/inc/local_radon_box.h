#pragma once

#include "pattern.h"

#include "eigen3/Eigen/Core"
#include "CImg.h"
using namespace cimg_library;

class LocalRadonBox
{
public:
    LocalRadonBox(int box_size, double increment) : box_size(box_size), increment(increment) {}  // constructor to initialize the local Radon box with the given size and increment, which will be used to compute the Radon transform of the pattern based on the microscope, crystal and source point parameters, and which will be used as input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern
    ~LocalRadonBox() {} 

    int compute(const Eigen::Vector4d& projectorNormal, const SourcePoint& viewPoint, const Pattern& pattern); // compute the Radon transform of the pattern for this projector normal, which will be used as input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern
    CImg<double> get_boxImage() const {
        return boxImage;
    } // get the local Radon box image, which will be used to compute the Radon transform of the pattern based on the microscope, crystal and source point parameters, and which will be used as input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern         
    

private:
    Eigen::Vector3d center; 
    int box_size; // the size of the local Radon box, which will be used to compute the Radon transform of the pattern based on the microscope, crystal and source point parameters, and which will be used as input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction patterndou
    double increment; // the increment of the local Radon box, which will be used to compute the Radon transform of the pattern based on the microscope, crystal and source point parameters, and which will be used as input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern

    CImg<double> boxImage; // the image of the local Radon box, which will be used to compute the Radon transform of the pattern based on the microscope, crystal and source point parameters, and which will be used as input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern
};