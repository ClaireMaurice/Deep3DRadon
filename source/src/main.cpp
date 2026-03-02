#include <iostream>
#include <string>
// let's try to use CImg library to load and display an image

#define cimg_use_png
#define cimg_use_tiff
#include "CImg.h"

#include "orientation.h"
#include "microscope.h"
#include "detector.h"
#include "crystal.h"
#include "source_point.h"
#include "pattern.h"
#include "radon_transform.h"

extern int dev_test();

#define DEV_TEST 0

int main()
{
    if (DEV_TEST)
    {
        std::cout << "Running dev test..." << std::endl;
        dev_test();
        return 0;
    }
    else
    {
        std::cout << "Running main program..." << std::endl;
    }
    // dev_test();

    Microscope microscope;
    // microscope.dump();

    Crystal crystal;

    crystal.buildUnitCell("Copper", "FCC", 3.615);
    crystal.buildReflectors();
    // crystal.dump();

    // quick and dirty solution to apply a distortion to the unit cell parameters to test the effect of unit cell deformation on the diffraction pattern and the Radon transform, which will be used for training the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern
    // this is the F tensor - should be symmetric if no rotation... to be further developed properly to deal with continuum mechanics stress and strain tensors
    // And remember that elastic strains are of the order of 10⁻3

    // For no deformation, simply use the Identity matrix as the deformation gradient, which will not change the lattice parameters of the unit cell 
    Eigen::Matrix3d deformationGradient = Eigen::Matrix3d::Identity();

    deformationGradient << 1.005, 0.0, 0.0,
                           0.0, 1.0, 0.0,
                           0.0, 0.0, 1.0; 

    
    Eigen::Matrix3d latticeMatrix = crystal.getUnitCell().getLatticeMatrix();
    Eigen::Matrix3d deformedLatticeMatrix = deformationGradient * latticeMatrix; // apply the deformation gradient to the lattice matrix to get the deformed lattice matrix
    UnitCell deformedUnitCell = UnitCell::getFromLatticeMatrix(deformedLatticeMatrix); // create a new unit cell based on the deformed lattice matrix, which will be used to calculate the diffraction pattern and the Radon transform based on the deformed unit cell parameters, which will be used for training the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern
    
    SourcePoint sourcePoint(Euler(10,0,0), Eigen::Vector3d(0.5, 0.7, 0.5), deformedUnitCell);
    
    Pattern pattern;
    pattern.simulate(microscope, crystal, sourcePoint);
    pattern.save("simulated_pattern.png");

    // Here we use the initial undeformed unit cell for the view point, which is what we would have in a real experiment where we don't know the exact unit cell parameters of the crystal
    SourcePoint viewPoint(Euler(10, 0, 0), Eigen::Vector3d(0.5, 0.7, 0.5), crystal.getUnitCell());

    RadonTransform radonTransform;
    radonTransform.compute(microscope, crystal, viewPoint, pattern);

    radonTransform.save("radon_transform_(10,0,0)_(5,7,5).png", crystal); // save the Radon transform of the pattern to a file with a name that includes the parameters of the view point, which can be used for visualization or for input to the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern

    return 0;
}
