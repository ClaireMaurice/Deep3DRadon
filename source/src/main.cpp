#include <iostream>
#include <string>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
// let's try to use CImg library to load and display an image

#define cimg_use_png
#define cimg_use_tiff
#include "CImg.h"

#include "config_handler.h"
#include "orientation.h"
#include "microscope.h"
#include "detector.h"
#include "crystal.h"
#include "source_point.h"
#include "pattern.h"
#include "radon_transform.h"


int main(int argc, char *argv[])
{
    // Load configuration from JSON file
    std::string configFile = "../config.json";
    if (argc > 1)
        configFile = argv[1]; // Allow specifying a different config file
    else
        std::cout << "No config file specified, using default '../config.json'." << std::endl;


    // Load and parse configuration
    ConfigHandler config;
    if (!config.load(configFile)) {
        std::cerr << "Failed to load configuration file. Aborting." << std::endl;
        return 1;
    }

    // Create output folder
    if (!config.createOutputFolder()){
        std::cerr << "Failed to create output folder. Aborting." << std::endl;
        return 1;
    }

    // Print loaded configuration
    //config.printParameters();

    Microscope microscope;
    microscope.setAcceleratingVoltage(config.getAcceleratingVoltage());
    microscope.setTilt(config.getTiltAngle(), config.getTiltAxis());
    microscope.dump();
    
    Crystal crystal;
    crystal.buildUnitCell(config.getCrystalElement(), config.getCrystalStructure(), config.getLatticeParameter());
    crystal.buildReflectors();
    crystal.dump();

    // quick and dirty solution to apply a distortion to the unit cell parameters to test the effect of unit cell deformation on the diffraction pattern and the Radon transform, which will be used for training the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern
    // this is the F tensor - should be symmetric if no rotation... to be further developed properly to deal with continuum mechanics stress and strain tensors
    // And remember that elastic strains are of the order of 10⁻3

   
    Eigen::Matrix3d latticeMatrix = crystal.getUnitCell().getLatticeMatrix();
    Eigen::Matrix3d deformedLatticeMatrix = config.getDeformationGradient() * latticeMatrix; // apply the deformation gradient to the lattice matrix to get the deformed lattice matrix
    UnitCell deformedUnitCell = UnitCell::getFromLatticeMatrix(deformedLatticeMatrix);       // create a new unit cell based on the deformed lattice matrix, which will be used to calculate the diffraction pattern and the Radon transform based on the deformed unit cell parameters, which will be used for training the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern

    SourcePoint sourcePoint(config.getSourceEuler(), config.getSourcePosition(),deformedUnitCell);
    config.getSourceEuler();

    Pattern pattern;
    pattern.simulate(microscope, crystal, sourcePoint);
    pattern.save(config.getOutputFolder() + "/simulated_pattern.png");

    // Here we use the initial undeformed unit cell for the view point, which is what we would have in a real experiment where we don't know the exact unit cell parameters of the crystal
    SourcePoint viewPoint(config.getViewEuler(), config.getViewPosition(), crystal.getUnitCell());
    
//    double deso = viewPoint.getOrientation().desorientationAngle(sourcePoint.getOrientation(), 195); // calculate the desorientation angle between the view point orientation and the source point orientation, which will be used to determine the visibility of the reflectors on the screen based on the microscope and source point parameters, and which will be used for training the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern
//    std::cout << "Desorientation angle between view point and source point: " << deso << " degrees" << std::endl;

    RadonTransform radonTransform(config.getRadonPrecompute(), config.getRadonBoxSizePrecompute(), config.getRadonBoxSizeFinal(), config.getRadonIncrementPrecompute(), config.getRadonIncrementFinal());
    radonTransform.compute(microscope, crystal, viewPoint, pattern, config.getOutputFolder());

    // Generate basename based on actual parameters
    Euler viewEuler = config.getViewEuler();
    Eigen::Vector3d viewPosition = config.getViewPosition();    

    std::string basename = config.getOutputFolder() + "/radon_transform_(" + std::to_string((int)viewEuler.phi1()) + "," +
                           std::to_string((int)viewEuler.phi()) + "," + std::to_string((int)viewEuler.phi2()) + ")_(" +
                           std::to_string((int)(viewPosition.x() * 10)) + "," + std::to_string((int)(viewPosition.y() * 10)) + "," +
                           std::to_string((int)(viewPosition.z() * 10)) + ")";
    radonTransform.save(basename, crystal); // save the Radon transform of the pattern to a file with a name that includes the parameters of the view point, which can be used for visualization or for input to the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern

    return 0;
}
