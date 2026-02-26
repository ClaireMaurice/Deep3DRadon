#include <iostream>
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

int main()
{
    if (0)
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

    SourcePoint sourcePoint(Euler(65, 76, 24), Eigen::Vector3d(0.5, 0.7, 0.5));
    // sourcePoint.dump();

    Pattern pattern;
    pattern.simulate(microscope, crystal, sourcePoint);
    pattern.save("simulated_pattern.png");

    SourcePoint viewPoint(Euler(10, 0, 0), Eigen::Vector3d(0.5, 0.7, 0.5));

    RadonTransform radonTransform;
    radonTransform.compute(microscope, crystal, viewPoint, pattern);

    // // generate a string representation of the viewPoint parameters, which will be used to create the filename for saving the Radon transform of the pattern, and which will be used for visualization or for input to the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern
    radonTransform.save("radon_transform_(10,0,0)_(5,7,5)", crystal); // save the Radon transform of the pattern to a file with a name that includes the parameters of the view point, which can be used for visualization or for input to the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern

    return 0;
}
