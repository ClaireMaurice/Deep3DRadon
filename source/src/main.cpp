#include <iostream>
// let's try to use CImg library to load and display an image
#define cimg_use_png
#include "CImg.h"

#include "orientation.h"
#include "microscope.h"
#include "detector.h"
#include "crystal.h"
#include "source_point.h"
#include "pattern.h"
#include "radon_transform.h"

extern int dev_test();

int main() {
    //dev_test();

    Microscope microscope;
    //microscope.dump();

    Crystal crystal;

    crystal.buildUnitCell("Copper","FCC",3.615);
    crystal.buildReflectors();
    //crystal.dump();

    SourcePoint sourcePoint(Euler(10,0,0), Eigen::Vector3d(0.5,0.7,0.5));
    //sourcePoint.dump();

    Pattern pattern;
    pattern.simulate(microscope, crystal, sourcePoint);
    pattern.save("simulated_pattern.png");


    // SourcePoint viewPoint(Euler(10,0,0), Eigen::Vector3d(0.5,0.7,0.5));

    // RadonTransform radonTransform;
    // radonTransform.computeRadonTransform(microscope, crystal, viewPoint, pattern);


    return 0;
}

