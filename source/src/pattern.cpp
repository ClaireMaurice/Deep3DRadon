#include "pattern.h"

#include "eigen3/Eigen/Dense"

#include "microscope.h"
#include "crystal.h"
#include "source_point.h"
#include "k_lines.h"
#include "projector.h"


void Pattern::simulate(Microscope& microscope, Crystal& crystal, SourcePoint& source_point)
{
    // TODO : implement the simulation of the diffraction pattern based on the microscope, crystal and source point parameters
    bool verbose = false; // for debugging purposes, we can set this to true to print the parameters of the simulation and the intermediate results, and we can set it to false to disable the debug output for a cleaner output when we are confident that the simulation is working correctly

    std::cout << "Simulating diffraction pattern..." << std::endl;

    Projector projector;
    projector.buildProjectorList(microscope, crystal, source_point);

    // allocate a blank image for the pattern of size obtained from the microscope parameters
    int img_width, img_height;
    microscope.getDetector()->getNumPixels(img_width, img_height);
    
    if (verbose) {
        std::cout << "Image width: " << img_width << std::endl;
        std::cout << "Image height: " << img_height << std::endl;
    }

    pattern_img = CImg<unsigned char>(img_width, img_height, 1,1, 0); // create a blank black image

        const unsigned char white[] = {255,255,255};
//    pattern_img.draw_circle(PCpix(0), PCpix(1), 5, white); // draw the source point on the image for visualization purposes, we can later replace this with a more accurate representation of the source point, but for now this is a simple way to visualize the position of the source point on the screen

    // iterate over the reflectors of the crystal and simulate the diffraction spots based on the source point and microscope parameters
    for (const auto& reflector : projector.getProjectorList()) {
        
        Eigen::Vector3d g_hkl = reflector.head<3>(); // the first three components of the projector normal are the reciprocal lattice vector, which we can use to compute the diffraction spots on the screen based on the microscope and source point parameters
        double s2 = reflector(3); // the fourth component of the projector normal is the sin^2(theta) value, 
     
        KLines k_lines;
        int result = k_lines.buildKLines(g_hkl, s2, source_point.getPosition(), img_width, img_height);
       
        if(result == 0) {
            if (verbose)
                std::cout << "Not visible on the screen, skipping this reflector." << std::endl;
            continue;
        } else {
            if (verbose)
                std::cout << "Visible on the screen, drawing the Kikuchi lines for this reflector." << std::endl;

            k_lines.draw(pattern_img); 
        }
    }
}


