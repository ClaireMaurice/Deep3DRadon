#include "pattern.h"

#include "eigen3/Eigen/Dense"

#include "microscope.h"
#include "crystal.h"
#include "source_point.h"
#include "k_lines.h"

void Pattern::simulate(Microscope& microscope, Crystal& crystal, SourcePoint& source_point)
{
    // TODO : implement the simulation of the diffraction pattern based on the microscope, crystal and source point parameters
    bool verbose = false; // for debugging purposes, we can set this to true to print the parameters of the simulation and the intermediate results, and we can set it to false to disable the debug output for a cleaner output when we are confident that the simulation is working correctly

    std::cout << "Simulating diffraction pattern..." << std::endl;

    // allocate a blank image for the pattern of size obtained from the microscope parameters
    int img_width, img_height;
    microscope.getDetector()->getNumPixels(img_width, img_height);
    
    if (verbose) {
        std::cout << "Image width: " << img_width << std::endl;
        std::cout << "Image height: " << img_height << std::endl;
    }

    pattern_img = CImg<unsigned char>(img_width, img_height, 1,1, 0); // create a blank black image

    Quat q_DM = microscope.getDetector()->getOrientation();
    Quat q_MS = microscope.getOrientation();
    Quat q_SC = source_point.getOrientation().Conjugate(); // we need the conjugate of the source point orientation to get the correct transformation from the crystal frame to the detector frame
    Quat q_DC = q_SC * q_MS * q_DM; // the order of multiplication is important here, we need to apply the transformations in the correct order to get the correct orientation of the crystal with respect to the detector
    
    Eigen::Matrix3d R_DC = q_DC.toRotationMatrix();

    if(verbose) {
        std::cout << "Transport matrix (DC) : " << std::endl << R_DC << std::endl;
    }


    // TODO : verification of the transport matrix, we can print the quaternions and check that they are correct, and we can also check that the resulting orientation of the crystal with respect to the detector is correct by comparing it with the expected orientation based on the microscope and source point parameters
    // // need the transport matrice DC = DM * MS * SC
    // Eigen::Matrix3d DM = microscope.getDetector()->getOrientation().toRotationMatrix();
    // Eigen::Matrix3d MS = microscope.getOrientation().toRotationMatrix();
    // Eigen::Matrix3d SC = source_point.getOrientation().toRotationMatrix().transpose();
    
    // Eigen::Matrix3d DC = DM * MS * SC; // the order of multiplication is important here, we need to apply the transformations in the correct order to get the correct orientation of the crystal with respect to the detector
    // std::cout << "Transport matrix (DC) : " << std::endl << DC << std::endl;


    
    // to calculate the reciprocal lattice vectors, we need the lattice parameters of the crystal, which we can get from the unit cell parameters
    // for simplicity we assume a cubic lattice for now
    double a0 = crystal.getUnitCell().getLatticeParameters()[0];  //in Angstroms

    // source point position
    Eigen::Vector3d PC = source_point.getPosition(); 
    
    // convert to PC coordinates in pixels (TODO : check that the conversion is correct)
    Eigen::Vector3d PCpix(PC(0)*img_width,(1-PC(1))*img_height,PC(2)*img_width); // in pixels, with the origin at the top left corner of the image, and the y axis pointing downwards
    
    const unsigned char white[] = {255,255,255};
//    pattern_img.draw_circle(PCpix(0), PCpix(1), 5, white); // draw the source point on the image for visualization purposes, we can later replace this with a more accurate representation of the source point, but for now this is a simple way to visualize the position of the source point on the screen

    // iterate over the reflectors of the crystal and simulate the diffraction spots based on the source point and microscope parameters
    for (const auto& reflector : crystal.getReflectors()) {
        
        Eigen::Vector3i hkl = reflector.getHKL();

        if (verbose) {
            std::cout << "Simulating Kikuchi Lines for reflector with Miller indices (" << hkl.transpose() << ")" << std::endl;   
        }
        
        // compute the reciprocal lattice vector for the given Miller indices, for simplicty a cubic lattice is assumed for now, but this can be easily extended to other lattice types by using the appropriate formulas for the reciprocal lattice vectors based on the lattice parameters of the crystal
        Eigen::Vector3d g_hkl = 1/a0 * hkl.cast<double>(); // in reciprocal angstroms, since a0 is in angstroms and hkl is dimensionless
        
        // transport to detector frame
        g_hkl = R_DC * g_hkl; // apply the transport matrix to get the orientation of the reciprocal lattice vector with respect to the detector frame

        // compute the interplanar spacing d_hkl = 1 / |g_hkl|, 
        //and then compute the diffraction angle theta using Bragg's law: n*lambda = 2*d_hkl*sin(theta), 
        //where n is the order of the reflection (we can assume n=1 for now), 
        //lambda is the electron wavelength obtained from the microscope parameters, 
        //and d_hkl is the interplanar spacing calculated from the reciprocal lattice vector.
    
        // beware of dimensions : 
        double d_hkl = 1.0 / g_hkl.norm();  // in angstroms, since g_hkl is in reciprocal angstroms
        double lambda = microscope.getElectronWavelength(); // in Angstroms
        double sin_theta = lambda / (2 * d_hkl); // Bragg's law
        double s2 = sin_theta * sin_theta;

        if (verbose) {
            std::cout << "Reciprocal lattice vector g_hkl: " << g_hkl.transpose() << std::endl;
            std::cout << "Interplanar spacing d_hkl: " << d_hkl << " Angstroms" << std::endl;
            std::cout << "Electron wavelength lambda: " << lambda << " Angstroms" << std::endl;
            std::cout << "sin(theta): " << sin_theta << std::endl;
        }

     
        KLines k_lines;
        int result = k_lines.buildKLines(g_hkl, s2, PCpix, img_width, img_height);

        if (verbose) {
            std::cout << "KLines parameters for reflector (" << hkl.transpose() << "): " << std::endl;
            k_lines.dump();
        }

        
        if(result == 0) {
            if (verbose)
                std::cout << "Not visible on the screen, skipping this reflector." << std::endl;
            continue;
        } else {
            if (verbose)
                std::cout << "Visible on the screen, drawing the Kikuchi lines for this reflector." << std::endl;

            // let's do the work:
            if(k_lines.isHorizontal()) {
                // we can iterate over the x coordinates of the image and compute the corresponding y coordinates of the conic section defined by the k_lines parameters, and then draw a line between the two intersection points to visualize the spot on the screen
                for (int x = 0; x < img_width ; x++) {
                    double y1, y2;
                    if(k_lines.getY(x, y1, y2)) { // get the corresponding y coordinates of the conic section defined by the k_lines parameters for the given x coordinate
                        if(y1 >= 0 && y1 < img_height) {
                            pattern_img.draw_point(x, y1, white); 
                        }
                        if(y2 >= 0 && y2 < img_height) {
                        pattern_img.draw_point(x, y2, white); 
                        }
                   }
                }
            } else {
                // we can iterate over the y coordinates of the image and compute the corresponding x coordinates of the conic section defined by the k_lines parameters, and then draw a line between the two intersection points to visualize the spot on the screen
                for (int y = 0; y < img_height ; y++) {
                    double x1, x2;
                    if(k_lines.getX(y, x1, x2)) {  // get the corresponding x coordinates of the conic section defined by the k_lines parameters for the given y coordinate
                        if(x1 >= 0 && x1 < img_width) {
                            pattern_img.draw_point(x1, y, white); 
                        }
                        if(x2 >= 0 && x2 < img_width) {
                            pattern_img.draw_point(x2, y, white); 
                        }
                    }
                }
            }
        }
    }
}


