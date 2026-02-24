#include "pattern.h"

#include "eigen3/Eigen/Dense"

#include "microscope.h"
#include "crystal.h"
#include "source_point.h"

void Pattern::simulate(Microscope& microscope, Crystal& crystal, SourcePoint& source_point)
{
    // TODO : implement the simulation of the diffraction pattern based on the microscope, crystal and source point parameters
    // This is a placeholder implementation and should be replaced with actual simulation logic
    std::cout << "Simulating diffraction pattern..." << std::endl;

    // allocate a blank image for the pattern of size obtained from the microscope parameters
    int img_width, img_height;
    microscope.getDetector()->getNumPixels(img_width, img_height);
    std::cout << "Creating a blank image of size " << img_width << " x " << img_height << std::endl;

    pattern_img = CImg<unsigned char>(img_width, img_height, 1,1, 0); // create a blank black image


    // need the transport matrice DC = DM * MS * SC
    Eigen::Matrix3d DM = microscope.getDetector()->getOrientation().toRotationMatrix();
    Eigen::Matrix3d MS = microscope.getOrientation().toRotationMatrix();
    Eigen::Matrix3d SC = source_point.getOrientation().toRotationMatrix().transpose();
    
    Eigen::Matrix3d DC = DM * MS * SC; // the order of multiplication is important here, we need to apply the transformations in the correct order to get the correct orientation of the crystal with respect to the detector
    std::cout << "Transport matrix (DC) : " << std::endl << DC << std::endl;

    Quat q_DM = microscope.getDetector()->getOrientation();
    Quat q_MS = microscope.getOrientation();
    Quat q_SC = source_point.getOrientation().Conjugate(); // we need the conjugate of the source point orientation to get the correct transformation from the crystal frame to the detector frame
    Quat q_DC = q_SC * q_MS * q_DM; // the order of multiplication is important here, we need to apply the transformations in the correct order to get the correct orientation of the crystal with respect to the detector
    
    // TODO : verification of the transport matrix, we can print the quaternions and check that they are correct, and we can also check that the resulting orientation of the crystal with respect to the detector is correct by comparing it with the expected orientation based on the microscope and source point parameters
    
    Eigen::Matrix3d R_DC = q_DC.toRotationMatrix();
    std::cout << "Transport matrix (DC) : " << std::endl << R_DC << std::endl;

    // to calculate the reciprocal lattice vectors, we need the lattice parameters of the crystal, which we can get from the unit cell parameters, and then we can calculate the reciprocal lattice vectors using the formula g_hkl = h*b1 + k*b2 + l*b3, where b1, b2 and b3 are the reciprocal lattice vectors corresponding to the lattice parameters a, b and c, and h, k and l are the Miller indices of the reflector
    // for simplicity we assume a cubic lattice for now
    double a0 = crystal.getUnitCell().getLatticeParameters()[0];  //in Angstroms

    // source point position
    Eigen::Vector3d PC = source_point.getPosition(); 
    
    // convert to PC coordinates in pixels
    Eigen::Vector3d PCpix(PC(0)*img_width,(1-PC(1))*img_height,PC(2)*img_width); // in pixels, with the origin at the top left corner of the image, and the y axis pointing downwards
    
    // short hand for the source position coordinates
    double sx = PCpix(0);
    double sy = PCpix(1);
    double sz = PCpix(2);
    
    const unsigned char white[] = {255,255,255};
    pattern_img.draw_circle(sx, sy, 5, white); // draw the source point on the image for visualization purposes, we can later replace this with a more accurate representation of the source point, but for now this is a simple way to visualize the position of the source point on the screen

    // iterate over the reflectors of the crystal and simulate the diffraction spots based on the source point and microscope parameters
    for (const auto& reflector : crystal.getReflectors()) {
        
        Eigen::Vector3i hkl = reflector.getHKL();

        std::cout << "Simulating spot for reflector with Miller indices (" << hkl.transpose() << ")" << std::endl;   
        
        // compute the reciprocal lattice vector for the given Miller indices, 
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

        std::cout << "sin(theta) : " << sin_theta << std::endl;

      
        // need the unit vector of the reciprocal lattice vector in the detector frame
        Eigen::Vector3d n_hkl = g_hkl.normalized();
        if(n_hkl(1) < 0) {
            n_hkl = -n_hkl;     
        }
        // short hand for the conic coefficients
        double nx = n_hkl(0);
        double ny = n_hkl(1);
        double nz = n_hkl(2);

        double F =  (nx*sx - ny*sy - nz*sz);

        bool is_horizontal = fabs(ny) > fabs(nx); // we can determine if the conic is a horizontal or vertical ellipse based on the coefficients of the conic, if the coefficient of y^2 is larger than the coefficient of x^2, then it is a horizontal ellipse, otherwise it is a vertical ellipse
        double a, b; // the coefficients for the plane trace

        if(is_horizontal) {
            a = nx/ny;
            b = -F/ny;
        } else {
            a = ny/nx;
            b = F/nx;
        }

        // check screen intersection
        Eigen::Vector2d TopLeft(0,0);
        Eigen::Vector2d BottomRight(img_width,img_height);

        Eigen::Vector2d P1 = TopLeft;
        Eigen::Vector2d P2 = BottomRight;

        if(is_horizontal) {
            P1(1) = b;
            P2(1) = a*P2(0) + b;

            if(P1(1) < 0 && P2(1) < 0) {
                //std::cout << "Spot is not visible on the screen (below the screen)" << std::endl;
                continue; // skip this spot if it is not visible on the screen
            }
            if(P1(1) > img_height && P2(1) > img_height) {
                //std::cout << "Spot is not visible on the screen (above the screen)" << std::endl;
                continue; // skip this spot if it is not visible on the screen
            }
        

        } else {
            P1(0) = b;
            P2(0) = a*P2(1) + b;

            if(P1(0) < 0 && P2(0) < 0) {
                //std::cout << "Spot is not visible on the screen (below the screen)" << std::endl;
                continue; // skip this spot if it is not visible on the screen
            }
            if(P1(0) > img_height && P2(0) > img_height) {
                //std::cout << "Spot is not visible on the screen (above the screen)" << std::endl;
                continue; // skip this spot if it is not visible on the screen
            }
        }

        
        //pattern_img.draw_line(P1(0),P1(1), P2(0), P2(1), white); // draw a straight line between the two points, we can later replace this with a more accurate drawing of the conic section, but for now this is a simple approximation to visualize the spot on the screen

        // compute the conic coefficients
        double A = nx*nx - s2;
        double B = -nx*ny;
        double C = ny*ny - s2;  
        double D = -A*sx - B*sy + nx*nz*sz;
        double E = -C*sy - B*sx - ny*nz*sz;
        F = F * F - s2 * (sx*sx + sy*sy + sz*sz);

        

        if(is_horizontal) {
            for(int x = 0; x < img_width ; x++) {
                double b = B*x + E;
                double delta = b*b - C*(F+x*(A*x + 2*D));
                if(delta < 0) {
                    continue; // no intersection with this vertical line
                }else{
                    double sqrt_delta = sqrt(delta);
                    double y1 = (-b + sqrt_delta) / C;
                    double y2 = (-b - sqrt_delta) / C;

                    
                    pattern_img.draw_point(x, y1, white); // draw a line between the two intersection points, we can later replace this with a more accurate drawing of the conic section, but for now this is a simple way to visualize the spot on the screen
                    pattern_img.draw_point(x, y2, white);
                    
                }
            }
        } else {
            for(int y = 0; y < img_height ; y++) {
                double b = B*y + D;
                double delta = b*b - A*(F+y*(C*y + 2*E));
                if(delta < 0) {
                    continue; // no intersection with this vertical line
                }else{
                    double sqrt_delta = sqrt(delta);
                    double x1 = (-b + sqrt_delta) / A;
                    double x2 = (-b - sqrt_delta) / A;

                    
                    pattern_img.draw_point(x1, y, white); // draw a line between the two intersection points, we can later replace this with a more accurate drawing of the conic section, but for now this is a simple way to visualize the spot on the screen
                    pattern_img.draw_point(x2, y, white);
                    
                }
            }

        }

         
    }

}


