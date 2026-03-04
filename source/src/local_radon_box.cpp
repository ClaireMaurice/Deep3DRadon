#include <iostream>

#include "local_radon_box.h"
#include "k_lines.h"
#include "source_point.h"
#include "pattern.h"

int LocalRadonBox::precompute(Eigen::Vector4d &projectorNormal, const SourcePoint &viewPoint, const Pattern &pattern, const int box_size, const double increment) {
      // this should check that the k-lines are valid, long enough and well equilibrated.
    Eigen::Vector3d g_hkl = projectorNormal.head<3>(); // the first three components of the projector vector are the unit reciprocal lattice vector
    double s = sqrt(projectorNormal(3));

    int img_width = pattern.getWidth();
    int img_height = pattern.getHeight();

    // quick and dirty way to recenter the box around the detected normal vector
    // build a large 2D slice of the local Radon box at the theoretical Bragg angle
    // defines the local frame
    Eigen::Vector3d e3 = g_hkl;

    // e1 is perpendicular to e3 and lies on the screen (e.g. parallel to the trace of the Kikuchi lines on the screen)
    Eigen::Vector3d e1 = Eigen::Vector3d(-e3(1), e3(0), 0).normalized();
    // e2 is perpendicular to both e3 and e1, and completes the right-handed local frame
    Eigen::Vector3d e2 = e3.cross(e1).normalized();


    CImg<double> boxSlice = CImg<double>(box_size, box_size, 1, 2, 0);
    const Eigen::Vector3d viewPos = viewPoint.getPositionInPixels(pattern.getWidth(), pattern.getHeight());

    int nbPlus_max = 0, nbMinus_max = 0;  // initialize the counts of pixels in the plus and minus branches of the KLines, which will be used to check that the KLines are long enough and well equilibrated, and which can be used as additional input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern
    double Iplus_max = 0, Iminus_max = 0; // initialize the maximum intensity values in the plus and minus branches of the KLines, which will be used to check that the KLines are well equilibrated, and which can be used as additional input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern
    for (int i = 0; i < box_size; i++)
    {
        for (int j = 0; j < box_size; j++)
        {
            // Create the normal vector once per (i,j) pair instead of for each k
            double p1 = (i - box_size / 2) * increment; // the p1 coordinate of the normal vector, which we can use to compute the angles for the Radon transform based on the microscope and source point parameters
            double p2 = (j - box_size / 2) * increment; // the p2 coordinate of the normal vector, which we can use to compute the angles for the Radon transform based on the microscope and source point parameters
            Eigen::Vector3d n = p1 * e1 + p2 * e2 + e3; // the normal vector of the plane defined by (p1, p2), which we can use to compute the angles for the Radon transform based on the microscope and source point parameters
            n.normalize();                              // normalize the normal vector to get the unit normal vector

            int nbPlus = 0, nbMinus = 0; // initialize the counts of pixels in the plus and minus branches of the KLines for this (i,j) pair, which will be used to check that the KLines are long enough and well equilibrated, and which can be used as additional input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern
            double Iplus = 0, Iminus = 0;
            KLines k_lines;
            if (k_lines.buildKLines(n, s * s, viewPos, img_width, img_height))
            {
                k_lines.getIntegral(pattern, &Iplus, &Iminus, &nbPlus, &nbMinus); // compute the integral of the pattern along the KLines, which will be used as input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern; we can also compute the integrals along the plus and minus branches of the KLines, which can be used as additional input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern
                boxSlice(i, j, 0) = Iplus;
                boxSlice(i, j, 1) = Iminus;
                if (Iplus > Iplus_max)
                {
                    Iplus_max = Iplus;   // update the maximum intensity value in the plus branch of the KLines, which will be used to check that the KLines are well equilibrated, and which can be used as additional input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern
                    nbPlus_max = nbPlus; // update the count of pixels in the plus branch of the KLines, which will be used to check that the KLines are long enough and well equilibrated, and which can be used as additional input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern
                }
                if (Iminus > Iminus_max)
                {
                    Iminus_max = Iminus;   // update the maximum intensity value in the minus branch of the KLines, which will be used to check that the KLines are well equilibrated, and which can be used as additional input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern
                    nbMinus_max = nbMinus; // update the count of pixels in the minus branch of the KLines, which will be used to check that the KLines are long enough and well equilibrated, and which can be used as additional input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern
                }
            }
            else
            {
                boxSlice(i, j, 0) = 0;
                boxSlice(i, j, 1) = 0;
            }
        }
    }

    // check that the Klines are long enough (i.e. have a minimum number of pixels with non-zero intensity - should define a threshold for this in the future)
    // if ((nbPlus_max + nbMinus_max) / 2 < 800) // check that the two branches have similar length (i.e. mas intensity ratio close to 1 - should define a threshold for this in the future)
    // {
    //     std::cout << "KLines are too short, skipping this box. nbPlus_max: " << nbPlus_max << ", nbMinus_max: " << nbMinus_max << std::endl;
    //     return 0; // return 0 to indicate that the precomputation was not successful
    // }

    double intensity_ratio = boxSlice.get_channel(0).max() / (boxSlice.get_channel(1).max() + 1e-6); // add a small value to avoid division by zero
    std::cout << "Intensity ratio " << intensity_ratio << std::endl;                                 // print the intensity ratio of the plus and minus branches of the KLines, which can be used as a quick check to see if the KLines are valid and well equilibrated, and which can be used as additional input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern

    if (intensity_ratio < 0.5 || intensity_ratio > 2.0) // if the intensity ratio is too far from 1, we can consider that the KLines are not valid or not well equilibrated, and we can skip this box
    {
        std::cout << "KLines are not well equilibrated, skipping this box." << std::endl;
        return 0; // return 0 to indicate that the precomputation was not successful
    }

    // I want the pixel coordinates of the maximum intensity in the plus and minus branches of the KLines, which can be used to recenter the smaller box
    double max_val = boxSlice.get_channel(0).max();
    int xmax = 0, ymax = 0;
    cimg_forXY(boxSlice, x, y)
    {
        if (boxSlice(x, y, 0) == max_val)
        {
            xmax += x;
            ymax += y;
            break;
        }
    }

    max_val = boxSlice.get_channel(1).max();

    cimg_forXY(boxSlice, x, y)
    {
        if (boxSlice(x, y, 1) == max_val)
        {
            xmax += x;
            ymax += y;
            break;
        }
    }
    

    //update the projector normal to recenter the box around the detected normal vector
    double p1 = (xmax/2.0 - box_size / 2) * increment; // the p1 coordinate of the normal vector, which we can use to compute the angles for the Radon transform based on the microscope and source point parameters
    double p2 = (ymax/2.0 - box_size / 2) * increment; // the p2 coordinate of the normal vector, which we can use to compute the angles for the Radon transform based on the microscope and source point parameters
    
    std::cout << "Recentered p1: " << p1 << ", p2: " << p2 << std::endl; // print the new p1 and p2 coordinates of the normal vector, which can be used to compute the angles for the Radon transform based on the microscope and source point parameters, and which can be used as additional input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern
    
    Eigen::Vector3d n = p1 * e1 + p2 * e2 + e3; // the normal vector of the plane defined by (p1, p2), which we can use to compute the angles for the Radon transform based on the microscope and source point parameters
    n.normalize();                              // normalize the normal vector to get the unit normal vector
    std::cout << "Recentered normal vector: " << n.transpose() << std::endl; // print the new normal vector, which can be used to compute the angles for the Radon transform based on the microscope and source point parameters, and which can be used as additional input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern
    projectorNormal.head<3>() = n;              // update the first three components of the projector normal vector with the new normal vector, which will be used to compute the angles for the Rad

    return 1; // return 1 to indicate that the precomputation was successful
}

int LocalRadonBox::compute(const Eigen::Vector4d &projectorVector, const SourcePoint &viewPoint, const Pattern &pattern)
{
    Eigen::Vector3d g_hkl = projectorVector.head<3>(); // the first three components of the projector vector are the unit reciprocal lattice vector
    double s = sqrt(projectorVector(3));               // the fourth component of the projector vector is the sin^2(theta) value

    // defines the local frame
    Eigen::Vector3d e3 = g_hkl;

    // e1 is perpendicular to e3 and lies on the screen (e.g. parallel to the trace of the Kikuchi lines on the screen)
    Eigen::Vector3d e1 = Eigen::Vector3d(-e3(1), e3(0), 0).normalized();
    // e2 is perpendicular to both e3 and e1, and completes the right-handed local frame
    Eigen::Vector3d e2 = e3.cross(e1).normalized();

    // set up the local grid for the local Radon box
    CImg<double> basalPlane(box_size, box_size, 1, 3, 0); // the basal plane contains a set of normal vectors - the center is e3 and n = p1*e1 + p2*e2 + e3
    for (int i = 0; i < box_size; i++)
    {
        for (int j = 0; j < box_size; j++)
        {
            double p1 = (i - box_size / 2) * increment; // the p1 coordinate of the normal vector, which we can use to compute the angles for the Radon transform based on the microscope and source point parameters
            double p2 = (j - box_size / 2) * increment; // the p2 coordinate of the normal vector, which we can use to compute the angles for the Radon transform based on the microscope and source point parameters
            Eigen::Vector3d n = p1 * e1 + p2 * e2 + e3; // the normal vector of the plane defined by (p1, p2), which we can use to compute the angles for the Radon transform based on the microscope and source point parameters
            n.normalize();                              // normalize the normal vector to get the unit normal vector

            basalPlane(i, j, 0) = n(0);
            basalPlane(i, j, 1) = n(1);
            basalPlane(i, j, 2) = n(2);
        }
    }
    // the z_axis contains the brag angle, centered on s
    CImg<double> z_axis(box_size, 1, 1, 1, 0);
    for (int i = 0; i < box_size; i++)
    {
        double p3 = (i - box_size / 2) * increment; // the p3 coordinate of the normal vector, which we can use to compute the angles for the Radon transform based on the microscope and source point parameters
        z_axis(i) = s + p3;
    }

    // let's compute the local Radon box image, which will be used as input for the neural network
    boxImage = CImg<double>(box_size, box_size, box_size, 1, 0);

    // Cache values that don't change in the loops to avoid repeated function calls
    const int img_width = pattern.getWidth();
    const int img_height = pattern.getHeight();
    const Eigen::Vector3d viewPos = viewPoint.getPositionInPixels(img_width, img_height);
    
    for (int i = 0; i < box_size; i++)
    {
        for (int j = 0; j < box_size; j++)
        {
            // Read the normal vector from the basal plane for this (i,j) pair
            Eigen::Vector3d n(basalPlane(i, j, 0), basalPlane(i, j, 1), basalPlane(i, j, 2));

            for (int k = 0; k < box_size; k++)
            {
                double s = z_axis(k);

                KLines k_lines;
                if (k_lines.buildKLines(n, s * s, viewPos, img_width, img_height))
                {
                    boxImage(i, j, k) = k_lines.getIntegral(pattern);
                }
                else
                {
                    boxImage(i, j, k) = 0;
                }
            }
        }
    }
    return 1; // return 1 to indicate that the computation was successful
}