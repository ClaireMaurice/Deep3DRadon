#include <iostream>

#include "local_radon_box.h"
#include "k_lines.h"
#include "source_point.h"
#include "pattern.h"

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
    const Eigen::Vector3d viewPos = viewPoint.getPosition();
    const double patWidth = pattern.getWidth();
    const double patHeight = pattern.getHeight();

    for (int i = 0; i < box_size; i++)
    {
        for (int j = 0; j < box_size; j++)
        {
            // Create the normal vector once per (i,j) pair instead of for each k
            Eigen::Vector3d n(basalPlane(i, j, 0), basalPlane(i, j, 1), basalPlane(i, j, 2));

            for (int k = 0; k < box_size; k++)
            {
                double s = z_axis(k);

                KLines k_lines;
                if (k_lines.buildKLines(n, s * s, viewPos, patWidth, patHeight))
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