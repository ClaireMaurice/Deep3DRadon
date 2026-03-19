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

    int nbPlus_max = 0, nbMinus_max = 0;  // initialize the counts of pixels in the plus and minus branches of the KLines
    double Iplus_max = 0, Iminus_max = 0; // initialize the maximum intensity values
    for (int i = 0; i < box_size; i++)
    {
        for (int j = 0; j < box_size; j++)
        {
            // Correction : division par 2.0 pour éviter l'arrondi entier
            double p1 = (i - box_size / 2.0) * increment; 
            double p2 = (j - box_size / 2.0) * increment; 
            Eigen::Vector3d n = p1 * e1 + p2 * e2 + e3; 
            n.normalize();                               

            int nbPlus = 0, nbMinus = 0; 
            double Iplus = 0, Iminus = 0;
            KLines k_lines;
            if (k_lines.buildKLines(n, s * s, viewPos, img_width, img_height))
            {
                k_lines.getIntegral(pattern, &Iplus, &Iminus, &nbPlus, &nbMinus); 
                boxSlice(i, j, 0) = Iplus;
                boxSlice(i, j, 1) = Iminus;
                if (Iplus > Iplus_max)
                {
                    Iplus_max = Iplus;   
                    nbPlus_max = nbPlus; 
                }
                if (Iminus > Iminus_max)
                {
                    Iminus_max = Iminus;   
                    nbMinus_max = nbMinus; 
                }
            }
            else
            {
                boxSlice(i, j, 0) = 0;
                boxSlice(i, j, 1) = 0;
            }
        }
    }

    double intensity_ratio = boxSlice.get_channel(0).max() / (boxSlice.get_channel(1).max() + 1e-6); 
    std::cout << "Intensity ratio " << intensity_ratio << std::endl;                                

    if (intensity_ratio < 0.5 || intensity_ratio > 2.0) 
    {
        std::cout << "KLines are not well equilibrated, skipping this box." << std::endl;
        return 0; 
    }

    // Extraction des coordonnées du maximum (Ta logique originale conservée)
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
    
    // Update the projector normal (Recentrage de Claire)
    // Utilisation de 2.0 partout pour garantir la précision flottante
    double p1 = (xmax/2.0 - box_size / 2.0) * increment; 
    double p2 = (ymax/2.0 - box_size / 2.0) * increment; 
    
    std::cout << "Recentered p1: " << p1 << ", p2: " << p2 << std::endl; 
    
    Eigen::Vector3d n = p1 * e1 + p2 * e2 + e3; 
    n.normalize();                              
    std::cout << "Recentered normal vector: " << n.transpose() << std::endl; 
    projectorNormal.head<3>() = n;              

    return 1; 
}

int LocalRadonBox::compute(const Eigen::Vector4d &projectorVector, const SourcePoint &viewPoint, const Pattern &pattern)
{
    // SAUVEGARDE DE LA NORMALE POUR LE CSV
    m_finalProjectorNormal = projectorVector; 

    Eigen::Vector3d g_hkl = projectorVector.head<3>(); 
    double s_center = sqrt(projectorVector(3));               

    // defines the local frame
    Eigen::Vector3d e3 = g_hkl;
    Eigen::Vector3d e1 = Eigen::Vector3d(-e3(1), e3(0), 0).normalized();
    Eigen::Vector3d e2 = e3.cross(e1).normalized();

    // set up the local grid for the local Radon box
    CImg<double> basalPlane(box_size, box_size, 1, 3, 0); 
    for (int i = 0; i < box_size; i++)
    {
        for (int j = 0; j < box_size; j++)
        {
            double p1 = (i - box_size / 2.0) * increment; 
            double p2 = (j - box_size / 2.0) * increment; 
            Eigen::Vector3d n = (p1 * e1 + p2 * e2 + e3).normalized(); 

            basalPlane(i, j, 0) = n(0);
            basalPlane(i, j, 1) = n(1);
            basalPlane(i, j, 2) = n(2);
        }
    }

    CImg<double> z_axis(box_size, 1, 1, 1, 0);
    for (int i = 0; i < box_size; i++)
    {
        double p3 = (i - box_size / 2.0) * increment; 
        z_axis(i) = s_center + p3;
    }

    boxImage = CImg<double>(box_size, box_size, box_size, 1, 0);

    const int img_width = pattern.getWidth();
    const int img_height = pattern.getHeight();
    const Eigen::Vector3d viewPos = viewPoint.getPositionInPixels(img_width, img_height);
    
    for (int i = 0; i < box_size; i++)
    {
        for (int j = 0; j < box_size; j++)
        {
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
    return 1; 
}