#pragma once

#include "pattern.h"
#include "source_point.h"
#include "eigen3/Eigen/Core"
#include "CImg.h"

using namespace cimg_library;

class LocalRadonBox
{
public:
    LocalRadonBox(int box_size, double increment) : box_size(box_size), increment(increment) {}
    ~LocalRadonBox() {}

    // Precompute pour trouver le centre optimal (le recentrage de Claire)
    int precompute(Eigen::Vector4d &projectorNormal, const SourcePoint &viewPoint, const Pattern &pattern, const int box_size, const double increment);

    // Compute pour générer le volume 3D final
    int compute(const Eigen::Vector4d &projectorNormal, const SourcePoint &viewPoint, const Pattern &pattern);
    
    // Getters pour l'IA et le CSV
    CImg<double> get_boxImage() const { return boxImage; }
    Eigen::Vector4d get_finalProjectorNormal() const { return m_finalProjectorNormal; }
    double get_increment() const { return increment; } // Indispensable pour l'échelle de l'IA

private:
    Eigen::Vector3d center;
    int box_size;
    double increment;
    CImg<double> boxImage;
    Eigen::Vector4d m_finalProjectorNormal; // Stocke la normale après recentrage
};