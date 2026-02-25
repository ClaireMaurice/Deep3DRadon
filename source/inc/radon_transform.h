#pragma once

#include "pattern.h"

#include "eigen3/Eigen/Core"

class RadonTransform
{
public:
    RadonTransform();
    ~RadonTransform();  

    void setViewPoint(const Eigen::Vector3d& viewPoint) {
        this->viewPoint = viewPoint;
    }

    void dump() const; // for debugging purposes, to print the parameters of the Radon transform, such as the angles and the distances, etc.

    void apply(const Pattern& pattern);

    private:
        Eigen::Vector3d viewPoint;
        
};
