#pragma once

#include "eigen3/Eigen/Dense"

#include "pattern.h"

class KLines {

public:
    KLines();
    ~KLines();

    int buildKLines(Eigen::Vector3d n, double s2, Eigen::Vector3d sp, int img_width, int img_height);
    bool isHorizontal() const { return is_horizontal; }

    int getY(double x, double& y1, double& y2) const;
    int getX(double y, double& x1, double& x2) const;

    double getIntegral(const Pattern& pattern) const; // compute the integral of the pattern along the KLines, which will be used as input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern
    
    // I want to add a function to draw the KLine of the pattern
    void draw(CImg<unsigned char>& img) const; // draw the KLine on the pattern, which can be used for visualization purposes


    void dump() const; // for debugging purposes, to print the parameters of the KLines, such as the coefficients of the conic section, the orientation of the line, etc.
    
 private:
    double a,b;
    double A,B,C,D,E,F;
    bool is_horizontal;
};