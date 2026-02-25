#pragma once

#include "eigen3/Eigen/Dense"

class KLines {

public:
    KLines();
    ~KLines();

    int buildKLines(Eigen::Vector3d n, double s2, Eigen::Vector3d sp, int img_width, int img_height);
    bool isHorizontal() const { return is_horizontal; }

    int getY(double x, double& y1, double& y2) const;
    int getX(double y, double& x1, double& x2) const;

    void dump() const; // for debugging purposes, to print the parameters of the KLines, such as the coefficients of the conic section, the orientation of the line, etc.
    
 private:
    double a,b;
    double A,B,C,D,E,F;
    bool is_horizontal;
};