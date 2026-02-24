#pragma once

#include "orientation.h"


class Detector
{
public:
    // Default constructor
    Detector();
    ~Detector();

    void setNumPixels(int numPixelsX, int numPixelsY);
    void setPixelSize(double pixelSize);
    void setOrientation(const Quat& qDM);

    void getNumPixels(int& numPixelsX, int& numPixelsY) const {
        numPixelsX = m_numPixelsX;
        numPixelsY = m_numPixelsY;
    }
    void getPixelSize(double& pixelSize) const {
        pixelSize = m_pixelSize;
    }

    Quat getOrientation() const { return m_qDM; }


    void dump() const; // for debugging purposes, print the detector parameters

private:
    Quat m_qDM; // the orientation of the detector frame with respect to the microscope frame, expressed as a quaternion
    
    int m_numPixelsX; // number of pixels in the X direction
    int m_numPixelsY; // number of pixels in the Y direction
    double m_pixelSize; // in nm, size of the pixels
};