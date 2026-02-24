#include <iostream>
#include "detector.h"

Detector::Detector() {
    // default constructor corresponding to my old camera Nordlys II, 
    //m_qDM = Quat(Euler(359.38, 95.25, 0.51)); 

        m_qDM = Quat(Euler(0, 90, 0)); // for now we will assume that the detector is aligned with the microscope frame, but we can change this later to test different orientations of the detector

        setNumPixels(1344, 1024); // default number of pixels
    setPixelSize(23.47); // default pixel size in microns
}

void Detector::setNumPixels(int numPixelsX, int numPixelsY) {
    m_numPixelsX = numPixelsX;
    m_numPixelsY = numPixelsY;
}

void Detector::setPixelSize(double pixelSize) {
    m_pixelSize = pixelSize;
}

void Detector::setOrientation(const Quat& qDM) {
    m_qDM = qDM;
}

void Detector::dump() const {
    std::cout << "Detector parameters:" << std::endl;
    std::cout << "Orientation (quaternion) : " << m_qDM << std::endl;
    std::cout << "Number of pixels : " << m_numPixelsX << " x " << m_numPixelsY << std::endl;
    std::cout << "Pixel size : " << m_pixelSize << " microns" << std::endl;
}

Detector::~Detector() {
    // nothing to do for now, but if we had allocated resources we should release them here
}