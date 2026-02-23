#include <iostream>
// let's try to use CImg library to load and display an image
#define cimg_use_png
#include "CImg.h"

#include "orientation.h"
#include "microscope.h"
#include "detector.h"


extern int dev_test();

void simulateWiredPattern();
    // TODO : simulate a wired pattern for testing the orientation estim


int main() {
    //dev_test();

    Microscope microscope;
    //Crystal crystal;
    
    std::cout << "Microscope initialized with accelerating voltage : " << microscope.getAcceleratingVoltage() << " kV" << std::endl;
    std::cout << "Electron wavelength : " << microscope.getElectronWavelength() << " nm" << std::endl;
    
    microscope.getDetector()->dump();


    return 0;
}

