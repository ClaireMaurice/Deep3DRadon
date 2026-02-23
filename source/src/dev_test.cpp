#include "orientation.h"
#include <iostream>

#define cimg_use_png
#include "CImg.h"

int dev_test() {
    std::cout << "Hello, World!" << std::endl;

    // let's try to use CImg library to load and display an image

    cimg_library::CImg<unsigned char> image("../source/data/Screenshot_20260223_111115.png"); // replace with your image path


    // lets's try to use the orientation class
    Miller m1; // default constructor
    std::cout << "Miller m1 : " << m1 << std::endl;

    // can I use Eigen vectors ?
    Eigen::Vector3i hkl(1,2,3);
    Eigen::Vector3i uvw(4,5,6);
    Miller m2(hkl, uvw);
    std::cout << "Miller m2 : " << m2 << std::endl;
}