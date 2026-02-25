#include <iostream>

#define cimg_use_tiff
#define cimg_use_png
#include "CImg.h"

#include "orientation.h"


int dev_test() {
    std::cout << "Hello, World!" << std::endl;

    // lets's try to use the orientation class
    Miller m1; // default constructor
    std::cout << "Miller m1 : " << m1 << std::endl;

    // can I use Eigen vectors ?
    Eigen::Vector3i hkl(1,2,3);
    Eigen::Vector3i uvw(4,5,6);
    Miller m2(hkl, uvw);
    std::cout << "Miller m2 : " << m2 << std::endl;

  
    // let's try to use CImg library to generate an image
    cimg_library::CImg<unsigned char> image(256,256,1,3); // replace with your image path
    // fill the image with a gradient
    cimg_forXY(image,x,y) {
        image(x,y,0) = x; // red channel
        image(x,y,1) = y; // green channel
        image(x,y,2) = 128; // blue channel
    }
    // draw a white circle in the center
    const unsigned char white[] = {255,255,255};
    image.draw_circle(128,128,50, white);
    // save the image
    image.save("gradient.png");
    image.save_tiff("gradient.tiff");
    // works ok on 2D images

    // let's try to use CImg library to generate a 3D image
    cimg_library::CImg<double> image3D(16,16,16,3); // replace with your image path
    // fill the image with a gradient
    cimg_forXYZ(image3D,x,y,z) {
        image3D(x,y,z,0) = x; // red channel
        image3D(x,y,z,1) = y; // green channel
        image3D(x,y,z,2) = z; // blue channel
    }
    // save the image
    image3D.save("gradient3D.tiff");








    return 0;

}