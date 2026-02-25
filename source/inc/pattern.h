#pragma once
#include <iostream>

#define cimg_use_png
#include "CImg.h"
using namespace cimg_library;

class Microscope;
class Crystal;
class SourcePoint;

class Pattern
{
public:
    Pattern() {};
    ~Pattern() {}; 

    void simulate(Microscope &microscope, Crystal& crystal, SourcePoint& source_point);
    void save(const std::string& filename) {
        pattern_img.save(filename.c_str());
    }

    int getWidth() const {
        return pattern_img.width();
    }
    int getHeight() const {
        return pattern_img.height();
    }
    const CImg<unsigned char>& getPattern() const {
        return pattern_img;
    }

private:
    CImg<unsigned char> pattern_img;

};