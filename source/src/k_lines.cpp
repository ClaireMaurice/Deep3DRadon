#include <iostream>

#include "k_lines.h"

KLines::KLines() {};

KLines::~KLines() {};

void KLines::dump() const
{
    std::cout << "KLines parameters: " << std::endl;
    std::cout << "a: " << a << std::endl;
    std::cout << "b: " << b << std::endl;
    std::cout << "A: " << A << std::endl;
    std::cout << "B: " << B << std::endl;
    std::cout << "C: " << C << std::endl;
    std::cout << "D: " << D << std::endl;
    std::cout << "E: " << E << std::endl;
    std::cout << "F: " << F << std::endl;
    std::cout << "is_horizontal: " << is_horizontal << std::endl;
}

int KLines::buildKLines(Eigen::Vector3d n, double s2, Eigen::Vector3d sp, int img_width, int img_height)
{

    bool verbose = false; // set to true to enable verbose output for debugging purposes

    if (verbose)
    {
        std::cout << "Building KLines for normal vector n: " << n.transpose() << std::endl;
        std::cout << "s2: " << s2 << std::endl;
        std::cout << "Source point sp: " << sp.transpose() << std::endl;
        std::cout << "Image width: " << img_width << std::endl;
        std::cout << "Image height: " << img_height << std::endl;
    }
    
    // convert to PC coordinates in pixels (TODO : check that the conversion is correct)
    Eigen::Vector3d PCpix(sp(0)*img_width,(1-sp(1))*img_height,sp(2)*img_width); // in pixels, with the origin at the top left corner of the image, and the y axis pointing downwards
 


    n = n.normalized();
    if (n(1) < 0)
    {
        n = -n;
    }

    // short hand notation
    double nx = n(0);
    double ny = n(1);
    double nz = n(2);
    double sx = PCpix(0);
    double sy = PCpix(1);
    double sz = PCpix(2);

    F = (nx * sx - ny * sy - nz * sz); // precompute the F term for later use in the conic coefficients

    is_horizontal = fabs(ny) > fabs(nx); // we can determine if the conic is a horizontal or vertical

    if (is_horizontal)
    {
        a = nx / ny;
        b = -F / ny;
    }
    else
    {
        a = ny / nx;
        b = F / nx;
    }

    // check screen intersection
    Eigen::Vector2d TopLeft(0, 0);
    Eigen::Vector2d BottomRight(img_width, img_height);

    Eigen::Vector2d P1 = TopLeft;
    Eigen::Vector2d P2 = BottomRight;

    if (is_horizontal)
    {
        P1(1) = b;
        P2(1) = a * P2(0) + b;

        if (P1(1) < 0 && P2(1) < 0)
        {
            return 0; // skip this reflector if it is not visible on the screen
        }
        if (P1(1) > img_height && P2(1) > img_height)
        {
            // std::cout << "Spot is not visible on the screen (above the screen)" << std::endl;
            return 0; // skip this reflector if it is not visible on the screen
        }
    }
    else
    {
        P1(0) = b;
        P2(0) = a * P2(1) + b;

        if (P1(0) < 0 && P2(0) < 0)
        {
            return 0; // skip this reflector if it is not visible on the screen
        }
        if (P1(0) > img_width && P2(0) > img_width)
        {
            return 0; // skip this reflector if it is not visible on the screen
        }
    }

    if (verbose)
    {
        std::cout << "Screen intersection points: " << std::endl;
        std::cout << "P1: " << P1.transpose() << std::endl;
        std::cout << "P2: " << P2.transpose() << std::endl;
    }

    // compute the conic coefficients
    A = nx * nx - s2;
    B = -nx * ny;
    C = ny * ny - s2;
    D = -A * sx - B * sy + nx * nz * sz;
    E = -C * sy - B * sx - ny * nz * sz;
    F = F * F - s2 * (sx * sx + sy * sy + sz * sz);

    return 1; // this reflector is visible on the screen, we can proceed to compute the intersection points with the screen and draw the spot on the screen
}

int KLines::getY(double x, double &y1, double &y2) const
{
    double b = B * x + E;
    double delta = b * b - C * (F + x * (A * x + 2 * D));

    if (delta < 0)
        return 0; // no intersection

    double sqrt_delta = sqrt(delta);
    y1 = (-b + sqrt_delta) / C;
    y2 = (-b - sqrt_delta) / C;
    return 1; // intersection exists
}

int KLines::getX(double y, double &x1, double &x2) const
{
    double b = B * y + D;
    double delta = b * b - A * (F + y * (C * y + 2 * E));

    if (delta < 0)
        return 0; // no intersection

    double sqrt_delta = sqrt(delta);
    x1 = (-b + sqrt_delta) / A;
    x2 = (-b - sqrt_delta) / A;
    return 1; // intersection exists
}

void KLines::draw(CImg<unsigned char> &img) const
{
    unsigned char white[] = {255, 255, 255};
     // let's do the work:
    if(isHorizontal()) {
        // we can iterate over the x coordinates of the image and compute the corresponding y coordinates of the conic section defined by the k_lines parameters, and then draw a line between the two intersection points to visualize the spot on the screen
        for (int x = 0; x < img.width() ; x++) {
            double y1, y2;
            if(getY(x, y1, y2)) { // get the corresponding y coordinates of the conic section defined by the k_lines parameters for the given x coordinate
                if(y1 >= 0 && y1 < img.height()) {
                    img.draw_point(x, y1, white); 
                }
                if(y2 >= 0 && y2 < img.height()) {
                    img.draw_point(x, y2, white); 
                }
            }
        }
    } else {
        // we can iterate over the y coordinates of the image and compute the corresponding x coordinates of the conic section defined by the k_lines parameters, and then draw a line between the two intersection points to visualize the spot on the screen
        for (int y = 0; y < img.height() ; y++) {
            double x1, x2;
            if(getX(y, x1, x2)) {  // get the corresponding x coordinates of the conic section defined by the k_lines parameters for the given y coordinate
                if(x1 >= 0 && x1 < img.width()) {
                    img.draw_point(x1, y, white); 
                }
                if(x2 >= 0 && x2 < img.width()) {
                    img.draw_point(x2, y, white); 
                }
            }
        }
    }
}



double KLines::getIntegral(const Pattern &pattern) const
{
    // TODO : implement the function to compute the integral of the pattern along the KLines, which will be used as input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern
    double integral = 0.0;

    if (is_horizontal)
    {
        for (int x = 0; x < pattern.getWidth(); x++)
        {
            double y1, y2;
            if (getY(x, y1, y2))
            { // get the corresponding y coordinates of the conic section defined by the k_lines parameters for the given x coordinate
                if (y1 >= 0 && y1 < pattern.getHeight())
                {
                    integral += pattern.getPattern().linear_atXY(x, y1);
                }
                if (y2 >= 0 && y2 < pattern.getHeight())
                {
                    integral += pattern.getPattern().linear_atXY(x, y2);
                }
            }
        }
    }
    else
    {
        for (int y = 0; y < pattern.getHeight(); y++)
        {
            double x1, x2;
            if (getX(y, x1, x2))
            { // get the corresponding x coordinates of the conic section defined by the k_lines parameters for the given y coordinate
                if (x1 >= 0 && x1 < pattern.getWidth())
                {
                    integral += pattern.getPattern().linear_atXY(x1, y);
                }
                if (x2 >= 0 && x2 < pattern.getWidth())
                {
                    integral += pattern.getPattern().linear_atXY(x2, y);
                }
            }
        }
    }

    return integral; // placeholder value, replace with actual integral computation
}