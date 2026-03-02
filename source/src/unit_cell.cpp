#include "unit_cell.h"
#include <eigen3/Eigen/Dense>   

namespace
{
    constexpr double kDegToRad = 3.14159265358979323846 / 180.0;
    constexpr double kRadToDeg = 180.0 / 3.14159265358979323846;
    constexpr double kEpsilon = 1e-12;

    // Clamps a double value to a specified range, ensuring the result is no less than minVal and no greater than maxVal. It achieves this by combining std::min and std::max to constrain value within the bounds [minVal, maxVal].
    double clamp(double value, double minVal, double maxVal)
    {
        return std::max(minVal, std::min(maxVal, value));
    }
}

Eigen::Matrix3d UnitCell::getLatticeMatrix() const
{
    const double a = m_latticeParameters[0];
    const double b = m_latticeParameters[1];
    const double c = m_latticeParameters[2];
    const double alpha = m_latticeParameters[3] * kDegToRad;
    const double beta = m_latticeParameters[4] * kDegToRad;
    const double gamma = m_latticeParameters[5] * kDegToRad;

    const double cosAlpha = std::cos(alpha);
    const double cosBeta = std::cos(beta);
    const double cosGamma = std::cos(gamma);
    const double sinAlpha = std::sin(alpha);
    

    // cVec is collinear to c (along z-axis)
    const Eigen::Vector3d cVec{0.0, 0.0, c};

    // bVec lies in the yz-plane, making angle alpha with cVec
    const Eigen::Vector3d bVec{0.0, b * sinAlpha, b * cosAlpha};

    // aVec is computed to satisfy angle constraints
    // aVec makes angle beta with cVec: aVec·cVec = |aVec|*|cVec|*cos(beta)
    // aVec makes angle gamma with bVec: aVec·bVec = |aVec|*|bVec|*cos(gamma)
    const double aZ = a * cosBeta;
    const double aY = a*(cosGamma-cosAlpha*cosBeta)/sinAlpha;
    const double aXSquared = a * a - aY * aY - aZ * aZ;
    const double aX = std::sqrt(std::max(0.0, aXSquared));
    const Eigen::Vector3d aVec{aX, aY, aZ};

    Eigen::Matrix3d matrix;
    matrix.col(0) = aVec;
    matrix.col(1) = bVec;
    matrix.col(2) = cVec;
    return matrix;
}


Eigen::Matrix3d UnitCell::getReciprocalLatticeMatrix() const
{
    Eigen::Matrix3d latticeMatrix = getLatticeMatrix();
    Eigen::Matrix3d reciprocalLatticeMatrix = latticeMatrix.inverse().transpose(); // the reciprocal lattice matrix is given by 2*pi times the inverse of the lattice matrix, transposed to get the correct orientation of the reciprocal lattice vectors with respect to the crystal frame
    return reciprocalLatticeMatrix;
}


UnitCell UnitCell::getFromLatticeMatrix(const Eigen::Matrix3d& latticeMatrix)
{
    // we can calculate the lattice parameters from the lattice matrix using the following formulas:
    // a = |aVec|, b = |bVec|, c = |cVec|
    // alpha = angle between bVec and cVec, beta = angle between aVec and cVec, gamma = angle between aVec and bVec
    const Eigen::Vector3d aVec = latticeMatrix.col(0);
    const Eigen::Vector3d bVec = latticeMatrix.col(1);
    const Eigen::Vector3d cVec = latticeMatrix.col(2);

    const double a = aVec.norm();
    const double b = bVec.norm();
    const double c = cVec.norm();

    const double alpha = std::acos(clamp(bVec.dot(cVec) / (b * c), -1.0, 1.0)) * kRadToDeg;
    const double beta = std::acos(clamp(aVec.dot(cVec) / (a * c), -1.0, 1.0)) * kRadToDeg;
    const double gamma = std::acos(clamp(aVec.dot(bVec) / (a * b), -1.0, 1.0)) * kRadToDeg;

    std::vector<double> latticeParameters{a, b, c, alpha, beta, gamma};
    return UnitCell(latticeParameters, {});
}