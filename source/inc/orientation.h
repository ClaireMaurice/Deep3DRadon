#ifndef ORIENTATION_H
#define ORIENTATION_H

#include "eigen3/Eigen/Core"
//#include "CImg.h"

/* This file and its companion orientation.cpp define and implement several
 * Rotation/Orientation descriptors used in EBSD and materials science in
 * general.
 * Each descriptor is implemented as a separate class
 * This is code legacy from earlier versions of StrainCorrelator
 *
 * Linear algebra entities are taken care of through Eigen3 objects. This
 * might be a little bit overkill...
 * */



/*!
 * \brief The Miller class
 * \details TODO
 */
class Miller
{
public:
    Miller();
    Miller(const Eigen::Vector3i &plan, const Eigen::Vector3i &dir );
    ~Miller() = default;

    // accessors :
    const Eigen::Vector3i * getHkl_p() const {return &hkl;}
    const Eigen::Vector3i * getUvw_p() const {return &uvw;}

    // temporary :
    void printInternalMemory();


private:
    Eigen::Vector3i hkl; //!< the crystallographic plane normal to Zs
    Eigen::Vector3i uvw; //!< the crystallographic direction parallel to Xs
};

std::ostream &operator << (std::ostream &s, const Miller &m);

#endif // ORIENTATION_H
