#include "orientation.h"
#include <iostream>

//! default constructor
Miller::Miller() : hkl(0,0,1), uvw(1,0,0) {}

//! constructor from (hkl) and (uvw) references
Miller::Miller(const Eigen::Vector3i &plan, const Eigen::Vector3i &dir ) :
    hkl(plan),
    uvw(dir)
{}

void Miller::printInternalMemory()
{
    std::cout << "address de hkl : " << &hkl << "\n";
    std::cout << "address de uvw : " << &uvw << "\n";
}


std::ostream &operator << (std::ostream &s, const Miller &m) {
    const Eigen::Vector3i *u = m.getHkl_p();
    const Eigen::Vector3i *v = m.getUvw_p();

    return(s << '(' << (*u)(0) << (*u)(1) << (*u)(2) << ")[" << (*v)(0) << (*v)(1) << (*v)(2) << ']');
}
