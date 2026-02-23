#include <iostream>
#include "crystal.h"

Crystal::Crystal() : m_unitCell(std::vector<double>{1,1,1,90,90,90}, std::vector<Atom>{}) {
    // default constructor corresponding to a simple cubic unit cell with one atom at the origin
}

Crystal::~Crystal() {
    // nothing to do for now, but if we had allocated resources we should release them here
}
