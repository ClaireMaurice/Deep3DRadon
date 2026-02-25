#include <iostream>

#include "radon_transform.h"

RadonTransform::RadonTransform() {
    // default constructor
}

RadonTransform::~RadonTransform() {
    // nothing to do for now, but if we had allocated resources we should release them here
}

void RadonTransform::dump() const {
    std::cout << "RadonTransform parameters:" << std::endl;
    std::cout << "View point: " << viewPoint.transpose() << std::endl;
    // nothing to do for now, but if we had allocated resources we should release them here
}

void RadonTransform::apply(const Pattern& pattern) {
    // TODO : implement the Radon transform on the given pattern, we can use the view point to determine the angles and distances for the Radon transform, and then we can apply the Radon transform to the pattern to get the transformed pattern, which we can then use for further processing or visualization
}

