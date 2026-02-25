#include <iostream>

#define cimg_use_tiff
#define cimg_use_png
#include "CImg.h"

#include "radon_transform.h"
#include "projector.h"

#define BOX_SIZE 16
#define INCREMENT 0.001


RadonTransform::RadonTransform() {
    // default constructor
}

RadonTransform::~RadonTransform() {
    // nothing to do for now, but if we had allocated resources we should release them here
}

void RadonTransform::compute(Microscope& microscope, Crystal& crystal, SourcePoint& viewPoint, Pattern& pattern) {

    Projector projector;
    projector.buildProjectorList(microscope, crystal, viewPoint);

    m_radonBoxes.clear(); // clear the list of local Radon boxes before computing the new ones for this view point, which will be used to compute the Radon transform of the pattern based on the microscope, crystal and source point parameters
    
    for (const auto& projectorNormal : projector.getProjectorList()) {

        LocalRadonBox radonBox(BOX_SIZE, INCREMENT); // create a local Radon box for this projector normal, which will be used to compute the Radon transform of the pattern based on the microscope, crystal and source point parameters, and which will be used as input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern
        
        int result = radonBox.compute(projectorNormal, viewPoint, pattern); // compute the Radon transform of the pattern for this projector normal

        if(result) 
            m_radonBoxes.push_back(radonBox); // add the local Radon box for this projector normal to the list of local Radon boxes, which will be used to compute the Radon transform of the pattern based on the microscope, crystal and source point parameters, and which will be used as input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern
    }   
}

void RadonTransform::dump() const {
    std::cout << "RadonTransform parameters:" << std::endl;
    std::cout << "Number of local Radon boxes: " << m_radonBoxes.size() << std::endl;
    // we can also print the parameters of each local Radon box if needed, but for now we will just print the number of local Radon boxes for simplicity
}


// Temporary function to save the Radon transform of the pattern 
// Should work more on the format of the saved file : 
void RadonTransform::save(const std::string& basename, Crystal& crystal) {

    // I want to go through the list of local Radon boxes and save each one as a separate file, 
    //with a name that includes the miller indices of the reflector handled by the local Radon box 
    std::vector<Reflector> reflectors = crystal.getReflectors(); // get the list of reflectors from the crystal, which will be used to generate the filenames for saving the local Radon box images, and which will be used to compute the Radon transform of the pattern based on the microscope, crystal and source point parameters, and which will be used as input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern
    for (size_t i = 0; i < m_radonBoxes.size(); i++) {
        const auto& radonBox = m_radonBoxes[i];
        const auto& reflector = reflectors[i];
        
        // the vector containing the hkl parameters of the reflector
        Eigen::Vector3i hkl = reflector.getHKL();

        // generate a string representation of the hkl parameters of the reflector, which will be used to create the filename for saving the local Radon box image, and which will be used to compute the Radon transform of the pattern based on the microscope, crystal and source point parameters, and which will be used as input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern
        std::string hkl_str = "(" + std::to_string(hkl.x()) +  std::to_string(hkl.y())  + std::to_string(hkl.z()) + ")";
        std::string filename = basename + "_radon_box_" + std::to_string(i) + "_" + hkl_str + ".tiff"; // create a filename that includes the index of the local Radon box and the parameters of the corresponding projector normal, which will be used to compute the Radon transform of the pattern based on the microscope, crystal and source point parameters, and which will be used as input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern 
              
        
        radonBox.get_boxImage().save_tiff(filename.c_str()); // save the local Radon box image to a file with the specified name, which can be used for visualization or for input to the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern
    }
}

