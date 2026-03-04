#include <iostream>

#define cimg_use_tiff
#define cimg_use_png
#include "CImg.h"

#include "radon_transform.h"
#include "projector.h"



RadonTransform::~RadonTransform()
{
    // nothing to do for now, but if we had allocated resources we should release them here
}

void RadonTransform::compute(Microscope &microscope, Crystal &crystal, SourcePoint &viewPoint, Pattern &pattern, const std::string &output_folder)
{
    Projector projector;
    projector.buildProjectorList(microscope, crystal, viewPoint);

    m_radonBoxes.clear(); // clear the list of local Radon boxes before computing the new ones for this view point, which will be used to compute the Radon transform of the pattern based on the microscope, crystal and source point parameters

    activeReflectors.clear();      // clear the list of active reflectors before computing the new ones for this view point, which will be used to compute the Radon transform of the pattern based on the microscope, crystal and source point parameters, and which can be used as additional input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern
    int currentReflectorIndex = 0; // this variable will keep track of the index of the current reflector being processed, which will be used to compute the Radon transform of the pattern based on the microscope, crystal and source point parameters, and which can be used as additional input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern

    for (auto projectorNormal : projector.getProjectorList())
    {

        LocalRadonBox radonBox(m_radon_box_size_final, m_radon_increment_final); // create a local Radon box for this projector normal, which will be used to compute the Radon transform of the pattern based on the microscope, crystal and source point parameters, and which will be used as input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern

        std::cout << "\n\nProcessing projector normal " << projectorNormal.transpose() << " for reflector index " << currentReflectorIndex << "." << std::endl;

        if(m_radon_precompute)
        {
            std::cout << "Precomputing local Radon box for projector normal " << projectorNormal.transpose() << " with box size " << m_radon_box_size_precompute << " and increment " << m_radon_increment_precompute << "." << std::endl;
            if(radonBox.precompute(projectorNormal, viewPoint, pattern, m_radon_box_size_precompute, m_radon_increment_precompute))
            {
                std::cout << "Precomputation successful for projector normal " << projectorNormal.transpose() << ", proceeding to compute the Radon transform for this projector normal." << std::endl;
                if (radonBox.compute(projectorNormal, viewPoint, pattern))
                {                                                      // compute the Radon transform of the pattern for this projector normal
                    activeReflectors.push_back(currentReflectorIndex); // add the index of the current reflector to the list of active reflectors for this view point
                    m_radonBoxes.push_back(radonBox);
                    // add the local Radon box for this projector normal to the list of local Radon boxes, which will be used to compute the Radon transform of the pattern based on the microscope, crystal and source point parameters, and which will be used as input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern
                }
            }
            else
            {
                std::cerr << "Warning: precomputation failed for projector normal " << projectorNormal.transpose() << ", skipping this projector normal for the Radon transform computation." << std::endl;
            }
            currentReflectorIndex++; // increment the index of the current reflector after processing each projector normal, whether precomputation was successful or not, to ensure that the indices of the active reflectors correspond to the correct projector normals in the list of local Radon boxes
        }
        else {
            std::cout << "Skipping precomputation for projector normal " << projectorNormal.transpose() << " and using default local Radon box parameters." << std::endl;   
            if (radonBox.compute(projectorNormal, viewPoint, pattern))
            {                                                      // compute the Radon transform of the pattern for this projector normal
                activeReflectors.push_back(currentReflectorIndex); // add the index of the current reflector to the list of active reflectors for this view point
                m_radonBoxes.push_back(radonBox);
                    // add the local Radon box for this projector normal to the list of local Radon boxes, which will be used to compute the Radon transform of the pattern based on the microscope, crystal and source point parameters, and which will be used as input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern
            }
             else
            {
                std::cerr << "Warning: computation failed for projector normal " << projectorNormal.transpose() << ", skipping this projector normal for the Radon transform computation." << std::endl;
            }
            currentReflectorIndex++; // increment the index of the current reflector after processing each projector normal, whether computation was successful or not, to ensure that the indices of the active reflectors correspond to the correct projector normals in the list of local Radon boxes
        }

        
        // increment the index of the current reflector after processing each projector normal

        // std::cout << "Press Enter to continue..." << std::endl; // wait for user input before proceeding, to allow time to visualize the 2D slice of the local Radon box
        // std::cin.get();
    }
}

void RadonTransform::dump() const
{
    std::cout << "RadonTransform parameters:" << std::endl;
    std::cout << "Number of local Radon boxes: " << m_radonBoxes.size() << std::endl;
    std::cout << "Active reflectors indices: ";
    for (const auto &index : activeReflectors)
    {
        std::cout << index << " ";
    }
    std::cout << std::endl;
    // we can also print the parameters of each local Radon box if needed, but for now we will just print the number of local Radon boxes for simplicity
}

// Temporary function to save the Radon transform of the pattern
// Should work more on the format of the saved file :
void RadonTransform::save(const std::string &basename, Crystal &crystal)
{

    // I want to go through the list of local Radon boxes and save each one as a separate file,
    // with a name that includes the miller indices of the reflector handled by the local Radon box
    std::vector<Reflector> reflectors = crystal.getReflectors(); // get the list of reflectors from the crystal, which will be used to generate the filenames for saving the local Radon box images, and which will be used to compute the Radon transform of the pattern based on the microscope, crystal and source point parameters, and which will be used as input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern
    for (size_t i = 0; i < m_radonBoxes.size(); i++)
    {
        const auto &radonBox = m_radonBoxes[i];
        const auto &index = activeReflectors[i];   // get the index of the reflector corresponding to this local Radon box, which will be used to generate the filename for saving the local Radon box image, and which will be used to compute the Radon transform of the pattern based on the microscope, crystal and source point parameters, and which will be used as input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern
        const auto &reflector = reflectors[index]; // get the reflector corresponding to this local Radon box using the index from the list of active reflectors, which will be used to generate the

        // the vector containing the hkl parameters of the reflector
        Eigen::Vector3i hkl = reflector.getHKL();

        // generate a string representation of the hkl parameters of the reflector, which will be used to create the filename for saving the local Radon box image, and which will be used to compute the Radon transform of the pattern based on the microscope, crystal and source point parameters, and which will be used as input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern
        std::string hkl_str = "(" + std::to_string(hkl.x()) + std::to_string(hkl.y()) + std::to_string(hkl.z()) + ")";
        std::string filename = basename + "_radon_box_" + std::to_string(i) + "_" + hkl_str + ".tiff"; // create a filename that includes the index of the local Radon box and the parameters of the corresponding projector normal, which will be used to compute the Radon transform of the pattern based on the microscope, crystal and source point parameters, and which will be used as input for the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern

        radonBox.get_boxImage().save_tiff(filename.c_str()); // save the local Radon box image to a file with the specified name, which can be used for visualization or for input to the neural network to predict the 3D structure of the crystal based on the simulated diffraction pattern
    }
}
