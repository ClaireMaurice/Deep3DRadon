#include <iostream>
#include <fstream>

#define cimg_use_tiff
#define cimg_use_png
#include "CImg.h"

#include "radon_transform.h"
#include "projector.h"

RadonTransform::~RadonTransform()
{
    // nothing to do
}

void RadonTransform::compute(Microscope &microscope, Crystal &crystal, SourcePoint &viewPoint, Pattern &pattern, const std::string &output_folder)
{
    Projector projector;
    projector.buildProjectorList(microscope, crystal, viewPoint);

    m_radonBoxes.clear(); 
    activeReflectors.clear(); 
    int currentReflectorIndex = 0; 

    for (auto projectorNormal : projector.getProjectorList())
    {
        LocalRadonBox radonBox(m_radon_box_size_final, m_radon_increment_final); 

        if(m_radon_precompute)
        {
            if(radonBox.precompute(projectorNormal, viewPoint, pattern, m_radon_box_size_precompute, m_radon_increment_precompute))
            {
                if (radonBox.compute(projectorNormal, viewPoint, pattern))
                {                                                                                                                    
                    activeReflectors.push_back(currentReflectorIndex); 
                    m_radonBoxes.push_back(radonBox);
                }
            }
        }
        else {
            if (radonBox.compute(projectorNormal, viewPoint, pattern))
            {                                                                                                                    
                activeReflectors.push_back(currentReflectorIndex); 
                m_radonBoxes.push_back(radonBox);
            }
        }
        currentReflectorIndex++; 
    }
}

void RadonTransform::dump() const
{
    std::cout << "RadonTransform parameters:" << std::endl;
    std::cout << "Number of local Radon boxes: " << m_radonBoxes.size() << std::endl;
}

// NOTE: On a SUPPRIMÉ RadonTransform::getFinalNormal d'ici car elle est déjà définie dans le .h

void RadonTransform::save(const std::string &basename, Crystal &crystal)
{
    std::string bin_filename = basename + ".bin";
    std::ofstream out(bin_filename, std::ios::binary);

    if (!out) {
        std::cerr << "Erreur : Impossible de creer le fichier " << bin_filename << std::endl;
        return;
    }

    std::vector<Reflector> reflectors = crystal.getReflectors();
    size_t total_expected_plans = reflectors.size();
    size_t box_data_points = m_radon_box_size_final * m_radon_box_size_final * m_radon_box_size_final;

    size_t active_ptr = 0;

    for (size_t i = 0; i < total_expected_plans; ++i) 
    {
        if (active_ptr < activeReflectors.size() && activeReflectors[active_ptr] == (int)i) 
        {
            const auto &box = m_radonBoxes[active_ptr].get_boxImage();
            out.write(reinterpret_cast<const char*>(box.data()), box.size() * sizeof(double));
            active_ptr++;
        } 
        else 
        {
            std::vector<double> empty_box(box_data_points, 0.0);
            out.write(reinterpret_cast<const char*>(empty_box.data()), empty_box.size() * sizeof(double));
        }
    }

    out.close();
}