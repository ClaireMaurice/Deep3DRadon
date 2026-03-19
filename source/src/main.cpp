#include <iostream>
#include <string>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <random>
#include <iomanip>
#include <sstream>
#include <cmath>

#define cimg_use_png
#define cimg_use_tiff
#include "CImg.h"

#include "config_handler.h"
#include "orientation.h"
#include "microscope.h"
#include "detector.h"
#include "crystal.h"
#include "source_point.h"
#include "pattern.h"
#include "radon_transform.h"

int main(int argc, char *argv[])
{
    // --- 1. Chargement de la configuration ---
    std::string configFile = "../config.json";
    if (argc > 1) configFile = argv[1]; 
    else std::cout << "No config file specified, using default '../config.json'." << std::endl;

    ConfigHandler config;
    if (!config.load(configFile)) {
        std::cerr << "Failed to load configuration file. Aborting." << std::endl;
        return 1;
    }

    if (!config.createOutputFolder()){
        std::cerr << "Failed to create output folder. Aborting." << std::endl;
        return 1;
    }

    // --- 2. Initialisation Microscope & Cristal (Fixes) ---
    Microscope microscope;
    microscope.setAcceleratingVoltage(config.getAcceleratingVoltage());
    microscope.setTilt(config.getTiltAngle(), config.getTiltAxis()); 
    
    Crystal crystal;
    crystal.buildUnitCell(config.getCrystalElement(), config.getCrystalStructure(), config.getLatticeParameter());
    crystal.buildReflectors();

    // Déformation de la maille (Logique de Claire)
    Eigen::Matrix3d latticeMatrix = crystal.getUnitCell().getLatticeMatrix();
    Eigen::Matrix3d deformedLatticeMatrix = config.getDeformationGradient() * latticeMatrix; 
    UnitCell deformedUnitCell = UnitCell::getFromLatticeMatrix(deformedLatticeMatrix); 

    // --- 3. Paramètres de la base de données ---
    const int NB_ORIENTATIONS = 100;
    const int NB_VP_PER_ORIENT = 10;
    
    // MODIFICATION 1 : Passage à une variation de 1/100.
    // Si la base_vp est à ~0.5, 1% de variation correspond à +/- 0.005.
    const double VP_RANGE = 0.005; 

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis_euler(0.0, 360.0);
    std::uniform_real_distribution<> dis_vp(-VP_RANGE, VP_RANGE);

    // --- 4. Préparation du fichier CSV pour l'IA ---
    std::ofstream csvLabels(config.getOutputFolder() + "/labels.csv");
    
    // MODIFICATION 2 : Ajout de la colonne "increment" dans l'en-tête
    csvLabels << "id,vp_x,vp_y,vp_z,euler_1,euler_2,euler_3,increment";
    
    int NB_PLANS = crystal.getReflectors().size();
    for(int n=0; n < NB_PLANS; ++n) {
        csvLabels << ",n" << n << "_x,n" << n << "_y,n" << n << "_z";
    }
    csvLabels << std::endl;

    int global_id = 0;
    Eigen::Vector3d base_vp = config.getViewPosition();

    std::cout << "\n=== Debut de la generation du dataset ( " << (NB_ORIENTATIONS * NB_VP_PER_ORIENT) << " samples ) ===" << std::endl;

    // On initialise RadonTransform
    RadonTransform radonTransform(config.getRadonPrecompute(), config.getRadonBoxSizePrecompute(), config.getRadonBoxSizeFinal(), config.getRadonIncrementPrecompute(), config.getRadonIncrementFinal());

    // --- 5. Boucle Principale ---
    for (int i = 0; i < NB_ORIENTATIONS; ++i) {
        
        double e1 = dis_euler(gen), e2 = dis_euler(gen), e3 = dis_euler(gen);
        Euler random_orient(e1 * M_PI / 180.0, e2 * M_PI / 180.0, e3 * M_PI / 180.0);

        SourcePoint sourcePoint(random_orient, config.getSourcePosition(), deformedUnitCell);
        Pattern pattern;
        pattern.simulate(microscope, crystal, sourcePoint);

        for (int j = 0; j < NB_VP_PER_ORIENT; ++j) {
            
            // Calcul du ViewPoint avec la nouvelle plage de variation (1/100)
            Eigen::Vector3d VP(base_vp.x() + dis_vp(gen), base_vp.y() + dis_vp(gen), base_vp.z() + dis_vp(gen));
            SourcePoint viewPoint(random_orient, VP, crystal.getUnitCell());

            // 1. Calcul de Radon (effectue le precompute/recentrage en interne)
            radonTransform.compute(microscope, crystal, viewPoint, pattern, config.getOutputFolder());

            std::stringstream ss;
            ss << std::setw(6) << std::setfill('0') << global_id;
            std::string id_str = ss.str();
            
            // 2. Sauvegarde du binaire (.bin)
            radonTransform.save(config.getOutputFolder() + "/sample_" + id_str, crystal); 

            // 3. Sauvegarde dans le CSV
            // On écrit l'ID, le VP, l'orientation et l'incrément (échelle physique)
            csvLabels << id_str << "," << VP.x() << "," << VP.y() << "," << VP.z()
                      << "," << e1 << "," << e2 << "," << e3
                      << "," << config.getRadonIncrementFinal(); // <--- Ajout de l'incrément ici

            // MODIFICATION 3 : Sauvegarde des normales RECENTRÉES
            // Ce sont les normales finales utilisées pour poser la boîte Radon.
            for (int n = 0; n < NB_PLANS; ++n) {
                Eigen::Vector3d finalNormal = radonTransform.getFinalNormal(n);
                csvLabels << "," << finalNormal.x() << "," << finalNormal.y() << "," << finalNormal.z();
            }
            csvLabels << std::endl;

            global_id++;
            std::cout << "Progression : " << global_id << " / " << (NB_ORIENTATIONS * NB_VP_PER_ORIENT) << "\r" << std::flush;
        }
    }

    std::cout << "\nGeneration terminee ! Le fichier CSV contient desormais l'increment et les normales de capture." << std::endl;
    csvLabels.close();
    return 0;
}