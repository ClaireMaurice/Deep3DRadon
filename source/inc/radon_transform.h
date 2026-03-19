#pragma once

#include <vector>
#include <string>
#include "pattern.h"
#include "microscope.h"
#include "crystal.h"
#include "source_point.h"
#include "local_radon_box.h"
#include "eigen3/Eigen/Core"

class RadonTransform
{
public:
    // Constructeur avec paramètres par défaut
    RadonTransform(bool precompute = true, 
                   int box_size_precompute = 64, 
                   int box_size_final = 16, 
                   double increment_precompute = 0.001, 
                   double increment_final = 0.0005) 
        : m_radon_precompute(precompute), 
          m_radon_box_size_precompute(box_size_precompute), 
          m_radon_box_size_final(box_size_final), 
          m_radon_increment_precompute(increment_precompute), 
          m_radon_increment_final(increment_final) {}

    ~RadonTransform();  

    // Métrologie et Sauvegarde
    void compute(Microscope& microscope, Crystal& crystal, SourcePoint& viewPoint, Pattern& pattern, const std::string& output_folder = ".");
    void save(const std::string &basename, Crystal &crystal);
    void dump() const;

    // --- RÉCUPÉRATION DE LA NORMALE RECENTRÉE ---
    // Cette méthode parcourt les réflecteurs actifs pour trouver la normale 
    // réellement utilisée lors de la capture de la boîte.
    Eigen::Vector3d getFinalNormal(int reflectorIndex) const {
        for (size_t i = 0; i < activeReflectors.size(); ++i) {
            if (activeReflectors[i] == reflectorIndex) {
                // On renvoie la normale stockée dans la RadonBox correspondante
                return m_radonBoxes[i].get_finalProjectorNormal().head<3>();
            }
        }
        // Si le plan n'était pas assez intense (rejeté par Claire), on renvoie 0.
        // Cela permet à l'IA d'identifier les plans "vides" dans le binaire.
        return Eigen::Vector3d::Zero();
    }

    // Getter pour l'incrément (utile pour le CSV)
    double getIncrementFinal() const { return m_radon_increment_final; }

private:
    std::vector<LocalRadonBox> m_radonBoxes; // Boîtes Radon calculées
    std::vector<int> activeReflectors;       // Indices des plans ayant survécu au precompute
    
    bool m_radon_precompute;
    int m_radon_box_size_precompute;
    int m_radon_box_size_final;
    double m_radon_increment_precompute;
    double m_radon_increment_final;
};