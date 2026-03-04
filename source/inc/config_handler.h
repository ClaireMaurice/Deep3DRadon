#pragma once

#include <string>
#include <iostream>
#include "eigen3/Eigen/Core"

#include "orientation.h"

class ConfigHandler
{
public:
    ConfigHandler();
    ~ConfigHandler() = default;

    // Load and parse configuration from JSON file
    bool load(const std::string &configFile);

    // Print configuration to console
    void printParameters() const;

    // Create output folder
    bool createOutputFolder() const;

    // Getters for output folder
    const std::string &getOutputFolder() const { return m_output_folder; }

    // Getters for microscope parameters
    double getAcceleratingVoltage() const { return m_accelerating_voltage_kV; }
    double getTiltAngle() const { return m_tilt_angle_deg; }
    int getTiltAxis() const { return m_tilt_axis; }

    // Getters for detector parameters
    int getNumPixelsX() const { return m_num_pixels_x; }
    int getNumPixelsY() const { return m_num_pixels_y; }
    double getPixelSize() const { return m_pixel_size; }

    // Getters for crystal parameters
    const std::string &getCrystalElement() const { return m_crystal_element; }
    const std::string &getCrystalStructure() const { return m_crystal_structure; }
    double getLatticeParameter() const { return m_lattice_parameter; }

    // Getters for source point parameters
    const Euler &getSourceEuler() const;
    const Eigen::Vector3d &getSourcePosition() const;
    const Eigen::Matrix3d &getDeformationGradient() const;

    // Getters for view point parameters
    const Euler &getViewEuler() const;
    const Eigen::Vector3d &getViewPosition() const;

    // Getters for Radon transform parameters
    bool getRadonPrecompute() const { return m_radon_precompute; }
    int getRadonBoxSizePrecompute() const { return m_radon_box_size_precompute; }
    int getRadonBoxSizeFinal() const { return m_radon_box_size_final; }
    double getRadonIncrementPrecompute() const { return m_radon_increment_precompute; }
    double getRadonIncrementFinal() const { return m_radon_increment_final; }

private:
    // Output folder
    std::string m_output_folder;

    // Microscope parameters
    double m_accelerating_voltage_kV;
    double m_tilt_angle_deg;
    int m_tilt_axis;

    // Detector parameters
    int m_num_pixels_x;
    int m_num_pixels_y;
    double m_pixel_size;

    // Crystal parameters
    std::string m_crystal_element;
    std::string m_crystal_structure;
    double m_lattice_parameter;

    // Source point parameters
    Euler m_source_euler;
    Eigen::Vector3d m_source_position;
    Eigen::Matrix3d m_deformationGradient;

    // View point parameters
    Euler m_view_euler;
    Eigen::Vector3d m_view_position;

    // Radon transform parameters
    bool m_radon_precompute;
    int m_radon_box_size_precompute;
    int m_radon_box_size_final;
    double m_radon_increment_precompute;
    double m_radon_increment_final;
};
