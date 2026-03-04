#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include "config_handler.h"
#include "json.hpp"

using json = nlohmann::json;

ConfigHandler::ConfigHandler()
    : m_output_folder("../output"),
      m_accelerating_voltage_kV(20.0),
      m_tilt_angle_deg(70.0),
      m_tilt_axis(0),
      m_num_pixels_x(1344),
      m_num_pixels_y(1024),
      m_pixel_size(23.47),
      m_crystal_element("Cu"),
      m_crystal_structure("FCC"),
      m_lattice_parameter(3.615),
      m_source_euler(0.0, 0.0, 0.0),
      m_source_position(0.5, 0.7, 0.5),
      m_deformationGradient(Eigen::Matrix3d::Identity()),
      m_view_euler(0.0, 0.0, 0.0),
      m_view_position(0.5, 0.7, 0.5),
      m_radon_precompute(true),
      m_radon_box_size_precompute(64),
      m_radon_box_size_final(16),
      m_radon_increment_precompute(0.001),
      m_radon_increment_final(0.0005)
{
}

bool ConfigHandler::load(const std::string &configFile)
{
    std::ifstream file(configFile);
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open config file '" << configFile << "'." << std::endl;
        return false;
    }

    try
    {
        json config = json::parse(file);

        // Parse output folder
        if (config.contains("output_folder"))
        {
            m_output_folder = config["output_folder"];
        }

        // Parse microscope parameters
        if (config.contains("microscope"))
        {
            auto &mic = config["microscope"];
            if (mic.contains("accelerating_voltage_kV"))
                m_accelerating_voltage_kV = mic["accelerating_voltage_kV"];
            if (mic.contains("tilt_angle_deg"))
                m_tilt_angle_deg = mic["tilt_angle_deg"];
            if (mic.contains("tilt_axis"))
                m_tilt_axis = mic["tilt_axis"];
        }

        // Parse detector parameters
        if (config.contains("detector"))
        {
            auto &det = config["detector"];
            if (det.contains("num_pixels_x"))
                m_num_pixels_x = det["num_pixels_x"];
            if (det.contains("num_pixels_y"))
                m_num_pixels_y = det["num_pixels_y"];
            if (det.contains("pixel_size"))
                m_pixel_size = det["pixel_size"];
        }

        // Parse crystal parameters
        if (config.contains("crystal"))
        {
            auto &crys = config["crystal"];
            if (crys.contains("element"))
                m_crystal_element = crys["element"];
            if (crys.contains("structure_type"))
                m_crystal_structure = crys["structure_type"];
            if (crys.contains("lattice_parameter"))
                m_lattice_parameter = crys["lattice_parameter"];
        }

        // Parse source point parameters
        if (config.contains("source_point"))
        {
            auto &src = config["source_point"];
            if (src.contains("euler_angles"))
            {
                double phi1 = 0.0, Phi = 0.0, phi2 = 0.0; // default values
                auto &euler = src["euler_angles"];
                if (euler.contains("phi1"))
                   phi1 = euler["phi1"];
                if (euler.contains("Phi"))
                    Phi = euler["Phi"];
                if (euler.contains("phi2"))
                    phi2 = euler["phi2"];
                m_source_euler = Euler(phi1, Phi, phi2);
            }
            if (src.contains("position"))
            {
                auto &pos = src["position"];
                if (pos.contains("x"))
                    m_source_position.x() = pos["x"];
                if (pos.contains("y"))
                    m_source_position.y() = pos["y"];
                if (pos.contains("z"))
                    m_source_position.z() = pos["z"];
            }
            if (src.contains("deformation_gradient"))
            {
                auto &def = src["deformation_gradient"];
                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        m_deformationGradient(i, j) = def[i][j];
                    }
                }
            }
        }

        // Parse view point parameters
        if (config.contains("view_point"))
        {
            auto &view = config["view_point"];
            if (view.contains("euler_angles"))
            {
                auto &euler = view["euler_angles"];
                double phi1 = 0.0, Phi = 0.0, phi2 = 0.0; // default values
                if (euler.contains("phi1"))
                    phi1 = euler["phi1"];
                if (euler.contains("Phi"))
                    Phi = euler["Phi"];
                if (euler.contains("phi2"))
                    phi2 = euler["phi2"];
                m_view_euler = Euler(phi1, Phi, phi2);
            }
            if (view.contains("position"))
            {
                auto &pos = view["position"];
                if (pos.contains("x"))
                    m_view_position.x() = pos["x"];
                if (pos.contains("y"))
                    m_view_position.y() = pos["y"];
                if (pos.contains("z"))
                    m_view_position.z() = pos["z"];
            }
        }

        // Parse Radon transform parameters
        if (config.contains("radon_transform"))
        {
            auto &radon = config["radon_transform"];
            if (radon.contains("precompute"))
                m_radon_precompute = radon["precompute"];
            if (radon.contains("box_size_precompute"))
                m_radon_box_size_precompute = radon["box_size_precompute"];
            if (radon.contains("box_size_final"))
                m_radon_box_size_final = radon["box_size_final"];
            if (radon.contains("increment_precompute"))
                m_radon_increment_precompute = radon["increment_precompute"];
            if (radon.contains("increment_final"))
                m_radon_increment_final = radon["increment_final"];
        }

        file.close();
        std::cout << "Config file '" << configFile << "' loaded successfully." << std::endl;
        return true;
    }
    catch (json::exception &e)
    {
        std::cerr << "Error parsing config file: " << e.what() << std::endl;
        return false;
    }
}

const Euler &ConfigHandler::getSourceEuler() const { return m_source_euler; }
const Eigen::Vector3d &ConfigHandler::getSourcePosition() const { return m_source_position; }
const Eigen::Matrix3d &ConfigHandler::getDeformationGradient() const { return m_deformationGradient; }
const Euler &ConfigHandler::getViewEuler() const { return m_view_euler; }
const Eigen::Vector3d &ConfigHandler::getViewPosition() const { return m_view_position; }

void ConfigHandler::printParameters() const
{
    std::cout << "\n=== Configuration Parameters ===" << std::endl;
    std::cout << "Output folder: " << m_output_folder << std::endl;

    std::cout << "\nMicroscope:" << std::endl;
    std::cout << "  Accelerating voltage: " << m_accelerating_voltage_kV << " kV" << std::endl;
    std::cout << "  Tilt angle: " << m_tilt_angle_deg << " deg, axis: " << m_tilt_axis << std::endl;

    std::cout << "Detector:" << std::endl;
    std::cout << "  Pixels: " << m_num_pixels_x << "x" << m_num_pixels_y << std::endl;
    std::cout << "  Pixel size: " << m_pixel_size << " nm" << std::endl;

    std::cout << "Crystal:" << std::endl;
    std::cout << "  Element: " << m_crystal_element << ", Structure: " << m_crystal_structure << std::endl;
    std::cout << "  Lattice parameter: " << m_lattice_parameter << " Å" << std::endl;

    std::cout << "Source point:" << std::endl;
    std::cout << "  Euler angles: (" << m_source_euler.phi1() << ", " << m_source_euler.phi() << ", " << m_source_euler.phi2() << ")" << std::endl;
    std::cout << "  Position: (" << m_source_position.x() << ", " << m_source_position.y() << ", " << m_source_position.z() << ")" << std::endl;
    std::cout << "  Deformation gradient:" << std::endl;
    for (int i = 0; i < 3; i++)
    {
        std::cout << "    ";
        for (int j = 0; j < 3; j++)
        {
            std::cout << m_deformationGradient(i, j) << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "View point:" << std::endl;
    std::cout << "  Euler angles: (" << m_view_euler.phi1() << ", " << m_view_euler.phi() << ", " << m_view_euler.phi2() << ")" << std::endl;
    std::cout << "  Position: (" << m_view_position.x() << ", " << m_view_position.y() << ", " << m_view_position.z() << ")" << std::endl;

    std::cout << "Radon transform:" << std::endl;
    std::cout << "  Precompute: " << (m_radon_precompute ? "Yes" : "No") << std::endl;
    std::cout << "  Box size for precompute: " << m_radon_box_size_precompute << std::endl;
    std::cout << "  Box size for final: " << m_radon_box_size_final << std::endl;
    std::cout << "  Increment for precompute: " << m_radon_increment_precompute << std::endl;
    std::cout << "  Increment for final: " << m_radon_increment_final << std::endl;

    std::cout << "==============================\n"
              << std::endl;
}

bool ConfigHandler::createOutputFolder() const
{
    if (mkdir(m_output_folder.c_str(), 0755) == 0 || errno == EEXIST)
    {
        return true;
    }
    std::cerr << "Error: Could not create output folder '" << m_output_folder << "'." << std::endl;
    return false;
}
