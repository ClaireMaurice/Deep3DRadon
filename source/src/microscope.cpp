#include <iostream>

#include "microscope.h"
#include "detector.h"

Microscope::Microscope() : m_tiltAngle(70), m_TiltAxis(0)
{
    // default accelerating voltage is 20 kV, corresponding to an electron wavelength of 0.0251 nm
    setAcceleratingVoltage(20);
    m_qMS = Quat(-m_tiltAngle, m_TiltAxis);
    m_detector = new Detector();
}


void Microscope::setAcceleratingVoltage(double voltage_kV) {
    m_acceleratingVoltage = voltage_kV;
    // we can compute the electron wavelength using the accelerating voltage
    // using the relativistic formula: lambda = h / sqrt(2*m*e*V*(1+e*V/(2*m*c^2)))
    // where h is Planck's constant, m is the electron mass, e is the electron charge, V is the accelerating voltage, and c is the speed of light
    const double h = 6.62607015e-34; // Planck's constant in J*s
    const double m = 9.10938356e-31; // electron mass in kg
    const double e = 1.602176634e-19; // electron charge in C
    const double c = 299792458; // speed of light in m/s

    double V = m_acceleratingVoltage * 1e3; // convert kV to V
    m_lambda = h / sqrt(2*m*e*V*(1+e*V/(2*m*c*c))) * 1e10; // convert m to A

    // checking the calculated wavelength
    // double expected_lambda = 12.2643/sqrt(V*(1+0.97845e-6*V)); // expected wavelength in A for 20 kV
    // std::cout << "Calculated electron wavelength: " << m_lambda << " A" << std::endl;
    // std::cout << "Expected electron wavelength: " << expected_lambda  << " A" << std::endl; 
}

void Microscope::dump() const {
    std::cout << "Microscope parameters:" << std::endl;
    std::cout << "Accelerating voltage: " << m_acceleratingVoltage << " kV" << std::endl;
    std::cout << "Electron wavelength: " << m_lambda << " nm" << std::endl;
    std::cout << "Tilt angle: " << m_tiltAngle << " degrees" << std::endl;
    std::cout << "Tilt axis: " << (m_TiltAxis == 0 ? "Xm" : "Ym") << std::endl;
}

Microscope::~Microscope() {
    // if the detector was allocated, we should delete it
    if (m_detector) {
        delete m_detector;
    }
}

