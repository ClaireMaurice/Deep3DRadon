#pragma once

#include "orientation.h"

class Detector;

class Microscope
{
public:
    // Default constructor
    Microscope();
    ~Microscope();    

    void setAcceleratingVoltage(double voltage_kV);
    void setTiltAngle(double tiltAngle_deg, int tiltAxis);
    void setDetector(Detector* detector);

    double getAcceleratingVoltage() const { return m_acceleratingVoltage; }
    double getElectronWavelength() const { return m_lambda; }
    Detector* getDetector() const { return m_detector; }



private:
    double m_acceleratingVoltage; // in kV
    double m_lambda; // in nm, electron wavelength
    double m_tiltAngle; // in degrees, tilt angle of the sample
    int m_TiltAxis; // 0 for Xm, 1 for Ym

    Quat m_qMS; // the orientation of the microscope frame with respect to the sample frame, expressed as a quaternion
    Detector* m_detector; // pointer to the detector
};

