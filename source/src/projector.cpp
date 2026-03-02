#include "projector.h"

#include "orientation.h"


void Projector::buildProjectorList(Microscope& microscope, Crystal& crystal, SourcePoint& source_point) {

    // initialize the projector list to an empty list
    projector_list.clear();

    Quat q_DM = microscope.getDetector()->getOrientation();
    Quat q_MS = microscope.getOrientation();
    Quat q_SC = source_point.getOrientation().Conjugate(); // we need the conjugate of the source point orientation to get the correct transformation from the crystal frame to the detector frame
    Quat q_DC = q_SC * q_MS * q_DM; // the order of multiplication is important here, we need to apply the transformations in the correct order to get the correct orientation of the crystal with respect to the detector
    
    Eigen::Matrix3d R_DC = q_DC.toRotationMatrix();

        
    // to calculate the reciprocal lattice vectors, we need the "local" reciprocal lattice matrix defined for the source point
    Eigen::Matrix3d reciprocalLatticeMatrix = source_point.getUnitCell().getReciprocalLatticeMatrix(); 
    std::cout << "Reciprocal lattice matrix:\n" << reciprocalLatticeMatrix << std::endl;

    for (const auto& reflector : crystal.getReflectors()) {
        
        Eigen::Vector3i hkl = reflector.getHKL();
       
        // compute the reciprocal lattice vector for the given Miller indices, for simplicty a cubic lattice is assumed for now, but this can be easily extended to other lattice types by using the appropriate formulas for the reciprocal lattice vectors based on the lattice parameters of the crystal
        Eigen::Vector3d g_hkl = reciprocalLatticeMatrix * hkl.cast<double>(); // in reciprocal angstroms, since a0 is in angstroms and hkl is dimensionless
        
        // compute the interplanar spacing d_hkl = 1 / |g_hkl|, 
        //and then compute the diffraction angle theta using Bragg's law: n*lambda = 2*d_hkl*sin(theta), 
        //where n is the order of the reflection (we can assume n=1 for now), 
        //lambda is the electron wavelength obtained from the microscope parameters, 
        //and d_hkl is the interplanar spacing calculated from the reciprocal lattice vector.
    
        // beware of dimensions : 
        double d_hkl = 1.0 / g_hkl.norm();  // in angstroms, since g_hkl is in reciprocal angstroms
        double lambda = microscope.getElectronWavelength(); // in Angstroms
        double sin_theta = lambda / (2 * d_hkl); // Bragg's law
        double s2 = sin_theta * sin_theta;

        // transport to detector frame
        g_hkl = R_DC * g_hkl; // apply the transport matrix to get the orientation of the reciprocal lattice vector with respect to the detector frame

        g_hkl.normalize(); // we can normalize the reciprocal lattice vector to get the projection normal, which will be used to determine the visibility of the reflector on the screen based on the microscope and source point parameters
        Eigen::Vector4d projector_normal(g_hkl(0), g_hkl(1), g_hkl(2), s2); // we can store the projection normal as a 4D vector, where the first three components are the reciprocal lattice vector and the fourth component is the sin^2(theta) value, which will be used to determine the visibility of the reflector on the screen based on the microscope and source point parameters
        projector_list.push_back(projector_normal); // add the projection normal to the projector list, which will be used to simulate the diffraction pattern on the detector based on the microscope and source point parameters
    }
}
