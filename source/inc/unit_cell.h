#pragma once

#include <vector>
#include "atom.h"
#include "eigen3/Eigen/Core"

class UnitCell {
    public:
        UnitCell() : m_latticeParameters(6, 0.0) {} // default constructor with 6 zeros for a, b, c, alpha, beta, gamma

        UnitCell(const std::vector<double>& latticeParameters, const std::vector<Atom>& atoms)
            : m_latticeParameters(latticeParameters), m_atoms(atoms) {}

        void setLatticeParameters(const std::vector<double>& latticeParameters) { m_latticeParameters = latticeParameters; }
        void setAtoms(const std::vector<Atom>& atoms) { m_atoms = atoms; }    

        const std::vector<double>& getLatticeParameters() const { return m_latticeParameters; }
        const std::vector<Atom>& getAtoms() const { return m_atoms; }

        Eigen::Matrix3d getLatticeMatrix() const;
        Eigen::Matrix3d getReciprocalLatticeMatrix() const;

        static UnitCell getFromLatticeMatrix(const Eigen::Matrix3d& latticeMatrix);

    private:
        std::vector<double> m_latticeParameters; // a, b, c, alpha, beta, gamma
        std::vector<Atom> m_atoms; // list of atomic positions in fractional coordinates   
};