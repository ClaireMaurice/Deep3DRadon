#pragma once

#include <vector>
#include "atom.h"

class UnitCell {
    public:
        UnitCell(const std::vector<double>& latticeParameters, const std::vector<Atom>& atoms)
            : m_latticeParameters(latticeParameters), m_atoms(atoms) {}

        const std::vector<double>& getLatticeParameters() const { return m_latticeParameters; }
        const std::vector<Atom>& getAtoms() const { return m_atoms; }

    private:
        std::vector<double> m_latticeParameters; // a, b, c, alpha, beta, gamma
        std::vector<Atom> m_atoms; // list of atomic positions in fractional coordinates   
};