#pragma once
#include "unit_cell.h"

class Crystal
{
public:
    Crystal();
    ~Crystal(); 

    private:
        // TODO : add lattice parameters, unit cell, etc.
        UnitCell m_unitCell;
        
};