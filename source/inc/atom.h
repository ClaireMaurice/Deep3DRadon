#pragma once


class Atom {
public:
    Atom() : Z(0), x(0.0), y(0.0), z(0.0) {} // default constructor
    Atom(int atomicNumber, double xCoord, double yCoord, double zCoord)
        : Z(atomicNumber), x(xCoord), y(yCoord), z(zCoord) {}

    int getAtomicNumber() const { return Z; }
    double getX() const { return x; }
    double getY() const { return y; }
    double getZ() const { return z; }

    void setAtomicNumber(int atomicNumber) { Z = atomicNumber; }
    void setX(double xCoord) { x = xCoord; }
    void setY(double yCoord) { y = yCoord; }
    void setZ(double zCoord) { z = zCoord; }

    
private:
    int Z; // atomic number
    double x, y, z; // coordinates of the atom 
};