/*
 *      ABCluster
 *
 * Copyright (c) 2015 Jun Zhang
 * All rights reserved.
 *
 * The source code in this file is part of ABCluster, and is provided under 
 * a licence and may be used, copied, transmitted, or stored only in accord 
 * with that licence.
 *
 */

#ifndef    __RIGIDMOLECULE_H__
#define    __RIGIDMOLECULE_H__

#include "lbfgs.h"
#include <string>
#include <vector>

namespace abcluster {

using namespace std;

class ExtCoord {
public:    
    ExtCoord(void) { /* Nothing to do */ }
    ExtCoord(double tphiZ1, double tphiY2, double tphiZ3, double tmX, double tmY, double tmZ) : phiZ1(tphiZ1), phiY2(tphiY2), phiZ3(tphiZ3), mX(tmX), mY(tmY), mZ(tmZ) { /* Nothing to do */ }
    ~ExtCoord(void) { /* Nothing to do */ }
    double phiZ1;
    double phiY2;
    double phiZ3;
    double mX;
    double mY;
    double mZ;
};

class RigidMolecule {
public:
    RigidMolecule(void);   
    ~RigidMolecule(void);

    // natoms 
    // Comments
    // C1  x1  y1  z1
    // ...
    // Cna xna yna zna
    void LoadMolecule(const string& fn);

    int GetNAtoms(void) const;
    const string& GetComments(void) const;
    const vector<double>& GetInternalCoordinates(void) const;
    const vector<string>& GetSymbols(void) const;
    
    void SetExternalCoordinates(const ExtCoord& teCoord);
    void GetExternalCoordinates(ExtCoord& teCoord) const;
    
    RigidMolecule& operator = (const RigidMolecule& rm);

    // labCoord:      x0 y0 z0 x1 y1 z1 ...   
    void GetLabQuantity(vector<double>& labCoord,  vector<double>& labDerivative) const;
    void GetLabQuantity(const double* x, vector<double>& labCoord,  vector<double>& labDerivative) const;

    // Get eCoords from labCoords
    void SolveExternalCoordinates(const double* labCoords, ExtCoord& eCoords) const;

private:
    int natoms;
    string comments;
    vector<double> iCoord; // x0, y0, z0, x1, y1, z1 ...    
    vector<string> symbols;
    ExtCoord eCoord;

    void MoveToCenter(void);

private:
    static int progress(void *instance, const lbfgsfloatval_t *x, 
            const lbfgsfloatval_t *g, const lbfgsfloatval_t fx,
            const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm,
            const lbfgsfloatval_t step, int n, int k, int ls);
    static lbfgsfloatval_t evaluate(void *instance, 
            const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, 
            const lbfgsfloatval_t step);

};

}

#endif //  __RIGIDMOLECULE_H__
