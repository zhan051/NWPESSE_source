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

#ifndef    __CERESULT_H__
#define    __CERESULT_H__

#include "RMCluster.h"
#include <string>
#include <vector>

namespace abcluster {

using namespace std;

class CEResult {
public:
    // Constructor.
    // tn: number of atoms.
    CEResult(void);
    // Destructor.
    ~CEResult(void);
    
    // Tackle the atoms.
    int Save(const string& fn) const; // Save X to fn.
    int Show(void) const;             // Write X on the screen.
    
    void SetN(int tn);
    void SetEnergy(double tenergy);
    void SetX(const vector<double>& tx);
    void SetSymbols(const vector<string>& tsymbols);
    
    double GetEnergy(void) const;
    
    // Move the coordinates to the orginal-center
    void MoveSC(void);

    // Solve the extent coordinates from the optimzied XYZ coordinates.
    void SolveExtCoords(const RMCluster& rmc, vector<ExtCoord>& eCoords) const;

private:
    int n;                  // Number of atoms.
    int dof;                // DOF
    double energy;          // Energy
    vector<double> x;       // Coordinates
    vector<string> symbols; // Printed symbols.
 
    // Save solutions.
    int SaveGJF(const string& fn) const; // Save as gjf file.
    int SaveXYZ(const string& fn) const; // Save as XYZ file.

};

}

#endif //  __CERESULT_H__
