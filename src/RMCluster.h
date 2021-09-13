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

#ifndef    __RMCLUSTER_H__
#define    __RMCLUSTER_H__

#include "RigidMolecule.h"
#include "lbfgs.h"
#include <string>
#include <vector>
#include <utility>

namespace abcluster {

using namespace std;

class RMCluster {
public:
    RMCluster(void);
    ~RMCluster(void);

    void LoadRMCluster(const string& fn);
    void SaveAsXYZ(const string& fn) const;

    double GetEnergy(void) const;
    void SetEnergy(double tenergy);
    const vector<RigidMolecule>& GetRigidMolecules(void) const;
    int GetNRigidMolecules(void) const;
    const vector<pair<string, int> >& GetComponents(void) const;

    // Warning: eCoord will be filled without allocation!
    void GeteCoordsArray(double* eCoord) const;
    void SeteCoordsArray(const double* eCoord);

    // Operators
    RMCluster& operator = (const RMCluster& rmc);

    int CoarseOpt(const vector<int>& tfix_indices);

private:
    int nrms;
    vector<RigidMolecule> rms;
    double energy;
    vector<pair<string, int> > components;

    // Compute energy and grad.
    void ComputEnergyGrad(double& tenergy, vector<double>& tgrad);
   void ComputeEnergyGradCore(
           int naI, const vector<double>& CartI, const vector<double>& DerivI,
           int naJ, const vector<double>& CartJ, const vector<double>& DerivJ,
           double& energyIJ,
           double& dEdIphiZ1, double& dEdIphiY2, double& dEdIphiZ3, double& dEdImX, double& dEdImY, double& dEdImZ,
           double& dEdJphiZ1, double& dEdJphiY2, double& dEdJphiZ3, double& dEdJmX, double& dEdJmY, double& dEdJmZ
           ) const;

private:
    vector<int> fix_indices; // To fix some coordinates.

    // For L-BFGS optimization. See the information in "lbfgs.h".
    static int progress(void *instance, const lbfgsfloatval_t *x, 
            const lbfgsfloatval_t *g, const lbfgsfloatval_t fx,
            const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm,
            const lbfgsfloatval_t step, int n, int k, int ls);   
    // For L-BFGS optimization. See the information in "lbfgs.h".
    static lbfgsfloatval_t evaluate(void *instance, const lbfgsfloatval_t *x,
            lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);    

};

}

#endif //  __RMCLUSTER_H__
